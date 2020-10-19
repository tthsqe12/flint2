/*
    Copyright (C) 2018 - 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "n_poly.h"
#include "nmod_mpoly.h"
#include "fq_nmod_mpoly.h"

/*
    For each j, set out[j] to the evaluation of A at x_i = alpha[i] (i != j)
    i.e. if nvars = 3
        out[0] = A(x, alpha[1], alpha[2])
        out[1] = A(alpha[0], x, alpha[2])
        out[2] = A(alpha[0], alpha[1], x)

    If ignore[j] is nonzero, then out[j] need not be calculated, probably
    because we shouldn't calculate it in dense form.
*/
static void nmod_mpoly_evals(
    n_poly_struct * out,
    const int * ignore,
    const nmod_mpoly_t A,
    ulong * Amin_exp,
    ulong * Amax_exp,
    ulong * Astride,
    mp_limb_t * alpha,
    const nmod_mpoly_ctx_t ctx)
{
    slong i, j;
    slong nvars = ctx->minfo->nvars;
    ulong mask = (-UWORD(1)) >> (FLINT_BITS - A->bits);
    slong * offsets, * shifts;
    slong N = mpoly_words_per_exp_sp(A->bits, ctx->minfo);
    ulong * varexps;
    n_poly_struct * caches;

    offsets = FLINT_ARRAY_ALLOC(2*nvars, slong);
    shifts = offsets + nvars;
    varexps = FLINT_ARRAY_ALLOC(nvars, ulong);
    caches = FLINT_ARRAY_ALLOC(3*nvars, n_poly_struct);
    for (j = 0; j < nvars; j++)
    {
        n_poly_zero(out + j);
        mpoly_gen_offset_shift_sp(offsets + j, shifts + j, j, A->bits, ctx->minfo);
        n_poly_init(caches + 3*j + 0);
        n_poly_init(caches + 3*j + 1);
        n_poly_init(caches + 3*j + 2);
        nmod_pow_cache_start(alpha[j], caches + 3*j + 0,
                                           caches + 3*j + 1, caches + 3*j + 2);
    }

    for (i = 0; i < A->length; i++)
    {
        mp_limb_t meval;
        mp_limb_t t;
        ulong varexp;

        meval = A->coeffs[i];
        for (j = 0; j < nvars; j++)
        {
            varexp = ((A->exps + N*i)[offsets[j]]>>shifts[j])&mask;

            FLINT_ASSERT((Astride[j] == 0 && varexp == Amin_exp[j]) ||
                                     (varexp - Amin_exp[j]) % Astride[j] == 0);

            varexps[j] = Astride[j] < 2 ? varexp - Amin_exp[j] :
                                         (varexp - Amin_exp[j])/Astride[j];

            t = nmod_pow_cache_mulpow_ui(meval, varexps[j], caches + 3*j + 0,
                         caches + 3*j + 1, caches + 3*j + 2, ctx->ffinfo->mod);

            FLINT_ASSERT(t == nmod_mul(meval, nmod_pow_ui(alpha[j], varexps[j],
                                         ctx->ffinfo->mod), ctx->ffinfo->mod));
            meval = t;
        }

        for (j = 0; j < nvars; j++)
        {
            varexp = varexps[j];

            if (ignore[j])
                continue;

            n_poly_fit_length(out + j, varexp + 1);

            while (out[j].length <= varexp)
            {
                out[j].coeffs[out[j].length] = 0;
                out[j].length++;
            }

            t = nmod_pow_cache_mulpow_neg_ui(meval, varexp, caches + 3*j + 0,
                         caches + 3*j + 1, caches + 3*j + 2, ctx->ffinfo->mod);

            FLINT_ASSERT(t == nmod_mul(meval, nmod_pow_ui(nmod_inv(alpha[j],
              ctx->ffinfo->mod), varexp, ctx->ffinfo->mod), ctx->ffinfo->mod));

            out[j].coeffs[varexp] = nmod_add(out[j].coeffs[varexp], t,
                                                             ctx->ffinfo->mod);
        }
    }

    for (j = 0; j < nvars; j++)
        _n_poly_normalise(out + j);

    for (j = 0; j < 3*nvars; j++)
        n_poly_clear(caches + j);

    flint_free(offsets);
    flint_free(varexps);
    flint_free(caches);
}


static void nmod_mpoly_evals_lgprime(
    n_fq_poly_struct * out,
    const int * ignore,
    const nmod_mpoly_t A,
    ulong * Amin_exp,
    ulong * Amax_exp,
    ulong * Astride,
    const nmod_mpoly_ctx_t smctx,
    const fq_nmod_struct * alpha,
    const fq_nmod_ctx_t lgctx)
{
    slong d = fq_nmod_ctx_degree(lgctx);
    slong i, j;
    slong nvars = smctx->minfo->nvars;
    ulong mask = (-UWORD(1)) >> (FLINT_BITS - A->bits);
    slong * offsets, * shifts;
    slong N = mpoly_words_per_exp_sp(A->bits, smctx->minfo);
    ulong * varexps;
    n_poly_struct * caches;
    mp_limb_t * t = FLINT_ARRAY_ALLOC(2*d, mp_limb_t);
    mp_limb_t * meval = t + d;

    offsets = FLINT_ARRAY_ALLOC(2*nvars, slong);
    shifts = offsets + nvars;
    varexps = FLINT_ARRAY_ALLOC(nvars, ulong);
    caches = FLINT_ARRAY_ALLOC(3*nvars, n_poly_struct);
    for (j = 0; j < nvars; j++)
    {
        n_fq_poly_zero(out + j);
        mpoly_gen_offset_shift_sp(offsets + j, shifts + j, j, A->bits,
                                                                 smctx->minfo);
        n_poly_init(caches + 3*j + 0);
        n_poly_init(caches + 3*j + 1);
        n_poly_init(caches + 3*j + 2);
        n_fq_pow_cache_start_fq_nmod(alpha + j, caches + 3*j + 0,
                                    caches + 3*j + 1, caches + 3*j + 2, lgctx);
    }

    for (i = 0; i < A->length; i++)
    {
        ulong varexp;

        _n_fq_set_nmod(meval, A->coeffs[i], d);

        for (j = 0; j < nvars; j++)
        {
            varexp = ((A->exps + N*i)[offsets[j]]>>shifts[j])&mask;

            FLINT_ASSERT((Astride[j] == 0 && varexp == Amin_exp[j]) ||
                                     (varexp - Amin_exp[j]) % Astride[j] == 0);

            varexps[j] = Astride[j] < 2 ? varexp - Amin_exp[j] :
                                         (varexp - Amin_exp[j])/Astride[j];

            n_fq_pow_cache_mulpow_ui(meval, meval, varexps[j], caches + 3*j + 0,
                                    caches + 3*j + 1, caches + 3*j + 2, lgctx);
        }

        for (j = 0; j < nvars; j++)
        {
            varexp = varexps[j];

            if (ignore[j])
                continue;

            n_poly_fit_length(out + j, d*(varexp + 1));

            while (out[j].length <= varexp)
            {
                _n_fq_zero(out[j].coeffs + d*out[j].length, d);
                out[j].length++;
            }

            n_fq_pow_cache_mulpow_neg_ui(t, meval, varexp, caches + 3*j + 0,
                                    caches + 3*j + 1, caches + 3*j + 2, lgctx);

            n_fq_add(out[j].coeffs + d*varexp, out[j].coeffs + d*varexp, t, lgctx);
        }
    }

    for (j = 0; j < nvars; j++)
        _n_fq_poly_normalise(out + j, d);

    for (j = 0; j < 3*nvars; j++)
        n_poly_clear(caches + j);

    flint_free(offsets);
    flint_free(varexps);
    flint_free(caches);
    flint_free(t);
}


void mpoly_gcd_info_set_estimates_nmod_mpoly(
    mpoly_gcd_info_t I,
    const nmod_mpoly_t A,
    const nmod_mpoly_t B,
    const nmod_mpoly_ctx_t ctx,
    const thread_pool_handle * handles,
    slong num_handles)
{
    int try_count = 0;
    slong nvars = ctx->minfo->nvars;
    slong i, j;
    n_poly_t Geval;
    n_poly_struct * Aevals, * Bevals;
    mp_limb_t * alpha;
    flint_rand_t state;
    slong ignore_limit;
    int * ignore;

    flint_randinit(state);

    ignore = FLINT_ARRAY_ALLOC(nvars, int);
    alpha  = FLINT_ARRAY_ALLOC(nvars, mp_limb_t);
    Aevals = FLINT_ARRAY_ALLOC(nvars, n_poly_struct);
    Bevals = FLINT_ARRAY_ALLOC(nvars, n_poly_struct);

    n_poly_init(Geval);
    for (j = 0; j < nvars; j++)
    {
        n_poly_init(Aevals + j);
        n_poly_init(Bevals + j);
    }

    ignore_limit = A->length/4096 + B->length/4096;
    ignore_limit = FLINT_MAX(WORD(9999), ignore_limit);
    I->Gdeflate_deg_bounds_are_nice = 1;
    for (j = 0; j < nvars; j++)
    {
        if (I->Adeflate_deg[j] > ignore_limit ||
            I->Bdeflate_deg[j] > ignore_limit)
        {
            ignore[j] = 1;
            I->Gdeflate_deg_bounds_are_nice = 0;
        }
        else
        {
            ignore[j] = 0;
        }
    }

try_again:

    if (++try_count > 10)
    {
        I->Gdeflate_deg_bounds_are_nice = 0;
        for (j = 0; j < nvars; j++)
        {
            I->Gdeflate_deg_bound[j] = FLINT_MIN(I->Adeflate_deg[j],
                                                 I->Bdeflate_deg[j]);
            I->Gterm_count_est[j] = 1 + I->Gdeflate_deg_bound[j]/2;
        }

        goto cleanup;
    }

    for (j = 0; j < nvars; j++)
        alpha[j] = n_urandint(state, ctx->ffinfo->mod.n - 1) + 1;

    nmod_mpoly_evals(Aevals, ignore, A, I->Amin_exp, I->Amax_exp, I->Gstride, alpha, ctx);
    nmod_mpoly_evals(Bevals, ignore, B, I->Bmin_exp, I->Bmax_exp, I->Gstride, alpha, ctx);

    for (j = 0; j < nvars; j++)
    {
        if (ignore[j])
        {
            I->Gdeflate_deg_bound[j] = FLINT_MIN(I->Adeflate_deg[j],
                                                 I->Bdeflate_deg[j]);
            I->Gterm_count_est[j] = 1 + I->Gdeflate_deg_bound[j]/2;
        }
        else
        {
            if (I->Adeflate_deg[j] != n_poly_degree(Aevals + j) ||
                I->Bdeflate_deg[j] != n_poly_degree(Bevals + j))
            {
                goto try_again;
            }

            n_poly_mod_gcd(Geval, Aevals + j, Bevals + j, ctx->ffinfo->mod);

            I->Gterm_count_est[j] = 0;
            I->Gdeflate_deg_bound[j] = n_poly_degree(Geval);
            for (i = I->Gdeflate_deg_bound[j]; i >= 0; i--)
                I->Gterm_count_est[j] += (Geval->coeffs[i] != 0);
        }
    }

cleanup:

    n_poly_clear(Geval);
    for (j = 0; j < nvars; j++)
    {
        n_poly_clear(Aevals + j);
        n_poly_clear(Bevals + j);
    }

    flint_free(ignore);
    flint_free(alpha);
    flint_free(Aevals);
    flint_free(Bevals);

    flint_randclear(state);

    return;
}

void mpoly_gcd_info_set_estimates_nmod_mpoly_lgprime(
    mpoly_gcd_info_t I,
    const nmod_mpoly_t A,
    const nmod_mpoly_t B,
    const nmod_mpoly_ctx_t smctx,
    const thread_pool_handle * handles,
    slong num_handles)
{
    slong nvars = smctx->minfo->nvars;
    int try_count = 0;
    slong i, j;
    n_fq_poly_t Geval;
    n_fq_poly_struct * Aevals, * Bevals;
    fq_nmod_struct * alpha;
    flint_rand_t state;
    slong ignore_limit;
    int * ignore;
    fq_nmod_mpoly_ctx_t lgctx;
    slong lgd;

    flint_randinit(state);

    lgd = WORD(20)/(FLINT_BIT_COUNT(smctx->ffinfo->mod.n));
    lgd = FLINT_MAX(WORD(2), lgd);
    fq_nmod_mpoly_ctx_init_deg(lgctx, nvars, ORD_LEX, smctx->ffinfo->mod.n, lgd);

    ignore = FLINT_ARRAY_ALLOC(nvars, int);
    alpha = FLINT_ARRAY_ALLOC(nvars, fq_nmod_struct);
    Aevals = FLINT_ARRAY_ALLOC(nvars, n_fq_poly_struct);
    Bevals = FLINT_ARRAY_ALLOC(nvars, n_fq_poly_struct);
    for (j = 0; j < nvars; j++)
    {
        n_fq_poly_init(Aevals + j);
        n_fq_poly_init(Bevals + j);
        fq_nmod_init(alpha + j, lgctx->fqctx);
    }

    n_fq_poly_init(Geval);

    ignore_limit = A->length/4096 + B->length/4096;
    ignore_limit = FLINT_MAX(WORD(9999), ignore_limit);
    I->Gdeflate_deg_bounds_are_nice = 1;
    for (j = 0; j < nvars; j++)
    {
        if (I->Adeflate_deg[j] > ignore_limit ||
            I->Bdeflate_deg[j] > ignore_limit)
        {
            ignore[j] = 1;
            I->Gdeflate_deg_bounds_are_nice = 0;
        }
        else
        {
            ignore[j] = 0;
        }
    }

try_again:

    if (++try_count > 10)
    {
        I->Gdeflate_deg_bounds_are_nice = 0;
        for (j = 0; j < nvars; j++)
        {
            I->Gdeflate_deg_bound[j] = FLINT_MIN(I->Adeflate_deg[j],
                                                 I->Bdeflate_deg[j]);
            I->Gterm_count_est[j] = 1 + I->Gdeflate_deg_bound[j]/2;
        }

        goto cleanup;
    }

    for (j = 0; j < nvars; j++)
    {
        fq_nmod_rand(alpha + j, state, lgctx->fqctx);
        if (fq_nmod_is_zero(alpha + j, lgctx->fqctx))
            fq_nmod_one(alpha + j, lgctx->fqctx);
    }

    nmod_mpoly_evals_lgprime(Aevals, ignore, A, I->Amin_exp, I->Amax_exp,
                                       I->Gstride, smctx, alpha, lgctx->fqctx);
    nmod_mpoly_evals_lgprime(Bevals, ignore, B, I->Bmin_exp, I->Bmax_exp,
                                       I->Gstride, smctx, alpha, lgctx->fqctx);

    for (j = 0; j < nvars; j++)
    {
        if (ignore[j])
        {
            I->Gdeflate_deg_bound[j] = FLINT_MIN(I->Adeflate_deg[j],
                                                 I->Bdeflate_deg[j]);
            I->Gterm_count_est[j] = 1 + I->Gdeflate_deg_bound[j]/2;
        }
        else
        {
            if (I->Adeflate_deg[j] != n_fq_poly_degree(Aevals + j) ||
                I->Bdeflate_deg[j] != n_fq_poly_degree(Bevals + j))
            {
                lgd++;
                fq_nmod_mpoly_ctx_change_modulus(lgctx, lgd);
                goto try_again;
            }

            n_fq_poly_gcd(Geval, Aevals + j, Bevals + j, lgctx->fqctx);

            I->Gterm_count_est[j] = 0;
            I->Gdeflate_deg_bound[j] = n_fq_poly_degree(Geval);
            for (i = I->Gdeflate_deg_bound[j]; i >= 0; i--)
            {
                I->Gterm_count_est[j] += (Geval->coeffs[i] != 0);
            }
        }
    }

cleanup:

    n_fq_poly_clear(Geval);
    for (j = 0; j < nvars; j++)
    {
        n_fq_poly_clear(Aevals + j);
        n_fq_poly_clear(Bevals + j);
        fq_nmod_clear(alpha + j, lgctx->fqctx);
    }
    flint_free(alpha);
    flint_free(Aevals);
    flint_free(Bevals);
    flint_free(ignore);

    fq_nmod_mpoly_ctx_clear(lgctx);

    flint_randclear(state);

    return;
}


/*********************** Easy when B is a monomial ***************************/
static int _try_monomial_gcd(
    nmod_mpoly_t G, flint_bitcnt_t Gbits,
    const nmod_mpoly_t A,
    const nmod_mpoly_t B,
    const nmod_mpoly_ctx_t ctx)
{
    slong i;
    fmpz * minAfields, * minAdegs, * minBdegs;
    TMP_INIT;

    FLINT_ASSERT(A->length > 0);
    FLINT_ASSERT(B->length == 1);

    TMP_START;

    /* get the field-wise minimum of A */
    minAfields = (fmpz *) TMP_ALLOC(ctx->minfo->nfields*sizeof(fmpz));
    for (i = 0; i < ctx->minfo->nfields; i++)
        fmpz_init(minAfields + i);
    mpoly_min_fields_fmpz(minAfields, A->exps, A->length, A->bits, ctx->minfo);

    /* unpack to get the min degrees of each variable in A */
    minAdegs = (fmpz *) TMP_ALLOC(ctx->minfo->nvars*sizeof(fmpz));
    for (i = 0; i < ctx->minfo->nvars; i++)
        fmpz_init(minAdegs + i);
    mpoly_get_monomial_ffmpz_unpacked_ffmpz(minAdegs, minAfields, ctx->minfo);

    /* get the degree of each variable in B */
    minBdegs = (fmpz *) TMP_ALLOC(ctx->minfo->nvars*sizeof(fmpz));
    for (i = 0; i < ctx->minfo->nvars; i++)
        fmpz_init(minBdegs + i);
    mpoly_get_monomial_ffmpz(minBdegs, B->exps, B->bits, ctx->minfo);

    /* compute the degree of each variable in G */
    _fmpz_vec_min_inplace(minBdegs, minAdegs, ctx->minfo->nvars);

    nmod_mpoly_fit_length_reset_bits(G, 1, Gbits, ctx);
    mpoly_set_monomial_ffmpz(G->exps, minBdegs, Gbits, ctx->minfo);
    G->coeffs[0] = 1;
    G->length = 1;

    for (i = 0; i < ctx->minfo->nfields; i++)
    {
        fmpz_clear(minAfields + i);
    }
    for (i = 0; i < ctx->minfo->nvars; i++)
    {
        fmpz_clear(minAdegs + i);
        fmpz_clear(minBdegs + i);
    }

    TMP_END;

    return 1;
}


/********************** See if cofactors are monomials ***********************/
static int _try_monomial_cofactors(
    nmod_mpoly_t G, flint_bitcnt_t Gbits,
    const nmod_mpoly_t A,
    const nmod_mpoly_t B,
    const nmod_mpoly_ctx_t ctx)
{
    int success;
    slong i, j;
    slong NA, NG;
    slong nvars = ctx->minfo->nvars;
    fmpz * Abarexps, * Bbarexps, * Texps;
    mp_limb_t a0inv;
    nmod_mpoly_t T;
    TMP_INIT;

    FLINT_ASSERT(A->length > 0);
    FLINT_ASSERT(B->length > 0);

    if (A->length != B->length)
        return 0;

    for (i = A->length - 1; i > 0; i--)
    {
        success = (nmod_mul(A->coeffs[0], B->coeffs[i], ctx->ffinfo->mod)
                == nmod_mul(B->coeffs[0], A->coeffs[i], ctx->ffinfo->mod));
        if (!success)
            goto cleanup;
    }

    TMP_START;

    Abarexps = (fmpz *) TMP_ALLOC(3*nvars*sizeof(fmpz));
    Bbarexps = Abarexps + 1*nvars;
    Texps    = Abarexps + 2*nvars;
    for (j = 0; j < nvars; j++)
    {
        fmpz_init(Abarexps + j);
        fmpz_init(Bbarexps + j);
        fmpz_init(Texps + j);
    }

    success = mpoly_monomial_cofactors(Abarexps, Bbarexps, A->exps, A->bits,
                                      B->exps, B->bits, A->length, ctx->minfo);
    if (!success)
        goto cleanup_tmp;

    nmod_mpoly_init3(T, A->length, Gbits, ctx);
    NG = mpoly_words_per_exp(Gbits, ctx->minfo);
    NA = mpoly_words_per_exp(A->bits, ctx->minfo);
    a0inv = nmod_inv(A->coeffs[0], ctx->ffinfo->mod);
    T->length = A->length;
    for (i = 0; i < A->length; i++)
    {
        mpoly_get_monomial_ffmpz(Texps, A->exps + NA*i, A->bits, ctx->minfo);
        _fmpz_vec_sub(Texps, Texps, Abarexps, nvars);
        mpoly_set_monomial_ffmpz(T->exps + NG*i, Texps, Gbits, ctx->minfo);
        T->coeffs[i] = nmod_mul(A->coeffs[i], a0inv, ctx->ffinfo->mod);
    }
    nmod_mpoly_swap(G, T, ctx);
    nmod_mpoly_clear(T, ctx);

    success = 1;

cleanup_tmp:

    for (j = 0; j < nvars; j++)
    {
        fmpz_clear(Abarexps + j);
        fmpz_clear(Bbarexps + j);
        fmpz_clear(Texps + j);
    }

    TMP_END;

cleanup:

    return success;
}


/********* Assume B has length one when converted to univar format ***********/
static int _try_missing_var(
    nmod_mpoly_t G, flint_bitcnt_t Gbits,
    slong var,
    const nmod_mpoly_t A, ulong Ashift,
    const nmod_mpoly_t B, ulong Bshift,
    const nmod_mpoly_ctx_t ctx)
{
    int success;
    slong i;
    nmod_mpoly_t tG;
    nmod_mpoly_univar_t Ax;

    nmod_mpoly_init(tG, ctx);
    nmod_mpoly_univar_init(Ax, ctx);

    nmod_mpoly_to_univar(Ax, A, var, ctx);

    FLINT_ASSERT(Ax->length > 0);
    success = _nmod_mpoly_gcd_threaded_pool(tG, Gbits, B, Ax->coeffs + 0,
                                                                 ctx, NULL, 0);
    if (!success)
        goto cleanup;

    for (i = 1; i < Ax->length; i++)
    {
        success = _nmod_mpoly_gcd_threaded_pool(tG, Gbits, tG, Ax->coeffs + i,
                                                                 ctx, NULL, 0);
        if (!success)
            goto cleanup;
    }

    nmod_mpoly_swap(G, tG, ctx);
    _mpoly_gen_shift_left(G->exps, G->bits, G->length,
                                   var, FLINT_MIN(Ashift, Bshift), ctx->minfo);

cleanup:

    nmod_mpoly_clear(tG, ctx);
    nmod_mpoly_univar_clear(Ax, ctx);

    return success;
}


/************************ See if B divides A ********************************/
static int _try_divides(
    nmod_mpoly_t G,
    const nmod_mpoly_t A,
    const nmod_mpoly_t BB,
    const nmod_mpoly_ctx_t ctx,
    const thread_pool_handle * handles,
    slong num_handles)
{
    int success = 0;
    nmod_mpoly_t Q, B, M;

    nmod_mpoly_init(Q, ctx);
    nmod_mpoly_init(B, ctx);
    nmod_mpoly_init(M, ctx);

    /* BB = M*B */
    nmod_mpoly_term_content(M, BB, ctx);
    nmod_mpoly_divides(B, BB, M, ctx);

    if (_nmod_mpoly_divides_threaded_pool(Q, A, B, ctx, handles, num_handles))
    {
        /* gcd(Q*B, M*B) */
        _try_monomial_gcd(G, A->bits, Q, M, ctx);
        nmod_mpoly_mul(G, G, B, ctx);
        success = 1;
    }

    nmod_mpoly_clear(Q, ctx);
    nmod_mpoly_clear(B, ctx);
    nmod_mpoly_clear(M, ctx);

    return success;
}


/********************** Hit A and B with zippel ******************************/
static int _try_zippel(
    nmod_mpoly_t G,
    const nmod_mpoly_t A,
    const nmod_mpoly_t B,
    const mpoly_gcd_info_t I,
    const nmod_mpoly_ctx_t ctx)
{
    slong i, k;
    slong m = I->mvars;
    int success;
    mpoly_zipinfo_t zinfo;
    flint_bitcnt_t wbits;
    flint_rand_t randstate;
    nmod_mpoly_ctx_t uctx;
    nmod_mpolyu_t Au, Bu, Gu, Abaru, Bbaru;
    nmod_mpoly_t Ac, Bc, Gc;

    FLINT_ASSERT(A->bits <= FLINT_BITS);
    FLINT_ASSERT(B->bits <= FLINT_BITS);

    if (!I->can_use_zippel)
        return 0;

    FLINT_ASSERT(m >= WORD(2));
    FLINT_ASSERT(A->length > 0);
    FLINT_ASSERT(B->length > 0);

    flint_randinit(randstate);

    /* interpolation will continue in m variables */
    mpoly_zipinfo_init(zinfo, m);

    /* uctx is context for Z[y_1,...,y_{m-1}]*/
    nmod_mpoly_ctx_init(uctx, m - 1, ORD_LEX, ctx->ffinfo->mod.n);

    /* fill in a valid zinfo->perm and degrees */
    for (i = 0; i < m; i++)
    {
        k = I->zippel_perm[i];
        zinfo->perm[i] = k;
        zinfo->Adegs[i] = I->Adeflate_deg[k];
        zinfo->Bdegs[i] = I->Bdeflate_deg[k];
        FLINT_ASSERT(I->Adeflate_deg[k] != 0);
        FLINT_ASSERT(I->Bdeflate_deg[k] != 0);
    }

    wbits = FLINT_MAX(A->bits, B->bits);

    nmod_mpolyu_init(Au, wbits, uctx);
    nmod_mpolyu_init(Bu, wbits, uctx);
    nmod_mpolyu_init(Gu, wbits, uctx);
    nmod_mpolyu_init(Abaru, wbits, uctx);
    nmod_mpolyu_init(Bbaru, wbits, uctx);
    nmod_mpoly_init3(Ac, 0, wbits, uctx);
    nmod_mpoly_init3(Bc, 0, wbits, uctx);
    nmod_mpoly_init3(Gc, 0, wbits, uctx);

    nmod_mpoly_to_mpolyu_perm_deflate_threaded_pool(Au, uctx, A, ctx, zinfo->perm,
                                             I->Amin_exp, I->Gstride, NULL, 0);
    nmod_mpoly_to_mpolyu_perm_deflate_threaded_pool(Bu, uctx, B, ctx, zinfo->perm,
                                             I->Bmin_exp, I->Gstride, NULL, 0);

    success = nmod_mpolyu_content_mpoly_threaded_pool(Ac, Au, uctx, NULL, 0);
    success = success &&
              nmod_mpolyu_content_mpoly_threaded_pool(Bc, Bu, uctx, NULL, 0);
    if (!success)
        goto cleanup;

    nmod_mpolyu_divexact_mpoly_inplace(Au, Ac, uctx);
    nmod_mpolyu_divexact_mpoly_inplace(Bu, Bc, uctx);

    /* after removing content, degree bounds in zinfo are still valid bounds */
    success = nmod_mpolyu_gcdm_zippel(Gu, Abaru, Bbaru, Au, Bu,
                                                       uctx, zinfo, randstate);
    if (!success)
        goto cleanup;

    success = _nmod_mpoly_gcd_threaded_pool(Gc, wbits, Ac, Bc, uctx, NULL, 0);
    if (!success)
        goto cleanup;

    nmod_mpolyu_mul_mpoly_inplace(Gu, Gc, uctx);

    nmod_mpoly_from_mpolyu_perm_inflate(G, I->Gbits, ctx, Gu, uctx,
                                         zinfo->perm, I->Gmin_exp, I->Gstride);
    success = 1;

cleanup:

    nmod_mpolyu_clear(Au, uctx);
    nmod_mpolyu_clear(Bu, uctx);
    nmod_mpolyu_clear(Gu, uctx);
    nmod_mpolyu_clear(Abaru, uctx);
    nmod_mpolyu_clear(Bbaru, uctx);
    nmod_mpoly_clear(Ac, uctx);
    nmod_mpoly_clear(Bc, uctx);
    nmod_mpoly_clear(Gc, uctx);

    nmod_mpoly_ctx_clear(uctx);

    mpoly_zipinfo_clear(zinfo);

    flint_randclear(randstate);

    return success;
}


static int _try_zippel3(
    nmod_mpoly_t G,
    const nmod_mpoly_t A,
    const nmod_mpoly_t B,
    const mpoly_gcd_info_t I,
    const nmod_mpoly_ctx_t ctx)
{
    slong i, k;
    slong m = I->mvars;
    int success;
    flint_bitcnt_t wbits;
    nmod_mpoly_ctx_t lctx;
    nmod_mpoly_t Al, Bl, Gl, Abarl, Bbarl;
    nmod_mpoly_t Al_lc, Bl_lc, Ac, Bc, Gc, Gamma;
    slong * tmp, * Gl_degs, * Al_degs, * Bl_degs, * Gamma_degs, * Gguess;
    slong max_degree;

    FLINT_ASSERT(A->bits <= FLINT_BITS);
    FLINT_ASSERT(B->bits <= FLINT_BITS);
    FLINT_ASSERT(A->length > 0);
    FLINT_ASSERT(B->length > 0);

    if (!I->can_use_bma)
        return 0;

    FLINT_ASSERT(m >= WORD(3));

    tmp = FLINT_ARRAY_ALLOC(5*m, slong);
    Al_degs   = tmp + 1*m;
    Bl_degs   = tmp + 2*m;
    Gamma_degs = tmp + 3*m;
    Gl_degs   = tmp + 4*m;

    nmod_mpoly_ctx_init(lctx, m, ORD_LEX, ctx->ffinfo->mod.n);

    max_degree = 0;
    for (i = 0; i < m; i++)
    {
        k = I->bma_perm[i];

        Gl_degs[i] = I->Gdeflate_deg_bound[k];

        Al_degs[i] = I->Adeflate_deg[k];
        max_degree = FLINT_MAX(max_degree, Al_degs[i]);

        Bl_degs[i] = I->Bdeflate_deg[k];
        max_degree = FLINT_MAX(max_degree, Bl_degs[i]);
    }

    wbits = 1 + FLINT_BIT_COUNT(max_degree);
    wbits = mpoly_fix_bits(wbits, lctx->minfo);
    FLINT_ASSERT(wbits <= FLINT_BITS);

    nmod_mpoly_init3(Al, 0, wbits, lctx);
    nmod_mpoly_init3(Bl, 0, wbits, lctx);
    nmod_mpoly_init3(Gl, 0, wbits, lctx);
    nmod_mpoly_init3(Abarl, 0, wbits, lctx);
    nmod_mpoly_init3(Bbarl, 0, wbits, lctx);
    nmod_mpoly_init3(Ac, 0, wbits, lctx);
    nmod_mpoly_init3(Bc, 0, wbits, lctx);
    nmod_mpoly_init3(Gc, 0, wbits, lctx);
    nmod_mpoly_init3(Gamma, 0, wbits, lctx);
    nmod_mpoly_init3(Al_lc, 0, wbits, lctx);
    nmod_mpoly_init3(Bl_lc, 0, wbits, lctx);

    nmod_mpoly_to_mpolyl_perm_deflate(Al, lctx, A, ctx, I->bma_perm, I->Amin_exp, I->Gstride);
    nmod_mpoly_to_mpolyl_perm_deflate(Bl, lctx, B, ctx, I->bma_perm, I->Bmin_exp, I->Gstride);

    success = nmod_mpolyl_content(Ac, Al, 2, lctx) &&
              nmod_mpolyl_content(Bc, Bl, 2, lctx) && 
              nmod_mpoly_gcd(Gc, Ac, Bc, lctx);
    if (!success)
        goto cleanup;

    nmod_mpoly_degrees_si(tmp, Ac, lctx);
    for (i = 0; i < m; i++)
        Al_degs[i] -= tmp[i];

    success = nmod_mpoly_divides(Al, Al, Ac, lctx);
    FLINT_ASSERT(success);

    nmod_mpoly_degrees_si(tmp, Bc, lctx);
    for (i = 0; i < m; i++)
        Bl_degs[i] -= tmp[i];

    success = nmod_mpoly_divides(Bl, Bl, Bc, lctx);
    FLINT_ASSERT(success);

    nmod_mpoly_degrees_si(tmp, Gc, lctx);
    for (i = 0; i < m; i++)
        Gl_degs[i] -= tmp[i];

    nmod_mpoly_repack_bits_inplace(Al, wbits, lctx);
    nmod_mpoly_repack_bits_inplace(Bl, wbits, lctx);
    nmod_mpolyl_lead_coeff(Al_lc, Al, 2, lctx);
    nmod_mpolyl_lead_coeff(Bl_lc, Bl, 2, lctx);
    success = nmod_mpoly_gcd(Gamma, Al_lc, Bl_lc, lctx);
    if (!success)
        goto cleanup;
    nmod_mpoly_repack_bits_inplace(Gamma, wbits, lctx);

    nmod_mpoly_degrees_si(Gamma_degs, Gamma, lctx);

    Gguess = I->Gdeflate_deg_bounds_are_nice ? Gl_degs : NULL;

    FLINT_ASSERT(Gl->bits == wbits);
    FLINT_ASSERT(Abarl->bits == wbits);
    FLINT_ASSERT(Bbarl->bits == wbits);
    success = nmod_mpolyl_gcd_zippel_smprime(Gl, Gguess, Abarl, Bbarl,
                            Al, Al_degs, Bl, Bl_degs, Gamma, Gamma_degs, lctx);
    if (!success)
    {
        success = nmod_mpolyl_gcd_zippel_lgprime(Gl, Gguess, Abarl, Bbarl,
                            Al, Al_degs, Bl, Bl_degs, Gamma, Gamma_degs, lctx);
        if (!success)
            goto cleanup;
    }

    nmod_mpoly_mul(Gl, Gl, Gc, lctx);

    nmod_mpoly_from_mpolyl_perm_inflate(G, I->Gbits, ctx, Gl, lctx,
                                         I->bma_perm, I->Gmin_exp, I->Gstride);
    success = 1;

cleanup:

    nmod_mpoly_clear(Al, lctx);
    nmod_mpoly_clear(Bl, lctx);
    nmod_mpoly_clear(Gl, lctx);
    nmod_mpoly_clear(Abarl, lctx);
    nmod_mpoly_clear(Bbarl, lctx);
    nmod_mpoly_clear(Ac, lctx);
    nmod_mpoly_clear(Bc, lctx);
    nmod_mpoly_clear(Gc, lctx);
    nmod_mpoly_clear(Gamma, lctx);
    nmod_mpoly_clear(Al_lc, lctx);
    nmod_mpoly_clear(Bl_lc, lctx);

    nmod_mpoly_ctx_clear(lctx);

    flint_free(tmp);

    return success;
}


static int _try_zippel2(
    nmod_mpoly_t G,
    const nmod_mpoly_t A,
    const nmod_mpoly_t B,
    const mpoly_gcd_info_t I,
    const nmod_mpoly_ctx_t ctx)
{
    slong i, k;
    slong m = I->mvars;
    int success;
    flint_bitcnt_t wbits;
    nmod_mpoly_ctx_t uctx;
    nmod_mpolyu_t Auu, Buu, Guu, Abaruu, Bbaruu;
    nmod_mpoly_t Ac, Bc, Gc, Gamma;
    slong * tmp, * Guu_degs, * Auu_degs, * Buu_degs, * Gamma_degs, * Gguess;
    slong max_minor_degree;

    FLINT_ASSERT(A->bits <= FLINT_BITS);
    FLINT_ASSERT(B->bits <= FLINT_BITS);
    FLINT_ASSERT(A->length > 0);
    FLINT_ASSERT(B->length > 0);

    if (!I->can_use_bma)
        return 0;

    FLINT_ASSERT(m >= WORD(3));

    tmp = FLINT_ARRAY_ALLOC(5*m, slong);
    Auu_degs   = tmp + 1*m;
    Buu_degs   = tmp + 2*m;
    Gamma_degs = tmp + 3*m;
    Guu_degs   = tmp + 4*m;

    /* uctx is context for Fp[y_2,...,y_{m - 1}]*/
    nmod_mpoly_ctx_init(uctx, m - 2, ORD_LEX, ctx->ffinfo->mod.n);

    max_minor_degree = 0;
    for (i = 2; i < m; i++)
    {
        k = I->bma_perm[i];

        Guu_degs[i - 2] = I->Gdeflate_deg_bound[k];

        Auu_degs[i - 2] = I->Adeflate_deg[k];
        max_minor_degree = FLINT_MAX(max_minor_degree, Auu_degs[i - 2]);

        Buu_degs[i - 2] = I->Bdeflate_deg[k];
        max_minor_degree = FLINT_MAX(max_minor_degree, Buu_degs[i - 2]);
    }

    wbits = 1 + FLINT_BIT_COUNT(max_minor_degree);
    wbits = mpoly_fix_bits(wbits, uctx->minfo);
    FLINT_ASSERT(wbits <= FLINT_BITS);

    nmod_mpolyu_init(Auu, wbits, uctx);
    nmod_mpolyu_init(Buu, wbits, uctx);
    nmod_mpolyu_init(Guu, wbits, uctx);
    nmod_mpolyu_init(Abaruu, wbits, uctx);
    nmod_mpolyu_init(Bbaruu, wbits, uctx);
    nmod_mpoly_init3(Ac, 0, wbits, uctx);
    nmod_mpoly_init3(Bc, 0, wbits, uctx);
    nmod_mpoly_init3(Gc, 0, wbits, uctx);
    nmod_mpoly_init3(Gamma, 0, wbits, uctx);

    nmod_mpoly_to_mpolyuu_perm_deflate_threaded_pool(Auu, uctx, A, ctx,
                   I->bma_perm, I->Amin_exp, I->Gstride, I->Amax_exp, NULL, 0);
    nmod_mpoly_to_mpolyuu_perm_deflate_threaded_pool(Buu, uctx, B, ctx,
                   I->bma_perm, I->Bmin_exp, I->Gstride, I->Bmax_exp, NULL, 0);

    success = nmod_mpolyu_content_mpoly_threaded_pool(Ac, Auu, uctx, NULL, 0);
    success = success &&
              nmod_mpolyu_content_mpoly_threaded_pool(Bc, Buu, uctx, NULL, 0);
    success = success &&
              _nmod_mpoly_gcd_threaded_pool(Gc, wbits, Ac, Bc, uctx, NULL, 0);
    if (!success)
        goto cleanup;

    nmod_mpoly_degrees_si(tmp, Ac, uctx);
    for (i = 0; i < m - 2; i++)
        Auu_degs[i] -= tmp[i];

    nmod_mpolyu_divexact_mpoly_inplace(Auu, Ac, uctx);

    nmod_mpoly_degrees_si(tmp, Bc, uctx);
    for (i = 0; i < m - 2; i++)
        Buu_degs[i] -= tmp[i];

    nmod_mpolyu_divexact_mpoly_inplace(Buu, Bc, uctx);

    nmod_mpoly_degrees_si(tmp, Gc, uctx);
    for (i = 0; i < m - 2; i++)
        Guu_degs[i] -= tmp[i];

    FLINT_ASSERT(Auu->bits == wbits);
    FLINT_ASSERT(Buu->bits == wbits);
    FLINT_ASSERT(Auu->length > 1);
    FLINT_ASSERT(Buu->length > 1);
    FLINT_ASSERT(Ac->bits == wbits);
    FLINT_ASSERT(Bc->bits == wbits);

    success = _nmod_mpoly_gcd_threaded_pool(Gamma, wbits, Auu->coeffs + 0,
                                               Buu->coeffs + 0, uctx, NULL, 0);
    if (!success)
        goto cleanup;

    nmod_mpoly_degrees_si(Gamma_degs, Gamma, uctx);

    /* Auu*Gamma and Buu*Gamma must not overflow */
    for (i = 0; i < m - 2; i++)
    {
        ulong this_deg = (ulong) Gamma_degs[i] + FLINT_MAX(Auu_degs[i], Buu_degs[i]);
        if ((FLINT_BIT_COUNT(this_deg)) >= wbits)
        {
            wbits = mpoly_fix_bits(wbits + 1, uctx->minfo);
            if (wbits > FLINT_BITS)
            {
                success = 0;
                goto cleanup;
            }

            nmod_mpolyu_repack_bits_inplace(Guu, wbits, uctx);
            nmod_mpolyu_repack_bits_inplace(Abaruu, wbits, uctx);
            nmod_mpolyu_repack_bits_inplace(Bbaruu, wbits, uctx);
            nmod_mpolyu_repack_bits_inplace(Auu, wbits, uctx);
            nmod_mpolyu_repack_bits_inplace(Buu, wbits, uctx);
            nmod_mpoly_repack_bits_inplace(Gamma, wbits, uctx);
            nmod_mpoly_repack_bits_inplace(Gc, wbits, uctx);
            break;
        }
    }

    Gguess = I->Gdeflate_deg_bounds_are_nice ? Guu_degs : NULL;

    success = nmod_mpolyuu_gcd_zippel_smprime(Guu, Gguess, Abaruu, Bbaruu,
                        Auu, Auu_degs, Buu, Buu_degs, Gamma, Gamma_degs, uctx);
    if (!success)
    {
        success = nmod_mpolyuu_gcd_zippel_lgprime(Guu, Gguess, Abaruu, Bbaruu,
                        Auu, Auu_degs, Buu, Buu_degs, Gamma, Gamma_degs, uctx);
        if (!success)
            goto cleanup;
    }

    nmod_mpolyu_mul_mpoly_inplace(Guu, Gc, uctx);

    nmod_mpoly_from_mpolyuu_perm_inflate(G, I->Gbits, ctx, Guu, uctx,
                                         I->bma_perm, I->Gmin_exp, I->Gstride);
    success = 1;

cleanup:

    nmod_mpolyu_clear(Auu, uctx);
    nmod_mpolyu_clear(Buu, uctx);
    nmod_mpolyu_clear(Guu, uctx);
    nmod_mpolyu_clear(Abaruu, uctx);
    nmod_mpolyu_clear(Bbaruu, uctx);
    nmod_mpoly_clear(Ac, uctx);
    nmod_mpoly_clear(Bc, uctx);
    nmod_mpoly_clear(Gc, uctx);
    nmod_mpoly_clear(Gamma, uctx);

    nmod_mpoly_ctx_clear(uctx);

    flint_free(tmp);

    return success;
}


/*********************** Hit A and B with brown ******************************/
typedef struct
{
    nmod_mpolyn_struct * Pn;
    const nmod_mpoly_ctx_struct * nctx;
    const nmod_mpoly_struct * P;
    const nmod_mpoly_ctx_struct * ctx;
    const slong * perm;
    const ulong * shift, * stride;
    const thread_pool_handle * handles;
    slong num_handles;
}
_convertn_arg_struct;

typedef _convertn_arg_struct _convertn_arg_t[1];

static void _worker_convertn(void * varg)
{
    _convertn_arg_struct * arg = (_convertn_arg_struct *) varg;

    nmod_mpoly_to_mpolyn_perm_deflate_threaded_pool(arg->Pn, arg->nctx, arg->P, arg->ctx,
           arg->perm, arg->shift, arg->stride, arg->handles, arg->num_handles);
}

static int _try_brown(
    nmod_mpoly_t G,
    const nmod_mpoly_t A,
    const nmod_mpoly_t B,
    mpoly_gcd_info_t I,
    const nmod_mpoly_ctx_t ctx,
    const thread_pool_handle * handles,
    slong num_handles)
{
    int success;
    slong m = I->mvars;
    flint_bitcnt_t wbits;
    nmod_mpoly_ctx_t nctx;
    nmod_mpolyn_t An, Bn, Gn, Abarn, Bbarn;
    nmod_poly_stack_t Sp;

    if (!I->can_use_brown)
        return 0;

    FLINT_ASSERT(m >= 2);
    FLINT_ASSERT(A->bits <= FLINT_BITS);
    FLINT_ASSERT(B->bits <= FLINT_BITS);
    FLINT_ASSERT(A->length > 0);
    FLINT_ASSERT(B->length > 0);

    wbits = FLINT_MAX(A->bits, B->bits);

    nmod_mpoly_ctx_init(nctx, m, ORD_LEX, ctx->ffinfo->mod.n);
    nmod_poly_stack_init(Sp, wbits, nctx);
    nmod_mpolyn_init(An, wbits, nctx);
    nmod_mpolyn_init(Bn, wbits, nctx);
    nmod_mpolyn_init(Gn, wbits, nctx);
    nmod_mpolyn_init(Abarn, wbits, nctx);
    nmod_mpolyn_init(Bbarn, wbits, nctx);

    if (num_handles > 0)
    {
        slong s = mpoly_divide_threads(num_handles, A->length, B->length);
        _convertn_arg_t arg;

        FLINT_ASSERT(s >= 0);
        FLINT_ASSERT(s < num_handles);

        arg->Pn = Bn;
        arg->nctx = nctx;
        arg->P = B;
        arg->ctx = ctx;
        arg->perm = I->brown_perm;
        arg->shift = I->Bmin_exp;
        arg->stride = I->Gstride;
        arg->handles = handles + (s + 1);
        arg->num_handles = num_handles - (s + 1);

        thread_pool_wake(global_thread_pool, handles[s], 0, _worker_convertn, arg);

        nmod_mpoly_to_mpolyn_perm_deflate_threaded_pool(An, nctx, A, ctx,
                       I->brown_perm, I->Amin_exp, I->Gstride, handles + 0, s);

        thread_pool_wait(global_thread_pool, handles[s]);
    }
    else
    {
        nmod_mpoly_to_mpolyn_perm_deflate_threaded_pool(An, nctx, A, ctx,
                              I->brown_perm, I->Amin_exp, I->Gstride, NULL, 0);
        nmod_mpoly_to_mpolyn_perm_deflate_threaded_pool(Bn, nctx, B, ctx,
                              I->brown_perm, I->Bmin_exp, I->Gstride, NULL, 0);
    }

    FLINT_ASSERT(An->bits == wbits);
    FLINT_ASSERT(Bn->bits == wbits);
    FLINT_ASSERT(An->length > 1);
    FLINT_ASSERT(Bn->length > 1);

    success = (num_handles > 0)
        ? nmod_mpolyn_gcd_brown_smprime_threaded_pool(Gn, Abarn, Bbarn, An, Bn,
                                          m - 1, nctx, I, handles, num_handles)
        : nmod_mpolyn_gcd_brown_smprime(Gn, Abarn, Bbarn, An, Bn,
                                                           m - 1, nctx, I, Sp);

    if (!success)
    {
        nmod_mpoly_to_mpolyn_perm_deflate_threaded_pool(An, nctx, A, ctx,
                              I->brown_perm, I->Amin_exp, I->Gstride, NULL, 0);
        nmod_mpoly_to_mpolyn_perm_deflate_threaded_pool(Bn, nctx, B, ctx,
                              I->brown_perm, I->Bmin_exp, I->Gstride, NULL, 0);
        success = nmod_mpolyn_gcd_brown_lgprime(Gn, Abarn, Bbarn, An, Bn,
                                                                  m - 1, nctx);
    }

    if (!success)
        goto cleanup;

    nmod_mpoly_from_mpolyn_perm_inflate(G, I->Gbits, ctx, Gn, nctx,
                                       I->brown_perm, I->Gmin_exp, I->Gstride);
    success = 1;

cleanup:

    nmod_mpolyn_clear(An, nctx);
    nmod_mpolyn_clear(Bn, nctx);
    nmod_mpolyn_clear(Gn, nctx);
    nmod_mpolyn_clear(Abarn, nctx);
    nmod_mpolyn_clear(Bbarn, nctx);
    nmod_poly_stack_clear(Sp);
    nmod_mpoly_ctx_clear(nctx);

    return success;
}


/*
    The function must pack its answer into bits = Gbits <= FLINT_BITS
    Both A and B have to be packed into bits <= FLINT_BITS

    return is 1 for success, 0 for failure.
*/
int _nmod_mpoly_gcd_threaded_pool(
    nmod_mpoly_t G, flint_bitcnt_t Gbits,
    const nmod_mpoly_t A,
    const nmod_mpoly_t B,
    const nmod_mpoly_ctx_t ctx,
    const thread_pool_handle * handles,
    slong num_handles)
{
    int success;
    slong v_in_both;
    slong v_in_either;
    slong v_in_A_only;
    slong v_in_B_only;
    slong j;
    slong nvars = ctx->minfo->nvars;
    mpoly_gcd_info_t I;
#if FLINT_WANT_ASSERT
    nmod_mpoly_t T, Asave, Bsave;
#endif

    if (A->length == 1)
    {
        return _try_monomial_gcd(G, Gbits, B, A, ctx);
    }
    else if (B->length == 1)
    {
        return _try_monomial_gcd(G, Gbits, A, B, ctx);
    }

#if FLINT_WANT_ASSERT
    nmod_mpoly_init(T, ctx);
    nmod_mpoly_init(Asave, ctx);
    nmod_mpoly_init(Bsave, ctx);
    nmod_mpoly_set(Asave, A, ctx);
    nmod_mpoly_set(Bsave, B, ctx);
#endif

    mpoly_gcd_info_init(I, nvars);

    /* entries of I are all now invalid */

    I->Gbits = Gbits;

    mpoly_gcd_info_limits(I->Amax_exp, I->Amin_exp, I->Alead_count,
                      I->Atail_count, A->exps, A->bits, A->length, ctx->minfo);
    mpoly_gcd_info_limits(I->Bmax_exp, I->Bmin_exp, I->Blead_count,
                      I->Btail_count, B->exps, B->bits, B->length, ctx->minfo);

    /* set ess(p) := p/term_content(p) */

    /* check if the cofactors could be monomials, i.e. ess(A) == ess(B) */
    for (j = 0; j < nvars; j++)
    {
        if (I->Amax_exp[j] - I->Amin_exp[j] != I->Bmax_exp[j] - I->Bmin_exp[j])
            goto skip_monomial_cofactors;
    }

    if (_try_monomial_cofactors(G, I->Gbits, A, B, ctx))
        goto successful;

skip_monomial_cofactors:

    mpoly_gcd_info_stride(I->Gstride,
            A->exps, A->bits, A->length, I->Amax_exp, I->Amin_exp,
            B->exps, B->bits, B->length, I->Bmax_exp, I->Bmin_exp, ctx->minfo);

    for (j = 0; j < nvars; j++)
    {
        ulong t = I->Gstride[j];

        if (t == 0)
        {
            FLINT_ASSERT(I->Amax_exp[j] == I->Amin_exp[j] ||
                         I->Bmax_exp[j] == I->Bmin_exp[j]);
        }
        else
        {
            FLINT_ASSERT((I->Amax_exp[j] - I->Amin_exp[j]) % t == 0);
            FLINT_ASSERT((I->Bmax_exp[j] - I->Bmin_exp[j]) % t == 0);
        }

        I->Adeflate_deg[j] = t == 0 ? 0 : (I->Amax_exp[j] - I->Amin_exp[j])/t;
        I->Bdeflate_deg[j] = t == 0 ? 0 : (I->Bmax_exp[j] - I->Bmin_exp[j])/t;
        I->Gmin_exp[j] = FLINT_MIN(I->Amin_exp[j], I->Bmin_exp[j]);
    }

    /*
        The following are now valid:
            I->Amax_exp, I->Amin_exp, I->Alead_count, I->Atail_count,
            I->Bmax_exp, I->Bmin_exp, I->Blead_count, I->Btail_count,
            I->Gstride
            I->Adeflate_deg
            I->Bdeflate_deg
            I->Gmin_exp
    */

    /* check if ess(A) and ess(B) have a variable v_in_both in common */
    v_in_both = -WORD(1);
    for (j = 0; j < nvars; j++)
    {
        if (I->Amax_exp[j] > I->Amin_exp[j] && I->Bmax_exp[j] > I->Bmin_exp[j])
        {
            v_in_both = j;
            break;
        }
    }
    if (v_in_both == -WORD(1))
    {
        /*
            The variables in ess(A) and ess(B) are disjoint.
            gcd is trivial to compute.
        */

calculate_trivial_gcd:

        nmod_mpoly_fit_length_reset_bits(G, 1, Gbits, ctx);
        mpoly_set_monomial_ui(G->exps, I->Gmin_exp, Gbits, ctx->minfo);
        G->coeffs[0] = 1;
        _nmod_mpoly_set_length(G, 1, ctx);

        goto successful;
    }

    /* check if ess(A) and ess(B) depend on another variable v_in_either */
    FLINT_ASSERT(0 <= v_in_both);
    FLINT_ASSERT(v_in_both < nvars);

    v_in_either = -WORD(1);
    for (j = 0; j < nvars; j++)
    {
        if (j == v_in_both)
            continue;

        if (I->Amax_exp[j] > I->Amin_exp[j] || I->Bmax_exp[j] > I->Bmin_exp[j])
        {
            v_in_either = j;
            break;
        }
    }

    if (v_in_either == -WORD(1))
    {
        /*
            The ess(A) and ess(B) depend on only one variable v_in_both
            Calculate gcd using univariates
        */
        nmod_poly_t a, b, g;

        nmod_poly_init_mod(a, ctx->ffinfo->mod);
        nmod_poly_init_mod(b, ctx->ffinfo->mod);
        nmod_poly_init_mod(g, ctx->ffinfo->mod);
        _nmod_mpoly_to_nmod_poly_deflate(a, A, v_in_both,
                                                 I->Amin_exp, I->Gstride, ctx);
        _nmod_mpoly_to_nmod_poly_deflate(b, B, v_in_both,
                                                 I->Bmin_exp, I->Gstride, ctx);
        nmod_poly_gcd(g, a, b);
        _nmod_mpoly_from_nmod_poly_inflate(G, Gbits, g, v_in_both,
                                                 I->Gmin_exp, I->Gstride, ctx);
        nmod_poly_clear(a);
        nmod_poly_clear(b);
        nmod_poly_clear(g);

        goto successful;
    }

    /* check if there is a variable in ess(A) that is not in ess(B) */
    v_in_A_only = -WORD(1);
    v_in_B_only = -WORD(1);
    for (j = 0; j < nvars; j++)
    {
        if (I->Amax_exp[j] > I->Amin_exp[j] && I->Bmax_exp[j] == I->Bmin_exp[j])
        {
            v_in_A_only = j;
            break;
        }
        if (I->Bmax_exp[j] > I->Bmin_exp[j] && I->Amax_exp[j] == I->Amin_exp[j])
        {
            v_in_B_only = j;
            break;
        }
    }
    if (v_in_A_only != -WORD(1))
    {
        success = _try_missing_var(G, I->Gbits,
                                   v_in_A_only,
                                   A, I->Amin_exp[v_in_A_only],
                                   B, I->Bmin_exp[v_in_A_only],
                                   ctx);
        goto cleanup;
    }
    if (v_in_B_only != -WORD(1))
    {
        success = _try_missing_var(G, I->Gbits,
                                   v_in_B_only,
                                   B, I->Bmin_exp[v_in_B_only],
                                   A, I->Amin_exp[v_in_B_only],
                                   ctx);
        goto cleanup;
    }

    /* all variable are now either
            missing from both ess(A) and ess(B), or
            present in both ess(A) and ess(B)
        and there are at least two in the latter case
    */

    mpoly_gcd_info_set_estimates_nmod_mpoly(I, A, B, ctx, handles, num_handles);
    if (!I->Gdeflate_deg_bounds_are_nice)
        mpoly_gcd_info_set_estimates_nmod_mpoly_lgprime(I, A, B, ctx, handles, num_handles);
    mpoly_gcd_info_set_perm(I, A->length, B->length, ctx->minfo);

    /* everything in I is valid now */

    /* check divisibility A/B and B/A */
    {
        int gcd_is_trivial = 1;
        int try_a = I->Gdeflate_deg_bounds_are_nice;
        int try_b = I->Gdeflate_deg_bounds_are_nice;
        for (j = 0; j < nvars; j++)
        {
            if (I->Gdeflate_deg_bound[j] != 0)
                gcd_is_trivial = 0;

            if (I->Adeflate_deg[j] != I->Gdeflate_deg_bound[j])
                try_a = 0;

            if (I->Bdeflate_deg[j] != I->Gdeflate_deg_bound[j])
                try_b = 0;
        }

        if (gcd_is_trivial)
            goto calculate_trivial_gcd;

        if (try_a && _try_divides(G, B, A, ctx, handles, num_handles))
            goto successful;

        if (try_b && _try_divides(G, A, B, ctx, handles, num_handles))
            goto successful;
    }

    mpoly_gcd_info_measure_brown(I, A->length, B->length, ctx->minfo);

    if (I->mvars == 2)
    {
        if (_try_brown(G, A, B, I, ctx, handles, num_handles))
            goto successful;
    }
    else
    {
        mpoly_gcd_info_measure_zippel(I, A->length, B->length, ctx->minfo);
        mpoly_gcd_info_measure_zippel2(I, A->length, B->length, ctx->minfo);

        if (I->zippel_time_est < I->brown_time_est)
        {
            if (_try_zippel3(G, A, B, I, ctx))
                goto successful;
            if (_try_zippel2(G, A, B, I, ctx))
                goto successful;
            if (_try_zippel(G, A, B, I, ctx))
                goto successful;
            if (_try_brown(G, A, B, I, ctx, handles, num_handles))
                goto successful;
        }
        else
        {
            if (_try_zippel3(G, A, B, I, ctx))
                goto successful;
            if (_try_zippel2(G, A, B, I, ctx))
                goto successful;
            if (_try_zippel(G, A, B, I, ctx))
                goto successful;
            if (_try_brown(G, A, B, I, ctx, handles, num_handles))
                goto successful;

        }
    }

    success = 0;
    goto cleanup;

successful:

    success = 1;

cleanup:

    mpoly_gcd_info_clear(I);

    if (success)
    {
        nmod_mpoly_repack_bits_inplace(G, Gbits, ctx);
        nmod_mpoly_make_monic(G, G, ctx);
        FLINT_ASSERT(nmod_mpoly_divides(T, Asave, G, ctx));
        FLINT_ASSERT(nmod_mpoly_divides(T, Bsave, G, ctx));
    }

#if FLINT_WANT_ASSERT
    nmod_mpoly_clear(T, ctx);
    nmod_mpoly_clear(Asave, ctx);
    nmod_mpoly_clear(Bsave, ctx);
#endif

    return success;
}


int nmod_mpoly_gcd(
    nmod_mpoly_t G,
    const nmod_mpoly_t A,
    const nmod_mpoly_t B,
    const nmod_mpoly_ctx_t ctx)
{
    flint_bitcnt_t Gbits;
    int success;
    thread_pool_handle * handles;
    slong num_handles;
    slong thread_limit;

    thread_limit = FLINT_MIN(A->length, B->length)/256;

    if (nmod_mpoly_is_zero(A, ctx))
    {
        if (nmod_mpoly_is_zero(B, ctx))
            nmod_mpoly_zero(G, ctx);
        else
            nmod_mpoly_make_monic(G, B, ctx);
        return 1;
    }

    if (nmod_mpoly_is_zero(B, ctx))
    {
        nmod_mpoly_make_monic(G, A, ctx);
        return 1;
    }

    Gbits = FLINT_MIN(A->bits, B->bits);

    if (A->bits <= FLINT_BITS && B->bits <= FLINT_BITS)
    {
        num_handles = flint_request_threads(&handles, thread_limit);
        success = _nmod_mpoly_gcd_threaded_pool(G, Gbits, A, B, ctx,
                                                         handles, num_handles);
        flint_give_back_threads(handles, num_handles);
        return success;
    }

    if (A->length == 1)
    {
        return _try_monomial_gcd(G, Gbits, B, A, ctx);
    }
    else if (B->length == 1)
    {
        return _try_monomial_gcd(G, Gbits, A, B, ctx);
    }
    else if (_try_monomial_cofactors(G, Gbits, A, B, ctx))
    {
        return 1;
    }
    else
    {
        /*
            The gcd calculation is unusual.
            First see if both inputs fit into FLINT_BITS.
            Then, try deflation as a last resort.
        */

        int success;
        slong k;
        fmpz * Ashift, * Astride;
        fmpz * Bshift, * Bstride;
        fmpz * Gshift, * Gstride;
        nmod_mpoly_t Anew, Bnew;
        const nmod_mpoly_struct * Ause, * Buse;

        nmod_mpoly_init(Anew, ctx);
        nmod_mpoly_init(Bnew, ctx);

        Ause = A;
        if (A->bits > FLINT_BITS)
        {
            if (!nmod_mpoly_repack_bits(Anew, A, FLINT_BITS, ctx))
                goto could_not_repack;
            Ause = Anew;
        }

        Buse = B;
        if (B->bits > FLINT_BITS)
        {
            if (!nmod_mpoly_repack_bits(Bnew, B, FLINT_BITS, ctx))
                goto could_not_repack;
            Buse = Bnew;
        }

        num_handles = flint_request_threads(&handles, thread_limit);
        Gbits = FLINT_MIN(Ause->bits, Buse->bits);
        success = _nmod_mpoly_gcd_threaded_pool(G, Gbits, Ause, Buse, ctx,
                                                         handles, num_handles);
        flint_give_back_threads(handles, num_handles);

        goto cleanup;

could_not_repack:

        /*
            One of A or B could not be repacked into FLINT_BITS. See if
            they both fit into FLINT_BITS after deflation.
        */

        Ashift  = _fmpz_vec_init(ctx->minfo->nvars);
        Astride = _fmpz_vec_init(ctx->minfo->nvars);
        Bshift  = _fmpz_vec_init(ctx->minfo->nvars);
        Bstride = _fmpz_vec_init(ctx->minfo->nvars);
        Gshift  = _fmpz_vec_init(ctx->minfo->nvars);
        Gstride = _fmpz_vec_init(ctx->minfo->nvars);

        nmod_mpoly_deflation(Ashift, Astride, A, ctx);
        nmod_mpoly_deflation(Bshift, Bstride, B, ctx);
        _fmpz_vec_min(Gshift, Ashift, Bshift, ctx->minfo->nvars);
        for (k = 0; k < ctx->minfo->nvars; k++)
        {
            fmpz_gcd(Gstride + k, Astride + k, Bstride + k);
        }

        success = 0;

        nmod_mpoly_deflate(Anew, A, Ashift, Gstride, ctx);
        if (Anew->bits > FLINT_BITS)
        {
            if (!nmod_mpoly_repack_bits(Anew, Anew, FLINT_BITS, ctx))
                goto deflate_cleanup;
        }

        nmod_mpoly_deflate(Bnew, B, Bshift, Gstride, ctx);
        if (Bnew->bits > FLINT_BITS)
        {
            if (!nmod_mpoly_repack_bits(Bnew, Bnew, FLINT_BITS, ctx))
                goto deflate_cleanup;
        }

        num_handles = flint_request_threads(&handles, thread_limit);
        Gbits = FLINT_MIN(Anew->bits, Bnew->bits);
        success = _nmod_mpoly_gcd_threaded_pool(G, Gbits, Anew, Bnew, ctx,
                                                         handles, num_handles);
        flint_give_back_threads(handles, num_handles);

        if (success)
        {
            nmod_mpoly_inflate(G, G, Gshift, Gstride, ctx);
            nmod_mpoly_make_monic(G, G, ctx);
        }

deflate_cleanup:

        _fmpz_vec_clear(Ashift, ctx->minfo->nvars);
        _fmpz_vec_clear(Astride, ctx->minfo->nvars);
        _fmpz_vec_clear(Bshift, ctx->minfo->nvars);
        _fmpz_vec_clear(Bstride, ctx->minfo->nvars);
        _fmpz_vec_clear(Gshift, ctx->minfo->nvars);
        _fmpz_vec_clear(Gstride, ctx->minfo->nvars);

cleanup:

        nmod_mpoly_clear(Anew, ctx);
        nmod_mpoly_clear(Bnew, ctx);

        return success;
    }
}

