/*
    Copyright (C) 2019-2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fq_nmod_mpoly.h"
#include "fq_nmod_mpoly_factor.h"

/*
    For each j, set out[j] to the evaluation of A at x_i = alpha[i] (i != j)
    i.e. if nvars = 3
        out[0] = A(x, alpha[1], alpha[2])
        out[1] = A(alpha[0], x, alpha[2])
        out[2] = A(alpha[0], alpha[1], x)

    If ignore[j] is nonzero, then out[j] need not be calculated, probably
    because we shouldn't calculate it in dense form.
*/
void fq_nmod_mpoly_evals(
    n_fq_poly_struct * out,
    const int * ignore,
    const fq_nmod_mpoly_t A,
    ulong * Amin_exp,
    ulong * Amax_exp,
    ulong * Astride,
    fq_nmod_struct * alpha,
    const fq_nmod_mpoly_ctx_t ctx)
{
    slong d = fq_nmod_ctx_degree(ctx->fqctx);
    slong i, j;
    slong nvars = ctx->minfo->nvars;
    ulong mask = (-UWORD(1)) >> (FLINT_BITS - A->bits);
    slong N = mpoly_words_per_exp_sp(A->bits, ctx->minfo);
    slong * offsets = FLINT_ARRAY_ALLOC(2*nvars, slong);
    slong * shifts = offsets + nvars;
    ulong * varexps = FLINT_ARRAY_ALLOC(nvars, ulong);
    n_poly_struct * caches = FLINT_ARRAY_ALLOC(3*nvars, n_poly_struct);
    mp_limb_t * t = FLINT_ARRAY_ALLOC(2*d, mp_limb_t);
    mp_limb_t * meval = t + d;

    for (j = 0; j < nvars; j++)
    {
        mpoly_gen_offset_shift_sp(offsets + j, shifts + j, j, A->bits,
                                                                   ctx->minfo);
        n_poly_init(caches + 3*j + 0);
        n_poly_init(caches + 3*j + 1);
        n_poly_init(caches + 3*j + 2);
        n_fq_pow_cache_start_fq_nmod(alpha + j, caches + 3*j + 0,
                               caches + 3*j + 1, caches + 3*j + 2, ctx->fqctx);

        if (ignore[j])
            continue;

        i = (Astride[j] < 2) ? Amax_exp[j] - Amin_exp[j] :
                            (Amax_exp[j] - Amin_exp[j])/Astride[j];

        n_poly_fit_length(out + j, d*(i + 1));
        _nmod_vec_zero(out[j].coeffs, d*(i + 1));
        out[j].length = i + 1;
    }

    for (i = 0; i < A->length; i++)
    {
        ulong varexp;
        mp_limb_t * s = A->coeffs + d*i; /* source */

        for (j = 0; j < nvars; j++)
        {
            varexp = ((A->exps + N*i)[offsets[j]]>>shifts[j])&mask;

            FLINT_ASSERT((Astride[j] == 0 && varexp == Amin_exp[j]) ||
                                     (varexp - Amin_exp[j]) % Astride[j] == 0);

            varexps[j] = Astride[j] < 2 ? varexp - Amin_exp[j] :
                                         (varexp - Amin_exp[j])/Astride[j];

            n_fq_pow_cache_mulpow_ui(meval, s, varexps[j], caches + 3*j + 0,
                               caches + 3*j + 1, caches + 3*j + 2, ctx->fqctx);
            s = meval;
        }

        for (j = 0; j < nvars; j++)
        {
            varexp = varexps[j];

            if (ignore[j])
                continue;

            FLINT_ASSERT(out[j].alloc >= d*(varexp + 1));

            n_fq_pow_cache_mulpow_neg_ui(t, meval, varexp, caches + 3*j + 0,
                               caches + 3*j + 1, caches + 3*j + 2, ctx->fqctx);

            n_fq_add(out[j].coeffs + d*varexp, out[j].coeffs + d*varexp,
                                                                t, ctx->fqctx);
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

void fq_nmod_mpoly_evals_lgprime(
    n_fq_poly_struct * out,
    const int * ignore,
    const fq_nmod_mpoly_t A,
    ulong * Amin_exp,
    ulong * Amax_exp,
    ulong * Astride,
    const fq_nmod_mpoly_ctx_t smctx,
    fq_nmod_struct * alpha,
    const fq_nmod_mpoly_ctx_t lgctx,
    const bad_fq_nmod_embed_t emb)
{
    slong smd = fq_nmod_ctx_degree(smctx->fqctx);
    slong lgd = fq_nmod_ctx_degree(lgctx->fqctx);
    slong i, j;
    slong nvars = smctx->minfo->nvars;
    ulong mask = (-UWORD(1)) >> (FLINT_BITS - A->bits);
    slong N = mpoly_words_per_exp_sp(A->bits, smctx->minfo);
    slong * offsets = FLINT_ARRAY_ALLOC(2*nvars, slong);
    slong * shifts = offsets + nvars;
    ulong * varexps = FLINT_ARRAY_ALLOC(nvars, ulong);
    n_poly_struct * caches = FLINT_ARRAY_ALLOC(3*nvars, n_poly_struct);
    mp_limb_t * t = FLINT_ARRAY_ALLOC(2*lgd, mp_limb_t);
    mp_limb_t * meval = t + lgd;

    for (j = 0; j < nvars; j++)
    {
        mpoly_gen_offset_shift_sp(offsets + j, shifts + j, j, A->bits,
                                                                 smctx->minfo);
        n_poly_init(caches + 3*j + 0);
        n_poly_init(caches + 3*j + 1);
        n_poly_init(caches + 3*j + 2);
        n_fq_pow_cache_start_fq_nmod(alpha + j, caches + 3*j + 0,
                             caches + 3*j + 1, caches + 3*j + 2, lgctx->fqctx);
        if (ignore[j])
            continue;

        i = (Astride[j] < 2) ? Amax_exp[j] - Amin_exp[j] :
                              (Amax_exp[j] - Amin_exp[j])/Astride[j];

        n_poly_fit_length(out + j, lgd*(i + 1));
        _nmod_vec_zero(out[j].coeffs, lgd*(i + 1));
        out[j].length = i + 1;
    }

    for (i = 0; i < A->length; i++)
    {
        ulong varexp;

        bad_n_fq_embed_sm_elem_to_lg(meval, A->coeffs + smd*i, emb);

        for (j = 0; j < nvars; j++)
        {
            varexp = ((A->exps + N*i)[offsets[j]]>>shifts[j])&mask;

            FLINT_ASSERT((Astride[j] == 0 && varexp == Amin_exp[j]) ||
                                     (varexp - Amin_exp[j]) % Astride[j] == 0);

            varexps[j] = Astride[j] < 2 ? varexp - Amin_exp[j] :
                                         (varexp - Amin_exp[j])/Astride[j];

            n_fq_pow_cache_mulpow_ui(meval, meval, varexps[j], caches + 3*j + 0,
                             caches + 3*j + 1, caches + 3*j + 2, lgctx->fqctx);
        }

        for (j = 0; j < nvars; j++)
        {
            varexp = varexps[j];

            if (ignore[j])
                continue;

            FLINT_ASSERT(out[j].alloc >= lgd*(varexp + 1));

            n_fq_pow_cache_mulpow_neg_ui(t, meval, varexp, caches + 3*j + 0,
                             caches + 3*j + 1, caches + 3*j + 2, lgctx->fqctx);

            n_fq_add(out[j].coeffs + lgd*varexp, out[j].coeffs + lgd*varexp,
                                                              t, lgctx->fqctx);
        }
    }

    for (j = 0; j < nvars; j++)
        _n_fq_poly_normalise(out + j, lgd);

    for (j = 0; j < 3*nvars; j++)
        n_poly_clear(caches + j);

    flint_free(offsets);
    flint_free(varexps);
    flint_free(caches);
    flint_free(t);
}


void mpoly_gcd_info_set_estimates_fq_nmod_mpoly(
    mpoly_gcd_info_t I,
    const fq_nmod_mpoly_t A,
    const fq_nmod_mpoly_t B,
    const fq_nmod_mpoly_ctx_t ctx)
{
    slong d = fq_nmod_ctx_degree(ctx->fqctx);
    int tries_left = 10;
    slong nvars = ctx->minfo->nvars;
    slong i, j;
    n_fq_poly_t Geval;
    n_fq_poly_struct * Aevals, * Bevals;
    fq_nmod_struct * alpha;
    flint_rand_t state;
    slong ignore_limit;
    int * ignore;

    flint_randinit(state);

    ignore = FLINT_ARRAY_ALLOC(nvars, int);
    alpha  = FLINT_ARRAY_ALLOC(nvars, fq_nmod_struct);
    Aevals = FLINT_ARRAY_ALLOC(2*nvars, n_fq_poly_struct);
    Bevals = Aevals + nvars;

    n_fq_poly_init(Geval);
    for (j = 0; j < nvars; j++)
    {
        fq_nmod_init(alpha + j, ctx->fqctx);
        n_fq_poly_init(Aevals + j);
        n_fq_poly_init(Bevals + j);
    }

    ignore_limit = FLINT_MAX(WORD(9999), (A->length + B->length)/4096);
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

    if (--tries_left < 0)
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
        fq_nmod_rand(alpha + j, state, ctx->fqctx);
        if (fq_nmod_is_zero(alpha + j, ctx->fqctx))
            fq_nmod_one(alpha + j, ctx->fqctx);
    }

    fq_nmod_mpoly_evals(Aevals, ignore, A, I->Amin_exp, I->Amax_exp, I->Gstride, alpha, ctx);
    fq_nmod_mpoly_evals(Bevals, ignore, B, I->Bmin_exp, I->Bmax_exp, I->Gstride, alpha, ctx);

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
                goto try_again;
            }

            n_fq_poly_gcd(Geval, Aevals + j, Bevals + j, ctx->fqctx);

            I->Gterm_count_est[j] = 0;
            I->Gdeflate_deg_bound[j] = n_fq_poly_degree(Geval);
            for (i = I->Gdeflate_deg_bound[j]; i >= 0; i--)
                I->Gterm_count_est[j] += _n_fq_is_zero(Geval->coeffs + d*i, d);
        }
    }

cleanup:

    n_fq_poly_clear(Geval);
    for (j = 0; j < nvars; j++)
    {
        fq_nmod_clear(alpha + j, ctx->fqctx);
        n_fq_poly_clear(Aevals + j);
        n_fq_poly_clear(Bevals + j);
    }

    flint_free(ignore);
    flint_free(alpha);
    flint_free(Aevals);

    flint_randclear(state);

    return;
}


void mpoly_gcd_info_set_estimates_fq_nmod_mpoly_lgprime(
    mpoly_gcd_info_t I,
    const fq_nmod_mpoly_t A,
    const fq_nmod_mpoly_t B,
    const fq_nmod_mpoly_ctx_t smctx)
{
    int tries_left = 10;
    slong nvars = smctx->minfo->nvars;
    slong i, j;
    n_fq_poly_t Geval;
    n_fq_poly_struct * Aevals, * Bevals;
    fq_nmod_struct * alpha;
    flint_rand_t state;
    slong ignore_limit;
    int * ignore;
    fq_nmod_mpoly_ctx_t lgctx;
    bad_fq_nmod_mpoly_embed_chooser_t embc;
    bad_fq_nmod_embed_struct * cur_emb;

    flint_randinit(state);

    cur_emb = bad_fq_nmod_mpoly_embed_chooser_init(embc, lgctx, smctx, state);

    ignore = FLINT_ARRAY_ALLOC(nvars, int);
    alpha  = FLINT_ARRAY_ALLOC(nvars, fq_nmod_struct);
    Aevals = FLINT_ARRAY_ALLOC(2*nvars, n_fq_poly_struct);
    Bevals = Aevals + nvars;

    n_fq_poly_init(Geval);
    for (j = 0; j < nvars; j++)
    {
        fq_nmod_init(alpha + j, lgctx->fqctx);
        n_fq_poly_init(Aevals + j);
        n_fq_poly_init(Bevals + j);
    }

    ignore_limit = FLINT_MAX(WORD(9999), (A->length + B->length)/4096);
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

    if (--tries_left < 0 || cur_emb == NULL)
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

    fq_nmod_mpoly_evals_lgprime(Aevals, ignore, A, I->Amin_exp, I->Amax_exp,
                                     I->Gstride, smctx, alpha, lgctx, cur_emb);
    fq_nmod_mpoly_evals_lgprime(Bevals, ignore, B, I->Bmin_exp, I->Bmax_exp,
                                     I->Gstride, smctx, alpha, lgctx, cur_emb);

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
            slong lgd = fq_nmod_ctx_degree(lgctx->fqctx);

            if (I->Adeflate_deg[j] != n_fq_poly_degree(Aevals + j) ||
                I->Bdeflate_deg[j] != n_fq_poly_degree(Bevals + j))
            {
                cur_emb = bad_fq_nmod_mpoly_embed_chooser_next(embc, lgctx, smctx, state);
                goto try_again;
            }

            n_fq_poly_gcd(Geval, Aevals + j, Bevals + j, lgctx->fqctx);

            I->Gterm_count_est[j] = 0;
            I->Gdeflate_deg_bound[j] = n_fq_poly_degree(Geval);
            for (i = I->Gdeflate_deg_bound[j]; i >= 0; i--)
                I->Gterm_count_est[j] += _n_fq_is_zero(Geval->coeffs + lgd*i, lgd);
        }
    }

cleanup:

    n_fq_poly_clear(Geval);
    for (j = 0; j < nvars; j++)
    {
        fq_nmod_clear(alpha + j, lgctx->fqctx);
        n_fq_poly_clear(Aevals + j);
        n_fq_poly_clear(Bevals + j);
    }

    flint_free(ignore);
    flint_free(alpha);
    flint_free(Aevals);

    bad_fq_nmod_mpoly_embed_chooser_clear(embc, lgctx, smctx, state);

    flint_randclear(state);

    return;
}


/* (Abar, Bbar) = (A, B) */
static void _parallel_set(
    fq_nmod_mpoly_t Abar, /* could be NULL */
    fq_nmod_mpoly_t Bbar, /* could be NULL */
    const fq_nmod_mpoly_t A,
    const fq_nmod_mpoly_t B,
    const fq_nmod_mpoly_ctx_t ctx)
{
    if (Abar == B && Bbar == A)
    {
        FLINT_ASSERT(Abar != NULL && Bbar != NULL);
        fq_nmod_mpoly_set(Abar, B, ctx);
        fq_nmod_mpoly_set(Bbar, A, ctx);
        fq_nmod_mpoly_swap(Abar, Bbar, ctx);
    }
    else if (Abar == B && Bbar != A)
    {
        FLINT_ASSERT(Abar != NULL);
        if (Bbar != NULL)
            fq_nmod_mpoly_set(Bbar, B, ctx);
        fq_nmod_mpoly_set(Abar, A, ctx);
    }
    else
    {
        if (Abar != NULL)
            fq_nmod_mpoly_set(Abar, A, ctx);
        if (Bbar != NULL)
            fq_nmod_mpoly_set(Bbar, B, ctx);
    }
}


/* The variables in ess(A) and ess(B) are disjoint. gcd is trivial to compute */
static int _do_trivial(
    fq_nmod_mpoly_t G,
    fq_nmod_mpoly_t Abar,  /* could be NULL */
    fq_nmod_mpoly_t Bbar,  /* could be NULL */
    const fq_nmod_mpoly_t A,
    const fq_nmod_mpoly_t B,
    const mpoly_gcd_info_t I,
    const fq_nmod_mpoly_ctx_t ctx)
{
    _parallel_set(Abar, Bbar, A, B, ctx);

    if (Abar != NULL)
        mpoly_monomials_shift_right_ui(Abar->exps, Abar->bits, Abar->length,
                                                      I->Gmin_exp, ctx->minfo);

    if (Bbar != NULL)
        mpoly_monomials_shift_right_ui(Bbar->exps, Bbar->bits, Bbar->length,
                                                      I->Gmin_exp, ctx->minfo);

    fq_nmod_mpoly_fit_length_reset_bits(G, 1, I->Gbits, ctx);
    mpoly_set_monomial_ui(G->exps, I->Gmin_exp, I->Gbits, ctx->minfo);
    _n_fq_one(G->coeffs, fq_nmod_ctx_degree(ctx->fqctx));
    _fq_nmod_mpoly_set_length(G, 1, ctx);

    return 1;
}

/*********************** Easy when B is a monomial ***************************/
static int _do_monomial_gcd(
    fq_nmod_mpoly_t G,
    fq_nmod_mpoly_t Abar,  /* could be NULL */
    fq_nmod_mpoly_t Bbar,  /* could be NULL */
    const fq_nmod_mpoly_t A,
    const fq_nmod_mpoly_t B,
    const fq_nmod_mpoly_ctx_t ctx)
{
    slong i;
    flint_bitcnt_t Gbits = FLINT_MIN(A->bits, B->bits);
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

    _parallel_set(Abar, Bbar, A, B, ctx);

    if (Abar != NULL)
        mpoly_monomials_shift_right_ffmpz(Abar->exps, Abar->bits, Abar->length,
                                                         minBdegs, ctx->minfo);

    if (Bbar != NULL)
        mpoly_monomials_shift_right_ffmpz(Bbar->exps, Bbar->bits, Bbar->length,
                                                         minBdegs, ctx->minfo);

    fq_nmod_mpoly_fit_length_reset_bits(G, 1, Gbits, ctx);
    mpoly_set_monomial_ffmpz(G->exps, minBdegs, Gbits, ctx->minfo);
    _n_fq_one(G->coeffs + 0, fq_nmod_ctx_degree(ctx->fqctx));
    _fq_nmod_mpoly_set_length(G, 1, ctx);

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
    fq_nmod_mpoly_t G,
    fq_nmod_mpoly_t Abar,  /* could be NULL */
    fq_nmod_mpoly_t Bbar,  /* could be NULL */
    const fq_nmod_mpoly_t A,
    const fq_nmod_mpoly_t B,
    const fq_nmod_mpoly_ctx_t ctx)
{
    slong d = fq_nmod_ctx_degree(ctx->fqctx);
    int success;
    slong i, j;
    slong NA, NG;
    slong nvars = ctx->minfo->nvars;
    fmpz * Abarexps, * Bbarexps, * Texps;
    mp_limb_t * tmp, * t1, * t2, * a0, * b0;
    fq_nmod_mpoly_t T;
    flint_bitcnt_t Gbits = FLINT_MIN(A->bits, B->bits);
    flint_bitcnt_t Abarbits = A->bits;
    flint_bitcnt_t Bbarbits = B->bits;
    TMP_INIT;

    FLINT_ASSERT(A->length > 0);
    FLINT_ASSERT(B->length > 0);

    if (A->length != B->length)
        return 0;

    TMP_START;

    tmp = (mp_limb_t *) TMP_ALLOC(d*(4 + FLINT_MAX(N_FQ_MUL_ITCH,
                                            N_FQ_INV_ITCH))*sizeof(mp_limb_t));
    t1 = tmp + d*FLINT_MAX(N_FQ_MUL_ITCH, N_FQ_INV_ITCH);
    t2 = t1 + d;
    a0 = t2 + d;
    b0 = a0 + d;

    for (i = A->length - 1; i > 0; i--)
    {
        _n_fq_mul(t1, A->coeffs + d*0, B->coeffs + d*i, ctx->fqctx, tmp);
        _n_fq_mul(t2, B->coeffs + d*0, A->coeffs + d*i, ctx->fqctx, tmp);
        success = _n_fq_equal(t1, t2, d);
        if (!success)
            goto cleanup_less;
    }

    _n_fq_set(a0, A->coeffs + d*0, d);
    _n_fq_set(b0, B->coeffs + d*0, d);

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
        goto cleanup_more;

    fq_nmod_mpoly_init3(T, A->length, Gbits, ctx);
    NG = mpoly_words_per_exp(Gbits, ctx->minfo);
    NA = mpoly_words_per_exp(A->bits, ctx->minfo);
    _n_fq_inv(t1, A->coeffs + d*0, ctx->fqctx, tmp);
    T->length = A->length;
    for (i = 0; i < A->length; i++)
    {
        mpoly_get_monomial_ffmpz(Texps, A->exps + NA*i, A->bits, ctx->minfo);
        _fmpz_vec_sub(Texps, Texps, Abarexps, nvars);
        mpoly_set_monomial_ffmpz(T->exps + NG*i, Texps, Gbits, ctx->minfo);
        n_fq_mul(T->coeffs + d*i, A->coeffs + d*i, t1, ctx->fqctx);
    }
    fq_nmod_mpoly_swap(G, T, ctx);
    fq_nmod_mpoly_clear(T, ctx);

    if (Abar != NULL)
    {
        fq_nmod_mpoly_fit_length_reset_bits(Abar, 1, Abarbits, ctx);
        mpoly_set_monomial_ffmpz(Abar->exps, Abarexps, Abarbits, ctx->minfo);
        _n_fq_set(Abar->coeffs + d*0, a0, d);
        _fq_nmod_mpoly_set_length(Abar, 1, ctx);
    }

    if (Bbar != NULL)
    {
        fq_nmod_mpoly_fit_length_reset_bits(Bbar, 1, Bbarbits, ctx);
        mpoly_set_monomial_ffmpz(Bbar->exps, Bbarexps, Bbarbits, ctx->minfo);
        _n_fq_set(Bbar->coeffs + d*0, b0, d);
        _fq_nmod_mpoly_set_length(Bbar, 1, ctx);
    }

    success = 1;

cleanup_more:

    for (j = 0; j < nvars; j++)
    {
        fmpz_clear(Abarexps + j);
        fmpz_clear(Bbarexps + j);
        fmpz_clear(Texps + j);
    }

cleanup_less:

    TMP_END;

    return success;
}


/***  ess(A) and ess(B) depend on only one variable v_in_both ****************/
static int _do_univar(
    fq_nmod_mpoly_t G,
    fq_nmod_mpoly_t Abar,
    fq_nmod_mpoly_t Bbar,
    const fq_nmod_mpoly_t A,
    const fq_nmod_mpoly_t B,
    slong v_in_both,
    const mpoly_gcd_info_t I,
    const fq_nmod_mpoly_ctx_t ctx)
{
    fq_nmod_poly_t a, b, g, t, r;

    fq_nmod_poly_init(a, ctx->fqctx);
    fq_nmod_poly_init(b, ctx->fqctx);
    fq_nmod_poly_init(g, ctx->fqctx);
    fq_nmod_poly_init(t, ctx->fqctx);
    fq_nmod_poly_init(r, ctx->fqctx);

    _fq_nmod_mpoly_to_fq_nmod_poly_deflate(a, A, v_in_both,
                                                 I->Amin_exp, I->Gstride, ctx);
    _fq_nmod_mpoly_to_fq_nmod_poly_deflate(b, B, v_in_both,
                                                 I->Bmin_exp, I->Gstride, ctx);

    fq_nmod_poly_gcd(g, a, b, ctx->fqctx);
    _fq_nmod_mpoly_from_fq_nmod_poly_inflate(G, I->Gbits, g, v_in_both,
                                                 I->Gmin_exp, I->Gstride, ctx);
    if (Abar != NULL)
    {
        fq_nmod_poly_divrem(t, r, a, g, ctx->fqctx);
        FLINT_ASSERT(fq_nmod_poly_is_zero(r, ctx->fqctx));
        _fq_nmod_mpoly_from_fq_nmod_poly_inflate(Abar, I->Abarbits, t,
                                   v_in_both, I->Abarmin_exp, I->Gstride, ctx);
    }

    if (Bbar != NULL)
    {
        fq_nmod_poly_divrem(t, r, b, g, ctx->fqctx);
        FLINT_ASSERT(fq_nmod_poly_is_zero(r, ctx->fqctx));
        _fq_nmod_mpoly_from_fq_nmod_poly_inflate(Bbar, I->Bbarbits, t, v_in_both,
                                              I->Bbarmin_exp, I->Gstride, ctx);
    }

    fq_nmod_poly_clear(a, ctx->fqctx);
    fq_nmod_poly_clear(b, ctx->fqctx);
    fq_nmod_poly_clear(g, ctx->fqctx);
    fq_nmod_poly_clear(t, ctx->fqctx);
    fq_nmod_poly_clear(r, ctx->fqctx);

    return 1;
}

/********* Assume B has length one when converted to univar format ***********/
static int _try_missing_var(
    fq_nmod_mpoly_t G, flint_bitcnt_t Gbits,
    fq_nmod_mpoly_t Abar,
    fq_nmod_mpoly_t Bbar,
    slong var,
    const fq_nmod_mpoly_t A, ulong Ashift,
    const fq_nmod_mpoly_t B, ulong Bshift,
    const fq_nmod_mpoly_ctx_t ctx)
{
    int success;
    fq_nmod_mpoly_univar_t Au;

    fq_nmod_mpoly_univar_init(Au, ctx);

#if FLINT_WANT_ASSERT
    fq_nmod_mpoly_to_univar(Au, B, var, ctx);
    FLINT_ASSERT(Au->length == 1);
#endif
    fq_nmod_mpoly_to_univar(Au, A, var, ctx);

    fq_nmod_mpoly_univar_fit_length(Au, Au->length + 1, ctx);
    fq_nmod_mpoly_set(Au->coeffs + Au->length, B, ctx);
    Au->length++;

    if (Abar == NULL && Bbar == NULL)
    {
        success = _fq_nmod_mpoly_vec_content_mpoly(G, Au->coeffs, Au->length, ctx);
        if (!success)
            goto cleanup;

        fq_nmod_mpoly_repack_bits_inplace(G, Gbits, ctx);
        _mpoly_gen_shift_left(G->exps, G->bits, G->length,
                                   var, FLINT_MIN(Ashift, Bshift), ctx->minfo);
    }
    else
    {
        fq_nmod_mpoly_t tG, tAbar, tBbar;

        fq_nmod_mpoly_init(tG, ctx);
        fq_nmod_mpoly_init(tAbar, ctx);
        fq_nmod_mpoly_init(tBbar, ctx);

        success = _fq_nmod_mpoly_vec_content_mpoly(G, Au->coeffs, Au->length, ctx);
        if (!success)
            goto cleanup;

        fq_nmod_mpoly_repack_bits_inplace(tG, Gbits, ctx);
        _mpoly_gen_shift_left(tG->exps, tG->bits, tG->length,
                                   var, FLINT_MIN(Ashift, Bshift), ctx->minfo);

        if (Abar != NULL)
        {
            success = fq_nmod_mpoly_divides(tAbar, A, tG, ctx);
            FLINT_ASSERT(success);
        }

        if (Bbar != NULL)
        {
            success = fq_nmod_mpoly_divides(tBbar, B, tG, ctx);
            FLINT_ASSERT(success);
        }

        fq_nmod_mpoly_swap(G, tG, ctx);

        if (Abar != NULL)
            fq_nmod_mpoly_swap(Abar, tAbar, ctx);

        if (Bbar != NULL)
            fq_nmod_mpoly_swap(Bbar, tBbar, ctx);

        fq_nmod_mpoly_clear(tG, ctx);
        fq_nmod_mpoly_clear(tAbar, ctx);
        fq_nmod_mpoly_clear(tBbar, ctx);        
    }

    success = 1;

cleanup:

    fq_nmod_mpoly_univar_clear(Au, ctx);

    return success;
}


/************************ See if B divides A ********************************/
static int _try_divides(
    fq_nmod_mpoly_t G,
    fq_nmod_mpoly_t Abar,
    fq_nmod_mpoly_t Bbar,
    const fq_nmod_mpoly_t A,
    const fq_nmod_mpoly_t BB,
    const fq_nmod_mpoly_ctx_t ctx)
{
    int success = 0;
    fq_nmod_mpoly_t Q, B, M;

    fq_nmod_mpoly_init(Q, ctx);
    fq_nmod_mpoly_init(B, ctx);
    fq_nmod_mpoly_init(M, ctx);

    /* BB = M*B */
    fq_nmod_mpoly_term_content(M, BB, ctx);
    fq_nmod_mpoly_divides(B, BB, M, ctx);

    if (fq_nmod_mpoly_divides(Q, A, B, ctx))
    {
        /* gcd(Q*B, M*B) */
        _do_monomial_gcd(G, Abar, Bbar, Q, M, ctx);
        fq_nmod_mpoly_mul(G, G, B, ctx);
        success = 1;
    }

    fq_nmod_mpoly_clear(Q, ctx);
    fq_nmod_mpoly_clear(B, ctx);
    fq_nmod_mpoly_clear(M, ctx);

    return success;
}

/********************** Hit A and B with zippel ******************************/
static int _try_zippel(
    fq_nmod_mpoly_t G,
    fq_nmod_mpoly_t Abar,
    fq_nmod_mpoly_t Bbar,
    const fq_nmod_mpoly_t A,
    const fq_nmod_mpoly_t B,
    const mpoly_gcd_info_t I,
    const fq_nmod_mpoly_ctx_t ctx)
{
    slong i, k;
    slong m = I->mvars;
    int success;
    mpoly_zipinfo_t zinfo;
    flint_bitcnt_t wbits;
    flint_rand_t randstate;
    fq_nmod_mpoly_ctx_t uctx;
    fq_nmod_mpolyu_t Au, Bu, Gu, Abaru, Bbaru;
    fq_nmod_mpoly_t Ac, Bc, Gc, Abarc, Bbarc;

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

    /* uctx is context for Fq[y_1,...,y_{m-1}]*/
    fq_nmod_mpoly_ctx_init(uctx, m - 1, ORD_LEX, ctx->fqctx);

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

    fq_nmod_mpolyu_init(Au, wbits, uctx);
    fq_nmod_mpolyu_init(Bu, wbits, uctx);
    fq_nmod_mpolyu_init(Gu, wbits, uctx);
    fq_nmod_mpolyu_init(Abaru, wbits, uctx);
    fq_nmod_mpolyu_init(Bbaru, wbits, uctx);
    fq_nmod_mpoly_init3(Ac, 0, wbits, uctx);
    fq_nmod_mpoly_init3(Bc, 0, wbits, uctx);
    fq_nmod_mpoly_init3(Gc, 0, wbits, uctx);
    fq_nmod_mpoly_init3(Abarc, 0, wbits, uctx);
    fq_nmod_mpoly_init3(Bbarc, 0, wbits, uctx);

    fq_nmod_mpoly_to_mpolyu_perm_deflate(Au, uctx, A, ctx, zinfo->perm,
                                                      I->Amin_exp, I->Gstride);
    fq_nmod_mpoly_to_mpolyu_perm_deflate(Bu, uctx, B, ctx, zinfo->perm,
                                                      I->Bmin_exp, I->Gstride);

    success = fq_nmod_mpolyu_content_mpoly(Ac, Au, uctx) &&
              fq_nmod_mpolyu_content_mpoly(Bc, Bu, uctx);
    if (!success)
        goto cleanup;

    fq_nmod_mpolyu_divexact_mpoly_inplace(Au, Ac, uctx);
    fq_nmod_mpolyu_divexact_mpoly_inplace(Bu, Bc, uctx);

    /* after removing content, degree bounds in zinfo are still valid bounds */
    success = fq_nmod_mpolyu_gcdm_zippel(Gu, Abaru, Bbaru, Au, Bu,
                                                       uctx, zinfo, randstate);
    if (!success)
        goto cleanup;

    if (Abar == NULL && Bbar == NULL)
    {
        success = fq_nmod_mpoly_gcd(Gc, Ac, Bc, uctx);
        if (!success)
            goto cleanup;

        fq_nmod_mpoly_repack_bits_inplace(Gc, wbits, uctx);
        fq_nmod_mpolyu_mul_mpoly_inplace(Gu, Gc, uctx);

        fq_nmod_mpoly_from_mpolyu_perm_inflate(G, I->Gbits, ctx, Gu, uctx,
                                         zinfo->perm, I->Gmin_exp, I->Gstride);
    }
    else
    {
        success = fq_nmod_mpoly_gcd_cofactors(Gc, Abarc, Bbarc, Ac, Bc, uctx);
        if (!success)
            goto cleanup;

        fq_nmod_mpoly_repack_bits_inplace(Gc, wbits, uctx);
        fq_nmod_mpoly_repack_bits_inplace(Abarc, wbits, uctx);
        fq_nmod_mpoly_repack_bits_inplace(Bbarc, wbits, uctx);

        fq_nmod_mpolyu_mul_mpoly_inplace(Gu, Gc, uctx);
        fq_nmod_mpolyu_mul_mpoly_inplace(Abaru, Abarc, uctx);
        fq_nmod_mpolyu_mul_mpoly_inplace(Bbaru, Bbarc, uctx);

        fq_nmod_mpoly_from_mpolyu_perm_inflate(G, I->Gbits, ctx, Gu, uctx,
                                        zinfo->perm, I->Gmin_exp, I->Gstride);

        if (Abar != NULL)
            fq_nmod_mpoly_from_mpolyu_perm_inflate(Abar, I->Abarbits, ctx,
                         Abaru, uctx, zinfo->perm, I->Abarmin_exp, I->Gstride);

        if (Bbar != NULL)
            fq_nmod_mpoly_from_mpolyu_perm_inflate(Bbar, I->Bbarbits, ctx,
                         Bbaru, uctx, zinfo->perm, I->Bbarmin_exp, I->Gstride);
    }

    success = 1;

cleanup:

    fq_nmod_mpolyu_clear(Au, uctx);
    fq_nmod_mpolyu_clear(Bu, uctx);
    fq_nmod_mpolyu_clear(Gu, uctx);
    fq_nmod_mpolyu_clear(Abaru, uctx);
    fq_nmod_mpolyu_clear(Bbaru, uctx);
    fq_nmod_mpoly_clear(Ac, uctx);
    fq_nmod_mpoly_clear(Bc, uctx);
    fq_nmod_mpoly_clear(Gc, uctx);
    fq_nmod_mpoly_clear(Abarc, uctx);
    fq_nmod_mpoly_clear(Bbarc, uctx);

    fq_nmod_mpoly_ctx_clear(uctx);

    mpoly_zipinfo_clear(zinfo);

    flint_randclear(randstate);

    return success;
}

static int _try_zippel2_old(
    fq_nmod_mpoly_t G,
    const fq_nmod_mpoly_t A,
    const fq_nmod_mpoly_t B,
    const mpoly_gcd_info_t I,
    const fq_nmod_mpoly_ctx_t ctx)
{
    slong i, k;
    slong m = I->mvars;
    int success;
    flint_bitcnt_t wbits;
    fq_nmod_mpoly_ctx_t uctx;
    fq_nmod_mpolyu_t Auu, Buu, Guu, Abaruu, Bbaruu;
    fq_nmod_mpoly_t Ac, Bc, Gc, Gamma;
    slong * tmp, * Guu_degs, * Auu_degs, * Buu_degs, * Gamma_degs, * Gguess;
    slong max_minor_degree;

    FLINT_ASSERT(A->bits <= FLINT_BITS);
    FLINT_ASSERT(B->bits <= FLINT_BITS);
    FLINT_ASSERT(A->length > 0);
    FLINT_ASSERT(B->length > 0);

    if (!I->can_use_bma)
        return 0;

flint_printf("_try_zippel2\n");

    FLINT_ASSERT(m >= WORD(3));

    tmp = FLINT_ARRAY_ALLOC(5*m, slong);
    Auu_degs   = tmp + 1*m;
    Buu_degs   = tmp + 2*m;
    Gamma_degs = tmp + 3*m;
    Guu_degs   = tmp + 4*m;

    /* uctx is context for Z[y_2,...,y_{m - 1}]*/
    fq_nmod_mpoly_ctx_init(uctx, m - 2, ORD_LEX, ctx->fqctx);

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

    fq_nmod_mpolyu_init(Auu, wbits, uctx);
    fq_nmod_mpolyu_init(Buu, wbits, uctx);
    fq_nmod_mpolyu_init(Guu, wbits, uctx);
    fq_nmod_mpolyu_init(Abaruu, wbits, uctx);
    fq_nmod_mpolyu_init(Bbaruu, wbits, uctx);
    fq_nmod_mpoly_init3(Ac, 0, wbits, uctx);
    fq_nmod_mpoly_init3(Bc, 0, wbits, uctx);
    fq_nmod_mpoly_init3(Gc, 0, wbits, uctx);
    fq_nmod_mpoly_init3(Gamma, 0, wbits, uctx);

    fq_nmod_mpoly_to_mpolyuu_perm_deflate(Auu, uctx, A, ctx,
                                         I->bma_perm, I->Amin_exp, I->Gstride);
    fq_nmod_mpoly_to_mpolyuu_perm_deflate(Buu, uctx, B, ctx,
                                         I->bma_perm, I->Bmin_exp, I->Gstride);

    success = fq_nmod_mpolyu_content_mpoly(Ac, Auu, uctx) &&
              fq_nmod_mpolyu_content_mpoly(Bc, Buu, uctx) &&
              fq_nmod_mpoly_gcd(Gc, Ac, Bc, uctx);
    if (!success)
        goto cleanup;

    fq_nmod_mpoly_repack_bits_inplace(Gc, wbits, uctx);

    fq_nmod_mpoly_degrees_si(tmp, Ac, uctx);
    for (i = 0; i < m - 2; i++)
        Auu_degs[i] -= tmp[i];

    fq_nmod_mpolyu_divexact_mpoly_inplace(Auu, Ac, uctx);

    fq_nmod_mpoly_degrees_si(tmp, Bc, uctx);
    for (i = 0; i < m - 2; i++)
        Buu_degs[i] -= tmp[i];

    fq_nmod_mpolyu_divexact_mpoly_inplace(Buu, Bc, uctx);

    fq_nmod_mpoly_degrees_si(tmp, Gc, uctx);
    for (i = 0; i < m - 2; i++)
        Guu_degs[i] -= tmp[i];

    FLINT_ASSERT(Auu->bits == wbits);
    FLINT_ASSERT(Buu->bits == wbits);
    FLINT_ASSERT(Auu->length > 1);
    FLINT_ASSERT(Buu->length > 1);
    FLINT_ASSERT(Ac->bits == wbits);
    FLINT_ASSERT(Bc->bits == wbits);

    success = fq_nmod_mpoly_gcd(Gamma, Auu->coeffs + 0, Buu->coeffs + 0, uctx);
    if (!success)
        goto cleanup;

    fq_nmod_mpoly_repack_bits_inplace(Gamma, wbits, uctx);

    fq_nmod_mpoly_degrees_si(Gamma_degs, Gamma, uctx);

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

            fq_nmod_mpolyu_repack_bits_inplace(Guu, wbits, uctx);
            fq_nmod_mpolyu_repack_bits_inplace(Abaruu, wbits, uctx);
            fq_nmod_mpolyu_repack_bits_inplace(Bbaruu, wbits, uctx);
            fq_nmod_mpolyu_repack_bits_inplace(Auu, wbits, uctx);
            fq_nmod_mpolyu_repack_bits_inplace(Buu, wbits, uctx);
            fq_nmod_mpoly_repack_bits_inplace(Gamma, wbits, uctx);
            fq_nmod_mpoly_repack_bits_inplace(Gc, wbits, uctx);
            break;
        }
    }

    Gguess = I->Gdeflate_deg_bounds_are_nice ? Guu_degs : NULL;

    success = fq_nmod_mpolyuu_gcd_zippel_smprime(Guu, Gguess, Abaruu, Bbaruu,
                        Auu, Auu_degs, Buu, Buu_degs, Gamma, Gamma_degs, uctx);
    if (!success)
    {
        success = fq_nmod_mpolyuu_gcd_zippel_lgprime(Guu, Gguess, Abaruu, Bbaruu,
                        Auu, Auu_degs, Buu, Buu_degs, Gamma, Gamma_degs, uctx);
        if (!success)
            goto cleanup;
    }

    fq_nmod_mpolyu_mul_mpoly_inplace(Guu, Gc, uctx);

    fq_nmod_mpoly_from_mpolyuu_perm_inflate(G, I->Gbits, ctx, Guu, uctx,
                                         I->bma_perm, I->Gmin_exp, I->Gstride);
    success = 1;

cleanup:

    fq_nmod_mpolyu_clear(Auu, uctx);
    fq_nmod_mpolyu_clear(Buu, uctx);
    fq_nmod_mpolyu_clear(Guu, uctx);
    fq_nmod_mpolyu_clear(Abaruu, uctx);
    fq_nmod_mpolyu_clear(Bbaruu, uctx);
    fq_nmod_mpoly_clear(Ac, uctx);
    fq_nmod_mpoly_clear(Bc, uctx);
    fq_nmod_mpoly_clear(Gc, uctx);
    fq_nmod_mpoly_clear(Gamma, uctx);

    fq_nmod_mpoly_ctx_clear(uctx);

    flint_free(tmp);

flint_printf("_try_zippel2 returning %d\n", success);

    return success;
}

static int _try_zippel2(
    fq_nmod_mpoly_t G,
    fq_nmod_mpoly_t Abar,
    fq_nmod_mpoly_t Bbar,
    const fq_nmod_mpoly_t A,
    const fq_nmod_mpoly_t B,
    const mpoly_gcd_info_t I,
    const fq_nmod_mpoly_ctx_t ctx)
{
    slong i, k;
    slong m = I->mvars;
    int success;
    flint_bitcnt_t wbits;
    fq_nmod_mpoly_ctx_t lctx;
    fq_nmod_mpoly_t Al, Bl, Gl, Abarl, Bbarl;
    fq_nmod_mpoly_t Al_lc, Bl_lc, Ac, Bc, Gc, Abarc, Bbarc, Gamma;
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
    Gl_degs   = tmp + 3*m;
    Gamma_degs = tmp + 4*m;

    fq_nmod_mpoly_ctx_init(lctx, m, ORD_LEX, ctx->fqctx);

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

    fq_nmod_mpoly_init3(Al, 0, wbits, lctx);
    fq_nmod_mpoly_init3(Bl, 0, wbits, lctx);
    fq_nmod_mpoly_init3(Gl, 0, wbits, lctx);
    fq_nmod_mpoly_init3(Abarl, 0, wbits, lctx);
    fq_nmod_mpoly_init3(Bbarl, 0, wbits, lctx);
    fq_nmod_mpoly_init3(Ac, 0, wbits, lctx);
    fq_nmod_mpoly_init3(Bc, 0, wbits, lctx);
    fq_nmod_mpoly_init3(Gc, 0, wbits, lctx);
    fq_nmod_mpoly_init3(Abarc, 0, wbits, lctx);
    fq_nmod_mpoly_init3(Bbarc, 0, wbits, lctx);
    fq_nmod_mpoly_init3(Gamma, 0, wbits, lctx);
    fq_nmod_mpoly_init3(Al_lc, 0, wbits, lctx);
    fq_nmod_mpoly_init3(Bl_lc, 0, wbits, lctx);

    fq_nmod_mpoly_to_mpolyl_perm_deflate(Al, lctx, A, ctx, I->bma_perm, I->Amin_exp, I->Gstride);
    fq_nmod_mpoly_to_mpolyl_perm_deflate(Bl, lctx, B, ctx, I->bma_perm, I->Bmin_exp, I->Gstride);

    success = fq_nmod_mpolyl_content(Ac, Al, 2, lctx) &&
              fq_nmod_mpolyl_content(Bc, Bl, 2, lctx);
    if (!success)
        goto cleanup;

    if (Abar == NULL && Bbar == NULL)
        success = fq_nmod_mpoly_gcd(Gc, Ac, Bc, lctx);
    else
        success = fq_nmod_mpoly_gcd_cofactors(Gc, Abarc, Bbarc, Ac, Bc, lctx);
    if (!success)
        goto cleanup;

    fq_nmod_mpoly_degrees_si(tmp, Ac, lctx);
    for (i = 0; i < m; i++)
        Al_degs[i] -= tmp[i];

    success = fq_nmod_mpoly_divides(Al, Al, Ac, lctx);
    FLINT_ASSERT(success);

    fq_nmod_mpoly_degrees_si(tmp, Bc, lctx);
    for (i = 0; i < m; i++)
        Bl_degs[i] -= tmp[i];

    success = fq_nmod_mpoly_divides(Bl, Bl, Bc, lctx);
    FLINT_ASSERT(success);

    fq_nmod_mpoly_degrees_si(tmp, Gc, lctx);
    for (i = 0; i < m; i++)
        Gl_degs[i] -= tmp[i];

    fq_nmod_mpoly_repack_bits_inplace(Al, wbits, lctx);
    fq_nmod_mpoly_repack_bits_inplace(Bl, wbits, lctx);
    fq_nmod_mpolyl_lead_coeff(Al_lc, Al, 2, lctx);
    fq_nmod_mpolyl_lead_coeff(Bl_lc, Bl, 2, lctx);
    success = fq_nmod_mpoly_gcd(Gamma, Al_lc, Bl_lc, lctx);
    if (!success)
        goto cleanup;
    fq_nmod_mpoly_repack_bits_inplace(Gamma, wbits, lctx);

    fq_nmod_mpoly_degrees_si(Gamma_degs, Gamma, lctx);

    Gguess = I->Gdeflate_deg_bounds_are_nice ? Gl_degs : NULL;

    success = fq_nmod_mpolyl_gcd_zippel_smprime(Gl, Gguess, Abarl, Bbarl,
                            Al, Al_degs, Bl, Bl_degs, Gamma, Gamma_degs, lctx);
    if (!success)
    {
        success = fq_nmod_mpolyl_gcd_zippel_lgprime(Gl, Gguess, Abarl, Bbarl,
                            Al, Al_degs, Bl, Bl_degs, Gamma, Gamma_degs, lctx);
        if (!success)
            goto cleanup;
    }

    fq_nmod_mpoly_mul(Gl, Gl, Gc, lctx);
    fq_nmod_mpoly_from_mpolyl_perm_inflate(G, I->Gbits, ctx, Gl, lctx,
                                         I->bma_perm, I->Gmin_exp, I->Gstride);
    if (Abar != NULL)
    {
        fq_nmod_mpoly_mul(Abarl, Abarl, Abarc, lctx);
        fq_nmod_mpoly_from_mpolyl_perm_inflate(Abar, I->Abarbits, ctx, Abarl, lctx,
                                      I->bma_perm, I->Abarmin_exp, I->Gstride);
    }

    if (Bbar != NULL)
    {
        fq_nmod_mpoly_mul(Bbarl, Bbarl, Bbarc, lctx);
        fq_nmod_mpoly_from_mpolyl_perm_inflate(Bbar, I->Bbarbits, ctx, Bbarl, lctx,
                                      I->bma_perm, I->Bbarmin_exp, I->Gstride);
    }

    success = 1;

cleanup:

    fq_nmod_mpoly_clear(Al, lctx);
    fq_nmod_mpoly_clear(Bl, lctx);
    fq_nmod_mpoly_clear(Gl, lctx);
    fq_nmod_mpoly_clear(Abarl, lctx);
    fq_nmod_mpoly_clear(Bbarl, lctx);
    fq_nmod_mpoly_clear(Ac, lctx);
    fq_nmod_mpoly_clear(Bc, lctx);
    fq_nmod_mpoly_clear(Gc, lctx);
    fq_nmod_mpoly_clear(Abarc, lctx);
    fq_nmod_mpoly_clear(Bbarc, lctx);
    fq_nmod_mpoly_clear(Gamma, lctx);
    fq_nmod_mpoly_clear(Al_lc, lctx);
    fq_nmod_mpoly_clear(Bl_lc, lctx);

    fq_nmod_mpoly_ctx_clear(lctx);

    flint_free(tmp);

    return success;
}

/******************** Hit A and B with hensel lifting ************************/
static int _try_hensel(
    fq_nmod_mpoly_t G,
    fq_nmod_mpoly_t Abar,
    fq_nmod_mpoly_t Bbar,
    const fq_nmod_mpoly_t A,
    const fq_nmod_mpoly_t B,
    const mpoly_gcd_info_t I,
    const fq_nmod_mpoly_ctx_t ctx)
{
    slong i, k;
    slong m = I->mvars;
    int success;
    flint_bitcnt_t wbits;
    fq_nmod_mpoly_ctx_t lctx;
    fq_nmod_mpoly_t Al, Bl, Gl, Abarl, Bbarl;
    fq_nmod_mpoly_t Ac, Bc, Gc, Abarc, Bbarc;
    slong max_degree;

    FLINT_ASSERT(A->bits <= FLINT_BITS);
    FLINT_ASSERT(B->bits <= FLINT_BITS);
    FLINT_ASSERT(A->length > 0);
    FLINT_ASSERT(B->length > 0);

    if (!I->can_use_brown)
        return 0;

    FLINT_ASSERT(m >= WORD(2));

    fq_nmod_mpoly_ctx_init(lctx, m, ORD_LEX, ctx->fqctx);

    max_degree = 0;
    for (i = 0; i < m; i++)
    {
        k = I->bma_perm[i];
        max_degree = FLINT_MAX(max_degree, I->Adeflate_deg[k]);
        max_degree = FLINT_MAX(max_degree, I->Bdeflate_deg[k]);
    }

    wbits = 1 + FLINT_BIT_COUNT(max_degree);
    wbits = mpoly_fix_bits(wbits, lctx->minfo);
    FLINT_ASSERT(wbits <= FLINT_BITS);

    fq_nmod_mpoly_init3(Al, 0, wbits, lctx);
    fq_nmod_mpoly_init3(Bl, 0, wbits, lctx);
    fq_nmod_mpoly_init3(Gl, 0, wbits, lctx);
    fq_nmod_mpoly_init3(Abarl, 0, wbits, lctx);
    fq_nmod_mpoly_init3(Bbarl, 0, wbits, lctx);
    fq_nmod_mpoly_init3(Ac, 0, wbits, lctx);
    fq_nmod_mpoly_init3(Bc, 0, wbits, lctx);
    fq_nmod_mpoly_init3(Gc, 0, wbits, lctx);
    fq_nmod_mpoly_init3(Abarc, 0, wbits, lctx);
    fq_nmod_mpoly_init3(Bbarc, 0, wbits, lctx);

    fq_nmod_mpoly_to_mpolyl_perm_deflate(Al, lctx, A, ctx, I->bma_perm, I->Amin_exp, I->Gstride);
    fq_nmod_mpoly_to_mpolyl_perm_deflate(Bl, lctx, B, ctx, I->bma_perm, I->Bmin_exp, I->Gstride);

    success = fq_nmod_mpolyl_content(Ac, Al, 1, lctx) &&
              fq_nmod_mpolyl_content(Bc, Bl, 1, lctx);
    if (!success)
        goto cleanup;

    if (Abar == NULL && Bbar == NULL)
        success = fq_nmod_mpoly_gcd(Gc, Ac, Bc, lctx);
    else
        success = fq_nmod_mpoly_gcd_cofactors(Gc, Abarc, Bbarc, Ac, Bc, lctx);
    if (!success)
        goto cleanup;

    success = fq_nmod_mpoly_divides(Al, Al, Ac, lctx);
    FLINT_ASSERT(success);

    success = fq_nmod_mpoly_divides(Bl, Bl, Bc, lctx);
    FLINT_ASSERT(success);

    fq_nmod_mpoly_repack_bits_inplace(Al, wbits, lctx);
    fq_nmod_mpoly_repack_bits_inplace(Bl, wbits, lctx);
    success = fq_nmod_mpolyl_gcd_hensel_smprime(Gl, Abarl, Bbarl, Al, Bl, lctx);
    if (!success)
        goto cleanup;

    fq_nmod_mpoly_mul(Gl, Gl, Gc, lctx);
    fq_nmod_mpoly_from_mpolyl_perm_inflate(G, I->Gbits, ctx, Gl, lctx,
                                         I->bma_perm, I->Gmin_exp, I->Gstride);
    if (Abar != NULL)
    {
        fq_nmod_mpoly_mul(Abarl, Abarl, Abarc, lctx);
        fq_nmod_mpoly_from_mpolyl_perm_inflate(Abar, I->Abarbits, ctx, Abarl, lctx,
                                      I->bma_perm, I->Abarmin_exp, I->Gstride);
    }

    if (Bbar != NULL)
    {
        fq_nmod_mpoly_mul(Bbarl, Bbarl, Bbarc, lctx);
        fq_nmod_mpoly_from_mpolyl_perm_inflate(Bbar, I->Bbarbits, ctx, Bbarl, lctx,
                                      I->bma_perm, I->Bbarmin_exp, I->Gstride);
    }

    success = 1;

cleanup:

    fq_nmod_mpoly_clear(Al, lctx);
    fq_nmod_mpoly_clear(Bl, lctx);
    fq_nmod_mpoly_clear(Gl, lctx);
    fq_nmod_mpoly_clear(Abarl, lctx);
    fq_nmod_mpoly_clear(Bbarl, lctx);
    fq_nmod_mpoly_clear(Ac, lctx);
    fq_nmod_mpoly_clear(Bc, lctx);
    fq_nmod_mpoly_clear(Gc, lctx);

    fq_nmod_mpoly_ctx_clear(lctx);

    return success;
}


/*********************** Hit A and B with brown ******************************/
static int _try_brown(
    fq_nmod_mpoly_t G,
    fq_nmod_mpoly_t Abar,
    fq_nmod_mpoly_t Bbar,
    const fq_nmod_mpoly_t A,
    const fq_nmod_mpoly_t B,
    mpoly_gcd_info_t I,
    const fq_nmod_mpoly_ctx_t ctx)
{
    int success;
    slong m = I->mvars;
    flint_bitcnt_t wbits;
    fq_nmod_mpoly_ctx_t nctx;
    fq_nmod_mpolyn_t An, Bn, Gn, Abarn, Bbarn;

    if (!I->can_use_brown)
        return 0;

    FLINT_ASSERT(m >= 2);
    FLINT_ASSERT(A->bits <= FLINT_BITS);
    FLINT_ASSERT(B->bits <= FLINT_BITS);
    FLINT_ASSERT(A->length > 0);
    FLINT_ASSERT(B->length > 0);

    wbits = FLINT_MAX(A->bits, B->bits);

    fq_nmod_mpoly_ctx_init(nctx, m, ORD_LEX, ctx->fqctx);
    fq_nmod_mpolyn_init(An, wbits, nctx);
    fq_nmod_mpolyn_init(Bn, wbits, nctx);
    fq_nmod_mpolyn_init(Gn, wbits, nctx);
    fq_nmod_mpolyn_init(Abarn, wbits, nctx);
    fq_nmod_mpolyn_init(Bbarn, wbits, nctx);

    fq_nmod_mpoly_to_mpolyn_perm_deflate(An, nctx, A, ctx,
                                       I->brown_perm, I->Amin_exp, I->Gstride);
    fq_nmod_mpoly_to_mpolyn_perm_deflate(Bn, nctx, B, ctx,
                                       I->brown_perm, I->Bmin_exp, I->Gstride);

    FLINT_ASSERT(An->bits == wbits);
    FLINT_ASSERT(Bn->bits == wbits);
    FLINT_ASSERT(An->length > 1);
    FLINT_ASSERT(Bn->length > 1);

    success = fq_nmod_mpolyn_gcd_brown_smprime(Gn, Abarn, Bbarn, An, Bn,
                                                                  m - 1, nctx);
    if (!success)
    {
        fq_nmod_mpoly_to_mpolyn_perm_deflate(An, nctx, A, ctx,
                                       I->brown_perm, I->Amin_exp, I->Gstride);
        fq_nmod_mpoly_to_mpolyn_perm_deflate(Bn, nctx, B, ctx,
                                       I->brown_perm, I->Bmin_exp, I->Gstride);
        success = fq_nmod_mpolyn_gcd_brown_lgprime(Gn, Abarn, Bbarn, An, Bn,
                                                                  m - 1, nctx);
    }

    if (!success)
        goto cleanup;

    fq_nmod_mpoly_from_mpolyn_perm_inflate(G, I->Gbits, ctx, Gn, nctx,
                                       I->brown_perm, I->Gmin_exp, I->Gstride);

    if (Abar != NULL)
        fq_nmod_mpoly_from_mpolyn_perm_inflate(Abar, I->Abarbits, ctx, Abarn, nctx,
                                    I->brown_perm, I->Abarmin_exp, I->Gstride);

    if (Bbar != NULL)
        fq_nmod_mpoly_from_mpolyn_perm_inflate(Bbar, I->Bbarbits, ctx, Bbarn, nctx,
                                    I->brown_perm, I->Bbarmin_exp, I->Gstride);

    success = 1;

cleanup:

    fq_nmod_mpolyn_clear(An, nctx);
    fq_nmod_mpolyn_clear(Bn, nctx);
    fq_nmod_mpolyn_clear(Gn, nctx);
    fq_nmod_mpolyn_clear(Abarn, nctx);
    fq_nmod_mpolyn_clear(Bbarn, nctx);
    fq_nmod_mpoly_ctx_clear(nctx);

    return success;
}


/*
    Both A and B have to be packed into bits <= FLINT_BITS
    return is 1 for success, 0 for failure.
*/
static int _fq_nmod_mpoly_gcd_algo_small(
    fq_nmod_mpoly_t G,
    fq_nmod_mpoly_t Abar, /* could be NULL */
    fq_nmod_mpoly_t Bbar, /* could be NULL */
    const fq_nmod_mpoly_t A,
    const fq_nmod_mpoly_t B,
    const fq_nmod_mpoly_ctx_t ctx,
    unsigned int algo)
{
    int success;
    flint_bitcnt_t Gbits = FLINT_MIN(A->bits, B->bits);
    flint_bitcnt_t Abarbits = A->bits;
    flint_bitcnt_t Bbarbits = B->bits;
    slong v_in_both;
    slong v_in_either;
    slong v_in_A_only;
    slong v_in_B_only;
    slong j;
    slong nvars = ctx->minfo->nvars;
    mpoly_gcd_info_t I;
#if FLINT_WANT_ASSERT
    fq_nmod_mpoly_t T, Asave, Bsave;
#endif

    if (A->length == 1)
        return _do_monomial_gcd(G, Bbar, Abar, B, A, ctx);
    else if (B->length == 1)
        return _do_monomial_gcd(G, Abar, Bbar, A, B, ctx);

#if FLINT_WANT_ASSERT
    fq_nmod_mpoly_init(T, ctx);
    fq_nmod_mpoly_init(Asave, ctx);
    fq_nmod_mpoly_init(Bsave, ctx);
    fq_nmod_mpoly_set(Asave, A, ctx);
    fq_nmod_mpoly_set(Bsave, B, ctx);
#endif

    mpoly_gcd_info_init(I, nvars);

    /* entries of I are all now invalid */

    I->Gbits = Gbits;
    I->Abarbits = Abarbits;
    I->Bbarbits = Bbarbits;

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

    if (_try_monomial_cofactors(G, Abar, Bbar, A, B, ctx))
    {
        goto successful;
    }

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

        t = FLINT_MIN(I->Amin_exp[j], I->Bmin_exp[j]);
        I->Gmin_exp[j] = t;
        I->Abarmin_exp[j] = I->Amin_exp[j] - t;
        I->Bbarmin_exp[j] = I->Bmin_exp[j] - t;
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
        _do_trivial(G, Abar, Bbar, A, B, I, ctx);
        goto successful;
    }

    /* check if ess(A) and ess(B) depend on another variable v_in_either */
    FLINT_ASSERT(0 <= v_in_both && v_in_both < nvars);

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
        _do_univar(G, Abar, Bbar, A, B, v_in_both, I, ctx);
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
        success = _try_missing_var(G, I->Gbits, Abar, Bbar, v_in_A_only,
                                   A, I->Amin_exp[v_in_A_only],
                                   B, I->Bmin_exp[v_in_A_only],
                                   ctx);
        goto cleanup;
    }
    if (v_in_B_only != -WORD(1))
    {
        success = _try_missing_var(G, I->Gbits, Bbar, Abar, v_in_B_only,
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

    mpoly_gcd_info_set_estimates_fq_nmod_mpoly(I, A, B, ctx);
    if (!I->Gdeflate_deg_bounds_are_nice)
        mpoly_gcd_info_set_estimates_fq_nmod_mpoly_lgprime(I, A, B, ctx);
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
        {
            _do_trivial(G, Abar, Bbar, A, B, I, ctx);
            goto successful;
        }

        if (try_a && _try_divides(G, Bbar, Abar, B, A, ctx))
            goto successful;

        if (try_b && _try_divides(G, Abar, Bbar, A, B, ctx))
            goto successful;
    }

    if (I->mvars < 3)
    {
        mpoly_gcd_info_measure_brown(I, A->length, B->length, ctx->minfo);

        if (_try_hensel(G, Abar, Bbar, A, B, I, ctx))
            goto successful;

        if (_try_brown(G, Abar, Bbar, A, B, I, ctx))
            goto successful;
    }
    else
    {
        mpoly_gcd_info_measure_brown(I, A->length, B->length, ctx->minfo);
        mpoly_gcd_info_measure_zippel(I, A->length, B->length, ctx->minfo);
        mpoly_gcd_info_measure_zippel2(I, A->length, B->length, ctx->minfo);

        if (I->zippel_time_est < I->brown_time_est)
        {
            if (_try_zippel2(G, Abar, Bbar, A, B, I, ctx))
                goto successful;
            if (_try_zippel2_old(G, A, B, I, ctx))
                goto successful;
            if (_try_zippel(G, Abar, Bbar, A, B, I, ctx))
                goto successful;
            if (_try_brown(G, Abar, Bbar, A, B, I, ctx))
                goto successful;
        }
        else
        {
            if (_try_zippel2(G, Abar, Bbar, A, B, I, ctx))
                goto successful;
            if (_try_zippel2_old(G, A, B, I, ctx))
                goto successful;
            if (_try_brown(G, Abar, Bbar, A, B, I, ctx))
                goto successful;
            if (_try_zippel(G, Abar, Bbar, A, B, I, ctx))
                goto successful;
        }
    }

    if (_try_hensel(G, Abar, Bbar, A, B, I, ctx))
        goto successful;

    success = 0;
    goto cleanup;

successful:

    success = 1;

cleanup:

    mpoly_gcd_info_clear(I);

    if (success)
    {
        FLINT_ASSERT(G->length > 0);

        if (!n_fq_is_one(G->coeffs + 0, ctx->fqctx))
        {
            if (Abar != NULL)
                fq_nmod_mpoly_scalar_mul_n_fq(Abar, Abar, G->coeffs + 0, ctx);

            if (Bbar != NULL)
                fq_nmod_mpoly_scalar_mul_n_fq(Bbar, Bbar, G->coeffs + 0, ctx);

            fq_nmod_mpoly_make_monic(G, G, ctx);
        }

        FLINT_ASSERT(fq_nmod_mpoly_divides(T, Asave, G, ctx));
        FLINT_ASSERT(Abar == NULL || fq_nmod_mpoly_equal(T, Abar, ctx));
            
        FLINT_ASSERT(fq_nmod_mpoly_divides(T, Bsave, G, ctx));
        FLINT_ASSERT(Bbar == NULL || fq_nmod_mpoly_equal(T, Bbar, ctx));
    }

#if FLINT_WANT_ASSERT
    fq_nmod_mpoly_clear(T, ctx);
    fq_nmod_mpoly_clear(Asave, ctx);
    fq_nmod_mpoly_clear(Bsave, ctx);
#endif

    return success;
}


/*
    The gcd calculation is unusual.
    First see if both inputs fit into FLINT_BITS.
    Then, try deflation as a last resort.
*/
static int _fq_nmod_mpoly_gcd_algo_large(
    fq_nmod_mpoly_t G,
    fq_nmod_mpoly_t Abar,
    fq_nmod_mpoly_t Bbar,
    const fq_nmod_mpoly_t A,
    const fq_nmod_mpoly_t B,
    const fq_nmod_mpoly_ctx_t ctx,
    unsigned int algo)
{
    int success;
    slong k;
    fmpz * Ashift, * Astride;
    fmpz * Bshift, * Bstride;
    fmpz * Gshift, * Gstride;
    fq_nmod_mpoly_t Anew, Bnew;
    const fq_nmod_mpoly_struct * Ause, * Buse;

    if (A->length == 1)
        return _do_monomial_gcd(G, Bbar, Abar, B, A, ctx);

    if (B->length == 1)
        return _do_monomial_gcd(G, Abar, Bbar, A, B, ctx);

    if (_try_monomial_cofactors(G, Abar, Bbar, A, B, ctx))
        return 1;

    fq_nmod_mpoly_init(Anew, ctx);
    fq_nmod_mpoly_init(Bnew, ctx);

    Ause = A;
    if (A->bits > FLINT_BITS)
    {
       if (!fq_nmod_mpoly_repack_bits(Anew, A, FLINT_BITS, ctx))
            goto could_not_repack;
        Ause = Anew;
    }

    Buse = B;
    if (B->bits > FLINT_BITS)
    {
        if (!fq_nmod_mpoly_repack_bits(Bnew, B, FLINT_BITS, ctx))
            goto could_not_repack;
        Buse = Bnew;
    }

    success = _fq_nmod_mpoly_gcd_algo_small(G, Abar, Bbar, Ause, Buse, ctx, algo);

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

    fq_nmod_mpoly_deflation(Ashift, Astride, A, ctx);
    fq_nmod_mpoly_deflation(Bshift, Bstride, B, ctx);
    _fmpz_vec_min(Gshift, Ashift, Bshift, ctx->minfo->nvars);
    for (k = 0; k < ctx->minfo->nvars; k++)
    {
        fmpz_gcd(Gstride + k, Astride + k, Bstride + k);
    }

    fq_nmod_mpoly_deflate(Anew, A, Ashift, Gstride, ctx);
    if (Anew->bits > FLINT_BITS)
    {
        success = fq_nmod_mpoly_repack_bits(Anew, Anew, FLINT_BITS, ctx);
        if (!success)
            goto deflate_cleanup;
    }

    fq_nmod_mpoly_deflate(Bnew, B, Bshift, Gstride, ctx);
    if (Bnew->bits > FLINT_BITS)
    {
        success = fq_nmod_mpoly_repack_bits(Bnew, Bnew, FLINT_BITS, ctx);
        if (!success)
            goto deflate_cleanup;
    }

    success = _fq_nmod_mpoly_gcd_algo(G, Abar, Bbar, Anew, Bnew, ctx, algo);
    if (!success)
        goto deflate_cleanup;

    for (k = 0; k < ctx->minfo->nvars; k++)
    {
        fmpz_sub(Ashift + k, Ashift + k, Gshift + k);
        fmpz_sub(Bshift + k, Bshift + k, Gshift + k);
        FLINT_ASSERT(fmpz_sgn(Ashift + k) >= 0);
        FLINT_ASSERT(fmpz_sgn(Bshift + k) >= 0);
    }

    fq_nmod_mpoly_inflate(G, G, Gshift, Gstride, ctx);
    if (Abar != NULL)
        fq_nmod_mpoly_inflate(Abar, Abar, Ashift, Gstride, ctx);
    if (Bbar != NULL)
        fq_nmod_mpoly_inflate(Bbar, Bbar, Bshift, Gstride, ctx);

    FLINT_ASSERT(G->length > 0);
    if (!n_fq_is_one(G->coeffs + 0, ctx->fqctx))
    {
        if (Abar != NULL)
            fq_nmod_mpoly_scalar_mul_n_fq(Abar, Abar, G->coeffs + 0, ctx);

        if (Bbar != NULL)
            fq_nmod_mpoly_scalar_mul_n_fq(Bbar, Bbar, G->coeffs + 0, ctx);

        fq_nmod_mpoly_make_monic(G, G, ctx);
    }
    fq_nmod_mpoly_make_monic(G, G, ctx);

deflate_cleanup:

    _fmpz_vec_clear(Ashift, ctx->minfo->nvars);
    _fmpz_vec_clear(Astride, ctx->minfo->nvars);
    _fmpz_vec_clear(Bshift, ctx->minfo->nvars);
    _fmpz_vec_clear(Bstride, ctx->minfo->nvars);
    _fmpz_vec_clear(Gshift, ctx->minfo->nvars);
    _fmpz_vec_clear(Gstride, ctx->minfo->nvars);

cleanup:

    fq_nmod_mpoly_clear(Anew, ctx);
    fq_nmod_mpoly_clear(Bnew, ctx);

    return success;
}

int _fq_nmod_mpoly_gcd_algo(
    fq_nmod_mpoly_t G,
    fq_nmod_mpoly_t Abar,
    fq_nmod_mpoly_t Bbar,
    const fq_nmod_mpoly_t A,
    const fq_nmod_mpoly_t B,
    const fq_nmod_mpoly_ctx_t ctx,
    unsigned int algo)
{
    if (A->bits <= FLINT_BITS && B->bits <= FLINT_BITS)
        return _fq_nmod_mpoly_gcd_algo_small(G, Abar, Bbar, A, B, ctx, algo);
    else
        return _fq_nmod_mpoly_gcd_algo_large(G, Abar, Bbar, A, B, ctx, algo);
}


int fq_nmod_mpoly_gcd(
    fq_nmod_mpoly_t G,
    const fq_nmod_mpoly_t A,
    const fq_nmod_mpoly_t B,
    const fq_nmod_mpoly_ctx_t ctx)
{
    if (fq_nmod_mpoly_is_zero(A, ctx))
    {
        if (fq_nmod_mpoly_is_zero(B, ctx))
            fq_nmod_mpoly_zero(G, ctx);
        else
            fq_nmod_mpoly_make_monic(G, B, ctx);
        return 1;
    }

    if (fq_nmod_mpoly_is_zero(B, ctx))
    {
        fq_nmod_mpoly_make_monic(G, A, ctx);
        return 1;
    }

    return _fq_nmod_mpoly_gcd_algo(G, NULL, NULL, A, B, ctx, MPOLY_GCD_USE_ALL);
}

