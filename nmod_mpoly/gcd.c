/*
    Copyright (C) 2018, 2019 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "nmod_mpoly.h"

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
    nmod_poly_struct * out,
    const int * ignore,
    const nmod_mpoly_t A,
    ulong * Amin_exp,
    ulong * Amax_exp,
    ulong * Astride,
    mp_limb_t * alpha,
    const nmod_mpoly_ctx_t ctx,
    const thread_pool_handle * handles,
    slong num_handles)
{
    slong i, j;
    slong nvars = ctx->minfo->nvars;
    slong total_limit, total_length;
    int use_direct_LUT;
    ulong varexp;
    ulong mask;
    slong * offsets, * shifts;
    slong N = mpoly_words_per_exp_sp(A->bits, ctx->minfo);
    ulong * Aexp = A->exps;
    mp_limb_t * Acoeff = A->coeffs;
    mp_limb_t meval;
    mp_limb_t t;

    FLINT_ASSERT(A->bits <= FLINT_BITS);

    mask = (-UWORD(1)) >> (FLINT_BITS - A->bits);
    offsets = (slong *) flint_malloc(ctx->minfo->nvars*sizeof(slong));
    shifts = (slong *) flint_malloc(ctx->minfo->nvars*sizeof(slong));

    for (j = 0; j < ctx->minfo->nvars; j++)
    {
        nmod_poly_zero(out + j);
        mpoly_gen_offset_shift_sp(offsets + j, shifts + j, j, A->bits, ctx->minfo);
    }

    /*
        two cases:
        (1) the Amax_exp[j] are small enough to calculate a direct LUT
        (2) use a LUT for exponents that are powers of two
    */

    total_limit = A->length/256;
    total_limit = FLINT_MAX(WORD(9999), total_limit);
    total_length = 0;
    use_direct_LUT = 1;
    for (j = 0; j < ctx->minfo->nvars; j++)
    {
        total_length += Amax_exp[j] + 1;
        if ((ulong) total_length > (ulong) total_limit)
            use_direct_LUT = 0;
    }

    if (use_direct_LUT)
    {
        slong off;
        mp_limb_t * LUT, ** LUTvalue, ** LUTvalueinv;

        /* value of powers of alpha[j] */
        LUT = (mp_limb_t *) flint_malloc(2*total_length*sizeof(mp_limb_t));

        /* pointers into LUT */
        LUTvalue    = (mp_limb_t **) flint_malloc(nvars*sizeof(mp_limb_t *));
        LUTvalueinv = (mp_limb_t **) flint_malloc(nvars*sizeof(mp_limb_t *));

        off = 0;
        for (j = 0; j < nvars; j++)
        {
            ulong k;
            mp_limb_t alphainvj = nmod_inv(alpha[j], (out + 0)->mod);

            LUTvalue[j] = LUT + off;
            LUTvalueinv[j] = LUT + total_length + off;
            LUTvalue[j][0] = 1;
            LUTvalueinv[j][0] = 1;
            for (k = 0; k < Amax_exp[j]; k++)
            {
                LUTvalue[j][k + 1] = nmod_mul(LUTvalue[j][k], alpha[j],
                                                               (out + 0)->mod);
                LUTvalueinv[j][k + 1] = nmod_mul(LUTvalueinv[j][k], alphainvj,
                                                               (out + 0)->mod);
            }

            off += Amax_exp[j] + 1;
        }
        FLINT_ASSERT(off == total_length);

        for (i = 0; i < A->length; i++)
        {
            meval = Acoeff[i];

            for (j = 0; j < nvars; j++)
            {
                varexp = ((Aexp + N*i)[offsets[j]]>>shifts[j])&mask;
                FLINT_ASSERT(varexp <= Amax_exp[j]);
                meval = nmod_mul(meval, LUTvalue[j][varexp], (out + 0)->mod);
            }

            for (j = 0; j < nvars; j++)
            {
                varexp = ((Aexp + N*i)[offsets[j]]>>shifts[j])&mask;

                if (ignore[j])
                    continue;

                t = nmod_mul(meval, LUTvalueinv[j][varexp], (out + j)->mod);

                FLINT_ASSERT((Astride[j] == 0 && varexp == Amin_exp[j])
                                  || (varexp - Amin_exp[j]) % Astride[j] == 0);

                varexp = Astride[j] < 2 ? varexp - Amin_exp[j] :
                                           (varexp - Amin_exp[j])/Astride[j];

                t = nmod_add(t, nmod_poly_get_coeff_ui(out + j, varexp),
                                                               (out + j)->mod);
                nmod_poly_set_coeff_ui(out + j, varexp, t);
            }
        }

        flint_free(LUT);
        flint_free(LUTvalue);
        flint_free(LUTvalueinv);
    }
    else
    {
        slong LUTlen;
        ulong * LUTmask;
        slong * LUToffset, * LUTvar;
        mp_limb_t * LUTvalue, * LUTvalueinv;
        mp_limb_t * vieval;
        mp_limb_t t, xpoweval, xinvpoweval;

        LUToffset   = (slong *) flint_malloc(N*FLINT_BITS*sizeof(slong));
        LUTmask     = (ulong *) flint_malloc(N*FLINT_BITS*sizeof(ulong));
        LUTvalue    = (mp_limb_t *) flint_malloc(N*FLINT_BITS*sizeof(mp_limb_t));
        LUTvar      = (slong *) flint_malloc(N*FLINT_BITS*sizeof(slong));
        LUTvalueinv = (mp_limb_t *) flint_malloc(N*FLINT_BITS*sizeof(mp_limb_t));

        vieval = (mp_limb_t *) flint_malloc(nvars*sizeof(mp_limb_t));

        LUTlen = 0;
        for (j = nvars - 1; j >= 0; j--)
        {
            flint_bitcnt_t bits = FLINT_BIT_COUNT(Amax_exp[j]);
            xpoweval = alpha[j]; /* xpoweval = alpha[j]^(2^i) */
            xinvpoweval = nmod_inv(xpoweval, (out + 0)->mod); /* alpha[j]^(-2^i) */
            for (i = 0; i < bits; i++)
            {
                LUToffset[LUTlen] = offsets[j];
                LUTmask[LUTlen] = (UWORD(1) << (shifts[j] + i));
                LUTvalue[LUTlen] = xpoweval;
                LUTvalueinv[LUTlen] = xinvpoweval;
                LUTvar[LUTlen] = j;
                LUTlen++;
                xpoweval = nmod_mul(xpoweval, xpoweval, (out + 0)->mod);
                xinvpoweval = nmod_mul(xinvpoweval, xinvpoweval, (out + 0)->mod);
            }

            vieval[j] = 1;
        }
        FLINT_ASSERT(LUTlen < N*FLINT_BITS);

        for (i = 0; i < A->length; i++)
        {
            meval = Acoeff[i];

            for (j = 0; j < LUTlen; j++)
            {
                if (((Aexp + N*i)[LUToffset[j]] & LUTmask[j]) != 0)
                {
                    meval = nmod_mul(meval, LUTvalue[j], (out + 0)->mod);
                    vieval[LUTvar[j]] = nmod_mul(vieval[LUTvar[j]],
                                               LUTvalueinv[j], (out + 0)->mod);
                }
            }

            for (j = 0; j < nvars; j++)
            {
                varexp = ((Aexp + N*i)[offsets[j]]>>shifts[j])&mask;

                FLINT_ASSERT((Astride[j] == 0 && varexp == Amin_exp[j])
                                  || (varexp - Amin_exp[j]) % Astride[j] == 0);

                varexp = Astride[j] < 2 ? varexp - Amin_exp[j] :
                                           (varexp - Amin_exp[j])/Astride[j];

                t = nmod_mul(meval, vieval[j], (out + j)->mod);
                t = nmod_add(t, nmod_poly_get_coeff_ui(out + j, varexp),
                                                               (out + j)->mod);
                nmod_poly_set_coeff_ui(out + j, varexp, t);
                vieval[j] = 1;
            }
        }

        flint_free(LUToffset);
        flint_free(LUTmask);
        flint_free(LUTvalue);
        flint_free(LUTvar);
        flint_free(LUTvalueinv);

        flint_free(vieval);
    }

    flint_free(offsets);
    flint_free(shifts);
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
    slong i, j;
    nmod_poly_t Geval;
    nmod_poly_struct * Aevals, * Bevals;
    mp_limb_t * alpha;
    flint_rand_t randstate;
    slong ignore_limit;
    int * ignore;

    flint_randinit(randstate);

    ignore = (int *) flint_malloc(ctx->minfo->nvars*sizeof(int));
    alpha = (mp_limb_t *) flint_malloc(ctx->minfo->nvars*sizeof(mp_limb_t));
    Aevals = (nmod_poly_struct *) flint_malloc(
                                   ctx->minfo->nvars*sizeof(nmod_poly_struct));
    Bevals = (nmod_poly_struct *) flint_malloc(
                                   ctx->minfo->nvars*sizeof(nmod_poly_struct));

    nmod_poly_init_mod(Geval, ctx->ffinfo->mod);
    for (j = 0; j < ctx->minfo->nvars; j++)
    {
        nmod_poly_init_mod(Aevals + j, ctx->ffinfo->mod);
        nmod_poly_init_mod(Bevals + j, ctx->ffinfo->mod);
    }

    ignore_limit = A->length/4096 + B->length/4096;
    ignore_limit = FLINT_MAX(WORD(9999), ignore_limit);
    I->Gdeflate_deg_bounds_are_nice = 1;
    for (j = 0; j < ctx->minfo->nvars; j++)
    {
        if (   I->Adeflate_deg[j] > ignore_limit
            || I->Bdeflate_deg[j] > ignore_limit)
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
        for (j = 0; j < ctx->minfo->nvars; j++)
        {
            I->Gdeflate_deg_bound[j] = FLINT_MIN(I->Adeflate_deg[j],
                                                 I->Bdeflate_deg[j]);
            I->Gterm_count_est[j] = (I->Gdeflate_deg_bound[j] + 1)/2;
        }

        goto cleanup;
    }

    for (j = 0; j < ctx->minfo->nvars; j++)
    {
        alpha[j] = n_urandint(randstate, ctx->ffinfo->mod.n - 1) + 1;
    }


    nmod_mpoly_evals(Aevals, ignore, A, I->Amin_exp, I->Amax_exp, I->Gstride,
                                             alpha, ctx, handles, num_handles);
    nmod_mpoly_evals(Bevals, ignore, B, I->Bmin_exp, I->Bmax_exp, I->Gstride,
                                             alpha, ctx, handles, num_handles);

    for (j = 0; j < ctx->minfo->nvars; j++)
    {
        if (ignore[j])
        {
            I->Gdeflate_deg_bound[j] = FLINT_MIN(I->Adeflate_deg[j],
                                                 I->Bdeflate_deg[j]);
            I->Gterm_count_est[j] = (I->Gdeflate_deg_bound[j] + 1)/2;
        }
        else
        {
            if (   I->Adeflate_deg[j] != nmod_poly_degree(Aevals + j)
                || I->Bdeflate_deg[j] != nmod_poly_degree(Bevals + j))
            {
                goto try_again;
            }

            nmod_poly_gcd(Geval, Aevals + j, Bevals + j);

            I->Gterm_count_est[j] = 0;
            I->Gdeflate_deg_bound[j] = nmod_poly_degree(Geval);
            for (i = I->Gdeflate_deg_bound[j]; i >= 0; i--)
            {
                I->Gterm_count_est[j] += (Geval->coeffs[i] != 0);
            }
        }
    }

cleanup:

    nmod_poly_clear(Geval);
    for (j = 0; j < ctx->minfo->nvars; j++)
    {
        nmod_poly_clear(Aevals + j);
        nmod_poly_clear(Bevals + j);
    }

    flint_free(ignore);
    flint_free(alpha);
    flint_free(Aevals);
    flint_free(Bevals);

    flint_randclear(randstate);

    return;
}

static int _try_zippel(
    nmod_mpoly_t G,
    flint_bitcnt_t Gbits,
    const nmod_mpoly_t A,
    const nmod_mpoly_t B,
    const mpoly_gcd_info_t I,
    const nmod_mpoly_ctx_t ctx)
{
    slong i, k;
    slong m = I->mvars;
    int success;
    mpoly_zipinfo_t zinfo;
    flint_bitcnt_t ABbits;
    flint_rand_t randstate;
    nmod_mpoly_ctx_t uctx;
    nmod_mpolyu_t Au, Bu, Gu;
    nmod_mpoly_t Acontent, Bcontent;
    nmod_mpolyu_t Abar, Bbar, Gbar;

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

    ABbits = FLINT_MAX(A->bits, B->bits);

    nmod_mpolyu_init(Au, ABbits, uctx);
    nmod_mpolyu_init(Bu, ABbits, uctx);
    nmod_mpolyu_init(Gu, ABbits, uctx);

    nmod_mpoly_to_mpolyu_perm_deflate(Au, uctx, A, ctx,
                        zinfo->perm, I->Amin_exp, I->Gstride, NULL, 0);
    nmod_mpoly_to_mpolyu_perm_deflate(Bu, uctx, B, ctx,
                        zinfo->perm, I->Bmin_exp, I->Gstride, NULL, 0);

    FLINT_ASSERT(Au->bits == ABbits);
    FLINT_ASSERT(Bu->bits == ABbits);
    FLINT_ASSERT(Au->length > 1);
    FLINT_ASSERT(Bu->length > 1);

    nmod_mpoly_init3(Acontent, 0, ABbits, uctx);
    nmod_mpoly_init3(Bcontent, 0, ABbits, uctx);
    nmod_mpolyu_init(Abar, ABbits, uctx);
    nmod_mpolyu_init(Bbar, ABbits, uctx);
    nmod_mpolyu_init(Gbar, ABbits, uctx);

    /* remove content from A and B */
    success = nmod_mpolyu_content_mpoly(Acontent, Au, uctx, NULL, 0);
    success = success
           && nmod_mpolyu_content_mpoly(Bcontent, Bu, uctx, NULL, 0);
    if (!success)
        goto cleanup;

    FLINT_ASSERT(Acontent->bits == ABbits);
    FLINT_ASSERT(Bcontent->bits == ABbits);
    nmod_mpolyu_divexact_mpoly(Abar, Au, Acontent, uctx);
    nmod_mpolyu_divexact_mpoly(Bbar, Bu, Bcontent, uctx);

    /* after removing content, degree bounds in zinfo are still valid bounds */

    /* compute GCD */
    success = nmod_mpolyu_gcdm_zippel(Gbar, Abar, Bbar, uctx, zinfo, randstate);
    if (!success)
        goto cleanup;

    /* put back content */
    success = _nmod_mpoly_gcd(Acontent, ABbits, Acontent, Bcontent, uctx,
                                                                      NULL, 0);
    if (!success)
        goto cleanup;

    FLINT_ASSERT(Acontent->bits == ABbits);
    nmod_mpolyu_mul_mpoly(Gu, Gbar, Acontent, uctx);
    nmod_mpoly_from_mpolyu_perm_inflate(G, Gbits, ctx, Gu, uctx,
                                         zinfo->perm, I->Gmin_exp, I->Gstride);
    success = 1;

cleanup:

    nmod_mpolyu_clear(Abar, uctx);
    nmod_mpolyu_clear(Bbar, uctx);
    nmod_mpolyu_clear(Gbar, uctx);
    nmod_mpoly_clear(Acontent, uctx);
    nmod_mpoly_clear(Bcontent, uctx);

    nmod_mpolyu_clear(Au, uctx);
    nmod_mpolyu_clear(Bu, uctx);
    nmod_mpolyu_clear(Gu, uctx);
    nmod_mpoly_ctx_clear(uctx);

    mpoly_zipinfo_clear(zinfo);

    flint_randclear(randstate);

    return success;
}


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

    nmod_mpoly_to_mpolyn_perm_deflate(arg->Pn, arg->nctx, arg->P, arg->ctx,
           arg->perm, arg->shift, arg->stride, arg->handles, arg->num_handles);
}

static int _try_brown(
    nmod_mpoly_t G,
    flint_bitcnt_t Gbits,
    const nmod_mpoly_t A,
    const nmod_mpoly_t B,
    mpoly_gcd_info_t I,
    const nmod_mpoly_ctx_t ctx,
    const thread_pool_handle * handles,
    slong num_handles)
{
    int success;
    slong m = I->mvars;
    flint_bitcnt_t ABbits;
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

    ABbits = FLINT_MAX(A->bits, B->bits);

    nmod_mpoly_ctx_init(nctx, m, ORD_LEX, ctx->ffinfo->mod.n);
    nmod_poly_stack_init(Sp, ABbits, nctx);
    nmod_mpolyn_init(An, ABbits, nctx);
    nmod_mpolyn_init(Bn, ABbits, nctx);
    nmod_mpolyn_init(Gn, ABbits, nctx);
    nmod_mpolyn_init(Abarn, ABbits, nctx);
    nmod_mpolyn_init(Bbarn, ABbits, nctx);

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

        thread_pool_wake(global_thread_pool, handles[s], _worker_convertn, arg);

        nmod_mpoly_to_mpolyn_perm_deflate(An, nctx, A, ctx,
                                   I->brown_perm, I->Amin_exp, I->Gstride,
                                                               handles + 0, s);

        thread_pool_wait(global_thread_pool, handles[s]);
    }
    else
    {
        nmod_mpoly_to_mpolyn_perm_deflate(An, nctx, A, ctx,
                              I->brown_perm, I->Amin_exp, I->Gstride, NULL, 0);
        nmod_mpoly_to_mpolyn_perm_deflate(Bn, nctx, B, ctx,
                              I->brown_perm, I->Bmin_exp, I->Gstride, NULL, 0);
    }

    FLINT_ASSERT(An->bits == ABbits);
    FLINT_ASSERT(Bn->bits == ABbits);
    FLINT_ASSERT(An->length > 1);
    FLINT_ASSERT(Bn->length > 1);

    success = (num_handles > 0)
        ? nmod_mpolyn_gcd_brown_smprime_threaded(Gn, Abarn, Bbarn, An, Bn,
                                          m - 1, nctx, I, handles, num_handles)
        : nmod_mpolyn_gcd_brown_smprime(Gn, Abarn, Bbarn, An, Bn,
                                                           m - 1, nctx, I, Sp);

    if (!success)
    {
        nmod_mpoly_to_mpolyn_perm_deflate(An, nctx, A, ctx,
                              I->brown_perm, I->Amin_exp, I->Gstride, NULL, 0);
        nmod_mpoly_to_mpolyn_perm_deflate(Bn, nctx, B, ctx,
                              I->brown_perm, I->Bmin_exp, I->Gstride, NULL, 0);
        success = nmod_mpolyn_gcd_brown_lgprime(Gn, Abarn, Bbarn, An, Bn,
                                                                  m - 1, nctx);
    }

    if (!success)
        goto cleanup;

    nmod_mpoly_from_mpolyn_perm_inflate(G, Gbits, ctx, Gn, nctx,
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
    Assume when B is converted to univar format, its length would be one.
    Gcd is gcd of coefficients of univar(A) and B (modulo some shifts).
*/
static int _try_missing_var(
    nmod_mpoly_t G,
    flint_bitcnt_t Gbits,
    slong var,
    const nmod_mpoly_t A,
    ulong Ashift,
    const nmod_mpoly_t B,
    ulong Bshift,
    const nmod_mpoly_ctx_t ctx)
{
    int success;
    slong i;
    nmod_mpoly_t tG;
    nmod_mpoly_univar_t Ax;

    nmod_mpoly_init(tG, ctx);
    nmod_mpoly_univar_init(Ax, ctx);

    success = nmod_mpoly_to_univar(Ax, A, var, ctx);
    if (!success)
        goto cleanup;

    FLINT_ASSERT(Ax->length > 0);
    success = _nmod_mpoly_gcd(tG, Gbits, B, Ax->coeffs + 0, ctx, NULL, 0);
    if (!success)
        goto cleanup;

    for (i = 1; i < Ax->length; i++)
    {
        success = _nmod_mpoly_gcd(tG, Gbits, tG, Ax->coeffs + i, ctx, NULL, 0);
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

/*
    Test if B divides A or A divides B
        TODO: incorporate deflation
*/
static int _try_divides(
    nmod_mpoly_t G,
    const nmod_mpoly_t A,
    int try_a,
    const nmod_mpoly_t B,
    int try_b,
    const nmod_mpoly_ctx_t ctx)
{
    int success;
    nmod_mpoly_t Q;

    nmod_mpoly_init(Q, ctx);

    if (try_b && nmod_mpoly_divides_threaded(Q, A, B, ctx, 1))
    {
        nmod_mpoly_set(G, B, ctx);
        success = 1;
        goto cleanup;
    }

    if (try_a && nmod_mpoly_divides_threaded(Q, B, A, ctx, 1))
    {
        nmod_mpoly_set(G, A, ctx);
        success = 1;
        goto cleanup;
    }

    success = 0;

cleanup:

    nmod_mpoly_clear(Q, ctx);

    return success;
}

/*
    The function must pack its answer into bits = Gbits <= FLINT_BITS
    Both A and B have to be packed into bits <= FLINT_BITS

    return is 1 for success, 0 for failure.
*/
int _nmod_mpoly_gcd(
    nmod_mpoly_t G,
    flint_bitcnt_t Gbits,
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

    if (A->length == 1)
    {
        return _nmod_mpoly_gcd_monomial(G, Gbits, B, A, ctx);
    }
    else if (B->length == 1)
    {
        return _nmod_mpoly_gcd_monomial(G, Gbits, A, B, ctx);
    }

    mpoly_gcd_info_init(I, nvars);

    /* entries of I are all now invalid */

    mpoly_gcd_info_limits(I->Amax_exp, I->Amin_exp, I->Alead_count,
                      I->Atail_count, A->exps, A->bits, A->length, ctx->minfo);
    mpoly_gcd_info_limits(I->Bmax_exp, I->Bmin_exp, I->Blead_count,
                      I->Btail_count, B->exps, B->bits, B->length, ctx->minfo);

    /* set ess(p) := p/term_content(p) */

    /* check if the cofactors could be monomials, i.e. ess(A) == ess(B) */
    if (A->length == B->length)
    {
        if (_nmod_mpoly_gcd_monomial_cofactors_sp(G, Gbits,
                                             A, I->Amax_exp, I->Amin_exp,
                                             B, I->Bmax_exp, I->Bmin_exp, ctx))
        {
            goto successful;
        }
    }

    mpoly_gcd_info_stride(I->Gstride,
            A->exps, A->bits, A->length, I->Amax_exp, I->Amin_exp,
            B->exps, B->bits, B->length, I->Bmax_exp, I->Bmin_exp, ctx->minfo);

    for (j = 0; j < nvars; j++)
    {
        if (I->Gstride[j] == 0)
        {
            FLINT_ASSERT(  I->Amax_exp[j] == I->Amin_exp[j]
                        || I->Bmax_exp[j] == I->Bmin_exp[j]);
        }
        else
        {
            FLINT_ASSERT((I->Amax_exp[j] - I->Amin_exp[j]) % I->Gstride[j] == 0);
            FLINT_ASSERT((I->Bmax_exp[j] - I->Bmin_exp[j]) % I->Gstride[j] == 0);
        }

        I->Adeflate_deg[j] = I->Gstride[j] == 0 ? 0
                           : (I->Amax_exp[j] - I->Amin_exp[j]) / I->Gstride[j];
        I->Bdeflate_deg[j] = I->Gstride[j] == 0 ? 0
                           : (I->Bmax_exp[j] - I->Bmin_exp[j]) / I->Gstride[j];
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

        nmod_mpoly_fit_length(G, 1, ctx);
        nmod_mpoly_fit_bits(G, Gbits, ctx);
        G->bits = Gbits;
        mpoly_set_monomial_ui(G->exps, I->Gmin_exp, Gbits, ctx->minfo);
        G->coeffs[0] = UWORD(1);
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
        success = _try_missing_var(G, Gbits, v_in_A_only,
                                             A, I->Amin_exp[v_in_A_only],
                                             B, I->Bmin_exp[v_in_A_only], ctx);
        goto cleanup;
    }
    if (v_in_B_only != -WORD(1))
    {
        success = _try_missing_var(G, Gbits, v_in_B_only,
                                             B, I->Bmin_exp[v_in_B_only],
                                             A, I->Amin_exp[v_in_B_only], ctx);
        goto cleanup;
    }

    /* all variable are now either
            missing from both ess(A) and ess(B), or
            present in both ess(A) and ess(B)
        and there are at least two in the latter case
    */

    mpoly_gcd_info_set_estimates_nmod_mpoly(I, A, B, ctx, handles, num_handles);
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
            {
                gcd_is_trivial = 0;
            }

            if (I->Adeflate_deg[j] != I->Gdeflate_deg_bound[j]
                || I->Amin_exp[j] > I->Bmin_exp[j])
            {
                try_a = 0;
            }

            if (I->Bdeflate_deg[j] != I->Gdeflate_deg_bound[j]
                || I->Bmin_exp[j] > I->Amin_exp[j])
            {
                try_b = 0;
            }
        }

        if (gcd_is_trivial)
            goto calculate_trivial_gcd;

        if ((try_a || try_b) && _try_divides(G, A, try_a, B, try_b, ctx))
            goto successful;
    }

    mpoly_gcd_info_measure_brown(I, A->length, B->length, ctx->minfo);
    mpoly_gcd_info_measure_zippel(I, A->length, B->length, ctx->minfo);

    if (I->zippel_time_est < I->brown_time_est)
    {
        if (_try_zippel(G, Gbits, A, B, I, ctx))
            goto successful;

        if (_try_brown(G, Gbits, A, B, I, ctx, handles, num_handles))
            goto successful;
    }
    else
    {
        if (_try_brown(G, Gbits, A, B, I, ctx, handles, num_handles))
            goto successful;

        if (_try_zippel(G, Gbits, A, B, I, ctx))
            goto successful;
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
    }

    return success;
}


int nmod_mpoly_gcd_threaded(
    nmod_mpoly_t G,
    const nmod_mpoly_t A,
    const nmod_mpoly_t B,
    const nmod_mpoly_ctx_t ctx,
    slong thread_limit)
{
    slong i;
    flint_bitcnt_t Gbits;
    int success;
    thread_pool_handle * handles;
    slong num_handles;

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
        /* usual gcd's go right down here */

        /* get workers */
        handles = NULL;
        num_handles = 0;
        if (global_thread_pool_initialized)
        {
            slong max_num_handles = thread_pool_get_size(global_thread_pool);
            max_num_handles = FLINT_MIN(thread_limit - 1, max_num_handles);
            if (max_num_handles > 0)
            {
                handles = (thread_pool_handle *) flint_malloc(
                                   max_num_handles*sizeof(thread_pool_handle));
                num_handles = thread_pool_request(global_thread_pool,
                                                     handles, max_num_handles);
            }
        }

        success = _nmod_mpoly_gcd(G, Gbits, A, B, ctx, handles, num_handles);

        for (i = 0; i < num_handles; i++)
        {
            thread_pool_give_back(global_thread_pool, handles[i]);
        }
        if (handles)
        {
            flint_free(handles);
        }

        return success;
    }

    if (A->length == 1)
    {
        return _nmod_mpoly_gcd_monomial(G, Gbits, B, A, ctx);
    }
    else if (B->length == 1)
    {
        return _nmod_mpoly_gcd_monomial(G, Gbits, A, B, ctx);
    }
    else if (_nmod_mpoly_gcd_monomial_cofactors(G, A, B, ctx))
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
        int useAnew = 0;
        int useBnew = 0;
        slong k;
        fmpz * Ashift, * Astride;
        fmpz * Bshift, * Bstride;
        fmpz * Gshift, * Gstride;
        nmod_mpoly_t Anew;
        nmod_mpoly_t Bnew;

        nmod_mpoly_init(Anew, ctx);
        nmod_mpoly_init(Bnew, ctx);

        if (A->bits > FLINT_BITS)
        {
            useAnew = nmod_mpoly_repack_bits(Anew, A, FLINT_BITS, ctx);
            if (!useAnew)
                goto could_not_repack;
        }

        if (B->bits > FLINT_BITS)
        {
            useBnew = nmod_mpoly_repack_bits(Bnew, B, FLINT_BITS, ctx);
            if (!useBnew)
                goto could_not_repack;
        }

        success = _nmod_mpoly_gcd(G, FLINT_BITS, useAnew ? Anew : A,
                                             useBnew ? Bnew : B, ctx, NULL, 0);
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

        success = _nmod_mpoly_gcd(G, FLINT_BITS, Anew, Bnew, ctx, NULL, 0);

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

int nmod_mpoly_gcd(
    nmod_mpoly_t G,
    const nmod_mpoly_t A,
    const nmod_mpoly_t B,
    const nmod_mpoly_ctx_t ctx)
{
    return nmod_mpoly_gcd_threaded(G, A, B, ctx, MPOLY_DEFAULT_THREAD_LIMIT);
}

