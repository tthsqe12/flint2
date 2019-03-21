/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "nmod_mpoly.h"
#include "thread_pool.h"
#include "profiler.h"

/*
    Try to set G to the gcd of A and B using Brown's alogrithm M.
    This function switches to a big primes version if needed.
    It should only really fail if the dense size of the inputs is too large.
*/
int nmod_mpoly_gcd_brown(nmod_mpoly_t G, const nmod_mpoly_t A,
                              const nmod_mpoly_t B, const nmod_mpoly_ctx_t ctx)
{
    int success;
    nmod_mpolyd_t Ad, Bd, Gd, Abar, Bbar;
    nmod_mpolyd_ctx_t dctx;
    slong nvars = ctx->minfo->nvars;
/*
timeit_t time;
*/
    success = 1;

    if (nmod_mpoly_is_zero(A, ctx)) {
        nmod_mpoly_set(G, B, ctx);
        goto cleanup_stage0;
    }
    if (nmod_mpoly_is_zero(B, ctx)) {
        nmod_mpoly_set(G, A, ctx);
        goto cleanup_stage0;
    }

    nmod_mpolyd_ctx_init(dctx, nvars);
    success = nmod_mpolyd_ctx_set_for_gcd(dctx, A, B, ctx);
    if (!success)
    {
        nmod_mpoly_zero(G, ctx);
        goto cleanup_stage1;
    }

    nmod_mpolyd_init(Ad, nvars);
    nmod_mpolyd_init(Bd, nvars);
    nmod_mpolyd_init(Gd, nvars);
    nmod_mpolyd_init(Abar, nvars);
    nmod_mpolyd_init(Bbar, nvars);
/*
timeit_start(time);
*/
    nmod_mpoly_convert_to_nmod_mpolyd(Ad, dctx, A, ctx);
    nmod_mpoly_convert_to_nmod_mpolyd(Bd, dctx, B, ctx);
/*
timeit_stop(time);
flint_printf(", (* convert in time *) %wd\n", time->wall);
*/
/*
timeit_start(time);
*/
    success = nmod_mpolyd_gcd_brown_smprime(Gd, Abar, Bbar, Ad, Bd, ctx->ffinfo);
/*
timeit_stop(time);
flint_printf(", (* d gcd time *) %wd\n", time->wall);
*/


    if (!success)
    {
printf("trying lgprime\n");
        nmod_mpoly_convert_to_nmod_mpolyd(Ad, dctx, A, ctx);
        nmod_mpoly_convert_to_nmod_mpolyd(Bd, dctx, B, ctx);
        success = nmod_mpolyd_gcd_brown_lgprime(Gd, Abar, Bbar, Ad, Bd, ctx->ffinfo);
        if (!success) {
            nmod_mpoly_zero(G, ctx);
        } else {
            nmod_mpoly_convert_from_nmod_mpolyd(G, ctx, Gd, dctx);
        }
    } else
    {
/*
timeit_start(time);
*/
        nmod_mpoly_convert_from_nmod_mpolyd(G, ctx, Gd, dctx);
/*
timeit_stop(time);
flint_printf(", (* convert out time *) %wd\n", time->wall);
*/
    }

    nmod_mpolyd_clear(Bbar);
    nmod_mpolyd_clear(Abar);
    nmod_mpolyd_clear(Gd);
    nmod_mpolyd_clear(Bd);
    nmod_mpolyd_clear(Ad);

cleanup_stage1:

    nmod_mpolyd_ctx_clear(dctx);

cleanup_stage0:

    if (!nmod_mpoly_is_zero(G, ctx))
        nmod_mpoly_make_monic(G, G, ctx);

    return success;
}





void nmod_mpoly_to_nmod_poly_keepbits(nmod_poly_t A, slong * Ashift,
               const nmod_mpoly_t B, slong var, const nmod_mpoly_ctx_t ctx);

void nmod_mpoly_from_nmod_poly_keepbits(nmod_mpoly_t A, const nmod_poly_t B,
                           slong Bshift, slong var, mp_bitcnt_t bits, const nmod_mpoly_ctx_t ctx);




/* if the coefficient doesn't exist, a new one is created (and set to zero) */
nmod_polydr_struct * _nmod_mpolyn_get_coeff(nmod_mpolyn_t A,
                                      ulong * pow, const nmod_mpoly_ctx_t uctx)
{
    slong i, j, a, b;
    nmod_polydr_struct * xk;
    slong N = mpoly_words_per_exp_sp(A->bits, uctx->minfo);
    int cmp;

    a = 0;
    b = A->length;

    if (b == 0 || mpoly_monomial_gt_nomask(pow, A->exps + N*0, N))
    {
        i = 0;
        goto create_new;
    }

    if (mpoly_monomial_equal(pow , A->exps + N*(b - 1), N))
    {
        return A->coeffs + b - 1;
    }

try_again:

    if (b - a < 4)
    {
        for (i = a; i < b && (cmp = mpoly_monomial_cmp_nomask(A->exps + N*i, pow, N)) >= 0; i++)
        {
            if (cmp == 0) 
            {
                return A->coeffs + i;
            }
        }
        goto create_new;
    }
    else
    {
        i = a + (b - a)/2;
        cmp = mpoly_monomial_cmp_nomask(A->exps + N*i, pow, N);
        if (cmp == 0) 
        {
            return A->coeffs + i;
        }
        else if (cmp > 0)
        {
            a = i;
        }
        else
        {
            b = i;
        }
        goto try_again;
    }

create_new:

    nmod_mpolyn_fit_length(A, A->length + 1, uctx);

    for (j = A->length; j > i; j--)
    {
        mpoly_monomial_set(A->exps + N*j, A->exps + N*(j - 1), N);
        nmod_polydr_swap(A->coeffs + j, A->coeffs + j - 1, uctx->ffinfo);
    }

    mpoly_monomial_set(A->exps + N*i, pow, N);
    A->length++;
    xk = A->coeffs + i;
    xk->length = 0;

    return xk;
}



/* if the coefficient doesn't exist, a new one is created (and set to zero) */
nmod_mpolyn_struct * _nmod_mpolyun_get_coeff(nmod_mpolyun_t A,
                                        ulong pow, const nmod_mpoly_ctx_t uctx)
{
    slong i, j, a, b;
    nmod_mpolyn_struct * xk;

    a = 0;
    b = A->length;

    if (b == 0 || pow > A->exps[0])
    {
        i = 0;
        goto create_new;
    }

    if (pow == A->exps[b - 1])
    {
        return A->coeffs + b - 1;
    }

try_again:

    if (b - a < 8)
    {
        for (i = a; i < b && A->exps[i] >= pow; i++)
        {
            if (A->exps[i] == pow) 
            {
                return A->coeffs + i;
            }
        }
        goto create_new;
    }
    else
    {
        i = a + (b - a)/2;
        if (A->exps[i] == pow) 
        {
            return A->coeffs + i;
        }
        else if (A->exps[i] > pow)
        {
            a = i;
        }
        else
        {
            b = i;
        }
        goto try_again;
    }

create_new:

    nmod_mpolyun_fit_length(A, A->length + 1, uctx);

    for (j = A->length; j > i; j--)
    {
        A->exps[j] = A->exps[j - 1];
        nmod_mpolyn_swap(A->coeffs + j, A->coeffs + j - 1);
    }
    
    A->length++;
    A->exps[i] = pow;
    xk = A->coeffs + i;
    xk->length = 0;
    FLINT_ASSERT(xk->bits == A->bits);

    return xk;
}



void nmod_mpoly_to_mpolyun_perm_deflate_bivar(nmod_mpolyun_t A, const nmod_mpoly_t B,
               const slong * perm, const ulong * shift, const ulong * stride,
                       const nmod_mpoly_ctx_t uctx, const nmod_mpoly_ctx_t ctx)
{
    slong j;
    slong m = 1;
    slong NB, NA = 1;
    slong p0 = perm[0], p1 = perm[1];
    ulong shift0 = shift[p0], shift1 = shift[p1];
    ulong stride0 = stride[p0], stride1 = stride[p1];
    ulong Bexp0, Bexp1;
    slong Boff0, Bshift0, Boff1, Bshift1;
    ulong mask;
    nmod_mpolyn_struct * Ac;

    FLINT_ASSERT(A->bits <= FLINT_BITS);
    FLINT_ASSERT(B->bits <= FLINT_BITS);
    FLINT_ASSERT(m == uctx->minfo->nvars);
    FLINT_ASSERT(NA == mpoly_words_per_exp_sp(A->bits, uctx->minfo));
    NB = mpoly_words_per_exp_sp(B->bits, ctx->minfo);

    mpoly_gen_offset_shift_sp(&Boff0, &Bshift0, p0, B->bits, ctx->minfo);
    mpoly_gen_offset_shift_sp(&Boff1, &Bshift1, p1, B->bits, ctx->minfo);
    mask = (-UWORD(1)) >> (FLINT_BITS - B->bits);

    for (j = 0; j < B->length; j++)
    {
        Bexp0 = ((B->exps + NB*j)[Boff0] >> Bshift0) & mask;
        Bexp1 = ((B->exps + NB*j)[Boff1] >> Bshift1) & mask;

        Ac = _nmod_mpolyun_get_coeff(A, stride1 == 1 ? (Bexp1 - shift1)
                                           : (Bexp1 - shift1) / stride1, uctx);
        FLINT_ASSERT(Ac->bits == A->bits);

        nmod_mpolyn_fit_length(Ac, 1, uctx);
        nmod_polydr_set_coeff_ui(Ac->coeffs + 0, stride0 == 1 ? (Bexp0 - shift0)
                                   : (Bexp0 - shift0) / stride0, B->coeffs[j], uctx->ffinfo);
        mpoly_monomial_zero(Ac->exps + NA*0, NA);
        Ac->length = 1;
    }
}

/*
    Convert B to A using the variable permutation perm.
    The uctx should be the context of the coefficients of A.
    The ctx should be the context of B.

    operation on each term:

    for 0 <= k <= m
        l = perm[k]
        Aexp[k] = (Bexp[l] - shift[l])/stride[l]

    the variable of index m - 1 in uctx is moved to dense storage in nmod_poly
*/
void nmod_mpoly_to_mpolyun_perm_deflate(nmod_mpolyun_t A, const nmod_mpoly_t B,
               const slong * perm, const ulong * shift, const ulong * stride,
                       const nmod_mpoly_ctx_t uctx, const nmod_mpoly_ctx_t ctx)
{
    slong j, k, l;
    slong NA = mpoly_words_per_exp_sp(A->bits, uctx->minfo);
    slong NB = mpoly_words_per_exp_sp(B->bits, ctx->minfo);
    slong n = ctx->minfo->nvars;
    slong m = uctx->minfo->nvars;
    ulong * Bexps;
    ulong * texp;
    slong * offs, * shifts;
    nmod_mpolyn_struct * Ac;
    nmod_polydr_struct * Acc;
    TMP_INIT;

    FLINT_ASSERT(m > 0);

    A->length = 0;

    if (m == 1)
    {
        nmod_mpoly_to_mpolyun_perm_deflate_bivar(A, B, perm, shift, stride, uctx, ctx);
        return;
    }

    if (m > 2)
    {
        nmod_mpolyu_t Au;
        nmod_mpolyu_init(Au, A->bits, uctx);
        nmod_mpoly_to_mpolyu_perm_deflate(Au, B, perm, shift, stride, uctx, ctx);
        nmod_mpolyu_cvtto_mpolyun(A, Au, m - 1, uctx);
        nmod_mpolyu_clear(Au, uctx);
        return;
    }

    TMP_START;
    Bexps = (ulong *) TMP_ALLOC(n*sizeof(ulong));
    texp = (ulong *) TMP_ALLOC(NA*sizeof(ulong));

    offs   = (slong *) TMP_ALLOC(m*sizeof(ulong));
    shifts = (slong *) TMP_ALLOC(m*sizeof(ulong));
    for (k = 0; k < m; k++)
    {
        mpoly_gen_offset_shift_sp(offs + k, shifts + k, k, A->bits, uctx->minfo);
    }

    for (j = 0; j < B->length; j++)
    {
        mpoly_get_monomial_ui(Bexps, B->exps + NB*j, B->bits, ctx->minfo);
        l = perm[m];
        Ac = _nmod_mpolyun_get_coeff(A, stride[l] == 1 ? (Bexps[l] - shift[l]) : (Bexps[l] - shift[l]) / stride[l], uctx);
        FLINT_ASSERT(Ac->bits == A->bits);

        mpoly_monomial_zero(texp, NA);
        for (k = 0; k + 1 < m; k++)
        {
            l = perm[k];
/*
            FLINT_ASSERT(stride[l] != UWORD(0));
            FLINT_ASSERT(((Bexps[l] - shift[l]) % stride[l]) == UWORD(0));
*/
            texp[offs[k]] += (stride[l] == 1 ? (Bexps[l] - shift[l]) : (Bexps[l] - shift[l]) / stride[l]) << shifts[k];
        }

        Acc = _nmod_mpolyn_get_coeff(Ac, texp, uctx);
        l = perm[m - 1];
        nmod_polydr_set_coeff_ui(Acc, stride[l] == 1 ? (Bexps[l] - shift[l]) : (Bexps[l] - shift[l]) / stride[l], B->coeffs[j], uctx->ffinfo);
    }
/*
printf("nmod_mpoly_to_mpolyun_perm_deflate returning\n");
printf("B: "); nmod_mpoly_print_pretty(B, NULL, ctx); printf("\n");
printf("A: "); nmod_mpolyun_print_pretty(A, NULL, uctx); printf("\n");
*/
    TMP_END;
}




/*
    Convert B to A using the variable permutation vector perm.
    A must be constructed with bits = Abits.

    operation on each term:

        for 0 <= l < n
            Aexp[l] = shift[l]

        for 0 <= k <= m
            l = perm[k]
            Aexp[l] += scale[l]*Bexp[k]
*/
void nmod_mpoly_from_mpolyun_perm_inflate(nmod_mpoly_t A, mp_bitcnt_t Abits,
                                                        const nmod_mpolyun_t B,
                const slong * perm, const ulong * shift, const ulong * stride,
                       const nmod_mpoly_ctx_t uctx, const nmod_mpoly_ctx_t ctx)
{
    slong n = ctx->minfo->nvars;
    slong m = uctx->minfo->nvars;
    slong i, j, h, k, l;
    slong NA, NB;
    slong Alen;
    mp_limb_t * Acoeff;
    ulong * Aexp;
    slong Aalloc;
    ulong * uexps;
    ulong * Aexps, * tAexp, * tAgexp;
    TMP_INIT;

    FLINT_ASSERT(B->length > 0);
    FLINT_ASSERT(Abits <= FLINT_BITS);
    FLINT_ASSERT(B->bits <= FLINT_BITS);
    FLINT_ASSERT(m + 1 <= n);
    TMP_START;

    uexps = (ulong *) TMP_ALLOC((m + 1)*sizeof(ulong));
    Aexps = (ulong *) TMP_ALLOC(n*sizeof(ulong));

    NA = mpoly_words_per_exp(Abits, ctx->minfo);
    NB = mpoly_words_per_exp(B->bits, uctx->minfo);

    tAexp = (ulong *) TMP_ALLOC(NA*sizeof(ulong));
    tAgexp = (ulong *) TMP_ALLOC(NA*sizeof(ulong));
    mpoly_gen_monomial_sp(tAgexp, perm[m - 1], Abits, ctx->minfo);
    for (i = 0; i < NA; i++)
    {
        tAgexp[i] *= stride[perm[m - 1]];
    }

    nmod_mpoly_fit_bits(A, Abits, ctx);
    A->bits = Abits;

    Acoeff = A->coeffs;
    Aexp = A->exps;
    Aalloc = A->alloc;
    Alen = 0;
    for (i = 0; i < B->length; i++)
    {
        nmod_mpolyn_struct * Bc = B->coeffs + i;
        FLINT_ASSERT(Bc->bits == B->bits);

        for (j = 0; j < Bc->length; j++)
        {
            mpoly_get_monomial_ui(uexps, Bc->exps + NB*j, Bc->bits, uctx->minfo);
            uexps[m] = B->exps[i];
            FLINT_ASSERT(uexps[m - 1] == 0);
            for (l = 0; l < n; l++)
            {
                Aexps[l] = shift[l];
            }
            for (k = 0; k <= m; k++)
            {
                l = perm[k];
                Aexps[l] += stride[l]*uexps[k];
            }

            mpoly_set_monomial_ui(tAexp, Aexps, Abits, ctx->minfo);

            l = perm[m - 1];
            h = (Bc->coeffs + j)->length;
            _nmod_mpoly_fit_length(&Acoeff, &Aexp, &Aalloc, Alen + h, NA);
            for (h--; h >= 0; h--)
            {
                mp_limb_t c = (Bc->coeffs + j)->coeffs[h];
                if (c == 0)
                {
                    continue;
                }
                mpoly_monomial_madd(Aexp + NA*Alen, tAexp, h, tAgexp, NA);
                Acoeff[Alen] = c;
                Alen++;
            }
        }
    }
    A->coeffs = Acoeff;
    A->exps = Aexp;
    A->alloc = Aalloc;
    _nmod_mpoly_set_length(A, Alen, ctx);
/*
printf("inflate returning\n");
printf("B: "); nmod_mpolyun_print_pretty(B,NULL,uctx); printf("\n");
printf("A: "); nmod_mpoly_print_pretty(A,NULL,ctx); printf("\n");
*/
    nmod_mpoly_sort_terms(A, ctx);

    TMP_END;
}




typedef struct
{
    nmod_mpolyun_struct * output;
    const nmod_mpoly_struct * input;
    const slong * perm;
    const ulong * shift;
    const ulong * stride;
    const nmod_mpoly_ctx_struct * uctx;
    const nmod_mpoly_ctx_struct * ctx;
}
_convertworker_arg_struct;

typedef _convertworker_arg_struct _convertworker_arg_t[1];

static void _convertworker(void * varg)
{
    _convertworker_arg_struct * arg = (_convertworker_arg_struct *) varg;

    nmod_mpoly_to_mpolyun_perm_deflate(arg->output, arg->input,
                      arg->perm, arg->shift, arg->stride, arg->uctx, arg->ctx);
}


int nmod_mpolyun_gcd_brown_lgprime(nmod_mpolyun_t G,
     nmod_mpolyun_t Abar, nmod_mpolyun_t Bbar, nmod_mpolyun_t A, nmod_mpolyun_t B,
                                         slong var, const nmod_mpoly_ctx_t ctx);

int nmod_mpoly_gcd_brownnew(nmod_mpoly_t G, const nmod_mpoly_t A,
                              const nmod_mpoly_t B, const nmod_mpoly_ctx_t ctx)
{
    int success;
    slong * perm;
    ulong * shift, * stride;
    slong i;
    mp_bitcnt_t new_bits;
    nmod_mpoly_ctx_t uctx;
    nmod_mpolyun_t An, Bn, Gn, Abarn, Bbarn;
/*
    timeit_t time;
*/

    if (nmod_mpoly_is_zero(A, ctx))
    {
        if (nmod_mpoly_is_zero(B, ctx))
        {
            nmod_mpoly_zero(G, ctx);
        }
        else
        {
            nmod_mpoly_make_monic(G, B, ctx);
        }
        return 1;
    }

    if (nmod_mpoly_is_zero(B, ctx))
    {
        nmod_mpoly_make_monic(G, A, ctx);
        return 1;
    }

    if (A->bits > FLINT_BITS || B->bits > FLINT_BITS)
    {
        return 0;
    }

    if (ctx->minfo->nvars == 1)
    {
        slong shiftA, shiftB;
        nmod_poly_t a, b, g;
        nmod_poly_init(a, ctx->ffinfo->mod.n);
        nmod_poly_init(b, ctx->ffinfo->mod.n);
        nmod_poly_init(g, ctx->ffinfo->mod.n);
        nmod_mpoly_to_nmod_poly_keepbits(a, &shiftA, A, 0, ctx);
        nmod_mpoly_to_nmod_poly_keepbits(b, &shiftB, B, 0, ctx);
        nmod_poly_gcd(g, a, b);
        nmod_mpoly_from_nmod_poly_keepbits(G, g, FLINT_MIN(shiftA, shiftB), 0, A->bits, ctx);
        nmod_poly_clear(a);
        nmod_poly_clear(b);
        nmod_poly_clear(g);
        return 1;
    }

    perm = (slong *) flint_malloc(ctx->minfo->nvars*sizeof(slong));
    shift = (ulong *) flint_malloc(ctx->minfo->nvars*sizeof(ulong));
    stride = (ulong *) flint_malloc(ctx->minfo->nvars*sizeof(ulong));
    for (i = 0; i < ctx->minfo->nvars; i++)
    {
        perm[i] = i + 1 < ctx->minfo->nvars ? i + 1 : 0;
        shift[i] = 0;
        stride[i] = 1;
    }

    new_bits = FLINT_MAX(A->bits, B->bits);

    nmod_mpoly_ctx_init(uctx, ctx->minfo->nvars - 1, ORD_LEX, ctx->ffinfo->mod.n);
    nmod_mpolyun_init(An, new_bits, uctx);
    nmod_mpolyun_init(Bn, new_bits, uctx);
    nmod_mpolyun_init(Gn, new_bits, uctx);
    nmod_mpolyun_init(Abarn, new_bits, uctx);
    nmod_mpolyun_init(Bbarn, new_bits, uctx);

    success = 0;
    {
        _convertworker_arg_t convertargs;
        thread_pool_handle handles[8];
/*
timeit_start(time);
*/
        if (global_thread_pool_initialized
                    && thread_pool_request(global_thread_pool, handles, 1) > 0)
        {
            convertargs->output = An;
            convertargs->input = A;
            convertargs->perm = perm;
            convertargs->shift = shift;
            convertargs->stride = stride;
            convertargs->uctx = uctx;
            convertargs->ctx = ctx;
            thread_pool_wake(global_thread_pool, handles[0], _convertworker, &convertargs);
            nmod_mpoly_to_mpolyun_perm_deflate(Bn, B, perm, shift, stride, uctx, ctx);
            thread_pool_wait(global_thread_pool, handles[0]);
            thread_pool_give_back(global_thread_pool, handles[0]);
        }
        else
        {
            nmod_mpoly_to_mpolyun_perm_deflate(An, A, perm, shift, stride, uctx, ctx);
            nmod_mpoly_to_mpolyun_perm_deflate(Bn, B, perm, shift, stride, uctx, ctx);
        }
/*
timeit_stop(time);
flint_printf(", (* convert in time *) %wd\n", time->wall);
*/

/*
timeit_start(time);
*/
        success = nmod_mpolyun_gcd_brown_smprime_threaded(Gn, Abarn, Bbarn, An, Bn, uctx->minfo->nvars - 1, uctx);
/*
timeit_stop(time);
flint_printf(", (*un gcd time *) %wd\n", time->wall);
*/

        if (success)
        {
/*
timeit_start(time);
*/
            nmod_mpoly_from_mpolyun_perm_inflate(G, new_bits, Gn, perm, shift, stride, uctx, ctx);
            nmod_mpoly_make_monic(G, G, ctx);
/*
timeit_stop(time);
flint_printf(", (*  convert out time *) %wd\n", time->wall);
*/
        }
    }

    if (!success)
    {
printf("trying lgprime\n");

        nmod_mpoly_to_mpolyun_perm_deflate(An, A, perm, shift, stride, uctx, ctx);
        nmod_mpoly_to_mpolyun_perm_deflate(Bn, B, perm, shift, stride, uctx, ctx);

        success = nmod_mpolyun_gcd_brown_lgprime(Gn, Abarn, Bbarn, An, Bn, uctx->minfo->nvars - 1, uctx);
        if (success)
        {
            nmod_mpoly_from_mpolyun_perm_inflate(G, new_bits, Gn, perm, shift, stride, uctx, ctx);
            nmod_mpoly_make_monic(G, G, ctx);
            success = 1;
        }

    }
    nmod_mpolyun_clear(An, uctx);
    nmod_mpolyun_clear(Bn, uctx);
    nmod_mpolyun_clear(Gn, uctx);
    nmod_mpolyun_clear(Abarn, uctx);
    nmod_mpolyun_clear(Bbarn, uctx);
    nmod_mpoly_ctx_clear(uctx);

    flint_free(perm);
    flint_free(shift);
    flint_free(stride);

    return success;
}

