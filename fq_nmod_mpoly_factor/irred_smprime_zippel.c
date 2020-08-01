/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "nmod_mpoly_factor.h"
#include "fq_nmod_mpoly_factor.h"
#include "ui_factor.h"

n_polyun_term_struct * n_polyun_get_term(n_polyun_t A, ulong k);

fq_nmod_mpoly_struct * _fq_nmod_mpolyu_get_coeff(fq_nmod_mpolyu_t A,
                                     ulong pow, const fq_nmod_mpoly_ctx_t uctx);

void fq_nmod_poly_product_roots(
    fq_nmod_poly_t master,
    const fq_nmod_struct * monomials,
    slong mlength,
    const fq_nmod_ctx_t ctx);

void _fq_nmod_eval_to_bpoly(
    fq_nmod_bpoly_t E,
    const fq_nmod_mpoly_t A,
    const fq_nmod_poly_struct * alphabetas,
    const fq_nmod_mpoly_ctx_t ctx);

void _fq_nmod_mpoly_set_bpoly_var1_zero(
    fq_nmod_mpoly_t A,
    flint_bitcnt_t Abits,
    const fq_nmod_bpoly_t B,
    slong var,
    const fq_nmod_mpoly_ctx_t ctx);

#if 0
static void fq_nmod_mpoly_get_mpolyu2(
    fq_nmod_mpolyu_t A,
    const fq_nmod_mpoly_t B,
    slong var0,
    slong var1,
    const fq_nmod_mpoly_ctx_t ctx)
{
    slong i;
    fq_nmod_mpoly_struct * Ac;
    ulong * Bexps;
    slong NA, NB;

    FLINT_ASSERT(ctx->minfo->ord == ORD_LEX);

    FLINT_ASSERT(var0 < ctx->minfo->nvars);
    FLINT_ASSERT(var1 < ctx->minfo->nvars);

    Bexps = FLINT_ARRAY_ALLOC(ctx->minfo->nvars, ulong);

    NA = mpoly_words_per_exp_sp(A->bits, ctx->minfo);
    NB = mpoly_words_per_exp(B->bits, ctx->minfo);

    A->length = 0;
    for (i = 0; i < B->length; i++)
    {
        mpoly_get_monomial_ui(Bexps, B->exps + NB*i, B->bits, ctx->minfo);
        FLINT_ASSERT(FLINT_BIT_COUNT(Bexps[var0]) < FLINT_BITS/2);
        FLINT_ASSERT(FLINT_BIT_COUNT(Bexps[var1]) < FLINT_BITS/2);
        Ac = _fq_nmod_mpolyu_get_coeff(A, pack_exp2(Bexps[var0], Bexps[var1]), ctx);
        FLINT_ASSERT(Ac->bits == A->bits);
        fq_nmod_mpoly_fit_length(Ac, Ac->length + 1, ctx);
        fq_nmod_set(Ac->coeffs + Ac->length, B->coeffs + i, ctx->fqctx);
        Bexps[var0] = 0;
        Bexps[var1] = 0;
        mpoly_set_monomial_ui(Ac->exps + NA*Ac->length, Bexps, A->bits, ctx->minfo);
        Ac->length++;
    }

    flint_free(Bexps);
/*
printf("nmod_mpoly_get_mpolyu3 returning: "); nmod_mpolyu3_print_pretty(A, ourvars[var0], ourvars[var1], ourvars[var2], ourvars, ctx); printf("\n");
*/
    FLINT_ASSERT(fq_nmod_mpolyu_is_canonical(A, ctx));
}
#endif

/*
    return 
        -1: singular vandermonde matrix
        0:  inconsistent system
        1:  success
*/
int fq_nmod_zip_find_coeffs_new(
    fq_nmod_struct * coeffs,             /* length mlength */
    const fq_nmod_struct * monomials,    /* length mlength */
    slong mlength,
    const fq_nmod_struct * evals,        /* length elength */
    slong elength,    
    fq_nmod_poly_t master,
    const fq_nmod_ctx_t ctx)
{
    int success;
    slong i, j;
    fq_nmod_t V, V0, T, S, r, p0;

    FLINT_ASSERT(elength >= mlength);

    fq_nmod_init(V, ctx);
    fq_nmod_init(V0, ctx);
    fq_nmod_init(T, ctx);
    fq_nmod_init(S, ctx);
    fq_nmod_init(r, ctx);
    fq_nmod_init(p0, ctx);

    fq_nmod_poly_product_roots(master, monomials, mlength, ctx);

    for (i = 0; i < mlength; i++)
    {
        /* coeffs[i] is (coeffs(P).values)/P(roots[i]) =: V/S
            where P(x) = master(x)/(x-roots[i])     */
        fq_nmod_zero(V0, ctx);
        fq_nmod_zero(T, ctx);
        fq_nmod_zero(S, ctx);
        fq_nmod_set(r, monomials + i, ctx);
        for (j = mlength; j > 0; j--)
        {
            fq_nmod_mul(T, r, T, ctx);
            fq_nmod_add(T, T, master->coeffs + j, ctx);

            fq_nmod_mul(S, r, S, ctx);
            fq_nmod_add(S, S, T, ctx);

            fq_nmod_mul(p0, evals + j - 1, T, ctx);
            fq_nmod_add(V0, V0, p0, ctx);
        }
        /* roots[i] should be a root of master */
#if WANT_ASSERT
        fq_nmod_mul(p0, r, T, ctx);
        fq_nmod_add(p0, p0, master->coeffs + 0, ctx);
        FLINT_ASSERT(fq_nmod_is_zero(p0, ctx));
#endif
        fq_nmod_set(V, V0, ctx);
        fq_nmod_mul(S, S, r, ctx);
        if (fq_nmod_is_zero(S, ctx))
        {
            success = -1;
            goto cleanup;
        }

        fq_nmod_inv(p0, S, ctx);
        fq_nmod_mul(coeffs + i, V, p0, ctx);
    }

    /* use the coefficients of master as temp work space */
    for (i = 0; i < mlength; i++)
    {
        fq_nmod_pow_ui(master->coeffs + i, monomials + i, mlength, ctx);
    }

    /* check that the remaining points match */
    for (i = mlength; i < elength; i++)
    {
        fq_nmod_zero(V0, ctx);
        fq_nmod_zero(S, ctx);
        for (j = 0; j < mlength; j++)
        {
            fq_nmod_mul(master->coeffs + j, master->coeffs + j, monomials + j, ctx);
            fq_nmod_mul(p0, coeffs + j, master->coeffs + j, ctx);
            fq_nmod_add(V0, V0, p0, ctx);
        }
        fq_nmod_set(V, V0, ctx);
        if (!fq_nmod_equal(V, evals + i, ctx))
        {
            success = 0;
            goto cleanup;
        }
    }

    success = 1;

cleanup:

    fq_nmod_init(V, ctx);
    fq_nmod_init(V0, ctx);
    fq_nmod_init(T, ctx);
    fq_nmod_init(S, ctx);
    fq_nmod_init(r, ctx);
    fq_nmod_init(p0, ctx);

    return success;
}

int fq_nmod_zip_find_coeffs_new2(
    fq_nmod_struct * coeffs,             /* length mlength */
    const fq_nmod_struct * monomials,    /* length mlength */
    slong mlength,
    const fq_nmod_struct * evals,        /* length elength */
    slong elength,
    const fq_nmod_struct * master,       /* length mlength + 1 */
    fq_nmod_struct * temp,               /* length mlength */
    const fq_nmod_ctx_t ctx)
{
    int success;
    slong i, j;
    fq_nmod_t V, V0, T, S, r, p0;

    FLINT_ASSERT(elength >= mlength);

    fq_nmod_init(V, ctx);
    fq_nmod_init(V0, ctx);
    fq_nmod_init(T, ctx);
    fq_nmod_init(S, ctx);
    fq_nmod_init(r, ctx);
    fq_nmod_init(p0, ctx);

    for (i = 0; i < mlength; i++)
    {
        /* coeffs[i] is (coeffs(P).values)/P(roots[i]) =: V/S
            where P(x) = master(x)/(x-roots[i])     */
        fq_nmod_zero(V0, ctx);
        fq_nmod_zero(T, ctx);
        fq_nmod_zero(S, ctx);
        fq_nmod_set(r, monomials + i, ctx);
        for (j = mlength; j > 0; j--)
        {
            fq_nmod_mul(T, r, T, ctx);
            fq_nmod_add(T, T, master + j, ctx);

            fq_nmod_mul(S, r, S, ctx);
            fq_nmod_add(S, S, T, ctx);

            fq_nmod_mul(p0, evals + j - 1, T, ctx);
            fq_nmod_add(V0, V0, p0, ctx);
        }
        /* roots[i] should be a root of master */
#if WANT_ASSERT
        fq_nmod_mul(p0, r, T, ctx);
        fq_nmod_add(p0, p0, master + 0, ctx);
        FLINT_ASSERT(fq_nmod_is_zero(p0, ctx));
#endif
        fq_nmod_set(V, V0, ctx);
        fq_nmod_mul(S, S, r, ctx);
        if (fq_nmod_is_zero(S, ctx))
        {
            success = -1;
            goto cleanup;
        }

        fq_nmod_inv(p0, S, ctx);
        fq_nmod_mul(coeffs + i, V, p0, ctx);
    }

    /* check that the remaining points match */
    for (j = 0; j < mlength; j++)
        fq_nmod_pow_ui(temp + j, monomials + j, mlength, ctx);

    for (i = mlength; i < elength; i++)
    {
        fq_nmod_zero(V0, ctx);
        fq_nmod_zero(S, ctx);
        for (j = 0; j < mlength; j++)
        {
            fq_nmod_mul(temp + j, temp + j, monomials + j, ctx);
            fq_nmod_mul(p0, coeffs + j, temp + j, ctx);
            fq_nmod_add(V0, V0, p0, ctx);
        }
        fq_nmod_set(V, V0, ctx);
        if (!fq_nmod_equal(V, evals + i, ctx))
        {
            success = 0;
            goto cleanup;
        }
    }

    success = 1;

cleanup:

    fq_nmod_init(V, ctx);
    fq_nmod_init(V0, ctx);
    fq_nmod_init(T, ctx);
    fq_nmod_init(S, ctx);
    fq_nmod_init(r, ctx);
    fq_nmod_init(p0, ctx);

    return success;
}

/*
    B vars: x0 x1 x2 x3 x4 xv           2 < v
    A vars: xv x0 x1 : 0 0 x2 x3 x4 0
*/

static void fq_nmod_mpoly_get_mpolyu3(
    fq_nmod_mpolyu_t A,
    const fq_nmod_mpoly_t B,
    slong var0,
    slong var1,
    slong var2,
    const fq_nmod_mpoly_ctx_t ctx)
{
    slong i;
    fq_nmod_mpoly_struct * Ac;
    ulong * Bexps;
    slong NA, NB;

    FLINT_ASSERT(ctx->minfo->ord == ORD_LEX);

    FLINT_ASSERT(var0 < ctx->minfo->nvars);
    FLINT_ASSERT(var1 < ctx->minfo->nvars);
    FLINT_ASSERT(var2 < ctx->minfo->nvars);

    Bexps = FLINT_ARRAY_ALLOC(ctx->minfo->nvars, ulong);

    NA = mpoly_words_per_exp_sp(A->bits, ctx->minfo);
    NB = mpoly_words_per_exp(B->bits, ctx->minfo);

    A->length = 0;
    for (i = 0; i < B->length; i++)
    {
        mpoly_get_monomial_ui(Bexps, B->exps + NB*i, B->bits, ctx->minfo);
        FLINT_ASSERT(FLINT_BIT_COUNT(Bexps[var0]) < FLINT_BITS/3);
        FLINT_ASSERT(FLINT_BIT_COUNT(Bexps[var1]) < FLINT_BITS/3);
        FLINT_ASSERT(FLINT_BIT_COUNT(Bexps[var2]) < FLINT_BITS/3);
        Ac = _fq_nmod_mpolyu_get_coeff(A, pack_exp3(Bexps[var0], Bexps[var1], Bexps[var2]), ctx);
        FLINT_ASSERT(Ac->bits == A->bits);
        fq_nmod_mpoly_fit_length(Ac, Ac->length + 1, ctx);
        fq_nmod_set(Ac->coeffs + Ac->length, B->coeffs + i, ctx->fqctx);
        Bexps[var0] = 0;
        Bexps[var1] = 0;
        Bexps[var2] = 0;
        mpoly_set_monomial_ui(Ac->exps + NA*Ac->length, Bexps, A->bits, ctx->minfo);
        Ac->length++;
    }

    flint_free(Bexps);
/*
printf("nmod_mpoly_get_mpolyu3 returning: "); nmod_mpolyu3_print_pretty(A, ourvars[var0], ourvars[var1], ourvars[var2], ourvars, ctx); printf("\n");
*/
    FLINT_ASSERT(fq_nmod_mpolyu_is_canonical(A, ctx));
}

void fq_nmod_mpoly_monomial_evals(
    fq_nmod_poly_t E,
    const fq_nmod_mpoly_t A,
    const fq_nmod_struct * alpha,
    const fq_nmod_mpoly_ctx_t ctx)
{
    slong i, j;
    slong offset, shift;
    slong N = mpoly_words_per_exp_sp(A->bits, ctx->minfo);
    slong nvars = ctx->minfo->nvars;
    ulong * Aexp;
    slong * LUToffset;
    ulong * LUTmask;
    fq_nmod_struct * LUTvalue;
    slong LUTlen;
    fq_nmod_t xpoweval;
    ulong * inputexpmask;
    TMP_INIT;

    FLINT_ASSERT(A->bits <= FLINT_BITS);

    TMP_START;

    LUToffset = (slong *) TMP_ALLOC(N*FLINT_BITS*sizeof(slong));
    LUTmask   = (ulong *) TMP_ALLOC(N*FLINT_BITS*sizeof(ulong));
    LUTvalue  = (fq_nmod_struct *) TMP_ALLOC(N*FLINT_BITS*sizeof(fq_nmod_struct));

    for (i = 0; i < N*FLINT_BITS; i++)
        fq_nmod_init(LUTvalue + i, ctx->fqctx);

    fq_nmod_init(xpoweval, ctx->fqctx);

    Aexp = A->exps;

    inputexpmask = (ulong *) TMP_ALLOC(N*sizeof(ulong));
    mpoly_monomial_zero(inputexpmask, N);
    for (i = 0; i < A->length; i++)
    {
        for (j = 0; j < N; j++)
        {
            inputexpmask[j] |= (Aexp + N*i)[j];
        }
    }

    LUTlen = 0;
    for (j = nvars - 1; j >= 0; j--)
    {
        mpoly_gen_offset_shift_sp(&offset, &shift, j, A->bits, ctx->minfo);

        fq_nmod_set(xpoweval, alpha + j, ctx->fqctx); /* xpoweval = alpha[i]^(2^i) */
        for (i = 0; i < A->bits; i++)
        {
            LUToffset[LUTlen] = offset;
            LUTmask[LUTlen] = (UWORD(1) << (shift + i));
            fq_nmod_set(LUTvalue + LUTlen, xpoweval, ctx->fqctx);
            if ((inputexpmask[offset] & LUTmask[LUTlen]) != 0)
            {
                LUTlen++;
            }
            fq_nmod_mul(xpoweval, xpoweval, xpoweval, ctx->fqctx);
        }
    }
    FLINT_ASSERT(LUTlen < N*FLINT_BITS);

    fq_nmod_poly_fit_length(E, A->length, ctx->fqctx);
    E->length = A->length;
    for (i = 0; i < A->length; i++)
    {
        fq_nmod_one(xpoweval, ctx->fqctx);
        for (j = 0; j < LUTlen; j++)
        {
            if (((Aexp + N*i)[LUToffset[j]] & LUTmask[j]) != 0)
            {
                fq_nmod_mul(xpoweval, xpoweval, LUTvalue + j, ctx->fqctx);
            }
        }
        fq_nmod_set(E->coeffs + i, xpoweval, ctx->fqctx);
    }

    for (i = 0; i < N*FLINT_BITS; i++)
        fq_nmod_clear(LUTvalue + i, ctx->fqctx);

    fq_nmod_clear(xpoweval, ctx->fqctx);

    TMP_END;
}

void _fq_nmod_mpoly_monomial_evals(
    fq_nmod_struct * E,
    const ulong * Aexps,
    flint_bitcnt_t Abits,
    slong Alen,
    const fq_nmod_struct * alpha,
    slong vstart,
    const fq_nmod_mpoly_ctx_t ctx)
{
    slong i, j;
    slong offset, shift;
    slong N = mpoly_words_per_exp_sp(Abits, ctx->minfo);
    slong nvars = ctx->minfo->nvars;
    slong * LUToffset;
    ulong * LUTmask;
    fq_nmod_struct * LUTvalue;
    slong LUTlen;
    fq_nmod_t xpoweval;
    ulong * inputexpmask;
    TMP_INIT;

    TMP_START;
/*
flint_printf("_nmod_mpoly_monomial_evals called Alen: %wd\n", Alen);
*/
    LUToffset = (slong *) TMP_ALLOC(N*FLINT_BITS*sizeof(slong));
    LUTmask   = (ulong *) TMP_ALLOC(N*FLINT_BITS*sizeof(ulong));
    LUTvalue  = (fq_nmod_struct *) TMP_ALLOC(N*FLINT_BITS*sizeof(fq_nmod_struct));

    for (i = 0; i < N*FLINT_BITS; i++)
        fq_nmod_init(LUTvalue + i, ctx->fqctx);
    fq_nmod_init(xpoweval, ctx->fqctx);

    inputexpmask = (ulong *) TMP_ALLOC(N*sizeof(ulong));
    mpoly_monomial_zero(inputexpmask, N);
    for (i = 0; i < Alen; i++)
    {
        for (j = 0; j < N; j++)
        {
            inputexpmask[j] |= (Aexps + N*i)[j];
        }
    }

    LUTlen = 0;
    for (j = nvars - 1; j >= vstart; j--)
    {
        mpoly_gen_offset_shift_sp(&offset, &shift, j, Abits, ctx->minfo);

        fq_nmod_set(xpoweval, alpha + j, ctx->fqctx); /* xpoweval = alpha[i]^(2^i) */
        for (i = 0; i < Abits; i++)
        {
            LUToffset[LUTlen] = offset;
            LUTmask[LUTlen] = (UWORD(1) << (shift + i));
            fq_nmod_set(LUTvalue + LUTlen, xpoweval, ctx->fqctx);
            if ((inputexpmask[offset] & LUTmask[LUTlen]) != 0)
                LUTlen++;
            fq_nmod_mul(xpoweval, xpoweval, xpoweval, ctx->fqctx);
        }
    }
    FLINT_ASSERT(LUTlen < N*FLINT_BITS);

    for (i = 0; i < Alen; i++)
    {
        fq_nmod_one(xpoweval, ctx->fqctx);
        for (j = 0; j < LUTlen; j++)
        {
            if (((Aexps + N*i)[LUToffset[j]] & LUTmask[j]) != 0)
            {
                fq_nmod_mul(xpoweval, xpoweval, LUTvalue + j, ctx->fqctx);
            }
        }
        fq_nmod_set(E + i, xpoweval, ctx->fqctx);
    }

    for (i = 0; i < N*FLINT_BITS; i++)
        fq_nmod_clear(LUTvalue + i, ctx->fqctx);
    fq_nmod_clear(xpoweval, ctx->fqctx);

    TMP_END;
}

void _fq_nmod_mpoly_monomial_evals_indirect(
    fq_nmod_struct * E,
    const ulong * Aexps,
    flint_bitcnt_t Abits,
    ulong * Aind,
    slong Alen,
    const fq_nmod_struct * alpha,
    slong vstart,
    slong vstop,
    const fq_nmod_mpoly_ctx_t ctx)
{
    slong i, j;
    slong offset, shift;
    slong N = mpoly_words_per_exp_sp(Abits, ctx->minfo);
    slong * LUToffset;
    ulong * LUTmask;
    fq_nmod_struct * LUTvalue;
    slong LUTlen;
    fq_nmod_t xpoweval;
    ulong * inputexpmask;
    const ulong * thisAexp;
    TMP_INIT;

    FLINT_ASSERT(0 <= vstart);
    FLINT_ASSERT(vstart < vstop);
    FLINT_ASSERT(vstop <= ctx->minfo->nvars);

    TMP_START;

    LUToffset = (slong *) TMP_ALLOC(N*FLINT_BITS*sizeof(slong));
    LUTmask   = (ulong *) TMP_ALLOC(N*FLINT_BITS*sizeof(ulong));
    LUTvalue  = (fq_nmod_struct *) TMP_ALLOC(N*FLINT_BITS*sizeof(fq_nmod_struct));

    for (i = 0; i < N*FLINT_BITS; i++)
        fq_nmod_init(LUTvalue + i, ctx->fqctx);
    fq_nmod_init(xpoweval, ctx->fqctx);

    inputexpmask = (ulong *) TMP_ALLOC(N*sizeof(ulong));
    mpoly_monomial_zero(inputexpmask, N);
    for (i = 0; i < Alen; i++)
    {
        thisAexp = Aexps + N*Aind[i];
        for (j = 0; j < N; j++)
            inputexpmask[j] |= thisAexp[j];
    }

    LUTlen = 0;
    for (j = vstop - 1; j >= vstart; j--)
    {
        mpoly_gen_offset_shift_sp(&offset, &shift, j, Abits, ctx->minfo);

        fq_nmod_set(xpoweval, alpha + j, ctx->fqctx); /* xpoweval = alpha[i]^(2^i) */
        for (i = 0; i < Abits; i++)
        {
            LUToffset[LUTlen] = offset;
            LUTmask[LUTlen] = (UWORD(1) << (shift + i));
            fq_nmod_set(LUTvalue + LUTlen, xpoweval, ctx->fqctx);
            if ((inputexpmask[offset] & LUTmask[LUTlen]) != 0)
                LUTlen++;
            fq_nmod_mul(xpoweval, xpoweval, xpoweval, ctx->fqctx);
        }
    }
    FLINT_ASSERT(LUTlen < N*FLINT_BITS);

    for (i = 0; i < Alen; i++)
    {
        thisAexp = Aexps + N*Aind[i];
        fq_nmod_one(xpoweval, ctx->fqctx);
        for (j = 0; j < LUTlen; j++)
        {
            if ((thisAexp[LUToffset[j]] & LUTmask[j]) != 0)
            {
                fq_nmod_mul(xpoweval, xpoweval, LUTvalue + j, ctx->fqctx);
            }
        }
        fq_nmod_set(E + i, xpoweval, ctx->fqctx);
    }

    for (i = 0; i < N*FLINT_BITS; i++)
        fq_nmod_clear(LUTvalue + i, ctx->fqctx);
    fq_nmod_clear(xpoweval, ctx->fqctx);

    TMP_END;
}

void fq_nmod_mpolyu_set_eval_helper(
    fq_nmod_polyun_t EH,
    const fq_nmod_mpolyu_t A,
    const fq_nmod_struct * alpha,
    const fq_nmod_mpoly_ctx_t ctx)
{
    slong i, j, n;
    fq_nmod_polyun_term_struct * EHterms;
    fq_nmod_struct * p, * q;

    fq_nmod_polyun_fit_length(EH, A->length, ctx->fqctx);
    EH->length = A->length;
    EHterms = EH->terms;

    for (i = 0; i < A->length; i++)
    {
        EHterms[i].exp = A->exps[i];
        n = A->coeffs[i].length;
        fq_nmod_poly_fit_length(EHterms[i].coeff, 3*n, ctx->fqctx);
        fq_nmod_mpoly_monomial_evals(EHterms[i].coeff, A->coeffs + i, alpha, ctx);
        FLINT_ASSERT(n == EHterms[i].coeff->length);
        p = EHterms[i].coeff->coeffs;
        q = A->coeffs[i].coeffs;
        for (j = n - 1; j >= 0; j--)
        {
            fq_nmod_t t1, t2;
            fq_nmod_init(t1, ctx->fqctx);
            fq_nmod_init(t2, ctx->fqctx);
            fq_nmod_set(t1, p + j, ctx->fqctx);
            fq_nmod_set(t2, q + j, ctx->fqctx);
            fq_nmod_set(p + 3*j + 0, t1, ctx->fqctx);
            fq_nmod_set(p + 3*j + 1, t2, ctx->fqctx);
            fq_nmod_set(p + 3*j + 2, t1, ctx->fqctx);
            fq_nmod_clear(t1, ctx->fqctx);
            fq_nmod_clear(t2, ctx->fqctx);
        }
    }
}

void fq_nmod_mpoly_delete_duplicate_terms(
    fq_nmod_mpoly_t A,
    const fq_nmod_mpoly_ctx_t ctx)
{
    slong i, j;
    slong N = mpoly_words_per_exp(A->bits, ctx->minfo);
    j = -1;
    for (i = 0; i < A->length; i++)
    {
        if (j >= 0 && mpoly_monomial_equal(A->exps + N*j, A->exps + N*i, N))
        {
            FLINT_ASSERT(fq_nmod_equal(A->coeffs + j, A->coeffs + i, ctx->fqctx));
            continue;
        }
        j++;
        fq_nmod_set(A->coeffs + j, A->coeffs + i, ctx->fqctx);
        mpoly_monomial_set(A->exps + N*j, A->exps + N*i, N);
    }
    j++;
    A->length = j;
}

slong fq_nmod_mpolyu_set_eval_helper_and_zip_form(
    fq_nmod_polyun_t EH,
    fq_nmod_mpolyu_t H,
    ulong deg,
    const fq_nmod_mpolyu_t B,
    const fq_nmod_struct * alpha,
    const fq_nmod_mpoly_ctx_t ctx)
{
    slong i, j, n;
    ulong x, y, z;
    fq_nmod_polyun_term_struct * EHterms;
    fq_nmod_struct * p, * q;
    fq_nmod_mpoly_struct * Hc;
    slong old_len, zip_length = 0;
    slong N = mpoly_words_per_exp(B->bits, ctx->minfo);
/*
flint_printf("nmod_mpolyu_set_eval_helper_and_zip_form called\n");
flint_printf("deg = %wu\n", deg);
flint_printf("B: ");
nmod_mpolyu3_print_pretty(B, "Y", "X", "Z", NULL, ctx);
flint_printf("\n");
*/
    fq_nmod_polyun_fit_length(EH, B->length, ctx->fqctx);
    EH->length = B->length;
    EHterms = EH->terms;

    H->length = 0;

    for (i = 0; i < B->length; i++)
    {
        EHterms[i].exp = B->exps[i];
        y = extract_exp(EHterms[i].exp, 2, 3);
        x = extract_exp(EHterms[i].exp, 1, 3);
        z = extract_exp(EHterms[i].exp, 0, 3);
        n = B->coeffs[i].length;
        fq_nmod_poly_fit_length(EHterms[i].coeff, 3*n, ctx->fqctx);
        fq_nmod_mpoly_monomial_evals(EHterms[i].coeff, B->coeffs + i, alpha, ctx);
        FLINT_ASSERT(n == EHterms[i].coeff->length);
        p = EHterms[i].coeff->coeffs;
        q = B->coeffs[i].coeffs;

        if (x < deg)
        {
            FLINT_ASSERT(y == 0 && "strange but ok");
            Hc = _fq_nmod_mpolyu_get_coeff(H, pack_exp3(0, x, z), ctx);
            fq_nmod_mpoly_fit_length(Hc, n, ctx);
            old_len = Hc->length;
            for (j = 0; j < n; j++)
                fq_nmod_set(Hc->coeffs + old_len + j, p + j, ctx->fqctx);
            mpoly_copy_monomials(Hc->exps + N*old_len, B->coeffs[i].exps, n, N);
            Hc->length += n;
            zip_length = FLINT_MAX(zip_length, Hc->length);
            if (old_len > 0)
            {
                FLINT_ASSERT(0 && "strange but ok");
                fq_nmod_mpoly_sort_terms(Hc, ctx);
                fq_nmod_mpoly_delete_duplicate_terms(Hc, ctx);
            }
        }

        for (j = n - 1; j >= 0; j--)
        {
            fq_nmod_t t1, t2;
            fq_nmod_init(t1, ctx->fqctx);
            fq_nmod_init(t2, ctx->fqctx);
            fq_nmod_set(t1, p + j, ctx->fqctx);
            fq_nmod_set(t2, q + j, ctx->fqctx);
            fq_nmod_set(p + 3*j + 0, t1, ctx->fqctx);
            fq_nmod_set(p + 3*j + 1, t2, ctx->fqctx);
            fq_nmod_set(p + 3*j + 2, t1, ctx->fqctx);
            fq_nmod_clear(t1, ctx->fqctx);
            fq_nmod_clear(t2, ctx->fqctx);
        }
    }

    return zip_length;
}


slong fq_nmod_mpoly_set_eval_helper_and_zip_form2(
    slong * deg1_, /* degree of B wrt main var 1 */
    fq_nmod_polyun_t EH,
    fq_nmod_polyun_t H,
    fq_nmod_polyun_t M,
    const fq_nmod_mpoly_t B,
    const fq_nmod_struct * betas,
    const fq_nmod_mpoly_ctx_t ctx)
{
    slong start, Bi, j, n;
    slong e0, e1, Hi, EHi;
    fq_nmod_polyun_term_struct * EHterms, * Hterms, * Mterms;
    fq_nmod_struct * p;
    slong zip_length = 0;
    flint_bitcnt_t Bbits = B->bits;
    const fq_nmod_struct * Bcoeffs = B->coeffs;
    slong Blen = B->length;
    const ulong * Bexps = B->exps;
    ulong mask = (-UWORD(1)) >> (FLINT_BITS - Bbits);
    slong N = mpoly_words_per_exp_sp(Bbits, ctx->minfo);
    slong off0, off1, shift0, shift1;
    slong deg0, deg1 = -1;
/*
flint_printf("nmod_mpoly_set_eval_helper_and_zip_form2 called\n");
flint_printf("B: "); nmod_mpoly_print_pretty(B, NULL, ctx); flint_printf("\n");
*/
    FLINT_ASSERT(Blen > 0);

    mpoly_gen_offset_shift_sp(&off0, &shift0, 0, Bbits, ctx->minfo);
    mpoly_gen_offset_shift_sp(&off1, &shift1, 1, Bbits, ctx->minfo);

    Bi = 0;
    deg0 = (Bexps[N*Bi + off0] >> shift0) & mask;

    EHi = 0;
    Hi = 0;

    while (Bi < Blen)
    {
        start = Bi;
        e0 = (Bexps[N*Bi + off0] >> shift0) & mask;
        e1 = (Bexps[N*Bi + off1] >> shift1) & mask;
        deg1 = FLINT_MAX(deg1, e1);
        while (1)
        {
            Bi++;
            if (Bi >= Blen)
                break;
            if (((Bexps[N*Bi + off0] >> shift0) & mask) != e0)
                break;
            if (((Bexps[N*Bi + off1] >> shift1) & mask) != e1)
                break;
        }

        n = Bi - start;

        fq_nmod_polyun_fit_length(EH, EHi + 1, ctx->fqctx);
        EHterms = EH->terms;
        EHterms[EHi].exp = pack_exp2(e0, e1);
        fq_nmod_poly_fit_length(EHterms[EHi].coeff, 3*n, ctx->fqctx);
        EHterms[EHi].coeff->length = n;
        p = EHterms[EHi].coeff->coeffs;
        EHi++;

        _fq_nmod_mpoly_monomial_evals(p, Bexps + N*start, Bbits, n, betas, 2, ctx);

        if (e0 < deg0)
        {
            fq_nmod_polyun_fit_length(H, Hi + 1, ctx->fqctx);
            fq_nmod_polyun_fit_length(M, Hi + 1, ctx->fqctx);
            Hterms = H->terms;
            Mterms = M->terms;
            Hterms[Hi].exp = pack_exp2(e0, e1);
            Mterms[Hi].exp = pack_exp2(e0, e1);
            fq_nmod_poly_fit_length(Hterms[Hi].coeff, n, ctx->fqctx);
            zip_length = FLINT_MAX(zip_length, n);
            Hterms[Hi].coeff->length = n;
            for (j = 0; j < n; j++)
                fq_nmod_set(Hterms[Hi].coeff->coeffs + j, p + j, ctx->fqctx);
            fq_nmod_poly_product_roots(Mterms[Hi].coeff, p, n, ctx->fqctx);
            Hi++;
        }

        for (j = n - 1; j >= 0; j--)
        {
            fq_nmod_t t1, t2;
            fq_nmod_init(t1, ctx->fqctx);
            fq_nmod_init(t2, ctx->fqctx);
            fq_nmod_set(t2, Bcoeffs + start + j, ctx->fqctx);
            fq_nmod_set(t1, p + j, ctx->fqctx);
            fq_nmod_set(p + 3*j + 0, t1, ctx->fqctx);
            fq_nmod_set(p + 3*j + 1, t2, ctx->fqctx);
            fq_nmod_set(p + 3*j + 2, t1, ctx->fqctx);
        }
    }

    EH->length = EHi;
    H->length = Hi;
    M->length = Hi;
/*
flint_printf("nmod_mpoly_set_eval_helper_and_zip_form2 returning deg1 = %wd\n", deg1);
*/
    *deg1_ = deg1;
    return zip_length;
}



static slong fq_nmod_mpoly_set_eval_helper_and_zip_form(
    ulong * deg_,       /* deg_X(B), output */
    fq_nmod_polyun_t EH,
    fq_nmod_mpolyu_t H,
    const fq_nmod_mpoly_t B,
    const fq_nmod_struct * alpha,
    slong yvar,         /* Y = gen(yvar) (X = gen(0), Z = gen(1))*/
    const fq_nmod_mpoly_ctx_t ctx)
{
    slong xvar = 0;
    slong zvar = 1;
    slong i, j, n;
    ulong y, x, z;
    slong yoff, xoff, zoff;
    slong yshift, xshift, zshift;
    fq_nmod_polyun_term_struct * EHterms;
    fq_nmod_struct * p;
    fq_nmod_mpoly_struct * Hc;
    slong old_len, zip_length = 0;
    flint_bitcnt_t bits = B->bits;
    slong Blen = B->length;
    const ulong * Bexps = B->exps;
    const fq_nmod_struct * Bcoeffs = B->coeffs;
    slong N = mpoly_words_per_exp(bits, ctx->minfo);
    ulong mask = (-UWORD(1)) >> (FLINT_BITS - bits);
    ulong * ind;
    n_polyun_t T;
    n_polyun_term_struct * Tt;
    ulong deg;
/*
flint_printf("nmod_mpoly_set_eval_helper_and_zip_form called\n");
flint_printf("deg = %wu\n", deg);
flint_printf("B: ");
nmod_mpoly_print_pretty(B, NULL, ctx);
flint_printf("\n");
*/
    FLINT_ASSERT(bits <= FLINT_BITS);
    FLINT_ASSERT(bits == B->bits);
    FLINT_ASSERT(bits == H->bits);
    FLINT_ASSERT(Blen > 0);

    n_polyun_init(T);

    mpoly_gen_offset_shift_sp(&yoff, &yshift, yvar, bits, ctx->minfo);
    mpoly_gen_offset_shift_sp(&xoff, &xshift, xvar, bits, ctx->minfo);
    mpoly_gen_offset_shift_sp(&zoff, &zshift, zvar, bits, ctx->minfo);

    deg = (Bexps[N*0 + xoff] >> xshift) & mask;
    FLINT_ASSERT(deg == fq_nmod_mpoly_degree_si(B, 0, ctx));

    /* TODO use a map here instead of this shit */
    for (i = 0; i < Blen; i++)
    {
        y = (Bexps[N*i + yoff] >> yshift) & mask;
        x = (Bexps[N*i + xoff] >> xshift) & mask;
        z = (Bexps[N*i + zoff] >> zshift) & mask;
        Tt = n_polyun_get_term(T, pack_exp3(y, x, z));
        FLINT_ASSERT(Tt->exp == pack_exp3(y, x, z));
        n_poly_fit_length(Tt->coeff, Tt->coeff->length + 1);
        Tt->coeff->coeffs[Tt->coeff->length] = i;
        Tt->coeff->length++;
    }

/*flint_printf("T:"); n_polyu3n_print_pretty(T, "Y", "X", "Z", "_"); flint_printf("\n");*/

    fq_nmod_polyun_fit_length(EH, T->length, ctx->fqctx);
    EH->length = T->length;
    EHterms = EH->terms;

    H->length = 0;

    for (i = 0; i < T->length; i++)
    {
        EHterms[i].exp = T->terms[i].exp;
        y = extract_exp(EHterms[i].exp, 2, 3);
        x = extract_exp(EHterms[i].exp, 1, 3);
        z = extract_exp(EHterms[i].exp, 0, 3);
        n = T->terms[i].coeff->length;
        fq_nmod_poly_fit_length(EHterms[i].coeff, 3*n, ctx->fqctx);
        EHterms[i].coeff->length = n;
        p = EHterms[i].coeff->coeffs;
        ind = T->terms[i].coeff->coeffs;
        _fq_nmod_mpoly_monomial_evals_indirect(p, Bexps, bits, ind, n, alpha,
                                                                2, yvar, ctx);
        if (x < deg)
        {
            FLINT_ASSERT(y == 0 && "strange but ok");
            Hc = _fq_nmod_mpolyu_get_coeff(H, pack_exp3(0, x, z), ctx);
            fq_nmod_mpoly_fit_length(Hc, n, ctx);
            old_len = Hc->length;
            for (j = 0; j < n; j++)
                fq_nmod_set(Hc->coeffs + old_len + j, p + j, ctx->fqctx);
            for (j = 0; j < n; j++)
            {
                mpoly_monomial_set(Hc->exps + N*(old_len + j),
                                   Bexps + N*ind[j], N);
            }
            Hc->length += n;
            zip_length = FLINT_MAX(zip_length, Hc->length);
            if (old_len > 0)
            {
                FLINT_ASSERT(0 && "strange but ok");
                fq_nmod_mpoly_sort_terms(Hc, ctx);
                fq_nmod_mpoly_delete_duplicate_terms(Hc, ctx);
            }
        }

        for (j = n - 1; j >= 0; j--)
        {
            fq_nmod_t t1, t2;
            fq_nmod_init(t1, ctx->fqctx);
            fq_nmod_init(t2, ctx->fqctx);
            fq_nmod_set(t1, p + j, ctx->fqctx);
            fq_nmod_set(t2, Bcoeffs + ind[j], ctx->fqctx);
            fq_nmod_set(p + 3*j + 0, t1, ctx->fqctx);
            fq_nmod_set(p + 3*j + 1, t2, ctx->fqctx);
            fq_nmod_set(p + 3*j + 2, t1, ctx->fqctx);
        }
    }

    n_polyun_clear(T);

    *deg_ = deg;
/*
flint_printf("nmod_mpoly_set_eval_helper_and_zip_form returning\n");
*/
    return zip_length;
}


static void fq_nmod_poly_eval_step(
    fq_nmod_t res,
    fq_nmod_poly_t A,
    const fq_nmod_ctx_t ctx)
{
    slong i, Alen = A->length;
    fq_nmod_struct * Acoeffs = A->coeffs;
    fq_nmod_t p;

    FLINT_ASSERT(3*Alen <= A->alloc);

    fq_nmod_init(p, ctx);

    fq_nmod_zero(res, ctx);
    for (i = 0; i < Alen; i++)
    {
        fq_nmod_mul(p, Acoeffs + 3*i + 0, Acoeffs + 3*i + 1, ctx);
        fq_nmod_add(res, res, p, ctx);
        fq_nmod_mul(Acoeffs + 3*i + 0, Acoeffs + 3*i + 0, Acoeffs + 3*i + 2, ctx);
    }

    fq_nmod_clear(p, ctx);
}

void fq_nmod_polyu_eval_step(
    fq_nmod_polyu_t E,
    fq_nmod_polyun_t A,
    const fq_nmod_ctx_t ctx)
{
    slong Ai, Ei;
    fq_nmod_polyun_term_struct * Aterms = A->terms;

    fq_nmod_polyu_fit_length(E, A->length, ctx);

    Ei = 0;
    for (Ai = 0; Ai < A->length; Ai++)
    {
        FLINT_ASSERT(Ei < E->alloc);
        E->exps[Ei] = Aterms[Ai].exp;
        fq_nmod_poly_eval_step(E->coeffs + Ei, Aterms[Ai].coeff, ctx);
        Ei += !fq_nmod_is_zero(E->coeffs + Ei, ctx);
    }
    E->length = Ei;
}

void fq_nmod_bpoly_eval_step(
    fq_nmod_bpoly_t E,
    fq_nmod_polyun_t A,
    fq_nmod_ctx_t ctx)
{
    slong Ai;
    fq_nmod_t c;
    ulong e0, e1;
    fq_nmod_polyun_term_struct * Aterms = A->terms;

    fq_nmod_init(c, ctx);

    fq_nmod_bpoly_zero(E, ctx);
    for (Ai = 0; Ai < A->length; Ai++)
    {
        fq_nmod_poly_eval_step(c, Aterms[Ai].coeff, ctx);
        e0 = extract_exp(Aterms[Ai].exp, 1, 2);
        e1 = extract_exp(Aterms[Ai].exp, 0, 2);
        if (fq_nmod_is_zero(c, ctx))
            continue;
        fq_nmod_bpoly_set_coeff(E, e0, e1, c, ctx);
    }

    fq_nmod_clear(c, ctx);
}


void fq_nmod_poly_eval_reset(fq_nmod_poly_t A, const fq_nmod_ctx_t ctx)
{
    slong i, Alen = A->length;
    fq_nmod_struct * Acoeffs = A->coeffs;

    FLINT_ASSERT(3*Alen <= A->alloc);

    for (i = 0; i < Alen; i++)
        fq_nmod_set(Acoeffs + 3*i + 0, Acoeffs + 3*i + 2, ctx);
}

void fq_nmod_polyun_eval_reset(fq_nmod_polyun_t A, const fq_nmod_ctx_t ctx)
{
    slong Ai;
    for (Ai = 0; Ai < A->length; Ai++)
        fq_nmod_poly_eval_reset(A->terms[Ai].coeff, ctx);
}


int fq_nmod_polyu2_add_zip_must_match(
    fq_nmod_polyun_t Z,
    const fq_nmod_bpoly_t A,
    slong cur_length,
    const fq_nmod_ctx_t ctx)
{
    slong i, Ai, ai;
    fq_nmod_polyun_term_struct * Zt = Z->terms;
    const fq_nmod_poly_struct * Acoeffs = A->coeffs;

    Ai = A->length - 1;
    ai = (Ai < 0) ? 0 : fq_nmod_poly_degree(A->coeffs + Ai, ctx);

    for (i = 0; i < Z->length; i++)
    {
        if (Ai >= 0 && Zt[i].exp == pack_exp2(Ai, ai))
        {
            /* Z present, A present */
            fq_nmod_set(Zt[i].coeff->coeffs + cur_length, Acoeffs[Ai].coeffs + ai, ctx);
            Zt[i].coeff->length = cur_length + 1;
            do {
                ai--;
            } while (ai >= 0 && fq_nmod_is_zero(Acoeffs[Ai].coeffs + ai, ctx));
            if (ai < 0)
            {
                do {
                    Ai--;
                } while (Ai >= 0 && Acoeffs[Ai].length == 0);
                if (Ai >= 0)
                    ai = fq_nmod_poly_degree(Acoeffs + Ai, ctx);
            }
        }
        else if (Ai < 0 || Zt[i].exp > pack_exp2(Ai, ai))
        {
            /* Z present, A missing */
            fq_nmod_zero(Zt[i].coeff->coeffs + cur_length, ctx);
            Zt[i].coeff->length = cur_length + 1;
        }
        else
        {
            /* Z missing, A present */
            return 0;
        }
    }

    return 1;
}

void fq_nmod_polyu3_add_zip_limit1(
    fq_nmod_polyun_t Z,
    const fq_nmod_polyun_t A,
    const ulong deg1,
    slong cur_length,
    slong fit_length,
    const fq_nmod_ctx_t ctx)
{
    const fq_nmod_polyun_term_struct * At = A->terms;
    const fq_nmod_polyun_term_struct * Ait;
    fq_nmod_polyun_term_struct * Zit;
    slong Ai, ai, Zi, j;

    Ai = -1;
    ai = -1;
    do {
        Ai++;
        Ait = At + Ai;
    } while (Ai < A->length && extract_exp(Ait->exp, 1, 3) >= deg1);
    if (Ai < A->length)
        ai = fq_nmod_poly_degree(Ait->coeff, ctx);

    Zi = 0;

    while (Ai < A->length && Zi < Z->length)
    {
        Zit = Z->terms + Zi;
        Ait = At + Ai;
        if (Ait->exp + ai > Zit->exp)
        {
            /* missing from Z */
            fq_nmod_polyun_fit_length(Z, Z->length + 1, ctx);
            for (j = Z->length; j > Zi; j--)
                fq_nmod_polyun_term_swap(Z->terms + j, Z->terms + j - 1);
            Z->length++;
            Zit = Z->terms + Zi;
            Zit->exp = Ait->exp + ai;
            fq_nmod_poly_fit_length(Zit->coeff, fit_length, ctx);
            Zit->coeff->length = cur_length;
            for (j = 0; j < cur_length; j++)
                fq_nmod_zero(Zit->coeff->coeffs + j, ctx);
            goto in_both;            
        }
        else if (Ait->exp + ai < Zit->exp)
        {
            /* missing from A */
            FLINT_ASSERT(cur_length + 1 <= Zit->coeff->alloc);
            fq_nmod_zero(Zit->coeff->coeffs + cur_length, ctx);
            Zit->coeff->length = cur_length + 1;
            Zi++;
        }
        else
        {
in_both:
            FLINT_ASSERT(cur_length == Zit->coeff->length);
            FLINT_ASSERT(cur_length + 1 <= Zit->coeff->alloc);
            fq_nmod_set(Zit->coeff->coeffs + cur_length, Ait->coeff->coeffs + ai, ctx);
            Zit->coeff->length = cur_length + 1;
            Zi++;
            do {
                ai--;
            } while (ai >= 0 && fq_nmod_is_zero(Ait->coeff->coeffs + ai, ctx));
            if (ai < 0)
            {
                do {
                    Ai++;
                    Ait = At + Ai;
                } while (Ai < A->length && extract_exp(Ait->exp, 1, 3) >= deg1);
                if (Ai < A->length)
                    ai = fq_nmod_poly_degree(Ait->coeff, ctx);
            }
        }
    }

    /* everything in A must be put on the end of Z */
    while (Ai < A->length)
    {
        Zi = Z->length;
        fq_nmod_polyun_fit_length(Z, Zi + A->length - Ai, ctx);
        Zit = Z->terms + Zi;
        Zit->exp = Ait->exp + ai;
        fq_nmod_poly_fit_length(Zit->coeff, fit_length, ctx);
        Zit->coeff->length = cur_length;
        for (j = 0; j < cur_length; j++)
            fq_nmod_zero(Zit->coeff->coeffs + j, ctx);
        Z->length = ++Zi;
        FLINT_ASSERT(cur_length == Zit->coeff->length);
        FLINT_ASSERT(cur_length + 1 <= Zit->coeff->alloc);
        fq_nmod_set(Zit->coeff->coeffs + cur_length, Ait->coeff->coeffs + ai, ctx);
        Zit->coeff->length = cur_length + 1;
        do {
            ai--;
        } while (ai >= 0 && fq_nmod_is_zero(Ait->coeff->coeffs + ai, ctx));
        if (ai < 0)
        {
            do {
                Ai++;
                Ait = At + Ai;
            } while (Ai < A->length && extract_exp(Ait->exp, 1, 3) >= deg1);
            if (Ai < A->length)
                ai = fq_nmod_poly_degree(Ait->coeff, ctx);
        }
    }

    /* everything in Z must have a zero appended */
    while (Zi < Z->length)
    {
        Zit = Z->terms + Zi;
        FLINT_ASSERT(cur_length == Zit->coeff->length);
        FLINT_ASSERT(cur_length + 1 <= Zit->coeff->alloc);
        fq_nmod_zero(Zit->coeff->coeffs + cur_length, ctx);
        Zit->coeff->length = cur_length + 1;
        Zi++;
    }

    for (Zi = 0; Zi < Z->length; Zi++)
    {
        FLINT_ASSERT(Z->terms[Zi].coeff->length == cur_length + 1);
    }
}

slong fq_nmod_mpolyu_find_term(const fq_nmod_mpolyu_t A, ulong e)
{
    slong i;
    for (i = 0; i < A->length; i++)
        if (A->exps[i] == e)
            return i;
    return -1;
}

/*
    for each Y^y*X^x*Z^z in B with x = deg,
        keep the Y^y*X^x*Z^z*poly(x1,...) in B
    for each Y^y*X^x*Z^z in Z,
        assert that x < deg
        if there is no Y^0*X^x*Z^y in H, fail
        find coefficients of poly using this entry in H
        output Y^y*X^x*Z^z*poly(x1,...) to A
    sort A

    return
        -1: singular vandermonde matrix encountered
        0:  inconsistent system encountered
        1:  success
*/
static int fq_nmod_mpoly_from_zip(
    fq_nmod_mpoly_t B,
    const fq_nmod_polyun_t Z,
    fq_nmod_mpolyu_t H,
    ulong deg,
    slong yvar,     /* Y = gen(yvar) */
    const fq_nmod_mpoly_ctx_t ctx)
{
    int success;
    slong Hi, Zi, Bi, i, j;
    slong xvar = 0;
    slong zvar = 1;
    ulong x, y, z;
    flint_bitcnt_t bits = B->bits;
    fq_nmod_struct * Bcoeffs;
    ulong * Bexps;
    slong N = mpoly_words_per_exp_sp(bits, ctx->minfo);
    ulong mask = (-UWORD(1)) >> (FLINT_BITS - bits);
    slong xoff, xshift, yoff, yshift, zoff, zshift;
    fq_nmod_polyun_term_struct * Zt = Z->terms;
    fq_nmod_mpoly_struct * Hc;
    slong Hlen = H->length;

    fq_nmod_polyun_t M;  /* temp */
    fq_nmod_polyun_init(M, ctx->fqctx);
/*
flint_printf("-----------------");
flint_printf("nmod_mpoly_from_zip called vars %wd, %wd, %wd\n", yvar, xvar, zvar);
flint_printf("Z: "); n_polyu3n_print_pretty(Z, "Y", "X", "Z", "_"); printf("\n");
flint_printf("H: "); nmod_mpolyu3_print_pretty(H, "Y", "X", "Z", NULL, ctx); printf("\n");
flint_printf("deg: %wd\n", deg);
*/
    FLINT_ASSERT(bits == H->bits);

    fq_nmod_polyun_fit_length(M, Hlen + 1, ctx->fqctx);
    for (i = 0; i <= Hlen; i++)
        M->terms[i].coeff->length = 0;

    mpoly_gen_offset_shift_sp(&yoff, &yshift, yvar, bits, ctx->minfo);
    mpoly_gen_offset_shift_sp(&xoff, &xshift, xvar, bits, ctx->minfo);
    mpoly_gen_offset_shift_sp(&zoff, &zshift, zvar, bits, ctx->minfo);

    /* x is most significant in ctx, so keeping the lc_x in B is easy */
    FLINT_ASSERT(xvar == 0);

    for (Bi = 0; Bi < B->length; Bi++)
    {
        x = (((B->exps + N*Bi)[xoff] >> xshift) & mask);
        FLINT_ASSERT(x <= deg);
        if (x != deg)
            break;
    }

    for (Zi = 0; Zi < Z->length; Zi++)
    {
        y = extract_exp(Zt[Zi].exp, 2, 3);
        x = extract_exp(Zt[Zi].exp, 1, 3);
        z = extract_exp(Zt[Zi].exp, 0, 3);
        FLINT_ASSERT(x < deg);
        Hi = fq_nmod_mpolyu_find_term(H, pack_exp3(0, x, z));
        if (Hi < 0)
            return 0;

        FLINT_ASSERT(Hi < Hlen);
        FLINT_ASSERT(H->exps[Hi] == pack_exp3(0, x, z));

        Hc = H->coeffs + Hi;
        FLINT_ASSERT(bits == Hc->bits);
        FLINT_ASSERT(Hc->length > 0);
        fq_nmod_mpoly_fit_length(B, Bi + Hc->length, ctx);
        Bcoeffs = B->coeffs;

        if (M->terms[Hi].coeff->length < 1)
        {
            fq_nmod_poly_product_roots(M->terms[Hi].coeff,
                                           Hc->coeffs, Hc->length, ctx->fqctx);
        }

        fq_nmod_poly_fit_length(M->terms[Hlen].coeff, Hc->length, ctx->fqctx);

        success = fq_nmod_zip_find_coeffs_new2(Bcoeffs + Bi, Hc->coeffs,
                    Hc->length, Zt[Zi].coeff->coeffs, Zt[Zi].coeff->length,
                    M->terms[Hi].coeff->coeffs, M->terms[Hlen].coeff->coeffs,
                                                             ctx->fqctx);
        if (success < 1)
            return success;

        Bexps = B->exps;
        for (j = Bi, i = 0; i < Hc->length; j++, i++)
        {
            if (fq_nmod_is_zero(Bcoeffs + j, ctx->fqctx))
                continue;
            fq_nmod_set(Bcoeffs + Bi, Bcoeffs + j, ctx->fqctx);
            FLINT_ASSERT(Bi < B->alloc);
            mpoly_monomial_set(Bexps + N*Bi, Hc->exps + N*i, N);
            (Bexps + N*Bi)[yoff] += y << yshift;
            Bi++;
        }
    }
    B->length = Bi;
    fq_nmod_mpoly_sort_terms(B, ctx);
    FLINT_ASSERT(fq_nmod_mpoly_is_canonical(B, ctx));
/*
flint_printf("nmod_mpoly_from_zip returning good\n");
flint_printf("B: "); nmod_mpoly_print_pretty(B, NULL, ctx); flint_printf("\n");
*/
    fq_nmod_polyun_clear(M, ctx->fqctx);

    return 1;
}


int fq_nmod_mpoly_hlift_zippel(
    slong m,
    fq_nmod_mpoly_struct * B,
    slong r,
    const fq_nmod_struct * alpha,
    const fq_nmod_mpoly_t A,
    const slong * degs,
    const fq_nmod_mpoly_ctx_t ctx,
    flint_rand_t state)
{
    int success;
    slong i;
    slong zip_fails_remaining;
    slong req_zip_images, cur_zip_image;
    fq_nmod_mpolyu_struct Au[1], * H;
    fq_nmod_polyun_struct Aeh[1], * Beh;
    fq_nmod_polyu_struct Aeval[1], * Beval;
    fq_nmod_polyun_struct * BBeval, * Z;
    fq_nmod_struct * beta;
    flint_bitcnt_t bits = A->bits;
    fq_nmod_mpoly_t T1, T2;
    ulong * Bdegs;
    const slong degs0 = degs[0];

    FLINT_ASSERT(m > 2);
    FLINT_ASSERT(r > 1);
    FLINT_ASSERT(bits <= FLINT_BITS);
/*
flint_printf("fq_nmod_mpoly_hlift_zippel called m = %wd\n", m);

flint_printf("alpha[m-1]: "); fq_nmod_print_pretty(alpha + m - 1, ctx->fqctx); flint_printf("\n");

for (i = 0; i < r; i++)
{
flint_printf("B[%wd]: ", i);
fq_nmod_mpoly_print_pretty(B + i, NULL, ctx);
flint_printf("\n");
}

flint_printf("A: ", i);
fq_nmod_mpoly_print_pretty(A, NULL, ctx);
flint_printf("\n");
*/

#if WANT_ASSERT
    {
        fq_nmod_mpoly_t T;
        fq_nmod_mpoly_init(T, ctx);
        fq_nmod_mpoly_one(T, ctx);
        for (i = 0; i < r; i++)
            fq_nmod_mpoly_mul(T, T, B + i, ctx);
        fq_nmod_mpoly_sub(T, A, T, ctx);
        fq_nmod_mpoly_evaluate_one_fq_nmod(T, T, m, alpha + m - 1, ctx);
        FLINT_ASSERT(fq_nmod_mpoly_is_zero(T, ctx));
        fq_nmod_mpoly_clear(T, ctx);
    }
#endif

    beta = FLINT_ARRAY_ALLOC(ctx->minfo->nvars, fq_nmod_struct);
    for (i = 0; i < ctx->minfo->nvars; i++)
        fq_nmod_init(beta + i, ctx->fqctx);

    Bdegs = FLINT_ARRAY_ALLOC(r, ulong);
    H = FLINT_ARRAY_ALLOC(r, fq_nmod_mpolyu_struct);
    Beh = FLINT_ARRAY_ALLOC(r, fq_nmod_polyun_struct);
    Beval = FLINT_ARRAY_ALLOC(r, fq_nmod_polyu_struct);
    BBeval = FLINT_ARRAY_ALLOC(r, fq_nmod_polyun_struct);
    Z = FLINT_ARRAY_ALLOC(r, fq_nmod_polyun_struct);

    fq_nmod_mpolyu_init(Au, bits, ctx);
    fq_nmod_polyun_init(Aeh, ctx->fqctx);
    fq_nmod_polyu_init(Aeval, ctx->fqctx);
    for (i = 0; i < r; i++)
    {
        fq_nmod_mpolyu_init(H + i, bits, ctx);
        fq_nmod_polyun_init(Beh + i, ctx->fqctx);
        fq_nmod_polyu_init(Beval + i, ctx->fqctx);
        fq_nmod_polyun_init(BBeval + i, ctx->fqctx);
        fq_nmod_polyun_init(Z + i, ctx->fqctx);
    }

    /* init done */

    for (i = 0; i < r; i++)
    {
        success = fq_nmod_mpoly_repack_bits_inplace(B + i, bits, ctx);
        if (!success)
            goto cleanup;
    }

    zip_fails_remaining = 3;

    fq_nmod_mpoly_get_mpolyu3(Au, A, m, 0, 1, ctx);

choose_betas:

    /* only beta[2], beta[3], ..., beta[m - 1] will be used */
    for (i = 0; i < ctx->minfo->nvars; i++)
    {
        fq_nmod_rand(beta + i, state, ctx->fqctx);
        if (fq_nmod_is_zero(beta + i, ctx->fqctx))
            fq_nmod_one(beta + i, ctx->fqctx);
/*
if (i >= 2 && i < m)
{
flint_printf("beta[%wd]: ", i);
fq_nmod_print_pretty(beta + i, ctx->fqctx);
flint_printf("\n");
}
*/
    }

    fq_nmod_mpolyu_set_eval_helper(Aeh, Au, beta, ctx);

    req_zip_images = 1;
    for (i = 0; i < r; i++)
    {
        slong this_zip_images;
        this_zip_images = fq_nmod_mpoly_set_eval_helper_and_zip_form(Bdegs + i,
                                          Beh + i, H + i, B + i, beta, m, ctx);
        req_zip_images = FLINT_MAX(req_zip_images, this_zip_images);
        FLINT_ASSERT(Bdegs[i] > 0);
    }

    cur_zip_image = 0;

next_zip_image:

    fq_nmod_polyu_eval_step(Aeval, Aeh, ctx->fqctx);
/*
flint_printf("Aeval: ");
fq_nmod_polyu3_print_pretty(Aeval, "Y", "X", "Z", ctx->fqctx);
flint_printf("\n");
*/
    for (i = 0; i < r; i++)
    {
        fq_nmod_polyu_eval_step(Beval + i, Beh + i, ctx->fqctx);
/*
flint_printf("Beval[%wd]: ", i);
fq_nmod_polyu3_print_pretty(Beval + i, "Y", "X", "Z", ctx->fqctx);
flint_printf("\n");
*/
    }


    success = fq_nmod_polyu3_hlift(r, BBeval, Aeval, Beval,
                                             alpha + m - 1, degs0, ctx->fqctx);
    if (success < 1)
    {
        if (--zip_fails_remaining >= 0)
            goto choose_betas;

flint_printf("fq_nmod_mpoly_hlift_zippel fail 1\n");

        success = 0;
        goto cleanup;
    }

    for (i = 0; i < r; i++)
    {
        fq_nmod_polyu3_add_zip_limit1(Z + i, BBeval + i, Bdegs[i],
                                    cur_zip_image, req_zip_images, ctx->fqctx);
    }

    cur_zip_image++;
    if (cur_zip_image < req_zip_images)
        goto next_zip_image;

    for (i = 0; i < r; i++)
    {
        success = fq_nmod_mpoly_from_zip(B + i, Z + i, H + i, Bdegs[i], m, ctx);
        if (success < 1)
        {
flint_printf("fq_nmod_mpoly_hlift_zippel fail 2\n");
            success = 0;
            goto cleanup;
        }
    }

    fq_nmod_mpoly_init3(T1, A->length, bits, ctx);
    fq_nmod_mpoly_init3(T2, A->length, bits, ctx);
    fq_nmod_mpoly_mul(T1, B + 0, B + 1, ctx);
    for (i = 2; i < r; i++)
    {
        fq_nmod_mpoly_mul(T2, T1, B + i, ctx);
        fq_nmod_mpoly_swap(T1, T2, ctx);
    }

    success = fq_nmod_mpoly_equal(T1, A, ctx);
    fq_nmod_mpoly_clear(T1, ctx);
    fq_nmod_mpoly_clear(T2, ctx);

cleanup:

    fq_nmod_mpolyu_clear(Au, ctx);
    fq_nmod_polyun_clear(Aeh, ctx->fqctx);
    fq_nmod_polyu_clear(Aeval, ctx->fqctx);
    for (i = 0; i < r; i++)
    {
        fq_nmod_mpolyu_clear(H + i, ctx);
        fq_nmod_polyun_clear(Beh + i, ctx->fqctx);
        fq_nmod_polyu_clear(Beval + i, ctx->fqctx);
        fq_nmod_polyun_clear(BBeval + i, ctx->fqctx);
        fq_nmod_polyun_clear(Z + i, ctx->fqctx);
    }

    for (i = 0; i < ctx->minfo->nvars; i++)
        fq_nmod_clear(beta + i, ctx->fqctx);
    flint_free(beta);

    flint_free(Bdegs);
    flint_free(H);
    flint_free(Beh);
    flint_free(Beval);
    flint_free(BBeval);
    flint_free(Z);

    return success;
}



int fq_nmod_mpoly_factor_irred_smprime_zippel(
    fq_nmod_mpolyv_t fac,
    const fq_nmod_mpoly_t A,
    const fq_nmod_mpoly_factor_t lcAfac,
    const fq_nmod_mpoly_t lcA,
    const fq_nmod_mpoly_ctx_t ctx,
    flint_rand_t state)
{
    int success;
    int alphas_tries_remaining, alphabetas_tries_remaining, alphabetas_length;
    const slong n = ctx->minfo->nvars - 1;
    slong i, j, k, r;
    fq_nmod_struct * alpha;
    fq_nmod_poly_struct * alphabetas;
    fq_nmod_mpoly_struct * Aevals;
    slong * degs, * degeval;
    fq_nmod_mpolyv_t tfac;
    fq_nmod_mpoly_t t, Acopy;
    fq_nmod_mpoly_struct * newA;
    fq_nmod_poly_t Abfc;
    fq_nmod_bpoly_t Ab;
    fq_nmod_tpoly_t Abfp;
    fq_nmod_mpoly_t m, mpow;
    fq_nmod_mpolyv_t new_lcs, lc_divs;
/*
flint_printf("fq_nmod_mpoly_factor_irred_smprime_zippel called\n");
flint_printf("     A: "); fq_nmod_mpoly_print_pretty(A, NULL, ctx); flint_printf("\n");
flint_printf("lcAfac: "); fq_nmod_mpoly_factor_print_pretty(lcAfac, NULL, ctx); flint_printf("\n");
flint_printf("   lcA: "); fq_nmod_mpoly_print_pretty(lcA, NULL, ctx); flint_printf("\n");
flint_printf("ctx:\n");
fq_nmod_ctx_print(ctx->fqctx);
*/
    FLINT_ASSERT(n > 1);
    FLINT_ASSERT(A->length > 1);
    FLINT_ASSERT(fq_nmod_is_one(A->coeffs + 0, ctx->fqctx));
    FLINT_ASSERT(A->bits <= FLINT_BITS);

    if (ctx->fqctx->modulus->length < n_clog(A->length, ctx->fqctx->modulus->mod.n))
        return 0;

    fq_nmod_mpoly_init(Acopy, ctx);
    fq_nmod_mpoly_init(m, ctx);
    fq_nmod_mpoly_init(mpow, ctx);

    fq_nmod_mpolyv_init(new_lcs, ctx);
    fq_nmod_mpolyv_init(lc_divs, ctx);

    fq_nmod_poly_init(Abfc, ctx->fqctx);
    fq_nmod_tpoly_init(Abfp, ctx->fqctx);
    fq_nmod_bpoly_init(Ab, ctx->fqctx);

    degs    = FLINT_ARRAY_ALLOC(n + 1, slong);
    degeval = FLINT_ARRAY_ALLOC(n + 1, slong);
	alpha   = FLINT_ARRAY_ALLOC(n, fq_nmod_struct);
    alphabetas = FLINT_ARRAY_ALLOC(n, fq_nmod_poly_struct);
    Aevals  = FLINT_ARRAY_ALLOC(n, fq_nmod_mpoly_struct);
	for (i = 0; i < n; i++)
    {
        fq_nmod_init(alpha + i, ctx->fqctx);
        fq_nmod_poly_init(alphabetas + i, ctx->fqctx);
		fq_nmod_mpoly_init(Aevals + i, ctx);
    }
    fq_nmod_mpolyv_init(tfac, ctx);
	fq_nmod_mpoly_init(t, ctx);

    /* init done */

    alphabetas_length = 2;
    alphas_tries_remaining = 10;
	fq_nmod_mpoly_degrees_si(degs, A, ctx);

next_alpha:

    if (--alphas_tries_remaining < 0)
	{
		success = 0;
        goto cleanup;
	}

    for (i = 0; i < n; i++)
    {
        fq_nmod_rand(alpha + i, state, ctx->fqctx);
/*
flint_printf("alpha[%wd]: ", i);
fq_nmod_print_pretty(alpha + i, ctx->fqctx);
flint_printf("\n");
*/
    }

    /* ensure degrees do not drop under evaluation */
	for (i = n - 1; i >= 0; i--)
	{
        fq_nmod_mpoly_evaluate_one_fq_nmod(Aevals + i,
                       i == n - 1 ? A : Aevals + i + 1, i + 1, alpha + i, ctx);
		fq_nmod_mpoly_degrees_si(degeval, Aevals + i, ctx);
		for (j = 0; j <= i; j++)
			if (degeval[j] != degs[j])
				goto next_alpha;
	}

    /* make sure univar is squarefree */
	fq_nmod_mpoly_derivative(t, Aevals + 0, 0, ctx);
	fq_nmod_mpoly_gcd(t, t, Aevals + 0, ctx);
	if (!fq_nmod_mpoly_is_one(t, ctx))
		goto next_alpha;

    alphabetas_tries_remaining = 2 + alphabetas_length;

next_alphabetas:

    if (--alphabetas_tries_remaining < 0)
    {
        if (++alphabetas_length > 5)
        {
            success = 0;
            goto cleanup;
        }
        goto next_alpha;
    }

    for (i = 0; i < n; i++)
    {
        fq_nmod_poly_fit_length(alphabetas + i, alphabetas_length, ctx->fqctx);
        fq_nmod_set(alphabetas[i].coeffs + 0, alpha + i, ctx->fqctx);
        for (j = 1; j < alphabetas_length; j++)
            fq_nmod_rand(alphabetas[i].coeffs + j, state, ctx->fqctx);
        alphabetas[i].length = alphabetas_length;
        _fq_nmod_poly_normalise(alphabetas + i, ctx->fqctx);
/*
flint_printf("alphabetas[%wd]: ", i);
fq_nmod_poly_print_pretty(alphabetas + i, "Y", ctx->fqctx);
flint_printf("\n");
*/
    }

    _fq_nmod_eval_to_bpoly(Ab, A, alphabetas, ctx);
    success = fq_nmod_bpoly_factor_smprime(Abfc, Abfp, Ab, 0, ctx->fqctx);
    if (!success)
    {
        FLINT_ASSERT(0 && "this should not happen");
        goto next_alpha;
    }

    r = Abfp->length;

    if (r < 2)
    {
        fq_nmod_mpolyv_fit_length(fac, 1, ctx);
        fac->length = 1;
        fq_nmod_mpoly_set(fac->coeffs + 0, A, ctx);
        success = 1;
        goto cleanup;
    }

    fq_nmod_mpolyv_fit_length(lc_divs, r, ctx);
    lc_divs->length = r;
    if (lcAfac->num > 0)
    {
        success = fq_nmod_mpoly_factor_lcc_wang(lc_divs->coeffs, lcAfac,
                                       Abfc, Abfp->coeffs, r, alphabetas, ctx);
        if (!success)
            goto next_alphabetas;
    }
    else
    {
        for (i = 0; i < r; i++)
            fq_nmod_mpoly_one(lc_divs->coeffs + i, ctx);
    }
/*
flint_printf("lc_divs:\n");
fq_nmod_mpolyv_print_pretty(lc_divs, NULL, ctx);
*/

    success = fq_nmod_mpoly_divides(m, lcA, lc_divs->coeffs + 0, ctx);
    FLINT_ASSERT(success);
    for (i = 1; i < r; i++)
    {
        success = fq_nmod_mpoly_divides(m, m, lc_divs->coeffs + i, ctx);
        FLINT_ASSERT(success);
    }
/*
flint_printf("m: ");
fq_nmod_mpoly_print_pretty(m, NULL, ctx);
flint_printf("\n");
*/
    fq_nmod_mpoly_pow_ui(mpow, m, r - 1, ctx);
    if (fq_nmod_mpoly_is_one(mpow, ctx))
    {
        newA = (fq_nmod_mpoly_struct *) A;
    }
    else
    {
        newA = Acopy;
        fq_nmod_mpoly_mul(newA, A, mpow, ctx);
    }

    if (newA->bits > FLINT_BITS)
    {
        success = 0;
        goto cleanup;
    }

    fq_nmod_mpoly_degrees_si(degs, newA, ctx);

    fq_nmod_mpoly_set(t, mpow, ctx);
    for (i = n - 1; i >= 0; i--)
    {
        fq_nmod_mpoly_evaluate_one_fq_nmod(t, mpow, i + 1, alpha + i, ctx);
        fq_nmod_mpoly_swap(t, mpow, ctx);
        fq_nmod_mpoly_mul(Aevals + i, Aevals + i, mpow, ctx);
    }

    fq_nmod_mpolyv_fit_length(new_lcs, (n + 1)*r, ctx);
    i = n;
    for (j = 0; j < r; j++)
    {
        fq_nmod_mpoly_mul(new_lcs->coeffs + i*r + j, lc_divs->coeffs + j, m, ctx);
    }
    for (i = n - 1; i >= 0; i--)
    {
        for (j = 0; j < r; j++)
        {
            fq_nmod_mpoly_evaluate_one_fq_nmod(new_lcs->coeffs + i*r + j,
                       new_lcs->coeffs + (i + 1)*r + j, i + 1, alpha + i, ctx);
        }
    }

    fq_nmod_mpolyv_fit_length(fac, r, ctx);
    fac->length = r;
    for (i = 0; i < r; i++)
    {
        fq_nmod_t q;
        fq_nmod_init(q, ctx->fqctx);
        FLINT_ASSERT(fq_nmod_mpoly_is_fq_nmod(new_lcs->coeffs + 0*r + i, ctx));
        FLINT_ASSERT(fq_nmod_mpoly_length(new_lcs->coeffs + 0*r + i, ctx) == 1);
        _fq_nmod_mpoly_set_bpoly_var1_zero(fac->coeffs + i, newA->bits, Abfp->coeffs + i, 0, ctx);
        FLINT_ASSERT(fac->coeffs[i].length > 0);
        fq_nmod_inv(q, fac->coeffs[i].coeffs + 0, ctx->fqctx);
        fq_nmod_mul(q, q, new_lcs->coeffs[0*r + i].coeffs + 0, ctx->fqctx);
        fq_nmod_mpoly_scalar_mul_fq_nmod(fac->coeffs + i, fac->coeffs + i, q, ctx);
        fq_nmod_clear(q, ctx->fqctx);
    }

    fq_nmod_mpolyv_fit_length(tfac, r, ctx);
    tfac->length = r;
    for (k = 1; k <= n; k++)
    {
        for (i = 0; i < r; i++)
        {
            _fq_nmod_mpoly_set_lead0(tfac->coeffs + i, fac->coeffs + i,
                                               new_lcs->coeffs + k*r + i, ctx);
        }

        if (k > 2)
        {
/*
flint_printf("starting zippel lift k = %wd\n", k);
*/
            success = fq_nmod_mpoly_hlift_zippel(k, tfac->coeffs, r, alpha,
                                  k < n ? Aevals + k : newA, degs, ctx, state);
/*
flint_printf("finished zippel lift k = %wd, success = %d\n", k, success);
*/
        }
        else
        {
            success = fq_nmod_mpoly_hlift(k, tfac->coeffs, r, alpha,
                                         k < n ? Aevals + k : newA, degs, ctx);
        }

        if (!success)
            goto next_alphabetas;

        fq_nmod_mpolyv_swap(tfac, fac, ctx);
    }

    if (!fq_nmod_mpoly_is_fq_nmod(m, ctx))
    {
        fq_nmod_mpoly_univar_t u;
        fq_nmod_mpoly_univar_init(u, ctx);
        for (i = 0; i < r; i++)
        {
            fq_nmod_mpoly_to_univar(u, fac->coeffs + i, 0, ctx);
            success = fq_nmod_mpoly_univar_content_mpoly(t, u, ctx);
            if (!success)
            {
                fq_nmod_mpoly_univar_clear(u, ctx);
                goto cleanup;
            }
            success = fq_nmod_mpoly_divides(fac->coeffs + i,
                                            fac->coeffs + i, t, ctx);
            FLINT_ASSERT(success);
        }
        fq_nmod_mpoly_univar_clear(u, ctx);
    }

    for (i = 0; i < r; i++)
        fq_nmod_mpoly_make_monic(fac->coeffs + i, fac->coeffs + i, ctx);

    success = 1;

cleanup:

    fq_nmod_mpolyv_clear(new_lcs, ctx);
    fq_nmod_mpolyv_clear(lc_divs, ctx);

    fq_nmod_poly_clear(Abfc, ctx->fqctx);
    fq_nmod_tpoly_clear(Abfp, ctx->fqctx);

	for (i = 0; i < n; i++)
    {
		fq_nmod_mpoly_clear(Aevals + i, ctx);
        fq_nmod_poly_clear(alphabetas + i, ctx->fqctx);
    }
    flint_free(alphabetas);
    flint_free(alpha);
    flint_free(Aevals);
    flint_free(degs);
    flint_free(degeval);

    fq_nmod_mpolyv_clear(tfac, ctx);
    fq_nmod_mpoly_clear(t, ctx);

    fq_nmod_mpoly_clear(Acopy, ctx);
    fq_nmod_mpoly_clear(m, ctx);
    fq_nmod_mpoly_clear(mpow, ctx);

#if WANT_ASSERT
    if (success)
    {
        fq_nmod_mpoly_t prod;
        fq_nmod_mpoly_init(prod, ctx);
        fq_nmod_mpoly_one(prod, ctx);
        for (i = 0; i < fac->length; i++)
            fq_nmod_mpoly_mul(prod, prod, fac->coeffs + i, ctx);
        FLINT_ASSERT(fq_nmod_mpoly_equal(prod, A, ctx));
        fq_nmod_mpoly_clear(prod, ctx);
    }
#endif
/*
flint_printf("fq_nmod_mpoly_factor_irred_smprime_wang returning %d\n", success);
*/
	return success;
}
