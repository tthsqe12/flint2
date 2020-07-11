/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fmpz_mpoly_factor.h"
#include "nmod_mpoly_factor.h"


static void nmod_mpoly_get_mpolyu2(
    nmod_mpolyu_t A,
    const nmod_mpoly_t B,
    slong var0,
    slong var1,
    const nmod_mpoly_ctx_t ctx)
{
    slong i;
    nmod_mpoly_struct * Ac;
    ulong * Bexps;
    slong NA, NB;

    FLINT_ASSERT(ctx->minfo->ord == ORD_LEX);

    FLINT_ASSERT(var0 < ctx->minfo->nvars);
    FLINT_ASSERT(var1 < ctx->minfo->nvars);

    Bexps = (ulong *) flint_malloc(ctx->minfo->nvars*sizeof(ulong));

    NA = mpoly_words_per_exp_sp(A->bits, ctx->minfo);
    NB = mpoly_words_per_exp(B->bits, ctx->minfo);

    A->length = 0;
    for (i = 0; i < B->length; i++)
    {
        mpoly_get_monomial_ui(Bexps, B->exps + NB*i, B->bits, ctx->minfo);
        FLINT_ASSERT(FLINT_BIT_COUNT(Bexps[var0]) < FLINT_BITS/2);
        FLINT_ASSERT(FLINT_BIT_COUNT(Bexps[var1]) < FLINT_BITS/2);
        Ac = _nmod_mpolyu_get_coeff(A, pack_exp2(Bexps[var0], Bexps[var1]), ctx);
        FLINT_ASSERT(Ac->bits == A->bits);
        nmod_mpoly_fit_length(Ac, Ac->length + 1, ctx);
        Ac->coeffs[Ac->length] = B->coeffs[i];
        Bexps[var0] = 0;
        Bexps[var1] = 0;
        mpoly_set_monomial_ui(Ac->exps + NA*Ac->length, Bexps, A->bits, ctx->minfo);
        Ac->length++;
    }

    flint_free(Bexps);
/*
printf("nmod_mpoly_get_mpolyu3 returning: "); nmod_mpolyu3_print_pretty(A, ourvars[var0], ourvars[var1], ourvars[var2], ourvars, ctx); printf("\n");
*/
    FLINT_ASSERT(nmod_mpolyu_is_canonical(A, ctx));
}


/*
    return 
        -1: singular vandermonde matrix
        0:  inconsistent system
        1:  success
*/
int nmod_zip_find_coeffs_new(
    mp_limb_t * coeffs,             /* length mlength */
    const mp_limb_t * monomials,    /* length mlength */
    slong mlength,
    const mp_limb_t * evals,        /* length elength */
    slong elength,    
    nmod_poly_t master,
    const nmodf_ctx_t ffinfo)
{
    slong i, j;
    mp_limb_t V, V0, V1, V2, T, S, r, p0, p1;

    FLINT_ASSERT(elength >= mlength);

    nmod_poly_product_roots_nmod_vec(master, monomials, mlength);

    for (i = 0; i < mlength; i++)
    {
        /* coeffs[i] is (coeffs(P).values)/P(roots[i]) =: V/S
            where P(x) = master(x)/(x-roots[i])     */
        V0 = V1 = V2 = T = S = 0;
        r = monomials[i];
        for (j = mlength; j > 0; j--)
        {
            T = nmod_add(nmod_mul(r, T, ffinfo->mod), master->coeffs[j], ffinfo->mod);
            S = nmod_add(nmod_mul(r, S, ffinfo->mod), T, ffinfo->mod);
            umul_ppmm(p1, p0, evals[j - 1], T);
            add_sssaaaaaa(V2, V1, V0, V2, V1, V0, WORD(0), p1, p0);
        }
        /* roots[i] should be a root of master */
        FLINT_ASSERT(nmod_add(nmod_mul(r, T, ffinfo->mod), master->coeffs[0], ffinfo->mod) == 0);
        NMOD_RED3(V, V2, V1, V0, ffinfo->mod);
        S = nmod_mul(S, r, ffinfo->mod); /* shift is one */
        if (S == 0)
            return -1;
        coeffs[i] = nmod_mul(V, nmod_inv(S, ffinfo->mod), ffinfo->mod);
    }

    /* use the coefficients of master as temp work space */
    for (i = 0; i < mlength; i++)
    {
        master->coeffs[i] = nmod_pow_ui(monomials[i], mlength, ffinfo->mod);
    }

    /* check that the remaining points match */
    for (i = mlength; i < elength; i++)
    {
        V0 = V1 = V2 = S = 0;
        for (j = 0; j < mlength; j++)
        {
            master->coeffs[j] = nmod_mul(master->coeffs[j], monomials[j], ffinfo->mod);
            umul_ppmm(p1, p0, coeffs[j], master->coeffs[j]);
            add_sssaaaaaa(V2, V1, V0, V2, V1, V0, WORD(0), p1, p0);
        }
        NMOD_RED3(V, V2, V1, V0, ffinfo->mod);
        if (V != evals[i])
            return 0;
    }

    return 1;
}

int nmod_zip_find_coeffs_new2(
    mp_limb_t * coeffs,             /* length mlength */
    const mp_limb_t * monomials,    /* length mlength */
    slong mlength,
    const mp_limb_t * evals,        /* length elength */
    slong elength,
    const mp_limb_t * master,       /* length mlength + 1 */
    mp_limb_t * temp,               /* length mlength */
    nmod_t mod)
{
    slong i, j;
    mp_limb_t V, V0, V1, V2, T, S, r, p0, p1;

    FLINT_ASSERT(elength >= mlength);

    for (i = 0; i < mlength; i++)
    {
        /* coeffs[i] is (coeffs(P).values)/P(roots[i]) =: V/S
            where P(x) = master(x)/(x-roots[i])     */
        V0 = V1 = V2 = T = S = 0;
        r = monomials[i];
        for (j = mlength; j > 0; j--)
        {
            T = nmod_add(nmod_mul(r, T, mod), master[j], mod);
            S = nmod_add(nmod_mul(r, S, mod), T, mod);
            umul_ppmm(p1, p0, evals[j - 1], T);
            add_sssaaaaaa(V2, V1, V0, V2, V1, V0, 0, p1, p0);
        }
        /* roots[i] should be a root of master */
        FLINT_ASSERT(nmod_add(nmod_mul(r, T, mod), master[0], mod) == 0);
        NMOD_RED3(V, V2, V1, V0, mod);
        S = nmod_mul(S, r, mod); /* shift is one */
        if (S == 0)
            return -1;
        coeffs[i] = nmod_mul(V, nmod_inv(S, mod), mod);
    }

    /* check that the remaining points match */
    for (j = 0; j < mlength; j++)
        temp[j] = nmod_pow_ui(monomials[j], mlength, mod);

    for (i = mlength; i < elength; i++)
    {
        V0 = V1 = V2 = S = 0;
        for (j = 0; j < mlength; j++)
        {
            temp[j] = nmod_mul(temp[j], monomials[j], mod);
            umul_ppmm(p1, p0, coeffs[j], temp[j]);
            add_sssaaaaaa(V2, V1, V0, V2, V1, V0, 0, p1, p0);
        }
        NMOD_RED3(V, V2, V1, V0, mod);
        if (V != evals[i])
            return 0;
    }
    return 1;
}


/*
    B vars: x0 x1 x2 x3 x4 xv           2 < v
    A vars: xv x0 x1 : 0 0 x2 x3 x4 0
*/
static void nmod_mpoly_get_mpolyu3(
    nmod_mpolyu_t A,
    const nmod_mpoly_t B,
    slong var0,
    slong var1,
    slong var2,
    const nmod_mpoly_ctx_t ctx)
{
    slong i;
    nmod_mpoly_struct * Ac;
    ulong * Bexps;
    slong NA, NB;

    FLINT_ASSERT(ctx->minfo->ord == ORD_LEX);

    FLINT_ASSERT(var0 < ctx->minfo->nvars);
    FLINT_ASSERT(var1 < ctx->minfo->nvars);
    FLINT_ASSERT(var2 < ctx->minfo->nvars);

    Bexps = (ulong *) flint_malloc(ctx->minfo->nvars*sizeof(ulong));

    NA = mpoly_words_per_exp_sp(A->bits, ctx->minfo);
    NB = mpoly_words_per_exp(B->bits, ctx->minfo);

    A->length = 0;
    for (i = 0; i < B->length; i++)
    {
        mpoly_get_monomial_ui(Bexps, B->exps + NB*i, B->bits, ctx->minfo);
        FLINT_ASSERT(FLINT_BIT_COUNT(Bexps[var0]) < FLINT_BITS/3);
        FLINT_ASSERT(FLINT_BIT_COUNT(Bexps[var1]) < FLINT_BITS/3);
        FLINT_ASSERT(FLINT_BIT_COUNT(Bexps[var2]) < FLINT_BITS/3);
        Ac = _nmod_mpolyu_get_coeff(A, pack_exp3(Bexps[var0], Bexps[var1], Bexps[var2]), ctx);
        FLINT_ASSERT(Ac->bits == A->bits);
        nmod_mpoly_fit_length(Ac, Ac->length + 1, ctx);
        Ac->coeffs[Ac->length] = B->coeffs[i];
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
    FLINT_ASSERT(nmod_mpolyu_is_canonical(A, ctx));
}

void nmod_mpoly_monomial_evals(
    n_poly_t E,
    const nmod_mpoly_t A,
    const mp_limb_t * alpha,
    const nmod_mpoly_ctx_t ctx)
{
    slong i, j;
    slong offset, shift;
    slong N = mpoly_words_per_exp_sp(A->bits, ctx->minfo);
    slong nvars = ctx->minfo->nvars;
    ulong * Aexp;
    slong * LUToffset;
    ulong * LUTmask;
    mp_limb_t * LUTvalue;
    slong LUTlen;
    mp_limb_t xpoweval;
    ulong * inputexpmask;
    TMP_INIT;

    FLINT_ASSERT(A->bits <= FLINT_BITS);

    TMP_START;

    LUToffset = (slong *) TMP_ALLOC(N*FLINT_BITS*sizeof(slong));
    LUTmask   = (ulong *) TMP_ALLOC(N*FLINT_BITS*sizeof(ulong));
    LUTvalue  = (mp_limb_t *) TMP_ALLOC(N*FLINT_BITS*sizeof(mp_limb_t));

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

        xpoweval = alpha[j]; /* xpoweval = alpha[i]^(2^i) */
        for (i = 0; i < A->bits; i++)
        {
            LUToffset[LUTlen] = offset;
            LUTmask[LUTlen] = (UWORD(1) << (shift + i));
            LUTvalue[LUTlen] = xpoweval;
            if ((inputexpmask[offset] & LUTmask[LUTlen]) != 0)
            {
                LUTlen++;
            }
            xpoweval = nmod_mul(xpoweval, xpoweval, ctx->ffinfo->mod);
        }
    }
    FLINT_ASSERT(LUTlen < N*FLINT_BITS);

    n_poly_fit_length(E, A->length);
    E->length = A->length;
    for (i = 0; i < A->length; i++)
    {
        xpoweval = 1;
        for (j = 0; j < LUTlen; j++)
        {
            if (((Aexp + N*i)[LUToffset[j]] & LUTmask[j]) != 0)
            {
                xpoweval = nmod_mul(xpoweval, LUTvalue[j], ctx->ffinfo->mod);
            }
        }
        E->coeffs[i] = xpoweval;
    }

    TMP_END;
}

void _nmod_mpoly_monomial_evals(
    mp_limb_t * E,
    const ulong * Aexps,
    flint_bitcnt_t Abits,
    slong Alen,
    const mp_limb_t * alpha,
    slong vstart,
    const nmod_mpoly_ctx_t ctx)
{
    slong i, j;
    slong offset, shift;
    slong N = mpoly_words_per_exp_sp(Abits, ctx->minfo);
    slong nvars = ctx->minfo->nvars;
    slong * LUToffset;
    ulong * LUTmask;
    mp_limb_t * LUTvalue;
    slong LUTlen;
    mp_limb_t xpoweval;
    ulong * inputexpmask;
    TMP_INIT;

    TMP_START;
/*
flint_printf("_nmod_mpoly_monomial_evals called Alen: %wd\n", Alen);
*/
    LUToffset = (slong *) TMP_ALLOC(N*FLINT_BITS*sizeof(slong));
    LUTmask   = (ulong *) TMP_ALLOC(N*FLINT_BITS*sizeof(ulong));
    LUTvalue  = (mp_limb_t *) TMP_ALLOC(N*FLINT_BITS*sizeof(mp_limb_t));

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

        xpoweval = alpha[j]; /* xpoweval = alpha[i]^(2^i) */
        for (i = 0; i < Abits; i++)
        {
            LUToffset[LUTlen] = offset;
            LUTmask[LUTlen] = (UWORD(1) << (shift + i));
            LUTvalue[LUTlen] = xpoweval;
            if ((inputexpmask[offset] & LUTmask[LUTlen]) != 0)
                LUTlen++;
            xpoweval = nmod_mul(xpoweval, xpoweval, ctx->ffinfo->mod);
        }
    }
    FLINT_ASSERT(LUTlen < N*FLINT_BITS);

    for (i = 0; i < Alen; i++)
    {
        xpoweval = 1;
        for (j = 0; j < LUTlen; j++)
        {
            if (((Aexps + N*i)[LUToffset[j]] & LUTmask[j]) != 0)
            {
                xpoweval = nmod_mul(xpoweval, LUTvalue[j], ctx->ffinfo->mod);
            }
        }
        E[i] = xpoweval;
    }

    TMP_END;
}

void _nmod_mpoly_monomial_evals_indirect(
    mp_limb_t * E,
    const ulong * Aexps,
    flint_bitcnt_t Abits,
    ulong * Aind,
    slong Alen,
    const mp_limb_t * alpha,
    slong vstart,
    slong vstop,
    const nmod_mpoly_ctx_t ctx)
{
    slong i, j;
    slong offset, shift;
    slong N = mpoly_words_per_exp_sp(Abits, ctx->minfo);
    slong * LUToffset;
    ulong * LUTmask;
    mp_limb_t * LUTvalue;
    slong LUTlen;
    mp_limb_t xpoweval;
    ulong * inputexpmask;
    const ulong * thisAexp;
    TMP_INIT;

    FLINT_ASSERT(0 <= vstart);
    FLINT_ASSERT(vstart < vstop);
    FLINT_ASSERT(vstop <= ctx->minfo->nvars);

    TMP_START;

    LUToffset = (slong *) TMP_ALLOC(N*FLINT_BITS*sizeof(slong));
    LUTmask   = (ulong *) TMP_ALLOC(N*FLINT_BITS*sizeof(ulong));
    LUTvalue  = (mp_limb_t *) TMP_ALLOC(N*FLINT_BITS*sizeof(mp_limb_t));

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

        xpoweval = alpha[j]; /* xpoweval = alpha[i]^(2^i) */
        for (i = 0; i < Abits; i++)
        {
            LUToffset[LUTlen] = offset;
            LUTmask[LUTlen] = (UWORD(1) << (shift + i));
            LUTvalue[LUTlen] = xpoweval;
            if ((inputexpmask[offset] & LUTmask[LUTlen]) != 0)
                LUTlen++;
            xpoweval = nmod_mul(xpoweval, xpoweval, ctx->ffinfo->mod);
        }
    }
    FLINT_ASSERT(LUTlen < N*FLINT_BITS);

    for (i = 0; i < Alen; i++)
    {
        thisAexp = Aexps + N*Aind[i];
        xpoweval = 1;
        for (j = 0; j < LUTlen; j++)
        {
            if ((thisAexp[LUToffset[j]] & LUTmask[j]) != 0)
            {
                xpoweval = nmod_mul(xpoweval, LUTvalue[j], ctx->ffinfo->mod);
            }
        }
        E[i] = xpoweval;
    }

    TMP_END;
}




/*
    evaluation helper:
        3*i + 0 monomial evaluated at current power (changing)
        3*i + 1 coefficient (constant)
        3*i + 2 monomial evaluated at first power (constant)
*/

void nmod_mpolyu_set_eval_helper(
    n_polyun_t EH,
    const nmod_mpolyu_t A,
    const mp_limb_t * alpha,
    const nmod_mpoly_ctx_t ctx)
{
    slong i, j, n;
    n_polyun_term_struct * EHterms;
    mp_limb_t * p, * q;

    n_polyun_fit_length(EH, A->length);
    EH->length = A->length;
    EHterms = EH->terms;

    for (i = 0; i < A->length; i++)
    {
        EHterms[i].exp = A->exps[i];
        n = A->coeffs[i].length;
        n_poly_fit_length(EHterms[i].coeff, 3*n);
        nmod_mpoly_monomial_evals(EHterms[i].coeff, A->coeffs + i, alpha, ctx);
        FLINT_ASSERT(n == EHterms[i].coeff->length);
        p = EHterms[i].coeff->coeffs;
        q = A->coeffs[i].coeffs;
        for (j = n - 1; j >= 0; j--)
        {
            mp_limb_t t1 = p[j];
            mp_limb_t t2 = q[j];
            p[3*j + 0] = t1;
            p[3*j + 1] = t2;
            p[3*j + 2] = t1;
        }
    }
}


void nmod_mpoly_delete_duplicate_terms(
    nmod_mpoly_t A,
    const nmod_mpoly_ctx_t ctx)
{
    slong i, j;
    slong N = mpoly_words_per_exp(A->bits, ctx->minfo);
    j = -1;
    for (i = 0; i < A->length; i++)
    {
        if (j >= 0 && mpoly_monomial_equal(A->exps + N*j, A->exps + N*i, N))
        {
            FLINT_ASSERT(A->coeffs[j] == A->coeffs[i]);
            continue;
        }
        j++;
        A->coeffs[j] = A->coeffs[i];
        mpoly_monomial_set(A->exps + N*j, A->exps + N*i, N);
    }
    j++;
    A->length = j;
}

/*
    for each term Y^y*X^x*Z^z * pol(x1,...) in B with j < deg
    set Y^0*X^x*Z^z in H as the monomials with the monomial evals as coeffs
        merge monomial sets comming from different y's (shouldn't happen)
*/
slong nmod_mpolyu_set_eval_helper_and_zip_form(
    n_polyun_t EH,
    nmod_mpolyu_t H,
    ulong deg,
    const nmod_mpolyu_t B,
    const mp_limb_t * alpha,
    const nmod_mpoly_ctx_t ctx)
{
    slong i, j, n;
    ulong x, y, z;
    n_polyun_term_struct * EHterms;
    mp_limb_t * p, * q;
    nmod_mpoly_struct * Hc;
    slong old_len, zip_length = 0;
    slong N = mpoly_words_per_exp(B->bits, ctx->minfo);
/*
flint_printf("nmod_mpolyu_set_eval_helper_and_zip_form called\n");
flint_printf("deg = %wu\n", deg);
flint_printf("B: ");
nmod_mpolyu3_print_pretty(B, "Y", "X", "Z", NULL, ctx);
flint_printf("\n");
*/
    n_polyun_fit_length(EH, B->length);
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
        n_poly_fit_length(EHterms[i].coeff, 3*n);
        nmod_mpoly_monomial_evals(EHterms[i].coeff, B->coeffs + i, alpha, ctx);
        FLINT_ASSERT(n == EHterms[i].coeff->length);
        p = EHterms[i].coeff->coeffs;
        q = B->coeffs[i].coeffs;

        if (x < deg)
        {
            FLINT_ASSERT(y == 0 && "strange but ok");
            Hc = _nmod_mpolyu_get_coeff(H, pack_exp3(0, x, z), ctx);
            nmod_mpoly_fit_length(Hc, n, ctx);
            old_len = Hc->length;
            flint_mpn_copyi(Hc->coeffs + old_len, p, n);
            mpoly_copy_monomials(Hc->exps + N*old_len, B->coeffs[i].exps, n, N);
            Hc->length += n;
            zip_length = FLINT_MAX(zip_length, Hc->length);
            if (old_len > 0)
            {
                FLINT_ASSERT(0 && "strange but ok");
                nmod_mpoly_sort_terms(Hc, ctx);
                nmod_mpoly_delete_duplicate_terms(Hc, ctx);
            }
        }

        for (j = n - 1; j >= 0; j--)
        {
            mp_limb_t t1 = p[j];
            mp_limb_t t2 = q[j];
            p[3*j + 0] = t1;
            p[3*j + 1] = t2;
            p[3*j + 2] = t1;
        }
    }

    return zip_length;
}



slong nmod_mpoly_set_eval_helper_and_zip_form2(
    slong * deg1_, /* degree of B wrt main var 1 */
    n_polyun_t EH,
    n_polyun_t H,
    n_polyun_t M,
    const nmod_mpoly_t B,
    const mp_limb_t * betas,
    const nmod_mpoly_ctx_t ctx)
{
    slong start, Bi, j, n;
    slong e0, e1, Hi, EHi;
    n_polyun_term_struct * EHterms, * Hterms, * Mterms;
    mp_limb_t * p;
    slong zip_length = 0;
    flint_bitcnt_t Bbits = B->bits;
    const mp_limb_t * Bcoeffs = B->coeffs;
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

        n_polyun_fit_length(EH, EHi + 1);
        EHterms = EH->terms;
        EHterms[EHi].exp = pack_exp2(e0, e1);
        n_poly_fit_length(EHterms[EHi].coeff, 3*n);
        EHterms[EHi].coeff->length = n;
        p = EHterms[EHi].coeff->coeffs;
        EHi++;

        _nmod_mpoly_monomial_evals(p, Bexps + N*start, Bbits, n, betas, 2, ctx);

        if (e0 < deg0)
        {
            n_polyun_fit_length(H, Hi + 1);
            n_polyun_fit_length(M, Hi + 1);
            Hterms = H->terms;
            Mterms = M->terms;
            Hterms[Hi].exp = pack_exp2(e0, e1);
            Mterms[Hi].exp = pack_exp2(e0, e1);
            n_poly_fit_length(Hterms[Hi].coeff, n);
            zip_length = FLINT_MAX(zip_length, n);
            Hterms[Hi].coeff->length = n;
            flint_mpn_copyi(Hterms[Hi].coeff->coeffs, p, n);
            n_poly_mod_product_roots_nmod_vec(Mterms[Hi].coeff, p, n, ctx->ffinfo->mod);
            Hi++;
        }

        for (j = n - 1; j >= 0; j--)
        {
            mp_limb_t t2 = Bcoeffs[start + j];
            mp_limb_t t1 = p[j];
            p[3*j + 0] = t1;
            p[3*j + 1] = t2;
            p[3*j + 2] = t1;
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


int _fmpz_mpoly_modpk_update_zip(
    fmpz_t pk,
    fmpz_mpoly_t A,
    const fmpz_mpoly_ctx_t ctx,
    n_polyun_t Z,
    const n_polyun_t H,
    const n_polyun_t M,
    const nmod_mpoly_ctx_t ctxp)
{
    slong i, j, Ai, n;
    int success;
    slong off, shift;
    slong N = mpoly_words_per_exp_sp(A->bits, ctx->minfo);
    ulong start, mask = (-UWORD(1)) >> (FLINT_BITS - A->bits);
    n_poly_t c, t;
    mp_limb_t * ccoeffs;

    mpoly_gen_offset_shift_sp(&off, &shift, 0, A->bits, ctx->minfo);

    mask = mask << shift;

    Ai = 1;
    start = (A->exps + N*0)[off] & mask;

    while (Ai < A->length && start == ((A->exps + N*Ai)[off] & mask))
    {
        Ai++;
    }

    FLINT_ASSERT(Ai < A->length);

    if (Ai >= A->length)
        return 1;

    n_poly_init(c);
    n_poly_init(t);

    FLINT_ASSERT(Z->length == H->length);
    FLINT_ASSERT(Z->length == M->length);

    for (i = 0; i < Z->length; i++)
    {
        n = H->terms[i].coeff->length;
        FLINT_ASSERT(M->terms[i].coeff->length == n + 1);
        FLINT_ASSERT(Z->terms[i].coeff->length >= n);

        n_poly_fit_length(c, n);
        n_poly_fit_length(t, n);

        ccoeffs = c->coeffs;

        success = nmod_zip_find_coeffs_new2(c->coeffs,
                        H->terms[i].coeff->coeffs, n,
                        Z->terms[i].coeff->coeffs, Z->terms[i].coeff->length,
                        M->terms[i].coeff->coeffs, t->coeffs,
                                                            ctxp->ffinfo->mod);
        if (success <= 0)
        {
            n_poly_clear(t);
            n_poly_clear(c);
            return success;
        }

        FLINT_ASSERT(Ai + n <= A->length);

        for (j = 0; j < n; j++)
        {
            if (ctxp->ffinfo->mod.n - ccoeffs[j] < ccoeffs[j])
                fmpz_submul_ui(A->coeffs + Ai + j, pk, ctxp->ffinfo->mod.n - ccoeffs[j]);
            else
                fmpz_addmul_ui(A->coeffs + Ai + j, pk, ccoeffs[j]);
        }

        Ai += n;
    }

    FLINT_ASSERT(Ai == A->length);

    n_poly_clear(t);
    n_poly_clear(c);

    return 1;
}

/* if the coefficient doesn't exist, a new one is created (and set to zero) */
n_polyun_term_struct * n_polyun_get_term(n_polyun_t A, ulong k)
{
    slong i, j;
    n_polyun_term_struct * xk;

    for (i = 0; i < A->length && A->terms[i].exp >= k; i++)
    {
        if (A->terms[i].exp == k) 
        {
            return A->terms + i;
        }
    }

    n_polyun_fit_length(A, A->length + 1);

    for (j = A->length; j > i; j--)
        n_polyun_term_swap(A->terms + j, A->terms + j - 1);
    
    A->length++;

    xk = A->terms + i;
    xk->exp = k;
    xk->coeff->length = 0;
    return xk;
}

/*
    for each term Y^y*X^x*Z^z * pol(x1,...) in B with j < deg
    set Y^0*X^x*Z^z in H as the monomials with the monomial evals as coeffs
        merge monomial sets comming from different y's (shouldn't happen)
*/
static slong nmod_mpoly_set_eval_helper_and_zip_form(
    ulong * deg_,       /* deg_X(B), output */
    n_polyun_t EH,
    nmod_mpolyu_t H,
    const nmod_mpoly_t B,
    const mp_limb_t * alpha,
    slong yvar,         /* Y = gen(yvar) (X = gen(0), Z = gen(1))*/
    const nmod_mpoly_ctx_t ctx)
{
    slong xvar = 0;
    slong zvar = 1;
    slong i, j, n;
    ulong y, x, z;
    slong yoff, xoff, zoff;
    slong yshift, xshift, zshift;
    n_polyun_term_struct * EHterms;
    mp_limb_t * p;
    nmod_mpoly_struct * Hc;
    slong old_len, zip_length = 0;
    flint_bitcnt_t bits = B->bits;
    slong Blen = B->length;
    const ulong * Bexps = B->exps;
    const mp_limb_t * Bcoeffs = B->coeffs;
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
    FLINT_ASSERT(deg == nmod_mpoly_degree_si(B, 0, ctx));

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

    n_polyun_fit_length(EH, T->length);
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
        n_poly_fit_length(EHterms[i].coeff, 3*n);
        EHterms[i].coeff->length = n;
        p = EHterms[i].coeff->coeffs;
        ind = T->terms[i].coeff->coeffs;
        _nmod_mpoly_monomial_evals_indirect(p, Bexps, bits, ind, n, alpha,
                                                                2, yvar, ctx);
        if (x < deg)
        {
            FLINT_ASSERT(y == 0 && "strange but ok");
            Hc = _nmod_mpolyu_get_coeff(H, pack_exp3(0, x, z), ctx);
            nmod_mpoly_fit_length(Hc, n, ctx);
            old_len = Hc->length;
            flint_mpn_copyi(Hc->coeffs + old_len, p, n);
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
                nmod_mpoly_sort_terms(Hc, ctx);
                nmod_mpoly_delete_duplicate_terms(Hc, ctx);
            }
        }

        for (j = n - 1; j >= 0; j--)
        {
            mp_limb_t t1 = p[j];
            mp_limb_t t2 = Bcoeffs[ind[j]];
            p[3*j + 0] = t1;
            p[3*j + 1] = t2;
            p[3*j + 2] = t1;
        }
    }

    n_polyun_clear(T);

    *deg_ = deg;
/*
flint_printf("nmod_mpoly_set_eval_helper_and_zip_form returning\n");
*/
    return zip_length;
}



static mp_limb_t n_poly_mod_eval_step(n_poly_t A, nmod_t mod)
{
    slong i, Alen = A->length;
    mp_limb_t * Acoeffs = A->coeffs;
    ulong t0, t1, t2, p0, p1;

    FLINT_ASSERT(3*Alen <= A->alloc);

    t2 = t1 = t0 = 0;
    for (i = 0; i < Alen; i++)
    {
        umul_ppmm(p1, p0, Acoeffs[3*i + 0], Acoeffs[3*i + 1]);
        add_sssaaaaaa(t2, t1, t0, t2, t1, t0, 0, p1, p0);
        Acoeffs[3*i + 0] = nmod_mul(Acoeffs[3*i + 0], Acoeffs[3*i + 2], mod);
    }
    NMOD_RED3(t0, t2, t1, t0, mod);
    return t0;
}

void n_polyu_mod_eval_step(n_polyu_t E, n_polyun_t A, nmod_t mod)
{
    slong Ai, Ei;
    n_polyu_term_struct * Eterms;
    n_polyun_term_struct * Aterms;

    n_polyu_fit_length(E, A->length);

    Eterms = E->terms;
    Aterms = A->terms;
    Ei = 0;
    for (Ai = 0; Ai < A->length; Ai++)
    {
        FLINT_ASSERT(Ei < E->alloc);
        Eterms[Ei].exp = Aterms[Ai].exp;
        Eterms[Ei].coeff = n_poly_mod_eval_step(Aterms[Ai].coeff, mod);
        Ei += (Eterms[Ei].coeff != 0);
    }
    E->length = Ei;
}

void n_bpoly_mod_eval_step(n_bpoly_t E, n_polyun_t A, nmod_t mod)
{
    slong Ai;
    mp_limb_t c;
    ulong e0, e1;
    n_polyun_term_struct * Aterms = A->terms;

    n_bpoly_zero(E);
    for (Ai = 0; Ai < A->length; Ai++)
    {
        c = n_poly_mod_eval_step(Aterms[Ai].coeff, mod);
        e0 = extract_exp(Aterms[Ai].exp, 1, 2);
        e1 = extract_exp(Aterms[Ai].exp, 0, 2);
        if (c == 0)
            continue;
        n_bpoly_set_coeff_nonzero(E, e0, e1, c);
    }
}

void n_poly_eval_reset(n_poly_t A)
{
    slong i, Alen = A->length;
    mp_limb_t * Acoeffs = A->coeffs;

    FLINT_ASSERT(3*Alen <= A->alloc);

    for (i = 0; i < Alen; i++)
        Acoeffs[3*i + 0] = Acoeffs[3*i + 2];
}

void n_polyun_eval_reset(n_polyun_t A)
{
    slong Ai;
    for (Ai = 0; Ai < A->length; Ai++)
        n_poly_eval_reset(A->terms[Ai].coeff);
}

int n_polyu2_add_zip_must_match(
    n_polyun_t Z,
    const n_bpoly_t A,
    slong cur_length)
{
    slong i, Ai, ai;
    n_polyun_term_struct * Zt = Z->terms;
    const n_poly_struct * Acoeffs = A->coeffs;

    Ai = A->length - 1;
    ai = (Ai < 0) ? 0 : n_poly_degree(A->coeffs + Ai);

    for (i = 0; i < Z->length; i++)
    {
        if (Ai >= 0 && Zt[i].exp == pack_exp2(Ai, ai))
        {
            /* Z present, A present */
            Zt[i].coeff->coeffs[cur_length] = Acoeffs[Ai].coeffs[ai];
            Zt[i].coeff->length = cur_length + 1;
            do {
                ai--;
            } while (ai >= 0 && Acoeffs[Ai].coeffs[ai] == 0);
            if (ai < 0)
            {
                do {
                    Ai--;
                } while (Ai >= 0 && Acoeffs[Ai].length == 0);
                if (Ai >= 0)
                    ai = n_poly_degree(Acoeffs + Ai);
            }
        }
        else if (Ai < 0 || Zt[i].exp > pack_exp2(Ai, ai))
        {
            /* Z present, A missing */
            Zt[i].coeff->coeffs[cur_length] = 0;
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

void n_polyu3_add_zip_limit1(
    n_polyun_t Z,
    const n_polyun_t A,
    const ulong deg1,
    slong cur_length,
    slong fit_length)
{
    const n_polyun_term_struct * At = A->terms;
    const n_polyun_term_struct * Ait;
    n_polyun_term_struct * Zit;
    slong Ai, ai, Zi, j;

    Ai = -1;
    ai = -1;
    do {
        Ai++;
        Ait = At + Ai;
    } while (Ai < A->length && extract_exp(Ait->exp, 1, 3) >= deg1);
    if (Ai < A->length)
        ai = n_poly_degree(Ait->coeff);

    Zi = 0;

    while (Ai < A->length && Zi < Z->length)
    {
        Zit = Z->terms + Zi;
        Ait = At + Ai;
        if (Ait->exp + ai > Zit->exp)
        {
            /* missing from Z */
            n_polyun_fit_length(Z, Z->length + 1);
            for (j = Z->length; j > Zi; j--)
                n_polyun_term_swap(Z->terms + j, Z->terms + j - 1);
            Z->length++;
            Zit = Z->terms + Zi;
            Zit->exp = Ait->exp + ai;
            n_poly_fit_length(Zit->coeff, fit_length);
            Zit->coeff->length = cur_length;
            mpn_zero(Zit->coeff->coeffs, cur_length);
            goto in_both;            
        }
        else if (Ait->exp + ai < Zit->exp)
        {
            /* missing from A */
            FLINT_ASSERT(cur_length + 1 <= Zit->coeff->alloc);
            Zit->coeff->coeffs[cur_length] = 0;
            Zit->coeff->length = cur_length + 1;
            Zi++;
        }
        else
        {
in_both:
            FLINT_ASSERT(cur_length == Zit->coeff->length);
            FLINT_ASSERT(cur_length + 1 <= Zit->coeff->alloc);
            Zit->coeff->coeffs[cur_length] = Ait->coeff->coeffs[ai];
            Zit->coeff->length = cur_length + 1;
            Zi++;
            do {
                ai--;
            } while (ai >= 0 && Ait->coeff->coeffs[ai] == 0);
            if (ai < 0)
            {
                do {
                    Ai++;
                    Ait = At + Ai;
                } while (Ai < A->length && extract_exp(Ait->exp, 1, 3) >= deg1);
                if (Ai < A->length)
                    ai = n_poly_degree(Ait->coeff);
            }
        }
    }

    /* everything in A must be put on the end of Z */
    while (Ai < A->length)
    {
        Zi = Z->length;
        n_polyun_fit_length(Z, Zi + A->length - Ai);
        Zit = Z->terms + Zi;
        Zit->exp = Ait->exp + ai;
        n_poly_fit_length(Zit->coeff, fit_length);
        Zit->coeff->length = cur_length;
        mpn_zero(Zit->coeff->coeffs, cur_length);
        Z->length = ++Zi;
        FLINT_ASSERT(cur_length == Zit->coeff->length);
        FLINT_ASSERT(cur_length + 1 <= Zit->coeff->alloc);
        Zit->coeff->coeffs[cur_length] = Ait->coeff->coeffs[ai];
        Zit->coeff->length = cur_length + 1;
        do {
            ai--;
        } while (ai >= 0 && Ait->coeff->coeffs[ai] == 0);
        if (ai < 0)
        {
            do {
                Ai++;
                Ait = At + Ai;
            } while (Ai < A->length && extract_exp(Ait->exp, 1, 3) >= deg1);
            if (Ai < A->length)
                ai = n_poly_degree(Ait->coeff);
        }
    }

    /* everything in Z must have a zero appended */
    while (Zi < Z->length)
    {
        Zit = Z->terms + Zi;
        FLINT_ASSERT(cur_length == Zit->coeff->length);
        FLINT_ASSERT(cur_length + 1 <= Zit->coeff->alloc);
        Zit->coeff->coeffs[cur_length] = 0;
        Zit->coeff->length = cur_length + 1;
        Zi++;
    }

    for (Zi = 0; Zi < Z->length; Zi++)
    {
        FLINT_ASSERT(Z->terms[Zi].coeff->length == cur_length + 1);
    }
}

slong nmod_mpolyu_find_term(const nmod_mpolyu_t A, ulong e)
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
static int nmod_mpoly_from_zip(
    nmod_mpoly_t B,
    const n_polyun_t Z,
    nmod_mpolyu_t H,
    ulong deg,
    slong yvar,     /* Y = gen(yvar) */
    const nmod_mpoly_ctx_t ctx)
{
    int success;
    slong Hi, Zi, Bi, i, j;
    slong xvar = 0;
    slong zvar = 1;
    ulong x, y, z;
    flint_bitcnt_t bits = B->bits;
    mp_limb_t * Bcoeffs;
    ulong * Bexps;
    slong N = mpoly_words_per_exp_sp(bits, ctx->minfo);
    ulong mask = (-UWORD(1)) >> (FLINT_BITS - bits);
    slong xoff, xshift, yoff, yshift, zoff, zshift;
    n_polyun_term_struct * Zt = Z->terms;
    nmod_mpoly_struct * Hc;
    slong Hlen = H->length;

    n_polyun_t M;  /* temp */
    n_polyun_init(M);
/*
flint_printf("-----------------");
flint_printf("nmod_mpoly_from_zip called vars %wd, %wd, %wd\n", yvar, xvar, zvar);
flint_printf("Z: "); n_polyu3n_print_pretty(Z, "Y", "X", "Z", "_"); printf("\n");
flint_printf("H: "); nmod_mpolyu3_print_pretty(H, "Y", "X", "Z", NULL, ctx); printf("\n");
flint_printf("deg: %wd\n", deg);
*/
    FLINT_ASSERT(bits == H->bits);

    n_polyun_fit_length(M, Hlen + 1);
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
        Hi = nmod_mpolyu_find_term(H, pack_exp3(0, x, z));
        if (Hi < 0)
            return 0;

        FLINT_ASSERT(Hi < Hlen);
        FLINT_ASSERT(H->exps[Hi] == pack_exp3(0, x, z));

        Hc = H->coeffs + Hi;
        FLINT_ASSERT(bits == Hc->bits);
        FLINT_ASSERT(Hc->length > 0);
        nmod_mpoly_fit_length(B, Bi + Hc->length, ctx);
        Bcoeffs = B->coeffs;

        if (M->terms[Hi].coeff->length < 1)
        {
            n_poly_mod_product_roots_nmod_vec(M->terms[Hi].coeff,
                                     Hc->coeffs, Hc->length, ctx->ffinfo->mod);
        }

        n_poly_fit_length(M->terms[Hlen].coeff, Hc->length);

        success = nmod_zip_find_coeffs_new2(Bcoeffs + Bi, Hc->coeffs,
                    Hc->length, Zt[Zi].coeff->coeffs, Zt[Zi].coeff->length,
                    M->terms[Hi].coeff->coeffs, M->terms[Hlen].coeff->coeffs,
                                                             ctx->ffinfo->mod);
        if (success <= 0)
            return success;

        Bexps = B->exps;
        for (j = Bi, i = 0; i < Hc->length; j++, i++)
        {
            if (Bcoeffs[j] == 0)
                continue;
            Bcoeffs[Bi] = Bcoeffs[j];
            FLINT_ASSERT(Bi < B->alloc);
            mpoly_monomial_set(Bexps + N*Bi, Hc->exps + N*i, N);
            (Bexps + N*Bi)[yoff] += y << yshift;
            Bi++;
        }
    }
    B->length = Bi;
    nmod_mpoly_sort_terms(B, ctx);
    FLINT_ASSERT(nmod_mpoly_is_canonical(B, ctx));
/*
flint_printf("nmod_mpoly_from_zip returning good\n");
flint_printf("B: "); nmod_mpoly_print_pretty(B, NULL, ctx); flint_printf("\n");
*/
    n_polyun_clear(M);

    return 1;
}


/*
    ctx is for Fp[x_0, ..., x_n] where n = ctx->minfo->nvars - 1

    We are lifting a factorization mod
        (x_m - alpha[m-1], ..., x_n - alpha[n-1])
    to a factorization mod
        (x_(m+1) - alpha[m], ..., x-n - alpha[n-1]),
    We are lifting variable x_m, and all polynomials live in Fp[x_0, ..., x_m].

    It is required that m >= 3. Rewrite polynomials as nmod_mpolyu's to
    live in Fp[x_2, ..., x_(m-1)][x_m, x_0, x_1] still using ctx for
    the Fp[x_2, ..., x_(m-1)] part. Evaluate inputs at
        (x_2, ..., x_(m-1)) = (beta_2^i, ..., beta_(m-1)^i) for i = 1, 2, ..., s
    The evaluations are trivariate n_polyu's in Fp[x_m, x_0, x_1], aka n_polyu3.
    Lift the evaluations in Fp[x_m, x_0, x_1] mod (x_m - alpha[m-1]) using
    n_polyu3_mod_factor_lift. The lifted factors come out as n_polyun's in
    Fp[x_1][x_m, x_0], i.e. x_1 is in dense storage. Interpolate the factors in
    Fp[x_2, ..., x_(m-1)][x_m, x_0, x_1] using Zippel's assumption.

    return:
       -1: function failed with no conclusion
        0: lifting is impossible without changing the lc_(x_0)
        1: B was successfully updated with the lifted factors
*/
int nmod_mpoly_hlift_zippel(
    slong r,
    nmod_mpoly_struct * B,
    flint_rand_t state,
    const mp_limb_t * alpha,
    const nmod_mpoly_t A,
    const slong * degs,
    slong m,
    const nmod_mpoly_ctx_t ctx)
{
    int success;
    slong i;
    slong zip_fails_remaining;
    slong req_zip_images, cur_zip_image;
    nmod_mpolyu_struct Au[1], * H;
    n_polyun_struct Aeh[1], * Beh;
    n_polyu_struct Aeval[1], * Beval;
    n_polyun_struct * BBeval, * Z;
    mp_limb_t * beta;
    flint_bitcnt_t bits = A->bits;
    nmod_mpoly_t T1, T2;
    ulong * Bdegs;
    const slong degs0 = degs[0];
/*
    const char * vars[] = {"x", "y", "z", "w", "u", "v"};
flint_printf("nmod_mpoly_factor_lift_tuncer called (m = %wd, bits = %wu)\n", m, bits); fflush(stdout);
for (i = 0; i < n; i++)
{
flint_printf("alpha[%wd]: %wu\n", i, alpha[i]);
}
flint_printf("A: "); nmod_mpoly_print_pretty(A, vars, ctx); flint_printf("\n");
for (i = 0; i < r; i++)
{
flint_printf("B[%wd]: ", i); nmod_mpoly_print_pretty(B + i, vars, ctx); flint_printf("\n");
}
*/
    FLINT_ASSERT(m > 2);
    FLINT_ASSERT(r > 1);
    FLINT_ASSERT(bits <= FLINT_BITS);

    beta = (mp_limb_t *) flint_malloc(ctx->minfo->nvars*sizeof(mp_limb_t));
    Bdegs = (ulong *) flint_malloc(r*sizeof(ulong));
    H = (nmod_mpolyu_struct *) flint_malloc(r*sizeof(nmod_mpolyu_struct));
    Beh = (n_polyun_struct *) flint_malloc(r*sizeof(n_polyun_struct));
    Beval = (n_polyu_struct *) flint_malloc(r*sizeof(n_polyu_struct));
    BBeval = (n_polyun_struct *) flint_malloc(r*sizeof(n_polyun_struct));
    Z = (n_polyun_struct *) flint_malloc(r*sizeof(n_polyun_struct));

    nmod_mpolyu_init(Au, bits, ctx);
    n_polyun_init(Aeh);
    n_polyu_init(Aeval);
    for (i = 0; i < r; i++)
    {
        nmod_mpolyu_init(H + i, bits, ctx);
        n_polyun_init(Beh + i);
        n_polyu_init(Beval + i);
        n_polyun_init(BBeval + i);
        n_polyun_init(Z + i);
    }

    /* init done */

    for (i = 0; i < r; i++)
    {
        success = nmod_mpoly_repack_bits_inplace(B + i, bits, ctx);
        if (!success)
            goto cleanup;
    }

    zip_fails_remaining = 3;

    nmod_mpoly_get_mpolyu3(Au, A, m, 0, 1, ctx);

choose_betas:

    /* only beta[2], beta[3], ..., beta[m - 1] will be used */
    FLINT_ASSERT(ctx->ffinfo->mod.n > 3);
    for (i = 0; i < ctx->minfo->nvars; i++)
        beta[i] = n_urandint(state, ctx->ffinfo->mod.n - 3) + 2;

    nmod_mpolyu_set_eval_helper(Aeh, Au, beta, ctx);

    req_zip_images = 1;
    for (i = 0; i < r; i++)
    {
        slong this_zip_images;
        this_zip_images = nmod_mpoly_set_eval_helper_and_zip_form(Bdegs + i,
                                          Beh + i, H + i, B + i, beta, m, ctx);
        req_zip_images = FLINT_MAX(req_zip_images, this_zip_images);
        FLINT_ASSERT(Bdegs[i] > 0);
    }

    cur_zip_image = 0;

next_zip_image:

    n_polyu_mod_eval_step(Aeval, Aeh, ctx->ffinfo->mod);
    for (i = 0; i < r; i++)
    {
        n_polyu_mod_eval_step(Beval + i, Beh + i, ctx->ffinfo->mod);
    }

    success = n_polyu3_mod_hensel_lift(r, BBeval, Aeval, Beval,
                                             alpha[m - 1], degs0, ctx->ffinfo);
    if (success < 1)
    {
        FLINT_ASSERT(0 && "spurious failure");
        if (--zip_fails_remaining >= 0)
            goto choose_betas;
        success = 0;
        goto cleanup;
    }

    for (i = 0; i < r; i++)
    {
        n_polyu3_add_zip_limit1(Z + i, BBeval + i, Bdegs[i],
                                                cur_zip_image, req_zip_images); 
    }

    cur_zip_image++;
    if (cur_zip_image < req_zip_images)
        goto next_zip_image;

    for (i = 0; i < r; i++)
    {
        success = nmod_mpoly_from_zip(B + i, Z + i, H + i, Bdegs[i], m, ctx);
        if (success < 1)
        {
            FLINT_ASSERT(0 && "spurious failure");
            success = 0;
            goto cleanup;
        }
    }

    nmod_mpoly_init3(T1, A->length, bits, ctx);
    nmod_mpoly_init3(T2, A->length, bits, ctx);
    nmod_mpoly_mul(T1, B + 0, B + 1, ctx);
    for (i = 2; i < r; i++)
    {
        nmod_mpoly_mul(T2, T1, B + i, ctx);
        nmod_mpoly_swap(T1, T2, ctx);
    }

    success = nmod_mpoly_equal(T1, A, ctx);
    nmod_mpoly_clear(T1, ctx);
    nmod_mpoly_clear(T2, ctx);

cleanup:

    nmod_mpolyu_clear(Au, ctx);
    n_polyun_clear(Aeh);
    n_polyu_clear(Aeval);
    for (i = 0; i < r; i++)
    {
        nmod_mpolyu_clear(H + i, ctx);
        n_polyun_clear(Beh + i);
        n_polyu_clear(Beval + i);
        n_polyun_clear(BBeval + i);
        n_polyun_clear(Z + i);
    }

    flint_free(beta);
    flint_free(Bdegs);
    flint_free(H);
    flint_free(Beh);
    flint_free(Beval);
    flint_free(BBeval);
    flint_free(Z);

    return success;
}


static void _nmod_mpoly_set_fmpz_mpoly(
    nmod_mpoly_t Ap,
    const nmod_mpoly_ctx_t ctxp,
    const fmpz_mpoly_t A,
    const fmpz_mpoly_ctx_t ctx)
{
    slong N = mpoly_words_per_exp(A->bits, ctx->minfo);
    slong Ap_len, i;
    FLINT_ASSERT(ctx->minfo->nvars == ctx->minfo->nvars);
    FLINT_ASSERT(ctx->minfo->ord == ctx->minfo->ord);
    nmod_mpoly_fit_bits(Ap, A->bits, ctxp);
    Ap->bits = A->bits;
    nmod_mpoly_fit_length(Ap, A->length, ctxp);
    Ap_len = 0;
    for (i = 0; i < A->length; i++)
    {
        Ap->coeffs[Ap_len] = fmpz_fdiv_ui(A->coeffs + i, ctxp->ffinfo->mod.n);
        if (Ap->coeffs[Ap_len] == 0)
            continue;
        mpoly_monomial_set(Ap->exps + N*Ap_len, A->exps + N*i, N);
        Ap_len++;
    }
    Ap->length = Ap_len;
}

static void _fmpz_mpoly_modpk_taylor_coeff(
    const fmpz_t pk,
    nmod_mpoly_t T,
    const nmod_mpoly_ctx_t ctxp,
    const fmpz_mpoly_t E,
    const fmpz_mpoly_ctx_t ctx)
{
    slong i, Tlen;
    slong N = mpoly_words_per_exp(E->bits, ctx->minfo);
    fmpz_t t;

    fmpz_init(t);

    nmod_mpoly_fit_bits(T, E->bits, ctxp);
    T->bits = E->bits;
    nmod_mpoly_fit_length(T, E->length, ctxp);
    Tlen = 0;
    for (i = 0; i < E->length; i++)
    {
        FLINT_ASSERT(fmpz_divisible(E->coeffs + i, pk)); /* TODO !!! */
        fmpz_divexact(t, E->coeffs + i, pk);
        T->coeffs[Tlen] = fmpz_fdiv_ui(t, ctxp->ffinfo->mod.n);
        if (T->coeffs[Tlen] == 0)
            continue;
        mpoly_monomial_set(T->exps + N*Tlen, E->exps + N*i, N);
        Tlen++;
    }
    T->length = Tlen;

    fmpz_clear(t);
}

void fmpz_set_nmods(fmpz_t a, ulong n, ulong p)
{
    if (p - n < n)
        fmpz_neg_ui(a, p - n);
    else
        fmpz_set_ui(a, n);
}

/*
    Operation on each coeff:
        A = (A < 0 ? A + pk : A) + pk*mods(Ap, p)
*/
void _fmpz_mpoly_modpk_update(
    const fmpz_t pk,
    fmpz_mpoly_t A,
    const fmpz_mpoly_ctx_t ctx,
    const nmod_mpoly_t Ap,
    const nmod_mpoly_ctx_t ctxp)
{
    slong i;
    fmpz_mpoly_t T;
    slong N = mpoly_words_per_exp(Ap->bits, ctx->minfo);
/*
    for (i = 0; i < A->length; i++)
    {
        if (fmpz_sgn(A->coeffs + i) < 0)
            fmpz_add(A->coeffs + i, A->coeffs + i, pk);
    }

printf("positive A: "); fmpz_mpoly_print_pretty(A, NULL, ctx); printf("\n");
*/

    fmpz_mpoly_init3(T, Ap->length, Ap->bits, ctx);
    T->length = Ap->length;
    for (i = 0; i < Ap->length; i++)
    {
        fmpz_set_nmods(T->coeffs + i, Ap->coeffs[i], ctxp->ffinfo->mod.n);
        fmpz_mul(T->coeffs + i, T->coeffs + i, pk);
        mpoly_monomial_set(T->exps + N*i, Ap->exps + N*i, N);
    }

    fmpz_mpoly_add(A, A, T, ctx);
    fmpz_mpoly_clear(T, ctx);
}


void _fmpz_mpoly_set_nmod_mpoly_smod(
    fmpz_mpoly_t A,
    const fmpz_mpoly_ctx_t ctx,
    const nmod_mpoly_t Ap,
    const nmod_mpoly_ctx_t ctxp)
{
    slong N = mpoly_words_per_exp(Ap->bits, ctxp->minfo);
/*
flint_printf("_fmpz_mpoly_set_nmod_mpoly_smod called bits = %wu\n", Ap->bits);
*/
    FLINT_ASSERT(ctx->minfo->ord == ctxp->minfo->ord);
    FLINT_ASSERT(ctx->minfo->nvars == ctxp->minfo->nvars);

    A->length = 0;
    fmpz_mpoly_fit_bits(A, Ap->bits, ctx);
    A->bits = Ap->bits;

    fmpz_mpoly_fit_length(A, Ap->length, ctx);
    A->length = Ap->length;

    mpoly_copy_monomials(A->exps, Ap->exps, Ap->length, N);
    _fmpz_vec_set_nmod_vec(A->coeffs, Ap->coeffs, Ap->length,
                                                      ctxp->ffinfo->mod);
/*
flint_printf("_fmpz_mpoly_set_nmod_mpoly_smod returning\n");
*/
}



/*
    The B[i] have the correct lc_x1 over ZZ.
    A = prod_i Bp[i] mod p
*/
int fmpz_mfactor_lift_prime_power(
    slong r,
    fmpz_mpoly_struct * B,
    const nmod_mpoly_struct * Bp,
    const fmpz_mpoly_t A,
    const slong * degs,
    const mp_limb_t * alphap,
    const fmpz_mpoly_ctx_t ctx,
    const nmod_mpoly_ctx_t ctxp,
    slong L)
{
    int success;
    slong i, k;
    nmod_mpoly_pfrac_t I;
    slong n = ctx->minfo->nvars - 1;
    fmpz_t pk;
    nmod_mpoly_struct * deltas;
    fmpz_mpoly_t e, t1, t2;
    nmod_mpoly_t tk;
    const char * vars [] = {"x", "y", "z", "w", "t", "u" ,"v"};

flint_printf("fmpz_mfactor_lift_prime_power\n");

    for (i = 0; i < r; i++)
    {
flint_printf("initial B[%wd]: ", i); fmpz_mpoly_print_pretty(B + i, vars, ctx); printf("\n");
    }


    FLINT_ASSERT(r > 1);

    nmod_mpoly_pfrac_init(I, r, n, Bp, alphap, ctxp);
    deltas = I->deltas + n*I->l;

    fmpz_init(pk);
    fmpz_one(pk);

    fmpz_mpoly_init(e, ctx);
    fmpz_mpoly_init(t1, ctx);
    fmpz_mpoly_init(t2, ctx);

    nmod_mpoly_init(tk, ctxp);

    k = 1;
    while (1)
    {
        fmpz_mul_ui(pk, pk, ctxp->ffinfo->mod.n);
        fmpz_mpoly_mul(t1, B + 0, B + 1, ctx);
        for (i = 2; i < r; i++)
        {
            fmpz_mpoly_mul(t2, t1, B + i, ctx);
            fmpz_mpoly_swap(t1, t2, ctx);
        }
        fmpz_mpoly_sub(e, A, t1, ctx);

flint_printf("e: "); fmpz_mpoly_print_pretty(e, vars, ctx); flint_printf("\n");

        if (fmpz_mpoly_is_zero(e, ctx))
        {
            success = 1;
            goto cleanup;
        }

        if (k > L)
        {
            success = 0;
            goto cleanup;
        }

        _fmpz_mpoly_modpk_taylor_coeff(pk, tk, ctxp, e, ctx);

flint_printf("calling nmod disolve\n");
        success = nmod_mpoly_pfrac(A->bits, n, r, tk, alphap, degs, I, ctxp);
        FLINT_ASSERT(success);

        for (i = 0; i < r; i++)
        {
flint_printf("delta[%wd]: ", i); nmod_mpoly_print_pretty(deltas + i, vars, ctxp); flint_printf("\n");
            _fmpz_mpoly_modpk_update(pk, B + i, ctx, deltas + i, ctxp);
flint_printf("updated B[%wd]: ", i); fmpz_mpoly_print_pretty(B + i, vars, ctx); flint_printf("\n");

        }
    }

    FLINT_ASSERT(0 && "unreachable");

cleanup:

    nmod_mpoly_pfrac_clear(I, ctxp);

    return success;
}


/*
    n_polyun_t Beh has all x0^i*x1^j*poly(x2, ...) with
        coeffs triples suitable for sequential eval at x2,... = beta^i

    n_polyun_t H has x0^i*x1^j*poly(x2, ...) with i < deg_x0(B[i]) with
        the monomial evals at beta in the coeffs

    n_polyun_t M has x0^i*x1^j*poly(x2, ...) with i < deg_x0(B[i])
        the master poly of H[x0^i*x1^j] in the coeff

    n_polyun_t Z has x0^i*x1^j with i < deg_x0(B[i])
        ready to collect images
*/
int fmpz_mfactor_lift_prime_power_tuncer(
    slong r,
    fmpz_mpoly_struct * B,
    flint_rand_t state,
    const nmod_mpoly_struct * Bp,
    const fmpz_mpoly_t A,
    const mp_limb_t * alphap,
    const fmpz_mpoly_ctx_t ctx,
    const nmod_mpoly_ctx_t ctxp,
    slong L)
{
    slong req_zip_images, cur_zip_image;
    flint_bitcnt_t bits = A->bits;
    slong n = ctxp->minfo->nvars;
    int success;
    slong i, j, k;
    n_polyun_struct * H, * M, * Z, * Beh;
    n_bpoly_struct * Ceval, * Beval, Teval[1];
    slong * Cdegs1;
    nmod_mpoly_t Tp;
    nmod_mpolyu_t Tu;
    n_polyun_t Teh;
    fmpz_t pk;
    mp_limb_t * betas;
    fmpz_mpoly_t e, t1, t2;
/*
    const char * vars [] = {"x", "y", "z", "w", "t", "u" ,"v"};

flint_printf("fmpz_mfactor_lift_prime_power_tuncer called\n");
    for (i = 0; i < r; i++)
    {
flint_printf("initial B[%wd]: ", i); fmpz_mpoly_print_pretty(B + i, vars, ctx); flint_printf("\n");

    }
*/
    FLINT_ASSERT(r > 1);

    fmpz_init(pk);

    fmpz_mpoly_init(e, ctx);
    fmpz_mpoly_init(t1, ctx);
    fmpz_mpoly_init(t2, ctx);

    nmod_mpoly_init(Tp, ctxp);
    nmod_mpolyu_init(Tu, bits, ctxp);
    n_polyun_init(Teh);

    Cdegs1 = (slong *) flint_malloc(r*sizeof(slong));

    betas = (mp_limb_t *) flint_malloc(n*sizeof(mp_limb_t));

    FLINT_ASSERT(ctxp->ffinfo->mod.n > 3);
    for (i = 0; i < ctx->minfo->nvars; i++)
        betas[i] = n_urandint(state, ctxp->ffinfo->mod.n - 3) + 2;

    Beh = (n_polyun_struct *) flint_malloc(r*sizeof(n_polyun_struct));
    H = (n_polyun_struct *) flint_malloc(r*sizeof(n_polyun_struct));
    M = (n_polyun_struct *) flint_malloc(r*sizeof(n_polyun_struct));
    Z = (n_polyun_struct *) flint_malloc(r*sizeof(n_polyun_struct));
    Ceval = (n_bpoly_struct *) flint_malloc(r*sizeof(n_bpoly_struct));
    Beval = (n_bpoly_struct *) flint_malloc(r*sizeof(n_bpoly_struct));
    for (i = 0; i < r; i++)
    {
        n_polyun_init(Beh + i);
        n_polyun_init(H + i);
        n_polyun_init(M + i);
        n_polyun_init(Z + i);
        n_bpoly_init(Ceval + i);
        n_bpoly_init(Beval + i);
    }

    n_bpoly_init(Teval);

    /* init done */

    req_zip_images = 1;
    for (i = 0; i < r; i++)
    {
        slong this_images;
        this_images = nmod_mpoly_set_eval_helper_and_zip_form2(Cdegs1 + i,
                                   Beh + i, H + i, M + i, Bp + i, betas, ctxp);
        req_zip_images = FLINT_MAX(req_zip_images, this_images);
    }

    for (i = 0; i < r; i++)
    {
        n_polyun_fit_length(Z + i, H[i].length);
        Z[i].length = H[i].length;
        for (j = 0; j < H[i].length; j++)
        {
            Z[i].terms[j].exp = H[i].terms[j].exp;
            n_poly_fit_length(Z[i].terms[j].coeff, req_zip_images);
            Z[i].terms[j].coeff->length = 0;
        }
    }

    fmpz_one(pk);

    k = 1;

next_power:

    fmpz_mul_ui(pk, pk, ctxp->ffinfo->mod.n);

    fmpz_mpoly_mul(t1, B + 0, B + 1, ctx);
    for (i = 2; i < r; i++)
    {
        fmpz_mpoly_mul(t2, t1, B + i, ctx);
        fmpz_mpoly_swap(t1, t2, ctx);
    }
    fmpz_mpoly_sub(e, A, t1, ctx);

    if (fmpz_mpoly_is_zero(e, ctx))
    {
        success = 1;
        goto cleanup;
    }

    if (k > L)
    {
FLINT_ASSERT(0 && "spurious failure");
        success = 0;
        goto cleanup;
    }

    _fmpz_mpoly_modpk_taylor_coeff(pk, Tp, ctxp, e, ctx);

    nmod_mpoly_get_mpolyu2(Tu, Tp, 0, 1, ctxp);

    nmod_mpolyu_set_eval_helper(Teh, Tu, betas, ctxp);
    if (fmpz_cmp_ui(pk, ctxp->ffinfo->mod.n) > 0)
    {
        for (i = 0; i < r; i++)
            n_polyun_eval_reset(Beh + i);
    }

    cur_zip_image = 0;

next_zip_image:

    n_bpoly_mod_eval_step(Teval, Teh, ctxp->ffinfo->mod);
    for (i = 0; i < r; i++)
        n_bpoly_mod_eval_step(Beval + i, Beh + i, ctxp->ffinfo->mod);

    success = n_bpoly_mod_disolve(r, Ceval, Cdegs1, Teval, Beval,
                                                       ctxp->ffinfo->mod);
    if (success <= 0)
    {
        FLINT_ASSERT(0 && "spurious failure");
        success = 0;
        goto cleanup;
    }

    for (i = 0; i < r; i++)
    {
        success = n_polyu2_add_zip_must_match(Z + i, Ceval + i, cur_zip_image);
        if (!success)
        {
            FLINT_ASSERT(0 && "spurious failure");
            success = 0;
            goto cleanup;
        }
    }

    cur_zip_image++;
    if (cur_zip_image < req_zip_images)
        goto next_zip_image;

    for (i = 0; i < r; i++)
    {
        success = _fmpz_mpoly_modpk_update_zip(pk, B + i, ctx,
                                                    Z + i, H + i, M + i, ctxp);
        if (success <= 0)
        {
            FLINT_ASSERT(0 && "spurious failure");
            success = 0;
            goto cleanup;
        }
    }

    goto next_power;

cleanup:

    fmpz_clear(pk);

    fmpz_mpoly_clear(e, ctx);
    fmpz_mpoly_clear(t1, ctx);
    fmpz_mpoly_clear(t2, ctx);

    nmod_mpoly_clear(Tp, ctxp);
    nmod_mpolyu_clear(Tu, ctxp);
    n_polyun_clear(Teh);

    flint_free(Cdegs1);

    flint_free(betas);

    for (i = 0; i < r; i++)
    {
        n_polyun_clear(Beh + i);
        n_polyun_clear(H + i);
        n_polyun_clear(M + i);
        n_polyun_clear(Z + i);
        n_bpoly_clear(Ceval + i);
        n_bpoly_clear(Beval + i);
    }
    flint_free(Beh);
    flint_free(H);
    flint_free(M);
    flint_free(Z);
    flint_free(Ceval);
    flint_free(Beval);

    n_bpoly_clear(Teval);

    FLINT_ASSERT(success == 0 || success == 1);
    return success;
}


void nmod_poly_set_fmpz_poly(nmod_poly_t a, const fmpz_poly_t b)
{
    slong i;
    nmod_poly_fit_length(a, b->length);
    for (i = 0; i < b->length; i++)
        a->coeffs[i] = fmpz_fdiv_ui(b->coeffs + i, a->mod.n);
    a->length = b->length;
    _nmod_poly_normalise(a);
}

int fmpz_mpoly_factor_irred_tuncer(
    fmpz_mpoly_factor_t fac,
    const fmpz_mpoly_t A,
    const fmpz_mpoly_factor_t lcAfac,
    const fmpz_mpoly_t lcA,
    const fmpz_mpoly_ctx_t ctx)
{
    int success;
    const slong n = ctx->minfo->nvars - 1;
    slong i, j, k;
    fmpz * alpha;
    fmpz_mpoly_struct * Aevals;
    slong * degs, * tdegs;
    fmpz_mpoly_factor_t tfac;
    fmpz_mpoly_t t, Acopy;
    fmpz_mpoly_struct * newA;
    fmpz_poly_t Au;
    fmpz_poly_factor_t Aufac;
    slong alpha_bits, alpha_count;
    flint_rand_t state;
    fmpz_mpoly_t m, mpow;
    fmpz_mpolyv_t Alc, lc_divs;
    fmpz_t q, facBound;
    mp_limb_t p;
    nmod_mpoly_ctx_t ctxp;
    nmod_mpoly_factor_t facp, tfacp;
    nmod_mpolyv_t Aevalp, Alcp;
    nmod_poly_t Aup;
    mp_limb_t * alphap;
    slong r, L;

    FLINT_ASSERT(n > 1);
    FLINT_ASSERT(ctx->minfo->ord == ORD_LEX);
    FLINT_ASSERT(A->bits <= FLINT_BITS);
    FLINT_ASSERT(fmpz_mpoly_factor_matches(lcA, lcAfac, ctx));

/*
flint_printf("_try_tuncer called n = %wd\n", n);
flint_printf("A: "); fmpz_mpoly_print_pretty(A, NULL, ctx); printf("\n");
*/
    flint_randinit(state);

    fmpz_init(facBound);
    fmpz_init(q);

    fmpz_mpoly_init(Acopy, ctx);
    fmpz_mpoly_init(m, ctx);
    fmpz_mpoly_init(mpow, ctx);

    fmpz_mpolyv_init(Alc, ctx);
    fmpz_mpolyv_init(lc_divs, ctx);

    fmpz_poly_factor_init(Aufac);
    fmpz_poly_init(Au);

    alpha = _fmpz_vec_init(n);
    alphap = (mp_limb_t *) flint_malloc(n*sizeof(mp_limb_t));

    degs  = (slong *) flint_malloc(2*(n + 1)*sizeof(slong));
    tdegs = degs + (n + 1);

    Aevals = (fmpz_mpoly_struct *) flint_malloc(n*sizeof(fmpz_mpoly_struct));
    for (i = 0; i < n; i++)
        fmpz_mpoly_init(Aevals + i, ctx);

    fmpz_mpoly_factor_init(tfac, ctx);
    fmpz_mpoly_init(t, ctx);

    nmod_mpoly_ctx_init(ctxp, n + 1, ORD_LEX, 2); /* modulus no care */
    nmod_mpoly_factor_init(facp, ctxp);
    nmod_mpoly_factor_init(tfacp, ctxp);
    nmod_mpolyv_init(Aevalp, ctxp);
    nmod_mpolyv_init(Alcp, ctxp);
    nmod_poly_init_mod(Aup, ctxp->ffinfo->mod);

    /* init done */

    fmpz_mpoly_degrees_si(degs, A, ctx);

    alpha_count = 0;
    alpha_bits = 10;

next_alpha:

    /*
        choose the αi from +-[1, 2^alpha_bits].
        alpha_bits is incremented every so often until it gets too high.
    */

    alpha_count++;
    if (alpha_count >= alpha_bits)
    {
        alpha_count = 0;
        alpha_bits++;
        if (alpha_bits > FLINT_BITS/2)
        {
            success = 0;
            goto cleanup;
        }
    }

    for (i = 0; i < n; i++)
    {
        ulong l = n_randlimb(state);
        ulong mask = UWORD(1) << alpha_bits;
        if (l & mask)
            fmpz_neg_ui(alpha + i, 1 + (l & (mask - 1)));
        else
            fmpz_set_ui(alpha + i, 1 + (l & (mask - 1)));
    }
/*
usleep(1000);
printf("---------------------------\n");
printf("alpha = "); tuple_print(alpha, n);
*/
    /* ensure degrees do not drop under evaluation */
    for (i = n - 1; i >= 0; i--)
    {
        fmpz_mpoly_evaluate_one_fmpz(Aevals + i, i == n - 1 ? A : Aevals + i + 1,
                                                        i + 1, alpha + i, ctx);
        fmpz_mpoly_degrees_si(tdegs, Aevals + i, ctx);
        for (j = 0; j <= i; j++)
        {
            if (tdegs[j] != degs[j])
                goto next_alpha;
        }
    }

    /* make sure our univar is squarefree */
    success = fmpz_mpoly_get_fmpz_poly(Au, Aevals + 0, 0, ctx);
    FLINT_ASSERT(success);
    fmpz_poly_factor(Aufac, Au);
    r = Aufac->num;
    FLINT_ASSERT(r >= 1);
    for (j = 0; j < r; j++)
    {
        if (Aufac->exp[j] != 1)
            goto next_alpha;
    }

    /* if the univariate is irreducible, then A is irreducible */
    if (r < 2)
    {
        fmpz_mpoly_factor_one(fac, ctx);
        fmpz_mpoly_factor_append_ui(fac, A, 1, ctx);
        success = 1;
        goto cleanup;
    }

    if (lcAfac->num > 0)
    {
        success = fmpz_mpoly_factor_lcc_wang(lc_divs, lcAfac, Aufac, alpha, ctx);
        if (!success)
            goto next_alpha;
    }
    else
    {
        /* lcA is constant */
        fmpz_mpolyv_fit_length(lc_divs, r, ctx);
        lc_divs->length = r;
        for (i = 0; i < r; i++)
        {
            FLINT_ASSERT(Aufac->p[i].length > 0);
            fmpz_mpoly_set_fmpz(lc_divs->coeffs + i,
                          Aufac->p[i].coeffs + Aufac->p[i].length - 1, ctx);
        }
    }

    /*
        Assuming no extraneous factors, we have

            A(X, x1, ..., xn) = F1 * ... * Fr  where  r = Aufac->num

        and lead_divisor[i] is a divisor of lc_X(Fi). We also have the
        univariate factorization

            A(X, α1, ..., αn) = (c1 X^? + ... ) * ... * (cr X^? + ... )

        Set c(x1, ..., xn) = lc_X(A) and

            m(x1, ..., xn) = c/(prod_i lc_divs[i])   division is exact

        Lift the univariate factorization

            m(α)^(r-1) f(X, α) = (m(α)*lc_divs[0](α) X^? + ...) * ... *
                                 (m(α)*lc_divs[r-1](α) X^? + ...)

        against

            m(x1, ..., xn)^(r-1) f(X, x1, ..., xn)

        Note m(x1, ..., xn) is usually constant here, but it does not have to be.
    */

    FLINT_ASSERT(r > 1);
    success = fmpz_mpoly_divides(m, lcA, lc_divs->coeffs + 0, ctx);
    FLINT_ASSERT(success);
    for (i = 1; i < r; i++)
    {
        success = fmpz_mpoly_divides(m, m, lc_divs->coeffs + i, ctx);
        FLINT_ASSERT(success);
    }

    fmpz_mpoly_pow_ui(mpow, m, r - 1, ctx);
    if (fmpz_mpoly_is_one(mpow, ctx))
    {
        newA = (fmpz_mpoly_struct *) A;
    }
    else
    {
        newA = Acopy;
        fmpz_mpoly_mul(newA, A, mpow, ctx);
    }

    if (newA->bits > FLINT_BITS)
    {
        success = 0;
        goto cleanup;
    }

    fmpz_mpoly_degrees_si(degs, newA, ctx);
    _fmpz_vec_height(facBound, newA->coeffs, newA->length);
    if (!fmpz_mpoly_factor_bound_si(facBound, facBound, degs, n + 1))
    {
        success = 0;
        goto cleanup;
    }
    fmpz_mul_2exp(facBound, facBound, 1);
    fmpz_add_ui(facBound, facBound, 1);

/*
printf("mpow: "); fmpz_mpoly_print_pretty(mpow, NULL, ctx); printf("\n");

flint_printf("modified A: "); fmpz_mpoly_print_pretty(newA, NULL, ctx); printf("\n");
*/
    fmpz_mpoly_set(t, mpow, ctx);
    for (i = n - 1; i >= 0; i--)
    {
        fmpz_mpoly_evaluate_one_fmpz(t, mpow, i + 1, alpha + i, ctx);
        fmpz_mpoly_swap(t, mpow, ctx);
        fmpz_mpoly_mul(Aevals + i, Aevals + i, mpow, ctx);
    }

    fmpz_mpolyv_fit_length(Alc, (n + 1)*r, ctx);

    i = n;
    for (j = 0; j < r; j++)
    {
        fmpz_mpoly_mul(Alc->coeffs + i*r + j, lc_divs->coeffs + j, m, ctx);
    }
    for (i = n - 1; i >= 0; i--)
    {
        for (j = 0; j < r; j++)
        {
            fmpz_mpoly_evaluate_one_fmpz(Alc->coeffs + i*r + j,
                 Alc->coeffs + (i + 1)*r + j, i + 1, alpha + i, ctx);
        }
    }

    fmpz_mpoly_factor_fit_length(fac, r, ctx);
    fac->num = r;
    fmpz_one(fac->constant);
    for (i = 0; i < r; i++)
    {
        FLINT_ASSERT(fmpz_mpoly_is_fmpz(Alc->coeffs + 0*r + i, ctx));
        FLINT_ASSERT(fmpz_mpoly_length(Alc->coeffs + 0*r + i, ctx) == 1);
        FLINT_ASSERT(fmpz_divisible(Alc->coeffs[i].coeffs + 0,
                                 Aufac->p[i].coeffs + Aufac->p[i].length - 1));
        fmpz_divexact(q, Alc->coeffs[i].coeffs + 0,
                                  Aufac->p[i].coeffs + Aufac->p[i].length - 1);
        _fmpz_mpoly_set_fmpz_poly(fac->poly + i, newA->bits,
                               Aufac->p[i].coeffs, Aufac->p[i].length, 0, ctx);
        fmpz_mpoly_scalar_mul_fmpz(fac->poly + i, fac->poly + i, q, ctx);
        fmpz_one(fac->exp + i);
    }

    /*
        At this point, we have 

            Anew = m(x_1, ..., x_n)^(r-1)*A(X, x_1, ...., x_n)
            Aeval[n]   = Anew(X, x_1, ..., x_{n-1}, x_n)  (not stored)
            Aeval[n-1] = Anew(X, x_1, ..., x_{n-1}, α_n)
            ...
            Aeval[1]   = Anew(X, x_1, α_2, ..., α_n)
            Aeval[0]   = Anew(X, α_1, α_2, ..., α_n)

        and the leading coefficients for the factorization of
        Aeval[k] are in Alc[k][i] for 0 <= i < r, 0 <= k <= n

        Choose a prime p that
            1. Divides no coefficient of any Aeval[k], 0 <= k <= n
            2. Divides no coefficient of any Alc[k][i]
            3. Keeps Aeval[0] squarefree
            4. Divides no coefficient of any factor of Aeval[0]

        Try several p before going back an picking a new α.
    */

    p = UWORD(1) << (FLINT_BITS - 1);

next_tuncer_prime:

    if (p >= UWORD_MAX_PRIME)
    {
        success = 0;
        goto cleanup;
    }

    p = n_nextprime(p, 1);
/*
flint_printf("---------- next_tuncer_prime p = %wu -----------\n", p);
    usleep(1000);
*/
    nmod_mpoly_ctx_set_modulus(ctxp, p);
    nmod_mpoly_factor_fit_length(facp, r, ctxp);
    nmod_mpoly_factor_fit_length(tfacp, r, ctxp);
    facp->num = r;
    tfacp->num = r;

    nmod_mpolyv_fit_length(Aevalp, n + 1, ctxp);
    nmod_mpolyv_fit_length(Alcp, (n + 1)*r, ctxp);

    facp->constant = 1;
    for (i = 0; i < r; i++)
    {
        /* check condition 4 */
        _nmod_mpoly_set_fmpz_mpoly(facp->poly + i, ctxp, fac->poly + i, ctx);
        if (facp->poly[i].length != fac->poly[i].length)
            goto next_tuncer_prime;

        fmpz_one(facp->exp + i);
    }

    /* check condition 3 */
    Aup->mod = ctxp->ffinfo->mod;
    nmod_poly_set_fmpz_poly(Aup, Au);
    if (Aup->length != Au->length || !nmod_poly_is_squarefree(Aup))
        goto next_tuncer_prime;

    for (k = 0; k <= n; k++)
    {
        /* check condition 1 */
        _nmod_mpoly_set_fmpz_mpoly(Aevalp->coeffs + k, ctxp,
                                               k < n ? Aevals + k : newA, ctx);
        if (Aevalp->coeffs[k].length != (k < n ? Aevals + k : newA)->length)
            goto next_tuncer_prime;

        /* check condition 2 */
        for (i = 0; i < r; i++)
        {
            _nmod_mpoly_set_fmpz_mpoly(Alcp->coeffs + k*r + i, ctxp,
                                                   Alc->coeffs + k*r + i, ctx);
            if (Alcp->coeffs[k*r + i].length != Alc->coeffs[k*r + i].length)
                goto next_tuncer_prime;
        }
    }

    /* set alpha's mod p */
    for (i = 0; i < n; i++)
        alphap[i] = fmpz_fdiv_ui(alpha + i, ctxp->ffinfo->mod.n);

    for (k = 1; k <= n; k++)
    {
        for (i = 0; i < r; i++)
        {
            fmpz_one(tfacp->exp + i);
            _nmod_mpoly_set_lead0(tfacp->poly + i, facp->poly + i,
                                         Alcp->coeffs + k*r + i, ctxp);
        }

        if (k > 2)
        {
            success = nmod_mpoly_hlift_zippel(r, tfacp->poly, state,
                                    alphap, Aevalp->coeffs + k, degs, k, ctxp);
            if (!success)
                goto next_tuncer_prime;
        }
        else
        {
            success = nmod_mpoly_hlift(k, tfacp->poly, r, alphap,
                                               Aevalp->coeffs + k, degs, ctxp);
            if (!success)
                goto next_tuncer_prime;
        }

        nmod_mpoly_factor_swap(tfacp, facp, ctxp);
    }

    fmpz_mpoly_factor_fit_length(fac, r, ctx);
    fac->num = r;
    for (i = 0; i < r; i++)
    {
        fmpz_one(fac->exp + i);
        _fmpz_mpoly_set_nmod_mpoly_smod(fac->poly + i, ctx,
                                                         facp->poly + i, ctxp);
        _fmpz_mpoly_set_lead0(fac->poly + i, fac->poly + i,
                                                   Alc->coeffs + n*r + i, ctx);
    }

    L = fmpz_clog_ui(facBound, ctxp->ffinfo->mod.n);

    if (0)
    {
        /* recursive approach */
        success = fmpz_mfactor_lift_prime_power(r, fac->poly,
                                 facp->poly, newA, degs, alphap, ctx, ctxp, L);
        FLINT_ASSERT(success == 1 && "don't know want to do");
    }
    else
    {
        /* more zippel :-( */
        success = fmpz_mfactor_lift_prime_power_tuncer(r, fac->poly, state,
                                       facp->poly, newA, alphap, ctx, ctxp, L);
        FLINT_ASSERT(success == 1 && "don't know want to do");
    }

    success = 1;

    if (fmpz_mpoly_is_fmpz(m, ctx))
    {
        for (i = 0; i < Aufac->num; i++)
        {
            _fmpz_vec_content(q, fac->poly[i].coeffs, fac->poly[i].length);
            if (fmpz_sgn(fac->poly[i].coeffs + 0) < 0)
                fmpz_neg(q, q);
            fmpz_mpoly_scalar_divexact_fmpz(fac->poly + i, fac->poly + i, q, ctx);
        }
    }
    else
    {
        fmpz_mpoly_univar_t u;
        fmpz_mpoly_univar_init(u, ctx);
        for (i = 0; i < Aufac->num; i++)
        {
            fmpz_mpoly_to_univar(u, fac->poly + i, 0, ctx);
            success = fmpz_mpoly_univar_content_mpoly(t, u, ctx);
            if (!success)
            {
                fmpz_mpoly_univar_clear(u, ctx);
                goto cleanup;
            }
            success = fmpz_mpoly_divides(fac->poly + i, fac->poly + i, t, ctx);
            FLINT_ASSERT(success);
            FLINT_ASSERT(fac->poly[i].length > 0);
            if (fmpz_sgn(fac->poly[i].coeffs + 0) < 0)
                fmpz_mpoly_neg(fac->poly + i, fac->poly + i, ctx);
        }
        fmpz_mpoly_univar_clear(u, ctx);
    }

cleanup:

    flint_randclear(state);

    fmpz_clear(facBound);
    fmpz_clear(q);

    fmpz_mpoly_clear(Acopy, ctx);
    fmpz_mpoly_clear(m, ctx);
    fmpz_mpoly_clear(mpow, ctx);

    fmpz_mpolyv_clear(Alc, ctx);
    fmpz_mpolyv_clear(lc_divs, ctx);

    fmpz_poly_factor_clear(Aufac);
    fmpz_poly_clear(Au);

    _fmpz_vec_clear(alpha, n);
    flint_free(alphap);
    flint_free(degs);

    for (i = 0; i < n; i++)
        fmpz_mpoly_clear(Aevals + i, ctx);
    flint_free(Aevals);

    fmpz_mpoly_factor_clear(tfac, ctx);
    fmpz_mpoly_clear(t, ctx);

    nmod_mpoly_factor_clear(facp, ctxp);
    nmod_mpoly_factor_clear(tfacp, ctxp);
    nmod_mpolyv_clear(Aevalp, ctxp);
    nmod_mpolyv_clear(Alcp, ctxp);
    nmod_poly_clear(Aup);
    nmod_mpoly_ctx_clear(ctxp);

/*
flint_printf("_try_tuncer returning %d\n", success);
*/
    return success;
}

