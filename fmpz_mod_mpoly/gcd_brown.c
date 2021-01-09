/*
    Copyright (C) 2019 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mod_mpoly.h"
#include "fmpz_mod_mpoly_factor.h"


int fmpz_mod_mpolyn_gcd_brown_bivar(
    fmpz_mod_mpolyn_t G,
    fmpz_mod_mpolyn_t Abar,
    fmpz_mod_mpolyn_t Bbar,
    fmpz_mod_mpolyn_t A,
    fmpz_mod_mpolyn_t B,
    const fmpz_mod_mpoly_ctx_t ctx)
{
    int success;
    slong bound;
    fmpz_t alpha, temp, gammaeval;
    fmpz_mod_poly_t Aeval, Beval, Geval, Abareval, Bbareval;
    fmpz_mod_mpolyn_t T;
    slong deggamma, ldegG, ldegAbar, ldegBbar, ldegA, ldegB;
    fmpz_mod_poly_t cA, cB, cG, cAbar, cBbar, gamma, r;
    fmpz_mod_poly_t modulus, modulus2;
    slong N, off, shift;
    flint_bitcnt_t bits = A->bits;
#if FLINT_WANT_ASSERT
    fmpz_mod_poly_t leadA, leadB;
#endif

    FLINT_ASSERT(A->bits == bits);
    FLINT_ASSERT(B->bits == bits);
    FLINT_ASSERT(G->bits == bits);
    FLINT_ASSERT(Abar->bits == bits);
    FLINT_ASSERT(Bbar->bits == bits);

    fmpz_init(alpha);
    fmpz_init(temp);
    fmpz_init(gammaeval);

#if FLINT_WANT_ASSERT
    fmpz_mod_poly_init(leadA, ctx->ffinfo);
    fmpz_mod_poly_init(leadB, ctx->ffinfo);
    fmpz_mod_poly_set(leadA, fmpz_mod_mpolyn_leadcoeff_poly(A), ctx->ffinfo);
    fmpz_mod_poly_set(leadB, fmpz_mod_mpolyn_leadcoeff_poly(B), ctx->ffinfo);
#endif

    N = mpoly_words_per_exp_sp(bits, ctx->minfo);
    mpoly_gen_offset_shift_sp(&off, &shift, 0, bits, ctx->minfo);

    fmpz_mod_mpolyn_init(T, bits, ctx);
    fmpz_mod_poly_init(r, ctx->ffinfo);
    fmpz_mod_poly_init(cA, ctx->ffinfo);
    fmpz_mod_poly_init(cB, ctx->ffinfo);
    fmpz_mod_poly_init(cG, ctx->ffinfo);
    fmpz_mod_poly_init(cAbar, ctx->ffinfo);
    fmpz_mod_poly_init(cBbar, ctx->ffinfo);
    fmpz_mod_poly_init(gamma, ctx->ffinfo);
    fmpz_mod_poly_init(Aeval, ctx->ffinfo);
    fmpz_mod_poly_init(Beval, ctx->ffinfo);
    fmpz_mod_poly_init(Geval, ctx->ffinfo);
    fmpz_mod_poly_init(Abareval, ctx->ffinfo);
    fmpz_mod_poly_init(Bbareval, ctx->ffinfo);
    fmpz_mod_poly_init(modulus, ctx->ffinfo);
    fmpz_mod_poly_init(modulus2, ctx->ffinfo);

    fmpz_mod_mpolyn_content_poly(cA, A, ctx);
    fmpz_mod_mpolyn_content_poly(cB, B, ctx);
    fmpz_mod_mpolyn_divexact_poly(A, cA, ctx);
    fmpz_mod_mpolyn_divexact_poly(B, cB, ctx);

    fmpz_mod_poly_gcd(cG, cA, cB, ctx->ffinfo);

    fmpz_mod_poly_divrem(cAbar, r, cA, cG, ctx->ffinfo);
    FLINT_ASSERT(fmpz_mod_poly_is_zero(r, ctx->ffinfo));
    fmpz_mod_poly_divrem(cBbar, r, cB, cG, ctx->ffinfo);
    FLINT_ASSERT(fmpz_mod_poly_is_zero(r, ctx->ffinfo));

    fmpz_mod_poly_gcd(gamma, fmpz_mod_mpolyn_leadcoeff_poly(A),
                             fmpz_mod_mpolyn_leadcoeff_poly(B), ctx->ffinfo);

    ldegA = fmpz_mod_mpolyn_lastdeg(A, ctx);
    ldegB = fmpz_mod_mpolyn_lastdeg(B, ctx);
    deggamma = fmpz_mod_poly_degree(gamma, ctx->ffinfo);
    bound = 1 + deggamma + FLINT_MAX(ldegA, ldegB);

    fmpz_mod_poly_set_ui(modulus, 1, ctx->ffinfo);

    fmpz_sub_ui(alpha, fmpz_mod_ctx_modulus(ctx->ffinfo), 1);

choose_prime:

    fmpz_sub_ui(alpha, alpha, 1);
    if (fmpz_sgn(alpha) <= 0)
    {
        success = 0;
        goto cleanup;
    }

    /* make sure evaluation point does not kill both lc(A) and lc(B) */
    fmpz_mod_poly_evaluate_fmpz(gammaeval, gamma, alpha, ctx->ffinfo);
    if (fmpz_is_zero(gammaeval))
    {
        goto choose_prime;
    }

    /* evaluation point should kill neither A nor B */
    fmpz_mod_mpolyn_intp_reduce_sm_poly(Aeval, A, alpha, ctx);
    fmpz_mod_mpolyn_intp_reduce_sm_poly(Beval, B, alpha, ctx);
    FLINT_ASSERT(Aeval->length > 0);
    FLINT_ASSERT(Beval->length > 0);

    fmpz_mod_poly_gcd(Geval, Aeval, Beval, ctx->ffinfo);
    fmpz_mod_poly_divrem(Abareval, r, Aeval, Geval, ctx->ffinfo);
    FLINT_ASSERT(fmpz_mod_poly_is_zero(r, ctx->ffinfo));
    fmpz_mod_poly_divrem(Bbareval, r, Beval, Geval, ctx->ffinfo);
    FLINT_ASSERT(fmpz_mod_poly_is_zero(r, ctx->ffinfo));

    FLINT_ASSERT(Geval->length > 0);
    FLINT_ASSERT(Abareval->length > 0);
    FLINT_ASSERT(Bbareval->length > 0);

    if (fmpz_mod_poly_degree(Geval, ctx->ffinfo) == 0)
    {
        fmpz_mod_mpolyn_one(G, ctx);
        fmpz_mod_mpolyn_swap(Abar, A, ctx);
        fmpz_mod_mpolyn_swap(Bbar, B, ctx);
        goto successful_put_content;    
    }

    if (fmpz_mod_poly_degree(modulus, ctx->ffinfo) > 0)
    {
        FLINT_ASSERT(G->length > 0);
        if (fmpz_mod_poly_degree(Geval, ctx->ffinfo) > ((G->exps + N*0)[off]>>shift))
        {
            goto choose_prime;
        }
        else if (fmpz_mod_poly_degree(Geval, ctx->ffinfo) < ((G->exps + N*0)[off]>>shift))
        {
            fmpz_mod_poly_set_ui(modulus, 1, ctx->ffinfo);
        }
    }

    fmpz_mod_poly_scalar_mul_fmpz(Geval, Geval, gammaeval, ctx->ffinfo);

    if (fmpz_mod_poly_degree(modulus, ctx->ffinfo) > 0)
    {
        fmpz_mod_poly_evaluate_fmpz(temp, modulus, alpha, ctx->ffinfo);
        fmpz_mod_inv(temp, temp, ctx->ffinfo);
        fmpz_mod_poly_scalar_mul_fmpz(modulus, modulus, temp, ctx->ffinfo);
        fmpz_mod_mpolyn_intp_crt_sm_poly(&ldegG, G, T, Geval, modulus, alpha, ctx);
        fmpz_mod_mpolyn_intp_crt_sm_poly(&ldegAbar, Abar, T, Abareval, modulus, alpha, ctx);
        fmpz_mod_mpolyn_intp_crt_sm_poly(&ldegBbar, Bbar, T, Bbareval, modulus, alpha, ctx);
    }
    else
    {
        fmpz_mod_mpolyn_intp_lift_sm_poly(G, Geval, ctx);
        fmpz_mod_mpolyn_intp_lift_sm_poly(Abar, Abareval, ctx);
        fmpz_mod_mpolyn_intp_lift_sm_poly(Bbar, Bbareval, ctx);
        ldegG = 0;
        ldegAbar = 0;
        ldegBbar = 0;
    }

    fmpz_mod_poly_scalar_mul_fmpz(modulus2, modulus, alpha, ctx->ffinfo);
    fmpz_mod_poly_shift_left(modulus, modulus, 1, ctx->ffinfo);
    fmpz_mod_poly_sub(modulus, modulus, modulus2, ctx->ffinfo);

    if (fmpz_mod_poly_degree(modulus, ctx->ffinfo) < bound)
    {
        goto choose_prime;
    }


    FLINT_ASSERT(ldegG >= 0);
    FLINT_ASSERT(ldegAbar >= 0);
    FLINT_ASSERT(ldegBbar >= 0);

    if (deggamma + ldegA == ldegG + ldegAbar &&
        deggamma + ldegB == ldegG + ldegBbar)
    {
        goto successful;
    }

    fmpz_mod_poly_one(modulus, ctx->ffinfo);
    goto choose_prime;

successful:

    fmpz_mod_mpolyn_content_poly(modulus, G, ctx);
    fmpz_mod_mpolyn_divexact_poly(G, modulus, ctx);
    fmpz_mod_mpolyn_divexact_poly(Abar, fmpz_mod_mpolyn_leadcoeff_poly(G), ctx);
    fmpz_mod_mpolyn_divexact_poly(Bbar, fmpz_mod_mpolyn_leadcoeff_poly(G), ctx);

successful_put_content:

    fmpz_mod_mpolyn_mul_poly(G, cG, ctx);
    fmpz_mod_mpolyn_mul_poly(Abar, cAbar, ctx);
    fmpz_mod_mpolyn_mul_poly(Bbar, cBbar, ctx);

    success = 1;

cleanup:

#if FLINT_WANT_ASSERT
    if (success)
    {
        FLINT_ASSERT(fmpz_is_one(fmpz_mod_mpolyn_leadcoeff(G)));
        fmpz_mod_poly_mul(modulus, fmpz_mod_mpolyn_leadcoeff_poly(G),
                            fmpz_mod_mpolyn_leadcoeff_poly(Abar), ctx->ffinfo);
        FLINT_ASSERT(fmpz_mod_poly_equal(modulus, leadA, ctx->ffinfo));
        fmpz_mod_poly_mul(modulus, fmpz_mod_mpolyn_leadcoeff_poly(G),
                            fmpz_mod_mpolyn_leadcoeff_poly(Bbar), ctx->ffinfo);
        FLINT_ASSERT(fmpz_mod_poly_equal(modulus, leadB, ctx->ffinfo));
    }
    fmpz_mod_poly_clear(leadA, ctx->ffinfo);
    fmpz_mod_poly_clear(leadB, ctx->ffinfo);
#endif

    fmpz_mod_poly_clear(r, ctx->ffinfo);
    fmpz_mod_poly_clear(cA, ctx->ffinfo);
    fmpz_mod_poly_clear(cB, ctx->ffinfo);
    fmpz_mod_poly_clear(cG, ctx->ffinfo);
    fmpz_mod_poly_clear(cAbar, ctx->ffinfo);
    fmpz_mod_poly_clear(cBbar, ctx->ffinfo);

    fmpz_mod_poly_clear(gamma, ctx->ffinfo);

    fmpz_mod_poly_clear(Aeval, ctx->ffinfo);
    fmpz_mod_poly_clear(Beval, ctx->ffinfo);
    fmpz_mod_poly_clear(Geval, ctx->ffinfo);
    fmpz_mod_poly_clear(Abareval, ctx->ffinfo);
    fmpz_mod_poly_clear(Bbareval, ctx->ffinfo);

    fmpz_mod_mpolyn_clear(T, ctx);

    fmpz_mod_poly_clear(modulus, ctx->ffinfo);
    fmpz_mod_poly_clear(modulus2, ctx->ffinfo);

    fmpz_clear(alpha);
    fmpz_clear(temp);
    fmpz_clear(gammaeval);

    return success;
}


void fmpz_mod_mpolyn_set(
    fmpz_mod_mpolyn_t A,
    const fmpz_mod_mpolyn_t B,
    const fmpz_mod_mpoly_ctx_t ctx)
{
    slong i;
    slong Blen = B->length;
    slong N = mpoly_words_per_exp(B->bits, ctx->minfo);

    fmpz_mod_mpolyn_fit_bits(A, B->bits, ctx);
    A->bits = B->bits;

    fmpz_mod_mpolyn_fit_length(A, Blen, ctx);

    mpoly_copy_monomials(A->exps, B->exps, Blen, N);

    for (i = 0; i < Blen; i++)
        fmpz_mod_poly_set(A->coeffs + i, B->coeffs + i, ctx->ffinfo);

    A->length = Blen;
}

void fmpz_mod_mpolyn_content_last(
    fmpz_mod_poly_t g,
    const fmpz_mod_mpolyn_t A,
    const fmpz_mod_mpoly_ctx_t ctx)
{
    slong i;
    fmpz_mod_poly_zero(g, ctx->ffinfo);
    for (i = 0; i < A->length; i++)
        fmpz_mod_poly_gcd(g, g, A->coeffs + i, ctx->ffinfo);
}

void fmpz_mod_mpolyn_divexact_last(
    fmpz_mod_mpolyn_t A,
    const fmpz_mod_poly_t g,
    const fmpz_mod_mpoly_ctx_t ctx)
{
    slong i;
    for (i = 0; i < A->length; i++)
        fmpz_mod_poly_div(A->coeffs + i, A->coeffs + i, g, ctx->ffinfo);
}

void fmpz_mod_mpolyn_mul_last(
    fmpz_mod_mpolyn_t A,
    const fmpz_mod_poly_t g,
    const fmpz_mod_mpoly_ctx_t ctx)
{
    slong i;
    for (i = 0; i < A->length; i++)
        fmpz_mod_poly_mul(A->coeffs + i, A->coeffs + i, g, ctx->ffinfo);
}

int fmpz_mod_mpolyn_is_nonzero_fmpz(
    const fmpz_mod_mpolyn_t A,
    const fmpz_mod_mpoly_ctx_t ctx)
{
    slong N;

    if (A->length != 1 || A->coeffs[0].length != 1)
        return 0;

    N = mpoly_words_per_exp_sp(A->bits, ctx->minfo);
    return mpoly_monomial_is_zero(A->exps + N*0, N);
}

int fmpz_mod_mpolyn_is_canonical(
    const fmpz_mod_mpolyn_t A,
    const fmpz_mod_mpoly_ctx_t ctx)
{
    slong i;

    if (!mpoly_monomials_valid_test(A->exps, A->length, A->bits, ctx->minfo))
    {
        return 0;
    }

    if (mpoly_monomials_overflow_test(A->exps, A->length, A->bits, ctx->minfo))
        return 0;

    if (!mpoly_monomials_inorder_test(A->exps, A->length, A->bits, ctx->minfo))
        return 0;

    for (i = 0; i < A->length; i++)
    {
        if (!fmpz_mod_poly_is_canonical(A->coeffs + i, ctx->ffinfo) ||
            fmpz_mod_poly_is_zero(A->coeffs + i, ctx->ffinfo))
        {
            return 0;
        }
    }

    return 1;
}


void fmpz_mod_mpolyn_interp_reduce_2sm_mpolyn(
    fmpz_mod_mpolyn_t E,
    fmpz_mod_mpolyn_t F,
    fmpz_mod_mpolyn_t A,
    slong var,
    fmpz_mod_poly_t alphapow,
    const fmpz_mod_mpoly_ctx_t ctx)
{
    slong N = mpoly_words_per_exp(A->bits, ctx->minfo);
    slong offset, shift, k;
    ulong mask;
    fmpz_t e, f;
    fmpz_mod_poly_struct * Acoeff = A->coeffs;
    ulong * Aexp = A->exps;
    slong Alen = A->length;
    slong Ai;
    fmpz_mod_poly_struct * Ecoeff;
    ulong * Eexp;
    slong Ei;
    fmpz_mod_poly_struct * Fcoeff;
    ulong * Fexp;
    slong Fi;

    fmpz_init(e);
    fmpz_init(f);

    FLINT_ASSERT(var > 0);
    mpoly_gen_offset_shift_sp(&offset, &shift, var - 1, A->bits, ctx->minfo);
    mask = (-UWORD(1)) >> (FLINT_BITS - A->bits);

    Ecoeff = E->coeffs;
    Eexp = E->exps;
    Ei = 0;
    Fcoeff = F->coeffs;
    Fexp = F->exps;
    Fi = 0;
    for (Ai = 0; Ai < Alen; Ai++)
    {
        fmpz_mod_poly_eval2_pow(e, f, Acoeff + Ai, alphapow, ctx->ffinfo);
        k = ((Aexp + N*Ai)[offset] >> shift) & mask;

        if (!fmpz_is_zero(e))
        {
            if (Ei > 0 && mpoly_monomial_equal_extra(Eexp + N*(Ei - 1), Aexp + N*Ai, N, offset, -(k << shift)))
            {
                /* append to previous */
                fmpz_mod_poly_set_coeff_fmpz(Ecoeff + Ei - 1, k, e, ctx->ffinfo);
            }
            else
            {
                FLINT_ASSERT(Ei == 0 || mpoly_monomial_gt_nomask_extra(Eexp + N*(Ei - 1), Aexp + N*Ai, N, offset, -(k << shift)));

                /* create new */
                if (Ei >= E->alloc)
                {
                    fmpz_mod_mpolyn_fit_length(E, Ei + 1, ctx);
                    Ecoeff = E->coeffs;
                    Eexp = E->exps;
                }
                mpoly_monomial_set_extra(Eexp + N*Ei, Aexp + N*Ai, N, offset, -(k << shift));
                fmpz_mod_poly_zero(Ecoeff + Ei, ctx->ffinfo);
                fmpz_mod_poly_set_coeff_fmpz(Ecoeff + Ei, k, e, ctx->ffinfo);
                Ei++;
            }
        }

        if (!fmpz_is_zero(f))
        {
            if (Fi > 0 && mpoly_monomial_equal_extra(Fexp + N*(Fi - 1), Aexp + N*Ai, N, offset, -(k << shift)))
            {
                /* append to previous */
                fmpz_mod_poly_set_coeff_fmpz(Fcoeff + Fi - 1, k, f, ctx->ffinfo);
            }
            else
            {
                FLINT_ASSERT(Fi == 0 || mpoly_monomial_gt_nomask_extra(Fexp + N*(Fi - 1), Aexp + N*Ai, N, offset, -(k << shift)));

                /* create new */
                if (Fi >= F->alloc)
                {
                    fmpz_mod_mpolyn_fit_length(F, Fi + 1, ctx);
                    Fcoeff = F->coeffs;
                    Fexp = F->exps;
                }
                mpoly_monomial_set_extra(Fexp + N*Fi, Aexp + N*Ai, N, offset, -(k << shift));
                fmpz_mod_poly_zero(Fcoeff + Fi, ctx->ffinfo);
                fmpz_mod_poly_set_coeff_fmpz(Fcoeff + Fi, k, f, ctx->ffinfo);
                Fi++;
            }
        }
    }
    E->length = Ei;
    F->length = Fi;

    fmpz_clear(e);
    fmpz_clear(f);
}


void fmpz_mod_mpolyn_interp_lift_2sm_mpolyn(
    slong * lastdeg,
    fmpz_mod_mpolyn_t T,
    fmpz_mod_mpolyn_t A,
    fmpz_mod_mpolyn_t B,
    slong var,
    const fmpz_t alpha,
    const fmpz_mod_mpoly_ctx_t ctx)
{
    slong lastlen = 0;
    slong N = mpoly_words_per_exp_sp(A->bits, ctx->minfo);
    slong offset, shift;

    fmpz_mod_poly_struct * Tcoeff;
    ulong * Texp;
    slong Ti;

    fmpz_mod_poly_struct * Acoeff = A->coeffs;
    slong Alen = A->length;
    ulong * Aexp = A->exps;
    slong Ai, ai;

    fmpz_mod_poly_struct * Bcoeff = B->coeffs;
    slong Blen = B->length;
    ulong * Bexp = B->exps;
    slong Bi, bi;

    fmpz zero0 = 0;
    fmpz * Avalue, * Bvalue;
    fmpz_t u, v, FvalueA, FvalueB;
    int cmp;
    fmpz_t d0;

    fmpz_init(d0);
    fmpz_mod_add(d0, alpha, alpha, ctx->ffinfo);
    fmpz_mod_inv(d0, d0, ctx->ffinfo);

    fmpz_init(u);
    fmpz_init(v);
    fmpz_init(FvalueA);
    fmpz_init(FvalueB);

    FLINT_ASSERT(var > 0);
    FLINT_ASSERT(T->bits == A->bits);

    mpoly_gen_offset_shift_sp(&offset, &shift, var - 1, A->bits, ctx->minfo);

    fmpz_mod_mpolyn_fit_length(T, FLINT_MAX(Alen, Blen), ctx);
    Tcoeff = T->coeffs;
    Texp = T->exps;

    Ti = Ai = Bi = 0;
    ai = (Ai >= Alen) ? 0 : A->coeffs[Ai].length - 1;
    bi = (Bi >= Blen) ? 0 : B->coeffs[Bi].length - 1;

    while (Ai < Alen || Bi < Blen)
    {
        if (Ti >= T->alloc)
        {
            slong extra = FLINT_MAX(Alen - Ai, Blen - Bi);
            fmpz_mod_mpolyn_fit_length(T, Ti + extra + 1, ctx);
            Tcoeff = T->coeffs;
            Texp = T->exps;
        }

        FLINT_ASSERT(Ai >= Alen || !fmpz_is_zero((Acoeff + Ai)->coeffs + ai));
        FLINT_ASSERT(Bi >= Blen || !fmpz_is_zero((Bcoeff + Bi)->coeffs + bi));

        Avalue = &zero0;
        if (Ai < Alen)
        {
            Avalue = (Acoeff + Ai)->coeffs + ai;
            mpoly_monomial_set_extra(Texp + N*Ti,
                                     Aexp + N*Ai, N, offset, ai << shift);
        }

        Bvalue = &zero0;
        if (Bi < Blen)
        {
            cmp = (Avalue == &zero0) ? -1 : mpoly_monomial_cmp_nomask_extra(
                             Texp + N*Ti, Bexp + N*Bi, N, offset, bi << shift);
            if (cmp <= 0)
            {
                Bvalue = (Bcoeff + Bi)->coeffs + bi;
            }
            if (cmp < 0)
            {
                Avalue = &zero0;
                mpoly_monomial_set_extra(Texp + N*Ti,
                                         Bexp + N*Bi, N, offset, bi << shift);
            }
        }

        fmpz_mod_neg(FvalueA, Avalue, ctx->ffinfo);
        fmpz_mod_neg(FvalueB, Bvalue, ctx->ffinfo);
        fmpz_mod_sub(u, FvalueB, FvalueA, ctx->ffinfo);
        fmpz_mod_add(v, FvalueB, FvalueA, ctx->ffinfo);
        fmpz_mod_mul(v, alpha, v, ctx->ffinfo);
        fmpz_mod_neg(v, v, ctx->ffinfo);

        FLINT_ASSERT(!fmpz_is_zero(u) || !fmpz_is_zero(v));
        fmpz_mod_poly_zero(Tcoeff + Ti, ctx->ffinfo);
        fmpz_mod_mul(u, u, d0, ctx->ffinfo);
        fmpz_mod_mul(v, v, d0, ctx->ffinfo);
        fmpz_mod_poly_set_coeff_fmpz(Tcoeff + Ti, 0, v, ctx->ffinfo);
        fmpz_mod_poly_set_coeff_fmpz(Tcoeff + Ti, 1, u, ctx->ffinfo);

        if (Avalue != &zero0)
        {
            do {
                ai--;
            } while (ai >= 0 && fmpz_is_zero((Acoeff + Ai)->coeffs + ai));
            if (ai < 0)
            {
                Ai++;
                if (Ai < Alen)
                    ai = A->coeffs[Ai].length - 1;
            }
        }
        if (Bvalue != &zero0)
        {
            do {
                bi--;
            } while (bi >= 0 && fmpz_is_zero((Bcoeff + Bi)->coeffs + bi));
            if (bi < 0)
            {
                Bi++;
                if (Bi < Blen)
                    bi = B->coeffs[Bi].length - 1;
            }
        }

        FLINT_ASSERT(!fmpz_mod_poly_is_zero(Tcoeff + Ti, ctx->ffinfo));
        lastlen = FLINT_MAX(lastlen, Tcoeff[Ti].length);
        Ti++;
    }
    T->length = Ti;

    *lastdeg = lastlen - 1;

    FLINT_ASSERT(fmpz_mod_mpolyn_is_canonical(T, ctx));

    fmpz_clear(d0);
    fmpz_clear(u);
    fmpz_clear(v);
    fmpz_clear(FvalueA);
    fmpz_clear(FvalueB);

    return;
}


int fmpz_mod_mpolyn_interp_crt_2sm_mpolyn(
    slong * lastdeg,
    fmpz_mod_mpolyn_t F,
    fmpz_mod_mpolyn_t T,
    fmpz_mod_mpolyn_t A,
    fmpz_mod_mpolyn_t B,
    slong var,
    fmpz_mod_poly_t modulus,
    fmpz_mod_poly_t alphapow,
    const fmpz_mod_mpoly_ctx_t ctx)
{
    slong lastlen = 0;
    int changed = 0;
    slong N = mpoly_words_per_exp_sp(A->bits, ctx->minfo);
    slong offset, shift;
    fmpz_mod_poly_t zero;
    fmpz zero0 = 0;

    fmpz_mod_poly_struct * Tcoeff;
    ulong * Texp;
    slong Ti;

    fmpz_mod_poly_struct * Fcoeff = F->coeffs;
    slong Flen = F->length;
    ulong * Fexp = F->exps;
    slong Fi;

    fmpz_mod_poly_struct * Acoeff = A->coeffs;
    slong Alen = A->length;
    ulong * Aexp = A->exps;
    slong Ai, ai;

    fmpz_mod_poly_struct * Bcoeff = B->coeffs;
    slong Blen = B->length;
    ulong * Bexp = B->exps;
    slong Bi, bi;

    fmpz_mod_poly_struct * Fvalue;
    fmpz * Avalue, * Bvalue;
    fmpz_t u, v, FvalueA, FvalueB;
    int texp_set, cmp;

    zero->coeffs = NULL;
    zero->alloc = 0;
    zero->length = 0;

    fmpz_init(u);
    fmpz_init(v);
    fmpz_init(FvalueA);
    fmpz_init(FvalueB);

#if FLINT_WANT_ASSERT
    fmpz_mod_poly_evaluate_fmpz(u, modulus, alphapow->coeffs + 1, ctx->ffinfo);
    fmpz_mod_mul(u, u, alphapow->coeffs + 1, ctx->ffinfo);
    fmpz_mod_add(u, u, u, ctx->ffinfo);
    FLINT_ASSERT(fmpz_is_one(u));
    fmpz_mod_neg(v, alphapow->coeffs + 1, ctx->ffinfo);
    fmpz_mod_poly_evaluate_fmpz(u, modulus, v, ctx->ffinfo);
    fmpz_mod_mul(u, u, alphapow->coeffs + 1, ctx->ffinfo);
    fmpz_mod_add(u, u, u, ctx->ffinfo);
    FLINT_ASSERT(fmpz_is_one(u));
#endif

    FLINT_ASSERT(fmpz_mod_mpolyn_is_canonical(A, ctx));
    FLINT_ASSERT(fmpz_mod_mpolyn_is_canonical(B, ctx));
    FLINT_ASSERT(fmpz_mod_mpolyn_is_canonical(F, ctx));

    FLINT_ASSERT(var > 0);
    FLINT_ASSERT(T->bits == A->bits);
    FLINT_ASSERT(F->bits == A->bits);
    FLINT_ASSERT(A->bits <= FLINT_BITS);

    mpoly_gen_offset_shift_sp(&offset, &shift, var - 1, A->bits, ctx->minfo);

    Flen = F->length;

    fmpz_mod_mpolyn_fit_length(T, FLINT_MAX(Flen, Alen), ctx);
    Tcoeff = T->coeffs;
    Texp = T->exps;

    Ti = Fi = Ai = Bi = 0;
    ai = (Ai >= Alen) ? 0 : A->coeffs[Ai].length - 1;
    bi = (Bi >= Blen) ? 0 : B->coeffs[Bi].length - 1;

    while (Fi < Flen || Ai < Alen || Bi < Blen)
    {
        if (Ti >= T->alloc)
        {
            slong extra = Flen - Fi;
            extra = FLINT_MAX(extra, Alen - Ai);
            extra = FLINT_MAX(extra, Blen - Bi);
            fmpz_mod_mpolyn_fit_length(T, Ti + extra + 1, ctx);
            Tcoeff = T->coeffs;
            Texp = T->exps;
        }

        FLINT_ASSERT(Fi >= Flen || (Fcoeff + Fi)->length != 0);
        FLINT_ASSERT(Ai >= Alen || !fmpz_is_zero((Acoeff + Ai)->coeffs + ai));
        FLINT_ASSERT(Bi >= Blen || !fmpz_is_zero((Bcoeff + Bi)->coeffs + bi));

        Fvalue = zero;
        texp_set = 0;
        if (Fi < Flen)
        {
            Fvalue = Fcoeff + Fi;
            texp_set = 1;
            mpoly_monomial_set(Texp + N*Ti, Fexp + N*Fi, N);
        }

        Avalue = &zero0;
        if (Ai < Alen)
        {
            cmp = (!texp_set) ? -1
                     : mpoly_monomial_cmp_nomask_extra(Texp + N*Ti,
                                      Aexp + N*Ai, N, offset, ai << shift);

            if (cmp <= 0)
            {
                Avalue = (Acoeff + Ai)->coeffs + ai;
            }
            if (cmp < 0)
            {
                Fvalue = zero;
                texp_set = 1;
                mpoly_monomial_set_extra(Texp + N*Ti,
                                         Aexp + N*Ai, N, offset, ai << shift);
            }
        }

        Bvalue = 0;
        if (Bi < Blen)
        {
            cmp = (!texp_set) ? -1 : mpoly_monomial_cmp_nomask_extra(
                             Texp + N*Ti, Bexp + N*Bi, N, offset, bi << shift);

            if (cmp <= 0)
            {
                Bvalue = (Bcoeff + Bi)->coeffs + bi;
            }
            if (cmp < 0)
            {
                Fvalue = zero;
                Avalue = &zero0;
                texp_set = 1;
                mpoly_monomial_set_extra(Texp + N*Ti,
                                         Bexp + N*Bi, N, offset, bi << shift);
            }
        }

        FLINT_ASSERT(texp_set);

        fmpz_mod_poly_eval2_pow(FvalueA, FvalueB, Fvalue, alphapow, ctx->ffinfo);
        fmpz_mod_sub(FvalueA, FvalueA, Avalue, ctx->ffinfo);
        fmpz_mod_sub(FvalueB, FvalueB, Bvalue, ctx->ffinfo);
        fmpz_mod_sub(u, FvalueB, FvalueA, ctx->ffinfo);
        fmpz_mod_add(v, FvalueB, FvalueA, ctx->ffinfo);
        fmpz_mod_mul(v, alphapow->coeffs + 1, v, ctx->ffinfo);
        fmpz_mod_neg(v, v, ctx->ffinfo);
        changed |= !fmpz_is_zero(u) || !fmpz_is_zero(v);
        fmpz_mod_poly_addmul_linear(Tcoeff + Ti, Fvalue, modulus, u, v, ctx->ffinfo);

        Fi += (Fvalue != zero);
        if (Avalue != &zero0)
        {
            do {
                ai--;
            } while (ai >= 0 && fmpz_is_zero((Acoeff + Ai)->coeffs + ai));
            if (ai < 0)
            {
                Ai++;
                if (Ai < Alen)
                    ai = A->coeffs[Ai].length - 1;
            }
        }
        if (Bvalue != &zero0)
        {
            do {
                bi--;
            } while (bi >= 0 && fmpz_is_zero((Bcoeff + Bi)->coeffs + bi));
            if (bi < 0)
            {
                Bi++;
                if (Bi < Blen)
                    bi = B->coeffs[Bi].length - 1;
            }
        }

        FLINT_ASSERT(!fmpz_mod_poly_is_zero(Tcoeff + Ti, ctx->ffinfo));
        lastlen = FLINT_MAX(lastlen, Tcoeff[Ti].length);
        Ti++;
    }
    T->length = Ti;

    *lastdeg = lastlen - 1;

    if (changed)
        fmpz_mod_mpolyn_swap(T, F, ctx);

    fmpz_clear(u);
    fmpz_clear(v);
    fmpz_clear(FvalueA);
    fmpz_clear(FvalueB);

    FLINT_ASSERT(fmpz_mod_mpolyn_is_canonical(F, ctx));

    return changed;
}

int fmpz_mod_mpolyn_gcd_brown_smprime(
    fmpz_mod_mpolyn_t G,
    fmpz_mod_mpolyn_t Abar,
    fmpz_mod_mpolyn_t Bbar,
    fmpz_mod_mpolyn_t A,
    fmpz_mod_mpolyn_t B,
    slong var,
    const fmpz_mod_mpoly_ctx_t ctx,
    const mpoly_gcd_info_t I,
    fmpz_mod_poly_polyun_mpolyn_stack_t St)
{
    int success;
    slong bound;
    slong offset, shift;
    fmpz_t alpha, temp, gammaevalp, gammaevalm;
    fmpz_mod_mpolyn_struct * Aevalp, * Bevalp, * Gevalp, * Abarevalp, * Bbarevalp;
    fmpz_mod_mpolyn_struct * Aevalm, * Bevalm, * Gevalm, * Abarevalm, * Bbarevalm;
    fmpz_mod_mpolyn_struct * T1;
    slong deggamma, ldegG, ldegAbar, ldegBbar, ldegA, ldegB;
    fmpz_mod_poly_struct * cA, * cB, * cG, * cAbar, * cBbar, * gamma;
    fmpz_mod_poly_struct * modulus, * alphapow;
    flint_bitcnt_t bits = A->bits;
    slong N = mpoly_words_per_exp_sp(bits, ctx->minfo);
#if FLINT_WANT_ASSERT
    fmpz_mod_mpolyn_t Aorg, Borg;
    fmpz_mod_poly_t leadA, leadB;
#endif

    FLINT_ASSERT(bits == St->mpolyn_stack->bits);
    FLINT_ASSERT(bits == A->bits);
    FLINT_ASSERT(bits == B->bits);
    FLINT_ASSERT(bits == G->bits);
    FLINT_ASSERT(bits == Abar->bits);
    FLINT_ASSERT(bits == Bbar->bits);

    FLINT_ASSERT(var > 0);

    if (var == 1)
        return fmpz_mod_mpolyn_gcd_brown_bivar(G, Abar, Bbar, A, B, ctx);

    mpoly_gen_offset_shift_sp(&offset, &shift, var - 1, G->bits, ctx->minfo);

#if FLINT_WANT_ASSERT
    fmpz_mod_poly_init(leadA, ctx->ffinfo);
    fmpz_mod_poly_init(leadB, ctx->ffinfo);
    fmpz_mod_poly_set(leadA, fmpz_mod_mpolyn_leadcoeff_poly(A), ctx->ffinfo);
    fmpz_mod_poly_set(leadB, fmpz_mod_mpolyn_leadcoeff_poly(B), ctx->ffinfo);
    fmpz_mod_mpolyn_init(Aorg, bits, ctx);
    fmpz_mod_mpolyn_init(Borg, bits, ctx);
    fmpz_mod_mpolyn_set(Aorg, A, ctx);
    fmpz_mod_mpolyn_set(Borg, B, ctx);
#endif

    fmpz_init(alpha);
    fmpz_init(temp);
    fmpz_init(gammaevalp);
    fmpz_init(gammaevalm);

    fmpz_mod_poly_stack_fit_request(St->poly_stack, 8);
    cA          = fmpz_mod_poly_stack_take_top(St->poly_stack);
    cB          = fmpz_mod_poly_stack_take_top(St->poly_stack);
    cG          = fmpz_mod_poly_stack_take_top(St->poly_stack);
    cAbar       = fmpz_mod_poly_stack_take_top(St->poly_stack);
    cBbar       = fmpz_mod_poly_stack_take_top(St->poly_stack);
    gamma       = fmpz_mod_poly_stack_take_top(St->poly_stack);
    alphapow    = fmpz_mod_poly_stack_take_top(St->poly_stack);
    modulus     = fmpz_mod_poly_stack_take_top(St->poly_stack);

    fmpz_mod_mpolyn_stack_fit_request(St->mpolyn_stack, 11, ctx);
    Aevalp      = fmpz_mod_mpolyn_stack_take_top(St->mpolyn_stack);
    Bevalp      = fmpz_mod_mpolyn_stack_take_top(St->mpolyn_stack);
    Gevalp      = fmpz_mod_mpolyn_stack_take_top(St->mpolyn_stack);
    Abarevalp   = fmpz_mod_mpolyn_stack_take_top(St->mpolyn_stack);
    Bbarevalp   = fmpz_mod_mpolyn_stack_take_top(St->mpolyn_stack);
    Aevalm      = fmpz_mod_mpolyn_stack_take_top(St->mpolyn_stack);
    Bevalm      = fmpz_mod_mpolyn_stack_take_top(St->mpolyn_stack);
    Gevalm      = fmpz_mod_mpolyn_stack_take_top(St->mpolyn_stack);
    Abarevalm   = fmpz_mod_mpolyn_stack_take_top(St->mpolyn_stack);
    Bbarevalm   = fmpz_mod_mpolyn_stack_take_top(St->mpolyn_stack);
    T1          = fmpz_mod_mpolyn_stack_take_top(St->mpolyn_stack);

    fmpz_mod_mpolyn_content_last(cA, A, ctx);
    fmpz_mod_mpolyn_content_last(cB, B, ctx);
    fmpz_mod_mpolyn_divexact_last(A, cA, ctx);
    fmpz_mod_mpolyn_divexact_last(B, cB, ctx);

    fmpz_mod_poly_gcd(cG, cA, cB, ctx->ffinfo);
    fmpz_mod_poly_div(cAbar, cA, cG, ctx->ffinfo);
    fmpz_mod_poly_div(cBbar, cB, cG, ctx->ffinfo);

    fmpz_mod_poly_gcd(gamma, fmpz_mod_mpolyn_leadcoeff_poly(A),
                             fmpz_mod_mpolyn_leadcoeff_poly(B), ctx->ffinfo);

    ldegA = fmpz_mod_mpolyn_lastdeg(A, ctx);
    ldegB = fmpz_mod_mpolyn_lastdeg(B, ctx);
    deggamma = fmpz_mod_poly_degree(gamma, ctx->ffinfo);
    bound = 1 + deggamma + FLINT_MAX(ldegA, ldegB);

    fmpz_mod_poly_fit_length(alphapow, FLINT_MAX(WORD(3), bound + 1), ctx->ffinfo);

    fmpz_mod_poly_one(modulus, ctx->ffinfo);

    fmpz_fdiv_q_2exp(alpha, fmpz_mod_ctx_modulus(ctx->ffinfo), 1);

choose_prime:

    fmpz_sub_ui(alpha, alpha, 1);
    if (fmpz_sgn(alpha) <= 0)
    {
        success = 0;
        goto cleanup;
    }

    FLINT_ASSERT(alphapow->alloc >= 2);
    alphapow->length = 2;
    fmpz_one(alphapow->coeffs + 0);
    fmpz_set(alphapow->coeffs + 1, alpha);

    /* make sure evaluation point does not kill both lc(A) and lc(B) */
    fmpz_mod_poly_eval2_pow(gammaevalp, gammaevalm, gamma, alphapow, ctx->ffinfo);
    if (fmpz_is_zero(gammaevalp) || fmpz_is_zero(gammaevalm))
    {
        goto choose_prime;
    }

    /* evaluation point should kill neither A nor B */
    fmpz_mod_mpolyn_interp_reduce_2sm_mpolyn(Aevalp, Aevalm, A, var, alphapow, ctx);
    fmpz_mod_mpolyn_interp_reduce_2sm_mpolyn(Bevalp, Bevalm, B, var, alphapow, ctx);
    FLINT_ASSERT(Aevalp->length > 0);
    FLINT_ASSERT(Aevalm->length > 0);
    FLINT_ASSERT(Bevalp->length > 0);
    FLINT_ASSERT(Bevalm->length > 0);

    success = fmpz_mod_mpolyn_gcd_brown_smprime(Gevalp, Abarevalp, Bbarevalp,
                                          Aevalp, Bevalp, var - 1, ctx, I, St);
    success = success &&
              fmpz_mod_mpolyn_gcd_brown_smprime(Gevalm, Abarevalm, Bbarevalm,
                                          Aevalm, Bevalm, var - 1, ctx, I, St);
    if (success == 0)
    {
        goto choose_prime;
    }

    FLINT_ASSERT(Gevalp->length > 0);
    FLINT_ASSERT(Abarevalp->length > 0);
    FLINT_ASSERT(Bbarevalp->length > 0);
    FLINT_ASSERT(Gevalm->length > 0);
    FLINT_ASSERT(Abarevalm->length > 0);
    FLINT_ASSERT(Bbarevalm->length > 0);

    if (fmpz_mod_mpolyn_is_nonzero_fmpz(Gevalp, ctx) ||
        fmpz_mod_mpolyn_is_nonzero_fmpz(Gevalm, ctx))
    {
        fmpz_mod_mpolyn_one(G, ctx);
        fmpz_mod_mpolyn_swap(Abar, A, ctx);
        fmpz_mod_mpolyn_swap(Bbar, B, ctx);
        goto successful_put_content;
    }

    if (Gevalp->coeffs[0].length != Gevalm->coeffs[0].length ||
        !mpoly_monomial_equal(Gevalp->exps + N*0, Gevalm->exps + N*0, N))
    {
        goto choose_prime;
    }

    /* the Geval have matching degrees */
    if (fmpz_mod_poly_degree(modulus, ctx->ffinfo) > 0)
    {
        slong k = fmpz_mod_poly_degree(Gevalp->coeffs + 0, ctx->ffinfo);
        int cmp = mpoly_monomial_cmp_nomask_extra(G->exps + N*0,
                                    Gevalp->exps + N*0, N, offset, k << shift);
        if (cmp < 0)
        {
            goto choose_prime;
        }
        else if (cmp > 0)
        {
            fmpz_mod_poly_one(modulus, ctx->ffinfo);
        }
    }

    /* update interpolants */
    fmpz_mod_mpolyn_scalar_mul_fmpz_mod(Gevalp, gammaevalp, ctx);
    fmpz_mod_mpolyn_scalar_mul_fmpz_mod(Gevalm, gammaevalm, ctx);
    if (fmpz_mod_poly_degree(modulus, ctx->ffinfo) > 0)
    {
        fmpz_mod_poly_evaluate_fmpz(temp, modulus, alpha, ctx->ffinfo);
        fmpz_mod_mul(temp, temp, alpha, ctx->ffinfo);
        fmpz_mod_add(temp, temp, temp, ctx->ffinfo);
        fmpz_mod_poly_scalar_div_fmpz(modulus, modulus, temp, ctx->ffinfo);

        fmpz_mod_mpolyn_interp_crt_2sm_mpolyn(&ldegG, G, T1,
                                  Gevalp, Gevalm, var, modulus, alphapow, ctx);
        fmpz_mod_mpolyn_interp_crt_2sm_mpolyn(&ldegAbar, Abar, T1,
                            Abarevalp, Abarevalm, var, modulus, alphapow, ctx);
        fmpz_mod_mpolyn_interp_crt_2sm_mpolyn(&ldegBbar, Bbar, T1,
                            Bbarevalp, Bbarevalm, var, modulus, alphapow, ctx);
    }
    else
    {
        fmpz_mod_mpolyn_interp_lift_2sm_mpolyn(&ldegG, G,
                                              Gevalp, Gevalm, var, alpha, ctx);
        fmpz_mod_mpolyn_interp_lift_2sm_mpolyn(&ldegAbar, Abar,
                                        Abarevalp, Abarevalm, var, alpha, ctx);
        fmpz_mod_mpolyn_interp_lift_2sm_mpolyn(&ldegBbar, Bbar,
                                        Bbarevalp, Bbarevalm, var, alpha, ctx);
    }

    fmpz_mod_mul(temp, alpha, alpha, ctx->ffinfo);
    fmpz_mod_neg(temp, temp, ctx->ffinfo);
    fmpz_mod_poly_shift_left_scalar_addmul_fmpz_mod(modulus, 2, temp, ctx->ffinfo);

    if (fmpz_mod_poly_degree(modulus, ctx->ffinfo) < bound)
    {
        goto choose_prime;
    }

    FLINT_ASSERT(ldegG >= 0);
    FLINT_ASSERT(ldegAbar >= 0);
    FLINT_ASSERT(ldegBbar >= 0);

    if (deggamma + ldegA == ldegG + ldegAbar &&
        deggamma + ldegB == ldegG + ldegBbar)
    {
        goto successful;
    }

    fmpz_mod_poly_one(modulus, ctx->ffinfo);
    goto choose_prime;

successful:

    fmpz_mod_mpolyn_content_last(modulus, G, ctx);
    fmpz_mod_mpolyn_divexact_last(G, modulus, ctx);
    fmpz_mod_mpolyn_divexact_last(Abar, fmpz_mod_mpolyn_leadcoeff_poly(G), ctx);
    fmpz_mod_mpolyn_divexact_last(Bbar, fmpz_mod_mpolyn_leadcoeff_poly(G), ctx);

successful_put_content:

    fmpz_mod_mpolyn_mul_last(G, cG, ctx);
    fmpz_mod_mpolyn_mul_last(Abar, cAbar, ctx);
    fmpz_mod_mpolyn_mul_last(Bbar, cBbar, ctx);

    success = 1;

cleanup:

#if FLINT_WANT_ASSERT
    if (success)
    {
        FLINT_ASSERT(fmpz_is_one(fmpz_mod_mpolyn_leadcoeff(G)));
        fmpz_mod_poly_mul(modulus, fmpz_mod_mpolyn_leadcoeff_poly(G),
                            fmpz_mod_mpolyn_leadcoeff_poly(Abar), ctx->ffinfo);
        FLINT_ASSERT(fmpz_mod_poly_equal(modulus, leadA, ctx->ffinfo));
        fmpz_mod_poly_mul(modulus, fmpz_mod_mpolyn_leadcoeff_poly(G),
                            fmpz_mod_mpolyn_leadcoeff_poly(Bbar), ctx->ffinfo);
        FLINT_ASSERT(fmpz_mod_poly_equal(modulus, leadB, ctx->ffinfo));
    }
    fmpz_mod_poly_clear(leadA, ctx->ffinfo);
    fmpz_mod_poly_clear(leadB, ctx->ffinfo);
    fmpz_mod_mpolyn_clear(Aorg, ctx);
    fmpz_mod_mpolyn_clear(Borg, ctx);
#endif

    fmpz_clear(alpha);
    fmpz_clear(temp);
    fmpz_clear(gammaevalp);
    fmpz_clear(gammaevalm);

    fmpz_mod_poly_stack_give_back(St->poly_stack, 8);
    fmpz_mod_mpolyn_stack_give_back(St->mpolyn_stack, 11);

FLINT_ASSERT(success);

    return success;
}


/* should find its way back here in interesting cases */
int fmpz_mod_mpoly_gcd_brown(
    fmpz_mod_mpoly_t G,
    const fmpz_mod_mpoly_t A,
    const fmpz_mod_mpoly_t B,
    const fmpz_mod_mpoly_ctx_t ctx)
{
    if (fmpz_mod_mpoly_is_zero(A, ctx) || fmpz_mod_mpoly_is_zero(B, ctx))
        return fmpz_mod_mpoly_gcd(G, A, B, ctx);

    return _fmpz_mod_mpoly_gcd_algo(G, NULL, NULL, A, B, ctx, MPOLY_GCD_USE_BROWN);
}
