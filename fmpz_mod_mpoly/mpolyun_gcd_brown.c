/*
    Copyright (C) 2019 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fmpz_mod_mpoly.h"

/*
    E = A(v = alpha)
    A is in R[X][v]
    E is in R[X]
*/
void fmpz_mod_mpolyn_intp_reduce_sm_poly(
    fmpz_mod_poly_t E,
    const fmpz_mod_mpolyn_t A,
    const fmpz_t alpha,
    const fmpz_mod_mpoly_ctx_t ctx)
{
    fmpz_t v;
    slong Ai, Alen, k;
    fmpz_mod_poly_struct * Acoeff;
    ulong * Aexp;
    slong N, off, shift;

    fmpz_init(v);

    N = mpoly_words_per_exp_sp(A->bits, ctx->minfo);
    mpoly_gen_offset_shift_sp(&off, &shift, 0, A->bits, ctx->minfo);

    Acoeff = A->coeffs;
    Aexp = A->exps;
    Alen = A->length;
    Ai = 0;
    fmpz_mod_poly_zero(E);
    for (Ai = 0; Ai < Alen; Ai++)
    {
        fmpz_mod_poly_evaluate_fmpz(v, Acoeff + Ai, alpha);
        k = (Aexp + N*Ai)[off] >> shift;
        fmpz_mod_poly_set_coeff_fmpz(E, k, v);
    }

    fmpz_clear(v);
}

/*
    A = B
    A, B are in R[X]
*/
void fmpz_mod_mpolyn_intp_lift_sm_poly(
    fmpz_mod_mpolyn_t A,
    const fmpz_mod_poly_t B,
    const fmpz_mod_mpoly_ctx_t ctx)
{
    slong Bi;
    slong Blen = fmpz_mod_poly_length(B);
    fmpz * Bcoeff = B->coeffs;
    fmpz_mod_poly_struct * Acoeff;
    ulong * Aexp;
    slong Ai;
    slong N, off, shift;

    N = mpoly_words_per_exp_sp(A->bits, ctx->minfo);
    mpoly_gen_offset_shift_sp(&off, &shift, 0, A->bits, ctx->minfo);

    fmpz_mod_mpolyn_fit_length(A, Blen, ctx);
    Acoeff = A->coeffs;
    Aexp = A->exps;

    Ai = 0;
    for (Bi = Blen - 1; Bi >= 0; Bi--)
    {
        if (!fmpz_is_zero(Bcoeff + Bi))
        {
            FLINT_ASSERT(Ai < A->alloc);

            fmpz_mod_poly_set_fmpz(Acoeff + Ai, Bcoeff + Bi);
            mpoly_monomial_zero(Aexp + N*Ai, N);
            (Aexp + N*Ai)[off] = Bi << shift;
            Ai++;
        }
    }
    A->length = Ai;
}


/*
    F = F + modulus*(A - F(v = alpha))
    no assumptions about matching monomials
    F is in Fp[X][v]
    A is in Fp[X]
    it is expected that modulus(alpha) == 1
*/
int fmpz_mod_mpolyn_intp_crt_sm_poly(
    slong * lastdeg_,
    fmpz_mod_mpolyn_t F,
    fmpz_mod_mpolyn_t T,
    const fmpz_mod_poly_t A,
    const fmpz_mod_poly_t modulus,
    const fmpz_t alpha,
    const fmpz_mod_mpoly_ctx_t ctx)
{
    int changed = 0;
    slong lastdeg = -WORD(1);
    fmpz_t u, v;
    slong Fi, Ti, Ai;
    fmpz * Acoeff = A->coeffs;
    slong Flen = F->length;
    fmpz_mod_poly_struct * Fcoeff = F->coeffs;
    ulong * Fexp = F->exps;
    fmpz_mod_poly_struct * Tcoeff;
    ulong * Texp;
    fmpz_mod_poly_t tp;
    slong N, off, shift;

    N = mpoly_words_per_exp_sp(F->bits, ctx->minfo);
    mpoly_gen_offset_shift_sp(&off, &shift, 0, F->bits, ctx->minfo);

    Fi = 0;
    Ai = fmpz_mod_poly_degree(A);

    fmpz_init(u);
    fmpz_init(v);
    fmpz_mod_poly_init(tp, fmpz_mod_ctx_modulus(ctx->ffinfo));

    fmpz_mod_mpolyn_fit_length(T, Flen + Ai + 1, ctx);
    Tcoeff = T->coeffs;
    Texp = T->exps;
    Ti = 0;

    while (Fi < Flen || Ai >= 0)
    {
        FLINT_ASSERT(Ti < T->alloc);

        if (Fi < Flen)
        {
            FLINT_ASSERT(!fmpz_mod_poly_is_zero(Fcoeff + Fi));
            FLINT_ASSERT(fmpz_mod_poly_degree(Fcoeff + Fi) < fmpz_mod_poly_degree(modulus));
        }

        if (Ai >= 0)
        {
            FLINT_ASSERT(!fmpz_is_zero(Acoeff + Ai));
        }

        mpoly_monomial_zero(Texp + N*Ti, N);

        if (Fi < Flen && Ai >= 0 && ((Fexp + N*Fi)[off]>>shift) == Ai)
        {
            /* F term ok, A term ok */
            fmpz_mod_poly_evaluate_fmpz(u, Fcoeff + Fi, alpha);
            fmpz_mod_sub(v, Acoeff + Ai, u, ctx->ffinfo);
            if (!fmpz_is_zero(v))
            {
                changed = 1;
                fmpz_mod_poly_scalar_mul_fmpz(tp, modulus, v);
                fmpz_mod_poly_add(Tcoeff + Ti, Fcoeff + Fi, tp);
            }
            else
            {
                fmpz_mod_poly_set(Tcoeff + Ti, Fcoeff + Fi);
            }

            (Texp + N*Ti)[off] = Ai << shift;
            Fi++;
            do {
                Ai--;
            } while (Ai >= 0 && fmpz_is_zero(Acoeff + Ai));
        }
        else if (Fi < Flen && (Ai < 0 || ((Fexp + N*Fi)[off]>>shift) > Ai))
        {
            /* F term ok, A term missing */
            fmpz_mod_poly_evaluate_fmpz(v, Fcoeff + Fi, alpha);
            if (!fmpz_is_zero(v))
            {
                changed = 1;
                fmpz_mod_poly_scalar_mul_fmpz(tp, modulus, v);
                fmpz_mod_poly_sub(Tcoeff + Ti, Fcoeff + Fi, tp);
            }
            else
            {
                fmpz_mod_poly_set(Tcoeff + Ti, Fcoeff + Fi);                
            }

            (Texp + N*Ti)[off] = (Fexp + N*Fi)[off];
            Fi++;
        }
        else if (Ai >= 0 && (Fi >= Flen || ((Fexp + N*Fi)[off]>>shift) < Ai))
        {
            /* F term missing, A term ok */
            changed = 1;
            fmpz_mod_poly_scalar_mul_fmpz(Tcoeff + Ti, modulus, Acoeff + Ai);

            (Texp + N*Ti)[off] = Ai << shift;
            do {
                Ai--;
            } while (Ai >= 0 && fmpz_is_zero(Acoeff + Ai));
        }
        else
        {
            FLINT_ASSERT(0);
        }

        FLINT_ASSERT(!fmpz_mod_poly_is_zero(Tcoeff + Ti));
        lastdeg = FLINT_MAX(lastdeg, fmpz_mod_poly_degree(Tcoeff + Ti));

        Ti++;
    }
    T->length = Ti;

    fmpz_clear(u);
    fmpz_clear(v);
    fmpz_mod_poly_clear(tp);

    if (changed)
    {
        fmpz_mod_mpolyn_swap(T, F, ctx);
    }

    *lastdeg_ = lastdeg;
    return changed;
}


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
#if WANT_ASSERT
    fmpz_mod_poly_t leadA, leadB;
#endif

    fmpz_init(alpha);
    fmpz_init(temp);
    fmpz_init(gammaeval);
/*
printf("fmpz_mod_mpolyn_gcd_bivar called\n");
printf("A: "); fmpz_mod_mpolyn_print_pretty(A, NULL, ctx); printf("\n");
printf("B: "); fmpz_mod_mpolyn_print_pretty(B, NULL, ctx); printf("\n");
*/

#if WANT_ASSERT
    fmpz_mod_poly_init(leadA, fmpz_mod_ctx_modulus(ctx->ffinfo));
    fmpz_mod_poly_init(leadB, fmpz_mod_ctx_modulus(ctx->ffinfo));
    fmpz_mod_poly_set(leadA, fmpz_mod_mpolyn_leadcoeff_poly(A, ctx));
    fmpz_mod_poly_set(leadB, fmpz_mod_mpolyn_leadcoeff_poly(B, ctx));
#endif

    N = mpoly_words_per_exp_sp(A->bits, ctx->minfo);
    mpoly_gen_offset_shift_sp(&off, &shift, 0, A->bits, ctx->minfo);

    fmpz_mod_poly_init(r, fmpz_mod_ctx_modulus(ctx->ffinfo));
    fmpz_mod_poly_init(cA, fmpz_mod_ctx_modulus(ctx->ffinfo));
    fmpz_mod_poly_init(cB, fmpz_mod_ctx_modulus(ctx->ffinfo));
    fmpz_mod_poly_init(cG, fmpz_mod_ctx_modulus(ctx->ffinfo));
    fmpz_mod_poly_init(cAbar, fmpz_mod_ctx_modulus(ctx->ffinfo));
    fmpz_mod_poly_init(cBbar, fmpz_mod_ctx_modulus(ctx->ffinfo));
    fmpz_mod_poly_init(gamma, fmpz_mod_ctx_modulus(ctx->ffinfo));
    fmpz_mod_poly_init(Aeval, fmpz_mod_ctx_modulus(ctx->ffinfo));
    fmpz_mod_poly_init(Beval, fmpz_mod_ctx_modulus(ctx->ffinfo));
    fmpz_mod_poly_init(Geval, fmpz_mod_ctx_modulus(ctx->ffinfo));
    fmpz_mod_poly_init(Abareval, fmpz_mod_ctx_modulus(ctx->ffinfo));
    fmpz_mod_poly_init(Bbareval, fmpz_mod_ctx_modulus(ctx->ffinfo));
    fmpz_mod_poly_init(modulus, fmpz_mod_ctx_modulus(ctx->ffinfo));
    fmpz_mod_poly_init(modulus2, fmpz_mod_ctx_modulus(ctx->ffinfo));
    fmpz_mod_mpolyn_init(T, A->bits, ctx);

    fmpz_mod_mpolyn_content_poly(cA, A, ctx);
    fmpz_mod_mpolyn_content_poly(cB, B, ctx);
    fmpz_mod_mpolyn_divexact_poly(A, cA, ctx);
    fmpz_mod_mpolyn_divexact_poly(B, cB, ctx);

    fmpz_mod_poly_gcd(cG, cA, cB);

    fmpz_mod_poly_divrem(cAbar, r, cA, cG); FLINT_ASSERT(fmpz_mod_poly_is_zero(r));
    fmpz_mod_poly_divrem(cBbar, r, cB, cG); FLINT_ASSERT(fmpz_mod_poly_is_zero(r));

    fmpz_mod_poly_gcd(gamma, fmpz_mod_mpolyn_leadcoeff_poly(A, ctx),
                             fmpz_mod_mpolyn_leadcoeff_poly(B, ctx));

    ldegA = fmpz_mod_mpolyn_lastdeg(A, ctx);
    ldegB = fmpz_mod_mpolyn_lastdeg(B, ctx);
    deggamma = fmpz_mod_poly_degree(gamma);
    bound = 1 + deggamma + FLINT_MAX(ldegA, ldegB);

    fmpz_mod_poly_set_ui(modulus, 1);

    fmpz_sub_ui(alpha, fmpz_mod_ctx_modulus(ctx->ffinfo), 1);

choose_prime: /* prime is v - alpha */

    fmpz_sub_ui(alpha, alpha, 1);
    if (fmpz_sgn(alpha) <= 0)
    {
        success = 0;
        goto cleanup;
    }

    /* make sure evaluation point does not kill both lc(A) and lc(B) */
    fmpz_mod_poly_evaluate_fmpz(gammaeval, gamma, alpha);
    if (fmpz_is_zero(gammaeval))
    {
        goto choose_prime;
    }

    /* evaluation point should kill neither A nor B */
    fmpz_mod_mpolyn_intp_reduce_sm_poly(Aeval, A, alpha, ctx);
    fmpz_mod_mpolyn_intp_reduce_sm_poly(Beval, B, alpha, ctx);
    FLINT_ASSERT(Aeval->length > 0);
    FLINT_ASSERT(Beval->length > 0);

    fmpz_mod_poly_gcd(Geval, Aeval, Beval);
    fmpz_mod_poly_divrem(Abareval, r, Aeval, Geval); FLINT_ASSERT(fmpz_mod_poly_is_zero(r));
    fmpz_mod_poly_divrem(Bbareval, r, Beval, Geval); FLINT_ASSERT(fmpz_mod_poly_is_zero(r));

    FLINT_ASSERT(Geval->length > 0);
    FLINT_ASSERT(Abareval->length > 0);
    FLINT_ASSERT(Bbareval->length > 0);

    if (fmpz_mod_poly_degree(Geval) == 0)
    {
        fmpz_mod_mpolyn_one(G, ctx);
        fmpz_mod_mpolyn_swap(Abar, A, ctx);
        fmpz_mod_mpolyn_swap(Bbar, B, ctx);
        goto successful_put_content;    
    }

    if (fmpz_mod_poly_degree(modulus) > 0)
    {
        FLINT_ASSERT(G->length > 0);
        if (fmpz_mod_poly_degree(Geval) > ((G->exps + N*0)[off]>>shift))
        {
            goto choose_prime;
        }
        else if (fmpz_mod_poly_degree(Geval) < ((G->exps + N*0)[off]>>shift))
        {
            fmpz_mod_poly_set_ui(modulus, 1);
        }
    }

    fmpz_mod_poly_scalar_mul_fmpz(Geval, Geval, gammaeval);

    if (fmpz_mod_poly_degree(modulus) > 0)
    {
        fmpz_mod_poly_evaluate_fmpz(temp, modulus, alpha);
        fmpz_mod_inv(temp, temp, ctx->ffinfo);
        fmpz_mod_poly_scalar_mul_fmpz(modulus, modulus, temp);
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

    fmpz_mod_poly_scalar_mul_fmpz(modulus2, modulus, alpha);
    fmpz_mod_poly_shift_left(modulus, modulus, 1);
    fmpz_mod_poly_sub(modulus, modulus, modulus2);

    if (fmpz_mod_poly_degree(modulus) < bound)
    {
        goto choose_prime;
    }


    FLINT_ASSERT(ldegG >= 0);
    FLINT_ASSERT(ldegAbar >= 0);
    FLINT_ASSERT(ldegBbar >= 0);

    if (   deggamma + ldegA == ldegG + ldegAbar
        && deggamma + ldegB == ldegG + ldegBbar )
    {
        goto successful;
    }

    fmpz_mod_poly_set_ui(modulus, 1);
    goto choose_prime;

successful:

    fmpz_mod_mpolyn_content_poly(modulus, G, ctx);
    fmpz_mod_mpolyn_divexact_poly(G, modulus, ctx);
    fmpz_mod_mpolyn_divexact_poly(Abar, fmpz_mod_mpolyn_leadcoeff_poly(G, ctx), ctx);
    fmpz_mod_mpolyn_divexact_poly(Bbar, fmpz_mod_mpolyn_leadcoeff_poly(G, ctx), ctx);

successful_put_content:

    fmpz_mod_mpolyn_mul_poly(G, cG, ctx);
    fmpz_mod_mpolyn_mul_poly(Abar, cAbar, ctx);
    fmpz_mod_mpolyn_mul_poly(Bbar, cBbar, ctx);

    success = 1;

cleanup:

#if WANT_ASSERT
    if (success)
    {
        FLINT_ASSERT(fmpz_is_one(fmpz_mod_mpolyn_leadcoeff_last_ref(G, ctx)));
        fmpz_mod_poly_mul(modulus, fmpz_mod_mpolyn_leadcoeff_poly(G, ctx),
                                   fmpz_mod_mpolyn_leadcoeff_poly(Abar, ctx));
        FLINT_ASSERT(fmpz_mod_poly_equal(modulus, leadA));
        fmpz_mod_poly_mul(modulus, fmpz_mod_mpolyn_leadcoeff_poly(G, ctx),
                                   fmpz_mod_mpolyn_leadcoeff_poly(Bbar, ctx));
        FLINT_ASSERT(fmpz_mod_poly_equal(modulus, leadB));
    }
    fmpz_mod_poly_clear(leadA);
    fmpz_mod_poly_clear(leadB);
#endif

    fmpz_mod_poly_clear(r);
    fmpz_mod_poly_clear(cA);
    fmpz_mod_poly_clear(cB);
    fmpz_mod_poly_clear(cG);
    fmpz_mod_poly_clear(cAbar);
    fmpz_mod_poly_clear(cBbar);

    fmpz_mod_poly_clear(gamma);

    fmpz_mod_poly_clear(Aeval);
    fmpz_mod_poly_clear(Beval);
    fmpz_mod_poly_clear(Geval);
    fmpz_mod_poly_clear(Abareval);
    fmpz_mod_poly_clear(Bbareval);

    fmpz_mod_mpolyn_clear(T, ctx);

    fmpz_mod_poly_clear(modulus);
    fmpz_mod_poly_clear(modulus2);

    fmpz_clear(alpha);
    fmpz_clear(temp);
    fmpz_clear(gammaeval);

    return success;
}
