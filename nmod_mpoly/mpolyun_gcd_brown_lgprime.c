/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "nmod_mpoly.h"
#include "fq_nmod_mpoly.h"
int usleep(ulong usec);
int nmod_mpolyun_is_canonical(const nmod_mpolyun_t A, const nmod_mpoly_ctx_t ctx);
mp_limb_t nmod_mpolyun_leadcoeff_last(nmod_mpolyun_t A, const nmod_mpoly_ctx_t ctx);

/*
    E = A mod alpha(v)
    A is in Fp[X][][v]
    E is in (Fp/alpha(v))[X]
*/
void nmod_mpolyun_reduce_last_fq_nmod_poly(
                             fq_nmod_poly_t E, const fq_nmod_ctx_t fqctx,
                            const nmod_mpolyun_t A, const nmod_mpoly_ctx_t ctx)
{
    fq_nmod_t v;
    slong Ai, Alen;
    nmod_mpolyn_struct * Acoeff;
    ulong * Aexp;

    fq_nmod_init(v, fqctx);

    Acoeff = A->coeffs;
    Aexp = A->exps;
    Alen = A->length;

    fq_nmod_poly_zero(E, fqctx);
    for (Ai = 0; Ai < Alen; Ai++)
    {
        nmod_poly_rem(v, (Acoeff + Ai)->coeffs + 0, fqctx->modulus);
        fq_nmod_poly_set_coeff(E, Aexp[Ai], v, fqctx);
    }

    fq_nmod_clear(v, fqctx);
}



/*
    Convert B to A using the lowest degree representative
    A is in           Fp [X][][v]
    B is in (Fp/alpha(v))[X]
*/
void nmod_mpolyun_startinterp_fq_nmod_poly(slong * lastdeg_,
                             nmod_mpolyun_t A, const nmod_mpoly_ctx_t ctx,
                             const fq_nmod_poly_t B, const fq_nmod_ctx_t fqctx)
{
    slong N = mpoly_words_per_exp_sp(A->bits, ctx->minfo);
    slong Bexp;
    slong Blen = fq_nmod_poly_length(B, fqctx);
    fq_nmod_struct * Bcoeff = B->coeffs;
    nmod_mpolyn_struct * Acoeff;
    ulong * Aexp;
    slong Ai;
    slong lastdeg = -WORD(1);

    nmod_mpolyun_fit_length(A, Blen, ctx);
    Acoeff = A->coeffs;
    Aexp = A->exps;

    Ai = 0;
    for (Bexp = Blen - 1; Bexp >= 0; Bexp--)
    {
        if (!fq_nmod_is_zero(Bcoeff + Bexp, fqctx))
        {
            FLINT_ASSERT(Ai < A->alloc);

            nmod_mpolyn_fit_length(Acoeff + Ai, 1, ctx);
            mpoly_monomial_zero((Acoeff + Ai)->exps + N*0, N);
            nmod_poly_set((Acoeff + Ai)->coeffs + 0, Bcoeff + Bexp);
            lastdeg = FLINT_MAX(lastdeg, nmod_poly_degree(Bcoeff + Bexp));
            Aexp[Ai] = Bexp;
            (Acoeff + Ai)->length = 1;
            Ai++;
        }
    }
    A->length = Ai;

    *lastdeg_ = lastdeg;
}


/*
    F = F + modulus*((A - F(alpha))/(modulus(alpha)))
    no assumptions about matching monomials
    F is in Fp[X][v]
    A is in (Fp/alpha(v))[X]    alpha(v) is fqctx->modulus(v)
*/
int nmod_mpolyun_addinterp_fq_nmod_poly(slong * lastdeg_,
                         nmod_mpolyun_t F, nmod_mpolyun_t T, nmod_poly_t modulus,
                         const nmod_mpoly_ctx_t ctx,
                                   fq_nmod_poly_t A, const fq_nmod_ctx_t fqctx)


{
    int changed = 0;
    slong lastdeg = -WORD(1);
    slong N = mpoly_words_per_exp_sp(F->bits, ctx->minfo);
    fq_nmod_t u, v;
    nmod_poly_t w;
    slong Fi, Toff, Aexp;
    fq_nmod_struct * Acoeff = A->coeffs;
    slong Flen = F->length;
    nmod_mpolyn_struct * Fcoeff = F->coeffs;
    ulong * Fexp = F->exps;
    nmod_mpolyn_struct * Tcoeff;
    ulong * Texp;
    nmod_poly_t tp;
    fq_nmod_t inv_m_eval;

    fq_nmod_init(inv_m_eval, fqctx);
    nmod_poly_rem(inv_m_eval, modulus, fqctx->modulus);
    fq_nmod_inv(inv_m_eval, inv_m_eval, fqctx);



    Fi = 0;

/*
flint_printf("\nnmod_mpolyun_addinterp_bivar: (alpha = %wu)\n", alpha);
flint_printf("A: "); nmod_poly_print_pretty(A, "X"); printf("\n");
flint_printf("F: "); nmod_mpolyun_print_pretty(F, NULL, ctx); printf("\n");
*/

    fq_nmod_init(u, fqctx);
    fq_nmod_init(v, fqctx);
    nmod_poly_init(w, fqctx->modulus->mod.n);

    Aexp = fq_nmod_poly_degree(A, fqctx);

    nmod_poly_init(tp, ctx->ffinfo->mod.n);

    nmod_mpolyun_fit_length(T, Flen + Aexp + 1, ctx);
    Tcoeff = T->coeffs;
    Texp = T->exps;
    Toff = 0;

    while (Fi < Flen || Aexp >= 0)
    {
        FLINT_ASSERT(Toff < T->alloc);
/*
flint_printf("Fi: %wd\n",Fi);
flint_printf("Aexp: %wd\n",Aexp);
*/
        if (Fi < Flen)
        {
            FLINT_ASSERT(!nmod_poly_is_zero((Fcoeff + Fi)->coeffs + 0));
            FLINT_ASSERT(nmod_poly_degree((Fcoeff + Fi)->coeffs + 0) < nmod_poly_degree(modulus));
        }

        if (Aexp >= 0)
        {
            FLINT_ASSERT(!fq_nmod_is_zero(Acoeff + Aexp, fqctx));
        }

        nmod_mpolyn_fit_length(Tcoeff + Toff, 1, ctx);

        if (Fi < Flen && Aexp >= 0 && Fexp[Fi] == Aexp)
        {
            /* F term ok, A term ok */
            nmod_poly_rem(u, (Fcoeff + Fi)->coeffs + 0, fqctx->modulus);
            fq_nmod_sub(v, Acoeff + Aexp, u, fqctx);
            if (!fq_nmod_is_zero(v, fqctx))
            {
                changed = 1;
                fq_nmod_mul(u, v, inv_m_eval, fqctx);
                nmod_poly_mul(w, modulus, u);
                nmod_poly_add((Tcoeff + Toff)->coeffs + 0, (Fcoeff + Fi)->coeffs + 0, w);
            }
            else
            {
                nmod_poly_set((Tcoeff + Toff)->coeffs + 0, (Fcoeff + Fi)->coeffs + 0);                
            }
            Texp[Toff] = Aexp;
            Fi++;
            do {
                Aexp--;
            } while (Aexp >= 0 && fq_nmod_is_zero(Acoeff + Aexp, fqctx));
        }
        else if (Fi < Flen && (Aexp < 0 || Fexp[Fi] > Aexp))
        {
            /* F term ok, A term missing */
            nmod_poly_rem(v, (Fcoeff + Fi)->coeffs + 0, fqctx->modulus);
            if (!fq_nmod_is_zero(v, fqctx))
            {
                changed = 1;
                fq_nmod_mul(u, v, inv_m_eval, fqctx);
                nmod_poly_mul(w, u, modulus);
                nmod_poly_sub((Tcoeff + Toff)->coeffs + 0, (Fcoeff + Fi)->coeffs + 0, w);
            }
            else
            {
                nmod_poly_set((Tcoeff + Toff)->coeffs + 0, (Fcoeff + Fi)->coeffs + 0);                
            }
            Texp[Toff] = Fexp[Fi];
            Fi++;
        }
        else if (Aexp >= 0 && (Fi >= Flen || Fexp[Fi] < Aexp))
        {
            /* F term missing, A term ok */
            changed = 1;
            fq_nmod_mul(u, Acoeff + Aexp, inv_m_eval, fqctx);
            nmod_poly_mul((Tcoeff + Toff)->coeffs + 0, modulus, u);
            Texp[Toff] = Aexp;
            do {
                Aexp--;
            } while (Aexp >= 0 && fq_nmod_is_zero(Acoeff + Aexp, fqctx));
        }
        else
        {
            FLINT_ASSERT(0);
        }

        lastdeg = FLINT_MAX(lastdeg, nmod_poly_degree((Tcoeff + Toff)->coeffs + 0));
        mpoly_monomial_zero((Tcoeff + Toff)->exps + N*0, N);
        FLINT_ASSERT(!nmod_poly_is_zero((Tcoeff + Toff)->coeffs + 0));
        (Tcoeff + Toff)->length = 1;
        Toff++;
    }
    T->length = Toff;

    if (changed)
    {
        nmod_mpolyun_swap(T, F);
    }

    fq_nmod_clear(u, fqctx);
    fq_nmod_clear(v, fqctx);
    nmod_poly_clear(w);

    fq_nmod_clear(inv_m_eval, fqctx);


    *lastdeg_ = lastdeg;
    return changed;
}




int nmod_mpolyun_gcd_brown_lgprime_bivar(nmod_mpolyun_t G,
  nmod_mpolyun_t Abar, nmod_mpolyun_t Bbar, nmod_mpolyun_t A, nmod_mpolyun_t B,
                                         const nmod_mpoly_ctx_t ctx)
{
    int success;
    slong bound;
    fq_nmod_t temp, gammaeval;
    fq_nmod_poly_t Aeval, Beval, Geval, Abareval, Bbareval;
    nmod_mpolyun_t T;
    slong deggamma, ldegG, ldegAbar, ldegBbar, ldegA, ldegB;
    nmod_poly_t cA, cB, cG, cAbar, cBbar, gamma;
    nmod_poly_t modulus;
    slong deg;
    fq_nmod_mpoly_ctx_t ectx;
#if WANT_ASSERT
    nmod_poly_t leadA, leadB;
#endif

#if WANT_ASSERT
    nmod_poly_init(leadA, ctx->ffinfo->mod.n);
    nmod_poly_init(leadB, ctx->ffinfo->mod.n);
    nmod_poly_set(leadA, nmod_mpolyun_leadcoeff_ref(A, ctx));
    nmod_poly_set(leadB, nmod_mpolyun_leadcoeff_ref(B, ctx));
#endif

    nmod_poly_init(cA, ctx->ffinfo->mod.n);
    nmod_poly_init(cB, ctx->ffinfo->mod.n);
    nmod_mpolyun_content_last(cA, A, ctx);
    nmod_mpolyun_content_last(cB, B, ctx);
    nmod_mpolyun_divexact_last(A, cA, ctx);
    nmod_mpolyun_divexact_last(B, cB, ctx);

    nmod_poly_init(cG, ctx->ffinfo->mod.n);
    nmod_poly_gcd_euclidean(cG, cA, cB);

    nmod_poly_init(cAbar, ctx->ffinfo->mod.n);
    nmod_poly_init(cBbar, ctx->ffinfo->mod.n);
    nmod_poly_div(cAbar, cA, cG);
    nmod_poly_div(cBbar, cB, cG);

    nmod_poly_init(gamma, ctx->ffinfo->mod.n);
    nmod_poly_gcd(gamma, nmod_mpolyun_leadcoeff_ref(A, ctx),
                         nmod_mpolyun_leadcoeff_ref(B, ctx));

    ldegA = nmod_mpolyun_lastdeg(A, ctx);
    ldegB = nmod_mpolyun_lastdeg(B, ctx);
    deggamma = nmod_poly_degree(gamma);
    bound = 1 + deggamma + FLINT_MAX(ldegA, ldegB);

    nmod_mpolyun_init(T, A->bits, ctx);

    nmod_poly_init(modulus, ctx->ffinfo->mod.n);
    nmod_poly_one(modulus);

    deg = WORD(20)/(FLINT_BIT_COUNT(ctx->ffinfo->mod.n));
    deg = FLINT_MAX(WORD(2), deg);
    fq_nmod_mpoly_ctx_init_deg(ectx, ctx->minfo->nvars, ORD_LEX,
                                                      ctx->ffinfo->mod.n, deg);

    fq_nmod_poly_init(Aeval, ectx->fqctx);
    fq_nmod_poly_init(Beval, ectx->fqctx);
    fq_nmod_poly_init(Geval, ectx->fqctx);
    fq_nmod_poly_init(Abareval, ectx->fqctx);
    fq_nmod_poly_init(Bbareval, ectx->fqctx);
    fq_nmod_init(gammaeval, ectx->fqctx);
    fq_nmod_init(temp, ectx->fqctx);

    /* initialization already picked a prime */
    goto have_prime;

choose_prime: /* prime will be irreducible element of Fp[v] */

    /* same TODO */
    deg++;
    if (deg > 10000)
    {
        /* ran out of primes */
        success = 0;
        goto cleanup;
    }

    fq_nmod_mpoly_ctx_change_modulus(ectx, deg);

have_prime:

    /* make sure reduction does not kill both lc(A) and lc(B) */
    nmod_poly_rem(gammaeval, gamma, ectx->fqctx->modulus);
    if (fq_nmod_is_zero(gammaeval, ectx->fqctx))
    {
        goto choose_prime;
    }

    /* reduction should kill neither A nor B */
    nmod_mpolyun_reduce_last_fq_nmod_poly(Aeval, ectx->fqctx, A, ctx);
    nmod_mpolyun_reduce_last_fq_nmod_poly(Beval, ectx->fqctx, B, ctx);
    FLINT_ASSERT(Aeval->length > 0);
    FLINT_ASSERT(Beval->length > 0);

    fq_nmod_poly_gcd(Geval, Aeval, Beval, ectx->fqctx);
    success = fq_nmod_poly_divides(Abareval, Aeval, Geval, ectx->fqctx);
    FLINT_ASSERT(success);
    success = fq_nmod_poly_divides(Bbareval, Beval, Geval, ectx->fqctx);
    FLINT_ASSERT(success);

    FLINT_ASSERT(Geval->length > 0);
    FLINT_ASSERT(Abareval->length > 0);
    FLINT_ASSERT(Bbareval->length > 0);

    if (fq_nmod_poly_degree(Geval, ectx->fqctx) == 0)
    {
        nmod_mpolyun_one(G, ctx);
        nmod_mpolyun_swap(Abar, A);
        nmod_mpolyun_swap(Bbar, B);
        goto successful_put_content;    
    }

    if (nmod_poly_degree(modulus) > 0)
    {
        FLINT_ASSERT(G->length > 0);
        if (fq_nmod_poly_degree(Geval, ectx->fqctx) > G->exps[0])
        {
            goto choose_prime;
        }
        else if (fq_nmod_poly_degree(Geval, ectx->fqctx) < G->exps[0])
        {
            nmod_poly_one(modulus);
        }
    }

    fq_nmod_poly_scalar_mul_fq_nmod(Geval, Geval, gammaeval, ectx->fqctx);

    if (nmod_poly_degree(modulus) > 0)
    {
        nmod_mpolyun_addinterp_fq_nmod_poly(&ldegG, G, T, modulus, ctx, Geval, ectx->fqctx);
        nmod_mpolyun_addinterp_fq_nmod_poly(&ldegAbar, Abar, T, modulus, ctx, Abareval, ectx->fqctx);
        nmod_mpolyun_addinterp_fq_nmod_poly(&ldegBbar, Bbar, T, modulus, ctx, Bbareval, ectx->fqctx);
    }
    else
    {
        nmod_mpolyun_startinterp_fq_nmod_poly(&ldegG, G, ctx, Geval, ectx->fqctx);
        nmod_mpolyun_startinterp_fq_nmod_poly(&ldegAbar, Abar, ctx, Abareval, ectx->fqctx);
        nmod_mpolyun_startinterp_fq_nmod_poly(&ldegBbar, Bbar, ctx, Bbareval, ectx->fqctx);
    }

    nmod_poly_mul(modulus, modulus, ectx->fqctx->modulus);

    if (nmod_poly_degree(modulus) < bound)
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

    nmod_poly_one(modulus);
    goto choose_prime;

successful:

    nmod_mpolyun_content_last(modulus, G, ctx);
    nmod_mpolyun_divexact_last(G, modulus, ctx);
    nmod_mpolyun_divexact_last(Abar, nmod_mpolyun_leadcoeff_ref(G, ctx), ctx);
    nmod_mpolyun_divexact_last(Bbar, nmod_mpolyun_leadcoeff_ref(G, ctx), ctx);

successful_put_content:

    nmod_mpolyun_mul_last(G, cG, ctx);
    nmod_mpolyun_mul_last(Abar, cAbar, ctx);
    nmod_mpolyun_mul_last(Bbar, cBbar, ctx);

    success = 1;

cleanup:

#if WANT_ASSERT
    if (success)
    {
        nmod_poly_mul(modulus, nmod_mpolyun_leadcoeff_ref(G, ctx),
                               nmod_mpolyun_leadcoeff_ref(Abar, ctx));
        FLINT_ASSERT(nmod_poly_equal(modulus, leadA));
        nmod_poly_mul(modulus, nmod_mpolyun_leadcoeff_ref(G, ctx),
                               nmod_mpolyun_leadcoeff_ref(Bbar, ctx));
        FLINT_ASSERT(nmod_poly_equal(modulus, leadB));
    }
    nmod_poly_clear(leadA);
    nmod_poly_clear(leadB);
#endif

    nmod_poly_clear(cA);
    nmod_poly_clear(cB);
    nmod_poly_clear(cG);
    nmod_poly_clear(cAbar);
    nmod_poly_clear(cBbar);

    nmod_poly_clear(gamma);

    nmod_mpolyun_clear(T, ctx);

    nmod_poly_clear(modulus);

    fq_nmod_poly_clear(Aeval, ectx->fqctx);
    fq_nmod_poly_clear(Beval, ectx->fqctx);
    fq_nmod_poly_clear(Geval, ectx->fqctx);
    fq_nmod_poly_clear(Abareval, ectx->fqctx);
    fq_nmod_poly_clear(Bbareval, ectx->fqctx);
    fq_nmod_clear(gammaeval, ectx->fqctx);
    fq_nmod_clear(temp, ectx->fqctx);

    fq_nmod_mpoly_ctx_clear(ectx);

    return success;
}



/**************************
    fq_nmod stuff
***************************/




int fq_nmod_mpolyn_is_canonical(const fq_nmod_mpolyn_t A, const fq_nmod_mpoly_ctx_t ctx)
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
        slong l = (A->coeffs + i)->length;

        if (l == 0)
        {
            printf("\n mpolyn coeff length zero!!!!!!!\n");
            return 0;
        }

        if (fq_nmod_is_zero((A->coeffs + i)->coeffs + l - 1, ctx->fqctx))
        {
            printf("\n mpolyn coeff top coeff zero!!!!!!!\n");
            return 0;
        }
    }

    return 1;
}

int fq_nmod_mpolyun_is_canonical(const fq_nmod_mpolyun_t A, const fq_nmod_mpoly_ctx_t ctx)
{
    slong i;

    if (A->length > A->alloc)
    {
        return 0;
    }

    for (i = 0; i < A->length; i++)
    {
        if (!fq_nmod_mpolyn_is_canonical(A->coeffs + i, ctx))
        {
            return 0;
        }

        if (i > 0 && A->exps[i - 1] <= A->exps[i])
        {
            return 0;
        }
    }

    return 1;
}

int fq_nmod_mpolyn_is_nonzero_fq_nmod(const fq_nmod_mpolyn_t A, const fq_nmod_mpoly_ctx_t ctx)
{
    slong N;

    if (A->length != WORD(1))
    {
        return 0;
    }

    if (fq_nmod_poly_degree(A->coeffs + 0, ctx->fqctx) != 0)
    {
        return 0;
    }

    N = mpoly_words_per_exp(A->bits, ctx->minfo);
    return mpoly_monomial_is_zero(A->exps + N*0, N);
}

int fq_nmod_mpolyun_is_nonzero_fq_nmod(const fq_nmod_mpolyun_t A, const fq_nmod_mpoly_ctx_t ctx)
{
    if (A->length != 1 || A->exps[0] != 0)
    {
        return 0;
    }

    return fq_nmod_mpolyn_is_nonzero_fq_nmod(A->coeffs + 0, ctx);
}

void fq_nmod_mpolyn_scalar_mul_fq_nmod(fq_nmod_mpolyn_t A, const fq_nmod_t c, const fq_nmod_mpoly_ctx_t ctx)
{
    slong i;
    for (i = 0; i < A->length; i++)
    {
        fq_nmod_poly_scalar_mul_fq_nmod(A->coeffs + i, A->coeffs + i, c, ctx->fqctx);
    }
}

void fq_nmod_mpolyun_scalar_mul_fq_nmod(fq_nmod_mpolyun_t A, const fq_nmod_t c, const fq_nmod_mpoly_ctx_t ctx)
{
    slong i;
    FLINT_ASSERT(!fq_nmod_is_zero(c, ctx->fqctx));
    for (i = 0; i < A->length; i++)
    {
        fq_nmod_mpolyn_scalar_mul_fq_nmod(A->coeffs + i, c, ctx);
    }
}

/*
    get the leading coefficient in x_0,...,x_var
    A is in Fq[x_0, ... x_(var-1)][x_var]
    return is in Fq
*/
fq_nmod_struct * fq_nmod_mpolyn_leadcoeff_last_ref(fq_nmod_mpolyn_t A, const fq_nmod_mpoly_ctx_t ctx)
{
    fq_nmod_poly_struct * leadpoly;
    FLINT_ASSERT(A->length > 0);
    leadpoly = A->coeffs + 0;
    FLINT_ASSERT(leadpoly->length > 0);
    return leadpoly->coeffs + leadpoly->length - 1;
}






/*
    E = A(v = alpha)
    A is in Fq[X][v]
    E is in Fq[X]
*/
void fq_nmod_mpolyun_eval_last_bivar(fq_nmod_poly_t E, const fq_nmod_mpolyun_t A,
                             const fq_nmod_t alpha, const fq_nmod_mpoly_ctx_t ctx)
{
    fq_nmod_t v;
    slong Ai, Alen;
    fq_nmod_mpolyn_struct * Acoeff;
    ulong * Aexp;

    fq_nmod_init(v, ctx->fqctx);

    Acoeff = A->coeffs;
    Aexp = A->exps;
    Alen = A->length;
    Ai = 0;
    fq_nmod_poly_zero(E, ctx->fqctx);
    for (Ai = 0; Ai < Alen; Ai++)
    {
        fq_nmod_poly_evaluate_fq_nmod(v, (Acoeff + Ai)->coeffs + 0, alpha, ctx->fqctx);
        fq_nmod_poly_set_coeff(E, Aexp[Ai], v, ctx->fqctx);
    }

    fq_nmod_clear(v, ctx->fqctx);
}

/*
    A = B
    A is in Fq[X][v]  (no v appears)
    B is in Fq[X]
*/
void fq_nmod_mpolyun_startinterp_bivar(fq_nmod_mpolyun_t A, const fq_nmod_poly_t B,
                                                   const fq_nmod_mpoly_ctx_t ctx)
{
    slong N = mpoly_words_per_exp_sp(A->bits, ctx->minfo);
    slong Bexp;
    slong Blen = fq_nmod_poly_length(B, ctx->fqctx);
    fq_nmod_struct * Bcoeff = B->coeffs;
    fq_nmod_mpolyn_struct * Acoeff;
    ulong * Aexp;
    slong Ai;

    fq_nmod_mpolyun_fit_length(A, Blen, ctx);
    Acoeff = A->coeffs;
    Aexp = A->exps;

    Ai = 0;
    for (Bexp = Blen - 1; Bexp >= 0; Bexp--)
    {
        if (!fq_nmod_is_zero(Bcoeff + Bexp, ctx->fqctx))
        {
            FLINT_ASSERT(Ai < A->alloc);

            fq_nmod_mpolyn_fit_length(Acoeff + Ai, 1, ctx);
            mpoly_monomial_zero((Acoeff + Ai)->exps + N*0, N);
            fq_nmod_poly_set_fq_nmod((Acoeff + Ai)->coeffs + 0, Bcoeff + Bexp, ctx->fqctx);
            Aexp[Ai] = Bexp;
            (Acoeff + Ai)->length = 1;
            Ai++;
        }
    }
    A->length = Ai;
}

/*
    F = F + modulus*(A - F(v = alpha))
    no assumptions about matching monomials
    F is in Fq[X][v]
    A is in Fq[X]
    it is expected that modulus(alpha) == 1
*/
int fq_nmod_mpolyun_addinterp_bivar(slong * lastdeg_,
                   fq_nmod_mpolyun_t F, fq_nmod_mpolyun_t T, const fq_nmod_poly_t A,
      const fq_nmod_poly_t modulus, const fq_nmod_t alpha,  const fq_nmod_mpoly_ctx_t ctx)
{
    int changed = 0;
    slong lastdeg = -WORD(1);
    slong N = mpoly_words_per_exp_sp(F->bits, ctx->minfo);
    fq_nmod_t u, v;
    slong Fi, Toff, Aexp;
    fq_nmod_struct * Acoeff = A->coeffs;
    slong Flen = F->length;
    fq_nmod_mpolyn_struct * Fcoeff = F->coeffs;
    ulong * Fexp = F->exps;
    fq_nmod_mpolyn_struct * Tcoeff;
    ulong * Texp;
    fq_nmod_poly_t tp;

/*
flint_printf("\nnmod_mpolyun_addinterp_bivar: (alpha = %wu)\n", alpha);
flint_printf("A: "); nmod_poly_print_pretty(A, "X"); printf("\n");
flint_printf("F: "); nmod_mpolyun_print_pretty(F, NULL, ctx); printf("\n");
*/

    Fi = 0;
    Aexp = fq_nmod_poly_degree(A, ctx->fqctx);

    fq_nmod_init(u, ctx->fqctx);
    fq_nmod_init(v, ctx->fqctx);
    fq_nmod_poly_init(tp, ctx->fqctx);

    fq_nmod_mpolyun_fit_length(T, Flen + Aexp + 1, ctx);
    Tcoeff = T->coeffs;
    Texp = T->exps;
    Toff = 0;

    while (Fi < Flen || Aexp >= 0)
    {
        FLINT_ASSERT(Toff < T->alloc);
/*
flint_printf("Fi: %wd\n",Fi);
flint_printf("Aexp: %wd\n",Aexp);
*/
        if (Fi < Flen)
        {
            FLINT_ASSERT(!fq_nmod_poly_is_zero((Fcoeff + Fi)->coeffs + 0, ctx->fqctx));
            FLINT_ASSERT(fq_nmod_poly_degree((Fcoeff + Fi)->coeffs + 0, ctx->fqctx) < fq_nmod_poly_degree(modulus, ctx->fqctx));
        }

        if (Aexp >= 0)
        {
            FLINT_ASSERT(!fq_nmod_is_zero(Acoeff + Aexp, ctx->fqctx));
        }

        fq_nmod_mpolyn_fit_length(Tcoeff + Toff, 1, ctx);

        if (Fi < Flen && Aexp >= 0 && Fexp[Fi] == Aexp)
        {
            /* F term ok, A term ok */
            fq_nmod_poly_evaluate_fq_nmod(u, (Fcoeff + Fi)->coeffs + 0, alpha, ctx->fqctx);
            fq_nmod_sub(v, Acoeff + Aexp, u, ctx->fqctx);
            if (!fq_nmod_is_zero(v, ctx->fqctx))
            {
                changed = 1;
                fq_nmod_poly_scalar_mul_fq_nmod(tp, modulus, v, ctx->fqctx);
                fq_nmod_poly_add((Tcoeff + Toff)->coeffs + 0, (Fcoeff + Fi)->coeffs + 0, tp, ctx->fqctx);
            }
            else
            {
                fq_nmod_poly_set((Tcoeff + Toff)->coeffs + 0, (Fcoeff + Fi)->coeffs + 0, ctx->fqctx);
            }
            Texp[Toff] = Aexp;
            Fi++;
            do {
                Aexp--;
            } while (Aexp >= 0 && fq_nmod_is_zero(Acoeff + Aexp, ctx->fqctx));
        }
        else if (Fi < Flen && (Aexp < 0 || Fexp[Fi] > Aexp))
        {
            /* F term ok, A term missing */
            fq_nmod_poly_evaluate_fq_nmod(v, (Fcoeff + Fi)->coeffs + 0, alpha, ctx->fqctx);
            if (!fq_nmod_is_zero(v, ctx->fqctx))
            {
                changed = 1;
                fq_nmod_poly_scalar_mul_fq_nmod(tp, modulus, v, ctx->fqctx);
                fq_nmod_poly_sub((Tcoeff + Toff)->coeffs + 0, (Fcoeff + Fi)->coeffs + 0, tp, ctx->fqctx);
            }
            else
            {
                fq_nmod_poly_set((Tcoeff + Toff)->coeffs + 0, (Fcoeff + Fi)->coeffs + 0, ctx->fqctx);
            }
            Texp[Toff] = Fexp[Fi];
            Fi++;
        }
        else if (Aexp >= 0 && (Fi >= Flen || Fexp[Fi] < Aexp))
        {
            /* F term missing, A term ok */
            changed = 1;
            fq_nmod_poly_scalar_mul_fq_nmod((Tcoeff + Toff)->coeffs + 0, modulus, Acoeff + Aexp, ctx->fqctx);
            Texp[Toff] = Aexp;
            do {
                Aexp--;
            } while (Aexp >= 0 && fq_nmod_is_zero(Acoeff + Aexp, ctx->fqctx));
        }
        else
        {
            FLINT_ASSERT(0);
        }

        lastdeg = FLINT_MAX(lastdeg, fq_nmod_poly_degree((Tcoeff + Toff)->coeffs + 0, ctx->fqctx));
        mpoly_monomial_zero((Tcoeff + Toff)->exps + N*0, N);
        FLINT_ASSERT(!fq_nmod_poly_is_zero((Tcoeff + Toff)->coeffs + 0, ctx->fqctx));
        (Tcoeff + Toff)->length = 1;
        Toff++;
    }
    T->length = Toff;

    fq_nmod_clear(u, ctx->fqctx);
    fq_nmod_clear(v, ctx->fqctx);
    fq_nmod_poly_clear(tp, ctx->fqctx);

    if (changed)
    {
        fq_nmod_mpolyun_swap(T, F);
    }

    *lastdeg_ = lastdeg;
    return changed;
}










int fq_nmod_mpolyun_gcd_brown_smprime_bivar(fq_nmod_mpolyun_t G,
          fq_nmod_mpolyun_t Abar, fq_nmod_mpolyun_t Bbar, fq_nmod_mpolyun_t A,
                            fq_nmod_mpolyun_t B, const fq_nmod_mpoly_ctx_t ctx)
{
    int success;
    slong bound;
    fq_nmod_t alpha, temp, gammaeval;
    fq_nmod_poly_t Aeval, Beval, Geval, Abareval, Bbareval;
    fq_nmod_mpolyun_t T;
    slong deggamma, ldegG, ldegAbar, ldegBbar, ldegA, ldegB;
    fq_nmod_poly_t cA, cB, cG, cAbar, cBbar, gamma, trem;
    fq_nmod_poly_t modulus, tempmod;
    mp_bitcnt_t bits = A->bits;
#if WANT_ASSERT
    fq_nmod_poly_t leadA, leadB;
#endif

#if WANT_ASSERT
    fq_nmod_poly_init(leadA, ctx->fqctx);
    fq_nmod_poly_init(leadB, ctx->fqctx);
    fq_nmod_poly_set(leadA, fq_nmod_mpolyun_leadcoeff_ref(A, ctx), ctx->fqctx);
    fq_nmod_poly_set(leadB, fq_nmod_mpolyun_leadcoeff_ref(B, ctx), ctx->fqctx);
#endif

    fq_nmod_poly_init(cA, ctx->fqctx);
    fq_nmod_poly_init(cB, ctx->fqctx);
    fq_nmod_mpolyun_content_last(cA, A, ctx);
    fq_nmod_mpolyun_content_last(cB, B, ctx);
    fq_nmod_mpolyun_divexact_last(A, cA, ctx);
    fq_nmod_mpolyun_divexact_last(B, cB, ctx);

    fq_nmod_poly_init(cG, ctx->fqctx);
    fq_nmod_poly_gcd(cG, cA, cB, ctx->fqctx);

    fq_nmod_poly_init(cAbar, ctx->fqctx);
    fq_nmod_poly_init(cBbar, ctx->fqctx);
    fq_nmod_poly_init(trem, ctx->fqctx);
    fq_nmod_poly_divrem(cAbar, trem, cA, cG, ctx->fqctx);
    FLINT_ASSERT(fq_nmod_poly_is_zero(trem, ctx->fqctx));
    fq_nmod_poly_divrem(cBbar, trem, cB, cG, ctx->fqctx);
    FLINT_ASSERT(fq_nmod_poly_is_zero(trem, ctx->fqctx));

    fq_nmod_poly_init(gamma, ctx->fqctx);
    fq_nmod_poly_gcd(gamma, fq_nmod_mpolyun_leadcoeff_ref(A, ctx),
                            fq_nmod_mpolyun_leadcoeff_ref(B, ctx), ctx->fqctx);

    ldegA = fq_nmod_mpolyun_lastdeg(A, ctx);
    ldegB = fq_nmod_mpolyun_lastdeg(B, ctx);
    deggamma = fq_nmod_poly_degree(gamma, ctx->fqctx);

    bound = 1 + deggamma + FLINT_MAX(ldegA, ldegB);

    fq_nmod_mpolyun_init(T, bits, ctx);

    fq_nmod_poly_init(modulus, ctx->fqctx);
    fq_nmod_poly_one(modulus, ctx->fqctx);
    fq_nmod_poly_init(tempmod, ctx->fqctx);
    fq_nmod_poly_gen(tempmod, ctx->fqctx);
    fq_nmod_poly_neg(tempmod, tempmod, ctx->fqctx);

    fq_nmod_poly_init(Aeval, ctx->fqctx);
    fq_nmod_poly_init(Beval, ctx->fqctx);
    fq_nmod_poly_init(Geval, ctx->fqctx);
    fq_nmod_poly_init(Abareval, ctx->fqctx);
    fq_nmod_poly_init(Bbareval, ctx->fqctx);

    fq_nmod_init(gammaeval, ctx->fqctx);
    fq_nmod_init(alpha, ctx->fqctx);
    fq_nmod_init(temp, ctx->fqctx);

    fq_nmod_set_ui(alpha, 0, ctx->fqctx);

choose_prime:   /* prime is v - alpha */

    if (fq_nmod_next(alpha, ctx->fqctx) == 0)
    {
        success = 0;
        goto cleanup;
    }

    /* make sure evaluation point does not kill both lc(A) and lc(B) */
    fq_nmod_poly_evaluate_fq_nmod(gammaeval, gamma, alpha, ctx->fqctx);
    if (fq_nmod_is_zero(gammaeval, ctx->fqctx))
    {
        goto choose_prime;
    }

    /* evaluation point should kill neither A nor B */
    fq_nmod_mpolyun_eval_last_bivar(Aeval, A, alpha, ctx);
    fq_nmod_mpolyun_eval_last_bivar(Beval, B, alpha, ctx);
    FLINT_ASSERT(Aeval->length > 0);
    FLINT_ASSERT(Beval->length > 0);

    fq_nmod_poly_gcd(Geval, Aeval, Beval, ctx->fqctx);
    fq_nmod_poly_divrem(Abareval, trem, Aeval, Geval, ctx->fqctx);
    FLINT_ASSERT(fq_nmod_poly_is_zero(trem, ctx->fqctx));
    fq_nmod_poly_divrem(Bbareval, trem, Beval, Geval, ctx->fqctx);
    FLINT_ASSERT(fq_nmod_poly_is_zero(trem, ctx->fqctx));

    FLINT_ASSERT(Geval->length > 0);
    FLINT_ASSERT(Abareval->length);
    FLINT_ASSERT(Bbareval->length > 0);

    if (fq_nmod_poly_degree(Geval, ctx->fqctx) == 0)
    {
        fq_nmod_mpolyun_one(G, ctx);
        fq_nmod_mpolyun_swap(Abar, A);
        fq_nmod_mpolyun_swap(Bbar, B);
        goto successful_put_content;
    }

    if (fq_nmod_poly_degree(modulus, ctx->fqctx) > 0)
    {
        FLINT_ASSERT(G->length > 0);
        if (fq_nmod_poly_degree(Geval, ctx->fqctx) > G->exps[0])
        {
            goto choose_prime;
        }
        else if (fq_nmod_poly_degree(Geval, ctx->fqctx) < G->exps[0])
        {
            fq_nmod_poly_one(modulus, ctx->fqctx);
        }
    }

    fq_nmod_poly_scalar_mul_fq_nmod(Geval, Geval, gammaeval, ctx->fqctx);

    if (fq_nmod_poly_degree(modulus, ctx->fqctx) > 0)
    {
        fq_nmod_poly_evaluate_fq_nmod(temp, modulus, alpha, ctx->fqctx);
        fq_nmod_inv(temp, temp, ctx->fqctx);
        fq_nmod_poly_scalar_mul_fq_nmod(modulus, modulus, temp, ctx->fqctx);
        fq_nmod_mpolyun_addinterp_bivar(&ldegG, G, T, Geval, modulus, alpha, ctx);
        fq_nmod_mpolyun_addinterp_bivar(&ldegAbar, Abar, T, Abareval, modulus, alpha, ctx);
        fq_nmod_mpolyun_addinterp_bivar(&ldegBbar, Bbar, T, Bbareval, modulus, alpha, ctx);
    }
    else
    {
        fq_nmod_mpolyun_startinterp_bivar(G, Geval, ctx);
        fq_nmod_mpolyun_startinterp_bivar(Abar, Abareval, ctx);
        fq_nmod_mpolyun_startinterp_bivar(Bbar, Bbareval, ctx);
        ldegG = 0;
        ldegAbar = 0;
        ldegBbar = 0;
    }
    fq_nmod_poly_set_coeff(tempmod, 0, alpha, ctx->fqctx);
    fq_nmod_poly_mul(modulus, modulus, tempmod, ctx->fqctx);

    if (fq_nmod_poly_degree(modulus, ctx->fqctx) < bound)
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

    fq_nmod_poly_one(modulus, ctx->fqctx);
    goto choose_prime;

successful:

    fq_nmod_mpolyun_content_last(modulus, G, ctx);
    fq_nmod_mpolyun_divexact_last(G, modulus, ctx);
    fq_nmod_mpolyun_divexact_last(Abar, fq_nmod_mpolyun_leadcoeff_ref(G, ctx), ctx);
    fq_nmod_mpolyun_divexact_last(Bbar, fq_nmod_mpolyun_leadcoeff_ref(G, ctx), ctx);

successful_put_content:

    fq_nmod_mpolyun_mul_last(G, cG, ctx);
    fq_nmod_mpolyun_mul_last(Abar, cAbar, ctx);
    fq_nmod_mpolyun_mul_last(Bbar, cBbar, ctx);

    success = 1;

cleanup:

#if WANT_ASSERT
    if (success)
    {
        fq_nmod_poly_mul(modulus, fq_nmod_mpolyun_leadcoeff_ref(G, ctx),
                                  fq_nmod_mpolyun_leadcoeff_ref(Abar, ctx), ctx->fqctx);
        FLINT_ASSERT(fq_nmod_poly_equal(modulus, leadA, ctx->fqctx));
        fq_nmod_poly_mul(modulus, fq_nmod_mpolyun_leadcoeff_ref(G, ctx),
                                  fq_nmod_mpolyun_leadcoeff_ref(Bbar, ctx), ctx->fqctx);
        FLINT_ASSERT(fq_nmod_poly_equal(modulus, leadB, ctx->fqctx));
    }
    fq_nmod_poly_clear(leadA, ctx->fqctx);
    fq_nmod_poly_clear(leadB, ctx->fqctx);
#endif

    fq_nmod_poly_clear(cA, ctx->fqctx);
    fq_nmod_poly_clear(cB, ctx->fqctx);
    fq_nmod_poly_clear(cG, ctx->fqctx);
    fq_nmod_poly_clear(cAbar, ctx->fqctx);
    fq_nmod_poly_clear(cBbar, ctx->fqctx);
    fq_nmod_poly_clear(trem, ctx->fqctx);
    fq_nmod_poly_clear(gamma, ctx->fqctx);

    fq_nmod_poly_clear(Aeval, ctx->fqctx);
    fq_nmod_poly_clear(Beval, ctx->fqctx);
    fq_nmod_poly_clear(Geval, ctx->fqctx);
    fq_nmod_poly_clear(Abareval, ctx->fqctx);
    fq_nmod_poly_clear(Bbareval, ctx->fqctx);

    fq_nmod_mpolyun_clear(T, ctx);

    fq_nmod_clear(gammaeval, ctx->fqctx);
    fq_nmod_clear(alpha, ctx->fqctx);
    fq_nmod_clear(temp, ctx->fqctx);

    fq_nmod_poly_clear(modulus, ctx->fqctx);
    fq_nmod_poly_clear(tempmod, ctx->fqctx);

    return success;
}


/*
    E = A(x_var = alpha)
    A is in Fq[x_0, ..., x_(var-2), x_(var-1)][x_var]
    E is in Fq[x_0, ..., x_(var-2)][x_(var-1)]
*/
void fq_nmod_mpolyn_eval_last_n(fq_nmod_mpolyn_t E, fq_nmod_mpolyn_t A, slong var,
                                fq_nmod_t alpha, const fq_nmod_mpoly_ctx_t ctx)
{
    slong N = mpoly_words_per_exp_sp(A->bits, ctx->minfo);
    slong offset, shift, k;
    ulong mask;
    fq_nmod_t v;
    fq_nmod_poly_struct * Acoeff = A->coeffs;
    ulong * Aexp = A->exps;
    slong Alen = A->length;
    slong Ai;
    fq_nmod_poly_struct * Ecoeff;
    ulong * Eexp;
    slong Ei;

    fq_nmod_init(v, ctx->fqctx);

    FLINT_ASSERT(var > 0);
    mpoly_gen_offset_shift_sp(&offset, &shift, var - 1, A->bits, ctx->minfo);
    mask = (-UWORD(1)) >> (FLINT_BITS - A->bits);

    Ecoeff = E->coeffs;
    Eexp = E->exps;
    Ei = 0;
    for (Ai = 0; Ai < Alen; Ai++)
    {
        fq_nmod_poly_evaluate_fq_nmod(v, Acoeff + Ai, alpha, ctx->fqctx);
        k = ((Aexp + N*Ai)[offset] >> shift) & mask;
        if (fq_nmod_is_zero(v, ctx->fqctx))
        {
            continue;
        }

        if (Ei > 0 && mpoly_monomial_equal_extra(Eexp + N*(Ei - 1), Aexp + N*Ai, N, offset, -(k << shift)))
        {
            /* append to previous */
            fq_nmod_poly_set_coeff(Ecoeff + Ei - 1, k, v, ctx->fqctx);
        }
        else
        {
            FLINT_ASSERT(Ei == 0 || mpoly_monomial_gt_nomask_extra(Eexp + N*(Ei - 1), Aexp + N*Ai, N, offset, -(k << shift)));

            /* create new */
            if (Ei >= E->alloc)
            {
                fq_nmod_mpolyn_fit_length(E, Ei + 1, ctx);
                Ecoeff = E->coeffs;
                Eexp = E->exps;
            }
            mpoly_monomial_set_extra(Eexp + N*Ei, Aexp + N*Ai, N, offset, -(k << shift));
            fq_nmod_poly_zero(Ecoeff + Ei, ctx->fqctx);
            fq_nmod_poly_set_coeff(Ecoeff + Ei, k, v, ctx->fqctx);
            Ei++;
        }
    }
    E->length = Ei;

    fq_nmod_clear(v, ctx->fqctx);
}


/*
    E = A(x_var = alpha)
    A is in Fq[X][x_0, ..., x_(var-2), x_(var-1)][x_var]
    E is in Fq[X][x_0, ..., x_(var-2)][x_(var-1)]
*/
void fq_nmod_mpolyun_eval_last_un(fq_nmod_mpolyun_t E, fq_nmod_mpolyun_t A, slong var,
                                   fq_nmod_t alpha, const fq_nmod_mpoly_ctx_t ctx)
{
    fq_nmod_mpolyn_struct * Acoeff = A->coeffs;
    ulong * Aexp = A->exps;
    slong Alen = A->length;
    slong Ai;
    fq_nmod_mpolyn_struct * Ecoeff;
    ulong * Eexp;
    slong Ei;

    FLINT_ASSERT(var > 0);

    fq_nmod_mpolyun_fit_length(E, Alen, ctx);
    Ecoeff = E->coeffs;
    Eexp = E->exps;

    Ei = 0;
    for (Ai = 0; Ai < Alen; Ai++)
    {
        fq_nmod_mpolyn_eval_last_n(Ecoeff + Ei, Acoeff + Ai, var, alpha, ctx);
        Eexp[Ei] = Aexp[Ai];
        Ei += ((Ecoeff + Ei)->length != 0);
    }
    E->length = Ei;

    FLINT_ASSERT(fq_nmod_mpolyun_is_canonical(E, ctx));
}

/*
    T = F + modulus*(A - F(x_var = alpha))
    no assumptions about matching monomials
    F is in Fp[x_0, ..., x_(var-1), x_(var-1)][x_var]
    A is in Fp[x_0, ..., x_(var-2)][x_(var-1)]
    in order to fxn correctly, modulus(alpha) should be 1
*/
int fq_nmod_mpolyn_addinterp_n(slong * lastdeg_,
             fq_nmod_mpolyn_t T, fq_nmod_mpolyn_t F, fq_nmod_mpolyn_t A, slong var,
              fq_nmod_poly_t modulus, const fq_nmod_t alpha, const fq_nmod_mpoly_ctx_t ctx)
{
    int changed = 0;
    slong N = mpoly_words_per_exp_sp(A->bits, ctx->minfo);
    slong lastdeg = -WORD(1);
    slong offset, shift;
    slong vi;
    fq_nmod_t v;
    fq_nmod_poly_t tp;

    fq_nmod_poly_struct * Tcoeff;
    ulong * Texp;
    slong Ti;

    fq_nmod_poly_struct * Acoeff = A->coeffs;
    slong Alen = A->length;
    ulong * Aexp = A->exps;
    slong Ai;

    fq_nmod_poly_struct * Fcoeff = F->coeffs;
    slong Flen = F->length;
    ulong * Fexp = F->exps;
    slong Fi;

    fq_nmod_poly_init(tp, ctx->fqctx);
    fq_nmod_init(v, ctx->fqctx);

    FLINT_ASSERT(var > 0);
    FLINT_ASSERT(T->bits == A->bits);
    FLINT_ASSERT(F->bits == A->bits);

    mpoly_gen_offset_shift_sp(&offset, &shift, var - 1, A->bits, ctx->minfo);

    Flen = F->length;

    fq_nmod_mpolyn_fit_length(T, FLINT_MAX(Flen, Alen), ctx);
    Tcoeff = T->coeffs;
    Texp = T->exps;
    Ti = 0;

    Fi = Ai = vi = 0;
    if (Ai < Alen)
    {
        vi = fq_nmod_poly_degree(A->coeffs + Ai, ctx->fqctx);
    }
    while (Fi < Flen || Ai < Alen)
    {
        if (Ti >= T->alloc)
        {
            fq_nmod_mpolyn_fit_length(T, Ti + FLINT_MAX(Flen - Fi, Alen - Ai), ctx);
            Tcoeff = T->coeffs;
            Texp = T->exps;
        }

        if (Ai < Alen)
        {
            FLINT_ASSERT(!fq_nmod_is_zero((Acoeff + Ai)->coeffs + vi, ctx->fqctx));
        }

        if (Fi < Flen && Ai < Alen && mpoly_monomial_equal_extra(Fexp + N*Fi, Aexp + N*Ai, N, offset, vi << shift))
        {
            /* F term ok, A term ok */
            fq_nmod_poly_evaluate_fq_nmod(v, Fcoeff + Fi, alpha, ctx->fqctx);
            fq_nmod_sub(v, (Acoeff + Ai)->coeffs + vi, v, ctx->fqctx);
            if (!fq_nmod_is_zero(v, ctx->fqctx))
            {
                changed = 1;
                fq_nmod_poly_scalar_mul_fq_nmod(tp, modulus, v, ctx->fqctx);
                fq_nmod_poly_add(Tcoeff + Ti, Fcoeff + Fi, tp, ctx->fqctx);
            }
            else
            {
                fq_nmod_poly_set(Tcoeff + Ti, Fcoeff + Fi, ctx->fqctx);
            }
            mpoly_monomial_set(Texp + N*Ti, Fexp + N*Fi, N);

            Fi++;
            do {
                vi--;
            } while (vi >= 0 && fq_nmod_is_zero((Acoeff + Ai)->coeffs + vi, ctx->fqctx));
            if (vi < 0)
            {
                Ai++;
                if (Ai < Alen)
                {
                    vi = fq_nmod_poly_degree(A->coeffs + Ai, ctx->fqctx);
                }
            }
        }
        else if (Fi < Flen && (Ai >= Alen || mpoly_monomial_gt_nomask_extra(Fexp + N*Fi, Aexp + N*Ai, N, offset, vi << shift)))
        {
            /* F term ok, A term missing */
            fq_nmod_poly_evaluate_fq_nmod(v, Fcoeff + Fi, alpha, ctx->fqctx);
            if (!fq_nmod_is_zero(v, ctx->fqctx))
            {
                changed = 1;
                fq_nmod_poly_scalar_mul_fq_nmod(tp, modulus, v, ctx->fqctx);
                fq_nmod_poly_sub(Tcoeff + Ti, Fcoeff + Fi, tp, ctx->fqctx);
            }
            else
            {
                fq_nmod_poly_set(Tcoeff + Ti, Fcoeff + Fi, ctx->fqctx);
            }
            mpoly_monomial_set(Texp + N*Ti, Fexp + N*Fi, N);

            Fi++;
        }
        else
        {
            FLINT_ASSERT(Ai < Alen && (Fi >= Flen || mpoly_monomial_lt_nomask_extra(Fexp + N*Fi, Aexp + N*Ai, N, offset, vi << shift)));

            /* F term missing, A term ok */
            changed = 1;
            fq_nmod_poly_scalar_mul_fq_nmod(Tcoeff + Ti, modulus, (Acoeff + Ai)->coeffs + vi, ctx->fqctx);
            mpoly_monomial_set_extra(Texp + N*Ti, Aexp + N*Ai, N, offset, vi << shift);

            do {
                vi--;
            } while (vi >= 0 && fq_nmod_is_zero((Acoeff + Ai)->coeffs + vi, ctx->fqctx));
            if (vi < 0)
            {
                Ai++;
                if (Ai < Alen)
                {
                    vi = fq_nmod_poly_degree(A->coeffs + Ai, ctx->fqctx);
                }
            }
        }

        FLINT_ASSERT(!fq_nmod_poly_is_zero(Tcoeff + Ti, ctx->fqctx));
        lastdeg = FLINT_MAX(lastdeg, fq_nmod_poly_degree(Tcoeff + Ti, ctx->fqctx));
        Ti++;
    }
    T->length = Ti;

    fq_nmod_poly_clear(tp, ctx->fqctx);
    fq_nmod_clear(v, ctx->fqctx);

    *lastdeg_ = FLINT_MAX(*lastdeg_, lastdeg);
    return changed;
}

/*
    F = F + modulus*(A - F(x_var = alpha))
    no assumptions about matching monomials
    F is in Fp[X][x_0, ..., x_(var-2), x_(var-1)][x_var]
    A is in Fp[X][x_0, ..., x_(var-2)][x_(var-1)]
    in order to fxn correctly, modulus(alpha) should be 1
*/
int fq_nmod_mpolyun_addinterp_un(slong * lastdeg,
             fq_nmod_mpolyun_t F, fq_nmod_mpolyun_t T, fq_nmod_mpolyun_t A, slong var,
              fq_nmod_poly_t modulus, const fq_nmod_t alpha, const fq_nmod_mpoly_ctx_t ctx)
{
    int changed = 0;
    slong i, j, k;
    ulong * Texp;
    ulong * Fexp;
    ulong * Aexp;
    slong Flen;
    slong Alen;
    fq_nmod_mpolyn_struct * Tcoeff;
    fq_nmod_mpolyn_struct * Fcoeff;
    fq_nmod_mpolyn_struct  * Acoeff;
    fq_nmod_mpolyn_t zero;

    FLINT_ASSERT(var > 0);

    *lastdeg = -WORD(1);

    FLINT_ASSERT(F->bits == T->bits);
    FLINT_ASSERT(T->bits == A->bits);

    Flen = F->length;
    Alen = A->length;
    fq_nmod_mpolyun_fit_length(T, Flen + Alen, ctx);

    Tcoeff = T->coeffs;
    Fcoeff = F->coeffs;
    Acoeff = A->coeffs;
    Texp = T->exps;
    Fexp = F->exps;
    Aexp = A->exps;   

    fq_nmod_mpolyn_init(zero, A->bits, ctx);
    zero->length = 0;

    i = j = k = 0;
    while (i < Flen || j < Alen)
    {
        if (i < Flen && j < Alen && (Fexp[i] == Aexp[j]))
        {
            /* F term ok, A term ok */
            changed |= fq_nmod_mpolyn_addinterp_n(lastdeg, Tcoeff + k, Fcoeff + i,
                                         Acoeff + j, var, modulus, alpha, ctx);
            Texp[k] = Aexp[j];
            i++;
            j++;
        }
        else if (i < Flen && (j >= Alen || Fexp[i] > Aexp[j]))
        {
            /* F term ok, A term missing */
            changed |= fq_nmod_mpolyn_addinterp_n(lastdeg, Tcoeff + k, Fcoeff + i,
                                               zero, var, modulus, alpha, ctx);
            Texp[k] = Fexp[i];
            i++;
        }
        else
        {
            FLINT_ASSERT(j < Alen && (i >= Flen || Aexp[j] > Fexp[i]));

            /* F term missing, A term ok */
            changed |= fq_nmod_mpolyn_addinterp_n(lastdeg, Tcoeff + k, zero,
                                         Acoeff + j, var, modulus, alpha, ctx);
            Texp[k] = Aexp[j];
            j++;
        }

        FLINT_ASSERT(!fq_nmod_mpolyn_is_zero(Tcoeff + k, ctx));
        k++;
    }
    T->length = k;

    if (changed)
    {
        fq_nmod_mpolyun_swap(T, F);
    }

    fq_nmod_mpolyn_clear(zero, ctx);

    FLINT_ASSERT(fq_nmod_mpolyun_is_canonical(F, ctx));

    return changed;    
}


/*
    A = B
    A is in Fq[x_0, ..., x_(var-1), x_(var-1)][x_var]
    B is in Fq[x_0, ..., x_(var-2)][x_(var-1)]
*/
void fq_nmod_mpolyn_startinterp_n(fq_nmod_mpolyn_t A, fq_nmod_mpolyn_t B, slong var,
                                                    const fq_nmod_mpoly_ctx_t ctx)
{
    slong N = mpoly_words_per_exp_sp(B->bits, ctx->minfo);
    slong offset, shift;
    slong vi;

    fq_nmod_poly_struct * Bcoeff = B->coeffs;
    ulong * Bexp = B->exps;
    slong Blen = B->length;
    slong Bi;

    fq_nmod_poly_struct * Acoeff;
    ulong * Aexp;
    slong Ai;

    fq_nmod_mpolyn_fit_length(A, Blen, ctx);
    Acoeff = A->coeffs;
    Aexp = A->exps;

    mpoly_gen_offset_shift_sp(&offset, &shift, var - 1, A->bits, ctx->minfo);

    Ai = 0;
    for (Bi = 0; Bi < Blen; Bi++)
    {
        if (Ai + (Bcoeff + Bi)->length >= A->alloc)
        {
            fq_nmod_mpolyn_fit_length(A, Ai + (Bcoeff + Bi)->length, ctx);
            Acoeff = A->coeffs;
            Aexp = A->exps;
        }
        for (vi = (Bcoeff + Bi)->length - 1; vi >= 0; vi--)
        {
            if (!fq_nmod_is_zero((Bcoeff + Bi)->coeffs + vi, ctx->fqctx))
            {
                mpoly_monomial_set_extra(Aexp + N*Ai, Bexp + N*Bi, N, offset, vi << shift);
                fq_nmod_poly_zero(Acoeff + Ai, ctx->fqctx);
                fq_nmod_poly_set_coeff(Acoeff + Ai, 0, (Bcoeff + Bi)->coeffs + vi, ctx->fqctx);
                Ai++;
            }
        }
    }
    A->length = Ai;
}


/*
    A = B
    A is in Fq[X][x_0, ..., x_(var-2), x_(var-1)][x_var]
    B is in Fq[X][x_0, ..., x_(var-2)][x_(var-1)]
*/
void fq_nmod_mpolyun_startinterp_un(fq_nmod_mpolyun_t A, fq_nmod_mpolyun_t B, slong var,
                                                    const fq_nmod_mpoly_ctx_t ctx)
{
    slong i;
    fq_nmod_mpolyn_struct * Bcoeff = B->coeffs;
    ulong * Bexp = B->exps;
    slong Blen = B->length;

    fq_nmod_mpolyn_struct * Acoeff;
    ulong * Aexp;

    fq_nmod_mpolyun_fit_length(A, Blen, ctx);
    Acoeff = A->coeffs;
    Aexp = A->exps;

    for (i = 0; i < Blen; i++)
    {
        Aexp[i] = Bexp[i];
        fq_nmod_mpolyn_startinterp_n(Acoeff + i, Bcoeff + i, var, ctx);
    }
    A->length = Blen;

    FLINT_ASSERT(fq_nmod_mpolyun_is_canonical(A, ctx));
}




int fq_nmod_mpolyun_gcd_brown_smprime(fq_nmod_mpolyun_t G,
          fq_nmod_mpolyun_t Abar, fq_nmod_mpolyun_t Bbar, fq_nmod_mpolyun_t A,
                  fq_nmod_mpolyun_t B, slong var, const fq_nmod_mpoly_ctx_t ctx)
{
    int success;
    slong bound;
    slong offset, shift;
    fq_nmod_t alpha, temp, gammaeval;
    fq_nmod_mpolyun_t Aeval, Beval, Geval, Abareval, Bbareval;
    fq_nmod_mpolyun_t T;
    slong deggamma, ldegG, ldegAbar, ldegBbar, ldegA, ldegB;
    fq_nmod_poly_t cA, cB, cG, cAbar, cBbar, gamma, trem;
    fq_nmod_poly_t modulus, tempmod;
    mp_bitcnt_t bits = A->bits;
    slong N = mpoly_words_per_exp(bits, ctx->minfo);
#if WANT_ASSERT
    fq_nmod_poly_t leadA, leadB;
#endif

    FLINT_ASSERT(var >= 0);
    if (var == WORD(0))
    {
        /* bivariate is more comfortable separated */
        return fq_nmod_mpolyun_gcd_brown_smprime_bivar(G, Abar, Bbar, A, B, ctx);
    }

    mpoly_gen_offset_shift_sp(&offset, &shift, var - 1, G->bits, ctx->minfo);

#if WANT_ASSERT
    fq_nmod_poly_init(leadA, ctx->fqctx);
    fq_nmod_poly_init(leadB, ctx->fqctx);
    fq_nmod_poly_set(leadA, fq_nmod_mpolyun_leadcoeff_ref(A, ctx), ctx->fqctx);
    fq_nmod_poly_set(leadB, fq_nmod_mpolyun_leadcoeff_ref(B, ctx), ctx->fqctx);
#endif

    fq_nmod_poly_init(cA, ctx->fqctx);
    fq_nmod_poly_init(cB, ctx->fqctx);
    fq_nmod_mpolyun_content_last(cA, A, ctx);
    fq_nmod_mpolyun_content_last(cB, B, ctx);
    fq_nmod_mpolyun_divexact_last(A, cA, ctx);
    fq_nmod_mpolyun_divexact_last(B, cB, ctx);

    fq_nmod_poly_init(cG, ctx->fqctx);
    fq_nmod_poly_gcd(cG, cA, cB, ctx->fqctx);

    fq_nmod_poly_init(cAbar, ctx->fqctx);
    fq_nmod_poly_init(cBbar, ctx->fqctx);
    fq_nmod_poly_init(trem, ctx->fqctx);
    fq_nmod_poly_divrem(cAbar, trem, cA, cG, ctx->fqctx);
    FLINT_ASSERT(fq_nmod_poly_is_zero(trem, ctx->fqctx));
    fq_nmod_poly_divrem(cBbar, trem, cB, cG, ctx->fqctx);
    FLINT_ASSERT(fq_nmod_poly_is_zero(trem, ctx->fqctx));

    fq_nmod_poly_init(gamma, ctx->fqctx);
    fq_nmod_poly_gcd(gamma, fq_nmod_mpolyun_leadcoeff_ref(A, ctx),
                            fq_nmod_mpolyun_leadcoeff_ref(B, ctx), ctx->fqctx);

    ldegA = fq_nmod_mpolyun_lastdeg(A, ctx);
    ldegB = fq_nmod_mpolyun_lastdeg(B, ctx);
    deggamma = fq_nmod_poly_degree(gamma, ctx->fqctx);

    bound = 1 + deggamma + FLINT_MAX(ldegA, ldegB);

    fq_nmod_mpolyun_init(T, bits, ctx);
    fq_nmod_poly_init(modulus, ctx->fqctx);
    fq_nmod_poly_one(modulus, ctx->fqctx);
    fq_nmod_poly_init(tempmod, ctx->fqctx);
    fq_nmod_poly_gen(tempmod, ctx->fqctx);
    fq_nmod_poly_neg(tempmod, tempmod, ctx->fqctx);

    fq_nmod_mpolyun_init(Aeval, bits, ctx);
    fq_nmod_mpolyun_init(Beval, bits, ctx);
    fq_nmod_mpolyun_init(Geval, bits, ctx);
    fq_nmod_mpolyun_init(Abareval, bits, ctx);
    fq_nmod_mpolyun_init(Bbareval, bits, ctx);

    fq_nmod_init(gammaeval, ctx->fqctx);
    fq_nmod_init(alpha, ctx->fqctx);
    fq_nmod_init(temp, ctx->fqctx);

    fq_nmod_set_ui(alpha, 0, ctx->fqctx);

choose_prime:

    if (fq_nmod_next(alpha, ctx->fqctx) == 0)
    {
        success = 0;
        goto cleanup;
    }

    /* make sure evaluation point does not kill both lc(A) and lc(B) */
    fq_nmod_poly_evaluate_fq_nmod(gammaeval, gamma, alpha, ctx->fqctx);
    if (fq_nmod_is_zero(gammaeval, ctx->fqctx))
    {
        goto choose_prime;
    }

    /* evaluation point should kill neither A nor B */
    fq_nmod_mpolyun_eval_last_un(Aeval, A, var, alpha, ctx);
    fq_nmod_mpolyun_eval_last_un(Beval, B, var, alpha, ctx);
    FLINT_ASSERT(Aeval->length > 0);
    FLINT_ASSERT(Beval->length > 0);

    success = fq_nmod_mpolyun_gcd_brown_smprime(Geval, Abareval, Bbareval,
                                                   Aeval, Beval, var - 1, ctx);
    if (success == 0)
    {
        goto choose_prime;
    }

    FLINT_ASSERT(Geval->length > 0);
    FLINT_ASSERT(Abareval->length > 0);
    FLINT_ASSERT(Bbareval->length > 0);

    if (fq_nmod_mpolyun_is_nonzero_fq_nmod(Geval, ctx))
    {
        fq_nmod_mpolyun_one(G, ctx);
        fq_nmod_mpolyun_swap(Abar, A);
        fq_nmod_mpolyun_swap(Bbar, B);
        goto successful_put_content;    
    }

    if (fq_nmod_poly_degree(modulus, ctx->fqctx) > 0)
    {
        int cmp = 0;
        FLINT_ASSERT(G->length > 0);
        if (G->exps[0] != Geval->exps[0])
        {
            cmp = G->exps[0] > Geval->exps[0] ? 1 : -1;
        }
        if (cmp == 0)
        {
            slong k = fq_nmod_poly_degree((Geval->coeffs + 0)->coeffs + 0, ctx->fqctx);
            FLINT_ASSERT(k >= 0);
            cmp = mpoly_monomial_cmp_nomask_extra(
                        (G->coeffs + 0)->exps + N*0,
                    (Geval->coeffs + 0)->exps + N*0, N, offset, k << shift);
        }

        if (cmp < 0)
        {
            goto choose_prime;
        }
        else if (cmp > 0)
        {
            fq_nmod_poly_one(modulus, ctx->fqctx);
        }
    }

    fq_nmod_inv(temp, fq_nmod_mpolyn_leadcoeff_last_ref(Geval->coeffs + 0, ctx), ctx->fqctx);
    fq_nmod_mul(temp, temp, gammaeval, ctx->fqctx);
    fq_nmod_mpolyun_scalar_mul_fq_nmod(Geval, temp, ctx);

    if (fq_nmod_poly_degree(modulus, ctx->fqctx) > 0)
    {
        fq_nmod_poly_evaluate_fq_nmod(temp, modulus, alpha, ctx->fqctx);
        fq_nmod_inv(temp, temp, ctx->fqctx);
        fq_nmod_poly_scalar_mul_fq_nmod(modulus, modulus, temp, ctx->fqctx);
        fq_nmod_mpolyun_addinterp_un(&ldegG, G, T, Geval, var, modulus, alpha, ctx);
        fq_nmod_mpolyun_addinterp_un(&ldegAbar, Abar, T, Abareval, var, modulus, alpha, ctx);
        fq_nmod_mpolyun_addinterp_un(&ldegBbar, Bbar, T, Bbareval, var, modulus, alpha, ctx);
    }
    else
    {
        fq_nmod_mpolyun_startinterp_un(G, Geval, var, ctx);
        fq_nmod_mpolyun_startinterp_un(Abar, Abareval, var, ctx);
        fq_nmod_mpolyun_startinterp_un(Bbar, Bbareval, var, ctx);
        ldegG = 0;
        ldegAbar = 0;
        ldegBbar = 0;
    }
    fq_nmod_poly_set_coeff(tempmod, 0, alpha, ctx->fqctx);
    fq_nmod_poly_mul(modulus, modulus, tempmod, ctx->fqctx);

    if (fq_nmod_poly_degree(modulus, ctx->fqctx) < bound)
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

    fq_nmod_poly_one(modulus, ctx->fqctx);
    goto choose_prime;

successful:

    fq_nmod_mpolyun_content_last(modulus, G, ctx);
    fq_nmod_mpolyun_divexact_last(G, modulus, ctx);
    fq_nmod_mpolyun_divexact_last(Abar, fq_nmod_mpolyun_leadcoeff_ref(G, ctx), ctx);
    fq_nmod_mpolyun_divexact_last(Bbar, fq_nmod_mpolyun_leadcoeff_ref(G, ctx), ctx);

successful_put_content:

    fq_nmod_mpolyun_mul_last(G, cG, ctx);
    fq_nmod_mpolyun_mul_last(Abar, cAbar, ctx);
    fq_nmod_mpolyun_mul_last(Bbar, cBbar, ctx);

    success = 1;

cleanup:

#if WANT_ASSERT
    if (success)
    {
        fq_nmod_poly_mul(modulus, fq_nmod_mpolyun_leadcoeff_ref(G, ctx),
                                  fq_nmod_mpolyun_leadcoeff_ref(Abar, ctx), ctx->fqctx);
        FLINT_ASSERT(fq_nmod_poly_equal(modulus, leadA, ctx->fqctx));
        fq_nmod_poly_mul(modulus, fq_nmod_mpolyun_leadcoeff_ref(G, ctx),
                                  fq_nmod_mpolyun_leadcoeff_ref(Bbar, ctx), ctx->fqctx);
        FLINT_ASSERT(fq_nmod_poly_equal(modulus, leadB, ctx->fqctx));
    }
    fq_nmod_poly_clear(leadA, ctx->fqctx);
    fq_nmod_poly_clear(leadB, ctx->fqctx);
#endif

    fq_nmod_poly_clear(cA, ctx->fqctx);
    fq_nmod_poly_clear(cB, ctx->fqctx);
    fq_nmod_poly_clear(cG, ctx->fqctx);
    fq_nmod_poly_clear(cAbar, ctx->fqctx);
    fq_nmod_poly_clear(cBbar, ctx->fqctx);
    fq_nmod_poly_clear(trem, ctx->fqctx);
    fq_nmod_poly_clear(gamma, ctx->fqctx);

    fq_nmod_mpolyun_clear(Aeval, ctx);
    fq_nmod_mpolyun_clear(Beval, ctx);
    fq_nmod_mpolyun_clear(Geval, ctx);
    fq_nmod_mpolyun_clear(Abareval, ctx);
    fq_nmod_mpolyun_clear(Bbareval, ctx);

    fq_nmod_mpolyun_clear(T, ctx);

    fq_nmod_clear(gammaeval, ctx->fqctx);
    fq_nmod_clear(alpha, ctx->fqctx);
    fq_nmod_clear(temp, ctx->fqctx);

    fq_nmod_poly_clear(modulus, ctx->fqctx);
    fq_nmod_poly_clear(tempmod, ctx->fqctx);

    return success;
}


/*
    E = A mod alpha(x_var)
    A is in                Fp[x_0, ..., x_(var-2), x_(var-1)][x_var]
    E is in (Fp/alpha(x_var))[x_0, ..., x_(var-2)][x_(var-1)]
    alpha is ectx->modulus
*/

void nmod_mpolyn_redto_fq_nmod_mpolyn(fq_nmod_mpolyn_t E, fq_nmod_mpoly_ctx_t ectx,
                        nmod_mpolyn_t A, slong var, const nmod_mpoly_ctx_t ctx)
{
    slong N = mpoly_words_per_exp_sp(A->bits, ctx->minfo);
    slong offset, shift, k;
    ulong mask;
    fq_nmod_t v;
    nmod_poly_struct * Acoeff = A->coeffs;
    ulong * Aexp = A->exps;
    slong Alen = A->length;
    slong Ai;
    fq_nmod_poly_struct * Ecoeff;
    ulong * Eexp;
    slong Ei;

    fq_nmod_init(v, ectx->fqctx);

    FLINT_ASSERT(var > 0);
    mpoly_gen_offset_shift_sp(&offset, &shift, var - 1, A->bits, ctx->minfo);
    mask = (-UWORD(1)) >> (FLINT_BITS - A->bits);

    Ecoeff = E->coeffs;
    Eexp = E->exps;
    Ei = 0;
    for (Ai = 0; Ai < Alen; Ai++)
    {
        nmod_poly_rem(v, Acoeff + Ai, ectx->fqctx->modulus);
        k = ((Aexp + N*Ai)[offset] >> shift) & mask;
        if (fq_nmod_is_zero(v, ectx->fqctx))
        {
            continue;
        }

        if (Ei > 0 && mpoly_monomial_equal_extra(Eexp + N*(Ei - 1), Aexp + N*Ai, N, offset, -(k << shift)))
        {
            /* append to previous */
            fq_nmod_poly_set_coeff(Ecoeff + Ei - 1, k, v, ectx->fqctx);
        }
        else
        {
            FLINT_ASSERT(Ei == 0 || mpoly_monomial_gt_nomask_extra(Eexp + N*(Ei - 1), Aexp + N*Ai, N, offset, -(k << shift)));

            /* create new */
            if (Ei >= E->alloc)
            {
                fq_nmod_mpolyn_fit_length(E, Ei + 1, ectx);
                Ecoeff = E->coeffs;
                Eexp = E->exps;
            }
            mpoly_monomial_set_extra(Eexp + N*Ei, Aexp + N*Ai, N, offset, -(k << shift));
            fq_nmod_poly_zero(Ecoeff + Ei, ectx->fqctx);
            fq_nmod_poly_set_coeff(Ecoeff + Ei, k, v, ectx->fqctx);
            Ei++;
        }
    }
    E->length = Ei;

    fq_nmod_clear(v, ectx->fqctx);
}

/*
    E = A mod alpha(x_var)
    A is in               Fp [X][x_0, ..., x_(var-2), x_(var-1)][x_var]
    E is in (Fp/alpha(x_var))[X][x_0, ..., x_(var-2)][x_(var-1)]
    alpha is ectx->modulus
*/
void nmod_mpolyun_redto_fq_nmod_mpolyun(fq_nmod_mpolyun_t E, fq_nmod_mpoly_ctx_t ectx,
                               nmod_mpolyun_t A, slong var, const nmod_mpoly_ctx_t ctx)
{
    nmod_mpolyn_struct * Acoeff = A->coeffs;
    ulong * Aexp = A->exps;
    slong Alen = A->length;
    slong Ai;
    fq_nmod_mpolyn_struct * Ecoeff;
    ulong * Eexp;
    slong Ei;

    fq_nmod_mpolyun_fit_length(E, Alen, ectx);
    Ecoeff = E->coeffs;
    Eexp = E->exps;

    Ei = 0;
    for (Ai = 0; Ai < Alen; Ai++)
    {
        nmod_mpolyn_redto_fq_nmod_mpolyn(Ecoeff + Ei, ectx, Acoeff + Ai, var, ctx);
        Eexp[Ei] = Aexp[Ai];
        Ei += ((Ecoeff + Ei)->length != 0);
    }
    E->length = Ei;

    FLINT_ASSERT(fq_nmod_mpolyun_is_canonical(E, ectx));
}




/*
    A = B using lowest degree representative
    A is in                      Fp [x_0, ..., x_(var-1), x_(var-1)][x_var]
    B is in (Fp[x_var]/alpha(x_var))[x_0, ..., x_(var-2)][x_(var-1)]
*/
void nmod_mpolyn_startinterp_fq_nmod_mpolyn(slong * lastdeg_,
                      nmod_mpolyn_t A, slong var, const nmod_mpoly_ctx_t ctx,
                            fq_nmod_mpolyn_t B, const fq_nmod_mpoly_ctx_t ectx)
{
    slong N = mpoly_words_per_exp_sp(B->bits, ctx->minfo);
    slong offset, shift;
    slong vi;
    fq_nmod_poly_struct * Bcoeff = B->coeffs;
    ulong * Bexp = B->exps;
    slong Blen = B->length;
    slong Bi;
    nmod_poly_struct * Acoeff;
    ulong * Aexp;
    slong Ai;
    slong lastdeg = -WORD(1);

    FLINT_ASSERT(A->bits == B->bits);

    nmod_mpolyn_fit_length(A, Blen, ctx);
    Acoeff = A->coeffs;
    Aexp = A->exps;

    mpoly_gen_offset_shift_sp(&offset, &shift, var - 1, A->bits, ctx->minfo);

    Ai = 0;
    for (Bi = 0; Bi < Blen; Bi++)
    {
        if (Ai + (Bcoeff + Bi)->length >= A->alloc)
        {
            nmod_mpolyn_fit_length(A, Ai + (Bcoeff + Bi)->length, ctx);
            Acoeff = A->coeffs;
            Aexp = A->exps;
        }
        for (vi = (Bcoeff + Bi)->length - 1; vi >= 0; vi--)
        {
            if (!fq_nmod_is_zero((Bcoeff + Bi)->coeffs + vi, ectx->fqctx))
            {
                mpoly_monomial_set_extra(Aexp + N*Ai, Bexp + N*Bi, N, offset, vi << shift);
                nmod_poly_set(Acoeff + Ai, (Bcoeff + Bi)->coeffs + vi);
                lastdeg = FLINT_MAX(lastdeg, nmod_poly_degree(Acoeff + Ai));
                Ai++;
            }
        }
    }
    A->length = Ai;

    *lastdeg_ = lastdeg;
}


/*
    A = B using lowest degree representative
    A is in                       Fp[X][x_0, ..., x_(var-2), x_(var-1)][x_var]
    B is in (Fp[x_var]/alpha(x_var))[X][x_0, ..., x_(var-2)][x_(var-1)]
*/
void nmod_mpolyun_startinterp_fq_nmod_mpolyun(slong * lastdeg,
                   nmod_mpolyun_t A, slong var, const nmod_mpoly_ctx_t ctx,
                            fq_nmod_mpolyun_t B, const fq_nmod_mpoly_ctx_t ectx)
{
    slong i;
    fq_nmod_mpolyn_struct * Bcoeff = B->coeffs;
    ulong * Bexp = B->exps;
    slong Blen = B->length;
    nmod_mpolyn_struct * Acoeff;
    ulong * Aexp;

    FLINT_ASSERT(var > 0);

    *lastdeg = -WORD(1);

    nmod_mpolyun_fit_length(A, Blen, ctx);
    Acoeff = A->coeffs;
    Aexp = A->exps;

    for (i = 0; i < Blen; i++)
    {
        Aexp[i] = Bexp[i];
        nmod_mpolyn_startinterp_fq_nmod_mpolyn(lastdeg, Acoeff + i, var, ctx, Bcoeff + i, ectx);
    }
    A->length = Blen;

    FLINT_ASSERT(nmod_mpolyun_is_canonical(A, ctx));
}




/*
    T = F + modulus*(A - F(x_var = alpha))
    no assumptions about matching monomials
    F is in Fp[x_0, ..., x_(var-1), x_(var-1)][x_var]
    A is in Fp[x_0, ..., x_(var-2)][x_(var-1)]
    in order to fxn correctly, modulus(alpha) should be 1
*/
int nmod_mpolyn_addinterp_fq_nmod_mpolyn(slong * lastdeg_,
             nmod_mpolyn_t T, nmod_mpolyn_t F, nmod_poly_t modulus, slong var, const nmod_mpoly_ctx_t ctx,
                              fq_nmod_mpolyn_t A, const fq_nmod_mpoly_ctx_t ectx)
{
    int changed = 0;
    slong N = mpoly_words_per_exp_sp(A->bits, ctx->minfo);
    slong offset, shift;
    slong lastdeg = -WORD(1);
    slong vi;
    fq_nmod_t u, v;
    nmod_poly_t w;

    nmod_poly_struct * Tcoeff;
    ulong * Texp;
    slong Ti;

    fq_nmod_poly_struct * Acoeff = A->coeffs;
    slong Alen = A->length;
    ulong * Aexp = A->exps;
    slong Ai;

    nmod_poly_struct * Fcoeff = F->coeffs;
    slong Flen = F->length;
    ulong * Fexp = F->exps;
    slong Fi;
    fq_nmod_t inv_m_eval;

    fq_nmod_init(inv_m_eval, ectx->fqctx);
    nmod_poly_rem(inv_m_eval, modulus, ectx->fqctx->modulus);
    fq_nmod_inv(inv_m_eval, inv_m_eval, ectx->fqctx);

    fq_nmod_init(u, ectx->fqctx);
    fq_nmod_init(v, ectx->fqctx);
    nmod_poly_init(w, ctx->ffinfo->mod.n);

    FLINT_ASSERT(var > 0);
    FLINT_ASSERT(T->bits == A->bits);
    FLINT_ASSERT(F->bits == A->bits);

    mpoly_gen_offset_shift_sp(&offset, &shift, var - 1, A->bits, ctx->minfo);

    Flen = F->length;

    nmod_mpolyn_fit_length(T, FLINT_MAX(Flen, Alen), ctx);
    Tcoeff = T->coeffs;
    Texp = T->exps;
    Ti = 0;

    Fi = Ai = vi = 0;
    if (Ai < Alen)
    {
        vi = fq_nmod_poly_degree(A->coeffs + Ai, ectx->fqctx);
    }
    while (Fi < Flen || Ai < Alen)
    {
        if (Ti >= T->alloc)
        {
            nmod_mpolyn_fit_length(T, Ti + FLINT_MAX(Flen - Fi, Alen - Ai), ctx);
            Tcoeff = T->coeffs;
            Texp = T->exps;
        }

        if (Fi < Flen && Ai < Alen && mpoly_monomial_equal_extra(Fexp + N*Fi, Aexp + N*Ai, N, offset, vi << shift))
        {
            /* F term ok, A term ok */
            nmod_poly_rem(u, Fcoeff + Fi, ectx->fqctx->modulus);
            fq_nmod_sub(v, (Acoeff + Ai)->coeffs + vi, u, ectx->fqctx);
            if (!fq_nmod_is_zero(v, ectx->fqctx))
            {
                changed = 1;
                fq_nmod_mul(u, v, inv_m_eval, ectx->fqctx);
                nmod_poly_mul(w, modulus, u);
                nmod_poly_add(Tcoeff + Ti, Fcoeff + Fi, w);
            }
            else
            {
                nmod_poly_set(Tcoeff + Ti, Fcoeff + Fi);
            }
            mpoly_monomial_set(Texp + N*Ti, Fexp + N*Fi, N);

            Fi++;
            do {
                vi--;
            } while (vi >= 0 && fq_nmod_is_zero((Acoeff + Ai)->coeffs + vi, ectx->fqctx));
            if (vi < 0)
            {
                Ai++;
                if (Ai < Alen)
                {
                    vi = fq_nmod_poly_degree(A->coeffs + Ai, ectx->fqctx);
                }
            }
        }
        else if (Fi < Flen && (Ai >= Alen || mpoly_monomial_gt_nomask_extra(Fexp + N*Fi, Aexp + N*Ai, N, offset, vi << shift)))
        {
            /* F term ok, A term missing */
            nmod_poly_rem(v, Fcoeff + Fi, ectx->fqctx->modulus);
            if (!fq_nmod_is_zero(v, ectx->fqctx))
            {
                changed = 1;
                fq_nmod_mul(u, v, inv_m_eval, ectx->fqctx);
                nmod_poly_mul(w, u, modulus);
                nmod_poly_sub(Tcoeff + Ti, Fcoeff + Fi, w);
            }
            else
            {
                nmod_poly_set(Tcoeff + Ti, Fcoeff + Fi);
            }
            mpoly_monomial_set(Texp + N*Ti, Fexp + N*Fi, N);

            Fi++;
        }
        else
        {
            FLINT_ASSERT(Ai < Alen && (Fi >= Flen || mpoly_monomial_lt_nomask_extra(Fexp + N*Fi, Aexp + N*Ai, N, offset, vi << shift)));

            /* F term missing, A term ok */
            changed = 1;
            fq_nmod_mul(u, (Acoeff + Ai)->coeffs + vi, inv_m_eval, ectx->fqctx);
            nmod_poly_mul(Tcoeff + Ti, modulus, u);
            mpoly_monomial_set_extra(Texp + N*Ti, Aexp + N*Ai, N, offset, vi << shift);

            do {
                vi--;
            } while (vi >= 0 && fq_nmod_is_zero((Acoeff + Ai)->coeffs + vi, ectx->fqctx));
            if (vi < 0)
            {
                Ai++;
                if (Ai < Alen)
                {
                    vi = fq_nmod_poly_degree(A->coeffs + Ai, ectx->fqctx);
                }
            }
        }

        FLINT_ASSERT(!nmod_poly_is_zero(Tcoeff + Ti));
        lastdeg = FLINT_MAX(lastdeg, nmod_poly_degree(Tcoeff + Ti));
        Ti++;
    }
    T->length = Ti;

    fq_nmod_clear(inv_m_eval, ectx->fqctx);
    fq_nmod_clear(u, ectx->fqctx);
    fq_nmod_clear(v, ectx->fqctx);
    nmod_poly_clear(w);

    *lastdeg_ = FLINT_MAX(*lastdeg_, lastdeg);
    return changed;
}

/*
    F = F + modulus*((A - F(alpha))/modulus(alpha))
    no assumptions about matching monomials
    F is in               Fp [X][x_0, ..., x_(var-1), x_(var-1)][x_var]
    A is in (Fp/alpha(x_var))[X][x_0, ..., x_(var-2)][x_(var-1)]
    alpha is ectx->fqctx->modulus
*/
int nmod_mpolyun_addinterp_fq_nmod_mpolyun(slong * lastdeg,
             nmod_mpolyun_t F, nmod_mpolyun_t T, nmod_poly_t modulus, slong var, const nmod_mpoly_ctx_t ctx,
                                 fq_nmod_mpolyun_t A, const fq_nmod_mpoly_ctx_t ectx)
{
    int changed = 0;
    slong i, j, k;
    ulong * Texp;
    ulong * Fexp;
    ulong * Aexp;
    slong Flen;
    slong Alen;
    nmod_mpolyn_struct * Tcoeff;
    nmod_mpolyn_struct * Fcoeff;
    nmod_mpolyn_t zero;
    fq_nmod_mpolyn_struct  * Acoeff;
    fq_nmod_mpolyn_t ezero;

    FLINT_ASSERT(var > 0);

    *lastdeg = -WORD(1);

    FLINT_ASSERT(F->bits == T->bits);
    FLINT_ASSERT(T->bits == A->bits);

    Flen = F->length;
    Alen = A->length;
    nmod_mpolyun_fit_length(T, Flen + Alen, ctx);

    Tcoeff = T->coeffs;
    Fcoeff = F->coeffs;
    Acoeff = A->coeffs;
    Texp = T->exps;
    Fexp = F->exps;
    Aexp = A->exps;   

    nmod_mpolyn_init(zero, A->bits, ctx);
    zero->length = 0;
    fq_nmod_mpolyn_init(ezero, A->bits, ectx);
    ezero->length = 0;

    i = j = k = 0;
    while (i < Flen || j < Alen)
    {
        if (i < Flen && j < Alen && (Fexp[i] == Aexp[j]))
        {
            /* F term ok, A term ok */
            changed |= nmod_mpolyn_addinterp_fq_nmod_mpolyn(lastdeg, Tcoeff + k, Fcoeff + i, modulus, var, ctx, Acoeff + j, ectx);
            Texp[k] = Aexp[j];
            i++;
            j++;
        }
        else if (i < Flen && (j >= Alen || Fexp[i] > Aexp[j]))
        {
            /* F term ok, A term missing */
            changed |= nmod_mpolyn_addinterp_fq_nmod_mpolyn(lastdeg, Tcoeff + k, Fcoeff + i, modulus, var, ctx, ezero, ectx);
            Texp[k] = Fexp[i];
            i++;
        }
        else
        {
            FLINT_ASSERT(j < Alen && (i >= Flen || Aexp[j] > Fexp[i]));

            /* F term missing, A term ok */
            changed |= nmod_mpolyn_addinterp_fq_nmod_mpolyn(lastdeg, Tcoeff + k, zero, modulus, var, ctx, Acoeff + j, ectx);
            Texp[k] = Aexp[j];
            j++;
        }

        FLINT_ASSERT(!nmod_mpolyn_is_zero(Tcoeff + k, ctx));
        k++;
    }
    T->length = k;

    if (changed)
    {
        nmod_mpolyun_swap(T, F);
    }

    nmod_mpolyn_clear(zero, ctx);
    fq_nmod_mpolyn_clear(ezero, ectx);

    FLINT_ASSERT(nmod_mpolyun_is_canonical(F, ctx));

    return changed;
}







int nmod_mpolyun_gcd_brown_lgprime(nmod_mpolyun_t G,
     nmod_mpolyun_t Abar, nmod_mpolyun_t Bbar, nmod_mpolyun_t A, nmod_mpolyun_t B,
                                         slong var, const nmod_mpoly_ctx_t ctx)
{
    int success;
    slong bound;
    slong offset, shift;
    fq_nmod_t temp, gammaeval;
    fq_nmod_mpolyun_t Aeval, Beval, Geval, Abareval, Bbareval;
    nmod_mpolyun_t T;
    slong deggamma, ldegG, ldegAbar, ldegBbar, ldegA, ldegB;
    nmod_poly_t cA, cB, cG, cAbar, cBbar, gamma;
    nmod_poly_t modulus;
    mp_bitcnt_t bits = A->bits;
    slong N = mpoly_words_per_exp_sp(bits, ctx->minfo);
    slong deg;
    fq_nmod_mpoly_ctx_t ectx;
#if WANT_ASSERT
    nmod_poly_t leadA, leadB;
#endif

    if (var == WORD(0))
    {
        return nmod_mpolyun_gcd_brown_lgprime_bivar(G, Abar, Bbar, A, B, ctx);
    }

    mpoly_gen_offset_shift_sp(&offset, &shift, var - 1, G->bits, ctx->minfo);

#if WANT_ASSERT
    nmod_poly_init(leadA, ctx->ffinfo->mod.n);
    nmod_poly_init(leadB, ctx->ffinfo->mod.n);
    nmod_poly_set(leadA, nmod_mpolyun_leadcoeff_ref(A, ctx));
    nmod_poly_set(leadB, nmod_mpolyun_leadcoeff_ref(B, ctx));
#endif

    nmod_poly_init(cA, ctx->ffinfo->mod.n);
    nmod_poly_init(cB, ctx->ffinfo->mod.n);
    nmod_mpolyun_content_last(cA, A, ctx);
    nmod_mpolyun_content_last(cB, B, ctx);
    nmod_mpolyun_divexact_last(A, cA, ctx);
    nmod_mpolyun_divexact_last(B, cB, ctx);

    nmod_poly_init(cG, ctx->ffinfo->mod.n);
    nmod_poly_gcd(cG, cA, cB);

    nmod_poly_init(cAbar, ctx->ffinfo->mod.n);
    nmod_poly_init(cBbar, ctx->ffinfo->mod.n);
    nmod_poly_div(cAbar, cA, cG);
    nmod_poly_div(cBbar, cB, cG);

    nmod_poly_init(gamma, ctx->ffinfo->mod.n);
    nmod_poly_gcd(gamma, nmod_mpolyun_leadcoeff_ref(A, ctx),
                         nmod_mpolyun_leadcoeff_ref(B, ctx));

    ldegA = nmod_mpolyun_lastdeg(A, ctx);
    ldegB = nmod_mpolyun_lastdeg(B, ctx);
    deggamma = nmod_poly_degree(gamma);

    bound = 1 + deggamma + FLINT_MAX(ldegA, ldegB);

    nmod_mpolyun_init(T, bits, ctx);

    nmod_poly_init(modulus, ctx->ffinfo->mod.n);
    nmod_poly_one(modulus);

    deg = WORD(20)/(FLINT_BIT_COUNT(ctx->ffinfo->mod.n));
    deg = FLINT_MAX(WORD(2), deg);

    fq_nmod_mpoly_ctx_init_deg(ectx, ctx->minfo->nvars, ORD_LEX,
                                                      ctx->ffinfo->mod.n, deg);

    fq_nmod_mpolyun_init(Aeval, bits, ectx);
    fq_nmod_mpolyun_init(Beval, bits, ectx);
    fq_nmod_mpolyun_init(Geval, bits, ectx);
    fq_nmod_mpolyun_init(Abareval, bits, ectx);
    fq_nmod_mpolyun_init(Bbareval, bits, ectx);
    fq_nmod_init(gammaeval, ectx->fqctx);
    fq_nmod_init(temp, ectx->fqctx);

    /* initialization already picked a prime */
    goto have_prime;

choose_prime: /* prime will be irreducible element of Fp[v] */

    /* same TODO */
    deg++;
    if (deg > 10000)
    {
        /* ran out of primes */
        success = 0;
        goto cleanup;
    }

    fq_nmod_mpoly_ctx_change_modulus(ectx, deg);

have_prime:

    /* make sure reduction does not kill both lc(A) and lc(B) */
    nmod_poly_rem(gammaeval, gamma, ectx->fqctx->modulus);
    if (fq_nmod_is_zero(gammaeval, ectx->fqctx))
    {
        goto choose_prime;
    }

    /* make sure reduction does not kill either A or B */
    nmod_mpolyun_redto_fq_nmod_mpolyun(Aeval, ectx, A, var, ctx);
    nmod_mpolyun_redto_fq_nmod_mpolyun(Beval, ectx, B, var, ctx);
    FLINT_ASSERT(fq_nmod_mpolyun_is_canonical(Aeval, ectx));
    FLINT_ASSERT(fq_nmod_mpolyun_is_canonical(Beval, ectx));
    if (Aeval->length == 0 || Beval->length == 0)
    {
        goto choose_prime;
    }

    success = fq_nmod_mpolyun_gcd_brown_smprime(Geval, Abareval, Bbareval,
                                                  Aeval, Beval, var - 1, ectx);
    if (success == 0)
    {
        goto choose_prime;
    }

    if (fq_nmod_mpolyun_is_nonzero_fq_nmod(Geval, ectx))
    {
        nmod_mpolyun_one(G, ctx);
        nmod_mpolyun_swap(Abar, A);
        nmod_mpolyun_swap(Bbar, B);
        goto successful_put_content;    
    }

    if (nmod_poly_degree(modulus) > 0)
    {
        /* compare leading monomials of Geval and G */
        int cmp = 0;
        FLINT_ASSERT(G->length > 0);
        if (G->exps[0] != Geval->exps[0])
        {
            cmp = G->exps[0] > Geval->exps[0] ? 1 : -1;
        }
        if (cmp == 0)
        {
            slong k = fq_nmod_poly_degree((Geval->coeffs + 0)->coeffs + 0, ectx->fqctx);
            cmp = mpoly_monomial_cmp_nomask_extra(
                        (G->coeffs + 0)->exps + N*0,
                    (Geval->coeffs + 0)->exps + N*0, N, offset, k << shift);
        }

        if (cmp < 0)
        {
            goto choose_prime;
        }
        else if (cmp > 0)
        {
            nmod_poly_one(modulus);
        }
    }

    fq_nmod_inv(temp, fq_nmod_mpolyn_leadcoeff_last_ref(Geval->coeffs + 0, ectx), ectx->fqctx);
    fq_nmod_mul(temp, temp, gammaeval, ectx->fqctx);
    fq_nmod_mpolyun_scalar_mul_fq_nmod(Geval, temp, ectx);

    if (nmod_poly_degree(modulus) > 0)
    {
        nmod_mpolyun_addinterp_fq_nmod_mpolyun(&ldegG, G, T, modulus, var, ctx, Geval, ectx);
        nmod_mpolyun_addinterp_fq_nmod_mpolyun(&ldegAbar, Abar, T, modulus, var, ctx, Abareval, ectx);
        nmod_mpolyun_addinterp_fq_nmod_mpolyun(&ldegBbar, Bbar, T, modulus, var, ctx, Bbareval, ectx);
    }
    else
    {
        nmod_mpolyun_startinterp_fq_nmod_mpolyun(&ldegG, G, var, ctx, Geval, ectx);
        nmod_mpolyun_startinterp_fq_nmod_mpolyun(&ldegAbar, Abar, var, ctx, Abareval, ectx);
        nmod_mpolyun_startinterp_fq_nmod_mpolyun(&ldegBbar, Bbar, var, ctx, Bbareval, ectx);
    }
    nmod_poly_mul(modulus, modulus, ectx->fqctx->modulus);

    if (nmod_poly_degree(modulus) < bound)
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

    nmod_poly_one(modulus);
    goto choose_prime;

successful:

    nmod_mpolyun_content_last(modulus, G, ctx);
    nmod_mpolyun_divexact_last(G, modulus, ctx);
    nmod_mpolyun_divexact_last(Abar, nmod_mpolyun_leadcoeff_ref(G, ctx), ctx);
    nmod_mpolyun_divexact_last(Bbar, nmod_mpolyun_leadcoeff_ref(G, ctx), ctx);

successful_put_content:

    nmod_mpolyun_mul_last(G, cG, ctx);
    nmod_mpolyun_mul_last(Abar, cAbar, ctx);
    nmod_mpolyun_mul_last(Bbar, cBbar, ctx);

    success = 1;

cleanup:

#if WANT_ASSERT
    if (success)
    {
        FLINT_ASSERT(1 == nmod_mpolyun_leadcoeff_last(G, ctx));
        nmod_poly_mul(modulus, nmod_mpolyun_leadcoeff_ref(G, ctx),
                               nmod_mpolyun_leadcoeff_ref(Abar, ctx));
        FLINT_ASSERT(nmod_poly_equal(modulus, leadA));
        nmod_poly_mul(modulus, nmod_mpolyun_leadcoeff_ref(G, ctx),
                               nmod_mpolyun_leadcoeff_ref(Bbar, ctx));
        FLINT_ASSERT(nmod_poly_equal(modulus, leadB));
    }
    nmod_poly_clear(leadA);
    nmod_poly_clear(leadB);
#endif

    nmod_poly_clear(cA);
    nmod_poly_clear(cB);
    nmod_poly_clear(cG);
    nmod_poly_clear(cAbar);
    nmod_poly_clear(cBbar);

    nmod_poly_clear(gamma);

    nmod_mpolyun_clear(T, ctx);

    nmod_poly_clear(modulus);

    fq_nmod_mpolyun_clear(Aeval, ectx);
    fq_nmod_mpolyun_clear(Beval, ectx);
    fq_nmod_mpolyun_clear(Geval, ectx);
    fq_nmod_mpolyun_clear(Abareval, ectx);
    fq_nmod_mpolyun_clear(Bbareval, ectx);
    fq_nmod_clear(gammaeval, ectx->fqctx);
    fq_nmod_clear(temp, ectx->fqctx);

    fq_nmod_mpoly_ctx_clear(ectx);

    return success;
}
