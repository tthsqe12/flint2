/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "nmod_mpoly.h"


int nmod_mpolyu_gcd_brown_univar(nmod_mpolyu_t G,
     nmod_mpolyu_t Abar, nmod_mpolyu_t Bbar, nmod_mpolyu_t A, nmod_mpolyu_t B,
                                                    const nmod_mpoly_ctx_t ctx)
{
    nmod_poly_t a, b, g, abar, bbar;
    FLINT_ASSERT(A->bits == B->bits);
    nmod_poly_init(a, ctx->ffinfo->mod.n);
    nmod_poly_init(b, ctx->ffinfo->mod.n);
    nmod_poly_init(g, ctx->ffinfo->mod.n);
    nmod_poly_init(abar, ctx->ffinfo->mod.n);
    nmod_poly_init(bbar, ctx->ffinfo->mod.n);
    nmod_mpolyu_cvtto_poly(a, A, ctx);
    nmod_mpolyu_cvtto_poly(b, B, ctx);
    nmod_poly_gcd(g, a, b);
    nmod_poly_div(abar, a, g);
    nmod_poly_div(bbar, b, g);
    nmod_mpolyu_cvtfrom_poly(G, g, ctx);
    nmod_mpolyu_cvtfrom_poly(Abar, abar, ctx);
    nmod_mpolyu_cvtfrom_poly(Bbar, bbar, ctx);
    nmod_poly_clear(a);
    nmod_poly_clear(b);
    nmod_poly_clear(g);
    nmod_poly_clear(abar);
    nmod_poly_clear(bbar);
    return 1;
    
}


int nmod_mpolyu_gcd_brown_bivar(nmod_mpolyu_t G,
     nmod_mpolyu_t Abar, nmod_mpolyu_t Bbar, nmod_mpolyu_t A, nmod_mpolyu_t B,
                                                    const nmod_mpoly_ctx_t ctx)
{
    return 0;
}

int nmod_mpolyu_gcd_brown_smprime(nmod_mpolyu_t G,
     nmod_mpolyu_t Abar, nmod_mpolyu_t Bbar, nmod_mpolyu_t A, nmod_mpolyu_t B,
                                         slong var, const nmod_mpoly_ctx_t ctx)
{
    int success;
    slong bound;
    mp_limb_t alpha, temp, gammaeval;
    nmod_mpolyu_t Aeval, Beval, Geval, Abareval, Bbareval;
    nmod_mpolyun_t An, Bn, Gs, Abars, Bbars, T;
    slong deggamma, ldegGs, ldegAbars, ldegBbars, ldegA, ldegB;
    nmod_poly_t cA, cB, cG, cAbar, cBbar, gamma;
    nmod_poly_t cGs, cAbars, cBbars, modulus, modulus2;

    mp_bitcnt_t bits = A->bits;
    slong N = mpoly_words_per_exp(bits, ctx->minfo);
/*
printf("nmod_mpolyu_gcd_brown_smprime called\n");
printf("A: "); nmod_mpolyu_print_pretty(A, NULL, ctx); printf("\n");
printf("B: "); nmod_mpolyu_print_pretty(B, NULL, ctx); printf("\n");
*/
    if (var == -WORD(1))
    {
        /* no more variables left to interpolate */
        return nmod_mpolyu_gcd_brown_univar(G, Abar, Bbar, A, B, ctx);
    }

    if (0 && var == WORD(0))
    {
        /* bivariate is more comfortable separated */
        return nmod_mpolyu_gcd_brown_bivar(G, Abar, Bbar, A, B, ctx);
    }

    nmod_mpolyun_init(An, bits, ctx);
    nmod_mpolyun_init(Bn, bits, ctx);
    nmod_mpolyu_cvtto_mpolyun(An, A, var, ctx);
    nmod_mpolyu_cvtto_mpolyun(Bn, B, var, ctx);

    nmod_poly_init(cA, ctx->ffinfo->mod.n);
    nmod_poly_init(cB, ctx->ffinfo->mod.n);
    nmod_mpolyun_content_last(cA, An, ctx);
    nmod_mpolyun_content_last(cB, Bn, ctx);
    nmod_mpolyun_divexact_last(An, cA, ctx);
    nmod_mpolyun_divexact_last(Bn, cB, ctx);

    nmod_poly_init(cG, ctx->ffinfo->mod.n);
    nmod_poly_gcd_euclidean(cG, cA, cB);

    nmod_poly_init(cAbar, ctx->ffinfo->mod.n);
    nmod_poly_init(cBbar, ctx->ffinfo->mod.n);
    nmod_poly_div(cAbar, cA, cG);
    nmod_poly_div(cBbar, cB, cG);

    nmod_poly_init(gamma, ctx->ffinfo->mod.n);
    nmod_poly_gcd(gamma, nmod_mpolyun_leadcoeff_ref(An, ctx),
                         nmod_mpolyun_leadcoeff_ref(Bn, ctx));

    ldegA = nmod_mpolyun_lastdeg(An, ctx);
    ldegB = nmod_mpolyun_lastdeg(Bn, ctx);
    deggamma = nmod_poly_degree(gamma);
    bound = 1 + deggamma + FLINT_MAX(ldegA, ldegB);

    nmod_mpolyu_init(Aeval, bits, ctx);
    nmod_mpolyu_init(Beval, bits, ctx);
    nmod_mpolyu_init(Geval, bits, ctx);
    nmod_mpolyu_init(Abareval, bits, ctx);
    nmod_mpolyu_init(Bbareval, bits, ctx);

    nmod_mpolyun_init(Gs, bits, ctx);
    nmod_mpolyun_init(Abars, bits, ctx);
    nmod_mpolyun_init(Bbars, bits, ctx);
    nmod_mpolyun_init(T, bits, ctx);

    nmod_poly_init(modulus, ctx->ffinfo->mod.n);
    nmod_poly_init(modulus2, ctx->ffinfo->mod.n);
    nmod_poly_one(modulus);

    if ((ctx->ffinfo->mod.n & UWORD(1)) == UWORD(0))
    {
        success = 0;
        goto cleanup;
    }

    for (alpha = (ctx->ffinfo->mod.n - UWORD(1))/UWORD(2); alpha != UWORD(0); alpha--)
    {
        /* make sure evaluation point does not kill both lc(A) and lc(B) */
        gammaeval = nmod_poly_evaluate_nmod(gamma, alpha);
        if (gammaeval == 0)
        {
            goto break_continue;
        }

        /* make sure evaluation point does not kill either A or B */
        nmod_mpolyun_eval_last(Aeval, An, alpha, ctx);
        nmod_mpolyun_eval_last(Beval, Bn, alpha, ctx);
        if (Aeval->length == 0 || Beval->length == 0)
        {
            goto break_continue;
        }

        success = nmod_mpolyu_gcd_brown_smprime(Geval, Abareval, Bbareval,
                                                   Aeval, Beval, var - 1, ctx);

        if (success == 0)
            goto break_continue;

        FLINT_ASSERT(Geval->length > 0);
        FLINT_ASSERT(Abareval->length);
        FLINT_ASSERT(Bbareval->length > 0);

        if (nmod_mpolyu_is_nonzero_nmod(Geval, ctx))
        {
            nmod_mpolyun_one(Gs, ctx);
            nmod_mpolyun_set(Abars, An, ctx);
            nmod_mpolyun_set(Bbars, Bn, ctx);
            goto successful;    
        }

        if (nmod_poly_degree(modulus) > 0)
        {
            int cmp = 0;
            FLINT_ASSERT(Gs->length > 0);
            if (Geval->exps[0] != Gs->exps[0])
            {
                cmp = Geval->exps[0] > Gs->exps[0] ? 1 : -1;
            }
            if (cmp == 0)
            {
                cmp = mpoly_monomial_cmp_nomask((Geval->coeffs + 0)->exps + N*0,
                                                (Gs->coeffs + 0)->exps + N*0, N);
            }

            if (cmp > 0)
            {
                goto break_continue;
            }
            else if (cmp < 0)
            {
                nmod_poly_one(modulus);
            }
        }

        temp = n_invmod(nmod_mpolyu_leadcoeff(Geval, ctx), ctx->ffinfo->mod.n);
        nmod_mpolyu_scalar_mul_nmod(Geval, nmod_mul(gammaeval, temp,
                                                   ctx->ffinfo->mod), ctx);

        if (nmod_poly_degree(modulus) > 0)
        {
            temp = nmod_poly_evaluate_nmod(modulus, alpha);
            temp = n_invmod(temp, ctx->ffinfo->mod.n);
            nmod_poly_scalar_mul_nmod(modulus, modulus, temp);
            nmod_mpolyun_addinterp(&ldegGs, Gs, T, Geval, modulus, alpha, ctx);
            nmod_mpolyun_addinterp(&ldegAbars, Abars, T, Abareval, modulus, alpha, ctx);
            nmod_mpolyun_addinterp(&ldegBbars, Bbars, T, Bbareval, modulus, alpha, ctx);
        }
        else
        {
            nmod_mpolyun_set_mpolyu(Gs, Geval, ctx);
            nmod_mpolyun_set_mpolyu(Abars, Abareval, ctx);
            nmod_mpolyun_set_mpolyu(Bbars, Bbareval, ctx);
            ldegGs = 0;
            ldegAbars = 0;
            ldegBbars = 0;
        }

        nmod_poly_scalar_mul_nmod(modulus2, modulus, alpha);
        nmod_poly_shift_left(modulus, modulus, 1);
        nmod_poly_sub(modulus, modulus, modulus2);

        if (nmod_poly_degree(modulus) < bound)
            continue;

        FLINT_ASSERT(ldegGs >= 0);
        FLINT_ASSERT(ldegAbars >= 0);
        FLINT_ASSERT(ldegBbars >= 0);

        if (   deggamma + ldegA == ldegGs + ldegAbars
            && deggamma + ldegB == ldegGs + ldegBbars
           )
        {
            goto successful;
        }
        else
        {
            nmod_poly_one(modulus);
            continue;
        }

break_continue:
        (void)(NULL);
    }

    success = 0;

cleanup:

    nmod_mpolyun_clear(An, ctx);
    nmod_mpolyun_clear(Bn, ctx);

    nmod_poly_clear(cA);
    nmod_poly_clear(cB);

    nmod_poly_clear(cG);

    nmod_poly_clear(cAbar);
    nmod_poly_clear(cBbar);

    nmod_poly_clear(gamma);

    nmod_mpolyu_clear(Aeval, ctx);
    nmod_mpolyu_clear(Beval, ctx);
    nmod_mpolyu_clear(Geval, ctx);
    nmod_mpolyu_clear(Abareval, ctx);
    nmod_mpolyu_clear(Bbareval, ctx);

    nmod_mpolyun_clear(Gs, ctx);
    nmod_mpolyun_clear(Abars, ctx);
    nmod_mpolyun_clear(Bbars, ctx);
    nmod_mpolyun_clear(T, ctx);

    nmod_poly_clear(modulus);
    nmod_poly_clear(modulus2);

    return success;

successful:

    nmod_poly_init(cGs, ctx->ffinfo->mod.n);
    nmod_poly_init(cAbars, ctx->ffinfo->mod.n);
    nmod_poly_init(cBbars, ctx->ffinfo->mod.n);

    nmod_mpolyun_content_last(cGs, Gs, ctx);
    nmod_mpolyun_content_last(cAbars, Abars, ctx);
    nmod_mpolyun_content_last(cBbars, Bbars, ctx);

    nmod_mpolyun_divexact_last(Gs, cGs, ctx);
    nmod_mpolyun_divexact_last(Abars, cAbars, ctx);
    nmod_mpolyun_divexact_last(Bbars, cBbars, ctx);

    nmod_mpolyun_mul_last(Gs, cG, ctx);
    nmod_mpolyun_mul_last(Abars, cAbar, ctx);
    nmod_mpolyun_mul_last(Bbars, cBbar, ctx);

    nmod_mpolyu_cvtfrom_mpolyun(G, Gs, var, ctx);
    nmod_mpolyu_cvtfrom_mpolyun(Abar, Abars, var, ctx);
    nmod_mpolyu_cvtfrom_mpolyun(Bbar, Bbars, var, ctx);

    nmod_poly_clear(cGs);
    nmod_poly_clear(cAbars);
    nmod_poly_clear(cBbars);

    success = 1;
    goto cleanup;
}
