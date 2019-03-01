/*
    Copyright (C) 2019 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fq_nmod_mpoly.h"


int fq_nmod_mpolyu_gcd_brown_univar(fq_nmod_mpolyu_t G,
               fq_nmod_mpolyu_t Abar, fq_nmod_mpolyu_t Bbar,
         fq_nmod_mpolyu_t A, fq_nmod_mpolyu_t B, const fq_nmod_mpoly_ctx_t ctx)
{
    fq_nmod_poly_t a, b, g, abar, bbar;
    fq_nmod_poly_t r;
    fq_nmod_poly_init(r, ctx->fqctx);
    FLINT_ASSERT(A->bits == B->bits);
    fq_nmod_poly_init(a, ctx->fqctx);
    fq_nmod_poly_init(b, ctx->fqctx);
    fq_nmod_poly_init(g, ctx->fqctx);
    fq_nmod_poly_init(abar, ctx->fqctx);
    fq_nmod_poly_init(bbar, ctx->fqctx);
    fq_nmod_mpolyu_cvtto_poly(a, A, ctx);
    fq_nmod_mpolyu_cvtto_poly(b, B, ctx);
    fq_nmod_poly_gcd(g, a, b, ctx->fqctx);
    fq_nmod_poly_divrem(abar, r, a, g, ctx->fqctx);
    FLINT_ASSERT(fq_nmod_poly_is_zero(r, ctx->fqctx));
    fq_nmod_poly_divrem(bbar, r, b, g, ctx->fqctx);
    FLINT_ASSERT(fq_nmod_poly_is_zero(r, ctx->fqctx));
    fq_nmod_mpolyu_cvtfrom_poly(G, g, ctx);
    fq_nmod_mpolyu_cvtfrom_poly(Abar, abar, ctx);
    fq_nmod_mpolyu_cvtfrom_poly(Bbar, bbar, ctx);
    fq_nmod_poly_clear(a, ctx->fqctx);
    fq_nmod_poly_clear(b, ctx->fqctx);
    fq_nmod_poly_clear(g, ctx->fqctx);
    fq_nmod_poly_clear(abar, ctx->fqctx);
    fq_nmod_poly_clear(bbar, ctx->fqctx);
    fq_nmod_poly_clear(r, ctx->fqctx);
    return 1;
}


int fq_nmod_mpolyu_gcd_brown_bivar(fq_nmod_mpolyu_t G,
                fq_nmod_mpolyu_t Abar, fq_nmod_mpolyu_t Bbar,
         fq_nmod_mpolyu_t A, fq_nmod_mpolyu_t B, const fq_nmod_mpoly_ctx_t ctx)
{
    return 0;
}

int fq_nmod_mpolyu_gcd_brown_smprime(fq_nmod_mpolyu_t G,
          fq_nmod_mpolyu_t Abar, fq_nmod_mpolyu_t Bbar, fq_nmod_mpolyu_t A,
                  fq_nmod_mpolyu_t B, slong var, const fq_nmod_mpoly_ctx_t ctx)
{
    int success;
    slong bound;
    fq_nmod_t alpha, temp, gammaeval;
    fq_nmod_mpolyu_t Aeval, Beval, Geval, Abareval, Bbareval;
    fq_nmod_mpolyun_t An, Bn, Gs, Abars, Bbars, T;
    slong deggamma, ldegGs, ldegAbars, ldegBbars, ldegA, ldegB;
    fq_nmod_poly_t cA, cB, cG, cAbar, cBbar, gamma, trem;
    fq_nmod_poly_t cGs, cAbars, cBbars, modulus, tempmod;

    mp_bitcnt_t bits = A->bits;
    slong N = mpoly_words_per_exp(bits, ctx->minfo);
    if (var == -WORD(1))
    {
        /* no more variables left to interpolate */
        return fq_nmod_mpolyu_gcd_brown_univar(G, Abar, Bbar, A, B, ctx);
    }

    if (0 && var == WORD(0))
    {
        /* bivariate is more comfortable separated */
        return fq_nmod_mpolyu_gcd_brown_bivar(G, Abar, Bbar, A, B, ctx);
    }

    fq_nmod_mpolyun_init(An, bits, ctx);
    fq_nmod_mpolyun_init(Bn, bits, ctx);
    fq_nmod_mpolyu_cvtto_mpolyun(An, A, var, ctx);
    fq_nmod_mpolyu_cvtto_mpolyun(Bn, B, var, ctx);

    fq_nmod_poly_init(cA, ctx->fqctx);
    fq_nmod_poly_init(cB, ctx->fqctx);
    fq_nmod_mpolyun_content_last(cA, An, ctx);
    fq_nmod_mpolyun_content_last(cB, Bn, ctx);
    fq_nmod_mpolyun_divexact_last(An, cA, ctx);
    fq_nmod_mpolyun_divexact_last(Bn, cB, ctx);

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
    fq_nmod_poly_gcd(gamma, fq_nmod_mpolyun_leadcoeff_ref(An, ctx),
                            fq_nmod_mpolyun_leadcoeff_ref(Bn, ctx), ctx->fqctx);

    ldegA = fq_nmod_mpolyun_lastdeg(An, ctx);
    ldegB = fq_nmod_mpolyun_lastdeg(Bn, ctx);
    deggamma = fq_nmod_poly_degree(gamma, ctx->fqctx);

    bound = 1 + deggamma + FLINT_MAX(ldegA, ldegB);

    fq_nmod_mpolyun_init(Gs, bits, ctx);
    fq_nmod_mpolyun_init(Abars, bits, ctx);
    fq_nmod_mpolyun_init(Bbars, bits, ctx);
    fq_nmod_mpolyun_init(T, bits, ctx);

    fq_nmod_mpolyu_init(Aeval, bits, ctx);
    fq_nmod_mpolyu_init(Beval, bits, ctx);
    fq_nmod_mpolyu_init(Geval, bits, ctx);
    fq_nmod_mpolyu_init(Abareval, bits, ctx);
    fq_nmod_mpolyu_init(Bbareval, bits, ctx);

    fq_nmod_init(gammaeval, ctx->fqctx);
    fq_nmod_init(alpha, ctx->fqctx);
    fq_nmod_init(temp, ctx->fqctx);

    fq_nmod_poly_init(modulus, ctx->fqctx);
    fq_nmod_poly_one(modulus, ctx->fqctx);
    fq_nmod_poly_init(tempmod, ctx->fqctx);
    fq_nmod_poly_gen(tempmod, ctx->fqctx);
    fq_nmod_poly_neg(tempmod, tempmod, ctx->fqctx);

    fq_nmod_set_ui(alpha, 0, ctx->fqctx);
    while (fq_nmod_next(alpha, ctx->fqctx) != 0)
    {
        /* make sure evaluation point does not kill both lc(A) and lc(B) */
        fq_nmod_poly_evaluate_fq_nmod(gammaeval, gamma, alpha, ctx->fqctx);
        if (fq_nmod_is_zero(gammaeval, ctx->fqctx))
        {
            goto break_continue;
        }

        /* make sure evaluation point does not kill either A or B */
        fq_nmod_mpolyun_eval_last(Aeval, An, alpha, ctx);
        fq_nmod_mpolyun_eval_last(Beval, Bn, alpha, ctx);
        if (Aeval->length == 0 || Beval->length == 0)
        {
            goto break_continue;
        }

        success = fq_nmod_mpolyu_gcd_brown_smprime(Geval, Abareval, Bbareval,
                                                   Aeval, Beval, var - 1, ctx);

        if (success == 0)
        {
            goto break_continue;
        }

        FLINT_ASSERT(Geval->length > 0);
        FLINT_ASSERT(Abareval->length);
        FLINT_ASSERT(Bbareval->length > 0);

        if (fq_nmod_mpolyu_is_nonzero_fq_nmod(Geval, ctx))
        {
            fq_nmod_mpolyun_one(Gs, ctx);
            fq_nmod_mpolyun_set(Abars, An, ctx);
            fq_nmod_mpolyun_set(Bbars, Bn, ctx);
            goto successful;    
        }

        if (fq_nmod_poly_degree(modulus, ctx->fqctx) > 0)
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
                fq_nmod_poly_one(modulus, ctx->fqctx);
            }
        }

        fq_nmod_inv(temp, fq_nmod_mpolyu_leadcoeff_ref(Geval, ctx), ctx->fqctx);
        fq_nmod_mul(temp, temp, gammaeval, ctx->fqctx);
        fq_nmod_mpolyu_scalar_mul_fq_nmod(Geval, temp, ctx);

        if (fq_nmod_poly_degree(modulus, ctx->fqctx) > 0)
        {
            fq_nmod_poly_evaluate_fq_nmod(temp, modulus, alpha, ctx->fqctx);
            fq_nmod_inv(temp, temp, ctx->fqctx);
            fq_nmod_poly_scalar_mul_fq_nmod(modulus, modulus, temp, ctx->fqctx);
            fq_nmod_mpolyun_addinterp(&ldegGs, Gs, T, Geval, modulus, alpha, ctx);
            fq_nmod_mpolyun_addinterp(&ldegAbars, Abars, T, Abareval, modulus, alpha, ctx);
            fq_nmod_mpolyun_addinterp(&ldegBbars, Bbars, T, Bbareval, modulus, alpha, ctx);
        }
        else
        {
            fq_nmod_mpolyun_set_mpolyu(Gs, Geval, ctx);
            fq_nmod_mpolyun_set_mpolyu(Abars, Abareval, ctx);
            fq_nmod_mpolyun_set_mpolyu(Bbars, Bbareval, ctx);
            ldegGs = 0;
            ldegAbars = 0;
            ldegBbars = 0;
        }
        fq_nmod_poly_set_coeff(tempmod, 0, alpha, ctx->fqctx);
        fq_nmod_poly_mul(modulus, modulus, tempmod, ctx->fqctx);

        if (fq_nmod_poly_degree(modulus, ctx->fqctx) < bound)
        {
            goto break_continue;
        }

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
            fq_nmod_poly_one(modulus, ctx->fqctx);
            goto break_continue;
        }

break_continue:
        (void)(NULL);
    }

    success = 0;

cleanup:

    fq_nmod_mpolyun_clear(An, ctx);
    fq_nmod_mpolyun_clear(Bn, ctx);

    fq_nmod_poly_clear(cA, ctx->fqctx);
    fq_nmod_poly_clear(cB, ctx->fqctx);

    fq_nmod_poly_clear(cG, ctx->fqctx);

    fq_nmod_poly_clear(cAbar, ctx->fqctx);
    fq_nmod_poly_clear(cBbar, ctx->fqctx);
    fq_nmod_poly_clear(trem, ctx->fqctx);
    fq_nmod_poly_clear(gamma, ctx->fqctx);

    fq_nmod_mpolyu_clear(Aeval, ctx);
    fq_nmod_mpolyu_clear(Beval, ctx);
    fq_nmod_mpolyu_clear(Geval, ctx);
    fq_nmod_mpolyu_clear(Abareval, ctx);
    fq_nmod_mpolyu_clear(Bbareval, ctx);

    fq_nmod_mpolyun_clear(Gs, ctx);
    fq_nmod_mpolyun_clear(Abars, ctx);
    fq_nmod_mpolyun_clear(Bbars, ctx);
    fq_nmod_mpolyun_clear(T, ctx);

    fq_nmod_clear(gammaeval, ctx->fqctx);
    fq_nmod_clear(alpha, ctx->fqctx);
    fq_nmod_clear(temp, ctx->fqctx);

    fq_nmod_poly_clear(modulus, ctx->fqctx);
    fq_nmod_poly_clear(tempmod, ctx->fqctx);

    return success;

successful:

    fq_nmod_poly_init(cGs, ctx->fqctx);
    fq_nmod_poly_init(cAbars, ctx->fqctx);
    fq_nmod_poly_init(cBbars, ctx->fqctx);

    fq_nmod_mpolyun_content_last(cGs, Gs, ctx);
    fq_nmod_mpolyun_content_last(cAbars, Abars, ctx);
    fq_nmod_mpolyun_content_last(cBbars, Bbars, ctx);

    fq_nmod_mpolyun_divexact_last(Gs, cGs, ctx);
    fq_nmod_mpolyun_divexact_last(Abars, cAbars, ctx);
    fq_nmod_mpolyun_divexact_last(Bbars, cBbars, ctx);

    fq_nmod_mpolyun_mul_last(Gs, cG, ctx);
    fq_nmod_mpolyun_mul_last(Abars, cAbar, ctx);
    fq_nmod_mpolyun_mul_last(Bbars, cBbar, ctx);

    fq_nmod_mpolyu_cvtfrom_mpolyun(G, Gs, var, ctx);
    fq_nmod_mpolyu_cvtfrom_mpolyun(Abar, Abars, var, ctx);
    fq_nmod_mpolyu_cvtfrom_mpolyun(Bbar, Bbars, var, ctx);

    fq_nmod_poly_clear(cGs, ctx->fqctx);
    fq_nmod_poly_clear(cAbars, ctx->fqctx);
    fq_nmod_poly_clear(cBbars, ctx->fqctx);

    success = 1;
    goto cleanup;
}
