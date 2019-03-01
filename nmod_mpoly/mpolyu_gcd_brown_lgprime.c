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

int nmod_mpolyu_gcd_brown_lgprime(nmod_mpolyu_t G,
     nmod_mpolyu_t Abar, nmod_mpolyu_t Bbar, nmod_mpolyu_t A, nmod_mpolyu_t B,
                                         slong var, const nmod_mpoly_ctx_t ctx)
{
    int success;
    slong bound;
    fq_nmod_t temp, gammaeval;
    fq_nmod_mpolyu_t Aeval, Beval, Geval, Abareval, Bbareval;
    nmod_mpolyun_t An, Bn, Gs, Abars, Bbars, T;
    slong deggamma, ldegGs, ldegAbars, ldegBbars, ldegA, ldegB;
    nmod_poly_t cA, cB, cG, cAbar, cBbar, gamma;
    nmod_poly_t cGs, cAbars, cBbars, modulus;
    mp_bitcnt_t bits = A->bits;
    slong N = mpoly_words_per_exp(bits, ctx->minfo);
    slong deg;
    fq_nmod_mpoly_ctx_t ectx;

    if (var == -WORD(1))
    {
        /* no more variables left to interpolate */
        return nmod_mpolyu_gcd_brown_univar(G, Abar, Bbar, A, B, ctx);
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
    nmod_poly_gcd(cG, cA, cB);

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

    nmod_mpolyun_init(Gs, bits, ctx);
    nmod_mpolyun_init(Abars, bits, ctx);
    nmod_mpolyun_init(Bbars, bits, ctx);
    nmod_mpolyun_init(T, bits, ctx);

    nmod_poly_init(modulus, ctx->ffinfo->mod.n);
    nmod_poly_one(modulus);

    deg = WORD(20)/(FLINT_BIT_COUNT(ctx->ffinfo->mod.n));
    deg = FLINT_MAX(WORD(2), deg);

    fq_nmod_mpoly_ctx_init_deg(ectx, ctx->minfo->nvars, ORD_LEX,
                                                      ctx->ffinfo->mod.n, deg);

    fq_nmod_mpolyu_init(Aeval, bits, ectx);
    fq_nmod_mpolyu_init(Beval, bits, ectx);
    fq_nmod_mpolyu_init(Geval, bits, ectx);
    fq_nmod_mpolyu_init(Abareval, bits, ectx);
    fq_nmod_mpolyu_init(Bbareval, bits, ectx);
    fq_nmod_init(gammaeval, ectx->fqctx);
    fq_nmod_init(temp, ectx->fqctx);

    /* initialization already picked a prime */
    goto have_prime;

choose_prime:

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
    nmod_mpolyun_redto_fq_nmod_mpolyu(Aeval, An, ectx, ctx);
    nmod_mpolyun_redto_fq_nmod_mpolyu(Beval, Bn, ectx, ctx);
    if (Aeval->length == 0 || Beval->length == 0)
    {
        goto choose_prime;
    }

    FLINT_ASSERT(fq_nmod_mpolyu_is_canonical(Aeval, ectx));
    FLINT_ASSERT(fq_nmod_mpolyu_is_canonical(Beval, ectx));

    success = fq_nmod_mpolyu_gcd_brown_smprime(Geval, Abareval, Bbareval,
                                                  Aeval, Beval, var - 1, ectx);

    if (success == 0)
    {
        goto choose_prime;
    }

    if (fq_nmod_mpolyu_is_nonzero_fq_nmod(Geval, ectx))
    {
        nmod_mpolyun_one(Gs, ctx);
        nmod_mpolyun_set(Abars, An, ctx);
        nmod_mpolyun_set(Bbars, Bn, ctx);
        goto successful;    
    }

    if (nmod_poly_degree(modulus) > 0)
    {
        /* compare leading monomials of Geval and Gs */
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
            goto choose_prime;
        }
        else if (cmp < 0)
        {
            nmod_poly_one(modulus);
        }
    }

    fq_nmod_inv(temp, fq_nmod_mpolyu_leadcoeff_ref(Geval, ectx), ectx->fqctx);
    fq_nmod_mul(temp, temp, gammaeval, ectx->fqctx);
    fq_nmod_mpolyu_scalar_mul_fq_nmod(Geval, temp, ectx);

    if (nmod_poly_degree(modulus) > 0)
    {
        nmod_mpolyun_addinterp_fq_nmod_mpolyu(&ldegGs, Gs, T, modulus, ctx, Geval, ectx);
        nmod_mpolyun_addinterp_fq_nmod_mpolyu(&ldegAbars, Abars, T, modulus, ctx, Abareval, ectx);
        nmod_mpolyun_addinterp_fq_nmod_mpolyu(&ldegBbars, Bbars, T, modulus, ctx, Bbareval, ectx);
    }
    else
    {
        nmod_mpolyun_set_fq_nmod_mpolyu(&ldegGs, Gs, ctx, Geval, ectx);
        nmod_mpolyun_set_fq_nmod_mpolyu(&ldegAbars, Abars, ctx, Abareval, ectx);
        nmod_mpolyun_set_fq_nmod_mpolyu(&ldegBbars, Bbars, ctx, Bbareval, ectx);
    }
    nmod_poly_mul(modulus, modulus, ectx->fqctx->modulus);

    if (nmod_poly_degree(modulus) < bound)
    {
        goto choose_prime;
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

    nmod_poly_one(modulus);
    goto choose_prime;

cleanup:

    nmod_mpolyun_clear(An, ctx);
    nmod_mpolyun_clear(Bn, ctx);

    nmod_poly_clear(cA);
    nmod_poly_clear(cB);

    nmod_poly_clear(cG);

    nmod_poly_clear(cAbar);
    nmod_poly_clear(cBbar);

    nmod_poly_clear(gamma);

    nmod_mpolyun_clear(Gs, ctx);
    nmod_mpolyun_clear(Abars, ctx);
    nmod_mpolyun_clear(Bbars, ctx);
    nmod_mpolyun_clear(T, ctx);

    nmod_poly_clear(modulus);

    fq_nmod_mpolyu_clear(Aeval, ectx);
    fq_nmod_mpolyu_clear(Beval, ectx);
    fq_nmod_mpolyu_clear(Geval, ectx);
    fq_nmod_mpolyu_clear(Abareval, ectx);
    fq_nmod_mpolyu_clear(Bbareval, ectx);
    fq_nmod_clear(gammaeval, ectx->fqctx);
    fq_nmod_clear(temp, ectx->fqctx);

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
