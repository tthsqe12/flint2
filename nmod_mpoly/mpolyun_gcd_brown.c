/*
    Copyright (C) 2019 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "nmod_mpoly.h"
#include "fq_nmod_mpoly.h"


int nmod_mpolyun_gcd_brown_smprime_bivar_ref(
    nmod_mpolyun_t G,
    nmod_mpolyun_t Abar,
    nmod_mpolyun_t Bbar,
    nmod_mpolyun_t A,
    nmod_mpolyun_t B,
    const nmod_mpoly_ctx_t ctx)
{
    int success;
    slong bound;
    mp_limb_t alpha, temp, gammaeval;
    nmod_poly_t Aeval, Beval, Geval, Abareval, Bbareval;
    nmod_mpolyun_t T;
    slong deggamma, ldegG, ldegAbar, ldegBbar, ldegA, ldegB;
    nmod_poly_t cA, cB, cG, cAbar, cBbar, gamma;
    nmod_poly_t modulus, modulus2;
#if WANT_ASSERT
    nmod_poly_t leadA, leadB;
#endif

#if WANT_ASSERT
    nmod_poly_init(leadA, ctx->ffinfo->mod.n);
    nmod_poly_init(leadB, ctx->ffinfo->mod.n);
    nmod_poly_set(leadA, nmod_mpolyun_leadcoeff_poly(A, ctx));
    nmod_poly_set(leadB, nmod_mpolyun_leadcoeff_poly(B, ctx));
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
    nmod_poly_gcd(gamma, nmod_mpolyun_leadcoeff_poly(A, ctx),
                         nmod_mpolyun_leadcoeff_poly(B, ctx));

    ldegA = nmod_mpolyun_lastdeg(A, ctx);
    ldegB = nmod_mpolyun_lastdeg(B, ctx);
    deggamma = nmod_poly_degree(gamma);
    bound = 1 + deggamma + FLINT_MAX(ldegA, ldegB);

    nmod_poly_init(Aeval, ctx->ffinfo->mod.n);
    nmod_poly_init(Beval, ctx->ffinfo->mod.n);
    nmod_poly_init(Geval, ctx->ffinfo->mod.n);
    nmod_poly_init(Abareval, ctx->ffinfo->mod.n);
    nmod_poly_init(Bbareval, ctx->ffinfo->mod.n);

    nmod_mpolyun_init(T, A->bits, ctx);

    nmod_poly_init(modulus, ctx->ffinfo->mod.n);
    nmod_poly_init(modulus2, ctx->ffinfo->mod.n);
    nmod_poly_one(modulus);

    if ((ctx->ffinfo->mod.n & UWORD(1)) == UWORD(0))
    {
        success = 0;
        goto cleanup;
    }

    alpha = (ctx->ffinfo->mod.n - UWORD(1))/UWORD(2);

choose_prime: /* prime is v - alpha */

    if (alpha < 2)
    {
        success = 0;
        goto cleanup;
    }

    alpha--;

    /* make sure evaluation point does not kill both lc(A) and lc(B) */
    gammaeval = nmod_poly_evaluate_nmod(gamma, alpha);
    if (gammaeval == 0)
    {
        goto choose_prime;
    }

    /* evaluation point should kill neither A nor B */
    nmod_mpolyun_eval_last_bivar(Aeval, A, alpha, ctx);
    nmod_mpolyun_eval_last_bivar(Beval, B, alpha, ctx);
    FLINT_ASSERT(Aeval->length > 0);
    FLINT_ASSERT(Beval->length > 0);

    nmod_poly_gcd(Geval, Aeval, Beval);
    nmod_poly_div(Abareval, Aeval, Geval);
    nmod_poly_div(Bbareval, Beval, Geval);

    FLINT_ASSERT(Geval->length > 0);
    FLINT_ASSERT(Abareval->length > 0);
    FLINT_ASSERT(Bbareval->length > 0);

    if (nmod_poly_degree(Geval) == 0)
    {
        nmod_mpolyun_one(G, ctx);
        nmod_mpolyun_swap(Abar, A);
        nmod_mpolyun_swap(Bbar, B);
        goto successful_put_content;    
    }

    if (nmod_poly_degree(modulus) > 0)
    {
        FLINT_ASSERT(G->length > 0);
        if (nmod_poly_degree(Geval) > G->exps[0])
        {
            goto choose_prime;
        }
        else if (nmod_poly_degree(Geval) < G->exps[0])
        {
            nmod_poly_one(modulus);
        }
    }

    nmod_poly_scalar_mul_nmod(Geval, Geval, gammaeval);

    if (nmod_poly_degree(modulus) > 0)
    {
        temp = nmod_poly_evaluate_nmod(modulus, alpha);
        temp = n_invmod(temp, ctx->ffinfo->mod.n);
        nmod_poly_scalar_mul_nmod(modulus, modulus, temp);
        nmod_mpolyun_addinterp_bivar(&ldegG, G, T, Geval, modulus, alpha, ctx);
        nmod_mpolyun_addinterp_bivar(&ldegAbar, Abar, T, Abareval, modulus, alpha, ctx);
        nmod_mpolyun_addinterp_bivar(&ldegBbar, Bbar, T, Bbareval, modulus, alpha, ctx);
    }
    else
    {
        nmod_mpolyun_startinterp_bivar(G, Geval, ctx);
        nmod_mpolyun_startinterp_bivar(Abar, Abareval, ctx);
        nmod_mpolyun_startinterp_bivar(Bbar, Bbareval, ctx);
        ldegG = 0;
        ldegAbar = 0;
        ldegBbar = 0;
    }

    nmod_poly_scalar_mul_nmod(modulus2, modulus, alpha);
    nmod_poly_shift_left(modulus, modulus, 1);
    nmod_poly_sub(modulus, modulus, modulus2);

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
    nmod_mpolyun_divexact_last(Abar, nmod_mpolyun_leadcoeff_poly(G, ctx), ctx);
    nmod_mpolyun_divexact_last(Bbar, nmod_mpolyun_leadcoeff_poly(G, ctx), ctx);

successful_put_content:

    nmod_mpolyun_mul_last(G, cG, ctx);
    nmod_mpolyun_mul_last(Abar, cAbar, ctx);
    nmod_mpolyun_mul_last(Bbar, cBbar, ctx);

    success = 1;

cleanup:

#if WANT_ASSERT
    if (success)
    {
        FLINT_ASSERT(1 == nmod_mpolyun_leadcoeff(G, ctx));
        nmod_poly_mul(modulus, nmod_mpolyun_leadcoeff_poly(G, ctx),
                               nmod_mpolyun_leadcoeff_poly(Abar, ctx));
        FLINT_ASSERT(nmod_poly_equal(modulus, leadA));
        nmod_poly_mul(modulus, nmod_mpolyun_leadcoeff_poly(G, ctx),
                               nmod_mpolyun_leadcoeff_poly(Bbar, ctx));
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

    nmod_poly_clear(Aeval);
    nmod_poly_clear(Beval);
    nmod_poly_clear(Geval);
    nmod_poly_clear(Abareval);
    nmod_poly_clear(Bbareval);

    nmod_mpolyun_clear(T, ctx);

    nmod_poly_clear(modulus);
    nmod_poly_clear(modulus2);

    return success;
}




int nmod_mpolyun_gcd_brown_smprime_bivar(nmod_mpolyun_t G,
  nmod_mpolyun_t Abar, nmod_mpolyun_t Bbar, nmod_mpolyun_t A, nmod_mpolyun_t B,
                                                    const nmod_mpoly_ctx_t ctx)
{
    int success;
    slong bound;
    mp_limb_t alpha, temp, gammaevalp, gammaevalm;
    nmod_poly_t Aevalp, Bevalp, Gevalp, Abarevalp, Bbarevalp;
    nmod_poly_t Aevalm, Bevalm, Gevalm, Abarevalm, Bbarevalm;
    nmod_mpolyun_t T;
    slong deggamma, ldegG, ldegAbar, ldegBbar, ldegA, ldegB;
    nmod_poly_t cA, cB, cG, cAbar, cBbar, gamma;
    nmod_poly_t modulus, modulus2, alphapow, r;
    int gstab, astab, bstab, use_stab;
#if WANT_ASSERT
    nmod_poly_t leadA, leadB;
#endif

#if WANT_ASSERT
    nmod_poly_init(leadA, ctx->ffinfo->mod.n);
    nmod_poly_init(leadB, ctx->ffinfo->mod.n);
    nmod_poly_set(leadA, nmod_mpolyun_leadcoeff_poly(A, ctx));
    nmod_poly_set(leadB, nmod_mpolyun_leadcoeff_poly(B, ctx));
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
    nmod_poly_gcd(gamma, nmod_mpolyun_leadcoeff_poly(A, ctx),
                         nmod_mpolyun_leadcoeff_poly(B, ctx));

    ldegA = nmod_mpolyun_lastdeg(A, ctx);
    ldegB = nmod_mpolyun_lastdeg(B, ctx);
    deggamma = nmod_poly_degree(gamma);
    bound = 1 + deggamma + FLINT_MAX(ldegA, ldegB);

    nmod_poly_init(Aevalp, ctx->ffinfo->mod.n);
    nmod_poly_init(Bevalp, ctx->ffinfo->mod.n);
    nmod_poly_init(Gevalp, ctx->ffinfo->mod.n);
    nmod_poly_init(Abarevalp, ctx->ffinfo->mod.n);
    nmod_poly_init(Bbarevalp, ctx->ffinfo->mod.n);
    nmod_poly_init(Aevalm, ctx->ffinfo->mod.n);
    nmod_poly_init(Bevalm, ctx->ffinfo->mod.n);
    nmod_poly_init(Gevalm, ctx->ffinfo->mod.n);
    nmod_poly_init(Abarevalm, ctx->ffinfo->mod.n);
    nmod_poly_init(Bbarevalm, ctx->ffinfo->mod.n);

    nmod_mpolyun_init(T, A->bits, ctx);

    nmod_poly_init(r, ctx->ffinfo->mod.n);
    nmod_poly_init(alphapow, ctx->ffinfo->mod.n);
    nmod_poly_fit_length(alphapow, FLINT_MAX(WORD(3), bound + 1));
    nmod_poly_init(modulus, ctx->ffinfo->mod.n);
    nmod_poly_init(modulus2, ctx->ffinfo->mod.n);
    nmod_poly_one(modulus);

    if ((ctx->ffinfo->mod.n & UWORD(1)) == UWORD(0))
    {
        success = 0;
        goto cleanup;
    }

    use_stab = 1;
    gstab = bstab = astab = 0;

    alpha = (ctx->ffinfo->mod.n - UWORD(1))/UWORD(2);

choose_prime: /* primes are v - alpha, v + alpha */

    if (alpha < 2)
    {
        success = 0;
        goto cleanup;
    }

    alpha--;

    FLINT_ASSERT(0 < alpha && alpha <= ctx->ffinfo->mod.n/2);
    FLINT_ASSERT(alphapow->alloc >= 2);
    alphapow->length = 2;
    alphapow->coeffs[0] = 1;
    alphapow->coeffs[1] = alpha;

    /* make sure evaluation point does not kill both lc(A) and lc(B) */
    _nmod_poly_eval2_pow(&gammaevalp, &gammaevalm, gamma, alphapow, ctx->ffinfo);
    if (gammaevalp == 0 || gammaevalm == 0)
    {
        goto choose_prime;
    }

    /* evaluation point should kill neither A nor B */
    nmod_mpolyun_eval2_last_bivar(Aevalp, Aevalm, A, alphapow, ctx);
    nmod_mpolyun_eval2_last_bivar(Bevalp, Bevalm, B, alphapow, ctx);
    FLINT_ASSERT(Aevalp->length > 0);
    FLINT_ASSERT(Aevalm->length > 0);
    FLINT_ASSERT(Bevalp->length > 0);
    FLINT_ASSERT(Bevalm->length > 0);

    if (use_stab && gstab)
    {
        slong Gdeg;
        nmod_mpolyun_eval2_last_bivar(Gevalp, Gevalm, G, alphapow, ctx);
        Gdeg = G->exps[0];
        success = 1;
        success = success && nmod_poly_degree(Gevalp) == Gdeg;
        success = success && nmod_poly_degree(Gevalm) == Gdeg;
        success = success && Gevalp->coeffs[Gdeg] == gammaevalp;
        success = success && Gevalm->coeffs[Gdeg] == gammaevalm;
        nmod_poly_divrem_basecase(Abarevalp, r, Aevalp, Gevalp);
        success = success && (r->length == 0);
        nmod_poly_divrem_basecase(Abarevalm, r, Aevalm, Gevalm);
        success = success && (r->length == 0);
        nmod_poly_divrem_basecase(Bbarevalp, r, Bevalp, Gevalp);
        success = success && (r->length == 0);
        nmod_poly_divrem_basecase(Bbarevalm, r, Bevalm, Gevalm);
        success = success && (r->length == 0);

        if (!success)
        {
            use_stab = 0;
            nmod_poly_one(modulus);
            alpha = (ctx->ffinfo->mod.n - UWORD(1))/UWORD(2);
            goto choose_prime;
        }

        nmod_poly_scalar_mul_nmod(Abarevalp, Abarevalp, gammaevalp);
        nmod_poly_scalar_mul_nmod(Abarevalm, Abarevalm, gammaevalm);
        nmod_poly_scalar_mul_nmod(Bbarevalp, Bbarevalp, gammaevalp);
        nmod_poly_scalar_mul_nmod(Bbarevalm, Bbarevalm, gammaevalm);
    }
    else
    {
        nmod_poly_gcd(Gevalp, Aevalp, Bevalp);
        nmod_poly_div(Abarevalp, Aevalp, Gevalp);
        nmod_poly_div(Bbarevalp, Bevalp, Gevalp);
        nmod_poly_gcd(Gevalm, Aevalm, Bevalm);
        nmod_poly_div(Abarevalm, Aevalm, Gevalm);
        nmod_poly_div(Bbarevalm, Bevalm, Gevalm);
        gstab = astab = bstab = 0;
    }

    FLINT_ASSERT(Gevalp->length > 0);
    FLINT_ASSERT(Abarevalp->length > 0);
    FLINT_ASSERT(Bbarevalp->length > 0);
    FLINT_ASSERT(Gevalm->length > 0);
    FLINT_ASSERT(Abarevalm->length > 0);
    FLINT_ASSERT(Bbarevalm->length > 0);

    if (nmod_poly_degree(Gevalp) == 0 || nmod_poly_degree(Gevalm) == 0)
    {
        nmod_mpolyun_one(G, ctx);
        nmod_mpolyun_swap(Abar, A);
        nmod_mpolyun_swap(Bbar, B);
        goto successful_put_content;    
    }

    if (nmod_poly_degree(Gevalp) != nmod_poly_degree(Gevalm))
    {
        goto choose_prime;
    }

    /* the Geval have matching degrees */
    if (nmod_poly_degree(modulus) > 0)
    {
        FLINT_ASSERT(G->length > 0);
        if (nmod_poly_degree(Gevalp) > G->exps[0])
        {
            goto choose_prime;
        }
        else if (nmod_poly_degree(Gevalp) < G->exps[0])
        {
            nmod_poly_one(modulus);
        }
    }

    /* update interpolants */
    nmod_poly_scalar_mul_nmod(Gevalp, Gevalp, gammaevalp);
    nmod_poly_scalar_mul_nmod(Gevalm, Gevalm, gammaevalm);
    if (nmod_poly_degree(modulus) > 0)
    {
        temp = nmod_poly_evaluate_nmod(modulus, alpha);
        FLINT_ASSERT(temp == nmod_poly_evaluate_nmod(modulus, ctx->ffinfo->mod.n - alpha));
        temp = nmod_mul(temp, alpha, ctx->ffinfo->mod);
        temp = nmod_add(temp, temp, ctx->ffinfo->mod);
        nmod_poly_scalar_mul_nmod(modulus, modulus, n_invmod(temp, ctx->ffinfo->mod.n));
        if (!gstab)
        {
            gstab = !nmod_mpolyun_addinterp2_bivar(&ldegG, G, T, Gevalp, Gevalm, modulus, alphapow, ctx);
        }
        nmod_mpolyun_addinterp2_bivar(&ldegAbar, Abar, T, Abarevalp, Abarevalm, modulus, alphapow, ctx);
        nmod_mpolyun_addinterp2_bivar(&ldegBbar, Bbar, T, Bbarevalp, Bbarevalm, modulus, alphapow, ctx);
    }
    else
    {
        nmod_mpolyun_startinterp2_bivar(&ldegG, G, Gevalp, Gevalm, alpha, ctx);
        nmod_mpolyun_startinterp2_bivar(&ldegAbar, Abar, Abarevalp, Abarevalm, alpha, ctx);
        nmod_mpolyun_startinterp2_bivar(&ldegBbar, Bbar, Bbarevalp, Bbarevalm, alpha, ctx);
        gstab = astab = bstab = 0;
    }
    temp = nmod_mul(alpha, alpha, ctx->ffinfo->mod);
    nmod_poly_scalar_mul_nmod(modulus2, modulus, temp);
    nmod_poly_shift_left(modulus, modulus, 2);
    nmod_poly_sub(modulus, modulus, modulus2);

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
    nmod_mpolyun_divexact_last(Abar, nmod_mpolyun_leadcoeff_poly(G, ctx), ctx);
    nmod_mpolyun_divexact_last(Bbar, nmod_mpolyun_leadcoeff_poly(G, ctx), ctx);

successful_put_content:

    nmod_mpolyun_mul_last(G, cG, ctx);
    nmod_mpolyun_mul_last(Abar, cAbar, ctx);
    nmod_mpolyun_mul_last(Bbar, cBbar, ctx);

    success = 1;

cleanup:

#if WANT_ASSERT
    if (success)
    {
        FLINT_ASSERT(1 == nmod_mpolyun_leadcoeff(G, ctx));
        nmod_poly_mul(modulus, nmod_mpolyun_leadcoeff_poly(G, ctx),
                               nmod_mpolyun_leadcoeff_poly(Abar, ctx));
        FLINT_ASSERT(nmod_poly_equal(modulus, leadA));
        nmod_poly_mul(modulus, nmod_mpolyun_leadcoeff_poly(G, ctx),
                               nmod_mpolyun_leadcoeff_poly(Bbar, ctx));
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

    nmod_poly_clear(Aevalp);
    nmod_poly_clear(Bevalp);
    nmod_poly_clear(Gevalp);
    nmod_poly_clear(Abarevalp);
    nmod_poly_clear(Bbarevalp);
    nmod_poly_clear(Aevalm);
    nmod_poly_clear(Bevalm);
    nmod_poly_clear(Gevalm);
    nmod_poly_clear(Abarevalm);
    nmod_poly_clear(Bbarevalm);

    nmod_mpolyun_clear(T, ctx);

    nmod_poly_clear(r);
    nmod_poly_clear(alphapow);
    nmod_poly_clear(modulus);
    nmod_poly_clear(modulus2);

    return success;
}



/*
    G, Abar, Bbar, A, B are all in
        Fp[X][x_0, ..., x_(var-1)][x_var],   var < ctx->minfo->nvars

    The inputs A and B might be modified.
    If the return is 1, then G is set to the gcd and Abar, Bbar the cofactors.
*/
int nmod_mpolyun_gcd_brown_smprime_ref(
    nmod_mpolyun_t G,
    nmod_mpolyun_t Abar,
    nmod_mpolyun_t Bbar,
    nmod_mpolyun_t A,
    nmod_mpolyun_t B,
    slong var,
    const nmod_mpoly_ctx_t ctx)
{
    int success;
    slong bound;
    slong offset, shift;
    mp_limb_t alpha, temp, gammaeval;
    nmod_mpolyun_t Aeval, Beval, Geval, Abareval, Bbareval;
    nmod_mpolyun_t T;
    slong deggamma, ldegG, ldegAbar, ldegBbar, ldegA, ldegB;
    nmod_poly_t cA, cB, cG, cAbar, cBbar, gamma;
    nmod_poly_t modulus, modulus2;
    mp_bitcnt_t bits = A->bits;
    slong N = mpoly_words_per_exp(bits, ctx->minfo);
#if WANT_ASSERT
    nmod_poly_t leadA, leadB;
#endif

    FLINT_ASSERT(var >= 0);
    if (var == 0)
    {
        /* bivariate is more comfortable separated */
        return nmod_mpolyun_gcd_brown_smprime_bivar_ref(G, Abar, Bbar, A, B, ctx);
    }

    mpoly_gen_offset_shift_sp(&offset, &shift, var - 1, G->bits, ctx->minfo);

#if WANT_ASSERT
    nmod_poly_init(leadA, ctx->ffinfo->mod.n);
    nmod_poly_init(leadB, ctx->ffinfo->mod.n);
    nmod_poly_set(leadA, nmod_mpolyun_leadcoeff_poly(A, ctx));
    nmod_poly_set(leadB, nmod_mpolyun_leadcoeff_poly(B, ctx));
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
    nmod_poly_gcd(gamma, nmod_mpolyun_leadcoeff_poly(A, ctx),
                         nmod_mpolyun_leadcoeff_poly(B, ctx));

    ldegA = nmod_mpolyun_lastdeg(A, ctx);
    ldegB = nmod_mpolyun_lastdeg(B, ctx);
    deggamma = nmod_poly_degree(gamma);
    bound = 1 + deggamma + FLINT_MAX(ldegA, ldegB);

    nmod_mpolyun_init(T, bits, ctx);
    nmod_poly_init(modulus, ctx->ffinfo->mod.n);
    nmod_poly_init(modulus2, ctx->ffinfo->mod.n);
    nmod_poly_one(modulus);

    nmod_mpolyun_init(Aeval, bits, ctx);
    nmod_mpolyun_init(Beval, bits, ctx);
    nmod_mpolyun_init(Geval, bits, ctx);
    nmod_mpolyun_init(Abareval, bits, ctx);
    nmod_mpolyun_init(Bbareval, bits, ctx);

    if ((ctx->ffinfo->mod.n & UWORD(1)) == UWORD(0))
    {
        success = 0;
        goto cleanup;
    }

    alpha = (ctx->ffinfo->mod.n - UWORD(1))/UWORD(2);

choose_prime: /* prime is v - alpha */

    if (alpha < 2)
    {
        success = 0;
        goto cleanup;
    }

    alpha--;

    /* make sure evaluation point does not kill both lc(A) and lc(B) */
    gammaeval = nmod_poly_evaluate_nmod(gamma, alpha);
    if (gammaeval == 0)
    {
        goto choose_prime;
    }

    /* evaluation point should kill neither A bor B */
    nmod_mpolyun_eval_last_un(Aeval, A, var, alpha, ctx);
    nmod_mpolyun_eval_last_un(Beval, B, var, alpha, ctx);
    FLINT_ASSERT(Aeval->length > 0);
    FLINT_ASSERT(Beval->length > 0);

    success = nmod_mpolyun_gcd_brown_smprime_ref(Geval, Abareval, Bbareval,
                                                   Aeval, Beval, var - 1, ctx);
    if (success == 0)
    {
        goto choose_prime;
    }

    FLINT_ASSERT(Geval->length > 0);
    FLINT_ASSERT(Abareval->length);
    FLINT_ASSERT(Bbareval->length > 0);

    if (nmod_mpolyun_is_nonzero_nmod(Geval, ctx))
    {
        nmod_mpolyun_one(G, ctx);
        nmod_mpolyun_swap(Abar, A);
        nmod_mpolyun_swap(Bbar, B);
        goto successful_put_content;    
    }

    if (nmod_poly_degree(modulus) > 0)
    {
        int cmp = 0;
        FLINT_ASSERT(G->length > 0);
        if (G->exps[0] != Geval->exps[0])
        {
            cmp = G->exps[0] > Geval->exps[0] ? 1 : -1;
        }
        if (cmp == 0)
        {
            slong k = nmod_poly_degree((Geval->coeffs + 0)->coeffs + 0);
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

    temp = nmod_mpolyn_leadcoeff(Geval->coeffs + 0, ctx);
    temp = n_invmod(temp, ctx->ffinfo->mod.n);
    temp = nmod_mul(gammaeval, temp, ctx->ffinfo->mod);
    nmod_mpolyun_scalar_mul_nmod(Geval, temp, ctx);

    if (nmod_poly_degree(modulus) > 0)
    {
        temp = nmod_poly_evaluate_nmod(modulus, alpha);
        temp = n_invmod(temp, ctx->ffinfo->mod.n);
        nmod_poly_scalar_mul_nmod(modulus, modulus, temp);
        nmod_mpolyun_addinterp_un(&ldegG, G, T, Geval, var, modulus, alpha, ctx);
        nmod_mpolyun_addinterp_un(&ldegAbar, Abar, T, Abareval, var, modulus, alpha, ctx);
        nmod_mpolyun_addinterp_un(&ldegBbar, Bbar, T, Bbareval, var, modulus, alpha, ctx);
    }
    else
    {
        nmod_mpolyun_set_popup(G, Geval, var, ctx);
        nmod_mpolyun_set_popup(Abar, Abareval, var, ctx);
        nmod_mpolyun_set_popup(Bbar, Bbareval, var, ctx);
        ldegG = 0;
        ldegAbar = 0;
        ldegBbar = 0;
    }

    nmod_poly_scalar_mul_nmod(modulus2, modulus, alpha);
    nmod_poly_shift_left(modulus, modulus, 1);
    nmod_poly_sub(modulus, modulus, modulus2);

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
    nmod_mpolyun_divexact_last(Abar, nmod_mpolyun_leadcoeff_poly(G, ctx), ctx);
    nmod_mpolyun_divexact_last(Bbar, nmod_mpolyun_leadcoeff_poly(G, ctx), ctx);

successful_put_content:

    nmod_mpolyun_mul_last(G, cG, ctx);
    nmod_mpolyun_mul_last(Abar, cAbar, ctx);
    nmod_mpolyun_mul_last(Bbar, cBbar, ctx);

    success = 1;

cleanup:

#if WANT_ASSERT
    if (success)
    {
        FLINT_ASSERT(1 == nmod_mpolyun_leadcoeff(G, ctx));
        nmod_poly_mul(modulus, nmod_mpolyun_leadcoeff_poly(G, ctx),
                               nmod_mpolyun_leadcoeff_poly(Abar, ctx));
        FLINT_ASSERT(nmod_poly_equal(modulus, leadA));
        nmod_poly_mul(modulus, nmod_mpolyun_leadcoeff_poly(G, ctx),
                               nmod_mpolyun_leadcoeff_poly(Bbar, ctx));
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

    nmod_mpolyun_clear(Aeval, ctx);
    nmod_mpolyun_clear(Beval, ctx);
    nmod_mpolyun_clear(Geval, ctx);
    nmod_mpolyun_clear(Abareval, ctx);
    nmod_mpolyun_clear(Bbareval, ctx);

    nmod_mpolyun_clear(T, ctx);

    nmod_poly_clear(modulus);
    nmod_poly_clear(modulus2);

    return success;
}

/* optimized version of nmod_mpolyun_gcd_brown_smprime_ref */
int nmod_mpolyun_gcd_brown_smprime(
    nmod_mpolyun_t G,
    nmod_mpolyun_t Abar,
    nmod_mpolyun_t Bbar,
    nmod_mpolyun_t A,
    nmod_mpolyun_t B,
    slong var,
    const nmod_mpoly_ctx_t ctx)
{
    int success;
    slong bound;
    slong offset, shift;
    mp_limb_t alpha, temp, gammaevalp, gammaevalm;
    nmod_mpolyun_t Aevalp, Bevalp, Gevalp, Abarevalp, Bbarevalp;
    nmod_mpolyun_t Aevalm, Bevalm, Gevalm, Abarevalm, Bbarevalm;
    nmod_mpolyun_t T;
    slong deggamma, ldegG, ldegAbar, ldegBbar, ldegA, ldegB;
    nmod_poly_t cA, cB, cG, cAbar, cBbar, gamma;
    nmod_poly_t modulus, modulus2, alphapow;
    mp_bitcnt_t bits = A->bits;
    slong N = mpoly_words_per_exp(bits, ctx->minfo);
#if WANT_ASSERT
    nmod_poly_t leadA, leadB;
#endif

    if (var == 0)
    {
        /* bivariate is more comfortable separated */
        return nmod_mpolyun_gcd_brown_smprime_bivar(G, Abar, Bbar, A, B, ctx);
    }

    mpoly_gen_offset_shift_sp(&offset, &shift, var - 1, G->bits, ctx->minfo);

#if WANT_ASSERT
    nmod_poly_init(leadA, ctx->ffinfo->mod.n);
    nmod_poly_init(leadB, ctx->ffinfo->mod.n);
    nmod_poly_set(leadA, nmod_mpolyun_leadcoeff_poly(A, ctx));
    nmod_poly_set(leadB, nmod_mpolyun_leadcoeff_poly(B, ctx));
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
    nmod_poly_gcd(gamma, nmod_mpolyun_leadcoeff_poly(A, ctx),
                         nmod_mpolyun_leadcoeff_poly(B, ctx));

    ldegA = nmod_mpolyun_lastdeg(A, ctx);
    ldegB = nmod_mpolyun_lastdeg(B, ctx);
    deggamma = nmod_poly_degree(gamma);
    bound = 1 + deggamma + FLINT_MAX(ldegA, ldegB);

    nmod_mpolyun_init(Aevalp, bits, ctx);
    nmod_mpolyun_init(Bevalp, bits, ctx);
    nmod_mpolyun_init(Gevalp, bits, ctx);
    nmod_mpolyun_init(Abarevalp, bits, ctx);
    nmod_mpolyun_init(Bbarevalp, bits, ctx);
    nmod_mpolyun_init(Aevalm, bits, ctx);
    nmod_mpolyun_init(Bevalm, bits, ctx);
    nmod_mpolyun_init(Gevalm, bits, ctx);
    nmod_mpolyun_init(Abarevalm, bits, ctx);
    nmod_mpolyun_init(Bbarevalm, bits, ctx);
    nmod_mpolyun_init(T, bits, ctx);

    nmod_poly_init(alphapow, ctx->ffinfo->mod.n);
    nmod_poly_fit_length(alphapow, FLINT_MAX(WORD(3), bound + 1));

    nmod_poly_init(modulus, ctx->ffinfo->mod.n);
    nmod_poly_init(modulus2, ctx->ffinfo->mod.n);
    nmod_poly_one(modulus);

    if ((ctx->ffinfo->mod.n & UWORD(1)) == UWORD(0))
    {
        success = 0;
        goto cleanup;
    }

    alpha = (ctx->ffinfo->mod.n - UWORD(1))/UWORD(2);

choose_prime:

    if (alpha < 2)
    {
        success = 0;
        goto cleanup;
    }

    alpha--;

    FLINT_ASSERT(0 < alpha && alpha <= ctx->ffinfo->mod.n/2);
    FLINT_ASSERT(alphapow->alloc >= 2);
    alphapow->length = 2;
    alphapow->coeffs[0] = 1;
    alphapow->coeffs[1] = alpha;

    /* make sure evaluation point does not kill both lc(A) and lc(B) */
    _nmod_poly_eval2_pow(&gammaevalp, &gammaevalm, gamma, alphapow, ctx->ffinfo);
    if (gammaevalp == 0 || gammaevalm == 0)
    {
        goto choose_prime;
    }

    /* evaluation point should kill neither A nor B */
    nmod_mpolyun_eval2_last_un(Aevalp, Aevalm, A, var, alphapow, ctx);
    nmod_mpolyun_eval2_last_un(Bevalp, Bevalm, B, var, alphapow, ctx);
    FLINT_ASSERT(Aevalp->length > 0);
    FLINT_ASSERT(Aevalm->length > 0);
    FLINT_ASSERT(Bevalp->length > 0);
    FLINT_ASSERT(Bevalm->length > 0);

    success = nmod_mpolyun_gcd_brown_smprime(Gevalp, Abarevalp, Bbarevalp,
                                               Aevalp, Bevalp, var - 1, ctx);
    success = success && nmod_mpolyun_gcd_brown_smprime(Gevalm, Abarevalm, Bbarevalm,
                                               Aevalm, Bevalm, var - 1, ctx);
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

    if (   nmod_mpolyun_is_nonzero_nmod(Gevalp, ctx)
        || nmod_mpolyun_is_nonzero_nmod(Gevalm, ctx))
    {
        nmod_mpolyun_one(G, ctx);
        nmod_mpolyun_swap(Abar, A);
        nmod_mpolyun_swap(Bbar, B);
        goto successful_put_content;    
    }

    if (Gevalp->exps[0] != Gevalm->exps[0])
    {
        goto choose_prime;
    }
    if (   nmod_poly_degree((Gevalp->coeffs + 0)->coeffs + 0)
        != nmod_poly_degree((Gevalm->coeffs + 0)->coeffs + 0))
    {
        goto choose_prime;
    }
    if (!mpoly_monomial_equal((Gevalp->coeffs + 0)->exps + N*0,
                              (Gevalm->coeffs + 0)->exps + N*0, N))
    {
        goto choose_prime;
    }

    /* the Geval have matching degrees */
    if (nmod_poly_degree(modulus) > 0)
    {
        int cmp = 0;
        FLINT_ASSERT(G->length > 0);
        if (G->exps[0] != Gevalp->exps[0])
        {
            cmp = G->exps[0] > Gevalp->exps[0] ? 1 : -1;
        }
        if (cmp == 0)
        {
            slong k = nmod_poly_degree((Gevalp->coeffs + 0)->coeffs + 0);
            cmp = mpoly_monomial_cmp_nomask_extra(
                   (G->coeffs + 0)->exps + N*0,
                   (Gevalp->coeffs + 0)->exps + N*0, N, offset, k << shift);
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

    /* update interpolants */
    temp = nmod_mpolyn_leadcoeff(Gevalp->coeffs + 0, ctx);
    temp = n_invmod(temp, ctx->ffinfo->mod.n);
    temp = nmod_mul(gammaevalp, temp, ctx->ffinfo->mod);
    nmod_mpolyun_scalar_mul_nmod(Gevalp, temp, ctx);
    temp = nmod_mpolyn_leadcoeff(Gevalm->coeffs + 0, ctx);
    temp = n_invmod(temp, ctx->ffinfo->mod.n);
    temp = nmod_mul(gammaevalm, temp, ctx->ffinfo->mod);
    nmod_mpolyun_scalar_mul_nmod(Gevalm, temp, ctx);
    if (nmod_poly_degree(modulus) > 0)
    {
        temp = nmod_poly_evaluate_nmod(modulus, alpha);
        FLINT_ASSERT(temp == nmod_poly_evaluate_nmod(modulus, ctx->ffinfo->mod.n - alpha));
        temp = nmod_mul(temp, alpha, ctx->ffinfo->mod);
        temp = nmod_add(temp, temp, ctx->ffinfo->mod);
        nmod_poly_scalar_mul_nmod(modulus, modulus, n_invmod(temp, ctx->ffinfo->mod.n));
        nmod_mpolyun_addinterp2_un(&ldegG, G, T, Gevalp, Gevalm, var, modulus, alphapow, ctx);
        nmod_mpolyun_addinterp2_un(&ldegAbar, Abar, T, Abarevalp, Abarevalm, var, modulus, alphapow, ctx);
        nmod_mpolyun_addinterp2_un(&ldegBbar, Bbar, T, Bbarevalp, Bbarevalm, var, modulus, alphapow, ctx);
    }
    else
    {
        nmod_mpolyun_startinterp2_un(&ldegG, G, Gevalp, Gevalm, var, alpha, ctx);
        nmod_mpolyun_startinterp2_un(&ldegAbar, Abar, Abarevalp, Abarevalm, var, alpha, ctx);
        nmod_mpolyun_startinterp2_un(&ldegBbar, Bbar, Bbarevalp, Bbarevalm, var, alpha, ctx);
    }
    temp = nmod_mul(alpha, alpha, ctx->ffinfo->mod);
    nmod_poly_scalar_mul_nmod(modulus2, modulus, temp);
    nmod_poly_shift_left(modulus, modulus, 2);
    nmod_poly_sub(modulus, modulus, modulus2);

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
    nmod_mpolyun_divexact_last(Abar, nmod_mpolyun_leadcoeff_poly(G, ctx), ctx);
    nmod_mpolyun_divexact_last(Bbar, nmod_mpolyun_leadcoeff_poly(G, ctx), ctx);

successful_put_content:

    nmod_mpolyun_mul_last(G, cG, ctx);
    nmod_mpolyun_mul_last(Abar, cAbar, ctx);
    nmod_mpolyun_mul_last(Bbar, cBbar, ctx);

    success = 1;

cleanup:

#if WANT_ASSERT
    if (success)
    {
        FLINT_ASSERT(1 == nmod_mpolyun_leadcoeff(G, ctx));
        nmod_poly_mul(modulus, nmod_mpolyun_leadcoeff_poly(G, ctx),
                               nmod_mpolyun_leadcoeff_poly(Abar, ctx));
        FLINT_ASSERT(nmod_poly_equal(modulus, leadA));
        nmod_poly_mul(modulus, nmod_mpolyun_leadcoeff_poly(G, ctx),
                               nmod_mpolyun_leadcoeff_poly(Bbar, ctx));
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

    nmod_mpolyun_clear(Aevalp, ctx);
    nmod_mpolyun_clear(Bevalp, ctx);
    nmod_mpolyun_clear(Gevalp, ctx);
    nmod_mpolyun_clear(Abarevalp, ctx);
    nmod_mpolyun_clear(Bbarevalp, ctx);
    nmod_mpolyun_clear(Aevalm, ctx);
    nmod_mpolyun_clear(Bevalm, ctx);
    nmod_mpolyun_clear(Gevalm, ctx);
    nmod_mpolyun_clear(Abarevalm, ctx);
    nmod_mpolyun_clear(Bbarevalm, ctx);
    nmod_mpolyun_clear(T, ctx);

    nmod_poly_clear(alphapow);
    nmod_poly_clear(modulus);
    nmod_poly_clear(modulus2);

    return success;
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
    nmod_poly_set(leadA, nmod_mpolyun_leadcoeff_poly(A, ctx));
    nmod_poly_set(leadB, nmod_mpolyun_leadcoeff_poly(B, ctx));
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
    nmod_poly_gcd(gamma, nmod_mpolyun_leadcoeff_poly(A, ctx),
                         nmod_mpolyun_leadcoeff_poly(B, ctx));

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
    nmod_mpolyun_divexact_last(Abar, nmod_mpolyun_leadcoeff_poly(G, ctx), ctx);
    nmod_mpolyun_divexact_last(Bbar, nmod_mpolyun_leadcoeff_poly(G, ctx), ctx);

successful_put_content:

    nmod_mpolyun_mul_last(G, cG, ctx);
    nmod_mpolyun_mul_last(Abar, cAbar, ctx);
    nmod_mpolyun_mul_last(Bbar, cBbar, ctx);

    success = 1;

cleanup:

#if WANT_ASSERT
    if (success)
    {
        nmod_poly_mul(modulus, nmod_mpolyun_leadcoeff_poly(G, ctx),
                               nmod_mpolyun_leadcoeff_poly(Abar, ctx));
        FLINT_ASSERT(nmod_poly_equal(modulus, leadA));
        nmod_poly_mul(modulus, nmod_mpolyun_leadcoeff_poly(G, ctx),
                               nmod_mpolyun_leadcoeff_poly(Bbar, ctx));
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


int nmod_mpolyun_gcd_brown_lgprime(
    nmod_mpolyun_t G,
    nmod_mpolyun_t Abar,
    nmod_mpolyun_t Bbar,
    nmod_mpolyun_t A,
    nmod_mpolyun_t B,
    slong var,
    const nmod_mpoly_ctx_t ctx)
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
    nmod_poly_set(leadA, nmod_mpolyun_leadcoeff_poly(A, ctx));
    nmod_poly_set(leadB, nmod_mpolyun_leadcoeff_poly(B, ctx));
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
    nmod_poly_gcd(gamma, nmod_mpolyun_leadcoeff_poly(A, ctx),
                         nmod_mpolyun_leadcoeff_poly(B, ctx));

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

    fq_nmod_inv(temp, fq_nmod_mpolyn_leadcoeff(Geval->coeffs + 0, ectx), ectx->fqctx);
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
    nmod_mpolyun_divexact_last(Abar, nmod_mpolyun_leadcoeff_poly(G, ctx), ctx);
    nmod_mpolyun_divexact_last(Bbar, nmod_mpolyun_leadcoeff_poly(G, ctx), ctx);

successful_put_content:

    nmod_mpolyun_mul_last(G, cG, ctx);
    nmod_mpolyun_mul_last(Abar, cAbar, ctx);
    nmod_mpolyun_mul_last(Bbar, cBbar, ctx);

    success = 1;

cleanup:

#if WANT_ASSERT
    if (success)
    {
        FLINT_ASSERT(1 == nmod_mpolyun_leadcoeff(G, ctx));
        nmod_poly_mul(modulus, nmod_mpolyun_leadcoeff_poly(G, ctx),
                               nmod_mpolyun_leadcoeff_poly(Abar, ctx));
        FLINT_ASSERT(nmod_poly_equal(modulus, leadA));
        nmod_poly_mul(modulus, nmod_mpolyun_leadcoeff_poly(G, ctx),
                               nmod_mpolyun_leadcoeff_poly(Bbar, ctx));
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
