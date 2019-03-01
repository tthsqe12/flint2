/*
    Copyright (C) 2019 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fq_nmod_mpoly.h"

int fq_nmod_mpolyu_gcd_brown_lgprime(fq_nmod_mpolyu_t G,
          fq_nmod_mpolyu_t Abar, fq_nmod_mpolyu_t Bbar, fq_nmod_mpolyu_t A,
                  fq_nmod_mpolyu_t B, slong var, const fq_nmod_mpoly_ctx_t ctx)
{
    int success;
    slong bound;
    fq_nmod_t temp, gammaeval;
    fq_nmod_mpolyu_t Aeval, Beval, Geval, Abareval, Bbareval;
    fq_nmod_mpolyun_t An, Bn, Gs, Abars, Bbars, T;
    slong deggamma, ldegGs, ldegAbars, ldegBbars, ldegA, ldegB;
    fq_nmod_poly_t cA, cB, cG, cAbar, cBbar, gamma, trem;
    fq_nmod_poly_t cGs, cAbars, cBbars, modulus;
    mp_bitcnt_t bits = A->bits;
    slong N = mpoly_words_per_exp(bits, ctx->minfo);
    flint_rand_t randstate;
    _fq_nmod_mpoly_embed_chooser_t embc;
    _fq_nmod_embed_struct * cur_emb;
    fq_nmod_mpoly_ctx_t ectx;

    if (var == -WORD(1))
    {
        /* no more variables left to interpolate */
        return fq_nmod_mpolyu_gcd_brown_univar(G, Abar, Bbar, A, B, ctx);
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

    fq_nmod_poly_init(modulus, ctx->fqctx);
    fq_nmod_poly_one(modulus, ctx->fqctx);

    flint_randinit(randstate);
    cur_emb = _fq_nmod_mpoly_embed_chooser_init(embc, ectx, ctx, randstate);

    /*
        Once Aeval, Beval, ..., t are inited in ectx->fqctx, they do not need
        to be cleared and reinited when ectx->fqctx changes.
    */
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

    cur_emb = _fq_nmod_mpoly_embed_chooser_next(embc, ectx, ctx, randstate);
    if (cur_emb == NULL)
    {
        /* ran out of primes */
        success = 0;
        goto cleanup;
    }

have_prime:

    /* make sure reduction does not kill both lc */
    _fq_nmod_embed_sm_to_lg(gammaeval, gamma, cur_emb);
    if (fq_nmod_is_zero(gammaeval, ectx->fqctx))
    {
        goto choose_prime;
    }

    /* make sure reduction does not kill either A or B */
    fq_nmod_mpolyun_redto_fq_nmod_mpolyu(Aeval, An, ectx, ctx, cur_emb);
    fq_nmod_mpolyun_redto_fq_nmod_mpolyu(Beval, Bn, ectx, ctx, cur_emb);
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
        fq_nmod_mpolyun_one(Gs, ctx);
        fq_nmod_mpolyun_set(Abars, An, ctx);
        fq_nmod_mpolyun_set(Bbars, Bn, ctx);
        goto successful;    
    }

    if (fq_nmod_poly_degree(modulus, ctx->fqctx) > 0)
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
            fq_nmod_poly_one(modulus, ctx->fqctx);
        }
    }

    fq_nmod_inv(temp, fq_nmod_mpolyu_leadcoeff_ref(Geval, ectx), ectx->fqctx);
    fq_nmod_mul(temp, temp, gammaeval, ectx->fqctx);
    fq_nmod_mpolyu_scalar_mul_fq_nmod(Geval, temp, ectx);

    if (fq_nmod_poly_degree(modulus, ctx->fqctx) > 0)
    {
        fq_nmod_mpolyun_addinterp_lgprime(&ldegGs, Gs, T, modulus, ctx, Geval, ectx, cur_emb);
        fq_nmod_mpolyun_addinterp_lgprime(&ldegAbars, Abars, T, modulus, ctx, Abareval, ectx, cur_emb);
        fq_nmod_mpolyun_addinterp_lgprime(&ldegBbars, Bbars, T, modulus, ctx, Bbareval, ectx, cur_emb);
    }
    else
    {
        fq_nmod_mpolyun_startinterp_lgprime(&ldegGs, Gs, ctx, Geval, ectx, cur_emb);
        fq_nmod_mpolyun_startinterp_lgprime(&ldegAbars, Abars, ctx, Abareval, ectx, cur_emb);
        fq_nmod_mpolyun_startinterp_lgprime(&ldegBbars, Bbars, ctx, Bbareval, ectx, cur_emb);
    }
    fq_nmod_poly_mul(modulus, modulus, cur_emb->h, ctx->fqctx);

    if (fq_nmod_poly_degree(modulus, ctx->fqctx) < bound)
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

    fq_nmod_poly_one(modulus, ctx->fqctx);
    goto choose_prime;

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

    fq_nmod_mpolyun_clear(Gs, ctx);
    fq_nmod_mpolyun_clear(Abars, ctx);
    fq_nmod_mpolyun_clear(Bbars, ctx);
    fq_nmod_mpolyun_clear(T, ctx);

    fq_nmod_poly_clear(modulus, ctx->fqctx);

    fq_nmod_mpolyu_clear(Aeval, ectx);
    fq_nmod_mpolyu_clear(Beval, ectx);
    fq_nmod_mpolyu_clear(Geval, ectx);
    fq_nmod_mpolyu_clear(Abareval, ectx);
    fq_nmod_mpolyu_clear(Bbareval, ectx);
    fq_nmod_clear(gammaeval, ectx->fqctx);
    fq_nmod_clear(temp, ectx->fqctx);

    _fq_nmod_mpoly_embed_chooser_clear(embc, ectx, ctx, randstate);

    flint_randclear(randstate);

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
