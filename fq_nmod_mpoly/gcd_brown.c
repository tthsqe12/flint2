/*
    Copyright (C) 2019 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fq_nmod_mpoly.h"

/*
    Try to set G to the gcd of A and B using Brown's alogrithm M.
    This function switches to a big primes version if needed.
    It should only really fail if the dense size of the inputs is too large.
*/
int fq_nmod_mpoly_gcd_brown(fq_nmod_mpoly_t G,
                            const fq_nmod_mpoly_t A, const fq_nmod_mpoly_t B,
                                                 const fq_nmod_mpoly_ctx_t ctx)
{
    int success;
    fq_nmod_mpolyd_t Ad, Bd, Gd, Abar, Bbar;
    fq_nmod_mpolyd_ctx_t dctx;
    slong nvars = ctx->minfo->nvars;

    success = 1;

    if (fq_nmod_mpoly_is_zero(A, ctx))
    {
        fq_nmod_mpoly_set(G, B, ctx);
        goto cleanup_stage0;
    }
    if (fq_nmod_mpoly_is_zero(B, ctx))
    {
        fq_nmod_mpoly_set(G, A, ctx);
        goto cleanup_stage0;
    }

    fq_nmod_mpolyd_ctx_init2(dctx, nvars, ctx->fqctx);
    success = fq_nmod_mpolyd_ctx_set_for_gcd(dctx, A, B, ctx);
    if (!success)
    {
        fq_nmod_mpoly_zero(G, ctx);
        goto cleanup_stage1;
    }

    fq_nmod_mpolyd_init(Ad, nvars, ctx->fqctx);
    fq_nmod_mpolyd_init(Bd, nvars, ctx->fqctx);
    fq_nmod_mpolyd_init(Gd, nvars, ctx->fqctx);
    fq_nmod_mpolyd_init(Abar, nvars, ctx->fqctx);
    fq_nmod_mpolyd_init(Bbar, nvars, ctx->fqctx);

    fq_nmod_mpoly_convert_to_fq_nmod_mpolyd(Ad, dctx, A, ctx);
    fq_nmod_mpoly_convert_to_fq_nmod_mpolyd(Bd, dctx, B, ctx);
    success = fq_nmod_mpolyd_gcd_brown_smprime(Gd, Abar, Bbar, Ad, Bd, dctx);
    if (!success)
    {
        fq_nmod_mpoly_convert_to_fq_nmod_mpolyd(Ad, dctx, A, ctx);
        fq_nmod_mpoly_convert_to_fq_nmod_mpolyd(Bd, dctx, B, ctx);
        success = fq_nmod_mpolyd_gcd_brown_lgprime(Gd, Abar, Bbar, Ad, Bd, dctx);
    }
    if (success)
    {
        fq_nmod_mpoly_convert_from_fq_nmod_mpolyd(G, ctx, Gd, dctx);
    }

    fq_nmod_mpolyd_clear(Bbar, ctx->fqctx);
    fq_nmod_mpolyd_clear(Abar, ctx->fqctx);
    fq_nmod_mpolyd_clear(Gd, ctx->fqctx);
    fq_nmod_mpolyd_clear(Bd, ctx->fqctx);
    fq_nmod_mpolyd_clear(Ad, ctx->fqctx);

cleanup_stage1:

    fq_nmod_mpolyd_ctx_clear(dctx);

cleanup_stage0:

    if (success && !fq_nmod_mpoly_is_zero(G, ctx))
    {
        fq_nmod_mpoly_make_monic(G, G, ctx);
    }

    return success;
}




void fq_nmod_mpoly_to_fq_nmod_poly_keepbits(fq_nmod_poly_t A, slong * Ashift,
            const fq_nmod_mpoly_t B, slong var, const fq_nmod_mpoly_ctx_t ctx);

void fq_nmod_mpoly_from_fq_nmod_poly_keepbits(fq_nmod_mpoly_t A,
         const fq_nmod_poly_t B, slong Bshift, slong var, mp_bitcnt_t bits,
                                                const fq_nmod_mpoly_ctx_t ctx);


int fq_nmod_mpoly_gcd_brownnew(fq_nmod_mpoly_t G, const fq_nmod_mpoly_t A,
                        const fq_nmod_mpoly_t B, const fq_nmod_mpoly_ctx_t ctx)
{
    int success;
    slong * perm;
    slong i;
    mp_bitcnt_t new_bits;
    fq_nmod_mpoly_ctx_t uctx;
    fq_nmod_mpolyu_t Au, Bu, Gu, Abaru, Bbaru;

    if (fq_nmod_mpoly_is_zero(A, ctx))
    {
        if (fq_nmod_mpoly_is_zero(B, ctx))
        {
            fq_nmod_mpoly_zero(G, ctx);
        }
        else
        {
            fq_nmod_mpoly_make_monic(G, B, ctx);
        }
        return 1;
    }

    if (fq_nmod_mpoly_is_zero(B, ctx))
    {
        fq_nmod_mpoly_make_monic(G, A, ctx);
        return 1;
    }

    if (A->bits > FLINT_BITS || B->bits > FLINT_BITS)
    {
        return 0;
    }

    if (ctx->minfo->nvars == 1)
    {
        slong shiftA, shiftB;
        fq_nmod_poly_t a, b, g;
        fq_nmod_poly_init(a, ctx->fqctx);
        fq_nmod_poly_init(b, ctx->fqctx);
        fq_nmod_poly_init(g, ctx->fqctx);
        fq_nmod_mpoly_to_fq_nmod_poly_keepbits(a, &shiftA, A, 0, ctx);
        fq_nmod_mpoly_to_fq_nmod_poly_keepbits(b, &shiftB, B, 0, ctx);
        fq_nmod_poly_gcd(g, a, b, ctx->fqctx);
        fq_nmod_mpoly_from_fq_nmod_poly_keepbits(G, g, FLINT_MIN(shiftA, shiftB), 0, A->bits, ctx);
        fq_nmod_poly_clear(a, ctx->fqctx);
        fq_nmod_poly_clear(b, ctx->fqctx);
        fq_nmod_poly_clear(g, ctx->fqctx);
        return 1;
    }

    perm = (slong *) flint_malloc(ctx->minfo->nvars*sizeof(slong));
    for (i = 0; i < ctx->minfo->nvars; i++)
    {
        perm[i] = i;
    }

    new_bits = FLINT_MAX(A->bits, B->bits);

    fq_nmod_mpoly_ctx_init(uctx, ctx->minfo->nvars - 1, ORD_LEX, ctx->fqctx);
    fq_nmod_mpolyu_init(Au, new_bits, uctx);
    fq_nmod_mpolyu_init(Bu, new_bits, uctx);
    fq_nmod_mpolyu_init(Gu, new_bits, uctx);
    fq_nmod_mpolyu_init(Abaru, new_bits, uctx);
    fq_nmod_mpolyu_init(Bbaru, new_bits, uctx);

    fq_nmod_mpoly_to_mpolyu_perm(Au, A, perm, uctx, ctx);
    fq_nmod_mpoly_to_mpolyu_perm(Bu, B, perm, uctx, ctx);

/*
printf("gcd_brownnew:\n");
printf("Au: "); fq_nmod_mpolyu_print_pretty(Au, NULL, uctx); printf("\n");
printf("Bu: "); fq_nmod_mpolyu_print_pretty(Bu, NULL, uctx); printf("\n");
*/
    success = fq_nmod_mpolyu_gcd_brown_smprime(Gu, Abaru, Bbaru, Au, Bu, uctx->minfo->nvars - 1, uctx);
    if (!success)
    {
        success = fq_nmod_mpolyu_gcd_brown_lgprime(Gu, Abaru, Bbaru, Au, Bu, uctx->minfo->nvars - 1, uctx);
    }
    if (success)
    {
        fq_nmod_mpoly_from_mpolyu_perm(G, Gu, 1, perm, uctx, ctx);
        fq_nmod_mpoly_make_monic(G, G, ctx);
        success = 1;
    }

    fq_nmod_mpolyu_clear(Au, uctx);
    fq_nmod_mpolyu_clear(Bu, uctx);
    fq_nmod_mpolyu_clear(Gu, uctx);
    fq_nmod_mpolyu_clear(Abaru, uctx);
    fq_nmod_mpolyu_clear(Bbaru, uctx);
    fq_nmod_mpoly_ctx_clear(uctx);

    flint_free(perm);

    return success;
}
