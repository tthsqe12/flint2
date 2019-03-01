/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "nmod_mpoly.h"

/*
    Try to set G to the gcd of A and B using Brown's alogrithm M.
    This function switches to a big primes version if needed.
    It should only really fail if the dense size of the inputs is too large.
*/
int nmod_mpoly_gcd_brown(nmod_mpoly_t G, const nmod_mpoly_t A,
                              const nmod_mpoly_t B, const nmod_mpoly_ctx_t ctx)
{
    int success;
    nmod_mpolyd_t Ad, Bd, Gd, Abar, Bbar;
    nmod_mpolyd_ctx_t dctx;
    slong nvars = ctx->minfo->nvars;

    success = 1;

    if (nmod_mpoly_is_zero(A, ctx)) {
        nmod_mpoly_set(G, B, ctx);
        goto cleanup_stage0;
    }
    if (nmod_mpoly_is_zero(B, ctx)) {
        nmod_mpoly_set(G, A, ctx);
        goto cleanup_stage0;
    }

    nmod_mpolyd_ctx_init(dctx, nvars);
    success = nmod_mpolyd_ctx_set_for_gcd(dctx, A, B, ctx);
    if (!success)
    {
        nmod_mpoly_zero(G, ctx);
        goto cleanup_stage1;
    }

    nmod_mpolyd_init(Ad, nvars);
    nmod_mpolyd_init(Bd, nvars);
    nmod_mpolyd_init(Gd, nvars);
    nmod_mpolyd_init(Abar, nvars);
    nmod_mpolyd_init(Bbar, nvars);

    nmod_mpoly_convert_to_nmod_mpolyd(Ad, dctx, A, ctx);
    nmod_mpoly_convert_to_nmod_mpolyd(Bd, dctx, B, ctx);
    success = nmod_mpolyd_gcd_brown_smprime(Gd, Abar, Bbar, Ad, Bd, ctx->ffinfo);
    if (!success)
    {
        nmod_mpoly_convert_to_nmod_mpolyd(Ad, dctx, A, ctx);
        nmod_mpoly_convert_to_nmod_mpolyd(Bd, dctx, B, ctx);
        success = nmod_mpolyd_gcd_brown_lgprime(Gd, Abar, Bbar, Ad, Bd, ctx->ffinfo);
        if (!success) {
            nmod_mpoly_zero(G, ctx);
        } else {
            nmod_mpoly_convert_from_nmod_mpolyd(G, ctx, Gd, dctx);
        }
    } else
    {
        nmod_mpoly_convert_from_nmod_mpolyd(G, ctx, Gd, dctx);
    }

    nmod_mpolyd_clear(Bbar);
    nmod_mpolyd_clear(Abar);
    nmod_mpolyd_clear(Gd);
    nmod_mpolyd_clear(Bd);
    nmod_mpolyd_clear(Ad);

cleanup_stage1:

    nmod_mpolyd_ctx_clear(dctx);

cleanup_stage0:

    if (!nmod_mpoly_is_zero(G, ctx))
        nmod_mpoly_make_monic(G, G, ctx);

    return success;
}





void nmod_mpoly_to_nmod_poly_keepbits(nmod_poly_t A, slong * Ashift,
               const nmod_mpoly_t B, slong var, const nmod_mpoly_ctx_t ctx);

void nmod_mpoly_from_nmod_poly_keepbits(nmod_mpoly_t A, const nmod_poly_t B,
                           slong Bshift, slong var, mp_bitcnt_t bits, const nmod_mpoly_ctx_t ctx);


int nmod_mpoly_gcd_brownnew(nmod_mpoly_t G, const nmod_mpoly_t A,
                              const nmod_mpoly_t B, const nmod_mpoly_ctx_t ctx)
{
    int success;
    slong * perm;
    slong i;
    mp_bitcnt_t new_bits;
    nmod_mpoly_ctx_t uctx;
    nmod_mpolyu_t Au, Bu, Gu, Abaru, Bbaru;

    if (nmod_mpoly_is_zero(A, ctx))
    {
        if (nmod_mpoly_is_zero(B, ctx))
        {
            nmod_mpoly_zero(G, ctx);
        }
        else
        {
            nmod_mpoly_make_monic(G, B, ctx);
        }
        return 1;
    }

    if (nmod_mpoly_is_zero(B, ctx))
    {
        nmod_mpoly_make_monic(G, A, ctx);
        return 1;
    }

    if (A->bits > FLINT_BITS || B->bits > FLINT_BITS)
    {
        return 0;
    }

    if (ctx->minfo->nvars == 1)
    {
        slong shiftA, shiftB;
        nmod_poly_t a, b, g;
        nmod_poly_init(a, ctx->ffinfo->mod.n);
        nmod_poly_init(b, ctx->ffinfo->mod.n);
        nmod_poly_init(g, ctx->ffinfo->mod.n);
        nmod_mpoly_to_nmod_poly_keepbits(a, &shiftA, A, 0, ctx);
        nmod_mpoly_to_nmod_poly_keepbits(b, &shiftB, B, 0, ctx);
        nmod_poly_gcd(g, a, b);
        nmod_mpoly_from_nmod_poly_keepbits(G, g, FLINT_MIN(shiftA, shiftB), 0, A->bits, ctx);
        nmod_poly_clear(a);
        nmod_poly_clear(b);
        nmod_poly_clear(g);
        return 1;
    }

    perm = (slong *) flint_malloc(ctx->minfo->nvars*sizeof(slong));
    for (i = 0; i < ctx->minfo->nvars; i++)
    {
        perm[i] = i;
    }

    new_bits = FLINT_MAX(A->bits, B->bits);

    nmod_mpoly_ctx_init(uctx, ctx->minfo->nvars - 1, ORD_LEX, ctx->ffinfo->mod.n);
    nmod_mpolyu_init(Au, new_bits, uctx);
    nmod_mpolyu_init(Bu, new_bits, uctx);
    nmod_mpolyu_init(Gu, new_bits, uctx);
    nmod_mpolyu_init(Abaru, new_bits, uctx);
    nmod_mpolyu_init(Bbaru, new_bits, uctx);

    nmod_mpoly_to_mpolyu_perm(Au, A, perm, uctx, ctx);
    nmod_mpoly_to_mpolyu_perm(Bu, B, perm, uctx, ctx);

    success = nmod_mpolyu_gcd_brown_smprime(Gu, Abaru, Bbaru, Au, Bu, uctx->minfo->nvars - 1, uctx);
    if (!success)
    {
        success = nmod_mpolyu_gcd_brown_lgprime(Gu, Abaru, Bbaru, Au, Bu, uctx->minfo->nvars - 1, uctx);
    }
    if (success)
    {
        nmod_mpoly_from_mpolyu_perm(G, Gu, 1, perm, uctx, ctx);
        nmod_mpoly_make_monic(G, G, ctx);
        success = 1;
    }

    nmod_mpolyu_clear(Au, uctx);
    nmod_mpolyu_clear(Bu, uctx);
    nmod_mpolyu_clear(Gu, uctx);
    nmod_mpolyu_clear(Abaru, uctx);
    nmod_mpolyu_clear(Bbaru, uctx);
    nmod_mpoly_ctx_clear(uctx);

    flint_free(perm);

    return success;
}

