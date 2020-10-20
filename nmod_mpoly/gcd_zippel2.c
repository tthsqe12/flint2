/*
    Copyright (C) 2019 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "nmod_mpoly.h"


int nmod_mpoly_gcd_zippel2(
    nmod_mpoly_t G,
    const nmod_mpoly_t A,
    const nmod_mpoly_t B,
    const nmod_mpoly_ctx_t ctx)
{
    int success;
    flint_bitcnt_t wbits = FLINT_MAX(A->bits, B->bits);
    slong i, m = ctx->minfo->nvars;
    nmod_mpoly_ctx_t lctx;
    nmod_mpoly_t Al, Bl, Gl, Abarl, Bbarl;
    nmod_mpoly_t Ac, Bc, Gc, Gamma, Al_lc, Bl_lc;
    slong * perm, * A_degs, * B_degs, * Gamma_degs;
    ulong * shift, * stride;

    if (nmod_mpoly_is_zero(A, ctx))
    {
        if (nmod_mpoly_is_zero(B, ctx))
            nmod_mpoly_zero(G, ctx);
        else
            nmod_mpoly_make_monic(G, B, ctx);
        return 1;
    }

    if (nmod_mpoly_is_zero(B, ctx))
    {
        nmod_mpoly_make_monic(G, A, ctx);
        return 1;
    }

    if (m < 3)
    {
        return nmod_mpoly_gcd_brown(G, A, B, ctx);
    }

    if (wbits > FLINT_BITS)
    {
        return 0;
    }

    perm = FLINT_ARRAY_ALLOC(m, slong);
    shift = FLINT_ARRAY_ALLOC(m, ulong);
    stride = FLINT_ARRAY_ALLOC(m, ulong);
    A_degs = FLINT_ARRAY_ALLOC(m, slong);
    B_degs = FLINT_ARRAY_ALLOC(m, slong);
    Gamma_degs = FLINT_ARRAY_ALLOC(m, slong);

    nmod_mpoly_ctx_init(lctx, m, ORD_LEX, ctx->ffinfo->mod.n);
    nmod_mpoly_init3(Al, 0, wbits, lctx);
    nmod_mpoly_init3(Bl, 0, wbits, lctx);
    nmod_mpoly_init3(Gl, 0, wbits, lctx);
    nmod_mpoly_init3(Abarl, 0, wbits, lctx);
    nmod_mpoly_init3(Bbarl, 0, wbits, lctx);
    nmod_mpoly_init3(Ac, 0, wbits, lctx);
    nmod_mpoly_init3(Bc, 0, wbits, lctx);
    nmod_mpoly_init3(Gc, 0, wbits, lctx);
    nmod_mpoly_init3(Gamma, 0, wbits, lctx);
    nmod_mpoly_init3(Al_lc, 0, wbits, lctx);
    nmod_mpoly_init3(Bl_lc, 0, wbits, lctx);

    nmod_mpoly_degrees_si(A_degs, A, ctx);
    nmod_mpoly_degrees_si(B_degs, B, ctx);
    if (FLINT_BIT_COUNT(A_degs[0]) >= FLINT_BITS/2 ||
        FLINT_BIT_COUNT(A_degs[1]) >= FLINT_BITS/2 ||
        FLINT_BIT_COUNT(B_degs[0]) >= FLINT_BITS/2 ||
        FLINT_BIT_COUNT(B_degs[1]) >= FLINT_BITS/2)
    {
        success = 0;
        goto cleanup;
    }

    for (i = 0; i < m; i++)
    {
        perm[i] = i;
        shift[i] = 0;
        stride[i] = 1;
    }

    nmod_mpoly_to_mpolyl_perm_deflate(Al, lctx, A, ctx, perm, shift, stride);
    nmod_mpoly_to_mpolyl_perm_deflate(Bl, lctx, B, ctx, perm, shift, stride);

    success = nmod_mpolyl_content(Ac, Al, 2, lctx) &&
              nmod_mpolyl_content(Bc, Bl, 2, lctx) &&
              nmod_mpoly_gcd(Gc, Ac, Bc, lctx);
    if (!success)
        goto cleanup;

    nmod_mpoly_divides(Al, Al, Ac, lctx);
    nmod_mpoly_divides(Bl, Bl, Bc, lctx);

    nmod_mpoly_repack_bits_inplace(Al, wbits, lctx);
    nmod_mpoly_repack_bits_inplace(Bl, wbits, lctx);
    nmod_mpolyl_lead_coeff(Al_lc, Al, 2, lctx);
    nmod_mpolyl_lead_coeff(Bl_lc, Bl, 2, lctx);
    success = nmod_mpoly_gcd(Gamma, Al_lc, Bl_lc, lctx);
    if (!success)
        goto cleanup;
    nmod_mpoly_repack_bits_inplace(Gamma, wbits, lctx);

    nmod_mpoly_degrees_si(Gamma_degs, Gamma, lctx);
    nmod_mpoly_degrees_si(A_degs, Al, lctx);
    nmod_mpoly_degrees_si(B_degs, Bl, lctx);

    success = nmod_mpolyl_gcd_zippel_smprime(Gl, NULL, Abarl, Bbarl,
                              Al, A_degs, Bl, B_degs, Gamma, Gamma_degs, lctx);
    if (!success)
    {
        success = nmod_mpolyl_gcd_zippel_lgprime(Gl, NULL, Abarl, Bbarl,
                              Al, A_degs, Bl, B_degs, Gamma, Gamma_degs, lctx);
        if (!success)
            goto cleanup;
    }

    nmod_mpoly_mul(Gl, Gl, Gc, lctx);

    nmod_mpoly_from_mpolyl_perm_inflate(G, wbits, ctx, Gl, lctx, perm, shift, stride);
    success = 1;

    nmod_mpoly_make_monic(G, G, ctx);

cleanup:

    nmod_mpoly_clear(Al, lctx);
    nmod_mpoly_clear(Bl, lctx);
    nmod_mpoly_clear(Gl, lctx);
    nmod_mpoly_clear(Abarl, lctx);
    nmod_mpoly_clear(Bbarl, lctx);
    nmod_mpoly_clear(Ac, lctx);
    nmod_mpoly_clear(Bc, lctx);
    nmod_mpoly_clear(Gc, lctx);
    nmod_mpoly_clear(Gamma, lctx);
    nmod_mpoly_clear(Al_lc, lctx);
    nmod_mpoly_clear(Bl_lc, lctx);
    nmod_mpoly_ctx_clear(lctx);
    flint_free(perm);
    flint_free(shift);
    flint_free(stride);
    flint_free(A_degs);
    flint_free(B_degs);
    flint_free(Gamma_degs);

    return success;
}
