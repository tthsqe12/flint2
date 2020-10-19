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
    nmod_mpoly_ctx_t uctx;
    nmod_mpolyu_t Auu, Buu, Guu, Abaruu, Bbaruu;
    nmod_mpoly_t Ac, Bc, Gc, Gamma;
    slong * perm, * A_degs, * B_degs, * Gamma_degs;
    ulong * shift, * stride;

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

    nmod_mpoly_ctx_init(uctx, m - 2, ORD_LEX, ctx->ffinfo->mod.n);
    nmod_mpolyu_init(Auu, wbits, uctx);
    nmod_mpolyu_init(Buu, wbits, uctx);
    nmod_mpolyu_init(Guu, wbits, uctx);
    nmod_mpolyu_init(Abaruu, wbits, uctx);
    nmod_mpolyu_init(Bbaruu, wbits, uctx);
    nmod_mpoly_init3(Ac, 0, wbits, uctx);
    nmod_mpoly_init3(Bc, 0, wbits, uctx);
    nmod_mpoly_init3(Gc, 0, wbits, uctx);
    nmod_mpoly_init3(Gamma, 0, wbits, uctx);

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

    nmod_mpoly_to_mpolyuu_perm_deflate_threaded_pool(Auu, uctx, A, ctx,
                                           perm, shift, stride, NULL, NULL, 0);
    nmod_mpoly_to_mpolyuu_perm_deflate_threaded_pool(Buu, uctx, B, ctx,
                                           perm, shift, stride, NULL, NULL, 0);

    success = nmod_mpolyu_content_mpoly_threaded_pool(Ac, Auu, uctx, NULL, 0) &&
              nmod_mpolyu_content_mpoly_threaded_pool(Bc, Buu, uctx, NULL, 0) &&
              _nmod_mpoly_gcd_threaded_pool(Gc, wbits, Ac, Bc, uctx, NULL, 0);
    if (!success)
        goto cleanup;

    nmod_mpolyu_divexact_mpoly_inplace(Auu, Ac, uctx);
    nmod_mpolyu_divexact_mpoly_inplace(Buu, Bc, uctx);

    FLINT_ASSERT(Auu->bits == wbits);
    FLINT_ASSERT(Buu->bits == wbits);
    FLINT_ASSERT(Auu->length > 0);
    FLINT_ASSERT(Buu->length > 0);
    FLINT_ASSERT(Ac->bits == wbits);
    FLINT_ASSERT(Bc->bits == wbits);

    success = _nmod_mpoly_gcd_threaded_pool(Gamma, wbits, Auu->coeffs + 0,
                                               Buu->coeffs + 0, uctx, NULL, 0);
    if (!success)
        goto cleanup;

    nmod_mpoly_degrees_si(Gamma_degs, Gamma, uctx);
    nmod_mpolyu_degrees_si(A_degs, Auu, uctx);
    nmod_mpolyu_degrees_si(B_degs, Buu, uctx);

    for (i = 0; i < m - 2; i++)
    {
        ulong this_deg = (ulong) Gamma_degs[i] + FLINT_MAX(A_degs[i], B_degs[i]);
        if ((FLINT_BIT_COUNT(this_deg)) >= wbits)
        {
            wbits = mpoly_fix_bits(wbits + 1, uctx->minfo);
            if (wbits > FLINT_BITS)
            {
                success = 0;
                goto cleanup;
            }

            nmod_mpolyu_repack_bits_inplace(Guu, wbits, uctx);
            nmod_mpolyu_repack_bits_inplace(Abaruu, wbits, uctx);
            nmod_mpolyu_repack_bits_inplace(Bbaruu, wbits, uctx);
            nmod_mpolyu_repack_bits_inplace(Auu, wbits, uctx);
            nmod_mpolyu_repack_bits_inplace(Buu, wbits, uctx);
            nmod_mpoly_repack_bits_inplace(Gamma, wbits, uctx);
            nmod_mpoly_repack_bits_inplace(Gc, wbits, uctx);

            break;
        }
    }

    success = nmod_mpolyuu_gcd_zippel_smprime(Guu, NULL, Abaruu, Bbaruu,
                            Auu, A_degs, Buu, B_degs, Gamma, Gamma_degs, uctx);
    if (!success)
    {
        success = nmod_mpolyuu_gcd_zippel_lgprime(Guu, NULL, Abaruu, Bbaruu,
                            Auu, A_degs, Buu, B_degs, Gamma, Gamma_degs, uctx);
        if (!success)
            goto cleanup;
    }

    nmod_mpolyu_mul_mpoly_inplace(Guu, Gc, uctx);

    nmod_mpoly_from_mpolyuu_perm_inflate(G, wbits, ctx, Guu, uctx,
                                                          perm, shift, stride);
    success = 1;

    nmod_mpoly_make_monic(G, G, ctx);

cleanup:

    nmod_mpolyu_clear(Auu, uctx);
    nmod_mpolyu_clear(Buu, uctx);
    nmod_mpolyu_clear(Guu, uctx);
    nmod_mpolyu_clear(Abaruu, uctx);
    nmod_mpolyu_clear(Bbaruu, uctx);
    nmod_mpoly_clear(Ac, uctx);
    nmod_mpoly_clear(Bc, uctx);
    nmod_mpoly_clear(Gc, uctx);
    nmod_mpoly_clear(Gamma, uctx);
    nmod_mpoly_ctx_clear(uctx);
    flint_free(perm);
    flint_free(shift);
    flint_free(stride);
    flint_free(A_degs);
    flint_free(B_degs);
    flint_free(Gamma_degs);

    return success;
}
