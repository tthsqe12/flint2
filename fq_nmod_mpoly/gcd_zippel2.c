/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fq_nmod_mpoly.h"


int fq_nmod_mpoly_gcd_zippel2(
    fq_nmod_mpoly_t G,
    const fq_nmod_mpoly_t A,
    const fq_nmod_mpoly_t B,
    const fq_nmod_mpoly_ctx_t ctx)
{
    int success;
    flint_bitcnt_t wbits = FLINT_MAX(A->bits, B->bits);
    slong i, m = ctx->minfo->nvars;
    fq_nmod_mpoly_ctx_t uctx;
    fq_nmod_mpolyu_t Auu, Buu, Guu, Abaruu, Bbaruu;
    fq_nmod_mpoly_t Ac, Bc, Gc, Gamma;
    slong * perm, * A_degs, * B_degs, * Gamma_degs;
    ulong * shift, * stride;

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

    if (m < 3)
    {
        return fq_nmod_mpoly_gcd_brown(G, A, B, ctx);
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

    fq_nmod_mpoly_ctx_init(uctx, m - 2, ORD_LEX, ctx->fqctx);
    fq_nmod_mpolyu_init(Auu, wbits, uctx);
    fq_nmod_mpolyu_init(Buu, wbits, uctx);
    fq_nmod_mpolyu_init(Guu, wbits, uctx);
    fq_nmod_mpolyu_init(Abaruu, wbits, uctx);
    fq_nmod_mpolyu_init(Bbaruu, wbits, uctx);
    fq_nmod_mpoly_init3(Ac, 0, wbits, uctx);
    fq_nmod_mpoly_init3(Bc, 0, wbits, uctx);
    fq_nmod_mpoly_init3(Gc, 0, wbits, uctx);
    fq_nmod_mpoly_init3(Gamma, 0, wbits, uctx);

    fq_nmod_mpoly_degrees_si(A_degs, A, ctx);
    fq_nmod_mpoly_degrees_si(B_degs, B, ctx);
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

    fq_nmod_mpoly_to_mpolyuu_perm_deflate(Auu, uctx, A, ctx, perm, shift, stride);
    fq_nmod_mpoly_to_mpolyuu_perm_deflate(Buu, uctx, B, ctx, perm, shift, stride);

    success = fq_nmod_mpolyu_content_mpoly(Ac, Auu, uctx) &&
              fq_nmod_mpolyu_content_mpoly(Bc, Buu, uctx) &&
              _fq_nmod_mpoly_gcd(Gc, wbits, Ac, Bc, uctx);
    if (!success)
        goto cleanup;

    fq_nmod_mpolyu_divexact_mpoly_inplace(Auu, Ac, uctx);
    fq_nmod_mpolyu_divexact_mpoly_inplace(Buu, Bc, uctx);

    FLINT_ASSERT(Auu->bits == wbits);
    FLINT_ASSERT(Buu->bits == wbits);
    FLINT_ASSERT(Auu->length > 0);
    FLINT_ASSERT(Buu->length > 0);
    FLINT_ASSERT(Ac->bits == wbits);
    FLINT_ASSERT(Bc->bits == wbits);

    success = _fq_nmod_mpoly_gcd(Gamma, wbits, Auu->coeffs + 0, Buu->coeffs + 0, uctx);
    if (!success)
        goto cleanup;

    fq_nmod_mpoly_degrees_si(Gamma_degs, Gamma, uctx);
    fq_nmod_mpolyu_degrees_si(A_degs, Auu, uctx);
    fq_nmod_mpolyu_degrees_si(B_degs, Buu, uctx);

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

            fq_nmod_mpolyu_repack_bits_inplace(Guu, wbits, uctx);
            fq_nmod_mpolyu_repack_bits_inplace(Abaruu, wbits, uctx);
            fq_nmod_mpolyu_repack_bits_inplace(Bbaruu, wbits, uctx);
            fq_nmod_mpolyu_repack_bits_inplace(Auu, wbits, uctx);
            fq_nmod_mpolyu_repack_bits_inplace(Buu, wbits, uctx);
            fq_nmod_mpoly_repack_bits_inplace(Gamma, wbits, uctx);
            fq_nmod_mpoly_repack_bits_inplace(Gc, wbits, uctx);

            break;
        }
    }

    success = fq_nmod_mpolyuu_gcd_zippel_smprime(Guu, NULL, Abaruu, Bbaruu,
                            Auu, A_degs, Buu, B_degs, Gamma, Gamma_degs, uctx);
    if (!success)
    {
        success = fq_nmod_mpolyuu_gcd_zippel_lgprime(Guu, NULL, Abaruu, Bbaruu,
                            Auu, A_degs, Buu, B_degs, Gamma, Gamma_degs, uctx);
        if (!success)
            goto cleanup;
    }

    fq_nmod_mpolyu_mul_mpoly_inplace(Guu, Gc, uctx);

    fq_nmod_mpoly_from_mpolyuu_perm_inflate(G, wbits, ctx, Guu, uctx,
                                                          perm, shift, stride);
    success = 1;

    fq_nmod_mpoly_make_monic(G, G, ctx);

cleanup:

    fq_nmod_mpolyu_clear(Auu, uctx);
    fq_nmod_mpolyu_clear(Buu, uctx);
    fq_nmod_mpolyu_clear(Guu, uctx);
    fq_nmod_mpolyu_clear(Abaruu, uctx);
    fq_nmod_mpolyu_clear(Bbaruu, uctx);
    fq_nmod_mpoly_clear(Ac, uctx);
    fq_nmod_mpoly_clear(Bc, uctx);
    fq_nmod_mpoly_clear(Gc, uctx);
    fq_nmod_mpoly_clear(Gamma, uctx);
    fq_nmod_mpoly_ctx_clear(uctx);
    flint_free(perm);
    flint_free(shift);
    flint_free(stride);
    flint_free(A_degs);
    flint_free(B_degs);
    flint_free(Gamma_degs);

    return success;
}
