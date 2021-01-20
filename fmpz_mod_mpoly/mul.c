/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mod_mpoly.h"
#include "long_extras.h"


static int _try_dense(
    slong * Bdegs,
    slong * Cdegs,
    slong Blen,
    slong Clen,
    slong nvars)
{
    const int max_bit_size = FLINT_MIN(FLINT_BITS/3 + 16, FLINT_BITS - 4);
    slong i, t, product_count, dense_size;

    FLINT_ASSERT(Blen > 0);
    FLINT_ASSERT(Clen > 0);

    dense_size = 1;
    for (i = 0; i < nvars; i++)
    {
        if (z_add_checked(&t, Bdegs[i], Cdegs[i]) ||
            z_add_checked(&t, t, 1) ||
            z_mul_checked(&dense_size, dense_size, t))
        {
            return 0;
        }
    }

    if (dense_size >= WORD(1) << max_bit_size)
        return 0;

    if (z_mul_checked(&product_count, Blen, Clen))
        return 1;

    return dense_size < product_count/32;
}

void fmpz_mod_mpoly_mul(
    fmpz_mod_mpoly_t A,
    const fmpz_mod_mpoly_t B,
    const fmpz_mod_mpoly_t C,
    const fmpz_mod_mpoly_ctx_t ctx)
{
    slong i;
    slong nvars = ctx->minfo->nvars;
    slong * Bdegs, * Cdegs;
    fmpz * maxBfields, * maxCfields;
    slong min_length, max_length;
    TMP_INIT;

    if (B->length < 1 || C->length < 1)
    {
        fmpz_mod_mpoly_zero(A, ctx);
        return;
    }

    TMP_START;

    maxBfields = TMP_ARRAY_ALLOC(2*ctx->minfo->nfields, fmpz);
    maxCfields = maxBfields + ctx->minfo->nfields;
    for (i = 0; i < 2*ctx->minfo->nfields; i++)
        fmpz_init(maxBfields + i);

    mpoly_max_fields_fmpz(maxBfields, B->exps, B->length, B->bits, ctx->minfo);
    mpoly_max_fields_fmpz(maxCfields, C->exps, C->length, C->bits, ctx->minfo);

    min_length = FLINT_MIN(B->length, C->length);
    max_length = FLINT_MAX(B->length, C->length);

    if (min_length < 20 || max_length < 50 ||
        B->bits > FLINT_BITS || C->bits > FLINT_BITS)
    {
        goto do_heap;
    }

    /*
        The multiplication is not trivial and each packed field fits
        into one word. In particular, the degrees must fit an slong.
    */
    Bdegs = TMP_ARRAY_ALLOC(2*nvars, slong);
    Cdegs = Bdegs + nvars;
    mpoly_get_monomial_ui_unpacked_ffmpz((ulong *)Bdegs, maxBfields, ctx->minfo);
    mpoly_get_monomial_ui_unpacked_ffmpz((ulong *)Cdegs, maxCfields, ctx->minfo);

    if (_try_dense(Bdegs, Cdegs, B->length, C->length, nvars))
    {
        if (_fmpz_mod_mpoly_mul_dense_maxfields(A, B, maxBfields, C, maxCfields, ctx))
            goto cleanup;
    }

do_heap:

    _fmpz_mod_mpoly_mul_johnson_maxfields(A, B, maxBfields, C, maxCfields, ctx);

cleanup:

    for (i = 0; i < 2*ctx->minfo->nfields; i++)
        fmpz_clear(maxBfields + i);

    TMP_END;
}

