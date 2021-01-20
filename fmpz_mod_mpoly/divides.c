/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mod_mpoly.h"
#include "ulong_extras.h"

static int _fmpz_mod_mpoly_divides_try_dense(
    slong * Adegs,
    slong * Bdegs,
    slong nvars,
    slong Alen,
    slong Blen)
{
    const int max_bit_size = FLINT_MIN(FLINT_BITS/3 + 16, FLINT_BITS - 4);
    slong i;
    ulong total_dense_size;

    FLINT_ASSERT(Alen > 0);
    FLINT_ASSERT(Blen > 0);

    total_dense_size = 1;
    for (i = 0; i < nvars; i++)
    {
        if (n_mul_checked(&total_dense_size, total_dense_size, Adegs[i] + 1))
            return 0;
    }

    return total_dense_size < (UWORD(1) << max_bit_size) &&
           total_dense_size/16 < Alen;
}


int fmpz_mod_mpoly_divides(
    fmpz_mod_mpoly_t Q,
    const fmpz_mod_mpoly_t A,
    const fmpz_mod_mpoly_t B,
    const fmpz_mod_mpoly_ctx_t ctx)
{
    int success;
    slong i, * Adegs, * Bdegs;
    TMP_INIT;

    if (fmpz_mod_mpoly_is_zero(B, ctx))
    {
        if (!fmpz_mod_mpoly_is_zero(A, ctx) ||
            !fmpz_is_one(fmpz_mod_mpoly_ctx_modulus(ctx)))
        {
            flint_throw(FLINT_DIVZERO, "fmpz_mod_mpoly_divides: divide by zero");
        }

        fmpz_mod_mpoly_zero(Q, ctx);
        return 1;
    }

    if (fmpz_mod_mpoly_is_zero(A, ctx))
    {
        fmpz_mod_mpoly_zero(Q, ctx);
        return 1;
    }

    if (A->length < 30)
    {
        return fmpz_mod_mpoly_divides_monagan_pearce(Q, A, B, ctx);
    }

    TMP_START;

    success = -1;

    if (A->bits <= FLINT_BITS &&
        B->bits <= FLINT_BITS &&
        A->length > 50)
    {
        Adegs = (slong *) TMP_ALLOC(ctx->minfo->nvars*sizeof(slong));
        Bdegs = (slong *) TMP_ALLOC(ctx->minfo->nvars*sizeof(slong));

        mpoly_degrees_si(Adegs, A->exps, A->length, A->bits, ctx->minfo);
        mpoly_degrees_si(Bdegs, B->exps, B->length, B->bits, ctx->minfo);

        /* quick degree check */
        for (i = 0; i < ctx->minfo->nvars; i++)
        {
            if (Adegs[i] < Bdegs[i])
            {
                fmpz_mod_mpoly_zero(Q, ctx);
                success = 0;
                goto cleanup;
            }
        }

        if (_fmpz_mod_mpoly_divides_try_dense(Adegs, Bdegs, ctx->minfo->nvars,
                                                         A->length, B->length))
        {
            success = fmpz_mod_mpoly_divides_dense(Q, A, B, ctx);
        }
    }

    if (success < 0)
        success = fmpz_mod_mpoly_divides_monagan_pearce(Q, A, B, ctx);

cleanup:

    TMP_END;
    return success;
}
