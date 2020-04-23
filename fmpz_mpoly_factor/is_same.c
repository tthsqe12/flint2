/*
    Copyright (C) 2019 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fmpz_mpoly_factor.h"


int fmpz_mpoly_factor_is_same(
    const fmpz_mpoly_factor_t A,
    const fmpz_mpoly_factor_t B,
    const fmpz_mpoly_ctx_t ctx)
{
    slong i;

    if (!fmpz_equal(A->content, B->content))
        return 0;

    if (A->length != B->length)
        return 0;

    for (i = 0; i < A->length; i++)
    {
        if (!fmpz_equal(A->exp + i, B->exp + i))
            return 0;

        if (!fmpz_mpoly_equal(A->poly + i, B->poly + i, ctx))
            return 0;
    }

    return 1;
}
