/*
    Copyright (C) 2019 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fmpz_mpoly_factor.h"


int fmpz_mpoly_factor_pow_fmpz(
    fmpz_mpoly_factor_t A,
    const fmpz_mpoly_factor_t B,
    const fmpz_t e,
    const fmpz_mpoly_ctx_t ctx)
{
    slong i;

    fmpz_mpoly_factor_set(A, B, ctx);

    if (!fmpz_pow_fmpz(A->content, A->content, e))
        return 0;

    for (i = 0; i < A->length; i++)
        fmpz_mul(A->exp + i, A->exp + i, e);

    return 1;
}

