/*
    Copyright (C) 2019 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fmpz_mpoly_factor.h"


int fmpz_mpoly_factor_mul(
    fmpz_mpoly_factor_t A,
    const fmpz_mpoly_factor_t B,
    const fmpz_mpoly_factor_t C,
    const fmpz_mpoly_ctx_t ctx)
{
    slong i, Ti;
    fmpz_mpoly_factor_t T;

    fmpz_mpoly_factor_init(T, ctx);

    fmpz_mul(T->constant, B->constant, C->constant);

    FLINT_ASSERT(B->num >= 0);
    FLINT_ASSERT(C->num >= 0);
    fmpz_mpoly_factor_fit_length(T, B->num + C->num, ctx);

    Ti = 0;

    for (i = 0; i < B->num; i++)
    {
        fmpz_mpoly_set(T->poly + Ti, B->poly + i, ctx);
        fmpz_set(T->exp + Ti, B->exp + i);
        Ti++;
    }

    for (i = 0; i < C->num; i++)
    {
        fmpz_mpoly_set(T->poly + Ti, C->poly + i, ctx);
        fmpz_set(T->exp + Ti, C->exp + i);
        Ti++;
    }

    T->num = Ti;

    fmpz_mpoly_factor_swap(A, T, ctx);
    fmpz_mpoly_factor_clear(T, ctx);

    return 1;
}

