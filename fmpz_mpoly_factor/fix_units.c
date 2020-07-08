/*
    Copyright (C) 2019 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fmpz_mpoly_factor.h"


int fmpz_mpoly_factor_fix_units(
    fmpz_mpoly_factor_t A,
    const fmpz_mpoly_ctx_t ctx)
{
    int changed = 0;
    slong i;

    for (i = 0; i < A->num; i++)
    {
        fmpz_mpoly_struct * Ai = A->poly + i;
        if (Ai->length > 0 && fmpz_sgn(Ai->coeffs + 0) < 0)
        {
            changed = 1;
            fmpz_mpoly_neg(Ai, Ai, ctx);
            if (fmpz_is_odd(A->exp + i))
                fmpz_neg(A->constant, A->constant);
        }
    }

    return changed;
}

