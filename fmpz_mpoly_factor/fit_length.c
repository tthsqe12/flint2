/*
    Copyright (C) 2019 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fmpz_mpoly_factor.h"

void fmpz_mpoly_factor_fit_length(fmpz_mpoly_factor_t fac, slong len, const fmpz_mpoly_ctx_t ctx)
{
    if (len > fac->alloc)
    {
        /* At least double number of allocated coeffs */
        if (len < 2 * fac->alloc)
            len = 2 * fac->alloc;
        fmpz_mpoly_factor_realloc(fac, len, ctx);
    }
}

