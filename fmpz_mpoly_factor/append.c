/*
    Copyright (C) 2019 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fmpz_mpoly_factor.h"

void fmpz_mpoly_factor_append(fmpz_mpoly_factor_t fac, const fmpz_mpoly_t p, slong exp, const fmpz_mpoly_ctx_t ctx)
{
    slong i = fac->length;

    fmpz_mpoly_factor_fit_length(fac, i + 1, ctx);

    fmpz_mpoly_set(fac->poly + i, p, ctx);
    fac->exp[i] = exp;
    fac->length = i + 1;
}
