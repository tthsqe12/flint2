/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fmpz_mpoly_factor.h"


void fmpz_mpoly_factor_scalar_mul_si(
    fmpz_mpoly_factor_t A,
    const fmpz_mpoly_factor_t B,
    slong c,
    const fmpz_mpoly_ctx_t ctx)
{
    fmpz_mpoly_factor_set(A, B, ctx);
    fmpz_mul_si(A->constant, A->constant, c);
}
