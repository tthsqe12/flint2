/*
    Copyright (C) 2019 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fmpz_mpoly_factor.h"


void fmpz_mpoly_factor_gen(
    fmpz_mpoly_factor_t A,
    slong v,
    const fmpz_mpoly_ctx_t ctx)
{
    fmpz_one(A->content);
    fmpz_mpoly_factor_fit_length(A, 1, ctx);
    fmpz_mpoly_gen(A->poly + 0, v, ctx);
    fmpz_one(A->exp + 0);
    A->length = 1;
}

