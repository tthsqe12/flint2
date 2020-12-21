/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mod_mpoly.h"


int fmpz_mod_mpoly_divides(
    fmpz_mod_mpoly_t Q,
    const fmpz_mod_mpoly_t A,
    const fmpz_mod_mpoly_t B,
    const fmpz_mod_mpoly_ctx_t ctx)
{
    return fmpz_mod_mpoly_divides_monagan_pearce(Q, A, B, ctx);
}
