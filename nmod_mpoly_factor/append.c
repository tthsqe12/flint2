/*
    Copyright (C) 2019 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "nmod_mpoly_factor.h"

void nmod_mpoly_factor_append_ui(nmod_mpoly_factor_t fac, const nmod_mpoly_t p, ulong e, const nmod_mpoly_ctx_t ctx)
{
    slong i = fac->length;
    nmod_mpoly_factor_fit_length(fac, i + 1, ctx);
    nmod_mpoly_set(fac->poly + i, p, ctx);
    fmpz_set_ui(fac->exp + i, e);
    fac->length = i + 1;
}

void nmod_mpoly_factor_append_fmpz(nmod_mpoly_factor_t fac, const nmod_mpoly_t p, const fmpz_t e, const nmod_mpoly_ctx_t ctx)
{
    slong i = fac->length;
    nmod_mpoly_factor_fit_length(fac, i + 1, ctx);
    nmod_mpoly_set(fac->poly + i, p, ctx);
    fmpz_set(fac->exp + i, e);
    fac->length = i + 1;
}
