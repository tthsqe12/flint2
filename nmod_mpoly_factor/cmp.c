/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "nmod_mpoly_factor.h"


int nmod_mpoly_factor_cmp(
    const nmod_mpoly_factor_t A,
    const nmod_mpoly_factor_t B,
    const nmod_mpoly_ctx_t ctx)
{
    int cmp;
    slong i;

    if (A->constant != B->constant)
        return A->constant > B->constant ? 1 : -1;

    if (A->num != B->num)
        return A->num > B->num ? 1 : -1;

    for (i = 0; i < A->num; i++)
    {
        cmp = fmpz_cmp(A->exp + i, B->exp + i);
        if (cmp != 0)
            return cmp;

        cmp = nmod_mpoly_cmp(A->poly + i, B->poly + i, ctx);
        if (cmp != 0)
            return cmp;
    }

    return 0;
}