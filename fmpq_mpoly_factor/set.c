/*
    Copyright (C) 2019 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fmpq_mpoly_factor.h"

void fmpq_mpoly_factor_set(fmpq_mpoly_factor_t A, const fmpq_mpoly_factor_t B,
                                                    const fmpq_mpoly_ctx_t ctx)
{
    slong i;

    if (A == B)
        return;

    fmpq_mpoly_factor_fit_length(A, B->num, ctx);
    fmpq_set(A->constant, B->constant);
    for (i = 0; i < B->num; i++)
    {
        fmpq_mpoly_set(A->poly + i, B->poly + i, ctx);
        fmpz_set(A->exp + i, B->exp + i);
    }
    A->num = B->num;
}
