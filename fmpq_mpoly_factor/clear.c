/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fmpq_mpoly_factor.h"

void fmpq_mpoly_factor_clear(fmpq_mpoly_factor_t f, const fmpq_mpoly_ctx_t ctx)
{
    if (f->alloc > 0)
    {
        slong i;

        for (i = 0; i < f->alloc; i++)
        {
            fmpq_mpoly_clear(f->poly + i, ctx);
			fmpz_clear(f->exp + i);
        }

        flint_free(f->poly);
        flint_free(f->exp);
    }

    fmpq_clear(f->constant);
}

