/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fmpq_mpoly_factor.h"


void fmpq_mpoly_factor_init(fmpq_mpoly_factor_t f, const fmpq_mpoly_ctx_t ctx)
{
    fmpq_init(f->constant);
	fmpq_one(f->constant);
    f->poly  = NULL;
    f->exp   = NULL;
    f->num   = 0;
    f->alloc = 0;
}


void fmpq_mpoly_factor_init2(fmpq_mpoly_factor_t f, slong alloc, const fmpq_mpoly_ctx_t ctx)
{
    fmpq_init(f->constant);
	fmpq_one(f->constant);

    if (alloc > 0)
    {
        slong i;

        f->exp = (fmpz *) flint_calloc(alloc, sizeof(fmpz));
        f->poly = (fmpq_mpoly_struct *) flint_malloc(alloc *
                                                    sizeof(fmpq_mpoly_struct));
        for (i = 0; i < alloc; i++)
            fmpq_mpoly_init(f->poly + i, ctx);

        f->alloc = alloc;
    }
    else
    {
        f->exp = NULL;
        f->poly = NULL;
        f->alloc = 0;
    }

    f->num = 0;
}

