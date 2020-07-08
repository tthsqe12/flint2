/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fmpz_mpoly_factor.h"
#include "fmpq_mpoly_factor.h"

int fmpq_mpoly_factor_squarefree(fmpq_mpoly_factor_t f, const fmpq_mpoly_t A,
                                                    const fmpq_mpoly_ctx_t ctx)
{
	slong i;
    int success;
    fmpz_mpoly_factor_t zf;

	fmpz_mpoly_factor_init(zf, ctx->zctx);
	success = fmpz_mpoly_factor_squarefree(zf, A->zpoly, ctx->zctx);

	fmpq_mpoly_factor_fit_length(f, zf->num, ctx);
	fmpq_mul_fmpz(f->constant, A->content, zf->constant);
	for (i = 0; i < zf->num; i++)
	{
		fmpz_swap(f->exp + i, zf->exp + i);
		fmpq_one(f->poly[i].content);
		fmpz_mpoly_swap(f->poly[i].zpoly, zf->poly + i, ctx->zctx);
        fmpq_mpoly_reduce(f->poly + i, ctx); /* just in case */
	}
	f->num = zf->num;

	fmpz_mpoly_factor_clear(zf, ctx->zctx);

	return success;
}

