/*
    Copyright (C) 2019 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fq_nmod_mpoly_factor.h"


int fq_nmod_mpoly_univar_content_mpoly(
    fq_nmod_mpoly_t g,
    const fq_nmod_mpoly_univar_t A,
    const fq_nmod_mpoly_ctx_t ctx)
{
    slong i;

    fq_nmod_mpoly_zero(g, ctx);
    for (i = 0; i < A->length; i++)
    {
		if (!fq_nmod_mpoly_gcd(g, g, A->coeffs + i, ctx))
			return 0;
    }

    return 1;
}

void fq_nmod_mpoly_univar_divexact_mpoly(
    fq_nmod_mpoly_univar_t A,
    const fq_nmod_mpoly_t b,
    const fq_nmod_mpoly_ctx_t ctx)
{
    int success;
    slong i;

    for (i = 0; i < A->length; i++)
    {
        success = fq_nmod_mpoly_divides(A->coeffs + i, A->coeffs + i, b, ctx);
        FLINT_ASSERT(success);
    }
}
