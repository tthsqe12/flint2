/*
    Copyright (C) 2019 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fq_zech_mpoly.h"


int fq_zech_mpoly_gcd(fq_zech_mpoly_t G, const fq_zech_mpoly_t A,
                       const fq_zech_mpoly_t B, const fq_zech_mpoly_ctx_t ctx)
{
    if (fq_zech_mpoly_is_zero(A, ctx))
    {
        if (fq_zech_mpoly_is_zero(B, ctx))
            fq_zech_mpoly_zero(G, ctx);
        else
            fq_zech_mpoly_make_monic(G, B, ctx);
        return 1;
    }

    if (fq_zech_mpoly_is_zero(B, ctx))
    {
        fq_zech_mpoly_make_monic(G, A, ctx);
        return 1;
    }

    flint_printf("fq_zech_mpoly_gcd not implemented\n");
    flint_abort();

    return 0;
}
