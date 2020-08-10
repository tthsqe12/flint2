/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fq_zech_mpoly.h"

int fq_zech_mpoly_is_fq_zech(const fq_zech_mpoly_t A,
                                                 const fq_zech_mpoly_ctx_t ctx)
{
    slong N;

    if (A->length > 1)
        return 0;

    if (A->length < 1)
        return 1;

    N = mpoly_words_per_exp(A->bits, ctx->minfo);
    return mpoly_monomial_is_zero(A->exps + N*0, N);
}
