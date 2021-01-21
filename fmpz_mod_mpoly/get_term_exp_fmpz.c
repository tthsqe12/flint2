/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mod_mpoly.h"

void fmpz_mod_mpoly_get_term_exp_fmpz(fmpz ** exp, const fmpz_mod_mpoly_t A, 
                                           slong i, const fmpz_mod_mpoly_ctx_t ctx)
{
    slong N;

    if (i >= (ulong) A->length)
        flint_throw(FLINT_ERROR,
                       "fmpz_mod_mpoly_get_term_exp_fmpz: index out of range");
 
    N = mpoly_words_per_exp(A->bits, ctx->minfo);
    mpoly_get_monomial_pfmpz(exp, A->exps + N*i, A->bits, ctx->minfo);
}