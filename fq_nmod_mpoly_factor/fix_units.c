/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fq_nmod_mpoly_factor.h"


int fq_nmod_mpoly_factor_fix_units(
    fq_nmod_mpoly_factor_t A,
    const fq_nmod_mpoly_ctx_t ctx)
{
    int changed = 0;
    slong i;
    fq_nmod_t c;

    fq_nmod_init(c, ctx->fqctx);

    for (i = 0; i < A->num; i++)
    {
        fq_nmod_mpoly_struct * Ai = A->poly + i;
        if (Ai->length > 0 && !fq_nmod_is_one(Ai->coeffs + 0, ctx->fqctx))
        {
            changed = 1;
            fq_nmod_pow(c, Ai->coeffs + 0, A->exp + i, ctx->fqctx);
            fq_nmod_mul(A->constant, A->constant, c, ctx->fqctx);
            fq_nmod_inv(c, c, ctx->fqctx);
            fq_nmod_mpoly_scalar_mul_fq_nmod(Ai, Ai, c, ctx);
        }
    }

    fq_nmod_clear(c, ctx->fqctx);

    return changed;
}

