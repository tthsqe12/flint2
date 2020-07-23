/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "nmod_mpoly_factor.h"


int nmod_mpoly_factor_fix_units(
    nmod_mpoly_factor_t A,
    const nmod_mpoly_ctx_t ctx)
{
    int changed = 0;
    slong i;
    mp_limb_t c;

    for (i = 0; i < A->num; i++)
    {
        nmod_mpoly_struct * Ai = A->poly + i;
        if (Ai->length > 0 && Ai->coeffs[0] != 1)
        {
            changed = 1;
            c = Ai->coeffs[0];
            nmod_mpoly_scalar_mul_ui(Ai, Ai, nmod_inv(c, ctx->ffinfo->mod), ctx);
            c = nmod_pow_fmpz(c, A->exp + i, ctx->ffinfo->mod);
            A->constant = nmod_mul(A->constant, c, ctx->ffinfo->mod);
        }
    }

    return changed;
}

