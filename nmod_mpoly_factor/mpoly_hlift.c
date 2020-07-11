/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "nmod_mpoly_factor.h"


int _hlift_quintic(
    slong m,
    nmod_mpoly_struct * f,
    slong r,
    const mp_limb_t * alpha,
    const nmod_mpoly_t A,
    const slong * degs,
    const nmod_mpoly_ctx_t ctx)
{
    int success;
    slong i, j;
    nmod_mpoly_t e, t, pow, xalpha, q;
    nmod_mpoly_struct * betas, * deltas;
    nmod_mpoly_pfrac_t I;
    flint_bitcnt_t bits = A->bits;

    FLINT_ASSERT(r > 1);

    nmod_mpoly_init(e, ctx);
    nmod_mpoly_init(t, ctx);
    nmod_mpoly_init(pow, ctx);
    nmod_mpoly_init(xalpha, ctx);
    nmod_mpoly_init(q, ctx);

    betas  = (nmod_mpoly_struct * ) flint_malloc(r*sizeof(nmod_mpoly_struct));
    for (i = 0; i < r; i++)
    {
        nmod_mpoly_init(betas + i, ctx);
        nmod_mpoly_repack_bits_inplace(f + i, bits, ctx);
        nmod_mpoly_evaluate_one_ui(betas + i, f + i, m, alpha[m - 1], ctx);
    }

    nmod_mpoly_mul(t, f + 0, f + 1, ctx);
    for (i = 2; i < r; i++)
        nmod_mpoly_mul(t, t, f + i, ctx);
    nmod_mpoly_sub(e, A, t, ctx);

    nmod_mpoly_one(pow, ctx);
    nmod_mpoly_repack_bits_inplace(pow, bits, ctx);

    nmod_mpoly_gen(xalpha, m, ctx);
    nmod_mpoly_sub_ui(xalpha, xalpha, alpha[m - 1], ctx);
    nmod_mpoly_repack_bits_inplace(xalpha, bits, ctx);

    nmod_mpoly_pfrac_init(I, r, m - 1, betas, alpha, ctx);
    deltas = I->deltas + (m - 1)*I->l;

    for (j = 1; j <= degs[m]; j++)
    {
        if (nmod_mpoly_is_zero(e, ctx))
        {
            success = 1;
            goto cleanup;
        }

        nmod_mpoly_mul(pow, pow, xalpha, ctx);
        success = nmod_mpoly_divides(q, e, pow, ctx);
        FLINT_ASSERT(success);
        nmod_mpoly_evaluate_one_ui(t, q, m, alpha[m - 1], ctx);

        success = nmod_mpoly_pfrac(A->bits, m - 1, r, t, alpha, degs, I, ctx);
        if (success < 1)
        {
            success = 0;
            goto cleanup;
        }

        for (i = 0; i < r; i++)
        {
            nmod_mpoly_mul(t, deltas + i, pow, ctx);
            nmod_mpoly_add(f + i, f + i, t, ctx);
        }

        nmod_mpoly_mul(t, f + 0, f + 1, ctx);
        for (i = 2; i < r; i++)
            nmod_mpoly_mul(t, t, f + i, ctx);
        nmod_mpoly_sub(e, A, t, ctx);
    }

    success = nmod_mpoly_is_zero(e, ctx);

cleanup:

    nmod_mpoly_pfrac_clear(I, ctx);

    nmod_mpoly_clear(e, ctx);
    nmod_mpoly_clear(t, ctx);
    nmod_mpoly_clear(pow, ctx);
    nmod_mpoly_clear(xalpha, ctx);
    nmod_mpoly_clear(q, ctx);

    for (i = 0; i < r; i++)
        nmod_mpoly_clear(betas + i, ctx);

    flint_free(betas);

    return success;
}


/* should have A = prod_i f[i] mod (gen(m) - alpha[m-1]) */
int nmod_mpoly_hlift(
    slong m,
    nmod_mpoly_struct * f,  /* length r */
    slong r,
    const mp_limb_t * alpha,
    const nmod_mpoly_t A,
    const slong * degs,
    const nmod_mpoly_ctx_t ctx)
{
    FLINT_ASSERT(r >= 2);
    FLINT_ASSERT(A->bits <= FLINT_BITS);
    FLINT_ASSERT(ctx->minfo->ord == ORD_LEX);

    return _hlift_quintic(m, f, r, alpha, A, degs, ctx);
}
