/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fq_nmod_mpoly_factor.h"


static int _hlift_quintic(
    slong m,
    fq_nmod_mpoly_struct * f,
    slong r,
    const fq_nmod_struct * alpha,
    const fq_nmod_mpoly_t A,
    const slong * degs,
    const fq_nmod_mpoly_ctx_t ctx)
{
    int success;
    slong i, j;
    fq_nmod_mpoly_t e, t, pow, xalpha, q;
    fq_nmod_mpoly_struct * betas, * deltas;
    fq_nmod_mpoly_pfrac_t I;
    flint_bitcnt_t bits = A->bits;

    FLINT_ASSERT(r > 1);

    fq_nmod_mpoly_init(e, ctx);
    fq_nmod_mpoly_init(t, ctx);
    fq_nmod_mpoly_init(pow, ctx);
    fq_nmod_mpoly_init(xalpha, ctx);
    fq_nmod_mpoly_init(q, ctx);

    betas  = (fq_nmod_mpoly_struct *) flint_malloc(r*sizeof(fq_nmod_mpoly_struct));
    for (i = 0; i < r; i++)
    {
        fq_nmod_mpoly_init(betas + i, ctx);
        fq_nmod_mpoly_evaluate_one_fq_nmod(betas + i, f + i, m, alpha + m - 1, ctx);
    }

    fq_nmod_mpoly_mul(t, f + 0, f + 1, ctx);
    for (i = 2; i < r; i++)
        fq_nmod_mpoly_mul(t, t, f + i, ctx);
    fq_nmod_mpoly_sub(e, A, t, ctx);

    fq_nmod_mpoly_one(pow, ctx);
    fq_nmod_mpoly_repack_bits_inplace(pow, bits, ctx);

    fq_nmod_mpoly_gen(xalpha, m, ctx);
    fq_nmod_mpoly_sub_fq_nmod(xalpha, xalpha, alpha + m - 1, ctx);
    fq_nmod_mpoly_repack_bits_inplace(xalpha, bits, ctx);

    fq_nmod_mpoly_pfrac_init(I, r, m - 1, betas, alpha, ctx);
    deltas = I->deltas + (m - 1)*I->l;

    for (j = 1; j <= degs[m]; j++)
    {
        if (fq_nmod_mpoly_is_zero(e, ctx))
        {
            success = 1;
            goto cleanup;
        }

        fq_nmod_mpoly_mul(pow, pow, xalpha, ctx);
        success = fq_nmod_mpoly_divides(q, e, pow, ctx);
        FLINT_ASSERT(success);
        fq_nmod_mpoly_evaluate_one_fq_nmod(t, q, m, alpha + m - 1, ctx);

        success = fq_nmod_mpoly_pfrac(bits, m - 1, r, t, alpha, degs, I, ctx);
        if (!success)
            goto cleanup;

        for (i = 0; i < r; i++)
        {
            fq_nmod_mpoly_mul(t, deltas + i, pow, ctx);
            fq_nmod_mpoly_add(f + i, f + i, t, ctx);
        }

        fq_nmod_mpoly_mul(t, f + 0, f + 1, ctx);
        for (i = 2; i < r; i++)
            fq_nmod_mpoly_mul(t, t, f + i, ctx);
        fq_nmod_mpoly_sub(e, A, t, ctx);
    }

    success = fq_nmod_mpoly_is_zero(e, ctx);

cleanup:

    fq_nmod_mpoly_pfrac_clear(I, ctx);

    fq_nmod_mpoly_clear(e, ctx);
    fq_nmod_mpoly_clear(t, ctx);
    fq_nmod_mpoly_clear(pow, ctx);
    fq_nmod_mpoly_clear(xalpha, ctx);
    fq_nmod_mpoly_clear(q, ctx);

    for (i = 0; i < r; i++)
        fq_nmod_mpoly_clear(betas + i, ctx);

    flint_free(betas);

    return success;
}

/* should have A = prod_i f[i] mod (gen(m) - alpha[m-1]) */
int fq_nmod_mpoly_hlift(
    slong m,
    fq_nmod_mpoly_struct * f, /* length r */
    slong r,
    const fq_nmod_struct * alpha,
    const fq_nmod_mpoly_t A,
    const slong * degs,
    const fq_nmod_mpoly_ctx_t ctx)
{
    FLINT_ASSERT(r >= 2);
    FLINT_ASSERT(A->bits <= FLINT_BITS);
    FLINT_ASSERT(ctx->minfo->ord == ORD_LEX);

    return _hlift_quintic(m, f, r, alpha, A, degs, ctx);
}
