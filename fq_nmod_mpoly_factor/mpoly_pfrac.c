/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fq_nmod_mpoly_factor.h"


static void _fq_nmod_to_poly(
    fq_nmod_poly_t A,
    const fq_nmod_mpoly_t B,
    const fq_nmod_mpoly_ctx_t ctx)
{
    slong i;
    ulong * shift, * stride;

    if (B->length == 0)
    {
        fq_nmod_poly_zero(A, ctx->fqctx);
        return;
    }

    shift = (ulong *) flint_malloc(ctx->minfo->nvars*sizeof(ulong));
    stride = (ulong *) flint_malloc(ctx->minfo->nvars*sizeof(ulong));
    for (i = 0; i < ctx->minfo->nvars; i++)
    {
        shift[i] = 0;
        stride[i] = 1;
    }

    _fq_nmod_mpoly_to_fq_nmod_poly_deflate(A, B, 0, shift, stride, ctx);

    flint_free(shift);
    flint_free(stride);
}

static void _fq_nmod_from_poly(
    fq_nmod_mpoly_t A,
    flint_bitcnt_t Abits,
    const fq_nmod_poly_t B,
    const fq_nmod_mpoly_ctx_t ctx)
{
    slong i;
    ulong * shift, * stride;

    if (B->length == 0)
    {
        fq_nmod_mpoly_zero(A, ctx);
        return;
    }

    shift = (ulong *) flint_malloc(ctx->minfo->nvars*sizeof(ulong));
    stride = (ulong *) flint_malloc(ctx->minfo->nvars*sizeof(ulong));
    for (i = 0; i < ctx->minfo->nvars; i++)
    {
        shift[i] = 0;
        stride[i] = 1;
    }

    _fq_nmod_mpoly_from_fq_nmod_poly_inflate(A, Abits, B, 0, shift, stride, ctx);

    flint_free(shift);
    flint_free(stride);
}


void fq_nmod_mpoly_pfrac_init(
    fq_nmod_mpoly_pfrac_t I,
    slong l, slong r,
    const fq_nmod_mpoly_struct * betas,
    const fq_nmod_struct * alpha,
    const fq_nmod_mpoly_ctx_t ctx)
{
    slong i, j, k;
    fq_nmod_poly_t p;
    fq_nmod_poly_t G, S, pq;

    I->l = l;
    I->r = r;

    FLINT_ASSERT(l > 0);

    I->dbetas = (fq_nmod_poly_struct *) flint_malloc(
                                                l*sizeof(fq_nmod_poly_struct));

    I->inv_prod_dbetas = (fq_nmod_poly_struct *) flint_malloc(
                                               l* sizeof(fq_nmod_poly_struct));

    I->prod_mbetas = (fq_nmod_mpoly_struct *) flint_malloc(
                                       (r + 1)*l*sizeof(fq_nmod_mpoly_struct));

    I->mbetas = (fq_nmod_mpoly_struct *) flint_malloc(
                                       (r + 1)*l*sizeof(fq_nmod_mpoly_struct));

    I->deltas = (fq_nmod_mpoly_struct *) flint_malloc(
                                       (r + 1)*l*sizeof(fq_nmod_mpoly_struct));

    fq_nmod_poly_init(p, ctx->fqctx);
    fq_nmod_poly_init(G, ctx->fqctx);
    fq_nmod_poly_init(S, ctx->fqctx);
    fq_nmod_poly_init(pq, ctx->fqctx);

    for (i = r; i >= 0; i--)
        for (j = 0; j < l; j++)
            fq_nmod_mpoly_init(I->deltas + i*l + j, ctx);

    /* set betas */
    i = r;
    for (j = 0; j < l; j++)
    {
        fq_nmod_mpoly_init(I->mbetas + i*l + j, ctx);
        fq_nmod_mpoly_set(I->mbetas + i*l + j, betas + j, ctx);
    }
    for (i--; i >= 0; i--)
    {
        for (j = 0; j < l; j++)
        {
            fq_nmod_mpoly_init(I->mbetas + i*l + j, ctx);
            fq_nmod_mpoly_evaluate_one_fq_nmod(I->mbetas + i*l + j, I->mbetas + (i + 1)*l + j, i + 1, alpha + i, ctx);
        }
    }
    for (j = 0; j < l; j++)
    {
        _fq_nmod_to_poly(p, I->mbetas + 0*l + j, ctx);
        fq_nmod_poly_init(I->dbetas + j, ctx->fqctx);
        fq_nmod_poly_set(I->dbetas + j, p, ctx->fqctx);
    }

    /* set product of betas */
    for (i = r; i >= 0; i--)
    {
        for (j = 0; j < l; j++)
        {
            fq_nmod_mpoly_init(I->prod_mbetas + i*l + j, ctx);
            fq_nmod_mpoly_one(I->prod_mbetas + i*l + j, ctx);
            for (k = 0; k < l; k++)
            {
                if (k == j)
                    continue;
                fq_nmod_mpoly_mul(I->prod_mbetas + i*l + j, I->prod_mbetas + i*l + j, I->mbetas + i*l + k, ctx);
            }
        }        
    }
    for (j = 0; j < l; j++)
    {
        fq_nmod_poly_one(pq, ctx->fqctx);
        for (k = 0; k < l; k++)
        {
            if (k == j)
                continue;
            fq_nmod_poly_mul(pq, pq, I->dbetas + k, ctx->fqctx);
        }
        fq_nmod_poly_init(I->inv_prod_dbetas + j, ctx->fqctx);
        fq_nmod_poly_xgcd(G, S, I->inv_prod_dbetas + j, I->dbetas + j, pq, ctx->fqctx);
        FLINT_ASSERT(fq_nmod_poly_is_one(G, ctx->fqctx));
    }

    fq_nmod_poly_clear(p, ctx->fqctx);
    fq_nmod_poly_clear(G, ctx->fqctx);
    fq_nmod_poly_clear(S, ctx->fqctx);
    fq_nmod_poly_clear(pq, ctx->fqctx);
}


void fq_nmod_mpoly_pfrac_clear(fq_nmod_mpoly_pfrac_t I, const fq_nmod_mpoly_ctx_t ctx)
{
    slong i, j;

    for (j = 0; j < I->l; j++)
    {
        fq_nmod_poly_clear(I->inv_prod_dbetas + j, ctx->fqctx);
        fq_nmod_poly_clear(I->dbetas + j, ctx->fqctx);
        for (i = 0; i <= I->r; i++)
        {
            fq_nmod_mpoly_clear(I->prod_mbetas + i*I->l + j, ctx);
            fq_nmod_mpoly_clear(I->mbetas + i*I->l + j, ctx);
            fq_nmod_mpoly_clear(I->deltas + i*I->l + j, ctx);
        }
    }

    flint_free(I->inv_prod_dbetas);
    flint_free(I->dbetas);
    flint_free(I->prod_mbetas);
    flint_free(I->mbetas);
    flint_free(I->deltas);
}


int fq_nmod_mpoly_pfrac(
    flint_bitcnt_t bits,
    slong r,
    slong num,
    const fq_nmod_mpoly_t t,
    const fq_nmod_struct * alpha,
    const slong * deg,
    const fq_nmod_mpoly_pfrac_t I,
    const fq_nmod_mpoly_ctx_t ctx)
{
    slong i, j;
    int success;
    fq_nmod_mpoly_struct * deltas = I->deltas + r*I->l;
    fq_nmod_mpoly_struct * newdeltas = I->deltas + (r - 1)*I->l;

    if (r < 1)
    {
        fq_nmod_poly_t dtq, S, R;

        fq_nmod_poly_init(dtq, ctx->fqctx);
        fq_nmod_poly_init(S, ctx->fqctx);
        fq_nmod_poly_init(R, ctx->fqctx);

        _fq_nmod_to_poly(dtq, t, ctx);

        success = 1;
        for (i = 0; i < num; i++)
        {
            fq_nmod_poly_mul(S, dtq, I->inv_prod_dbetas + i, ctx->fqctx);
            fq_nmod_poly_rem(R, S, I->dbetas + i, ctx->fqctx);
            _fq_nmod_from_poly(deltas + i, bits, R, ctx);
        }

        fq_nmod_poly_clear(dtq, ctx->fqctx);
        fq_nmod_poly_clear(S, ctx->fqctx);
        fq_nmod_poly_clear(R, ctx->fqctx);
    }
    else
    {
        fq_nmod_mpoly_t e, tt, pow, g, q, newt;

FLINT_ASSERT(num == I->l);

        fq_nmod_mpoly_init(e, ctx);
        fq_nmod_mpoly_init(tt, ctx);
        fq_nmod_mpoly_init(pow, ctx);
        fq_nmod_mpoly_init(g, ctx);
        fq_nmod_mpoly_init(q, ctx);
        fq_nmod_mpoly_init(newt, ctx);
        for (i = 0; i < num; i++)
            fq_nmod_mpoly_zero(deltas + i, ctx);

        fq_nmod_mpoly_set(e, t, ctx);

        fq_nmod_mpoly_one(pow, ctx);
        fq_nmod_mpoly_gen(g, r, ctx);
        fq_nmod_mpoly_sub_fq_nmod(g, g, alpha + r - 1, ctx);
        for (j = 0; j <= deg[r]; j++)
        {
            success = fq_nmod_mpoly_divides(q, e, pow, ctx);
            FLINT_ASSERT(success);
            fq_nmod_mpoly_evaluate_one_fq_nmod(newt, q, r, alpha + r - 1, ctx);
            success = fq_nmod_mpoly_pfrac(bits, r - 1, num, newt, alpha, deg, I, ctx);
            if (!success)
                goto cleanup;

            fq_nmod_mpoly_set(e, t, ctx);
            for (i = 0; i < num; i++)
            {
                fq_nmod_mpoly_mul(tt, newdeltas + i, pow, ctx);
                fq_nmod_mpoly_add(deltas + i, deltas + i, tt, ctx);
                fq_nmod_mpoly_mul(tt, deltas + i, I->prod_mbetas + r*I->l + i, ctx);
                fq_nmod_mpoly_sub(e, e, tt, ctx);
            }

            fq_nmod_mpoly_mul(pow, pow, g, ctx);
        }

        success = fq_nmod_mpoly_is_zero(e, ctx);

cleanup:

        fq_nmod_mpoly_clear(e, ctx);
        fq_nmod_mpoly_clear(tt, ctx);
        fq_nmod_mpoly_clear(pow, ctx);
        fq_nmod_mpoly_clear(g, ctx);
        fq_nmod_mpoly_clear(q, ctx);
        fq_nmod_mpoly_clear(newt, ctx);
    }

    return success;
}

