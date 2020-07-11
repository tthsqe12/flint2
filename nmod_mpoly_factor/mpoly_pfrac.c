/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "nmod_mpoly_factor.h"


static void _to_poly(
    nmod_poly_t A,
    const nmod_mpoly_t B,
    const nmod_mpoly_ctx_t ctx)
{
    slong i;
    ulong * shift, * stride;

    if (B->length == 0)
    {
        nmod_poly_zero(A);
        return;
    }

    shift = (ulong *) flint_malloc(ctx->minfo->nvars*sizeof(ulong));
    stride = (ulong *) flint_malloc(ctx->minfo->nvars*sizeof(ulong));
    for (i = 0; i < ctx->minfo->nvars; i++)
    {
        shift[i] = 0;
        stride[i] = 1;
    }

    _nmod_mpoly_to_nmod_poly_deflate(A, B, 0, shift, stride, ctx);

    flint_free(shift);
    flint_free(stride);
}


static void _from_poly(
    nmod_mpoly_t A,
    flint_bitcnt_t Abits,
    const nmod_poly_t B,
    const nmod_mpoly_ctx_t ctx)
{
    slong i;
    ulong * shift, * stride;

    if (B->length == 0)
    {
        nmod_mpoly_zero(A, ctx);
        return;
    }

    shift = (ulong *) flint_malloc(ctx->minfo->nvars*sizeof(ulong));
    stride = (ulong *) flint_malloc(ctx->minfo->nvars*sizeof(ulong));
    for (i = 0; i < ctx->minfo->nvars; i++)
    {
        shift[i] = 0;
        stride[i] = 1;
    }

    _nmod_mpoly_from_nmod_poly_inflate(A, Abits, B, 0, shift, stride, ctx);

    flint_free(shift);
    flint_free(stride);
}



void nmod_mpoly_pfrac_init(
    nmod_mpoly_pfrac_t I,
    slong l, slong r,
    const nmod_mpoly_struct * betas,
    const mp_limb_t * alpha,
    const nmod_mpoly_ctx_t ctx)
{
    slong i, j, k;
    nmod_poly_t p;
    nmod_poly_t G, S, pq;
/*
flint_printf("_mfactor_disolve_init called(l = %wd, r = %wd)\n",l,r);
*/
    I->l = l;
    I->r = r;

    FLINT_ASSERT(l > 0);

    I->inv_prod_dbetas = (nmod_poly_struct *) flint_malloc(l*sizeof(nmod_poly_struct));
    I->dbetas = (nmod_poly_struct *) flint_malloc(l*sizeof(nmod_poly_struct));
    I->prod_mbetas = (nmod_mpoly_struct *) flint_malloc((r + 1)*l*sizeof(nmod_mpoly_struct));
    I->mbetas = (nmod_mpoly_struct *) flint_malloc((r + 1)*l*sizeof(nmod_mpoly_struct));
    I->deltas = (nmod_mpoly_struct *) flint_malloc((r + 1)*l*sizeof(nmod_mpoly_struct));

    nmod_poly_init_mod(p, ctx->ffinfo->mod);
    nmod_poly_init_mod(G, ctx->ffinfo->mod);
    nmod_poly_init_mod(S, ctx->ffinfo->mod);
    nmod_poly_init_mod(pq, ctx->ffinfo->mod);

    /* initialize deltas */
    for (i = r; i >= 0; i--)
        for (j = 0; j < l; j++)
            nmod_mpoly_init(I->deltas + i*l + j, ctx);

    /* set betas */
    i = r;
    for (j = 0; j < l; j++)
    {
        nmod_mpoly_init(I->mbetas + i*l + j, ctx);
        nmod_mpoly_set(I->mbetas + i*l + j, betas + j, ctx);
/*
flint_printf("mbetas[%wd][%wd]: ",i,j); nmod_mpoly_print_pretty(I->mbetas + i*l + j, NULL, ctx); printf("\n");
*/
    }
    for (i--; i >= 0; i--)
    {
        for (j = 0; j < l; j++)
        {
            nmod_mpoly_init(I->mbetas + i*l + j, ctx);
            nmod_mpoly_evaluate_one_ui(I->mbetas + i*l + j, I->mbetas + (i + 1)*l + j, i + 1, alpha[i], ctx);
/*
flint_printf("mbetas[%wd][%wd]: ",i,j); nmod_mpoly_print_pretty(I->mbetas + i*l + j, NULL, ctx); printf("\n");
*/
        }
    }
    for (j = 0; j < l; j++)
    {
        _to_poly(p, I->mbetas + 0*l + j, ctx);
        nmod_poly_init_mod(I->dbetas + j, ctx->ffinfo->mod);
        nmod_poly_set(I->dbetas + j, p);
/*
flint_printf("dbetas[%wd]: ",j); nmod_poly_print_pretty(I->dbetas + j, "x1"); printf("\n");
*/
    }

    /* set product of betas */
    for (i = r; i >= 0; i--)
    {
        for (j = 0; j < l; j++)
        {
            nmod_mpoly_init(I->prod_mbetas + i*l + j, ctx);
            nmod_mpoly_one(I->prod_mbetas + i*l + j, ctx);
            for (k = 0; k < l; k++)
            {
                if (k == j)
                    continue;
                nmod_mpoly_mul(I->prod_mbetas + i*l + j, I->prod_mbetas + i*l + j, I->mbetas + i*l + k, ctx);
            }
/*
flint_printf("prod_mbetas[%wd][%wd]: ",i,j); nmod_mpoly_print_pretty(I->prod_mbetas + i*l + j, NULL, ctx); printf("\n");
*/
        }        
    }
    for (j = 0; j < l; j++)
    {
        nmod_poly_one(pq);
        for (k = 0; k < l; k++)
        {
            if (k == j)
                continue;
            nmod_poly_mul(pq, pq, I->dbetas + k);
        }
        nmod_poly_init_mod(I->inv_prod_dbetas + j, ctx->ffinfo->mod);
        nmod_poly_xgcd(G, S, I->inv_prod_dbetas + j, I->dbetas + j, pq);
        FLINT_ASSERT(nmod_poly_is_one(G));
/*
flint_printf("inv_prod_dbetas[%wd]: ",j); nmod_poly_print_pretty(I->inv_prod_dbetas + j, "x1"); printf("\n");
*/
    }

    nmod_poly_clear(p);
    nmod_poly_clear(G);
    nmod_poly_clear(S);
    nmod_poly_clear(pq);
/*
printf("_mfactor_disolve_init returning\n");
*/
}

void nmod_mpoly_pfrac_clear(nmod_mpoly_pfrac_t I, const nmod_mpoly_ctx_t ctx)
{
    slong i, j;

    for (j = 0; j < I->l; j++)
    {
        nmod_poly_clear(I->inv_prod_dbetas + j);
        nmod_poly_clear(I->dbetas + j);
        for (i = 0; i <= I->r; i++)
        {
            nmod_mpoly_clear(I->prod_mbetas + i*I->l + j, ctx);
            nmod_mpoly_clear(I->mbetas + i*I->l + j, ctx);
            nmod_mpoly_clear(I->deltas + i*I->l + j, ctx);
        }
    }

    flint_free(I->deltas);
    flint_free(I->inv_prod_dbetas);
    flint_free(I->dbetas);
    flint_free(I->prod_mbetas);
    flint_free(I->mbetas);
}


int nmod_mpoly_pfrac(
    flint_bitcnt_t bits,
    slong r,
    slong num,
    const nmod_mpoly_t t,
    const mp_limb_t * alpha,
    const slong * deg,
    const nmod_mpoly_pfrac_t I,
    const nmod_mpoly_ctx_t ctx)
{
    slong i, j;
    int success;
    nmod_mpoly_struct * deltas = I->deltas + r*I->l;
    nmod_mpoly_struct * newdeltas = I->deltas + (r - 1)*I->l;

    if (r < 1)
    {
        nmod_poly_t dtq, S, R;

        nmod_poly_init_mod(dtq, ctx->ffinfo->mod);
        nmod_poly_init_mod(S, ctx->ffinfo->mod);
        nmod_poly_init_mod(R, ctx->ffinfo->mod);

        _to_poly(dtq, t, ctx);

        success = 1;
        for (i = 0; i < num; i++)
        {
            nmod_poly_mul(S, dtq, I->inv_prod_dbetas + i);
            nmod_poly_rem(R, S, I->dbetas + i);
            _from_poly(deltas + i, bits, R, ctx);
        }

        nmod_poly_clear(dtq);
        nmod_poly_clear(S);
        nmod_poly_clear(R);
    }
    else
    {
        nmod_mpoly_t e, tt, pow, g, q, newt;

FLINT_ASSERT(num == I->l);

        nmod_mpoly_init(e, ctx);
        nmod_mpoly_init(tt, ctx);
        nmod_mpoly_init(pow, ctx);
        nmod_mpoly_init(g, ctx);
        nmod_mpoly_init(q, ctx);
        nmod_mpoly_init(newt, ctx);
        for (i = 0; i < num; i++)
            nmod_mpoly_zero(deltas + i, ctx);

        nmod_mpoly_set(e, t, ctx);

        nmod_mpoly_one(pow, ctx);
        nmod_mpoly_gen(g, r, ctx);
        nmod_mpoly_sub_ui(g, g, alpha[r - 1], ctx);
        for (j = 0; j <= deg[r]; j++)
        {
            success = nmod_mpoly_divides(q, e, pow, ctx);
            FLINT_ASSERT(success);
            nmod_mpoly_evaluate_one_ui(newt, q, r, alpha[r - 1], ctx);
            success = nmod_mpoly_pfrac(bits, r - 1, num, newt, alpha, deg, I, ctx);
            if (!success)
                goto cleanup;

            nmod_mpoly_set(e, t, ctx);
            for (i = 0; i < num; i++)
            {
                nmod_mpoly_mul(tt, newdeltas + i, pow, ctx);
                nmod_mpoly_add(deltas + i, deltas + i, tt, ctx);
                nmod_mpoly_mul(tt, deltas + i, I->prod_mbetas + r*I->l + i, ctx);
                nmod_mpoly_sub(e, e, tt, ctx);
            }

            nmod_mpoly_mul(pow, pow, g, ctx);
        }

        success = nmod_mpoly_is_zero(e, ctx);

cleanup:

        nmod_mpoly_clear(e, ctx);
        nmod_mpoly_clear(tt, ctx);
        nmod_mpoly_clear(pow, ctx);
        nmod_mpoly_clear(g, ctx);
        nmod_mpoly_clear(q, ctx);
        nmod_mpoly_clear(newt, ctx);
    }

    return success;
}

