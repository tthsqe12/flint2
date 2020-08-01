/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fq_nmod_mpoly_factor.h"
#include "ui_factor.h"


int fq_nmod_mpoly_factor_lcc_wang(
    fq_nmod_mpoly_struct * lc_divs,
    const fq_nmod_mpoly_factor_t lcAfac,
    const fq_nmod_poly_t Auc,
    const fq_nmod_bpoly_struct * Auf,
    slong r,
    const fq_nmod_poly_struct * alpha,
    const fq_nmod_mpoly_ctx_t ctx)
{
    int success;
    slong i, j, k;
    const slong n = ctx->minfo->nvars - 1;
    fq_nmod_poly_struct * lcAfaceval;
    fq_nmod_poly_struct ** salpha;
    fq_nmod_poly_struct * d;
    fq_nmod_poly_t Q, R, T;
    fq_nmod_mpoly_t t;

    fq_nmod_poly_init(Q, ctx->fqctx);
    fq_nmod_poly_init(R, ctx->fqctx);
    fq_nmod_poly_init(T, ctx->fqctx);

    salpha = FLINT_ARRAY_ALLOC((n + 1), fq_nmod_poly_struct *);
    salpha[0] = T;
    for (i = 0; i < n; i++)
        salpha[i + 1] = (fq_nmod_poly_struct *) alpha + i;

    lcAfaceval = FLINT_ARRAY_ALLOC(lcAfac->num, fq_nmod_poly_struct);
    for (i = 0; i < lcAfac->num; i++)
        fq_nmod_poly_init(lcAfaceval + i, ctx->fqctx);

    d = FLINT_ARRAY_ALLOC(lcAfac->num + 1, fq_nmod_poly_struct);
    for (i = 0; i < lcAfac->num + 1; i++)
        fq_nmod_poly_init(d + i, ctx->fqctx);

    fq_nmod_mpoly_init(t, ctx);

    for (j = 0; j < lcAfac->num; j++)
        fq_nmod_mpoly_compose_fq_nmod_poly(lcAfaceval + j, lcAfac->poly + j, salpha, ctx);

    fq_nmod_poly_set(d + 0, Auc, ctx->fqctx);
    for (i = 0; i < lcAfac->num; i++)
    {
        fq_nmod_poly_make_monic(Q, lcAfaceval + i, ctx->fqctx);
        if (fq_nmod_poly_degree(Q, ctx->fqctx) < 1)
        {
            success = 0;
            goto cleanup;
        }
        for (j = i; j >= 0; j--)
        {
            fq_nmod_poly_set(R, d + j, ctx->fqctx);
            while (fq_nmod_poly_degree(R, ctx->fqctx) > 0)
            {
                fq_nmod_poly_gcd(R, R, Q, ctx->fqctx);
                fq_nmod_poly_divrem(Q, T, Q, R, ctx->fqctx);
                if (fq_nmod_poly_degree(Q, ctx->fqctx) < 1)
                {
                    success = 0;
                    goto cleanup;
                }
            }
        }
        fq_nmod_poly_set(d + i + 1, Q, ctx->fqctx);
    }

    for (j = 0; j < r; j++)
    {
        fq_nmod_mpoly_one(lc_divs + j, ctx);
        fq_nmod_poly_mul(R, Auf[j].coeffs + Auf[j].length - 1, Auc, ctx->fqctx);
        for (i = lcAfac->num - 1; i >= 0; i--)
        {
            fq_nmod_poly_make_monic(Q, lcAfaceval + i, ctx->fqctx);
            if (fq_nmod_poly_degree(Q, ctx->fqctx) < 1)
                continue;
            k = fq_nmod_poly_remove(R, Q, ctx->fqctx);
            fq_nmod_mpoly_pow_ui(t, lcAfac->poly + i, k, ctx);
            fq_nmod_mpoly_mul(lc_divs + j, lc_divs + j, t, ctx);
        }
    }

    success = 1;

cleanup:

    fq_nmod_poly_clear(Q, ctx->fqctx);
    fq_nmod_poly_clear(R, ctx->fqctx);
    fq_nmod_poly_clear(T, ctx->fqctx);
    fq_nmod_mpoly_clear(t, ctx);

    for (i = 0; i < lcAfac->num; i++)
        fq_nmod_poly_clear(lcAfaceval + i, ctx->fqctx);
    flint_free(lcAfaceval);

    for (i = 0; i < lcAfac->num + 1; i++)
        fq_nmod_poly_clear(d + i, ctx->fqctx);
    flint_free(d);

    flint_free(salpha);

    return success;
}
