/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fq_zech_mpoly_factor.h"


void fq_zech_mpoly_evaluate_all_fq_zech_poly(
    fq_zech_poly_t E,
    const fq_zech_mpoly_t A,
    const fq_zech_poly_struct * alphabetas,
    const fq_zech_mpoly_ctx_t ctx)
{
    slong n = ctx->minfo->nvars;
    slong i, N = mpoly_words_per_exp_sp(A->bits, ctx->minfo);
    slong * offsets, * shifts;
    ulong mask = (-UWORD(1)) >> (FLINT_BITS - A->bits);
    slong * starts, * ends, * stops;
    ulong * es;
    fq_zech_poly_struct * realE;

    if (A->length < 1)
    {
        fq_zech_poly_zero(E, ctx->fqctx);
        return;
    }

    starts = FLINT_ARRAY_ALLOC(n, slong);
    ends   = FLINT_ARRAY_ALLOC(n, slong);
    stops  = FLINT_ARRAY_ALLOC(n, slong);
    es     = FLINT_ARRAY_ALLOC(n, ulong);
    realE  = FLINT_ARRAY_ALLOC(n + 1, fq_zech_poly_struct);
    for (i = 0; i < n + 1; i++)
        fq_zech_poly_init(realE + i, ctx->fqctx);

    offsets = FLINT_ARRAY_ALLOC(n, slong);
    shifts  = FLINT_ARRAY_ALLOC(n, slong);
    for (i = 0; i < ctx->minfo->nvars; i++)
        mpoly_gen_offset_shift_sp(offsets + i, shifts + i, i, A->bits, ctx->minfo);

    _fq_zech_mpoly_evaluate_rest_n_poly_fq(realE, starts, ends, stops, es,
                    A->coeffs + 0, A->exps, A->length, 0,
                                        alphabetas, offsets, shifts, N, mask,
                                                ctx->minfo->nvars, ctx->fqctx);
    fq_zech_poly_swap(E, realE + 0, ctx->fqctx);

    for (i = 0; i < n + 1; i++)
        fq_zech_poly_clear(realE + i, ctx->fqctx);
    flint_free(realE);
    flint_free(starts);
    flint_free(ends);
    flint_free(stops);
    flint_free(es);

    flint_free(offsets);
    flint_free(shifts);
}


int fq_zech_mpoly_factor_lcc_wang(
    fq_zech_mpoly_struct * lc_divs,
    const fq_zech_mpoly_factor_t lcAfac,
    const fq_zech_poly_t Auc,
    const fq_zech_bpoly_struct * Auf,
    slong r,
    const fq_zech_poly_struct * alpha,
    const fq_zech_mpoly_ctx_t ctx)
{
    int success;
    slong i, j, k;
    const slong n = ctx->minfo->nvars - 1;
    fq_zech_poly_struct * lcAfaceval;
    fq_zech_poly_struct * salpha;
    fq_zech_poly_struct * d;
    fq_zech_poly_t Q, R;
    fq_zech_mpoly_t t;

    fq_zech_poly_init(Q, ctx->fqctx);
    fq_zech_poly_init(R, ctx->fqctx);
    fq_zech_mpoly_init(t, ctx);

    salpha = FLINT_ARRAY_ALLOC((n + 1), fq_zech_poly_struct);
    fq_zech_poly_init(salpha + 0, ctx->fqctx);
    for (i = 0; i < n; i++)
        salpha[i + 1] = alpha[i];

    lcAfaceval = FLINT_ARRAY_ALLOC(lcAfac->num, fq_zech_poly_struct);
    for (i = 0; i < lcAfac->num; i++)
        fq_zech_poly_init(lcAfaceval + i, ctx->fqctx);

    d = FLINT_ARRAY_ALLOC(lcAfac->num + 1, fq_zech_poly_struct);
    for (i = 0; i < lcAfac->num + 1; i++)
        fq_zech_poly_init(d + i, ctx->fqctx);

    /* init done */

    for (j = 0; j < lcAfac->num; j++)
        fq_zech_mpoly_evaluate_all_fq_zech_poly(lcAfaceval + j, lcAfac->poly + j, salpha, ctx);

    fq_zech_poly_set(d + 0, Auc, ctx->fqctx);
    for (i = 0; i < lcAfac->num; i++)
    {
        fq_zech_poly_make_monic(Q, lcAfaceval + i, ctx->fqctx);
        if (fq_zech_poly_degree(Q, ctx->fqctx) < 1)
        {
            success = 0;
            goto cleanup;
        }
        for (j = i; j >= 0; j--)
        {
            fq_zech_poly_set(R, d + j, ctx->fqctx);
            while (fq_zech_poly_degree(R, ctx->fqctx) > 0)
            {
                fq_zech_poly_gcd(R, R, Q, ctx->fqctx);
                fq_zech_poly_divrem(Q, salpha + 0, Q, R, ctx->fqctx);
                FLINT_ASSERT(fq_zech_poly_is_zero(salpha + 0, ctx->fqctx));
                if (fq_zech_poly_degree(Q, ctx->fqctx) < 1)
                {
                    success = 0;
                    goto cleanup;
                }
            }
        }
        fq_zech_poly_set(d + i + 1, Q, ctx->fqctx);
    }

    for (j = 0; j < r; j++)
    {
        fq_zech_mpoly_one(lc_divs + j, ctx);
        fq_zech_poly_mul(R, Auf[j].coeffs + Auf[j].length - 1, Auc, ctx->fqctx);
        for (i = lcAfac->num - 1; i >= 0; i--)
        {
            fq_zech_poly_make_monic(Q, lcAfaceval + i, ctx->fqctx);
            if (fq_zech_poly_degree(Q, ctx->fqctx) < 1)
                continue;
            k = fq_zech_poly_remove(R, Q, ctx->fqctx);
            fq_zech_mpoly_pow_ui(t, lcAfac->poly + i, k, ctx);
            fq_zech_mpoly_mul(lc_divs + j, lc_divs + j, t, ctx);
        }
    }

    success = 1;

cleanup:

    fq_zech_poly_clear(Q, ctx->fqctx);
    fq_zech_poly_clear(R, ctx->fqctx);
    fq_zech_mpoly_clear(t, ctx);

    for (i = 0; i < lcAfac->num; i++)
        fq_zech_poly_clear(lcAfaceval + i, ctx->fqctx);
    flint_free(lcAfaceval);

    for (i = 0; i < lcAfac->num + 1; i++)
        fq_zech_poly_clear(d + i, ctx->fqctx);
    flint_free(d);

    fq_zech_poly_clear(salpha + 0, ctx->fqctx);
    flint_free(salpha);

    return success;
}
