/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fmpz_mpoly_factor.h"


int fmpz_mpoly_factor_lcc_wang(
    fmpz_mpolyv_t lc_divs,
    const fmpz_mpoly_factor_t lcAfac,
    const fmpz_t Auc,
    const fmpz_poly_struct * Auf,
    slong r,
    const fmpz * alpha,
    const fmpz_mpoly_ctx_t ctx)
{
    int success;
    slong i, j, k;
    const slong n = ctx->minfo->nvars - 1;
    fmpz * lcAfaceval = _fmpz_vec_init(lcAfac->num);
    const fmpz ** salpha = (const fmpz **) flint_malloc((n + 1)*sizeof(fmpz *));
    fmpz * d = _fmpz_vec_init(1 + lcAfac->num);
    fmpz zero = 0;
    fmpz * dtilde = _fmpz_vec_init(r);
    fmpz_t delta, Q, R;
    fmpz_mpoly_t t;

    fmpz_init(delta);
    fmpz_init(Q);
    fmpz_init(R);

    fmpz_mpoly_init(t, ctx);

    salpha[0] = &zero;
    for (i = 0; i < n; i++)
        salpha[i + 1] = alpha + i;

    for (j = 0; j < lcAfac->num; j++)
        fmpz_mpoly_evaluate_all_fmpz(lcAfaceval + j, lcAfac->poly + j,
                                                 (fmpz * const *) salpha, ctx);

    fmpz_mul(d + 0, Auc, lcAfac->constant);
    for (i = 0; i < lcAfac->num; i++)
    {
        fmpz_abs(Q, lcAfaceval + i);
        if (fmpz_is_one(Q) || fmpz_is_zero(Q))
        {
            success = 0;
            goto cleanup;
        }
        for (j = i; j >= 0; j--)
        {
            fmpz_set(R, d + j);
            while (!fmpz_is_one(R))
            {
                fmpz_gcd(R, R, Q);
                fmpz_divexact(Q, Q, R);
                if (fmpz_is_one(Q))
                {
                    success = 0;
                    goto cleanup;
                }
            }
        }
        fmpz_set(d + i + 1, Q);
    }

    fmpz_mpolyv_fit_length(lc_divs, r, ctx);
    lc_divs->length = r;

    for (j = 0; j < r; j++)
    {
        fmpz_mpoly_one(lc_divs->coeffs + j, ctx);
        fmpz_mul(R, Auf[j].coeffs + Auf[j].length - 1, Auc);
        fmpz_one(dtilde + j);
        for (i = lcAfac->num - 1; i >= 0; i--)
        {
            fmpz_abs(Q, lcAfaceval + i);
            if (fmpz_cmp_ui(Q, 2) < 0)
                continue;
            k = fmpz_remove(R, R, Q);
            fmpz_mpoly_pow_ui(t, lcAfac->poly + i, k, ctx);
            fmpz_mpoly_mul(lc_divs->coeffs + j, lc_divs->coeffs + j, t, ctx);
            fmpz_pow_ui(Q, lcAfaceval + i, k);
            fmpz_mul(dtilde + j, dtilde + j, Q);
        }
    }

    fmpz_set(delta, Auc);
    for (j = 0; j < r; j++)
    {
        FLINT_ASSERT(Auf[j].length > 0);
        fmpz_gcd(R, Auf[j].coeffs + Auf[j].length - 1, dtilde + j);
        FLINT_ASSERT(fmpz_divisible(Auf[j].coeffs + Auf[j].length - 1, R));
        fmpz_divexact(Q, Auf[j].coeffs + Auf[j].length - 1, R);
        fmpz_mpoly_scalar_mul_fmpz(lc_divs->coeffs + j, lc_divs->coeffs + j, Q, ctx);
    }

    success = 1;

cleanup:

    fmpz_clear(delta);
    fmpz_clear(Q);
    fmpz_clear(R);
    fmpz_mpoly_clear(t, ctx);
    _fmpz_vec_clear(lcAfaceval, lcAfac->num);
    _fmpz_vec_clear(d, 1 + lcAfac->num);
    _fmpz_vec_clear(dtilde, r);
    flint_free(salpha);

    return success;
}
