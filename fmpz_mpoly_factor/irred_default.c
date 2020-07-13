/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fmpz_mpoly_factor.h"
#include "nmod_mpoly_factor.h"


static int _try_lift(
    fmpz_mpolyv_t qfac,
    const fmpz_mpoly_t q,
    const fmpz_mpolyv_t pfac,
    const fmpz_mpoly_t p,
    slong m,
    fmpz * alpha,
    slong n,
    const fmpz_mpoly_ctx_t ctx)
{
    int success;
    slong i;
    slong * newdeg;
    fmpz_mpoly_t lcq, lcp, t, newq;
    fmpz_mpoly_univar_t u;

    FLINT_ASSERT(pfac->length > 1);

    newdeg = (slong *) flint_malloc((n + 1)*sizeof(slong));
    fmpz_mpoly_init(lcq, ctx);
    fmpz_mpoly_init(lcp, ctx);
    fmpz_mpoly_init(t, ctx);
    fmpz_mpoly_init(newq, ctx);
    fmpz_mpoly_univar_init(u, ctx);

#if WANT_ASSERT
    fmpz_mpoly_one(t, ctx);
    for (i = 0; i < pfac->length; i++)
        fmpz_mpoly_mul(t, t, pfac->coeffs + i, ctx);
    FLINT_ASSERT(fmpz_mpoly_equal(t, p, ctx));
#endif

    _fmpz_mpoly_get_lead0(lcq, q, ctx);
    fmpz_mpoly_evaluate_one_fmpz(lcp, lcq, m, alpha + m - 1, ctx);

    FLINT_ASSERT(lcp->length > 0);

    fmpz_mpoly_pow_ui(t, lcq, pfac->length - 1, ctx);
    fmpz_mpoly_mul(newq, q, t, ctx);

    if (newq->bits > FLINT_BITS)
    {
        success = -1;
        goto cleanup;
    }

    fmpz_mpoly_degrees_si(newdeg, newq, ctx);

    fmpz_mpolyv_fit_length(qfac, pfac->length, ctx);
    qfac->length = pfac->length;
    for (i = 0; i < pfac->length; i++)
    {
        _fmpz_mpoly_get_lead0(t, pfac->coeffs + i, ctx);
        success = fmpz_mpoly_divides(t, lcp, t, ctx);
        FLINT_ASSERT(success);
        fmpz_mpoly_mul(qfac->coeffs + i, pfac->coeffs + i, t, ctx);
        _fmpz_mpoly_set_lead0(qfac->coeffs + i, qfac->coeffs + i, lcq, ctx);
    }

    success = fmpz_mpoly_hlift(m, qfac->coeffs, qfac->length,
                                                     alpha, newq, newdeg, ctx);
    if (!success)
        goto cleanup;

    for (i = 0; i < qfac->length; i++)
    {
        fmpz_mpoly_to_univar(u, qfac->coeffs + i, 0, ctx);
        success = fmpz_mpoly_univar_content_mpoly(t, u, ctx);
        if (!success)
        {
            success = -1;
            goto cleanup;
        }
        success = fmpz_mpoly_divides(qfac->coeffs + i,
                                     qfac->coeffs + i, t, ctx);
        FLINT_ASSERT(success);
        fmpz_mpoly_unit_normalize(qfac->coeffs + i, ctx);
    }

    success = 1;

cleanup:

    flint_free(newdeg);
    fmpz_mpoly_clear(lcq, ctx);
    fmpz_mpoly_clear(lcp, ctx);
    fmpz_mpoly_clear(t, ctx);
    fmpz_mpoly_clear(newq, ctx);
    fmpz_mpoly_univar_clear(u, ctx);

#if WANT_ASSERT
    if (success)
    {
        fmpz_mpoly_init(t, ctx);
        fmpz_mpoly_one(t, ctx);
        for (i = 0; i < qfac->length; i++)
            fmpz_mpoly_mul(t, t, qfac->coeffs + i, ctx);
        FLINT_ASSERT(fmpz_mpoly_equal(t, q, ctx));
        fmpz_mpoly_clear(t, ctx);
    }
#endif

    return success;
}


/* A is square free and primitive w.r.t all variables */
int fmpz_mpoly_factor_irred_default(
    fmpz_mpolyv_t fac,
    const fmpz_mpoly_t A,
    const fmpz_mpoly_ctx_t ctx)
{
    int success;
    const slong n = ctx->minfo->nvars - 1;
    slong i, j, m, r;
    fmpz_t subset;
    fmpz * alpha, * alphait;
    fmpz_mpoly_struct * Aevals;
    slong * deg, * degeval;
    fmpz_mpolyv_t qfac, pfac, tfac, dfac;
    fmpz_mpoly_t t, p, q;
    fmpz_mpoly_univar_t u;
    fmpz_poly_t c;
    fmpz_bpoly_t B;
    fmpz_tpoly_t F;

    FLINT_ASSERT(n > 1);
    FLINT_ASSERT(A->length > 0);
    FLINT_ASSERT(fmpz_sgn(A->coeffs + 0) > 0);

    fmpz_init(subset);
    alphait = _fmpz_vec_init(n);
    alpha   = _fmpz_vec_init(n);
    Aevals  = (fmpz_mpoly_struct *) flint_malloc(n*sizeof(fmpz_mpoly_struct));
    deg     = (slong *) flint_malloc((n + 1)*sizeof(slong));
    degeval = (slong *) flint_malloc((n + 1)*sizeof(slong));
    for (i = 0; i < n; i++)
        fmpz_mpoly_init(Aevals + i, ctx);
    fmpz_mpolyv_init(pfac, ctx);
    fmpz_mpolyv_init(qfac, ctx);
    fmpz_mpolyv_init(tfac, ctx);
    fmpz_mpolyv_init(dfac, ctx);
    fmpz_mpoly_init(t, ctx);
    fmpz_mpoly_init(p, ctx);
    fmpz_mpoly_init(q, ctx);
    fmpz_mpoly_univar_init(u, ctx);
    fmpz_poly_init(c);
    fmpz_bpoly_init(B);
    fmpz_tpoly_init(F);

    fmpz_mpoly_degrees_si(deg, A, ctx);
    goto got_alpha;

next_alpha:

    tuple_next(alphait, n);
    for (i = 0; i < n; i++)
    {
        j = n - 1 - i;
        fmpz_cdiv_q_2exp(alpha + j, alphait + i, 1);
        if (fmpz_is_even(alphait + i))
            fmpz_neg(alpha + j, alpha + j);
    }

got_alpha:

    /* ensure degrees do not drop under evaluation */
    for (i = n - 1; i >= 0; i--)
    {
        fmpz_mpoly_evaluate_one_fmpz(Aevals + i,
                       i == n - 1 ? A : Aevals + i + 1, i + 1, alpha + i, ctx);
        fmpz_mpoly_degrees_si(degeval, Aevals + i, ctx);
        for (j = 0; j <= i; j++)
        {
            if (degeval[j] != deg[j])
            {
                tuple_saturate(alphait, n, n - i);
                goto next_alpha;
            }
        }
    }

    /* make sure univar is squarefree */
    fmpz_mpoly_derivative(t, Aevals + 0, 0, ctx);
    fmpz_mpoly_gcd(t, t, Aevals + 0, ctx);
    if (!fmpz_mpoly_is_fmpz(t, ctx))
        goto next_alpha;

    /* make evaluations primitive */
    for (i = n - 1; i > 0; i--)
    {
        fmpz_mpoly_to_univar(u, Aevals + i, 0, ctx);
        success = fmpz_mpoly_univar_content_mpoly(t, u, ctx);
        if (!success)
            goto cleanup;
        success = fmpz_mpoly_divides(Aevals + i, Aevals + i, t, ctx);
        FLINT_ASSERT(success);
        fmpz_mpoly_unit_normalize(Aevals + i, ctx);
    }

    fmpz_mpoly_get_bpoly(B, Aevals + 1, 0, 1, ctx);
    fmpz_bpoly_factor(c, F, B);

    /* c = +-1 */
    FLINT_ASSERT(c->length == 1 && fmpz_is_pm1(c->coeffs + 0));

    fmpz_mpolyv_fit_length(pfac, F->length, ctx);
    pfac->length = F->length;
    for (i = 0; i < F->length; i++)
    {
        fmpz_mpoly_set_fmpz_bpoly(pfac->coeffs + i, A->bits,
                                                 F->coeffs + i, 0, 1, ctx);
        fmpz_mpoly_unit_normalize(pfac->coeffs + i, ctx);
    }

    for (m = 2; m <= n; m++)
    {
        fmpz_mpoly_set(q, m < n ? Aevals + m : A, ctx);
        fmpz_mpoly_set(p, Aevals + m - 1, ctx);

    #if WANT_ASSERT
        fmpz_mpoly_one(t, ctx);
        for (i = 0; i < pfac->length; i++)
            fmpz_mpoly_mul(t, t, pfac->coeffs + i, ctx);
        FLINT_ASSERT(fmpz_mpoly_equal(t, p, ctx));
    #endif

        /* if p has only one factor, A must be irreducible */
        if (pfac->length < 2)
        {
            fmpz_mpolyv_fit_length(fac, 1, ctx);
            fac->length = 1;
            fmpz_mpoly_set(fac->coeffs + 0, A, ctx);
            success = 1;
            goto cleanup;
        }

        success = _try_lift(qfac, q, pfac, p, m, alpha, n, ctx);
        if (success > 0)
        {
            fmpz_mpolyv_swap(qfac, pfac, ctx);
            continue;
        }
        else if (success < 0)
        {
            success = 0;
            goto cleanup;
        }

        /* if we couldn't lift two local factors, A must be irreducible */
        if (pfac->length == 2)
        {
            fmpz_mpolyv_fit_length(fac, 1, ctx);
            fac->length = 1;
            fmpz_mpoly_set(fac->coeffs + 0, A, ctx);
            success = 1;
            goto cleanup;
        }

        qfac->length = 0;

try_again:

        for (r = 1; r <= pfac->length/2; r++)
        {
            subset_first(subset, pfac->length, r);

        #if WANT_ASSERT
            fmpz_mpoly_one(t, ctx);
            for (i = 0; i < pfac->length; i++)
                fmpz_mpoly_mul(t, t, pfac->coeffs + i, ctx);
            FLINT_ASSERT(fmpz_mpoly_equal(t, p, ctx));
        #endif

            do {
                fmpz_mpolyv_fit_length(dfac, 2, ctx);
                dfac->length = 2;
                fmpz_mpoly_one(dfac->coeffs + 0, ctx);
                fmpz_mpoly_one(dfac->coeffs + 1, ctx);
                for (i = 0; i < pfac->length; i++)
                {
                    j = fmpz_tstbit(subset, i);
                    fmpz_mpoly_mul(dfac->coeffs + j,
                                   dfac->coeffs + j, pfac->coeffs + i, ctx);
                }

                success = _try_lift(tfac, q, dfac, p, m, alpha, n, ctx);
                if (success > 0)
                {
                    for (i = pfac->length - 1; i >= 0; i--)
                    {
                        if (fmpz_tstbit(subset, i))
                        {
                            pfac->length--;
                            fmpz_mpoly_swap(pfac->coeffs + i,
                                            pfac->coeffs + pfac->length, ctx);
                        }
                    }

                    fmpz_mpolyv_fit_length(qfac, qfac->length + 1, ctx);
                    fmpz_mpoly_swap(qfac->coeffs + qfac->length,
                                                        tfac->coeffs + 1, ctx);
                    qfac->length++;

                    fmpz_mpoly_swap(q, tfac->coeffs + 0, ctx);
                    fmpz_mpoly_swap(p, dfac->coeffs + 0, ctx);
                    goto try_again;
                }
                else if (success < 0)
                {
                    success = 0;
                    goto cleanup;
                }
            }
            while (subset_next(subset, subset, pfac->length));
        }

        /* if pfac could not be combined, q must be irreducible */
        fmpz_mpolyv_fit_length(qfac, qfac->length + 1, ctx);
        fmpz_mpoly_swap(qfac->coeffs + qfac->length, q, ctx);
        qfac->length++;
        fmpz_mpolyv_swap(qfac, pfac, ctx);
    }

    success = 1;

    fmpz_mpolyv_swap(fac, pfac, ctx);

cleanup:

    fmpz_clear(subset);
    _fmpz_vec_clear(alphait, n);
    _fmpz_vec_clear(alpha, n);
    for (i = 0; i < n; i++)
        fmpz_mpoly_clear(Aevals + i, ctx);
    flint_free(Aevals);
    flint_free(deg);
    flint_free(degeval);
    fmpz_mpolyv_clear(pfac, ctx);
    fmpz_mpolyv_clear(qfac, ctx);
    fmpz_mpolyv_clear(tfac, ctx);
    fmpz_mpolyv_clear(dfac, ctx);
    fmpz_mpoly_clear(t, ctx);
    fmpz_mpoly_clear(p, ctx);
    fmpz_mpoly_clear(q, ctx);
    fmpz_mpoly_univar_clear(u, ctx);
    fmpz_poly_clear(c);
    fmpz_bpoly_clear(B);
    fmpz_tpoly_clear(F);

#if WANT_ASSERT
    if (success)
    {
        fmpz_mpoly_init(t, ctx);
        fmpz_mpoly_one(t, ctx);
        for (i = 0; i < fac->length; i++)
            fmpz_mpoly_mul(t, t, fac->coeffs + i, ctx);
        FLINT_ASSERT(fmpz_mpoly_equal(t, A, ctx));
        fmpz_mpoly_clear(t, ctx);
    }
#endif

    return success;
}
