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
    fmpz_mpoly_factor_t qfac,
    const fmpz_mpoly_t q,
    const fmpz_mpoly_factor_t pfac,
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

    FLINT_ASSERT(pfac->num > 1);

    newdeg = (slong *) flint_malloc((n + 1)*sizeof(slong));
    fmpz_mpoly_init(lcq, ctx);
    fmpz_mpoly_init(lcp, ctx);
    fmpz_mpoly_init(t, ctx);
    fmpz_mpoly_init(newq, ctx);
    fmpz_mpoly_univar_init(u, ctx);

    FLINT_ASSERT(fmpz_is_one(pfac->constant));
    FLINT_ASSERT(fmpz_mpoly_factor_matches(p, pfac, ctx));

    _fmpz_mpoly_get_lead0(lcq, q, ctx);
    fmpz_mpoly_evaluate_one_fmpz(lcp, lcq, m, alpha + m - 1, ctx);

    FLINT_ASSERT(lcp->length > 0);

    fmpz_mpoly_pow_ui(t, lcq, pfac->num - 1, ctx);
    fmpz_mpoly_mul(newq, q, t, ctx);
    fmpz_mpoly_degrees_si(newdeg, newq, ctx);

    fmpz_set(qfac->constant, pfac->constant);
    fmpz_mpoly_factor_fit_length(qfac, pfac->num, ctx);
    qfac->num = pfac->num;
    for (i = 0; i < pfac->num; i++)
    {
        _fmpz_mpoly_get_lead0(t, pfac->poly + i, ctx);
        success = fmpz_mpoly_divides(t, lcp, t, ctx);
        FLINT_ASSERT(success);
        fmpz_mpoly_mul(qfac->poly + i, pfac->poly + i, t, ctx);
        _fmpz_mpoly_set_lead0(qfac->poly + i, qfac->poly + i, lcq, ctx);
        fmpz_one(qfac->exp + i);
    }

    success = fmpz_mpoly_hlift(m, qfac->poly, qfac->num,
                                                     alpha, newq, newdeg, ctx);
    if (!success)
        goto cleanup;

    for (i = 0; i < qfac->num; i++)
    {
        fmpz_mpoly_to_univar(u, qfac->poly + i, 0, ctx);
        success = fmpz_mpoly_univar_content_mpoly(t, u, ctx);
        if (!success)
        {
            success = -1;
            goto cleanup;
        }
        success = fmpz_mpoly_divides(qfac->poly + i, qfac->poly + i, t, ctx);
        FLINT_ASSERT(success);
        FLINT_ASSERT(qfac->poly[i].length > 0);
        if (fmpz_sgn(qfac->poly[i].coeffs + 0) < 0)
            fmpz_mpoly_neg(qfac->poly + i, qfac->poly + i, ctx);
    }

    success = 1;

cleanup:

    flint_free(newdeg);
    fmpz_mpoly_clear(lcq, ctx);
    fmpz_mpoly_clear(lcp, ctx);
    fmpz_mpoly_clear(t, ctx);
    fmpz_mpoly_clear(newq, ctx);
    fmpz_mpoly_univar_clear(u, ctx);

    /* q and its factors are primitive with positive lc */
    FLINT_ASSERT(!success || (fmpz_is_one(qfac->constant) &&
                              fmpz_mpoly_factor_matches(q, qfac, ctx)));
    return success;
}


/* A is square free and primitive w.r.t all variables */
int fmpz_mpoly_factor_irred_default(
    fmpz_mpoly_factor_t fac,
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
    fmpz_mpoly_factor_t qfac, pfac, tfac, dfac;
    fmpz_mpoly_t t, p, q;
    fmpz_mpoly_univar_t u;

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
    fmpz_mpoly_factor_init(pfac, ctx);
    fmpz_mpoly_factor_init(qfac, ctx);
    fmpz_mpoly_factor_init(tfac, ctx);
    fmpz_mpoly_factor_init(dfac, ctx);
    fmpz_mpoly_init(t, ctx);
    fmpz_mpoly_init(p, ctx);
    fmpz_mpoly_init(q, ctx);
    fmpz_mpoly_univar_init(u, ctx);

    fmpz_mpoly_factor_one(fac, ctx);

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

/*
usleep(1000000);
printf("alpha = "); tuple_print(alpha, n);
*/

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

    /* make sure our univar is squarefree */
    fmpz_mpoly_derivative(t, Aevals + 0, 0, ctx);
    fmpz_mpoly_gcd(t, t, Aevals + 0, ctx);
    if (!fmpz_mpoly_is_fmpz(t, ctx))
        goto next_alpha;

    /* make our evaluations primitive */
    for (i = n - 1; i > 0; i--)
    {
        fmpz_mpoly_to_univar(u, Aevals + i, 0, ctx);
        fmpz_mpoly_univar_content_mpoly(t, u, ctx);
        success = fmpz_mpoly_divides(Aevals + i, Aevals + i, t, ctx);
        FLINT_ASSERT(success);
        FLINT_ASSERT(Aevals[i].length > 0);
        if (fmpz_sgn(Aevals[i].coeffs + 0) < 0)
            fmpz_mpoly_neg(Aevals + i, Aevals + i, ctx);
    }

    fmpz_mpoly_factor_irred_bivar(pfac, Aevals + 1, 0, 1, ctx);

    for (m = 2; m <= n; m++)
    {
        fmpz_mpoly_set(q, m < n ? Aevals + m : A, ctx);
        fmpz_mpoly_set(p, Aevals + m - 1, ctx);

        FLINT_ASSERT(fmpz_is_one(pfac->constant));
        FLINT_ASSERT(fmpz_mpoly_factor_matches(p, pfac, ctx));

        /* if p has only one factor, A must be irreducible */
        if (pfac->num < 2)
        {
            fmpz_mpoly_factor_append_ui(fac, A, 1, ctx);
            success = 1;
            goto cleanup;
        }

        success = _try_lift(qfac, q, pfac, p, m, alpha, n, ctx);
        if (success > 0)
        {
            fmpz_mpoly_factor_swap(qfac, pfac, ctx);
            continue;
        }
        else if (success < 0)
        {
            success = 0;
            goto cleanup;
        }

        /* if we couldn't lift two local factors, A must be irreducible */
        if (pfac->num == 2)
        {
            fmpz_mpoly_factor_append_ui(fac, A, 1, ctx);
            success = 1;
            goto cleanup;
        }

        qfac->num = 0;

try_again:

        for (r = 1; r <= pfac->num/2; r++)
        {
            subset_first(subset, pfac->num, r);

            FLINT_ASSERT(fmpz_is_one(pfac->constant));
            FLINT_ASSERT(fmpz_mpoly_factor_matches(p, pfac, ctx));

            do {
                fmpz_mpoly_factor_fit_length(dfac, 2, ctx);
                dfac->num = 2;
                fmpz_one(dfac->constant);
                fmpz_mpoly_one(dfac->poly + 0, ctx);
                fmpz_mpoly_one(dfac->poly + 1, ctx);
                fmpz_one(dfac->exp + 0);
                fmpz_one(dfac->exp + 1);
                for (i = 0; i < pfac->num; i++)
                {
                    j = fmpz_tstbit(subset, i);
                    fmpz_mpoly_mul(dfac->poly + j, dfac->poly + j,
                                                          pfac->poly + i, ctx);
                }

                success = _try_lift(tfac, q, dfac, p, m, alpha, n, ctx);
                if (success > 0)
                {
                    for (i = pfac->num - 1; i >= 0; i--)
                    {
                        if (fmpz_tstbit(subset, i))
                        {
                            fmpz_mpoly_swap(pfac->poly + i,
                                              pfac->poly + pfac->num - 1, ctx);
                            pfac->num--;
                        }
                    }
                    fmpz_mpoly_factor_append_ui(qfac, tfac->poly + 1, 1, ctx);
                    fmpz_mpoly_swap(q, tfac->poly + 0, ctx);
                    fmpz_mpoly_swap(p, dfac->poly + 0, ctx);
                    goto try_again;
                }
                else if (success < 0)
                {
                    success = 0;
                    goto cleanup;
                }
            }
            while (subset_next(subset, subset, pfac->num));
        }
        /* if pfac could not be combined, p must be irreducible */
        fmpz_mpoly_factor_append_ui(qfac, q, 1, ctx);
        fmpz_mpoly_factor_swap(qfac, pfac, ctx);
    }

    success = 1;

    for (i = 0; i < pfac->num; i++)
        fmpz_mpoly_factor_append_ui(fac, pfac->poly + i, 1, ctx);

cleanup:

    fmpz_clear(subset);
    _fmpz_vec_clear(alphait, n);
    _fmpz_vec_clear(alpha, n);
    for (i = 0; i < n; i++)
        fmpz_mpoly_clear(Aevals + i, ctx);
    flint_free(Aevals);
    flint_free(deg);
    flint_free(degeval);
    fmpz_mpoly_factor_clear(pfac, ctx);
    fmpz_mpoly_factor_clear(qfac, ctx);
    fmpz_mpoly_factor_clear(tfac, ctx);
    fmpz_mpoly_factor_clear(dfac, ctx);
    fmpz_mpoly_clear(t, ctx);
    fmpz_mpoly_clear(p, ctx);
    fmpz_mpoly_clear(q, ctx);
    fmpz_mpoly_univar_clear(u, ctx);

    FLINT_ASSERT(!success || (fmpz_is_one(fac->constant) &&
                              fmpz_mpoly_factor_matches(A, fac, ctx)));

    return success;
}
