/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fmpz_mpoly_factor.h"


int fmpz_mpoly_factor_irred_wang(
    fmpz_mpolyv_t fac,
    const fmpz_mpoly_t A,
    const fmpz_mpoly_factor_t lcAfac,
    const fmpz_mpoly_t lcA,
    const fmpz_mpoly_ctx_t ctx)
{
    int success;
    const slong n = ctx->minfo->nvars - 1;
    slong i, j, k;
    fmpz * alpha;
    fmpz_mpoly_struct * Aevals;
    slong * degs, * degeval;
    fmpz_mpolyv_t tfac;
    fmpz_mpoly_t t, Acopy;
    fmpz_mpoly_struct * newA;
    fmpz_poly_t Au;
    fmpz_poly_factor_t Aufac;
    slong alpha_modulus, alpha_count;
    flint_rand_t state;
    fmpz_mpoly_t m, mpow;
    fmpz_mpolyv_t new_lcs, lc_divs;
    fmpz_t q;
    slong r;

    FLINT_ASSERT(n > 1);
    FLINT_ASSERT(ctx->minfo->ord == ORD_LEX);
    FLINT_ASSERT(A->bits <= FLINT_BITS);
    FLINT_ASSERT(fmpz_mpoly_factor_matches(lcA, lcAfac, ctx));

/*
flint_printf("_try_wang called n = %wd\n", n);
flint_printf("A: "); fmpz_mpoly_print_pretty(A, NULL, ctx); flint_printf("\n");
*/
    fmpz_init(q);

    flint_randinit(state);

    fmpz_mpoly_init(Acopy, ctx);
    fmpz_mpoly_init(m, ctx);
    fmpz_mpoly_init(mpow, ctx);

    fmpz_mpolyv_init(new_lcs, ctx);
    fmpz_mpolyv_init(lc_divs, ctx);

    fmpz_poly_factor_init(Aufac);
    fmpz_poly_init(Au);

    alpha = _fmpz_vec_init(n);

    degs    = (slong *) flint_malloc((n + 1)*sizeof(slong));
    degeval = (slong *) flint_malloc((n + 1)*sizeof(slong));
    Aevals    = (fmpz_mpoly_struct *) flint_malloc(n*sizeof(fmpz_mpoly_struct));
    for (i = 0; i < n; i++)
        fmpz_mpoly_init(Aevals + i, ctx);

    fmpz_mpolyv_init(tfac, ctx);
    fmpz_mpoly_init(t, ctx);

    fmpz_mpoly_degrees_si(degs, A, ctx);

    alpha_count = 0;
    alpha_modulus = 1;
    goto got_alpha;

next_alpha:

    alpha_count++;
    if (alpha_count >= alpha_modulus)
    {
        alpha_count = 0;
        alpha_modulus++;
        if (alpha_modulus > 100)
        {
            success = 0;
            goto cleanup;
        }
    }

    for (i = 0; i < n; i++)
        fmpz_set_si(alpha + i, n_urandint(state, alpha_modulus) - alpha_modulus/2);

got_alpha:
/*
usleep(1000);
flint_printf("---------------------------\n");
flint_printf("(mod %wd) alpha = ", alpha_modulus); tuple_print(alpha, n);
*/
    /* ensure degrees do not drop under evaluation */
    for (i = n - 1; i >= 0; i--)
    {
        fmpz_mpoly_evaluate_one_fmpz(Aevals + i, i == n - 1 ? A : Aevals + i + 1, i + 1, alpha + i, ctx);
        fmpz_mpoly_degrees_si(degeval, Aevals + i, ctx);
        for (j = 0; j <= i; j++)
        {
            if (degeval[j] != degs[j])
                goto next_alpha;
        }
    }

    /* make sure our univar is squarefree */
    FLINT_ASSERT(fmpz_mpoly_is_fmpz_poly(Aevals + 0, 0, ctx));
    success = fmpz_mpoly_get_fmpz_poly(Au, Aevals + 0, 0, ctx);
    FLINT_ASSERT(success);
    fmpz_poly_factor(Aufac, Au);
    r = Aufac->num;
    for (j = 0; j < r; j++)
    {
        if (Aufac->exp[j] != 1)
            goto next_alpha;
    }

    /* if the univariate is irreducible, then A irreducible */
    if (r < 2)
    {
        fmpz_mpolyv_fit_length(fac, 1, ctx);
        fac->length = 1;
        fmpz_mpoly_set(fac->coeffs + 0, A, ctx);
        success = 1;
        goto cleanup;
    }

    if (lcAfac->num > 0)
    {
        if (!fmpz_mpoly_factor_lcc_wang(lc_divs, lcAfac, Aufac, alpha, ctx))
            goto next_alpha;
    }
    else
    {
        /* lcA is constant */
        fmpz_mpolyv_fit_length(lc_divs, r, ctx);
        lc_divs->length = r;
        for (i = 0; i < r; i++)
        {
            FLINT_ASSERT(Aufac->p[i].length > 0);
            fmpz_mpoly_set_fmpz(lc_divs->coeffs + i,
                             Aufac->p[i].coeffs + Aufac->p[i].length - 1, ctx);
        }
    }

    /*
        Assuming no extraneous divisors, we have

            A(X, x1, ..., xn) = F1 * ... * Fr

        and lead_divisor[i] is a divisor of lc_X(Fi). We also have the
        univariate factorization

            A(X, α1, ..., αn) = (c1 X^? + ... ) * ... * (cr X^? + ... )

        Set c(x1, ..., xn) = lc_X(A) and

            m(x1, ..., xn) = c/(prod_i lc_divs[i])   division is exact

        Lift the univariate factorization

            m(α)^(r-1) f(X, α) = (m(α)*lc_divs[0](α) X^? + ...) * ... * (m(α)*lc_divs[r-1](α) X^? + ...)

        against

            m(x1, ..., xn)^(r-1) f(X, x1, ..., xn)

        Note m(x1, ..., xn) is usually constant here, but it certainly does not have to be.
    */

    FLINT_ASSERT(r > 1);
    success = fmpz_mpoly_divides(m, lcA, lc_divs->coeffs + 0, ctx);
    FLINT_ASSERT(success);
    for (i = 1; i < r; i++)
    {
        success = fmpz_mpoly_divides(m, m, lc_divs->coeffs + i, ctx);
        FLINT_ASSERT(success);
    }
/*
printf("m: "); fmpz_mpoly_print_pretty(m, NULL, ctx); printf("\n");
*/
    fmpz_mpoly_pow_ui(mpow, m, r - 1, ctx);
    if (fmpz_mpoly_is_one(mpow, ctx))
    {
        newA = (fmpz_mpoly_struct *) A;
    }
    else
    {
        newA = Acopy;
        fmpz_mpoly_mul(newA, A, mpow, ctx);
    }

    fmpz_mpoly_degrees_si(degs, newA, ctx);

/*
printf("mpow: "); fmpz_mpoly_print_pretty(mpow, NULL, ctx); printf("\n");
flint_printf("modified A: "); fmpz_mpoly_print_pretty(newA, NULL, ctx); printf("\n");
*/
    fmpz_mpoly_set(t, mpow, ctx);
    for (i = n - 1; i >= 0; i--)
    {
        fmpz_mpoly_evaluate_one_fmpz(t, mpow, i + 1, alpha + i, ctx);
        fmpz_mpoly_swap(t, mpow, ctx);
        fmpz_mpoly_mul(Aevals + i, Aevals + i, mpow, ctx);
/*
flint_printf("Aeval[%wd]: ", i);  fmpz_mpoly_print_pretty(Aevals + i, NULL, ctx); printf("\n");
*/
    }

    fmpz_mpolyv_fit_length(new_lcs, (n + 1)*r, ctx);

    i = n;
    for (j = 0; j < r; j++)
    {
        fmpz_mpoly_mul(new_lcs->coeffs + i*r + j, lc_divs->coeffs + j, m, ctx);
    }
    for (i = n - 1; i >= 0; i--)
    {
        for (j = 0; j < r; j++)
        {
            fmpz_mpoly_evaluate_one_fmpz(new_lcs->coeffs + i*r + j,
                       new_lcs->coeffs + (i + 1)*r + j, i + 1, alpha + i, ctx);
        }
    }

    fmpz_mpolyv_fit_length(fac, r, ctx);
    fac->length = r;
    for (i = 0; i < r; i++)
    {
        FLINT_ASSERT(fmpz_mpoly_is_fmpz(new_lcs->coeffs + 0*r + i, ctx));
        FLINT_ASSERT(fmpz_mpoly_length(new_lcs->coeffs + 0*r + i, ctx) == 1);
        FLINT_ASSERT(fmpz_divisible(new_lcs->coeffs[i].coeffs + 0, Aufac->p[i].coeffs + Aufac->p[i].length - 1));
        fmpz_divexact(q, new_lcs->coeffs[i].coeffs + 0,
                                  Aufac->p[i].coeffs + Aufac->p[i].length - 1);
        _fmpz_mpoly_set_fmpz_poly(fac->coeffs + i, newA->bits,
                               Aufac->p[i].coeffs, Aufac->p[i].length, 0, ctx);
        fmpz_mpoly_scalar_mul_fmpz(fac->coeffs + i, fac->coeffs + i, q, ctx);
    }

    fmpz_mpolyv_fit_length(tfac, r, ctx);
    tfac->length = r;
    for (k = 1; k <= n; k++)
    {
        for (i = 0; i < r; i++)
        {
            _fmpz_mpoly_set_lead0(tfac->coeffs + i, fac->coeffs + i,
                                               new_lcs->coeffs + k*r + i, ctx);
        }

        success = fmpz_mpoly_hlift(k, tfac->coeffs, r, alpha,
                                         k < n ? Aevals + k : newA, degs, ctx);
        if (!success)
            goto next_alpha;

        fmpz_mpolyv_swap(tfac, fac, ctx);
    }

    if (fmpz_mpoly_is_fmpz(m, ctx))
    {
        for (i = 0; i < r; i++)
        {
            _fmpz_vec_content(q, fac->coeffs[i].coeffs, fac->coeffs[i].length);
            if (fmpz_sgn(fac->coeffs[i].coeffs + 0) < 0)
                fmpz_neg(q, q);
            fmpz_mpoly_scalar_divexact_fmpz(fac->coeffs + i,
                                            fac->coeffs + i, q, ctx);
        }
    }
    else
    {
        fmpz_mpoly_univar_t u;
        fmpz_mpoly_univar_init(u, ctx);
        for (i = 0; i < r; i++)
        {
            fmpz_mpoly_to_univar(u, fac->coeffs + i, 0, ctx);
            success = fmpz_mpoly_univar_content_mpoly(t, u, ctx);
            if (!success)
            {
                fmpz_mpoly_univar_clear(u, ctx);
                goto cleanup;
            }
            success = fmpz_mpoly_divides(fac->coeffs + i,
                                         fac->coeffs + i, t, ctx);
            FLINT_ASSERT(success);
            fmpz_mpoly_unit_normalize(fac->coeffs + i, ctx);
        }
        fmpz_mpoly_univar_clear(u, ctx);
    }

    success = 1;

cleanup:

    fmpz_clear(q);

    fmpz_mpolyv_clear(new_lcs, ctx);
    fmpz_mpolyv_clear(lc_divs, ctx);

    fmpz_poly_factor_clear(Aufac);
        
    _fmpz_vec_clear(alpha, n);

    for (i = 0; i < n; i++)
        fmpz_mpoly_clear(Aevals + i, ctx);
    flint_free(Aevals);

    flint_free(degs);
    flint_free(degeval);
    fmpz_mpolyv_clear(tfac, ctx);
    fmpz_mpoly_clear(t, ctx);

    fmpz_mpoly_clear(Acopy, ctx);
    fmpz_mpoly_clear(m, ctx);
    fmpz_mpoly_clear(mpow, ctx);

    fmpz_poly_clear(Au);

#if WANT_ASSERT
    if (success)
    {
        fmpz_mpoly_t prod;
        fmpz_mpoly_init(prod, ctx);
        fmpz_mpoly_one(prod, ctx);
        for (i = 0; i < fac->length; i++)
            fmpz_mpoly_mul(prod, prod, fac->coeffs + i, ctx);
        FLINT_ASSERT(fmpz_mpoly_equal(prod, A, ctx));
        fmpz_mpoly_clear(prod, ctx);
    }
#endif

    return success;
}