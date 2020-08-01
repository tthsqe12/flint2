/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fq_nmod_mpoly_factor.h"

/*
    input A, B0, B1 with A(y,x) = B0(y,x) * B1(y,x) mod (y-alpha)
    return
       -1: B0(alpha,x) & B1(alpha,x) are not pairwise prime, or
           A(alpha,x) has wrong degree w.r.t x
        0: lift of B0 and B1 to true factors is impossible
        1: successfully lifted B0 and B1 to true factors without changing lc_x
*/
int fq_nmod_bpoly_hlift2(
    fq_nmod_bpoly_t A, /* clobbered (shifted by alpha) */
    fq_nmod_bpoly_t B0,
    fq_nmod_bpoly_t B1,
    const fq_nmod_t alpha,
    slong degree_inner, /* required degree in x */
    const fq_nmod_ctx_t ctx)
{
    int success;
    slong i, j;
    fq_nmod_poly_t c, s, t, u, v, g;
    fq_nmod_t malpha;

/*
flint_printf("fq_nmod_bpoly_hlift2 called: degree_inner = %wd\n", degree_inner);
flint_printf("alpha: "); fq_nmod_print_pretty(alpha, ctx); flint_printf("\n");
printf(" A: "); fq_nmod_bpoly_print_pretty(A, "Y", "X", ctx); flint_printf("\n");
printf("B0: "); fq_nmod_bpoly_print_pretty(B0, "Y", "X", ctx); flint_printf("\n");
printf("B1: "); fq_nmod_bpoly_print_pretty(B1, "Y", "X", ctx); flint_printf("\n");
*/
    FLINT_ASSERT(fq_nmod_bpoly_is_canonical(A, ctx));
    FLINT_ASSERT(fq_nmod_bpoly_is_canonical(B0, ctx));
    FLINT_ASSERT(fq_nmod_bpoly_is_canonical(B1, ctx));

    fq_nmod_poly_init(c, ctx);
    fq_nmod_poly_init(s, ctx);
    fq_nmod_poly_init(t, ctx);
    fq_nmod_poly_init(u, ctx);
    fq_nmod_poly_init(v, ctx);
    fq_nmod_poly_init(g, ctx);
    fq_nmod_init(malpha, ctx);

    fq_nmod_neg(malpha, alpha, ctx);

    fq_nmod_bpoly_taylor_shift_var0(A, alpha, ctx);
    fq_nmod_bpoly_taylor_shift_var0(B0, alpha, ctx);
    fq_nmod_bpoly_taylor_shift_var0(B1, alpha, ctx);

#if WANT_ASSERT
    {
        fq_nmod_poly_t T;
        fq_nmod_poly_init(T, ctx);
        fq_nmod_poly_mul(T, B0->coeffs + 0, B1->coeffs + 0, ctx);
        FLINT_ASSERT(fq_nmod_poly_equal(A->coeffs + 0, T, ctx));
        fq_nmod_poly_clear(T, ctx);
    }
#endif

    if (fq_nmod_poly_degree(A->coeffs + 0, ctx) != degree_inner)
    {
flint_printf("fq_nmod_bpoly_hlift2 fail 1\n");
        success = -1;
        goto cleanup;
    }

    /* the required degree in x is supposed to be deg_x(A) */
    FLINT_ASSERT(fq_nmod_bpoly_degree1(A, ctx) == fq_nmod_poly_degree(A->coeffs + 0, ctx));

    fq_nmod_poly_xgcd(g, s, t, B1->coeffs + 0, B0->coeffs + 0, ctx);
    if (!fq_nmod_poly_is_one(g, ctx))
    {
flint_printf("fq_nmod_bpoly_hlift2 fail 2\n");
        success = -2;
        goto cleanup;
    }

    fq_nmod_bpoly_fit_length(B0, A->length, ctx);
    fq_nmod_bpoly_fit_length(B1, A->length, ctx);
    for (j = 1; j < A->length; j++)
    {
        fq_nmod_poly_set(c, A->coeffs + j, ctx);
        for (i = 0; i <= j; i++)
        {
            if (i < B0->length && j - i < B1->length)
            {
                fq_nmod_poly_mul(t, B0->coeffs + i, B1->coeffs + j - i, ctx);
                fq_nmod_poly_sub(c, c, t, ctx);
            }
        }

        if (fq_nmod_poly_is_zero(c, ctx))
            continue;

        fq_nmod_poly_mul(t, s, c, ctx);
        fq_nmod_poly_divrem(g, u, t, B0->coeffs + 0, ctx);
        fq_nmod_poly_mul(t, u, B1->coeffs + 0, ctx);
        fq_nmod_poly_sub(c, c, t, ctx);
        fq_nmod_poly_divrem(v, g, c, B0->coeffs + 0, ctx);

        if (j < B0->length)
            fq_nmod_poly_add(B0->coeffs + j, B0->coeffs + j, u, ctx);
        else
            fq_nmod_poly_set(B0->coeffs + j, u, ctx);

        if (j < B1->length)
            fq_nmod_poly_add(B1->coeffs + j, B1->coeffs + j, v, ctx);
        else
            fq_nmod_poly_set(B1->coeffs + j, v, ctx);

        if (!fq_nmod_poly_is_zero(B0->coeffs + j, ctx))
            B0->length = FLINT_MAX(B0->length, j + 1);
        if (!fq_nmod_poly_is_zero(B1->coeffs + j, ctx))
            B1->length = FLINT_MAX(B1->length, j + 1);

        if (B0->length - 1 + B1->length - 1 > A->length - 1)
        {
flint_printf("fq_nmod_bpoly_hlift2 fail 3\n");
            success = 0;
            goto cleanup;
        }
    }

    fq_nmod_bpoly_taylor_shift_var0(B0, malpha, ctx);
    fq_nmod_bpoly_taylor_shift_var0(B1, malpha, ctx);

    success = 1;

cleanup:

/*
flint_printf("n_bpoly_hlift2 returning %d\n", success);
flint_printf("B0: "); n_bpoly_print_pretty(B0, "Y", "X"); flint_printf("\n");
flint_printf("B1: "); n_bpoly_print_pretty(B1, "Y", "X"); flint_printf("\n");
*/

    if (success > 0)
    {
        fq_nmod_bpoly_t tp1, tp2;
        fq_nmod_bpoly_init(tp1, ctx);
        fq_nmod_bpoly_init(tp2, ctx);

        fq_nmod_bpoly_taylor_shift_var0(A, malpha, ctx);
        fq_nmod_bpoly_mul(tp1, B0, B1, ctx);
        FLINT_ASSERT(fq_nmod_bpoly_equal(tp1, A, ctx));

        fq_nmod_bpoly_clear(tp1, ctx);
        fq_nmod_bpoly_clear(tp2, ctx);
    }

    fq_nmod_poly_clear(c, ctx);
    fq_nmod_poly_clear(s, ctx);
    fq_nmod_poly_clear(t, ctx);
    fq_nmod_poly_clear(u, ctx);
    fq_nmod_poly_clear(v, ctx);
    fq_nmod_poly_clear(g, ctx);
    fq_nmod_clear(malpha, ctx);

    return success;
}


/* r factor version */
int fq_nmod_bpoly_hlift(
    slong r,
    fq_nmod_bpoly_t A, /* clobbered (shifted by alpha) */
    fq_nmod_bpoly_struct * B,
    const fq_nmod_t alpha,
    slong degree_inner, /* required degree in x */
    const fq_nmod_ctx_t ctx)
{
    int success;
    slong i, j, k, tdeg;
    fq_nmod_poly_struct * s, * v;
    fq_nmod_poly_t c, t, u, g1, g2;
    fq_nmod_bpoly_struct * U;
    fq_nmod_t malpha;
/*
flint_printf("n_bpoly_hlift(%wd) called: degree_inner = %wd, alpha = %wu\n", r, degree_inner, alpha);
flint_printf(" A: ");
n_bpoly_print_pretty(A, "y", "x");
flint_printf("\n");
for (i = 0; i < r; i++)
{
flint_printf("B[%wd]: ", i);
n_bpoly_print_pretty(B + i, "y", "x");
flint_printf("\n");
}
*/
    FLINT_ASSERT(r > 2);
    FLINT_ASSERT(fq_nmod_bpoly_is_canonical(A, ctx));

    for (i = 0; i < r; i++)
    {
        FLINT_ASSERT(B[i].length > 0);
        FLINT_ASSERT(fq_nmod_bpoly_is_canonical(B + i, ctx));
    }

    U = (fq_nmod_bpoly_struct *) flint_malloc(r*sizeof(fq_nmod_bpoly_struct));
    for (i = 0; i < r; i++)
    {
        fq_nmod_bpoly_init(U + i, ctx);
        fq_nmod_bpoly_fit_length(U + i, A->length, ctx);
        for (j = 0; j < A->length; j++)
            fq_nmod_poly_zero(U[i].coeffs + j, ctx);
        U[i].length = A->length;
        fq_nmod_bpoly_fit_length(B + i, A->length, ctx);
    }

    s = (fq_nmod_poly_struct *) flint_malloc(r*sizeof(fq_nmod_poly_struct));
    v = (fq_nmod_poly_struct *) flint_malloc(r*sizeof(fq_nmod_poly_struct));
    for (i = 0; i < r; i++)
    {
        fq_nmod_poly_init(s + i, ctx);
        fq_nmod_poly_init(v + i, ctx);
    }

    fq_nmod_poly_init(c, ctx);
    fq_nmod_poly_init(t, ctx);
    fq_nmod_poly_init(u, ctx);
    fq_nmod_poly_init(g1, ctx);
    fq_nmod_poly_init(g2, ctx);
    fq_nmod_init(malpha, ctx);

    fq_nmod_neg(malpha, alpha, ctx);

    fq_nmod_bpoly_taylor_shift_var0(A, alpha, ctx);
    for (i = 0; i < r; i++)
        fq_nmod_bpoly_taylor_shift_var0(B + i, alpha, ctx);

    /* supposed to have A(alpha,x) = B0(alpha,x) * B1(alpha,x) * ... */
#if WANT_ASSERT
    {
        fq_nmod_poly_t T;
        fq_nmod_poly_init(T, ctx);
        fq_nmod_poly_mul(T, B[0].coeffs + 0, B[1].coeffs + 0, ctx);
        for (i = 2; i < r; i++)
            fq_nmod_poly_mul(T, T, B[i].coeffs + 0, ctx);
        FLINT_ASSERT(fq_nmod_poly_equal(A->coeffs + 0, T, ctx));
        fq_nmod_poly_clear(T, ctx);
    }
#endif

    if (fq_nmod_poly_degree(A->coeffs + 0, ctx) != degree_inner)
    {
flint_printf("fq_nmod_bpoly_hlift fail 1\n");
        success = -1;
        goto cleanup;
    }

    /* the required degree in x is supposed to be deg_x(A) */
    FLINT_ASSERT(fq_nmod_bpoly_degree1(A, ctx) == fq_nmod_poly_degree(A->coeffs + 0, ctx));

    for (i = 0; i < r; i++)
    {
        fq_nmod_poly_one(t, ctx);
        for (j = 0; j < r; j++)
        {
            if (j != i)
                fq_nmod_poly_mul(t, t, B[j].coeffs + 0, ctx);
        }

        fq_nmod_poly_xgcd(g1, s + i, g2, t, B[i].coeffs + 0, ctx);
        if (!fq_nmod_poly_is_one(g1, ctx))
        {
flint_printf("fq_nmod_bpoly_hlift fail 2\n");
            success = -1;
            goto cleanup;
        }
    }

    k = r - 2;
    fq_nmod_poly_mul(U[k].coeffs + 0, B[k].coeffs + 0, B[k + 1].coeffs + 0, ctx);
    for (k--; k > 0; k--)
        fq_nmod_poly_mul(U[k].coeffs + 0, B[k].coeffs + 0, U[k + 1].coeffs + 0, ctx);

    for (j = 1; j < A->length; j++)
    {
        for (k = 0; k < r; k++)
            fq_nmod_poly_zero(U[k].coeffs + j, ctx);

        k = r - 2;
        fq_nmod_poly_zero(U[k].coeffs + j, ctx);
        for (i = 0; i <= j; i++)
        {
            if (i < B[k].length && j - i < B[k + 1].length)
            {
                fq_nmod_poly_mul(t, B[k].coeffs + i, B[k + 1].coeffs + j - i, ctx);
                fq_nmod_poly_add(U[k].coeffs + j, U[k].coeffs + j, t, ctx);
            }
        }
        for (k--; k > 0; k--)
        {
            fq_nmod_poly_zero(U[k].coeffs + j, ctx);
            for (i = 0; i <= j; i++)
            {
                if (i < B[k].length)
                {
                    fq_nmod_poly_mul(t, B[k].coeffs + i, U[k + 1].coeffs + j - i, ctx);
                    fq_nmod_poly_add(U[k].coeffs + j, U[k].coeffs + j, t, ctx);
                }
            }
        }

        fq_nmod_poly_set(c, A->coeffs + j, ctx);

        for (i = 0; i <= j; i++)
        {
            if (i < B[0].length)
            {
                fq_nmod_poly_mul(t, B[0].coeffs + i, U[1].coeffs + j - i, ctx);
                fq_nmod_poly_sub(c, c, t, ctx);
            }
        }

        if (fq_nmod_poly_is_zero(c, ctx))
            continue;

        tdeg = 0;
        for (i = 0; i < r; i++)
        {
            fq_nmod_poly_mul(t, s + i, c, ctx);
            fq_nmod_poly_divrem(g1, v + i, t, B[i].coeffs + 0, ctx);
            while (j >= B[i].length)
            {
                fq_nmod_poly_zero(B[i].coeffs + B[i].length, ctx);
                B[i].length++;
            }
            fq_nmod_poly_add(B[i].coeffs + j, B[i].coeffs + j, v + i, ctx);
            fq_nmod_bpoly_normalise(B + i, ctx);
            tdeg += B[i].length - 1;
        }

        if (tdeg >= A->length)
        {
flint_printf("fq_nmod_bpoly_hlift fail 3\n");
            success = 0;
            goto cleanup;
        }

        k = r - 2;
        fq_nmod_poly_mul(t, B[k].coeffs + 0, v + k + 1, ctx);
        fq_nmod_poly_mul(u, B[k + 1].coeffs + 0, v + k, ctx);
        fq_nmod_poly_add(t, t, u, ctx);
        fq_nmod_poly_add(U[k].coeffs + j, U[k].coeffs + j, t, ctx);
        for (k--; k > 0; k--)
        {
            fq_nmod_poly_mul(u, B[k].coeffs + 0, t, ctx);
            fq_nmod_poly_mul(t, U[k + 1].coeffs + 0, v + k, ctx);
            fq_nmod_poly_add(t, t, u, ctx);
            fq_nmod_poly_add(U[k].coeffs + j, U[k].coeffs + j, t, ctx);
        }
    }

    for (i = 0; i < r; i++)
        fq_nmod_bpoly_taylor_shift_var0(B + i, malpha, ctx);

    success = 1;

cleanup:
/*
flint_printf("n_bpoly_hlift(%wd) returning %d\n", r, success);
for (i = 0; i < r; i++)
{
flint_printf("B[%wd]: ", i);
n_bpoly_print_pretty(B + i, "y", "x");
flint_printf("\n");
}

FLINT_ASSERT(success);
*/
    if (success > 0)
    {
        fq_nmod_bpoly_t tp1, tp2;
        fq_nmod_bpoly_init(tp1, ctx);
        fq_nmod_bpoly_init(tp2, ctx);

        fq_nmod_bpoly_taylor_shift_var0(A, malpha, ctx);
        fq_nmod_bpoly_mul(tp1, B + 0, B + 1, ctx);
        for (i = 2; i < r; i++)
        {
            fq_nmod_bpoly_mul(tp2, tp1, B + i, ctx);
            fq_nmod_bpoly_swap(tp1, tp2, ctx);
        }
        FLINT_ASSERT(fq_nmod_bpoly_equal(tp1, A, ctx));

        fq_nmod_bpoly_clear(tp1, ctx);
        fq_nmod_bpoly_clear(tp2, ctx);
    }

    for (i = 0; i < r; i++)
    {
        fq_nmod_bpoly_clear(U + i, ctx);
        fq_nmod_poly_clear(s + i, ctx);
        fq_nmod_poly_clear(v + i, ctx);
    }
    flint_free(U);
    flint_free(s);
    flint_free(v);

    fq_nmod_poly_clear(c, ctx);
    fq_nmod_poly_clear(t, ctx);
    fq_nmod_poly_clear(u, ctx);
    fq_nmod_poly_clear(g1, ctx);
    fq_nmod_poly_clear(g2, ctx);

    fq_nmod_clear(malpha, ctx);

    return success;
}

