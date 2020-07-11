/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fmpz_mpoly_factor.h"


void fmpz_bpoly_clear(fmpz_bpoly_t A)
{
    if (A->alloc > 0)
    {
        slong i;
        for (i = 0; i < A->alloc; i++)
            fmpz_poly_clear(A->coeffs + i);
        flint_free(A->coeffs);
    }
}

void fmpz_bpoly_print_pretty(fmpz_bpoly_t A, const char * var0, const char * var1)
{
    slong i;
    int first = 1;

    for (i = A->length - 1; i >= 0; i--)
    {
        if (fmpz_poly_is_zero(A->coeffs + i))
            continue;

        if (!first)
            flint_printf(" + ");
        first = 0;

        flint_printf("(");
        fmpz_poly_print_pretty(A->coeffs + i, var1);
        flint_printf(")*%s^%wd", var0, i);
    }

    if (first)
        flint_printf("0");
}

void fmpz_bpoly_realloc(fmpz_bpoly_t A, slong len)
{
    slong i;
    slong old_alloc = A->alloc;

    if (len <= old_alloc)
        return;

    len = FLINT_MAX(len, 2*old_alloc);

    if (A->alloc == 0)
        A->coeffs = (fmpz_poly_struct *) flint_malloc(
                                               len * sizeof(fmpz_poly_struct));
    else
        A->coeffs = (fmpz_poly_struct *) flint_realloc(A->coeffs,
                                               len * sizeof(fmpz_poly_struct));

    for (i = old_alloc; i < len; i++)
        fmpz_poly_init(A->coeffs + i);

    A->alloc = len;
}

void fmpz_bpoly_set_coeff(fmpz_bpoly_t A, slong xi, slong yi, const fmpz_t c)
{
    slong i;

    FLINT_ASSERT(!fmpz_is_zero(c));

    if (xi >= A->length)
    {
        fmpz_bpoly_fit_length(A, xi + 1);
        for (i = A->length; i <= xi; i++)
            fmpz_poly_zero(A->coeffs + i);
        A->length = xi + 1;
    }

    fmpz_poly_set_coeff_fmpz(A->coeffs + xi, yi, c);
}
