/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "nmod_mpoly_factor.h"


slong n_bpoly_degree1(const n_bpoly_t A)
{
    slong i, len = 0;
    for (i = 0; i < A->length; i++)
        len = FLINT_MAX(len, A->coeffs[i].length);
    return len - 1;    
}

int n_bpoly_mod_is_canonical(const n_bpoly_t A, nmod_t mod)
{
    slong i;

    if (A->length <= 0)
        return A->length == 0;

    for (i = 0; i < A->length; i++)
    {
        if (!n_poly_mod_is_canonical(A->coeffs + i, mod))
        {
flint_printf("n_bpoly coeff i = %wd is bad\n", i);
            return 0;
        }
    }

    return !n_poly_is_zero(A->coeffs + A->length - 1);
}



int n_bpoly_equal(const n_bpoly_t A, const n_bpoly_t B)
{
    slong i;

    if (A->length != B->length)
        return 0;

    for (i = 0; i < A->length; i++)
    {
        if (!n_poly_equal(A->coeffs + i, B->coeffs + i))
            return 0;
    }

    return 1;
}


void n_bpoly_mod_mul(n_bpoly_t A, const n_bpoly_t B, const n_bpoly_t C,
                                                                    nmod_t mod)
{
    slong i, j;
    n_poly_struct * t;

    FLINT_ASSERT(A != B);
    FLINT_ASSERT(A != C);

    if (B->length <= 0 || C->length <= 0)
    {
        A->length = 0;
        return;
    }

    n_bpoly_fit_length(A, B->length + C->length);
    for (i = 0; i < B->length + C->length - 1; i++)
        n_poly_zero(A->coeffs + i);

    t = A->coeffs + B->length + C->length - 1;

    for (i = 0; i < B->length; i++)
    {
        for (j = 0; j < C->length; j++)
        {
            _n_poly_mod_mul(t, B->coeffs + i, C->coeffs + j, mod);
            n_poly_mod_add(A->coeffs + i + j, A->coeffs + i + j, t, mod);
        }
    }

    A->length = B->length + C->length - 1;
}
