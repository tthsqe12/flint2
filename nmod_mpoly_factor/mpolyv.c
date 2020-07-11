/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "nmod_mpoly_factor.h"


void nmod_mpolyv_clear(nmod_mpolyv_t A, const nmod_mpoly_ctx_t ctx)
{
    slong i;
    for (i = 0; i < A->alloc; i++)
        nmod_mpoly_clear(A->coeffs + i, ctx);
    flint_free(A->coeffs);
}


void nmod_mpolyv_print_pretty(
    const nmod_mpolyv_t poly,
    const char ** x,
    const nmod_mpoly_ctx_t ctx)
{
    slong i;
    for (i = 0; i < poly->length; i++)
    {
        flint_printf("coeff[%wd]: ", i);
        nmod_mpoly_print_pretty(poly->coeffs + i, x, ctx);
        flint_printf("\n");
    }
}


void nmod_mpolyv_fit_length(
    nmod_mpolyv_t A,
    slong length,
    const nmod_mpoly_ctx_t ctx)
{
    slong i;
    slong old_alloc = A->alloc;
    slong new_alloc = FLINT_MAX(length, 2*A->alloc);

    if (length <= old_alloc)
        return;

    if (old_alloc > 0)
    {
        A->coeffs = (nmod_mpoly_struct *) flint_realloc(A->coeffs,
                                          new_alloc*sizeof(nmod_mpoly_struct));
    }
    else
    {
        A->coeffs = (nmod_mpoly_struct *) flint_malloc(
                                          new_alloc*sizeof(nmod_mpoly_struct));
    }

    for (i = old_alloc; i < new_alloc; i++)
        nmod_mpoly_init(A->coeffs + i, ctx);

    A->alloc = new_alloc;
}
