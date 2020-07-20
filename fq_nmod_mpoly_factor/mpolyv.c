/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fq_nmod_mpoly_factor.h"


void fq_nmod_mpolyv_clear(fq_nmod_mpolyv_t A, const fq_nmod_mpoly_ctx_t ctx)
{
    slong i;
    if (A->alloc > 0)
    {
        for (i = 0; i < A->alloc; i++)
            fq_nmod_mpoly_clear(A->coeffs + i, ctx);
        flint_free(A->coeffs);
    }
}


void fq_nmod_mpolyv_print_pretty(
    const fq_nmod_mpolyv_t poly,
    const char ** x,
    const fq_nmod_mpoly_ctx_t ctx)
{
    slong i;
    for (i = 0; i < poly->length; i++)
    {
        flint_printf("coeff[%wd]: ", i);
        fq_nmod_mpoly_print_pretty(poly->coeffs + i, x, ctx);
        flint_printf("\n");
    }
}


void fq_nmod_mpolyv_fit_length(
    fq_nmod_mpolyv_t A,
    slong length,
    const fq_nmod_mpoly_ctx_t ctx)
{
    slong i;
    slong old_alloc = A->alloc;
    slong new_alloc = FLINT_MAX(length, 2*A->alloc);

    if (length <= old_alloc)
        return;

    if (old_alloc > 0)
    {
        A->coeffs = (fq_nmod_mpoly_struct *) flint_realloc(A->coeffs,
                                       new_alloc*sizeof(fq_nmod_mpoly_struct));
    }
    else
    {
        A->coeffs = (fq_nmod_mpoly_struct *) flint_malloc(
                                       new_alloc*sizeof(fq_nmod_mpoly_struct));
    }

    for (i = old_alloc; i < new_alloc; i++)
        fq_nmod_mpoly_init(A->coeffs + i, ctx);

    A->alloc = new_alloc;
}
