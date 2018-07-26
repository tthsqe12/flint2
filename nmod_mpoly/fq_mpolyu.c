/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "nmod_mpoly.h"


void fq_nmod_mpolyu_init(fq_nmod_mpolyu_t A, mp_bitcnt_t bits, const fq_nmod_mpoly_ctx_t ctx)
{
    A->coeffs = NULL;
    A->exps = NULL;
    A->alloc = 0;
    A->length = 0;
    A->bits = bits;
}


void fq_nmod_mpolyu_clear(fq_nmod_mpolyu_t A, const fq_nmod_mpoly_ctx_t uctx)
{
    slong i;
    for (i = 0; i < A->alloc; i++)
        fq_nmod_mpoly_clear(A->coeffs + i, uctx);
    flint_free(A->coeffs);
    flint_free(A->exps);
}

void fq_nmod_mpolyu_print_pretty(const fq_nmod_mpolyu_t poly,
                                const char ** x, const fq_nmod_mpoly_ctx_t ctx)
{
    slong i;
    if (poly->length == 0)
        flint_printf("0");
    for (i = 0; i < poly->length; i++)
    {
        if (i != 0)
            flint_printf(" + ");
        flint_printf("(");
        fq_nmod_mpoly_print_pretty(poly->coeffs + i,x,ctx);
        flint_printf(")*X^%wd",poly->exps[i]);
    }
}

void fq_nmod_mpolyu_fit_length(fq_nmod_mpolyu_t A, slong length,
                                                const fq_nmod_mpoly_ctx_t uctx)
{
    slong i;
    slong old_alloc = A->alloc;
    slong new_alloc = FLINT_MAX(length, 2*A->alloc);

    if (length > old_alloc)
    {
        if (old_alloc == 0)
        {
            A->exps = (ulong *) flint_malloc(new_alloc*sizeof(ulong));
            A->coeffs = (fq_nmod_mpoly_struct *) flint_malloc(
                                       new_alloc*sizeof(fq_nmod_mpoly_struct));
        } else
        {
            A->exps = (ulong *) flint_realloc(A->exps,
                                                      new_alloc*sizeof(ulong));
            A->coeffs = (fq_nmod_mpoly_struct *) flint_realloc(A->coeffs,
                                       new_alloc*sizeof(fq_nmod_mpoly_struct));
        }

        for (i = old_alloc; i < new_alloc; i++)
        {
            fq_nmod_mpoly_init(A->coeffs + i, A->bits, uctx);
        }
        A->alloc = new_alloc;
    }
}
