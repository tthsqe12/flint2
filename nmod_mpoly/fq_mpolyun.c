/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "nmod_mpoly.h"

void fq_nmod_mpolyun_init(fq_nmod_mpolyun_t A, mp_bitcnt_t bits, const fq_nmod_mpoly_ctx_t ctx)
{
    A->coeffs = NULL;
    A->exps = NULL;
    A->alloc = 0;
    A->length = 0;
    A->bits = bits;
}

void fq_nmod_mpolyun_clear(fq_nmod_mpolyun_t A, const fq_nmod_mpoly_ctx_t ctx)
{
    slong i;
    for (i = 0; i < A->alloc; i++)
        fq_nmod_mpolyn_clear(A->coeffs + i, ctx);
    flint_free(A->coeffs);
    flint_free(A->exps);
}


void fq_nmod_mpolyun_swap(fq_nmod_mpolyun_t A, fq_nmod_mpolyun_t B)
{
   fq_nmod_mpolyun_struct t = *A;
   *A = *B;
   *B = t;
}

void fq_nmod_mpolyun_zero(fq_nmod_mpolyun_t A, const fq_nmod_mpoly_ctx_t ctx)
{
    slong i;
    for (i = 0; i < A->alloc; i++) {
        fq_nmod_mpolyn_clear(A->coeffs + i, ctx);
        fq_nmod_mpolyn_init(A->coeffs + i, A->bits, ctx);
    }
    A->length = 0;
}

void fq_nmod_mpolyun_print_pretty(const fq_nmod_mpolyun_t poly,
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
        FLINT_ASSERT((poly->coeffs + i)->bits == poly->bits);
        fq_nmod_mpolyn_print_pretty(poly->coeffs + i,x,ctx);
        flint_printf(")*X^%wd",poly->exps[i]);
    }
}

void fq_nmod_mpolyun_fit_length(fq_nmod_mpolyun_t A, slong length, const fq_nmod_mpoly_ctx_t ctx)
{
    slong i;
    slong old_alloc = A->alloc;
    slong new_alloc = FLINT_MAX(length, 2*A->alloc);

    if (length > old_alloc)
    {
        if (old_alloc == 0)
        {
            A->exps = (ulong *) flint_malloc(new_alloc*sizeof(ulong));
            A->coeffs = (fq_nmod_mpolyn_struct *) flint_malloc(
                                      new_alloc*sizeof(fq_nmod_mpolyn_struct));
        } else
        {
            A->exps = (ulong *) flint_realloc(A->exps,
                                                      new_alloc*sizeof(ulong));
            A->coeffs = (fq_nmod_mpolyn_struct *) flint_realloc(A->coeffs,
                                      new_alloc*sizeof(fq_nmod_mpolyn_struct));
        }

        for (i = old_alloc; i < new_alloc; i++)
        {
            fq_nmod_mpolyn_init(A->coeffs + i, A->bits, ctx);
        }
        A->alloc = new_alloc;
    }
}
