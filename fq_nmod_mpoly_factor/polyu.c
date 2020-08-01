/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fq_nmod_mpoly_factor.h"


void fq_nmod_polyu_clear(
    fq_nmod_polyu_t A,
    const fq_nmod_ctx_t ctx)
{
    if (A->alloc > 0)
    {
        slong i;
        for (i = 0; i < A->alloc; i++)
            fq_nmod_clear(A->coeffs + i, ctx);
        flint_free(A->exps);
        flint_free(A->coeffs);
    }
    else
    {
        FLINT_ASSERT(A->coeffs == NULL);
        FLINT_ASSERT(A->exps == NULL);
    }
}

void fq_nmod_polyu_realloc(fq_nmod_polyu_t A, slong len, const fq_nmod_ctx_t ctx)
{
    slong i;
    slong old_alloc = A->alloc;
    slong new_alloc = FLINT_MAX(len, old_alloc + 1 + old_alloc/2);

    FLINT_ASSERT(A->alloc >= 0);
    if (len <= A->alloc)
        return;

    if (old_alloc > 0)
    {
        A->coeffs = (fq_nmod_struct *) flint_realloc(A->coeffs,
                                             new_alloc*sizeof(fq_nmod_struct));
        A->exps = (ulong *) flint_realloc(A->exps, new_alloc*sizeof(ulong));
    }
    else
    {
        FLINT_ASSERT(A->coeffs == NULL);
        FLINT_ASSERT(A->exps == NULL);
        A->coeffs = (fq_nmod_struct *) flint_malloc(
                                             new_alloc*sizeof(fq_nmod_struct));
        A->exps = (ulong *) flint_malloc(new_alloc*sizeof(ulong));
    }

    for (i = old_alloc; i < new_alloc; i++)
        fq_nmod_init(A->coeffs + i, ctx);

    A->alloc = new_alloc;
}

void fq_nmod_polyu3_print_pretty(
    const fq_nmod_polyu_t A,
    const char * var0,
    const char * var1,
    const char * var2,
    const fq_nmod_ctx_t ctx)
{
    slong i;
    int first = 1;

    for (i = 0; i < A->length; i++)
    {
        if (!first)
            printf(" + ");
        first = 0;
        flint_printf("(");
        fq_nmod_print_pretty(A->coeffs + i, ctx);
        flint_printf(")*%s^%wu*%s^%wu*%s^%wu",
            var0, extract_exp(A->exps[i], 2, 3),
            var1, extract_exp(A->exps[i], 1, 3),
            var2, extract_exp(A->exps[i], 0, 3));
    }

    if (first)
        flint_printf("0");
}

void fq_nmod_polyu3_degrees(
    slong * deg0,
    slong * deg1,
    slong * deg2,
    const fq_nmod_polyu_t A)
{
    slong i;
    ulong m;
    ulong mask = mpoly_overflow_mask_sp(FLINT_BITS/3);

    if (A->length < 1)
    {
        *deg0 = *deg1 = *deg2 = -1;
        return;
    }

    m = A->exps[0];
    for (i = 1; i < A->length; i++)
        m = mpoly_monomial_max1(m, A->exps[i], FLINT_BITS/3, mask);

    *deg0 = extract_exp(m, 2, 3);
    *deg1 = extract_exp(m, 1, 3);
    *deg2 = extract_exp(m, 0, 3);
}

int fq_nmod_polyu_is_canonical(
    const fq_nmod_polyu_t A,
    const fq_nmod_ctx_t ctx)
{
    slong i;
    if (A->length < 0)
        return 0;
    for (i = 0; i < A->length; i++)
    {
        if (fq_nmod_is_zero(A->coeffs + i, ctx))
            return 0;
        if (i > 0 && A->exps[i] >= A->exps[i - 1])
            return 0;
    }
    return 1;
}
