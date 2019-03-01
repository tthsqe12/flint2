/*
    Copyright (C) 2019 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "nmod_poly.h"
#include "nmod_mpoly.h"

void nmod_poly_stack_init(nmod_poly_stack_t S, slong size, mp_limb_t n)
{
    slong alloc, i;

    /* lets NOT deal with null pointers and zero alloc's */
    alloc = FLINT_MAX(size, WORD(1));
    S->array = (nmod_poly_struct **) flint_malloc(alloc*sizeof(nmod_poly_struct *));
    for (i = 0; i < alloc; i++)
    {
        S->array[i] = (nmod_poly_struct *) flint_malloc(sizeof(nmod_poly_struct));
        nmod_poly_init(S->array[i], n);
    }
    S->alloc = alloc;
    S->top = 0;
    S->modulus = n;
}

void nmod_poly_stack_clear(nmod_poly_stack_t S)
{
    slong i;

    FLINT_ASSERT(S->top == 0);

    for (i = 0; i < S->alloc; i++)
    {
        nmod_poly_clear(S->array[i]);
        flint_free(S->array[i]);
    }
    flint_free(S->array);
    S->array = NULL;
}

/* insure that k slots are available after top and return pointer to top */
nmod_poly_struct ** nmod_poly_stack_fit_request(nmod_poly_stack_t S, slong k)
{
    slong newalloc, i;

    FLINT_ASSERT(S->alloc >= S->top);

    if (S->top + k > S->alloc)
    {
        newalloc = S->top + k;
        S->array = (nmod_poly_struct **) flint_realloc(S->array,
                                           newalloc*sizeof(nmod_poly_struct*));
        for (i = S->alloc; i < newalloc; i++)
        {
            S->array[i] = (nmod_poly_struct *) flint_malloc(
                                                     sizeof(nmod_poly_struct));
            nmod_poly_init(S->array[i], S->modulus);
        }
        S->alloc = newalloc;
    }

    return S->array + S->top;
}
