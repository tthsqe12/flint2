/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "ui_factor.h"


void ui_factor_stack_init(ui_factor_stack_t S)
{
    S->mpz_alloc = 0;
    S->mpz_array = NULL;
    S->mpz_top = 0;

    S->factor_alloc = 0;
    S->factor_array = NULL;
    S->factor_top = 0;
}

void ui_factor_stack_clear(ui_factor_stack_t S)
{
    slong i;

    FLINT_ASSERT(S->mpz_top == 0);

    for (i = 0; i < S->mpz_alloc; i++)
    {
        mpz_clear(S->mpz_array[i]);
        flint_free(S->mpz_array[i]);
    }

    if (S->mpz_array)
        flint_free(S->mpz_array);

    S->mpz_array = NULL;
    S->mpz_alloc = 0;

    FLINT_ASSERT(S->factor_top == 0);

    for (i = 0; i < S->factor_alloc; i++)
    {
        ui_factor_clear(S->factor_array[i]);
        flint_free(S->factor_array[i]);
    }

    if (S->factor_array)
        flint_free(S->factor_array);

    S->factor_array = NULL;
    S->factor_alloc = 0;
}

/* insure that k slots are available after top and return pointer to top */
__mpz_struct ** ui_factor_stack_fit_request_mpz(ui_factor_stack_t S, slong k)
{
    slong newalloc, i;

    FLINT_ASSERT(S->mpz_alloc >= S->mpz_top);

    if (S->mpz_top + k > S->mpz_alloc)
    {
        newalloc = FLINT_MAX(WORD(1), S->mpz_top + k);

        if (S->mpz_array)
        {
            S->mpz_array = (__mpz_struct **) flint_realloc(S->mpz_array,
                                           newalloc*sizeof(__mpz_struct*));
        }
        else
        {
            S->mpz_array = (__mpz_struct **) flint_malloc(
                                           newalloc*sizeof(__mpz_struct*));
        }

        for (i = S->mpz_alloc; i < newalloc; i++)
        {
            S->mpz_array[i] = (__mpz_struct *) flint_malloc(
                                                     sizeof(__mpz_struct));
            mpz_init(S->mpz_array[i]);
        }
        S->mpz_alloc = newalloc;
    }

    return S->mpz_array + S->mpz_top;
}

ui_factor_struct ** ui_factor_stack_fit_request_factor(ui_factor_stack_t S, slong k)
{
    slong newalloc, i;

    FLINT_ASSERT(S->factor_alloc >= S->factor_top);

    if (S->factor_top + k > S->factor_alloc)
    {
        newalloc = FLINT_MAX(WORD(1), S->factor_top + k);

        if (S->factor_array)
        {
            S->factor_array = (ui_factor_struct **) flint_realloc(S->factor_array,
                                           newalloc*sizeof(ui_factor_struct*));
        }
        else
        {
            S->factor_array = (ui_factor_struct **) flint_malloc(
                                           newalloc*sizeof(ui_factor_struct*));
        }

        for (i = S->factor_alloc; i < newalloc; i++)
        {
            S->factor_array[i] = (ui_factor_struct *) flint_malloc(
                                                     sizeof(ui_factor_struct));
            ui_factor_init(S->factor_array[i]);
        }
        S->factor_alloc = newalloc;
    }

    return S->factor_array + S->factor_top;
}
