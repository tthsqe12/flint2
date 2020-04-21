/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include <limits.h>
#include "ui_factor.h"


void ui_factor_sieve_init(ui_factor_sieve_t S)
{
    S->array = NULL;
    S->alloc = 0;
    S->max = 0;
}

void ui_factor_sieve_clear(ui_factor_sieve_t S)
{
    if (S->array)
        flint_free(S->array);
}

void ui_factor_sieve_fit_length(ui_factor_sieve_t S, slong length)
{
    if (S->alloc >= length)
        return;

    if (S->array)
    {
        S->array = (ui_factor_sieve_entry *) flint_realloc(S->array, length *
                                                sizeof(ui_factor_sieve_entry));
    }
    else
    {
        S->array = (ui_factor_sieve_entry *) flint_malloc(length *
                                                sizeof(ui_factor_sieve_entry));
    }

    S->alloc = length;
}

#define IDX(i) ((i)-1)/2

/* try to ensure that s can be used to factor numbers <= m */
void ui_factor_sieve_build(ui_factor_sieve_t S, ulong m)
{
    ulong i, k, n;
    ui_factor_sieve_entry * s;

    m = FLINT_MIN(m, INT_MAX);
    m |= 1;

    if (m <= S->max)
        return;

    ui_factor_sieve_fit_length(S, (m + 1)/2);
    s = S->array;

    s[IDX(1)].pminus = 0;
    s[IDX(1)].cofactor = 0;

    /*
        S->max = 0  =>  n = 1
        S->max = 1  =>  n = 1
        S->max = 2  =>  n = 1
        S->max = 3  =>  n = 3
        S->max = 4  =>  n = 3
    */
    n = FLINT_MAX(S->max, 1);
    n -= 1;
    n |= 1;

    for (i = n + 2; i <= m; i += 2)
    {
        s[IDX(i)].pminus = 0;
        s[IDX(i)].cofactor = 0;
    }

    for (i = 3; i <= n; i += 2)
    {
        FLINT_ASSERT(s[IDX(i)].pminus >= 3);
        if (s[IDX(i)].cofactor != 0)
            continue;

        /* i is prime. scan odd k*i in (n, m] */
        for (k = m/i, k -= !(k&1); k*i > n; k -= 2)
        {
            if (s[IDX(k*i)].pminus != 0)
                continue;
            s[IDX(k*i)].pminus = i;
            s[IDX(k*i)].cofactor = IDX(k);
        }
    }

    for (i = n + 2; i <= m; i += 2)
    {
        if (s[IDX(i)].pminus != 0)
            continue;

        /* i is prime. scan odd k*i >= i^2 with k*i in (n, m] */
        s[IDX(i)].pminus = i;
        s[IDX(i)].cofactor = IDX(1);
        for (k = i; k <= m/i; k += 2)
        {
            if (s[IDX(k*i)].pminus != 0)
                continue;
            s[IDX(k*i)].pminus = i;
            s[IDX(k*i)].cofactor = IDX(k);
        }
    }

    for (i = 3; i <= m; i += 2)
    {
        FLINT_ASSERT(s[IDX(i)].pminus >= 3);
        FLINT_ASSERT(i == s[IDX(i)].pminus * (2*s[IDX(i)].cofactor + 1));
    }

    S->max = m;
}

