/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "nmod_mpoly_factor.h"


void subset_first(fmpz_t a, slong n, slong r)
{
    FLINT_ASSERT(r >= 0);
    FLINT_ASSERT(r <= n);
    fmpz_one(a);
    fmpz_mul_2exp(a, a, r); 
    fmpz_sub_ui(a, a, 1);
}

int subset_next(fmpz_t a, const fmpz_t b, slong n)
{
    slong t1, t2, i;
    int r;
    if (a == b)
    {
        fmpz_t t;
        fmpz_init(t);
        r = subset_next(t, b, n);
        fmpz_swap(a, t);
        fmpz_clear(t);
        return r;
    }

    i = 0;
    while (i < n && fmpz_tstbit(b, i) == 0)
        i++;
    t1 = i;
/*flint_printf("t1: %wd\n", t1);*/
    while (i<n && fmpz_tstbit(b,i) == 1)
        i++;
    t2 = i;
/*flint_printf("t2: %wd\n", t2);*/
    if (t2 < n)
    {
        fmpz_t t;
        fmpz_init_set_ui(t, 1);
        fmpz_one(a);
        fmpz_mul_2exp(a, a, n - t2);
        fmpz_sub_ui(a, a, 1);
        fmpz_mul_2exp(a, a, t2);
        fmpz_and(a, b, a);
        fmpz_setbit(a, t2);
        if (t2 > t1)
            fmpz_mul_2exp(t, t, t2 - t1 - 1);
        fmpz_sub_ui(t, t, 1);
        fmpz_add(a, a, t);
        fmpz_clear(t);
        return 1;
    }
    else
    {
        return 0;
    }
}

void subset_print(const fmpz_t a, slong n)
{
    slong i;
    for (i = n - 1; i >= 0; i --)
    {
        flint_printf("%d",fmpz_tstbit(a, i));
    }
}

void subset_map_down(fmpz_t a, const fmpz_t b, const fmpz_t m)
{
    ulong i, j, bbits = fmpz_bits(b);

    FLINT_ASSERT(a != b);
    FLINT_ASSERT(a != m);

    j = 0;
    fmpz_zero(a);
    for (i = 0; i < bbits; i++)
    {
        if (fmpz_tstbit(b, i))
        {
            FLINT_ASSERT(!fmpz_tstbit(m, i));
            fmpz_setbit(a, j++);
        }
        else
        {
            j += !fmpz_tstbit(m, i);
        }
    }
}
