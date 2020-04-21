/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "ui_factor.h"


static void _ui_factor_sort_terms(ui_factor_entry * a,
                        slong length, ulong posmask, ulong totalmask)
{
    slong i, j;
    slong cur, mid;

    if (length < 20)
    {
        for (i = 1; i < length; i++)
        {
            for (j = i; j > 0 && a[j-1].base > a[j].base; j--)
            {
                ulong t1 = a[j].base;
                ulong t2 = a[j].pow;
                a[j].base =  a[j-1].base;
                a[j].pow = a[j-1].pow;
                a[j-1].base = t1;
                a[j-1].pow = t2;                
            }
        }
        return;
    }

    if ((posmask & totalmask) == 0)
    {
        posmask >>= 1;
        if (posmask != 0)
        {
            _ui_factor_sort_terms(a, length, posmask, totalmask);
        }

        return;
    }

    mid = 0;
    while (mid < length && (a[mid].base & posmask) == 0)
    {
        mid++;
    }

    cur = mid;
    while (++cur < length)
    {
        if ((a[cur].base & posmask) == 0)
        {
            ulong t1 = a[cur].base;
            ulong t2 = a[cur].pow;
            a[cur].base =  a[mid].base;
            a[cur].pow = a[mid].pow;
            a[mid].base = t1;
            a[mid].pow = t2;
            mid++;
        }
    }

    posmask >>= 1;
    if (posmask != 0)
    {
        _ui_factor_sort_terms(a, mid, posmask, totalmask);
        _ui_factor_sort_terms(a + mid, length - mid, posmask, totalmask);
    }
}

void ui_factor_canonicalise(ui_factor_t f)
{
    ui_factor_entry * fd = f->data;
    slong fn = f->length;
    slong i, j;
    ulong smallest_idx, smallest, mask;

    if (unlikely(fn < 2))
        return;

    smallest_idx = 0;
    smallest = fd[smallest_idx].base;

    mask = smallest;
    j = 1;
    for (i = 1; i < fn; i++)
    {
        ulong t1 = fd[i].base;
        ulong t2 = fd[i].pow;
        mask |= t1;
        if (t1 == smallest)
        {
            fd[smallest_idx].pow += t2;
        }
        else
        {
            if (t1 < smallest)
            {
                smallest = t1;
                smallest_idx = j;
            }
            fd[j].base = t1;
            fd[j].pow = t2;
            j++;
        }
    }

    fn = j;

    _ui_factor_sort_terms(fd, fn, UWORD(1) << (FLINT_BIT_COUNT(mask)-1), mask);

    j = 0;
    for (i = 1; i < fn; i++)
    {
        FLINT_ASSERT(j < i);

        if (fd[j].base == fd[i].base)
        {
            fd[j].pow += fd[i].pow;
        }
        else
        {
            j++;
            fd[j].base = fd[i].base;
            fd[j].pow = fd[i].pow;
        }
    }

    j++;
    f->length = j;
}

