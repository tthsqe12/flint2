/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "ui_factor.h"


void ui_factor_push_ui_without_sieve(ui_factor_t f, ulong b)
{
    int i;
    slong fi = f->length;
    ui_factor_entry * fd;
    n_factor_t fac;

    n_factor_init(&fac);
    n_factor(&fac, b, 1);

    ui_factor_fit_length(f, 1 + FLINT_MAX_FACTORS_IN_LIMB + fi);

    fd = f->data;
    for (i = 0; i < fac.num; i++)
    {
        fd[fi].base = fac.p[i];
        fd[fi].pow = fac.exp[i];
        fi++;
    }

    f->length = fi;
}

/* f *= base */
void ui_factor_push_ui_with_sieve(ui_factor_t f, ulong b,
                                                     const ui_factor_sieve_t S)
{
    slong fi = f->length;
    ui_factor_entry * fd;
    ui_factor_sieve_entry * s;

    ui_factor_fit_length(f, fi + 16);

    fd = f->data;

    if ((b % 2) == 0)
    {
        ulong e = 0;
        do {
            e += 1;
            b /= 2;
        } while ((b % 2) == 0);

        fd[fi].base = 2;
        fd[fi].pow = e;
        fi++;
    }

    if (unlikely(b > S->max))
    {
        f->length = fi;
        ui_factor_push_ui_without_sieve(f, b);
        return;
    }

    s = S->array;

    b = b/2;
    if (b > 0)
    {
        fd[fi].base = s[b].pminus;
        b = s[b].cofactor;
        fd[fi].pow = 1;
        fi++;

        while (b > 0)
        {
            ulong p = s[b].pminus;
            b = s[b].cofactor;
            if (fd[fi - 1].base == p)
            {
                fd[fi - 1].pow++;
            }
            else
            {
                fd[fi].base = p;
                fd[fi].pow = 1;
                fi++;
            }
        }
    }

    f->length = fi;
}
