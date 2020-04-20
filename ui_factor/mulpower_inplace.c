/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "ui_factor.h"

/* f *= g^p, like addmul on the exponents */
void ui_factor_mulpow_inplace(ui_factor_t f, const ui_factor_t g, ulong p)
{
    ui_factor_entry * fd = f->data;
    const ui_factor_entry * gd = g->data;
    slong fn = f->length;
    slong gn = g->length;
    slong gi;
    slong i;
    slong fi = 0;
    for (gi = 0; gi < gn; gi++)
    {
        while (1)
        {
            if ((fi >= fn))
            {
                ui_factor_fit_length(f, fi + gn - gi);
                fd = f->data;
                for ( ; gi < gn; gi++, fi++)
                {
                    fd[fi].base = gd[gi].base;
                    fd[fi].pow = p*gd[gi].pow;
                }
                f->length = fi;
                return;
            }

            if (fd[fi].base >= gd[gi].base)
                break;

            fi++;
        }

        if ((fd[fi].base > gd[gi].base))
        {
            ui_factor_entry t1, t2;
            ui_factor_fit_length(f, fn + 1);
            fd = f->data;
            t2.base = gd[gi].base;
            t2.pow = p*gd[gi].pow;
            for (i = fi; i <= fn; i++)
            {
                t1 = fd[fi];
                fd[fi] = t2;
                t2 = t1;
            }
            fn++;
            fi++;
            f->length = fn;
        }
        else
        {
            fd[fi].pow += p*gd[gi].pow;
            fi++;
        }
    }
}

