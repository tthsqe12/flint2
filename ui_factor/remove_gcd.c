/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "ui_factor.h"


void ui_factor_remove_gcd(ui_factor_t f, ui_factor_t g)
{
    slong fi, gi, fj, gj;
    slong fn = f->length;
    slong gn = g->length;
    ui_factor_entry * fd = f->data;
    ui_factor_entry * gd = g->data;

    FLINT_ASSERT(ui_factor_is_canonical(f));
    FLINT_ASSERT(ui_factor_is_canonical(g));

    fi = gi = 0;
    fj = gj = 0;
    while (fj < fn && gj < gn)
    {
        if (fd[fj].base == gd[gj].base)
        {
            if (fd[fj].pow < gd[gj].pow)
            {
                FLINT_ASSERT(gi <= gj);
                gd[gi].base = gd[gj].base;
                gd[gi].pow = gd[gj].pow - fd[fj].pow;
                gi++;
            }
            else if (fd[fj].pow > gd[gj].pow)
            {
                FLINT_ASSERT(fi <= fj);
                fd[fi].base = fd[fj].base;
                fd[fi].pow = fd[fj].pow - gd[gj].pow;
                fi++;
            }
            fj++;
            gj++;
        }
        else if (fd[fj].base < gd[gj].base)
        {
            FLINT_ASSERT(fi <= fj);
            fd[fi].base = fd[fj].base;
            fd[fi].pow = fd[fj].pow;
            fi++;
            fj++;
        }
        else
        {
            FLINT_ASSERT(gi <= gj);
            gd[gi].base = gd[gj].base;
            gd[gi].pow = gd[gj].pow;
            gi++;
            gj++;
        }
    }

    while (fj < fn)
    {
        FLINT_ASSERT(fi <= fj);
        fd[fi].base = fd[fj].base;
        fd[fi].pow = fd[fj].pow;
        fi++;
        fj++;
    }

    while (gj < gn)
    {
        FLINT_ASSERT(gi <= gj);
        gd[gi].base = gd[gj].base;
        gd[gi].pow = gd[gj].pow;
        gi++;
        gj++;
    }

    f->length = fi;
    g->length = gi;

    FLINT_ASSERT(ui_factor_is_canonical(f));
    FLINT_ASSERT(ui_factor_is_canonical(g));
}
