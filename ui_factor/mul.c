/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "ui_factor.h"


void ui_factor_mul(ui_factor_t z, const ui_factor_t f, const ui_factor_t g)
{
    slong fi, gi, zi;
    slong fn = f->length;
    slong gn = g->length;
    const ui_factor_entry * fd, * gd;
    ui_factor_entry * zd;

    FLINT_ASSERT(z != f);
    FLINT_ASSERT(z != f);

    ui_factor_fit_length(z, fn + gn);

    zd = z->data;
    fd = f->data;
    gd = g->data;

    fi = gi = zi = 0;

    while (fi < fn && gi < gn)
    {
        if (fd[fi].base == gd[gi].base)
        {
            zd[zi].base = fd[fi].base;
            zd[zi].pow = fd[fi].pow + gd[gi].pow;
            fi++;
            gi++;
        }
        else if (fd[fi].base < gd[gi].base)
        {
            zd[zi].base = fd[fi].base;
            zd[zi].pow = fd[fi].pow;
            fi++;
        }
        else
        {
            zd[zi].base = gd[gi].base;
            zd[zi].pow = gd[gi].pow;
            gi++;
        }

        zi++;
    }

    while (fi < fn)
    {
        zd[zi].base = fd[fi].base;
        zd[zi].pow = fd[fi].pow;
        fi++;
        zi++;
    }

    while (gi < gn)
    {
        zd[zi].base = gd[gi].base;
        zd[zi].pow = gd[gi].pow;
        gi++;
        zi++;
    }

    FLINT_ASSERT(zi <= z->alloc);
    z->length = zi;

    FLINT_ASSERT(!(ui_factor_is_canonical(f) && ui_factor_is_canonical(g)) ||
                 ui_factor_is_canonical(z));
}

