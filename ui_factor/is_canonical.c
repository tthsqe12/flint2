/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "ui_factor.h"

/* bases strictly increasing, exponents > 0 */
int ui_factor_is_canonical(const ui_factor_t f)
{
    const ui_factor_entry * fd = f->data;
    ulong fi;
    for (fi = 0; fi < f->length; fi++)
    {
        if (fd[fi].pow == 0)
            return 0;

        if (fi > 0 && fd[fi - 1].base >= fd[fi].base)
            return 0;
    }
    return 1;
}
