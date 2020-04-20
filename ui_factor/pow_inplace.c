/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "ui_factor.h"

/* f = f^pow */
void ui_factor_pow_inplace(ui_factor_t f, ulong pow)
{
    ui_factor_entry * fd = f->data;
    slong fi, fn = f->length;
    for (fi = 0; fi < fn; fi++)
        fd[fi].pow *= pow;
}
