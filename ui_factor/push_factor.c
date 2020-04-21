/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "ui_factor.h"


void ui_factor_push_factor(ui_factor_t f, ulong base, ulong pow)
{
    ui_factor_fit_length(f, f->length + 1);
    f->data[f->length].base = base;
    f->data[f->length].pow = pow;
    f->length++;
}

