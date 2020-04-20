/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "ui_factor.h"


void ui_factor_fit_length(ui_factor_t f, ulong length)
{
    if (f->alloc >= length)
        return;

    length = FLINT_MAX(length, f->alloc + f->alloc/2);

    if (f->data)
    {
        f->data = (ui_factor_entry *) flint_realloc(f->data, length *
                                                      sizeof(ui_factor_entry));
    }
    else
    {
        f->data = (ui_factor_entry *) flint_malloc(length *
                                                      sizeof(ui_factor_entry));
    }

    f->alloc = length;
}

