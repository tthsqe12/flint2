/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "ui_factor.h"


mp_size_t flint_mpn_mul_111(mp_limb_t * y, const mp_limb_t * x,
                         mp_size_t n, mp_limb_t a1, mp_limb_t a2, mp_limb_t a3)
{
    mp_limb_t hi, lo, p1, p2;

    umul_ppmm(hi, lo, a1, a2);
    if (hi == 0)
    {
        p1 = lo;
        umul_ppmm(hi, lo, p1, a3);
        if (hi == 0)
        {
            y[n] = mpn_mul_1(y, x, n, lo); n += (y[n] != 0);
            return n;
        }
        p2 = a3;
        goto do_two;
    }

    umul_ppmm(hi, lo, a1, a3);
    if (hi == 0)
    {
        p1 = lo;
        p2 = a2;
        goto do_two;
    }

    umul_ppmm(hi, lo, a2, a3);
    if (hi == 0)
    {
        p1 = lo;
        p2 = a1;
        goto do_two;
    }    

    y[n] = mpn_mul_1(y, x, n, a1); n += (y[n] != 0);
    y[n] = mpn_mul_1(y, y, n, a2); n += (y[n] != 0);
    y[n] = mpn_mul_1(y, y, n, a3); n += (y[n] != 0);
    return n;

do_two:

    y[n] = mpn_mul_1(y, x, n, p1); n += (y[n] != 0);
    y[n] = mpn_mul_1(y, y, n, p2); n += (y[n] != 0);
    return n;
}
