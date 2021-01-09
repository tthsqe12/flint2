/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_vec.h"
#include "fmpz_mod_mpoly.h"

void _fmpz_mod_vec_neg(fmpz * A, const fmpz * B, slong len,
                                                      const fmpz_mod_ctx_t ctx)
{
    for (len--; len >= 0; len--)
        fmpz_mod_neg(A + len, B + len, ctx);
}

void _fmpz_mod_vec_scalar_mul_fmpz_mod(
    fmpz * A,
    const fmpz * B,
    slong len,
    const fmpz_t c,
    const fmpz_mod_ctx_t ctx)
{
    if (fmpz_is_one(c))
    {
        _fmpz_vec_set(A, B, len);
    }
    else
    {
        for (len--; len >= 0; len--)
            fmpz_mod_mul(A + len, B + len, c, ctx);
    }
}

void _fmpz_mod_vec_scalar_div_fmpz_mod(
    fmpz * A,
    const fmpz * B,
    slong len,
    const fmpz_t c,
    const fmpz_mod_ctx_t ctx)
{
    fmpz_t d;
    fmpz_init(d);
    fmpz_mod_inv(d, c, ctx);
    for (len--; len >= 0; len--)
        fmpz_mod_mul(A + len, B + len, d, ctx);
    fmpz_clear(d);
}

void _fmpz_mod_vec_dot(
    fmpz_t d,
    const fmpz * a,
    const fmpz * b,
    slong len,
    const fmpz_mod_ctx_t ctx)
{
    fmpz_zero(d);
    for (len--; len >= 0; len--)
        fmpz_addmul(d, a + len, b + len);
    fmpz_mod(d, d, fmpz_mod_ctx_modulus(ctx));
}