/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mod_mpoly_factor.h"
#include "mpn_extras.h"
#include "fmpz_vec.h"

void fmpz_mod_polyun_content_poly(
    fmpz_mod_poly_t g,
    const fmpz_mod_polyun_t A,
    const fmpz_mod_ctx_t ctx)
{
    slong i;
    fmpz_mod_poly_zero(g, ctx);
    for (i = 0; i < A->length; i++)
        fmpz_mod_poly_gcd(g, g, A->terms[i].coeff, ctx);
}

void fmpz_mod_polyun_divexact_poly(
    fmpz_mod_polyun_t A,
    const fmpz_mod_poly_t g,
    const fmpz_mod_ctx_t ctx)
{
    slong i;
    for (i = 0; i < A->length; i++)
        fmpz_mod_poly_div(A->terms[i].coeff, A->terms[i].coeff, g, ctx);
}

void fmpz_mod_polyun_mul_poly(
    fmpz_mod_polyun_t A,
    const fmpz_mod_poly_t g,
    const fmpz_mod_ctx_t ctx)
{
    slong i;
    for (i = 0; i < A->length; i++)
        fmpz_mod_poly_mul(A->terms[i].coeff, A->terms[i].coeff, g, ctx);
}


slong fmpz_mod_polyun_lastdeg(
    const fmpz_mod_polyun_t A,
    const fmpz_mod_ctx_t ctx)
{
    slong i, len = 0;
    for (i = 0; i < A->length; i++)
        len = FLINT_MAX(len, A->terms[i].coeff->length);
    return len - 1;
}

void fmpz_mod_polyun_one(fmpz_mod_polyun_t A, const fmpz_mod_ctx_t ctx)
{
    fmpz_mod_polyun_fit_length(A, 1, ctx);
    fmpz_mod_poly_one(A->terms[0].coeff, ctx);
    A->terms[0].exp = 0;
    A->length = 1;
}

