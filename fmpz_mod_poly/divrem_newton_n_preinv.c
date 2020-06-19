/*
    Copyright (C) 2013 Martin Lee

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#undef ulong
#define ulong ulongxx/* interferes with system includes */

#include <stdlib.h>

#undef ulong

#include <gmp.h>

#define ulong mp_limb_t

#include "fmpz_vec.h"
#include "fmpz_mod_poly.h"

void _fmpz_mod_poly_divrem_newton_n_preinv (fmpz* Q, fmpz* R, const fmpz* A,
                              slong lenA, const fmpz* B, slong lenB,
                              const fmpz* Binv, slong lenBinv, const fmpz_t mod)
{
    const slong lenQ = lenA - lenB + 1;

    _fmpz_mod_poly_div_newton_n_preinv(Q, A, lenA, B, lenB, Binv, lenBinv, mod);

    if (lenB > 1)
    {
        if (lenQ >= lenB - 1)
            _fmpz_mod_poly_mullow(R, Q, lenQ, B, lenB - 1, mod, lenB - 1);
        else
            _fmpz_mod_poly_mullow(R, B, lenB - 1, Q, lenQ, mod, lenB - 1);

        _fmpz_vec_sub(R, A, R, lenB - 1);
        _fmpz_vec_scalar_mod_fmpz(R, R, lenB - 1, mod);
    }
}

void fmpz_mod_poly_divrem_newton_n_preinv(fmpz_mod_poly_t Q, fmpz_mod_poly_t R,
                              const fmpz_mod_poly_t A, const fmpz_mod_poly_t B,
                          const fmpz_mod_poly_t Binv, const fmpz_mod_ctx_t ctx)
{
    const slong lenA = A->length, lenB = B->length, lenBinv= Binv->length,
                       lenQ = lenA - lenB + 1;
    fmpz* q, * r;

    if (lenB == 0)
    {
        if (fmpz_is_one(fmpz_mod_ctx_modulus(ctx)))
        {
            fmpz_mod_poly_set(Q, A);
            fmpz_mod_poly_zero(R);
            return;
        } else
        {
            flint_printf("Exception (fmpz_mod_poly_divrem_newton_n_preinv)."
                         " Division by zero.\n");
            flint_abort();
        }
    }

    if (lenA < lenB)
    {
        fmpz_mod_poly_set(R, A);
        fmpz_mod_poly_zero(Q);
        return;
    }

    if (lenA > 2 * lenB - 2)
    {
        flint_printf ("Exception (fmpz_mod_poly_divrem_newton_n_preinv).\n");
        flint_abort();
    }

    if (Q == A || Q == B || Q == Binv)
    {
        q = _fmpz_vec_init(lenQ);
    }
    else
    {
        fmpz_mod_poly_fit_length(Q, lenQ);
        q = Q->coeffs;
    }
    if (R == A || R == B || R == Binv)
    {
        r = _fmpz_vec_init(lenB - 1);
    }
    else
    {
        fmpz_mod_poly_fit_length(R, lenB - 1);
        r = R->coeffs;
    }

    _fmpz_mod_poly_divrem_newton_n_preinv (q, r, A->coeffs, lenA,
                                           B->coeffs, lenB, Binv->coeffs,
                                           lenBinv, fmpz_mod_ctx_modulus(ctx));

    if (Q == A || Q == B || Q == Binv)
    {
        _fmpz_vec_clear(Q->coeffs, Q->alloc);
        Q->coeffs = q;
        Q->alloc  = lenQ;
        Q->length = lenQ;
    }
    else
    {
        _fmpz_mod_poly_set_length(Q, lenQ);
    }
    if (R == A || R == B || R == Binv)
    {
        _fmpz_vec_clear(R->coeffs, R->alloc);
        R->coeffs = r;
        R->alloc  = lenB - 1;
        R->length = lenB - 1;
    }
    else
    {
      _fmpz_mod_poly_set_length(R, lenB - 1);
    }

    _fmpz_mod_poly_normalise(R);
}
