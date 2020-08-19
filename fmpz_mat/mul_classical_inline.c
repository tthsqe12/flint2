/*
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mat.h"
#include "longlong.h"

void
fmpz_mat_mul_classical_inline(fmpz_mat_t C, const fmpz_mat_t A,
                                                const fmpz_mat_t B)
{
    slong ar, bc, br;
    slong i, j, k;
    fmpz a, b;
    mpz_t t1, t2;
    __mpz_struct * cij;
    mp_limb_t sign, hi, lo, acc[3];
    __mpz_struct r;

    r._mp_d = acc;
    r._mp_alloc = 3;

    ar = A->r;
    br = B->r;
    bc = B->c;

    mpz_init(t1);
    mpz_init(t2);

    for (i = 0; i < ar; i++)
    {
        for (j = 0; j < bc; j++)
        {
            fmpz * cij_ref = fmpz_mat_entry(C, i, j);
            int cij_was_big = COEFF_IS_MPZ(*cij_ref);

            if (cij_was_big)
            {
                cij = COEFF_TO_PTR(*cij_ref);
            }
            else
            {
                cij = t2;
            }

            flint_mpz_set_ui(cij, UWORD(0));

            acc[2] = acc[1] = acc[0] = UWORD(0);

            for (k = 0; k < br; k++)
            {
                a = A->rows[i][k];
                b = B->rows[k][j];

                if (a == 0 || b == 0)
                    continue;

                if (!COEFF_IS_MPZ(a))   /* a is small */
                {
                    if (!COEFF_IS_MPZ(b))  /* both are small */
                    {
                        smul_ppmm(hi, lo, a, b);
                        sign = FLINT_SIGN_EXT(hi);
                        add_sssaaaaaa(acc[2], acc[1], acc[0],
                                      acc[2], acc[1], acc[0], sign, hi, lo);
                    }
                    else
                    {
                        if (a >= 0)
                            flint_mpz_addmul_ui(cij, COEFF_TO_PTR(b), a);
                        else
                            flint_mpz_submul_ui(cij, COEFF_TO_PTR(b), -a);
                    }
                }
                else if (!COEFF_IS_MPZ(b))  /* b is small */
                {
                    if (b >= 0)
                        flint_mpz_addmul_ui(cij, COEFF_TO_PTR(a), b);
                    else
                        flint_mpz_submul_ui(cij, COEFF_TO_PTR(a), -b);
                }
                else if (mpz_sgn(cij) == 0)
                {
                    mpz_mul(cij, COEFF_TO_PTR(a), COEFF_TO_PTR(b));
                }
                else if (COEFF_TO_PTR(a)->_mp_size == +1 ||
                         COEFF_TO_PTR(a)->_mp_size == -1 ||
                         COEFF_TO_PTR(b)->_mp_size == +1 ||
                         COEFF_TO_PTR(b)->_mp_size == -1)
                {
                    mpz_addmul(cij, COEFF_TO_PTR(a), COEFF_TO_PTR(b));
                }
                else
                {
                    mpz_mul(t1, COEFF_TO_PTR(a), COEFF_TO_PTR(b));
                    mpz_add(cij, cij, t1);
                }
            }

            sign = FLINT_SIGN_EXT(acc[2]);

            acc[0] ^= sign;
            acc[1] ^= sign;
            acc[2] ^= sign;

            r._mp_size = acc[2] ? 3 : (acc[1] ? 2 : acc[0] != 0);

            if (sign)
            {
                r._mp_size = -r._mp_size;
                sub_dddmmmsss(acc[2], acc[1], acc[0],
                              acc[2], acc[1], acc[0], sign, sign, sign);
            }

            if (mpz_sgn(cij) == 0)
            {
                if (cij_was_big)
                {
                    mpz_set(cij, &r);
                    _fmpz_demote_val(cij_ref);
                }
                else
                {
                    fmpz_set_mpz(cij_ref, &r);
                }
            }
            else
            {
                mpz_add(cij, cij, &r);            

                if (cij_was_big)
                    _fmpz_demote_val(cij_ref);
                else
                    fmpz_set_mpz(cij_ref, cij);
            }
        }
    }

    mpz_clear(t1);
    mpz_clear(t2);
}
