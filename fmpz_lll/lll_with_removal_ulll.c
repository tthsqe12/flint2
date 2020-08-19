/*
    Copyright (C) 2009, 2010 William Hart
    Copyright (C) 2009, 2010 Andy Novocin
    Copyright (C) 2014 Abhinav Baid

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_lll.h"


void
my_fmpz_mat_mul_classical_inline(fmpz_mat_t C, const fmpz_mat_t A,
                                         const fmpz_mat_t B, mpz_t tt)
{
    slong ar, bc, br;
    slong i, j, k;
    __mpz_struct * t;

    fmpz a, b;

    mp_limb_t au, bu;
    mp_limb_t pos[3];
    mp_limb_t neg[3];

    ar = A->r;
    br = B->r;
    bc = B->c;

    for (i = 0; i < ar; i++)
    {
        for (j = 0; j < bc; j++)
        {
            t = _fmpz_promote(fmpz_mat_entry(C, i, j));

            flint_mpz_set_ui(t, UWORD(0));

            pos[2] = pos[1] = pos[0] = neg[2] = neg[1] = neg[0] = UWORD(0);

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
                        au = FLINT_ABS(a);
                        bu = FLINT_ABS(b);

                        umul_ppmm(au, bu, au, bu);

                        if ((a ^ b) >= WORD(0))
                            add_sssaaaaaa(pos[2], pos[1], pos[0],
                                          pos[2], pos[1], pos[0], 0, au, bu);
                        else
                            add_sssaaaaaa(neg[2], neg[1], neg[0],
                                          neg[2], neg[1], neg[0], 0, au, bu);
                    }
                    else
                    {
                        if (a >= 0)
                            flint_mpz_addmul_ui(t, COEFF_TO_PTR(b), a);
                        else
                            flint_mpz_submul_ui(t, COEFF_TO_PTR(b), -a);
                    }
                }
                else if (!COEFF_IS_MPZ(b))  /* b is small */
                {
                    if (b >= 0)
                        flint_mpz_addmul_ui(t, COEFF_TO_PTR(a), b);
                    else
                        flint_mpz_submul_ui(t, COEFF_TO_PTR(a), -b);
                }
                else
                {
                    if (mpz_sgn(t) == 0)
                    {
                        mpz_mul(t, COEFF_TO_PTR(a), COEFF_TO_PTR(b));
                    }
                    else
                    {
                        mpz_mul(tt, COEFF_TO_PTR(a), COEFF_TO_PTR(b));
                        mpz_add(t, t, tt);
/*
                        mpz_addmul(t, COEFF_TO_PTR(a), COEFF_TO_PTR(b));
*/
                    }
                }
            }

            if (mpz_sgn(t) != 0 || pos[2] || neg[2] || pos[1] || neg[1])
            {
                __mpz_struct r;

                r._mp_size = pos[2] ? 3 : (pos[1] ? 2 : pos[0] != 0);
                r._mp_alloc = r._mp_size;
                r._mp_d = pos;

                mpz_add(t, t, &r);

                r._mp_size = neg[2] ? 3 : (neg[1] ? 2 : neg[0] != 0);
                r._mp_alloc = r._mp_size;
                r._mp_d = neg;

                mpz_sub(t, t, &r);
            }
            else
            {
                if (neg[0] > pos[0])
                {
                    flint_mpz_set_ui(t, neg[0] - pos[0]);
                    mpz_neg(t, t);
                }
                else
                {
                    flint_mpz_set_ui(t, pos[0] - neg[0]);
                }
            }

            _fmpz_demote_val(fmpz_mat_entry(C, i, j));

        }
    }
}


void
my_fmpz_mat_mul(fmpz_mat_t C, const fmpz_mat_t A, const fmpz_mat_t B,
mpz_t t1)
{
    slong ar, br, bc;
    slong abits, bbits, bits;
    slong i, j, dim;

    ar = fmpz_mat_nrows(A);
    br = fmpz_mat_nrows(B);
    bc = fmpz_mat_ncols(B);

    if (ar == 0 || br == 0 || bc == 0)
    {
        fmpz_mat_zero(C);
        return;
    }

    if (C == A || C == B)
    {
        fmpz_mat_t T;
        fmpz_mat_init(T, ar, bc);
        fmpz_mat_mul(T, A, B);
        fmpz_mat_swap(C, T);
        fmpz_mat_clear(T);
        return;
    }

    if (br == 1)
    {
        for (i = 0; i < ar; i++)
            for (j = 0; j < bc; j++)
                fmpz_mul(fmpz_mat_entry(C, i, j),
                    fmpz_mat_entry(A, i, 0), fmpz_mat_entry(B, 0, j));
        return;
    }

    if (br == 2)
    {
        for (i = 0; i < ar; i++)
            for (j = 0; j < bc; j++)
                fmpz_fmma(fmpz_mat_entry(C, i, j),
                    fmpz_mat_entry(A, i, 0), fmpz_mat_entry(B, 0, j),
                    fmpz_mat_entry(A, i, 1), fmpz_mat_entry(B, 1, j));
        return;
    }

    abits = fmpz_mat_max_bits(A);
    bbits = fmpz_mat_max_bits(B);
    abits = FLINT_ABS(abits);
    bbits = FLINT_ABS(bbits);

    if (abits == 0 || bbits == 0)
    {
        fmpz_mat_zero(C);
        return;
    }

    bits = abits + bbits + FLINT_BIT_COUNT(br) + 1;
    dim = FLINT_MIN(ar, bc);
    dim = FLINT_MIN(dim, br);

    if (bits <= FLINT_BITS - 2)
    {
        if ((dim > 160 && abits + bbits <= 20) || dim > 600) /* tuning param */
            _fmpz_mat_mul_multi_mod(C, A, B, bits);
        else if (dim > 160) /* tuning param */
            fmpz_mat_mul_strassen(C, A, B);
        else
            fmpz_mat_mul_1(C, A, B);
    }
    else if (abits <= FLINT_BITS - 2 && bbits <= FLINT_BITS - 2)
    {
        if (dim > 400) /* tuning param */
            _fmpz_mat_mul_multi_mod(C, A, B, bits);
        else if (bits <= 2 * FLINT_BITS - 1)
            fmpz_mat_mul_2a(C, A, B);
        else
            fmpz_mat_mul_2b(C, A, B);
    }
    else if (abits <= 2 * FLINT_BITS && bbits <= 2 * FLINT_BITS
                                     && bits <= 4 * FLINT_BITS - 1)
    {
        if (dim > 40) /* tuning param */
            _fmpz_mat_mul_multi_mod(C, A, B, bits);
        else
            fmpz_mat_mul_4(C, A, B);
    }
    else
    {
        if (dim >= 3 * FLINT_BIT_COUNT(bits))  /* tuning param */
            _fmpz_mat_mul_multi_mod(C, A, B, bits);
        else if (abits >= 500 && bbits >= 500 && dim >= 8)  /* tuning param */
            fmpz_mat_mul_strassen(C, A, B);
        else
            my_fmpz_mat_mul_classical_inline(C, A, B, t1);
    }
}


int
fmpz_lll_with_removal_ulll(fmpz_mat_t FM, fmpz_mat_t UM, slong new_size,
                           const fmpz_t gs_B, const fmpz_lll_t fl)
{

slong loop_count = 0;

    int newd;
    if (fl->rt == Z_BASIS)
    {
        slong r, c, mbits, prev_mbits, i, j;
        int full_prec = 1, done = 0, is_U_I;
        fmpz_mat_t U, big_td, trunc_data;
mpz_t t1, t2;
fmpz_mat_t T;

flint_printf("doing fl->rt == Z_BASIS\n");

        r = FM->r;
        c = FM->c;
        mbits = FLINT_ABS(fmpz_mat_max_bits(FM));
        prev_mbits = mbits;

mpz_init(t1);
mpz_init(t2);
fmpz_mat_init(T, r, c);

        fmpz_mat_init(big_td, r, c + r);
        fmpz_mat_init(trunc_data, r, c);

flint_printf("mbits: %wd, new_size: %wd\n", mbits, new_size);

        if (mbits > new_size)
        {
            full_prec = 0;

            /* do some truncating */
            fmpz_mat_scalar_tdiv_q_2exp(trunc_data, FM,
                                        (ulong) (mbits - new_size));

            /* make a large lattice which has identity in one corner and trunc_data in the other */
            for (i = 0; i < r; i++)
            {
                fmpz_one(fmpz_mat_entry(big_td, i, i));

                for (j = r; j < r + c; j++)
                    fmpz_set(fmpz_mat_entry(big_td, i, j),
                             fmpz_mat_entry(trunc_data, i, j - r));
            }
        }
        else
        {
            full_prec = 1;
        }


        while (done == 0)
        {
loop_count++;
/*
flint_printf("loop_count: %wd\n", loop_count);
*/
            if (full_prec == 0)
            {
                fmpz_lll_wrapper_with_removal_knapsack(big_td, UM, gs_B, fl);
            }
            else
            {
                newd =
                    fmpz_lll_wrapper_with_removal_knapsack(FM, UM, gs_B, fl);
            }

            if (full_prec == 1)
                done = 1;
            else
            {
                /* get U and compare it to the identity */
                fmpz_mat_window_init(U, big_td, 0, 0, r, r);
                is_U_I = fmpz_mat_is_one(U);

                if (is_U_I == 0)
                {
                    my_fmpz_mat_mul(T, U, FM, t1);
                    fmpz_mat_swap(FM, T);
                }

                mbits = FLINT_ABS(fmpz_mat_max_bits(FM));
                /* make this condition better? */
                if ((mbits - new_size > 0)
                    && (mbits <= prev_mbits - (slong) (new_size / 4))
                    && is_U_I == 0)
                {
                    /* do some truncating */
                    fmpz_mat_scalar_tdiv_q_2exp(trunc_data, FM,
                                                (ulong) (mbits - new_size));

                    /* keep with the big_td concept */
                    for (i = 0; i < r; i++)
                    {
                        for (j = 0; j < i; j++)
                            fmpz_zero(fmpz_mat_entry(big_td, i, j));
                        fmpz_one(fmpz_mat_entry(big_td, i, i));

                        for (j = i + 1; j < r; j++)
                            fmpz_zero(fmpz_mat_entry(big_td, i, j));

                        for (j = r; j < r + c; j++)
                            fmpz_set(fmpz_mat_entry
                                     (big_td, i, j),
                                     fmpz_mat_entry(trunc_data, i, j - r));
                    }
                }
                else
                {
flint_printf("switching to full_prec = 1 at loop_count = %wd\n", loop_count);
                    /* can switch to FM, no need for a new identity */
                    full_prec = 1;
                }

                prev_mbits = mbits;
                fmpz_mat_window_clear(U);
            }
        }

flint_printf("loop_count: %wd\n", loop_count);


flint_printf("fmpz_mat_max_bits(FM): %wd\n", fmpz_mat_max_bits(FM));

        fmpz_mat_clear(trunc_data);
        fmpz_mat_clear(big_td);

fmpz_mat_clear(T);
mpz_clear(t1);
mpz_clear(t2);

    }
    else
    {
flint_printf("calling fmpz_lll_wrapper_with_removal_knapsack\n");


        newd = fmpz_lll_wrapper_with_removal_knapsack(FM, UM, gs_B, fl);
    }

    return newd;
}
