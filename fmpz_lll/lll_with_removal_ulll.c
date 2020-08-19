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


int
fmpz_lll_with_removal_ulll(fmpz_mat_t FM, fmpz_mat_t UM, slong new_size,
                           const fmpz_t gs_B, const fmpz_lll_t fl)
{
    int newd;
    if (fl->rt == Z_BASIS)
    {
        slong r, c, mbits, prev_mbits, i, j;
        int full_prec = 1, done = 0, is_U_I;
        fmpz_mat_t U, big_td, T;

        r = FM->r;
        c = FM->c;
        mbits = FLINT_ABS(fmpz_mat_max_bits(FM));
        prev_mbits = mbits;

        fmpz_mat_init(T, r, c);
        fmpz_mat_init(big_td, r, c + r);

        if (mbits > new_size)
        {
            full_prec = 0;

            /* make a large lattice which has identity on the left
               and truncated FM on the right
            */
            for (i = 0; i < r; i++)
            {
                fmpz_one(fmpz_mat_entry(big_td, i, i));

                for (j = r; j < r + c; j++)
                    fmpz_tdiv_q_2exp(fmpz_mat_entry(big_td, i, j),
                               fmpz_mat_entry(FM, i, j - r), mbits - new_size);
            }
        }
        else
        {
            full_prec = 1;
        }

        while (done == 0)
        {
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
                    fmpz_mat_mul(T, U, FM);
                    fmpz_mat_swap(FM, T);
                }

                mbits = FLINT_ABS(fmpz_mat_max_bits(FM));
                /* make this condition better? */
                if ((mbits - new_size > 0)
                    && (mbits <= prev_mbits - (slong) (new_size / 4))
                    && is_U_I == 0)
                {
                    /* keep with the big_td concept */
                    for (i = 0; i < r; i++)
                    {
                        for (j = 0; j < r; j++)
                            fmpz_set_si(fmpz_mat_entry(big_td, i, j), i == j);

                        for (j = 0; j < c; j++)
                            fmpz_tdiv_q_2exp(fmpz_mat_entry(big_td, i, j + r),
                               fmpz_mat_entry(FM, i, j), mbits - new_size);
                    }
                }
                else
                {
                    /* can switch to FM, no need for a new identity */
                    full_prec = 1;
                }

                prev_mbits = mbits;
                fmpz_mat_window_clear(U);
            }
        }

        fmpz_mat_clear(big_td);
        fmpz_mat_clear(T);
    }
    else
    {
        newd = fmpz_lll_wrapper_with_removal_knapsack(FM, UM, gs_B, fl);
    }

    return newd;
}
