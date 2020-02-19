/*
    Copyright (C) 2011 Andy Novocin
    Copyright (C) 2011 William Hart
    Copyright (C) 2011 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include <stdlib.h>
#include "fmpz_poly.h"

#define TRACE 0

/*
    A degree k factor should not be considered if
        possible_degs[k] != possible_degs[0]
*/
void fmpz_poly_factor_zassenhaus_recombination(
    fmpz_poly_factor_t final_fac,
    const unsigned int * possible_degs,
	const fmpz_poly_factor_t lifted_fac, 
    const fmpz_poly_t F, const fmpz_t P, slong exp)
{
    const slong r = lifted_fac->num;
    slong tryme_deg;

    slong k, *used_arr, *sub_arr;
    fmpz_poly_t f, Q, tryme;
    fmpz *leadF;

    used_arr = flint_calloc(2 * r, sizeof(slong));
    sub_arr  = used_arr + r;

    fmpz_poly_init(f);
    fmpz_poly_init(Q);
    fmpz_poly_init(tryme);
    fmpz_poly_set(f, F);

#if TRACE == 1
    fmpz_poly_factor_print(lifted_fac); flint_printf(" lifted_fac\n");
#endif

    leadF = fmpz_poly_lead(F);

    for (k = 1; k < r; k++)
    {
        slong count = 0, indx = k - 1, l;

        for(l = 0; l < k; l++)
            sub_arr[l] = l;

        sub_arr[indx]--;
        while ((indx >= 0))
        {
            sub_arr[indx] = sub_arr[indx] + 1;

            for (l = indx + 1; l < k; l++)
                sub_arr[l] = sub_arr[l - 1] + 1;

            if (sub_arr[k - 1] > r - 1)
                indx--;
            else
            {
                for(l = 0; l < k; l++)
                {
                    if (used_arr[sub_arr[l]] == 1)
                        break;
                }

                /* try this degree only if it is possible */
                tryme_deg = 0;
                for (l = 0; l < k; l++)
                    tryme_deg += fmpz_poly_degree(lifted_fac->p + sub_arr[l]);

                if (possible_degs[tryme_deg] != possible_degs[0])
                    goto skip;

             /* Need to involve leadF, perhaps set coeff 0 to leadF and do 
                leadF * rest and check if under M_bits... here I'm using a 
                trial division... */
                fmpz_poly_set_fmpz(tryme, leadF);

                for(l = 0; l < k; l++)
                    fmpz_poly_mul(tryme, tryme, lifted_fac->p + (sub_arr[l]));

                fmpz_poly_scalar_smod_fmpz(tryme, tryme, P);
                fmpz_poly_primitive_part(tryme, tryme);

                if (fmpz_poly_divides(Q, f, tryme))
                {
                    fmpz_poly_factor_insert(final_fac, tryme, exp);

                    for(l = 0; l < k; l++)
                    {
                        used_arr[sub_arr[l]] = 1;
                        count++;
                    }

                    fmpz_poly_set(f, Q);
                    leadF = fmpz_poly_lead(f);
                 /* If r - count = k then the rest are irreducible.  
                    TODO: Add a test for that case */
                }
skip:
                indx = k - 1;
            }
        }

     /* This is where we switch to the next loop for k.  So we will have 
        found all factors using <= k local factors.  We should/could update 
        f to be the rest divided away (or multiply the remaining), could 
        also adjust r.  It is the number of remaining factors so if you 
        update then test if r = k or k+1 in which case the remaining f is 
        irreducible. */
    }

    {
        slong test = 0;

        for (k = 0; k < r; k++)
            test = test + used_arr[k];

        if (test == 0)
            fmpz_poly_factor_insert(final_fac, f, exp);
    }

    fmpz_poly_clear(f);
    fmpz_poly_clear(tryme);
    fmpz_poly_clear(Q);
    flint_free(used_arr);
}

#undef TRACE

