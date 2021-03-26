/*
    Copyright 2021 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include <stdlib.h>
#include "flint.h"
#include "fmpz_mat.h"
#include "profiler.h"

int main(void)
{
    slong dim, i, reps, total, mint, maxt;
    flint_bitcnt_t Abits, Bbits, Cbits;
    timeit_t timer;
    FLINT_TEST_INIT(state);

    flint_set_num_threads(1);

    for (dim = 16; dim <= 200; dim += 5 + dim/8)
    {
        fmpz_mat_t A, B, C, D, E;

        fmpz_mat_init(A, dim, dim);
        fmpz_mat_init(B, dim, dim);
        fmpz_mat_init(C, dim, dim);
        fmpz_mat_init(D, dim, dim);
        fmpz_mat_init(E, dim, dim);

        total = 0;
        mint = 10000000000;
        maxt = 0;

        for (Abits = 20; Abits < 2*FLINT_BITS; Abits += 15)
        for (Bbits = Abits; Bbits < 2*FLINT_BITS; Bbits += 15)
        {
            Cbits = Abits + Bbits + FLINT_BIT_COUNT(dim) + 1;

            if (Cbits >= 4*FLINT_BITS)
                continue;

            fmpz_mat_randtest(A, state, Abits);
            fmpz_mat_randtest(B, state, Bbits);

            fmpz_mat_mul_classical_inline(D, A, B);

            reps = 1 + 50000000/dim/dim/dim;

            timeit_start(timer);
            for (i = reps; i > 0; i--)
                fmpz_mat_mul_multi_mod(E, A, B);
            timeit_stop(timer);
/*
            flint_printf("dim %3wd, Abits %3wu, Bbits %3wu: %5wd\n",
                         dim, Abits, Bbits, timer->wall);
*/
            total += timer->wall;
            mint = FLINT_MIN(mint, timer->wall);
            maxt = FLINT_MAX(maxt, timer->wall);

            if (!fmpz_mat_equal(D, E))
            {
                flint_printf("oops\n");
                flint_abort();
            }
        }

        flint_printf("dim %3wd: min %4wd  max %4wd  total %6wd\n", dim, mint, maxt, total);

        fmpz_mat_clear(A);
        fmpz_mat_clear(B);
        fmpz_mat_clear(C);
        fmpz_mat_clear(D);
        fmpz_mat_clear(E);
    }

    FLINT_TEST_CLEANUP(state);
    return 0;
}
