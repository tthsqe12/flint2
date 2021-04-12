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
#include "profiler.h"
#include "flint.h"
#include "fmpz_mat.h"
#include "fmpz.h"
#include "ulong_extras.h"



static void fmpz_mat_pascal(fmpz_mat_t mat, int triangular)
{
    slong R, C, i, j;

    R = fmpz_mat_nrows(mat);
    C = fmpz_mat_ncols(mat);

    if (R == 0 || C == 0)
        return;

    if (triangular > 0)
    {
        for (i = 1; i < R; i++)
            for (j = 0; j < i && j < C; j++)
                fmpz_zero(fmpz_mat_entry(mat, i, j));

        for (j = 0; j < C; j++)
            fmpz_one(fmpz_mat_entry(mat, 0, j));

        for (i = 1; i < R && i < C; i++)
            fmpz_one(fmpz_mat_entry(mat, i, i));

        for (i = 1; i < R; i++)
            for (j = i + 1; j < C; j++)
                fmpz_add(fmpz_mat_entry(mat, i, j),
                    fmpz_mat_entry(mat, i, j - 1),
                    fmpz_mat_entry(mat, i - 1, j - 1));
    }
    else if (triangular < 0)
    {
        for (i = 0; i < R; i++)
            for (j = i + 1; j < C; j++)
                fmpz_zero(fmpz_mat_entry(mat, i, j));

        for (i = 0; i < R; i++)
            fmpz_one(fmpz_mat_entry(mat, i, 0));

        for (i = 1; i < R && i < C; i++)
            fmpz_one(fmpz_mat_entry(mat, i, i));

        for (i = 2; i < R; i++)
            for (j = 1; j < i && j < C; j++)
                fmpz_add(fmpz_mat_entry(mat, i, j),
                    fmpz_mat_entry(mat, i - 1, j - 1),
                    fmpz_mat_entry(mat, i - 1, j));
    }
    else
    {
        for (j = 0; j < C; j++)
            fmpz_one(fmpz_mat_entry(mat, 0, j));

        for (i = 1; i < R; i++)
            fmpz_one(fmpz_mat_entry(mat, i, 0));

        for (i = 1; i < R; i++)
            for (j = 1; j < C; j++)
                fmpz_add(fmpz_mat_entry(mat, i, j),
                    fmpz_mat_entry(mat, i, j - 1),
                    fmpz_mat_entry(mat, i - 1, j));
    }
}

slong rep_time(timeit_t timer, fmpz_mat_t C, fmpz_mat_t A, fmpz_mat_t B)
{
    slong i, j, reps = 1;

    timeit_start(timer);
    fmpz_mat_mul(C, A, B);
    if (timeit_query_wall(timer) < 300)
    {
        j = 2 + 300/(1 + timeit_query_wall(timer));
        do {
            reps += j;
            for (i = j; i > 0; i--)
                fmpz_mat_mul(C, A, B);
        } while (timeit_query_wall(timer) < 300);
    }
    timeit_stop(timer);

    return reps;
}

int main(void)
{
    slong t, tmul = 300;
    slong m, k, n, reps;
    flint_bitcnt_t Abits, Bbits;
    fmpz_mat_t A, B, C;
    timeit_t timer;
    FLINT_TEST_INIT(state);

    flint_set_num_threads(1);

    for (t = 0; t < 0*tmul; t++)
    {
        if (t % 10 == 0)
            fprintf(stderr, "#randtest %d/%d\n", (int)t, (int)tmul);

        m = 2 + n_randint(state, 500);
        k = 3 + n_randint(state, 500);
        n = 2 + n_randint(state, 500);
        Abits = 5 + n_randint(state, 1000);
        Bbits = 5 + n_randint(state, 1000);

        fmpz_mat_init(A, m, k);
        fmpz_mat_init(B, k, n);
        fmpz_mat_init(C, m, n);

        fmpz_mat_randtest(A, state, Abits);
        fmpz_mat_randtest(B, state, Bbits);

        reps = rep_time(timer, C, A, B);

        flint_printf("randtest[%wd, %wd, %wd, %wd, %wd] = %.3f #ns %wd reps\n",
                      m, k, n, Abits, Bbits, timer->wall*1000000.0/reps/m/k/n, reps);

        fmpz_mat_clear(A);
        fmpz_mat_clear(B);
        fmpz_mat_clear(C);
    }

    for (t = 0; t < 1+0*tmul; t++)
    {
        if (t % 10 == 0)
            fprintf(stderr, "#randtest square %d/%d\n", (int)t, (int)tmul);

        m = k = n = 3 + n_randint(state, 1500);
        Abits = Bbits = 5 + n_randint(state, 1 + 100000/(10 + m));

m = k = n = 773;
Abits = Bbits = 7;

        fmpz_mat_init(A, m, k);
        fmpz_mat_init(B, k, n);
        fmpz_mat_init(C, m, n);

        fmpz_mat_randtest(A, state, Abits);
        fmpz_mat_randtest(B, state, Bbits);

        reps = rep_time(timer, C, A, B);

        flint_printf("randtest[%wd, %wd, %wd, %wd, %wd] = %.3f #ns %wd reps\n",
                      m, k, n, Abits, Bbits, timer->wall*1000000.0/reps/m/k/n, reps);

        fmpz_mat_clear(A);
        fmpz_mat_clear(B);
        fmpz_mat_clear(C);
    }

    for (t = 0; t < 0*tmul; t++)
    {
        int atri, btri;

        if (t % 10 == 0)
            fprintf(stderr, "#pascal %d/%d\n", (int)t, (int)tmul);

        m = 2 + n_randint(state, 600);
        k = 3 + n_randint(state, 600);
        n = 2 + n_randint(state, 600);

        fmpz_mat_init(A, m, k);
        fmpz_mat_init(B, k, n);
        fmpz_mat_init(C, m, n);

        atri = -1 + (int)n_randint(state, 3);
        btri = -1 + (int)n_randint(state, 3);

        fmpz_mat_pascal(A, atri);
        fmpz_mat_pascal(B, btri);

        reps = rep_time(timer, C, A, B);

        flint_printf("pascal[%wd, %wd, %wd, %d, %d] = %.3f #ns %wd reps\n",
                      m, k, n, atri, btri, timer->wall*1000000.0/reps/m/k/n, reps);

        fmpz_mat_clear(A);
        fmpz_mat_clear(B);
        fmpz_mat_clear(C);
    }

    FLINT_TEST_CLEANUP(state);
    return 0;
}
