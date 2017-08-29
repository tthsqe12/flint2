/*
    Copyright (C) 2017 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "profiler.h"
#include "fmpz.h"
#include "fmpz_mpoly.h"
#include "ulong_extras.h"

int
main(void)
{
    int i, j, result, max_threads = 5;
    FLINT_TEST_INIT(state);

    flint_printf("mul_heap_threaded....");
    fflush(stdout);

/*
    slong power;
    fmpz_mpoly_ctx_t ctx;
    fmpz_mpoly_t f, g, h, h1, X, Y, Z, T, U;
    timeit_t time;


    fmpz_mpoly_ctx_init(ctx, 5, ORD_LEX);

    fmpz_mpoly_init(f, ctx);
    fmpz_mpoly_init(g, ctx);
    fmpz_mpoly_init(h, ctx);
    fmpz_mpoly_init(h1, ctx);

    fmpz_mpoly_init(X, ctx);
    fmpz_mpoly_init(Y, ctx);
    fmpz_mpoly_init(Z, ctx);
    fmpz_mpoly_init(T, ctx);
    fmpz_mpoly_init(U, ctx);

    fmpz_mpoly_gen(X, 0, ctx);
    fmpz_mpoly_gen(Y, 1, ctx);
    fmpz_mpoly_gen(Z, 2, ctx);
    fmpz_mpoly_gen(T, 3, ctx);
    fmpz_mpoly_gen(U, 4, ctx);

    fmpz_mpoly_set_si(f, WORD(1), ctx);

    fmpz_mpoly_pow_fps(h, X, WORD(1), ctx);
    fmpz_mpoly_scalar_mul_si(h, h, WORD(1), ctx);
    fmpz_mpoly_add(f, f, h, ctx);

    fmpz_mpoly_pow_fps(h, Y, WORD(1), ctx);
    fmpz_mpoly_scalar_mul_si(h, h, WORD(1), ctx);
    fmpz_mpoly_add(f, f, h, ctx);

    fmpz_mpoly_pow_fps(h, Z, WORD(2), ctx);
    fmpz_mpoly_scalar_mul_si(h, h, WORD(2), ctx);
    fmpz_mpoly_add(f, f, h, ctx);

    fmpz_mpoly_pow_fps(h, T, WORD(3), ctx);
    fmpz_mpoly_scalar_mul_si(h, h, WORD(3), ctx);
    fmpz_mpoly_add(f, f, h, ctx);

    fmpz_mpoly_pow_fps(h, U, WORD(5), ctx);
    fmpz_mpoly_scalar_mul_si(h, h, WORD(5), ctx);
    fmpz_mpoly_add(f, f, h, ctx);


    fmpz_mpoly_set_si(g, WORD(1), ctx);

    fmpz_mpoly_pow_fps(h, U, WORD(1), ctx);
    fmpz_mpoly_scalar_mul_si(h, h, WORD(1), ctx);
    fmpz_mpoly_add(g, g, h, ctx);

    fmpz_mpoly_pow_fps(h, T, WORD(1), ctx);
    fmpz_mpoly_scalar_mul_si(h, h, WORD(1), ctx);
    fmpz_mpoly_add(g, g, h, ctx);

    fmpz_mpoly_pow_fps(h, Z, WORD(2), ctx);
    fmpz_mpoly_scalar_mul_si(h, h, WORD(2), ctx);
    fmpz_mpoly_add(g, g, h, ctx);

    fmpz_mpoly_pow_fps(h, Y, WORD(3), ctx);
    fmpz_mpoly_scalar_mul_si(h, h, WORD(3), ctx);
    fmpz_mpoly_add(g, g, h, ctx);

    fmpz_mpoly_pow_fps(h, X, WORD(5), ctx);
    fmpz_mpoly_scalar_mul_si(h, h, WORD(5), ctx);
    fmpz_mpoly_add(g, g, h, ctx);

    printf("f = "); fmpz_mpoly_print_pretty(f, NULL, ctx); printf("\n");
    printf("g = "); fmpz_mpoly_print_pretty(g, NULL, ctx); printf("\n");

    power = 2;
    fmpz_mpoly_pow_fps(f, f, power, ctx);
    fmpz_mpoly_pow_fps(g, g, power, ctx);
    flint_printf("power: %wd\n", power);

    flint_printf("mul_johnson          ... ");
    timeit_start(time);
    fmpz_mpoly_mul_johnson(h1, f, g, ctx);
    timeit_stop(time);
    flint_printf("wall %wd\n", time->wall);

    flint_printf("mul_johnson          ... ");
    timeit_start(time);
    fmpz_mpoly_mul_johnson(h, f, g, ctx);
    timeit_stop(time);
    flint_printf("wall %wd\n", time->wall);

    for (i = 1; i <= max_threads; i++)
    {
        flint_set_num_threads(i);
        flint_printf("mul_heap_threadedC %wd ... ",i);
        fflush(stdout);
        timeit_start(time);
        fmpz_mpoly_mul_heap_threaded(h, f, g, ctx);
        timeit_stop(time);
        flint_printf("wall %wd\n", time->wall);
        fmpz_mpoly_test(h, ctx);
        if (!fmpz_mpoly_equal(h,h1,ctx))
            printf("h != h1\n");
    }

    fmpz_mpoly_clear(U, ctx);
    fmpz_mpoly_clear(T, ctx);
    fmpz_mpoly_clear(Z, ctx);
    fmpz_mpoly_clear(Y, ctx);
    fmpz_mpoly_clear(X, ctx);
    fmpz_mpoly_clear(h1, ctx);
    fmpz_mpoly_clear(h, ctx);
    fmpz_mpoly_clear(g, ctx);
    fmpz_mpoly_clear(f, ctx);
*/

    /* Check mul_heap_threaded matches mul_johnson */
    for (i = 0; i < 10 * flint_test_multiplier(); i++)
    {
        fmpz_mpoly_ctx_t ctx;
        fmpz_mpoly_t f, g, h, k;
        ordering_t ord;
        slong nvars, len, len1, len2, exp_bound, exp_bound1, exp_bound2;
        slong coeff_bits, exp_bits, exp_bits1, exp_bits2;

        ord = mpoly_ordering_randtest(state);
        nvars = n_randint(state, 10) + 1;

        fmpz_mpoly_ctx_init(ctx, nvars, ord);

        fmpz_mpoly_init(f, ctx);
        fmpz_mpoly_init(g, ctx);
        fmpz_mpoly_init(h, ctx);
        fmpz_mpoly_init(k, ctx);

        len = n_randint(state, 100);
        len1 = n_randint(state, 100);
        len2 = n_randint(state, 100);

        exp_bits = n_randint(state, 20/(nvars + 
                            mpoly_ordering_isdeg(ord) + (nvars == 1)) + 1) + 1;
        exp_bits1 = n_randint(state, 20/(nvars + 
                            mpoly_ordering_isdeg(ord) + (nvars == 1)) + 1) + 1;
        exp_bits2 = n_randint(state, 20/(nvars + 
                            mpoly_ordering_isdeg(ord) + (nvars == 1)) + 1) + 1;
        exp_bound = n_randbits(state, exp_bits);
        exp_bound1 = n_randbits(state, exp_bits1);
        exp_bound2 = n_randbits(state, exp_bits2);

        coeff_bits = n_randint(state, 200);

        for (j = 0; j < 4; j++)
        {
            fmpz_mpoly_randtest(f, state, len1, exp_bound1, coeff_bits, ctx);
            fmpz_mpoly_randtest(g, state, len2, exp_bound2, coeff_bits, ctx);
            fmpz_mpoly_randtest(h, state, len, exp_bound, coeff_bits, ctx);
            fmpz_mpoly_randtest(k, state, len, exp_bound, coeff_bits, ctx);

            flint_set_num_threads(n_randint(state, max_threads) + 1);

            fmpz_mpoly_mul_johnson(h, f, g, ctx);
            fmpz_mpoly_test(h, ctx);

            fmpz_mpoly_mul_heap_threaded(k, f, g, ctx);
            fmpz_mpoly_test(k, ctx);

            result = fmpz_mpoly_equal(h, k, ctx);

            if (!result)
            {
                printf("FAIL\n");

                printf("ord = "); mpoly_ordering_print(ord);
                printf(", len = %ld, exp_bits = %ld, exp_bound = %lx, "
                    "len1 = %ld, exp_bits1 = %ld, exp_bound1 = %lx, "
                    "len2 = %ld, exp_bits2 = %ld, exp_bound2 = %lx, "
                                      "coeff_bits = %ld, nvars = %ld\n\n",
                       len, exp_bits, exp_bound, len1, exp_bits1, exp_bound1,
                               len2, exp_bits2, exp_bound2, coeff_bits, nvars);

                fmpz_mpoly_print_pretty(f, NULL, ctx); printf("\n\n");
                fmpz_mpoly_print_pretty(g, NULL, ctx); printf("\n\n");
                fmpz_mpoly_print_pretty(h, NULL, ctx); printf("\n\n");
                fmpz_mpoly_print_pretty(k, NULL, ctx); printf("\n\n");

                flint_abort();
            }

        }

        fmpz_mpoly_clear(f, ctx);  
        fmpz_mpoly_clear(g, ctx);  
        fmpz_mpoly_clear(h, ctx);  
        fmpz_mpoly_clear(k, ctx);  
    }


    /* Check aliasing first argument */
    for (i = 0; i < 10 * flint_test_multiplier(); i++)
    {
        fmpz_mpoly_ctx_t ctx;
        fmpz_mpoly_t f, g, h;
        ordering_t ord;
        slong nvars, len, len1, len2, exp_bound, exp_bound1, exp_bound2;
        slong coeff_bits, exp_bits, exp_bits1, exp_bits2;

        ord = mpoly_ordering_randtest(state);
        nvars = n_randint(state, 10) + 1;

        fmpz_mpoly_ctx_init(ctx, nvars, ord);

        fmpz_mpoly_init(f, ctx);
        fmpz_mpoly_init(g, ctx);
        fmpz_mpoly_init(h, ctx);

        len = n_randint(state, 100);
        len1 = n_randint(state, 100);
        len2 = n_randint(state, 100);

        exp_bits = n_randint(state, FLINT_BITS - 1 - 
                  mpoly_ordering_isdeg(ord)*FLINT_BIT_COUNT(nvars)) + 1;
        exp_bits1 = n_randint(state, FLINT_BITS - 2 -
                  mpoly_ordering_isdeg(ord)*FLINT_BIT_COUNT(nvars)) + 1;
        exp_bits2 = n_randint(state, FLINT_BITS - 2 -
                  mpoly_ordering_isdeg(ord)*FLINT_BIT_COUNT(nvars)) + 1;

        exp_bound = n_randbits(state, exp_bits);
        exp_bound1 = n_randbits(state, exp_bits1);
        exp_bound2 = n_randbits(state, exp_bits2);

        coeff_bits = n_randint(state, 200);

        for (j = 0; j < 4; j++)
        {
            fmpz_mpoly_randtest(f, state, len1, exp_bound1, coeff_bits, ctx);
            fmpz_mpoly_randtest(g, state, len2, exp_bound2, coeff_bits, ctx);
            fmpz_mpoly_randtest(h, state, len, exp_bound, coeff_bits, ctx);

            flint_set_num_threads(n_randint(state, max_threads) + 1);

            fmpz_mpoly_mul_johnson(h, f, g, ctx);
            fmpz_mpoly_test(h, ctx);
  
            fmpz_mpoly_mul_heap_threaded(f, f, g, ctx);
            fmpz_mpoly_test(f, ctx);

            result = fmpz_mpoly_equal(h, f, ctx);

            if (!result)
            {
                printf("FAIL\n");
                printf("Aliasing test1\n");

                printf("ord = "); mpoly_ordering_print(ord);
                printf(", len = %ld, exp_bits = %ld, exp_bound = %lx, "
                    "len1 = %ld, exp_bits1 = %ld, exp_bound1 = %lx, "
                    "len2 = %ld, exp_bits2 = %ld, exp_bound2 = %lx, "
                                      "coeff_bits = %ld, nvars = %ld\n\n",
                       len, exp_bits, exp_bound, len1, exp_bits1, exp_bound1,
                               len2, exp_bits2, exp_bound2, coeff_bits, nvars);

                fmpz_mpoly_print_pretty(f, NULL, ctx); printf("\n\n");
                fmpz_mpoly_print_pretty(g, NULL, ctx); printf("\n\n");
                fmpz_mpoly_print_pretty(h, NULL, ctx); printf("\n\n");

                flint_abort();
            }
        }

        fmpz_mpoly_clear(f, ctx);  
        fmpz_mpoly_clear(g, ctx);  
        fmpz_mpoly_clear(h, ctx);  
    }

    /* Check aliasing second argument */
    for (i = 0; i < 10 * flint_test_multiplier(); i++)
    {
        fmpz_mpoly_ctx_t ctx;
        fmpz_mpoly_t f, g, h;
        ordering_t ord;
        slong nvars, len, len1, len2, exp_bound, exp_bound1, exp_bound2;
        slong coeff_bits, exp_bits, exp_bits1, exp_bits2;

        ord = mpoly_ordering_randtest(state);
        nvars = n_randint(state, 10) + 1;

        fmpz_mpoly_ctx_init(ctx, nvars, ord);

        fmpz_mpoly_init(f, ctx);
        fmpz_mpoly_init(g, ctx);
        fmpz_mpoly_init(h, ctx);

        len = n_randint(state, 100);
        len1 = n_randint(state, 100);
        len2 = n_randint(state, 100);

        exp_bits = n_randint(state, FLINT_BITS - 1 - 
                  mpoly_ordering_isdeg(ord)*FLINT_BIT_COUNT(nvars)) + 1;
        exp_bits1 = n_randint(state, FLINT_BITS - 2 -
                  mpoly_ordering_isdeg(ord)*FLINT_BIT_COUNT(nvars)) + 1;
        exp_bits2 = n_randint(state, FLINT_BITS - 2 -
                  mpoly_ordering_isdeg(ord)*FLINT_BIT_COUNT(nvars)) + 1;

        exp_bound = n_randbits(state, exp_bits);
        exp_bound1 = n_randbits(state, exp_bits1);
        exp_bound2 = n_randbits(state, exp_bits2);

        coeff_bits = n_randint(state, 200);

        for (j = 0; j < 4; j++)
        {
            fmpz_mpoly_randtest(f, state, len1, exp_bound1, coeff_bits, ctx);
            fmpz_mpoly_randtest(g, state, len2, exp_bound2, coeff_bits, ctx);
            fmpz_mpoly_randtest(h, state, len, exp_bound, coeff_bits, ctx);

            flint_set_num_threads(n_randint(state, max_threads) + 1);

            fmpz_mpoly_mul_johnson(h, f, g, ctx);
            fmpz_mpoly_test(h, ctx);
  
            fmpz_mpoly_mul_heap_threaded(g, f, g, ctx);
            fmpz_mpoly_test(g, ctx);

            result = fmpz_mpoly_equal(h, g, ctx);

            if (!result)
            {
                printf("FAIL\n");
                printf("Aliasing test1\n");

                printf("ord = "); mpoly_ordering_print(ord);
                printf(", len = %ld, exp_bits = %ld, exp_bound = %lx, "
                    "len1 = %ld, exp_bits1 = %ld, exp_bound1 = %lx, "
                    "len2 = %ld, exp_bits2 = %ld, exp_bound2 = %lx, "
                                      "coeff_bits = %ld, nvars = %ld\n\n",
                       len, exp_bits, exp_bound, len1, exp_bits1, exp_bound1,
                               len2, exp_bits2, exp_bound2, coeff_bits, nvars);

                fmpz_mpoly_print_pretty(f, NULL, ctx); printf("\n\n");
                fmpz_mpoly_print_pretty(g, NULL, ctx); printf("\n\n");
                fmpz_mpoly_print_pretty(h, NULL, ctx); printf("\n\n");

                flint_abort();
            }
        }

        fmpz_mpoly_clear(f, ctx);  
        fmpz_mpoly_clear(g, ctx);  
        fmpz_mpoly_clear(h, ctx);  
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}

