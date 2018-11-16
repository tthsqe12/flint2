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
#include "nmod_mpoly.h"
#include "profiler.h"

#define multiplier 0
int
main(void)
{
    int i, j, result, max_threads = 5;
    FLINT_TEST_INIT(state);

    flint_printf("mul_heap_threaded....");
    fflush(stdout);

    {
        timeit_t time;
        nmod_mpoly_ctx_t ctx;
        nmod_mpoly_t p, f, g, h;
        const char * vars[] = {"x","y","z","t","u"};

flint_printf("\n******************\n");

        nmod_mpoly_ctx_init(ctx, 5, ORD_DEGLEX, UWORD(2147483659));
        nmod_mpoly_init(f, ctx);
        nmod_mpoly_init(g, ctx);
        nmod_mpoly_init(p, ctx);
        nmod_mpoly_init(h, ctx);
        nmod_mpoly_set_str_pretty(g, "(1+x+y+2*z^2+3*t^3+5*u^5)^14", vars, ctx);
        nmod_mpoly_set_str_pretty(f, "(1+u+t+2*z^2+3*y^3+5*x^5)^14", vars, ctx);
timeit_start(time);
        nmod_mpoly_mul(p, f, g, ctx);
timeit_stop(time);
flint_printf(" mul time: %wd\n", time->wall);

    for (i = 1; i <= 4; i++)
    {
        flint_set_num_threads(i);
timeit_start(time);
        nmod_mpoly_mul_heap_threaded(h, f, g, ctx);
timeit_stop(time);
flint_printf("pmul time(%wd): %wd\n", i, time->wall);

        if (!nmod_mpoly_equal(p, h, ctx))
        {
            printf("FAIL\n");
            flint_abort();
        }
    }

        nmod_mpoly_clear(f, ctx);
        nmod_mpoly_clear(g, ctx);
        nmod_mpoly_clear(p, ctx);
        nmod_mpoly_clear(h, ctx);
        nmod_mpoly_ctx_clear(ctx);
    }

    /* Check mul_heap_threaded matches mul_johnson */
    for (i = 0; i < multiplier * flint_test_multiplier(); i++)
    {
        nmod_mpoly_ctx_t ctx;
        nmod_mpoly_t f, g, h, k;
        mp_limb_t modulus;
        slong len, len1, len2;
        mp_bitcnt_t exp_bits, exp_bits1, exp_bits2;

        modulus = n_randint(state, FLINT_BITS - 1) + 1;
        modulus = n_randbits(state, modulus);
        modulus = n_nextprime(modulus, 1);
        nmod_mpoly_ctx_init_rand(ctx, state, 10, modulus);

        nmod_mpoly_init(f, ctx);
        nmod_mpoly_init(g, ctx);
        nmod_mpoly_init(h, ctx);
        nmod_mpoly_init(k, ctx);

        len = n_randint(state, 100);
        len1 = n_randint(state, 100);
        len2 = n_randint(state, 100);

        exp_bits = n_randint(state, 200) + 2;
        exp_bits1 = n_randint(state, 200) + 2;
        exp_bits2 = n_randint(state, 200) + 2;

        for (j = 0; j < 4; j++)
        {
            nmod_mpoly_randtest_bits(f, state, len1, exp_bits1, ctx);
            nmod_mpoly_randtest_bits(g, state, len2, exp_bits2, ctx);
            nmod_mpoly_randtest_bits(h, state, len, exp_bits, ctx);
            nmod_mpoly_randtest_bits(k, state, len, exp_bits, ctx);

            flint_set_num_threads(n_randint(state, max_threads) + 1);

            nmod_mpoly_mul_johnson(h, f, g, ctx);
            nmod_mpoly_test(h, ctx);
            nmod_mpoly_mul_heap_threaded(k, f, g, ctx);
            nmod_mpoly_test(k, ctx);
            result = nmod_mpoly_equal(h, k, ctx);

            if (!result)
            {
                printf("FAIL\n");
                flint_printf("Check mul_heap_threaded matches mul_johnson\ni = %wd, j = %wd\n", i ,j);
                flint_abort();
            }
        }

        nmod_mpoly_clear(f, ctx);
        nmod_mpoly_clear(g, ctx);
        nmod_mpoly_clear(h, ctx);
        nmod_mpoly_clear(k, ctx);

        nmod_mpoly_ctx_clear(ctx);
    }


    /* Check aliasing first argument */
    for (i = 0; i < multiplier * flint_test_multiplier(); i++)
    {
        nmod_mpoly_ctx_t ctx;
        nmod_mpoly_t f, g, h;
        mp_limb_t modulus;
        slong len, len1, len2;
        slong exp_bits, exp_bits1, exp_bits2;

        modulus = n_randint(state, FLINT_BITS - 1) + 1;
        modulus = n_randbits(state, modulus);
        modulus = n_nextprime(modulus, 1);
        nmod_mpoly_ctx_init_rand(ctx, state, 10, modulus);

        nmod_mpoly_init(f, ctx);
        nmod_mpoly_init(g, ctx);
        nmod_mpoly_init(h, ctx);

        len = n_randint(state, 100);
        len1 = n_randint(state, 100);
        len2 = n_randint(state, 100);

        exp_bits = n_randint(state, 200) + 2;
        exp_bits1 = n_randint(state, 200) + 2;
        exp_bits2 = n_randint(state, 200) + 2;

        for (j = 0; j < 4; j++)
        {
            nmod_mpoly_randtest_bits(f, state, len1, exp_bits1, ctx);
            nmod_mpoly_randtest_bits(g, state, len2, exp_bits2, ctx);
            nmod_mpoly_randtest_bits(h, state, len, exp_bits, ctx);

            flint_set_num_threads(n_randint(state, max_threads) + 1);

            nmod_mpoly_mul_johnson(h, f, g, ctx);
            nmod_mpoly_test(h, ctx);
            nmod_mpoly_mul_heap_threaded(f, f, g, ctx);
            nmod_mpoly_test(f, ctx);
            result = nmod_mpoly_equal(h, f, ctx);

            if (!result)
            {
                printf("FAIL\n");
                flint_printf("Check aliasing first argument\ni = %wd, j = %wd\n", i ,j);
                flint_abort();
            }
        }

        nmod_mpoly_clear(f, ctx);
        nmod_mpoly_clear(g, ctx);
        nmod_mpoly_clear(h, ctx);

        nmod_mpoly_ctx_clear(ctx);
    }

    /* Check aliasing second argument */
    for (i = 0; i < multiplier * flint_test_multiplier(); i++)
    {
        nmod_mpoly_ctx_t ctx;
        nmod_mpoly_t f, g, h;
        mp_limb_t modulus;
        slong len, len1, len2;
        slong exp_bits, exp_bits1, exp_bits2;

        modulus = n_randint(state, FLINT_BITS - 1) + 1;
        modulus = n_randbits(state, modulus);
        modulus = n_nextprime(modulus, 1);
        nmod_mpoly_ctx_init_rand(ctx, state, 10, modulus);

        nmod_mpoly_init(f, ctx);
        nmod_mpoly_init(g, ctx);
        nmod_mpoly_init(h, ctx);

        len = n_randint(state, 100);
        len1 = n_randint(state, 100);
        len2 = n_randint(state, 100);

        exp_bits = n_randint(state, 200) + 2;
        exp_bits1 = n_randint(state, 200) + 2;
        exp_bits2 = n_randint(state, 200) + 2;

        for (j = 0; j < 4; j++)
        {
            nmod_mpoly_randtest_bits(f, state, len1, exp_bits1, ctx);
            nmod_mpoly_randtest_bits(g, state, len2, exp_bits2, ctx);
            nmod_mpoly_randtest_bits(h, state, len, exp_bits, ctx);

            flint_set_num_threads(n_randint(state, max_threads) + 1);

            nmod_mpoly_mul_johnson(h, f, g, ctx);
            nmod_mpoly_test(h, ctx);
            nmod_mpoly_mul_heap_threaded(g, f, g, ctx);
            nmod_mpoly_test(g, ctx);
            result = nmod_mpoly_equal(h, g, ctx);

            if (!result)
            {
                printf("FAIL\n");
                flint_printf("Check aliasing first argument\ni = %wd, j = %wd\n", i ,j);
                flint_abort();
            }
        }

        nmod_mpoly_clear(f, ctx);
        nmod_mpoly_clear(g, ctx);
        nmod_mpoly_clear(h, ctx);

        nmod_mpoly_ctx_clear(ctx);
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}

