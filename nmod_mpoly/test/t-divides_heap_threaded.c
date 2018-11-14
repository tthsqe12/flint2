/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include <stdlib.h>
#include "nmod_mpoly.h"
#include "fmpz_mpoly.h"
#include "profiler.h"

int
main(void)
{
    slong test_count = 10;
    int result, result2, i;
    FLINT_TEST_INIT(state);

    flint_printf("divides_heap_threaded....\n");
    fflush(stdout);

    {
        timeit_t time;
        fmpz_mpoly_ctx_t ctx;
        fmpz_mpoly_t p, f, g, h;
        const char * vars[] = {"x","y","z","t","u"};

        fmpz_mpoly_ctx_init(ctx, 5, ORD_DEGREVLEX);
        fmpz_mpoly_init(f, ctx);
        fmpz_mpoly_init(g, ctx);
        fmpz_mpoly_init(p, ctx);
        fmpz_mpoly_init(h, ctx);
        fmpz_mpoly_set_str_pretty(f, "(1+x+y+2*z^2+3*t^3+5*u^5)^1", vars, ctx);
        fmpz_mpoly_set_str_pretty(g, "(1+u+t+2*z^2+3*y^3+5*x^5)^1", vars, ctx);
timeit_start(time);
        fmpz_mpoly_mul_johnson(p, f, g, ctx);
timeit_stop(time);
flint_printf("z mul time: %wd\n", time->wall);


timeit_start(time);
        result = fmpz_mpoly_divides_monagan_pearce(h, p, f, ctx);
timeit_stop(time);
flint_printf("z div time: %wd\n", time->wall);

        if (result == 0 || !fmpz_mpoly_equal(h, g, ctx))
        {
            printf("FAIL 0\n");
            flint_abort();
        }

        fmpz_mpoly_clear(f, ctx);
        fmpz_mpoly_clear(g, ctx);
        fmpz_mpoly_clear(p, ctx);
        fmpz_mpoly_clear(h, ctx);
        fmpz_mpoly_ctx_clear(ctx);
    }

    {
        timeit_t time;
        nmod_mpoly_ctx_t ctx;
        nmod_mpoly_t p, f, g, h, h2;
        const char * vars[] = {"x","y","z","t","u"};

        nmod_mpoly_ctx_init(ctx, 5, ORD_DEGLEX, 179424691);
        nmod_mpoly_init(f, ctx);
        nmod_mpoly_init(g, ctx);
        nmod_mpoly_init(p, ctx);
        nmod_mpoly_init(h, ctx);
        nmod_mpoly_init(h2, ctx);
        nmod_mpoly_set_str_pretty(g, "(1+x+y+2*z^2+3*t^3+5*u^5)^12", vars, ctx);
        nmod_mpoly_set_str_pretty(f, "(1+u+t+2*z^2+3*y^3+5*x^5)^12", vars, ctx);
timeit_start(time);
        nmod_mpoly_mul(p, f, g, ctx);
timeit_stop(time);
flint_printf(" mul time: %wd\n", time->wall);


timeit_start(time);
        result = nmod_mpoly_divides_monagan_pearce(h2, p, f, ctx);
timeit_stop(time);
flint_printf(" div time: %wd\n", time->wall);

/*
flint_printf("h: "); nmod_mpoly_print_pretty(h, vars, ctx); flint_printf("\n");
*/

timeit_start(time);
        result = nmod_mpoly_divides_heap_threaded(h, p, f, ctx);
timeit_stop(time);
flint_printf("pdiv time: %wd\n", time->wall);

flint_printf("result: %wd\n", result);
/*
flint_printf("h: "); nmod_mpoly_print_pretty(h, vars, ctx); flint_printf("\n");
*/
        if (result == 0 || !nmod_mpoly_equal(h, g, ctx))
        {
            printf("FAIL 1\n");
            flint_abort();
        }

printf("*******************************************\n");


        nmod_mpoly_add_ui(p, p, 1, ctx);
        result = nmod_mpoly_divides_heap_threaded(h, p, f, ctx);
flint_printf("h: "); nmod_mpoly_print_pretty(h, vars, ctx); flint_printf("\n");

        if (result != 0)
        {
            printf("FAIL 2\n");
            flint_abort();
        }

        for (i = 0; i < test_count*flint_test_multiplier(); i++)
        {
            slong len;
/*
flint_printf("first i = %wd\n", i);
*/
            len = n_randint(state, 40);
            nmod_mpoly_randtest_bound(f, state, len, 10, ctx);
            nmod_mpoly_randtest_bound(g, state, len, 10, ctx);
            nmod_mpoly_randtest_bound(p, state, len, 10, ctx);
            nmod_mpoly_randtest_bound(h, state, len, 10, ctx);

            if (nmod_mpoly_is_zero(g, ctx))
                continue;

            result = nmod_mpoly_divides_heap_threaded(h, f, g, ctx);
            result2 = nmod_mpoly_divides_monagan_pearce(p, f, g, ctx);

            if (result != result2)
            {
                flint_printf("FAIL 3 i = %wd\n", i);
                flint_abort();
            }

            if (result == 1 && !nmod_mpoly_equal(h, p, ctx))
            {
                flint_printf("FAIL 4 i = %wd\n", i);
                flint_abort();
            }
        }

        for (i = 0; i < test_count*flint_test_multiplier(); i++)
        {
            slong len;
/*
flint_printf("second i = %wd\n", i);
*/
            len = n_randint(state, 50);
            nmod_mpoly_randtest_bound(f, state, len, 20, ctx);
            nmod_mpoly_randtest_bound(g, state, len, 20, ctx);
            nmod_mpoly_randtest_bound(p, state, len, 20, ctx);
            nmod_mpoly_randtest_bound(h, state, len, 20, ctx);

            if (nmod_mpoly_is_zero(g, ctx))
                continue;

            nmod_mpoly_mul(f, f, g, ctx);
/*
flint_printf("calling threaded\n");
*/
            result = nmod_mpoly_divides_heap_threaded(h, f, g, ctx);
/*
flint_printf("returned from from threaded\n");
*/
            result2 = nmod_mpoly_divides_monagan_pearce(p, f, g, ctx);
/*
flint_printf("returned from single thread\n");
printf("h: "); nmod_mpoly_print_pretty(h, vars, ctx); printf("\n");
printf("p: "); nmod_mpoly_print_pretty(p, vars, ctx); printf("\n");
*/
            if (result != result2)
            {
                flint_printf("FAIL 3 i = %wd\n", i);
                flint_abort();
            }

            if (result == 1 && !nmod_mpoly_equal(h, p, ctx))
            {
                flint_printf("FAIL 4 i = %wd\n", i);
                flint_abort();
            }
        }

        for (i = 0; i < test_count*flint_test_multiplier(); i++)
        {
            slong len;
/*
flint_printf("third i = %wd\n", i);
*/
            len = n_randint(state, 30);
            nmod_mpoly_randtest_bound(f, state, len, 10, ctx);
            nmod_mpoly_randtest_bound(g, state, len, 10, ctx);
            nmod_mpoly_randtest_bound(p, state, len/2, 5, ctx);
            nmod_mpoly_randtest_bound(h, state, len/2, 5, ctx);

            if (nmod_mpoly_is_zero(g, ctx))
                continue;

            nmod_mpoly_mul(f, f, g, ctx);
            nmod_mpoly_add(f, f, p, ctx);
/*
flint_printf("calling threaded\n");
*/
            result = nmod_mpoly_divides_heap_threaded(h, f, g, ctx);
/*
flint_printf("returned from from threaded\n");
*/
            result2 = nmod_mpoly_divides_monagan_pearce(p, f, g, ctx);
/*
flint_printf("returned from single thread\n");
printf("h: "); nmod_mpoly_print_pretty(h, vars, ctx); printf("\n");
printf("p: "); nmod_mpoly_print_pretty(p, vars, ctx); printf("\n");
*/

            if (result != result2)
            {
                flint_printf("FAIL 3 i = %wd\n", i);
                flint_abort();
            }

            if (result == 1 && !nmod_mpoly_equal(h, p, ctx))
            {
                flint_printf("FAIL 4 i = %wd\n", i);
                flint_abort();
            }
        }

        nmod_mpoly_clear(f, ctx);
        nmod_mpoly_clear(g, ctx);
        nmod_mpoly_clear(p, ctx);
        nmod_mpoly_clear(h, ctx);
        nmod_mpoly_clear(h2, ctx);

        nmod_mpoly_ctx_clear(ctx);
    }

    FLINT_TEST_CLEANUP(state);

    printf("PASS\n");
    return 0;
}

