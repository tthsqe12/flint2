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


int
main(void)
{
    int result, result2, i;
    FLINT_TEST_INIT(state);

    flint_printf("divides_heap_threaded....\n");
    fflush(stdout);

    {
        nmod_mpoly_ctx_t ctx;
        nmod_mpoly_t p, f, g, h;
        const char * vars[] = {"x","y","z","t","u"};

        nmod_mpoly_ctx_init(ctx, 5, ORD_DEGLEX, 101);
        nmod_mpoly_init(f, ctx);
        nmod_mpoly_init(g, ctx);
        nmod_mpoly_init(p, ctx);
        nmod_mpoly_init(h, ctx);
        nmod_mpoly_set_str_pretty(f, "(1+x+y+2*z^2+3*t^3+5*u^5)", vars, ctx);
        nmod_mpoly_set_str_pretty(g, "(1+u+t+2*z^2+3*y^3+5*x^5)", vars, ctx);
        nmod_mpoly_mul(p, f, g, ctx);

flint_printf("f: "); nmod_mpoly_print_pretty(f, vars, ctx); flint_printf("\n");
flint_printf("g: "); nmod_mpoly_print_pretty(g, vars, ctx); flint_printf("\n");
flint_printf("p: "); nmod_mpoly_print_pretty(p, vars, ctx); flint_printf("\n");


        result = nmod_mpoly_divides_heap_threaded(h, p, f, ctx);
flint_printf("h: "); nmod_mpoly_print_pretty(h, vars, ctx); flint_printf("\n");

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

        for (i = 0; i < 10*flint_test_multiplier(); i++)
        {
            slong len;

flint_printf("first i = %wd\n", i);

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

        for (i = 0; i < 10*flint_test_multiplier(); i++)
        {
            slong len;

flint_printf("second i = %wd\n", i);

            len = n_randint(state, 40);
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

        for (i = 0; i < 10*flint_test_multiplier(); i++)
        {
            slong len;

flint_printf("third i = %wd\n", i);

            len = n_randint(state, 20);
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

        nmod_mpoly_ctx_clear(ctx);
    }

    FLINT_TEST_CLEANUP(state);

    printf("PASS\n");
    return 0;
}

