/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "nmod_mpoly.h"
#include "profiler.h"

int
main(void)
{
    slong i, j;
    FLINT_TEST_INIT(state);

    flint_printf("gcd_brown....\n");
    fflush(stdout);

flint_set_num_threads(1);

    if (1) {
        nmod_mpoly_ctx_t ctx;
        nmod_mpoly_t g, a, b;
        const char * vars[] = {"x","y","z","t","u"};
        timeit_t time;

        nmod_mpoly_ctx_init(ctx, 3, ORD_LEX, 179424691);
        nmod_mpoly_init(g, ctx);
        nmod_mpoly_init(a, ctx);
        nmod_mpoly_init(b, ctx);
flint_printf("{%wd\n", 60);
timeit_start(time);
        nmod_mpoly_set_str_pretty(a, "(1+x^2+y^3+z)^5", vars, ctx);
        nmod_mpoly_set_str_pretty(b, "(1+x+y+z^2)^5", vars, ctx);
        nmod_mpoly_set_str_pretty(g, "(1+x+z^2+z)^5", vars, ctx);
timeit_stop(time);
flint_printf(", (*pow time*) %wd\n", time->wall);

timeit_start(time);
        nmod_mpoly_mul(a,a,g,ctx);
        nmod_mpoly_mul(b,b,g,ctx);
timeit_stop(time);
flint_printf(", (*mul time*) %wd\n", time->wall);
        nmod_mpoly_set_str_pretty(g, "1+x^200", vars, ctx);
/*
printf("a: "); nmod_mpoly_print_pretty(a, vars, ctx); printf("\n");
printf("b: "); nmod_mpoly_print_pretty(b, vars, ctx); printf("\n");
*/
timeit_start(time);
        nmod_mpoly_gcd_brownnew(g, a, b, ctx);
timeit_stop(time);
flint_printf(", (* tot gcd time *) %wd\n", time->wall);

/*
printf("first g: "); nmod_mpoly_print_pretty(g, vars, ctx); printf("\n");
*/

timeit_start(time);
        nmod_mpoly_divides(a, a, g, ctx);
        nmod_mpoly_divides(b, b, g, ctx);
timeit_stop(time);
flint_printf(", (* div time *) %wd}\n", time->wall);

/*
printf("a: "); nmod_mpoly_print_pretty(a, vars, ctx); printf("\n");
printf("b: "); nmod_mpoly_print_pretty(b, vars, ctx); printf("\n");
*/
        nmod_mpoly_gcd_brown(g, a, b, ctx);

printf("cofactor g: "); nmod_mpoly_print_pretty(g, vars, ctx); printf("\n");


        nmod_mpoly_clear(g, ctx);
        nmod_mpoly_clear(a, ctx);
        nmod_mpoly_clear(b, ctx);
        nmod_mpoly_ctx_clear(ctx);
    }

    if (0) {
        nmod_mpoly_ctx_t ctx;
        nmod_mpoly_t g, a, b;
        const char * vars[] = {"x","y","z","t","u"};
        timeit_t time;

        nmod_mpoly_ctx_init(ctx, 4, ORD_DEGLEX, 179424691);
        nmod_mpoly_init(g, ctx);
        nmod_mpoly_init(a, ctx);
        nmod_mpoly_init(b, ctx);
        nmod_mpoly_set_str_pretty(a, "(1+x+y-x*y+z+t+t*z)^20", vars, ctx);
        nmod_mpoly_set_str_pretty(b, "(1+x+y-x*y+z-t+t*z)^20", vars, ctx);
        nmod_mpoly_set_str_pretty(g, "(1+x+2*x*t+y+z+t)^20", vars, ctx);
        nmod_mpoly_mul(a,a,g,ctx);
        nmod_mpoly_mul(b,b,g,ctx);
        nmod_mpoly_set_str_pretty(g, "1+x^200", vars, ctx);
/*
printf("a: "); nmod_mpoly_print_pretty(a, vars, ctx); printf("\n");
printf("b: "); nmod_mpoly_print_pretty(b, vars, ctx); printf("\n");
*/
timeit_start(time);
        nmod_mpoly_gcd_brownnew(g, a, b, ctx);
timeit_stop(time);
flint_printf("gcd time: %wd\n", time->wall);
/*
printf("g: "); nmod_mpoly_print_pretty(g, vars, ctx); printf("\n");
*/

        nmod_mpoly_divides(a, a, g, ctx);
        nmod_mpoly_divides(b, b, g, ctx);
/*
printf("a: "); nmod_mpoly_print_pretty(a, vars, ctx); printf("\n");
printf("b: "); nmod_mpoly_print_pretty(b, vars, ctx); printf("\n");
*/
        nmod_mpoly_gcd_brownnew(g, a, b, ctx);

printf("g: "); nmod_mpoly_print_pretty(g, vars, ctx); printf("\n");


        nmod_mpoly_clear(g, ctx);
        nmod_mpoly_clear(a, ctx);
        nmod_mpoly_clear(b, ctx);
        nmod_mpoly_ctx_clear(ctx);
    }

    for (i = 0; i < 10 * flint_test_multiplier(); i++)
    {
        nmod_mpoly_ctx_t ctx;
        nmod_mpoly_t a, b, g, ca, cb, cg, t;
        slong len, len1, len2;
        slong degbound;
        mp_limb_t modulus;
        int res;

        modulus = n_randint(state, FLINT_BITS - 1) + 1;
        modulus = n_randbits(state, modulus);
        modulus = n_nextprime(modulus, 1);

        nmod_mpoly_ctx_init(ctx, 4, ORD_LEX, modulus);

        nmod_mpoly_init(g, ctx);
        nmod_mpoly_init(a, ctx);
        nmod_mpoly_init(b, ctx);
        nmod_mpoly_init(ca, ctx);
        nmod_mpoly_init(cb, ctx);
        nmod_mpoly_init(cg, ctx);
        nmod_mpoly_init(t, ctx);

        len = n_randint(state, 100) + 1;
        len1 = n_randint(state, 100);
        len2 = n_randint(state, 100);

        degbound = 1 + 60/ctx->minfo->nvars/ctx->minfo->nvars;

flint_printf("i = %wd\n", i);

        for (j = 0; j < 4; j++)
        {
            do {
                nmod_mpoly_randtest_bound(t, state, len, degbound, ctx);
            } while (t->length == 0);
            nmod_mpoly_randtest_bound(a, state, len1, degbound, ctx);
            nmod_mpoly_randtest_bound(b, state, len2, degbound, ctx);
/*
printf("a: "); nmod_mpoly_print_pretty(a,NULL,ctx); printf("\n");
printf("b: "); nmod_mpoly_print_pretty(b,NULL,ctx); printf("\n");
printf("t: "); nmod_mpoly_print_pretty(t,NULL,ctx); printf("\n");
*/
            nmod_mpoly_mul_johnson(a, a, t, ctx);
            nmod_mpoly_mul_johnson(b, b, t, ctx);

            nmod_mpoly_randtest_bits(g, state, len, FLINT_BITS, ctx);

            res = nmod_mpoly_gcd_brownnew(g, a, b, ctx);
            if (!res)
            {
                printf("FAIL\n");
                flint_printf("Check gcd could be computed\ni = %wd, j = %wd\n", i ,j);
                flint_abort();
            }
            nmod_mpoly_assert_canonical(g, ctx);

            if (nmod_mpoly_is_zero(g, ctx))
            {
                if (!nmod_mpoly_is_zero(a, ctx) || !nmod_mpoly_is_zero(b, ctx))
                {
                    printf("FAIL\n");
                    flint_printf("Check zero gcd only results from zero inputs\ni = %wd, j = %wd\n", i ,j);
                    flint_abort();
                }
                continue;
            }

            if (g->coeffs[0] != UWORD(1))
            {
                printf("FAIL\n");
                flint_printf("Check gcd is monic\ni = %wd, j = %wd\n", i ,j);
                flint_abort();
            }

            res = 1;
            res = res && nmod_mpoly_divides_monagan_pearce(ca, a, g, ctx);
            res = res && nmod_mpoly_divides_monagan_pearce(cb, b, g, ctx);
            if (!res)
            {
                printf("FAIL\n");
                flint_printf("Check divisibility\ni = %wd, j = %wd\n", i ,j);
                flint_abort();
            }

            res = nmod_mpoly_gcd_brownnew(cg, ca, cb, ctx);

            if (!res)
            {
                printf("FAIL\n");
                flint_printf("Check cofactor gcd could be computed\ni = %wd, j = %wd\n", i ,j);
                flint_abort();
            }

            if (!nmod_mpoly_equal_ui(cg, UWORD(1), ctx))
            {
                printf("FAIL\n");
                flint_printf("Check cofactors are relatively prime\ni = %wd, j = %wd\n", i ,j);                
                flint_abort();
            }
        }

        nmod_mpoly_clear(g, ctx);
        nmod_mpoly_clear(a, ctx);
        nmod_mpoly_clear(b, ctx);
        nmod_mpoly_clear(ca, ctx);
        nmod_mpoly_clear(cb, ctx);
        nmod_mpoly_clear(cg, ctx);
        nmod_mpoly_clear(t, ctx);
        nmod_mpoly_ctx_clear(ctx);
    }


    printf("PASS\n");
    FLINT_TEST_CLEANUP(state);

    return 0;
}

