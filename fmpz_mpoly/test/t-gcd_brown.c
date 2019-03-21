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
#include "fmpz.h"
#include "fmpz_mpoly.h"
#include "ulong_extras.h"
#include "profiler.h"

int
main(void)
{
    int i, j;
    FLINT_TEST_INIT(state);

    flint_printf("gcd_brown....");
    fflush(stdout);

    flint_set_num_threads(2);
    if (1) {
        fmpz_mpoly_ctx_t ctx;
        fmpz_mpoly_t g, a, b;
        const char * vars[] = {"x","y","z","t","u"};
        timeit_t time;

        fmpz_mpoly_ctx_init(ctx, 3, ORD_LEX);
        fmpz_mpoly_init(g, ctx);
        fmpz_mpoly_init(a, ctx);
        fmpz_mpoly_init(b, ctx);
flint_printf("{%wd\n", 60);
timeit_start(time);
        fmpz_mpoly_set_str_pretty(a, "3*(1-17*x^6+5555*y^5+64*z^3)", vars, ctx);
        fmpz_mpoly_set_str_pretty(b, "4*(1+19*x^2+66666*y^9+81*z^2)", vars, ctx);
        fmpz_mpoly_set_str_pretty(g, "2*(7777777-3*5*7*x^9+y^3+z^9)", vars, ctx);
timeit_stop(time);
flint_printf(", (*pow time*) %wd\n", time->wall);

timeit_start(time);
        fmpz_mpoly_mul(a,a,g,ctx);
        fmpz_mpoly_mul(b,b,g,ctx);
timeit_stop(time);
flint_printf(", (*mul time*) %wd\n", time->wall);
        fmpz_mpoly_set_str_pretty(g, "1+x^200", vars, ctx);
/*
printf("a: "); nmod_mpoly_print_pretty(a, vars, ctx); printf("\n");
printf("b: "); nmod_mpoly_print_pretty(b, vars, ctx); printf("\n");
*/
timeit_start(time);
        fmpz_mpoly_gcd_brownnew(g, a, b, ctx);
timeit_stop(time);
flint_printf(", (* tot gcd time *) %wd\n", time->wall);


printf("g: "); fmpz_mpoly_print_pretty(g, vars, ctx); printf("\n");


timeit_start(time);
        fmpz_mpoly_divides(a, a, g, ctx);
        fmpz_mpoly_divides(b, b, g, ctx);
timeit_stop(time);
flint_printf(", (* div time *) %wd}\n", time->wall);

/*
printf("a: "); nmod_mpoly_print_pretty(a, vars, ctx); printf("\n");
printf("b: "); nmod_mpoly_print_pretty(b, vars, ctx); printf("\n");

        fmpz_mpoly_gcd_brownnew(g, a, b, ctx);

printf("g: "); fmpz_mpoly_print_pretty(g, vars, ctx); printf("\n");
*/

        fmpz_mpoly_clear(g, ctx);
        fmpz_mpoly_clear(a, ctx);
        fmpz_mpoly_clear(b, ctx);
        fmpz_mpoly_ctx_clear(ctx);
    }



    for (i = 0; i < 0 * flint_test_multiplier(); i++)
    {
        fmpz_mpoly_ctx_t ctx;
        fmpz_mpoly_t a, b, g, ca, cb, cg, t;
        mp_bitcnt_t coeff_bits;
        slong len, len1, len2;
        slong degbound;
        int res;

        fmpz_mpoly_ctx_init_rand(ctx, state, 4);

        fmpz_mpoly_init(g, ctx);
        fmpz_mpoly_init(a, ctx);
        fmpz_mpoly_init(b, ctx);
        fmpz_mpoly_init(ca, ctx);
        fmpz_mpoly_init(cb, ctx);
        fmpz_mpoly_init(cg, ctx);
        fmpz_mpoly_init(t, ctx);

        len = n_randint(state, 25) + 1;
        len1 = n_randint(state, 50);
        len2 = n_randint(state, 50);

        degbound = 1 + 25/(2*ctx->minfo->nvars - 1);

        coeff_bits = n_randint(state, 400);

        for (j = 0; j < 4; j++)
        {

flint_printf("i = %wd, j = %wd, nvars = %wd\n", i, j, ctx->minfo->nvars);

            do {
                fmpz_mpoly_randtest_bound(t, state, len, coeff_bits + 1, degbound, ctx);
            } while (t->length == 0);
            fmpz_mpoly_randtest_bound(a, state, len1, coeff_bits, degbound, ctx);
            fmpz_mpoly_randtest_bound(b, state, len2, coeff_bits, degbound, ctx);
            fmpz_mpoly_mul(a, a, t, ctx);
            fmpz_mpoly_mul(b, b, t, ctx);

            fmpz_mpoly_randtest_bits(g, state, len, coeff_bits, FLINT_BITS, ctx);

            res = fmpz_mpoly_gcd_brownnew(g, a, b, ctx);
            fmpz_mpoly_assert_canonical(g, ctx);

            if (!res) {
                continue;
            }

            if (fmpz_mpoly_is_zero(g, ctx))
            {
                if (!fmpz_mpoly_is_zero(a, ctx) || !fmpz_mpoly_is_zero(b, ctx)) {
                    printf("FAIL\n");
                    flint_printf("Check zero gcd only results from zero inputs\ni = %wd, j = %wd\n", i ,j);
                    flint_abort();
                }
                continue;
            }

            if (fmpz_sgn(g->coeffs + 0) <= 0)
            {
                printf("FAIL\n");
                flint_printf("Check gcd has positive lc\ni = %wd, j = %wd\n", i ,j);
                flint_abort();
            }

            res = 1;
            res = res && fmpz_mpoly_divides(ca, a, g, ctx);
            res = res && fmpz_mpoly_divides(cb, b, g, ctx);
            if (!res)
            {
                printf("FAIL\n");
                flint_printf("Check divisibility\ni = %wd, j = %wd\n", i ,j);
                flint_abort();
            }

            res = fmpz_mpoly_gcd_brownnew(cg, ca, cb, ctx);

            if (!res)
                continue;

            if (!fmpz_mpoly_equal_ui(cg, UWORD(1), ctx))
            {
                printf("FAIL\n");
                flint_printf("Check cofactors are relatively prime\ni = %wd, j = %wd\n", i ,j);                
                flint_abort();
            }
        }

        fmpz_mpoly_clear(g, ctx);
        fmpz_mpoly_clear(a, ctx);
        fmpz_mpoly_clear(b, ctx);
        fmpz_mpoly_clear(ca, ctx);
        fmpz_mpoly_clear(cb, ctx);
        fmpz_mpoly_clear(cg, ctx);
        fmpz_mpoly_clear(t, ctx);
        fmpz_mpoly_ctx_clear(ctx);
    }


    FLINT_TEST_CLEANUP(state);

    printf("PASS\n");
    return 0;
}

