/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "nmod_mpoly.h"

int
main(void)
{
    slong i, j;
    FLINT_TEST_INIT(state);

    flint_printf("gcd_zippel....\n");
    fflush(stdout);

    {
        nmod_mpoly_ctx_t ctx;
        nmod_mpoly_t a, b, g;
        const char* vars[] = {"X","x1","x0"};

        nmod_mpoly_ctx_init(ctx, 3, ORD_DEGREVLEX, 101);
        nmod_mpoly_init(a, ctx);
        nmod_mpoly_init(b, ctx);
        nmod_mpoly_init(g, ctx);
        nmod_mpoly_set_str_pretty(a, "(x0+x1)*((x0^2+x0+x1)*X^2 + (x1+2)*x0^2 + x1^5+15)*(x0^3 + X + 2*x1 + 3*x0)", vars, ctx);
        nmod_mpoly_set_str_pretty(b, "(x0+x1)*((x0^2+x0+x1)*X^2 + (x1+2)*x0^2 + x1^5+15)*(x0^2 + X + 4*x1 + 5*x0)", vars, ctx);

        nmod_mpoly_gcd_zippel(g, a, b, ctx);
        nmod_mpoly_assert_canonical(g, ctx);

printf("g: "); nmod_mpoly_print_pretty(g, vars, ctx); printf("\n");

        nmod_mpoly_gcd_brown(g, a, b, ctx);
        nmod_mpoly_assert_canonical(g, ctx);

printf("g: "); nmod_mpoly_print_pretty(g, vars, ctx); printf("\n");


        nmod_mpoly_clear(g, ctx);
        nmod_mpoly_clear(a, ctx);
        nmod_mpoly_clear(b, ctx);
        nmod_mpoly_ctx_clear(ctx);
    }

    for (i = 0; i < 0*10 * flint_test_multiplier(); i++)
    {
        nmod_mpoly_ctx_t ctx;
        nmod_mpoly_t a, b, g, ca, cb, cg, t;
        slong len, len1, len2;
        slong degbound, degbound1, degbound2;
        mp_limb_t modulus;
        int res;

        modulus = n_randint(state, (i % 10 == 0) ? 4: FLINT_BITS - 1) + 1;
        modulus = n_randbits(state, modulus);
        modulus = n_nextprime(modulus, 1);

        nmod_mpoly_ctx_init_rand(ctx, state, WORD(2), modulus);
        if (ctx->minfo->nvars < 2)
            goto continue_loop;

flint_printf("modulus: %wu\n", modulus);

        nmod_mpoly_init(g, ctx);
        nmod_mpoly_init(a, ctx);
        nmod_mpoly_init(b, ctx);
        nmod_mpoly_init(ca, ctx);
        nmod_mpoly_init(cb, ctx);
        nmod_mpoly_init(cg, ctx);
        nmod_mpoly_init(t, ctx);

        len = n_randint(state, 100) + 1;
        len1 = n_randint(state, 200) + 1;
        len2 = n_randint(state, 200) + 1;

        degbound  = (10 + n_randint(state, 10))/ctx->minfo->nvars;
        degbound1 = (10 + n_randint(state, 10))/ctx->minfo->nvars;
        degbound2 = (10 + n_randint(state, 10))/ctx->minfo->nvars;

        for (j = 0; j < 1; j++)
        {
            do {
                nmod_mpoly_randtest_bound(t, state, len, degbound, ctx);
            } while (t->length == 0);
            do {
                nmod_mpoly_randtest_bound(a, state, len1, degbound1, ctx);
            } while (a->length == 0);
            do {
                nmod_mpoly_randtest_bound(b, state, len2, degbound2, ctx);
            } while (b->length == 0);
            nmod_mpoly_mul_johnson(a, a, t, ctx);
            nmod_mpoly_mul_johnson(b, b, t, ctx);

            nmod_mpoly_randtest_bits(g, state, len, FLINT_BITS, ctx);

            res = nmod_mpoly_gcd_zippel(g, a, b, ctx);

            if (!res)
            {
printf("!!!!!!!!!!!!!could not compute!!!!!!!!!!!!!!!!\n");
                continue;
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

            res = nmod_mpoly_gcd_zippel(cg, ca, cb, ctx);
            if (!res)
            {
printf("!!!!!!!!!!!!!could not compute!!!!!!!!!!!!!!!!\n");

                continue;
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
continue_loop:
        nmod_mpoly_ctx_clear(ctx);
    }

    printf("PASS\n");
    FLINT_TEST_CLEANUP(state);

    return 0;
}

