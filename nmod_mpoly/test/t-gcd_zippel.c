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

    if (0) {
        nmod_mpoly_ctx_t ctx;
        nmod_mpoly_t a, b, g;
        const char* vars[] = {"X","x1","x0"};

        nmod_mpoly_ctx_init(ctx, 3, ORD_DEGREVLEX, 101);
        nmod_mpoly_init(a, ctx);
        nmod_mpoly_init(b, ctx);
        nmod_mpoly_init(g, ctx);
        nmod_mpoly_set_str_pretty(a, "(x0^2+x1+1)*(x1+x0+1)*((x0^2+x0+x1)*X^2 + (x1+2)*x0^2 + x1^2+15)*(1 + X + 2*x1 + 3*x0)", vars, ctx);
        nmod_mpoly_set_str_pretty(b, "(x0^2+x1+1)*(x1+x0+1)*((x0^2+x0+x1)*X^2 + (x1+2)*x0^2 + x1^2+15)*(1 + X + 4*x1 + 5*x0)", vars, ctx);

        nmod_mpoly_gcd_zippel(g, a, b, ctx);
        nmod_mpoly_assert_canonical(g, ctx);

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
        ulong degbound;
        ulong * degbounds, * degbounds1, * degbounds2;
        mp_limb_t modulus;
        int res;

        modulus = n_randint(state, (i % 10 == 0) ? 4: FLINT_BITS - 1) + 1;
        modulus = n_randbits(state, modulus);
        modulus = n_nextprime(modulus, 1);

        nmod_mpoly_ctx_init_rand(ctx, state, WORD(10), modulus);
/*
flint_printf("trying example %wd modulus: %wu nvars: %wd\n", i,modulus,ctx->minfo->nvars);
*/
        nmod_mpoly_init(g, ctx);
        nmod_mpoly_init(a, ctx);
        nmod_mpoly_init(b, ctx);
        nmod_mpoly_init(ca, ctx);
        nmod_mpoly_init(cb, ctx);
        nmod_mpoly_init(cg, ctx);
        nmod_mpoly_init(t, ctx);

        len = n_randint(state, 15) + 1;
        len1 = n_randint(state, 15);
        len2 = n_randint(state, 15);

        degbound = 100/(2*ctx->minfo->nvars - 1);
        degbounds = (ulong * ) flint_malloc(ctx->minfo->nvars*sizeof(ulong));
        degbounds1 = (ulong * ) flint_malloc(ctx->minfo->nvars*sizeof(ulong));
        degbounds2 = (ulong * ) flint_malloc(ctx->minfo->nvars*sizeof(ulong));
        for (j = 0; j < ctx->minfo->nvars; j++)
        {
            degbounds[j] = n_randint(state, degbound + UWORD(1)) + UWORD(1);
            degbounds1[j] = n_randint(state, degbound + UWORD(1)) + UWORD(1);
            degbounds2[j] = n_randint(state, degbound + UWORD(1)) + UWORD(1);
        }

        for (j = 0; j < 4; j++)
        {
            do {
                nmod_mpoly_randtest_bounds(t, state, len, degbounds, ctx);
            } while (t->length == 0);
            nmod_mpoly_randtest_bounds(a, state, len1, degbounds1, ctx);
            nmod_mpoly_randtest_bounds(b, state, len2, degbounds2, ctx);

/*
printf(" --------------- t: "); nmod_mpoly_print_pretty(t, NULL, ctx); printf("\n");
printf(" --------------- a: "); nmod_mpoly_print_pretty(a, NULL, ctx); printf("\n");
printf(" --------------- b: "); nmod_mpoly_print_pretty(b, NULL, ctx); printf("\n");
*/

            nmod_mpoly_mul_johnson(a, a, t, ctx);
            nmod_mpoly_mul_johnson(b, b, t, ctx);




            nmod_mpoly_randtest_bits(g, state, len, FLINT_BITS, ctx);

            res = nmod_mpoly_gcd_zippel(g, a, b, ctx);
            if (!res)
            {
flint_printf("!could not compute modulus = %wu !\n", modulus);
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
flint_printf("!could not compute modulus = %wu !\n", modulus);
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
        nmod_mpoly_ctx_clear(ctx);
    }

    printf("PASS\n");
    FLINT_TEST_CLEANUP(state);

    return 0;
}

