/*
    Copyright (C) 2019 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fq_nmod_mpoly.h"

int
main(void)
{
    slong i;
    FLINT_TEST_INIT(state);

    flint_printf("gcd_brown....\n");
    fflush(stdout);

    /* check simple example */
    {
        fq_nmod_mpoly_ctx_t ctx;
        fq_nmod_mpoly_t g, a, b;
        /* get the coefficient of y^1*x^2*/
        const char * vars[] = {"x", "y", "z", "w"};

        fq_nmod_mpoly_ctx_init_deg(ctx, 3, ORD_LEX, 2, 2);
        fq_nmod_mpoly_init(g, ctx);
        fq_nmod_mpoly_init(a, ctx);
        fq_nmod_mpoly_init(b, ctx);
        fq_nmod_mpoly_set_str_pretty(a, "(x+y+z^2)*(x-y^9+z^3)", vars, ctx);
        fq_nmod_mpoly_set_str_pretty(b, "(x+y+z^9)*(x^9+y+z^2)", vars, ctx);

printf("a: "); fq_nmod_mpoly_print_pretty(a, vars, ctx); printf("\n");
printf("b: "); fq_nmod_mpoly_print_pretty(b, vars, ctx); printf("\n");

        fq_nmod_mpoly_gcd_brownnew(g, a, b, ctx);

printf("g: "); fq_nmod_mpoly_print_pretty(g, vars, ctx); printf("\n");

        fq_nmod_mpoly_clear(g, ctx);
        fq_nmod_mpoly_clear(a, ctx);
        fq_nmod_mpoly_clear(b, ctx);
        fq_nmod_mpoly_ctx_clear(ctx);
    }

    for (i = 0; i < 10*flint_test_multiplier(); i++)
    {
        fq_nmod_mpoly_ctx_t ctx;
        fq_nmod_mpoly_t a, b, g, ca, cb, cg, t;
        slong len, len1, len2;
        slong degbound;
        int res;
        mp_bitcnt_t pbits;
        slong deg;

        pbits = 1 + n_randint(state, FLINT_BITS);
        pbits = 1 + n_randint(state, pbits);
        deg = 1 + n_randint(state, 4);
        fq_nmod_mpoly_ctx_init_rand(ctx, state, 4, pbits, deg);

        fq_nmod_mpoly_init(g, ctx);
        fq_nmod_mpoly_init(a, ctx);
        fq_nmod_mpoly_init(b, ctx);
        fq_nmod_mpoly_init(ca, ctx);
        fq_nmod_mpoly_init(cb, ctx);
        fq_nmod_mpoly_init(cg, ctx);
        fq_nmod_mpoly_init(t, ctx);

        len = n_randint(state, 100) + 1;
        len1 = n_randint(state, 150);
        len2 = n_randint(state, 150);

flint_printf("i: %wd\n", i);

        degbound = 1 + 60/ctx->minfo->nvars/ctx->minfo->nvars;

        do {
            fq_nmod_mpoly_randtest_bound(t, state, len, degbound, ctx);
        } while (t->length == 0);
        fq_nmod_mpoly_randtest_bound(a, state, len1, degbound, ctx);
        fq_nmod_mpoly_randtest_bound(b, state, len2, degbound, ctx);
        fq_nmod_mpoly_mul(a, a, t, ctx);
        fq_nmod_mpoly_mul(b, b, t, ctx);

        fq_nmod_mpoly_randtest_bits(g, state, len, FLINT_BITS, ctx);

        res = fq_nmod_mpoly_gcd_brownnew(g, a, b, ctx);
        if (!res)
        {
            printf("FAIL\n");
            flint_printf("Could not compute gcd\n"
                            "i = %wd, p = %wu, d = %wd\n", i, ctx->fqctx->modulus->mod.n, ctx->fqctx->modulus->length - 1);
            flint_abort();
        }
        fq_nmod_mpoly_assert_canonical(g, ctx);

        if (fq_nmod_mpoly_is_zero(g, ctx))
        {
            if (!fq_nmod_mpoly_is_zero(a, ctx) || !fq_nmod_mpoly_is_zero(b, ctx))
            {
                printf("FAIL\n");
                flint_printf("Check zero gcd only results from zero inputs\ni = %wd\n", i);
                flint_abort();
            }
            goto break_continue;
        }

        if (!fq_nmod_is_one(g->coeffs + 0, ctx->fqctx))
        {
            printf("FAIL\n");
            flint_printf("Check gcd is monic\ni = %wd\n", i);
            flint_abort();
        }

        res = 1;
        res = res && fq_nmod_mpoly_divides(ca, a, g, ctx);
        res = res && fq_nmod_mpoly_divides(cb, b, g, ctx);
        if (!res)
        {
            printf("FAIL\n");
            flint_printf("Check divisibility\ni = %wd\n", i);
            flint_abort();
        }

        res = fq_nmod_mpoly_gcd_brownnew(cg, ca, cb, ctx);
        if (!res)
        {
            printf("FAIL\n");
            flint_printf("Could not compute cofactor gcd\ni = %wd\n", i);
            flint_abort();
        }

        if (!fq_nmod_mpoly_is_one(cg, ctx))
        {
            printf("FAIL\n");
            flint_printf("Check cofactors are relatively prime\ni = %wd\n", i);
            flint_abort();
        }

break_continue:

        fq_nmod_mpoly_clear(g, ctx);
        fq_nmod_mpoly_clear(a, ctx);
        fq_nmod_mpoly_clear(b, ctx);
        fq_nmod_mpoly_clear(ca, ctx);
        fq_nmod_mpoly_clear(cb, ctx);
        fq_nmod_mpoly_clear(cg, ctx);
        fq_nmod_mpoly_clear(t, ctx);
        fq_nmod_mpoly_ctx_clear(ctx);
    }


    printf("PASS\n");
    FLINT_TEST_CLEANUP(state);

    return 0;
}
