/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fq_poly.h"


void fq_poly_rand_monic(fq_poly_t f, flint_rand_t state, slong len, fq_ctx_t ctx)
{
    slong i;

    fq_poly_fit_length(f, len, ctx);
    f->length = len;
    for (i = 0; i < len; i++)
    {
        if (i + 1 < len)
            fq_rand(f->coeffs + i, state, ctx);
        else
            fq_one(f->coeffs + i, ctx);
    }
    f->length = len;
    _fq_poly_normalise(f, ctx);
}


int
main(void)
{
    slong i, j;
    FLINT_TEST_INIT(state);

    flint_printf("xgcd_euclidean_f_composite_characteristic....");
    fflush(stdout);

    for (i = 0; i < 10 * flint_test_multiplier(); i++)
    {
        fmpz_t p;
        fmpz_mod_poly_t h;
        fq_ctx_t ctx;
        fq_poly_t a, b, g, s, t, t1, t2;
        fq_t f;

        fmpz_init(p);
        fmpz_randtest_unsigned(p, state, 200);
        fmpz_add_ui(p, p, 2 + n_randint(state, 10));
        fmpz_nextprime(p, p, 1);

        fmpz_mod_poly_init(h, p);
        fmpz_mod_poly_randtest_monic(h, state, 2 + n_randint(state, 7));

flint_printf("h: "); fmpz_mod_poly_print_pretty(h, "t"); flint_printf(" mod "); fmpz_print(p); flint_printf("\n");

        fq_ctx_init_modulus(ctx, h, "t");
        fq_poly_init(a, ctx);
        fq_poly_init(b, ctx);
        fq_poly_init(g, ctx);
        fq_poly_init(s, ctx);
        fq_poly_init(t, ctx);
        fq_poly_init(t1, ctx);
        fq_poly_init(t2, ctx);
        fq_init(f, ctx);

        for (j = 0; j < 10 * flint_test_multiplier(); j++)
        {
flint_printf("i = %wd, j = %wd\n", i, j);

            fq_poly_rand_monic(a, state, n_randint(state, 10), ctx);
            fq_poly_rand_monic(b, state, n_randint(state, 10), ctx);

            fq_poly_xgcd_euclidean_f(f, g, s, t, a, b, ctx);
            if (fq_is_one(f, ctx))
            {
                fq_poly_mul(t1, s, a, ctx);
                fq_poly_mul(t2, t, b, ctx);
                fq_poly_add(t1, t1, t2, ctx);
                if (!fq_poly_equal(t1, g, ctx))
                {
                    flint_printf("FAIL\n");
                    flint_printf("f: "); fq_print_pretty(f, ctx); flint_printf("\n");
                    flint_printf("a: "); fq_poly_print_pretty(a, "x", ctx); flint_printf("\n");
                    flint_printf("b: "); fq_poly_print_pretty(b, "x", ctx); flint_printf("\n");
                    flint_printf("s: "); fq_poly_print_pretty(s, "x", ctx); flint_printf("\n");
                    flint_printf("t: "); fq_poly_print_pretty(t, "x", ctx); flint_printf("\n");
                    flint_printf("g: "); fq_poly_print_pretty(g, "x", ctx); flint_printf("\n");
                    flint_printf("s*a+t*b: "); fq_poly_print_pretty(t1, "x", ctx); flint_printf("\n");
                    flint_abort();
                }
            }
        }

        fq_poly_clear(a, ctx);
        fq_poly_clear(b, ctx);
        fq_poly_clear(g, ctx);
        fq_poly_clear(s, ctx);
        fq_poly_clear(t, ctx);
        fq_poly_clear(t1, ctx);
        fq_poly_clear(t2, ctx);
        fq_clear(f, ctx);
        fq_ctx_clear(ctx);
        fmpz_mod_poly_clear(h);
        fmpz_clear(p);
    }

    FLINT_TEST_CLEANUP(state);

    flint_printf("PASS\n");
    return 0;
}

