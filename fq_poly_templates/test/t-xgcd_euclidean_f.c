/*
    Copyright (C) 2012 Sebastian Pancratz
    Copyright (C) 2013 Mike Hansen

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#ifdef T

#include "templates.h"

int
main(void)
{
    int i, result;
    FLINT_TEST_INIT(state);

    flint_printf("xgcd_euclidean_f....");
    fflush(stdout);

    /* Generic case, most likely co-prime arguments ***************************** */

    /* Compare with result from GCD and check correctness */
    for (i = 0; i < 50 * flint_test_multiplier(); i++)
    {
        TEMPLATE(T, ctx_t) ctx;
        TEMPLATE(T, t) f1, f2;
        TEMPLATE(T, poly_t) a, b, d, g, s, t, v, w;

        TEMPLATE(T, ctx_randtest_reducible) (ctx, state);

        TEMPLATE(T, init) (f1, ctx);
        TEMPLATE(T, init) (f2, ctx);
        TEMPLATE(T, poly_init) (a, ctx);
        TEMPLATE(T, poly_init) (b, ctx);
        TEMPLATE(T, poly_init) (d, ctx);
        TEMPLATE(T, poly_init) (g, ctx);
        TEMPLATE(T, poly_init) (s, ctx);
        TEMPLATE(T, poly_init) (t, ctx);
        TEMPLATE(T, poly_init) (v, ctx);
        TEMPLATE(T, poly_init) (w, ctx);
        TEMPLATE(T, poly_randtest) (a, state, n_randint(state, 10), ctx);
        TEMPLATE(T, poly_randtest) (b, state, n_randint(state, 10), ctx);

        TEMPLATE(T, poly_gcd_euclidean_f) (f1, d, a, b, ctx);
        TEMPLATE(T, poly_xgcd_euclidean_f) (f2, g, s, t, a, b, ctx);

        if (!TEMPLATE(T, is_one) (f2, ctx) || !TEMPLATE(T, is_one) (f1, ctx))
        {
            result = TEMPLATE(T, equal) (f1, f2, ctx);
        }
        else
        {
            TEMPLATE(T, poly_mul) (v, s, a, ctx);
            TEMPLATE(T, poly_mul) (w, t, b, ctx);
            TEMPLATE(T, poly_add) (w, v, w, ctx);

            result = (TEMPLATE(T, poly_equal) (d, g, ctx) &&
                      TEMPLATE(T, poly_equal) (g, w, ctx));
        }

        if (!result)
        {
            flint_printf("FAIL:\n");
            TEMPLATE(T, ctx_print) (ctx);
            TEMPLATE(T, print_pretty) (f1, ctx), flint_printf("\n\n");
            TEMPLATE(T, print_pretty) (f2, ctx), flint_printf("\n\n");
            TEMPLATE(T, poly_print_pretty) (a, "x", ctx), flint_printf("\n\n");
            TEMPLATE(T, poly_print_pretty) (b, "x", ctx), flint_printf("\n\n");
            TEMPLATE(T, poly_print_pretty) (d, "x", ctx), flint_printf("\n\n");
            TEMPLATE(T, poly_print_pretty) (g, "x", ctx), flint_printf("\n\n");
            TEMPLATE(T, poly_print_pretty) (s, "x", ctx), flint_printf("\n\n");
            TEMPLATE(T, poly_print_pretty) (t, "x", ctx), flint_printf("\n\n");
            TEMPLATE(T, poly_print_pretty) (v, "x", ctx), flint_printf("\n\n");
            TEMPLATE(T, poly_print_pretty) (w, "x", ctx), flint_printf("\n\n");
            abort();
        }

        TEMPLATE(T, clear) (f1, ctx);
        TEMPLATE(T, clear) (f2, ctx);
        TEMPLATE(T, poly_clear) (a, ctx);
        TEMPLATE(T, poly_clear) (b, ctx);
        TEMPLATE(T, poly_clear) (d, ctx);
        TEMPLATE(T, poly_clear) (g, ctx);
        TEMPLATE(T, poly_clear) (s, ctx);
        TEMPLATE(T, poly_clear) (t, ctx);
        TEMPLATE(T, poly_clear) (v, ctx);
        TEMPLATE(T, poly_clear) (w, ctx);
        TEMPLATE(T, ctx_clear) (ctx);
    }

    /* Special case, arguments share a factor ******************************* */

    /* Compare with result from GCD and check correctness */
    for (i = 0; i < 50 * flint_test_multiplier(); i++)
    {
        TEMPLATE(T, ctx_t) ctx;
        TEMPLATE(T, t) f1, f2;
        TEMPLATE(T, poly_t) a, b, d, f, g, s, t, v, w;

        TEMPLATE(T, ctx_randtest_reducible) (ctx, state);

        TEMPLATE(T, init) (f1, ctx);
        TEMPLATE(T, init) (f2, ctx);

        TEMPLATE(T, poly_init) (a, ctx);
        TEMPLATE(T, poly_init) (b, ctx);
        TEMPLATE(T, poly_init) (d, ctx);
        TEMPLATE(T, poly_init) (f, ctx);
        TEMPLATE(T, poly_init) (g, ctx);
        TEMPLATE(T, poly_init) (s, ctx);
        TEMPLATE(T, poly_init) (t, ctx);
        TEMPLATE(T, poly_init) (v, ctx);
        TEMPLATE(T, poly_init) (w, ctx);
        TEMPLATE(T, poly_randtest) (a, state, n_randint(state, 10), ctx);
        TEMPLATE(T, poly_randtest) (b, state, n_randint(state, 10), ctx);
        TEMPLATE(T, poly_randtest) (f, state, n_randint(state, 20), ctx);
        TEMPLATE(T, poly_mul) (a, a, f, ctx);
        TEMPLATE(T, poly_mul) (b, b, f, ctx);

        TEMPLATE(T, poly_gcd_euclidean_f) (f1, d, a, b, ctx);
        TEMPLATE(T, poly_xgcd_euclidean_f) (f2, g, s, t, a, b, ctx);

        if (!TEMPLATE(T, is_one) (f2, ctx))
        {
            result = TEMPLATE(T, equal) (f1, f2, ctx);
        }
        else
        {
            TEMPLATE(T, poly_mul) (v, s, a, ctx);
            TEMPLATE(T, poly_mul) (w, t, b, ctx);
            TEMPLATE(T, poly_add) (w, v, w, ctx);

            result = (TEMPLATE(T, poly_equal) (d, g, ctx) &&
                      TEMPLATE(T, poly_equal) (g, w, ctx));
        }

        if (!result)
        {
            flint_printf("FAIL:\n");
            TEMPLATE(T, print_pretty) (f1, ctx), flint_printf("\n\n");
            TEMPLATE(T, print_pretty) (f2, ctx), flint_printf("\n\n");
            TEMPLATE(T, poly_print_pretty) (a, "x", ctx);
            flint_printf("\n\n");
            TEMPLATE(T, poly_print_pretty) (b, "x", ctx);
            flint_printf("\n\n");
            TEMPLATE(T, poly_print_pretty) (d, "x", ctx);
            flint_printf("\n\n");
            TEMPLATE(T, poly_print_pretty) (f, "x", ctx);
            flint_printf("\n\n");
            TEMPLATE(T, poly_print_pretty) (g, "x", ctx);
            flint_printf("\n\n");
            TEMPLATE(T, poly_print_pretty) (s, "x", ctx);
            flint_printf("\n\n");
            TEMPLATE(T, poly_print_pretty) (t, "x", ctx);
            flint_printf("\n\n");
            TEMPLATE(T, poly_print_pretty) (v, "x", ctx);
            flint_printf("\n\n");
            TEMPLATE(T, poly_print_pretty) (w, "x", ctx);
            flint_printf("\n\n");
            abort();
        }

        TEMPLATE(T, clear) (f1, ctx);
        TEMPLATE(T, clear) (f2, ctx);
        TEMPLATE(T, poly_clear) (a, ctx);
        TEMPLATE(T, poly_clear) (b, ctx);
        TEMPLATE(T, poly_clear) (d, ctx);
        TEMPLATE(T, poly_clear) (f, ctx);
        TEMPLATE(T, poly_clear) (g, ctx);
        TEMPLATE(T, poly_clear) (s, ctx);
        TEMPLATE(T, poly_clear) (t, ctx);
        TEMPLATE(T, poly_clear) (v, ctx);
        TEMPLATE(T, poly_clear) (w, ctx);
        TEMPLATE(T, ctx_clear) (ctx);
    }

    FLINT_TEST_CLEANUP(state);

    flint_printf("PASS\n");
    return 0;
}

#endif
