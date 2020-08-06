/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include <stdlib.h>
#include "n_poly.h"
#include "profiler.h"

int main(void)
{
    slong i, j;
    FLINT_TEST_INIT(state);

    flint_printf("n_poly_fq_mul....");
    fflush(stdout);

    if (0) {
        slong t1, t2;
        slong deg = 20;
        slong n = 5000000/deg;
        fq_nmod_ctx_t ctx;
        n_poly_struct * a;
        n_poly_t b;
        fq_nmod_poly_struct * A;
        fq_nmod_poly_t B;
        n_poly_stack_t St;
        nmod_poly_t modulus;
        timeit_t timer;

        n_poly_stack_init(St);

        nmod_poly_init(modulus, UWORD(9223372036854775837));

        nmod_poly_randtest_irreducible(modulus, state, 8 + 1);
        nmod_poly_make_monic(modulus, modulus);

flint_printf("\nmodulus: "); nmod_poly_print_pretty(modulus, "#"); flint_printf("\n");

        fq_nmod_ctx_init_modulus(ctx, modulus, "#");

        n_poly_init(b);
        fq_nmod_poly_init(B, ctx);

        a = (n_poly_struct *) flint_malloc(n*sizeof(n_poly_struct));
        A = (fq_nmod_poly_struct *) flint_malloc(n*sizeof(fq_nmod_poly_struct));
        for (i = 0; i < n; i++)
        {
            n_poly_init(a + i);
            fq_nmod_poly_init(A + i, ctx);
            n_poly_fq_randtest(a + i, state,
                         n_randint(state, n_randint(state, deg) + 1) + 1, ctx);
            n_poly_fq_get_fq_nmod_poly(A + i, a + i, ctx);
        }

        for (j = 0; j < 4; j++)
        {
            timeit_start(timer);
            for (i = 0; i + 1 < n; i++)
                n_poly_fq_mul_(b, a + i, a + i + 1, ctx, St);
            timeit_stop(timer);
            t1 = timer->wall;
            flint_printf("\n  new: %wd\n", t1);

            timeit_start(timer);
            for (i = 0; i + 1 < n; i++)
                fq_nmod_poly_mul(B, A + i, A + i + 1, ctx);
            timeit_stop(timer);
            t2 = timer->wall;
            flint_printf("  old: %wd\n", t2);
            flint_printf("ratio: %f\n", (double)(t2)/(double)(t1));
        }

        n_poly_stack_clear(St);
    }

    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fq_nmod_ctx_t ctx;
        n_poly_t a, b, c, d, e;
        fq_nmod_poly_t A, B, C;

        fq_nmod_ctx_randtest(ctx, state);

        n_poly_init(a);
        n_poly_init(b);
        n_poly_init(c);
        n_poly_init(d);
        n_poly_init(e);
        fq_nmod_poly_init(A, ctx);
        fq_nmod_poly_init(B, ctx);
        fq_nmod_poly_init(C, ctx);

        for (j = 0; j < 10; j++)
        {
            n_poly_fq_randtest(a, state, n_randint(state, 20), ctx);
            n_poly_fq_randtest(b, state, n_randint(state, 20), ctx);
            n_poly_fq_randtest(c, state, n_randint(state, 20), ctx);
            n_poly_fq_randtest(d, state, n_randint(state, 20), ctx);

            n_poly_fq_get_fq_nmod_poly(B, b, ctx);
            n_poly_fq_get_fq_nmod_poly(C, c, ctx);
            fq_nmod_poly_mul(A, B, C, ctx);
            n_poly_fq_set_fq_nmod_poly(a, A, ctx);

            n_poly_fq_mul(d, b, c, ctx);
            n_poly_fq_set(e, b, ctx);
            n_poly_fq_mul(b, b, c, ctx);
            n_poly_fq_mul(c, e, c, ctx);

            if (!n_poly_fq_equal(a, d, ctx) ||
                !n_poly_fq_equal(a, b, ctx) ||
                !n_poly_fq_equal(a, c, ctx))
            {
                flint_printf("FAIL\n i = %wd, j = %wd\n", i, j);
                flint_abort();
            }
        }

        n_poly_clear(a);
        n_poly_clear(b);
        n_poly_clear(c);
        n_poly_clear(d);
        n_poly_clear(e);
        fq_nmod_poly_clear(A, ctx);
        fq_nmod_poly_clear(B, ctx);
        fq_nmod_poly_clear(C, ctx);
        fq_nmod_ctx_clear(ctx);
    }

    FLINT_TEST_CLEANUP(state);

    flint_printf("PASS\n");
    return 0;
}
