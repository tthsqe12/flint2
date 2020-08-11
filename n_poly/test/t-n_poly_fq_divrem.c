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

int
main(void)
{
    slong i, j;
    FLINT_TEST_INIT(state);

    flint_printf("n_poly_fq_divrem....");
    fflush(stdout);

    if (0) {
        ulong t1, t2, t3, t4;
        slong deg = 40;
        slong n = 1000000/deg;
        fmpz_t p;
        fq_nmod_ctx_t ctx;
        n_poly_struct * a, * b;
        n_poly_t q, r;
        fq_nmod_poly_struct * A, * B;
        fq_nmod_poly_t Q, R;
        n_poly_stack_t St;
        nmod_poly_t modulus;
        timeit_t timer;

        n_poly_stack_init(St);

        nmod_poly_init(modulus, UWORD(9223372036854775837));

        nmod_poly_randtest_irreducible(modulus, state, 2 + 1);
        nmod_poly_make_monic(modulus, modulus);

flint_printf("\nmodulus: "); nmod_poly_print_pretty(modulus, "#"); flint_printf("\n");

        fq_nmod_ctx_init_modulus(ctx, modulus, "#");

        n_poly_init(q);
        n_poly_init(r);
        fq_nmod_poly_init(Q, ctx);
        fq_nmod_poly_init(R, ctx);

        a = (n_poly_struct *) flint_malloc(n*sizeof(n_poly_struct));
        b = (n_poly_struct *) flint_malloc(n*sizeof(n_poly_struct));
        A = (fq_nmod_poly_struct *) flint_malloc(n*sizeof(fq_nmod_poly_struct));
        B = (fq_nmod_poly_struct *) flint_malloc(n*sizeof(fq_nmod_poly_struct));
        for (i = 0; i < n; i++)
        {
            slong k = n_randint(state, n_randint(state, deg) + 1) + 1;
            n_poly_init(a + i);
            n_poly_init(b + i);
            fq_nmod_poly_init(A + i, ctx);
            fq_nmod_poly_init(B + i, ctx);
            n_poly_fq_randtest(a + i, state, k + n_randint(state, k) + 1, ctx);
            n_poly_fq_randtest(b + i, state, k, ctx);
            n_poly_fq_get_fq_nmod_poly(A + i, a + i, ctx);
            n_poly_fq_get_fq_nmod_poly(B + i, b + i, ctx);
        }

        for (j = 0; j < 4; j++)
        {
            timeit_start(timer);
            for (i = 0; i + 1 < n; i++)
                n_poly_fq_divrem_divconquer_(q, r, a + i, b + i, ctx, St);
            timeit_stop(timer);
            t1 = timer->wall;
            flint_printf("\n  new: %wd\n", t1);

            timeit_start(timer);
            for (i = 0; i + 1 < n; i++)
                fq_nmod_poly_divrem(Q, R, A + i, B + i, ctx);
            timeit_stop(timer);
            t2 = timer->wall;
            flint_printf("  old: %wd\n", t2);
            flint_printf("ratio: %f\n", (double)(t2)/(double)(t1));
        }

        n_poly_stack_clear(St);

    }

    for (i = 0; i < 2000 * flint_test_multiplier(); i++)
    {
        fq_nmod_ctx_t ctx;
        n_poly_t a, b, c, d, e, f, q, r;
        fq_nmod_poly_t A, B, Q, R;

        fq_nmod_ctx_randtest(ctx, state);

        n_poly_init(a);
        n_poly_init(b);
        n_poly_init(c);
        n_poly_init(d);
        n_poly_init(e);
        n_poly_init(f);
        n_poly_init(q);
        n_poly_init(r);
        fq_nmod_poly_init(A, ctx);
        fq_nmod_poly_init(B, ctx);
        fq_nmod_poly_init(Q, ctx);
        fq_nmod_poly_init(R, ctx);

        for (j = 0; j < 10; j++)
        {
            n_poly_fq_randtest(a, state, n_randint(state, 20), ctx);
            n_poly_fq_randtest(b, state, n_randint(state, 15), ctx);
            n_poly_fq_randtest(c, state, n_randint(state, 20), ctx);
            n_poly_fq_randtest(d, state, n_randint(state, 20), ctx);
            n_poly_fq_randtest(e, state, n_randint(state, 20), ctx);
            n_poly_fq_randtest(f, state, n_randint(state, 20), ctx);

            if (n_poly_is_zero(b))
                n_poly_fq_one(b, ctx);

            n_poly_fq_get_fq_nmod_poly(A, a, ctx);
            n_poly_fq_get_fq_nmod_poly(B, b, ctx);
            fq_nmod_poly_divrem(Q, R, A, B, ctx);
            n_poly_fq_set_fq_nmod_poly(q, Q, ctx);
            n_poly_fq_set_fq_nmod_poly(r, R, ctx);

            n_poly_fq_divrem(c, d, a, b, ctx);
            n_poly_fq_divrem(c, d, a, b, ctx);
            if (!n_poly_fq_equal(c, q, ctx) ||
                !n_poly_fq_equal(d, r, ctx))
            {
flint_printf("\n");
flint_printf("q: "); n_poly_fq_print_pretty(q, "x", ctx); flint_printf("\n");
flint_printf("r: "); n_poly_fq_print_pretty(r, "x", ctx); flint_printf("\n");
flint_printf("\n");
flint_printf("c: "); n_poly_fq_print_pretty(c, "x", ctx); flint_printf("\n");
flint_printf("d: "); n_poly_fq_print_pretty(d, "x", ctx); flint_printf("\n");

                flint_printf("FAIL\n i = %wd, j = %wd\n", i, j);
                flint_abort();
            }
        }

        n_poly_clear(a);
        n_poly_clear(b);
        n_poly_clear(c);
        n_poly_clear(d);
        n_poly_clear(e);
        n_poly_clear(f);
        n_poly_clear(q);
        n_poly_clear(r);
        fq_nmod_poly_clear(A, ctx);
        fq_nmod_poly_clear(B, ctx);
        fq_nmod_poly_clear(Q, ctx);
        fq_nmod_poly_clear(R, ctx);
        fq_nmod_ctx_clear(ctx);
    }

    FLINT_TEST_CLEANUP(state);

    flint_printf("PASS\n");
    return 0;
}