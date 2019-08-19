/*
    Copyright (C) 2019 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "nmod_mpoly_factor.h"
#include "profiler.h"

int
main(void)
{
    slong i, j, k;
    FLINT_TEST_INIT(state);

    flint_printf("factor....");
    fflush(stdout);

    {
timeit_t timer;
        nmod_mpoly_ctx_t ctx;
        nmod_mpoly_t a;
        nmod_mpoly_factor_t fac;
        const char * vars[] = {"x", "y", "z", "t"};

        nmod_mpoly_ctx_init(ctx, 2, ORD_LEX, 2);
        nmod_mpoly_init(a, ctx);
        nmod_mpoly_set_str_pretty(a, "(x+1)^16385", vars, ctx);

printf("****nmod_mpoly****\n");
printf("    factoring: "); nmod_mpoly_print_pretty(a, vars, ctx); printf("\n");

        nmod_mpoly_factor_init(fac, ctx);
timeit_start(timer);
        nmod_mpoly_factor(fac, a, 1, ctx);
timeit_stop(timer);
printf("factorization: "); nmod_mpoly_factor_print_pretty(fac, vars, ctx); printf("\n");
flint_printf("         time: %wd\n", a->length, timer->wall);

        nmod_mpoly_factor_clear(fac, ctx);
        nmod_mpoly_clear(a, ctx);
        nmod_mpoly_ctx_clear(ctx);
    }

	{
timeit_t timer;
		nmod_poly_factor_t fac;
		nmod_poly_t a;

		nmod_poly_factor_init(fac);
		nmod_poly_init(a, 2);

		nmod_poly_set_coeff_ui(a, 1, 1);
		nmod_poly_set_coeff_ui(a, 0, 1);
		nmod_poly_pow(a, a, 16385);

printf("****nmod_poly****\n");
printf("    factoring: "); nmod_poly_print_pretty(a, "x"); printf("\n");
timeit_start(timer);
        nmod_poly_factor(fac, a);
timeit_stop(timer);
printf("factorization: "); nmod_poly_factor_print(fac);
flint_printf("         time: %wd\n", timer->wall);

		nmod_poly_factor_clear(fac);
		nmod_poly_clear(a);
	}

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}
