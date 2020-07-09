/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "nmod_mpoly_factor.h"
#include "profiler.h"


/* check total number of factors with multiplicity is between lower and upper */
void check_omega(slong lower, slong upper, const nmod_mpoly_t p, const nmod_mpoly_ctx_t ctx)
{
    slong i;
    nmod_mpoly_t q;
    nmod_mpoly_factor_t g, g2;
    fmpz_t omega;

    flint_printf("checking %wd <= # <= %wd: length %wd in %wd vars\n",
                                   lower, upper, p->length, ctx->minfo->nvars);

    fmpz_init(omega);
    nmod_mpoly_factor_init(g, ctx);
    nmod_mpoly_factor_init(g2, ctx);
    nmod_mpoly_init(q, ctx);

    nmod_mpoly_factor(g, p, ctx);

    if (nmod_mpoly_factor_fix_units(g, ctx))
    {
        flint_printf("FAIL:\nfactorization is not unit normal\n");
        flint_abort();
    }

    fmpz_zero(omega);
    for (i = 0; i < g->num; i++)
        fmpz_add(omega, omega, g->exp + i);

    if (fmpz_cmp_si(omega, lower) < 0 || fmpz_cmp_si(omega, upper) > 0)
    {
        flint_printf("FAIL:\nfactorization has wrong number of factors\n");
        flint_abort();        
    }

    nmod_mpoly_factor_expand(q, g, ctx);
    if (!nmod_mpoly_equal(q, p, ctx))
    {
        flint_printf("FAIL:\nfactorization does not match original polynomial\n");
        flint_abort();        
    }

    for (i = 0; i < g->num; i++)
    {
        nmod_mpoly_factor(g2, g->poly + i, ctx);
        if (g2->num != 1 || !fmpz_is_one(g2->exp + 0))
        {
            flint_printf("FAIL:\nfactor is reducible\n");
            flint_abort();
        }
    }

    nmod_mpoly_clear(q, ctx);
    nmod_mpoly_factor_clear(g, ctx);
    nmod_mpoly_factor_clear(g2, ctx);
    fmpz_clear(omega);
}

int
main(void)
{
    slong i, j, tmul = 10;
    FLINT_TEST_INIT(state);

    flint_printf("factor....");
    fflush(stdout);

    /* check random factors */
    for (i = 0; i < tmul * flint_test_multiplier(); i++)
    {
        slong lower;
        nmod_mpoly_ctx_t ctx;
        nmod_mpoly_t a, t;
        slong nfacs, len;
        ulong expbound;
        mp_limb_t p;

        p = n_randint(state, (i % 10 == 0) ? 4 : FLINT_BITS - 1) + 1;
        p = n_randbits(state, p);
        p = n_nextprime(p, 1);

        nmod_mpoly_ctx_init_rand(ctx, state, 6, p);

        nmod_mpoly_init(a, ctx);
        nmod_mpoly_init(t, ctx);

        nfacs = 1 + (5 + n_randint(state, 5))/ctx->minfo->nvars;
        expbound = 3 + 40/nfacs/ctx->minfo->nvars;

        lower = 0;
        nmod_mpoly_one(a, ctx);
        for (j = 0; j < nfacs; j++)
        {
            len = 1 + n_randint(state, 7);
            nmod_mpoly_randtest_bound(t, state, len, expbound, ctx);
            if (nmod_mpoly_is_zero(t, ctx))
                nmod_mpoly_one(t, ctx);
            lower += !nmod_mpoly_is_ui(t, ctx);
            nmod_mpoly_mul(a, a, t, ctx);
        }

flint_printf("%wd ", i);
        check_omega(lower, WORD_MAX, a, ctx);

        nmod_mpoly_clear(t, ctx);
        nmod_mpoly_clear(a, ctx);
        nmod_mpoly_ctx_clear(ctx);
    }

#if 0
    {
        nmod_mpoly_ctx_t ctx;
        nmod_mpoly_t a;
        nmod_mpoly_factor_t fac;
        const char * vars[] = {"x", "y", "z", "t"};
        nmod_mpoly_ctx_init(ctx, 3, ORD_LEX, 2);
        nmod_mpoly_init(a, ctx);
        nmod_mpoly_set_str_pretty(a,
"(z^8*x^8+x^1+y^16+y^1+z^8+z^3)*((y^4+z^3+z)*x^8+x^1+y^16+y^1+z^8+z^3)"
, vars, ctx);
printf("\n******** starting example 3 vars ********\n");
printf("       factoring: "); nmod_mpoly_print_pretty(a, vars, ctx); printf("\n");
        nmod_mpoly_factor_init(fac, ctx);
        result = nmod_mpoly_factor(fac, a, 1, ctx);
        nmod_mpoly_factor_sort(fac, ctx);
printf("factorization(%d): ", result); nmod_mpoly_factor_print_pretty(fac, vars, ctx); printf("\n");
        nmod_mpoly_factor_clear(fac, ctx);
        nmod_mpoly_clear(a, ctx);
        nmod_mpoly_ctx_clear(ctx);
    }
    {
timeit_t timer;
        nmod_mpoly_ctx_t ctx;
        nmod_mpoly_t a;
        nmod_mpoly_factor_t fac;
        const char * vars[] = {"x", "y", "z", "t"};

        nmod_mpoly_ctx_init(ctx, 2, ORD_LEX, 2);
        nmod_mpoly_init(a, ctx);
        nmod_mpoly_set_str_pretty(a, "(x^8+x+y^16+y)*(x^8+x+y^4+y)", vars, ctx);

printf("\n******** starting univar nmod_mpoly ********\n");
printf("    factoring: "); nmod_mpoly_print_pretty(a, vars, ctx); printf("\n");

        nmod_mpoly_factor_init(fac, ctx);
timeit_start(timer);
        nmod_mpoly_factor(fac, a, 1, ctx);
        nmod_mpoly_factor_sort(fac, ctx);
timeit_stop(timer);
printf("factorization: "); nmod_mpoly_factor_print_pretty(fac, vars, ctx); printf("\n");
flint_printf("         time: %wd\n", timer->wall);

        nmod_mpoly_factor_clear(fac, ctx);
        nmod_mpoly_clear(a, ctx);
        nmod_mpoly_ctx_clear(ctx);
    }

    for (i = 0; i < 4; i++)
    {
timeit_t timer;
        nmod_mpoly_ctx_t ctx;
        nmod_mpoly_t a;
        nmod_mpoly_factor_t fac;
        const char * vars[] = {"x", "y", "z", "t"};
        const mp_limb_t moduli[] = {3, 7, 13, 8191};

        nmod_mpoly_ctx_init(ctx, 2, ORD_LEX, moduli[i]);
        nmod_mpoly_init(a, ctx);

flint_printf("\n******* starting bernardin example mod %wu *********\n", ctx->ffinfo->mod.n);


        nmod_mpoly_set_str_pretty(a, "(x^2+y+1)*(x+y^2+1)^2*(x*y+1)^13*(x*y+2)^14*(x+y+1)^113", vars, ctx);
/*
printf("    factoring: "); nmod_mpoly_print_pretty(a, vars, ctx); printf("\n");
*/
        nmod_mpoly_factor_init(fac, ctx);
timeit_start(timer);
        nmod_mpoly_factor(fac, a, 1, ctx);
        nmod_mpoly_factor_sort(fac, ctx);
timeit_stop(timer);
printf("factorization: "); nmod_mpoly_factor_print_pretty(fac, vars, ctx); printf("\n");
flint_printf("         time: %wd\n", timer->wall);

        nmod_mpoly_factor_clear(fac, ctx);
        nmod_mpoly_clear(a, ctx);
        nmod_mpoly_ctx_clear(ctx);
    }

    {
timeit_t timer;
        nmod_mpoly_ctx_t ctx;
        nmod_mpoly_t a;
        nmod_mpoly_factor_t fac;
        const char * vars[] = {"x", "y", "z", "t"};

        nmod_mpoly_ctx_init(ctx, 2, ORD_LEX, 2);
        nmod_mpoly_init(a, ctx);
        nmod_mpoly_set_str_pretty(a, "(x+1)^4097", vars, ctx);

printf("\n******** starting univar nmod_mpoly ********\n");
printf("    factoring: "); nmod_mpoly_print_pretty(a, vars, ctx); printf("\n");

        nmod_mpoly_factor_init(fac, ctx);
timeit_start(timer);
        nmod_mpoly_factor(fac, a, 1, ctx);
        nmod_mpoly_factor_sort(fac, ctx);
timeit_stop(timer);
printf("factorization: "); nmod_mpoly_factor_print_pretty(fac, vars, ctx); printf("\n");
flint_printf("         time: %wd\n", timer->wall);

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
		nmod_poly_pow(a, a, 4097);

printf("\n******** starting univar nmod_poly ********\n");
printf("    factoring: "); nmod_poly_print_pretty(a, "x"); printf("\n");
timeit_start(timer);
        nmod_poly_factor(fac, a);
timeit_stop(timer);
printf("factorization: "); nmod_poly_factor_print(fac);
flint_printf("         time: %wd\n", timer->wall);

		nmod_poly_factor_clear(fac);
		nmod_poly_clear(a);
	}
#endif
    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}
