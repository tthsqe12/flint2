/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fq_nmod_mpoly_factor.h"


/* check total number of factors with multiplicity is between lower and upper */
void check_omega(slong lower, slong upper, const fq_nmod_mpoly_t p, const fq_nmod_mpoly_ctx_t ctx)
{
    slong i;
    fq_nmod_mpoly_t q;
    fq_nmod_mpoly_factor_t g, h;
    fmpz_t omega;

    fmpz_init(omega);
    fq_nmod_mpoly_factor_init(g, ctx);
    fq_nmod_mpoly_factor_init(h, ctx);
    fq_nmod_mpoly_init(q, ctx);

flint_printf("p: ");
fq_nmod_mpoly_print_pretty(p, NULL, ctx);
flint_printf("\n");

flint_printf("doing factor_wang\n");
    if (!fq_nmod_mpoly_factor_wang(g, p, ctx))
    {
        flint_printf("FAIL:\ncheck factorization 1 could be computed\n");
        flint_abort();        
    }


flint_printf("g: ");
fq_nmod_mpoly_factor_print_pretty(g, NULL, ctx);
flint_printf("\n");



flint_printf("doing factor\n");
    if (!fq_nmod_mpoly_factor(h, p, ctx))
    {
        flint_printf("FAIL:\ncheck factorization 2 could be computed\n");
        flint_abort();        
    }

flint_printf("done factor\n");

    for (i = 0; i < g->num; i++)
    {
        if (g->poly[i].length < 1 || !fq_nmod_is_one(g->poly[i].coeffs + 0, ctx->fqctx))
        {
            flint_printf("FAIL:\nfactorization is not unit normal\n");
            flint_abort();
        }
    }

    fmpz_zero(omega);
    for (i = 0; i < g->num; i++)
        fmpz_add(omega, omega, g->exp + i);

    if (fmpz_cmp_si(omega, lower) < 0 || fmpz_cmp_si(omega, upper) > 0)
    {
        flint_printf("FAIL:\nfactorization has wrong number of factors\n");
        flint_abort();        
    }

    fq_nmod_mpoly_factor_expand(q, g, ctx);

flint_printf("q: ");
fq_nmod_mpoly_print_pretty(q, NULL, ctx);
flint_printf("\n");


    if (!fq_nmod_mpoly_equal(q, p, ctx))
    {
        flint_printf("FAIL:\nfactorization does not match original polynomial\n");
        flint_abort();        
    }

    fq_nmod_mpoly_factor_sort(g, ctx);
    fq_nmod_mpoly_factor_sort(h, ctx);
    if (fq_nmod_mpoly_factor_cmp(g, h, ctx) != 0)
    {
        flint_printf("factorizations do not match\n");
        flint_abort();        
    }

    for (i = 0; i < g->num; i++)
    {
        fq_nmod_mpoly_factor(h, g->poly + i, ctx);
        if (h->num != 1 || !fmpz_is_one(h->exp + 0))
        {
            flint_printf("FAIL:\nfactor is reducible\n");
            flint_abort();
        }
    }

    fq_nmod_mpoly_clear(q, ctx);
    fq_nmod_mpoly_factor_clear(g, ctx);
    fq_nmod_mpoly_factor_clear(h, ctx);
    fmpz_clear(omega);
}


int
main(void)
{
    slong i, j, tmul = 15;
    FLINT_TEST_INIT(state);

    flint_printf("factor_wang....");
    fflush(stdout);

    for (i = 0; i < tmul * flint_test_multiplier(); i++)
    {
        slong lower;
        fq_nmod_mpoly_ctx_t ctx;
        fq_nmod_mpoly_t a, t;
        slong nfacs, len;
        ulong expbound, powbound, pow;

        fq_nmod_mpoly_ctx_init_rand(ctx, state, 5, FLINT_BITS, 5);
        fq_nmod_mpoly_init(a, ctx);
        fq_nmod_mpoly_init(t, ctx);

        nfacs = 1 + n_randint(state, 6);
        powbound = 1 + n_randint(state, 3);
        powbound = 1 + n_randint(state, powbound);
        expbound = 2 + 25/nfacs/ctx->minfo->nvars;

flint_printf("i: %wd, powbound: %wu, expbound: %wu\n", i, powbound, expbound);
fflush(stdout);


        lower = 0;
        fq_nmod_mpoly_one(a, ctx);
        for (j = 0; j < nfacs; j++)
        {
            pow = 1 + n_randint(state, powbound);
            len = 1 + n_randint(state, 35/pow/nfacs);
            fq_nmod_mpoly_randtest_bound(t, state, len, expbound, ctx);
            if (fq_nmod_mpoly_is_zero(t, ctx))
                fq_nmod_mpoly_one(t, ctx);
            if (!fq_nmod_mpoly_is_fq_nmod(t, ctx))
                lower += pow;
flint_printf("[ length %wd ]\n", t->length);
fflush(stdout);

flint_printf("*(");
fq_nmod_mpoly_print_pretty(t, NULL, ctx);
flint_printf(")^%wu\n", pow);
fflush(stdout);

            fq_nmod_mpoly_pow_ui(t, t, pow, ctx);
            fq_nmod_mpoly_mul(a, a, t, ctx);
        }

        check_omega(lower, WORD_MAX, a, ctx);

        fq_nmod_mpoly_clear(t, ctx);
        fq_nmod_mpoly_clear(a, ctx);
        fq_nmod_mpoly_ctx_clear(ctx);
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}
