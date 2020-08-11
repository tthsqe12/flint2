/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fq_nmod_mpoly_factor.h"
#include "profiler.h"


/* check total number of factors with multiplicity is between lower and upper */
slong check_omega(slong lower, slong upper, const fq_nmod_mpoly_t p, const fq_nmod_mpoly_ctx_t ctx)
{
    slong i, res;
    fq_nmod_mpoly_t q;
    fq_nmod_mpoly_factor_t g, g2;
    fmpz_t omega;

    flint_printf("checking %wd <= # <= %wd: length %wd in %wd vars\n",
                  lower, FLINT_MIN(9999, upper), p->length, ctx->minfo->nvars);

    fmpz_init(omega);
    fq_nmod_mpoly_factor_init(g, ctx);
    fq_nmod_mpoly_factor_init(g2, ctx);
    fq_nmod_mpoly_init(q, ctx);

    if (!fq_nmod_mpoly_factor(g, p, ctx))
    {
        flint_printf("FAIL:\ncheck factorization could be computed\n");
        flint_abort();        
    }

/*
flint_printf("p: "); nmod_mpoly_print_pretty(p, NULL, ctx); flint_printf("\n");
flint_printf("g: "); nmod_mpoly_factor_print_pretty(g, NULL, ctx); flint_printf("\n");
*/

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
    if (!fq_nmod_mpoly_equal(q, p, ctx))
    {
        flint_printf("FAIL:\nfactorization does not match original polynomial\n");
        flint_abort();        
    }

    for (i = 0; i < g->num; i++)
    {
        fq_nmod_mpoly_factor(g2, g->poly + i, ctx);
        if (g2->num != 1 || !fmpz_is_one(g2->exp + 0))
        {
            flint_printf("FAIL:\nfactor is reducible\n");
            flint_abort();
        }
    }

    fq_nmod_mpoly_clear(q, ctx);
    fq_nmod_mpoly_factor_clear(g, ctx);
    fq_nmod_mpoly_factor_clear(g2, ctx);
    res = fmpz_get_si(omega);
    fmpz_clear(omega);
    return res;
}


void init_vars(char *** vars, slong * nvars, const char * s)
{
    slong i, sn = strlen(s);

    (*vars) = (char **) flint_malloc(sn * sizeof(char *));
    for (i = 0; i < sn; i++)
    {
        (*vars)[i] = flint_calloc(sn + 1, sizeof(char));
    }

    (*nvars) = 0;
    for (i = 0; i < sn; i++)
    {
        if (s[i] == ' ')
        {
            if (strlen((*vars)[(*nvars)]) != 0)
            {
               (*nvars)++;
            }
        }
        else
        {
            strncat((*vars)[(*nvars)], s + i, 1);
        }
    }

    (*nvars)++;
}

void clear_vars(char *** vars, slong * nvars, const char * s)
{
    slong i, sn = strlen(s);
    for (i = 0; i < sn; i++)
        flint_free((*vars)[i]);
    flint_free((*vars));
}


int
main(void)
{
int check = 0;
    slong i, j, tmul = 30;
    slong total;
    FLINT_TEST_INIT(state);

    flint_printf("factor....");
    fflush(stdout);

if (check)
{
    total = 0;
    for (i = 0; i < tmul * flint_test_multiplier(); i++)
    {
        slong lower;
        fq_nmod_mpoly_ctx_t ctx;
        fq_nmod_mpoly_t a, t;
        slong nfacs, len;
        ulong expbounds[2];
        ulong expbound;

        fq_nmod_mpoly_ctx_init_rand(ctx, state, 2, FLINT_BITS, 5);
        fq_nmod_mpoly_init(a, ctx);
        fq_nmod_mpoly_init(t, ctx);

        nfacs = 2 + n_randint(state, 10);
        expbound = 3 + n_randint(state, 3 + 40/nfacs);

        lower = 0;
        fq_nmod_mpoly_one(a, ctx);
        for (j = 0; j < nfacs; j++)
        {
            expbounds[0] = 2 + n_randint(state, expbound);
            expbounds[1] = 2 + n_randint(state, expbound);
            len = expbounds[0] + expbounds[1] + 50;
            fq_nmod_mpoly_randtest_bounds(t, state, len, expbounds, ctx);
            if (fq_nmod_mpoly_is_zero(t, ctx))
                fq_nmod_mpoly_one(t, ctx);
            lower += !fq_nmod_mpoly_is_fq_nmod(t, ctx);
            fq_nmod_mpoly_mul(a, a, t, ctx);
        }

flint_printf("1:%wd ", i);
        total += check_omega(lower, WORD_MAX, a, ctx);

        fq_nmod_mpoly_clear(t, ctx);
        fq_nmod_mpoly_clear(a, ctx);
        fq_nmod_mpoly_ctx_clear(ctx);
    }
flint_printf("**********total number of bvar factors: %wd ******\n", total);

    if (total != 2805)
    {
        flint_printf("total number of bvar factors should be 2805\n");
        flint_abort();
    }
}

    total = 0;
    for (i = 0; i < tmul * flint_test_multiplier(); i++)
    {
        slong lower;
        fq_nmod_mpoly_ctx_t ctx;
        fq_nmod_mpoly_t a, t;
        slong nfacs, len;
        ulong expbound;

        fq_nmod_mpoly_ctx_init_rand(ctx, state, 8, FLINT_BITS, 5);
        fq_nmod_mpoly_init(a, ctx);
        fq_nmod_mpoly_init(t, ctx);

        nfacs = 1 + (5 + n_randint(state, 5))/ctx->minfo->nvars;
        expbound = 2 + 60/nfacs/ctx->minfo->nvars;

        lower = 0;
        fq_nmod_mpoly_one(a, ctx);
        for (j = 0; j < nfacs; j++)
        {
            len = 1 + n_randint(state, 20);
            fq_nmod_mpoly_randtest_bound(t, state, len, expbound, ctx);
            if (fq_nmod_mpoly_is_zero(t, ctx))
                fq_nmod_mpoly_one(t, ctx);
            lower += !fq_nmod_mpoly_is_fq_nmod(t, ctx);
            fq_nmod_mpoly_mul(a, a, t, ctx);
        }

flint_printf("2:%wd ", i);
        total += check_omega(lower, WORD_MAX, a, ctx);

        fq_nmod_mpoly_clear(t, ctx);
        fq_nmod_mpoly_clear(a, ctx);
        fq_nmod_mpoly_ctx_clear(ctx);
    }
flint_printf("**********total number of mvar factors: %wd ******\n", total);

if (check)
{
    if (total != 3443)
    {
        flint_printf("total number of mvar factors should be 3443\n");
        flint_abort();
    }
}

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}
