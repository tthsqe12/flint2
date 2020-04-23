/*
    Copyright (C) 2019 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fmpz_mpoly_factor.h"
#include "profiler.h"




/* check factoration matches f within reason */
void check_same(fmpz_mpoly_factor_t f, const fmpz_mpoly_ctx_t ctx)
{
    fmpz_mpoly_factor_t g;
    fmpz_mpoly_t p;

    flint_printf("checking same: "); fmpz_mpoly_factor_print_pretty(f, NULL, ctx); printf("\n");

    fmpz_mpoly_factor_init(g, ctx);
    fmpz_mpoly_init(p, ctx);

    if (!fmpz_mpoly_factor_expand(p, f, ctx))
    {
        flint_printf("correct factorization could not be expanded\n");
        flint_abort();            
    }

    if (!fmpz_mpoly_factor(g, p, 1, ctx))
    {
        flint_printf("factorization could not be computed\n");
        flint_abort();            

    }

    if (fmpz_mpoly_factor_fix_units(g, ctx))
    {
        flint_printf("factorization is not unit normal\n");
        flint_abort();
    }

    fmpz_mpoly_factor_sort(g, ctx);
    fmpz_mpoly_factor_fix_units(f, ctx);
    fmpz_mpoly_factor_sort(f, ctx);

    if (!fmpz_mpoly_factor_is_same(f, g, ctx))
    {
        flint_printf("factorization does not match\n");
        flint_printf(" correct: "); fmpz_mpoly_factor_print_pretty(f, NULL, ctx); printf("\n");
        flint_printf("computed: "); fmpz_mpoly_factor_print_pretty(g, NULL, ctx); printf("\n");
        flint_abort();        
    }

    fmpz_mpoly_clear(p, ctx);
    fmpz_mpoly_factor_clear(g, ctx);
}

/* check total number of factors with multiplicity is between lower and upper */
void check_omega(slong lower, slong upper, const fmpz_mpoly_t p, const fmpz_mpoly_ctx_t ctx)
{
    slong i;
    fmpz_mpoly_t q;
    fmpz_mpoly_factor_t g;
    fmpz_t omega;

    flint_printf("checking %wd <= # <= %wd: ", lower, upper);
    if (p->length < 20)
        fmpz_mpoly_print_pretty(p, NULL, ctx);
    else
        flint_printf("length %wd in %wd vars", p->length, ctx->minfo->nvars);
    flint_printf("\n");

    fmpz_init(omega);
    fmpz_mpoly_factor_init(g, ctx);
    fmpz_mpoly_init(q, ctx);

    fmpz_mpoly_factor(g, p, 1, ctx);

    if (fmpz_mpoly_factor_fix_units(g, ctx))
    {
        flint_printf("factorization is not unit normal\n");
        flint_abort();
    }

    fmpz_zero(omega);
    for (i = 0; i < g->length; i++)
        fmpz_add(omega, omega, g->exp + i);

    if (fmpz_cmp_si(omega, lower) < 0 ||
        fmpz_cmp_si(omega, upper) > 0)
    {
        flint_printf("factorization has wrong number of factors\n");
        flint_abort();        
    }

    fmpz_mpoly_factor_expand(q, g, ctx);
    if (!fmpz_mpoly_equal(q, p, ctx))
    {
        flint_printf("factorization does not match original polynomial\n");
        flint_abort();        
    }

    fmpz_mpoly_clear(q, ctx);
    fmpz_mpoly_factor_clear(g, ctx);
    fmpz_clear(omega);
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

void check_same_str(const char * s, const char * polys)
{
    slong nvars;
    char ** vars;
    fmpz_mpoly_ctx_t ctx;
    fmpz_mpoly_factor_t f;

    init_vars(&vars, &nvars, s);

    fmpz_mpoly_ctx_init(ctx, nvars, ORD_LEX);
    fmpz_mpoly_factor_init(f, ctx);

    if (fmpz_mpoly_factor_set_str_pretty(f, polys, (const char **) vars, ctx))
    {
        flint_printf("bad string input\n");
        flint_abort();
    }

    clear_vars(&vars, &nvars, s);

    check_same(f, ctx);

    fmpz_mpoly_factor_clear(f, ctx);
    fmpz_mpoly_ctx_clear(ctx);
}


void check_omega_str(slong lower, slong upper, const char * s, const char * polys)
{
    slong nvars;
    char ** vars;
    fmpz_mpoly_ctx_t ctx;
    fmpz_mpoly_t p;

    init_vars(&vars, &nvars, s);

    fmpz_mpoly_ctx_init(ctx, nvars, ORD_LEX);
    fmpz_mpoly_init(p, ctx);

    if (fmpz_mpoly_set_str_pretty(p, polys, (const char **) vars, ctx))
    {
        flint_printf("bad string input\n");
        flint_abort();
    }

    clear_vars(&vars, &nvars, s);

    check_omega(lower, upper, p, ctx);

    fmpz_mpoly_clear(p, ctx);
    fmpz_mpoly_ctx_clear(ctx);
}




int
main(void)
{
    slong i, j, k;

    FLINT_TEST_INIT(state);

    flint_printf("factor....");

    check_same_str("x y", "(x+y)*(x-y)*(x^2-x*y+y^2)*(x^2+x*y+y^2)");
    check_same_str("x y z", "((y*z+1)*x+y+z+1)*(x^2+y^2+z^2+2)*(x^3+y^3+z^3+3)");
    check_same_str("x y z", "(x^2+y+1)*(x+y^2+1)^2*(x*y+1)^13*(x*y+2)^14*(x+y+1)^113");
    check_same_str("x y", "-2*x^2*y^4*(2*x^2*y^2-9)*(3*x^2*y-1)*(10*x*y^3-6*x*y^2+1)*(4*x^2+5*y-9)*(7*x^3 + y^2)");

    check_omega_str(3, 3, "x y", "(1+y)^9*x^9-y^9");
    check_omega_str(2, 2, "x y z", "x^9-y^9*z^3");

    /* check random factors */
    for (i = 0; i < 20 * flint_test_multiplier(); i++)
    {
        slong lower;
        fmpz_mpoly_ctx_t ctx;
        fmpz_mpoly_t a, t;
        flint_bitcnt_t coeff_bits;
        slong nfacs, len;
        ulong expbound;

        fmpz_mpoly_ctx_init_rand(ctx, state, 4);

        fmpz_mpoly_init(a, ctx);
        fmpz_mpoly_init(t, ctx);

        nfacs = 2 + (5 + n_randint(state, 5))/ctx->minfo->nvars;
        expbound = 3 + 20/ctx->minfo->nvars/nfacs;

        lower = 0;
        fmpz_mpoly_one(a, ctx);
        for (j = 0; j < nfacs; j++)
        {
            do {
                len = 1 + n_randint(state, 10);
                coeff_bits = 2 + n_randint(state, 100)/nfacs;
                fmpz_mpoly_randtest_bound(t, state, len, coeff_bits, expbound, ctx);
            } while (t->length == 0);
            lower += !fmpz_mpoly_is_fmpz(t, ctx);
            fmpz_mpoly_mul(a, a, t, ctx);
        }

        check_omega(lower, 999, a, ctx);

        fmpz_mpoly_clear(t, ctx);
        fmpz_mpoly_clear(a, ctx);
        fmpz_mpoly_ctx_clear(ctx);
    }

    /* fateman */
    for (i = 0; i <= 15; i++)
    {
        fmpz_mpoly_ctx_t ctx;
        fmpz_mpoly_t a, b;
        const char * vars[] = {"x", "y", "z", "t"};

        fmpz_mpoly_ctx_init(ctx, 4, ORD_LEX);
        fmpz_mpoly_init(a, ctx);
        fmpz_mpoly_init(b, ctx);
        fmpz_mpoly_set_str_pretty(b, "1+x+y+z+t", vars, ctx);
        fmpz_mpoly_pow_ui(b, b, i, ctx);
        fmpz_mpoly_add_ui(a, b, 1, ctx);
        fmpz_mpoly_add_ui(b, b, 2, ctx);
        fmpz_mpoly_mul(a, a, b, ctx);

        k = (i > 0);
        for (j = 1; j <= i; j++)
            if ((j%2) != 0 && (i%j) == 0)
                k++;

        check_omega(k, k, a, ctx);

        fmpz_mpoly_clear(a, ctx);
        fmpz_mpoly_clear(b, ctx);
        fmpz_mpoly_ctx_clear(ctx);
    }

    /* bivariate examples */
    {
        fmpz_mpoly_ctx_t ctx;
        fmpz_mpoly_t a, b;
        fmpz_mpoly_factor_t fac;
        const char * vars[] = {"x", "y", "z1", "z2", "z3"};
        fmpz * shift  = _fmpz_vec_init(5);
        fmpz * stride = _fmpz_vec_init(5);
        fmpz_mpoly_struct * sub[5];

        fmpz_mpoly_ctx_init(ctx, 5, ORD_LEX);
        fmpz_mpoly_init(a, ctx);
        fmpz_mpoly_init(b, ctx);
        fmpz_mpoly_factor_init(fac, ctx);

        fmpz_mpoly_set_str_pretty(b, "(x+z1+z2+z3)"
                                    "*(x+z1+z2-z3)"
                                    "*(x+z1-z2+z3)"
                                    "*(x+z1-z2-z3)"
                                    "*(x-z1+z2+z3)"
                                    "*(x-z1+z2-z3)"
                                    "*(x-z1-z2+z3)"
                                    "*(x-z1-z2-z3)", vars, ctx);

        for (i = 0; i < 5; i++)
        {
            sub[i] = (fmpz_mpoly_struct *) flint_malloc(sizeof(fmpz_mpoly_struct));
            fmpz_mpoly_init(sub[i], ctx);
            fmpz_zero(shift + i);
            fmpz_one(stride + i);
        }
        fmpz_set_ui(stride + 2, 2);
        fmpz_set_ui(stride + 3, 2);
        fmpz_set_ui(stride + 4, 2);
        fmpz_mpoly_deflate(b, b, shift, stride, ctx);

        fmpz_mpoly_set_str_pretty(sub[0], "x", vars, ctx);
        fmpz_mpoly_set_str_pretty(sub[1], "y", vars, ctx);
        fmpz_mpoly_set_str_pretty(sub[2], "y+1", vars, ctx);
        fmpz_mpoly_set_str_pretty(sub[3], "y+2", vars, ctx);
        fmpz_mpoly_set_str_pretty(sub[4], "y+3", vars, ctx);
        fmpz_mpoly_compose_fmpz_mpoly(a, b, sub, ctx, ctx);
        check_omega(1, 1, a, ctx);

        fmpz_mpoly_set_str_pretty(sub[0], "x", vars, ctx);
        fmpz_mpoly_set_str_pretty(sub[1], "y", vars, ctx);
        fmpz_mpoly_set_str_pretty(sub[2], "y", vars, ctx);
        fmpz_mpoly_set_str_pretty(sub[3], "y+1", vars, ctx);
        fmpz_mpoly_set_str_pretty(sub[4], "y+2", vars, ctx);
        fmpz_mpoly_compose_fmpz_mpoly(a, b, sub, ctx, ctx);
        check_omega(1, 1, a, ctx);

        fmpz_mpoly_set_str_pretty(sub[0], "x + y + 2", vars, ctx);
        fmpz_mpoly_set_str_pretty(sub[1], "y", vars, ctx);
        fmpz_mpoly_set_str_pretty(sub[2], "y+1", vars, ctx);
        fmpz_mpoly_set_str_pretty(sub[3], "y+4", vars, ctx);
        fmpz_mpoly_set_str_pretty(sub[4], "y+9", vars, ctx);
        fmpz_mpoly_compose_fmpz_mpoly(a, b, sub, ctx, ctx);
        check_omega(1, 1, a, ctx);

        fmpz_mpoly_factor_clear(fac, ctx);
        fmpz_mpoly_clear(a, ctx);
        fmpz_mpoly_ctx_clear(ctx);
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}
