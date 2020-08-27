/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mpoly_factor.h"

#define FLINT_ARRAY_ALLOC(n, T) (T *) flint_malloc((n)*sizeof(T))


int
main(void)
{
    slong i, j, k, l;
    FLINT_TEST_INIT(state);

    flint_printf("lcc_kaltofen....");
    fflush(stdout);

    {
        fmpz_mpoly_ctx_t ctx;
        fmpz_mpoly_struct divs[10];
        fmpz_poly_struct Auf[10];
        fmpz_mpoly_factor_t Af;
        fmpz_mpoly_t A;
        fmpz alphas[10];


        fmpz_mpoly_ctx_init(ctx, 3, ORD_LEX);
        fmpz_init_set_ui(alphas + 0, 1);
        fmpz_init_set_ui(alphas + 1, 1);
        fmpz_init_set_ui(alphas + 2, 1);
        fmpz_init_set_ui(alphas + 3, 1);
        fmpz_mpoly_factor_init(Af, ctx);
        fmpz_mpoly_init(A, ctx);
        for (i = 0; i < 10; i++)
        {
            fmpz_mpoly_init(divs + i, ctx);
            fmpz_mpoly_one(divs + i, ctx);
            fmpz_poly_init(Auf + i);
        }


        fmpz_mpoly_set_str_pretty(A, "(x2+1)", NULL, ctx);
        fmpz_mpoly_get_fmpz_poly(Auf + 0, A, 1, ctx);
        fmpz_mpoly_set_str_pretty(A, "(x2+2)^2", NULL, ctx);
        fmpz_mpoly_get_fmpz_poly(Auf + 1, A, 1, ctx);
        fmpz_mpoly_set_str_pretty(A, "(x2+1)^3", NULL, ctx);
        fmpz_mpoly_get_fmpz_poly(Auf + 2, A, 1, ctx);

        fmpz_mpoly_set_str_pretty(A, "(x2+x3)^4*(x2+x3+1)^2*x3^3", NULL, ctx);
        fmpz_mpoly_factor_fit_length(Af, 3, ctx);
        Af->num = 3;
        fmpz_set_ui(Af->exp + 0, 4); fmpz_mpoly_set_str_pretty(Af->poly + 0, "x2+x3", NULL, ctx);
        fmpz_set_ui(Af->exp + 1, 2); fmpz_mpoly_set_str_pretty(Af->poly + 1, "x2+x3+1", NULL, ctx);
        fmpz_set_ui(Af->exp + 2, 3); fmpz_mpoly_set_str_pretty(Af->poly + 2, "x3", NULL, ctx);

        fmpz_mpoly_factor_lcc_kaltofen(divs, 3, Af, Auf, 1, alphas, ctx);
    }

    for (i = 0; i < 10*flint_test_multiplier(); i++)
    {
        slong r, v, nvars;
        ulong * bounds;
        fmpz_mpoly_ctx_t ctx;
        fmpz_mpoly_struct A[1], t1[1], t2[1], * divs, * lcs;
        fmpz_mpoly_factor_t Af;
        fmpz_poly_struct ut1[1], ut2[1], * ulcs;
        fmpz * alphas;
        

flint_printf("--------------------------------\n------------ i = %wd -------------\n", i);

        nvars = 3 + n_randint(state, 3);

        fmpz_mpoly_ctx_init(ctx, nvars, ORD_LEX);
        fmpz_mpoly_init(A, ctx);
        fmpz_mpoly_init(t1, ctx);
        fmpz_mpoly_init(t2, ctx);
        fmpz_mpoly_factor_init(Af, ctx);
        fmpz_poly_init(ut1);
        fmpz_poly_init(ut2);

        bounds = FLINT_ARRAY_ALLOC(nvars, ulong);
        alphas = _fmpz_vec_init(nvars - 1);

        r = 2 + n_randint(state, 5);
        divs = FLINT_ARRAY_ALLOC(r, fmpz_mpoly_struct);
        lcs = FLINT_ARRAY_ALLOC(r, fmpz_mpoly_struct);
        ulcs = FLINT_ARRAY_ALLOC(r, fmpz_poly_struct);
        for (j = 0; j < r; j++)
        {
            fmpz_mpoly_init(divs + j, ctx);
            fmpz_mpoly_one(divs + j, ctx);
            fmpz_mpoly_init(lcs + j, ctx);
            fmpz_mpoly_set_ui(lcs + j, 1 + n_randint(state, 100), ctx);
            fmpz_poly_init(ulcs + j);
        }

        for (j = 0; j < r; j++)
        {
            bounds[0] = 1;
            for (k = 1; k < nvars; k++)
                bounds[k] = 1 + n_randint(state, 4);
            fmpz_mpoly_randtest_bounds(t1, state, 4, 20, bounds, ctx);
            for (k = n_randint(state, 3); k >= 0; k--)
            {
                l = n_randint(state, r);
                fmpz_mpoly_mul(lcs + l, lcs + l, t1, ctx);
            }
        }

        fmpz_mpoly_one(A, ctx);
        for (j = 0; j < r; j++)
            fmpz_mpoly_mul(A, A, lcs + j, ctx);
/*
flint_printf("A: ");
fmpz_mpoly_print_pretty(A, NULL, ctx);
flint_printf("\n");
*/
        if (!fmpz_mpoly_factor_squarefree(Af, A, ctx))
        {
            flint_printf("FAIL:\ncheck factor_squarefree could be computed\n");
            flint_abort();
        }

flint_printf("Af: ");
fmpz_mpoly_factor_print_pretty(Af, NULL, ctx);
flint_printf("\n");

        for (j = 0; j < nvars - 1; j++)
        {
            fmpz_set_ui(alphas + j, n_urandint(state, 100));
            if (n_randint(state, 2))
                fmpz_neg(alphas + j, alphas + j);
flint_printf("alphas[%wd]: ", j); fmpz_print(alphas + j); flint_printf("\n");
        }

        for (v = 1; v < nvars; v++)
        {
            int have_zero = 0;

            for (j = 0; j < r; j++)
            {
                /* evaluation of lcs[j] and divs[j] down to univar */

                fmpz_mpoly_set(t1, lcs + j, ctx);
                fmpz_mpoly_set(t2, divs + j, ctx);
                for (k = 1; k < nvars; k++)
                {
                    if (k == v)
                        continue;
                    fmpz_mpoly_evaluate_one_fmpz(t1, t1, k, alphas + k - 1, ctx);
                    fmpz_mpoly_evaluate_one_fmpz(t2, t2, k, alphas + k - 1, ctx);
                }

                if (fmpz_mpoly_is_zero(t1, ctx) || fmpz_mpoly_is_zero(t2, ctx))
                {
                    have_zero = 1;
                }
                else
                {
                    FLINT_ASSERT(fmpz_mpoly_is_fmpz_poly(t1, v, ctx));
                    FLINT_ASSERT(fmpz_mpoly_is_fmpz_poly(t2, v, ctx));
                    fmpz_mpoly_get_fmpz_poly(ut1, t1, v, ctx);
                    fmpz_mpoly_get_fmpz_poly(ut2, t2, v, ctx);
                    fmpz_poly_primitive_part(ut2, ut2);
                    if (!fmpz_poly_divides(ulcs + j, ut1, ut2))
                    {
                        flint_printf("FAIL:\nbad divisor\n");
                        flint_abort();
                    }
                }
            }

            if (have_zero)
            {
flint_printf("had zero at v = %wd\n", v);
                continue;
            }

            fmpz_mpoly_factor_lcc_kaltofen(divs, r, Af, ulcs, v, alphas, ctx);
        }

flint_printf("Af: ");
fmpz_mpoly_factor_print_pretty(Af, NULL, ctx);
flint_printf("\n");

if (Af->num > 0)
{
    flint_printf("oops\n");
usleep(2000000);
}



        flint_free(bounds);
        _fmpz_vec_clear(alphas, nvars - 1);

        for (j = 0; j < r; j++)
        {
            fmpz_mpoly_clear(divs + j, ctx);
            fmpz_mpoly_clear(lcs + j, ctx);
            fmpz_poly_clear(ulcs + j);
        }
        flint_free(divs);
        flint_free(lcs);
        flint_free(ulcs);

        fmpz_mpoly_clear(A, ctx);
        fmpz_mpoly_clear(t1, ctx);
        fmpz_mpoly_clear(t2, ctx);
        fmpz_mpoly_factor_clear(Af, ctx);
        fmpz_poly_clear(ut1);
        fmpz_poly_clear(ut2);
        fmpz_mpoly_ctx_clear(ctx);


    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}
