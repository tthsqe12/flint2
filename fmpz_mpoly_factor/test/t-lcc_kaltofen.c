/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mpoly_factor.h"


int
main(void)
{
    slong i;
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

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}
