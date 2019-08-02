/*
    Copyright (C) 2019 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fmpz_mpoly_factor.h"

int
main(void)
{
    FLINT_TEST_INIT(state);

    flint_printf("factor....");
    fflush(stdout);

    {
        fmpz_mpoly_ctx_t ctx;
        fmpz_mpoly_t a;
        fmpz_mpoly_factor_t fac;
        const char * vars[] = {"x", "y", "z"};

        fmpz_mpoly_ctx_init(ctx, 3, ORD_DEGLEX);
        fmpz_mpoly_init(a, ctx);
        fmpz_mpoly_set_str_pretty(a, "x^9-y^9", vars, ctx);

        printf("a: "); fmpz_mpoly_print_pretty(a, vars, ctx); printf("\n");

        fmpz_mpoly_factor_init(fac, ctx);
        fmpz_mpoly_factor(fac, a, 1, ctx);

        printf("fac:\n");
        fmpz_mpoly_factor_print(fac, vars, ctx);

        fmpz_mpoly_factor_clear(fac, ctx);
        fmpz_mpoly_clear(a, ctx);
        fmpz_mpoly_ctx_clear(ctx);
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}
