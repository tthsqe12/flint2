/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fmpz_mpoly.h"

int
main(void)
{
    slong i, j;
    FLINT_TEST_INIT(state);

    flint_printf("gcd_zippel....\n");
    fflush(stdout);

    {
        fmpz_mpoly_ctx_t ctx;
        fmpz_mpoly_t a, b, g;
        const char* vars[] = {"X","x1","x0"};

        fmpz_mpoly_ctx_init(ctx, 3, ORD_DEGREVLEX);
        fmpz_mpoly_init(a, ctx);
        fmpz_mpoly_init(b, ctx);
        fmpz_mpoly_init(g, ctx);
        fmpz_mpoly_set_str_pretty(a, "-6*(1+x0)*(1+x1)*((x1+x0^2)*X^2-3*x0^5*x1^6+x1+20000000000*x0)*(x0^3 + X + 2*x1 + 3*x0)", vars, ctx);
        fmpz_mpoly_set_str_pretty(b, "-4*(1+x0)*(1+x1)*((x1+x0^2)*X^2-3*x0^5*x1^6+x1+20000000000*x0)*(x0^2 + X + 4*x1 + 5*x0)", vars, ctx);

        fmpz_mpoly_gcd_zippel(g, a, b, ctx);
        fmpz_mpoly_assert_canonical(g, ctx);

printf("g: "); fmpz_mpoly_print_pretty(g, vars, ctx); printf("\n");

        fmpz_mpoly_gcd_brown(g, a, b, ctx);
        fmpz_mpoly_assert_canonical(g, ctx);

printf("g: "); fmpz_mpoly_print_pretty(g, vars, ctx); printf("\n");


        fmpz_mpoly_clear(g, ctx);
        fmpz_mpoly_clear(a, ctx);
        fmpz_mpoly_clear(b, ctx);
        fmpz_mpoly_ctx_clear(ctx);
    }


    printf("PASS\n");
    FLINT_TEST_CLEANUP(state);

    return 0;
}

