/*
    Copyright 2019 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include <stdlib.h>
#include "profiler.h"
#include "fmpz_mpoly_factor.h"


slong check_omega(slong om, const fmpz_mpoly_t p, const fmpz_mpoly_ctx_t ctx)
{
    slong i;
    fmpz_mpoly_factor_t g;
    fmpz_t omega;
    timeit_t timer;

    fmpz_init(omega);
    fmpz_mpoly_factor_init(g, ctx);

    timeit_start(timer);
    fmpz_mpoly_factor_wang(g, p, ctx);
    timeit_stop(timer);

    fmpz_zero(omega);
    for (i = 0; i < g->num; i++)
        fmpz_add(omega, omega, g->exp + i);

    if (fmpz_cmp_si(omega, om) != 0)
    {
        flint_printf("factorization has wrong number of factors\n");
        flint_abort();        
    }

    fmpz_mpoly_factor_clear(g, ctx);
    fmpz_clear(omega);

    return timer->wall;
}


int main(int argc, char *argv[])
{
    slong i, j, k;
    slong time, total_time = 0;

    {
        fmpz_mpoly_ctx_t ctx;
        fmpz_mpoly_t a, b, c;
        const char * vars[] = {"x", "y", "z", "t"};

        fmpz_mpoly_ctx_init(ctx, 4, ORD_LEX);
        fmpz_mpoly_init(a, ctx);
        fmpz_mpoly_init(b, ctx);
        fmpz_mpoly_init(c, ctx);

        for (i = 0; i <= 20; i++)
        {
            fmpz_mpoly_set_str_pretty(a, "x", vars, ctx);
            fmpz_mpoly_set_str_pretty(b, "y", vars, ctx);
            fmpz_mpoly_set_str_pretty(c, "1+x+y+z+t", vars, ctx);
            fmpz_mpoly_pow_ui(c, c, i, ctx);
            fmpz_mpoly_add(a, a, c, ctx);
            fmpz_mpoly_add(b, b, c, ctx);
            fmpz_mpoly_mul(a, a, b, ctx);

            k = (i > 0);
            for (j = 1; j <= i; j++)
                if ((j%2) != 0 && (i%j) == 0)
                    k++;
            k = 2;

            time = check_omega(k, a, ctx);
            flint_printf("%wd: %wd\n", i, time);
            total_time += time;
        }

        fmpz_mpoly_clear(a, ctx);
        fmpz_mpoly_clear(b, ctx);
        fmpz_mpoly_clear(c, ctx);
        fmpz_mpoly_ctx_clear(ctx);
    }

    flint_printf("total_time: %wd\n", total_time);

    flint_cleanup_master();
    return 0;
}
