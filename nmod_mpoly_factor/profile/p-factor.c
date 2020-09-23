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
#include "nmod_mpoly_factor.h"
#include "profiler.h"


slong check_omega(slong om, const nmod_mpoly_t p, const nmod_mpoly_ctx_t ctx)
{
    slong i;
    nmod_mpoly_factor_t g;
    fmpz_t omega;
    timeit_t timer;

    fmpz_init(omega);
    nmod_mpoly_factor_init(g, ctx);

    timeit_start(timer);
    nmod_mpoly_factor_wang(g, p, ctx);
    timeit_stop(timer);

    fmpz_zero(omega);
    for (i = 0; i < g->num; i++)
        fmpz_add(omega, omega, g->exp + i);

    if (fmpz_cmp_si(omega, om) < 0)
    {
        flint_printf("factorization has wrong number of factors\n");
        flint_abort();        
    }

    nmod_mpoly_factor_clear(g, ctx);
    fmpz_clear(omega);

    return timer->wall;
}


int main(int argc, char *argv[])
{
    slong i, time, total_time = 0;

    {
        nmod_mpoly_ctx_t ctx;
        nmod_mpoly_t a, b;
        const char * vars[] = {"x", "y", "z", "t"};

        nmod_mpoly_ctx_init(ctx, 4, ORD_LEX, 104891);
        nmod_mpoly_init(a, ctx);
        nmod_mpoly_init(b, ctx);

        for (i = 1; i <= 20; i++)
        {
            nmod_mpoly_set_str_pretty(b, "1+x+y+z+t", vars, ctx);
            nmod_mpoly_pow_ui(b, b, i, ctx);
            nmod_mpoly_add_ui(a, b, 3, ctx);
            nmod_mpoly_add_ui(b, b, 7, ctx);
            nmod_mpoly_mul(a, a, b, ctx);

            time = check_omega(2, a, ctx);
            flint_printf("power %wd: %wd\n", i, time);
            total_time += time;
        }

        nmod_mpoly_clear(a, ctx);
        nmod_mpoly_clear(b, ctx);
        nmod_mpoly_ctx_clear(ctx);
    }

    flint_printf("total_time: %wd\n", total_time);

    flint_cleanup_master();
    return 0;
}

