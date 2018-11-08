/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include <stdlib.h>
#include "nmod_mpoly.h"


int
main(void)
{
    FLINT_TEST_INIT(state);

    flint_printf("divides_heap_threaded....\n");
    fflush(stdout);

    {
        nmod_mpoly_ctx_t ctx;
        nmod_mpoly_t p, f, g, h;
        const char * vars[] = {"x","y","z","t","u"};

        nmod_mpoly_ctx_init(ctx, 5, ORD_DEGLEX, 101);
        nmod_mpoly_init(f, ctx);
        nmod_mpoly_init(g, ctx);
        nmod_mpoly_init(p, ctx);
        nmod_mpoly_init(h, ctx);
        nmod_mpoly_set_str_pretty(f, "(1+x+y+2*z^2+3*t^3+5*u^5)^2", vars, ctx);
        nmod_mpoly_set_str_pretty(g, "(1+u+t+2*z^2+3*y^3+5*x^5)^2", vars, ctx);
        nmod_mpoly_mul(p, f, g, ctx);

flint_printf("f: "); nmod_mpoly_print_pretty(f, vars, ctx); flint_printf("\n");
flint_printf("g: "); nmod_mpoly_print_pretty(g, vars, ctx); flint_printf("\n");
flint_printf("p: "); nmod_mpoly_print_pretty(p, vars, ctx); flint_printf("\n");


        nmod_mpoly_divides_heap_threaded(h, p, f, ctx);
flint_printf("h: "); nmod_mpoly_print_pretty(h, vars, ctx); flint_printf("\n");


        nmod_mpoly_clear(f, ctx);
        nmod_mpoly_clear(g, ctx);
        nmod_mpoly_clear(p, ctx);
        nmod_mpoly_clear(h, ctx);
        nmod_mpoly_ctx_clear(ctx);
    }

    FLINT_TEST_CLEANUP(state);

    printf("PASS\n");
    return 0;
}

