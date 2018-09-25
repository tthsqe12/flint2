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
/*
    int i, j, result;
*/
    FLINT_TEST_INIT(state);

    flint_printf("evaluate_one....");
    fflush(stdout);

    {
        nmod_mpoly_ctx_t ctx, ctx2;
        nmod_mpoly_t f, g, h;
        ulong res, vals[10];
        nmod_poly_struct * subs[10];
        nmod_mpoly_struct * comp[10];
        nmod_poly_t p;
        const char * vars[] = {"x", "y", "z"};
        const char * vars2[] = {"X", "Y"};

        nmod_mpoly_ctx_init(ctx, 3, ORD_DEGLEX, 100);
        nmod_mpoly_init(f, ctx);
        nmod_mpoly_init(g, ctx);
        nmod_mpoly_set_str_pretty(f, "(1+y^9)*x*z+(2*y+y^1)*z^2", vars, ctx);
        printf("f: "); nmod_mpoly_print_pretty(f, vars, ctx); printf("\n");
        nmod_mpoly_evaluate_one_ui(g, f, 1, 11, ctx);
        nmod_mpoly_assert_canonical(g, ctx);
        printf("g: "); nmod_mpoly_print_pretty(g, vars, ctx); printf("\n");

        vals[0] = 2;
        vals[1] = 5;
        vals[2] = 7;
        res = nmod_mpoly_evaluate_all_ui(f, vals, ctx);
        flint_printf("res: %wu\n", res);
/*
        nmod_mpoly_fit_bits(f, 128, ctx);
*/
        subs[0] = (nmod_poly_struct*) flint_malloc(sizeof(nmod_poly_struct));
        subs[1] = (nmod_poly_struct*) flint_malloc(sizeof(nmod_poly_struct));
        subs[2] = (nmod_poly_struct*) flint_malloc(sizeof(nmod_poly_struct));
        nmod_poly_init(subs[0], ctx->ffinfo->mod.n);
        nmod_poly_init(subs[1], ctx->ffinfo->mod.n);
        nmod_poly_init(subs[2], ctx->ffinfo->mod.n);
        nmod_poly_init(p, ctx->ffinfo->mod.n);
        nmod_poly_set_coeff_ui(subs[0], 0, 2);
        nmod_poly_set_coeff_ui(subs[1], 0, 5);
        nmod_poly_set_coeff_ui(subs[2], 0, 7);
printf("x: "); nmod_poly_print_pretty(subs[0], "X"); printf("\n");
printf("y: "); nmod_poly_print_pretty(subs[1], "X"); printf("\n");
printf("z: "); nmod_poly_print_pretty(subs[2], "X"); printf("\n");
        nmod_mpoly_compose_nmod_poly(p, f, subs, ctx);
printf("r: "); nmod_poly_print_pretty(p, "X"); printf("\n");
        nmod_poly_clear(p);
        nmod_poly_clear(subs[0]);
        nmod_poly_clear(subs[1]);
        nmod_poly_clear(subs[2]);
        flint_free(subs[0]);
        flint_free(subs[1]);
        flint_free(subs[2]);

        nmod_mpoly_ctx_init(ctx2, 2, ORD_DEGLEX, 100);
        comp[0] = (nmod_mpoly_struct*) flint_malloc(sizeof(nmod_mpoly_struct));
        comp[1] = (nmod_mpoly_struct*) flint_malloc(sizeof(nmod_mpoly_struct));
        comp[2] = (nmod_mpoly_struct*) flint_malloc(sizeof(nmod_mpoly_struct));
        nmod_mpoly_init(comp[0], ctx2);
        nmod_mpoly_init(comp[1], ctx2);
        nmod_mpoly_init(comp[2], ctx2);
        nmod_mpoly_set_str_pretty(comp[0], "X+1*Y", vars2, ctx2);
        nmod_mpoly_set_str_pretty(comp[1], "X+2*Y", vars2, ctx2);
        nmod_mpoly_set_str_pretty(comp[2], "X+3*Y", vars2, ctx2);
printf("f: "); nmod_mpoly_print_pretty(f, vars, ctx); printf("\n");
printf("x: "); nmod_mpoly_print_pretty(comp[0], vars2, ctx2); printf("\n");
printf("y: "); nmod_mpoly_print_pretty(comp[1], vars2, ctx2); printf("\n");
printf("z: "); nmod_mpoly_print_pretty(comp[2], vars2, ctx2); printf("\n");
        nmod_mpoly_init(h, ctx2);
        nmod_mpoly_compose_nmod_mpoly(h, f, comp, ctx, ctx2);
printf("h: "); nmod_mpoly_print_pretty(h, vars2, ctx2); printf("\n");


        nmod_mpoly_clear(comp[0], ctx2);
        nmod_mpoly_clear(comp[1], ctx2);
        nmod_mpoly_clear(comp[2], ctx2);
        nmod_mpoly_ctx_clear(ctx2);

        nmod_mpoly_clear(f, ctx);
        nmod_mpoly_clear(g, ctx);
        nmod_mpoly_ctx_clear(ctx);
    }

    FLINT_TEST_CLEANUP(state);
    flint_printf("PASS\n");
    return 0;
}

