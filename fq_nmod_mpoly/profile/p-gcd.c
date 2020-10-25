/*
    Copyright 2019 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include <stdlib.h>
#include "fq_nmod_mpoly.h"
#include "fq_nmod_mpoly_factor.h"
#include "profiler.h"

int main(int argc, char *argv[])
{
    {
        fq_nmod_mpoly_ctx_t ctx;
        fq_nmod_mpoly_t a, b, t, g;
        timeit_t timer;

        fq_nmod_mpoly_ctx_init_deg(ctx, 7, ORD_LEX, UWORD(4611686018427388039), 2);
        fq_nmod_mpoly_init(a, ctx);
        fq_nmod_mpoly_init(b, ctx);
        fq_nmod_mpoly_init(t, ctx);
        fq_nmod_mpoly_init(g, ctx);

        fq_nmod_mpoly_set_str_pretty(a, "(1 + #*x1 + x2 + x3 + #*x4 + #*x5 + 2*x6 + x7)^7 + x1", NULL, ctx);
        fq_nmod_mpoly_set_str_pretty(b, "(1 + x1 + #*x2 + x3 + #*x4 + x5 + 3*x6 + #*x7)^7 + x2", NULL, ctx);
        fq_nmod_mpoly_set_str_pretty(t, "(1 + #*x1 + x2 + #*x3 + #*x4 + x5 + 4*x6 + x7)^7 + x3", NULL, ctx);
        fq_nmod_mpoly_mul(a, a, t, ctx);
        fq_nmod_mpoly_mul(b, b, t, ctx);

flint_printf("a->length: %wd\n", a->length);
flint_printf("b->length: %wd\n", b->length);
flint_printf("t->length: %wd\n", t->length);

        {
            fq_nmod_mpoly_factor_t f;
            fq_nmod_mpoly_factor_init(f, ctx);
            timeit_start(timer);
            fq_nmod_mpoly_factor_wang(f, a, ctx);
            timeit_stop(timer);
            flint_printf("factor time: %wd\n", timer->wall);
            if (f->num != 2)
                flint_printf("factor oops!!!\n");
            fq_nmod_mpoly_factor_clear(f, ctx);            
        }

        timeit_start(timer);
        fq_nmod_mpoly_gcd(g, a, b, ctx);
        timeit_stop(timer);

        flint_printf("time: %wd\n", timer->wall);
        if (g->length != t->length)
            flint_printf("oops!!!\n");

timeit_start(timer);
fq_nmod_mpoly_divides(t, a, g, ctx);
fq_nmod_mpoly_divides(t, b, g, ctx);
timeit_stop(timer);
flint_printf("divides time: %wd\n", timer->wall);


        fq_nmod_mpoly_clear(a, ctx);
        fq_nmod_mpoly_clear(b, ctx);
        fq_nmod_mpoly_clear(t, ctx);
        fq_nmod_mpoly_clear(g, ctx);
        fq_nmod_mpoly_ctx_clear(ctx);
    }

    {
        fq_nmod_mpoly_ctx_t ctx;
        fq_nmod_mpoly_t a, b, t, g;
        const char * vars[] = {"x", "y", "z", "t" ,"u", "v", "w"};
        timeit_t timer;

        fq_nmod_mpoly_ctx_init_deg(ctx, 7, ORD_LEX, UWORD(4611686018427388039), 2);
        fq_nmod_mpoly_init(a, ctx);
        fq_nmod_mpoly_init(b, ctx);
        fq_nmod_mpoly_init(t, ctx);
        fq_nmod_mpoly_init(g, ctx);

        fq_nmod_mpoly_set_str_pretty(a, "(1+#*x+y+z+#*t+u+2*v+w)^7+x", vars, ctx);
        fq_nmod_mpoly_set_str_pretty(b, "(1+#-2*x-#*y+z+t+#*u+v+w)^7+y", vars, ctx);
        fq_nmod_mpoly_set_str_pretty(t, "(#+x+y+#*z-t+u+v+3*#*w)^7+w", vars, ctx);
        fq_nmod_mpoly_mul(a, a, t, ctx);
        fq_nmod_mpoly_mul(b, b, t, ctx);

flint_printf("a->length: %wd\n", a->length);
flint_printf("b->length: %wd\n", b->length);
flint_printf("t->length: %wd\n", t->length);

        timeit_start(timer);
        fq_nmod_mpoly_gcd(g, a, b, ctx);
        timeit_stop(timer);

        flint_printf("time: %wd\n", timer->wall);
        if (g->length != t->length)
            flint_printf("oops!!!\n");

        fq_nmod_mpoly_clear(a, ctx);
        fq_nmod_mpoly_clear(b, ctx);
        fq_nmod_mpoly_clear(t, ctx);
        fq_nmod_mpoly_clear(g, ctx);
        fq_nmod_mpoly_ctx_clear(ctx);
    }

    flint_cleanup_master();
    return 0;
}

