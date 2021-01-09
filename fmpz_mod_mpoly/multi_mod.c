/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mod_mpoly.h"

void fmpz_mod_poly_multi_mod_init(
    fmpz_mod_poly_multi_mod_t P,
    const fmpz_mod_ctx_t ctx)
{
    P->prog = NULL;
    P->alloc = 0;
    P->length = 0;
    P->localsize = 1;
    P->temp1loc = 0;
}

static void _fmpz_mod_poly_multi_mod_fit_length(
    fmpz_mod_poly_multi_mod_t P,
    slong k,
    const fmpz_mod_ctx_t ctx)
{
    slong i;

    k = FLINT_MAX(WORD(1), k);

    for (i = k; i < P->alloc; i++)
    {
        fmpz_mod_poly_clear(P->prog[i].modulus, ctx);
    }

    P->prog = FLINT_ARRAY_REALLOC(P->prog, k, _fmpz_mod_poly_multi_mod_prog_instr);

    for (i = P->alloc; i < k; i++)
    {
        fmpz_mod_poly_init(P->prog[i].modulus, ctx);
    }

    P->alloc = k;
}


void fmpz_mod_poly_multi_mod_clear(
    fmpz_mod_poly_multi_mod_t P,
    const fmpz_mod_ctx_t ctx)
{
    slong i;

    for (i = 0; i < P->alloc; i++)
    {
        fmpz_mod_poly_clear(P->prog[i].modulus, ctx);
    }

    flint_free(P->prog);
}


static void _fill_prog(
    fmpz_mod_poly_multi_mod_t P,
    slong * link,
    fmpz_mod_poly_struct * v,
    slong j,
    slong input_idx,    /* < 0 for original input, >= 0 for tmp */
    const fmpz_mod_ctx_t ctx)
{
    slong i, k;

    FLINT_ASSERT(j >= 0);

    P->localsize = FLINT_MAX(P->localsize, 1 + input_idx);

    for (k = 0; k < 2; k++)
    {
        i = P->length;
        FLINT_ASSERT(i < P->alloc);

        P->prog[i].b_idx = input_idx;
        fmpz_mod_poly_set(P->prog[i].modulus, v + j + k, ctx);

        if (link[j + k] >= 0)
        {
            P->prog[i].a_idx = input_idx + 1;
            P->length = i + 1;
            _fill_prog(P, link, v, link[j + k], input_idx + 1, ctx);
        }
        else
        {
            P->prog[i].a_idx = link[j + k];
            P->length = i + 1;
        }
    }
}


void fmpz_mod_poly_multi_mod_precompute(
    fmpz_mod_poly_multi_mod_t P,
    const fmpz_mod_poly_struct * f,
    slong r,
    const fmpz_mod_ctx_t ctx)
{
    slong i, j;
    slong * link;
    fmpz_mod_poly_struct * v;

    FLINT_ASSERT(r > 0);

    _fmpz_mod_poly_multi_mod_fit_length(P, 2*r - 1, ctx);
    P->length = 0;
    P->localsize = 1;

    if (r < 2)
    {
        /*
            There is only one modulus. Let's compute as
                output[0] = ((input[0]*1) mod f[0])*1
                           +(input[0]*1) mod f[0])*0
        */

        i = 0;
        fmpz_mod_poly_set(P->prog[i].modulus, f + 0, ctx);
        P->prog[i].a_idx = -1;
        P->prog[i].b_idx = -1;
        P->length = i + 1;

        P->temp1loc = P->localsize++;
        return;
    }

    link = FLINT_ARRAY_ALLOC(2*r - 2, slong);
    v = FLINT_ARRAY_ALLOC(2*r - 2, fmpz_mod_poly_struct);

    for (i = 0; i < 2*r - 2; i++)
        fmpz_mod_poly_init(v + i, ctx);

    for (i = 0; i < r; i++)
    {
        fmpz_mod_poly_set(v + i, f + i, ctx);
        link[i] = -i - 1;
    }

    for (i = r, j = 0; j < 2*r - 4; i++, j += 2)
    {
        slong s, minp, mind;

        minp = j;
        mind = fmpz_mod_poly_degree(v + j, ctx);
        for (s = j + 1; s < i; s++)
        {
            if (fmpz_mod_poly_degree(v + s, ctx) < mind)
            {
                mind = fmpz_mod_poly_degree(v + s, ctx);
                minp = s;
            }
        }
        fmpz_mod_poly_swap(v + j, v + minp, ctx);
        SLONG_SWAP(link[j], link[minp]);

        minp = j + 1;
        mind = fmpz_mod_poly_degree(v + j + 1, ctx);
        for (s = j + 2; s < i; s++)
        {
            if (fmpz_mod_poly_degree(v + s, ctx) < mind)
            {
                mind = fmpz_mod_poly_degree(v + s, ctx);
                minp = s;
            }
        }
        fmpz_mod_poly_swap(v + j + 1, v + minp, ctx);
        SLONG_SWAP(link[j + 1], link[minp]);

        fmpz_mod_poly_mul(v + i, v + j, v + j + 1, ctx);
        link[i] = j;
    }

    _fill_prog(P, link, v, 2*r - 4, -1, ctx);
    P->temp1loc = P->localsize++;

    for (i = 0; i < 2*r - 2; i++)
        fmpz_mod_poly_clear(v + i, ctx);

    flint_free(link);
    flint_free(v);
}


void fmpz_mod_poly_multi_mod_precomp(
    fmpz_mod_poly_struct * outputs,
    const fmpz_mod_poly_multi_mod_t P,
    const fmpz_mod_poly_t input,
    const fmpz_mod_ctx_t ctx)
{
    slong i;
    fmpz_mod_poly_struct * tmps;
    TMP_INIT;

    TMP_START;
    tmps = TMP_ARRAY_ALLOC(P->localsize, fmpz_mod_poly_struct);
    for (i = 0; i < P->localsize; i++)
        fmpz_mod_poly_init(tmps + i, ctx);

    _fmpz_mod_poly_multi_mod_run(outputs, P, input, tmps, ctx);

    for (i = 0; i < P->localsize; i++)
        fmpz_mod_poly_clear(tmps + i, ctx);

    TMP_END;
}

void fmpz_mod_poly_multi_mod(
    fmpz_mod_poly_struct * outputs,
    const fmpz_mod_poly_struct * moduli,
    const fmpz_mod_poly_t input,
    slong len,
    const fmpz_mod_ctx_t ctx)
{
    fmpz_mod_poly_multi_mod_t P;
    fmpz_mod_poly_multi_mod_init(P, ctx);
    fmpz_mod_poly_multi_mod_precompute(P, moduli, len, ctx);
    fmpz_mod_poly_multi_mod_precomp(outputs, P, input, ctx);
    fmpz_mod_poly_multi_mod_clear(P, ctx);
}

void _fmpz_mod_poly_multi_mod_run(
    fmpz_mod_poly_struct * outputs,
    const fmpz_mod_poly_multi_mod_t P,
    const fmpz_mod_poly_t input,
    fmpz_mod_poly_struct * tmps,
    const fmpz_mod_ctx_t ctx)
{
    slong i, a, b;
    fmpz_mod_poly_struct * A, * t1;
    const fmpz_mod_poly_struct * B;

    t1 = tmps + P->temp1loc;

    for (i = 0; i < P->length; i++)
    {
        a = P->prog[i].a_idx;
        b = P->prog[i].b_idx;

        if (a < 0)
            A = outputs + (-1 - a);
        else
            A = tmps + a;

        if (b < 0)
            B = input;
        else
            B = tmps + b;

        fmpz_mod_poly_divrem(t1, A, B, P->prog[i].modulus, ctx);
    }
}
