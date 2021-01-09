/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mod_mpoly.h"

void fmpz_mod_poly_multi_crt_init(
    fmpz_mod_poly_multi_crt_t P,
    const fmpz_mod_ctx_t ctx)
{
    P->prog = NULL;
    P->moduli = NULL;
    P->invmoduli = NULL;
    P->fracmoduli = NULL;
    P->alloc = 0;
    P->length = 0;
    P->localsize = 1;
    P->temp1loc = 0;
    P->temp2loc = 0;
    P->temp3loc = 0;
    P->temp4loc = 0;
    P->good = 0;
}

static void _fmpz_mod_poly_multi_crt_fit_length(
    fmpz_mod_poly_multi_crt_t P,
    slong k,
    const fmpz_mod_ctx_t ctx)
{
    slong i;

    k = FLINT_MAX(WORD(1), k);

    for (i = k; i < P->alloc; i++)
    {
        fmpz_mod_poly_clear(P->prog[i].b_modulus, ctx);
        fmpz_mod_poly_clear(P->prog[i].c_modulus, ctx);
        fmpz_mod_poly_clear(P->moduli + i, ctx);
        fmpz_mod_poly_clear(P->invmoduli + i, ctx);
        fmpz_mod_poly_clear(P->fracmoduli + i, ctx);
    }

    P->prog = FLINT_ARRAY_REALLOC(P->prog, k, _fmpz_mod_poly_multi_crt_prog_instr);
    P->moduli = FLINT_ARRAY_REALLOC(P->moduli, k, fmpz_mod_poly_struct);
    P->invmoduli = FLINT_ARRAY_REALLOC(P->invmoduli, k, fmpz_mod_poly_struct);
    P->fracmoduli = FLINT_ARRAY_REALLOC(P->fracmoduli, k, fmpz_mod_poly_struct);

    for (i = P->alloc; i < k; i++)
    {
        fmpz_mod_poly_init(P->prog[i].b_modulus, ctx);
        fmpz_mod_poly_init(P->prog[i].c_modulus, ctx);
        fmpz_mod_poly_init(P->moduli + i, ctx);
        fmpz_mod_poly_init(P->invmoduli + i, ctx);
        fmpz_mod_poly_init(P->fracmoduli + i, ctx);
    }

    P->alloc = k;
}


void fmpz_mod_poly_multi_crt_clear(
    fmpz_mod_poly_multi_crt_t P,
    const fmpz_mod_ctx_t ctx)
{
    slong i;

    for (i = 0; i < P->alloc; i++)
    {
        fmpz_mod_poly_clear(P->prog[i].b_modulus, ctx);
        fmpz_mod_poly_clear(P->prog[i].c_modulus, ctx);
        fmpz_mod_poly_clear(P->moduli + i, ctx);
        fmpz_mod_poly_clear(P->invmoduli + i, ctx);
        fmpz_mod_poly_clear(P->fracmoduli + i, ctx);
    }

    flint_free(P->prog);
    flint_free(P->moduli);
    flint_free(P->invmoduli);
    flint_free(P->fracmoduli);
}

static int _fill_pfrac(
    slong * link,
    const fmpz_mod_poly_struct * v,
    fmpz_mod_poly_struct * w,
    slong j,
    const fmpz_mod_poly_t W,
    const fmpz_mod_ctx_t ctx,
    fmpz_mod_poly_t g, /* temps */
    fmpz_mod_poly_t s,
    fmpz_mod_poly_t t)
{
    if (j < 0)
        return 1;

    /* W/(v[j]*v[j+1]) = w[j]/v[j] + w[j+1]/v[j+1] */

    if (fmpz_mod_poly_degree(v + j, ctx) < 1)
        return 0;
    if (fmpz_mod_poly_degree(v + j + 1, ctx) < 1)
        return 0;

    fmpz_mod_poly_xgcd(g, s, t, v + j, v + j + 1, ctx);
    if (!fmpz_mod_poly_is_one(g, ctx))
        return 0;

    fmpz_mod_poly_mulmod(w + j, W, t, v + j, ctx);
    fmpz_mod_poly_mulmod(w + j + 1, W, s, v + j + 1, ctx);

#if FLINT_WANT_ASSERT
    fmpz_mod_poly_mul(s, w + j, v + j + 1, ctx);
    fmpz_mod_poly_mul(t, w + j + 1, v + j, ctx);
    fmpz_mod_poly_add(s, s, t, ctx);
    FLINT_ASSERT(fmpz_mod_poly_equal(s, W, ctx));
#endif

    return _fill_pfrac(link, v, w, link[j], w + j, ctx, g, s, t) &&
           _fill_pfrac(link, v, w, link[j + 1], w + j + 1, ctx, g, s, t);
}


static void _invert(
    fmpz_mod_poly_t ainv,
    const fmpz_mod_poly_t a,
    const fmpz_mod_ctx_t ctx)
{
    fmpz_mod_poly_reverse(ainv, a, a->length, ctx);
    fmpz_mod_poly_inv_series(ainv, ainv, a->length, ctx);
}


static void _fill_prog(
    fmpz_mod_poly_multi_crt_t P,
    slong * link,
    fmpz_mod_poly_struct * v,
    fmpz_mod_poly_struct * w,
    slong j,
    slong ret_idx,
    const fmpz_mod_ctx_t ctx)
{
    slong i, b_idx, c_idx;
    slong next_ret_idx = ret_idx;

    FLINT_ASSERT(j >= 0);

    if (link[j] >= 0)
    {
        b_idx = ++next_ret_idx;
        _fill_prog(P, link, v, w, link[j], b_idx, ctx);
    }
    else
    {
        b_idx = -1 - link[j];
        FLINT_ASSERT(b_idx < P->alloc);
        fmpz_mod_poly_set(P->moduli + b_idx, v + j, ctx);
        fmpz_mod_poly_set(P->fracmoduli + b_idx, w + j, ctx);
        _invert(P->invmoduli + b_idx, P->moduli + b_idx, ctx);
        b_idx = -1 - b_idx;
    }

    if (link[j + 1] >= 0)
    {
        c_idx = ++next_ret_idx;
        _fill_prog(P, link, v, w, link[j + 1], c_idx, ctx);
    }
    else
    {
        c_idx = -1 - link[j + 1];
        FLINT_ASSERT(c_idx < P->alloc);
        fmpz_mod_poly_set(P->moduli + c_idx, v + j + 1, ctx);
        fmpz_mod_poly_set(P->fracmoduli + c_idx, w + j + 1, ctx);
        _invert(P->invmoduli + c_idx, P->moduli + c_idx, ctx);
        c_idx = -1 - c_idx;
    }

    i = P->length;
    FLINT_ASSERT(i < P->alloc);
    P->prog[i].a_idx = ret_idx;
    P->prog[i].b_idx = b_idx;
    P->prog[i].c_idx = c_idx;
    fmpz_mod_poly_set(P->prog[i].b_modulus, v + j, ctx);
    fmpz_mod_poly_set(P->prog[i].c_modulus, v + j + 1, ctx);
    P->length = i + 1;

    P->localsize = FLINT_MAX(P->localsize, 1 + next_ret_idx);
}


int fmpz_mod_poly_multi_crt_precompute(
    fmpz_mod_poly_multi_crt_t P,
    const fmpz_mod_poly_struct * f,
    slong r,
    const fmpz_mod_ctx_t ctx)
{
    slong i, j;
    slong * link;
    fmpz_mod_poly_struct * v;
    fmpz_mod_poly_struct * w;
    fmpz_mod_poly_t one, g, s, t;

    FLINT_ASSERT(r > 0);

    _fmpz_mod_poly_multi_crt_fit_length(P, r, ctx);
    P->length = 0;
    P->localsize = 1;

    if (r < 2)
    {
        /*
            There is only one modulus. Let's compute as
                output[0] = ((input[0]*1) mod f[0])*1
                           +(input[0]*1) mod f[0])*0
        */

        P->good = !fmpz_mod_poly_is_zero(f + 0, ctx);
        
        if (P->good)
        {
            fmpz_mod_poly_set(P->moduli + 0, f + 0, ctx);
            fmpz_mod_poly_one(P->fracmoduli + 0, ctx);
            _invert(P->invmoduli + 0, P->moduli + 0, ctx);

            i = 0;
            fmpz_mod_poly_one(P->prog[i].b_modulus, ctx);
            fmpz_mod_poly_zero(P->prog[i].c_modulus, ctx);
            P->prog[i].a_idx = 0;
            P->prog[i].b_idx = -WORD(1);
            P->prog[i].c_idx = -WORD(1);
            P->length = i + 1;

            P->temp1loc = P->localsize++;
            P->temp2loc = P->localsize++;
            P->temp3loc = P->localsize++;
            P->temp4loc = P->localsize++;
        }

        return P->good;
    }

    fmpz_mod_poly_init(one, ctx);
    fmpz_mod_poly_init(g, ctx);
    fmpz_mod_poly_init(s, ctx);
    fmpz_mod_poly_init(t, ctx);

    link = FLINT_ARRAY_ALLOC(2*r - 2, slong);
    v = FLINT_ARRAY_ALLOC(2*(2*r - 2), fmpz_mod_poly_struct);
    w = v + 2*r - 2;

    for (i = 0; i < 2*(2*r - 2); i++)
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

    fmpz_mod_poly_one(one, ctx);
    P->good = _fill_pfrac(link, v, w, 2*r - 4, one, ctx, g, s, t);

    if (P->good)
    {
        _fill_prog(P, link, v, w, 2*r - 4, 0, ctx);
        P->temp1loc = P->localsize++;
        P->temp2loc = P->localsize++;
        P->temp3loc = P->localsize++;
        P->temp4loc = P->localsize++;
    }

    fmpz_mod_poly_clear(one, ctx);
    fmpz_mod_poly_clear(g, ctx);
    fmpz_mod_poly_clear(s, ctx);
    fmpz_mod_poly_clear(t, ctx);

    for (i = 0; i < 2*(2*r - 2); i++)
        fmpz_mod_poly_clear(v + i, ctx);

    flint_free(link);
    flint_free(v);

    return P->good;
}


void fmpz_mod_poly_multi_crt_precomp(
    fmpz_mod_poly_t output,
    const fmpz_mod_poly_multi_crt_t P,
    const fmpz_mod_poly_struct * inputs,
    const fmpz_mod_ctx_t ctx)
{
    slong i;
    fmpz_mod_poly_struct * out;
    TMP_INIT;

    TMP_START;
    out = TMP_ARRAY_ALLOC(P->localsize, fmpz_mod_poly_struct);
    for (i = 0; i < P->localsize; i++)
        fmpz_mod_poly_init(out + i, ctx);

    fmpz_mod_poly_swap(out + 0, output, ctx);
    _fmpz_mod_poly_multi_crt_run(out, P, inputs, ctx);
    fmpz_mod_poly_swap(out + 0, output, ctx);

    for (i = 0; i < P->localsize; i++)
        fmpz_mod_poly_clear(out + i, ctx);

    TMP_END;
}

int fmpz_mod_poly_multi_crt(
    fmpz_mod_poly_t output,
    const fmpz_mod_poly_struct * moduli,
    const fmpz_mod_poly_struct * inputs,
    slong len,
    const fmpz_mod_ctx_t ctx)
{
    int success;
    fmpz_mod_poly_multi_crt_t P;
    fmpz_mod_poly_multi_crt_init(P, ctx);
    success = fmpz_mod_poly_multi_crt_precompute(P, moduli, len, ctx);
    fmpz_mod_poly_multi_crt_precomp(output, P, inputs, ctx);
    fmpz_mod_poly_multi_crt_clear(P, ctx);
    return success;
}


void _fmpz_mod_poly_multi_crt_run(
    fmpz_mod_poly_struct * outputs,
    const fmpz_mod_poly_multi_crt_t P,
    const fmpz_mod_poly_struct * inputs,
    const fmpz_mod_ctx_t ctx)
{
    slong i, a, b, c;
    const fmpz_mod_poly_struct * m = P->moduli;
    const fmpz_mod_poly_struct * mi = P->invmoduli;
    const fmpz_mod_poly_struct * mf = P->fracmoduli;
    fmpz_mod_poly_struct * A, * B, * C, * t1, * t2, * t3, * t4;

    t1 = outputs + P->temp1loc;
    t2 = outputs + P->temp2loc;
    t3 = outputs + P->temp3loc;
    t4 = outputs + P->temp4loc;

    for (i = 0; i < P->length; i++)
    {
        a = P->prog[i].a_idx;
        b = P->prog[i].b_idx;
        c = P->prog[i].c_idx;

        A = outputs + a;
        B = outputs + b;
        C = outputs + c;

        FLINT_ASSERT(a >= 0);

        /* leaves need a premultiplication by (prod_{i!=b} m[i])^-1 mod m[i] */
        if (b < 0)
        {
            b = -b - 1;
            B = t1;

            fmpz_mod_poly_mul(t3, inputs + b, mf + b, ctx);

            if (t3->length <= 2*m[b].length - 2)
                fmpz_mod_poly_divrem_newton_n_preinv(t4, B, t3, m + b, mi + b, ctx);
            else
                fmpz_mod_poly_rem(B, t3, m + b, ctx);
        }

        /* ditto */
        if (c < 0)
        {
            c = -c - 1;
            C = t2;

            fmpz_mod_poly_mul(t3, inputs + c, mf + c, ctx);

            if (t3->length <= 2*m[c].length - 2)
                fmpz_mod_poly_divrem_newton_n_preinv(t4, C, t3, m + c, mi + c, ctx);
            else
                fmpz_mod_poly_rem(C, t3, m + c, ctx);
        }

        /* A = B*c_m + C*b_m */
        fmpz_mod_poly_mul(A, B, P->prog[i].c_modulus, ctx);
        fmpz_mod_poly_mul(t3, C, P->prog[i].b_modulus, ctx);
        fmpz_mod_poly_add(A, A, t3, ctx);

        /* last calculation should write answer to outputs[0] */
        FLINT_ASSERT(i + 1 < P->length || A == outputs + 0);
    }
}
