/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mod_mpoly_factor.h"


void _fmpz_mod_poly_div_root(
    fmpz_t E,
    fmpz * Q,
    const fmpz * A,
    slong len,
    const fmpz_t c,
    const fmpz_mod_ctx_t ctx)
{
    fmpz_t r, t, t2;
    slong i;

    if (len < 2)
    {
        fmpz_zero(E);
        return;
    }

    fmpz_init(r);
    fmpz_init(t);
    fmpz_init(t2);

    fmpz_set(t, A + len - 2);
    fmpz_set(Q + len - 2, A + len - 1);
    fmpz_set(r, Q + len - 2);

    for (i = len - 2; i > 0; i--)
    {
        fmpz_mod_mul(t2, r, c, ctx);
        fmpz_mod_add(r, t2, t, ctx);
        fmpz_set(t, A + i-1);
        fmpz_set(Q + i-1, r);
    }

    fmpz_mod_mul(t2, r, c, ctx);
    fmpz_mod_add(E, t2, t, ctx);

    fmpz_clear(r);
    fmpz_clear(t);
    fmpz_clear(t2);
    return;
}

void fmpz_mod_poly_div_root(fmpz_t E, fmpz_mod_poly_t Q, 
             const fmpz_mod_poly_t A, const fmpz_t c, const fmpz_mod_ctx_t ctx)
{
    slong len = A->length;

    if (len == 0)
    {
        fmpz_mod_poly_zero(Q, ctx);
        fmpz_zero(E);
        return;
    }

    if (len == 1)
    {
        fmpz_set(E, A->coeffs + 0);
        fmpz_mod_poly_zero(Q, ctx);
        return;
    }

    if (fmpz_is_zero(c))
    {
        fmpz_set(E, A->coeffs + 0);
        fmpz_mod_poly_shift_right(Q, A, 1, ctx);
        return;
    }

    fmpz_mod_poly_fit_length(Q, len - 1, ctx);
    _fmpz_mod_poly_div_root(E, Q->coeffs, A->coeffs, len, c, ctx);
    _fmpz_mod_poly_set_length(Q, len - 1);
}


/*
    conversion between polynomials in coefficient form and point-value form
    and artihmetic in point-value form
*/

/**************** conversion ************************************************/

/* p = 1 mod 4 */

static slong _find_eval_points4(
    fmpz * list,
    slong d,
    const fmpz_mod_ctx_t ctx)
{
    slong i, len;
    const fmpz * p = fmpz_mod_ctx_modulus(ctx);
    fmpz_t n, halfp, mn2, t;

    FLINT_ASSERT(d > 0);
    FLINT_ASSERT((fmpz_get_ui(p) & UWORD(3)) == 1);

    fmpz_init(n);
    fmpz_init(halfp);
    fmpz_init(mn2);
    fmpz_init(t);

    fmpz_sub_ui(halfp, p, 1);
    fmpz_fdiv_q_2exp(halfp, halfp, 1);

    fmpz_one(list + 0);
    len = 1;

    for (fmpz_set_ui(n, 2); len < d && fmpz_cmp(n, halfp) <= 0; fmpz_add_ui(n, n, 1))
    {
        int ok = 1;

        fmpz_mod_mul(mn2, n, n, ctx);
        fmpz_mod_neg(mn2, mn2, ctx);

        for (i = 0; ok && i < len; i++)
        {
            fmpz_mod_mul(t, list + i, list + i, ctx);
            ok = !fmpz_equal(t, mn2);
        }

        if (ok)
        {
            fmpz_set(list + len, n);
            len++;
        }
    }
    return len;
}

static int _fill_matrices4(
    fmpz * M,          /* length d by 4d */
    fmpz * Q,          /* length d by 4d+1 */
    slong d,
    const fmpz_mod_ctx_t ctx)
{
    slong i, j;
    fmpz_mod_poly_t g, h;
    fmpz * list;
    fmpz_t g0i, c, t;

    list = _fmpz_vec_init(d);
    if (d != _find_eval_points4(list, d, ctx))
    {
        _fmpz_vec_clear(list, d);
        return 0;
    }

    fmpz_init(g0i);
    fmpz_init(c);
    fmpz_init(t);
    fmpz_mod_poly_init2(g, 4*d + 4, ctx);
    fmpz_mod_poly_init2(h, 4*d + 4, ctx);

    fmpz_mod_poly_one(g, ctx);
    for (i = 0; i < d; i++)
    {
        fmpz_mod_pow_ui(t, list + i, 4, ctx);
        fmpz_mod_neg(t, t, ctx);
        fmpz_mod_poly_shift_left_scalar_addmul_fmpz_mod(g, 4, t, ctx);
    }

    fmpz_mod_inv(g0i, g->coeffs + 0, ctx);
    for (i = 0; i < d; i++)
    {
        FLINT_ASSERT(4*i+4 < g->length);
        fmpz_mod_mul(Q + i*(4*d+1) + 0, g0i, g->coeffs + 4*i+4, ctx);
        fmpz_mod_poly_div_root(c, h, g, list + i, ctx);
        fmpz_mod_poly_evaluate_fmpz(c, h, list + i, ctx);
        fmpz_mod_mul(c, list + i, c, ctx);
        fmpz_mod_inv(c, c, ctx);
        for (j = 0; j < 4*d; j++)
        {
            fmpz_mod_pow_ui(M + i*(4*d) + j, list + i, 1+j, ctx);
            fmpz_mod_mul(Q + (j/4)*(4*d+1) + 4*i + (j%4) + 1, h->coeffs + j, c, ctx);
        }
    }

    fmpz_clear(g0i);
    fmpz_clear(c);
    fmpz_clear(t);
    fmpz_mod_poly_clear(g, ctx);
    fmpz_mod_poly_clear(h, ctx);
    _fmpz_vec_clear(list, d);
    return 1;
}


static void _from_coeffs4(
    fmpz * v,          /* length 4d+1 */
    const fmpz * a,
    slong alen,
    const fmpz * M,    /* length d by 4d */
    slong d,
    const fmpz_t w,
    const fmpz_mod_ctx_t ctx)
{
    slong i, j;
    fmpz_t t1, t2, t3, t4;
    fmpz_t c1, c2, c3, c4;

    FLINT_ASSERT(0 <= alen);
    FLINT_ASSERT(alen <= 1 + 4*d);

    if (alen < 1)
    {
        for (i = 0; i < 4*d+1; i++)
            fmpz_zero(v + i);
        return;
    }
    else if (alen == 1)
    {
        for (i = 0; i < 4*d+1; i++)
            fmpz_set(v + i, a + 0);
        return;
    }

    fmpz_init(t1);
    fmpz_init(t2);
    fmpz_init(t3);
    fmpz_init(t4);
    fmpz_init(c1);
    fmpz_init(c2);
    fmpz_init(c3);
    fmpz_init(c4);

    fmpz_set(v + 0, a + 0);
    for (i = 0; i < d; i++)
    {
        fmpz_zero(c1);
        fmpz_zero(c2);
        fmpz_zero(c3);
        fmpz_zero(c4);

        for (j = 0; j + 4 < alen; j += 4)
        {
            FLINT_ASSERT(j + 3 < 4*d);
            fmpz_addmul(c1, a + j + 1, M + j + 0);
            fmpz_addmul(c2, a + j + 2, M + j + 1);
            fmpz_addmul(c3, a + j + 3, M + j + 2);
            fmpz_addmul(c4, a + j + 4, M + j + 3);
        }

        if (j + 1 < alen)
            fmpz_addmul(c1, a + j + 1, M + j + 0);

        if (j + 2 < alen)
            fmpz_addmul(c1, a + j + 2, M + j + 1);

        if (j + 3 < alen)
            fmpz_addmul(c1, a + j + 3, M + j + 2);

        FLINT_ASSERT(j + 4 >= alen);

        fmpz_mod_set_fmpz(c4, c4, ctx);
        fmpz_mod_set_fmpz(c1, c1, ctx);
        fmpz_mod_set_fmpz(c2, c2, ctx);
        fmpz_mod_set_fmpz(c3, c3, ctx);

        M += 4*d;

        fmpz_mod_add(c4, c4, a + 0, ctx);
        fmpz_mod_add(t1, c4, c2, ctx);
        fmpz_mod_sub(t2, c4, c2, ctx);
        fmpz_mod_add(t3, c1, c3, ctx);
        fmpz_mod_sub(t4, c1, c3, ctx);
        fmpz_mod_mul(t4, t4, w, ctx);

        fmpz_mod_add(v + 4*i + 1, t1, t3, ctx);
        fmpz_mod_add(v + 4*i + 2, t2, t4, ctx);
        fmpz_mod_sub(v + 4*i + 3, t1, t3, ctx);
        fmpz_mod_sub(v + 4*i + 4, t2, t4, ctx);
    }

    fmpz_clear(t1);
    fmpz_clear(t2);
    fmpz_clear(t3);
    fmpz_clear(t4);
    fmpz_clear(c1);
    fmpz_clear(c2);
    fmpz_clear(c3);
    fmpz_clear(c4);
}

static void _to_coeffs4(
    fmpz * a,          /* length 4d+1 */
    const fmpz * v,    /* length 4d+1 */
    fmpz * t,          /* length 4d   */
    const fmpz * Q,    /* length d by 4d+1 */
    slong d,
    const fmpz_t w,
    const fmpz_mod_ctx_t ctx)
{
    slong i, j;
    fmpz_t t1, t2, t3, t4;
    fmpz_t c1, c2, c3, c4;

    fmpz_init(t1);
    fmpz_init(t2);
    fmpz_init(t3);
    fmpz_init(t4);
    fmpz_init(c1);
    fmpz_init(c2);
    fmpz_init(c3);
    fmpz_init(c4);

    fmpz_set(a + 0, v + 0);

    for (i = 0; i < d; i++)
    {
        fmpz_mod_add(t2, v + 1+4*i+0, v + 1+4*i+2, ctx);
        fmpz_mod_sub(t1, v + 1+4*i+0, v + 1+4*i+2, ctx);
        fmpz_mod_add(t3, v + 1+4*i+1, v + 1+4*i+3, ctx);
        fmpz_mod_sub(t4, v + 1+4*i+1, v + 1+4*i+3, ctx);
        fmpz_mod_mul(t4, t4, w, ctx);
        fmpz_mod_sub(t + 4*i+0, t1, t4, ctx);
        fmpz_mod_sub(t + 4*i+1, t2, t3, ctx);
        fmpz_mod_add(t + 4*i+2, t1, t4, ctx);
        fmpz_mod_add(t + 4*i+3, t2, t3, ctx);
    }

    for (i = 0; i < d; i++)
    {
        fmpz_zero(c1);
        fmpz_zero(c2);
        fmpz_zero(c3);
        fmpz_mul(c4, Q + 0, v + 0);

        for (j = 0; j < d; j++)
        {
            fmpz_addmul(c1, t + 4*j + 0, Q + 4*j + 1);
            fmpz_addmul(c2, t + 4*j + 1, Q + 4*j + 2);
            fmpz_addmul(c3, t + 4*j + 2, Q + 4*j + 3);
            fmpz_addmul(c4, t + 4*j + 3, Q + 4*j + 4);
        }

        Q += 4*d + 1;

        fmpz_mod_set_fmpz(a + 4*i + 1, c1, ctx);
        fmpz_mod_set_fmpz(a + 4*i + 2, c2, ctx);
        fmpz_mod_set_fmpz(a + 4*i + 3, c3, ctx);
        fmpz_mod_set_fmpz(a + 4*i + 4, c4, ctx);
    }
}


static int _fill_matrices2(
    fmpz * M,          /* length d by 2d */
    fmpz * Q,          /* length d by 2d+1 */
    slong d,
    const fmpz_mod_ctx_t ctx)
{
    slong i, j;
    fmpz_mod_poly_t g, h;
    fmpz_t g0i, c, t;

    if (fmpz_cmp_ui(fmpz_mod_ctx_modulus(ctx), 2*d) <= 0)
        return 0;

    fmpz_init(g0i);
    fmpz_init(c);
    fmpz_init(t);
    fmpz_mod_poly_init2(g, 2*d + 2, ctx);
    fmpz_mod_poly_init2(h, 2*d + 2, ctx);

    fmpz_mod_poly_one(g, ctx);
    for (i = 0; i < d; i++)
    {
        fmpz_set_ui(t, i + 1);
        fmpz_mod_mul(t, t, t, ctx);
        fmpz_mod_neg(t, t, ctx);
        fmpz_mod_poly_shift_left_scalar_addmul_fmpz_mod(g, 2, t, ctx);
    }

    fmpz_mod_inv(g0i, g->coeffs + 0, ctx);
    for (i = 0; i < d; i++)
    {
        FLINT_ASSERT(2*(i+1) < g->length);
        fmpz_mod_mul(Q + i*(2*d+1) + 0, g0i, g->coeffs + 2*(i+1), ctx);
        fmpz_set_ui(t, i + 1);
        fmpz_mod_poly_div_root(c, h, g, t, ctx);
        FLINT_ASSERT(fmpz_is_zero(c));
        fmpz_mod_poly_evaluate_fmpz(c, h, t, ctx);
        fmpz_mod_mul(c, t, c, ctx);
        fmpz_mod_inv(c, c, ctx);
        for (j = 0; j < 2*d; j++)
        {
            fmpz_mod_pow_ui(M + i*(2*d) + j, t, 1+j, ctx);
            fmpz_mod_mul(Q + (j/2)*(2*d+1) + 2*i + (j%2) + 1, h->coeffs + j, c, ctx);
        }
    }

    fmpz_clear(g0i);
    fmpz_clear(c);
    fmpz_clear(t);
    fmpz_mod_poly_clear(g, ctx);
    fmpz_mod_poly_clear(h, ctx);

    return 1;
}


static void _from_coeffs2(
    fmpz * v,          /* length 2d+1 */
    const fmpz * a,    /* length alen <= 2d+1 */
    slong alen,
    const fmpz * M,    /* length d by 2d */
    slong d,
    const fmpz_mod_ctx_t ctx)
{
    slong i, j;
    fmpz_t c1, c2;

    FLINT_ASSERT(0 <= alen);
    FLINT_ASSERT(alen <= 1 + 2*d);

    if (alen < 1)
    {
        for (i = 0; i < 2*d+1; i++)
            fmpz_zero(v + i);
        return;
    }
    else if (alen == 1)
    {
        for (i = 0; i < 2*d+1; i++)
            fmpz_set(v + i, a + 0);
        return;
    }

    fmpz_init(c1);
    fmpz_init(c2);

    fmpz_set(v + 0, a + 0);
    for (i = 0; i < d; i++)
    {
        fmpz_zero(c1);
        fmpz_zero(c2);

        for (j = 0; j + 2 < alen; j += 2)
        {
            FLINT_ASSERT(j + 1 < 2*d);
            fmpz_addmul(c1, a + j + 1, M + j + 0);
            fmpz_addmul(c2, a + j + 2, M + j + 1);
        }

        if (j + 1 < alen)
            fmpz_addmul(c1, a + j + 1, M + j + 0);

        FLINT_ASSERT(j + 2 >= alen);

        fmpz_mod_set_fmpz(c2, c2, ctx);
        fmpz_mod_set_fmpz(c1, c1, ctx);

        M += 2*d;

        fmpz_mod_add(c2, c2, a + 0, ctx);

        fmpz_mod_add(v + 2*i + 1, c2, c1, ctx);
        fmpz_mod_sub(v + 2*i + 2, c2, c1, ctx);
    }

    fmpz_clear(c1);
    fmpz_clear(c2);
}


static void _to_coeffs2(
    fmpz * a,          /* length 2d+1 */
    const fmpz * v,    /* length 2d+1 */
    fmpz * t,          /* length 2d   */
    const fmpz * Q,    /* length d by 2d+1 */
    slong d,
    const fmpz_mod_ctx_t ctx)
{
    slong i, j;
    fmpz_t c1, c2;

    fmpz_init(c1);
    fmpz_init(c2);

    fmpz_set(a + 0, v + 0);

    for (i = 0; i < d; i++)
    {
        fmpz_mod_sub(t + 2*i+0, v + 1+2*i+0, v + 1+2*i+1, ctx);
        fmpz_mod_add(t + 2*i+1, v + 1+2*i+0, v + 1+2*i+1, ctx);
    }

    for (i = 0; i < d; i++)
    {
        fmpz_zero(c1);
        fmpz_mul(c2, Q + 0, v + 0);
        for (j = 0; j < d; j++)
        {
            fmpz_addmul(c1, t + 2*j + 0, Q + 2*j + 1);
            fmpz_addmul(c2, t + 2*j + 1, Q + 2*j + 2);
        }

        Q += 2*d + 1;

        fmpz_mod_set_fmpz(a + 2*i + 1, c1, ctx);
        fmpz_mod_set_fmpz(a + 2*i + 2, c2, ctx);
    }

    fmpz_clear(c1);
    fmpz_clear(c2);
}


void fmpz_mod_eval_interp_init(fmpz_mod_eval_interp_t E)
{
    E->M = NULL;
    E->T = NULL;
    E->Q = NULL;
    E->array = NULL;
    E->alloc = 0;
    E->d = 0;
    E->radix = 0;
    fmpz_init(E->w);
}

void fmpz_mod_eval_interp_clear(fmpz_mod_eval_interp_t E)
{
    slong i;
    for (i = 0; i < E->alloc; i++)
        fmpz_clear(E->array + i);
    flint_free(E->array);
}

int fmpz_mod_eval_interp_set_degree_modulus(
    fmpz_mod_eval_interp_t E,
    slong deg,
    const fmpz_mod_ctx_t ctx)
{
    slong i, d, new_alloc;
    const fmpz * p = fmpz_mod_ctx_modulus(ctx);

    FLINT_ASSERT(deg >= 0);

    if (fmpz_cmp_ui(p, 3) < 0 || (fmpz_get_ui(p) % 2) == 0 || fmpz_cmp_ui(p, deg) <= 0)
        return 0;

    if ((fmpz_get_ui(p) % 4) == 1)
    {
        d = (deg + 3)/4;
        d = FLINT_MAX(d, 1);

        new_alloc = d*(4*d) + 4*d + d*(4*d + 1);

        for (i = new_alloc; i < E->alloc; i++)
            fmpz_clear(E->array + i);
        E->array = FLINT_ARRAY_REALLOC(E->array, new_alloc, fmpz);
        for (i = E->alloc; i < new_alloc; i++)
            fmpz_init(E->array + i);

        E->radix = 4;
        E->alloc = new_alloc;
        E->d = d;
        E->M = E->array;
        E->T = E->M + d*(4*d);
        E->Q = E->T + 4*d;
        fmpz_mod_set_si(E->w, -1, ctx);
        fmpz_sqrtmod(E->w, E->w, fmpz_mod_ctx_modulus(ctx));

        return _fill_matrices4(E->M, E->Q, d, ctx);
    }
    else
    {
        d = (deg + 1)/2;
        d = FLINT_MAX(d, 1);

        new_alloc = d*(2*d) + 2*d + d*(2*d + 1);

        for (i = new_alloc; i < E->alloc; i++)
            fmpz_clear(E->array + i);
        E->array = FLINT_ARRAY_REALLOC(E->array, new_alloc, fmpz);
        for (i = E->alloc; i < new_alloc; i++)
            fmpz_init(E->array + i);

        E->radix = 2;
        E->alloc = new_alloc;
        E->d = d;
        E->M = E->array;
        E->T = E->M + d*(2*d);
        E->Q = E->T + 2*d;
        fmpz_mod_set_si(E->w, -1, ctx);

        return _fill_matrices2(E->M, E->Q, d, ctx);
    }
}

static void fmpz_mod_eval_interp_to_coeffs(
    fmpz * a,
    const fmpz * v,
    fmpz_mod_eval_interp_t E,
    const fmpz_mod_ctx_t ctx)
{
    if (E->radix == 4)
        _to_coeffs4(a, v, E->T, E->Q, E->d, E->w, ctx);
    else
        _to_coeffs2(a, v, E->T, E->Q, E->d, ctx);
}

static void fmpz_mod_eval_interp_from_coeffs(
    fmpz * v,
    const fmpz * a,
    slong alen,
    fmpz_mod_eval_interp_t E,
    const fmpz_mod_ctx_t ctx)
{
    if (E->radix == 4)
        _from_coeffs4(v, a, alen, E->M, E->d, E->w, ctx);
    else
        _from_coeffs2(v, a, alen, E->M, E->d, ctx);
}


/********** conversion over Fp **********/

void fmpz_mod_eval_interp_to_coeffs_poly(
    fmpz_mod_poly_t a,
    const fmpz_mod_poly_t v,
    fmpz_mod_eval_interp_t E,
    const fmpz_mod_ctx_t ctx)
{
    slong l = fmpz_mod_eval_interp_eval_length(E);
    if (v->length == 0)
    {
        a->length = 0;
        return;
    }
    FLINT_ASSERT(v->length == l);
    fmpz_mod_poly_fit_length(a, l, ctx);
    fmpz_mod_eval_interp_to_coeffs(a->coeffs, v->coeffs, E, ctx);
    a->length = l;
    _fmpz_mod_poly_normalise(a);
}

void fmpz_mod_eval_interp_from_coeffs_poly(
    fmpz_mod_poly_t v,
    const fmpz_mod_poly_t a,
    fmpz_mod_eval_interp_t E,
    const fmpz_mod_ctx_t ctx)
{
    slong l = fmpz_mod_eval_interp_eval_length(E);
    if (a->length == 0)
    {
        v->length = 0;
        return;
    }
    fmpz_mod_poly_fit_length(v, l, ctx);
    v->length = l;
    fmpz_mod_eval_interp_from_coeffs(v->coeffs, a->coeffs, a->length, E, ctx);
}


/*********** arithmetic *****************************************************/

/* a += b */
void fmpz_mod_evals_add_inplace(
    fmpz_mod_poly_t a,
    fmpz_mod_poly_t b,
    slong len,
    const fmpz_mod_ctx_t ctx)
{
    slong i;

    if (b->length == 0)
        return;

    fmpz_mod_poly_fit_length(a, len, ctx);

    if (a->length == 0)
    {
        _fmpz_vec_set(a->coeffs, b->coeffs, len);
        a->length = len;
        return;
    }

    for (i = 0; i < len; i++)
        fmpz_mod_add(a->coeffs + i, a->coeffs + i, b->coeffs + i, ctx);

    a->length = _fmpz_vec_is_zero(a->coeffs, len) ? 0 : len;
}


/* a = b*c */
void fmpz_mod_evals_mul(
    fmpz_mod_poly_t a,
    fmpz_mod_poly_t b,
    fmpz_mod_poly_t c,
    slong len,
    const fmpz_mod_ctx_t ctx)
{
    slong i;

    if (b->length == 0 || c->length == 0)
    {
        a->length = 0;
        return;
    }

    fmpz_mod_poly_fit_length(a, len, ctx);

    for (i = 0; i < len; i++)
        fmpz_mod_mul(a->coeffs + i, b->coeffs + i, c->coeffs + i, ctx);

    a->length = _fmpz_vec_is_zero(a->coeffs, len) ? 0 : len;
}

/* a += b*c */
void fmpz_mod_evals_addmul(
    fmpz_mod_poly_t a,
    fmpz_mod_poly_t b,
    fmpz_mod_poly_t c,
    slong len,
    const fmpz_mod_ctx_t ctx)
{
    slong i;

    if (b->length == 0 || c->length == 0)
        return;

    if (a->length == 0)
    {
        fmpz_mod_evals_mul(a, b, c, len, ctx);
        return;
    }

    for (i = 0; i < len; i++)
    {
        fmpz_addmul(a->coeffs + i, b->coeffs + i, c->coeffs + i);
        fmpz_mod_set_fmpz(a->coeffs + i, a->coeffs + i, ctx);
    }

    a->length = _fmpz_vec_is_zero(a->coeffs, len) ? 0 : len;
}

void fmpz_mod_evals_fmma(
    fmpz_mod_poly_t a,
    fmpz_mod_poly_t b,
    fmpz_mod_poly_t c,
    fmpz_mod_poly_t d,
    fmpz_mod_poly_t e,
    slong len,
    const fmpz_mod_ctx_t ctx)
{
    fmpz_t t1, t2;
    slong i;

    if (b->length == 0 || c->length == 0)
    {
        fmpz_mod_evals_mul(a, d, e, len, ctx);
        return;
    }

    if (d->length == 0 || e->length == 0)
    {
        fmpz_mod_evals_mul(a, b, c, len, ctx);
        return;
    }

    fmpz_mod_poly_fit_length(a, len, ctx);

    fmpz_init(t1);
    fmpz_init(t2);

    for (i = 0; i < len; i++)
    {
        fmpz_mul(t1, b->coeffs + i, c->coeffs + i);
        fmpz_mul(t2, d->coeffs + i, e->coeffs + i);
        fmpz_add(t1, t1, t2);
        fmpz_mod_set_fmpz(a->coeffs + i, t1, ctx);
    }

    fmpz_clear(t1);
    fmpz_clear(t2);

    a->length = _fmpz_vec_is_zero(a->coeffs, len) ? 0 : len;
}

