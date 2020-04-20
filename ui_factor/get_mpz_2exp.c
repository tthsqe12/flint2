/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "ui_factor.h"


/* x should have at least 3*len limbs allocated */
static mp_size_t ui_product_get_mpn(
    mp_limb_t * x,
    const ulong * v,
    slong len,
    slong stride)
{
    mp_size_t xlen, ylen, zlen;
    mp_limb_t out;

    FLINT_ASSERT(0 < len);

    if (len <= 16)
    {
        slong i;
        xlen = 1;
        x[0] = v[0];
        for (i = 1; i < len; i++)
        {
            out = mpn_mul_1(x, x, xlen, v[stride*i]);
            if (out != 0)
            {
                FLINT_ASSERT(xlen < len);
                x[xlen] = out;
                xlen++;
            }
        }
    }
    else
    {
        ylen = ui_product_get_mpn(x + len, v, len/2 + (len % 2), 2*stride);
        zlen = ui_product_get_mpn(x + len + ylen, v + stride, len/2, 2*stride);
        out = (ylen >= zlen) ? mpn_mul(x, x + len, ylen, x + len + ylen, zlen)
                             : mpn_mul(x, x + len + ylen, zlen, x + len, ylen);
        xlen = ylen + zlen - (out == 0);
    }
    FLINT_ASSERT(xlen <= len);
    FLINT_ASSERT(xlen > 0);
    FLINT_ASSERT(x[xlen - 1] != 0);
    return xlen;
}


static __inline__ slong ui_product_one(ulong * terms)
{
    terms[0] = 1;
    return 0;
}


static __inline__ slong ui_product_push(ulong * terms, slong top, ulong b)
{
    ulong hi, lo;
    umul_ppmm(hi, lo, terms[top], b);
    if (hi == 0)
        terms[top] = lo;
    else
        terms[++top] = b;
    return top;
}


/* x*2^e = f, return is e */
ulong ui_factor_get_mpz_2exp(mpz_t x, const ui_factor_t f, ui_factor_stack_t S)
{
    slong fi, ti, tj;
    slong tn, fn = f->length;
    ui_factor_struct * t;
    ui_factor_entry * td;
    const ui_factor_entry * fd = f->data;
    ulong b, p, e = 0;
    mp_limb_t * xd;
    mp_size_t xn;
    mp_limb_t * zd;
    mp_size_t zn;
    mpz_ptr z;
    mpz_ptr y[FLINT_BITS];
    mp_size_t l;
    slong i;
    mp_limb_t out;
    slong top;
    ulong * terms;

    fi = 0;

    if (fi < fn && fd[fi].base == 2)
    {
        e = fd[fi].pow;
        fi++;
    }

    if (fi >= fn)
    {
        mpz_set_ui(x, 1);
        return e;
    }

    ui_factor_stack_fit_request_mpz(S, FLINT_BITS + 2);
    z = ui_factor_stack_take_top_mpz(S);
    terms = flint_mpz_fit_length(ui_factor_stack_take_top_mpz(S), fn);
    for (i = 0; i < FLINT_BITS; i++)
        y[i] = ui_factor_stack_take_top_mpz(S);

    ui_factor_stack_fit_request_factor(S, 1);
    t = ui_factor_stack_take_top_factor(S);
    ui_factor_fit_length(t, fn);
    td = t->data;

    top = ui_product_one(terms);
    tj = 0;
    for ( ; fi < fn; fi++)
    {
        b = fd[fi].base;
        p = fd[fi].pow;

        if (p % 2)
            top = ui_product_push(terms, top, b);

        p = p/2;
        if (p > 0)
        {
            td[tj].base = b;
            td[tj].pow = p;
            tj++;
        }
    }
    tn = tj;

    i = 0;
    y[i]->_mp_size = ui_product_get_mpn(
                   flint_mpz_fit_length(y[i], 3*(top + 1)), terms, top + 1, 1);
    l = mpz_sizeinbase(y[i], 2);

    while (tn > 0)
    {
        top = ui_product_one(terms);
        tj = 0;
        for (ti = 0; ti < tn; ti++)
        {
            b = td[ti].base;
            p = td[ti].pow;

            if (p % 2)
                top = ui_product_push(terms, top, b);

            p = p/2;
            if (p > 0)
            {
                td[tj].base = b;
                td[tj].pow = p;
                tj++;
            }
        }
        tn = tj;

        i++;
        y[i]->_mp_size = ui_product_get_mpn(
                   flint_mpz_fit_length(y[i], 3*(top + 1)), terms, top + 1, 1);
        l += mpz_sizeinbase(y[i], 2) << i;
    }

    /* +1 for extra bits in top limbs, and +1 for extra mpn carry-out */
    l = 1 + 1 + l/FLINT_BITS;

    if (i > 0)
    {
        xd = flint_mpz_fit_length(x, l);
        zd = flint_mpz_fit_length(z, l);

        zn = 2*y[i]->_mp_size;
        FLINT_ASSERT(zn <= l);
        mpn_sqr(zd, y[i]->_mp_d, y[i]->_mp_size);
        zn -= (zd[zn - 1] == 0);

        while (1)
        {
            i--;
            xn = zn + y[i]->_mp_size;
            FLINT_ASSERT(xn <= l);
            out = (zn >= y[i]->_mp_size) ?
                    mpn_mul(xd, zd, zn, y[i]->_mp_d, y[i]->_mp_size) :
                    mpn_mul(xd, y[i]->_mp_d, y[i]->_mp_size, zd, zn);
            xn -= (out == 0);

            if (i <= 0)
                break;

            FLINT_ASSERT(2*xn <= l);
            mpn_sqr(zd, xd, xn);
            zn = 2*xn - (zd[2*xn - 1] == 0);
        }

        x->_mp_size = xn;
    }
    else
    {
        FLINT_ASSERT(i == 0);
        mpz_set(x, y[i]);
    }

    ui_factor_stack_give_back_mpz(S, FLINT_BITS + 2);
    ui_factor_stack_give_back_factor(S, 1);

    return e;
}

