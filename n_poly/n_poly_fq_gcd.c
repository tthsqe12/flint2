/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "n_poly.h"
#include "ui_factor.h"

slong _n_poly_fq_gcd_euclidean_inplace_(    
    mp_limb_t * A, slong Alen,
    mp_limb_t * B, slong Blen,
    const fq_nmod_ctx_t ctx,
    mp_limb_t * tmp)
{
    slong d = fq_nmod_ctx_degree(ctx);
    slong i;
    mp_limb_t * u = tmp;
    mp_limb_t * q0 = u + d;
    mp_limb_t * q1 = q0 + d;
    mp_limb_t * t = q1 + d;

again:

    if (Alen < 2 || Blen < 2)
    {
        if (Alen < 1)
        {
            if (Blen < 1)
                return 0;

            _n_fq_inv(u, B + d*(Blen - 1), ctx, t);
            for (i = 0; i + 1 < Blen; i++)
                _n_fq_mul(B + d*i, B + d*i, u, ctx, t);
            _n_fq_one(B + d*(Blen - 1), d);
            return -Blen - 1;
        }

        if (Blen < 1)
        {
            _n_fq_inv(u, A + d*(Alen - 1), ctx, t);
            for (i = 0; i + 1 < Alen; i++)
                _n_fq_mul(A + d*i, A + d*i, u, ctx, t);
            _n_fq_one(A + d*(Alen - 1), d);
            return Alen;
        }

        if (Blen < 2)
        {
            _n_fq_one(B + d*0, d);
            return -1 - 1;        
        }

        FLINT_ASSERT(Alen < 2);

        _n_fq_one(A + d*0, d);
        return 1;        
    }

    if (Alen > Blen)
    {
        /* Q = A[Alen-1]/B[Blen-1] x + */
        _n_fq_inv(u, B + d*(Blen - 1), ctx, t);
        _n_fq_mul(q1, A + d*(Alen - 1), u, ctx, t);
        _n_fq_mul(q0, q1, B + d*(Blen - 2), ctx, t);
        _n_fq_sub(q0, q0, A + d*(Alen - 2), ctx);
        _n_fq_mul(q0, q0, u, ctx, t);

        _nmod_vec_neg(q1, q1, d, ctx->mod);

        _n_fq_mul(u, q0, B + d*0, ctx, t);
        _n_fq_add(A + d*(-1 + Alen - Blen),
                  A + d*(-1 + Alen - Blen), u, ctx);

        for (i = 0; i < Blen - 1; i++)
        {
            _n_fq_mul2(t, q1, B + d*i, ctx);
            _n_fq_madd2(t, q0, B + d*(i + 1), ctx, t + 2*d);
            _n_fq_reduce2(u, t, ctx, t + 2*d);
            _n_fq_add(A + d*(i + Alen - Blen),
                      A + d*(i + Alen - Blen), u, ctx);
        }

        Alen -= 2;
        while (Alen > 0 && _n_fq_is_zero(A + d*(Alen - 1), d))
            Alen--;

        goto again;
    }
    else if (Blen > Alen)
    {
        _n_fq_inv(u, A + d*(Alen - 1), ctx, t);
        _n_fq_mul(q1, B + d*(Blen - 1), u, ctx, t);
        _n_fq_mul(q0, q1, A + d*(Alen - 2), ctx, t);
        _n_fq_sub(q0, q0, B + d*(Blen - 2), ctx);
        _n_fq_mul(q0, q0, u, ctx, t);

        _nmod_vec_neg(q1, q1, d, ctx->mod);

        i = -1;
        _n_fq_mul(u, q0, A + d*(i + 1), ctx, t);
        _n_fq_add(B + d*(i + Blen - Alen), B + d*(i + Blen - Alen), u, ctx);

        for (i = 0; i < Alen - 2; i++)
        {
            _n_fq_mul2(t, q1, A + d*i, ctx);
            _n_fq_madd2(t, q0, A + d*(i + 1), ctx, t + 2*d);
            _n_fq_reduce2(u, t, ctx, t + 2*d);
            _n_fq_add(B + d*(i + Blen - Alen), B + d*(i + Blen - Alen), u, ctx);
        }

        Blen -= 2;
        while (Blen > 0 && _n_fq_is_zero(B + d*(Blen - 1), d))
            Blen--;

        goto again;
    }
    else
    {
        _n_fq_inv(u, B + d*(Blen - 1), ctx, t);
        _n_fq_mul(q0, A + d*(Alen - 1), u, ctx, t);

        for (i = 0; i < Blen - 1; i++)
        {
            _n_fq_mul(u, q0, B + d*i, ctx, t);
            _n_fq_sub(A + d*i, A + d*i, u, ctx);
        }

        Alen -= 1;
        while (Alen > 0 && _n_fq_is_zero(A + d*(Alen - 1), d))
            Alen--;

        goto again;
    }
}


void n_poly_fq_gcd_(
    n_poly_t G,
    const n_poly_t A,
    const n_poly_t B,
    const fq_nmod_ctx_t ctx,
    n_poly_stack_t St)
{
    slong d = fq_nmod_ctx_degree(ctx);
    slong n;
    mp_limb_t * a, * b, * t;

    n_poly_stack_fit_request(St, 3);

    t = n_poly_stack_vec_init(St, 8*d);
    a = n_poly_stack_vec_init(St, d*A->length + 1);
    b = n_poly_stack_vec_init(St, d*B->length + 1);

    _nmod_vec_set(a, A->coeffs, d*A->length);
    _nmod_vec_set(b, B->coeffs, d*B->length);

    n = _n_poly_fq_gcd_euclidean_inplace_(a, A->length, b, B->length, ctx, t);

    if (n < 0)
    {
        n = -n - 1;
        n_poly_fit_length(G, d*n);
        _nmod_vec_set(G->coeffs, b, d*n);
        G->length = n;
    }
    else
    {
        n_poly_fit_length(G, d*n);
        _nmod_vec_set(G->coeffs, a, d*n);
        G->length = n;
    }

    n_poly_stack_vec_clear(St);
    n_poly_stack_vec_clear(St);
    n_poly_stack_vec_clear(St);
}

void n_poly_fq_gcd(
    n_poly_t G,
    const n_poly_t A,
    const n_poly_t B,
    const fq_nmod_ctx_t ctx)
{
#if 0
    n_poly_stack_t St;
    n_poly_stack_init(St);
    n_poly_fq_gcd_(G, A, B, ctx, St);
    n_poly_stack_clear(St);
#else
    fq_nmod_poly_t g, a, b;
    fq_nmod_poly_init(g, ctx);
    fq_nmod_poly_init(a, ctx);
    fq_nmod_poly_init(b, ctx);
    n_poly_fq_get_fq_nmod_poly(g, G, ctx);
    n_poly_fq_get_fq_nmod_poly(a, A, ctx);
    n_poly_fq_get_fq_nmod_poly(b, B, ctx);
    fq_nmod_poly_gcd(g, a, b, ctx);
    n_poly_fq_set_fq_nmod_poly(G, g, ctx);
    fq_nmod_poly_clear(g, ctx);
    fq_nmod_poly_clear(a, ctx);
    fq_nmod_poly_clear(b, ctx);
#endif
}

