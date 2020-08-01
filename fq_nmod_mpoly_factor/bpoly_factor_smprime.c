/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "nmod_mpoly_factor.h"
#include "fq_nmod_mpoly_factor.h"
#include "ui_factor.h"


void fq_nmod_bpoly_eval(
    fq_nmod_poly_t E,
    const fq_nmod_bpoly_t A,
    const fq_nmod_t alpha,
    const fq_nmod_ctx_t ctx)
{
    slong i;
    fq_nmod_t t;
    fq_nmod_init(t, ctx);
    fq_nmod_poly_zero(E, ctx);
    for (i = A->length - 1; i >= 0; i--)
    {
        fq_nmod_poly_evaluate_fq_nmod(t, A->coeffs + i, alpha, ctx);
        fq_nmod_poly_set_coeff(E, i, t, ctx);
    }
    fq_nmod_clear(t, ctx);
}

void fq_nmod_bpoly_make_monic_series(
    fq_nmod_bpoly_t A,
    const fq_nmod_bpoly_t B,
    slong order,
    const fq_nmod_ctx_t ctx)
{
    slong i;
    fq_nmod_poly_t lcinv;

    FLINT_ASSERT(A != B);
    FLINT_ASSERT(B->length > 0);
    FLINT_ASSERT(fq_nmod_bpoly_is_canonical(B, ctx));

    fq_nmod_poly_init(lcinv, ctx);
    fq_nmod_poly_inv_series(lcinv, B->coeffs + B->length - 1, order, ctx);

    fq_nmod_bpoly_fit_length(A, B->length, ctx);
    for (i = 0; i < B->length; i++)
        fq_nmod_poly_mullow(A->coeffs + i, B->coeffs + i, lcinv, order, ctx);

    A->length = B->length;
    fq_nmod_bpoly_normalise(A, ctx);

    fq_nmod_poly_clear(lcinv, ctx);
}


static void _lattice(
    nmod_mat_t N,
    fq_nmod_bpoly_struct * const * g,
    slong r,
    slong lift_order,
    slong * starts,
    const fq_nmod_bpoly_t f,
    const fq_nmod_ctx_t ctx)
{
    slong i, j, k, l;
    slong deg = ctx->modulus->length;
    fq_nmod_bpoly_t Q, R, dg;
    fq_nmod_bpoly_struct * ld;
    nmod_mat_t M, T1, T2;
    int nlimbs;
    mp_limb_t * trow;

    nlimbs = _nmod_vec_dot_bound_limbs(r, ctx->modulus->mod);
    trow = (mp_limb_t *) flint_malloc(r*sizeof(mp_limb_t));
    fq_nmod_bpoly_init(Q, ctx);
    fq_nmod_bpoly_init(R, ctx);
    fq_nmod_bpoly_init(dg, ctx);
    ld = (fq_nmod_bpoly_struct *) flint_malloc(r*sizeof(fq_nmod_bpoly_struct));
    for (i = 0; i < r; i++)
        fq_nmod_bpoly_init(ld + i, ctx);

    /* init done */

    for (i = 0; i < r; i++)
    {
        fq_nmod_bpoly_divrem_series(Q, R, f, g[i], lift_order, ctx);
        FLINT_ASSERT(R->length == 0);
        fq_nmod_bpoly_derivative(R, g[i], ctx);
        fq_nmod_bpoly_mul_series(ld + i, Q, R, lift_order, ctx);
    }

    for (k = 0; k + 1 < f->length; k++)
    {
        slong d = nmod_mat_nrows(N);

        if (lift_order <= starts[k])
            continue;

        nmod_mat_init(M, deg*(lift_order - starts[k]), d, ctx->modulus->mod.n);

        for (j = starts[k]; j < lift_order; j++)
        for (l = 0; l < deg; l++)
        {
            for (i = 0; i < r; i++)
            {
                if (k >= ld[i].length ||
                    j >= ld[i].coeffs[k].length ||
                    l >= ld[i].coeffs[k].coeffs[j].length)
                {
                    trow[i] = 0;
                }
                else
                {
                    trow[i] = ld[i].coeffs[k].coeffs[j].coeffs[l];
                }
            }

            for (i = 0; i < d; i++)
                nmod_mat_entry(M, (j - starts[k])*deg + l, i) =
                 _nmod_vec_dot(trow, N->rows[i], r, ctx->modulus->mod, nlimbs);
        }

        nmod_mat_init_nullspace_tr(T1, M);

        nmod_mat_init(T2, nmod_mat_nrows(T1), nmod_mat_ncols(N),
                                                          ctx->modulus->mod.n);
        nmod_mat_mul(T2, T1, N);
        nmod_mat_swap(T2, N);
        nmod_mat_rref(N);

        nmod_mat_clear(M);
        nmod_mat_clear(T1);
        nmod_mat_clear(T2);
    }

    flint_free(trow);
    fq_nmod_bpoly_clear(Q, ctx);
    fq_nmod_bpoly_clear(R, ctx);
    fq_nmod_bpoly_clear(dg, ctx);
    for (i = 0; i < r; i++)
        fq_nmod_bpoly_clear(ld + i, ctx);
    flint_free(ld);
}


static int _zassenhaus(
    slong limit,
    fq_nmod_tpoly_t F,
    const fq_nmod_t malpha,
    const nmod_mat_t N,
    fq_nmod_bpoly_struct * const * loc_fac_org,
    slong r,
    slong order,
    const fq_nmod_bpoly_t B,
    const fq_nmod_ctx_t ctx)
{
    int success;
    fq_nmod_bpoly_t Q, R, t1, t2;
    fq_nmod_poly_t g;
    slong * subset;
    slong i, j, s, len, d = nmod_mat_nrows(N);
    fq_nmod_bpoly_struct * loc_fac;
    fq_nmod_bpoly_struct * f;
    fq_nmod_bpoly_t B_copy;

    FLINT_ASSERT(nmod_mat_ncols(N) == r);

    fq_nmod_poly_init(g, ctx);
    fq_nmod_bpoly_init(Q, ctx);
    fq_nmod_bpoly_init(R, ctx);
    fq_nmod_bpoly_init(t1, ctx);
    fq_nmod_bpoly_init(t2, ctx);
    fq_nmod_bpoly_init(B_copy, ctx);

    subset = FLINT_ARRAY_ALLOC(d, slong);
    loc_fac = FLINT_ARRAY_ALLOC(d, fq_nmod_bpoly_struct);
    for (i = 0; i < d; i++)
    {
        subset[i] = i;
        fq_nmod_bpoly_init(loc_fac + i, ctx);
    }

    for (i = 0; i < d; i++)
    {
        fq_nmod_bpoly_one(loc_fac + i, ctx);
        for (j = 0; j < r; j++)
        {
            if (nmod_mat_entry(N, i, j) == 0)
                continue;
            FLINT_ASSERT(nmod_mat_entry(N, i, j) == 1);
            fq_nmod_bpoly_mul_series(t1, loc_fac + i, loc_fac_org[j], order, ctx);
            fq_nmod_bpoly_swap(t1, loc_fac + i, ctx);
        }
    }

    f = (fq_nmod_bpoly_struct *) B;

    len = d;
    for (s = 1; s <= len/2; s++)
    {
        if (s > limit)
        {
            success = 0;
            goto cleanup;
        }

        zassenhaus_subset_first(subset, len, s);
        while (1)
        {
            FLINT_ASSERT(f->length > 0);
            fq_nmod_bpoly_set_poly_var1(t1, f->coeffs + f->length - 1, ctx);
            for (i = 0; i < len; i++)
            {
                if (subset[i] >= 0)
                {
                    fq_nmod_bpoly_mul_series(t2, t1, loc_fac + subset[i], order, ctx);
                    fq_nmod_bpoly_swap(t1, t2, ctx);
                }
            }

            fq_nmod_bpoly_make_primitive(g, t1, ctx);
            if (fq_nmod_bpoly_divides(Q, f, t1, ctx))
            {
                fq_nmod_bpoly_taylor_shift_var1(t1, t1, malpha, ctx);
                fq_nmod_tpoly_fit_length(F, F->length + 1, ctx);
                fq_nmod_bpoly_swap(F->coeffs + F->length, t1, ctx);
                F->length++;
                f = B_copy;
                fq_nmod_bpoly_swap(f, Q, ctx);
                len -= s;
                if (!zassenhaus_subset_next_disjoint(subset, len + s))
                    break;
            }
            else
            {
                if (!zassenhaus_subset_next(subset, len))
                    break;
            }
        }
    }

    if (f->length > 1)
    {
        fq_nmod_tpoly_fit_length(F, F->length + 1, ctx);
        fq_nmod_bpoly_taylor_shift_var1(F->coeffs + F->length, f, malpha, ctx);
        F->length++;
    }
    else
    {
        FLINT_ASSERT(f->length == 1);
        FLINT_ASSERT(fq_nmod_poly_is_one(f->coeffs + 0, ctx));
    }

    success = 1;

cleanup:

    fq_nmod_poly_clear(g, ctx);
    fq_nmod_bpoly_clear(Q, ctx);
    fq_nmod_bpoly_clear(R, ctx);
    fq_nmod_bpoly_clear(t1, ctx);
    fq_nmod_bpoly_clear(t2, ctx);
    fq_nmod_bpoly_clear(B_copy, ctx);

    for (i = 0; i < d; i++)
        fq_nmod_bpoly_clear(loc_fac + i, ctx);
    flint_free(loc_fac);
    flint_free(subset);

    return success;
}

static void _hensel_build_tree(
    slong * link,
    fq_nmod_bpoly_struct * v,
    fq_nmod_bpoly_struct * w,
    const fq_nmod_poly_struct * local_facs,
    slong r,
    const fq_nmod_ctx_t ctx)
{
    slong i, j;
    fq_nmod_poly_t d;
    fq_nmod_poly_struct * V;
    fq_nmod_poly_struct * W;

    V = (fq_nmod_poly_struct *) flint_malloc((2*r - 2)*sizeof(fq_nmod_poly_struct));
    W = (fq_nmod_poly_struct *) flint_malloc((2*r - 2)*sizeof(fq_nmod_poly_struct));

    fq_nmod_poly_init(d, ctx);
    for (i = 0; i < 2*r - 2; i++)
    {
        fq_nmod_poly_init(V + i, ctx);
        fq_nmod_poly_init(W + i, ctx);
    }

    for (i = 0; i < r; i++)
    {
        fq_nmod_poly_set(V + i, local_facs + i, ctx);
        link[i] = -i - 1;
    }

    for (i = r, j = 0; j < 2*r - 4; i++, j += 2)
    {
        slong s;
        slong minp, mind;
        slong tmp;

        minp = j;
        mind = fq_nmod_poly_degree(V + j, ctx);
        for (s = j + 1; s < i; s++)
        {
            if (fq_nmod_poly_degree(V + s, ctx) < mind)
            {
                minp = s;
                mind = fq_nmod_poly_degree(V + s, ctx);
            }
        }
        fq_nmod_poly_swap(V + j, V + minp, ctx);
        tmp = link[j]; link[j] = link[minp]; link[minp] = tmp;

        minp = j + 1;
        mind = fq_nmod_poly_degree(V + j + 1, ctx);
        for (s = j + 2; s < i; s++)
        {
            if (fq_nmod_poly_degree(V + s, ctx) < mind)
            {
                minp = s;
                mind = fq_nmod_poly_degree(V + s, ctx);
            }
        }
        fq_nmod_poly_swap(V + j + 1, V + minp, ctx);
        tmp = link[j + 1]; link[j + 1] = link[minp]; link[minp] = tmp;

        fq_nmod_poly_mul(V + i, V + j, V + j + 1, ctx);
        link[i] = j;
    }

    for (j = 0; j < 2*r - 2; j += 2)
    {
        fq_nmod_poly_xgcd(d, W + j, W + j + 1, V + j, V + j + 1, ctx);
        FLINT_ASSERT(fq_nmod_poly_is_one(d, ctx));
    }

    for (j = 0; j < 2*r - 2; j++)
    {
        fq_nmod_bpoly_set_poly_var0(v + j, V + j, ctx);
        fq_nmod_bpoly_set_poly_var0(w + j, W + j, ctx);
    }

    fq_nmod_poly_clear(d, ctx);
    for (i = 0; i < 2*r - 2; i++)
    {
        fq_nmod_poly_clear(V + i, ctx);
        fq_nmod_poly_clear(W + i, ctx);
    }
    flint_free(V);
    flint_free(W);
}

static void _hensel_lift_fac(
    fq_nmod_bpoly_t G,
    fq_nmod_bpoly_t H,
    const fq_nmod_bpoly_t f,
    fq_nmod_bpoly_t g,
    fq_nmod_bpoly_t h,
    const fq_nmod_bpoly_t a,
    const fq_nmod_bpoly_t b,
    slong p0,
    slong p1,
    const fq_nmod_ctx_t ctx)
{
    slong i, j;
    fq_nmod_bpoly_t c, t1, t2, q, r;

    fq_nmod_bpoly_init(c, ctx);
    fq_nmod_bpoly_init(t1, ctx);
    fq_nmod_bpoly_init(t2, ctx);
    fq_nmod_bpoly_init(q, ctx);
    fq_nmod_bpoly_init(r, ctx);

    fq_nmod_bpoly_mul(t1, g, h, ctx);
    fq_nmod_bpoly_sub(c, f, t1, ctx);
    for (i = 0; i < c->length; i++)
    {
    #if WANT_ASSERT
        {
            fq_nmod_t cc;
            fq_nmod_init(cc, ctx);
            for (j = 0; j < p0; j++)
            {
                fq_nmod_poly_get_coeff(cc, c->coeffs + i, j, ctx);
                FLINT_ASSERT(fq_nmod_is_zero(cc, ctx));
            }
            fq_nmod_clear(cc, ctx);
        }
    #endif
        fq_nmod_poly_shift_right(c->coeffs + i, c->coeffs + i, p0, ctx);
        fq_nmod_poly_truncate(c->coeffs + i, p1, ctx);
    }

    fq_nmod_bpoly_mul_series(t1, c, b, p1, ctx);
    fq_nmod_bpoly_divrem_series(q, r, t1, g, p1, ctx);
    for (i = 0; i < r->length; i++)
        fq_nmod_poly_shift_left(r->coeffs + i, r->coeffs + i, p0, ctx);
    for (i = 0; i < g->length; i++)
        fq_nmod_poly_truncate(g->coeffs + i, p0, ctx);
    fq_nmod_bpoly_add(t1, r, g, ctx);

    fq_nmod_bpoly_mul_series(t2, c, a, p1, ctx);
    fq_nmod_bpoly_divrem_series(q, r, t2, h, p1, ctx);
    for (i = 0; i < r->length; i++)
        fq_nmod_poly_shift_left(r->coeffs + i, r->coeffs + i, p0, ctx);
    for (i = 0; i < h->length; i++)
        fq_nmod_poly_truncate(h->coeffs + i, p0, ctx);
    fq_nmod_bpoly_add(t2, r, h, ctx);

    fq_nmod_bpoly_swap(G, t1, ctx);
    fq_nmod_bpoly_swap(H, t2, ctx);

#if WANT_ASSERT
    {
        fq_nmod_t cc;
        fq_nmod_init(cc, ctx);
        fq_nmod_bpoly_mul(t1, G, H, ctx);
        fq_nmod_bpoly_sub(c, f, t1, ctx);
        for (i = 0; i < c->length; i++)        
        for (j = 0; j < p0 + p1; j++)
        {
            fq_nmod_poly_get_coeff(cc, c->coeffs + i, j, ctx);
            FLINT_ASSERT(fq_nmod_is_zero(cc, ctx));
        }
        fq_nmod_clear(cc, ctx);
    }
#endif

    fq_nmod_bpoly_clear(c, ctx);
    fq_nmod_bpoly_clear(t1, ctx);
    fq_nmod_bpoly_clear(t2, ctx);
    fq_nmod_bpoly_clear(q, ctx);
    fq_nmod_bpoly_clear(r, ctx);
}

static void _hensel_lift_inv(
    fq_nmod_bpoly_t A,
    fq_nmod_bpoly_t B,
    const fq_nmod_bpoly_t G,
    const fq_nmod_bpoly_t H,
    fq_nmod_bpoly_t a,
    fq_nmod_bpoly_t b,
    slong p0,
    slong p1,
    const fq_nmod_ctx_t ctx)
{
    slong i, j;
    fq_nmod_bpoly_t c, t1, t2, q, r;

/*
flint_printf("lift_inv p0: %wd, p1: %wd\n", p0, p1);
flint_printf("G: "); fq_nmod_bpoly_print_pretty(G, "X", "Y", ctx); flint_printf("\n");
flint_printf("H: "); fq_nmod_bpoly_print_pretty(H, "X", "Y", ctx); flint_printf("\n");

flint_printf("a: "); fq_nmod_bpoly_print_pretty(a, "X", "Y", ctx); flint_printf("\n");
flint_printf("b: "); fq_nmod_bpoly_print_pretty(b, "X", "Y", ctx); flint_printf("\n");
*/


    fq_nmod_bpoly_init(c, ctx);
    fq_nmod_bpoly_init(t1, ctx);
    fq_nmod_bpoly_init(t2, ctx);
    fq_nmod_bpoly_init(q, ctx);
    fq_nmod_bpoly_init(r, ctx);

    for (i = 0; i < b->length; i++)
        fq_nmod_poly_truncate(b->coeffs + i, p0, ctx);
    for (i = 0; i < a->length; i++)
        fq_nmod_poly_truncate(a->coeffs + i, p0, ctx);

    fq_nmod_bpoly_mul(t1, G, a, ctx);
    fq_nmod_bpoly_mul(t2, H, b, ctx);
    fq_nmod_bpoly_add(c, t1, t2, ctx);

    FLINT_ASSERT(c->length > 0);
    for (i = 0; i < c->length; i++)
        fq_nmod_poly_neg(c->coeffs + i, c->coeffs + i, ctx);
    fq_nmod_poly_add_si(c->coeffs + 0, c->coeffs + 0, 1, ctx);
    fq_nmod_bpoly_normalise(c, ctx);

    for (i = 0; i < c->length; i++)
    {
    #if WANT_ASSERT
        {
            fq_nmod_t cc;
            fq_nmod_init(cc, ctx);
            for (j = 0; j < p0; j++)
            {
                fq_nmod_poly_get_coeff(cc, c->coeffs + i, j, ctx);
                FLINT_ASSERT(fq_nmod_is_zero(cc, ctx));
            }
            fq_nmod_clear(cc, ctx);
        }
    #endif
        fq_nmod_poly_shift_right(c->coeffs + i, c->coeffs + i, p0, ctx);
        fq_nmod_poly_truncate(c->coeffs + i, p1, ctx);
    }

    fq_nmod_bpoly_mul_series(t1, c, b, p1, ctx);
    fq_nmod_bpoly_divrem_series(q, r, t1, G, p1, ctx);
/*
flint_printf("before shift r: ");
fq_nmod_bpoly_print_pretty(r, "X", "Y", ctx);
flint_printf("\n");
*/

    for (i = 0; i < r->length; i++)
        fq_nmod_poly_shift_left(r->coeffs + i, r->coeffs + i, p0, ctx);
/*
flint_printf("after shift r: ");
fq_nmod_bpoly_print_pretty(r, "X", "Y", ctx);
flint_printf("\n");
*/

    fq_nmod_bpoly_add(t1, r, b, ctx);

    fq_nmod_bpoly_mul_series(t2, c, a, p1, ctx);
    fq_nmod_bpoly_divrem_series(q, r, t2, H, p1, ctx);
/*
flint_printf("before shift r: ");
fq_nmod_bpoly_print_pretty(r, "X", "Y", ctx);
flint_printf("\n");
*/

    for (i = 0; i < r->length; i++)
        fq_nmod_poly_shift_left(r->coeffs + i, r->coeffs + i, p0, ctx);
/*
flint_printf("after shift r: ");
fq_nmod_bpoly_print_pretty(r, "X", "Y", ctx);
flint_printf("\n");
*/


    fq_nmod_bpoly_add(t2, r, a, ctx);

    fq_nmod_bpoly_swap(t1, B, ctx);
    fq_nmod_bpoly_swap(t2, A, ctx);
/*
flint_printf("A: "); fq_nmod_bpoly_print_pretty(A, "X", "Y", ctx); flint_printf("\n");
flint_printf("B: "); fq_nmod_bpoly_print_pretty(B, "X", "Y", ctx); flint_printf("\n");
*/

#if WANT_ASSERT
    fq_nmod_bpoly_mul(t1, G, A, ctx);
    fq_nmod_bpoly_mul(t2, H, B, ctx);
    fq_nmod_bpoly_add(c, t1, t2, ctx);

    FLINT_ASSERT(c->length > 0);
    for (i = 0; i < c->length; i++)
        fq_nmod_poly_neg(c->coeffs + i, c->coeffs + i, ctx);
    fq_nmod_poly_add_si(c->coeffs + 0, c->coeffs + 0, 1, ctx);
    fq_nmod_bpoly_normalise(c, ctx);

/*
flint_printf("p0: %wd, p1: %wd, c: ", p0, p1);
fq_nmod_bpoly_print_pretty(c, "X", "Y", ctx);
flint_printf("\n");
*/
    {
        fq_nmod_t cc;
        fq_nmod_init(cc, ctx);
        for (i = 0; i < c->length; i++)        
        for (j = 0; j < p0 + p1; j++)
        {
            fq_nmod_poly_get_coeff(cc, c->coeffs + i, j, ctx);
            FLINT_ASSERT(fq_nmod_is_zero(cc, ctx));
        }
        fq_nmod_clear(cc, ctx);
    }
#endif

    fq_nmod_bpoly_clear(c, ctx);
    fq_nmod_bpoly_clear(t1, ctx);
    fq_nmod_bpoly_clear(t2, ctx);
    fq_nmod_bpoly_clear(q, ctx);
    fq_nmod_bpoly_clear(r, ctx); 
}


static void _hensel_lift_tree(
    int opt,
    slong * link,
    fq_nmod_bpoly_struct * v,
    fq_nmod_bpoly_struct * w,
    const fq_nmod_bpoly_t f,
    slong j,
    slong p0,
    slong p1,
    const fq_nmod_ctx_t ctx)
{
    FLINT_ASSERT(p1 <= p0);

    if (j < 0)
        return;

    if (opt >= 0)
        _hensel_lift_fac(v + j, v + j + 1,
                           f, v + j, v + j + 1, w + j, w + j + 1, p0, p1, ctx);

    if (opt <= 0)
        _hensel_lift_inv(w + j, w + j + 1,
                              v + j, v + j + 1, w + j, w + j + 1, p0, p1, ctx);

    _hensel_lift_tree(opt, link, v, w, v + j, link[j], p0, p1, ctx);
    _hensel_lift_tree(opt, link, v, w, v + j + 1, link[j + 1], p0, p1, ctx);
}


int fq_nmod_bpoly_factor_smprime(
    fq_nmod_poly_t c,
    fq_nmod_tpoly_t F,
    fq_nmod_bpoly_t B,
    int allow_shift,
    const fq_nmod_ctx_t ctx)
{
    int success;
    slong i, r;
    slong Blenx = B->length;
    slong Bleny;
    slong final_order, curr_lift_order, prev_lift_order, next_lift_order;
    slong * starts;
    fq_nmod_t alpha, Blc;
    fq_nmod_poly_t Beval;
    fq_nmod_poly_factor_t local_fac;
    fq_nmod_bpoly_t monicB;
    nmod_mat_t N;
    slong * link;
    fq_nmod_bpoly_struct * v, * w, ** lift_fac;
    fq_nmod_tpoly_t tmp;
    slong e[FLINT_BITS];
    slong old_nrows;
    slong zas_limit;
/*
flint_printf("fq_nmod_bpoly_factor_smprime called allow_shift = %d\n", allow_shift);
*/
    FLINT_ASSERT(Blenx > 1);

    fq_nmod_init(alpha, ctx);
    fq_nmod_init(Blc, ctx);
    fq_nmod_poly_init(Beval, ctx);
    fq_nmod_poly_factor_init(local_fac, ctx);
    fq_nmod_bpoly_init(monicB, ctx);
    fq_nmod_tpoly_init(tmp, ctx);
    nmod_mat_init(N, 0, 0, ctx->modulus->mod.n);
    starts = (slong *) flint_malloc(Blenx*sizeof(slong));
    link = (slong *) flint_malloc(sizeof(slong));
    lift_fac = (fq_nmod_bpoly_struct **) flint_malloc(
                                               sizeof(fq_nmod_bpoly_struct *));

    /* init done */

    fq_nmod_bpoly_make_primitive(c, B, ctx);

    /* deg_y(B) + 1 */
    Bleny = 0;
    for (i = 0; i < B->length; i++)
        Bleny = FLINT_MAX(Bleny, (B->coeffs + i)->length);

    /* CLD bounds */
    for (i = 0; i < Blenx; i++)
        starts[i] = Bleny;

    fq_nmod_zero(alpha, ctx);
    goto got_alpha;

next_alpha:

    if (!allow_shift || !fq_nmod_next(alpha, ctx))
    {
        success = 0;
        goto cleanup;
    }

got_alpha:

    fq_nmod_bpoly_eval(Beval, B, alpha, ctx);

    /* if killed leading coeff, get new alpha */
    if (Beval->length != Blenx)
        goto next_alpha;

    local_fac->num = 0;
    fq_nmod_poly_factor(local_fac, Blc, Beval, ctx);

    r = local_fac->num;

    /* if multiple factors, get new alpha */
    for (i = 0; i < r; i++)
    {
        if (local_fac->exp[i] != 1)
            goto next_alpha;
    }

	/* if one factor, A is irreducible */
	if (r < 2)
	{
        fq_nmod_tpoly_fit_length(F, 1, ctx);
        fq_nmod_bpoly_swap(F->coeffs + 0, B, ctx);
        F->length = 1;
        success = 1;
        goto cleanup;
	}

    /* done if B is constant in y */
    if (Bleny < 2)
    {
        fq_nmod_tpoly_fit_length(F, r, ctx);
        F->length = r;
        for (i = 0; i < r; i++)
            fq_nmod_bpoly_set_poly_var0(F->coeffs + i, local_fac->poly + i, ctx);
        success = 1;
        goto cleanup;
    }

    /* precision for constructing true factors */
    final_order = Bleny;

    fq_nmod_bpoly_taylor_shift_var1(B, B, alpha, ctx);

    nmod_mat_clear(N);
    nmod_mat_init(N, r, r, ctx->modulus->mod.n);
    for (i = 0; i < r; i++)
        nmod_mat_entry(N, i, i) = 1;

    link = (slong *) flint_realloc(link, (2*r - 2)*sizeof(slong));
    lift_fac = (fq_nmod_bpoly_struct **) flint_realloc(lift_fac,
                                               r*sizeof(fq_nmod_bpoly_struct));

    fq_nmod_tpoly_fit_length(tmp, 2*(2*r - 2), ctx);
    v = tmp->coeffs + 0;
    w = tmp->coeffs + (2*r - 2);

    curr_lift_order = final_order + r;

    fq_nmod_bpoly_make_monic_series(monicB, B, curr_lift_order, ctx);

    _hensel_build_tree(link, v, w, local_fac->poly, r, ctx);
    for (i = 0; i < 2*r - 2; i++)
        if (-link[i] - 1 >= 0)
            lift_fac[-link[i] - 1] = v + i;

    FLINT_ASSERT(curr_lift_order > 1);
    for (i = 0, e[i] = curr_lift_order; e[i] > 1; i++)
        e[i+1] = (e[i] + 1) / 2;

    for (i--; i > 0; i--)
        _hensel_lift_tree(0, link, v, w, monicB, 2*r-4, e[i+1], e[i]-e[i+1], ctx);

    prev_lift_order = e[1];
    _hensel_lift_tree(1, link, v, w, monicB, 2*r-4, e[1], e[0]-e[1], ctx);

    zas_limit = 1;

try_zas:

    F->length = 0;
    fq_nmod_neg(alpha, alpha, ctx);
    success = _zassenhaus(zas_limit, F, alpha, N, lift_fac, r, final_order, B, ctx);
    fq_nmod_neg(alpha, alpha, ctx);
    if (success)
        goto cleanup;

    zas_limit = 2;

more:

    old_nrows = nmod_mat_nrows(N);
    _lattice(N, lift_fac, r, curr_lift_order, starts, B, ctx);
    if (nmod_mat_nrows(N) < old_nrows && nmod_mat_is_reduced(N))
        goto try_zas;

    next_lift_order = curr_lift_order + r;
    next_lift_order = FLINT_MIN(next_lift_order, 2*curr_lift_order);

    _hensel_lift_tree(-1, link, v, w, monicB, 2*r-4, prev_lift_order,
                                       curr_lift_order - prev_lift_order, ctx);

    fq_nmod_bpoly_make_monic_series(monicB, B, next_lift_order, ctx);

    _hensel_lift_tree(0, link, v, w, monicB, 2*r-4, curr_lift_order,
                                       next_lift_order - curr_lift_order, ctx);

    prev_lift_order = curr_lift_order;
    curr_lift_order = next_lift_order;

    goto more;

cleanup:

    flint_free(starts);
    flint_free(link);
    flint_free(lift_fac);

    nmod_mat_clear(N);
    fq_nmod_poly_clear(Beval, ctx);
    fq_nmod_poly_factor_clear(local_fac, ctx);
    fq_nmod_bpoly_clear(monicB, ctx);
    fq_nmod_tpoly_clear(tmp, ctx);

    fq_nmod_clear(alpha, ctx);
    fq_nmod_clear(Blc, ctx);
/*
flint_printf("fq_nmod_bpoly_factor_smprime returning %d\n", success);
flint_printf("F->length: %wd\n", F->length);
*/
    return success;
}

