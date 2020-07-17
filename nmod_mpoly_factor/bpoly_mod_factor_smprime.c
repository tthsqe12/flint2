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
#include "profiler.h"

int nmod_mat_is_reduced(const nmod_mat_t N)
{
    slong i, j, k = 0;
    slong r = nmod_mat_ncols(N);
    slong d = nmod_mat_nrows(N);
    
    for (i = 0; i < d; i++)
    for (j = 0; j < r; j++)
    {
        if (nmod_mat_entry(N, i, j) != 0)
        {
            if (nmod_mat_entry(N, i, j) == 1)
                k++;
            else
                return 0;
        }   
    }
    return k == r;
}

void nmod_mat_init_nullspace_tr(nmod_mat_t X, nmod_mat_t tmp)
{
    slong i, j, k, m, n, rank, nullity;
    slong * p;
    slong * pivots;
    slong * nonpivots;

    m = tmp->r;
    n = tmp->c;

    p = flint_malloc(sizeof(slong) * FLINT_MAX(m, n));

    rank = nmod_mat_rref(tmp);

    nullity = n - rank;

    nmod_mat_init(X, nullity, n, tmp->mod.n);

    if (rank == 0)
    {
        for (i = 0; i < nullity; i++)
            nmod_mat_entry(X, i, i) = UWORD(1);
    }
    else if (nullity)
    {
        pivots = p;            /* length = rank */
        nonpivots = p + rank;  /* length = nullity */

        for (i = j = k = 0; i < rank; i++)
        {
            while (nmod_mat_entry(tmp, i, j) == UWORD(0))
            {
                nonpivots[k] = j;
                k++;
                j++;
            }
            pivots[i] = j;
            j++;
        }
        while (k < nullity)
        {
            nonpivots[k] = j;
            k++;
            j++;
        }

        for (i = 0; i < nullity; i++)
        {
            for (j = 0; j < rank; j++)
            {
                mp_limb_t c = nmod_mat_entry(tmp, j, nonpivots[i]);
                nmod_mat_entry(X, i, pivots[j]) = nmod_neg(c, tmp->mod);
            }

            nmod_mat_entry(X, i, nonpivots[i]) = UWORD(1);
        }
    }

    flint_free(p);
}



static void n_bpoly_mod_add_poly_shift(
    n_bpoly_t A,
    const n_poly_t B,
    slong yshift,
    const nmod_t ctx)
{
    slong i;

    FLINT_ASSERT(A->length > B->length);

    for (i = 0; i < B->length; i++)
    {
        FLINT_ASSERT(0 == n_poly_get_coeff(A->coeffs + i, yshift));
        n_poly_set_coeff(A->coeffs + i, yshift, B->coeffs[i]);
    }
}

void n_bpoly_mod_make_monic_series(
    n_bpoly_t A,
    const n_bpoly_t B,
    slong order,
    nmod_t ctx)
{
    slong i;
    n_poly_t lcinv;

    FLINT_ASSERT(A != B);
    FLINT_ASSERT(B->length > 0);
    FLINT_ASSERT(n_bpoly_mod_is_canonical(B, ctx));

    n_poly_init(lcinv);
    n_poly_mod_inv_series(lcinv, B->coeffs + B->length - 1, order, ctx);

    n_bpoly_fit_length(A, B->length);
    A->length = B->length;
    for (i = 0; i < B->length; i++)
        n_poly_mod_mullow(A->coeffs + i, B->coeffs + i, lcinv, order, ctx);

    n_poly_clear(lcinv);
}

int nmod_partial_fraction_coeffs(
    slong r,
    n_poly_struct * out,
    const n_poly_struct * f,
    nmod_t ctx)
{
    slong i;
    n_poly_t num, den, a, b, g, t;

    FLINT_ASSERT(r >= 2);

    n_poly_init(num);
    n_poly_init(den);
    n_poly_init(a);
    n_poly_init(b);
    n_poly_init(g);
    n_poly_init(t);

    n_poly_one(num);
    _n_poly_mod_mul(den, f + 0, f + 1, ctx);
    for (i = 2; i < r; i++)
    {
        _n_poly_mod_mul(g, den, f + i, ctx);
        n_poly_swap(g, den);
    }

    for (i = 0; i < r; i++)
    {
        n_poly_mod_divrem(den, t, den, f + i, ctx);
        FLINT_ASSERT(n_poly_is_zero(t));
        n_poly_mod_xgcd(g, a, b, f + i, den, ctx);
        if (n_poly_degree(g) != 0)
            return 0; /* TODO leak */
        FLINT_ASSERT(n_poly_is_one(g));
        n_poly_mod_mul(t, b, num, ctx);
        n_poly_mod_rem(out + i, t, f + i, ctx);
        n_poly_mod_mul(t, a, num, ctx);
        n_poly_mod_rem(num, t, den, ctx);
    }

    n_poly_clear(num);
    n_poly_clear(den);
    n_poly_clear(a);
    n_poly_clear(b);
    n_poly_clear(g);
    n_poly_clear(t);
    return 1;
}

void n_bpoly_mod_eval(
    n_poly_t E,
    const n_bpoly_t A,
    mp_limb_t alpha,
    nmod_t ctx)
{
    slong i;
    n_poly_zero(E);
    for (i = A->length - 1; i >= 0; i--)
    {
        mp_limb_t c = n_poly_mod_evaluate_nmod(A->coeffs + i, alpha, ctx);
        n_poly_set_coeff(E, i, c);
    }
}


static void _recombine_lattice(
    nmod_mat_t N,
    const n_bpoly_struct * g,
    slong r,
    slong lift_order,
    slong * starts,
    const n_bpoly_t f,
    nmod_t ctx)
{
    slong i, j, k;
    n_bpoly_t Q, R, dg;
    n_bpoly_struct * ld;
    nmod_mat_t M, T1, T2;
    int nlimbs;
    mp_limb_t * trow;

    nlimbs = _nmod_vec_dot_bound_limbs(r, ctx);
    trow = (mp_limb_t *) flint_malloc(r*sizeof(mp_limb_t));
    n_bpoly_init(Q);
    n_bpoly_init(R);
    n_bpoly_init(dg);
    ld = (n_bpoly_struct *) flint_malloc(r*sizeof(n_bpoly_struct));
    for (i = 0; i < r; i++)
        n_bpoly_init(ld + i);

    /* init done */

    for (i = 0; i < r; i++)
    {
        n_bpoly_mod_divrem_series(Q, R, f, g + i, lift_order, ctx);
        FLINT_ASSERT(R->length == 0);
        n_bpoly_mod_derivative(R, g + i, ctx);
        n_bpoly_mod_mul_series(ld + i, Q, R, lift_order, ctx);
    }

    for (k = 0; k + 1 < f->length; k++)
    {
        slong d = nmod_mat_nrows(N);

        if (lift_order <= starts[k])
            continue;

        nmod_mat_init(M, lift_order - starts[k], d, ctx.n);

        for (j = starts[k]; j < lift_order; j++)
        {
            for (i = 0; i < r; i++)
                trow[i] = n_bpoly_get_coeff(ld + i, k, j);

            for (i = 0; i < d; i++)
                nmod_mat_entry(M, j - starts[k], i) =
                             _nmod_vec_dot(trow, N->rows[i], r, ctx, nlimbs);
        }

        nmod_mat_init_nullspace_tr(T1, M);

        nmod_mat_init(T2, nmod_mat_nrows(T1), nmod_mat_ncols(N), ctx.n);
        nmod_mat_mul(T2, T1, N);
        nmod_mat_swap(T2, N);
        nmod_mat_rref(N);

        nmod_mat_clear(M);
        nmod_mat_clear(T1);
        nmod_mat_clear(T2);
    }

    flint_free(trow);
    n_bpoly_clear(Q);
    n_bpoly_clear(R);
    n_bpoly_clear(dg);
    for (i = 0; i < r; i++)
        n_bpoly_clear(ld + i);
    flint_free(ld);
}


static void _recombine_zassenhaus(
    n_tpoly_t F,
    mp_limb_t malpha,
    const nmod_mat_t N,
    const n_bpoly_struct * loc_fac_org,
    slong r,
    slong order,
    n_bpoly_t f,
    nmod_t ctx)
{
    n_bpoly_t Q, R, t1, t2;
    n_poly_t leadf, g;
    slong * idx;
    slong i, j, s, len, d = nmod_mat_nrows(N);
    fmpz_t subset;
    n_bpoly_struct * loc_fac;

    FLINT_ASSERT(nmod_mat_ncols(N) == r);

    loc_fac = (n_bpoly_struct *) flint_malloc(d*sizeof(n_bpoly_struct));
    for (i = 0; i < d; i++)
        n_bpoly_init(loc_fac + i);

    idx = (slong *) flint_malloc(r * sizeof(slong));
    for (i = 0; i < r; i++)
        idx[i] = i;

    fmpz_init(subset);

    n_poly_init(g);
    n_bpoly_init(Q);
    n_bpoly_init(R);
    n_bpoly_init(t1);
    n_bpoly_init(t2);
    n_poly_init(leadf);

    for (i = 0; i < d; i++)
    {
        n_bpoly_one(loc_fac + i);
        for (j = 0; j < r; j++)
        {
            if (nmod_mat_entry(N, i, j) == 0)
                continue;
            n_bpoly_mod_mul_series(t1, loc_fac + i, loc_fac_org + j, order, ctx);
            n_bpoly_swap(t1, loc_fac + i);
        }
    }

    FLINT_ASSERT(f->length > 0);
    n_poly_set(leadf, f->coeffs + f->length - 1);

    len = d;
    for (s = 1; s <= len/2; s++)
    {
        subset_first(subset, len, s);
        do {
try_subset:
            n_bpoly_set_poly_var1(t1, leadf);
            for (i = 0; i < len; i++)
            {
                if (fmpz_tstbit(subset, i))
                {
                    n_bpoly_mod_mul_series(t2, t1, loc_fac + idx[i], order, ctx);
                    n_bpoly_swap(t1, t2);
                }
            }

            n_bpoly_mod_make_primitive(g, t1, ctx);
            if (n_bpoly_mod_divides(Q, f, t1, ctx))
            {
                n_bpoly_mod_taylor_shift_var1(t1, malpha, ctx);
                n_tpoly_fit_length(F, F->length + 1);
                n_bpoly_swap(F->coeffs + F->length, t1);
                F->length++;
                n_bpoly_swap(f, Q);
                FLINT_ASSERT(f->length > 0);
                n_poly_set(leadf, f->coeffs + f->length - 1);

                for (j = 0, i = 0; i < len; i++)
                    if (!fmpz_tstbit(subset, i))
                        idx[j++] = idx[i];
                len -= s;

                if (!subset_fix(subset, len + s))
                    goto sloop_continue;
                goto try_subset;
            }
        }
        while (subset_next(subset, subset, len));
sloop_continue:
        (void)(NULL);
    }

    n_bpoly_mod_taylor_shift_var1(f, malpha, ctx);
    if (f->length > 1)
    {
        n_tpoly_fit_length(F, F->length + 1);
        n_bpoly_swap(F->coeffs + F->length, f);
        F->length++;
    }
    else
    {
        FLINT_ASSERT(f->length == 1);
        FLINT_ASSERT(n_poly_is_one(f->coeffs + 0));
    }

    fmpz_clear(subset);

    n_poly_clear(g);
    n_bpoly_clear(Q);
    n_bpoly_clear(R);
    n_bpoly_clear(t1);
    n_bpoly_clear(t2);
    n_poly_clear(leadf);

    for (i = 0; i < d; i++)
        n_bpoly_clear(loc_fac + i);
    flint_free(loc_fac);

    flint_free(idx);
}


static void _lift_quintic(
    n_bpoly_struct * lift_facs,
    n_poly_struct * local_facs,
    slong r,
    const n_bpoly_t IBtilde,
    slong lift_order,
    nmod_t ctx)
{
    slong i, j, k;
    n_bpoly_t Id;
    n_bpoly_t tp, tp1, error;
    n_poly_t ss, tt;

    n_bpoly_init(Id);
    n_poly_init(ss);
    n_poly_init(tt);
    n_bpoly_init(tp);
    n_bpoly_init(tp1);
    n_bpoly_init(error);

    n_bpoly_fit_length(Id, r);
    Id->length = r;

    for (i = 0; i < r; i++)
    {
        n_bpoly_set_poly_var0(lift_facs + i, local_facs + i);
    }

    nmod_partial_fraction_coeffs(r, Id->coeffs, local_facs, ctx);

    n_bpoly_mod_mul_series(tp, lift_facs + 0, lift_facs + 1, lift_order, ctx);
    for (i = 2; i < r; i++)
    {
        n_bpoly_mod_mul_series(tp1, tp, lift_facs + i, lift_order, ctx);
        n_bpoly_swap(tp1, tp);
    }

    n_bpoly_mod_sub(error, IBtilde, tp, ctx);

    for (j = 1; j < lift_order; j++)
    {
        n_poly_zero(ss);
        for (i = error->length - 1; i >= 0; i--)
        {
            n_poly_set_coeff(ss, i, n_bpoly_get_coeff(error, i, j));
            for (k = 0; k < j; k++)
            {
                FLINT_ASSERT(0 == n_bpoly_get_coeff(error, i, k));
            }
        }

        for (i = 0; i < r; i++)
        {
            n_poly_mod_mul(tt, ss, Id->coeffs + i, ctx);
            n_poly_mod_rem(tt, tt, local_facs + i, ctx);
            n_bpoly_mod_add_poly_shift(lift_facs + i, tt, j, ctx);
        }

        n_bpoly_mod_mul_series(tp, lift_facs + 0, lift_facs + 1, lift_order, ctx);
        for (i = 2; i < r; i++)
        {
            n_bpoly_mod_mul_series(tp1, tp, lift_facs + i, lift_order, ctx);
            n_bpoly_swap(tp1, tp);
        }
        n_bpoly_mod_sub(error, IBtilde, tp, ctx);
    }

    n_bpoly_clear(Id);
    n_poly_clear(ss);
    n_poly_clear(tt);
    n_bpoly_clear(tp);
    n_bpoly_clear(tp1);
    n_bpoly_clear(error);
}


static void _hensel_build_tree(
    slong * link,
    n_bpoly_struct * v,
    n_bpoly_struct * w,
    const n_poly_struct * local_facs,
    slong r,
    nmod_t ctx)
{
    slong i, j;
    n_poly_t d;
    n_poly_struct * V;
    n_poly_struct * W;

    V = flint_malloc((2*r - 2)*sizeof(nmod_poly_t));
    W = flint_malloc((2*r - 2)*sizeof(nmod_poly_t));

    n_poly_init(d);
    for (i = 0; i < 2*r - 2; i++)
    {
        n_poly_init(V + i);
        n_poly_init(W + i);
    }

    for (i = 0; i < r; i++)
    {
        n_poly_set(V + i, local_facs + i);
        link[i] = -i - 1;
    }

    for (i = r, j = 0; j < 2*r - 4; i++, j += 2)
    {
        slong s;
        slong minp, mind;
        slong tmp;

        minp = j;
        mind = n_poly_degree(V + j);
        for (s = j + 1; s < i; s++)
        {
            if (n_poly_degree(V + s) < mind)
            {
                minp = s;
                mind = n_poly_degree(V + s);
            }
        }

        n_poly_swap(V + j, V + minp);
        tmp = link[j]; link[j] = link[minp]; link[minp] = tmp;

        minp = j + 1;
        mind = n_poly_degree(V + j + 1);

        for (s = j + 2; s < i; s++)
        {
            if (n_poly_degree(V + s) < mind)
            {
                minp = s;
                mind = n_poly_degree(V + s);
            }
        }

        n_poly_swap(V + j + 1, V + minp);
        tmp = link[j + 1]; link[j + 1] = link[minp]; link[minp] = tmp;

        n_poly_mod_mul(V + i, V + j, V + j + 1, ctx);
        link[i] = j;
    }

    for (j = 0; j < 2*r - 2; j += 2)
    {
        n_poly_mod_xgcd(d, W + j, W + j + 1, V + j, V + j + 1, ctx);
        FLINT_ASSERT(n_poly_is_one(d));
    }

    for (j = 0; j < 2*r - 2; j++)
    {
        n_bpoly_set_poly_var0(v + j, V + j);
        n_bpoly_set_poly_var0(w + j, W + j);
    }

    n_poly_clear(d);
    for (i = 0; i < 2*r - 2; i++)
    {
        n_poly_clear(V + i);
        n_poly_clear(W + i);
    }
    flint_free(V);
    flint_free(W);
}

static void _hensel_lift_fac(
    n_bpoly_t G,
    n_bpoly_t H,
    const n_bpoly_t f,
    n_bpoly_t g,
    n_bpoly_t h,
    const n_bpoly_t a,
    const n_bpoly_t b,
    slong p0,
    slong p1,
    nmod_t ctx)
{
    slong i, j;
    n_bpoly_t c, t1, t2, q, r;

    n_bpoly_init(c);
    n_bpoly_init(t1);
    n_bpoly_init(t2);
    n_bpoly_init(q);
    n_bpoly_init(r);

    n_bpoly_mod_mul(t1, g, h, ctx);
    n_bpoly_mod_sub(c, f, t1, ctx);
    for (i = 0; i < c->length; i++)
    {
        for (j = 0; j < p0; j++)
            FLINT_ASSERT(n_poly_get_coeff(c->coeffs + i, j) == 0);
        n_poly_shift_right(c->coeffs + i, c->coeffs + i, p0);
        n_poly_truncate(c->coeffs + i, p1);
    }

    n_bpoly_mod_mul_series(t1, c, b, p1, ctx);
    n_bpoly_mod_divrem_series(q, r, t1, g, p1, ctx);
    for (i = 0; i < r->length; i++)
        n_poly_shift_left(r->coeffs + i, r->coeffs + i, p0);
    for (i = 0; i < g->length; i++)
        n_poly_truncate(g->coeffs + i, p0);
    n_bpoly_mod_add(t1, r, g, ctx);

    n_bpoly_mod_mul_series(t2, c, a, p1, ctx);
    n_bpoly_mod_divrem_series(q, r, t2, h, p1, ctx);
    for (i = 0; i < r->length; i++)
        n_poly_shift_left(r->coeffs + i, r->coeffs + i, p0);
    for (i = 0; i < h->length; i++)
        n_poly_truncate(h->coeffs + i, p0);
    n_bpoly_mod_add(t2, r, h, ctx);

    n_bpoly_swap(G, t1);
    n_bpoly_swap(H, t2);

    n_bpoly_clear(c);
    n_bpoly_clear(t1);
    n_bpoly_clear(t2);
    n_bpoly_clear(q);
    n_bpoly_clear(r);
}

static void _hensel_lift_inv(
    n_bpoly_t A,
    n_bpoly_t B,
    const n_bpoly_t G,
    const n_bpoly_t H,
    n_bpoly_t a,
    n_bpoly_t b,
    slong p0,
    slong p1,
    nmod_t ctx)
{
    slong i, j;
    n_bpoly_t c, t1, t2, q, r;

    n_bpoly_init(c);
    n_bpoly_init(t1);
    n_bpoly_init(t2);
    n_bpoly_init(q);
    n_bpoly_init(r);

    n_bpoly_mod_mul(t1, G, a, ctx);
    n_bpoly_mod_mul(t2, H, b, ctx);
    n_bpoly_mod_add(c, t1, t2, ctx);

    FLINT_ASSERT(c->length > 0);
    for (i = 0; i < c->length; i++)
        n_poly_mod_neg(c->coeffs + i, c->coeffs + i, ctx);
    n_poly_mod_add_ui(c->coeffs + 0, c->coeffs + 0, 1, ctx);
    n_bpoly_normalise(c);

    for (i = 0; i < c->length; i++)
    {
        for (j = 0; j < p0; j++)
            FLINT_ASSERT(n_poly_get_coeff(c->coeffs + i, j) == 0);
        n_poly_shift_right(c->coeffs + i, c->coeffs + i, p0);
        n_poly_truncate(c->coeffs + i, p1);
    }

    n_bpoly_mod_mul_series(t1, c, b, p1, ctx);
    n_bpoly_mod_divrem_series(q, r, t1, G, p1, ctx);
    for (i = 0; i < r->length; i++)
        n_poly_shift_left(r->coeffs + i, r->coeffs + i, p0);
    for (i = 0; i < b->length; i++)
        n_poly_truncate(b->coeffs + i, p0);
    n_bpoly_mod_add(t1, r, b, ctx);

    n_bpoly_mod_mul_series(t2, c, a, p1, ctx);
    n_bpoly_mod_divrem_series(q, r, t2, H, p1, ctx);
    for (i = 0; i < r->length; i++)
        n_poly_shift_left(r->coeffs + i, r->coeffs + i, p0);
    for (i = 0; i < a->length; i++)
        n_poly_truncate(a->coeffs + i, p0);
    n_bpoly_mod_add(t2, r, a, ctx);

    n_bpoly_swap(t1, B);
    n_bpoly_swap(t2, A);

    n_bpoly_clear(c);
    n_bpoly_clear(t1);
    n_bpoly_clear(t2);
    n_bpoly_clear(q);
    n_bpoly_clear(r); 
}

static void _hensel_lift_tree(
    int opt,
    slong * link,
    n_bpoly_struct * v,
    n_bpoly_struct * w,
    n_bpoly_t f,
    slong j,
    slong p0,
    slong p1,
    nmod_t ctx)
{
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

int n_bpoly_mod_factor_smprime(
    n_poly_t c,
    n_tpoly_t F,
    n_bpoly_t B,
    int allow_shift,
    nmod_t ctx)
{
    int success;
    slong i, r;
    slong Blenx = B->length;
    slong Bleny;
    slong final_order, lift_order;
    slong * starts;
    mp_limb_t alpha;
    n_poly_t Beval;
    n_poly_factor_t local_fac;
    n_bpoly_t monicB;
    n_tpoly_t lift_fac;
    nmod_mat_t N;
timeit_t timer;

    FLINT_ASSERT(Blenx > 1);

    n_poly_init(Beval);
    n_poly_factor_init(local_fac);
    n_bpoly_init(monicB);
    n_tpoly_init(lift_fac);

    starts = (slong *) flint_malloc(Blenx*sizeof(slong));

    /* init done */

    n_bpoly_mod_make_primitive(c, B, ctx);

    alpha = 0;
    goto got_alpha;

next_alpha:

    alpha++;
    if (!allow_shift || alpha >= ctx.n)
    {
        success = 0;
        goto cleanup;
    }

got_alpha:

    n_bpoly_mod_eval(Beval, B, alpha, ctx);

    /* if killed leading coeff, get new alpha */
    if (Beval->length != Blenx)
        goto next_alpha;

    n_poly_mod_factor(local_fac, Beval, ctx);

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
        n_tpoly_fit_length(F, 1);
        F->length = 1;
        n_bpoly_swap(F->coeffs + 0, B);
        success = 1;
        goto cleanup;
	}

    n_bpoly_mod_taylor_shift_var1(B, alpha, ctx);

    Bleny = 0;
    for (i = 0; i < B->length; i++)
        Bleny = FLINT_MAX(Bleny, (B->coeffs + i)->length);

    lift_order = Bleny + 10;
    final_order = Bleny;

    n_bpoly_mod_make_monic_series(monicB, B, lift_order, ctx);

    n_tpoly_fit_length(lift_fac, r);
    lift_fac->length = r;

    for (i = 0; i < r; i++)
    {
        FLINT_ASSERT(local_fac->poly[i].length > 1);
        FLINT_ASSERT(local_fac->poly[i].coeffs[local_fac->poly[i].length - 1] == 1);
    }

timeit_start(timer);
    {
        slong * link = flint_malloc((2*r - 2)*sizeof(slong));
        n_bpoly_struct * v = flint_malloc(2*(2*r - 2) * sizeof(n_bpoly_struct));
        n_bpoly_struct * w = v + (2*r - 2);

        for (i = 0; i < 2*r - 2; i++)
        {
            n_bpoly_init(v + i);
            n_bpoly_init(w + i);
        }

        _hensel_build_tree(link, v, w, local_fac->poly, r, ctx);

        {
            slong * e, n = FLINT_CLOG2(lift_order) + 1;
            e = flint_malloc(n * sizeof(slong));
            for (e[i = 0] = lift_order; e[i] > 1; i++)
                e[i + 1] = (e[i] + 1) / 2;

            for (i--; i > 0; i--)
                _hensel_lift_tree(0, link, v, w, monicB, 2*r-4, e[i+1], e[i]-e[i+1], ctx);

            if (lift_order > 1)
                _hensel_lift_tree(1, link, v, w, monicB, 2*r-4, e[i+1], e[i]-e[i+1], ctx);

            for (i = 0; i < 2*r - 2; i++)
            { 
                if (link[i] < 0)
                {
/*
flint_printf("lifted factor: ");
n_bpoly_print_pretty(v + i, "x", "y");
flint_printf("\n");
*/
                }
            }

        }

        for (i = 0; i < 2*r - 2; i++)
        {
            n_bpoly_clear(v + i);
            n_bpoly_clear(w + i);
        }
        flint_free(link);
        flint_free(v);
    }
timeit_stop(timer);
flint_printf("   tree: %wd\n", timer->wall);


timeit_start(timer);
    _lift_quintic(lift_fac->coeffs, local_fac->poly, r, monicB, lift_order, ctx);
timeit_stop(timer);
flint_printf("quintic: %wd\n", timer->wall);

    for (i = 0; i < Blenx; i++)
        starts[i] = Bleny;

    nmod_mat_init(N, r, r, ctx.n);
    for (i = 0; i < r; i++)
        nmod_mat_entry(N, i, i) = 1;

    _recombine_lattice(N, lift_fac->coeffs, r, lift_order, starts, B, ctx);

    if (!nmod_mat_is_reduced(N))
    {
flint_printf("it wasn't reduced\n");
usleep(1000000);
        nmod_mat_clear(N);
        nmod_mat_init(N, r, r, ctx.n);
        for (i = 0; i < r; i++)
            nmod_mat_entry(N, i, i) = 1;
    }

    F->length = 0;
    _recombine_zassenhaus(F, nmod_neg(alpha, ctx), N,
                                     lift_fac->coeffs, r, final_order, B, ctx);

    nmod_mat_clear(N);
    success = 1;

cleanup:

    n_poly_clear(Beval);
    n_poly_factor_clear(local_fac);
    n_bpoly_clear(monicB);
    n_tpoly_clear(lift_fac);
    flint_free(starts);

    return success;
}

