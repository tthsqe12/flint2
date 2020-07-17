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


void n_bpoly_eval_fq_nmod_poly(
    fq_nmod_poly_t A,
    const fq_nmod_ctx_t ectx,
    const n_bpoly_t B)
{
    slong i;
    n_poly_t t;
    n_poly_t mock;
    nmod_poly_t mock2;

    n_poly_init(t);

    fq_nmod_poly_zero(A, ectx);
    for (i = B->length - 1; i >= 0; i--)
    {
        n_poly_mock(mock, ectx->modulus);
        n_poly_mod_rem(t, B->coeffs + i, mock, ectx->modulus->mod);
        nmod_poly_mock(mock2, t, ectx->modulus->mod);
        fq_nmod_poly_set_coeff(A, i, mock2, ectx);
    }

    n_poly_clear(t);
}


void n_bpoly_mod_make_monic_mod(n_bpoly_t A, n_poly_t mk, nmod_t mod)
{
    slong i;
    int success;
    n_poly_t t, lcinv;

    FLINT_ASSERT(A->length > 0);
    FLINT_ASSERT(n_bpoly_mod_is_canonical(A, mod));

    n_poly_init(t);
    n_poly_init(lcinv);
    success = n_poly_mod_invmod(lcinv, A->coeffs + A->length - 1, mk, mod);
    FLINT_ASSERT(success);

    for (i = 0; i < A->length; i++)
    {
        n_poly_mod_mulmod(t, A->coeffs + i, lcinv, mk, mod);
        n_poly_swap(A->coeffs + i, t);
    }

    n_poly_clear(t);
    n_poly_clear(lcinv);
}


void n_bpoly_set_fq_nmod_poly(
    n_bpoly_t A,
    const fq_nmod_poly_t B,
    const fq_nmod_ctx_t ectx)
{
    slong i;

    n_bpoly_fit_length(A, B->length);
    A->length = B->length;

    for (i = 0; i < B->length; i++)
        n_poly_set_nmod_poly(A->coeffs + i, B->coeffs + i);
}


void n_bpoly_mod_mul_mod_poly(
    n_bpoly_t A,
    const n_bpoly_t B,
    const n_bpoly_t C,
    const n_poly_t m,
    nmod_t mod)
{
    slong i, j;
    n_poly_t t;

    FLINT_ASSERT(A != B);
    FLINT_ASSERT(A != C);

    n_poly_init(t);

    n_bpoly_fit_length(A, B->length + C->length - 1);
    for (i = 0; i < B->length + C->length - 1; i++)
        n_poly_zero(A->coeffs + i);

    for (i = 0; i < B->length; i++)
    {
        for (j = 0; j < C->length; j++)
        {
            n_poly_mod_mul(t, B->coeffs + i, C->coeffs + j, mod);
            n_poly_mod_add(A->coeffs + i + j, A->coeffs + i + j, t, mod);
            n_poly_mod_rem(A->coeffs + i + j, A->coeffs + i + j, m, mod);
        }
    }

    A->length = B->length + C->length - 1;

    while (A->length > 0 && n_poly_is_zero(A->coeffs + A->length - 1))
        A->length--;

    n_poly_clear(t);
}

/* division in ((Z/nZ)[y]/m(y))[x] */
void n_bpoly_mod_divrem_mod_poly(
    n_bpoly_t Q,
    n_bpoly_t R,
    const n_bpoly_t A,
    const n_bpoly_t B,
    const n_poly_t m,
    nmod_t ctx)
{
    int success;
    slong i, qoff;
    n_poly_t q, t, Binv;

    FLINT_ASSERT(R != A);
    FLINT_ASSERT(R != B);
    FLINT_ASSERT(Q != A);
    FLINT_ASSERT(Q != B);
    FLINT_ASSERT(B->length > 0);

    n_poly_init(q);
    n_poly_init(t);
    n_poly_init(Binv);

    n_bpoly_set(R, A);

    Q->length = 0;

    success = n_poly_mod_invmod(Binv, B->coeffs + B->length - 1, m, ctx);
    FLINT_ASSERT(success);

    while (R->length >= B->length)
    {
        n_poly_mod_mulmod(q, R->coeffs + R->length - 1, Binv, m, ctx);

        for (i = 0; i < B->length; i++)
        {
            n_poly_mod_mulmod(t, B->coeffs + i, q, m, ctx);
            n_poly_mod_sub(R->coeffs + i + R->length - B->length,
                            R->coeffs + i + R->length - B->length, t, ctx);
        }

        qoff = R->length - B->length;

        FLINT_ASSERT(qoff >= 0);
        if (qoff >= Q->length)
        {
            n_bpoly_fit_length(Q, qoff + 1);
            for (i = Q->length; i <= qoff; i++)
                n_poly_zero(Q->coeffs + i);
            Q->length = qoff + 1;
        }

        n_poly_set(Q->coeffs + qoff, q);

        FLINT_ASSERT(n_poly_is_zero(R->coeffs + R->length - 1));

        while (R->length > 0 && n_poly_is_zero(R->coeffs + R->length - 1))
            R->length--;
    }

    n_poly_clear(q);
    n_poly_clear(t);
    n_poly_clear(Binv);
}


void n_bpoly_add_fq_nmod_poly_mul(
    n_bpoly_t A,
    const fq_nmod_poly_t B,
    const n_poly_t mpow,
    nmod_t mod)
{
    slong i;
    n_poly_t t, mock;

    n_poly_init(t);

    FLINT_ASSERT(A->length > B->length);

    for (i = 0; i < B->length; i++)
    {
        n_poly_mock(mock, B->coeffs + i);
        n_poly_mod_mul(t, mock, mpow, mod);
        n_poly_mod_add(A->coeffs + i, A->coeffs + i, t, mod);
    }

    n_poly_clear(t);
}

static void _recombine_zassenhaus(
    n_tpoly_t F,
    const n_poly_t finalmpow,
    const nmod_mat_t N,
    const n_bpoly_struct * loc_fac_,
    slong r,
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

            n_bpoly_mod_mul_mod_poly(t1, loc_fac + i, loc_fac_ + j, finalmpow, ctx);
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
                    n_bpoly_mod_mul_mod_poly(t2, t1, loc_fac + idx[i], finalmpow, ctx);
                    n_bpoly_swap(t1, t2);
                }
            }

            n_bpoly_mod_make_primitive(g, t1, ctx);
            if (n_bpoly_mod_divides(Q, f, t1, ctx))
            {
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
    fq_nmod_poly_struct * local_facs,
    slong r,
    const n_bpoly_t IBtilde,
    slong liftpow,
    const n_poly_t liftmpow,
    const fq_nmod_ctx_t ectx)
{
    const nmod_t ctx = ectx->modulus->mod;
    slong i, j, k;
    fq_nmod_bpoly_t Id;
    fq_nmod_poly_t ss, tt;
    n_bpoly_t tp, tp1, error;
    n_poly_t mpow;
    n_poly_t mock;
    nmod_poly_t mock2;

    fq_nmod_bpoly_init(Id, ectx);
    fq_nmod_poly_init(ss, ectx);
    fq_nmod_poly_init(tt, ectx);
    n_bpoly_init(tp);
    n_bpoly_init(tp1);
    n_bpoly_init(error);
    n_poly_init(mpow);

    fq_nmod_bpoly_fit_length(Id, r, ectx);
    Id->length = r;

    for (i = 0; i < r; i++)
    {
        n_bpoly_set_fq_nmod_poly(lift_facs + i, local_facs + i, ectx);
    }

    fq_nmod_partial_fraction_coeffs(r, Id->coeffs, local_facs, ectx);

    n_bpoly_mod_mul_mod_poly(tp, lift_facs + 0, lift_facs + 1, liftmpow, ctx);
    for (i = 2; i < r; i++)
    {
        n_bpoly_mod_mul_mod_poly(tp1, tp, lift_facs + i, liftmpow, ctx);
        n_bpoly_swap(tp1, tp);
    }

    n_bpoly_mod_sub(error, IBtilde, tp, ctx);

    n_poly_one(mpow);

    for (j = 1; j < liftpow; j++)
    {
        n_poly_mock(mock, ectx->modulus);
        n_poly_mod_mul(mpow, mpow, mock, ctx);

        fq_nmod_poly_zero(ss, ectx);

        for (k = error->length - 1; k >= 0; k--)
        {
            n_poly_t q, r;
            n_poly_init(q);
            n_poly_init(r);
            n_poly_mod_divrem(q, r, error->coeffs + k, mpow, ctx);
            FLINT_ASSERT(n_poly_is_zero(r));
            n_poly_mock(mock, ectx->modulus);
            n_poly_mod_rem(r, q, mock, ctx);
            nmod_poly_mock(mock2, r, ctx);
            fq_nmod_poly_set_coeff(ss, k, mock2, ectx);
            n_poly_clear(q);
            n_poly_clear(r);
        }

        for (i = 0; i < r; i++)
        {
            fq_nmod_poly_mul(tt, ss, Id->coeffs + i, ectx);
            fq_nmod_poly_rem(tt, tt, local_facs + i, ectx);
            n_bpoly_add_fq_nmod_poly_mul(lift_facs + i, tt, mpow, ctx);
        }

        n_bpoly_mod_mul_mod_poly(tp, lift_facs + 0, lift_facs + 1, liftmpow, ctx);
        for (i = 2; i < r; i++)
        {
            n_bpoly_mod_mul_mod_poly(tp1, tp, lift_facs + i, liftmpow, ctx);
            n_bpoly_swap(tp1, tp);
        }
        n_bpoly_mod_sub(error, IBtilde, tp, ctx);
    }

    fq_nmod_bpoly_clear(Id, ectx);
    fq_nmod_poly_clear(ss, ectx);
    fq_nmod_poly_clear(tt, ectx);
    n_bpoly_clear(tp);
    n_bpoly_clear(tp1);
    n_bpoly_clear(error);
    n_poly_clear(mpow);
}


static void _recombine_lattice(
    nmod_mat_t N,
    const n_bpoly_struct * g,
    slong r,
    const n_poly_t lift_alpha_pow,
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
    slong lift_order = lift_alpha_pow->length - 1;

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
        n_bpoly_mod_divrem_mod_poly(Q, R, f, g + i, lift_alpha_pow, ctx);
        FLINT_ASSERT(R->length == 0);
        n_bpoly_mod_derivative(R, g + i, ctx);
        n_bpoly_mod_mul_mod_poly(ld + i, Q, R, lift_alpha_pow, ctx);
    }

    for (k = 0; k + 1 < f->length; k++)
    {
        slong d = nmod_mat_nrows(N);

        if (d < 2)
            break;

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
/*
flint_printf("k = %wd, N:\n", k); nmod_mat_print_pretty(N); flint_printf("\n");
*/
        nmod_mat_clear(M);
        nmod_mat_clear(T1);
        nmod_mat_clear(T2);

        if (nmod_mat_is_reduced(N))
            break;
    }

    flint_free(trow);
    n_bpoly_clear(Q);
    n_bpoly_clear(R);
    n_bpoly_clear(dg);
    for (i = 0; i < r; i++)
        n_bpoly_clear(ld + i);
    flint_free(ld);
}

void n_bpoly_mod_factor_lgprime(
    n_poly_t c,
    n_tpoly_t F,
    n_bpoly_t B,
    nmod_t ctx)
{
    slong i, r, deg;
    fq_nmod_poly_t Beval;
    fq_nmod_poly_factor_t local_fac;
    fq_nmod_t Blc;
    slong Blenx = B->length;
    slong Bleny;
    n_bpoly_t monicB;
    n_tpoly_t lift_fac;
    fq_nmod_bpoly_t Id;
    fq_nmod_ctx_t ectx;
    fmpz_t P;
    slong final_pow, lift_pow;
    n_poly_t final_alpha_pow, lift_alpha_pow;
    n_poly_t mock;
    nmod_mat_t N;
    slong * starts;

timeit_t timer;

    deg = 2;
    fmpz_init_set_ui(P, ctx.n);
    fq_nmod_ctx_init(ectx, P, deg, "y");

    FLINT_ASSERT(Blenx > 1);

    n_poly_init(final_alpha_pow);
    n_poly_init(lift_alpha_pow);

    fq_nmod_poly_init(Beval, ectx);
    fq_nmod_poly_factor_init(local_fac, ectx);
    fq_nmod_init(Blc, ectx);

    n_bpoly_init(monicB);
    n_tpoly_init(lift_fac);
    fq_nmod_bpoly_init(Id, ectx);

    starts = (slong *) flint_malloc(Blenx*sizeof(slong));

    /* init done */

    n_bpoly_mod_make_primitive(c, B, ctx);

    goto got_alpha;

next_alpha:

    deg++;

	fq_nmod_ctx_clear(ectx);
	fq_nmod_ctx_init(ectx, P, deg, "y");

got_alpha:

    n_bpoly_eval_fq_nmod_poly(Beval, ectx, B);

    /* if killed leading coeff, get new alpha */
    if (Beval->length != Blenx)
        goto next_alpha;

    local_fac->num = 0;
    fq_nmod_poly_factor(local_fac, Blc, Beval, ectx);

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
        goto cleanup;
	}

    Bleny = 0;
    for (i = 0; i < B->length; i++)
        Bleny = FLINT_MAX(Bleny, (B->coeffs + i)->length);

    final_pow = (Bleny - 1 + deg)/deg;
    lift_pow = final_pow + 4;

    n_poly_mock(mock, ectx->modulus);
    n_poly_mod_pow(final_alpha_pow, mock, final_pow, ctx);
    n_poly_mod_pow(lift_alpha_pow, mock, lift_pow, ctx);

    n_bpoly_set(monicB, B);
    n_bpoly_mod_make_monic_mod(monicB, lift_alpha_pow, ctx);

    n_tpoly_fit_length(lift_fac, r);
    lift_fac->length = r;

    for (i = 0; i < r; i++)
    {
        FLINT_ASSERT(local_fac->poly[i].length > 1);
        FLINT_ASSERT(fq_nmod_is_one(local_fac->poly[i].coeffs + local_fac->poly[i].length - 1, ectx));
    }

timeit_start(timer);
    _lift_quintic(lift_fac->coeffs, local_fac->poly, r, monicB,
                                               lift_pow, lift_alpha_pow, ectx);
timeit_stop(timer);
flint_printf("lift: %wd\n", timer->wall);

    for (i = 0; i < Blenx; i++)
        starts[i] = Bleny;

    nmod_mat_init(N, r, r, ctx.n);
    for (i = 0; i < r; i++)
        nmod_mat_entry(N, i, i) = 1;

timeit_start(timer);
    _recombine_lattice(N, lift_fac->coeffs, r, lift_alpha_pow, starts, B, ctx);
timeit_stop(timer);
flint_printf("latt: %wd\n", timer->wall);

    if (!nmod_mat_is_reduced(N))
    {
flint_printf("it wasn't reduced\n");
usleep(1000000);
        nmod_mat_clear(N);
        nmod_mat_init(N, r, r, ctx.n);
        for (i = 0; i < r; i++)
            nmod_mat_entry(N, i, i) = 1;
    }

timeit_start(timer);
    F->length = 0;
    _recombine_zassenhaus(F, final_alpha_pow, N, lift_fac->coeffs, r, B, ctx);
timeit_stop(timer);
flint_printf("zass: %wd\n", timer->wall);

    nmod_mat_clear(N);

cleanup:

    n_poly_clear(final_alpha_pow);
    n_poly_clear(lift_alpha_pow);

    fq_nmod_poly_clear(Beval, ectx);
    fq_nmod_poly_factor_clear(local_fac, ectx);
    fq_nmod_clear(Blc, ectx);
    n_bpoly_clear(monicB);
    n_tpoly_clear(lift_fac);
    fq_nmod_bpoly_clear(Id, ectx);

    fq_nmod_ctx_clear(ectx);

    return;
}
