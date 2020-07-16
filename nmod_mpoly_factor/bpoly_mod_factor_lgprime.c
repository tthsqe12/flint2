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


void n_bpoly_mod_mulmod_poly(
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

            n_bpoly_mod_mulmod_poly(t1, loc_fac + i, loc_fac_ + j, finalmpow, ctx);
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
                    n_bpoly_mod_mulmod_poly(t2, t1, loc_fac + idx[i], finalmpow, ctx);
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


void n_bpoly_mod_factor_lgprime(
    n_poly_t c,
    n_tpoly_t F,
    n_bpoly_t B,
    nmod_t ctx)
{
    slong i, j, k, finalpow;
    fq_nmod_poly_t Beval;
    fq_nmod_poly_factor_t Bevalfac;
    fq_nmod_t Bfaclc;
    slong Blenx = B->length;
    slong Bleny;
    n_bpoly_t IBtilde;
    n_tpoly_t InewBitilde;
    fq_nmod_bpoly_t Id;
    fq_nmod_bpoly_t IBitilde;
    n_bpoly_t tp, tp1, error;
    fq_nmod_poly_t ss, tt;
	slong deg, r;
    fq_nmod_ctx_t ectx;
    fmpz_t P;
    n_poly_t finalmpow, mpow;
    n_poly_t mock;
    nmod_poly_t mock2;
    n_poly_t g;
    nmod_mat_t Ntr;
timeit_t timer;

flint_printf("n_bpoly_mod_factor_lgprime called\n");

    deg = 2;
    fmpz_init_set_ui(P, ctx.n);
    fq_nmod_ctx_init(ectx, P, deg, "y");

    FLINT_ASSERT(Blenx > 1);

    n_poly_init(g);
    n_poly_init(mpow);
    n_poly_init(finalmpow);
    fq_nmod_poly_init(Beval, ectx);
    fq_nmod_poly_factor_init(Bevalfac, ectx);
    fq_nmod_init(Bfaclc, ectx);

    n_bpoly_init(IBtilde);
    n_tpoly_init(InewBitilde);
    fq_nmod_bpoly_init(Id, ectx);
    fq_nmod_bpoly_init(IBitilde, ectx);

    fq_nmod_poly_init(ss, ectx);
    fq_nmod_poly_init(tt, ectx);
    n_bpoly_init(tp);
    n_bpoly_init(tp1);
    n_bpoly_init(error);

    F->length = 0;

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

    Bevalfac->num = 0;
    fq_nmod_poly_factor(Bevalfac, Bfaclc, Beval, ectx);

    r = Bevalfac->num;
flint_printf("r: %wd\n", r);

    /* if multiple factors, get new alpha */
    for (i = 0; i < r; i++)
    {
        if (Bevalfac->exp[i] != 1)
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

    finalpow = (Bleny - 1 + deg)/deg;

    n_poly_mock(mock, ectx->modulus);
    n_poly_mod_pow(finalmpow, mock, finalpow, ctx);

    n_tpoly_fit_length(InewBitilde, r);
    InewBitilde->length = r;
    fq_nmod_bpoly_fit_length(Id, r, ectx);
    Id->length = r;
    fq_nmod_bpoly_fit_length(IBitilde, r, ectx);
    IBitilde->length = r;

    n_bpoly_set(IBtilde, B);
    n_bpoly_mod_make_monic_mod(IBtilde, finalmpow, ctx);

    for (i = 0; i < r; i++)
    {
        fq_nmod_poly_make_monic(IBitilde->coeffs + i, Bevalfac->poly + i, ectx);
        n_bpoly_set_fq_nmod_poly(InewBitilde->coeffs + i, IBitilde->coeffs + i, ectx);
    }

    fq_nmod_partial_fraction_coeffs(r, Id->coeffs, IBitilde->coeffs, ectx);

    n_bpoly_mod_mulmod_poly(tp, InewBitilde->coeffs + 0, InewBitilde->coeffs + 1, finalmpow, ctx);
    for (i = 2; i < r; i++)
    {
        n_bpoly_mod_mulmod_poly(tp1, tp, InewBitilde->coeffs + i, finalmpow, ctx);
        n_bpoly_swap(tp1, tp);
    }

    n_bpoly_mod_sub(error, IBtilde, tp, ctx);

    n_poly_one(mpow);

timeit_start(timer);

    for (j = 1; j < finalpow; j++)
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
            fq_nmod_poly_rem(tt, tt, IBitilde->coeffs + i, ectx);
            n_bpoly_add_fq_nmod_poly_mul(InewBitilde->coeffs + i, tt, mpow, ctx);
        }

        n_bpoly_mod_mulmod_poly(tp, InewBitilde->coeffs + 0, InewBitilde->coeffs + 1, finalmpow, ctx);
        for (i = 2; i < r; i++)
        {
            n_bpoly_mod_mulmod_poly(tp1, tp, InewBitilde->coeffs + i, finalmpow, ctx);
            n_bpoly_swap(tp1, tp);
        }
        n_bpoly_mod_sub(error, IBtilde, tp, ctx);
    }

timeit_stop(timer);
flint_printf("lift: %wd\n", timer->wall);

timeit_start(timer);

    nmod_mat_init(Ntr, r, r, ctx.n);
    for (i = 0; i < r; i++)
        nmod_mat_entry(Ntr, i, i) = 1;

    F->length = 0;
    _recombine_zassenhaus(F, finalmpow, Ntr, InewBitilde->coeffs, r, B, ctx);

timeit_stop(timer);
flint_printf("zass: %wd\n", timer->wall);

    nmod_mat_clear(Ntr);

cleanup:

    n_poly_clear(g);
    n_poly_clear(mpow);
    n_poly_clear(finalmpow);

    fq_nmod_poly_clear(Beval, ectx);
    fq_nmod_poly_factor_clear(Bevalfac, ectx);
    fq_nmod_clear(Bfaclc, ectx);
    n_bpoly_clear(IBtilde);
    n_tpoly_clear(InewBitilde);
    fq_nmod_bpoly_clear(Id, ectx);
    fq_nmod_bpoly_clear(IBitilde, ectx);

    fq_nmod_poly_clear(ss, ectx);
    fq_nmod_poly_clear(tt, ectx);
    n_bpoly_clear(tp);
    n_bpoly_clear(tp1);
    n_bpoly_clear(error);

    fq_nmod_ctx_clear(ectx);
/*
flint_printf("n_bpoly_mod_factor_lgprime returning\n");
for (i = 0; i < F->length; i++)
{
flint_printf("F[%wd]: ", i); n_bpoly_print_pretty(F->coeffs + i, "x", "y"); flint_printf("\n");
}
*/

flint_printf("n_bpoly_mod_factor_lgprime called\n");

    return;
}
