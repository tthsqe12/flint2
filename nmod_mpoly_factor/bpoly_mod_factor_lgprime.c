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


typedef struct {
    slong r; /* number of local factors */
    n_bpoly_t Btilde;
    n_bpoly_struct * newBitilde;
    fq_nmod_poly_struct * d;
    fq_nmod_poly_struct * Bitilde;
} nmod_lgprime_info_struct;

typedef nmod_lgprime_info_struct nmod_lgprime_info_t[1];

void nmod_lgprime_info_init(
    nmod_lgprime_info_t I,
    slong r,
    const fq_nmod_ctx_t ectx)
{
    slong i;

    FLINT_ASSERT(r >= 2);

    I->r = r;

    n_bpoly_init(I->Btilde);

    I->newBitilde = (n_bpoly_struct *) flint_malloc(r*sizeof(n_bpoly_struct));
    I->d          = (fq_nmod_poly_struct *) flint_malloc(r*sizeof(fq_nmod_poly_struct));
    I->Bitilde    = (fq_nmod_poly_struct *) flint_malloc(r*sizeof(fq_nmod_poly_struct));

    for (i = 0; i < r; i++)
    {
        n_bpoly_init(I->newBitilde + i);
        fq_nmod_poly_init(I->d + i, ectx);
        fq_nmod_poly_init(I->Bitilde + i, ectx);
    }
}

void nmod_lgprime_info_clear(
    nmod_lgprime_info_t I,
    const fq_nmod_ctx_t ectx)
{
    slong i;

    n_bpoly_clear(I->Btilde);

    for (i = 0; i < I->r; i++)
    {
        n_bpoly_clear(I->newBitilde + i);
        fq_nmod_poly_clear(I->d + i, ectx);
        fq_nmod_poly_clear(I->Bitilde + i, ectx);
    }

    flint_free(I->newBitilde);
    flint_free(I->d);
    flint_free(I->Bitilde);
}

int nmod_lgprime_info_disolve(
    nmod_lgprime_info_t I,
    const fq_nmod_ctx_t ectx)
{
    return fq_nmod_partial_fraction_coeffs(I->r, I->d, I->Bitilde, ectx);
}

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
    const n_poly_t alpha,
    const n_bpoly_struct * loc_fac,
    slong r,
    n_bpoly_t f,
    nmod_t ctx)
{
    n_poly_t g;
    n_bpoly_t Q, R, trymez;
    n_bpoly_t tryme, trymet;
    n_poly_t leadf;
    slong kk, *used_arr, *sub_arr;

    used_arr = (slong *) flint_calloc(2*r, sizeof(slong));
    sub_arr  = used_arr + r;

    n_poly_init(g);
    n_bpoly_init(Q);
    n_bpoly_init(R);
    n_bpoly_init(trymez);
    n_bpoly_init(tryme);
    n_bpoly_init(trymet);
    n_poly_init(leadf);

    FLINT_ASSERT(f->length > 0);
    n_poly_set(leadf, f->coeffs + f->length - 1);

    for (kk = 1; kk < r; kk++)
    {
        slong count = 0, indx = kk - 1, l;

        for(l = 0; l < kk; l++)
            sub_arr[l] = l;

        sub_arr[indx]--;
        while ((indx >= 0))
        {
            sub_arr[indx] = sub_arr[indx] + 1;

            for (l = indx + 1; l < kk; l++)
                sub_arr[l] = sub_arr[l - 1] + 1;

            if (sub_arr[kk - 1] > r - 1)
                indx--;
            else
            {
                for (l = 0; l < kk; l++)
                {
                    if (used_arr[sub_arr[l]] == 1)
                        break;
                }

                n_bpoly_set_poly_var1(tryme, leadf);
                for (l = 0; l < kk; l++)
                {
                    n_bpoly_mod_mulmod_poly(trymet, tryme, loc_fac + sub_arr[l], alpha, ctx);
                    n_bpoly_swap(trymet, tryme);
                }

                n_bpoly_set(trymez, tryme);
                n_bpoly_mod_make_primitive(g, trymez, ctx);

                if (n_bpoly_mod_divides(Q, f, trymez, ctx))
                {
                    n_tpoly_fit_length(F, F->length + 1);

                    n_bpoly_swap(F->coeffs + F->length, trymez);
                    F->length++;

                    for (l = 0; l < kk; l++)
                    {
                        used_arr[sub_arr[l]] = 1;
                        count++;
                    }

                    n_bpoly_set(f, Q);
                    FLINT_ASSERT(f->length > 0);
                    n_poly_set(leadf, f->coeffs + f->length - 1);

                 /* If r - count = kk then the rest are irreducible.  
                    TODO: Add a test for that case */
                }

                indx = kk - 1;
            }
        }
    }

    {
        slong test = 0;

        for (kk = 0; kk < r; kk++)
            test = test + used_arr[kk];

        if (test == 0)
        {
            n_tpoly_fit_length(F, F->length + 1);
            n_bpoly_swap(F->coeffs + F->length, f);
            F->length++;
        }
    }

    n_poly_clear(g);
    n_bpoly_clear(Q);
    n_bpoly_clear(R);
    n_bpoly_clear(trymez);
    n_bpoly_clear(tryme);
    n_bpoly_clear(trymet);
    n_poly_clear(leadf);

    flint_free(used_arr);
}


void n_bpoly_mod_factor_lgprime(
    n_poly_t c,
    n_tpoly_t F,
    n_bpoly_t B,
    nmod_t ctx)
{
    slong i, j, k, w;
    fq_nmod_poly_t Beval;
    fq_nmod_poly_factor_t Bevalfac;
    fq_nmod_t Bfaclc;
    slong Blengthx, Blengthy;
    nmod_lgprime_info_t I;
    n_bpoly_t tp, tp1, error;
    fq_nmod_poly_t ss, tt;
	slong deg;
    fq_nmod_ctx_t ectx;
    fmpz_t P;
    n_poly_t finalmpow, mpow;
    n_poly_t mock;
    nmod_poly_t mock2;
    n_poly_t g;

    deg = 2;
    fmpz_init_set_ui(P, ctx.n);
    fq_nmod_ctx_init(ectx, P, deg, "y");

    n_poly_init(g);
    n_poly_init(mpow);
    n_poly_init(finalmpow);
    fq_nmod_poly_init(Beval, ectx);
    fq_nmod_poly_factor_init(Bevalfac, ectx);
    fq_nmod_init(Bfaclc, ectx);
    nmod_lgprime_info_init(I, 2, ectx);
    fq_nmod_poly_init(ss, ectx);
    fq_nmod_poly_init(tt, ectx);
    n_bpoly_init(tp);
    n_bpoly_init(tp1);
    n_bpoly_init(error);

    F->length = 0;

    Blengthx = B->length;
    FLINT_ASSERT(Blengthx > 1);

    goto got_alpha;

next_alpha:

    deg++;

	fq_nmod_ctx_clear(ectx);
	fq_nmod_ctx_init(ectx, P, deg, "y");

got_alpha:

    n_bpoly_eval_fq_nmod_poly(Beval, ectx, B);

    /* if killed leading coeff, get new alpha */
    if (Beval->length != Blengthx)
        goto next_alpha;

    Bevalfac->num = 0;
    fq_nmod_poly_factor(Bevalfac, Bfaclc, Beval, ectx);

    /* if multiple factors, get new alpha */
    for (i = 0; i < Bevalfac->num; i++)
    {
        if (Bevalfac->exp[i] != 1)
            goto next_alpha;
    }

	/* if one factor, A is irreducible */
	if (Bevalfac->num < 2)
	{
        n_tpoly_fit_length(F, 1);
        F->length = 1;
        n_bpoly_swap(F->coeffs + 0, B);
        goto cleanup;
	}

    Blengthy = 0;
    for (i = 0; i < B->length; i++)
        Blengthy = FLINT_MAX(Blengthy, (B->coeffs + i)->length);

    n_poly_one(finalmpow);
    w = 0;
    while (finalmpow->length <= Blengthy)
    {
        w++;
        n_poly_mock(mock, ectx->modulus);
        n_poly_mod_mul(finalmpow, finalmpow, mock, ctx);
    }

    nmod_lgprime_info_clear(I, ectx);
    nmod_lgprime_info_init(I, Bevalfac->num, ectx);

    n_bpoly_set(I->Btilde, B);
    n_bpoly_mod_make_monic_mod(I->Btilde, finalmpow, ctx);

    for (i = 0; i < I->r; i++)
    {
        fq_nmod_poly_make_monic(I->Bitilde + i, Bevalfac->poly + i, ectx);
        n_bpoly_set_fq_nmod_poly(I->newBitilde + i, I->Bitilde + i, ectx);
    }

    nmod_lgprime_info_disolve(I, ectx);

    FLINT_ASSERT(I->r > 1);

    n_bpoly_mod_mulmod_poly(tp, I->newBitilde + 0, I->newBitilde + 1, finalmpow, ctx);
    for (i = 2; i < I->r; i++)
    {
        n_bpoly_mod_mulmod_poly(tp1, tp, I->newBitilde + i, finalmpow, ctx);
        n_bpoly_swap(tp1, tp);
    }

    n_bpoly_mod_sub(error, I->Btilde, tp, ctx);

    n_poly_one(mpow);

    for (j = 1; j < w; j++)
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

        for (i = 0; i < I->r; i++)
        {
            fq_nmod_poly_mul(tt, ss, I->d + i, ectx);
            fq_nmod_poly_rem(tt, tt, I->Bitilde + i, ectx);
            n_bpoly_add_fq_nmod_poly_mul(I->newBitilde + i, tt, mpow, ctx);
        }

        n_bpoly_mod_mulmod_poly(tp, I->newBitilde + 0, I->newBitilde + 1, finalmpow, ctx);
        for (i = 2; i < I->r; i++)
        {
            n_bpoly_mod_mulmod_poly(tp1, tp, I->newBitilde + i, finalmpow, ctx);
            n_bpoly_swap(tp1, tp);
        }
        n_bpoly_mod_sub(error, I->Btilde, tp, ctx);
    }

    _recombine_zassenhaus(F, finalmpow, I->newBitilde, I->r, B, ctx);

cleanup:

    n_poly_clear(g);
    n_poly_clear(mpow);
    n_poly_clear(finalmpow);

    fq_nmod_poly_clear(Beval, ectx);
    fq_nmod_poly_factor_clear(Bevalfac, ectx);
    fq_nmod_clear(Bfaclc, ectx);
    nmod_lgprime_info_clear(I, ectx);
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
    return;
}
