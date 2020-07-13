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

void n_bpoly_mod_make_monic_series(n_bpoly_t A, slong order, nmod_t mod)
{
    slong i;
    n_poly_t t, lcinv;

    FLINT_ASSERT(A->length > 0);
    FLINT_ASSERT(n_bpoly_mod_is_canonical(A, mod));

    n_poly_init(t);
    n_poly_init(lcinv);
    n_poly_mod_inv_series(lcinv, A->coeffs + A->length - 1, order, mod);

    for (i = 0; i < A->length; i++)
    {
        n_poly_mod_mullow(t, A->coeffs + i, lcinv, order, mod);
        n_poly_swap(A->coeffs + i, t);
    }

    n_poly_clear(t);
    n_poly_clear(lcinv);
}

typedef struct {
    slong r; /* number of local factors */
    nmod_t modulus;
    n_bpoly_t Btilde;
    n_bpoly_struct * newBitilde;
    n_poly_struct * P;
    n_poly_struct * d;
    n_poly_struct * Bitilde;
} nmod_bpoly_info_struct;

typedef nmod_bpoly_info_struct nmod_bpoly_info_t[1];

void nmod_bpoly_info_init(nmod_bpoly_info_t I, slong r)
{
    slong i;

    FLINT_ASSERT(r >= 2);

    I->r = r;

    n_bpoly_init(I->Btilde);

    I->newBitilde = (n_bpoly_struct *) flint_malloc(r*sizeof(n_bpoly_struct));
    I->P          = (n_poly_struct *) flint_malloc(r*sizeof(n_poly_struct));
    I->d          = (n_poly_struct *) flint_malloc(r*sizeof(n_poly_struct));
    I->Bitilde    = (n_poly_struct *) flint_malloc(r*sizeof(n_poly_struct));

    for (i = 0; i < r; i++)
    {
        n_bpoly_init(I->newBitilde + i);
        n_poly_init(I->P + i);
        n_poly_init(I->d + i);
        n_poly_init(I->Bitilde + i);
    }
}

void nmod_bpoly_info_clear(nmod_bpoly_info_t I)
{
    slong i;

    n_bpoly_clear(I->Btilde);

    for (i = 0; i < I->r; i++)
    {
        n_bpoly_clear(I->newBitilde + i);
        n_poly_clear(I->P + i);
        n_poly_clear(I->d + i);
        n_poly_clear(I->Bitilde + i);
    }

    flint_free(I->newBitilde);
    flint_free(I->P);
    flint_free(I->d);
    flint_free(I->Bitilde);
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

int nmod_bpoly_info_disolve(nmod_bpoly_info_t I, nmod_t ctx)
{
    return nmod_partial_fraction_coeffs(I->r, I->d, I->Bitilde, ctx);
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


int n_bpoly_mod_factor_smprime(
    n_poly_t c,
    n_tpoly_t F,
    n_bpoly_t B,
    int allow_shift,
    nmod_t ctx)
{
    int success;
    slong i, j;
    ulong k;
    mp_limb_t alpha;
    n_poly_t Beval;
    n_poly_factor_t Bevalfac;
    slong Blengthx, Blengthy;
    nmod_bpoly_info_t I;
    n_bpoly_t tp, tp1, error;
    n_poly_t ss, tt;
    n_poly_t g;

    n_poly_init(g);
    n_poly_init(Beval);
    n_poly_factor_init(Bevalfac);
    nmod_bpoly_info_init(I, 2);
    n_poly_init(ss);
    n_poly_init(tt);
    n_bpoly_init(tp);
    n_bpoly_init(tp1);
    n_bpoly_init(error);

    /* init done */

    n_bpoly_mod_make_primitive(c, B, ctx);

    Blengthx = B->length;
    FLINT_ASSERT(Blengthx > 1);

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
    if (Beval->length != Blengthx)
        goto next_alpha;

    n_poly_mod_factor(Bevalfac, Beval, ctx);

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
        success = 1;
        goto cleanup;
	}

    n_bpoly_mod_taylor_shift_var1(B, alpha, ctx);

    Blengthy = 0;
    for (i = 0; i < B->length; i++)
        Blengthy = FLINT_MAX(Blengthy, (B->coeffs + i)->length);

    nmod_bpoly_info_clear(I);
    nmod_bpoly_info_init(I, Bevalfac->num);

    n_bpoly_set(I->Btilde, B);
    n_bpoly_mod_make_monic_series(I->Btilde, Blengthy, ctx);

    for (i = 0; i < I->r; i++)
    {
        n_poly_mod_make_monic(I->Bitilde + i, Bevalfac->poly + i, ctx);
        n_bpoly_set_poly_var0(I->newBitilde + i, I->Bitilde + i);
    }

    nmod_bpoly_info_disolve(I, ctx);

    FLINT_ASSERT(I->r >= 2);

    n_bpoly_mod_mullow(tp, I->newBitilde + 0, I->newBitilde + 1, Blengthy, ctx);
    for (i = 2; i < I->r; i++)
    {
        n_bpoly_mod_mullow(tp1, tp, I->newBitilde + i, Blengthy, ctx);
        n_bpoly_swap(tp1, tp);
    }

    n_bpoly_mod_sub(error, I->Btilde, tp, ctx);

    for (j = 1; j < Blengthy; j++)
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

        for (i = 0; i < I->r; i++)
        {
            n_poly_mod_mul(tt, ss, I->d + i, ctx);
            n_poly_mod_rem(tt, tt, I->Bitilde + i, ctx);
            n_bpoly_mod_add_poly_shift(I->newBitilde + i, tt, j, ctx);
        }

        n_bpoly_mod_mullow(tp, I->newBitilde + 0, I->newBitilde + 1, Blengthy, ctx);
        for (i = 2; i < I->r; i++)
        {
            n_bpoly_mod_mullow(tp1, tp, I->newBitilde + i, Blengthy, ctx);
            n_bpoly_swap(tp1, tp);
        }
        n_bpoly_mod_sub(error, I->Btilde, tp, ctx);
    }

    {
        n_bpoly_t f, Q, R, trymez;
        n_bpoly_t tryme, trymet;
        n_poly_t leadf;
        slong kk, *used_arr, *sub_arr;

        used_arr = (slong *) flint_calloc(2 * I->r, sizeof(slong));
        sub_arr  = used_arr + I->r;

        n_bpoly_init(f);
        n_bpoly_init(Q);
        n_bpoly_init(R);
        n_bpoly_init(trymez);
        n_bpoly_init(tryme);
        n_bpoly_init(trymet);
        n_poly_init(leadf);

        n_bpoly_set(f, B);
        FLINT_ASSERT(f->length > 0);
        n_poly_set(leadf, f->coeffs + f->length - 1);

        for (kk = 1; kk < I->r; kk++)
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

                if (sub_arr[kk - 1] > I->r - 1)
                    indx--;
                else
                {
                    for(l = 0; l < kk; l++)
                    {
                        if (used_arr[sub_arr[l]] == 1)
                            break;
                    }

                    n_bpoly_set_poly_var1(tryme, leadf);
                    for (l = 0; l < kk; l++)
                    {
                        n_bpoly_mod_mullow(trymet, tryme, I->newBitilde + sub_arr[l], Blengthy, ctx);
                        n_bpoly_swap(trymet, tryme);
                    }

                    n_bpoly_set(trymez, tryme);
                    n_bpoly_mod_make_primitive(g, trymez, ctx);
                    if (n_bpoly_mod_divides(Q, f, trymez, ctx))
                    {
                        n_bpoly_mod_taylor_shift_var1(trymez, nmod_neg(alpha, ctx), ctx);
                        n_tpoly_fit_length(F, F->length + 1);
                        n_bpoly_swap(F->coeffs + F->length, trymez);
                        F->length++;

                        for(l = 0; l < kk; l++)
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

            for (kk = 0; kk < I->r; kk++)
                test = test + used_arr[kk];

            if (test == 0)
            {
                n_bpoly_mod_taylor_shift_var1(f, nmod_neg(alpha, ctx), ctx);
                n_tpoly_fit_length(F, F->length + 1);
                n_bpoly_swap(F->coeffs + F->length, f);
                F->length++;
            }
        }

        n_bpoly_clear(f);
        n_bpoly_clear(Q);
        n_bpoly_clear(R);
        n_bpoly_clear(trymez);
        n_bpoly_clear(tryme);
        n_bpoly_clear(trymet);
        n_poly_clear(leadf);

        flint_free(used_arr);
    }

    success = 1;

cleanup:

    n_poly_clear(g);
    n_poly_clear(Beval);
    n_poly_factor_clear(Bevalfac);
    nmod_bpoly_info_clear(I);
    n_poly_clear(ss);
    n_poly_clear(tt);
    n_bpoly_clear(tp);
    n_bpoly_clear(tp1);
    n_bpoly_clear(error);

    return success;
}


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


void n_bpoly_mod_factor_lgprime(
    n_poly_t c,
    n_tpoly_t F,
    n_bpoly_t B,
    nmod_t ctx)
{
    slong i, j, k;
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

    Blengthx = B->length;
    FLINT_ASSERT(Blengthx > 1);

    goto got_alpha;

next_alpha:

    deg++;

	fq_nmod_ctx_clear(ectx);
	fq_nmod_ctx_init(ectx, P, deg, "y");

got_alpha:
/*
printf("***alpha: finite field m: "); nmod_poly_print_pretty(ectx->modulus, "y"); printf("\n");
*/
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
    {
        Blengthy = FLINT_MAX(Blengthy, (B->coeffs + i)->length);
    }

    n_poly_one(finalmpow);
    while (finalmpow->length <= Blengthy)
    {
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

    for (j = 1; j < Blengthy; j++)
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

    {
        n_bpoly_t f, Q, R, trymez;
        n_bpoly_t tryme, trymet;
        n_poly_t leadf;
        slong kk, *used_arr, *sub_arr;

        used_arr = (slong *) flint_calloc(2 * I->r, sizeof(slong));
        sub_arr  = used_arr + I->r;

        n_bpoly_init(f);
        n_bpoly_init(Q);
        n_bpoly_init(R);
        n_bpoly_init(trymez);
        n_bpoly_init(tryme);
        n_bpoly_init(trymet);
        n_poly_init(leadf);

        n_bpoly_set(f, B);
        FLINT_ASSERT(f->length > 0);
        n_poly_set(leadf, f->coeffs + f->length - 1);

        for (kk = 1; kk < I->r; kk++)
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

                if (sub_arr[kk - 1] > I->r - 1)
                    indx--;
                else
                {
                    for(l = 0; l < kk; l++)
                    {
                        if (used_arr[sub_arr[l]] == 1)
                            break;
                    }

                    n_bpoly_set_poly_var1(tryme, leadf);
                    for (l = 0; l < kk; l++)
                    {
                        n_bpoly_mod_mulmod_poly(trymet, tryme, I->newBitilde + sub_arr[l], finalmpow, ctx);
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

            for (kk = 0; kk < I->r; kk++)
                test = test + used_arr[kk];

            if (test == 0)
            {
                n_tpoly_fit_length(F, F->length + 1);
                n_bpoly_swap(F->coeffs + F->length, f);
                F->length++;
            }
        }

        n_bpoly_clear(f);
        n_bpoly_clear(Q);
        n_bpoly_clear(R);
        n_bpoly_clear(trymez);
        n_bpoly_clear(tryme);
        n_bpoly_clear(trymet);
        n_poly_clear(leadf);

        flint_free(used_arr);
    }

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

    return;
}
