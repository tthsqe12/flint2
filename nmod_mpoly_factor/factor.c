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





/* a *= b^e */
void nmod_mpoly_factor_mul_factor_fmpz(
    nmod_mpoly_factor_t a,
    const nmod_mpoly_factor_t b,
    const fmpz_t e,
    const nmod_mpoly_ctx_t ctx)
{
    slong i;
    fmpz_t t;

    fmpz_init(t);

    a->constant = nmod_mul(a->constant,
            nmod_pow_fmpz(b->constant, e, ctx->ffinfo->mod), ctx->ffinfo->mod);

    for (i = 0; i < b->num; i++)
    {
        fmpz_mul(t, b->exp + i, e);
        nmod_mpoly_factor_append_fmpz(a, b->poly + i, t, ctx);
    }

    fmpz_clear(t);
}

void nmod_mpoly_factor_mul_factor_ui(
    nmod_mpoly_factor_t a,
    const nmod_mpoly_factor_t b,
    ulong e,
    const nmod_mpoly_ctx_t ctx)
{
    slong i;
    fmpz_t t;

    fmpz_init(t);

    a->constant = nmod_mul(a->constant,
              nmod_pow_ui(b->constant, e, ctx->ffinfo->mod), ctx->ffinfo->mod);

    for (i = 0; i < b->num; i++)
    {
        fmpz_mul_ui(t, b->exp + i, e);
        nmod_mpoly_factor_append_fmpz(a, b->poly + i, t, ctx);
    }

    fmpz_clear(t);
}



/* f *= g^e, no combination */
void nmod_mpoly_factor_append_factor_fmpz(
	nmod_mpoly_factor_t f,
	nmod_mpoly_factor_t g,
	fmpz_t e,
	const nmod_mpoly_ctx_t ctx)
{
	slong i;
	fmpz_t c;

	fmpz_init(c);

	f->constant = nmod_mul(f->constant, g->constant, ctx->ffinfo->mod);

	for (i = 0; i < g->num; i++)
	{
		fmpz_mul(c, g->exp + i, e);
		nmod_mpoly_factor_append_fmpz(f, g->poly + i, c, ctx);
	}

	fmpz_clear(c);
}

void nmod_mpoly_factor_append_factor_ui(
	nmod_mpoly_factor_t f,
	nmod_mpoly_factor_t g,
	ulong e,
	const nmod_mpoly_ctx_t ctx)
{
	slong i;
	fmpz_t c;

	fmpz_init(c);

	f->constant = nmod_mul(f->constant,
              nmod_pow_ui(g->constant, e, ctx->ffinfo->mod), ctx->ffinfo->mod);

	for (i = 0; i < g->num; i++)
	{
		fmpz_mul_si(c, g->exp + i, e);
		nmod_mpoly_factor_append_fmpz(f, g->poly + i, c, ctx);
	}

	fmpz_clear(c);
}


void n_bpoly_clear(n_bpoly_t A)
{
    slong i;
    if (A->alloc > 0)
    {
        FLINT_ASSERT(A->coeffs != NULL);
        for (i = 0; i < A->alloc; i++)
            n_poly_clear(A->coeffs + i);
        flint_free(A->coeffs);
    }
    else
    {
        FLINT_ASSERT(A->coeffs == NULL);
    }
}

void fq_nmod_bpoly_clear(fq_nmod_bpoly_t A, const fq_nmod_ctx_t ctx)
{
    slong i;
    if (A->alloc > 0)
    {
        FLINT_ASSERT(A->coeffs != NULL);
        for (i = 0; i < A->alloc; i++)
            fq_nmod_poly_clear(A->coeffs + i, ctx);
        flint_free(A->coeffs);
    }
    else
    {
        FLINT_ASSERT(A->coeffs == NULL);
    }
}

void n_bpoly_print_pretty(
    const n_bpoly_t A,
    const char * xvar,
    const char * yvar)
{
    slong i;
    int first;

    first = 1;
    for (i = A->length - 1; i >= 0; i--)
    {
        if (n_poly_is_zero(A->coeffs + i))
            continue;

        if (!first)
            flint_printf(" + ");

        first = 0;

        flint_printf("(");
        n_poly_print_pretty(A->coeffs + i, yvar);
        flint_printf(")*%s^%wd", xvar, i);
    }

    if (first)
        flint_printf("0");
}

void fq_nmod_bpoly_print_pretty(
    const fq_nmod_bpoly_t A,
    const char * xvar,
    const char * yvar,
    const fq_nmod_ctx_t ectx)
{
    slong i;
    int first;

    first = 1;
    for (i = A->length - 1; i >= 0; i--)
    {
        if (fq_nmod_poly_is_zero(A->coeffs + i, ectx))
        {
            FLINT_ASSERT(i < A->length - 1);
            continue;
        }

        if (!first)
            flint_printf(" + ");

        first = 0;

        flint_printf("(");
        fq_nmod_poly_print_pretty(A->coeffs + i, yvar, ectx);
        flint_printf(")*%s^%wd", xvar, i);
    }

    if (first)
        flint_printf("0");
}

void n_bpoly_realloc(n_bpoly_t A, slong len)
{
    slong i;
    slong old_alloc = A->alloc;
    slong new_alloc = FLINT_MAX(len, old_alloc + 1 + old_alloc/2);

    FLINT_ASSERT(A->alloc >= 0);

    if (len <= A->alloc)
        return;

    if (A->alloc > 0)
    {
        FLINT_ASSERT(A->coeffs != NULL);
        A->coeffs = (n_poly_struct *) flint_realloc(A->coeffs,
                                       new_alloc * sizeof(n_poly_struct));
    }
    else
    {
        FLINT_ASSERT(A->coeffs == NULL);
        A->coeffs = (n_poly_struct *) flint_malloc(
                                       new_alloc * sizeof(n_poly_struct));
    }

    for (i = old_alloc; i < new_alloc; i++)
        n_poly_init(A->coeffs + i);

    A->alloc = len;
}

void fq_nmod_bpoly_realloc(fq_nmod_bpoly_t A, slong len, const fq_nmod_ctx_t ctx)
{
    slong i;
    slong old_alloc = A->alloc;
    slong new_alloc = FLINT_MAX(len, old_alloc + 1 + old_alloc/2);

    FLINT_ASSERT(A->alloc >= 0);

    if (len <= A->alloc)
        return;

    if (A->alloc > 0)
    {
        FLINT_ASSERT(A->coeffs != NULL);
        A->coeffs = (fq_nmod_poly_struct *) flint_realloc(A->coeffs,
                                      new_alloc * sizeof(fq_nmod_poly_struct));
    }
    else
    {
        FLINT_ASSERT(A->coeffs == NULL);
        A->coeffs = (fq_nmod_poly_struct *) flint_malloc(
                                      new_alloc * sizeof(fq_nmod_poly_struct));
    }

    for (i = old_alloc; i < new_alloc; i++)
        fq_nmod_poly_init(A->coeffs + i, ctx);

    A->alloc = len;
}

void fq_nmod_bpoly_set_polyx(
    fq_nmod_bpoly_t A,
    const fq_nmod_poly_t B,
    const fq_nmod_ctx_t ectx)
{
    slong i;
    fq_nmod_bpoly_fit_length(A, B->length, ectx);
    A->length = 0;
    for (i = 0; i < B->length; i++)
    {
        fq_nmod_poly_set_fq_nmod(A->coeffs + i, B->coeffs + i, ectx);
        if (!fq_nmod_poly_is_zero(A->coeffs + i, ectx))
            A->length = i + 1;
    }    
}

void fq_nmod_bpoly_set_polyy(
    fq_nmod_bpoly_t A,
    const fq_nmod_poly_t B,
    const fq_nmod_ctx_t ectx)
{
    fq_nmod_bpoly_fit_length(A, 1, ectx);
	fq_nmod_poly_set(A->coeffs + 0, B, ectx);
	A->length = !fq_nmod_poly_is_zero(A->coeffs + 0, ectx);
}

void fq_nmod_bpoly_get_coeff(fq_nmod_t c, const fq_nmod_bpoly_t A, slong xi, slong yi, const fq_nmod_ctx_t ectx)
{
    if (xi >= A->length)
        fq_nmod_zero(c, ectx);
    else
        fq_nmod_poly_get_coeff(c, A->coeffs + xi, yi, ectx);
}


void fq_nmod_bpoly_make_monic(fq_nmod_bpoly_t A, slong order, const fq_nmod_ctx_t ctx)
{
    slong i;
    fq_nmod_poly_t t, lcinv;

    FLINT_ASSERT(A->length > 0);

    fq_nmod_poly_init(t, ctx);
    fq_nmod_poly_init(lcinv, ctx);
    fq_nmod_poly_inv_series(lcinv, A->coeffs + A->length - 1, order, ctx);

    for (i = 0; i < A->length; i++)
    {
        fq_nmod_poly_mullow(t, A->coeffs + i, lcinv, order, ctx);
        fq_nmod_poly_swap(A->coeffs + i, t, ctx);
    }

    fq_nmod_poly_clear(t, ctx);
    fq_nmod_poly_clear(lcinv, ctx);
}

void fq_nmod_bpoly_mullow(
    fq_nmod_bpoly_t A,
    const fq_nmod_bpoly_t B,
    const fq_nmod_bpoly_t C,
    slong order,
    const fq_nmod_ctx_t ectx)
{
    slong i, j;
    fq_nmod_poly_t t;

    FLINT_ASSERT(A != B);
    FLINT_ASSERT(A != C);

    fq_nmod_poly_init(t, ectx);

    fq_nmod_bpoly_fit_length(A, B->length + C->length - 1, ectx);
    for (i = 0; i < B->length + C->length - 1; i++)
        fq_nmod_poly_zero(A->coeffs + i, ectx);

    for (i = 0; i < B->length; i++)
    {
        for (j = 0; j < C->length; j++)
        {
            fq_nmod_poly_mullow(t, B->coeffs + i, C->coeffs + j, order, ectx);
            fq_nmod_poly_add(A->coeffs + i + j, A->coeffs + i + j, t, ectx);
        }
    }

    A->length = B->length + C->length - 1;
    while (A->length > 0 && fq_nmod_poly_is_zero(A->coeffs + A->length - 1, ectx))
        A->length--;

    fq_nmod_poly_clear(t, ectx);
}


void fq_nmod_bpoly_add_poly_shift(
    fq_nmod_bpoly_t A,
    const fq_nmod_poly_t B,
    slong yshift,
    const fq_nmod_ctx_t ctx)
{
    slong i;

    FLINT_ASSERT(A->length > B->length);

    for (i = 0; i < B->length; i++)
    {
        fq_nmod_t ccc;
        fq_nmod_init(ccc, ctx);
        fq_nmod_poly_get_coeff(ccc, A->coeffs + i, yshift, ctx);
        FLINT_ASSERT(fq_nmod_is_zero(ccc, ctx));
        fq_nmod_poly_set_coeff(A->coeffs + i, yshift, B->coeffs + i, ctx);
        fq_nmod_clear(ccc, ctx);
    }
}


void fq_nmod_bpoly_sub(
    fq_nmod_bpoly_t A,
    const fq_nmod_bpoly_t B,
    const fq_nmod_bpoly_t C,
    const fq_nmod_ctx_t ectx)
{
    slong i;

    FLINT_ASSERT(A != B);
    FLINT_ASSERT(A != C);

    fq_nmod_bpoly_fit_length(A, FLINT_MAX(B->length, C->length), ectx);
    A->length = 0;
    for (i = 0; i < FLINT_MAX(B->length, C->length) - 1; i++)
    {
        if (i < B->length)
        {
            if (i < C->length)
            {
                fq_nmod_poly_sub(A->coeffs + i, B->coeffs + i, C->coeffs + i, ectx);
            }
            else
            {
                fq_nmod_poly_set(A->coeffs + i, B->coeffs + i, ectx);
            }
        }
        else
        {
            FLINT_ASSERT(i < C->length);
            fq_nmod_poly_neg(A->coeffs + i, C->coeffs + i, ectx);
        }

        if (!fq_nmod_poly_is_zero(A->coeffs + i, ectx))
            A->length = i + 1;
    }
}

void n_bpoly_set(n_bpoly_t A, const n_bpoly_t B)
{
    slong i;

    FLINT_ASSERT(A != B);

    n_bpoly_fit_length(A, B->length);
    A->length = B->length;

    for (i = 0; i < B->length; i++)
        n_poly_set(A->coeffs + i, B->coeffs + i);
}

void fq_nmod_bpoly_set(fq_nmod_bpoly_t A, const fq_nmod_bpoly_t B, const fq_nmod_ctx_t ctx)
{
    slong i;

    FLINT_ASSERT(A != B);

    fq_nmod_bpoly_fit_length(A, B->length, ctx);
    A->length = B->length;

    for (i = 0; i < B->length; i++)
        fq_nmod_poly_set(A->coeffs + i, B->coeffs + i, ctx);
}

int fq_nmod_bpoly_divides(
    fq_nmod_bpoly_t Q,
    const fq_nmod_bpoly_t A,
    const fq_nmod_bpoly_t B,
    const fq_nmod_ctx_t ctx)
{
    slong i, qoff;
    int divides;
    fq_nmod_poly_t q, t;
    fq_nmod_bpoly_t R;

    FLINT_ASSERT(Q != A);
    FLINT_ASSERT(Q != B);
    FLINT_ASSERT(B->length > 0);

    fq_nmod_poly_init(q, ctx);
    fq_nmod_poly_init(t, ctx);
    fq_nmod_bpoly_init(R, ctx);
    fq_nmod_bpoly_set(R, A, ctx);

    Q->length = 0;

    while (R->length >= B->length)
    {
        fq_nmod_poly_divrem(q, t, R->coeffs + R->length - 1, B->coeffs + B->length - 1, ctx);
        if (!fq_nmod_poly_is_zero(t, ctx))
        {
            divides = 0;
            goto cleanup;
        }

        for (i = 0; i < B->length; i++)
        {
            fq_nmod_poly_mul(t, B->coeffs + i, q, ctx);
            fq_nmod_poly_sub(R->coeffs + i + R->length - B->length, R->coeffs + i + R->length - B->length, t, ctx);
        }

        qoff = R->length - B->length;

        FLINT_ASSERT(qoff >= 0);
        if (qoff >= Q->length)
        {
            fq_nmod_bpoly_fit_length(Q, qoff + 1, ctx);
            for (i = Q->length; i <= qoff; i++)
                fq_nmod_poly_zero(Q->coeffs + i, ctx);
            Q->length = qoff + 1;
        }

        fq_nmod_poly_set(Q->coeffs + qoff, q, ctx);

        while (R->length > 0 && fq_nmod_poly_is_zero(R->coeffs + R->length - 1, ctx))
            R->length--;
    }

    divides = (R->length == 0);

cleanup:

    fq_nmod_poly_clear(q, ctx);
    fq_nmod_poly_clear(t, ctx);
    fq_nmod_bpoly_clear(R, ctx);

    return divides;
}

void fq_nmod_bpoly_make_primitive(fq_nmod_bpoly_t A, const fq_nmod_ctx_t ctx)
{
    slong i;
    fq_nmod_poly_t g, q, r;

    fq_nmod_poly_init(g, ctx);
    fq_nmod_poly_init(q, ctx);
    fq_nmod_poly_init(r, ctx);

    for (i = 0; i < A->length; i++)
	{
        fq_nmod_poly_gcd(q, g, A->coeffs + i, ctx);
		fq_nmod_poly_swap(g, q, ctx);
	}

    for (i = 0; i < A->length; i++)
    {
        fq_nmod_poly_divrem(q, r, A->coeffs + i, g, ctx);
		fq_nmod_poly_swap(A->coeffs + i, q, ctx);
    }

    fq_nmod_poly_clear(g, ctx);
    fq_nmod_poly_clear(q, ctx);
    fq_nmod_poly_clear(r, ctx);
}


void n_bpoly_set_coeff_nonzero(n_bpoly_t A, slong xi, slong yi, mp_limb_t c)
{
    slong i;

    FLINT_ASSERT(c != 0);

    if (xi >= A->length)
    {
        n_bpoly_fit_length(A, xi + 1);
        for (i = A->length; i <= xi; i++)
            n_poly_zero(A->coeffs + i);
        A->length = xi + 1;
    }

    n_poly_set_coeff_nonzero(A->coeffs + xi, yi, c);
    FLINT_ASSERT(!n_poly_is_zero(A->coeffs + A->length - 1));
}

void n_bpoly_set_coeff(n_bpoly_t A, slong xi, slong yi, mp_limb_t c)
{
    slong i;

    if (xi >= A->length)
    {
        n_bpoly_fit_length(A, xi + 1);
        for (i = A->length; i <= xi; i++)
            n_poly_zero(A->coeffs + i);
        A->length = xi + 1;
    }

    n_poly_set_coeff(A->coeffs + xi, yi, c);
    while (A->length > 0 && n_poly_is_zero(A->coeffs + A->length - 1))
        A->length--;
}

void fq_nmod_bpoly_set_coeff(
    fq_nmod_bpoly_t A,
    slong xi, slong yi,
    const fq_nmod_t c,
    const fq_nmod_ctx_t ctx)
{
    slong i;

    if (xi >= A->length)
    {
        fq_nmod_bpoly_fit_length(A, xi + 1, ctx);
        for (i = A->length; i <= xi; i++)
            fq_nmod_poly_zero(A->coeffs + i, ctx);
        A->length = xi + 1;
    }

    fq_nmod_poly_set_coeff(A->coeffs + xi, yi, c, ctx);
    while (A->length > 0 && fq_nmod_poly_is_zero(A->coeffs + A->length - 1, ctx))
        A->length--;
}


void nmod_bpoly_zero(n_bpoly_t A)
{
    A->length = 0;
}

void fq_nmod_bpoly_zero(fq_nmod_bpoly_t A, const fq_nmod_ctx_t fqctx)
{
    A->length = 0;
}

void fq_nmod_mpoly_to_bpoly(
    fq_nmod_bpoly_t A,
    const fq_nmod_mpoly_t B,
    slong varx,
    slong vary,
    const fq_nmod_mpoly_ctx_t ctx)
{
    slong j;
    slong NB;
    ulong Bexpx, Bexpy;
    slong Boffx, Bshiftx, Boffy, Bshifty;
    ulong mask;

    FLINT_ASSERT(B->bits <= FLINT_BITS);
    NB = mpoly_words_per_exp_sp(B->bits, ctx->minfo);

    mpoly_gen_offset_shift_sp(&Boffx, &Bshiftx, varx, B->bits, ctx->minfo);
    mpoly_gen_offset_shift_sp(&Boffy, &Bshifty, vary, B->bits, ctx->minfo);
    mask = (-UWORD(1)) >> (FLINT_BITS - B->bits);

    fq_nmod_bpoly_zero(A, ctx->fqctx);

    for (j = 0; j < B->length; j++)
    {
        Bexpx = ((B->exps + NB*j)[Boffx] >> Bshiftx) & mask;
        Bexpy = ((B->exps + NB*j)[Boffy] >> Bshifty) & mask;
        fq_nmod_bpoly_set_coeff(A, Bexpx, Bexpy, B->coeffs + j, ctx->fqctx);
    }
}


void fq_nmod_mpoly_from_fq_nmod_bpoly(
    fq_nmod_mpoly_t A,
    flint_bitcnt_t Abits,
    const fq_nmod_bpoly_t B,
    slong varx,
    slong vary,
    const fq_nmod_mpoly_ctx_t ctx)
{
    slong n = ctx->minfo->nvars;
    slong i, j;
    slong NA;
    slong Alen;
    fq_nmod_struct * Acoeff;
    ulong * Aexp;
    slong Aalloc;
    ulong * Aexps;
    TMP_INIT;

    FLINT_ASSERT(B->length > 0);
    FLINT_ASSERT(Abits <= FLINT_BITS);
    TMP_START;

    Aexps = (ulong *) TMP_ALLOC(n*sizeof(ulong));
    for (i = 0; i < n; i++)
        Aexps[i] = 0;

    NA = mpoly_words_per_exp(Abits, ctx->minfo);
    fq_nmod_mpoly_fit_bits(A, Abits, ctx);
    A->bits = Abits;

    Acoeff = A->coeffs;
    Aexp = A->exps;
    Aalloc = A->alloc;
    Alen = 0;
    for (i = 0; i < B->length; i++)
    {
        fq_nmod_poly_struct * Bc = B->coeffs + i;
        _fq_nmod_mpoly_fit_length(&Acoeff, &Aexp, &Aalloc, Alen + Bc->length, NA, ctx->fqctx);

        for (j = 0; j < Bc->length; j++)
        {
            if (fq_nmod_is_zero(Bc->coeffs + j, ctx->fqctx))
                continue;

            Aexps[varx] = i;
            Aexps[vary] = j;
            fq_nmod_set(Acoeff + Alen, Bc->coeffs + j, ctx->fqctx);
            mpoly_set_monomial_ui(Aexp + NA*Alen, Aexps, Abits, ctx->minfo);
            Alen++;
        }
    }
    A->coeffs = Acoeff;
    A->exps = Aexp;
    A->alloc = Aalloc;
    A->length = Alen;

    TMP_END;

    fq_nmod_mpoly_sort_terms(A, ctx);
}


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

void _fq_nmod_poly_taylor_shift_horner(
    fq_nmod_struct * poly,
    const fq_nmod_t c,
    slong n,
    const fq_nmod_ctx_t ctx)
{
    slong i, j;
    fq_nmod_t p;
    fq_nmod_init(p, ctx);
    for (i = n - 2; i >= 0; i--)
    {
        for (j = i; j < n - 1; j++)
        {
            fq_nmod_mul(p, poly + j + 1, c, ctx);
            fq_nmod_add(poly + j, poly + j, p, ctx);
        }
    }
    fq_nmod_clear(p, ctx);
}

void fq_nmod_poly_taylor_shift_horner(
    fq_nmod_poly_t g,
    const fq_nmod_poly_t f,
    const fq_nmod_t c,
    const fq_nmod_ctx_t ctx)
{
    if (f != g)
        fq_nmod_poly_set(g, f, ctx);

    _fq_nmod_poly_taylor_shift_horner(g->coeffs, c, g->length, ctx);
}


void fq_nmod_bpoly_taylor_shift(
    fq_nmod_bpoly_t A,
    const fq_nmod_t alpha,
    const fq_nmod_ctx_t ctx)
{
    slong i;
    for (i = A->length - 1; i >= 0; i--)
        fq_nmod_poly_taylor_shift_horner(A->coeffs + i, A->coeffs + i, alpha, ctx);    
}




/*****************************************************************************/


typedef struct {
    slong r; /* number of local factors */
    fq_nmod_bpoly_t Btilde;
    fq_nmod_bpoly_struct * newBitilde;
    fq_nmod_poly_struct * P;
    fq_nmod_poly_struct * d;
    fq_nmod_poly_struct * Bitilde;
} fq_nmod_smprime_info_struct;

typedef fq_nmod_smprime_info_struct fq_nmod_smprime_info_t[1];


void fq_nmod_smprime_info_init(
    fq_nmod_smprime_info_t I,
    slong r,
    const fq_nmod_ctx_t fqctx)
{
    slong i;

    FLINT_ASSERT(r >= 2);

    I->r = r;

    fq_nmod_bpoly_init(I->Btilde, fqctx);

    I->newBitilde = (fq_nmod_bpoly_struct *) flint_malloc(I->r*sizeof(fq_nmod_bpoly_struct));
    I->d          = (fq_nmod_poly_struct *) flint_malloc(I->r*sizeof(fq_nmod_poly_struct));
    I->Bitilde    = (fq_nmod_poly_struct *) flint_malloc(I->r*sizeof(fq_nmod_poly_struct));

    for (i = 0; i < I->r; i++)
    {
        fq_nmod_bpoly_init(I->newBitilde + i, fqctx);
        fq_nmod_poly_init(I->d + i, fqctx);
        fq_nmod_poly_init(I->Bitilde + i, fqctx);
    }
}

void fq_nmod_smprime_info_clear(fq_nmod_smprime_info_t I,
    const fq_nmod_ctx_t fqctx)
{
    slong i;

    fq_nmod_bpoly_clear(I->Btilde, fqctx);

    for (i = 0; i < I->r; i++)
    {
        fq_nmod_bpoly_clear(I->newBitilde + i, fqctx);
        fq_nmod_poly_clear(I->d + i, fqctx);
        fq_nmod_poly_clear(I->Bitilde + i, fqctx);
    }

    flint_free(I->newBitilde);
    flint_free(I->d);
    flint_free(I->Bitilde);
}




/*
    set out[i] so that
    1/(f[0]*f[1]*...*f[n-1]) = out[0]/f[0] + ... + out[n-1]/f[n-1]
*/
int fq_nmod_partial_fraction_coeffs(
    slong r,
    fq_nmod_poly_struct * out,
    const fq_nmod_poly_struct * f,
    const fq_nmod_ctx_t ectx)
{
    slong i;
    fq_nmod_poly_t num, den, a, b, g, t;

    FLINT_ASSERT(r >= 2);

    fq_nmod_poly_init(num, ectx);
    fq_nmod_poly_init(den, ectx);
    fq_nmod_poly_init(a, ectx);
    fq_nmod_poly_init(b, ectx);
    fq_nmod_poly_init(g, ectx);
    fq_nmod_poly_init(t, ectx);

    fq_nmod_poly_one(num, ectx);
    fq_nmod_poly_mul(den, f + 0, f + 1, ectx);
    for (i = 2; i < r; i++)
        fq_nmod_poly_mul(den, den, f + i, ectx);

    for (i = 0; i < r; i++)
    {
        fq_nmod_poly_divrem(den, t, den, f + i, ectx);
        FLINT_ASSERT(fq_nmod_poly_is_zero(t, ectx));
        fq_nmod_poly_xgcd(g, a, b, f + i, den, ectx);
        if (fq_nmod_poly_degree(g, ectx) != 0)
            return 0;
        FLINT_ASSERT(fq_nmod_poly_is_one(g, ectx));
        fq_nmod_poly_mul(t, b, num, ectx);
        fq_nmod_poly_rem(out + i, t, f + i, ectx);
        fq_nmod_poly_mul(t, a, num, ectx);
        fq_nmod_poly_rem(num, t, den, ectx);
    }

    fq_nmod_poly_clear(num, ectx);
    fq_nmod_poly_clear(den, ectx);
    fq_nmod_poly_clear(a, ectx);
    fq_nmod_poly_clear(b, ectx);
    fq_nmod_poly_clear(g, ectx);
    fq_nmod_poly_clear(t, ectx);
    return 1;
}

int fq_nmod_smprime_info_disolve(fq_nmod_smprime_info_t I, const fq_nmod_ctx_t ctx)
{
    return fq_nmod_partial_fraction_coeffs(I->r, I->d, I->Bitilde, ctx);
}




static int _fq_nmod_irreducible_bivar_factors_smprime(
    fq_nmod_mpoly_factor_t fac,
    const fq_nmod_mpoly_t A,
    slong xvar,
    slong yvar,
    const fq_nmod_mpoly_ctx_t ctx)
{
    int success;
    slong i, j;
    ulong k;
    fq_nmod_t alpha;
    fq_nmod_poly_t Beval;
    fq_nmod_bpoly_t B;
    fq_nmod_poly_factor_t Bevalfac;
    fq_nmod_t Bevalfaclc;
    slong Blengthx, Blengthy;
    fq_nmod_smprime_info_t I;
    fq_nmod_bpoly_t tp, tp1, error;
    fq_nmod_poly_t ss, tt;

    fq_nmod_init(alpha, ctx->fqctx);
    fq_nmod_init(Bevalfaclc, ctx->fqctx);
    fq_nmod_poly_init(ss, ctx->fqctx);
    fq_nmod_poly_init(tt, ctx->fqctx);
    fq_nmod_bpoly_init(tp, ctx->fqctx);
    fq_nmod_bpoly_init(tp1, ctx->fqctx);
    fq_nmod_bpoly_init(error, ctx->fqctx);

    fq_nmod_poly_init(Beval, ctx->fqctx);
    fq_nmod_poly_factor_init(Bevalfac, ctx->fqctx);
    fq_nmod_bpoly_init(B, ctx->fqctx);
    fq_nmod_smprime_info_init(I, 2, ctx->fqctx);

    /* init done */

    fq_nmod_mpoly_factor_one(fac, ctx);

    fq_nmod_mpoly_to_bpoly(B, A, xvar, yvar, ctx);


    Blengthx = B->length;
    FLINT_ASSERT(Blengthx > 1);

    fq_nmod_set_ui(alpha, 0, ctx->fqctx);
    goto got_alpha;

next_alpha:

    if (fq_nmod_next(alpha, ctx->fqctx) == 0)
    {
        /* TODO leak */
        success = 0;
        goto cleanup;
    }

got_alpha:

    fq_nmod_bpoly_eval(Beval, B, alpha, ctx->fqctx);

    /* if killed leading coeff, get new alpha */
    if (Beval->length != Blengthx)
        goto next_alpha;

    Bevalfac->num = 0;
    fq_nmod_poly_factor(Bevalfac, Bevalfaclc, Beval, ctx->fqctx);

    /* if multiple factors, get new alpha */
    for (i = 0; i < Bevalfac->num; i++)
    {
        if (Bevalfac->exp[i] != 1)
            goto next_alpha;
    }

	/* if one factor, A is irreducible */
	if (Bevalfac->num == 1)
	{
		fq_nmod_mpoly_factor_append_ui(fac, A, 1, ctx);
        success = 1;
        goto cleanup;
	}

    fq_nmod_bpoly_taylor_shift(B, alpha, ctx->fqctx);

    Blengthy = 0;
    for (i = 0; i < B->length; i++)
    {
        Blengthy = FLINT_MAX(Blengthy, (B->coeffs + i)->length);
    }
    fq_nmod_smprime_info_clear(I, ctx->fqctx);
    fq_nmod_smprime_info_init(I, Bevalfac->num, ctx->fqctx);

    fq_nmod_bpoly_set(I->Btilde, B, ctx->fqctx);
    fq_nmod_bpoly_make_monic(I->Btilde, Blengthy, ctx->fqctx);

    for (i = 0; i < I->r; i++)
    {
        fq_nmod_poly_make_monic(I->Bitilde + i, Bevalfac->poly + i, ctx->fqctx);
        fq_nmod_bpoly_set_polyx(I->newBitilde + i, I->Bitilde + i, ctx->fqctx);
    }

    fq_nmod_smprime_info_disolve(I, ctx->fqctx);

    FLINT_ASSERT(I->r > 1);

    fq_nmod_bpoly_mullow(tp, I->newBitilde + 0, I->newBitilde + 1, Blengthy, ctx->fqctx);
    for (i = 2; i < I->r; i++)
    {
        fq_nmod_bpoly_mullow(tp1, tp, I->newBitilde + i, Blengthy, ctx->fqctx);
        fq_nmod_bpoly_swap(tp1, tp);
    }

    fq_nmod_bpoly_sub(error, I->Btilde, tp, ctx->fqctx);

    for (j = 1; j < Blengthy; j++)
    {
        fq_nmod_poly_zero(ss, ctx->fqctx);
        for (i = error->length - 1; i >= 0; i--)
        {
            fq_nmod_t ccc;
            fq_nmod_init(ccc, ctx->fqctx);
            fq_nmod_bpoly_get_coeff(ccc, error, i, j, ctx->fqctx);
            fq_nmod_poly_set_coeff(ss, i, ccc, ctx->fqctx);
            for (k = 0; k < j; k++)
            {
                fq_nmod_bpoly_get_coeff(ccc, error, i, k, ctx->fqctx);
                FLINT_ASSERT(fq_nmod_is_zero(ccc, ctx->fqctx));
            }
            fq_nmod_clear(ccc, ctx->fqctx);
        }

        for (i = 0; i < I->r; i++)
        {
            fq_nmod_poly_mul(tt, ss, I->d + i, ctx->fqctx);
            fq_nmod_poly_rem(tt, tt, I->Bitilde + i, ctx->fqctx);
            fq_nmod_bpoly_add_poly_shift(I->newBitilde + i, tt, j, ctx->fqctx);
        }

        fq_nmod_bpoly_mullow(tp, I->newBitilde + 0, I->newBitilde + 1, Blengthy, ctx->fqctx);
        for (i = 2; i < I->r; i++)
        {
            fq_nmod_bpoly_mullow(tp1, tp, I->newBitilde + i, Blengthy, ctx->fqctx);
            fq_nmod_bpoly_swap(tp1, tp);
        }
        fq_nmod_bpoly_sub(error, I->Btilde, tp, ctx->fqctx);
    }

    {
        fq_nmod_bpoly_t f, Q, R, trymez;
        fq_nmod_bpoly_t tryme, trymet;
        fq_nmod_poly_t leadf;
        slong kk, *used_arr, *sub_arr;

        used_arr = (slong *) flint_calloc(2 * I->r, sizeof(slong));
        sub_arr  = used_arr + I->r;

        fq_nmod_bpoly_init(f, ctx->fqctx);
        fq_nmod_bpoly_init(Q, ctx->fqctx);
        fq_nmod_bpoly_init(R, ctx->fqctx);
        fq_nmod_bpoly_init(trymez, ctx->fqctx);
        fq_nmod_bpoly_init(tryme, ctx->fqctx);
        fq_nmod_bpoly_init(trymet, ctx->fqctx);
        fq_nmod_poly_init(leadf, ctx->fqctx);

        fq_nmod_bpoly_set(f, B, ctx->fqctx);
        FLINT_ASSERT(f->length > 0);
        fq_nmod_poly_set(leadf, f->coeffs + f->length - 1, ctx->fqctx);

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

                    fq_nmod_bpoly_set_polyy(tryme, leadf, ctx->fqctx);
                    for (l = 0; l < kk; l++)
                    {
                        fq_nmod_bpoly_mullow(trymet, tryme, I->newBitilde + sub_arr[l], Blengthy, ctx->fqctx);
                        fq_nmod_bpoly_swap(trymet, tryme);
                    }

                    fq_nmod_bpoly_set(trymez, tryme, ctx->fqctx);
                    fq_nmod_bpoly_make_primitive(trymez, ctx->fqctx);
                    if (fq_nmod_bpoly_divides(Q, f, trymez, ctx->fqctx))
                    {
                        fq_nmod_mpoly_t goodtry;
                        fq_nmod_neg(alpha, alpha, ctx->fqctx);
                        fq_nmod_bpoly_taylor_shift(trymez, alpha, ctx->fqctx);
                        fq_nmod_neg(alpha, alpha, ctx->fqctx);
                        fq_nmod_mpoly_init(goodtry, ctx);
                        fq_nmod_mpoly_from_fq_nmod_bpoly(goodtry, A->bits, trymez, xvar, yvar, ctx);
                        fq_nmod_mpoly_make_monic(goodtry, goodtry, ctx);
                        fq_nmod_mpoly_factor_append_ui(fac, goodtry, 1, ctx);
                        fq_nmod_mpoly_clear(goodtry, ctx);

                        for(l = 0; l < kk; l++)
                        {
                            used_arr[sub_arr[l]] = 1;
                            count++;
                        }

                        fq_nmod_bpoly_set(f, Q, ctx->fqctx);
                        FLINT_ASSERT(f->length > 0);
                        fq_nmod_poly_set(leadf, f->coeffs + f->length - 1, ctx->fqctx);

                     /* If r - count = k then the rest are irreducible.  
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
                fq_nmod_mpoly_t goodtry;
                fq_nmod_neg(alpha, alpha, ctx->fqctx);
                fq_nmod_bpoly_taylor_shift(f, alpha, ctx->fqctx);
                fq_nmod_neg(alpha, alpha, ctx->fqctx);
                fq_nmod_mpoly_init(goodtry, ctx);
                fq_nmod_mpoly_from_fq_nmod_bpoly(goodtry, A->bits, f, xvar, yvar, ctx);
                fq_nmod_mpoly_make_monic(goodtry, goodtry, ctx);
                fq_nmod_mpoly_factor_append_ui(fac, goodtry, 1, ctx);
                fq_nmod_mpoly_clear(goodtry, ctx);
            }
        }

        flint_free(used_arr);
        fq_nmod_bpoly_clear(f, ctx->fqctx);
        fq_nmod_bpoly_clear(Q, ctx->fqctx);
        fq_nmod_bpoly_clear(R, ctx->fqctx);
        fq_nmod_bpoly_clear(trymez, ctx->fqctx);
        fq_nmod_bpoly_clear(tryme, ctx->fqctx);
        fq_nmod_bpoly_clear(trymet, ctx->fqctx);
        fq_nmod_poly_clear(leadf, ctx->fqctx);
    }

    success = 1;

cleanup:

    fq_nmod_clear(alpha, ctx->fqctx);
    fq_nmod_clear(Bevalfaclc, ctx->fqctx);
    fq_nmod_poly_clear(ss, ctx->fqctx);
    fq_nmod_poly_clear(tt, ctx->fqctx);
    fq_nmod_bpoly_clear(tp, ctx->fqctx);
    fq_nmod_bpoly_clear(tp1, ctx->fqctx);
    fq_nmod_bpoly_clear(error, ctx->fqctx);

    fq_nmod_poly_clear(Beval, ctx->fqctx);
    fq_nmod_poly_factor_clear(Bevalfac, ctx->fqctx);
    fq_nmod_bpoly_clear(B, ctx->fqctx);
    fq_nmod_smprime_info_clear(I, ctx->fqctx);

    return success;
}




/*
void fq_nmod_poly_shift_nmod_poly(
    fq_nmod_poly_t A,
    nmod_poly_t B,
    const fq_nmod_ctx_t ectx)
{
    slong i;
    nmod_poly_t q, b, r;
    nmod_poly_init_mod(q, B->mod);
    nmod_poly_init_mod(b, B->mod);
    nmod_poly_init_mod(r, B->mod);
    nmod_poly_set(b, B);
    fq_nmod_poly_zero(A, ectx);
    for (i = 0; i < B->length; i++)
    {
        nmod_poly_divrem(q, r, b, ectx->modulus);
        fq_nmod_poly_set_coeff(A, i, r, ectx);
        nmod_poly_swap(b, q);
    }
    FLINT_ASSERT(nmod_poly_is_zero(b));
    nmod_poly_clear(b);
    nmod_poly_clear(q);
    nmod_poly_clear(r);
}

void nmod_poly_unshift_fq_nmod_poly(
    nmod_poly_t A,
    fq_nmod_poly_t B,
    const fq_nmod_ctx_t ectx)
{
    slong i;
    nmod_poly_zero(A);
    for (i = B->length - 1; i >= 0; i++)
    {
        nmod_poly_mul(A, A, ectx->modulus);
        nmod_poly_add(A, A, B->coeffs + i);
    }
}

void fq_nmod_bpoly_shift_nmod_bpoly(
    fq_nmod_bpoly_t A,
    nmod_bpoly_t B,
    const fq_nmod_ctx_t ectx)
{
    slong i;
    fq_nmod_bpoly_fit_length(A, B->length, ectx);
    A->length = 0;
    for (i = 0; i < B->length; i++)
    {
        fq_nmod_poly_shift_nmod_poly(A->coeffs + i, B->coeffs + i, ectx);
        if (!fq_nmod_poly_is_zero(A->coeffs + i, ectx))
            A->length = i + 1;
    }
}

void nmod_bpoly_unshift_fq_nmod_bpoly(
    nmod_bpoly_t A,
    fq_nmod_bpoly_t B,
    const fq_nmod_ctx_t ectx)
{
    slong i;
    nmod_bpoly_fit_length(A, B->length);
    A->length = 0;
    for (i = 0; i < B->length; i++)
    {
        nmod_poly_unshift_fq_nmod_poly(A->coeffs + i, B->coeffs + i, ectx);
        if (!nmod_poly_is_zero(A->coeffs + i))
            A->length = i + 1;
    }
}
*/



void fq_nmod_mpoly_set_nmod_mpoly(
    fq_nmod_mpoly_t A,
    const fq_nmod_mpoly_ctx_t Actx,
    const nmod_mpoly_t B,
    const nmod_mpoly_ctx_t Bctx)
{
    slong i, N;

    FLINT_ASSERT(Actx->minfo->ord == Bctx->minfo->ord);
    FLINT_ASSERT(Actx->minfo->nvars == Bctx->minfo->nvars);

    fq_nmod_mpoly_fit_bits(A, B->bits, Actx);
    A->bits = B->bits;

    N = mpoly_words_per_exp(B->bits, Bctx->minfo);

    fq_nmod_mpoly_fit_length(A, B->length, Actx);
    A->length = B->length;

    mpoly_copy_monomials(A->exps, B->exps, B->length, N);

    for (i = 0; i < B->length; i++)
        fq_nmod_set_ui(A->coeffs + i, B->coeffs[i], Actx->fqctx);
}



int nmod_mpoly_get_fq_nmod_mpoly(
    nmod_mpoly_t A,
    const nmod_mpoly_ctx_t Actx,
    const fq_nmod_mpoly_t B,
    const fq_nmod_mpoly_ctx_t Bctx)
{
    slong i, N;

    FLINT_ASSERT(Actx->minfo->ord == Bctx->minfo->ord);
    FLINT_ASSERT(Actx->minfo->nvars == Bctx->minfo->nvars);

    nmod_mpoly_fit_bits(A, B->bits, Actx);
    A->bits = B->bits;

    N = mpoly_words_per_exp(B->bits, Bctx->minfo);

    nmod_mpoly_fit_length(A, B->length, Actx);
    A->length = B->length;

    mpoly_copy_monomials(A->exps, B->exps, B->length, N);

    for (i = 0; i < B->length; i++)
    {
        if ((B->coeffs + i)->length == 1)
        {
            A->coeffs[i] = (B->coeffs + i)->coeffs[0];
        }
        else
        {
            return 0;
        }
    }

    return 1;
}


void get_conjugates(fq_nmod_mpoly_struct * C, const fq_nmod_mpoly_t A, slong deg, const fq_nmod_mpoly_ctx_t ctx)
{
    slong i, j;
    for (i = 0; i < deg; i++)
    {
        fq_nmod_mpoly_struct * Ci = C + i;
        fq_nmod_mpoly_set(Ci, A, ctx);
        for (j = 0; j < Ci->length; j++)
        {
            fq_nmod_frobenius(Ci->coeffs + j, Ci->coeffs + j, i, ctx->fqctx);
        }
    }
}


static int _fq_nmod_try_lift(
    fq_nmod_mpoly_factor_t qfac,
    const fq_nmod_mpoly_t q,
    const fq_nmod_mpoly_factor_t pfac,
    const fq_nmod_mpoly_t p,
    slong m,
    fq_nmod_struct * alpha,
    slong n,
    const fq_nmod_mpoly_ctx_t ctx)
{
    int success;
    slong i;
    slong * newdeg;
    fq_nmod_mpoly_t lcq, lcp, t, newq;
    fq_nmod_mpoly_univar_t u;

    FLINT_ASSERT(pfac->num > 1);

    newdeg = (slong *) flint_malloc((n + 1)*sizeof(slong));
    fq_nmod_mpoly_init(lcq, ctx);
    fq_nmod_mpoly_init(lcp, ctx);
    fq_nmod_mpoly_init(t, ctx);
    fq_nmod_mpoly_init(newq, ctx);
    fq_nmod_mpoly_univar_init(u, ctx);

    FLINT_ASSERT(fq_nmod_is_one(pfac->constant, ctx->fqctx));
    FLINT_ASSERT(fq_nmod_mpoly_factor_matches(p, pfac, ctx));

    _fq_nmod_mpoly_get_lc(lcq, q, ctx);
    fq_nmod_mpoly_evaluate_one_fq_nmod(lcp, lcq, m, alpha + m - 1, ctx);

    FLINT_ASSERT(lcp->length > 0);

    fq_nmod_mpoly_pow_ui(t, lcq, pfac->num - 1, ctx);
    fq_nmod_mpoly_mul(newq, q, t, ctx);
    fq_nmod_mpoly_degrees_si(newdeg, newq, ctx);

    fq_nmod_set(qfac->constant, pfac->constant, ctx->fqctx);
    fq_nmod_mpoly_factor_fit_length(qfac, pfac->num, ctx);
    qfac->num = pfac->num;
    for (i = 0; i < pfac->num; i++)
    {
        _fq_nmod_mpoly_get_lc(t, pfac->poly + i, ctx);
        success = fq_nmod_mpoly_divides(t, lcp, t, ctx);
        FLINT_ASSERT(success);
        fq_nmod_mpoly_mul(qfac->poly + i, pfac->poly + i, t, ctx);
        _fq_nmod_mpoly_set_lc(qfac->poly + i, qfac->poly + i, lcq, ctx);
        fmpz_one(qfac->exp + i);
    }

    success = fq_nmod_mpoly_hlift(m, qfac->poly, qfac->num, alpha, newq, newdeg, ctx);
    if (!success)
        goto cleanup;

    for (i = 0; i < qfac->num; i++)
    {
        fq_nmod_mpoly_to_univar(u, qfac->poly + i, 0, ctx);
        success = fq_nmod_mpoly_univar_content_mpoly(t, u, ctx);
        if (!success)
        {
            success = -1;
            goto cleanup;
        }
        success = fq_nmod_mpoly_divides(qfac->poly + i, qfac->poly + i, t, ctx);
        FLINT_ASSERT(success);
        fq_nmod_mpoly_make_monic(qfac->poly + i, qfac->poly + i, ctx);
    }

    success = 1;

cleanup:

    flint_free(newdeg);
    fq_nmod_mpoly_clear(lcq, ctx);
    fq_nmod_mpoly_clear(lcp, ctx);
    fq_nmod_mpoly_clear(t, ctx);
    fq_nmod_mpoly_clear(newq, ctx);
    fq_nmod_mpoly_univar_clear(u, ctx);

    /* q and its factors are primitive with positive lc */
    FLINT_ASSERT(!success || (fq_nmod_is_one(qfac->constant, ctx->fqctx) &&
                              fq_nmod_mpoly_factor_matches(q, qfac, ctx)));
    return success;
}



int nmod_mpoly_factor_irred_lgprime_default(
    nmod_mpolyv_t fac_,
    const nmod_mpoly_t A_,
    const nmod_mpoly_ctx_t ctx_)
{
    fq_nmod_mpoly_factor_t fac;
    fq_nmod_mpoly_t A;
    fq_nmod_mpoly_ctx_t ctx;
    slong edeg;
    int success;
    const slong n = ctx_->minfo->nvars - 1;
    slong i, j, m, r;
	fmpz_t subset;
    flint_rand_t randstate;
    fq_nmod_struct * alpha;
    fq_nmod_mpoly_struct * Aevals;
    slong * deg, * degeval;
    fq_nmod_mpoly_factor_t qfac, pfac, tfac, dfac;
    fq_nmod_mpoly_t t, p, q;
	fq_nmod_mpoly_univar_t u;
/*
flint_printf("_irreducible_mvar_factors_lgprime(n = %wd) called\n", n);
flint_printf("A_: "); nmod_mpoly_print_pretty(A_, NULL, ctx_); printf("\n");
*/
    FLINT_ASSERT(A_->length > 0);
    FLINT_ASSERT(A_->coeffs[0] == 1);
    FLINT_ASSERT(ctx_->minfo->ord == ORD_LEX);

    edeg = 2;
    fq_nmod_mpoly_ctx_init_deg(ctx, n + 1, ORD_LEX, ctx_->ffinfo->mod.n, edeg);
    fq_nmod_mpoly_init(A, ctx);
    fq_nmod_mpoly_set_nmod_mpoly(A, ctx, A_, ctx_);

    fq_nmod_mpoly_factor_init(fac, ctx);
    fq_nmod_mpoly_factor_one(fac, ctx);

	fmpz_init(subset);
    flint_randinit(randstate);
	alpha = (fq_nmod_struct *) flint_malloc(n*sizeof(fq_nmod_struct));
    for (i = 0; i < n; i++)
        fq_nmod_init(alpha + i, ctx->fqctx);

    deg     = (slong *) flint_malloc((n + 1)*sizeof(slong));
    degeval = (slong *) flint_malloc((n + 1)*sizeof(slong));

    Aevals    = (fq_nmod_mpoly_struct *) flint_malloc(n*sizeof(fq_nmod_mpoly_struct));
	for (i = 0; i < n; i++)
		fq_nmod_mpoly_init(Aevals + i, ctx);

    fq_nmod_mpoly_factor_init(pfac, ctx);
	fq_nmod_mpoly_factor_init(qfac, ctx);
    fq_nmod_mpoly_factor_init(tfac, ctx);
	fq_nmod_mpoly_factor_init(dfac, ctx);
	fq_nmod_mpoly_init(t, ctx);
	fq_nmod_mpoly_init(p, ctx);
	fq_nmod_mpoly_init(q, ctx);
	fq_nmod_mpoly_univar_init(u, ctx);

	fq_nmod_mpoly_degrees_si(deg, A, ctx);
    goto got_alpha;

next_alpha:

    edeg++;
    if (edeg > 1000)
    {
        success = -1;
        goto cleanup;
    }

    fq_nmod_mpoly_ctx_change_modulus(ctx, edeg);
    fq_nmod_mpoly_set_nmod_mpoly(A, ctx, A_, ctx_);

    for (i = 0; i < n; i++)
        fq_nmod_rand(alpha + i, randstate, ctx->fqctx);

got_alpha:

	/* ensure degrees do not drop under evalutaion */
    for (i = n - 1; i >= 0; i--)
    {
    	fq_nmod_mpoly_evaluate_one_fq_nmod(Aevals + i,
                       i == n - 1 ? A : Aevals + i + 1, i + 1, alpha + i, ctx);
    	fq_nmod_mpoly_degrees_si(degeval, Aevals + i, ctx);
	    for (j = 0; j <= i; j++)
        {
	    	if (degeval[j] != deg[j])
			    goto next_alpha;
        }
    }

    /* make sure our univar is squarefree */
	fq_nmod_mpoly_derivative(t, Aevals + 0, 0, ctx);
	fq_nmod_mpoly_gcd(t, t, Aevals + 0, ctx);
	if (!fq_nmod_mpoly_is_one(t, ctx))
		goto next_alpha;

    /* make our evaluations primitive */
    for (i = n - 1; i > 0; i--)
    {
        fq_nmod_mpoly_to_univar(u, Aevals + i, 0, ctx);
        fq_nmod_mpoly_univar_content_mpoly(t, u, ctx);
        success = fq_nmod_mpoly_divides(Aevals + i, Aevals + i, t, ctx);
        FLINT_ASSERT(success);
        fq_nmod_mpoly_make_monic(Aevals + i, Aevals + i, ctx);
    }

    success = _fq_nmod_irreducible_bivar_factors_smprime(pfac, Aevals + 1,
                                                                    0, 1, ctx);
    if (!success)
        goto next_alpha;

    for (m = 2; m <= n; m++)
    {
        fq_nmod_mpoly_set(q, m < n ? Aevals + m : A, ctx);
        fq_nmod_mpoly_set(p, Aevals + m - 1, ctx);

        FLINT_ASSERT(fq_nmod_is_one(pfac->constant, ctx->fqctx));
        FLINT_ASSERT(fq_nmod_mpoly_factor_matches(p, pfac, ctx));

        /* if p has only one factor, A must be irreducible */
        if (pfac->num < 2)
        {
	        fq_nmod_mpoly_factor_append_ui(fac, A, 1, ctx);
		    success = 1;
		    goto cleanup;
        }

        success = _fq_nmod_try_lift(qfac, q, pfac, p, m, alpha, n, ctx);
        if (success)
        {
            fq_nmod_mpoly_factor_swap(qfac, pfac, ctx);
            continue;
        }

        if (pfac->num == 2)
        {
            fq_nmod_mpoly_factor_append_ui(fac, A, 1, ctx);
            success = 1;
            goto cleanup;
        }

        qfac->num = 0;

try_again:

        for (r = 1; r <= pfac->num/2; r++)
        {
            subset_first(subset, pfac->num, r);
            FLINT_ASSERT(fq_nmod_is_one(pfac->constant, ctx->fqctx));
            FLINT_ASSERT(fq_nmod_mpoly_factor_matches(p, pfac, ctx));

            do {
                fq_nmod_mpoly_factor_fit_length(dfac, 2, ctx);
                dfac->num = 2;
                fq_nmod_one(dfac->constant, ctx->fqctx);
                fq_nmod_mpoly_one(dfac->poly + 0, ctx);
                fq_nmod_mpoly_one(dfac->poly + 1, ctx);
                fmpz_one(dfac->exp + 0);
                fmpz_one(dfac->exp + 1);
                for (i = 0; i < pfac->num; i++)
                {
                    j = fmpz_tstbit(subset, i);
                    fq_nmod_mpoly_mul(dfac->poly + j, dfac->poly + j,
                                                          pfac->poly + i, ctx);
                }

                success = _fq_nmod_try_lift(tfac, q, dfac, p, m, alpha, n, ctx);
                if (success)
                {
                    for (i = pfac->num - 1; i >= 0; i--)
                    {
                        if (fmpz_tstbit(subset, i))
                        {
                            fq_nmod_mpoly_swap(pfac->poly + i,
                                              pfac->poly + pfac->num - 1, ctx);
                            pfac->num--;                            
                        }
                    }
                    fq_nmod_mpoly_factor_append_ui(qfac, tfac->poly + 1, 1, ctx);
                    fq_nmod_mpoly_swap(q, tfac->poly + 0, ctx);
                    fq_nmod_mpoly_swap(p, dfac->poly + 0, ctx);
                    goto try_again;
                }
            }
            while (subset_next(subset, subset, pfac->num));
        }
        /* if pfac could not be combined, p must be irreducible */
        fq_nmod_mpoly_factor_append_ui(qfac, q, 1, ctx);
        fq_nmod_mpoly_factor_swap(qfac, pfac, ctx);
    }

    success = 1;

    for (i = 0; i < pfac->num; i++)
        fq_nmod_mpoly_factor_append_ui(fac, pfac->poly + i, 1, ctx);

cleanup:

    fac_->length = 0;
    if (success)
	{
        nmod_mpoly_t truefactor_;
        fq_nmod_mpoly_t truefactor;
        fq_nmod_mpoly_struct * conjugates;
/*
printf("now must frob combine\n");
*/
	    fq_nmod_mpoly_set_nmod_mpoly(A, ctx, A_, ctx_);
/*
flint_printf("A: "); fq_nmod_mpoly_print_pretty(A, NULL, ctx); printf("\n");
flint_printf("fac: "); fq_nmod_mpoly_factor_print_pretty(fac, NULL, ctx); printf("\n");
*/
		FLINT_ASSERT(fq_nmod_mpoly_factor_matches(A, fac, ctx));

        conjugates = (fq_nmod_mpoly_struct *) flint_malloc(edeg*sizeof(fq_nmod_mpoly_struct));
        for (i = 0; i < edeg; i++)
            fq_nmod_mpoly_init(conjugates + i, ctx);

        nmod_mpoly_init(truefactor_, ctx_);
        fq_nmod_mpoly_init(truefactor, ctx);

        while (fac->num > 0)
        {
            fq_nmod_mpoly_one(truefactor, ctx);
            get_conjugates(conjugates, fac->poly + 0, edeg, ctx);
/*
flint_printf("fac->poly[0]: "); fq_nmod_mpoly_print_pretty(fac->poly + 0, NULL, ctx);printf("\n");
for (i = 0; i < edeg; i++)
{
flint_printf("conjugate[%wd]: ", i); fq_nmod_mpoly_print_pretty(conjugates + i, NULL, ctx);printf("\n");
}
*/
            for (i = 0; i < fac->num; i++)
            {
                for (j = 0; j < edeg; j++)
                {
                    if (fq_nmod_mpoly_equal(fac->poly + i, conjugates + j, ctx))
                    {
/*
flint_printf("match i = %wd, j = %wd\n", i, j);
*/
                        fq_nmod_mpoly_mul(truefactor, truefactor, fac->poly + i, ctx);
                        fq_nmod_mpoly_swap(fac->poly + i, fac->poly + fac->num - 1, ctx);
                        fac->num--;
                        i--;
						break;
                    }
                }
            }

            success = nmod_mpoly_get_fq_nmod_mpoly(truefactor_, ctx_, truefactor, ctx);
            FLINT_ASSERT(success);

            nmod_mpolyv_fit_length(fac_, fac_->length + 1, ctx_);
            nmod_mpoly_swap(fac_->coeffs + fac_->length, truefactor_, ctx_);
            fac_->length++;
        }

        for (i = 0; i < edeg; i++)
            fq_nmod_mpoly_clear(conjugates + i, ctx);
        flint_free(conjugates);

        fq_nmod_mpoly_clear(truefactor, ctx);
        nmod_mpoly_clear(truefactor_, ctx_);

        success = 1;
    }

    fq_nmod_mpoly_clear(A, ctx);

    fq_nmod_mpoly_factor_clear(fac, ctx);

	fmpz_clear(subset);
    flint_randclear(randstate);

    for (i = 0; i < n; i++)
        fq_nmod_clear(alpha + i, ctx->fqctx);
    flint_free(alpha);

    flint_free(deg);
    flint_free(degeval);

	for (i = 0; i < n; i++)
		fq_nmod_mpoly_clear(Aevals + i, ctx);
    flint_free(Aevals);

    fq_nmod_mpoly_factor_clear(pfac, ctx);
	fq_nmod_mpoly_factor_clear(qfac, ctx);
    fq_nmod_mpoly_factor_clear(tfac, ctx);
	fq_nmod_mpoly_factor_clear(dfac, ctx);
	fq_nmod_mpoly_clear(t, ctx);
	fq_nmod_mpoly_clear(p, ctx);
	fq_nmod_mpoly_clear(q, ctx);
	fq_nmod_mpoly_univar_clear(u, ctx);

    fq_nmod_mpoly_ctx_clear(ctx);

    return success;
}


void nmod_mpoly_convert_perm(
    nmod_mpoly_t A,
    flint_bitcnt_t Abits,
    const nmod_mpoly_ctx_t lctx,
    const nmod_mpoly_t B,
    const nmod_mpoly_ctx_t ctx,
    const slong * perm)
{
    slong n = ctx->minfo->nvars;
    slong m = lctx->minfo->nvars;
    slong i, k, l;
    slong NA, NB;
    ulong * Aexps;
    ulong * Bexps;
    TMP_INIT;

    FLINT_ASSERT(B->length > 0);
    FLINT_ASSERT(Abits <= FLINT_BITS);
    FLINT_ASSERT(B->bits <= FLINT_BITS);
    TMP_START;

    Aexps = (ulong *) TMP_ALLOC(m*sizeof(ulong));
    Bexps = (ulong *) TMP_ALLOC(n*sizeof(ulong));

    NA = mpoly_words_per_exp(Abits, lctx->minfo);
    NB = mpoly_words_per_exp(B->bits, ctx->minfo);

    nmod_mpoly_fit_bits(A, Abits, ctx);
    A->bits = Abits;
    nmod_mpoly_fit_length(A, B->length, lctx);
    A->length = B->length;
    for (i = 0; i < B->length; i++)
    {        
	    A->coeffs[i] = B->coeffs[i];
	    mpoly_get_monomial_ui(Bexps, B->exps + NB*i, B->bits, ctx->minfo);
	    for (k = 0; k < m; k++)
	    {
	        l = perm[k];
	        Aexps[k] = l < 0 ? 0 : Bexps[l];
	    }
	    mpoly_set_monomial_ui(A->exps + NA*(i), Aexps, Abits, lctx->minfo);
     }  
    TMP_END;
    nmod_mpoly_sort_terms(A, lctx);
}

/*
    A is primitive w.r.t to any variable appearing in A.
    A is squarefree and monic.

    return 1 for success, 0 for failure
*/
static int _irreducible_factors(
    nmod_mpolyv_t Af,
    nmod_mpoly_t A,
    const nmod_mpoly_ctx_t ctx)
{
    int success;
    slong i, t, nzdvar, mvars;
    slong nvars = ctx->minfo->nvars;
    slong * Adegs, * perm, * iperm;
    nmod_mpoly_t G, Abar, Bbar, nzdpoly;
    flint_bitcnt_t Lbits, Abits;
    int perm_is_id;
#if WANT_ASSERT
    nmod_mpoly_t Aorg;

    nmod_mpoly_init(Aorg, ctx);
    nmod_mpoly_set(Aorg, A, ctx);
#endif

    nmod_mpoly_init(G, ctx);
    nmod_mpoly_init(Abar, ctx);
    nmod_mpoly_init(Bbar, ctx);
    nmod_mpoly_init(nzdpoly, ctx);

    Adegs = (slong *) flint_malloc(3*nvars*sizeof(slong));
    perm = Adegs + nvars;
    iperm = perm + nvars;

    if (A->length < 2)
    {
        FLINT_ASSERT(A->length == 1);
        FLINT_ASSERT(!nmod_mpoly_is_ui(A, ctx));
        nmod_mpolyv_fit_length(Af, 1, ctx);
        Af->length = 1;
        nmod_mpoly_swap(Af->coeffs + 0, A, ctx);
        success = 1;
        goto cleanup;
    }

    if (!nmod_mpoly_degrees_fit_si(A, ctx))
    {
        success = 0;
        goto cleanup;
    }

    if (A->bits > FLINT_BITS &&
        !nmod_mpoly_repack_bits_inplace(A, FLINT_BITS, ctx))
    {
        success = 0;
        goto cleanup;
    }

    nmod_mpoly_degrees_si(Adegs, A, ctx);

    Abits = A->bits;

    mvars = 0;
    Lbits = 0;
    nzdvar = -1;
    for (i = 0; i < nvars; i++)
    {
		iperm[i] = -1;
        if (Adegs[i] > 0)
        {
            flint_bitcnt_t this_bits = FLINT_BIT_COUNT(Adegs[i]);
            Lbits = FLINT_MAX(Lbits, this_bits);
            perm[mvars] = i;
            mvars++;

            if (nzdvar < 0)
            {
                nmod_mpoly_derivative(nzdpoly, A, i, ctx);
                if (!nmod_mpoly_is_zero(nzdpoly, ctx))
                {
                    nzdvar = i;
                    t = perm[mvars - 1];
                    perm[mvars - 1] = perm[0];
                    perm[0] = t;
                }
            }
        }
    }

    FLINT_ASSERT(nzdvar >= 0);
    FLINT_ASSERT(perm[0] == nzdvar);
    FLINT_ASSERT(!nmod_mpoly_is_zero(nzdpoly, ctx));

    /*
        Check for annoying things like (x^2+y)(x^p+y). The squarefree
        requirement should already have ruled out (x^2+y)(x^p+y^p).
    */
    if (!nmod_mpoly_gcd_cofactors(G, Abar, Bbar, A, nzdpoly, ctx))
    {
        success = 0;
        goto cleanup;
    }
    if (!nmod_mpoly_is_one(G, ctx))
    {
        nmod_mpolyv_t newf;
        nmod_mpolyv_init(newf, ctx);
        success = _irreducible_factors(Af, Abar, ctx);
        success = success && _irreducible_factors(newf, Bbar, ctx);
        if (success)
        {
            nmod_mpolyv_fit_length(Af, Af->length + newf->length, ctx);
            for (i = 0; i < newf->length; i++)
            {
                nmod_mpoly_swap(Af->coeffs + Af->length, newf->coeffs + i, ctx);
                Af->length++;
            }
        }
        nmod_mpolyv_clear(newf, ctx);
        goto cleanup;
    }

	/* A is separable wrt gen(perm[0]) now */

    if (Lbits > FLINT_BITS - 10)
    {
        success = 0;
        goto cleanup;
    }

    /* TODO nice permutation */

    /* invert perm */
    perm_is_id = (mvars == nvars);
    for (i = 0; i < mvars; i++)
    {
        perm_is_id = perm_is_id && (perm[i] == i);
        iperm[perm[i]] = i;
    }

    if (mvars < 2)
    {
        nmod_poly_t Au;
        nmod_poly_factor_t Auf;

        FLINT_ASSERT(mvars == 1);

        nmod_poly_init_mod(Au, ctx->ffinfo->mod);
        nmod_poly_factor_init(Auf);

        FLINT_ASSERT(nmod_mpoly_is_nmod_poly(A, perm[0], ctx));
        success = nmod_mpoly_get_nmod_poly(Au, A, perm[0], ctx);
        FLINT_ASSERT(success);
        nmod_poly_factor(Auf, Au);

        nmod_mpolyv_fit_length(Af, Auf->num, ctx);
        Af->length = Auf->num; 
        for (i = 0; i < Auf->num; i++)
        {
            FLINT_ASSERT(Auf->exp[i] == 1);
            _nmod_mpoly_set_nmod_poly(Af->coeffs + i, Abits,
                             Auf->p[i].coeffs, Auf->p[i].length, perm[0], ctx);
        }

        nmod_poly_clear(Au);
        nmod_poly_factor_clear(Auf);

        success = 1;
    }
    else if (mvars == 2)
    {
        n_poly_t c;
        n_bpoly_t Ab;
        n_tpoly_t Abf;

        n_poly_init(c);
        n_bpoly_init(Ab);
        n_tpoly_init(Abf);

        nmod_mpoly_get_bpoly(Ab, A, perm[0], perm[1], ctx);
        success = n_bpoly_mod_factor_smprime(c, Abf, Ab, 1, ctx->ffinfo->mod);
        if (!success)
        {
            nmod_mpoly_get_bpoly(Ab, A, perm[0], perm[1], ctx);
            n_bpoly_mod_factor_lgprime(c, Abf, Ab, ctx->ffinfo->mod);
        }

        FLINT_ASSERT(n_poly_degree(c) == 0);

        nmod_mpolyv_fit_length(Af, Abf->length, ctx);
        Af->length = Abf->length;
        for (i = 0; i < Abf->length; i++)
        {
            nmod_mpoly_set_bpoly(Af->coeffs + i, Abits, Abf->coeffs + i,
                                                        perm[0], perm[1], ctx);
            nmod_mpoly_make_monic(Af->coeffs + i, Af->coeffs + i, ctx);
        }

        n_poly_clear(c);
        n_bpoly_clear(Ab);
        n_tpoly_clear(Abf);

        success = 1;
    }
    else
    {
		nmod_mpoly_ctx_t Lctx;
		nmod_mpoly_t L;
		nmod_mpolyv_t Lf;

		nmod_mpoly_ctx_init(Lctx, mvars, ORD_LEX, ctx->ffinfo->mod.n);
		nmod_mpoly_init(L, Lctx);
		nmod_mpolyv_init(Lf, Lctx);

        Lbits = mpoly_fix_bits(Lbits + 1, Lctx->minfo);

		nmod_mpoly_convert_perm(L, Lbits, Lctx, A, ctx, perm);
        nmod_mpoly_make_monic(L, L, ctx);

		success = nmod_mpoly_factor_irred_smprime_default(Lf, L, Lctx);
        if (!success)
        {
    		success = nmod_mpoly_factor_irred_lgprime_default(Lf, L, Lctx);
        }

		if (success)
        {
            nmod_mpolyv_fit_length(Af, Lf->length, ctx);
            Af->length = Lf->length;
		    for (i = 0; i < Lf->length; i++)
		    {
                nmod_mpoly_convert_perm(Af->coeffs + i, Abits, ctx,
                                                  Lf->coeffs + i, Lctx, iperm);
                nmod_mpoly_make_monic(Af->coeffs + i, Af->coeffs + i, ctx);
		    }
        }

	    nmod_mpolyv_clear(Lf, Lctx);
	    nmod_mpoly_clear(L, Lctx);
	    nmod_mpoly_ctx_clear(Lctx);
    }

cleanup:

    nmod_mpoly_clear(G, ctx);
    nmod_mpoly_clear(Abar, ctx);
    nmod_mpoly_clear(Bbar, ctx);
    nmod_mpoly_clear(nzdpoly, ctx);

    flint_free(Adegs);

#if WANT_ASSERT
    if (success)
    {
        nmod_mpoly_t prod;
        nmod_mpoly_init(prod, ctx);
        nmod_mpoly_one(prod, ctx);
        for (i = 0; i < Af->length; i++)
            nmod_mpoly_mul(prod, prod, Af->coeffs + i, ctx);
        FLINT_ASSERT(nmod_mpoly_equal(prod, Aorg, ctx));
        nmod_mpoly_clear(prod, ctx);
        nmod_mpoly_clear(Aorg, ctx);
    }
#endif

    return success;
}


int nmod_mpoly_factor(
    nmod_mpoly_factor_t f,
    const nmod_mpoly_t A,
    const nmod_mpoly_ctx_t ctx)
{
    int success;
    slong i, j;
    nmod_mpolyv_t t;
    nmod_mpoly_factor_t g;

    nmod_mpolyv_init(t, ctx);
    nmod_mpoly_factor_init(g, ctx);

    success = nmod_mpoly_factor_squarefree(f, A, ctx);
    if (!success)
        goto cleanup;

    g->constant = f->constant;
    g->num = 0;
    for (j = 0; j < f->num; j++)
    {
        success = _irreducible_factors(t, f->poly + j, ctx);
        if (!success)
        {
            flint_printf("it failed\n");
            goto cleanup;
        }

        nmod_mpoly_factor_fit_length(g, g->num + t->length, ctx);
        for (i = 0; i < t->length; i++)
        {
            fmpz_set(g->exp + g->num, f->exp + j);
            nmod_mpoly_swap(g->poly + g->num, t->coeffs + i, ctx);
            g->num++;
        }
    }

    nmod_mpoly_factor_swap(f, g, ctx);

    success = 1;

cleanup:

    nmod_mpolyv_clear(t, ctx);
    nmod_mpoly_factor_clear(g, ctx);

    FLINT_ASSERT(!success || nmod_mpoly_factor_matches(A, f, ctx));

    return success;
}
