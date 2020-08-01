/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fq_nmod_mpoly_factor.h"


void fq_nmod_bpoly_clear(fq_nmod_bpoly_t A, const fq_nmod_ctx_t ctx)
{
    slong i;
    if (A->alloc > 0)
    {
        for (i = 0; i < A->alloc; i++)
            fq_nmod_poly_clear(A->coeffs + i, ctx);
        flint_free(A->coeffs);
    }
    else
    {
        FLINT_ASSERT(A->coeffs == NULL);
    }
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


void fq_nmod_bpoly_print_pretty(
    const fq_nmod_bpoly_t A,
    const char * var0,
    const char * var1,
    const fq_nmod_ctx_t ctx)
{
    slong i;
    int first;

    first = 1;
    for (i = A->length - 1; i >= 0; i--)
    {
        if (fq_nmod_poly_is_zero(A->coeffs + i, ctx))
            continue;

        if (!first)
            flint_printf(" + ");

        first = 0;

        flint_printf("(");
        fq_nmod_poly_print_pretty(A->coeffs + i, var1, ctx);
        flint_printf(")*%s^%wd", var0, i);
    }

    if (first)
        flint_printf("0");
}

slong fq_nmod_bpoly_degree1(const fq_nmod_bpoly_t A, const fq_nmod_ctx_t ctx)
{
    slong i, len = 0;
    for (i = 0; i < A->length; i++)
        len = FLINT_MAX(len, A->coeffs[i].length);
    return len - 1;    
}

int fq_nmod_bpoly_is_canonical(const fq_nmod_bpoly_t A, const fq_nmod_ctx_t ctx)
{
    if (A->length <= 0)
        return A->length == 0;
/*
    slong i;

    for (i = 0; i < A->length; i++)
    {
        if (!fq_nmod_poly_is_canonical(A->coeffs + i, ctx))
            return 0;
    }
*/
    return !fq_nmod_poly_is_zero(A->coeffs + A->length - 1, ctx);
}

int fq_nmod_bpoly_equal(
    const fq_nmod_bpoly_t A,
    const fq_nmod_bpoly_t B,
    const fq_nmod_ctx_t ctx)
{
    slong i;

    if (A->length != B->length)
        return 0;

    for (i = 0; i < A->length; i++)
    {
        if (!fq_nmod_poly_equal(A->coeffs + i, B->coeffs + i, ctx))
            return 0;
    }

    return 1;
}

void fq_nmod_bpoly_set_coeff(
    fq_nmod_bpoly_t A,
    slong xi,
    slong yi,
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
    fq_nmod_bpoly_normalise(A, ctx);
}


void fq_nmod_bpoly_set_poly_var0(
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

void fq_nmod_bpoly_set_poly_var1(
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


void fq_nmod_bpoly_mul(
    fq_nmod_bpoly_t A,
    const fq_nmod_bpoly_t B,
    const fq_nmod_bpoly_t C,
    const fq_nmod_ctx_t ctx)
{
    slong i, j;
    fq_nmod_poly_struct * t;

    FLINT_ASSERT(A != B);
    FLINT_ASSERT(A != C);

    if (B->length < 1 || C->length < 1)
    {
        A->length = 0;
        return;
    }

    fq_nmod_bpoly_fit_length(A, B->length + C->length, ctx);
    for (i = 0; i < B->length + C->length - 1; i++)
        fq_nmod_poly_zero(A->coeffs + i, ctx);

    t = A->coeffs + B->length + C->length - 1;

    for (i = 0; i < B->length; i++)
    {
        for (j = 0; j < C->length; j++)
        {
            fq_nmod_poly_mul(t, B->coeffs + i, C->coeffs + j, ctx);
            fq_nmod_poly_add(A->coeffs + i + j, A->coeffs + i + j, t, ctx);
        }
    }

    A->length = B->length + C->length - 1;
    fq_nmod_bpoly_normalise(A, ctx);
}


void fq_nmod_bpoly_mul_series(
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


void fq_nmod_bpoly_add(
    fq_nmod_bpoly_t A,
    const fq_nmod_bpoly_t B,
    const fq_nmod_bpoly_t C,
    const fq_nmod_ctx_t ctx)
{
    slong i;
    slong Alen = FLINT_MAX(B->length, C->length);

    FLINT_ASSERT(A != B);
    FLINT_ASSERT(A != C);

    fq_nmod_bpoly_fit_length(A, Alen, ctx);

    A->length = 0;
    for (i = 0; i < Alen; i++)
    {
        if (i < B->length)
        {
            if (i < C->length)
            {
                fq_nmod_poly_add(A->coeffs + i, B->coeffs + i, C->coeffs + i, ctx);
            }
            else
            {
                fq_nmod_poly_set(A->coeffs + i, B->coeffs + i, ctx);
            }
        }
        else
        {
            FLINT_ASSERT(i < C->length);
            fq_nmod_poly_set(A->coeffs + i, C->coeffs + i, ctx);
        }

        if (!fq_nmod_poly_is_zero(A->coeffs + i, ctx))
            A->length = i + 1;
    }
}

void fq_nmod_bpoly_one(fq_nmod_bpoly_t A, const fq_nmod_ctx_t ctx)
{
    fq_nmod_bpoly_fit_length(A, 1, ctx);
    A->length = 1;
    fq_nmod_poly_one(A->coeffs + 0, ctx);
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

void fq_nmod_bpoly_derivative(
    fq_nmod_bpoly_t A,
    const fq_nmod_bpoly_t B,
    const fq_nmod_ctx_t ctx)
{
    slong i;
    fq_nmod_t c;

    if (B->length < 2)
    {
        A->length = 0;
        return;
    }

    fq_nmod_init(c, ctx);

    fq_nmod_bpoly_fit_length(A, B->length - 1, ctx);

    for (i = 1; i < B->length; i++)
    {
        fq_nmod_set_si(c, i, ctx);
        fq_nmod_poly_scalar_mul_fq_nmod(A->coeffs + i - 1, B->coeffs + i, c, ctx);
    }

    A->length = B->length - 1;
    fq_nmod_bpoly_normalise(A, ctx);

    fq_nmod_clear(c, ctx);
}

/*
    division in (Fq[y]/y^order)[x]
    inputs need not be reduced mod y^order
*/
void fq_nmod_bpoly_divrem_series(
    fq_nmod_bpoly_t Q,
    fq_nmod_bpoly_t R,
    const fq_nmod_bpoly_t A,
    const fq_nmod_bpoly_t B,
    slong order,
    const fq_nmod_ctx_t ctx)
{
    slong i, qoff;
    fq_nmod_poly_t q, t;

    FLINT_ASSERT(R != A);
    FLINT_ASSERT(R != B);
    FLINT_ASSERT(Q != A);
    FLINT_ASSERT(Q != B);
    FLINT_ASSERT(B->length > 0);

    fq_nmod_poly_init(q, ctx);
    fq_nmod_poly_init(t, ctx);

    fq_nmod_bpoly_set(R, A, ctx);
    for (i = 0; i < R->length; i++)
        fq_nmod_poly_truncate(R->coeffs + i, order, ctx);
    fq_nmod_bpoly_normalise(R, ctx);

    Q->length = 0;

    while (R->length >= B->length)
    {
        fq_nmod_poly_div_series(q, R->coeffs + R->length - 1,
                                        B->coeffs + B->length - 1, order, ctx);

        for (i = 0; i < B->length; i++)
        {
            fq_nmod_poly_mullow(t, B->coeffs + i, q, order, ctx);
            fq_nmod_poly_sub(R->coeffs + i + R->length - B->length,
                                R->coeffs + i + R->length - B->length, t, ctx);
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

        FLINT_ASSERT(fq_nmod_poly_is_zero(R->coeffs + R->length - 1, ctx));

        fq_nmod_bpoly_normalise(R, ctx);
    }

    fq_nmod_poly_clear(q, ctx);
    fq_nmod_poly_clear(t, ctx);
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
        fq_nmod_poly_divrem(q, t, R->coeffs + R->length - 1,
                                               B->coeffs + B->length - 1, ctx);
        if (!fq_nmod_poly_is_zero(t, ctx))
        {
            divides = 0;
            goto cleanup;
        }

        for (i = 0; i < B->length; i++)
        {
            fq_nmod_poly_mul(t, B->coeffs + i, q, ctx);
            fq_nmod_poly_sub(R->coeffs + i + R->length - B->length,
                                R->coeffs + i + R->length - B->length, t, ctx);
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

        fq_nmod_bpoly_normalise(R, ctx);
    }

    divides = (R->length == 0);

cleanup:

    fq_nmod_poly_clear(q, ctx);
    fq_nmod_poly_clear(t, ctx);
    fq_nmod_bpoly_clear(R, ctx);

    return divides;
}

void fq_nmod_bpoly_set(fq_nmod_bpoly_t A, const fq_nmod_bpoly_t B,
                                                       const fq_nmod_ctx_t ctx)
{
    slong i;

    if (A == B)
        return;

    fq_nmod_bpoly_fit_length(A, B->length, ctx);
    A->length = B->length;

    for (i = 0; i < B->length; i++)
        fq_nmod_poly_set(A->coeffs + i, B->coeffs + i, ctx);
}

void fq_nmod_bpoly_make_primitive(
    fq_nmod_poly_t g,
    fq_nmod_bpoly_t A,
    const fq_nmod_ctx_t ctx)
{
    slong Alen = A->length;
    slong i;
    fq_nmod_poly_t q, r;
    fq_nmod_t c;

    fq_nmod_poly_init(q, ctx);
    fq_nmod_poly_init(r, ctx);
    fq_nmod_init(c, ctx);

    fq_nmod_poly_zero(g, ctx);
    for (i = 0; i < Alen; i++)
	{
        fq_nmod_poly_gcd(q, g, A->coeffs + i, ctx);
		fq_nmod_poly_swap(g, q, ctx);
	}

    for (i = 0; i < A->length; i++)
    {
        fq_nmod_poly_divrem(q, r, A->coeffs + i, g, ctx);
		fq_nmod_poly_swap(A->coeffs + i, q, ctx);
    }

    /* make lc_xy(A) one */
    if (Alen > 0)
    {
        fq_nmod_struct * c_ = A->coeffs[Alen - 1].coeffs
                                              + A->coeffs[Alen - 1].length - 1;
        if (!fq_nmod_is_one(c_, ctx))
        {
            fq_nmod_poly_scalar_mul_fq_nmod(g, g, c_, ctx);
            fq_nmod_inv(c, c_, ctx);
            for (i = 0; i < Alen; i++)
                fq_nmod_poly_scalar_mul_fq_nmod(A->coeffs + i,
                                                A->coeffs + i, c, ctx);
        }
    }


    fq_nmod_poly_clear(q, ctx);
    fq_nmod_poly_clear(r, ctx);
    fq_nmod_clear(c, ctx);
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


void fq_nmod_bpoly_taylor_shift_var1(
    fq_nmod_bpoly_t A,
    const fq_nmod_bpoly_t B,
    const fq_nmod_t alpha,
    const fq_nmod_ctx_t ctx)
{
    slong i;
    fq_nmod_bpoly_set(A, B, ctx);
    for (i = A->length - 1; i >= 0; i--)
        fq_nmod_poly_taylor_shift_horner(A->coeffs + i, A->coeffs + i, alpha, ctx);    
}

void fq_nmod_bpoly_taylor_shift_var0(
    fq_nmod_bpoly_t A,
    const fq_nmod_t alpha,
    const fq_nmod_ctx_t ctx)
{
    slong n, i, j;
    fq_nmod_poly_t t;

    if (fq_nmod_is_zero(alpha, ctx))
        return;

    fq_nmod_poly_init(t, ctx);
    n = A->length;

    for (i = n - 2; i >= 0; i--)
    {
        for (j = i; j < n - 1; j++)
        {
            fq_nmod_poly_scalar_mul_fq_nmod(t, A->coeffs + j + 1, alpha, ctx);
            fq_nmod_poly_add(A->coeffs + j, A->coeffs + j, t, ctx);
        }
    }

    fq_nmod_poly_clear(t, ctx);
}


void fq_nmod_mpoly_get_fq_nmod_bpoly(
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


void fq_nmod_mpoly_set_fq_nmod_bpoly(
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


void fq_nmod_mpoly_get_bpoly(
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


void fq_nmod_mpoly_set_bpoly(
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

    fq_nmod_mpoly_fit_length_set_bits(A, 8, Abits, ctx);

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

