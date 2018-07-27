/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "nmod_mpoly.h"

void fq_nmod_mpolyun_init(fq_nmod_mpolyun_t A, mp_bitcnt_t bits, const fq_nmod_mpoly_ctx_t ctx)
{
    A->coeffs = NULL;
    A->exps = NULL;
    A->alloc = 0;
    A->length = 0;
    A->bits = bits;
}

void fq_nmod_mpolyun_clear(fq_nmod_mpolyun_t A, const fq_nmod_mpoly_ctx_t ctx)
{
    slong i;
    for (i = 0; i < A->alloc; i++)
        fq_nmod_mpolyn_clear(A->coeffs + i, ctx);
    flint_free(A->coeffs);
    flint_free(A->exps);
}


void fq_nmod_mpolyun_swap(fq_nmod_mpolyun_t A, fq_nmod_mpolyun_t B)
{
   fq_nmod_mpolyun_struct t = *A;
   *A = *B;
   *B = t;
}

void fq_nmod_mpolyun_zero(fq_nmod_mpolyun_t A, const fq_nmod_mpoly_ctx_t ctx)
{
    slong i;
    for (i = 0; i < A->alloc; i++) {
        fq_nmod_mpolyn_clear(A->coeffs + i, ctx);
        fq_nmod_mpolyn_init(A->coeffs + i, A->bits, ctx);
    }
    A->length = 0;
}

void fq_nmod_mpolyun_print_pretty(const fq_nmod_mpolyun_t poly,
                                   const char ** x, const fq_nmod_mpoly_ctx_t ctx)
{
    slong i;
    if (poly->length == 0)
        flint_printf("0");
    for (i = 0; i < poly->length; i++)
    {
        if (i != 0)
            flint_printf(" + ");
        flint_printf("(");
        FLINT_ASSERT((poly->coeffs + i)->bits == poly->bits);
        fq_nmod_mpolyn_print_pretty(poly->coeffs + i,x,ctx);
        flint_printf(")*X^%wd",poly->exps[i]);
    }
}

void fq_nmod_mpolyun_fit_length(fq_nmod_mpolyun_t A, slong length, const fq_nmod_mpoly_ctx_t ctx)
{
    slong i;
    slong old_alloc = A->alloc;
    slong new_alloc = FLINT_MAX(length, 2*A->alloc);

    if (length > old_alloc)
    {
        if (old_alloc == 0)
        {
            A->exps = (ulong *) flint_malloc(new_alloc*sizeof(ulong));
            A->coeffs = (fq_nmod_mpolyn_struct *) flint_malloc(
                                      new_alloc*sizeof(fq_nmod_mpolyn_struct));
        } else
        {
            A->exps = (ulong *) flint_realloc(A->exps,
                                                      new_alloc*sizeof(ulong));
            A->coeffs = (fq_nmod_mpolyn_struct *) flint_realloc(A->coeffs,
                                      new_alloc*sizeof(fq_nmod_mpolyn_struct));
        }

        for (i = old_alloc; i < new_alloc; i++)
        {
            fq_nmod_mpolyn_init(A->coeffs + i, A->bits, ctx);
        }
        A->alloc = new_alloc;
    }
}

void fq_nmod_mpolyun_shift_right(fq_nmod_mpolyun_t A, ulong s)
{
    slong i;
    for (i = 0; i < A->length; i++)
    {
        FLINT_ASSERT(A->exps[i] >= s);
        A->exps[i] -= s;
    }
}

void fq_nmod_mpolyun_shift_left(fq_nmod_mpolyun_t A, ulong s)
{
    slong i;
    for (i = 0; i < A->length; i++)
    {
        A->exps[i] += s;
    }
}

slong fq_nmod_mpolyun_lastdeg(fq_nmod_mpolyun_t A, fq_nmod_mpoly_ctx_t ctx)
{
    slong i, j;
    slong deg = -WORD(1);

    for (i = 0; i < A->length; i++)
    {
        for (j = 0; j < (A->coeffs + i)->length; j++)
        {
            deg = FLINT_MAX(deg, fq_nmod_poly_degree((A->coeffs + i)->coeffs + j, ctx->fqctx));
        }
    }
    FLINT_ASSERT(deg >= 0);
    return deg;
}

void fq_nmod_mpolyun_set(fq_nmod_mpolyun_t A, const fq_nmod_mpolyun_t B, const fq_nmod_mpoly_ctx_t ctx)
{
    slong i, Blen;
    fq_nmod_mpolyn_struct * Acoeff, * Bcoeff;
    ulong * Aexp, * Bexp;

    Blen = B->length;
    fq_nmod_mpolyun_fit_length(A, Blen, ctx);
    Acoeff = A->coeffs;
    Bcoeff = B->coeffs;
    Aexp = A->exps;
    Bexp = B->exps;

    for (i = 0; i < Blen; i++)
    {
        fq_nmod_mpolyn_set(Acoeff + i, Bcoeff + i, ctx);
        Aexp[i] = Bexp[i];
    }

    /* demote remaining coefficients */
    for (i = Blen; i < A->length; i++)
    {
        fq_nmod_mpolyn_clear(Acoeff + i, ctx);
        fq_nmod_mpolyn_init(Acoeff + i, A->bits, ctx);
    }
    A->length = Blen;
}


void fq_nmod_mpolyun_mul_poly(
    fq_nmod_mpolyun_t A,
    const fq_nmod_mpolyun_t B,
    const fq_nmod_poly_t c,
    const fq_nmod_mpoly_ctx_t ctx)
{
    slong i, Blen;
    fq_nmod_mpolyn_struct * Acoeff, * Bcoeff;
    ulong * Aexp, * Bexp;

    Blen = B->length;
    fq_nmod_mpolyun_fit_length(A, Blen, ctx);
    Acoeff = A->coeffs;
    Bcoeff = B->coeffs;
    Aexp = A->exps;
    Bexp = B->exps;

    for (i = 0; i < Blen; i++)
    {
        fq_nmod_mpolyn_mul_poly(Acoeff + i, Bcoeff + i, c, ctx);
        Aexp[i] = Bexp[i];
    }

    /* demote remaining coefficients */
    for (i = Blen; i < A->length; i++)
    {
        fq_nmod_mpolyn_clear(Acoeff + i, ctx);
        fq_nmod_mpolyn_init(Acoeff + i, A->bits, ctx);
    }
    A->length = Blen;
}


void fq_nmod_mpolyun_content_last(fq_nmod_poly_t a, fq_nmod_mpolyun_t B,
                                                          fq_nmod_mpoly_ctx_t ctx)
{
    slong i, j;

    fq_nmod_poly_zero(a, ctx->fqctx);
    for (i = 0; i < B->length; i++)
    {
        for (j = 0; j < (B->coeffs + i)->length; j++)
        {
            fq_nmod_poly_gcd(a, a, (B->coeffs + i)->coeffs + j, ctx->fqctx);
            if (fq_nmod_poly_degree(a, ctx->fqctx) == 0)
                break;
        }
    }
}

void fq_nmod_mpolyun_divexact_last(fq_nmod_mpolyun_t A, fq_nmod_poly_t b,
                                                       fq_nmod_mpoly_ctx_t ctx)
{
    slong i, j;

    fq_nmod_poly_t r;
    fq_nmod_poly_init(r, ctx->fqctx);
    for (i = 0; i < A->length; i++)
    {
        for (j = 0; j < (A->coeffs + i)->length; j++)
        {
            fq_nmod_poly_divrem((A->coeffs + i)->coeffs + j, r, (A->coeffs + i)->coeffs + j, b, ctx->fqctx);
            FLINT_ASSERT(fq_nmod_poly_is_zero(r, ctx->fqctx));
        }
    }
    fq_nmod_poly_clear(r, ctx->fqctx);
}
