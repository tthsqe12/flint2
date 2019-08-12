/*
    Copyright (C) 2019 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fmpz_mpoly_factor.h"
#include "fmpq_poly.h"

void fmpz_poly_factor_print_pretty(const fmpz_poly_factor_t fac, const char * var)
{
    slong i;

    flint_printf("content: ");
    fmpz_print(&(fac->c));
    flint_printf("\n");
    for (i = 0; i < fac->num; i++)
    {
        flint_printf("factor[%wd]: ", i);
        fmpz_poly_print_pretty(fac->p + i, var);
        flint_printf(" ^ %wd\n", fac->exp[i]);
    }
}


/******* fmpz_mod_bpoly - (Z/nZ)[x,y] ****************************************/

typedef struct
{
    fmpz_mod_poly_struct * coeffs;
    slong alloc;
    slong length;
    fmpz_t modulus;
} fmpz_mod_bpoly_struct;

typedef fmpz_mod_bpoly_struct fmpz_mod_bpoly_t[1];


/******* fmpz_bpoly - Z[x,y] *************************************************/

typedef struct
{
    fmpz_poly_struct * coeffs;
    slong alloc;
    slong length;
} fmpz_bpoly_struct;

typedef fmpz_bpoly_struct fmpz_bpoly_t[1];


void fmpz_mod_bpoly_init(fmpz_mod_bpoly_t A, const fmpz_t modulus)
{
    A->coeffs = NULL;
    A->alloc = 0;
    A->length = 0;
    fmpz_init_set(A->modulus, modulus);
}

void fmpz_mod_bpoly_clear(fmpz_mod_bpoly_t A)
{
    fmpz_clear(A->modulus);
    if (A->coeffs)
    {
        slong i;
        for (i = 0; i < A->alloc; i++)
            fmpz_mod_poly_clear(A->coeffs + i);
        flint_free(A->coeffs);
    }
}

void fmpz_mod_bpoly_swap(fmpz_mod_bpoly_t A, fmpz_mod_bpoly_t B)
{
    fmpz_mod_bpoly_struct t = *A;
    *A = *B;
    *B = t;
}

void fmpz_mod_bpoly_print(fmpz_mod_bpoly_t A, const char * xvar, const char * yvar)
{
    slong i;
    int first;

    first = 1;
    for (i = A->length - 1; i >= 0; i--)
    {
        if (fmpz_mod_poly_is_zero(A->coeffs + i))
            continue;

        if (!first)
            flint_printf(" + ");

        first = 0;

        flint_printf("(");
        fmpz_mod_poly_print_pretty(A->coeffs + i, yvar);
        flint_printf(")*%s^%wd", xvar, i);
    }

    if (first)
        flint_printf("0");
}

void fmpz_mod_bpoly_fit_length(fmpz_mod_bpoly_t A, slong len)
{
    slong i;

    if (len <= A->alloc)
        return;

    if (len < 2 * A->alloc)
        len = 2 * A->alloc;

    if (A->alloc == 0)
        A->coeffs = (fmpz_mod_poly_struct *) flint_malloc(len * sizeof(fmpz_mod_poly_struct));
    else
        A->coeffs = (fmpz_mod_poly_struct *) flint_realloc(A->coeffs, len * sizeof(fmpz_mod_poly_struct));

    for (i = A->alloc; i < len; i++)
        fmpz_mod_poly_init(A->coeffs + i, A->modulus);

    A->alloc = len;
}

void fmpz_mod_bpoly_set_fmpz_bpoly(fmpz_mod_bpoly_t A, const fmpz_bpoly_t B)
{
    slong i;
    fmpz_mod_bpoly_fit_length(A, B->length);
    for (i = 0; i < B->length; i++)
    {
        fmpz_mod_poly_set_fmpz_poly(A->coeffs + i, B->coeffs + i);
        if (!fmpz_mod_poly_is_zero(A->coeffs + i))
            A->length = i + 1;
    }
}

void fmpz_mod_bpoly_set_polyx(fmpz_mod_bpoly_t A, const fmpz_mod_poly_t B)
{
    slong i;
    fmpz_mod_bpoly_fit_length(A, B->length);
    for (i = 0; i < B->length; i++)
    {
        fmpz_mod_poly_set_fmpz(A->coeffs + i, B->coeffs + i);
        if (!fmpz_mod_poly_is_zero(A->coeffs + i))
            A->length = i + 1;
    }    
}

void fmpz_mod_bpoly_set_polyy(fmpz_mod_bpoly_t A, const fmpz_mod_poly_t B)
{
    fmpz_mod_bpoly_fit_length(A, 1);
	fmpz_mod_poly_set(A->coeffs + 0, B);
	A->length = !fmpz_mod_poly_is_zero(A->coeffs + 0);
}


void fmpz_mod_bpoly_get_coeff(fmpz_t c, const fmpz_mod_bpoly_t A, slong xi, slong yi)
{
    if (xi >= A->length)
        fmpz_zero(c);

    fmpz_mod_poly_get_coeff_fmpz(c, A->coeffs + xi, yi);
}

void fmpz_mod_bpoly_make_monic(fmpz_mod_bpoly_t A, slong order)
{
    slong i;
    fmpz_mod_poly_t t, lcinv;

    FLINT_ASSERT(A->length > 0);

    fmpz_mod_poly_init(t, A->modulus);
    fmpz_mod_poly_init(lcinv, A->modulus);
    fmpz_mod_poly_inv_series(lcinv, A->coeffs + A->length - 1, order);

    for (i = 0; i < A->length; i++)
    {
        fmpz_mod_poly_mullow(t, A->coeffs + i, lcinv, order);
        fmpz_mod_poly_swap(A->coeffs + i, t);
    }

    fmpz_mod_poly_clear(t);
    fmpz_mod_poly_clear(lcinv);
}

void fmpz_mod_bpoly_mul(fmpz_mod_bpoly_t A, const fmpz_mod_bpoly_t B, const fmpz_mod_bpoly_t C, slong order)
{
    slong i, j;
    fmpz_mod_poly_t t;

    FLINT_ASSERT(A != B);
    FLINT_ASSERT(A != C);

    fmpz_mod_poly_init(t, B->modulus);

    fmpz_mod_bpoly_fit_length(A, B->length + C->length - 1);
    for (i = 0; i < B->length + C->length - 1; i++)
        fmpz_mod_poly_zero(A->coeffs + i);

    for (i = 0; i < B->length; i++)
    {
        for (j = 0; j < C->length; j++)
        {
            fmpz_mod_poly_mullow(t, B->coeffs + i, C->coeffs + j, order);
            fmpz_mod_poly_add(A->coeffs + i + j, A->coeffs + i + j, t);
        }
    }

    A->length = B->length + C->length - 1;

    fmpz_mod_poly_clear(t);
}

void fmpz_mod_bpoly_add_poly_shift(fmpz_mod_bpoly_t A, const fmpz_mod_poly_t B, slong yshift)
{
    slong i;
    fmpz_t c;

    FLINT_ASSERT(A->length > B->length);

    fmpz_init(c);

    for (i = 0; i < B->length; i++)
    {
        fmpz_mod_poly_get_coeff_fmpz(c, A->coeffs + i, yshift);
        FLINT_ASSERT(fmpz_is_zero(c));
        fmpz_mod_poly_set_coeff_fmpz(A->coeffs + i, yshift, B->coeffs + i);
    }

    fmpz_clear(c);
}

void fmpz_mod_bpoly_sub(fmpz_mod_bpoly_t A, const fmpz_mod_bpoly_t B, const fmpz_mod_bpoly_t C)
{
    slong i;

    FLINT_ASSERT(A != B);
    FLINT_ASSERT(A != C);

    fmpz_mod_bpoly_fit_length(A, FLINT_MAX(B->length, C->length));

    for (i = 0; i < FLINT_MAX(B->length, C->length) - 1; i++)
    {
        if (i < B->length)
        {
            if (i < C->length)
            {
                fmpz_mod_poly_sub(A->coeffs + i, B->coeffs + i, C->coeffs + i);
            }
            else
            {
                fmpz_mod_poly_set(A->coeffs + i, B->coeffs + i);
            }
        }
        else
        {
            FLINT_ASSERT(i < C->length);
            fmpz_mod_poly_neg(A->coeffs + i, C->coeffs + i);
        }

        if (!fmpz_mod_poly_is_zero(A->coeffs + i))
            A->length = i + 1;
    }
}

void fmpz_bpoly_init(fmpz_bpoly_t A)
{
    A->coeffs = NULL;
    A->alloc = 0;
    A->length = 0;
}

void fmpz_bpoly_clear(fmpz_bpoly_t A)
{
    if (A->coeffs)
    {
        slong i;
        for (i = 0; i < A->alloc; i++)
            fmpz_poly_clear(A->coeffs + i);
        flint_free(A->coeffs);
    }
}

void fmpz_bpoly_print(fmpz_bpoly_t A, const char * xvar, const char * yvar)
{
    slong i;
    int first;

    first = 1;
    for (i = A->length - 1; i >= 0; i--)
    {
        if (fmpz_poly_is_zero(A->coeffs + i))
            continue;

        if (!first)
            flint_printf(" + ");

        first = 0;

        flint_printf("(");
        fmpz_poly_print_pretty(A->coeffs + i, yvar);
        flint_printf(")*%s^%wd", xvar, i);
    }

    if (first)
        flint_printf("0");
}

void fmpz_bpoly_fit_length(fmpz_bpoly_t A, slong len)
{
    slong i;

    if (len <= A->alloc)
        return;

    if (len < 2 * A->alloc)
        len = 2 * A->alloc;

    if (A->alloc == 0)
        A->coeffs = (fmpz_poly_struct *) flint_malloc(len * sizeof(fmpz_poly_struct));
    else
        A->coeffs = (fmpz_poly_struct *) flint_realloc(A->coeffs, len * sizeof(fmpz_poly_struct));

    for (i = A->alloc; i < len; i++)
        fmpz_poly_init(A->coeffs + i);

    A->alloc = len;
}

void fmpz_bpoly_set(fmpz_bpoly_t A, const fmpz_bpoly_t B)
{
    slong i;

    FLINT_ASSERT(A != B);

    fmpz_bpoly_fit_length(A, B->length);
    A->length = B->length;

    for (i = 0; i < B->length; i++)
        fmpz_poly_set(A->coeffs + i, B->coeffs + i);
}

void fmpz_bpoly_make_primitive(fmpz_bpoly_t A)
{
    slong i;
    fmpz_poly_t g, q;

    fmpz_poly_init(g);
    fmpz_poly_init(q);

    for (i = 0; i < A->length; i++)
	{
        fmpz_poly_gcd(q, g, A->coeffs + i);
		fmpz_poly_swap(g, q);
	}

    for (i = 0; i < A->length; i++)
    {
        fmpz_poly_div(q, A->coeffs + i, g);
		fmpz_poly_swap(A->coeffs + i, q);
    }

    fmpz_poly_clear(g);
    fmpz_poly_clear(q);
}

void fmpz_bpoly_set_coeff(fmpz_bpoly_t A, slong xi, slong yi, const fmpz_t c)
{
    slong i;

    FLINT_ASSERT(!fmpz_is_zero(c));

    if (xi >= A->length)
    {
        fmpz_bpoly_fit_length(A, xi + 1);
        for (i = A->length; i <= xi; i++)
            fmpz_poly_zero(A->coeffs + i);
        A->length = xi + 1;
    }

    fmpz_poly_set_coeff_fmpz(A->coeffs + xi, yi, c);
}

void fmpz_bpoly_zero(fmpz_bpoly_t A)
{
    A->length = 0;
}


int fmpz_bpoly_divides(fmpz_bpoly_t Q, fmpz_bpoly_t A, fmpz_bpoly_t B)
{
    slong i, qoff;
    int divides;
    fmpz_poly_t q, t;
    fmpz_bpoly_t R;

    FLINT_ASSERT(Q != A);
    FLINT_ASSERT(Q != B);
    FLINT_ASSERT(B->length > 0);

    fmpz_poly_init(q);
    fmpz_poly_init(t);
    fmpz_bpoly_init(R);
    fmpz_bpoly_set(R, A);

    Q->length = 0;

    while (R->length >= B->length)
    {
        divides = fmpz_poly_divides(q, R->coeffs + R->length - 1, B->coeffs + B->length - 1);
        if (!divides)
            goto cleanup;

        for (i = 0; i < B->length; i++)
        {
            fmpz_poly_mul(t, B->coeffs + i, q);
            fmpz_poly_sub(R->coeffs + i + R->length - B->length, R->coeffs + i + R->length - B->length, t);
        }

        qoff = R->length - B->length;

        FLINT_ASSERT(qoff >= 0);
        if (qoff >= Q->length)
        {
            fmpz_bpoly_fit_length(Q, qoff + 1);
            for (i = Q->length; i <= qoff; i++)
                fmpz_poly_zero(Q->coeffs + i);
            Q->length = qoff + 1;
        }

        fmpz_poly_set(Q->coeffs + qoff, q);

        while (R->length > 0 && fmpz_poly_is_zero(R->coeffs + R->length - 1))
            R->length--;
    }

    divides = (R->length == 0);

cleanup:

    fmpz_poly_clear(q);
    fmpz_poly_clear(t);
    fmpz_bpoly_clear(R);

    return divides;
}


void fmpz_mpoly_from_fmpz_bpoly(
    fmpz_mpoly_t A,
    flint_bitcnt_t Abits,
    const fmpz_bpoly_t B,
    slong varx,
    slong vary,
    const fmpz_mpoly_ctx_t ctx)
{
    slong n = ctx->minfo->nvars;
    slong i, j;
    slong NA;
    slong Alen;
    fmpz * Acoeff;
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
    fmpz_mpoly_fit_bits(A, Abits, ctx);
    A->bits = Abits;

    Acoeff = A->coeffs;
    Aexp = A->exps;
    Aalloc = A->alloc;
    Alen = 0;
    for (i = 0; i < B->length; i++)
    {
        fmpz_poly_struct * Bc = B->coeffs + i;
        _fmpz_mpoly_fit_length(&Acoeff, &Aexp, &Aalloc, Alen + Bc->length, NA);

        for (j = 0; j < Bc->length; j++)
        {
            if (fmpz_is_zero(Bc->coeffs + j))
                continue;
            Aexps[varx] = i;
            Aexps[vary] = j;
            fmpz_set(Acoeff + Alen, Bc->coeffs + j);
            mpoly_set_monomial_ui(Aexp + NA*Alen, Aexps, Abits, ctx->minfo);
            Alen++;
        }
    }
    A->coeffs = Acoeff;
    A->exps = Aexp;
    A->alloc = Aalloc;
    _fmpz_mpoly_set_length(A, Alen, ctx);

    fmpz_mpoly_sort_terms(A, ctx);
    TMP_END;
}

void fmpz_bpoly_set_fmpz_mod_bpoly(fmpz_bpoly_t A, const fmpz_mod_bpoly_t B)
{
    slong i;

    fmpz_bpoly_fit_length(A, B->length);
    A->length = B->length;

    for (i = 0; i < B->length; i++)
    {
        fmpz_poly_fit_length(A->coeffs + i, (B->coeffs + i)->length);
        (A->coeffs + i)->length = (B->coeffs + i)->length;
        _fmpz_vec_scalar_smod_fmpz((A->coeffs + i)->coeffs, (B->coeffs + i)->coeffs, (B->coeffs + i)->length, B->modulus);
    }
}

void fmpz_mpoly_to_bpoly(fmpz_bpoly_t A, const fmpz_mpoly_t B, slong varx, slong vary, const fmpz_mpoly_ctx_t ctx)
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

    fmpz_bpoly_zero(A);

    for (j = 0; j < B->length; j++)
    {
        Bexpx = ((B->exps + NB*j)[Boffx] >> Bshiftx) & mask;
        Bexpy = ((B->exps + NB*j)[Boffy] >> Bshifty) & mask;
        fmpz_bpoly_set_coeff(A, Bexpx, Bexpy, B->coeffs + j);
    }
}

void fmpz_bpoly_eval(fmpz_poly_t E, const fmpz_bpoly_t A, const fmpz_t alpha)
{
    slong i;
    fmpz_t t;

    fmpz_init(t);

    fmpz_poly_zero(E);
    for (i = A->length - 1; i >= 0; i--)
    {
        fmpz_poly_evaluate_fmpz(t, A->coeffs + i, alpha);
        fmpz_poly_set_coeff_fmpz(E, i, t);
    }

    fmpz_clear(t);
}

void fmpz_bpoly_taylor_shift(fmpz_bpoly_t A, const fmpz_t alpha)
{
    slong i;
    for (i = A->length - 1; i >= 0; i--)
        _fmpz_poly_taylor_shift((A->coeffs + i)->coeffs, alpha, (A->coeffs + i)->length);    
}

typedef struct {
    slong r; /* number of local factors */
    ulong k;
    fmpz_t p;
    fmpz_t pk;
    fmpz_mod_bpoly_t Btilde;                /* mod p^k */
    fmpz_mod_bpoly_struct * newBitilde;     /* mod p^k */
    fmpz_mod_poly_struct * P;               /* mod p^k */
    fmpz_mod_poly_struct * d;               /* mod p^k */
    fmpz_mod_poly_struct * Bitilde;         /* mod p^k */
    fmpz_mod_poly_struct * d1;              /* mod p */
    fmpz_mod_poly_struct * Bitilde1;        /* mod p */
} bpoly_info_struct;

typedef bpoly_info_struct bpoly_info_t[1];

void bpoly_info_init(bpoly_info_t I, slong r, const fmpz_t p, ulong k)
{
    slong i;

    FLINT_ASSERT(r >= 2);

    I->r = r;

    I->k = k;
    fmpz_init_set(I->p, p);
    fmpz_init(I->pk);
    fmpz_pow_ui(I->pk, p, k);

    fmpz_mod_bpoly_init(I->Btilde, I->pk);

    I->newBitilde = (fmpz_mod_bpoly_struct *) flint_malloc(I->r*sizeof(fmpz_mod_bpoly_struct));
    I->P          = (fmpz_mod_poly_struct *) flint_malloc(I->r*sizeof(fmpz_mod_poly_struct));
    I->d          = (fmpz_mod_poly_struct *) flint_malloc(I->r*sizeof(fmpz_mod_poly_struct));
    I->Bitilde    = (fmpz_mod_poly_struct *) flint_malloc(I->r*sizeof(fmpz_mod_poly_struct));
    I->d1         = (fmpz_mod_poly_struct *) flint_malloc(I->r*sizeof(fmpz_mod_poly_struct));
    I->Bitilde1   = (fmpz_mod_poly_struct *) flint_malloc(I->r*sizeof(fmpz_mod_poly_struct));

    for (i = 0; i < I->r; i++)
    {
        fmpz_mod_bpoly_init(I->newBitilde + i, I->pk);
        fmpz_mod_poly_init(I->P + i, I->pk);
        fmpz_mod_poly_init(I->d + i, I->pk);
        fmpz_mod_poly_init(I->Bitilde + i, I->pk);
        fmpz_mod_poly_init(I->d1 + i, I->p);
        fmpz_mod_poly_init(I->Bitilde1 + i, I->p);
    }
}

void bpoly_info_clear(bpoly_info_t I)
{
    slong i;

    fmpz_clear(I->p);
    fmpz_clear(I->pk);

    fmpz_mod_bpoly_clear(I->Btilde);

    for (i = 0; i < I->r; i++)
    {
        fmpz_mod_bpoly_clear(I->newBitilde + i);
        fmpz_mod_poly_clear(I->P + i);
        fmpz_mod_poly_clear(I->d + i);
        fmpz_mod_poly_clear(I->Bitilde + i);
        fmpz_mod_poly_clear(I->d1 + i);
        fmpz_mod_poly_clear(I->Bitilde1 + i);
    }

    flint_free(I->newBitilde);
    flint_free(I->P);
    flint_free(I->d);
    flint_free(I->Bitilde);
    flint_free(I->d1);
    flint_free(I->Bitilde1);
}

int fmpz_mpoly_univar_content_mpoly(
    fmpz_mpoly_t g,
    const fmpz_mpoly_univar_t A,
    const fmpz_mpoly_ctx_t ctx)
{
    slong i;
    int success;

    fmpz_mpoly_zero(g, ctx);
    for (i = 0; i < A->length; i++)
    {
        success = fmpz_mpoly_gcd(g, g, A->coeffs + i, ctx);
        FLINT_ASSERT(success);
    }

    return 1;
}

void fmpz_mpoly_univar_divexact_mpoly(
    fmpz_mpoly_univar_t A,
    const fmpz_mpoly_t b,
    const fmpz_mpoly_ctx_t ctx)
{
    int success;
    slong i;

    for (i = 0; i < A->length; i++)
    {
        success = fmpz_mpoly_divides(A->coeffs + i, A->coeffs + i, b, ctx);
        FLINT_ASSERT(success);
    }
}

/*
    set out[i] so that
    1/(f[0]*f[1]*...*f[n-1]) = out[0]/f[0] + ... + out[n-1]/f[n-1]
*/
int partial_fraction_coeffs(fmpz_mod_poly_struct * out, const fmpz_mod_poly_struct * f, const fmpz_t p, slong n)
{
    slong i;
    fmpz_mod_poly_t num, den, a, b, g, t;

    FLINT_ASSERT(n > 1);

    fmpz_mod_poly_init(num, p);
    fmpz_mod_poly_init(den, p);
    fmpz_mod_poly_init(a, p);
    fmpz_mod_poly_init(b, p);
    fmpz_mod_poly_init(g, p);
    fmpz_mod_poly_init(t, p);

    fmpz_mod_poly_set_ui(num, 1);
    fmpz_mod_poly_mul(den, f + 0, f + 1);
    for (i = 2; i < n; i++)
        fmpz_mod_poly_mul(den, den, f + i);

    for (i = 0; i < n; i++)
    {
        fmpz_mod_poly_divrem(den, t, den, f + i);
        FLINT_ASSERT(fmpz_mod_poly_is_zero(t));
        fmpz_mod_poly_xgcd(g, a, b, f + i, den);
        if (fmpz_mod_poly_degree(g) != 0)
            return 0;
        FLINT_ASSERT(fmpz_is_one(g->coeffs + 0));
        fmpz_mod_poly_mul(t, b, num);
        fmpz_mod_poly_rem(out + i, t, f + i);
        fmpz_mod_poly_mul(t, a, num);
        fmpz_mod_poly_rem(num, t, den);
    }

    fmpz_mod_poly_clear(num);
    fmpz_mod_poly_clear(den);
    fmpz_mod_poly_clear(a);
    fmpz_mod_poly_clear(b);
    fmpz_mod_poly_clear(g);
    fmpz_mod_poly_clear(t);
    return 1;
}


int bpoly_info_disolve(bpoly_info_t I)
{
    slong i, j;
    fmpz_t pj, t1;
    fmpz_mod_poly_t error, t, s, s1, s2;

    if (!partial_fraction_coeffs(I->d1, I->Bitilde1, I->p, I->r))
        return 0;

    fmpz_init(pj);
    fmpz_init(t1);
    fmpz_mod_poly_init(error, I->pk);
    fmpz_mod_poly_init(t, I->pk);
    fmpz_mod_poly_init(s, I->p);
    fmpz_mod_poly_init(s1, I->p);
    fmpz_mod_poly_init(s2, I->pk);

    for (i = 0; i < I->r; i++)
    {
        fmpz_mod_poly_set_ui(I->P + i, 1);
        for (j = 0; j < I->r; j++)
        {
            if (i == j)
                continue;
            fmpz_mod_poly_mul(I->P + i, I->P + i, I->Bitilde + j);
        }
    }

    fmpz_mod_poly_set_ui(error, 1);
    for (i = 0; i < I->r; i++)
    {
        fmpz_mod_poly_set(I->d + i, I->d1 + i); /* slight abuse because moduli are different */
        fmpz_mod_poly_mul(t, I->d + i, I->P + i);
        fmpz_mod_poly_sub(error, error, t);
    }

    fmpz_one(pj);
    for (j = 1; j < I->k; j++)
    {
        fmpz_mul(pj, pj, I->p);
        fmpz_mod_poly_zero(s);
        for (i = error->length - 1; i >= 0; i--)
        {
            FLINT_ASSERT(fmpz_divisible(error->coeffs + i, pj));
            fmpz_divexact(t1, error->coeffs + i, pj);
            fmpz_mod(t1, t1, I->p);
            fmpz_mod_poly_set_coeff_fmpz(s, i, t1);
        }

        for (i = 0; i < I->r; i++)
        {
            fmpz_mod_poly_mul(s1, s, I->d1 + i);
            fmpz_set(&s2->p, I->p);
            fmpz_mod_poly_rem(s2, s1, I->Bitilde1 + i);
            fmpz_set(&s2->p, I->pk);
            fmpz_mod_poly_scalar_mul_fmpz(s2, s2, pj);
            fmpz_mod_poly_add(I->d + i, I->d + i, s2);
        }

        fmpz_mod_poly_set_ui(error, 1);
        for (i = 0; i < I->r; i++)
        {
            fmpz_mod_poly_mul(t, I->d + i, I->P + i);
            fmpz_mod_poly_sub(error, error, t);
        }
    }

    FLINT_ASSERT(fmpz_mod_poly_is_zero(error));

    fmpz_clear(pj);
    fmpz_clear(t1);
    fmpz_mod_poly_clear(error);
    fmpz_mod_poly_clear(s);
    fmpz_mod_poly_clear(s1);
    fmpz_mod_poly_clear(s2);

    return 1;
}



void fmpz_mpoly_convert_perm(
    fmpz_mpoly_t A,
    flint_bitcnt_t Abits,
    const fmpz_mpoly_ctx_t lctx,
    const fmpz_mpoly_t B,
    const fmpz_mpoly_ctx_t ctx,
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

    Aexps = (ulong *) TMP_ALLOC((m)*sizeof(ulong));
    Bexps = (ulong *) TMP_ALLOC(n*sizeof(ulong));

    NA = mpoly_words_per_exp(Abits, lctx->minfo);
    NB = mpoly_words_per_exp(B->bits, ctx->minfo);

    fmpz_mpoly_fit_bits(A, Abits, ctx);
    A->bits = Abits;
    fmpz_mpoly_fit_length(A, B->length, lctx);
    A->length = B->length;
    for (i = 0; i < B->length; i++)
    {        
	    fmpz_set(A->coeffs + i, B->coeffs + i);
	    mpoly_get_monomial_ui(Bexps, B->exps + NB*i, B->bits, ctx->minfo);
	    for (k = 0; k < m; k++)
	    {
	        l = perm[k];
	        Aexps[k] = l < 0 ? 0 : Bexps[l];
	    }
	    mpoly_set_monomial_ui(A->exps + NA*(i), Aexps, Abits, lctx->minfo);
     }  
    fmpz_mpoly_sort_terms(A, lctx);
    TMP_END;
}


void _to_poly(fmpz_poly_t A, const fmpz_mpoly_t B, const fmpz_mpoly_ctx_t ctx)
{
    slong i;
    ulong * shift, * stride;

    if (B->length == 0)
    {
        fmpz_poly_zero(A);
        return;
    }

    shift = (ulong *) flint_malloc(ctx->minfo->nvars*sizeof(ulong));
    stride = (ulong *) flint_malloc(ctx->minfo->nvars*sizeof(ulong));
    for (i = 0; i < ctx->minfo->nvars; i++)
    {
        shift[i] = 0;
        stride[i] = 1;
    }

    _fmpz_mpoly_to_fmpz_poly_deflate(A, B, 0, shift, stride, ctx);

    flint_free(shift);
    flint_free(stride);
}

void _to_polyq(fmpq_poly_t A, const fmpz_mpoly_t B, const fmpz_mpoly_ctx_t ctx)
{
    slong i;
    ulong mask;
    slong shift, off, N;
    slong Blen = B->length;
    fmpz * Bcoeff = B->coeffs;
    ulong * Bexp = B->exps;

    FLINT_ASSERT(B->bits <= FLINT_BITS);

    fmpq_poly_zero(A);

    N = mpoly_words_per_exp_sp(B->bits, ctx->minfo);
    mpoly_gen_offset_shift_sp(&off, &shift, 0, B->bits, ctx->minfo);

    mask = (-UWORD(1)) >> (FLINT_BITS - B->bits);
    for (i = 0; i < Blen; i++)
        fmpq_poly_set_coeff_fmpz(A, (Bexp[N*i + off] >> shift) & mask, Bcoeff + i);
}


void _from_poly(fmpz_mpoly_t A, flint_bitcnt_t Abits, const fmpz_poly_t B, const fmpz_mpoly_ctx_t ctx)
{
    slong i;
    ulong * shift, * stride;

    if (B->length == 0)
    {
        fmpz_mpoly_zero(A, ctx);
        return;
    }

    shift = (ulong *) flint_malloc(ctx->minfo->nvars*sizeof(ulong));
    stride = (ulong *) flint_malloc(ctx->minfo->nvars*sizeof(ulong));
    for (i = 0; i < ctx->minfo->nvars; i++)
    {
        shift[i] = 0;
        stride[i] = 1;
    }

    _fmpz_mpoly_from_fmpz_poly_inflate(A, Abits, B, 0, shift, stride, ctx);

    flint_free(shift);
    flint_free(stride);
}

int _from_polyq(fmpz_mpoly_t A, flint_bitcnt_t Abits, const fmpq_poly_t B, const fmpz_mpoly_ctx_t ctx)
{
    slong var = 0;
    slong N;
    slong k;
    slong Alen;
    fmpz * Acoeff;
    ulong * Aexp;
    slong Aalloc;
    ulong * strideexp;
    TMP_INIT;

    if (B->length == 0)
    {
        fmpz_mpoly_zero(A, ctx);
        return 1;
    }

    if (!fmpz_is_one(fmpq_poly_denref(B)))
    {
        return 0;
    }

    TMP_START;

    FLINT_ASSERT(Abits <= FLINT_BITS);

    N = mpoly_words_per_exp(Abits, ctx->minfo);
    strideexp = (ulong*) TMP_ALLOC(N*sizeof(ulong));
    mpoly_gen_monomial_sp(strideexp, var, Abits, ctx->minfo);

    fmpz_mpoly_fit_bits(A, Abits, ctx);
    A->bits = Abits;

    Acoeff = A->coeffs;
    Aexp = A->exps;
    Aalloc = A->alloc;
    Alen = 0;
    for (k = B->length - 1; k >= 0; k--)
    {
        _fmpz_mpoly_fit_length(&Acoeff, &Aexp, &Aalloc, Alen + 1, N);
        if (!fmpz_is_zero(B->coeffs + k))
        {
            fmpz_swap(Acoeff + Alen, B->coeffs + k);
            mpoly_monomial_mul_ui(Aexp + N*Alen, strideexp, N, k);
            Alen++;
        }
    }
    A->coeffs = Acoeff;
    A->exps = Aexp;
    A->alloc = Aalloc;
    _fmpz_mpoly_set_length(A, Alen, ctx);

    TMP_END;

    return 1;
}



typedef struct {
    slong n;
    slong r;
    slong l;
    fmpq_poly_struct * inv_prod_dbetas;
    fmpq_poly_struct * dbetas;
    fmpz_mpoly_struct * prod_mbetas;
    fmpz_mpoly_struct * mbetas;
    fmpz_mpoly_struct * deltas;
} mfactor_disolve_struct;

typedef mfactor_disolve_struct mfactor_disolve_t[1];


void _mfactor_disolve_init(
    mfactor_disolve_t I,
    slong l, slong r,
    const fmpz_mpoly_struct * betas,
    const fmpz * alpha,
    const fmpz_mpoly_ctx_t ctx)
{
    slong i, j, k;
    fmpz_poly_t p;
    fmpq_poly_t G, S, pq;
/*
flint_printf("_mfactor_disolve_init called(l = %wd, r = %wd)\n",l,r);
*/
    I->l = l;
    I->r = r;

    FLINT_ASSERT(l > 0);

    I->inv_prod_dbetas = (fmpq_poly_struct *) flint_malloc(l*sizeof(fmpq_poly_struct));
    I->dbetas = (fmpq_poly_struct *) flint_malloc(l*sizeof(fmpq_poly_struct));
    I->prod_mbetas = (fmpz_mpoly_struct *) flint_malloc((r + 1)*l*sizeof(fmpz_mpoly_struct));
    I->mbetas = (fmpz_mpoly_struct *) flint_malloc((r + 1)*l*sizeof(fmpz_mpoly_struct));
    I->deltas = (fmpz_mpoly_struct *) flint_malloc((r + 1)*l*sizeof(fmpz_mpoly_struct));

    fmpz_poly_init(p);
    fmpq_poly_init(G);
    fmpq_poly_init(S);
    fmpq_poly_init(pq);

    /* initialize deltas */
    for (i = r; i >= 0; i--)
        for (j = 0; j < l; j++)
            fmpz_mpoly_init(I->deltas + i*l + j, ctx);

    /* set betas */
    i = r;
    for (j = 0; j < l; j++)
    {
        fmpz_mpoly_init(I->mbetas + i*l + j, ctx);
        fmpz_mpoly_set(I->mbetas + i*l + j, betas + j, ctx);
/*
flint_printf("mbetas[%wd][%wd]: ",i,j); fmpz_mpoly_print_pretty(I->mbetas + i*l + j, NULL, ctx); printf("\n");
*/
    }
    for (i--; i >= 0; i--)
    {
        for (j = 0; j < l; j++)
        {
            fmpz_mpoly_init(I->mbetas + i*l + j, ctx);
            fmpz_mpoly_evaluate_one_fmpz(I->mbetas + i*l + j, I->mbetas + (i + 1)*l + j, i + 1, alpha + i, ctx);
/*
flint_printf("mbetas[%wd][%wd]: ",i,j); fmpz_mpoly_print_pretty(I->mbetas + i*l + j, NULL, ctx); printf("\n");
*/
        }
    }
    for (j = 0; j < l; j++)
    {
        _to_poly(p, I->mbetas + 0*l + j, ctx);
        fmpq_poly_init(I->dbetas + j);
        fmpq_poly_set_fmpz_poly(I->dbetas + j, p);
/*
flint_printf("dbetas[%wd]: ",j); fmpq_poly_print_pretty(I->dbetas + j, "x1"); printf("\n");
*/
    }

    /* set product of betas */
    for (i = r; i >= 0; i--)
    {
        for (j = 0; j < l; j++)
        {
            fmpz_mpoly_init(I->prod_mbetas + i*l + j, ctx);
            fmpz_mpoly_one(I->prod_mbetas + i*l + j, ctx);
            for (k = 0; k < l; k++)
            {
                if (k == j)
                    continue;
                fmpz_mpoly_mul(I->prod_mbetas + i*l + j, I->prod_mbetas + i*l + j, I->mbetas + i*l + k, ctx);
            }
/*
flint_printf("prod_mbetas[%wd][%wd]: ",i,j); fmpz_mpoly_print_pretty(I->prod_mbetas + i*l + j, NULL, ctx); printf("\n");
*/
        }        
    }
    for (j = 0; j < l; j++)
    {
        fmpq_poly_one(pq);
        for (k = 0; k < l; k++)
        {
            if (k == j)
                continue;
            fmpq_poly_mul(pq, pq, I->dbetas + k);
        }
        fmpq_poly_init(I->inv_prod_dbetas + j);
        fmpq_poly_xgcd(G, S, I->inv_prod_dbetas + j, I->dbetas + j, pq);
        FLINT_ASSERT(fmpq_poly_is_one(G));
/*
flint_printf("inv_prod_dbetas[%wd]: ",j); fmpq_poly_print_pretty(I->inv_prod_dbetas + j, "x1"); printf("\n");
*/
    }

    fmpz_poly_clear(p);
    fmpq_poly_clear(G);
    fmpq_poly_clear(S);
    fmpq_poly_clear(pq);
/*
printf("_mfactor_disolve_init returning\n");
*/
}


void _mfactor_disolve_clear(mfactor_disolve_t I, const fmpz_mpoly_ctx_t ctx)
{
    slong i, j;

    for (j = 0; j < I->l; j++)
    {
        fmpq_poly_clear(I->inv_prod_dbetas + j);
        fmpq_poly_clear(I->dbetas + j);
        for (i = 0; i <= I->r; i++)
        {
            fmpz_mpoly_clear(I->prod_mbetas + i*I->l + j, ctx);
            fmpz_mpoly_clear(I->mbetas + i*I->l + j, ctx);
            fmpz_mpoly_clear(I->deltas + i*I->l + j, ctx);
        }
    }

    flint_free(I->inv_prod_dbetas);
    flint_free(I->dbetas);
    flint_free(I->prod_mbetas);
    flint_free(I->mbetas);
}


int _mfactor_disolve(
    flint_bitcnt_t bits,
    slong r,
    slong num,
    const fmpz_mpoly_t t,
    const fmpz * alpha,
    const slong * deg,
    const mfactor_disolve_t I,
    const fmpz_mpoly_ctx_t ctx)
{
    slong i, j;
    int success;
    fmpz_mpoly_struct * deltas = I->deltas + r*I->l;
    fmpz_mpoly_struct * newdeltas = I->deltas + (r - 1)*I->l;

/*
flint_printf("_mfactor_disolve(r = %wd, num = %wd) called:\n", r, num);
flint_printf("t: "); fmpz_mpoly_print_pretty(t, NULL, ctx); printf("\n");
for (i = 0; i < num; i++)
{
flint_printf("betas[%wd]: ", i); fmpz_mpoly_print_pretty(betas + i, NULL, ctx); printf("\n");
}
*/

    if (r == 0)
    {
        fmpq_poly_t dtq, S, R;

        fmpq_poly_init(dtq);
        fmpq_poly_init(S);
        fmpq_poly_init(R);

        _to_polyq(dtq, t, ctx);

        success = 1;
        for (i = 0; i < num; i++)
        {
            fmpq_poly_mul(S, dtq, I->inv_prod_dbetas + i);
            fmpq_poly_rem(R, S, I->dbetas + i);
            success = _from_polyq(deltas + i, bits, R, ctx);
            if (!success)
                break;
        }

        fmpq_poly_clear(dtq);
        fmpq_poly_clear(S);
        fmpq_poly_clear(R);
    }
    else
    {
        fmpz_mpoly_t e, tt, pow, g, q, newt;

FLINT_ASSERT(num == I->l);

        fmpz_mpoly_init(e, ctx);
        fmpz_mpoly_init(tt, ctx);
        fmpz_mpoly_init(pow, ctx);
        fmpz_mpoly_init(g, ctx);
        fmpz_mpoly_init(q, ctx);
        fmpz_mpoly_init(newt, ctx);
        for (i = 0; i < num; i++)
            fmpz_mpoly_zero(deltas + i, ctx);

        fmpz_mpoly_set(e, t, ctx);

        fmpz_mpoly_one(pow, ctx);
        fmpz_mpoly_gen(g, r, ctx);
        fmpz_mpoly_sub_fmpz(g, g, alpha + r - 1, ctx);
        for (j = 0; j <= deg[r]; j++)
        {
            success = fmpz_mpoly_divides(q, e, pow, ctx);
            FLINT_ASSERT(success);
            fmpz_mpoly_evaluate_one_fmpz(newt, q, r, alpha + r - 1, ctx);
            success = _mfactor_disolve(bits, r - 1, num, newt, alpha, deg, I, ctx);
            if (!success)
                goto cleanup;

            fmpz_mpoly_set(e, t, ctx);
            for (i = 0; i < num; i++)
            {
                fmpz_mpoly_mul(tt, newdeltas + i, pow, ctx);
                fmpz_mpoly_add(deltas + i, deltas + i, tt, ctx);
                fmpz_mpoly_mul(tt, deltas + i, I->prod_mbetas + r*I->l + i, ctx);
                fmpz_mpoly_sub(e, e, tt, ctx);
            }

            fmpz_mpoly_mul(pow, pow, g, ctx);
        }

        success = fmpz_mpoly_is_zero(e, ctx);

cleanup:

        fmpz_mpoly_clear(e, ctx);
        fmpz_mpoly_clear(tt, ctx);
        fmpz_mpoly_clear(pow, ctx);
        fmpz_mpoly_clear(g, ctx);
        fmpz_mpoly_clear(q, ctx);
        fmpz_mpoly_clear(newt, ctx);
    }

/*
flint_printf("_mfactor_disolve(r = %wd, num = %wd) returning %d:\n", r, num, success);
flint_printf("t: "); fmpz_mpoly_print_pretty(t, NULL, ctx); printf("\n");
for (i = 0; i < num; i++)
{
flint_printf("deltas[%wd]: ", i); fmpz_mpoly_print_pretty(deltas + i, NULL, ctx); printf("\n");
}
*/
    return success;
}

int _mfactor_lift(
    slong m,
    fmpz_mpoly_factor_t lfac,
    const fmpz * alpha,
    const fmpz_mpoly_t A,
    const slong * degs,
    const fmpz_mpoly_ctx_t ctx)
{
    int success;
    slong i, j;
    slong r = lfac->length;
    fmpz_mpoly_t e, t, pow, g, q;
    fmpz_mpoly_struct * betas, * deltas;
    mfactor_disolve_t I;
/*
flint_printf("_mfactor_lift called (m = %wd)\n", m);
flint_printf("lfac: \n"); fmpz_mpoly_factor_print(lfac, NULL, ctx); printf("\n");
flint_printf("A: "); fmpz_mpoly_print_pretty(A, NULL, ctx); printf("\n");
*/
    FLINT_ASSERT(r > 1);

    fmpz_mpoly_init(e, ctx);
    fmpz_mpoly_init(t, ctx);
    fmpz_mpoly_init(pow, ctx);
    fmpz_mpoly_init(g, ctx);
    fmpz_mpoly_init(q, ctx);

    betas  = (fmpz_mpoly_struct * ) flint_malloc(r*sizeof(fmpz_mpoly_struct));
    for (i = 0; i < r; i++)
    {
        fmpz_mpoly_init(betas + i, ctx);
        fmpz_mpoly_evaluate_one_fmpz(betas + i, lfac->poly + i, m, alpha + m - 1, ctx);
    }

    fmpz_mpoly_mul(t, lfac->poly + 0, lfac->poly + 1, ctx);
    for (i = 2; i < r; i++)
        fmpz_mpoly_mul(t, t, lfac->poly + i, ctx);
    fmpz_mpoly_sub(e, A, t, ctx);

    fmpz_mpoly_one(pow, ctx);
    fmpz_mpoly_gen(g, m, ctx);
    fmpz_mpoly_sub_fmpz(g, g, alpha + m - 1, ctx);

    _mfactor_disolve_init(I, lfac->length, m - 1, betas, alpha, ctx);
    deltas = I->deltas + (m - 1)*I->l;

    for (j = 1; j <= degs[m]; j++)
    {
/*
flint_printf("<_mfactor_lift> j = %wd, error (length: %wd)\n", j, e->length); fmpz_mpoly_print_pretty(e, NULL, ctx); printf("\n");
*/
        if (fmpz_mpoly_is_zero(e, ctx))
        {
            success = 1;
            goto cleanup;
        }

        fmpz_mpoly_mul(pow, pow, g, ctx);
        success = fmpz_mpoly_divides(q, e, pow, ctx);
        FLINT_ASSERT(success);
        fmpz_mpoly_evaluate_one_fmpz(t, q, m, alpha + m - 1, ctx);

        success = _mfactor_disolve(A->bits, m - 1, r, t, alpha, degs, I, ctx);
        if (!success)
            goto cleanup;

        for (i = 0; i < r; i++)
        {
            fmpz_mpoly_mul(t, deltas + i, pow, ctx);
            fmpz_mpoly_add(lfac->poly + i, lfac->poly + i, t, ctx);
        }

        fmpz_mpoly_mul(t, lfac->poly + 0, lfac->poly + 1, ctx);
        for (i = 2; i < r; i++)
            fmpz_mpoly_mul(t, t, lfac->poly + i, ctx);
        fmpz_mpoly_sub(e, A, t, ctx);
    }

    success = fmpz_mpoly_is_zero(e, ctx);

cleanup:

    _mfactor_disolve_clear(I, ctx);

    fmpz_mpoly_clear(e, ctx);
    fmpz_mpoly_clear(t, ctx);
    fmpz_mpoly_clear(pow, ctx);
    fmpz_mpoly_clear(g, ctx);
    fmpz_mpoly_clear(q, ctx);

    for (i = 0; i < r; i++)
    {
        fmpz_mpoly_clear(betas + i, ctx);
    }

    flint_free(betas);

    return success;
}




int _append_irreducible_bivar_factors(
    fmpz_mpoly_factor_t fac,
    fmpz_mpoly_t A,
    slong xvar,
    slong yvar,
    slong Apow,
    const fmpz_mpoly_ctx_t ctx)
{
    slong i, j;
        fmpz_t alpha;
        fmpz_poly_t Beval;
        fmpz_bpoly_t B;
        fmpz_poly_factor_t Bevalfac;
        slong Blengthx, Blengthy;
        flint_bitcnt_t Bbits;
        ulong pkbits;
        ulong k;
        fmpz_t p;
        bpoly_info_t I;
        fmpz_mod_bpoly_t tp, tp1, error;
        fmpz_mod_poly_t ss, tt;

        k = 1;
        fmpz_init_set_ui(p, UWORD(1) << (FLINT_BITS - 1));
        fmpz_init(alpha);
        fmpz_poly_init(Beval);
        fmpz_poly_factor_init(Bevalfac);
        fmpz_bpoly_init(B);
        bpoly_info_init(I, 2, p, k);

        fmpz_mpoly_to_bpoly(B, A, xvar, yvar, ctx);
/*
printf("B: "); fmpz_bpoly_print(B, "x", "y"); printf("\n");
*/
        Blengthx = B->length;
        FLINT_ASSERT(Blengthx > 1);
/*
flint_printf("Blengthx: %wd\n");
*/
        fmpz_zero(alpha);
        goto got_alpha;

next_alpha:

        fmpz_neg(alpha, alpha);
        if (fmpz_sgn(alpha) >= 0)
            fmpz_add_ui(alpha, alpha, 1);

got_alpha:
/*
printf("<_append_irreducible_bivar_factors> trying alpha = "); fmpz_print(alpha); printf("\n");
*/
        fmpz_bpoly_eval(Beval, B, alpha);
/*
printf("Beval: "); fmpz_poly_print_pretty(Beval, "x"); printf("\n");
*/

        /* if killed leading coeff, get new alpha */
        if (Beval->length != Blengthx)
            goto next_alpha;

        fmpz_one(&Bevalfac->c);
        Bevalfac->num = 0;
        fmpz_poly_factor(Bevalfac, Beval);
/*
printf("<_append_irreducible_bivar_factors> Bevalfac:\n");
fmpz_poly_factor_print_pretty(Bevalfac, "x");
*/
        /* if multiple factors, get new alpha */
        for (i = 0; i < Bevalfac->num; i++)
        {
            if (Bevalfac->exp[i] != 1)
                goto next_alpha;
        }

		/* if one factor, A is irreducible */
		if (Bevalfac->num == 1)
		{
			fmpz_mpoly_factor_append(fac, A, Apow, ctx);
            /* TODO leak */
			return 1;
		}

/*
printf("**************\ngood alpha: "); fmpz_print(alpha); printf("\n");
*/
        fmpz_bpoly_taylor_shift(B, alpha);
/*
printf("B: "); fmpz_bpoly_print(B, "x", "y"); printf("\n");
*/

        Blengthy = 0;
        Bbits = 0;
        for (i = 0; i < B->length; i++)
        {
            slong this_bits;
            Blengthy = FLINT_MAX(Blengthy, (B->coeffs + i)->length);
            this_bits = _fmpz_vec_max_bits((B->coeffs + i)->coeffs,
                                           (B->coeffs + i)->length);
            Bbits = FLINT_MAX(Bbits, FLINT_ABS(this_bits));
        }
/*
flint_printf("Blengthy: %wd\n", Blengthy);
flint_printf("Bbits: %wu\n", Bbits);
*/
        pkbits = (FLINT_BIT_COUNT(Blengthx*Blengthy) + 1)/2;
        pkbits += Blengthx + Blengthy + Bbits - 3;
/*
flint_printf("pkbits: %wu\n", pkbits);
*/
next_prime:

        fmpz_nextprime(p, p, 1);

        FLINT_ASSERT(B->length > 0);
        FLINT_ASSERT((B->coeffs + B->length - 1)->length > 0);
        FLINT_ASSERT(!fmpz_is_zero((B->coeffs + B->length - 1)->coeffs + 0));

        if (fmpz_divisible((B->coeffs + B->length - 1)->coeffs + 0, p))
            goto next_prime;

        k = (pkbits + fmpz_bits(p))/ fmpz_bits(p);

        bpoly_info_clear(I);
        bpoly_info_init(I, Bevalfac->num, p, k);
/*
flint_printf("I->r : %wd\n", I->r);
flint_printf("I->p : "); fmpz_print(I->p); printf("\n");
flint_printf("I->pk: "); fmpz_print(I->pk); printf("\n");
*/
        fmpz_mod_bpoly_set_fmpz_bpoly(I->Btilde, B);
/*
printf("I->Btilde: "); fmpz_mod_bpoly_print(I->Btilde, "x", "y"); printf("\n");
*/
        fmpz_mod_bpoly_make_monic(I->Btilde, Blengthy);
/*
printf("I->Btilde: "); fmpz_mod_bpoly_print(I->Btilde, "x", "y"); printf("\n");
*/


        for (i = 0; i < I->r; i++)
        {
            fmpz_mod_poly_set_fmpz_poly(I->Bitilde1 + i, Bevalfac->p + i);
            fmpz_mod_poly_make_monic(I->Bitilde1 + i, I->Bitilde1 + i);
/*
flint_printf("  I->Bitilde1[%wd]: ", i); fmpz_mod_poly_print_pretty(I->Bitilde1 + i, "x"); printf("\n");
*/

            fmpz_mod_poly_set_fmpz_poly(I->Bitilde + i, Bevalfac->p + i);
            fmpz_mod_poly_make_monic(I->Bitilde + i, I->Bitilde + i);
/*
flint_printf("   I->Bitilde[%wd]: ", i); fmpz_mod_poly_print_pretty(I->Bitilde + i, "x"); printf("\n");
*/

            fmpz_mod_bpoly_set_polyx(I->newBitilde + i, I->Bitilde + i);
/*
flint_printf("I->newBitilde[%wd]: ", i); fmpz_mod_bpoly_print(I->newBitilde + i, "x", "y"); printf("\n");
*/

        }

        if (!bpoly_info_disolve(I))
            goto next_prime;

        FLINT_ASSERT(I->r > 1);

        fmpz_mod_poly_init(ss, I->pk);
        fmpz_mod_poly_init(tt, I->pk);
        fmpz_mod_bpoly_init(tp, I->pk);
        fmpz_mod_bpoly_init(tp1, I->pk);
        fmpz_mod_bpoly_init(error, I->pk);

        fmpz_mod_bpoly_mul(tp, I->newBitilde + 0, I->newBitilde + 1, Blengthy);
/*
flint_printf("       tp: "); fmpz_mod_bpoly_print(tp, "x", "y"); printf("\n");
*/
        for (i = 2; i < I->r; i++)
        {
            fmpz_mod_bpoly_mul(tp1, tp, I->newBitilde + i, Blengthy);
            fmpz_mod_bpoly_swap(tp1, tp);
/*
flint_printf("       tp: "); fmpz_mod_bpoly_print(tp, "x", "y"); printf("\n");
*/
        }
/*
flint_printf("I->Btilde: "); fmpz_mod_bpoly_print(I->Btilde, "x", "y"); printf("\n");
*/

        fmpz_mod_bpoly_sub(error, I->Btilde, tp);
/*
flint_printf("error: "); fmpz_mod_bpoly_print(error, "x", "y"); printf("\n");
*/
        for (j = 1; j < Blengthy; j++)
        {
            fmpz_mod_poly_zero(ss);
            for (i = error->length - 1; i >= 0; i--)
            {
                fmpz_t ct;
                fmpz_init(ct);

                fmpz_mod_bpoly_get_coeff(ct, error, i, j);
                fmpz_mod_poly_set_coeff_fmpz(ss, i, ct);
                for (k = 0; k < j; k++)
                {
                    fmpz_mod_bpoly_get_coeff(ct, error, i, k);
                    FLINT_ASSERT(fmpz_is_zero(ct));
                }

                fmpz_clear(ct);
            }
/*
flint_printf("j = %wd lift: ss: ", j); fmpz_mod_poly_print_pretty(ss, "x"); printf("\n");
*/
            for (i = 0; i < I->r; i++)
            {
                fmpz_mod_poly_mul(tt, ss, I->d + i);
                fmpz_mod_poly_rem(tt, tt, I->Bitilde + i);
                fmpz_mod_bpoly_add_poly_shift(I->newBitilde + i, tt, j);
            }

            fmpz_mod_bpoly_mul(tp, I->newBitilde + 0, I->newBitilde + 1, Blengthy);
            for (i = 2; i < I->r; i++)
            {
                fmpz_mod_bpoly_mul(tp1, tp, I->newBitilde + i, Blengthy);
                fmpz_mod_bpoly_swap(tp1, tp);
            }
            fmpz_mod_bpoly_sub(error, I->Btilde, tp);
/*
flint_printf("j = %wd lift: error: ", j); fmpz_mod_bpoly_print(error, "x", "y"); printf("\n");
*/
        }


        fmpz_mod_poly_clear(ss);
        fmpz_mod_poly_clear(tt);
        fmpz_mod_bpoly_clear(tp);
        fmpz_mod_bpoly_clear(tp1);
        fmpz_mod_bpoly_clear(error);

/*
printf("error: "); fmpz_mod_bpoly_print(error, "x", "y"); printf("\n");
for (i = 0; i < I->r; i++)
{
flint_printf("I->newBitilde[%wd]: ", i); fmpz_mod_bpoly_print(I->newBitilde + i, "x", "y"); printf("\n");
}
*/

        {
            fmpz_bpoly_t f, Q, R, trymez;
            fmpz_mod_bpoly_t tryme, trymet;
            fmpz_mod_poly_t leadf;
            slong k, *used_arr, *sub_arr;

            used_arr = flint_calloc(2 * I->r, sizeof(slong));
            sub_arr  = used_arr + I->r;

            fmpz_bpoly_init(f);
            fmpz_bpoly_init(Q);
            fmpz_bpoly_init(R);
            fmpz_bpoly_init(trymez);
            fmpz_mod_bpoly_init(tryme, I->pk);
            fmpz_mod_bpoly_init(trymet, I->pk);
            fmpz_mod_poly_init(leadf, I->pk);

            fmpz_bpoly_set(f, B);
            FLINT_ASSERT(f->length > 0);
            fmpz_mod_poly_set_fmpz_poly(leadf, f->coeffs + f->length - 1);
/*
printf("++++f: "); fmpz_bpoly_print(f, "x", "y"); printf("\n");
printf("leadf: "); fmpz_mod_poly_print_pretty(leadf, "y"); printf("\n");
*/
            for (k = 1; k < I->r; k++)
            {
                slong count = 0, indx = k - 1, l;

                for(l = 0; l < k; l++)
                    sub_arr[l] = l;

                sub_arr[indx]--;
                while ((indx >= 0))
                {
                    sub_arr[indx] = sub_arr[indx] + 1;

                    for (l = indx + 1; l < k; l++)
                        sub_arr[l] = sub_arr[l - 1] + 1;

                    if (sub_arr[k - 1] > I->r - 1)
                        indx--;
                    else
                    {
                        for(l = 0; l < k; l++)
                        {
                            if (used_arr[sub_arr[l]] == 1)
                                break;
                        }

                        fmpz_mod_bpoly_set_polyy(tryme, leadf);
                        for (l = 0; l < k; l++)
                        {
                            fmpz_mod_bpoly_mul(trymet, tryme, I->newBitilde + sub_arr[l], Blengthy);
                            fmpz_mod_bpoly_swap(trymet, tryme);
                        }
/*
printf("tryme: "); fmpz_mod_bpoly_print(tryme, "x", "y"); printf("\n");
*/
                        fmpz_bpoly_set_fmpz_mod_bpoly(trymez, tryme);
/*
printf("trymez: "); fmpz_bpoly_print(trymez, "x", "y"); printf("\n");
*/
                        fmpz_bpoly_make_primitive(trymez);
/*
printf("trymez: "); fmpz_bpoly_print(trymez, "x", "y"); printf("\n");
*/

                        if (fmpz_bpoly_divides(Q, f, trymez))
                        {
                            fmpz_mpoly_t goodtry;

                            fmpz_neg(alpha, alpha);
                            fmpz_bpoly_taylor_shift(trymez, alpha);
                            fmpz_neg(alpha, alpha);

                            fmpz_mpoly_init(goodtry, ctx);
/*
printf("good trymez: "); fmpz_bpoly_print(trymez, "x1", "x2"); printf("\n");
*/

                            fmpz_mpoly_from_fmpz_bpoly(goodtry, A->bits, trymez, xvar, yvar, ctx);
                            if (fmpz_sgn(goodtry->coeffs) < 0)
                                fmpz_mpoly_neg(goodtry, goodtry, ctx);
                            fmpz_mpoly_factor_append(fac, goodtry, Apow, ctx);
/*
printf("appended "); fmpz_mpoly_print_pretty(goodtry, NULL, ctx); printf("\n");
*/
                            fmpz_mpoly_clear(goodtry, ctx);

                            for(l = 0; l < k; l++)
                            {
                                used_arr[sub_arr[l]] = 1;
                                count++;
                            }

                            fmpz_bpoly_set(f, Q);
                            FLINT_ASSERT(f->length > 0);
                            fmpz_mod_poly_set_fmpz_poly(leadf, f->coeffs + f->length - 1);
/*
printf("++++f: "); fmpz_bpoly_print(f, "x", "y"); printf("\n");
printf("leadf: "); fmpz_mod_poly_print_pretty(leadf, "y"); printf("\n");
*/

                         /* If r - count = k then the rest are irreducible.  
                            TODO: Add a test for that case */
                        }

                        indx = k - 1;
                    }
                }

             /* This is where we switch to the next loop for k.  So we will have 
                found all factors using <= k local factors.  We should/could update 
                f to be the rest divided away (or multiply the remaining), could 
                also adjust r.  It is the number of remaining factors so if you 
                update then test if r = k or k+1 in which case the remaining f is 
                irreducible. */
            }

            {
                slong test = 0;

                for (k = 0; k < I->r; k++)
                    test = test + used_arr[k];

                if (test == 0)
                {
/*                fmpz_poly_factor_insert(final_fac, f, exp);*/

                    fmpz_mpoly_t goodtry;
/*
printf("----f: "); fmpz_bpoly_print(f, "x", "y"); printf("\n");
printf("leadf: "); fmpz_mod_poly_print_pretty(leadf, "y"); printf("\n");
*/

                    fmpz_neg(alpha, alpha);
                    fmpz_bpoly_taylor_shift(f, alpha);
                    fmpz_neg(alpha, alpha);

                    fmpz_mpoly_init(goodtry, ctx);

                    fmpz_mpoly_from_fmpz_bpoly(goodtry, A->bits, f, xvar, yvar, ctx);

                    if (fmpz_sgn(goodtry->coeffs) < 0)
                        fmpz_mpoly_neg(goodtry, goodtry, ctx);
                    fmpz_mpoly_factor_append(fac, goodtry, Apow, ctx);
/*
printf("appended "); fmpz_mpoly_print_pretty(goodtry, NULL, ctx); printf("\n");
*/

                    fmpz_mpoly_clear(goodtry, ctx);
                }
            }

            fmpz_bpoly_clear(f);
            fmpz_bpoly_clear(Q);
            fmpz_bpoly_clear(R);
            fmpz_bpoly_clear(trymez);
            fmpz_mod_bpoly_clear(tryme);
            fmpz_mod_bpoly_clear(trymet);
            fmpz_mod_poly_clear(leadf);
        }

        bpoly_info_clear(I);

    return 1;
}

int _append_irreducible_factors(fmpz_mpoly_factor_t fac, fmpz_mpoly_t A, slong Apow, const fmpz_mpoly_ctx_t ctx);


/* */
int fmpz_mpoly_mfactor(fmpz_mpoly_factor_t fac, fmpz_mpoly_t A, const fmpz_mpoly_ctx_t ctx)
{
    int success;
	const slong n = ctx->minfo->nvars - 1;
    slong i, j, m;
    fmpz * alpha;
	fmpz_mpoly_struct * Aevals, * newAevals, * lcAevals;
	slong * deg, * degeval, * newdeg;
	fmpz_mpoly_t t, g, lcA, newA;
	fmpz_mpoly_univar_t u;
    fmpz_mpoly_factor_t lfac;
    slong dummyvars[] = {0};
    ulong dummydegs[] = {0};

/*
flint_printf("fmpz_mpoly_mfactor(n = %wd) called\n       A: ", n); fmpz_mpoly_print_pretty(A, NULL, ctx); printf("\n");
*/
	alpha = _fmpz_vec_init(n);
    Aevals    = (fmpz_mpoly_struct *) flint_malloc(n*sizeof(fmpz_mpoly_struct));
    newAevals = (fmpz_mpoly_struct *) flint_malloc(n*sizeof(fmpz_mpoly_struct));
    lcAevals  = (fmpz_mpoly_struct *) flint_malloc(n*sizeof(fmpz_mpoly_struct));

    deg     = (slong *) flint_malloc((n + 1)*sizeof(slong));
    degeval = (slong *) flint_malloc((n + 1)*sizeof(slong));
    newdeg  = (slong *) flint_malloc((n + 1)*sizeof(slong));

	for (i = 0; i < n; i++)
	{
		fmpz_zero(alpha + i);
		fmpz_mpoly_init(Aevals + i, ctx);
		fmpz_mpoly_init(newAevals + i, ctx);
		fmpz_mpoly_init(lcAevals + i, ctx);
	}

	fmpz_mpoly_init(newA, ctx);
	fmpz_mpoly_init(lcA, ctx);
	fmpz_mpoly_init(t, ctx);
	fmpz_mpoly_init(g, ctx);

	fmpz_mpoly_univar_init(u, ctx);

    fmpz_mpoly_factor_init(lfac, ctx);


	fmpz_mpoly_degrees_si(deg, A, ctx);
	
	goto got_alpha;

next_alpha:

    for (i = 0; i < n; i++)
        fmpz_add_ui(alpha + i, alpha + i, 1);
/*
printf("!!!!!!!!!!!!!!!!!!!!!! chose alphas = "); fmpz_print(alpha + 0); printf("\n");
*/
got_alpha:

	/* ensure degrees do not drop under evalutaion */
	i = n - 1;
	fmpz_mpoly_evaluate_one_fmpz(Aevals + i, A, i + 1, alpha + i, ctx);
/*
flint_printf("Aeval[%wd]: ", i ); fmpz_mpoly_print_pretty(Aevals+ i, NULL, ctx); printf("\n");
*/
	fmpz_mpoly_degrees_si(degeval, Aevals + i, ctx);
	for (j = 0; j <= i; j++)
		if (degeval[j] != deg[j])
			goto next_alpha;
	for (i--; i >= 0; i--)
	{
		fmpz_mpoly_evaluate_one_fmpz(Aevals + i, Aevals + i + 1, i + 1, alpha + i, ctx);
/*
flint_printf("Aeval[%wd]: ", i ); fmpz_mpoly_print_pretty(Aevals + i, NULL, ctx); printf("\n");
*/
		fmpz_mpoly_degrees_si(degeval, Aevals + i, ctx);
		for (j = 0; j <= i; j++)
			if (degeval[j] != deg[j])
				goto next_alpha;
	}

	fmpz_mpoly_derivative(t, Aevals + 0, 0, ctx);
	fmpz_mpoly_gcd(t, t, Aevals + 0, ctx);
/*
flint_printf("squarefree check gcd: " ); fmpz_mpoly_print_pretty(t, NULL, ctx); printf("\n");
*/
	if (!fmpz_mpoly_is_fmpz(t, ctx))
		goto next_alpha;

	fmpz_mpoly_to_univar(u, Aevals + 1, 0, ctx);
    fmpz_mpoly_univar_content_mpoly(t, u, ctx);
/*
flint_printf("content    check gcd: " ); fmpz_mpoly_print_pretty(t, NULL, ctx); printf("\n");
*/
	if (!fmpz_mpoly_is_fmpz(t, ctx))
		goto next_alpha;

    lfac->length = 0;
    _append_irreducible_bivar_factors(lfac, Aevals + 1, 0, 1, 1, ctx);
    FLINT_ASSERT(lfac->length > 0);
/*
flint_printf("lfac:\n"); fmpz_mpoly_factor_print(lfac, NULL, ctx);
*/
    if (lfac->length == 1)
    {
/*
flint_printf("mfactor irreducible\n");
*/
	    fmpz_mpoly_factor_append(fac, A, 1, ctx);
        success = 1;
        goto cleanup;
    }

    dummyvars[0] = 0;
    dummydegs[0] = deg[0];
    fmpz_mpoly_get_coeff_vars_ui(lcA, A, dummyvars, dummydegs, 1, ctx);
/*
flint_printf("lcA: "); fmpz_mpoly_print_pretty(lcA, NULL, ctx); printf("\n");
*/

    /* rescale A */
    fmpz_mpoly_pow_ui(t, lcA, lfac->length - 1, ctx);
    fmpz_mpoly_mul(newA, A, t, ctx);
/*
flint_printf("newA: "); fmpz_mpoly_print_pretty(newA, NULL, ctx); printf("\n");
*/
    fmpz_mpoly_degrees_si(newdeg, newA, ctx);

    /* evaluations of leading coefficient */
    i = n - 1;
    fmpz_mpoly_evaluate_one_fmpz(lcAevals + i, lcA, i + 1, alpha + i, ctx);
/*
flint_printf("lcAevals[%wd]: ", i); fmpz_mpoly_print_pretty(lcAevals + i, NULL, ctx); printf("\n");
*/
    for (i--; i >= 0; i--)
    {
	    fmpz_mpoly_evaluate_one_fmpz(lcAevals + i, lcAevals + i + 1, i + 1, alpha + i, ctx);
/*
flint_printf("lcAevals[%wd]: ", i); fmpz_mpoly_print_pretty(lcAevals + i, NULL, ctx); printf("\n");
*/
    }

    /* evaluations */
    i = n - 1;
	fmpz_mpoly_evaluate_one_fmpz(newAevals + i, newA, i + 1, alpha + i, ctx);
/*
flint_printf("newAevals[%wd]: ", i ); fmpz_mpoly_print_pretty(newAevals + i, NULL, ctx); printf("\n");
*/
	for (i--; i >= 0; i--)
	{

		fmpz_mpoly_evaluate_one_fmpz(newAevals + i, newAevals + i + 1, i + 1, alpha + i, ctx);
/*
flint_printf("newAevals[%wd]: ", i ); fmpz_mpoly_print_pretty(newAevals + i, NULL, ctx); printf("\n");
*/
	}

    /* rescale factors */
    for (i = 0; i < lfac->length; i++)
    {
        dummyvars[0] = 0;
        dummydegs[0] = fmpz_mpoly_degree_si(lfac->poly + i, 0, ctx);
        fmpz_mpoly_get_coeff_vars_ui(g, lfac->poly + i, dummyvars, dummydegs, 1, ctx);
/*
flint_printf("lfac[%wd]lc: ", i); fmpz_mpoly_print_pretty(g, NULL, ctx); printf("\n");
*/
        success = fmpz_mpoly_divides(t, lcAevals + 1, g, ctx);
        FLINT_ASSERT(success);
/*
flint_printf("lfac[%wd] t: ", i); fmpz_mpoly_print_pretty(g, NULL, ctx); printf("\n");
*/
        fmpz_mpoly_mul(lfac->poly + i, lfac->poly + i, t, ctx);
/*
flint_printf("lfac[%wd]: ", i); fmpz_mpoly_print_pretty(lfac->poly + i, NULL, ctx); printf("\n");
*/
    }
/*
flint_printf("lift m = %wd, lfac:\n", 1); fmpz_mpoly_factor_print(lfac, NULL, ctx);
*/

    for (m = 2; m <= n; m++)
    {
        for (i = 0; i < lfac->length; i++)
        {
            slong dd = fmpz_mpoly_degree_si(lfac->poly + i, 0, ctx);
            fmpz_mpoly_gen(g, 0, ctx);
            fmpz_mpoly_pow_ui(g, g, dd, ctx);
            fmpz_mpoly_sub(t, m < n ? lcAevals + m : lcA, lcAevals + m - 1, ctx);
            fmpz_mpoly_mul(t, t, g, ctx);
            fmpz_mpoly_add(lfac->poly + i, lfac->poly + i, t, ctx);
        }

        success = _mfactor_lift(m, lfac, alpha, m < n ? newAevals + m : newA, newdeg, ctx);
        if (!success)
            goto next_alpha;
/*
flint_printf("lift m = %wd, lfac:\n", m); fmpz_mpoly_factor_print(lfac, NULL, ctx);
*/
    }

    for (i = 0; i < lfac->length; i++)
    {
        fmpz_mpoly_to_univar(u, lfac->poly + i, 0, ctx);
        success = fmpz_mpoly_univar_content_mpoly(g, u, ctx);
        if (!success)
            goto cleanup;
        success = fmpz_mpoly_divides(t, lfac->poly + i, g, ctx);
        FLINT_ASSERT(success);
        fmpz_mpoly_factor_append(fac, t, 1, ctx);
    }

    success = 1;

cleanup:

    fmpz_mpoly_factor_clear(lfac, ctx);
	fmpz_mpoly_univar_clear(u, ctx);
	fmpz_mpoly_clear(newA, ctx);
	fmpz_mpoly_clear(lcA, ctx);
	fmpz_mpoly_clear(t, ctx);
	fmpz_mpoly_clear(g, ctx);
	for (i = 0; i < n; i++)
	{
		fmpz_mpoly_clear(Aevals + i, ctx);
		fmpz_mpoly_clear(newAevals + i, ctx);
		fmpz_mpoly_clear(lcAevals + i, ctx);
	}
    flint_free(Aevals);
    flint_free(newAevals);
    flint_free(lcAevals);
    flint_free(deg);
    flint_free(degeval);
    flint_free(newdeg);
	_fmpz_vec_clear(alpha, n);

	return success;
}



/*
    A is primitive w.r.t to any variable appearing in A.
    A is square free.

    return 1 for success, 0 for failure
*/
int _append_irreducible_factors(fmpz_mpoly_factor_t fac, fmpz_mpoly_t A, slong Apow, const fmpz_mpoly_ctx_t ctx)
{
    slong i, mvars;
    slong * Adegs, * perm, * iperm;
    ulong * shift, * stride;
    TMP_INIT;
/*
printf("************_append_irreducible_factors:\n");
fmpz_mpoly_print_pretty(A, NULL, ctx); printf("\n");
*/
    if (A->bits > FLINT_BITS)
        return 0;

    TMP_START;

    Adegs = (slong *) TMP_ALLOC(ctx->minfo->nvars*sizeof(slong));
    perm = (slong *) TMP_ALLOC(ctx->minfo->nvars*sizeof(slong));
    iperm = (slong *) TMP_ALLOC(ctx->minfo->nvars*sizeof(slong));
	shift = (ulong *) TMP_ALLOC(ctx->minfo->nvars*sizeof(slong));
    stride = (ulong *) TMP_ALLOC(ctx->minfo->nvars*sizeof(slong));
    fmpz_mpoly_degrees_si(Adegs, A, ctx);

    mvars = 0;

    for (i = 0; i < ctx->minfo->nvars; i++)
    {
        shift[i] = 0;
        stride[i] = 1;
		iperm[i] = -1;
        if (Adegs[i] > 0)
        {
            perm[mvars] = i;
            mvars++;
        }
    }
	/* TODO: figure out nice perm */

	for (i = 0; i < mvars; i++)
		iperm[perm[i]] = i;

    if (mvars == 1)
    {
        fmpz_mpoly_t t;
        fmpz_poly_t Au;
        fmpz_poly_factor_t facu;

        fmpz_mpoly_init(t, ctx);
        fmpz_poly_init(Au);
        fmpz_poly_factor_init(facu);

        _fmpz_mpoly_to_fmpz_poly_deflate(Au, A, perm[0], shift, stride, ctx);
        fmpz_poly_factor(facu, Au);

        fmpz_mul(fac->content, fac->content, &facu->c); /* facu->c should be 1 */
        for (i = 0; i < facu->num; i++)
        {
            _fmpz_mpoly_from_fmpz_poly_inflate(t, A->bits, facu->p + i, perm[0], shift, stride, ctx);
            fmpz_mpoly_factor_append(fac, t, facu->exp[i]*Apow, ctx); /* facu->exp[i] should be 1 */
        }

        fmpz_mpoly_clear(t, ctx);
        fmpz_poly_clear(Au);
        fmpz_poly_factor_clear(facu);
        goto cleanup;
    }
    else if (mvars == 2)
    {
        _append_irreducible_bivar_factors(fac, A, perm[0], perm[1], Apow, ctx);
    }
    else
    {
		fmpz_mpoly_ctx_t lctx;
		fmpz_mpoly_t Al, B;
		fmpz_mpoly_factor_t Amfactors;
		int success;

		fmpz_mpoly_ctx_init(lctx, mvars, ORD_LEX);
		fmpz_mpoly_init(Al, lctx);
        fmpz_mpoly_init(B, ctx);
		fmpz_mpoly_convert_perm(Al, A->bits, lctx, A, ctx, perm);
		fmpz_mpoly_factor_init(Amfactors, lctx);
		success = fmpz_mpoly_mfactor(Amfactors, Al, lctx);
		if (!success)
			return 0;
		for (i = 0; i < Amfactors->length; i++)
		{
			fmpz_mpoly_convert_perm(B, A->bits, ctx, Amfactors->poly + i, lctx, iperm);
            FLINT_ASSERT(B->length > 0);
            if (fmpz_sgn(B->coeffs + 0) < 0)
                fmpz_mpoly_neg(B, B, ctx);
			fmpz_mpoly_factor_append(fac, B, Apow, ctx);
		}
		fmpz_mpoly_factor_clear(Amfactors, lctx);
		fmpz_mpoly_clear(B, ctx);
		fmpz_mpoly_clear(Al, lctx);
		fmpz_mpoly_ctx_clear(lctx);
		return 1;
    }

cleanup:

    TMP_END;

    return 1;
}



void fmpz_mpoly_factor_swap(fmpz_mpoly_factor_t A, fmpz_mpoly_factor_t B, const fmpz_mpoly_ctx_t ctx)
{
   fmpz_mpoly_factor_struct t = *A;
   *A = *B;
   *B = t;
}




void fmpz_mpoly_univar_shift_right(
    fmpz_mpoly_univar_t A,
    ulong shift,
    const fmpz_mpoly_ctx_t ctx)
{
    slong i;
    for (i = 0; i < A->length; i++)
    {
        FLINT_ASSERT(A->exps[i] >= shift);
        A->exps[i] -= shift;
    }
}

void fmpz_mpoly_factor_expand(fmpz_mpoly_t a, const fmpz_mpoly_factor_t f, const fmpz_mpoly_ctx_t ctx)
{
    slong i;
    fmpz_mpoly_t t1, t2;

    fmpz_mpoly_init(t1, ctx);
    fmpz_mpoly_init(t2, ctx);

    fmpz_mpoly_set_fmpz(a, f->content, ctx);

    for (i = 0; i < f->length; i++)
    {
        fmpz_mpoly_pow_ui(t1, f->poly + i, f->exp[i], ctx);
        fmpz_mpoly_mul(t2, a, t1, ctx);
        fmpz_mpoly_swap(a, t2, ctx);
    }

    fmpz_mpoly_clear(t1, ctx);
    fmpz_mpoly_clear(t2, ctx);
}

int fmpz_mpoly_factor_matches(const fmpz_mpoly_t a, const fmpz_mpoly_factor_t f, const fmpz_mpoly_ctx_t ctx)
{
    int matches;
    fmpz_mpoly_t t;
    fmpz_mpoly_init(t, ctx);
    fmpz_mpoly_factor_expand(t, f, ctx);
    matches = fmpz_mpoly_equal(t, a, ctx);
    fmpz_mpoly_clear(t, ctx);
    return matches;
}

/*
    S is nonzero and primitve w.r.t variable var
    Sp is the derivative of S w.r.t variable var

    return 1 for success, 0 for failure
*/
int _append_squarefree_factors(fmpz_mpoly_factor_t f, fmpz_mpoly_t S, fmpz_mpoly_t Sp, slong var, slong pow, const fmpz_mpoly_ctx_t ctx)
{
    slong k;
    int success;
    fmpz_mpoly_t Sm, Ss, Y, Z;

    fmpz_mpoly_init(Sm, ctx);
    fmpz_mpoly_init(Ss, ctx);
    fmpz_mpoly_init(Y, ctx);
    fmpz_mpoly_init(Z, ctx);

    success = fmpz_mpoly_gcd(Sm, S, Sp, ctx);
    if (!success)
        goto cleanup;
    fmpz_mpoly_divides(Ss, S, Sm, ctx);
    fmpz_mpoly_divides(Y, Sp, Sm, ctx);
    FLINT_ASSERT(!fmpz_mpoly_is_zero(Ss, ctx));
    FLINT_ASSERT(!fmpz_mpoly_is_zero(Y, ctx));

    k = 1;

    while (fmpz_mpoly_derivative(Sp, Ss, var, ctx),
           fmpz_mpoly_sub(Z, Y, Sp, ctx),
           !fmpz_mpoly_is_zero(Z, ctx))
    {
        fmpz_mpoly_gcd(S, Ss, Z, ctx);
        if (!success)
            goto cleanup;
        fmpz_mpoly_divides(Sp, Ss, S, ctx);
        fmpz_mpoly_swap(Ss, Sp, ctx);
        fmpz_mpoly_divides(Y, Z, S, ctx);
        FLINT_ASSERT(!fmpz_mpoly_is_zero(Ss, ctx));
        FLINT_ASSERT(!fmpz_mpoly_is_zero(Y, ctx));

        if (!fmpz_mpoly_is_fmpz(S, ctx))
        {
            fmpz_mpoly_factor_append(f, S, k*pow, ctx);
        }
        else
        {
            FLINT_ASSERT(fmpz_mpoly_is_one(S, ctx));
        }

        k++;
    }

    if (!fmpz_mpoly_is_fmpz(Ss, ctx))
    {
        fmpz_mpoly_factor_append(f, Ss, k*pow, ctx);
    }
    else
    {
        FLINT_ASSERT(fmpz_mpoly_is_one(Ss, ctx));
    }

cleanup:

    fmpz_mpoly_clear(Sm, ctx);
    fmpz_mpoly_clear(Ss, ctx);
    fmpz_mpoly_clear(Y, ctx);
    fmpz_mpoly_clear(Z, ctx);

    return success;
}

int fmpz_mpoly_factor(fmpz_mpoly_factor_t f, const fmpz_mpoly_t A, int full, const fmpz_mpoly_ctx_t ctx)
{
    int success;
    slong j, v;
    fmpz_mpoly_factor_t newf;
    fmpz_t cc;
    fmpz_mpoly_t c;
    fmpz_mpoly_univar_t u;

    /* 0. set trivial factorization */
    f->length = 0;
    if (fmpz_mpoly_is_fmpz(A, ctx))
    {
        fmpz_mpoly_get_fmpz(f->content, A, ctx);
        return 1;
    }
    else
    {
        fmpz_one(f->content);
        fmpz_mpoly_factor_append(f, A, 1, ctx);

        if (fmpz_sgn(f->poly->coeffs + 0) < 0)
        {
            fmpz_neg(f->content, f->content);
            fmpz_mpoly_neg(f->poly + 0, f->poly + 0, ctx);
        }
    }

    if (A->bits > FLINT_BITS)
    {
        return 0;
    }

    fmpz_mpoly_factor_init(newf, ctx);
    fmpz_mpoly_univar_init(u, ctx);
    fmpz_mpoly_init(c, ctx);
    fmpz_init(cc);

    /* 1. ensure factors are primitive w.r.t any variable */
    FLINT_ASSERT(fmpz_mpoly_factor_matches(A, f, ctx));

    for (v = 0; v < ctx->minfo->nvars; v++)
    {
        fmpz_set(newf->content, f->content);
        newf->length = 0;

        for (j = 0; j < f->length; j++)
        {
            if (f->exp[j] == 1)
            {
                fmpz_mpoly_to_univar(u, f->poly + j, v, ctx);
                FLINT_ASSERT(u->length > 0);

                fmpz_mpoly_univar_content_mpoly(c, u, ctx);
                fmpz_mpoly_univar_divexact_mpoly(u, c, ctx);

                if (fmpz_mpoly_is_fmpz(c, ctx))
                {
                    FLINT_ASSERT(c->length == 1);
                    fmpz_mul(newf->content, newf->content, c->coeffs + 0);
                }
                else
                {
                    fmpz_mpoly_factor_append(newf, c, 1, ctx);
                }

                FLINT_ASSERT(u->length > 0);

                if (u->exps[u->length - 1] != 0)
                {
                    fmpz_mpoly_gen(c, v, ctx);
                    fmpz_mpoly_factor_append(newf, c, u->exps[u->length - 1], ctx);
                    fmpz_mpoly_univar_shift_right(u, u->exps[u->length - 1], ctx);
                }

                if (u->length > 1)
                {
                    fmpz_mpoly_from_univar_bits(c, A->bits, u, ctx);
                    fmpz_mpoly_factor_append(newf, c, 1, ctx);
                }
                else
                {
                    FLINT_ASSERT(fmpz_mpoly_is_one(u->coeffs + 0, ctx));
                }
            }
            else
            {
                FLINT_ASSERT(f->exp[j] > 1);
                fmpz_mpoly_factor_append(newf, f->poly + j, f->exp[j], ctx);
            }
        }

        fmpz_mpoly_factor_swap(f, newf, ctx);
    }

    /* 2. ensure factors are squarefree */
    FLINT_ASSERT(fmpz_mpoly_factor_matches(A, f, ctx));

    fmpz_set(newf->content, f->content);
    newf->length = 0;

    for (j = 0; j < f->length; j++)
    {
        for (v = 0; v < ctx->minfo->nvars; v++)
        {
            fmpz_mpoly_derivative(c, f->poly + j, v, ctx);
            if (!fmpz_mpoly_is_zero(c, ctx))
                break;
        }

        if (v < ctx->minfo->nvars)
        {
            success = _append_squarefree_factors(newf, f->poly + j, c, v, f->exp[j], ctx);
            if (!success)
                goto cleanup;
        }
        else
        {
            /*
                The factors in f should all be nonconstant,
                but we can correct a constant factor here.
            */
            FLINT_ASSERT(fmpz_mpoly_is_fmpz(f->poly + j, ctx));
            fmpz_mpoly_get_fmpz(cc, f->poly + j, ctx);
            fmpz_mul(newf->content, newf->content, cc);
        }
    }

    fmpz_mpoly_factor_swap(f, newf, ctx);

    /* 3. ensure factors are irreducible */
    if (!full)
        goto cleanup;

    FLINT_ASSERT(fmpz_mpoly_factor_matches(A, f, ctx));

    fmpz_set(newf->content, f->content);
    newf->length = 0;

    for (j = 0; j < f->length; j++)
    {
        _append_irreducible_factors(newf, f->poly + j, f->exp[j], ctx);
    }

    fmpz_mpoly_factor_swap(f, newf, ctx);

    success = 1;

cleanup:

    FLINT_ASSERT(fmpz_mpoly_factor_matches(A, f, ctx));

    fmpz_mpoly_factor_clear(newf, ctx);
    fmpz_mpoly_univar_clear(u, ctx);
    fmpz_mpoly_clear(c, ctx);
    fmpz_clear(cc);

    return 1;
}

