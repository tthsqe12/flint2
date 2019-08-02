/*
    Copyright (C) 2019 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fmpz_mpoly_factor.h"


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

void fmpz_mod_bpoly_print(fmpz_mod_bpoly_t A)
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
        fmpz_mod_poly_print_pretty(A->coeffs + i, "y");
        flint_printf(")*x^%wd", i);
    }

    if (first)
        flint_printf("0");
}



/******* fmpz_bpoly - Z[x,y] *************************************************/

typedef struct
{
    fmpz_poly_struct * coeffs;
    slong alloc;
    slong length;
} fmpz_bpoly_struct;

typedef fmpz_bpoly_struct fmpz_bpoly_t[1];

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

/*
    A is primitive w.r.t to any variable appearing in A.
    A is square free.

    return 1 for success, 0 for failure
*/
int _append_irreducible_factors(fmpz_mpoly_factor_t f, fmpz_mpoly_t A, slong Apow, const fmpz_mpoly_ctx_t ctx)
{
    slong i, mvars;
    slong * Adegs, * perm;
    ulong * shift, * stride;
    TMP_INIT;

    if (A->bits > FLINT_BITS)
        return 0;

    TMP_START;

    Adegs = (slong *) TMP_ALLOC(ctx->minfo->nvars*sizeof(slong));
    perm = (slong *) TMP_ALLOC(ctx->minfo->nvars*sizeof(slong));
    shift = (ulong *) TMP_ALLOC(ctx->minfo->nvars*sizeof(slong));
    stride = (ulong *) TMP_ALLOC(ctx->minfo->nvars*sizeof(slong));
    fmpz_mpoly_degrees_si(Adegs, A, ctx);

    mvars = 0;

    for (i = 0; i < ctx->minfo->nvars; i++)
    {
        shift[i] = 0;
        stride[i] = 1;
        if (Adegs[i] > 0)
        {
            perm[mvars] = i;
            mvars++;
        }
    }

    if (mvars == 1)
    {
        fmpz_mpoly_t t;
        fmpz_poly_t Au;
        fmpz_poly_factor_t fu;

        fmpz_mpoly_init(t, ctx);
        fmpz_poly_init(Au);
        fmpz_poly_factor_init(fu);

        _fmpz_mpoly_to_fmpz_poly_deflate(Au, A, perm[0], shift, stride, ctx);
        fmpz_poly_factor(fu, Au);

        fmpz_mul(f->content, f->content, &fu->c); /* fu->c should be 1 */
        for (i = 0; i < fu->num; i++)
        {
            _fmpz_mpoly_from_fmpz_poly_inflate(t, A->bits, fu->p + i, perm[0], shift, stride, ctx);
            fmpz_mpoly_factor_append(f, t, fu->exp[i]*Apow, ctx); /* fu->exp[i] should be 1 */
        }

        fmpz_mpoly_clear(t, ctx);
        fmpz_poly_clear(Au);
        fmpz_poly_factor_clear(fu);
        goto cleanup;
    }
    else if (mvars == 2)
    {
        fmpz_t alpha;
        fmpz_poly_t Beval;
        fmpz_bpoly_t B;
        fmpz_poly_factor_t Bevalfac;
        slong Blengthx, Blengthy;
        flint_bitcnt_t Bbits;


        fmpz_poly_init(Beval);
        fmpz_poly_factor_init(Bevalfac);
        fmpz_bpoly_init(B);



        fmpz_mpoly_to_bpoly(B, A, perm[0], perm[1], ctx);
printf("B: "); fmpz_bpoly_print(B, "x", "y"); printf("\n");

        Blengthx = B->length;
        FLINT_ASSERT(Blengthx > 1);
flint_printf("Blengthx: %wd\n");


        fmpz_zero(alpha);
        goto got_alpha;

next_alpha:

        fmpz_neg(alpha, alpha);
        if (fmpz_sgn(alpha) >= 0)
            fmpz_add_ui(alpha, alpha, 1);

got_alpha:

printf("trying alpha = "); fmpz_print(alpha); printf("\n");
usleep(1000000);

        fmpz_bpoly_eval(Beval, B, alpha);

printf("Beval: "); fmpz_poly_print_pretty(Beval, "x"); printf("\n");


        /* if killed leading coeff, get new alpha */
        if (Beval->length != Blengthx)
            goto next_alpha;

        fmpz_poly_factor(Bevalfac, Beval);

printf("Bevalfac:\n");
fmpz_poly_factor_print_pretty(Bevalfac, "x");

        /* if multiple factors, get new alpha */
        for (i = 0; i < Bevalfac->num; i++)
        {
            if (Bevalfac->exp[i] != 1)
                goto next_alpha;
        }

printf("good alpha: "); fmpz_print(alpha); printf("\n");

        fmpz_bpoly_taylor_shift(B, alpha);

printf("B: "); fmpz_bpoly_print(B, "x", "y"); printf("\n");


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
flint_printf("Blengthy: %wd\n", Blengthy);
flint_printf("Bbits: %wu\n", Bbits);


        
    }

    fmpz_mpoly_factor_append(f, A, Apow, ctx); /* TODO */


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
int _append_squarefree_factors(fmpz_mpoly_factor_t f, fmpz_mpoly_t S, fmpz_mpoly_t Sp, slong var, const fmpz_mpoly_ctx_t ctx)
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
            fmpz_mpoly_factor_append(f, S, k, ctx);
        }
        else
        {
            FLINT_ASSERT(fmpz_mpoly_is_one(S, ctx));
        }

        k++;
    }

    if (!fmpz_mpoly_is_fmpz(Ss, ctx))
    {
        fmpz_mpoly_factor_append(f, Ss, k, ctx);
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
            success = _append_squarefree_factors(newf, f->poly + j, c, v, ctx);
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

