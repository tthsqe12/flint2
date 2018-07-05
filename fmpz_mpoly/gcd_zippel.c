/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "nmod_mpoly.h"
#include "fmpz_mpoly.h"

void nmod_mpolyu_setform(nmod_mpolyu_t A, nmod_mpolyu_t B, nmod_mpoly_ctx_t ctx);

void nmod_mpolyun_init(nmod_mpolyun_t A, mp_bitcnt_t bits, const nmod_mpoly_ctx_t ctx);

void nmod_mpolyu_cvtto_mpolyun(
    nmod_mpolyun_t A,
    nmod_mpolyu_t B,
    slong k,
    nmod_mpoly_ctx_t ctx);

int nmod_mpolyu_pgcd_zippel_univar(
    nmod_mpolyu_t G,
    nmod_mpolyu_t A,
    nmod_mpolyu_t B,
    nmod_mpoly_ctx_t ctx);

int nmod_mpolyu_pgcd_zippel_bivar(
    nmod_mpolyu_t G,
    nmod_mpolyu_t A,
    nmod_mpolyu_t B,
    nmod_mpoly_ctx_t ctx,
    mpoly_zipinfo_t zinfo);

void nmod_mpolyun_print_pretty(const nmod_mpolyun_t poly,
                                   const char ** x, const nmod_mpoly_ctx_t ctx);
void nmod_mpolyu_print_pretty(const nmod_mpolyu_t poly,
                                   const char ** x, const nmod_mpoly_ctx_t ctx);

void nmod_mpolyun_shift_right(nmod_mpolyun_t A, ulong s);
void nmod_mpolyun_shift_left(nmod_mpolyun_t A, ulong s);

void nmod_mpolyun_content_last(nmod_poly_t a, nmod_mpolyun_t B, nmod_mpoly_ctx_t ctx);

nmod_poly_struct * nmod_mpolyn_leadcoeff_ref(nmod_mpolyn_t A, nmod_mpoly_ctx_t ctx);
nmod_poly_struct * nmod_mpolyun_leadcoeff_ref(nmod_mpolyun_t A, nmod_mpoly_ctx_t ctx);
mp_limb_t nmod_mpoly_leadcoeff(nmod_mpoly_t A, nmod_mpoly_ctx_t ctx);
mp_limb_t nmod_mpolyu_leadcoeff(nmod_mpolyu_t A, nmod_mpoly_ctx_t ctx);

void nmod_mpolyu_init(nmod_mpolyu_t A, mp_bitcnt_t bits, const nmod_mpoly_ctx_t ctx);

void nmod_mpolyun_eval_last(nmod_mpolyu_t B, nmod_mpolyun_t A, mp_limb_t alpha, nmod_mpoly_ctx_t ctx);

void nmod_mpolyun_zero(nmod_mpolyun_t A, const nmod_mpoly_ctx_t ctx);

int nmod_mpolyu_is_one(nmod_mpolyu_t A, const nmod_mpoly_ctx_t uctx);

void nmod_mpolyu_cvtfrom_poly_notmain(nmod_mpolyu_t A, nmod_poly_t a, slong var, nmod_mpoly_ctx_t ctx);

void nmod_mpolyu_scalar_mul_nmod(nmod_mpolyu_t A, mp_limb_t c, nmod_mpoly_ctx_t ctx);

int nmod_mpolyun_addinterp(
    slong * lastdeg,
    nmod_mpolyun_t F,
    nmod_mpolyun_t T,
    nmod_mpolyu_t A,
    nmod_poly_t modulus,
    mp_limb_t alpha,
    const nmod_mpoly_ctx_t ctx);

void nmod_mpolyun_mul_poly(
    nmod_mpolyun_t A,
    const nmod_mpolyun_t B,
    const nmod_poly_t c,
    const nmod_mpoly_ctx_t ctx);


int nmod_mpolyu_gcd_zippel_linzipp(
    nmod_mpolyu_t G,
    nmod_mpolyu_t A,
    nmod_mpolyu_t B,
    slong var,
    nmod_mpoly_ctx_t ctx,
    mpoly_zipinfo_t zinfo,
    flint_rand_t randstate);


int nmod_mpolyu_sgcd_zippel(
    nmod_mpolyu_t G,
    nmod_mpolyu_t A,
    nmod_mpolyu_t B,
    nmod_mpolyu_t f,
    slong var,
    nmod_mpoly_ctx_t ctx,
    flint_rand_t randstate);


void nmod_mpolyu_fit_length(nmod_mpolyu_t A, slong length, const nmod_mpoly_ctx_t uctx);

void nmod_mpolyu_clear(nmod_mpolyu_t A, const nmod_mpoly_ctx_t uctx);

/*
    fmpz_mpolyu_t
    sparse univariates with fmpz_mpoly_t coefficients
*/
typedef struct
{
   fmpz_mpoly_struct * coeffs;
   ulong * exps;
   slong alloc;
   slong length;
   mp_bitcnt_t bits;    /* default bits to construct coeffs */
} fmpz_mpolyu_struct;
typedef fmpz_mpolyu_struct fmpz_mpolyu_t[1];





void fmpz_mpolyu_init(fmpz_mpolyu_t A, mp_bitcnt_t bits, const fmpz_mpoly_ctx_t ctx)
{
    A->coeffs = NULL;
    A->exps = NULL;
    A->alloc = 0;
    A->length = 0;
    A->bits = bits;
}


void fmpz_mpolyu_clear(fmpz_mpolyu_t A, const fmpz_mpoly_ctx_t uctx)
{
    slong i;
    for (i = 0; i < A->alloc; i++)
        fmpz_mpoly_clear(A->coeffs + i, uctx);
    flint_free(A->coeffs);
    flint_free(A->exps);
}


void fmpz_mpolyu_swap(fmpz_mpolyu_t A, fmpz_mpolyu_t B, const fmpz_mpoly_ctx_t uctx)
{
   fmpz_mpolyu_struct t = *A;
   *A = *B;
   *B = t;
}

void fmpz_mpolyu_zero(fmpz_mpolyu_t A, const fmpz_mpoly_ctx_t uctx)
{
    slong i;
    for (i = 0; i < A->alloc; i++) {
        fmpz_mpoly_clear(A->coeffs + i, uctx);
        fmpz_mpoly_init(A->coeffs + i, uctx);
    }
    A->length = 0;
}

/*
int nmod_mpolyu_is_one(nmod_mpolyu_t A, const nmod_mpoly_ctx_t uctx)
{
    if (A->length != 1 || A->exps[0] != UWORD(0))
        return 0;

    return nmod_mpoly_is_one(A->coeffs + 0, uctx);
}
*/



void fmpz_mpolyu_print_pretty(const fmpz_mpolyu_t poly,
                                   const char ** x, const fmpz_mpoly_ctx_t ctx)
{
    slong i;
    if (poly->length == 0)
        flint_printf("0");
    for (i = 0; i < poly->length; i++)
    {
        if (i != 0)
            flint_printf(" + ");
        flint_printf("(");
        fmpz_mpoly_print_pretty(poly->coeffs + i, x, ctx);
        flint_printf(")*X^%wd", poly->exps[i]);
    }
}



void fmpz_mpolyu_fit_length(fmpz_mpolyu_t A, slong length, const fmpz_mpoly_ctx_t uctx)
{
    slong i;
    slong old_alloc = A->alloc;
    slong new_alloc = FLINT_MAX(length, 2*A->alloc);

    if (length > old_alloc)
    {
        if (old_alloc == 0)
        {
            A->exps = (ulong *) flint_malloc(new_alloc*sizeof(ulong));
            A->coeffs = (fmpz_mpoly_struct *) flint_malloc(
                                          new_alloc*sizeof(fmpz_mpoly_struct));
        } else
        {
            A->exps = (ulong *) flint_realloc(A->exps,
                                                      new_alloc*sizeof(ulong));
            A->coeffs = (fmpz_mpoly_struct *) flint_realloc(A->coeffs,
                                          new_alloc*sizeof(fmpz_mpoly_struct));
        }

        for (i = old_alloc; i < new_alloc; i++)
        {
            fmpz_mpoly_init(A->coeffs + i, uctx);
            fmpz_mpoly_fit_bits(A->coeffs + i, A->bits, uctx);
            (A->coeffs + i)->bits = A->bits;
        }
        A->alloc = new_alloc;
    }
}

void fmpz_mpolyu_one(fmpz_mpolyu_t A, const fmpz_mpoly_ctx_t uctx)
{
    slong i;
    fmpz_mpolyu_fit_length(A, WORD(1), uctx);
    A->exps[0] = UWORD(0);
    fmpz_mpoly_one(A->coeffs + 0, uctx);
    A->length = WORD(1);
}


void fmpz_mpolyu_set(fmpz_mpolyu_t A, const fmpz_mpolyu_t B, const fmpz_mpoly_ctx_t uctx)
{
    slong i;
    fmpz_mpoly_struct * Acoeff, * Bcoeff;
    ulong * Aexp, * Bexp;
    slong Alen, Blen;

    Alen = 0;
    Blen = B->length;
    fmpz_mpolyu_fit_length(A, Blen, uctx);
    Acoeff = A->coeffs;
    Bcoeff = B->coeffs;
    Aexp = A->exps;
    Bexp = B->exps;

    for (i = 0; i < Blen; i++)
    {
        fmpz_mpoly_set(Acoeff + Alen, Bcoeff + i, uctx);
        Aexp[Alen++] = Bexp[i];
    }
    Alen = Blen;

    /* demote remaining coefficients */
    for (i = Alen; i < A->length; i++)
    {
        fmpz_mpoly_clear(Acoeff + i, uctx);
        fmpz_mpoly_init(Acoeff + i, uctx);
    }
    A->length = Alen;
}


/* if the coefficient doesn't exist, a new one is created (and set to zero) */
fmpz_mpoly_struct * _fmpz_mpolyu_get_coeff(fmpz_mpolyu_t A,
                             ulong pow, const fmpz_mpoly_ctx_t uctx)
{
    slong i, j;
    fmpz_mpoly_struct * xk;

    for (i = 0; i < A->length && A->exps[i] >= pow; i++)
    {
//flint_printf(" get coeff i = %wd\n",i);
        if (A->exps[i] == pow) 
        {
            return A->coeffs + i;
        }
    }

    fmpz_mpolyu_fit_length(A, A->length + 1, uctx);

//flint_printf("bits = %wd\n",bits);


    for (j = A->length; j > i; j--)
    {
        A->exps[j] = A->exps[j - 1];
        fmpz_mpoly_swap(A->coeffs + j, A->coeffs + j - 1, uctx);
    }

//flint_printf("bits = %wd\n",bits);
    
    A->length++;
    A->exps[i] = pow;
    xk = A->coeffs + i;
    xk->length = 0;
    FLINT_ASSERT(xk->bits == A->bits);
//    nmod_mpoly_fit_bits(xk, A->bits, uctx);
//    xk->bits = A->bits;

//flint_printf("xk->bits = %wd\n",xk->bits);

    return xk;
}


void fmpz_mpoly_to_mpolyu_perm(fmpz_mpolyu_t A, const fmpz_mpoly_t B,
                                                    const slong * perm,
                                                    const fmpz_mpoly_ctx_t uctx,
                                                    const fmpz_mpoly_ctx_t ctx)
{
    slong i, j;
    TMP_INIT;

    slong n = ctx->minfo->nvars;
    FLINT_ASSERT(uctx->minfo->nvars == n - 1);

    TMP_START;

    fmpz * uexps = (fmpz *) TMP_ALLOC(n*sizeof(fmpz));
    fmpz * exps  = (fmpz *) TMP_ALLOC(n*sizeof(fmpz));
    for (i = 0; i < n; i++)
    {
        fmpz_init(uexps + i);
        fmpz_init(exps + i);
    }

    A->bits = B->bits;
    fmpz_mpolyu_zero(A, uctx);

    slong N = mpoly_words_per_exp(B->bits, ctx->minfo);
    slong NA = mpoly_words_per_exp(A->bits, uctx->minfo);

    for (j = 0; j < B->length; j++)
    {
        mpoly_get_monomial_ffmpz(exps, B->exps + N*j, B->bits, ctx->minfo);

//flint_printf("original exp:");for (i=0; i<n; i++){flint_printf(", %wd", exps[i]);}printf("\n");

        for (i = 0; i < n; i++)
        {
            fmpz_swap(uexps + i, exps + perm[i]);
        }
        fmpz_mpoly_struct * Ac = _fmpz_mpolyu_get_coeff(A, uexps[n - 1], uctx);

        FLINT_ASSERT(Ac->bits == B->bits);


//flint_printf("        uexp:");for (i=0; i<n; i++){flint_printf(", %wd", uexps[i]);}printf("\n");

        fmpz_mpoly_fit_length(Ac, Ac->length + 1, uctx);
        fmpz_set(Ac->coeffs + Ac->length, B->coeffs + j);
        mpoly_set_monomial_ffmpz(Ac->exps + NA*Ac->length, uexps, A->bits, uctx->minfo);
        Ac->length++;

//flint_printf("A after: "); nmod_mpolyu_print_pretty(A, NULL, uctx); printf("\n");

    }



    for (i = 0; i < A->length; i++)
    {
        fmpz_mpoly_sort(A->coeffs + i, uctx);
    }

    for (i = 0; i < n; i++)
    {
        fmpz_clear(uexps + i);
        fmpz_clear(exps + i);
    }

    TMP_END;
}

/*
    
*/
void fmpz_mpoly_from_mpolyu_perm(fmpz_mpoly_t A, const fmpz_mpolyu_t B,
                                                    const slong * perm,
                                                    const fmpz_mpoly_ctx_t uctx,
                                                    const fmpz_mpoly_ctx_t ctx)
{
    slong i, j, k, bits, N;
    slong Alen;
    fmpz * Acoeff;
    ulong * Aexp;
    slong Aalloc;
    slong n = ctx->minfo->nvars;
    TMP_INIT;

    FLINT_ASSERT(uctx->minfo->nvars == n - 1);

    if (B->length == 0)
    {
        fmpz_mpoly_zero(A, ctx);
        return;        
    }

    TMP_START;

    /* find bits required to represent result */
    fmpz * texps = (fmpz *) TMP_ALLOC(n*sizeof(fmpz));
    fmpz * uexps = (fmpz *) TMP_ALLOC(n*sizeof(fmpz));
    fmpz * exps  = (fmpz *) TMP_ALLOC(n*sizeof(fmpz));
    for (i = 0; i < n; i++)
    {
        fmpz_init(texps + i);
        fmpz_init(uexps + i);
        fmpz_init(exps + i);
    }
    for (i = 0; i < B->length; i++)
    {
        fmpz_mpoly_struct * pi = B->coeffs + i;
        mpoly_degrees_ffmpz(texps, pi->exps, pi->length, pi->bits, uctx->minfo);
        _fmpz_vec_max_inplace(uexps, texps, n - 1);
    }
    fmpz_set_ui(uexps + n - 1, B->exps[0]);
    for (i = 0; i < n; i++)
    {
        fmpz_swap(uexps + i, exps + perm[i]);
    }
    bits = mpoly_exp_bits_required_ffmpz(exps, ctx->minfo);
    bits = FLINT_MAX(MPOLY_MIN_BITS, bits);
    bits = mpoly_fix_bits(bits, ctx->minfo);
/*
flint_printf("bits: %wd\n", bits);
*/
    N = mpoly_words_per_exp(bits, ctx->minfo);

    fmpz_mpoly_fit_bits(A, bits, ctx);
    A->bits = bits;

    Acoeff = A->coeffs;
    Aexp = A->exps;
    Aalloc = A->alloc;
    Alen = 0;
    for (i = 0; i < B->length; i++)
    {
        fmpz_mpoly_struct * Bc = B->coeffs + i;
        _fmpz_mpoly_fit_length(&Acoeff, &Aexp, &Aalloc, Alen + Bc->length, N);
        slong NB = mpoly_words_per_exp(Bc->bits, uctx->minfo);
        for (j = 0; j < Bc->length; j++)
        {
            fmpz_set(Acoeff + Alen + j, Bc->coeffs + j);
            mpoly_get_monomial_ffmpz(uexps, Bc->exps + NB*j, Bc->bits, uctx->minfo);
            fmpz_set_ui(uexps + n - 1, B->exps[i]);

            for (k = 0; k < n; k++)
            {
                fmpz_swap(uexps + k, exps + perm[k]);
            }

            mpoly_set_monomial_ffmpz(Aexp + N*(Alen + j), exps, bits, ctx->minfo);
        }
        Alen += (B->coeffs + i)->length;
    }
    A->coeffs = Acoeff;
    A->exps = Aexp;
    A->alloc = Aalloc;
    _fmpz_mpoly_set_length(A, Alen, ctx);
/*
printf("poly: "); nmod_mpoly_print_pretty(A, NULL, ctx); printf("\n");
*/
    for (i = 0; i < n; i++)
    {
        fmpz_clear(texps + i);
        fmpz_clear(uexps + i);
        fmpz_clear(exps + i);
    }

    fmpz_mpoly_sort(A, ctx);
    TMP_END;
}




void fmpz_mpoly_from_mpolyu_perm_keepbits(fmpz_mpoly_t A, const fmpz_mpolyu_t B,
                                                    const slong * perm,
                                                    const fmpz_mpoly_ctx_t uctx,
                                                    const fmpz_mpoly_ctx_t ctx)
{
    slong i, j, k, bits, N;
    slong Alen;
    fmpz * Acoeff;
    ulong * Aexp;
    slong Aalloc;
    slong n = ctx->minfo->nvars;
    TMP_INIT;

    FLINT_ASSERT(uctx->minfo->nvars == n - 1);

    if (B->length == 0)
    {
        fmpz_mpoly_zero(A, ctx);
        return;        
    }

    TMP_START;

    /* find bits required to represent result */
    fmpz * texps = (fmpz *) TMP_ALLOC(n*sizeof(fmpz));
    fmpz * uexps = (fmpz *) TMP_ALLOC(n*sizeof(fmpz));
    fmpz * exps  = (fmpz *) TMP_ALLOC(n*sizeof(fmpz));
    for (i = 0; i < n; i++)
    {
        fmpz_init(texps + i);
        fmpz_init(uexps + i);
        fmpz_init(exps + i);
    }
    for (i = 0; i < B->length; i++)
    {
        fmpz_mpoly_struct * pi = B->coeffs + i;
        mpoly_degrees_ffmpz(texps, pi->exps, pi->length, pi->bits, uctx->minfo);
        _fmpz_vec_max_inplace(uexps, texps, n - 1);
    }
    fmpz_set_ui(uexps + n - 1, B->exps[0]);
    for (i = 0; i < n; i++)
    {
        fmpz_swap(uexps + i, exps + perm[i]);
    }
    bits = mpoly_exp_bits_required_ffmpz(exps, ctx->minfo);

    FLINT_ASSERT(bits <= B->bits);
    bits = B->bits;
/*
flint_printf("bits: %wd\n", bits);
*/
    N = mpoly_words_per_exp(bits, ctx->minfo);

    fmpz_mpoly_fit_bits(A, bits, ctx);
    A->bits = bits;

    Acoeff = A->coeffs;
    Aexp = A->exps;
    Aalloc = A->alloc;
    Alen = 0;
    for (i = 0; i < B->length; i++)
    {
        fmpz_mpoly_struct * Bc = B->coeffs + i;
        _fmpz_mpoly_fit_length(&Acoeff, &Aexp, &Aalloc, Alen + Bc->length, N);
        slong NB = mpoly_words_per_exp(Bc->bits, uctx->minfo);
        for (j = 0; j < Bc->length; j++)
        {
            fmpz_set(Acoeff + Alen + j, Bc->coeffs + j);
            mpoly_get_monomial_ffmpz(uexps, Bc->exps + NB*j, Bc->bits, uctx->minfo);
            fmpz_set_ui(uexps + n - 1, B->exps[i]);

            for (k = 0; k < n; k++)
            {
                fmpz_swap(uexps + k, exps + perm[k]);
            }

            mpoly_set_monomial_ffmpz(Aexp + N*(Alen + j), exps, bits, ctx->minfo);
        }
        Alen += (B->coeffs + i)->length;
    }
    A->coeffs = Acoeff;
    A->exps = Aexp;
    A->alloc = Aalloc;
    _fmpz_mpoly_set_length(A, Alen, ctx);
/*
printf("poly: "); nmod_mpoly_print_pretty(A, NULL, ctx); printf("\n");
*/
    for (i = 0; i < n; i++)
    {
        fmpz_clear(texps + i);
        fmpz_clear(uexps + i);
        fmpz_clear(exps + i);
    }

    fmpz_mpoly_sort(A, ctx);
    TMP_END;
}







void fmpz_mpoly_to_nmod_mpoly(
    nmod_mpoly_t Ap,
    nmod_mpoly_ctx_t ctxp,
    fmpz_mpoly_t A,
    fmpz_mpoly_ctx_t ctx)
{
    fmpz_t t;
    fmpz_init(t);
    FLINT_ASSERT(Ap->bits == A->bits);
    slong N = mpoly_words_per_exp(A->bits, ctx->minfo);
    nmod_mpoly_fit_length(Ap, A->length, ctxp);
    slong k = 0;
    for (slong i = 0; i < A->length; i++)
    {
        mpoly_monomial_set(Ap->exps + N*k, A->exps + N*i, N);
        fmpz_set_ui(t, ctxp->ffinfo->mod.n);
        fmpz_mod(t, A->coeffs + i, t);
        Ap->coeffs[k] = fmpz_get_ui(t);
        k += (Ap->coeffs[k] != UWORD(0));
    }
    Ap->length = k;
    fmpz_clear(t);
}
void fmpz_mpolyu_to_nmod_mpolyu(
    nmod_mpolyu_t Ap,
    nmod_mpoly_ctx_t ctxp,
    fmpz_mpolyu_t A,
    fmpz_mpoly_ctx_t ctx)
{
    FLINT_ASSERT(Ap->bits == A->bits);
    nmod_mpolyu_fit_length(Ap, A->length, ctxp);
    for (slong i = 0; i < A->length; i++)
    {
        Ap->exps[i] = A->exps[i];
        fmpz_mpoly_to_nmod_mpoly(Ap->coeffs + i, ctxp, A->coeffs + i, ctx);
    }
    Ap->length = A->length;
}






void fmpz_mpoly_set_nmod_mpoly(
    fmpz_mpoly_t A,
    fmpz_mpoly_ctx_t ctx,
    nmod_mpoly_t Ap,
    nmod_mpoly_ctx_t ctxp)
{
    slong i, N;
    FLINT_ASSERT(Ap->bits == A->bits);
    fmpz_mpoly_fit_length(A, Ap->length, ctx);
    N = mpoly_words_per_exp(A->bits, ctx->minfo);
    for (i = 0; i < Ap->length*N; i++)
         A->exps[i] = Ap->exps[i];
    _fmpz_vec_set_nmod_vec(A->coeffs, Ap->coeffs, Ap->length, ctxp->ffinfo->mod);
    A->length = Ap->length;
}
void fmpz_mpolyu_set_nmod_mpolyu(
    fmpz_mpolyu_t A,
    fmpz_mpoly_ctx_t ctx,
    nmod_mpolyu_t Ap,
    nmod_mpoly_ctx_t ctxp)
{
    slong i;
    FLINT_ASSERT(Ap->bits == A->bits);
    fmpz_mpolyu_fit_length(A, Ap->length, ctx);
    for (i = 0; i < Ap->length; i++)
    {
        A->exps[i] = Ap->exps[i];
        fmpz_mpoly_set_nmod_mpoly(A->coeffs + i, ctx, Ap->coeffs + i, ctxp);
    }
    A->length = Ap->length;
}


fmpz * fmpz_mpoly_leadcoeff_ref(fmpz_mpoly_t A)
{
    FLINT_ASSERT(A->length > 0);
    return A->coeffs + 0;
}

fmpz * fmpz_mpolyu_leadcoeff_ref(fmpz_mpolyu_t A)
{
    FLINT_ASSERT(A->length > 0);
    return fmpz_mpoly_leadcoeff_ref(A->coeffs + 0);
}




int fmpz_mpoly_CRT_nmod_mpoly(
    fmpz_mpoly_t H,
    fmpz_mpoly_ctx_t ctx,
    fmpz_t m,
    nmod_mpoly_t A,
    nmod_mpoly_ctx_t ctxp)
{
    int changed = 0;
    slong i, N;
    fmpz_t t;
    FLINT_ASSERT(H->length == A->length);
    FLINT_ASSERT(H->bits == A->bits);
    N = mpoly_words_per_exp(A->bits, ctx->minfo);
    fmpz_init(t);
    for (i = 0; i < A->length; i++)
    {
        FLINT_ASSERT(mpoly_monomial_equal(H->exps + N*i, A->exps + N*i, N));
        fmpz_CRT_ui(t, H->coeffs + i, m, A->coeffs[i], ctxp->ffinfo->mod.n, 1);
        changed |= !fmpz_equal(t, H->coeffs + i);
        fmpz_swap(t, H->coeffs + i);
    }
    fmpz_clear(t);
    return changed;
}
int fmpz_mpolyu_CRT_nmod_mpolyu(
    fmpz_mpolyu_t H,
    fmpz_mpoly_ctx_t ctx,
    fmpz_t m,
    nmod_mpolyu_t A,
    nmod_mpoly_ctx_t ctxp)
{
    int changed = 0;
    FLINT_ASSERT(H->bits == A->bits);
    FLINT_ASSERT(H->length == A->length);;
    for (slong i = 0; i < A->length; i++)
    {
        FLINT_ASSERT(H->exps[i] == A->exps[i]);
        changed |= fmpz_mpoly_CRT_nmod_mpoly(H->coeffs + i, ctx, m, A->coeffs + i, ctxp);
    }
    H->length = A->length;
    return changed;
}








void fmpz_mpolyu_msub(
    fmpz_mpolyu_t R,
    fmpz_mpolyu_t A,
    fmpz_mpolyu_t B,
    fmpz_mpoly_t c,
    slong e,
    const fmpz_mpoly_ctx_t ctx)
{
    fmpz_mpolyu_fit_length(R, A->length + B->length, ctx);
    fmpz_mpoly_t T;
    fmpz_mpoly_init(T, ctx);

    slong i = 0, j = 0, k = 0;
    while (i < A->length || j < B->length)
    {
        if (i < A->length && (j >= B->length || A->exps[i] > B->exps[j] + e))
        {
            /* only A ok */
            fmpz_mpoly_set(R->coeffs + k, A->coeffs + i, ctx);
            R->exps[k] = A->exps[i];
            k++;
            i++;
        }
        else if (j < B->length && (i >= A->length || B->exps[j] + e > A->exps[i]))
        {
            /* only B ok */
            fmpz_mpoly_mul_johnson(R->coeffs + k, B->coeffs + j, c, ctx);
            fmpz_mpoly_neg(R->coeffs + k, R->coeffs + k, ctx);
            R->exps[k] = B->exps[j] + e;
            k++;
            j++;
        }
        else if (i < A->length && j < B->length && (A->exps[i] == B->exps[j] + e))
        {
            fmpz_mpoly_mul_johnson(T, B->coeffs + j, c, ctx);

//printf("T: "); nmod_mpoly_print_pretty(T, NULL, ctx); printf("\n");
//printf("A(%d): ",i); nmod_mpoly_print_pretty(A->coeffs + i, NULL, ctx); printf("\n");

            fmpz_mpoly_sub(R->coeffs + k, A->coeffs + i, T, ctx);

//printf("R(%d): ",k); nmod_mpoly_print_pretty(R->coeffs + k, NULL, ctx); printf("\n");


            R->exps[k] = A->exps[i];
            k += !fmpz_mpoly_is_zero(R->coeffs + k, ctx);
            i++;
            j++;
        } else 
        {
            FLINT_ASSERT(0);
        }
    }

    fmpz_mpoly_clear(T, ctx);
    R->length = k;
}

int fmpz_mpolyu_divides(fmpz_mpolyu_t A, fmpz_mpolyu_t B, const fmpz_mpoly_ctx_t ctx)
{
    int ret = 0;
    fmpz_mpolyu_t P, R;
    fmpz_mpoly_t t;
    fmpz_mpoly_init(t, ctx);
    fmpz_mpolyu_init(P, A->bits, ctx);
    fmpz_mpolyu_init(R, A->bits, ctx);
    fmpz_mpolyu_set(R, A, ctx);

//flint_printf("testing divisiblity:\n");
//printf("A: "); nmod_mpolyu_print_pretty(A, NULL, ctx); printf("\n");
//printf("B: "); nmod_mpolyu_print_pretty(B, NULL, ctx); printf("\n");

    FLINT_ASSERT(B->length > 0);

    while (R->length > 0)
    {
        if (R->exps[0] < B->exps[0])
            goto done;

        if (!fmpz_mpoly_divides_monagan_pearce(t, R->coeffs + 0, B->coeffs + 0, ctx))
            goto done;

//printf("t: "); nmod_mpoly_print_pretty(t, NULL, ctx); printf("\n");


        fmpz_mpolyu_msub(P, R, B, t, R->exps[0] - B->exps[0], ctx);
        fmpz_mpolyu_swap(P, R, ctx);

//printf("R: "); nmod_mpolyu_print_pretty(R, NULL, ctx); printf("\n");

    }
    ret = 1;

done:
    fmpz_mpoly_clear(t, ctx);
    fmpz_mpolyu_clear(P, ctx);
    fmpz_mpolyu_clear(R, ctx);

//flint_printf("returning divisiblity: %d\n", ret);

    return ret;
}



void fmpz_mpoly_remove_fmpz_content(fmpz_mpoly_t A, fmpz_mpoly_t B, fmpz_t g, fmpz_mpoly_ctx_t ctx)
{
    slong i;
    slong N;

    FLINT_ASSERT(B->bits == A->bits);

    fmpz_mpoly_fit_length(A, B->length, ctx);
    N = mpoly_words_per_exp(B->bits, ctx->minfo);
    for (i = 0; i < B->length; i++)
    {
        mpoly_monomial_set(A->exps + N*i, B->exps + N*i, N);
        fmpz_divexact(A->coeffs + i, B->coeffs + i, g);
    }
    A->length = B->length;
}

void fmpz_mpolyu_remove_fmpz_content(fmpz_mpolyu_t A, fmpz_mpolyu_t B, fmpz_mpoly_ctx_t ctx)
{
    slong i, j;
    fmpz_t g;

    FLINT_ASSERT(A->bits == B->bits);
    fmpz_mpolyu_fit_length(A, B->length, ctx);

    fmpz_init_set_si(g, WORD(0));

    for (i = 0; i < B->length; i++)
    {
        for (j = 0; j < (B->coeffs + i)->length; j++)
        {
            fmpz_gcd(g, g, (B->coeffs + i)->coeffs + j);
            if (fmpz_is_one(g))
                break;
        }
    }

    FLINT_ASSERT(!fmpz_is_zero(g));

    for (i = 0; i < B->length; i++)
    {
        A->exps[i] = B->exps[i];
        fmpz_mpoly_remove_fmpz_content(A->coeffs + i, B->coeffs + i, g, ctx);
    }
    A->length = B->length;

    fmpz_clear(g);
}


void nmod_mpoly_ctx_change_modulus(nmod_mpoly_ctx_t ctx, mp_limb_t modulus)
{
    nmodf_ctx_clear(ctx->ffinfo);
    nmodf_ctx_init(ctx->ffinfo, modulus);
}



int fmpz_mpolyu_gcd_zippel_linzipm(
    fmpz_mpolyu_t G,
    fmpz_mpolyu_t A,
    fmpz_mpolyu_t B,
    fmpz_mpoly_ctx_t ctx,
    mpoly_zipinfo_t zinfo,
    flint_rand_t randstate)
{
    int success, changed;
    mp_limb_t p = UWORD(1) << (FLINT_BITS - 1), old_p, t, gammap;
    fmpz_t gamma, pp, gammapp, m;
    nmod_mpolyu_t Ap, Bp, Gp, Gform;
    fmpz_mpolyu_t H;
    nmod_mpoly_ctx_t ctxp;

    p=12;
/*
flint_printf("**********fmpz_mpolyu_gcd_zippel_linzipm *******\n");
flint_printf("A: "); fmpz_mpolyu_print_pretty(A, NULL, ctx); flint_printf("\n");
flint_printf("B: "); fmpz_mpolyu_print_pretty(B, NULL, ctx); flint_printf("\n");
*/
    fmpz_init(pp);
    fmpz_init(gammapp);
    fmpz_init_set_si(m, 1);
    fmpz_init(gamma);
    fmpz_gcd(gamma, fmpz_mpolyu_leadcoeff_ref(A), fmpz_mpolyu_leadcoeff_ref(B));

    FLINT_ASSERT(A->bits == B->bits);
    FLINT_ASSERT(A->bits == G->bits);
    FLINT_ASSERT(A->length > 0);
    FLINT_ASSERT(B->length > 0);
    FLINT_ASSERT(A->exps[A->length - 1] == 0);
    FLINT_ASSERT(B->exps[B->length - 1] == 0);

    nmod_mpoly_ctx_init(ctxp, ctx->minfo->nvars, ORD_LEX, 2);

    nmod_mpolyu_init(Ap, A->bits, ctxp);
    nmod_mpolyu_init(Bp, A->bits, ctxp);
    nmod_mpolyu_init(Gp, A->bits, ctxp);
    nmod_mpolyu_init(Gform, A->bits, ctxp);

    fmpz_mpolyu_init(H, A->bits, ctx);

choose_prime:
    old_p = p;
    p = n_nextprime(p, 1);
    if (p <= old_p) {
        /* ran out of primes */
        success = 0;
        goto finished;
    }


//flint_printf("outer p: %wd\n", p);


    fmpz_set_ui(pp, p);
    fmpz_mod(gammapp, gamma, pp);
    gammap = fmpz_get_ui(gammapp);
    if (gammap == UWORD(0))
        goto choose_prime;

//flint_printf("gamma = "); fmpz_print(gamma); printf("\n");
//flint_printf("gamma mod p = %d\n", gammap);

    nmod_mpoly_ctx_change_modulus(ctxp, p);

    fmpz_mpolyu_to_nmod_mpolyu(Ap, ctxp, A, ctx);
    fmpz_mpolyu_to_nmod_mpolyu(Bp, ctxp, B, ctx);
    success = nmod_mpolyu_gcd_zippel_linzipp(Gp, Ap, Bp, ctx->minfo->nvars - 1, ctxp, zinfo, randstate);
    if (success != 1)
        goto choose_prime;
    t = nmod_mpolyu_leadcoeff(Gp, ctxp);
    t = nmod_inv(t, ctxp->ffinfo->mod);
    t = nmod_mul(t, gammap, ctxp->ffinfo->mod);
    nmod_mpolyu_scalar_mul_nmod(Gp, t, ctxp);

//    printf("           Gp: "); nmod_mpolyu_print_pretty(Gp, NULL, ctxp); printf("\n");

    if (Gp->length == 1) {
        FLINT_ASSERT(Gp->exps[0] == 0);
        FLINT_ASSERT(nmod_mpoly_is_nmod(Gp->coeffs + 0, ctxp));
        FLINT_ASSERT(!nmod_mpoly_is_zero(Gp->coeffs + 0, ctxp));
        fmpz_mpolyu_one(G, ctx);
        success = 1;
        goto finished;
    }


    nmod_mpolyu_setform(Gform, Gp, ctxp);
    fmpz_mpolyu_set_nmod_mpolyu(H, ctx, Gp, ctxp);
    fmpz_set_ui(m, p);

/*
    printf("fmpz H: "); fmpz_mpolyu_print_pretty(H, NULL, ctx); printf("\n");

flint_printf("**********fmpz_mpolyu_gcd_zippel_linzipm entering inner loop *******\n");
flint_printf("A: "); fmpz_mpolyu_print_pretty(A, NULL, ctx); flint_printf("\n");
flint_printf("B: "); fmpz_mpolyu_print_pretty(B, NULL, ctx); flint_printf("\n");
*/

choose_prime_inner:
    old_p = p;
    p = n_nextprime(p, 1);
    if (p <= old_p) {
        /* ran out of primes */
        success = 0;
        goto finished;
    }

//flint_printf("inner p: %wd\n", p);

    fmpz_set_ui(pp, p);
    fmpz_mod(gammapp, gamma, pp);
    gammap = fmpz_get_ui(gammapp);
    if (gammap == UWORD(0))
        goto choose_prime_inner;

//flint_printf("gamma = "); fmpz_print(gamma); printf("\n");
//flint_printf("gamma mod p = %d\n", gammap);


    nmod_mpoly_ctx_change_modulus(ctxp, p);

    fmpz_mpolyu_to_nmod_mpolyu(Ap, ctxp, A, ctx);
    fmpz_mpolyu_to_nmod_mpolyu(Bp, ctxp, B, ctx);

    success = nmod_mpolyu_sgcd_zippel(Gp, Ap, Bp, Gform, ctx->minfo->nvars, ctxp, randstate);
    if (success == 0)
        goto choose_prime;
    if (success == -1)
        goto choose_prime_inner;

    t = nmod_mpolyu_leadcoeff(Gp, ctxp);
    t = nmod_inv(t, ctxp->ffinfo->mod);
    t = nmod_mul(t, gammap, ctxp->ffinfo->mod);
    nmod_mpolyu_scalar_mul_nmod(Gp, t, ctxp);
/*
    printf("           Gp: "); nmod_mpolyu_print_pretty(Gp, NULL, ctxp); printf("\n");

    printf("before fmpz H: "); fmpz_mpolyu_print_pretty(H, NULL, ctx); printf("\n");
*/

    changed = fmpz_mpolyu_CRT_nmod_mpolyu(H, ctx, m, Gp, ctxp);

    fmpz_mul_ui(m, m, p);
//    printf(" after fmpz H: "); fmpz_mpolyu_print_pretty(H, NULL, ctx); printf("\n");


    if (changed)
        goto choose_prime_inner;

/*
    flint_printf("nvars: %d\n", ctx->minfo->nvars);
    flint_printf("          A: "); fmpz_mpolyu_print_pretty(A, NULL, ctx); flint_printf("\n");
    flint_printf("          B: "); fmpz_mpolyu_print_pretty(B, NULL, ctx); flint_printf("\n");
    flint_printf("candidate G: "); fmpz_mpolyu_print_pretty(H, NULL, ctx); printf("\n");
*/
    fmpz_mpolyu_remove_fmpz_content(G, H, ctx);

    if (!fmpz_mpolyu_divides(A, G, ctx) || !fmpz_mpolyu_divides(B, G, ctx))
    {
//        flint_printf("fmpz division failed\n");
        goto choose_prime_inner;
    }
    success = 1;
    goto finished;




finished:

    nmod_mpolyu_clear(Ap, ctxp);
    nmod_mpolyu_clear(Bp, ctxp);
    nmod_mpolyu_clear(Gp, ctxp);
    nmod_mpolyu_clear(Gform, ctxp);

    fmpz_mpolyu_clear(H, ctx);

    nmod_mpoly_ctx_clear(ctxp);

    fmpz_clear(gammapp);
    fmpz_clear(gamma);
    fmpz_clear(m);
    fmpz_clear(pp);

/*
flint_printf("**********fmpz_mpolyu_gcd_zippel_linzipm returning %d *******\n", success);
flint_printf("A: "); fmpz_mpolyu_print_pretty(A, NULL, ctx); flint_printf("\n");
flint_printf("B: "); fmpz_mpolyu_print_pretty(B, NULL, ctx); flint_printf("\n");
flint_printf("G: "); fmpz_mpolyu_print_pretty(G, NULL, ctx); flint_printf("\n");
*/

    return success;
}







slong _fmpz_mpoly_divides_monagan_pearce(fmpz ** poly1, ulong ** exp1,
         slong * alloc, const fmpz * poly2, const ulong * exp2, slong len2,
       const fmpz * poly3, const ulong * exp3, slong len3, mp_bitcnt_t bits, slong N,
                                                         const ulong * cmpmask);

void fmpz_mpolyu_divexact_mpoly(fmpz_mpolyu_t A, fmpz_mpolyu_t B, fmpz_mpoly_t c, fmpz_mpoly_ctx_t ctx)
{
    slong i;
    slong len;
    slong N;
    mp_bitcnt_t exp_bits;
    ulong * cmpmask;
    TMP_INIT;

    TMP_START;

    exp_bits = B->bits;
    FLINT_ASSERT(A->bits == B->bits);

    fmpz_mpolyu_fit_length(A, B->length, ctx);

    N = mpoly_words_per_exp(exp_bits, ctx->minfo);
    cmpmask = (ulong*) TMP_ALLOC(N*sizeof(ulong));
    mpoly_get_cmpmask(cmpmask, N, exp_bits, ctx->minfo);

    for (i = 0; i < B->length; i++)
    {
        fmpz_mpoly_struct * poly1 = A->coeffs + i;
        fmpz_mpoly_struct * poly2 = B->coeffs + i;
        fmpz_mpoly_struct * poly3 = c;
        A->exps[i] = B->exps[i];

        fmpz_mpoly_fit_length(poly1, poly2->length/poly3->length + 1, ctx);
        fmpz_mpoly_fit_bits(poly1, exp_bits, ctx);
        poly1->bits = exp_bits;

        len = _fmpz_mpoly_divides_monagan_pearce(&poly1->coeffs, &poly1->exps,
                            &poly1->alloc, poly2->coeffs, poly2->exps, poly2->length,
                              poly3->coeffs, poly3->exps, poly3->length, exp_bits, N,
                                                  cmpmask);
        FLINT_ASSERT(len > 0);
        _fmpz_mpoly_set_length(poly1, len, ctx);
    }
    A->length = B->length;

    TMP_END;
}


slong _fmpz_mpoly_mul_johnson(fmpz ** poly1, ulong ** exp1, slong * alloc,
                 const fmpz * poly2, const ulong * exp2, slong len2,
                 const fmpz * poly3, const ulong * exp3, slong len3,
                              mp_bitcnt_t bits, slong N, const ulong * cmpmask);

void fmpz_mpolyu_mul_mpoly(fmpz_mpolyu_t A, fmpz_mpolyu_t B, fmpz_mpoly_t c, fmpz_mpoly_ctx_t ctx)
{
    slong i;
    slong len;
    slong N;
    mp_bitcnt_t exp_bits;
    ulong * cmpmask;
    TMP_INIT;

    TMP_START;

    exp_bits = B->bits;
    FLINT_ASSERT(A->bits == B->bits);
    FLINT_ASSERT(A->bits == c->bits);

    fmpz_mpolyu_fit_length(A, B->length, ctx);

    N = mpoly_words_per_exp(exp_bits, ctx->minfo);
    cmpmask = (ulong*) TMP_ALLOC(N*sizeof(ulong));
    mpoly_get_cmpmask(cmpmask, N, exp_bits, ctx->minfo);

    for (i = 0; i < B->length; i++)
    {
        fmpz_mpoly_struct * poly1 = A->coeffs + i;
        fmpz_mpoly_struct * poly2 = B->coeffs + i;
        fmpz_mpoly_struct * poly3 = c;
        A->exps[i] = B->exps[i];

        fmpz_mpoly_fit_length(poly1, poly2->length/poly3->length + 1, ctx);
        fmpz_mpoly_fit_bits(poly1, exp_bits, ctx);
        poly1->bits = exp_bits;

        len = _fmpz_mpoly_mul_johnson(&poly1->coeffs, &poly1->exps,
                            &poly1->alloc, poly2->coeffs, poly2->exps, poly2->length,
                              poly3->coeffs, poly3->exps, poly3->length, exp_bits, N,
                                                  cmpmask);

        _fmpz_mpoly_set_length(poly1, len, ctx);

    }
    A->length = B->length;

    TMP_END;
}

int fmpz_mpoly_gcd_zippel_keepbits(fmpz_mpoly_t G, fmpz_mpoly_t A, fmpz_mpoly_t B, const fmpz_mpoly_ctx_t ctx);



void fmpz_mpolyu_shift_right(fmpz_mpolyu_t A, ulong s)
{
    slong i;
    for (i = 0; i < A->length; i++)
    {
        FLINT_ASSERT(A->exps[i] >= s);
        A->exps[i] -= s;
    }
}

void fmpz_mpolyu_shift_left(fmpz_mpolyu_t A, ulong s)
{
    slong i;
    for (i = 0; i < A->length; i++)
    {
        FLINT_ASSERT((slong)(A->exps[i] + s) >= 0);
        A->exps[i] += s;
    }
}

int fmpz_mpolyu_gcd_zippel(
    fmpz_mpolyu_t G,
    fmpz_mpolyu_t A,
    fmpz_mpolyu_t B,
    fmpz_mpoly_ctx_t ctx,
    mpoly_zipinfo_t zinfo,
    flint_rand_t randstate)
{
    int success;
    slong i;
    slong ABminshift;
    fmpz_mpoly_t content;
    fmpz_mpolyu_t Abar, Bbar, Gbar;

    FLINT_ASSERT(A->bits == B->bits);
    FLINT_ASSERT(A->length > 0);
    FLINT_ASSERT(B->length > 0);

    fmpz_mpoly_init(content, ctx);
    fmpz_mpolyu_init(Abar, A->bits, ctx);
    fmpz_mpolyu_init(Bbar, A->bits, ctx);
    fmpz_mpolyu_init(Gbar, A->bits, ctx);

    fmpz_mpoly_set(content, A->coeffs + 0, ctx);
    for (i = 1; i < A->length; i++)
    {
        success = fmpz_mpoly_gcd_zippel_keepbits(content, content, A->coeffs + i, ctx);
        if (!success)
            return 0;
    }
    for (i = 0; i < B->length; i++)
    {
        success = fmpz_mpoly_gcd_zippel_keepbits(content, content, B->coeffs + i, ctx);
        if (!success)
            return 0;
    }
/*
    flint_printf("(%d)linzipm A: ",ctx->minfo->nvars); fmpz_mpolyu_print_pretty(A, NULL, ctx); printf("\n");
    flint_printf("(%d)linzipm B: ",ctx->minfo->nvars); fmpz_mpolyu_print_pretty(B, NULL, ctx); printf("\n");
    flint_printf("(%d)  content: ",ctx->minfo->nvars); fmpz_mpoly_print_pretty(content, NULL, ctx); printf("\n");
*/
    fmpz_mpolyu_divexact_mpoly(Abar, A, content, ctx);
    fmpz_mpolyu_divexact_mpoly(Bbar, B, content, ctx);
/*
    flint_printf("(%d)linzipm Abar: ",ctx->minfo->nvars); fmpz_mpolyu_print_pretty(Abar, NULL, ctx); printf("\n");
    flint_printf("(%d)linzipm Bbar: ",ctx->minfo->nvars); fmpz_mpolyu_print_pretty(Bbar, NULL, ctx); printf("\n");
*/
    ABminshift = FLINT_MIN(Abar->exps[Abar->length - 1], Bbar->exps[Bbar->length - 1]);
    fmpz_mpolyu_shift_right(Abar, Abar->exps[Abar->length - 1]);
    fmpz_mpolyu_shift_right(Bbar, Bbar->exps[Bbar->length - 1]);

    success = fmpz_mpolyu_gcd_zippel_linzipm(Gbar, Abar, Bbar, ctx, zinfo, randstate);
    if (!success)
        return 0;

    fmpz_mpolyu_shift_left(Gbar, ABminshift);
    fmpz_mpolyu_mul_mpoly(G, Gbar, content, ctx);

    fmpz_mpolyu_clear(Abar, ctx);
    fmpz_mpolyu_clear(Bbar, ctx);
    fmpz_mpolyu_clear(Gbar, ctx);
    fmpz_mpoly_clear(content, ctx);

    return 1;
}










void fmpz_mpoly_to_fmpz_poly_keepbits(fmpz_poly_t A, slong * Ashift,
               const fmpz_mpoly_t B, slong var, const fmpz_mpoly_ctx_t ctx)
{
    slong i, shift, off, N;
    slong _Ashift = 0, len = B->length;
    fmpz * coeff = B->coeffs;
    ulong * exp = B->exps;
    mp_bitcnt_t bits = B->bits;

    N = mpoly_words_per_exp(bits, ctx->minfo);
    mpoly_gen_offset_shift(&off, &shift, var, N, bits, ctx->minfo);

    fmpz_poly_zero(A);
    if (len > 0)
    {
        ulong mask = (-UWORD(1)) >> (FLINT_BITS - bits);
        _Ashift = (exp[N*(len - 1)] >> shift) & mask;
        for (i = 0; i < len; i++)
        {
            ulong k = ((exp[N*i + off] >> shift) & mask) - _Ashift;
            FLINT_ASSERT(((slong)k) >= 0);
            fmpz_poly_set_coeff_fmpz(A, k, coeff + i);
        }
    }

    *Ashift = _Ashift;
}

void fmpz_mpoly_from_fmpz_poly_keepbits(fmpz_mpoly_t A, const fmpz_poly_t B,
                           slong Bshift, slong var, mp_bitcnt_t bits, const fmpz_mpoly_ctx_t ctx)
{
    slong shift, off, N;
    slong k;
    slong Alen;
    fmpz * Acoeff;
    ulong * Aexp;
    slong Aalloc;
    ulong * one;
    TMP_INIT;

    TMP_START;

    FLINT_ASSERT(!fmpz_poly_is_zero(B));
    FLINT_ASSERT(Bshift >= 0);
    FLINT_ASSERT(Bshift + fmpz_poly_degree(B) >= 0);
    FLINT_ASSERT(1 + FLINT_BIT_COUNT(Bshift + fmpz_poly_degree(B)) <= bits);
    
    N = mpoly_words_per_exp(bits, ctx->minfo);
    one = (ulong*) TMP_ALLOC(N*sizeof(ulong));
    mpoly_gen_oneexp_offset_shift(one, &off, &shift, var, N, bits, ctx->minfo);

    fmpz_mpoly_fit_bits(A, bits, ctx);
    A->bits = bits;

    Acoeff = A->coeffs;
    Aexp = A->exps;
    Aalloc = A->alloc;
    Alen = 0;
    for (k = fmpz_poly_degree(B); k >= 0; k--)
    {
        _fmpz_mpoly_fit_length(&Acoeff, &Aexp, &Aalloc, Alen + 1, N);
        mpoly_monomial_mul_si(Aexp + N*Alen, one, N, k + Bshift);
        fmpz_poly_get_coeff_fmpz(Acoeff + Alen, B, k);
        Alen += !fmpz_is_zero(Acoeff + Alen);
    }

    A->coeffs = Acoeff;
    A->exps = Aexp;
    A->alloc = Aalloc;
    _fmpz_mpoly_set_length(A, Alen, ctx);

    TMP_END;
}


int fmpz_mpoly_gcd_zippel_keepbits(fmpz_mpoly_t G, fmpz_mpoly_t A, fmpz_mpoly_t B, const fmpz_mpoly_ctx_t ctx)
{
    flint_rand_t randstate;
    int ret, success = 0;
    mpoly_zipinfo_t zinfo;
    fmpz_mpoly_ctx_t uctx;
    fmpz_mpolyu_t Au, Bu, Gu;

    fmpz_mpoly_t Acopy, Bcopy;
    fmpz_mpoly_init(Acopy, ctx);
    fmpz_mpoly_init(Bcopy, ctx);
    fmpz_mpoly_set(Acopy, A, ctx);
    fmpz_mpoly_set(Bcopy, B, ctx);

    FLINT_ASSERT(A->bits <= FLINT_BITS);
    FLINT_ASSERT(B->bits <= FLINT_BITS);
    FLINT_ASSERT(A->bits == B->bits);

//flint_printf("(%d) **** fmpz_mpoly_gcd_zippel_keepbits: \n",ctx->minfo->nvars);
//flint_printf("(%d) A(%wd): ",ctx->minfo->nvars, A->bits); fmpz_mpoly_print_pretty(A, NULL, ctx); printf("\n");
//flint_printf("(%d) B(%wd): ",ctx->minfo->nvars, A->bits); fmpz_mpoly_print_pretty(B, NULL, ctx); printf("\n");

    FLINT_ASSERT(!fmpz_mpoly_is_zero(A, ctx));
    FLINT_ASSERT(!fmpz_mpoly_is_zero(B, ctx));

    if (ctx->minfo->nvars == 1) {
        slong shiftA, shiftB;
        fmpz_poly_t a, b, g;
        fmpz_poly_init(a);
        fmpz_poly_init(b);
        fmpz_poly_init(g);
        fmpz_mpoly_to_fmpz_poly_keepbits(a, &shiftA, A, 0, ctx);
        fmpz_mpoly_to_fmpz_poly_keepbits(b, &shiftB, B, 0, ctx);
        fmpz_poly_gcd(g, a, b);
        fmpz_mpoly_from_fmpz_poly_keepbits(G, g, FLINT_MIN(shiftA, shiftB), 0, A->bits, ctx);
        fmpz_poly_clear(a);
        fmpz_poly_clear(b);
        fmpz_poly_clear(g);

//flint_printf("(%d) **** fmpz_mpoly_gcd_zippel_keepbits returning: \n",ctx->minfo->nvars);
//flint_printf("(%d) +uu+ A(%wd): ",ctx->minfo->nvars, Acopy->bits); fmpz_mpoly_print_pretty(Acopy, NULL, ctx); printf("\n");
//flint_printf("(%d) +uu+ B(%wd): ",ctx->minfo->nvars, Acopy->bits); fmpz_mpoly_print_pretty(Bcopy, NULL, ctx); printf("\n");
//flint_printf("(%d) +uu+ G(%wd): ",ctx->minfo->nvars, G->bits); fmpz_mpoly_print_pretty(G, NULL, ctx); printf("\n");

        fmpz_mpoly_clear(Acopy, ctx);
        fmpz_mpoly_clear(Bcopy, ctx);
        return 1;
    }


    flint_randinit(randstate);

//flint_printf("1A bits: %d\n", A->bits);

    mpoly_zipinfo_init(zinfo, ctx->minfo->nvars);
    fmpz_mpoly_degrees_si(zinfo->Adegs, A, ctx);
    fmpz_mpoly_degrees_si(zinfo->Bdegs, B, ctx);
    for (slong i = 0; i < ctx->minfo->nvars; i++)
    {
        zinfo->perm[i] = ctx->minfo->nvars - 1 - i;
        zinfo->perm[i] = i;
    }


    slong new_bits = FLINT_MAX(A->bits, B->bits);


//flint_printf("new_bits: %d\n", new_bits);


    fmpz_mpoly_ctx_init(uctx, ctx->minfo->nvars - 1, ORD_LEX);
    fmpz_mpolyu_init(Au, new_bits, uctx);
    fmpz_mpolyu_init(Bu, new_bits, uctx);
    fmpz_mpolyu_init(Gu, new_bits, uctx);


    fmpz_mpoly_to_mpolyu_perm(Au, A, zinfo->perm, uctx, ctx);
    fmpz_mpoly_to_mpolyu_perm(Bu, B, zinfo->perm, uctx, ctx);

//flint_printf("Au: "); fmpz_mpolyu_print_pretty(Au, NULL, uctx); flint_printf("\n");
//flint_printf("Bu: "); fmpz_mpolyu_print_pretty(Bu, NULL, uctx); flint_printf("\n");

    ret = fmpz_mpolyu_gcd_zippel(Gu, Au, Bu, uctx, zinfo, randstate);
    if (ret) {
        fmpz_mpoly_from_mpolyu_perm_keepbits(G, Gu, zinfo->perm, uctx, ctx);
        if (fmpz_sgn(G->coeffs + 0) < 0)
            fmpz_mpoly_neg(G, G, ctx);
        success = 1;
    }

//flint_printf("Gu: "); fmpz_mpolyu_print_pretty(Gu, NULL, uctx); flint_printf("\n");

    fmpz_mpolyu_clear(Au, uctx);
    fmpz_mpolyu_clear(Bu, uctx);
    fmpz_mpolyu_clear(Gu, uctx);
    fmpz_mpoly_ctx_clear(uctx);

    mpoly_zipinfo_clear(zinfo);

    flint_randclear(randstate);

//flint_printf("(%d) **** fmpz_mpoly_gcd_zippel_keepbits returning: \n",ctx->minfo->nvars);
//flint_printf("(%d) ++++ A(%wd): ",ctx->minfo->nvars, Acopy->bits); fmpz_mpoly_print_pretty(Acopy, NULL, ctx); printf("\n");
//flint_printf("(%d) ++++ B(%wd): ",ctx->minfo->nvars, Acopy->bits); fmpz_mpoly_print_pretty(Bcopy, NULL, ctx); printf("\n");
//flint_printf("(%d) ++++ G(%wd): ",ctx->minfo->nvars, G->bits); fmpz_mpoly_print_pretty(G, NULL, ctx); printf("\n");

    fmpz_mpoly_clear(Acopy, ctx);
    fmpz_mpoly_clear(Bcopy, ctx);

//    TMP_END;
    return success;
}




int fmpz_mpoly_gcd_zippel(fmpz_mpoly_t G, fmpz_mpoly_t A, fmpz_mpoly_t B, const fmpz_mpoly_ctx_t ctx)
{
    flint_rand_t randstate;
    int ret, success = 0;
    mpoly_zipinfo_t zinfo;
    fmpz_mpoly_ctx_t uctx;
    fmpz_mpolyu_t Au, Bu, Gu;
//    const char * vars[] = {"x","y","z"};

//flint_printf("\n\n **** fmpz_mpoly_gcd_zippel\n");
//flint_printf("  A: "); fmpz_mpoly_print_pretty(A, NULL, ctx); printf("\n");
//flint_printf("  B: "); fmpz_mpoly_print_pretty(B, NULL, ctx); printf("\n");


    if (fmpz_mpoly_is_zero(A, ctx))
    {
        if (fmpz_mpoly_is_zero(B, ctx))
        {
            fmpz_mpoly_zero(G, ctx);
            return 1;
        }
        if (fmpz_sgn(B->coeffs + 0) < 0)
        {
            fmpz_mpoly_neg(G, B, ctx);
            return 1;
        } else
        {
            fmpz_mpoly_set(G, B, ctx);
            return 1;
        }
    }

    if (fmpz_mpoly_is_zero(B, ctx))
    {
        if (fmpz_sgn(A->coeffs + 0) < 0)
        {
            fmpz_mpoly_neg(G, A, ctx);
            return 1;
        } else
        {
            fmpz_mpoly_set(G, A, ctx);
            return 1;
        }
    }

    if (ctx->minfo->nvars == 1) {
        slong shiftA, shiftB;
        fmpz_poly_t a, b, g;
        fmpz_poly_init(a);
        fmpz_poly_init(b);
        fmpz_poly_init(g);
        fmpz_mpoly_to_fmpz_poly_keepbits(a, &shiftA, A, 0, ctx);
        fmpz_mpoly_to_fmpz_poly_keepbits(b, &shiftB, B, 0, ctx);
        fmpz_poly_gcd(g, a, b);
        fmpz_mpoly_from_fmpz_poly_keepbits(G, g, FLINT_MIN(shiftA, shiftB), 0, A->bits, ctx);
        fmpz_poly_clear(a);
        fmpz_poly_clear(b);
        fmpz_poly_clear(g);
        return 1;
    }

    FLINT_ASSERT(A->bits <= FLINT_BITS);
    FLINT_ASSERT(B->bits <= FLINT_BITS);

//flint_printf("**************\nfmpz_mpoly_gcd_zippel: \n");
//flint_printf("A(%wd): ", A->bits); fmpz_mpoly_print_pretty(A, NULL, ctx); printf("\n");
//flint_printf("B(%wd): ", A->bits); fmpz_mpoly_print_pretty(B, NULL, ctx); printf("\n");

    FLINT_ASSERT(ctx->minfo->nvars > WORD(1));
    FLINT_ASSERT(!fmpz_mpoly_is_zero(A, ctx));
    FLINT_ASSERT(!fmpz_mpoly_is_zero(B, ctx));

    flint_randinit(randstate);

//flint_printf("1A bits: %d\n", A->bits);

    mpoly_zipinfo_init(zinfo, ctx->minfo->nvars);
    fmpz_mpoly_degrees_si(zinfo->Adegs, A, ctx);
    fmpz_mpoly_degrees_si(zinfo->Bdegs, B, ctx);
    for (slong i = 0; i < ctx->minfo->nvars; i++)
    {
        zinfo->perm[i] = ctx->minfo->nvars - 1 - i;
        zinfo->perm[i] = i;
    }

    slong new_bits = FLINT_MAX(A->bits, B->bits);

    fmpz_mpoly_ctx_init(uctx, ctx->minfo->nvars - 1, ORD_LEX);
    fmpz_mpolyu_init(Au, new_bits, uctx);
    fmpz_mpolyu_init(Bu, new_bits, uctx);
    fmpz_mpolyu_init(Gu, new_bits, uctx);


    fmpz_mpoly_to_mpolyu_perm(Au, A, zinfo->perm, uctx, ctx);
    fmpz_mpoly_to_mpolyu_perm(Bu, B, zinfo->perm, uctx, ctx);

//flint_printf("Au: "); fmpz_mpolyu_print_pretty(Au, NULL, uctx); flint_printf("\n");
//flint_printf("Bu: "); fmpz_mpolyu_print_pretty(Bu, NULL, uctx); flint_printf("\n");

    ret = fmpz_mpolyu_gcd_zippel(Gu, Au, Bu, uctx, zinfo, randstate);
    if (ret) {
        fmpz_mpoly_from_mpolyu_perm(G, Gu, zinfo->perm, uctx, ctx);
        if (fmpz_sgn(G->coeffs + 0) < 0)
            fmpz_mpoly_neg(G, G, ctx);
        success = 1;
    }

//flint_printf("Gu: "); fmpz_mpolyu_print_pretty(Gu, NULL, uctx); flint_printf("\n");

    fmpz_mpolyu_clear(Au, uctx);
    fmpz_mpolyu_clear(Bu, uctx);
    fmpz_mpolyu_clear(Gu, uctx);
    fmpz_mpoly_ctx_clear(uctx);

    mpoly_zipinfo_clear(zinfo);

    flint_randclear(randstate);
/*
flint_printf("\n\n **** fmpz_mpoly_gcd_zippel returning\n");
flint_printf("++++ G(%wd): ", G->bits); fmpz_mpoly_print_pretty(G, NULL, ctx); printf("\n");
*/
//    TMP_END;
    return success;
}
