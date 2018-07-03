/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "nmod_mpoly.h"




/************ zipinfo ********************/
void mpoly_zipinfo_init(mpoly_zipinfo_t zinfo, slong nvars)
{
    zinfo->nvars = nvars;
    zinfo->Adegs = (slong *) flint_malloc(nvars*sizeof(slong));
    zinfo->Bdegs = (slong *) flint_malloc(nvars*sizeof(slong));
    zinfo->perm  = (slong *) flint_malloc(nvars*sizeof(slong));
}

void mpoly_zipinfo_clear(mpoly_zipinfo_t zinfo)
{
    flint_free(zinfo->Adegs);
    flint_free(zinfo->Bdegs);
    flint_free(zinfo->perm);
}


/*********** mpolyu **************************/

void nmod_mpolyu_init(nmod_mpolyu_t A, mp_bitcnt_t bits, const nmod_mpoly_ctx_t ctx)
{
    A->coeffs = NULL;
    A->exps = NULL;
    A->alloc = 0;
    A->length = 0;
    A->bits = bits;
}


void nmod_mpolyu_clear(nmod_mpolyu_t A, const nmod_mpoly_ctx_t uctx)
{
    slong i;
    for (i = 0; i < A->alloc; i++)
        nmod_mpoly_clear(A->coeffs + i, uctx);
    flint_free(A->coeffs);
    flint_free(A->exps);
}


void nmod_mpolyu_swap(nmod_mpolyu_t A, nmod_mpolyu_t B, const nmod_mpoly_ctx_t uctx)
{
   nmod_mpolyu_struct t = *A;
   *A = *B;
   *B = t;
}

void nmod_mpolyu_zero(nmod_mpolyu_t A, const nmod_mpoly_ctx_t uctx)
{
    slong i;
    for (i = 0; i < A->alloc; i++) {
        nmod_mpoly_clear(A->coeffs + i, uctx);
        nmod_mpoly_init(A->coeffs + i, uctx);
    }
    A->length = 0;
}

int nmod_mpolyu_is_one(nmod_mpolyu_t A, const nmod_mpoly_ctx_t uctx)
{
    if (A->length != 1 || A->exps[0] != UWORD(0))
        return 0;

    return nmod_mpoly_is_one(A->coeffs + 0, uctx);
}

void nmod_mpolyu_print_pretty(const nmod_mpolyu_t poly,
                                   const char ** x, const nmod_mpoly_ctx_t ctx)
{
    slong i;
    if (poly->length == 0)
        flint_printf("0");
    for (i = 0; i < poly->length; i++)
    {
        if (i != 0)
            flint_printf(" + ");
        flint_printf("(");
        nmod_mpoly_print_pretty(poly->coeffs + i,x,ctx);
        flint_printf(")*X^%wd",poly->exps[i]);
    }
}

void nmod_mpolyu_fit_length(nmod_mpolyu_t A, slong length, const nmod_mpoly_ctx_t uctx)
{
    slong i;
    slong old_alloc = A->alloc;
    slong new_alloc = FLINT_MAX(length, 2*A->alloc);

    if (length > old_alloc)
    {
        if (old_alloc == 0)
        {
            A->exps = (ulong *) flint_malloc(new_alloc*sizeof(ulong));
            A->coeffs = (nmod_mpoly_struct *) flint_malloc(
                                          new_alloc*sizeof(nmod_mpoly_struct));
        } else
        {
            A->exps = (ulong *) flint_realloc(A->exps,
                                                      new_alloc*sizeof(ulong));
            A->coeffs = (nmod_mpoly_struct *) flint_realloc(A->coeffs,
                                          new_alloc*sizeof(nmod_mpoly_struct));
        }

        for (i = old_alloc; i < new_alloc; i++)
        {
            nmod_mpoly_init(A->coeffs + i, uctx);
            nmod_mpoly_fit_bits(A->coeffs + i, A->bits, uctx);
            (A->coeffs + i)->bits = A->bits;
        }
        A->alloc = new_alloc;
    }
}

void nmod_mpolyu_set(nmod_mpolyu_t A, const nmod_mpolyu_t B, const nmod_mpoly_ctx_t uctx)
{
    slong i;
    nmod_mpoly_struct * Acoeff, * Bcoeff;
    ulong * Aexp, * Bexp;
    slong Alen, Blen;

    Alen = 0;
    Blen = B->length;
    nmod_mpolyu_fit_length(A, Blen, uctx);
    Acoeff = A->coeffs;
    Bcoeff = B->coeffs;
    Aexp = A->exps;
    Bexp = B->exps;

    for (i = 0; i < Blen; i++)
    {
        nmod_mpoly_set(Acoeff + Alen, Bcoeff + i, uctx);
        Aexp[Alen++] = Bexp[i];
    }
    Alen = Blen;

    /* demote remaining coefficients */
    for (i = Alen; i < A->length; i++)
    {
        nmod_mpoly_clear(Acoeff + i, uctx);
        nmod_mpoly_init(Acoeff + i, uctx);
    }
    A->length = Alen;
}

void nmod_mpolyu_msub(
    nmod_mpolyu_t R,
    nmod_mpolyu_t A,
    nmod_mpolyu_t B,
    nmod_mpoly_t c,
    slong e,
    const nmod_mpoly_ctx_t ctx)
{
    nmod_mpolyu_fit_length(R, A->length + B->length, ctx);
    nmod_mpoly_t T;
    nmod_mpoly_init(T, ctx);

    slong i = 0, j = 0, k = 0;
    while (i < A->length || j < B->length)
    {
        if (i < A->length && (j >= B->length || A->exps[i] > B->exps[j] + e))
        {
            /* only A ok */
            nmod_mpoly_set(R->coeffs + k, A->coeffs + i, ctx);
            R->exps[k] = A->exps[i];
            k++;
            i++;
        }
        else if (j < B->length && (i >= A->length || B->exps[j] + e > A->exps[i]))
        {
            /* only B ok */
            nmod_mpoly_mul(R->coeffs + k, B->coeffs + j, c, ctx);
            nmod_mpoly_neg(R->coeffs + k, R->coeffs + k, ctx);
            R->exps[k] = B->exps[j] + e;
            k++;
            j++;
        }
        else if (i < A->length && j < B->length && (A->exps[i] == B->exps[j] + e))
        {
            nmod_mpoly_mul(T, B->coeffs + j, c, ctx);

//printf("T: "); nmod_mpoly_print_pretty(T, NULL, ctx); printf("\n");
//printf("A(%d): ",i); nmod_mpoly_print_pretty(A->coeffs + i, NULL, ctx); printf("\n");

            nmod_mpoly_sub(R->coeffs + k, A->coeffs + i, T, ctx);

//printf("R(%d): ",k); nmod_mpoly_print_pretty(R->coeffs + k, NULL, ctx); printf("\n");


            R->exps[k] = A->exps[i];
            k += !nmod_mpoly_is_zero(R->coeffs + k, ctx);
            i++;
            j++;
        } else 
        {
            FLINT_ASSERT(0);
        }
    }

    nmod_mpoly_clear(T, ctx);
    R->length = k;
}

int nmod_mpolyu_divides(nmod_mpolyu_t A, nmod_mpolyu_t B, const nmod_mpoly_ctx_t ctx)
{
    int ret = 0;
    nmod_mpolyu_t P, R;
    nmod_mpoly_t t;
    nmod_mpoly_init(t, ctx);
    nmod_mpolyu_init(P, A->bits, ctx);
    nmod_mpolyu_init(R, A->bits, ctx);
    nmod_mpolyu_set(R, A, ctx);

//flint_printf("testing divisiblity:\n");
//printf("A: "); nmod_mpolyu_print_pretty(A, NULL, ctx); printf("\n");
//printf("B: "); nmod_mpolyu_print_pretty(B, NULL, ctx); printf("\n");

    FLINT_ASSERT(B->length > 0);

    while (R->length > 0)
    {
        if (R->exps[0] < B->exps[0])
            goto done;

        if (!nmod_mpoly_divides(t, R->coeffs + 0, B->coeffs + 0, ctx))
            goto done;

//printf("t: "); nmod_mpoly_print_pretty(t, NULL, ctx); printf("\n");


        nmod_mpolyu_msub(P, R, B, t, R->exps[0] - B->exps[0], ctx);
        nmod_mpolyu_swap(P, R, ctx);

//printf("R: "); nmod_mpolyu_print_pretty(R, NULL, ctx); printf("\n");

    }
    ret = 1;

done:
    nmod_mpoly_clear(t, ctx);
    nmod_mpolyu_clear(P, ctx);
    nmod_mpolyu_clear(R, ctx);

//flint_printf("returning divisiblity: %d\n", ret);

    return ret;
}





/*********** mpolyn **************************/

void nmod_mpolyn_init(nmod_mpolyn_t A, mp_bitcnt_t bits, const nmod_mpoly_ctx_t ctx)
{
    A->coeffs = NULL;
    A->exps = NULL;
    A->alloc = 0;
    A->length = 0;
    A->bits = bits;
}

void nmod_mpolyn_clear(nmod_mpolyn_t A, const nmod_mpoly_ctx_t ctx)
{
    slong i;
    for (i = 0; i < A->alloc; i++)
        nmod_poly_clear(A->coeffs + i);
    flint_free(A->coeffs);
    flint_free(A->exps);
}

inline void nmod_mpolyn_swap(nmod_mpolyn_t A, nmod_mpolyn_t B)
{
   nmod_mpolyn_struct t = *A;
   *A = *B;
   *B = t;
}

void nmod_mpolyn_zero(nmod_mpolyn_t A, const nmod_mpoly_ctx_t ctx)
{
    slong i;
    for (i = 0; i < A->alloc; i++) {
        nmod_poly_clear(A->coeffs + i);
        nmod_poly_init(A->coeffs + i, ctx->ffinfo->mod.n);
    }
    A->length = 0;
}

int nmod_mpolyn_is_zero(nmod_mpolyn_t A, const nmod_mpoly_ctx_t ctx)
{
    return A->length == 0;
}

void nmod_mpolyn_print_pretty(const nmod_mpolyn_t A,
                                   const char ** x_in, const nmod_mpoly_ctx_t ctx)
{
    nmod_poly_struct * coeff = A->coeffs;
    slong len = A->length;
    ulong * exp = A->exps;
    slong bits = A->bits;
    slong i, j, N;
    fmpz * exponents;
    char ** x = (char **) x_in;
    TMP_INIT;

    if (len == 0)
    {
        flint_printf("0");
        return;
    }

    N = mpoly_words_per_exp(bits, ctx->minfo);

    TMP_START;

    if (x == NULL)
    {
        x = (char **) TMP_ALLOC(ctx->minfo->nvars*sizeof(char *));
        for (i = 0; i < ctx->minfo->nvars; i++)
        {
            x[i] = (char *) TMP_ALLOC(((FLINT_BITS+4)/3)*sizeof(char));
            flint_sprintf(x[i], "x%wd", i);
        }
    }

    exponents = (fmpz *) TMP_ALLOC(ctx->minfo->nvars*sizeof(ulong));
    for (i = 0; i < ctx->minfo->nvars; i++)
        fmpz_init(exponents + i);
   
    for (i = 0; i < len; i++)
    {
        if (i > 0)
        {
            printf(" + ");
        }

        printf("(");
        nmod_poly_print_pretty(coeff + i, "v");
        printf(")");

        mpoly_get_monomial_ffmpz(exponents, exp + N*i, bits, ctx->minfo);

        for (j = 0; j < ctx->minfo->nvars; j++)
        {
            int cmp = fmpz_cmp_ui(exponents + j, WORD(1));

            if (cmp > 0)
            {
                printf("*%s^", x[j]);
                fmpz_print(exponents + j);
            } else if (cmp == 0)
            {
                printf("*%s", x[j]);
            }
        }
    }

    for (i = 0; i < ctx->minfo->nvars; i++)
        fmpz_clear(exponents + i);

    TMP_END;
}

void nmod_mpolyn_fit_length(nmod_mpolyn_t A, slong length, const nmod_mpoly_ctx_t ctx)
{
    slong i;
    slong old_alloc = A->alloc;
    slong new_alloc = FLINT_MAX(length, 2*A->alloc);

    if (length > old_alloc)
    {
        slong N = mpoly_words_per_exp(A->bits, ctx->minfo);

        if (old_alloc == 0)
        {
            A->exps = (ulong *) flint_malloc(new_alloc*N*sizeof(ulong));
            A->coeffs = (nmod_poly_struct *) flint_malloc(new_alloc*sizeof(nmod_poly_struct));
        } else
        {
            A->exps = (ulong *) flint_realloc(A->exps, new_alloc*N*sizeof(ulong));
            A->coeffs = (nmod_poly_struct *) flint_realloc(A->coeffs, new_alloc*sizeof(nmod_poly_struct));
        }

        for (i = old_alloc; i < new_alloc; i++)
        {
            nmod_poly_init(A->coeffs + i, ctx->ffinfo->mod.n);
        }
        A->alloc = new_alloc;
    }
}

void nmod_mpolyn_set_length(nmod_mpolyn_t A, slong newlen, const nmod_mpoly_ctx_t ctx)
{
    if (A->length > newlen)
    {
        slong i;
        for (i = newlen; i < A->length; i++)
        {
            nmod_poly_clear(A->coeffs + i);
            nmod_poly_init(A->coeffs + i, ctx->ffinfo->mod.n);
        }
    }
    A->length = newlen;
}

void nmod_mpolyn_fit_bits(nmod_mpolyn_t A, slong bits, const nmod_mpoly_ctx_t ctx)
{
   slong N;
   ulong * t;

   if (A->bits < bits)
   {
      if (A->alloc != 0)
      {
         N = mpoly_words_per_exp(bits, ctx->minfo);
         t = flint_malloc(N*A->alloc*sizeof(ulong));
         mpoly_repack_monomials(t, bits, A->exps, A->bits, A->length, ctx->minfo);
         flint_free(A->exps);
         A->exps = t;
      }

      A->bits = bits;
   }
}

void nmod_mpolyn_set(nmod_mpolyn_t A, const nmod_mpolyn_t B, const nmod_mpoly_ctx_t ctx)
{
    slong i;
    nmod_poly_struct * Acoeff, * Bcoeff;
    ulong * Aexp, * Bexp;
    slong Blen;

    nmod_mpolyn_fit_bits(A, B->bits, ctx);
    A->bits = B->bits;

    Blen = B->length;
    nmod_mpolyn_fit_length(A, Blen, ctx);
    Acoeff = A->coeffs;
    Bcoeff = B->coeffs;
    Aexp = A->exps;
    Bexp = B->exps;

    slong N = mpoly_words_per_exp(B->bits, ctx->minfo);

    for (i = 0; i < Blen; i++)
    {
        nmod_poly_set(Acoeff + i, Bcoeff + i);
        mpoly_monomial_set(Aexp + N*i, Bexp + N*i, N);
    }

    /* demote remaining coefficients */
    for (i = Blen; i < A->length; i++)
    {
        nmod_poly_clear(Acoeff + i);
        nmod_poly_init(Acoeff + i, ctx->ffinfo->mod.n);
    }
    A->length = Blen;
}

void nmod_mpolyn_set_mpoly(nmod_mpolyn_t A, const nmod_mpoly_t B, const nmod_mpoly_ctx_t ctx)
{
    slong i;
    nmod_poly_struct * Acoeff;
    mp_limb_t * Bcoeff;
    ulong * Aexp, * Bexp;
    slong Blen;

    nmod_mpolyn_fit_bits(A, B->bits, ctx);
    A->bits = B->bits;

    Blen = B->length;
    nmod_mpolyn_fit_length(A, Blen, ctx);
    Acoeff = A->coeffs;
    Bcoeff = B->coeffs;
    Aexp = A->exps;
    Bexp = B->exps;

    slong N = mpoly_words_per_exp(B->bits, ctx->minfo);

    for (i = 0; i < Blen; i++)
    {
        nmod_poly_zero(Acoeff + i);
        nmod_poly_set_coeff_ui(Acoeff + i, 0, Bcoeff[i]);
        mpoly_monomial_set(Aexp + N*i, Bexp + N*i, N);
    }

    /* demote remaining coefficients */
    for (i = Blen; i < A->length; i++)
    {
        nmod_poly_clear(Acoeff + i);
        nmod_poly_init(Acoeff + i, ctx->ffinfo->mod.n);
    }
    A->length = Blen;
}



void nmod_mpolyn_mul_poly(
    nmod_mpolyn_t A,
    const nmod_mpolyn_t B,
    const nmod_poly_t c,
    const nmod_mpoly_ctx_t ctx)
{
    slong i;
    nmod_poly_struct * Acoeff, * Bcoeff;
    ulong * Aexp, * Bexp;
    slong Blen;

    nmod_mpolyn_fit_bits(A, B->bits, ctx);
    A->bits = B->bits;

    Blen = B->length;
    nmod_mpolyn_fit_length(A, Blen, ctx);
    Acoeff = A->coeffs;
    Bcoeff = B->coeffs;
    Aexp = A->exps;
    Bexp = B->exps;

    slong N = mpoly_words_per_exp(B->bits, ctx->minfo);

    FLINT_ASSERT(!nmod_poly_is_zero(c));

    for (i = 0; i < Blen; i++)
    {
        nmod_poly_mul(Acoeff + i, Bcoeff + i, c);
        mpoly_monomial_set(Aexp + N*i, Bexp + N*i, N);
    }

    /* demote remaining coefficients */
    for (i = Blen; i < A->length; i++)
    {
        nmod_poly_clear(Acoeff + i);
        nmod_poly_init(Acoeff + i, ctx->ffinfo->mod.n);
    }
    A->length = Blen;
}


/*********** mpolyun **************************/

void nmod_mpolyun_init(nmod_mpolyun_t A, mp_bitcnt_t bits, const nmod_mpoly_ctx_t ctx)
{
    A->coeffs = NULL;
    A->exps = NULL;
    A->alloc = 0;
    A->length = 0;
    A->bits = bits;
}


void nmod_mpolyun_clear(nmod_mpolyun_t A, const nmod_mpoly_ctx_t ctx)
{
    slong i;
    for (i = 0; i < A->alloc; i++)
        nmod_mpolyn_clear(A->coeffs + i, ctx);
    flint_free(A->coeffs);
    flint_free(A->exps);
}


inline void nmod_mpolyun_swap(nmod_mpolyun_t A, nmod_mpolyun_t B)
{
   nmod_mpolyun_struct t = *A;
   *A = *B;
   *B = t;
}

void nmod_mpolyun_zero(nmod_mpolyun_t A, const nmod_mpoly_ctx_t ctx)
{
    slong i;
    for (i = 0; i < A->alloc; i++) {
        nmod_mpolyn_clear(A->coeffs + i, ctx);
        nmod_mpolyn_init(A->coeffs + i, A->bits, ctx);
    }
    A->length = 0;
}

void nmod_mpolyun_print_pretty(const nmod_mpolyun_t poly,
                                   const char ** x, const nmod_mpoly_ctx_t ctx)
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
        nmod_mpolyn_print_pretty(poly->coeffs + i,x,ctx);
        flint_printf(")*X^%wd",poly->exps[i]);
    }
}

void nmod_mpolyun_fit_length(nmod_mpolyun_t A, slong length, const nmod_mpoly_ctx_t ctx)
{
    slong i;
    slong old_alloc = A->alloc;
    slong new_alloc = FLINT_MAX(length, 2*A->alloc);

    if (length > old_alloc)
    {
        if (old_alloc == 0)
        {
            A->exps = (ulong *) flint_malloc(new_alloc*sizeof(ulong));
            A->coeffs = (nmod_mpolyn_struct *) flint_malloc(
                                          new_alloc*sizeof(nmod_mpolyn_struct));
        } else
        {
            A->exps = (ulong *) flint_realloc(A->exps,
                                                      new_alloc*sizeof(ulong));
            A->coeffs = (nmod_mpolyn_struct *) flint_realloc(A->coeffs,
                                          new_alloc*sizeof(nmod_mpolyn_struct));
        }

        for (i = old_alloc; i < new_alloc; i++)
        {
            nmod_mpolyn_init(A->coeffs + i, A->bits, ctx);
        }
        A->alloc = new_alloc;
    }
}

void nmod_mpolyun_set(nmod_mpolyun_t A, const nmod_mpolyun_t B, const nmod_mpoly_ctx_t ctx)
{
    slong i, Blen;
    nmod_mpolyn_struct * Acoeff, * Bcoeff;
    ulong * Aexp, * Bexp;

    Blen = B->length;
    nmod_mpolyun_fit_length(A, Blen, ctx);
    Acoeff = A->coeffs;
    Bcoeff = B->coeffs;
    Aexp = A->exps;
    Bexp = B->exps;

    for (i = 0; i < Blen; i++)
    {
        nmod_mpolyn_set(Acoeff + i, Bcoeff + i, ctx);
        Aexp[i] = Bexp[i];
    }

    /* demote remaining coefficients */
    for (i = Blen; i < A->length; i++)
    {
        nmod_mpolyn_clear(Acoeff + i, ctx);
        nmod_mpolyn_init(Acoeff + i, A->bits, ctx);
    }
    A->length = Blen;
}


void nmod_mpolyun_mul_poly(
    nmod_mpolyun_t A,
    const nmod_mpolyun_t B,
    const nmod_poly_t c,
    const nmod_mpoly_ctx_t ctx)
{
    slong i, Blen;
    nmod_mpolyn_struct * Acoeff, * Bcoeff;
    ulong * Aexp, * Bexp;

    Blen = B->length;
    nmod_mpolyun_fit_length(A, Blen, ctx);
    Acoeff = A->coeffs;
    Bcoeff = B->coeffs;
    Aexp = A->exps;
    Bexp = B->exps;

    for (i = 0; i < Blen; i++)
    {
        nmod_mpolyn_mul_poly(Acoeff + i, Bcoeff + i, c, ctx);
        Aexp[i] = Bexp[i];
    }

    /* demote remaining coefficients */
    for (i = Blen; i < A->length; i++)
    {
        nmod_mpolyn_clear(Acoeff + i, ctx);
        nmod_mpolyn_init(Acoeff + i, A->bits, ctx);
    }
    A->length = Blen;
}


nmod_poly_struct * nmod_mpolyn_leadcoeff_ref(nmod_mpolyn_t A, nmod_mpoly_ctx_t ctx)
{
    FLINT_ASSERT(A->length > 0);
    return A->coeffs + 0;
}

nmod_poly_struct * nmod_mpolyun_leadcoeff_ref(nmod_mpolyun_t A, nmod_mpoly_ctx_t ctx)
{
    FLINT_ASSERT(A->length > 0);
    return nmod_mpolyn_leadcoeff_ref(A->coeffs + 0, ctx);
}


mp_limb_t nmod_mpoly_leadcoeff(nmod_mpoly_t A, nmod_mpoly_ctx_t ctx)
{
    FLINT_ASSERT(A->length > 0);
    return A->coeffs[0];
}

mp_limb_t nmod_mpolyu_leadcoeff(nmod_mpolyu_t A, nmod_mpoly_ctx_t ctx)
{
    FLINT_ASSERT(A->length > 0);
    return nmod_mpoly_leadcoeff(A->coeffs + 0, ctx);
}






/* if the coefficient doesn't exist, a new one is created (and set to zero) */
nmod_mpoly_struct * _nmod_mpolyu_get_coeff(nmod_mpolyu_t A,
                             ulong pow, const nmod_mpoly_ctx_t uctx)
{
    slong i, j;
    nmod_mpoly_struct * xk;

    for (i = 0; i < A->length && A->exps[i] >= pow; i++)
    {
//flint_printf(" get coeff i = %wd\n",i);
        if (A->exps[i] == pow) 
        {
            return A->coeffs + i;
        }
    }

    nmod_mpolyu_fit_length(A, A->length + 1, uctx);

//flint_printf("bits = %wd\n",bits);


    for (j = A->length; j > i; j--)
    {
        A->exps[j] = A->exps[j - 1];
        nmod_mpoly_swap(A->coeffs + j, A->coeffs + j - 1, uctx);
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


/*
    emplaceterm assumes that c is valid modulo ctx->ffinfo->mod.n
*/
void _nmod_mpoly_emplacebackterm_ui_ffmpz(nmod_mpoly_t poly,
                    mp_limb_t c, const fmpz * exp, const nmod_mpoly_ctx_t ctx)
{
    mp_bitcnt_t exp_bits;
    slong N;
    FLINT_ASSERT(c < ctx->ffinfo->mod.n);

    exp_bits = mpoly_exp_bits_required_ffmpz(exp, ctx->minfo);
    exp_bits = mpoly_fix_bits(exp_bits, ctx->minfo);
    nmod_mpoly_fit_bits(poly, exp_bits, ctx);

    N = mpoly_words_per_exp(poly->bits, ctx->minfo);

    nmod_mpoly_fit_length(poly, poly->length + 1, ctx);
    poly->coeffs[poly->length] = c;
    mpoly_set_monomial_ffmpz(poly->exps + N*poly->length, exp, poly->bits, ctx->minfo);
    poly->length++; /* safe because length is increasing */
}

void nmod_mpoly_to_mpolyu_perm(nmod_mpolyu_t A, const nmod_mpoly_t B,
                                                    const slong * perm,
                                                    const nmod_mpoly_ctx_t uctx,
                                                    const nmod_mpoly_ctx_t ctx)
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
    nmod_mpolyu_zero(A, uctx);

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
        nmod_mpoly_struct * Ac = _nmod_mpolyu_get_coeff(A, uexps[n - 1], uctx);

        FLINT_ASSERT(Ac->bits == B->bits);


//flint_printf("        uexp:");for (i=0; i<n; i++){flint_printf(", %wd", uexps[i]);}printf("\n");

        nmod_mpoly_fit_length(Ac, Ac->length + 1, uctx);
        Ac->coeffs[Ac->length] = B->coeffs[j];
        mpoly_set_monomial_ffmpz(Ac->exps + NA*Ac->length, uexps, A->bits, uctx->minfo);
        Ac->length++;

//flint_printf("A after: "); nmod_mpolyu_print_pretty(A, NULL, uctx); printf("\n");

    }



    for (i = 0; i < A->length; i++)
    {
        nmod_mpoly_sort_terms(A->coeffs + i, uctx);
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
void nmod_mpoly_from_mpolyu_perm(nmod_mpoly_t poly1, const nmod_mpolyu_t poly2,
                                                    const slong * perm,
                                                    const nmod_mpoly_ctx_t uctx,
                                                    const nmod_mpoly_ctx_t ctx)
{
    slong i, j, k, bits, N;
    slong p_len;
    mp_limb_t * p_coeff;
    ulong * p_exp;
    slong p_alloc;
    slong n = ctx->minfo->nvars;
    TMP_INIT;

/*
flint_printf("nmod_mpoly_from_mpolyu_perm called\n");
nmod_mpolyu_print_pretty(poly2, NULL, ctx); printf("\n");
flint_printf("ctx->minfo->ord: %d\n", ctx->minfo->ord);
flint_printf("n: %wd\n", n);
*/

    FLINT_ASSERT(uctx->minfo->nvars == n - 1);

    if (poly2->length == 0)
    {
        nmod_mpoly_zero(poly1, ctx);
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
    for (i = 0; i < poly2->length; i++)
    {
        nmod_mpoly_struct * pi = poly2->coeffs + i;
        mpoly_degrees_ffmpz(texps, pi->exps, pi->length, pi->bits, uctx->minfo);
        _fmpz_vec_max_inplace(uexps, texps, n - 1);
    }
    fmpz_set_ui(uexps + n - 1, poly2->exps[0]);
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

    nmod_mpoly_fit_bits(poly1, bits, ctx);
    poly1->bits = bits;

    p_coeff = poly1->coeffs;
    p_exp = poly1->exps;
    p_alloc = poly1->alloc;
    p_len = 0;
    for (i = 0; i < poly2->length; i++)
    {
        nmod_mpoly_struct * Bc = poly2->coeffs + i;
        _nmod_mpoly_fit_length(&p_coeff, &p_exp, &p_alloc, p_len + Bc->length, N);
        slong NB = mpoly_words_per_exp(Bc->bits, uctx->minfo);
        for (j = 0; j < Bc->length; j++)
        {
            p_coeff[p_len + j] = Bc->coeffs[j];
            mpoly_get_monomial_ffmpz(uexps, Bc->exps + NB*j, Bc->bits, uctx->minfo);
            fmpz_set_ui(uexps + n - 1, poly2->exps[i]);
/*
flint_printf("(%wd, %wd) uexp",i,j);
for (k=0; k<n; k++){
flint_printf(", %wd", uexps[k]);
}
printf("\n");
*/

            for (k = 0; k < n; k++)
            {
                fmpz_swap(uexps + k, exps + perm[k]);
            }
/*
flint_printf("        exp",i,j);
for (k=0; k<n; k++){
flint_printf(", %wd", exps[k]);
}
printf("\n");
*/

            mpoly_set_monomial_ffmpz(p_exp + N*(p_len + j), exps, bits, ctx->minfo);
/*
printf("packed exp: %016llx\n", (p_exp + N*(p_len + j))[0]);
*/
        }
        p_len += (poly2->coeffs + i)->length;
    }
    poly1->coeffs = p_coeff;
    poly1->exps = p_exp;
    poly1->alloc = p_alloc;
    _nmod_mpoly_set_length(poly1, p_len, ctx);
/*
printf("poly: "); nmod_mpoly_print_pretty(poly1, NULL, ctx); printf("\n");
*/
    for (i = 0; i < n; i++)
    {
        fmpz_clear(texps + i);
        fmpz_clear(uexps + i);
        fmpz_clear(exps + i);
    }

    nmod_mpoly_sort_terms(poly1, ctx);
    TMP_END;
}



void nmod_mpoly_from_mpolyu_keepbits(nmod_mpoly_t poly1, const nmod_mpolyu_t poly2,
                                                    const slong * perm,
                                                    const nmod_mpoly_ctx_t uctx,
                                                    const nmod_mpoly_ctx_t ctx)
{
    slong i, j, k, bits, N;
    slong p_len;
    mp_limb_t * p_coeff;
    ulong * p_exp;
    slong p_alloc;
    slong n = ctx->minfo->nvars;
    TMP_INIT;

/*
flint_printf("nmod_mpoly_from_mpolyu_perm called\n");
nmod_mpolyu_print_pretty(poly2, NULL, ctx); printf("\n");
flint_printf("ctx->minfo->ord: %d\n", ctx->minfo->ord);
flint_printf("n: %wd\n", n);
*/

    FLINT_ASSERT(uctx->minfo->nvars == n - 1);

    if (poly2->length == 0)
    {
        nmod_mpoly_zero(poly1, ctx);
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
    for (i = 0; i < poly2->length; i++)
    {
        nmod_mpoly_struct * pi = poly2->coeffs + i;
        mpoly_degrees_ffmpz(texps, pi->exps, pi->length, pi->bits, uctx->minfo);
        _fmpz_vec_max_inplace(uexps, texps, n - 1);
    }
    fmpz_set_ui(uexps + n - 1, poly2->exps[0]);
    for (i = 0; i < n; i++)
    {
        fmpz_swap(uexps + i, exps + perm[i]);
    }
    bits = mpoly_exp_bits_required_ffmpz(exps, ctx->minfo);

    FLINT_ASSERT(bits <= poly2->bits);
    bits = poly2->bits;
/*
flint_printf("bits: %wd\n", bits);
*/
    N = mpoly_words_per_exp(bits, ctx->minfo);

    nmod_mpoly_fit_bits(poly1, bits, ctx);
    poly1->bits = bits;

    p_coeff = poly1->coeffs;
    p_exp = poly1->exps;
    p_alloc = poly1->alloc;
    p_len = 0;
    for (i = 0; i < poly2->length; i++)
    {
        nmod_mpoly_struct * Bc = poly2->coeffs + i;
        _nmod_mpoly_fit_length(&p_coeff, &p_exp, &p_alloc, p_len + Bc->length, N);
        slong NB = mpoly_words_per_exp(Bc->bits, uctx->minfo);
        for (j = 0; j < Bc->length; j++)
        {
            p_coeff[p_len + j] = Bc->coeffs[j];
            mpoly_get_monomial_ffmpz(uexps, Bc->exps + NB*j, Bc->bits, uctx->minfo);
            fmpz_set_ui(uexps + n - 1, poly2->exps[i]);
/*
flint_printf("(%wd, %wd) uexp",i,j);
for (k=0; k<n; k++){
flint_printf(", %wd", uexps[k]);
}
printf("\n");
*/

            for (k = 0; k < n; k++)
            {
                fmpz_swap(uexps + k, exps + perm[k]);
            }
/*
flint_printf("        exp",i,j);
for (k=0; k<n; k++){
flint_printf(", %wd", exps[k]);
}
printf("\n");
*/

            mpoly_set_monomial_ffmpz(p_exp + N*(p_len + j), exps, bits, ctx->minfo);
/*
printf("packed exp: %016llx\n", (p_exp + N*(p_len + j))[0]);
*/
        }
        p_len += (poly2->coeffs + i)->length;
    }
    poly1->coeffs = p_coeff;
    poly1->exps = p_exp;
    poly1->alloc = p_alloc;
    _nmod_mpoly_set_length(poly1, p_len, ctx);
/*
printf("poly: "); nmod_mpoly_print_pretty(poly1, NULL, ctx); printf("\n");
*/
    for (i = 0; i < n; i++)
    {
        fmpz_clear(texps + i);
        fmpz_clear(uexps + i);
        fmpz_clear(exps + i);
    }

    nmod_mpoly_sort_terms(poly1, ctx);
    TMP_END;
}




/*
    F = F + modulus*(A - F(alpha))
*/
int
nmod_mpolyn_addinterp(
    slong * lastdeg,
    nmod_mpolyn_t F,
    nmod_mpolyn_t T,
    nmod_mpoly_t A,
    nmod_poly_t modulus,
    mp_limb_t alpha,
    const nmod_mpoly_ctx_t ctx)
{
    int changed = 0;
    slong i, j, k;
    mp_limb_t v;
    mp_bitcnt_t bits = A->bits;
    slong Flen = F->length, Alen = A->length;
    ulong * Fexp = F->exps, * Aexp = A->exps;
    mp_limb_t * Acoeff = A->coeffs;
    nmod_poly_struct * Fcoeff = F->coeffs;
    nmod_poly_t tp;
/*
flint_printf("********** nmod_mpolyn_addinterp:\n");
flint_printf("F(%d): ",F->bits); nmod_mpolyn_print_pretty(F, NULL, ctx); printf("\n");
flint_printf("A(%d): ",A->bits); nmod_mpoly_print_pretty(A, NULL, ctx); printf("\n");
flint_printf("modulus: "); nmod_poly_print_pretty(modulus, "v"); printf("\n");
flint_printf("alpha: %d\n", alpha);
*/
    FLINT_ASSERT(F->bits == bits);
    FLINT_ASSERT(ctx->minfo->ord == ORD_LEX);

    nmod_poly_init(tp, ctx->ffinfo->mod.n);

    nmod_mpolyn_fit_length(T, Flen + Alen, ctx);
    ulong * Texp = T->exps;
    nmod_poly_struct * Tcoeff = T->coeffs;

    slong N = mpoly_words_per_exp(bits, ctx->minfo);

    i = j = k = 0;
    while (i < Flen || j < Alen)
    {
        if (i < Flen && (j >= Alen || mpoly_monomial_gt_nomask(Fexp + N*i, Aexp + N*j, N)))
        {
            FLINT_ASSERT(!nmod_poly_is_zero(Fcoeff + i));
            FLINT_ASSERT(nmod_poly_degree(Fcoeff + i) < nmod_poly_degree(modulus));
            /* F term ok, A term missing */
            v = nmod_poly_evaluate_nmod(Fcoeff + i, alpha);
            if (v != UWORD(0))
            {
                changed = 1;
                nmod_poly_scalar_mul_nmod(tp, modulus, v);
                nmod_poly_sub(Tcoeff + k, Fcoeff + i, tp);
            } else {
                nmod_poly_set(Tcoeff + k, Fcoeff + i);                
            }
            lastdeg[0] = FLINT_MAX(lastdeg[0], nmod_poly_degree(Tcoeff + k));

            mpoly_monomial_set(Texp + N*k, Fexp + N*i, N);
            FLINT_ASSERT(!nmod_poly_is_zero(Tcoeff + k));
            k++;
            i++;
        }
        else if (j < Alen && (i >= Flen || mpoly_monomial_lt_nomask(Fexp + N*i, Aexp + N*j, N)))
        {
            /* F term missing, A term ok */
            if (Acoeff[j] != UWORD(0))
            {
                changed = 1;
                nmod_poly_zero(Tcoeff + k);
                nmod_poly_scalar_mul_nmod(Tcoeff + k, modulus, Acoeff[j]);
                lastdeg[0] = FLINT_MAX(lastdeg[0], nmod_poly_degree(Tcoeff + k));
                mpoly_monomial_set(Texp + N*k, Aexp + N*j, N);
                k++;
            }
            j++;
        }
        else if (i < Flen && j < Alen && mpoly_monomial_equal(Fexp + N*i, Aexp + N*j, N))
        {
            FLINT_ASSERT(!nmod_poly_is_zero(Fcoeff + i));
            FLINT_ASSERT(nmod_poly_degree(Fcoeff + i) < nmod_poly_degree(modulus));

            /* F term ok, A term ok */
            v = nmod_poly_evaluate_nmod(Fcoeff + i, alpha);
            v = nmod_sub(Acoeff[j], v, ctx->ffinfo->mod);
            if (v != UWORD(0))
            {
                changed = 1;
                nmod_poly_scalar_mul_nmod(tp, modulus, v);
                nmod_poly_add(Tcoeff + k, Fcoeff + i, tp);
            } else {
                nmod_poly_set(Tcoeff + k, Fcoeff + i);                
            }
            lastdeg[0] = FLINT_MAX(lastdeg[0], nmod_poly_degree(Tcoeff + k));
            mpoly_monomial_set(Texp + N*k, Aexp + N*j, N);
            FLINT_ASSERT(!nmod_poly_is_zero(Tcoeff + k));
            k++;
            i++;
            j++;
        } else 
        {
            FLINT_ASSERT(0);
        }
    }

    nmod_mpolyn_set_length(T, k, ctx);

    if (changed)
    {
        nmod_mpolyn_swap(T, F);
    }

    nmod_poly_clear(tp);
    return changed;
}


int nmod_mpolyun_addinterp(
    slong * lastdeg,
    nmod_mpolyun_t F,
    nmod_mpolyun_t T,
    nmod_mpolyu_t A,
    nmod_poly_t modulus,
    mp_limb_t alpha,
    const nmod_mpoly_ctx_t ctx)
{
    int changed = 0;
    slong i, j, k;
/*
flint_printf("nmod_mpolyun_addinterp:\n");
flint_printf("F(%d): ",F->bits); nmod_mpolyun_print_pretty(F, NULL, ctx); printf("\n");
flint_printf("A(%d): ",A->bits); nmod_mpolyu_print_pretty(A, NULL, ctx); printf("\n");
flint_printf("modulus: "); nmod_poly_print_pretty(modulus, "v"); printf("\n");
flint_printf("alpha: %d\n", alpha);
*/
    lastdeg[0] = -WORD(1);

    FLINT_ASSERT(F->bits == T->bits);
    FLINT_ASSERT(T->bits == A->bits);

    nmod_mpolyn_t S;
    nmod_mpolyn_init(S, F->bits, ctx);

    slong Flen = F->length;
    slong Alen = A->length;
    nmod_mpolyun_fit_length(T, Flen + Alen, ctx);

    nmod_mpolyn_struct * Tcoeff = T->coeffs;
    nmod_mpolyn_struct * Fcoeff = F->coeffs;
    nmod_mpoly_struct  * Acoeff = A->coeffs;
    ulong * Texp = T->exps;
    ulong * Fexp = F->exps;
    ulong * Aexp = A->exps;   

    nmod_mpoly_t zero;
    nmod_mpoly_init(zero, ctx);
    nmod_mpoly_fit_bits(zero, A->bits, ctx);
    zero->bits = A->bits;

    i = j = k = 0;
    while (i < Flen || j < Alen)
    {
//flint_printf("i: %d  j: %d\n",i,j);

        if (i < Flen && (j >= Alen || Fexp[i] > Aexp[j]))
        {
//flint_printf("case 1  i: %d  j: %d\n",i,j);

            /* F term ok, A term missing */
            nmod_mpolyn_set(Tcoeff + k, Fcoeff + i, ctx);
            changed |= nmod_mpolyn_addinterp(lastdeg, Tcoeff + k, S, zero, modulus, alpha, ctx);
            Texp[k] = Fexp[i];
            k++;
            i++;
        }
        else if (j < Alen && (i >= Flen || Aexp[j] > Fexp[i]))
        {
//flint_printf("case 2  i: %d  j: %d\n",i,j);

            /* F term missing, A term ok */
            nmod_mpolyn_zero(Tcoeff + k, ctx);
            changed |= nmod_mpolyn_addinterp(lastdeg, Tcoeff + k, S, Acoeff + j, modulus, alpha, ctx);
            Texp[k] = Aexp[j];
            k++;
            j++;
        }
        else if (i < Flen && j < Alen && (Fexp[i] == Aexp[j]))
        {
//flint_printf("case 3  i: %d  j: %d\n",i,j);

            /* F term ok, A term ok */
            nmod_mpolyn_set(Tcoeff + k, Fcoeff + i, ctx);
            changed |= nmod_mpolyn_addinterp(lastdeg, Tcoeff + k, S, Acoeff + j, modulus, alpha, ctx);
            Texp[k] = Aexp[j];
            FLINT_ASSERT(!nmod_mpolyn_is_zero(Tcoeff + k, ctx));
            k++;
            i++;
            j++;
        } else 
        {

            FLINT_ASSERT(0);
        }
    }

//printf("3\n");

    /* demote remaining coefficients */
    for (i = k; i < T->length; i++)
    {
        nmod_mpolyn_clear(Tcoeff + i, ctx);
        nmod_mpolyn_init(Tcoeff + i, T->bits, ctx);
    }
    T->length = k;

//printf("4\n");

    if (changed)
    {
        nmod_mpolyun_swap(T, F);
    }


//flint_printf("nmod_mpolyun_addinterp done!\n");

    nmod_mpolyn_clear(S, ctx);
    nmod_mpoly_clear(zero, ctx);
    return changed;    
}



/*
void nmod_mpolyn_set_mpoly(
    nmod_mpolyn_t A,
    nmod_mpoly_t B,
    const nmod_mpoly_ctx_t ctx)
{
    nmod_mpolyn_fit_length(A, B->length, ctx);
    FLINT_ASSERT(A->bits == B->bits);

    slong N = mpoly_words_per_exp(B->bits, ctx);

    for (slong i = 0; i < B->length; i++)
    {
        nmod_poly_zero(A->coeffs + i);
        nmod_poly_set_coeff_ui(A->coeffs + i, 0, B->coeffs[i]);
        mpoly_monomial_set(A->exps + N*i, B->exps + N*i, N);
    }
}
*/

void nmod_mpolyun_set_mpolyu(
    nmod_mpolyun_t A,
    const nmod_mpolyu_t B,
    const nmod_mpoly_ctx_t ctx)
{
    FLINT_ASSERT(A->bits == B->bits);
    nmod_mpolyun_fit_length(A, B->length, ctx);
    for (slong i = 0; i < B->length; i++)
    {
        A->exps[i] = B->exps[i];
        nmod_mpolyn_set_mpoly(A->coeffs + i, B->coeffs + i, ctx);

        FLINT_ASSERT((A->coeffs + i)->bits == B->bits);

    }
    A->length = B->length;
}


void nmod_mpoly_cvtto_mpolyn(
    nmod_mpolyn_t A,
    nmod_mpoly_t B,
    slong var,
    nmod_mpoly_ctx_t ctx)
{
    TMP_INIT;

    FLINT_ASSERT(B->bits <= FLINT_BITS);
    FLINT_ASSERT(ctx->minfo->ord == ORD_LEX);

    TMP_START;

    slong N = mpoly_words_per_exp(B->bits, ctx->minfo);
    ulong * oneexp = TMP_ALLOC(N*sizeof(ulong));
    slong offset;
    slong shift;
    ulong mask = (-UWORD(1)) >> (FLINT_BITS - B->bits);
    mpoly_gen_oneexp_offset_shift(oneexp, &offset, &shift, var, N, B->bits, ctx->minfo);

    nmod_mpolyn_fit_bits(A, B->bits, ctx);
    A->bits = B->bits;

    slong k = 0;
    nmod_mpolyn_fit_length(A, k + 1, ctx);
    for (slong i = 0; i < B->length; i++)
    {
        ulong c = (B->exps[N*i + offset] >> shift) & mask;
        mpoly_monomial_msub(A->exps + N*k, B->exps + N*i, c, oneexp, N);

        if (k > 0 && mpoly_monomial_equal(A->exps + N*k, A->exps + N*(k - 1), N))
        {
            nmod_poly_set_coeff_ui(A->coeffs + k - 1, c, B->coeffs[i]);
        } else
        {
            nmod_poly_zero(A->coeffs + k);
            nmod_poly_set_coeff_ui(A->coeffs + k, c, B->coeffs[i]);
            k++;
            nmod_mpolyn_fit_length(A, k + 1, ctx);
        }
    }

    nmod_mpolyn_set_length(A, k, ctx);
    TMP_END;
}
void nmod_mpolyu_cvtto_mpolyun(
    nmod_mpolyun_t A,
    nmod_mpolyu_t B,
    slong k,
    nmod_mpoly_ctx_t ctx)
{
    slong i, Blen;
    nmod_mpolyn_struct * Acoeff;
    nmod_mpoly_struct * Bcoeff;
    ulong * Aexp, * Bexp;

    Blen = B->length;
    nmod_mpolyun_fit_length(A, Blen, ctx);
    Acoeff = A->coeffs;
    Bcoeff = B->coeffs;
    Aexp = A->exps;
    Bexp = B->exps;

    for (i = 0; i < Blen; i++)
    {
        nmod_mpoly_cvtto_mpolyn(Acoeff + i, Bcoeff + i, k, ctx);
        Aexp[i] = Bexp[i];
    }

    /* demote remaining coefficients */
    for (i = Blen; i < A->length; i++)
    {
        nmod_mpolyn_clear(Acoeff + i, ctx);
        nmod_mpolyn_init(Acoeff + i, A->bits, ctx);
    }
    A->length = Blen;  
}





void nmod_mpoly_cvtfrom_mpolyn(
    nmod_mpoly_t A,
    nmod_mpolyn_t B,
    slong var,
    nmod_mpoly_ctx_t ctx)
{
    TMP_INIT;

    FLINT_ASSERT(B->bits == A->bits);
    FLINT_ASSERT(B->bits <= FLINT_BITS);
    FLINT_ASSERT(ctx->minfo->ord == ORD_LEX);

    TMP_START;

    slong N = mpoly_words_per_exp(B->bits, ctx->minfo);
    ulong * oneexp = TMP_ALLOC(N*sizeof(ulong));
    slong offset;
    slong shift;
    mpoly_gen_oneexp_offset_shift(oneexp, &offset, &shift, var, N, B->bits, ctx->minfo);

    nmod_mpoly_fit_length(A, B->length, ctx);

    slong k = 0;
    for (slong i = 0; i < B->length; i++)
    {
        for (slong j = (B->coeffs + i)->length - 1; j >= 0; j--)
        {
            mp_limb_t c = (B->coeffs + i)->coeffs[j];
            if (c != UWORD(0))
            {
                nmod_mpoly_fit_length(A, k + 1, ctx);
                A->coeffs[k] = c;
                mpoly_monomial_madd(A->exps + N*k, B->exps + N*i, j, oneexp, N);                
                k++;
            }
        }
    }

    A->length = k;
    TMP_END;
}
void nmod_mpolyu_cvtfrom_mpolyun(
    nmod_mpolyu_t A,
    nmod_mpolyun_t B,
    slong var,
    nmod_mpoly_ctx_t ctx)
{
    nmod_mpolyu_fit_length(A, B->length, ctx);

    for (slong i = 0; i < B->length; i++)
    {
        nmod_mpoly_cvtfrom_mpolyn(A->coeffs + i, B->coeffs + i, var, ctx);
        A->exps[i] = B->exps[i];
    }
    A->length = B->length;
}




void nmod_mpolyun_content_last(nmod_poly_t a, nmod_mpolyun_t B, nmod_mpoly_ctx_t ctx)
{
    nmod_poly_zero(a);
    for (slong i = 0; i < B->length; i++)
    {
        for (slong j = 0; j < (B->coeffs + i)->length; j++)
        {
            nmod_poly_gcd(a, a, (B->coeffs + i)->coeffs + j);
        }
    }
}

void nmod_mpolyun_divexact_last(nmod_mpolyun_t A, nmod_poly_t b, nmod_mpoly_ctx_t ctx)
{
    nmod_poly_t r;
    nmod_poly_init(r, ctx->ffinfo->mod.n);
    for (slong i = 0; i < A->length; i++)
    {
        for (slong j = 0; j < (A->coeffs + i)->length; j++)
        {
            nmod_poly_divrem((A->coeffs + i)->coeffs + j, r, (A->coeffs + i)->coeffs + j, b);
            FLINT_ASSERT(nmod_poly_is_zero(r));
        }
    }
    nmod_poly_clear(r);
}


void nmod_mpoly_evalsk(
    nmod_mpoly_t A,
    nmod_mpoly_t B,
    slong entries,
    slong * offs,
    ulong * masks,
    mp_limb_t * powers,
    nmod_mpoly_ctx_t ctx)
{
    FLINT_ASSERT(A->bits == B->bits);
    nmod_mpoly_fit_length(A, B->length, ctx);
    slong N = mpoly_words_per_exp(B->bits, ctx->minfo);
    for (slong i = 0; i < B->length; i++)
    {
        mp_limb_t prod = UWORD(1);

        for (slong j = 0; j < entries; j++)
        {
            if ((B->exps + N*i)[offs[j]] & masks[j])
            {
                prod = nmod_mul(prod, powers[j], ctx->ffinfo->mod);
            }
        }

        A->coeffs[i] = prod;
        mpoly_monomial_set(A->exps + N*i, B->exps + N*i, N);
    }
    A->length = B->length;
}

void nmod_mpolyu_evalsk(
    nmod_mpolyu_t A,
    nmod_mpolyu_t B,
    slong entries,
    slong * offs,
    ulong * masks,
    mp_limb_t * powers,
    nmod_mpoly_ctx_t ctx)
{
    nmod_mpolyu_fit_length(A, B->length, ctx);
    for (slong i = 0; i < B->length; i++)
    {
        A->exps[i] = B->exps[i];
        nmod_mpoly_evalsk(A->coeffs + i, B->coeffs + i, entries, offs, masks, powers, ctx);
    }
    A->length = B->length;
}

void nmod_mpolyu_mulsk(nmod_mpolyu_t A, nmod_mpolyu_t B, nmod_mpoly_ctx_t ctx)
{
    FLINT_ASSERT(A->length == B->length);
    for (slong i = 0; i < A->length; i++)
    {
        FLINT_ASSERT(A->exps[i] == B->exps[i]);

        FLINT_ASSERT((A->coeffs + i)->length == (B->coeffs + i)->length);
        for (slong j = 0; j < (A->coeffs + i)->length; j++)
        {
            (A->coeffs + i)->coeffs[j] = nmod_mul((A->coeffs + i)->coeffs[j], (B->coeffs + i)->coeffs[j], ctx->ffinfo->mod);
        }
    }
}

void nmod_mpolyu_cvtto_poly(nmod_poly_t a, nmod_mpolyu_t A, nmod_mpoly_ctx_t ctx)
{
    nmod_poly_zero(a);
    for (slong i = 0; i < A->length; i++)
    {
        slong N = mpoly_words_per_exp((A->coeffs + i)->bits, ctx->minfo);
        FLINT_ASSERT((A->coeffs + i)->length == 1);
        FLINT_ASSERT(mpoly_monomial_is_zero((A->coeffs + i)->exps + N*0, N));
        nmod_poly_set_coeff_ui(a, A->exps[i], (A->coeffs + i)->coeffs[0]);
    }
}

void nmod_mpolyu_cvtfrom_poly(
    nmod_mpolyu_t A,
    nmod_poly_t a,
    nmod_mpoly_ctx_t ctx)
{
    nmod_mpolyu_zero(A, ctx);
    slong k = 0;
//flint_printf("nmod_mpolyu_cvtfrom_poly:\n");
//flint_printf("A bits: %d\n", A->bits);
//flint_printf("ctx minfo nfields: %d\n", ctx->minfo->nfields);

    slong N = mpoly_words_per_exp(A->bits, ctx->minfo);
//flint_printf("N: %d\n", N);

    for (slong i = nmod_poly_length(a) - 1; i >= 0; i--)
    {
        mp_limb_t c = nmod_poly_get_coeff_ui(a, i);

//flint_printf("i = %d  c = %d\n", i, c);

        if (c != UWORD(0))
        {

//flint_printf("ctx minfo nfields: %d\n", ctx->minfo->nfields);

            nmod_mpolyu_fit_length(A, k + 1, ctx);
            A->exps[k] = i;
            nmod_mpoly_fit_length(A->coeffs + k, 1, ctx);
            nmod_mpoly_fit_bits(A->coeffs + k, A->bits, ctx);
            (A->coeffs + k)->bits = A->bits;
            (A->coeffs + k)->coeffs[0] = c;
            (A->coeffs + k)->length = 1;
            mpoly_monomial_zero((A->coeffs + k)->exps + N*0, N);
//flint_printf("coeff: "); nmod_mpoly_print_pretty(A->coeffs + k, NULL, ctx); flint_printf("\n");
            k++;
        }
    }
    A->length = k;

//flint_printf("ret: "); nmod_mpolyu_print_pretty(A, NULL, ctx); flint_printf("\n");
}


/*
    return 0 if the leading coeff of A vanishes
    else return 1
*/
int nmod_mpolyu_evalfromsk(nmod_poly_t e, nmod_mpolyu_t A, nmod_mpolyu_t SK, nmod_mpoly_ctx_t ctx)
{
    int ret = 0;

//printf("nmod_mpolyu_evalfromsk:\n");
//printf(" A: "); nmod_mpolyu_print_pretty(A, NULL, ctx); printf("\n");
//printf("SK: "); nmod_mpolyu_print_pretty(SK, NULL, ctx); printf("\n");

    FLINT_ASSERT(A->length == SK->length);

    nmod_poly_zero(e);
    for (slong i = 0; i < A->length; i++)
    {
        FLINT_ASSERT((A->coeffs + i)->length == (SK->coeffs + i)->length);

        mp_limb_t v, pp0, pp1, ac0 = 0, ac1 = 0, ac2 = 0;
        for (slong j = 0; j < (A->coeffs + i)->length; j++)
        {
            umul_ppmm(pp1, pp0, (A->coeffs + i)->coeffs[j], (SK->coeffs + i)->coeffs[j]);
            add_sssaaaaaa(ac2, ac1, ac0, ac2, ac1, ac0, WORD(0), pp1, pp0);
        }
        NMOD_RED3(v, ac2, ac1, ac0, ctx->ffinfo->mod);

//flint_printf("coeff of S^%d is %d\n", A->exps[i], v);

        nmod_poly_set_coeff_ui(e, A->exps[i], v);
        ret |= (i == 0 && v != 0);
    }

    return ret;
}




int nmod_vandsolve(mp_limb_t * x, mp_limb_t * a, mp_limb_t * b, slong n, nmod_t mod)
{
//printf("nmod_vandsolve\n");

    int success = 0;
    slong i, j;
    mp_limb_t t;


//for(i=0; i<n; i++) {
//printf("a[%d]: %d  b[%d]: %d\n", i, a[i], i, b[i]);
//}

    for (i = 0; i < n; i++)
        x[i] = 0;

    nmod_poly_t Q, P, R, u;
    nmod_poly_init(Q, mod.n);
    nmod_poly_init(P, mod.n);
    nmod_poly_init(R, mod.n);
    nmod_poly_init(u, mod.n);
    nmod_poly_set_coeff_ui(u, 1, 1);
    nmod_poly_product_roots_nmod_vec(P, a, n);
    for (i = 0; i < n; i++)
    {
        if (a[i] == UWORD(0))
            goto cleanup;

        nmod_poly_set_coeff_ui(u, 0, mod.n - a[i]);
        nmod_poly_divrem(Q, R, P, u);
        t = nmod_mul(a[i], nmod_poly_evaluate_nmod(Q, a[i]), mod);
        if (t == UWORD(0))
            goto cleanup;

        mp_limb_t Dinv = nmod_inv(t, mod);
        for (j = 0; j < n; j++)
        {
            t = nmod_mul(b[j], Dinv, mod);
            t = nmod_mul(t, nmod_poly_get_coeff_ui(Q, j), mod);
            x[i] = nmod_add(x[i], t, mod);
        }
    }
    success = 1;

cleanup:
    nmod_poly_clear(Q);
    nmod_poly_clear(P);
    nmod_poly_clear(R);
    nmod_poly_clear(u);


//for(i=0; i<n; i++) {
//printf("x[%d]: %d\n", i, x[i]);
//}

//printf("success: %d\n", success);

    return success;
}


/*
    return -1: failed due to inability to find scale factors
    return  0: failed due to wrong assumed form
    return  1: success ("G" is correct assuming assumed form "f" is correct)
*/
int nmod_mpolyu_sgcd_zippel(
    nmod_mpolyu_t G,
    nmod_mpolyu_t A,
    nmod_mpolyu_t B,
    nmod_mpolyu_t f,
    slong var,
    nmod_mpoly_ctx_t ctx,
    flint_rand_t randstate)
{
    nmod_mpolyu_t Aevalsk1, Bevalsk1, fevalsk1, Aevalski, Bevalski, fevalski;
    nmod_poly_t Aeval, Beval, Geval;
    mp_limb_t * alpha, * b;
    nmod_mat_struct * M, * ML;
    nmod_mat_t MF, Mwindow, Msol;
    int lc_ok;
    int * ML_is_initialized;
    slong i, j, k, s, S, nullity, success;
    TMP_INIT;

flint_printf("sgcd:\n");
flint_printf("A: "); nmod_mpolyu_print_pretty(A, NULL, ctx); flint_printf("\n");
flint_printf("B: "); nmod_mpolyu_print_pretty(B, NULL, ctx); flint_printf("\n");
flint_printf("f: "); nmod_mpolyu_print_pretty(f, NULL, ctx); flint_printf("\n");


    FLINT_ASSERT(A->bits == B->bits);
    FLINT_ASSERT(A->bits == G->bits);
    FLINT_ASSERT(A->bits == f->bits);
    FLINT_ASSERT(var >= 0);
    if (var == 0)
    {
        nmod_poly_t a, b, g;
        nmod_poly_init(a, ctx->ffinfo->mod.n);
        nmod_poly_init(b, ctx->ffinfo->mod.n);
        nmod_poly_init(g, ctx->ffinfo->mod.n);
        nmod_mpolyu_cvtto_poly(a, A, ctx);
        nmod_mpolyu_cvtto_poly(b, B, ctx);
        nmod_poly_gcd(g, a, b);
        nmod_mpolyu_cvtfrom_poly(G, g, ctx);
        nmod_poly_clear(a);
        nmod_poly_clear(b);
        nmod_poly_clear(g);
        return 1;
    }

    if (f->length == 1)
    {
        return -1;
    }

    TMP_START;

    nmod_mpolyu_init(Aevalsk1, f->bits, ctx);
    nmod_mpolyu_init(Bevalsk1, f->bits, ctx);
    nmod_mpolyu_init(fevalsk1, f->bits, ctx);
    nmod_mpolyu_init(Aevalski, f->bits, ctx);
    nmod_mpolyu_init(Bevalski, f->bits, ctx);
    nmod_mpolyu_init(fevalski, f->bits, ctx);
    nmod_poly_init(Aeval, ctx->ffinfo->mod.n);
    nmod_poly_init(Beval, ctx->ffinfo->mod.n);
    nmod_poly_init(Geval, ctx->ffinfo->mod.n);



    slong * d = (slong *) TMP_ALLOC(f->length*sizeof(slong));
    for (i = 0; i < f->length; i++)
    {
        d[i] = i;
    }

    /*
        make d sort the coeffs so that
        (f->coeffs + d[j-1])->length <= (f->coeffs + d[j-0])->length
        for all j
    */
    for (i = 1; i<f->length; i++)
    {
        for (j=i; j > 1 && (f->coeffs + d[j-1])->length 
                         > (f->coeffs + d[j-0])->length; j--)
        {
            slong temp = d[j-1];
            d[j-1] = d[j-0];
            d[j-0] = temp;
        }
    }

    slong l = f->length - 3;
    for (i = 0; i < f->length; i++)
    {
        l += (f->coeffs + i)->length;
    }
    l = l / (f->length - 1);
    l = FLINT_MAX(l, (f->coeffs + f->length - 1)->length);

//flint_printf("l: %d\n", l);


    alpha = (mp_limb_t *) TMP_ALLOC(var*sizeof(mp_limb_t));
    ML = (nmod_mat_struct *) TMP_ALLOC(f->length*sizeof(nmod_mat_struct));
    b = (mp_limb_t *) TMP_ALLOC((f->coeffs + d[f->length - 1])->length*sizeof(mp_limb_t));

    nmod_mat_init(MF, 0, l, ctx->ffinfo->mod.n);

    M = (nmod_mat_struct *) TMP_ALLOC(f->length*sizeof(nmod_mat_struct));
    ML_is_initialized = (int *) TMP_ALLOC(f->length*sizeof(int));
    for (i = 0; i < f->length; i++)
    {
        nmod_mat_init(M + i, l, (f->coeffs + i)->length, ctx->ffinfo->mod.n);
        ML_is_initialized[i] = 0;
    }

    mp_limb_t * W = (mp_limb_t *) flint_malloc(l*f->length*sizeof(mp_limb_t));
    for (i = 0; i < l*f->length; i++)
    {
        W[i] = 0;
    }

    nmod_mat_init(Msol, l, 1, ctx->ffinfo->mod.n);

    /* compute how many masks are needed */
    slong entries = f->bits * var;
    slong * offs = (slong *) TMP_ALLOC(entries*sizeof(slong));
    ulong * masks = (ulong *) TMP_ALLOC(entries*sizeof(slong));
    mp_limb_t * powers = (mp_limb_t *) TMP_ALLOC(entries*sizeof(mp_limb_t));

    slong N = mpoly_words_per_exp(f->bits, ctx->minfo);




pick_evaluation_point:

    FLINT_ASSERT(ctx->ffinfo->mod.n > UWORD(2));

    for (i = 0; i < var; i++)
    {
        alpha[i] = UWORD(2) + n_randint(randstate, ctx->ffinfo->mod.n - UWORD(2));
    }

//for (i = 0; i <var; i++) {
//flint_printf("-------------------------evaluation_point  x%d: %d\n",i,alpha[i]);
//}


    /* store bit masks for each power of two of the non-main variables */
    for (i = 0; i < var; i++)
    {
        slong shift, off;
        mpoly_gen_offset_shift(&off, &shift, i, N, f->bits, ctx->minfo);
        for (j = 0; j < f->bits; j++)
        {
            offs[f->bits*i + j] = off;
            masks[f->bits*i + j] = UWORD(1) << (j + shift);
            if (j == 0)
                powers[f->bits*i + j] = alpha[i];
            else
                powers[f->bits*i + j] = nmod_mul(powers[f->bits*i + j-1], powers[f->bits*i + j-1], ctx->ffinfo->mod);
        }
    }
/*
for (i = 0; i < entries; i++)
{
printf("entry[%d]:  off = %d  mask = %016llx  power = %d\n", i, offs[i], masks[i], powers[i]);
}
*/
    nmod_mpolyu_evalsk(Aevalsk1, A, entries, offs, masks, powers, ctx);
    nmod_mpolyu_evalsk(Bevalsk1, B, entries, offs, masks, powers, ctx);
    nmod_mpolyu_evalsk(fevalsk1, f, entries, offs, masks, powers, ctx);

//printf("Aevalsk1: "); nmod_mpolyu_print_pretty(Aevalsk1, NULL, ctx); printf("\n");
//printf("Bevalsk1: "); nmod_mpolyu_print_pretty(Bevalsk1, NULL, ctx); printf("\n");
//printf("fevalsk1: "); nmod_mpolyu_print_pretty(fevalsk1, NULL, ctx); printf("\n");


    for (i = 0; i < l; i++)
    {
        if (i == 0)
        {
            nmod_mpolyu_set(Aevalski, Aevalsk1, ctx);
            nmod_mpolyu_set(Bevalski, Bevalsk1, ctx);
            nmod_mpolyu_set(fevalski, fevalsk1, ctx);
        } else
        {
            nmod_mpolyu_mulsk(Aevalski, Aevalsk1, ctx);
            nmod_mpolyu_mulsk(Bevalski, Bevalsk1, ctx);
            nmod_mpolyu_mulsk(fevalski, fevalsk1, ctx);
        }

        for (j = 0; j < f->length; j++)
        {
            for (k = 0; k < (f->coeffs + j)->length; k++)
            {
                (M + j)->rows[i][k] = (fevalski->coeffs + j)->coeffs[k];
            }
        }

        lc_ok = 1;
        lc_ok = lc_ok && nmod_mpolyu_evalfromsk(Aeval, A, Aevalski, ctx);
        lc_ok = lc_ok && nmod_mpolyu_evalfromsk(Beval, B, Bevalski, ctx);

//printf("Aeval: "); nmod_poly_print_pretty(Aeval, "X"); printf("\n");
//printf("Beval: "); nmod_poly_print_pretty(Beval, "X"); printf("\n");

        nmod_poly_gcd(Geval, Aeval, Beval);
//printf("Geval: "); nmod_poly_print_pretty(Geval, "X"); printf("\n");


        k = nmod_poly_length(Geval);
        j = WORD(0);
        while ((--k) >= 0)
        {
            mp_limb_t ck = nmod_poly_get_coeff_ui(Geval, k);
            if (ck != UWORD(0))
            {
                while (j < f->length && f->exps[j] > k)
                {
                    j++;
                }
                if (j >= f->length || f->exps[j] != k)
                {
                    /* Geval does not fit the form f */
                    success = 0;
                    goto finished;
                }
                W[l*j + i] = ck;
            }
        }
    }
/*
for (j = 0; j < f->length; j++) {
flint_printf("M[%d]:\n", j);
nmod_mat_print_pretty(M+j);
}
for (j = 0; j < l*f->length; j++) {
flint_printf("W[%d]: %d\n", j, W[j]);
}
*/

    nullity = -1;

    for (S = 0; S < f->length; S++)
    {
//flint_printf("starting S = %d\n",S);

        s = d[S];
        if (!ML_is_initialized[s])
        {
            nmod_mat_init(ML + s, l, (f->coeffs + s)->length + l, ctx->ffinfo->mod.n);
            for (i = 0; i < l; i++)
            {
                for (j = 0; j < (f->coeffs + s)->length; j++)
                {
                    (ML + s)->rows[i][j] = (M + s)->rows[i][j];
                }
                (ML + s)->rows[i][(f->coeffs + s)->length + i] = W[l*s + i];
            }
        } else {
            for (i = 0; i < l; i++)
            {
                for (j = 0; j < (f->coeffs + s)->length; j++)
                {
                    (ML + s)->rows[i][j] = (M + s)->rows[i][j];
                }
                for (j = 0; j < l; j++) {
                    (ML + s)->rows[i][(f->coeffs + s)->length + j]
                                             = (j==i ? W[l*s + i] : UWORD(0));
                }
            }

        }

//flint_printf("ML[%d]:\n", s);
//nmod_mat_print_pretty(ML+s);
//printf("\n");



        nmod_mat_rref(ML + s);

//printf("after rref ");
//flint_printf("ML[%d]:\n", s);
//nmod_mat_print_pretty(ML+s);
//printf("\n");


        for (i = 0; i < (f->coeffs + s)->length; i++)
        {
            if ((ML + s)->rows[i][i] != UWORD(1))
            {
                /* evaluation points produced a singular vandermonde matrix */
                goto pick_evaluation_point;
            }
        }

        nmod_mat_window_init(Mwindow, ML + s,
                (f->coeffs + s)->length, (f->coeffs + s)->length,
                 l, (f->coeffs + s)->length + l);
        nmod_mat_t MFtemp;
        nmod_mat_init(MFtemp, nmod_mat_nrows(MF) + l - (f->coeffs + s)->length, l, ctx->ffinfo->mod.n);
        nmod_mat_concat_vertical(MFtemp, MF, Mwindow);
        nmod_mat_swap(MFtemp, MF);
        nmod_mat_clear(MFtemp);
        nmod_mat_window_clear(Mwindow);


        nullity = l - nmod_mat_rref(MF);

//printf("rref MF:\n");
//nmod_mat_print_pretty(MF);
//printf("\n");
//flint_printf("nullity: %d\n",nullity);


        if (nullity == 0)
        {
            /* There is no solution for scale factors. Form f must be wrong */
            success = 0;
            goto finished;
        }
        if (nullity == 1)
        {
            /*
                There is one solution for scale factors based on equations
                considered thus far. Accept this as a solution and perform
                some checks of the remaining equations at the end.
            */
            break;
        }
    }

    if (nullity != 1)
    {
        /* Gcd might have content that needs factoring out */
        success = -1;
        //assert(0 && "oh man");
        goto finished;
    }

    nullity = nmod_mat_nullspace(Msol, MF);
    FLINT_ASSERT(nullity == 1);

//printf("Msol:\n");
//nmod_mat_print_pretty(Msol);
//printf("\n");

    nmod_mpolyu_set(G, f, ctx);

    for (i = 0; i < f->length; i++)
    {
        for (j = 0; j < (f->coeffs + i)->length; j++)
        {
//flint_printf("Msol[%d] = %d\n",j,nmod_mat_get_entry(Msol, j, 0));

            b[j] = nmod_mul(W[l*i + j], nmod_mat_get_entry(Msol, j, 0), ctx->ffinfo->mod);
        }
        success = nmod_vandsolve((G->coeffs + i)->coeffs,
                                 (fevalsk1->coeffs + i)->coeffs, b, 
                                    (f->coeffs + i)->length, ctx->ffinfo->mod);
        if (!success)
        {
            /* evaluation points produced a singular vandermonde matrix */
            goto pick_evaluation_point;
        }
    }

    /* check solution */
    for (s = 0; s < f->length; s++)
    {
        mp_limb_t pp0, pp1, ac0, ac1, ac2, u, v;

        for (i = 0; i < l; i++)
        {
            ac0 = ac1 = ac2 = 0;
            for (j = 0; j < (f->coeffs + s)->length; j++)
            {
                umul_ppmm(pp1, pp0, (M + s)->rows[i][j], (G->coeffs + s)->coeffs[j]);
                add_sssaaaaaa(ac2, ac1, ac0, ac2, ac1, ac0, WORD(0), pp1, pp0);
            }

            NMOD_RED3(v, ac2, ac1, ac0, ctx->ffinfo->mod);
            u = nmod_mul(W[l*s + i], nmod_mat_get_entry(Msol, i, 0), ctx->ffinfo->mod);
            if (v != u)
            {
                success = 0;
                printf("solution not a solution\n");
                goto finished;
            }
        }
    }

    success = 1;

finished:

    nmod_mat_clear(MF);
    nmod_mat_clear(Msol);
    for (i = 0; i < f->length; i++)
    {
        nmod_mat_clear(M + i);
        if (ML_is_initialized[i])
            nmod_mat_clear(ML + i);
    }
    nmod_mpolyu_clear(Aevalsk1, ctx);
    nmod_mpolyu_clear(Bevalsk1, ctx);
    nmod_mpolyu_clear(fevalsk1, ctx);
    nmod_mpolyu_clear(Aevalski, ctx);
    nmod_mpolyu_clear(Bevalski, ctx);
    nmod_mpolyu_clear(fevalski, ctx);
    nmod_poly_clear(Aeval);
    nmod_poly_clear(Beval);
    nmod_poly_clear(Geval);

printf("sgcd returning G: "); nmod_mpolyu_print_pretty(G, NULL, ctx); printf("\n");

    TMP_END;
    return success;
}









void nmod_mpolyn_eval_last(nmod_mpoly_t B, nmod_mpolyn_t A, mp_limb_t alpha, nmod_mpoly_ctx_t ctx)
{
//flint_printf("B bits = %d  A bits = %d\n",B->bits, A->bits);

    FLINT_ASSERT(B->bits == A->bits);

    nmod_mpoly_fit_length(B, A->length, ctx);
    slong N = mpoly_words_per_exp(A->bits, ctx->minfo);
    slong k = 0;
    for (slong i = 0; i < A->length; i++)
    {
        mpoly_monomial_set(B->exps + N*k, A->exps + N*i, N);
        B->coeffs[k] = nmod_poly_evaluate_nmod(A->coeffs + i, alpha);
        if (B->coeffs[k] != UWORD(0))
        {
            k++;
        }
    }
    B->length = k;
}

void nmod_mpolyun_eval_last(nmod_mpolyu_t B, nmod_mpolyun_t A, mp_limb_t alpha, nmod_mpoly_ctx_t ctx)
{
    FLINT_ASSERT(B->bits == A->bits);

//flint_printf("un_eval_last\n");

    nmod_mpolyu_fit_length(B, A->length, ctx);
    slong k = 0;
    for (slong i = 0; i < A->length; i++)
    {
        B->exps[k] = A->exps[i];
        nmod_mpolyn_eval_last(B->coeffs + k, A->coeffs + i, alpha, ctx);
        if (!nmod_mpoly_is_zero(B->coeffs + k, ctx))
        {
            k++;
        }
    }
    B->length = k;
}



void nmod_mpoly_setform(nmod_mpoly_t A, nmod_mpoly_t B, nmod_mpoly_ctx_t ctx)
{
    nmod_mpoly_set(A, B, ctx);
    for (slong i = 0; i < A->length; i++)
    {
        A->coeffs[i] = UWORD(0);
    }
}

void nmod_mpolyu_setform(nmod_mpolyu_t A, nmod_mpolyu_t B, nmod_mpoly_ctx_t ctx)
{
    nmod_mpolyu_fit_length(A, B->length, ctx);
    for (slong i = 0; i < B->length; i++)
    {
        nmod_mpoly_setform(A->coeffs + i, B->coeffs + i, ctx);
        A->exps[i] = B->exps[i];
    }
    A->length = B->length;
}


void nmod_mpoly_setform_mpolyn(nmod_mpoly_t A, nmod_mpolyn_t B, nmod_mpoly_ctx_t ctx)
{
    slong i, N;
    FLINT_ASSERT(A->bits == B->bits);

    nmod_mpoly_fit_length(A, B->length, ctx);
    N = mpoly_words_per_exp(B->bits, ctx->minfo);
    for (i = 0; i < B->length; i++)
    {
        A->coeffs[i] = UWORD(0);
        mpoly_monomial_set(A->exps + N*i, B->exps + N*i, N);
    }
    A->length = B->length;
}

void nmod_mpolyu_setform_mpolyun(nmod_mpolyu_t A, nmod_mpolyun_t B, nmod_mpoly_ctx_t ctx)
{
    slong i;

    nmod_mpolyu_fit_length(A, B->length, ctx);
    for (i = 0; i < B->length; i++)
    {
        FLINT_ASSERT((B->coeffs + i)->bits == B->bits);
        nmod_mpoly_setform_mpolyn(A->coeffs + i, B->coeffs + i, ctx);
        A->exps[i] = B->exps[i];
    }
    A->length = B->length;
}


void nmod_mpoly_cvtfrom_poly_notmain(nmod_mpoly_t A, nmod_poly_t a, slong var, nmod_mpoly_ctx_t ctx)
{
    TMP_INIT;
    TMP_START;

//printf("nmod_mpoly_cvtfrom_poly_notmain:\n");
//printf("a: "); nmod_poly_print_pretty(a, "v"); printf("\n");

    slong N = mpoly_words_per_exp(A->bits, ctx->minfo);

    ulong * oneexp = (ulong *)TMP_ALLOC(N*sizeof(ulong));
    slong offset, shift;
    mpoly_gen_oneexp_offset_shift(oneexp, &offset, &shift, var, N, A->bits, ctx->minfo);

    nmod_mpoly_fit_length(A, nmod_poly_length(a), ctx);

    slong k = 0;
    for (slong i = nmod_poly_length(a) - 1; i >= 0; i--)
    {
        mp_limb_t c = nmod_poly_get_coeff_ui(a, i);
        if (c != UWORD(0))
        {
            A->coeffs[k] = c;
            mpoly_monomial_mul_si(A->exps + N*k, oneexp, N, i);
            k++;
        }
    }
    A->length = k;

//printf("A: "); nmod_mpoly_print_pretty(A, NULL, ctx); printf("\n");

    TMP_END;
}

void nmod_mpolyu_cvtfrom_poly_notmain(nmod_mpolyu_t A, nmod_poly_t a, slong var, nmod_mpoly_ctx_t ctx)
{
//printf("nmod_mpolyu_cvtfrom_poly_notmain:\n");

    nmod_mpolyu_fit_length(A, 1, ctx);
    A->exps[0] = 0;
    nmod_mpoly_cvtfrom_poly_notmain(A->coeffs + 0, a, var, ctx);
    A->length = !nmod_mpoly_is_zero(A->coeffs + 0, ctx);

//printf("A: "); nmod_mpoly_print_pretty(A, NULL, ctx); printf("\n");

}



void nmod_mpolyu_scalar_mul_nmod(nmod_mpolyu_t A, mp_limb_t c, nmod_mpoly_ctx_t ctx)
{
    for (slong i = 0; i < A->length; i++)
    {
        for (slong j = 0; j < (A->coeffs + i)->length; j++)
        {
            (A->coeffs + i)->coeffs[j] = nmod_mul((A->coeffs + i)->coeffs[j], c, ctx->ffinfo->mod);
        }
    }
}


int nmod_mpolyu_pgcd_zippel_univar(
    nmod_mpolyu_t G,
    nmod_mpolyu_t A,
    nmod_mpolyu_t B,
    nmod_mpoly_ctx_t ctx)
{
    FLINT_ASSERT(A->bits == B->bits);
    nmod_poly_t a, b, g;
    nmod_poly_init(a, ctx->ffinfo->mod.n);
    nmod_poly_init(b, ctx->ffinfo->mod.n);
    nmod_poly_init(g, ctx->ffinfo->mod.n);
    nmod_mpolyu_cvtto_poly(a, A, ctx);
    nmod_mpolyu_cvtto_poly(b, B, ctx);
    nmod_poly_gcd(g, a, b);
    nmod_mpolyu_cvtfrom_poly(G, g, ctx);
    nmod_poly_clear(a);
    nmod_poly_clear(b);
    nmod_poly_clear(g);
    return 1;
}


int nmod_mpolyu_pgcd_zippel_bivar(
    nmod_mpolyu_t G,
    nmod_mpolyu_t A,
    nmod_mpolyu_t B,
    nmod_mpoly_ctx_t ctx,
    mpoly_zipinfo_t zinfo)
{
    nmod_poly_t a, b, c, g, modulus, tempmod;
    nmod_mpolyu_t Aeval, Beval, Geval;
    nmod_mpolyun_t An, Bn, H, Ht;
    mp_limb_t geval, temp, alpha;
    int success = 0, changed;
    slong lastdeg;
    slong var = 0;

    FLINT_ASSERT(ctx->minfo->ord == ORD_LEX);
    FLINT_ASSERT(var >= -WORD(1));

    FLINT_ASSERT(G->bits == A->bits);
    FLINT_ASSERT(G->bits == B->bits);

//flint_printf("**********nmod_mpolyu_pgcd_zippel_bivar var = %d *********\n", var);
//flint_printf("A: "); nmod_mpolyu_print_pretty(A, NULL, ctx); flint_printf("\n");
//flint_printf("B: "); nmod_mpolyu_print_pretty(B, NULL, ctx); flint_printf("\n");

    nmod_mpolyun_init(An, A->bits, ctx);
    nmod_mpolyun_init(Bn, A->bits, ctx);

    nmod_mpolyu_cvtto_mpolyun(An, A, var, ctx);
    nmod_mpolyu_cvtto_mpolyun(Bn, B, var, ctx);

    nmod_poly_init(a, ctx->ffinfo->mod.n);
    nmod_poly_init(b, ctx->ffinfo->mod.n);
    nmod_poly_init(c, ctx->ffinfo->mod.n);
    nmod_poly_init(g, ctx->ffinfo->mod.n);
    nmod_mpolyun_content_last(a, An, ctx);
    nmod_mpolyun_content_last(b, Bn, ctx);
    nmod_mpolyun_divexact_last(An, a, ctx);
    nmod_mpolyun_divexact_last(Bn, b, ctx);
    nmod_poly_gcd(c, a, b);
    nmod_poly_gcd(g, nmod_mpolyun_leadcoeff_ref(An, ctx),
                     nmod_mpolyun_leadcoeff_ref(Bn, ctx));

//printf("c: "); nmod_poly_print_pretty(c, "v"); printf("\n");
//printf("g: "); nmod_poly_print_pretty(g, "v"); printf("\n");

    nmod_poly_init(modulus, ctx->ffinfo->mod.n);
    nmod_poly_init(tempmod, ctx->ffinfo->mod.n);
    nmod_poly_set_coeff_ui(tempmod, 1, UWORD(1));
    nmod_mpolyu_init(Aeval, A->bits, ctx);
    nmod_mpolyu_init(Beval, A->bits, ctx);
    nmod_mpolyu_init(Geval, A->bits, ctx);
    nmod_mpolyun_init(H, A->bits, ctx);
    nmod_mpolyun_init(Ht, A->bits, ctx);

    nmod_poly_one(modulus);
    nmod_mpolyun_zero(H, ctx);

    alpha = ctx->ffinfo->mod.n;
    while (1)
    {
        if (alpha == UWORD(0))
            goto ret_fail;
        alpha--;

//flint_printf("starting loop x%d = %wu\n", var, alpha);
//usleep(100000);

        geval = nmod_poly_evaluate_nmod(g, alpha);
        if (geval == WORD(0))
            goto outer_continue;

        nmod_mpolyun_eval_last(Aeval, An, alpha, ctx);
        nmod_mpolyun_eval_last(Beval, Bn, alpha, ctx);

        nmod_mpolyu_pgcd_zippel_univar(Geval, Aeval, Beval, ctx);


        if (nmod_mpolyu_is_one(Geval, ctx))
        {
            nmod_mpolyu_cvtfrom_poly_notmain(G, c, var, ctx);
            goto ret_success;
        }

//flint_printf("Geval: "); nmod_mpolyu_print_pretty(Geval, NULL, ctx); flint_printf("\n");

        FLINT_ASSERT(Geval->length > 0);

        if (nmod_poly_degree(modulus) > 0)
        {
            if (Geval->exps[0] > H->exps[0])
            {
                goto outer_continue;

            } else if (Geval->exps[0] < H->exps[0])
            {
                nmod_poly_one(modulus);                
            }
        }

        temp = n_invmod(nmod_mpolyu_leadcoeff(Geval, ctx), ctx->ffinfo->mod.n);
        nmod_mpolyu_scalar_mul_nmod(Geval, nmod_mul(geval, temp, ctx->ffinfo->mod), ctx);

        /* update interpolant H */
//flint_printf("updating H\n");

        if (nmod_poly_degree(modulus) > 0)
        {
            nmod_poly_scalar_mul_nmod(modulus, modulus, n_invmod(nmod_poly_evaluate_nmod(modulus, alpha), ctx->ffinfo->mod.n));
//flint_printf("before modulus: "); nmod_poly_print_pretty(modulus, "v"); printf("\n");
            changed = nmod_mpolyun_addinterp(&lastdeg, H, Ht, Geval, modulus, alpha, ctx);
            if (!changed)
            {
                nmod_mpolyun_content_last(a, H, ctx);
                nmod_mpolyun_mul_poly(Ht, H, c, ctx);
                nmod_mpolyun_divexact_last(Ht, a, ctx);

                nmod_mpolyu_cvtfrom_mpolyun(G, Ht, var, ctx);

//flint_printf("checking: "); nmod_mpolyu_print_pretty(G, NULL, ctx); printf("\n");

                if (!nmod_mpolyu_divides(A, G, ctx))
                    goto outer_continue;
                if (!nmod_mpolyu_divides(B, G, ctx))
                    goto outer_continue;        
                goto ret_success;
            }

        } else
        {
            nmod_mpolyun_set_mpolyu(H, Geval, ctx);
            lastdeg = WORD(0);
        }
        nmod_poly_set_coeff_ui(tempmod, 0, ctx->ffinfo->mod.n - alpha);
        nmod_poly_mul(modulus, modulus, tempmod);

//flint_printf("H: "); nmod_mpolyun_print_pretty(H, NULL, ctx); flint_printf("\n");

/*
        if (lastdeg > zinfo->Adegs[var] || lastdeg > zinfo->Bdegs[var])
        {
//            flint_printf("lastdeg: %d, %d %d\n",lastdeg, zinfo->Adegs[var], zinfo->Bdegs[var]);
            nmod_poly_one(modulus);
            goto outer_continue;
        }
*/

//flint_printf("after modulus: "); nmod_poly_print_pretty(modulus, "v"); printf("\n");


outer_continue:
        NULL;
    }

ret_success:
    success = 1;

ret_fail:
    nmod_poly_clear(a);
    nmod_poly_clear(b);
    nmod_poly_clear(c);
    nmod_poly_clear(g);
    nmod_poly_clear(modulus);
    nmod_poly_clear(tempmod);
    nmod_mpolyu_clear(Aeval, ctx);
    nmod_mpolyu_clear(Beval, ctx);
    nmod_mpolyu_clear(Geval, ctx);
    nmod_mpolyun_clear(An, ctx);
    nmod_mpolyun_clear(Bn, ctx);
    nmod_mpolyun_clear(H, ctx);
    nmod_mpolyun_clear(Ht, ctx);

    return success;
}



void nmod_mpolyun_shift_right(nmod_mpolyun_t A, ulong s)
{
    slong i;
    for (i = 0; i < A->length; i++)
    {
        FLINT_ASSERT(A->exps[i] >= s);
        A->exps[i] -= s;
    }
}

void nmod_mpolyun_shift_left(nmod_mpolyun_t A, ulong s)
{
    slong i;
    for (i = 0; i < A->length; i++)
    {
        A->exps[i] += s;
    }
}

int nmod_mpolyu_pgcd_zippel(
    nmod_mpolyu_t G,
    nmod_mpolyu_t A,
    nmod_mpolyu_t B,
    slong var,
    nmod_mpoly_ctx_t ctx,
    mpoly_zipinfo_t zinfo,
    flint_rand_t randstate)
{
    int divcheck_fail_count, inner_gcd_added;
    slong lastdeg;
    ulong ABminshift;
    int success = 0, changed;

    FLINT_ASSERT(ctx->minfo->ord == ORD_LEX);
    FLINT_ASSERT(var >= -WORD(1));

    FLINT_ASSERT(G->bits == A->bits);
    FLINT_ASSERT(G->bits == B->bits);

//flint_printf("**********nmod_mpolyu_pgcd_zippel var = %d *********\n", var);
//flint_printf("A: "); nmod_mpolyu_print_pretty(A, NULL, ctx); flint_printf("\n");
//flint_printf("B: "); nmod_mpolyu_print_pretty(B, NULL, ctx); flint_printf("\n");

    FLINT_ASSERT(A->length > 0);
    FLINT_ASSERT(B->length > 0);

    if (var == -WORD(1))
    {
        /* no more variables left to interpolate */
        return nmod_mpolyu_pgcd_zippel_univar(G, A, B, ctx);
    }

    if (var == WORD(0))
    {
        /* bivariate is more comfortable separated */
        return nmod_mpolyu_pgcd_zippel_bivar(G, A, B, ctx, zinfo);
    }

    nmod_mpolyun_t An, Bn;
    nmod_mpolyun_init(An, A->bits, ctx);
    nmod_mpolyun_init(Bn, A->bits, ctx);

    nmod_mpolyu_cvtto_mpolyun(An, A, var, ctx);
    nmod_mpolyu_cvtto_mpolyun(Bn, B, var, ctx);
    ABminshift = FLINT_MIN(A->exps[A->length - 1], B->exps[B->length - 1]);
    nmod_mpolyun_shift_right(An, ABminshift);
    nmod_mpolyun_shift_right(Bn, ABminshift);

flint_printf("An: "); nmod_mpolyun_print_pretty(An, NULL, ctx); flint_printf("\n");
flint_printf("Bn: "); nmod_mpolyun_print_pretty(Bn, NULL, ctx); flint_printf("\n");


    nmod_poly_t a, b, c, g;
    nmod_poly_init(a, ctx->ffinfo->mod.n);
    nmod_poly_init(b, ctx->ffinfo->mod.n);
    nmod_poly_init(c, ctx->ffinfo->mod.n);
    nmod_poly_init(g, ctx->ffinfo->mod.n);

    nmod_mpolyun_content_last(a, An, ctx);
    nmod_mpolyun_content_last(b, Bn, ctx);

flint_printf("a: "); nmod_poly_print_pretty(a, "v"); flint_printf("\n");
flint_printf("b: "); nmod_poly_print_pretty(b, "v"); flint_printf("\n");

    nmod_mpolyun_divexact_last(An, a, ctx);
    nmod_mpolyun_divexact_last(Bn, b, ctx);

flint_printf("An: "); nmod_mpolyun_print_pretty(An, NULL, ctx); flint_printf("\n");
flint_printf("Bn: "); nmod_mpolyun_print_pretty(Bn, NULL, ctx); flint_printf("\n");


    nmod_poly_gcd(c, a, b);

flint_printf("c: "); nmod_poly_print_pretty(c, "v"); flint_printf("\n");

    nmod_poly_gcd(g, nmod_mpolyun_leadcoeff_ref(An, ctx),
                     nmod_mpolyun_leadcoeff_ref(Bn, ctx));

flint_printf("g: "); nmod_poly_print_pretty(c, "v"); flint_printf("\n");

    nmod_poly_t modulus, tempmod;
    nmod_poly_init(modulus, ctx->ffinfo->mod.n);
    nmod_poly_init(tempmod, ctx->ffinfo->mod.n);
    nmod_poly_set_coeff_ui(tempmod, 1, UWORD(1));
    nmod_mpolyu_t Aeval, Beval, Geval, Gform;
    nmod_mpolyu_init(Aeval, A->bits, ctx);
    nmod_mpolyu_init(Beval, A->bits, ctx);
    nmod_mpolyu_init(Geval, A->bits, ctx);
    nmod_mpolyu_init(Gform, A->bits, ctx);
    nmod_mpolyun_t H, Ht;
    nmod_mpolyun_init(H, A->bits, ctx);
    nmod_mpolyun_init(Ht, A->bits, ctx);
    mp_limb_t geval, temp;


    nmod_poly_one(modulus);
    nmod_mpolyun_zero(H, ctx);

    mp_limb_t alpha = ctx->ffinfo->mod.n;


    while (1)
    {
        if (alpha == UWORD(0))
        {
            success = 0;
            goto finished;
        }
        alpha--;

flint_printf("------starting outer loop x%d = %d\n", var, alpha);


        geval = nmod_poly_evaluate_nmod(g, alpha);
        if (geval == WORD(0))
            goto outer_continue;

//flint_printf("Aeval bits: %d\n", Aeval->bits);
//flint_printf("   An bits: %d\n",    An->bits);
        nmod_mpolyun_eval_last(Aeval, An, alpha, ctx);

        nmod_mpolyun_eval_last(Beval, Bn, alpha, ctx);
        success = nmod_mpolyu_pgcd_zippel(Geval, Aeval, Beval, var - 1, ctx, zinfo, randstate);
flint_printf("********** pgcd return(%d) var = %d ", var-1, success); nmod_mpolyu_print_pretty(Geval, NULL, ctx); flint_printf("\n");
//flint_printf("Geval bits: %d\n", Geval->bits);
        if (!success)
            goto outer_continue;

        if (nmod_mpolyu_is_one(Geval, ctx))
        {
            nmod_mpolyu_cvtfrom_poly_notmain(G, c, var, ctx);
            return 1;
        }

        if (nmod_poly_degree(modulus) > 0)
        {
            if (Geval->exps[0] > H->exps[0])
            {
                goto outer_continue;

            } else if (Geval->exps[0] < H->exps[0])
            {
                nmod_poly_one(modulus);                
            }
        }

        /* update interpolant H */
flint_printf("------updating H\n");

        temp = n_invmod(nmod_mpolyu_leadcoeff(Geval, ctx), ctx->ffinfo->mod.n);
        nmod_mpolyu_scalar_mul_nmod(Geval, nmod_mul(geval, temp, ctx->ffinfo->mod), ctx);

        if (nmod_poly_degree(modulus) > 0)
        {
            nmod_poly_scalar_mul_nmod(modulus, modulus, n_invmod(nmod_poly_evaluate_nmod(modulus, alpha), ctx->ffinfo->mod.n));
            changed = nmod_mpolyun_addinterp(&lastdeg, H, Ht, Geval, modulus, alpha, ctx);
            if (!changed)
            {
                nmod_mpolyun_content_last(a, H, ctx);
                nmod_mpolyun_mul_poly(Ht, H, c, ctx);
                nmod_mpolyun_divexact_last(Ht, a, ctx);
                nmod_mpolyun_shift_left(Ht, ABminshift);
                nmod_mpolyu_cvtfrom_mpolyun(G, Ht, var, ctx);
                if (!nmod_mpolyu_divides(A, G, ctx))
                    goto outer_continue;
                if (!nmod_mpolyu_divides(B, G, ctx))
                    goto outer_continue;
                return 1;
            }

        } else
        {
            nmod_mpolyun_set_mpolyu(H, Geval, ctx);
        }
        nmod_poly_set_coeff_ui(tempmod, 0, ctx->ffinfo->mod.n - alpha);
        nmod_poly_mul(modulus, modulus, tempmod);

//flint_printf("after modulus: "); nmod_poly_print_pretty(modulus, "v"); printf("\n");


flint_printf("------H: "); nmod_mpolyun_print_pretty(H, NULL, ctx); flint_printf("\n");


        nmod_mpolyu_setform_mpolyun(Gform, H, ctx);

        divcheck_fail_count = 0;
        inner_gcd_added = 0;
        while (1)
        {
//flint_printf("asdf\n");
            if (alpha == UWORD(0))
            {
                success = 0;
                goto finished;
            }
            alpha--;

flint_printf("------starting inner loop x%d = %d\n", var, alpha);
        assert(alpha >= 80);

            geval = nmod_poly_evaluate_nmod(g, alpha);
            if (geval == WORD(0))
                goto inner_continue;

            nmod_mpolyun_eval_last(Aeval, An, alpha, ctx);
            nmod_mpolyun_eval_last(Beval, Bn, alpha, ctx);
            success = nmod_mpolyu_sgcd_zippel(Geval, Aeval, Beval, Gform, var, ctx, randstate);
flint_printf("sgcd return(%d) var = %d : ", success, var); nmod_mpolyu_print_pretty(Geval, NULL, ctx); flint_printf("\n");
//flint_printf("Geval bits: %d\n", Geval->bits);

            if (success == -1)
            {
                success = nmod_mpolyu_pgcd_zippel_rmcontent(G, A, B, var, ctx, zinfo, randstate);
                goto finished;
            }

            if (success != 1)
            {
                if (inner_gcd_added)
                    nmod_poly_one(modulus);
                goto outer_continue;
            }

            /* update interpolant H */
            temp = n_invmod(nmod_mpolyu_leadcoeff(Geval, ctx), ctx->ffinfo->mod.n);
            nmod_mpolyu_scalar_mul_nmod(Geval, nmod_mul(geval, temp, ctx->ffinfo->mod), ctx);
            FLINT_ASSERT(nmod_poly_degree(modulus) > 0);
            nmod_poly_scalar_mul_nmod(modulus, modulus, n_invmod(nmod_poly_evaluate_nmod(modulus, alpha), ctx->ffinfo->mod.n));
//flint_printf("Geval bits: %d\n", Geval->bits);
//flint_printf("before modulus: "); nmod_poly_print_pretty(modulus, "v"); printf("\n");
//flint_printf("before add     H: "); nmod_mpolyun_print_pretty(H, NULL, ctx); flint_printf("\n");
//flint_printf("before add Geval: "); nmod_mpolyu_print_pretty(Geval, NULL, ctx); flint_printf("\n");

            changed = nmod_mpolyun_addinterp(&lastdeg, H, Ht, Geval, modulus, alpha, ctx);
            nmod_poly_set_coeff_ui(tempmod, 0, ctx->ffinfo->mod.n - alpha);
            nmod_poly_mul(modulus, modulus, tempmod);
flint_printf("------updated H: "); nmod_mpolyun_print_pretty(H, NULL, ctx); flint_printf("\n");

            inner_gcd_added = changed;

            if (!changed)
            {
printf("didnt change\n");

                nmod_mpolyun_content_last(a, H, ctx);
                nmod_mpolyun_mul_poly(Ht, H, c, ctx);
                nmod_mpolyun_divexact_last(Ht, a, ctx);

                nmod_mpolyu_cvtfrom_mpolyun(G, Ht, var, ctx);

printf("candidate G: "); nmod_mpolyu_print_pretty(G, NULL, ctx); printf("\n");

                if (!nmod_mpolyu_divides(A, G, ctx)
                             || !nmod_mpolyu_divides(B, G, ctx))
                {
                    ++divcheck_fail_count;
                    flint_printf("division failed %d times\n", divcheck_fail_count);
                    if (divcheck_fail_count >= 2)
                    {
                        nmod_poly_one(modulus);
                        goto outer_continue;
                    }
                    goto inner_continue;
                }
                success = 1;
                goto finished;
            }

//flint_printf("after modulus: "); nmod_poly_print_pretty(modulus, "v"); printf("\n");

inner_continue:
            NULL;
        }
        FLINT_ASSERT(0 && "not reachable");

outer_continue:
        NULL;
    }
    FLINT_ASSERT(0 && "not reachable");


finished:

    nmod_poly_clear(a);
    nmod_poly_clear(b);
    nmod_poly_clear(c);
    nmod_poly_clear(g);
    nmod_poly_clear(modulus);
    nmod_poly_clear(tempmod);
    nmod_mpolyu_clear(Aeval, ctx);
    nmod_mpolyu_clear(Beval, ctx);
    nmod_mpolyu_clear(Geval, ctx);
    nmod_mpolyun_clear(An, ctx);
    nmod_mpolyun_clear(Bn, ctx);
    nmod_mpolyun_clear(H, ctx);
    nmod_mpolyun_clear(Ht, ctx);

    return success;
}


int nmod_mpolyu_content(
    nmod_mpoly_t c,
    nmod_mpolyu_t A,
    nmod_mpoly_ctx_t ctx)
{
    slong i;
    FLINT_ASSERT(A->length > 0);
flint_printf("starting content. bits = %wd\n", A->bits);

    nmod_mpoly_set(c, A->coeffs + 0, ctx);
flint_printf("c->bits = %wd\n", c->bits);

    FLINT_ASSERT(c->bits == A->bits);
    for (i = 1; i < A->length; i++)
    {
        int success;
        success = nmod_mpoly_gcd_zippel_keepbits(c, c, A->coeffs + i, ctx);
        if (!success)
            return 0;
        FLINT_ASSERT(c->bits == A->bits);
    }
    return 1;
}


slong _nmod_mpoly_divides_monagan_pearce(
                     mp_limb_t ** coeff1,      ulong ** exp1, slong * alloc,
                const mp_limb_t * coeff2, const ulong * exp2, slong len2,
                const mp_limb_t * coeff3, const ulong * exp3, slong len3,
     mp_bitcnt_t bits, slong N, const ulong * cmpmask, const nmodf_ctx_t fctx);

void nmod_mpolyu_divexact_mpoly(nmod_mpolyu_t A, nmod_mpolyu_t B, nmod_mpoly_t c, nmod_mpoly_ctx_t ctx)
{
    slong len, N;
    mp_bitcnt_t exp_bits;
    nmod_mpoly_struct * poly1, * poly2, * poly3;
    ulong * cmpmask;
    TMP_INIT;

    TMP_START;

    exp_bits = B->bits;
    FLINT_ASSERT(A->bits == B->bits);

    nmod_mpolyu_fit_length(A, B->length, ctx);

    N = mpoly_words_per_exp(exp_bits, ctx->minfo);
    cmpmask = (ulong*) TMP_ALLOC(N*sizeof(ulong));
    mpoly_get_cmpmask(cmpmask, N, exp_bits, ctx->minfo);


    for (slong i = 0; i < B->length; i++)
    {

        flint_printf("dividing i = %d\n",i);
        flint_printf("(B->coeffs + i)->bits: %wd\n", (B->coeffs + i)->bits);
        flint_printf("              c->bits: %wd\n", c->bits);

        poly1 = A->coeffs + i;
        poly2 = B->coeffs + i;
        poly3 = c;

      nmod_mpoly_fit_length(poly1, poly2->length/poly3->length + 1, ctx);
      nmod_mpoly_fit_bits(poly1, exp_bits, ctx);
      poly1->bits = exp_bits;

      len = _nmod_mpoly_divides_monagan_pearce(&poly1->coeffs, &poly1->exps,
                            &poly1->alloc, poly2->coeffs, poly2->exps, poly2->length,
                              poly3->coeffs, poly3->exps, poly3->length, exp_bits, N,
                                                  cmpmask, ctx->ffinfo);

    _nmod_mpoly_set_length(poly1, len, ctx);

        FLINT_ASSERT(len > 0);

        flint_printf("(A->coeffs + i)->bits: %wd\n", (A->coeffs + i)->bits);
        FLINT_ASSERT((A->coeffs + i)->bits == A->bits);


        A->exps[i] = B->exps[i];
    }
    A->length = B->length;

    TMP_END;
}


slong _nmod_mpoly_mul_johnson(mp_limb_t ** coeff1, ulong ** exp1, slong * alloc,
                 const mp_limb_t * coeff2, const ulong * exp2, slong len2,
                 const mp_limb_t * coeff3, const ulong * exp3, slong len3,
      mp_bitcnt_t bits, slong N, const ulong * cmpmask, const nmodf_ctx_t fctx);

void nmod_mpolyu_mul_mpoly(nmod_mpolyu_t A, nmod_mpolyu_t B, nmod_mpoly_t c, nmod_mpoly_ctx_t ctx)
{
    slong len, N;
    mp_bitcnt_t exp_bits;
    nmod_mpoly_struct * poly1, * poly2, * poly3;
    ulong * cmpmask;
    TMP_INIT;

    TMP_START;

    exp_bits = B->bits;
    FLINT_ASSERT(A->bits == B->bits);

    nmod_mpolyu_fit_length(A, B->length, ctx);

    N = mpoly_words_per_exp(exp_bits, ctx->minfo);
    cmpmask = (ulong*) TMP_ALLOC(N*sizeof(ulong));
    mpoly_get_cmpmask(cmpmask, N, exp_bits, ctx->minfo);


    for (slong i = 0; i < B->length; i++)
    {

        flint_printf("dividing i = %d\n",i);
        flint_printf("(B->coeffs + i)->bits: %wd\n", (B->coeffs + i)->bits);
        flint_printf("              c->bits: %wd\n", c->bits);

        poly1 = A->coeffs + i;
        poly2 = B->coeffs + i;
        poly3 = c;

      nmod_mpoly_fit_length(poly1, poly2->length/poly3->length + 1, ctx);
      nmod_mpoly_fit_bits(poly1, exp_bits, ctx);
      poly1->bits = exp_bits;

      len = _nmod_mpoly_mul_johnson(&poly1->coeffs, &poly1->exps,
                            &poly1->alloc, poly2->coeffs, poly2->exps, poly2->length,
                              poly3->coeffs, poly3->exps, poly3->length, exp_bits, N,
                                                  cmpmask, ctx->ffinfo);

    _nmod_mpoly_set_length(poly1, len, ctx);

        FLINT_ASSERT(len > 0);

        flint_printf("(A->coeffs + i)->bits: %wd\n", (A->coeffs + i)->bits);
        FLINT_ASSERT((A->coeffs + i)->bits == A->bits);


        A->exps[i] = B->exps[i];
    }
    A->length = B->length;

    TMP_END;
}






int nmod_mpolyu_pgcd_zippel_rmcontent(
    nmod_mpolyu_t G,
    nmod_mpolyu_t A,
    nmod_mpolyu_t B,
    slong var,
    nmod_mpoly_ctx_t ctx,
    mpoly_zipinfo_t zinfo,
    flint_rand_t randstate)
{
    int success;
    nmod_mpoly_t Ac, Bc, Gc;
    nmod_mpolyu_t Abar, Bbar, Gbar;

    FLINT_ASSERT(A->bits == B->bits);
    FLINT_ASSERT(A->bits == G->bits);

    nmod_mpoly_init(Ac, ctx);
    nmod_mpoly_init(Bc, ctx);
    nmod_mpoly_init(Gc, ctx);

    nmod_mpolyu_init(Abar, A->bits, ctx);
    nmod_mpolyu_init(Bbar, A->bits, ctx);
    nmod_mpolyu_init(Gbar, A->bits, ctx);

printf("!!!!!!!!!!!!!!!!!!!1 removing content\n");


    success = 1;
    success = success && nmod_mpolyu_content(Ac, A, ctx);
    success = success && nmod_mpolyu_content(Bc, B, ctx);
    if (!success)
    {
        assert(0 && "content failed");
        goto finished;
    }

printf("Ac: "); nmod_mpoly_print_pretty(Ac, NULL, ctx); printf("\n");
printf("Bc: "); nmod_mpoly_print_pretty(Bc, NULL, ctx); printf("\n");



    nmod_mpolyu_divexact_mpoly(Abar, A, Ac, ctx);
    nmod_mpolyu_divexact_mpoly(Bbar, B, Bc, ctx);


printf("Abar: "); nmod_mpolyu_print_pretty(Abar, NULL, ctx); printf("\n");
printf("Bbar: "); nmod_mpolyu_print_pretty(Bbar, NULL, ctx); printf("\n");




    success = nmod_mpoly_gcd_zippel_keepbits(Gc, Ac, Bc, ctx);
    if (!success)
    {
        assert(0 && "content gcd failed");
        goto finished;
    }

printf("Gc: "); nmod_mpoly_print_pretty(Gc, NULL, ctx); printf("\n");



    success = nmod_mpolyu_pgcd_zippel(Gbar, Abar, Bbar, var, ctx, zinfo, randstate);
    if (!success)
    {
        assert(0 && "contentless gcd failed");
        goto finished;
    }

printf("Gbar: "); nmod_mpolyu_print_pretty(Gbar, NULL, ctx); printf("\n");


    nmod_mpolyu_mul_mpoly(G, Gbar, Gc, ctx);

finished:

    nmod_mpoly_clear(Ac, ctx);
    nmod_mpoly_clear(Bc, ctx);
    nmod_mpoly_clear(Gc, ctx);
    nmod_mpolyu_clear(Abar, ctx);
    nmod_mpolyu_clear(Bbar, ctx);
    nmod_mpolyu_clear(Gbar, ctx);

    return success;
}




int nmod_mpoly_gcd_zippel_keepbits(nmod_mpoly_t G, nmod_mpoly_t A, nmod_mpoly_t B, const nmod_mpoly_ctx_t ctx)
{
    flint_rand_t randstate;
    int ret, success = 0;
    mpoly_zipinfo_t zinfo;
    nmod_mpoly_ctx_t uctx;
    nmod_mpolyu_t Au, Bu, Gu;
    const char * vars[] = {"x","y","z"};
//    TMP_INIT;

//    TMP_START;

    FLINT_ASSERT(A->bits <= FLINT_BITS);
    FLINT_ASSERT(B->bits <= FLINT_BITS);
    FLINT_ASSERT(A->bits == B->bits);

//flint_printf("**************\nnmod_mpoly_gcd_zippel: \nA(%wd): ", A->bits); nmod_mpoly_print_pretty(A, vars, ctx);
//flint_printf("\nB(%wd): ",B->bits); nmod_mpoly_print_pretty(B, vars, ctx); flint_printf("\n");

    FLINT_ASSERT(ctx->minfo->nvars > WORD(1));
    FLINT_ASSERT(!nmod_mpoly_is_zero(A, ctx));
    FLINT_ASSERT(!nmod_mpoly_is_zero(B, ctx));

    flint_randinit(randstate);

//flint_printf("1A bits: %d\n", A->bits);

    mpoly_zipinfo_init(zinfo, ctx->minfo->nvars);
    nmod_mpoly_degrees_si(zinfo->Adegs, A, ctx);
    nmod_mpoly_degrees_si(zinfo->Bdegs, B, ctx);
    for (slong i = 0; i < ctx->minfo->nvars; i++)
    {
        zinfo->perm[i] = i;
    }

//flint_printf("2A bits: %d\n", A->bits);


    slong new_bits = FLINT_MAX(A->bits, B->bits);


//flint_printf("new_bits: %d\n", new_bits);


    nmod_mpoly_ctx_init(uctx, ctx->minfo->nvars - 1, ORD_LEX, ctx->ffinfo->mod.n);
    nmod_mpolyu_init(Au, new_bits, uctx);
    nmod_mpolyu_init(Bu, new_bits, uctx);
    nmod_mpolyu_init(Gu, new_bits, uctx);

//flint_printf("3A bits: %d\n", A->bits);

    nmod_mpoly_to_mpolyu_perm(Au, A, zinfo->perm, uctx, ctx);
    nmod_mpoly_to_mpolyu_perm(Bu, B, zinfo->perm, uctx, ctx);

//flint_printf("Au: "); nmod_mpolyu_print_pretty(Au, NULL, uctx); flint_printf("\n");
//flint_printf("Bu: "); nmod_mpolyu_print_pretty(Bu, NULL, uctx); flint_printf("\n");
//flint_printf("4A bits: %d\n", A->bits);

    ret = nmod_mpolyu_pgcd_zippel(Gu, Au, Bu, uctx->minfo->nvars - 1, uctx, zinfo, randstate);
//flint_printf("5A bits: %d\n", A->bits);
    if (ret) {
        nmod_mpoly_from_mpolyu_keepbits(G, Gu, zinfo->perm, uctx, ctx);
        nmod_mpoly_make_monic(G, G, ctx);
        success = 1;
    }

//flint_printf("6A bits: %d\n", A->bits);

    nmod_mpolyu_clear(Au, uctx);
    nmod_mpolyu_clear(Bu, uctx);
    nmod_mpolyu_clear(Gu, uctx);
    nmod_mpoly_ctx_clear(uctx);

    mpoly_zipinfo_clear(zinfo);

    flint_randclear(randstate);

flint_printf("++++ A(%wd): ", A->bits); nmod_mpoly_print_pretty(A, NULL, ctx); printf("\n");
flint_printf("++++ B(%wd): ", B->bits); nmod_mpoly_print_pretty(B, NULL, ctx); printf("\n");
flint_printf("++++ G(%wd): ", G->bits); nmod_mpoly_print_pretty(G, NULL, ctx); printf("\n");

//    TMP_END;
    return success;
}

int nmod_mpoly_gcd_zippel(nmod_mpoly_t G, nmod_mpoly_t A, nmod_mpoly_t B, const nmod_mpoly_ctx_t ctx)
{
    flint_rand_t randstate;
    int ret, success = 0;
    mpoly_zipinfo_t zinfo;
    nmod_mpoly_ctx_t uctx;
    nmod_mpolyu_t Au, Bu, Gu;
    const char * vars[] = {"x","y","z"};
//    TMP_INIT;

//    TMP_START;

    FLINT_ASSERT(A->bits <= FLINT_BITS);
    FLINT_ASSERT(B->bits <= FLINT_BITS);

flint_printf("**************\nnmod_mpoly_gcd_zippel: \nA(%wd): ", A->bits); nmod_mpoly_print_pretty(A, vars, ctx);
flint_printf("\nB(%wd): ",B->bits); nmod_mpoly_print_pretty(B, vars, ctx); flint_printf("\n");

    FLINT_ASSERT(ctx->minfo->nvars > WORD(1));
    FLINT_ASSERT(!nmod_mpoly_is_zero(A, ctx));
    FLINT_ASSERT(!nmod_mpoly_is_zero(B, ctx));

    flint_randinit(randstate);

//flint_printf("1A bits: %d\n", A->bits);

    mpoly_zipinfo_init(zinfo, ctx->minfo->nvars);
    nmod_mpoly_degrees_si(zinfo->Adegs, A, ctx);
    nmod_mpoly_degrees_si(zinfo->Bdegs, B, ctx);
    for (slong i = 0; i < ctx->minfo->nvars; i++)
    {
        zinfo->perm[i] = ctx->minfo->nvars - 1 - i;
    }


    slong new_bits = FLINT_MAX(A->bits, B->bits);


flint_printf("new_bits: %d\n", new_bits);


    nmod_mpoly_ctx_init(uctx, ctx->minfo->nvars - 1, ORD_LEX, ctx->ffinfo->mod.n);
    nmod_mpolyu_init(Au, new_bits, uctx);
    nmod_mpolyu_init(Bu, new_bits, uctx);
    nmod_mpolyu_init(Gu, new_bits, uctx);


    nmod_mpoly_to_mpolyu_perm(Au, A, zinfo->perm, uctx, ctx);
    nmod_mpoly_to_mpolyu_perm(Bu, B, zinfo->perm, uctx, ctx);

//flint_printf("Au: "); nmod_mpolyu_print_pretty(Au, NULL, uctx); flint_printf("\n");
//flint_printf("Bu: "); nmod_mpolyu_print_pretty(Bu, NULL, uctx); flint_printf("\n");

    ret = nmod_mpolyu_pgcd_zippel(Gu, Au, Bu, uctx->minfo->nvars - 1, uctx, zinfo, randstate);
    if (ret) {
        nmod_mpoly_from_mpolyu_perm(G, Gu, zinfo->perm, uctx, ctx);
        nmod_mpoly_make_monic(G, G, ctx);
        success = 1;
    }


    nmod_mpolyu_clear(Au, uctx);
    nmod_mpolyu_clear(Bu, uctx);
    nmod_mpolyu_clear(Gu, uctx);
    nmod_mpoly_ctx_clear(uctx);

    mpoly_zipinfo_clear(zinfo);

    flint_randclear(randstate);

flint_printf("++++ A(%wd): ", A->bits); nmod_mpoly_print_pretty(A, NULL, ctx); printf("\n");
flint_printf("++++ B(%wd): ", B->bits); nmod_mpoly_print_pretty(B, NULL, ctx); printf("\n");
flint_printf("++++ G(%wd): ", G->bits); nmod_mpoly_print_pretty(G, NULL, ctx); printf("\n");

//    TMP_END;
    return success;
}











int nmod_mpolyu_gcd_zippel_linzipp(
    nmod_mpolyu_t G,
    nmod_mpolyu_t A,
    nmod_mpolyu_t B,
    slong var,
    nmod_mpoly_ctx_t ctx,
    mpoly_zipinfo_t zinfo,
    flint_rand_t randstate)
{
    int divcheck_fail_count, inner_gcd_added;
    slong lastdeg;
    ulong ABminshift;
    int success = 0, changed;

    FLINT_ASSERT(ctx->minfo->ord == ORD_LEX);
    FLINT_ASSERT(var >= -WORD(1));

    FLINT_ASSERT(G->bits == A->bits);
    FLINT_ASSERT(G->bits == B->bits);

flint_printf("**********nmod_mpolyu_gcd_zippel_linzipp var = %d *********\n", var);
flint_printf("A: "); nmod_mpolyu_print_pretty(A, NULL, ctx); flint_printf("\n");
flint_printf("B: "); nmod_mpolyu_print_pretty(B, NULL, ctx); flint_printf("\n");

    FLINT_ASSERT(A->length > 0);
    FLINT_ASSERT(B->length > 0);

    if (var == -WORD(1))
    {
        /* no more variables left to interpolate */
        return nmod_mpolyu_pgcd_zippel_univar(G, A, B, ctx);
    }

    if (var == WORD(0))
    {
        /* bivariate is more comfortable separated */
        return nmod_mpolyu_pgcd_zippel_bivar(G, A, B, ctx, zinfo);
    }

    nmod_mpolyun_t An, Bn;
    nmod_mpolyun_init(An, A->bits, ctx);
    nmod_mpolyun_init(Bn, A->bits, ctx);

    nmod_mpolyu_cvtto_mpolyun(An, A, var, ctx);
    nmod_mpolyu_cvtto_mpolyun(Bn, B, var, ctx);
    ABminshift = FLINT_MIN(A->exps[A->length - 1], B->exps[B->length - 1]);
    nmod_mpolyun_shift_right(An, ABminshift);
    nmod_mpolyun_shift_right(Bn, ABminshift);

flint_printf("An: "); nmod_mpolyun_print_pretty(An, NULL, ctx); flint_printf("\n");
flint_printf("Bn: "); nmod_mpolyun_print_pretty(Bn, NULL, ctx); flint_printf("\n");


    nmod_poly_t a, b, c, g;
    nmod_poly_init(a, ctx->ffinfo->mod.n);
    nmod_poly_init(b, ctx->ffinfo->mod.n);
    nmod_poly_init(c, ctx->ffinfo->mod.n);
    nmod_poly_init(g, ctx->ffinfo->mod.n);

    nmod_mpolyun_content_last(a, An, ctx);
    nmod_mpolyun_content_last(b, Bn, ctx);

flint_printf("a: "); nmod_poly_print_pretty(a, "v"); flint_printf("\n");
flint_printf("b: "); nmod_poly_print_pretty(b, "v"); flint_printf("\n");

//    nmod_mpolyun_divexact_last(An, a, ctx);
//    nmod_mpolyun_divexact_last(Bn, b, ctx);
//
//flint_printf("An: "); nmod_mpolyun_print_pretty(An, NULL, ctx); flint_printf("\n");
//flint_printf("Bn: "); nmod_mpolyun_print_pretty(Bn, NULL, ctx); flint_printf("\n");

    nmod_poly_gcd(c, a, b);

flint_printf("c: "); nmod_poly_print_pretty(c, "v"); flint_printf("\n");

    nmod_poly_gcd(g, nmod_mpolyun_leadcoeff_ref(An, ctx),
                     nmod_mpolyun_leadcoeff_ref(Bn, ctx));

flint_printf("g: "); nmod_poly_print_pretty(c, "v"); flint_printf("\n");

    nmod_poly_t modulus, tempmod;
    nmod_poly_init(modulus, ctx->ffinfo->mod.n);
    nmod_poly_init(tempmod, ctx->ffinfo->mod.n);
    nmod_poly_set_coeff_ui(tempmod, 1, UWORD(1));
    nmod_mpolyu_t Aeval, Beval, Geval, Gform;
    nmod_mpolyu_init(Aeval, A->bits, ctx);
    nmod_mpolyu_init(Beval, A->bits, ctx);
    nmod_mpolyu_init(Geval, A->bits, ctx);
    nmod_mpolyu_init(Gform, A->bits, ctx);
    nmod_mpolyun_t H, Ht;
    nmod_mpolyun_init(H, A->bits, ctx);
    nmod_mpolyun_init(Ht, A->bits, ctx);
    mp_limb_t geval, temp;

    if (nmod_poly_degree(c) > 0)
    {
        success = 0;
        goto finished;
    }

    nmod_poly_one(modulus);
    nmod_mpolyun_zero(H, ctx);

    mp_limb_t alpha = ctx->ffinfo->mod.n;


    while (1)
    {
        if (alpha == UWORD(0))
        {
            success = 0;
            goto finished;
        }
        alpha--;

flint_printf("------starting outer loop x%d = %d\n", var, alpha);


        geval = nmod_poly_evaluate_nmod(g, alpha);
        if (geval == WORD(0))
            goto outer_continue;

//flint_printf("Aeval bits: %d\n", Aeval->bits);
//flint_printf("   An bits: %d\n",    An->bits);
        nmod_mpolyun_eval_last(Aeval, An, alpha, ctx);

        nmod_mpolyun_eval_last(Beval, Bn, alpha, ctx);
        success = nmod_mpolyu_gcd_zippel_linzipp(Geval, Aeval, Beval, var - 1, ctx, zinfo, randstate);
flint_printf("********** pgcd return(%d) var = %d ", var-1, success); nmod_mpolyu_print_pretty(Geval, NULL, ctx); flint_printf("\n");
//flint_printf("Geval bits: %d\n", Geval->bits);
        if (!success)
            goto outer_continue;

        if (nmod_mpolyu_is_one(Geval, ctx))
        {
            nmod_mpolyu_cvtfrom_poly_notmain(G, c, var, ctx);
            return 1;
        }

        if (nmod_poly_degree(modulus) > 0)
        {
            if (Geval->exps[0] > H->exps[0])
            {
                goto outer_continue;

            } else if (Geval->exps[0] < H->exps[0])
            {
                nmod_poly_one(modulus);                
            }
        }

        /* update interpolant H */
flint_printf("------updating H\n");

        temp = n_invmod(nmod_mpolyu_leadcoeff(Geval, ctx), ctx->ffinfo->mod.n);
        nmod_mpolyu_scalar_mul_nmod(Geval, nmod_mul(geval, temp, ctx->ffinfo->mod), ctx);

        if (nmod_poly_degree(modulus) > 0)
        {
            nmod_poly_scalar_mul_nmod(modulus, modulus, n_invmod(nmod_poly_evaluate_nmod(modulus, alpha), ctx->ffinfo->mod.n));
            changed = nmod_mpolyun_addinterp(&lastdeg, H, Ht, Geval, modulus, alpha, ctx);
            if (!changed)
            {
                nmod_mpolyun_content_last(a, H, ctx);
                nmod_mpolyun_mul_poly(Ht, H, c, ctx);
                nmod_mpolyun_divexact_last(Ht, a, ctx);
                nmod_mpolyun_shift_left(Ht, ABminshift);
                nmod_mpolyu_cvtfrom_mpolyun(G, Ht, var, ctx);
                if (!nmod_mpolyu_divides(A, G, ctx))
                    goto outer_continue;
                if (!nmod_mpolyu_divides(B, G, ctx))
                    goto outer_continue;
                return 1;
            }

        } else
        {
            nmod_mpolyun_set_mpolyu(H, Geval, ctx);
        }
        nmod_poly_set_coeff_ui(tempmod, 0, ctx->ffinfo->mod.n - alpha);
        nmod_poly_mul(modulus, modulus, tempmod);

//flint_printf("after modulus: "); nmod_poly_print_pretty(modulus, "v"); printf("\n");


flint_printf("------H: "); nmod_mpolyun_print_pretty(H, NULL, ctx); flint_printf("\n");


        nmod_mpolyu_setform_mpolyun(Gform, H, ctx);

        divcheck_fail_count = 0;
        inner_gcd_added = 0;
        while (1)
        {
//flint_printf("asdf\n");
            if (alpha == UWORD(0))
            {
                success = 0;
                goto finished;
            }
            alpha--;

flint_printf("------starting inner loop x%d = %d\n", var, alpha);
        assert(alpha >= 80);

            geval = nmod_poly_evaluate_nmod(g, alpha);
            if (geval == WORD(0))
                goto inner_continue;

            nmod_mpolyun_eval_last(Aeval, An, alpha, ctx);
            nmod_mpolyun_eval_last(Beval, Bn, alpha, ctx);
            success = nmod_mpolyu_sgcd_zippel(Geval, Aeval, Beval, Gform, var, ctx, randstate);
flint_printf("sgcd return(%d) var = %d : ", success, var); nmod_mpolyu_print_pretty(Geval, NULL, ctx); flint_printf("\n");
//flint_printf("Geval bits: %d\n", Geval->bits);

            if (success == -1)
            {
                success = nmod_mpolyu_pgcd_zippel_rmcontent(G, A, B, var, ctx, zinfo, randstate);
                goto finished;
            }

            if (success != 1)
            {
                if (inner_gcd_added)
                    nmod_poly_one(modulus);
                goto outer_continue;
            }

            /* update interpolant H */
            temp = n_invmod(nmod_mpolyu_leadcoeff(Geval, ctx), ctx->ffinfo->mod.n);
            nmod_mpolyu_scalar_mul_nmod(Geval, nmod_mul(geval, temp, ctx->ffinfo->mod), ctx);
            FLINT_ASSERT(nmod_poly_degree(modulus) > 0);
            nmod_poly_scalar_mul_nmod(modulus, modulus, n_invmod(nmod_poly_evaluate_nmod(modulus, alpha), ctx->ffinfo->mod.n));
//flint_printf("Geval bits: %d\n", Geval->bits);
//flint_printf("before modulus: "); nmod_poly_print_pretty(modulus, "v"); printf("\n");
//flint_printf("before add     H: "); nmod_mpolyun_print_pretty(H, NULL, ctx); flint_printf("\n");
//flint_printf("before add Geval: "); nmod_mpolyu_print_pretty(Geval, NULL, ctx); flint_printf("\n");

            changed = nmod_mpolyun_addinterp(&lastdeg, H, Ht, Geval, modulus, alpha, ctx);
            nmod_poly_set_coeff_ui(tempmod, 0, ctx->ffinfo->mod.n - alpha);
            nmod_poly_mul(modulus, modulus, tempmod);
flint_printf("------updated H: "); nmod_mpolyun_print_pretty(H, NULL, ctx); flint_printf("\n");

            inner_gcd_added = changed;

            if (!changed)
            {
printf("didnt change\n");

                nmod_mpolyun_content_last(a, H, ctx);
                nmod_mpolyun_mul_poly(Ht, H, c, ctx);
                nmod_mpolyun_divexact_last(Ht, a, ctx);

                nmod_mpolyu_cvtfrom_mpolyun(G, Ht, var, ctx);

printf("candidate G: "); nmod_mpolyu_print_pretty(G, NULL, ctx); printf("\n");

                if (!nmod_mpolyu_divides(A, G, ctx)
                             || !nmod_mpolyu_divides(B, G, ctx))
                {
                    ++divcheck_fail_count;
                    flint_printf("division failed %d times\n", divcheck_fail_count);
                    if (divcheck_fail_count >= 2)
                    {
                        nmod_poly_one(modulus);
                        goto outer_continue;
                    }
                    goto inner_continue;
                }
                success = 1;
                goto finished;
            }

//flint_printf("after modulus: "); nmod_poly_print_pretty(modulus, "v"); printf("\n");

inner_continue:
            NULL;
        }
        FLINT_ASSERT(0 && "not reachable");

outer_continue:
        NULL;
    }
    FLINT_ASSERT(0 && "not reachable");


finished:

    nmod_poly_clear(a);
    nmod_poly_clear(b);
    nmod_poly_clear(c);
    nmod_poly_clear(g);
    nmod_poly_clear(modulus);
    nmod_poly_clear(tempmod);
    nmod_mpolyu_clear(Aeval, ctx);
    nmod_mpolyu_clear(Beval, ctx);
    nmod_mpolyu_clear(Geval, ctx);
    nmod_mpolyun_clear(An, ctx);
    nmod_mpolyun_clear(Bn, ctx);
    nmod_mpolyun_clear(H, ctx);
    nmod_mpolyun_clear(Ht, ctx);

    return success;
}
