/*
    Copyright (C) 2019-2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef FMPZ_MOD_MPOLY_H
#define FMPZ_MOD_MPOLY_H

#ifdef FMPZ_MOD_MPOLY_INLINES_C
#define FMPZ_MOD_MPOLY_INLINE FLINT_DLL
#else
#define FMPZ_MOD_MPOLY_INLINE static __inline__
#endif

#undef ulong
#define ulong ulongxx /* interferes with system includes */
#include <stdio.h>
#undef ulong

#include <gmp.h>
#define ulong mp_limb_t

#include "flint.h"
#include "fmpz_mod.h"
#include "fmpz_mpoly.h"
#include "mpoly.h"
#include "n_poly.h"

#ifdef __cplusplus
 extern "C" {
#endif

FLINT_DLL void fmpz_mod_ctx_init_rand_bits(fmpz_mod_ctx_t ctx,
                                  flint_rand_t state, flint_bitcnt_t max_bits);

FLINT_DLL void fmpz_mod_ctx_init_rand_bits_prime(fmpz_mod_ctx_t ctx,
                                  flint_rand_t state, flint_bitcnt_t max_bits);

FLINT_DLL void _fmpz_mod_vec_neg(fmpz * A, const fmpz * B, slong len,
                                                     const fmpz_mod_ctx_t ctx);

FLINT_DLL void _fmpz_mod_vec_scalar_mul_fmpz_mod(fmpz * A, const fmpz * B,
                          slong len, const fmpz_t c, const fmpz_mod_ctx_t ctx);

FLINT_DLL void _fmpz_mod_vec_scalar_div_fmpz_mod(fmpz * A, const fmpz * B,
                          slong len, const fmpz_t c, const fmpz_mod_ctx_t ctx);

FLINT_DLL void _fmpz_mod_vec_dot(fmpz_t d, const fmpz * a, const fmpz * b,
                                          slong len, const fmpz_mod_ctx_t ctx);

FMPZ_MOD_MPOLY_INLINE
void flint_mpz_add_uiui(mpz_ptr a, mpz_srcptr b, ulong c1, ulong c0)
{
    ulong d[2];
    mpz_t c;

    c->_mp_d = d;
    c->_mp_alloc = 2;

    d[1] = c1;
    d[0] = c0;

    c->_mp_size = d[1] != 0 ? 2 :
                  d[0] != 0;

    mpz_add(a, b, c);
}

FMPZ_MOD_MPOLY_INLINE
void flint_mpz_add_uiuiui(mpz_ptr a, mpz_srcptr b, ulong c2, ulong c1, ulong c0)
{
    ulong d[3];
    mpz_t c;

    c->_mp_d = d;
    c->_mp_alloc = 3;

    d[2] = c2;
    d[1] = c1;
    d[0] = c0;

    c->_mp_size = d[2] != 0 ? 3 :
                  d[1] != 0 ? 2 :
                  d[0] != 0;

    mpz_add(a, b, c);
}

typedef struct
{
    mpoly_ctx_t minfo;
    fmpz_mod_ctx_t ffinfo;
} fmpz_mod_mpoly_ctx_struct;

typedef fmpz_mod_mpoly_ctx_struct fmpz_mod_mpoly_ctx_t[1];

/*
    fmpz_mod_mpoly_t
    sparse multivariates with fmpz_mod coeffs
*/
typedef struct
{
    fmpz * coeffs;
    ulong * exps;
    slong length;
    flint_bitcnt_t bits;    /* number of bits per exponent */
    slong coeffs_alloc;     /* abs size in mp_limb_t units */
    slong exps_alloc;       /* abs size in ulong units */
} fmpz_mod_mpoly_struct;

typedef fmpz_mod_mpoly_struct fmpz_mod_mpoly_t[1];

/*
    fmpz_mod_mpoly_univar_t
    sparse univariates with multivariate coefficients
*/
typedef struct
{
   fmpz_mod_mpoly_struct * coeffs; /* multivariate coefficients */
   fmpz * exps;
   slong alloc;
   slong length;
} fmpz_mod_mpoly_univar_struct;

typedef fmpz_mod_mpoly_univar_struct fmpz_mod_mpoly_univar_t[1];



/* Internal type definitions *************************************************/

typedef struct
{
    fmpz_mod_poly_struct * coeffs;
    slong alloc;
    slong length;
} fmpz_mod_bpoly_struct;

typedef fmpz_mod_bpoly_struct fmpz_mod_bpoly_t[1];


typedef struct
{
    fmpz_mod_bpoly_struct * coeffs;
    slong alloc;
    slong length;
} fmpz_mod_tpoly_struct;

typedef fmpz_mod_tpoly_struct fmpz_mod_tpoly_t[1];


typedef struct
{
    ulong exp;
    fmpz_t coeff;
} fmpz_mod_polyu_term_struct;

typedef struct
{
    ulong * exps;
    fmpz * coeffs;
    slong length;
    slong alloc;
} fmpz_mod_polyu_struct;

typedef fmpz_mod_polyu_struct fmpz_mod_polyu_t[1];


typedef struct
{
    ulong exp;
    fmpz_mod_poly_t coeff;
} fmpz_mod_polyun_term_struct;

typedef struct
{
    fmpz_mod_polyun_term_struct * terms;
    slong length;
    slong alloc;
} fmpz_mod_polyun_struct;

typedef fmpz_mod_polyun_struct fmpz_mod_polyun_t[1];

/*
    fmpz_mod_mpolyu_t
    sparse univariates with fmpz_mod_mpoly_t coefficients
        with uniform bits and LEX ordering
*/
typedef struct
{
   fmpz_mod_mpoly_struct * coeffs;
   ulong * exps;
   slong alloc;
   slong length;
   flint_bitcnt_t bits;    /* default bits to construct coeffs */
} fmpz_mod_mpolyu_struct;
typedef fmpz_mod_mpolyu_struct fmpz_mod_mpolyu_t[1];

/*
    fmpz_mod_mpolyn_t
    sparse multivariates with fmpz_mod_poly_t coefficients
        with LEX ordering
*/
typedef struct
{
   fmpz_mod_poly_struct * coeffs;
   ulong * exps;
   slong alloc;
   slong length;
   slong bits;
} fmpz_mod_mpolyn_struct;

typedef fmpz_mod_mpolyn_struct fmpz_mod_mpolyn_t[1];

typedef struct
{
    fmpz_mod_mpolyn_struct * coeffs;
    ulong * exps;
    slong alloc;
    slong length;
    flint_bitcnt_t bits;
} fmpz_mod_mpolyun_struct;

typedef fmpz_mod_mpolyun_struct fmpz_mod_mpolyun_t[1];


typedef struct
{
    fmpz_mod_poly_struct ** array;
    slong alloc;
    slong top;
} fmpz_mod_poly_stack_struct;

typedef fmpz_mod_poly_stack_struct fmpz_mod_poly_stack_t[1];


typedef struct
{
    fmpz_mod_bpoly_struct ** array;
    slong alloc;
    slong top;
} fmpz_mod_bpoly_stack_struct;

typedef fmpz_mod_bpoly_stack_struct fmpz_mod_bpoly_stack_t[1];

typedef struct
{
    fmpz_mod_polyun_struct ** array;
    slong alloc;
    slong top;
} fmpz_mod_polyun_stack_struct;

typedef fmpz_mod_polyun_stack_struct fmpz_mod_polyun_stack_t[1];

typedef struct
{
    fmpz_mod_mpolyn_struct ** array;
    slong alloc;
    slong top;
    flint_bitcnt_t bits;
} fmpz_mod_mpolyn_stack_struct;

typedef fmpz_mod_mpolyn_stack_struct fmpz_mod_mpolyn_stack_t[1];

typedef struct {
    fmpz_mod_poly_stack_t poly_stack;
    fmpz_mod_bpoly_stack_t bpoly_stack;
} fmpz_mod_poly_bpoly_stack_struct;

typedef fmpz_mod_poly_bpoly_stack_struct fmpz_mod_poly_bpoly_stack_t[1];

typedef struct {
    fmpz_mod_poly_stack_t poly_stack;
    fmpz_mod_polyun_stack_t polyun_stack;
} fmpz_mod_poly_polyun_stack_struct;

typedef fmpz_mod_poly_polyun_stack_struct fmpz_mod_poly_polyun_stack_t[1];

typedef struct {
    fmpz_mod_poly_stack_t poly_stack;
    fmpz_mod_polyun_stack_t polyun_stack;
    fmpz_mod_mpolyn_stack_t mpolyn_stack;
} fmpz_mod_poly_polyun_mpolyn_stack_struct;

typedef fmpz_mod_poly_polyun_mpolyn_stack_struct fmpz_mod_poly_polyun_mpolyn_stack_t[1];


/*****************************************************************************/

FLINT_DLL void fmpz_mod_poly_stack_init(fmpz_mod_poly_stack_t S);

FLINT_DLL void fmpz_mod_poly_stack_clear(fmpz_mod_poly_stack_t S);

FLINT_DLL fmpz_mod_poly_struct ** fmpz_mod_poly_stack_fit_request(fmpz_mod_poly_stack_t S, slong k);

FMPZ_MOD_MPOLY_INLINE
fmpz_mod_poly_struct ** fmpz_mod_poly_stack_request(fmpz_mod_poly_stack_t S, slong k)
{
    fmpz_mod_poly_struct ** poly_top;
    poly_top = fmpz_mod_poly_stack_fit_request(S, k);
    S->top += k;
    return poly_top;
}

FMPZ_MOD_MPOLY_INLINE
fmpz_mod_poly_struct * fmpz_mod_poly_stack_take_top(fmpz_mod_poly_stack_t S)
{
    /* assume the request for 1 has already been fitted */
    fmpz_mod_poly_struct ** poly_top;
    FLINT_ASSERT(S->top + 1 <= S->alloc);
    poly_top = S->array + S->top;
    S->top += 1;
    return poly_top[0];
}

FMPZ_MOD_MPOLY_INLINE
void fmpz_mod_poly_stack_give_back(fmpz_mod_poly_stack_t S, slong k)
{
    FLINT_ASSERT(S->top >= k);
    S->top -= k;
}

FMPZ_MOD_MPOLY_INLINE
slong fmpz_mod_poly_stack_size(const fmpz_mod_poly_stack_t S)
{
    return S->top;
}


FLINT_DLL void fmpz_mod_bpoly_stack_init(fmpz_mod_bpoly_stack_t S);

FLINT_DLL void fmpz_mod_bpoly_stack_clear(fmpz_mod_bpoly_stack_t S);

FLINT_DLL fmpz_mod_bpoly_struct ** fmpz_mod_bpoly_stack_fit_request(fmpz_mod_bpoly_stack_t S, slong k);

FMPZ_MOD_MPOLY_INLINE
fmpz_mod_bpoly_struct ** fmpz_mod_bpoly_stack_request(fmpz_mod_bpoly_stack_t S, slong k)
{
    fmpz_mod_bpoly_struct ** bpoly_top;
    bpoly_top = fmpz_mod_bpoly_stack_fit_request(S, k);
    S->top += k;
    return bpoly_top;
}

FMPZ_MOD_MPOLY_INLINE
fmpz_mod_bpoly_struct * fmpz_mod_bpoly_stack_take_top(fmpz_mod_bpoly_stack_t S)
{
    /* assume the request for 1 has already been fitted */
    fmpz_mod_bpoly_struct ** bpoly_top;
    FLINT_ASSERT(S->top + 1 <= S->alloc);
    bpoly_top = S->array + S->top;
    S->top += 1;
    return bpoly_top[0];
}

FMPZ_MOD_MPOLY_INLINE
void fmpz_mod_bpoly_stack_give_back(fmpz_mod_bpoly_stack_t S, slong k)
{
    FLINT_ASSERT(S->top >= k);
    S->top -= k;
}

FMPZ_MOD_MPOLY_INLINE
slong fmpz_mod_bpoly_stack_size(const fmpz_mod_bpoly_stack_t S)
{
    return S->top;
}


FLINT_DLL void fmpz_mod_polyun_stack_init(fmpz_mod_polyun_stack_t S);

FLINT_DLL void fmpz_mod_polyun_stack_clear(fmpz_mod_polyun_stack_t S);

FLINT_DLL fmpz_mod_polyun_struct ** fmpz_mod_polyun_stack_fit_request(fmpz_mod_polyun_stack_t S, slong k);

FMPZ_MOD_MPOLY_INLINE
fmpz_mod_polyun_struct ** fmpz_mod_polyun_stack_request(fmpz_mod_polyun_stack_t S, slong k)
{
    fmpz_mod_polyun_struct ** polyun_top;
    polyun_top = fmpz_mod_polyun_stack_fit_request(S, k);
    S->top += k;
    return polyun_top;
}

FMPZ_MOD_MPOLY_INLINE
fmpz_mod_polyun_struct * fmpz_mod_polyun_stack_take_top(fmpz_mod_polyun_stack_t S)
{
    /* assume the request for 1 has already been fitted */
    fmpz_mod_polyun_struct ** polyun_top;
    FLINT_ASSERT(S->top + 1 <= S->alloc);
    polyun_top = S->array + S->top;
    S->top += 1;
    return polyun_top[0];
}

FMPZ_MOD_MPOLY_INLINE
void fmpz_mod_polyun_stack_give_back(fmpz_mod_polyun_stack_t S, slong k)
{
    FLINT_ASSERT(S->top >= k);
    S->top -= k;
}

FMPZ_MOD_MPOLY_INLINE
slong fmpz_mod_polyun_stack_size(const fmpz_mod_polyun_stack_t S)
{
    return S->top;
}



FLINT_DLL void fmpz_mod_mpolyn_stack_init(fmpz_mod_mpolyn_stack_t S,
                          flint_bitcnt_t bits, const fmpz_mod_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mod_mpolyn_stack_clear(fmpz_mod_mpolyn_stack_t S,
                                               const fmpz_mod_mpoly_ctx_t ctx);

FLINT_DLL fmpz_mod_mpolyn_struct ** fmpz_mod_mpolyn_stack_fit_request(
           fmpz_mod_mpolyn_stack_t S, slong k, const fmpz_mod_mpoly_ctx_t ctx);

FMPZ_MOD_MPOLY_INLINE
fmpz_mod_mpolyn_struct ** fmpz_mod_mpolyn_stack_request(
            fmpz_mod_mpolyn_stack_t S, slong k, const fmpz_mod_mpoly_ctx_t ctx)
{
    fmpz_mod_mpolyn_struct ** mpolyn_top;
    mpolyn_top = fmpz_mod_mpolyn_stack_fit_request(S, k, ctx);
    S->top += k;
    return mpolyn_top;
}

FMPZ_MOD_MPOLY_INLINE
fmpz_mod_mpolyn_struct * fmpz_mod_mpolyn_stack_take_top(fmpz_mod_mpolyn_stack_t S)
{
    /* assume the request for 1 has already been fitted */
    fmpz_mod_mpolyn_struct ** mpolyn_top;
    FLINT_ASSERT(S->top + 1 <= S->alloc);
    mpolyn_top = S->array + S->top;
    S->top += 1;
    return mpolyn_top[0];
}

FMPZ_MOD_MPOLY_INLINE
void fmpz_mod_mpolyn_stack_give_back(fmpz_mod_mpolyn_stack_t S, slong k)
{
    FLINT_ASSERT(S->top >= k);
    S->top -= k;
}

FMPZ_MOD_MPOLY_INLINE
slong fmpz_mod_mpolyn_stack_size(const fmpz_mod_mpolyn_stack_t S)
{
    return S->top;
}




/* Context object ************************************************************/

FLINT_DLL void fmpz_mod_mpoly_ctx_init(fmpz_mod_mpoly_ctx_t ctx, 
                      slong nvars, const ordering_t ord, const fmpz_t modulus);

FLINT_DLL void fmpz_mod_mpoly_ctx_init_rand(fmpz_mod_mpoly_ctx_t ctx,
                    flint_rand_t state, slong max_nvars, const fmpz_t modulus);

FLINT_DLL void fmpz_mod_mpoly_ctx_init_rand_bits_prime(fmpz_mod_mpoly_ctx_t ctx,
                 flint_rand_t state, slong max_nvars, flint_bitcnt_t max_bits);

FLINT_DLL void fmpz_mod_mpoly_ctx_init_rand_bits(fmpz_mod_mpoly_ctx_t ctx,
                 flint_rand_t state, slong max_nvars, flint_bitcnt_t max_bits);

FLINT_DLL void fmpz_mod_mpoly_ctx_clear(fmpz_mod_mpoly_ctx_t ctx);

FMPZ_MOD_MPOLY_INLINE
slong fmpz_mod_mpoly_ctx_nvars(const fmpz_mod_mpoly_ctx_t ctx)
{
    return ctx->minfo->nvars;
}

FMPZ_MOD_MPOLY_INLINE
ordering_t fmpz_mod_mpoly_ctx_ord(const fmpz_mod_mpoly_ctx_t ctx)
{
    return ctx->minfo->ord;
}

FMPZ_MOD_MPOLY_INLINE
const fmpz * fmpz_mod_mpoly_ctx_modulus(const fmpz_mod_mpoly_ctx_t ctx)
{
    return fmpz_mod_ctx_modulus(ctx->ffinfo);
}

FMPZ_MOD_MPOLY_INLINE
void fmpz_mod_mpoly_ctx_get_modulus(fmpz_t m, const fmpz_mod_mpoly_ctx_t ctx)
{
    fmpz_set(m, fmpz_mod_mpoly_ctx_modulus(ctx));
}

/*  Memory management ********************************************************/


FMPZ_MOD_MPOLY_INLINE
void fmpz_mod_mpoly_init(fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_ctx_t ctx)
{
    A->coeffs = NULL;
    A->exps = NULL;
    A->length = 0;
    A->bits = MPOLY_MIN_BITS;
    A->coeffs_alloc = 0;
    A->exps_alloc = 0;
}

FLINT_DLL void fmpz_mod_mpoly_clear(fmpz_mod_mpoly_t A,
                                               const fmpz_mod_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mod_mpoly_init2(fmpz_mod_mpoly_t A, slong alloc,
                                               const fmpz_mod_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mod_mpoly_init3(fmpz_mod_mpoly_t A, slong alloc,
                          flint_bitcnt_t bits, const fmpz_mod_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mod_mpoly_realloc(fmpz_mod_mpoly_t A,
                                  slong alloc, const fmpz_mod_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mod_mpoly_fit_length(fmpz_mod_mpoly_t A, slong length,
                                               const fmpz_mod_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mod_mpoly_fit_length_fit_bits(fmpz_mod_mpoly_t A,
               slong len, flint_bitcnt_t bits, const fmpz_mod_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mod_mpoly_fit_length_reset_bits(fmpz_mod_mpoly_t A,
               slong len, flint_bitcnt_t bits, const fmpz_mod_mpoly_ctx_t ctx);

FMPZ_MOD_MPOLY_INLINE
void _fmpz_mod_mpoly_fit_length(
    fmpz ** coeffs,
    slong * coeffs_alloc,
    ulong ** exps,
    slong * exps_alloc,
    slong N,
    slong length)
{
    if (length > *coeffs_alloc)
    {
        slong i, old_alloc = *coeffs_alloc;
        slong new_alloc = FLINT_MAX(length, old_alloc*2);
        *coeffs_alloc = new_alloc;
        *coeffs = (fmpz *) flint_realloc(*coeffs, new_alloc*sizeof(fmpz));
        for (i = old_alloc; i < new_alloc; i++)
            fmpz_init(*coeffs + i);
    }

    if (N*length > *exps_alloc)
    {
        *exps_alloc = FLINT_MAX(N*length, *exps_alloc*2);
        *exps = (ulong *) flint_realloc(*exps, *exps_alloc*sizeof(ulong));
    }
}

FMPZ_MOD_MPOLY_INLINE
void _fmpz_mod_mpoly_set_length(fmpz_mod_mpoly_t A, slong newlen, 
                                                const fmpz_mod_mpoly_ctx_t ctx)
{
/*
    slong i;
*/

    FLINT_ASSERT(newlen <= A->coeffs_alloc);
    FLINT_ASSERT(mpoly_words_per_exp(A->bits, ctx->minfo)*newlen <= A->exps_alloc);
/*
    for (i = A->length - 1; i >= newlen; i--)
       _fmpz_demote(A->coeffs + i);
*/
    A->length = newlen;
}

FMPZ_MOD_MPOLY_INLINE
void fmpz_mod_mpoly_truncate(fmpz_mod_mpoly_t A, slong newlen, 
                                                const fmpz_mod_mpoly_ctx_t ctx)
{
    if (A->length > newlen)
    {
        slong i;

        for (i = newlen; i < A->length; i++)
            _fmpz_demote(A->coeffs + i);

        A->length = newlen;
    }
}


/* Input/output **************************************************************/

FLINT_DLL int fmpz_mod_mpoly_set_str_pretty(fmpz_mod_mpoly_t A, const char * str,
                                  const char ** x, const fmpz_mod_mpoly_ctx_t ctx);

FLINT_DLL char * fmpz_mod_mpoly_get_str_pretty(const fmpz_mod_mpoly_t A,
                                  const char ** x, const fmpz_mod_mpoly_ctx_t ctx);

FLINT_DLL int fmpz_mod_mpoly_fprint_pretty(FILE * file, 
            const fmpz_mod_mpoly_t A, const char ** x, const fmpz_mod_mpoly_ctx_t ctx);

FMPZ_MOD_MPOLY_INLINE
int fmpz_mod_mpoly_print_pretty(const fmpz_mod_mpoly_t A,
                                   const char ** x, const fmpz_mod_mpoly_ctx_t ctx)
{
   return fmpz_mod_mpoly_fprint_pretty(stdout, A, x, ctx);
}


/*  Basic manipulation *******************************************************/

FLINT_DLL void fmpz_mod_mpoly_gen(fmpz_mod_mpoly_t A, slong var,
                                               const fmpz_mod_mpoly_ctx_t ctx);

FLINT_DLL int fmpz_mod_mpoly_is_gen(const fmpz_mod_mpoly_t A,
                                    slong var, const fmpz_mod_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mod_mpoly_set(fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_t B,
                                               const fmpz_mod_mpoly_ctx_t ctx);

FLINT_DLL int fmpz_mod_mpoly_equal(const fmpz_mod_mpoly_t A,
                     const fmpz_mod_mpoly_t B, const fmpz_mod_mpoly_ctx_t ctx);

FMPZ_MOD_MPOLY_INLINE
void fmpz_mod_mpoly_swap(fmpz_mod_mpoly_t A, fmpz_mod_mpoly_t B,
                                                const fmpz_mod_mpoly_ctx_t ctx)
{
    fmpz_mod_mpoly_struct t = *A;
    *A = *B;
    *B = t;
}

/* Constants *****************************************************************/

FLINT_DLL int fmpz_mod_mpoly_is_fmpz(const fmpz_mod_mpoly_t A,
                                               const fmpz_mod_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mod_mpoly_get_fmpz(fmpz_t c, const fmpz_mod_mpoly_t A,
                                               const fmpz_mod_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mod_mpoly_set_fmpz_mod(fmpz_mod_mpoly_t A,
                               const fmpz_t c, const fmpz_mod_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mod_mpoly_set_fmpz(fmpz_mod_mpoly_t A,
                               const fmpz_t c, const fmpz_mod_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mod_mpoly_set_si(fmpz_mod_mpoly_t A,
                                      slong c, const fmpz_mod_mpoly_ctx_t ctx);

FMPZ_MOD_MPOLY_INLINE
void fmpz_mod_mpoly_zero(fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_ctx_t ctx)
{
   _fmpz_mod_mpoly_set_length(A, 0, ctx);
}

FMPZ_MOD_MPOLY_INLINE
void fmpz_mod_mpoly_one(fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_ctx_t ctx)
{
    fmpz_mod_mpoly_set_si(A, 1, ctx);
}

FLINT_DLL int fmpz_mod_mpoly_equal_fmpz(const fmpz_mod_mpoly_t A,
                               const fmpz_t c, const fmpz_mod_mpoly_ctx_t ctx);

FLINT_DLL int fmpz_mod_mpoly_equal_si(const fmpz_mod_mpoly_t A,
                                      slong c, const fmpz_mod_mpoly_ctx_t ctx);

FMPZ_MOD_MPOLY_INLINE
int fmpz_mod_mpoly_is_zero(const fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_ctx_t ctx)
{
   return A->length < 1;
}

FMPZ_MOD_MPOLY_INLINE
int fmpz_mod_mpoly_is_one(const fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_ctx_t ctx)
{
   return fmpz_mod_mpoly_equal_si(A, 1, ctx);
}

/* Degrees *******************************************************************/

FMPZ_MOD_MPOLY_INLINE
int fmpz_mod_mpoly_degrees_fit_si(const fmpz_mod_mpoly_t A,
                                                const fmpz_mod_mpoly_ctx_t ctx)
{
    return A->bits <= FLINT_BITS ? 1
               : mpoly_degrees_fit_si(A->exps, A->length, A->bits, ctx->minfo);
}

FMPZ_MOD_MPOLY_INLINE
void fmpz_mod_mpoly_degrees_fmpz(fmpz ** degs, const fmpz_mod_mpoly_t A,
                                                const fmpz_mod_mpoly_ctx_t ctx)
{
    mpoly_degrees_pfmpz(degs, A->exps, A->length, A->bits, ctx->minfo);
}

FMPZ_MOD_MPOLY_INLINE
void fmpz_mod_mpoly_degrees_si(slong * degs, const fmpz_mod_mpoly_t A,
                                                const fmpz_mod_mpoly_ctx_t ctx)
{
    mpoly_degrees_si(degs, A->exps, A->length, A->bits, ctx->minfo);
}

FMPZ_MOD_MPOLY_INLINE
void fmpz_mod_mpoly_degree_fmpz(fmpz_t deg, const fmpz_mod_mpoly_t A, slong var,
                                                const fmpz_mod_mpoly_ctx_t ctx)
{
    mpoly_degree_fmpz(deg, A->exps, A->length, A->bits, var, ctx->minfo);
}

FMPZ_MOD_MPOLY_INLINE
slong fmpz_mod_mpoly_degree_si(const fmpz_mod_mpoly_t A, slong var,
                                                const fmpz_mod_mpoly_ctx_t ctx)
{
    return mpoly_degree_si(A->exps, A->length, A->bits, var, ctx->minfo);
}

FMPZ_MOD_MPOLY_INLINE
int fmpz_mod_mpoly_total_degree_fits_si(const fmpz_mod_mpoly_t A,
                                                const fmpz_mod_mpoly_ctx_t ctx)
{
    return mpoly_total_degree_fits_si(A->exps, A->length, A->bits, ctx->minfo);
}

FMPZ_MOD_MPOLY_INLINE
void fmpz_mod_mpoly_total_degree_fmpz(fmpz_t td, const fmpz_mod_mpoly_t A,
                                                const fmpz_mod_mpoly_ctx_t ctx)
{
    mpoly_total_degree_fmpz(td, A->exps, A->length, A->bits, ctx->minfo);
}

FMPZ_MOD_MPOLY_INLINE
slong fmpz_mod_mpoly_total_degree_si(const fmpz_mod_mpoly_t A,
                                                const fmpz_mod_mpoly_ctx_t ctx)
{
    return mpoly_total_degree_si(A->exps, A->length, A->bits, ctx->minfo);
}


/* Coefficients **************************************************************/

FLINT_DLL void fmpz_mod_mpoly_get_coeff_fmpz_monomial(fmpz_t c,
                        const fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_t M,
                                               const fmpz_mod_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mod_mpoly_set_coeff_fmpz_monomial(fmpz_mod_mpoly_t A,
                                    const fmpz_t c, const fmpz_mod_mpoly_t M,
                                               const fmpz_mod_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mod_mpoly_get_coeff_fmpz_fmpz(fmpz_t c,
                                const fmpz_mod_mpoly_t A, fmpz * const * exp,
                                               const fmpz_mod_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mod_mpoly_get_coeff_fmpz_ui(fmpz_t c,
                                const fmpz_mod_mpoly_t A, const ulong * exp,
                                               const fmpz_mod_mpoly_ctx_t ctx);

FLINT_DLL void _fmpz_mod_mpoly_set_coeff_fmpz_fmpz(fmpz_mod_mpoly_t A,
             const fmpz_t c, const fmpz * exp, const fmpz_mod_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mod_mpoly_set_coeff_fmpz_fmpz(fmpz_mod_mpoly_t A,
           const fmpz_t c, fmpz * const * exp, const fmpz_mod_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mod_mpoly_set_coeff_fmpz_ui(fmpz_mod_mpoly_t A,
            const fmpz_t c, const ulong * exp, const fmpz_mod_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mod_mpoly_get_coeff_vars_ui(fmpz_mod_mpoly_t C,
             const fmpz_mod_mpoly_t A, const slong * vars, const ulong * exps,
                                 slong length, const fmpz_mod_mpoly_ctx_t ctx);

FMPZ_MOD_MPOLY_INLINE fmpz * fmpz_mod_mpoly_leadcoeff(fmpz_mod_mpoly_t A)
{
    FLINT_ASSERT(A->length > 0);
    return A->coeffs + 0;
}

/* comparison ****************************************************************/

FLINT_DLL int fmpz_mod_mpoly_cmp(const fmpz_mod_mpoly_t A,
                     const fmpz_mod_mpoly_t B, const fmpz_mod_mpoly_ctx_t ctx);


/* container operations ******************************************************/

FLINT_DLL int fmpz_mod_mpoly_is_canonical(const fmpz_mod_mpoly_t A,
                                               const fmpz_mod_mpoly_ctx_t ctx);

FMPZ_MOD_MPOLY_INLINE
slong fmpz_mod_mpoly_length(const fmpz_mod_mpoly_t A,
                                                const fmpz_mod_mpoly_ctx_t ctx)
{
    return A->length;
}

FLINT_DLL void fmpz_mod_mpoly_resize(fmpz_mod_mpoly_t A, slong new_length,
                                               const fmpz_mod_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mod_mpoly_get_term_coeff_fmpz(fmpz_t c,
            const fmpz_mod_mpoly_t A, slong i, const fmpz_mod_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mod_mpoly_set_term_coeff_fmpz(fmpz_mod_mpoly_t A, slong i,
                               const fmpz_t c, const fmpz_mod_mpoly_ctx_t ctx);

FMPZ_MOD_MPOLY_INLINE
int fmpz_mod_mpoly_term_exp_fits_ui(const fmpz_mod_mpoly_t A, slong i,
                                                const fmpz_mod_mpoly_ctx_t ctx)
{
    return A->bits <= FLINT_BITS ? 1
                     : mpoly_term_exp_fits_ui(A->exps, A->bits, i, ctx->minfo);
}

FMPZ_MOD_MPOLY_INLINE
int fmpz_mod_mpoly_term_exp_fits_si(const fmpz_mod_mpoly_t A, slong i,
                                                const fmpz_mod_mpoly_ctx_t ctx)
{
    return A->bits <= FLINT_BITS ? 1
                     : mpoly_term_exp_fits_si(A->exps, A->bits, i, ctx->minfo);
}

FLINT_DLL void fmpz_mod_mpoly_get_term_exp_fmpz(fmpz ** exp,
            const fmpz_mod_mpoly_t A, slong i, const fmpz_mod_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mod_mpoly_get_term_exp_ui(ulong * exp,
            const fmpz_mod_mpoly_t A, slong i, const fmpz_mod_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mod_mpoly_get_term_exp_si(slong * exp,
            const fmpz_mod_mpoly_t A, slong i, const fmpz_mod_mpoly_ctx_t ctx);

FLINT_DLL ulong fmpz_mod_mpoly_get_term_var_exp_ui(const fmpz_mod_mpoly_t A,
                           slong i, slong var, const fmpz_mod_mpoly_ctx_t ctx);

FLINT_DLL slong fmpz_mod_mpoly_get_term_var_exp_si(const fmpz_mod_mpoly_t A,
                           slong i, slong var, const fmpz_mod_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mod_mpoly_set_term_exp_fmpz(fmpz_mod_mpoly_t A, slong i,
                           fmpz * const * exp, const fmpz_mod_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mod_mpoly_set_term_exp_ui(fmpz_mod_mpoly_t A, slong i,
                            const ulong * exp, const fmpz_mod_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mod_mpoly_get_term(fmpz_mod_mpoly_t M,
            const fmpz_mod_mpoly_t A, slong i, const fmpz_mod_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mod_mpoly_get_term_monomial(fmpz_mod_mpoly_t M,
            const fmpz_mod_mpoly_t A, slong i, const fmpz_mod_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mod_mpoly_push_term_fmpz_fmpz(fmpz_mod_mpoly_t A,
           const fmpz_t c, fmpz * const * exp, const fmpz_mod_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mod_mpoly_push_term_fmpz_ui(fmpz_mod_mpoly_t A,
           const fmpz_t c, const ulong * exp, const fmpz_mod_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mod_mpoly_sort_terms(fmpz_mod_mpoly_t A,
                                               const fmpz_mod_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mod_mpoly_combine_like_terms(fmpz_mod_mpoly_t A,
                                               const fmpz_mod_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mod_mpoly_reverse(fmpz_mod_mpoly_t A,
                     const fmpz_mod_mpoly_t B, const fmpz_mod_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mod_mpoly_assert_canonical(const fmpz_mod_mpoly_t A,
                                               const fmpz_mod_mpoly_ctx_t ctx);

FLINT_DLL void _fmpz_mod_mpoly_radix_sort1(fmpz_mod_mpoly_t A, slong left,
              slong right, flint_bitcnt_t pos, ulong cmpmask, ulong totalmask);

FLINT_DLL void _fmpz_mod_mpoly_radix_sort(fmpz_mod_mpoly_t A, slong left,
                    slong right, flint_bitcnt_t pos, slong N, ulong * cmpmask);

FLINT_DLL void _fmpz_mod_mpoly_push_exp_ffmpz(fmpz_mod_mpoly_t A,
                             const fmpz * exp, const fmpz_mod_mpoly_ctx_t ctx);

FLINT_DLL void _fmpz_mod_mpoly_push_exp_pfmpz(fmpz_mod_mpoly_t A,
                           fmpz * const * exp, const fmpz_mod_mpoly_ctx_t ctx);

FLINT_DLL void _fmpz_mod_mpoly_push_exp_ui(fmpz_mod_mpoly_t A,
                            const ulong * exp, const fmpz_mod_mpoly_ctx_t ctx);


/* Random generation *********************************************************/

FLINT_DLL void fmpz_mod_mpoly_randtest_bounds(fmpz_mod_mpoly_t A,
                        flint_rand_t state, slong length, ulong * exp_bounds,
                                               const fmpz_mod_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mod_mpoly_randtest_bound(fmpz_mod_mpoly_t A,
                            flint_rand_t state, slong length, ulong exp_bound,
                                               const fmpz_mod_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mod_mpoly_randtest_bits(fmpz_mod_mpoly_t A,
                    flint_rand_t state, slong length, flint_bitcnt_t exp_bits,
                                               const fmpz_mod_mpoly_ctx_t ctx);


/* Addition/Subtraction ******************************************************/

FLINT_DLL void fmpz_mod_mpoly_add_fmpz_mod(fmpz_mod_mpoly_t A,
                                    const fmpz_mod_mpoly_t B, const fmpz_t c,
                                               const fmpz_mod_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mod_mpoly_add_fmpz(fmpz_mod_mpoly_t A,
                                    const fmpz_mod_mpoly_t B, const fmpz_t c,
                                               const fmpz_mod_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mod_mpoly_add_si(fmpz_mod_mpoly_t A,
                                    const fmpz_mod_mpoly_t B, slong c,
                                               const fmpz_mod_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mod_mpoly_sub_fmpz(fmpz_mod_mpoly_t A,
                                    const fmpz_mod_mpoly_t B, const fmpz_t c,
                                               const fmpz_mod_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mod_mpoly_sub_si(fmpz_mod_mpoly_t A,
                                    const fmpz_mod_mpoly_t B, slong c,
                                               const fmpz_mod_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mod_mpoly_add(fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_t B,
                     const fmpz_mod_mpoly_t C, const fmpz_mod_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mod_mpoly_sub(fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_t B,
                     const fmpz_mod_mpoly_t C, const fmpz_mod_mpoly_ctx_t ctx);

/* Scalar operations *********************************************************/

FLINT_DLL void fmpz_mod_mpoly_neg(fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_t B,
                                               const fmpz_mod_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mod_mpoly_scalar_mul_fmpz_mod_invertible(fmpz_mod_mpoly_t A,
                                    const fmpz_mod_mpoly_t B, const fmpz_t c,
                                               const fmpz_mod_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mod_mpoly_scalar_mul_fmpz(fmpz_mod_mpoly_t A,
                                    const fmpz_mod_mpoly_t B, const fmpz_t c,
                                               const fmpz_mod_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mod_mpoly_make_monic(fmpz_mod_mpoly_t A,
                     const fmpz_mod_mpoly_t B, const fmpz_mod_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mod_mpoly_scalar_mul_fmpz_mod_invertible(fmpz_mod_mpoly_t A,
                                      const fmpz_mod_mpoly_t B, const fmpz_t c,
                                               const fmpz_mod_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mod_mpoly_scalar_addmul_fmpz(fmpz_mod_mpoly_t A,
                        const fmpz_mod_mpoly_t B, const fmpz_mod_mpoly_t C,
                               const fmpz_t d, const fmpz_mod_mpoly_ctx_t ctx);

/* Differention **************************************************************/

FLINT_DLL void fmpz_mod_mpoly_derivative(fmpz_mod_mpoly_t A,
          const fmpz_mod_mpoly_t B, slong var, const fmpz_mod_mpoly_ctx_t ctx);


/* Evaluation ****************************************************************/

FLINT_DLL void _fmpz_mod_mpoly_eval_all_fmpz_mod(fmpz_t eval,
                        const fmpz * Acoeffs, const ulong * Aexps, slong Alen,
                            flint_bitcnt_t Abits, const fmpz * alphas,
                            const mpoly_ctx_t mctx, const fmpz_mod_ctx_t fctx);

FLINT_DLL void fmpz_mod_mpoly_evaluate_all_fmpz(fmpz_t eval,
                            const fmpz_mod_mpoly_t A, fmpz * const * alphas,
                                               const fmpz_mod_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mod_mpoly_evaluate_one_fmpz(fmpz_mod_mpoly_t A,
                        const fmpz_mod_mpoly_t B, slong var, const fmpz_t val,
                                               const fmpz_mod_mpoly_ctx_t ctx);

/* Multiplication ************************************************************/

FLINT_DLL void fmpz_mod_mpoly_mul(fmpz_mod_mpoly_t A,
                          const fmpz_mod_mpoly_t B, const fmpz_mod_mpoly_t C,
                                               const fmpz_mod_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mod_mpoly_mul_johnson(fmpz_mod_mpoly_t A,
                          const fmpz_mod_mpoly_t B, const fmpz_mod_mpoly_t C,
                                               const fmpz_mod_mpoly_ctx_t ctx);

/* Powering ******************************************************************/

FLINT_DLL int fmpz_mod_mpoly_pow_fmpz(fmpz_mod_mpoly_t A,
     const fmpz_mod_mpoly_t B, const fmpz_t k, const fmpz_mod_mpoly_ctx_t ctx);

FLINT_DLL int fmpz_mod_mpoly_pow_ui(fmpz_mod_mpoly_t A,
            const fmpz_mod_mpoly_t B, ulong k, const fmpz_mod_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mod_mpoly_pow_rmul(fmpz_mod_mpoly_t A,
            const fmpz_mod_mpoly_t B, ulong k, const fmpz_mod_mpoly_ctx_t ctx);

/* Division ******************************************************************/

FLINT_DLL int fmpz_mod_mpoly_divides(fmpz_mod_mpoly_t Q,
                        const fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_t B,
                                               const fmpz_mod_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mod_mpoly_div(fmpz_mod_mpoly_t Q,
                        const fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_t B,
                                               const fmpz_mod_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mod_mpoly_divrem(fmpz_mod_mpoly_t Q, fmpz_mod_mpoly_t R,
                        const fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_t B,
                                               const fmpz_mod_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mod_mpoly_divrem_ideal(fmpz_mod_mpoly_struct ** Q,
                            fmpz_mod_mpoly_t R, const fmpz_mod_mpoly_t A,
                            fmpz_mod_mpoly_struct * const * B, slong len,
                                               const fmpz_mod_mpoly_ctx_t ctx);

FLINT_DLL int fmpz_mod_mpoly_divides_monagan_pearce(fmpz_mod_mpoly_t Q,
                        const fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_t B,
                                               const fmpz_mod_mpoly_ctx_t ctx);
FLINT_DLL void fmpz_mod_mpoly_div_monagan_pearce(fmpz_mod_mpoly_t Q,
                          const fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_t B,
                                               const fmpz_mod_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mod_mpoly_divrem_monagan_pearce(fmpz_mod_mpoly_t Q,
        fmpz_mod_mpoly_t R, const fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_t B,
                                               const fmpz_mod_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mod_mpoly_divrem_ideal_monagan_pearce(
        fmpz_mod_mpoly_struct ** Q, fmpz_mod_mpoly_t R,
        const fmpz_mod_mpoly_t A, fmpz_mod_mpoly_struct * const * B, slong len,
                                               const fmpz_mod_mpoly_ctx_t ctx);

/* Square root ***************************************************************/

FLINT_DLL int fmpz_mod_mpoly_sqrt_heap(fmpz_mod_mpoly_t Q,
                     const fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_ctx_t ctx);

FMPZ_MOD_MPOLY_INLINE
int fmpz_mod_mpoly_sqrt(fmpz_mod_mpoly_t Q, const fmpz_mod_mpoly_t A,
                                                const fmpz_mod_mpoly_ctx_t ctx)
{
    return fmpz_mod_mpoly_sqrt_heap(Q, A, ctx);
}

FMPZ_MOD_MPOLY_INLINE
int fmpz_mod_mpoly_is_square(const fmpz_mod_mpoly_t A,
                                                const fmpz_mod_mpoly_ctx_t ctx)
{
    int res;
    fmpz_mod_mpoly_t Q;
    fmpz_mod_mpoly_init(Q, ctx);
    res = fmpz_mod_mpoly_sqrt_heap(Q, A, ctx);
    fmpz_mod_mpoly_clear(Q, ctx);
    return res;
}

FLINT_DLL int fmpz_mod_mpoly_quadratic_root(fmpz_mod_mpoly_t Q,
                    const fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_t B,
                                               const fmpz_mod_mpoly_ctx_t ctx);

/* GCD ***********************************************************************/

FLINT_DLL void fmpz_mod_mpoly_term_content(fmpz_mod_mpoly_t M,
                     const fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_ctx_t ctx);

FLINT_DLL int fmpz_mod_mpoly_content_vars(fmpz_mod_mpoly_t g,
                    const fmpz_mod_mpoly_t A, slong * vars, slong vars_length,
                                               const fmpz_mod_mpoly_ctx_t ctx);

FLINT_DLL int fmpz_mod_mpoly_gcd(fmpz_mod_mpoly_t G, const fmpz_mod_mpoly_t A,
                     const fmpz_mod_mpoly_t B, const fmpz_mod_mpoly_ctx_t ctx);

FLINT_DLL int _fmpz_mod_mpoly_gcd_algo(fmpz_mod_mpoly_t G,
                            fmpz_mod_mpoly_t Abar, fmpz_mod_mpoly_t Bbar,
                            const fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_t B,
                            const fmpz_mod_mpoly_ctx_t ctx, unsigned int algo);

FLINT_DLL int fmpz_mod_mpoly_gcd_cofactors(fmpz_mod_mpoly_t G,
                        fmpz_mod_mpoly_t Abar, fmpz_mod_mpoly_t Bbar,
                        const fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_t B,
                                               const fmpz_mod_mpoly_ctx_t ctx);

FLINT_DLL int fmpz_mod_mpoly_gcd_brown(fmpz_mod_mpoly_t G,
                        const fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_t B,
                                               const fmpz_mod_mpoly_ctx_t ctx);

FLINT_DLL int fmpz_mod_mpoly_gcd_hensel(fmpz_mod_mpoly_t G,
                        const fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_t B,
                                               const fmpz_mod_mpoly_ctx_t ctx);

FLINT_DLL int fmpz_mod_mpoly_gcd_zippel(fmpz_mod_mpoly_t G,
                        const fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_t B,
                                               const fmpz_mod_mpoly_ctx_t ctx);

FLINT_DLL int fmpz_mod_mpoly_gcd_zippel2(fmpz_mod_mpoly_t G,
                        const fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_t B,
                                               const fmpz_mod_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mod_mpoly_deflation(fmpz * shift, fmpz * stride,
                     const fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mod_mpoly_deflate(fmpz_mod_mpoly_t A,
            const fmpz_mod_mpoly_t B, const fmpz * shift, const fmpz * stride,
                                               const fmpz_mod_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mod_mpoly_inflate(fmpz_mod_mpoly_t A,
            const fmpz_mod_mpoly_t B, const fmpz * shift, const fmpz * stride,
                                               const fmpz_mod_mpoly_ctx_t ctx);

/******************************************************************************

   Internal functions (guaranteed to change without notice)

******************************************************************************/

FLINT_DLL void fmpz_mod_mpolyl_lead_coeff(fmpz_mod_mpoly_t c,
     const fmpz_mod_mpoly_t A, slong num_vars, const fmpz_mod_mpoly_ctx_t ctx);

FLINT_DLL int fmpz_mod_mpolyl_content(fmpz_mod_mpoly_t g,
     const fmpz_mod_mpoly_t A, slong num_vars, const fmpz_mod_mpoly_ctx_t ctx);

FLINT_DLL void _fmpz_mod_mpoly_to_fmpz_mod_poly_deflate(
    fmpz_mod_poly_t A,
    const fmpz_mod_mpoly_t B,
    slong var,
    const ulong * Bshift,
    const ulong * Bstride,
    const fmpz_mod_mpoly_ctx_t ctx);

FLINT_DLL void _fmpz_mod_mpoly_from_fmpz_mod_poly_inflate(
    fmpz_mod_mpoly_t A,
    flint_bitcnt_t Abits,
    const fmpz_mod_poly_t B,
    slong var,
    const ulong * Ashift,
    const ulong * Astride,
    const fmpz_mod_mpoly_ctx_t ctx);

FLINT_DLL void _fmpz_mod_mpoly_set_nmod_mpoly(
    fmpz_mod_mpoly_t A,
    const fmpz_mod_mpoly_ctx_t ctx,
    const nmod_mpoly_t nA,
    const nmod_mpoly_ctx_t nctx);

FLINT_DLL void _fmpz_mod_mpoly_get_nmod_mpoly(
    nmod_mpoly_t nA,
    const nmod_mpoly_ctx_t nctx,
    const fmpz_mod_mpoly_t A,
    const fmpz_mod_mpoly_ctx_t ctx);

FLINT_DLL int fmpz_mod_mpoly_repack_bits(fmpz_mod_mpoly_t A,
                            const fmpz_mod_mpoly_t B, flint_bitcnt_t Abits,
                                               const fmpz_mod_mpoly_ctx_t ctx);

FLINT_DLL int fmpz_mod_mpoly_repack_bits_inplace(fmpz_mod_mpoly_t A,
                         flint_bitcnt_t Abits, const fmpz_mod_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mod_mpoly_to_mpolyl_perm_deflate(
    fmpz_mod_mpoly_t A,
    const fmpz_mod_mpoly_ctx_t lctx,
    const fmpz_mod_mpoly_t B,
    const fmpz_mod_mpoly_ctx_t ctx,
    const slong * perm,
    const ulong * shift,
    const ulong * stride);

FLINT_DLL void fmpz_mod_mpoly_from_mpolyl_perm_inflate(
    fmpz_mod_mpoly_t A,
    flint_bitcnt_t Abits,
    const fmpz_mod_mpoly_ctx_t ctx,
    const fmpz_mod_mpoly_t B,
    const fmpz_mod_mpoly_ctx_t lctx,
    const slong * perm,
    const ulong * shift,
    const ulong * stride);


/* CRT ***********************************************************************/

/* instructions do A = B*c_modulus + C*b_modulus */
typedef struct
{
    slong a_idx; /* index of A */
    slong b_idx; /* index of B */
    slong c_idx; /* index of C */
    fmpz_mod_poly_t b_modulus;
    fmpz_mod_poly_t c_modulus;
} _fmpz_mod_poly_multi_crt_prog_instr;

typedef struct
{
    _fmpz_mod_poly_multi_crt_prog_instr * prog; /* straight line program */
    fmpz_mod_poly_struct * moduli, * invmoduli, * fracmoduli;
    slong length; /* length of prog */
    slong alloc;  /* alloc of prog */
    slong localsize; /* length of outputs required in poly_multi_mod_run */
    slong temp1loc, temp2loc, temp3loc, temp4loc;
    int good;   /* the moduli are good for CRT, essentially relatively prime */
} fmpz_mod_poly_multi_crt_struct;

typedef fmpz_mod_poly_multi_crt_struct fmpz_mod_poly_multi_crt_t[1];

FLINT_DLL void fmpz_mod_poly_multi_crt_init(fmpz_mod_poly_multi_crt_t CRT,
                                                     const fmpz_mod_ctx_t ctx);

FLINT_DLL int fmpz_mod_poly_multi_crt_precompute(fmpz_mod_poly_multi_crt_t CRT,
                               const fmpz_mod_poly_struct * moduli, slong len,
                                                     const fmpz_mod_ctx_t ctx);

FLINT_DLL void fmpz_mod_poly_multi_crt_precomp(fmpz_mod_poly_t output,
                const fmpz_mod_poly_multi_crt_t CRT,
                const fmpz_mod_poly_struct * values, const fmpz_mod_ctx_t ctx);

FLINT_DLL int fmpz_mod_poly_multi_crt(fmpz_mod_poly_t output,
                                const fmpz_mod_poly_struct * moduli,
                                const fmpz_mod_poly_struct * values, slong len,
                                                     const fmpz_mod_ctx_t ctx);

FLINT_DLL void fmpz_mod_poly_multi_crt_clear(fmpz_mod_poly_multi_crt_t CRT,
                                                     const fmpz_mod_ctx_t ctx);

FMPZ_MOD_MPOLY_INLINE
slong _fmpz_mod_poly_multi_crt_local_size(const fmpz_mod_poly_multi_crt_t CRT)
{
    return CRT->localsize;
}

FLINT_DLL void _fmpz_mod_poly_multi_crt_run(fmpz_mod_poly_struct * outputs,
                                        const fmpz_mod_poly_multi_crt_t CRT,
                const fmpz_mod_poly_struct * inputs, const fmpz_mod_ctx_t ctx);


/* instructions do A = B mod modulus */
typedef struct
{
    slong a_idx; /* index of A */
    slong b_idx; /* index of B */
    fmpz_mod_poly_t modulus;
} _fmpz_mod_poly_multi_mod_prog_instr;

typedef struct
{
    _fmpz_mod_poly_multi_mod_prog_instr * prog; /* straight line program */
    slong length; /* length of prog */
    slong alloc;  /* alloc of prog */
    slong localsize; /* length of temps required in poly_multi_mod_run */
    slong temp1loc;
} fmpz_mod_poly_multi_mod_struct;

typedef fmpz_mod_poly_multi_mod_struct fmpz_mod_poly_multi_mod_t[1];

FLINT_DLL void fmpz_mod_poly_multi_mod_init(
    fmpz_mod_poly_multi_mod_t P,
    const fmpz_mod_ctx_t ctx);

FLINT_DLL void fmpz_mod_poly_multi_mod_clear(
    fmpz_mod_poly_multi_mod_t P,
    const fmpz_mod_ctx_t ctx);

FLINT_DLL void fmpz_mod_poly_multi_mod_precompute(
    fmpz_mod_poly_multi_mod_t P,
    const fmpz_mod_poly_struct * f,
    slong r,
    const fmpz_mod_ctx_t ctx);

FLINT_DLL void fmpz_mod_poly_multi_mod_precomp(
    fmpz_mod_poly_t outputs,
    const fmpz_mod_poly_multi_mod_t P,
    const fmpz_mod_poly_struct * moduli,
    const fmpz_mod_ctx_t ctx);

FLINT_DLL void fmpz_mod_poly_multi_mod(
    fmpz_mod_poly_struct * outputs,
    const fmpz_mod_poly_struct * moduli,
    const fmpz_mod_poly_t input,
    slong len,
    const fmpz_mod_ctx_t ctx);

FMPZ_MOD_MPOLY_INLINE
slong _fmpz_mod_poly_multi_mod_local_size(const fmpz_mod_poly_multi_mod_t P)
{
    return P->localsize;
}

FLINT_DLL void _fmpz_mod_poly_multi_mod_run(
    fmpz_mod_poly_struct * outputs,
    const fmpz_mod_poly_multi_mod_t P,
    const fmpz_mod_poly_t input,
    fmpz_mod_poly_struct * tmps,
    const fmpz_mod_ctx_t ctx);


/* Univariates ***************************************************************/

FLINT_DLL void fmpz_mod_mpoly_univar_init(fmpz_mod_mpoly_univar_t A,
                                               const fmpz_mod_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mod_mpoly_univar_clear(fmpz_mod_mpoly_univar_t A,
                                               const fmpz_mod_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mod_mpoly_univar_fit_length(fmpz_mod_mpoly_univar_t A,
                                 slong length, const fmpz_mod_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mod_mpoly_univar_print_pretty(const fmpz_mod_mpoly_univar_t A,
                              const char ** x, const fmpz_mod_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mod_mpoly_univar_assert_canonical(fmpz_mod_mpoly_univar_t A,
                                               const fmpz_mod_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mod_mpoly_to_univar(fmpz_mod_mpoly_univar_t A,
          const fmpz_mod_mpoly_t B, slong var, const fmpz_mod_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mod_mpoly_from_univar_bits(fmpz_mod_mpoly_t A,
            flint_bitcnt_t Abits, const fmpz_mod_mpoly_univar_t B, slong var,
                                               const fmpz_mod_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mod_mpoly_from_univar(fmpz_mod_mpoly_t A,
                                  const fmpz_mod_mpoly_univar_t B, slong var,
                                               const fmpz_mod_mpoly_ctx_t ctx);

FMPZ_MOD_MPOLY_INLINE
void fmpz_mod_mpoly_univar_swap(fmpz_mod_mpoly_univar_t A,
                     fmpz_mod_mpoly_univar_t B, const fmpz_mod_mpoly_ctx_t ctx)
{
    fmpz_mod_mpoly_univar_struct t = *A;
    *A = *B;
    *B = t;
}

FMPZ_MOD_MPOLY_INLINE
int fmpz_mod_mpoly_univar_degree_fits_si(const fmpz_mod_mpoly_univar_t A,
                                                const fmpz_mod_mpoly_ctx_t ctx)
{
    return A->length == 0 || fmpz_fits_si(A->exps + 0);
}

FMPZ_MOD_MPOLY_INLINE
slong fmpz_mod_mpoly_univar_length(const fmpz_mod_mpoly_univar_t A,
                                                const fmpz_mod_mpoly_ctx_t ctx)
{
    return A->length;
}

FMPZ_MOD_MPOLY_INLINE
slong fmpz_mod_mpoly_univar_get_term_exp_si(fmpz_mod_mpoly_univar_t A, slong i,
                                                const fmpz_mod_mpoly_ctx_t ctx)
{
    FLINT_ASSERT(i < (ulong)A->length);
    return fmpz_get_si(A->exps + i);
}

FMPZ_MOD_MPOLY_INLINE
void fmpz_mod_mpoly_univar_get_term_coeff(fmpz_mod_mpoly_t c,
      const fmpz_mod_mpoly_univar_t A, slong i, const fmpz_mod_mpoly_ctx_t ctx)
{
    FLINT_ASSERT((ulong)i < (ulong)A->length);
    fmpz_mod_mpoly_set(c, A->coeffs + i, ctx);
}

FMPZ_MOD_MPOLY_INLINE
void fmpz_mod_mpoly_univar_swap_term_coeff(fmpz_mod_mpoly_t c,
            fmpz_mod_mpoly_univar_t A, slong i, const fmpz_mod_mpoly_ctx_t ctx)
{
    FLINT_ASSERT((ulong)i < (ulong)A->length);
    fmpz_mod_mpoly_swap(c, A->coeffs + i, ctx);
}


/* mpolyn and mpolyun ********************************************************/

FLINT_DLL void fmpz_mod_mpolyn_init(fmpz_mod_mpolyn_t A, flint_bitcnt_t bits,
                                               const fmpz_mod_mpoly_ctx_t ctx);

FMPZ_MOD_MPOLY_INLINE void fmpz_mod_mpolyn_swap(fmpz_mod_mpolyn_t A,
                           fmpz_mod_mpolyn_t B, const fmpz_mod_mpoly_ctx_t ctx)
{
   fmpz_mod_mpolyn_struct t = *A;
   *A = *B;
   *B = t;    
}

FMPZ_MOD_MPOLY_INLINE void fmpz_mod_mpolyn_zero(fmpz_mod_mpolyn_t A,
                                                const fmpz_mod_mpoly_ctx_t ctx)
{
    A->length = 0;
}

FMPZ_MOD_MPOLY_INLINE int fmpz_mod_mpolyn_is_zero(fmpz_mod_mpolyn_t A,
                                                const fmpz_mod_mpoly_ctx_t ctx)
{
    return A->length == 0;
}

FMPZ_MOD_MPOLY_INLINE
const fmpz_mod_poly_struct * fmpz_mod_mpolyn_leadcoeff_poly(
                                                     const fmpz_mod_mpolyn_t A)
{
    FLINT_ASSERT(A->length > 0);
    return A->coeffs + 0;
}


FLINT_DLL void fmpz_mod_mpolyn_fit_length(fmpz_mod_mpolyn_t A,
                                 slong length, const fmpz_mod_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mod_mpolyun_fit_length(fmpz_mod_mpolyun_t A,
                                 slong length, const fmpz_mod_mpoly_ctx_t ctx);

FMPZ_MOD_MPOLY_INLINE
const fmpz * fmpz_mod_mpolyn_leadcoeff(const fmpz_mod_mpolyn_t A)
{
    FLINT_ASSERT(A->length > 0);
    FLINT_ASSERT(A->coeffs[0].length > 0);
    return A->coeffs[0].coeffs + A->coeffs[0].length - 1;
}

FLINT_DLL void fmpz_mod_mpolyn_fit_bits(
    fmpz_mod_mpolyn_t A,
    slong bits,
    const fmpz_mod_mpoly_ctx_t ctx);

FLINT_DLL fmpz_mod_poly_struct * fmpz_mod_mpolyun_leadcoeff_ref(
                         fmpz_mod_mpolyun_t A, const fmpz_mod_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mod_mpolyn_content_poly(fmpz_mod_poly_t a,
                    const fmpz_mod_mpolyn_t B, const fmpz_mod_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mod_mpolyun_content_last(fmpz_mod_poly_t a,
                   const fmpz_mod_mpolyun_t B, const fmpz_mod_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mod_mpolyn_divexact_poly(fmpz_mod_mpolyn_t A,
                      const fmpz_mod_poly_t b, const fmpz_mod_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mod_mpolyun_divexact_last(fmpz_mod_mpolyun_t A,
                      const fmpz_mod_poly_t b, const fmpz_mod_mpoly_ctx_t ctx);

FLINT_DLL ulong nmod_mpolyn_bidegree(const nmod_mpolyn_t A);

FLINT_DLL ulong fmpz_mod_mpolyn_bidegree(const fmpz_mod_mpolyn_t A);

FLINT_DLL slong fmpz_mod_mpolyn_lastdeg(const fmpz_mod_mpolyn_t A,
                                               const fmpz_mod_mpoly_ctx_t ctx);

FLINT_DLL slong fmpz_mod_mpolyun_lastdeg(const fmpz_mod_mpolyun_t A,
                                               const fmpz_mod_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mod_mpolyun_init(fmpz_mod_mpolyun_t A, flint_bitcnt_t bits,
                                               const fmpz_mod_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mod_mpolyn_clear(fmpz_mod_mpolyn_t A,
                                               const fmpz_mod_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mod_mpolyun_clear(fmpz_mod_mpolyun_t A,
                                               const fmpz_mod_mpoly_ctx_t ctx);

FLINT_DLL fmpz * fmpz_mod_mpolyun_leadcoeff_last_ref(const fmpz_mod_mpolyun_t A,
                                               const fmpz_mod_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mod_mpolyn_one(fmpz_mod_mpolyn_t A,
                                               const fmpz_mod_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mod_mpolyun_one(fmpz_mod_mpolyun_t A,
                                               const fmpz_mod_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mod_mpolyn_mul_poly(fmpz_mod_mpolyn_t A,
                            fmpz_mod_poly_t b, const fmpz_mod_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mod_mpolyun_mul_last(fmpz_mod_mpolyun_t A, fmpz_mod_poly_t b,
                                               const fmpz_mod_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mod_mpolyn_scalar_mul_fmpz_mod(fmpz_mod_mpolyn_t A,
                               const fmpz_t c, const fmpz_mod_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mod_mpolyun_scalar_mul_fmpz_mod(fmpz_mod_mpolyun_t A,
                               const fmpz_t c, const fmpz_mod_mpoly_ctx_t ctx);

FLINT_DLL int fmpz_mod_mpolyn_equal(const fmpz_mod_mpolyn_t A,
                   const fmpz_mod_mpolyn_t B, const fmpz_mod_mpoly_ctx_t ctx);

FLINT_DLL int fmpz_mod_mpolyun_equal(const fmpz_mod_mpolyun_t A,
                   const fmpz_mod_mpolyun_t B, const fmpz_mod_mpoly_ctx_t ctx);

FLINT_DLL int fmpz_mod_mpolyn_gcd_brown_bivar(
     fmpz_mod_mpolyn_t G, fmpz_mod_mpolyn_t Abar, fmpz_mod_mpolyn_t Bbar,
     fmpz_mod_mpolyn_t A, fmpz_mod_mpolyn_t B, const fmpz_mod_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mod_mpolyn_print_pretty(const fmpz_mod_mpolyn_t poly,
                              const char ** x, const fmpz_mod_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mod_mpolyun_print_pretty(const fmpz_mod_mpolyun_t poly,
                              const char ** x, const fmpz_mod_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mod_mpoly_cvtfrom_mpolyn(
    fmpz_mod_mpoly_t A,
    const fmpz_mod_mpolyn_t B,
    slong var,
    const fmpz_mod_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mod_mpoly_cvtto_mpolyn(
    fmpz_mod_mpolyn_t A,
    const fmpz_mod_mpoly_t B,
    slong var,
    const fmpz_mod_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mod_mpoly_to_mpolyn_perm_deflate(
    fmpz_mod_mpolyn_t A,
    const fmpz_mod_mpoly_ctx_t nctx,
    const fmpz_mod_mpoly_t B,
    const fmpz_mod_mpoly_ctx_t ctx,
    const slong * perm,
    const ulong * shift,
    const ulong * stride);

FLINT_DLL void fmpz_mod_mpoly_from_mpolyn_perm_inflate(
    fmpz_mod_mpoly_t A,
    flint_bitcnt_t Abits,
    const fmpz_mod_mpoly_ctx_t ctx,
    const fmpz_mod_mpolyn_t B,
    const fmpz_mod_mpoly_ctx_t nctx,
    const slong * perm,
    const ulong * shift,
    const ulong * stride);

/* interp ********************************************************************/

FLINT_DLL void fmpz_mod_mpolyn_intp_reduce_sm_poly(fmpz_mod_poly_t E,
                            const fmpz_mod_mpolyn_t A, const fmpz_t alpha,
                                               const fmpz_mod_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mod_mpolyn_intp_lift_sm_poly(fmpz_mod_mpolyn_t A,
                      const fmpz_mod_poly_t B, const fmpz_mod_mpoly_ctx_t ctx);

FLINT_DLL int fmpz_mod_mpolyn_intp_crt_sm_poly(slong * lastdeg_,
             fmpz_mod_mpolyn_t F, fmpz_mod_mpolyn_t T, const fmpz_mod_poly_t A,
                            const fmpz_mod_poly_t modulus, const fmpz_t alpha,
                                               const fmpz_mod_mpoly_ctx_t ctx);

/* helpers *******************************************************************/

typedef struct {
    nmod_berlekamp_massey_struct * coeffs;
    ulong * exps;
    slong length;
    slong alloc;
    slong pointcount;
} nmod_bma_mpoly_struct;

typedef nmod_bma_mpoly_struct nmod_bma_mpoly_t[1];

typedef struct {
    fmpz_mod_berlekamp_massey_struct * coeffs;
    ulong * exps;
    slong length;
    slong alloc;
    slong pointcount;
} fmpz_mod_bma_mpoly_struct;

typedef fmpz_mod_bma_mpoly_struct fmpz_mod_bma_mpoly_t[1];

typedef struct
{
    slong * degbounds;
    ulong * subdegs;
    fmpz_mod_discrete_log_pohlig_hellman_t dlogenv;
    nmod_discrete_log_pohlig_hellman_t dlogenv_sp;
} mpoly_bma_interpolate_ctx_struct;
typedef mpoly_bma_interpolate_ctx_struct mpoly_bma_interpolate_ctx_t[1];


FLINT_DLL void nmod_bma_mpoly_init(nmod_bma_mpoly_t A);

FLINT_DLL void nmod_bma_mpoly_reset_prime(nmod_bma_mpoly_t A, nmod_t fpctx);

FLINT_DLL void nmod_bma_mpoly_clear(nmod_bma_mpoly_t A);

FLINT_DLL void nmod_bma_mpoly_print(const nmod_bma_mpoly_t A);

FLINT_DLL void nmod_bma_mpoly_fit_length(nmod_bma_mpoly_t A, slong length,
                                                                 nmod_t fpctx);

FLINT_DLL void nmod_bma_mpoly_zero(nmod_bma_mpoly_t L);

FLINT_DLL int nmod_bma_mpoly_reduce(nmod_bma_mpoly_t L);

FLINT_DLL void nmod_bma_mpoly_add_point(
    nmod_bma_mpoly_t L,
    const n_bpoly_t A,
    const nmod_mpoly_ctx_t ctx_sp);

FLINT_DLL int nmod_bma_mpoly_get_fmpz_mpolyu(fmpz_mpolyu_t A,
      const fmpz_mpoly_ctx_t ctx, ulong alphashift, const nmod_bma_mpoly_t L,
                         const mpoly_bma_interpolate_ctx_t Ictx, nmod_t fpctx);



FLINT_DLL void fmpz_mod_bma_mpoly_init(fmpz_mod_bma_mpoly_t A);

FLINT_DLL void fmpz_mod_bma_mpoly_clear(fmpz_mod_bma_mpoly_t A,
                                                   const fmpz_mod_ctx_t fpctx);

FLINT_DLL void fmpz_mod_bma_mpoly_fit_length(fmpz_mod_bma_mpoly_t A,
                                     slong length, const fmpz_mod_ctx_t fpctx);

FLINT_DLL void fmpz_mod_bma_mpoly_zero(fmpz_mod_bma_mpoly_t L);

FLINT_DLL void fmpz_mod_bma_mpoly_add_point(fmpz_mod_bma_mpoly_t L,
                  const fmpz_mod_mpolyn_t A, const fmpz_mod_mpoly_ctx_t ctx_mp);

FLINT_DLL int fmpz_mod_bma_mpoly_get_fmpz_mpolyu(fmpz_mpolyu_t A,
    const fmpz_mpoly_ctx_t ctx, const fmpz_t alphashift, const fmpz_mod_bma_mpoly_t L,
           const mpoly_bma_interpolate_ctx_t Ictx, const fmpz_mod_ctx_t fpctx);

FMPZ_MOD_MPOLY_INLINE
void mpoly_bma_interpolate_ctx_init(mpoly_bma_interpolate_ctx_t I, slong nvars)
{
    I->degbounds = (slong *) flint_malloc(nvars*sizeof(slong));
    I->subdegs   = (ulong *) flint_malloc(nvars*sizeof(ulong));
    fmpz_mod_discrete_log_pohlig_hellman_init(I->dlogenv);
    nmod_discrete_log_pohlig_hellman_init(I->dlogenv_sp);
}

FMPZ_MOD_MPOLY_INLINE
void mpoly_bma_interpolate_ctx_clear(mpoly_bma_interpolate_ctx_t I)
{
    flint_free(I->degbounds);
    flint_free(I->subdegs);
    fmpz_mod_discrete_log_pohlig_hellman_clear(I->dlogenv);
    nmod_discrete_log_pohlig_hellman_clear(I->dlogenv_sp);
}

FLINT_DLL int nmod_mpoly_bma_get_fmpz_mpoly(fmpz_mpoly_t A,
     const fmpz_mpoly_ctx_t ctx, ulong alphashift, nmod_berlekamp_massey_t I,
                         const mpoly_bma_interpolate_ctx_t Ictx, nmod_t fpctx);

FLINT_DLL int fmpz_mod_bma_get_fmpz_mpoly(fmpz_mpoly_t A,
     const fmpz_mpoly_ctx_t ctx, const fmpz_t alphashift, fmpz_mod_berlekamp_massey_t I,
           const mpoly_bma_interpolate_ctx_t Ictx, const fmpz_mod_ctx_t fpctx);


FLINT_DLL void nmod_mpoly_bma_interpolate_alpha_powers(mp_limb_t * out,
                            ulong w, const mpoly_bma_interpolate_ctx_t Ictx,
                                     const fmpz_mpoly_ctx_t ctx, nmod_t fpctx);


FLINT_DLL void fmpz_mod_mpoly_bma_interpolate_alpha_powers(fmpz * out,
                     const fmpz_t w, const mpoly_bma_interpolate_ctx_t Ictx,
                       const fmpz_mpoly_ctx_t ctx, const fmpz_mod_ctx_t fpctx);

/* skel */

FLINT_DLL void fmpz_mod_mpoly_red_skel(fmpz_mpolyc_t Ared, const fmpz_mpoly_t A,
                                                      const fmpz_mod_ctx_t fpctx);

FLINT_DLL void fmpz_mod_mpolyu_red_skel(fmpz_mpolycu_t Ared, const fmpz_mpolyu_t A,
                                                      const fmpz_mod_ctx_t fpctx);

FLINT_DLL void fmpz_mod_mpoly_copy_skel(fmpz_mpolyc_t M, const fmpz_mpolyc_t S);

FLINT_DLL void fmpz_mod_mpolyu_copy_skel(fmpz_mpolycu_t M, const fmpz_mpolycu_t S);

FLINT_DLL void fmpz_mod_mpoly_pow_skel(fmpz_mpolyc_t M, const fmpz_mpolyc_t S,
                                          ulong k, const fmpz_mod_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mod_mpolyu_pow_skel(fmpz_mpolycu_t M, const fmpz_mpolycu_t S,
                                          ulong k, const fmpz_mod_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mod_mpoly_set_skel(fmpz_mpolyc_t S,
                      const fmpz_mod_mpoly_ctx_t ctx_mp, const fmpz_mpoly_t A,
                               const fmpz * alpha, const fmpz_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mod_mpolyu_set_skel(fmpz_mpolycu_t S,
                     const fmpz_mod_mpoly_ctx_t ctx_mp, const fmpz_mpolyu_t A,
                               const fmpz * alpha, const fmpz_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mod_mpoly_use_skel_mul(fmpz_t eval, fmpz_mpolyc_t Ared,
                               fmpz_mpolyc_t Avar, const fmpz_mpolyc_t Ainc,
                                                   const fmpz_mod_ctx_t fpctx);

FLINT_DLL void fmpz_mod_mpolyuu_use_skel_mul(fmpz_mod_mpolyn_t E,
        const fmpz_mpolyu_t A, const fmpz_mpolycu_t Ared, fmpz_mpolycu_t Acur,
                 const fmpz_mpolycu_t Ainc, const fmpz_mod_mpoly_ctx_t ctx_mp);


/* eval */

FLINT_DLL mp_limb_t fmpz_mpoly_eval_nmod(nmod_t fpctx,
    const fmpz_mpoly_t A, const mp_limb_t * alpha, const fmpz_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mpoly_eval_fmpz_mod(fmpz_t eval,
                        const fmpz_mod_ctx_t fpctx, const fmpz_mpoly_t A,
                               const fmpz * alpha, const fmpz_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mpolyuu_eval_nmod(
    n_bpoly_t E,
    const nmod_mpoly_ctx_t ctx_sp,
    const fmpz_mpolyu_t A,
    const mp_limb_t * alpha,
    const fmpz_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mpolyuu_eval_fmpz_mod(fmpz_mod_mpolyn_t E,
                 const fmpz_mod_mpoly_ctx_t ctx_mp, const fmpz_mpolyu_t A,
                               const fmpz * alpha, const fmpz_mpoly_ctx_t ctx);

FLINT_DLL void _mpoly_monomial_evals_nmod(
    mp_limb_t * E,
    const ulong * Aexps,
    flint_bitcnt_t Abits,
    slong Alen,
    const mp_limb_t * alpha,
    slong vstart,
    const mpoly_ctx_t mctx,
    nmod_t fctx);

FLINT_DLL void nmod_mpolyuu_eval_step2(
    n_bpoly_t E,
    n_bpoly_t Acur,
    const n_polyun_t Ainc,
    const nmod_mpoly_ctx_t ctx_sp);

/* gcd ***********************************************************************/

FLINT_DLL int fmpz_mod_mpolyn_gcd_brown_smprime(
    fmpz_mod_mpolyn_t G,
    fmpz_mod_mpolyn_t Abar,
    fmpz_mod_mpolyn_t Bbar,
    fmpz_mod_mpolyn_t A,
    fmpz_mod_mpolyn_t B,
    slong var,
    const fmpz_mod_mpoly_ctx_t ctx,
    const mpoly_gcd_info_t I,
    fmpz_mod_poly_polyun_mpolyn_stack_t St);

FLINT_DLL int fmpz_mod_mpolyl_gcd_zippel_smprime(
    fmpz_mod_mpoly_t rG, const slong * rGdegs, /* guess at rG degrees, could be NULL */
    fmpz_mod_mpoly_t rAbar,
    fmpz_mod_mpoly_t rBbar,
    const fmpz_mod_mpoly_t A, const slong * Adegs,
    const fmpz_mod_mpoly_t B, const slong * Bdegs,
    const fmpz_mod_mpoly_t gamma, const slong * gammadegs,
    const fmpz_mod_mpoly_ctx_t ctx);

FLINT_DLL int fmpz_mod_mpolyl_gcd_hensel_smprime(
    fmpz_mod_mpoly_t G, slong Gdeg,
    fmpz_mod_mpoly_t Abar,
    fmpz_mod_mpoly_t Bbar,
    const fmpz_mod_mpoly_t A,
    const fmpz_mod_mpoly_t B,
    const fmpz_mod_mpoly_ctx_t ctx);

/* geobuckets ****************************************************************/

typedef struct fmpz_mod_mpoly_geobucket
{
    fmpz_mod_mpoly_struct polys[FLINT_BITS/2];
    fmpz_mod_mpoly_struct temps[FLINT_BITS/2];
    slong length;
} fmpz_mod_mpoly_geobucket_struct;

typedef fmpz_mod_mpoly_geobucket_struct fmpz_mod_mpoly_geobucket_t[1];

FLINT_DLL void fmpz_mod_mpoly_geobucket_init(fmpz_mod_mpoly_geobucket_t B,
                                               const fmpz_mod_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mod_mpoly_geobucket_clear(fmpz_mod_mpoly_geobucket_t B,
                                               const fmpz_mod_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mod_mpoly_geobucket_empty(fmpz_mod_mpoly_t p,
                 fmpz_mod_mpoly_geobucket_t B, const fmpz_mod_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mod_mpoly_geobucket_fit_length(fmpz_mod_mpoly_geobucket_t B,
                                      slong i, const fmpz_mod_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mod_mpoly_geobucket_set(fmpz_mod_mpoly_geobucket_t B,
                           fmpz_mod_mpoly_t p, const fmpz_mod_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mod_mpoly_geobucket_add(fmpz_mod_mpoly_geobucket_t B,
                           fmpz_mod_mpoly_t p, const fmpz_mod_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mod_mpoly_geobucket_sub(fmpz_mod_mpoly_geobucket_t B,
                           fmpz_mod_mpoly_t p, const fmpz_mod_mpoly_ctx_t ctx);

/******************************************************************************

   Internal consistency checks

******************************************************************************/

/*
   test that r is a valid remainder upon division by g
   this means that no monomial of r is divisible by lm(g)
*/
FMPZ_MOD_MPOLY_INLINE
void fmpz_mod_mpoly_remainder_strongtest(const fmpz_mod_mpoly_t r,
                      const fmpz_mod_mpoly_t g, const fmpz_mod_mpoly_ctx_t ctx)
{
   slong i, N, bits;
   ulong mask = 0;
   ulong * rexp, * gexp;

   bits = FLINT_MAX(r->bits, g->bits);
   N = mpoly_words_per_exp(bits, ctx->minfo);

    if (g->length == 0)
        flint_throw(FLINT_ERROR, "Zero denominator in remainder test");

    if (r->length == 0)
        return;

    rexp = (ulong *) flint_malloc(N*r->length*sizeof(ulong));
    gexp = (ulong *) flint_malloc(N*1        *sizeof(ulong));
    mpoly_repack_monomials(rexp, bits, r->exps, r->bits, r->length, ctx->minfo);
    mpoly_repack_monomials(gexp, bits, g->exps, g->bits, 1,         ctx->minfo);

    mask = (bits <= FLINT_BITS) ? mpoly_overflow_mask_sp(bits) : 0;

    for (i = 0; i < r->length; i++)
    {
        int divides;

        if (bits <= FLINT_BITS)
            divides = mpoly_monomial_divides_test(rexp + i*N, gexp + 0*N, N, mask);
        else
            divides = mpoly_monomial_divides_mp_test(rexp + i*N, gexp + 0*N, N, bits);

        if (divides)
        {
            flint_printf("fmpz_mod_mpoly_remainder_strongtest FAILED i = %wd\n", i);
            flint_printf("rem ");fmpz_mod_mpoly_print_pretty(r, NULL, ctx); printf("\n\n");
            flint_printf("den ");fmpz_mod_mpoly_print_pretty(g, NULL, ctx); printf("\n\n");
            flint_abort();
        }
    }

   flint_free(rexp);
   flint_free(gexp);
}

#ifdef __cplusplus
}
#endif

#endif
