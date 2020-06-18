/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#ifndef NMOD_MPOLY_FACTOR_H
#define NMOD_MPOLY_FACTOR_H

#ifdef NMOD_MPOLY_FACTOR_INLINES_C
#define NMOD_MPOLY_FACTOR_INLINE FLINT_DLL
#else
#define NMOD_MPOLY_FACTOR_INLINE static __inline__
#endif

#undef ulong
#define ulong ulongxx /* interferes with system includes */
#include <stdio.h>
#undef ulong
#include <gmp.h>
#define ulong mp_limb_t

#include "nmod_mpoly.h"

#ifdef __cplusplus
 extern "C" {
#endif

typedef struct {
    mp_limb_t content;
    nmod_mpoly_struct * poly;
    fmpz * exp;
    slong length;
    slong alloc;
} nmod_mpoly_factor_struct;

typedef nmod_mpoly_factor_struct nmod_mpoly_factor_t[1];

FLINT_DLL void nmod_mpoly_factor_init(nmod_mpoly_factor_t f, const nmod_mpoly_ctx_t ctx);

FLINT_DLL void nmod_mpoly_factor_init2(nmod_mpoly_factor_t f, slong alloc, const nmod_mpoly_ctx_t ctx);

FLINT_DLL void nmod_mpoly_factor_realloc(nmod_mpoly_factor_t f, slong alloc, const nmod_mpoly_ctx_t ctx);

FLINT_DLL void nmod_mpoly_factor_fit_length(nmod_mpoly_factor_t f, slong len, const nmod_mpoly_ctx_t ctx);

FLINT_DLL void nmod_mpoly_factor_clear(nmod_mpoly_factor_t f, const nmod_mpoly_ctx_t ctx);

FLINT_DLL void nmod_mpoly_factor_append_ui(nmod_mpoly_factor_t f, const nmod_mpoly_t A, ulong e, const nmod_mpoly_ctx_t ctx);

FLINT_DLL void nmod_mpoly_factor_append_fmpz(nmod_mpoly_factor_t f, const nmod_mpoly_t A, const fmpz_t e, const nmod_mpoly_ctx_t ctx);

FLINT_DLL void nmod_mpoly_factor_set(nmod_mpoly_factor_t a, const nmod_mpoly_factor_t b, const nmod_mpoly_ctx_t ctx);

FLINT_DLL void nmod_mpoly_factor_print_pretty(const nmod_mpoly_factor_t f, const char ** vars, const nmod_mpoly_ctx_t ctx);

FLINT_DLL int nmod_mpoly_factor(nmod_mpoly_factor_t f, const nmod_mpoly_t A, int full, const nmod_mpoly_ctx_t ctx);

FLINT_DLL void nmod_mpoly_factor_sort(nmod_mpoly_factor_t f, const nmod_mpoly_ctx_t ctx);

NMOD_MPOLY_FACTOR_INLINE
void nmod_mpoly_factor_swap(nmod_mpoly_factor_t A, nmod_mpoly_factor_t B,
                                                    const nmod_mpoly_ctx_t ctx)
{
   nmod_mpoly_factor_struct t = *A;
   *A = *B;
   *B = t;
}

FLINT_DLL void subset_first(fmpz_t a, slong n, slong r);

FLINT_DLL int subset_next(fmpz_t a, const fmpz_t b, slong n);

FLINT_DLL void subset_print(const fmpz_t a, slong n);





/**************** dense univarates (Z/nZ)[x] *********************************/

typedef struct
{
    mp_limb_t * coeffs;
    slong length;
    slong alloc;
} nmod_polydr_struct;

typedef nmod_polydr_struct nmod_polydr_t[1];

NMOD_MPOLY_FACTOR_INLINE
void nmod_polydr_init(nmod_polydr_t A)
{
    A->length = 0;
    A->alloc = 0;
    A->coeffs = NULL;
}

NMOD_MPOLY_FACTOR_INLINE
void nmod_polydr_init2(nmod_polydr_t A, slong alloc)
{
    FLINT_ASSERT(alloc >= 0);
    A->length = 0;
    A->alloc = alloc;
    A->coeffs = alloc > 0 ? flint_malloc(alloc*sizeof(mp_limb_t)) : NULL;
}

FLINT_DLL void nmod_polydr_clear(nmod_polydr_t A);

FLINT_DLL void nmod_polydr_realloc(nmod_polydr_t A, slong len);

FLINT_DLL void nmod_polydr_print_pretty(const nmod_polydr_t A, const char * x);

NMOD_MPOLY_FACTOR_INLINE
void nmod_polydr_fit_length(nmod_polydr_t A, slong len)
{
    if (len > A->alloc)
        nmod_polydr_realloc(A, len);
}


NMOD_MPOLY_FACTOR_INLINE
void nmod_poly_mock(nmod_poly_t a, const nmod_polydr_t b, nmod_t ctx)
{
    a->coeffs = b->coeffs;
    a->length = b->length;
    a->alloc = b->alloc;
    a->mod = ctx;
}

NMOD_MPOLY_FACTOR_INLINE
void nmod_polydr_mock(nmod_polydr_t a, const nmod_poly_t b)
{
    a->coeffs = b->coeffs;
    a->length = b->length;
    a->alloc = b->alloc;
}


NMOD_MPOLY_FACTOR_INLINE
void nmod_polydr_set(nmod_polydr_t A, const nmod_polydr_t B)
{
    nmod_polydr_fit_length(A, B->length);
    flint_mpn_copyi(A->coeffs, B->coeffs, B->length);
    A->length = B->length;
}


NMOD_MPOLY_FACTOR_INLINE
void nmod_polydr_swap(nmod_polydr_t A, nmod_polydr_t B)
{
    nmod_polydr_struct t = *B;
    *B = *A;
    *A = t;
}

NMOD_MPOLY_FACTOR_INLINE
void _nmod_polydr_normalise(nmod_polydr_t A)
{
    while (A->length > 0 && A->coeffs[A->length - 1] == 0)
        A->length--;
}

NMOD_MPOLY_FACTOR_INLINE
slong nmod_polydr_degree(const nmod_polydr_t A)
{
    FLINT_ASSERT(A->length >= 0);
    return A->length - 1;
}

NMOD_MPOLY_FACTOR_INLINE
int nmod_polydr_is_one(const nmod_polydr_t A)
{
    return A->length == 1 && A->coeffs[0] == 1;
}

NMOD_MPOLY_FACTOR_INLINE
void nmod_polydr_one(nmod_polydr_t A, const nmodf_ctx_t ctx)
{
    nmod_polydr_fit_length(A, 1);
    A->length = (ctx->mod.n != 1);
    A->coeffs[0] = 1;
}

NMOD_MPOLY_FACTOR_INLINE
void nmod_polydr_set_nmod(nmod_polydr_t A, mp_limb_t c, const nmodf_ctx_t ctx)
{
    FLINT_ASSERT(c < ctx->mod.n);
    nmod_polydr_fit_length(A, 1);
    A->coeffs[0] = c;
    A->length = (c != 0);
}

NMOD_MPOLY_FACTOR_INLINE
int nmod_polydr_is_zero(const nmod_polydr_t poly)
{
    return poly->length == 0;
}

NMOD_MPOLY_FACTOR_INLINE
void nmod_polydr_zero(nmod_polydr_t res)
{
    res->length = 0;
}

NMOD_MPOLY_FACTOR_INLINE
int nmod_polydr_equal(const nmod_polydr_t a, const nmod_polydr_t b,
                                                         const nmodf_ctx_t ctx)
{
    if (a->length != b->length)
        return 0;

    if (a != b)
    {
        if (!_nmod_vec_equal(a->coeffs, b->coeffs, a->length))
            return 0;
    }

    return 1;
}

NMOD_MPOLY_FACTOR_INLINE
void nmod_polydr_make_monic(nmod_polydr_t A, const nmod_polydr_t B,
                                                         const nmodf_ctx_t ctx)
{
    FLINT_ASSERT(B->length > 0);
    nmod_polydr_fit_length(A, B->length);
    _nmod_poly_make_monic(A->coeffs, B->coeffs, B->length, ctx->mod);
    A->length = B->length;
}

NMOD_MPOLY_FACTOR_INLINE
void nmod_polydr_taylor_shift(nmod_polydr_t g, const nmod_polydr_t f,
                                            mp_limb_t c, const nmodf_ctx_t ctx)
{
    if (f != g)
        nmod_polydr_set(g, f);

    _nmod_poly_taylor_shift(g->coeffs, c, g->length, ctx->mod);
}

NMOD_MPOLY_FACTOR_INLINE
ulong nmod_polydr_get_coeff_ui(const nmod_polydr_t poly, slong j)
{
    return (j >= poly->length) ? 0 : poly->coeffs[j];
}

FLINT_DLL void nmod_polydr_set_coeff_ui(nmod_polydr_t A, slong j, ulong c,
                                                        const nmodf_ctx_t ctx);


NMOD_MPOLY_FACTOR_INLINE
void nmod_polydr_set_nmod_poly(nmod_polydr_t a, const nmod_poly_t b)
{
    nmod_polydr_fit_length(a, b->length);
    flint_mpn_copyi(a->coeffs, b->coeffs, b->length);
    a->length = b->length;
}

NMOD_MPOLY_FACTOR_INLINE
void nmod_poly_set_nmod_polydr(nmod_poly_t a, const nmod_polydr_t b)
{
    nmod_poly_fit_length(a, b->length);
    flint_mpn_copyi(a->coeffs, b->coeffs, b->length);
    a->length = b->length;
}

NMOD_MPOLY_FACTOR_INLINE
void nmod_polydr_shift_left(nmod_polydr_t A, const nmod_polydr_t B,
                                                slong k, const nmodf_ctx_t ctx)
{
    nmod_polydr_fit_length(A, B->length + k);
    _nmod_poly_shift_left(A->coeffs, B->coeffs, B->length, k);
    A->length = B->length + k;
}

FLINT_DLL void nmod_polydr_scalar_mul_nmod(nmod_polydr_t A,
                    const nmod_polydr_t B, mp_limb_t c, const nmodf_ctx_t ctx);

NMOD_MPOLY_FACTOR_INLINE
mp_limb_t nmod_polydr_evaluate_nmod(const nmod_polydr_t A, mp_limb_t c,
                                                         const nmodf_ctx_t ctx)
{
    return _nmod_poly_evaluate_nmod(A->coeffs, A->length, c, ctx->mod);
}

NMOD_MPOLY_FACTOR_INLINE
void nmod_polydr_neg(nmod_polydr_t A, const nmod_polydr_t B, const nmodf_ctx_t ctx)
{
    nmod_polydr_fit_length(A, B->length);
    _nmod_vec_neg(A->coeffs, B->coeffs, B->length, ctx->mod);
    A->length = B->length;
}


FLINT_DLL void nmod_polydr_add(nmod_polydr_t A, const nmod_polydr_t B,
                                 const nmod_polydr_t C, const nmodf_ctx_t ctx);

FLINT_DLL void nmod_polydr_sub(nmod_polydr_t A, const nmod_polydr_t B,
                                 const nmod_polydr_t C, const nmodf_ctx_t ctx);

FLINT_DLL void nmod_polydr_mul(nmod_polydr_t A, const nmod_polydr_t B,
                                 const nmod_polydr_t C, const nmodf_ctx_t ctx);

FLINT_DLL void nmod_polydr_mullow(nmod_polydr_t res, const nmod_polydr_t poly1,
                const nmod_polydr_t poly2, slong trunc, const nmodf_ctx_t ctx);

FLINT_DLL void nmod_polydr_div(nmod_polydr_t Q, const nmod_polydr_t A,
                                 const nmod_polydr_t B, const nmodf_ctx_t ctx);

FLINT_DLL void nmod_polydr_rem(nmod_polydr_t R, const nmod_polydr_t A,
                                 const nmod_polydr_t B, const nmodf_ctx_t ctx);

FLINT_DLL void nmod_polydr_divrem(nmod_polydr_t Q, nmod_polydr_t R,
          const nmod_polydr_t A, const nmod_polydr_t B, const nmodf_ctx_t ctx);

FLINT_DLL int nmod_polydr_invmod(nmod_polydr_t A, const nmod_polydr_t B,
                                 const nmod_polydr_t P, const nmodf_ctx_t ctx);

FLINT_DLL void nmod_polydr_gcd(nmod_polydr_t G, const nmod_polydr_t A,
                                 const nmod_polydr_t B, const nmodf_ctx_t ctx);

FLINT_DLL void nmod_polydr_xgcd(nmod_polydr_t G, nmod_polydr_t S,
                nmod_polydr_t T, const nmod_polydr_t A, const nmod_polydr_t B,
                                                        const nmodf_ctx_t ctx);

FLINT_DLL void nmod_polydr_inv_series(nmod_polydr_t Qinv,
                        const nmod_polydr_t Q, slong n, const nmodf_ctx_t ctx);


/************** sparse poly with dense coefficients ************************/

typedef struct
{
    ulong exp;
    nmod_polydr_t vec;
} nmod_polydru_entry_struct;

typedef struct
{
    nmod_polydru_entry_struct * coeffs;
    slong length;
    slong alloc;
} nmod_polydru_struct;

typedef nmod_polydru_struct nmod_polydru_t[1];

NMOD_MPOLY_FACTOR_INLINE
void nmod_polydru_init(nmod_polydru_t A)
{
    A->coeffs = NULL;
    A->length = 0;
    A->alloc = 0;
}

FLINT_DLL void nmod_polydru_clear(nmod_polydru_t A);

FLINT_DLL void nmod_polydru_realloc(nmod_polydru_t A, slong len);

NMOD_MPOLY_FACTOR_INLINE
void nmod_polydru_fit_length(nmod_polydru_t A, slong len)
{
    if (len > A->alloc)
        nmod_polydru_realloc(A, len);
}

NMOD_MPOLY_FACTOR_INLINE
void nmod_polydru_swap(nmod_polydru_t A, nmod_polydru_t B)
{
    nmod_polydru_struct t = *B;
    *B = *A;
    *A = t;
}


/************** sparse poly with pair{polydr, polydr} coefficients ***********/

typedef struct
{
    ulong exp;
    nmod_polydr_t vec1;
    nmod_polydr_t vec2;
} nmod_polydr2u_entry_struct;

typedef struct
{
    nmod_polydr2u_entry_struct * coeffs;
    slong length;
    slong alloc;
} nmod_polydr2u_struct;

typedef nmod_polydr2u_struct nmod_polydr2u_t[1];

NMOD_MPOLY_FACTOR_INLINE
void nmod_polydr2u_init(nmod_polydr2u_t A)
{
    A->coeffs = NULL;
    A->length = 0;
    A->alloc = 0;
}

FLINT_DLL void nmod_polydr2u_clear(nmod_polydr2u_t A);

FLINT_DLL void nmod_polydr2u_realloc(nmod_polydr2u_t A, slong len);

NMOD_MPOLY_FACTOR_INLINE
void nmod_polydr2u_fit_length(nmod_polydr2u_t A, slong len)
{
    if (len > A->alloc)
        nmod_polydr2u_realloc(A, len);
}

NMOD_MPOLY_FACTOR_INLINE
void nmod_polydr2u_swap(nmod_polydr2u_t A, nmod_polydr2u_t B)
{
    nmod_polydr2u_struct t = *B;
    *B = *A;
    *A = t;
}


/******* dense bivariates (Z/nZ)[x,y] ****************************************/

typedef struct
{
    nmod_polydr_struct * coeffs;
    slong alloc;
    slong length;
} nmod_bpoly_struct;

typedef nmod_bpoly_struct nmod_bpoly_t[1];

NMOD_MPOLY_FACTOR_INLINE
void nmod_bpoly_init(nmod_bpoly_t A)
{
    A->coeffs = NULL;
    A->alloc = 0;
    A->length = 0;
}

FLINT_DLL void nmod_bpoly_clear(nmod_bpoly_t A);

NMOD_MPOLY_FACTOR_INLINE
void nmod_bpoly_swap(nmod_bpoly_t A, nmod_bpoly_t B)
{
    nmod_bpoly_struct t = *A;
    *A = *B;
    *B = t;
}

FLINT_DLL void nmod_bpoly_print_pretty(const nmod_bpoly_t A,
                                         const char * xvar, const char * yvar);


FLINT_DLL void nmod_bpoly_realloc(nmod_bpoly_t A, slong len);

NMOD_MPOLY_FACTOR_INLINE
void nmod_bpoly_fit_length(nmod_bpoly_t A, slong len)
{
    if (len > A->alloc)
        nmod_bpoly_realloc(A, len);
}


/*********** dense bivariates Fq[x,y] ****************************************/

typedef struct
{
    fq_nmod_poly_struct * coeffs;
    slong alloc;
    slong length;
} fq_nmod_bpoly_struct;

typedef fq_nmod_bpoly_struct fq_nmod_bpoly_t[1];

NMOD_MPOLY_FACTOR_INLINE
void fq_nmod_bpoly_init(fq_nmod_bpoly_t A, const fq_nmod_ctx_t ectx)
{
    A->coeffs = NULL;
    A->alloc = 0;
    A->length = 0;
}

FLINT_DLL void fq_nmod_bpoly_clear(fq_nmod_bpoly_t A, const fq_nmod_ctx_t ectx);

NMOD_MPOLY_FACTOR_INLINE
void fq_nmod_bpoly_swap(fq_nmod_bpoly_t A, fq_nmod_bpoly_t B)
{
    fq_nmod_bpoly_struct t = *A;
    *A = *B;
    *B = t;
}

FLINT_DLL void fq_nmod_bpoly_realloc(fq_nmod_bpoly_t A, slong len,
                                                      const fq_nmod_ctx_t ctx);

NMOD_MPOLY_FACTOR_INLINE
void fq_nmod_bpoly_fit_length(fq_nmod_bpoly_t A, slong len,
                                                       const fq_nmod_ctx_t ctx)
{
    if (len > A->alloc)
        fq_nmod_bpoly_realloc(A, len, ctx);
}


FLINT_DLL void fq_nmod_bpoly_print_pretty(const fq_nmod_bpoly_t A,
               const char * xvar, const char * yvar, const fq_nmod_ctx_t ectx);


/***************** dio solver **********************************************/

typedef struct {
    slong n;
    slong r;
    slong l;
    nmod_poly_struct * inv_prod_dbetas;
    nmod_poly_struct * dbetas;
    nmod_mpoly_struct * prod_mbetas;
    nmod_mpoly_struct * mbetas;
    nmod_mpoly_struct * deltas;
} nmod_disolve_struct;

typedef nmod_disolve_struct nmod_disolve_t[1];


void nmod_disolve_init(
    nmod_disolve_t I,
    slong l, slong r,
    const nmod_mpoly_struct * betas,
    const mp_limb_t * alpha,
    const nmod_mpoly_ctx_t ctx);

void nmod_disolve_clear(
    nmod_disolve_t I,
    const nmod_mpoly_ctx_t ctx);


int nmod_mfactor_disolve(
    flint_bitcnt_t bits,
    slong r,
    slong num,
    const nmod_mpoly_t t,
    const mp_limb_t * alpha,
    const slong * deg,
    const nmod_disolve_t I,
    const nmod_mpoly_ctx_t ctx);


int nmod_mfactor_lift(
    slong m,
    nmod_mpoly_factor_t lfac,
    const mp_limb_t * alpha,
    const nmod_mpoly_t A,
    const slong * degs,
    const nmod_mpoly_ctx_t ctx);




#ifdef __cplusplus
}
#endif

#endif

