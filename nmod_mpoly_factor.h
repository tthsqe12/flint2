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
    slong alloc;
    slong length;
} n_poly_struct;

typedef n_poly_struct n_poly_t[1];

NMOD_MPOLY_FACTOR_INLINE
void n_poly_init(n_poly_t A)
{
    A->length = 0;
    A->alloc = 0;
    A->coeffs = NULL;
}

NMOD_MPOLY_FACTOR_INLINE
void n_poly_init2(n_poly_t A, slong alloc)
{
    FLINT_ASSERT(alloc >= 0);
    A->length = 0;
    A->alloc = alloc;
    A->coeffs = alloc > 0 ? flint_malloc(alloc*sizeof(mp_limb_t)) : NULL;
}

FLINT_DLL void n_poly_clear(n_poly_t A);

FLINT_DLL void n_poly_realloc(n_poly_t A, slong len);

FLINT_DLL void n_poly_print_pretty(const n_poly_t A, const char * x);

NMOD_MPOLY_FACTOR_INLINE
void n_poly_fit_length(n_poly_t A, slong len)
{
    if (len > A->alloc)
        n_poly_realloc(A, len);
}


NMOD_MPOLY_FACTOR_INLINE
void nmod_poly_mock(nmod_poly_t a, const n_poly_t b, nmod_t mod)
{
    a->coeffs = b->coeffs;
    a->length = b->length;
    a->alloc = b->alloc;
    a->mod = mod;
}

NMOD_MPOLY_FACTOR_INLINE
void n_poly_mock(n_poly_t a, const nmod_poly_t b)
{
    a->coeffs = b->coeffs;
    a->length = b->length;
    a->alloc = b->alloc;
}


NMOD_MPOLY_FACTOR_INLINE
void n_poly_set(n_poly_t A, const n_poly_t B)
{
    n_poly_fit_length(A, B->length);
    flint_mpn_copyi(A->coeffs, B->coeffs, B->length);
    A->length = B->length;
}


NMOD_MPOLY_FACTOR_INLINE
void n_poly_swap(n_poly_t A, n_poly_t B)
{
    n_poly_struct t = *B;
    *B = *A;
    *A = t;
}

NMOD_MPOLY_FACTOR_INLINE
void _n_poly_normalise(n_poly_t A)
{
    while (A->length > 0 && A->coeffs[A->length - 1] == 0)
        A->length--;
}

NMOD_MPOLY_FACTOR_INLINE
slong n_poly_degree(const n_poly_t A)
{
    FLINT_ASSERT(A->length >= 0);
    return A->length - 1;
}

NMOD_MPOLY_FACTOR_INLINE
int n_poly_is_one(const n_poly_t A)
{
    return A->length == 1 && A->coeffs[0] == 1;
}

NMOD_MPOLY_FACTOR_INLINE
void n_poly_one(n_poly_t A)
{
    n_poly_fit_length(A, 1);
    A->length = 1;
    A->coeffs[0] = 1;
}

NMOD_MPOLY_FACTOR_INLINE
void n_poly_set_ui(n_poly_t A, mp_limb_t c)
{
    n_poly_fit_length(A, 1);
    A->coeffs[0] = c;
    A->length = (c != 0);
}

NMOD_MPOLY_FACTOR_INLINE
int n_poly_is_zero(const n_poly_t poly)
{
    return poly->length == 0;
}

NMOD_MPOLY_FACTOR_INLINE
void n_poly_zero(n_poly_t res)
{
    res->length = 0;
}

NMOD_MPOLY_FACTOR_INLINE
int n_poly_equal(const n_poly_t a, const n_poly_t b)
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
void n_poly_mod_make_monic(n_poly_t A, const n_poly_t B, nmod_t mod)
{
    FLINT_ASSERT(B->length > 0);
    n_poly_fit_length(A, B->length);
    A->length = B->length;
    _nmod_poly_make_monic(A->coeffs, B->coeffs, B->length, mod);
}

NMOD_MPOLY_FACTOR_INLINE
void n_poly_mod_taylor_shift(n_poly_t g, mp_limb_t c, nmod_t mod)
{
    _nmod_poly_taylor_shift(g->coeffs, c, g->length, mod);
}

NMOD_MPOLY_FACTOR_INLINE
ulong n_poly_get_coeff(const n_poly_t poly, slong j)
{
    return (j >= poly->length) ? 0 : poly->coeffs[j];
}

FLINT_DLL void n_poly_set_coeff(n_poly_t A, slong e, ulong c);

FLINT_DLL void n_poly_mod_set_coeff_ui(n_poly_t A, slong j, ulong c, nmod_t mod);

NMOD_MPOLY_FACTOR_INLINE
void n_poly_set_nmod_poly(n_poly_t a, const nmod_poly_t b)
{
    n_poly_fit_length(a, b->length);
    flint_mpn_copyi(a->coeffs, b->coeffs, b->length);
    a->length = b->length;
}

NMOD_MPOLY_FACTOR_INLINE
void nmod_poly_set_n_poly(nmod_poly_t a, const n_poly_t b)
{
    nmod_poly_fit_length(a, b->length);
    flint_mpn_copyi(a->coeffs, b->coeffs, b->length);
    a->length = b->length;
}

NMOD_MPOLY_FACTOR_INLINE
void n_poly_shift_left(n_poly_t A, const n_poly_t B, slong k)
{
    n_poly_fit_length(A, B->length + k);
    _nmod_poly_shift_left(A->coeffs, B->coeffs, B->length, k);
    A->length = B->length + k;
}

/* basic linear arithmetic */

NMOD_MPOLY_FACTOR_INLINE
void n_poly_mod_scalar_mul_nmod(n_poly_t A, const n_poly_t B, mp_limb_t c,
                                                                    nmod_t mod)
{
    n_poly_fit_length(A, B->length);
    _nmod_vec_scalar_mul_nmod(A->coeffs, B->coeffs, B->length, c, mod);
    A->length = B->length;
}

NMOD_MPOLY_FACTOR_INLINE
mp_limb_t n_poly_mod_evaluate_nmod(const n_poly_t A, mp_limb_t c, nmod_t mod)
{
    return _nmod_poly_evaluate_nmod(A->coeffs, A->length, c, mod);
}

NMOD_MPOLY_FACTOR_INLINE
void n_poly_mod_neg(n_poly_t A, const n_poly_t B, nmod_t mod)
{
    n_poly_fit_length(A, B->length);
    _nmod_vec_neg(A->coeffs, B->coeffs, B->length, mod);
    A->length = B->length;
}

NMOD_MPOLY_FACTOR_INLINE
void n_poly_mod_add(n_poly_t A, const n_poly_t B, const n_poly_t C, nmod_t mod)
{
    slong Alen = FLINT_MAX(B->length, C->length);
    n_poly_fit_length(A, Alen);
    _nmod_poly_add(A->coeffs, B->coeffs, B->length, C->coeffs, C->length, mod);
    A->length = Alen;
    _n_poly_normalise(A);
}

NMOD_MPOLY_FACTOR_INLINE
void n_poly_mod_sub(n_poly_t A, const n_poly_t B, const n_poly_t C, nmod_t mod)
{
    slong Alen = FLINT_MAX(B->length, C->length);
    n_poly_fit_length(A, Alen);
    _nmod_poly_sub(A->coeffs, B->coeffs, B->length, C->coeffs, C->length, mod);
    A->length = Alen;
    _n_poly_normalise(A);
}

/* no-nonsense quadratic arithmetic: no aliasing, no mod 1, ... */
NMOD_MPOLY_FACTOR_INLINE
void _n_poly_mod_mul(n_poly_t A, const n_poly_t B, const n_poly_t C,
                                                                    nmod_t mod)
{
    slong Blen = B->length;
    slong Clen = C->length;
    slong Alen = Blen + Clen - 1;

    FLINT_ASSERT(A != B);
    FLINT_ASSERT(A != C);

    if (Clen <= 0 || Blen <= 0)
    {
        A->length = 0;
        return;
    }

    n_poly_fit_length(A, Alen);
    A->length = Alen;

    if (Blen >= Clen)
        _nmod_poly_mul(A->coeffs, B->coeffs, Blen, C->coeffs, Clen, mod);
    else
        _nmod_poly_mul(A->coeffs, C->coeffs, Clen, B->coeffs, Blen, mod);
}

NMOD_MPOLY_FACTOR_INLINE
void _n_poly_mod_div(n_poly_t Q, const n_poly_t A, const n_poly_t B, nmod_t mod)
{
    const slong lenA = A->length, lenB = B->length;
    FLINT_ASSERT(lenB > 0);
    FLINT_ASSERT(Q != A && Q != B);
    if (lenA < lenB)
    {
        n_poly_zero(Q);
        return;
    }
    n_poly_fit_length(Q, lenA - lenB + 1);
    _nmod_poly_div(Q->coeffs, A->coeffs, lenA, B->coeffs, lenB, mod);
    Q->length = lenA - lenB + 1;
}

NMOD_MPOLY_FACTOR_INLINE
void _n_poly_mod_rem(n_poly_t R, const n_poly_t A, const n_poly_t B, nmod_t mod)
{
    const slong lenA = A->length, lenB = B->length;
    FLINT_ASSERT(R != A && R != B);
    FLINT_ASSERT(lenB > 0);
    if (lenA < lenB)
    {
        n_poly_set(R, A);
        return;
    }
    n_poly_fit_length(R, lenB - 1);
    _nmod_poly_rem(R->coeffs, A->coeffs, lenA, B->coeffs, lenB, mod);
    R->length = lenB - 1;
    _n_poly_normalise(R);
}

NMOD_MPOLY_FACTOR_INLINE
void _n_poly_mod_divrem(n_poly_t Q, n_poly_t R, const n_poly_t A,
                                                  const n_poly_t B, nmod_t mod)
{
    const slong lenA = A->length, lenB = B->length;

    FLINT_ASSERT(lenB > 0);
    FLINT_ASSERT(Q != A && Q != B);
    FLINT_ASSERT(R != A && R != B);

    if (lenA < lenB)
    {
        n_poly_set(R, A);
        n_poly_zero(Q);
        return;
    }

    n_poly_fit_length(Q, lenA - lenB + 1);
    n_poly_fit_length(R, lenB - 1);
    _nmod_poly_divrem(Q->coeffs, R->coeffs, A->coeffs, lenA, B->coeffs, lenB, mod);
    Q->length = lenA - lenB + 1;
    R->length = lenB - 1;
    _n_poly_normalise(R);
}

/* quadratic arithmetic dealing with aliasing, mod 1, and illegal input */

FLINT_DLL void n_poly_mod_mul(n_poly_t A, const n_poly_t B, const n_poly_t C,
                                                                   nmod_t mod);

FLINT_DLL void n_poly_mod_mullow(n_poly_t A, const n_poly_t B,
                                        const n_poly_t C, slong n, nmod_t mod);

FLINT_DLL void n_poly_mod_div(n_poly_t Q, const n_poly_t A, const n_poly_t B,
                                                                   nmod_t mod);

FLINT_DLL void n_poly_mod_rem(n_poly_t R, const n_poly_t A, const n_poly_t B,
                                                                   nmod_t mod);

FLINT_DLL void n_poly_mod_divrem(n_poly_t Q, n_poly_t R, const n_poly_t A,
                                                 const n_poly_t B, nmod_t mod);

FLINT_DLL int n_poly_mod_invmod(n_poly_t A, const n_poly_t B, const n_poly_t P,
                                                                   nmod_t mod);

FLINT_DLL void n_poly_mod_gcd(n_poly_t G, const n_poly_t A, const n_poly_t B,
                                                                   nmod_t mod);

FLINT_DLL void n_poly_mod_xgcd(n_poly_t G, n_poly_t S, n_poly_t T,
                               const n_poly_t A, const n_poly_t B, nmod_t mod);

FLINT_DLL void n_poly_mod_inv_series(n_poly_t Qinv, const n_poly_t Q, slong n,
                                                                   nmod_t mod);


/************** sparse univariate poly ************************/

typedef struct
{
    ulong exp;
    mp_limb_t coeff;
} n_polyu_term_struct;

typedef struct
{
    n_polyu_term_struct * terms;
    slong length;
    slong alloc;
} n_polyu_struct;

typedef n_polyu_struct n_polyu_t[1];

NMOD_MPOLY_FACTOR_INLINE
void n_polyu_init(n_polyu_t A)
{
    A->terms = NULL;
    A->length = 0;
    A->alloc = 0;
}

FLINT_DLL void n_polyu_clear(n_polyu_t A);

FLINT_DLL void n_polyu_realloc(n_polyu_t A, slong len);

NMOD_MPOLY_FACTOR_INLINE
void n_polyu_fit_length(n_polyu_t A, slong len)
{
    FLINT_ASSERT(A->alloc >= 0);
    if (len > A->alloc)
        n_polyu_realloc(A, len);
}

NMOD_MPOLY_FACTOR_INLINE
void n_polyu_swap(n_polyu_t A, n_polyu_t B)
{
    n_polyu_struct t = *B;
    *B = *A;
    *A = t;
}

/************** sparse poly with dense coefficients ************************/

typedef struct
{
    ulong exp;
    n_poly_t coeff;
} n_polyun_term_struct;

typedef struct
{
    n_polyun_term_struct * terms;
    slong length;
    slong alloc;
} n_polyun_struct;

typedef n_polyun_struct n_polyun_t[1];

NMOD_MPOLY_FACTOR_INLINE
void n_polyun_init(n_polyun_t A)
{
    A->terms = NULL;
    A->length = 0;
    A->alloc = 0;
}

FLINT_DLL void n_polyun_clear(n_polyun_t A);

FLINT_DLL void n_polyun_realloc(n_polyun_t A, slong len);

NMOD_MPOLY_FACTOR_INLINE
void n_polyun_fit_length(n_polyun_t A, slong len)
{
    if (len > A->alloc)
        n_polyun_realloc(A, len);
}

NMOD_MPOLY_FACTOR_INLINE
void n_polyun_swap(n_polyun_t A, n_polyun_t B)
{
    n_polyun_struct t = *B;
    *B = *A;
    *A = t;
}


/************** sparse poly with pair{polydr, polydr} coefficients ***********/

typedef struct
{
    ulong exp;
    n_poly_t vec1;
    n_poly_t vec2;
} n_poly2u_entry_struct;

typedef struct
{
    n_poly2u_entry_struct * coeffs;
    slong length;
    slong alloc;
} n_poly2u_struct;

typedef n_poly2u_struct n_poly2u_t[1];

NMOD_MPOLY_FACTOR_INLINE
void n_poly2u_init(n_poly2u_t A)
{
    A->coeffs = NULL;
    A->length = 0;
    A->alloc = 0;
}

FLINT_DLL void n_poly2u_clear(n_poly2u_t A);

FLINT_DLL void n_poly2u_realloc(n_poly2u_t A, slong len);

NMOD_MPOLY_FACTOR_INLINE
void n_poly2u_fit_length(n_poly2u_t A, slong len)
{
    if (len > A->alloc)
        n_poly2u_realloc(A, len);
}

NMOD_MPOLY_FACTOR_INLINE
void n_poly2u_swap(n_poly2u_t A, n_poly2u_t B)
{
    n_poly2u_struct t = *B;
    *B = *A;
    *A = t;
}


/******* dense bivariates (Z/nZ)[x,y] ****************************************/

typedef struct
{
    n_poly_struct * coeffs;
    slong alloc;
    slong length;
} n_bpoly_struct;

typedef n_bpoly_struct n_bpoly_t[1];

NMOD_MPOLY_FACTOR_INLINE
void n_bpoly_init(n_bpoly_t A)
{
    A->coeffs = NULL;
    A->alloc = 0;
    A->length = 0;
}

FLINT_DLL void n_bpoly_clear(n_bpoly_t A);

NMOD_MPOLY_FACTOR_INLINE
void n_bpoly_swap(n_bpoly_t A, n_bpoly_t B)
{
    n_bpoly_struct t = *A;
    *A = *B;
    *B = t;
}

FLINT_DLL void n_bpoly_print_pretty(const n_bpoly_t A,
                                        const char * xvar, const char * yvar);


FLINT_DLL void n_bpoly_realloc(n_bpoly_t A, slong len);

NMOD_MPOLY_FACTOR_INLINE
void n_bpoly_fit_length(n_bpoly_t A, slong len)
{
    if (len > A->alloc)
        n_bpoly_realloc(A, len);
}

NMOD_MPOLY_FACTOR_INLINE
void n_bpoly_zero(n_bpoly_t A)
{
    A->length = 0;
}

FLINT_DLL void n_bpoly_set_coeff(n_bpoly_t A, slong e0, slong e1, mp_limb_t c);


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

