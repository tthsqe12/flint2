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


/*****************************************************************************/

FLINT_DLL int nmod_mat_is_reduced(const nmod_mat_t N);

FLINT_DLL void nmod_mat_init_nullspace_tr(nmod_mat_t X, nmod_mat_t tmp);

FLINT_DLL int mpoly_is_poly(
    const ulong * Aexps,
    slong Alen,
    flint_bitcnt_t Abits,
    slong var,
    const mpoly_ctx_t mctx);

FLINT_DLL int nmod_mpoly_is_nmod_poly(
    const nmod_mpoly_t A,
    slong var,
    const nmod_mpoly_ctx_t ctx);

FLINT_DLL int nmod_mpoly_get_nmod_poly(
    nmod_poly_t A,
    const nmod_mpoly_t B,
    slong var,
    const nmod_mpoly_ctx_t ctx);

FLINT_DLL void nmod_mpoly_fit_length_set_bits(
    nmod_mpoly_t A,
    slong len,
    flint_bitcnt_t bits,
    const nmod_mpoly_ctx_t ctx);

FLINT_DLL void _nmod_mpoly_set_nmod_poly(
    nmod_mpoly_t A,
    flint_bitcnt_t Abits,
    const mp_limb_t * Bcoeffs,
    slong Blen,
    slong var,
    const nmod_mpoly_ctx_t ctx);

FLINT_DLL void nmod_mpoly_set_nmod_poly(
    nmod_mpoly_t A,
    const nmod_poly_t B,
    slong var,
    const nmod_mpoly_ctx_t ctx);

/*****************************************************************************/

typedef struct {
    mp_limb_t constant;
    nmod_mpoly_struct * poly;
    fmpz * exp;
    slong num;
    slong alloc;
} nmod_mpoly_factor_struct;

typedef nmod_mpoly_factor_struct nmod_mpoly_factor_t[1];

NMOD_MPOLY_FACTOR_INLINE
void nmod_mpoly_factor_init(nmod_mpoly_factor_t f, const nmod_mpoly_ctx_t ctx)
{
	f->constant = 1;
    f->poly  = NULL;
    f->exp   = NULL;
    f->num   = 0;
    f->alloc = 0;
}

FLINT_DLL void nmod_mpoly_factor_init2(nmod_mpoly_factor_t f, slong alloc,
                                                   const nmod_mpoly_ctx_t ctx);

FLINT_DLL void nmod_mpoly_factor_realloc(nmod_mpoly_factor_t f, slong alloc,
                                                   const nmod_mpoly_ctx_t ctx);

FLINT_DLL void nmod_mpoly_factor_fit_length(nmod_mpoly_factor_t f, slong len,
                                                   const nmod_mpoly_ctx_t ctx);

FLINT_DLL void nmod_mpoly_factor_clear(nmod_mpoly_factor_t f,
                                                   const nmod_mpoly_ctx_t ctx);

FLINT_DLL void nmod_mpoly_factor_append_ui(nmod_mpoly_factor_t f,
                    const nmod_mpoly_t A, ulong e, const nmod_mpoly_ctx_t ctx);

FLINT_DLL void nmod_mpoly_factor_append_fmpz(nmod_mpoly_factor_t f,
             const nmod_mpoly_t A, const fmpz_t e, const nmod_mpoly_ctx_t ctx);

FLINT_DLL void nmod_mpoly_factor_set(nmod_mpoly_factor_t f,
                      const nmod_mpoly_factor_t g, const nmod_mpoly_ctx_t ctx);

FLINT_DLL void nmod_mpoly_factor_print_pretty(const nmod_mpoly_factor_t f,
                               const char ** vars, const nmod_mpoly_ctx_t ctx);

FLINT_DLL int nmod_mpoly_factor_squarefree(nmod_mpoly_factor_t f,
                             const nmod_mpoly_t A, const nmod_mpoly_ctx_t ctx);

FLINT_DLL int nmod_mpoly_factor(nmod_mpoly_factor_t f, const nmod_mpoly_t A,
                                                   const nmod_mpoly_ctx_t ctx);

FLINT_DLL void nmod_mpoly_factor_sort(nmod_mpoly_factor_t f,
                                                   const nmod_mpoly_ctx_t ctx);

FLINT_DLL int nmod_mpoly_factor_expand(nmod_mpoly_t A,
                      const nmod_mpoly_factor_t f, const nmod_mpoly_ctx_t ctx);

NMOD_MPOLY_FACTOR_INLINE
int nmod_mpoly_factor_matches(const nmod_mpoly_t a, const nmod_mpoly_factor_t f, const nmod_mpoly_ctx_t ctx)
{
    int matches;
    nmod_mpoly_t t;
    nmod_mpoly_init(t, ctx);
    nmod_mpoly_factor_expand(t, f, ctx);
    matches = nmod_mpoly_equal(t, a, ctx);
    nmod_mpoly_clear(t, ctx);
    return matches;
}

FLINT_DLL int nmod_mpoly_factor_fix_units(nmod_mpoly_factor_t f,
                                                   const nmod_mpoly_ctx_t ctx);

NMOD_MPOLY_FACTOR_INLINE
void nmod_mpoly_factor_swap(nmod_mpoly_factor_t f, nmod_mpoly_factor_t g,
                                                    const nmod_mpoly_ctx_t ctx)
{
   nmod_mpoly_factor_struct t = *f;
   *f = *g;
   *g = t;
}

NMOD_MPOLY_FACTOR_INLINE
void nmod_mpoly_factor_one(nmod_mpoly_factor_t f, const nmod_mpoly_ctx_t ctx)
{
	f->constant = 1;
	f->num = 0;
}

FLINT_DLL void _nmod_mpoly_get_lead0(
    nmod_mpoly_t c,
    const nmod_mpoly_t A,
    const nmod_mpoly_ctx_t ctx);

FLINT_DLL void _nmod_mpoly_set_lead0(
    nmod_mpoly_t A,
    const nmod_mpoly_t B,
    const nmod_mpoly_t c,
    const nmod_mpoly_ctx_t ctx);

/*****************************************************************************/

FLINT_DLL void subset_first(fmpz_t a, slong n, slong r);

FLINT_DLL int subset_next(fmpz_t a, const fmpz_t b, slong n);

FLINT_DLL void subset_print(const fmpz_t a, slong n);

FLINT_DLL void subset_map_down(fmpz_t a, const fmpz_t b, const fmpz_t m);

FLINT_DLL int subset_fix(fmpz_t subset, slong len);

FLINT_DLL void tuple_print(fmpz * alpha, slong n);

FLINT_DLL void tuple_saturate(fmpz * alpha, slong n, slong m);

FLINT_DLL void tuple_next(fmpz * alpha, slong n);

/*****************************************************************************/

typedef struct
{
    nmod_mpoly_struct * coeffs;
    slong alloc;
    slong length;
} nmod_mpolyv_struct;

typedef nmod_mpolyv_struct nmod_mpolyv_t[1];

NMOD_MPOLY_FACTOR_INLINE
void nmod_mpolyv_init(nmod_mpolyv_t A, const nmod_mpoly_ctx_t ctx)
{
    A->coeffs = NULL;
    A->alloc = 0;
    A->length = 0;
}

NMOD_MPOLY_FACTOR_INLINE
void nmod_mpolyv_swap(nmod_mpolyv_t A, nmod_mpolyv_t B,
                                                    const nmod_mpoly_ctx_t ctx)
{
   nmod_mpolyv_struct t = *A;
   *A = *B;
   *B = t;
}

FLINT_DLL void nmod_mpolyv_clear(nmod_mpolyv_t A, const nmod_mpoly_ctx_t ctx);

FLINT_DLL void nmod_mpolyv_print_pretty(const nmod_mpolyv_t poly,
                                  const char ** x, const nmod_mpoly_ctx_t ctx);

FLINT_DLL void nmod_mpolyv_fit_length(nmod_mpolyv_t A, slong length,
                                                   const nmod_mpoly_ctx_t ctx);

/*****************************************************************************/

FLINT_DLL int nmod_mpoly_univar_content_mpoly(
    nmod_mpoly_t g,
    const nmod_mpoly_univar_t A,
    const nmod_mpoly_ctx_t ctx);

FLINT_DLL void nmod_mpoly_univar_divexact_mpoly(
    nmod_mpoly_univar_t A,
    const nmod_mpoly_t b,
    const nmod_mpoly_ctx_t ctx);

/*****************************************************************************/

FLINT_DLL int nmod_mpoly_factor_irred_smprime_default(
    nmod_mpolyv_t fac,
    const nmod_mpoly_t A,
    const nmod_mpoly_ctx_t ctx);

FLINT_DLL int nmod_mpoly_factor_irred_lgprime_default(
    nmod_mpolyv_t fac,
    const nmod_mpoly_t A,
    const nmod_mpoly_ctx_t ctx);

FLINT_DLL int nmod_irreducible_bivar_factors_smprime(
    nmod_mpoly_factor_t fac,
    const nmod_mpoly_t A,
    slong xvar,
    slong yvar,
    const nmod_mpoly_ctx_t ctx);

FLINT_DLL int nmod_irreducible_bivar_factors_lgprime(
    nmod_mpoly_factor_t fac,
    const nmod_mpoly_t A,
    slong xvar,
    slong yvar,
    const nmod_mpoly_ctx_t ctx);

/*****************************************************************************/

NMOD_MPOLY_FACTOR_INLINE
flint_bitcnt_t mpoly_gen_pow_exp_bits_required(
    slong v,
    ulong e,
    const mpoly_ctx_t mctx)
{
    return 1 + FLINT_BIT_COUNT(e); /* only lex and deg supported */
}


NMOD_MPOLY_FACTOR_INLINE
int z_add_checked(slong * a, slong b, slong c)
{
	ulong ahi, alo;
	add_ssaaaa(ahi, alo, 0, b, 0, c);
	*a = alo;
	return FLINT_SIGN_EXT(alo) != ahi;
}

NMOD_MPOLY_FACTOR_INLINE
mp_limb_t nmod_addmul(mp_limb_t a, mp_limb_t b, mp_limb_t c, nmod_t mod)
{
    NMOD_ADDMUL(a, b, c, mod);
    return a;
}


NMOD_MPOLY_FACTOR_INLINE
ulong pack_exp2(ulong e0, ulong e1)
{
    return (e0 << (1*(FLINT_BITS/2))) +
           (e1 << (0*(FLINT_BITS/2)));
}

NMOD_MPOLY_FACTOR_INLINE
ulong pack_exp3(ulong e0, ulong e1, ulong e2)
{
    return (e0 << (2*(FLINT_BITS/3))) +
           (e1 << (1*(FLINT_BITS/3))) +
           (e2 << (0*(FLINT_BITS/3)));
}

NMOD_MPOLY_FACTOR_INLINE
ulong extract_exp(ulong e, int idx, int nvars)
{
    return (e >> (idx*(FLINT_BITS/nvars))) &
            ((-UWORD(1)) >> (FLINT_BITS - FLINT_BITS/nvars));

}


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
    A->alloc = 0;       /* alloc >= 0 */
    A->coeffs = NULL;   /* alloc == 0  =>  coeffs == NULL */
}

NMOD_MPOLY_FACTOR_INLINE
void n_poly_init2(n_poly_t A, slong alloc)
{
    FLINT_ASSERT(A->alloc >= 0);
    A->length = 0;
    A->alloc = alloc;
    A->coeffs = NULL;
    if (alloc > 0)
        A->coeffs = (mp_limb_t *) flint_malloc(alloc*sizeof(mp_limb_t));
}

NMOD_MPOLY_FACTOR_INLINE
void n_poly_clear(n_poly_t A)
{
    FLINT_ASSERT(A->alloc != 0 || A->coeffs == NULL);
    if (A->alloc > 0)
        flint_free(A->coeffs);
}

FLINT_DLL void n_poly_realloc(n_poly_t A, slong len);

FLINT_DLL void n_poly_print_pretty(const n_poly_t A, const char * x);

FLINT_DLL int n_poly_mod_is_canonical(const n_poly_t A, nmod_t mod);

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

NMOD_MPOLY_FACTOR_INLINE
void n_poly_set_coeff_nonzero(n_poly_t A, slong j, ulong c)
{
    FLINT_ASSERT(c != 0);
    if (j >= A->length)
    {
        n_poly_fit_length(A, j + 1);
        flint_mpn_zero(A->coeffs + A->length, j - A->length);
        A->length = j + 1;
    }
    A->coeffs[j] = c;
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

NMOD_MPOLY_FACTOR_INLINE
void n_poly_shift_right(n_poly_t res, const n_poly_t poly, slong k)
{
    if (k >= poly->length)
    {
        res->length = 0;
    }
    else
    {
        const slong len = poly->length - k;
        n_poly_fit_length(res, len);
        _nmod_poly_shift_right(res->coeffs, poly->coeffs, len, k);
        res->length = len;
    }
}

NMOD_MPOLY_FACTOR_INLINE
void n_poly_truncate(n_poly_t poly, slong len)
{
    if (poly->length > len)
    {
        poly->length = len;
        _n_poly_normalise(poly);
    }
}

/* basic linear arithmetic */

NMOD_MPOLY_FACTOR_INLINE
void _n_poly_mod_scalar_mul_nmod(n_poly_t A, const n_poly_t B, mp_limb_t c,
                                                                    nmod_t mod)
{
    slong i;
    FLINT_ASSERT(B->length <= B->alloc);
    n_poly_fit_length(A, B->length);
    for (i = 0; i < B->length; i++)
        A->coeffs[i] = nmod_mul(B->coeffs[i], c, mod);
/*
    _nmod_vec_scalar_mul_nmod(A->coeffs, B->coeffs, B->length, c, mod);
*/
    A->length = B->length;
}

FLINT_DLL void n_poly_mod_scalar_mul_ui(n_poly_t A, const n_poly_t B,
                                                      mp_limb_t c, nmod_t ctx);

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

FLINT_DLL void n_poly_mod_add_ui(n_poly_t res, const n_poly_t poly, ulong c, nmod_t ctx);


NMOD_MPOLY_FACTOR_INLINE
void n_poly_mod_sub(n_poly_t A, const n_poly_t B, const n_poly_t C, nmod_t mod)
{
    slong Alen = FLINT_MAX(B->length, C->length);
    n_poly_fit_length(A, Alen);
    _nmod_poly_sub(A->coeffs, B->coeffs, B->length, C->coeffs, C->length, mod);
    A->length = Alen;
    _n_poly_normalise(A);
}

NMOD_MPOLY_FACTOR_INLINE
void n_poly_mod_product_roots_nmod_vec(n_poly_t A, mp_srcptr r, slong n, nmod_t mod)
{
    n_poly_fit_length(A, n + 1);
    A->length = n + 1;
    _nmod_poly_product_roots_nmod_vec(A->coeffs, r, n, mod);
}

FLINT_DLL void n_poly_mod_shift_left_scalar_addmul(n_poly_t A, slong k,
                                                      mp_limb_t c, nmod_t mod);

FLINT_DLL void n_poly_mod_addmul_linear(n_poly_t A, const n_poly_t B,
                     const n_poly_t C, mp_limb_t d1, mp_limb_t d0, nmod_t mod);

FLINT_DLL void n_poly_mod_eval2_pow(mp_limb_t * vp, mp_limb_t * vm,
                              const n_poly_t P, n_poly_t alphapow, nmod_t mod);

/* quadratic arithmetic without the bullshit: no aliasing, no mod 1, no throwing */

NMOD_MPOLY_FACTOR_INLINE
void _n_poly_mod_mul(n_poly_t A, const n_poly_t B, const n_poly_t C, nmod_t mod)
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

FLINT_DLL void n_poly_mod_pow(n_poly_t res, const n_poly_t poly, ulong e,
                                                                   nmod_t ctx);

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

FLINT_DLL void n_poly_mod_mulmod(n_poly_t res, const n_poly_t poly1,
                           const n_poly_t poly2, const n_poly_t f, nmod_t mod);

FLINT_DLL int n_poly_mod_invmod(n_poly_t A, const n_poly_t B, const n_poly_t P,
                                                                   nmod_t mod);

FLINT_DLL void n_poly_mod_gcd(n_poly_t G, const n_poly_t A, const n_poly_t B,
                                                                   nmod_t mod);

FLINT_DLL void n_poly_mod_xgcd(n_poly_t G, n_poly_t S, n_poly_t T,
                               const n_poly_t A, const n_poly_t B, nmod_t mod);

FLINT_DLL void n_poly_mod_inv_series(n_poly_t Qinv, const n_poly_t Q, slong n,
                                                                   nmod_t mod);

FLINT_DLL void n_poly_mod_div_series(n_poly_t Q, const n_poly_t A,
                                    const n_poly_t B, slong order, nmod_t ctx);

typedef struct {
    mp_limb_t constant;
    n_poly_struct * poly;
    slong * exp;
    slong num;
    slong alloc;
} n_poly_factor_struct;

typedef n_poly_factor_struct n_poly_factor_t[1];

FLINT_DLL void n_poly_factor_init(n_poly_factor_t f);

FLINT_DLL void n_poly_factor_fit_length(n_poly_factor_t f, slong len);

FLINT_DLL void n_poly_factor_clear(n_poly_factor_t f);

FLINT_DLL void n_poly_mod_factor(n_poly_factor_t f, const n_poly_t A, nmod_t mod);

FLINT_DLL void n_poly_factor_print_pretty(const n_poly_factor_t f, const char * x);

/*****************************************************************************/

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
void n_polyu_term_swap(n_polyu_term_struct * A, n_polyu_term_struct * B)
{
    n_polyu_term_struct T = *A;
    *A = *B;
    *B = T;
}
NMOD_MPOLY_FACTOR_INLINE
void n_polyu_swap(n_polyu_t A, n_polyu_t B)
{
    n_polyu_struct T = *B;
    *B = *A;
    *A = T;
}


FLINT_DLL void n_polyu_clear(n_polyu_t A);

FLINT_DLL void n_polyu_realloc(n_polyu_t A, slong len);

FLINT_DLL void n_polyu3_print_pretty(const n_polyu_t A, const char * var0,
                                         const char * var1, const char * var2);

FLINT_DLL void n_polyu3_degrees(slong * deg0, slong * deg1, slong * deg2,
                                                            const n_polyu_t A);

FLINT_DLL int n_polyu_mod_is_canonical(const n_polyu_t A, nmod_t mod);


/*****************************************************************************/

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
void n_polyun_term_swap(n_polyun_term_struct * A, n_polyun_term_struct * B)
{
    n_polyun_term_struct T = *A;
    *A = *B;
    *B = T;
}

NMOD_MPOLY_FACTOR_INLINE
void n_polyun_swap(n_polyun_t A, n_polyun_t B)
{
    n_polyun_struct t = *B;
    *B = *A;
    *A = t;
}

FLINT_DLL void n_polyun_clear(n_polyun_t A);

FLINT_DLL void n_polyun_realloc(n_polyun_t A, slong len);

FLINT_DLL void n_polyu2n_print_pretty(const n_polyun_t A, const char * var0,
                                      const char * var1, const char * varlast);

FLINT_DLL void n_polyu3n_print_pretty(const n_polyun_t A, const char * var0,
                   const char * var1, const char * var2, const char * varlast);

FLINT_DLL int n_polyun_mod_is_canonical(const n_polyun_t A, nmod_t mod);


/*****************************************************************************/

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

NMOD_MPOLY_FACTOR_INLINE
void n_bpoly_normalise(n_bpoly_t A)
{
    while (A->length > 0 && n_poly_is_zero(A->coeffs + A->length - 1))
        A->length--;
}

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

NMOD_MPOLY_FACTOR_INLINE
int n_bpoly_is_zero(const n_bpoly_t A)
{
    return A->length == 0;
}

FLINT_DLL void nmod_mpoly_get_bpoly(n_bpoly_t A, const nmod_mpoly_t B,
                           slong var0, slong var1, const nmod_mpoly_ctx_t ctx);

FLINT_DLL void nmod_mpoly_set_bpoly(nmod_mpoly_t A, flint_bitcnt_t Abits,
        const n_bpoly_t B, slong var0, slong var1, const nmod_mpoly_ctx_t ctx);

FLINT_DLL void _n_bpoly_set(n_bpoly_t A, const n_bpoly_t B);

NMOD_MPOLY_FACTOR_INLINE
void n_bpoly_set(n_bpoly_t A, const n_bpoly_t B)
{
    if (A != B)
        _n_bpoly_set(A, B);
}

FLINT_DLL void n_bpoly_one(n_bpoly_t A);

FLINT_DLL int n_bpoly_equal(const n_bpoly_t A, const n_bpoly_t B);

FLINT_DLL int n_bpoly_mod_is_canonical(const n_bpoly_t A, nmod_t mod);

FLINT_DLL void n_bpoly_set_coeff(n_bpoly_t A, slong e0, slong e1, mp_limb_t c);

FLINT_DLL void n_bpoly_set_coeff_nonzero(n_bpoly_t A, slong e0, slong e1,
                                                                  mp_limb_t c);

FLINT_DLL void n_bpoly_mod_derivative(n_bpoly_t A, const n_bpoly_t B,
                                                                   nmod_t ctx);

NMOD_MPOLY_FACTOR_INLINE
mp_limb_t n_bpoly_get_coeff(const n_bpoly_t A, slong e0, slong e1)
{
    if (e0 >= A->length)
        return 0;
    else
        return n_poly_get_coeff(A->coeffs + e0, e1);
}

NMOD_MPOLY_FACTOR_INLINE
slong n_bpoly_degree0(const n_bpoly_t A)
{
    return A->length - 1;
}

FLINT_DLL slong n_bpoly_degree1(const n_bpoly_t A);

FLINT_DLL void n_bpoly_set_poly_var1(n_bpoly_t A, const n_poly_t B);

FLINT_DLL void n_bpoly_set_poly_var0(n_bpoly_t A, const n_poly_t B);

FLINT_DLL void n_bpoly_mod_taylor_shift_var1(n_bpoly_t A, const n_bpoly_t B,
                                                      mp_limb_t c, nmod_t ctx);

FLINT_DLL void n_bpoly_mod_taylor_shift_var0(n_bpoly_t A, mp_limb_t c,
                                                                   nmod_t ctx);

FLINT_DLL void n_bpoly_mod_add(n_bpoly_t A, const n_bpoly_t B,
                                                const n_bpoly_t C, nmod_t ctx);

FLINT_DLL void n_bpoly_mod_sub(n_bpoly_t A, const n_bpoly_t B,
                                                const n_bpoly_t C, nmod_t ctx);

FLINT_DLL void n_bpoly_mod_make_primitive(n_poly_t g, n_bpoly_t A, nmod_t ctx);

FLINT_DLL void n_bpoly_mod_mul(n_bpoly_t A, const n_bpoly_t B,
                                                const n_bpoly_t C, nmod_t ctx);

FLINT_DLL int n_bpoly_mod_divides(n_bpoly_t Q, const n_bpoly_t A,
                                                const n_bpoly_t B, nmod_t ctx);

FLINT_DLL void n_bpoly_mod_mul_series(n_bpoly_t A, const n_bpoly_t B,
                                   const n_bpoly_t C, slong order, nmod_t ctx);

FLINT_DLL void n_bpoly_mod_divrem_series(n_bpoly_t Q, n_bpoly_t R,
                const n_bpoly_t A, const n_bpoly_t B, slong order, nmod_t ctx);

typedef struct
{
    n_bpoly_struct * coeffs;
    slong alloc;
    slong length;
} n_tpoly_struct;

typedef n_tpoly_struct n_tpoly_t[1];

NMOD_MPOLY_FACTOR_INLINE
void n_tpoly_init(n_tpoly_t A)
{
    A->coeffs = NULL;
    A->alloc = 0;
    A->length = 0;
}

NMOD_MPOLY_FACTOR_INLINE
void n_tpoly_swap(n_tpoly_t A, n_tpoly_t B)
{
    n_tpoly_struct t = *A;
    *A = *B;
    *B = t;
}

FLINT_DLL void n_tpoly_fit_length(n_tpoly_t A, slong len);

FLINT_DLL void n_tpoly_clear(n_tpoly_t A);

FLINT_DLL int n_bpoly_mod_factor_smprime(
    n_poly_t c,
    n_tpoly_t F,
    n_bpoly_t B,
    int allow_shift,
    nmod_t ctx);

FLINT_DLL void n_bpoly_mod_factor_lgprime(
    n_poly_t c,
    n_tpoly_t F,
    n_bpoly_t B,
    nmod_t ctx);


/*****************************************************************************/

FLINT_DLL int nmod_mpolyu_is_canonical(const nmod_mpolyu_t A,
                                                   const nmod_mpoly_ctx_t ctx);

FLINT_DLL void nmod_mpolyu3_print_pretty(
    const nmod_mpolyu_t A,
    const char * var0,
    const char * var1,
    const char * var2,
    const char ** vars,
    const nmod_mpoly_ctx_t ctx);

/*****************************************************************************/

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


/*****************************************************************************/

typedef struct {
    slong n;
    slong r;
    slong l;
    nmod_poly_struct * inv_prod_dbetas;
    nmod_poly_struct * dbetas;
    nmod_mpoly_struct * prod_mbetas;
    nmod_mpoly_struct * mbetas;
    nmod_mpoly_struct * deltas;
} nmod_mpoly_pfrac_struct;

typedef nmod_mpoly_pfrac_struct nmod_mpoly_pfrac_t[1];


FLINT_DLL void nmod_mpoly_pfrac_init(
    nmod_mpoly_pfrac_t I,
    slong l, slong r,
    const nmod_mpoly_struct * betas,
    const mp_limb_t * alpha,
    const nmod_mpoly_ctx_t ctx);

FLINT_DLL void nmod_mpoly_pfrac_clear(
    nmod_mpoly_pfrac_t I,
    const nmod_mpoly_ctx_t ctx);

FLINT_DLL int nmod_mpoly_pfrac(
    flint_bitcnt_t bits,
    slong r,
    slong num,
    const nmod_mpoly_t t,
    const mp_limb_t * alpha,
    const slong * deg,
    const nmod_mpoly_pfrac_t I,
    const nmod_mpoly_ctx_t ctx);

FLINT_DLL int nmod_mpoly_hlift(
    slong m,
    nmod_mpoly_struct * f,
    slong r,
    const mp_limb_t * alpha,
    const nmod_mpoly_t A,
    const slong * degs,
    const nmod_mpoly_ctx_t ctx);

FLINT_DLL int n_bpoly_mod_disolve(
    slong r,
    n_bpoly_struct * C,
    slong * C_deg1_bound,
    n_bpoly_t A,
    n_bpoly_struct * B,
    nmod_t mod);



FLINT_DLL int n_bpoly_mod_hensel_lift2(
    n_bpoly_t A, /* clobbered (shifted by alpha) */
    n_bpoly_t B0,
    n_bpoly_t B1,
    mp_limb_t alpha,
    slong degree_inner, /* required degree in x */
    nmod_t mod);

FLINT_DLL int n_bpoly_mod_hensel_lift(
    slong r,
    n_bpoly_t A, /* clobbered (shifted by alpha) */
    n_bpoly_struct * B,
    mp_limb_t alpha,
    slong degree_inner, /* required degree in x */
    nmod_t mod);

FLINT_DLL int n_polyu3_mod_hensel_lift(
    slong r,
    n_polyun_struct * BB,
    n_polyu_t A,
    n_polyu_struct * B,
    mp_limb_t beta,
    slong degree_inner, /* required degree in x */
    const nmodf_ctx_t ctx);



#ifdef __cplusplus
}
#endif

#endif

