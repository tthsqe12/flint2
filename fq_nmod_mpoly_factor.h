/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#ifndef FQ_NMOD_MPOLY_FACTOR_H
#define FQ_NMOD_MPOLY_FACTOR_H

#ifdef FQ_NMOD_MPOLY_FACTOR_INLINES_C
#define FQ_NMOD_MPOLY_FACTOR_INLINE FLINT_DLL
#else
#define FQ_NMOD_MPOLY_FACTOR_INLINE static __inline__
#endif

#undef ulong
#define ulong ulongxx /* interferes with system includes */
#include <stdio.h>
#undef ulong
#include <gmp.h>
#define ulong mp_limb_t

#include "fq_nmod_mpoly.h"
#include "nmod_mpoly_factor.h"

#ifdef __cplusplus
 extern "C" {
#endif

/*****************************************************************************/

FQ_NMOD_MPOLY_FACTOR_INLINE
void _fq_nmod_add(
    mp_limb_t * a,          /* length d */
    const mp_limb_t * b,    /* length d */
    const mp_limb_t * c,    /* length d */
    const fq_nmod_ctx_t ctx)
{
    slong d = ctx->modulus->length - 1;
    FLINT_ASSERT(d > 0);
    _nmod_vec_add(a, b, c, d, ctx->modulus->mod);
}

FQ_NMOD_MPOLY_FACTOR_INLINE
void _fq_nmod_sub(
    mp_limb_t * a,          /* length d */
    const mp_limb_t * b,    /* length d */
    const mp_limb_t * c,    /* length d */
    const fq_nmod_ctx_t ctx)
{
    slong d = ctx->modulus->length - 1;
    FLINT_ASSERT(d > 0);
    _nmod_vec_sub(a, b, c, d, ctx->modulus->mod);
}

FQ_NMOD_MPOLY_FACTOR_INLINE
void _fq_nmod_reduce2(
    mp_limb_t * a,          /* length d */
    mp_limb_t * b,          /* length 2*d - 1 */
    const fq_nmod_ctx_t ctx,
    mp_limb_t * t)          /* length d */
{
    slong i, k, d = ctx->modulus->length - 1;
    FLINT_ASSERT(d > 0);

    FLINT_ASSERT(a != b);

    if (ctx->sparse_modulus)
    {
        for (i = 2*d - 2; i >= d; i--)
        {
            if (b[i] == 0)
                continue;

            for (k = d - 1; k >= 0; k--)
            {
                FLINT_ASSERT(ctx->a[k] > 0);
                b[ctx->j[k] + i - d] = nmod_addmul(b[ctx->j[k] + i - d],
                                       b[i], ctx->mod.n - ctx->a[k], ctx->mod);
            }

            b[i] = 0;
        }

        for (i = 0; i < d; i++)
            a[i] = b[i];
    }
    else
    {
        slong blen = 2*d - 1;

        while (blen > d && b[blen - 1] == 0)
            blen--;

        _nmod_poly_divrem_newton_n_preinv(t, a, b, blen,
                                 ctx->modulus->coeffs, ctx->modulus->length,
                                 ctx->inv->coeffs, ctx->inv->length, ctx->mod);
    }
}

FQ_NMOD_MPOLY_FACTOR_INLINE
void _fq_nmod_mul(
    mp_limb_t * a,          /* length d */
    const mp_limb_t * b,    /* length d */
    const mp_limb_t * c,    /* length d */
    const fq_nmod_ctx_t ctx,
    mp_limb_t * t)          /* length 4*d */
{
    slong d = ctx->modulus->length - 1;
    FLINT_ASSERT(d > 0);
    _nmod_poly_mul(t + 0, b, d, c, d, ctx->modulus->mod);
    _fq_nmod_reduce2(a, t, ctx, t + 2*d);
}

FQ_NMOD_MPOLY_FACTOR_INLINE
void _fq_nmod_inv(
    mp_limb_t * a,
    const mp_limb_t * b,
    const fq_nmod_ctx_t ctx,
    mp_limb_t * t)  /* length d */
{
    slong blen;
    slong d = ctx->modulus->length - 1;
    FLINT_ASSERT(d > 0);

    while (blen > 0 && b[blen - 1] == 0)
        blen--;

    if (blen < 1)
    {
        flint_throw(FLINT_ERROR, "impossible inverse in _fq_nmod_inv");
    }
    else if (blen == 1)
    {
        a[0] = n_invmod(b[0], ctx->mod.n);
        _nmod_vec_zero(a + 1, d - 1);
    }
    else
    {
        if (1 != _nmod_poly_gcdinv(t, a, b, blen, ctx->modulus->coeffs, d + 1, ctx->mod))
        {
            flint_throw(FLINT_ERROR, "impossible inverse in _fq_nmod_inv");
        }

        if (t[0] != 1)
        {
            _nmod_vec_scalar_mul_nmod(a, a, d, n_invmod(t[0], ctx->mod.n), ctx->mod);
        }
    }
}


/*****************************************************************************/

FLINT_DLL int fq_nmod_mpoly_is_fq_nmod_poly(
    const fq_nmod_mpoly_t A,
    slong var,
    const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL int fq_nmod_mpoly_get_fq_nmod_poly(
    fq_nmod_poly_t A,
    const fq_nmod_mpoly_t B,
    slong var,
    const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL void fq_nmod_mpoly_fit_length_set_bits(
    fq_nmod_mpoly_t A,
    slong len,
    flint_bitcnt_t bits,
    const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL void _fq_nmod_mpoly_set_fq_nmod_poly(
    fq_nmod_mpoly_t A,
    flint_bitcnt_t Abits,
    const fq_nmod_struct * Bcoeffs,
    slong Blen,
    slong var,
    const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL void fq_nmod_mpoly_set_fq_nmod_poly(
    fq_nmod_mpoly_t A,
    const fq_nmod_poly_t B,
    slong var,
    const fq_nmod_mpoly_ctx_t ctx);

/*****************************************************************************/

void fq_nmod_mpoly_evaluate_one_fq_nmod(
    fq_nmod_mpoly_t A,
    const fq_nmod_mpoly_t B,
    slong var,
    const fq_nmod_t val,
    const fq_nmod_mpoly_ctx_t ctx);

typedef struct {
    fq_nmod_t constant;
    fq_nmod_mpoly_struct * poly;
    fmpz * exp;
    slong num;
    slong alloc;
} fq_nmod_mpoly_factor_struct;

typedef fq_nmod_mpoly_factor_struct fq_nmod_mpoly_factor_t[1];

FLINT_DLL void fq_nmod_mpoly_factor_init(fq_nmod_mpoly_factor_t f,
                                                const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL void fq_nmod_mpoly_factor_realloc(fq_nmod_mpoly_factor_t f,
                                   slong alloc, const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL void fq_nmod_mpoly_factor_fit_length(fq_nmod_mpoly_factor_t f,
                                     slong len, const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL void fq_nmod_mpoly_factor_clear(fq_nmod_mpoly_factor_t f,
                                                const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL void fq_nmod_mpoly_factor_set(fq_nmod_mpoly_factor_t a,
                const fq_nmod_mpoly_factor_t b, const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL void fq_nmod_mpoly_factor_print_pretty(const fq_nmod_mpoly_factor_t f,
                            const char ** vars, const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL void fq_nmod_mpoly_factor_append_ui(fq_nmod_mpoly_factor_t f,
              const fq_nmod_mpoly_t A, ulong e, const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL void fq_nmod_mpoly_factor_append_fmpz(fq_nmod_mpoly_factor_t f,
       const fq_nmod_mpoly_t A, const fmpz_t e, const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL int fq_nmod_mpoly_factor_squarefree(fq_nmod_mpoly_factor_t f,
                       const fq_nmod_mpoly_t A, const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL int fq_nmod_mpoly_factor(fq_nmod_mpoly_factor_t f,
                       const fq_nmod_mpoly_t A, const fq_nmod_mpoly_ctx_t ctx);

FQ_NMOD_MPOLY_FACTOR_INLINE
void fq_nmod_mpoly_factor_swap(fq_nmod_mpoly_factor_t A,
                       fq_nmod_mpoly_factor_t B, const fq_nmod_mpoly_ctx_t ctx)
{
   fq_nmod_mpoly_factor_struct t = *A;
   *A = *B;
   *B = t;
}

FQ_NMOD_MPOLY_FACTOR_INLINE
void fq_nmod_mpoly_factor_one(fq_nmod_mpoly_factor_t a,
                                                 const fq_nmod_mpoly_ctx_t ctx)
{
	fq_nmod_one(a->constant, ctx->fqctx);
	a->num = 0;
}

FLINT_DLL int fq_nmod_mpoly_factor_expand(fq_nmod_mpoly_t A,
                const fq_nmod_mpoly_factor_t f, const fq_nmod_mpoly_ctx_t ctx);


FQ_NMOD_MPOLY_FACTOR_INLINE
int fq_nmod_mpoly_factor_matches(const fq_nmod_mpoly_t a, const fq_nmod_mpoly_factor_t f, const fq_nmod_mpoly_ctx_t ctx)
{
    int matches;
    fq_nmod_mpoly_t t;
    fq_nmod_mpoly_init(t, ctx);
    fq_nmod_mpoly_factor_expand(t, f, ctx);
    matches = fq_nmod_mpoly_equal(t, a, ctx);
    fq_nmod_mpoly_clear(t, ctx);
    return matches;
}

FLINT_DLL int fq_nmod_mpoly_factor_fix_units(fq_nmod_mpoly_factor_t A,
                                                const fq_nmod_mpoly_ctx_t ctx);


FLINT_DLL void _fq_nmod_mpoly_get_lead0(
    fq_nmod_mpoly_t c,
    const fq_nmod_mpoly_t A,
    const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL void _fq_nmod_mpoly_set_lead0(
    fq_nmod_mpoly_t A,
    const fq_nmod_mpoly_t B,
    const fq_nmod_mpoly_t c,
    const fq_nmod_mpoly_ctx_t ctx);

/*****************************************************************************/

typedef struct
{
    fq_nmod_mpoly_struct * coeffs;
    slong alloc;
    slong length;
} fq_nmod_mpolyv_struct;

typedef fq_nmod_mpolyv_struct fq_nmod_mpolyv_t[1];

FQ_NMOD_MPOLY_FACTOR_INLINE
void fq_nmod_mpolyv_init(fq_nmod_mpolyv_t A, const fq_nmod_mpoly_ctx_t ctx)
{
    A->coeffs = NULL;
    A->alloc = 0;
    A->length = 0;
}

FQ_NMOD_MPOLY_FACTOR_INLINE
void fq_nmod_mpolyv_swap(fq_nmod_mpolyv_t A, fq_nmod_mpolyv_t B,
                                                 const fq_nmod_mpoly_ctx_t ctx)
{
   fq_nmod_mpolyv_struct t = *A;
   *A = *B;
   *B = t;
}

FLINT_DLL void fq_nmod_mpolyv_clear(fq_nmod_mpolyv_t A,
                                                const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL void fq_nmod_mpolyv_print_pretty(const fq_nmod_mpolyv_t poly,
                               const char ** x, const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL void fq_nmod_mpolyv_fit_length(fq_nmod_mpolyv_t A, slong length,
                                                const fq_nmod_mpoly_ctx_t ctx);

/*****************************************************************************/

FLINT_DLL int fq_nmod_mpoly_univar_content_mpoly(
    fq_nmod_mpoly_t g,
    const fq_nmod_mpoly_univar_t A,
    const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL void fq_nmod_mpoly_univar_divexact_mpoly(
    fq_nmod_mpoly_univar_t A,
    const fq_nmod_mpoly_t b,
    const fq_nmod_mpoly_ctx_t ctx);

/*****************************************************************************/

FLINT_DLL int fq_nmod_mpoly_factor_irred_smprime_default(
    fq_nmod_mpolyv_t fac,
    const fq_nmod_mpoly_t A,
    const fq_nmod_mpoly_ctx_t ctx,
    flint_rand_t state);

FLINT_DLL int fq_nmod_mpoly_factor_irred_lgprime_default(
    fq_nmod_mpolyv_t fac,
    const fq_nmod_mpoly_t A,
    const fq_nmod_mpoly_ctx_t ctx,
    flint_rand_t state);

/*****************************************************************************/

typedef struct {
    slong n;
    slong r;
    slong l;
    fq_nmod_poly_struct * inv_prod_dbetas;
    fq_nmod_poly_struct * dbetas;
    fq_nmod_mpoly_struct * prod_mbetas;
    fq_nmod_mpoly_struct * mbetas;
    fq_nmod_mpoly_struct * deltas;
} fq_nmod_mpoly_pfrac_struct;

typedef fq_nmod_mpoly_pfrac_struct fq_nmod_mpoly_pfrac_t[1];


FLINT_DLL void fq_nmod_mpoly_pfrac_init(
    fq_nmod_mpoly_pfrac_t I,
    slong l, slong r,
    const fq_nmod_mpoly_struct * betas,
    const fq_nmod_struct * alpha,
    const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL void fq_nmod_mpoly_pfrac_clear(
    fq_nmod_mpoly_pfrac_t I,
    const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL int fq_nmod_mpoly_pfrac(
    flint_bitcnt_t bits,
    slong r,
    slong num,
    const fq_nmod_mpoly_t t,
    const fq_nmod_struct * alpha,
    const slong * deg,
    const fq_nmod_mpoly_pfrac_t I,
    const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL int fq_nmod_mpoly_hlift(
    slong m,
    fq_nmod_mpoly_struct * f, /* length r */
    slong r,
    const fq_nmod_struct * alpha,
    const fq_nmod_mpoly_t A,
    const slong * degs,
    const fq_nmod_mpoly_ctx_t ctx);

/*****************************************************************************/

typedef struct
{
    fq_nmod_poly_struct * coeffs;
    slong alloc;
    slong length;
} fq_nmod_bpoly_struct;

typedef fq_nmod_bpoly_struct fq_nmod_bpoly_t[1];

FQ_NMOD_MPOLY_FACTOR_INLINE
void fq_nmod_bpoly_init(fq_nmod_bpoly_t A, const fq_nmod_ctx_t ectx)
{
    A->coeffs = NULL;
    A->alloc = 0;
    A->length = 0;
}

FLINT_DLL void fq_nmod_bpoly_clear(fq_nmod_bpoly_t A, const fq_nmod_ctx_t ctx);

FQ_NMOD_MPOLY_FACTOR_INLINE
void fq_nmod_bpoly_swap(fq_nmod_bpoly_t A, fq_nmod_bpoly_t B,
                                                       const fq_nmod_ctx_t ctx)
{
    fq_nmod_bpoly_struct t = *A;
    *A = *B;
    *B = t;
}

FLINT_DLL void fq_nmod_bpoly_realloc(fq_nmod_bpoly_t A, slong len,
                                                      const fq_nmod_ctx_t ctx);

FQ_NMOD_MPOLY_FACTOR_INLINE
void fq_nmod_bpoly_fit_length(fq_nmod_bpoly_t A, slong len,
                                                       const fq_nmod_ctx_t ctx)
{
    if (len > A->alloc)
        fq_nmod_bpoly_realloc(A, len, ctx);
}


FLINT_DLL void fq_nmod_bpoly_print_pretty(const fq_nmod_bpoly_t A,
               const char * xvar, const char * yvar, const fq_nmod_ctx_t ectx);

FQ_NMOD_MPOLY_FACTOR_INLINE
void fq_nmod_bpoly_zero(fq_nmod_bpoly_t A, const fq_nmod_ctx_t fqctx)
{
    A->length = 0;
}

FLINT_DLL void fq_nmod_bpoly_one(fq_nmod_bpoly_t A, const fq_nmod_ctx_t ctx);

FQ_NMOD_MPOLY_FACTOR_INLINE
void fq_nmod_bpoly_normalise(fq_nmod_bpoly_t A, const fq_nmod_ctx_t ctx)
{
    while (A->length > 0 && fq_nmod_poly_is_zero(A->coeffs + A->length - 1, ctx))
        A->length--;
}

FLINT_DLL int fq_nmod_bpoly_is_canonical(const fq_nmod_bpoly_t A,
                                                      const fq_nmod_ctx_t ctx);

FLINT_DLL void fq_nmod_bpoly_set(fq_nmod_bpoly_t A, const fq_nmod_bpoly_t B,
                                                      const fq_nmod_ctx_t ctx);

FLINT_DLL void fq_nmod_bpoly_set_coeff(fq_nmod_bpoly_t A, slong xi,
                         slong yi, const fq_nmod_t c, const fq_nmod_ctx_t ctx);

FLINT_DLL void fq_nmod_bpoly_set_poly_var0(fq_nmod_bpoly_t A,
                             const fq_nmod_poly_t B, const fq_nmod_ctx_t ectx);

FLINT_DLL void fq_nmod_bpoly_set_poly_var1(fq_nmod_bpoly_t A,
                             const fq_nmod_poly_t B, const fq_nmod_ctx_t ectx);

FLINT_DLL void fq_nmod_bpoly_get_coeff(fq_nmod_t c, const fq_nmod_bpoly_t A,
                                 slong xi, slong yi, const fq_nmod_ctx_t ectx);

FLINT_DLL void fq_nmod_bpoly_make_monic(fq_nmod_bpoly_t A, slong order,
                                                      const fq_nmod_ctx_t ctx);

FLINT_DLL void fq_nmod_bpoly_mul(
    fq_nmod_bpoly_t A,
    const fq_nmod_bpoly_t B,
    const fq_nmod_bpoly_t C,
    const fq_nmod_ctx_t ctx);

FLINT_DLL void fq_nmod_bpoly_mul_series(
    fq_nmod_bpoly_t A,
    const fq_nmod_bpoly_t B,
    const fq_nmod_bpoly_t C,
    slong order,
    const fq_nmod_ctx_t ectx);

FLINT_DLL void fq_nmod_bpoly_add_poly_shift(
    fq_nmod_bpoly_t A,
    const fq_nmod_poly_t B,
    slong yshift,
    const fq_nmod_ctx_t ctx);

FLINT_DLL void fq_nmod_bpoly_add(
    fq_nmod_bpoly_t A,
    const fq_nmod_bpoly_t B,
    const fq_nmod_bpoly_t C,
    const fq_nmod_ctx_t ctx);

FLINT_DLL void fq_nmod_bpoly_sub(
    fq_nmod_bpoly_t A,
    const fq_nmod_bpoly_t B,
    const fq_nmod_bpoly_t C,
    const fq_nmod_ctx_t ectx);

FLINT_DLL void fq_nmod_bpoly_derivative(
    fq_nmod_bpoly_t A,
    const fq_nmod_bpoly_t B,
    const fq_nmod_ctx_t ctx);

FLINT_DLL void fq_nmod_bpoly_divrem_series(
    fq_nmod_bpoly_t Q,
    fq_nmod_bpoly_t R,
    const fq_nmod_bpoly_t A,
    const fq_nmod_bpoly_t B,
    slong order,
    const fq_nmod_ctx_t ctx);

FLINT_DLL int fq_nmod_bpoly_divides(
    fq_nmod_bpoly_t Q,
    const fq_nmod_bpoly_t A,
    const fq_nmod_bpoly_t B,
    const fq_nmod_ctx_t ctx);

FLINT_DLL void fq_nmod_bpoly_make_primitive(
    fq_nmod_poly_t c,
    fq_nmod_bpoly_t A,
    const fq_nmod_ctx_t ctx);

FLINT_DLL void _fq_nmod_poly_taylor_shift_horner(
    fq_nmod_struct * poly,
    const fq_nmod_t c,
    slong n,
    const fq_nmod_ctx_t ctx);

FLINT_DLL void fq_nmod_poly_taylor_shift_horner(
    fq_nmod_poly_t g,
    const fq_nmod_poly_t f,
    const fq_nmod_t c,
    const fq_nmod_ctx_t ctx);

FLINT_DLL void fq_nmod_bpoly_taylor_shift_var1(
    fq_nmod_bpoly_t A,
    const fq_nmod_bpoly_t B,
    const fq_nmod_t alpha,
    const fq_nmod_ctx_t ctx);

FLINT_DLL void fq_nmod_mpoly_get_bpoly(
    fq_nmod_bpoly_t A,
    const fq_nmod_mpoly_t B,
    slong varx,
    slong vary,
    const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL void fq_nmod_mpoly_set_bpoly(
    fq_nmod_mpoly_t A,
    flint_bitcnt_t Abits,
    const fq_nmod_bpoly_t B,
    slong varx,
    slong vary,
    const fq_nmod_mpoly_ctx_t ctx);


typedef struct
{
    fq_nmod_bpoly_struct * coeffs;
    slong alloc;
    slong length;
} fq_nmod_tpoly_struct;

typedef fq_nmod_tpoly_struct fq_nmod_tpoly_t[1];

FQ_NMOD_MPOLY_FACTOR_INLINE
void fq_nmod_tpoly_init(fq_nmod_tpoly_t A, const fq_nmod_ctx_t ctx)
{
    A->coeffs = NULL;
    A->alloc = 0;
    A->length = 0;
}

FQ_NMOD_MPOLY_FACTOR_INLINE
void fq_nmod_tpoly_swap(fq_nmod_tpoly_t A, fq_nmod_tpoly_t B,
                                                       const fq_nmod_ctx_t ctx)
{
    fq_nmod_tpoly_struct t = *A;
    *A = *B;
    *B = t;
}

FLINT_DLL void fq_nmod_tpoly_fit_length(fq_nmod_tpoly_t A, slong len,
                                                      const fq_nmod_ctx_t ctx);

FLINT_DLL void fq_nmod_tpoly_clear(fq_nmod_tpoly_t A, const fq_nmod_ctx_t ctx);

FLINT_DLL int fq_nmod_bpoly_factor_smprime(
    fq_nmod_poly_t c,
    fq_nmod_tpoly_t F,
    fq_nmod_bpoly_t B,
    int allow_shift,
    const fq_nmod_ctx_t ctx);

FLINT_DLL int fq_nmod_bpoly_factor_lgprime(
    fq_nmod_poly_t c,
    fq_nmod_tpoly_t F,
    fq_nmod_bpoly_t B,
    const fq_nmod_ctx_t ctx,
    flint_rand_t state);

#ifdef __cplusplus
}
#endif

#endif

