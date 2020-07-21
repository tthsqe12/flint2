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

#ifdef __cplusplus
 extern "C" {
#endif

FLINT_DLL int fq_nmod_partial_fraction_coeffs(
    slong r,
    fq_nmod_poly_struct * out,
    const fq_nmod_poly_struct * f,
    const fq_nmod_ctx_t ectx);

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

FLINT_DLL void _fq_nmod_mpoly_get_lc(
    fq_nmod_mpoly_t c,
    const fq_nmod_mpoly_t A,
    const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL void _fq_nmod_mpoly_set_lc(
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
    const fq_nmod_mpoly_ctx_t ctx);

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

FLINT_DLL void fq_nmod_mpoly_get_fq_nmod_bpoly(
    fq_nmod_bpoly_t A,
    const fq_nmod_mpoly_t B,
    slong varx,
    slong vary,
    const fq_nmod_mpoly_ctx_t ctx);


FLINT_DLL void fq_nmod_mpoly_set_fq_nmod_bpoly(
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

FLINT_DLL void fq_nmod_bpoly_factor_lgprime(
    fq_nmod_poly_t c,
    fq_nmod_tpoly_t F,
    fq_nmod_bpoly_t B,
    const fq_nmod_ctx_t ctx);

#ifdef __cplusplus
}
#endif

#endif

