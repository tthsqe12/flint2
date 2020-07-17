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

FLINT_DLL int fq_nmod_mpoly_univar_content_mpoly(
    fq_nmod_mpoly_t g,
    const fq_nmod_mpoly_univar_t A,
    const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL void fq_nmod_mpoly_univar_divexact_mpoly(
    fq_nmod_mpoly_univar_t A,
    const fq_nmod_mpoly_t b,
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


#ifdef __cplusplus
}
#endif

#endif
