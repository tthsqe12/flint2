/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#ifndef FMPZ_MPOLY_FACTOR_H
#define FMPZ_MPOLY_FACTOR_H

#ifdef FMPZ_MPOLY_FACTOR_INLINES_C
#define FMPZ_MPOLY_FACTOR_INLINE FLINT_DLL
#else
#define FMPZ_MPOLY_FACTOR_INLINE static __inline__
#endif

#undef ulong
#define ulong ulongxx /* interferes with system includes */
#include <stdio.h>
#undef ulong
#include <gmp.h>
#define ulong mp_limb_t

#include "fmpz_mpoly.h"

#ifdef __cplusplus
 extern "C" {
#endif

typedef struct {
    fmpz_t content;
    fmpz_mpoly_struct * poly;
    fmpz * exp;
    slong length;
    slong alloc;
} fmpz_mpoly_factor_struct;

typedef fmpz_mpoly_factor_struct fmpz_mpoly_factor_t[1];

FLINT_DLL void fmpz_mpoly_factor_init(fmpz_mpoly_factor_t fac, const fmpz_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mpoly_factor_init2(fmpz_mpoly_factor_t fac, slong alloc, const fmpz_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mpoly_factor_realloc(fmpz_mpoly_factor_t fac, slong alloc, const fmpz_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mpoly_factor_fit_length(fmpz_mpoly_factor_t fac, slong len, const fmpz_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mpoly_factor_clear(fmpz_mpoly_factor_t fac, const fmpz_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mpoly_factor_set(fmpz_mpoly_factor_t res, const fmpz_mpoly_factor_t fac, const fmpz_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mpoly_factor_append_ui(fmpz_mpoly_factor_t f, const fmpz_mpoly_t A, ulong e, const fmpz_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mpoly_factor_append_fmpz(fmpz_mpoly_factor_t f, const fmpz_mpoly_t A, const fmpz_t e, const fmpz_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mpoly_factor_print_pretty(const fmpz_mpoly_factor_t f, const char ** vars, const fmpz_mpoly_ctx_t ctx);

FLINT_DLL int fmpz_mpoly_factor(fmpz_mpoly_factor_t fac, const fmpz_mpoly_t A, int full, const fmpz_mpoly_ctx_t ctx);

FMPZ_MPOLY_FACTOR_INLINE void fmpz_mpoly_factor_swap(fmpz_mpoly_factor_t A, fmpz_mpoly_factor_t B, const fmpz_mpoly_ctx_t ctx)
{
   fmpz_mpoly_factor_struct t = *A;
   *A = *B;
   *B = t;
}

FMPZ_MPOLY_FACTOR_INLINE void fmpz_mpoly_factor_set_fmpz(fmpz_mpoly_factor_t A,
                                    const fmpz_t b, const fmpz_mpoly_ctx_t ctx)
{
    A->length = 0;
    fmpz_set(A->content, b);
}

FMPZ_MPOLY_FACTOR_INLINE void fmpz_mpoly_factor_zero(fmpz_mpoly_factor_t A,
                                                    const fmpz_mpoly_ctx_t ctx)
{
    fmpz_zero(A->content);
    A->length = 0;
}

FLINT_DLL void fmpz_mpoly_factor_sort(fmpz_mpoly_factor_t f, const fmpz_mpoly_ctx_t ctx);

FLINT_DLL int fmpz_mpoly_factor_expand(fmpz_mpoly_t A,
                      const fmpz_mpoly_factor_t f, const fmpz_mpoly_ctx_t ctx);

FLINT_DLL int fmpz_mpoly_factor_fix_units(fmpz_mpoly_factor_t A,
                                                   const fmpz_mpoly_ctx_t ctx);

FLINT_DLL int fmpz_mpoly_factor_add(fmpz_mpoly_factor_t A,
                      const fmpz_mpoly_factor_t B, const fmpz_mpoly_factor_t C,
                                                   const fmpz_mpoly_ctx_t ctx);

FLINT_DLL int fmpz_mpoly_factor_sub(fmpz_mpoly_factor_t A,
                      const fmpz_mpoly_factor_t B, const fmpz_mpoly_factor_t C,
                                                   const fmpz_mpoly_ctx_t ctx);

FLINT_DLL int fmpz_mpoly_factor_mul(fmpz_mpoly_factor_t A,
                      const fmpz_mpoly_factor_t B, const fmpz_mpoly_factor_t C,
                                                   const fmpz_mpoly_ctx_t ctx);

FLINT_DLL int fmpz_mpoly_factor_div(fmpz_mpoly_factor_t A,
                      const fmpz_mpoly_factor_t B, const fmpz_mpoly_factor_t C,
                                                   const fmpz_mpoly_ctx_t ctx);

FLINT_DLL int fmpz_mpoly_factor_pow_fmpz(fmpz_mpoly_factor_t A,
      const fmpz_mpoly_factor_t B, const fmpz_t e, const fmpz_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mpoly_factor_scalar_mul_si(fmpz_mpoly_factor_t A,
             const fmpz_mpoly_factor_t B, slong c, const fmpz_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mpoly_factor_gen(fmpz_mpoly_factor_t A, slong v,
                                                   const fmpz_mpoly_ctx_t ctx);

FLINT_DLL int fmpz_mpoly_factor_is_same(const fmpz_mpoly_factor_t A,
                      const fmpz_mpoly_factor_t B, const fmpz_mpoly_ctx_t ctx);

FLINT_DLL int fmpz_mpoly_factor_set_str_pretty(fmpz_mpoly_factor_t poly,
              const char * str, const char** x_in, const fmpz_mpoly_ctx_t ctx);

FLINT_DLL int fmpz_mpoly_factor_bound_si(fmpz_t B, const fmpz_t A,
                                              const slong * degs, slong nvars);


FLINT_DLL void subset_first(fmpz_t a, slong n, slong r);

FLINT_DLL int subset_next(fmpz_t a, const fmpz_t b, slong n);

FLINT_DLL void subset_print(const fmpz_t a, slong n);

#ifdef __cplusplus
}
#endif

#endif

