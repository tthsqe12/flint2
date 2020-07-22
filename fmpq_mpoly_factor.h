/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#ifndef FMPQ_MPOLY_FACTOR_H
#define FMPQ_MPOLY_FACTOR_H

#ifdef FMPQ_MPOLY_FACTOR_INLINES_C
#define FMPQ_MPOLY_FACTOR_INLINE FLINT_DLL
#else
#define FMPQ_MPOLY_FACTOR_INLINE static __inline__
#endif

#undef ulong
#define ulong ulongxx /* interferes with system includes */
#include <stdio.h>
#undef ulong
#include <gmp.h>
#define ulong mp_limb_t

#include "flint/fmpq_mpoly.h"

#ifdef __cplusplus
 extern "C" {
#endif

typedef struct {
    fmpq_t constant;
    fmpq_mpoly_struct * poly;
    fmpz * exp;
    slong num;
    slong alloc;
} fmpq_mpoly_factor_struct;

typedef fmpq_mpoly_factor_struct fmpq_mpoly_factor_t[1];

FLINT_DLL void fmpq_mpoly_factor_init(fmpq_mpoly_factor_t f,
                                                   const fmpq_mpoly_ctx_t ctx);

FLINT_DLL void fmpq_mpoly_factor_realloc(fmpq_mpoly_factor_t f,
                                      slong alloc, const fmpq_mpoly_ctx_t ctx);

FLINT_DLL void fmpq_mpoly_factor_fit_length(fmpq_mpoly_factor_t f,
                                        slong len, const fmpq_mpoly_ctx_t ctx);

FLINT_DLL void fmpq_mpoly_factor_clear(fmpq_mpoly_factor_t f,
                                                   const fmpq_mpoly_ctx_t ctx);

FLINT_DLL int fmpq_mpoly_factor_squarefree(fmpq_mpoly_factor_t f,
                             const fmpq_mpoly_t A, const fmpq_mpoly_ctx_t ctx);

FLINT_DLL int fmpq_mpoly_factor(fmpq_mpoly_factor_t f, const fmpq_mpoly_t A,
                                                   const fmpq_mpoly_ctx_t ctx);

FLINT_DLL void _fmpq_mpoly_factor_swap_fmpz_mpoly_factor(fmpq_mpoly_factor_t f,
            fmpz_mpoly_factor_t g, const fmpq_t c, const fmpq_mpoly_ctx_t ctx);

#ifdef __cplusplus
}
#endif

#endif

