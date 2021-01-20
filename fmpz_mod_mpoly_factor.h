/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef FMPZ_MOD_MPOLY_FACTOR_H
#define FMPZ_MOD_MPOLY_FACTOR_H

#ifdef FMPZ_MOD_MPOLY_FACTOR_INLINES_C
#define FMPZ_MOD_MPOLY_FACTOR_INLINE FLINT_DLL
#else
#define FMPZ_MOD_MPOLY_FACTOR_INLINE static __inline__
#endif

#undef ulong
#define ulong ulongxx /* interferes with system includes */
#include <stdio.h>
#undef ulong
#include <gmp.h>
#define ulong mp_limb_t

#include "n_poly.h"
#include "nmod_mpoly_factor.h"
#include "fmpz_mod_mpoly.h"

#ifdef __cplusplus
 extern "C" {
#endif

/* fmpz_mod_poly extras ******************************************************/

FLINT_DLL void fmpz_mod_mpoly_set_fmpz_mod_poly(
    fmpz_mod_mpoly_t A,
    const fmpz_mod_poly_t B,
    slong var,
    const fmpz_mod_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mod_poly_eval2_pow(
    fmpz_t vp,
    fmpz_t vm,
    const fmpz_mod_poly_t P, 
    fmpz_mod_poly_t alphapow,
    const fmpz_mod_ctx_t ctx);

FLINT_DLL void fmpz_mod_poly_scalar_addmul_fmpz_mod(
    fmpz_mod_poly_t A,
    const fmpz_mod_poly_t B,
    const fmpz_mod_poly_t C,
    const fmpz_t d0,
    const fmpz_mod_ctx_t ctx);

FLINT_DLL void fmpz_mod_poly_addmul_linear(
    fmpz_mod_poly_t A,
    const fmpz_mod_poly_t B,
    const fmpz_mod_poly_t C,
    const fmpz_t d1, const fmpz_t d0,
    const fmpz_mod_ctx_t ctx);

FLINT_DLL void fmpz_mod_poly_shift_left_scalar_addmul_fmpz_mod(
    fmpz_mod_poly_t A,
    slong k,
    const fmpz_t c,
    const fmpz_mod_ctx_t ctx);

/* fmpz_mod_mat extras *******************************************************/

FLINT_DLL int fmpz_mod_mat_is_reduced(const fmpz_mod_mat_t N);

FLINT_DLL void fmpz_mod_mat_init_nullspace_tr(fmpz_mod_mat_t X, fmpz_mod_mat_t tmp, const fmpz_mod_ctx_t ctx);


/* type definitions **********************************************************/

typedef struct {
    fmpz_t constant;
    fmpz_mod_mpoly_struct * poly;
    fmpz * exp;
    slong num;
    slong alloc;
} fmpz_mod_mpoly_factor_struct;

typedef fmpz_mod_mpoly_factor_struct fmpz_mod_mpoly_factor_t[1];

FMPZ_MOD_MPOLY_FACTOR_INLINE
void fmpz_mod_mpoly_factor_init(fmpz_mod_mpoly_factor_t f, const fmpz_mod_mpoly_ctx_t ctx)
{
	fmpz_init_set_ui(f->constant, 1);
    f->poly  = NULL;
    f->exp   = NULL;
    f->num   = 0;
    f->alloc = 0;
}

FLINT_DLL void fmpz_mod_mpoly_factor_init2(fmpz_mod_mpoly_factor_t f, slong alloc,
                                                   const fmpz_mod_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mod_mpoly_factor_realloc(fmpz_mod_mpoly_factor_t f, slong alloc,
                                                   const fmpz_mod_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mod_mpoly_factor_fit_length(fmpz_mod_mpoly_factor_t f, slong len,
                                                   const fmpz_mod_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mod_mpoly_factor_clear(fmpz_mod_mpoly_factor_t f,
                                                   const fmpz_mod_mpoly_ctx_t ctx);

FMPZ_MOD_MPOLY_FACTOR_INLINE
slong fmpz_mod_mpoly_factor_length(const fmpz_mod_mpoly_factor_t f,
                                                    const fmpz_mod_mpoly_ctx_t ctx)
{
    return f->num;
}

FMPZ_MOD_MPOLY_FACTOR_INLINE
void fmpz_mod_mpoly_factor_get_constant_fmpz(fmpz_t c,
               const fmpz_mod_mpoly_factor_t f, const fmpz_mod_mpoly_ctx_t ctx)
{
    fmpz_set(c, f->constant);
}

FMPZ_MOD_MPOLY_FACTOR_INLINE
void fmpz_mod_mpoly_factor_get_base(fmpz_mod_mpoly_t p, const fmpz_mod_mpoly_factor_t f,
                                           slong i, const fmpz_mod_mpoly_ctx_t ctx)
{
    FLINT_ASSERT(i < (ulong) f->num);
    fmpz_mod_mpoly_set(p, f->poly + i, ctx);
}

FMPZ_MOD_MPOLY_FACTOR_INLINE
void fmpz_mod_mpoly_factor_swap_base(fmpz_mod_mpoly_t p, fmpz_mod_mpoly_factor_t f,
                                           slong i, const fmpz_mod_mpoly_ctx_t ctx)
{
    FLINT_ASSERT(i < (ulong) f->num);
    fmpz_mod_mpoly_swap(p, f->poly + i, ctx);
}

FMPZ_MOD_MPOLY_FACTOR_INLINE
slong fmpz_mod_mpoly_factor_get_exp_si(fmpz_mod_mpoly_factor_t f,
                                           slong i, const fmpz_mod_mpoly_ctx_t ctx)
{
    FLINT_ASSERT(i < (ulong) f->num);
    return fmpz_get_si(f->exp + i);
}

FLINT_DLL void fmpz_mod_mpoly_factor_set(fmpz_mod_mpoly_factor_t f,
              const fmpz_mod_mpoly_factor_t g, const fmpz_mod_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mod_mpoly_factor_print_pretty(const fmpz_mod_mpoly_factor_t f,
                           const char ** vars, const fmpz_mod_mpoly_ctx_t ctx);

FLINT_DLL int fmpz_mod_mpoly_factor_content(fmpz_mod_mpoly_factor_t f,
                     const fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_ctx_t ctx);

FLINT_DLL int fmpz_mod_mpoly_factor_squarefree(fmpz_mod_mpoly_factor_t f,
                     const fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_ctx_t ctx);

FLINT_DLL int fmpz_mod_mpoly_factor(fmpz_mod_mpoly_factor_t f,
                     const fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mod_mpoly_factor_sort(fmpz_mod_mpoly_factor_t f,
                                               const fmpz_mod_mpoly_ctx_t ctx);

FLINT_DLL int fmpz_mod_mpoly_factor_cmp(const fmpz_mod_mpoly_factor_t A,
              const fmpz_mod_mpoly_factor_t B, const fmpz_mod_mpoly_ctx_t ctx);

/******************************************************************************

   Internal functions (guaranteed to change without notice)

******************************************************************************/

FMPZ_MOD_MPOLY_FACTOR_INLINE slong _fmpz_mod_poly_degree(const fmpz_mod_poly_t a)
{
    return a->length - 1;
}

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

/*****************************************************************************/

FLINT_DLL int fmpz_mod_mpoly_factor_separable(fmpz_mod_mpoly_factor_t f,
            const fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_ctx_t ctx, int sep);

FLINT_DLL int fmpz_mod_mpoly_factor_expand(fmpz_mod_mpoly_t A,
              const fmpz_mod_mpoly_factor_t f, const fmpz_mod_mpoly_ctx_t ctx);

FMPZ_MOD_MPOLY_FACTOR_INLINE
int fmpz_mod_mpoly_factor_matches(const fmpz_mod_mpoly_t a, const fmpz_mod_mpoly_factor_t f, const fmpz_mod_mpoly_ctx_t ctx)
{
    int matches;
    fmpz_mod_mpoly_t t;
    fmpz_mod_mpoly_init(t, ctx);
    fmpz_mod_mpoly_factor_expand(t, f, ctx);
    matches = fmpz_mod_mpoly_equal(t, a, ctx);
    fmpz_mod_mpoly_clear(t, ctx);
    return matches;
}

FMPZ_MOD_MPOLY_FACTOR_INLINE
void fmpz_mod_mpoly_factor_append_fmpz_swap(fmpz_mod_mpoly_factor_t f,
            fmpz_mod_mpoly_t A, const fmpz_t e, const fmpz_mod_mpoly_ctx_t ctx)
{
    slong i = f->num;
    fmpz_mod_mpoly_factor_fit_length(f, i + 1, ctx);
    fmpz_mod_mpoly_swap(f->poly + i, A, ctx);
    fmpz_set(f->exp + i, e);
    f->num = i + 1;
}

FMPZ_MOD_MPOLY_FACTOR_INLINE
void fmpz_mod_mpoly_factor_append_fmpz(fmpz_mod_mpoly_factor_t f,
            fmpz_mod_mpoly_t A, const fmpz_t e, const fmpz_mod_mpoly_ctx_t ctx)
{
    slong i = f->num;
    fmpz_mod_mpoly_factor_fit_length(f, i + 1, ctx);
    fmpz_mod_mpoly_set(f->poly + i, A, ctx);
    fmpz_set(f->exp + i, e);
    f->num = i + 1;
}

FMPZ_MOD_MPOLY_FACTOR_INLINE
void fmpz_mod_mpoly_factor_append_ui(fmpz_mod_mpoly_factor_t f,
             const fmpz_mod_mpoly_t A, ulong e, const fmpz_mod_mpoly_ctx_t ctx)
{
    slong i = f->num;
    fmpz_mod_mpoly_factor_fit_length(f, i + 1, ctx);
    fmpz_mod_mpoly_set(f->poly + i, A, ctx);
    fmpz_set_ui(f->exp + i, e);
    f->num = i + 1;
}

FLINT_DLL int fmpz_mod_mpoly_factor_fix_units(fmpz_mod_mpoly_factor_t f,
                                               const fmpz_mod_mpoly_ctx_t ctx);

FMPZ_MOD_MPOLY_FACTOR_INLINE
void fmpz_mod_mpoly_factor_swap(fmpz_mod_mpoly_factor_t f,
                     fmpz_mod_mpoly_factor_t g, const fmpz_mod_mpoly_ctx_t ctx)
{
   fmpz_mod_mpoly_factor_struct t = *f;
   *f = *g;
   *g = t;
}

FMPZ_MOD_MPOLY_FACTOR_INLINE
void fmpz_mod_mpoly_factor_one(fmpz_mod_mpoly_factor_t f, const fmpz_mod_mpoly_ctx_t ctx)
{
	fmpz_one(f->constant);
	f->num = 0;
}

FLINT_DLL void _fmpz_mod_mpoly_get_lead0(
    fmpz_mod_mpoly_t c,
    const fmpz_mod_mpoly_t A,
    const fmpz_mod_mpoly_ctx_t ctx);

FLINT_DLL void _fmpz_mod_mpoly_set_lead0(
    fmpz_mod_mpoly_t A,
    const fmpz_mod_mpoly_t B,
    const fmpz_mod_mpoly_t c,
    const fmpz_mod_mpoly_ctx_t ctx);

/*****************************************************************************/

void _fmpz_mod_mpoly_factor_set_nmod_mpoly_factor(
    fmpz_mod_mpoly_factor_t f,
    const fmpz_mod_mpoly_ctx_t ctx,
    const nmod_mpoly_factor_t nf,
    const nmod_mpoly_ctx_t nctx);


/* stack *********************************************************************/

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

/* polyun ********************************************************************/

FMPZ_MOD_MPOLY_INLINE
ulong fmpz_mod_polyu1n_bidegree(const fmpz_mod_polyun_t A)
{
    ulong x, y;
    FLINT_ASSERT(A->length > 0);
    x = A->terms[0].exp;
    y = A->terms[0].coeff->length - 1;
    return (x << (FLINT_BITS/2)) + y;
}

FMPZ_MOD_MPOLY_INLINE
void fmpz_mod_polyun_swap(
    fmpz_mod_polyun_t A,
    fmpz_mod_polyun_t B)
{
    fmpz_mod_polyun_struct t = *A;
    *A = *B;
    *B = t;
}

FLINT_DLL int fmpz_mod_polyun_is_canonical(
    const fmpz_mod_polyun_t A,
    const fmpz_mod_ctx_t ctx);

FLINT_DLL void fmpz_mod_polyun_init(
    fmpz_mod_polyun_t A,
    const fmpz_mod_ctx_t ctx);

FLINT_DLL void fmpz_mod_polyun_clear(
    fmpz_mod_polyun_t A,
    const fmpz_mod_ctx_t ctx);

FLINT_DLL void fmpz_mod_polyun_realloc(
    fmpz_mod_polyun_t A,
    slong len,
    const fmpz_mod_ctx_t ctx);

FLINT_DLL void fmpz_mod_polyu2n_print_pretty(
    const fmpz_mod_polyun_t A,
    const char * var0,
    const char * var1,
    const char * varlast,
    const fmpz_mod_ctx_t ctx);

FLINT_DLL void fmpz_mod_polyun_set(
    fmpz_mod_polyun_t A,
    const fmpz_mod_polyun_t B,
    const fmpz_mod_ctx_t ctx);

FLINT_DLL void fmpz_mod_polyu3n_print_pretty(
    const fmpz_mod_polyun_t A,
    const char * var0,
    const char * var1,
    const char * var2,
    const char * varlast,
    const fmpz_mod_ctx_t ctx);

FLINT_DLL void fmpz_mod_polyu1n_print_pretty(
    const fmpz_mod_polyun_t A,
    const char * var0,
    const char * varlast,
    const fmpz_mod_ctx_t ctx);


FMPZ_MOD_MPOLY_INLINE
void fmpz_mod_polyun_fit_length(fmpz_mod_polyun_t A, slong len, const fmpz_mod_ctx_t ctx)
{
    if (len > A->alloc)
        fmpz_mod_polyun_realloc(A, len, ctx);
}

/* mpolyn ********************************************************************/

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

FMPZ_MOD_MPOLY_INLINE
const fmpz * fmpz_mod_mpolyn_leadcoeff(const fmpz_mod_mpolyn_t A)
{
    FLINT_ASSERT(A->length > 0);
    FLINT_ASSERT(A->coeffs[0].length > 0);
    return A->coeffs[0].coeffs + A->coeffs[0].length - 1;
}

FLINT_DLL int fmpz_mod_mpolyn_is_canonical(const fmpz_mod_mpolyn_t A,
                                               const fmpz_mod_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mod_mpolyn_content_poly(fmpz_mod_poly_t a,
                    const fmpz_mod_mpolyn_t B, const fmpz_mod_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mod_mpolyn_divexact_poly(fmpz_mod_mpolyn_t A,
                      const fmpz_mod_poly_t b, const fmpz_mod_mpoly_ctx_t ctx);

FLINT_DLL ulong nmod_mpolyn_bidegree(const nmod_mpolyn_t A);

FLINT_DLL ulong fmpz_mod_mpolyn_bidegree(const fmpz_mod_mpolyn_t A);

FLINT_DLL slong fmpz_mod_mpolyn_lastdeg(const fmpz_mod_mpolyn_t A,
                                               const fmpz_mod_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mod_mpolyn_clear(fmpz_mod_mpolyn_t A,
                                               const fmpz_mod_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mod_mpolyn_one(fmpz_mod_mpolyn_t A,
                                               const fmpz_mod_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mod_mpolyn_mul_poly(fmpz_mod_mpolyn_t A,
                            fmpz_mod_poly_t b, const fmpz_mod_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mod_mpolyn_scalar_mul_fmpz_mod(fmpz_mod_mpolyn_t A,
                               const fmpz_t c, const fmpz_mod_mpoly_ctx_t ctx);

FLINT_DLL int fmpz_mod_mpolyn_equal(const fmpz_mod_mpolyn_t A,
                   const fmpz_mod_mpolyn_t B, const fmpz_mod_mpoly_ctx_t ctx);

FLINT_DLL int fmpz_mod_mpolyn_gcd_brown_bivar(
     fmpz_mod_mpolyn_t G, fmpz_mod_mpolyn_t Abar, fmpz_mod_mpolyn_t Bbar,
     fmpz_mod_mpolyn_t A, fmpz_mod_mpolyn_t B, const fmpz_mod_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mod_mpolyn_print_pretty(const fmpz_mod_mpolyn_t poly,
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

FLINT_DLL void fmpz_mod_polyu1n_interp_reduce_2sm_poly(
    fmpz_mod_poly_t E,
    fmpz_mod_poly_t F,
    const fmpz_mod_polyun_t A,
    fmpz_mod_poly_t alphapow,
    const fmpz_mod_ctx_t ctx);

FLINT_DLL void fmpz_mod_polyu1n_interp_lift_2sm_poly(
    slong * lastdeg,
    fmpz_mod_polyun_t F,
    const fmpz_mod_poly_t A,
    const fmpz_mod_poly_t B,
    const fmpz_t alpha,
    const fmpz_mod_ctx_t ctx);

FLINT_DLL int fmpz_mod_polyu1n_interp_crt_2sm_poly(
    slong * lastdeg,
    fmpz_mod_polyun_t F,
    fmpz_mod_polyun_t T,
    const fmpz_mod_poly_t A,
    const fmpz_mod_poly_t B,
    const fmpz_mod_poly_t modulus,
    fmpz_mod_poly_t alphapow,
    const fmpz_mod_ctx_t ctx);

FLINT_DLL void fmpz_mod_mpolyn_interp_reduce_sm_poly(fmpz_mod_poly_t E,
                            const fmpz_mod_mpolyn_t A, const fmpz_t alpha,
                                               const fmpz_mod_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mod_mpolyn_interp_lift_sm_poly(fmpz_mod_mpolyn_t A,
                      const fmpz_mod_poly_t B, const fmpz_mod_mpoly_ctx_t ctx);

FLINT_DLL int fmpz_mod_mpolyn_interp_crt_sm_poly(slong * lastdeg_,
             fmpz_mod_mpolyn_t F, fmpz_mod_mpolyn_t T, const fmpz_mod_poly_t A,
                            const fmpz_mod_poly_t modulus, const fmpz_t alpha,
                                               const fmpz_mod_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mod_mpolyn_interp_reduce_2sm_mpolyn(
    fmpz_mod_mpolyn_t E,
    fmpz_mod_mpolyn_t F,
    fmpz_mod_mpolyn_t A,
    slong var,
    fmpz_mod_poly_t alphapow,
    const fmpz_mod_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mod_mpolyn_interp_lift_2sm_mpolyn(
    slong * lastdeg,
    fmpz_mod_mpolyn_t T,
    fmpz_mod_mpolyn_t A,
    fmpz_mod_mpolyn_t B,
    slong var,
    const fmpz_t alpha,
    const fmpz_mod_mpoly_ctx_t ctx);

FLINT_DLL int fmpz_mod_mpolyn_interp_crt_2sm_mpolyn(
    slong * lastdeg,
    fmpz_mod_mpolyn_t F,
    fmpz_mod_mpolyn_t T,
    fmpz_mod_mpolyn_t A,
    fmpz_mod_mpolyn_t B,
    slong var,
    fmpz_mod_poly_t modulus,
    fmpz_mod_poly_t alphapow,
    const fmpz_mod_mpoly_ctx_t ctx);


/*****************************************************************************/


FMPZ_MOD_MPOLY_FACTOR_INLINE
void fmpz_mod_polyu_swap(
    fmpz_mod_polyu_t A,
    fmpz_mod_polyu_t B)
{
    fmpz_mod_polyu_struct t = *A;
    *A = *B;
    *B = t;
}

FLINT_DLL void fmpz_mod_polyu_init(fmpz_mod_polyu_t A);

FLINT_DLL void fmpz_mod_polyu_clear(fmpz_mod_polyu_t A);

FLINT_DLL void fmpz_mod_polyu_realloc(fmpz_mod_polyu_t A, slong len);

FMPZ_MOD_MPOLY_FACTOR_INLINE
void fmpz_mod_polyu_fit_length(
    fmpz_mod_polyu_t a,
    slong len,
    const fmpz_mod_ctx_t ctx)
{
    if (len > a->alloc)
        fmpz_mod_polyu_realloc(a, len);
}

FLINT_DLL void fmpz_mod_polyu3_degrees(
    slong * deg0,
    slong * deg1,
    slong * deg2,
    const fmpz_mod_polyu_t A);

/*****************************************************************************/



/*****************************************************************************/

FMPZ_MOD_MPOLY_FACTOR_INLINE
void fmpz_mod_bpoly_init(fmpz_mod_bpoly_t A, const fmpz_mod_ctx_t ctx)
{
    A->coeffs = NULL;
    A->alloc = 0;
    A->length = 0;
}

FLINT_DLL void fmpz_mod_bpoly_clear(fmpz_mod_bpoly_t A, const fmpz_mod_ctx_t ctx);

FMPZ_MOD_MPOLY_FACTOR_INLINE
void fmpz_mod_bpoly_swap(fmpz_mod_bpoly_t A, fmpz_mod_bpoly_t B,
                                                      const fmpz_mod_ctx_t ctx)
{
    fmpz_mod_bpoly_struct t = *A;
    *A = *B;
    *B = t;
}


FMPZ_MOD_MPOLY_FACTOR_INLINE
void fmpz_mod_bpoly_get_coeff(fmpz_t c, const fmpz_mod_bpoly_t A,
                                  slong e0, slong e1, const fmpz_mod_ctx_t ctx)
{
    if (e0 >= A->length)
        fmpz_zero(c);
    else
        fmpz_mod_poly_get_coeff_fmpz(c, A->coeffs + e0, e1, ctx);
}


FMPZ_MOD_MPOLY_FACTOR_INLINE
slong fmpz_mod_bpoly_degree0(const fmpz_mod_bpoly_t A, const fmpz_mod_ctx_t ctx)
{
    return A->length - 1;
}

FMPZ_MOD_MPOLY_FACTOR_INLINE
void fmpz_mod_bpoly_normalise(fmpz_mod_bpoly_t A, const fmpz_mod_ctx_t ctx)
{
    while (A->length > 0 && fmpz_mod_poly_is_zero(A->coeffs + A->length - 1, ctx))
        A->length--;
}

FLINT_DLL int fmpz_mod_bpoly_equal(
    const fmpz_mod_bpoly_t A,
    const fmpz_mod_bpoly_t B,
    const fmpz_mod_ctx_t ctx);

FLINT_DLL void fmpz_mod_bpoly_set(
    fmpz_mod_bpoly_t A,
    const fmpz_mod_bpoly_t B,
    const fmpz_mod_ctx_t ctx);


FLINT_DLL void fmpz_mod_bpoly_set_poly_gen1(
    fmpz_mod_bpoly_t A,
    const fmpz_mod_poly_t B,
    const fmpz_mod_ctx_t ctx);


FLINT_DLL void fmpz_mod_bpoly_set_poly_gen0(
    fmpz_mod_bpoly_t A,
    const fmpz_mod_poly_t B,
    const fmpz_mod_ctx_t ctx);

FLINT_DLL void fmpz_mod_bpoly_one(fmpz_mod_bpoly_t A, const fmpz_mod_ctx_t ctx);

FMPZ_MOD_MPOLY_FACTOR_INLINE
int fmpz_mod_bpoly_is_one(const fmpz_mod_bpoly_t A, const fmpz_mod_ctx_t ctx)
{
    return A->length == 1 && fmpz_mod_poly_is_one(A->coeffs + 0, ctx);
}

FLINT_DLL slong fmpz_mod_bpoly_degree1(const fmpz_mod_bpoly_t A, const fmpz_mod_ctx_t ctx);

FLINT_DLL void fmpz_mod_bpoly_set_poly_gen1(fmpz_mod_bpoly_t A, const fmpz_mod_poly_t B, const fmpz_mod_ctx_t ctx);

FLINT_DLL void fmpz_mod_bpoly_set_poly_gen0(fmpz_mod_bpoly_t A, const fmpz_mod_poly_t B, const fmpz_mod_ctx_t ctx);

FLINT_DLL void fmpz_mod_bpoly_print_pretty(const fmpz_mod_bpoly_t A,
               const char * xvar, const char * yvar, const fmpz_mod_ctx_t ctx);

FLINT_DLL void fmpz_mod_bpoly_fit_length(fmpz_mod_bpoly_t A, slong len,
                                                     const fmpz_mod_ctx_t ctx);

FLINT_DLL void fmpz_mod_bpoly_set_coeff(fmpz_mod_bpoly_t A, slong xi, slong yi,
                                     const fmpz_t c, const fmpz_mod_ctx_t ctx);

FLINT_DLL void fmpz_mod_bpoly_zero(fmpz_mod_bpoly_t A, const fmpz_mod_ctx_t ctx);

FLINT_DLL void fmpz_mod_bpoly_reverse_vars(fmpz_mod_bpoly_t A,
                           const fmpz_mod_bpoly_t B, const fmpz_mod_ctx_t ctx);

FLINT_DLL void fmpz_mod_bpoly_taylor_shift_gen1(
    fmpz_mod_bpoly_t A,
    const fmpz_mod_bpoly_t B,
    const fmpz_t c,
    const fmpz_mod_ctx_t ctx);

FLINT_DLL void fmpz_mod_bpoly_sub(
    fmpz_mod_bpoly_t A,
    const fmpz_mod_bpoly_t B,
    const fmpz_mod_bpoly_t C,
    const fmpz_mod_ctx_t ctx);

FLINT_DLL void fmpz_mod_bpoly_add(
    fmpz_mod_bpoly_t A,
    const fmpz_mod_bpoly_t B,
    const fmpz_mod_bpoly_t C,
    const fmpz_mod_ctx_t ctx);

FLINT_DLL void fmpz_mod_bpoly_make_primitive(
    fmpz_mod_poly_t g,
    fmpz_mod_bpoly_t A,
    const fmpz_mod_ctx_t ctx);

FLINT_DLL void fmpz_mod_bpoly_mul(
    fmpz_mod_bpoly_t A,
    const fmpz_mod_bpoly_t B,
    const fmpz_mod_bpoly_t C,
    const fmpz_mod_ctx_t ctx);

FLINT_DLL void fmpz_mod_bpoly_mul_series(
    fmpz_mod_bpoly_t A,
    const fmpz_mod_bpoly_t B,
    const fmpz_mod_bpoly_t C,
    slong order,
    const fmpz_mod_ctx_t ctx);

FLINT_DLL void fmpz_mod_bpoly_divrem_series(
    fmpz_mod_bpoly_t Q,
    fmpz_mod_bpoly_t R,
    const fmpz_mod_bpoly_t A,
    const fmpz_mod_bpoly_t B,
    slong order,
    const fmpz_mod_ctx_t ctx);

FLINT_DLL int fmpz_mod_bpoly_divides(
    fmpz_mod_bpoly_t Q,
    const fmpz_mod_bpoly_t A,
    const fmpz_mod_bpoly_t B,
    const fmpz_mod_ctx_t ctx);

FLINT_DLL void fmpz_mod_bpoly_taylor_shift_gen0(
    fmpz_mod_bpoly_t A,
    const fmpz_t alpha,
    const fmpz_mod_ctx_t ctx);

FLINT_DLL void fmpz_mod_bpoly_derivative_gen0(
    fmpz_mod_bpoly_t A,
    const fmpz_mod_bpoly_t B,
    const fmpz_mod_ctx_t ctx);

FLINT_DLL void fmpz_mod_bpoly_make_monic_series(
    fmpz_mod_bpoly_t A,
    const fmpz_mod_bpoly_t B,
    slong order,
    const fmpz_mod_ctx_t ctx);


FMPZ_MOD_MPOLY_INLINE
void fmpz_mod_tpoly_init(fmpz_mod_tpoly_t A, const fmpz_mod_ctx_t ctx)
{
    A->coeffs = NULL;
    A->alloc = 0;
    A->length = 0;
}

FMPZ_MOD_MPOLY_INLINE
void fmpz_mod_tpoly_swap(fmpz_mod_tpoly_t A, fmpz_mod_tpoly_t B,
                                                      const fmpz_mod_ctx_t ctx)
{
    fmpz_mod_tpoly_struct t = *A;
    *A = *B;
    *B = t;
}

FLINT_DLL void fmpz_mod_tpoly_fit_length(fmpz_mod_tpoly_t A, slong len,
                                                     const fmpz_mod_ctx_t ctx);

FLINT_DLL void fmpz_mod_tpoly_clear(fmpz_mod_tpoly_t A, const fmpz_mod_ctx_t ctx);


FLINT_DLL void fmpz_mod_mpoly_get_fmpz_mod_bpoly(fmpz_mod_bpoly_t A,
                            const fmpz_mod_mpoly_t B, slong var0, slong var1,
                                               const fmpz_mod_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mod_mpoly_set_fmpz_mod_bpoly(fmpz_mod_mpoly_t A,
        flint_bitcnt_t Abits, const fmpz_mod_bpoly_t B, slong var0, slong var1,
                                               const fmpz_mod_mpoly_ctx_t ctx);

FLINT_DLL int fmpz_mod_bpoly_factor_smprime(fmpz_mod_poly_t c,
                    fmpz_mod_tpoly_t F, fmpz_mod_bpoly_t B, int allow_shift,
                                                     const fmpz_mod_ctx_t ctx);

/*****************************************************************************/

FLINT_DLL int fmpz_mod_bpoly_is_canonical(const fmpz_mod_bpoly_t A, const fmpz_mod_ctx_t ctx);

/*****************************************************************************/

FLINT_DLL int _fmpz_mod_zip_vand_solve(
    fmpz * coeffs,             /* in Fp: size mlength */
    const fmpz * monomials,    /* in Fp: size mlength */
    slong mlength,
    const fmpz * evals,        /* in Fp: size elength */
    slong elength,
    const fmpz * master,       /* in Fp: size mlength + 1 */
    fmpz * scratch,            /* in Fp: size mlength */
    const fmpz_mod_ctx_t ctx);

FLINT_DLL void _fmpz_mod_zip_eval_step(
    fmpz_t ev,
    fmpz * cur,            /* in Fp */
    const fmpz * inc,      /* in Fp */
    const fmpz * coeffs,   /* in Fp */
    slong length,
    const fmpz_mod_ctx_t ctx);


/*****************************************************************************/

typedef struct
{
    fmpz_mod_mpoly_struct * coeffs;
    slong alloc;
    slong length;
} fmpz_mod_mpolyv_struct;

typedef fmpz_mod_mpolyv_struct fmpz_mod_mpolyv_t[1];

FMPZ_MOD_MPOLY_FACTOR_INLINE
void fmpz_mod_mpolyv_init(fmpz_mod_mpolyv_t A, const fmpz_mod_mpoly_ctx_t ctx)
{
    A->coeffs = NULL;
    A->alloc = 0;
    A->length = 0;
}

FMPZ_MOD_MPOLY_FACTOR_INLINE
void fmpz_mod_mpolyv_swap(fmpz_mod_mpolyv_t A, fmpz_mod_mpolyv_t B,
                                                const fmpz_mod_mpoly_ctx_t ctx)
{
   fmpz_mod_mpolyv_struct t = *A;
   *A = *B;
   *B = t;
}

FLINT_DLL void fmpz_mod_mpolyv_clear(fmpz_mod_mpolyv_t A,
                                               const fmpz_mod_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mod_mpolyv_print_pretty(const fmpz_mod_mpolyv_t poly,
                              const char ** x, const fmpz_mod_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mod_mpolyv_fit_length(fmpz_mod_mpolyv_t A, slong length,
                                               const fmpz_mod_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mod_mpolyv_set_coeff(fmpz_mod_mpolyv_t A, slong i,
                                   fmpz_mod_mpoly_t c, const fmpz_mod_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mod_mpoly_to_mpolyv(fmpz_mod_mpolyv_t A, const fmpz_mod_mpoly_t B,
                        const fmpz_mod_mpoly_t xalpha, const fmpz_mod_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mod_mpoly_from_mpolyv(fmpz_mod_mpoly_t A, flint_bitcnt_t Abits,
                            const fmpz_mod_mpolyv_t B, const fmpz_mod_mpoly_t xalpha,
                                                   const fmpz_mod_mpoly_ctx_t ctx);

FLINT_DLL int _fmpz_mod_mpoly_vec_content_mpoly(fmpz_mod_mpoly_t g,
          const fmpz_mod_mpoly_struct * A, slong Alen, const fmpz_mod_mpoly_ctx_t ctx);

/*****************************************************************************/

FLINT_DLL int _fmpz_mod_mpoly_factor_separable(fmpz_mod_mpoly_factor_t f,
                    const fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_ctx_t ctx, int sep);

FLINT_DLL int fmpz_mod_mpoly_factor_lcc_wang(fmpz_mod_mpoly_struct * lc_divs,
             const fmpz_mod_mpoly_factor_t lcAfac, const fmpz_mod_poly_t Auc,
             const fmpz_mod_bpoly_struct * Auf, slong r,
           const fmpz_mod_poly_struct * alpha, const fmpz_mod_mpoly_ctx_t ctx);

FLINT_DLL int fmpz_mod_mpoly_factor_irred_smprime_zassenhaus(fmpz_mod_mpolyv_t fac,
         const fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_ctx_t ctx, flint_rand_t state);

FLINT_DLL int fmpz_mod_mpoly_factor_irred_smprime_wang(fmpz_mod_mpolyv_t fac,
       const fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_factor_t lcAfac,
       const fmpz_mod_mpoly_t lcA, const fmpz_mod_mpoly_ctx_t ctx, flint_rand_t state);

FLINT_DLL int fmpz_mod_mpoly_factor_irred_smprime_zippel(fmpz_mod_mpolyv_t fac,
       const fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_factor_t lcAfac,
       const fmpz_mod_mpoly_t lcA, const fmpz_mod_mpoly_ctx_t ctx, flint_rand_t state);

/*****************************************************************************/

FLINT_DLL void fmpz_mod_mpoly_compression_do(fmpz_mod_mpoly_t L,
                 const fmpz_mod_mpoly_ctx_t Lctx, fmpz * Acoeffs, slong Alen,
                                                        mpoly_compression_t M);

FLINT_DLL void fmpz_mod_mpoly_compression_undo(fmpz_mod_mpoly_t A,
             flint_bitcnt_t Abits, const fmpz_mod_mpoly_ctx_t Actx, fmpz_mod_mpoly_t L,
                           const fmpz_mod_mpoly_ctx_t Lctx, mpoly_compression_t M);

/*****************************************************************************/

FLINT_DLL int fmpz_mod_mpolyu_is_canonical(const fmpz_mod_mpolyu_t A,
                                               const fmpz_mod_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mod_mpolyu3_print_pretty(const fmpz_mod_mpolyu_t A,
                    const char * var0, const char * var1, const char * var2,
                               const char ** vars, const fmpz_mod_mpoly_ctx_t ctx);

/*****************************************************************************/

typedef struct {
    fmpz * M;
    fmpz * T;
    fmpz * Q;
    fmpz * array;
    slong alloc;
    slong d;
    slong radix;
    fmpz_t w;
} fmpz_mod_eval_interp_struct;

typedef fmpz_mod_eval_interp_struct fmpz_mod_eval_interp_t[1];


FLINT_DLL void fmpz_mod_eval_interp_init(fmpz_mod_eval_interp_t E);

FLINT_DLL void fmpz_mod_eval_interp_clear(fmpz_mod_eval_interp_t E);

FLINT_DLL int fmpz_mod_eval_interp_set_degree_modulus(
    fmpz_mod_eval_interp_t E,
    slong deg,
    const fmpz_mod_ctx_t ctx);

FMPZ_MOD_MPOLY_FACTOR_INLINE
slong fmpz_mod_eval_interp_eval_length(fmpz_mod_eval_interp_t E)
{
    return 1 + E->radix*E->d;
}

FLINT_DLL void fmpz_mod_eval_interp_to_coeffs_poly(
    fmpz_mod_poly_t a,
    const fmpz_mod_poly_t v,
    fmpz_mod_eval_interp_t E,
    const fmpz_mod_ctx_t ctx);

FLINT_DLL void fmpz_mod_eval_interp_from_coeffs_poly(
    fmpz_mod_poly_t v,
    const fmpz_mod_poly_t a,
    fmpz_mod_eval_interp_t E,
    const fmpz_mod_ctx_t ctx);

FMPZ_MOD_MPOLY_FACTOR_INLINE void
fmpz_mod_evals_zero(fmpz_mod_poly_t a)
{
    a->length = 0;
}

FLINT_DLL void fmpz_mod_evals_add_inplace(fmpz_mod_poly_t a, fmpz_mod_poly_t b, slong len,
                                                                   const fmpz_mod_ctx_t ctx);

FLINT_DLL void fmpz_mod_evals_mul(fmpz_mod_poly_t a, fmpz_mod_poly_t b, fmpz_mod_poly_t c, slong len,
                                                                   const fmpz_mod_ctx_t ctx);

FLINT_DLL void fmpz_mod_evals_addmul(fmpz_mod_poly_t a, fmpz_mod_poly_t b, fmpz_mod_poly_t c, slong len,
                                                                   const fmpz_mod_ctx_t ctx);

FLINT_DLL void fmpz_mod_evals_fmma(fmpz_mod_poly_t a, fmpz_mod_poly_t b, fmpz_mod_poly_t c,
                                fmpz_mod_poly_t d, fmpz_mod_poly_t e, slong len, const fmpz_mod_ctx_t ctx);

/*****************************************************************************/

typedef struct {
    flint_bitcnt_t bits;
    slong w;
    slong r;
    fmpz_mod_poly_struct * inv_prod_dbetas;
    fmpz_mod_mpoly_struct * inv_prod_dbetas_mvar;
    fmpz_mod_poly_struct * dbetas;
    fmpz_mod_mpoly_struct * dbetas_mvar;
    fmpz_mod_mpoly_struct * prod_mbetas;
    fmpz_mod_mpolyv_struct * prod_mbetas_coeffs;
    fmpz_mod_mpoly_struct * mbetas;
    fmpz_mod_mpoly_struct * deltas;
    fmpz_mod_mpoly_struct * xalpha;
    fmpz_mod_mpoly_struct * q;
    fmpz_mod_mpoly_geobucket_struct * G;
    fmpz_mod_mpoly_struct * qt;
    fmpz_mod_mpoly_struct * newt;
    fmpz_mod_mpolyv_struct * delta_coeffs;
    fmpz_mod_mpoly_t T;
    fmpz_mod_mpoly_t Q;
    fmpz_mod_mpoly_t R;
} fmpz_mod_mpoly_pfrac_struct;

typedef fmpz_mod_mpoly_pfrac_struct fmpz_mod_mpoly_pfrac_t[1];


FLINT_DLL int fmpz_mod_mpoly_pfrac_init(fmpz_mod_mpoly_pfrac_t I, flint_bitcnt_t bits,
                         slong l, slong r, const fmpz_mod_mpoly_struct * betas, 
                          const fmpz * alpha, const fmpz_mod_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mod_mpoly_pfrac_clear(fmpz_mod_mpoly_pfrac_t I,
                                                const fmpz_mod_mpoly_ctx_t ctx);

FLINT_DLL int fmpz_mod_mpoly_pfrac(slong r, fmpz_mod_mpoly_t t, const slong * deg,
                     fmpz_mod_mpoly_pfrac_t I, const fmpz_mod_mpoly_ctx_t ctx);

FLINT_DLL int fmpz_mod_mpoly_hlift(slong m, fmpz_mod_mpoly_struct * f, slong r,
            const fmpz * alpha, const fmpz_mod_mpoly_t A, const slong * degs,
                                               const fmpz_mod_mpoly_ctx_t ctx);

FLINT_DLL int fmpz_mod_bpoly_pfrac(slong r, fmpz_mod_bpoly_struct * C,
          slong * C_deg1_bound, fmpz_mod_bpoly_t A, fmpz_mod_bpoly_struct * B,
                                                     const fmpz_mod_ctx_t ctx);

FLINT_DLL int fmpz_mod_bpoly_hlift2_cubic(
    fmpz_mod_bpoly_t A, /* clobbered (shifted by alpha) */
    fmpz_mod_bpoly_t B0,
    fmpz_mod_bpoly_t B1,
    const fmpz_t alpha,
    slong degree_inner, /* required degree in x */
    const fmpz_mod_ctx_t ctx,
    fmpz_mod_eval_interp_t E,
    fmpz_mod_poly_bpoly_stack_t St);

FLINT_DLL int fmpz_mod_bpoly_hlift2(fmpz_mod_bpoly_t A, fmpz_mod_bpoly_t B0,
                  fmpz_mod_bpoly_t B1, const fmpz_t alpha, slong degree_inner,
                     const fmpz_mod_ctx_t ctx, fmpz_mod_poly_bpoly_stack_t St);

FLINT_DLL int fmpz_mod_bpoly_hlift(slong r, fmpz_mod_bpoly_t A,
           fmpz_mod_bpoly_struct * B, const fmpz_t alpha, slong degree_inner,
                     const fmpz_mod_ctx_t ctx, fmpz_mod_poly_bpoly_stack_t St);

FLINT_DLL int fmpz_mod_polyu3_hlift(slong r, fmpz_mod_polyun_struct * BB,
              fmpz_mod_polyu_t A, fmpz_mod_polyu_struct * B, const fmpz_t beta,
                                 slong degree_inner, const fmpz_mod_ctx_t ctx);

FLINT_DLL int fmpz_mod_mpoly_hlift_zippel(slong m, fmpz_mod_mpoly_struct * B, slong r,
            const fmpz * alpha, const fmpz_mod_mpoly_t A, const slong * degs,
                           const fmpz_mod_mpoly_ctx_t ctx, flint_rand_t state);

FLINT_DLL int fmpz_mod_mpoly_factor_algo(fmpz_mod_mpoly_factor_t f,
                    const fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_ctx_t ctx,
                                                            unsigned int algo);

FLINT_DLL int fmpz_mod_mpoly_factor_zassenhaus(fmpz_mod_mpoly_factor_t f,
                     const fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_ctx_t ctx);

FLINT_DLL int fmpz_mod_mpoly_factor_wang(fmpz_mod_mpoly_factor_t f,
                     const fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_ctx_t ctx);

FLINT_DLL int fmpz_mod_mpoly_factor_zippel(fmpz_mod_mpoly_factor_t f,
                             const fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_ctx_t ctx);

FLINT_DLL int _fmpz_mod_mpoly_evaluate_rest_fmpz_mod_poly(fmpz_mod_poly_struct * E,
                slong * starts, slong * ends, slong * stops, ulong * es,
                const fmpz * Acoeffs, const ulong * Aexps, slong Alen, slong var,
                const fmpz_mod_poly_struct * alphas,
                const slong * offsets, const slong * shifts,
                   slong N, ulong mask, slong nvars, const fmpz_mod_ctx_t ctx);

FLINT_DLL void _fmpz_mod_mpoly_eval_rest_to_fmpz_mod_bpoly(fmpz_mod_bpoly_t E,
             const fmpz_mod_mpoly_t A, const fmpz_mod_poly_struct * alphabetas,
                                                const fmpz_mod_mpoly_ctx_t ctx);

FLINT_DLL void _fmpz_mod_mpoly_set_fmpz_mod_bpoly_var1_zero(fmpz_mod_mpoly_t A,
                     flint_bitcnt_t Abits, const fmpz_mod_bpoly_t B, slong var,
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

FLINT_DLL void fmpz_mod_poly_eval_pow(
    fmpz_t eval,
    const fmpz_mod_poly_t P,
    fmpz_mod_poly_t alphapow,
    const fmpz_mod_ctx_t ctx);

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

FLINT_DLL int fmpz_mod_polyu1n_gcd_brown_smprime(
    fmpz_mod_polyun_t G,
    fmpz_mod_polyun_t Abar,
    fmpz_mod_polyun_t Bbar,
    fmpz_mod_polyun_t A,
    fmpz_mod_polyun_t B,
    const fmpz_mod_ctx_t ctx,
    fmpz_mod_poly_stack_t St_poly,
    fmpz_mod_polyun_stack_t St_polyun);

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

#ifdef __cplusplus
}
#endif

#endif

