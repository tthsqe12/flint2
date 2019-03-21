

#ifndef NMOD_POLYDR_H
#define NMOD_POLYDR_H

#ifdef NMOD_POLYDR_INLINES_C
#define NMOD_POLYDR_INLINE FLINT_DLL
#else
#define NMOD_POLYDR_INLINE static __inline__
#endif

#undef ulong
#define ulong ulongxx /* interferes with system includes */
#include <stdio.h>
#undef ulong
#include <gmp.h>
#define ulong mp_limb_t

#include "flint.h"
#include "nmod_vec.h"
#include "nmod_mat.h"
#include "ulong_extras.h"
#include "fmpz.h"

#ifdef __cplusplus
    extern "C" {
#endif

/* nmod_polydr_t - nmod_poly_t done right */




typedef struct
{
    mp_ptr coeffs;
    slong alloc;
    slong length;
} nmod_polydr_struct;

typedef nmod_polydr_struct nmod_polydr_t[1];


FLINT_DLL void nmod_polydr_init(nmod_polydr_t poly, const nmod_ctx_t ctx);

FLINT_DLL void nmod_polydr_init2(nmod_polydr_t poly, slong alloc, const nmod_ctx_t ctx);

FLINT_DLL void nmod_polydr_realloc(nmod_polydr_t poly, slong alloc, const nmod_ctx_t ctx);

FLINT_DLL void nmod_polydr_clear(nmod_polydr_t poly, const nmod_ctx_t ctx);

FLINT_DLL void nmod_polydr_fit_length(nmod_polydr_t poly, slong alloc, const nmod_ctx_t ctx);

NMOD_POLYDR_INLINE
void _nmod_polydr_set_length(nmod_polydr_t poly, slong len)
{
    poly->length = len;
}

NMOD_POLYDR_INLINE
void _nmod_polydr_normalise(nmod_polydr_t poly)
{
    while (poly->length > 0 && (poly->coeffs[poly->length - 1] == WORD(0)))
        poly->length--;
}

NMOD_POLYDR_INLINE
slong nmod_polydr_length(const nmod_polydr_t poly, const nmod_ctx_t ctx)
{
    return poly->length;
}

NMOD_POLYDR_INLINE
slong nmod_polydr_degree(const nmod_polydr_t poly, const nmod_ctx_t ctx)
{
    return poly->length - 1;
}

NMOD_POLYDR_INLINE
mp_bitcnt_t nmod_polydr_max_bits(const nmod_polydr_t poly)
{
    return _nmod_vec_max_bits(poly->coeffs, poly->length);
}

NMOD_POLYDR_INLINE
mp_ptr nmod_polydr_lead(const nmod_polydr_t poly)
{
    if (poly->length > 0)
        return poly->coeffs + (poly->length - 1);
    else
        return NULL;
}

NMOD_POLYDR_INLINE
void nmod_polydr_set(nmod_polydr_t a, const nmod_polydr_t b, const nmod_ctx_t ctx)
{
    if (a != b)
    {
        nmod_polydr_fit_length(a, b->length, ctx);
        flint_mpn_copyi(a->coeffs, b->coeffs, b->length);
        a->length = b->length;
    }
}


NMOD_POLYDR_INLINE
void nmod_polydr_set_nmod_poly(nmod_polydr_t a, const nmod_poly_t b, const nmod_ctx_t ctx)
{
    nmod_polydr_fit_length(a, b->length, ctx);
    flint_mpn_copyi(a->coeffs, b->coeffs, b->length);
    a->length = b->length;
}

NMOD_POLYDR_INLINE
void nmod_polydr_swap(nmod_polydr_t poly1, nmod_polydr_t poly2, const nmod_ctx_t ctx)
{
    slong t;
    mp_ptr tp;

    t = poly1->alloc;
    poly1->alloc = poly2->alloc;
    poly2->alloc = t;

    t = poly1->length;
    poly1->length = poly2->length;
    poly2->length = t;

    tp = poly1->coeffs;
    poly1->coeffs = poly2->coeffs;
    poly2->coeffs = tp;
}

NMOD_POLYDR_INLINE
void nmod_polydr_zero(nmod_polydr_t res, const nmod_ctx_t ctx)
{
    res->length = 0;
}

NMOD_POLYDR_INLINE
void nmod_polydr_one(nmod_polydr_t res, const nmod_ctx_t ctx)
{
    nmod_polydr_fit_length(res, 1, ctx);
    res->length = (ctx->mod.n != UWORD(1));
    res->coeffs[0] = UWORD(1);
}

NMOD_POLYDR_INLINE
int nmod_polydr_is_zero(const nmod_polydr_t poly, const nmod_ctx_t ctx)
{
    return (poly->length == 0);
}

NMOD_POLYDR_INLINE
int nmod_polydr_is_one(const nmod_polydr_t poly, const nmod_ctx_t ctx)
{
    return (poly->length == 1) && (poly->coeffs[0] == 1);
}

FLINT_DLL void nmod_polydr_set_coeff_ui(nmod_polydr_t poly, slong j, ulong c, const nmod_ctx_t ctx);


FLINT_DLL void nmod_polydr_make_monic(nmod_polydr_t output,
                              const nmod_polydr_t input, const nmod_ctx_t ctx);

FLINT_DLL void nmod_polydr_add(nmod_polydr_t A,
           const nmod_polydr_t B, const nmod_polydr_t C, const nmod_ctx_t ctx);

FLINT_DLL void nmod_polydr_sub(nmod_polydr_t A,
           const nmod_polydr_t B, const nmod_polydr_t C, const nmod_ctx_t ctx);

FLINT_DLL void nmod_polydr_mul(nmod_polydr_t A,
           const nmod_polydr_t B, const nmod_polydr_t C, const nmod_ctx_t ctx);

FLINT_DLL void nmod_polydr_div(nmod_polydr_t Q,
           const nmod_polydr_t A, const nmod_polydr_t B, const nmod_ctx_t ctx);

FLINT_DLL void nmod_polydr_rem(nmod_polydr_t R, 
           const nmod_polydr_t A, const nmod_polydr_t B, const nmod_ctx_t ctx);

FLINT_DLL void nmod_polydr_divrem(nmod_polydr_t Q, nmod_polydr_t R,
           const nmod_polydr_t A, const nmod_polydr_t B, const nmod_ctx_t ctx);

FLINT_DLL void nmod_polydr_gcd(nmod_polydr_t G,
           const nmod_polydr_t A, const nmod_polydr_t B, const nmod_ctx_t ctx);

#endif
