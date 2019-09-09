/*
    Copyright (C) 2019 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#ifndef FMPZ_MOD_MPOLY_H
#define FMPZ_MOD_MPOLY_H

#ifdef FMPZ_MOD_MPOLY_INLINES_C
#define FMPZ_MOD_MPOLY_INLINE FLINT_DLL
#else
#define FMPZ_MOD_MPOLY_INLINE static __inline__
#endif

#undef ulong
#define ulong ulongxx /* interferes with system includes */
#include <stdio.h>
#undef ulong

#include <gmp.h>
#define ulong mp_limb_t

#include "flint.h"
#include "fmpz_mod.h"
#include "fmpz_mpoly.h"
#include "mpoly.h"

#ifdef __cplusplus
 extern "C" {
#endif

typedef struct
{
    fmpz_mod_ctx_t ffinfo;
    mpoly_ctx_t minfo;
} fmpz_mod_mpoly_ctx_struct;

typedef fmpz_mod_mpoly_ctx_struct fmpz_mod_mpoly_ctx_t[1];

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

FLINT_DLL void fmpz_mod_mpolyn_fit_length(fmpz_mod_mpolyn_t A,
                                 slong length, const fmpz_mod_mpoly_ctx_t ctx);

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

FLINT_DLL fmpz_mod_poly_struct * fmpz_mod_mpolyun_leadcoeff_ref(
                         fmpz_mod_mpolyun_t A, const fmpz_mod_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mod_mpolyn_content_last(fmpz_mod_poly_t a,
                    const fmpz_mod_mpolyn_t B, const fmpz_mod_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mod_mpolyun_content_last(fmpz_mod_poly_t a,
                   const fmpz_mod_mpolyun_t B, const fmpz_mod_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mod_mpolyn_divexact_last(fmpz_mod_mpolyn_t A,
                      const fmpz_mod_poly_t b, const fmpz_mod_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mod_mpolyun_divexact_last(fmpz_mod_mpolyun_t A,
                      const fmpz_mod_poly_t b, const fmpz_mod_mpoly_ctx_t ctx);

FLINT_DLL slong fmpz_mod_mpolyun_lastdeg(const fmpz_mod_mpolyun_t A,
                       const fmpz_mpoly_ctx_t ctx, const fmpz_mod_ctx_t fpctx);

FLINT_DLL void fmpz_mod_mpolyun_init(fmpz_mod_mpolyun_t A, flint_bitcnt_t bits,
                       const fmpz_mpoly_ctx_t ctx, const fmpz_mod_ctx_t fpctx);

FLINT_DLL void fmpz_mod_mpolyun_clear(fmpz_mod_mpolyun_t A,
                       const fmpz_mpoly_ctx_t ctx, const fmpz_mod_ctx_t fpctx);

FLINT_DLL fmpz * fmpz_mod_mpolyun_leadcoeff_last_ref(const fmpz_mod_mpolyun_t A,
                       const fmpz_mpoly_ctx_t ctx, const fmpz_mod_ctx_t fpctx);

FLINT_DLL void fmpz_mod_mpolyun_one(fmpz_mod_mpolyun_t A,
                       const fmpz_mpoly_ctx_t ctx, const fmpz_mod_ctx_t fpctx);

FLINT_DLL void fmpz_mod_mpolyun_mul_last(fmpz_mod_mpolyun_t A, fmpz_mod_poly_t b,
                       const fmpz_mpoly_ctx_t ctx, const fmpz_mod_ctx_t fpctx);

FLINT_DLL void fmpz_mod_mpolyun_set_modulus(fmpz_mod_mpolyun_t A,
                                                   const fmpz_mod_ctx_t fpctx);

FLINT_DLL void fmpz_mod_mpolyun_scalar_mul_fmpz_mod(fmpz_mod_mpolyun_t A,
       const fmpz_t c, const fmpz_mpoly_ctx_t ctx, const fmpz_mod_ctx_t fpctx);

FLINT_DLL int fmpz_mod_mpolyun_equal(
                    const fmpz_mod_mpolyun_t A, const fmpz_mod_mpolyun_t B,
                       const fmpz_mpoly_ctx_t ctx, const fmpz_mod_ctx_t fpctx);

FLINT_DLL int fmpz_mod_mpolyun_gcd_brown_bivar(
    fmpz_mod_mpolyun_t G, fmpz_mod_mpolyun_t Abar, fmpz_mod_mpolyun_t Bbar,
    fmpz_mod_mpolyun_t A, fmpz_mod_mpolyun_t B,
                       const fmpz_mpoly_ctx_t ctx, const fmpz_mod_ctx_t fpctx);

FLINT_DLL void fmpz_mod_mpolyun_print_pretty(const fmpz_mod_mpolyun_t poly,
      const char ** x, const fmpz_mpoly_ctx_t ctx, const fmpz_mod_ctx_t fpctx);

/* geobuckets ****************************************************************/
typedef struct fmpz_mpoly_geobucket
{
    fmpz_mpoly_struct polys[FLINT_BITS/2];
    slong length;
} fmpz_mpoly_geobucket_struct;

typedef fmpz_mpoly_geobucket_struct fmpz_mpoly_geobucket_t[1];

FLINT_DLL void fmpz_mpoly_geobucket_init(fmpz_mpoly_geobucket_t B,
                                                   const fmpz_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mpoly_geobucket_clear(fmpz_mpoly_geobucket_t B,
                                                   const fmpz_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mpoly_geobucket_empty(fmpz_mpoly_t p,
                         fmpz_mpoly_geobucket_t B, const fmpz_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mpoly_geobucket_print(fmpz_mpoly_geobucket_t B,
                                  const char ** x, const fmpz_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mpoly_geobucket_fit_length(fmpz_mpoly_geobucket_t B,
                                          slong i, const fmpz_mpoly_ctx_t ctx);

FLINT_DLL void _fmpz_mpoly_geobucket_fix(fmpz_mpoly_geobucket_t B, slong i,
                                                   const fmpz_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mpoly_geobucket_set(fmpz_mpoly_geobucket_t B,
                                   fmpz_mpoly_t p, const fmpz_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mpoly_geobucket_add(fmpz_mpoly_geobucket_t B,
                                   fmpz_mpoly_t p, const fmpz_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mpoly_geobucket_sub(fmpz_mpoly_geobucket_t B,
                                   fmpz_mpoly_t p, const fmpz_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mpoly_geobucket_set_fmpz(fmpz_mpoly_geobucket_t B,
                                         fmpz_t c, const fmpz_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mpoly_geobucket_gen(fmpz_mpoly_geobucket_t B, slong var,
                                                   const fmpz_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mpoly_geobucket_add_inplace(fmpz_mpoly_geobucket_t B1,
                        fmpz_mpoly_geobucket_t B2, const fmpz_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mpoly_geobucket_sub_inplace(fmpz_mpoly_geobucket_t B1,
                        fmpz_mpoly_geobucket_t B2, const fmpz_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mpoly_geobucket_neg_inplace(fmpz_mpoly_geobucket_t B1,
                                                   const fmpz_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mpoly_geobucket_mul_inplace(fmpz_mpoly_geobucket_t B1,
                        fmpz_mpoly_geobucket_t B2, const fmpz_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mpoly_geobucket_pow_ui_inplace(fmpz_mpoly_geobucket_t B1,
                                          ulong k, const fmpz_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mpoly_geobucket_pow_fmpz_inplace(fmpz_mpoly_geobucket_t B1,
                                   const fmpz_t k, const fmpz_mpoly_ctx_t ctx);

FLINT_DLL int fmpz_mpoly_geobucket_divides_inplace(fmpz_mpoly_geobucket_t B1,
                        fmpz_mpoly_geobucket_t B2, const fmpz_mpoly_ctx_t ctx);


/* Helpers for gcd_berlekamp_massey ******************************************/

typedef struct
{
    slong * degbounds;
    ulong * subdegs;
    fmpz_mod_discrete_log_pohlig_hellman_t dlogenv;
    nmod_discrete_log_pohlig_hellman_t dlogenv_sp;
} mpoly_bma_interpolate_ctx_struct;
typedef mpoly_bma_interpolate_ctx_struct mpoly_bma_interpolate_ctx_t[1];

FLINT_DLL void mpoly_bma_interpolate_ctx_init(mpoly_bma_interpolate_ctx_t Ictx,
                                                                  slong nvars);

FLINT_DLL void mpoly_bma_interpolate_ctx_clear(mpoly_bma_interpolate_ctx_t Ictx);

FLINT_DLL int nmod_mpoly_bma_get_fmpz_mpoly(fmpz_mpoly_t A,
     const fmpz_mpoly_ctx_t ctx, ulong alphashift, nmod_berlekamp_massey_t I,
              const mpoly_bma_interpolate_ctx_t Ictx, const nmodf_ctx_t fpctx);

/*
    nmod_mpoly "skeletons" - just the coefficients
*/
typedef struct
{
   mp_limb_t * coeffs;
   slong alloc;
   slong length;
} nmod_mpolyc_struct;

typedef nmod_mpolyc_struct nmod_mpolyc_t[1];

FLINT_DLL void nmod_mpolyc_init(nmod_mpolyc_t A);

FLINT_DLL void nmod_mpolyc_clear(nmod_mpolyc_t A);

FLINT_DLL void nmod_mpolyc_fit_length(nmod_mpolyc_t A, slong length);

typedef struct
{
   nmod_mpolyc_struct * coeffs;
   slong alloc;
   slong length;
} nmod_mpolycu_struct;

typedef nmod_mpolycu_struct nmod_mpolycu_t[1];

FLINT_DLL void nmod_mpolycu_init(nmod_mpolycu_t A);

FLINT_DLL void nmod_mpolycu_clear(nmod_mpolycu_t A);

FLINT_DLL void nmod_mpolycu_fit_length(nmod_mpolycu_t A, slong length);

/*
    fmpz_mod_mpoly "skeletons" - just the coefficients
*/

typedef struct
{
   fmpz * coeffs;
   slong alloc;
   slong length;
} fmpz_mpolyc_struct;

typedef fmpz_mpolyc_struct fmpz_mpolyc_t[1];

FLINT_DLL void fmpz_mpolyc_init(fmpz_mpolyc_t A);

FLINT_DLL void fmpz_mpolyc_clear(fmpz_mpolyc_t A);

FLINT_DLL void fmpz_mpolyc_fit_length(fmpz_mpolyc_t A, slong length);

typedef struct
{
   fmpz_mpolyc_struct * coeffs;
   slong alloc;
   slong length;
} fmpz_mpolycu_struct;

typedef fmpz_mpolycu_struct fmpz_mpolycu_t[1];

FLINT_DLL void fmpz_mpolycu_init(fmpz_mpolycu_t A);

FLINT_DLL void fmpz_mpolycu_clear(fmpz_mpolycu_t A);

FLINT_DLL void fmpz_mpolycu_fit_length(fmpz_mpolycu_t A, slong length);

FLINT_DLL void nmod_mpoly_bma_interpolate_alpha_powers(mp_limb_t * out,
  ulong w, const mpoly_bma_interpolate_ctx_t Ictx, const nmod_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mod_mpoly_bma_interpolate_alpha_powers(fmpz * out,
                     const fmpz_t w, const mpoly_bma_interpolate_ctx_t Ictx,
                       const fmpz_mpoly_ctx_t ctx, const fmpz_mod_ctx_t fpctx);

FLINT_DLL void fmpz_mpoly_set_skel(fmpz_mpolyc_t M, const fmpz_mpoly_t A,
   const fmpz * alpha, const fmpz_mpoly_ctx_t ctx, const fmpz_mod_ctx_t fpctx);

FLINT_DLL void fmpz_mpolyu_set_skel(fmpz_mpolycu_t M, const fmpz_mpolyu_t A,
   const fmpz * alpha, const fmpz_mpoly_ctx_t ctx, const fmpz_mod_ctx_t fpctx);

FLINT_DLL void fmpz_mpoly_copy_skel(fmpz_mpolyc_t M, const fmpz_mpolyc_t S);

FLINT_DLL void fmpz_mpolyu_copy_skel(fmpz_mpolycu_t M, const fmpz_mpolycu_t S);

FLINT_DLL void fmpz_mpoly_red_skel(fmpz_mpolyc_t Ared, const fmpz_mpoly_t A,
                                                   const fmpz_mod_ctx_t fpctx);

FLINT_DLL void fmpz_mpolyu_red_skel(fmpz_mpolycu_t Ared, const fmpz_mpolyu_t A,
                                                   const fmpz_mod_ctx_t fpctx);

FLINT_DLL void fmpz_mpoly_use_skel_mul(fmpz_t eval, fmpz_mpolyc_t Ared,
           fmpz_mpolyc_t M, const fmpz_mpolyc_t S, const fmpz_mod_ctx_t fpctx);

FLINT_DLL void fmpz_mpolyuu_use_skel_mul(fmpz_mod_mpolyun_t E,
                           const fmpz_mpolyu_t A, fmpz_mpolycu_t Ared,
                           fmpz_mpolycu_t M, const fmpz_mpolycu_t S,
                       const fmpz_mpoly_ctx_t ctx, const fmpz_mod_ctx_t fpctx);

FLINT_DLL int fmpz_mod_bma_get_fmpz_mpoly(fmpz_mpoly_t A,
          const fmpz_t alphashift, fmpz_mod_berlekamp_massey_t I,
          const mpoly_bma_interpolate_ctx_t Ictx, const fmpz_mpoly_ctx_t ctx,
                                                   const fmpz_mod_ctx_t fpctx);

typedef struct {
    nmod_berlekamp_massey_struct * coeffs;
    ulong * exps;
    slong length;
    slong alloc;
    slong pointcount;
} nmod_bma_mpoly_struct;

typedef nmod_bma_mpoly_struct nmod_bma_mpoly_t[1];

FLINT_DLL void nmod_bma_mpoly_init(nmod_bma_mpoly_t A);

FLINT_DLL void nmod_bma_mpoly_reset_prime(nmod_bma_mpoly_t A,
                                                      const nmodf_ctx_t fpctx);

FLINT_DLL void nmod_bma_mpoly_clear(nmod_bma_mpoly_t A);

FLINT_DLL void nmod_bma_mpoly_print(const nmod_bma_mpoly_t A);

FLINT_DLL void nmod_bma_mpoly_fit_length(nmod_bma_mpoly_t A, slong length,
                                                      const nmodf_ctx_t fpctx);

FLINT_DLL void nmod_bma_mpoly_zero(nmod_bma_mpoly_t L);

FLINT_DLL int nmod_bma_mpoly_reduce(nmod_bma_mpoly_t L);

FLINT_DLL void nmod_bma_mpoly_add_point(nmod_bma_mpoly_t L,
                              const nmod_mpolyun_t A, const nmodf_ctx_t fpctx);

FLINT_DLL int nmod_bma_mpoly_get_fmpz_mpolyu(fmpz_mpolyu_t A,
      const fmpz_mpoly_ctx_t ctx, ulong alphashift, const nmod_bma_mpoly_t L,
              const mpoly_bma_interpolate_ctx_t Ictx, const nmodf_ctx_t fpctx);

typedef struct {
    fmpz_mod_berlekamp_massey_struct * coeffs;
    ulong * exps;
    slong length;
    slong alloc;
    slong pointcount;
} fmpz_mod_bma_mpoly_struct;

typedef fmpz_mod_bma_mpoly_struct fmpz_mod_bma_mpoly_t[1];

FLINT_DLL void fmpz_mod_bma_mpoly_init(fmpz_mod_bma_mpoly_t A);

FLINT_DLL void fmpz_mod_bma_mpoly_reset_prime(
    fmpz_mod_bma_mpoly_t A,
    const fmpz_mod_ctx_t fpctx);

FLINT_DLL void fmpz_mod_bma_mpoly_clear(fmpz_mod_bma_mpoly_t A);

FLINT_DLL void fmpz_mod_bma_mpoly_print(
    fmpz_mod_bma_mpoly_t A,
    const mpoly_bma_interpolate_ctx_t Ictx);

FLINT_DLL void fmpz_mod_bma_mpoly_fit_length(
    fmpz_mod_bma_mpoly_t A,
    slong length,
    const fmpz_mod_ctx_t fpctx);

FLINT_DLL void fmpz_mod_bma_mpoly_zero(fmpz_mod_bma_mpoly_t L);

FLINT_DLL int fmpz_mod_bma_mpoly_reduce(fmpz_mod_bma_mpoly_t L);

FLINT_DLL void fmpz_mod_bma_mpoly_add_point(fmpz_mod_bma_mpoly_t L,
                    const fmpz_mod_mpolyun_t A, const fmpz_mpoly_ctx_t ctx,
                                                   const fmpz_mod_ctx_t fpctx);

FLINT_DLL int fmpz_mod_bma_mpoly_get_fmpz_mpolyu(fmpz_mpolyu_t A,
            const fmpz_t alphashift, const fmpz_mod_bma_mpoly_t L,
          const mpoly_bma_interpolate_ctx_t Ictx, const fmpz_mpoly_ctx_t ctx,
                                                   const fmpz_mod_ctx_t fpctx);

FLINT_DLL ulong fmpz_mod_mpolyn_bidegree(const fmpz_mod_mpolyn_t A);

FLINT_DLL ulong nmod_mpolyn_bidegree(const nmod_mpolyn_t A);

FLINT_DLL void fmpz_mpoly_eval_fmpz_mod(fmpz_t eval, const fmpz_mpoly_t A,
   const fmpz * alpha, const fmpz_mpoly_ctx_t ctx, const fmpz_mod_ctx_t fpctx);

FLINT_DLL void fmpz_mpolyu_symmetrize_coeffs(fmpz_mpolyu_t A,
                       const fmpz_mpoly_ctx_t ctx, const fmpz_mod_ctx_t fpctx);

FLINT_DLL void fmpz_mpolyuu_eval_fmpz_mod(fmpz_mod_mpolyn_t E,
                         const fmpz_mpoly_ctx_t ctx_mp, const fmpz_mpolyu_t A,
                             const fmpz * alpha, const fmpz_mpoly_ctx_t ctx,
                                                   const fmpz_mod_ctx_t fpctx);

FLINT_DLL mp_limb_t fmpz_mpoly_eval_nmod(const nmodf_ctx_t fpctx,
                            const fmpz_mpoly_t A, const mp_limb_t * alpha,
                                                   const fmpz_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mpolyuu_eval_nmod(nmod_mpolyn_t E,
                     const nmod_mpoly_ctx_t ctx_sp, const fmpz_mpolyu_t A,
                          const mp_limb_t * alpha, const fmpz_mpoly_ctx_t ctx);

typedef struct {
    slong mlength;
    slong malloc;
    mp_limb_t * coeffs;
    mp_limb_t * monomials;
    slong ealloc;
    mp_limb_t * evals;
} nmod_zip_struct;

typedef nmod_zip_struct nmod_zip_t[1];

FLINT_DLL void nmod_zip_init(nmod_zip_t Z);

FLINT_DLL void nmod_zip_clear(nmod_zip_t Z);

FLINT_DLL void nmod_zip_set_lengths(nmod_zip_t A, slong mlength, slong elength);

FLINT_DLL void nmod_zip_print(const nmod_zip_t Z, slong elength);

typedef struct {
    nmod_zip_struct * coeffs;
    ulong * exps;
    slong length;
    slong alloc;
    slong pointcount;
} nmod_zip_mpolyu_struct;

typedef nmod_zip_mpolyu_struct nmod_zip_mpolyu_t[1];

FLINT_DLL void nmod_zip_mpolyu_init(nmod_zip_mpolyu_t Z);

FLINT_DLL void nmod_zip_mpolyu_clear(nmod_zip_mpolyu_t Z);

FLINT_DLL void nmod_zip_mpolyu_fit_length(nmod_zip_mpolyu_t A, slong length);

FLINT_DLL void nmod_zip_mpolyu_fit_poly(nmod_zip_mpolyu_t Z, fmpz_mpolyu_t H,
                                                            slong eval_length);

FLINT_DLL void nmod_mpoly_set_skel(nmod_mpolyc_t S,
                     const nmod_mpoly_ctx_t ctx_sp, const fmpz_mpoly_t A,
                          const mp_limb_t * alpha, const fmpz_mpoly_ctx_t ctx);

FLINT_DLL void nmod_zip_mpolyu_set_skel(nmod_zip_mpolyu_t Z,
                       const nmod_mpoly_ctx_t ctx_sp, const fmpz_mpolyu_t A,
                          const mp_limb_t * alpha, const fmpz_mpoly_ctx_t ctx);

FLINT_DLL void nmod_zip_mpolyuu_print(const nmod_zip_mpolyu_t A);

FLINT_DLL int nmod_zip_mpolyuu_add_point(nmod_zip_mpolyu_t L,
                                                       const nmod_mpolyun_t A);

FLINT_DLL void nmod_mpolyu_set_skel(nmod_mpolycu_t S,
                        const nmod_mpoly_ctx_t ctx_sp, const fmpz_mpolyu_t A,
                          const mp_limb_t * alpha, const fmpz_mpoly_ctx_t ctx);

FLINT_DLL void nmod_mpoly_red_skel(nmod_mpolyc_t Ared, const fmpz_mpoly_t A,
                                                      const nmodf_ctx_t fpctx);

FLINT_DLL void nmod_mpolyu_red_skel(nmod_mpolycu_t Ared, const fmpz_mpolyu_t A,
                                                      const nmodf_ctx_t fpctx);

FLINT_DLL void nmod_mpoly_copy_skel(nmod_mpolyc_t M, const nmod_mpolyc_t S);

FLINT_DLL void nmod_mpolyu_copy_skel(nmod_mpolycu_t M, const nmod_mpolycu_t S);

FLINT_DLL void nmod_mpoly_pow_skel(nmod_mpolyc_t M, const nmod_mpolyc_t S,
                                          ulong k, const nmod_mpoly_ctx_t ctx);

FLINT_DLL void nmod_mpolyu_pow_skel(nmod_mpolycu_t M, const nmod_mpolycu_t S,
                                          ulong k, const nmod_mpoly_ctx_t ctx);

FLINT_DLL mp_limb_t nmod_mpoly_use_skel_mul(const nmod_mpolyc_t Ared,
  nmod_mpolyc_t Acur, const nmod_mpolyc_t Ainc, const nmod_mpoly_ctx_t ctx_sp);

FLINT_DLL void nmod_mpolyuu_use_skel_mul(nmod_mpolyun_t E,
     const fmpz_mpolyu_t A, const nmod_mpolycu_t Ared, nmod_mpolycu_t Acur,
                     const nmod_mpolycu_t Ainc, const nmod_mpoly_ctx_t ctx_sp);

typedef enum {
    nmod_zip_find_coeffs_good,
    nmod_zip_find_coeffs_no_match,
    nmod_zip_find_coeffs_non_invertible
} nmod_zip_find_coeffs_ret_t;

FLINT_DLL nmod_zip_find_coeffs_ret_t nmod_zip_find_coeffs(nmod_zip_t Z,
               nmod_poly_t master, slong pointcount, const nmodf_ctx_t ffinfo);

FLINT_DLL nmod_zip_find_coeffs_ret_t nmod_mpolyu_zip_find_coeffs(
                           nmod_zip_mpolyu_t Z, const nmod_mpoly_ctx_t ctx_sp);

FLINT_DLL int fmpz_mpolyu_addinterp_zip(fmpz_mpolyu_t H, const fmpz_t Hmodulus,
                          const nmod_zip_mpolyu_t Z, const nmodf_ctx_t ffinfo);

FLINT_DLL int fmpz_mpoly_repack_bits_inplace(fmpz_mpoly_t A, flint_bitcnt_t Abits,
                                                   const fmpz_mpoly_ctx_t ctx);

FLINT_DLL int fmpz_mpolyu_repack_bits(fmpz_mpolyu_t A, flint_bitcnt_t Abits,
                                                   const fmpz_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mpoly_set_bits(fmpz_mpoly_t A, flint_bitcnt_t Abits,
                                                   const fmpz_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mpolyu_set_bits(fmpz_mpolyu_t A, flint_bitcnt_t Abits,
                                                   const fmpz_mpoly_ctx_t ctx);

FLINT_DLL int fmpz_mpolyuu_eval_all_but_one_nmod(slong * out_deg,
                        nmod_poly_t out, const fmpz_mpolyu_t A, slong var,
                               mp_limb_t * values, const fmpz_mpoly_ctx_t ctx);

FLINT_DLL slong fmpz_mpolyuu_gcd_degree_bound_minor(slong * Adeg, slong * Bdeg,
                      const fmpz_mpolyu_t A, const fmpz_mpolyu_t B, slong var,
                               const fmpz_mpoly_ctx_t ctx, flint_rand_t state);

FLINT_DLL void fmpz_mpoly_ksub_content(fmpz_t content, const fmpz_mpoly_t A,
                            const ulong * subdegs, const fmpz_mpoly_ctx_t ctx);

/* Helpers for array methods *************************************************/

FLINT_DLL void _fmpz_mpoly_mul_array_chunked_DEG(fmpz_mpoly_t P,
                             const fmpz_mpoly_t A, const fmpz_mpoly_t B, 
                                       ulong degb, const fmpz_mpoly_ctx_t ctx);

FLINT_DLL void _fmpz_mpoly_mul_array_chunked_LEX(fmpz_mpoly_t P,
                             const fmpz_mpoly_t A, const fmpz_mpoly_t B, 
                              const ulong * mults, const fmpz_mpoly_ctx_t ctx);

FLINT_DLL void _fmpz_mpoly_addmul_array1_slong1(ulong * poly1, 
                 const slong * poly2, const ulong * exp2, slong len2,
                          const slong * poly3, const ulong * exp3, slong len3);

FLINT_DLL void _fmpz_mpoly_addmul_array1_slong(ulong * poly1, 
                 const slong * poly2, const ulong * exp2, slong len2,
                          const slong * poly3, const ulong * exp3, slong len3);

FLINT_DLL void _fmpz_mpoly_addmul_array1_slong2(ulong * poly1, 
                 const slong * poly2, const ulong * exp2, slong len2,
                          const slong * poly3, const ulong * exp3, slong len3);

FLINT_DLL void _fmpz_mpoly_addmul_array1_fmpz(fmpz * poly1, 
                 const fmpz * poly2, const ulong * exp2, slong len2,
                           const fmpz * poly3, const ulong * exp3, slong len3);

FLINT_DLL void _fmpz_mpoly_submul_array1_slong(ulong * poly1, 
                  const slong * poly2, const ulong * exp2, slong len2,
                          const slong * poly3, const ulong * exp3, slong len3);

FLINT_DLL void _fmpz_mpoly_submul_array1_slong2(ulong * poly1, 
                  const slong * poly2, const ulong * exp2, slong len2,
                          const slong * poly3, const ulong * exp3, slong len3);

FLINT_DLL void _fmpz_mpoly_submul_array1_slong1(ulong * poly1, 
                 const slong * poly2, const ulong * exp2, slong len2,
                          const slong * poly3, const ulong * exp3, slong len3);

FLINT_DLL void _fmpz_mpoly_submul_array1_fmpz(fmpz * poly1, 
                 const fmpz * poly2, const ulong * exp2, slong len2,
                           const fmpz * poly3, const ulong * exp3, slong len3);

FLINT_DLL void _fmpz_mpoly_submul_array1_slong_1(ulong * poly1, 
                          slong d, const ulong exp2,
                          const slong * poly3, const ulong * exp3, slong len3);

FLINT_DLL void _fmpz_mpoly_submul_array1_slong2_1(ulong * poly1, 
                           slong d, const ulong exp2,
                          const slong * poly3, const ulong * exp3, slong len3);

FLINT_DLL void _fmpz_mpoly_submul_array1_fmpz_1(fmpz * poly1, 
                          const fmpz_t d, ulong exp2,
                           const fmpz * poly3, const ulong * exp3, slong len3);

FLINT_DLL void mpoly_main_variable_split_LEX(slong * ind, ulong * pexp,
             const ulong * Aexp,
             slong l1, slong Alen, const ulong * mults, slong num, slong Abits);

FLINT_DLL void mpoly_main_variable_split_DEG(slong * ind, ulong * pexp,
             const ulong * Aexp,
             slong l1, slong Alen, ulong deg, slong num, slong Abits);

FLINT_DLL slong fmpz_mpoly_append_array_sm1_LEX(fmpz_mpoly_t P,
                        slong Plen, ulong * coeff_array,
                  const ulong * mults, slong num, slong array_size, slong top);
FLINT_DLL slong fmpz_mpoly_append_array_sm2_LEX(fmpz_mpoly_t P,
                        slong Plen, ulong * coeff_array,
                  const ulong * mults, slong num, slong array_size, slong top);
FLINT_DLL slong fmpz_mpoly_append_array_sm3_LEX(fmpz_mpoly_t P,
                         slong Plen, ulong * coeff_array,
                  const ulong * mults, slong num, slong array_size, slong top);
FLINT_DLL slong fmpz_mpoly_append_array_fmpz_LEX(fmpz_mpoly_t P,
                        slong Plen, fmpz * coeff_array,
                  const ulong * mults, slong num, slong array_size, slong top);

FLINT_DLL slong fmpz_mpoly_append_array_sm1_DEGLEX(fmpz_mpoly_t P,
          slong Plen, ulong * coeff_array, slong top, slong nvars, slong degb);
FLINT_DLL slong fmpz_mpoly_append_array_sm2_DEGLEX(fmpz_mpoly_t P,
          slong Plen, ulong * coeff_array, slong top, slong nvars, slong degb);
FLINT_DLL slong fmpz_mpoly_append_array_sm3_DEGLEX(fmpz_mpoly_t P,
          slong Plen, ulong * coeff_array, slong top, slong nvars, slong degb);
FLINT_DLL slong fmpz_mpoly_append_array_fmpz_DEGLEX(fmpz_mpoly_t P,
           slong Plen, fmpz * coeff_array, slong top, slong nvars, slong degb);

FLINT_DLL slong fmpz_mpoly_append_array_sm1_DEGREVLEX(fmpz_mpoly_t P,
          slong Plen, ulong * coeff_array, slong top, slong nvars, slong degb);
FLINT_DLL slong fmpz_mpoly_append_array_sm2_DEGREVLEX(fmpz_mpoly_t P,
          slong Plen, ulong * coeff_array, slong top, slong nvars, slong degb);
FLINT_DLL slong fmpz_mpoly_append_array_sm3_DEGREVLEX(fmpz_mpoly_t P,
          slong Plen, ulong * coeff_array, slong top, slong nvars, slong degb);
FLINT_DLL slong fmpz_mpoly_append_array_fmpz_DEGREVLEX(fmpz_mpoly_t P,
           slong Plen, fmpz * coeff_array, slong top, slong nvars, slong degb);

FLINT_DLL slong _fmpz_mpoly_from_ulong_array(fmpz ** poly1,
                         ulong ** exp1, slong * alloc, ulong * poly2,
                          const slong * mults, slong num, slong bits, slong k);

FLINT_DLL slong _fmpz_mpoly_from_ulong_array2(fmpz ** poly1,
                         ulong ** exp1, slong * alloc, ulong * poly2, 
                          const slong * mults, slong num, slong bits, slong k);

FLINT_DLL slong _fmpz_mpoly_from_ulong_array1(fmpz ** poly1,
                         ulong ** exp1, slong * alloc, ulong * poly2,
                          const slong * mults, slong num, slong bits, slong k);

FLINT_DLL slong _fmpz_mpoly_from_fmpz_array(fmpz ** poly1,
                         ulong ** exp1, slong * alloc, fmpz * poly2,
                          const slong * mults, slong num, slong bits, slong k);

FLINT_DLL void _fmpz_mpoly_to_ulong_array2(ulong * p, const fmpz * coeffs,
                                                const ulong * exps, slong len);

FLINT_DLL void _fmpz_mpoly_to_ulong_array1(ulong * p, const fmpz * coeffs,
                                                const ulong * exps, slong len);

FLINT_DLL void _fmpz_mpoly_to_ulong_array(ulong * p, const fmpz * coeffs,
                                                const ulong * exps, slong len);

FLINT_DLL void _fmpz_mpoly_to_fmpz_array(fmpz * p, const fmpz * coeffs,
                                                const ulong * exps, slong len);


/* Misc arithmetic - has nothing to do with mpoly, should be moved out *******/

FMPZ_MPOLY_INLINE
void _fmpz_mpoly_sub_uiuiui_fmpz(ulong * c, const fmpz_t d)
{
   fmpz fc = *d;

   if (!COEFF_IS_MPZ(fc))
   {
        ulong f0, f1, f2;
        f0 = fc;
        f1 = f2 = FLINT_SIGN_EXT(f0);
        sub_dddmmmsss(c[2], c[1], c[0], c[2], c[1], c[0], f2, f1, f0);
   } else
   {
      slong size = fmpz_size(d);
      __mpz_struct * m = COEFF_TO_PTR(fc);
      if (fmpz_sgn(d) < 0)
         mpn_add(c, c, 3, m->_mp_d, size);
      else
         mpn_sub(c, c, 3, m->_mp_d, size);
   }
}

FMPZ_MPOLY_INLINE
void _fmpz_mpoly_add_uiuiui_fmpz(ulong * c, const fmpz_t d)
{
    fmpz fc = *d;

    if (!COEFF_IS_MPZ(fc))
    {
        ulong f0, f1, f2;
        f0 = fc;
        f1 = f2 = FLINT_SIGN_EXT(f0);
        add_sssaaaaaa(c[2], c[1], c[0], c[2], c[1], c[0], f2, f1, f0);
   } else
   {
      slong size = fmpz_size(d);
      __mpz_struct * m = COEFF_TO_PTR(fc);
      if (fmpz_sgn(d) < 0)
         mpn_sub(c, c, 3, m->_mp_d, size);
      else
         mpn_add(c, c, 3, m->_mp_d, size);
   }
}

FMPZ_MPOLY_INLINE
void _fmpz_mpoly_submul_uiuiui_fmpz(ulong * c, slong d1, slong d2)
{
    ulong p[2], p2;
    smul_ppmm(p[1], p[0], d1, d2);
    p2 = FLINT_SIGN_EXT(p[1]);
    sub_dddmmmsss(c[2], c[1], c[0], c[2], c[1], c[0], p2, p[1], p[0]);
}

FMPZ_MPOLY_INLINE
void _fmpz_mpoly_addmul_uiuiui_fmpz(ulong * c, slong d1, slong d2)
{
    ulong p[2], p2;
    smul_ppmm(p[1], p[0], d1, d2);
    p2 = FLINT_SIGN_EXT(p[1]);
    add_sssaaaaaa(c[2], c[1], c[0], c[2], c[1], c[0], p2, p[1], p[0]);
}



/******************************************************************************

   Internal consistency checks

******************************************************************************/

/*
   test that r is a valid remainder upon division by g
   this means that if c*x^a is a term of r and x^a is divisible by the leading
   monomial of g, then |c| < |leading coefficient of g|
*/
FMPZ_MPOLY_INLINE
void fmpz_mpoly_remainder_test(const fmpz_mpoly_t r, const fmpz_mpoly_t g,
                                                    const fmpz_mpoly_ctx_t ctx)
{
   slong i, N, bits;
   ulong mask = 0;
   ulong * rexp, * gexp;

   bits = FLINT_MAX(r->bits, g->bits);
   N = mpoly_words_per_exp(bits, ctx->minfo);

   if (g->length == 0 )
      flint_throw(FLINT_ERROR, "Zero denominator in remainder test");

   if (r->length == 0 )
      return;

   rexp = (ulong *) flint_malloc(N*r->length*sizeof(ulong));
   gexp = (ulong *) flint_malloc(N*1        *sizeof(ulong));
   mpoly_repack_monomials(rexp, bits, r->exps, r->bits, r->length, ctx->minfo);
   mpoly_repack_monomials(gexp, bits, g->exps, g->bits, 1,         ctx->minfo);

   /* mask with high bit set in each field of exponent vector */
   for (i = 0; i < FLINT_BITS/bits; i++)
      mask = (mask << bits) + (UWORD(1) << (bits - 1));

    for (i = 0; i < r->length; i++)
    {
        int divides;

        if (bits <= FLINT_BITS)
            divides = mpoly_monomial_divides_test(rexp + i*N, gexp + 0*N, N, mask);
        else
            divides = mpoly_monomial_divides_mp_test(rexp + i*N, gexp + 0*N, N, bits);

        if (divides && fmpz_cmpabs(g->coeffs + 0, r->coeffs + i) <= 0)
        {
            flint_printf("fmpz_mpoly_remainder_test FAILED i = %wd\n", i);
            flint_printf("rem ");fmpz_mpoly_print_pretty(r, NULL, ctx); printf("\n\n");
            flint_printf("den ");fmpz_mpoly_print_pretty(g, NULL, ctx); printf("\n\n");
            flint_abort();
        }
    }

   flint_free(rexp);
   flint_free(gexp);
}


/*
   test that r is a valid remainder upon division by g over Q
   this means that no term of r is divisible by lt(g)
*/
FMPZ_MPOLY_INLINE
void fmpz_mpoly_remainder_strongtest(const fmpz_mpoly_t r, const fmpz_mpoly_t g,
                                                    const fmpz_mpoly_ctx_t ctx)
{
    slong i, N, bits;
    ulong mask = 0;
    ulong * rexp, * gexp;

    bits = FLINT_MAX(r->bits, g->bits);
    N = mpoly_words_per_exp(bits, ctx->minfo);

    if (g->length == 0 )
        flint_throw(FLINT_ERROR, "Zero denominator in remainder test");

    if (r->length == 0 )
        return;

    rexp = (ulong *) flint_malloc(N*r->length*sizeof(ulong));
    gexp = (ulong *) flint_malloc(N*1        *sizeof(ulong));
    mpoly_repack_monomials(rexp, bits, r->exps, r->bits, r->length, ctx->minfo);
    mpoly_repack_monomials(gexp, bits, g->exps, g->bits, 1,         ctx->minfo);

    /* mask with high bit set in each field of exponent vector */
    for (i = 0; i < FLINT_BITS/bits; i++)
        mask = (mask << bits) + (UWORD(1) << (bits - 1));

    for (i = 0; i < r->length; i++)
    {
        int divides;

        if (bits <= FLINT_BITS)
            divides = mpoly_monomial_divides_test(rexp + i*N, gexp + 0*N, N, mask);
        else
            divides = mpoly_monomial_divides_mp_test(rexp + i*N, gexp + 0*N, N, bits);

        if (divides)
        {
            flint_printf("fmpz_mpoly_remainder_strongtest FAILED i = %wd\n", i);
            flint_printf("rem ");fmpz_mpoly_print_pretty(r, NULL, ctx); printf("\n\n");
            flint_printf("den ");fmpz_mpoly_print_pretty(g, NULL, ctx); printf("\n\n");
            flint_abort();
        }
    }

    flint_free(rexp);
    flint_free(gexp);
}

#ifdef __cplusplus
}
#endif

#endif
