/*
    Copyright (C) 2019 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#ifndef UI_FACTOR_H
#define UI_FACTOR_H

#ifdef UI_FACTOR_INLINES_C
#define UI_FACTOR_INLINE FLINT_DLL
#else
#define UI_FACTOR_INLINE static __inline__
#endif

#undef ulong
#define ulong ulongxx /* interferes with system includes */
#include <stdio.h>
#undef ulong
#include <gmp.h>
#define ulong mp_limb_t

#include "flint/flint.h"
#include "flint/fmpz.h"


#define likely(x)       __builtin_expect((x),1)
#define unlikely(x)     __builtin_expect((x),0)

/*
#define likely(x)       (x)
#define unlikely(x)     (x)
*/

#ifdef __cplusplus
 extern "C" {
#endif

/*********** representation of an integer as prod_i base[i]^pow[i] ***********/
typedef struct {
    ulong base, pow;
} ui_factor_entry;

typedef struct {
    ui_factor_entry * data;
    slong alloc;
    slong length;
} ui_factor_struct;

typedef ui_factor_struct ui_factor_t[1];


/******************** vectors of mpz's and ui_factor's ***********************/
typedef struct {
    __mpz_struct ** mpz_array;
    slong mpz_alloc;
    slong mpz_top;
    ui_factor_struct ** factor_array;
    slong factor_alloc;
    slong factor_top;
} ui_factor_stack_struct;

typedef ui_factor_stack_struct ui_factor_stack_t[1];


/******** sieve type that will hopefully eventually be done right ************/
/* the index of this struct represents the odd number n = 2*index + 1 */
/* index = 0 (n = 1) is special */
typedef struct {
  unsigned int pminus;      /* smallest prime factor p of n */
  unsigned int cofactor;    /* index (n/p - 1)/2 of the cofactor n/p */
} ui_factor_sieve_entry;

typedef struct {
    ui_factor_sieve_entry * array;
    ulong max;     /* the biggest number we can factor */
    slong alloc;
} ui_factor_sieve_struct;

typedef ui_factor_sieve_struct ui_factor_sieve_t[1];


/* ui_factor_sieve_t *********************************************************/

FLINT_DLL void ui_factor_sieve_init(ui_factor_sieve_t S);

FLINT_DLL void ui_factor_sieve_clear(ui_factor_sieve_t S);

FLINT_DLL void ui_factor_sieve_fit_length(ui_factor_sieve_t S, slong length);

FLINT_DLL void ui_factor_sieve_build(ui_factor_sieve_t S, ulong m);


/* ui_factor_stack_t *********************************************************/

FLINT_DLL void ui_factor_stack_init(ui_factor_stack_t S);

FLINT_DLL void ui_factor_stack_clear(ui_factor_stack_t S);

FLINT_DLL mpz_ptr * ui_factor_stack_fit_request_mpz(ui_factor_stack_t S, slong k);

FLINT_DLL ui_factor_struct ** ui_factor_stack_fit_request_factor(
                                                 ui_factor_stack_t S, slong k);

UI_FACTOR_INLINE mpz_ptr ui_factor_stack_take_top_mpz(ui_factor_stack_t S)
{
    /* assume the request for 1 has already been fitted */
    mpz_ptr * mpz_top = S->mpz_array + S->mpz_top;
    FLINT_ASSERT(S->mpz_top + 1 <= S->mpz_alloc);
    S->mpz_top += 1;
    return mpz_top[0];
}

UI_FACTOR_INLINE ui_factor_struct * ui_factor_stack_take_top_factor(ui_factor_stack_t S)
{
    /* assume the request for 1 has already been fitted */
    ui_factor_struct ** factor_top = S->factor_array + S->factor_top;
    FLINT_ASSERT(S->factor_top + 1 <= S->factor_alloc);
    S->factor_top += 1;
    return factor_top[0];
}

UI_FACTOR_INLINE void ui_factor_stack_give_back_mpz(ui_factor_stack_t S, slong k)
{
    FLINT_ASSERT(S->mpz_top >= k);
    S->mpz_top -= k;
}

UI_FACTOR_INLINE void ui_factor_stack_give_back_factor(ui_factor_stack_t S, slong k)
{
    FLINT_ASSERT(S->factor_top >= k);
    S->factor_top -= k;
}



/* ui_factor_t ***************************************************************/

UI_FACTOR_INLINE void ui_factor_init(ui_factor_t f)
{
    f->data = NULL;
    f->alloc = 0;
    f->length = 0;
}

FLINT_DLL void ui_factor_init2(ui_factor_t f, slong length);

UI_FACTOR_INLINE void ui_factor_clear(ui_factor_t f)
{
    if (f->data)
        flint_free(f->data);
}

UI_FACTOR_INLINE void ui_factor_swap(ui_factor_t f, ui_factor_t g)
{
    ui_factor_t t;
    *t = *f;
    *f = *g;
    *g = *t;
}

FLINT_DLL void ui_factor_fit_length(ui_factor_t f, ulong length);

FLINT_DLL void ui_factor_print(const ui_factor_t f);

FLINT_DLL int ui_factor_is_canonical(const ui_factor_t f);

UI_FACTOR_INLINE void ui_factor_one(ui_factor_t f)
{
    f->length = 0;
}

FLINT_DLL void ui_factor_push_factor(ui_factor_t f, ulong base, ulong pow);

FLINT_DLL void ui_factor_push_ui_without_sieve(ui_factor_t f, ulong b);

FLINT_DLL void ui_factor_push_ui_with_sieve(ui_factor_t f, ulong b,
                                                    const ui_factor_sieve_t S);

FLINT_DLL ulong ui_factor_get_mpz_2exp(mpz_t x, const ui_factor_t f,
                                                         ui_factor_stack_t St);

FLINT_DLL int ui_factor_equal_fmpz(const ui_factor_t f, const fmpz_t x);

FLINT_DLL int ui_factor_equal_mpz(const ui_factor_t f, const mpz_t x);

FLINT_DLL void ui_factor_remove_gcd(ui_factor_t f, ui_factor_t g);

FLINT_DLL void ui_factor_mul(ui_factor_t z, const ui_factor_t f,
                                                          const ui_factor_t g);

FLINT_DLL void ui_factor_pow_inplace(ui_factor_t f, ulong pow);

FLINT_DLL void ui_factor_mulpow_inplace(ui_factor_t f, const ui_factor_t g,
                                                                      ulong p);

FLINT_DLL void ui_factor_canonicalise(ui_factor_t f);

/****** useful mpn stuff *****************************************************/

FLINT_DLL mp_size_t flint_mpn_mul_11(mp_limb_t * y,
                 const mp_limb_t * x, mp_size_t n, mp_limb_t a1, mp_limb_t a2);

FLINT_DLL mp_size_t flint_mpn_mul_111(mp_limb_t * y, const mp_limb_t * x,
                        mp_size_t n, mp_limb_t a1, mp_limb_t a2, mp_limb_t a3);

UI_FACTOR_INLINE mp_limb_t * flint_mpz_fit_length(mpz_ptr z, mp_size_t n)
{
    if (n > z->_mp_alloc)
    {
        n = FLINT_MAX(z->_mp_alloc + z->_mp_alloc/2, n);
        return (mp_ptr) _mpz_realloc(z, n);
    }
    else
    {
        return z->_mp_d;
    }
}


/* zassenhaus ****************************************************************/

FLINT_DLL void zassenhaus_subset_first(slong * s, slong r, slong m);

FLINT_DLL int zassenhaus_subset_next(slong * s, slong r);

FLINT_DLL slong zassenhaus_subset_next_disjoint(slong * s, slong r);

typedef struct {
    slong deg;
    unsigned char * pos_degs;   /* possible degrees: entries are 0 or 1*/
    slong new_length;
    slong new_total;
    slong * new_degs;
    slong alloc;
} zassenhaus_prune_struct;

typedef zassenhaus_prune_struct zassenhaus_prune_t[1];

UI_FACTOR_INLINE
void zassenhaus_prune_init(zassenhaus_prune_t Z)
{
    Z->deg = 0;
    Z->pos_degs = NULL;
    Z->new_length = 0;
    Z->new_total = 0;
    Z->new_degs = NULL;
    Z->alloc = 0;
}

FLINT_DLL void zassenhaus_prune_clear(zassenhaus_prune_t Z);

FLINT_DLL void zassenhaus_prune_set_degree(zassenhaus_prune_t Z, slong d);

UI_FACTOR_INLINE
void zassenhaus_prune_start_add_factors(zassenhaus_prune_t Z)
{
    Z->new_length = 0;
    Z->new_total = 0;
}

FLINT_DLL void zassenhaus_prune_add_factor(zassenhaus_prune_t Z,
                                                         slong deg, slong exp);

FLINT_DLL void zassenhaus_prune_end_add_factors(zassenhaus_prune_t Z);

FLINT_DLL int zassenhaus_prune_must_be_irreducible(const zassenhaus_prune_t Z);

UI_FACTOR_INLINE
int zassenhaus_prune_degree_is_possible(const zassenhaus_prune_t Z, slong d)
{
    if (d <= 0)
        return d == 0;

    if (d >= Z->deg)
        return d == Z->deg;

    return Z->pos_degs[d];
}

#ifdef __cplusplus
}
#endif

#endif

