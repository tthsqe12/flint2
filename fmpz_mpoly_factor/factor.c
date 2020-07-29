/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fmpz_mpoly_factor.h"
#include "fmpq_poly.h"
#include "fmpz_mod_mpoly.h"
#include "nmod_mpoly_factor.h"
#include "mpn_extras.h"
#include "nmod_vec.h"

/*
    fmpz_mod_mpoly_t
    sparse multivariates with fmpz_mod coeffs
*/
typedef struct
{
   fmpz * coeffs;
   ulong * exps;  
   slong alloc;
   slong length;
   flint_bitcnt_t bits;     /* number of bits per exponent */
} fmpz_mod_mpoly_struct;

typedef fmpz_mod_mpoly_struct fmpz_mod_mpoly_t[1];

void fmpz_mod_mpoly_clear(
    fmpz_mod_mpoly_t poly,
    const fmpz_mod_mpoly_ctx_t ctx)
{
   if (poly->coeffs != NULL)
   {
      slong i;

      for (i = 0; i < poly->alloc; i++)
         _fmpz_demote(poly->coeffs + i);

      flint_free(poly->coeffs);
      flint_free(poly->exps);
   }
}


void fmpz_mod_mpoly_init(
    fmpz_mod_mpoly_t A,
    const fmpz_mod_mpoly_ctx_t ctx)
{
    A->coeffs = NULL;
    A->exps = NULL;
    A->alloc = 0;
    A->length = 0;
    A->bits = MPOLY_MIN_BITS;
}

void fmpz_mod_mpoly_init3(
    fmpz_mod_mpoly_t A,
    slong alloc,
    flint_bitcnt_t bits,
    const fmpz_mod_mpoly_ctx_t ctx)
{
    if (alloc > 0)
    {
        slong N = mpoly_words_per_exp(bits, ctx->minfo);
        A->exps   = (ulong *) flint_malloc(alloc*N*sizeof(ulong));
        A->coeffs = (fmpz *) flint_calloc(alloc, sizeof(fmpz));
    }
    else
    {
        alloc = 0;
        A->coeffs = NULL;
        A->exps = NULL;
    }
    A->alloc = alloc;
    A->length = 0;
    A->bits = bits;
}

void fmpz_mod_mpoly_truncate(
    fmpz_mod_mpoly_t A,
    slong newlen,
    const fmpz_mod_mpoly_ctx_t ctx)
{
    if (A->length > newlen)
    {
        slong i;

        for (i = newlen; i < A->length; i++)
            _fmpz_demote(A->coeffs + i);

        A->length = newlen;
    }  
}

void fmpz_mod_mpoly_realloc(
    fmpz_mod_mpoly_t poly,
    slong alloc,
    const fmpz_mod_mpoly_ctx_t ctx)
{
    slong N;

    if (alloc <= 0)             /* Clear up, reinitialise */
    {
        fmpz_mod_mpoly_clear(poly, ctx);
        fmpz_mod_mpoly_init(poly, ctx);
        return;
    }

    N = mpoly_words_per_exp(poly->bits, ctx->minfo);

    if (poly->alloc != 0)            /* Realloc */
    {
        fmpz_mod_mpoly_truncate(poly, alloc, ctx);

        poly->coeffs = (fmpz *) flint_realloc(poly->coeffs, alloc*sizeof(fmpz));
        poly->exps = (ulong *) flint_realloc(poly->exps, alloc*N*sizeof(ulong));

        if (alloc > poly->alloc)
            memset(poly->coeffs + poly->alloc, 0,
                                           (alloc - poly->alloc)*sizeof(fmpz));
    }
    else                        /* Nothing allocated already so do it now */
    {
        poly->coeffs = (fmpz *) flint_calloc(alloc, sizeof(fmpz));
        poly->exps   = (ulong *) flint_malloc(alloc*N*sizeof(ulong));
    }

    poly->alloc = alloc;
}

void fmpz_mod_mpoly_fit_length(
    fmpz_mod_mpoly_t poly,
    slong len,
    const fmpz_mod_mpoly_ctx_t ctx)
{
    if (len > poly->alloc)
    {
        /* At least double number of allocated coeffs */
        if (len < 2 * poly->alloc)
            len = 2 * poly->alloc;
        fmpz_mod_mpoly_realloc(poly, len, ctx);
    }
}



void fmpz_mpoly_convert_perm(
    fmpz_mpoly_t A,
    flint_bitcnt_t Abits,
    const fmpz_mpoly_ctx_t Actx,
    const fmpz_mpoly_t B,
    const fmpz_mpoly_ctx_t Bctx,
    const slong * perm)
{
    slong n = Bctx->minfo->nvars;
    slong m = Actx->minfo->nvars;
    slong i, k, l;
    slong NA, NB;
    ulong * Aexps;
    ulong * Bexps;
    TMP_INIT;

    FLINT_ASSERT(A->bits <= FLINT_BITS);
    FLINT_ASSERT(B->bits <= FLINT_BITS);
    FLINT_ASSERT(B->length > 0);
    TMP_START;

    Aexps = (ulong *) TMP_ALLOC(m*sizeof(ulong));
    Bexps = (ulong *) TMP_ALLOC(n*sizeof(ulong));

    NA = mpoly_words_per_exp(Abits, Actx->minfo);
    NB = mpoly_words_per_exp(B->bits, Bctx->minfo);

    fmpz_mpoly_fit_length_set_bits(A, B->length, Abits, Actx);
    A->length = B->length;
    for (i = 0; i < B->length; i++)
    {        
        fmpz_set(A->coeffs + i, B->coeffs + i);
        mpoly_get_monomial_ui(Bexps, B->exps + NB*i, B->bits, Bctx->minfo);
        for (k = 0; k < m; k++)
        {
            l = perm[k];
            Aexps[k] = l < 0 ? 0 : Bexps[l];
        }
        mpoly_set_monomial_ui(A->exps + NA*i, Aexps, Abits, Actx->minfo);
     }  
    fmpz_mpoly_sort_terms(A, Actx);
    TMP_END;
}


/*
    A is primitive w.r.t to any variable appearing in A.
    A is square free with positive lead coeff.

    return 1 for success, 0 for failure
*/
static int _irreducible_factors(
    fmpz_mpolyv_t Af,
    fmpz_mpoly_t A,
    const fmpz_mpoly_ctx_t ctx)
{
    int success;
    slong i, mvars;
    slong * Adegs, * perm, * iperm;
    slong nvars = ctx->minfo->nvars;
    flint_bitcnt_t Lbits, Abits;
    int perm_is_id;
    flint_rand_t state;
#if WANT_ASSERT
    fmpz_mpoly_t Aorg;

    fmpz_mpoly_init(Aorg, ctx);
    fmpz_mpoly_set(Aorg, A, ctx);
#endif

    FLINT_ASSERT(A->length == 0 || fmpz_sgn(A->coeffs + 0) > 0);

    flint_randinit(state);
    Adegs = (slong *) flint_malloc(3*nvars*sizeof(slong));
    perm = Adegs + nvars;
    iperm = perm + nvars;

    if (A->length < 2)
    {
        FLINT_ASSERT(A->length == 1);
        FLINT_ASSERT(!fmpz_mpoly_is_fmpz(A, ctx));
        fmpz_mpolyv_fit_length(Af, 1, ctx);
        Af->length = 1;
        fmpz_mpoly_swap(Af->coeffs + 0, A, ctx);
        success = 1;
        goto cleanup;
    }

    if (!fmpz_mpoly_degrees_fit_si(A, ctx))
    {
        success = 0;
        goto cleanup;
    }

    if (A->bits > FLINT_BITS &&
        !fmpz_mpoly_repack_bits_inplace(A, FLINT_BITS, ctx))
    {
        success = 0;
        goto cleanup;
    }

    fmpz_mpoly_degrees_si(Adegs, A, ctx);

    Abits = A->bits;

    mvars = 0;
    Lbits = 0;
    for (i = 0; i < nvars; i++)
    {
        iperm[i] = -1;
        if (Adegs[i] > 0)
        {
            flint_bitcnt_t this_bits = FLINT_BIT_COUNT(Adegs[i]);
            Lbits = FLINT_MAX(Lbits, this_bits);
            perm[mvars] = i;
            mvars++;
        }
    }

    if (Lbits > FLINT_BITS - 10)
    {
        success = 0;
        goto cleanup;
    }

    /* TODO nice permutation */

    /* invert perm */
    perm_is_id = (mvars == nvars);
    for (i = 0; i < mvars; i++)
    {
        perm_is_id = perm_is_id && (perm[i] == i);
        iperm[perm[i]] = i;
    }

    if (mvars < 2)
    {
        fmpz_poly_t Au;
        fmpz_poly_factor_t Auf;

        FLINT_ASSERT(mvars == 1);

        fmpz_poly_init(Au);
        fmpz_poly_factor_init(Auf);

        FLINT_ASSERT(fmpz_mpoly_is_fmpz_poly(A, perm[0], ctx));
        success = fmpz_mpoly_get_fmpz_poly(Au, A, perm[0], ctx);
        FLINT_ASSERT(success);
        fmpz_poly_factor(Auf, Au);

        /* c = +-1 */
        FLINT_ASSERT(fmpz_is_pm1(&Auf->c));

        fmpz_mpolyv_fit_length(Af, Auf->num, ctx);
        Af->length = Auf->num;
        for (i = 0; i < Auf->num; i++)
        {
            FLINT_ASSERT(Auf->exp[i] == 1);
            _fmpz_mpoly_set_fmpz_poly(Af->coeffs + i, Abits,
                             Auf->p[i].coeffs, Auf->p[i].length, perm[0], ctx);
        }

        fmpz_poly_clear(Au);
        fmpz_poly_factor_clear(Auf);

        success = 1;
    }
    else if (mvars == 2)
    {
        fmpz_poly_t c;
        fmpz_bpoly_t Ab;
        fmpz_tpoly_t Abf;

        fmpz_poly_init(c);
        fmpz_bpoly_init(Ab);
        fmpz_tpoly_init(Abf);

        fmpz_mpoly_get_bpoly(Ab, A, perm[0], perm[1], ctx);
        fmpz_bpoly_factor(c, Abf, Ab);

        FLINT_ASSERT(c->length == 1 && fmpz_is_pm1(c->coeffs + 0));

        fmpz_mpolyv_fit_length(Af, Abf->length, ctx);
        Af->length = Abf->length;
        for (i = 0; i < Abf->length; i++)
        {
            fmpz_mpoly_set_fmpz_bpoly(Af->coeffs + i, Abits, Abf->coeffs + i,
                                                        perm[0], perm[1], ctx);
            fmpz_mpoly_unit_normalize(Af->coeffs + i, ctx);
        }

        fmpz_poly_clear(c);
        fmpz_bpoly_clear(Ab);
        fmpz_tpoly_clear(Abf);

        success = 1;
    }
    else
    {
        fmpz_mpoly_ctx_t Lctx;
        fmpz_mpoly_t L, lcL;
        fmpz_mpolyv_t Lf;
        fmpz_mpoly_factor_t lcLf;

        fmpz_mpoly_ctx_init(Lctx, mvars, ORD_LEX);
        fmpz_mpoly_init(L, Lctx);
        fmpz_mpoly_init(lcL, Lctx);
        fmpz_mpolyv_init(Lf, Lctx);
        fmpz_mpoly_factor_init(lcLf, Lctx);

        Lbits = mpoly_fix_bits(Lbits + 1, Lctx->minfo);

        fmpz_mpoly_convert_perm(L, Lbits, Lctx, A, ctx, perm);
        fmpz_mpoly_unit_normalize(L, ctx);

        _fmpz_mpoly_get_lead0(lcL, L, Lctx);
        success = fmpz_mpoly_factor(lcLf, lcL, Lctx);
        if (success)
        {
            success = fmpz_mpoly_factor_irred_zippel(Lf, L, lcLf, lcL, Lctx, state);
            if (!success)
                success = fmpz_mpoly_factor_irred_wang(Lf, L, lcLf, lcL, Lctx, state);
        }
        if (!success)
            success = fmpz_mpoly_factor_irred_default(Lf, L, Lctx);

        if (success)
        {
            fmpz_mpolyv_fit_length(Af, Lf->length, ctx);
            Af->length = Lf->length;
            for (i = 0; i < Lf->length; i++)
            {
                fmpz_mpoly_convert_perm(Af->coeffs + i, Abits, ctx,
                                                  Lf->coeffs + i, Lctx, iperm);
                fmpz_mpoly_unit_normalize(Af->coeffs + i, ctx);
            }
        }

        fmpz_mpoly_clear(L, Lctx);
        fmpz_mpoly_clear(lcL, Lctx);
        fmpz_mpolyv_clear(Lf, Lctx);
        fmpz_mpoly_factor_clear(lcLf, Lctx);
        fmpz_mpoly_ctx_clear(Lctx);
    }

cleanup:

    flint_randclear(state);
    flint_free(Adegs);

#if WANT_ASSERT
    if (success)
    {
        fmpz_mpoly_t prod;
        fmpz_mpoly_init(prod, ctx);
        fmpz_mpoly_one(prod, ctx);
        for (i = 0; i < Af->length; i++)
            fmpz_mpoly_mul(prod, prod, Af->coeffs + i, ctx);
        FLINT_ASSERT(fmpz_mpoly_equal(prod, Aorg, ctx));
        fmpz_mpoly_clear(prod, ctx);
        fmpz_mpoly_clear(Aorg, ctx);
    }
#endif

    return success;
}


int fmpz_mpoly_factor(
    fmpz_mpoly_factor_t f,
    const fmpz_mpoly_t A,
    const fmpz_mpoly_ctx_t ctx)
{
    int success;
    slong i, j;
    fmpz_mpolyv_t t;
    fmpz_mpoly_factor_t g;

    fmpz_mpolyv_init(t, ctx);
    fmpz_mpoly_factor_init(g, ctx);

    success = fmpz_mpoly_factor_squarefree(f, A, ctx);
    if (!success)
        goto cleanup;

    FLINT_ASSERT(fmpz_mpoly_factor_matches(A, f, ctx));

    fmpz_swap(g->constant, f->constant);
    g->num = 0;
    for (j = 0; j < f->num; j++)
    {
        success = _irreducible_factors(t, f->poly + j, ctx);
        if (!success)
            goto cleanup;

        fmpz_mpoly_factor_fit_length(g, g->num + t->length, ctx);
        for (i = 0; i < t->length; i++)
        {
            fmpz_set(g->exp + g->num, f->exp + j);
            fmpz_mpoly_swap(g->poly + g->num, t->coeffs + i, ctx);
            g->num++;
        }
    }
    fmpz_mpoly_factor_swap(f, g, ctx);

    success = 1;

cleanup:

    fmpz_mpolyv_clear(t, ctx);
    fmpz_mpoly_factor_clear(g, ctx);

    FLINT_ASSERT(!success || fmpz_mpoly_factor_matches(A, f, ctx));

    return success;
}
