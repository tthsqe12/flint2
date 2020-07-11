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
    fmpz_mpoly_factor_t Af,
    fmpz_mpoly_t A,
    const fmpz_mpoly_ctx_t ctx)
{
    int success;
    slong i, mvars;
    slong * Adegs, * perm, * iperm;
    slong nvars = ctx->minfo->nvars;
    flint_bitcnt_t Lbits, Abits = A->bits;
    int same;
/*
flint_printf("_irreducible_factors called\n");
flint_printf("A: "); fmpz_mpoly_print_pretty(A, NULL, ctx); flint_printf("\n");
*/

    FLINT_ASSERT(A->length == 0 || fmpz_sgn(A->coeffs + 0) > 0);

    Adegs = (slong *) flint_malloc(3*nvars*sizeof(slong));
    perm = Adegs + nvars;
    iperm = perm + nvars;

    if (A->length <= 1)
    {
        if (A->length == 1)
        {
            FLINT_ASSERT(!fmpz_mpoly_is_fmpz(A, ctx));
            fmpz_mpoly_factor_fit_length(Af, 1, ctx);
            Af->num = 1;
            fmpz_one(Af->constant);
            fmpz_mpoly_set(Af->poly + 0, A, ctx);
            fmpz_one(Af->exp + 0);
        }
        else
        {
            fmpz_mpoly_factor_zero(Af, ctx);
        }

        success = 1;
        goto cleanup;
    }

    if (!fmpz_mpoly_degrees_fit_si(A, ctx))
    {
        success = 0;
        goto cleanup;
    }

    fmpz_mpoly_degrees_si(Adegs, A, ctx);

    mvars = 0;
    Lbits = 0;
    for (i = 0; i < ctx->minfo->nvars; i++)
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

    if (Lbits + 10 > FLINT_BITS)
    {
        success = 0;
        goto cleanup;
    }

    /* TODO: figure out nice perm */

    /* invert perm */
    same = (mvars == nvars);
    for (i = 0; i < mvars; i++)
    {
        same = same && (perm[i] == i);
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

        fmpz_mpoly_factor_fit_length(Af, Auf->num, ctx);
        fmpz_swap(Af->constant, &Auf->c); /* Auf->c should be 1 */
        for (i = 0; i < Auf->num; i++)
        {
            fmpz_set_ui(Af->exp + i, Auf->exp[i]); /* Aufac->exp[i] should be 1 */
            _fmpz_mpoly_set_fmpz_poly(Af->poly + i, Abits,
                             Auf->p[i].coeffs, Auf->p[i].length, perm[0], ctx);
        }
        Af->num = Auf->num;

        fmpz_poly_clear(Au);
        fmpz_poly_factor_clear(Auf);

        success = 1;
    }
    else if (mvars == 2)
    {
        fmpz_mpoly_factor_irred_bivar(Af, A, perm[0], perm[1], ctx);
        success = 1;
    }
    else
    {
        fmpz_mpoly_ctx_t Lctx;
        fmpz_mpoly_t L, lcL;
        fmpz_mpoly_factor_t Lf, lcLf;

        fmpz_mpoly_ctx_init(Lctx, mvars, ORD_LEX);
        fmpz_mpoly_init(L, Lctx);
        fmpz_mpoly_init(lcL, Lctx);
        fmpz_mpoly_factor_init(Lf, Lctx);
        fmpz_mpoly_factor_init(lcLf, Lctx);

        Lbits = mpoly_fix_bits(Lbits + 1, Lctx->minfo);

        fmpz_mpoly_convert_perm(L, Lbits, Lctx, A, ctx, perm);
        fmpz_mpoly_unit_normalize(L, ctx);

        _fmpz_mpoly_get_lead0(lcL, L, Lctx);
        success = fmpz_mpoly_factor(lcLf, lcL, Lctx);
        if (success)
        {
            success = fmpz_mpoly_factor_irred_tuncer(Lf, L, lcLf, lcL, Lctx);
            if (!success)
                success = fmpz_mpoly_factor_irred_wang(Lf, L, lcLf, lcL, Lctx);
        }
        if (!success)
            success = fmpz_mpoly_factor_irred_default(Lf, L, Lctx);

        if (success)
        {
            fmpz_mpoly_factor_fit_length(Af, Lf->num, ctx);
            fmpz_one(Af->constant);
            Af->num = Lf->num;
            for (i = 0; i < Lf->num; i++)
            {
                fmpz_one(Af->exp + i);
                fmpz_mpoly_convert_perm(Af->poly + i, Abits, ctx,
                                                    Lf->poly + i, Lctx, iperm);
                fmpz_mpoly_unit_normalize(Af->poly + i, ctx);
            }
        }

        fmpz_mpoly_clear(L, Lctx);
        fmpz_mpoly_clear(lcL, Lctx);
        fmpz_mpoly_factor_clear(Lf, Lctx);
        fmpz_mpoly_factor_clear(lcLf, Lctx);
        fmpz_mpoly_ctx_clear(Lctx);
    }

cleanup:

    flint_free(Adegs);
/*
flint_printf("_irreducible_factors returning %d\n", success);
flint_printf("A: "); fmpz_mpoly_print_pretty(A, NULL, ctx); flint_printf("\n");
flint_printf("f: "); fmpz_mpoly_factor_print_pretty(f, NULL, ctx); flint_printf("\n");
*/
    FLINT_ASSERT(!success || fmpz_mpoly_factor_matches(A, Af, ctx));

    return success;
}


int fmpz_mpoly_factor(
    fmpz_mpoly_factor_t f,
    const fmpz_mpoly_t A,
    const fmpz_mpoly_ctx_t ctx)
{
    int success;
    slong i, j;
    fmpz_mpoly_factor_t newf, tempf;

    fmpz_mpoly_factor_init(newf, ctx);
    fmpz_mpoly_factor_init(tempf, ctx);

    success = fmpz_mpoly_factor_squarefree(f, A, ctx);
    if (!success)
        goto cleanup;

    /* ensure factors are irreducible */
    FLINT_ASSERT(fmpz_mpoly_factor_matches(A, f, ctx));

    fmpz_swap(newf->constant, f->constant);
    newf->num = 0;
    for (j = 0; j < f->num; j++)
    {
        success = _irreducible_factors(tempf, f->poly + j, ctx);
        if (!success)
            goto cleanup;

        FLINT_ASSERT(fmpz_is_one(tempf->constant));
        for (i = 0; i < tempf->num; i++)
        {
            FLINT_ASSERT(fmpz_is_one(tempf->exp + i));
            fmpz_mpoly_factor_append_fmpz_swap(newf,
                                             tempf->poly + i, f->exp + j, ctx);
        }
    }
    fmpz_mpoly_factor_swap(f, newf, ctx);

    success = 1;

cleanup:

    fmpz_mpoly_factor_clear(newf, ctx);
    fmpz_mpoly_factor_clear(tempf, ctx);

    FLINT_ASSERT(fmpz_mpoly_factor_matches(A, f, ctx));

    return success;
}

