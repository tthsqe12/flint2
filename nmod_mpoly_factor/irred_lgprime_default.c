/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "nmod_mpoly_factor.h"
#include "fq_nmod_mpoly_factor.h"


static void fq_nmod_mpoly_set_nmod_mpoly(
    fq_nmod_mpoly_t A,
    const fq_nmod_mpoly_ctx_t Actx,
    const nmod_mpoly_t B,
    const nmod_mpoly_ctx_t Bctx)
{
    slong i, N;

    FLINT_ASSERT(Actx->minfo->ord == Bctx->minfo->ord);
    FLINT_ASSERT(Actx->minfo->nvars == Bctx->minfo->nvars);

    fq_nmod_mpoly_fit_bits(A, B->bits, Actx);
    A->bits = B->bits;

    N = mpoly_words_per_exp(B->bits, Bctx->minfo);

    fq_nmod_mpoly_fit_length(A, B->length, Actx);
    A->length = B->length;

    mpoly_copy_monomials(A->exps, B->exps, B->length, N);

    for (i = 0; i < B->length; i++)
        fq_nmod_set_ui(A->coeffs + i, B->coeffs[i], Actx->fqctx);
}


static int nmod_mpoly_get_fq_nmod_mpoly(
    nmod_mpoly_t A,
    const nmod_mpoly_ctx_t Actx,
    const fq_nmod_mpoly_t B,
    const fq_nmod_mpoly_ctx_t Bctx)
{
    slong i;
    slong N = mpoly_words_per_exp(B->bits, Bctx->minfo);

    FLINT_ASSERT(Actx->minfo->ord == Bctx->minfo->ord);
    FLINT_ASSERT(Actx->minfo->nvars == Bctx->minfo->nvars);

    nmod_mpoly_fit_length_set_bits(A, B->length, B->bits, Actx);
    A->length = B->length;

    mpoly_copy_monomials(A->exps, B->exps, B->length, N);

    for (i = 0; i < B->length; i++)
    {
        if ((B->coeffs + i)->length == 1)
        {
            A->coeffs[i] = (B->coeffs + i)->coeffs[0];
        }
        else
        {
            return 0;
        }
    }

    return 1;
}


static void get_conjugates(
    fq_nmod_mpoly_struct * C,
    const fq_nmod_mpoly_t A,
    slong deg,
    const fq_nmod_mpoly_ctx_t ctx)
{
    slong i, j;
    for (i = 0; i < deg; i++)
    {
        fq_nmod_mpoly_struct * Ci = C + i;
        fq_nmod_mpoly_set(Ci, A, ctx);
        for (j = 0; j < Ci->length; j++)
        {
            fq_nmod_frobenius(Ci->coeffs + j, Ci->coeffs + j, i, ctx->fqctx);
        }
    }
}


int nmod_mpoly_factor_irred_lgprime_default(
    nmod_mpolyv_t fac_,
    const nmod_mpoly_t A_,
    const nmod_mpoly_ctx_t ctx_)
{
    fq_nmod_mpolyv_t fac;
    fq_nmod_mpoly_t A;
    fq_nmod_mpoly_ctx_t ctx;
    slong edeg;
    int success;
    const slong n = ctx_->minfo->nvars - 1;
    slong i, j;

    FLINT_ASSERT(A_->length > 0);
    FLINT_ASSERT(A_->coeffs[0] == 1);
    FLINT_ASSERT(ctx_->minfo->ord == ORD_LEX);

    edeg = 2;
    fq_nmod_mpoly_ctx_init_deg(ctx, n + 1, ORD_LEX, ctx_->ffinfo->mod.n, edeg);
    fq_nmod_mpoly_init(A, ctx);
    fq_nmod_mpoly_set_nmod_mpoly(A, ctx, A_, ctx_);
    fq_nmod_mpolyv_init(fac, ctx);

try_again:

flint_printf("trying extension of degree %wd\n", edeg);

    success = fq_nmod_mpoly_factor_irred_smprime_default(fac, A, ctx);
    if (success < 1)
    {
        if (success < 0)
        {
            success = 0;
            goto cleanup;
        }

        edeg++;
        fq_nmod_mpoly_ctx_change_modulus(ctx, edeg);
        fq_nmod_mpoly_set_nmod_mpoly(A, ctx, A_, ctx_);
        goto try_again;
    }

cleanup:

    fac_->length = 0;
    if (success)
	{
        nmod_mpoly_t truefactor_;
        fq_nmod_mpoly_t truefactor;
        fq_nmod_mpoly_struct * conjugates;
/*
printf("now must frob combine\n");
*/
	    fq_nmod_mpoly_set_nmod_mpoly(A, ctx, A_, ctx_);
/*
flint_printf("A: "); fq_nmod_mpoly_print_pretty(A, NULL, ctx); printf("\n");
flint_printf("fac: "); fq_nmod_mpoly_factor_print_pretty(fac, NULL, ctx); printf("\n");
*/
    #if WANT_ASSERT
        {
            fq_nmod_mpoly_t t;
            fq_nmod_mpoly_init(t, ctx);
            fq_nmod_mpoly_one(t, ctx);
            for (i = 0; i < fac->length; i++)
                fq_nmod_mpoly_mul(t, t, fac->coeffs + i, ctx);
		    FLINT_ASSERT(fq_nmod_mpoly_equal(A, t, ctx));
            fq_nmod_mpoly_clear(t, ctx);
        }
    #endif

        conjugates = (fq_nmod_mpoly_struct *) flint_malloc(edeg*sizeof(fq_nmod_mpoly_struct));
        for (i = 0; i < edeg; i++)
            fq_nmod_mpoly_init(conjugates + i, ctx);

        nmod_mpoly_init(truefactor_, ctx_);
        fq_nmod_mpoly_init(truefactor, ctx);

        while (fac->length > 0)
        {
            fq_nmod_mpoly_one(truefactor, ctx);
            get_conjugates(conjugates, fac->coeffs + 0, edeg, ctx);

            for (i = 0; i < fac->length; i++)
            for (j = 0; j < edeg; j++)
            {
                if (!fq_nmod_mpoly_equal(fac->coeffs + i, conjugates + j, ctx))
                    continue;

                fq_nmod_mpoly_mul(truefactor, truefactor, fac->coeffs + i, ctx);
                fac->length--;
                fq_nmod_mpoly_swap(fac->coeffs + i, fac->coeffs + fac->length, ctx);
                i--;
				break;
            }

            success = nmod_mpoly_get_fq_nmod_mpoly(truefactor_, ctx_, truefactor, ctx);
            FLINT_ASSERT(success);

            nmod_mpolyv_fit_length(fac_, fac_->length + 1, ctx_);
            nmod_mpoly_swap(fac_->coeffs + fac_->length, truefactor_, ctx_);
            fac_->length++;
        }

        for (i = 0; i < edeg; i++)
            fq_nmod_mpoly_clear(conjugates + i, ctx);
        flint_free(conjugates);

        fq_nmod_mpoly_clear(truefactor, ctx);
        nmod_mpoly_clear(truefactor_, ctx_);

        success = 1;
    }

    fq_nmod_mpoly_clear(A, ctx);
    fq_nmod_mpolyv_clear(fac, ctx);
    fq_nmod_mpoly_ctx_clear(ctx);

    return success;
}

