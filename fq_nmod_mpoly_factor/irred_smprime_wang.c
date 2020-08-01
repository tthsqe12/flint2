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
#include "ui_factor.h"


static void _eval_rest(
    fq_nmod_poly_t E,
    const fq_nmod_struct * Acoeffs,
    const ulong * Aexps,
    slong Alen,
    slong var,
    const fq_nmod_poly_struct * alphabetas,
    const slong * offsets,
    const slong * shifts,
    slong N,
    ulong mask,
    slong nvars,
    const fq_nmod_ctx_t ctx)
{
    slong offset = offsets[var];
    slong shift = shifts[var];
    slong start, stop;
    ulong e, next_e;
    fq_nmod_poly_t T;

    fq_nmod_poly_zero(E, ctx);

    if (Alen < 1)
        return;

    if (var >= nvars)
    {
        FLINT_ASSERT(Alen == 1);
        fq_nmod_poly_set_fq_nmod(E, Acoeffs + 0, ctx);
        return;
    }

    fq_nmod_poly_init(T, ctx);

    start = 0;
    e = mask & (Aexps[N*start + offset] >> shift);

next:

    FLINT_ASSERT(start < Alen);
    FLINT_ASSERT(e == (mask & (Aexps[N*start + offset] >> shift)));

    stop = start + 1;
    while (stop < Alen && (mask & (Aexps[N*stop + offset] >> shift)) == e)
        stop++;

    _eval_rest(T, Acoeffs + start, Aexps + N*start, stop - start, var + 1,
                         alphabetas + 1, offsets, shifts, N, mask, nvars, ctx);
    fq_nmod_poly_add(E, E, T, ctx);

    if (stop < Alen)
    {
        next_e = (mask & (Aexps[N*stop + offset] >> shift));
        FLINT_ASSERT(next_e < e);
        fq_nmod_poly_pow(T, alphabetas, e - next_e, ctx);
        fq_nmod_poly_mul(E, E, T, ctx);
        e = next_e;
        start = stop;
        goto next;
    }
    else
    {
        fq_nmod_poly_pow(T, alphabetas, e, ctx);
        fq_nmod_poly_mul(E, E, T, ctx);
    }

    fq_nmod_poly_clear(T, ctx);
}

void _fq_nmod_eval_to_bpoly(
    fq_nmod_bpoly_t E,
    const fq_nmod_mpoly_t A,
    const fq_nmod_poly_struct * alphabetas,
    const fq_nmod_mpoly_ctx_t ctx)
{
    slong i, N = mpoly_words_per_exp_sp(A->bits, ctx->minfo);
    slong * offsets, * shifts;
    slong offset, shift;
    slong start, stop;
    ulong e, mask = (-UWORD(1)) >> (FLINT_BITS - A->bits);

    E->length = 0;
    if (A->length < 1)
        return;

    offsets = FLINT_ARRAY_ALLOC(ctx->minfo->nvars, slong);
    shifts  = FLINT_ARRAY_ALLOC(ctx->minfo->nvars, slong);
    for (i = 0; i < ctx->minfo->nvars; i++)
        mpoly_gen_offset_shift_sp(offsets + i, shifts + i, i, A->bits, ctx->minfo);

    offset = offsets[0];
    shift = shifts[0];

    start = 0;
    e = mask & (A->exps[N*start + offset] >> shift);

next:

    FLINT_ASSERT(start < A->length);
    FLINT_ASSERT(e == (mask & (A->exps[N*start + offset] >> shift)));

    stop = start + 1;
    while (stop < A->length && (mask & (A->exps[N*stop + offset] >> shift)) == e)
        stop++;

    fq_nmod_bpoly_fit_length(E, e + 1, ctx->fqctx);
    while (E->length <= e)
    {
        fq_nmod_poly_zero(E->coeffs + E->length, ctx->fqctx);
        E->length++;
    }

    _eval_rest(E->coeffs + e, A->coeffs + start, A->exps + N*start, stop - start, 1,
          alphabetas, offsets, shifts, N, mask, ctx->minfo->nvars, ctx->fqctx);

    if (stop < A->length)
    {
        FLINT_ASSERT(e > (mask & (A->exps[N*stop + offset] >> shift)));
        e = (mask & (A->exps[N*stop + offset] >> shift));
        start = stop;
        goto next;
    }

    fq_nmod_bpoly_normalise(E, ctx->fqctx);

    flint_free(offsets);
    flint_free(shifts);
}


/* A = B(gen(var), 0) */
void _fq_nmod_mpoly_set_bpoly_var1_zero(
    fq_nmod_mpoly_t A,
    flint_bitcnt_t Abits,
    const fq_nmod_bpoly_t B,
    slong var,
    const fq_nmod_mpoly_ctx_t ctx)
{
    slong N = mpoly_words_per_exp(Abits, ctx->minfo);
    slong i, Alen;
    slong Blen = B->length;
    ulong * genexp;
    TMP_INIT;

    TMP_START;

    genexp = (ulong *) TMP_ALLOC(N*sizeof(ulong));
    if (Abits <= FLINT_BITS)
        mpoly_gen_monomial_sp(genexp, var, Abits, ctx->minfo);
    else
        mpoly_gen_monomial_offset_mp(genexp, var, Abits, ctx->minfo);

    Alen = 2;
    for (i = 0; i < Blen; i++)
        Alen += (B->coeffs[i].length > 0);

    fq_nmod_mpoly_fit_length_set_bits(A, Alen, Abits, ctx);

    Alen = 0;
    for (i = Blen - 1; i >= 0; i--)
    {
        FLINT_ASSERT(Alen < A->alloc);
        fq_nmod_bpoly_get_coeff(A->coeffs + Alen, B, i, 0, ctx->fqctx);
        if (fq_nmod_is_zero(A->coeffs + Alen, ctx->fqctx))
            continue;

        if (Abits <= FLINT_BITS)
            mpoly_monomial_mul_ui(A->exps + N*Alen, genexp, N, i);
        else
            mpoly_monomial_mul_ui_mp(A->exps + N*Alen, genexp, N, i);
        Alen++;
    }
    A->length = Alen;

    TMP_END;
}

int fq_nmod_mpoly_factor_irred_smprime_wang(
    fq_nmod_mpolyv_t fac,
    const fq_nmod_mpoly_t A,
    const fq_nmod_mpoly_factor_t lcAfac,
    const fq_nmod_mpoly_t lcA,
    const fq_nmod_mpoly_ctx_t ctx,
    flint_rand_t state)
{
    int success;
    int alphas_tries_remaining, alphabetas_tries_remaining, alphabetas_length;
    const slong n = ctx->minfo->nvars - 1;
    slong i, j, k, r;
    fq_nmod_struct * alpha;
    fq_nmod_poly_struct * alphabetas;
    fq_nmod_mpoly_struct * Aevals;
    slong * degs, * degeval;
    fq_nmod_mpolyv_t tfac;
    fq_nmod_mpoly_t t, Acopy;
    fq_nmod_mpoly_struct * newA;
    fq_nmod_poly_t Abfc;
    fq_nmod_bpoly_t Ab;
    fq_nmod_tpoly_t Abfp;
    fq_nmod_mpoly_t m, mpow;
    fq_nmod_mpolyv_t new_lcs, lc_divs;
/*
flint_printf("fq_nmod_mpoly_factor_irred_smprime_wang called p = %wu\n", ctx->fqctx->modulus->mod.n);
flint_printf("     A: "); fq_nmod_mpoly_print_pretty(A, NULL, ctx); flint_printf("\n");
flint_printf("lcAfac: "); fq_nmod_mpoly_factor_print_pretty(lcAfac, NULL, ctx); flint_printf("\n");
flint_printf("   lcA: "); fq_nmod_mpoly_print_pretty(lcA, NULL, ctx); flint_printf("\n");
*/
    FLINT_ASSERT(n > 1);
    FLINT_ASSERT(A->length > 1);
    FLINT_ASSERT(fq_nmod_is_one(A->coeffs + 0, ctx->fqctx));
    FLINT_ASSERT(A->bits <= FLINT_BITS);

    fq_nmod_mpoly_init(Acopy, ctx);
    fq_nmod_mpoly_init(m, ctx);
    fq_nmod_mpoly_init(mpow, ctx);

    fq_nmod_mpolyv_init(new_lcs, ctx);
    fq_nmod_mpolyv_init(lc_divs, ctx);

    fq_nmod_poly_init(Abfc, ctx->fqctx);
    fq_nmod_tpoly_init(Abfp, ctx->fqctx);
    fq_nmod_bpoly_init(Ab, ctx->fqctx);

    degs    = FLINT_ARRAY_ALLOC(n + 1, slong);
    degeval = FLINT_ARRAY_ALLOC(n + 1, slong);
	alpha   = FLINT_ARRAY_ALLOC(n, fq_nmod_struct);
    alphabetas = FLINT_ARRAY_ALLOC(n, fq_nmod_poly_struct);
    Aevals  = FLINT_ARRAY_ALLOC(n, fq_nmod_mpoly_struct);
	for (i = 0; i < n; i++)
    {
        fq_nmod_init(alpha + i, ctx->fqctx);
        fq_nmod_poly_init(alphabetas + i, ctx->fqctx);
		fq_nmod_mpoly_init(Aevals + i, ctx);
    }
    fq_nmod_mpolyv_init(tfac, ctx);
	fq_nmod_mpoly_init(t, ctx);

    /* init done */

    alphabetas_length = 2;
    alphas_tries_remaining = 10;
	fq_nmod_mpoly_degrees_si(degs, A, ctx);

next_alpha:

    if (--alphas_tries_remaining < 0)
	{
		success = 0;
        goto cleanup;
	}

    for (i = 0; i < n; i++)
        fq_nmod_rand(alpha + i, state, ctx->fqctx);

    /* ensure degrees do not drop under evaluation */
	for (i = n - 1; i >= 0; i--)
	{
        fq_nmod_mpoly_evaluate_one_fq_nmod(Aevals + i,
                       i == n - 1 ? A : Aevals + i + 1, i + 1, alpha + i, ctx);
		fq_nmod_mpoly_degrees_si(degeval, Aevals + i, ctx);
		for (j = 0; j <= i; j++)
			if (degeval[j] != degs[j])
				goto next_alpha;
	}

    /* make sure univar is squarefree */
	fq_nmod_mpoly_derivative(t, Aevals + 0, 0, ctx);
	fq_nmod_mpoly_gcd(t, t, Aevals + 0, ctx);
	if (!fq_nmod_mpoly_is_one(t, ctx))
		goto next_alpha;

    alphabetas_tries_remaining = 2 + alphabetas_length;

next_alphabetas:

    if (--alphabetas_tries_remaining < 0)
    {
        if (++alphabetas_length > 10)
        {
            success = 0;
            goto cleanup;
        }
        goto next_alpha;
    }

    for (i = 0; i < n; i++)
    {
        fq_nmod_poly_fit_length(alphabetas + i, alphabetas_length, ctx->fqctx);
        fq_nmod_set(alphabetas[i].coeffs + 0, alpha + i, ctx->fqctx);
        for (j = 1; j < alphabetas_length; j++)
            fq_nmod_rand(alphabetas[i].coeffs + j, state, ctx->fqctx);
        alphabetas[i].length = alphabetas_length;
        _fq_nmod_poly_normalise(alphabetas + i, ctx->fqctx);
/*
flint_printf("alphabetas[%wd]: ", i);
fq_nmod_poly_print_pretty(alphabetas + i, "Y", ctx->fqctx);
flint_printf("\n");
*/
    }

    _fq_nmod_eval_to_bpoly(Ab, A, alphabetas, ctx);
    success = fq_nmod_bpoly_factor_smprime(Abfc, Abfp, Ab, 0, ctx->fqctx);
    if (!success)
    {
        FLINT_ASSERT(0 && "this should not happen");
        goto next_alpha;
    }

    r = Abfp->length;

    if (r < 2)
    {
        fq_nmod_mpolyv_fit_length(fac, 1, ctx);
        fac->length = 1;
        fq_nmod_mpoly_set(fac->coeffs + 0, A, ctx);
        success = 1;
        goto cleanup;
    }

    fq_nmod_mpolyv_fit_length(lc_divs, r, ctx);
    lc_divs->length = r;
    if (lcAfac->num > 0)
    {
        success = fq_nmod_mpoly_factor_lcc_wang(lc_divs->coeffs, lcAfac,
                                       Abfc, Abfp->coeffs, r, alphabetas, ctx);
        if (!success)
            goto next_alphabetas;
    }
    else
    {
        for (i = 0; i < r; i++)
            fq_nmod_mpoly_one(lc_divs->coeffs + i, ctx);
    }
/*
flint_printf("lc_divs:\n");
fq_nmod_mpolyv_print_pretty(lc_divs, NULL, ctx);
*/

    success = fq_nmod_mpoly_divides(m, lcA, lc_divs->coeffs + 0, ctx);
    FLINT_ASSERT(success);
    for (i = 1; i < r; i++)
    {
        success = fq_nmod_mpoly_divides(m, m, lc_divs->coeffs + i, ctx);
        FLINT_ASSERT(success);
    }

    fq_nmod_mpoly_pow_ui(mpow, m, r - 1, ctx);
    if (fq_nmod_mpoly_is_one(mpow, ctx))
    {
        newA = (fq_nmod_mpoly_struct *) A;
    }
    else
    {
        newA = Acopy;
        fq_nmod_mpoly_mul(newA, A, mpow, ctx);
    }

    if (newA->bits > FLINT_BITS)
    {
        success = 0;
        goto cleanup;
    }

    fq_nmod_mpoly_degrees_si(degs, newA, ctx);

    fq_nmod_mpoly_set(t, mpow, ctx);
    for (i = n - 1; i >= 0; i--)
    {
        fq_nmod_mpoly_evaluate_one_fq_nmod(t, mpow, i + 1, alpha + i, ctx);
        fq_nmod_mpoly_swap(t, mpow, ctx);
        fq_nmod_mpoly_mul(Aevals + i, Aevals + i, mpow, ctx);
    }

    fq_nmod_mpolyv_fit_length(new_lcs, (n + 1)*r, ctx);
    i = n;
    for (j = 0; j < r; j++)
    {
        fq_nmod_mpoly_mul(new_lcs->coeffs + i*r + j, lc_divs->coeffs + j, m, ctx);
    }
    for (i = n - 1; i >= 0; i--)
    {
        for (j = 0; j < r; j++)
        {
            fq_nmod_mpoly_evaluate_one_fq_nmod(new_lcs->coeffs + i*r + j,
                       new_lcs->coeffs + (i + 1)*r + j, i + 1, alpha + i, ctx);
        }
    }

    fq_nmod_mpolyv_fit_length(fac, r, ctx);
    fac->length = r;
    for (i = 0; i < r; i++)
    {
        fq_nmod_t q;
        fq_nmod_init(q, ctx->fqctx);
        FLINT_ASSERT(fq_nmod_mpoly_is_fq_nmod(new_lcs->coeffs + 0*r + i, ctx));
        FLINT_ASSERT(fq_nmod_mpoly_length(new_lcs->coeffs + 0*r + i, ctx) == 1);
        _fq_nmod_mpoly_set_bpoly_var1_zero(fac->coeffs + i, newA->bits, Abfp->coeffs + i, 0, ctx);
        FLINT_ASSERT(fac->coeffs[i].length > 0);
        fq_nmod_inv(q, fac->coeffs[i].coeffs + 0, ctx->fqctx);
        fq_nmod_mul(q, q, new_lcs->coeffs[0*r + i].coeffs + 0, ctx->fqctx);
        fq_nmod_mpoly_scalar_mul_fq_nmod(fac->coeffs + i, fac->coeffs + i, q, ctx);
        fq_nmod_clear(q, ctx->fqctx);
    }

    fq_nmod_mpolyv_fit_length(tfac, r, ctx);
    tfac->length = r;
    for (k = 1; k <= n; k++)
    {
        for (i = 0; i < r; i++)
        {
            _fq_nmod_mpoly_set_lead0(tfac->coeffs + i, fac->coeffs + i,
                                               new_lcs->coeffs + k*r + i, ctx);
        }

        success = fq_nmod_mpoly_hlift(k, tfac->coeffs, r, alpha,
                                         k < n ? Aevals + k : newA, degs, ctx);

        if (!success)
            goto next_alphabetas;

        fq_nmod_mpolyv_swap(tfac, fac, ctx);
    }

    if (!fq_nmod_mpoly_is_fq_nmod(m, ctx))
    {
        fq_nmod_mpoly_univar_t u;
        fq_nmod_mpoly_univar_init(u, ctx);
        for (i = 0; i < r; i++)
        {
            fq_nmod_mpoly_to_univar(u, fac->coeffs + i, 0, ctx);
            success = fq_nmod_mpoly_univar_content_mpoly(t, u, ctx);
            if (!success)
            {
                fq_nmod_mpoly_univar_clear(u, ctx);
                goto cleanup;
            }
            success = fq_nmod_mpoly_divides(fac->coeffs + i,
                                            fac->coeffs + i, t, ctx);
            FLINT_ASSERT(success);
        }
        fq_nmod_mpoly_univar_clear(u, ctx);
    }

    for (i = 0; i < r; i++)
        fq_nmod_mpoly_make_monic(fac->coeffs + i, fac->coeffs + i, ctx);

    success = 1;

cleanup:

    fq_nmod_mpolyv_clear(new_lcs, ctx);
    fq_nmod_mpolyv_clear(lc_divs, ctx);

    fq_nmod_poly_clear(Abfc, ctx->fqctx);
    fq_nmod_tpoly_clear(Abfp, ctx->fqctx);

	for (i = 0; i < n; i++)
    {
		fq_nmod_mpoly_clear(Aevals + i, ctx);
        fq_nmod_poly_clear(alphabetas + i, ctx->fqctx);
    }
    flint_free(alphabetas);
    flint_free(alpha);
    flint_free(Aevals);
    flint_free(degs);
    flint_free(degeval);

    fq_nmod_mpolyv_clear(tfac, ctx);
    fq_nmod_mpoly_clear(t, ctx);

    fq_nmod_mpoly_clear(Acopy, ctx);
    fq_nmod_mpoly_clear(m, ctx);
    fq_nmod_mpoly_clear(mpow, ctx);

#if WANT_ASSERT
    if (success)
    {
        fq_nmod_mpoly_t prod;
        fq_nmod_mpoly_init(prod, ctx);
        fq_nmod_mpoly_one(prod, ctx);
        for (i = 0; i < fac->length; i++)
            fq_nmod_mpoly_mul(prod, prod, fac->coeffs + i, ctx);
        FLINT_ASSERT(fq_nmod_mpoly_equal(prod, A, ctx));
        fq_nmod_mpoly_clear(prod, ctx);
    }
#endif
/*
flint_printf("fq_nmod_mpoly_factor_irred_smprime_wang returning %d\n", success);
*/
	return success;
}
