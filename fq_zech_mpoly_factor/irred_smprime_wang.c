/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fq_zech_mpoly_factor.h"

/*
    only E and alphas are shifted by "var"
    so output is in E[0]
*/
int _fq_zech_mpoly_evaluate_rest_n_poly_fq(
    fq_zech_poly_struct * E,
    slong * starts,
    slong * ends,
    slong * stops,
    ulong * es,
    const fq_zech_struct * Acoeffs,
    const ulong * Aexps,
    slong Alen,
    slong var,
    const fq_zech_poly_struct * alphas,
    const slong * offsets,
    const slong * shifts,
    slong N,
    ulong mask,
    slong nvars,
    const fq_zech_ctx_t ctx)
{
    slong v, stop;
    ulong next_e;

    FLINT_ASSERT(var < nvars);

    E -= var;
    alphas -= var;

    v = var;
    starts[v] = 0;
    ends[v] = Alen;
    fq_zech_poly_zero(E + v, ctx);

    if (Alen < 1)
        return 1;

calculate:
/*
    input:
        v
        starts[v]
        ends[v]
*/
    FLINT_ASSERT(ends[v] > starts[v]);
    es[v] = mask & (Aexps[N*starts[v] + offsets[v]] >> shifts[v]);

    fq_zech_poly_zero(E + v, ctx);

next:

    FLINT_ASSERT(es[v] == (mask & (Aexps[N*starts[v] + offsets[v]] >> shifts[v])));

    stop = starts[v] + 1;
    while (stop < ends[v] &&
           (mask & (Aexps[N*stop + offsets[v]] >> shifts[v])) == es[v])
    {
        stop++;
    }
    stops[v] = stop;

    if (v + 1 < nvars)
    {
        starts[v + 1] = starts[v];
        ends[v + 1] = stops[v];
        v++;
        goto calculate;
calculate_return:
        fq_zech_poly_add(E + v, E + v, E + v + 1, ctx);
    }
    else
    {
        fq_zech_poly_set_fq_zech(E + v + 1, Acoeffs + starts[v], ctx);
        fq_zech_poly_add(E + v, E + v, E + v + 1, ctx);
    }

    if (stops[v] < ends[v])
    {
        next_e = mask & (Aexps[N*stops[v] + offsets[v]] >> shifts[v]);
        FLINT_ASSERT(next_e < es[v]);
        fq_zech_poly_pow(E + v + 1, alphas + v, es[v] - next_e, ctx);
        fq_zech_poly_mul(E + v, E + v, E + v + 1, ctx);
        es[v] = next_e;
        starts[v] = stops[v];
        goto next;
    }
    else
    {
        fq_zech_poly_pow(E + v + 1, alphas + v, es[v], ctx);
        fq_zech_poly_mul(E + v, E + v, E + v + 1, ctx);
    }

    if (v > var)
    {
        v--;
        goto calculate_return;
    }

    return 1;
}


void _fq_zech_eval_to_bpoly(
    fq_zech_bpoly_t E,
    const fq_zech_mpoly_t A,
    const fq_zech_poly_struct * alphabetas,
    const fq_zech_mpoly_ctx_t ctx)
{
    slong n = ctx->minfo->nvars;
    slong i, N = mpoly_words_per_exp_sp(A->bits, ctx->minfo);
    slong * offsets, * shifts;
    slong offset, shift;
    slong start, stop;
    ulong e, mask = (-UWORD(1)) >> (FLINT_BITS - A->bits);
    slong * starts, * ends, * stops;
    ulong * es;
    fq_zech_poly_struct * realE;

    E->length = 0;
    if (A->length < 1)
        return;

    starts = FLINT_ARRAY_ALLOC(n, slong);
    ends   = FLINT_ARRAY_ALLOC(n, slong);
    stops  = FLINT_ARRAY_ALLOC(n, slong);
    es     = FLINT_ARRAY_ALLOC(n, ulong);
    realE  = FLINT_ARRAY_ALLOC(n + 1, fq_zech_poly_struct);
    for (i = 0; i < n + 1; i++)
        fq_zech_poly_init(realE + i, ctx->fqctx);

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

    fq_zech_bpoly_fit_length(E, e + 1, ctx->fqctx);
    while (E->length <= e)
    {
        fq_zech_poly_zero(E->coeffs + E->length, ctx->fqctx);
        E->length++;
    }

    _fq_zech_mpoly_evaluate_rest_n_poly_fq(realE, starts, ends, stops, es,
                    A->coeffs + start, A->exps + N*start, stop - start, 1,
                                        alphabetas, offsets, shifts, N, mask,
                                                ctx->minfo->nvars, ctx->fqctx);
    fq_zech_poly_set(E->coeffs + e, realE + 0, ctx->fqctx);

    if (stop < A->length)
    {
        FLINT_ASSERT(e > (mask & (A->exps[N*stop + offset] >> shift)));
        e = (mask & (A->exps[N*stop + offset] >> shift));
        start = stop;
        goto next;
    }

    fq_zech_bpoly_normalise(E, ctx->fqctx);

    for (i = 0; i < n + 1; i++)
        fq_zech_poly_clear(realE + i, ctx->fqctx);
    flint_free(realE);
    flint_free(starts);
    flint_free(ends);
    flint_free(stops);
    flint_free(es);

    flint_free(offsets);
    flint_free(shifts);
}


/* A = B(gen(var), 0) */
void _fq_zech_mpoly_set_fq_zech_bpoly_var1_zero(
    fq_zech_mpoly_t A,
    flint_bitcnt_t Abits,
    const fq_zech_bpoly_t B,
    slong var,
    const fq_zech_mpoly_ctx_t ctx)
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

    fq_zech_mpoly_fit_length_set_bits(A, Alen, Abits, ctx);

    Alen = 0;
    for (i = Blen - 1; i >= 0; i--)
    {
        FLINT_ASSERT(Alen < A->alloc);
        fq_zech_bpoly_get_coeff(A->coeffs + Alen, B, i, 0, ctx->fqctx);
        if (fq_zech_is_zero(A->coeffs + Alen, ctx->fqctx))
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

int fq_zech_mpoly_factor_irred_smprime_wang(
    fq_zech_mpolyv_t fac,
    const fq_zech_mpoly_t A,
    const fq_zech_mpoly_factor_t lcAfac,
    const fq_zech_mpoly_t lcA,
    const fq_zech_mpoly_ctx_t ctx,
    flint_rand_t state)
{
    int success;
    int alphas_tries_remaining, alphabetas_tries_remaining, alphabetas_length;
    const slong n = ctx->minfo->nvars - 1;
    slong i, j, k, r;
    fq_zech_struct * alpha;
    fq_zech_poly_struct * alphabetas;
    fq_zech_mpoly_struct * Aevals;
    slong * degs, * degeval;
    fq_zech_mpolyv_t tfac;
    fq_zech_mpoly_t t, Acopy;
    fq_zech_mpoly_struct * newA;
    fq_zech_poly_t Abfc;
    fq_zech_bpoly_t Ab;
    fq_zech_tpoly_t Abfp;
    fq_zech_mpoly_t m, mpow;
    fq_zech_mpolyv_t new_lcs, lc_divs;
/*
flint_printf("fq_nmod_mpoly_factor_irred_smprime_wang called p = %wu\n", ctx->fqctx->modulus->mod.n);
flint_printf("     A: "); fq_nmod_mpoly_print_pretty(A, NULL, ctx); flint_printf("\n");
flint_printf("lcAfac: "); fq_nmod_mpoly_factor_print_pretty(lcAfac, NULL, ctx); flint_printf("\n");
flint_printf("   lcA: "); fq_nmod_mpoly_print_pretty(lcA, NULL, ctx); flint_printf("\n");
*/
    FLINT_ASSERT(n > 1);
    FLINT_ASSERT(A->length > 1);
    FLINT_ASSERT(fq_zech_is_one(A->coeffs + 0, ctx->fqctx));
    FLINT_ASSERT(A->bits <= FLINT_BITS);

    fq_zech_mpoly_init(Acopy, ctx);
    fq_zech_mpoly_init(m, ctx);
    fq_zech_mpoly_init(mpow, ctx);

    fq_zech_mpolyv_init(new_lcs, ctx);
    fq_zech_mpolyv_init(lc_divs, ctx);

    fq_zech_poly_init(Abfc, ctx->fqctx);
    fq_zech_tpoly_init(Abfp, ctx->fqctx);
    fq_zech_bpoly_init(Ab, ctx->fqctx);

    degs    = FLINT_ARRAY_ALLOC(n + 1, slong);
    degeval = FLINT_ARRAY_ALLOC(n + 1, slong);
	alpha   = FLINT_ARRAY_ALLOC(n, fq_zech_struct);
    alphabetas = FLINT_ARRAY_ALLOC(n, fq_zech_poly_struct);
    Aevals  = FLINT_ARRAY_ALLOC(n, fq_zech_mpoly_struct);
	for (i = 0; i < n; i++)
    {
        fq_zech_init(alpha + i, ctx->fqctx);
        fq_zech_poly_init(alphabetas + i, ctx->fqctx);
		fq_zech_mpoly_init(Aevals + i, ctx);
    }
    fq_zech_mpolyv_init(tfac, ctx);
	fq_zech_mpoly_init(t, ctx);

    /* init done */

    alphabetas_length = 2;
    alphas_tries_remaining = 10;
	fq_zech_mpoly_degrees_si(degs, A, ctx);

next_alpha:

    if (--alphas_tries_remaining < 0)
	{
		success = 0;
        goto cleanup;
	}

    for (i = 0; i < n; i++)
    {
        fq_zech_rand(alpha + i, state, ctx->fqctx);
    }

	for (i = n - 1; i >= 0; i--)
	{
        fq_zech_mpoly_evaluate_one_fq_zech(Aevals + i,
                       i == n - 1 ? A : Aevals + i + 1, i + 1, alpha + i, ctx);
		fq_zech_mpoly_degrees_si(degeval, Aevals + i, ctx);
		for (j = 0; j <= i; j++)
			if (degeval[j] != degs[j])
				goto next_alpha;
	}

    fq_zech_mpoly_get_fq_zech_poly(Abfc, Aevals + 0, 0, ctx);
    if (!fq_zech_poly_is_squarefree(Abfc, ctx->fqctx))
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
        fq_zech_poly_fit_length(alphabetas + i, alphabetas_length, ctx->fqctx);
        fq_zech_set(alphabetas[i].coeffs + 0, alpha + i, ctx->fqctx);
        for (j = 1; j < alphabetas_length; j++)
            fq_zech_rand(alphabetas[i].coeffs + j, state, ctx->fqctx);
        alphabetas[i].length = alphabetas_length;
        _fq_zech_poly_normalise(alphabetas + i, ctx->fqctx);
    }

    _fq_zech_eval_to_bpoly(Ab, A, alphabetas, ctx);
    success = fq_zech_bpoly_factor_smprime(Abfc, Abfp, Ab, 0, ctx->fqctx);
    if (!success)
    {
        FLINT_ASSERT(0 && "this should not happen");
        goto next_alpha;
    }

    r = Abfp->length;

    if (r < 2)
    {
        fq_zech_mpolyv_fit_length(fac, 1, ctx);
        fac->length = 1;
        fq_zech_mpoly_set(fac->coeffs + 0, A, ctx);
        success = 1;
        goto cleanup;
    }

    fq_zech_mpolyv_fit_length(lc_divs, r, ctx);
    lc_divs->length = r;
    if (lcAfac->num > 0)
    {
        success = fq_zech_mpoly_factor_lcc_wang(lc_divs->coeffs, lcAfac,
                                       Abfc, Abfp->coeffs, r, alphabetas, ctx);
        if (!success)
            goto next_alphabetas;
    }
    else
    {
        for (i = 0; i < r; i++)
            fq_zech_mpoly_one(lc_divs->coeffs + i, ctx);
    }

    success = fq_zech_mpoly_divides(m, lcA, lc_divs->coeffs + 0, ctx);
    FLINT_ASSERT(success);
    for (i = 1; i < r; i++)
    {
        success = fq_zech_mpoly_divides(m, m, lc_divs->coeffs + i, ctx);
        FLINT_ASSERT(success);
    }

    fq_zech_mpoly_pow_ui(mpow, m, r - 1, ctx);
    if (fq_zech_mpoly_is_one(mpow, ctx))
    {
        newA = (fq_zech_mpoly_struct *) A;
    }
    else
    {
        newA = Acopy;
        fq_zech_mpoly_mul(newA, A, mpow, ctx);
    }

    if (newA->bits > FLINT_BITS)
    {
        success = 0;
        goto cleanup;
    }

    fq_zech_mpoly_degrees_si(degs, newA, ctx);

    fq_zech_mpoly_set(t, mpow, ctx);
    for (i = n - 1; i >= 0; i--)
    {
        fq_zech_mpoly_evaluate_one_fq_zech(t, mpow, i + 1, alpha + i, ctx);
        fq_zech_mpoly_swap(t, mpow, ctx);
        fq_zech_mpoly_mul(Aevals + i, Aevals + i, mpow, ctx);
    }

    fq_zech_mpolyv_fit_length(new_lcs, (n + 1)*r, ctx);
    i = n;
    for (j = 0; j < r; j++)
    {
        fq_zech_mpoly_mul(new_lcs->coeffs + i*r + j, lc_divs->coeffs + j, m, ctx);
    }
    for (i = n - 1; i >= 0; i--)
    {
        for (j = 0; j < r; j++)
        {
            fq_zech_mpoly_evaluate_one_fq_zech(new_lcs->coeffs + i*r + j,
                       new_lcs->coeffs + (i + 1)*r + j, i + 1, alpha + i, ctx);
        }
    }

    fq_zech_mpolyv_fit_length(fac, r, ctx);
    fac->length = r;
    for (i = 0; i < r; i++)
    {
        fq_zech_t q;
        fq_zech_init(q, ctx->fqctx);
        FLINT_ASSERT(fq_zech_mpoly_is_fq_zech(new_lcs->coeffs + 0*r + i, ctx));
        FLINT_ASSERT(fq_zech_mpoly_length(new_lcs->coeffs + 0*r + i, ctx) == 1);
        _fq_zech_mpoly_set_fq_zech_bpoly_var1_zero(fac->coeffs + i, newA->bits, Abfp->coeffs + i, 0, ctx);
        FLINT_ASSERT(fac->coeffs[i].length > 0);
        fq_zech_inv(q, fac->coeffs[i].coeffs + 0, ctx->fqctx);
        fq_zech_mul(q, q, new_lcs->coeffs[0*r + i].coeffs + 0, ctx->fqctx);
        fq_zech_mpoly_scalar_mul_fq_zech(fac->coeffs + i, fac->coeffs + i, q, ctx);
        fq_zech_clear(q, ctx->fqctx);
    }

    fq_zech_mpolyv_fit_length(tfac, r, ctx);
    tfac->length = r;
    for (k = 1; k <= n; k++)
    {
        for (i = 0; i < r; i++)
        {
            _fq_zech_mpoly_set_lead0(tfac->coeffs + i, fac->coeffs + i,
                                               new_lcs->coeffs + k*r + i, ctx);
        }

        success = fq_zech_mpoly_hlift(k, tfac->coeffs, r, alpha,
                                         k < n ? Aevals + k : newA, degs, ctx);

        if (!success)
            goto next_alphabetas;

        fq_zech_mpolyv_swap(tfac, fac, ctx);
    }

    if (!fq_zech_mpoly_is_fq_zech(m, ctx))
    {
        success = 0;
        goto cleanup;
#if 0
        fq_zech_mpoly_univar_t u;
        fq_zech_mpoly_univar_init(u, ctx);
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
#endif
    }

    for (i = 0; i < r; i++)
        fq_zech_mpoly_make_monic(fac->coeffs + i, fac->coeffs + i, ctx);

    success = 1;

cleanup:

    fq_zech_mpolyv_clear(new_lcs, ctx);
    fq_zech_mpolyv_clear(lc_divs, ctx);

    fq_zech_bpoly_clear(Ab, ctx->fqctx);
    fq_zech_poly_clear(Abfc, ctx->fqctx);
    fq_zech_tpoly_clear(Abfp, ctx->fqctx);

	for (i = 0; i < n; i++)
    {
		fq_zech_mpoly_clear(Aevals + i, ctx);
        fq_zech_poly_clear(alphabetas + i, ctx->fqctx);
    }
    flint_free(alphabetas);
    flint_free(alpha);
    flint_free(Aevals);
    flint_free(degs);
    flint_free(degeval);

    fq_zech_mpolyv_clear(tfac, ctx);
    fq_zech_mpoly_clear(t, ctx);

    fq_zech_mpoly_clear(Acopy, ctx);
    fq_zech_mpoly_clear(m, ctx);
    fq_zech_mpoly_clear(mpow, ctx);

#if WANT_ASSERT
    if (success)
    {
        fq_zech_mpoly_t prod;
        fq_zech_mpoly_init(prod, ctx);
        fq_zech_mpoly_one(prod, ctx);
        for (i = 0; i < fac->length; i++)
            fq_zech_mpoly_mul(prod, prod, fac->coeffs + i, ctx);
        FLINT_ASSERT(fq_zech_mpoly_equal(prod, A, ctx));
        fq_zech_mpoly_clear(prod, ctx);
    }
#endif

flint_printf("fq_zech_mpoly_factor_irred_smprime_wang returning %d\n", success);
usleep(1000000);


	return success;
}
