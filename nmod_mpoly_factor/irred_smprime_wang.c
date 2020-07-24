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

static void _eval_rest(
    n_poly_t E,
    const mp_limb_t * Acoeffs,
    const ulong * Aexps,
    slong Alen,
    slong var,
    const n_poly_struct * alphabetas,
    const slong * offsets,
    const slong * shifts,
    slong N,
    ulong mask,
    slong nvars,
    nmod_t ctx)
{
    slong offset = offsets[var];
    slong shift = shifts[var];
    slong start, stop;
    ulong e, next_e;
    n_poly_t T;

    n_poly_zero(E);

    if (Alen < 1)
        return;

    if (var >= nvars)
    {
        FLINT_ASSERT(Alen == 1);
        n_poly_set_ui(E, Acoeffs[0]);
        return;
    }

    n_poly_init(T);

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
    n_poly_mod_add(E, E, T, ctx);

    if (stop < Alen)
    {
        next_e = (mask & (Aexps[N*stop + offset] >> shift));
        FLINT_ASSERT(next_e < e);
        n_poly_mod_pow(T, alphabetas, e - next_e, ctx);
        n_poly_mod_mul(E, E, T, ctx);
        e = next_e;
        start = stop;
        goto next;
    }
    else
    {
        n_poly_mod_pow(T, alphabetas, e, ctx);
        n_poly_mod_mul(E, E, T, ctx);
    }

    n_poly_clear(T);
}

static void _eval_to_bpoly(
    n_bpoly_t E,
    const nmod_mpoly_t A,
    const n_poly_struct * alphabetas,
    const nmod_mpoly_ctx_t ctx)
{
    slong i, N = mpoly_words_per_exp_sp(A->bits, ctx->minfo);
    slong * offsets, * shifts;
    slong offset, shift;
    slong start, stop;
    ulong e, mask = (-UWORD(1)) >> (FLINT_BITS - A->bits);

    E->length = 0;
    if (A->length < 1)
        return;

    offsets = (slong *) flint_malloc(ctx->minfo->nvars*sizeof(slong));
    shifts = (slong *) flint_malloc(ctx->minfo->nvars*sizeof(slong));
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

    n_bpoly_fit_length(E, e + 1);
    while (E->length <= e)
    {
        n_poly_zero(E->coeffs + E->length);
        E->length++;
    }

    _eval_rest(E->coeffs + e, A->coeffs + start, A->exps + N*start, stop - start, 1,
                             alphabetas, offsets, shifts, N, mask,
                                          ctx->minfo->nvars, ctx->ffinfo->mod);

    if (stop < A->length)
    {
        FLINT_ASSERT(e > (mask & (A->exps[N*stop + offset] >> shift)));
        e = (mask & (A->exps[N*stop + offset] >> shift));
        start = stop;
        goto next;
    }

    n_bpoly_normalise(E);

    flint_free(offsets);
    flint_free(shifts);
}


int nmod_mpoly_factor_irred_smprime_wang(
    nmod_mpolyv_t fac,
    const nmod_mpoly_t A,
    const nmod_mpoly_factor_t lcAfac,
    const nmod_mpoly_t lcA,
    const nmod_mpoly_ctx_t ctx,
    flint_rand_t state)
{
    int success;
    int tries_remaining = 10;
    const slong n = ctx->minfo->nvars - 1;
    slong i, j, m, r;
	fmpz_t subset;
    mp_limb_t * alpha;
    nmod_mpoly_struct * Aevals;
    slong * deg, * degeval;
    nmod_mpolyv_t qfac, pfac, tfac, dfac;
    nmod_mpoly_t t, p, q;
    nmod_mpoly_univar_t u;
    n_poly_t c;
    n_bpoly_t B;
    n_tpoly_t F;
    n_poly_struct * alphabetas;
    slong alphabetas_length;

flint_printf("nmod_mpoly_factor_irred_smprime_wang called p = %wu\n", ctx->ffinfo->mod.n);
flint_printf("     A: "); nmod_mpoly_print_pretty(A, NULL, ctx); flint_printf("\n");
flint_printf("lcAfac: "); nmod_mpoly_factor_print_pretty(lcAfac, NULL, ctx); flint_printf("\n");
flint_printf("   lcA: "); nmod_mpoly_print_pretty(lcA, NULL, ctx); flint_printf("\n");

    FLINT_ASSERT(n > 1);
    FLINT_ASSERT(A->length > 1);
    FLINT_ASSERT(A->coeffs[0] == 1);
    FLINT_ASSERT(A->bits <= FLINT_BITS);

	fmpz_init(subset);
	alpha = (mp_limb_t *) flint_malloc(n*sizeof(mp_limb_t));
    Aevals    = (nmod_mpoly_struct *) flint_malloc(n*sizeof(nmod_mpoly_struct));
    deg     = (slong *) flint_malloc((n + 1)*sizeof(slong));
    degeval = (slong *) flint_malloc((n + 1)*sizeof(slong));
    alphabetas = (n_poly_struct *) flint_malloc(n*sizeof(n_poly_struct));
	for (i = 0; i < n; i++)
    {
		nmod_mpoly_init(Aevals + i, ctx);
        n_poly_init(alphabetas + i);
    }
    nmod_mpolyv_init(pfac, ctx);
    nmod_mpolyv_init(qfac, ctx);
    nmod_mpolyv_init(tfac, ctx);
    nmod_mpolyv_init(dfac, ctx);
	nmod_mpoly_init(t, ctx);
	nmod_mpoly_init(p, ctx);
	nmod_mpoly_init(q, ctx);
	nmod_mpoly_univar_init(u, ctx);
    n_poly_init(c);
    n_bpoly_init(B);
    n_tpoly_init(F);

    /* init done */

    alphabetas_length = 2;
	nmod_mpoly_degrees_si(deg, A, ctx);

next_alpha:

    if (--tries_remaining < 0)
	{
		success = 0;
        goto cleanup;
	}

    for (i = 0; i < n; i++)
        alpha[i] = n_urandint(state, ctx->ffinfo->mod.n - 1) + 1;

    /* ensure degrees do not drop under evaluation */
	for (i = n - 1; i >= 0; i--)
	{
		nmod_mpoly_evaluate_one_ui(Aevals + i,
                        i == n - 1 ? A : Aevals + i + 1, i + 1, alpha[i], ctx);
		nmod_mpoly_degrees_si(degeval, Aevals + i, ctx);
		for (j = 0; j <= i; j++)
			if (degeval[j] != deg[j])
				goto next_alpha;
	}

    /* make sure univar is squarefree */
	nmod_mpoly_derivative(t, Aevals + 0, 0, ctx);
	nmod_mpoly_gcd(t, t, Aevals + 0, ctx);
	if (!nmod_mpoly_is_one(t, ctx))
		goto next_alpha;

    for (i = 0; i < n; i++)
    {
        n_poly_fit_length(alphabetas + i, alphabetas_length);
        alphabetas[i].coeffs[0] = alpha[i];
        for (j = 1; j < alphabetas_length; j++)
            alphabetas[i].coeffs[j] = n_urandint(state, ctx->ffinfo->mod.n);
        alphabetas[i].length = alphabetas_length;
        _n_poly_normalise(alphabetas + i);
flint_printf("trying alhpabetas[%wd]: ", i);
n_poly_print_pretty(alphabetas + i, "Y");
flint_printf("\n");
    }

    _eval_to_bpoly(B, A, alphabetas, ctx);

flint_printf("B: "); n_bpoly_print_pretty(B, "X", "Y"); flint_printf("\n");

    success = n_bpoly_mod_factor_smprime(c, F, B, 0, ctx->ffinfo->mod);
    FLINT_ASSERT(success);

flint_printf("c: "); n_poly_print_pretty(c, "Y"); flint_printf("\n");
for (i = 0; i < F->length; i++)
{
flint_printf("F[%wd]: ", i);
n_bpoly_print_pretty(F->coeffs + i, "X", "Y");
flint_printf("\n");
}


FLINT_ASSERT(0);

cleanup:

    

	return success;
}

