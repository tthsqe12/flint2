/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "nmod_mpoly_factor.h"
#include "profiler.h"

/*
    return 1: success
           0: failed
          -1: exception (large exps)
*/
int nmod_mpolyl_gcd_hensel(
    nmod_mpoly_t G,
    nmod_mpoly_t Abar,
    nmod_mpoly_t Bbar,
    const nmod_mpoly_t A,
    const nmod_mpoly_t B,
    const nmod_mpoly_ctx_t ctx)
{
    int success;
    int alphas_tries_remaining;
    const slong n = ctx->minfo->nvars - 1;
    slong i, k;
    mp_limb_t * alphas;
    nmod_mpoly_struct * Aevals, * Bevals;
    slong * Adegs, * Bdegs;
    nmod_mpoly_t t, Acopy, Bcopy, gamma, g, abar, bbar;
    nmod_mpoly_struct * newA, * newB;
    nmod_mpoly_struct * Glcs, * Abarlcs, * Bbarlcs;
    nmod_mpoly_struct Afac[2], Bfac[2], Atfac[2], Btfac[2];
    mp_limb_t q;
    slong Adegx, Bdegx;
    flint_rand_t state;
timeit_t timer;

timeit_start(timer);

    FLINT_ASSERT(n > 0);
    FLINT_ASSERT(A->length > 0);
    FLINT_ASSERT(B->length > 0);
    FLINT_ASSERT(A->bits <= FLINT_BITS);
    FLINT_ASSERT(B->bits <= FLINT_BITS);

    flint_randinit(state);

    Adegs  = FLINT_ARRAY_ALLOC(n + 1, slong);
    Bdegs  = FLINT_ARRAY_ALLOC(n + 1, slong);

    Glcs    = FLINT_ARRAY_ALLOC(n + 1, nmod_mpoly_struct);
    Abarlcs = FLINT_ARRAY_ALLOC(n + 1, nmod_mpoly_struct);
    Bbarlcs = FLINT_ARRAY_ALLOC(n + 1, nmod_mpoly_struct);
    for (i = 0; i < n + 1; i++)
    {
        nmod_mpoly_init(Glcs + i, ctx);
        nmod_mpoly_init(Abarlcs + i, ctx);
        nmod_mpoly_init(Bbarlcs + i, ctx);
    }

	alphas  = FLINT_ARRAY_ALLOC(n, mp_limb_t);
    Aevals = FLINT_ARRAY_ALLOC(n, nmod_mpoly_struct);
    Bevals = FLINT_ARRAY_ALLOC(n, nmod_mpoly_struct);
	for (i = 0; i < n; i++)
    {
		nmod_mpoly_init(Aevals + i, ctx);
		nmod_mpoly_init(Bevals + i, ctx);
    }

    nmod_mpoly_init(t, ctx);
    nmod_mpoly_init(gamma, ctx);
    nmod_mpoly_init(Acopy, ctx);
    nmod_mpoly_init(Bcopy, ctx);
    nmod_mpoly_init(g, ctx);
    nmod_mpoly_init(abar, ctx);
    nmod_mpoly_init(bbar, ctx);

    for (i = 0; i < 2; i++)
    {
        nmod_mpoly_init(Afac + i, ctx);
        nmod_mpoly_init(Bfac + i, ctx);
        nmod_mpoly_init(Atfac + i, ctx);
        nmod_mpoly_init(Btfac + i, ctx);
    }

    nmod_mpoly_init(gamma, ctx);
    nmod_mpoly_one(gamma, ctx); /* haha */

    /* init done */

    alphas_tries_remaining = 10;

    for (i = 0; i < n; i++)
    {
        alphas[i] = 0;
    }

    goto got_alpha;

next_alpha:

flint_printf("next_alpha\n");

    if (--alphas_tries_remaining < 0)
	{
		success = 0;
        goto cleanup;
	}

    for (i = 0; i < n; i++)
    {
        alphas[i] = n_urandint(state, ctx->ffinfo->mod.n - 1) + 1;
    }

got_alpha:


    /* ensure deg_gen(0) do not drop under evaluation */
    Adegx = nmod_mpoly_degree_si(A, 0, ctx);
    Bdegx = nmod_mpoly_degree_si(B, 0, ctx);
	for (i = n - 1; i >= 0; i--)
	{
		nmod_mpoly_evaluate_one_ui(Aevals + i, i == n - 1 ? A : Aevals + i + 1, i + 1, alphas[i], ctx);
		nmod_mpoly_evaluate_one_ui(Bevals + i, i == n - 1 ? B : Bevals + i + 1, i + 1, alphas[i], ctx);
		if (Adegx != nmod_mpoly_degree_si(Aevals + i, 0, ctx) ||
            Bdegx != nmod_mpoly_degree_si(Bevals + i, 0, ctx))
        {
    		goto next_alpha;
        }
	}

    /* univariate gcd */
	if (!nmod_mpoly_gcd_cofactors(g, abar, bbar, Aevals + 0, Bevals + 0, ctx))
    {
        success = -1;
        goto cleanup;
    }

    if (nmod_mpoly_is_one(gamma, ctx))
    {
        newA = (nmod_mpoly_struct *) A;
        newB = (nmod_mpoly_struct *) B;
    }
    else
    {
        newA = Acopy;
        nmod_mpoly_mul(newA, A, gamma, ctx);
        newB = Bcopy;
        nmod_mpoly_mul(newB, B, gamma, ctx);
    }

    if (newA->bits > FLINT_BITS || newB->bits > FLINT_BITS)
    {
        success = -1;
        goto cleanup;
    }

    nmod_mpoly_degrees_si(Adegs, newA, ctx);
    nmod_mpoly_degrees_si(Bdegs, newB, ctx);

    nmod_mpoly_set(t, gamma, ctx);
    for (i = n - 1; i >= 0; i--)
    {
        nmod_mpoly_evaluate_one_ui(t, t, i + 1, alphas[i], ctx);
        nmod_mpoly_mul(Aevals + i, Aevals + i, t, ctx);
        nmod_mpoly_mul(Bevals + i, Bevals + i, t, ctx);
    }

    /* lcs have length n+1 */
    i = n;
    nmod_mpoly_set(Glcs + i, gamma, ctx);
    nmod_mpolyl_lead_coeff(Abarlcs + i, A, 1, ctx);
    nmod_mpolyl_lead_coeff(Bbarlcs + i, B, 1, ctx);
    for (i = n - 1; i >= 0; i--)
    {
        nmod_mpoly_evaluate_one_ui(Glcs + i, Glcs + i + 1, i + 1, alphas[i], ctx);
        nmod_mpoly_evaluate_one_ui(Abarlcs + i, Abarlcs + i + 1, i + 1, alphas[i], ctx);
        nmod_mpoly_evaluate_one_ui(Bbarlcs + i, Bbarlcs + i + 1, i + 1, alphas[i], ctx);

        /* evaluation could have killed gamma */
        if (nmod_mpoly_is_zero(Glcs + i, ctx) ||
            nmod_mpoly_is_zero(Abarlcs + i, ctx) ||
            nmod_mpoly_is_zero(Bbarlcs + i, ctx))
        {
            goto next_alpha;
        }
    }

    /* make the leading coefficients match Glcs[0], Abarlcs[0], Bbarlcs[0] */

    FLINT_ASSERT(nmod_mpoly_is_ui(Glcs + 0, ctx) && Glcs[0].length == 1);
    FLINT_ASSERT(nmod_mpoly_is_ui(Abarlcs + 0, ctx) && Abarlcs[0].length == 1);
    FLINT_ASSERT(nmod_mpoly_is_ui(Bbarlcs + 0, ctx) && Bbarlcs[0].length == 1);

    q = nmod_inv(g->coeffs[0], ctx->ffinfo->mod);
    q = nmod_mul(q, Glcs[0].coeffs[0], ctx->ffinfo->mod);
    nmod_mpoly_scalar_mul_nmod_invertible(Afac + 0, g, q, ctx);
    nmod_mpoly_scalar_mul_nmod_invertible(Bfac + 0, g, q, ctx);

    FLINT_ASSERT(abar->coeffs[0] == Abarlcs[0].coeffs[0]);
    FLINT_ASSERT(bbar->coeffs[0] == Bbarlcs[0].coeffs[0]);
    nmod_mpoly_set(Afac + 1, abar, ctx);
    nmod_mpoly_set(Bfac + 1, bbar, ctx);

    for (k = 1; k <= n; k++)
    {

flint_printf("starting k = %wd\n", k);
/*
flint_printf("before Afac[0]: ");
nmod_mpoly_print_pretty(Afac + 0, NULL, ctx);
flint_printf("\n");
flint_printf("before Afac[1]: ");
nmod_mpoly_print_pretty(Afac + 1, NULL, ctx);
flint_printf("\n");
flint_printf("Glcs[%wd]: ", k);
nmod_mpoly_print_pretty(Glcs + k, NULL, ctx);
flint_printf("\n");
*/
        _nmod_mpoly_set_lead0(Atfac + 0, Afac + 0, Glcs + k, ctx);
        _nmod_mpoly_set_lead0(Atfac + 1, Afac + 1, Abarlcs + k, ctx);
        success = nmod_mpoly_hlift(k, Atfac, 2, alphas,
                                        k < n ? Aevals + k : newA, Adegs, ctx);
        if (!success)
{
FLINT_ASSERT(0);
            goto next_alpha;

}

        nmod_mpoly_swap(Afac + 0, Atfac + 0, ctx);
        nmod_mpoly_swap(Afac + 1, Atfac + 1, ctx);
/*
flint_printf(" after Afac[0]: ");
nmod_mpoly_print_pretty(Afac + 0, NULL, ctx);
flint_printf("\n");
flint_printf(" after Afac[1]: ");
nmod_mpoly_print_pretty(Afac + 1, NULL, ctx);
flint_printf("\n");
flint_printf("before Bfac[0]: ");
nmod_mpoly_print_pretty(Bfac + 0, NULL, ctx);
flint_printf("\n");
flint_printf("before Bfac[1]: ");
nmod_mpoly_print_pretty(Bfac + 1, NULL, ctx);
flint_printf("\n");
*/

        _nmod_mpoly_set_lead0(Btfac + 0, Bfac + 0, Glcs + k, ctx);
        _nmod_mpoly_set_lead0(Btfac + 1, Bfac + 1, Bbarlcs + k, ctx);
        success = nmod_mpoly_hlift(k, Btfac, 2, alphas,
                                        k < n ? Bevals + k : newB, Bdegs, ctx);
        if (!success)
{
FLINT_ASSERT(0);
            goto next_alpha;
}

        nmod_mpoly_swap(Bfac + 0, Btfac + 0, ctx);
        nmod_mpoly_swap(Bfac + 1, Btfac + 1, ctx);
/*
flint_printf(" after Bfac[0]: ");
nmod_mpoly_print_pretty(Bfac + 0, NULL, ctx);
flint_printf("\n");
flint_printf(" after Bfac[1]: ");
nmod_mpoly_print_pretty(Bfac + 1, NULL, ctx);
flint_printf("\n");
*/
    }

timeit_stop(timer);

flint_printf("Afac[0] length: %wd\n", Afac[0].length);
flint_printf("Afac[1] length: %wd\n", Afac[1].length);

flint_printf("Bfac[0] length: %wd\n", Bfac[0].length);
flint_printf("Bfac[1] length: %wd\n", Bfac[1].length);

flint_printf("time: %wd\n", timer->wall);

FLINT_ASSERT(0);

    success = 1;

cleanup:

FLINT_ASSERT(0);

	return success;
}
