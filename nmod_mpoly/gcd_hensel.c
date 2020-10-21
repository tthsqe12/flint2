/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "nmod_mpoly_factor.h"


/*
    return 1: success
           0: failed
          -1: exception (large exps)

    in F[x_1, ..., x_n][X]:

    cont_X(A) = cont_X(B) = 1
    compute gamma = gcd(lc_X(A), lc_X(B))

    try to find evaluation point x_i -> alpha_i
    a = A(x_i = alpha_i) in F[X]
    b = B(x_i = alpha_i) in F[X]

    compute univariate a = g*abar, b = g*bbar

    try to find mu1, mu2 in F with gcd(g, mu1*abar + mu2*bbar) = 1

    set H = mu1*A + mu2*B
    lift the univariate factorization g * (mu1*abar + mu2*bbar) against gamma*H
    imposing (gamma, lc_X H) as the leading coefficients

    remove content from the lift of g to get G and test divisibility to get
    Abar and Bbar
*/
int nmod_mpolyl_gcd_hensel(
    nmod_mpoly_t G,
    nmod_mpoly_t Abar,
    nmod_mpoly_t Bbar,
    const nmod_mpoly_t A,
    const nmod_mpoly_t B,
    const nmod_mpoly_ctx_t ctx)
{
    int success, alphas_tries_remaining, gamma_is_one;
    const slong n = ctx->minfo->nvars - 1;
    slong i, k;
    flint_bitcnt_t bits = A->bits;
    mp_limb_t * alphas;
    mp_limb_t q, mu;
    nmod_mpoly_struct * Aevals, * Bevals, * Hevals;
    nmod_mpoly_struct * H; /* points to A, B, or Hevals + n */
    nmod_mpoly_struct * Glcs, * Hlcs;
    nmod_mpoly_struct Hfac[2], Htfac[2];
    slong * Hdegs;
    slong Adegx, Bdegx;
    nmod_mpoly_t t1, t2, g, abar, bbar, hbar;
    flint_rand_t state;

    FLINT_ASSERT(n > 0);
    FLINT_ASSERT(A->length > 0);
    FLINT_ASSERT(B->length > 0);
    FLINT_ASSERT(bits <= FLINT_BITS);
    FLINT_ASSERT(A->bits == bits);
    FLINT_ASSERT(B->bits == bits);
    FLINT_ASSERT(ctx->minfo->ord == ORD_LEX);

    flint_randinit(state);

    Hdegs  = FLINT_ARRAY_ALLOC(n + 1, slong);

    Glcs   = FLINT_ARRAY_ALLOC(3*(n + 1), nmod_mpoly_struct);
    Hlcs   = Glcs + (n + 1);
    Hevals = Hlcs + (n + 1);
    for (i = 0; i < n + 1; i++)
    {
        nmod_mpoly_init(Glcs + i, ctx);
        nmod_mpoly_init(Hlcs + i, ctx);
        nmod_mpoly_init(Hevals + i, ctx);
    }

	alphas = FLINT_ARRAY_ALLOC(n, mp_limb_t);
    Aevals = FLINT_ARRAY_ALLOC(2*(n + 1), nmod_mpoly_struct);
    Bevals = Aevals + (n + 1);
	for (i = 0; i < n; i++)
    {
		nmod_mpoly_init(Aevals + i, ctx);
		nmod_mpoly_init(Bevals + i, ctx);
    }

    nmod_mpoly_init(t1, ctx);
    nmod_mpoly_init(t2, ctx);
    nmod_mpoly_init(g, ctx);
    nmod_mpoly_init(abar, ctx);
    nmod_mpoly_init(bbar, ctx);
    nmod_mpoly_init(hbar, ctx);

    nmod_mpoly_init(Hfac + 0, ctx);
    nmod_mpoly_init(Hfac + 1, ctx);
    nmod_mpoly_init(Htfac + 0, ctx);
    nmod_mpoly_init(Htfac + 1, ctx);

    /* init done */

    alphas_tries_remaining = 10;

    for (i = 0; i < n; i++)
    {
        alphas[i] = 0;
    }

    goto got_alpha;

next_alpha:

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

    /* ensure deg_X do not drop under evaluation */
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
	if (!nmod_mpoly_gcd_cofactors(g, abar, bbar, Aevals + 0, Bevals + 0, ctx) ||
        !nmod_mpoly_gcd(t1, g, abar, ctx) ||
        !nmod_mpoly_gcd(t2, g, bbar, ctx))
    {
        goto fail_hard;
    }

    /* set Hlcs[n], Glcs[n] (gamma), H, and Hevals */
    if (nmod_mpoly_is_one(t1, ctx))
    {
        mu = 0; /* mu = mu2/mu1 = 0 */

        nmod_mpoly_swap(hbar, abar, ctx);

        nmod_mpolyl_lead_coeff(Hlcs + n, A, 1, ctx);
        nmod_mpolyl_lead_coeff(t2, B, 1, ctx);
        if (!nmod_mpoly_gcd(Glcs + n, Hlcs + n, t2, ctx))
            goto fail_hard;

        H = (nmod_mpoly_struct *) A;

        gamma_is_one = nmod_mpoly_is_one(Glcs + n, ctx);
        if (gamma_is_one)
            for (i = 0; i < n; i++)
                nmod_mpoly_swap(Hevals + i, Aevals + i, ctx);
    }
    else if (nmod_mpoly_is_one(t2, ctx))
    {
        mu = ctx->ffinfo->mod.n; /* mu = mu2/mu1 = infinity */

        nmod_mpoly_swap(hbar, bbar, ctx);

        nmod_mpolyl_lead_coeff(Hlcs + n, B, 1, ctx);
        nmod_mpolyl_lead_coeff(t2, A, 1, ctx);
        if (!nmod_mpoly_gcd(Glcs + n, Hlcs + n, t2, ctx))
            goto fail_hard;

        H = (nmod_mpoly_struct *) B;

        gamma_is_one = nmod_mpoly_is_one(Glcs + n, ctx);
        if (gamma_is_one)
            for (i = 0; i < n; i++)
                nmod_mpoly_swap(Hevals + i, Bevals + i, ctx);
    }
    else
    {
        int mu_tries_remaining = 10;

    next_mu:

        if (--mu_tries_remaining < 0)
        {
            success = 0;
            goto cleanup;
        }

        mu = n_urandint(state, ctx->ffinfo->mod.n - 1) + 1;
        nmod_mpoly_scalar_addmul_ui(hbar, abar, bbar, mu, ctx);
        if (!nmod_mpoly_gcd(t1, hbar, g, ctx))
        {
            success = -1;
            goto cleanup;
        }

        if (!nmod_mpoly_is_one(t1, ctx))
            goto next_mu;

        nmod_mpolyl_lead_coeff(t1, A, 1, ctx);
        nmod_mpolyl_lead_coeff(t2, B, 1, ctx);
        if (!nmod_mpoly_gcd(Glcs + n, t1, t2, ctx))
        {
            success = -1;
            goto cleanup;
        }

        H = Hevals + n;
        nmod_mpoly_scalar_addmul_ui(H, A, B, mu, ctx);
        nmod_mpolyl_lead_coeff(Hlcs + n, H, 1, ctx);

        gamma_is_one = nmod_mpoly_is_one(Glcs + n, ctx);
        if (gamma_is_one)
            for (i = 0; i < n; i++)
                nmod_mpoly_scalar_addmul_ui(Hevals + i, Aevals + i,
                                                        Bevals + i, mu, ctx);
    }

    if (!gamma_is_one)
    {
        nmod_mpoly_mul(Hevals + n, H, Glcs + n, ctx);
        H = Hevals + n;
        for (i = n - 1; i >= 0; i--)
            nmod_mpoly_evaluate_one_ui(Hlcs + i, Hevals + i + 1, i + 1, alphas[i], ctx);
    }

    if (H->bits > FLINT_BITS &&
        !nmod_mpoly_repack_bits_inplace(H, FLINT_BITS, ctx))
    {
        goto fail_hard;
    }

    /* the evals should all fit in H->bits */
    for (i = 0; i < n; i++)
        nmod_mpoly_repack_bits_inplace(Hevals + i, H->bits, ctx);

    nmod_mpoly_degrees_si(Hdegs, H, ctx);

    /* computed evaluated leading coeffs */
    for (i = n - 1; i >= 0; i--)
    {
        nmod_mpoly_evaluate_one_ui(Glcs + i, Glcs + i + 1, i + 1, alphas[i], ctx);
        nmod_mpoly_evaluate_one_ui(Hlcs + i, Hlcs + i + 1, i + 1, alphas[i], ctx);
        /* evaluation could have killed gamma */
        if (nmod_mpoly_is_zero(Glcs + i, ctx) ||
            nmod_mpoly_is_zero(Hlcs + i, ctx))
        {
            goto next_alpha;
        }
    }

    /* make the leading coefficients match Glcs[0], Hlcs[0] */

    FLINT_ASSERT(nmod_mpoly_is_ui(Glcs + 0, ctx) && Glcs[0].length == 1);
    FLINT_ASSERT(nmod_mpoly_is_ui(Hlcs + 0, ctx) && Hlcs[0].length == 1);

    q = nmod_inv(g->coeffs[0], ctx->ffinfo->mod);
    q = nmod_mul(q, Glcs[0].coeffs[0], ctx->ffinfo->mod);
    nmod_mpoly_scalar_mul_nmod_invertible(Hfac + 0, g, q, ctx);

    q = nmod_inv(hbar->coeffs[0], ctx->ffinfo->mod);
    q = nmod_mul(q, Hlcs[0].coeffs[0], ctx->ffinfo->mod);
    nmod_mpoly_scalar_mul_nmod_invertible(Hfac + 1, hbar, q, ctx);

    for (k = 1; k <= n; k++)
    {
        _nmod_mpoly_set_lead0(Htfac + 0, Hfac + 0, Glcs + k, ctx);
        _nmod_mpoly_set_lead0(Htfac + 1, Hfac + 1, Hlcs + k, ctx);
        success = nmod_mpoly_hlift(k, Htfac, 2, alphas,
                                           k < n ? Hevals + k : H, Hdegs, ctx);
        if (!success)
            goto next_alpha;

        nmod_mpoly_swap(Hfac + 0, Htfac + 0, ctx);
        nmod_mpoly_swap(Hfac + 1, Htfac + 1, ctx);
    }

    if (!nmod_mpolyl_content(t1, Hfac + 0, 1, ctx))
        goto fail_hard;

    success = nmod_mpoly_divides(G, Hfac + 0, t1, ctx);
    FLINT_ASSERT(success);

    if (mu == 0)
    {
        /* the division by t1 should succeed, but let's be careful */
        nmod_mpolyl_lead_coeff(t1, G, 1, ctx);
        success = nmod_mpoly_divides(Abar, Hfac + 1, t1, ctx) &&
                  nmod_mpoly_divides(Bbar, B, G, ctx);
    }
    else if (mu == ctx->ffinfo->mod.n)
    {
        /* ditto */
        nmod_mpolyl_lead_coeff(t1, G, 1, ctx);
        success = nmod_mpoly_divides(Bbar, Hfac + 1, t1, ctx) &&
                  nmod_mpoly_divides(Abar, A, G, ctx);
    }
    else
    {
        success = nmod_mpoly_divides(Abar, A, G, ctx) &&
                  nmod_mpoly_divides(Bbar, B, G, ctx);
    }

    if (!success)
        goto next_alpha;

    nmod_mpoly_repack_bits_inplace(G, bits, ctx);
    nmod_mpoly_repack_bits_inplace(Abar, bits, ctx);
    nmod_mpoly_repack_bits_inplace(Bbar, bits, ctx);

    success = 1;

cleanup:

	return success;

fail_hard:

    success = -1;
    goto cleanup;
}

