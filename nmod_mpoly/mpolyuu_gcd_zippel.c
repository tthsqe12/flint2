/*
    Copyright (C) 2019 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "nmod_mpoly.h"
#include "nmod_mpoly_factor.h"


void nmod_mpolyuu_print_pretty(
    const nmod_mpolyu_t poly,
    const char ** x,
    slong nmainvars,
    const nmod_mpoly_ctx_t ctx)
{
    slong i, j;
    ulong mask = (-UWORD(1)) >> (FLINT_BITS - FLINT_BITS/nmainvars);

    if (poly->length == 0)
        flint_printf("0");

    for (i = 0; i < poly->length; i++)
    {
        if (i != 0)
            flint_printf(" + ");
        flint_printf("(");
        nmod_mpoly_print_pretty(poly->coeffs + i, x, ctx);
        flint_printf(")");
        for (j = nmainvars - 1; j >= 0; j--)
        {
            flint_printf("*X%wd^%wd", nmainvars - 1 - j,
                    mask & (poly->exps[i] >> (FLINT_BITS/nmainvars*j)));
        }
    }
}

/*
    for 0 <= i < mvars
        gen(i) -> alpha[i]
*/
static void nmod_mpolyu_set_eval_helper(
    n_polyun_t EH,
    const nmod_mpolyu_t A,
    const mp_limb_t * betas,
    slong mvars,
    const nmod_mpoly_ctx_t ctx)
{
    slong i, j, n;
    nmod_mpoly_struct * Acoeffs = A->coeffs;
    ulong * Aexps = A->exps;
    flint_bitcnt_t Abits = Abits;
    n_polyun_term_struct * EHterms;
    mp_limb_t * p, * q;

    n_polyun_fit_length(EH, A->length);
    EH->length = A->length;
    EHterms = EH->terms;

    for (i = 0; i < A->length; i++)
    {
        EHterms[i].exp = Aexps[i];
        n = Acoeffs[i].length;
        n_poly_fit_length(EHterms[i].coeff, 3*n);
        EHterms[i].coeff->length = n;
        p = EHterms[i].coeff->coeffs;
        _nmod_mpoly_monomial_evals(p, Acoeffs[i].exps, Abits, n, betas,
                                                                0, mvars, ctx);
        q = Acoeffs[i].coeffs;
        for (j = n - 1; j >= 0; j--)
        {
            mp_limb_t t1 = p[j];
            mp_limb_t t2 = q[j];
            p[3*j + 0] = t1;
            p[3*j + 1] = t2;
            p[3*j + 2] = t1;
        }
    }
}


void nmod_mpolyu_evaluate_one_ui(
    nmod_mpolyu_t E,
    nmod_mpolyu_t A,
    slong var,
    mp_limb_t alpha,
    const nmod_mpoly_ctx_t ctx)
{
    slong i, Elen;

    FLINT_ASSERT(E != A);
    FLINT_ASSERT(E->bits == A->bits);

    nmod_mpolyu_fit_length(E, A->length, ctx);
    Elen = 0;
    for (i = 0; i < A->length; i++)
    {
        E->exps[i] = A->exps[i];
        nmod_mpoly_evaluate_one_ui(E->coeffs + i, A->coeffs + i, var, alpha, ctx);
        nmod_mpoly_repack_bits_inplace(E->coeffs + i, E->bits, ctx);
        Elen += !nmod_mpoly_is_zero(E->coeffs + i, ctx);
    }

    E->length = Elen;
}

static void zip_start(n_polyun_t Z, n_polyun_t H, slong req_images)
{
    slong j;
    n_polyun_fit_length(Z, H->length);
    Z->length = H->length;
    for (j = 0; j < H->length; j++)
    {
        Z->terms[j].exp = H->terms[j].exp;
        n_poly_fit_length(Z->terms[j].coeff, req_images);
        Z->terms[j].coeff->length = 0;
    }
}


static void nmod_mpolyuu_get_n_bpoly(
    n_bpoly_t A,
    nmod_mpolyu_t B,
    const nmod_mpoly_ctx_t ctx)
{
    slong i;
    n_bpoly_zero(A);
    for (i = 0; i < B->length; i++)
    {
        slong x = extract_exp(B->exps[i], 0, 2);
        slong y = extract_exp(B->exps[i], 1, 2);
        n_bpoly_set_coeff(A, x, y, nmod_mpoly_get_ui(B->coeffs + i, ctx));
    }
}


static void nmod_mpolyuun_interp_lift_sm_bpoly(
    nmod_mpolyun_t An,
    n_bpoly_t A,
    const nmod_mpoly_ctx_t ctx)
{
    slong N = mpoly_words_per_exp_sp(An->bits, ctx->minfo);
    slong i, j;

    An->length = 0;

    for (i = A->length - 1; i >= 0; i--)
    {
        n_poly_struct * Ai = A->coeffs + i;
        for (j = Ai->length - 1; j >= 0; j--)
        {
            if (Ai->coeffs[j] == 0)
                continue;

            nmod_mpolyun_fit_length(An, An->length + 1, ctx);
            An->exps[An->length] = pack_exp2(i, j);
            nmod_mpolyn_fit_length(An->coeffs + 0, 1, ctx);
            mpoly_monomial_zero(An->coeffs[0].exps + N*0, N);
            nmod_poly_zero(An->coeffs[0].coeffs + 0);
            nmod_poly_set_coeff_ui(An->coeffs[0].coeffs + 0, 0, Ai->coeffs[j]);
            An->length++;
        }
    }
}

static int nmod_mpolyuun_interp_crt_sm_bpoly(
    slong * lastdeg,
    nmod_mpolyun_t F,
    nmod_mpolyun_t T,
    n_bpoly_t A,
    nmod_poly_t modulus,
    mp_limb_t alpha,
    const nmod_mpoly_ctx_t ctx)
{
    FLINT_ASSERT(0);
    return 0;
}    



int nmod_mpolyuu_gcd_zippel(
    nmod_mpolyu_t G,
    nmod_mpolyu_t Abar,
    nmod_mpolyu_t Bbar,
    nmod_mpolyu_t A,
    nmod_mpolyu_t B,
    const nmod_mpoly_t Gamma,
    const nmod_mpoly_ctx_t ctx)
{
    int success;
    nmod_t mod = ctx->ffinfo->mod;
    slong i, j, m;
    slong nvars = ctx->minfo->nvars;
    flint_bitcnt_t bits = A->bits;
    mp_limb_t * alphas, * betas;
    flint_rand_t state;
    n_polyun_t HG, HAbar, HBbar, MG, MAbar, MBbar, ZG, ZAbar, ZBbar;
    n_bpoly_t Aeval, Beval, Geval, Abareval, Bbareval;
    nmod_mpolyun_t Tn, Gn, Abarn, Bbarn;
    slong lastdegG, lastdegAbar, lastdegBbar;
    slong abar_images, bbar_images, g_images;
    n_polyun_t Aeh, Beh;
    n_poly_t Gammaeh;
    nmod_mpolyu_struct * Aevals, * Bevals;
    nmod_mpoly_struct * Gammaevals;
    slong * Adegs, * Bdegs, * Gammadegs;
    nmod_poly_t modulus;
    n_poly_bpoly_stack_t St;
    mp_limb_t c, start_alpha, gamma;

flint_printf("nmod_mpolyuu_gcd_zippel called nvars = %wd\n", nvars);
flint_printf("A: "); nmod_mpolyuu_print_pretty(A, NULL, 2, ctx); flint_printf("\n");
flint_printf("B: "); nmod_mpolyuu_print_pretty(B, NULL, 2, ctx); flint_printf("\n");
flint_printf("Gamma: "); nmod_mpoly_print_pretty(Gamma, NULL, ctx); flint_printf("\n");

    FLINT_ASSERT(bits == G->bits);
    FLINT_ASSERT(bits == Abar->bits);
    FLINT_ASSERT(bits == Bbar->bits);
    FLINT_ASSERT(bits == A->bits);
    FLINT_ASSERT(bits == B->bits);
    FLINT_ASSERT(bits == Gamma->bits);

    if (ctx->ffinfo->mod.n < 5)
        return 0;

    n_polyun_init(HG);
    n_polyun_init(HAbar);
    n_polyun_init(HBbar);
    n_polyun_init(MG);
    n_polyun_init(MAbar);
    n_polyun_init(MBbar);
    n_polyun_init(ZG);
    n_polyun_init(ZAbar);
    n_polyun_init(ZBbar);
    n_bpoly_init(Aeval);
    n_bpoly_init(Beval);
    n_bpoly_init(Geval);
    n_bpoly_init(Abareval);
    n_bpoly_init(Bbareval);
    nmod_mpolyun_init(Tn, bits, ctx);
    nmod_mpolyun_init(Gn, bits, ctx);
    nmod_mpolyun_init(Abarn, bits, ctx);
    nmod_mpolyun_init(Bbarn, bits, ctx);
    n_polyun_init(Aeh);
    n_polyun_init(Beh);
    n_poly_init(Gammaeh);
    nmod_poly_init_mod(modulus, mod);
    n_poly_stack_init(St->poly_stack);
    n_bpoly_stack_init(St->bpoly_stack);

    betas = FLINT_ARRAY_ALLOC(nvars, mp_limb_t);
    alphas = FLINT_ARRAY_ALLOC(nvars, mp_limb_t);
    flint_randinit(state);

    Aevals = FLINT_ARRAY_ALLOC(nvars + 1, nmod_mpolyu_struct);
    Bevals = FLINT_ARRAY_ALLOC(nvars + 1, nmod_mpolyu_struct);
    Gammaevals = FLINT_ARRAY_ALLOC(nvars + 1, nmod_mpoly_struct);
    for (i = 0; i < nvars; i++)
    {
        nmod_mpolyu_init(Aevals + i, bits, ctx);
        nmod_mpolyu_init(Bevals + i, bits, ctx);
        nmod_mpoly_init3(Gammaevals + i, 0, bits, ctx);
    }
    Aevals[nvars] = *A;
    Bevals[nvars] = *B;
    Gammaevals[nvars] = *Gamma;

    Gammadegs = FLINT_ARRAY_ALLOC(nvars, slong);
    Adegs = FLINT_ARRAY_ALLOC(nvars, slong);
    Bdegs = FLINT_ARRAY_ALLOC(nvars, slong);

    for (j = 0; j < nvars; j++)
        Adegs[j] = Bdegs[j] = 0;

    for (i = 0; i < A->length; i++)
    {
        nmod_mpoly_degrees_si(Gammadegs, A->coeffs + i, ctx);
        for (j = 0; j < nvars; j++)
            Adegs[j] = FLINT_MAX(Adegs[j], Gammadegs[j]);
    }

    for (i = 0; i < B->length; i++)
    {
        nmod_mpoly_degrees_si(Gammadegs, B->coeffs + i, ctx);
        for (j = 0; j < nvars; j++)
            Bdegs[j] = FLINT_MAX(Bdegs[j], Gammadegs[j]);
    }

    nmod_mpoly_degrees_si(Gammadegs, Gamma, ctx);

choose_alphas:

    for (i = 0; i < ctx->minfo->nvars; i++)
        alphas[i] = n_urandint(state, ctx->ffinfo->mod.n - 2) + 1;

    for (i = nvars - 1; i >= 0; i--)
    {
        nmod_mpolyu_evaluate_one_ui(Aevals + i, Aevals + i + 1, i, alphas[i], ctx);
        nmod_mpolyu_evaluate_one_ui(Bevals + i, Bevals + i + 1, i, alphas[i], ctx);
        nmod_mpoly_evaluate_one_ui(Gammaevals + i, Gammaevals + i, i, alphas[i], ctx);
        if (nmod_mpoly_is_zero(Gammaevals + i, ctx))
            goto choose_alphas;
        nmod_mpoly_repack_bits_inplace(Gammaevals + i, bits, ctx);
    }


    m = 0;

    nmod_mpolyuu_get_n_bpoly(Aeval, Aevals + m, ctx);
    nmod_mpolyuu_get_n_bpoly(Beval, Bevals + m, ctx);
    success = n_bpoly_mod_gcd_brown_smprime(Geval, Abareval, Bbareval, Aeval, Beval, mod, St);
    FLINT_ASSERT(success);
    n_bpoly_scalar_mul_nmod(Geval, nmod_mpoly_get_ui(Gammaevals + m, ctx), mod);

    nmod_mpolyuun_interp_lift_sm_bpoly(Gn, Geval, ctx);
    nmod_mpolyuun_interp_lift_sm_bpoly(Abarn, Abareval, ctx);
    nmod_mpolyuun_interp_lift_sm_bpoly(Bbarn, Bbareval, ctx);

    n_poly_one((n_poly_struct *)modulus);
    c = nmod_neg(alphas[m], mod);
    n_poly_mod_shift_left_scalar_addmul((n_poly_struct *)modulus, 1, c, mod);

    start_alpha = alphas[m];
    while (1)
    {
        alphas[m] = (alphas[m] < 1) ? mod.n - 1 : alphas[m] - 1;
        if (alphas[m] == start_alpha)
        {
            FLINT_ASSERT(0);
            success = 0;
            goto cleanup;
        }

        nmod_mpolyu_evaluate_one_ui(Aevals + m, Aevals + m + 1, m, alphas[m], ctx);
        nmod_mpolyu_evaluate_one_ui(Bevals + m, Bevals + m + 1, m, alphas[m], ctx);
        nmod_mpoly_evaluate_one_ui(Gammaevals + m, Gammaevals + m + 1, m, alphas[m], ctx);
        
        nmod_mpolyuu_get_n_bpoly(Aeval, Aevals + m, ctx);
        nmod_mpolyuu_get_n_bpoly(Beval, Bevals + m, ctx);
        gamma = nmod_mpoly_get_ui(Gammaevals + m, ctx);
        FLINT_ASSERT(gamma != 0);
        success = n_bpoly_mod_gcd_brown_smprime(Geval, Abareval, Bbareval, Aeval, Beval, mod, St);        
        FLINT_ASSERT(success);
        n_bpoly_scalar_mul_nmod(Geval, gamma, mod);

        c = n_poly_mod_evaluate_nmod((n_poly_struct *)modulus, alphas[m], mod);
        c = nmod_inv(c, mod);
        _n_poly_mod_scalar_mul_nmod((n_poly_struct *)modulus, (n_poly_struct *)modulus, c, mod);
        nmod_mpolyuun_interp_crt_sm_bpoly(&lastdegG, Gn, Tn, Geval, modulus, alphas[m], ctx);
        nmod_mpolyuun_interp_crt_sm_bpoly(&lastdegAbar, Abarn, Tn, Abareval, modulus, alphas[m], ctx);
        nmod_mpolyuun_interp_crt_sm_bpoly(&lastdegBbar, Bbarn, Tn, Bbareval, modulus, alphas[m], ctx);
        c = nmod_neg(alphas[m], mod);
        n_poly_mod_shift_left_scalar_addmul((n_poly_struct*)modulus, 1, c, mod);

        if (n_poly_degree((n_poly_struct*)modulus) > 1 + Gammadegs[m] + Adegs[m] &&
            n_poly_degree((n_poly_struct*)modulus) > 1 + Gammadegs[m] + Bdegs[m])
        {
            break;
        }
    }

    nmod_mpolyu_cvtfrom_mpolyun(G, Gn, m, ctx);
    nmod_mpolyu_cvtfrom_mpolyun(Abar, Abarn, m, ctx);
    nmod_mpolyu_cvtfrom_mpolyun(Bbar, Bbarn, m, ctx);

    for (m = 1; m < nvars; m++)
    {
        /* G, Abar, Bbar are in gen(0), ..., gen(m - 1) */

        nmod_mpolyun_interp_lift_sm_mpolyu(Gn, G, ctx);
        nmod_mpolyun_interp_lift_sm_mpolyu(Abarn, Abar, ctx);
        nmod_mpolyun_interp_lift_sm_mpolyu(Bbarn, Bbar, ctx);

        n_poly_one((n_poly_struct*)modulus);
        c = nmod_neg(alphas[m], mod);
        n_poly_mod_shift_left_scalar_addmul((n_poly_struct*)modulus, 1, c, mod);

    choose_betas:

        /* only beta[0], beta[1], ..., beta[m - 1] will be used */
        for (i = 0; i < ctx->minfo->nvars; i++)
            betas[i] = n_urandint(state, ctx->ffinfo->mod.n - 3) + 2;

        abar_images = nmod_mpolyu_set_zip_form(HAbar, MAbar, Abar, betas, m, ctxp);
        bbar_images = nmod_mpolyu_set_zip_form(HAbar, MAbar, Bbar, betas, m, ctxp);
        g_images    = nmod_mpolyu_set_zip_form(HG, MG, G, betas, m, ctxp);
        req_zip_images = FLINT_MAX(abar_images, bbar_images);
        req_zip_images = FLINT_MAX(req_zip_images, g_images);

        start_alpha = alphas[m];
        while (1)
        {
            alphas[m] = (alphas[m] < 1) ? mod.n - 1 : alphas[m] - 1;
            if (alphas[m] == start_alpha)
            {
                FLINT_ASSERT(0);
                success = 0;
                goto cleanup;
            }

            nmod_mpolyu_evaluate_one(Aevals + m, Aevals + m + 1, m, alpha[m], ctx);
            nmod_mpolyu_evaluate_one(Bevals + m, Bevals + m + 1, m, alpha[m], ctx);
            nmod_mpoly_evaluate_one(Gammaevals + m, Gammaevals + m + 1, m, alpha[m], ctx);

            nmod_mpolyu_set_eval_helper(Aeh, Aevals + m, ctx);
            nmod_mpolyu_set_eval_helper(Beh, Aevals + m, ctx);
            nmod_mpoly_set_eval_helper(Gammaeh, Gammaevals + m, ctx);

            zip_start(ZAbar, HAbar);
            zip_start(ZBbar, HBbar);
            zip_start(ZG, HG);
            for (cur_zip_image = 0; cur_zip_image < req_zip_images; cur_zip_image++)
            {
                n_polyu_eval_step_to_bpoly(Ab, Aeh, mod);
                n_polyu_eval_step_to_bpoly(Bb, Beh, mod);
                gamma = nmod_mpoly_get_ui(Gammaevals + m);
                FLINT_ASSERT(gamma != 0);
                success = n_bpoly_mod_gcd_brown_smprime(Gb, Abarb, Bbarb, Ab, Bb, mod, St);        
                FLINT_ASSERT(success);
                n_bpoly_scalar_mul_nmod(Gb, gamma);

                success = n_polyu2_add_zip_must_match(ZG, Gb, cur_zip_image);
                FLINT_ASSERT(success);
                success = n_polyu2_add_zip_must_match(ZAbar, Abarb, cur_zip_image);
                FLINT_ASSERT(success);
                success = n_polyu2_add_zip_must_match(ZBbar, Bbarb, cur_zip_image);
                FLINT_ASSERT(success);
            }

            success = zip_solve(G, ZG, HG, MG, ctx);
            FLINT_ASSERT(success);
            success = zip_solve(Abar, ZAbar, HAbar, MAbar, ctx);
            FLINT_ASSERT(success);
            success = zip_solve(Bbar, ZBbar, HBbar, MBbar, ctx);
            FLINT_ASSERT(success);

            FLINT_ASSERT(nmod_poly_degree(modulus) > 0);
            c = nmod_poly_evaluate_nmod(modulus, alphas[m]);
            c = n_invmod(temp, ctx->ffinfo->mod.n);
            nmod_poly_scalar_mul_nmod(modulus, modulus, c);

            changed = nmod_mpolyun_interp_crt_sm_mpolyu(&lastdegG, Gn, Tn,
                                                   G, modulus, alphas[m], ctx);
            changed = nmod_mpolyun_interp_crt_sm_mpolyu(&lastdegAbar, Gn, Tn,
                                                Abar, modulus, alphas[m], ctx);
            changed = nmod_mpolyun_interp_crt_sm_mpolyu(&lastdegBbar, Gn, Tn,
                                                Bbar, modulus, alphas[m], ctx);

            if (n_poly_degree(modulus) > 1 + Gammadegs[m] + Adegs[m] &&
                n_poly_degree(modulus) > 1 + Gammadegs[m] + Bdegs[m])
            {
                break;
            }
        }

        nmod_mpolyu_cvtfrom_mpolyun(G, Gn, m, ctx);
        nmod_mpolyu_cvtfrom_mpolyun(Abar, Abarn, m, ctx);
        nmod_mpolyu_cvtfrom_mpolyun(Bbar, Bbarn, m, ctx);
    }

    success = 1;

cleanup:
    
    


    FLINT_ASSERT(0);
    return success;
}

