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
#include "profiler.h"

void nmod_mpolyuun_print_pretty(
    const nmod_mpolyun_t poly,
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
        nmod_mpolyn_print_pretty(poly->coeffs + i, x, ctx);
        flint_printf(")");
        for (j = nmainvars - 1; j >= 0; j--)
        {
            flint_printf("*X%wd^%wd", nmainvars - 1 - j,
                    mask & (poly->exps[i] >> (FLINT_BITS/nmainvars*j)));
        }
    }
}

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

static void nmod_mpoly_set_eval_helper(
    n_poly_t EH,
    const nmod_mpoly_t A,
    const mp_limb_t * betas,
    slong mvars,
    const nmod_mpoly_ctx_t ctx)
{
    slong j, n;
    mp_limb_t * p, * q;

    n = A->length;
    n_poly_fit_length(EH, 3*n);
    EH->length = n;
    p = EH->coeffs;
    _nmod_mpoly_monomial_evals(p, A->exps, A->bits, n, betas, 0, mvars, ctx);
    q = A->coeffs;
    for (j = n - 1; j >= 0; j--)
    {
        mp_limb_t t1 = p[j];
        mp_limb_t t2 = q[j];
        p[3*j + 0] = t1;
        p[3*j + 1] = t2;
        p[3*j + 2] = t1;
    }
}


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
    flint_bitcnt_t Abits = A->bits;
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

static void nmod_mpolyuu_get_n_bpoly(
    n_bpoly_t A,
    nmod_mpolyu_t B,
    const nmod_mpoly_ctx_t ctx)
{
    slong i;
    n_bpoly_zero(A);
    for (i = 0; i < B->length; i++)
    {
        slong x = extract_exp(B->exps[i], 1, 2);
        slong y = extract_exp(B->exps[i], 0, 2);
        n_bpoly_set_coeff(A, x, y, nmod_mpoly_get_ui(B->coeffs + i, ctx));
    }
}


static void nmod_mpolyuun_interp_lift_sm_bpoly(
    nmod_mpolyun_t An,
    n_bpoly_t A,
    const nmod_mpoly_ctx_t ctx)
{
    slong N = mpoly_words_per_exp_sp(An->bits, ctx->minfo);
    slong i, j, Ani;

    Ani = 0;
    for (i = A->length - 1; i >= 0; i--)
    {
        n_poly_struct * Ai = A->coeffs + i;
        for (j = Ai->length - 1; j >= 0; j--)
        {
            if (Ai->coeffs[j] == 0)
                continue;

            nmod_mpolyun_fit_length(An, Ani + 1, ctx);
            An->exps[Ani] = pack_exp2(i, j);
            nmod_mpolyn_fit_length(An->coeffs + Ani, 1, ctx);
            An->coeffs[Ani].length = 1;
            mpoly_monomial_zero(An->coeffs[Ani].exps + N*0, N);
            nmod_poly_zero(An->coeffs[Ani].coeffs + 0);
            nmod_poly_set_coeff_ui(An->coeffs[Ani].coeffs + 0, 0, Ai->coeffs[j]);
            Ani++;
        }
    }

    An->length = Ani;
}

static nmod_poly_struct * getcf(nmod_mpolyn_t A)
{
    FLINT_ASSERT(A->length == 1);
    FLINT_ASSERT(A->exps[0] == 0);
    return A->coeffs + 0;
}


mp_limb_t _n_poly_eval_pow(n_poly_t P, n_poly_t alphapow, int nlimbs, nmod_t ctx)
{
    mp_limb_t * Pcoeffs = P->coeffs;
    slong Plen = P->length;
    mp_limb_t * alpha_powers = alphapow->coeffs;
    mp_limb_t res;
    slong k;

    if (Plen > alphapow->length)
    {
        slong oldlength = alphapow->length;
        FLINT_ASSERT(2 <= oldlength);
        n_poly_fit_length(alphapow, Plen);
        alphapow->length = Plen;
        alpha_powers = alphapow->coeffs;
        for (k = oldlength; k < Plen; k++)
            alpha_powers[k] = nmod_mul(alpha_powers[k - 1], alpha_powers[1], ctx);
    }

    NMOD_VEC_DOT(res, k, Plen, Pcoeffs[k], alpha_powers[k], ctx, nlimbs);

    return res;

}

static int nmod_mpolyuun_interp_crt_sm_bpoly(
    slong * lastdeg,
    nmod_mpolyun_t F,
    nmod_mpolyun_t T,
    n_bpoly_t A,
    nmod_poly_t modulus,
    n_poly_t alphapow,
    const nmod_mpoly_ctx_t ctx)
{
    int changed = 0;
    int nlimbs = _nmod_vec_dot_bound_limbs(modulus->length, ctx->ffinfo->mod);
    slong N = mpoly_words_per_exp(T->bits, ctx->minfo);
    n_poly_struct * Acoeffs = A->coeffs;
    slong Fi, Ti, Ai, ai;
    slong Flen = F->length;
    ulong * Fexps = F->exps;
    nmod_mpolyn_struct * Fcoeffs = F->coeffs;
    ulong * Texps = T->exps;
    nmod_mpolyn_struct * Tcoeffs = T->coeffs;
    mp_limb_t v;

    FLINT_ASSERT(T->bits == F->bits);

    *lastdeg = -1;

    Ti = Fi = 0;
    Ai = A->length - 1;
    ai = (Ai < 0) ? 0 : n_poly_degree(Acoeffs + Ai);

    while (Fi < Flen || Ai >= 0)
    {
        if (Ti >= T->alloc)
        {
            slong extra = Flen - Fi;
            extra = FLINT_MAX(extra, Ai);
            nmod_mpolyun_fit_length(T, Ti + extra + 1, ctx);
            Tcoeffs = T->coeffs;
            Texps = T->exps;
        }

        nmod_mpolyn_fit_length(Tcoeffs + Ti, 1, ctx);
        Tcoeffs[Ti].length = 1;
        mpoly_monomial_zero(Tcoeffs[Ti].exps, N);

        if (Fi < Flen)
        {
            FLINT_ASSERT((Fcoeffs + Fi)->length == 1);
            FLINT_ASSERT((Fcoeffs + Fi)->exps[0] == 0);
        }

        if (Fi < Flen && Ai >= 0 && Fexps[Fi] == pack_exp2(Ai, ai))
        {
            /* F term ok, A term ok */
            Texps[Ti] = pack_exp2(Ai, ai);

            v = _n_poly_eval_pow((n_poly_struct *)getcf(Fcoeffs + Fi), alphapow, nlimbs, ctx->ffinfo->mod);
            v = nmod_sub(Acoeffs[Ai].coeffs[ai], v, ctx->ffinfo->mod);
            if (v != 0)
            {
                changed = 1;
                n_poly_mod_scalar_addmul_nmod((n_poly_struct*)getcf(Tcoeffs + Ti), (n_poly_struct*)getcf(Fcoeffs + Fi), (n_poly_struct*)modulus, v, ctx->ffinfo->mod);
            }
            else
            {
                nmod_poly_set(getcf(Tcoeffs + Ti), getcf(Fcoeffs + Fi));
            }

            Fi++;
            do {
                ai--;
            } while (ai >= 0 && Acoeffs[Ai].coeffs[ai] == 0);
            if (ai < 0)
            {
                do {
                    Ai--;
                } while (Ai >= 0 && Acoeffs[Ai].length == 0);
                if (Ai >= 0)
                    ai = n_poly_degree(Acoeffs + Ai);
            }
        }
        else if (Ai >= 0 && (Fi >= Flen || Fexps[Fi] < pack_exp2(Ai, ai)))
        {
            /* F term missing, A term ok */

            Texps[Ti] = pack_exp2(Ai, ai);

            changed = 1;
            nmod_poly_scalar_mul_nmod(getcf(Tcoeffs + Ti), modulus, Acoeffs[Ai].coeffs[ai]);

            do {
                ai--;
            } while (ai >= 0 && Acoeffs[Ai].coeffs[ai] == 0);
            if (ai < 0)
            {
                do {
                    Ai--;
                } while (Ai >= 0 && Acoeffs[Ai].length == 0);
                if (Ai >= 0)
                    ai = n_poly_degree(Acoeffs + Ai);
            }
        }
        else
        {
            FLINT_ASSERT(Fi < Flen && (Ai < 0 || Fexps[Fi] > pack_exp2(Ai, ai)));
            /* F term ok, Aterm missing */

            Texps[Ti] = Fexps[Fi];

            v = _n_poly_eval_pow((n_poly_struct *)getcf(Fcoeffs + Fi), alphapow, nlimbs, ctx->ffinfo->mod);
            if (v != 0)
            {
                changed = 1;
                n_poly_mod_scalar_addmul_nmod((n_poly_struct*)getcf(Tcoeffs + Ti), (n_poly_struct*)getcf(Fcoeffs + Fi), (n_poly_struct*)modulus, v, ctx->ffinfo->mod);
            }
            else
            {
                nmod_poly_set(getcf(Tcoeffs + Ti), getcf(Fcoeffs + Fi));                
            }

            Fi++;
        }

        FLINT_ASSERT(!nmod_poly_is_zero(getcf(Tcoeffs + Ti)));
        *lastdeg = FLINT_MAX(*lastdeg, nmod_poly_degree(getcf(Tcoeffs + Ti)));
        Ti++;
    }

    T->length = Ti;

    if (changed)
        nmod_mpolyun_swap(T, F);

    return changed;
}





static slong nmod_mpolyu_set_zip_form(
    n_polyun_t H, /* monomial evals */
    n_polyun_t M, /* master polys */
    const nmod_mpolyu_t A,
    const mp_limb_t * betas,
    slong mvars,
    const nmod_mpoly_ctx_t ctx)
{
    slong i, n, zip_length = 0;
    n_polyun_term_struct * Hterms, * Mterms;

    n_polyun_fit_length(H, A->length);
    H->length = A->length;
    Hterms = H->terms;

    n_polyun_fit_length(M, A->length);
    M->length = A->length;
    Mterms = M->terms;

    for (i = 0; i < A->length; i++)
    {
        n = A->coeffs[i].length;
        zip_length = FLINT_MAX(zip_length, n);
        Hterms[i].exp = A->exps[i];
        Mterms[i].exp = A->exps[i];
        n_poly_fit_length(Hterms[i].coeff, n);
        Hterms[i].coeff->length = n;
        _nmod_mpoly_monomial_evals(Hterms[i].coeff->coeffs,
                          A->coeffs[i].exps, A->bits, n, betas, 0, mvars, ctx);
        n_poly_mod_product_roots_nmod_vec(Mterms[i].coeff,
                                 Hterms[i].coeff->coeffs, n, ctx->ffinfo->mod);
    }

    return zip_length;
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

static int zip_solve(
    nmod_mpolyu_t A,
    n_polyun_t Z,
    n_polyun_t H,
    n_polyun_t M,
    const nmod_mpoly_ctx_t ctx)
{
    int success;
    slong i, n;
    n_poly_t t;
    n_polyun_term_struct * Zterms = Z->terms;
    n_polyun_term_struct * Hterms = H->terms;
    n_polyun_term_struct * Mterms = M->terms;

    n_poly_init(t);

    FLINT_ASSERT(A->length == H->length);
    FLINT_ASSERT(A->length == M->length);
    FLINT_ASSERT(A->length == Z->length);

    for (i = 0; i < A->length; i++)
    {
        n = A->coeffs[i].length;
        FLINT_ASSERT(Mterms[i].coeff->length == n + 1);
        FLINT_ASSERT(Zterms[i].coeff->length >= n);
        FLINT_ASSERT(Hterms[i].coeff->length == n);

        n_poly_fit_length(t, n);

        success = nmod_zip_find_coeffs_new(A->coeffs[i].coeffs,
                         Hterms[i].coeff->coeffs, n,
                         Zterms[i].coeff->coeffs, Zterms[i].coeff->length,
                         Mterms[i].coeff->coeffs, t->coeffs, ctx->ffinfo->mod);
        if (success < 1)
        {
            n_poly_clear(t);
            return success;
        }
    }

    n_poly_clear(t);
    return 1;
}


#define USE_G    1
#define USE_ABAR 2
#define USE_BBAR 4

/*
    gamma = gcd(lc(A), lc(B))

    fake answers G, Abar, Bbar with

        G*Abar = gamma*A
        G*Bbar = gamma*B
        lc(G) = gamma
        lc(Abar) = lc(A)
        lc(Bbar) = lc(B)

    real answers

        rG = pp(G)
        rAbar = Abar/lc(G)
        rBbar = Bbar/lc(G)

    deg(G) + deg


*/
int nmod_mpolyuu_gcd_zippel(
    nmod_mpolyu_t rG,
    nmod_mpolyu_t rAbar,
    nmod_mpolyu_t rBbar,
    const nmod_mpolyu_t A,
    const nmod_mpolyu_t B,
    const nmod_mpoly_t gamma,
    const nmod_mpoly_ctx_t ctx)
{
    int success, changed, use = USE_G;
    nmod_t mod = ctx->ffinfo->mod;
    slong i, j, m;
    slong nvars = ctx->minfo->nvars;
    flint_bitcnt_t bits = A->bits;
    mp_limb_t * alphas, * betas;
    flint_rand_t state;
    nmod_mpoly_t cont;
    nmod_mpolyu_t T, G, Abar, Bbar;
    n_polyun_t HG, HAbar, HBbar, MG, MAbar, MBbar, ZG, ZAbar, ZBbar;
    n_bpoly_t Aev, Bev, Gev, Abarev, Bbarev;
    mp_limb_t gammaev;
    nmod_mpolyun_t Tn, Gn, Abarn, Bbarn;
    slong lastdeg, alpha_tries_remaining;
    slong cur_zip_image, req_zip_images, abar_images, bbar_images, g_images;
    n_polyun_t Aeh, Beh;
    n_poly_t gammaeh;
    n_poly_t alphapow;
    nmod_mpolyu_struct * Aevals, * Bevals;
    nmod_mpoly_struct * gammaevals;
    slong * Adegs, * Bdegs, * gammadegs;
    nmod_poly_t modulus;
    n_poly_bpoly_stack_t St;
    mp_limb_t c, start_alpha;
    ulong GdegboundXY, newdegXY;
/*
timeit_t timer;
timeit_start(timer);

flint_printf("nmod_mpolyuu_gcd_zippel called nvars = %wd\n", nvars);
flint_printf("mod.n = %wu\n", mod.n);
flint_printf("A: "); nmod_mpolyuu_print_pretty(A, NULL, 2, ctx); flint_printf("\n");
flint_printf("B: "); nmod_mpolyuu_print_pretty(B, NULL, 2, ctx); flint_printf("\n");
flint_printf("gamma: "); nmod_mpoly_print_pretty(gamma, NULL, ctx); flint_printf("\n");
*/
    FLINT_ASSERT(bits <= FLINT_BITS);
    FLINT_ASSERT(bits == rG->bits);
    FLINT_ASSERT(bits == rAbar->bits);
    FLINT_ASSERT(bits == rBbar->bits);
    FLINT_ASSERT(bits == A->bits);
    FLINT_ASSERT(bits == B->bits);
    FLINT_ASSERT(bits == gamma->bits);

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
    n_bpoly_init(Aev);
    n_bpoly_init(Bev);
    n_bpoly_init(Gev);
    n_bpoly_init(Abarev);
    n_bpoly_init(Bbarev);
    n_poly_init2(alphapow, 4);
    nmod_mpoly_init3(cont, 1, bits, ctx);
    nmod_mpolyu_init(T, bits, ctx);
    nmod_mpolyu_init(G, bits, ctx);
    nmod_mpolyu_init(Abar, bits, ctx);
    nmod_mpolyu_init(Bbar, bits, ctx);
    nmod_mpolyun_init(Tn, bits, ctx);
    nmod_mpolyun_init(Gn, bits, ctx);
    nmod_mpolyun_init(Abarn, bits, ctx);
    nmod_mpolyun_init(Bbarn, bits, ctx);
    n_polyun_init(Aeh);
    n_polyun_init(Beh);
    n_poly_init(gammaeh);
    nmod_poly_init_mod(modulus, mod);
    n_poly_stack_init(St->poly_stack);
    n_bpoly_stack_init(St->bpoly_stack);

    betas = FLINT_ARRAY_ALLOC(nvars, mp_limb_t);
    alphas = FLINT_ARRAY_ALLOC(nvars, mp_limb_t);
    flint_randinit(state);

    Aevals = FLINT_ARRAY_ALLOC(nvars + 1, nmod_mpolyu_struct);
    Bevals = FLINT_ARRAY_ALLOC(nvars + 1, nmod_mpolyu_struct);
    gammaevals = FLINT_ARRAY_ALLOC(nvars + 1, nmod_mpoly_struct);
    for (i = 0; i < nvars; i++)
    {
        nmod_mpolyu_init(Aevals + i, bits, ctx);
        nmod_mpolyu_init(Bevals + i, bits, ctx);
        nmod_mpoly_init3(gammaevals + i, 0, bits, ctx);
    }
    Aevals[nvars] = *A;
    Bevals[nvars] = *B;
    gammaevals[nvars] = *gamma;

    gammadegs = FLINT_ARRAY_ALLOC(nvars, slong);
    Adegs = FLINT_ARRAY_ALLOC(nvars, slong);
    Bdegs = FLINT_ARRAY_ALLOC(nvars, slong);

    for (j = 0; j < nvars; j++)
        Adegs[j] = Bdegs[j] = 0;

    for (i = 0; i < A->length; i++)
    {
        nmod_mpoly_degrees_si(gammadegs, A->coeffs + i, ctx);
        for (j = 0; j < nvars; j++)
            Adegs[j] = FLINT_MAX(Adegs[j], gammadegs[j]);
    }

    for (i = 0; i < B->length; i++)
    {
        nmod_mpoly_degrees_si(gammadegs, B->coeffs + i, ctx);
        for (j = 0; j < nvars; j++)
            Bdegs[j] = FLINT_MAX(Bdegs[j], gammadegs[j]);
    }

    nmod_mpoly_degrees_si(gammadegs, gamma, ctx);

    GdegboundXY = FLINT_MIN(A->exps[0], B->exps[0]);
    if (GdegboundXY == 0)
        goto gcd_is_trivial;

    alpha_tries_remaining = 100;

choose_alphas:

    if (--alpha_tries_remaining < 0)
    {
        success = 0;
        goto cleanup;
    }

    for (i = 0; i < ctx->minfo->nvars; i++)
        alphas[i] = n_urandint(state, ctx->ffinfo->mod.n - 2) + 1;

    for (i = nvars - 1; i >= 0; i--)
    {
        nmod_mpolyu_evaluate_one_ui(Aevals + i, Aevals + i + 1, i, alphas[i], ctx);
        nmod_mpolyu_evaluate_one_ui(Bevals + i, Bevals + i + 1, i, alphas[i], ctx);
        nmod_mpoly_evaluate_one_ui(gammaevals + i, gammaevals + i + 1, i, alphas[i], ctx);
        if (nmod_mpoly_is_zero(gammaevals + i, ctx))
            goto choose_alphas;
        nmod_mpoly_repack_bits_inplace(gammaevals + i, bits, ctx);
    }

    m = 0;

    use = USE_G | USE_ABAR | USE_BBAR;
    use = USE_G;

    FLINT_ASSERT(alphapow->alloc > 1);
    alphapow->coeffs[0] = 1;
    alphapow->coeffs[1] = alphas[m];
    alphapow->length = 2;

    nmod_mpolyuu_get_n_bpoly(Aev, Aevals + m, ctx);
    nmod_mpolyuu_get_n_bpoly(Bev, Bevals + m, ctx);
    success = n_bpoly_mod_gcd_brown_smprime(Gev, Abarev, Bbarev, Aev, Bev, mod, St);
    if (!success)
        goto cleanup;

    newdegXY = n_bpoly_bidegree(Gev);
    if (newdegXY > GdegboundXY)
        goto choose_alphas;
    if (newdegXY < GdegboundXY)
    {
        GdegboundXY = newdegXY;
        if (GdegboundXY == 0)
            goto gcd_is_trivial;
    }

    gammaev = nmod_mpoly_get_ui(gammaevals + m, ctx);
    n_bpoly_scalar_mul_nmod(Gev, gammaev, mod);

    nmod_mpolyuun_interp_lift_sm_bpoly(Gn, Gev, ctx);
    nmod_mpolyuun_interp_lift_sm_bpoly(Abarn, Abarev, ctx);
    nmod_mpolyuun_interp_lift_sm_bpoly(Bbarn, Bbarev, ctx);

    n_poly_one((n_poly_struct *)modulus);
    c = nmod_neg(alphas[m], mod);
    n_poly_mod_shift_left_scalar_addmul((n_poly_struct *)modulus, 1, c, mod);

    start_alpha = alphas[m];
    while (1)
    {
choose_alpha_0:

        alphas[m] = (alphas[m] < 1) ? mod.n - 1 : alphas[m] - 1;
        if (alphas[m] == start_alpha)
            goto choose_alphas;

        FLINT_ASSERT(alphapow->alloc > 1);
        alphapow->coeffs[0] = 1;
        alphapow->coeffs[1] = alphas[m];
        alphapow->length = 2;

        nmod_mpolyu_evaluate_one_ui(Aevals + m, Aevals + m + 1, m, alphas[m], ctx);
        nmod_mpolyu_evaluate_one_ui(Bevals + m, Bevals + m + 1, m, alphas[m], ctx);
        nmod_mpoly_evaluate_one_ui(gammaevals + m, gammaevals + m + 1, m, alphas[m], ctx);

        nmod_mpolyuu_get_n_bpoly(Aev, Aevals + m, ctx);
        nmod_mpolyuu_get_n_bpoly(Bev, Bevals + m, ctx);
        gammaev = nmod_mpoly_get_ui(gammaevals + m, ctx);
        if (gammaev == 0)
            goto choose_alpha_0;

        success = n_bpoly_mod_gcd_brown_smprime(Gev, Abarev, Bbarev, Aev, Bev, mod, St);
        if (!success)
            goto cleanup;

        newdegXY = n_bpoly_bidegree(Gev);
        if (newdegXY > GdegboundXY)
            goto choose_alpha_0;
        if (newdegXY < GdegboundXY)
        {
            GdegboundXY = newdegXY;
            if (GdegboundXY == 0)
                goto gcd_is_trivial;
            goto choose_alphas;
        }

        n_bpoly_scalar_mul_nmod(Gev, gammaev, mod);

        c = n_poly_mod_evaluate_nmod((n_poly_struct *)modulus, alphas[m], mod);
        c = nmod_inv(c, mod);
        _n_poly_mod_scalar_mul_nmod((n_poly_struct *)modulus, (n_poly_struct *)modulus, c, mod);

        if ((use & USE_G) && !nmod_mpolyuun_interp_crt_sm_bpoly(
                                &lastdeg, Gn, Tn, Gev, modulus, alphapow, ctx))
        {
            if (m == nvars - 1)
            {
                nmod_mpolyu_cvtfrom_mpolyun(rG, Gn, m, ctx);
                success = nmod_mpolyu_content_mpoly_threaded_pool(cont, rG, ctx, NULL, 0);
                if (!success)
                    goto cleanup;
                nmod_mpolyu_divexact_mpoly_inplace(rG, cont, ctx);
                if (nmod_mpolyuu_divides(rAbar, A, rG, 2, ctx) &&
                    nmod_mpolyuu_divides(rBbar, B, rG, 2, ctx))
                {
                    break;
                }                
            }
            else
            {
                nmod_mpolyu_cvtfrom_mpolyun(G, Gn, m, ctx);
                nmod_mpolyu_mul_mpoly(T, Aevals + m + 1, gammaevals + m + 1, ctx);
                if (nmod_mpolyuu_divides(Abar, T, G, 2, ctx))
                {
                    nmod_mpolyu_mul_mpoly(T, Bevals + m + 1, gammaevals + m + 1, ctx);
                    if (nmod_mpolyuu_divides(Bbar, T, G, 2, ctx))
                        break;
                }
            }
        }

        if ((use & USE_ABAR) && !nmod_mpolyuun_interp_crt_sm_bpoly(
                          &lastdeg, Abarn, Tn, Abarev, modulus, alphapow, ctx))
        {
            if (m == nvars - 1)
            {
                nmod_mpolyu_cvtfrom_mpolyun(rAbar, Abarn, m, ctx);
                success = nmod_mpolyu_content_mpoly_threaded_pool(cont, rAbar, ctx, NULL, 0);
                if (!success)
                    goto cleanup;
                nmod_mpolyu_divexact_mpoly_inplace(rAbar, cont, ctx);
                if (nmod_mpolyuu_divides(rG, A, rAbar, 2, ctx) &&
                    nmod_mpolyuu_divides(rBbar, B, rG, 2, ctx))
                {
                    break;
                }                
            }
            else
            {
                nmod_mpolyu_cvtfrom_mpolyun(Abar, Abarn, m, ctx);
                nmod_mpolyu_mul_mpoly(T, Aevals + m + 1, gammaevals + m + 1, ctx);
                if (nmod_mpolyuu_divides(G, T, Abar, 2, ctx))
                {
                    nmod_mpolyu_mul_mpoly(T, Bevals + m + 1, gammaevals + m + 1, ctx);
                    if (nmod_mpolyuu_divides(Bbar, T, G, 2, ctx))
                        break;
                }
            }
        }

        if ((use & USE_BBAR) && !nmod_mpolyuun_interp_crt_sm_bpoly(
                          &lastdeg, Bbarn, Tn, Bbarev, modulus, alphapow, ctx))
        {
            if (m == nvars - 1)
            {
                nmod_mpolyu_cvtfrom_mpolyun(rBbar, Bbarn, m, ctx);
                success = nmod_mpolyu_content_mpoly_threaded_pool(cont, rBbar, ctx, NULL, 0);
                if (!success)
                    goto cleanup;
                nmod_mpolyu_divexact_mpoly_inplace(rBbar, cont, ctx);
                if (nmod_mpolyuu_divides(rG, A, rBbar, 2, ctx) &&
                    nmod_mpolyuu_divides(rAbar, B, rG, 2, ctx))
                {
                    break;
                }                
            }
            else
            {
                nmod_mpolyu_cvtfrom_mpolyun(Bbar, Bbarn, m, ctx);
                nmod_mpolyu_mul_mpoly(T, Bevals + m + 1, gammaevals + m + 1, ctx);
                if (nmod_mpolyuu_divides(G, T, Bbar, 2, ctx))
                {
                    nmod_mpolyu_mul_mpoly(T, Aevals + m + 1, gammaevals + m + 1, ctx);
                    if (nmod_mpolyuu_divides(Abar, T, G, 2, ctx))
                        break;
                }
            }
        }

        c = nmod_neg(alphas[m], mod);
        n_poly_mod_shift_left_scalar_addmul((n_poly_struct*)modulus, 1, c, mod);

        if (nmod_poly_degree(modulus) - 1 > gammadegs[m] + Adegs[m] &&
            nmod_poly_degree(modulus) - 1 > gammadegs[m] + Bdegs[m])
        {
            goto choose_alphas;
        }
    }

    for (m = 1; m < nvars; m++)
    {
        /* G, Abar, Bbar are in Fp[gen(0), ..., gen(m - 1)] */

        nmod_mpolyun_interp_lift_sm_mpolyu(Gn, G, ctx);
        nmod_mpolyun_interp_lift_sm_mpolyu(Abarn, Abar, ctx);
        nmod_mpolyun_interp_lift_sm_mpolyu(Bbarn, Bbar, ctx);

        n_poly_one((n_poly_struct*)modulus);
        c = nmod_neg(alphas[m], mod);
        n_poly_mod_shift_left_scalar_addmul((n_poly_struct*)modulus, 1, c, mod);

        start_alpha = alphas[m];

        use = USE_G | USE_ABAR | USE_BBAR;
        use = USE_G;

    choose_betas:

        /* only beta[0], beta[1], ..., beta[m - 1] will be used */
        for (i = 0; i < ctx->minfo->nvars; i++)
            betas[i] = n_urandint(state, ctx->ffinfo->mod.n - 3) + 2;

        abar_images = nmod_mpolyu_set_zip_form(HAbar, MAbar, Abar, betas, m, ctx);
        bbar_images = nmod_mpolyu_set_zip_form(HBbar, MBbar, Bbar, betas, m, ctx);
        g_images    = nmod_mpolyu_set_zip_form(HG, MG, G, betas, m, ctx);
        req_zip_images = FLINT_MAX(abar_images, bbar_images);
        req_zip_images = FLINT_MAX(req_zip_images, g_images);

        while (1)
        {
choose_alpha_m:

            alphas[m] = (alphas[m] < 1) ? mod.n - 1 : alphas[m] - 1;
            if (alphas[m] == start_alpha)
                goto choose_alphas;

            nmod_mpolyu_evaluate_one_ui(Aevals + m, Aevals + m + 1, m, alphas[m], ctx);
            nmod_mpolyu_evaluate_one_ui(Bevals + m, Bevals + m + 1, m, alphas[m], ctx);
            nmod_mpoly_evaluate_one_ui(gammaevals + m, gammaevals + m + 1, m, alphas[m], ctx);
            if (nmod_mpoly_is_zero(gammaevals, ctx))
                goto choose_alpha_m;

            nmod_mpolyu_set_eval_helper(Aeh, Aevals + m, betas, m, ctx);
            nmod_mpolyu_set_eval_helper(Beh, Bevals + m, betas, m, ctx);
            nmod_mpoly_set_eval_helper(gammaeh, gammaevals + m, betas, m, ctx);

            zip_start(ZG, HG, req_zip_images);
            zip_start(ZAbar, HAbar, req_zip_images);
            zip_start(ZBbar, HBbar, req_zip_images);
            for (cur_zip_image = 0; cur_zip_image < req_zip_images; cur_zip_image++)
            {
                n_bpoly_mod_eval_step(Aev, Aeh, mod);
                n_bpoly_mod_eval_step(Bev, Beh, mod);
                gammaev = n_poly_mod_eval_step(gammaeh, mod);
                if (gammaev == 0)
                    goto choose_betas;

                success = n_bpoly_mod_gcd_brown_smprime(Gev, Abarev, Bbarev, Aev, Bev, mod, St);        
                if (!success)
                    goto cleanup;

                newdegXY = n_bpoly_bidegree(Gev);
                if (newdegXY > GdegboundXY)
                    goto choose_betas;
                if (newdegXY < GdegboundXY)
                {
                    GdegboundXY = newdegXY;
                    if (GdegboundXY == 0)
                        goto gcd_is_trivial;
                    goto choose_alphas;
                }

                n_bpoly_scalar_mul_nmod(Gev, gammaev, mod);
                if ((use & USE_G) && !n_polyu2_add_zip_must_match(ZG, Gev, cur_zip_image))
                    goto choose_alphas;
                if ((use & USE_ABAR) && !n_polyu2_add_zip_must_match(ZAbar, Abarev, cur_zip_image))
                    goto choose_alphas;
                if ((use & USE_BBAR) && !n_polyu2_add_zip_must_match(ZBbar, Bbarev, cur_zip_image))
                    goto choose_alphas;
            }

            if ((use & USE_G) && zip_solve(G, ZG, HG, MG, ctx) < 1)
                goto choose_alphas;
            if ((use & USE_ABAR) && zip_solve(Abar, ZAbar, HAbar, MAbar, ctx) < 1)
                goto choose_alphas;
            if ((use & USE_BBAR) && zip_solve(Bbar, ZBbar, HBbar, MBbar, ctx) < 1)
                goto choose_alphas;

            FLINT_ASSERT(nmod_poly_degree(modulus) > 0);
            c = nmod_poly_evaluate_nmod(modulus, alphas[m]);
            c = nmod_inv(c, mod);
            nmod_poly_scalar_mul_nmod(modulus, modulus, c);

            if ((use & USE_G) && !nmod_mpolyun_interp_crt_sm_mpolyu(
                                 &lastdeg, Gn, Tn, G, modulus, alphas[m], ctx))
            {
                nmod_mpolyu_cvtfrom_mpolyun(rG, Gn, m, ctx);
                if (m == nvars - 1)
                {
                    success = nmod_mpolyu_content_mpoly_threaded_pool(cont, rG, ctx, NULL, 0);
                    if (!success)
                        goto cleanup;
                    nmod_mpolyu_divexact_mpoly_inplace(rG, cont, ctx);
                    if (nmod_mpolyuu_divides(rAbar, A, rG, 2, ctx) &&
                        nmod_mpolyuu_divides(rBbar, B, rG, 2, ctx))
                    {
                        break;
                    }
                }
                else
                {
                    nmod_mpolyu_mul_mpoly(T, Aevals + m + 1, gammaevals + m + 1, ctx);
                    if (nmod_mpolyuu_divides(rAbar, T, rG, 2, ctx))
                    {
                        nmod_mpolyu_mul_mpoly(T, Bevals + m + 1, gammaevals + m + 1, ctx);
                        if (nmod_mpolyuu_divides(rBbar, T, rG, 2, ctx))
                        {
                            nmod_mpolyu_swap(G, rG, ctx);
                            nmod_mpolyu_swap(Abar, rAbar, ctx);
                            nmod_mpolyu_swap(Bbar, rBbar, ctx);
                            break;
                        }
                    }
                }
            }
            if ((use & USE_ABAR) && !nmod_mpolyun_interp_crt_sm_mpolyu(
                           &lastdeg, Abarn, Tn, Abar, modulus, alphas[m], ctx))
            {
                nmod_mpolyu_cvtfrom_mpolyun(rAbar, Abarn, m, ctx);
                if (m == nvars - 1)
                {
                    success = nmod_mpolyu_content_mpoly_threaded_pool(cont, rAbar, ctx, NULL, 0);
                    if (!success)
                        goto cleanup;
                    nmod_mpolyu_divexact_mpoly_inplace(rAbar, cont, ctx);
                    if (nmod_mpolyuu_divides(rG, A, rAbar, 2, ctx) &&
                        nmod_mpolyuu_divides(rBbar, B, rG, 2, ctx))
                    {
                        break;
                    }
                }
                else
                {
                    nmod_mpolyu_mul_mpoly(T, Aevals + m + 1, gammaevals + m + 1, ctx);
                    if (nmod_mpolyuu_divides(rG, T, rAbar, 2, ctx))
                    {
                        nmod_mpolyu_mul_mpoly(T, Bevals + m + 1, gammaevals + m + 1, ctx);
                        if (nmod_mpolyuu_divides(rBbar, T, rG, 2, ctx))
                        {
                            nmod_mpolyu_swap(G, rG, ctx);
                            nmod_mpolyu_swap(Abar, rAbar, ctx);
                            nmod_mpolyu_swap(Bbar, rBbar, ctx);
                            break;
                        }
                    }
                }
            }
            if ((use & USE_BBAR) && !nmod_mpolyun_interp_crt_sm_mpolyu(
                           &lastdeg, Bbarn, Tn, Bbar, modulus, alphas[m], ctx))
            {
                nmod_mpolyu_cvtfrom_mpolyun(rBbar, Bbarn, m, ctx);
                if (m == nvars - 1)
                {
                    success = nmod_mpolyu_content_mpoly_threaded_pool(cont, rBbar, ctx, NULL, 0);
                    if (!success)
                        goto cleanup;
                    nmod_mpolyu_divexact_mpoly_inplace(rBbar, cont, ctx);
                    if (nmod_mpolyuu_divides(rG, B, rBbar, 2, ctx) &&
                        nmod_mpolyuu_divides(rAbar, A, rG, 2, ctx))
                    {
                        break;
                    }
                }
                else
                {
                    nmod_mpolyu_mul_mpoly(T, Bevals + m + 1, gammaevals + m + 1, ctx);
                    if (nmod_mpolyuu_divides(rG, T, rBbar, 2, ctx))
                    {
                        nmod_mpolyu_mul_mpoly(T, Aevals + m + 1, gammaevals + m + 1, ctx);
                        if (nmod_mpolyuu_divides(rAbar, T, rG, 2, ctx))
                        {
                            nmod_mpolyu_swap(G, rG, ctx);
                            nmod_mpolyu_swap(Abar, rAbar, ctx);
                            nmod_mpolyu_swap(Bbar, rBbar, ctx);
                            break;
                        }
                    }
                }
            }

            c = nmod_neg(alphas[m], mod);
            n_poly_mod_shift_left_scalar_addmul((n_poly_struct*)modulus, 1, c, mod);

            if (nmod_poly_degree(modulus) - 1 > gammadegs[m] + Adegs[m] &&
                nmod_poly_degree(modulus) - 1 > gammadegs[m] + Bdegs[m])
            {
                goto choose_alphas;
            }
        }
    }

    success = 1;

cleanup:

    n_polyun_clear(HG);
    n_polyun_clear(HAbar);
    n_polyun_clear(HBbar);
    n_polyun_clear(MG);
    n_polyun_clear(MAbar);
    n_polyun_clear(MBbar);
    n_polyun_clear(ZG);
    n_polyun_clear(ZAbar);
    n_polyun_clear(ZBbar);
    n_bpoly_clear(Aev);
    n_bpoly_clear(Bev);
    n_bpoly_clear(Gev);
    n_bpoly_clear(Abarev);
    n_bpoly_clear(Bbarev);
    n_poly_clear(alphapow);
    nmod_mpoly_clear(cont, ctx);
    nmod_mpolyu_clear(T, ctx);
    nmod_mpolyu_clear(G, ctx);
    nmod_mpolyu_clear(Abar, ctx);
    nmod_mpolyu_clear(Bbar, ctx);
    nmod_mpolyun_clear(Tn, ctx);
    nmod_mpolyun_clear(Gn, ctx);
    nmod_mpolyun_clear(Abarn, ctx);
    nmod_mpolyun_clear(Bbarn, ctx);
    n_polyun_clear(Aeh);
    n_polyun_clear(Beh);
    n_poly_clear(gammaeh);
    nmod_poly_clear(modulus);
    n_poly_stack_clear(St->poly_stack);
    n_bpoly_stack_clear(St->bpoly_stack);

    flint_free(betas);
    flint_free(alphas);
    flint_randclear(state);

    for (i = 0; i < nvars; i++)
    {
        nmod_mpolyu_clear(Aevals + i, ctx);
        nmod_mpolyu_clear(Bevals + i, ctx);
        nmod_mpoly_clear(gammaevals + i, ctx);
    }
    flint_free(Aevals);
    flint_free(Bevals);
    flint_free(gammaevals);

    flint_free(gammadegs);
    flint_free(Adegs);
    flint_free(Bdegs);
/*
timeit_stop(timer);
flint_printf("(%wd ms) mpolyuu_gcd returning %d\n", timer->wall, success);

flint_printf("rG: ");
nmod_mpolyuu_print_pretty(rG, NULL, 2, ctx);
flint_printf("\n");
flint_printf("rAbar: ");
nmod_mpolyuu_print_pretty(rAbar, NULL, 2, ctx);
flint_printf("\n");
flint_printf("rBbar: ");
nmod_mpolyuu_print_pretty(rBbar, NULL, 2, ctx);
flint_printf("\n");
*/
    return success;

gcd_is_trivial:

    nmod_mpolyu_one(rG, ctx);
    nmod_mpolyu_set(rAbar, A, ctx);
    nmod_mpolyu_set(rBbar, B, ctx);

    success = 1;
    
    goto cleanup;
}

