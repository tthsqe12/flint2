/*
    Copyright (C) 2019 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mod.h"
#include "nmod_mpoly.h"
#include "fmpz_mpoly.h"
#include "fmpz_mod_mpoly.h"
#include "n_poly.h"
#include "nmod_mpoly_factor.h"


void n_bpoly_mod_interp_reduce_2sm_poly(
    n_poly_t Ap,
    n_poly_t Am,
    const n_bpoly_t A,
    n_poly_t alphapow,
    nmod_t mod);

void n_bpoly_mod_interp_lift_2sm_poly(
    slong * deg1,
    n_bpoly_t T,
    const n_poly_t A,
    const n_poly_t B,    
    mp_limb_t alpha,
    nmod_t mod);

int n_bpoly_mod_interp_crt_2sm_poly(
    slong * deg1,
    n_bpoly_t F,
    n_bpoly_t T,
    n_poly_t A,
    n_poly_t B,
    const n_poly_t modulus,
    n_poly_t alphapow,
    nmod_t mod);

typedef struct {
    n_poly_stack_t poly_stack;
    n_bpoly_stack_t bpoly_stack;
} n_poly_bpoly_stack_struct;

typedef n_poly_bpoly_stack_struct n_poly_bpoly_stack_t[1];


ulong n_bpoly_bidegree(const n_bpoly_t A)
{
    ulong x, y;

    FLINT_ASSERT(A->length > 0);

    x = A->length - 1;
    y = A->coeffs[x].length - 1;

    return (x << (FLINT_BITS/2)) + y;
}


void n_bpoly_scalar_mul_nmod(
    n_bpoly_t A,
    mp_limb_t c,
    nmod_t ctx)
{
    slong i;

    FLINT_ASSERT(c != 0);

    if (c == 1)
        return;

    for (i = 0; i < A->length; i++)
        _n_poly_mod_scalar_mul_nmod_inplace(A->coeffs + i, c, ctx);
}


void n_bpoly_mod_content_last(n_poly_t g, const n_bpoly_t A, nmod_t ctx)
{
    slong i;

    n_poly_zero(g);
    for (i = 0; i < A->length; i++)
    {
        n_poly_mod_gcd(g, g, A->coeffs + i, ctx);
        if (n_poly_is_one(g))
            break;
    }
}

void n_bpoly_mod_divexact_last(n_bpoly_t A, const n_poly_t b, nmod_t ctx)
{
    slong i;
    n_poly_struct * t;

    if (b->length == 1)
    {
        if (b->coeffs[0] != 1)
            n_bpoly_scalar_mul_nmod(A, nmod_inv(b->coeffs[0], ctx), ctx);
        return;
    }

    n_bpoly_fit_length(A, A->length + 1);

    t = A->coeffs + A->length;

    for (i = 0; i < A->length; i++)
    {
        if (A->coeffs[i].length < 1)
            continue;

        n_poly_mod_div(t, A->coeffs + i, b, ctx);
        n_poly_swap(A->coeffs + i, t);
    }
}


void n_bpoly_mod_mul_last(n_bpoly_t A, const n_poly_t b, nmod_t ctx)
{
    slong i;
    n_poly_struct * t;

    if (n_poly_is_one(b))
        return;

    n_bpoly_fit_length(A, A->length + 1);

    t = A->coeffs + A->length;

    for (i = 0; i < A->length; i++)
    {
        if (A->coeffs[i].length < 1)
            continue;

        n_poly_mod_mul(t, A->coeffs + i, b, ctx);
        n_poly_swap(A->coeffs + i, t);
    }
}



int n_bpoly_mod_gcd_brown_smprime(
    n_bpoly_t G,
    n_bpoly_t Abar,
    n_bpoly_t Bbar,
    n_bpoly_t A,
    n_bpoly_t B,
    nmod_t ctx,
    n_poly_bpoly_stack_t Sp)
{
    int success;
    slong bound;
    mp_limb_t alpha, temp, gammaevalp, gammaevalm;
    n_poly_struct * Aevalp, * Bevalp, * Gevalp, * Abarevalp, * Bbarevalp;
    n_poly_struct * Aevalm, * Bevalm, * Gevalm, * Abarevalm, * Bbarevalm;
    n_bpoly_struct * T;
    slong deggamma, ldegG, ldegAbar, ldegBbar, ldegA, ldegB;
    n_poly_struct * cA, * cB, * cG, * cAbar, * cBbar, * gamma;
    n_poly_struct * modulus, * alphapow, * r;
    int gstab, astab, bstab, use_stab;
#if WANT_ASSERT
    n_poly_t leadA, leadB;
    const slong Sp_size_poly = n_poly_stack_size(Sp->poly_stack);
    const slong Sp_size_bpoly = n_bpoly_stack_size(Sp->bpoly_stack);
#endif

    FLINT_ASSERT(A->length > 0);
    FLINT_ASSERT(B->length > 0);

#if WANT_ASSERT
    n_poly_init(leadA);
    n_poly_init(leadB);
    n_poly_set(leadA, A->coeffs + A->length - 1);
    n_poly_set(leadB, B->coeffs + B->length - 1);
#endif

    n_poly_stack_fit_request(Sp->poly_stack, 19);
    cA          = n_poly_stack_take_top(Sp->poly_stack);
    cB          = n_poly_stack_take_top(Sp->poly_stack);
    cG          = n_poly_stack_take_top(Sp->poly_stack);
    cAbar       = n_poly_stack_take_top(Sp->poly_stack);
    cBbar       = n_poly_stack_take_top(Sp->poly_stack);
    gamma       = n_poly_stack_take_top(Sp->poly_stack);
    Aevalp      = n_poly_stack_take_top(Sp->poly_stack);
    Bevalp      = n_poly_stack_take_top(Sp->poly_stack);
    Gevalp      = n_poly_stack_take_top(Sp->poly_stack);
    Abarevalp   = n_poly_stack_take_top(Sp->poly_stack);
    Bbarevalp   = n_poly_stack_take_top(Sp->poly_stack);
    Aevalm      = n_poly_stack_take_top(Sp->poly_stack);
    Bevalm      = n_poly_stack_take_top(Sp->poly_stack);
    Gevalm      = n_poly_stack_take_top(Sp->poly_stack);
    Abarevalm   = n_poly_stack_take_top(Sp->poly_stack);
    Bbarevalm   = n_poly_stack_take_top(Sp->poly_stack);
    r           = n_poly_stack_take_top(Sp->poly_stack);
    alphapow    = n_poly_stack_take_top(Sp->poly_stack);
    modulus     = n_poly_stack_take_top(Sp->poly_stack);

    n_bpoly_stack_fit_request(Sp->bpoly_stack, 1);
    T           = n_bpoly_stack_take_top(Sp->bpoly_stack);

    n_bpoly_mod_content_last(cA, A, ctx);
    n_bpoly_mod_content_last(cB, B, ctx);
    n_bpoly_mod_divexact_last(A, cA, ctx);
    n_bpoly_mod_divexact_last(B, cB, ctx);

    n_poly_mod_gcd(cG, cA, cB, ctx);

    n_poly_mod_div(cAbar, cA, cG, ctx);
    n_poly_mod_div(cBbar, cB, cG, ctx);

    n_poly_mod_gcd(gamma, A->coeffs + A->length - 1, B->coeffs + B->length - 1, ctx);

    ldegA = n_bpoly_degree1(A);
    ldegB = n_bpoly_degree1(B);
    deggamma = n_poly_degree(gamma);
    bound = 1 + deggamma + FLINT_MAX(ldegA, ldegB);

    n_poly_fit_length(alphapow, FLINT_MAX(WORD(3), bound + 1));
    n_poly_one(modulus);

    use_stab = 1;
    gstab = bstab = astab = 0;

    if ((ctx.n & UWORD(1)) == UWORD(0))
    {
        success = 0;
        goto cleanup;
    }

    alpha = (ctx.n - UWORD(1))/UWORD(2);

choose_prime: /* primes are v - alpha, v + alpha */

    if (alpha < 2)
    {
        success = 0;
        goto cleanup;
    }

    alpha--;

    FLINT_ASSERT(0 < alpha && alpha <= ctx.n/2);
    FLINT_ASSERT(alphapow->alloc >= 2);
    alphapow->length = 2;
    alphapow->coeffs[0] = 1;
    alphapow->coeffs[1] = alpha;

    /* make sure evaluation point does not kill both lc(A) and lc(B) */
    n_poly_mod_eval2_pow(&gammaevalp, &gammaevalm, gamma, alphapow, ctx);
    if (gammaevalp == 0 || gammaevalm == 0)
        goto choose_prime;

    /* evaluation point should kill neither A nor B */
    n_bpoly_mod_interp_reduce_2sm_poly(Aevalp, Aevalm, A, alphapow, ctx);
    n_bpoly_mod_interp_reduce_2sm_poly(Bevalp, Bevalm, B, alphapow, ctx);
    FLINT_ASSERT(Aevalp->length > 0);
    FLINT_ASSERT(Aevalm->length > 0);
    FLINT_ASSERT(Bevalp->length > 0);
    FLINT_ASSERT(Bevalm->length > 0);

    if (use_stab && gstab)
    {
        slong Gdeg;
        n_bpoly_mod_interp_reduce_2sm_poly(Gevalp, Gevalm, G, alphapow, ctx);
        Gdeg = n_bpoly_degree0(G);
        success = 1;
        success = success && n_poly_degree(Gevalp) == Gdeg;
        success = success && n_poly_degree(Gevalm) == Gdeg;
        success = success && Gevalp->coeffs[Gdeg] == gammaevalp;
        success = success && Gevalm->coeffs[Gdeg] == gammaevalm;
        n_poly_mod_divrem(Abarevalp, r, Aevalp, Gevalp, ctx);
        success = success && (r->length == 0);
        n_poly_mod_divrem(Abarevalm, r, Aevalm, Gevalm, ctx);
        success = success && (r->length == 0);
        n_poly_mod_divrem(Bbarevalp, r, Bevalp, Gevalp, ctx);
        success = success && (r->length == 0);
        n_poly_mod_divrem(Bbarevalm, r, Bevalm, Gevalm, ctx);
        success = success && (r->length == 0);

        if (!success)
        {
            use_stab = 0;
            n_poly_one(modulus);
            alpha = (ctx.n - UWORD(1))/UWORD(2);
            goto choose_prime;
        }

        _n_poly_mod_scalar_mul_nmod_inplace(Abarevalp, gammaevalp, ctx);
        _n_poly_mod_scalar_mul_nmod_inplace(Abarevalm, gammaevalm, ctx);
        _n_poly_mod_scalar_mul_nmod_inplace(Bbarevalp, gammaevalp, ctx);
        _n_poly_mod_scalar_mul_nmod_inplace(Bbarevalm, gammaevalm, ctx);
    }
    else
    {
        n_poly_mod_gcd(Gevalp, Aevalp, Bevalp, ctx);
        n_poly_mod_div(Abarevalp, Aevalp, Gevalp, ctx);
        n_poly_mod_div(Bbarevalp, Bevalp, Gevalp, ctx);
        n_poly_mod_gcd(Gevalm, Aevalm, Bevalm, ctx);
        n_poly_mod_div(Abarevalm, Aevalm, Gevalm, ctx);
        n_poly_mod_div(Bbarevalm, Bevalm, Gevalm, ctx);
        gstab = astab = bstab = 0;
    }

    FLINT_ASSERT(Gevalp->length > 0);
    FLINT_ASSERT(Abarevalp->length > 0);
    FLINT_ASSERT(Bbarevalp->length > 0);
    FLINT_ASSERT(Gevalm->length > 0);
    FLINT_ASSERT(Abarevalm->length > 0);
    FLINT_ASSERT(Bbarevalm->length > 0);

    if (n_poly_degree(Gevalp) == 0 || n_poly_degree(Gevalm) == 0)
    {
        n_bpoly_one(G);
        n_bpoly_swap(Abar, A);
        n_bpoly_swap(Bbar, B);
        goto successful_put_content;    
    }

    if (n_poly_degree(Gevalp) != n_poly_degree(Gevalm))
    {
        goto choose_prime;
    }

    /* the Geval have matching degrees */
    if (n_poly_degree(modulus) > 0)
    {
        FLINT_ASSERT(G->length > 0);
        if (n_poly_degree(Gevalp) > n_bpoly_degree0(G))
        {
            goto choose_prime;
        }
        else if (n_poly_degree(Gevalp) < n_bpoly_degree0(G))
        {
            n_poly_one(modulus);
        }
    }

    /* update interpolants */
    _n_poly_mod_scalar_mul_nmod_inplace(Gevalp, gammaevalp, ctx);
    _n_poly_mod_scalar_mul_nmod_inplace(Gevalm, gammaevalm, ctx);
    if (n_poly_degree(modulus) > 0)
    {
        temp = n_poly_mod_evaluate_nmod(modulus, alpha, ctx);
        FLINT_ASSERT(temp == n_poly_mod_evaluate_nmod(modulus, ctx.n - alpha, ctx));
        temp = nmod_mul(temp, alpha, ctx);
        temp = nmod_add(temp, temp, ctx);
        temp = nmod_inv(temp, ctx);
        _n_poly_mod_scalar_mul_nmod_inplace(modulus, temp, ctx);
        if (!gstab)
        {
            gstab = !n_bpoly_mod_interp_crt_2sm_poly(&ldegG, G, T,
                                       Gevalp, Gevalm, modulus, alphapow, ctx);
        }
        n_bpoly_mod_interp_crt_2sm_poly(&ldegAbar, Abar, T,
                                 Abarevalp, Abarevalm, modulus, alphapow, ctx);
        n_bpoly_mod_interp_crt_2sm_poly(&ldegBbar, Bbar, T,
                                 Bbarevalp, Bbarevalm, modulus, alphapow, ctx);
    }
    else
    {
        n_bpoly_mod_interp_lift_2sm_poly(&ldegG, G, Gevalp, Gevalm, alpha, ctx);
        n_bpoly_mod_interp_lift_2sm_poly(&ldegAbar, Abar,
                                             Abarevalp, Abarevalm, alpha, ctx);
        n_bpoly_mod_interp_lift_2sm_poly(&ldegBbar, Bbar,
                                             Bbarevalp, Bbarevalm, alpha, ctx);
        gstab = astab = bstab = 0;
    }

    temp = ctx.n - nmod_mul(alpha, alpha, ctx);
    n_poly_mod_shift_left_scalar_addmul(modulus, 2, temp, ctx);

    if (n_poly_degree(modulus) < bound)
        goto choose_prime;

    FLINT_ASSERT(ldegG >= 0);
    FLINT_ASSERT(ldegAbar >= 0);
    FLINT_ASSERT(ldegBbar >= 0);

    if (   deggamma + ldegA == ldegG + ldegAbar
        && deggamma + ldegB == ldegG + ldegBbar )
    {
        goto successful;
    }

    n_poly_one(modulus);
    goto choose_prime;

successful:

    n_bpoly_mod_content_last(modulus, G, ctx);
    n_bpoly_mod_divexact_last(G, modulus, ctx);
    n_bpoly_mod_divexact_last(Abar, G->coeffs + G->length - 1, ctx);
    n_bpoly_mod_divexact_last(Bbar, G->coeffs + G->length - 1, ctx);

successful_put_content:

    n_bpoly_mod_mul_last(G, cG, ctx);
    n_bpoly_mod_mul_last(Abar, cAbar, ctx);
    n_bpoly_mod_mul_last(Bbar, cBbar, ctx);

    success = 1;

cleanup:

#if WANT_ASSERT
    if (success)
    {
        n_poly_struct * Glead = G->coeffs + G->length - 1;
        FLINT_ASSERT(1 == Glead->coeffs[Glead->length - 1]);
        n_poly_mod_mul(modulus, G->coeffs + G->length - 1, Abar->coeffs + Abar->length - 1, ctx);
        FLINT_ASSERT(n_poly_equal(modulus, leadA));
        n_poly_mod_mul(modulus, G->coeffs + G->length - 1, Bbar->coeffs + Bbar->length - 1, ctx);
        FLINT_ASSERT(n_poly_equal(modulus, leadB));
    }
    n_poly_clear(leadA);
    n_poly_clear(leadB);
#endif

    n_poly_stack_give_back(Sp->poly_stack, 19);
    n_bpoly_stack_give_back(Sp->bpoly_stack, 1);
    FLINT_ASSERT(Sp_size_poly == n_poly_stack_size(Sp->poly_stack));
    FLINT_ASSERT(Sp_size_bpoly == n_bpoly_stack_size(Sp->bpoly_stack));

    return success;
}



#define LOW_HALF_MASK ((-UWORD(1)) >> (FLINT_BITS - FLINT_BITS/2))

/*
    If I was formed from evaluations at
        alpha^alphashift, alpha^(alphashift + 1), ...
    construct the corresponding mpoly if possible with coeffs in (-p/2, p/2]
    The substitution degrees and degree bounds in Ictx are used.
*/
int nmod_mpoly_bma_get_fmpz_mpoly(
    fmpz_mpoly_t A,
    const fmpz_mpoly_ctx_t ctx,
    ulong alphashift,
    nmod_berlekamp_massey_t I,
    const mpoly_bma_interpolate_ctx_t Ictx,
    const nmodf_ctx_t fpctx)
{
    slong i, j, t, N;
    int success;
    ulong new_exp, this_exp;
    slong * shifts, * offsets;
    mp_limb_t * values, * roots;
    fmpz * Acoeff;
    ulong * Aexp;
    slong Alen;
    mp_limb_t T, S, V, V0, V1, V2, p0, p1, r;
    TMP_INIT;

    t = nmod_poly_degree(I->V1);
    FLINT_ASSERT(I->points->length >= t);

    if (t <= 0)
    {
        return 0;
    }

    TMP_START;

    /* use the rt member of I as temp space for roots - slightly dirty */
    nmod_poly_fit_length(I->rt, t);
    I->rt->length = t;
    roots = I->rt->coeffs;
    values = I->points->coeffs;

    success = nmod_poly_find_distinct_nonzero_roots(roots, I->V1);
    if (!success)
    {
        goto cleanup;
    }

    FLINT_ASSERT(ctx->minfo->ord == ORD_LEX);
    FLINT_ASSERT(A->bits <= FLINT_BITS);
    N = mpoly_words_per_exp_sp(A->bits, ctx->minfo);
    fmpz_mpoly_fit_length(A, t, ctx);
    Acoeff = A->coeffs;
    Aexp = A->exps;
    A->length = 0;

    shifts = (slong *) TMP_ALLOC(ctx->minfo->nvars*sizeof(slong));
    offsets = (slong *) TMP_ALLOC(ctx->minfo->nvars*sizeof(slong));
    for (j = 0; j < ctx->minfo->nvars; j++)
    {
        mpoly_gen_offset_shift_sp(offsets + j, shifts + j, j, A->bits, ctx->minfo);
    }

    Alen = 0;
    for (i = 0; i < t; i++)
    {
        /* coeffs[i] is (coeffs(P).values)/P(roots[i]) =: V/S
            where P(x) = master(x)/(x-roots[i])     */
        V0 = V1 = V2 = T = S = 0;
        r = roots[i];
        for (j = t; j > 0; j--)
        {
            T = nmod_add(nmod_mul(r, T, fpctx->mod), I->V1->coeffs[j], fpctx->mod);
            S = nmod_add(nmod_mul(r, S, fpctx->mod), T, fpctx->mod);
            umul_ppmm(p1, p0, values[j - 1], T);
            add_sssaaaaaa(V2, V1, V0, V2, V1, V0, WORD(0), p1, p0);
        }
        /* roots[i] should be a root of master */
        FLINT_ASSERT(nmod_add(nmod_mul(r, T, fpctx->mod), I->V1->coeffs[0], fpctx->mod) == 0);
        NMOD_RED3(V, V2, V1, V0, fpctx->mod);
        S = nmod_mul(S, nmod_pow_ui(r, alphashift, fpctx->mod), fpctx->mod);
        V0 = nmod_mul(V, nmod_inv(S, fpctx->mod), fpctx->mod);
        if (V0 == 0)
        {
            /* hmmm */
            continue;
        }

        fmpz_set_ui(Acoeff + Alen, V0);
        if (fpctx->mod.n - V0 < V0)
        {
            fmpz_sub_ui(Acoeff + Alen, Acoeff + Alen, fpctx->mod.n);
        }

        mpoly_monomial_zero(Aexp + N*Alen, N);
        new_exp = nmod_discrete_log_pohlig_hellman_run(Ictx->dlogenv_sp, roots[i]);
        for (j = ctx->minfo->nvars - 1; j >= 0; j--)
        {
            this_exp = new_exp % Ictx->subdegs[j];
            new_exp = new_exp / Ictx->subdegs[j];
            if (this_exp >= Ictx->degbounds[j])
            {
                success = 0;
                goto cleanup;
            }
            (Aexp + N*Alen)[offsets[j]] |= this_exp << shifts[j];
        }
        if (new_exp != 0)
        {
            success = 0;
            goto cleanup;
        }
        Alen++;
    }
    A->length = Alen;

    fmpz_mpoly_sort_terms(A, ctx);

    success = 1;

cleanup:

    TMP_END;
    return success;
}

/*
    nmod_mpoly "skeletons" - just the coefficients
*/

void nmod_mpolyc_init(nmod_mpolyc_t A)
{
    A->coeffs = NULL;
    A->alloc = 0;
    A->length = 0;
}

void nmod_mpolyc_clear(nmod_mpolyc_t A)
{
    if (A->coeffs)
        flint_free(A->coeffs);
}

void nmod_mpolyc_fit_length(nmod_mpolyc_t A, slong length)
{
    slong old_alloc = A->alloc;
    slong new_alloc = FLINT_MAX(length, A->alloc + A->alloc/2);

    if (length > old_alloc)
    {
        if (old_alloc == 0)
        {
            A->coeffs = (mp_limb_t *) flint_malloc(new_alloc*sizeof(mp_limb_t));
        }
        else
        {
            A->coeffs = (mp_limb_t *) flint_realloc(A->coeffs, new_alloc*sizeof(mp_limb_t));
        }

        A->alloc = new_alloc;
    }
}

void nmod_mpolycu_init(nmod_mpolycu_t A)
{
    A->coeffs = NULL;
    A->alloc = 0;
    A->length = 0;
}

void nmod_mpolycu_clear(nmod_mpolycu_t A)
{
    slong i;
    for (i = 0; i < A->alloc; i++)
    {
        nmod_mpolyc_clear(A->coeffs + i);
    }
    if (A->coeffs)
        flint_free(A->coeffs);
}

void nmod_mpolycu_fit_length(nmod_mpolycu_t A, slong length)
{
    slong i;
    slong old_alloc = A->alloc;
    slong new_alloc = FLINT_MAX(length, A->alloc + A->alloc/2);

    if (length > old_alloc)
    {
        if (old_alloc == 0)
        {
            A->coeffs = (nmod_mpolyc_struct *) flint_malloc(
                                         new_alloc*sizeof(nmod_mpolyc_struct));
        }
        else
        {
            A->coeffs = (nmod_mpolyc_struct *) flint_realloc(A->coeffs,
                                         new_alloc*sizeof(nmod_mpolyc_struct));
        }

        for (i = old_alloc; i < new_alloc; i++)
        {
            nmod_mpolyc_init(A->coeffs + i);
        }
        A->alloc = new_alloc;
    }
}

/*
    fmpz_mod_mpoly "skeletons" - just the coefficients
*/

void fmpz_mpolyc_init(fmpz_mpolyc_t A)
{
    A->coeffs = NULL;
    A->alloc = 0;
    A->length = 0;
}

void fmpz_mpolyc_clear(fmpz_mpolyc_t A)
{
    slong i;
    for (i = 0; i < A->alloc; i++)
    {
        fmpz_clear(A->coeffs + i);
    }
    if (A->coeffs)
        flint_free(A->coeffs);
}

void fmpz_mpolyc_fit_length(fmpz_mpolyc_t A, slong length)
{
    slong i;
    slong old_alloc = A->alloc;
    slong new_alloc = FLINT_MAX(length, A->alloc + A->alloc/2);

    if (length > old_alloc)
    {
        if (old_alloc == 0)
        {
            A->coeffs = (fmpz *) flint_malloc(new_alloc*sizeof(fmpz));
        }
        else
        {
            A->coeffs = (fmpz *) flint_realloc(A->coeffs, new_alloc*sizeof(fmpz));
        }

        for (i = old_alloc; i < new_alloc; i++)
        {
            fmpz_init(A->coeffs + i);
        }
        A->alloc = new_alloc;
    }
}

void fmpz_mpolycu_init(fmpz_mpolycu_t A)
{
    A->coeffs = NULL;
    A->alloc = 0;
    A->length = 0;
}

void fmpz_mpolycu_clear(fmpz_mpolycu_t A)
{
    slong i;
    for (i = 0; i < A->alloc; i++)
    {
        fmpz_mpolyc_clear(A->coeffs + i);
    }
    if (A->coeffs)
        flint_free(A->coeffs);
}

void fmpz_mpolycu_fit_length(fmpz_mpolycu_t A, slong length)
{
    slong i;
    slong old_alloc = A->alloc;
    slong new_alloc = FLINT_MAX(length, A->alloc + A->alloc/2);

    if (length > old_alloc)
    {
        if (old_alloc == 0)
        {
            A->coeffs = (fmpz_mpolyc_struct *) flint_malloc(
                                         new_alloc*sizeof(fmpz_mpolyc_struct));
        }
        else
        {
            A->coeffs = (fmpz_mpolyc_struct *) flint_realloc(A->coeffs,
                                         new_alloc*sizeof(fmpz_mpolyc_struct));
        }

        for (i = old_alloc; i < new_alloc; i++)
        {
            fmpz_mpolyc_init(A->coeffs + i);
        }
        A->alloc = new_alloc;
    }
}


/*
    set out to the evaluation of variables after ksub at alpha^w

    out[0]   = alpha ^ (w * subdegs[n-1] * subdegs[n-2] * ... * * subdegs[1])
      ...
    out[n-3] = alpha ^ (w * subdegs[n-1] * subdegs[n-2])
    out[n-2] = alpha ^ (w * subdegs[n-1])
    out[n-1] = alpha ^ (w)

    secret: subdegs[0] is not used
*/

void nmod_mpoly_bma_interpolate_alpha_powers(
    mp_limb_t * out,
    ulong w,
    const mpoly_bma_interpolate_ctx_t Ictx,
    const fmpz_mpoly_ctx_t ctx,
    const nmodf_ctx_t fpctx)
{
    slong j = ctx->minfo->nvars - 1;
    out[j] = nmod_pow_ui(Ictx->dlogenv_sp->alpha, w, fpctx->mod);
    for (; j > 0; j--)
    {
        out[j - 1] = nmod_pow_ui(out[j], Ictx->subdegs[j], fpctx->mod);
    }
}

void fmpz_mod_mpoly_bma_interpolate_alpha_powers(
    fmpz * out,
    const fmpz_t w,
    const mpoly_bma_interpolate_ctx_t Ictx,
    const fmpz_mpoly_ctx_t ctx,
    const fmpz_mod_ctx_t fpctx)
{
    slong j = ctx->minfo->nvars - 1;
    fmpz_mod_pow_fmpz(out + j, Ictx->dlogenv->alpha, w, fpctx);
    for (; j > 0; j--)
    {
        fmpz_mod_pow_ui(out + j - 1, out + j, Ictx->subdegs[j], fpctx);
    }
}


/*
    Set the coefficients of E to the evaluation of cooresponding monomials of A
    evaluation at x_i = alpha[i]
*/
void fmpz_mod_mpoly_set_skel(
    fmpz_mpolyc_t M,
    const fmpz_mod_mpoly_ctx_t ctx_mp,
    const fmpz_mpoly_t A,
    const fmpz * alpha,
    const fmpz_mpoly_ctx_t ctx)
{
    slong i, j;
    slong offset, shift;
    slong N = mpoly_words_per_exp_sp(A->bits, ctx->minfo);
    slong nvars = ctx->minfo->nvars;
    ulong * Aexp;
    fmpz * Mcoeff;
    slong * LUToffset;
    ulong * LUTmask;
    fmpz * LUTvalue;
    slong LUTlen;
    fmpz_t xpoweval;
    ulong * inputexpmask;
    TMP_INIT;

    TMP_START;

    fmpz_init(xpoweval);

    LUToffset = (slong *) TMP_ALLOC(N*FLINT_BITS*sizeof(slong));
    LUTmask   = (ulong *) TMP_ALLOC(N*FLINT_BITS*sizeof(ulong));
    LUTvalue  = (fmpz *) TMP_ALLOC(N*FLINT_BITS*sizeof(fmpz));
    for (i = 0; i < N*FLINT_BITS; i++)
    {
        fmpz_init(LUTvalue + i);
    }

    fmpz_mpolyc_fit_length(M, A->length);
    M->length = A->length;

    Mcoeff = M->coeffs;
    Aexp = A->exps;

    inputexpmask = (ulong *) TMP_ALLOC(N*sizeof(ulong));
    mpoly_monomial_zero(inputexpmask, N);
    for (i = 0; i < A->length; i++)
    {
        for (j = 0; j < N; j++)
        {
            inputexpmask[j] |= (Aexp + N*i)[j];
        }
    }

    LUTlen = 0;
    for (j = nvars - 1; j >= 0; j--)
    {
        mpoly_gen_offset_shift_sp(&offset, &shift, j, A->bits, ctx->minfo);

        fmpz_set(xpoweval, alpha + j); /* xpoweval = alpha[j]^(2^i) */
        for (i = 0; i < A->bits; i++)
        {
            LUToffset[LUTlen] = offset;
            LUTmask[LUTlen] = (UWORD(1) << (shift + i));
            fmpz_set(LUTvalue + LUTlen, xpoweval);
            if ((inputexpmask[offset] & LUTmask[LUTlen]) != 0)
            {
                LUTlen++;
            }
            fmpz_mod_mul(xpoweval, xpoweval, xpoweval, ctx_mp->ffinfo);
        }
    }
    FLINT_ASSERT(LUTlen < N*FLINT_BITS);

    for (i = 0; i < A->length; i++)
    {
        fmpz_one(xpoweval);
        for (j = 0; j < LUTlen; j++)
        {
            if (((Aexp + N*i)[LUToffset[j]] & LUTmask[j]) != 0)
            {
                fmpz_mod_mul(xpoweval, xpoweval, LUTvalue + j, ctx_mp->ffinfo);
            }
        }
        fmpz_set(Mcoeff + i, xpoweval);
    }

    fmpz_clear(xpoweval);
    for (i = 0; i < N*FLINT_BITS; i++)
    {
        fmpz_clear(LUTvalue + i);
    }
    TMP_END;
}

void fmpz_mod_mpolyu_set_skel(
    fmpz_mpolycu_t M,
    const fmpz_mod_mpoly_ctx_t ctx_mp,
    const fmpz_mpolyu_t A,
    const fmpz * alpha,
    const fmpz_mpoly_ctx_t ctx)
{
    slong i;
    fmpz_mpolycu_fit_length(M, A->length);
    M->length = A->length;
    for (i = 0; i < A->length; i++)
    {
        fmpz_mod_mpoly_set_skel(M->coeffs + i, ctx_mp, A->coeffs + i, alpha, ctx);
    }
}


/* M = S */
void fmpz_mod_mpoly_copy_skel(fmpz_mpolyc_t M, const fmpz_mpolyc_t S)
{
    fmpz_mpolyc_fit_length(M, S->length);
    M->length = S->length;
    _fmpz_vec_set(M->coeffs, S->coeffs, S->length);
}

void fmpz_mod_mpolyu_copy_skel(fmpz_mpolycu_t M, const fmpz_mpolycu_t S)
{
    slong i;
    fmpz_mpolycu_fit_length(M, S->length);
    M->length = S->length;
    for (i = 0; i < S->length; i++)
    {
        fmpz_mod_mpoly_copy_skel(M->coeffs + i, S->coeffs + i);
    }
}


/* Ared = A mod p */
void fmpz_mod_mpoly_red_skel(
    fmpz_mpolyc_t Ared,
    const fmpz_mpoly_t A,
    const fmpz_mod_ctx_t fpctx)
{
    fmpz_mpolyc_fit_length(Ared, A->length);
    Ared->length = A->length;
    _fmpz_vec_scalar_mod_fmpz(Ared->coeffs, A->coeffs, A->length,
                                                  fmpz_mod_ctx_modulus(fpctx));
}

void fmpz_mod_mpolyu_red_skel(
    fmpz_mpolycu_t Ared,
    const fmpz_mpolyu_t A,
    const fmpz_mod_ctx_t fpctx)
{
    slong i;
    fmpz_mpolycu_fit_length(Ared, A->length);
    Ared->length = A->length;
    for (i = 0; i < A->length; i++)
    {
        fmpz_mod_mpoly_red_skel(Ared->coeffs + i, A->coeffs + i, fpctx);
    }
}


/*
    return Ared.Avar and multiply Avar by Ainc
    the coefficients of Ared must be reduced mod fpctx->p
*/
void fmpz_mod_mpoly_use_skel_mul(
    fmpz_t eval,
    fmpz_mpolyc_t Ared,
    fmpz_mpolyc_t Avar,
    const fmpz_mpolyc_t Ainc,
    const fmpz_mod_ctx_t fpctx)
{
    slong i;
    fmpz_zero(eval);
    FLINT_ASSERT(Avar->length == Ared->length);
    for (i = 0; i < Ared->length; i++)
    {
        fmpz_addmul(eval, Ared->coeffs + i, Avar->coeffs + i);
        fmpz_mod_mul(Avar->coeffs + i, Avar->coeffs + i, Ainc->coeffs + i, fpctx);
    }
    fmpz_mod(eval, eval, fmpz_mod_ctx_modulus(fpctx));
}


/*
    A is in ZZ[X,Y][x_0, ..., x_(n-1)]
    E is in Fp[X][Y]

    Set E to the evaluation of A. The evaluations of A's monomials are in M.
    M is then multiplied by S.
*/
void fmpz_mod_mpolyuu_use_skel_mul(
    fmpz_mod_mpolyn_t E,
    const fmpz_mpolyu_t A,
    const fmpz_mpolycu_t Ared,
    fmpz_mpolycu_t M,
    const fmpz_mpolycu_t S,
    const fmpz_mod_mpoly_ctx_t ctx_mp)
{
    slong xexp, yexp;
    slong i;
    fmpz_t eval;

    FLINT_ASSERT(E->bits == FLINT_BITS/2);
    FLINT_ASSERT(1 == mpoly_words_per_exp_sp(E->bits, ctx_mp->minfo));

    FLINT_ASSERT(A->length == Ared->length);
    FLINT_ASSERT(A->length == M->length);
    FLINT_ASSERT(A->length == S->length);

    fmpz_init(eval);

    E->length = 0;
    for (i = 0; i < A->length; i++)
    {
        fmpz_mod_mpoly_use_skel_mul(eval, Ared->coeffs + i, M->coeffs + i,
                                                S->coeffs + i, ctx_mp->ffinfo);
        if (fmpz_is_zero(eval))
        {
            continue;
        }

        xexp = A->exps[i] >> (FLINT_BITS/2);
        yexp = A->exps[i] & LOW_HALF_MASK;

        if (E->length > 0 && (E->exps[E->length - 1] >> (FLINT_BITS/2)) == xexp)
        {
            fmpz_mod_poly_set_coeff_fmpz(E->coeffs + E->length - 1, yexp, eval);
        }
        else
        {
            fmpz_mod_mpolyn_fit_length(E, E->length + 1, ctx_mp);
            fmpz_mod_poly_zero(E->coeffs + E->length);
            fmpz_mod_poly_set_coeff_fmpz(E->coeffs + E->length, yexp, eval);
            E->exps[E->length] = xexp << (FLINT_BITS/2);
            E->length++;
        }
    }

    fmpz_clear(eval);
}



int fmpz_mod_bma_get_fmpz_mpoly(
    fmpz_mpoly_t A,
    const fmpz_mpoly_ctx_t ctx,
    const fmpz_t alphashift,
    fmpz_mod_berlekamp_massey_t I,
    const mpoly_bma_interpolate_ctx_t Ictx,
    const fmpz_mod_ctx_t fpctx)
{
    slong i, j, t, N;
    int success;
    ulong this_exp;
    fmpz_t new_exp;
    slong * shifts, * offsets;
    fmpz * values, * roots;
    fmpz * Acoeff;
    ulong * Aexp;
    slong Alen;
    fmpz_t T, S, V, temp, halfp;
    TMP_INIT;

    TMP_START;

    fmpz_init(halfp);
    fmpz_init(T);
    fmpz_init(S);
    fmpz_init(V);
    fmpz_init(temp);
    fmpz_init(new_exp);

    fmpz_tdiv_q_2exp(halfp, fmpz_mod_ctx_modulus(fpctx), 1);

    fmpz_mod_berlekamp_massey_reduce(I);
    t = fmpz_mod_poly_degree(I->V1);
    FLINT_ASSERT(t > 0);
    FLINT_ASSERT(I->points->length >= t);

    fmpz_mod_poly_fit_length(I->rt, t);
    I->rt->length = t;
    success = fmpz_mod_poly_find_distinct_nonzero_roots(I->rt->coeffs, I->V1);
    if (!success)
    {
        goto cleanup;
    }

    roots = I->rt->coeffs;
    values = I->points->coeffs;

    FLINT_ASSERT(ctx->minfo->ord == ORD_LEX);
    FLINT_ASSERT(A->bits <= FLINT_BITS);
    N = mpoly_words_per_exp_sp(A->bits, ctx->minfo);
    fmpz_mpoly_fit_length(A, t, ctx);
    Acoeff = A->coeffs;
    Aexp = A->exps;
    A->length = 0;

    shifts = (slong *) TMP_ALLOC(ctx->minfo->nvars*sizeof(slong));
    offsets = (slong *) TMP_ALLOC(ctx->minfo->nvars*sizeof(slong));
    for (j = 0; j < ctx->minfo->nvars; j++)
    {
        mpoly_gen_offset_shift_sp(offsets + j, shifts + j, j, A->bits, ctx->minfo);
    }

    Alen = 0;
    for (i = 0; i < t; i++)
    {
        /* coeffs[i] is (coeffs(P).values)/P(roots[i]) =: V/S
            where P(x) = master(x)/(x-roots[i])     */
        fmpz_zero(V);
        fmpz_zero(T);
        fmpz_zero(S);
        for (j = t; j > 0; j--)
        {
            fmpz_mod_mul(temp, roots + i, T, fpctx);
            fmpz_mod_add(T, temp, I->V1->coeffs + j, fpctx);
            fmpz_mod_mul(temp, roots + i, S, fpctx);
            fmpz_mod_add(S, temp, T, fpctx);
            fmpz_mod_mul(temp, values + j - 1, T, fpctx);
            fmpz_mod_add(V, V, temp, fpctx);
        }
        /* roots[i] should be a root of master */
#if WANT_ASSERT
        fmpz_mod_mul(temp, roots + i, T, fpctx);
        fmpz_mod_add(temp, temp, I->V1->coeffs + 0, fpctx);
        FLINT_ASSERT(fmpz_is_zero(temp));
#endif
        fmpz_mod_pow_fmpz(temp, roots + i, alphashift, fpctx);
        fmpz_mod_mul(S, S, temp, fpctx);
        fmpz_mod_inv(temp, S, fpctx);
        fmpz_mod_mul(Acoeff + Alen, V, temp, fpctx);
        if (fmpz_is_zero(Acoeff + Alen))
        {
            /* hmmm */
            continue;
        }

        if (fmpz_cmp(Acoeff + Alen, halfp) > 0)
        {
            fmpz_sub(Acoeff + Alen, Acoeff + Alen, fmpz_mod_ctx_modulus(fpctx));
        }


        mpoly_monomial_zero(Aexp + N*Alen, N);
        fmpz_mod_discrete_log_pohlig_hellman_run(new_exp, Ictx->dlogenv, roots + i);
        for (j = ctx->minfo->nvars - 1; j >= 0; j--)
        {
            this_exp = fmpz_fdiv_ui(new_exp, Ictx->subdegs[j]);
            fmpz_fdiv_q_ui(new_exp, new_exp, Ictx->subdegs[j]);
            if (this_exp >= Ictx->degbounds[j])
            {
                success = 0;
                goto cleanup;
            }
            (Aexp + N*Alen)[offsets[j]] |= this_exp << shifts[j];
        }
        if (!fmpz_is_zero(new_exp))
        {
            success = 0;
            goto cleanup;
        }
        Alen++;
    }
    A->length = Alen;

    fmpz_mpoly_sort_terms(A, ctx);

    success = 1;

cleanup:

    fmpz_clear(T);
    fmpz_clear(S);
    fmpz_clear(V);
    fmpz_clear(temp);
    fmpz_clear(halfp);

    TMP_END;
    return success;
}



void nmod_bma_mpoly_init(nmod_bma_mpoly_t A)
{
    A->length = 0;
    A->alloc = 0;
    A->exps = NULL;
    A->coeffs = NULL;
    A->pointcount = 0;
}

void nmod_bma_mpoly_reset_prime(
    nmod_bma_mpoly_t A,
    const nmodf_ctx_t fpctx)
{
    slong i;
    for (i = 0; i < A->alloc; i++)
    {
        nmod_berlekamp_massey_set_prime(A->coeffs + i, fpctx->mod.n);
    }
}


void nmod_bma_mpoly_clear(nmod_bma_mpoly_t A)
{
    slong i;
    for (i = 0; i < A->alloc; i++)
    {
        nmod_berlekamp_massey_clear(A->coeffs + i);
    }

    if (A->exps)
        flint_free(A->exps);
    if (A->coeffs)
        flint_free(A->coeffs);
}

void nmod_bma_mpoly_print(const nmod_bma_mpoly_t A)
{
    slong i;
    flint_printf("0");
    for (i = 0; i < A->length; i++)
    {
        flint_printf(" + [");
        nmod_berlekamp_massey_print(A->coeffs + i);
        flint_printf("]*X^%wd*Y^%wd",
                 A->exps[i] >> (FLINT_BITS/2),
                 A->exps[i] & ((-UWORD(1)) >> (FLINT_BITS/2)));
    }
}


void nmod_bma_mpoly_fit_length(
    nmod_bma_mpoly_t A,
    slong length,
    const nmodf_ctx_t fpctx)
{
    slong i;
    slong old_alloc = A->alloc;
    slong new_alloc = FLINT_MAX(length, A->alloc + A->alloc/2);

    if (length > old_alloc)
    {
        if (old_alloc == 0)
        {
            A->exps = (ulong *) flint_malloc(new_alloc*sizeof(ulong));
            A->coeffs = (nmod_berlekamp_massey_struct *) flint_malloc(
                               new_alloc*sizeof(nmod_berlekamp_massey_struct));
        }
        else
        {
            A->exps = (ulong *) flint_realloc(A->exps, new_alloc*sizeof(ulong));
            A->coeffs = (nmod_berlekamp_massey_struct *) flint_realloc(A->coeffs,
                               new_alloc*sizeof(nmod_berlekamp_massey_struct));
        }

        for (i = old_alloc; i < new_alloc; i++)
        {
            nmod_berlekamp_massey_init(A->coeffs + i, fpctx->mod.n);
        }
        A->alloc = new_alloc;
    }
}

void nmod_bma_mpoly_zero(nmod_bma_mpoly_t L)
{
    L->length = 0;
    L->pointcount = 0;
}


int nmod_bma_mpoly_reduce(nmod_bma_mpoly_t L)
{
    slong i;
    int changed;

    changed = 0;

    for (i = 0; i < L->length; i++)
    {
        FLINT_ASSERT(L->pointcount == nmod_berlekamp_massey_point_count(L->coeffs + i));
        changed |= nmod_berlekamp_massey_reduce(L->coeffs + i);
    }

    return changed;
}

void nmod_bma_mpoly_add_point(
    nmod_bma_mpoly_t L,
    const nmod_mpolyn_t A,
    const nmod_mpoly_ctx_t ctx_sp)
{
    slong j;
    slong Alen = A->length;
    nmod_poly_struct * Acoeff = A->coeffs;
    slong Li, Ai, ai;
    ulong Aexp;
    nmod_berlekamp_massey_struct * Lcoeff;
    slong Llen;
    ulong * Lexp;

    if (L->length == 0)
    {
        slong tot = 0;
        for (Ai = 0; Ai < Alen; Ai++)
        {
            for (ai = (Acoeff + Ai)->length - 1; ai >= 0; ai--)
            {
                tot += (0 != (Acoeff + Ai)->coeffs[ai]);
            }
        }
        nmod_bma_mpoly_fit_length(L, tot, ctx_sp->ffinfo);
    }

    Lcoeff = L->coeffs;
    Llen = L->length;
    Lexp = L->exps;

    Li = 0;
    Ai = ai = 0;
    Aexp = 0;
    if (Ai < Alen)
    {
        ai = nmod_poly_degree(A->coeffs + Ai);
        FLINT_ASSERT(FLINT_BIT_COUNT(ai) < FLINT_BITS/2);
        FLINT_ASSERT(0 == (A->exps[Ai] & LOW_HALF_MASK));
        Aexp = A->exps[Ai] + ai;
    }

    while (Li < Llen || Ai < Alen)
    {
        if (Li < Llen && Ai < Alen && Lexp[Li] == Aexp)
        {
            /* L term present, A term present */
add_same_exp:
            nmod_berlekamp_massey_add_point(Lcoeff + Li, (Acoeff + Ai)->coeffs[ai]);
            Li++;
            do {
                ai--;
            } while (ai >= 0 && (0 == (Acoeff + Ai)->coeffs[ai]));
            if (ai < 0)
            {
                Ai++;
                if (Ai < Alen)
                {
                    ai = nmod_poly_degree(A->coeffs + Ai);
                    FLINT_ASSERT(FLINT_BIT_COUNT(ai) < FLINT_BITS/2);
                    FLINT_ASSERT(0 == (A->exps[Ai] & LOW_HALF_MASK));
                    Aexp = A->exps[Ai] + ai;
                }
            }
            else
            {
                FLINT_ASSERT(Ai < Alen);
                FLINT_ASSERT(FLINT_BIT_COUNT(ai) < FLINT_BITS/2);
                FLINT_ASSERT(0 == (A->exps[Ai] & LOW_HALF_MASK));
                Aexp = A->exps[Ai] + ai;
            }
        }
        else if (Li < Llen && (Ai >= Alen || Lexp[Li] > Aexp))
        {
            /* L term present, A term missing */
            nmod_berlekamp_massey_add_zeros(Lcoeff + Li, 1);
            Li++;
        }
        else
        {
            /* L term missing, A term present */
            FLINT_ASSERT(Ai < Alen && (Li >= Llen || Lexp[Li] < Aexp));
            {
                ulong texp;
                nmod_berlekamp_massey_struct tcoeff;

                nmod_bma_mpoly_fit_length(L, Llen + 1, ctx_sp->ffinfo);
                Lcoeff = L->coeffs;
                Lexp = L->exps;

                texp = Lexp[Llen];
                tcoeff = Lcoeff[Llen];
                for (j = Llen - 1; j >= Li; j--)
                {
                    Lexp[j + 1] = Lexp[j];
                    Lcoeff[j + 1] = Lcoeff[j];
                }
                Lexp[Li] = texp;
                Lcoeff[Li] = tcoeff;
            }

            nmod_berlekamp_massey_start_over(Lcoeff + Li);
            nmod_berlekamp_massey_add_zeros(Lcoeff + Li, L->pointcount);
            Lexp[Li] = Aexp;
            Llen++;
            L->length = Llen;

            goto add_same_exp;
        }
    }

    L->pointcount++;
}



void new_nmod_bma_mpoly_add_point(
    nmod_bma_mpoly_t L,
    const n_bpoly_t A,
    const nmod_mpoly_ctx_t ctx_sp)
{
    slong j;
    slong Alen = A->length;
    n_poly_struct * Acoeffs = A->coeffs;
    slong Li, Ai, ai;
    nmod_berlekamp_massey_struct * Lcoeff;
    slong Llen;
    ulong * Lexp;

    if (L->length == 0)
    {
        slong tot = 0;
        for (Ai = 0; Ai < Alen; Ai++)
        for (ai = n_poly_degree(Acoeffs + Ai); ai >= 0; ai--)
            tot += (0 != Acoeffs[Ai].coeffs[ai]);
        nmod_bma_mpoly_fit_length(L, tot, ctx_sp->ffinfo);
    }

    Lcoeff = L->coeffs;
    Llen = L->length;
    Lexp = L->exps;

    Li = 0;
    Ai = A->length - 1;
    ai = 0;
    if (Ai >= 0)
        ai = n_poly_degree(Acoeffs + Ai);

    while (Li < Llen || Ai >= 0)
    {
        if (Li < Llen && Ai >= 0 && Lexp[Li] == pack_exp2(Ai, ai))
        {
            /* L term present, A term present */
add_same_exp:
            nmod_berlekamp_massey_add_point(Lcoeff + Li, Acoeffs[Ai].coeffs[ai]);
            Li++;

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
        else if (Li < Llen && (Ai < 0 || Lexp[Li] > pack_exp2(Ai, ai)))
        {
            /* L term present, A term missing */
            nmod_berlekamp_massey_add_zeros(Lcoeff + Li, 1);
            Li++;
        }
        else
        {
            /* L term missing, A term present */
            FLINT_ASSERT(Ai >= 0 && (Li >= Llen || Lexp[Li] < pack_exp2(Ai, ai)));
            {
                ulong texp;
                nmod_berlekamp_massey_struct tcoeff;

                nmod_bma_mpoly_fit_length(L, Llen + 1, ctx_sp->ffinfo);
                Lcoeff = L->coeffs;
                Lexp = L->exps;

                texp = Lexp[Llen];
                tcoeff = Lcoeff[Llen];
                for (j = Llen - 1; j >= Li; j--)
                {
                    Lexp[j + 1] = Lexp[j];
                    Lcoeff[j + 1] = Lcoeff[j];
                }
                Lexp[Li] = texp;
                Lcoeff[Li] = tcoeff;
            }

            nmod_berlekamp_massey_start_over(Lcoeff + Li);
            nmod_berlekamp_massey_add_zeros(Lcoeff + Li, L->pointcount);
            Lexp[Li] = pack_exp2(Ai, ai);
            Llen++;
            L->length = Llen;

            goto add_same_exp;
        }
    }

    L->pointcount++;
}



int new_nmod_zip_mpolyuu_add_point(
    nmod_zip_mpolyu_t L,
    const n_bpoly_t A)
{
    n_poly_struct * Acoeffs = A->coeffs;
    slong Li, Ai, ai;
    nmod_zip_struct * Lcoeff;
    slong Llen;
    ulong * Lexp;
    slong pointcount = L->pointcount;

    Lcoeff = L->coeffs;
    Llen = L->length;
    Lexp = L->exps;

    Li = 0;
    Ai = A->length - 1;
    ai = 0;
    if (Ai >= 0)
        ai = n_poly_degree(Acoeffs + Ai);

    for (Li = 0; Li < Llen; Li++)
    {
        nmod_zip_struct * Lc = Lcoeff + Li;

        if (Ai >= 0 && Lexp[Li] == pack_exp2(Ai, ai))
        {
            /* L present A present */
            FLINT_ASSERT(pointcount <= Lc->ealloc);
            Lc->evals[pointcount] = Acoeffs[Ai].coeffs[ai];
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
        else if (Ai < 0 || Lexp[Li] > pack_exp2(Ai, ai))
        {
            /* L present A missing */
            FLINT_ASSERT(pointcount <= Lc->ealloc);
            Lc->evals[pointcount] = 0;
        }
        else
        {
            /* L missing A present */
            return 0;
        }
    }

    L->pointcount = pointcount + 1;
    return 1;
}




int nmod_bma_mpoly_get_fmpz_mpolyu(
    fmpz_mpolyu_t A,
    const fmpz_mpoly_ctx_t ctx,
    ulong alphashift,
    const nmod_bma_mpoly_t L,
    const mpoly_bma_interpolate_ctx_t Ictx,
    const nmodf_ctx_t fpctx)
{
    int success;
    slong i;

    fmpz_mpolyu_fit_length(A, L->length, ctx);
    A->length = 0;
    for (i = 0; i < L->length; i++)
    {
        A->exps[A->length] = L->exps[i];
        success = nmod_mpoly_bma_get_fmpz_mpoly(A->coeffs + A->length, ctx,
                                       alphashift, L->coeffs + i, Ictx, fpctx);
        if (!success)
        {
            return 0;
        }
        A->length += !fmpz_mpoly_is_zero(A->coeffs + A->length, ctx);
    }
    return 1;
}



void fmpz_mod_bma_mpoly_init(fmpz_mod_bma_mpoly_t A)
{
    A->length = 0;
    A->alloc = 0;
    A->exps = NULL;
    A->coeffs = NULL;
    A->pointcount = 0;
}

void fmpz_mod_bma_mpoly_reset_prime(
    fmpz_mod_bma_mpoly_t A,
    const fmpz_mod_ctx_t fpctx)
{
    slong i;
    for (i = 0; i < A->alloc; i++)
    {
        fmpz_mod_berlekamp_massey_set_prime(A->coeffs + i, fmpz_mod_ctx_modulus(fpctx));
    }
}


void fmpz_mod_bma_mpoly_clear(fmpz_mod_bma_mpoly_t A)
{
    slong i;
    for (i = 0; i < A->alloc; i++)
    {
        fmpz_mod_berlekamp_massey_clear(A->coeffs + i);
    }

    if (A->exps)
        flint_free(A->exps);
    if (A->coeffs)
        flint_free(A->coeffs);
}

void fmpz_mod_bma_mpoly_print(
    fmpz_mod_bma_mpoly_t A,
    const mpoly_bma_interpolate_ctx_t Ictx)
{
    slong i;
    flint_printf("0");
    for (i = 0; i < A->length; i++)
    {
        flint_printf(" + [");
        fmpz_mod_berlekamp_massey_print(A->coeffs + i);
        flint_printf("]*X^%wd*Y^%wd",
                 A->exps[i] >> (FLINT_BITS/2),
                 A->exps[i] & ((-UWORD(1)) >> (FLINT_BITS - FLINT_BITS/2)));
    }
}


void fmpz_mod_bma_mpoly_fit_length(
    fmpz_mod_bma_mpoly_t A,
    slong length,
    const fmpz_mod_ctx_t fpctx)
{
    slong i;
    slong old_alloc = A->alloc;
    slong new_alloc = FLINT_MAX(length, A->alloc + A->alloc/2);

    if (length > old_alloc)
    {
        if (old_alloc == 0)
        {
            A->exps = (ulong *) flint_malloc(new_alloc*sizeof(ulong));
            A->coeffs = (fmpz_mod_berlekamp_massey_struct *) flint_malloc(
                           new_alloc*sizeof(fmpz_mod_berlekamp_massey_struct));
        }
        else
        {
            A->exps = (ulong *) flint_realloc(A->exps, new_alloc*sizeof(ulong));
            A->coeffs = (fmpz_mod_berlekamp_massey_struct *) flint_realloc(A->coeffs,
                           new_alloc*sizeof(fmpz_mod_berlekamp_massey_struct));
        }

        for (i = old_alloc; i < new_alloc; i++)
        {
            fmpz_mod_berlekamp_massey_init(A->coeffs + i, fmpz_mod_ctx_modulus(fpctx));
        }
        A->alloc = new_alloc;
    }
}

void fmpz_mod_bma_mpoly_zero(fmpz_mod_bma_mpoly_t L)
{
    L->length = 0;
    L->pointcount = 0;
}

int fmpz_mod_bma_mpoly_reduce(fmpz_mod_bma_mpoly_t L)
{
    slong i;
    int changed;

    changed = 0;

    for (i = 0; i < L->length; i++)
    {
        FLINT_ASSERT(L->pointcount == fmpz_mod_berlekamp_massey_point_count(L->coeffs + i));
        changed |= fmpz_mod_berlekamp_massey_reduce(L->coeffs + i);
    }

    return changed;
}

void fmpz_mod_bma_mpoly_add_point(
    fmpz_mod_bma_mpoly_t L,
    const fmpz_mod_mpolyn_t A,
    const fmpz_mod_mpoly_ctx_t ctx_mp)
{
    slong j;
    slong Alen = A->length;
    fmpz_mod_poly_struct * Acoeff = A->coeffs;
    slong Li, Ai, ai;
    ulong Aexp;
    fmpz_mod_berlekamp_massey_struct * Lcoeff;
    slong Llen;
    ulong * Lexp;

    if (L->length == 0)
    {
        slong tot = 0;
        for (Ai = 0; Ai < Alen; Ai++)
        {
            for (ai = (Acoeff + Ai)->length - 1; ai >= 0; ai--)
            {
                tot += !fmpz_is_zero((Acoeff + Ai)->coeffs + ai);
            }
        }
        fmpz_mod_bma_mpoly_fit_length(L, tot, ctx_mp->ffinfo);
    }

    Lcoeff = L->coeffs;
    Llen = L->length;
    Lexp = L->exps;

    Li = 0;
    Ai = ai = 0;
    Aexp = 0;
    if (Ai < Alen)
    {
        ai = fmpz_mod_poly_degree(A->coeffs + Ai);
        FLINT_ASSERT(FLINT_BIT_COUNT(ai) < FLINT_BITS/2);
        FLINT_ASSERT(0 == (A->exps[Ai] & LOW_HALF_MASK));
        Aexp = A->exps[Ai] + ai;
    }

    while (Li < Llen || Ai < Alen)
    {
        if (Li < Llen && Ai < Alen && Lexp[Li] == Aexp)
        {
            /* L term present, A term present */
add_same_exp:
            fmpz_mod_berlekamp_massey_add_point(Lcoeff + Li,
                                                   (Acoeff + Ai)->coeffs + ai);
            Li++;
            do {
                ai--;
            } while (ai >= 0 && fmpz_is_zero((Acoeff + Ai)->coeffs + ai));
            if (ai < 0)
            {
                Ai++;
                if (Ai < Alen)
                {
                    ai = fmpz_mod_poly_degree(A->coeffs + Ai);
                    FLINT_ASSERT(FLINT_BIT_COUNT(ai) < FLINT_BITS/2);
                    FLINT_ASSERT(0 == (A->exps[Ai] & LOW_HALF_MASK));
                    Aexp = A->exps[Ai] + ai;
                }
            }
            else
            {
                FLINT_ASSERT(Ai < Alen);
                FLINT_ASSERT(FLINT_BIT_COUNT(ai) < FLINT_BITS/2);
                FLINT_ASSERT(0 == (A->exps[Ai] & LOW_HALF_MASK));
                Aexp = A->exps[Ai] + ai;
            }
        }
        else if (Li < Llen && (Ai >= Alen || Lexp[Li] > Aexp))
        {
            /* L term present, A term missing */
            fmpz_mod_berlekamp_massey_add_zeros(Lcoeff + Li, 1);
            Li++;
        }
        else
        {
            /* L term missing, A term present */
            FLINT_ASSERT(Ai < Alen && (Li >= Llen || Lexp[Li] < Aexp));
            {
                ulong texp;
                fmpz_mod_berlekamp_massey_struct tcoeff;

                fmpz_mod_bma_mpoly_fit_length(L, Llen + 1, ctx_mp->ffinfo);
                Lcoeff = L->coeffs;
                Lexp = L->exps;

                texp = Lexp[Llen];
                tcoeff = Lcoeff[Llen];
                for (j = Llen - 1; j >= Li; j--)
                {
                    Lexp[j + 1] = Lexp[j];
                    Lcoeff[j + 1] = Lcoeff[j];
                }
                Lexp[Li] = texp;
                Lcoeff[Li] = tcoeff;
            }

            fmpz_mod_berlekamp_massey_start_over(Lcoeff + Li);
            fmpz_mod_berlekamp_massey_add_zeros(Lcoeff + Li, L->pointcount);
            Lexp[Li] = Aexp;
            Llen++;
            L->length = Llen;
            goto add_same_exp;
        }
    }

    L->pointcount++;
}

int fmpz_mod_bma_mpoly_get_fmpz_mpolyu(
    fmpz_mpolyu_t A,
    const fmpz_mpoly_ctx_t ctx,
    const fmpz_t alphashift,
    const fmpz_mod_bma_mpoly_t L,
    const mpoly_bma_interpolate_ctx_t Ictx,
    const fmpz_mod_ctx_t fpctx)
{
    int success;
    slong i;

    fmpz_mpolyu_fit_length(A, L->length, ctx);
    A->length = 0;
    for (i = 0; i < L->length; i++)
    {
        A->exps[A->length] = L->exps[i];
        success = fmpz_mod_bma_get_fmpz_mpoly(A->coeffs + A->length, ctx,
                                       alphashift, L->coeffs + i, Ictx, fpctx);
        if (!success)
        {
            return 0;
        }
        A->length += !fmpz_mpoly_is_zero(A->coeffs + A->length, ctx);
    }
    return 1;
}



ulong fmpz_mod_mpolyn_bidegree(const fmpz_mod_mpolyn_t A)
{
    ulong degx, degy;

    FLINT_ASSERT(A->length > 0);
    FLINT_ASSERT(A->bits == FLINT_BITS/2);

    degx = (A->exps + 1*0)[0] >> (FLINT_BITS/2);
    degy = fmpz_mod_poly_degree(A->coeffs + 0);

    FLINT_ASSERT(FLINT_BIT_COUNT(degx) < FLINT_BITS/2);
    FLINT_ASSERT(FLINT_BIT_COUNT(degy) < FLINT_BITS/2);

    return (degx << (FLINT_BITS/2)) + degy;
}

ulong nmod_mpolyn_bidegree(const nmod_mpolyn_t A)
{
    ulong degx, degy;

    FLINT_ASSERT(A->length > 0);
    FLINT_ASSERT(A->bits == FLINT_BITS/2);

    degx = A->exps[0] >> (FLINT_BITS/2);
    degy = nmod_poly_degree(A->coeffs + 0);

    FLINT_ASSERT(FLINT_BIT_COUNT(degx) < FLINT_BITS/2);
    FLINT_ASSERT(FLINT_BIT_COUNT(degy) < FLINT_BITS/2);

    return (degx << (FLINT_BITS/2)) + degy;
}


void fmpz_mpoly_eval_fmpz_mod(
    fmpz_t eval,
    const fmpz_mod_ctx_t fpctx,
    const fmpz_mpoly_t A,
    const fmpz * alpha,
    const fmpz_mpoly_ctx_t ctx)
{
    slong i, j;
    slong offset, shift;
    slong N = mpoly_words_per_exp_sp(A->bits, ctx->minfo);
    slong nvars = ctx->minfo->nvars;
    ulong * Aexp;
    slong * LUToffset;
    ulong * LUTmask;
    fmpz * LUTvalue;
    slong LUTlen;
    fmpz_t xpoweval;
    ulong * inputexpmask;
    TMP_INIT;

    TMP_START;

    fmpz_init(xpoweval);

    LUToffset = (slong *) TMP_ALLOC(N*FLINT_BITS*sizeof(slong));
    LUTmask   = (ulong *) TMP_ALLOC(N*FLINT_BITS*sizeof(ulong));
    LUTvalue  = (fmpz *) TMP_ALLOC(N*FLINT_BITS*sizeof(fmpz));
    for (i = 0; i < N*FLINT_BITS; i++)
    {
        fmpz_init(LUTvalue + i);
    }

    Aexp = A->exps;

    inputexpmask = (ulong *) TMP_ALLOC(N*sizeof(ulong));
    mpoly_monomial_zero(inputexpmask, N);
    for (i = 0; i < A->length; i++)
    {
        for (j = 0; j < N; j++)
        {
            inputexpmask[j] |= (Aexp + N*i)[j];
        }
    }

    LUTlen = 0;
    for (j = nvars - 1; j >= 0; j--)
    {
        mpoly_gen_offset_shift_sp(&offset, &shift, j, A->bits, ctx->minfo);

        fmpz_set(xpoweval, alpha + j); /* xpoweval = alpha[i]^(2^i) */
        for (i = 0; i < A->bits; i++)
        {
            LUToffset[LUTlen] = offset;
            LUTmask[LUTlen] = (UWORD(1) << (shift + i));
            fmpz_set(LUTvalue + LUTlen, xpoweval);
            if ((inputexpmask[offset] & LUTmask[LUTlen]) != 0)
            {
                LUTlen++;
            }
            fmpz_mod_mul(xpoweval, xpoweval, xpoweval, fpctx);
        }
    }
    FLINT_ASSERT(LUTlen < N*FLINT_BITS);

    fmpz_zero(eval);
    for (i = 0; i < A->length; i++)
    {
        fmpz_mod(xpoweval, A->coeffs + i, fmpz_mod_ctx_modulus(fpctx));
        for (j = 0; j < LUTlen; j++)
        {
            if (((Aexp + N*i)[LUToffset[j]] & LUTmask[j]) != 0)
            {
                fmpz_mod_mul(xpoweval, xpoweval, LUTvalue + j, fpctx);
            }
        }
        fmpz_mod_add(eval, eval, xpoweval, fpctx);
    }

    fmpz_clear(xpoweval);
    for (i = 0; i < N*FLINT_BITS; i++)
    {
        fmpz_clear(LUTvalue + i);
    }
    TMP_END;
}

/* take coeffs from [0, p) to (-p/2, p/2]
*/
void fmpz_mpolyu_symmetrize_coeffs(
    fmpz_mpolyu_t A,
    const fmpz_mpoly_ctx_t ctx,
    const fmpz_mod_ctx_t fpctx)
{
    slong i, j;
    fmpz_mpoly_struct * Ac;

    for (i = 0; i < A->length; i++)
    {
        Ac = A->coeffs + i;
        for (j = 0; j < Ac->length; j++)
        {
            fmpz_smod(Ac->coeffs + j, Ac->coeffs + j,
                                                  fmpz_mod_ctx_modulus(fpctx));
        }
    }
}



void fmpz_mpolyuu_eval_nmod(
    nmod_mpolyn_t E,
    const nmod_mpoly_ctx_t ctx_sp,
    const fmpz_mpolyu_t A,
    const mp_limb_t * alpha,
    const fmpz_mpoly_ctx_t ctx)
{
    slong i;
    slong xexp, yexp;
    mp_limb_t eval;

    FLINT_ASSERT(E->bits == FLINT_BITS/2);
    FLINT_ASSERT(1 == mpoly_words_per_exp_sp(E->bits, ctx_sp->minfo));

    E->length = 0;
    for (i = 0; i < A->length; i++)
    {
        eval = fmpz_mpoly_eval_nmod(ctx_sp->ffinfo, A->coeffs + i, alpha, ctx);
        if (eval == 0)
        {
            continue;
        }

        xexp = A->exps[i] >> (FLINT_BITS/2);
        yexp = A->exps[i] & (-UWORD(1) >> (FLINT_BITS - FLINT_BITS/2));

        if (E->length > 0 && (E->exps[E->length - 1] >> (FLINT_BITS/2)) == xexp)
        {
            nmod_poly_set_coeff_ui(E->coeffs + E->length - 1, yexp, eval);
        }
        else
        {
            nmod_mpolyn_fit_length(E, E->length + 1, ctx_sp);
            nmod_poly_zero(E->coeffs + E->length);
            nmod_poly_set_coeff_ui(E->coeffs + E->length, yexp, eval);
            E->exps[E->length] = xexp << (FLINT_BITS/2);
            E->length++;
        }
    }
}


void new_fmpz_mpolyuu_eval_nmod(
    n_bpoly_t E,
    const nmod_mpoly_ctx_t ctx_sp,
    const fmpz_mpolyu_t A,
    const mp_limb_t * alpha,
    const fmpz_mpoly_ctx_t ctx)
{
    slong i;
    slong xexp, yexp;
    mp_limb_t eval;

    E->length = 0;
    for (i = 0; i < A->length; i++)
    {
        eval = fmpz_mpoly_eval_nmod(ctx_sp->ffinfo, A->coeffs + i, alpha, ctx);
        if (eval == 0)
        {
            continue;
        }

        xexp = A->exps[i] >> (FLINT_BITS/2);
        yexp = A->exps[i] & (-UWORD(1) >> (FLINT_BITS - FLINT_BITS/2));

#if 0
        if (E->length > 0 && (E->exps[E->length - 1] >> (FLINT_BITS/2)) == xexp)
        {
            nmod_poly_set_coeff_ui(E->coeffs + E->length - 1, yexp, eval);
        }
        else
        {
            nmod_mpolyn_fit_length(E, E->length + 1, ctx_sp);
            nmod_poly_zero(E->coeffs + E->length);
            nmod_poly_set_coeff_ui(E->coeffs + E->length, yexp, eval);
            E->exps[E->length] = xexp << (FLINT_BITS/2);
            E->length++;
        }
#endif
        n_bpoly_set_coeff(E, xexp, yexp, eval);
    }
}


void fmpz_mpolyuu_eval_fmpz_mod(
    fmpz_mod_mpolyn_t E,
    const fmpz_mod_mpoly_ctx_t ctx_mp,
    const fmpz_mpolyu_t A,
    const fmpz * alpha,
    const fmpz_mpoly_ctx_t ctx)
{
    slong i;
    slong xexp, yexp;
    fmpz_t eval;

    FLINT_ASSERT(E->bits == FLINT_BITS/2);
    FLINT_ASSERT(1 == mpoly_words_per_exp_sp(E->bits, ctx_mp->minfo));

    fmpz_init(eval);

    E->length = 0;
    for (i = 0; i < A->length; i++)
    {
        fmpz_mpoly_eval_fmpz_mod(eval, ctx_mp->ffinfo, A->coeffs + i, alpha, ctx);
        if (fmpz_is_zero(eval))
        {
            continue;
        }

        xexp = A->exps[i] >> (FLINT_BITS/2);
        yexp = A->exps[i] & (-UWORD(1) >> (FLINT_BITS - FLINT_BITS/2));

        if (E->length > 0 && (E->exps[E->length - 1] >> (FLINT_BITS/2)) == xexp)
        {
            fmpz_mod_poly_set_coeff_fmpz(E->coeffs + E->length - 1, yexp, eval);
        }
        else
        {
            fmpz_mod_mpolyn_fit_length(E, E->length + 1, ctx_mp);
            fmpz_mod_poly_zero(E->coeffs + E->length);
            fmpz_mod_poly_set_coeff_fmpz(E->coeffs + E->length, yexp, eval);
            E->exps[E->length] = xexp << (FLINT_BITS/2);
            E->length++;
        }
    }

    fmpz_clear(eval);
}


mp_limb_t fmpz_mpoly_eval_nmod(
    const nmodf_ctx_t fpctx,
    const fmpz_mpoly_t A,
    const mp_limb_t * alpha,
    const fmpz_mpoly_ctx_t ctx)
{
    mp_limb_t eval;
    slong i, j;
    slong offset, shift;
    slong N = mpoly_words_per_exp_sp(A->bits, ctx->minfo);
    slong nvars = ctx->minfo->nvars;
    ulong * Aexp;
    slong * LUToffset;
    ulong * LUTmask;
    mp_limb_t * LUTvalue;
    slong LUTlen;
    mp_limb_t xpoweval;
    ulong * inputexpmask;
    TMP_INIT;

    FLINT_ASSERT(A->bits <= FLINT_BITS);

    TMP_START;

    LUToffset = (slong *) TMP_ALLOC(N*FLINT_BITS*sizeof(slong));
    LUTmask   = (ulong *) TMP_ALLOC(N*FLINT_BITS*sizeof(ulong));
    LUTvalue  = (mp_limb_t *) TMP_ALLOC(N*FLINT_BITS*sizeof(mp_limb_t));

    Aexp = A->exps;

    inputexpmask = (ulong *) TMP_ALLOC(N*sizeof(ulong));
    mpoly_monomial_zero(inputexpmask, N);
    for (i = 0; i < A->length; i++)
    {
        for (j = 0; j < N; j++)
        {
            inputexpmask[j] |= (Aexp + N*i)[j];
        }
    }

    LUTlen = 0;
    for (j = nvars - 1; j >= 0; j--)
    {
        mpoly_gen_offset_shift_sp(&offset, &shift, j, A->bits, ctx->minfo);

        xpoweval = alpha[j]; /* xpoweval = alpha[i]^(2^i) */
        for (i = 0; i < A->bits; i++)
        {
            LUToffset[LUTlen] = offset;
            LUTmask[LUTlen] = (UWORD(1) << (shift + i));
            LUTvalue[LUTlen] = xpoweval;
            if ((inputexpmask[offset] & LUTmask[LUTlen]) != 0)
            {
                LUTlen++;
            }
            xpoweval = nmod_mul(xpoweval, xpoweval, fpctx->mod);
        }
    }
    FLINT_ASSERT(LUTlen < N*FLINT_BITS);

    eval = 0;
    for (i = 0; i < A->length; i++)
    {
        xpoweval = fmpz_fdiv_ui(A->coeffs + i, fpctx->mod.n);
        for (j = 0; j < LUTlen; j++)
        {
            if (((Aexp + N*i)[LUToffset[j]] & LUTmask[j]) != 0)
            {
                xpoweval = nmod_mul(xpoweval, LUTvalue[j], fpctx->mod);
            }
        }
        eval = nmod_add(eval, xpoweval, fpctx->mod);
    }

    TMP_END;

    return eval;
}





void nmod_zip_init(nmod_zip_t Z)
{
    Z->mlength = 0;
    Z->malloc = 0;
    Z->coeffs = NULL;
    Z->monomials = NULL;
    Z->ealloc = 0;
    Z->evals = NULL;
}

void nmod_zip_clear(nmod_zip_t Z)
{
    if (Z->coeffs)
        flint_free(Z->coeffs);
    if (Z->monomials)
        flint_free(Z->monomials);
    if (Z->evals)
        flint_free(Z->evals);
}

void nmod_zip_set_lengths(nmod_zip_t A, slong mlength, slong elength)
{
    slong old_malloc = A->malloc;
    slong new_malloc = FLINT_MAX(mlength, A->malloc + A->malloc/2);
    slong old_ealloc = A->ealloc;
    slong new_ealloc = FLINT_MAX(elength, A->ealloc + A->ealloc/2);

    FLINT_ASSERT(mlength > 0);
    FLINT_ASSERT(elength > 0);

    if (mlength > old_malloc)
    {
        if (old_malloc == 0)
        {
            A->coeffs    = (mp_limb_t *) flint_malloc(
                                                 new_malloc*sizeof(mp_limb_t));
            A->monomials = (mp_limb_t *) flint_malloc(
                                                 new_malloc*sizeof(mp_limb_t));
        }
        else
        {
            A->coeffs    = (mp_limb_t *) flint_realloc(A->coeffs,
                                                 new_malloc*sizeof(mp_limb_t));
            A->monomials = (mp_limb_t *) flint_realloc(A->monomials,
                                                 new_malloc*sizeof(mp_limb_t));
        }
        A->malloc = new_malloc;
    }

    A->mlength = mlength;

    if (elength > old_ealloc)
    {
        if (old_ealloc == 0)
        {
            A->evals = (mp_limb_t *) flint_malloc(new_ealloc*sizeof(mp_limb_t));
        }
        else
        {
            A->evals = (mp_limb_t *) flint_realloc(A->evals,
                                                 new_ealloc*sizeof(mp_limb_t));
        }
        A->ealloc = new_ealloc;
    }
}

void nmod_zip_print(const nmod_zip_t Z, slong elength)
{
    slong i;

    printf("m ");
    for (i = 0; i < Z->mlength; i++)
    {
        flint_printf("(%wu %wu) ", Z->coeffs[i], Z->monomials[i]);
    }
    printf("e ");
    for (i = 0; i < elength; i++)
    {
        flint_printf("%wu ", Z->evals[i]);
    }
}



void nmod_zip_mpolyu_init(nmod_zip_mpolyu_t Z)
{
    Z->alloc = 0;
    Z->length = 0;
    Z->exps = NULL;
    Z->coeffs = NULL;
    Z->pointcount = 0;
}

void nmod_zip_mpolyu_clear(nmod_zip_mpolyu_t Z)
{
    slong i;

    for (i = 0; i < Z->alloc; i++)
    {
        nmod_zip_clear(Z->coeffs + i);
    }

    if (Z->exps)
        flint_free(Z->exps);
    if (Z->coeffs)
        flint_free(Z->coeffs);
}

void nmod_zip_mpolyu_fit_length(
    nmod_zip_mpolyu_t A,
    slong length)
{
    slong i;
    slong old_alloc = A->alloc;
    slong new_alloc = FLINT_MAX(length, A->alloc + A->alloc/2);

    if (length > old_alloc)
    {
        if (old_alloc == 0)
        {
            A->exps = (ulong *) flint_malloc(new_alloc*sizeof(ulong));
            A->coeffs = (nmod_zip_struct *) flint_malloc(
                                            new_alloc*sizeof(nmod_zip_struct));
        }
        else
        {
            A->exps = (ulong *) flint_realloc(A->exps, new_alloc*sizeof(ulong));
            A->coeffs = (nmod_zip_struct *) flint_realloc(A->coeffs,
                                            new_alloc*sizeof(nmod_zip_struct));
        }

        for (i = old_alloc; i < new_alloc; i++)
        {
            nmod_zip_init(A->coeffs + i);
        }
        A->alloc = new_alloc;
    }
}

void nmod_zip_mpolyu_fit_poly(
    nmod_zip_mpolyu_t Z,
    fmpz_mpolyu_t H,
    slong eval_length)
{
    slong i;
    FLINT_ASSERT(H->length > 0);

    nmod_zip_mpolyu_fit_length(Z, H->length);

    for (i = 0; i < H->length; i++)
    {
        Z->exps[i] = H->exps[i];
        nmod_zip_set_lengths(Z->coeffs + i, (H->coeffs + i)->length, eval_length);
    }

    Z->length = H->length;
    Z->pointcount = 0;
}


void nmod_mpoly_set_skel(
    nmod_mpolyc_t S,
    const nmod_mpoly_ctx_t ctx_sp,
    const fmpz_mpoly_t A,
    const mp_limb_t * alpha,
    const fmpz_mpoly_ctx_t ctx)
{
    slong i, j;
    slong offset, shift;
    slong N = mpoly_words_per_exp_sp(A->bits, ctx->minfo);
    slong nvars = ctx->minfo->nvars;
    ulong * Aexp;
    slong * LUToffset;
    ulong * LUTmask;
    mp_limb_t * LUTvalue;
    slong LUTlen;
    mp_limb_t xpoweval;
    ulong * inputexpmask;
    TMP_INIT;

    FLINT_ASSERT(A->bits <= FLINT_BITS);

    TMP_START;

    LUToffset = (slong *) TMP_ALLOC(N*FLINT_BITS*sizeof(slong));
    LUTmask   = (ulong *) TMP_ALLOC(N*FLINT_BITS*sizeof(ulong));
    LUTvalue  = (mp_limb_t *) TMP_ALLOC(N*FLINT_BITS*sizeof(mp_limb_t));

    Aexp = A->exps;

    inputexpmask = (ulong *) TMP_ALLOC(N*sizeof(ulong));
    mpoly_monomial_zero(inputexpmask, N);
    for (i = 0; i < A->length; i++)
    {
        for (j = 0; j < N; j++)
        {
            inputexpmask[j] |= (Aexp + N*i)[j];
        }
    }

    LUTlen = 0;
    for (j = nvars - 1; j >= 0; j--)
    {
        mpoly_gen_offset_shift_sp(&offset, &shift, j, A->bits, ctx->minfo);

        xpoweval = alpha[j]; /* xpoweval = alpha[i]^(2^i) */
        for (i = 0; i < A->bits; i++)
        {
            LUToffset[LUTlen] = offset;
            LUTmask[LUTlen] = (UWORD(1) << (shift + i));
            LUTvalue[LUTlen] = xpoweval;
            if ((inputexpmask[offset] & LUTmask[LUTlen]) != 0)
            {
                LUTlen++;
            }
            xpoweval = nmod_mul(xpoweval, xpoweval, ctx_sp->ffinfo->mod);
        }
    }
    FLINT_ASSERT(LUTlen < N*FLINT_BITS);

    nmod_mpolyc_fit_length(S, A->length);
    for (i = 0; i < A->length; i++)
    {
        xpoweval = 1;
        for (j = 0; j < LUTlen; j++)
        {
            if (((Aexp + N*i)[LUToffset[j]] & LUTmask[j]) != 0)
            {
                xpoweval = nmod_mul(xpoweval, LUTvalue[j], ctx_sp->ffinfo->mod);
            }
        }
        S->coeffs[i] = xpoweval;
    }
    S->length = A->length;

    TMP_END;
}

void nmod_zip_mpolyu_set_skel(
    nmod_zip_mpolyu_t Z,
    const nmod_mpoly_ctx_t ctx_sp,
    const fmpz_mpolyu_t A,
    const mp_limb_t * alpha,
    const fmpz_mpoly_ctx_t ctx)
{
    slong i, j;

    nmod_mpolyc_t T;
    nmod_mpolyc_init(T);

    FLINT_ASSERT(Z->length == A->length);
    for (i = 0; i < Z->length; i++)
    {
        nmod_zip_struct * Zc = Z->coeffs + i;

        nmod_mpoly_set_skel(T, ctx_sp, A->coeffs + i, alpha, ctx);

        Z->exps[i] = A->exps[i];
        FLINT_ASSERT(Zc->mlength == T->length);
        for (j = 0; j < Zc->mlength; j++)
        {
            Zc->coeffs[j] = 0;
            Zc->monomials[j] = T->coeffs[j];
        }
    }
    Z->pointcount = 0;

    nmod_mpolyc_clear(T);
}

void nmod_zip_mpolyuu_print(const nmod_zip_mpolyu_t A)
{
    slong i;
    flint_printf("0");
    for (i = 0; i < A->length; i++)
    {
        flint_printf(" + [");
        nmod_zip_print(A->coeffs + i, A->pointcount);
        flint_printf("]*X^%wd*Y^%wd",
                 A->exps[i] >> (FLINT_BITS/2),
                 A->exps[i] & ((-UWORD(1)) >> (FLINT_BITS/2)));
    }
}

int nmod_zip_mpolyuu_add_point(
    nmod_zip_mpolyu_t L,
    const nmod_mpolyn_t A)
{
    slong Alen = A->length;
    nmod_poly_struct * Acoeff = A->coeffs;
    slong Li, Ai, ai;
    ulong Aexp;
    nmod_zip_struct * Lcoeff;
    slong Llen;
    ulong * Lexp;
    slong pointcount = L->pointcount;

    Lcoeff = L->coeffs;
    Llen = L->length;
    Lexp = L->exps;

    Li = 0;
    Ai = ai = 0;
    Aexp = 0;
    if (Ai < Alen)
    {
        ai = nmod_poly_degree(A->coeffs + Ai);
        FLINT_ASSERT(FLINT_BIT_COUNT(ai) < FLINT_BITS/2);
        FLINT_ASSERT(0 == (A->exps[Ai] & LOW_HALF_MASK));
        Aexp = A->exps[Ai] + ai;
    }

    for (Li = 0; Li < Llen; Li++)
    {
        nmod_zip_struct * Lc = Lcoeff + Li;

        if (Ai < Alen && Lexp[Li] == Aexp)
        {
            /* L present A present */
            FLINT_ASSERT(pointcount <= Lc->ealloc);
            Lc->evals[pointcount] = (Acoeff + Ai)->coeffs[ai];
            do {
                ai--;
            } while (ai >= 0 && (0 == (Acoeff + Ai)->coeffs[ai]));
            if (ai < 0)
            {
                Ai++;
                if (Ai < Alen)
                {
                    ai = nmod_poly_degree(A->coeffs + Ai);        
                    FLINT_ASSERT(FLINT_BIT_COUNT(ai) < FLINT_BITS/2);
                    FLINT_ASSERT(0 == (A->exps[Ai] & LOW_HALF_MASK));
                    Aexp = A->exps[Ai] + ai;
                }
            }
            else
            {
                FLINT_ASSERT(Ai < Alen);
                FLINT_ASSERT(FLINT_BIT_COUNT(ai) < FLINT_BITS/2);
                FLINT_ASSERT(0 == (A->exps[Ai] & LOW_HALF_MASK));
                Aexp = A->exps[Ai] + ai;
            }
        }
        else if (Ai >= Alen || Lexp[Li] > Aexp)
        {
            /* L present A missing */
            FLINT_ASSERT(pointcount <= Lc->ealloc);
            Lc->evals[pointcount] = 0;
        }
        else
        {
            /* L missing A present */
            return 0;
        }
    }

    L->pointcount = pointcount + 1;
    return 1;
}



void nmod_mpolyu_set_skel(
    nmod_mpolycu_t S,
    const nmod_mpoly_ctx_t ctx_sp,
    const fmpz_mpolyu_t A,
    const mp_limb_t * alpha,
    const fmpz_mpoly_ctx_t ctx)
{
    slong i;
    nmod_mpolycu_fit_length(S, A->length);
    for (i = 0; i < A->length; i++)
    {
        nmod_mpoly_set_skel(S->coeffs + i, ctx_sp, A->coeffs + i, alpha, ctx);
    }
    S->length = A->length;
}



/* Ared = A mod p */
void nmod_mpoly_red_skel(
    nmod_mpolyc_t Ared,
    const fmpz_mpoly_t A,
    const nmodf_ctx_t fpctx)
{
    nmod_mpolyc_fit_length(Ared, A->length);
    Ared->length = A->length;
    _fmpz_vec_get_nmod_vec(Ared->coeffs, A->coeffs, A->length, fpctx->mod);
}

void nmod_mpolyu_red_skel(
    nmod_mpolycu_t Ared,
    const fmpz_mpolyu_t A,
    const nmodf_ctx_t fpctx)
{
    slong i;
    nmod_mpolycu_fit_length(Ared, A->length);
    Ared->length = A->length;
    for (i = 0; i < A->length; i++)
    {
        nmod_mpoly_red_skel(Ared->coeffs + i, A->coeffs + i, fpctx);
    }
}

/* M = S */
void nmod_mpoly_copy_skel(nmod_mpolyc_t M, const nmod_mpolyc_t S)
{
    slong i;
    nmod_mpolyc_fit_length(M, S->length);
    M->length = S->length;
    for (i = 0; i < S->length; i++)
    {
        M->coeffs[i] = S->coeffs[i];
    }
}

void nmod_mpolyu_copy_skel(nmod_mpolycu_t M, const nmod_mpolycu_t S)
{
    slong i;
    nmod_mpolycu_fit_length(M, S->length);
    M->length = S->length;
    for (i = 0; i < S->length; i++)
    {
        nmod_mpoly_copy_skel(M->coeffs + i, S->coeffs + i);
    }
}

/* return Ared.Acur and multiply Acur by Ainc */
mp_limb_t nmod_mpoly_use_skel_mul(
    const nmod_mpolyc_t Ared,
    nmod_mpolyc_t Acur,
    const nmod_mpolyc_t Ainc,
    const nmodf_ctx_t fpctx)
{
    slong i;
    mp_limb_t V, V0, V1, V2, p1, p0;
    FLINT_ASSERT(Ared->length == Acur->length);
    FLINT_ASSERT(Ared->length == Ainc->length);

    V0 = V1 = V2 = 0;
    for (i = 0; i < Ared->length; i++)
    {
        umul_ppmm(p1, p0, Ared->coeffs[i], Acur->coeffs[i]);
        add_sssaaaaaa(V2, V1, V0, V2, V1, V0, WORD(0), p1, p0);
        Acur->coeffs[i] = nmod_mul(Acur->coeffs[i], Ainc->coeffs[i], fpctx->mod);
    }
    
    NMOD_RED3(V, V2, V1, V0, fpctx->mod);
    return V;
}


/* M = S^k */
void nmod_mpoly_pow_skel(
    nmod_mpolyc_t M,
    const nmod_mpolyc_t S,
    ulong k,
    const nmod_mpoly_ctx_t ctx)
{
    slong i;
    nmod_mpolyc_fit_length(M, S->length);
    M->length = S->length;
    for (i = 0; i < S->length; i++)
    {
        M->coeffs[i] = nmod_pow_ui(S->coeffs[i], k, ctx->ffinfo->mod);
    }
}

void nmod_mpolyu_pow_skel(
    nmod_mpolycu_t M,
    const nmod_mpolycu_t S,
    ulong k,
    const nmod_mpoly_ctx_t ctx)
{
    slong i;
    nmod_mpolycu_fit_length(M, S->length);
    M->length = S->length;
    for (i = 0; i < S->length; i++)
    {
        nmod_mpoly_pow_skel(M->coeffs + i, S->coeffs + i, k, ctx);
    }
}





mp_limb_t n_poly_mod_eval_step2(
    n_poly_t Acur,
    const n_poly_t Ainc,
    nmod_t mod)
{
    slong i, Alen = Acur->length;
    mp_limb_t * cur = Acur->coeffs;
    const mp_limb_t * inc = Ainc->coeffs;
    ulong t0, t1, t2, p0, p1;

    FLINT_ASSERT(2*Alen == Ainc->length);

    t2 = t1 = t0 = 0;
    for (i = 0; i < Alen; i++)
    {
        umul_ppmm(p1, p0, cur[i], inc[2*i + 0]);
        add_sssaaaaaa(t2, t1, t0, t2, t1, t0, 0, p1, p0);
        cur[i] = nmod_mul(cur[i], inc[2*i + 1], mod);
    }
    NMOD_RED3(t0, t2, t1, t0, mod);
    return t0;
}

void nmod_mpolyuu_eval_step2(
    nmod_mpolyn_t E,
    n_bpoly_t Acur,
    const n_polyun_t Ainc,
    const nmod_mpoly_ctx_t ctx_sp)
{
    slong xexp, yexp;
    slong i;
    slong Alen = Acur->length;
    n_polyun_term_struct * Aterms = Ainc->terms;
    mp_limb_t eval;

    FLINT_ASSERT(Alen == Ainc->length);
    FLINT_ASSERT(E->bits == FLINT_BITS/2);
    FLINT_ASSERT(1 == mpoly_words_per_exp_sp(E->bits, ctx_sp->minfo));

    E->length = 0;
    for (i = 0; i < Alen; i++)
    {
        eval = n_poly_mod_eval_step2(Acur->coeffs + i, Aterms[i].coeff,
                                                          ctx_sp->ffinfo->mod);
        if (eval == 0)
        {
            continue;
        }

        xexp = Aterms[i].exp >> (FLINT_BITS/2);
        yexp = Aterms[i].exp & LOW_HALF_MASK;

        if (E->length > 0 && (E->exps[E->length - 1] >> (FLINT_BITS/2)) == xexp)
        {
            nmod_poly_set_coeff_ui(E->coeffs + E->length - 1, yexp, eval);
        }
        else
        {
            nmod_mpolyn_fit_length(E, E->length + 1, ctx_sp);
            nmod_poly_zero(E->coeffs + E->length);
            nmod_poly_set_coeff_ui(E->coeffs + E->length, yexp, eval);
            E->exps[E->length] = xexp << (FLINT_BITS/2);
            E->length++;
        }
    }
}


void new_nmod_mpolyuu_eval_step2(
    n_bpoly_t E,
    n_bpoly_t Acur,
    const n_polyun_t Ainc,
    const nmod_mpoly_ctx_t ctx_sp)
{
    slong xexp, yexp;
    slong i;
    slong Alen = Acur->length;
    n_polyun_term_struct * Aterms = Ainc->terms;
    mp_limb_t eval;

    FLINT_ASSERT(Alen == Ainc->length);

    E->length = 0;
    for (i = 0; i < Alen; i++)
    {
        eval = n_poly_mod_eval_step2(Acur->coeffs + i, Aterms[i].coeff,
                                                          ctx_sp->ffinfo->mod);
        if (eval == 0)
        {
            continue;
        }

        xexp = Aterms[i].exp >> (FLINT_BITS/2);
        yexp = Aterms[i].exp & LOW_HALF_MASK;

        n_bpoly_set_coeff(E, xexp, yexp, eval);
    }

    FLINT_ASSERT(n_bpoly_mod_is_canonical(E, ctx_sp->ffinfo->mod));
}



void n_bpoly_set_mpolyn2(
    n_bpoly_t A,
    const nmod_mpolyn_t B,
    const nmod_mpoly_ctx_t ctx_sp)
{
    slong i, x;

    FLINT_ASSERT(B->bits == FLINT_BITS/2);
    FLINT_ASSERT(1 == mpoly_words_per_exp_sp(B->bits, ctx_sp->minfo));

    A->length = 0;

    for (i = 0; i < B->length; i++)
    {
        x = B->exps[i] >> (FLINT_BITS/2);
        n_bpoly_fit_length(A, x + 1);
        while (A->length <= x)
        {
            n_poly_zero(A->coeffs + A->length);
            A->length++;
        }

        n_poly_set_nmod_poly(A->coeffs + x, B->coeffs + i);
    }

    FLINT_ASSERT(n_bpoly_mod_is_canonical(A, ctx_sp->ffinfo->mod));
}


void nmod_mpolyuu_use_skel_mul(
    nmod_mpolyn_t E,
    const fmpz_mpolyu_t A,
    const nmod_mpolycu_t Ared,
    nmod_mpolycu_t Acur,
    const nmod_mpolycu_t Ainc,
    const nmod_mpoly_ctx_t ctx_sp)
{
    slong xexp, yexp;
    slong i;
    mp_limb_t eval;

    FLINT_ASSERT(E->bits == FLINT_BITS/2);
    FLINT_ASSERT(1 == mpoly_words_per_exp_sp(E->bits, ctx_sp->minfo));

    FLINT_ASSERT(A->length == Ared->length);
    FLINT_ASSERT(A->length == Acur->length);
    FLINT_ASSERT(A->length == Ainc->length);

    E->length = 0;
    for (i = 0; i < A->length; i++)
    {
        eval = nmod_mpoly_use_skel_mul(Ared->coeffs + i, Acur->coeffs + i,
                                             Ainc->coeffs + i, ctx_sp->ffinfo);
        if (eval == 0)
        {
            continue;
        }

        xexp = A->exps[i] >> (FLINT_BITS/2);
        yexp = A->exps[i] & ((-UWORD(1)) >> (FLINT_BITS - FLINT_BITS/2));

        if (E->length > 0 && (E->exps[E->length - 1] >> (FLINT_BITS/2)) == xexp)
        {
            nmod_poly_set_coeff_ui(E->coeffs + E->length - 1, yexp, eval);
        }
        else
        {
            nmod_mpolyn_fit_length(E, E->length + 1, ctx_sp);
            nmod_poly_zero(E->coeffs + E->length);
            nmod_poly_set_coeff_ui(E->coeffs + E->length, yexp, eval);
            E->exps[E->length] = xexp << (FLINT_BITS/2);
            E->length++;
        }
    }
}



nmod_zip_find_coeffs_ret_t nmod_zip_find_coeffs(
    nmod_zip_t Z,
    nmod_poly_t master,
    slong pointcount,
    const nmodf_ctx_t ffinfo)
{
    slong i, j;
    mp_limb_t V, V0, V1, V2, T, S, r, p0, p1;

    FLINT_ASSERT(pointcount >= Z->mlength);

    nmod_poly_product_roots_nmod_vec(master, Z->monomials, Z->mlength);

    for (i = 0; i < Z->mlength; i++)
    {
        /* coeffs[i] is (coeffs(P).values)/P(roots[i]) =: V/S
            where P(x) = master(x)/(x-roots[i])     */
        V0 = V1 = V2 = T = S = 0;
        r = Z->monomials[i];
        for (j = Z->mlength; j > 0; j--)
        {
            T = nmod_add(nmod_mul(r, T, ffinfo->mod), master->coeffs[j], ffinfo->mod);
            S = nmod_add(nmod_mul(r, S, ffinfo->mod), T, ffinfo->mod);
            umul_ppmm(p1, p0, Z->evals[j - 1], T);
            add_sssaaaaaa(V2, V1, V0, V2, V1, V0, WORD(0), p1, p0);
        }
        /* roots[i] should be a root of master */
        FLINT_ASSERT(nmod_add(nmod_mul(r, T, ffinfo->mod), master->coeffs[0], ffinfo->mod) == 0);
        NMOD_RED3(V, V2, V1, V0, ffinfo->mod);
        S = nmod_mul(S, r, ffinfo->mod); /* shift is one */
        if (S == 0)
        {
            return nmod_zip_find_coeffs_non_invertible;
        }
        Z->coeffs[i] = nmod_mul(V, nmod_inv(S, ffinfo->mod), ffinfo->mod);
    }

    /* use the coefficients of master as temp work space */
    for (i = 0; i < Z->mlength; i++)
    {
        master->coeffs[i] = nmod_pow_ui(Z->monomials[i], Z->mlength, ffinfo->mod);
    }

    /* check that the remaining points match */
    for (i = Z->mlength; i < pointcount; i++)
    {
        V0 = V1 = V2 = S = 0;
        for (j = 0; j < Z->mlength; j++)
        {
            master->coeffs[j] = nmod_mul(master->coeffs[j], Z->monomials[j], ffinfo->mod);
            umul_ppmm(p1, p0, Z->coeffs[j], master->coeffs[j]);
            add_sssaaaaaa(V2, V1, V0, V2, V1, V0, WORD(0), p1, p0);
        }
        NMOD_RED3(V, V2, V1, V0, ffinfo->mod);
        if (V != Z->evals[i])
        {
            return nmod_zip_find_coeffs_no_match;
        }
    }

    return nmod_zip_find_coeffs_good;
}

nmod_zip_find_coeffs_ret_t nmod_mpolyu_zip_find_coeffs(
    nmod_zip_mpolyu_t Z,
    const nmod_mpoly_ctx_t ctx_sp)
{
    slong i;
    nmod_zip_find_coeffs_ret_t ret;
    nmod_poly_t T;

    nmod_poly_init_mod(T, ctx_sp->ffinfo->mod);

    for (i = 0; i < Z->length; i++)
    {
        ret = nmod_zip_find_coeffs(Z->coeffs + i, T, Z->pointcount, ctx_sp->ffinfo);
        if (ret != nmod_zip_find_coeffs_good)
        {
            goto cleanup;
        }
    }

    ret = nmod_zip_find_coeffs_good;

cleanup:

    nmod_poly_clear(T);

    return ret;
}


int fmpz_mpolyu_addinterp_zip(
    fmpz_mpolyu_t H,
    const fmpz_t Hmodulus,
    const nmod_zip_mpolyu_t Z,
    const nmodf_ctx_t ffinfo)
{
    int changed = 0;
    slong i, j;
    fmpz_t t;

    fmpz_init(t);

    FLINT_ASSERT(H->length == Z->length);
    for (i = 0; i < H->length; i++)
    {
        fmpz_mpoly_struct * Hc = H->coeffs + i;
        nmod_zip_struct * Zc = Z->coeffs + i;

        FLINT_ASSERT(Hc->length == Zc->mlength);
        for (j = 0; j < Hc->length; j++)
        {
            fmpz_CRT_ui(t, Hc->coeffs + j, Hmodulus, Zc->coeffs[j], ffinfo->mod.n, 1);
            changed |= !fmpz_equal(t, Hc->coeffs + j);
            fmpz_swap(t, Hc->coeffs + j);
        }
    }

    fmpz_clear(t);
    return changed;
}


int fmpz_mpolyu_repack_bits(
    fmpz_mpolyu_t A,
    flint_bitcnt_t Abits,
    const fmpz_mpoly_ctx_t ctx)
{
    flint_bitcnt_t org_bits = A->bits;
    int success;
    slong i, j;

    for (i = 0; i < A->length; i++)
    {
        FLINT_ASSERT((A->coeffs + i)->bits == A->bits);
        success = fmpz_mpoly_repack_bits_inplace(A->coeffs + i, Abits, ctx);
        if (!success)
        {
            /* repack changed coeffs */
            for (j = 0; j < i; j++)
            {
                success = fmpz_mpoly_repack_bits_inplace(A->coeffs + j, org_bits, ctx);
                FLINT_ASSERT(success);
            }
            return 0;
        }
    }
    return 1;
}


/* fit bits is not the best fxn to accomplish this */
void fmpz_mpoly_set_bits(
    fmpz_mpoly_t A,
    flint_bitcnt_t Abits,
    const fmpz_mpoly_ctx_t ctx)
{
    fmpz_mpoly_fit_bits(A, Abits, ctx);
    A->bits = Abits;
}

void fmpz_mpolyu_set_bits(
    fmpz_mpolyu_t A,
    flint_bitcnt_t Abits,
    const fmpz_mpoly_ctx_t ctx)
{
    slong i;
    for (i = 0; i < A->length; i++)
    {
        fmpz_mpoly_set_bits(A->coeffs + i, Abits, ctx);
    }
    A->bits = Abits;
}


/*
    A is in ZZ[x_0,..., x_(n-1)]
    evaluation at x_i = values[i]
    out is in Fp[x_var]   (p is stored in out :( )

    out_deg is deg of out.

    This function fails to produce out if it will have too large degree:
    If return is 0, then out is invalid. Otherwise, out is valid.
    out_deg is always valid.
*/
int fmpz_mpoly_eval_all_but_one_nmod(
    slong * out_deg,
    nmod_poly_t out,
    const fmpz_mpoly_t A,
    slong var,
    mp_limb_t * values,
    const fmpz_mpoly_ctx_t ctx)
{
    slong i, j;
    slong deg;
    const slong deg_limit = 9999;
    ulong varexp, thisexp;
    mp_limb_t t, v;
    ulong mask;
    slong * offsets, * shifts;
    slong N = mpoly_words_per_exp_sp(A->bits, ctx->minfo);
    ulong * Aexp = A->exps;
    fmpz * Acoeff = A->coeffs;
    TMP_INIT;

    FLINT_ASSERT(A->bits <= FLINT_BITS);

    TMP_START;

    mask = (-UWORD(1)) >> (FLINT_BITS - A->bits);
    offsets = (slong *) TMP_ALLOC(ctx->minfo->nvars*sizeof(slong));
    shifts = (slong *) TMP_ALLOC(ctx->minfo->nvars*sizeof(slong));
    for (j = 0; j < ctx->minfo->nvars; j++)
    {
        mpoly_gen_offset_shift_sp(offsets + j, shifts + j, j, A->bits, ctx->minfo);
    }

    nmod_poly_zero(out);
    deg = -WORD(1);
    for (i = 0; i < A->length; i++)
    {
        varexp = ((Aexp + N*i)[offsets[var]]>>shifts[var])&mask;
        deg = FLINT_MAX(deg, (slong)(varexp));
        v = fmpz_fdiv_ui(Acoeff + i, out->mod.n);
        for (j = 0; j < ctx->minfo->nvars; j++)
        {
            thisexp = ((Aexp + N*i)[offsets[j]]>>shifts[j])&mask;
            if (j != var)
            {
                v = nmod_mul(v, nmod_pow_ui(values[j], thisexp, out->mod), out->mod);
            }
        }
        t = nmod_poly_get_coeff_ui(out, varexp);
        if (deg <= deg_limit)
        {
            nmod_poly_set_coeff_ui(out, varexp, nmod_add(t, v, out->mod));
        }
    }

    TMP_END;

    *out_deg = deg;
    return deg <= deg_limit;
}

/*
    A is in ZZ[X,Y][x_0,..., x_(n-1)]
    out is in Fp[x_var]   (p is stored in out :( )

    X = values[0]
    Y = values[1]
    x_0 = values[2]
    ...
*/
int fmpz_mpolyuu_eval_all_but_one_nmod(
    slong * out_deg,
    nmod_poly_t out,
    const fmpz_mpolyu_t A,
    slong var,
    mp_limb_t * values,
    const fmpz_mpoly_ctx_t ctx)
{
    slong i;
    slong deg, this_deg;
    ulong * Aexp = A->exps;
    const fmpz_mpoly_struct * Acoeff = A->coeffs;
    nmod_poly_t temp;
    mp_limb_t t, t1;
    int success, this_success;

    FLINT_ASSERT(A->bits <= FLINT_BITS);

    nmod_poly_zero(out);
    nmod_poly_init(temp, out->mod.n);

    deg = -WORD(1);
    success = 1;
    for (i = 0; i < A->length; i++)
    {
        t = nmod_pow_ui(values[0], Aexp[i] >> (FLINT_BITS/2), out->mod);
        t1 = nmod_pow_ui(values[1], Aexp[i] & LOW_HALF_MASK, out->mod);
        t = nmod_mul(t, t1, out->mod);
        this_success = fmpz_mpoly_eval_all_but_one_nmod(&this_deg, temp, Acoeff + i, var, values + 2, ctx);
        deg = FLINT_MAX(deg, this_deg);
        success = success && this_success;
        if (success)
        {
            nmod_poly_scalar_mul_nmod(temp, temp, t);
            nmod_poly_add(out, out, temp);
        }
    }

    nmod_poly_clear(temp);  

    *out_deg = deg;
 
    return success; 
}

/*
    A, B are in ZZ[X,Y][x_0, ..., x_(n-1)]

    Set Adeg (resp Bdeg) to the degree of A (resp B) wrt x_var.
    Return an upper bound on the degree of gcd(A,B) wrt x_var.
*/
slong fmpz_mpolyuu_gcd_degree_bound_minor(
    slong * Adeg,
    slong * Bdeg,
    const fmpz_mpolyu_t A,
    const fmpz_mpolyu_t B,
    slong var,
    const fmpz_mpoly_ctx_t ctx,
    flint_rand_t state)
{
    slong i;
    int tries = 0;
    slong degA, degB, degRet;
    mp_limb_t p = UWORD(1) << (FLINT_BITS - 2);
    nmod_poly_t Geval, Aeval, Beval;
    mp_limb_t * values;
    int Asuccess, Bsuccess;
    TMP_INIT;

    TMP_START;

    values = (mp_limb_t *) TMP_ALLOC((ctx->minfo->nvars + 2)*sizeof(mp_limb_t));

    p = n_nextprime(p, 1);
    nmod_poly_init(Geval, p);
    nmod_poly_init(Aeval, p);
    nmod_poly_init(Beval, p);

try_again:

    for (i = 0; i < ctx->minfo->nvars + 2; i++)
    {
        values[i] = n_urandint(state, p);
    }

    Asuccess = fmpz_mpolyuu_eval_all_but_one_nmod(&degA, Aeval, A, var, values, ctx);
    Bsuccess = fmpz_mpolyuu_eval_all_but_one_nmod(&degB, Beval, B, var, values, ctx);
    *Adeg = degA;
    *Bdeg = degB;

    if (!Asuccess || !Bsuccess)
    {
        degRet = FLINT_MIN(degA, degB);
        goto cleanup;
    }

    if (degA != nmod_poly_degree(Aeval) || degB != nmod_poly_degree(Beval))
    {
        if (++tries > 100)
        {
            degRet = FLINT_MIN(degA, degB);
            goto cleanup;
        }
        p = n_nextprime(p, 1);
        nmod_poly_clear(Geval);
        nmod_poly_clear(Aeval);
        nmod_poly_clear(Beval);
        nmod_poly_init(Geval, p);
        nmod_poly_init(Aeval, p);
        nmod_poly_init(Beval, p);
        goto try_again;
    }

    nmod_poly_gcd(Geval, Aeval, Beval);
    degRet = nmod_poly_degree(Geval);

cleanup:

    nmod_poly_clear(Geval);
    nmod_poly_clear(Aeval);
    nmod_poly_clear(Beval);
    TMP_END;
    return degRet;
}

/*
    A is in ZZ[x_0, ..., x_(n-1)]

    After the substitutions
        x_0     = x ^ (sub[1] * sub[2] * ... * sub[n-1])

        x_(n-2) = x ^ (sub[n-1])
        x_(n-1) = x ^ (1)
    a univariate in ZZ[x] remains. Return the content of this poly.
*/
void fmpz_mpoly_ksub_content(
    fmpz_t content,
    const fmpz_mpoly_t A,
    const ulong * subdegs,
    const fmpz_mpoly_ctx_t ctx)
{
    slong i, j;
    slong N = mpoly_words_per_exp_sp(A->bits, ctx->minfo);
    fmpz_mpoly_t T;
    fmpz_mpoly_ctx_t Tctx;
    fmpz_t e;
    fmpz * Acoeff = A->coeffs;
    ulong * Aexp = A->exps;
    ulong mask;
    slong * offsets, * shifts;
    TMP_INIT;

    TMP_START;
    fmpz_init(e);

    fmpz_mpoly_ctx_init(Tctx, 1, ORD_LEX);
    fmpz_mpoly_init(T, Tctx);

    mask = (-UWORD(1)) >> (FLINT_BITS - A->bits);
    offsets = (slong *) TMP_ALLOC(ctx->minfo->nvars*sizeof(slong));
    shifts = (slong *) TMP_ALLOC(ctx->minfo->nvars*sizeof(slong));
    for (j = 0; j < ctx->minfo->nvars; j++)
    {
        mpoly_gen_offset_shift_sp(offsets + j, shifts + j, j, A->bits, ctx->minfo);
    }

    for (i = 0; i < A->length; i++)
    {
        fmpz_zero(e);
        for (j = 0; j < ctx->minfo->nvars; j++)
        {
            fmpz_mul_ui(e, e, subdegs[j]);
            fmpz_add_ui(e, e, ((Aexp + N*i)[offsets[j]]>>shifts[j])&mask);
        }
        _fmpz_mpoly_push_exp_ffmpz(T, e, Tctx);
        fmpz_set(T->coeffs + T->length - 1, Acoeff + i);
    }
    fmpz_mpoly_sort_terms(T, Tctx);
    fmpz_mpoly_combine_like_terms(T, Tctx);
    _fmpz_vec_content(content, T->coeffs, T->length);
    fmpz_mpoly_clear(T, Tctx);
    fmpz_mpoly_ctx_clear(Tctx);

    fmpz_clear(e);
    TMP_END;
}

/*
    Check if H evaluated a random point (eval(H)) matches gcd(eval(A), eval(B))
    with leading coefficient eval(Gamma). The evaluations are in x_i.

    H, A, B are in Z[X,Y][x_0, ..., x_(n-1)]
    Gamma is in Z[x_0, ..., x_(n-1)]
*/
typedef enum {
    random_check_good,
    random_check_point_not_found,
    random_check_image_degree_high,
    random_check_image_degree_low,
    random_check_image_no_match    
} random_check_ret_t;

random_check_ret_t static _random_check_sp(
    ulong * GevaldegXY,
    ulong GdegboundXY,
    int which_check,
    n_bpoly_t new_Aeval_sp,
    n_bpoly_t new_Beval_sp,
    n_bpoly_t new_Geval_sp,
    n_bpoly_t new_Abareval_sp,
    n_bpoly_t new_Bbareval_sp,
    mp_limb_t * checkalpha_sp,
    const fmpz_mpolyu_t H,
    const fmpz_mpolyu_t A,
    const fmpz_mpolyu_t B,
    const fmpz_mpoly_t Gamma,
    const fmpz_mpoly_ctx_t ctx,
    const nmod_mpoly_ctx_t ctx_sp,
    flint_rand_t randstate,
    n_poly_bpoly_stack_t new_Sp)
{
    mp_limb_t Gammaeval_sp;
    int success;
    int point_try_count;
    slong i;

    /* try to test H at a random evaluation point */
    for (point_try_count = 0; point_try_count < 10; point_try_count++)
    {
        /* evaluate Gamma, A, and B at random point */
        for (i = 0; i < ctx->minfo->nvars; i++)
            checkalpha_sp[i] = n_urandint(randstate, ctx_sp->ffinfo->mod.n);
        new_fmpz_mpolyuu_eval_nmod(new_Aeval_sp, ctx_sp, A, checkalpha_sp, ctx);
        new_fmpz_mpolyuu_eval_nmod(new_Beval_sp, ctx_sp, B, checkalpha_sp, ctx);

        /* make sure that evaluation did not kill either lc(A) or lc(B) */
        if (new_Aeval_sp->length < 1 || new_Beval_sp->length < 1 ||
            n_bpoly_bidegree(new_Aeval_sp) != A->exps[0] ||
            n_bpoly_bidegree(new_Beval_sp) != B->exps[0])
        {
            continue;
        }

        /* Gamma is gcd(lc(A), lc(B)) so it evaluation should not be zero */
        Gammaeval_sp = fmpz_mpoly_eval_nmod(ctx_sp->ffinfo, Gamma,
                                                           checkalpha_sp, ctx);
        FLINT_ASSERT(Gammaeval_sp != 0);

        success = n_bpoly_mod_gcd_brown_smprime(new_Geval_sp,
                  new_Abareval_sp, new_Bbareval_sp, new_Aeval_sp, new_Beval_sp, ctx_sp->ffinfo->mod, new_Sp);
        if (!success)
            continue;

        n_bpoly_scalar_mul_nmod(new_Geval_sp, Gammaeval_sp, ctx_sp->ffinfo->mod);

        FLINT_ASSERT(new_Geval_sp->length > 0);
        *GevaldegXY = n_bpoly_bidegree(new_Geval_sp);

        if (GdegboundXY < *GevaldegXY)
        {
            return random_check_image_degree_high;
        }
        else if (GdegboundXY > *GevaldegXY)
        {
            return random_check_image_degree_low;
        }

        if (which_check == 0)
        {
            new_fmpz_mpolyuu_eval_nmod(new_Bbareval_sp, ctx_sp, H, checkalpha_sp, ctx);
            if (!n_bpoly_equal(new_Bbareval_sp, new_Geval_sp))
                return random_check_image_no_match;
        }
        else
        {
            new_fmpz_mpolyuu_eval_nmod(new_Geval_sp, ctx_sp, H, checkalpha_sp, ctx);

            if (!n_bpoly_equal(new_Geval_sp,
                         which_check == 1 ? new_Abareval_sp : new_Bbareval_sp))
            {
                return random_check_image_no_match;
            }
        }

        return random_check_good;
    }

    return random_check_point_not_found;
}

random_check_ret_t static _random_check(
    ulong * GevaldegXY,
    ulong GdegboundXY,
    fmpz_mod_mpolyn_t Aeval,
    fmpz_mod_mpolyn_t Beval,
    fmpz_mod_mpolyn_t Geval,
    fmpz_mod_mpolyn_t Abareval,
    fmpz_mod_mpolyn_t Bbareval,
    fmpz_t Gammaeval,
    fmpz * checkalpha,
    const fmpz_mpolyu_t H,
    const fmpz_mpolyu_t A,
    const fmpz_mpolyu_t B,
    const fmpz_mpoly_t Gamma,
    const fmpz_mpoly_ctx_t ctx,
    const fmpz_mod_mpoly_ctx_t ctx_mp,
    flint_rand_t randstate)
{
    int success;
    int point_try_count;
    slong i;

    /* try to test H at a random evaluation point */
    for (point_try_count = 0; point_try_count < 10; point_try_count++)
    {
        /* evaluate Gamma, A, and B at random point */
        for (i = 0; i < ctx->minfo->nvars; i++)
        {
            fmpz_randm(checkalpha + i, randstate, fmpz_mod_ctx_modulus(ctx_mp->ffinfo));
        }
        fmpz_mpolyuu_eval_fmpz_mod(Aeval, ctx_mp, A, checkalpha, ctx);
        fmpz_mpolyuu_eval_fmpz_mod(Beval, ctx_mp, B, checkalpha, ctx);

        /* make sure that evaluation did not kill either lc(A) or lc(B) */
        if ( Aeval->length == 0 || Beval->length == 0 
            || fmpz_mod_mpolyn_bidegree(Aeval) != A->exps[0]
            || fmpz_mod_mpolyn_bidegree(Beval) != B->exps[0])
        {
            continue;
        }

        /* Gamma is gcd(lc(A), lc(B)) so it evaluation should not be zero */
        fmpz_mpoly_eval_fmpz_mod(Gammaeval, ctx_mp->ffinfo, Gamma, checkalpha, ctx);
        FLINT_ASSERT(!fmpz_is_zero(Gammaeval));

        success = fmpz_mod_mpolyn_gcd_brown_bivar(Geval, Abareval, Bbareval,
                                                         Aeval, Beval, ctx_mp);
        if (!success)
        {
            continue;
        }
        fmpz_mod_mpolyn_scalar_mul_fmpz_mod(Geval, Gammaeval, ctx_mp);

        FLINT_ASSERT(Geval->length > 0);
        *GevaldegXY = fmpz_mod_mpolyn_bidegree(Geval);

        if (GdegboundXY < *GevaldegXY)
        {
            return random_check_image_degree_high;
        }
        else if (GdegboundXY > *GevaldegXY)
        {
            return random_check_image_degree_low;
        }

        /* reuse Bbareval for Heval */
        fmpz_mpolyuu_eval_fmpz_mod(Bbareval, ctx_mp, H, checkalpha, ctx);
        if (!fmpz_mod_mpolyn_equal(Bbareval, Geval, ctx_mp))
        {
            return random_check_image_no_match;
        }

        return random_check_good;
    }

    return random_check_point_not_found;
}


void _mpoly_monomial_evals_nmod(
    mp_limb_t * E,
    const ulong * Aexps,
    flint_bitcnt_t Abits,
    slong Alen,
    const mp_limb_t * alpha,
    slong vstart,
    const mpoly_ctx_t mctx,
    nmod_t fctx)
{
    slong i, j;
    slong offset, shift;
    slong N = mpoly_words_per_exp_sp(Abits, mctx);
    slong nvars = mctx->nvars;
    slong * LUToffset;
    ulong * LUTmask;
    mp_limb_t * LUTvalue;
    slong LUTlen;
    mp_limb_t xpoweval;
    ulong * inputexpmask;
    TMP_INIT;

    TMP_START;

    LUToffset = (slong *) TMP_ALLOC(N*FLINT_BITS*sizeof(slong));
    LUTmask   = (ulong *) TMP_ALLOC(N*FLINT_BITS*sizeof(ulong));
    LUTvalue  = (mp_limb_t *) TMP_ALLOC(N*FLINT_BITS*sizeof(mp_limb_t));

    inputexpmask = (ulong *) TMP_ALLOC(N*sizeof(ulong));
    mpoly_monomial_zero(inputexpmask, N);
    for (i = 0; i < Alen; i++)
    {
        for (j = 0; j < N; j++)
        {
            inputexpmask[j] |= (Aexps + N*i)[j];
        }
    }

    LUTlen = 0;
    for (j = nvars - 1; j >= vstart; j--)
    {
        mpoly_gen_offset_shift_sp(&offset, &shift, j, Abits, mctx);

        xpoweval = alpha[j]; /* xpoweval = alpha[i]^(2^i) */
        for (i = 0; i < Abits; i++)
        {
            LUToffset[LUTlen] = offset;
            LUTmask[LUTlen] = (UWORD(1) << (shift + i));
            LUTvalue[LUTlen] = xpoweval;
            if ((inputexpmask[offset] & LUTmask[LUTlen]) != 0)
                LUTlen++;
            xpoweval = nmod_mul(xpoweval, xpoweval, fctx);
        }
    }
    FLINT_ASSERT(LUTlen < N*FLINT_BITS);

    for (i = 0; i < Alen; i++)
    {
        xpoweval = 1;
        for (j = 0; j < LUTlen; j++)
        {
            if (((Aexps + N*i)[LUToffset[j]] & LUTmask[j]) != 0)
            {
                xpoweval = nmod_mul(xpoweval, LUTvalue[j], fctx);
            }
        }
        E[i] = xpoweval;
    }

    TMP_END;
}


static void nmod_mpoly_get_eval_helper(
    n_poly_t Acur,
    n_poly_t Ainc,
    const nmod_mpoly_ctx_t ctx_sp,
    const fmpz_mpoly_t A,
    const mp_limb_t * alphas,
    const fmpz_mpoly_ctx_t ctx)
{
    mp_limb_t * p, * q;
    const fmpz * Acoeffs = A->coeffs;
    slong j, Alen = A->length;

    n_poly_fit_length(Acur, Alen);
    Acur->length = Alen;

    n_poly_fit_length(Ainc, 2*Alen);
    Ainc->length = 2*Alen;

    p = Acur->coeffs;
    _mpoly_monomial_evals_nmod(p, A->exps, A->bits, Alen, alphas, 0,
                                              ctx->minfo, ctx_sp->ffinfo->mod);
    q = Ainc->coeffs;
    for (j = 0; j < Alen; j++)
    {
        q[2*j + 0] = fmpz_fdiv_ui(Acoeffs + j, ctx_sp->ffinfo->mod.n);
        q[2*j + 1] = p[j];
    }
}

static void nmod_mpolyu_get_eval_helper(
    n_bpoly_t Acur,
    n_polyun_t Ainc,
    const nmod_mpoly_ctx_t ctx_sp,
    const fmpz_mpolyu_t A,
    const mp_limb_t * alphas,
    const fmpz_mpoly_ctx_t ctx)
{
    slong i;

    n_bpoly_fit_length(Acur, A->length);
    Acur->length = A->length;

    n_polyun_fit_length(Ainc, A->length);
    Ainc->length = A->length;

    for (i = 0; i < A->length; i++)
    {
        Ainc->terms[i].exp = A->exps[i];
        nmod_mpoly_get_eval_helper(Acur->coeffs + i, Ainc->terms[i].coeff,
                                           ctx_sp, A->coeffs + i, alphas, ctx);
    }
}


void fmpz_mpolyuu_degrees_si(
    slong * degs,
    const fmpz_mpolyu_t A,
    const fmpz_mpoly_ctx_t ctx)
{
    slong i, j, * t;
    TMP_INIT;

    TMP_START;

    t = TMP_ALLOC(ctx->minfo->nvars*sizeof(slong));
    
    fmpz_mpoly_degrees_si(degs, A->coeffs + 0, ctx);

    for (i = 0; i < A->length; i++)
    {
        fmpz_mpoly_degrees_si(t, A->coeffs + i, ctx);
        for (j = 0; j < ctx->minfo->nvars; j++)
            degs[j] = FLINT_MAX(degs[j], t[j]);
    }

    TMP_END;
}
    

int fmpz_mpolyuu_gcd_berlekamp_massey(
    fmpz_mpolyu_t G,
    fmpz_mpolyu_t Abar,
    fmpz_mpolyu_t Bbar,
    fmpz_mpolyu_t A,
    fmpz_mpolyu_t B,
    const fmpz_mpoly_t Gamma,
    const fmpz_mpoly_ctx_t ctx)
{
    int which_check, changed, success, point_try_count;
    flint_bitcnt_t bits = A->bits;
    mpoly_bma_interpolate_ctx_t Ictx;
    nmod_mpoly_ctx_t ctx_sp;
    fmpz_mod_mpoly_ctx_t ctx_mp;
    fmpz_mpolyu_t H;
    fmpz_mpoly_t Hcontent;
    /* multi precision workspace */
    fmpz_mod_bma_mpoly_t Lambda;
    fmpz_mod_mpolyn_t Aeval, Beval, Geval, Abareval, Bbareval;
    fmpz_mpolycu_t Ainc, Acur, Binc, Bcur, Ared, Bred;
    fmpz_mpolyc_t Gammainc, Gammacur, Gammared;
    fmpz_t p, pm1, sshift, last_unlucky_sshift_plus_1, image_count;
    fmpz * checkalpha;
    /* single precision workspace */
    n_poly_bpoly_stack_t new_Sp;
    nmod_bma_mpoly_t GLambda_sp, AbarLambda_sp, BbarLambda_sp;
    n_bpoly_t new_Aeval_sp, new_Beval_sp, new_Geval_sp, new_Abareval_sp, new_Bbareval_sp;

    n_poly_t Gammacur_sp;
    n_bpoly_t Acur_sp, Bcur_sp;
    n_poly_t Gammainc_sp;
    n_polyun_t Ainc_sp, Binc_sp;

    mp_limb_t p_sp, sshift_sp, last_unlucky_sshift_plus_1_sp, image_count_sp;
    fmpz_t Gammaeval;
    mp_limb_t Gammaeval_sp;
    mp_limb_t * checkalpha_sp;
    /* misc */
    slong i, j;
    ulong GdegboundXY, GevaldegXY;
    slong * Adegs, * Bdegs;
    flint_rand_t randstate;
    fmpz_t subprod, cAksub, cBksub;
    int unlucky_count;
    fmpz_t Hmodulus;
    nmod_zip_mpolyu_t newZ;

    slong zip_evals;
    ulong ABtotal_length;

    FLINT_ASSERT(bits == A->bits);
    FLINT_ASSERT(bits == B->bits);
    FLINT_ASSERT(bits == G->bits);
    FLINT_ASSERT(bits == Abar->bits);
    FLINT_ASSERT(bits == Bbar->bits);
    FLINT_ASSERT(bits == Gamma->bits);

    /* let's initialize everything at once to avoid complicated cleanup */

    ABtotal_length = 0;
    for (i = 0; i < A->length; i++)
        ABtotal_length += (A->coeffs + i)->length;
    for (i = 0; i < B->length; i++)
        ABtotal_length += (B->coeffs + i)->length;

    flint_randinit(randstate);
    fmpz_init(p);
    fmpz_init(pm1); /* p - 1 */
    fmpz_init(image_count);
    fmpz_init(subprod);
    fmpz_init(cAksub);
    fmpz_init(cBksub);
    fmpz_init(sshift);
    fmpz_init(last_unlucky_sshift_plus_1);
    fmpz_init(Gammaeval);
    fmpz_init(Hmodulus);

    mpoly_bma_interpolate_ctx_init(Ictx, ctx->minfo->nvars);

    /* multiprecision workspace */
    fmpz_set_ui(p, 2);    /* modulus no care */
    fmpz_mod_mpoly_ctx_init(ctx_mp, 2, ORD_LEX, p); /* modulus no care */
    fmpz_mod_bma_mpoly_init(Lambda);

    fmpz_mod_mpolyn_init(Aeval, FLINT_BITS/2, ctx_mp);
    fmpz_mod_mpolyn_init(Beval, FLINT_BITS/2, ctx_mp);
    fmpz_mod_mpolyn_init(Geval, FLINT_BITS/2, ctx_mp);
    fmpz_mod_mpolyn_init(Abareval, FLINT_BITS/2, ctx_mp);
    fmpz_mod_mpolyn_init(Bbareval, FLINT_BITS/2, ctx_mp);
    fmpz_mpolyc_init(Gammainc);
    fmpz_mpolyc_init(Gammacur);
    fmpz_mpolyc_init(Gammared);
    fmpz_mpolycu_init(Ainc);
    fmpz_mpolycu_init(Acur);
    fmpz_mpolycu_init(Ared);
    fmpz_mpolycu_init(Binc);
    fmpz_mpolycu_init(Bcur);
    fmpz_mpolycu_init(Bred);
    checkalpha = (fmpz *) flint_malloc(ctx->minfo->nvars*sizeof(fmpz));
    for (i = 0; i < ctx->minfo->nvars; i++)
        fmpz_init(checkalpha + i);

    /* machine precision workspace "sp" */
    nmod_mpoly_ctx_init(ctx_sp, 2, ORD_LEX, 2); /* modulus no care */

    n_poly_stack_init(new_Sp->poly_stack);
    n_bpoly_stack_init(new_Sp->bpoly_stack);

    nmod_bma_mpoly_init(GLambda_sp);
    nmod_bma_mpoly_init(AbarLambda_sp);
    nmod_bma_mpoly_init(BbarLambda_sp);

    n_bpoly_init(new_Aeval_sp);
    n_bpoly_init(new_Beval_sp);
    n_bpoly_init(new_Geval_sp);
    n_bpoly_init(new_Abareval_sp);
    n_bpoly_init(new_Bbareval_sp);

    n_poly_init(Gammacur_sp);
    n_bpoly_init(Acur_sp);
    n_bpoly_init(Bcur_sp);
    n_poly_init(Gammainc_sp);
    n_polyun_init(Ainc_sp);
    n_polyun_init(Binc_sp);

    checkalpha_sp = (mp_limb_t *) flint_malloc(ctx->minfo->nvars*sizeof(mp_limb_t));

    /* the zippler */

    nmod_zip_mpolyu_init(newZ);

    Adegs      = (slong *) flint_malloc(ctx->minfo->nvars*sizeof(slong));
    Bdegs      = (slong *) flint_malloc(ctx->minfo->nvars*sizeof(slong));

    /* find a degree bound on G in the two main variables */
    GdegboundXY = FLINT_MIN(A->exps[0], B->exps[0]);
    p_sp = UWORD(1) << (FLINT_BITS - 2);
    for (point_try_count = 0; point_try_count < 10; point_try_count++)
    {
        p_sp = n_nextprime(p_sp, 1);
        nmod_mpoly_ctx_set_modulus(ctx_sp, p_sp);
        /* unfortunate nmod_poly's need mod set */

        for (i = 0; i < ctx->minfo->nvars; i++)
            checkalpha_sp[i] = n_urandint(randstate, p_sp);
        new_fmpz_mpolyuu_eval_nmod(new_Aeval_sp, ctx_sp, A, checkalpha_sp, ctx);
        new_fmpz_mpolyuu_eval_nmod(new_Beval_sp, ctx_sp, B, checkalpha_sp, ctx);

        if (new_Aeval_sp->length < 1 || new_Beval_sp->length < 1 ||
            n_bpoly_bidegree(new_Aeval_sp) != A->exps[0] ||
            n_bpoly_bidegree(new_Beval_sp) != B->exps[0])
        {
            /* evaluation killed at least one of lc(A) or lc(B) */
            continue;
        }
        success = n_bpoly_mod_gcd_brown_smprime(new_Geval_sp,
                  new_Abareval_sp, new_Bbareval_sp, new_Aeval_sp, new_Beval_sp, ctx_sp->ffinfo->mod, new_Sp);
        if (success)
        {
            FLINT_ASSERT(new_Geval_sp->length > 0);
            GdegboundXY = n_bpoly_bidegree(new_Geval_sp);
            break;
        }
    }

    fmpz_mpolyuu_degrees_si(Adegs, A, ctx);
    fmpz_mpolyuu_degrees_si(Bdegs, A, ctx);

    /*
        Find bits into which H can be packed. The degrees satsify
            deg_(x_i)(H) <= deg_(x_i)(A)
            deg_(x_i)(H) <= deg_(x_i)(B)
            deg_(x_i)(H) <= deg_(x_i)(Gamma) + deg_(x_i)(G)
    */
    for (i = 0; i < ctx->minfo->nvars; i++)
        Ictx->degbounds[i] = 1 + FLINT_MAX(Adegs[i], Bdegs[i]);

    fmpz_mpolyu_init(Abar, bits, ctx);
    fmpz_mpolyu_init(Bbar, bits, ctx);
    fmpz_mpolyu_init(H, bits, ctx);
    fmpz_mpoly_init3(Hcontent, 0, bits, ctx);

    /* initialization done! */

    if (GdegboundXY == 0)
    {
        fmpz_mpolyu_one(G, ctx);
        fmpz_mpolyu_swap(Abar, A, ctx);
        fmpz_mpolyu_swap(Bbar, B, ctx);
        success = 1;
        goto cleanup;
    }

    for (i = 0; i < ctx->minfo->nvars; i++)
        Ictx->subdegs[i] = Ictx->degbounds[i];

    goto got_ksub;

pick_ksub:

    if (ctx->minfo->nvars > 1)
    {
        /* just increment the smallest subdegs[j] */
        j = 1;
        for (i = 2; i < ctx->minfo->nvars; i++)
        {
            if (Ictx->subdegs[i] < Ictx->subdegs[j])
            {
                j = i;
            }
        }
        Ictx->subdegs[j]++;
    }

got_ksub:

    fmpz_one(subprod);
    for (i = 0; i < ctx->minfo->nvars; i++)
    {
        if ((slong)(Ictx->subdegs[i]) < 0)
        {
            /* ksub has overflown - absolute falure */
            success = 0;
            goto cleanup;
        }
        fmpz_mul_ui(subprod, subprod, Ictx->subdegs[i]);
    }

    /* see if the ksub killed either lc(A) or lc(B) */
    fmpz_mpoly_ksub_content(cAksub, A->coeffs + 0, Ictx->subdegs, ctx);
    fmpz_mpoly_ksub_content(cBksub, B->coeffs + 0, Ictx->subdegs, ctx);
    if (fmpz_is_zero(cAksub) || fmpz_is_zero(cBksub))
    {
        /* try a new substitution if we killed either leading coefficient */
        goto pick_ksub;
    }

pick_bma_prime:

    if (fmpz_cmp_ui(p, ABtotal_length) < 0)
        fmpz_set_ui(p, ABtotal_length);

    if (fmpz_cmp(p, subprod) < 0)
        fmpz_set(p, subprod);

    success = fmpz_next_smooth_prime(p, p);
    fmpz_sub_ui(pm1, p, 1);
    if (!success)
    {
        /* ran out of smooth primes - absolute falure */
        success = 0;
        goto cleanup;
    }

    /* make sure reduction does not kill either leading coeff after ksub */
    if (fmpz_divisible(cAksub, p) || fmpz_divisible(cBksub, p))
    {
        goto pick_bma_prime;
    }

    /* make sure p does not divide any coefficient of Gamma */
    for (i = 0; i < Gamma->length; i++)
    {
        if (fmpz_divisible(Gamma->coeffs + i, p))
        {
            goto pick_bma_prime;
        }
    }

    if (fmpz_abs_fits_ui(p))
    {
        p_sp = fmpz_get_ui(p);
        sshift_sp = 1;

        unlucky_count = 0;
        last_unlucky_sshift_plus_1_sp = 0;

        nmod_mpoly_ctx_set_modulus(ctx_sp, p_sp);
        nmod_discrete_log_pohlig_hellman_precompute_prime(Ictx->dlogenv_sp, p_sp);

        nmod_bma_mpoly_reset_prime(GLambda_sp, ctx_sp->ffinfo);
        nmod_bma_mpoly_reset_prime(AbarLambda_sp, ctx_sp->ffinfo);
        nmod_bma_mpoly_reset_prime(BbarLambda_sp, ctx_sp->ffinfo);
        nmod_bma_mpoly_zero(GLambda_sp);
        nmod_bma_mpoly_zero(AbarLambda_sp);
        nmod_bma_mpoly_zero(BbarLambda_sp);

        FLINT_ASSERT(sshift_sp == 1);
        nmod_mpoly_bma_interpolate_alpha_powers(checkalpha_sp, sshift_sp,
                                                    Ictx, ctx, ctx_sp->ffinfo);

        nmod_mpoly_get_eval_helper(Gammacur_sp, Gammainc_sp, ctx_sp, Gamma, checkalpha_sp, ctx);
        nmod_mpolyu_get_eval_helper(Acur_sp, Ainc_sp, ctx_sp, A, checkalpha_sp, ctx);
        nmod_mpolyu_get_eval_helper(Bcur_sp, Binc_sp, ctx_sp, B, checkalpha_sp, ctx);

        image_count_sp = 0;

    next_bma_image_sp:

        /* image count is also the current power of alpha we are evaluating */
        image_count_sp++;

        FLINT_ASSERT(GLambda_sp->pointcount == AbarLambda_sp->pointcount);
        FLINT_ASSERT(GLambda_sp->pointcount == BbarLambda_sp->pointcount);
        FLINT_ASSERT(sshift_sp + GLambda_sp->pointcount == image_count_sp);

        if (image_count_sp >= p_sp - 1)
        {
            /* out of evaluation points alpha^image_count in Fp* */
            goto pick_bma_prime;
        }

        Gammaeval_sp = n_poly_mod_eval_step2(Gammacur_sp, Gammainc_sp, ctx_sp->ffinfo->mod);
        new_nmod_mpolyuu_eval_step2(new_Aeval_sp, Acur_sp, Ainc_sp, ctx_sp);
        new_nmod_mpolyuu_eval_step2(new_Beval_sp, Bcur_sp, Binc_sp, ctx_sp);

        if (new_Aeval_sp->length < 1 || new_Beval_sp->length < 1 ||
            n_bpoly_bidegree(new_Aeval_sp) != A->exps[0] ||
            n_bpoly_bidegree(new_Beval_sp) != B->exps[0])
        {
            /* evaluation killed either lc(A) or lc(B) */
            sshift_sp += GLambda_sp->pointcount + 1;
            nmod_bma_mpoly_zero(GLambda_sp);
            nmod_bma_mpoly_zero(AbarLambda_sp);
            nmod_bma_mpoly_zero(BbarLambda_sp);
            goto next_bma_image_sp;
        }

        /* the evaluation killed neither lc(A) nor lc(B) */
        FLINT_ASSERT(Gammaeval_sp != 0);

        success = n_bpoly_mod_gcd_brown_smprime(new_Geval_sp,
              new_Abareval_sp, new_Bbareval_sp, new_Aeval_sp, new_Beval_sp,
                                               ctx_sp->ffinfo->mod, new_Sp);

        if (!success)
        {
            sshift_sp += GLambda_sp->pointcount + 1;
            nmod_bma_mpoly_zero(GLambda_sp);
            nmod_bma_mpoly_zero(AbarLambda_sp);
            nmod_bma_mpoly_zero(BbarLambda_sp);
            goto next_bma_image_sp;
        }

        FLINT_ASSERT(new_Geval_sp->length > 0);
        GevaldegXY = n_bpoly_bidegree(new_Geval_sp);
        n_bpoly_scalar_mul_nmod(new_Geval_sp, Gammaeval_sp, ctx_sp->ffinfo->mod);

        if (GdegboundXY < GevaldegXY)
        {
            /* this image in Fp[X,Y] was unlucky */
            if (sshift_sp == last_unlucky_sshift_plus_1_sp)
            {
                /* this ksub is probably unlucky */
                goto pick_ksub;
            }
            if (++unlucky_count > 2)
            {
                goto pick_bma_prime;
            }
            last_unlucky_sshift_plus_1_sp = sshift_sp + 1;
            sshift_sp += GLambda_sp->pointcount + 1;
            nmod_bma_mpoly_zero(GLambda_sp);
            nmod_bma_mpoly_zero(AbarLambda_sp);
            nmod_bma_mpoly_zero(BbarLambda_sp);
            goto next_bma_image_sp;        
        }
        else if (GdegboundXY > GevaldegXY)
        {
            /* new bound on deg_XY(G) */
            GdegboundXY = GevaldegXY;
            if (GdegboundXY == 0)
            {
                fmpz_mpolyu_one(G, ctx);
                fmpz_mpolyu_swap(Abar, A, ctx);
                fmpz_mpolyu_swap(Bbar, B, ctx);
                success = 1;
                goto cleanup;
            }
            sshift_sp += GLambda_sp->pointcount;
            nmod_bma_mpoly_zero(GLambda_sp);
            nmod_bma_mpoly_zero(AbarLambda_sp);
            nmod_bma_mpoly_zero(BbarLambda_sp);
            new_nmod_bma_mpoly_add_point(GLambda_sp, new_Geval_sp, ctx_sp);
            new_nmod_bma_mpoly_add_point(AbarLambda_sp, new_Abareval_sp, ctx_sp);
            new_nmod_bma_mpoly_add_point(BbarLambda_sp, new_Bbareval_sp, ctx_sp);
            goto next_bma_image_sp;
        }

        new_nmod_bma_mpoly_add_point(GLambda_sp, new_Geval_sp, ctx_sp);
        new_nmod_bma_mpoly_add_point(AbarLambda_sp, new_Abareval_sp, ctx_sp);
        new_nmod_bma_mpoly_add_point(BbarLambda_sp, new_Bbareval_sp, ctx_sp);

        if ((GLambda_sp->pointcount & 7) != 0)
        {
            goto next_bma_image_sp;
        }

        if (Gamma->length <= GLambda_sp->pointcount/2 &&
            !nmod_bma_mpoly_reduce(GLambda_sp))
        {
            flint_printf("G stabilized at %wd\n", GLambda_sp->pointcount);

            success = nmod_bma_mpoly_get_fmpz_mpolyu(H, ctx, sshift_sp,
                                              GLambda_sp, Ictx, ctx_sp->ffinfo);
            if (!success || H->length < 1 ||
                H->coeffs[0].length!= Gamma->length)
            {
                goto next_bma_image_sp;
            }

            which_check = 0;
            goto check_sp;
        }

        if (A->coeffs[0].length <= AbarLambda_sp->pointcount/2 &&
            !nmod_bma_mpoly_reduce(AbarLambda_sp))
        {
            flint_printf("Abar stabilized at %wd\n", AbarLambda_sp->pointcount);

            success = nmod_bma_mpoly_get_fmpz_mpolyu(H, ctx, sshift_sp,
                                          AbarLambda_sp, Ictx, ctx_sp->ffinfo);
            if (!success || H->length < 1 ||
                H->coeffs[0].length != A->coeffs[0].length)
            {
                goto next_bma_image_sp;
            }

            which_check = 1;
            goto check_sp;
        }

        if (B->coeffs[0].length <= BbarLambda_sp->pointcount/2 &&
            !nmod_bma_mpoly_reduce(BbarLambda_sp))
        {
            flint_printf("Bbar stabilized at %wd\n", BbarLambda_sp->pointcount);

            success = nmod_bma_mpoly_get_fmpz_mpolyu(H, ctx, sshift_sp,
                                          BbarLambda_sp, Ictx, ctx_sp->ffinfo);
            if (!success || H->length < 1 ||
                H->coeffs[0].length != B->coeffs[0].length)
            {
                goto next_bma_image_sp;
            }

            which_check = 2;
            goto check_sp;
        }

        goto next_bma_image_sp;

check_sp:

        if (which_check == 0)
        {
            /* GdegboundXY should be the bidegree of H */
            FLINT_ASSERT(GdegboundXY == H->exps[0]);
        }

        switch (_random_check_sp(&GevaldegXY, GdegboundXY, which_check,
                    new_Aeval_sp, new_Beval_sp, new_Geval_sp, new_Abareval_sp, new_Bbareval_sp,
                 checkalpha_sp, H, A, B, Gamma, ctx, ctx_sp, randstate, new_Sp))
        {
            default:
                FLINT_ASSERT(0);
            case random_check_image_no_match:
            case random_check_image_degree_high:
                goto next_bma_image_sp;
            case random_check_image_degree_low:
                /* the random evaluation point gave us a better degree bound */
                sshift_sp += GLambda_sp->pointcount;
                GdegboundXY = GevaldegXY;
                if (GdegboundXY == 0)
                {
                    fmpz_mpolyu_one(G, ctx);
                    fmpz_mpolyu_swap(Abar, A, ctx);
                    fmpz_mpolyu_swap(Bbar, B, ctx);
                    success = 1;
                    goto cleanup;
                }
                nmod_bma_mpoly_zero(GLambda_sp);
                nmod_bma_mpoly_zero(AbarLambda_sp);
                nmod_bma_mpoly_zero(BbarLambda_sp);
                goto next_bma_image_sp;
            case random_check_point_not_found:
                /* hmmm */
            case random_check_good:
                break;
        }
    }
    else
    {
        fmpz_one(sshift);

flint_printf("not implemented\n");
flint_abort();
which_check = -1;


        unlucky_count = 0;
        fmpz_zero(last_unlucky_sshift_plus_1);

        fmpz_mod_ctx_set_modulus(ctx_mp->ffinfo, p);
        fmpz_mod_discrete_log_pohlig_hellman_precompute_prime(Ictx->dlogenv, p);
        fmpz_mod_bma_mpoly_reset_prime(Lambda, ctx_mp->ffinfo);
        fmpz_mod_bma_mpoly_zero(Lambda);

        /* unfortunate fmpz_mod_poly's store their own ctx :( */
        fmpz_mod_mpolyn_set_modulus(Aeval, ctx_mp->ffinfo);
        fmpz_mod_mpolyn_set_modulus(Beval, ctx_mp->ffinfo);
        fmpz_mod_mpolyn_set_modulus(Geval, ctx_mp->ffinfo);
        fmpz_mod_mpolyn_set_modulus(Abareval, ctx_mp->ffinfo);
        fmpz_mod_mpolyn_set_modulus(Bbareval, ctx_mp->ffinfo);

        FLINT_ASSERT(fmpz_is_one(sshift));
        fmpz_mod_mpoly_bma_interpolate_alpha_powers(checkalpha, sshift,
                                                    Ictx, ctx, ctx_mp->ffinfo);

        /* set evaluation of monomials */
        fmpz_mod_mpoly_set_skel(Gammainc, ctx_mp, Gamma, checkalpha, ctx);
        fmpz_mod_mpolyu_set_skel(Ainc, ctx_mp, A, checkalpha, ctx);
        fmpz_mod_mpolyu_set_skel(Binc, ctx_mp, B, checkalpha, ctx);

        /* set reduction of coeffs */
        fmpz_mod_mpoly_red_skel(Gammared, Gamma, ctx_mp->ffinfo);
        fmpz_mod_mpolyu_red_skel(Ared, A, ctx_mp->ffinfo);
        fmpz_mod_mpolyu_red_skel(Bred, B, ctx_mp->ffinfo);

        /* copy evaluation of monomials */
        fmpz_mod_mpoly_copy_skel(Gammacur, Gammainc);
        fmpz_mod_mpolyu_copy_skel(Acur, Ainc);
        fmpz_mod_mpolyu_copy_skel(Bcur, Binc);

        fmpz_zero(image_count);

    next_bma_image:

        /* image count is also the current power of alpha we are evaluating */
        fmpz_add_ui(image_count, image_count, 1);

    #if WANT_ASSERT
        /* image_count == sshift + Lambda->pointcount */
        {
            fmpz_t t;
            fmpz_init(t);
            fmpz_add_ui(t, sshift, Lambda->pointcount);
            FLINT_ASSERT(fmpz_equal(t, image_count));
            fmpz_clear(t);
        }
    #endif

        if (fmpz_cmp(image_count, pm1) >= 0)
        {
            /* out of evaluation points alpha^image_count in Fp* */
            goto pick_bma_prime;
        }

        fmpz_mod_mpoly_use_skel_mul(Gammaeval, Gammared, Gammacur, Gammainc,
                                                               ctx_mp->ffinfo);
        fmpz_mod_mpolyuu_use_skel_mul(Aeval, A, Ared, Acur, Ainc, ctx_mp);
        fmpz_mod_mpolyuu_use_skel_mul(Beval, B, Bred, Bcur, Binc, ctx_mp);
        if (Aeval->length == 0 || Beval->length == 0
            || fmpz_mod_mpolyn_bidegree(Aeval) != A->exps[0]
            || fmpz_mod_mpolyn_bidegree(Beval) != B->exps[0])
        {
            /* evaluation killed either lc(A) or lc(B) */
            fmpz_add_ui(sshift, sshift, Lambda->pointcount + 1);
            fmpz_mod_bma_mpoly_zero(Lambda);
            goto next_bma_image;
        }

        /* the evaluation killed neither lc(A) nor lc(B) */
        FLINT_ASSERT(!fmpz_is_zero(Gammaeval));

        success = fmpz_mod_mpolyn_gcd_brown_bivar(Geval, Abareval, Bbareval,
                                                         Aeval, Beval, ctx_mp);
        if (!success)
        {
            fmpz_add_ui(sshift, sshift, Lambda->pointcount + 1);
            fmpz_mod_bma_mpoly_zero(Lambda);
            goto next_bma_image;
        }

        FLINT_ASSERT(Geval->length > 0);
        GevaldegXY = fmpz_mod_mpolyn_bidegree(Geval);
        fmpz_mod_mpolyn_scalar_mul_fmpz_mod(Geval, Gammaeval, ctx_mp);

        FLINT_ASSERT(fmpz_equal(Gammaeval, fmpz_mod_mpolyn_leadcoeff_last_ref(Geval, ctx_mp)));

        if (GdegboundXY < GevaldegXY)
        {
            /* this image in Fp[X,Y] was unlucky */
            if (fmpz_equal(sshift, last_unlucky_sshift_plus_1))
            {
                /* this ksub is probably unlucky */
                goto pick_ksub;
            }
            if (++unlucky_count > 2)
            {
                goto pick_bma_prime;
            }
            fmpz_add_ui(last_unlucky_sshift_plus_1, sshift, 1);
            fmpz_add_ui(sshift, sshift, Lambda->pointcount + 1);
            fmpz_mod_bma_mpoly_zero(Lambda);
            goto next_bma_image;        
        }
        else if (GdegboundXY > GevaldegXY)
        {
            /* new bound on deg_XY(G) */
            fmpz_add_ui(sshift, sshift, Lambda->pointcount);
            fmpz_mod_bma_mpoly_zero(Lambda);
            fmpz_mod_bma_mpoly_add_point(Lambda, Geval, ctx_mp);
            GdegboundXY = GevaldegXY;
            if (GdegboundXY == 0)
            {
                fmpz_mpolyu_one(G, ctx);
                fmpz_mpolyu_swap(Abar, A, ctx);
                fmpz_mpolyu_swap(Bbar, B, ctx);
                success = 1;
                goto cleanup;
            }
            goto next_bma_image;
        }

        fmpz_mod_bma_mpoly_add_point(Lambda, Geval, ctx_mp);
        if (   (Lambda->pointcount & 1) != 0
            || Gamma->length > Lambda->pointcount/2)
        {
            goto next_bma_image;
        }

        changed = fmpz_mod_bma_mpoly_reduce(Lambda);
        if (changed)
        {
            goto next_bma_image;
        }

        success = fmpz_mod_bma_mpoly_get_fmpz_mpolyu(H, ctx, sshift,
                                                 Lambda, Ictx, ctx_mp->ffinfo);
        if (!success
            || H->length == 0
            || (H->coeffs + 0)->length != Gamma->length)
        {
            goto next_bma_image;
        }

        /* GdegboundXY should be the bidegree of H */
        FLINT_ASSERT(GdegboundXY == H->exps[0]);

        switch (_random_check(&GevaldegXY, GdegboundXY,
                    Aeval, Beval, Geval, Abareval, Bbareval, Gammaeval,
                           checkalpha, H, A, B, Gamma, ctx, ctx_mp, randstate))
        {
            default:
                FLINT_ASSERT(0);
            case random_check_image_no_match:
            case random_check_image_degree_high:
                goto next_bma_image;
            case random_check_image_degree_low:
                /* the random evaluation point gave us a better degree bound */
                fmpz_add_ui(sshift, sshift, Lambda->pointcount);
                fmpz_mod_bma_mpoly_zero(Lambda);
                GdegboundXY = GevaldegXY;
                if (GdegboundXY == 0)
                {
                    fmpz_mpolyu_one(G, ctx);
                    fmpz_mpolyu_swap(Abar, A, ctx);
                    fmpz_mpolyu_swap(Bbar, B, ctx);
                    success = 1;
                    goto cleanup;
                }
                goto next_bma_image;
            case random_check_point_not_found:
                /* hmmm */
            case random_check_good:
                break;
        }
    }

    /* assume that H is correct modulo Hmodulus = p */
    FLINT_ASSERT(0 <= which_check && which_check <= 2);
    fmpz_set(Hmodulus, p);

    /* find number of evals for zip interp */
    FLINT_ASSERT(H->length > 0);
    zip_evals = H->coeffs[0].length;
    for (i = 1; i < H->length; i++)
        zip_evals = FLINT_MAX(zip_evals, H->coeffs[i].length);
    zip_evals += 1; /* one extra check eval */
    nmod_zip_mpolyu_fit_poly(newZ, H, zip_evals);

    p_sp = UWORD(1) << (FLINT_BITS - 2);

pick_zip_prime:

    if (p_sp >= UWORD_MAX_PRIME)
    {
        success = 0;
        goto cleanup;
    }
    p_sp = n_nextprime(p_sp, 1);

    if (0 == fmpz_fdiv_ui(Hmodulus, p_sp))
        goto pick_zip_prime;

    nmod_mpoly_ctx_set_modulus(ctx_sp, p_sp);

    FLINT_ASSERT(p_sp > 3);
    for (i = 0; i < ctx->minfo->nvars; i++)
        checkalpha_sp[i] = n_urandint(randstate, p_sp - 3) + 2;

    nmod_zip_mpolyu_set_skel(newZ, ctx_sp, H, checkalpha_sp, ctx);
    nmod_mpoly_get_eval_helper(Gammacur_sp, Gammainc_sp, ctx_sp, Gamma, checkalpha_sp, ctx);
    nmod_mpolyu_get_eval_helper(Acur_sp, Ainc_sp, ctx_sp, A, checkalpha_sp, ctx);
    nmod_mpolyu_get_eval_helper(Bcur_sp, Binc_sp, ctx_sp, B, checkalpha_sp, ctx);

next_zip_image:

    Gammaeval_sp = n_poly_mod_eval_step2(Gammacur_sp, Gammainc_sp, ctx_sp->ffinfo->mod);
    new_nmod_mpolyuu_eval_step2(new_Aeval_sp, Acur_sp, Ainc_sp, ctx_sp);
    new_nmod_mpolyuu_eval_step2(new_Beval_sp, Bcur_sp, Binc_sp, ctx_sp);

    if (new_Aeval_sp->length < 1 || new_Beval_sp->length < 1 ||
        n_bpoly_bidegree(new_Aeval_sp) != A->exps[0] ||
        n_bpoly_bidegree(new_Beval_sp) != B->exps[0])
    {
        goto pick_zip_prime;
    }

    /* the evaluation killed neither lc(A) nor lc(B) */
    FLINT_ASSERT(Gammaeval_sp != 0);

    success = n_bpoly_mod_gcd_brown_smprime(new_Geval_sp,
              new_Abareval_sp, new_Bbareval_sp, new_Aeval_sp, new_Beval_sp,
                                               ctx_sp->ffinfo->mod, new_Sp);
    if (!success)
        goto pick_zip_prime;

    FLINT_ASSERT(new_Geval_sp->length > 0);
    GevaldegXY = n_bpoly_bidegree(new_Geval_sp);

    if (GevaldegXY > GdegboundXY)
    {
        goto pick_zip_prime;        
    }
    else if (GevaldegXY < GdegboundXY)
    {
        GdegboundXY = GevaldegXY;
        if (GdegboundXY == 0)
        {
            fmpz_mpolyu_one(G, ctx);
            fmpz_mpolyu_swap(Abar, A, ctx);
            fmpz_mpolyu_swap(Bbar, B, ctx);
            success = 1;
            goto cleanup;
        }
        goto pick_bma_prime;
    }

    if (which_check == 1)
    {
        success = new_nmod_zip_mpolyuu_add_point(newZ, new_Abareval_sp);
    }
    else if (which_check == 2)
    {
        success = new_nmod_zip_mpolyuu_add_point(newZ, new_Bbareval_sp);
    }
    else
    {
        FLINT_ASSERT(which_check == 0);
        n_bpoly_scalar_mul_nmod(new_Geval_sp, Gammaeval_sp, ctx_sp->ffinfo->mod);
        success = new_nmod_zip_mpolyuu_add_point(newZ, new_Geval_sp);
    }

    if (!success)
        goto pick_bma_prime;

    if (newZ->pointcount < zip_evals)
        goto next_zip_image;

    switch (nmod_mpolyu_zip_find_coeffs(newZ, ctx_sp))
    {
        default:
            FLINT_ASSERT(0);
        case nmod_zip_find_coeffs_no_match:
            goto pick_bma_prime;
        case nmod_zip_find_coeffs_non_invertible:
            goto pick_zip_prime;
        case nmod_zip_find_coeffs_good:
            break;
    }

    FLINT_ASSERT(bits == H->bits);
    changed = fmpz_mpolyu_addinterp_zip(H, Hmodulus, newZ, ctx_sp->ffinfo);
    fmpz_mul_ui(Hmodulus, Hmodulus, ctx_sp->ffinfo->mod.n);

    if (changed)
    {
        /* TODO if the coefficients of H are getting to large? */
        goto pick_zip_prime;
    }

    success = fmpz_mpolyu_content_mpoly_threaded_pool(Hcontent, H, ctx, NULL, 0);
    FLINT_ASSERT(Hcontent->bits == bits);
    if (!success)
    {
        success = 0;
        goto cleanup;
    }

    if (which_check == 1)
    {
        fmpz_mpolyu_divexact_mpoly(Abar, H, 1, Hcontent, ctx);
        if (!fmpz_mpolyuu_divides(G, A, Abar, 2, ctx) ||
            !fmpz_mpolyuu_divides(Bbar, B, G, 2, ctx))
        {
            goto pick_zip_prime;
        }
    }
    else if (which_check == 2)
    {
        fmpz_mpolyu_divexact_mpoly(Bbar, H, 1, Hcontent, ctx);
        if (!fmpz_mpolyuu_divides(G, B, Bbar, 2, ctx) ||
            !fmpz_mpolyuu_divides(Abar, A, G, 2, ctx))
        {
            goto pick_zip_prime;
        }
    }
    else
    {
        FLINT_ASSERT(which_check == 0);
        fmpz_mpolyu_divexact_mpoly(G, H, 1, Hcontent, ctx);
        if (!fmpz_mpolyuu_divides(Abar, A, G, 2, ctx) ||
            !fmpz_mpolyuu_divides(Bbar, B, G, 2, ctx))
        {
            goto pick_zip_prime;
        }
    }

    success = 1;

cleanup:

    fmpz_mpoly_clear(Hcontent, ctx);
    fmpz_mpolyu_clear(H, ctx);

    flint_free(Adegs);
    flint_free(Bdegs);

    /* the zippler */
    nmod_zip_mpolyu_clear(newZ);

    /* machine precision workspace */
    flint_free(checkalpha_sp);

    n_bpoly_clear(new_Aeval_sp);
    n_bpoly_clear(new_Beval_sp);
    n_bpoly_clear(new_Geval_sp);
    n_bpoly_clear(new_Abareval_sp);
    n_bpoly_clear(new_Bbareval_sp);

    n_poly_clear(Gammacur_sp);
    n_bpoly_clear(Acur_sp);
    n_bpoly_clear(Bcur_sp);

    n_poly_clear(Gammainc_sp);
    n_polyun_clear(Ainc_sp);
    n_polyun_clear(Binc_sp);

    nmod_bma_mpoly_clear(GLambda_sp);
    nmod_bma_mpoly_clear(AbarLambda_sp);
    nmod_bma_mpoly_clear(BbarLambda_sp);

    n_poly_stack_clear(new_Sp->poly_stack);
    n_bpoly_stack_clear(new_Sp->bpoly_stack);
    nmod_mpoly_ctx_clear(ctx_sp);

    /* multiprecision workspace */
    for (i = 0; i < ctx->minfo->nvars; i++)
        fmpz_clear(checkalpha + i);

    flint_free(checkalpha);
    fmpz_mod_mpolyn_clear(Aeval, ctx_mp);
    fmpz_mod_mpolyn_clear(Beval, ctx_mp);
    fmpz_mod_mpolyn_clear(Geval, ctx_mp);
    fmpz_mod_mpolyn_clear(Abareval, ctx_mp);
    fmpz_mod_mpolyn_clear(Bbareval, ctx_mp);
    fmpz_mpolyc_clear(Gammainc);
    fmpz_mpolyc_clear(Gammacur);
    fmpz_mpolyc_clear(Gammared);
    fmpz_mpolycu_clear(Ainc);
    fmpz_mpolycu_clear(Acur);
    fmpz_mpolycu_clear(Ared);
    fmpz_mpolycu_clear(Binc);
    fmpz_mpolycu_clear(Bcur);
    fmpz_mpolycu_clear(Bred);
    fmpz_mod_bma_mpoly_clear(Lambda);
    fmpz_mod_mpoly_ctx_clear(ctx_mp);

    mpoly_bma_interpolate_ctx_clear(Ictx);

    fmpz_clear(Hmodulus);
    fmpz_clear(Gammaeval);
    fmpz_clear(last_unlucky_sshift_plus_1);
    fmpz_clear(sshift);
    fmpz_clear(cBksub);
    fmpz_clear(cAksub);
    fmpz_clear(subprod);
    fmpz_clear(image_count);
    fmpz_clear(pm1);
    fmpz_clear(p);
    flint_randclear(randstate);

    if (success)
    {
        FLINT_ASSERT(G->bits == bits);
        FLINT_ASSERT(Abar->bits == bits);
        FLINT_ASSERT(Bbar->bits == bits);
    }
    else
    {
        fmpz_mpolyu_set_bits(G, bits, ctx);
        fmpz_mpolyu_set_bits(Abar, bits, ctx);
        fmpz_mpolyu_set_bits(Bbar, bits, ctx);
    }

    return success;
}


int fmpz_mpoly_gcd_berlekamp_massey(
    fmpz_mpoly_t G,
    const fmpz_mpoly_t A,
    const fmpz_mpoly_t B,
    const fmpz_mpoly_ctx_t ctx)
{
    slong i;
    flint_bitcnt_t wbits;
    int success = 0;
    fmpz_mpoly_ctx_t uctx;
    fmpz_mpolyu_t Auu, Buu, Guu, Abaruu, Bbaruu;
    fmpz_mpoly_t Ac, Bc, Gc, Gamma;
    slong * Adegs, * Bdegs, * perm;
    ulong * shift, * stride;
    ulong max_main_degree, max_minor_degree;

    if (fmpz_mpoly_is_zero(A, ctx))
    {
        if (fmpz_mpoly_is_zero(B, ctx))
        {
            fmpz_mpoly_zero(G, ctx);
            return 1;
        }
        if (fmpz_sgn(B->coeffs + 0) < 0)
        {
            fmpz_mpoly_neg(G, B, ctx);
            return 1;
        }
        else
        {
            fmpz_mpoly_set(G, B, ctx);
            return 1;
        }
    }

    if (fmpz_mpoly_is_zero(B, ctx))
    {
        if (fmpz_sgn(A->coeffs + 0) < 0)
        {
            fmpz_mpoly_neg(G, A, ctx);
            return 1;
        }
        else
        {
            fmpz_mpoly_set(G, A, ctx);
            return 1;
        }
    }

    if (A->bits > FLINT_BITS || B->bits > FLINT_BITS)
    {
        return 0;
    }

    if (ctx->minfo->nvars < 3)
    {
        return fmpz_mpoly_gcd_zippel(G, A, B, ctx);
    }

    FLINT_ASSERT(A->bits <= FLINT_BITS);
    FLINT_ASSERT(B->bits <= FLINT_BITS);
    FLINT_ASSERT(ctx->minfo->nvars >= 3);
    FLINT_ASSERT(!fmpz_mpoly_is_zero(A, ctx));
    FLINT_ASSERT(!fmpz_mpoly_is_zero(B, ctx));

    Adegs = (slong *) flint_malloc(ctx->minfo->nvars*sizeof(slong));
    Bdegs = (slong *) flint_malloc(ctx->minfo->nvars*sizeof(slong));
    perm = (slong *) flint_malloc((ctx->minfo->nvars)*sizeof(slong));
    shift = (ulong *) flint_malloc((ctx->minfo->nvars)*sizeof(ulong));
    stride = (ulong *) flint_malloc((ctx->minfo->nvars)*sizeof(ulong));

    mpoly_degrees_si(Adegs, A->exps, A->length, A->bits, ctx->minfo);
    mpoly_degrees_si(Bdegs, B->exps, B->length, B->bits, ctx->minfo);

    max_main_degree = 0;
    max_minor_degree = 0;
    for (i = 0; i < ctx->minfo->nvars; i++)
    {
        perm[i] = i;
        shift[i] = 0;
        stride[i] = 1;
        FLINT_ASSERT(Adegs[i] >= 0);
        FLINT_ASSERT(Bdegs[i] >= 0);
        if (i < 2)
        {
            max_main_degree = FLINT_MAX(max_main_degree, Adegs[i]);
            max_main_degree = FLINT_MAX(max_main_degree, Bdegs[i]);
        }
        else
        {
            max_minor_degree = FLINT_MAX(max_minor_degree, Adegs[i]);
            max_minor_degree = FLINT_MAX(max_minor_degree, Bdegs[i]);
        }
    }

    fmpz_mpoly_ctx_init(uctx, ctx->minfo->nvars - 2, ORD_LEX);

    /* wbits is bits for intermediates in ZZ[x_0,x_1][x_2,...,x_(n-1)] */
    wbits = 1 + FLINT_BIT_COUNT(max_minor_degree);
    wbits = FLINT_MAX(MPOLY_MIN_BITS, wbits);
    wbits = mpoly_fix_bits(wbits, uctx->minfo);
    FLINT_ASSERT(wbits <= FLINT_BITS);

    fmpz_mpolyu_init(Auu, wbits, uctx);
    fmpz_mpolyu_init(Buu, wbits, uctx);
    fmpz_mpolyu_init(Guu, wbits, uctx);
    fmpz_mpolyu_init(Abaruu, wbits, uctx);
    fmpz_mpolyu_init(Bbaruu, wbits, uctx);
    fmpz_mpoly_init3(Ac, 0, wbits, uctx);
    fmpz_mpoly_init3(Bc, 0, wbits, uctx);
    fmpz_mpoly_init3(Gc, 0, wbits, uctx);
    fmpz_mpoly_init3(Gamma, 0, wbits, uctx);

    /* two main variables must be packed into bits = FLINT_BITS/2 */
    if (FLINT_BIT_COUNT(max_main_degree) >= FLINT_BITS/2)
    {
        success = 0;
        goto cleanup;
    }

    fmpz_mpoly_to_mpolyuu_perm_deflate_threaded_pool(Auu, uctx, A, ctx,
                                           perm, shift, stride, NULL, NULL, 0);
    fmpz_mpoly_to_mpolyuu_perm_deflate_threaded_pool(Buu, uctx, B, ctx,
                                           perm, shift, stride, NULL, NULL, 0);

    success = fmpz_mpolyu_content_mpoly_threaded_pool(Ac, Auu, uctx, NULL, 0);
    success = success &&
              fmpz_mpolyu_content_mpoly_threaded_pool(Bc, Buu, uctx, NULL, 0);
    if (!success)
        goto cleanup;

    fmpz_mpolyu_divexact_mpoly_inplace(Auu, Ac, uctx);
    fmpz_mpolyu_divexact_mpoly_inplace(Buu, Bc, uctx);

    success = _fmpz_mpoly_gcd_threaded_pool(Gamma, wbits, Auu->coeffs + 0,
                                            Buu->coeffs + 0, uctx, NULL, 0);
    if (!success)
        goto cleanup;

    success = fmpz_mpolyuu_gcd_berlekamp_massey(Guu, Abaruu, Bbaruu,
                                                        Auu, Buu, Gamma, uctx);
    if (!success)
        goto cleanup;

    success = _fmpz_mpoly_gcd_threaded_pool(Gc, wbits, Ac, Bc, uctx, NULL, 0);
    if (!success)
        goto cleanup;

    fmpz_mpolyu_mul_mpoly_inplace(Guu, Gc, uctx);

    fmpz_mpoly_from_mpolyuu_perm_inflate(G, FLINT_MIN(A->bits, B->bits), ctx,
                                               Guu, uctx, perm, shift, stride);
    if (fmpz_sgn(G->coeffs + 0) < 0)
        fmpz_mpoly_neg(G, G, ctx);

    success = 1;

cleanup:

    flint_free(Adegs);
    flint_free(Bdegs);
    flint_free(perm);
    flint_free(shift);
    flint_free(stride);

    fmpz_mpolyu_clear(Auu, uctx);
    fmpz_mpolyu_clear(Buu, uctx);
    fmpz_mpolyu_clear(Guu, uctx);
    fmpz_mpolyu_clear(Abaruu, uctx);
    fmpz_mpolyu_clear(Bbaruu, uctx);
    fmpz_mpoly_clear(Ac, uctx);
    fmpz_mpoly_clear(Bc, uctx);
    fmpz_mpoly_clear(Gc, uctx);
    fmpz_mpoly_clear(Gamma, uctx);

    fmpz_mpoly_ctx_clear(uctx);

    return success;
}

