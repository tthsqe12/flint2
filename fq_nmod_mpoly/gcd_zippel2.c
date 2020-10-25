/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fq_nmod_mpoly.h"
#include "fq_nmod_mpoly_factor.h"
#include "profiler.h"

void usleep(ulong);

void n_fq_polyun_set(n_fq_polyun_t A, const n_fq_polyun_t B, const fq_nmod_ctx_t ctx);

ulong _mpoly_bidegree(
    const ulong * Aexps,
    flint_bitcnt_t Abits,
    const mpoly_ctx_t mctx);

void nmod_mpoly_monomial_evals2_new(
    n_polyun_t E,
    const ulong * Aexps,
    flint_bitcnt_t Abits,
    slong Alen,
    const mp_limb_t * betas,
    slong m,
    const mpoly_ctx_t ctx,
    nmod_t mod);

void fq_nmod_mpoly_monomial_evals2_new(
    n_fq_polyun_t E,
    const ulong * Aexps,
    flint_bitcnt_t Abits,
    slong Alen,
    const fq_nmod_struct * betas,
    slong m,
    const fq_nmod_mpoly_ctx_t ctx);

void nmod_mpoly_monomial_evals_new(
    n_poly_t E,
    const ulong * Aexps,
    flint_bitcnt_t Abits,
    slong Alen,
    const mp_limb_t * betas,
    slong start,
    slong stop,
    const mpoly_ctx_t ctx,
    nmod_t mod);

void fq_nmod_mpoly_monomial_evals_new(
    n_fq_poly_t E,
    const ulong * Aexps,
    flint_bitcnt_t Abits,
    slong Alen,
    const fq_nmod_struct * betas,
    slong start,
    slong stop,
    const fq_nmod_mpoly_ctx_t ctx);

void n_fq_bpoly_eval_step_new(
    n_fq_bpoly_t E,
    n_fq_polyun_t cur,
    const n_fq_polyun_t inc,
    const fq_nmod_mpoly_t A,
    const fq_nmod_ctx_t ctx);


void n_fq_poly_eval_step_new(
    mp_limb_t * res,
    n_fq_poly_t cur,
    const n_fq_poly_t inc,
    const fq_nmod_mpoly_t A,
    const fq_nmod_ctx_t ctx);

void _n_fq_poly_evalp_step_new(
    mp_limb_t * res,            /* length d */
    mp_limb_t * cur,            /* length length */
    const mp_limb_t * inc,      /* length length */
    const mp_limb_t * coeffs,   /* length d*length */
    slong length,
    slong d,
    nmod_t mod)
{
    slong i, j;
    mp_limb_t p0, p1;
    mp_limb_t * tmp;
    TMP_INIT;

    if (length < 1)
    {
        _n_fq_zero(res, d);
        return;
    }

    TMP_START;
    tmp = (mp_limb_t *) TMP_ALLOC(3*d*sizeof(mp_limb_t));

    i = 0;

    for (j = 0; j < d; j++)
    {
        umul_ppmm(tmp[3*j+1], tmp[3*j+0], cur[i], (coeffs + d*i)[j]);
        tmp[3*j+2] = 0;
    }
    cur[i] = nmod_mul(cur[i], inc[i], mod);

    for (i = 1; i < length; i++)
    {
        for (j = 0; j < d; j++)
        {
            umul_ppmm(p1, p0, cur[i], (coeffs + d*i)[j]);
            add_sssaaaaaa(tmp[3*j+2], tmp[3*j+1], tmp[3*j+0],
                          tmp[3*j+2], tmp[3*j+1], tmp[3*j+0], 0, p1, p0);
        }
        cur[i] = nmod_mul(cur[i], inc[i], mod);
    }

    for (j = 0; j < d; j++)
        NMOD_RED3(res[j], tmp[3*j+2], tmp[3*j+1], tmp[3*j+0], mod);

    TMP_END;
}


void n_fq_bpoly_evalp_step_new(
    n_fq_bpoly_t E,
    n_polyun_t cur,
    const n_polyun_t inc,
    const fq_nmod_mpoly_t A,
    const fq_nmod_ctx_t ctx)
{
    slong d = fq_nmod_ctx_degree(ctx);
    slong i, Ai;
    slong e0, e1;
    mp_limb_t * c = FLINT_ARRAY_ALLOC(d, mp_limb_t);

    n_bpoly_zero(E);

    Ai = 0;
    for (i = 0; i < cur->length; i++)
    {
        slong this_len = cur->terms[i].coeff->length;

        _n_fq_poly_evalp_step_new(c, cur->terms[i].coeff->coeffs,
                                  inc->terms[i].coeff->coeffs,
                                  A->coeffs + d*Ai, this_len, d, ctx->mod);

        Ai += this_len;

        e0 = extract_exp(cur->terms[i].exp, 1, 2);
        e1 = extract_exp(cur->terms[i].exp, 0, 2);
        if (_n_fq_is_zero(c, d))
            continue;

        n_fq_bpoly_set_coeff_n_fq(E, e0, e1, c, ctx);
    }

    FLINT_ASSERT(Ai == A->length);

    flint_free(c);
}

void n_fq_poly_evalp_step_new(
    mp_limb_t * res,
    n_poly_t cur,
    const n_poly_t inc,
    const fq_nmod_mpoly_t A,
    const fq_nmod_ctx_t ctx)
{
    FLINT_ASSERT(A->length == cur->length);
    FLINT_ASSERT(A->length == inc->length);
    _n_fq_poly_evalp_step_new(res, cur->coeffs, inc->coeffs, A->coeffs,
                                 A->length, fq_nmod_ctx_degree(ctx), ctx->mod);
}


int nmod_mpoly_gcd_get_use_new(
    slong rGdeg,
    slong Adeg,
    slong Bdeg,
    slong gammadeg,
    slong degxAB,
    slong degyAB,
    slong numABgamma,
    const n_polyun_t G,
    const n_polyun_t Abar,
    const n_polyun_t Bbar);

slong n_polyun_product_roots(
    n_polyun_t M,
    const n_polyun_t H,
    nmod_t ctx);

slong n_fq_polyun_product_roots(
    n_fq_polyun_t M,
    const n_fq_polyun_t H,
    const fq_nmod_ctx_t ctx);

static ulong _fq_nmod_mpoly_bidegree(
    const fq_nmod_mpoly_t A,
    const fq_nmod_mpoly_ctx_t ctx)
{
    FLINT_ASSERT(A->length > 0);
    return _mpoly_bidegree(A->exps, A->bits, ctx->minfo);
}

static void fq_nmod_mpoly_monomial_evals2(
    n_fq_polyun_t E,
    const fq_nmod_mpoly_t A,
    const fq_nmod_struct * betas,
    slong m,
    const fq_nmod_mpoly_ctx_t ctx)
{
    fq_nmod_mpoly_monomial_evals2_new(E, A->exps, A->bits, A->length, betas, m, ctx);
}

static void fq_nmod_mpoly_monomial_evalsp2(
    n_polyun_t E,
    const fq_nmod_mpoly_t A,
    const mp_limb_t * betas,
    slong m,
    const fq_nmod_mpoly_ctx_t ctx)
{
    nmod_mpoly_monomial_evals2_new(E, A->exps, A->bits, A->length, betas, m,
                                                  ctx->minfo, ctx->fqctx->mod);
}

static void fq_nmod_mpoly_monomial_evalsp(
    n_poly_t E,
    const fq_nmod_mpoly_t A,
    const mp_limb_t * betas,
    slong start,
    slong stop,
    const fq_nmod_mpoly_ctx_t ctx)
{
    nmod_mpoly_monomial_evals_new(E, A->exps, A->bits, A->length, betas,
                                     start, stop, ctx->minfo, ctx->fqctx->mod);
}

static void fq_nmod_mpoly_monomial_evals(
    n_fq_poly_t E,
    const fq_nmod_mpoly_t A,
    const fq_nmod_struct * betas,
    slong start,
    slong stop,
    const fq_nmod_mpoly_ctx_t ctx)
{
    fq_nmod_mpoly_monomial_evals_new(E, A->exps, A->bits, A->length,
                                                      betas, start, stop, ctx);
}



int qzip_solvel(
    fq_nmod_mpoly_t A,
    n_fq_polyun_t Z,
    n_fq_polyun_t H,
    n_fq_polyun_t M,
    const fq_nmod_mpoly_ctx_t ctx);


int fq_nmod_zip_findp_coeffs(
    mp_limb_t * coeffs,             /* length mlength */
    const mp_limb_t * monomials,    /* length mlength */
    slong mlength,
    const mp_limb_t * evals,        /* length d*elength */
    slong elength,
    const mp_limb_t * master,       /* length d*(mlength + 1) */
    mp_limb_t * scratch,            /* length d*mlength */
    const fq_nmod_ctx_t ctx)
{
    slong d = fq_nmod_ctx_degree(ctx);
    nmod_t mod = fq_nmod_ctx_mod(ctx);
    int success;
    slong i, j, k;
    mp_limb_t * tmp = FLINT_ARRAY_ALLOC(d*20, mp_limb_t);
    mp_limb_t * V = tmp + 6*d;
    mp_limb_t * V0 = V + d;
    mp_limb_t * T = V0 + d;
    mp_limb_t * S = T + d;
    mp_limb_t * r = S + d;
    mp_limb_t * p0 = r + d;
    mp_limb_t * V_p = p0 + d;
    mp_limb_t r_p, T_p, S_p;

    FLINT_ASSERT(elength >= mlength);

    for (i = 0; i < mlength; i++)
    {
        /* coeffs[i] is (coeffs(P).values)/P(roots[i]) =: V/S
            where P(x) = master(x)/(x-roots[i])     */
        _n_fq_zero(V0, d);
        _n_fq_zero(T, d);
        _n_fq_zero(S, d);

        _nmod_vec_zero(V_p, 3*d);
        T_p = S_p = 0;
        r_p = (monomials + d*i)[0];
        for (j = mlength; j > 0; j--)
        {
            T_p = nmod_add(nmod_mul(r_p, T_p, mod), (master + d*j)[0], mod);
            S_p = nmod_add(nmod_mul(r_p, S_p, mod), T_p, mod);
            for (k = 0; k < d; k++)
            {
                mp_limb_t hi, lo;
                umul_ppmm(hi, lo, T_p, (evals + d*(j - 1))[k]);
                add_sssaaaaaa(V_p[3*k+2], V_p[3*k+1], V_p[3*k+0],
                              V_p[3*k+2], V_p[3*k+1], V_p[3*k+0], 0, hi, lo);

            }
        }

        /* roots[i] should be a root of master */
        FLINT_ASSERT(nmod_add(nmod_mul(r_p, T_p, mod), master[0], mod) == 0);

        S_p = nmod_mul(S_p, r_p, mod); /* shift is one */
        if (S == 0)
            return -1;

        S_p = nmod_inv(S_p, mod);

        for (k = 0; k < d; k++)
        {
            mp_limb_t vk;
            NMOD_RED3(vk, V_p[3*k+2], V_p[3*k+1], V_p[3*k+0], mod);
            (coeffs + d*i)[k] = nmod_mul(vk, S_p, mod);
        }
    }

    /* check that the remaining points match */
    for (j = 0; j < mlength; j++)
        scratch[j] = nmod_pow_ui((monomials + d*j)[0], mlength, mod);

    for (i = mlength; i < elength; i++)
    {
        _nmod_vec_zero(V_p, 3*d);
        S_p = 0;
        for (j = 0; j < mlength; j++)
        {
            scratch[j] = nmod_mul(scratch[j], (monomials + d*j)[0], mod);
            for (k = 0; k < d; k++)
            {
                mp_limb_t hi, lo;
                umul_ppmm(hi, lo, scratch[j], (coeffs + d*j)[k]);
                add_sssaaaaaa(V_p[3*k+2], V_p[3*k+1], V_p[3*k+0],
                              V_p[3*k+2], V_p[3*k+1], V_p[3*k+0], 0, hi, lo);
            }
        }

        for (k = 0; k < d; k++)
        {
            mp_limb_t vk;
            NMOD_RED3(vk, V_p[3*k+2], V_p[3*k+1], V_p[3*k+0], mod);
            if (vk != (evals + d*i)[k])
            {
                success = 0;
                goto cleanup;
            }
        }
    }

    success = 1;

cleanup:

    flint_free(tmp);

    return success;
}

int fq_nmod_zip_findp_coeffs_new(
    mp_limb_t * coeffs,             /* length d*mlength */
    const mp_limb_t * monomials,    /* length mlength */
    slong mlength,
    const mp_limb_t * evals,        /* length d*elength */
    slong elength,
    const mp_limb_t * master,       /* length (mlength + 1) */
    mp_limb_t * scratch,            /* length mlength */
    const fq_nmod_ctx_t ctx)
{
    slong d = fq_nmod_ctx_degree(ctx);
    nmod_t mod = fq_nmod_ctx_mod(ctx);
    int success;
    slong i, j, k;
    mp_limb_t * tmp = FLINT_ARRAY_ALLOC(d*20, mp_limb_t);
    mp_limb_t * V = tmp + 6*d;
    mp_limb_t * V0 = V + d;
    mp_limb_t * T = V0 + d;
    mp_limb_t * S = T + d;
    mp_limb_t * r = S + d;
    mp_limb_t * p0 = r + d;
    mp_limb_t * V_p = p0 + d;
    mp_limb_t r_p, T_p, S_p;

    FLINT_ASSERT(elength >= mlength);

    for (i = 0; i < mlength; i++)
    {
        /* coeffs[i] is (coeffs(P).values)/P(roots[i]) =: V/S
            where P(x) = master(x)/(x-roots[i])     */
        _n_fq_zero(V0, d);
        _n_fq_zero(T, d);
        _n_fq_zero(S, d);

        _nmod_vec_zero(V_p, 3*d);
        T_p = S_p = 0;
        r_p = monomials[i];
        for (j = mlength; j > 0; j--)
        {
            T_p = nmod_add(nmod_mul(r_p, T_p, mod), master[j], mod);
            S_p = nmod_add(nmod_mul(r_p, S_p, mod), T_p, mod);
            for (k = 0; k < d; k++)
            {
                mp_limb_t hi, lo;
                umul_ppmm(hi, lo, T_p, (evals + d*(j - 1))[k]);
                add_sssaaaaaa(V_p[3*k+2], V_p[3*k+1], V_p[3*k+0],
                              V_p[3*k+2], V_p[3*k+1], V_p[3*k+0], 0, hi, lo);

            }
        }

        /* roots[i] should be a root of master */
        FLINT_ASSERT(nmod_add(nmod_mul(r_p, T_p, mod), master[0], mod) == 0);

        S_p = nmod_mul(S_p, r_p, mod); /* shift is one */
        if (S == 0)
            return -1;

        S_p = nmod_inv(S_p, mod);

        for (k = 0; k < d; k++)
        {
            mp_limb_t vk;
            NMOD_RED3(vk, V_p[3*k+2], V_p[3*k+1], V_p[3*k+0], mod);
            (coeffs + d*i)[k] = nmod_mul(vk, S_p, mod);
        }
    }

    /* check that the remaining points match */
    for (j = 0; j < mlength; j++)
        scratch[j] = nmod_pow_ui(monomials[j], mlength, mod);

    for (i = mlength; i < elength; i++)
    {
        _nmod_vec_zero(V_p, 3*d);
        S_p = 0;
        for (j = 0; j < mlength; j++)
        {
            scratch[j] = nmod_mul(scratch[j], monomials[j], mod);
            for (k = 0; k < d; k++)
            {
                mp_limb_t hi, lo;
                umul_ppmm(hi, lo, scratch[j], (coeffs + d*j)[k]);
                add_sssaaaaaa(V_p[3*k+2], V_p[3*k+1], V_p[3*k+0],
                              V_p[3*k+2], V_p[3*k+1], V_p[3*k+0], 0, hi, lo);
            }
        }

        for (k = 0; k < d; k++)
        {
            mp_limb_t vk;
            NMOD_RED3(vk, V_p[3*k+2], V_p[3*k+1], V_p[3*k+0], mod);
            if (vk != (evals + d*i)[k])
            {
                success = 0;
                goto cleanup;
            }
        }
    }

    success = 1;

cleanup:

    flint_free(tmp);

    return success;
}



static int qzip_solvelp(
    fq_nmod_mpoly_t A,
    n_fq_polyun_t Z,
    n_fq_polyun_t H,
    n_fq_polyun_t M,
    const fq_nmod_mpoly_ctx_t ctx)
{
    slong d = fq_nmod_ctx_degree(ctx->fqctx);
    int success;
    slong Ai, i, n;
    n_poly_t t;

    n_poly_init(t);

    FLINT_ASSERT(Z->length == H->length);
    FLINT_ASSERT(Z->length == M->length);

    /* A has the right length, but possibly from a smaller ctx */
    if (A->length*d > A->coeffs_alloc)
    {
        slong new_alloc = FLINT_MAX(A->coeffs_alloc + A->coeffs_alloc/2, A->length*d);
        A->coeffs = (mp_limb_t *) flint_realloc(A->coeffs, new_alloc*sizeof(mp_limb_t));
        A->coeffs_alloc = new_alloc;
    }

    Ai = 0;
    for (i = 0; i < H->length; i++)
    {
        n = H->terms[i].coeff->length;
        FLINT_ASSERT(M->terms[i].coeff->length == n + 1);
        FLINT_ASSERT(Z->terms[i].coeff->length >= n);
        FLINT_ASSERT(Ai + n <= A->length);

        n_poly_fit_length(t, n);

        success = fq_nmod_zip_findp_coeffs_new(A->coeffs + d*Ai,
                         H->terms[i].coeff->coeffs, n,
                         Z->terms[i].coeff->coeffs, Z->terms[i].coeff->length,
                         M->terms[i].coeff->coeffs, t->coeffs, ctx->fqctx);
        if (success < 1)
        {
            n_poly_clear(t);
            return success;
        }

        Ai += n;
        FLINT_ASSERT(Ai <= A->length);
    }

    FLINT_ASSERT(Ai == A->length);

    n_poly_clear(t);
    return 1;
}





void fq_nmod_mpolyn_interp_lift_sm_bpoly(
    fq_nmod_mpolyn_t F,
    n_fq_bpoly_t A,
    const fq_nmod_mpoly_ctx_t ctx);

int fq_nmod_mpolyn_interp_crt_sm_bpoly(
    slong * lastdeg,
    fq_nmod_mpolyn_t F,
    fq_nmod_mpolyn_t T,
    const n_fq_bpoly_t A,
    const n_fq_poly_t modulus,
    n_fq_poly_t alphapow,
    const fq_nmod_mpoly_ctx_t ctx);


void fq_nmod_next_not_zero(fq_nmod_t alpha, const fq_nmod_ctx_t fqctx);

int fq_nmod_zip_find_coeffs_new_fq_nmod(
    mp_limb_t * coeffs,             /* length mlength */
    const mp_limb_t * monomials,    /* length mlength */
    slong mlength,
    const mp_limb_t * evals,        /* length d*elength */
    slong elength,
    const mp_limb_t * master,       /* length d*(mlength + 1) */
    mp_limb_t * temp,               /* length d*mlength */
    const fq_nmod_ctx_t ctx);

void fq_nmod_mpolyuun_print_pretty(
    const fq_nmod_mpolyun_t poly,
    const char ** x,
    slong nmainvars,
    const fq_nmod_mpoly_ctx_t ctx)
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
        fq_nmod_mpolyn_print_pretty(poly->coeffs + i, x, ctx);
        flint_printf(")");
        for (j = nmainvars - 1; j >= 0; j--)
        {
            flint_printf("*X%wd^%wd", nmainvars - 1 - j,
                    mask & (poly->exps[i] >> (FLINT_BITS/nmainvars*j)));
        }
    }
}

void fq_nmod_mpolyuu_print_pretty(
    const fq_nmod_mpolyu_t poly,
    const char ** x,
    slong nmainvars,
    const fq_nmod_mpoly_ctx_t ctx)
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
        fq_nmod_mpoly_print_pretty(poly->coeffs + i, x, ctx);
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

void fq_nmod_mpoly_set_eval_helper(
    n_fq_poly_t EH,
    const fq_nmod_mpoly_t A,
    const fq_nmod_struct * betas,
    slong mvars,
    const fq_nmod_mpoly_ctx_t ctx)
{
    slong d = fq_nmod_ctx_degree(ctx->fqctx);
    slong j, n;
    mp_limb_t * p, * q;

    n = A->length;
    n_poly_fit_length(EH, 3*d*n);
    EH->length = n;
    p = EH->coeffs;
    _fq_nmod_mpoly_monomial_evals(p, A->exps, A->bits, n, betas, 0, mvars, ctx);
    q = A->coeffs;
    for (j = n - 1; j >= 0; j--)
    {
        _n_fq_set(p + d*(3*j + 2), p + d*j, d);
        _n_fq_set(p + d*(3*j + 0), p + d*(3*j + 2), d);
        _n_fq_set(p + d*(3*j + 1), q + d*j, d);
    }
}


void fq_nmod_mpoly_set_evalp_helper(
    n_fq_poly_t EH,
    const fq_nmod_mpoly_t A,
    const fq_nmod_struct * betas,
    slong mvars,
    const fq_nmod_mpoly_ctx_t ctx)
{
    slong d = fq_nmod_ctx_degree(ctx->fqctx);
    slong j, n;
    mp_limb_t * p, * q;

    n = A->length;
    n_poly_fit_length(EH, 3*d*n);
    EH->length = n;
    p = EH->coeffs;
    _fq_nmod_mpoly_monomial_evals(p, A->exps, A->bits, n, betas, 0, mvars, ctx);
    q = A->coeffs;
    for (j = n - 1; j >= 0; j--)
    {
        _n_fq_set(p + d*(3*j + 2), p + d*j, d);
        _n_fq_set(p + d*(3*j + 0), p + d*(3*j + 2), d);
        _n_fq_set(p + d*(3*j + 1), q + d*j, d);
    }
}


void fq_nmod_mpolyu_set_eval_helper(
    n_fq_polyun_t EH,
    const fq_nmod_mpolyu_t A,
    const fq_nmod_struct * betas,
    slong mvars,
    const fq_nmod_mpoly_ctx_t ctx)
{
    slong d = fq_nmod_ctx_degree(ctx->fqctx);
    slong i, j, n;
    fq_nmod_mpoly_struct * Acoeffs = A->coeffs;
    ulong * Aexps = A->exps;
    flint_bitcnt_t Abits = A->bits;
    n_fq_polyun_term_struct * EHterms;
    mp_limb_t * p, * q;

    n_polyun_fit_length(EH, A->length);
    EH->length = A->length;
    EHterms = EH->terms;

    for (i = 0; i < A->length; i++)
    {
        EHterms[i].exp = Aexps[i];
        n = Acoeffs[i].length;
        n_poly_fit_length(EHterms[i].coeff, d*3*n);
        EHterms[i].coeff->length = n;
        p = EHterms[i].coeff->coeffs;
        _fq_nmod_mpoly_monomial_evals(p, Acoeffs[i].exps, Abits, n, betas,
                                                                0, mvars, ctx);
        q = Acoeffs[i].coeffs;
        for (j = n - 1; j >= 0; j--)
        {
            _n_fq_set(p + d*(3*j + 2), p + d*j, d);
            _n_fq_set(p + d*(3*j + 0), p + d*(3*j + 2), d);
            _n_fq_set(p + d*(3*j + 1), q + d*j, d);
        }
    }
}


void fq_nmod_mpolyu_evaluate_one_fq_nmod(
    fq_nmod_mpolyu_t E,
    fq_nmod_mpolyu_t A,
    slong var,
    const fq_nmod_t alpha,
    const fq_nmod_mpoly_ctx_t ctx)
{
    slong i, Ei;

    FLINT_ASSERT(E != A);
    FLINT_ASSERT(E->bits == A->bits);
    FLINT_ASSERT(fq_nmod_mpolyu_is_canonical(A, ctx));

    fq_nmod_mpolyu_fit_length(E, A->length, ctx);
    Ei = 0;
    for (i = 0; i < A->length; i++)
    {
        E->exps[Ei] = A->exps[i];
        fq_nmod_mpoly_evaluate_one_fq_nmod(E->coeffs + Ei, A->coeffs + i, var, alpha, ctx);
        fq_nmod_mpoly_repack_bits_inplace(E->coeffs + Ei, E->bits, ctx);
        Ei += !fq_nmod_mpoly_is_zero(E->coeffs + Ei, ctx);
    }

    E->length = Ei;
}

void fq_nmod_mpolyuu_get_n_fq_bpoly(
    n_fq_bpoly_t A,
    fq_nmod_mpolyu_t B,
    const fq_nmod_mpoly_ctx_t ctx)
{
    slong i;
    n_fq_bpoly_zero(A);
    for (i = 0; i < B->length; i++)
    {
        slong x = extract_exp(B->exps[i], 1, 2);
        slong y = extract_exp(B->exps[i], 0, 2);
        n_fq_bpoly_set_coeff_n_fq(A, x, y,
                fq_nmod_mpoly_get_nonzero_n_fq(B->coeffs + i, ctx), ctx->fqctx);
    }
}


void fq_nmod_mpolyuun_interp_lift_sm_bpoly(
    fq_nmod_mpolyun_t C,
    n_fq_bpoly_t B,
    const fq_nmod_mpoly_ctx_t ctx)
{
    slong d = fq_nmod_ctx_degree(ctx->fqctx);
    slong N = mpoly_words_per_exp_sp(C->bits, ctx->minfo);
    slong i, j, Ci;

    Ci = 0;
    for (i = B->length - 1; i >= 0; i--)
    {
        n_poly_struct * Bi = B->coeffs + i;
        for (j = Bi->length - 1; j >= 0; j--)
        {
            if (_n_fq_is_zero(Bi->coeffs + d*j, d))
                continue;

            fq_nmod_mpolyun_fit_length(C, Ci + 1, ctx);
            C->exps[Ci] = pack_exp2(i, j);
            fq_nmod_mpolyn_fit_length(C->coeffs + Ci, 1, ctx);
            C->coeffs[Ci].length = 1;
            mpoly_monomial_zero(C->coeffs[Ci].exps + N*0, N);
            n_fq_poly_set_n_fq(C->coeffs[Ci].coeffs + 0, Bi->coeffs + d*j, ctx->fqctx);
            Ci++;
        }
    }

    C->length = Ci;
}

static void fq_nmod_mpolyn_interp_lift_lg_bpoly(
    slong * lastdeg_,
    fq_nmod_mpolyn_t F,
    const fq_nmod_mpoly_ctx_t smctx,
    n_fq_bpoly_t A,
    const fq_nmod_mpoly_ctx_t lgctx,
    const bad_fq_nmod_embed_t emb)
{
    slong lgd = fq_nmod_ctx_degree(lgctx->fqctx);
    slong N = mpoly_words_per_exp_sp(F->bits, smctx->minfo);
    slong i, j, Fi;
    slong lastdeg = -WORD(1);
    slong off0, shift0, off1, shift1;

    mpoly_gen_offset_shift_sp(&off0, &shift0, 0, F->bits, smctx->minfo);
    mpoly_gen_offset_shift_sp(&off1, &shift1, 1, F->bits, smctx->minfo);

    Fi = 0;
    for (i = A->length - 1; i >= 0; i--)
    {
        n_fq_poly_struct * Ai = A->coeffs + i;
        for (j = Ai->length - 1; j >= 0; j--)
        {
            if (_n_fq_is_zero(Ai->coeffs + lgd*j, lgd))
                continue;

            fq_nmod_mpolyn_fit_length(F, Fi + 1, smctx);
            mpoly_monomial_zero(F->exps + N*Fi, N);
            (F->exps + N*Fi)[off0] += (i << shift0);
            (F->exps + N*Fi)[off1] += (j << shift1);
            bad_n_fq_embed_lg_to_sm(F->coeffs + Fi, Ai->coeffs + lgd*j, emb);
            lastdeg = FLINT_MAX(lastdeg, n_fq_poly_degree(F->coeffs + Fi));

            Fi++;
        }
    }

    F->length = Fi;

    *lastdeg_ = lastdeg;
}


static n_fq_poly_struct * getcf(fq_nmod_mpolyn_t A)
{
    FLINT_ASSERT(A->length == 1);
    FLINT_ASSERT(A->exps[0] == 0);
    return A->coeffs + 0;
}

int fq_nmod_mpolyuun_interp_crt_sm_bpoly(
    slong * lastdeg,
    fq_nmod_mpolyun_t F,
    fq_nmod_mpolyun_t T,
    n_fq_bpoly_t A,
    n_fq_poly_t modulus,
    n_fq_poly_t alphapow,
    const fq_nmod_mpoly_ctx_t ctx)
{
    int changed = 0;
    slong d = fq_nmod_ctx_degree(ctx->fqctx);
    slong N = mpoly_words_per_exp(T->bits, ctx->minfo);
    n_fq_poly_struct * Acoeffs = A->coeffs;
    slong Fi, Ti, Ai, ai;
    slong Flen = F->length;
    ulong * Fexps = F->exps;
    fq_nmod_mpolyn_struct * Fcoeffs = F->coeffs;
    ulong * Texps = T->exps;
    fq_nmod_mpolyn_struct * Tcoeffs = T->coeffs;
    mp_limb_t * v = FLINT_ARRAY_ALLOC(d, mp_limb_t);

    FLINT_ASSERT(fq_nmod_mpolyun_is_canonical(F, ctx));
    FLINT_ASSERT(n_fq_bpoly_is_canonical(A, ctx->fqctx));

    FLINT_ASSERT(T->bits == F->bits);

    *lastdeg = -1;

    Ti = Fi = 0;
    Ai = A->length - 1;
    ai = (Ai < 0) ? 0 : n_fq_poly_degree(Acoeffs + Ai);

    while (Fi < Flen || Ai >= 0)
    {
        if (Ti >= T->alloc)
        {
            slong extra = Flen - Fi;
            extra = FLINT_MAX(extra, Ai);
            fq_nmod_mpolyun_fit_length(T, Ti + extra + 1, ctx);
            Tcoeffs = T->coeffs;
            Texps = T->exps;
        }

        fq_nmod_mpolyn_fit_length(Tcoeffs + Ti, 1, ctx);
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

            n_fq_poly_eval_pow(v, getcf(Fcoeffs + Fi), alphapow, ctx->fqctx);
            n_fq_sub(v, Acoeffs[Ai].coeffs + d*ai, v, ctx->fqctx);
            if (!_n_fq_is_zero(v, d))
            {
                changed = 1;
                n_fq_poly_scalar_addmul_n_fq(getcf(Tcoeffs + Ti),
                                  getcf(Fcoeffs + Fi), modulus, v, ctx->fqctx);
            }
            else
            {
                n_fq_poly_set(getcf(Tcoeffs + Ti), getcf(Fcoeffs + Fi), ctx->fqctx);
            }

            Fi++;
            do {
                ai--;
            } while (ai >= 0 && _n_fq_is_zero(Acoeffs[Ai].coeffs + d*ai, d));
            if (ai < 0)
            {
                do {
                    Ai--;
                } while (Ai >= 0 && Acoeffs[Ai].length == 0);
                if (Ai >= 0)
                    ai = n_fq_poly_degree(Acoeffs + Ai);
            }
        }
        else if (Ai >= 0 && (Fi >= Flen || Fexps[Fi] < pack_exp2(Ai, ai)))
        {
            /* F term missing, A term ok */

            Texps[Ti] = pack_exp2(Ai, ai);

            changed = 1;
            n_fq_poly_scalar_mul_n_fq(getcf(Tcoeffs + Ti), modulus,
                                        Acoeffs[Ai].coeffs + d*ai, ctx->fqctx);
            do {
                ai--;
            } while (ai >= 0 && _n_fq_is_zero(Acoeffs[Ai].coeffs + d*ai, d));
            if (ai < 0)
            {
                do {
                    Ai--;
                } while (Ai >= 0 && Acoeffs[Ai].length == 0);
                if (Ai >= 0)
                    ai = n_fq_poly_degree(Acoeffs + Ai);
            }
        }
        else
        {
            FLINT_ASSERT(Fi < Flen && (Ai < 0 || Fexps[Fi] > pack_exp2(Ai, ai)));

            /* F term ok, Aterm missing */

            Texps[Ti] = Fexps[Fi];

            n_fq_poly_eval_pow(v, getcf(Fcoeffs + Fi), alphapow, ctx->fqctx);
            if (!_n_fq_is_zero(v, d))
            {
                changed = 1;
flint_printf("mpolyuu_gcd_zippel.c line 451: this code looks wrong");
flint_abort();
                n_fq_poly_scalar_addmul_n_fq(getcf(Tcoeffs + Ti), 
                                  getcf(Fcoeffs + Fi), modulus, v, ctx->fqctx);
            }
            else
            {
                n_fq_poly_set(getcf(Tcoeffs + Ti), getcf(Fcoeffs + Fi), ctx->fqctx);                
            }

            Fi++;
        }

        FLINT_ASSERT(!n_fq_poly_is_zero(getcf(Tcoeffs + Ti)));
        *lastdeg = FLINT_MAX(*lastdeg, n_fq_poly_degree(getcf(Tcoeffs + Ti)));
        Ti++;
    }

    T->length = Ti;

    if (changed)
        fq_nmod_mpolyun_swap(T, F);

    FLINT_ASSERT(fq_nmod_mpolyun_is_canonical(F, ctx));

    flint_free(v);

    return changed;
}

int fq_nmod_mpolyn_interp_crt_lg_bpoly(
    slong * lastdeg,
    fq_nmod_mpolyn_t F,
    fq_nmod_mpolyn_t T,
    n_fq_poly_t modulus,
    const fq_nmod_mpoly_ctx_t smctx,
    n_fq_bpoly_t A,
    const fq_nmod_mpoly_ctx_t lgctx,
    const bad_fq_nmod_embed_t emb)
{
    slong lgd = fq_nmod_ctx_degree(lgctx->fqctx);
    int changed = 0;
    slong N = mpoly_words_per_exp(T->bits, smctx->minfo);
    slong off0, shift0, off1, shift1;
    n_fq_poly_struct * Acoeffs = A->coeffs;
    slong Fi, Ti, Ai, ai;
    slong Flen = F->length;
    ulong * Fexps = F->exps;
    n_fq_poly_struct * Fcoeffs = F->coeffs;
    ulong * Texps = T->exps;
    n_fq_poly_struct * Tcoeffs = T->coeffs;
    mp_limb_t * u = FLINT_ARRAY_ALLOC(3*lgd, mp_limb_t);
    mp_limb_t * v = u + lgd;
    mp_limb_t * inv_m_eval = v + lgd;
    n_fq_poly_t u_sm;
    ulong Fexpi, mask;

    mask = (-UWORD(1)) >> (FLINT_BITS - F->bits);
    mpoly_gen_offset_shift_sp(&off0, &shift0, 0, F->bits, smctx->minfo);
    mpoly_gen_offset_shift_sp(&off1, &shift1, 1, F->bits, smctx->minfo);

    n_fq_poly_init(u_sm);

    bad_n_fq_embed_sm_to_lg(u, modulus, emb);
    n_fq_inv(inv_m_eval, u, lgctx->fqctx);

    FLINT_ASSERT(fq_nmod_mpolyn_is_canonical(F, smctx));
    FLINT_ASSERT(n_fq_bpoly_is_canonical(A, lgctx->fqctx));

    FLINT_ASSERT(T->bits == F->bits);

    *lastdeg = -1;

    Ti = Fi = 0;
    Ai = A->length - 1;
    ai = (Ai < 0) ? 0 : n_fq_poly_degree(Acoeffs + Ai);

    while (Fi < Flen || Ai >= 0)
    {
        if (Ti >= T->alloc)
        {
            slong extra = FLINT_MAX(Flen - Fi, Ai);
            fq_nmod_mpolyn_fit_length(T, Ti + extra + 1, smctx);
            Tcoeffs = T->coeffs;
            Texps = T->exps;
        }

        if (Fi < Flen)
        {
            Fexpi = pack_exp2(((Fexps + N*Fi)[off0]>>shift0)&mask,
                              ((Fexps + N*Fi)[off1]>>shift1)&mask);
        }
        else
        {
            Fexpi = 0;
        }

        if (Fi < Flen && Ai >= 0 && Fexpi == pack_exp2(Ai, ai))
        {
            /* F term ok, A term ok */
            mpoly_monomial_set(Texps + N*Ti, Fexps + N*Fi, N);

            bad_n_fq_embed_sm_to_lg(u, Fcoeffs + Fi, emb);
            n_fq_sub(v, Acoeffs[Ai].coeffs + lgd*ai, u, lgctx->fqctx);
            if (!_n_fq_is_zero(v, lgd))
            {
                changed = 1;
                n_fq_mul(u, v, inv_m_eval, lgctx->fqctx);
                bad_n_fq_embed_lg_to_sm(u_sm, u, emb);
                n_fq_poly_mul(Tcoeffs + Ti, modulus, u_sm, smctx->fqctx);
                n_fq_poly_add(Tcoeffs + Ti, Tcoeffs + Ti, Fcoeffs + Fi, smctx->fqctx);
            }
            else
            {
                n_fq_poly_set(Tcoeffs + Ti, Fcoeffs + Fi, smctx->fqctx);
            }

            Fi++;
            do {
                ai--;
            } while (ai >= 0 && _n_fq_is_zero(Acoeffs[Ai].coeffs + lgd*ai, lgd));
            if (ai < 0)
            {
                do {
                    Ai--;
                } while (Ai >= 0 && Acoeffs[Ai].length == 0);
                if (Ai >= 0)
                    ai = n_fq_poly_degree(Acoeffs + Ai);
            }
        }
        else if (Fi < Flen && (Ai < 0 || Fexpi > pack_exp2(Ai, ai)))
        {
            /* F term ok, Aterm missing */
            mpoly_monomial_set(Texps + N*Ti, Fexps + N*Fi, N);

            bad_n_fq_embed_sm_to_lg(v, Fcoeffs + Fi, emb);
            if (!_n_fq_is_zero(v, lgd))
            {
                changed = 1;
                n_fq_mul(u, v, inv_m_eval, lgctx->fqctx);
                bad_n_fq_embed_lg_to_sm(u_sm, u, emb);
                n_fq_poly_mul(Tcoeffs + Ti, modulus, u_sm, smctx->fqctx);
                n_fq_poly_sub(Tcoeffs + Ti, Fcoeffs + Fi, Tcoeffs + Ti, smctx->fqctx);
            }
            else
            {
                n_fq_poly_set(Tcoeffs + Ti, Fcoeffs + Fi, smctx->fqctx);
            }

            Fi++;
        }
        else
        {
            FLINT_ASSERT(Ai >= 0 && (Fi >= Flen || Fexpi < pack_exp2(Ai, ai)));

            /* F term missing, A term ok */
            mpoly_monomial_zero(Texps + N*Ti, N);
            (Texps + N*Ti)[off0] += (Ai << shift0);
            (Texps + N*Ti)[off1] += (ai << shift1);

            changed = 1;
            n_fq_mul(u, Acoeffs[Ai].coeffs + lgd*ai, inv_m_eval, lgctx->fqctx);
            bad_n_fq_embed_lg_to_sm(u_sm, u, emb);
            n_fq_poly_mul(Tcoeffs + Ti, modulus, u_sm, smctx->fqctx);

            do {
                ai--;
            } while (ai >= 0 && _n_fq_is_zero(Acoeffs[Ai].coeffs + lgd*ai, lgd));
            if (ai < 0)
            {
                do {
                    Ai--;
                } while (Ai >= 0 && Acoeffs[Ai].length == 0);
                if (Ai >= 0)
                    ai = n_fq_poly_degree(Acoeffs + Ai);
            }
        }

        FLINT_ASSERT(!n_fq_poly_is_zero(Tcoeffs + Ti));
        *lastdeg = FLINT_MAX(*lastdeg, n_fq_poly_degree(Tcoeffs + Ti));
        Ti++;
    }

    T->length = Ti;

    if (changed)
        fq_nmod_mpolyn_swap(T, F);

    FLINT_ASSERT(fq_nmod_mpolyn_is_canonical(F, smctx));

    n_fq_poly_clear(u_sm);
    flint_free(u);

    return changed;
}

slong fq_nmod_mpolyu_set_zip_form(
    n_fq_polyun_t H, /* monomial evals */
    n_fq_polyun_t M, /* master polys */
    const fq_nmod_mpolyu_t A,
    const fq_nmod_struct * betas,
    slong mvars,
    const fq_nmod_mpoly_ctx_t ctx)
{
    slong d = fq_nmod_ctx_degree(ctx->fqctx);
    slong i, n, zip_length = 0;
    n_fq_polyun_term_struct * Hterms, * Mterms;

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
        n_poly_fit_length(Hterms[i].coeff, d*n);
        Hterms[i].coeff->length = n;
        _fq_nmod_mpoly_monomial_evals(Hterms[i].coeff->coeffs,
                          A->coeffs[i].exps, A->bits, n, betas, 0, mvars, ctx);
        n_fq_poly_product_roots_n_fq(Mterms[i].coeff,
                                       Hterms[i].coeff->coeffs, n, ctx->fqctx);
    }

    return zip_length;
}


void qzip_start(
    n_fq_polyun_t Z,
    n_fq_polyun_t H,
    slong req_images,
    const fq_nmod_ctx_t ctx)
{
    slong d = fq_nmod_ctx_degree(ctx);
    slong j;
    n_polyun_fit_length(Z, H->length);
    Z->length = H->length;
    for (j = 0; j < H->length; j++)
    {
        Z->terms[j].exp = H->terms[j].exp;
        n_poly_fit_length(Z->terms[j].coeff, d*req_images);
        Z->terms[j].coeff->length = 0;
    }
}

int qzip_solve(
    fq_nmod_mpolyu_t A,
    n_fq_polyun_t Z,
    n_fq_polyun_t H,
    n_fq_polyun_t M,
    const fq_nmod_mpoly_ctx_t ctx)
{
    int success;
    slong i, n;
    n_fq_poly_t t;
    n_fq_polyun_term_struct * Zterms = Z->terms;
    n_polyun_term_struct * Hterms = H->terms;
    n_polyun_term_struct * Mterms = M->terms;

    n_fq_poly_init(t);

    FLINT_ASSERT(A->length == H->length);
    FLINT_ASSERT(A->length == M->length);
    FLINT_ASSERT(A->length == Z->length);

    for (i = 0; i < A->length; i++)
    {
        slong d = fq_nmod_ctx_degree(ctx->fqctx);
        n = A->coeffs[i].length;
        FLINT_ASSERT(Mterms[i].coeff->length == n + 1);
        FLINT_ASSERT(Zterms[i].coeff->length >= n);
        FLINT_ASSERT(Hterms[i].coeff->length == n);

        n_poly_fit_length(t, d*n);

        if (A->coeffs[i].coeffs_alloc < d*n)
        {
            slong new_alloc = FLINT_MAX(d*n, A->coeffs[i].coeffs_alloc);
            A->coeffs[i].coeffs_alloc = new_alloc;
            A->coeffs[i].coeffs = flint_realloc(A->coeffs[i].coeffs,
                                                  new_alloc*sizeof(mp_limb_t));
        }

        success = fq_nmod_zip_find_coeffs_new_fq_nmod(A->coeffs[i].coeffs,
                             Hterms[i].coeff->coeffs, n,
                             Zterms[i].coeff->coeffs, Zterms[i].coeff->length,
                             Mterms[i].coeff->coeffs, t->coeffs, ctx->fqctx);
        if (success < 1)
        {
            n_fq_poly_clear(t);
            return success;
        }
    }

    n_fq_poly_clear(t);
    return 1;
}

#if 0
static int qzip_solvep(
    fq_nmod_mpolyu_t A,
    n_fq_polyun_t Z,
    n_fq_polyun_t H,
    n_fq_polyun_t M,
    const fq_nmod_mpoly_ctx_t ctx)
{
    int success;
    slong i, n;
    n_fq_poly_t t;
    n_fq_polyun_term_struct * Zterms = Z->terms;
    n_polyun_term_struct * Hterms = H->terms;
    n_polyun_term_struct * Mterms = M->terms;

    n_fq_poly_init(t);

    FLINT_ASSERT(A->length == H->length);
    FLINT_ASSERT(A->length == M->length);
    FLINT_ASSERT(A->length == Z->length);

    for (i = 0; i < A->length; i++)
    {
        slong d = fq_nmod_ctx_degree(ctx->fqctx);
        n = A->coeffs[i].length;
        FLINT_ASSERT(Mterms[i].coeff->length == n + 1);
        FLINT_ASSERT(Zterms[i].coeff->length >= n);
        FLINT_ASSERT(Hterms[i].coeff->length == n);

        n_poly_fit_length(t, d*n);

        success = fq_nmod_zip_findp_coeffs(A->coeffs[i].coeffs,
                             Hterms[i].coeff->coeffs, n,
                             Zterms[i].coeff->coeffs, Zterms[i].coeff->length,
                             Mterms[i].coeff->coeffs, t->coeffs, ctx->fqctx);
        if (success < 1)
        {
            n_fq_poly_clear(t);
            return success;
        }
    }

    n_fq_poly_clear(t);
    return 1;
}
#endif

void n_fq_bpoly_scalar_mul_n_fq(
    n_fq_bpoly_t A,
    const mp_limb_t * c,
    const fq_nmod_ctx_t ctx)
{
    slong d = fq_nmod_ctx_degree(ctx);
    slong i;

    if (_n_fq_is_zero(c, d))
    {
        A->length = 0;
        return;
    }

    if (_n_fq_is_one(c, d))
    {
        return;
    }

    for (i = 0; i < A->length; i++)
        n_fq_poly_scalar_mul_n_fq(A->coeffs + i, A->coeffs + i, c, ctx);
}

void n_fq_poly_eval_step(
    mp_limb_t * res,
    n_fq_poly_t A,
    const fq_nmod_ctx_t ctx)
{
    slong d = fq_nmod_ctx_degree(ctx);
    slong i, Alen = A->length;
    mp_limb_t * Acoeffs = A->coeffs;
    mp_limb_t * tmp, * sum;
    TMP_INIT;

    FLINT_ASSERT(d*3*Alen <= A->alloc);

    if (Alen < 1)
    {
        _n_fq_zero(res, d);
        return;
    }

    TMP_START;
    tmp = (mp_limb_t *) TMP_ALLOC(8*d*sizeof(mp_limb_t));
    sum = tmp + 4*d;

    i = 0;
    _n_fq_mul2(sum, Acoeffs + d*(3*i + 0), Acoeffs + d*(3*i + 1), ctx);
    _n_fq_mul(Acoeffs + d*(3*i + 0), Acoeffs + d*(3*i + 0), Acoeffs + d*(3*i + 2), ctx, tmp);
    for (i = 1; i < Alen; i++)
    {
        _n_fq_madd2(sum, Acoeffs + d*(3*i + 0), Acoeffs + d*(3*i + 1), ctx, tmp);
        _n_fq_mul(Acoeffs + d*(3*i + 0), Acoeffs + d*(3*i + 0), Acoeffs + d*(3*i + 2), ctx, tmp);
    }
    _n_fq_reduce2(res, sum, ctx, tmp);

    TMP_END;
}

void n_fq_poly_evalp_step(
    mp_limb_t * res,
    n_fq_poly_t A,
    const fq_nmod_ctx_t ctx)
{
    slong d = fq_nmod_ctx_degree(ctx);
    nmod_t mod = fq_nmod_ctx_mod(ctx);
    slong i, j, Alen = A->length;
    mp_limb_t * Acoeffs = A->coeffs;
    mp_limb_t * tmp, p1, p0;
    TMP_INIT;

    FLINT_ASSERT(d*3*Alen <= A->alloc);

    if (Alen < 1)
    {
        _n_fq_zero(res, d);
        return;
    }

    /*
        3i+0 current power Fp
        3i+1 coeff         Fq
        3i+2 incr power    Fp
    */

    if (d == 2)
    {
        mp_limb_t tmp[6];

        i = 0;

        for (j = 0; j < d; j++)
        {
            umul_ppmm(tmp[3*j + 1], tmp[3*j + 0], (Acoeffs + d*(3*i + 0))[0],
                                                  (Acoeffs + d*(3*i + 1))[j]);
            tmp[3*j + 2] = 0;
        }
        (Acoeffs + d*(3*i + 0))[0] = nmod_mul((Acoeffs + d*(3*i + 0))[0],
                                              (Acoeffs + d*(3*i + 2))[0], mod);
        for (i = 1; i < Alen; i++)
        {
            for (j = 0; j < d; j++)
            {
                umul_ppmm(p1, p0, (Acoeffs + d*(3*i + 0))[0],
                                  (Acoeffs + d*(3*i + 1))[j]);
                add_sssaaaaaa(tmp[3*j + 2], tmp[3*j + 1], tmp[3*j + 0],
                              tmp[3*j + 2], tmp[3*j + 1], tmp[3*j + 0], 0, p1, p0);
            }
            (Acoeffs + d*(3*i + 0))[0] = nmod_mul((Acoeffs + d*(3*i + 0))[0],
                                                  (Acoeffs + d*(3*i + 2))[0], mod);
        }

        for (j = 0; j < d; j++)
        {
            NMOD_RED3(res[j], tmp[3*j + 2], tmp[3*j + 1], tmp[3*j + 0], mod);
        }
    }
    else
    {
        TMP_START;
        tmp = (mp_limb_t *) TMP_ALLOC(3*d*sizeof(mp_limb_t));

        i = 0;

        for (j = 0; j < d; j++)
        {
            umul_ppmm(tmp[3*j + 1], tmp[3*j + 0], (Acoeffs + d*(3*i + 0))[0],
                                                  (Acoeffs + d*(3*i + 1))[j]);
            tmp[3*j + 2] = 0;
        }
        (Acoeffs + d*(3*i + 0))[0] = nmod_mul((Acoeffs + d*(3*i + 0))[0],
                                              (Acoeffs + d*(3*i + 2))[0], mod);
        for (i = 1; i < Alen; i++)
        {
            for (j = 0; j < d; j++)
            {
                umul_ppmm(p1, p0, (Acoeffs + d*(3*i + 0))[0],
                                  (Acoeffs + d*(3*i + 1))[j]);
                add_sssaaaaaa(tmp[3*j + 2], tmp[3*j + 1], tmp[3*j + 0],
                              tmp[3*j + 2], tmp[3*j + 1], tmp[3*j + 0], 0, p1, p0);
            }
            (Acoeffs + d*(3*i + 0))[0] = nmod_mul((Acoeffs + d*(3*i + 0))[0],
                                                  (Acoeffs + d*(3*i + 2))[0], mod);
        }

        for (j = 0; j < d; j++)
        {
            NMOD_RED3(res[j], tmp[3*j + 2], tmp[3*j + 1], tmp[3*j + 0], mod);
        }

        TMP_END;
    }
}


void n_fq_bpoly_eval_step(
    n_fq_bpoly_t E,
    n_fq_polyun_t A,
    const fq_nmod_ctx_t ctx)
{
    slong d = fq_nmod_ctx_degree(ctx);
    slong Ai;
    ulong e0, e1;
    n_fq_polyun_term_struct * Aterms = A->terms;
    mp_limb_t * c = FLINT_ARRAY_ALLOC(d, mp_limb_t);

    n_fq_bpoly_zero(E);
    for (Ai = 0; Ai < A->length; Ai++)
    {
        n_fq_poly_eval_step(c, Aterms[Ai].coeff, ctx);
        e0 = extract_exp(Aterms[Ai].exp, 1, 2);
        e1 = extract_exp(Aterms[Ai].exp, 0, 2);
        if (_n_fq_is_zero(c, d))
            continue;
        n_fq_bpoly_set_coeff_n_fq(E, e0, e1, c, ctx);
    }

    flint_free(c);
}


void n_fq_bpoly_evalp_step(
    n_fq_bpoly_t E,
    n_fq_polyun_t A,
    const fq_nmod_ctx_t ctx)
{
    slong d = fq_nmod_ctx_degree(ctx);
    slong Ai;
    ulong e0, e1;
    n_fq_polyun_term_struct * Aterms = A->terms;
    mp_limb_t * c = FLINT_ARRAY_ALLOC(d, mp_limb_t);

    n_fq_bpoly_zero(E);
    for (Ai = 0; Ai < A->length; Ai++)
    {
        n_fq_poly_evalp_step(c, Aterms[Ai].coeff, ctx);
        e0 = extract_exp(Aterms[Ai].exp, 1, 2);
        e1 = extract_exp(Aterms[Ai].exp, 0, 2);
        if (_n_fq_is_zero(c, d))
            continue;
        n_fq_bpoly_set_coeff_n_fq(E, e0, e1, c, ctx);
    }

    flint_free(c);
}

int n_fq_polyu2_add_zip_must_match(
    n_fq_polyun_t Z,
    const n_fq_bpoly_t A,
    slong cur_length,
    const fq_nmod_ctx_t ctx)
{
    slong d = fq_nmod_ctx_degree(ctx);
    slong i, Ai, ai;
    n_polyun_term_struct * Zt = Z->terms;
    const n_poly_struct * Acoeffs = A->coeffs;

    Ai = A->length - 1;
    ai = (Ai < 0) ? 0 : n_fq_poly_degree(A->coeffs + Ai);

    for (i = 0; i < Z->length; i++)
    {
        if (Ai >= 0 && Zt[i].exp == pack_exp2(Ai, ai))
        {
            /* Z present, A present */
            _n_fq_set(Zt[i].coeff->coeffs + d*cur_length, Acoeffs[Ai].coeffs + d*ai, d);
            Zt[i].coeff->length = cur_length + 1;
            do {
                ai--;
            } while (ai >= 0 && _n_fq_is_zero(Acoeffs[Ai].coeffs + d*ai, d));
            if (ai < 0)
            {
                do {
                    Ai--;
                } while (Ai >= 0 && Acoeffs[Ai].length == 0);
                if (Ai >= 0)
                    ai = n_fq_poly_degree(Acoeffs + Ai);
            }
        }
        else if (Ai < 0 || Zt[i].exp > pack_exp2(Ai, ai))
        {
            /* Z present, A missing */
            _n_fq_zero(Zt[i].coeff->coeffs + d*cur_length, d);
            Zt[i].coeff->length = cur_length + 1;
        }
        else
        {
            /* Z missing, A present */
            return 0;
        }
    }

    return 1;
}

/* F = F + modulus*(A - F(alpha))*/
int newfq_nmod_mpolyn_interp_crt_sm_mpoly(
    slong * lastdeg_,
    fq_nmod_mpolyn_t F,
    fq_nmod_mpolyn_t T,
    fq_nmod_mpoly_t A,  /* could have zero coeffs */
    const n_fq_poly_t modulus,
    n_fq_poly_t alphapow,
    const fq_nmod_mpoly_ctx_t ctx)
{
    slong lastdeg = *lastdeg_;
    slong d = fq_nmod_ctx_degree(ctx->fqctx);
    flint_bitcnt_t bits = A->bits;
    slong N = mpoly_words_per_exp(bits, ctx->minfo);
    int changed = 0;
    slong i, j, k;
    mp_limb_t * v = FLINT_ARRAY_ALLOC(d, mp_limb_t);
    slong Flen = F->length, Alen = A->length;
    ulong * Fexp = F->exps, * Aexp = A->exps;
    ulong * Texp;
    mp_limb_t * Acoeff = A->coeffs;
    n_fq_poly_struct * Fcoeff = F->coeffs;
    n_fq_poly_struct * Tcoeff;

    FLINT_ASSERT(F->bits == bits);
    FLINT_ASSERT(ctx->minfo->ord == ORD_LEX);

    fq_nmod_mpolyn_fit_length(T, Flen + Alen, ctx);
    Texp = T->exps;
    Tcoeff = T->coeffs;

    i = j = k = 0;
    while (i < Flen || j < Alen)
    {
        if (i < Flen && j < Alen &&
                               mpoly_monomial_equal(Fexp + N*i, Aexp + N*j, N))
        {
            FLINT_ASSERT(!n_fq_poly_is_zero(Fcoeff + i));
            FLINT_ASSERT(n_fq_poly_degree(Fcoeff + i) < n_fq_poly_degree(modulus));

            /* F term ok, A term ok */
            n_fq_poly_eval_pow(v, Fcoeff + i, alphapow, ctx->fqctx);
            n_fq_sub(v, Acoeff + d*j, v, ctx->fqctx);
            if (!_n_fq_is_zero(v, d))
            {
                changed = 1;
                n_fq_poly_scalar_addmul_n_fq(Tcoeff + k, Fcoeff + i,
                                                       modulus, v, ctx->fqctx);
            }
            else
            {
                n_fq_poly_set(Tcoeff + k, Fcoeff + i, ctx->fqctx);
            }

            FLINT_ASSERT(!n_fq_poly_is_zero(Tcoeff + k));
            lastdeg = FLINT_MAX(lastdeg, n_fq_poly_degree(Tcoeff + k));

            mpoly_monomial_set(Texp + N*k, Aexp + N*j, N);
            k++;
            i++;
            j++;
        }
        else if (i < Flen && (j >= Alen ||
                          mpoly_monomial_gt_nomask(Fexp + N*i, Aexp + N*j, N)))
        {
            FLINT_ASSERT(!n_fq_poly_is_zero(Fcoeff + i));
            FLINT_ASSERT(n_fq_poly_degree(Fcoeff + i) < n_fq_poly_degree(modulus));

            /* F term ok, A term missing */
            n_fq_poly_eval_pow(v, Fcoeff + i, alphapow, ctx->fqctx);
            if (!_n_fq_is_zero(v, d))
            {
                changed = 1;
                _n_fq_neg(v, v, d, ctx->fqctx->mod);
                n_fq_poly_scalar_addmul_n_fq(Tcoeff + k, Fcoeff + i, modulus, v, ctx->fqctx);
            }
            else
            {
                n_fq_poly_set(Tcoeff + k, Fcoeff + i, ctx->fqctx);                
            }

            FLINT_ASSERT(!n_fq_poly_is_zero(Tcoeff + k));
            lastdeg = FLINT_MAX(lastdeg, n_fq_poly_degree(Tcoeff + k));

            mpoly_monomial_set(Texp + N*k, Fexp + N*i, N);
            k++;
            i++;
        }
        else
        {
            FLINT_ASSERT(j < Alen && (i >= Flen ||
                         mpoly_monomial_lt_nomask(Fexp + N*i, Aexp + N*j, N)));

            /* F term missing, A term ok */
            if (!_n_fq_is_zero(Acoeff + d*j, d))
            {
                changed = 1;
                n_fq_poly_scalar_mul_n_fq(Tcoeff + k,
                                            modulus, Acoeff + d*j, ctx->fqctx);

                FLINT_ASSERT(!n_fq_poly_is_zero(Tcoeff + k));
                lastdeg = FLINT_MAX(lastdeg, n_fq_poly_degree(Tcoeff + k));
                mpoly_monomial_set(Texp + N*k, Aexp + N*j, N);
                k++;
            }
            j++;
        }
    }
    T->length = k;

    if (changed)
        fq_nmod_mpolyn_swap(T, F);

    flint_free(v);

    *lastdeg_ = lastdeg;
    return changed;
}


int newfq_nmod_mpolyun_interp_crt_sm_mpolyu(
    slong * lastdeg,
    fq_nmod_mpolyun_t F,
    fq_nmod_mpolyun_t T,
    fq_nmod_mpolyu_t A, /* could have zero coeffs */
    const n_fq_poly_t modulus,
    n_fq_poly_t alphapow,
    const fq_nmod_mpoly_ctx_t ctx)
{
    int changed = 0;
    slong i, j, k;
    ulong * Texp;
    ulong * Fexp;
    ulong * Aexp;
    slong Flen;
    slong Alen;
    fq_nmod_mpolyn_t S;
    fq_nmod_mpolyn_struct * Tcoeff;
    fq_nmod_mpolyn_struct * Fcoeff;
    fq_nmod_mpoly_struct  * Acoeff;
    fq_nmod_mpoly_t zero;

    *lastdeg = -WORD(1);

    FLINT_ASSERT(F->bits == T->bits);
    FLINT_ASSERT(T->bits == A->bits);

    fq_nmod_mpolyn_init(S, F->bits, ctx);

    Flen = F->length;
    Alen = A->length;
    fq_nmod_mpolyun_fit_length(T, Flen + Alen, ctx);

    Tcoeff = T->coeffs;
    Fcoeff = F->coeffs;
    Acoeff = A->coeffs;
    Texp = T->exps;
    Fexp = F->exps;
    Aexp = A->exps;   

    fq_nmod_mpoly_init(zero, ctx);
    fq_nmod_mpoly_fit_length_reset_bits(zero, 0, A->bits, ctx);

    i = j = k = 0;
    while (i < Flen || j < Alen)
    {
        if (i < Flen && (j >= Alen || Fexp[i] > Aexp[j]))
        {
            /* F term ok, A term missing */
            fq_nmod_mpolyn_set(Tcoeff + k, Fcoeff + i, ctx);
            changed |= newfq_nmod_mpolyn_interp_crt_sm_mpoly(lastdeg,
                                  Tcoeff + k, S, zero, modulus, alphapow, ctx);
            Texp[k] = Fexp[i];
            k++;
            i++;
        }
        else if (j < Alen && (i >= Flen || Aexp[j] > Fexp[i]))
        {
            /* F term missing, A term ok */
            fq_nmod_mpolyn_zero(Tcoeff + k, ctx);
            changed |= newfq_nmod_mpolyn_interp_crt_sm_mpoly(lastdeg,
                            Tcoeff + k, S, Acoeff + j, modulus, alphapow, ctx);
            Texp[k] = Aexp[j];
            k++;
            j++;
        }
        else if (i < Flen && j < Alen && (Fexp[i] == Aexp[j]))
        {
            /* F term ok, A term ok */
            fq_nmod_mpolyn_set(Tcoeff + k, Fcoeff + i, ctx);
            changed |= newfq_nmod_mpolyn_interp_crt_sm_mpoly(lastdeg,
                            Tcoeff + k, S, Acoeff + j, modulus, alphapow, ctx);
            Texp[k] = Aexp[j];
            FLINT_ASSERT(!fq_nmod_mpolyn_is_zero(Tcoeff + k, ctx));
            k++;
            i++;
            j++;
        }
        else 
        {
            FLINT_ASSERT(0);
        }
    }
    T->length = k;

    if (changed)
        fq_nmod_mpolyun_swap(T, F);

    fq_nmod_mpolyn_clear(S, ctx);
    fq_nmod_mpoly_clear(zero, ctx);

    return changed;
}



/* F = F + modulus*(A - F(alpha)) with matching monomials */
int fq_nmod_mpolyn_interp_mcrt_sm_mpoly(
    slong * lastdeg_,
    fq_nmod_mpolyn_t F,
    fq_nmod_mpoly_t A,  /* could have zero coeffs */
    const n_fq_poly_t modulus,
    n_fq_poly_t alphapow,
    const fq_nmod_mpoly_ctx_t ctx)
{
    slong d = fq_nmod_ctx_degree(ctx->fqctx);
#if FLINT_WANT_ASSERT
    slong N = mpoly_words_per_exp(A->bits, ctx->minfo);
#endif
    slong lastdeg = *lastdeg_;
    int changed = 0;
    mp_limb_t * v = FLINT_ARRAY_ALLOC(d, mp_limb_t);
    slong i, Alen = A->length;
    mp_limb_t * Acoeffs = A->coeffs;
    n_fq_poly_struct * Fcoeffs = F->coeffs;

    FLINT_ASSERT(F->bits == A->bits);
    FLINT_ASSERT(ctx->minfo->ord == ORD_LEX);

    FLINT_ASSERT(F->length == Alen);
    for (i = 0; i < Alen; i++)
    {
        FLINT_ASSERT(mpoly_monomial_equal(F->exps + N*i, A->exps + N*i, N));
        FLINT_ASSERT(!n_fq_poly_is_zero(Fcoeffs + i));
        FLINT_ASSERT(n_fq_poly_degree(Fcoeffs + i) < n_fq_poly_degree(modulus));

        n_fq_poly_eval_pow(v, Fcoeffs + i, alphapow, ctx->fqctx);
        n_fq_sub(v, Acoeffs + d*i, v, ctx->fqctx);
        if (!_n_fq_is_zero(v, d))
        {
            changed = 1;
            n_fq_poly_scalar_addmul_n_fq(Fcoeffs + i, Fcoeffs + i,
                                                       modulus, v, ctx->fqctx);
        }

        FLINT_ASSERT(!n_fq_poly_is_zero(Fcoeffs + i));
        lastdeg = FLINT_MAX(lastdeg, n_fq_poly_degree(Fcoeffs + i));
    }

    flint_free(v);

    *lastdeg_ = lastdeg;
    return changed;
}


int fq_nmod_mpolyun_interp_mcrt_sm_mpolyu(
    slong * lastdeg,
    fq_nmod_mpolyun_t F,
    fq_nmod_mpolyu_t A, /* could have zero coeffs */
    const n_fq_poly_t modulus,
    n_fq_poly_t alphapow,
    const fq_nmod_mpoly_ctx_t ctx)
{
    int changed = 0;
    slong i, Alen = A->length;
    fq_nmod_mpolyn_struct * Fcoeffs = F->coeffs;
    fq_nmod_mpoly_struct  * Acoeffs = A->coeffs;

    FLINT_ASSERT(F->bits == A->bits);

    *lastdeg = -WORD(1);

    FLINT_ASSERT(F->length == Alen);
    for (i = 0; i < Alen; i++)
    {
        FLINT_ASSERT(F->exps[i] == A->exps[i]);
        changed |= fq_nmod_mpolyn_interp_mcrt_sm_mpoly(lastdeg,
                             Fcoeffs + i, Acoeffs + i, modulus, alphapow, ctx);
    }
/*
flint_printf("fq_nmod_mpolyun_interp_mcrt_sm_mpolyu returning %d\n", changed);
*/

    return changed;
}


#define USE_G    1
#define USE_ABAR 2
#define USE_BBAR 4

static double interp_cost(
    double degG,
    double numABgamma,
    double maxnumci,
    double totnumci,
    double degxAB,
    double degyAB)
{
    return degG*(degG*totnumci + numABgamma + 0.01*maxnumci*(
                     numABgamma + totnumci + (degxAB*degyAB)*(degxAB*degyAB)));
}

int mpoly_gcd_get_use_first(
    slong rGdeg,
    slong Adeg,
    slong Bdeg,
    slong gammadeg);


int fq_nmod_mpoly_gcd_get_use(
    slong rGdeg,
    slong Adeg,
    slong Bdeg,
    slong gammadeg,
    slong degxAB,
    slong degyAB,
    slong numABgamma,
    const fq_nmod_mpolyu_t G,
    const fq_nmod_mpolyu_t Abar,
    const fq_nmod_mpolyu_t Bbar)
{
    int use = 0;
    slong i, lower = FLINT_MAX(gammadeg, rGdeg);
    slong upper = gammadeg + FLINT_MIN(FLINT_MIN(Adeg, Bdeg), rGdeg);
    if (lower <= upper)
    {
        slong Gdeg = ((ulong)upper + (ulong)lower)/2;
        slong maxnumci, totnumci;
        double Gcost, Abarcost, Bbarcost;

        maxnumci = totnumci = 0;
        for (i = 0; i < G->length; i++)
        {
            maxnumci = FLINT_MAX(maxnumci, G->coeffs[i].length);
            totnumci += G->coeffs[i].length;
        }
        FLINT_ASSERT(Gdeg >= 0);
        Gcost = interp_cost(Gdeg,
                       numABgamma, maxnumci, totnumci, degxAB, degyAB);

        maxnumci = totnumci = 0;
        for (i = 0; i < Abar->length; i++)
        {
            maxnumci = FLINT_MAX(maxnumci, Abar->coeffs[i].length);
            totnumci += Abar->coeffs[i].length;
        }
        FLINT_ASSERT(gammadeg + Adeg - Gdeg >= 0);
        Abarcost = interp_cost(gammadeg + Adeg - Gdeg,
                       numABgamma, maxnumci, totnumci, degxAB, degyAB);

        maxnumci = totnumci = 0;
        for (i = 0; i < Bbar->length; i++)
        {
            maxnumci = FLINT_MAX(maxnumci, Bbar->coeffs[i].length);
            totnumci += Bbar->coeffs[i].length;
        }
        FLINT_ASSERT(gammadeg + Bdeg - Gdeg >= 0);
        Bbarcost = interp_cost(gammadeg + Bdeg - Gdeg,
                       numABgamma, maxnumci, totnumci, degxAB, degyAB);

        if (Gcost <= FLINT_MIN(Abarcost, Bbarcost)*1.125)
            use |= USE_G;

        if (Abarcost <= FLINT_MIN(Gcost, Bbarcost)*1.125)
            use |= USE_ABAR;

        if (Bbarcost <= FLINT_MIN(Gcost, Abarcost)*1.125)
            use |= USE_BBAR;
    }

    if (use == 0)
        use = USE_G | USE_ABAR | USE_BBAR;

    return use;
}


void n_fq_bpoly_print_pretty(
    const n_fq_bpoly_t A,
    const char * xvar,
    const char * yvar,
    const fq_nmod_ctx_t ctx);

ulong mpoly_bivar_degrees1_chained(ulong deg, ulong * exps, slong len)
{
    slong i;
    ulong mask = mpoly_overflow_mask_sp(FLINT_BITS/2);
    for (i = 0; i < len; i++)
        deg = mpoly_monomial_max1(deg, exps[i], FLINT_BITS/2, mask);
    return deg;
}

slong fq_nmod_mpolyu_total_length(const fq_nmod_mpolyu_t A)
{
    slong i, tot = 0;
    for (i = 0; i < A->length; i++)
        tot += A->coeffs[i].length;
    return tot;
}

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
        rAbar = Abar/lc(rG)
        rBbar = Bbar/lc(rG)
*/
int fq_nmod_mpolyl_gcd_zippel_smprime(
    fq_nmod_mpoly_t rG, const slong * rGdegs,
    fq_nmod_mpoly_t rAbar,
    fq_nmod_mpoly_t rBbar,
    const fq_nmod_mpoly_t A, const slong * Adegs,
    const fq_nmod_mpoly_t B, const slong * Bdegs,
    const fq_nmod_mpoly_t gamma, const slong * gammadegs,
    const fq_nmod_mpoly_ctx_t ctx)
{
    slong d = fq_nmod_ctx_degree(ctx->fqctx);
    int success, use, betas_in_fp, main_tries_left;
    slong i, m;
    slong nvars = ctx->minfo->nvars;
    flint_bitcnt_t bits = A->bits;
    fq_nmod_struct * alphas, * betas;
    mp_limb_t * betasp;
    flint_rand_t state;
    fq_nmod_mpoly_t cont;
    fq_nmod_mpoly_t T, G, Abar, Bbar;
    n_fq_polyun_t HG, HAbar, HBbar, MG, MAbar, MBbar, ZG, ZAbar, ZBbar;
    n_fq_bpoly_t Aev, Bev, Gev, Abarev, Bbarev;
    const mp_limb_t * gammaev;
    fq_nmod_mpolyn_t Tn, Gn, Abarn, Bbarn;
    slong lastdeg;
    slong cur_zip_image, req_zip_images, this_length;
    n_fq_polyun_t Aeh_cur, Aeh_inc, Beh_cur, Beh_inc;
    n_fq_poly_t gammaeh_cur, gammaeh_inc;
    n_fq_poly_t modulus, alphapow;
    fq_nmod_mpoly_struct * Aevals, * Bevals;
    fq_nmod_mpoly_struct * gammaevals;
    n_poly_bpoly_stack_t St;
    mp_limb_t * c;
    fq_nmod_t start_alpha;
    ulong GdegboundXY, newdegXY, Abideg, Bbideg;
    slong degxAB, degyAB;

flint_printf("fq_nmod_mpolyl_gcd_zippel_smprime called nvars = %wd\n", nvars);


/*
fq_nmod_ctx_print(ctx->fqctx);

flint_printf("A: ");
fq_nmod_mpoly_print_pretty(A, NULL, ctx);
flint_printf("\n");
flint_printf("B: ");
fq_nmod_mpoly_print_pretty(B, NULL, ctx);
flint_printf("\n");
*/

    FLINT_ASSERT(bits <= FLINT_BITS);
    FLINT_ASSERT(bits == A->bits);
    FLINT_ASSERT(bits == B->bits);
    FLINT_ASSERT(bits == gamma->bits);
    FLINT_ASSERT(A->length > 0);
    FLINT_ASSERT(B->length > 0);
    FLINT_ASSERT(gamma->length > 0);

    fq_nmod_mpoly_fit_length_reset_bits(rG, 1, bits, ctx);
    fq_nmod_mpoly_fit_length_reset_bits(rAbar, 1, bits, ctx);
    fq_nmod_mpoly_fit_length_reset_bits(rBbar, 1, bits, ctx);

#if FLINT_WANT_ASSERT
    {
        slong * tmp_degs = FLINT_ARRAY_ALLOC(nvars, slong);

        fq_nmod_mpoly_degrees_si(tmp_degs, A, ctx);
        for (i = 0; i < nvars; i++)
            FLINT_ASSERT(tmp_degs[i] == Adegs[i]);

        fq_nmod_mpoly_degrees_si(tmp_degs, B, ctx);
        for (i = 0; i < nvars; i++)
            FLINT_ASSERT(tmp_degs[i] == Bdegs[i]);

        fq_nmod_mpoly_degrees_si(tmp_degs, gamma, ctx);
        for (i = 0; i < nvars; i++)
            FLINT_ASSERT(tmp_degs[i] == gammadegs[i]);

        flint_free(tmp_degs);
    }
#endif

    fq_nmod_init(start_alpha, ctx->fqctx);
    c = FLINT_ARRAY_ALLOC(d, mp_limb_t);

    n_fq_polyun_init(HG);
    n_fq_polyun_init(HAbar);
    n_fq_polyun_init(HBbar);
    n_fq_polyun_init(MG);
    n_fq_polyun_init(MAbar);
    n_fq_polyun_init(MBbar);
    n_fq_polyun_init(ZG);
    n_fq_polyun_init(ZAbar);
    n_fq_polyun_init(ZBbar);
    n_fq_bpoly_init(Aev);
    n_fq_bpoly_init(Bev);
    n_fq_bpoly_init(Gev);
    n_fq_bpoly_init(Abarev);
    n_fq_bpoly_init(Bbarev);
    n_fq_poly_init2(alphapow, 4, ctx->fqctx);
    fq_nmod_mpoly_init3(cont, 1, bits, ctx);
    fq_nmod_mpoly_init3(T, 0, bits, ctx);
    fq_nmod_mpoly_init3(G, 0, bits, ctx);
    fq_nmod_mpoly_init3(Abar, 0, bits, ctx);
    fq_nmod_mpoly_init3(Bbar, 0, bits, ctx);
    fq_nmod_mpolyn_init(Tn, bits, ctx);
    fq_nmod_mpolyn_init(Gn, bits, ctx);
    fq_nmod_mpolyn_init(Abarn, bits, ctx);
    fq_nmod_mpolyn_init(Bbarn, bits, ctx);
    n_fq_polyun_init(Aeh_cur);
    n_fq_polyun_init(Aeh_inc);
    n_fq_polyun_init(Beh_cur);
    n_fq_polyun_init(Beh_inc);
    n_fq_poly_init(gammaeh_cur);
    n_fq_poly_init(gammaeh_inc);
    n_fq_poly_init(modulus);
    n_poly_stack_init(St->poly_stack);
    n_bpoly_stack_init(St->bpoly_stack);

    betasp = FLINT_ARRAY_ALLOC(nvars, mp_limb_t);
    betas = FLINT_ARRAY_ALLOC(nvars, fq_nmod_struct);
    alphas = FLINT_ARRAY_ALLOC(nvars, fq_nmod_struct);
    for (i = 0; i < nvars; i++)
    {
        fq_nmod_init(betas + i, ctx->fqctx);
        fq_nmod_init(alphas + i, ctx->fqctx);
    }

    flint_randinit(state);

    Aevals = FLINT_ARRAY_ALLOC(nvars + 1, fq_nmod_mpoly_struct);
    Bevals = FLINT_ARRAY_ALLOC(nvars + 1, fq_nmod_mpoly_struct);
    gammaevals = FLINT_ARRAY_ALLOC(nvars + 1, fq_nmod_mpoly_struct);
    for (i = 0; i < nvars; i++)
    {
        fq_nmod_mpoly_init3(Aevals + i, 0, bits, ctx);
        fq_nmod_mpoly_init3(Bevals + i, 0, bits, ctx);
        fq_nmod_mpoly_init3(gammaevals + i, 0, bits, ctx);
    }
    Aevals[nvars] = *A;
    Bevals[nvars] = *B;
    gammaevals[nvars] = *gamma;

    Abideg = _fq_nmod_mpoly_bidegree(A, ctx);
    Bbideg = _fq_nmod_mpoly_bidegree(B, ctx);

    degxAB = FLINT_MAX(Adegs[0], Bdegs[0]);
    degyAB = FLINT_MAX(Adegs[1], Bdegs[1]);

    GdegboundXY = pack_exp2(FLINT_MIN(Adegs[0], Bdegs[0]),
                            FLINT_MIN(Adegs[1], Bdegs[1]));
    if (GdegboundXY == 0)
        goto gcd_is_trivial;

    main_tries_left = 3;

choose_main:
/*
flint_printf("-------- choose_main ----------\n");
*/
    if (--main_tries_left < 0)
    {
        success = 0;
        goto cleanup;
    }

    for (i = 2; i < nvars; i++)
{
        fq_nmod_rand_not_zero(alphas + i, state, ctx->fqctx);
/*
flint_printf("alphas[%wd]: ", i);
fq_nmod_print_pretty(alphas + i, ctx->fqctx);
flint_printf("\n");
*/
}

    for (i = nvars - 1; i >= 2; i--)
    {
        fq_nmod_mpoly_evaluate_one_fq_nmod(Aevals + i, Aevals + i + 1, i, alphas + i, ctx);
        fq_nmod_mpoly_repack_bits_inplace(Aevals + i, bits, ctx);
        fq_nmod_mpoly_evaluate_one_fq_nmod(Bevals + i, Bevals + i + 1, i, alphas + i, ctx);
        fq_nmod_mpoly_repack_bits_inplace(Bevals + i, bits, ctx);
        fq_nmod_mpoly_evaluate_one_fq_nmod(gammaevals + i, gammaevals + i + 1, i, alphas + i, ctx);
        fq_nmod_mpoly_repack_bits_inplace(gammaevals + i, bits, ctx);
        if (fq_nmod_mpoly_is_zero(gammaevals + i, ctx))
            goto choose_main;
        if (Aevals[i].length < 1 || _fq_nmod_mpoly_bidegree(Aevals + i, ctx) != Abideg)
            goto choose_main;
        if (Bevals[i].length < 1 || _fq_nmod_mpoly_bidegree(Bevals + i, ctx) != Bbideg)
            goto choose_main;
    }

    m = 2;
    if (rGdegs == NULL)
        use = USE_G | USE_ABAR | USE_BBAR;
    else
        use = mpoly_gcd_get_use_first(rGdegs[m], Adegs[m], Bdegs[m], gammadegs[m]);

    fq_nmod_mpoly_get_n_fq_bpoly(Aev, Aevals + m, 0, 1, ctx);
    fq_nmod_mpoly_get_n_fq_bpoly(Bev, Bevals + m, 0, 1, ctx);
/*
flint_printf("Aev: ");
n_fq_bpoly_print_pretty(Aev, "x1", "x2", ctx->fqctx);
flint_printf("\n");
flint_printf("Bev: ");
n_fq_bpoly_print_pretty(Bev, "x1", "x2", ctx->fqctx);
flint_printf("\n");
*/
    success = n_fq_bpoly_gcd_brown_smprime(Gev, Abarev, Bbarev,
                                                     Aev, Bev, ctx->fqctx, St);
/*
flint_printf("Gev: ");
n_fq_bpoly_print_pretty(Gev, "x1", "x2", ctx->fqctx);
flint_printf("\n");
flint_printf("Abarev: ");
n_fq_bpoly_print_pretty(Abarev, "x1", "x2", ctx->fqctx);
flint_printf("\n");
flint_printf("Bbarev: ");
n_fq_bpoly_print_pretty(Bbarev, "x1", "x2", ctx->fqctx);
flint_printf("\n");
*/

    if (!success)
        goto cleanup;

    newdegXY = n_bpoly_bidegree(Gev);
    if (newdegXY > GdegboundXY)
        goto choose_main;
    if (newdegXY < GdegboundXY)
    {
        GdegboundXY = newdegXY;
        if (GdegboundXY == 0)
            goto gcd_is_trivial;
    }

    gammaev = fq_nmod_mpoly_get_nonzero_n_fq(gammaevals + m, ctx);
    n_fq_bpoly_scalar_mul_n_fq(Gev, gammaev, ctx->fqctx);

    fq_nmod_mpolyn_interp_lift_sm_bpoly(Gn, Gev, ctx);
    fq_nmod_mpolyn_interp_lift_sm_bpoly(Abarn, Abarev, ctx);
    fq_nmod_mpolyn_interp_lift_sm_bpoly(Bbarn, Bbarev, ctx);

    n_fq_poly_one(modulus, ctx->fqctx);
    n_fq_set_fq_nmod(c, alphas + m, ctx->fqctx);
    n_fq_poly_shift_left_scalar_submul(modulus, 1, c, ctx->fqctx);

    fq_nmod_set(start_alpha, alphas + m, ctx->fqctx);
    while (1)
    {
    choose_alpha_2:

        fq_nmod_next_not_zero(alphas + m, ctx->fqctx);
/*
flint_printf("choosing alphas[%wd]: ", m);
fq_nmod_print_pretty(alphas + m, ctx->fqctx);
flint_printf("\n");
*/

        if (fq_nmod_equal(alphas + m, start_alpha, ctx->fqctx))
            goto choose_main;

        FLINT_ASSERT(alphapow->alloc >= d*2);
        alphapow->length = 2;
        _n_fq_one(alphapow->coeffs + d*0, d);
        n_fq_set_fq_nmod(alphapow->coeffs + d*1, alphas + m, ctx->fqctx);

        fq_nmod_mpoly_evaluate_one_fq_nmod(Aevals + m, Aevals + m + 1, m, alphas + m, ctx);
        fq_nmod_mpoly_repack_bits_inplace(Aevals + m, bits, ctx);
        fq_nmod_mpoly_evaluate_one_fq_nmod(Bevals + m, Bevals + m + 1, m, alphas + m, ctx);
        fq_nmod_mpoly_repack_bits_inplace(Bevals + m, bits, ctx);
        fq_nmod_mpoly_evaluate_one_fq_nmod(gammaevals + m, gammaevals + m + 1, m, alphas + m, ctx);
        fq_nmod_mpoly_repack_bits_inplace(gammaevals + m, bits, ctx);
        if (fq_nmod_mpoly_is_zero(gammaevals + m, ctx))
            goto choose_alpha_2;
        if (Aevals[m].length < 1 || _fq_nmod_mpoly_bidegree(Aevals + m, ctx) != Abideg)
            goto choose_alpha_2;
        if (Bevals[m].length < 1 || _fq_nmod_mpoly_bidegree(Bevals + m, ctx) != Bbideg)
            goto choose_alpha_2;

        fq_nmod_mpoly_get_n_fq_bpoly(Aev, Aevals + m, 0, 1, ctx);
        fq_nmod_mpoly_get_n_fq_bpoly(Bev, Bevals + m, 0, 1, ctx);

        success = n_fq_bpoly_gcd_brown_smprime(Gev, Abarev, Bbarev,
                                                     Aev, Bev, ctx->fqctx, St);
        if (!success)
            goto cleanup;

        newdegXY = n_bpoly_bidegree(Gev);
        if (newdegXY > GdegboundXY)
            goto choose_alpha_2;
        if (newdegXY < GdegboundXY)
        {
            GdegboundXY = newdegXY;
            if (GdegboundXY == 0)
                goto gcd_is_trivial;
            goto choose_main;
        }

        gammaev = fq_nmod_mpoly_get_nonzero_n_fq(gammaevals + m, ctx);
        n_fq_bpoly_scalar_mul_n_fq(Gev, gammaev, ctx->fqctx);

/*
flint_printf("Gev: ");
n_fq_bpoly_print_pretty(Gev, "x1", "x2", ctx->fqctx);
flint_printf("\n");
flint_printf("Abarev: ");
n_fq_bpoly_print_pretty(Abarev, "x1", "x2", ctx->fqctx);
flint_printf("\n");
flint_printf("Bbarev: ");
n_fq_bpoly_print_pretty(Bbarev, "x1", "x2", ctx->fqctx);
flint_printf("\n");


flint_printf("addind interpolant\n");
*/
        n_fq_poly_eval_pow(c, modulus, alphapow, ctx->fqctx);
        n_fq_inv(c, c, ctx->fqctx);
        n_fq_poly_scalar_mul_n_fq(modulus, modulus, c, ctx->fqctx);

        if ((use & USE_G) && !fq_nmod_mpolyn_interp_crt_sm_bpoly(
                                &lastdeg, Gn, Tn, Gev, modulus, alphapow, ctx))
        {
/*
flint_printf("G stabilized\n");
*/
            if (m == nvars - 1)
            {
                fq_nmod_mpoly_cvtfrom_mpolyn(rG, Gn, m, ctx);
                success = fq_nmod_mpolyl_content(cont, rG, 2, ctx);
                if (!success)
                    goto cleanup;
                fq_nmod_mpoly_divides(rG, rG, cont, ctx);
                fq_nmod_mpoly_repack_bits_inplace(rG, bits, ctx);
                if (fq_nmod_mpoly_divides(rAbar, A, rG, ctx) &&
                    fq_nmod_mpoly_divides(rBbar, B, rG, ctx))
                {
                    fq_nmod_mpoly_repack_bits_inplace(rAbar, bits, ctx);
                    fq_nmod_mpoly_repack_bits_inplace(rBbar, bits, ctx);
                    break;
                }                
            }
            else
            {
                fq_nmod_mpoly_cvtfrom_mpolyn(G, Gn, m, ctx);
                fq_nmod_mpoly_mul(T, Aevals + m + 1, gammaevals + m + 1, ctx);
                if (fq_nmod_mpoly_divides(Abar, T, G, ctx))
                {
                    fq_nmod_mpoly_repack_bits_inplace(Abar, bits, ctx);
                    fq_nmod_mpoly_mul(T, Bevals + m + 1, gammaevals + m + 1, ctx);
                    if (fq_nmod_mpoly_divides(Bbar, T, G, ctx))
                    {
                        fq_nmod_mpoly_repack_bits_inplace(Bbar, bits, ctx);
                        break;
                    }
                }
            }
        }

        if ((use & USE_ABAR) && !fq_nmod_mpolyn_interp_crt_sm_bpoly(
                          &lastdeg, Abarn, Tn, Abarev, modulus, alphapow, ctx))
        {
/*
flint_printf("Abar stabilized\n");
*/
            if (m == nvars - 1)
            {
                fq_nmod_mpoly_cvtfrom_mpolyn(rAbar, Abarn, m, ctx);
                success = fq_nmod_mpolyl_content(cont, rAbar, 2, ctx);
                if (!success)
                    goto cleanup;
                fq_nmod_mpoly_divides(rAbar, rAbar, cont, ctx);
                fq_nmod_mpoly_repack_bits_inplace(rAbar, bits, ctx);
                if (fq_nmod_mpoly_divides(rG, A, rAbar, ctx) &&
                    fq_nmod_mpoly_divides(rBbar, B, rG, ctx))
                {
                    fq_nmod_mpoly_repack_bits_inplace(rG, bits, ctx);
                    fq_nmod_mpoly_repack_bits_inplace(rBbar, bits, ctx);
                    break;
                }
            }
            else
            {
                fq_nmod_mpoly_cvtfrom_mpolyn(Abar, Abarn, m, ctx);
                fq_nmod_mpoly_mul(T, Aevals + m + 1, gammaevals + m + 1, ctx);
                if (fq_nmod_mpoly_divides(G, T, Abar, ctx))
                {
                    fq_nmod_mpoly_repack_bits_inplace(G, bits, ctx);
                    fq_nmod_mpoly_mul(T, Bevals + m + 1, gammaevals + m + 1, ctx);
                    if (fq_nmod_mpoly_divides(Bbar, T, G, ctx))
                    {
                        fq_nmod_mpoly_repack_bits_inplace(Bbar, bits, ctx);
                        break;
                    }
                }
            }
        }

        if ((use & USE_BBAR) && !fq_nmod_mpolyn_interp_crt_sm_bpoly(
                          &lastdeg, Bbarn, Tn, Bbarev, modulus, alphapow, ctx))
        {
/*
flint_printf("Bbar stabilized\n");
*/
            if (m == nvars - 1)
            {
                fq_nmod_mpoly_cvtfrom_mpolyn(rBbar, Bbarn, m, ctx);
                success = fq_nmod_mpolyl_content(cont, rBbar, 2, ctx);
                if (!success)
                    goto cleanup;
                fq_nmod_mpoly_divides(rBbar, rBbar, cont, ctx);
                fq_nmod_mpoly_repack_bits_inplace(rBbar, bits, ctx);
                if (fq_nmod_mpoly_divides(rG, B, rBbar, ctx) &&
                    fq_nmod_mpoly_divides(rAbar, A, rG, ctx))
                {
                    fq_nmod_mpoly_repack_bits_inplace(rG, bits, ctx);
                    fq_nmod_mpoly_repack_bits_inplace(rAbar, bits, ctx);
                    break;
                }
            }
            else
            {
                fq_nmod_mpoly_cvtfrom_mpolyn(Bbar, Bbarn, m, ctx);
                fq_nmod_mpoly_mul(T, Bevals + m + 1, gammaevals + m + 1, ctx);
                if (fq_nmod_mpoly_divides(G, T, Bbar, ctx))
                {
                    fq_nmod_mpoly_repack_bits_inplace(G, bits, ctx);
                    fq_nmod_mpoly_mul(T, Aevals + m + 1, gammaevals + m + 1, ctx);
                    if (fq_nmod_mpoly_divides(Abar, T, G, ctx))
                    {
                        fq_nmod_mpoly_repack_bits_inplace(Abar, bits, ctx);
                        break;
                    }
                }
            }
        }
/*
flint_printf("modulus: ");
n_fq_poly_print_pretty(modulus, "v", ctx->fqctx);
flint_printf("\n");
*/
        if (n_fq_poly_degree(modulus) > gammadegs[m] + Adegs[m] &&
            n_fq_poly_degree(modulus) > gammadegs[m] + Bdegs[m])
        {
            goto choose_main;
        }

        FLINT_ASSERT(n_fq_equal_fq_nmod(alphapow->coeffs + d*1, alphas + m, ctx->fqctx));
        n_fq_poly_shift_left_scalar_submul(modulus, 1, alphapow->coeffs + d*1, ctx->fqctx);
    }

    for (m = 3; m < nvars; m++)
    {
/*
flint_printf("m = %wd\n", m);
*/
        /* G, Abar, Bbar are in Fq[gen(0), ..., gen(m - 1)] */
        fq_nmod_mpolyn_interp_lift_sm_mpoly(Gn, G, ctx);
        fq_nmod_mpolyn_interp_lift_sm_mpoly(Abarn, Abar, ctx);
        fq_nmod_mpolyn_interp_lift_sm_mpoly(Bbarn, Bbar, ctx);

        if (rGdegs == NULL)
            use = USE_G | USE_ABAR | USE_BBAR;
        else
            use = 0;

        n_fq_poly_one(modulus, ctx->fqctx);
        n_fq_set_fq_nmod(c, alphas + m, ctx->fqctx);
        n_fq_poly_shift_left_scalar_submul(modulus, 1, c, ctx->fqctx);

        fq_nmod_set(start_alpha, alphas + m, ctx->fqctx);

    choose_betas:

        /* only beta[0], beta[1], ..., beta[m - 1] will be used */
        betas_in_fp = (ctx->fqctx->mod.norm < FLINT_BITS/4);
        if (betas_in_fp)
        {
            for (i = 2; i < ctx->minfo->nvars; i++)
                betasp[i] = n_urandint(state, ctx->fqctx->mod.n - 3) + 2;

            fq_nmod_mpoly_monomial_evalsp2(HG, G, betasp + 2, m, ctx);
            fq_nmod_mpoly_monomial_evalsp2(HAbar, Abar, betasp + 2, m, ctx);
            fq_nmod_mpoly_monomial_evalsp2(HBbar, Bbar, betasp + 2, m, ctx);
        }
        else
        {
            for (i = 2; i < ctx->minfo->nvars; i++)
                fq_nmod_rand_not_zero(betas + i, state, ctx->fqctx);

            FLINT_ASSERT(G->bits == bits);
            FLINT_ASSERT(Abar->bits == bits);
            FLINT_ASSERT(Bbar->bits == bits);

            fq_nmod_mpoly_monomial_evals2(HG, G, betas + 2, m, ctx);
            fq_nmod_mpoly_monomial_evals2(HAbar, Abar, betas + 2, m, ctx);
            fq_nmod_mpoly_monomial_evals2(HBbar, Bbar, betas + 2, m, ctx);
        }

        if (use == 0)
        {
            this_length = gammaevals[m + 1].length + Aevals[m + 1].length +
                                                     Bevals[m + 1].length;

            use = nmod_mpoly_gcd_get_use_new(rGdegs[m], Adegs[m], Bdegs[m],
                  gammadegs[m], degxAB, degyAB, this_length, HG, HAbar, HBbar);
        }

        req_zip_images = 1;
        if (use & USE_G)
        {
            this_length = betas_in_fp ?
                         n_polyun_product_roots(MG, HG, ctx->fqctx->mod) :
                         n_fq_polyun_product_roots(MG, HG, ctx->fqctx);
            req_zip_images = FLINT_MAX(req_zip_images, this_length);
        }
        if (use & USE_ABAR)
        {
            this_length = betas_in_fp ?
                         n_polyun_product_roots(MAbar, HAbar, ctx->fqctx->mod) :
                         n_fq_polyun_product_roots(MAbar, HAbar, ctx->fqctx);
            req_zip_images = FLINT_MAX(req_zip_images, this_length);
        }
        if (use & USE_BBAR)
        {
            this_length = betas_in_fp ?
                         n_polyun_product_roots(MBbar, HBbar, ctx->fqctx->mod) :
                         n_fq_polyun_product_roots(MBbar, HBbar, ctx->fqctx);
            req_zip_images = FLINT_MAX(req_zip_images, this_length);
        }

        while (1)
        {
        choose_alpha_m:

            fq_nmod_next_not_zero(alphas + m, ctx->fqctx);
/*
flint_printf("choosing alphas[%wd]: ", m);
fq_nmod_print_pretty(alphas + m, ctx->fqctx);
flint_printf("\n");
*/
            if (fq_nmod_equal(alphas + m, start_alpha, ctx->fqctx))
                goto choose_main;

            FLINT_ASSERT(alphapow->alloc >= d*2);
            alphapow->length = 2;
            _n_fq_one(alphapow->coeffs + d*0, d);
            n_fq_set_fq_nmod(alphapow->coeffs + d*1, alphas + m, ctx->fqctx);

            fq_nmod_mpoly_evaluate_one_fq_nmod(Aevals + m, Aevals + m + 1, m, alphas + m, ctx);
            fq_nmod_mpoly_repack_bits_inplace(Aevals + m, bits, ctx);
            fq_nmod_mpoly_evaluate_one_fq_nmod(Bevals + m, Bevals + m + 1, m, alphas + m, ctx);
            fq_nmod_mpoly_repack_bits_inplace(Bevals + m, bits, ctx);
            fq_nmod_mpoly_evaluate_one_fq_nmod(gammaevals + m, gammaevals + m + 1, m, alphas + m, ctx);
            fq_nmod_mpoly_repack_bits_inplace(gammaevals + m, bits, ctx);
            if (fq_nmod_mpoly_is_zero(gammaevals + m, ctx))
                goto choose_alpha_m;
            if (Aevals[m].length < 1 || _fq_nmod_mpoly_bidegree(Aevals + m, ctx) != Abideg)
                goto choose_alpha_m;
            if (Bevals[m].length < 1 || _fq_nmod_mpoly_bidegree(Bevals + m, ctx) != Bbideg)
                goto choose_alpha_m;

            if (betas_in_fp)
            {
                fq_nmod_mpoly_monomial_evalsp2(Aeh_inc, Aevals + m, betasp + 2, m, ctx);
                fq_nmod_mpoly_monomial_evalsp2(Beh_inc, Bevals + m, betasp + 2, m, ctx);
                fq_nmod_mpoly_monomial_evalsp(gammaeh_inc, gammaevals + m, betasp + 2, 2, m, ctx);
                n_polyun_set(Aeh_cur, Aeh_inc);
                n_polyun_set(Beh_cur, Beh_inc);
                n_poly_set(gammaeh_cur, gammaeh_inc);
            }
            else
            {
                fq_nmod_mpoly_monomial_evals2(Aeh_inc, Aevals + m, betas + 2, m, ctx);
                fq_nmod_mpoly_monomial_evals2(Beh_inc, Bevals + m, betas + 2, m, ctx);
                fq_nmod_mpoly_monomial_evals(gammaeh_inc, gammaevals + m, betas + 2, 2, m, ctx);
                n_fq_polyun_set(Aeh_cur, Aeh_inc, ctx->fqctx);
                n_fq_polyun_set(Beh_cur, Beh_inc, ctx->fqctx);
                n_fq_poly_set(gammaeh_cur, gammaeh_inc, ctx->fqctx);
            }

            qzip_start(ZG, HG, req_zip_images, ctx->fqctx);
            qzip_start(ZAbar, HAbar, req_zip_images, ctx->fqctx);
            qzip_start(ZBbar, HBbar, req_zip_images, ctx->fqctx);
            for (cur_zip_image = 0; cur_zip_image < req_zip_images; cur_zip_image++)
            {
                if (betas_in_fp)
                {
                    n_fq_bpoly_evalp_step_new(Aev, Aeh_cur, Aeh_inc, Aevals + m, ctx->fqctx);
                    n_fq_bpoly_evalp_step_new(Bev, Beh_cur, Beh_inc, Bevals + m, ctx->fqctx);
                    n_fq_poly_evalp_step_new(c, gammaeh_cur, gammaeh_inc, gammaevals + m, ctx->fqctx);
                }
                else
                {
                    n_fq_bpoly_eval_step_new(Aev, Aeh_cur, Aeh_inc, Aevals + m, ctx->fqctx);
                    n_fq_bpoly_eval_step_new(Bev, Beh_cur, Beh_inc, Bevals + m, ctx->fqctx);
                    n_fq_poly_eval_step_new(c, gammaeh_cur, gammaeh_inc, gammaevals + m, ctx->fqctx);
                }

                if (_n_fq_is_zero(c, d))
                    goto choose_betas;
                if (Aev->length < 1 || n_bpoly_bidegree(Aev) != Abideg)
                    goto choose_betas;
                if (Bev->length < 1 || n_bpoly_bidegree(Bev) != Bbideg)
                    goto choose_betas;

                success = n_fq_bpoly_gcd_brown_smprime(Gev, Abarev, Bbarev,
                                                     Aev, Bev, ctx->fqctx, St);        
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
                    goto choose_main;
                }

                n_fq_bpoly_scalar_mul_n_fq(Gev, c, ctx->fqctx);
                if ((use & USE_G) && !n_fq_polyu2_add_zip_must_match(ZG, Gev, cur_zip_image, ctx->fqctx))
                    goto choose_main;
                if ((use & USE_ABAR) && !n_fq_polyu2_add_zip_must_match(ZAbar, Abarev, cur_zip_image, ctx->fqctx))
                    goto choose_main;
                if ((use & USE_BBAR) && !n_fq_polyu2_add_zip_must_match(ZBbar, Bbarev, cur_zip_image, ctx->fqctx))
                    goto choose_main;
            }

            if (betas_in_fp)
            {
                if ((use & USE_G) && qzip_solvelp(G, ZG, HG, MG, ctx) < 1)
                    goto choose_main;
                if ((use & USE_ABAR) && qzip_solvelp(Abar, ZAbar, HAbar, MAbar, ctx) < 1)
                    goto choose_main;
                if ((use & USE_BBAR) && qzip_solvelp(Bbar, ZBbar, HBbar, MBbar, ctx) < 1)
                    goto choose_main;
            }
            else
            {
                if ((use & USE_G) && qzip_solvel(G, ZG, HG, MG, ctx) < 1)
                    goto choose_main;
                if ((use & USE_ABAR) && qzip_solvel(Abar, ZAbar, HAbar, MAbar, ctx) < 1)
                    goto choose_main;
                if ((use & USE_BBAR) && qzip_solvel(Bbar, ZBbar, HBbar, MBbar, ctx) < 1)
                    goto choose_main;
            }
/*
flint_printf("adding interpolant\n");
*/
            n_fq_poly_eval_pow(c, modulus, alphapow, ctx->fqctx);
            n_fq_inv(c, c, ctx->fqctx);
            n_fq_poly_scalar_mul_n_fq(modulus, modulus, c, ctx->fqctx);

            if ((use & USE_G) && !fq_nmod_mpolyn_interp_mcrt_sm_mpoly(
                                      &lastdeg, Gn, G, modulus, alphapow, ctx))
            {
/*
flint_printf("G stabilized\n");
*/
                fq_nmod_mpoly_cvtfrom_mpolyn(rG, Gn, m, ctx);
                if (m == nvars - 1)
                {
                    success = fq_nmod_mpolyl_content(cont, rG, 2, ctx);
                    if (!success)
                        goto cleanup;
                    fq_nmod_mpoly_divides(rG, rG, cont, ctx);
                    fq_nmod_mpoly_repack_bits_inplace(rG, bits, ctx);
                    if (fq_nmod_mpoly_divides(rAbar, A, rG, ctx) &&
                        fq_nmod_mpoly_divides(rBbar, B, rG,  ctx))
                    {
                        fq_nmod_mpoly_repack_bits_inplace(rAbar, bits, ctx);
                        fq_nmod_mpoly_repack_bits_inplace(rBbar, bits, ctx);
                        break;
                    }
                }
                else
                {
                    fq_nmod_mpoly_mul(T, Aevals + m + 1, gammaevals + m + 1, ctx);
                    if (fq_nmod_mpoly_divides(rAbar, T, rG, ctx))
                    {
                        fq_nmod_mpoly_mul(T, Bevals + m + 1, gammaevals + m + 1, ctx);
                        if (fq_nmod_mpoly_divides(rBbar, T, rG, ctx))
                        {
                            fq_nmod_mpoly_repack_bits_inplace(rAbar, bits, ctx);
                            fq_nmod_mpoly_repack_bits_inplace(rBbar, bits, ctx);
                            fq_nmod_mpoly_swap(G, rG, ctx);
                            fq_nmod_mpoly_swap(Abar, rAbar, ctx);
                            fq_nmod_mpoly_swap(Bbar, rBbar, ctx);
                            break;
                        }
                    }
                }
            }
            if ((use & USE_ABAR) && !fq_nmod_mpolyn_interp_mcrt_sm_mpoly(
                                &lastdeg, Abarn, Abar, modulus, alphapow, ctx))
            {
/*
flint_printf("Abar stabilized\n");
*/
                fq_nmod_mpoly_cvtfrom_mpolyn(rAbar, Abarn, m, ctx);
                if (m == nvars - 1)
                {
                    success = fq_nmod_mpolyl_content(cont, rAbar, 2, ctx);
                    if (!success)
                        goto cleanup;
                    success = fq_nmod_mpoly_divides(rAbar, rAbar, cont, ctx);
                    FLINT_ASSERT(success);
                    fq_nmod_mpoly_repack_bits_inplace(rAbar, bits, ctx);
                    if (fq_nmod_mpoly_divides(rG, A, rAbar, ctx) &&
                        fq_nmod_mpoly_divides(rBbar, B, rG, ctx))
                    {
                        fq_nmod_mpoly_repack_bits_inplace(rG, bits, ctx);
                        fq_nmod_mpoly_repack_bits_inplace(rBbar, bits, ctx);
                        break;
                    }
                }
                else
                {
                    fq_nmod_mpoly_mul(T, Aevals + m + 1, gammaevals + m + 1, ctx);
                    if (fq_nmod_mpoly_divides(rG, T, rAbar, ctx))
                    {
                        fq_nmod_mpoly_mul(T, Bevals + m + 1, gammaevals + m + 1, ctx);
                        if (fq_nmod_mpoly_divides(rBbar, T, rG, ctx))
                        {
                            fq_nmod_mpoly_repack_bits_inplace(rG, bits, ctx);
                            fq_nmod_mpoly_repack_bits_inplace(rBbar, bits, ctx);
                            fq_nmod_mpoly_swap(G, rG, ctx);
                            fq_nmod_mpoly_swap(Abar, rAbar, ctx);
                            fq_nmod_mpoly_swap(Bbar, rBbar, ctx);
                            break;
                        }
                    }
                }
            }
            if ((use & USE_BBAR) && !fq_nmod_mpolyn_interp_mcrt_sm_mpoly(
                                &lastdeg, Bbarn, Bbar, modulus, alphapow, ctx))
            {
/*
flint_printf("Bbar stabilized\n");
*/
                fq_nmod_mpoly_cvtfrom_mpolyn(rBbar, Bbarn, m, ctx);
                if (m == nvars - 1)
                {
                    success = fq_nmod_mpolyl_content(cont, rBbar, 2, ctx);
                    if (!success)
                        goto cleanup;
                    fq_nmod_mpoly_divides(rBbar, rBbar, cont, ctx);
                    fq_nmod_mpoly_repack_bits_inplace(rBbar, bits, ctx);
                    if (fq_nmod_mpoly_divides(rG, B, rBbar, ctx) &&
                        fq_nmod_mpoly_divides(rAbar, A, rG, ctx))
                    {
                        fq_nmod_mpoly_repack_bits_inplace(rG, bits, ctx);
                        fq_nmod_mpoly_repack_bits_inplace(rAbar, bits, ctx);
                        break;
                    }
                }
                else
                {
                    fq_nmod_mpoly_mul(T, Bevals + m + 1, gammaevals + m + 1, ctx);
                    if (fq_nmod_mpoly_divides(rG, T, rBbar, ctx))
                    {
                        fq_nmod_mpoly_mul(T, Aevals + m + 1, gammaevals + m + 1, ctx);
                        if (fq_nmod_mpoly_divides(rAbar, T, rG, ctx))
                        {
                            fq_nmod_mpoly_repack_bits_inplace(rG, bits, ctx);
                            fq_nmod_mpoly_repack_bits_inplace(rAbar, bits, ctx);
                            fq_nmod_mpoly_swap(G, rG, ctx);
                            fq_nmod_mpoly_swap(Abar, rAbar, ctx);
                            fq_nmod_mpoly_swap(Bbar, rBbar, ctx);
                            break;
                        }
                    }
                }
            }

            if (n_fq_poly_degree(modulus) > gammadegs[m] + Adegs[m] &&
                n_fq_poly_degree(modulus) > gammadegs[m] + Bdegs[m])
            {
                goto choose_main;
            }

            FLINT_ASSERT(n_fq_equal_fq_nmod(alphapow->coeffs + d*1, alphas + m, ctx->fqctx));
            n_fq_poly_shift_left_scalar_submul(modulus, 1, alphapow->coeffs + d*1, ctx->fqctx);
        }
    }

    success = 1;

cleanup:

    n_fq_polyun_clear(HG);
    n_fq_polyun_clear(HAbar);
    n_fq_polyun_clear(HBbar);
    n_fq_polyun_clear(MG);
    n_fq_polyun_clear(MAbar);
    n_fq_polyun_clear(MBbar);
    n_fq_polyun_clear(ZG);
    n_fq_polyun_clear(ZAbar);
    n_fq_polyun_clear(ZBbar);
    n_fq_bpoly_clear(Aev);
    n_fq_bpoly_clear(Bev);
    n_fq_bpoly_clear(Gev);
    n_fq_bpoly_clear(Abarev);
    n_fq_bpoly_clear(Bbarev);
    n_fq_poly_clear(alphapow);
    fq_nmod_mpoly_clear(cont, ctx);
    fq_nmod_mpoly_clear(T, ctx);
    fq_nmod_mpoly_clear(G, ctx);
    fq_nmod_mpoly_clear(Abar, ctx);
    fq_nmod_mpoly_clear(Bbar, ctx);
    fq_nmod_mpolyn_clear(Tn, ctx);
    fq_nmod_mpolyn_clear(Gn, ctx);
    fq_nmod_mpolyn_clear(Abarn, ctx);
    fq_nmod_mpolyn_clear(Bbarn, ctx);
    n_fq_polyun_clear(Aeh_cur);
    n_fq_polyun_clear(Aeh_inc);
    n_fq_polyun_clear(Beh_cur);
    n_fq_polyun_clear(Beh_inc);
    n_fq_poly_clear(gammaeh_cur);
    n_fq_poly_clear(gammaeh_inc);
    n_fq_poly_clear(modulus);
    n_poly_stack_clear(St->poly_stack);
    n_bpoly_stack_clear(St->bpoly_stack);

    fq_nmod_clear(start_alpha, ctx->fqctx);
    flint_free(c);

    flint_free(betasp);

    for (i = 0; i < nvars; i++)
    {
        fq_nmod_clear(betas + i, ctx->fqctx);
        fq_nmod_clear(alphas + i, ctx->fqctx);
    }
    flint_free(betas);
    flint_free(alphas);
    flint_randclear(state);

    for (i = 0; i < nvars; i++)
    {
        fq_nmod_mpoly_clear(Aevals + i, ctx);
        fq_nmod_mpoly_clear(Bevals + i, ctx);
        fq_nmod_mpoly_clear(gammaevals + i, ctx);
    }
    flint_free(Aevals);
    flint_free(Bevals);
    flint_free(gammaevals);

    FLINT_ASSERT(!success || rG->bits == bits);
    FLINT_ASSERT(!success || rAbar->bits == bits);
    FLINT_ASSERT(!success || rBbar->bits == bits);

flint_printf("fq_nmod_mpolyl_gcd_zippel_smprime returning %d\n", success);
/*
flint_printf("rG: ");
fq_nmod_mpoly_print_pretty(rG, NULL, ctx);
flint_printf("\n");

flint_printf("rAbar: ");
fq_nmod_mpoly_print_pretty(rAbar, NULL, ctx);
flint_printf("\n");

flint_printf("rBbar: ");
fq_nmod_mpoly_print_pretty(rBbar, NULL, ctx);
flint_printf("\n");
*/
FLINT_ASSERT(ctx->fqctx->mod.n < 1000 || success);

    return success;

gcd_is_trivial:

    fq_nmod_mpoly_one(rG, ctx);
    fq_nmod_mpoly_set(rAbar, A, ctx);
    fq_nmod_mpoly_set(rBbar, B, ctx);

    success = 1;
    
    goto cleanup;
}



int newfq_nmod_mpolyn_interp_mcrt_lg_mpoly(
    slong * lastdeg,
    fq_nmod_mpolyn_t H,
    const fq_nmod_mpoly_ctx_t ctx,
    const n_fq_poly_t m,
    const mp_limb_t * inv_m_eval,
    fq_nmod_mpoly_t A,
    const fq_nmod_mpoly_ctx_t ectx,
    const bad_fq_nmod_embed_t emb)
{
    slong lgd = fq_nmod_ctx_degree(ectx->fqctx);
    slong i;
#if FLINT_WANT_ASSERT
    slong N = mpoly_words_per_exp(A->bits, ctx->minfo);
#endif
    int changed = 0;
    mp_limb_t * u, * v, * tmp;
    n_fq_poly_struct * w, * u_sm;
    n_poly_stack_t St;

    FLINT_ASSERT(H->length == A->length);
    FLINT_ASSERT(H->bits == A->bits);

    n_poly_stack_init(St);
    n_poly_stack_fit_request(St, 3);

    w = n_poly_stack_take_top(St);
    u_sm = n_poly_stack_take_top(St);

    tmp = n_poly_stack_vec_init(St, lgd*(2 + N_FQ_MUL_ITCH));
    u = tmp + lgd*N_FQ_MUL_ITCH;
    v = u + lgd;

    for (i = 0; i < A->length; i++)
    {
        FLINT_ASSERT(mpoly_monomial_equal(H->exps + N*i, A->exps + N*i, N));
        bad_n_fq_embed_sm_to_lg(u, H->coeffs + i, emb);
        _n_fq_sub(v, A->coeffs + lgd*i, u, lgd, ectx->fqctx->mod);
        if (!_n_fq_is_zero(v, lgd))
        {
            changed = 1;
            _n_fq_mul(u, v, inv_m_eval, ectx->fqctx, tmp);
            bad_n_fq_embed_lg_to_sm(u_sm, u, emb);
            n_fq_poly_mul_(w, u_sm, m, ctx->fqctx, St);
            n_fq_poly_add(H->coeffs + i, H->coeffs + i, w, ctx->fqctx);
        }

        *lastdeg = FLINT_MAX(*lastdeg, n_fq_poly_degree(H->coeffs + i));

        FLINT_ASSERT(n_fq_poly_degree(H->coeffs + i) < n_fq_poly_degree(m) +
                                      fq_nmod_poly_degree(emb->h, ctx->fqctx));
    }

    n_poly_stack_vec_clear(St);

    n_poly_stack_give_back(St, 2);

    n_poly_stack_clear(St);

    return changed;
}



int newfq_nmod_mpolyun_interp_mcrt_lg_mpolyu(
    slong * lastdeg,
    fq_nmod_mpolyun_t H,
    const fq_nmod_mpoly_ctx_t ctx,
    n_fq_poly_t m,
    fq_nmod_mpolyu_t A,
    const fq_nmod_mpoly_ctx_t ectx,
    bad_fq_nmod_embed_t emb)
{
    fq_nmod_poly_t m_;
    int changed;
    fq_nmod_poly_init(m_, ctx->fqctx);
    n_fq_poly_get_fq_nmod_poly(m_, m, ctx->fqctx);
    changed = fq_nmod_mpolyun_interp_mcrt_lg_mpolyu(lastdeg, H, ctx, m_, A, ectx, emb);
    fq_nmod_poly_clear(m_, ctx->fqctx);
    return changed;
}    


int fq_nmod_mpolyl_gcd_zippel_lgprime(
    fq_nmod_mpoly_t rG, const slong * rGdegs,
    fq_nmod_mpoly_t rAbar,
    fq_nmod_mpoly_t rBbar,
    const fq_nmod_mpoly_t A, const slong * Adegs,
    const fq_nmod_mpoly_t B, const slong * Bdegs,
    const fq_nmod_mpoly_t gamma, const slong * gammadegs,
    const fq_nmod_mpoly_ctx_t smctx)
{
    slong lgd;
    int success, use;
    slong i, m;
    slong nvars = smctx->minfo->nvars;
    flint_bitcnt_t bits = A->bits;
    fq_nmod_struct * alphas, * betas;
    flint_rand_t state;
    fq_nmod_mpoly_t cont;
    fq_nmod_mpoly_t T, G, Abar, Bbar;
    n_fq_polyun_t HG, HAbar, HBbar, MG, MAbar, MBbar, ZG, ZAbar, ZBbar;
    n_fq_bpoly_t Aev, Bev, Gev, Abarev, Bbarev;
    const mp_limb_t * gammaev;
    fq_nmod_mpolyn_t Tn, Gn, Abarn, Bbarn;
    slong lastdeg;
    slong cur_zip_image, req_zip_images, this_length;
    n_fq_polyun_t Aeh_cur, Aeh_inc, Beh_cur, Beh_inc;
    n_fq_poly_t gammaeh_cur, gammaeh_inc;
    fq_nmod_mpoly_struct * Aevals, * Bevals;
    fq_nmod_mpoly_struct * gammaevals;
    n_fq_poly_t modulus, alphapow;
    n_poly_bpoly_stack_t St;
    fq_nmod_t start_alpha;
    n_poly_t tmp;  /* tmp arithmetic space */
    ulong GdegboundXY, newdegXY, Abideg, Bbideg;
    slong degxAB, degyAB;
    bad_fq_nmod_mpoly_embed_chooser_t embc;
    bad_fq_nmod_embed_struct * cur_emb;
    fq_nmod_mpoly_ctx_t lgctx;
    fq_nmod_mpolyn_t gamman;
    fq_nmod_mpolyn_t An, Bn;

flint_printf("fq_nmod_mpolyl_gcd_zippel_lgprime called nvars = %wd\n", nvars);

    FLINT_ASSERT(bits <= FLINT_BITS);
    FLINT_ASSERT(bits == A->bits);
    FLINT_ASSERT(bits == B->bits);
    FLINT_ASSERT(bits == gamma->bits);
    FLINT_ASSERT(A->length > 0);
    FLINT_ASSERT(B->length > 0);
    FLINT_ASSERT(gamma->length > 0);

    fq_nmod_mpoly_fit_length_reset_bits(rG, 1, bits, smctx);
    fq_nmod_mpoly_fit_length_reset_bits(rAbar, 1, bits, smctx);
    fq_nmod_mpoly_fit_length_reset_bits(rBbar, 1, bits, smctx);

#if FLINT_WANT_ASSERT
    {
        slong * tmp_degs = FLINT_ARRAY_ALLOC(nvars, slong);

        fq_nmod_mpoly_degrees_si(tmp_degs, A, smctx);
        for (i = 0; i < nvars; i++)
            FLINT_ASSERT(tmp_degs[i] == Adegs[i]);

        fq_nmod_mpoly_degrees_si(tmp_degs, B, smctx);
        for (i = 0; i < nvars; i++)
            FLINT_ASSERT(tmp_degs[i] == Bdegs[i]);

        fq_nmod_mpoly_degrees_si(tmp_degs, gamma, smctx);
        for (i = 0; i < nvars; i++)
            FLINT_ASSERT(tmp_degs[i] == gammadegs[i]);

        flint_free(tmp_degs);
    }
#endif

    n_fq_polyun_init(HG);
    n_fq_polyun_init(HAbar);
    n_fq_polyun_init(HBbar);
    n_fq_polyun_init(MG);
    n_fq_polyun_init(MAbar);
    n_fq_polyun_init(MBbar);
    n_fq_polyun_init(ZG);
    n_fq_polyun_init(ZAbar);
    n_fq_polyun_init(ZBbar);
    n_fq_bpoly_init(Aev);
    n_fq_bpoly_init(Bev);
    n_fq_bpoly_init(Gev);
    n_fq_bpoly_init(Abarev);
    n_fq_bpoly_init(Bbarev);
    fq_nmod_mpoly_init3(cont, 1, bits, smctx);
    fq_nmod_mpoly_init3(T, 0, bits, smctx);
    fq_nmod_mpoly_init3(G, 0, bits, smctx);
    fq_nmod_mpoly_init3(Abar, 0, bits, smctx);
    fq_nmod_mpoly_init3(Bbar, 0, bits, smctx);
    fq_nmod_mpolyn_init(Tn, bits, smctx);
    fq_nmod_mpolyn_init(Gn, bits, smctx);
    fq_nmod_mpolyn_init(Abarn, bits, smctx);
    fq_nmod_mpolyn_init(Bbarn, bits, smctx);
    n_fq_polyun_init(Aeh_cur);
    n_fq_polyun_init(Aeh_inc);
    n_fq_polyun_init(Beh_cur);
    n_fq_polyun_init(Beh_inc);
    n_fq_poly_init(gammaeh_cur);
    n_fq_poly_init(gammaeh_inc);
    n_fq_poly_init(modulus);
    n_poly_stack_init(St->poly_stack);
    n_bpoly_stack_init(St->bpoly_stack);
    fq_nmod_mpolyn_init(An, bits, smctx);
    fq_nmod_mpolyn_init(Bn, bits, smctx);
    fq_nmod_mpolyn_init(gamman, bits, smctx);

    /* alphas[nvars - 1] not used - it is replaced by cur_emb */
    betas = FLINT_ARRAY_ALLOC(nvars, fq_nmod_struct);
    alphas = FLINT_ARRAY_ALLOC(nvars, fq_nmod_struct);
    for (i = 0; i < nvars; i++)
    {
        fq_nmod_init(betas + i, smctx->fqctx);
        fq_nmod_init(alphas + i, smctx->fqctx);
    }

    flint_randinit(state);
    cur_emb = bad_fq_nmod_mpoly_embed_chooser_init(embc, lgctx, smctx, state);
    lgd = fq_nmod_ctx_degree(lgctx->fqctx);
    n_poly_init2(tmp, lgd);
    n_poly_init2(alphapow, 2*lgd);
    fq_nmod_init(start_alpha, lgctx->fqctx);

    /* Aevals[nvars] does not exist - it is replaced by An */
    Aevals = FLINT_ARRAY_ALLOC(nvars, fq_nmod_mpoly_struct);
    Bevals = FLINT_ARRAY_ALLOC(nvars, fq_nmod_mpoly_struct);
    gammaevals = FLINT_ARRAY_ALLOC(nvars, fq_nmod_mpoly_struct);
    for (i = 0; i < nvars; i++)
    {
        fq_nmod_mpoly_init3(Aevals + i, 0, bits, smctx);
        fq_nmod_mpoly_init3(Bevals + i, 0, bits, smctx);
        fq_nmod_mpoly_init3(gammaevals + i, 0, bits, smctx);
    }
    fq_nmod_mpoly_cvtto_mpolyn(An, A, nvars - 1, smctx);
    fq_nmod_mpoly_cvtto_mpolyn(Bn, B, nvars - 1, smctx);
    fq_nmod_mpoly_cvtto_mpolyn(gamman, gamma, nvars - 1, smctx);

    Abideg = _fq_nmod_mpoly_bidegree(A, smctx);
    Bbideg = _fq_nmod_mpoly_bidegree(B, smctx);

    degxAB = FLINT_MAX(Adegs[0], Bdegs[0]);
    degyAB = FLINT_MAX(Adegs[1], Bdegs[1]);

    GdegboundXY = pack_exp2(FLINT_MIN(Adegs[0], Bdegs[0]),
                            FLINT_MIN(Adegs[1], Bdegs[1]));
    if (GdegboundXY == 0)
        goto gcd_is_trivial;

    goto got_alpha_m;

increase_degree:

    do {
        cur_emb = bad_fq_nmod_mpoly_embed_chooser_next(embc, lgctx, smctx, state);
        if (cur_emb == NULL)
        {
            /* ran out of primes */
            success = 0;
            goto cleanup;
        }
    } while (fq_nmod_ctx_degree(lgctx->fqctx) <= lgd);

    lgd = fq_nmod_ctx_degree(lgctx->fqctx);

    goto got_alpha_m;

choose_main:

    cur_emb = bad_fq_nmod_mpoly_embed_chooser_next(embc, lgctx, smctx, state);
    if (cur_emb == NULL)
    {
        /* ran out of primes */
        success = 0;
        goto cleanup;
    }

    lgd = fq_nmod_ctx_degree(lgctx->fqctx);

got_alpha_m:

    n_poly_fit_length(tmp, lgd);
    n_poly_fit_length(alphapow, 2*lgd);

    for (i = 2; i < nvars - 1; i++)
        fq_nmod_rand_not_zero(alphas + i, state, lgctx->fqctx);

    i = nvars - 1;
    fq_nmod_mpolyn_interp_reduce_lg_mpoly(Aevals + i, An, lgctx, smctx, cur_emb);
    fq_nmod_mpolyn_interp_reduce_lg_mpoly(Bevals + i, Bn, lgctx, smctx, cur_emb);
    fq_nmod_mpolyn_interp_reduce_lg_mpoly(gammaevals + i, gamman, lgctx, smctx, cur_emb);
    if (fq_nmod_mpoly_is_zero(gammaevals + i, lgctx))
        goto choose_main;
    if (Aevals[i].length < 1 || _fq_nmod_mpoly_bidegree(Aevals + i, lgctx) != Abideg)
        goto choose_main;
    if (Bevals[i].length < 1 || _fq_nmod_mpoly_bidegree(Bevals + i, lgctx) != Bbideg)
        goto choose_main;
    for (i--; i >= 2; i--)
    {
        fq_nmod_mpoly_evaluate_one_fq_nmod(Aevals + i, Aevals + i + 1, i, alphas + i, lgctx);
        fq_nmod_mpoly_repack_bits_inplace(Aevals + i, bits, lgctx);
        fq_nmod_mpoly_evaluate_one_fq_nmod(Bevals + i, Bevals + i + 1, i, alphas + i, lgctx);
        fq_nmod_mpoly_repack_bits_inplace(Bevals + i, bits, lgctx);
        fq_nmod_mpoly_evaluate_one_fq_nmod(gammaevals + i, gammaevals + i + 1, i, alphas + i, lgctx);
        fq_nmod_mpoly_repack_bits_inplace(gammaevals + i, bits, lgctx);
        if (fq_nmod_mpoly_is_zero(gammaevals + i, lgctx))
            goto choose_main;
        if (Aevals[i].length < 1 || _fq_nmod_mpoly_bidegree(Aevals + i, lgctx) != Abideg)
            goto choose_main;
        if (Bevals[i].length < 1 || _fq_nmod_mpoly_bidegree(Bevals + i, lgctx) != Bbideg)
            goto choose_main;
    }

    m = 2;

    if (rGdegs == NULL)
        use = USE_G | USE_ABAR | USE_BBAR;
    else
        use = mpoly_gcd_get_use_first(rGdegs[m], Adegs[m], Bdegs[m], gammadegs[m]);

    fq_nmod_mpoly_get_n_fq_bpoly(Aev, Aevals + m, 0, 1, lgctx);
    fq_nmod_mpoly_get_n_fq_bpoly(Bev, Bevals + m, 0, 1, lgctx);

    success = n_fq_bpoly_gcd_brown_smprime(Gev, Abarev, Bbarev,
                                                   Aev, Bev, lgctx->fqctx, St);
    if (!success)
        goto increase_degree;

    newdegXY = n_bpoly_bidegree(Gev);
    if (newdegXY > GdegboundXY)
        goto choose_main;
    if (newdegXY < GdegboundXY)
    {
        GdegboundXY = newdegXY;
        if (GdegboundXY == 0)
            goto gcd_is_trivial;
    }

    gammaev = fq_nmod_mpoly_get_nonzero_n_fq(gammaevals + m, lgctx);
    n_fq_bpoly_scalar_mul_n_fq(Gev, gammaev, lgctx->fqctx);

    if (nvars == 3)
    {
        fq_nmod_mpolyn_interp_lift_lg_bpoly(&lastdeg, Gn, smctx, Gev, lgctx, cur_emb);
        fq_nmod_mpolyn_interp_lift_lg_bpoly(&lastdeg, Abarn, smctx, Abarev, lgctx, cur_emb);
        fq_nmod_mpolyn_interp_lift_lg_bpoly(&lastdeg, Bbarn, smctx, Bbarev, lgctx, cur_emb);

        n_fq_poly_set_fq_nmod_poly(modulus, cur_emb->h, smctx->fqctx);

        while (1)
        {
        choose_alpha_2_last:

            cur_emb = bad_fq_nmod_mpoly_embed_chooser_next(embc, lgctx, smctx, state);
            if (cur_emb == NULL)
            {
                /* ran out of primes */
                success = 0;
                goto cleanup;
            }

            lgd = fq_nmod_ctx_degree(lgctx->fqctx);
            n_poly_fit_length(tmp, lgd);
            n_poly_fit_length(alphapow, 2*lgd);
/*
flint_printf("choosing alpha[%wd]:\n", m);
fq_nmod_ctx_print(lgctx->fqctx);
flint_printf("\n");
*/

            fq_nmod_mpolyn_interp_reduce_lg_mpoly(Aevals + m, An, lgctx, smctx, cur_emb);
            fq_nmod_mpolyn_interp_reduce_lg_mpoly(Bevals + m, Bn, lgctx, smctx, cur_emb);
            fq_nmod_mpolyn_interp_reduce_lg_mpoly(gammaevals + m, gamman, lgctx, smctx, cur_emb);
            if (fq_nmod_mpoly_is_zero(gammaevals + m, lgctx))
                goto choose_alpha_2_last;
            if (Aevals[m].length < 1 || _fq_nmod_mpoly_bidegree(Aevals + m, lgctx) != Abideg)
                goto choose_alpha_2_last;
            if (Bevals[m].length < 1 || _fq_nmod_mpoly_bidegree(Bevals + m, lgctx) != Bbideg)
                goto choose_alpha_2_last;

            fq_nmod_mpoly_get_n_fq_bpoly(Aev, Aevals + m, 0, 1, lgctx);
            fq_nmod_mpoly_get_n_fq_bpoly(Bev, Bevals + m, 0, 1, lgctx);

            success = n_fq_bpoly_gcd_brown_smprime(Gev, Abarev, Bbarev,
                                                   Aev, Bev, lgctx->fqctx, St);
            if (!success)
                goto increase_degree;

            newdegXY = n_bpoly_bidegree(Gev);
            if (newdegXY > GdegboundXY)
                goto choose_alpha_2_last;
            if (newdegXY < GdegboundXY)
            {
                GdegboundXY = newdegXY;
                if (GdegboundXY == 0)
                    goto gcd_is_trivial;
                goto choose_main;
            }

            gammaev = fq_nmod_mpoly_get_nonzero_n_fq(gammaevals + m, lgctx);
            n_fq_bpoly_scalar_mul_n_fq(Gev, gammaev, lgctx->fqctx);

            if ((use & USE_G) && !fq_nmod_mpolyn_interp_crt_lg_bpoly(
                        &lastdeg, Gn, Tn, modulus, smctx, Gev, lgctx, cur_emb))
            {
                fq_nmod_mpoly_cvtfrom_mpolyn(rG, Gn, m, smctx);
                success = fq_nmod_mpolyl_content(cont, rG, 2, smctx);
                if (!success)
                    goto cleanup;
                fq_nmod_mpoly_divides(rG, rG, cont, smctx);
                fq_nmod_mpoly_repack_bits_inplace(rG, bits, smctx);
                if (fq_nmod_mpoly_divides(rAbar, A, rG, smctx) &&
                    fq_nmod_mpoly_divides(rBbar, B, rG, smctx))
                {
                    fq_nmod_mpoly_repack_bits_inplace(rAbar, bits, smctx);
                    fq_nmod_mpoly_repack_bits_inplace(rBbar, bits, smctx);
                    break;
                }
            }

            if ((use & USE_ABAR) && !fq_nmod_mpolyn_interp_crt_lg_bpoly(
                  &lastdeg, Abarn, Tn, modulus, smctx, Abarev, lgctx, cur_emb))
            {
                fq_nmod_mpoly_cvtfrom_mpolyn(rAbar, Abarn, m, smctx);
                success = fq_nmod_mpolyl_content(cont, rAbar, 2, smctx);
                if (!success)
                    goto cleanup;
                fq_nmod_mpoly_divides(rAbar, rAbar, cont, smctx);
                fq_nmod_mpoly_repack_bits_inplace(rAbar, bits, smctx);
                if (fq_nmod_mpoly_divides(rG, A, rAbar, smctx) &&
                    fq_nmod_mpoly_divides(rBbar, B, rG, smctx))
                {
                    fq_nmod_mpoly_repack_bits_inplace(rG, bits, smctx);
                    fq_nmod_mpoly_repack_bits_inplace(rBbar, bits, smctx);
                    break;
                }
            }

            if ((use & USE_BBAR) && !fq_nmod_mpolyn_interp_crt_lg_bpoly(
                  &lastdeg, Bbarn, Tn, modulus, smctx, Bbarev, lgctx, cur_emb))
            {
                fq_nmod_mpoly_cvtfrom_mpolyn(rBbar, Bbarn, m, smctx);
                success = fq_nmod_mpolyl_content(cont, rBbar, 2, smctx);
                if (!success)
                    goto cleanup;
                fq_nmod_mpoly_divides(rBbar, rBbar, cont, smctx);
                fq_nmod_mpoly_repack_bits_inplace(rBbar, bits, smctx);
                if (fq_nmod_mpoly_divides(rG, B, rBbar, smctx) &&
                    fq_nmod_mpoly_divides(rAbar, A, rG, smctx))
                {
                    fq_nmod_mpoly_repack_bits_inplace(rG, bits, smctx);
                    fq_nmod_mpoly_repack_bits_inplace(rAbar, bits, smctx);
                    break;
                }
            }

            if (n_fq_poly_degree(modulus) > gammadegs[m] + Adegs[m] &&
                n_fq_poly_degree(modulus) > gammadegs[m] + Bdegs[m])
            {
                goto choose_main;
            }

            {
            n_fq_poly_t h_;
            n_fq_poly_init(h_);
            n_fq_poly_set_fq_nmod_poly(h_, cur_emb->h, smctx->fqctx);
            n_fq_poly_mul(modulus, modulus, h_, smctx->fqctx);
            n_fq_poly_clear(h_);
            }
        }

        success = 1;
        goto cleanup;
    }

    fq_nmod_mpolyn_interp_lift_sm_bpoly(Gn, Gev, lgctx);
    fq_nmod_mpolyn_interp_lift_sm_bpoly(Abarn, Abarev, lgctx);
    fq_nmod_mpolyn_interp_lift_sm_bpoly(Bbarn, Bbarev, lgctx);

    FLINT_ASSERT(tmp->alloc >= lgd);
    n_fq_poly_one(modulus, lgctx->fqctx);
    n_fq_set_fq_nmod(tmp->coeffs, alphas + m, lgctx->fqctx);
    n_fq_poly_shift_left_scalar_submul(modulus, 1, tmp->coeffs, lgctx->fqctx);

    fq_nmod_set(start_alpha, alphas + m, lgctx->fqctx);

    while (1)
    {
    choose_alpha_2:

        fq_nmod_next_not_zero(alphas + m, lgctx->fqctx);
/*
flint_printf("choosing alpha[%wd]: ", m);
fq_nmod_print_pretty(alphas + m, lgctx->fqctx);
flint_printf("\n");
*/
        if (fq_nmod_equal(alphas + m, start_alpha, lgctx->fqctx))
            goto increase_degree;

        FLINT_ASSERT(alphapow->alloc >= lgd*2);
        alphapow->length = 2;
        _n_fq_one(alphapow->coeffs + lgd*0, lgd);
        n_fq_set_fq_nmod(alphapow->coeffs + lgd*1, alphas + m, lgctx->fqctx);

        fq_nmod_mpoly_evaluate_one_fq_nmod(Aevals + m, Aevals + m + 1, m, alphas + m, lgctx);
        fq_nmod_mpoly_repack_bits_inplace(Aevals + m, bits, lgctx);
        fq_nmod_mpoly_evaluate_one_fq_nmod(Bevals + m, Bevals + m + 1, m, alphas + m, lgctx);
        fq_nmod_mpoly_repack_bits_inplace(Bevals + m, bits, lgctx);
        fq_nmod_mpoly_evaluate_one_fq_nmod(gammaevals + m, gammaevals + m + 1, m, alphas + m, lgctx);
        fq_nmod_mpoly_repack_bits_inplace(gammaevals + m, bits, lgctx);
        if (fq_nmod_mpoly_is_zero(gammaevals + m, lgctx))
            goto choose_alpha_2;
        if (Aevals[m].length < 1 || _fq_nmod_mpoly_bidegree(Aevals + m, lgctx) != Abideg)
            goto choose_alpha_2;
        if (Bevals[m].length < 1 || _fq_nmod_mpoly_bidegree(Bevals + m, lgctx) != Bbideg)
            goto choose_alpha_2;

        fq_nmod_mpoly_get_n_fq_bpoly(Aev, Aevals + m, 0, 1, lgctx);
        fq_nmod_mpoly_get_n_fq_bpoly(Bev, Bevals + m, 0, 1, lgctx);

        success = n_fq_bpoly_gcd_brown_smprime(Gev, Abarev, Bbarev,
                                                   Aev, Bev, lgctx->fqctx, St);
        if (!success)
            goto increase_degree;

        newdegXY = n_bpoly_bidegree(Gev);
        if (newdegXY > GdegboundXY)
            goto choose_alpha_2;
        if (newdegXY < GdegboundXY)
        {
            GdegboundXY = newdegXY;
            if (GdegboundXY == 0)
                goto gcd_is_trivial;
            goto choose_main;
        }

        gammaev = fq_nmod_mpoly_get_nonzero_n_fq(gammaevals + m, lgctx);
        n_fq_bpoly_scalar_mul_n_fq(Gev, gammaev, lgctx->fqctx);

        FLINT_ASSERT(tmp->alloc >= lgd);
        n_fq_poly_eval_pow(tmp->coeffs, modulus, alphapow, lgctx->fqctx);
        n_fq_inv(tmp->coeffs, tmp->coeffs, lgctx->fqctx);
        n_fq_poly_scalar_mul_n_fq(modulus, modulus, tmp->coeffs, lgctx->fqctx);

        if ((use & USE_G) && !fq_nmod_mpolyn_interp_crt_sm_bpoly(
                              &lastdeg, Gn, Tn, Gev, modulus, alphapow, lgctx))
        {
            fq_nmod_mpoly_cvtfrom_mpolyn(G, Gn, m, lgctx);
            fq_nmod_mpoly_mul(T, Aevals + m + 1, gammaevals + m + 1, lgctx);
            if (fq_nmod_mpoly_divides(Abar, T, G, lgctx))
            {
                fq_nmod_mpoly_mul(T, Bevals + m + 1, gammaevals + m + 1, lgctx);
                if (fq_nmod_mpoly_divides(Bbar, T, G, lgctx))
                {
                    fq_nmod_mpoly_repack_bits_inplace(Abar, bits, lgctx);
                    fq_nmod_mpoly_repack_bits_inplace(Bbar, bits, lgctx);
                    break;
                }
            }
        }

        if ((use & USE_ABAR) && !fq_nmod_mpolyn_interp_crt_sm_bpoly(
                        &lastdeg, Abarn, Tn, Abarev, modulus, alphapow, lgctx))
        {
            fq_nmod_mpoly_cvtfrom_mpolyn(Abar, Abarn, m, lgctx);
            fq_nmod_mpoly_mul(T, Aevals + m + 1, gammaevals + m + 1, lgctx);
            if (fq_nmod_mpoly_divides(G, T, Abar, lgctx))
            {
                fq_nmod_mpoly_mul(T, Bevals + m + 1, gammaevals + m + 1, lgctx);
                if (fq_nmod_mpoly_divides(Bbar, T, G, lgctx))
                {
                    fq_nmod_mpoly_repack_bits_inplace(G, bits, lgctx);
                    fq_nmod_mpoly_repack_bits_inplace(Bbar, bits, lgctx);
                    break;
                }
            }
        }

        if ((use & USE_BBAR) && !fq_nmod_mpolyn_interp_crt_sm_bpoly(
                        &lastdeg, Bbarn, Tn, Bbarev, modulus, alphapow, lgctx))
        {
            fq_nmod_mpoly_cvtfrom_mpolyn(Bbar, Bbarn, m, lgctx);
            fq_nmod_mpoly_mul(T, Bevals + m + 1, gammaevals + m + 1, lgctx);
            if (fq_nmod_mpoly_divides(G, T, Bbar, lgctx))
            {
                fq_nmod_mpoly_mul(T, Aevals + m + 1, gammaevals + m + 1, lgctx);
                if (fq_nmod_mpoly_divides(Abar, T, G, lgctx))
                {
                    fq_nmod_mpoly_repack_bits_inplace(G, bits, lgctx);
                    fq_nmod_mpoly_repack_bits_inplace(Abar, bits, lgctx);
                    break;
                }
            }
        }

        if (n_fq_poly_degree(modulus) > gammadegs[m] + Adegs[m] &&
            n_fq_poly_degree(modulus) > gammadegs[m] + Bdegs[m])
        {
            goto choose_main;
        }

        FLINT_ASSERT(n_fq_equal_fq_nmod(alphapow->coeffs + lgd*1, alphas + m, lgctx->fqctx));
        n_fq_poly_shift_left_scalar_submul(modulus, 1, alphapow->coeffs + lgd*1, lgctx->fqctx);
    }

    for (m = 3; m < nvars - 1; m++)
    {
        /* G, Abar, Bbar are in Fq[gen(0), ..., gen(m - 1)] */
        fq_nmod_mpolyn_interp_lift_sm_mpoly(Gn, G, lgctx);
        fq_nmod_mpolyn_interp_lift_sm_mpoly(Abarn, Abar, lgctx);
        fq_nmod_mpolyn_interp_lift_sm_mpoly(Bbarn, Bbar, lgctx);

        if (rGdegs == NULL)
            use = USE_G | USE_ABAR | USE_BBAR;
        else
            use = 0;

        FLINT_ASSERT(tmp->alloc >= lgd);
        n_fq_poly_one(modulus, lgctx->fqctx);
        n_fq_set_fq_nmod(tmp->coeffs, alphas + m, lgctx->fqctx);
        n_fq_poly_shift_left_scalar_submul(modulus, 1, tmp->coeffs, lgctx->fqctx);

        fq_nmod_set(start_alpha, alphas + m, lgctx->fqctx);

    choose_betas_m:

        /* only beta[2], beta[1], ..., beta[m - 1] will be used */
        for (i = 2; i < nvars; i++)
            fq_nmod_rand_not_zero(betas + i, state, lgctx->fqctx);

        FLINT_ASSERT(G->bits == bits);
        FLINT_ASSERT(Abar->bits == bits);
        FLINT_ASSERT(Bbar->bits == bits);

        fq_nmod_mpoly_monomial_evals2(HG, G, betas + 2, m, lgctx);
        fq_nmod_mpoly_monomial_evals2(HAbar, Abar, betas + 2, m, lgctx);
        fq_nmod_mpoly_monomial_evals2(HBbar, Bbar, betas + 2, m, lgctx);
        if (use == 0)
        {
            this_length = gammaevals[m + 1].length + Aevals[m + 1].length +
                                                     Bevals[m + 1].length;

            use = nmod_mpoly_gcd_get_use_new(rGdegs[m], Adegs[m], Bdegs[m],
                   gammadegs[m], degxAB, degyAB, this_length, HG, HAbar, HBbar);
        }

        req_zip_images = 1;
        if (use & USE_G)
        {
            this_length = n_fq_polyun_product_roots(MG, HG, lgctx->fqctx);
            req_zip_images = FLINT_MAX(req_zip_images, this_length);
        }
        if (use & USE_ABAR)
        {
            this_length = n_fq_polyun_product_roots(MAbar, HAbar, lgctx->fqctx);
            req_zip_images = FLINT_MAX(req_zip_images, this_length);
        }
        if (use & USE_BBAR)
        {
            this_length = n_fq_polyun_product_roots(MBbar, HBbar, lgctx->fqctx);
            req_zip_images = FLINT_MAX(req_zip_images, this_length);
        }

        while (1)
        {
        choose_alpha_m:

            fq_nmod_next_not_zero(alphas + m, lgctx->fqctx);
/*
flint_printf("choosing alpha[%wd]: ", m);
fq_nmod_print_pretty(alphas + m, lgctx->fqctx);
flint_printf("\n");
*/

            if (fq_nmod_equal(alphas + m, start_alpha, lgctx->fqctx))
                goto choose_main;

            FLINT_ASSERT(alphapow->alloc >= lgd*2);
            alphapow->length = 2;
            _n_fq_one(alphapow->coeffs + lgd*0, lgd);
            n_fq_set_fq_nmod(alphapow->coeffs + lgd*1, alphas + m, lgctx->fqctx);

            fq_nmod_mpoly_evaluate_one_fq_nmod(Aevals + m, Aevals + m + 1, m, alphas + m, lgctx);
            fq_nmod_mpoly_repack_bits_inplace(Aevals + m, bits, lgctx);
            fq_nmod_mpoly_evaluate_one_fq_nmod(Bevals + m, Bevals + m + 1, m, alphas + m, lgctx);
            fq_nmod_mpoly_repack_bits_inplace(Bevals + m, bits, lgctx);
            fq_nmod_mpoly_evaluate_one_fq_nmod(gammaevals + m, gammaevals + m + 1, m, alphas + m, lgctx);
            fq_nmod_mpoly_repack_bits_inplace(gammaevals + m, bits, lgctx);
            if (fq_nmod_mpoly_is_zero(gammaevals + m, lgctx))
                goto choose_alpha_m;
            if (Aevals[m].length < 1 || _fq_nmod_mpoly_bidegree(Aevals + m, lgctx) != Abideg)
                goto choose_alpha_m;
            if (Bevals[m].length < 1 || _fq_nmod_mpoly_bidegree(Bevals + m, lgctx) != Bbideg)
                goto choose_alpha_m;

            fq_nmod_mpoly_monomial_evals2(Aeh_inc, Aevals + m, betas + 2, m, lgctx);
            fq_nmod_mpoly_monomial_evals2(Beh_inc, Bevals + m, betas + 2, m, lgctx);
            fq_nmod_mpoly_monomial_evals(gammaeh_inc, gammaevals + m, betas + 2, 2, m, lgctx);
            n_fq_polyun_set(Aeh_cur, Aeh_inc, lgctx->fqctx);
            n_fq_polyun_set(Beh_cur, Beh_inc, lgctx->fqctx);
            n_fq_poly_set(gammaeh_cur, gammaeh_inc, lgctx->fqctx);

            qzip_start(ZG, HG, req_zip_images, lgctx->fqctx);
            qzip_start(ZAbar, HAbar, req_zip_images, lgctx->fqctx);
            qzip_start(ZBbar, HBbar, req_zip_images, lgctx->fqctx);
            for (cur_zip_image = 0; cur_zip_image < req_zip_images; cur_zip_image++)
            {
                n_fq_bpoly_eval_step_new(Aev, Aeh_cur, Aeh_inc, Aevals + m, lgctx->fqctx);
                n_fq_bpoly_eval_step_new(Bev, Beh_cur, Beh_inc, Bevals + m, lgctx->fqctx);
                FLINT_ASSERT(tmp->alloc >= lgd);
                n_fq_poly_eval_step_new(tmp->coeffs, gammaeh_cur, gammaeh_inc, gammaevals + m, lgctx->fqctx);
                if (_n_fq_is_zero(tmp->coeffs, lgd))
                    goto choose_betas_m;
                if (Aev->length < 1 || n_bpoly_bidegree(Aev) != Abideg)
                    goto choose_betas_m;
                if (Bev->length < 1 || n_bpoly_bidegree(Bev) != Bbideg)
                    goto choose_betas_m;

                success = n_fq_bpoly_gcd_brown_smprime(Gev, Abarev, Bbarev,
                                                   Aev, Bev, lgctx->fqctx, St);
                if (!success)
                    goto increase_degree;

                newdegXY = n_bpoly_bidegree(Gev);
                if (newdegXY > GdegboundXY)
                    goto choose_betas_m;
                if (newdegXY < GdegboundXY)
                {
                    GdegboundXY = newdegXY;
                    if (GdegboundXY == 0)
                        goto gcd_is_trivial;
                    goto choose_main;
                }

                n_fq_bpoly_scalar_mul_n_fq(Gev, tmp->coeffs, lgctx->fqctx);

                if ((use & USE_G) && !n_fq_polyu2_add_zip_must_match(ZG, Gev, cur_zip_image, lgctx->fqctx))
                    goto choose_main;
                if ((use & USE_ABAR) && !n_fq_polyu2_add_zip_must_match(ZAbar, Abarev, cur_zip_image, lgctx->fqctx))
                    goto choose_main;
                if ((use & USE_BBAR) && !n_fq_polyu2_add_zip_must_match(ZBbar, Bbarev, cur_zip_image, lgctx->fqctx))
                    goto choose_main;
            }

            if ((use & USE_G) && qzip_solvel(G, ZG, HG, MG, lgctx) < 1)
                goto choose_main;
            if ((use & USE_ABAR) && qzip_solvel(Abar, ZAbar, HAbar, MAbar, lgctx) < 1)
                goto choose_main;
            if ((use & USE_BBAR) && qzip_solvel(Bbar, ZBbar, HBbar, MBbar, lgctx) < 1)
                goto choose_main;

            FLINT_ASSERT(n_fq_poly_degree(modulus) > 0);

            FLINT_ASSERT(tmp->alloc >= lgd);
            n_fq_poly_eval_pow(tmp->coeffs, modulus, alphapow, lgctx->fqctx);
            n_fq_inv(tmp->coeffs, tmp->coeffs, lgctx->fqctx);
            n_fq_poly_scalar_mul_n_fq(modulus, modulus, tmp->coeffs, lgctx->fqctx);

            if ((use & USE_G) && !fq_nmod_mpolyn_interp_mcrt_sm_mpoly(
                                    &lastdeg, Gn, G, modulus, alphapow, lgctx))
            {
                fq_nmod_mpoly_cvtfrom_mpolyn(rG, Gn, m, lgctx);
                fq_nmod_mpoly_mul(T, Aevals + m + 1, gammaevals + m + 1, lgctx);
                if (fq_nmod_mpoly_divides(rAbar, T, rG, lgctx))
                {
                    fq_nmod_mpoly_mul(T, Bevals + m + 1, gammaevals + m + 1, lgctx);
                    if (fq_nmod_mpoly_divides(rBbar, T, rG, lgctx))
                    {
                        fq_nmod_mpoly_repack_bits_inplace(rAbar, bits, lgctx);
                        fq_nmod_mpoly_repack_bits_inplace(rBbar, bits, lgctx);
                        fq_nmod_mpoly_swap(G, rG, lgctx);
                        fq_nmod_mpoly_swap(Abar, rAbar, lgctx);
                        fq_nmod_mpoly_swap(Bbar, rBbar, lgctx);
                        break;
                    }
                }
            }

            if ((use & USE_ABAR) && !fq_nmod_mpolyn_interp_mcrt_sm_mpoly(
                              &lastdeg, Abarn, Abar, modulus, alphapow, lgctx))
            {
                fq_nmod_mpoly_cvtfrom_mpolyn(rAbar, Abarn, m, lgctx);
                fq_nmod_mpoly_mul(T, Aevals + m + 1, gammaevals + m + 1, lgctx);
                if (fq_nmod_mpoly_divides(rG, T, rAbar, lgctx))
                {
                    fq_nmod_mpoly_mul(T, Bevals + m + 1, gammaevals + m + 1, lgctx);
                    if (fq_nmod_mpoly_divides(rBbar, T, rG, lgctx))
                    {
                        fq_nmod_mpoly_repack_bits_inplace(rG, bits, lgctx);
                        fq_nmod_mpoly_repack_bits_inplace(rBbar, bits, lgctx);
                        fq_nmod_mpoly_swap(G, rG, lgctx);
                        fq_nmod_mpoly_swap(Abar, rAbar, lgctx);
                        fq_nmod_mpoly_swap(Bbar, rBbar, lgctx);
                        break;
                    }
                }
            }

            if ((use & USE_BBAR) && !fq_nmod_mpolyn_interp_mcrt_sm_mpoly(
                              &lastdeg, Bbarn, Bbar, modulus, alphapow, lgctx))
            {
                fq_nmod_mpoly_cvtfrom_mpolyn(rBbar, Bbarn, m, lgctx);
                fq_nmod_mpoly_mul(T, Bevals + m + 1, gammaevals + m + 1, lgctx);
                if (fq_nmod_mpoly_divides(rG, T, rBbar, lgctx))
                {
                    fq_nmod_mpoly_mul(T, Aevals + m + 1, gammaevals + m + 1, lgctx);
                    if (fq_nmod_mpoly_divides(rAbar, T, rG, lgctx))
                    {
                        fq_nmod_mpoly_repack_bits_inplace(rG, bits, lgctx);
                        fq_nmod_mpoly_repack_bits_inplace(rAbar, bits, lgctx);
                        fq_nmod_mpoly_swap(G, rG, lgctx);
                        fq_nmod_mpoly_swap(Abar, rAbar, lgctx);
                        fq_nmod_mpoly_swap(Bbar, rBbar, lgctx);
                        break;
                    }
                }
            }

            if (n_fq_poly_degree(modulus) > gammadegs[m] + Adegs[m] &&
                n_fq_poly_degree(modulus) > gammadegs[m] + Bdegs[m])
            {
                goto choose_main;
            }

            FLINT_ASSERT(n_fq_equal_fq_nmod(alphapow->coeffs + lgd*1, alphas + m, lgctx->fqctx));
            n_fq_poly_shift_left_scalar_submul(modulus, 1, alphapow->coeffs + lgd*1, lgctx->fqctx);
        }
    }

    m = nvars - 1;
    {
        /* G, Abar, Bbar are in Fq/alpha(gen(m-1))[gen(0), ..., gen(m - 1)] */

        fq_nmod_mpolyn_interp_lift_lg_mpoly(Gn, smctx, G, lgctx, cur_emb);
        fq_nmod_mpolyn_interp_lift_lg_mpoly(Abarn, smctx, Abar, lgctx, cur_emb);
        fq_nmod_mpolyn_interp_lift_lg_mpoly(Bbarn, smctx, Bbar, lgctx, cur_emb);

        if (rGdegs == NULL)
            use = USE_G | USE_ABAR | USE_BBAR;
        else
            use = 0;

        n_fq_poly_set_fq_nmod_poly(modulus, cur_emb->h, smctx->fqctx);

        while (1)
        {
        choose_alpha_last:

            cur_emb = bad_fq_nmod_mpoly_embed_chooser_next(embc, lgctx, smctx, state);
            if (cur_emb == NULL)
            {
                /* ran out of primes */
                success = 0;
                goto cleanup;
            }
/*
flint_printf("choosing alpha[%wd]:\n", m);
fq_nmod_ctx_print(lgctx->fqctx);
flint_printf("\n");
*/
            lgd = fq_nmod_ctx_degree(lgctx->fqctx);
            n_poly_fit_length(tmp, lgd);
            n_poly_fit_length(alphapow, 2*lgd);

            fq_nmod_mpolyn_interp_reduce_lg_mpoly(Aevals + m, An, lgctx, smctx, cur_emb);
            fq_nmod_mpolyn_interp_reduce_lg_mpoly(Bevals + m, Bn, lgctx, smctx, cur_emb);
            fq_nmod_mpolyn_interp_reduce_lg_mpoly(gammaevals + m, gamman, lgctx, smctx, cur_emb);
            if (fq_nmod_mpoly_is_zero(gammaevals + m, lgctx))
                goto choose_alpha_last;
            if (Aevals[m].length < 1 || _fq_nmod_mpoly_bidegree(Aevals + m, lgctx) != Abideg)
                goto choose_alpha_last;
            if (Bevals[m].length < 1 || _fq_nmod_mpoly_bidegree(Bevals + m, lgctx) != Bbideg)
                goto choose_alpha_last;

        choose_betas_last:

            /* only beta[2], ..., beta[m - 1] will be used */
            for (i = 2; i < nvars; i++)
                fq_nmod_rand_not_zero(betas + i, state, lgctx->fqctx);

            fq_nmod_mpoly_monomial_evals2(HG, G, betas + 2, m, lgctx);
            fq_nmod_mpoly_monomial_evals2(HAbar, Abar, betas + 2, m, lgctx);
            fq_nmod_mpoly_monomial_evals2(HBbar, Bbar, betas + 2, m, lgctx);

            if (use == 0)
            {
                this_length = gamma->length + A->length + B->length;
                use = nmod_mpoly_gcd_get_use_new(rGdegs[m], Adegs[m], Bdegs[m],
                       gammadegs[m], degxAB, degyAB, this_length, HG, HAbar, HBbar);
            }

            req_zip_images = 1;
            if (use & USE_G)
            {
                this_length = n_fq_polyun_product_roots(MG, HG, lgctx->fqctx);
                req_zip_images = FLINT_MAX(req_zip_images, this_length);
            }
            if (use & USE_ABAR)
            {
                this_length = n_fq_polyun_product_roots(MAbar, HAbar, lgctx->fqctx);
                req_zip_images = FLINT_MAX(req_zip_images, this_length);
            }
            if (use & USE_BBAR)
            {
                this_length = n_fq_polyun_product_roots(MBbar, HBbar, lgctx->fqctx);
                req_zip_images = FLINT_MAX(req_zip_images, this_length);
            }

            fq_nmod_mpoly_monomial_evals2(Aeh_inc, Aevals + m, betas + 2, m, lgctx);
            fq_nmod_mpoly_monomial_evals2(Beh_inc, Bevals + m, betas + 2, m, lgctx);
            fq_nmod_mpoly_monomial_evals(gammaeh_inc, gammaevals + m, betas + 2, 2, m, lgctx);
            n_fq_polyun_set(Aeh_cur, Aeh_inc, lgctx->fqctx);
            n_fq_polyun_set(Beh_cur, Beh_inc, lgctx->fqctx);
            n_fq_poly_set(gammaeh_cur, gammaeh_inc, lgctx->fqctx);

            qzip_start(ZG, HG, req_zip_images, lgctx->fqctx);
            qzip_start(ZAbar, HAbar, req_zip_images, lgctx->fqctx);
            qzip_start(ZBbar, HBbar, req_zip_images, lgctx->fqctx);
            for (cur_zip_image = 0; cur_zip_image < req_zip_images; cur_zip_image++)
            {
                n_fq_bpoly_eval_step_new(Aev, Aeh_cur, Aeh_inc, Aevals + m, lgctx->fqctx);
                n_fq_bpoly_eval_step_new(Bev, Beh_cur, Beh_inc, Bevals + m, lgctx->fqctx);
                FLINT_ASSERT(tmp->alloc >= lgd);
                n_fq_poly_eval_step_new(tmp->coeffs, gammaeh_cur, gammaeh_inc, gammaevals + m, lgctx->fqctx);
                if (_n_fq_is_zero(tmp->coeffs, lgd))
                    goto choose_betas_last;
                if (Aev->length < 1 || n_bpoly_bidegree(Aev) != Abideg)
                    goto choose_betas_last;
                if (Bev->length < 1 || n_bpoly_bidegree(Bev) != Bbideg)
                    goto choose_betas_last;

                success = n_fq_bpoly_gcd_brown_smprime(Gev, Abarev, Bbarev,
                                                   Aev, Bev, lgctx->fqctx, St);        
                if (!success)
                    goto cleanup;

                newdegXY = n_bpoly_bidegree(Gev);
                if (newdegXY > GdegboundXY)
                    goto choose_betas_last;
                if (newdegXY < GdegboundXY)
                {
                    GdegboundXY = newdegXY;
                    if (GdegboundXY == 0)
                        goto gcd_is_trivial;
                    goto choose_main;
                }

                n_fq_bpoly_scalar_mul_n_fq(Gev, tmp->coeffs, lgctx->fqctx);

                if ((use & USE_G) && !n_fq_polyu2_add_zip_must_match(ZG, Gev, cur_zip_image, lgctx->fqctx))
                    goto choose_main;
                if ((use & USE_ABAR) && !n_fq_polyu2_add_zip_must_match(ZAbar, Abarev, cur_zip_image, lgctx->fqctx))
                    goto choose_main;
                if ((use & USE_BBAR) && !n_fq_polyu2_add_zip_must_match(ZBbar, Bbarev, cur_zip_image, lgctx->fqctx))
                    goto choose_main;
            }
            if ((use & USE_G) && qzip_solvel(G, ZG, HG, MG, lgctx) < 1)
                goto choose_main;
            if ((use & USE_ABAR) && qzip_solvel(Abar, ZAbar, HAbar, MAbar, lgctx) < 1)
                goto choose_main;
            if ((use & USE_BBAR) && qzip_solvel(Bbar, ZBbar, HBbar, MBbar, lgctx) < 1)
                goto choose_main;

            FLINT_ASSERT(n_fq_poly_degree(modulus) > 0);

            FLINT_ASSERT(tmp->alloc >= lgd);
            bad_n_fq_embed_sm_to_lg(tmp->coeffs, modulus, cur_emb);
            n_fq_inv(tmp->coeffs, tmp->coeffs, lgctx->fqctx);

            if ((use & USE_G) && !newfq_nmod_mpolyn_interp_mcrt_lg_mpoly(
                                    &lastdeg, Gn, smctx, modulus, tmp->coeffs,
                                                            G, lgctx, cur_emb))
            {
                fq_nmod_mpoly_cvtfrom_mpolyn(rG, Gn, m, smctx);
                success = fq_nmod_mpolyl_content(cont, rG, 2, smctx);
                if (!success)
                    goto cleanup;
                fq_nmod_mpoly_divides(rG, rG, cont, smctx);
                fq_nmod_mpoly_repack_bits_inplace(rG, bits, smctx);
                if (fq_nmod_mpoly_divides(rAbar, A, rG, smctx) &&
                    fq_nmod_mpoly_divides(rBbar, B, rG, smctx))
                {
                    fq_nmod_mpoly_repack_bits_inplace(rAbar, bits, smctx);
                    fq_nmod_mpoly_repack_bits_inplace(rBbar, bits, smctx);
                    break;
                }
            }

            if ((use & USE_ABAR) && !newfq_nmod_mpolyn_interp_mcrt_lg_mpoly(
                                 &lastdeg, Abarn, smctx, modulus, tmp->coeffs,
                                                         Abar, lgctx, cur_emb))
            {
                fq_nmod_mpoly_cvtfrom_mpolyn(rAbar, Abarn, m, smctx);
                success = fq_nmod_mpolyl_content(cont, rAbar, 2, smctx);
                if (!success)
                    goto cleanup;
                fq_nmod_mpoly_divides(rAbar, rAbar, cont, smctx);
                fq_nmod_mpoly_repack_bits_inplace(rAbar, bits, smctx);
                if (fq_nmod_mpoly_divides(rG, A, rAbar, smctx) &&
                    fq_nmod_mpoly_divides(rBbar, B, rG, smctx))
                {
                    fq_nmod_mpoly_repack_bits_inplace(rG, bits, smctx);
                    fq_nmod_mpoly_repack_bits_inplace(rBbar, bits, smctx);
                    break;
                }
            }

            if ((use & USE_BBAR) && !newfq_nmod_mpolyn_interp_mcrt_lg_mpoly(
                                 &lastdeg, Bbarn, smctx, modulus, tmp->coeffs,
                                                         Bbar, lgctx, cur_emb))
            {
                fq_nmod_mpoly_cvtfrom_mpolyn(rBbar, Bbarn, m, smctx);
                success = fq_nmod_mpolyl_content(cont, rBbar, 2, smctx);
                if (!success)
                    goto cleanup;
                fq_nmod_mpoly_divides(rBbar, rBbar, cont, smctx);
                fq_nmod_mpoly_repack_bits_inplace(rBbar, bits, smctx);
                if (fq_nmod_mpoly_divides(rG, B, rBbar, smctx) &&
                    fq_nmod_mpoly_divides(rAbar, A, rG, smctx))
                {
                    fq_nmod_mpoly_repack_bits_inplace(rG, bits, smctx);
                    fq_nmod_mpoly_repack_bits_inplace(rAbar, bits, smctx);
                    break;
                }
            }

            if (n_fq_poly_degree(modulus) > gammadegs[m] + Adegs[m] &&
                n_fq_poly_degree(modulus) > gammadegs[m] + Bdegs[m])
            {
                goto choose_main;
            }

            {
            n_fq_poly_t h_;
            n_fq_poly_init(h_);
            n_fq_poly_set_fq_nmod_poly(h_, cur_emb->h, smctx->fqctx);
            n_fq_poly_mul(modulus, modulus, h_, smctx->fqctx);
            n_fq_poly_clear(h_);
            }
        }
    }

    success = 1;

cleanup:

    for (i = 0; i < nvars; i++)
    {
        fq_nmod_clear(betas + i, smctx->fqctx);
        fq_nmod_clear(alphas + i, smctx->fqctx);
    }
    flint_free(betas);
    flint_free(alphas);

    fq_nmod_mpolyn_clear(An, smctx);
    fq_nmod_mpolyn_clear(Bn, smctx);
    fq_nmod_mpolyn_clear(gamman, smctx);
    bad_fq_nmod_mpoly_embed_chooser_clear(embc, lgctx, smctx, state);

    n_fq_polyun_clear(HG);
    n_fq_polyun_clear(HAbar);
    n_fq_polyun_clear(HBbar);
    n_fq_polyun_clear(MG);
    n_fq_polyun_clear(MAbar);
    n_fq_polyun_clear(MBbar);
    n_fq_polyun_clear(ZG);
    n_fq_polyun_clear(ZAbar);
    n_fq_polyun_clear(ZBbar);
    n_fq_bpoly_clear(Aev);
    n_fq_bpoly_clear(Bev);
    n_fq_bpoly_clear(Gev);
    n_fq_bpoly_clear(Abarev);
    n_fq_bpoly_clear(Bbarev);
    fq_nmod_mpoly_clear(cont, smctx);
    fq_nmod_mpoly_clear(T, smctx);
    fq_nmod_mpoly_clear(G, smctx);
    fq_nmod_mpoly_clear(Abar, smctx);
    fq_nmod_mpoly_clear(Bbar, smctx);
    fq_nmod_mpolyn_clear(Tn, smctx);
    fq_nmod_mpolyn_clear(Gn, smctx);
    fq_nmod_mpolyn_clear(Abarn, smctx);
    fq_nmod_mpolyn_clear(Bbarn, smctx);
    n_fq_polyun_clear(Aeh_cur);
    n_fq_polyun_clear(Aeh_inc);
    n_fq_polyun_clear(Beh_cur);
    n_fq_polyun_clear(Beh_inc);
    n_fq_poly_clear(gammaeh_cur);
    n_fq_poly_clear(gammaeh_inc);
    n_fq_poly_clear(modulus);
    n_poly_stack_clear(St->poly_stack);
    n_bpoly_stack_clear(St->bpoly_stack);

    n_poly_clear(tmp);
    n_fq_poly_clear(alphapow);
    fq_nmod_clear(start_alpha, smctx->fqctx);

    flint_randclear(state);

    for (i = 0; i < nvars; i++)
    {
        fq_nmod_mpoly_clear(Aevals + i, smctx);
        fq_nmod_mpoly_clear(Bevals + i, smctx);
        fq_nmod_mpoly_clear(gammaevals + i, smctx);
    }
    flint_free(Aevals);
    flint_free(Bevals);
    flint_free(gammaevals);

flint_printf("fq_nmod_mpolyl_gcd_zippel_lgprime returning %d\n", success);
FLINT_ASSERT(success);


    return success;

gcd_is_trivial:

    fq_nmod_mpoly_one(rG, smctx);
    fq_nmod_mpoly_set(rAbar, A, smctx);
    fq_nmod_mpoly_set(rBbar, B, smctx);

    success = 1;
    
    goto cleanup;
}


/* should find its way back here in interesting cases */
int fq_nmod_mpoly_gcd_zippel2(
    fq_nmod_mpoly_t G,
    const fq_nmod_mpoly_t A,
    const fq_nmod_mpoly_t B,
    const fq_nmod_mpoly_ctx_t ctx)
{
    if (fq_nmod_mpoly_is_zero(A, ctx) || fq_nmod_mpoly_is_zero(B, ctx))
        return fq_nmod_mpoly_gcd(G, A, B, ctx);

    return _fq_nmod_mpoly_gcd_algo(G, NULL, NULL, A, B, ctx, MPOLY_GCD_USE_ZIPPEL2);
}

