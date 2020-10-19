/*
    Copyright (C) 2019 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fq_nmod_mpoly.h"
#include "fq_nmod_mpoly_factor.h"
#include "profiler.h"

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


static void fq_nmod_mpolyu_set_evalp_helper(
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

static void fq_nmod_mpolyuun_interp_lift_lg_bpoly(
    slong * lastdeg_,
    fq_nmod_mpolyun_t A,
    const fq_nmod_mpoly_ctx_t smctx,
    n_fq_bpoly_t B,
    const fq_nmod_mpoly_ctx_t lgctx,
    const bad_fq_nmod_embed_t emb)
{
    slong lgd = fq_nmod_ctx_degree(lgctx->fqctx);
    slong N = mpoly_words_per_exp_sp(A->bits, smctx->minfo);
    slong i, j, Ai;
    slong lastdeg = -WORD(1);

    Ai = 0;
    for (i = B->length - 1; i >= 0; i--)
    {
        n_poly_struct * Bi = B->coeffs + i;
        for (j = Bi->length - 1; j >= 0; j--)
        {
            if (_n_fq_is_zero(Bi->coeffs + lgd*j, lgd))
                continue;

            fq_nmod_mpolyun_fit_length(A, Ai + 1, smctx);
            A->exps[Ai] = pack_exp2(i, j);
            fq_nmod_mpolyn_fit_length(A->coeffs + Ai, 1, smctx);
            A->coeffs[Ai].length = 1;
            mpoly_monomial_zero(A->coeffs[Ai].exps + N*0, N);
            bad_n_fq_embed_lg_to_sm(A->coeffs[Ai].coeffs + 0, Bi->coeffs + lgd*j, emb);
            lastdeg = FLINT_MAX(lastdeg, n_fq_poly_degree(A->coeffs[Ai].coeffs + 0));
            Ai++;
        }
    }
    A->length = Ai;

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

int fq_nmod_mpolyuun_interp_crt_lg_bpoly(
    slong * lastdeg,
    fq_nmod_mpolyun_t F,
    fq_nmod_mpolyun_t T,
    n_fq_poly_t modulus,
    const fq_nmod_mpoly_ctx_t smctx,
    n_fq_bpoly_t A,
    const fq_nmod_mpoly_ctx_t lgctx,
    const bad_fq_nmod_embed_t emb)
{
    slong lgd = fq_nmod_ctx_degree(lgctx->fqctx);
    int changed = 0;
    slong N = mpoly_words_per_exp(T->bits, smctx->minfo);
    n_fq_poly_struct * Acoeffs = A->coeffs;
    slong Fi, Ti, Ai, ai;
    slong Flen = F->length;
    ulong * Fexps = F->exps;
    fq_nmod_mpolyn_struct * Fcoeffs = F->coeffs;
    ulong * Texps = T->exps;
    fq_nmod_mpolyn_struct * Tcoeffs = T->coeffs;
    mp_limb_t * u = FLINT_ARRAY_ALLOC(3*lgd, mp_limb_t);
    mp_limb_t * v = u + lgd;
    mp_limb_t * inv_m_eval = v + lgd;
    n_fq_poly_t u_sm;

    n_fq_poly_init(u_sm);

    bad_n_fq_embed_sm_to_lg(u, modulus, emb);
    n_fq_inv(inv_m_eval, u, lgctx->fqctx);

    FLINT_ASSERT(fq_nmod_mpolyun_is_canonical(F, smctx));
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
            slong extra = Flen - Fi;
            extra = FLINT_MAX(extra, Ai);
            fq_nmod_mpolyun_fit_length(T, Ti + extra + 1, smctx);
            Tcoeffs = T->coeffs;
            Texps = T->exps;
        }

        fq_nmod_mpolyn_fit_length(Tcoeffs + Ti, 1, smctx);
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

            bad_n_fq_embed_sm_to_lg(u, getcf(Fcoeffs + Fi), emb);
            n_fq_sub(v, Acoeffs[Ai].coeffs + lgd*ai, u, lgctx->fqctx);
            if (!_n_fq_is_zero(v, lgd))
            {
                changed = 1;
                n_fq_mul(u, v, inv_m_eval, lgctx->fqctx);
                bad_n_fq_embed_lg_to_sm(u_sm, u, emb);
                n_fq_poly_mul(getcf(Tcoeffs + Ti), modulus, u_sm, smctx->fqctx);
                n_fq_poly_add(getcf(Tcoeffs + Ti), getcf(Tcoeffs + Ti), getcf(Fcoeffs + Fi), smctx->fqctx);
            }
            else
            {
                n_fq_poly_set(getcf(Tcoeffs + Ti), getcf(Fcoeffs + Fi), smctx->fqctx);
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
        else if (Fi < Flen && (Ai < 0 || Fexps[Fi] > pack_exp2(Ai, ai)))
        {
            /* F term ok, Aterm missing */

            Texps[Ti] = Fexps[Fi];

            bad_n_fq_embed_sm_to_lg(v, getcf(Fcoeffs + Fi), emb);
            if (!_n_fq_is_zero(v, lgd))
            {
                changed = 1;
                n_fq_mul(u, v, inv_m_eval, lgctx->fqctx);
                bad_n_fq_embed_lg_to_sm(u_sm, u, emb);
                n_fq_poly_mul(getcf(Tcoeffs + Ti), modulus, u_sm, smctx->fqctx);
                n_fq_poly_add(getcf(Tcoeffs + Ti), getcf(Tcoeffs + Ti), getcf(Fcoeffs + Fi), smctx->fqctx);
            }
            else
            {
                n_fq_poly_set(getcf(Tcoeffs + Ti), getcf(Fcoeffs + Fi), smctx->fqctx);
            }

            Fi++;
        }
        else
        {
            FLINT_ASSERT(Ai >= 0 && (Fi >= Flen || Fexps[Fi] < pack_exp2(Ai, ai)));

            /* F term missing, A term ok */

            Texps[Ti] = pack_exp2(Ai, ai);

            changed = 1;
            n_fq_mul(u, Acoeffs[Ai].coeffs + lgd*ai, inv_m_eval, lgctx->fqctx);
            bad_n_fq_embed_lg_to_sm(u_sm, u, emb);
            n_fq_poly_mul(getcf(Tcoeffs + Ti), modulus, u_sm, smctx->fqctx);

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

        FLINT_ASSERT(!n_fq_poly_is_zero(getcf(Tcoeffs + Ti)));
        *lastdeg = FLINT_MAX(*lastdeg, n_fq_poly_degree(getcf(Tcoeffs + Ti)));
        Ti++;
    }

    T->length = Ti;

    if (changed)
        fq_nmod_mpolyun_swap(T, F);

    FLINT_ASSERT(fq_nmod_mpolyun_is_canonical(F, smctx));

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


static slong fq_nmod_mpolyu_set_zip_formp(
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
int fq_nmod_mpolyuu_gcd_zippel_smprime(
    fq_nmod_mpolyu_t rG, const slong * rGdegs,
    fq_nmod_mpolyu_t rAbar,
    fq_nmod_mpolyu_t rBbar,
    const fq_nmod_mpolyu_t A, const slong * Adegs,
    const fq_nmod_mpolyu_t B, const slong * Bdegs,
    const fq_nmod_mpoly_t gamma, const slong * gammadegs,
    const fq_nmod_mpoly_ctx_t ctx)
{
    slong d = fq_nmod_ctx_degree(ctx->fqctx);
    int success, use, betas_in_fp, main_tries_left;
    slong i, j, m;
    slong nvars = ctx->minfo->nvars;
    flint_bitcnt_t bits = A->bits;
    fq_nmod_struct * alphas, * betas;
    flint_rand_t state;
    fq_nmod_mpoly_t cont;
    fq_nmod_mpolyu_t T, G, Abar, Bbar;
    n_fq_polyun_t HG, HAbar, HBbar, MG, MAbar, MBbar, ZG, ZAbar, ZBbar;
    n_fq_bpoly_t Aev, Bev, Gev, Abarev, Bbarev;
    const mp_limb_t * gammaev;
    fq_nmod_mpolyun_t Tn, Gn, Abarn, Bbarn;
    slong lastdeg, cur_zip_image, req_zip_images, this_images;
    n_fq_polyun_t Aeh, Beh;
    n_fq_poly_t gammaeh;
    fq_nmod_mpolyu_struct * Aevals, * Bevals;
    fq_nmod_mpoly_struct * gammaevals;
    n_fq_poly_t modulus, alphapow;
    n_poly_bpoly_stack_t St;
    mp_limb_t * c;
    fq_nmod_t start_alpha;
    ulong GdegboundXY, newdegXY;
    slong degxAB, degyAB;
timeit_t timer;

flint_printf("fq_nmod_mpolyuu_gcd_zippel_smprime called\n");

    FLINT_ASSERT(bits <= FLINT_BITS);
    FLINT_ASSERT(bits == rG->bits);
    FLINT_ASSERT(bits == rAbar->bits);
    FLINT_ASSERT(bits == rBbar->bits);
    FLINT_ASSERT(bits == A->bits);
    FLINT_ASSERT(bits == B->bits);
    FLINT_ASSERT(bits == gamma->bits);

#if FLINT_WANT_ASSERT
    {
        slong * tmp_degs = FLINT_ARRAY_ALLOC(nvars, slong);

        fq_nmod_mpolyu_degrees_si(tmp_degs, A, ctx);
        for (j = 0; j < nvars; j++)
            FLINT_ASSERT(tmp_degs[j] == Adegs[j]);

        fq_nmod_mpolyu_degrees_si(tmp_degs, B, ctx);
        for (j = 0; j < nvars; j++)
            FLINT_ASSERT(tmp_degs[j] == Bdegs[j]);

        fq_nmod_mpoly_degrees_si(tmp_degs, gamma, ctx);
        for (j = 0; j < nvars; j++)
            FLINT_ASSERT(tmp_degs[j] == gammadegs[j]);

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
    fq_nmod_mpolyu_init(T, bits, ctx);
    fq_nmod_mpolyu_init(G, bits, ctx);
    fq_nmod_mpolyu_init(Abar, bits, ctx);
    fq_nmod_mpolyu_init(Bbar, bits, ctx);
    fq_nmod_mpolyun_init(Tn, bits, ctx);
    fq_nmod_mpolyun_init(Gn, bits, ctx);
    fq_nmod_mpolyun_init(Abarn, bits, ctx);
    fq_nmod_mpolyun_init(Bbarn, bits, ctx);
    n_fq_polyun_init(Aeh);
    n_fq_polyun_init(Beh);
    n_fq_poly_init(gammaeh);
    n_fq_poly_init(modulus);
    n_poly_stack_init(St->poly_stack);
    n_bpoly_stack_init(St->bpoly_stack);

    betas = FLINT_ARRAY_ALLOC(nvars, fq_nmod_struct);
    alphas = FLINT_ARRAY_ALLOC(nvars, fq_nmod_struct);
    for (i = 0; i < nvars; i++)
    {
        fq_nmod_init(betas + i, ctx->fqctx);
        fq_nmod_init(alphas + i, ctx->fqctx);
    }

    flint_randinit(state);

    Aevals = FLINT_ARRAY_ALLOC(nvars + 1, fq_nmod_mpolyu_struct);
    Bevals = FLINT_ARRAY_ALLOC(nvars + 1, fq_nmod_mpolyu_struct);
    gammaevals = FLINT_ARRAY_ALLOC(nvars + 1, fq_nmod_mpoly_struct);
    for (i = 0; i < nvars; i++)
    {
        fq_nmod_mpolyu_init(Aevals + i, bits, ctx);
        fq_nmod_mpolyu_init(Bevals + i, bits, ctx);
        fq_nmod_mpoly_init3(gammaevals + i, 0, bits, ctx);
    }
    Aevals[nvars] = *A;
    Bevals[nvars] = *B;
    gammaevals[nvars] = *gamma;

    /* the products gamma*A and gamma*B should not overflow */
    for (j = 0; j < nvars; j++)
    {
        FLINT_ASSERT(FLINT_BIT_COUNT(gammadegs[j] + Adegs[j]) < bits);
        FLINT_ASSERT(FLINT_BIT_COUNT(gammadegs[j] + Bdegs[j]) < bits);
    }

    GdegboundXY = mpoly_bivar_degrees1_chained(0, A->exps, A->length);
    GdegboundXY = mpoly_bivar_degrees1_chained(GdegboundXY, B->exps, B->length);
    degxAB = extract_exp(GdegboundXY, 0, 2);
    degyAB = extract_exp(GdegboundXY, 1, 2);

    GdegboundXY = FLINT_MIN(A->exps[0], B->exps[0]);
    if (GdegboundXY == 0)
        goto gcd_is_trivial;

    main_tries_left = 10;

choose_main:

    if (--main_tries_left < 0)
    {
        success = 0;
        goto cleanup;
    }

    for (i = 0; i < ctx->minfo->nvars; i++)
    {
        fq_nmod_rand(alphas + i, state, ctx->fqctx);
        if (fq_nmod_is_zero(alphas + i, ctx->fqctx))
            fq_nmod_one(alphas + i, ctx->fqctx);
    }

    for (i = nvars - 1; i >= 0; i--)
    {
        fq_nmod_mpolyu_evaluate_one_fq_nmod(Aevals + i, Aevals + i + 1, i, alphas + i, ctx);
        fq_nmod_mpolyu_evaluate_one_fq_nmod(Bevals + i, Bevals + i + 1, i, alphas + i, ctx);
        fq_nmod_mpoly_evaluate_one_fq_nmod(gammaevals + i, gammaevals + i + 1, i, alphas + i, ctx);
        fq_nmod_mpoly_repack_bits_inplace(gammaevals + i, bits, ctx);
        if (fq_nmod_mpoly_is_zero(gammaevals + i, ctx))
            goto choose_main;
        if (Aevals[i].length < 1 || Aevals[i].exps[0] != A->exps[0])
            goto choose_main;
        if (Bevals[i].length < 1 || Bevals[i].exps[0] != B->exps[0])
            goto choose_main;
    }

    m = 0;

    if (rGdegs == NULL)
    {
        use = USE_G | USE_ABAR | USE_BBAR;
    }
    else
    {
        use = mpoly_gcd_get_use_first(rGdegs[m], Adegs[m], Bdegs[m], gammadegs[m]);
    }

    FLINT_ASSERT(alphapow->alloc >= d*2);
    alphapow->length = 2;
    _n_fq_one(alphapow->coeffs + d*0, d);
    n_fq_set_fq_nmod(alphapow->coeffs + d*1, alphas + m, ctx->fqctx);

    fq_nmod_mpolyuu_get_n_fq_bpoly(Aev, Aevals + m, ctx);
    fq_nmod_mpolyuu_get_n_fq_bpoly(Bev, Bevals + m, ctx);

    success = n_fq_bpoly_gcd_brown_smprime(Gev, Abarev, Bbarev, Aev, Bev, ctx->fqctx, St);
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

    fq_nmod_mpolyuun_interp_lift_sm_bpoly(Gn, Gev, ctx);
    fq_nmod_mpolyuun_interp_lift_sm_bpoly(Abarn, Abarev, ctx);
    fq_nmod_mpolyuun_interp_lift_sm_bpoly(Bbarn, Bbarev, ctx);

    n_fq_poly_one(modulus, ctx->fqctx);
    n_fq_set_fq_nmod(c, alphas + m, ctx->fqctx);
    n_fq_poly_shift_left_scalar_submul(modulus, 1, c, ctx->fqctx);

    fq_nmod_set(start_alpha, alphas + m, ctx->fqctx);

    while (1)
    {
choose_alpha_0:

        fq_nmod_next_not_zero(alphas + m, ctx->fqctx);
        if (fq_nmod_equal(alphas + m, start_alpha, ctx->fqctx))
            goto choose_main;

        FLINT_ASSERT(alphapow->alloc >= d*2);
        alphapow->length = 2;
        _n_fq_one(alphapow->coeffs + d*0, d);
        n_fq_set_fq_nmod(alphapow->coeffs + d*1, alphas + m, ctx->fqctx);

        fq_nmod_mpolyu_evaluate_one_fq_nmod(Aevals + m, Aevals + m + 1, m, alphas + m, ctx);
        fq_nmod_mpolyu_evaluate_one_fq_nmod(Bevals + m, Bevals + m + 1, m, alphas + m, ctx);
        fq_nmod_mpoly_evaluate_one_fq_nmod(gammaevals + m, gammaevals + m + 1, m, alphas + m, ctx);
        if (fq_nmod_mpoly_is_zero(gammaevals + m, ctx))
            goto choose_alpha_0;
        if (Aevals[m].length < 1 || Aevals[m].exps[0] != A->exps[0])
            goto choose_alpha_0;
        if (Bevals[m].length < 1 || Bevals[m].exps[0] != B->exps[0])
            goto choose_alpha_0;

        fq_nmod_mpolyuu_get_n_fq_bpoly(Aev, Aevals + m, ctx);
        fq_nmod_mpolyuu_get_n_fq_bpoly(Bev, Bevals + m, ctx);

        success = n_fq_bpoly_gcd_brown_smprime(Gev, Abarev, Bbarev, Aev, Bev, ctx->fqctx, St);
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
            goto choose_main;
        }

        gammaev = fq_nmod_mpoly_get_nonzero_n_fq(gammaevals + m, ctx);
        n_fq_bpoly_scalar_mul_n_fq(Gev, gammaev, ctx->fqctx);

        n_fq_poly_eval_pow(c, modulus, alphapow, ctx->fqctx);
        n_fq_inv(c, c, ctx->fqctx);
        n_fq_poly_scalar_mul_n_fq(modulus, modulus, c, ctx->fqctx);

        if ((use & USE_G) && !fq_nmod_mpolyuun_interp_crt_sm_bpoly(
                                &lastdeg, Gn, Tn, Gev, modulus, alphapow, ctx))
        {
timeit_start(timer);
            if (m == nvars - 1)
            {
                fq_nmod_mpolyu_cvtfrom_mpolyun(rG, Gn, m, ctx);
                success = fq_nmod_mpolyu_content_mpoly(cont, rG, ctx);
                if (!success)
                    goto cleanup;
                fq_nmod_mpolyu_divexact_mpoly_inplace(rG, cont, ctx);
                if (fq_nmod_mpolyuu_divides(rAbar, A, rG, 2, ctx) &&
                    fq_nmod_mpolyuu_divides(rBbar, B, rG, 2, ctx))
                {
timeit_stop(timer);
flint_printf("m = %wd divisibility: %wd\n", m, timer->wall);
                    break;
                }                
            }
            else
            {
                fq_nmod_mpolyu_cvtfrom_mpolyun(G, Gn, m, ctx);
                fq_nmod_mpolyu_mul_mpoly(T, Aevals + m + 1, gammaevals + m + 1, ctx);
                if (fq_nmod_mpolyuu_divides(Abar, T, G, 2, ctx))
                {
                    fq_nmod_mpolyu_mul_mpoly(T, Bevals + m + 1, gammaevals + m + 1, ctx);
                    if (fq_nmod_mpolyuu_divides(Bbar, T, G, 2, ctx))
                    {
timeit_stop(timer);
flint_printf("m = %wd divisibility: %wd\n", m, timer->wall);
                        break;
                    }
                }
            }
        }

        if ((use & USE_ABAR) && !fq_nmod_mpolyuun_interp_crt_sm_bpoly(
                          &lastdeg, Abarn, Tn, Abarev, modulus, alphapow, ctx))
        {
            if (m == nvars - 1)
            {
                fq_nmod_mpolyu_cvtfrom_mpolyun(rAbar, Abarn, m, ctx);
                success = fq_nmod_mpolyu_content_mpoly(cont, rAbar, ctx);
                if (!success)
                    goto cleanup;
                fq_nmod_mpolyu_divexact_mpoly_inplace(rAbar, cont, ctx);
                if (fq_nmod_mpolyuu_divides(rG, A, rAbar, 2, ctx) &&
                    fq_nmod_mpolyuu_divides(rBbar, B, rG, 2, ctx))
                {
                    break;
                }
            }
            else
            {
                fq_nmod_mpolyu_cvtfrom_mpolyun(Abar, Abarn, m, ctx);
                fq_nmod_mpolyu_mul_mpoly(T, Aevals + m + 1, gammaevals + m + 1, ctx);
                if (fq_nmod_mpolyuu_divides(G, T, Abar, 2, ctx))
                {
                    fq_nmod_mpolyu_mul_mpoly(T, Bevals + m + 1, gammaevals + m + 1, ctx);
                    if (fq_nmod_mpolyuu_divides(Bbar, T, G, 2, ctx))
                        break;
                }
            }
        }

        if ((use & USE_BBAR) && !fq_nmod_mpolyuun_interp_crt_sm_bpoly(
                          &lastdeg, Bbarn, Tn, Bbarev, modulus, alphapow, ctx))
        {
            if (m == nvars - 1)
            {
                fq_nmod_mpolyu_cvtfrom_mpolyun(rBbar, Bbarn, m, ctx);
                success = fq_nmod_mpolyu_content_mpoly(cont, rBbar, ctx);
                if (!success)
                    goto cleanup;
                fq_nmod_mpolyu_divexact_mpoly_inplace(rBbar, cont, ctx);
                if (fq_nmod_mpolyuu_divides(rG, A, rBbar, 2, ctx) &&
                    fq_nmod_mpolyuu_divides(rAbar, B, rG, 2, ctx))
                {
                    break;
                }
            }
            else
            {
                fq_nmod_mpolyu_cvtfrom_mpolyun(Bbar, Bbarn, m, ctx);
                fq_nmod_mpolyu_mul_mpoly(T, Bevals + m + 1, gammaevals + m + 1, ctx);
                if (fq_nmod_mpolyuu_divides(G, T, Bbar, 2, ctx))
                {
                    fq_nmod_mpolyu_mul_mpoly(T, Aevals + m + 1, gammaevals + m + 1, ctx);
                    if (fq_nmod_mpolyuu_divides(Abar, T, G, 2, ctx))
                        break;
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

    for (m = 1; m < nvars; m++)
    {
        /* G, Abar, Bbar are in Fq[gen(0), ..., gen(m - 1)] */
        fq_nmod_mpolyun_interp_lift_sm_mpolyu(Gn, G, ctx);
        fq_nmod_mpolyun_interp_lift_sm_mpolyu(Abarn, Abar, ctx);
        fq_nmod_mpolyun_interp_lift_sm_mpolyu(Bbarn, Bbar, ctx);

        if (rGdegs == NULL)
        {
            use = USE_G | USE_ABAR | USE_BBAR;
        }
        else
        {
            slong numABgamma = gamma[m + 1].length +
                               fq_nmod_mpolyu_total_length(Aevals + m + 1) +
                               fq_nmod_mpolyu_total_length(Bevals + m + 1);

            use = fq_nmod_mpoly_gcd_get_use(rGdegs[m], Adegs[m], Bdegs[m],
                      gammadegs[m], degxAB, degyAB, numABgamma, G, Abar, Bbar);
        }

        n_fq_poly_one(modulus, ctx->fqctx);
        n_fq_set_fq_nmod(c, alphas + m, ctx->fqctx);
        n_fq_poly_shift_left_scalar_submul(modulus, 1, c, ctx->fqctx);

        fq_nmod_set(start_alpha, alphas + m, ctx->fqctx);

    choose_betas:

        /* only beta[0], beta[1], ..., beta[m - 1] will be used */
        betas_in_fp = (ctx->fqctx->mod.norm < FLINT_BITS/4);
        if (betas_in_fp)
        {
            for (i = 0; i < ctx->minfo->nvars; i++)
            {
                fq_nmod_set_ui(betas + i,
                      n_urandint(state, ctx->fqctx->mod.n - 3) + 2, ctx->fqctx);
            }
        }
        else
        {
            for (i = 0; i < ctx->minfo->nvars; i++)
            {
                fq_nmod_rand(betas + i, state, ctx->fqctx);
                if (fq_nmod_is_zero(betas + i, ctx->fqctx))
                    fq_nmod_one(betas + i, ctx->fqctx);
            }
        }

        req_zip_images = 2;
        if (use & USE_G)
        {
            this_images = betas_in_fp ?
                        fq_nmod_mpolyu_set_zip_formp(HG, MG, G, betas, m, ctx) :
                        fq_nmod_mpolyu_set_zip_form(HG, MG, G, betas, m, ctx);
            req_zip_images = FLINT_MAX(req_zip_images, this_images);
        }
        if (use & USE_ABAR)
        {
            this_images = betas_in_fp ?
                fq_nmod_mpolyu_set_zip_formp(HAbar, MAbar, Abar, betas, m, ctx) :
                fq_nmod_mpolyu_set_zip_form(HAbar, MAbar, Abar, betas, m, ctx);
            req_zip_images = FLINT_MAX(req_zip_images, this_images);            
        }
        if (use & USE_BBAR)
        {
            this_images = betas_in_fp ?
                fq_nmod_mpolyu_set_zip_formp(HBbar, MBbar, Bbar, betas, m, ctx) :
                fq_nmod_mpolyu_set_zip_form(HBbar, MBbar, Bbar, betas, m, ctx);                
            req_zip_images = FLINT_MAX(req_zip_images, this_images);            
        }

        while (1)
        {
choose_alpha_m:

            fq_nmod_next_not_zero(alphas + m, ctx->fqctx);
            if (fq_nmod_equal(alphas + m, start_alpha, ctx->fqctx))
                goto choose_main;

            FLINT_ASSERT(alphapow->alloc >= d*2);
            alphapow->length = 2;
            _n_fq_one(alphapow->coeffs + d*0, d);
            n_fq_set_fq_nmod(alphapow->coeffs + d*1, alphas + m, ctx->fqctx);

            fq_nmod_mpolyu_evaluate_one_fq_nmod(Aevals + m, Aevals + m + 1, m, alphas + m, ctx);
            fq_nmod_mpolyu_evaluate_one_fq_nmod(Bevals + m, Bevals + m + 1, m, alphas + m, ctx);
            fq_nmod_mpoly_evaluate_one_fq_nmod(gammaevals + m, gammaevals + m + 1, m, alphas + m, ctx);
            if (fq_nmod_mpoly_is_zero(gammaevals + m, ctx))
                goto choose_alpha_m;
            if (Aevals[m].length < 1 || Aevals[m].exps[0] != A->exps[0])
                goto choose_alpha_m;
            if (Bevals[m].length < 1 || Bevals[m].exps[0] != B->exps[0])
                goto choose_alpha_m;

            qzip_start(ZG, HG, req_zip_images, ctx->fqctx);
            qzip_start(ZAbar, HAbar, req_zip_images, ctx->fqctx);
            qzip_start(ZBbar, HBbar, req_zip_images, ctx->fqctx);

            if (betas_in_fp)
            {
                fq_nmod_mpolyu_set_evalp_helper(Aeh, Aevals + m, betas, m, ctx);
                fq_nmod_mpolyu_set_evalp_helper(Beh, Bevals + m, betas, m, ctx);
                fq_nmod_mpoly_set_evalp_helper(gammaeh, gammaevals + m, betas, m, ctx);
            }
            else
            {
                fq_nmod_mpolyu_set_eval_helper(Aeh, Aevals + m, betas, m, ctx);
                fq_nmod_mpolyu_set_eval_helper(Beh, Bevals + m, betas, m, ctx);
                fq_nmod_mpoly_set_eval_helper(gammaeh, gammaevals + m, betas, m, ctx);
            }

            for (cur_zip_image = 0; cur_zip_image < req_zip_images; cur_zip_image++)
            {
                if (betas_in_fp)
                {
                    n_fq_bpoly_evalp_step(Aev, Aeh, ctx->fqctx);
                    n_fq_bpoly_evalp_step(Bev, Beh, ctx->fqctx);
                    n_fq_poly_evalp_step(c, gammaeh, ctx->fqctx);
                }
                else
                {
                    n_fq_bpoly_eval_step(Aev, Aeh, ctx->fqctx);
                    n_fq_bpoly_eval_step(Bev, Beh, ctx->fqctx);
                    n_fq_poly_eval_step(c, gammaeh, ctx->fqctx);
                }

                if (_n_fq_is_zero(c, d))
                    goto choose_betas;
                if (Aev->length < 1 || n_bpoly_bidegree(Aev) != A->exps[0])
                    goto choose_betas;
                if (Bev->length < 1 || n_bpoly_bidegree(Bev) != B->exps[0])
                    goto choose_betas;

                success = n_fq_bpoly_gcd_brown_smprime(Gev, Abarev, Bbarev, Aev, Bev, ctx->fqctx, St);        
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
                if ((use & USE_G) && qzip_solvep(G, ZG, HG, MG, ctx) < 1)
                    goto choose_main;
                if ((use & USE_ABAR) && qzip_solvep(Abar, ZAbar, HAbar, MAbar, ctx) < 1)
                    goto choose_main;
                if ((use & USE_BBAR) && qzip_solvep(Bbar, ZBbar, HBbar, MBbar, ctx) < 1)
                    goto choose_main;
            }
            else
            {
                if ((use & USE_G) && qzip_solve(G, ZG, HG, MG, ctx) < 1)
                    goto choose_main;
                if ((use & USE_ABAR) && qzip_solve(Abar, ZAbar, HAbar, MAbar, ctx) < 1)
                    goto choose_main;
                if ((use & USE_BBAR) && qzip_solve(Bbar, ZBbar, HBbar, MBbar, ctx) < 1)
                    goto choose_main;
            }

            n_fq_poly_eval_pow(c, modulus, alphapow, ctx->fqctx);
            n_fq_inv(c, c, ctx->fqctx);
            n_fq_poly_scalar_mul_n_fq(modulus, modulus, c, ctx->fqctx);

            if ((use & USE_G) && !fq_nmod_mpolyun_interp_mcrt_sm_mpolyu(
                                      &lastdeg, Gn, G, modulus, alphapow, ctx))
            {
timeit_start(timer);
                fq_nmod_mpolyu_cvtfrom_mpolyun(rG, Gn, m, ctx);
                if (m == nvars - 1)
                {
                    success = fq_nmod_mpolyu_content_mpoly(cont, rG, ctx);
                    if (!success)
                        goto cleanup;
                    fq_nmod_mpolyu_divexact_mpoly_inplace(rG, cont, ctx);

                    if (fq_nmod_mpolyuu_divides(rAbar, A, rG, 2, ctx) &&
                        fq_nmod_mpolyuu_divides(rBbar, B, rG, 2, ctx))
                    {
timeit_stop(timer);
flint_printf("m = %wd divisibility: %wd\n", m, timer->wall);
                        break;
                    }
                }
                else
                {
                    fq_nmod_mpolyu_mul_mpoly(T, Aevals + m + 1, gammaevals + m + 1, ctx);
                    if (fq_nmod_mpolyuu_divides(rAbar, T, rG, 2, ctx))
                    {
                        fq_nmod_mpolyu_mul_mpoly(T, Bevals + m + 1, gammaevals + m + 1, ctx);
                        if (fq_nmod_mpolyuu_divides(rBbar, T, rG, 2, ctx))
                        {
                            fq_nmod_mpolyu_swap(G, rG, ctx);
                            fq_nmod_mpolyu_swap(Abar, rAbar, ctx);
                            fq_nmod_mpolyu_swap(Bbar, rBbar, ctx);
timeit_stop(timer);
flint_printf("m = %wd divisibility: %wd\n", m, timer->wall);
                            break;
                        }
                    }
                }
            }
            if ((use & USE_ABAR) && !fq_nmod_mpolyun_interp_mcrt_sm_mpolyu(
                                &lastdeg, Abarn, Abar, modulus, alphapow, ctx))
            {
                fq_nmod_mpolyu_cvtfrom_mpolyun(rAbar, Abarn, m, ctx);
                if (m == nvars - 1)
                {
                    success = fq_nmod_mpolyu_content_mpoly(cont, rAbar, ctx);
                    if (!success)
                        goto cleanup;
                    fq_nmod_mpolyu_divexact_mpoly_inplace(rAbar, cont, ctx);
                    if (fq_nmod_mpolyuu_divides(rG, A, rAbar, 2, ctx) &&
                        fq_nmod_mpolyuu_divides(rBbar, B, rG, 2, ctx))
                    {
                        break;
                    }
                }
                else
                {
                    fq_nmod_mpolyu_mul_mpoly(T, Aevals + m + 1, gammaevals + m + 1, ctx);
                    if (fq_nmod_mpolyuu_divides(rG, T, rAbar, 2, ctx))
                    {
                        fq_nmod_mpolyu_mul_mpoly(T, Bevals + m + 1, gammaevals + m + 1, ctx);
                        if (fq_nmod_mpolyuu_divides(rBbar, T, rG, 2, ctx))
                        {
                            fq_nmod_mpolyu_swap(G, rG, ctx);
                            fq_nmod_mpolyu_swap(Abar, rAbar, ctx);
                            fq_nmod_mpolyu_swap(Bbar, rBbar, ctx);
                            break;
                        }
                    }
                }
            }
            if ((use & USE_BBAR) && !fq_nmod_mpolyun_interp_mcrt_sm_mpolyu(
                            &lastdeg, Bbarn, Bbar, modulus, alphapow, ctx))
            {
                fq_nmod_mpolyu_cvtfrom_mpolyun(rBbar, Bbarn, m, ctx);
                if (m == nvars - 1)
                {
                    success = fq_nmod_mpolyu_content_mpoly(cont, rBbar, ctx);
                    if (!success)
                        goto cleanup;
                    fq_nmod_mpolyu_divexact_mpoly_inplace(rBbar, cont, ctx);
                    if (fq_nmod_mpolyuu_divides(rG, B, rBbar, 2, ctx) &&
                        fq_nmod_mpolyuu_divides(rAbar, A, rG, 2, ctx))
                    {
                        break;
                    }
                }
                else
                {
                    fq_nmod_mpolyu_mul_mpoly(T, Bevals + m + 1, gammaevals + m + 1, ctx);
                    if (fq_nmod_mpolyuu_divides(rG, T, rBbar, 2, ctx))
                    {
                        fq_nmod_mpolyu_mul_mpoly(T, Aevals + m + 1, gammaevals + m + 1, ctx);
                        if (fq_nmod_mpolyuu_divides(rAbar, T, rG, 2, ctx))
                        {
                            fq_nmod_mpolyu_swap(G, rG, ctx);
                            fq_nmod_mpolyu_swap(Abar, rAbar, ctx);
                            fq_nmod_mpolyu_swap(Bbar, rBbar, ctx);
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
    fq_nmod_mpoly_clear(cont, ctx);
    fq_nmod_mpolyu_clear(T, ctx);
    fq_nmod_mpolyu_clear(G, ctx);
    fq_nmod_mpolyu_clear(Abar, ctx);
    fq_nmod_mpolyu_clear(Bbar, ctx);
    fq_nmod_mpolyun_clear(Tn, ctx);
    fq_nmod_mpolyun_clear(Gn, ctx);
    fq_nmod_mpolyun_clear(Abarn, ctx);
    fq_nmod_mpolyun_clear(Bbarn, ctx);
    n_fq_polyun_clear(Aeh);
    n_fq_polyun_clear(Beh);
    n_fq_poly_clear(gammaeh);
    n_fq_poly_clear(modulus);
    n_fq_poly_clear(alphapow);
    n_poly_stack_clear(St->poly_stack);
    n_bpoly_stack_clear(St->bpoly_stack);

    fq_nmod_clear(start_alpha, ctx->fqctx);
    flint_free(c);

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
        fq_nmod_mpolyu_clear(Aevals + i, ctx);
        fq_nmod_mpolyu_clear(Bevals + i, ctx);
        fq_nmod_mpoly_clear(gammaevals + i, ctx);
    }
    flint_free(Aevals);
    flint_free(Bevals);
    flint_free(gammaevals);

    return success;

gcd_is_trivial:

    fq_nmod_mpolyu_one(rG, ctx);
    fq_nmod_mpolyu_set(rAbar, A, ctx);
    fq_nmod_mpolyu_set(rBbar, B, ctx);

    success = 1;
    
    goto cleanup;
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


int fq_nmod_mpolyuu_gcd_zippel_lgprime(
    fq_nmod_mpolyu_t rG, const slong * rGdegs,
    fq_nmod_mpolyu_t rAbar,
    fq_nmod_mpolyu_t rBbar,
    const fq_nmod_mpolyu_t A, const slong * Adegs,
    const fq_nmod_mpolyu_t B, const slong * Bdegs,
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
    fq_nmod_mpolyu_t T, G, Abar, Bbar;
    n_fq_polyun_t HG, HAbar, HBbar, MG, MAbar, MBbar, ZG, ZAbar, ZBbar;
    n_fq_bpoly_t Aev, Bev, Gev, Abarev, Bbarev;
    const mp_limb_t * gammaev;
    fq_nmod_mpolyun_t Tn, Gn, Abarn, Bbarn;
    slong lastdeg;
    slong cur_zip_image, req_zip_images, this_images;
    n_fq_polyun_t Aeh, Beh;
    n_fq_poly_t gammaeh;
    fq_nmod_mpolyu_struct * Aevals, * Bevals;
    fq_nmod_mpoly_struct * gammaevals;
    n_fq_poly_t modulus, alphapow;
    n_poly_bpoly_stack_t St;
    fq_nmod_t start_alpha;
    n_poly_t tmp;  /* tmp arithmetic space */
    ulong GdegboundXY, newdegXY;
    slong degxAB, degyAB;
    bad_fq_nmod_mpoly_embed_chooser_t embc;
    bad_fq_nmod_embed_struct * cur_emb;
    fq_nmod_mpoly_ctx_t lgctx;
    fq_nmod_mpolyn_t gamman;
    fq_nmod_mpolyun_t An, Bn;

    FLINT_ASSERT(bits <= FLINT_BITS);
    FLINT_ASSERT(bits == rG->bits);
    FLINT_ASSERT(bits == rAbar->bits);
    FLINT_ASSERT(bits == rBbar->bits);
    FLINT_ASSERT(bits == A->bits);
    FLINT_ASSERT(bits == B->bits);
    FLINT_ASSERT(bits == gamma->bits);

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
    fq_nmod_mpolyu_init(T, bits, smctx);
    fq_nmod_mpolyu_init(G, bits, smctx);
    fq_nmod_mpolyu_init(Abar, bits, smctx);
    fq_nmod_mpolyu_init(Bbar, bits, smctx);
    fq_nmod_mpolyun_init(Tn, bits, smctx);
    fq_nmod_mpolyun_init(Gn, bits, smctx);
    fq_nmod_mpolyun_init(Abarn, bits, smctx);
    fq_nmod_mpolyun_init(Bbarn, bits, smctx);
    n_fq_polyun_init(Aeh);
    n_fq_polyun_init(Beh);
    n_fq_poly_init(gammaeh);
    n_fq_poly_init(modulus);
    n_poly_stack_init(St->poly_stack);
    n_bpoly_stack_init(St->bpoly_stack);
    fq_nmod_mpolyun_init(An, bits, smctx);
    fq_nmod_mpolyun_init(Bn, bits, smctx);
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
    Aevals = FLINT_ARRAY_ALLOC(nvars, fq_nmod_mpolyu_struct);
    Bevals = FLINT_ARRAY_ALLOC(nvars, fq_nmod_mpolyu_struct);
    gammaevals = FLINT_ARRAY_ALLOC(nvars, fq_nmod_mpoly_struct);
    for (i = 0; i < nvars; i++)
    {
        fq_nmod_mpolyu_init(Aevals + i, bits, smctx);
        fq_nmod_mpolyu_init(Bevals + i, bits, smctx);
        fq_nmod_mpoly_init3(gammaevals + i, 0, bits, smctx);
    }
    fq_nmod_mpolyu_cvtto_mpolyun(An, A, nvars - 1, smctx);
    fq_nmod_mpolyu_cvtto_mpolyun(Bn, B, nvars - 1, smctx);
    fq_nmod_mpoly_cvtto_mpolyn(gamman, gamma, nvars - 1, smctx);

    GdegboundXY = mpoly_bivar_degrees1_chained(0, A->exps, A->length);
    GdegboundXY = mpoly_bivar_degrees1_chained(GdegboundXY, B->exps, B->length);
    degxAB = extract_exp(GdegboundXY, 0, 2);
    degyAB = extract_exp(GdegboundXY, 1, 2);

    GdegboundXY = FLINT_MIN(A->exps[0], B->exps[0]);
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

    for (i = 0; i < nvars - 1; i++)
    {
        fq_nmod_rand(alphas + i, state, lgctx->fqctx);
        if (fq_nmod_is_zero(alphas + i, lgctx->fqctx))
            fq_nmod_one(alphas + i, lgctx->fqctx);
    }

    i = nvars - 1;
    fq_nmod_mpolyun_interp_reduce_lg_mpolyu(Aevals + i, An, lgctx, smctx, cur_emb);
    fq_nmod_mpolyun_interp_reduce_lg_mpolyu(Bevals + i, Bn, lgctx, smctx, cur_emb);
    fq_nmod_mpolyn_interp_reduce_lg_mpoly(gammaevals + i, gamman, lgctx, smctx, cur_emb);
    if (fq_nmod_mpoly_is_zero(gammaevals + i, lgctx))
        goto choose_main;
    if (Aevals[i].length < 1 || Aevals[i].exps[0] != A->exps[0])
        goto choose_main;
    if (Bevals[i].length < 1 || Bevals[i].exps[0] != B->exps[0])
        goto choose_main;
    for (i--; i >= 0; i--)
    {
        fq_nmod_mpolyu_evaluate_one_fq_nmod(Aevals + i, Aevals + i + 1, i, alphas + i, lgctx);
        fq_nmod_mpolyu_evaluate_one_fq_nmod(Bevals + i, Bevals + i + 1, i, alphas + i, lgctx);
        fq_nmod_mpoly_evaluate_one_fq_nmod(gammaevals + i, gammaevals + i + 1, i, alphas + i, lgctx);
        fq_nmod_mpoly_repack_bits_inplace(gammaevals + i, bits, lgctx);
        if (fq_nmod_mpoly_is_zero(gammaevals + i, lgctx))
            goto choose_main;
        if (Aevals[i].length < 1 || Aevals[i].exps[0] != A->exps[0])
            goto choose_main;
        if (Bevals[i].length < 1 || Bevals[i].exps[0] != B->exps[0])
            goto choose_main;
    }

    m = 0;

    if (rGdegs == NULL)
    {
        use = USE_G | USE_ABAR | USE_BBAR;
    }
    else
    {
        use = mpoly_gcd_get_use_first(rGdegs[m], Adegs[m], Bdegs[m], gammadegs[m]);
    }

    fq_nmod_mpolyuu_get_n_fq_bpoly(Aev, Aevals + m, lgctx);
    fq_nmod_mpolyuu_get_n_fq_bpoly(Bev, Bevals + m, lgctx);

    success = n_fq_bpoly_gcd_brown_smprime(Gev, Abarev, Bbarev, Aev, Bev, lgctx->fqctx, St);
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

    if (nvars == 1)
    {
        fq_nmod_mpolyuun_interp_lift_lg_bpoly(&lastdeg, Gn, smctx, Gev, lgctx, cur_emb);
        fq_nmod_mpolyuun_interp_lift_lg_bpoly(&lastdeg, Abarn, smctx, Abarev, lgctx, cur_emb);
        fq_nmod_mpolyuun_interp_lift_lg_bpoly(&lastdeg, Bbarn, smctx, Bbarev, lgctx, cur_emb);

        n_fq_poly_set_fq_nmod_poly(modulus, cur_emb->h, smctx->fqctx);

        while (1)
        {
        choose_alpha_0_last:

            cur_emb = bad_fq_nmod_mpoly_embed_chooser_next(embc, lgctx, smctx, state);
            if (cur_emb == NULL)
            {
                /* ran out of primes */
                success = 0;
                goto cleanup;
            }

            lgd = fq_nmod_ctx_degree(lgctx->fqctx);

            fq_nmod_mpolyun_interp_reduce_lg_mpolyu(Aevals + m, An, lgctx, smctx, cur_emb);
            fq_nmod_mpolyun_interp_reduce_lg_mpolyu(Bevals + m, Bn, lgctx, smctx, cur_emb);
            fq_nmod_mpolyn_interp_reduce_lg_mpoly(gammaevals + m, gamman, lgctx, smctx, cur_emb);
            if (fq_nmod_mpoly_is_zero(gammaevals + m, lgctx))
                goto choose_alpha_0_last;
            if (Aevals[m].length < 1 || Aevals[m].exps[0] != A->exps[0])
                goto choose_alpha_0_last;
            if (Bevals[m].length < 1 || Bevals[m].exps[0] != B->exps[0])
                goto choose_alpha_0_last;

            fq_nmod_mpolyuu_get_n_fq_bpoly(Aev, Aevals + m, lgctx);
            fq_nmod_mpolyuu_get_n_fq_bpoly(Bev, Bevals + m, lgctx);

            success = n_fq_bpoly_gcd_brown_smprime(Gev, Abarev, Bbarev, Aev, Bev, lgctx->fqctx, St);
            if (!success)
                goto increase_degree;

            newdegXY = n_bpoly_bidegree(Gev);
            if (newdegXY > GdegboundXY)
                goto choose_alpha_0_last;
            if (newdegXY < GdegboundXY)
            {
                GdegboundXY = newdegXY;
                if (GdegboundXY == 0)
                    goto gcd_is_trivial;
                goto choose_main;
            }

            gammaev = fq_nmod_mpoly_get_nonzero_n_fq(gammaevals + m, lgctx);
            n_fq_bpoly_scalar_mul_n_fq(Gev, gammaev, lgctx->fqctx);

            if ((use & USE_G) && !fq_nmod_mpolyuun_interp_crt_lg_bpoly(
                        &lastdeg, Gn, Tn, modulus, smctx, Gev, lgctx, cur_emb))
            {
                fq_nmod_mpolyu_cvtfrom_mpolyun(rG, Gn, m, smctx);
                success = fq_nmod_mpolyu_content_mpoly(cont, rG, smctx);
                if (!success)
                    goto cleanup;
                fq_nmod_mpolyu_divexact_mpoly_inplace(rG, cont, smctx);
                if (fq_nmod_mpolyuu_divides(rAbar, A, rG, 2, smctx) &&
                    fq_nmod_mpolyuu_divides(rBbar, B, rG, 2, smctx))
                {
                    break;
                }
            }

            if ((use & USE_ABAR) && !fq_nmod_mpolyuun_interp_crt_lg_bpoly(
                  &lastdeg, Abarn, Tn, modulus, smctx, Abarev, lgctx, cur_emb))
            {
                fq_nmod_mpolyu_cvtfrom_mpolyun(rAbar, Abarn, m, smctx);
                success = fq_nmod_mpolyu_content_mpoly(cont, rAbar, smctx);
                if (!success)
                    goto cleanup;
                fq_nmod_mpolyu_divexact_mpoly_inplace(rAbar, cont, smctx);
                if (fq_nmod_mpolyuu_divides(rG, A, rAbar, 2, smctx) &&
                    fq_nmod_mpolyuu_divides(rBbar, B, rG, 2, smctx))
                {
                    break;
                }
            }

            if ((use & USE_BBAR) && !fq_nmod_mpolyuun_interp_crt_lg_bpoly(
                  &lastdeg, Bbarn, Tn, modulus, smctx, Bbarev, lgctx, cur_emb))
            {
                fq_nmod_mpolyu_cvtfrom_mpolyun(rBbar, Bbarn, m, smctx);
                success = fq_nmod_mpolyu_content_mpoly(cont, rBbar, smctx);
                if (!success)
                    goto cleanup;
                fq_nmod_mpolyu_divexact_mpoly_inplace(rBbar, cont, smctx);
                if (fq_nmod_mpolyuu_divides(rG, B, rBbar, 2, smctx) &&
                    fq_nmod_mpolyuu_divides(rAbar, A, rG, 2, smctx))
                {
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

    fq_nmod_mpolyuun_interp_lift_sm_bpoly(Gn, Gev, lgctx);
    fq_nmod_mpolyuun_interp_lift_sm_bpoly(Abarn, Abarev, lgctx);
    fq_nmod_mpolyuun_interp_lift_sm_bpoly(Bbarn, Bbarev, lgctx);

    FLINT_ASSERT(tmp->alloc >= lgd);
    n_fq_poly_one(modulus, lgctx->fqctx);
    n_fq_set_fq_nmod(tmp->coeffs, alphas + m, lgctx->fqctx);
    n_fq_poly_shift_left_scalar_submul(modulus, 1, tmp->coeffs, lgctx->fqctx);

    fq_nmod_set(start_alpha, alphas + m, lgctx->fqctx);

    while (1)
    {
choose_alpha_0:

        fq_nmod_next_not_zero(alphas + m, lgctx->fqctx);
        if (fq_nmod_equal(alphas + m, start_alpha, lgctx->fqctx))
            goto increase_degree;

        FLINT_ASSERT(alphapow->alloc >= lgd*2);
        _n_fq_one(alphapow->coeffs + lgd*0, lgd);
        n_fq_set_fq_nmod(alphapow->coeffs + lgd*1, alphas + m, lgctx->fqctx);
        alphapow->length = 2;

        fq_nmod_mpolyu_evaluate_one_fq_nmod(Aevals + m, Aevals + m + 1, m, alphas + m, lgctx);
        fq_nmod_mpolyu_evaluate_one_fq_nmod(Bevals + m, Bevals + m + 1, m, alphas + m, lgctx);
        fq_nmod_mpoly_evaluate_one_fq_nmod(gammaevals + m, gammaevals + m + 1, m, alphas + m, lgctx);
        if (fq_nmod_mpoly_is_zero(gammaevals + m, lgctx))
            goto choose_alpha_0;
        if (Aevals[m].length < 1 || Aevals[m].exps[0] != A->exps[0])
            goto choose_alpha_0;
        if (Bevals[m].length < 1 || Bevals[m].exps[0] != B->exps[0])
            goto choose_alpha_0;

        fq_nmod_mpolyuu_get_n_fq_bpoly(Aev, Aevals + m, lgctx);
        fq_nmod_mpolyuu_get_n_fq_bpoly(Bev, Bevals + m, lgctx);

        success = n_fq_bpoly_gcd_brown_smprime(Gev, Abarev, Bbarev, Aev, Bev, lgctx->fqctx, St);
        if (!success)
            goto increase_degree;;

        newdegXY = n_bpoly_bidegree(Gev);
        if (newdegXY > GdegboundXY)
            goto choose_alpha_0;
        if (newdegXY < GdegboundXY)
        {
            GdegboundXY = newdegXY;
            if (GdegboundXY == 0)
                goto gcd_is_trivial;
            goto choose_main;
        }

        gammaev = fq_nmod_mpoly_get_nonzero_n_fq(gammaevals + m, lgctx);
        n_fq_bpoly_scalar_mul_n_fq(Gev, gammaev, lgctx->fqctx);

        n_poly_fit_length(tmp, lgd);
        n_fq_poly_eval_pow(tmp->coeffs, modulus, alphapow, lgctx->fqctx);
        n_fq_inv(tmp->coeffs, tmp->coeffs, lgctx->fqctx);
        n_fq_poly_scalar_mul_n_fq(modulus, modulus, tmp->coeffs, lgctx->fqctx);

        if ((use & USE_G) && !fq_nmod_mpolyuun_interp_crt_sm_bpoly(
                              &lastdeg, Gn, Tn, Gev, modulus, alphapow, lgctx))
        {
            fq_nmod_mpolyu_cvtfrom_mpolyun(G, Gn, m, lgctx);
            fq_nmod_mpolyu_mul_mpoly(T, Aevals + m + 1, gammaevals + m + 1, lgctx);
            if (fq_nmod_mpolyuu_divides(Abar, T, G, 2, lgctx))
            {
                fq_nmod_mpolyu_mul_mpoly(T, Bevals + m + 1, gammaevals + m + 1, lgctx);
                if (fq_nmod_mpolyuu_divides(Bbar, T, G, 2, lgctx))
                {
                    break;
                }
            }
        }

        if ((use & USE_ABAR) && !fq_nmod_mpolyuun_interp_crt_sm_bpoly(
                        &lastdeg, Abarn, Tn, Abarev, modulus, alphapow, lgctx))
        {
            fq_nmod_mpolyu_cvtfrom_mpolyun(Abar, Abarn, m, lgctx);
            fq_nmod_mpolyu_mul_mpoly(T, Aevals + m + 1, gammaevals + m + 1, lgctx);
            if (fq_nmod_mpolyuu_divides(G, T, Abar, 2, lgctx))
            {
                fq_nmod_mpolyu_mul_mpoly(T, Bevals + m + 1, gammaevals + m + 1, lgctx);
                if (fq_nmod_mpolyuu_divides(Bbar, T, G, 2, lgctx))
                    break;
            }
        }

        if ((use & USE_BBAR) && !fq_nmod_mpolyuun_interp_crt_sm_bpoly(
                          &lastdeg, Bbarn, Tn, Bbarev, modulus, alphapow, lgctx))
        {
            fq_nmod_mpolyu_cvtfrom_mpolyun(Bbar, Bbarn, m, lgctx);
            fq_nmod_mpolyu_mul_mpoly(T, Bevals + m + 1, gammaevals + m + 1, lgctx);
            if (fq_nmod_mpolyuu_divides(G, T, Bbar, 2, lgctx))
            {
                fq_nmod_mpolyu_mul_mpoly(T, Aevals + m + 1, gammaevals + m + 1, lgctx);
                if (fq_nmod_mpolyuu_divides(Abar, T, G, 2, lgctx))
                    break;
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

    for (m = 1; m < nvars - 1; m++)
    {
        /* G, Abar, Bbar are in Fq[gen(0), ..., gen(m - 1)] */
        fq_nmod_mpolyun_interp_lift_sm_mpolyu(Gn, G, lgctx);
        fq_nmod_mpolyun_interp_lift_sm_mpolyu(Abarn, Abar, lgctx);
        fq_nmod_mpolyun_interp_lift_sm_mpolyu(Bbarn, Bbar, lgctx);

        if (rGdegs == NULL)
        {
            use = USE_G | USE_ABAR | USE_BBAR;
        }
        else
        {
            slong numABgamma = gamma[m + 1].length +
                               fq_nmod_mpolyu_total_length(Aevals + m + 1) +
                               fq_nmod_mpolyu_total_length(Bevals + m + 1);

            use = fq_nmod_mpoly_gcd_get_use(rGdegs[m], Adegs[m], Bdegs[m],
                      gammadegs[m], degxAB, degyAB, numABgamma, G, Abar, Bbar);
        }

        FLINT_ASSERT(tmp->alloc >= lgd);
        n_fq_poly_one(modulus, lgctx->fqctx);
        n_fq_set_fq_nmod(tmp->coeffs, alphas + m, lgctx->fqctx);
        n_fq_poly_shift_left_scalar_submul(modulus, 1, tmp->coeffs, lgctx->fqctx);

        fq_nmod_set(start_alpha, alphas + m, lgctx->fqctx);

choose_betas_m:

        /* only beta[0], beta[1], ..., beta[m - 1] will be used */
        for (i = 0; i < nvars; i++)
        {
            fq_nmod_rand(betas + i, state, lgctx->fqctx);
            if (fq_nmod_is_zero(betas + i, lgctx->fqctx))
                fq_nmod_one(betas + i, lgctx->fqctx);
        }

        req_zip_images = 2;
        if (use & USE_G)
        {
            this_images = fq_nmod_mpolyu_set_zip_form(HG, MG, G, betas, m, lgctx);
            req_zip_images = FLINT_MAX(req_zip_images, this_images);
        }
        if (use & USE_ABAR)
        {
            this_images =  fq_nmod_mpolyu_set_zip_form(HAbar, MAbar, Abar, betas, m, lgctx);
            req_zip_images = FLINT_MAX(req_zip_images, this_images);            
        }
        if (use & USE_BBAR)
        {
            this_images =  fq_nmod_mpolyu_set_zip_form(HBbar, MBbar, Bbar, betas, m, lgctx);
            req_zip_images = FLINT_MAX(req_zip_images, this_images);            
        }

        while (1)
        {
choose_alpha_m:

            fq_nmod_next_not_zero(alphas + m, lgctx->fqctx);
            if (fq_nmod_equal(alphas + m, start_alpha, lgctx->fqctx))
                goto choose_main;

            FLINT_ASSERT(alphapow->alloc >= lgd*2);
            alphapow->length = 2;
            _n_fq_one(alphapow->coeffs + lgd*0, lgd);
            n_fq_set_fq_nmod(alphapow->coeffs + lgd*1, alphas + m, lgctx->fqctx);

            fq_nmod_mpolyu_evaluate_one_fq_nmod(Aevals + m, Aevals + m + 1, m, alphas + m, lgctx);
            fq_nmod_mpolyu_evaluate_one_fq_nmod(Bevals + m, Bevals + m + 1, m, alphas + m, lgctx);
            fq_nmod_mpoly_evaluate_one_fq_nmod(gammaevals + m, gammaevals + m + 1, m, alphas + m, lgctx);
            if (fq_nmod_mpoly_is_zero(gammaevals + m, lgctx))
                goto choose_alpha_m;
            if (Aevals[m].length < 1 || Aevals[m].exps[0] != A->exps[0])
                goto choose_alpha_m;
            if (Bevals[m].length < 1 || Bevals[m].exps[0] != B->exps[0])
                goto choose_alpha_m;

            qzip_start(ZG, HG, req_zip_images, lgctx->fqctx);
            qzip_start(ZAbar, HAbar, req_zip_images, lgctx->fqctx);
            qzip_start(ZBbar, HBbar, req_zip_images, lgctx->fqctx);

            fq_nmod_mpolyu_set_eval_helper(Aeh, Aevals + m, betas, m, lgctx);
            fq_nmod_mpolyu_set_eval_helper(Beh, Bevals + m, betas, m, lgctx);
            fq_nmod_mpoly_set_eval_helper(gammaeh, gammaevals + m, betas, m, lgctx);

            for (cur_zip_image = 0; cur_zip_image < req_zip_images; cur_zip_image++)
            {
                n_fq_bpoly_eval_step(Aev, Aeh, lgctx->fqctx);
                n_fq_bpoly_eval_step(Bev, Beh, lgctx->fqctx);
                FLINT_ASSERT(tmp->alloc >= lgd);
                n_fq_poly_eval_step(tmp->coeffs, gammaeh, lgctx->fqctx);

                if (_n_fq_is_zero(tmp->coeffs, lgd))
                    goto choose_betas_m;
                if (Aev->length < 1 || n_bpoly_bidegree(Aev) != A->exps[0])
                    goto choose_betas_m;
                if (Bev->length < 1 || n_bpoly_bidegree(Bev) != B->exps[0])
                    goto choose_betas_m;

                success = n_fq_bpoly_gcd_brown_smprime(Gev, Abarev, Bbarev, Aev, Bev, lgctx->fqctx, St);
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

            if ((use & USE_G) && qzip_solve(G, ZG, HG, MG, lgctx) < 1)
                goto choose_main;
            if ((use & USE_ABAR) && qzip_solve(Abar, ZAbar, HAbar, MAbar, lgctx) < 1)
                goto choose_main;
            if ((use & USE_BBAR) && qzip_solve(Bbar, ZBbar, HBbar, MBbar, lgctx) < 1)
                goto choose_main;

            FLINT_ASSERT(n_fq_poly_degree(modulus) > 0);

            n_poly_fit_length(tmp, lgd);
            n_fq_poly_eval_pow(tmp->coeffs, modulus, alphapow, lgctx->fqctx);
            n_fq_inv(tmp->coeffs, tmp->coeffs, lgctx->fqctx);
            n_fq_poly_scalar_mul_n_fq(modulus, modulus, tmp->coeffs, lgctx->fqctx);

            if ((use & USE_G) && !fq_nmod_mpolyun_interp_mcrt_sm_mpolyu(
                                      &lastdeg, Gn, G, modulus, alphapow, lgctx))
            {
                fq_nmod_mpolyu_cvtfrom_mpolyun(rG, Gn, m, lgctx);
                fq_nmod_mpolyu_mul_mpoly(T, Aevals + m + 1, gammaevals + m + 1, lgctx);
                if (fq_nmod_mpolyuu_divides(rAbar, T, rG, 2, lgctx))
                {
                    fq_nmod_mpolyu_mul_mpoly(T, Bevals + m + 1, gammaevals + m + 1, lgctx);
                    if (fq_nmod_mpolyuu_divides(rBbar, T, rG, 2, lgctx))
                    {
                        fq_nmod_mpolyu_swap(G, rG, lgctx);
                        fq_nmod_mpolyu_swap(Abar, rAbar, lgctx);
                        fq_nmod_mpolyu_swap(Bbar, rBbar, lgctx);
                        break;
                    }
                }
            }
            if ((use & USE_ABAR) && !fq_nmod_mpolyun_interp_mcrt_sm_mpolyu(
                              &lastdeg, Abarn, Abar, modulus, alphapow, lgctx))
            {
                fq_nmod_mpolyu_cvtfrom_mpolyun(rAbar, Abarn, m, lgctx);
                fq_nmod_mpolyu_mul_mpoly(T, Aevals + m + 1, gammaevals + m + 1, lgctx);
                if (fq_nmod_mpolyuu_divides(rG, T, rAbar, 2, lgctx))
                {
                    fq_nmod_mpolyu_mul_mpoly(T, Bevals + m + 1, gammaevals + m + 1, lgctx);
                    if (fq_nmod_mpolyuu_divides(rBbar, T, rG, 2, lgctx))
                    {
                        fq_nmod_mpolyu_swap(G, rG, lgctx);
                        fq_nmod_mpolyu_swap(Abar, rAbar, lgctx);
                        fq_nmod_mpolyu_swap(Bbar, rBbar, lgctx);
                        break;
                    }
                }
            }
            if ((use & USE_BBAR) && !fq_nmod_mpolyun_interp_mcrt_sm_mpolyu(
                              &lastdeg, Bbarn, Bbar, modulus, alphapow, lgctx))
            {
                fq_nmod_mpolyu_cvtfrom_mpolyun(rBbar, Bbarn, m, lgctx);
                fq_nmod_mpolyu_mul_mpoly(T, Bevals + m + 1, gammaevals + m + 1, lgctx);
                if (fq_nmod_mpolyuu_divides(rG, T, rBbar, 2, lgctx))
                {
                    fq_nmod_mpolyu_mul_mpoly(T, Aevals + m + 1, gammaevals + m + 1, lgctx);
                    if (fq_nmod_mpolyuu_divides(rAbar, T, rG, 2, lgctx))
                    {
                        fq_nmod_mpolyu_swap(G, rG, lgctx);
                        fq_nmod_mpolyu_swap(Abar, rAbar, lgctx);
                        fq_nmod_mpolyu_swap(Bbar, rBbar, lgctx);
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

        fq_nmod_mpolyun_interp_lift_lg_mpolyu(Gn, smctx, G, lgctx, cur_emb);
        fq_nmod_mpolyun_interp_lift_lg_mpolyu(Abarn, smctx, Abar, lgctx, cur_emb);
        fq_nmod_mpolyun_interp_lift_lg_mpolyu(Bbarn, smctx, Bbar, lgctx, cur_emb);

        n_fq_poly_set_fq_nmod_poly(modulus, cur_emb->h, smctx->fqctx);

        if (rGdegs == NULL)
        {
            use = USE_G | USE_ABAR | USE_BBAR;
        }
        else
        {
            slong numABgamma = gamma->length + 
                               fq_nmod_mpolyu_total_length(A) + 
                               fq_nmod_mpolyu_total_length(B);

            use = fq_nmod_mpoly_gcd_get_use(rGdegs[m], Adegs[m], Bdegs[m],
                      gammadegs[m], degxAB, degyAB, numABgamma, G, Abar, Bbar);
        }

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

            lgd = fq_nmod_ctx_degree(lgctx->fqctx);
            n_poly_fit_length(tmp, lgd);
            n_poly_fit_length(alphapow, 2*lgd);

            fq_nmod_mpolyun_interp_reduce_lg_mpolyu(Aevals + m, An, lgctx, smctx, cur_emb);
            fq_nmod_mpolyun_interp_reduce_lg_mpolyu(Bevals + m, Bn, lgctx, smctx, cur_emb);
            fq_nmod_mpolyn_interp_reduce_lg_mpoly(gammaevals + m, gamman, lgctx, smctx, cur_emb);
            if (fq_nmod_mpoly_is_zero(gammaevals + m, lgctx))
                goto choose_alpha_last;
            if (Aevals[m].length < 1 || Aevals[m].exps[0] != A->exps[0])
                goto choose_alpha_last;
            if (Bevals[m].length < 1 || Bevals[m].exps[0] != B->exps[0])
                goto choose_alpha_last;

choose_betas_last:

            /* only beta[0], beta[1], ..., beta[m - 1] will be used */
            for (i = 0; i < nvars; i++)
            {
                fq_nmod_rand(betas + i, state, lgctx->fqctx);
                if (fq_nmod_is_zero(betas + i, lgctx->fqctx))
                    fq_nmod_one(betas + i, lgctx->fqctx);
            }

            req_zip_images = 2;
            if (use & USE_G)
            {
                this_images = fq_nmod_mpolyu_set_zip_form(HG, MG, G, betas, m, lgctx);
                req_zip_images = FLINT_MAX(req_zip_images, this_images);
            }
            if (use & USE_ABAR)
            {
                this_images =  fq_nmod_mpolyu_set_zip_form(HAbar, MAbar, Abar, betas, m, lgctx);
                req_zip_images = FLINT_MAX(req_zip_images, this_images);            
            }
            if (use & USE_BBAR)
            {
                this_images =  fq_nmod_mpolyu_set_zip_form(HBbar, MBbar, Bbar, betas, m, lgctx);
                req_zip_images = FLINT_MAX(req_zip_images, this_images);            
            }

            qzip_start(ZG, HG, req_zip_images, lgctx->fqctx);
            qzip_start(ZAbar, HAbar, req_zip_images, lgctx->fqctx);
            qzip_start(ZBbar, HBbar, req_zip_images, lgctx->fqctx);

            fq_nmod_mpolyu_set_eval_helper(Aeh, Aevals + m, betas, m, lgctx);
            fq_nmod_mpolyu_set_eval_helper(Beh, Bevals + m, betas, m, lgctx);
            fq_nmod_mpoly_set_eval_helper(gammaeh, gammaevals + m, betas, m, lgctx);

            for (cur_zip_image = 0; cur_zip_image < req_zip_images; cur_zip_image++)
            {
                n_fq_bpoly_eval_step(Aev, Aeh, lgctx->fqctx);
                n_fq_bpoly_eval_step(Bev, Beh, lgctx->fqctx);
                FLINT_ASSERT(tmp->alloc >= lgd);
                n_fq_poly_eval_step(tmp->coeffs, gammaeh, lgctx->fqctx);

                if (_n_fq_is_zero(tmp->coeffs, lgd))
                    goto choose_betas_last;
                if (Aev->length < 1 || n_bpoly_bidegree(Aev) != A->exps[0])
                    goto choose_betas_last;
                if (Bev->length < 1 || n_bpoly_bidegree(Bev) != B->exps[0])
                    goto choose_betas_last;


                success = n_fq_bpoly_gcd_brown_smprime(Gev, Abarev, Bbarev, Aev, Bev, lgctx->fqctx, St);        
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
            if ((use & USE_G) && qzip_solve(G, ZG, HG, MG, lgctx) < 1)
                goto choose_main;
            if ((use & USE_ABAR) && qzip_solve(Abar, ZAbar, HAbar, MAbar, lgctx) < 1)
                goto choose_main;
            if ((use & USE_BBAR) && qzip_solve(Bbar, ZBbar, HBbar, MBbar, lgctx) < 1)
                goto choose_main;

            FLINT_ASSERT(n_fq_poly_degree(modulus) > 0);

            if ((use & USE_G) && !newfq_nmod_mpolyun_interp_mcrt_lg_mpolyu(
                              &lastdeg, Gn, smctx, modulus, G, lgctx, cur_emb))
            {
                fq_nmod_mpolyu_cvtfrom_mpolyun(rG, Gn, m, smctx);
                success = fq_nmod_mpolyu_content_mpoly(cont, rG, smctx);
                if (!success)
                    goto cleanup;
                fq_nmod_mpolyu_divexact_mpoly_inplace(rG, cont, smctx);
                if (fq_nmod_mpolyuu_divides(rAbar, A, rG, 2, smctx) &&
                    fq_nmod_mpolyuu_divides(rBbar, B, rG, 2, smctx))
                {
                    break;
                }
            }
            if ((use & USE_ABAR) && !newfq_nmod_mpolyun_interp_mcrt_lg_mpolyu(
                        &lastdeg, Abarn, smctx, modulus, Abar, lgctx, cur_emb))
            {
                fq_nmod_mpolyu_cvtfrom_mpolyun(rAbar, Abarn, m, smctx);
                success = fq_nmod_mpolyu_content_mpoly(cont, rAbar, smctx);
                if (!success)
                    goto cleanup;
                fq_nmod_mpolyu_divexact_mpoly_inplace(rAbar, cont, smctx);
                if (fq_nmod_mpolyuu_divides(rG, A, rAbar, 2, smctx) &&
                    fq_nmod_mpolyuu_divides(rBbar, B, rG, 2, smctx))
                {
                    break;
                }
            }
            if ((use & USE_BBAR) && !newfq_nmod_mpolyun_interp_mcrt_lg_mpolyu(
                        &lastdeg, Bbarn, smctx, modulus, Bbar, lgctx, cur_emb))
            {
                fq_nmod_mpolyu_cvtfrom_mpolyun(rBbar, Bbarn, m, smctx);
                success = fq_nmod_mpolyu_content_mpoly(cont, rBbar, smctx);
                if (!success)
                    goto cleanup;
                fq_nmod_mpolyu_divexact_mpoly_inplace(rBbar, cont, smctx);
                if (fq_nmod_mpolyuu_divides(rG, B, rBbar, 2, smctx) &&
                    fq_nmod_mpolyuu_divides(rAbar, A, rG, 2, smctx))
                {
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

    fq_nmod_mpolyun_clear(An, smctx);
    fq_nmod_mpolyun_clear(Bn, smctx);
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
    fq_nmod_mpolyu_clear(T, smctx);
    fq_nmod_mpolyu_clear(G, smctx);
    fq_nmod_mpolyu_clear(Abar, smctx);
    fq_nmod_mpolyu_clear(Bbar, smctx);
    fq_nmod_mpolyun_clear(Tn, smctx);
    fq_nmod_mpolyun_clear(Gn, smctx);
    fq_nmod_mpolyun_clear(Abarn, smctx);
    fq_nmod_mpolyun_clear(Bbarn, smctx);
    n_fq_polyun_clear(Aeh);
    n_fq_polyun_clear(Beh);
    n_fq_poly_clear(gammaeh);
    n_fq_poly_clear(modulus);
    n_poly_stack_clear(St->poly_stack);
    n_bpoly_stack_clear(St->bpoly_stack);

    n_poly_clear(tmp);
    n_fq_poly_clear(alphapow);
    fq_nmod_clear(start_alpha, smctx->fqctx);

    flint_randclear(state);

    for (i = 0; i < nvars; i++)
    {
        fq_nmod_mpolyu_clear(Aevals + i, smctx);
        fq_nmod_mpolyu_clear(Bevals + i, smctx);
        fq_nmod_mpoly_clear(gammaevals + i, smctx);
    }
    flint_free(Aevals);
    flint_free(Bevals);
    flint_free(gammaevals);

    return success;

gcd_is_trivial:

    fq_nmod_mpolyu_one(rG, smctx);
    fq_nmod_mpolyu_set(rAbar, A, smctx);
    fq_nmod_mpolyu_set(rBbar, B, smctx);

    success = 1;
    
    goto cleanup;
}
