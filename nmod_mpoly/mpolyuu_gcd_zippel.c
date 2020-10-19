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
#include "fq_nmod_mpoly.h"
#include "fq_nmod_mpoly_factor.h"

void usleep(ulong);

void nmod_mpolyn_interp_lift_lg_mpoly(
    nmod_mpolyn_t A,
    const nmod_mpoly_ctx_t ctx,
    fq_nmod_mpoly_t Ap,
    const fq_nmod_mpoly_ctx_t ctxp);

int fq_nmod_mpolyn_interp_mcrt_sm_mpoly(
    slong * lastdeg_,
    fq_nmod_mpolyn_t F,
    fq_nmod_mpoly_t A,  /* could have zero coeffs */
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

int fq_nmod_mpolyuun_interp_crt_sm_bpoly(
    slong * lastdeg,
    fq_nmod_mpolyun_t F,
    fq_nmod_mpolyun_t T,
    n_fq_bpoly_t A,
    n_fq_poly_t modulus,
    n_fq_poly_t alphapow,
    const fq_nmod_mpoly_ctx_t ctx);


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
    slong i, Ei;

    FLINT_ASSERT(E != A);
    FLINT_ASSERT(E->bits == A->bits);

    nmod_mpolyu_fit_length(E, A->length, ctx);
    Ei = 0;
    for (i = 0; i < A->length; i++)
    {
        E->exps[Ei] = A->exps[i];
        nmod_mpoly_evaluate_one_ui(E->coeffs + Ei, A->coeffs + i, var, alpha, ctx);
        nmod_mpoly_repack_bits_inplace(E->coeffs + Ei, E->bits, ctx);
        Ei += !nmod_mpoly_is_zero(E->coeffs + Ei, ctx);
    }

    E->length = Ei;
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

static void nmod_mpolyn_interp_lift_sm_bpoly(
    nmod_mpolyn_t F,
    n_bpoly_t A,
    const nmod_mpoly_ctx_t ctx)
{
    slong N = mpoly_words_per_exp_sp(F->bits, ctx->minfo);
    slong i, j, Fi;
    slong off0, shift0, off1, shift1;

    mpoly_gen_offset_shift_sp(&off0, &shift0, 0, F->bits, ctx->minfo);
    mpoly_gen_offset_shift_sp(&off1, &shift1, 1, F->bits, ctx->minfo);

    Fi = 0;
    for (i = A->length - 1; i >= 0; i--)
    {
        n_poly_struct * Ai = A->coeffs + i;
        for (j = Ai->length - 1; j >= 0; j--)
        {
            if (Ai->coeffs[j] == 0)
                continue;

            nmod_mpolyn_fit_length(F, Fi + 1, ctx);
            mpoly_monomial_zero(F->exps + N*Fi, N);
            (F->exps + N*Fi)[off0] += (i << shift0);
            (F->exps + N*Fi)[off1] += (j << shift1);
            n_poly_set_ui((n_poly_struct *)(F->coeffs + Fi), Ai->coeffs[j]);
            Fi++;
        }
    }

    F->length = Fi;
}


static void fq_nmod_mpolyn_interp_lift_sm_bpoly(
    fq_nmod_mpolyn_t F,
    n_fq_bpoly_t A,
    const fq_nmod_mpoly_ctx_t ctx)
{
    slong d = fq_nmod_ctx_degree(ctx->fqctx);
    slong N = mpoly_words_per_exp_sp(F->bits, ctx->minfo);
    slong i, j, Fi;
    slong off0, shift0, off1, shift1;

    mpoly_gen_offset_shift_sp(&off0, &shift0, 0, F->bits, ctx->minfo);
    mpoly_gen_offset_shift_sp(&off1, &shift1, 1, F->bits, ctx->minfo);

    Fi = 0;
    for (i = A->length - 1; i >= 0; i--)
    {
        n_fq_poly_struct * Ai = A->coeffs + i;
        for (j = Ai->length - 1; j >= 0; j--)
        {
            if (_n_fq_is_zero(Ai->coeffs + d*j, d))
                continue;

            fq_nmod_mpolyn_fit_length(F, Fi + 1, ctx);
            mpoly_monomial_zero(F->exps + N*Fi, N);
            (F->exps + N*Fi)[off0] += (i << shift0);
            (F->exps + N*Fi)[off1] += (j << shift1);
            n_fq_poly_set_n_fq(F->coeffs + Fi, Ai->coeffs + d*j, ctx->fqctx);
            Fi++;
        }
    }

    F->length = Fi;
}


static void nmod_mpolyn_interp_lift_lg_bpoly(
    slong * lastdeg_,
    nmod_mpolyn_t F,
    const nmod_mpoly_ctx_t smctx,
    n_fq_bpoly_t A,
    const fq_nmod_mpoly_ctx_t lgctx)
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

            nmod_mpolyn_fit_length(F, Fi + 1, smctx);
            mpoly_monomial_zero(F->exps + N*Fi, N);
            (F->exps + N*Fi)[off0] += (i << shift0);
            (F->exps + N*Fi)[off1] += (j << shift1);
            n_fq_get_fq_nmod(F->coeffs + Fi, Ai->coeffs + lgd*j, lgctx->fqctx);
            lastdeg = FLINT_MAX(lastdeg, nmod_poly_degree(F->coeffs + Fi));

            Fi++;
        }
    }

    F->length = Fi;

    *lastdeg_ = lastdeg;
}




static ulong _mpoly_bidegree(
    const ulong * Aexps,
    flint_bitcnt_t Abits,
    const mpoly_ctx_t mctx)
{
    slong off0, shift0, off1, shift1;
    ulong mask = (-UWORD(1)) >> (FLINT_BITS - Abits);

    mpoly_gen_offset_shift_sp(&off0, &shift0, 0, Abits, mctx);
    mpoly_gen_offset_shift_sp(&off1, &shift1, 1, Abits, mctx);

    return pack_exp2((Aexps[off0] >> shift0) & mask,
                     (Aexps[off1] >> shift1) & mask);
}


static ulong _nmod_mpoly_bidegree(
    const nmod_mpoly_t A,
    const nmod_mpoly_ctx_t ctx)
{
    FLINT_ASSERT(A->length > 0);
    return _mpoly_bidegree(A->exps, A->bits, ctx->minfo);
}
/*
static ulong _nmod_mpolyn_bidegree(
    const nmod_mpolyn_t A,
    const nmod_mpoly_ctx_t ctx)
{
    FLINT_ASSERT(A->length > 0);
    return _mpoly_bidegree(A->exps, A->bits, ctx->minfo);
}
*/
static ulong _fq_nmod_mpoly_bidegree(
    const fq_nmod_mpoly_t A,
    const fq_nmod_mpoly_ctx_t ctx)
{
    FLINT_ASSERT(A->length > 0);
    return _mpoly_bidegree(A->exps, A->bits, ctx->minfo);
}



static nmod_poly_struct * getcf(nmod_mpolyn_t A)
{
    FLINT_ASSERT(A->length == 1);
    FLINT_ASSERT(A->exps[0] == 0);
    return A->coeffs + 0;
}

static n_poly_struct * getcfn(nmod_mpolyn_t A)
{
    FLINT_ASSERT(A->length == 1);
    FLINT_ASSERT(A->exps[0] == 0);
    return (n_poly_struct *) (A->coeffs + 0);
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

mp_limb_t n_poly_mod_eval_pow(n_poly_t P, n_poly_t alphapow, nmod_t ctx)
{
    int nlimbs = _nmod_vec_dot_bound_limbs(P->length, ctx);
    return _n_poly_eval_pow(P, alphapow, nlimbs, ctx);
}


static int nmod_mpolyn_interp_crt_sm_bpoly(
    slong * lastdeg,
    nmod_mpolyn_t F,
    nmod_mpolyn_t T,
    n_bpoly_t A,
    n_poly_t modulus,
    n_poly_t alphapow,
    const nmod_mpoly_ctx_t ctx)
{
    int changed = 0;
    int nlimbs = _nmod_vec_dot_bound_limbs(modulus->length, ctx->ffinfo->mod);
    slong N = mpoly_words_per_exp(F->bits, ctx->minfo);
    slong off0, shift0, off1, shift1;
    n_poly_struct * Acoeffs = A->coeffs;
    slong Fi, Ti, Ai, ai;
    slong Flen = F->length;
    ulong * Fexps = F->exps;
    nmod_poly_struct * Fcoeffs = F->coeffs;
    ulong * Texps = T->exps;
    nmod_poly_struct * Tcoeffs = T->coeffs;
    mp_limb_t v;
    ulong Fexpi, mask;

    mask = (-UWORD(1)) >> (FLINT_BITS - F->bits);
    mpoly_gen_offset_shift_sp(&off0, &shift0, 0, F->bits, ctx->minfo);
    mpoly_gen_offset_shift_sp(&off1, &shift1, 1, F->bits, ctx->minfo);

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
            nmod_mpolyn_fit_length(T, Ti + extra + 1, ctx);
            Tcoeffs = T->coeffs;
            Texps = T->exps;
        }

        if (Fi < Flen)
            Fexpi = pack_exp2(((Fexps + N*Fi)[off0]>>shift0)&mask,
                              ((Fexps + N*Fi)[off1]>>shift1)&mask);
        else
            Fexpi = 0;

        if (Fi < Flen && Ai >= 0 && Fexpi == pack_exp2(Ai, ai))
        {
            /* F term ok, A term ok */
            mpoly_monomial_set(Texps + N*Ti, Fexps + N*Fi, N);

            v = _n_poly_eval_pow((n_poly_struct *)(Fcoeffs + Fi), alphapow, nlimbs, ctx->ffinfo->mod);
            v = nmod_sub(Acoeffs[Ai].coeffs[ai], v, ctx->ffinfo->mod);
            if (v != 0)
            {
                changed = 1;
                n_poly_mod_scalar_addmul_nmod((n_poly_struct*)(Tcoeffs + Ti),(n_poly_struct *)(Fcoeffs + Fi), modulus, v, ctx->ffinfo->mod);
            }
            else
            {
                n_poly_set((n_poly_struct*)(Tcoeffs + Ti), (n_poly_struct*)(Fcoeffs + Fi));
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
        else if (Ai >= 0 && (Fi >= Flen || Fexpi < pack_exp2(Ai, ai)))
        {
            /* F term missing, A term ok */

            mpoly_monomial_zero(Texps + N*Ti, N);
            (Texps + N*Ti)[off0] += (Ai << shift0);
            (Texps + N*Ti)[off1] += (ai << shift1);

            changed = 1;
            _n_poly_mod_scalar_mul_nmod((n_poly_struct*)(Tcoeffs + Ti), modulus, Acoeffs[Ai].coeffs[ai], ctx->ffinfo->mod);

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
            FLINT_ASSERT(Fi < Flen && (Ai < 0 || Fexpi > pack_exp2(Ai, ai)));
            /* F term ok, Aterm missing */
            mpoly_monomial_set(Texps + N*Ti, Fexps + N*Fi, N);

            v = _n_poly_eval_pow((n_poly_struct *)(Fcoeffs + Fi), alphapow, nlimbs, ctx->ffinfo->mod);
            if (v != 0)
            {
                changed = 1;
                v = nmod_neg(v, ctx->ffinfo->mod);
                n_poly_mod_scalar_addmul_nmod((n_poly_struct*)(Tcoeffs + Ti), (n_poly_struct*)(Fcoeffs + Fi), modulus, v, ctx->ffinfo->mod);
            }
            else
            {
                nmod_poly_set(Tcoeffs + Ti, Fcoeffs + Fi);
            }

            Fi++;
        }

        FLINT_ASSERT(!n_poly_is_zero((n_poly_struct *)(Tcoeffs + Ti)));
        *lastdeg = FLINT_MAX(*lastdeg, n_poly_degree((n_poly_struct*)(Tcoeffs + Ti)));
        Ti++;
    }

    T->length = Ti;

    if (changed)
        nmod_mpolyn_swap(T, F);

    return changed;
}


static int fq_nmod_mpolyn_interp_crt_sm_bpoly(
    slong * lastdeg,
    fq_nmod_mpolyn_t F,
    fq_nmod_mpolyn_t T,
    const n_fq_bpoly_t A,
    const n_fq_poly_t modulus,
    n_fq_poly_t alphapow,
    const fq_nmod_mpoly_ctx_t ctx)
{
    int changed = 0;
    slong d = fq_nmod_ctx_degree(ctx->fqctx);
    slong N = mpoly_words_per_exp(F->bits, ctx->minfo);
    slong off0, shift0, off1, shift1;
    n_poly_struct * Acoeffs = A->coeffs;
    slong Fi, Ti, Ai, ai;
    slong Flen = F->length;
    ulong * Fexps = F->exps;
    n_fq_poly_struct * Fcoeffs = F->coeffs;
    ulong * Texps = T->exps;
    n_fq_poly_struct * Tcoeffs = T->coeffs;
    mp_limb_t * v = FLINT_ARRAY_ALLOC(d, mp_limb_t);
    ulong Fexpi, mask;

    FLINT_ASSERT(fq_nmod_mpolyn_is_canonical(F, ctx));
    FLINT_ASSERT(n_fq_bpoly_is_canonical(A, ctx->fqctx));

    mask = (-UWORD(1)) >> (FLINT_BITS - F->bits);
    mpoly_gen_offset_shift_sp(&off0, &shift0, 0, F->bits, ctx->minfo);
    mpoly_gen_offset_shift_sp(&off1, &shift1, 1, F->bits, ctx->minfo);

    FLINT_ASSERT(T->bits == F->bits);

    *lastdeg = -1;

    Ti = Fi = 0;
    Ai = A->length - 1;
    ai = (Ai < 0) ? 0 : n_poly_degree(Acoeffs + Ai);

    while (Fi < Flen || Ai >= 0)
    {
        if (Ti >= T->alloc)
        {
            slong extra = FLINT_MAX(Flen - Fi, Ai);
            fq_nmod_mpolyn_fit_length(T, Ti + extra + 1, ctx);
            Tcoeffs = T->coeffs;
            Texps = T->exps;
        }

        if (Fi < Flen)
            Fexpi = pack_exp2(((Fexps + N*Fi)[off0]>>shift0)&mask,
                              ((Fexps + N*Fi)[off1]>>shift1)&mask);
        else
            Fexpi = 0;

        if (Fi < Flen && Ai >= 0 && Fexpi == pack_exp2(Ai, ai))
        {
            /* F term ok, A term ok */
            mpoly_monomial_set(Texps + N*Ti, Fexps + N*Fi, N);

            n_fq_poly_eval_pow(v, Fcoeffs + Fi, alphapow, ctx->fqctx);
            n_fq_sub(v, Acoeffs[Ai].coeffs + d*ai, v, ctx->fqctx);
            if (!_n_fq_is_zero(v, d))
            {
                changed = 1;
                n_fq_poly_scalar_addmul_n_fq(Tcoeffs + Ti,
                                         Fcoeffs + Fi, modulus, v, ctx->fqctx);
            }
            else
            {
                n_fq_poly_set(Tcoeffs + Ti, Fcoeffs + Fi, ctx->fqctx);
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
        else if (Ai >= 0 && (Fi >= Flen || Fexpi < pack_exp2(Ai, ai)))
        {
            /* F term missing, A term ok */

            mpoly_monomial_zero(Texps + N*Ti, N);
            (Texps + N*Ti)[off0] += (Ai << shift0);
            (Texps + N*Ti)[off1] += (ai << shift1);

            changed = 1;
            n_fq_poly_scalar_mul_n_fq(Tcoeffs + Ti, modulus,
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
            FLINT_ASSERT(Fi < Flen && (Ai < 0 || Fexpi > pack_exp2(Ai, ai)));
            /* F term ok, Aterm missing */
            mpoly_monomial_set(Texps + N*Ti, Fexps + N*Fi, N);

            n_fq_poly_eval_pow(v, Fcoeffs + Fi, alphapow, ctx->fqctx);
            if (!_n_fq_is_zero(v, d))
            {
                changed = 1;
                _n_fq_neg(v, v, d, ctx->fqctx->mod);
                n_fq_poly_scalar_addmul_n_fq(Tcoeffs + Ti, 
                                         Fcoeffs + Fi, modulus, v, ctx->fqctx);
            }
            else
            {
                n_fq_poly_set(Tcoeffs + Ti, Fcoeffs + Fi, ctx->fqctx);                
            }

            Fi++;
        }

        FLINT_ASSERT(!n_fq_poly_is_zero(Tcoeffs + Ti));
        *lastdeg = FLINT_MAX(*lastdeg, n_fq_poly_degree(Tcoeffs + Ti));
        Ti++;
    }

    T->length = Ti;

    if (changed)
        fq_nmod_mpolyn_swap(T, F);

    flint_free(v);

    return changed;
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

void qzip_start(
    n_fq_polyun_t Z,
    n_fq_polyun_t H,
    slong req_images,
    const fq_nmod_ctx_t ctx);

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

static int zip_solvel(
    nmod_mpoly_t A,
    n_polyun_t Z,
    n_polyun_t H,
    n_polyun_t M,
    const nmod_mpoly_ctx_t ctx)
{
    int success;
    slong Ai, i, n;
    mp_limb_t * Acoeffs = A->coeffs;
    n_poly_t t;

    n_poly_init(t);

    FLINT_ASSERT(Z->length == H->length);
    FLINT_ASSERT(Z->length == M->length);

    Ai = 0;
    for (i = 0; i < H->length; i++)
    {
        n = H->terms[i].coeff->length;
        FLINT_ASSERT(M->terms[i].coeff->length == n + 1);
        FLINT_ASSERT(Z->terms[i].coeff->length >= n);
        FLINT_ASSERT(Ai + n <= A->length);

        n_poly_fit_length(t, n);

        success = nmod_zip_find_coeffs_new(Acoeffs + Ai,
                         H->terms[i].coeff->coeffs, n,
                         Z->terms[i].coeff->coeffs, Z->terms[i].coeff->length,
                         M->terms[i].coeff->coeffs, t->coeffs, ctx->ffinfo->mod);
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

int fq_nmod_zip_find_coeffs_new_fq_nmod(
    mp_limb_t * coeffs,             /* length mlength */
    const mp_limb_t * monomials,    /* length mlength */
    slong mlength,
    const mp_limb_t * evals,        /* length d*elength */
    slong elength,
    const mp_limb_t * master,       /* length d*(mlength + 1) */
    mp_limb_t * temp,               /* length d*mlength */
    const fq_nmod_ctx_t ctx);

static int qzip_solvel(
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
/*
flint_printf("i = %wd\n", i);
*/
        n = H->terms[i].coeff->length;
        FLINT_ASSERT(M->terms[i].coeff->length == n + 1);
        FLINT_ASSERT(Z->terms[i].coeff->length >= n);
        FLINT_ASSERT(Ai + n <= A->length);
/*
flint_printf("H[%wd]: ", i);
n_fq_poly_print_pretty(H->terms[i].coeff, "v", ctx->fqctx);
flint_printf("\n");

flint_printf("M[%wd]: ", i);
n_fq_poly_print_pretty(M->terms[i].coeff, "v", ctx->fqctx);
flint_printf("\n");

flint_printf("Z[%wd]: ", i);
n_fq_poly_print_pretty(Z->terms[i].coeff, "v", ctx->fqctx);
flint_printf("\n");
*/

        n_poly_fit_length(t, d*n);

        success = fq_nmod_zip_find_coeffs_new_fq_nmod(A->coeffs + d*Ai,
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
    The cost of interpolating gen(m) into
        G mod (gen(m) - alpha) = sum_i x^e_i*y^f_i c_i(gen(0), ..., gen(m-1))
    is (# = length, deg = deg_(gen(m))):
        (deg G)*(#A + #B + #gamma)                      eval setup
      + (deg G)*(max_i #c_i)*(#A + #B + #gamma)         zip eval
      + (deg G)*(max_i #c_i)*(deg_x AB)^2 (deg_y AB)^2  base gcd
      + (deg G)*(max_i #c_i)*(sum_i # c_i)              zip solve
      + (deg(G))^2*(sum_i #c_i)                         final interp
*/
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
    slong gammadeg)
{
    int use = 0;
    slong lower = FLINT_MAX(gammadeg, rGdeg);
    slong upper = gammadeg + FLINT_MIN(FLINT_MIN(Adeg, Bdeg), rGdeg);
    if (lower <= upper)
    {
        slong Gdeg = ((ulong)upper + (ulong)lower)/2;
        slong Abardeg = gammadeg + Adeg - Gdeg;
        slong Bbardeg = gammadeg + Bdeg - Gdeg;

        if (Gdeg <= Abardeg && Gdeg <= Bbardeg)
            use |= USE_G;

        if (Abardeg <= Gdeg && Abardeg <= Bbardeg)
            use |= USE_ABAR;

        if (Bbardeg <= Gdeg && Bbardeg <= Abardeg)
            use |= USE_BBAR;
    }

    if (use == 0)
        use = USE_G | USE_ABAR | USE_BBAR;

    return use;
}


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
    const fq_nmod_mpolyu_t Bbar);

static int nmod_mpoly_gcd_get_use(
    slong rGdeg,
    slong Adeg,
    slong Bdeg,
    slong gammadeg,
    slong degxAB,
    slong degyAB,
    slong numABgamma,
    const nmod_mpolyu_t G,
    const nmod_mpolyu_t Abar,
    const nmod_mpolyu_t Bbar)
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

static int nmod_mpoly_gcd_get_use_new(
    slong rGdeg,
    slong Adeg,
    slong Bdeg,
    slong gammadeg,
    slong degxAB,
    slong degyAB,
    slong numABgamma,
    const n_polyun_t G,
    const n_polyun_t Abar,
    const n_polyun_t Bbar)
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
            maxnumci = FLINT_MAX(maxnumci, G->terms[i].coeff->length);
            totnumci += G->terms[i].coeff->length;
        }
        FLINT_ASSERT(Gdeg >= 0);
        Gcost = interp_cost(Gdeg,
                       numABgamma, maxnumci, totnumci, degxAB, degyAB);

        maxnumci = totnumci = 0;
        for (i = 0; i < Abar->length; i++)
        {
            maxnumci = FLINT_MAX(maxnumci, Abar->terms[i].coeff->length);
            totnumci += Abar->terms[i].coeff->length;
        }
        FLINT_ASSERT(gammadeg + Adeg - Gdeg >= 0);
        Abarcost = interp_cost(gammadeg + Adeg - Gdeg,
                       numABgamma, maxnumci, totnumci, degxAB, degyAB);

        maxnumci = totnumci = 0;
        for (i = 0; i < Bbar->length; i++)
        {
            maxnumci = FLINT_MAX(maxnumci, Bbar->terms[i].coeff->length);
            totnumci += Bbar->terms[i].coeff->length;
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

slong fq_nmod_mpolyu_total_length(const fq_nmod_mpolyu_t A);

slong nmod_mpolyu_total_length(const nmod_mpolyu_t A)
{
    slong i, tot = 0;
    for (i = 0; i < A->length; i++)
        tot += A->coeffs[i].length;
    return tot;
}

void nmod_mpoly_cvtfrom_mpolyn(
    nmod_mpoly_t A,
    const nmod_mpolyn_t B,
    slong var,
    const nmod_mpoly_ctx_t ctx);

ulong mpoly_bivar_degrees1_chained(ulong deg, ulong * exps, slong len);


void nmod_mpolyn_interp_lift_sm_mpoly(
    nmod_mpolyn_t A,
    const nmod_mpoly_t B,
    const nmod_mpoly_ctx_t ctx);



mp_limb_t _n_poly_mod_eval_step_new(
    mp_limb_t * cur,
    const mp_limb_t * inc,
    const mp_limb_t * coeffs,
    slong length,
    nmod_t ctx)
{
    slong i;
    ulong t0, t1, t2, p0, p1;
    t2 = t1 = t0 = 0;
    for (i = 0; i < length; i++)
    {
        umul_ppmm(p1, p0, cur[i], coeffs[i]);
        add_sssaaaaaa(t2, t1, t0, t2, t1, t0, 0, p1, p0);
        cur[i] = nmod_mul(cur[i], inc[i], ctx);
    }
    NMOD_RED3(t0, t2, t1, t0, ctx);
    return t0;
}

void _n_fq_poly_eval_step_new(
    mp_limb_t * res,
    mp_limb_t * cur,
    const mp_limb_t * inc,
    const mp_limb_t * coeffs,
    slong length,
    const fq_nmod_ctx_t ctx)
{
    slong d = fq_nmod_ctx_degree(ctx);
    slong i;
    mp_limb_t * tmp, * sum;
    TMP_INIT;

    if (length < 1)
    {
        _n_fq_zero(res, d);
        return;
    }

    TMP_START;
    tmp = (mp_limb_t *) TMP_ALLOC(8*d*sizeof(mp_limb_t));
    sum = tmp + 4*d;

    i = 0;
    _n_fq_mul2(sum, cur + d*i, coeffs + d*i, ctx);
    _n_fq_mul(cur + d*i, cur + d*i, inc + d*i, ctx, tmp);
    for (i = 1; i < length; i++)
    {
        _n_fq_madd2(sum, cur + d*i, coeffs + d*i, ctx, tmp);
        _n_fq_mul(cur + d*i, cur + d*i, inc + d*i, ctx, tmp);
    }
    _n_fq_reduce2(res, sum, ctx, tmp);

    TMP_END;
}



mp_limb_t n_poly_mod_eval_step_new(
    n_poly_t cur,
    const n_poly_t inc,
    const nmod_mpoly_t A,
    const nmod_mpoly_ctx_t ctx)
{
    FLINT_ASSERT(A->length == cur->length);
    FLINT_ASSERT(A->length == inc->length);
    return _n_poly_mod_eval_step_new(cur->coeffs, inc->coeffs, A->coeffs,
                                                  A->length, ctx->ffinfo->mod);
}

void n_fq_poly_eval_step_new(
    mp_limb_t * res,
    n_fq_poly_t cur,
    const n_fq_poly_t inc,
    const fq_nmod_mpoly_t A,
    const fq_nmod_ctx_t ctx)
{
    FLINT_ASSERT(A->length == cur->length);
    FLINT_ASSERT(A->length == inc->length);
    _n_fq_poly_eval_step_new(res, cur->coeffs, inc->coeffs, A->coeffs,
                                                               A->length, ctx);
}

void n_bpoly_mod_eval_step_new(
    n_bpoly_t E,
    n_polyun_t cur,
    const n_polyun_t inc,
    const nmod_mpoly_t A,
    const nmod_mpoly_ctx_t ctx)
{
    slong i, Ai;
    slong e0, e1;
    mp_limb_t c;

    n_bpoly_zero(E);

    Ai = 0;
    for (i = 0; i < cur->length; i++)
    {
        slong this_len = cur->terms[i].coeff->length;

        c = _n_poly_mod_eval_step_new(cur->terms[i].coeff->coeffs,
                                  inc->terms[i].coeff->coeffs,
                                  A->coeffs + Ai, this_len, ctx->ffinfo->mod);

        Ai += this_len;

        e0 = extract_exp(cur->terms[i].exp, 1, 2);
        e1 = extract_exp(cur->terms[i].exp, 0, 2);
        if (c == 0)
            continue;

        n_bpoly_set_coeff_nonzero(E, e0, e1, c);
    }

    FLINT_ASSERT(Ai == A->length);
}

void n_fq_bpoly_eval_step_new(
    n_fq_bpoly_t E,
    n_fq_polyun_t cur,
    const n_fq_polyun_t inc,
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

        _n_fq_poly_eval_step_new(c, cur->terms[i].coeff->coeffs,
                                  inc->terms[i].coeff->coeffs,
                                  A->coeffs + d*Ai, this_len, ctx);

        Ai += this_len;

        e0 = extract_exp(cur->terms[i].exp, 1, 2);
        e1 = extract_exp(cur->terms[i].exp, 0, 2);
        if (c == 0)
            continue;

        n_fq_bpoly_set_coeff_n_fq(E, e0, e1, c, ctx);
    }

    FLINT_ASSERT(Ai == A->length);

    flint_free(c);
}



/*
    evaluation at
        gen(start) -> betas[0]
        gen(start+1) -> betas[1]
        ...
        gen(stop-1) -> betas[stop-start-1]

    the other gen are assumed to not appear in A
*/
void nmod_mpoly_get_monomial_evals(
    n_poly_t E,
    const nmod_mpoly_t A,
    const mp_limb_t * betas,
    slong start,
    slong stop,
    const nmod_mpoly_ctx_t ctx)
{
    slong i, Ai;
    flint_bitcnt_t bits = A->bits;
    slong Alen = A->length;
    const ulong * Aexps = A->exps;
    ulong mask = (-UWORD(1)) >> (FLINT_BITS - bits);
    slong N = mpoly_words_per_exp_sp(bits, ctx->minfo);
    slong * off, * shift;
    n_poly_struct * caches;
    mp_limb_t * c;
    slong num = stop - start;

    FLINT_ASSERT(bits <= FLINT_BITS);
    FLINT_ASSERT(Alen > 0);
    FLINT_ASSERT(num > 0);

    caches = FLINT_ARRAY_ALLOC(3*num, n_poly_struct);
    off = FLINT_ARRAY_ALLOC(2*num, slong);
    shift = off + num;
    for (i = 0; i < num; i++)
    {
        mpoly_gen_offset_shift_sp(&off[i], &shift[i], i + start, bits, ctx->minfo);
        n_poly_init(caches + 3*i + 0);
        n_poly_init(caches + 3*i + 1);
        n_poly_init(caches + 3*i + 2);
        nmod_pow_cache_start(betas[i], caches + 3*i + 0,
                                       caches + 3*i + 1,
                                       caches + 3*i + 2);
    }

    n_poly_fit_length(E, Alen);
    E->length = Alen;

    for (Ai = 0; Ai < Alen; Ai++)
    {
        c = E->coeffs + Ai;
        *c = 1;
        for (i = 0; i < num; i++)
        {
            ulong ei = (Aexps[N*Ai + off[i]] >> shift[i]) & mask;
            *c = nmod_pow_cache_mulpow_ui(*c, ei, caches + 3*i + 0,
                                                  caches + 3*i + 1,
                                                  caches + 3*i + 2, ctx->ffinfo->mod);
        }
    }

    for (i = 0; i < num; i++)
    {
        n_poly_clear(caches + 3*i + 0);
        n_poly_clear(caches + 3*i + 1);
        n_poly_clear(caches + 3*i + 2);
    }
    flint_free(caches);
    flint_free(off);
}
    
void fq_nmod_mpoly_get_monomial_evals(
    n_fq_poly_t E,
    const fq_nmod_mpoly_t A,
    const fq_nmod_struct * betas,
    slong start,
    slong stop,
    const fq_nmod_mpoly_ctx_t ctx)
{
    slong d = fq_nmod_ctx_degree(ctx->fqctx);
    slong i, Ai;
    flint_bitcnt_t bits = A->bits;
    slong Alen = A->length;
    const ulong * Aexps = A->exps;
    ulong mask = (-UWORD(1)) >> (FLINT_BITS - bits);
    slong N = mpoly_words_per_exp_sp(bits, ctx->minfo);
    slong * off, * shift;
    n_poly_struct * caches;
    mp_limb_t * c;
    slong num = stop - start;

    FLINT_ASSERT(bits <= FLINT_BITS);
    FLINT_ASSERT(Alen > 0);
    FLINT_ASSERT(num > 0);

    caches = FLINT_ARRAY_ALLOC(3*num, n_fq_poly_struct);
    off = FLINT_ARRAY_ALLOC(2*num, slong);
    shift = off + num;
    for (i = 0; i < num; i++)
    {
        mpoly_gen_offset_shift_sp(&off[i], &shift[i], i + start, bits, ctx->minfo);
        n_fq_poly_init(caches + 3*i + 0);
        n_fq_poly_init(caches + 3*i + 1);
        n_fq_poly_init(caches + 3*i + 2);
        n_fq_pow_cache_start_fq_nmod(betas + i, caches + 3*i + 0,
                                                caches + 3*i + 1,
                                                caches + 3*i + 2, ctx->fqctx);
    }

    n_poly_fit_length(E, d*Alen);
    E->length = Alen;

    for (Ai = 0; Ai < Alen; Ai++)
    {
        c = E->coeffs + d*Ai;
        _n_fq_one(c, d);
        for (i = 0; i < num; i++)
        {
            ulong ei = (Aexps[N*Ai + off[i]] >> shift[i]) & mask;
            n_fq_pow_cache_mulpow_ui(c, c, ei, caches + 3*i + 0,
                                               caches + 3*i + 1,
                                               caches + 3*i + 2, ctx->fqctx);
        }
    }

    for (i = 0; i < num; i++)
    {
        n_fq_poly_clear(caches + 3*i + 0);
        n_fq_poly_clear(caches + 3*i + 1);
        n_fq_poly_clear(caches + 3*i + 2);
    }
    flint_free(caches);
    flint_free(off);
}



/*
    evaluation at
        gen(0) -> x
        gen(1) -> y
        gen(2) -> betas[0]
        gen(3) -> betas[1]
        ...
        gen(m-1) -> betas[m-3]
*/
void nmod_mpoly_get_monomial_evals2(
    n_polyun_t E,
    const nmod_mpoly_t A,
    const mp_limb_t * betas,
    slong m,
    const nmod_mpoly_ctx_t ctx)
{
    slong i, Ai, Ei;
    ulong e0, e1, e01;
    flint_bitcnt_t bits = A->bits;
    slong Alen = A->length;
    const ulong * Aexps = A->exps;
    ulong mask = (-UWORD(1)) >> (FLINT_BITS - bits);
    slong N = mpoly_words_per_exp_sp(bits, ctx->minfo);
    slong * off, * shift;
    n_poly_struct * caches;
    mp_limb_t * c;

    FLINT_ASSERT(bits <= FLINT_BITS);
    FLINT_ASSERT(Alen > 0);
    FLINT_ASSERT(m > 2);

    caches = FLINT_ARRAY_ALLOC(3*(m - 2), n_poly_struct);
    off = FLINT_ARRAY_ALLOC(2*m, slong);
    shift = off + m;
    for (i = 0; i < m; i++)
    {
        mpoly_gen_offset_shift_sp(&off[i], &shift[i], i, bits, ctx->minfo);
        if (i >= 2)
        {
            n_poly_init(caches + 3*(i - 2) + 0);
            n_poly_init(caches + 3*(i - 2) + 1);
            n_poly_init(caches + 3*(i - 2) + 2);
            nmod_pow_cache_start(betas[i - 2], caches + 3*(i - 2) + 0,
                                           caches + 3*(i - 2) + 1,
                                           caches + 3*(i - 2) + 2);
        }
    }

    Ai = 0;
    Ei = 0;

    e0 = (Aexps[N*Ai + off[0]] >> shift[0]) & mask;
    e1 = (Aexps[N*Ai + off[1]] >> shift[1]) & mask;
    e01 = pack_exp2(e0, e1);
    n_polyun_fit_length(E, Ei + 1);
    E->terms[Ei].exp = e01;
    n_poly_fit_length(E->terms[Ei].coeff, 1);
    c = E->terms[Ei].coeff->coeffs + 0;
    E->terms[Ei].coeff->length = 1;
    Ei++;
    *c = 1;
    for (i = 2; i < m; i++)
    {
        ulong ei = (Aexps[N*Ai + off[i]] >> shift[i]) & mask;
        *c = nmod_pow_cache_mulpow_ui(*c, ei, caches + 3*(i - 2) + 0,
                                              caches + 3*(i - 2) + 1,
                                              caches + 3*(i - 2) + 2, ctx->ffinfo->mod);
    }

    for (Ai++; Ai < Alen; Ai++)
    {
        e0 = (Aexps[N*Ai + off[0]] >> shift[0]) & mask;
        e1 = (Aexps[N*Ai + off[1]] >> shift[1]) & mask;
        e01 = pack_exp2(e0, e1);
        if (e01 == E->terms[Ei-1].exp)
        {
            slong len = E->terms[Ei-1].coeff->length;
            n_poly_fit_length(E->terms[Ei-1].coeff, len + 1);
            c = E->terms[Ei-1].coeff->coeffs + len;
            E->terms[Ei-1].coeff->length = len + 1;
        }
        else
        {
            n_polyun_fit_length(E, Ei + 1);
            E->terms[Ei].exp = e01;
            n_poly_fit_length(E->terms[Ei].coeff, 1);
            c = E->terms[Ei].coeff->coeffs + 0;
            E->terms[Ei].coeff->length = 1;
            Ei++;
        }

        *c = 1;
        for (i = 2; i < m; i++)
        {
            ulong ei = (Aexps[N*Ai + off[i]] >> shift[i]) & mask;
            *c = nmod_pow_cache_mulpow_ui(*c, ei, caches + 3*(i - 2) + 0,
                                                  caches + 3*(i - 2) + 1,
                                                  caches + 3*(i - 2) + 2, ctx->ffinfo->mod);
        }
    }

    E->length = Ei;

    for (i = 0; i < m - 2; i++)
    {
        n_poly_clear(caches + 3*i + 0);
        n_poly_clear(caches + 3*i + 1);
        n_poly_clear(caches + 3*i + 2);
    }
    flint_free(caches);
    flint_free(off);

#if FLINT_WANT_ASSERT
    Ai = 0;
    for (i = 0; i < E->length; i++)
        Ai += E->terms[i].coeff->length;
    FLINT_ASSERT(Ai == A->length);
#endif
}

void fq_nmod_mpoly_get_monomial_evals2(
    n_fq_polyun_t E,
    const fq_nmod_mpoly_t A,
    const fq_nmod_struct * betas,
    slong m,
    const fq_nmod_mpoly_ctx_t ctx)
{
    slong d = fq_nmod_ctx_degree(ctx->fqctx);
    slong i, Ai, Ei;
    ulong e0, e1, e01;
    flint_bitcnt_t bits = A->bits;
    slong Alen = A->length;
    const ulong * Aexps = A->exps;
    ulong mask = (-UWORD(1)) >> (FLINT_BITS - bits);
    slong N = mpoly_words_per_exp_sp(bits, ctx->minfo);
    slong * off, * shift;
    n_fq_poly_struct * caches;
    mp_limb_t * c;

    FLINT_ASSERT(bits <= FLINT_BITS);
    FLINT_ASSERT(Alen > 0);
    FLINT_ASSERT(m > 2);

    caches = FLINT_ARRAY_ALLOC(3*(m - 2), n_fq_poly_struct);
    off = FLINT_ARRAY_ALLOC(2*m, slong);
    shift = off + m;
    for (i = 0; i < m; i++)
    {
        mpoly_gen_offset_shift_sp(&off[i], &shift[i], i, bits, ctx->minfo);
        if (i >= 2)
        {
            n_fq_poly_init(caches + 3*(i - 2) + 0);
            n_fq_poly_init(caches + 3*(i - 2) + 1);
            n_fq_poly_init(caches + 3*(i - 2) + 2);
            n_fq_pow_cache_start_fq_nmod(betas + i - 2, caches + 3*(i - 2) + 0,
                                           caches + 3*(i - 2) + 1,
                                           caches + 3*(i - 2) + 2, ctx->fqctx);
        }
    }

    Ai = 0;
    Ei = 0;

    e0 = (Aexps[N*Ai + off[0]] >> shift[0]) & mask;
    e1 = (Aexps[N*Ai + off[1]] >> shift[1]) & mask;
    e01 = pack_exp2(e0, e1);
    n_polyun_fit_length(E, Ei + 1);
    E->terms[Ei].exp = e01;
    n_poly_fit_length(E->terms[Ei].coeff, d*1);
    c = E->terms[Ei].coeff->coeffs + d*0;
    E->terms[Ei].coeff->length = 1;
    Ei++;
    _n_fq_one(c, d);
    for (i = 2; i < m; i++)
    {
        ulong ei = (Aexps[N*Ai + off[i]] >> shift[i]) & mask;
        n_fq_pow_cache_mulpow_ui(c, c, ei, caches + 3*(i - 2) + 0,
                                           caches + 3*(i - 2) + 1,
                                           caches + 3*(i - 2) + 2, ctx->fqctx);
    }

    for (Ai++; Ai < Alen; Ai++)
    {
        e0 = (Aexps[N*Ai + off[0]] >> shift[0]) & mask;
        e1 = (Aexps[N*Ai + off[1]] >> shift[1]) & mask;
        e01 = pack_exp2(e0, e1);
        if (e01 == E->terms[Ei-1].exp)
        {
            slong len = E->terms[Ei-1].coeff->length;
            n_poly_fit_length(E->terms[Ei-1].coeff, d*(len + 1));
            c = E->terms[Ei-1].coeff->coeffs + d*len;
            E->terms[Ei-1].coeff->length = len + 1;
        }
        else
        {
            n_polyun_fit_length(E, Ei + 1);
            E->terms[Ei].exp = e01;
            n_poly_fit_length(E->terms[Ei].coeff, d*1);
            c = E->terms[Ei].coeff->coeffs + d*0;
            E->terms[Ei].coeff->length = 1;
            Ei++;
        }

        _n_fq_one(c, d);
        for (i = 2; i < m; i++)
        {
            ulong ei = (Aexps[N*Ai + off[i]] >> shift[i]) & mask;
            n_fq_pow_cache_mulpow_ui(c, c, ei, caches + 3*(i - 2) + 0,
                                               caches + 3*(i - 2) + 1,
                                               caches + 3*(i - 2) + 2, ctx->fqctx);
        }
    }

    E->length = Ei;

    for (i = 0; i < m - 2; i++)
    {
        n_fq_poly_clear(caches + 3*i + 0);
        n_fq_poly_clear(caches + 3*i + 1);
        n_fq_poly_clear(caches + 3*i + 2);
    }
    flint_free(caches);
    flint_free(off);

#if FLINT_WANT_ASSERT
    Ai = 0;
    for (i = 0; i < E->length; i++)
        Ai += E->terms[i].coeff->length;
    FLINT_ASSERT(Ai == A->length);
#endif
}



/*
    evaluation at
        gen(0) -> x
        gen(1) -> y
        gen(2) -> betas[0]
        gen(3) -> betas[1]
        ...
        gen(m-1) -> betas[m - 3]
*/
slong nmod_mpoly_set_zip_form2(
    n_polyun_t H,
    n_polyun_t M,
    const nmod_mpoly_t A,
    const mp_limb_t * betas,
    slong m,
    const nmod_mpoly_ctx_t ctx)
{
    slong i, zip_length;

    nmod_mpoly_get_monomial_evals2(H, A, betas, m, ctx);

    n_polyun_fit_length(M, H->length);
    M->length = H->length;
    zip_length = 0;
    for (i = 0; i < H->length; i++)
    {
        slong len = H->terms[i].coeff->length;
        M->terms[i].exp = H->terms[i].exp;
        zip_length = FLINT_MAX(zip_length, len);
        n_poly_mod_product_roots_nmod_vec(M->terms[i].coeff,
                            H->terms[i].coeff->coeffs, len, ctx->ffinfo->mod);
    }

    return zip_length;
}


slong fq_nmod_mpoly_set_zip_form2(
    n_fq_polyun_t H,
    n_fq_polyun_t M,
    const fq_nmod_mpoly_t A,
    const fq_nmod_struct * betas,
    slong m,
    const fq_nmod_mpoly_ctx_t ctx)
{
    slong i, zip_length;

    fq_nmod_mpoly_get_monomial_evals2(H, A, betas, m, ctx);

    n_polyun_fit_length(M, H->length);
    M->length = H->length;
    zip_length = 0;
    for (i = 0; i < H->length; i++)
    {
        slong len = H->terms[i].coeff->length;
        M->terms[i].exp = H->terms[i].exp;
        zip_length = FLINT_MAX(zip_length, len);
        n_fq_poly_product_roots_n_fq(M->terms[i].coeff,
                                   H->terms[i].coeff->coeffs, len, ctx->fqctx);
    }

    return zip_length;
}



/*
    F = F + modulus*(A - F(alpha))
    monomials assumed to match
*/
int nmod_mpolyn_interp_mcrt_sm_mpoly(
    slong * lastdeg_,
    nmod_mpolyn_t F,
    const nmod_mpoly_t A,
    const n_poly_t modulus,
    n_poly_t alphapow,
    const nmod_mpoly_ctx_t ctx)
{
    slong lastdeg = -1;
    int changed = 0;
    slong i;
    mp_limb_t v;
    mp_limb_t * Acoeff = A->coeffs;
    slong Flen = F->length;

    FLINT_ASSERT(Flen == A->length);

    for (i = 0; i < Flen; i++)
    {
        /* F term ok, A term ok */
        v = n_poly_mod_eval_pow((n_poly_struct*)(F->coeffs + i), alphapow, ctx->ffinfo->mod);
        v = nmod_sub(Acoeff[i], v, ctx->ffinfo->mod);
        if (v != 0)
        {
            changed = 1;
            n_poly_mod_scalar_addmul_nmod((n_poly_struct*)(F->coeffs + i),
                (n_poly_struct*)(F->coeffs + i), modulus, v, ctx->ffinfo->mod);
        }
        lastdeg = FLINT_MAX(lastdeg, n_poly_degree((n_poly_struct*)(F->coeffs + i)));
    }

    *lastdeg_ = lastdeg;
    return changed;
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

    The degrees of A, B, and gamma wrt the minor vars must be passed in.
    A guess of the degrees of rG wrt the minor vars can be passed in.


    deg(G) = deg(gamma) - deg(lc(rG)) +  deg(rG)
    deg(G) <= deg(gamma) + deg(rG)
    deg(G) <= deg(gamma) + deg(A)
    deg(G) <= deg(gamma) + deg(B)
    deg(G) >= deg(gamma)
    deg(G) >= deg(rG)

    deg(A) = deg(gamma) + deg(A) - deg(G)
    deg(B) = deg(gamma) + deg(B) - deg(G)
*/
int nmod_mpolyl_gcd_zippel_smprime(
    nmod_mpoly_t rG, const slong * rGdegs, /* guess at rG degrees, could be NULL */
    nmod_mpoly_t rAbar,
    nmod_mpoly_t rBbar,
    const nmod_mpoly_t A, const slong * Adegs,
    const nmod_mpoly_t B, const slong * Bdegs,
    const nmod_mpoly_t gamma, const slong * gammadegs,
    const nmod_mpoly_ctx_t ctx)
{
    int success, use, alpha_tries_left;
    nmod_t mod = ctx->ffinfo->mod;
    slong i, j, m;
    slong nvars = ctx->minfo->nvars;
    flint_bitcnt_t bits = A->bits;
    mp_limb_t * alphas, * betas;
    flint_rand_t state;
    nmod_mpoly_t cont;
    nmod_mpoly_t T, G, Abar, Bbar;
    n_polyun_t HG, HAbar, HBbar, MG, MAbar, MBbar, ZG, ZAbar, ZBbar;
    n_bpoly_t Aev, Bev, Gev, Abarev, Bbarev;
    mp_limb_t gammaev;
    nmod_mpolyn_t Tn, Gn, Abarn, Bbarn;
    slong lastdeg;
    slong cur_zip_image, req_zip_images, this_length;
    n_polyun_t Aeh_cur, Aeh_inc, Beh_cur, Beh_inc;
    n_poly_t gammaeh_cur, gammaeh_inc;
    n_poly_t alphapow;
    nmod_mpoly_struct * Aevals, * Bevals;
    nmod_mpoly_struct * gammaevals;
    n_poly_t modulus;
    n_poly_bpoly_stack_t St;
    mp_limb_t c, start_alpha;
    ulong GdegboundXY, newdegXY, Abideg, Bbideg;
    slong degxAB, degyAB;
/*
flint_printf("nmod_mpolyl_gcd_zippel_smprime called mod.n = %wu\n", ctx->ffinfo->mod.n);
*/

    FLINT_ASSERT(bits <= FLINT_BITS);
    FLINT_ASSERT(bits == rG->bits);
    FLINT_ASSERT(bits == rAbar->bits);
    FLINT_ASSERT(bits == rBbar->bits);
    FLINT_ASSERT(bits == A->bits);
    FLINT_ASSERT(bits == B->bits);
    FLINT_ASSERT(bits == gamma->bits);
    FLINT_ASSERT(A->length > 0);
    FLINT_ASSERT(B->length > 0);
    FLINT_ASSERT(gamma->length > 0);

#if FLINT_WANT_ASSERT
    {
        slong * tmp_degs = FLINT_ARRAY_ALLOC(nvars, slong);

        nmod_mpoly_degrees_si(tmp_degs, A, ctx);
        for (j = 0; j < nvars; j++)
            FLINT_ASSERT(tmp_degs[j] == Adegs[j]);

        nmod_mpoly_degrees_si(tmp_degs, B, ctx);
        for (j = 0; j < nvars; j++)
            FLINT_ASSERT(tmp_degs[j] == Bdegs[j]);

        nmod_mpoly_degrees_si(tmp_degs, gamma, ctx);
        for (j = 0; j < nvars; j++)
            FLINT_ASSERT(tmp_degs[j] == gammadegs[j]);

        flint_free(tmp_degs);
    }
#endif

    if (ctx->ffinfo->mod.n < 7)
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
    nmod_mpoly_init3(T, 1, bits, ctx);
    nmod_mpoly_init3(G, 1, bits, ctx);
    nmod_mpoly_init3(Abar, 1, bits, ctx);
    nmod_mpoly_init3(Bbar, 1, bits, ctx);
    nmod_mpolyn_init(Tn, bits, ctx);
    nmod_mpolyn_init(Gn, bits, ctx);
    nmod_mpolyn_init(Abarn, bits, ctx);
    nmod_mpolyn_init(Bbarn, bits, ctx);
    n_polyun_init(Aeh_cur);
    n_polyun_init(Aeh_inc);
    n_polyun_init(Beh_cur);
    n_polyun_init(Beh_inc);
    n_poly_init(gammaeh_cur);
    n_poly_init(gammaeh_inc);
    n_poly_init(modulus);
    n_poly_stack_init(St->poly_stack);
    n_bpoly_stack_init(St->bpoly_stack);

    betas = FLINT_ARRAY_ALLOC(nvars, mp_limb_t);
    alphas = FLINT_ARRAY_ALLOC(nvars, mp_limb_t);
    flint_randinit(state);

    Aevals = FLINT_ARRAY_ALLOC(nvars + 1, nmod_mpoly_struct);
    Bevals = FLINT_ARRAY_ALLOC(nvars + 1, nmod_mpoly_struct);
    gammaevals = FLINT_ARRAY_ALLOC(nvars + 1, nmod_mpoly_struct);
    for (i = 0; i < nvars; i++)
    {
        nmod_mpoly_init3(Aevals + i, 0, bits, ctx);
        nmod_mpoly_init3(Bevals + i, 0, bits, ctx);
        nmod_mpoly_init3(gammaevals + i, 0, bits, ctx);
    }
    Aevals[nvars] = *A;
    Bevals[nvars] = *B;
    gammaevals[nvars] = *gamma;

    Abideg = _nmod_mpoly_bidegree(A, ctx);
    Bbideg = _nmod_mpoly_bidegree(B, ctx);

    degxAB = FLINT_MAX(Adegs[0], Bdegs[0]);
    degyAB = FLINT_MAX(Adegs[1], Bdegs[1]);

    GdegboundXY = pack_exp2(FLINT_MIN(Adegs[0], Bdegs[0]),
                            FLINT_MIN(Adegs[1], Bdegs[1]));
    if (GdegboundXY == 0)
        goto gcd_is_trivial;

    alpha_tries_left = 20;

choose_alphas:

    if (--alpha_tries_left < 0)
    {
        success = 0;
        goto cleanup;
    }

    for (i = 2; i < nvars; i++)
        alphas[i] = n_urandint(state, ctx->ffinfo->mod.n - 2) + 1;

    for (i = nvars - 1; i >= 2; i--)
    {
        nmod_mpoly_evaluate_one_ui(Aevals + i, Aevals + i + 1, i, alphas[i], ctx);
        nmod_mpoly_repack_bits_inplace(Aevals + i, bits, ctx);
        nmod_mpoly_evaluate_one_ui(Bevals + i, Bevals + i + 1, i, alphas[i], ctx);
        nmod_mpoly_repack_bits_inplace(Bevals + i, bits, ctx);
        nmod_mpoly_evaluate_one_ui(gammaevals + i, gammaevals + i + 1, i, alphas[i], ctx);
        nmod_mpoly_repack_bits_inplace(gammaevals + i, bits, ctx);
        if (nmod_mpoly_is_zero(gammaevals + i, ctx))
            goto choose_alphas;
        if (Aevals[i].length < 1 || _nmod_mpoly_bidegree(Aevals + i, ctx) != Abideg)
            goto choose_alphas;
        if (Bevals[i].length < 1 || _nmod_mpoly_bidegree(Bevals + i, ctx) != Bbideg)
            goto choose_alphas;
    }

    m = 2;

    if (rGdegs == NULL)
        use = USE_G | USE_ABAR | USE_BBAR;
    else
        use = mpoly_gcd_get_use_first(rGdegs[m], Adegs[m], Bdegs[m], gammadegs[m]);

    nmod_mpoly_get_bpoly(Aev, Aevals + m, 0, 1, ctx);
    nmod_mpoly_get_bpoly(Bev, Bevals + m, 0, 1, ctx);

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

    nmod_mpolyn_interp_lift_sm_bpoly(Gn, Gev, ctx);
    nmod_mpolyn_interp_lift_sm_bpoly(Abarn, Abarev, ctx);
    nmod_mpolyn_interp_lift_sm_bpoly(Bbarn, Bbarev, ctx);

    n_poly_one(modulus);
    c = nmod_neg(alphas[m], mod);
    n_poly_mod_shift_left_scalar_addmul(modulus, 1, c, mod);

    start_alpha = alphas[m];
    while (1)
    {
    choose_alpha_2:

        alphas[m] = (alphas[m] < 2) ? mod.n - 1 : alphas[m] - 1;

        if (alphas[m] == start_alpha)
            goto choose_alphas;

        FLINT_ASSERT(alphapow->alloc >= 2);
        alphapow->coeffs[0] = 1;
        alphapow->coeffs[1] = alphas[m];
        alphapow->length = 2;

        nmod_mpoly_evaluate_one_ui(Aevals + m, Aevals + m + 1, m, alphas[m], ctx);
        nmod_mpoly_repack_bits_inplace(Aevals + m, bits, ctx);
        nmod_mpoly_evaluate_one_ui(Bevals + m, Bevals + m + 1, m, alphas[m], ctx);
        nmod_mpoly_repack_bits_inplace(Bevals + m, bits, ctx);
        nmod_mpoly_evaluate_one_ui(gammaevals + m, gammaevals + m + 1, m, alphas[m], ctx);
        nmod_mpoly_repack_bits_inplace(gammaevals + m, bits, ctx);

        if (nmod_mpoly_is_zero(gammaevals + m, ctx))
            goto choose_alpha_2;
        if (Aevals[m].length < 1 || _nmod_mpoly_bidegree(Aevals + m, ctx) != Abideg)
            goto choose_alpha_2;
        if (Bevals[m].length < 1 || _nmod_mpoly_bidegree(Bevals + m, ctx) != Bbideg)
            goto choose_alpha_2;

        nmod_mpoly_get_bpoly(Aev, Aevals + m, 0, 1, ctx);
        nmod_mpoly_get_bpoly(Bev, Bevals + m, 0, 1, ctx);

        success = n_bpoly_mod_gcd_brown_smprime(Gev, Abarev, Bbarev, Aev, Bev, mod, St);
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
            goto choose_alphas;
        }

        gammaev = nmod_mpoly_get_ui(gammaevals + m, ctx);
        n_bpoly_scalar_mul_nmod(Gev, gammaev, mod);

        c = n_poly_mod_eval_pow(modulus, alphapow, mod);
        c = nmod_inv(c, mod);
        _n_poly_mod_scalar_mul_nmod(modulus, modulus, c, mod);

        if ((use & USE_G) && !nmod_mpolyn_interp_crt_sm_bpoly(
                                &lastdeg, Gn, Tn, Gev, modulus, alphapow, ctx))
        {
            if (m == nvars - 1)
            {
                nmod_mpoly_cvtfrom_mpolyn(rG, Gn, m, ctx);
                success = nmod_mpolyl_content(cont, rG, 2, ctx);
                if (!success)
                    goto cleanup;
                nmod_mpoly_divides(rG, rG, cont, ctx);
                nmod_mpoly_repack_bits_inplace(rG, bits, ctx);
                if (nmod_mpoly_divides(rAbar, A, rG, ctx) &&
                    nmod_mpoly_divides(rBbar, B, rG, ctx))
                {
                    nmod_mpoly_repack_bits_inplace(rAbar, bits, ctx);
                    nmod_mpoly_repack_bits_inplace(rBbar, bits, ctx);
                    break;
                }                
            }
            else
            {
                nmod_mpoly_cvtfrom_mpolyn(G, Gn, m, ctx);
                nmod_mpoly_mul(T, Aevals + m + 1, gammaevals + m + 1, ctx);
                if (nmod_mpoly_divides(Abar, T, G, ctx))
                {
                    nmod_mpoly_repack_bits_inplace(Abar, bits, ctx);
                    nmod_mpoly_mul(T, Bevals + m + 1, gammaevals + m + 1, ctx);
                    if (nmod_mpoly_divides(Bbar, T, G, ctx))
                    {
                        nmod_mpoly_repack_bits_inplace(Bbar, bits, ctx);
                        break;
                    }
                }
            }
        }

        if ((use & USE_ABAR) && !nmod_mpolyn_interp_crt_sm_bpoly(
                          &lastdeg, Abarn, Tn, Abarev, modulus, alphapow, ctx))
        {
            if (m == nvars - 1)
            {
                nmod_mpoly_cvtfrom_mpolyn(rAbar, Abarn, m, ctx);
                success = nmod_mpolyl_content(cont, rAbar, 2, ctx);
                if (!success)
                    goto cleanup;
                nmod_mpoly_divides(rAbar, rAbar, cont, ctx);
                nmod_mpoly_repack_bits_inplace(rAbar, bits, ctx);
                if (nmod_mpoly_divides(rG, A, rAbar, ctx) &&
                    nmod_mpoly_divides(rBbar, B, rG, ctx))
                {
                    nmod_mpoly_repack_bits_inplace(rG, bits, ctx);
                    nmod_mpoly_repack_bits_inplace(rBbar, bits, ctx);
                    break;
                }                
            }
            else
            {
                nmod_mpoly_cvtfrom_mpolyn(Abar, Abarn, m, ctx);
                nmod_mpoly_mul(T, Aevals + m + 1, gammaevals + m + 1, ctx);
                if (nmod_mpoly_divides(G, T, Abar, ctx))
                {
                    nmod_mpoly_repack_bits_inplace(G, bits, ctx);
                    nmod_mpoly_mul(T, Bevals + m + 1, gammaevals + m + 1, ctx);
                    if (nmod_mpoly_divides(Bbar, T, G, ctx))
                    {
                        nmod_mpoly_repack_bits_inplace(Bbar, bits, ctx);
                        break;
                    }
                }
            }
        }

        if ((use & USE_BBAR) && !nmod_mpolyn_interp_crt_sm_bpoly(
                          &lastdeg, Bbarn, Tn, Bbarev, modulus, alphapow, ctx))
        {
            if (m == nvars - 1)
            {
                nmod_mpoly_cvtfrom_mpolyn(rBbar, Bbarn, m, ctx);
                success = nmod_mpolyl_content(cont, rBbar, 2, ctx);
                if (!success)
                    goto cleanup;
                nmod_mpoly_divides(rBbar, rBbar, cont, ctx);
                nmod_mpoly_repack_bits_inplace(rBbar, bits, ctx);
                if (nmod_mpoly_divides(rG, B, rBbar, ctx) &&
                    nmod_mpoly_divides(rAbar, A, rG, ctx))
                {
                    nmod_mpoly_repack_bits_inplace(rG, bits, ctx);
                    nmod_mpoly_repack_bits_inplace(rAbar, bits, ctx);
                    break;
                }                
            }
            else
            {
                nmod_mpoly_cvtfrom_mpolyn(Bbar, Bbarn, m, ctx);
                nmod_mpoly_mul(T, Bevals + m + 1, gammaevals + m + 1, ctx);
                if (nmod_mpoly_divides(G, T, Bbar, ctx))
                {
                    nmod_mpoly_repack_bits_inplace(G, bits, ctx);
                    nmod_mpoly_mul(T, Aevals + m + 1, gammaevals + m + 1, ctx);
                    if (nmod_mpoly_divides(Abar, T, G, ctx))
                    {
                        nmod_mpoly_repack_bits_inplace(Abar, bits, ctx);
                        break;
                    }
                }
            }
        }

        if (n_poly_degree(modulus) > gammadegs[m] + Adegs[m] &&
            n_poly_degree(modulus) > gammadegs[m] + Bdegs[m])
        {
            goto choose_alphas;
        }

        c = nmod_neg(alphas[m], mod);
        n_poly_mod_shift_left_scalar_addmul(modulus, 1, c, mod);
    }

    for (m = 3; m < nvars; m++)
    {
        /* G, Abar, Bbar are in Fp[gen(0), ..., gen(m - 1)] */
        nmod_mpolyn_interp_lift_sm_mpoly(Gn, G, ctx);
        nmod_mpolyn_interp_lift_sm_mpoly(Abarn, Abar, ctx);
        nmod_mpolyn_interp_lift_sm_mpoly(Bbarn, Bbar, ctx);

        if (rGdegs == NULL)
            use = USE_G | USE_ABAR | USE_BBAR;
        else
            use = 0;

        n_poly_one(modulus);
        c = nmod_neg(alphas[m], mod);
        n_poly_mod_shift_left_scalar_addmul(modulus, 1, c, mod);

        start_alpha = alphas[m];

    choose_betas:

        /* only beta[2], beta[1], ..., beta[m - 1] will be used */
        for (i = 2; i < ctx->minfo->nvars; i++)
            betas[i] = n_urandint(state, ctx->ffinfo->mod.n - 3) + 2;

        req_zip_images = 1;

        /* TODO compute H first and then compute M only if needed */
        this_length = nmod_mpoly_set_zip_form2(HAbar, MAbar, Abar, betas + 2, m, ctx);
        req_zip_images = FLINT_MAX(req_zip_images, this_length);

        this_length = nmod_mpoly_set_zip_form2(HBbar, MBbar, Bbar, betas + 2, m, ctx);
        req_zip_images = FLINT_MAX(req_zip_images, this_length);

        this_length = nmod_mpoly_set_zip_form2(HG, MG, G, betas + 2, m, ctx);
        req_zip_images = FLINT_MAX(req_zip_images, this_length);

        if (use == 0)
        {
            this_length = gammaevals[m + 1].length + Aevals[m + 1].length +
                                                     Bevals[m + 1].length;

            use = nmod_mpoly_gcd_get_use_new(rGdegs[m], Adegs[m], Bdegs[m],
                  gammadegs[m], degxAB, degyAB, this_length, HG, HAbar, HBbar);
        }

        while (1)
        {
        choose_alpha_m:

            alphas[m] = (alphas[m] < 2) ? mod.n - 1 : alphas[m] - 1;
            if (alphas[m] == start_alpha)
                goto choose_alphas;

            FLINT_ASSERT(alphapow->alloc >= 2);
            alphapow->coeffs[0] = 1;
            alphapow->coeffs[1] = alphas[m];
            alphapow->length = 2;

            nmod_mpoly_evaluate_one_ui(Aevals + m, Aevals + m + 1, m, alphas[m], ctx);
            nmod_mpoly_evaluate_one_ui(Bevals + m, Bevals + m + 1, m, alphas[m], ctx);
            nmod_mpoly_evaluate_one_ui(gammaevals + m, gammaevals + m + 1, m, alphas[m], ctx);
            if (nmod_mpoly_is_zero(gammaevals + m, ctx))
                goto choose_alpha_m;
            if (Aevals[m].length < 1 || _nmod_mpoly_bidegree(Aevals + m, ctx) != Abideg)
                goto choose_alpha_m;
            if (Bevals[m].length < 1 || _nmod_mpoly_bidegree(Bevals + m, ctx) != Bbideg)
                goto choose_alpha_m;

            nmod_mpoly_get_monomial_evals2(Aeh_inc, Aevals + m, betas + 2, m, ctx);
            nmod_mpoly_get_monomial_evals2(Beh_inc, Bevals + m, betas + 2, m, ctx);
            nmod_mpoly_get_monomial_evals(gammaeh_inc, gammaevals + m, betas + 2, 2, m, ctx);

            n_polyun_set(Aeh_cur, Aeh_inc);
            n_polyun_set(Beh_cur, Beh_inc);
            n_poly_set(gammaeh_cur, gammaeh_inc);

            zip_start(ZG, HG, req_zip_images);
            zip_start(ZAbar, HAbar, req_zip_images);
            zip_start(ZBbar, HBbar, req_zip_images);
            for (cur_zip_image = 0; cur_zip_image < req_zip_images; cur_zip_image++)
            {
                n_bpoly_mod_eval_step_new(Aev, Aeh_cur, Aeh_inc, Aevals + m, ctx);
                n_bpoly_mod_eval_step_new(Bev, Beh_cur, Beh_inc, Bevals + m, ctx);
                gammaev = n_poly_mod_eval_step_new(gammaeh_cur, gammaeh_inc, gammaevals + m, ctx);
                if (gammaev == 0)
                    goto choose_betas;
                if (Aev->length < 1 || n_bpoly_bidegree(Aev) != Abideg)
                    goto choose_betas;
                if (Bev->length < 1 || n_bpoly_bidegree(Bev) != Bbideg)
                    goto choose_betas;

                success = n_bpoly_mod_gcd_brown_smprime(Gev, Abarev, Bbarev,
                                                            Aev, Bev, mod, St);        
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

            if ((use & USE_G) && zip_solvel(G, ZG, HG, MG, ctx) < 1)
                goto choose_alphas;
            if ((use & USE_ABAR) && zip_solvel(Abar, ZAbar, HAbar, MAbar, ctx) < 1)
                goto choose_alphas;
            if ((use & USE_BBAR) && zip_solvel(Bbar, ZBbar, HBbar, MBbar, ctx) < 1)
                goto choose_alphas;

            FLINT_ASSERT(n_poly_degree(modulus) > 0);
            c = n_poly_mod_eval_pow(modulus, alphapow, mod);
            c = nmod_inv(c, mod);
            _n_poly_mod_scalar_mul_nmod(modulus, modulus, c, mod);

            if ((use & USE_G) && !nmod_mpolyn_interp_mcrt_sm_mpoly(
                                      &lastdeg, Gn, G, modulus, alphapow, ctx))
            {
                nmod_mpoly_cvtfrom_mpolyn(rG, Gn, m, ctx);
                if (m == nvars - 1)
                {
                    success = nmod_mpolyl_content(cont, rG, 2, ctx);
                    if (!success)
                        goto cleanup;
                    nmod_mpoly_divides(rG, rG, cont, ctx);
                    nmod_mpoly_repack_bits_inplace(rG, bits, ctx);
                    if (nmod_mpoly_divides(rAbar, A, rG, ctx) &&
                        nmod_mpoly_divides(rBbar, B, rG, ctx))
                    {
                        nmod_mpoly_repack_bits_inplace(rAbar, bits, ctx);
                        nmod_mpoly_repack_bits_inplace(rBbar, bits, ctx);
                        break;
                    }
                }
                else
                {
                    nmod_mpoly_mul(T, Aevals + m + 1, gammaevals + m + 1, ctx);
                    if (nmod_mpoly_divides(rAbar, T, rG, ctx))
                    {
                        nmod_mpoly_mul(T, Bevals + m + 1, gammaevals + m + 1, ctx);
                        if (nmod_mpoly_divides(rBbar, T, rG, ctx))
                        {
                            nmod_mpoly_repack_bits_inplace(rAbar, bits, ctx);
                            nmod_mpoly_repack_bits_inplace(rBbar, bits, ctx);
                            nmod_mpoly_swap(G, rG, ctx);
                            nmod_mpoly_swap(Abar, rAbar, ctx);
                            nmod_mpoly_swap(Bbar, rBbar, ctx);
                            break;
                        }
                    }
                }
            }
            if ((use & USE_ABAR) && !nmod_mpolyn_interp_mcrt_sm_mpoly(
                                &lastdeg, Abarn, Abar, modulus, alphapow, ctx))
            {
                nmod_mpoly_cvtfrom_mpolyn(rAbar, Abarn, m, ctx);
                if (m == nvars - 1)
                {
                    success = nmod_mpolyl_content(cont, rAbar, 2, ctx);
                    if (!success)
                        goto cleanup;
                    nmod_mpoly_divides(rAbar, rAbar, cont, ctx);
                    nmod_mpoly_repack_bits_inplace(rAbar, bits, ctx);
                    if (nmod_mpoly_divides(rG, A, rAbar, ctx) &&
                        nmod_mpoly_divides(rBbar, B, rG, ctx))
                    {
                        nmod_mpoly_repack_bits_inplace(rG, bits, ctx);
                        nmod_mpoly_repack_bits_inplace(rBbar, bits, ctx);
                        break;
                    }
                }
                else
                {
                    nmod_mpoly_mul(T, Aevals + m + 1, gammaevals + m + 1, ctx);
                    if (nmod_mpoly_divides(rG, T, rAbar, ctx))
                    {
                        nmod_mpoly_mul(T, Bevals + m + 1, gammaevals + m + 1, ctx);
                        if (nmod_mpoly_divides(rBbar, T, rG, ctx))
                        {
                            nmod_mpoly_repack_bits_inplace(rG, bits, ctx);
                            nmod_mpoly_repack_bits_inplace(rBbar, bits, ctx);
                            nmod_mpoly_swap(G, rG, ctx);
                            nmod_mpoly_swap(Abar, rAbar, ctx);
                            nmod_mpoly_swap(Bbar, rBbar, ctx);
                            break;
                        }
                    }
                }
            }
            if ((use & USE_BBAR) && !nmod_mpolyn_interp_mcrt_sm_mpoly(
                                &lastdeg, Bbarn, Bbar, modulus, alphapow, ctx))
            {
                nmod_mpoly_cvtfrom_mpolyn(rBbar, Bbarn, m, ctx);
                if (m == nvars - 1)
                {
                    success = nmod_mpolyl_content(cont, rBbar, 2, ctx);
                    if (!success)
                        goto cleanup;
                    nmod_mpoly_divides(rBbar, rBbar, cont, ctx);
                    nmod_mpoly_repack_bits_inplace(rBbar, bits, ctx);
                    if (nmod_mpoly_divides(rG, B, rBbar, ctx) &&
                        nmod_mpoly_divides(rAbar, A, rG, ctx))
                    {
                        nmod_mpoly_repack_bits_inplace(rG, bits, ctx);
                        nmod_mpoly_repack_bits_inplace(rAbar, bits, ctx);
                        break;
                    }
                }
                else
                {
                    nmod_mpoly_mul(T, Bevals + m + 1, gammaevals + m + 1, ctx);
                    if (nmod_mpoly_divides(rG, T, rBbar, ctx))
                    {
                        nmod_mpoly_mul(T, Aevals + m + 1, gammaevals + m + 1, ctx);
                        if (nmod_mpoly_divides(rAbar, T, rG, ctx))
                        {
                            nmod_mpoly_repack_bits_inplace(rG, bits, ctx);
                            nmod_mpoly_repack_bits_inplace(rAbar, bits, ctx);
                            nmod_mpoly_swap(G, rG, ctx);
                            nmod_mpoly_swap(Abar, rAbar, ctx);
                            nmod_mpoly_swap(Bbar, rBbar, ctx);
                            break;
                        }
                    }
                }
            }

            if (n_poly_degree(modulus) > gammadegs[m] + Adegs[m] &&
                n_poly_degree(modulus) > gammadegs[m] + Bdegs[m])
            {
                goto choose_alphas;
            }

            c = nmod_neg(alphas[m], mod);
            n_poly_mod_shift_left_scalar_addmul(modulus, 1, c, mod);
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
    nmod_mpoly_clear(T, ctx);
    nmod_mpoly_clear(G, ctx);
    nmod_mpoly_clear(Abar, ctx);
    nmod_mpoly_clear(Bbar, ctx);
    nmod_mpolyn_clear(Tn, ctx);
    nmod_mpolyn_clear(Gn, ctx);
    nmod_mpolyn_clear(Abarn, ctx);
    nmod_mpolyn_clear(Bbarn, ctx);
    n_polyun_clear(Aeh_cur);
    n_polyun_clear(Aeh_inc);
    n_polyun_clear(Beh_cur);
    n_polyun_clear(Beh_inc);
    n_poly_clear(gammaeh_cur);
    n_poly_clear(gammaeh_inc);
    n_poly_clear(modulus);
    n_poly_stack_clear(St->poly_stack);
    n_bpoly_stack_clear(St->bpoly_stack);

    flint_free(betas);
    flint_free(alphas);
    flint_randclear(state);

    for (i = 0; i < nvars; i++)
    {
        nmod_mpoly_clear(Aevals + i, ctx);
        nmod_mpoly_clear(Bevals + i, ctx);
        nmod_mpoly_clear(gammaevals + i, ctx);
    }
    flint_free(Aevals);
    flint_free(Bevals);
    flint_free(gammaevals);
/*
flint_printf("nmod_mpolyl_gcd_zippel_smprime returing %d\n", success);
*/
    return success;

gcd_is_trivial:

    nmod_mpoly_one(rG, ctx);
    nmod_mpoly_set(rAbar, A, ctx);
    nmod_mpoly_set(rBbar, B, ctx);

    success = 1;
    
    goto cleanup;
}


int nmod_mpolyuu_gcd_zippel_smprime(
    nmod_mpolyu_t rG, const slong * rGdegs, /* guess at rG degrees, could be NULL */
    nmod_mpolyu_t rAbar,
    nmod_mpolyu_t rBbar,
    const nmod_mpolyu_t A, const slong * Adegs,
    const nmod_mpolyu_t B, const slong * Bdegs,
    const nmod_mpoly_t gamma, const slong * gammadegs,
    const nmod_mpoly_ctx_t ctx)
{
    int success, use, alpha_tries_left;
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
    slong lastdeg;
    slong cur_zip_image, req_zip_images, abar_images, bbar_images, g_images;
    n_polyun_t Aeh, Beh;
    n_poly_t gammaeh;
    n_poly_t alphapow;
    nmod_mpolyu_struct * Aevals, * Bevals;
    nmod_mpoly_struct * gammaevals;
    nmod_poly_t modulus;
    n_poly_bpoly_stack_t St;
    mp_limb_t c, start_alpha;
    ulong GdegboundXY, newdegXY;
    slong degxAB, degyAB;
/*
flint_printf("nmod_mpolyuu_gcd_zippel_smprime called\n");
*/
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

        nmod_mpolyu_degrees_si(tmp_degs, A, ctx);
        for (j = 0; j < nvars; j++)
            FLINT_ASSERT(tmp_degs[j] == Adegs[j]);

        nmod_mpolyu_degrees_si(tmp_degs, B, ctx);
        for (j = 0; j < nvars; j++)
            FLINT_ASSERT(tmp_degs[j] == Bdegs[j]);

        nmod_mpoly_degrees_si(tmp_degs, gamma, ctx);
        for (j = 0; j < nvars; j++)
            FLINT_ASSERT(tmp_degs[j] == gammadegs[j]);

        flint_free(tmp_degs);
    }
#endif

    if (ctx->ffinfo->mod.n < 7)
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

    alpha_tries_left = 20;

choose_alphas:

    if (--alpha_tries_left < 0)
    {
        success = 0;
        goto cleanup;
    }

    for (i = 0; i < nvars; i++)
        alphas[i] = n_urandint(state, ctx->ffinfo->mod.n - 2) + 1;

    for (i = nvars - 1; i >= 0; i--)
    {
        nmod_mpolyu_evaluate_one_ui(Aevals + i, Aevals + i + 1, i, alphas[i], ctx);
        nmod_mpolyu_evaluate_one_ui(Bevals + i, Bevals + i + 1, i, alphas[i], ctx);
        nmod_mpoly_evaluate_one_ui(gammaevals + i, gammaevals + i + 1, i, alphas[i], ctx);
        nmod_mpoly_repack_bits_inplace(gammaevals + i, bits, ctx);

        if (nmod_mpoly_is_zero(gammaevals + i, ctx))
            goto choose_alphas;
        if (Aevals[i].length < 1 || Aevals[i].exps[0] != A->exps[0])
            goto choose_alphas;
        if (Bevals[i].length < 1 || Bevals[i].exps[0] != B->exps[0])
            goto choose_alphas;
    }

    m = 0;

    if (rGdegs == NULL)
        use = USE_G | USE_ABAR | USE_BBAR;
    else
        use = mpoly_gcd_get_use_first(rGdegs[m], Adegs[m], Bdegs[m], gammadegs[m]);

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

        alphas[m] = (alphas[m] < 2) ? mod.n - 1 : alphas[m] - 1;
        if (alphas[m] == start_alpha)
            goto choose_alphas;

        FLINT_ASSERT(alphapow->alloc > 1);
        alphapow->coeffs[0] = 1;
        alphapow->coeffs[1] = alphas[m];
        alphapow->length = 2;

        nmod_mpolyu_evaluate_one_ui(Aevals + m, Aevals + m + 1, m, alphas[m], ctx);
        nmod_mpolyu_evaluate_one_ui(Bevals + m, Bevals + m + 1, m, alphas[m], ctx);
        nmod_mpoly_evaluate_one_ui(gammaevals + m, gammaevals + m + 1, m, alphas[m], ctx);

        if (nmod_mpoly_is_zero(gammaevals + m, ctx))
            goto choose_alpha_0;
        if (Aevals[m].length < 1 || Aevals[m].exps[0] != A->exps[0])
            goto choose_alpha_0;
        if (Bevals[m].length < 1 || Bevals[m].exps[0] != B->exps[0])
            goto choose_alpha_0;

        nmod_mpolyuu_get_n_bpoly(Aev, Aevals + m, ctx);
        nmod_mpolyuu_get_n_bpoly(Bev, Bevals + m, ctx);

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

        gammaev = nmod_mpoly_get_ui(gammaevals + m, ctx);
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
                if (nmod_mpolyuu_divides(rG, B, rBbar, 2, ctx) &&
                    nmod_mpolyuu_divides(rAbar, A, rG, 2, ctx))
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

        if (nmod_poly_degree(modulus) > gammadegs[m] + Adegs[m] &&
            nmod_poly_degree(modulus) > gammadegs[m] + Bdegs[m])
        {
            goto choose_alphas;
        }

        c = nmod_neg(alphas[m], mod);
        n_poly_mod_shift_left_scalar_addmul((n_poly_struct*)modulus, 1, c, mod);
    }

    for (m = 1; m < nvars; m++)
    {
        /* G, Abar, Bbar are in Fp[gen(0), ..., gen(m - 1)] */

        nmod_mpolyun_interp_lift_sm_mpolyu(Gn, G, ctx);
        nmod_mpolyun_interp_lift_sm_mpolyu(Abarn, Abar, ctx);
        nmod_mpolyun_interp_lift_sm_mpolyu(Bbarn, Bbar, ctx);

        if (rGdegs == NULL)
        {
            use = USE_G | USE_ABAR | USE_BBAR;
        }
        else
        {
            slong numABgamma = gammaevals[m + 1].length +
                               nmod_mpolyu_total_length(Aevals + m + 1) +
                               nmod_mpolyu_total_length(Bevals + m + 1);

            use = nmod_mpoly_gcd_get_use(rGdegs[m], Adegs[m], Bdegs[m],
                      gammadegs[m], degxAB, degyAB, numABgamma, G, Abar, Bbar);
        }

        n_poly_one((n_poly_struct*)modulus);
        c = nmod_neg(alphas[m], mod);
        n_poly_mod_shift_left_scalar_addmul((n_poly_struct*)modulus, 1, c, mod);

        start_alpha = alphas[m];

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

            alphas[m] = (alphas[m] < 2) ? mod.n - 1 : alphas[m] - 1;
            if (alphas[m] == start_alpha)
                goto choose_alphas;

            nmod_mpolyu_evaluate_one_ui(Aevals + m, Aevals + m + 1, m, alphas[m], ctx);
            nmod_mpolyu_evaluate_one_ui(Bevals + m, Bevals + m + 1, m, alphas[m], ctx);
            nmod_mpoly_evaluate_one_ui(gammaevals + m, gammaevals + m + 1, m, alphas[m], ctx);
            if (nmod_mpoly_is_zero(gammaevals + m, ctx))
                goto choose_alpha_m;
            if (Aevals[m].length < 1 || Aevals[m].exps[0] != A->exps[0])
                goto choose_alpha_m;
            if (Bevals[m].length < 1 || Bevals[m].exps[0] != B->exps[0])
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
                if (Aev->length < 1 || n_bpoly_bidegree(Aev) != A->exps[0])
                    goto choose_betas;
                if (Bev->length < 1 || n_bpoly_bidegree(Bev) != B->exps[0])
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
/*
flint_printf("nmod_mpolyuu_gcd_zippel_smprime returning %d\n", success);
*/
FLINT_ASSERT(!success);

    return success;

gcd_is_trivial:

    nmod_mpolyu_one(rG, ctx);
    nmod_mpolyu_set(rAbar, A, ctx);
    nmod_mpolyu_set(rBbar, B, ctx);

    success = 1;
    
    goto cleanup;
}


slong fq_nmod_mpolyu_set_zip_form(
    n_fq_polyun_t H, /* monomial evals */
    n_fq_polyun_t M, /* master polys */
    const fq_nmod_mpolyu_t A,
    const fq_nmod_struct * betas,
    slong mvars,
    const fq_nmod_mpoly_ctx_t ctx);


void _n_fq_set_n_poly(
    mp_limb_t * a,
    const mp_limb_t * bcoeffs, slong blen,
    const fq_nmod_ctx_t ctx)
{
    slong d = fq_nmod_ctx_degree(ctx);

    if (blen > d)
    {
        _nmod_poly_rem(a, bcoeffs, blen, ctx->modulus->coeffs, d + 1, ctx->mod);
    }
    else
    {
        slong i;
        for (i = 0; i < blen; i++)
            a[i] = bcoeffs[i];
        for (; i < d; i++)
            a[i] = 0;
    }
}

void _n_poly_mul_n_fq(
    n_poly_t a,
    const n_poly_t b,
    const mp_limb_t * c,
    const fq_nmod_ctx_t ctx)
{
    slong d = fq_nmod_ctx_degree(ctx);
    n_poly_t C;
    C->coeffs = (mp_limb_t *) c;
    C->length = d;
    C->alloc = d;
    _n_poly_normalise(C);
    n_poly_mod_mul(a, b, C, ctx->mod);
}


static void nmod_mpolyuun_interp_lift_lg_bpoly(
    slong * lastdeg_,
    nmod_mpolyun_t A,
    const nmod_mpoly_ctx_t smctx,
    n_fq_bpoly_t B,
    const fq_nmod_mpoly_ctx_t lgctx)
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

            nmod_mpolyun_fit_length(A, Ai + 1, smctx);
            A->exps[Ai] = pack_exp2(i, j);
            nmod_mpolyn_fit_length(A->coeffs + Ai, 1, smctx);
            A->coeffs[Ai].length = 1;
            mpoly_monomial_zero(A->coeffs[Ai].exps + N*0, N);
            n_fq_get_fq_nmod(A->coeffs[Ai].coeffs + 0, Bi->coeffs + lgd*j, lgctx->fqctx);
            lastdeg = FLINT_MAX(lastdeg, nmod_poly_degree(A->coeffs[Ai].coeffs + 0));
            Ai++;
        }
    }
    A->length = Ai;

    *lastdeg_ = lastdeg;
}





void fq_nmod_mpolyuu_get_n_fq_bpoly(
    n_fq_bpoly_t A,
    fq_nmod_mpolyu_t B,
    const fq_nmod_mpoly_ctx_t ctx);

void n_fq_bpoly_scalar_mul_n_fq(
    n_fq_bpoly_t A,
    const mp_limb_t * c,
    const fq_nmod_ctx_t ctx);

void fq_nmod_mpolyuun_interp_lift_sm_bpoly(
    fq_nmod_mpolyun_t C,
    n_fq_bpoly_t B,
    const fq_nmod_mpoly_ctx_t ctx);

int nmod_mpolyuun_interp_crt_lg_bpoly(
    slong * lastdeg,
    nmod_mpolyun_t F,
    nmod_mpolyun_t T,
    n_fq_poly_t modulus,
    const nmod_mpoly_ctx_t smctx,
    n_fq_bpoly_t A,
    const fq_nmod_mpoly_ctx_t lgctx)
{
    slong lgd = fq_nmod_ctx_degree(lgctx->fqctx);
    int changed = 0;
    slong N = mpoly_words_per_exp(T->bits, smctx->minfo);
    n_fq_poly_struct * Acoeffs = A->coeffs;
    slong Fi, Ti, Ai, ai;
    slong Flen = F->length;
    ulong * Fexps = F->exps;
    nmod_mpolyn_struct * Fcoeffs = F->coeffs;
    ulong * Texps = T->exps;
    nmod_mpolyn_struct * Tcoeffs = T->coeffs;
    mp_limb_t * u = FLINT_ARRAY_ALLOC(3*lgd, mp_limb_t);
    mp_limb_t * v = u + lgd;
    mp_limb_t * inv_m_eval = v + lgd;

    _n_fq_set_n_poly(u, modulus->coeffs, modulus->length, lgctx->fqctx);

    n_fq_inv(inv_m_eval, u, lgctx->fqctx);

    FLINT_ASSERT(nmod_mpolyun_is_canonical(F, smctx));
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
            nmod_mpolyun_fit_length(T, Ti + extra + 1, smctx);
            Tcoeffs = T->coeffs;
            Texps = T->exps;
        }

        nmod_mpolyn_fit_length(Tcoeffs + Ti, 1, smctx);
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

            _n_fq_set_n_poly(u, getcfn(Fcoeffs + Fi)->coeffs, getcfn(Fcoeffs + Fi)->length,lgctx->fqctx);
            n_fq_sub(v, Acoeffs[Ai].coeffs + lgd*ai, u, lgctx->fqctx);
            if (!_n_fq_is_zero(v, lgd))
            {
                changed = 1;
                n_fq_mul(u, v, inv_m_eval, lgctx->fqctx);
                _n_poly_mul_n_fq(getcfn(Tcoeffs + Ti), modulus, u, lgctx->fqctx);
                n_poly_mod_add(getcfn(Tcoeffs + Ti), getcfn(Tcoeffs + Ti), getcfn(Fcoeffs + Fi), smctx->ffinfo->mod);
            }
            else
            {
                n_poly_set(getcfn(Tcoeffs + Ti), getcfn(Fcoeffs + Fi));
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

            _n_fq_set_n_poly(v, getcfn(Fcoeffs + Fi)->coeffs, getcfn(Fcoeffs + Fi)->length,lgctx->fqctx);
            if (!_n_fq_is_zero(v, lgd))
            {
                changed = 1;
                n_fq_mul(u, v, inv_m_eval, lgctx->fqctx);
                _n_poly_mul_n_fq(getcfn(Tcoeffs + Ti), modulus, u, lgctx->fqctx);
                n_poly_mod_sub(getcfn(Tcoeffs + Ti), getcfn(Fcoeffs + Fi), getcfn(Tcoeffs + Ti), smctx->ffinfo->mod);
            }
            else
            {
                n_poly_set(getcfn(Tcoeffs + Ti), getcfn(Fcoeffs + Fi));
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
            _n_poly_mul_n_fq(getcfn(Tcoeffs + Ti), modulus, u, lgctx->fqctx);

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

        FLINT_ASSERT(!n_poly_is_zero(getcfn(Tcoeffs + Ti)));
        *lastdeg = FLINT_MAX(*lastdeg, n_poly_degree(getcfn(Tcoeffs + Ti)));
        Ti++;
    }

    T->length = Ti;

    if (changed)
        nmod_mpolyun_swap(T, F);

    FLINT_ASSERT(nmod_mpolyun_is_canonical(F, smctx));

    flint_free(u);

    return changed;
}



int nmod_mpolyn_interp_crt_lg_bpoly(
    slong * lastdeg,
    nmod_mpolyn_t F,
    nmod_mpolyn_t T,
    n_fq_poly_t modulus,
    const nmod_mpoly_ctx_t smctx,
    n_fq_bpoly_t A,
    const fq_nmod_mpoly_ctx_t lgctx)
{
    slong lgd = fq_nmod_ctx_degree(lgctx->fqctx);
    int changed = 0;
    slong N = mpoly_words_per_exp(T->bits, smctx->minfo);
    slong off0, shift0, off1, shift1;
    n_fq_poly_struct * Acoeffs = A->coeffs;
    slong Fi, Ti, Ai, ai;
    slong Flen = F->length;
    ulong * Fexps = F->exps;
    nmod_poly_struct * Fcoeffs = F->coeffs;
    ulong * Texps = T->exps;
    nmod_poly_struct * Tcoeffs = T->coeffs;
    mp_limb_t * u = FLINT_ARRAY_ALLOC(3*lgd, mp_limb_t);
    mp_limb_t * v = u + lgd;
    mp_limb_t * inv_m_eval = v + lgd;
    ulong Fexpi, mask;

    mask = (-UWORD(1)) >> (FLINT_BITS - F->bits);
    mpoly_gen_offset_shift_sp(&off0, &shift0, 0, F->bits, smctx->minfo);
    mpoly_gen_offset_shift_sp(&off1, &shift1, 1, F->bits, smctx->minfo);

    _n_fq_set_n_poly(u, modulus->coeffs, modulus->length, lgctx->fqctx);

    n_fq_inv(inv_m_eval, u, lgctx->fqctx);

    FLINT_ASSERT(nmod_mpolyn_is_canonical(F, smctx));
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
            nmod_mpolyn_fit_length(T, Ti + extra + 1, smctx);
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

            _n_fq_set_n_poly(u, Fcoeffs[Fi].coeffs, Fcoeffs[Fi].length, lgctx->fqctx);
            n_fq_sub(v, Acoeffs[Ai].coeffs + lgd*ai, u, lgctx->fqctx);
            if (!_n_fq_is_zero(v, lgd))
            {
                changed = 1;
                n_fq_mul(u, v, inv_m_eval, lgctx->fqctx);
                _n_poly_mul_n_fq((n_poly_struct*)(Tcoeffs + Ti), modulus, u, lgctx->fqctx);
                n_poly_mod_add((n_poly_struct*)(Tcoeffs + Ti), (n_poly_struct*)(Tcoeffs + Ti), (n_poly_struct*)(Fcoeffs + Fi), smctx->ffinfo->mod);
            }
            else
            {
                n_poly_set((n_poly_struct*)(Tcoeffs + Ti), (n_poly_struct*)(Fcoeffs + Fi));
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

            _n_fq_set_n_poly(v, Fcoeffs[Fi].coeffs, Fcoeffs[Fi].length, lgctx->fqctx);
            if (!_n_fq_is_zero(v, lgd))
            {
                changed = 1;
                n_fq_mul(u, v, inv_m_eval, lgctx->fqctx);
                _n_poly_mul_n_fq((n_poly_struct*)(Tcoeffs + Ti), modulus, u, lgctx->fqctx);
                n_poly_mod_sub((n_poly_struct*)(Tcoeffs + Ti), (n_poly_struct*)(Fcoeffs + Fi), (n_poly_struct*)(Tcoeffs + Ti), smctx->ffinfo->mod);
            }
            else
            {
                n_poly_set((n_poly_struct*)(Tcoeffs + Ti), (n_poly_struct*)(Fcoeffs + Fi));
            }

            Fi++;
        }
        else
        {
            FLINT_ASSERT(Ai >= 0 && (Fi >= Flen || Fexps[Fi] < pack_exp2(Ai, ai)));

            /* F term missing, A term ok */
            mpoly_monomial_zero(Texps + N*Ti, N);
            (Texps + N*Ti)[off0] += (Ai << shift0);
            (Texps + N*Ti)[off1] += (ai << shift1);

            changed = 1;
            n_fq_mul(u, Acoeffs[Ai].coeffs + lgd*ai, inv_m_eval, lgctx->fqctx);
            _n_poly_mul_n_fq((n_poly_struct*)(Tcoeffs + Ti), modulus, u, lgctx->fqctx);

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

        FLINT_ASSERT(!n_poly_is_zero((n_poly_struct*)(Tcoeffs + Ti)));
        *lastdeg = FLINT_MAX(*lastdeg, n_poly_degree((n_poly_struct*)(Tcoeffs + Ti)));
        Ti++;
    }

    T->length = Ti;

    if (changed)
        nmod_mpolyn_swap(T, F);

    FLINT_ASSERT(nmod_mpolyn_is_canonical(F, smctx));

    flint_free(u);

    return changed;
}



int newnmod_mpolyn_interp_mcrt_lg_mpoly(
    slong * lastdeg_,
    nmod_mpolyn_t H,
    const nmod_mpoly_ctx_t smctx,
    const n_poly_t m,
    const mp_limb_t * inv_m_eval,
    fq_nmod_mpoly_t A,
    const fq_nmod_mpoly_ctx_t lgctx)
{
    slong lastdeg = *lastdeg_;
    slong i, lgd = fq_nmod_ctx_degree(lgctx->fqctx);
#if FLINT_WANT_ASSERT
    slong N = mpoly_words_per_exp(A->bits, smctx->minfo);
#endif
    int changed = 0;
    mp_limb_t * u = FLINT_ARRAY_ALLOC(lgd, mp_limb_t);
    n_poly_t w;

    n_poly_init(w);

    FLINT_ASSERT(H->length == A->length);
    FLINT_ASSERT(H->bits == A->bits);

    for (i = 0; i < A->length; i++)
    {
        FLINT_ASSERT(mpoly_monomial_equal(H->exps + N*i, A->exps + N*i, N));
        _n_fq_set_n_poly(u, H->coeffs[i].coeffs, H->coeffs[i].length, lgctx->fqctx);
        n_fq_sub(u, A->coeffs + lgd*i, u, lgctx->fqctx);
        if (!_n_fq_is_zero(u, lgd))
        {
            changed = 1;
            n_fq_mul(u, u, inv_m_eval, lgctx->fqctx);
            _n_poly_mul_n_fq(w, m, u, lgctx->fqctx);
            {
                nmod_poly_t wmock;
                nmod_poly_mock(wmock, w, smctx->ffinfo->mod);
                nmod_poly_add(H->coeffs + i, H->coeffs + i, wmock);
            }
        }

        lastdeg = FLINT_MAX(lastdeg, nmod_poly_degree(H->coeffs + i));

        FLINT_ASSERT(nmod_poly_degree(H->coeffs + i) < n_poly_degree(m)
                                     + nmod_poly_degree(lgctx->fqctx->modulus));
    }

    *lastdeg_ = lastdeg;

    flint_free(u);
    n_poly_clear(w);
    
    return changed;
}

int newnmod_mpolyun_interp_mcrt_lg_mpolyu(
    slong * lastdeg,
    nmod_mpolyun_t H,
    const nmod_mpoly_ctx_t smctx,
    const n_poly_t m,
    fq_nmod_mpolyu_t A,
    const fq_nmod_mpoly_ctx_t lgctx)
{
    slong lgd = fq_nmod_ctx_degree(lgctx->fqctx);
    slong i;
    int changed = 0;
    mp_limb_t * inv_m_eval = FLINT_ARRAY_ALLOC(lgd, mp_limb_t);

    *lastdeg = -WORD(1);

    _n_fq_set_n_poly(inv_m_eval, m->coeffs, m->length, lgctx->fqctx);
    n_fq_inv(inv_m_eval, inv_m_eval, lgctx->fqctx);

    FLINT_ASSERT(H->bits == A->bits);
    FLINT_ASSERT(H->length == A->length);
    for (i = 0; i < A->length; i++)
    {
        FLINT_ASSERT(H->exps[i] == A->exps[i]);
        changed |= newnmod_mpolyn_interp_mcrt_lg_mpoly(lastdeg, H->coeffs + i,
                                   smctx, m, inv_m_eval, A->coeffs + i, lgctx);
    }
    H->length = A->length;

    flint_free(inv_m_eval);

    return changed;
}

void fq_nmod_mpoly_set_eval_helper(
    n_fq_poly_t EH,
    const fq_nmod_mpoly_t A,
    const fq_nmod_struct * betas,
    slong mvars,
    const fq_nmod_mpoly_ctx_t ctx);


void fq_nmod_mpolyu_set_eval_helper(
    n_fq_polyun_t EH,
    const fq_nmod_mpolyu_t A,
    const fq_nmod_struct * betas,
    slong mvars,
    const fq_nmod_mpoly_ctx_t ctx);

void n_fq_poly_eval_step(
    mp_limb_t * res,
    n_fq_poly_t A,
    const fq_nmod_ctx_t ctx);

void n_fq_poly_evalp_step(
    mp_limb_t * res,
    n_fq_poly_t A,
    const fq_nmod_ctx_t ctx);

void n_fq_bpoly_eval_step(
    n_fq_bpoly_t E,
    n_fq_polyun_t A,
    const fq_nmod_ctx_t ctx);

int n_fq_polyu2_add_zip_must_match(
    n_fq_polyun_t Z,
    const n_fq_bpoly_t A,
    slong cur_length,
    const fq_nmod_ctx_t ctx);

int qzip_solve(
    fq_nmod_mpolyu_t A,
    n_fq_polyun_t Z,
    n_fq_polyun_t H,
    n_fq_polyun_t M,
    const fq_nmod_mpoly_ctx_t ctx);

int fq_nmod_mpolyun_interp_mcrt_sm_mpolyu(
    slong * lastdeg,
    fq_nmod_mpolyun_t F,
    fq_nmod_mpolyu_t A, /* could have zero coeffs */
    const n_fq_poly_t modulus,
    n_fq_poly_t alphapow,
    const fq_nmod_mpoly_ctx_t ctx);


void fq_nmod_mpolyuu_print_pretty(
    const fq_nmod_mpolyu_t poly,
    const char ** x,
    slong nmainvars,
    const fq_nmod_mpoly_ctx_t ctx);

int nmod_mpolyuu_gcd_zippel_lgprime(
    nmod_mpolyu_t rG, const slong * rGdegs,
    nmod_mpolyu_t rAbar,
    nmod_mpolyu_t rBbar,
    const nmod_mpolyu_t A, const slong * Adegs,
    const nmod_mpolyu_t B, const slong * Bdegs,
    const nmod_mpoly_t gamma, const slong * gammadegs,
    const nmod_mpoly_ctx_t smctx)
{
    slong lgd;
    int success, use;
    slong i, j, m;
    slong nvars = smctx->minfo->nvars;
    flint_bitcnt_t bits = A->bits;
    fq_nmod_struct * alphas, * betas;
    flint_rand_t state;
    nmod_mpoly_t cont;
    nmod_mpolyu_t T, G, Abar, Bbar;
    nmod_mpolyun_t Tn, Gn, Abarn, Bbarn;
    fq_nmod_mpolyu_t qT, qG, qAbar, qBbar;
    fq_nmod_mpolyu_t qrG, qrAbar, qrBbar;
    fq_nmod_mpolyun_t qTn, qGn, qAbarn, qBbarn;
    n_fq_polyun_t HG, HAbar, HBbar, MG, MAbar, MBbar, ZG, ZAbar, ZBbar;
    n_fq_bpoly_t Aev, Bev, Gev, Abarev, Bbarev;
    const mp_limb_t * gammaev;
    slong lastdeg;
    slong cur_zip_image, req_zip_images, this_images;
    n_fq_polyun_t Aeh, Beh;
    n_fq_poly_t gammaeh;
    fq_nmod_mpolyu_struct * Aevals, * Bevals;
    fq_nmod_mpoly_struct * gammaevals;
    n_fq_poly_t alphapow;
    n_poly_t modulus;
    n_poly_bpoly_stack_t St;
    n_poly_t tmp;  /* tmp arithmetic space */
    fq_nmod_t start_alpha;
    ulong GdegboundXY, newdegXY;
    fq_nmod_mpoly_ctx_t lgctx;
    nmod_mpolyn_t gamman;
    nmod_mpolyun_t An, Bn;
    slong degxAB, degyAB;
/*
flint_printf("nmod_mpolyuu_gcd_zippel_lgprime called\n");
*/

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

        nmod_mpolyu_degrees_si(tmp_degs, A, smctx);
        for (j = 0; j < nvars; j++)
            FLINT_ASSERT(tmp_degs[j] == Adegs[j]);

        nmod_mpolyu_degrees_si(tmp_degs, B, smctx);
        for (j = 0; j < nvars; j++)
            FLINT_ASSERT(tmp_degs[j] == Bdegs[j]);

        nmod_mpoly_degrees_si(tmp_degs, gamma, smctx);
        for (j = 0; j < nvars; j++)
            FLINT_ASSERT(tmp_degs[j] == gammadegs[j]);

        flint_free(tmp_degs);
    }
#endif

    flint_randinit(state);

    lgd = WORD(20)/(FLINT_BIT_COUNT(smctx->ffinfo->mod.n));
    lgd = FLINT_MAX(WORD(2), lgd);
    fq_nmod_mpoly_ctx_init_deg(lgctx, nvars, ORD_LEX, smctx->ffinfo->mod.n, lgd);
    n_poly_init2(tmp, lgd);
    n_poly_init2(alphapow, 2*lgd);

    fq_nmod_init(start_alpha, lgctx->fqctx);
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
    nmod_mpoly_init3(cont, 1, bits, smctx);

    nmod_mpolyu_init(T, bits, smctx);
    nmod_mpolyu_init(G, bits, smctx);
    nmod_mpolyu_init(Abar, bits, smctx);
    nmod_mpolyu_init(Bbar, bits, smctx);
    nmod_mpolyun_init(Tn, bits, smctx);
    nmod_mpolyun_init(Gn, bits, smctx);
    nmod_mpolyun_init(Abarn, bits, smctx);
    nmod_mpolyun_init(Bbarn, bits, smctx);

    fq_nmod_mpolyu_init(qT, bits, lgctx);
    fq_nmod_mpolyu_init(qG, bits, lgctx);
    fq_nmod_mpolyu_init(qAbar, bits, lgctx);
    fq_nmod_mpolyu_init(qBbar, bits, lgctx);
    fq_nmod_mpolyu_init(qrG, bits, lgctx);
    fq_nmod_mpolyu_init(qrAbar, bits, lgctx);
    fq_nmod_mpolyu_init(qrBbar, bits, lgctx);
    fq_nmod_mpolyun_init(qTn, bits, lgctx);
    fq_nmod_mpolyun_init(qGn, bits, lgctx);
    fq_nmod_mpolyun_init(qAbarn, bits, lgctx);
    fq_nmod_mpolyun_init(qBbarn, bits, lgctx);

    n_fq_polyun_init(Aeh);
    n_fq_polyun_init(Beh);
    n_fq_poly_init(gammaeh);

    n_poly_init(modulus);

    n_poly_stack_init(St->poly_stack);
    n_bpoly_stack_init(St->bpoly_stack);

    nmod_mpolyun_init(An, bits, smctx);
    nmod_mpolyun_init(Bn, bits, smctx);
    nmod_mpolyn_init(gamman, bits, smctx);

    /* alphas[nvars - 1] not used - it is replaced lgctx->fqctx->modulus */
    betas = FLINT_ARRAY_ALLOC(nvars, fq_nmod_struct);
    alphas = FLINT_ARRAY_ALLOC(nvars, fq_nmod_struct);
    for (i = 0; i < nvars; i++)
    {
        fq_nmod_init(betas + i, lgctx->fqctx);
        fq_nmod_init(alphas + i, lgctx->fqctx);
    }

    /* Aevals[nvars] does not exist - it is replaced by An */
    Aevals = FLINT_ARRAY_ALLOC(nvars, fq_nmod_mpolyu_struct);
    Bevals = FLINT_ARRAY_ALLOC(nvars, fq_nmod_mpolyu_struct);
    gammaevals = FLINT_ARRAY_ALLOC(nvars, fq_nmod_mpoly_struct);
    for (i = 0; i < nvars; i++)
    {
        fq_nmod_mpolyu_init(Aevals + i, bits, lgctx);
        fq_nmod_mpolyu_init(Bevals + i, bits, lgctx);
        fq_nmod_mpoly_init3(gammaevals + i, 0, bits, lgctx);
    }
    nmod_mpolyu_cvtto_mpolyun(An, A, nvars - 1, smctx);
    nmod_mpolyu_cvtto_mpolyun(Bn, B, nvars - 1, smctx);
    nmod_mpoly_cvtto_mpolyn(gamman, gamma, nvars - 1, smctx);

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

    goto got_alpha_m;

increase_degree:

choose_alphas:

    /* TODO: don't necessarily increase degree here */
    lgd++;
    if (lgd > 10000)
    {
        /* ran out of primes */
        success = 0;
        goto cleanup;
    }

    fq_nmod_mpoly_ctx_change_modulus(lgctx, lgd);
    n_poly_fit_length(tmp, lgd);
    n_poly_fit_length(alphapow, 2*lgd);

got_alpha_m:

    for (i = 0; i < nvars - 1; i++)
    {
        fq_nmod_rand(alphas + i, state, lgctx->fqctx);
        if (fq_nmod_is_zero(alphas + i, lgctx->fqctx))
            fq_nmod_one(alphas + i, lgctx->fqctx);
    }

    i = nvars - 1;
    nmod_mpolyun_interp_reduce_lg_mpolyu(Aevals + i, An, lgctx, smctx);
    nmod_mpolyun_interp_reduce_lg_mpolyu(Bevals + i, Bn, lgctx, smctx);
    nmod_mpolyn_interp_reduce_lg_mpoly(gammaevals + i, gamman, lgctx, smctx);
    if (fq_nmod_mpoly_is_zero(gammaevals + i, lgctx))
        goto choose_alphas;
    if (Aevals[i].length < 1 || Aevals[i].exps[0] != A->exps[0])
        goto choose_alphas;
    if (Bevals[i].length < 1 || Bevals[i].exps[0] != B->exps[0])
        goto choose_alphas;
    for (i--; i >= 0; i--)
    {
        fq_nmod_mpolyu_evaluate_one_fq_nmod(Aevals + i, Aevals + i + 1, i, alphas + i, lgctx);
        fq_nmod_mpolyu_evaluate_one_fq_nmod(Bevals + i, Bevals + i + 1, i, alphas + i, lgctx);
        fq_nmod_mpoly_evaluate_one_fq_nmod(gammaevals + i, gammaevals + i + 1, i, alphas + i, lgctx);
        fq_nmod_mpoly_repack_bits_inplace(gammaevals + i, bits, lgctx);
        if (fq_nmod_mpoly_is_zero(gammaevals + i, lgctx))
            goto choose_alphas;
        if (Aevals[i].length < 1 || Aevals[i].exps[0] != A->exps[0])
            goto choose_alphas;
        if (Bevals[i].length < 1 || Bevals[i].exps[0] != B->exps[0])
            goto choose_alphas;
    }

    m = 0;

    if (rGdegs == NULL)
        use = USE_G | USE_ABAR | USE_BBAR;
    else
        use = mpoly_gcd_get_use_first(rGdegs[m], Adegs[m], Bdegs[m], gammadegs[m]);

    fq_nmod_mpolyuu_get_n_fq_bpoly(Aev, Aevals + m, lgctx);
    fq_nmod_mpolyuu_get_n_fq_bpoly(Bev, Bevals + m, lgctx);

    success = n_fq_bpoly_gcd_brown_smprime(Gev, Abarev, Bbarev, Aev, Bev, lgctx->fqctx, St);
    if (!success)
        goto increase_degree;

    newdegXY = n_bpoly_bidegree(Gev);
    if (newdegXY > GdegboundXY)
        goto choose_alphas;
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
        nmod_mpolyuun_interp_lift_lg_bpoly(&lastdeg, Gn, smctx, Gev, lgctx);
        nmod_mpolyuun_interp_lift_lg_bpoly(&lastdeg, Abarn, smctx, Abarev, lgctx);
        nmod_mpolyuun_interp_lift_lg_bpoly(&lastdeg, Bbarn, smctx, Bbarev, lgctx);

        n_poly_set_nmod_poly(modulus, lgctx->fqctx->modulus);

        while (1)
        {
        choose_alpha_0_last:
            lgd++;
            if (lgd > 10000)
            {
                /* ran out of primes */
                success = 0;
                goto cleanup;
            }

            fq_nmod_mpoly_ctx_change_modulus(lgctx, lgd);
            n_poly_fit_length(tmp, lgd);
            n_poly_fit_length(alphapow, 2*lgd);

            nmod_mpolyun_interp_reduce_lg_mpolyu(Aevals + m, An, lgctx, smctx);
            nmod_mpolyun_interp_reduce_lg_mpolyu(Bevals + m, Bn, lgctx, smctx);
            nmod_mpolyn_interp_reduce_lg_mpoly(gammaevals + m, gamman, lgctx, smctx);
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
                goto choose_alphas;
            }

            gammaev = fq_nmod_mpoly_get_nonzero_n_fq(gammaevals + m, lgctx);
            n_fq_bpoly_scalar_mul_n_fq(Gev, gammaev, lgctx->fqctx);

            if ((use & USE_G) && !nmod_mpolyuun_interp_crt_lg_bpoly(
                                 &lastdeg, Gn, Tn, modulus, smctx, Gev, lgctx))
            {
                nmod_mpolyu_cvtfrom_mpolyun(rG, Gn, m, smctx);
                success = nmod_mpolyu_content_mpoly_threaded_pool(cont, rG, smctx, NULL, 0);
                if (!success)
                    goto cleanup;
                nmod_mpolyu_divexact_mpoly_inplace(rG, cont, smctx);
                if (nmod_mpolyuu_divides(rAbar, A, rG, 2, smctx) &&
                    nmod_mpolyuu_divides(rBbar, B, rG, 2, smctx))
                {
                    break;
                }
            }

            if ((use & USE_ABAR) && !nmod_mpolyuun_interp_crt_lg_bpoly(
                           &lastdeg, Abarn, Tn, modulus, smctx, Abarev, lgctx))
            {
                nmod_mpolyu_cvtfrom_mpolyun(rAbar, Abarn, m, smctx);
                success = nmod_mpolyu_content_mpoly_threaded_pool(cont, rAbar, smctx, NULL, 0);
                if (!success)
                    goto cleanup;
                nmod_mpolyu_divexact_mpoly_inplace(rAbar, cont, smctx);
                if (nmod_mpolyuu_divides(rG, A, rAbar, 2, smctx) &&
                    nmod_mpolyuu_divides(rBbar, B, rG, 2, smctx))
                {
                    break;
                }
            }

            if ((use & USE_BBAR) && !nmod_mpolyuun_interp_crt_lg_bpoly(
                          &lastdeg, Bbarn, Tn, modulus, smctx, Bbarev, lgctx))
            {
                nmod_mpolyu_cvtfrom_mpolyun(rBbar, Bbarn, m, smctx);
                success = nmod_mpolyu_content_mpoly_threaded_pool(cont, rBbar, smctx, NULL, 0);
                if (!success)
                    goto cleanup;
                nmod_mpolyu_divexact_mpoly_inplace(rBbar, cont, smctx);
                if (nmod_mpolyuu_divides(rG, B, rBbar, 2, smctx) &&
                    nmod_mpolyuu_divides(rAbar, A, rG, 2, smctx))
                {
                    break;
                }
            }

            if (n_fq_poly_degree(modulus) > gammadegs[m] + Adegs[m] &&
                n_fq_poly_degree(modulus) > gammadegs[m] + Bdegs[m])
            {
                goto choose_alphas;
            }

            {
                n_poly_t h_mock;
                n_poly_mock(h_mock, lgctx->fqctx->modulus);
                n_poly_mod_mul(modulus, modulus, h_mock, smctx->ffinfo->mod);
            }
        }

        success = 1;
        goto cleanup;
    }

    fq_nmod_mpolyuun_interp_lift_sm_bpoly(qGn, Gev, lgctx);
    fq_nmod_mpolyuun_interp_lift_sm_bpoly(qAbarn, Abarev, lgctx);
    fq_nmod_mpolyuun_interp_lift_sm_bpoly(qBbarn, Bbarev, lgctx);

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
        alphapow->length = 2;
        _n_fq_one(alphapow->coeffs + lgd*0, lgd);
        n_fq_set_fq_nmod(alphapow->coeffs + lgd*1, alphas + m, lgctx->fqctx);

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
            goto choose_alphas;
        }

        gammaev = fq_nmod_mpoly_get_nonzero_n_fq(gammaevals + m, lgctx);
        n_fq_bpoly_scalar_mul_n_fq(Gev, gammaev, lgctx->fqctx);

        FLINT_ASSERT(tmp->alloc >= lgd);
        n_fq_poly_eval_pow(tmp->coeffs, modulus, alphapow, lgctx->fqctx);
        n_fq_inv(tmp->coeffs, tmp->coeffs, lgctx->fqctx);
        n_fq_poly_scalar_mul_n_fq(modulus, modulus, tmp->coeffs, lgctx->fqctx);

        if ((use & USE_G) && !fq_nmod_mpolyuun_interp_crt_sm_bpoly(
                              &lastdeg, qGn, qTn, Gev, modulus, alphapow, lgctx))
        {
            fq_nmod_mpolyu_cvtfrom_mpolyun(qG, qGn, m, lgctx);
            fq_nmod_mpolyu_mul_mpoly(qT, Aevals + m + 1, gammaevals + m + 1, lgctx);
            if (fq_nmod_mpolyuu_divides(qAbar, qT, qG, 2, lgctx))
            {
                fq_nmod_mpolyu_mul_mpoly(qT, Bevals + m + 1, gammaevals + m + 1, lgctx);
                if (fq_nmod_mpolyuu_divides(qBbar, qT, qG, 2, lgctx))
                {
                    break;
                }
            }
        }

        if ((use & USE_ABAR) && !fq_nmod_mpolyuun_interp_crt_sm_bpoly(
                      &lastdeg, qAbarn, qTn, Abarev, modulus, alphapow, lgctx))
        {
            fq_nmod_mpolyu_cvtfrom_mpolyun(qAbar, qAbarn, m, lgctx);
            fq_nmod_mpolyu_mul_mpoly(qT, Aevals + m + 1, gammaevals + m + 1, lgctx);
            if (fq_nmod_mpolyuu_divides(qG, qT, qAbar, 2, lgctx))
            {
                fq_nmod_mpolyu_mul_mpoly(qT, Bevals + m + 1, gammaevals + m + 1, lgctx);
                if (fq_nmod_mpolyuu_divides(qBbar, qT, qG, 2, lgctx))
                    break;
            }
        }

        if ((use & USE_BBAR) && !fq_nmod_mpolyuun_interp_crt_sm_bpoly(
                          &lastdeg, qBbarn, qTn, Bbarev, modulus, alphapow, lgctx))
        {
            fq_nmod_mpolyu_cvtfrom_mpolyun(qBbar, qBbarn, m, lgctx);
            fq_nmod_mpolyu_mul_mpoly(qT, Bevals + m + 1, gammaevals + m + 1, lgctx);
            if (fq_nmod_mpolyuu_divides(qG, qT, qBbar, 2, lgctx))
            {
                fq_nmod_mpolyu_mul_mpoly(qT, Aevals + m + 1, gammaevals + m + 1, lgctx);
                if (fq_nmod_mpolyuu_divides(qAbar, qT, qG, 2, lgctx))
                    break;
            }
        }

        if (n_fq_poly_degree(modulus) > gammadegs[m] + Adegs[m] &&
            n_fq_poly_degree(modulus) > gammadegs[m] + Bdegs[m])
        {
            goto choose_alphas;
        }

        FLINT_ASSERT(n_fq_equal_fq_nmod(alphapow->coeffs + lgd*1, alphas + m, lgctx->fqctx));
        n_fq_poly_shift_left_scalar_submul(modulus, 1, alphapow->coeffs + lgd*1, lgctx->fqctx);
    }

    for (m = 1; m < nvars - 1; m++)
    {
        /* qG, qAbar, qBbar are in Fq[gen(0), ..., gen(m - 1)] */
        fq_nmod_mpolyun_interp_lift_sm_mpolyu(qGn, qG, lgctx);
        fq_nmod_mpolyun_interp_lift_sm_mpolyu(qAbarn, qAbar, lgctx);
        fq_nmod_mpolyun_interp_lift_sm_mpolyu(qBbarn, qBbar, lgctx);

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
                   gammadegs[m], degxAB, degyAB, numABgamma, qG, qAbar, qBbar);
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
            this_images = fq_nmod_mpolyu_set_zip_form(HG, MG, qG, betas, m, lgctx);
            req_zip_images = FLINT_MAX(req_zip_images, this_images);
        }
        if (use & USE_ABAR)
        {
            this_images =  fq_nmod_mpolyu_set_zip_form(HAbar, MAbar, qAbar, betas, m, lgctx);
            req_zip_images = FLINT_MAX(req_zip_images, this_images);            
        }
        if (use & USE_BBAR)
        {
            this_images =  fq_nmod_mpolyu_set_zip_form(HBbar, MBbar, qBbar, betas, m, lgctx);
            req_zip_images = FLINT_MAX(req_zip_images, this_images);            
        }

        while (1)
        {
choose_alpha_m:

            fq_nmod_next_not_zero(alphas + m, lgctx->fqctx);
            if (fq_nmod_equal(alphas + m, start_alpha, lgctx->fqctx))
                goto increase_degree;

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
                n_poly_fit_length(tmp, lgd);
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
                    goto choose_alphas;
                }

                n_fq_bpoly_scalar_mul_n_fq(Gev, tmp->coeffs, lgctx->fqctx);
                if ((use & USE_G) && !n_fq_polyu2_add_zip_must_match(ZG, Gev, cur_zip_image, lgctx->fqctx))
                    goto choose_alphas;
                if ((use & USE_ABAR) && !n_fq_polyu2_add_zip_must_match(ZAbar, Abarev, cur_zip_image, lgctx->fqctx))
                    goto choose_alphas;
                if ((use & USE_BBAR) && !n_fq_polyu2_add_zip_must_match(ZBbar, Bbarev, cur_zip_image, lgctx->fqctx))
                    goto choose_alphas;
            }

            if ((use & USE_G) && qzip_solve(qG, ZG, HG, MG, lgctx) < 1)
                goto choose_alphas;
            if ((use & USE_ABAR) && qzip_solve(qAbar, ZAbar, HAbar, MAbar, lgctx) < 1)
                goto choose_alphas;
            if ((use & USE_BBAR) && qzip_solve(qBbar, ZBbar, HBbar, MBbar, lgctx) < 1)
                goto choose_alphas;

            FLINT_ASSERT(n_fq_poly_degree(modulus) > 0);

            FLINT_ASSERT(tmp->alloc >= lgd);
            n_fq_poly_eval_pow(tmp->coeffs, modulus, alphapow, lgctx->fqctx);
            n_fq_inv(tmp->coeffs, tmp->coeffs, lgctx->fqctx);
            n_fq_poly_scalar_mul_n_fq(modulus, modulus, tmp->coeffs, lgctx->fqctx);

            if ((use & USE_G) && !fq_nmod_mpolyun_interp_mcrt_sm_mpolyu(
                                      &lastdeg, qGn, qG, modulus, alphapow, lgctx))
            {
                fq_nmod_mpolyu_cvtfrom_mpolyun(qrG, qGn, m, lgctx);
                fq_nmod_mpolyu_mul_mpoly(qT, Aevals + m + 1, gammaevals + m + 1, lgctx);
                if (fq_nmod_mpolyuu_divides(qrAbar, qT, qrG, 2, lgctx))
                {
                    fq_nmod_mpolyu_mul_mpoly(qT, Bevals + m + 1, gammaevals + m + 1, lgctx);
                    if (fq_nmod_mpolyuu_divides(qrBbar, qT, qrG, 2, lgctx))
                    {
                        fq_nmod_mpolyu_swap(qG, qrG, lgctx);
                        fq_nmod_mpolyu_swap(qAbar, qrAbar, lgctx);
                        fq_nmod_mpolyu_swap(qBbar, qrBbar, lgctx);
                        break;
                    }
                }
            }

            if ((use & USE_ABAR) && !fq_nmod_mpolyun_interp_mcrt_sm_mpolyu(
                              &lastdeg, qAbarn, qAbar, modulus, alphapow, lgctx))
            {
                fq_nmod_mpolyu_cvtfrom_mpolyun(qrAbar, qAbarn, m, lgctx);
                fq_nmod_mpolyu_mul_mpoly(qT, Aevals + m + 1, gammaevals + m + 1, lgctx);
                if (fq_nmod_mpolyuu_divides(qrG, qT, qrAbar, 2, lgctx))
                {
                    fq_nmod_mpolyu_mul_mpoly(qT, Bevals + m + 1, gammaevals + m + 1, lgctx);
                    if (fq_nmod_mpolyuu_divides(qrBbar, qT, qrG, 2, lgctx))
                    {
                        fq_nmod_mpolyu_swap(qG, qrG, lgctx);
                        fq_nmod_mpolyu_swap(qAbar, qrAbar, lgctx);
                        fq_nmod_mpolyu_swap(qBbar, qrBbar, lgctx);
                        break;
                    }
                }
            }

            if ((use & USE_BBAR) && !fq_nmod_mpolyun_interp_mcrt_sm_mpolyu(
                              &lastdeg, qBbarn, qBbar, modulus, alphapow, lgctx))
            {
                fq_nmod_mpolyu_cvtfrom_mpolyun(qrBbar, qBbarn, m, lgctx);
                fq_nmod_mpolyu_mul_mpoly(qT, Bevals + m + 1, gammaevals + m + 1, lgctx);
                if (fq_nmod_mpolyuu_divides(qrG, qT, qrBbar, 2, lgctx))
                {
                    fq_nmod_mpolyu_mul_mpoly(qT, Aevals + m + 1, gammaevals + m + 1, lgctx);
                    if (fq_nmod_mpolyuu_divides(qrAbar, qT, qrG, 2, lgctx))
                    {
                        fq_nmod_mpolyu_swap(qG, qrG, lgctx);
                        fq_nmod_mpolyu_swap(qAbar, qrAbar, lgctx);
                        fq_nmod_mpolyu_swap(qBbar, qrBbar, lgctx);
                        break;
                    }
                }
            }

            if (n_fq_poly_degree(modulus) > gammadegs[m] + Adegs[m] &&
                n_fq_poly_degree(modulus) > gammadegs[m] + Bdegs[m])
            {
                goto choose_alphas;
            }

            FLINT_ASSERT(n_fq_equal_fq_nmod(alphapow->coeffs + lgd*1, alphas + m, lgctx->fqctx));
            n_fq_poly_shift_left_scalar_submul(modulus, 1, alphapow->coeffs + lgd*1, lgctx->fqctx);
        }
    }

    m = nvars - 1;
    {
        /* G, Abar, Bbar are in Fq/alpha(gen(m-1))[gen(0), ..., gen(m - 1)] */
        nmod_mpolyun_interp_lift_lg_mpolyu(Gn, smctx, qG, lgctx);
        nmod_mpolyun_interp_lift_lg_mpolyu(Abarn, smctx, qAbar, lgctx);
        nmod_mpolyun_interp_lift_lg_mpolyu(Bbarn, smctx, qBbar, lgctx);

        if (rGdegs == NULL)
        {
            use = USE_G | USE_ABAR | USE_BBAR;
        }
        else
        {
            slong numABgamma = gamma->length + 
                               nmod_mpolyu_total_length(A) + 
                               nmod_mpolyu_total_length(B);            

            use = fq_nmod_mpoly_gcd_get_use(rGdegs[m], Adegs[m], Bdegs[m],
                   gammadegs[m], degxAB, degyAB, numABgamma, qG, qAbar, qBbar);

        }

        n_poly_set_nmod_poly(modulus, lgctx->fqctx->modulus);

        while (1)
        {
choose_alpha_last:

            lgd++;
            if (lgd > 10000)
            {
                /* ran out of primes */
                success = 0;
                goto cleanup;
            }

            fq_nmod_mpoly_ctx_change_modulus(lgctx, lgd);

            nmod_mpolyun_interp_reduce_lg_mpolyu(Aevals + m, An, lgctx, smctx);
            nmod_mpolyun_interp_reduce_lg_mpolyu(Bevals + m, Bn, lgctx, smctx);
            nmod_mpolyn_interp_reduce_lg_mpoly(gammaevals + m, gamman, lgctx, smctx);
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
                this_images = fq_nmod_mpolyu_set_zip_form(HG, MG, qG, betas, m, lgctx);
                req_zip_images = FLINT_MAX(req_zip_images, this_images);
            }
            if (use & USE_ABAR)
            {
                this_images =  fq_nmod_mpolyu_set_zip_form(HAbar, MAbar, qAbar, betas, m, lgctx);
                req_zip_images = FLINT_MAX(req_zip_images, this_images);            
            }
            if (use & USE_BBAR)
            {
                this_images =  fq_nmod_mpolyu_set_zip_form(HBbar, MBbar, qBbar, betas, m, lgctx);
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
                n_poly_fit_length(tmp, lgd);
                n_fq_poly_eval_step(tmp->coeffs, gammaeh, lgctx->fqctx);

                if (_n_fq_is_zero(tmp->coeffs, lgd))
                    goto choose_betas_last;
                if (Aev->length < 1 || n_bpoly_bidegree(Aev) != A->exps[0])
                    goto choose_betas_last;
                if (Bev->length < 1 || n_bpoly_bidegree(Bev) != B->exps[0])
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
                    goto choose_alphas;
                }

                n_fq_bpoly_scalar_mul_n_fq(Gev, tmp->coeffs, lgctx->fqctx);

                if ((use & USE_G) && !n_fq_polyu2_add_zip_must_match(ZG, Gev, cur_zip_image, lgctx->fqctx))
                    goto choose_alphas;
                if ((use & USE_ABAR) && !n_fq_polyu2_add_zip_must_match(ZAbar, Abarev, cur_zip_image, lgctx->fqctx))
                    goto choose_alphas;
                if ((use & USE_BBAR) && !n_fq_polyu2_add_zip_must_match(ZBbar, Bbarev, cur_zip_image, lgctx->fqctx))
                    goto choose_alphas;
            }
            if ((use & USE_G) && qzip_solve(qG, ZG, HG, MG, lgctx) < 1)
                goto choose_alphas;
            if ((use & USE_ABAR) && qzip_solve(qAbar, ZAbar, HAbar, MAbar, lgctx) < 1)
                goto choose_alphas;
            if ((use & USE_BBAR) && qzip_solve(qBbar, ZBbar, HBbar, MBbar, lgctx) < 1)
                goto choose_alphas;

            FLINT_ASSERT(n_fq_poly_degree(modulus) > 0);

            if ((use & USE_G) && !newnmod_mpolyun_interp_mcrt_lg_mpolyu(
                                      &lastdeg, Gn, smctx, modulus, qG, lgctx))
            {
                nmod_mpolyu_cvtfrom_mpolyun(rG, Gn, m, smctx);
                success = nmod_mpolyu_content_mpoly_threaded_pool(cont, rG, smctx, NULL, 0);
                if (!success)
                    goto cleanup;
                nmod_mpolyu_divexact_mpoly_inplace(rG, cont, smctx);
                if (nmod_mpolyuu_divides(rAbar, A, rG, 2, smctx) &&
                    nmod_mpolyuu_divides(rBbar, B, rG, 2, smctx))
                {
                    break;
                }
            }
            if ((use & USE_ABAR) && !newnmod_mpolyun_interp_mcrt_lg_mpolyu(
                                &lastdeg, Abarn, smctx, modulus, qAbar, lgctx))
            {
                nmod_mpolyu_cvtfrom_mpolyun(rAbar, Abarn, m, smctx);
                success = nmod_mpolyu_content_mpoly_threaded_pool(cont, rAbar, smctx, NULL, 0);
                if (!success)
                    goto cleanup;
                nmod_mpolyu_divexact_mpoly_inplace(rAbar, cont, smctx);
                if (nmod_mpolyuu_divides(rG, A, rAbar, 2, smctx) &&
                    nmod_mpolyuu_divides(rBbar, B, rG, 2, smctx))
                {
                    break;
                }
            }
            if ((use & USE_BBAR) && !newnmod_mpolyun_interp_mcrt_lg_mpolyu(
                                &lastdeg, Bbarn, smctx, modulus, qBbar, lgctx))
            {
                nmod_mpolyu_cvtfrom_mpolyun(rBbar, Bbarn, m, smctx);
                success = nmod_mpolyu_content_mpoly_threaded_pool(cont, rBbar, smctx, NULL, 0);
                if (!success)
                    goto cleanup;
                nmod_mpolyu_divexact_mpoly_inplace(rBbar, cont, smctx);
                if (nmod_mpolyuu_divides(rG, B, rBbar, 2, smctx) &&
                    nmod_mpolyuu_divides(rAbar, A, rG, 2, smctx))
                {
                    break;
                }
            }

            if (n_fq_poly_degree(modulus) > gammadegs[m] + Adegs[m] &&
                n_fq_poly_degree(modulus) > gammadegs[m] + Bdegs[m])
            {
                goto choose_alphas;
            }

            {
                n_poly_t h_mock;
                n_poly_mock(h_mock, lgctx->fqctx->modulus);
                n_poly_mod_mul(modulus, modulus, h_mock, smctx->ffinfo->mod);
            }
        }
    }

    success = 1;

cleanup:

    nmod_mpolyun_clear(An, smctx);
    nmod_mpolyun_clear(Bn, smctx);
    nmod_mpolyn_clear(gamman, smctx);

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
    nmod_mpoly_clear(cont, smctx);

    nmod_mpolyu_clear(T, smctx);
    nmod_mpolyu_clear(G, smctx);
    nmod_mpolyu_clear(Abar, smctx);
    nmod_mpolyu_clear(Bbar, smctx);
    nmod_mpolyun_clear(Tn, smctx);
    nmod_mpolyun_clear(Gn, smctx);
    nmod_mpolyun_clear(Abarn, smctx);
    nmod_mpolyun_clear(Bbarn, smctx);

    fq_nmod_mpolyu_clear(qT, lgctx);
    fq_nmod_mpolyu_clear(qG, lgctx);
    fq_nmod_mpolyu_clear(qAbar, lgctx);
    fq_nmod_mpolyu_clear(qBbar, lgctx);
    fq_nmod_mpolyu_clear(qrG, lgctx);
    fq_nmod_mpolyu_clear(qrAbar, lgctx);
    fq_nmod_mpolyu_clear(qrBbar, lgctx);
    fq_nmod_mpolyun_clear(qTn, lgctx);
    fq_nmod_mpolyun_clear(qGn, lgctx);
    fq_nmod_mpolyun_clear(qAbarn, lgctx);
    fq_nmod_mpolyun_clear(qBbarn, lgctx);

    n_fq_polyun_clear(Aeh);
    n_fq_polyun_clear(Beh);
    n_fq_poly_clear(gammaeh);
    n_poly_clear(modulus);
    n_fq_poly_clear(alphapow);
    n_poly_stack_clear(St->poly_stack);
    n_bpoly_stack_clear(St->bpoly_stack);

    fq_nmod_clear(start_alpha, lgctx->fqctx);

    n_poly_clear(tmp);

    flint_free(betas);
    flint_free(alphas);
    flint_randclear(state);

    for (i = 0; i < nvars; i++)
    {
        fq_nmod_mpolyu_clear(Aevals + i, lgctx);
        fq_nmod_mpolyu_clear(Bevals + i, lgctx);
        fq_nmod_mpoly_clear(gammaevals + i, lgctx);
    }
    flint_free(Aevals);
    flint_free(Bevals);
    flint_free(gammaevals);

    fq_nmod_mpoly_ctx_clear(lgctx);
/*
flint_printf("nmod_mpolyuu_gcd_zippel_lgprime returning %d\n", success);
*/
    return success;

gcd_is_trivial:

    nmod_mpolyu_one(rG, smctx);
    nmod_mpolyu_set(rAbar, A, smctx);
    nmod_mpolyu_set(rBbar, B, smctx);

    success = 1;
    
    goto cleanup;
}


void n_fq_polyun_set(n_fq_polyun_t A, const n_fq_polyun_t B, const fq_nmod_ctx_t ctx)
{
    slong i;
    n_polyun_fit_length(A, B->length);
    for (i = 0; i < B->length; i++)
    {
        A->terms[i].exp = B->terms[i].exp;
        n_fq_poly_set(A->terms[i].coeff, B->terms[i].coeff, ctx);
    }
    A->length = B->length;
}


void n_fq_polyu2n_print_pretty(
    const n_polyun_t A,
    const char * var0,
    const char * var1,
    const char * varlast,
    const fq_nmod_ctx_t ctx)
{
    slong i;
    int first = 1;

    for (i = 0; i < A->length; i++)
    {
        if (!first)
            printf(" + ");
        first = 0;
        flint_printf("(");
        n_fq_poly_print_pretty(A->terms[i].coeff, varlast, ctx);
        flint_printf(")*%s^%wu*%s^%wu",
            var0, extract_exp(A->terms[i].exp, 1, 2),
            var1, extract_exp(A->terms[i].exp, 0, 2));
    }

    if (first)
        flint_printf("0");
}



int nmod_mpolyl_gcd_zippel_lgprime(
    nmod_mpoly_t rG, const slong * rGdegs, /* guess at rG degrees, could be NULL */
    nmod_mpoly_t rAbar,
    nmod_mpoly_t rBbar,
    const nmod_mpoly_t A, const slong * Adegs,
    const nmod_mpoly_t B, const slong * Bdegs,
    const nmod_mpoly_t gamma, const slong * gammadegs,
    const nmod_mpoly_ctx_t smctx)
{
    slong lgd;
    int success, use;
    slong i, j, m;
    slong nvars = smctx->minfo->nvars;
    flint_bitcnt_t bits = A->bits;
    fq_nmod_struct * alphas, * betas;
    flint_rand_t state;
    nmod_mpoly_t cont;
    nmod_mpoly_t T, G, Abar, Bbar;
    nmod_mpolyn_t Tn, Gn, Abarn, Bbarn;
    fq_nmod_mpoly_t qT, qG, qAbar, qBbar;
    fq_nmod_mpoly_t qrG, qrAbar, qrBbar;
    fq_nmod_mpolyn_t qTn, qGn, qAbarn, qBbarn;
    n_fq_polyun_t HG, HAbar, HBbar, MG, MAbar, MBbar, ZG, ZAbar, ZBbar;
    n_fq_bpoly_t Aev, Bev, Gev, Abarev, Bbarev;
    const mp_limb_t * gammaev;
    slong lastdeg;
    slong cur_zip_image, req_zip_images, this_length;
    n_polyun_t Aeh_cur, Aeh_inc, Beh_cur, Beh_inc;
    n_fq_poly_t gammaeh_cur, gammaeh_inc;
    n_fq_poly_t alphapow;
    fq_nmod_mpoly_struct * Aevals, * Bevals;
    fq_nmod_mpoly_struct * gammaevals;
    n_poly_t modulus;
    n_poly_bpoly_stack_t St;
    n_poly_t tmp;  /* tmp arithmetic space */
    fq_nmod_t start_alpha;
    ulong GdegboundXY, newdegXY, Abideg, Bbideg;
    slong degxAB, degyAB;
    fq_nmod_mpoly_ctx_t lgctx;
    nmod_mpolyn_t gamman;
    nmod_mpolyn_t An, Bn;

flint_printf("nmod_mpolyl_gcd_zippel_lgprime called  nvars = %wd\n", smctx->minfo->nvars);

    FLINT_ASSERT(bits <= FLINT_BITS);
    FLINT_ASSERT(bits == rG->bits);
    FLINT_ASSERT(bits == rAbar->bits);
    FLINT_ASSERT(bits == rBbar->bits);
    FLINT_ASSERT(bits == A->bits);
    FLINT_ASSERT(bits == B->bits);
    FLINT_ASSERT(bits == gamma->bits);
    FLINT_ASSERT(A->length > 0);
    FLINT_ASSERT(B->length > 0);
    FLINT_ASSERT(gamma->length > 0);

#if FLINT_WANT_ASSERT
    {
        slong * tmp_degs = FLINT_ARRAY_ALLOC(nvars, slong);

        nmod_mpoly_degrees_si(tmp_degs, A, smctx);
        for (j = 0; j < nvars; j++)
            FLINT_ASSERT(tmp_degs[j] == Adegs[j]);

        nmod_mpoly_degrees_si(tmp_degs, B, smctx);
        for (j = 0; j < nvars; j++)
            FLINT_ASSERT(tmp_degs[j] == Bdegs[j]);

        nmod_mpoly_degrees_si(tmp_degs, gamma, smctx);
        for (j = 0; j < nvars; j++)
            FLINT_ASSERT(tmp_degs[j] == gammadegs[j]);

        flint_free(tmp_degs);
    }
#endif

    flint_randinit(state);

    lgd = WORD(20)/(FLINT_BIT_COUNT(smctx->ffinfo->mod.n));
    lgd = FLINT_MAX(WORD(2), lgd);
    fq_nmod_mpoly_ctx_init_deg(lgctx, nvars, ORD_LEX, smctx->ffinfo->mod.n, lgd);
    n_poly_init2(tmp, lgd);
    n_poly_init2(alphapow, 2*lgd);

    fq_nmod_init(start_alpha, lgctx->fqctx);
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
    nmod_mpoly_init3(cont, 1, bits, smctx);

    nmod_mpoly_init3(T, 1, bits, smctx);
    nmod_mpoly_init3(G, 1, bits, smctx);
    nmod_mpoly_init3(Abar, 1, bits, smctx);
    nmod_mpoly_init3(Bbar, 1, bits, smctx);
    nmod_mpolyn_init(Tn, bits, smctx);
    nmod_mpolyn_init(Gn, bits, smctx);
    nmod_mpolyn_init(Abarn, bits, smctx);
    nmod_mpolyn_init(Bbarn, bits, smctx);

    fq_nmod_mpoly_init3(qT, 1, bits, lgctx);
    fq_nmod_mpoly_init3(qG, 1, bits, lgctx);
    fq_nmod_mpoly_init3(qAbar, 1, bits, lgctx);
    fq_nmod_mpoly_init3(qBbar, 1, bits, lgctx);
    fq_nmod_mpoly_init3(qrG, 1, bits, lgctx);
    fq_nmod_mpoly_init3(qrAbar, 1, bits, lgctx);
    fq_nmod_mpoly_init3(qrBbar, 1, bits, lgctx);
    fq_nmod_mpolyn_init(qTn, bits, lgctx);
    fq_nmod_mpolyn_init(qGn, bits, lgctx);
    fq_nmod_mpolyn_init(qAbarn, bits, lgctx);
    fq_nmod_mpolyn_init(qBbarn, bits, lgctx);

    n_fq_polyun_init(Aeh_cur);
    n_fq_polyun_init(Aeh_inc);
    n_fq_polyun_init(Beh_cur);
    n_fq_polyun_init(Beh_inc);
    n_fq_poly_init(gammaeh_cur);
    n_fq_poly_init(gammaeh_inc);

    n_poly_init(modulus);

    n_poly_stack_init(St->poly_stack);
    n_bpoly_stack_init(St->bpoly_stack);

    nmod_mpolyn_init(An, bits, smctx);
    nmod_mpolyn_init(Bn, bits, smctx);
    nmod_mpolyn_init(gamman, bits, smctx);

    /* alphas[nvars - 1] not used - it is replaced lgctx->fqctx->modulus */
    betas = FLINT_ARRAY_ALLOC(nvars, fq_nmod_struct);
    alphas = FLINT_ARRAY_ALLOC(nvars, fq_nmod_struct);
    for (i = 0; i < nvars; i++)
    {
        fq_nmod_init(betas + i, lgctx->fqctx);
        fq_nmod_init(alphas + i, lgctx->fqctx);
    }

    /* Aevals[nvars] does not exist - it is replaced by An */
    Aevals = FLINT_ARRAY_ALLOC(nvars, fq_nmod_mpoly_struct);
    Bevals = FLINT_ARRAY_ALLOC(nvars, fq_nmod_mpoly_struct);
    gammaevals = FLINT_ARRAY_ALLOC(nvars, fq_nmod_mpoly_struct);
    for (i = 0; i < nvars; i++)
    {
        fq_nmod_mpoly_init3(Aevals + i, 0, bits, lgctx);
        fq_nmod_mpoly_init3(Bevals + i, 0, bits, lgctx);
        fq_nmod_mpoly_init3(gammaevals + i, 0, bits, lgctx);
    }
    nmod_mpoly_cvtto_mpolyn(An, A, nvars - 1, smctx);
    nmod_mpoly_cvtto_mpolyn(Bn, B, nvars - 1, smctx);
    nmod_mpoly_cvtto_mpolyn(gamman, gamma, nvars - 1, smctx);

    Abideg = _nmod_mpoly_bidegree(A, smctx);
    Bbideg = _nmod_mpoly_bidegree(B, smctx);

    degxAB = FLINT_MAX(Adegs[0], Bdegs[0]);
    degyAB = FLINT_MAX(Adegs[1], Bdegs[1]);

    GdegboundXY = pack_exp2(FLINT_MIN(Adegs[0], Bdegs[0]),
                            FLINT_MIN(Adegs[1], Bdegs[1]));
    if (GdegboundXY == 0)
        goto gcd_is_trivial;

    goto got_alpha_last;

increase_degree:

choose_alphas:

    /* TODO: don't necessarily increase degree here */
    lgd++;
    if (lgd > 10000)
    {
        /* ran out of primes */
        success = 0;
        goto cleanup;
    }

    fq_nmod_mpoly_ctx_change_modulus(lgctx, lgd);
    n_poly_fit_length(tmp, lgd);
    n_poly_fit_length(alphapow, 2*lgd);

got_alpha_last:

    for (i = 2; i < nvars - 1; i++)
    {
        fq_nmod_rand(alphas + i, state, lgctx->fqctx);
        if (fq_nmod_is_zero(alphas + i, lgctx->fqctx))
            fq_nmod_one(alphas + i, lgctx->fqctx);
    }

flint_printf("---- got new alphas ----\n");

fq_nmod_ctx_print(lgctx->fqctx);


    i = nvars - 1;
    nmod_mpolyn_interp_reduce_lg_mpoly(Aevals + i, An, lgctx, smctx);
    nmod_mpolyn_interp_reduce_lg_mpoly(Bevals + i, Bn, lgctx, smctx);
    nmod_mpolyn_interp_reduce_lg_mpoly(gammaevals + i, gamman, lgctx, smctx);
    if (fq_nmod_mpoly_is_zero(gammaevals + i, lgctx))
        goto choose_alphas;
    if (Aevals[i].length < 1 || _fq_nmod_mpoly_bidegree(Aevals + i, lgctx) != Abideg)
        goto choose_alphas;
    if (Bevals[i].length < 1 || _fq_nmod_mpoly_bidegree(Bevals + i, lgctx) != Bbideg)
        goto choose_alphas;
    for (i--; i >= 2; i--)
    {
        fq_nmod_mpoly_evaluate_one_fq_nmod(Aevals + i, Aevals + i + 1, i, alphas + i, lgctx);
        fq_nmod_mpoly_repack_bits_inplace(Aevals + i, bits, lgctx);
        fq_nmod_mpoly_evaluate_one_fq_nmod(Bevals + i, Bevals + i + 1, i, alphas + i, lgctx);
        fq_nmod_mpoly_repack_bits_inplace(Bevals + i, bits, lgctx);
        fq_nmod_mpoly_evaluate_one_fq_nmod(gammaevals + i, gammaevals + i + 1, i, alphas + i, lgctx);
        fq_nmod_mpoly_repack_bits_inplace(gammaevals + i, bits, lgctx);
        if (fq_nmod_mpoly_is_zero(gammaevals + i, lgctx))
            goto choose_alphas;
        if (Aevals[i].length < 1 || _fq_nmod_mpoly_bidegree(Aevals + i, lgctx) != Abideg)
            goto choose_alphas;
        if (Bevals[i].length < 1 || _fq_nmod_mpoly_bidegree(Bevals + i, lgctx) != Bbideg)
            goto choose_alphas;
    }

    m = 2;

flint_printf("m = %wd\n", m);

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
        goto choose_alphas;
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
        nmod_mpolyn_interp_lift_lg_bpoly(&lastdeg, Gn, smctx, Gev, lgctx);
        nmod_mpolyn_interp_lift_lg_bpoly(&lastdeg, Abarn, smctx, Abarev, lgctx);
        nmod_mpolyn_interp_lift_lg_bpoly(&lastdeg, Bbarn, smctx, Bbarev, lgctx);

        n_poly_set_nmod_poly(modulus, lgctx->fqctx->modulus);

        while (1)
        {
        choose_alpha_2_last:
            lgd++;
            if (lgd > 10000)
            {
                /* ran out of primes */
                success = 0;
                goto cleanup;
            }

            fq_nmod_mpoly_ctx_change_modulus(lgctx, lgd);
            n_poly_fit_length(tmp, lgd);
            n_poly_fit_length(alphapow, 2*lgd);

            nmod_mpolyn_interp_reduce_lg_mpoly(Aevals + m, An, lgctx, smctx);
            nmod_mpolyn_interp_reduce_lg_mpoly(Bevals + m, Bn, lgctx, smctx);
            nmod_mpolyn_interp_reduce_lg_mpoly(gammaevals + m, gamman, lgctx, smctx);
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
                goto choose_alphas;
            }

            gammaev = fq_nmod_mpoly_get_nonzero_n_fq(gammaevals + m, lgctx);
            n_fq_bpoly_scalar_mul_n_fq(Gev, gammaev, lgctx->fqctx);

            if ((use & USE_G) && !nmod_mpolyn_interp_crt_lg_bpoly(
                                 &lastdeg, Gn, Tn, modulus, smctx, Gev, lgctx))
            {
                nmod_mpoly_cvtfrom_mpolyn(rG, Gn, m, smctx);
                success = nmod_mpolyl_content(cont, rG, 2, smctx);
                if (!success)
                    goto cleanup;
                nmod_mpoly_divides(rG, rG, cont, smctx);
                nmod_mpoly_repack_bits_inplace(rG, bits, smctx);
                if (nmod_mpoly_divides(rAbar, A, rG, smctx) &&
                    nmod_mpoly_divides(rBbar, B, rG, smctx))
                {
                    nmod_mpoly_repack_bits_inplace(rAbar, bits, smctx);
                    nmod_mpoly_repack_bits_inplace(rBbar, bits, smctx);
                    break;
                }
            }

            if ((use & USE_ABAR) && !nmod_mpolyn_interp_crt_lg_bpoly(
                           &lastdeg, Abarn, Tn, modulus, smctx, Abarev, lgctx))
            {
                nmod_mpoly_cvtfrom_mpolyn(rAbar, Abarn, m, smctx);
                success = nmod_mpolyl_content(cont, rAbar, 2, smctx);
                if (!success)
                    goto cleanup;
                nmod_mpoly_divides(rAbar, rAbar, cont, smctx);
                nmod_mpoly_repack_bits_inplace(rAbar, bits, smctx);
                if (nmod_mpoly_divides(rG, A, rAbar, smctx) &&
                    nmod_mpoly_divides(rBbar, B, rG, smctx))
                {
                    nmod_mpoly_repack_bits_inplace(rG, bits, smctx);
                    nmod_mpoly_repack_bits_inplace(rBbar, bits, smctx);
                    break;
                }
            }

            if ((use & USE_BBAR) && !nmod_mpolyn_interp_crt_lg_bpoly(
                          &lastdeg, Bbarn, Tn, modulus, smctx, Bbarev, lgctx))
            {
                nmod_mpoly_cvtfrom_mpolyn(rBbar, Bbarn, m, smctx);
                success = nmod_mpolyl_content(cont, rBbar, 2, smctx);
                if (!success)
                    goto cleanup;
                nmod_mpoly_divides(rBbar, rBbar, cont, smctx);
                nmod_mpoly_repack_bits_inplace(rBbar, bits, smctx);
                if (nmod_mpoly_divides(rG, B, rBbar, smctx) &&
                    nmod_mpoly_divides(rAbar, A, rG, smctx))
                {
                    nmod_mpoly_repack_bits_inplace(rG, bits, smctx);
                    nmod_mpoly_repack_bits_inplace(rAbar, bits, smctx);
                    break;
                }
            }

            if (n_fq_poly_degree(modulus) > gammadegs[m] + Adegs[m] &&
                n_fq_poly_degree(modulus) > gammadegs[m] + Bdegs[m])
            {
                goto choose_alphas;
            }

            {
                n_poly_t h_mock;
                n_poly_mock(h_mock, lgctx->fqctx->modulus);
                n_poly_mod_mul(modulus, modulus, h_mock, smctx->ffinfo->mod);
            }
        }

        success = 1;
        goto cleanup;
    }

    fq_nmod_mpolyn_interp_lift_sm_bpoly(qGn, Gev, lgctx);
    fq_nmod_mpolyn_interp_lift_sm_bpoly(qAbarn, Abarev, lgctx);
    fq_nmod_mpolyn_interp_lift_sm_bpoly(qBbarn, Bbarev, lgctx);

    FLINT_ASSERT(tmp->alloc >= lgd);
    n_fq_poly_one(modulus, lgctx->fqctx);
    n_fq_set_fq_nmod(tmp->coeffs, alphas + m, lgctx->fqctx);
    n_fq_poly_shift_left_scalar_submul(modulus, 1, tmp->coeffs, lgctx->fqctx);

    fq_nmod_set(start_alpha, alphas + m, lgctx->fqctx);

    while (1)
    {
choose_alpha_2:

        fq_nmod_next_not_zero(alphas + m, lgctx->fqctx);
        if (fq_nmod_equal(alphas + m, start_alpha, lgctx->fqctx))
            goto increase_degree;

flint_printf("alphas[%wd]: ", m);
fq_nmod_print_pretty(alphas + m, lgctx->fqctx);
flint_printf("\n");


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

flint_printf("gammaevals[%wd]: ");
fq_nmod_mpoly_print_pretty(gammaevals + m, NULL, lgctx);
flint_printf("\n");

        if (fq_nmod_mpoly_is_zero(gammaevals + m, lgctx))
            goto choose_alpha_2;
        if (Aevals[m].length < 1 || _fq_nmod_mpoly_bidegree(Aevals + m, lgctx) != Abideg)
            goto choose_alpha_2;
        if (Bevals[m].length < 1 || _fq_nmod_mpoly_bidegree(Bevals + m, lgctx) != Bbideg)
            goto choose_alpha_2;

flint_printf("good evals\n");

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
            goto choose_alphas;
        }

        gammaev = fq_nmod_mpoly_get_nonzero_n_fq(gammaevals + m, lgctx);
        n_fq_bpoly_scalar_mul_n_fq(Gev, gammaev, lgctx->fqctx);

        FLINT_ASSERT(tmp->alloc >= lgd);
        n_fq_poly_eval_pow(tmp->coeffs, modulus, alphapow, lgctx->fqctx);
        n_fq_inv(tmp->coeffs, tmp->coeffs, lgctx->fqctx);
        n_fq_poly_scalar_mul_n_fq(modulus, modulus, tmp->coeffs, lgctx->fqctx);

        if ((use & USE_G) && !fq_nmod_mpolyn_interp_crt_sm_bpoly(
                              &lastdeg, qGn, qTn, Gev, modulus, alphapow, lgctx))
        {
flint_printf("G stabilized\n");

            fq_nmod_mpoly_cvtfrom_mpolyn(qG, qGn, m, lgctx);
            fq_nmod_mpoly_mul(qT, Aevals + m + 1, gammaevals + m + 1, lgctx);
            if (fq_nmod_mpoly_divides(qAbar, qT, qG, lgctx))
            {
                fq_nmod_mpoly_mul(qT, Bevals + m + 1, gammaevals + m + 1, lgctx);
                if (fq_nmod_mpoly_divides(qBbar, qT, qG, lgctx))
                {
                    fq_nmod_mpoly_repack_bits_inplace(qAbar, bits, lgctx);
                    fq_nmod_mpoly_repack_bits_inplace(qBbar, bits, lgctx);
                    break;
                }
            }
        }

        if ((use & USE_ABAR) && !fq_nmod_mpolyn_interp_crt_sm_bpoly(
                      &lastdeg, qAbarn, qTn, Abarev, modulus, alphapow, lgctx))
        {

flint_printf("Abar stabilized\n");


            fq_nmod_mpoly_cvtfrom_mpolyn(qAbar, qAbarn, m, lgctx);
            fq_nmod_mpoly_mul(qT, Aevals + m + 1, gammaevals + m + 1, lgctx);
            if (fq_nmod_mpoly_divides(qG, qT, qAbar, lgctx))
            {
                fq_nmod_mpoly_mul(qT, Bevals + m + 1, gammaevals + m + 1, lgctx);
                if (fq_nmod_mpoly_divides(qBbar, qT, qG, lgctx))
                {
                    fq_nmod_mpoly_repack_bits_inplace(qG, bits, lgctx);
                    fq_nmod_mpoly_repack_bits_inplace(qBbar, bits, lgctx);
                    break;
                }
            }
        }

        if ((use & USE_BBAR) && !fq_nmod_mpolyn_interp_crt_sm_bpoly(
                      &lastdeg, qBbarn, qTn, Bbarev, modulus, alphapow, lgctx))
        {

flint_printf("Bbar stabilized\n");

            fq_nmod_mpoly_cvtfrom_mpolyn(qBbar, qBbarn, m, lgctx);
            fq_nmod_mpoly_mul(qT, Bevals + m + 1, gammaevals + m + 1, lgctx);
            if (fq_nmod_mpoly_divides(qG, qT, qBbar, lgctx))
            {
                fq_nmod_mpoly_mul(qT, Aevals + m + 1, gammaevals + m + 1, lgctx);
                if (fq_nmod_mpoly_divides(qAbar, qT, qG, lgctx))
                {
                    fq_nmod_mpoly_repack_bits_inplace(qG, bits, lgctx);
                    fq_nmod_mpoly_repack_bits_inplace(qAbar, bits, lgctx);
                    break;
                }
            }
        }

        if (n_fq_poly_degree(modulus) > gammadegs[m] + Adegs[m] &&
            n_fq_poly_degree(modulus) > gammadegs[m] + Bdegs[m])
        {
            goto choose_alphas;
        }

        FLINT_ASSERT(n_fq_equal_fq_nmod(alphapow->coeffs + lgd*1, alphas + m, lgctx->fqctx));
        n_fq_poly_shift_left_scalar_submul(modulus, 1, alphapow->coeffs + lgd*1, lgctx->fqctx);
    }

    for (m = 3; m < nvars - 1; m++)
    {
flint_printf("m = %wd\n", m);

        /* qG, qAbar, qBbar are in Fq[gen(0), ..., gen(m - 1)] */
        fq_nmod_mpolyn_interp_lift_sm_mpoly(qGn, qG, lgctx);
        fq_nmod_mpolyn_interp_lift_sm_mpoly(qAbarn, qAbar, lgctx);
        fq_nmod_mpolyn_interp_lift_sm_mpoly(qBbarn, qBbar, lgctx);

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

        /* only beta[0], beta[1], ..., beta[m - 1] will be used */
        for (i = 0; i < nvars; i++)
        {
            fq_nmod_rand(betas + i, state, lgctx->fqctx);
            if (fq_nmod_is_zero(betas + i, lgctx->fqctx))
                fq_nmod_one(betas + i, lgctx->fqctx);

flint_printf("betas[%wd]: ", i);
fq_nmod_print_pretty(betas + i, lgctx->fqctx);
flint_printf("\n");

        }

        /* TODO: ditto */
        req_zip_images = 2;

            this_length = fq_nmod_mpoly_set_zip_form2(HG, MG, qG, betas + 2, m, lgctx);
            req_zip_images = FLINT_MAX(req_zip_images, this_length);

            this_length =  fq_nmod_mpoly_set_zip_form2(HAbar, MAbar, qAbar, betas + 2, m, lgctx);
            req_zip_images = FLINT_MAX(req_zip_images, this_length);            

            this_length =  fq_nmod_mpoly_set_zip_form2(HBbar, MBbar, qBbar, betas + 2, m, lgctx);
            req_zip_images = FLINT_MAX(req_zip_images, this_length);            

        if (use == 0)
        {
            this_length = gammaevals[m + 1].length + Aevals[m + 1].length +
                                                     Bevals[m + 1].length;

            use = nmod_mpoly_gcd_get_use_new(rGdegs[m], Adegs[m], Bdegs[m],
                   gammadegs[m], degxAB, degyAB, this_length, HG, HAbar, HBbar);
        }

        while (1)
        {
        choose_alpha_m:

            fq_nmod_next_not_zero(alphas + m, lgctx->fqctx);
            if (fq_nmod_equal(alphas + m, start_alpha, lgctx->fqctx))
                goto increase_degree;
/*
flint_printf("alphas[%wd]: ", m);
fq_nmod_print_pretty(alphas + m, lgctx->fqctx);
flint_printf("\n");
usleep(100000);
*/


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

            fq_nmod_mpoly_get_monomial_evals2(Aeh_inc, Aevals + m, betas + 2, m, lgctx);
            fq_nmod_mpoly_get_monomial_evals2(Beh_inc, Bevals + m, betas + 2, m, lgctx);
            fq_nmod_mpoly_get_monomial_evals(gammaeh_inc, gammaevals + m, betas + 2, 2, m, lgctx);

            n_fq_polyun_set(Aeh_cur, Aeh_inc, lgctx->fqctx);
            n_fq_polyun_set(Beh_cur, Beh_inc, lgctx->fqctx);
            n_fq_poly_set(gammaeh_cur, gammaeh_inc, lgctx->fqctx);
/*
flint_printf("Aevals[%wd]: ",m);
fq_nmod_mpoly_print_pretty(Aevals + m, NULL, lgctx);
flint_printf("\n");

flint_printf("Bevals[%wd]: ",m);
fq_nmod_mpoly_print_pretty(Bevals + m, NULL, lgctx);
flint_printf("\n");

flint_printf("gammaevals[%wd]: ",m);
fq_nmod_mpoly_print_pretty(gammaevals + m, NULL, lgctx);
flint_printf("\n");


flint_printf("starting zip loop use = %d\n", use);
*/

            qzip_start(ZG, HG, req_zip_images, lgctx->fqctx);
            qzip_start(ZAbar, HAbar, req_zip_images, lgctx->fqctx);
            qzip_start(ZBbar, HBbar, req_zip_images, lgctx->fqctx);
            for (cur_zip_image = 0; cur_zip_image < req_zip_images; cur_zip_image++)
            {
                n_fq_bpoly_eval_step_new(Aev, Aeh_cur, Aeh_inc, Aevals + m, lgctx->fqctx);
                n_fq_bpoly_eval_step_new(Bev, Beh_cur, Beh_inc, Bevals + m, lgctx->fqctx);
                n_poly_fit_length(tmp, lgd);
                n_fq_poly_eval_step_new(tmp->coeffs, gammaeh_cur, gammaeh_inc, gammaevals + m, lgctx->fqctx);
/*
flint_printf("Aev: ");
n_fq_bpoly_print_pretty(Aev, "x1", "x2", lgctx->fqctx);
flint_printf("\n");

flint_printf("Bev: ");
n_fq_bpoly_print_pretty(Bev, "x1", "x2", lgctx->fqctx);
flint_printf("\n");

flint_printf("gammaev: ");
n_fq_print_pretty(tmp->coeffs, lgctx->fqctx);
flint_printf("\n");
*/



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
                    goto choose_alphas;
                }

                n_fq_bpoly_scalar_mul_n_fq(Gev, tmp->coeffs, lgctx->fqctx);
/*
flint_printf("Gev: ");
n_fq_bpoly_print_pretty(Gev, "x0", "x1", lgctx->fqctx);
flint_printf("\n");
*/

                if ((use & USE_G) && !n_fq_polyu2_add_zip_must_match(ZG, Gev, cur_zip_image, lgctx->fqctx))
                    goto choose_alphas;
                if ((use & USE_ABAR) && !n_fq_polyu2_add_zip_must_match(ZAbar, Abarev, cur_zip_image, lgctx->fqctx))
                    goto choose_alphas;
                if ((use & USE_BBAR) && !n_fq_polyu2_add_zip_must_match(ZBbar, Bbarev, cur_zip_image, lgctx->fqctx))
                    goto choose_alphas;
            }

            if ((use & USE_G) && qzip_solvel(qG, ZG, HG, MG, lgctx) < 1)
                goto choose_alphas;
            if ((use & USE_ABAR) && qzip_solvel(qAbar, ZAbar, HAbar, MAbar, lgctx) < 1)
                goto choose_alphas;
            if ((use & USE_BBAR) && qzip_solvel(qBbar, ZBbar, HBbar, MBbar, lgctx) < 1)
                goto choose_alphas;

            FLINT_ASSERT(n_fq_poly_degree(modulus) > 0);

            FLINT_ASSERT(tmp->alloc >= lgd);
            n_fq_poly_eval_pow(tmp->coeffs, modulus, alphapow, lgctx->fqctx);
            n_fq_inv(tmp->coeffs, tmp->coeffs, lgctx->fqctx);
            n_fq_poly_scalar_mul_n_fq(modulus, modulus, tmp->coeffs, lgctx->fqctx);

            if ((use & USE_G) && !fq_nmod_mpolyn_interp_mcrt_sm_mpoly(
                                  &lastdeg, qGn, qG, modulus, alphapow, lgctx))
            {
                fq_nmod_mpoly_cvtfrom_mpolyn(qrG, qGn, m, lgctx);
                fq_nmod_mpoly_mul(qT, Aevals + m + 1, gammaevals + m + 1, lgctx);
                if (fq_nmod_mpoly_divides(qrAbar, qT, qrG, lgctx))
                {
                    fq_nmod_mpoly_mul(qT, Bevals + m + 1, gammaevals + m + 1, lgctx);
                    if (fq_nmod_mpoly_divides(qrBbar, qT, qrG, lgctx))
                    {
                        fq_nmod_mpoly_repack_bits_inplace(qrAbar, bits, lgctx);
                        fq_nmod_mpoly_repack_bits_inplace(qrBbar, bits, lgctx);
                        fq_nmod_mpoly_swap(qG, qrG, lgctx);
                        fq_nmod_mpoly_swap(qAbar, qrAbar, lgctx);
                        fq_nmod_mpoly_swap(qBbar, qrBbar, lgctx);
                        break;
                    }
                }
            }

            if ((use & USE_ABAR) && !fq_nmod_mpolyn_interp_mcrt_sm_mpoly(
                              &lastdeg, qAbarn, qAbar, modulus, alphapow, lgctx))
            {
                fq_nmod_mpoly_cvtfrom_mpolyn(qrAbar, qAbarn, m, lgctx);
                fq_nmod_mpoly_mul(qT, Aevals + m + 1, gammaevals + m + 1, lgctx);
                if (fq_nmod_mpoly_divides(qrG, qT, qrAbar, lgctx))
                {
                    fq_nmod_mpoly_mul(qT, Bevals + m + 1, gammaevals + m + 1, lgctx);
                    if (fq_nmod_mpoly_divides(qrBbar, qT, qrG, lgctx))
                    {
                        fq_nmod_mpoly_repack_bits_inplace(qrG, bits, lgctx);
                        fq_nmod_mpoly_repack_bits_inplace(qrBbar, bits, lgctx);
                        fq_nmod_mpoly_swap(qG, qrG, lgctx);
                        fq_nmod_mpoly_swap(qAbar, qrAbar, lgctx);
                        fq_nmod_mpoly_swap(qBbar, qrBbar, lgctx);
                        break;
                    }
                }
            }

            if ((use & USE_BBAR) && !fq_nmod_mpolyn_interp_mcrt_sm_mpoly(
                            &lastdeg, qBbarn, qBbar, modulus, alphapow, lgctx))
            {
                fq_nmod_mpoly_cvtfrom_mpolyn(qrBbar, qBbarn, m, lgctx);
                fq_nmod_mpoly_mul(qT, Bevals + m + 1, gammaevals + m + 1, lgctx);
                if (fq_nmod_mpoly_divides(qrG, qT, qrBbar, lgctx))
                {
                    fq_nmod_mpoly_mul(qT, Aevals + m + 1, gammaevals + m + 1, lgctx);
                    if (fq_nmod_mpoly_divides(qrAbar, qT, qrG, lgctx))
                    {
                        fq_nmod_mpoly_repack_bits_inplace(qrG, bits, lgctx);
                        fq_nmod_mpoly_repack_bits_inplace(qrAbar, bits, lgctx);
                        fq_nmod_mpoly_swap(qG, qrG, lgctx);
                        fq_nmod_mpoly_swap(qAbar, qrAbar, lgctx);
                        fq_nmod_mpoly_swap(qBbar, qrBbar, lgctx);
                        break;
                    }
                }
            }

            if (n_fq_poly_degree(modulus) > gammadegs[m] + Adegs[m] &&
                n_fq_poly_degree(modulus) > gammadegs[m] + Bdegs[m])
            {
                goto choose_alphas;
            }

            FLINT_ASSERT(n_fq_equal_fq_nmod(alphapow->coeffs + lgd*1, alphas + m, lgctx->fqctx));
            n_fq_poly_shift_left_scalar_submul(modulus, 1, alphapow->coeffs + lgd*1, lgctx->fqctx);
        }
    }

    m = nvars - 1;
    {
        /* G, Abar, Bbar are in Fq/alpha(gen(m-1))[gen(0), ..., gen(m - 1)] */
        nmod_mpolyn_interp_lift_lg_mpoly(Gn, smctx, qG, lgctx);
        nmod_mpolyn_interp_lift_lg_mpoly(Abarn, smctx, qAbar, lgctx);
        nmod_mpolyn_interp_lift_lg_mpoly(Bbarn, smctx, qBbar, lgctx);

        if (rGdegs == NULL)
            use = USE_G | USE_ABAR | USE_BBAR;
        else
            use = 0;

        n_poly_set_nmod_poly(modulus, lgctx->fqctx->modulus);

        while (1)
        {
        choose_alpha_last:
/*
flint_printf("choosing alpha[%wd] last\n", m);
*/
            lgd++;
            if (lgd > 10000)
            {
                /* ran out of primes */
                success = 0;
                goto cleanup;
            }

            fq_nmod_mpoly_ctx_change_modulus(lgctx, lgd);

            nmod_mpolyn_interp_reduce_lg_mpoly(Aevals + m, An, lgctx, smctx);
            nmod_mpolyn_interp_reduce_lg_mpoly(Bevals + m, Bn, lgctx, smctx);
            nmod_mpolyn_interp_reduce_lg_mpoly(gammaevals + m, gamman, lgctx, smctx);
            if (fq_nmod_mpoly_is_zero(gammaevals + m, lgctx))
                goto choose_alpha_last;
            if (Aevals[m].length < 1 || _fq_nmod_mpoly_bidegree(Aevals + m, lgctx) != Abideg)
                goto choose_alpha_last;
            if (Bevals[m].length < 1 || _fq_nmod_mpoly_bidegree(Bevals + m, lgctx) != Bbideg)
                goto choose_alpha_last;

        choose_betas_last:

            /* only beta[2], beta[1], ..., beta[m - 1] will be used */
            for (i = 2; i < nvars; i++)
            {
                fq_nmod_rand(betas + i, state, lgctx->fqctx);
                if (fq_nmod_is_zero(betas + i, lgctx->fqctx))
                    fq_nmod_one(betas + i, lgctx->fqctx);
            }

            /* TODO ditto */

            req_zip_images = 1;

            if (use == 0)
            {
                this_length = fq_nmod_mpoly_set_zip_form2(HG, MG, qG, betas + 2, m, lgctx);
                req_zip_images = FLINT_MAX(req_zip_images, this_length);

                this_length =  fq_nmod_mpoly_set_zip_form2(HAbar, MAbar, qAbar, betas + 2, m, lgctx);
                req_zip_images = FLINT_MAX(req_zip_images, this_length);

                this_length =  fq_nmod_mpoly_set_zip_form2(HBbar, MBbar, qBbar, betas + 2, m, lgctx);
                req_zip_images = FLINT_MAX(req_zip_images, this_length);

                this_length = gamma->length + A->length + B->length;
                use = nmod_mpoly_gcd_get_use_new(rGdegs[m], Adegs[m], Bdegs[m],
                       gammadegs[m], degxAB, degyAB, this_length, HG, HAbar, HBbar);
            }
            else
            {
                if (use & USE_G)
                {
                    this_length = fq_nmod_mpoly_set_zip_form2(HG, MG, qG, betas + 2, m, lgctx);
                    req_zip_images = FLINT_MAX(req_zip_images, this_length);
                }
                if (use & USE_ABAR)
                {
                    this_length =  fq_nmod_mpoly_set_zip_form2(HAbar, MAbar, qAbar, betas + 2, m, lgctx);
                    req_zip_images = FLINT_MAX(req_zip_images, this_length);
                }
                if (use & USE_BBAR)
                {
                    this_length =  fq_nmod_mpoly_set_zip_form2(HBbar, MBbar, qBbar, betas + 2, m, lgctx);
                    req_zip_images = FLINT_MAX(req_zip_images, this_length);            
                }
            }

            fq_nmod_mpoly_get_monomial_evals2(Aeh_inc, Aevals + m, betas + 2, m, lgctx);
            fq_nmod_mpoly_get_monomial_evals2(Beh_inc, Bevals + m, betas + 2, m, lgctx);
            fq_nmod_mpoly_get_monomial_evals(gammaeh_inc, gammaevals + m, betas + 2, 2, m, lgctx);

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
                n_poly_fit_length(tmp, lgd);
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
                    goto choose_alphas;
                }

                n_fq_bpoly_scalar_mul_n_fq(Gev, tmp->coeffs, lgctx->fqctx);

                if ((use & USE_G) && !n_fq_polyu2_add_zip_must_match(ZG, Gev, cur_zip_image, lgctx->fqctx))
                    goto choose_alphas;
                if ((use & USE_ABAR) && !n_fq_polyu2_add_zip_must_match(ZAbar, Abarev, cur_zip_image, lgctx->fqctx))
                    goto choose_alphas;
                if ((use & USE_BBAR) && !n_fq_polyu2_add_zip_must_match(ZBbar, Bbarev, cur_zip_image, lgctx->fqctx))
                    goto choose_alphas;
            }

            if ((use & USE_G) && qzip_solvel(qG, ZG, HG, MG, lgctx) < 1)
                goto choose_alphas;
            if ((use & USE_ABAR) && qzip_solvel(qAbar, ZAbar, HAbar, MAbar, lgctx) < 1)
                goto choose_alphas;
            if ((use & USE_BBAR) && qzip_solvel(qBbar, ZBbar, HBbar, MBbar, lgctx) < 1)
                goto choose_alphas;

            FLINT_ASSERT(n_fq_poly_degree(modulus) > 0);

            n_poly_fit_length(tmp, lgd);
            _n_fq_set_n_poly(tmp->coeffs, modulus->coeffs, modulus->length, lgctx->fqctx);
            n_fq_inv(tmp->coeffs, tmp->coeffs, lgctx->fqctx);

            if ((use & USE_G) && !newnmod_mpolyn_interp_mcrt_lg_mpoly(
                         &lastdeg, Gn, smctx, modulus, tmp->coeffs, qG, lgctx))
            {
                nmod_mpoly_cvtfrom_mpolyn(rG, Gn, m, smctx);
                success = nmod_mpolyl_content(cont, rG, 2, smctx);
                if (!success)
                    goto cleanup;
                nmod_mpoly_divides(rG, rG, cont, smctx);
                nmod_mpoly_repack_bits_inplace(rG, bits, smctx);
                if (nmod_mpoly_divides(rAbar, A, rG, smctx) &&
                    nmod_mpoly_divides(rBbar, B, rG,  smctx))
                {
                    nmod_mpoly_repack_bits_inplace(rAbar, bits, smctx);
                    nmod_mpoly_repack_bits_inplace(rBbar, bits, smctx);
                    break;
                }
            }
            if ((use & USE_ABAR) && !newnmod_mpolyn_interp_mcrt_lg_mpoly(
                   &lastdeg, Abarn, smctx, modulus, tmp->coeffs, qAbar, lgctx))
            {
                nmod_mpoly_cvtfrom_mpolyn(rAbar, Abarn, m, smctx);
                success = nmod_mpolyl_content(cont, rAbar, 2, smctx);
                if (!success)
                    goto cleanup;
                nmod_mpoly_divides(rAbar, rAbar, cont, smctx);
                nmod_mpoly_repack_bits_inplace(rAbar, bits, smctx);
                if (nmod_mpoly_divides(rG, A, rAbar, smctx) &&
                    nmod_mpoly_divides(rBbar, B, rG, smctx))
                {
                    nmod_mpoly_repack_bits_inplace(rG, bits, smctx);
                    nmod_mpoly_repack_bits_inplace(rBbar, bits, smctx);
                    break;
                }
            }
            if ((use & USE_BBAR) && !newnmod_mpolyn_interp_mcrt_lg_mpoly(
                   &lastdeg, Bbarn, smctx, modulus, tmp->coeffs, qBbar, lgctx))
            {
                nmod_mpoly_cvtfrom_mpolyn(rBbar, Bbarn, m, smctx);
                success = nmod_mpolyl_content(cont, rBbar, 2, smctx);
                if (!success)
                    goto cleanup;
                nmod_mpoly_divides(rBbar, rBbar, cont, smctx);
                nmod_mpoly_repack_bits_inplace(rBbar, bits, smctx);
                if (nmod_mpoly_divides(rG, B, rBbar, smctx) &&
                    nmod_mpoly_divides(rAbar, A, rG, smctx))
                {
                    nmod_mpoly_repack_bits_inplace(rG, bits, smctx);
                    nmod_mpoly_repack_bits_inplace(rAbar, bits, smctx);
                    break;
                }
            }

            if (n_fq_poly_degree(modulus) > gammadegs[m] + Adegs[m] &&
                n_fq_poly_degree(modulus) > gammadegs[m] + Bdegs[m])
            {
                goto choose_alphas;
            }

            {
                n_poly_t h_mock;
                n_poly_mock(h_mock, lgctx->fqctx->modulus);
                n_poly_mod_mul(modulus, modulus, h_mock, smctx->ffinfo->mod);
            }
        }
    }

    success = 1;

cleanup:

    nmod_mpolyn_clear(An, smctx);
    nmod_mpolyn_clear(Bn, smctx);
    nmod_mpolyn_clear(gamman, smctx);

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
    nmod_mpoly_clear(cont, smctx);
    nmod_mpoly_clear(T, smctx);
    nmod_mpoly_clear(G, smctx);
    nmod_mpoly_clear(Abar, smctx);
    nmod_mpoly_clear(Bbar, smctx);
    nmod_mpolyn_clear(Tn, smctx);
    nmod_mpolyn_clear(Gn, smctx);
    nmod_mpolyn_clear(Abarn, smctx);
    nmod_mpolyn_clear(Bbarn, smctx);
    fq_nmod_mpoly_clear(qT, lgctx);
    fq_nmod_mpoly_clear(qG, lgctx);
    fq_nmod_mpoly_clear(qAbar, lgctx);
    fq_nmod_mpoly_clear(qBbar, lgctx);
    fq_nmod_mpoly_clear(qrG, lgctx);
    fq_nmod_mpoly_clear(qrAbar, lgctx);
    fq_nmod_mpoly_clear(qrBbar, lgctx);
    fq_nmod_mpolyn_clear(qTn, lgctx);
    fq_nmod_mpolyn_clear(qGn, lgctx);
    fq_nmod_mpolyn_clear(qAbarn, lgctx);
    fq_nmod_mpolyn_clear(qBbarn, lgctx);
    n_fq_polyun_clear(Aeh_cur);
    n_fq_polyun_clear(Aeh_inc);
    n_fq_polyun_clear(Beh_cur);
    n_fq_polyun_clear(Beh_inc);
    n_fq_poly_clear(gammaeh_cur);
    n_fq_poly_clear(gammaeh_inc);
    n_poly_clear(modulus);
    n_poly_stack_clear(St->poly_stack);
    n_bpoly_stack_clear(St->bpoly_stack);

    fq_nmod_clear(start_alpha, lgctx->fqctx);

    n_poly_clear(tmp);

    flint_free(betas);
    flint_free(alphas);
    flint_randclear(state);

    for (i = 0; i < nvars; i++)
    {
        fq_nmod_mpoly_clear(Aevals + i, lgctx);
        fq_nmod_mpoly_clear(Bevals + i, lgctx);
        fq_nmod_mpoly_clear(gammaevals + i, lgctx);
    }
    flint_free(Aevals);
    flint_free(Bevals);
    flint_free(gammaevals);

    fq_nmod_mpoly_ctx_clear(lgctx);

flint_printf("nmod_mpolyl_gcd_zippel_lgprime returning %d\n", success);
/*
flint_printf("rG: ");
nmod_mpoly_print_pretty(rG, NULL, smctx);
flint_printf("\n");

flint_printf("rAbar: ");
nmod_mpoly_print_pretty(rAbar, NULL, smctx);
flint_printf("\n");

flint_printf("rBbar: ");
nmod_mpoly_print_pretty(rBbar, NULL, smctx);
flint_printf("\n");
*/


    return success;

gcd_is_trivial:

    nmod_mpoly_one(rG, smctx);
    nmod_mpoly_set(rAbar, A, smctx);
    nmod_mpoly_set(rBbar, B, smctx);

    success = 1;
    
    goto cleanup;
}
