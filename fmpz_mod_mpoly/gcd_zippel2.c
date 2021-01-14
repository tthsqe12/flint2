/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mod_mpoly.h"
#include "fmpz_mod_mpoly_factor.h"


void fmpz_mod_poly_eval_pow(
    fmpz_t eval,
    const fmpz_mod_poly_t P,
    fmpz_mod_poly_t alphapow,
    const fmpz_mod_ctx_t ctx)
{
    slong Plen = P->length;
    if (Plen > alphapow->length)
    {
        slong i = alphapow->length;
        FLINT_ASSERT(2 <= i);
        fmpz_mod_poly_fit_length(alphapow, Plen, ctx);
        alphapow->length = Plen;
        for ( ; i < Plen; i++)
            fmpz_mod_mul(alphapow->coeffs + i, alphapow->coeffs + i - 1,
                                                    alphapow->coeffs + 1, ctx);
    }

    _fmpz_mod_vec_dot(eval, P->coeffs, alphapow->coeffs, Plen, ctx);
}


/*
    evaluation at
        gen(start) -> betas[0]
        gen(start+1) -> betas[1]
        ...
        gen(stop-1) -> betas[stop-start-1]

    the other gen are assumed to not appear in A
*/
void _fmpz_mod_mpoly_monomial_evals_cache(
    fmpz_mod_poly_t E,
    const ulong * Aexps,
    flint_bitcnt_t Abits,
    slong Alen,
    const fmpz * betas,
    slong start,
    slong stop,
    const mpoly_ctx_t mctx,
    const fmpz_mod_ctx_t fctx)
{
    slong i, Ai;
    ulong mask = (-UWORD(1)) >> (FLINT_BITS - Abits);
    slong N = mpoly_words_per_exp_sp(Abits, mctx);
    slong * off, * shift;
    fmpz * c;
    fmpz_t tt;
    slong num = stop - start;

    FLINT_ASSERT(Abits <= FLINT_BITS);
    FLINT_ASSERT(Alen > 0);
    FLINT_ASSERT(num > 0);

    fmpz_init(tt);
    off = FLINT_ARRAY_ALLOC(2*num, slong);
    shift = off + num;
    for (i = 0; i < num; i++)
    {
        mpoly_gen_offset_shift_sp(&off[i], &shift[i], i + start, Abits, mctx);
    }

    fmpz_mod_poly_fit_length(E, Alen, fctx);
    E->length = Alen;

    for (Ai = 0; Ai < Alen; Ai++)
    {
        c = E->coeffs + Ai;
        fmpz_one(c);
        for (i = 0; i < num; i++)
        {
            ulong ei = (Aexps[N*Ai + off[i]] >> shift[i]) & mask;
            fmpz_mod_pow_ui(tt, betas + i, ei, fctx);
            fmpz_mod_mul(c, c, tt, fctx);
        }
    }

    flint_free(off);
    fmpz_clear(tt);
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
void _fmpz_mod_mpoly_monomial_evals2_cache(
    fmpz_mod_polyun_t E,
    const ulong * Aexps,
    flint_bitcnt_t Abits,
    slong Alen,
    const fmpz * betas,
    slong m,
    const mpoly_ctx_t mctx,
    const fmpz_mod_ctx_t fctx)
{
    slong i, Ai, Ei;
    ulong e0, e1, e01;
    ulong mask = (-UWORD(1)) >> (FLINT_BITS - Abits);
    slong N = mpoly_words_per_exp_sp(Abits, mctx);
    slong * off, * shift;
    fmpz * c;
    fmpz_t tt;

    FLINT_ASSERT(Abits <= FLINT_BITS);
    FLINT_ASSERT(Alen > 0);
    FLINT_ASSERT(m > 2);

    fmpz_init(tt);
    off = FLINT_ARRAY_ALLOC(2*m, slong);
    shift = off + m;
    for (i = 0; i < m; i++)
    {
        mpoly_gen_offset_shift_sp(&off[i], &shift[i], i, Abits, mctx);
    }

    Ai = 0;
    Ei = 0;

    e0 = (Aexps[N*Ai + off[0]] >> shift[0]) & mask;
    e1 = (Aexps[N*Ai + off[1]] >> shift[1]) & mask;
    e01 = pack_exp2(e0, e1);
    fmpz_mod_polyun_fit_length(E, Ei + 1, fctx);
    E->terms[Ei].exp = e01;
    fmpz_mod_poly_fit_length(E->terms[Ei].coeff, 1, fctx);
    c = E->terms[Ei].coeff->coeffs + 0;
    E->terms[Ei].coeff->length = 1;
    Ei++;
    fmpz_one(c);
    for (i = 2; i < m; i++)
    {
        ulong ei = (Aexps[N*Ai + off[i]] >> shift[i]) & mask;
        fmpz_mod_pow_ui(tt, betas + i - 2, ei, fctx);
        fmpz_mod_mul(c, c, tt, fctx);
    }

    for (Ai++; Ai < Alen; Ai++)
    {
        e0 = (Aexps[N*Ai + off[0]] >> shift[0]) & mask;
        e1 = (Aexps[N*Ai + off[1]] >> shift[1]) & mask;
        e01 = pack_exp2(e0, e1);
        if (e01 == E->terms[Ei-1].exp)
        {
            slong len = E->terms[Ei-1].coeff->length;
            fmpz_mod_poly_fit_length(E->terms[Ei-1].coeff, len + 1, fctx);
            c = E->terms[Ei-1].coeff->coeffs + len;
            E->terms[Ei-1].coeff->length = len + 1;
        }
        else
        {
            fmpz_mod_polyun_fit_length(E, Ei + 1, fctx);
            E->terms[Ei].exp = e01;
            fmpz_mod_poly_fit_length(E->terms[Ei].coeff, 1, fctx);
            c = E->terms[Ei].coeff->coeffs + 0;
            E->terms[Ei].coeff->length = 1;
            Ei++;
        }

        fmpz_one(c);
        for (i = 2; i < m; i++)
        {
            ulong ei = (Aexps[N*Ai + off[i]] >> shift[i]) & mask;
            fmpz_mod_pow_ui(tt, betas + i - 2, ei, fctx);
            fmpz_mod_mul(c, c, tt, fctx);
        }
    }

    E->length = Ei;

    flint_free(off);
    fmpz_clear(tt);

#if FLINT_WANT_ASSERT
    Ai = 0;
    for (i = 0; i < E->length; i++)
        Ai += E->terms[i].coeff->length;
    FLINT_ASSERT(Ai == Alen);
#endif
}

static ulong _fmpz_mod_mpoly_bidegree(
    const fmpz_mod_mpoly_t A,
    const fmpz_mod_mpoly_ctx_t ctx)
{
    FLINT_ASSERT(A->length > 0);
    return _mpoly_bidegree(A->exps, A->bits, ctx->minfo);
}

static void fmpz_mod_polyun_zip_start(
    fmpz_mod_polyun_t Z,
    fmpz_mod_polyun_t H,
    slong req_images,
    const fmpz_mod_ctx_t fctx)
{
    slong j;
    fmpz_mod_polyun_fit_length(Z, H->length, fctx);
    Z->length = H->length;
    for (j = 0; j < H->length; j++)
    {
        Z->terms[j].exp = H->terms[j].exp;
        fmpz_mod_poly_fit_length(Z->terms[j].coeff, req_images, fctx);
        Z->terms[j].coeff->length = 0;
    }
}

static int fmpz_mod_polyun_zip_solve(
    fmpz_mod_mpoly_t A,
    fmpz_mod_polyun_t Z,
    fmpz_mod_polyun_t H,
    fmpz_mod_polyun_t M,
    const fmpz_mod_mpoly_ctx_t ctx)
{
    int success;
    slong Ai, i, n;
    fmpz * Acoeffs = A->coeffs;
    fmpz_mod_poly_t t;

    fmpz_mod_poly_init(t, ctx->ffinfo);

    FLINT_ASSERT(Z->length == H->length);
    FLINT_ASSERT(Z->length == M->length);

    Ai = 0;
    for (i = 0; i < H->length; i++)
    {
        n = H->terms[i].coeff->length;
        FLINT_ASSERT(M->terms[i].coeff->length == n + 1);
        FLINT_ASSERT(Z->terms[i].coeff->length >= n);
        FLINT_ASSERT(Ai + n <= A->length);

        fmpz_mod_poly_fit_length(t, n, ctx->ffinfo);

        success = _fmpz_mod_zip_vand_solve(Acoeffs + Ai,
                         H->terms[i].coeff->coeffs, n,
                         Z->terms[i].coeff->coeffs, Z->terms[i].coeff->length,
                         M->terms[i].coeff->coeffs, t->coeffs, ctx->ffinfo);
        if (success < 1)
        {
            fmpz_mod_poly_clear(t, ctx->ffinfo);
            return success;
        }

        Ai += n;
        FLINT_ASSERT(Ai <= A->length);

    }

    FLINT_ASSERT(Ai == A->length);

    fmpz_mod_poly_clear(t, ctx->ffinfo);
    return 1;
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

int fmpz_mod_mpoly_gcd_get_use_new(
    slong rGdeg,
    slong Adeg,
    slong Bdeg,
    slong gammadeg,
    slong degxAB,
    slong degyAB,
    slong numABgamma,
    const fmpz_mod_polyun_t G,
    const fmpz_mod_polyun_t Abar,
    const fmpz_mod_polyun_t Bbar)
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

void nmod_mpoly_cvtfrom_mpolyn(
    nmod_mpoly_t A,
    const nmod_mpolyn_t B,
    slong var,
    const nmod_mpoly_ctx_t ctx);

void nmod_mpolyn_interp_lift_sm_mpoly(
    nmod_mpolyn_t A,
    const nmod_mpoly_t B,
    const nmod_mpoly_ctx_t ctx);


void fmpz_mod_poly_eval_step_sep(
    fmpz_t eval,
    fmpz_mod_poly_t cur,
    const fmpz_mod_poly_t inc,
    const fmpz_mod_mpoly_t A,
    const fmpz_mod_mpoly_ctx_t ctx)
{
    FLINT_ASSERT(A->length == cur->length);
    FLINT_ASSERT(A->length == inc->length);
    _fmpz_mod_zip_eval_step(eval, cur->coeffs, inc->coeffs, A->coeffs, A->length, ctx->ffinfo);
}

void fmpz_mod_polyun_scalar_mul_fmpz(
    fmpz_mod_polyun_t A,
    const fmpz_t c,
    const fmpz_mod_ctx_t ctx)
{
    slong i;
    for (i = 0; i < A->length; i++)
        fmpz_mod_poly_scalar_mul_fmpz(A->terms[i].coeff, A->terms[i].coeff, c, ctx);
}

void static fmpz_mod_polyu2n_eval_step_sep(
    fmpz_mod_polyun_t E,
    fmpz_mod_polyun_t cur,
    const fmpz_mod_polyun_t inc,
    const fmpz_mod_mpoly_t A,
    const fmpz_mod_mpoly_ctx_t ctx)
{
    slong i, Ai, Ei, e0, e1, this_len;
    fmpz_mod_poly_struct * Ec;
    fmpz_t eval;

    FLINT_ASSERT(cur->length > 0);
    FLINT_ASSERT(cur->length == inc->length);

    fmpz_init(eval);

    e0 = extract_exp(cur->terms[0].exp, 1, 2);
    e1 = extract_exp(cur->terms[0].exp, 0, 2);

    fmpz_mod_polyun_fit_length(E, 4, ctx->ffinfo);
    Ei = 0;
    E->terms[Ei].exp = e1;
    Ec = E->terms[Ei].coeff;
    fmpz_mod_poly_zero(Ec, ctx->ffinfo);

    Ai = 0;
    for (i = 0; i < cur->length; i++)
    {
        this_len = cur->terms[i].coeff->length;
        _fmpz_mod_zip_eval_step(eval, cur->terms[i].coeff->coeffs,
                                  inc->terms[i].coeff->coeffs,
                                  A->coeffs + Ai, this_len, ctx->ffinfo);
        e0 = extract_exp(cur->terms[i].exp, 1, 2);
        e1 = extract_exp(cur->terms[i].exp, 0, 2);

        if (E->terms[Ei].exp != e0)
        {
            fmpz_mod_polyun_fit_length(E, Ei + 2, ctx->ffinfo);
            Ei += !fmpz_mod_poly_is_zero(E->terms[Ei].coeff, ctx->ffinfo);
            E->terms[Ei].exp = e0;
            Ec = E->terms[Ei].coeff;
            fmpz_mod_poly_zero(Ec, ctx->ffinfo);
        }

        fmpz_mod_poly_set_coeff_fmpz(Ec, e1, eval, ctx->ffinfo);
        Ai += this_len;
    }

    FLINT_ASSERT(Ai == A->length);

    Ei += !fmpz_mod_poly_is_zero(E->terms[Ei].coeff, ctx->ffinfo);
    E->length = Ei;

    FLINT_ASSERT(fmpz_mod_polyun_is_canonical(E, ctx->ffinfo));

    fmpz_clear(eval);
}


void fmpz_mod_mpoly_get_polyu1n(
    fmpz_mod_polyun_t A,
    const fmpz_mod_mpoly_t B,
    slong varx,
    slong vary,
    const fmpz_mod_mpoly_ctx_t ctx)
{
    slong j, Ai;
    ulong Bexpx, Bexpy;
    slong Boffx, Bshiftx, Boffy, Bshifty;
    ulong mask = (-UWORD(1)) >> (FLINT_BITS - B->bits);
    slong NB = mpoly_words_per_exp_sp(B->bits, ctx->minfo);
    
    mpoly_gen_offset_shift_sp(&Boffx, &Bshiftx, varx, B->bits, ctx->minfo);
    mpoly_gen_offset_shift_sp(&Boffy, &Bshifty, vary, B->bits, ctx->minfo);

    Ai = -1;
    for (j = 0; j < B->length; j++)
    {
        Bexpx = ((B->exps + NB*j)[Boffx] >> Bshiftx) & mask;
        Bexpy = ((B->exps + NB*j)[Boffy] >> Bshifty) & mask;

        if (Ai < 0 || A->terms[Ai].exp != Bexpx)
        {
            Ai++;
            fmpz_mod_polyun_fit_length(A, Ai + 1, ctx->ffinfo);
            A->terms[Ai].exp = Bexpx;
            fmpz_mod_poly_zero(A->terms[Ai].coeff, ctx->ffinfo);
        }

        fmpz_mod_poly_set_coeff_fmpz(A->terms[Ai].coeff, Bexpy, B->coeffs + j, ctx->ffinfo);
        if (fmpz_mod_poly_is_zero(A->terms[Ai].coeff, ctx->ffinfo))
            Ai--;
    }

    A->length = Ai + 1;

    FLINT_ASSERT(fmpz_mod_polyun_is_canonical(A, ctx->ffinfo));
}


void fmpz_mod_mpolyn_interp_lift_sm_polyu1n(
    fmpz_mod_mpolyn_t F,
    fmpz_mod_polyun_t A,
    const fmpz_mod_mpoly_ctx_t ctx)
{
    slong N = mpoly_words_per_exp_sp(F->bits, ctx->minfo);
    slong i, j, Fi;
    slong off0, shift0, off1, shift1;

    mpoly_gen_offset_shift_sp(&off0, &shift0, 0, F->bits, ctx->minfo);
    mpoly_gen_offset_shift_sp(&off1, &shift1, 1, F->bits, ctx->minfo);

    Fi = 0;
    for (i = 0; i < A->length; i++)
    {
        fmpz_mod_poly_struct * Ai = A->terms[i].coeff;
        ulong e0 = A->terms[i].exp << shift0;

        for (j = Ai->length - 1; j >= 0; j--)
        {
            if (fmpz_is_zero(Ai->coeffs + j))
                continue;

            fmpz_mod_mpolyn_fit_length(F, Fi + 1, ctx);
            mpoly_monomial_zero(F->exps + N*Fi, N);
            (F->exps + N*Fi)[off0] = e0;
            (F->exps + N*Fi)[off1] += (j << shift1);
            fmpz_mod_poly_set_fmpz(F->coeffs + Fi, Ai->coeffs + j, ctx->ffinfo);
            Fi++;
        }
    }

    F->length = Fi;
}

int fmpz_mod_mpolyn_interp_crt_sm_polyu1n(
    slong * lastdeg,
    fmpz_mod_mpolyn_t F,
    fmpz_mod_mpolyn_t T,
    fmpz_mod_polyun_t A,
    fmpz_mod_poly_t modulus,
    fmpz_mod_poly_t alphapow,
    const fmpz_mod_mpoly_ctx_t ctx)
{
    int changed = 0;
    slong N = mpoly_words_per_exp(F->bits, ctx->minfo);
    slong off0, shift0, off1, shift1;
    fmpz_mod_polyun_term_struct * Aterms = A->terms;
    slong Fi, Ti, Ai, ai;
    slong Alen = A->length;
    slong Flen = F->length;
    ulong * Fexps = F->exps;
    fmpz_mod_poly_struct * Fcoeffs = F->coeffs;
    ulong * Texps = T->exps;
    fmpz_mod_poly_struct * Tcoeffs = T->coeffs;
    fmpz_t v;
    ulong Fexpi, mask;

    fmpz_init(v);

    mask = (-UWORD(1)) >> (FLINT_BITS - F->bits);
    mpoly_gen_offset_shift_sp(&off0, &shift0, 0, F->bits, ctx->minfo);
    mpoly_gen_offset_shift_sp(&off1, &shift1, 1, F->bits, ctx->minfo);

    FLINT_ASSERT(T->bits == F->bits);

    *lastdeg = -1;

    Ti = Fi = 0;
    Ai = 0;
    ai = 0;
    if (Ai < Alen)
        ai = fmpz_mod_poly_degree(Aterms[Ai].coeff, ctx->ffinfo);

    while (Fi < Flen || Ai < Alen)
    {
        if (Ti >= T->alloc)
        {
            slong extra = FLINT_MAX(Flen - Fi, Alen - Ai);
            fmpz_mod_mpolyn_fit_length(T, Ti + extra + 1, ctx);
            Tcoeffs = T->coeffs;
            Texps = T->exps;
        }

        if (Fi < Flen)
            Fexpi = pack_exp2(((Fexps + N*Fi)[off0]>>shift0)&mask,
                              ((Fexps + N*Fi)[off1]>>shift1)&mask);
        else
            Fexpi = 0;

        if (Fi < Flen && Ai < Alen && Fexpi == pack_exp2(Aterms[Ai].exp, ai))
        {
            /* F term ok, A term ok */
            mpoly_monomial_set(Texps + N*Ti, Fexps + N*Fi, N);

            fmpz_mod_poly_eval_pow(v, Fcoeffs + Fi, alphapow, ctx->ffinfo);
            fmpz_mod_sub(v, Aterms[Ai].coeff->coeffs + ai, v, ctx->ffinfo);
            changed |= !fmpz_is_zero(v);
            fmpz_mod_poly_scalar_addmul_fmpz_mod(Tcoeffs + Ti,
                                        Fcoeffs + Fi, modulus, v, ctx->ffinfo);
            Fi++;
            do {
                ai--;
            } while (ai >= 0 && fmpz_is_zero(Aterms[Ai].coeff->coeffs + ai));
            if (ai < 0)
            {
                Ai++;
                if (Ai < Alen)
                    ai = fmpz_mod_poly_degree(Aterms[Ai].coeff, ctx->ffinfo);
            }
        }
        else if (Ai < Alen && (Fi >= Flen || Fexpi < pack_exp2(Aterms[Ai].exp, ai)))
        {
            /* F term missing, A term ok */
            mpoly_monomial_zero(Texps + N*Ti, N);
            (Texps + N*Ti)[off0] += (Ai << shift0);
            (Texps + N*Ti)[off1] += (ai << shift1);

            changed = 1;
            fmpz_mod_poly_scalar_mul_fmpz(Tcoeffs + Ti, modulus,
                                   Aterms[Ai].coeff->coeffs + ai, ctx->ffinfo);

            do {
                ai--;
            } while (ai >= 0 && fmpz_is_zero(Aterms[Ai].coeff->coeffs + ai));
            if (ai < 0)
            {
                Ai++;
                if (Ai < Alen)
                    ai = fmpz_mod_poly_degree(Aterms[Ai].coeff, ctx->ffinfo);
            }
        }
        else
        {
            FLINT_ASSERT(Fi < Flen && (Ai >= Alen || Fexpi > pack_exp2(Aterms[Ai].exp, ai)));
            /* F term ok, Aterm missing */
            mpoly_monomial_set(Texps + N*Ti, Fexps + N*Fi, N);

            fmpz_mod_poly_eval_pow(v, Fcoeffs + Fi, alphapow, ctx->ffinfo);
            fmpz_mod_neg(v, v, ctx->ffinfo);
            changed |= !fmpz_is_zero(v);
            fmpz_mod_poly_scalar_addmul_fmpz_mod(Tcoeffs + Ti,
                                        Fcoeffs + Fi, modulus, v, ctx->ffinfo);
            Fi++;
        }

        FLINT_ASSERT(!fmpz_mod_poly_is_zero(Tcoeffs + Ti, ctx->ffinfo));
        *lastdeg = FLINT_MAX(*lastdeg, fmpz_mod_poly_degree(Tcoeffs + Ti, ctx->ffinfo));
        Ti++;
    }

    T->length = Ti;

    if (changed)
        fmpz_mod_mpolyn_swap(T, F, ctx);

    fmpz_clear(v);

    return changed;
}

void fmpz_mod_polyun_content_poly(
    fmpz_mod_poly_t g,
    const fmpz_mod_polyun_t A,
    const fmpz_mod_ctx_t ctx)
{
    slong i;
    fmpz_mod_poly_zero(g, ctx);
    for (i = 0; i < A->length; i++)
        fmpz_mod_poly_gcd(g, g, A->terms[i].coeff, ctx);
}

void fmpz_mod_polyun_divexact_poly(
    fmpz_mod_polyun_t A,
    const fmpz_mod_poly_t g,
    const fmpz_mod_ctx_t ctx)
{
    slong i;
    for (i = 0; i < A->length; i++)
        fmpz_mod_poly_div(A->terms[i].coeff, A->terms[i].coeff, g, ctx);
}

void fmpz_mod_polyun_mul_poly(
    fmpz_mod_polyun_t A,
    const fmpz_mod_poly_t g,
    const fmpz_mod_ctx_t ctx)
{
    slong i;
    for (i = 0; i < A->length; i++)
        fmpz_mod_poly_mul(A->terms[i].coeff, A->terms[i].coeff, g, ctx);
}


slong fmpz_mod_polyun_lastdeg(
    const fmpz_mod_polyun_t A,
    const fmpz_mod_ctx_t ctx)
{
    slong i, len = 0;
    for (i = 0; i < A->length; i++)
        len = FLINT_MAX(len, A->terms[i].coeff->length);
    return len - 1;
}



/*
    E = A(y = alpha)
    A is in R[y][X]
    E is in R[X]
*/
void fmpz_mod_polyu1n_intp_reduce_sm_poly(
    fmpz_mod_poly_t E,
    const fmpz_mod_polyun_t A,
    const fmpz_t alpha,
    const fmpz_mod_ctx_t ctx)
{
    fmpz_t v;
    slong Ai;
    fmpz_mod_polyun_term_struct * Aterms = A->terms;
    slong Alen = A->length;

    fmpz_init(v);
    fmpz_mod_poly_zero(E, ctx);
    for (Ai = 0; Ai < Alen; Ai++)
    {
        fmpz_mod_poly_evaluate_fmpz(v, Aterms[Ai].coeff, alpha, ctx);
        fmpz_mod_poly_set_coeff_fmpz(E, Aterms[Ai].exp, v, ctx);
    }
    fmpz_clear(v);
}

/*
    A = B
    A, B are in R[X]
*/
void fmpz_mod_polyu1n_intp_lift_sm_poly(
    fmpz_mod_polyun_t A,
    const fmpz_mod_poly_t B,
    const fmpz_mod_ctx_t ctx)
{
    slong Bi;
    slong Blen = B->length;
    fmpz * Bcoeff = B->coeffs;
    fmpz_mod_polyun_term_struct * Aterms;
    slong Ai;

    fmpz_mod_polyun_fit_length(A, Blen, ctx);
    Aterms = A->terms;

    Ai = 0;
    for (Bi = Blen - 1; Bi >= 0; Bi--)
    {
        if (fmpz_is_zero(Bcoeff + Bi))
            continue;

        FLINT_ASSERT(Ai < A->alloc);

        fmpz_mod_poly_set_fmpz(Aterms[Ai].coeff, Bcoeff + Bi, ctx);
        Aterms[Ai].exp = Bi;
        Ai++;
    }
    A->length = Ai;
}


/*
    F = F + modulus*(A - F(v = alpha))
    no assumptions about matching monomials
    F is in Fp[X][v]
    A is in Fp[X]
    it is expected that modulus(alpha) == 1
*/
int fmpz_mod_polyu1n_intp_crt_sm_poly(
    slong * lastdeg,
    fmpz_mod_polyun_t F,
    fmpz_mod_polyun_t T,
    const fmpz_mod_poly_t A,
    const fmpz_mod_poly_t modulus,
    const fmpz_t alpha,
    const fmpz_mod_ctx_t ctx)
{
    int changed = 0;
    slong lastlen = 0;
    fmpz_t v;
    slong Fi, Ti, Ai;
    fmpz * Acoeffs = A->coeffs;
    slong Flen = F->length;
    fmpz_mod_polyun_term_struct * Fterms = F->terms;
    fmpz_mod_polyun_term_struct * Tterms;

    Fi = 0;
    Ai = fmpz_mod_poly_degree(A, ctx);

    fmpz_init(v);

    fmpz_mod_polyun_fit_length(T, Flen + Ai + 1, ctx);
    Tterms = T->terms;
    Ti = 0;

    while (Fi < Flen || Ai >= 0)
    {
        FLINT_ASSERT(Ti < T->alloc);

        if (Fi < Flen)
        {
            FLINT_ASSERT(!fmpz_mod_poly_is_zero(Fterms[Fi].coeff, ctx));
            FLINT_ASSERT(fmpz_mod_poly_degree(Fterms[Fi].coeff, ctx) <
                                           fmpz_mod_poly_degree(modulus, ctx));
        }

        if (Ai >= 0)
        {
            FLINT_ASSERT(!fmpz_is_zero(Acoeffs + Ai));
        }

        if (Fi < Flen && Ai >= 0 && Fterms[Fi].exp == Ai)
        {
            /* F term ok, A term ok */
            fmpz_mod_poly_evaluate_fmpz(v, Fterms[Fi].coeff, alpha, ctx);
            fmpz_mod_sub(v, Acoeffs + Ai, v, ctx);
            changed |= !fmpz_is_zero(v);
            fmpz_mod_poly_scalar_addmul_fmpz_mod(Tterms[Ti].coeff,
                                            Fterms[Fi].coeff, modulus, v, ctx);
            Tterms[Ti].exp = Ai;
            Fi++;
            do {
                Ai--;
            } while (Ai >= 0 && fmpz_is_zero(Acoeffs + Ai));
        }
        else if (Fi < Flen && (Ai < 0 || Fterms[Fi].exp > Ai))
        {
            /* F term ok, A term missing */
            fmpz_mod_poly_evaluate_fmpz(v, Fterms[Fi].coeff, alpha, ctx);
            fmpz_mod_neg(v, v, ctx);
            changed |= !fmpz_is_zero(v);
            fmpz_mod_poly_scalar_addmul_fmpz_mod(Tterms[Ti].coeff,
                                            Fterms[Fi].coeff, modulus, v, ctx);
            Tterms[Ti].exp = Fterms[Fi].exp;
            Fi++;
        }
        else if (Ai >= 0 && (Fi >= Flen || Fterms[Fi].exp < Ai))
        {
            /* F term missing, A term ok */
            changed = 1;
            fmpz_mod_poly_scalar_mul_fmpz(Tterms[Ti].coeff,
                                                   modulus, Acoeffs + Ai, ctx);
            Tterms[Ti].exp = Ai;
            do {
                Ai--;
            } while (Ai >= 0 && fmpz_is_zero(Acoeffs + Ai));
        }
        else
        {
            FLINT_ASSERT(0);
        }

        FLINT_ASSERT(!fmpz_mod_poly_is_zero(Tterms[Ti].coeff, ctx));
        lastlen = FLINT_MAX(lastlen, Tterms[Ti].coeff->length);

        Ti++;
    }
    T->length = Ti;

    fmpz_clear(v);

    if (changed)
        fmpz_mod_polyun_swap(T, F);

    *lastdeg = lastlen - 1;
    return changed;
}

void fmpz_mod_polyun_one(fmpz_mod_polyun_t A, const fmpz_mod_ctx_t ctx)
{
    fmpz_mod_polyun_fit_length(A, 1, ctx);
    fmpz_mod_poly_one(A->terms[0].coeff, ctx);
    A->terms[0].exp = 0;
    A->length = 1;
}


int fmpz_mod_polyu1n_gcd_brown_smprime(
    fmpz_mod_polyun_t G,
    fmpz_mod_polyun_t Abar,
    fmpz_mod_polyun_t Bbar,
    fmpz_mod_polyun_t A,
    fmpz_mod_polyun_t B,
    const fmpz_mod_ctx_t ctx,
    fmpz_mod_poly_polyun_stack_t St)
{
    int success;
    slong bound;
    fmpz_t alpha, temp, gammaeval;
    fmpz_mod_poly_t Aeval, Beval, Geval, Abareval, Bbareval;
    fmpz_mod_polyun_t T;
    slong deggamma, ldegG, ldegAbar, ldegBbar, ldegA, ldegB;
    fmpz_mod_poly_t cA, cB, cG, cAbar, cBbar, gamma, r;
    fmpz_mod_poly_t modulus;
#if FLINT_WANT_ASSERT
    fmpz_mod_poly_t leadA, leadB;
#endif

    fmpz_init(alpha);
    fmpz_init(temp);
    fmpz_init(gammaeval);

#if FLINT_WANT_ASSERT
    fmpz_mod_poly_init(leadA, ctx);
    fmpz_mod_poly_init(leadB, ctx);
    fmpz_mod_poly_set(leadA, A->terms[0].coeff, ctx);
    fmpz_mod_poly_set(leadB, B->terms[0].coeff, ctx);
#endif

    fmpz_mod_polyun_init(T, ctx);
    fmpz_mod_poly_init(r, ctx);
    fmpz_mod_poly_init(cA, ctx);
    fmpz_mod_poly_init(cB, ctx);
    fmpz_mod_poly_init(cG, ctx);
    fmpz_mod_poly_init(cAbar, ctx);
    fmpz_mod_poly_init(cBbar, ctx);
    fmpz_mod_poly_init(gamma, ctx);
    fmpz_mod_poly_init(Aeval, ctx);
    fmpz_mod_poly_init(Beval, ctx);
    fmpz_mod_poly_init(Geval, ctx);
    fmpz_mod_poly_init(Abareval, ctx);
    fmpz_mod_poly_init(Bbareval, ctx);
    fmpz_mod_poly_init(modulus, ctx);

    fmpz_mod_polyun_content_poly(cA, A, ctx);
    fmpz_mod_polyun_content_poly(cB, B, ctx);
    fmpz_mod_polyun_divexact_poly(A, cA, ctx);
    fmpz_mod_polyun_divexact_poly(B, cB, ctx);

    fmpz_mod_poly_gcd(cG, cA, cB, ctx);

    fmpz_mod_poly_divrem(cAbar, r, cA, cG, ctx);
    FLINT_ASSERT(fmpz_mod_poly_is_zero(r, ctx));
    fmpz_mod_poly_divrem(cBbar, r, cB, cG, ctx);
    FLINT_ASSERT(fmpz_mod_poly_is_zero(r, ctx));

    fmpz_mod_poly_gcd(gamma, A->terms[0].coeff, B->terms[0].coeff, ctx);

    ldegA = fmpz_mod_polyun_lastdeg(A, ctx);
    ldegB = fmpz_mod_polyun_lastdeg(B, ctx);
    deggamma = fmpz_mod_poly_degree(gamma, ctx);
    bound = 1 + deggamma + FLINT_MAX(ldegA, ldegB);

    fmpz_mod_poly_one(modulus, ctx);

    fmpz_sub_ui(alpha, fmpz_mod_ctx_modulus(ctx), 1);

choose_prime:

    fmpz_sub_ui(alpha, alpha, 1);
    if (fmpz_sgn(alpha) <= 0)
    {
        success = 0;
        goto cleanup;
    }

    /* make sure evaluation point does not kill both lc(A) and lc(B) */
    fmpz_mod_poly_evaluate_fmpz(gammaeval, gamma, alpha, ctx);
    if (fmpz_is_zero(gammaeval))
    {
        goto choose_prime;
    }

    /* evaluation point should kill neither A nor B */
    fmpz_mod_polyu1n_intp_reduce_sm_poly(Aeval, A, alpha, ctx);
    fmpz_mod_polyu1n_intp_reduce_sm_poly(Beval, B, alpha, ctx);
    FLINT_ASSERT(Aeval->length > 0);
    FLINT_ASSERT(Beval->length > 0);

    fmpz_mod_poly_gcd(Geval, Aeval, Beval, ctx);
    fmpz_mod_poly_divrem(Abareval, r, Aeval, Geval, ctx);
    FLINT_ASSERT(fmpz_mod_poly_is_zero(r, ctx));
    fmpz_mod_poly_divrem(Bbareval, r, Beval, Geval, ctx);
    FLINT_ASSERT(fmpz_mod_poly_is_zero(r, ctx));

    FLINT_ASSERT(Geval->length > 0);
    FLINT_ASSERT(Abareval->length > 0);
    FLINT_ASSERT(Bbareval->length > 0);

    if (fmpz_mod_poly_degree(Geval, ctx) == 0)
    {
        fmpz_mod_polyun_one(G, ctx);
        fmpz_mod_polyun_swap(Abar, A);
        fmpz_mod_polyun_swap(Bbar, B);
        goto successful_put_content;
    }

    if (fmpz_mod_poly_degree(modulus, ctx) > 0)
    {
        FLINT_ASSERT(G->length > 0);
        if (fmpz_mod_poly_degree(Geval, ctx) > G->terms[0].exp)
        {
            goto choose_prime;
        }
        else if (fmpz_mod_poly_degree(Geval, ctx) < G->terms[0].exp)
        {
            fmpz_mod_poly_one(modulus, ctx);
        }
    }

    fmpz_mod_poly_scalar_mul_fmpz(Geval, Geval, gammaeval, ctx);

    if (fmpz_mod_poly_degree(modulus, ctx) > 0)
    {
        fmpz_mod_poly_evaluate_fmpz(temp, modulus, alpha, ctx);
        fmpz_mod_inv(temp, temp, ctx);
        fmpz_mod_poly_scalar_mul_fmpz(modulus, modulus, temp, ctx);
        fmpz_mod_polyu1n_intp_crt_sm_poly(&ldegG, G, T, Geval, modulus, alpha, ctx);
        fmpz_mod_polyu1n_intp_crt_sm_poly(&ldegAbar, Abar, T, Abareval, modulus, alpha, ctx);
        fmpz_mod_polyu1n_intp_crt_sm_poly(&ldegBbar, Bbar, T, Bbareval, modulus, alpha, ctx);
    }
    else
    {
        fmpz_mod_polyu1n_intp_lift_sm_poly(G, Geval, ctx);
        fmpz_mod_polyu1n_intp_lift_sm_poly(Abar, Abareval, ctx);
        fmpz_mod_polyu1n_intp_lift_sm_poly(Bbar, Bbareval, ctx);
        ldegG = 0;
        ldegAbar = 0;
        ldegBbar = 0;
    }

    fmpz_mod_neg(temp, alpha, ctx);
    fmpz_mod_poly_shift_left_scalar_addmul_fmpz_mod(modulus, 1, temp, ctx);

    if (fmpz_mod_poly_degree(modulus, ctx) < bound)
    {
        goto choose_prime;
    }


    FLINT_ASSERT(ldegG >= 0);
    FLINT_ASSERT(ldegAbar >= 0);
    FLINT_ASSERT(ldegBbar >= 0);

    if (deggamma + ldegA == ldegG + ldegAbar &&
        deggamma + ldegB == ldegG + ldegBbar )
    {
        goto successful;
    }

    fmpz_mod_poly_one(modulus, ctx);
    goto choose_prime;

successful:

    fmpz_mod_polyun_content_poly(modulus, G, ctx);
    fmpz_mod_polyun_divexact_poly(G, modulus, ctx);
    fmpz_mod_polyun_divexact_poly(Abar, G->terms[0].coeff, ctx);
    fmpz_mod_polyun_divexact_poly(Bbar, G->terms[0].coeff, ctx);

successful_put_content:

    fmpz_mod_polyun_mul_poly(G, cG, ctx);
    fmpz_mod_polyun_mul_poly(Abar, cAbar, ctx);
    fmpz_mod_polyun_mul_poly(Bbar, cBbar, ctx);

    success = 1;

cleanup:

#if FLINT_WANT_ASSERT
    if (success)
    {
        FLINT_ASSERT(fmpz_is_one(fmpz_mod_poly_lead(G->terms[0].coeff, ctx)));
        fmpz_mod_poly_mul(modulus, G->terms[0].coeff, Abar->terms[0].coeff, ctx);
        FLINT_ASSERT(fmpz_mod_poly_equal(modulus, leadA, ctx));
        fmpz_mod_poly_mul(modulus, G->terms[0].coeff, Bbar->terms[0].coeff, ctx);
        FLINT_ASSERT(fmpz_mod_poly_equal(modulus, leadB, ctx));
    }
    fmpz_mod_poly_clear(leadA, ctx);
    fmpz_mod_poly_clear(leadB, ctx);
#endif

    fmpz_mod_poly_clear(r, ctx);
    fmpz_mod_poly_clear(cA, ctx);
    fmpz_mod_poly_clear(cB, ctx);
    fmpz_mod_poly_clear(cG, ctx);
    fmpz_mod_poly_clear(cAbar, ctx);
    fmpz_mod_poly_clear(cBbar, ctx);

    fmpz_mod_poly_clear(gamma, ctx);

    fmpz_mod_poly_clear(Aeval, ctx);
    fmpz_mod_poly_clear(Beval, ctx);
    fmpz_mod_poly_clear(Geval, ctx);
    fmpz_mod_poly_clear(Abareval, ctx);
    fmpz_mod_poly_clear(Bbareval, ctx);

    fmpz_mod_polyun_clear(T, ctx);

    fmpz_mod_poly_clear(modulus, ctx);

    fmpz_clear(alpha);
    fmpz_clear(temp);
    fmpz_clear(gammaeval);

    return success;
}


static void fmpz_mod_mpoly_monomial_evals(
    fmpz_mod_poly_t E,
    const fmpz_mod_mpoly_t A,
    const fmpz * betas,
    slong start,
    slong stop,
    const fmpz_mod_mpoly_ctx_t ctx)
{
    _fmpz_mod_mpoly_monomial_evals_cache(E, A->exps, A->bits, A->length,
                                  betas, start, stop, ctx->minfo, ctx->ffinfo);
}

static void fmpz_mod_mpoly_monomial_evals2(
    fmpz_mod_polyun_t E,
    const fmpz_mod_mpoly_t A,
    const fmpz * betas,
    slong m,
    const fmpz_mod_mpoly_ctx_t ctx)
{
    _fmpz_mod_mpoly_monomial_evals2_cache(E, A->exps, A->bits, A->length, betas, m,
                                                         ctx->minfo, ctx->ffinfo);
}


void fmpz_mod_mpolyn_interp_lift_sm_mpoly(
    fmpz_mod_mpolyn_t A,
    const fmpz_mod_mpoly_t B,
    const fmpz_mod_mpoly_ctx_t ctx)
{
    slong i, N;
    fmpz_mod_poly_struct * Acoeff;
    const fmpz * Bcoeff;
    ulong * Aexp, * Bexp;
    slong Blen;

    FLINT_ASSERT(A->bits == B->bits);

    Blen = B->length;
    fmpz_mod_mpolyn_fit_length(A, Blen, ctx);
    Acoeff = A->coeffs;
    Bcoeff = B->coeffs;
    Aexp = A->exps;
    Bexp = B->exps;

    N = mpoly_words_per_exp(B->bits, ctx->minfo);

    for (i = 0; i < Blen; i++)
    {
        fmpz_mod_poly_set_fmpz(Acoeff + i, Bcoeff + i, ctx->ffinfo);
        mpoly_monomial_set(Aexp + N*i, Bexp + N*i, N);
    }
    A->length = Blen;
}


int fmpz_mod_mpolyn_interp_crt_sm_mpoly(
    slong * lastdeg,
    fmpz_mod_mpolyn_t F,
    fmpz_mod_mpolyn_t T,
    fmpz_mod_mpoly_t A,
    fmpz_mod_poly_t modulus,
    const fmpz_t alpha,
    const fmpz_mod_mpoly_ctx_t ctx)
{
    slong lastlen = 0;
    int changed = 0;
    slong i, j, k;
    slong N;
    fmpz_t v;
    flint_bitcnt_t bits = A->bits;
    slong Flen = F->length, Alen = A->length;
    ulong * Fexp = F->exps, * Aexp = A->exps;
    ulong * Texp;
    fmpz * Acoeff = A->coeffs;
    fmpz_mod_poly_struct * Fcoeff = F->coeffs;
    fmpz_mod_poly_struct * Tcoeff;

    FLINT_ASSERT(F->bits == bits);
    FLINT_ASSERT(ctx->minfo->ord == ORD_LEX);

    fmpz_mod_mpolyn_fit_length(T, Flen + Alen, ctx);
    Texp = T->exps;
    Tcoeff = T->coeffs;

    N = mpoly_words_per_exp(bits, ctx->minfo);

    i = j = k = 0;
    while (i < Flen || j < Alen)
    {
        if (i < Flen && (j >= Alen
                       || mpoly_monomial_gt_nomask(Fexp + N*i, Aexp + N*j, N)))
        {
            FLINT_ASSERT(!fmpz_mod_poly_is_zero(Fcoeff + i, ctx->ffinfo));
            FLINT_ASSERT(fmpz_mod_poly_degree(Fcoeff + i, ctx->ffinfo) <
                                   fmpz_mod_poly_degree(modulus, ctx->ffinfo));

            /* F term ok, A term missing */
            fmpz_mod_poly_evaluate_fmpz(v, Fcoeff + i, alpha, ctx->ffinfo);
            if (!fmpz_is_zero(v))
            {
                changed = 1;
                fmpz_mod_neg(v, v, ctx->ffinfo);
                fmpz_mod_poly_scalar_addmul_fmpz_mod(Tcoeff + k,
                                          Fcoeff + i, modulus, v, ctx->ffinfo);
            }
            else
            {
                fmpz_mod_poly_set(Tcoeff + k, Fcoeff + i, ctx->ffinfo);                
            }

            FLINT_ASSERT(!fmpz_mod_poly_is_zero(Tcoeff + k, ctx->ffinfo));
            lastlen = FLINT_MAX(lastlen, Tcoeff[k].length);

            mpoly_monomial_set(Texp + N*k, Fexp + N*i, N);
            k++;
            i++;
        }
        else if (j < Alen && (i >= Flen
                        || mpoly_monomial_lt_nomask(Fexp + N*i, Aexp + N*j, N)))
        {
            /* F term missing, A term ok */
            FLINT_ASSERT(!fmpz_is_zero(Acoeff + j));

            changed = 1;
            fmpz_mod_poly_scalar_mul_fmpz(Tcoeff + k, modulus, Acoeff + j, ctx->ffinfo);
            lastlen = FLINT_MAX(lastlen, Tcoeff[k].length);
            mpoly_monomial_set(Texp + N*k, Aexp + N*j, N);
            k++;
            j++;
        }
        else if (i < Flen && j < Alen
                            && mpoly_monomial_equal(Fexp + N*i, Aexp + N*j, N))
        {
            FLINT_ASSERT(!fmpz_mod_poly_is_zero(Fcoeff + i, ctx->ffinfo));
            FLINT_ASSERT(fmpz_mod_poly_degree(Fcoeff + i, ctx->ffinfo) <
                                   fmpz_mod_poly_degree(modulus, ctx->ffinfo));

            /* F term ok, A term ok */
            fmpz_mod_poly_evaluate_fmpz(v, Fcoeff + i, alpha, ctx->ffinfo);
            fmpz_mod_sub(v, Acoeff + j, v, ctx->ffinfo);
            if (!fmpz_is_zero(v))
            {
                changed = 1;
                fmpz_mod_poly_scalar_addmul_fmpz_mod(Tcoeff + k,
                                          Fcoeff + i, modulus, v, ctx->ffinfo);
            }
            else
            {
                fmpz_mod_poly_set(Tcoeff + k, Fcoeff + i, ctx->ffinfo);            
            }
            FLINT_ASSERT(!fmpz_mod_poly_is_zero(Tcoeff + k, ctx->ffinfo));
            lastlen = FLINT_MAX(lastlen, Tcoeff[k].length);
            mpoly_monomial_set(Texp + N*k, Aexp + N*j, N);
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

    *lastdeg = lastlen - 1;

    if (changed)
        fmpz_mod_mpolyn_swap(T, F, ctx);

    return changed;
}


/* return the largest degree */
slong fmpz_mod_polyun_product_roots(
    fmpz_mod_polyun_t M,
    const fmpz_mod_polyun_t H,
    const fmpz_mod_ctx_t ctx)
{
    slong i, max_length = 0;

    fmpz_mod_polyun_fit_length(M, H->length, ctx);
    M->length = H->length;
    for (i = 0; i < H->length; i++)
    {
        slong len = H->terms[i].coeff->length;
        M->terms[i].exp = H->terms[i].exp;
        max_length = FLINT_MAX(max_length, len);
        fmpz_mod_poly_product_roots_fmpz_vec(M->terms[i].coeff,
                                          H->terms[i].coeff->coeffs, len, ctx);
    }

    return max_length;
}


int fmpz_mod_polyun_add_zip_must_match(
    fmpz_mod_polyun_t Z,
    const fmpz_mod_polyun_t A,
    slong cur_length)
{
    slong i, Ai, ai;
    slong Alen = A->length;
    fmpz_mod_polyun_term_struct * Zterms = Z->terms;
    const fmpz_mod_polyun_term_struct * Aterms = A->terms;

    Ai = 0;
    ai = 0;
    if (Ai < Alen)
        ai = Aterms[Ai].coeff->length - 1;

    for (i = 0; i < Z->length; i++)
    {
        if (Ai < Alen && Zterms[i].exp == pack_exp2(Aterms[Ai].exp, ai))
        {
            /* Z present, A present */
            fmpz_set(Zterms[i].coeff->coeffs + cur_length,
                     Aterms[Ai].coeff->coeffs + ai);
            Zterms[i].coeff->length = cur_length + 1;
            do {
                ai--;
            } while (ai >= 0 && fmpz_is_zero(Aterms[Ai].coeff->coeffs + ai));
            if (ai < 0)
            {
                Ai++;
                if (Ai < Alen)
                    ai = Aterms[Ai].coeff->length - 1;
            }
        }
        else if (Ai < 0 || Zterms[i].exp > pack_exp2(Aterms[Ai].exp, ai))
        {
            /* Z present, A missing */
            fmpz_zero(Zterms[i].coeff->coeffs + cur_length);
            Zterms[i].coeff->length = cur_length + 1;
        }
        else
        {
            /* Z missing, A present */
            return 0;
        }
    }

    return 1;
}


int fmpz_mod_mpolyn_interp_mcrt_sm_mpoly(
    slong * lastdeg,
    fmpz_mod_mpolyn_t F,
    const fmpz_mod_mpoly_t A,
    const fmpz_mod_poly_t modulus,
    fmpz_mod_poly_t alphapow,
    const fmpz_mod_mpoly_ctx_t ctx)
{
    slong lastlen = 0;
    int changed = 0;
    slong i;
    fmpz_t v;
    fmpz * Acoeff = A->coeffs;
    slong Flen = F->length;

    FLINT_ASSERT(Flen == A->length);

    fmpz_init(v);

    for (i = 0; i < Flen; i++)
    {
        /* F term ok, A term ok */
        fmpz_mod_poly_eval_pow(v, F->coeffs + i, alphapow, ctx->ffinfo);
        fmpz_mod_sub(v, Acoeff + i, v, ctx->ffinfo);
        if (!fmpz_is_zero(v))
        {
            changed = 1;
            fmpz_mod_poly_scalar_addmul_fmpz_mod(F->coeffs + i,
                                       F->coeffs + i, modulus, v, ctx->ffinfo);
        }
        lastlen = FLINT_MAX(lastlen, F->coeffs[i].length);
    }

    fmpz_clear(v);

    *lastdeg = lastlen - 1;
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
int fmpz_mod_mpolyl_gcd_zippel_smprime(
    fmpz_mod_mpoly_t rG, const slong * rGdegs, /* guess at rG degrees, could be NULL */
    fmpz_mod_mpoly_t rAbar,
    fmpz_mod_mpoly_t rBbar,
    const fmpz_mod_mpoly_t A, const slong * Adegs,
    const fmpz_mod_mpoly_t B, const slong * Bdegs,
    const fmpz_mod_mpoly_t gamma, const slong * gammadegs,
    const fmpz_mod_mpoly_ctx_t ctx)
{
    int success, use, alpha_tries_left;
    slong i, m;
    slong nvars = ctx->minfo->nvars;
    flint_bitcnt_t bits = A->bits;
    fmpz * alphas, * betas;
    flint_rand_t state;
    fmpz_mod_mpoly_t cont;
    fmpz_mod_mpoly_t T, G, Abar, Bbar;
    fmpz_mod_polyun_t HG, HAbar, HBbar, MG, MAbar, MBbar, ZG, ZAbar, ZBbar;
    fmpz_mod_polyun_t Aev, Bev, Gev, Abarev, Bbarev;
    fmpz_t gammaev;
    fmpz_mod_mpolyn_t Tn, Gn, Abarn, Bbarn;
    slong lastdeg;
    slong cur_zip_image, req_zip_images, this_length;
    fmpz_mod_polyun_t Aeh_cur, Aeh_inc, Beh_cur, Beh_inc;
    fmpz_mod_poly_t gammaeh_cur, gammaeh_inc;
    fmpz_mod_poly_t modulus, alphapow;
    fmpz_mod_mpoly_struct * Aevals, * Bevals;
    fmpz_mod_mpoly_struct * gammaevals;
    fmpz_mod_poly_polyun_stack_t St;
    fmpz_t c, start_alpha;
    ulong GdegboundXY, newdegXY, Abideg, Bbideg;
    slong degxAB, degyAB;

    FLINT_ASSERT(bits <= FLINT_BITS);
    FLINT_ASSERT(bits == A->bits);
    FLINT_ASSERT(bits == B->bits);
    FLINT_ASSERT(bits == gamma->bits);
    FLINT_ASSERT(A->length > 0);
    FLINT_ASSERT(B->length > 0);
    FLINT_ASSERT(gamma->length > 0);

    fmpz_mod_mpoly_fit_length_reset_bits(rG, 1, bits, ctx);
    fmpz_mod_mpoly_fit_length_reset_bits(rAbar, 1, bits, ctx);
    fmpz_mod_mpoly_fit_length_reset_bits(rBbar, 1, bits, ctx);

    fmpz_init(gammaev);
    fmpz_init(c);
    fmpz_init(start_alpha);

#if FLINT_WANT_ASSERT
    {
        slong * tmp_degs = FLINT_ARRAY_ALLOC(nvars, slong);

        fmpz_mod_mpoly_degrees_si(tmp_degs, A, ctx);
        for (i = 0; i < nvars; i++)
            FLINT_ASSERT(tmp_degs[i] == Adegs[i]);

        fmpz_mod_mpoly_degrees_si(tmp_degs, B, ctx);
        for (i = 0; i < nvars; i++)
            FLINT_ASSERT(tmp_degs[i] == Bdegs[i]);

        fmpz_mod_mpoly_degrees_si(tmp_degs, gamma, ctx);
        for (i = 0; i < nvars; i++)
            FLINT_ASSERT(tmp_degs[i] == gammadegs[i]);

        flint_free(tmp_degs);
    }
#endif

    FLINT_ASSERT(gammadegs[0] == 0);
    FLINT_ASSERT(gammadegs[1] == 0);

    fmpz_mod_polyun_init(HG, ctx->ffinfo);
    fmpz_mod_polyun_init(HAbar, ctx->ffinfo);
    fmpz_mod_polyun_init(HBbar, ctx->ffinfo);
    fmpz_mod_polyun_init(MG, ctx->ffinfo);
    fmpz_mod_polyun_init(MAbar, ctx->ffinfo);
    fmpz_mod_polyun_init(MBbar, ctx->ffinfo);
    fmpz_mod_polyun_init(ZG, ctx->ffinfo);
    fmpz_mod_polyun_init(ZAbar, ctx->ffinfo);
    fmpz_mod_polyun_init(ZBbar, ctx->ffinfo);
    fmpz_mod_polyun_init(Aev, ctx->ffinfo);
    fmpz_mod_polyun_init(Bev, ctx->ffinfo);
    fmpz_mod_polyun_init(Gev, ctx->ffinfo);
    fmpz_mod_polyun_init(Abarev, ctx->ffinfo);
    fmpz_mod_polyun_init(Bbarev, ctx->ffinfo);
    fmpz_mod_poly_init2(alphapow, 4, ctx->ffinfo);
    fmpz_mod_mpoly_init3(cont, 1, bits, ctx);
    fmpz_mod_mpoly_init3(T, 1, bits, ctx);
    fmpz_mod_mpoly_init3(G, 1, bits, ctx);
    fmpz_mod_mpoly_init3(Abar, 1, bits, ctx);
    fmpz_mod_mpoly_init3(Bbar, 1, bits, ctx);
    fmpz_mod_mpolyn_init(Tn, bits, ctx);
    fmpz_mod_mpolyn_init(Gn, bits, ctx);
    fmpz_mod_mpolyn_init(Abarn, bits, ctx);
    fmpz_mod_mpolyn_init(Bbarn, bits, ctx);
    fmpz_mod_polyun_init(Aeh_cur, ctx->ffinfo);
    fmpz_mod_polyun_init(Aeh_inc, ctx->ffinfo);
    fmpz_mod_polyun_init(Beh_cur, ctx->ffinfo);
    fmpz_mod_polyun_init(Beh_inc, ctx->ffinfo);
    fmpz_mod_poly_init(gammaeh_cur, ctx->ffinfo);
    fmpz_mod_poly_init(gammaeh_inc, ctx->ffinfo);
    fmpz_mod_poly_init(modulus, ctx->ffinfo);
    fmpz_mod_poly_stack_init(St->poly_stack);
    fmpz_mod_polyun_stack_init(St->polyun_stack);

    betas = _fmpz_vec_init(nvars);
    alphas = _fmpz_vec_init(nvars);
    flint_randinit(state);

    Aevals = FLINT_ARRAY_ALLOC(nvars + 1, fmpz_mod_mpoly_struct);
    Bevals = FLINT_ARRAY_ALLOC(nvars + 1, fmpz_mod_mpoly_struct);
    gammaevals = FLINT_ARRAY_ALLOC(nvars + 1, fmpz_mod_mpoly_struct);
    for (i = 0; i < nvars; i++)
    {
        fmpz_mod_mpoly_init3(Aevals + i, 0, bits, ctx);
        fmpz_mod_mpoly_init3(Bevals + i, 0, bits, ctx);
        fmpz_mod_mpoly_init3(gammaevals + i, 0, bits, ctx);
    }
    Aevals[nvars] = *A;
    Bevals[nvars] = *B;
    gammaevals[nvars] = *gamma;

    Abideg = _fmpz_mod_mpoly_bidegree(A, ctx);
    Bbideg = _fmpz_mod_mpoly_bidegree(B, ctx);

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
        fmpz_mod_rand_not_zero(alphas + i, state, ctx->ffinfo);

    for (i = nvars - 1; i >= 2; i--)
    {
        fmpz_mod_mpoly_evaluate_one_fmpz(Aevals + i, Aevals + i + 1, i, alphas + i, ctx);
        fmpz_mod_mpoly_evaluate_one_fmpz(Bevals + i, Bevals + i + 1, i, alphas + i, ctx);
        fmpz_mod_mpoly_evaluate_one_fmpz(gammaevals + i, gammaevals + i + 1, i, alphas + i, ctx);
        if (fmpz_mod_mpoly_is_zero(gammaevals + i, ctx))
            goto choose_alphas;
        if (Aevals[i].length < 1 || _fmpz_mod_mpoly_bidegree(Aevals + i, ctx) != Abideg)
            goto choose_alphas;
        if (Bevals[i].length < 1 || _fmpz_mod_mpoly_bidegree(Bevals + i, ctx) != Bbideg)
            goto choose_alphas;
    }

    m = 2;

    if (rGdegs == NULL)
        use = USE_G | USE_ABAR | USE_BBAR;
    else
        use = mpoly_gcd_get_use_first(rGdegs[m], Adegs[m], Bdegs[m], gammadegs[m]);

    fmpz_mod_mpoly_get_polyu1n(Aev, Aevals + m, 0, 1, ctx);
    fmpz_mod_mpoly_get_polyu1n(Bev, Bevals + m, 0, 1, ctx);

    success = fmpz_mod_polyu1n_gcd_brown_smprime(Gev, Abarev, Bbarev, Aev, Bev, ctx->ffinfo, St);
    if (!success)
        goto cleanup;

    newdegXY = fmpz_mod_polyu1n_bidegree(Gev);
    if (newdegXY > GdegboundXY)
        goto choose_alphas;
    if (newdegXY < GdegboundXY)
    {
        GdegboundXY = newdegXY;
        if (GdegboundXY == 0)
            goto gcd_is_trivial;
    }

    fmpz_mod_mpoly_get_fmpz(gammaev, gammaevals + m, ctx);
    fmpz_mod_polyun_scalar_mul_fmpz(Gev, gammaev, ctx->ffinfo);

    fmpz_mod_mpolyn_interp_lift_sm_polyu1n(Gn, Gev, ctx);
    fmpz_mod_mpolyn_interp_lift_sm_polyu1n(Abarn, Abarev, ctx);
    fmpz_mod_mpolyn_interp_lift_sm_polyu1n(Bbarn, Bbarev, ctx);

    fmpz_mod_poly_one(modulus, ctx->ffinfo);
    fmpz_mod_neg(c, alphas + m, ctx->ffinfo);
    fmpz_mod_poly_shift_left_scalar_addmul_fmpz_mod(modulus, 1, c, ctx->ffinfo);

    fmpz_set(start_alpha, alphas + m);
    while (1)
    {
    choose_alpha_2:

        if (fmpz_cmp_ui(alphas + m, 2) < 0)
            fmpz_set(alphas + m, fmpz_mod_ctx_modulus(ctx->ffinfo));
        fmpz_sub_ui(alphas + m, alphas + m, 1);
        if (fmpz_equal(alphas + m, start_alpha))
            goto choose_alphas;

        FLINT_ASSERT(alphapow->alloc >= 2);
        fmpz_one(alphapow->coeffs + 0);
        fmpz_set(alphapow->coeffs + 1, alphas + m);
        alphapow->length = 2;

        fmpz_mod_mpoly_evaluate_one_fmpz(Aevals + m, Aevals + m + 1, m, alphas + m, ctx);
        fmpz_mod_mpoly_evaluate_one_fmpz(Bevals + m, Bevals + m + 1, m, alphas + m, ctx);

        fmpz_mod_mpoly_evaluate_one_fmpz(gammaevals + m, gammaevals + m + 1, m, alphas + m, ctx);
        if (fmpz_mod_mpoly_is_zero(gammaevals + m, ctx))
            goto choose_alpha_2;
        if (Aevals[m].length < 1 || _fmpz_mod_mpoly_bidegree(Aevals + m, ctx) != Abideg)
            goto choose_alpha_2;
        if (Bevals[m].length < 1 || _fmpz_mod_mpoly_bidegree(Bevals + m, ctx) != Bbideg)
            goto choose_alpha_2;

        fmpz_mod_mpoly_get_polyu1n(Aev, Aevals + m, 0, 1, ctx);
        fmpz_mod_mpoly_get_polyu1n(Bev, Bevals + m, 0, 1, ctx);

        success = fmpz_mod_polyu1n_gcd_brown_smprime(Gev, Abarev, Bbarev, Aev, Bev, ctx->ffinfo, St);
        if (!success)
            goto cleanup;

        newdegXY = fmpz_mod_polyu1n_bidegree(Gev);
        if (newdegXY > GdegboundXY)
            goto choose_alpha_2;
        if (newdegXY < GdegboundXY)
        {
            GdegboundXY = newdegXY;
            if (GdegboundXY == 0)
                goto gcd_is_trivial;
            goto choose_alphas;
        }

        fmpz_mod_mpoly_get_fmpz(gammaev, gammaevals + m, ctx);
        fmpz_mod_polyun_scalar_mul_fmpz(Gev, gammaev, ctx->ffinfo);

        fmpz_mod_poly_eval_pow(c, modulus, alphapow, ctx->ffinfo);
        fmpz_mod_inv(c, c, ctx->ffinfo);
        fmpz_mod_poly_scalar_mul_fmpz(modulus, modulus, c, ctx->ffinfo);

        if ((use & USE_G) && !fmpz_mod_mpolyn_interp_crt_sm_polyu1n(
                                &lastdeg, Gn, Tn, Gev, modulus, alphapow, ctx))
        {
            if (m == nvars - 1)
            {
                fmpz_mod_mpoly_cvtfrom_mpolyn(rG, Gn, m, ctx);
                success = fmpz_mod_mpolyl_content(cont, rG, 2, ctx);
                if (!success)
                    goto cleanup;
                fmpz_mod_mpoly_divides(rG, rG, cont, ctx);
                fmpz_mod_mpoly_repack_bits_inplace(rG, bits, ctx);
                if (fmpz_mod_mpoly_divides(rAbar, A, rG, ctx) &&
                    fmpz_mod_mpoly_divides(rBbar, B, rG, ctx))
                {
                    fmpz_mod_mpoly_repack_bits_inplace(rAbar, bits, ctx);
                    fmpz_mod_mpoly_repack_bits_inplace(rBbar, bits, ctx);
                    break;
                }
            }
            else
            {
                fmpz_mod_mpoly_cvtfrom_mpolyn(G, Gn, m, ctx);
                fmpz_mod_mpoly_mul(T, Aevals + m + 1, gammaevals + m + 1, ctx);
                if (fmpz_mod_mpoly_divides(Abar, T, G, ctx))
                {
                    fmpz_mod_mpoly_repack_bits_inplace(Abar, bits, ctx);
                    fmpz_mod_mpoly_mul(T, Bevals + m + 1, gammaevals + m + 1, ctx);
                    if (fmpz_mod_mpoly_divides(Bbar, T, G, ctx))
                    {
                        fmpz_mod_mpoly_repack_bits_inplace(Bbar, bits, ctx);
                        break;
                    }
                }
            }
        }

        if ((use & USE_ABAR) && !fmpz_mod_mpolyn_interp_crt_sm_polyu1n(
                          &lastdeg, Abarn, Tn, Abarev, modulus, alphapow, ctx))
        {
            if (m == nvars - 1)
            {
                fmpz_mod_mpoly_cvtfrom_mpolyn(rAbar, Abarn, m, ctx);
                success = fmpz_mod_mpolyl_content(cont, rAbar, 2, ctx);
                if (!success)
                    goto cleanup;
                fmpz_mod_mpoly_divides(rAbar, rAbar, cont, ctx);
                fmpz_mod_mpoly_repack_bits_inplace(rAbar, bits, ctx);
                if (fmpz_mod_mpoly_divides(rG, A, rAbar, ctx) &&
                    fmpz_mod_mpoly_divides(rBbar, B, rG, ctx))
                {
                    fmpz_mod_mpoly_repack_bits_inplace(rG, bits, ctx);
                    fmpz_mod_mpoly_repack_bits_inplace(rBbar, bits, ctx);
                    break;
                }                
            }
            else
            {
                fmpz_mod_mpoly_cvtfrom_mpolyn(Abar, Abarn, m, ctx);
                fmpz_mod_mpoly_mul(T, Aevals + m + 1, gammaevals + m + 1, ctx);
                if (fmpz_mod_mpoly_divides(G, T, Abar, ctx))
                {
                    fmpz_mod_mpoly_repack_bits_inplace(G, bits, ctx);
                    fmpz_mod_mpoly_mul(T, Bevals + m + 1, gammaevals + m + 1, ctx);
                    if (fmpz_mod_mpoly_divides(Bbar, T, G, ctx))
                    {
                        fmpz_mod_mpoly_repack_bits_inplace(Bbar, bits, ctx);
                        break;
                    }
                }
            }
        }

        if ((use & USE_BBAR) && !fmpz_mod_mpolyn_interp_crt_sm_polyu1n(
                          &lastdeg, Bbarn, Tn, Bbarev, modulus, alphapow, ctx))
        {
            if (m == nvars - 1)
            {
                fmpz_mod_mpoly_cvtfrom_mpolyn(rBbar, Bbarn, m, ctx);
                success = fmpz_mod_mpolyl_content(cont, rBbar, 2, ctx);
                if (!success)
                    goto cleanup;
                fmpz_mod_mpoly_divides(rBbar, rBbar, cont, ctx);
                fmpz_mod_mpoly_repack_bits_inplace(rBbar, bits, ctx);
                if (fmpz_mod_mpoly_divides(rG, B, rBbar, ctx) &&
                    fmpz_mod_mpoly_divides(rAbar, A, rG, ctx))
                {
                    fmpz_mod_mpoly_repack_bits_inplace(rG, bits, ctx);
                    fmpz_mod_mpoly_repack_bits_inplace(rAbar, bits, ctx);
                    break;
                }                
            }
            else
            {
                fmpz_mod_mpoly_cvtfrom_mpolyn(Bbar, Bbarn, m, ctx);
                fmpz_mod_mpoly_mul(T, Bevals + m + 1, gammaevals + m + 1, ctx);
                if (fmpz_mod_mpoly_divides(G, T, Bbar, ctx))
                {
                    fmpz_mod_mpoly_repack_bits_inplace(G, bits, ctx);
                    fmpz_mod_mpoly_mul(T, Aevals + m + 1, gammaevals + m + 1, ctx);
                    if (fmpz_mod_mpoly_divides(Abar, T, G, ctx))
                    {
                        fmpz_mod_mpoly_repack_bits_inplace(Abar, bits, ctx);
                        break;
                    }
                }
            }
        }

        if (fmpz_mod_poly_degree(modulus, ctx->ffinfo) > gammadegs[m] + Adegs[m] &&
            fmpz_mod_poly_degree(modulus, ctx->ffinfo) > gammadegs[m] + Bdegs[m])
        {
            goto choose_alphas;
        }

        fmpz_mod_neg(c, alphas + m, ctx->ffinfo);
        fmpz_mod_poly_shift_left_scalar_addmul_fmpz_mod(modulus, 1, c, ctx->ffinfo);
    }

    for (m = 3; m < nvars; m++)
    {
        /* G, Abar, Bbar are in Fp[gen(0), ..., gen(m - 1)] */
        fmpz_mod_mpolyn_interp_lift_sm_mpoly(Gn, G, ctx);
        fmpz_mod_mpolyn_interp_lift_sm_mpoly(Abarn, Abar, ctx);
        fmpz_mod_mpolyn_interp_lift_sm_mpoly(Bbarn, Bbar, ctx);

        if (rGdegs == NULL)
            use = USE_G | USE_ABAR | USE_BBAR;
        else
            use = 0;

        fmpz_mod_poly_one(modulus, ctx->ffinfo);
        fmpz_mod_neg(c, alphas + m, ctx->ffinfo);
        fmpz_mod_poly_shift_left_scalar_addmul_fmpz_mod(modulus, 1, c, ctx->ffinfo);

        fmpz_set(start_alpha, alphas + m);

    choose_betas:

        /* only beta[2], beta[1], ..., beta[m - 1] will be used */
        for (i = 2; i < ctx->minfo->nvars; i++)
            fmpz_mod_rand_not_zero(betas + i, state, ctx->ffinfo);

        fmpz_mod_mpoly_monomial_evals2(HG, G, betas + 2, m, ctx);
        fmpz_mod_mpoly_monomial_evals2(HAbar, Abar, betas + 2, m, ctx);
        fmpz_mod_mpoly_monomial_evals2(HBbar, Bbar, betas + 2, m, ctx);

        if (use == 0)
        {
            this_length = gammaevals[m + 1].length + Aevals[m + 1].length +
                                                     Bevals[m + 1].length;

            use = fmpz_mod_mpoly_gcd_get_use_new(rGdegs[m], Adegs[m], Bdegs[m],
                  gammadegs[m], degxAB, degyAB, this_length, HG, HAbar, HBbar);
        }

        req_zip_images = 1;
        if (use & USE_G)
        {
            this_length = fmpz_mod_polyun_product_roots(MG, HG, ctx->ffinfo);
            req_zip_images = FLINT_MAX(req_zip_images, this_length);
        }
        if (use & USE_ABAR)
        {
            this_length = fmpz_mod_polyun_product_roots(MAbar, HAbar, ctx->ffinfo);
            req_zip_images = FLINT_MAX(req_zip_images, this_length);
        }
        if (use & USE_BBAR)
        {
            this_length = fmpz_mod_polyun_product_roots(MBbar, HBbar, ctx->ffinfo);
            req_zip_images = FLINT_MAX(req_zip_images, this_length);
        }

        while (1)
        {
        choose_alpha_m:

            if (fmpz_cmp_ui(alphas + m, 2) < 0)
                fmpz_set(alphas + m, fmpz_mod_ctx_modulus(ctx->ffinfo));
            fmpz_sub_ui(alphas + m, alphas + m, 1);
            if (fmpz_equal(alphas + m, start_alpha))
                goto choose_alphas;

            FLINT_ASSERT(alphapow->alloc >= 2);
            fmpz_one(alphapow->coeffs + 0);
            fmpz_set(alphapow->coeffs + 1, alphas + m);
            alphapow->length = 2;

            fmpz_mod_mpoly_evaluate_one_fmpz(Aevals + m, Aevals + m + 1, m, alphas + m, ctx);
            fmpz_mod_mpoly_evaluate_one_fmpz(Bevals + m, Bevals + m + 1, m, alphas + m, ctx);
            fmpz_mod_mpoly_evaluate_one_fmpz(gammaevals + m, gammaevals + m + 1, m, alphas + m, ctx);
            if (fmpz_mod_mpoly_is_zero(gammaevals + m, ctx))
                goto choose_alpha_m;
            if (Aevals[m].length < 1 || _fmpz_mod_mpoly_bidegree(Aevals + m, ctx) != Abideg)
                goto choose_alpha_m;
            if (Bevals[m].length < 1 || _fmpz_mod_mpoly_bidegree(Bevals + m, ctx) != Bbideg)
                goto choose_alpha_m;

            fmpz_mod_mpoly_monomial_evals2(Aeh_inc, Aevals + m, betas + 2, m, ctx);
            fmpz_mod_mpoly_monomial_evals2(Beh_inc, Bevals + m, betas + 2, m, ctx);
            fmpz_mod_mpoly_monomial_evals(gammaeh_inc, gammaevals + m, betas + 2, 2, m, ctx);
            fmpz_mod_polyun_set(Aeh_cur, Aeh_inc, ctx->ffinfo);
            fmpz_mod_polyun_set(Beh_cur, Beh_inc, ctx->ffinfo);
            fmpz_mod_poly_set(gammaeh_cur, gammaeh_inc, ctx->ffinfo);

            fmpz_mod_polyun_zip_start(ZG, HG, req_zip_images, ctx->ffinfo);
            fmpz_mod_polyun_zip_start(ZAbar, HAbar, req_zip_images, ctx->ffinfo);
            fmpz_mod_polyun_zip_start(ZBbar, HBbar, req_zip_images, ctx->ffinfo);
            for (cur_zip_image = 0; cur_zip_image < req_zip_images; cur_zip_image++)
            {
                fmpz_mod_polyu2n_eval_step_sep(Aev, Aeh_cur, Aeh_inc, Aevals + m, ctx);
                fmpz_mod_polyu2n_eval_step_sep(Bev, Beh_cur, Beh_inc, Bevals + m, ctx);
                fmpz_mod_poly_eval_step_sep(gammaev, gammaeh_cur, gammaeh_inc, gammaevals + m, ctx);
                if (fmpz_is_zero(gammaev))
                    goto choose_betas;
                if (Aev->length < 1 || fmpz_mod_polyu1n_bidegree(Aev) != Abideg)
                    goto choose_betas;
                if (Bev->length < 1 || fmpz_mod_polyu1n_bidegree(Bev) != Bbideg)
                    goto choose_betas;

                success = fmpz_mod_polyu1n_gcd_brown_smprime(Gev, Abarev, Bbarev,
                                                            Aev, Bev, ctx->ffinfo, St);        
                if (!success)
                    goto cleanup;

                newdegXY = fmpz_mod_polyu1n_bidegree(Gev);
                if (newdegXY > GdegboundXY)
                    goto choose_betas;
                if (newdegXY < GdegboundXY)
                {
                    GdegboundXY = newdegXY;
                    if (GdegboundXY == 0)
                        goto gcd_is_trivial;
                    goto choose_alphas;
                }

                fmpz_mod_polyun_scalar_mul_fmpz(Gev, gammaev, ctx->ffinfo);
                if ((use & USE_G) && !fmpz_mod_polyun_add_zip_must_match(ZG, Gev, cur_zip_image))
                    goto choose_alphas;
                if ((use & USE_ABAR) && !fmpz_mod_polyun_add_zip_must_match(ZAbar, Abarev, cur_zip_image))
                    goto choose_alphas;
                if ((use & USE_BBAR) && !fmpz_mod_polyun_add_zip_must_match(ZBbar, Bbarev, cur_zip_image))
                    goto choose_alphas;

            }

            if ((use & USE_G) && fmpz_mod_polyun_zip_solve(G, ZG, HG, MG, ctx) < 1)
                goto choose_alphas;
            if ((use & USE_ABAR) && fmpz_mod_polyun_zip_solve(Abar, ZAbar, HAbar, MAbar, ctx) < 1)
                goto choose_alphas;
            if ((use & USE_BBAR) && fmpz_mod_polyun_zip_solve(Bbar, ZBbar, HBbar, MBbar, ctx) < 1)
                goto choose_alphas;

            FLINT_ASSERT(fmpz_mod_poly_degree(modulus, ctx->ffinfo) > 0);
            fmpz_mod_poly_eval_pow(c, modulus, alphapow, ctx->ffinfo);
            fmpz_mod_inv(c, c, ctx->ffinfo);
            fmpz_mod_poly_scalar_mul_fmpz(modulus, modulus, c, ctx->ffinfo);

            if ((use & USE_G) && !fmpz_mod_mpolyn_interp_mcrt_sm_mpoly(
                                      &lastdeg, Gn, G, modulus, alphapow, ctx))
            {
                fmpz_mod_mpoly_cvtfrom_mpolyn(rG, Gn, m, ctx);
                if (m == nvars - 1)
                {
                    success = fmpz_mod_mpolyl_content(cont, rG, 2, ctx);
                    if (!success)
                        goto cleanup;
                    fmpz_mod_mpoly_divides(rG, rG, cont, ctx);
                    fmpz_mod_mpoly_repack_bits_inplace(rG, bits, ctx);
                    if (fmpz_mod_mpoly_divides(rAbar, A, rG, ctx) &&
                        fmpz_mod_mpoly_divides(rBbar, B, rG, ctx))
                    {
                        fmpz_mod_mpoly_repack_bits_inplace(rAbar, bits, ctx);
                        fmpz_mod_mpoly_repack_bits_inplace(rBbar, bits, ctx);
                        break;
                    }
                }
                else
                {
                    fmpz_mod_mpoly_mul(T, Aevals + m + 1, gammaevals + m + 1, ctx);
                    if (fmpz_mod_mpoly_divides(rAbar, T, rG, ctx))
                    {
                        fmpz_mod_mpoly_mul(T, Bevals + m + 1, gammaevals + m + 1, ctx);
                        if (fmpz_mod_mpoly_divides(rBbar, T, rG, ctx))
                        {
                            fmpz_mod_mpoly_repack_bits_inplace(rAbar, bits, ctx);
                            fmpz_mod_mpoly_repack_bits_inplace(rBbar, bits, ctx);
                            fmpz_mod_mpoly_swap(G, rG, ctx);
                            fmpz_mod_mpoly_swap(Abar, rAbar, ctx);
                            fmpz_mod_mpoly_swap(Bbar, rBbar, ctx);
                            break;
                        }
                    }
                }
            }
            if ((use & USE_ABAR) && !fmpz_mod_mpolyn_interp_mcrt_sm_mpoly(
                                &lastdeg, Abarn, Abar, modulus, alphapow, ctx))
            {
                fmpz_mod_mpoly_cvtfrom_mpolyn(rAbar, Abarn, m, ctx);
                if (m == nvars - 1)
                {
                    success = fmpz_mod_mpolyl_content(cont, rAbar, 2, ctx);
                    if (!success)
                        goto cleanup;
                    fmpz_mod_mpoly_divides(rAbar, rAbar, cont, ctx);
                    fmpz_mod_mpoly_repack_bits_inplace(rAbar, bits, ctx);
                    if (fmpz_mod_mpoly_divides(rG, A, rAbar, ctx) &&
                        fmpz_mod_mpoly_divides(rBbar, B, rG, ctx))
                    {
                        fmpz_mod_mpoly_repack_bits_inplace(rG, bits, ctx);
                        fmpz_mod_mpoly_repack_bits_inplace(rBbar, bits, ctx);
                        break;
                    }
                }
                else
                {
                    fmpz_mod_mpoly_mul(T, Aevals + m + 1, gammaevals + m + 1, ctx);
                    if (fmpz_mod_mpoly_divides(rG, T, rAbar, ctx))
                    {
                        fmpz_mod_mpoly_mul(T, Bevals + m + 1, gammaevals + m + 1, ctx);
                        if (fmpz_mod_mpoly_divides(rBbar, T, rG, ctx))
                        {
                            fmpz_mod_mpoly_repack_bits_inplace(rG, bits, ctx);
                            fmpz_mod_mpoly_repack_bits_inplace(rBbar, bits, ctx);
                            fmpz_mod_mpoly_swap(G, rG, ctx);
                            fmpz_mod_mpoly_swap(Abar, rAbar, ctx);
                            fmpz_mod_mpoly_swap(Bbar, rBbar, ctx);
                            break;
                        }
                    }
                }
            }
            if ((use & USE_BBAR) && !fmpz_mod_mpolyn_interp_mcrt_sm_mpoly(
                                &lastdeg, Bbarn, Bbar, modulus, alphapow, ctx))
            {
                fmpz_mod_mpoly_cvtfrom_mpolyn(rBbar, Bbarn, m, ctx);
                if (m == nvars - 1)
                {
                    success = fmpz_mod_mpolyl_content(cont, rBbar, 2, ctx);
                    if (!success)
                        goto cleanup;
                    fmpz_mod_mpoly_divides(rBbar, rBbar, cont, ctx);
                    fmpz_mod_mpoly_repack_bits_inplace(rBbar, bits, ctx);
                    if (fmpz_mod_mpoly_divides(rG, B, rBbar, ctx) &&
                        fmpz_mod_mpoly_divides(rAbar, A, rG, ctx))
                    {
                        fmpz_mod_mpoly_repack_bits_inplace(rG, bits, ctx);
                        fmpz_mod_mpoly_repack_bits_inplace(rAbar, bits, ctx);
                        break;
                    }
                }
                else
                {
                    fmpz_mod_mpoly_mul(T, Bevals + m + 1, gammaevals + m + 1, ctx);
                    if (fmpz_mod_mpoly_divides(rG, T, rBbar, ctx))
                    {
                        fmpz_mod_mpoly_mul(T, Aevals + m + 1, gammaevals + m + 1, ctx);
                        if (fmpz_mod_mpoly_divides(rAbar, T, rG, ctx))
                        {
                            fmpz_mod_mpoly_repack_bits_inplace(rG, bits, ctx);
                            fmpz_mod_mpoly_repack_bits_inplace(rAbar, bits, ctx);
                            fmpz_mod_mpoly_swap(G, rG, ctx);
                            fmpz_mod_mpoly_swap(Abar, rAbar, ctx);
                            fmpz_mod_mpoly_swap(Bbar, rBbar, ctx);
                            break;
                        }
                    }
                }
            }

            if (fmpz_mod_poly_degree(modulus, ctx->ffinfo) > gammadegs[m] + Adegs[m] &&
                fmpz_mod_poly_degree(modulus, ctx->ffinfo) > gammadegs[m] + Bdegs[m])
            {
                goto choose_alphas;
            }

            fmpz_mod_neg(c, alphas + m, ctx->ffinfo);
            fmpz_mod_poly_shift_left_scalar_addmul_fmpz_mod(modulus, 1, c, ctx->ffinfo);
        }
    }

    success = 1;

cleanup:

    fmpz_mod_polyun_clear(HG, ctx->ffinfo);
    fmpz_mod_polyun_clear(HAbar, ctx->ffinfo);
    fmpz_mod_polyun_clear(HBbar, ctx->ffinfo);
    fmpz_mod_polyun_clear(MG, ctx->ffinfo);
    fmpz_mod_polyun_clear(MAbar, ctx->ffinfo);
    fmpz_mod_polyun_clear(MBbar, ctx->ffinfo);
    fmpz_mod_polyun_clear(ZG, ctx->ffinfo);
    fmpz_mod_polyun_clear(ZAbar, ctx->ffinfo);
    fmpz_mod_polyun_clear(ZBbar, ctx->ffinfo);
    fmpz_mod_polyun_clear(Aev, ctx->ffinfo);
    fmpz_mod_polyun_clear(Bev, ctx->ffinfo);
    fmpz_mod_polyun_clear(Gev, ctx->ffinfo);
    fmpz_mod_polyun_clear(Abarev, ctx->ffinfo);
    fmpz_mod_polyun_clear(Bbarev, ctx->ffinfo);
    fmpz_mod_poly_clear(alphapow, ctx->ffinfo);
    fmpz_mod_mpoly_clear(cont, ctx);
    fmpz_mod_mpoly_clear(T, ctx);
    fmpz_mod_mpoly_clear(G, ctx);
    fmpz_mod_mpoly_clear(Abar, ctx);
    fmpz_mod_mpoly_clear(Bbar, ctx);
    fmpz_mod_mpolyn_clear(Tn, ctx);
    fmpz_mod_mpolyn_clear(Gn, ctx);
    fmpz_mod_mpolyn_clear(Abarn, ctx);
    fmpz_mod_mpolyn_clear(Bbarn, ctx);
    fmpz_mod_polyun_clear(Aeh_cur, ctx->ffinfo);
    fmpz_mod_polyun_clear(Aeh_inc, ctx->ffinfo);
    fmpz_mod_polyun_clear(Beh_cur, ctx->ffinfo);
    fmpz_mod_polyun_clear(Beh_inc, ctx->ffinfo);
    fmpz_mod_poly_clear(gammaeh_cur, ctx->ffinfo);
    fmpz_mod_poly_clear(gammaeh_inc, ctx->ffinfo);
    fmpz_mod_poly_clear(modulus, ctx->ffinfo);
    fmpz_mod_poly_stack_clear(St->poly_stack);
    fmpz_mod_polyun_stack_clear(St->polyun_stack);

    flint_free(betas);
    flint_free(alphas);
    flint_randclear(state);

    for (i = 0; i < nvars; i++)
    {
        fmpz_mod_mpoly_clear(Aevals + i, ctx);
        fmpz_mod_mpoly_clear(Bevals + i, ctx);
        fmpz_mod_mpoly_clear(gammaevals + i, ctx);
    }
    flint_free(Aevals);
    flint_free(Bevals);
    flint_free(gammaevals);

    fmpz_clear(gammaev);
    fmpz_clear(c);
    fmpz_clear(start_alpha);

    FLINT_ASSERT(!success || rG->bits == bits);
    FLINT_ASSERT(!success || rAbar->bits == bits);
    FLINT_ASSERT(!success || rBbar->bits == bits);

    return success;

gcd_is_trivial:

    fmpz_mod_mpoly_one(rG, ctx);
    fmpz_mod_mpoly_set(rAbar, A, ctx);
    fmpz_mod_mpoly_set(rBbar, B, ctx);

    success = 1;
    
    goto cleanup;
}


/* should find its way back here in interesting cases */
int fmpz_mod_mpoly_gcd_zippel2(
    fmpz_mod_mpoly_t G,
    const fmpz_mod_mpoly_t A,
    const fmpz_mod_mpoly_t B,
    const fmpz_mod_mpoly_ctx_t ctx)
{
    if (fmpz_mod_mpoly_is_zero(A, ctx) || fmpz_mod_mpoly_is_zero(B, ctx))
        return fmpz_mod_mpoly_gcd(G, A, B, ctx);

    return _fmpz_mod_mpoly_gcd_algo(G, NULL, NULL, A, B, ctx, MPOLY_GCD_USE_ZIPPEL2);
}

