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

/*
    evaluation at
        gen(start) -> caches[0]
        gen(start+1) -> caches[1]
        ...
        gen(stop-1) -> caches[stop-start-1]

    the other gen are assumed to not appear in A
*/
void mpoly_monomial_evals_fmpz_mod(
    fmpz_mod_poly_t EH,
    const ulong * Aexps, slong Alen, flint_bitcnt_t Abits,
    fmpz_mod_poly_struct * alpha_caches,
    slong start,
    slong stop,
    const mpoly_ctx_t mctx,
    const fmpz_mod_ctx_t fpctx)
{
    slong i, k;
    fmpz * p;
    ulong mask = (-UWORD(1)) >> (FLINT_BITS - Abits);
    slong N = mpoly_words_per_exp_sp(Abits, mctx);
    slong * off, * shift;
    slong num = stop - start;
    TMP_INIT;

    TMP_START;

    off = TMP_ARRAY_ALLOC(2*num, slong);
    shift = off + num;
    for (k = 0; k < num; k++)
        mpoly_gen_offset_shift_sp(&off[k], &shift[k], k + start, Abits, mctx);

    fmpz_mod_poly_fit_length(EH, Alen, fpctx);
    EH->length = Alen;
    p = EH->coeffs;

    for (i = 0; i < Alen; i++)
    {
        fmpz_one(p + i);
        for (k = 0; k < num; k++)
        {
            ulong ei = (Aexps[N*i + off[k]] >> shift[k]) & mask;
            fmpz_mod_pow_cache_mulpow_ui(p + i, p + i, ei, alpha_caches + k, fpctx);
        }
    }

    TMP_END;
}


/*
    evaluation at

    gen(0) -> x
    gen(1) -> y
    gen(2) -> alphas[0]
    gen(3) -> alphas[1]
    ...
    gen(m-1) -> alphas[m-3]

    the bivariate marks should be filled in by mpoly_monomials2_fill_marks
*/
void mpoly2_monomial_evals_fmpz_mod(
    fmpz_mod_polyun_t EH,
    const ulong * Aexps, flint_bitcnt_t Abits, ulong * Amarks, slong Amarkslen,
    fmpz_mod_poly_struct * alpha_caches,
    slong m,
    const mpoly_ctx_t mctx,
    const fmpz_mod_ctx_t fpctx)
{
    slong start, stop, i, j, k, n;
    ulong e0, e1;
    fmpz * p;
    ulong mask = (-UWORD(1)) >> (FLINT_BITS - Abits);
    slong N = mpoly_words_per_exp_sp(Abits, mctx);
    slong * off, * shift;
    TMP_INIT;

    FLINT_ASSERT(2 < m && m <= mctx->nvars);
    FLINT_ASSERT(Amarkslen > 0);

    TMP_START;

    off = TMP_ARRAY_ALLOC(2*m, slong);
    shift = off + m;
    for (k = 0; k < m; k++)
        mpoly_gen_offset_shift_sp(&off[k], &shift[k], k, Abits, mctx);

    fmpz_mod_polyun_fit_length(EH, Amarkslen, fpctx);

    for (i = 0; i < Amarkslen; i++)
    {
        start = Amarks[i];
        stop = Amarks[i + 1];
        FLINT_ASSERT(start < stop);
        n = stop - start;

        e0 = (Aexps[N*start + off[0]] >> shift[0]) & mask;
        e1 = (Aexps[N*start + off[1]] >> shift[1]) & mask;

        EH->exps[i] = pack_exp2(e0, e1);
        fmpz_mod_poly_fit_length(EH->coeffs + i, n, fpctx);
        EH->coeffs[i].length = n;
        p = EH->coeffs[i].coeffs;

        for (j = 0; j < n; j++)
        {
            fmpz_one(p + j);
            for (k = 2; k < m; k++)
            {
                ulong ei = (Aexps[N*(start + j) + off[k]] >> shift[k]) & mask;
                fmpz_mod_pow_cache_mulpow_ui(p + j, p + j, ei,
                                                  alpha_caches + k - 2, fpctx);
            }
        }
    }

    EH->length = Amarkslen;

    TMP_END;
}


void fmpz_mod_polyu2n_zip_eval_cur_inc_coeff(
    fmpz_mod_polyun_t E,
    fmpz_mod_polyun_t Acur,
    const fmpz_mod_polyun_t Ainc,
    const fmpz_mod_polyun_t Acoeff,
    const fmpz_mod_ctx_t ctx)
{
    slong i, Ei;
    slong e0, e1;
    fmpz_t c;

    FLINT_ASSERT(Acur->length > 0);
    FLINT_ASSERT(Acur->length == Ainc->length);
    FLINT_ASSERT(Acur->length == Acoeff->length);

    fmpz_init(c);

    e0 = extract_exp(Acur->exps[0], 1, 2);
    e1 = extract_exp(Acur->exps[0], 0, 2);

    fmpz_mod_polyun_fit_length(E, 4, ctx);
    Ei = 0;
    E->exps[Ei] = e1;
    fmpz_mod_poly_zero(E->coeffs + Ei, ctx);

    for (i = 0; i < Acur->length; i++)
    {
        slong this_len = Acur->coeffs[i].length;
        FLINT_ASSERT(this_len == Ainc->coeffs[i].length);
        FLINT_ASSERT(this_len == Acoeff->coeffs[i].length);

        _fmpz_mod_zip_eval_step(c, Acur->coeffs[i].coeffs,
              Ainc->coeffs[i].coeffs, Acoeff->coeffs[i].coeffs, this_len, ctx);

        e0 = extract_exp(Acur->exps[i], 1, 2);
        e1 = extract_exp(Acur->exps[i], 0, 2);

        if (E->exps[Ei] != e0)
        {
            fmpz_mod_polyun_fit_length(E, Ei + 2, ctx);
            Ei += !fmpz_mod_poly_is_zero(E->coeffs + Ei, ctx);
            E->exps[Ei] = e0;
            fmpz_mod_poly_zero(E->coeffs + Ei, ctx);
        }

        fmpz_mod_poly_set_coeff_fmpz(E->coeffs + Ei, e1, c, ctx);
    }

    Ei += !fmpz_mod_poly_is_zero(E->coeffs + Ei, ctx);
    E->length = Ei;

    FLINT_ASSERT(fmpz_mod_polyun_is_canonical(E, ctx));

    fmpz_clear(c);
}

