/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fmpz_mpoly_factor.h"


static void fmpz_mpoly_from_fmpz_bpoly(
    fmpz_mpoly_t A,
    flint_bitcnt_t Abits,
    const fmpz_bpoly_t B,
    slong varx,
    slong vary,
    const fmpz_mpoly_ctx_t ctx)
{
    slong n = ctx->minfo->nvars;
    slong i, j;
    slong NA;
    slong Alen;
    fmpz * Acoeff;
    ulong * Aexp;
    slong Aalloc;
    ulong * Aexps;
    TMP_INIT;

    FLINT_ASSERT(B->length > 0);
    FLINT_ASSERT(Abits <= FLINT_BITS);

    TMP_START;

    Aexps = (ulong *) TMP_ALLOC(n*sizeof(ulong));
    for (i = 0; i < n; i++)
        Aexps[i] = 0;

    NA = mpoly_words_per_exp(Abits, ctx->minfo);
    fmpz_mpoly_fit_bits(A, Abits, ctx);
    A->bits = Abits;

    Acoeff = A->coeffs;
    Aexp = A->exps;
    Aalloc = A->alloc;
    Alen = 0;
    for (i = 0; i < B->length; i++)
    {
        fmpz_poly_struct * Bc = B->coeffs + i;
        _fmpz_mpoly_fit_length(&Acoeff, &Aexp, &Aalloc, Alen + Bc->length, NA);

        for (j = 0; j < Bc->length; j++)
        {
            if (fmpz_is_zero(Bc->coeffs + j))
                continue;
            Aexps[varx] = i;
            Aexps[vary] = j;
            fmpz_set(Acoeff + Alen, Bc->coeffs + j);
            mpoly_set_monomial_ui(Aexp + NA*Alen, Aexps, Abits, ctx->minfo);
            Alen++;
        }
    }
    A->coeffs = Acoeff;
    A->exps = Aexp;
    A->alloc = Aalloc;
    _fmpz_mpoly_set_length(A, Alen, ctx);

    fmpz_mpoly_sort_terms(A, ctx);
    TMP_END;
}

static void fmpz_mpoly_to_bpoly(
    fmpz_bpoly_t A,
    const fmpz_mpoly_t B,
    slong varx,
    slong vary,
    const fmpz_mpoly_ctx_t ctx)
{
    slong j;
    slong NB;
    ulong Bexpx, Bexpy;
    slong Boffx, Bshiftx, Boffy, Bshifty;
    ulong mask;

    FLINT_ASSERT(B->bits <= FLINT_BITS);
    NB = mpoly_words_per_exp_sp(B->bits, ctx->minfo);

    mpoly_gen_offset_shift_sp(&Boffx, &Bshiftx, varx, B->bits, ctx->minfo);
    mpoly_gen_offset_shift_sp(&Boffy, &Bshifty, vary, B->bits, ctx->minfo);
    mask = (-UWORD(1)) >> (FLINT_BITS - B->bits);

    fmpz_bpoly_zero(A);

    for (j = 0; j < B->length; j++)
    {
        Bexpx = ((B->exps + NB*j)[Boffx] >> Bshiftx) & mask;
        Bexpy = ((B->exps + NB*j)[Boffy] >> Bshifty) & mask;
        fmpz_bpoly_set_coeff(A, Bexpx, Bexpy, B->coeffs + j);
    }
}


int fmpz_mpoly_factor_irred_bivar(
    fmpz_mpoly_factor_t fac,
    const fmpz_mpoly_t A,
    slong xvar,
    slong yvar,
    const fmpz_mpoly_ctx_t ctx)
{
    flint_bitcnt_t bits = A->bits;
    slong i;
    fmpz_poly_t c;
    fmpz_poly_factor_t cfac;
    fmpz_bpoly_t B;
    fmpz_tpoly_t F;
    fmpz_mpoly_t t;

    fmpz_poly_init(c);
    fmpz_poly_factor_init(cfac);
    fmpz_bpoly_init(B);
    fmpz_tpoly_init(F);
    fmpz_mpoly_init(t, ctx);

    fmpz_mpoly_to_bpoly(B, A, xvar, yvar, ctx);

    fmpz_bpoly_factor(c, F, B);
    fmpz_poly_factor(cfac, c);

    fac->num = 0;
    fmpz_swap(fac->constant, &cfac->c);
    for (i = 0; i < cfac->num; i++)
    {
        _fmpz_mpoly_set_fmpz_poly(t, bits, cfac->p[i].coeffs,
                                                 cfac->p[i].length, yvar, ctx);
        fmpz_mpoly_factor_append_ui_swap(fac, t, cfac->exp[i], ctx);
    }
    for (i = 0; i < F->length; i++)
    {
        fmpz_mpoly_from_fmpz_bpoly(t, bits, F->coeffs + i, xvar, yvar, ctx);
        fmpz_mpoly_factor_append_ui_swap(fac, t, 1, ctx);
    }
    fmpz_mpoly_factor_fix_units(fac, ctx);

    fmpz_poly_clear(c);
    fmpz_poly_factor_clear(cfac);
    fmpz_bpoly_clear(B);
    fmpz_tpoly_clear(F);
    fmpz_mpoly_clear(t, ctx);
/*
flint_printf("A: "); fmpz_mpoly_print_pretty(A, NULL, ctx); flint_printf("\n");
flint_printf("fac: "); fmpz_mpoly_factor_print_pretty(fac, NULL, ctx); flint_printf("\n");
*/
    FLINT_ASSERT(fmpz_mpoly_factor_matches(A, fac, ctx));

    return 1;
}
