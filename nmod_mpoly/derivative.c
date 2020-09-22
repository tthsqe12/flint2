/*
    Copyright (C) 2017 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "nmod_mpoly.h"


slong _nmod_mpoly_derivative(
    mp_limb_t * coeff1, ulong * exp1,
    const mp_limb_t * coeff2, const ulong * exp2, slong len2,
    flint_bitcnt_t bits,
    slong N,
    slong offset,
    slong shift,
    ulong * oneexp,
    const nmodf_ctx_t fctx)
{
    slong i, len1;

    /* x^c -> c*x^(c-1) */
    len1 = 0;
    for (i = 0; i < len2; i++)
    {
        mp_limb_t cr;
        ulong mask = (-UWORD(1)) >> (FLINT_BITS - bits);
        ulong c = (exp2[N*i + offset] >> shift) & mask;
        if (c == 0)
            continue;
        NMOD_RED(cr, c, fctx->mod);
        if (cr == 0)
            continue;
        coeff1[len1] = nmod_mul(coeff2[i], cr, fctx->mod);
        mpoly_monomial_sub(exp1 + N*len1, exp2 + N*i, oneexp, N);
        len1++;
    }

    return len1;
}


slong _nmod_mpoly_derivative_mp(
    mp_limb_t * coeff1, ulong * exp1,
    const mp_limb_t * coeff2, const ulong * exp2, slong len2,
    flint_bitcnt_t bits,
    slong N,
    slong offset,
    ulong * oneexp,
    const nmodf_ctx_t fctx)
{
    slong i, len1;
    fmpz_t c;
    fmpz_init(c);

    /* x^c -> c*x^(c-1) */
    len1 = 0;
    for (i = 0; i < len2; i++)
    {
        mp_limb_t cr;
        fmpz_set_ui_array(c, exp2 + N*i + offset, bits/FLINT_BITS);
        if (fmpz_is_zero(c))
            continue;
        cr = fmpz_fdiv_ui(c, fctx->mod.n);
        if (cr == 0)
            continue;
        coeff1[len1] = nmod_mul(coeff2[i], cr, fctx->mod);
        mpoly_monomial_sub_mp(exp1 + N*len1, exp2 + N*i, oneexp, N);
        len1++;
    }

    fmpz_clear(c);

    return len1;
}


void nmod_mpoly_derivative(
    nmod_mpoly_t poly1,
    const nmod_mpoly_t poly2,
    slong var,
    const nmod_mpoly_ctx_t ctx)
{
    flint_bitcnt_t bits = poly2->bits;
    slong N = mpoly_words_per_exp(bits, ctx->minfo);
    slong offset, shift;
    ulong * oneexp;
    slong len1;
    TMP_INIT;

    TMP_START;

    nmod_mpoly_fit_length_reset_bits(poly1, poly2->length, bits, ctx);

    oneexp = (ulong *) TMP_ALLOC(N*sizeof(ulong));

    if (bits <= FLINT_BITS)
    {
        mpoly_gen_monomial_offset_shift_sp(oneexp, &offset, &shift,
                                                        var, bits, ctx->minfo);

        len1 = _nmod_mpoly_derivative(poly1->coeffs, poly1->exps,
                                  poly2->coeffs, poly2->exps, poly2->length,
                                  bits, N, offset, shift, oneexp, ctx->ffinfo);
    }
    else
    {
        offset = mpoly_gen_monomial_offset_mp(oneexp, var, bits, ctx->minfo);

        len1 = _nmod_mpoly_derivative_mp(poly1->coeffs, poly1->exps,
                                  poly2->coeffs, poly2->exps, poly2->length,
                                  bits, N, offset,        oneexp, ctx->ffinfo);        
    }

    _nmod_mpoly_set_length(poly1, len1, ctx);

    TMP_END;
}
