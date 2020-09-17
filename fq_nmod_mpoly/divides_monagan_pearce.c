/*
    Copyright (C) 2019 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fq_nmod_mpoly.h"

int _fq_nmod_mpoly_divides_monagan_pearce(fq_nmod_mpoly_t Q,
                const mp_limb_t * coeff2, const ulong * exp2, slong len2,
                const mp_limb_t * coeff3, const ulong * exp3, slong len3,
     flint_bitcnt_t bits, slong N, const ulong * cmpmask, const fq_nmod_ctx_t fqctx)
{
    slong d = fq_nmod_ctx_degree(fqctx);
    int lt_divides;
    slong i, j, s;
    slong next_loc, heap_len;
    mpoly_heap_s * heap;
    mpoly_heap_t * chain;
    slong * store, * store_base;
    mpoly_heap_t * x;
    mp_limb_t * Qcoeffs = Q->coeffs;
    ulong * Qexps = Q->exps;
    slong Qlen;
    ulong * exp, * exps;
    ulong ** exp_list;
    slong exp_next;
    mp_limb_t * lc_minus_inv, * pp;
    ulong mask;
    slong * hind;
    TMP_INIT;

    TMP_START;

    pp = (mp_limb_t *) TMP_ALLOC(d*sizeof(mp_limb_t));
    lc_minus_inv = (mp_limb_t *) TMP_ALLOC(d*sizeof(mp_limb_t));

    next_loc = len3 + 4;   /* something bigger than heap can ever be */
    heap = (mpoly_heap_s *) TMP_ALLOC((len3 + 1)*sizeof(mpoly_heap_s));
    chain = (mpoly_heap_t *) TMP_ALLOC(len3*sizeof(mpoly_heap_t));
    store = store_base = (slong *) TMP_ALLOC(2*len3*sizeof(mpoly_heap_t *));
    exps = (ulong *) TMP_ALLOC(len3*N*sizeof(ulong));
    exp_list = (ulong **) TMP_ALLOC(len3*sizeof(ulong *));
    exp_next = 0;
    for (i = 0; i < len3; i++)
        exp_list[i] = exps + i*N;

    hind = (slong *) TMP_ALLOC(len3*sizeof(slong));
    for (i = 0; i < len3; i++)
        hind[i] = 1;

    /* mask with high bit set in each word of each field of exponent vector */
    mask = 0;
    for (i = 0; i < FLINT_BITS/bits; i++)
        mask = (mask << bits) + (UWORD(1) << (bits - 1));

    Qlen = 0;

    /* s is the number of terms * (latest quotient) we should put into heap */
    s = len3;

    /* insert (-1, 0, exp2[0]) into heap */
    heap_len = 2;
    x = chain + 0;
    x->i = -WORD(1);
    x->j = 0;
    x->next = NULL;
    heap[1].next = x;
    heap[1].exp = exp_list[exp_next++];
    mpoly_monomial_set(heap[1].exp, exp2, N);

    /* precompute leading cofficient info */
    n_fq_inv(lc_minus_inv, coeff3 + d*0, fqctx);
    _n_fq_neg(lc_minus_inv, lc_minus_inv, d, fqctx->mod);

    while (heap_len > 1)
    {
        exp = heap[1].exp;

        if (bits <= FLINT_BITS)
        {
            if (mpoly_monomial_overflows(exp, N, mask))
                goto not_exact_division;
        }
        else
        {
            if (mpoly_monomial_overflows_mp(exp, N, bits))
                goto not_exact_division;
        }

        _fq_nmod_mpoly_fit_length(&Qcoeffs, &Q->coeffs_alloc, d,
                                  &Qexps, &Q->exps_alloc, N, Qlen + 1);

        if (bits <= FLINT_BITS)
            lt_divides = mpoly_monomial_divides(Qexps + N*Qlen, exp, exp3, N, mask);
        else
            lt_divides = mpoly_monomial_divides_mp(Qexps + N*Qlen, exp, exp3, N, bits);

        _n_fq_zero(Qcoeffs + d*Qlen, d);
        do {
            exp_list[--exp_next] = heap[1].exp;
            x = _mpoly_heap_pop(heap, &heap_len, N, cmpmask);
            do
            {
                *store++ = x->i;
                *store++ = x->j;

                if (x->i == -WORD(1))
                {
                    n_fq_sub(Qcoeffs + d*Qlen, Qcoeffs + d*Qlen, coeff2 + d*x->j, fqctx);
                }
                else
                {
                    hind[x->i] |= WORD(1);
                    n_fq_mul(pp, coeff3 + d*x->i, Qcoeffs + d*x->j, fqctx);
                    n_fq_add(Qcoeffs + d*Qlen, Qcoeffs + d*Qlen, pp, fqctx);
                }
            } while ((x = x->next) != NULL);
        } while (heap_len > 1 && mpoly_monomial_equal(heap[1].exp, exp, N));

        /* process nodes taken from the heap */
        while (store > store_base)
        {
            j = *--store;
            i = *--store;

            if (i == -WORD(1))
            {
                /* take next dividend term */
                if (j + 1 < len2)
                {
                    x = chain + 0;
                    x->i = i;
                    x->j = j + 1;
                    x->next = NULL;
                    mpoly_monomial_set(exp_list[exp_next], exp2 + x->j*N, N);
                    exp_next += _mpoly_heap_insert(heap, exp_list[exp_next], x,
                                             &next_loc, &heap_len, N, cmpmask);
                }
            }
            else
            {
                /* should we go up */
                if (  (i + 1 < len3)
                   && (hind[i + 1] == 2*j + 1)
                   )
                {
                    x = chain + i + 1;
                    x->i = i + 1;
                    x->j = j;
                    x->next = NULL;
                    hind[x->i] = 2*(x->j + 1) + 0;

                    if (bits <= FLINT_BITS)
                        mpoly_monomial_add(exp_list[exp_next], exp3 + N*x->i,
                                                              Qexps + N*x->j, N);
                    else
                        mpoly_monomial_add_mp(exp_list[exp_next], exp3 + N*x->i,
                                                                 Qexps + N*x->j, N);

                    exp_next += _mpoly_heap_insert(heap, exp_list[exp_next], x,
                                             &next_loc, &heap_len, N, cmpmask);
                }
                /* should we go up? */
                if (j + 1 == Qlen)
                {
                    s++;
                }
                else if (((hind[i] & 1) == 1) &&
                         ((i == 1) || (hind[i - 1] >= 2*(j + 2) + 1)))
                {
                    x = chain + i;
                    x->i = i;
                    x->j = j + 1;
                    x->next = NULL;
                    hind[x->i] = 2*(x->j + 1) + 0;

                    if (bits <= FLINT_BITS)
                        mpoly_monomial_add(exp_list[exp_next], exp3 + N*x->i,
                                                              Qexps + N*x->j, N);
                    else
                        mpoly_monomial_add_mp(exp_list[exp_next], exp3 + N*x->i,
                                                                 Qexps + N*x->j, N);

                    exp_next += _mpoly_heap_insert(heap, exp_list[exp_next], x,
                                             &next_loc, &heap_len, N, cmpmask);
                }
            }
        }

        if (_n_fq_is_zero(Qcoeffs + d*Qlen, d))
        {
            continue;
        }

        n_fq_mul(Qcoeffs + d*Qlen, Qcoeffs + d*Qlen, lc_minus_inv, fqctx);

        if (!lt_divides ||
            mpoly_monomial_gt(exp2 + N*(len2 - 1), exp, N, cmpmask))
        {
            goto not_exact_division;
        }

        if (s > 1)
        {
            i = 1;
            x = chain + i;
            x->i = i;
            x->j = Qlen;
            x->next = NULL;
            hind[x->i] = 2*(x->j + 1) + 0;

            if (bits <= FLINT_BITS)
                mpoly_monomial_add(exp_list[exp_next], exp3 + x->i*N, Qexps + N*x->j, N);
            else
                mpoly_monomial_add_mp(exp_list[exp_next], exp3 + x->i*N, Qexps + N*x->j, N);

            exp_next += _mpoly_heap_insert(heap, exp_list[exp_next], x,
                                             &next_loc, &heap_len, N, cmpmask);
        }
        s = 1;      
        Qlen++;
    }

    Q->coeffs = Qcoeffs;
    Q->exps = Qexps;
    Q->length = Qlen;

    TMP_END;

    return 1;

not_exact_division:

    Q->coeffs = Qcoeffs;
    Q->exps = Qexps;
    Q->length = 0;

    TMP_END;

    return 0;
}

/* return 1 if quotient is exact */
int fq_nmod_mpoly_divides_monagan_pearce(
    fq_nmod_mpoly_t Q,
    const fq_nmod_mpoly_t A,
    const fq_nmod_mpoly_t B,
    const fq_nmod_mpoly_ctx_t ctx)
{
    slong i, N;
    flint_bitcnt_t Qbits;
    fmpz * Amaxfields, * Bmaxfields;
    ulong * cmpmask;
    ulong * Aexps = A->exps, * Bexps = B->exps, * expq;
    int success, easy_exit, freeAexps = 0, freeBexps = 0;
    ulong mask = 0;
    TMP_INIT;

    if (fq_nmod_mpoly_is_zero(B, ctx))
    {
        flint_throw(FLINT_DIVZERO, "Divide by zero in fq_nmod_mpoly_divides_monagan_pearce");
    }

    if (fq_nmod_mpoly_is_zero(A, ctx))
    {
        fq_nmod_mpoly_zero(Q, ctx);
        return 1;
    }

    TMP_START;

    Amaxfields = (fmpz *) TMP_ALLOC(ctx->minfo->nfields*sizeof(fmpz));
    Bmaxfields = (fmpz *) TMP_ALLOC(ctx->minfo->nfields*sizeof(fmpz));
    for (i = 0; i < ctx->minfo->nfields; i++)
    {
        fmpz_init(Amaxfields + i);
        fmpz_init(Bmaxfields + i);
    }

    mpoly_max_fields_fmpz(Amaxfields, Aexps, A->length, A->bits, ctx->minfo);
    mpoly_max_fields_fmpz(Bmaxfields, Bexps, B->length, B->bits, ctx->minfo);
    easy_exit = 0;
    for (i = 0; i < ctx->minfo->nfields; i++)
    {
        if (fmpz_cmp(Amaxfields + i, Bmaxfields + i) < 0)
            easy_exit = 1;
    }

    Qbits = 1 + _fmpz_vec_max_bits(Amaxfields, ctx->minfo->nfields);
    Qbits = FLINT_MAX(Qbits, A->bits);
    Qbits = FLINT_MAX(Qbits, B->bits);
    Qbits = mpoly_fix_bits(Qbits, ctx->minfo);

    for (i = 0; i < ctx->minfo->nfields; i++)
    {
        fmpz_clear(Amaxfields + i);
        fmpz_clear(Bmaxfields + i);
    }

    if (easy_exit)
    {
        success = 0;
        goto cleanup;
    }

    N = mpoly_words_per_exp(Qbits, ctx->minfo);
    cmpmask = (ulong *) TMP_ALLOC(N*sizeof(ulong));
    mpoly_get_cmpmask(cmpmask, N, Qbits, ctx->minfo);

    /* temporary space to check leading monomials divide */
    expq = (ulong *) TMP_ALLOC(N*sizeof(ulong));

    /* quick check for easy case of inexact division of leading monomials */
    if (Qbits == A->bits && Qbits == B->bits && A->exps[N - 1] < B->exps[N - 1])
    {
        success = 0;
        goto cleanup;
    }

    /* ensure input exponents packed to same size as output exponents */
    if (Qbits > A->bits)
    {
        freeAexps = 1;
        Aexps = (ulong *) flint_malloc(N*A->length*sizeof(ulong));
        mpoly_repack_monomials(Aexps, Qbits, A->exps, A->bits, A->length, ctx->minfo);
    }

    if (Qbits > B->bits)
    {
        freeBexps = 1;
        Bexps = (ulong *) flint_malloc(N*B->length*sizeof(ulong));
        mpoly_repack_monomials(Bexps, Qbits, B->exps, B->bits, B->length, ctx->minfo);
    }

    /* check leading monomial divides exactly */
    if (Qbits <= FLINT_BITS)
    {
        /* mask with high bit of each exponent vector field set */
        for (i = 0; i < FLINT_BITS/Qbits; i++)
            mask = (mask << Qbits) + (UWORD(1) << (Qbits - 1));

        if (!mpoly_monomial_divides(expq, Aexps, Bexps, N, mask))
        {
            success = 0;
            goto cleanup;
        }
    }
    else
    {
        if (!mpoly_monomial_divides_mp(expq, Aexps, Bexps, N, Qbits))
        {
            success = 0;
            goto cleanup;
        }
    }

    /* deal with aliasing and divide polynomials */
    if (Q == A || Q == B)
    {
        fq_nmod_mpoly_t T;
        fq_nmod_mpoly_init(T, ctx);
        fq_nmod_mpoly_fit_length_reset_bits(T, A->length/B->length + 1, Qbits, ctx);
        success = _fq_nmod_mpoly_divides_monagan_pearce(T,
                                          A->coeffs, Aexps, A->length,
                                          B->coeffs, Bexps, B->length,
                                                Qbits, N, cmpmask, ctx->fqctx);
        fq_nmod_mpoly_swap(Q, T, ctx);
        fq_nmod_mpoly_clear(T, ctx);
    }
    else
    {
        fq_nmod_mpoly_fit_length_reset_bits(Q, A->length/B->length + 1, Qbits, ctx);
        success = _fq_nmod_mpoly_divides_monagan_pearce(Q,
                                          A->coeffs, Aexps, A->length,
                                          B->coeffs, Bexps, B->length,
                                                Qbits, N, cmpmask, ctx->fqctx);
   }

cleanup:

    if (freeAexps)
        flint_free(Aexps);

    if (freeBexps)
        flint_free(Bexps);

    TMP_END;

    return success;
}
