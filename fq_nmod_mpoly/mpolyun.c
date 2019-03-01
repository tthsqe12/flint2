/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "nmod_mpoly.h"
#include "fq_nmod_mpoly.h"

void fq_nmod_mpolyun_init(fq_nmod_mpolyun_t A, mp_bitcnt_t bits,
                                                 const fq_nmod_mpoly_ctx_t ctx)
{
    A->coeffs = NULL;
    A->exps = NULL;
    A->alloc = 0;
    A->length = 0;
    A->bits = bits;
}

void fq_nmod_mpolyun_clear(fq_nmod_mpolyun_t A, const fq_nmod_mpoly_ctx_t ctx)
{
    slong i;
    for (i = 0; i < A->alloc; i++)
        fq_nmod_mpolyn_clear(A->coeffs + i, ctx);
    flint_free(A->coeffs);
    flint_free(A->exps);
}


void fq_nmod_mpolyun_swap(fq_nmod_mpolyun_t A, fq_nmod_mpolyun_t B)
{
   fq_nmod_mpolyun_struct t = *A;
   *A = *B;
   *B = t;
}

void fq_nmod_mpolyun_zero(fq_nmod_mpolyun_t A, const fq_nmod_mpoly_ctx_t ctx)
{
    slong i;
    for (i = 0; i < A->alloc; i++) {
        fq_nmod_mpolyn_clear(A->coeffs + i, ctx);
        fq_nmod_mpolyn_init(A->coeffs + i, A->bits, ctx);
    }
    A->length = 0;
}

void fq_nmod_mpolyun_print_pretty(const fq_nmod_mpolyun_t poly,
                                const char ** x, const fq_nmod_mpoly_ctx_t ctx)
{
    slong i;
    if (poly->length == 0)
        flint_printf("0");
    for (i = 0; i < poly->length; i++)
    {
        if (i != 0)
            flint_printf(" + ");
        flint_printf("(");
        FLINT_ASSERT((poly->coeffs + i)->bits == poly->bits);
        fq_nmod_mpolyn_print_pretty(poly->coeffs + i,x,ctx);
        flint_printf(")*X^%wd",poly->exps[i]);
    }
}

void fq_nmod_mpolyun_fit_length(fq_nmod_mpolyun_t A, slong length,
                                                 const fq_nmod_mpoly_ctx_t ctx)
{
    slong i;
    slong old_alloc = A->alloc;
    slong new_alloc = FLINT_MAX(length, 2*A->alloc);

    if (length > old_alloc)
    {
        if (old_alloc == 0)
        {
            A->exps = (ulong *) flint_malloc(new_alloc*sizeof(ulong));
            A->coeffs = (fq_nmod_mpolyn_struct *) flint_malloc(
                                      new_alloc*sizeof(fq_nmod_mpolyn_struct));
        } else
        {
            A->exps = (ulong *) flint_realloc(A->exps,
                                                      new_alloc*sizeof(ulong));
            A->coeffs = (fq_nmod_mpolyn_struct *) flint_realloc(A->coeffs,
                                      new_alloc*sizeof(fq_nmod_mpolyn_struct));
        }

        for (i = old_alloc; i < new_alloc; i++)
        {
            fq_nmod_mpolyn_init(A->coeffs + i, A->bits, ctx);
        }
        A->alloc = new_alloc;
    }
}

void fq_nmod_mpolyun_shift_right(fq_nmod_mpolyun_t A, ulong s)
{
    slong i;
    for (i = 0; i < A->length; i++)
    {
        FLINT_ASSERT(A->exps[i] >= s);
        A->exps[i] -= s;
    }
}

void fq_nmod_mpolyun_shift_left(fq_nmod_mpolyun_t A, ulong s)
{
    slong i;
    for (i = 0; i < A->length; i++)
    {
        A->exps[i] += s;
    }
}

slong fq_nmod_mpolyun_lastdeg(fq_nmod_mpolyun_t A, const fq_nmod_mpoly_ctx_t ctx)
{
    slong i, j;
    slong deg = -WORD(1);

    for (i = 0; i < A->length; i++)
    {
        for (j = 0; j < (A->coeffs + i)->length; j++)
        {
            deg = FLINT_MAX(deg,
                 fq_nmod_poly_degree((A->coeffs + i)->coeffs + j, ctx->fqctx));
        }
    }
    FLINT_ASSERT(deg >= 0);
    return deg;
}

void fq_nmod_mpolyun_set(fq_nmod_mpolyun_t A, const fq_nmod_mpolyun_t B,
                                                 const fq_nmod_mpoly_ctx_t ctx)
{
    slong i, Blen;
    fq_nmod_mpolyn_struct * Acoeff, * Bcoeff;
    ulong * Aexp, * Bexp;

    Blen = B->length;
    fq_nmod_mpolyun_fit_length(A, Blen, ctx);
    Acoeff = A->coeffs;
    Bcoeff = B->coeffs;
    Aexp = A->exps;
    Bexp = B->exps;

    for (i = 0; i < Blen; i++)
    {
        fq_nmod_mpolyn_set(Acoeff + i, Bcoeff + i, ctx);
        Aexp[i] = Bexp[i];
    }

    /* demote remaining coefficients */
    for (i = Blen; i < A->length; i++)
    {
        fq_nmod_mpolyn_clear(Acoeff + i, ctx);
        fq_nmod_mpolyn_init(Acoeff + i, A->bits, ctx);
    }
    A->length = Blen;
}

void fq_nmod_mpolyun_one(fq_nmod_mpolyun_t A, const fq_nmod_mpoly_ctx_t ctx)
{
    fq_nmod_mpolyun_fit_length(A, 1, ctx);
    fq_nmod_mpolyn_one(A->coeffs + 0, ctx);
    A->exps[0] = 0;
    A->length = 1;
}

/* take the last variable of B out */
void fq_nmod_mpoly_cvtto_mpolyn(fq_nmod_mpolyn_t A, fq_nmod_mpoly_t B,
                                      slong var, const fq_nmod_mpoly_ctx_t ctx)
{
    slong i;
    slong k;
    ulong * oneexp;
    slong offset;
    slong shift;
    ulong mask;
    slong N;
    TMP_INIT;

    FLINT_ASSERT(B->bits <= FLINT_BITS);
    FLINT_ASSERT(ctx->minfo->ord == ORD_LEX);

    TMP_START;

    N = mpoly_words_per_exp_sp(B->bits, ctx->minfo);
    oneexp = TMP_ALLOC(N*sizeof(ulong));
    mask = (-UWORD(1)) >> (FLINT_BITS - B->bits);
    mpoly_gen_monomial_offset_shift_sp(oneexp, &offset, &shift, var,
                                                          B->bits, ctx->minfo);

    fq_nmod_mpolyn_fit_bits(A, B->bits, ctx);
    A->bits = B->bits;

    k = 0;
    fq_nmod_mpolyn_fit_length(A, k + 1, ctx);
    for (i = 0; i < B->length; i++)
    {
        ulong c = (B->exps[N*i + offset] >> shift) & mask;
        mpoly_monomial_msub(A->exps + N*k, B->exps + N*i, c, oneexp, N);

        if (k > 0 && mpoly_monomial_equal(A->exps + N*k, A->exps + N*(k - 1), N))
        {
            fq_nmod_poly_set_coeff(A->coeffs + k - 1, c, B->coeffs + i,
                                                                   ctx->fqctx);
        } else
        {
            fq_nmod_poly_zero(A->coeffs + k, ctx->fqctx);
            fq_nmod_poly_set_coeff(A->coeffs + k, c, B->coeffs + i, ctx->fqctx);
            k++;
            fq_nmod_mpolyn_fit_length(A, k + 1, ctx);
        }
    }

    fq_nmod_mpolyn_set_length(A, k, ctx);
    TMP_END;
}

void fq_nmod_mpolyu_cvtto_mpolyun(fq_nmod_mpolyun_t A,
                    fq_nmod_mpolyu_t B, slong k, const fq_nmod_mpoly_ctx_t ctx)
{
    slong i, Blen;
    fq_nmod_mpolyn_struct * Acoeff;
    fq_nmod_mpoly_struct * Bcoeff;
    ulong * Aexp, * Bexp;

    Blen = B->length;
    fq_nmod_mpolyun_fit_length(A, Blen, ctx);
    Acoeff = A->coeffs;
    Bcoeff = B->coeffs;
    Aexp = A->exps;
    Bexp = B->exps;

    for (i = 0; i < Blen; i++)
    {
        fq_nmod_mpoly_cvtto_mpolyn(Acoeff + i, Bcoeff + i, k, ctx);
        Aexp[i] = Bexp[i];
    }

    /* demote remaining coefficients */
    for (i = Blen; i < A->length; i++)
    {
        fq_nmod_mpolyn_clear(Acoeff + i, ctx);
        fq_nmod_mpolyn_init(Acoeff + i, A->bits, ctx);
    }
    A->length = Blen;  
}


/* put the last variable of B back into A */
void fq_nmod_mpoly_cvtfrom_mpolyn(fq_nmod_mpoly_t A,
                  fq_nmod_mpolyn_t B, slong var, const fq_nmod_mpoly_ctx_t ctx)
{
    slong i, j;
    slong k;
    slong N;
    ulong * oneexp;
    TMP_INIT;

    FLINT_ASSERT(B->bits == A->bits);
    FLINT_ASSERT(B->bits <= FLINT_BITS);
    FLINT_ASSERT(ctx->minfo->ord == ORD_LEX);

    TMP_START;

    N = mpoly_words_per_exp_sp(B->bits, ctx->minfo);
    oneexp = TMP_ALLOC(N*sizeof(ulong));
    mpoly_gen_monomial_sp(oneexp, var, B->bits, ctx->minfo);

    fq_nmod_mpoly_fit_length(A, B->length, ctx);

    k = 0;
    for (i = 0; i < B->length; i++)
    {
        for (j = (B->coeffs + i)->length - 1; j >= 0; j--)
        {
            if (!fq_nmod_is_zero((B->coeffs + i)->coeffs + j, ctx->fqctx))
            {
                fq_nmod_mpoly_fit_length(A, k + 1, ctx);
                fq_nmod_set(A->coeffs + k, (B->coeffs + i)->coeffs + j,
                                                                   ctx->fqctx);
                mpoly_monomial_madd(A->exps + N*k, B->exps + N*i, j, oneexp, N);                
                k++;
            }
        }
    }

    A->length = k;
    TMP_END;
}

void fq_nmod_mpolyu_cvtfrom_mpolyun(fq_nmod_mpolyu_t A,
                 fq_nmod_mpolyun_t B, slong var, const fq_nmod_mpoly_ctx_t ctx)
{
    slong i;

    fq_nmod_mpolyu_fit_length(A, B->length, ctx);

    for (i = 0; i < B->length; i++)
    {
        fq_nmod_mpoly_cvtfrom_mpolyn(A->coeffs + i, B->coeffs + i, var, ctx);
        A->exps[i] = B->exps[i];
    }
    A->length = B->length;
}


void fq_nmod_mpolyn_set_mpoly(fq_nmod_mpolyn_t A, const fq_nmod_mpoly_t B,
                                                 const fq_nmod_mpoly_ctx_t ctx)
{
    slong i;
    slong N;
    fq_nmod_poly_struct * Acoeff;
    fq_nmod_struct * Bcoeff;
    ulong * Aexp, * Bexp;
    slong Blen;

    fq_nmod_mpolyn_fit_bits(A, B->bits, ctx);
    A->bits = B->bits;

    Blen = B->length;
    fq_nmod_mpolyn_fit_length(A, Blen, ctx);
    Acoeff = A->coeffs;
    Bcoeff = B->coeffs;
    Aexp = A->exps;
    Bexp = B->exps;

    N = mpoly_words_per_exp(B->bits, ctx->minfo);

    for (i = 0; i < Blen; i++)
    {
        fq_nmod_poly_zero(Acoeff + i, ctx->fqctx);
        fq_nmod_poly_set_coeff(Acoeff + i, 0, Bcoeff + i, ctx->fqctx);
        mpoly_monomial_set(Aexp + N*i, Bexp + N*i, N);
    }

    /* demote remaining coefficients */
    for (i = Blen; i < A->length; i++)
    {
        fq_nmod_poly_clear(Acoeff + i, ctx->fqctx);
        fq_nmod_poly_init(Acoeff + i, ctx->fqctx);
    }
    A->length = Blen;
}

void fq_nmod_mpolyun_set_mpolyu(fq_nmod_mpolyun_t A, const fq_nmod_mpolyu_t B,
                                                 const fq_nmod_mpoly_ctx_t ctx)
{
    slong i;

    FLINT_ASSERT(A->bits == B->bits);
    fq_nmod_mpolyun_fit_length(A, B->length, ctx);
    for (i = 0; i < B->length; i++)
    {
        A->exps[i] = B->exps[i];
        fq_nmod_mpolyn_set_mpoly(A->coeffs + i, B->coeffs + i, ctx);

        FLINT_ASSERT((A->coeffs + i)->bits == B->bits);

    }
    A->length = B->length;
}


void fq_nmod_mpolyun_mul_poly(fq_nmod_mpolyun_t A, const fq_nmod_mpolyun_t B,
                         const fq_nmod_poly_t c, const fq_nmod_mpoly_ctx_t ctx)
{
    slong i, Blen;
    fq_nmod_mpolyn_struct * Acoeff, * Bcoeff;
    ulong * Aexp, * Bexp;

    Blen = B->length;
    fq_nmod_mpolyun_fit_length(A, Blen, ctx);
    Acoeff = A->coeffs;
    Bcoeff = B->coeffs;
    Aexp = A->exps;
    Bexp = B->exps;

    for (i = 0; i < Blen; i++)
    {
        fq_nmod_mpolyn_mul_poly(Acoeff + i, Bcoeff + i, c, ctx);
        Aexp[i] = Bexp[i];
    }

    /* demote remaining coefficients */
    for (i = Blen; i < A->length; i++)
    {
        fq_nmod_mpolyn_clear(Acoeff + i, ctx);
        fq_nmod_mpolyn_init(Acoeff + i, A->bits, ctx);
    }
    A->length = Blen;
}


void fq_nmod_mpolyun_content_last(fq_nmod_poly_t a, fq_nmod_mpolyun_t B,
                                                 const fq_nmod_mpoly_ctx_t ctx)
{
    slong i, j;

    fq_nmod_poly_zero(a, ctx->fqctx);
    for (i = 0; i < B->length; i++)
    {
        for (j = 0; j < (B->coeffs + i)->length; j++)
        {
            fq_nmod_poly_gcd(a, a, (B->coeffs + i)->coeffs + j, ctx->fqctx);
            if (fq_nmod_poly_degree(a, ctx->fqctx) == 0)
                break;
        }
    }
}

void fq_nmod_mpolyun_divexact_last(fq_nmod_mpolyun_t A, fq_nmod_poly_t b,
                                                 const fq_nmod_mpoly_ctx_t ctx)
{
    slong i, j;

    fq_nmod_poly_t r;
    fq_nmod_poly_init(r, ctx->fqctx);
    for (i = 0; i < A->length; i++)
    {
        for (j = 0; j < (A->coeffs + i)->length; j++)
        {
            fq_nmod_poly_divrem((A->coeffs + i)->coeffs + j, r,
                                (A->coeffs + i)->coeffs + j, b, ctx->fqctx);
            FLINT_ASSERT(fq_nmod_poly_is_zero(r, ctx->fqctx));
        }
    }
    fq_nmod_poly_clear(r, ctx->fqctx);
}


void fq_nmod_mpolyun_mul_last(fq_nmod_mpolyun_t A, fq_nmod_poly_t b,
                                                 const fq_nmod_mpoly_ctx_t ctx)
{
    slong i, j;

    for (i = 0; i < A->length; i++)
    {
        for (j = 0; j < (A->coeffs + i)->length; j++)
        {
            fq_nmod_poly_mul((A->coeffs + i)->coeffs + j, 
                             (A->coeffs + i)->coeffs + j, b, ctx->fqctx);
        }
    }
}

/*
    Update H so that it does not change mod m, and is now A mod p
    It is asserted that the monomials in H and A match
*/
int fq_nmod_mpolyn_CRT_fq_nmod_mpoly(slong * lastdeg,
                          fq_nmod_mpolyn_t H, const fq_nmod_mpoly_ctx_t ctx,
                          fq_nmod_poly_t m, fq_nmod_t inv_m_eval,
                             fq_nmod_mpoly_t A, const fq_nmod_mpoly_ctx_t ectx,
                                                    const _fq_nmod_embed_t emb)
{
    slong i;
#if WANT_ASSERT
    slong N;
#endif
    int changed = 0;
    fq_nmod_t u, v;
    fq_nmod_poly_t w, u_sm;

    fq_nmod_init(u, ectx->fqctx);
    fq_nmod_init(v, ectx->fqctx);
    fq_nmod_poly_init(w, ctx->fqctx);
    fq_nmod_poly_init(u_sm, ctx->fqctx);

    FLINT_ASSERT(H->length == A->length);
    FLINT_ASSERT(H->bits == A->bits);
#if WANT_ASSERT
    N = mpoly_words_per_exp(A->bits, ctx->minfo);
#endif
    for (i = 0; i < A->length; i++)
    {
        FLINT_ASSERT(mpoly_monomial_equal(H->exps + N*i, A->exps + N*i, N));
        _fq_nmod_embed_sm_to_lg(u, H->coeffs + i, emb);
        fq_nmod_sub(v, A->coeffs + i, u, ectx->fqctx);
        if (!fq_nmod_is_zero(v, ectx->fqctx))
        {
            changed = 1;
            fq_nmod_mul(u, v, inv_m_eval, ectx->fqctx);
            _fq_nmod_embed_lg_to_sm(u_sm, u, emb);
            fq_nmod_poly_mul(w, u_sm, m, ctx->fqctx);
            fq_nmod_poly_add(H->coeffs + i, H->coeffs + i, w, ctx->fqctx);
        }

        lastdeg[0] = FLINT_MAX(lastdeg[0],
                               fq_nmod_poly_degree(H->coeffs + i, ctx->fqctx));

        FLINT_ASSERT(fq_nmod_poly_degree(H->coeffs + i, ctx->fqctx)
                         <  fq_nmod_poly_degree(m, ctx->fqctx)
                          + fq_nmod_poly_degree(emb->h, ctx->fqctx));
    }

    fq_nmod_clear(u, ectx->fqctx);
    fq_nmod_clear(v, ectx->fqctx);
    fq_nmod_poly_clear(w, ctx->fqctx);
    fq_nmod_poly_clear(u_sm, ctx->fqctx);

    return changed;
}

int fq_nmod_mpolyun_CRT_fq_nmod_mpolyu(slong * lastdeg,
         fq_nmod_mpolyun_t H, const fq_nmod_mpoly_ctx_t ctx, fq_nmod_poly_t m,
                        fq_nmod_mpolyu_t A, const fq_nmod_mpoly_ctx_t ectx,
                                                          _fq_nmod_embed_t emb)
{
    slong i;
    int changed = 0;
    fq_nmod_t inv_m_eval;

    lastdeg[0] = -WORD(1);

    fq_nmod_init(inv_m_eval, ectx->fqctx);
    _fq_nmod_embed_sm_to_lg(inv_m_eval, m, emb);
    fq_nmod_inv(inv_m_eval, inv_m_eval, ectx->fqctx);

    FLINT_ASSERT(H->bits == A->bits);
    FLINT_ASSERT(H->length == A->length);
    for (i = 0; i < A->length; i++)
    {
        FLINT_ASSERT(H->exps[i] == A->exps[i]);
        changed |= fq_nmod_mpolyn_CRT_fq_nmod_mpoly(lastdeg, H->coeffs + i,
                                 ctx, m, inv_m_eval, A->coeffs + i, ectx, emb);
    }
    H->length = A->length;
    fq_nmod_clear(inv_m_eval, ectx->fqctx);
    return changed;
}



/* reduce B via the map F_q[x] -> F_q^n */
void fq_nmod_mpolyn_redto_fq_nmod_mpoly(fq_nmod_mpoly_t A, fq_nmod_mpolyn_t B,
                 const fq_nmod_mpoly_ctx_t ectx, const fq_nmod_mpoly_ctx_t ctx,
                                                    const _fq_nmod_embed_t emb)
{
    slong i;
    slong k;
    slong N;

    FLINT_ASSERT(A->bits == B->bits);
    FLINT_ASSERT(B->bits <= FLINT_BITS);
    FLINT_ASSERT(ctx->minfo->ord == ORD_LEX);
    FLINT_ASSERT(ectx->minfo->ord == ORD_LEX);
    FLINT_ASSERT(ectx->minfo->nvars == ctx->minfo->nvars);

    N = mpoly_words_per_exp_sp(B->bits, ctx->minfo);

    k = 0;
    fq_nmod_mpoly_fit_length(A, k + 1, ectx);
    for (i = 0; i < B->length; i++)
    {
        fq_nmod_mpoly_fit_length(A, k + 1, ectx);
        mpoly_monomial_set(A->exps + N*k, B->exps + N*i, N);
        _fq_nmod_embed_sm_to_lg(A->coeffs + k, B->coeffs + i, emb);
        k += !fq_nmod_is_zero(A->coeffs + k, ectx->fqctx);
    }

    _fq_nmod_mpoly_set_length(A, k, ectx);
}

/* reduce B via the map from F_q[x] -> F_q^n */
void fq_nmod_mpolyun_redto_fq_nmod_mpolyu(fq_nmod_mpolyu_t A, fq_nmod_mpolyun_t B,
                 const fq_nmod_mpoly_ctx_t ectx, const fq_nmod_mpoly_ctx_t ctx,
                                                    const _fq_nmod_embed_t emb)
{
    slong i, k, Blen;
    fq_nmod_mpoly_struct * Acoeff;
    fq_nmod_mpolyn_struct * Bcoeff;
    ulong * Aexp, * Bexp;

    Blen = B->length;
    fq_nmod_mpolyu_fit_length(A, Blen, ectx);
    Acoeff = A->coeffs;
    Bcoeff = B->coeffs;
    Aexp = A->exps;
    Bexp = B->exps;

    k = 0;
    for (i = 0; i < Blen; i++)
    {
        fq_nmod_mpolyn_redto_fq_nmod_mpoly(Acoeff + k, Bcoeff + i, ectx, ctx, emb);
        Aexp[k] = Bexp[i];
        k += !fq_nmod_mpoly_is_zero(Acoeff + k, ectx);
    }
    A->length = k;  
}



/* evaluate A at lastvar = alpha */
void fq_nmod_mpolyn_eval_last(fq_nmod_mpoly_t B, fq_nmod_mpolyn_t A,
                                fq_nmod_t alpha, const fq_nmod_mpoly_ctx_t ctx)
{
    slong i;
    slong N;
    slong k;
    FLINT_ASSERT(B->bits == A->bits);

    fq_nmod_mpoly_fit_length(B, A->length, ctx);
    N = mpoly_words_per_exp(A->bits, ctx->minfo);
    k = 0;
    for (i = 0; i < A->length; i++)
    {
        mpoly_monomial_set(B->exps + N*k, A->exps + N*i, N);
        fq_nmod_poly_evaluate_fq_nmod(B->coeffs + k, A->coeffs + i, alpha,
                                                                   ctx->fqctx);
        if (!fq_nmod_is_zero(B->coeffs + k, ctx->fqctx))
        {
            k++;
        }
    }
    B->length = k;
}

void fq_nmod_mpolyun_eval_last(fq_nmod_mpolyu_t B, fq_nmod_mpolyun_t A,
                                fq_nmod_t alpha, const fq_nmod_mpoly_ctx_t ctx)
{
    slong i;
    slong k;

    FLINT_ASSERT(B->bits == A->bits);

    fq_nmod_mpolyu_fit_length(B, A->length, ctx);
    k = 0;
    for (i = 0; i < A->length; i++)
    {
        B->exps[k] = A->exps[i];
        fq_nmod_mpolyn_eval_last(B->coeffs + k, A->coeffs + i, alpha, ctx);
        if (!fq_nmod_mpoly_is_zero(B->coeffs + k, ctx))
        {
            k++;
        }
    }
    B->length = k;
}


/*
    F = F + modulus*(A - F(alpha))
*/
int
fq_nmod_mpolyn_addinterp(slong * lastdeg,
        fq_nmod_mpolyn_t F, fq_nmod_mpolyn_t T, fq_nmod_mpoly_t A,
        fq_nmod_poly_t modulus, fq_nmod_t alpha, const fq_nmod_mpoly_ctx_t ctx)
{
    int changed = 0;
    slong i, j, k;
    slong N;
    fq_nmod_t v;
    mp_bitcnt_t bits = A->bits;
    slong Flen = F->length, Alen = A->length;
    ulong * Fexp = F->exps, * Aexp = A->exps;
    ulong * Texp;
    fq_nmod_struct * Acoeff = A->coeffs;
    fq_nmod_poly_struct * Fcoeff = F->coeffs;
    fq_nmod_poly_struct * Tcoeff;
    fq_nmod_poly_t tp;

    FLINT_ASSERT(F->bits == bits);
    FLINT_ASSERT(ctx->minfo->ord == ORD_LEX);

    fq_nmod_init(v, ctx->fqctx);
    fq_nmod_poly_init(tp, ctx->fqctx);

    fq_nmod_mpolyn_fit_length(T, Flen + Alen, ctx);
    Texp = T->exps;
    Tcoeff = T->coeffs;

    N = mpoly_words_per_exp(bits, ctx->minfo);

    i = j = k = 0;
    while (i < Flen || j < Alen)
    {
        if (i < Flen && (j >= Alen
                        || mpoly_monomial_gt_nomask(Fexp + N*i, Aexp + N*j, N)))
        {
            FLINT_ASSERT(!fq_nmod_poly_is_zero(Fcoeff + i, ctx->fqctx));
            FLINT_ASSERT(fq_nmod_poly_degree(Fcoeff + i, ctx->fqctx)
                                   < fq_nmod_poly_degree(modulus, ctx->fqctx));

            /* F term ok, A term missing */
            fq_nmod_poly_evaluate_fq_nmod(v, Fcoeff + i, alpha, ctx->fqctx);
            if (!fq_nmod_is_zero(v, ctx->fqctx))
            {
                changed = 1;
                fq_nmod_poly_scalar_mul_fq_nmod(tp, modulus, v, ctx->fqctx);
                fq_nmod_poly_sub(Tcoeff + k, Fcoeff + i, tp, ctx->fqctx);
            } else {
                fq_nmod_poly_set(Tcoeff + k, Fcoeff + i, ctx->fqctx);                
            }
            lastdeg[0] = FLINT_MAX(lastdeg[0],
                                  fq_nmod_poly_degree(Tcoeff + k, ctx->fqctx));

            mpoly_monomial_set(Texp + N*k, Fexp + N*i, N);
            FLINT_ASSERT(!fq_nmod_poly_is_zero(Tcoeff + k, ctx->fqctx));
            k++;
            i++;
        }
        else if (j < Alen && (i >= Flen
                        || mpoly_monomial_lt_nomask(Fexp + N*i, Aexp + N*j, N)))
        {
            /* F term missing, A term ok */
            if (!fq_nmod_is_zero(Acoeff + j, ctx->fqctx))
            {
                changed = 1;
                fq_nmod_poly_zero(Tcoeff + k, ctx->fqctx);
                fq_nmod_poly_scalar_mul_fq_nmod(Tcoeff + k, modulus, Acoeff + j,
                                                                   ctx->fqctx);
                lastdeg[0] = FLINT_MAX(lastdeg[0],
                                  fq_nmod_poly_degree(Tcoeff + k, ctx->fqctx));
                mpoly_monomial_set(Texp + N*k, Aexp + N*j, N);
                k++;
            }
            j++;
        }
        else if (i < Flen && j < Alen
                             && mpoly_monomial_equal(Fexp + N*i, Aexp + N*j, N))
        {
            FLINT_ASSERT(!fq_nmod_poly_is_zero(Fcoeff + i, ctx->fqctx));
            FLINT_ASSERT(fq_nmod_poly_degree(Fcoeff + i, ctx->fqctx) < fq_nmod_poly_degree(modulus, ctx->fqctx));

            /* F term ok, A term ok */
            fq_nmod_poly_evaluate_fq_nmod(v, Fcoeff + i, alpha, ctx->fqctx);
            fq_nmod_sub(v, Acoeff + j, v, ctx->fqctx);
            if (!fq_nmod_is_zero(v, ctx->fqctx))
            {
                changed = 1;
                fq_nmod_poly_scalar_mul_fq_nmod(tp, modulus, v, ctx->fqctx);
                fq_nmod_poly_add(Tcoeff + k, Fcoeff + i, tp, ctx->fqctx);
            } else {
                fq_nmod_poly_set(Tcoeff + k, Fcoeff + i, ctx->fqctx);
            }
            lastdeg[0] = FLINT_MAX(lastdeg[0],
                                  fq_nmod_poly_degree(Tcoeff + k, ctx->fqctx));
            mpoly_monomial_set(Texp + N*k, Aexp + N*j, N);
            FLINT_ASSERT(!fq_nmod_poly_is_zero(Tcoeff + k, ctx->fqctx));
            k++;
            i++;
            j++;
        } else 
        {
            FLINT_ASSERT(0);
        }
    }

    fq_nmod_mpolyn_set_length(T, k, ctx);

    if (changed)
    {
        fq_nmod_mpolyn_swap(T, F);
    }

    fq_nmod_poly_clear(tp, ctx->fqctx);
    fq_nmod_clear(v, ctx->fqctx);

    return changed;
}


int fq_nmod_mpolyun_addinterp(
    slong * lastdeg,
    fq_nmod_mpolyun_t F,
    fq_nmod_mpolyun_t T,
    fq_nmod_mpolyu_t A,
    fq_nmod_poly_t modulus,
    fq_nmod_t alpha,
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

    lastdeg[0] = -WORD(1);

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
    fq_nmod_mpoly_fit_bits(zero, A->bits, ctx);
    zero->bits = A->bits;

    i = j = k = 0;
    while (i < Flen || j < Alen)
    {
        if (i < Flen && (j >= Alen || Fexp[i] > Aexp[j]))
        {
            /* F term ok, A term missing */
            fq_nmod_mpolyn_set(Tcoeff + k, Fcoeff + i, ctx);
            changed |= fq_nmod_mpolyn_addinterp(lastdeg, Tcoeff + k, S, zero, modulus, alpha, ctx);
            Texp[k] = Fexp[i];
            k++;
            i++;
        }
        else if (j < Alen && (i >= Flen || Aexp[j] > Fexp[i]))
        {
            /* F term missing, A term ok */
            fq_nmod_mpolyn_zero(Tcoeff + k, ctx);
            changed |= fq_nmod_mpolyn_addinterp(lastdeg, Tcoeff + k, S, Acoeff + j, modulus, alpha, ctx);
            Texp[k] = Aexp[j];
            k++;
            j++;
        }
        else if (i < Flen && j < Alen && (Fexp[i] == Aexp[j]))
        {
            /* F term ok, A term ok */
            fq_nmod_mpolyn_set(Tcoeff + k, Fcoeff + i, ctx);
            changed |= fq_nmod_mpolyn_addinterp(lastdeg, Tcoeff + k, S, Acoeff + j, modulus, alpha, ctx);
            Texp[k] = Aexp[j];
            FLINT_ASSERT(!fq_nmod_mpolyn_is_zero(Tcoeff + k, ctx));
            k++;
            i++;
            j++;
        } else 
        {
            FLINT_ASSERT(0);
        }
    }

    /* demote remaining coefficients */
    for (i = k; i < T->length; i++)
    {
        fq_nmod_mpolyn_clear(Tcoeff + i, ctx);
        fq_nmod_mpolyn_init(Tcoeff + i, T->bits, ctx);
    }
    T->length = k;

    if (changed)
    {
        fq_nmod_mpolyun_swap(T, F);
    }

    fq_nmod_mpolyn_clear(S, ctx);
    fq_nmod_mpoly_clear(zero, ctx);
    return changed;    
}



/* Convert B to A using the map  F_q^n -> F_q[x] */
void fq_nmod_mpolyn_startinterp_lgprime(slong * lastdeg,
                            fq_nmod_mpolyn_t A, const fq_nmod_mpoly_ctx_t ctx,
                            fq_nmod_mpoly_t B, const fq_nmod_mpoly_ctx_t ectx,
                                                    const _fq_nmod_embed_t emb)
{
    slong i, N;

    FLINT_ASSERT(B->bits == A->bits);
    fq_nmod_mpolyn_fit_length(A, B->length, ctx);
    N = mpoly_words_per_exp(A->bits, ctx->minfo);
    for (i = 0; i < B->length; i++)
    {
        mpoly_monomial_set(A->exps + N*i, B->exps + N*i, N);
        _fq_nmod_embed_lg_to_sm(A->coeffs + i, B->coeffs + i, emb);
        FLINT_ASSERT(!fq_nmod_poly_is_zero(A->coeffs + i, ctx->fqctx));
        lastdeg[0] = FLINT_MAX(lastdeg[0],
                               fq_nmod_poly_degree(A->coeffs + i, ctx->fqctx));

    }
    A->length = B->length;
}

void fq_nmod_mpolyun_startinterp_lgprime(slong * lastdeg,
                           fq_nmod_mpolyun_t A, const fq_nmod_mpoly_ctx_t ctx,
                           fq_nmod_mpolyu_t B, const fq_nmod_mpoly_ctx_t ectx,
                                                    const _fq_nmod_embed_t emb)
{
    slong i;

    lastdeg[0] = -WORD(1);

    FLINT_ASSERT(B->bits == A->bits);
    fq_nmod_mpolyun_fit_length(A, B->length, ctx);
    for (i = 0; i < B->length; i++)
    {
        A->exps[i] = B->exps[i];
        fq_nmod_mpolyn_startinterp_lgprime(lastdeg, A->coeffs + i, ctx, B->coeffs + i, ectx, emb);
        FLINT_ASSERT(!fq_nmod_mpolyn_is_zero(A->coeffs + i, ctx));
    }
    A->length = B->length;
}



/*
    update F so that it doesn't change mod m and is A mod emb->h

    F = F + m*((A - F(alpha))/(m(alpha)))

    no assumptions about matching monomials
*/
int fq_nmod_mpolyn_addinterp_lgprime(slong * lastdeg,
                     fq_nmod_mpolyn_t F, fq_nmod_mpolyn_t T, fq_nmod_poly_t m,
                         const fq_nmod_mpoly_ctx_t ctx, fq_nmod_mpoly_t A,
                         fq_nmod_t inv_m_eval, const fq_nmod_mpoly_ctx_t ectx,
                                                    const _fq_nmod_embed_t emb)
{
    int changed = 0;
    slong i, j, k;
    slong N;
    fq_nmod_t u, v;
    fq_nmod_poly_t u_sm, w;
    mp_bitcnt_t bits = A->bits;
    slong Flen = F->length, Alen = A->length;
    ulong * Fexp = F->exps, * Aexp = A->exps;
    ulong * Texp;
    fq_nmod_struct * Acoeff = A->coeffs;
    fq_nmod_poly_struct * Fcoeff = F->coeffs;
    fq_nmod_poly_struct * Tcoeff;

    FLINT_ASSERT(bits <= FLINT_BITS);
    FLINT_ASSERT(F->bits == bits);
    FLINT_ASSERT(ctx->minfo->ord == ORD_LEX);

    fq_nmod_init(u, ectx->fqctx);
    fq_nmod_init(v, ectx->fqctx);
    fq_nmod_poly_init(u_sm, ctx->fqctx);
    fq_nmod_poly_init(w, ctx->fqctx);

    fq_nmod_mpolyn_fit_length(T, Flen + Alen, ctx);
    Texp = T->exps;
    Tcoeff = T->coeffs;

    N = mpoly_words_per_exp(bits, ctx->minfo);

    i = j = k = 0;
    while (i < Flen || j < Alen)
    {
        if (i < Flen && (j >= Alen
                        || mpoly_monomial_gt_nomask(Fexp + N*i, Aexp + N*j, N)))
        {
            FLINT_ASSERT(!fq_nmod_poly_is_zero(Fcoeff + i, ctx->fqctx));
            FLINT_ASSERT(fq_nmod_poly_degree(Fcoeff + i, ctx->fqctx)
                                  < fq_nmod_poly_degree(m, ctx->fqctx));

            /* F term ok, A term missing */
            _fq_nmod_embed_sm_to_lg(v, Fcoeff + i, emb);
            if (!fq_nmod_is_zero(v, ectx->fqctx))
            {
                changed = 1;
                fq_nmod_mul(u, v, inv_m_eval, ectx->fqctx);
                _fq_nmod_embed_lg_to_sm(u_sm, u, emb);
                fq_nmod_poly_mul(w, u_sm, m, ctx->fqctx);
                fq_nmod_poly_sub(Tcoeff + k, Fcoeff + i, w, ctx->fqctx);
            }
            else
            {
                fq_nmod_poly_set(Tcoeff + k, Fcoeff + i, ctx->fqctx);
            }
            lastdeg[0] = FLINT_MAX(lastdeg[0], fq_nmod_poly_degree(Tcoeff + k, ctx->fqctx));

            mpoly_monomial_set(Texp + N*k, Fexp + N*i, N);
            FLINT_ASSERT(!fq_nmod_poly_is_zero(Tcoeff + k, ctx->fqctx));
            k++;
            i++;
        }
        else if (j < Alen && (i >= Flen
                        || mpoly_monomial_lt_nomask(Fexp + N*i, Aexp + N*j, N)))
        {
            /* F term missing, A term ok */
            if (!fq_nmod_is_zero(Acoeff + j, ectx->fqctx))
            {
                changed = 1;
                fq_nmod_mul(u, Acoeff + j, inv_m_eval, ectx->fqctx);
                _fq_nmod_embed_lg_to_sm(u_sm, u, emb);
                fq_nmod_poly_mul(Tcoeff + k, m, u_sm, ctx->fqctx);
                lastdeg[0] = FLINT_MAX(lastdeg[0], fq_nmod_poly_degree(Tcoeff + k, ctx->fqctx));
                mpoly_monomial_set(Texp + N*k, Aexp + N*j, N);
                k++;
            }
            j++;
        }
        else if (i < Flen && j < Alen
                             && mpoly_monomial_equal(Fexp + N*i, Aexp + N*j, N))
        {
            FLINT_ASSERT(!fq_nmod_poly_is_zero(Fcoeff + i, ctx->fqctx));
            FLINT_ASSERT(fq_nmod_poly_degree(Fcoeff + i, ctx->fqctx)
                                        < fq_nmod_poly_degree(m, ctx->fqctx));

            /* F term ok, A term ok */
            _fq_nmod_embed_sm_to_lg(u, Fcoeff + i, emb);
            fq_nmod_sub(v, Acoeff + j, u, ectx->fqctx);
            if (!fq_nmod_is_zero(v, ectx->fqctx))
            {
                changed = 1;
                fq_nmod_mul(u, v, inv_m_eval, ectx->fqctx);
                _fq_nmod_embed_lg_to_sm(u_sm, u, emb);
                fq_nmod_poly_mul(w, m, u_sm, ctx->fqctx);
                fq_nmod_poly_add(Tcoeff + k, Fcoeff + i, w, ctx->fqctx);
            }
            else
            {
                fq_nmod_poly_set(Tcoeff + k, Fcoeff + i, ctx->fqctx);
            }
            lastdeg[0] = FLINT_MAX(lastdeg[0], fq_nmod_poly_degree(Tcoeff + k, ctx->fqctx));
            mpoly_monomial_set(Texp + N*k, Aexp + N*j, N);
            FLINT_ASSERT(!fq_nmod_poly_is_zero(Tcoeff + k, ctx->fqctx));
            k++;
            i++;
            j++;
        }
        else 
        {
            FLINT_ASSERT(0);
        }
    }

    fq_nmod_mpolyn_set_length(T, k, ctx);

    if (changed)
    {
        fq_nmod_mpolyn_swap(T, F);
    }

    fq_nmod_clear(u, ectx->fqctx);
    fq_nmod_clear(v, ectx->fqctx);
    fq_nmod_poly_clear(u_sm, ctx->fqctx);
    fq_nmod_poly_clear(w, ctx->fqctx);

    return changed;
}


int fq_nmod_mpolyun_addinterp_lgprime(slong * lastdeg,
                   fq_nmod_mpolyun_t F, fq_nmod_mpolyun_t T, fq_nmod_poly_t m,
                       const fq_nmod_mpoly_ctx_t ctx, fq_nmod_mpolyu_t A,
                    const fq_nmod_mpoly_ctx_t ectx, const _fq_nmod_embed_t emb)
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
    fq_nmod_t inv_m_eval;

    lastdeg[0] = -WORD(1);

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

    fq_nmod_mpoly_init(zero, ectx);
    fq_nmod_mpoly_fit_bits(zero, A->bits, ectx);
    zero->bits = A->bits;

    fq_nmod_init(inv_m_eval, ectx->fqctx);
    _fq_nmod_embed_sm_to_lg(inv_m_eval, m, emb);
    fq_nmod_inv(inv_m_eval, inv_m_eval, ectx->fqctx);

    i = j = k = 0;
    while (i < Flen || j < Alen)
    {
        if (i < Flen && (j >= Alen || Fexp[i] > Aexp[j]))
        {
            /* F term ok, A term missing */
            fq_nmod_mpolyn_set(Tcoeff + k, Fcoeff + i, ctx);
            changed |= fq_nmod_mpolyn_addinterp_lgprime(lastdeg, Tcoeff + k,
                                       S, m, ctx, zero, inv_m_eval, ectx, emb);
            Texp[k] = Fexp[i];
            FLINT_ASSERT(!fq_nmod_mpolyn_is_zero(Tcoeff + k, ctx));
            k++;
            i++;
        }
        else if (j < Alen && (i >= Flen || Aexp[j] > Fexp[i]))
        {
            /* F term missing, A term ok */
            fq_nmod_mpolyn_zero(Tcoeff + k, ctx);
            changed |= fq_nmod_mpolyn_addinterp_lgprime(lastdeg, Tcoeff + k,
                                 S, m, ctx, Acoeff + j, inv_m_eval, ectx, emb);
            Texp[k] = Aexp[j];
            FLINT_ASSERT(!fq_nmod_mpolyn_is_zero(Tcoeff + k, ctx));
            k++;
            j++;
        }
        else if (i < Flen && j < Alen && (Fexp[i] == Aexp[j]))
        {
            /* F term ok, A term ok */
            fq_nmod_mpolyn_set(Tcoeff + k, Fcoeff + i, ctx);
            changed |= fq_nmod_mpolyn_addinterp_lgprime(lastdeg, Tcoeff + k,
                                 S, m, ctx, Acoeff + j, inv_m_eval, ectx, emb);
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
    {
        fq_nmod_mpolyun_swap(T, F);
    }

    fq_nmod_clear(inv_m_eval, ectx->fqctx);

    fq_nmod_mpolyn_clear(S, ctx);
    fq_nmod_mpoly_clear(zero, ectx);
    return changed;    
}
