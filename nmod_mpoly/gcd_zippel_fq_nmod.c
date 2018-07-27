/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "nmod_mpoly.h"


int fq_nmod_next_not_zero(fq_nmod_t a, fq_nmod_t b, fq_nmod_ctx_t ctx)
{
    assert(0);
    return 0;
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
    assert(0);
    return 0;
}


void fq_nmod_mpolyun_set_mpolyu(
    fq_nmod_mpolyun_t A,
    const fq_nmod_mpolyu_t B,
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








void fq_nmod_mpoly_cvtfrom_mpolyn(
    fq_nmod_mpoly_t A,
    fq_nmod_mpolyn_t B,
    slong var,
    fq_nmod_mpoly_ctx_t ctx)
{
    slong i, j;
    slong k;
    slong N;
    ulong * oneexp;
    slong offset;
    slong shift;
    TMP_INIT;

    FLINT_ASSERT(B->bits == A->bits);
    FLINT_ASSERT(B->bits <= FLINT_BITS);
    FLINT_ASSERT(ctx->minfo->ord == ORD_LEX);

    TMP_START;

    N = mpoly_words_per_exp(B->bits, ctx->minfo);
    oneexp = TMP_ALLOC(N*sizeof(ulong));
    mpoly_gen_oneexp_offset_shift(oneexp, &offset, &shift, var, N, B->bits, ctx->minfo);

    fq_nmod_mpoly_fit_length(A, B->length, ctx);

    k = 0;
    for (i = 0; i < B->length; i++)
    {
        for (j = (B->coeffs + i)->length - 1; j >= 0; j--)
        {
            if (!fq_nmod_is_zero((B->coeffs + i)->coeffs + j, ctx->fqctx))
            {
                fq_nmod_mpoly_fit_length(A, k + 1, ctx);
                fq_nmod_set(A->coeffs + k, (B->coeffs + i)->coeffs + j, ctx->fqctx);
                mpoly_monomial_madd(A->exps + N*k, B->exps + N*i, j, oneexp, N);                
                k++;
            }
        }
    }

    A->length = k;
    TMP_END;
}
void fq_nmod_mpolyu_cvtfrom_mpolyun(
    fq_nmod_mpolyu_t A,
    fq_nmod_mpolyun_t B,
    slong var,
    fq_nmod_mpoly_ctx_t ctx)
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




void fq_nmod_mpoly_cvtto_mpolyn(
    fq_nmod_mpolyn_t A,
    fq_nmod_mpoly_t B,
    slong var,
    fq_nmod_mpoly_ctx_t ctx)
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

    N = mpoly_words_per_exp(B->bits, ctx->minfo);
    oneexp = TMP_ALLOC(N*sizeof(ulong));
    mask = (-UWORD(1)) >> (FLINT_BITS - B->bits);
    mpoly_gen_oneexp_offset_shift(oneexp, &offset, &shift, var, N, B->bits, ctx->minfo);

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
            fq_nmod_poly_set_coeff(A->coeffs + k - 1, c, B->coeffs + i, ctx->fqctx);
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
void fq_nmod_mpolyu_cvtto_mpolyun(
    fq_nmod_mpolyun_t A,
    fq_nmod_mpolyu_t B,
    slong k,
    fq_nmod_mpoly_ctx_t ctx)
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








/*
    Try to set G to the gcd of A and B given the form f of G.
    return codes as enumerated in nmod_mpoly.h:

    nmod_mpoly_sgcd_success,
    nmod_mpoly_sgcd_form_wrong,
    nmod_mpoly_sgcd_no_solution,
    nmod_mpoly_sgcd_scales_not_found,
    nmod_mpoly_sgcd_eval_point_not_found,
    nmod_mpoly_sgcd_eval_gcd_deg_too_high
*/
nmod_sgcd_ret_t fq_nmod_mpolyu_sgcd_zippel(
    fq_nmod_mpolyu_t G,
    fq_nmod_mpolyu_t A,
    fq_nmod_mpolyu_t B,
    fq_nmod_mpolyu_t f,
    slong var,
    fq_nmod_mpoly_ctx_t ctx,
    flint_rand_t randstate,
    slong * degbound)
{
    assert(0);
    return 0;
}




void fq_nmod_mpolyn_eval_last(fq_nmod_mpoly_t B, fq_nmod_mpolyn_t A, fq_nmod_t alpha, fq_nmod_mpoly_ctx_t ctx)
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
        fq_nmod_poly_evaluate_fq_nmod(B->coeffs + k, A->coeffs + i, alpha, ctx->fqctx);
        if (!fq_nmod_is_zero(B->coeffs + k, ctx->fqctx))
        {
            k++;
        }
    }
    B->length = k;
}

void fq_nmod_mpolyun_eval_last(fq_nmod_mpolyu_t B, fq_nmod_mpolyun_t A, fq_nmod_t alpha, fq_nmod_mpoly_ctx_t ctx)
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










void fq_nmod_mpoly_setform(fq_nmod_mpoly_t A, fq_nmod_mpoly_t B, fq_nmod_mpoly_ctx_t ctx)
{
    slong i;

    fq_nmod_mpoly_set(A, B, ctx);
    for (i = 0; i < A->length; i++)
    {
        fq_nmod_zero(A->coeffs + i, ctx->fqctx);
    }
}

void fq_nmod_mpolyu_setform(fq_nmod_mpolyu_t A, fq_nmod_mpolyu_t B, fq_nmod_mpoly_ctx_t ctx)
{
    slong i;

    fq_nmod_mpolyu_fit_length(A, B->length, ctx);
    for (i = 0; i < B->length; i++)
    {
        fq_nmod_mpoly_setform(A->coeffs + i, B->coeffs + i, ctx);
        A->exps[i] = B->exps[i];
    }
    A->length = B->length;
}


void fq_nmod_mpoly_setform_mpolyn(fq_nmod_mpoly_t A, fq_nmod_mpolyn_t B, fq_nmod_mpoly_ctx_t ctx)
{
    slong i;
    slong N;

    FLINT_ASSERT(A->bits == B->bits);

    fq_nmod_mpoly_fit_length(A, B->length, ctx);
    N = mpoly_words_per_exp(B->bits, ctx->minfo);
    for (i = 0; i < B->length; i++)
    {
        fq_nmod_zero(A->coeffs + i, ctx->fqctx);
        mpoly_monomial_set(A->exps + N*i, B->exps + N*i, N);
    }
    A->length = B->length;
}

void fq_nmod_mpolyu_setform_mpolyun(fq_nmod_mpolyu_t A, fq_nmod_mpolyun_t B, fq_nmod_mpoly_ctx_t ctx)
{
    slong i;

    fq_nmod_mpolyu_fit_length(A, B->length, ctx);
    for (i = 0; i < B->length; i++)
    {
        FLINT_ASSERT((B->coeffs + i)->bits == B->bits);
        fq_nmod_mpoly_setform_mpolyn(A->coeffs + i, B->coeffs + i, ctx);
        A->exps[i] = B->exps[i];
    }
    A->length = B->length;
}

int fq_nmod_mpolyu_pgcd_zippel_univar(
    fq_nmod_mpolyu_t G,
    fq_nmod_mpolyu_t A,
    fq_nmod_mpolyu_t B,
    fq_nmod_mpoly_ctx_t ctx)
{
    fq_nmod_poly_t a, b, g;
    FLINT_ASSERT(A->bits == B->bits);
    fq_nmod_poly_init(a, ctx->fqctx);
    fq_nmod_poly_init(b, ctx->fqctx);
    fq_nmod_poly_init(g, ctx->fqctx);
    fq_nmod_mpolyu_cvtto_poly(a, A, ctx);
    fq_nmod_mpolyu_cvtto_poly(b, B, ctx);
    fq_nmod_poly_gcd(g, a, b, ctx->fqctx);
    fq_nmod_mpolyu_cvtfrom_poly(G, g, ctx);
    fq_nmod_poly_clear(a, ctx->fqctx);
    fq_nmod_poly_clear(b, ctx->fqctx);
    fq_nmod_poly_clear(g, ctx->fqctx);
    return 1;
}


int fq_nmod_mpolyu_pgcd_zippel_bivar(
    fq_nmod_mpolyu_t G,
    fq_nmod_mpolyu_t A,
    fq_nmod_mpolyu_t B,
    fq_nmod_mpoly_ctx_t ctx,
    mpoly_zipinfo_t zinfo,
    flint_rand_t randstate)
{
    fq_nmod_poly_t a, b, c, g, modulus, tempmod;
    fq_nmod_mpolyu_t Aeval, Beval, Geval;
    fq_nmod_mpolyun_t An, Bn, H, Ht;
    fq_nmod_t geval, temp, alpha, alphastart;
    fmpz_t minusone;
    int success = 0, changed;
    slong Alastdeg;
    slong Blastdeg;
    slong lastdeg;
    slong var = 0;

    FLINT_ASSERT(ctx->minfo->ord == ORD_LEX);
    FLINT_ASSERT(var >= -WORD(1));

    FLINT_ASSERT(G->bits == A->bits);
    FLINT_ASSERT(G->bits == B->bits);

    fmpz_init(minusone);
    fmpz_set_si(minusone, -WORD(1));

    fq_nmod_mpolyun_init(An, A->bits, ctx);
    fq_nmod_mpolyun_init(Bn, A->bits, ctx);

    fq_nmod_mpolyu_cvtto_mpolyun(An, A, var, ctx);
    fq_nmod_mpolyu_cvtto_mpolyun(Bn, B, var, ctx);

    fq_nmod_poly_init(a, ctx->fqctx);
    fq_nmod_poly_init(b, ctx->fqctx);
    fq_nmod_poly_init(c, ctx->fqctx);
    fq_nmod_poly_init(g, ctx->fqctx);
    fq_nmod_mpolyun_content_last(a, An, ctx);
    fq_nmod_mpolyun_content_last(b, Bn, ctx);
    fq_nmod_mpolyun_divexact_last(An, a, ctx);
    fq_nmod_mpolyun_divexact_last(Bn, b, ctx);
    fq_nmod_poly_gcd(c, a, b, ctx->fqctx);
    fq_nmod_poly_gcd(g, fq_nmod_mpolyun_leadcoeff_ref(An, ctx),
                        fq_nmod_mpolyun_leadcoeff_ref(Bn, ctx), ctx->fqctx);

    Alastdeg = fq_nmod_mpolyun_lastdeg(An, ctx);
    Blastdeg = fq_nmod_mpolyun_lastdeg(Bn, ctx);

    fq_nmod_poly_init(modulus, ctx->fqctx);
    fq_nmod_poly_init(tempmod, ctx->fqctx);
    fq_nmod_poly_set_coeff_fmpz(tempmod, 1, minusone, ctx->fqctx);
    fq_nmod_mpolyu_init(Aeval, A->bits, ctx);
    fq_nmod_mpolyu_init(Beval, A->bits, ctx);
    fq_nmod_mpolyu_init(Geval, A->bits, ctx);
    fq_nmod_mpolyun_init(H, A->bits, ctx);
    fq_nmod_mpolyun_init(Ht, A->bits, ctx);

    fq_nmod_poly_one(modulus, ctx->fqctx);
    fq_nmod_mpolyun_zero(H, ctx);

    fq_nmod_randtest_not_zero(alphastart, randstate, ctx->fqctx);
    fq_nmod_set(alpha, alphastart, ctx->fqctx);
    while (1)
    {
        /* get new evaluation point */
        if (fq_nmod_next_not_zero(alpha, alphastart, ctx->fqctx) == 0)
        {
            success = 0;
            goto finished;
        }

        /* make sure evaluation point does not kill both lc(A) and lc(B) */
        fq_nmod_poly_evaluate_fq_nmod(geval, g, alpha, ctx->fqctx);
        if (fq_nmod_is_zero(geval, ctx->fqctx))
            goto outer_continue;

        /* make sure evaluation point does not kill either A or B */
        fq_nmod_mpolyun_eval_last(Aeval, An, alpha, ctx);
        fq_nmod_mpolyun_eval_last(Beval, Bn, alpha, ctx);
        if (Aeval->length == 0 || Beval->length == 0)
            goto outer_continue;

        fq_nmod_mpolyu_pgcd_zippel_univar(Geval, Aeval, Beval, ctx);

        if (fq_nmod_mpolyu_is_one(Geval, ctx))
        {
            fq_nmod_mpolyu_cvtfrom_poly_notmain(G, c, var, ctx);
            success = 1;
            goto finished;
        }

        FLINT_ASSERT(Geval->length > 0);

        if (fq_nmod_poly_degree(modulus, ctx->fqctx) > 0)
        {
            if (Geval->exps[0] > H->exps[0])
            {
                goto outer_continue;

            } else if (Geval->exps[0] < H->exps[0])
            {
                fq_nmod_poly_one(modulus, ctx->fqctx);                
            }
        }

        /* update interpolant H */
        fq_nmod_inv(temp, fq_nmod_mpolyu_leadcoeff_ref(Geval, ctx), ctx->fqctx);
        fq_nmod_mul(temp, geval, temp, ctx->fqctx);
        fq_nmod_mpolyu_scalar_mul_fq_nmod(Geval, temp, ctx);
        if (fq_nmod_poly_degree(modulus, ctx->fqctx) > 0)
        {
            fq_nmod_poly_evaluate_fq_nmod(temp, modulus, alpha, ctx->fqctx);
            fq_nmod_inv(temp, temp, ctx->fqctx);
            fq_nmod_poly_scalar_mul_fq_nmod(modulus, modulus, temp, ctx->fqctx);
            changed = fq_nmod_mpolyun_addinterp(&lastdeg, H, Ht, Geval, modulus, alpha, ctx);
            if (!changed)
            {
                fq_nmod_mpolyun_content_last(a, H, ctx);
                fq_nmod_mpolyun_mul_poly(Ht, H, c, ctx);
                fq_nmod_mpolyun_divexact_last(Ht, a, ctx);

                fq_nmod_mpolyu_cvtfrom_mpolyun(G, Ht, var, ctx);

                if (   !fq_nmod_mpolyu_divides(A, G, ctx)
                    || !fq_nmod_mpolyu_divides(B, G, ctx))
                {
                    goto outer_continue;
                }
                success = 1;
                goto finished;
            }

        } else
        {
            fq_nmod_mpolyun_set_mpolyu(H, Geval, ctx);
            lastdeg = WORD(0);
        }
        fq_nmod_poly_set_coeff(tempmod, 0, alpha, ctx->fqctx);
        fq_nmod_poly_mul(modulus, modulus, tempmod, ctx->fqctx);

        /* something is wrong if the interpolation degree is too high */
        if (lastdeg > Alastdeg || lastdeg > Blastdeg)
        {
            fq_nmod_poly_one(modulus, ctx->fqctx);
            goto outer_continue;
        }

outer_continue:
        NULL;
    }

    success = 1;

finished:

    fmpz_clear(minusone);
    fq_nmod_clear(temp, ctx->fqctx);
    fq_nmod_clear(alpha, ctx->fqctx);
    fq_nmod_clear(alphastart, ctx->fqctx);
    fq_nmod_poly_clear(a, ctx->fqctx);
    fq_nmod_poly_clear(b, ctx->fqctx);
    fq_nmod_poly_clear(c, ctx->fqctx);
    fq_nmod_poly_clear(g, ctx->fqctx);
    fq_nmod_poly_clear(modulus, ctx->fqctx);
    fq_nmod_poly_clear(tempmod, ctx->fqctx);
    fq_nmod_mpolyu_clear(Aeval, ctx);
    fq_nmod_mpolyu_clear(Beval, ctx);
    fq_nmod_mpolyu_clear(Geval, ctx);
    fq_nmod_mpolyun_clear(An, ctx);
    fq_nmod_mpolyun_clear(Bn, ctx);
    fq_nmod_mpolyun_clear(H, ctx);
    fq_nmod_mpolyun_clear(Ht, ctx);

    return success;
}



int fq_nmod_mpolyu_pgcd_zippel(
    fq_nmod_mpolyu_t G,
    fq_nmod_mpolyu_t A,
    fq_nmod_mpolyu_t B,
    slong var,
    fq_nmod_mpoly_ctx_t ctx,
    mpoly_zipinfo_t zinfo,
    flint_rand_t randstate)
{
    int divcheck_fail_count;
    slong lastdeg;
    slong Alastdeg;
    slong Blastdeg;
    ulong ABminshift;
    slong degbound;
    int success = 0, changed;
    fq_nmod_mpolyun_t An, Bn;
    fq_nmod_poly_t a, b, c, g;
    fq_nmod_poly_t modulus, tempmod;
    fq_nmod_mpolyu_t Aeval, Beval, Geval, Gform;
    fq_nmod_mpolyun_t H, Ht;
    fq_nmod_t geval, temp;
    fq_nmod_t alpha, alphastart;
    fmpz_t minusone;

    FLINT_ASSERT(ctx->minfo->ord == ORD_LEX);
    FLINT_ASSERT(var >= -WORD(1));

    FLINT_ASSERT(G->bits == A->bits);
    FLINT_ASSERT(G->bits == B->bits);
    FLINT_ASSERT(A->length > 0);
    FLINT_ASSERT(B->length > 0);

    if (var == -WORD(1))
    {
        /* no more variables left to interpolate */
        return fq_nmod_mpolyu_pgcd_zippel_univar(G, A, B, ctx);
    }

    fq_nmod_mpolyun_init(An, A->bits, ctx);
    fq_nmod_mpolyun_init(Bn, A->bits, ctx);

    fq_nmod_mpolyu_cvtto_mpolyun(An, A, var, ctx);
    fq_nmod_mpolyu_cvtto_mpolyun(Bn, B, var, ctx);
    ABminshift = FLINT_MIN(A->exps[A->length - 1], B->exps[B->length - 1]);
    fq_nmod_mpolyun_shift_right(An, A->exps[A->length - 1]);
    fq_nmod_mpolyun_shift_right(Bn, B->exps[B->length - 1]);
    Alastdeg = fq_nmod_mpolyun_lastdeg(An, ctx);
    Blastdeg = fq_nmod_mpolyun_lastdeg(Bn, ctx);

    fq_nmod_poly_init(a, ctx->fqctx);
    fq_nmod_poly_init(b, ctx->fqctx);
    fq_nmod_poly_init(c, ctx->fqctx);
    fq_nmod_poly_init(g, ctx->fqctx);

    fq_nmod_mpolyun_content_last(a, An, ctx);
    fq_nmod_mpolyun_content_last(b, Bn, ctx);

    /* if the gcd has content wrt last variable, we are going to fail */
    /*nmod_mpolyun_divexact_last(An, a, ctx);*/
    /*nmod_mpolyun_divexact_last(Bn, b, ctx);*/

    fq_nmod_poly_gcd(c, a, b, ctx->fqctx);
    fq_nmod_poly_gcd(g, fq_nmod_mpolyun_leadcoeff_ref(An, ctx),
                        fq_nmod_mpolyun_leadcoeff_ref(Bn, ctx), ctx->fqctx);

    fq_nmod_poly_init(modulus, ctx->fqctx);
    fq_nmod_poly_init(tempmod, ctx->fqctx);
    fmpz_init(minusone);
    fmpz_set_si(minusone, WORD(-1));
    fq_nmod_poly_set_coeff_fmpz(tempmod, 1, minusone, ctx->fqctx);
    fq_nmod_mpolyu_init(Aeval, A->bits, ctx);
    fq_nmod_mpolyu_init(Beval, A->bits, ctx);
    fq_nmod_mpolyu_init(Geval, A->bits, ctx);
    fq_nmod_mpolyu_init(Gform, A->bits, ctx);
    fq_nmod_mpolyun_init(H, A->bits, ctx);
    fq_nmod_mpolyun_init(Ht, A->bits, ctx);

    if (fq_nmod_poly_degree(c, ctx->fqctx) > 0)
    {
        success = 0;
        goto finished;
    }

    if (var == WORD(0))
    {
        /* bivariate is more comfortable separated */
        success = fq_nmod_mpolyu_pgcd_zippel_bivar(G, A, B, ctx, zinfo, randstate);
        goto finished;
    }

    if (nmod_poly_degree(ctx->fqctx->modulus) < WORD(2))
    {
        success = 0;
        goto finished;
    }


    degbound = FLINT_MIN(A->exps[0], B->exps[0]);

    fq_nmod_poly_one(modulus, ctx->fqctx);
    fq_nmod_mpolyun_zero(H, ctx);

    fq_nmod_randtest_not_zero(alphastart, randstate, ctx->fqctx);
    fq_nmod_set(alpha, alphastart, ctx->fqctx);
    while (1)
    {
        /* get new evaluation point */
        if (fq_nmod_next_not_zero(alpha, alphastart, ctx->fqctx) == 0)
        {
            success = 0;
            goto finished;
        }

        /* make sure evaluation point does not kill both lc(A) and lc(B) */
        fq_nmod_poly_evaluate_fq_nmod(geval, g, alpha, ctx->fqctx);
        if (fq_nmod_is_zero(geval, ctx->fqctx))
            goto outer_continue;

        /* make sure evaluation point does not kill either A or B */
        fq_nmod_mpolyun_eval_last(Aeval, An, alpha, ctx);
        fq_nmod_mpolyun_eval_last(Beval, Bn, alpha, ctx);
        if (Aeval->length == 0 || Beval->length == 0)
            goto outer_continue;

        success = fq_nmod_mpolyu_pgcd_zippel(Geval, Aeval, Beval, var - 1, ctx, zinfo, randstate);
        if (!success || Geval->exps[0] > degbound)
        {
            success = 0;
            goto finished;
        }
        
        degbound = Geval->exps[0];

        if (fq_nmod_mpolyu_is_one(Geval, ctx))
        {
            fq_nmod_mpolyu_cvtfrom_poly_notmain(G, c, var, ctx);
            fq_nmod_mpolyu_shift_left(G, ABminshift);
            success = 1;
            goto finished;
        }

        if (fq_nmod_poly_degree(modulus, ctx->fqctx) > 0)
        {
            if (Geval->exps[0] > H->exps[0])
            {
                goto outer_continue;

            } else if (Geval->exps[0] < H->exps[0])
            {
                fq_nmod_poly_one(modulus, ctx->fqctx);                
            }
        }

        /* update interpolant H */
        fq_nmod_inv(temp, fq_nmod_mpolyu_leadcoeff_ref(Geval, ctx), ctx->fqctx);
        fq_nmod_mul(temp, geval, temp, ctx->fqctx);
        fq_nmod_mpolyu_scalar_mul_fq_nmod(Geval, temp, ctx);
        if (fq_nmod_poly_degree(modulus, ctx->fqctx) > 0)
        {
            fq_nmod_poly_evaluate_fq_nmod(temp, modulus, alpha, ctx->fqctx);
            fq_nmod_inv(temp, temp, ctx->fqctx);
            fq_nmod_poly_scalar_mul_fq_nmod(modulus, modulus, temp, ctx->fqctx);
            changed = fq_nmod_mpolyun_addinterp(&lastdeg, H, Ht, Geval, modulus, alpha, ctx);
            if (!changed)
            {
                fq_nmod_mpolyun_content_last(a, H, ctx);
                fq_nmod_mpolyun_mul_poly(Ht, H, c, ctx);
                fq_nmod_mpolyun_divexact_last(Ht, a, ctx);
                fq_nmod_mpolyun_shift_left(Ht, ABminshift);
                fq_nmod_mpolyu_cvtfrom_mpolyun(G, Ht, var, ctx);
                if (    !fq_nmod_mpolyu_divides(A, G, ctx)
                     || !fq_nmod_mpolyu_divides(B, G, ctx))
                {
                    fq_nmod_poly_one(modulus, ctx->fqctx);
                    goto outer_continue;
                }
                success = 1;
                goto finished;
            }

        } else
        {
            fq_nmod_mpolyun_set_mpolyu(H, Geval, ctx);
        }
        fq_nmod_poly_set_coeff(tempmod, 0, alpha, ctx->fqctx);
        fq_nmod_poly_mul(modulus, modulus, tempmod, ctx->fqctx);

        fq_nmod_mpolyu_setform_mpolyun(Gform, H, ctx);

        divcheck_fail_count = 0;
        while (1)
        {
            /* get new evaluation point */
            if (fq_nmod_next_not_zero(alpha, alphastart, ctx->fqctx) == 0)
            {
                success = 0;
                goto finished;
            }

            /* make sure evaluation point does not kill both lc(A) and lc(B) */
            fq_nmod_poly_evaluate_fq_nmod(geval, g, alpha, ctx->fqctx);
            if (fq_nmod_is_zero(geval, ctx->fqctx))
                goto inner_continue;

            /* make sure evaluation point does not kill either A or B */
            fq_nmod_mpolyun_eval_last(Aeval, An, alpha, ctx);
            fq_nmod_mpolyun_eval_last(Beval, Bn, alpha, ctx);
            if (Aeval->length == 0 || Beval->length == 0)
                goto inner_continue;

            switch (fq_nmod_mpolyu_sgcd_zippel(Geval, Aeval, Beval, Gform, var,
                                                    ctx, randstate, &degbound))
            {
                default:
                    FLINT_ASSERT(0);
                case nmod_sgcd_form_main_degree_too_high:
                    /* nmod_mpolyu_sgcd_zippel has updated degbound */
                    fq_nmod_poly_one(modulus, ctx->fqctx);
                    goto outer_continue;
                case nmod_sgcd_form_wrong:
                case nmod_sgcd_no_solution:
                    success = 0;
                    goto finished;
                case nmod_sgcd_scales_not_found:
                case nmod_sgcd_eval_point_not_found:
                case nmod_sgcd_eval_gcd_deg_too_high:
                    goto inner_continue;
                case nmod_sgcd_success:
                    (void)(NULL);
            }

            if (fq_nmod_is_zero(fq_nmod_mpolyu_leadcoeff_ref(Geval, ctx), ctx->fqctx))
                goto inner_continue;

            /* update interpolant H */
            FLINT_ASSERT(fq_nmod_poly_degree(modulus, ctx->fqctx) > 0);

            fq_nmod_poly_evaluate_fq_nmod(temp, modulus, alpha, ctx->fqctx);
            fq_nmod_inv(temp, temp, ctx->fqctx);
            fq_nmod_poly_scalar_mul_fq_nmod(modulus, modulus, temp, ctx->fqctx);
            changed = fq_nmod_mpolyun_addinterp(&lastdeg, H, Ht, Geval, modulus, alpha, ctx);

            fq_nmod_poly_set_coeff(tempmod, 0, alpha, ctx->fqctx);
            fq_nmod_poly_mul(modulus, modulus, tempmod, ctx->fqctx);

            if (!changed)
            {
                fq_nmod_mpolyun_content_last(a, H, ctx);
                fq_nmod_mpolyun_mul_poly(Ht, H, c, ctx);
                fq_nmod_mpolyun_divexact_last(Ht, a, ctx);
                fq_nmod_mpolyun_shift_left(Ht, ABminshift);
                fq_nmod_mpolyu_cvtfrom_mpolyun(G, Ht, var, ctx);
                if (    !fq_nmod_mpolyu_divides(A, G, ctx)
                     || !fq_nmod_mpolyu_divides(B, G, ctx))
                {
                    ++divcheck_fail_count;
                    if (divcheck_fail_count >= 2)
                    {
                        fq_nmod_poly_one(modulus, ctx->fqctx);
                        goto outer_continue;
                    }
                    goto inner_continue;
                }
                success = 1;
                goto finished;
            }

            /* something is wrong if the interpolation degree is too high */
            if (lastdeg > Alastdeg || lastdeg > Blastdeg)
            {
                fq_nmod_poly_one(modulus, ctx->fqctx);
                goto outer_continue;
            }
inner_continue:
            NULL;
        }
        FLINT_ASSERT(0 && "not reachable");

outer_continue:
        NULL;
    }
    FLINT_ASSERT(0 && "not reachable");


finished:

    fmpz_clear(minusone);
    fq_nmod_poly_clear(a, ctx->fqctx);
    fq_nmod_poly_clear(b, ctx->fqctx);
    fq_nmod_poly_clear(c, ctx->fqctx);
    fq_nmod_poly_clear(g, ctx->fqctx);
    fq_nmod_poly_clear(modulus, ctx->fqctx);
    fq_nmod_poly_clear(tempmod, ctx->fqctx);
    fq_nmod_mpolyu_clear(Aeval, ctx);
    fq_nmod_mpolyu_clear(Beval, ctx);
    fq_nmod_mpolyu_clear(Geval, ctx);
    fq_nmod_mpolyu_clear(Gform, ctx);
    fq_nmod_mpolyun_clear(An, ctx);
    fq_nmod_mpolyun_clear(Bn, ctx);
    fq_nmod_mpolyun_clear(H, ctx);
    fq_nmod_mpolyun_clear(Ht, ctx);

    return success;
}
