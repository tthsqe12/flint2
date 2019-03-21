/*
    Copyright (C) 2019 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fq_nmod_mpoly.h"
int usleep(ulong usec);

int fq_nmod_mpolyun_is_canonical(const fq_nmod_mpolyun_t A, const fq_nmod_mpoly_ctx_t ctx);
int fq_nmod_mpolyun_is_nonzero_fq_nmod(const fq_nmod_mpolyun_t A, const fq_nmod_mpoly_ctx_t ctx);
fq_nmod_struct * fq_nmod_mpolyn_leadcoeff_last_ref(fq_nmod_mpolyn_t A, const fq_nmod_mpoly_ctx_t ctx);
void fq_nmod_mpolyun_scalar_mul_fq_nmod(fq_nmod_mpolyun_t A, const fq_nmod_t c, const fq_nmod_mpoly_ctx_t ctx);
int fq_nmod_mpolyun_gcd_brown_smprime(fq_nmod_mpolyun_t G,
          fq_nmod_mpolyun_t Abar, fq_nmod_mpolyun_t Bbar, fq_nmod_mpolyun_t A,
                  fq_nmod_mpolyun_t B, slong var, const fq_nmod_mpoly_ctx_t ctx);

/*
    Try to set G to the gcd of A and B using Brown's alogrithm M.
    This function switches to a big primes version if needed.
    It should only really fail if the dense size of the inputs is too large.
*/
int fq_nmod_mpoly_gcd_brown(fq_nmod_mpoly_t G,
                            const fq_nmod_mpoly_t A, const fq_nmod_mpoly_t B,
                                                 const fq_nmod_mpoly_ctx_t ctx)
{
    int success;
    fq_nmod_mpolyd_t Ad, Bd, Gd, Abar, Bbar;
    fq_nmod_mpolyd_ctx_t dctx;
    slong nvars = ctx->minfo->nvars;

    success = 1;

    if (fq_nmod_mpoly_is_zero(A, ctx))
    {
        fq_nmod_mpoly_set(G, B, ctx);
        goto cleanup_stage0;
    }
    if (fq_nmod_mpoly_is_zero(B, ctx))
    {
        fq_nmod_mpoly_set(G, A, ctx);
        goto cleanup_stage0;
    }

    fq_nmod_mpolyd_ctx_init2(dctx, nvars, ctx->fqctx);
    success = fq_nmod_mpolyd_ctx_set_for_gcd(dctx, A, B, ctx);
    if (!success)
    {
        fq_nmod_mpoly_zero(G, ctx);
        goto cleanup_stage1;
    }

    fq_nmod_mpolyd_init(Ad, nvars, ctx->fqctx);
    fq_nmod_mpolyd_init(Bd, nvars, ctx->fqctx);
    fq_nmod_mpolyd_init(Gd, nvars, ctx->fqctx);
    fq_nmod_mpolyd_init(Abar, nvars, ctx->fqctx);
    fq_nmod_mpolyd_init(Bbar, nvars, ctx->fqctx);

    fq_nmod_mpoly_convert_to_fq_nmod_mpolyd(Ad, dctx, A, ctx);
    fq_nmod_mpoly_convert_to_fq_nmod_mpolyd(Bd, dctx, B, ctx);
    success = fq_nmod_mpolyd_gcd_brown_smprime(Gd, Abar, Bbar, Ad, Bd, dctx);
    if (!success)
    {
        fq_nmod_mpoly_convert_to_fq_nmod_mpolyd(Ad, dctx, A, ctx);
        fq_nmod_mpoly_convert_to_fq_nmod_mpolyd(Bd, dctx, B, ctx);
        success = fq_nmod_mpolyd_gcd_brown_lgprime(Gd, Abar, Bbar, Ad, Bd, dctx);
    }
    if (success)
    {
        fq_nmod_mpoly_convert_from_fq_nmod_mpolyd(G, ctx, Gd, dctx);
    }

    fq_nmod_mpolyd_clear(Bbar, ctx->fqctx);
    fq_nmod_mpolyd_clear(Abar, ctx->fqctx);
    fq_nmod_mpolyd_clear(Gd, ctx->fqctx);
    fq_nmod_mpolyd_clear(Bd, ctx->fqctx);
    fq_nmod_mpolyd_clear(Ad, ctx->fqctx);

cleanup_stage1:

    fq_nmod_mpolyd_ctx_clear(dctx);

cleanup_stage0:

    if (success && !fq_nmod_mpoly_is_zero(G, ctx))
    {
        fq_nmod_mpoly_make_monic(G, G, ctx);
    }

    return success;
}


void fq_nmod_mpoly_to_mpolyun_perm_deflate(fq_nmod_mpolyun_t A, const fq_nmod_mpoly_t B,
               const slong * perm, const ulong * shift, const ulong * stride,
                       const fq_nmod_mpoly_ctx_t uctx, const fq_nmod_mpoly_ctx_t ctx)
{
    slong m = uctx->minfo->nvars;
    fq_nmod_mpolyu_t Au;
    fq_nmod_mpolyu_init(Au, A->bits, uctx);
    fq_nmod_mpoly_to_mpolyu_perm_deflate(Au, B, perm, shift, stride, uctx, ctx);
    fq_nmod_mpolyu_cvtto_mpolyun(A, Au, m - 1, uctx);
    fq_nmod_mpolyu_clear(Au, uctx);
    return;
}



void fq_nmod_mpoly_from_mpolyun_perm_inflate(fq_nmod_mpoly_t A, mp_bitcnt_t Abits,
                                                         fq_nmod_mpolyun_t B,
                const slong * perm, const ulong * shift, const ulong * stride,
                       const fq_nmod_mpoly_ctx_t uctx, const fq_nmod_mpoly_ctx_t ctx)
{
    slong m = uctx->minfo->nvars;
    fq_nmod_mpolyu_t Au;
    fq_nmod_mpolyu_init(Au, B->bits, uctx);
    fq_nmod_mpolyu_cvtfrom_mpolyun(Au, B, m - 1, uctx);
    fq_nmod_mpoly_from_mpolyu_perm_inflate(A, Abits, Au, perm, shift, stride, uctx, ctx);
    fq_nmod_mpolyu_clear(Au, uctx);
    return;
}


/*
    A is in              Fq [X][][v]
    E is in (Fq[v]/alpha(v))[X]
*/

void fq_nmod_mpolyun_redto_fq_nmod_poly(fq_nmod_poly_t E, const fq_nmod_mpoly_ctx_t ectx,
          fq_nmod_mpolyun_t A, const fq_nmod_mpoly_ctx_t ctx, const _fq_nmod_embed_t emb)
{
    fq_nmod_t v;
    slong Ai, Alen;
    fq_nmod_mpolyn_struct * Acoeff;
    ulong * Aexp;
/*
printf("fq_nmod_mpolyun_redto_fq_nmod_poly called\n");
printf("A: "); fq_nmod_mpolyun_print_pretty(A, NULL, ctx); printf("\n");
*/

    fq_nmod_init(v, ectx->fqctx);

    Acoeff = A->coeffs;
    Aexp = A->exps;
    Alen = A->length;

    fq_nmod_poly_zero(E, ectx->fqctx);
    for (Ai = 0; Ai < Alen; Ai++)
    {
        FLINT_ASSERT((Acoeff + Ai)->length == 1);
        _fq_nmod_embed_sm_to_lg(v, (Acoeff + Ai)->coeffs + 0, emb);
        fq_nmod_poly_set_coeff(E, Aexp[Ai], v, ectx->fqctx);
    }

    fq_nmod_clear(v, ectx->fqctx);
/*
printf("fq_nmod_mpolyun_redto_fq_nmod_poly returning\n");
printf("E: "); fq_nmod_poly_print_pretty(E, "X", ectx->fqctx); printf("\n");
*/
}


/*
    A = B
    A is in              Fq [X][][v]
    B is in (Fq[v]/alpha(v))[X]
*/
void fq_nmod_mpolyun_startinterp_fq_nmod_poly(slong * lastdeg_,
       fq_nmod_mpolyun_t A, const fq_nmod_mpoly_ctx_t ctx,
              fq_nmod_poly_t B, const fq_nmod_mpoly_ctx_t ectx, const _fq_nmod_embed_t emb)
{
    slong N = mpoly_words_per_exp_sp(A->bits, ctx->minfo);
    slong Bexp;
    slong Blen = fq_nmod_poly_length(B, ectx->fqctx);
    fq_nmod_struct * Bcoeff = B->coeffs;
    fq_nmod_mpolyn_struct * Acoeff;
    ulong * Aexp;
    slong Ai;
    slong lastdeg = -WORD(1);
/*
printf("fq_nmod_mpolyun_startinterp_fq_nmod_poly called\n");
printf("B: "); fq_nmod_poly_print_pretty(B, "X", ectx->fqctx); printf("\n");
*/

    fq_nmod_mpolyun_fit_length(A, Blen, ctx);
    Acoeff = A->coeffs;
    Aexp = A->exps;

    Ai = 0;
    for (Bexp = Blen - 1; Bexp >= 0; Bexp--)
    {
        if (!fq_nmod_is_zero(Bcoeff + Bexp, ectx->fqctx))
        {
            FLINT_ASSERT(Ai < A->alloc);

            fq_nmod_mpolyn_fit_length(Acoeff + Ai, 1, ctx);
            mpoly_monomial_zero((Acoeff + Ai)->exps + N*0, N);
            _fq_nmod_embed_lg_to_sm((Acoeff + Ai)->coeffs + 0, Bcoeff + Bexp, emb);
            FLINT_ASSERT(!fq_nmod_poly_is_zero((Acoeff + Ai)->coeffs + 0, ctx->fqctx));
            lastdeg = FLINT_MAX(lastdeg, fq_nmod_poly_degree((Acoeff + Ai)->coeffs + 0, ctx->fqctx));
            Aexp[Ai] = Bexp;
            (Acoeff + Ai)->length = 1;
            Ai++;
        }
    }
    A->length = Ai;
/*
printf("fq_nmod_mpolyun_startinterp_fq_nmod_poly returning\n");
printf("A: "); fq_nmod_mpolyun_print_pretty(A, NULL, ctx); printf("\n");
*/
    *lastdeg_ = lastdeg;
}


/*
    F = F + modulus*((A - F(alpha))/(modulus(alpha)))
    F is in              Fq [X][][v]
    A is in (Fq[v]/alpha(v))[X]
    alpha(v) is ectx->fqctx->modulus(v)
*/
int fq_nmod_mpolyun_addinterp_fq_nmod_poly(slong * lastdeg_,
  fq_nmod_mpolyun_t F, fq_nmod_mpolyun_t T, fq_nmod_poly_t modulus, const fq_nmod_mpoly_ctx_t ctx,
              fq_nmod_poly_t A, const fq_nmod_mpoly_ctx_t ectx, const _fq_nmod_embed_t emb)
{
    int changed = 0;
    slong lastdeg = -WORD(1);
    slong N = mpoly_words_per_exp_sp(F->bits, ctx->minfo);
    fq_nmod_t u, v;
    fq_nmod_poly_t u_sm, w;
    slong Fi, Toff, Aexp;
    fq_nmod_struct * Acoeff = A->coeffs;
    slong Flen = F->length;
    fq_nmod_mpolyn_struct * Fcoeff = F->coeffs;
    ulong * Fexp = F->exps;
    fq_nmod_mpolyn_struct * Tcoeff;
    ulong * Texp;
    fq_nmod_poly_t tp;
    fq_nmod_t inv_m_eval;

    fq_nmod_init(inv_m_eval, ectx->fqctx);
    _fq_nmod_embed_sm_to_lg(inv_m_eval, modulus, emb);
    fq_nmod_inv(inv_m_eval, inv_m_eval, ectx->fqctx);

    fq_nmod_init(u, ectx->fqctx);
    fq_nmod_init(v, ectx->fqctx);
    fq_nmod_poly_init(u_sm, ctx->fqctx);
    fq_nmod_poly_init(w, ctx->fqctx);

    Fi = 0;

/*
flint_printf("\nnmod_mpolyun_addinterp_bivar: (alpha = %wu)\n", alpha);
flint_printf("A: "); nmod_poly_print_pretty(A, "X"); printf("\n");
flint_printf("F: "); nmod_mpolyun_print_pretty(F, NULL, ctx); printf("\n");
*/

    Aexp = fq_nmod_poly_degree(A, ectx->fqctx);

    fq_nmod_poly_init(tp, ctx->fqctx);

    fq_nmod_mpolyun_fit_length(T, Flen + Aexp + 1, ctx);
    Tcoeff = T->coeffs;
    Texp = T->exps;
    Toff = 0;

    while (Fi < Flen || Aexp >= 0)
    {
        FLINT_ASSERT(Toff < T->alloc);
/*
flint_printf("Fi: %wd\n",Fi);
flint_printf("Aexp: %wd\n",Aexp);
*/
        if (Fi < Flen)
        {
            FLINT_ASSERT(!fq_nmod_poly_is_zero((Fcoeff + Fi)->coeffs + 0, ctx->fqctx));
            FLINT_ASSERT(fq_nmod_poly_degree((Fcoeff + Fi)->coeffs + 0, ctx->fqctx) < fq_nmod_poly_degree(modulus, ctx->fqctx));
        }

        if (Aexp >= 0)
        {
            FLINT_ASSERT(!fq_nmod_is_zero(Acoeff + Aexp, ectx->fqctx));
        }

        fq_nmod_mpolyn_fit_length(Tcoeff + Toff, 1, ctx);

        if (Fi < Flen && Aexp >= 0 && Fexp[Fi] == Aexp)
        {
            /* F term ok, A term ok */
            _fq_nmod_embed_sm_to_lg(u, (Fcoeff + Fi)->coeffs + 0, emb);
            fq_nmod_sub(v, Acoeff + Aexp, u, ectx->fqctx);
            if (!fq_nmod_is_zero(v, ectx->fqctx))
            {
                changed = 1;
                fq_nmod_mul(u, v, inv_m_eval, ectx->fqctx);
                _fq_nmod_embed_lg_to_sm(u_sm, u, emb);
                fq_nmod_poly_mul(w, modulus, u_sm, ctx->fqctx);
                fq_nmod_poly_add((Tcoeff + Toff)->coeffs + 0, (Fcoeff + Fi)->coeffs + 0, w, ctx->fqctx);
            }
            else
            {
                fq_nmod_poly_set((Tcoeff + Toff)->coeffs + 0, (Fcoeff + Fi)->coeffs + 0, ctx->fqctx);
            }
            Texp[Toff] = Aexp;
            Fi++;
            do {
                Aexp--;
            } while (Aexp >= 0 && fq_nmod_is_zero(Acoeff + Aexp, ectx->fqctx));
        }
        else if (Fi < Flen && (Aexp < 0 || Fexp[Fi] > Aexp))
        {
            /* F term ok, A term missing */
            _fq_nmod_embed_sm_to_lg(v, (Fcoeff + Fi)->coeffs + 0, emb);
            if (!fq_nmod_is_zero(v, ectx->fqctx))
            {
                changed = 1;
                fq_nmod_mul(u, v, inv_m_eval, ectx->fqctx);
                _fq_nmod_embed_lg_to_sm(u_sm, u, emb);
                fq_nmod_poly_mul(w, u_sm, modulus, ctx->fqctx);
                fq_nmod_poly_sub((Tcoeff + Toff)->coeffs + 0, (Fcoeff + Fi)->coeffs + 0, w, ctx->fqctx);
            }
            else
            {
                fq_nmod_poly_set((Tcoeff + Toff)->coeffs + 0, (Fcoeff + Fi)->coeffs + 0, ctx->fqctx);
            }
            Texp[Toff] = Fexp[Fi];
            Fi++;
        }
        else if (Aexp >= 0 && (Fi >= Flen || Fexp[Fi] < Aexp))
        {
            /* F term missing, A term ok */
            changed = 1;
            fq_nmod_mul(u, Acoeff + Aexp, inv_m_eval, ectx->fqctx);
            _fq_nmod_embed_lg_to_sm(u_sm, u, emb);
            fq_nmod_poly_mul((Tcoeff + Toff)->coeffs + 0, modulus, u_sm, ctx->fqctx);
            Texp[Toff] = Aexp;
            do {
                Aexp--;
            } while (Aexp >= 0 && fq_nmod_is_zero(Acoeff + Aexp, ectx->fqctx));
        }
        else
        {
            FLINT_ASSERT(0);
        }

        lastdeg = FLINT_MAX(lastdeg, fq_nmod_poly_degree((Tcoeff + Toff)->coeffs + 0, ctx->fqctx));
        mpoly_monomial_zero((Tcoeff + Toff)->exps + N*0, N);
        FLINT_ASSERT(!fq_nmod_poly_is_zero((Tcoeff + Toff)->coeffs + 0, ctx->fqctx));
        (Tcoeff + Toff)->length = 1;
        Toff++;
    }
    T->length = Toff;

    if (changed)
    {
        fq_nmod_mpolyun_swap(T, F);
    }

    fq_nmod_clear(u, ectx->fqctx);
    fq_nmod_clear(v, ectx->fqctx);
    fq_nmod_poly_clear(u_sm, ctx->fqctx);
    fq_nmod_poly_clear(w, ctx->fqctx);

    fq_nmod_clear(inv_m_eval, ectx->fqctx);

    *lastdeg_ = lastdeg;
    return changed;
}







int fq_nmod_mpolyun_gcd_brown_lgprime_bivar(fq_nmod_mpolyun_t G,
          fq_nmod_mpolyun_t Abar, fq_nmod_mpolyun_t Bbar, fq_nmod_mpolyun_t A,
                           fq_nmod_mpolyun_t B, const fq_nmod_mpoly_ctx_t ctx)
{
    int success;
    slong bound;
    fq_nmod_t temp, gammaeval;
    fq_nmod_poly_t Aeval, Beval, Geval, Abareval, Bbareval;
    fq_nmod_mpolyun_t T;
    slong deggamma, ldegG, ldegAbar, ldegBbar, ldegA, ldegB;
    fq_nmod_poly_t cA, cB, cG, cAbar, cBbar, gamma, trem;
    fq_nmod_poly_t modulus;
    flint_rand_t randstate;
    _fq_nmod_mpoly_embed_chooser_t embc;
    _fq_nmod_embed_struct * cur_emb;
    fq_nmod_mpoly_ctx_t ectx;
#if WANT_ASSERT
    fq_nmod_poly_t leadA, leadB;
#endif

#if WANT_ASSERT
    fq_nmod_poly_init(leadA, ctx->fqctx);
    fq_nmod_poly_init(leadB, ctx->fqctx);
    fq_nmod_poly_set(leadA, fq_nmod_mpolyun_leadcoeff_ref(A, ctx), ctx->fqctx);
    fq_nmod_poly_set(leadB, fq_nmod_mpolyun_leadcoeff_ref(B, ctx), ctx->fqctx);
#endif

    fq_nmod_poly_init(cA, ctx->fqctx);
    fq_nmod_poly_init(cB, ctx->fqctx);
    fq_nmod_mpolyun_content_last(cA, A, ctx);
    fq_nmod_mpolyun_content_last(cB, B, ctx);
    fq_nmod_mpolyun_divexact_last(A, cA, ctx);
    fq_nmod_mpolyun_divexact_last(B, cB, ctx);

    fq_nmod_poly_init(cG, ctx->fqctx);
    fq_nmod_poly_gcd(cG, cA, cB, ctx->fqctx);

    fq_nmod_poly_init(cAbar, ctx->fqctx);
    fq_nmod_poly_init(cBbar, ctx->fqctx);
    fq_nmod_poly_init(trem, ctx->fqctx);
    fq_nmod_poly_divrem(cAbar, trem, cA, cG, ctx->fqctx);
    FLINT_ASSERT(fq_nmod_poly_is_zero(trem, ctx->fqctx));
    fq_nmod_poly_divrem(cBbar, trem, cB, cG, ctx->fqctx);
    FLINT_ASSERT(fq_nmod_poly_is_zero(trem, ctx->fqctx));

    fq_nmod_poly_init(gamma, ctx->fqctx);
    fq_nmod_poly_gcd(gamma, fq_nmod_mpolyun_leadcoeff_ref(A, ctx),
                            fq_nmod_mpolyun_leadcoeff_ref(B, ctx), ctx->fqctx);

    ldegA = fq_nmod_mpolyun_lastdeg(A, ctx);
    ldegB = fq_nmod_mpolyun_lastdeg(B, ctx);
    deggamma = fq_nmod_poly_degree(gamma, ctx->fqctx);

    bound = 1 + deggamma + FLINT_MAX(ldegA, ldegB);

    fq_nmod_mpolyun_init(T, A->bits, ctx);

    fq_nmod_poly_init(modulus, ctx->fqctx);
    fq_nmod_poly_one(modulus, ctx->fqctx);

    flint_randinit(randstate);
    cur_emb = _fq_nmod_mpoly_embed_chooser_init(embc, ectx, ctx, randstate);

    /*
        Once Aeval, Beval, ..., t are inited in ectx->fqctx, they do not need
        to be cleared and reinited when ectx->fqctx changes.
    */
    fq_nmod_poly_init(Aeval, ectx->fqctx);
    fq_nmod_poly_init(Beval, ectx->fqctx);
    fq_nmod_poly_init(Geval, ectx->fqctx);
    fq_nmod_poly_init(Abareval, ectx->fqctx);
    fq_nmod_poly_init(Bbareval, ectx->fqctx);
    fq_nmod_init(gammaeval, ectx->fqctx);
    fq_nmod_init(temp, ectx->fqctx);

    /* initialization already picked a prime */
    goto have_prime;

choose_prime:

    cur_emb = _fq_nmod_mpoly_embed_chooser_next(embc, ectx, ctx, randstate);
    if (cur_emb == NULL)
    {
        /* ran out of primes */
        success = 0;
        goto cleanup;
    }

have_prime:

    /* make sure reduction does not kill both lc */
    _fq_nmod_embed_sm_to_lg(gammaeval, gamma, cur_emb);
    if (fq_nmod_is_zero(gammaeval, ectx->fqctx))
    {
        goto choose_prime;
    }

    /* make sure reduction does not kill either A or B */
    fq_nmod_mpolyun_redto_fq_nmod_poly(Aeval, ectx, A, ctx, cur_emb);
    fq_nmod_mpolyun_redto_fq_nmod_poly(Beval, ectx, B, ctx, cur_emb);
    if (Aeval->length == 0 || Beval->length == 0)
    {
        goto choose_prime;
    }

    fq_nmod_poly_gcd(Geval, Aeval, Beval, ectx->fqctx);
    success = fq_nmod_poly_divides(Abareval, Aeval, Geval, ectx->fqctx);
    FLINT_ASSERT(success);
    success = fq_nmod_poly_divides(Bbareval, Beval, Geval, ectx->fqctx);
    FLINT_ASSERT(success);

    FLINT_ASSERT(Geval->length > 0);
    FLINT_ASSERT(Abareval->length > 0);
    FLINT_ASSERT(Bbareval->length > 0);

    if (fq_nmod_poly_degree(Geval, ectx->fqctx) == 0)
    {
        fq_nmod_mpolyun_one(G, ctx);
        fq_nmod_mpolyun_swap(Abar, A);
        fq_nmod_mpolyun_swap(Bbar, B);
        goto successful_put_content;
    }

    if (fq_nmod_poly_degree(modulus, ctx->fqctx) > 0)
    {
        FLINT_ASSERT(G->length > 0);
        if (fq_nmod_poly_degree(Geval, ectx->fqctx) > G->exps[0])
        {
            goto choose_prime;
        }
        else if (fq_nmod_poly_degree(Geval, ectx->fqctx) < G->exps[0])
        {
            fq_nmod_poly_one(modulus, ctx->fqctx);
        }
    }

    fq_nmod_poly_scalar_mul_fq_nmod(Geval, Geval, gammaeval, ectx->fqctx);

    if (fq_nmod_poly_degree(modulus, ctx->fqctx) > 0)
    {
        fq_nmod_mpolyun_addinterp_fq_nmod_poly(&ldegG, G, T, modulus, ctx, Geval, ectx, cur_emb);
        fq_nmod_mpolyun_addinterp_fq_nmod_poly(&ldegAbar, Abar, T, modulus, ctx, Abareval, ectx, cur_emb);
        fq_nmod_mpolyun_addinterp_fq_nmod_poly(&ldegBbar, Bbar, T, modulus, ctx, Bbareval, ectx, cur_emb);
    }
    else
    {
        fq_nmod_mpolyun_startinterp_fq_nmod_poly(&ldegG, G, ctx, Geval, ectx, cur_emb);
        fq_nmod_mpolyun_startinterp_fq_nmod_poly(&ldegAbar, Abar, ctx, Abareval, ectx, cur_emb);
        fq_nmod_mpolyun_startinterp_fq_nmod_poly(&ldegBbar, Bbar, ctx, Bbareval, ectx, cur_emb);
    }
    fq_nmod_poly_mul(modulus, modulus, cur_emb->h, ctx->fqctx);

    if (fq_nmod_poly_degree(modulus, ctx->fqctx) < bound)
    {
        goto choose_prime;
    }

    FLINT_ASSERT(ldegG >= 0);
    FLINT_ASSERT(ldegAbar >= 0);
    FLINT_ASSERT(ldegBbar >= 0);

    if (   deggamma + ldegA == ldegG + ldegAbar
        && deggamma + ldegB == ldegG + ldegBbar )
    {
        goto successful;
    }

    fq_nmod_poly_one(modulus, ctx->fqctx);
    goto choose_prime;

successful:

    fq_nmod_mpolyun_content_last(modulus, G, ctx);
    fq_nmod_mpolyun_divexact_last(G, modulus, ctx);
    fq_nmod_mpolyun_divexact_last(Abar, fq_nmod_mpolyun_leadcoeff_ref(G, ctx), ctx);
    fq_nmod_mpolyun_divexact_last(Bbar, fq_nmod_mpolyun_leadcoeff_ref(G, ctx), ctx);

successful_put_content:

    fq_nmod_mpolyun_mul_last(G, cG, ctx);
    fq_nmod_mpolyun_mul_last(Abar, cAbar, ctx);
    fq_nmod_mpolyun_mul_last(Bbar, cBbar, ctx);

    success = 1;

cleanup:

#if WANT_ASSERT
    if (success)
    {
        fq_nmod_poly_mul(modulus, fq_nmod_mpolyun_leadcoeff_ref(G, ctx),
                                  fq_nmod_mpolyun_leadcoeff_ref(Abar, ctx), ctx->fqctx);
        FLINT_ASSERT(fq_nmod_poly_equal(modulus, leadA, ctx->fqctx));
        fq_nmod_poly_mul(modulus, fq_nmod_mpolyun_leadcoeff_ref(G, ctx),
                                  fq_nmod_mpolyun_leadcoeff_ref(Bbar, ctx), ctx->fqctx);
        FLINT_ASSERT(fq_nmod_poly_equal(modulus, leadB, ctx->fqctx));
    }
    fq_nmod_poly_clear(leadA, ctx->fqctx);
    fq_nmod_poly_clear(leadB, ctx->fqctx);
#endif

    fq_nmod_poly_clear(cA, ctx->fqctx);
    fq_nmod_poly_clear(cB, ctx->fqctx);
    fq_nmod_poly_clear(cG, ctx->fqctx);
    fq_nmod_poly_clear(cAbar, ctx->fqctx);
    fq_nmod_poly_clear(cBbar, ctx->fqctx);
    fq_nmod_poly_clear(trem, ctx->fqctx);
    fq_nmod_poly_clear(gamma, ctx->fqctx);

    fq_nmod_mpolyun_clear(T, ctx);

    fq_nmod_poly_clear(modulus, ctx->fqctx);

    fq_nmod_poly_clear(Aeval, ectx->fqctx);
    fq_nmod_poly_clear(Beval, ectx->fqctx);
    fq_nmod_poly_clear(Geval, ectx->fqctx);
    fq_nmod_poly_clear(Abareval, ectx->fqctx);
    fq_nmod_poly_clear(Bbareval, ectx->fqctx);
    fq_nmod_clear(gammaeval, ectx->fqctx);
    fq_nmod_clear(temp, ectx->fqctx);

    _fq_nmod_mpoly_embed_chooser_clear(embc, ectx, ctx, randstate);

    flint_randclear(randstate);

    return success;
}




/*
    A is in                      Fq [x_0,...,x_(var-2), x_(var-1)][x_var]
    E is in (Fq[x_var]/alpha(x_var))[x_0,...,x_(var-2)][x_(var-1)]
*/

void fq_nmod_mpolyn_redto_fq_nmod_mpolyn(fq_nmod_mpolyn_t E, const fq_nmod_mpoly_ctx_t ectx,
 fq_nmod_mpolyn_t A, slong var, const fq_nmod_mpoly_ctx_t ctx, const _fq_nmod_embed_t emb)
{
    slong N = mpoly_words_per_exp_sp(A->bits, ctx->minfo);
    slong offset, shift, k;
    ulong mask;
    fq_nmod_t v;
    fq_nmod_poly_struct * Acoeff = A->coeffs;
    ulong * Aexp = A->exps;
    slong Alen = A->length;
    slong Ai;
    fq_nmod_poly_struct * Ecoeff;
    ulong * Eexp;
    slong Ei;

    fq_nmod_init(v, ectx->fqctx);

    FLINT_ASSERT(var > 0);
    mpoly_gen_offset_shift_sp(&offset, &shift, var - 1, A->bits, ctx->minfo);
    mask = (-UWORD(1)) >> (FLINT_BITS - A->bits);

    Ecoeff = E->coeffs;
    Eexp = E->exps;
    Ei = 0;
    for (Ai = 0; Ai < Alen; Ai++)
    {
        _fq_nmod_embed_sm_to_lg(v, Acoeff + Ai, emb);
        k = ((Aexp + N*Ai)[offset] >> shift) & mask;
        if (fq_nmod_is_zero(v, ectx->fqctx))
        {
            continue;
        }

        if (Ei > 0 && mpoly_monomial_equal_extra(Eexp + N*(Ei - 1), Aexp + N*Ai, N, offset, -(k << shift)))
        {
            /* append to previous */
            fq_nmod_poly_set_coeff(Ecoeff + Ei - 1, k, v, ectx->fqctx);
        }
        else
        {
            FLINT_ASSERT(Ei == 0 || mpoly_monomial_gt_nomask_extra(Eexp + N*(Ei - 1), Aexp + N*Ai, N, offset, -(k << shift)));

            /* create new */
            if (Ei >= E->alloc)
            {
                fq_nmod_mpolyn_fit_length(E, Ei + 1, ectx);
                Ecoeff = E->coeffs;
                Eexp = E->exps;
            }
            mpoly_monomial_set_extra(Eexp + N*Ei, Aexp + N*Ai, N, offset, -(k << shift));
            fq_nmod_poly_zero(Ecoeff + Ei, ectx->fqctx);
            fq_nmod_poly_set_coeff(Ecoeff + Ei, k, v, ectx->fqctx);
            Ei++;
        }
    }
    E->length = Ei;

    fq_nmod_clear(v, ectx->fqctx);
}

/*
    A is in                      Fq [X][x_0,...,x_(var-2), x_(var-1)][x_var]
    E is in (Fq[x_var]/alpha(x_var))[X][x_0,...,x_(var-2)][x_(var-1)]
*/
void fq_nmod_mpolyun_redto_fq_nmod_mpolyun(fq_nmod_mpolyun_t E, const fq_nmod_mpoly_ctx_t ectx,
 fq_nmod_mpolyun_t A, slong var, const fq_nmod_mpoly_ctx_t ctx, const _fq_nmod_embed_t emb)
{
    fq_nmod_mpolyn_struct * Acoeff = A->coeffs;
    ulong * Aexp = A->exps;
    slong Alen = A->length;
    slong Ai;
    fq_nmod_mpolyn_struct * Ecoeff;
    ulong * Eexp;
    slong Ei;

/*
printf("fq_nmod_mpolyun_redto_fq_nmod_mpolyun\n");
printf("A: "); fq_nmod_mpolyun_print_pretty(A, NULL, ctx); printf("\n");
*/

    fq_nmod_mpolyun_fit_length(E, Alen, ectx);
    Ecoeff = E->coeffs;
    Eexp = E->exps;

    Ei = 0;
    for (Ai = 0; Ai < Alen; Ai++)
    {
        fq_nmod_mpolyn_redto_fq_nmod_mpolyn(Ecoeff + Ei, ectx, Acoeff + Ai, var, ctx, emb);
        Eexp[Ei] = Aexp[Ai];
        Ei += ((Ecoeff + Ei)->length != 0);
    }
    E->length = Ei;

    FLINT_ASSERT(fq_nmod_mpolyun_is_canonical(E, ectx));
/*
printf("fq_nmod_mpolyun_redto_fq_nmod_mpolyun returning\n");
printf("E: "); fq_nmod_mpolyun_print_pretty(E, NULL, ectx); printf("\n");
*/
}


/*
    A = B using lowest degree representative
    A is in                      Fq [x_0,...,x_(var-2), x_(var-1)][x_var]
    B is in (Fq[x_var]/alpha(x_var))[x_0,...,x_(var-2)][x_(var-1)]
    alpha is emb->h
*/
void fq_nmod_mpolyn_startinterp_n_lgprime(slong * lastdeg_,
       fq_nmod_mpolyn_t A, slong var, const fq_nmod_mpoly_ctx_t ctx,
              fq_nmod_mpolyn_t B, const fq_nmod_mpoly_ctx_t ectx, const _fq_nmod_embed_t emb)
{
    slong N = mpoly_words_per_exp_sp(B->bits, ctx->minfo);
    slong offset, shift;
    slong vi;
    fq_nmod_poly_struct * Bcoeff = B->coeffs;
    ulong * Bexp = B->exps;
    slong Blen = B->length;
    slong Bi;
    fq_nmod_poly_struct * Acoeff;
    ulong * Aexp;
    slong Ai;
    slong lastdeg = -WORD(1);

    FLINT_ASSERT(A->bits == B->bits);

    fq_nmod_mpolyn_fit_length(A, Blen, ctx);
    Acoeff = A->coeffs;
    Aexp = A->exps;

    mpoly_gen_offset_shift_sp(&offset, &shift, var - 1, A->bits, ctx->minfo);

    Ai = 0;
    for (Bi = 0; Bi < Blen; Bi++)
    {
        if (Ai + (Bcoeff + Bi)->length >= A->alloc)
        {
            fq_nmod_mpolyn_fit_length(A, Ai + (Bcoeff + Bi)->length, ctx);
            Acoeff = A->coeffs;
            Aexp = A->exps;
        }
        for (vi = (Bcoeff + Bi)->length - 1; vi >= 0; vi--)
        {
            if (!fq_nmod_is_zero((Bcoeff + Bi)->coeffs + vi, ectx->fqctx))
            {
                mpoly_monomial_set_extra(Aexp + N*Ai, Bexp + N*Bi, N, offset, vi << shift);
                _fq_nmod_embed_lg_to_sm(Acoeff + Ai, (Bcoeff + Bi)->coeffs + vi, emb);
                FLINT_ASSERT(!fq_nmod_poly_is_zero(Acoeff + Ai, ctx->fqctx));
                lastdeg = FLINT_MAX(lastdeg, fq_nmod_poly_degree(Acoeff + Ai, ctx->fqctx));
                Ai++;
            }
        }
    }
    A->length = Ai;

    *lastdeg_ = lastdeg;
}

/*
    A = B using lowest degree representative
    A is in                      Fq [X][x_0,...,x_(var-2), x_(var-1)][x_var]
    B is in (Fq[x_var]/alpha(x_var))[X][x_0,...,x_(var-2)][x_(var-1)]
    alpha is emb->h
*/
void fq_nmod_mpolyun_startinterp_un_lgprime(slong * lastdeg,
       fq_nmod_mpolyun_t A, slong var, const fq_nmod_mpoly_ctx_t ctx,
              fq_nmod_mpolyun_t B, const fq_nmod_mpoly_ctx_t ectx, const _fq_nmod_embed_t emb)
{
    slong i;
    fq_nmod_mpolyn_struct * Bcoeff = B->coeffs;
    ulong * Bexp = B->exps;
    slong Blen = B->length;
    fq_nmod_mpolyn_struct * Acoeff;
    ulong * Aexp;

    FLINT_ASSERT(var > 0);
/*
flint_printf("fq_nmod_mpolyun_startinterp_un_lgprime called (var = %wd)\n", var);
printf("B: "); fq_nmod_mpolyun_print_pretty(B, NULL, ectx); printf("\n");
*/
    *lastdeg = -WORD(1);

    fq_nmod_mpolyun_fit_length(A, Blen, ctx);
    Acoeff = A->coeffs;
    Aexp = A->exps;

    for (i = 0; i < Blen; i++)
    {
        Aexp[i] = Bexp[i];
        fq_nmod_mpolyn_startinterp_n_lgprime(lastdeg, Acoeff + i, var, ctx, Bcoeff + i, ectx, emb);
    }
    A->length = Blen;

    FLINT_ASSERT(fq_nmod_mpolyun_is_canonical(A, ctx));
/*
printf("fq_nmod_mpolyun_startinterp_un_lgprime returning\n");
printf("A: "); fq_nmod_mpolyun_print_pretty(A, NULL, ectx); printf("\n");
*/
}


/*
    T = F + modulus*((A - F(alpha))/modulus(alpha))
    F is in                      Fq [x_0,...,x_(var-2), x_(var-1)][x_var]
    A is in (Fq[x_var]/alpha(x_var))[x_0,...,x_(var-2)][x_(var-1)]
*/
int fq_nmod_mpolyn_addinterp_n_lgprime(slong * lastdeg_,
  fq_nmod_mpolyn_t T, fq_nmod_mpolyn_t F, fq_nmod_poly_t modulus, slong var, const fq_nmod_mpoly_ctx_t ctx,
              fq_nmod_mpolyn_t A, const fq_nmod_mpoly_ctx_t ectx, const _fq_nmod_embed_t emb)
{
    int changed = 0;
    slong N = mpoly_words_per_exp_sp(A->bits, ctx->minfo);
    slong offset, shift;
    slong lastdeg = -WORD(1);
    slong vi;
    fq_nmod_t u, v;
    fq_nmod_poly_t u_sm, w;
    fq_nmod_t inv_m_eval;

    fq_nmod_poly_struct * Tcoeff;
    ulong * Texp;
    slong Ti;

    fq_nmod_poly_struct * Acoeff = A->coeffs;
    slong Alen = A->length;
    ulong * Aexp = A->exps;
    slong Ai;

    fq_nmod_poly_struct * Fcoeff = F->coeffs;
    slong Flen = F->length;
    ulong * Fexp = F->exps;
    slong Fi;

    fq_nmod_init(inv_m_eval, ectx->fqctx);
    _fq_nmod_embed_sm_to_lg(inv_m_eval, modulus, emb);
    fq_nmod_inv(inv_m_eval, inv_m_eval, ectx->fqctx);

    fq_nmod_init(u, ectx->fqctx);
    fq_nmod_init(v, ectx->fqctx);
    fq_nmod_poly_init(u_sm, ctx->fqctx);
    fq_nmod_poly_init(w, ctx->fqctx);

    FLINT_ASSERT(var > 0);
    FLINT_ASSERT(T->bits == A->bits);
    FLINT_ASSERT(F->bits == A->bits);

    mpoly_gen_offset_shift_sp(&offset, &shift, var - 1, A->bits, ctx->minfo);

    Flen = F->length;

    fq_nmod_mpolyn_fit_length(T, FLINT_MAX(Flen, Alen), ctx);
    Tcoeff = T->coeffs;
    Texp = T->exps;
    Ti = 0;

    Fi = Ai = vi = 0;
    if (Ai < Alen)
    {
        vi = fq_nmod_poly_degree(A->coeffs + Ai, ectx->fqctx);
    }
    while (Fi < Flen || Ai < Alen)
    {
        if (Ti >= T->alloc)
        {
            fq_nmod_mpolyn_fit_length(T, Ti + FLINT_MAX(Flen - Fi, Alen - Ai), ctx);
            Tcoeff = T->coeffs;
            Texp = T->exps;
        }

        if (Fi < Flen && Ai < Alen && mpoly_monomial_equal_extra(Fexp + N*Fi, Aexp + N*Ai, N, offset, vi << shift))
        {
            /* F term ok, A term ok */
            _fq_nmod_embed_sm_to_lg(u, Fcoeff + Fi, emb);
            fq_nmod_sub(v, (Acoeff + Ai)->coeffs + vi, u, ectx->fqctx);
            if (!fq_nmod_is_zero(v, ectx->fqctx))
            {
                changed = 1;
                fq_nmod_mul(u, v, inv_m_eval, ectx->fqctx);
                _fq_nmod_embed_lg_to_sm(u_sm, u, emb);
                fq_nmod_poly_mul(w, modulus, u_sm, ctx->fqctx);
                fq_nmod_poly_add(Tcoeff + Ti, Fcoeff + Fi, w, ctx->fqctx);
            }
            else
            {
                fq_nmod_poly_set(Tcoeff + Ti, Fcoeff + Fi, ctx->fqctx);
            }
            mpoly_monomial_set(Texp + N*Ti, Fexp + N*Fi, N);

            Fi++;
            do {
                vi--;
            } while (vi >= 0 && fq_nmod_is_zero((Acoeff + Ai)->coeffs + vi, ectx->fqctx));
            if (vi < 0)
            {
                Ai++;
                if (Ai < Alen)
                {
                    vi = fq_nmod_poly_degree(A->coeffs + Ai, ectx->fqctx);
                }
            }
        }
        else if (Fi < Flen && (Ai >= Alen || mpoly_monomial_gt_nomask_extra(Fexp + N*Fi, Aexp + N*Ai, N, offset, vi << shift)))
        {
            /* F term ok, A term missing */
            _fq_nmod_embed_sm_to_lg(v, Fcoeff + Fi, emb);
            if (!fq_nmod_is_zero(v, ectx->fqctx))
            {
                changed = 1;
                fq_nmod_mul(u, v, inv_m_eval, ectx->fqctx);
                _fq_nmod_embed_lg_to_sm(u_sm, u, emb);
                fq_nmod_poly_mul(w, u_sm, modulus, ctx->fqctx);
                fq_nmod_poly_sub(Tcoeff + Ti, Fcoeff + Fi, w, ctx->fqctx);
            }
            else
            {
                fq_nmod_poly_set(Tcoeff + Ti, Fcoeff + Fi, ctx->fqctx);
            }
            mpoly_monomial_set(Texp + N*Ti, Fexp + N*Fi, N);

            Fi++;
        }
        else
        {
            FLINT_ASSERT(Ai < Alen && (Fi >= Flen || mpoly_monomial_lt_nomask_extra(Fexp + N*Fi, Aexp + N*Ai, N, offset, vi << shift)));

            /* F term missing, A term ok */
            changed = 1;
            fq_nmod_mul(u, (Acoeff + Ai)->coeffs + vi, inv_m_eval, ectx->fqctx);
            _fq_nmod_embed_lg_to_sm(u_sm, u, emb);
            fq_nmod_poly_mul(Tcoeff + Ti, modulus, u_sm, ctx->fqctx);
            mpoly_monomial_set_extra(Texp + N*Ti, Aexp + N*Ai, N, offset, vi << shift);

            do {
                vi--;
            } while (vi >= 0 && fq_nmod_is_zero((Acoeff + Ai)->coeffs + vi, ectx->fqctx));
            if (vi < 0)
            {
                Ai++;
                if (Ai < Alen)
                {
                    vi = fq_nmod_poly_degree(A->coeffs + Ai, ectx->fqctx);
                }
            }
        }

        FLINT_ASSERT(!fq_nmod_poly_is_zero(Tcoeff + Ti, ctx->fqctx));
        lastdeg = FLINT_MAX(lastdeg, fq_nmod_poly_degree(Tcoeff + Ti, ctx->fqctx));
        Ti++;
    }
    T->length = Ti;


    fq_nmod_clear(inv_m_eval, ectx->fqctx);
    fq_nmod_clear(u, ectx->fqctx);
    fq_nmod_clear(v, ectx->fqctx);
    fq_nmod_poly_clear(u_sm, ctx->fqctx);
    fq_nmod_poly_clear(w, ctx->fqctx);

    *lastdeg_ = FLINT_MAX(*lastdeg_, lastdeg);
    return changed;
}

/*
    F = F + modulus*((A - F(alpha))/modulus(alpha))
    F is in                      Fq [X][x_0,...,x_(var-2), x_(var-1)][x_var]
    A is in (Fq[x_var]/alpha(x_var))[X][x_0,...,x_(var-2)][x_(var-1)]
*/
int fq_nmod_mpolyun_addinterp_un_lgprime(slong * lastdeg,
  fq_nmod_mpolyun_t F, fq_nmod_mpolyun_t T, fq_nmod_poly_t modulus, slong var, const fq_nmod_mpoly_ctx_t ctx,
              fq_nmod_mpolyun_t A, const fq_nmod_mpoly_ctx_t ectx, const _fq_nmod_embed_t emb)
{
    int changed = 0;
    slong i, j, k;
    ulong * Texp;
    ulong * Fexp;
    ulong * Aexp;
    slong Flen;
    slong Alen;
    fq_nmod_mpolyn_struct * Tcoeff;
    fq_nmod_mpolyn_struct * Fcoeff;
    fq_nmod_mpolyn_t zero;
    fq_nmod_mpolyn_struct  * Acoeff;
    fq_nmod_mpolyn_t ezero;

/*
flint_printf("fq_nmod_mpolyun_addinterp_un_lgprime called (var = %wd)\n", var);
printf("F: "); fq_nmod_mpolyun_print_pretty(F, NULL, ctx); printf("\n");
printf("A: "); fq_nmod_mpolyun_print_pretty(A, NULL, ectx); printf("\n");
*/

    FLINT_ASSERT(var > 0);

    *lastdeg = -WORD(1);

    FLINT_ASSERT(F->bits == T->bits);
    FLINT_ASSERT(T->bits == A->bits);

    Flen = F->length;
    Alen = A->length;
    fq_nmod_mpolyun_fit_length(T, Flen + Alen, ctx);

    Tcoeff = T->coeffs;
    Fcoeff = F->coeffs;
    Acoeff = A->coeffs;
    Texp = T->exps;
    Fexp = F->exps;
    Aexp = A->exps;   

    fq_nmod_mpolyn_init(zero, A->bits, ctx);
    zero->length = 0;
    fq_nmod_mpolyn_init(ezero, A->bits, ectx);
    ezero->length = 0;

    i = j = k = 0;
    while (i < Flen || j < Alen)
    {
        if (i < Flen && j < Alen && (Fexp[i] == Aexp[j]))
        {
            /* F term ok, A term ok */
            changed |= fq_nmod_mpolyn_addinterp_n_lgprime(lastdeg, Tcoeff + k,
                         Fcoeff + i, modulus, var, ctx, Acoeff + j, ectx, emb);
            Texp[k] = Aexp[j];
            i++;
            j++;
        }
        else if (i < Flen && (j >= Alen || Fexp[i] > Aexp[j]))
        {
            /* F term ok, A term missing */
            changed |= fq_nmod_mpolyn_addinterp_n_lgprime(lastdeg, Tcoeff + k,
                              Fcoeff + i, modulus, var, ctx, ezero, ectx, emb);
            Texp[k] = Fexp[i];
            i++;
        }
        else
        {
            FLINT_ASSERT(j < Alen && (i >= Flen || Aexp[j] > Fexp[i]));

            /* F term missing, A term ok */
            changed |= fq_nmod_mpolyn_addinterp_n_lgprime(lastdeg, Tcoeff + k,
                               zero, modulus, var, ctx, Acoeff + j, ectx, emb);
            Texp[k] = Aexp[j];
            j++;
        }

        FLINT_ASSERT(!fq_nmod_mpolyn_is_zero(Tcoeff + k, ctx));
        k++;
    }
    T->length = k;

    if (changed)
    {
        fq_nmod_mpolyun_swap(T, F);
    }

    fq_nmod_mpolyn_clear(zero, ctx);
    fq_nmod_mpolyn_clear(ezero, ectx);

    FLINT_ASSERT(fq_nmod_mpolyun_is_canonical(F, ctx));
/*
printf("fq_nmod_mpolyun_addinterp_un_lgprime returning %d\n", changed);
printf("F: "); fq_nmod_mpolyun_print_pretty(F, NULL, ctx); printf("\n");
*/

    return changed;
}



int fq_nmod_mpolyun_gcd_brown_lgprime(fq_nmod_mpolyun_t G,
          fq_nmod_mpolyun_t Abar, fq_nmod_mpolyun_t Bbar, fq_nmod_mpolyun_t A,
                  fq_nmod_mpolyun_t B, slong var, const fq_nmod_mpoly_ctx_t ctx)
{
    int success;
    slong bound;
    slong offset, shift;
    fq_nmod_t temp, gammaeval;
    fq_nmod_mpolyun_t Aeval, Beval, Geval, Abareval, Bbareval;
    fq_nmod_mpolyun_t T;
    slong deggamma, ldegG, ldegAbar, ldegBbar, ldegA, ldegB;
    fq_nmod_poly_t cA, cB, cG, cAbar, cBbar, gamma, trem;
    fq_nmod_poly_t modulus;
    mp_bitcnt_t bits = A->bits;
    slong N = mpoly_words_per_exp_sp(bits, ctx->minfo);
    flint_rand_t randstate;
    _fq_nmod_mpoly_embed_chooser_t embc;
    _fq_nmod_embed_struct * cur_emb;
    fq_nmod_mpoly_ctx_t ectx;
#if WANT_ASSERT
    fq_nmod_poly_t leadA, leadB;
#endif

    FLINT_ASSERT(var >= 0);
    if (var == WORD(0))
    {
        return fq_nmod_mpolyun_gcd_brown_lgprime_bivar(G, Abar, Bbar, A, B, ctx);
    }

    mpoly_gen_offset_shift_sp(&offset, &shift, var - 1, bits, ctx->minfo);

#if WANT_ASSERT
    fq_nmod_poly_init(leadA, ctx->fqctx);
    fq_nmod_poly_init(leadB, ctx->fqctx);
    fq_nmod_poly_set(leadA, fq_nmod_mpolyun_leadcoeff_ref(A, ctx), ctx->fqctx);
    fq_nmod_poly_set(leadB, fq_nmod_mpolyun_leadcoeff_ref(B, ctx), ctx->fqctx);
#endif

    fq_nmod_poly_init(cA, ctx->fqctx);
    fq_nmod_poly_init(cB, ctx->fqctx);
    fq_nmod_mpolyun_content_last(cA, A, ctx);
    fq_nmod_mpolyun_content_last(cB, B, ctx);
    fq_nmod_mpolyun_divexact_last(A, cA, ctx);
    fq_nmod_mpolyun_divexact_last(B, cB, ctx);

    fq_nmod_poly_init(cG, ctx->fqctx);
    fq_nmod_poly_gcd(cG, cA, cB, ctx->fqctx);

    fq_nmod_poly_init(cAbar, ctx->fqctx);
    fq_nmod_poly_init(cBbar, ctx->fqctx);
    fq_nmod_poly_init(trem, ctx->fqctx);
    fq_nmod_poly_divrem(cAbar, trem, cA, cG, ctx->fqctx);
    FLINT_ASSERT(fq_nmod_poly_is_zero(trem, ctx->fqctx));
    fq_nmod_poly_divrem(cBbar, trem, cB, cG, ctx->fqctx);
    FLINT_ASSERT(fq_nmod_poly_is_zero(trem, ctx->fqctx));

    fq_nmod_poly_init(gamma, ctx->fqctx);
    fq_nmod_poly_gcd(gamma, fq_nmod_mpolyun_leadcoeff_ref(A, ctx),
                            fq_nmod_mpolyun_leadcoeff_ref(B, ctx), ctx->fqctx);

    ldegA = fq_nmod_mpolyun_lastdeg(A, ctx);
    ldegB = fq_nmod_mpolyun_lastdeg(B, ctx);
    deggamma = fq_nmod_poly_degree(gamma, ctx->fqctx);

    bound = 1 + deggamma + FLINT_MAX(ldegA, ldegB);

    fq_nmod_mpolyun_init(T, bits, ctx);

    fq_nmod_poly_init(modulus, ctx->fqctx);
    fq_nmod_poly_one(modulus, ctx->fqctx);

    flint_randinit(randstate);
    cur_emb = _fq_nmod_mpoly_embed_chooser_init(embc, ectx, ctx, randstate);

    /*
        Once Aeval, Beval, ..., t are inited in ectx->fqctx, they do not need
        to be cleared and reinited when ectx->fqctx changes.
    */
    fq_nmod_mpolyun_init(Aeval, bits, ectx);
    fq_nmod_mpolyun_init(Beval, bits, ectx);
    fq_nmod_mpolyun_init(Geval, bits, ectx);
    fq_nmod_mpolyun_init(Abareval, bits, ectx);
    fq_nmod_mpolyun_init(Bbareval, bits, ectx);
    fq_nmod_init(gammaeval, ectx->fqctx);
    fq_nmod_init(temp, ectx->fqctx);

    /* initialization already picked a prime */
    goto have_prime;

choose_prime: /* prime is irreducible element of Fq[v] (cur_emb->h) */

    cur_emb = _fq_nmod_mpoly_embed_chooser_next(embc, ectx, ctx, randstate);
    if (cur_emb == NULL)
    {
        /* ran out of primes */
        success = 0;
        goto cleanup;
    }

have_prime:
    /* make sure reduction does not kill both lc */
    _fq_nmod_embed_sm_to_lg(gammaeval, gamma, cur_emb);
    if (fq_nmod_is_zero(gammaeval, ectx->fqctx))
    {
        goto choose_prime;
    }

    /* make sure reduction does not kill either A or B */
    fq_nmod_mpolyun_redto_fq_nmod_mpolyun(Aeval, ectx, A, var, ctx, cur_emb);
    fq_nmod_mpolyun_redto_fq_nmod_mpolyun(Beval, ectx, B, var, ctx, cur_emb);
    if (Aeval->length == 0 || Beval->length == 0)
    {
        goto choose_prime;
    }

    FLINT_ASSERT(fq_nmod_mpolyun_is_canonical(Aeval, ectx));
    FLINT_ASSERT(fq_nmod_mpolyun_is_canonical(Beval, ectx));

    success = fq_nmod_mpolyun_gcd_brown_smprime(Geval, Abareval, Bbareval,
                                                  Aeval, Beval, var - 1, ectx);
    if (success == 0)
    {
        goto choose_prime;
    }

    if (fq_nmod_mpolyun_is_nonzero_fq_nmod(Geval, ectx))
    {
        fq_nmod_mpolyun_one(G, ctx);
        fq_nmod_mpolyun_swap(Abar, A);
        fq_nmod_mpolyun_swap(Bbar, B);
        goto successful_put_content;
    }

    if (fq_nmod_poly_degree(modulus, ctx->fqctx) > 0)
    {
        /* compare leading monomials of Geval and G */
        int cmp = 0;
        FLINT_ASSERT(G->length > 0);
        if (G->exps[0] != Geval->exps[0])
        {
            cmp = G->exps[0] > Geval->exps[0] ? 1 : -1;
        }
        if (cmp == 0)
        {
            slong k = fq_nmod_poly_degree((Geval->coeffs + 0)->coeffs + 0, ectx->fqctx);
            cmp = mpoly_monomial_cmp_nomask_extra(
                        (G->coeffs + 0)->exps + N*0,
                    (Geval->coeffs + 0)->exps + N*0, N, offset, k << shift);
        }

        if (cmp < 0)
        {
            goto choose_prime;
        }
        else if (cmp > 0)
        {
            fq_nmod_poly_one(modulus, ctx->fqctx);
        }
    }

    fq_nmod_inv(temp, fq_nmod_mpolyn_leadcoeff_last_ref(Geval->coeffs + 0, ectx), ectx->fqctx);
    fq_nmod_mul(temp, temp, gammaeval, ectx->fqctx);
    fq_nmod_mpolyun_scalar_mul_fq_nmod(Geval, temp, ectx);

    if (fq_nmod_poly_degree(modulus, ctx->fqctx) > 0)
    {
        fq_nmod_mpolyun_addinterp_un_lgprime(&ldegG, G, T, modulus, var, ctx, Geval, ectx, cur_emb);
        fq_nmod_mpolyun_addinterp_un_lgprime(&ldegAbar, Abar, T, modulus, var, ctx, Abareval, ectx, cur_emb);
        fq_nmod_mpolyun_addinterp_un_lgprime(&ldegBbar, Bbar, T, modulus, var, ctx, Bbareval, ectx, cur_emb);
    }
    else
    {
        fq_nmod_mpolyun_startinterp_un_lgprime(&ldegG, G, var, ctx, Geval, ectx, cur_emb);
        fq_nmod_mpolyun_startinterp_un_lgprime(&ldegAbar, Abar, var, ctx, Abareval, ectx, cur_emb);
        fq_nmod_mpolyun_startinterp_un_lgprime(&ldegBbar, Bbar, var, ctx, Bbareval, ectx, cur_emb);
    }
    fq_nmod_poly_mul(modulus, modulus, cur_emb->h, ctx->fqctx);

    if (fq_nmod_poly_degree(modulus, ctx->fqctx) < bound)
    {
        goto choose_prime;
    }

    FLINT_ASSERT(ldegG >= 0);
    FLINT_ASSERT(ldegAbar >= 0);
    FLINT_ASSERT(ldegBbar >= 0);

    if (   deggamma + ldegA == ldegG + ldegAbar
        && deggamma + ldegB == ldegG + ldegBbar )
    {
        goto successful;
    }

    fq_nmod_poly_one(modulus, ctx->fqctx);
    goto choose_prime;

successful:

    fq_nmod_mpolyun_content_last(modulus, G, ctx);
    fq_nmod_mpolyun_divexact_last(G, modulus, ctx);
    fq_nmod_mpolyun_divexact_last(Abar, fq_nmod_mpolyun_leadcoeff_ref(G, ctx), ctx);
    fq_nmod_mpolyun_divexact_last(Bbar, fq_nmod_mpolyun_leadcoeff_ref(G, ctx), ctx);

successful_put_content:

    fq_nmod_mpolyun_mul_last(G, cG, ctx);
    fq_nmod_mpolyun_mul_last(Abar, cAbar, ctx);
    fq_nmod_mpolyun_mul_last(Bbar, cBbar, ctx);

    success = 1;

cleanup:

#if WANT_ASSERT
    if (success)
    {
        fq_nmod_poly_mul(modulus, fq_nmod_mpolyun_leadcoeff_ref(G, ctx),
                                  fq_nmod_mpolyun_leadcoeff_ref(Abar, ctx), ctx->fqctx);
        FLINT_ASSERT(fq_nmod_poly_equal(modulus, leadA, ctx->fqctx));
        fq_nmod_poly_mul(modulus, fq_nmod_mpolyun_leadcoeff_ref(G, ctx),
                                  fq_nmod_mpolyun_leadcoeff_ref(Bbar, ctx), ctx->fqctx);
        FLINT_ASSERT(fq_nmod_poly_equal(modulus, leadB, ctx->fqctx));
    }
    fq_nmod_poly_clear(leadA, ctx->fqctx);
    fq_nmod_poly_clear(leadB, ctx->fqctx);
#endif

    fq_nmod_poly_clear(cA, ctx->fqctx);
    fq_nmod_poly_clear(cB, ctx->fqctx);
    fq_nmod_poly_clear(cG, ctx->fqctx);
    fq_nmod_poly_clear(cAbar, ctx->fqctx);
    fq_nmod_poly_clear(cBbar, ctx->fqctx);
    fq_nmod_poly_clear(trem, ctx->fqctx);
    fq_nmod_poly_clear(gamma, ctx->fqctx);

    fq_nmod_mpolyun_clear(T, ctx);

    fq_nmod_poly_clear(modulus, ctx->fqctx);

    fq_nmod_mpolyun_clear(Aeval, ectx);
    fq_nmod_mpolyun_clear(Beval, ectx);
    fq_nmod_mpolyun_clear(Geval, ectx);
    fq_nmod_mpolyun_clear(Abareval, ectx);
    fq_nmod_mpolyun_clear(Bbareval, ectx);
    fq_nmod_clear(gammaeval, ectx->fqctx);
    fq_nmod_clear(temp, ectx->fqctx);

    _fq_nmod_mpoly_embed_chooser_clear(embc, ectx, ctx, randstate);

    flint_randclear(randstate);

    return success;
}





void fq_nmod_mpoly_to_fq_nmod_poly_keepbits(fq_nmod_poly_t A, slong * Ashift,
            const fq_nmod_mpoly_t B, slong var, const fq_nmod_mpoly_ctx_t ctx);

void fq_nmod_mpoly_from_fq_nmod_poly_keepbits(fq_nmod_mpoly_t A,
         const fq_nmod_poly_t B, slong Bshift, slong var, mp_bitcnt_t bits,
                                                const fq_nmod_mpoly_ctx_t ctx);


int fq_nmod_mpoly_gcd_brownnew(fq_nmod_mpoly_t G, const fq_nmod_mpoly_t A,
                        const fq_nmod_mpoly_t B, const fq_nmod_mpoly_ctx_t ctx)
{
    int success;
    slong * perm;
    ulong * shift, * stride;
    slong i;
    mp_bitcnt_t new_bits;
    fq_nmod_mpoly_ctx_t uctx;
    fq_nmod_mpolyun_t An, Bn, Gn, Abarn, Bbarn;

    if (fq_nmod_mpoly_is_zero(A, ctx))
    {
        if (fq_nmod_mpoly_is_zero(B, ctx))
        {
            fq_nmod_mpoly_zero(G, ctx);
        }
        else
        {
            fq_nmod_mpoly_make_monic(G, B, ctx);
        }
        return 1;
    }

    if (fq_nmod_mpoly_is_zero(B, ctx))
    {
        fq_nmod_mpoly_make_monic(G, A, ctx);
        return 1;
    }

    if (A->bits > FLINT_BITS || B->bits > FLINT_BITS)
    {
        return 0;
    }

    if (ctx->minfo->nvars == 1)
    {
        slong shiftA, shiftB;
        fq_nmod_poly_t a, b, g;
        fq_nmod_poly_init(a, ctx->fqctx);
        fq_nmod_poly_init(b, ctx->fqctx);
        fq_nmod_poly_init(g, ctx->fqctx);
        fq_nmod_mpoly_to_fq_nmod_poly_keepbits(a, &shiftA, A, 0, ctx);
        fq_nmod_mpoly_to_fq_nmod_poly_keepbits(b, &shiftB, B, 0, ctx);
        fq_nmod_poly_gcd(g, a, b, ctx->fqctx);
        fq_nmod_mpoly_from_fq_nmod_poly_keepbits(G, g, FLINT_MIN(shiftA, shiftB), 0, A->bits, ctx);
        fq_nmod_poly_clear(a, ctx->fqctx);
        fq_nmod_poly_clear(b, ctx->fqctx);
        fq_nmod_poly_clear(g, ctx->fqctx);
        return 1;
    }

    perm = (slong *) flint_malloc(ctx->minfo->nvars*sizeof(slong));
    shift = (ulong *) flint_malloc(ctx->minfo->nvars*sizeof(ulong));
    stride = (ulong *) flint_malloc(ctx->minfo->nvars*sizeof(ulong));
    for (i = 0; i < ctx->minfo->nvars; i++)
    {
        perm[i] = i + 1 < ctx->minfo->nvars ? i + 1 : 0;
        shift[i] = 0;
        stride[i] = 1;
    }

    new_bits = FLINT_MAX(A->bits, B->bits);

    fq_nmod_mpoly_ctx_init(uctx, ctx->minfo->nvars - 1, ORD_LEX, ctx->fqctx);
    fq_nmod_mpolyun_init(An, new_bits, uctx);
    fq_nmod_mpolyun_init(Bn, new_bits, uctx);
    fq_nmod_mpolyun_init(Gn, new_bits, uctx);
    fq_nmod_mpolyun_init(Abarn, new_bits, uctx);
    fq_nmod_mpolyun_init(Bbarn, new_bits, uctx);

    fq_nmod_mpoly_to_mpolyun_perm_deflate(An, A, perm, shift, stride, uctx, ctx);
    fq_nmod_mpoly_to_mpolyun_perm_deflate(Bn, B, perm, shift, stride, uctx, ctx);
    success = fq_nmod_mpolyun_gcd_brown_smprime(Gn, Abarn, Bbarn, An, Bn, uctx->minfo->nvars - 1, uctx);
    if (!success)
    {
        fq_nmod_mpoly_to_mpolyun_perm_deflate(An, A, perm, shift, stride, uctx, ctx);
        fq_nmod_mpoly_to_mpolyun_perm_deflate(Bn, B, perm, shift, stride, uctx, ctx);
        success = fq_nmod_mpolyun_gcd_brown_lgprime(Gn, Abarn, Bbarn, An, Bn, uctx->minfo->nvars - 1, uctx);
    }

    if (success)
    {
        fq_nmod_mpoly_from_mpolyun_perm_inflate(G, new_bits, Gn, perm, shift, stride, uctx, ctx);
        fq_nmod_mpoly_make_monic(G, G, ctx);        
    }
    else
    {
printf("gcd failing\n");
    }

    fq_nmod_mpolyun_clear(An, uctx);
    fq_nmod_mpolyun_clear(Bn, uctx);
    fq_nmod_mpolyun_clear(Gn, uctx);
    fq_nmod_mpolyun_clear(Abarn, uctx);
    fq_nmod_mpolyun_clear(Bbarn, uctx);
    fq_nmod_mpoly_ctx_clear(uctx);

    flint_free(perm);
    flint_free(shift);
    flint_free(stride);

    return success;
}
