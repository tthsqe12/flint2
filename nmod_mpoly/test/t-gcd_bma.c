/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "nmod_mpoly.h"
#include "fmpz_mpoly.h"
#include "profiler.h"
void usleep(ulong);





typedef struct {
    fmpz_t p;
    fmpz_preinvn_t pinv;
} fmpz_mod_ctx_struct;
typedef fmpz_mod_ctx_struct fmpz_mod_ctx_t[1];

void fmpz_mod_ctx_init(fmpz_mod_ctx_t ctx, const fmpz_t p)
{
    fmpz_init_set(ctx->p, p);
    fmpz_preinvn_init(ctx->pinv, p);
}

void fmpz_mod_ctx_clear(fmpz_mod_ctx_t ctx)
{
    fmpz_preinvn_clear(ctx->pinv);
    fmpz_clear(ctx->p);
}

void fmpz_mod_ctx_set_mod(fmpz_mod_ctx_t ctx, const fmpz_t p)
{
    fmpz_mod_ctx_clear(ctx);
    fmpz_mod_ctx_init(ctx, p);
}

int fmpz_mod_is_canonical(const fmpz_t a, const fmpz_mod_ctx_t ctx)
{
    return fmpz_sgn(a) >= 0 && fmpz_cmp(a, ctx->p) < 0;
}

void fmpz_mod_add(fmpz_t a, const fmpz_t b, const fmpz_t c, const fmpz_mod_ctx_t ctx)
{
    FLINT_ASSERT(fmpz_mod_is_canonical(b, ctx));
    FLINT_ASSERT(fmpz_mod_is_canonical(c, ctx));
    fmpz_add(a, b, c);
    if (fmpz_cmpabs(a, ctx->p) >= 0)
    {
        fmpz_sub(a, a, ctx->p);
    }
    FLINT_ASSERT(fmpz_mod_is_canonical(a, ctx));
}

void fmpz_mod_sub(fmpz_t a, const fmpz_t b, const fmpz_t c, const fmpz_mod_ctx_t ctx)
{
    FLINT_ASSERT(fmpz_mod_is_canonical(b, ctx));
    FLINT_ASSERT(fmpz_mod_is_canonical(c, ctx));
    fmpz_sub(a, b, c);
    if (fmpz_sgn(a) < 0)
    {
        fmpz_add(a, a, ctx->p);
    }
    FLINT_ASSERT(fmpz_mod_is_canonical(a, ctx));
}

void fmpz_mod_neg(fmpz_t a, const fmpz_t b, const fmpz_mod_ctx_t ctx)
{
    FLINT_ASSERT(fmpz_mod_is_canonical(b, ctx));
    fmpz_neg(a, b);
    if (fmpz_sgn(a) < 0)
    {
        fmpz_add(a, a, ctx->p);
    }
    FLINT_ASSERT(fmpz_mod_is_canonical(a, ctx));
}

void fmpz_mod_mul(fmpz_t a, const fmpz_t b, const fmpz_t c, const fmpz_mod_ctx_t ctx)
{
/*
    fmpz_t q, t;
    FLINT_ASSERT(fmpz_mod_is_canonical(b, ctx));
    FLINT_ASSERT(fmpz_mod_is_canonical(c, ctx));
    fmpz_init(q);
    fmpz_init(t);
    fmpz_mul(t, b, c);
    fmpz_fdiv_qr_preinvn(q, a, t, ctx->p, ctx->pinv);
    fmpz_clear(q);
    fmpz_clear(t);
    FLINT_ASSERT(fmpz_mod_is_canonical(a, ctx));
*/

    fmpz_t t;
    FLINT_ASSERT(fmpz_mod_is_canonical(b, ctx));
    FLINT_ASSERT(fmpz_mod_is_canonical(c, ctx));
    fmpz_init(t);
    fmpz_mul(t, b, c);
    fmpz_mod(a, t, ctx->p);
    fmpz_clear(t);
    FLINT_ASSERT(fmpz_mod_is_canonical(a, ctx));

}

void fmpz_mod_inv(fmpz_t a, const fmpz_t b, const fmpz_mod_ctx_t ctx)
{
    fmpz_t d;
    FLINT_ASSERT(fmpz_mod_is_canonical(b, ctx));
    fmpz_init(d);
    fmpz_gcdinv(d, a, b, ctx->p);
/*
printf("inv p: "); fmpz_print(ctx->p); printf("\n");
printf("inv b: "); fmpz_print(b); printf("\n");
printf("inv a: "); fmpz_print(a); printf("\n");
printf("inv d: "); fmpz_print(d); printf("\n");
*/
    if (!fmpz_is_one(d))
    {
        flint_throw(FLINT_IMPINV, "Cannot invert in fmpz_mod_inv\n");        
    }
    fmpz_clear(d);

    fmpz_mod(a, a, ctx->p);

    FLINT_ASSERT(fmpz_mod_is_canonical(a, ctx));
}

void fmpz_mod_pow_ui(fmpz_t a, const fmpz_t b, ulong pow, const fmpz_mod_ctx_t ctx)
{
    fmpz_t r, x;

    fmpz_init_set_ui(r, 1);
    fmpz_init_set(x, b);

    while (pow != 0)
    {
        if ((pow & 1) != 0)
        {
            fmpz_mod_mul(r, r, x, ctx);
        }
        pow = pow >> 1;
        if (pow != 0)
        {
            fmpz_mod_mul(x, x, x, ctx);
        }
    }

    fmpz_swap(a, r);
    fmpz_clear(r);
    fmpz_clear(x);
    FLINT_ASSERT(fmpz_mod_is_canonical(a, ctx));
}

void fmpz_mod_pow_fmpz(fmpz_t a, const fmpz_t b, const fmpz_t pow, const fmpz_mod_ctx_t ctx)
{
    mp_bitcnt_t i, bits;
    fmpz_t r, x;

    fmpz_init_set_ui(r, 1);
    if (fmpz_sgn(pow) < 0)
    {
        fmpz_init(x);
        fmpz_mod_inv(x, b, ctx);
    }
    else
    {
        fmpz_init_set(x, b);
    }

    bits = fmpz_bits(pow);

    for (i = 0; i < bits; i++)
    {
        if (fmpz_tstbit(pow, i) != 0)
        {
            fmpz_mod_mul(r, r, x, ctx);
        }
        fmpz_mod_mul(x, x, x, ctx);
    }

    fmpz_swap(a, r);

    fmpz_clear(r);
    fmpz_clear(x);
    FLINT_ASSERT(fmpz_mod_is_canonical(a, ctx));
}





void nmod_poly_init_mod(nmod_poly_t poly, const nmod_t mod)
{
    poly->coeffs = NULL;
    poly->alloc = 0;
    poly->length = 0;
    poly->mod = mod;
}

void nmod_poly_set_mod(nmod_poly_t poly, const nmod_t mod)
{
    poly->mod = mod;
}



void nmod_mpoly_ctx_init_mpoly_ctx(
    nmod_mpoly_ctx_t ctx,
    const mpoly_ctx_t mctx,
    mp_limb_t modulus)
{
    mpoly_ctx_init(ctx->minfo, mctx->nvars, mctx->ord);
    nmodf_ctx_init(ctx->ffinfo, modulus);
}

/*
void nmod_mpoly_ctx_set_mod(nmod_mpoly_ctx_t ctx, mp_limb_t modulus)
{
    nmodf_ctx_clear(ctx->ffinfo);
    nmodf_ctx_init(ctx->ffinfo, modulus);
}
*/
void nmod_mpoly_ctx_set_mod(
    nmod_mpoly_ctx_t ctx,
    mp_limb_t p)
{
    nmodf_ctx_reset(ctx->ffinfo, p);
}






void nmod_mpolyun_set_mod(nmod_mpolyun_t A, const nmod_t mod)
{
    slong i, j;

    for (i = 0; i < A->alloc; i++)
    {
        nmod_mpolyn_struct * Ac = A->coeffs + i;
        for (j = 0; j < Ac->alloc; j++)
        {
            (Ac->coeffs + j)->mod = mod;
        }
    }
}

void nmod_mpolyn_scalar_mul_nmod(nmod_mpolyn_t A, mp_limb_t c, const nmod_mpoly_ctx_t ctx)
{
    slong i;
    for (i = 0; i < A->length; i++)
    {
        nmod_poly_scalar_mul_nmod(A->coeffs + i, A->coeffs + i, c);
    }
}

void nmod_mpolyun_scalar_mul_nmod(nmod_mpolyun_t A, mp_limb_t c, const nmod_mpoly_ctx_t ctx)
{
    slong i;
    FLINT_ASSERT(c != 0);
    for (i = 0; i < A->length; i++)
    {
        nmod_mpolyn_scalar_mul_nmod(A->coeffs + i, c, ctx);
    }
}

void nmod_mpolyn_one(nmod_mpolyn_t A, const nmod_mpoly_ctx_t ctx)
{
    nmod_poly_struct * Acoeff;
    ulong * Aexp;
    slong N;

    nmod_mpolyn_fit_length(A, 1, ctx);
    Acoeff = A->coeffs;
    Aexp = A->exps;

    N = mpoly_words_per_exp(A->bits, ctx->minfo);

    nmod_poly_one(Acoeff + 0);
    mpoly_monomial_zero(Aexp + N*0, N);

    A->length = 1;
}


void nmod_mpolyun_one(nmod_mpolyun_t A, const nmod_mpoly_ctx_t ctx)
{
    nmod_mpolyun_fit_length(A, 1, ctx);
    nmod_mpolyn_one(A->coeffs + 0, ctx);
    A->exps[0] = 0;
    A->length = 1;
}


/*
    get the leading exponent in x_0,...,x_var
    A is in R[x_0, ... x_(var-1)][x_var]
*/
mp_limb_t nmod_mpolyn_leadcoeff_last(nmod_mpolyn_t A, const nmod_mpoly_ctx_t ctx)
{
    nmod_poly_struct * leadpoly;

    FLINT_ASSERT(A->length > 0);
    FLINT_ASSERT(nmod_poly_degree(A->coeffs + 0) >= 0);

    leadpoly = A->coeffs + 0;
    return leadpoly->coeffs[leadpoly->length - 1];
}

mp_limb_t nmod_mpolyun_leadcoeff_last(nmod_mpolyun_t A, const nmod_mpoly_ctx_t ctx)
{
    FLINT_ASSERT(A->length > 0);
    return nmod_mpolyn_leadcoeff_last(A->coeffs + 0, ctx);
}


void nmod_mpolyun_mul_last(nmod_mpolyun_t A, nmod_poly_t b,
                                                    const nmod_mpoly_ctx_t ctx)
{
    slong i, j;
    nmod_poly_t t;

    nmod_poly_init_mod(t, ctx->ffinfo->mod);

    for (i = 0; i < A->length; i++)
    {
        for (j = 0; j < (A->coeffs + i)->length; j++)
        {
            nmod_poly_mul(t, (A->coeffs + i)->coeffs + j, b);
            nmod_poly_swap(t, (A->coeffs + i)->coeffs + j);
        }
    }

    nmod_poly_clear(t);
}




/*
    E = A(v = alpha)
    A is in R[X][v]
    E is in R[X]
*/
void nmod_mpolyun_eval_last_bivar(nmod_poly_t E, const nmod_mpolyun_t A,
                                  mp_limb_t alpha, const nmod_mpoly_ctx_t ctx)
{
    mp_limb_t v;
    slong Ai, Alen;
    nmod_mpolyn_struct * Acoeff;
    ulong * Aexp;

    Acoeff = A->coeffs;
    Aexp = A->exps;
    Alen = A->length;
    Ai = 0;
    nmod_poly_zero(E);
    for (Ai = 0; Ai < Alen; Ai++)
    {
        v = nmod_poly_evaluate_nmod((Acoeff + Ai)->coeffs + 0, alpha);
        nmod_poly_set_coeff_ui(E, Aexp[Ai], v);
    }
}



/*
    A = B
    A, B are in R[X]
*/
void nmod_mpolyun_set_poly(nmod_mpolyun_t A, const nmod_poly_t B,
                                                   const nmod_mpoly_ctx_t ctx)
{
    slong N = mpoly_words_per_exp_sp(A->bits, ctx->minfo);
    slong Bexp;
    slong Blen = nmod_poly_length(B);
    mp_limb_t * Bcoeff = B->coeffs;
    nmod_mpolyn_struct * Acoeff;
    ulong * Aexp;
    slong Ai;

    nmod_mpolyun_fit_length(A, Blen, ctx);
    Acoeff = A->coeffs;
    Aexp = A->exps;

    Ai = 0;
    for (Bexp = Blen - 1; Bexp >= 0; Bexp--)
    {
        if (Bcoeff[Bexp] != UWORD(0))
        {
            FLINT_ASSERT(Ai < A->alloc);

            nmod_mpolyn_fit_length(Acoeff + Ai, 1, ctx);
            mpoly_monomial_zero((Acoeff + Ai)->exps + N*0, N);
            nmod_poly_zero((Acoeff + Ai)->coeffs + 0);
            nmod_poly_set_coeff_ui((Acoeff + Ai)->coeffs + 0, 0, Bcoeff[Bexp]);
            Aexp[Ai] = Bexp;
            (Acoeff + Ai)->length = 1;
            Ai++;
        }
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
int nmod_mpolyun_addinterp_bivar(slong * lastdeg_,
                   nmod_mpolyun_t F, nmod_mpolyun_t T, const nmod_poly_t A,
      const nmod_poly_t modulus,  mp_limb_t alpha,  const nmod_mpoly_ctx_t ctx)
{
    int changed = 0;
    slong lastdeg = -WORD(1);
    slong N = mpoly_words_per_exp_sp(F->bits, ctx->minfo);
    mp_limb_t v;
    slong Fi, Toff, Aexp;
    mp_limb_t * Acoeff = A->coeffs;
    slong Flen = F->length;
    nmod_mpolyn_struct * Fcoeff = F->coeffs;
    ulong * Fexp = F->exps;
    nmod_mpolyn_struct * Tcoeff;
    ulong * Texp;
    nmod_poly_t tp;
    
    Fi = 0;

/*
flint_printf("\nnmod_mpolyun_addinterp_bivar: (alpha = %wu)\n", alpha);
flint_printf("A: "); nmod_poly_print_pretty(A, "X"); printf("\n");
flint_printf("F: "); nmod_mpolyun_print_pretty(F, NULL, ctx); printf("\n");
*/

    Aexp = nmod_poly_degree(A);

    nmod_poly_init(tp, ctx->ffinfo->mod.n);

    nmod_mpolyun_fit_length(T, Flen + Aexp + 1, ctx);
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
            FLINT_ASSERT(!nmod_poly_is_zero((Fcoeff + Fi)->coeffs + 0));
            FLINT_ASSERT(nmod_poly_degree((Fcoeff + Fi)->coeffs + 0) < nmod_poly_degree(modulus));
        }

        if (Aexp >= 0)
        {
            FLINT_ASSERT(Acoeff[Aexp] != UWORD(0));
        }

        nmod_mpolyn_fit_length(Tcoeff + Toff, 1, ctx);

        if (Fi < Flen && Aexp >= 0 && Fexp[Fi] == Aexp)
        {
            /* F term ok, A term ok */
            v = nmod_poly_evaluate_nmod((Fcoeff + Fi)->coeffs + 0, alpha);
            v = nmod_sub(Acoeff[Aexp], v, ctx->ffinfo->mod);
            if (v != UWORD(0))
            {
                changed = 1;
                nmod_poly_scalar_mul_nmod(tp, modulus, v);
                nmod_poly_add((Tcoeff + Toff)->coeffs + 0, (Fcoeff + Fi)->coeffs + 0, tp);
            }
            else
            {
                nmod_poly_set((Tcoeff + Toff)->coeffs + 0, (Fcoeff + Fi)->coeffs + 0);                
            }
            Texp[Toff] = Aexp;
            Fi++;
            do {
                Aexp--;
            } while (Aexp >= 0 && Acoeff[Aexp] == UWORD(0));
        }
        else if (Fi < Flen && (Aexp < 0 || Fexp[Fi] > Aexp))
        {
            /* F term ok, A term missing */
            v = nmod_poly_evaluate_nmod((Fcoeff + Fi)->coeffs + 0, alpha);
            if (v != UWORD(0))
            {
                changed = 1;
                nmod_poly_scalar_mul_nmod(tp, modulus, v);
                nmod_poly_sub((Tcoeff + Toff)->coeffs + 0, (Fcoeff + Fi)->coeffs + 0, tp);
            }
            else
            {
                nmod_poly_set((Tcoeff + Toff)->coeffs + 0, (Fcoeff + Fi)->coeffs + 0);                
            }
            Texp[Toff] = Fexp[Fi];
            Fi++;
        }
        else if (Aexp >= 0 && (Fi >= Flen || Fexp[Fi] < Aexp))
        {
            /* F term missing, A term ok */
            changed = 1;
            nmod_poly_scalar_mul_nmod((Tcoeff + Toff)->coeffs + 0, modulus, Acoeff[Aexp]);
            Texp[Toff] = Aexp;
            do {
                Aexp--;
            } while (Aexp >= 0 && Acoeff[Aexp] == UWORD(0));
        }
        else
        {
            FLINT_ASSERT(0);
        }

        lastdeg = FLINT_MAX(lastdeg, nmod_poly_degree((Tcoeff + Toff)->coeffs + 0));
        mpoly_monomial_zero((Tcoeff + Toff)->exps + N*0, N);
        FLINT_ASSERT(!nmod_poly_is_zero((Tcoeff + Toff)->coeffs + 0));
        (Tcoeff + Toff)->length = 1;
        Toff++;
    }
    T->length = Toff;

    nmod_poly_clear(tp);

    if (changed)
    {
        nmod_mpolyun_swap(T, F);
    }

    *lastdeg_ = lastdeg;
    return changed;
}





/*
    set vp = P(alpha), vm = P(-alpha) given powers of alpha
*/
void _nmod_poly_eval2_pow(mp_limb_t * vp, mp_limb_t * vm, nmod_poly_t P, 
                                  nmod_poly_t alphapow, const nmodf_ctx_t fctx)
{
    mp_limb_t * Pcoeffs = P->coeffs;
    slong Plen = P->length;
    mp_limb_t * alpha_powers = alphapow->coeffs;
    mp_limb_t p1, p0, a0, a1, a2, q1, q0, b0, b1, b2;
    slong k;

    a0 = a1 = a2 = UWORD(0);
    b0 = b1 = b2 = UWORD(0);

    if (Plen > alphapow->length)
    {
        slong oldlength = alphapow->length;
        FLINT_ASSERT(2 <= oldlength);
        nmod_poly_fit_length(alphapow, Plen);
        for (k = oldlength; k < Plen; k++)
        {
            alphapow->coeffs[k] = nmod_mul(alphapow->coeffs[k - 1],
                                               alphapow->coeffs[1], fctx->mod);
        }
        alphapow->length = Plen;
    }

    for (k = 0; k + 2 <= Plen; k += 2)
    {
        umul_ppmm(p1, p0, Pcoeffs[k + 0], alpha_powers[k + 0]);
        umul_ppmm(q1, q0, Pcoeffs[k + 1], alpha_powers[k + 1]);
        add_sssaaaaaa(a2, a1, a0, a2, a1, a0, WORD(0), p1, p0);
        add_sssaaaaaa(b2, b1, b0, b2, b1, b0, WORD(0), q1, q0);
    }

    if (k < Plen)
    {
        umul_ppmm(p1, p0, Pcoeffs[k + 0], alpha_powers[k + 0]);
        add_sssaaaaaa(a2, a1, a0, a2, a1, a0, WORD(0), p1, p0);
        k++;
    }

    FLINT_ASSERT(k == Plen);

    NMOD_RED3(p0, a2, a1, a0, fctx->mod);
    NMOD_RED3(q0, b2, b1, b0, fctx->mod);

    vp[0] = nmod_add(p0, q0, fctx->mod);
    vm[0] = nmod_sub(p0, q0, fctx->mod);
}



/*
    E = A(v = alpha), F = A(v = -alpha)
    A is in R[X][v]
    E is in R[X]
    F is in R[X]
*/
void nmod_mpolyun_eval2_last_bivar(nmod_poly_t E, nmod_poly_t F,
                      const nmod_mpolyun_t A, nmod_poly_t alphapow,
                                                    const nmod_mpoly_ctx_t ctx)
{
    mp_limb_t u, v;
    slong Ai, Alen;
    nmod_mpolyn_struct * Acoeff;
    ulong * Aexp;

    Acoeff = A->coeffs;
    Aexp = A->exps;
    Alen = A->length;
    Ai = 0;
    nmod_poly_zero(E);
    nmod_poly_zero(F);
    for (Ai = 0; Ai < Alen; Ai++)
    {
        _nmod_poly_eval2_pow(&u, &v, (Acoeff + Ai)->coeffs + 0, alphapow, ctx->ffinfo);
        nmod_poly_set_coeff_ui(E, Aexp[Ai], u);
        nmod_poly_set_coeff_ui(F, Aexp[Ai], v);
    }
}

/*
    set F from its value A at v = alpha and its value B at v = -alpha
    no assumptions about matching monomials
    F is in R[X][v]
    A is in R[X]
    B is in R[X]
*/
void nmod_mpolyun_startinterp2_bivar(slong * lastdeg_,
                   nmod_mpolyun_t F, const nmod_poly_t A, const nmod_poly_t B,
                                  mp_limb_t alpha,  const nmod_mpoly_ctx_t ctx)
{
    slong lastdeg = -WORD(1);
    slong N = mpoly_words_per_exp_sp(F->bits, ctx->minfo);
    mp_limb_t u, v, d0, d1, Avalue, Bvalue;
    slong Fi, Aexp, Bexp;
    mp_limb_t * Acoeff = A->coeffs;
    mp_limb_t * Bcoeff = B->coeffs;
    nmod_mpolyn_struct * Fcoeff;
    ulong * Fexp;
    slong e;

    Aexp = nmod_poly_degree(A);
    Bexp = nmod_poly_degree(B);

    nmod_mpolyun_fit_length(F, FLINT_MAX(Aexp, Bexp) + 1, ctx);
    Fcoeff = F->coeffs;
    Fexp = F->exps;

    d0 = n_invmod(UWORD(2), ctx->ffinfo->mod.n);
    d1 = n_invmod(nmod_add(alpha, alpha, ctx->ffinfo->mod), ctx->ffinfo->mod.n);

    Fi = 0;
    while (Aexp >= 0 || Bexp >= 0)
    {
        e = Aexp;
        Avalue = 0;
        Bvalue = 0;
        if (Aexp == Bexp)
        {
            Avalue = Acoeff[Aexp];
            Bvalue = Bcoeff[Bexp];
        }
        else if (Aexp > Bexp)
        {
            Avalue = Acoeff[Aexp];
        }
        else
        {
            FLINT_ASSERT(Bexp > Aexp);
            e = Bexp;
            Bvalue = Bcoeff[Bexp];
        }
        FLINT_ASSERT(Avalue != 0 || Bvalue != 0);
        u = nmod_add(Avalue, Bvalue, ctx->ffinfo->mod);
        v = nmod_sub(Avalue, Bvalue, ctx->ffinfo->mod);
        u = nmod_mul(u, d0, ctx->ffinfo->mod);
        v = nmod_mul(v, d1, ctx->ffinfo->mod);

        FLINT_ASSERT(Fi < F->alloc);
        nmod_mpolyn_fit_length(Fcoeff + Fi, 1, ctx);
        mpoly_monomial_zero((Fcoeff + Fi)->exps + N*0, N);
        nmod_poly_zero((Fcoeff + Fi)->coeffs + 0);
        nmod_poly_set_coeff_ui((Fcoeff + Fi)->coeffs + 0, 0, u);
        nmod_poly_set_coeff_ui((Fcoeff + Fi)->coeffs + 0, 1, v);
        lastdeg = FLINT_MAX(lastdeg, nmod_poly_degree((Fcoeff + Fi)->coeffs + 0));
        Fexp[Fi] = e;
        (Fcoeff + Fi)->length = 1;
        Fi++;

        if (e == Aexp)
        {
            do {
                Aexp--;
            } while (Aexp >= 0 && Acoeff[Aexp] == 0);
        }
        if (e == Bexp)
        {
            do {
                Bexp--;
            } while (Bexp >= 0 && Bcoeff[Bexp] == 0);
        }
    }
    F->length = Fi;
/*
pthread_mutex_lock(&iomutex);
flint_printf("startinterp2_bivar returning  alpha = %wu\n", alpha);
flint_printf("A: "); nmod_poly_print_pretty(A, "X"); printf("\n");
flint_printf("B: "); nmod_poly_print_pretty(B, "X"); printf("\n");
flint_printf("F: "); nmod_mpolyun_print_pretty(F, NULL, ctx); printf("\n");
pthread_mutex_unlock(&iomutex);
*/
    *lastdeg_ = lastdeg;
    return;
}


/*
    update F from its value A at v = alpha and its value B at v = -alpha
    no assumptions about matching monomials
    F is in R[X][v]
    A is in R[X]
    B is in R[X]
    it is expected that modulus(alpha) == modulus(-alpha) == 1/(2*alpha)
*/
int nmod_mpolyun_addinterp2_bivar(slong * lastdeg_,
                   nmod_mpolyun_t F, nmod_mpolyun_t T, const nmod_poly_t A, const nmod_poly_t B,
      const nmod_poly_t modulus, nmod_poly_t alphapow, const nmod_mpoly_ctx_t ctx)
{
    int changed = 0, Finc;
    mp_limb_t alpha = nmod_poly_get_coeff_ui(alphapow, 1);
    slong lastdeg = -WORD(1);
    slong N = mpoly_words_per_exp_sp(F->bits, ctx->minfo);
    mp_limb_t u, v, FvalueA, FvalueB;
    slong Fi, Toff, Aexp, Bexp, e;
    mp_limb_t * Acoeff = A->coeffs;
    mp_limb_t * Bcoeff = B->coeffs;
    slong Flen = F->length;
    nmod_mpolyn_struct * Fcoeff = F->coeffs;
    ulong * Fexp = F->exps;
    nmod_mpolyn_struct * Tcoeff;
    ulong * Texp;
    nmod_poly_t tp;

    Fi = 0;
    Aexp = nmod_poly_degree(A);
    Bexp = nmod_poly_degree(B);

    nmod_poly_init(tp, ctx->ffinfo->mod.n);

    nmod_mpolyun_fit_length(T, Flen + FLINT_MAX(Aexp, Bexp) + 1, ctx);
    Tcoeff = T->coeffs;
    Texp = T->exps;
    Toff = 0;

#if WANT_ASSERT
/*
    u = nmod_poly_evaluate_nmod(modulus, alpha);
    u = nmod_mul(u, alpha, ctx->ffinfo->mod);
    u = nmod_mul(u, 2, ctx->ffinfo->mod);
    FLINT_ASSERT(u == 1);
    u = nmod_poly_evaluate_nmod(modulus, ctx->ffinfo->mod.n - alpha);
    u = nmod_mul(u, alpha, ctx->ffinfo->mod);
    u = nmod_mul(u, 2, ctx->ffinfo->mod);
    FLINT_ASSERT(u == 1);
*/
#endif

/*
pthread_mutex_lock(&iomutex);
flint_printf("addinterp2_bivar called  alpha = %wu\n", alpha);
flint_printf("A: "); nmod_poly_print_pretty(A, "X"); printf("\n");
flint_printf("B: "); nmod_poly_print_pretty(B, "X"); printf("\n");
flint_printf("F: "); nmod_mpolyun_print_pretty(F, NULL, ctx); printf("\n");
*/

    while (Fi < Flen || Aexp >= 0 || Bexp >= 0)
    {
        FLINT_ASSERT(Toff < T->alloc);

        e = -WORD(1);
        if (Fi < Flen)
        {
            e = Fexp[Fi];
            FLINT_ASSERT(!nmod_poly_is_zero((Fcoeff + Fi)->coeffs + 0));
            FLINT_ASSERT(nmod_poly_degree((Fcoeff + Fi)->coeffs + 0) < nmod_poly_degree(modulus));
        }
        if (Aexp >= 0)
        {
            e = FLINT_MAX(e, Aexp);
            FLINT_ASSERT(Acoeff[Aexp] != UWORD(0));
        }
        if (Bexp >= 0)
        {
            e = FLINT_MAX(e, Bexp);
            FLINT_ASSERT(Bcoeff[Bexp] != UWORD(0));
        }

        FLINT_ASSERT(e >= 0);
        nmod_mpolyn_fit_length(Tcoeff + Toff, 1, ctx);
        Texp[Toff] = e;

        FvalueA = FvalueB = 0;
        Finc = 0;
        if (Fi < Flen && e == Fexp[Fi])
        {
            Finc = 1;
            _nmod_poly_eval2_pow(&FvalueA, &FvalueB, (Fcoeff + Fi)->coeffs + 0, alphapow, ctx->ffinfo);
        }

        if (e == Aexp)
        {
            FvalueA = nmod_sub(FvalueA, Acoeff[Aexp], ctx->ffinfo->mod);
        }
        if (e == Bexp)
        {
            FvalueB = nmod_sub(FvalueB, Bcoeff[Bexp], ctx->ffinfo->mod);
        }

        u = nmod_sub(FvalueB, FvalueA, ctx->ffinfo->mod);
        v = nmod_mul(ctx->ffinfo->mod.n - alpha, nmod_add(FvalueB, FvalueA, ctx->ffinfo->mod), ctx->ffinfo->mod);

        if (u != 0 || v != 0)
        {
            changed = 1;
            nmod_poly_set_coeff_ui(tp, 0, v);
            nmod_poly_set_coeff_ui(tp, 1, u);
            nmod_poly_mul_classical((Tcoeff + Toff)->coeffs + 0, modulus, tp);
            if (Finc)
            {
                nmod_poly_add((Tcoeff + Toff)->coeffs + 0, (Tcoeff + Toff)->coeffs + 0, (Fcoeff + Fi)->coeffs + 0);
            }
        }
        else
        {
            FLINT_ASSERT(Finc == 1);
            nmod_poly_set((Tcoeff + Toff)->coeffs + 0, (Fcoeff + Fi)->coeffs + 0);
        }

        lastdeg = FLINT_MAX(lastdeg, nmod_poly_degree((Tcoeff + Toff)->coeffs + 0));
        mpoly_monomial_zero((Tcoeff + Toff)->exps + N*0, N);
        FLINT_ASSERT(!nmod_poly_is_zero((Tcoeff + Toff)->coeffs + 0));
        (Tcoeff + Toff)->length = 1;
        Toff++;

        Fi += Finc;
        if (e == Aexp)
        {
            do {
                Aexp--;
            } while (Aexp >= 0 && Acoeff[Aexp] == 0);
        }
        if (e == Bexp)
        {
            do {
                Bexp--;
            } while (Bexp >= 0 && Bcoeff[Bexp] == 0);
        }
    }
    T->length = Toff;

    nmod_poly_clear(tp);

    if (changed)
    {
        nmod_mpolyun_swap(T, F);
    }
/*
flint_printf("returning F: "); nmod_mpolyun_print_pretty(F, NULL, ctx); printf("\n");
pthread_mutex_unlock(&iomutex);
*/

    *lastdeg_ = lastdeg;
    return changed;
}





int nmod_mpolyun_gcd_brown_smprime_bivar_ref(nmod_mpolyun_t G,
  nmod_mpolyun_t Abar, nmod_mpolyun_t Bbar, nmod_mpolyun_t A, nmod_mpolyun_t B,
                       const nmod_mpoly_ctx_t ctx)
{
    int success;
    slong bound;
    mp_limb_t alpha, temp, gammaeval;
    nmod_poly_t Aeval, Beval, Geval, Abareval, Bbareval;
    nmod_mpolyun_t T;
    slong deggamma, ldegG, ldegAbar, ldegBbar, ldegA, ldegB;
    nmod_poly_t cA, cB, cG, cAbar, cBbar, gamma;
    nmod_poly_t modulus, modulus2;
#if WANT_ASSERT
    nmod_poly_t leadA, leadB;
#endif

#if WANT_ASSERT
    nmod_poly_init(leadA, ctx->ffinfo->mod.n);
    nmod_poly_init(leadB, ctx->ffinfo->mod.n);
    nmod_poly_set(leadA, nmod_mpolyun_leadcoeff_ref(A, ctx));
    nmod_poly_set(leadB, nmod_mpolyun_leadcoeff_ref(B, ctx));
#endif

    nmod_poly_init(cA, ctx->ffinfo->mod.n);
    nmod_poly_init(cB, ctx->ffinfo->mod.n);
    nmod_mpolyun_content_last(cA, A, ctx);
    nmod_mpolyun_content_last(cB, B, ctx);
    nmod_mpolyun_divexact_last(A, cA, ctx);
    nmod_mpolyun_divexact_last(B, cB, ctx);

    nmod_poly_init(cG, ctx->ffinfo->mod.n);
    nmod_poly_gcd_euclidean(cG, cA, cB);

    nmod_poly_init(cAbar, ctx->ffinfo->mod.n);
    nmod_poly_init(cBbar, ctx->ffinfo->mod.n);
    nmod_poly_div(cAbar, cA, cG);
    nmod_poly_div(cBbar, cB, cG);

    nmod_poly_init(gamma, ctx->ffinfo->mod.n);
    nmod_poly_gcd(gamma, nmod_mpolyun_leadcoeff_ref(A, ctx),
                         nmod_mpolyun_leadcoeff_ref(B, ctx));

    ldegA = nmod_mpolyun_lastdeg(A, ctx);
    ldegB = nmod_mpolyun_lastdeg(B, ctx);
    deggamma = nmod_poly_degree(gamma);
    bound = 1 + deggamma + FLINT_MAX(ldegA, ldegB);

    nmod_poly_init(Aeval, ctx->ffinfo->mod.n);
    nmod_poly_init(Beval, ctx->ffinfo->mod.n);
    nmod_poly_init(Geval, ctx->ffinfo->mod.n);
    nmod_poly_init(Abareval, ctx->ffinfo->mod.n);
    nmod_poly_init(Bbareval, ctx->ffinfo->mod.n);

    nmod_mpolyun_init(T, A->bits, ctx);

    nmod_poly_init(modulus, ctx->ffinfo->mod.n);
    nmod_poly_init(modulus2, ctx->ffinfo->mod.n);
    nmod_poly_one(modulus);

    if ((ctx->ffinfo->mod.n & UWORD(1)) == UWORD(0))
    {
        success = 0;
        goto cleanup;
    }

    alpha = (ctx->ffinfo->mod.n - UWORD(1))/UWORD(2);

choose_prime: /* prime is v - alpha */

    if (alpha < 2)
    {
        success = 0;
        goto cleanup;
    }

    alpha--;

    /* make sure evaluation point does not kill both lc(A) and lc(B) */
    gammaeval = nmod_poly_evaluate_nmod(gamma, alpha);
    if (gammaeval == 0)
    {
        goto choose_prime;
    }

    /* evaluation point should kill neither A nor B */
    nmod_mpolyun_eval_last_bivar(Aeval, A, alpha, ctx);
    nmod_mpolyun_eval_last_bivar(Beval, B, alpha, ctx);
    FLINT_ASSERT(Aeval->length > 0);
    FLINT_ASSERT(Beval->length > 0);

    nmod_poly_gcd(Geval, Aeval, Beval);
    nmod_poly_div(Abareval, Aeval, Geval);
    nmod_poly_div(Bbareval, Beval, Geval);

    FLINT_ASSERT(Geval->length > 0);
    FLINT_ASSERT(Abareval->length > 0);
    FLINT_ASSERT(Bbareval->length > 0);

    if (nmod_poly_degree(Geval) == 0)
    {
        nmod_mpolyun_one(G, ctx);
        nmod_mpolyun_swap(Abar, A);
        nmod_mpolyun_swap(Bbar, B);
        goto successful_put_content;    
    }

    if (nmod_poly_degree(modulus) > 0)
    {
        FLINT_ASSERT(G->length > 0);
        if (nmod_poly_degree(Geval) > G->exps[0])
        {
            goto choose_prime;
        }
        else if (nmod_poly_degree(Geval) < G->exps[0])
        {
            nmod_poly_one(modulus);
        }
    }

    nmod_poly_scalar_mul_nmod(Geval, Geval, gammaeval);

    if (nmod_poly_degree(modulus) > 0)
    {
        temp = nmod_poly_evaluate_nmod(modulus, alpha);
        temp = n_invmod(temp, ctx->ffinfo->mod.n);
        nmod_poly_scalar_mul_nmod(modulus, modulus, temp);
        nmod_mpolyun_addinterp_bivar(&ldegG, G, T, Geval, modulus, alpha, ctx);
        nmod_mpolyun_addinterp_bivar(&ldegAbar, Abar, T, Abareval, modulus, alpha, ctx);
        nmod_mpolyun_addinterp_bivar(&ldegBbar, Bbar, T, Bbareval, modulus, alpha, ctx);
    }
    else
    {
        nmod_mpolyun_set_poly(G, Geval, ctx);
        nmod_mpolyun_set_poly(Abar, Abareval, ctx);
        nmod_mpolyun_set_poly(Bbar, Bbareval, ctx);
        ldegG = 0;
        ldegAbar = 0;
        ldegBbar = 0;
    }

    nmod_poly_scalar_mul_nmod(modulus2, modulus, alpha);
    nmod_poly_shift_left(modulus, modulus, 1);
    nmod_poly_sub(modulus, modulus, modulus2);

    if (nmod_poly_degree(modulus) < bound)
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

    nmod_poly_one(modulus);
    goto choose_prime;

successful:

    nmod_mpolyun_content_last(modulus, G, ctx);
    nmod_mpolyun_divexact_last(G, modulus, ctx);
    nmod_mpolyun_divexact_last(Abar, nmod_mpolyun_leadcoeff_ref(G, ctx), ctx);
    nmod_mpolyun_divexact_last(Bbar, nmod_mpolyun_leadcoeff_ref(G, ctx), ctx);

successful_put_content:

    nmod_mpolyun_mul_last(G, cG, ctx);
    nmod_mpolyun_mul_last(Abar, cAbar, ctx);
    nmod_mpolyun_mul_last(Bbar, cBbar, ctx);

    success = 1;

cleanup:

#if WANT_ASSERT
    if (success)
    {
        FLINT_ASSERT(1 == nmod_mpolyun_leadcoeff_last(G, ctx));
        nmod_poly_mul(modulus, nmod_mpolyun_leadcoeff_ref(G, ctx),
                               nmod_mpolyun_leadcoeff_ref(Abar, ctx));
        FLINT_ASSERT(nmod_poly_equal(modulus, leadA));
        nmod_poly_mul(modulus, nmod_mpolyun_leadcoeff_ref(G, ctx),
                               nmod_mpolyun_leadcoeff_ref(Bbar, ctx));
        FLINT_ASSERT(nmod_poly_equal(modulus, leadB));
    }
    nmod_poly_clear(leadA);
    nmod_poly_clear(leadB);
#endif

    nmod_poly_clear(cA);
    nmod_poly_clear(cB);
    nmod_poly_clear(cG);
    nmod_poly_clear(cAbar);
    nmod_poly_clear(cBbar);

    nmod_poly_clear(gamma);

    nmod_poly_clear(Aeval);
    nmod_poly_clear(Beval);
    nmod_poly_clear(Geval);
    nmod_poly_clear(Abareval);
    nmod_poly_clear(Bbareval);

    nmod_mpolyun_clear(T, ctx);

    nmod_poly_clear(modulus);
    nmod_poly_clear(modulus2);

    return success;
}




int nmod_mpolyun_gcd_brown_smprime_bivar(nmod_mpolyun_t G,
  nmod_mpolyun_t Abar, nmod_mpolyun_t Bbar, nmod_mpolyun_t A, nmod_mpolyun_t B,
                                                    const nmod_mpoly_ctx_t ctx)
{
    int success;
    slong bound;
    mp_limb_t alpha, temp, gammaevalp, gammaevalm;
    nmod_poly_t Aevalp, Bevalp, Gevalp, Abarevalp, Bbarevalp;
    nmod_poly_t Aevalm, Bevalm, Gevalm, Abarevalm, Bbarevalm;
    nmod_mpolyun_t T;
    slong deggamma, ldegG, ldegAbar, ldegBbar, ldegA, ldegB;
    nmod_poly_t cA, cB, cG, cAbar, cBbar, gamma;
    nmod_poly_t modulus, modulus2, alphapow, r;
    int gstab, astab, bstab, use_stab;
#if WANT_ASSERT
/*
    nmod_poly_t leadA, leadB;
*/
#endif

#if WANT_ASSERT
/*
    nmod_poly_init(leadA, ctx->ffinfo->mod.n);
    nmod_poly_init(leadB, ctx->ffinfo->mod.n);
    nmod_poly_set(leadA, nmod_mpolyun_leadcoeff_ref(A, ctx));
    nmod_poly_set(leadB, nmod_mpolyun_leadcoeff_ref(B, ctx));
*/
#endif

    nmod_poly_init(cA, ctx->ffinfo->mod.n);
    nmod_poly_init(cB, ctx->ffinfo->mod.n);
    nmod_mpolyun_content_last(cA, A, ctx);
    nmod_mpolyun_content_last(cB, B, ctx);
    nmod_mpolyun_divexact_last(A, cA, ctx);
    nmod_mpolyun_divexact_last(B, cB, ctx);

    nmod_poly_init(cG, ctx->ffinfo->mod.n);
    nmod_poly_gcd_euclidean(cG, cA, cB);

    nmod_poly_init(cAbar, ctx->ffinfo->mod.n);
    nmod_poly_init(cBbar, ctx->ffinfo->mod.n);
    nmod_poly_div(cAbar, cA, cG);
    nmod_poly_div(cBbar, cB, cG);

    nmod_poly_init(gamma, ctx->ffinfo->mod.n);
    nmod_poly_gcd(gamma, nmod_mpolyun_leadcoeff_ref(A, ctx),
                         nmod_mpolyun_leadcoeff_ref(B, ctx));

    ldegA = nmod_mpolyun_lastdeg(A, ctx);
    ldegB = nmod_mpolyun_lastdeg(B, ctx);
    deggamma = nmod_poly_degree(gamma);
    bound = 1 + deggamma + FLINT_MAX(ldegA, ldegB);

    nmod_poly_init(Aevalp, ctx->ffinfo->mod.n);
    nmod_poly_init(Bevalp, ctx->ffinfo->mod.n);
    nmod_poly_init(Gevalp, ctx->ffinfo->mod.n);
    nmod_poly_init(Abarevalp, ctx->ffinfo->mod.n);
    nmod_poly_init(Bbarevalp, ctx->ffinfo->mod.n);
    nmod_poly_init(Aevalm, ctx->ffinfo->mod.n);
    nmod_poly_init(Bevalm, ctx->ffinfo->mod.n);
    nmod_poly_init(Gevalm, ctx->ffinfo->mod.n);
    nmod_poly_init(Abarevalm, ctx->ffinfo->mod.n);
    nmod_poly_init(Bbarevalm, ctx->ffinfo->mod.n);

    nmod_mpolyun_init(T, A->bits, ctx);

    nmod_poly_init(r, ctx->ffinfo->mod.n);
    nmod_poly_init(alphapow, ctx->ffinfo->mod.n);
    nmod_poly_fit_length(alphapow, FLINT_MAX(WORD(3), bound + 1));
    nmod_poly_init(modulus, ctx->ffinfo->mod.n);
    nmod_poly_init(modulus2, ctx->ffinfo->mod.n);
    nmod_poly_one(modulus);

    if ((ctx->ffinfo->mod.n & UWORD(1)) == UWORD(0))
    {
        success = 0;
        goto cleanup;
    }

    use_stab = 1;
    gstab = bstab = astab = 0;

    alpha = (ctx->ffinfo->mod.n - UWORD(1))/UWORD(2);

choose_prime: /* primes are v - alpha, v + alpha */

    if (alpha < 2)
    {
        success = 0;
        goto cleanup;
    }

    alpha--;

    FLINT_ASSERT(0 < alpha && alpha <= ctx->ffinfo->mod.n/2);
    FLINT_ASSERT(alphapow->alloc >= 2);
    alphapow->length = 2;
    alphapow->coeffs[0] = 1;
    alphapow->coeffs[1] = alpha;

    /* make sure evaluation point does not kill both lc(A) and lc(B) */
    _nmod_poly_eval2_pow(&gammaevalp, &gammaevalm, gamma, alphapow, ctx->ffinfo);
    if (gammaevalp == 0 || gammaevalm == 0)
    {
        goto choose_prime;
    }

    /* evaluation point should kill neither A nor B */
    nmod_mpolyun_eval2_last_bivar(Aevalp, Aevalm, A, alphapow, ctx);
    nmod_mpolyun_eval2_last_bivar(Bevalp, Bevalm, B, alphapow, ctx);
    FLINT_ASSERT(Aevalp->length > 0);
    FLINT_ASSERT(Aevalm->length > 0);
    FLINT_ASSERT(Bevalp->length > 0);
    FLINT_ASSERT(Bevalm->length > 0);

    if (use_stab && gstab)
    {
        slong Gdeg;
        nmod_mpolyun_eval2_last_bivar(Gevalp, Gevalm, G, alphapow, ctx);
        Gdeg = G->exps[0];
        success = 1;
        success = success && nmod_poly_degree(Gevalp) == Gdeg;
        success = success && nmod_poly_degree(Gevalm) == Gdeg;
        success = success && Gevalp->coeffs[Gdeg] == gammaevalp;
        success = success && Gevalm->coeffs[Gdeg] == gammaevalm;
        nmod_poly_divrem_basecase(Abarevalp, r, Aevalp, Gevalp);
        success = success && (r->length == 0);
        nmod_poly_divrem_basecase(Abarevalm, r, Aevalm, Gevalm);
        success = success && (r->length == 0);
        nmod_poly_divrem_basecase(Bbarevalp, r, Bevalp, Gevalp);
        success = success && (r->length == 0);
        nmod_poly_divrem_basecase(Bbarevalm, r, Bevalm, Gevalm);
        success = success && (r->length == 0);

        if (!success)
        {
            use_stab = 0;
            nmod_poly_one(modulus);
            alpha = (ctx->ffinfo->mod.n - UWORD(1))/UWORD(2);
            goto choose_prime;
        }

        nmod_poly_scalar_mul_nmod(Abarevalp, Abarevalp, gammaevalp);
        nmod_poly_scalar_mul_nmod(Abarevalm, Abarevalm, gammaevalm);
        nmod_poly_scalar_mul_nmod(Bbarevalp, Bbarevalp, gammaevalp);
        nmod_poly_scalar_mul_nmod(Bbarevalm, Bbarevalm, gammaevalm);
    }
    else
    {
        nmod_poly_gcd(Gevalp, Aevalp, Bevalp);
        nmod_poly_div(Abarevalp, Aevalp, Gevalp);
        nmod_poly_div(Bbarevalp, Bevalp, Gevalp);
        nmod_poly_gcd(Gevalm, Aevalm, Bevalm);
        nmod_poly_div(Abarevalm, Aevalm, Gevalm);
        nmod_poly_div(Bbarevalm, Bevalm, Gevalm);
        gstab = astab = bstab = 0;
    }

    FLINT_ASSERT(Gevalp->length > 0);
    FLINT_ASSERT(Abarevalp->length > 0);
    FLINT_ASSERT(Bbarevalp->length > 0);
    FLINT_ASSERT(Gevalm->length > 0);
    FLINT_ASSERT(Abarevalm->length > 0);
    FLINT_ASSERT(Bbarevalm->length > 0);

    if (nmod_poly_degree(Gevalp) == 0 || nmod_poly_degree(Gevalm) == 0)
    {
        nmod_mpolyun_one(G, ctx);
        nmod_mpolyun_swap(Abar, A);
        nmod_mpolyun_swap(Bbar, B);
        goto successful_put_content;    
    }

    if (nmod_poly_degree(Gevalp) != nmod_poly_degree(Gevalm))
    {
        goto choose_prime;
    }

    /* the Geval have matching degrees */
    if (nmod_poly_degree(modulus) > 0)
    {
        FLINT_ASSERT(G->length > 0);
        if (nmod_poly_degree(Gevalp) > G->exps[0])
        {
            goto choose_prime;
        }
        else if (nmod_poly_degree(Gevalp) < G->exps[0])
        {
            nmod_poly_one(modulus);
        }
    }

    /* update interpolants */
    nmod_poly_scalar_mul_nmod(Gevalp, Gevalp, gammaevalp);
    nmod_poly_scalar_mul_nmod(Gevalm, Gevalm, gammaevalm);
    if (nmod_poly_degree(modulus) > 0)
    {
        temp = nmod_poly_evaluate_nmod(modulus, alpha);
        FLINT_ASSERT(temp == nmod_poly_evaluate_nmod(modulus, ctx->ffinfo->mod.n - alpha));
        temp = nmod_mul(temp, alpha, ctx->ffinfo->mod);
        temp = nmod_add(temp, temp, ctx->ffinfo->mod);
        nmod_poly_scalar_mul_nmod(modulus, modulus, n_invmod(temp, ctx->ffinfo->mod.n));
        if (!gstab)
        {
            gstab = !nmod_mpolyun_addinterp2_bivar(&ldegG, G, T, Gevalp, Gevalm, modulus, alphapow, ctx);
        }
        nmod_mpolyun_addinterp2_bivar(&ldegAbar, Abar, T, Abarevalp, Abarevalm, modulus, alphapow, ctx);
        nmod_mpolyun_addinterp2_bivar(&ldegBbar, Bbar, T, Bbarevalp, Bbarevalm, modulus, alphapow, ctx);
    }
    else
    {
        nmod_mpolyun_startinterp2_bivar(&ldegG, G, Gevalp, Gevalm, alpha, ctx);
        nmod_mpolyun_startinterp2_bivar(&ldegAbar, Abar, Abarevalp, Abarevalm, alpha, ctx);
        nmod_mpolyun_startinterp2_bivar(&ldegBbar, Bbar, Bbarevalp, Bbarevalm, alpha, ctx);
        gstab = astab = bstab = 0;
    }
    temp = nmod_mul(alpha, alpha, ctx->ffinfo->mod);
    nmod_poly_scalar_mul_nmod(modulus2, modulus, temp);
    nmod_poly_shift_left(modulus, modulus, 2);
    nmod_poly_sub(modulus, modulus, modulus2);

    if (nmod_poly_degree(modulus) < bound)
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

    nmod_poly_one(modulus);
    goto choose_prime;

successful:

    nmod_mpolyun_content_last(modulus, G, ctx);
    nmod_mpolyun_divexact_last(G, modulus, ctx);
    nmod_mpolyun_divexact_last(Abar, nmod_mpolyun_leadcoeff_ref(G, ctx), ctx);
    nmod_mpolyun_divexact_last(Bbar, nmod_mpolyun_leadcoeff_ref(G, ctx), ctx);

successful_put_content:

    nmod_mpolyun_mul_last(G, cG, ctx);
    nmod_mpolyun_mul_last(Abar, cAbar, ctx);
    nmod_mpolyun_mul_last(Bbar, cBbar, ctx);

    success = 1;

cleanup:

#if WANT_ASSERT
/*
    if (success)
    {
        FLINT_ASSERT(1 == nmod_mpolyun_leadcoeff_last(G, ctx));
        nmod_poly_mul(modulus, nmod_mpolyun_leadcoeff_ref(G, ctx),
                               nmod_mpolyun_leadcoeff_ref(Abar, ctx));
        FLINT_ASSERT(nmod_poly_equal(modulus, leadA));
        nmod_poly_mul(modulus, nmod_mpolyun_leadcoeff_ref(G, ctx),
                               nmod_mpolyun_leadcoeff_ref(Bbar, ctx));
        FLINT_ASSERT(nmod_poly_equal(modulus, leadB));
    }
    nmod_poly_clear(leadA);
    nmod_poly_clear(leadB);
*/
#endif

    nmod_poly_clear(cA);
    nmod_poly_clear(cB);
    nmod_poly_clear(cG);
    nmod_poly_clear(cAbar);
    nmod_poly_clear(cBbar);

    nmod_poly_clear(gamma);

    nmod_poly_clear(Aevalp);
    nmod_poly_clear(Bevalp);
    nmod_poly_clear(Gevalp);
    nmod_poly_clear(Abarevalp);
    nmod_poly_clear(Bbarevalp);
    nmod_poly_clear(Aevalm);
    nmod_poly_clear(Bevalm);
    nmod_poly_clear(Gevalm);
    nmod_poly_clear(Abarevalm);
    nmod_poly_clear(Bbarevalm);

    nmod_mpolyun_clear(T, ctx);

    nmod_poly_clear(r);
    nmod_poly_clear(alphapow);
    nmod_poly_clear(modulus);
    nmod_poly_clear(modulus2);

    return success;
}






/************************************************

    discrete logs

***********************************************/


typedef struct {
    mp_limb_t gammapow;
    ulong cm;
} nmod_dlog_table_entry_struct;

static int nmod_dlog_table_entry_struct_cmp(const nmod_dlog_table_entry_struct * lhs,
                                       const nmod_dlog_table_entry_struct * rhs)
{
    return (lhs->gammapow < rhs->gammapow) ? -1 : (lhs->gammapow > rhs->gammapow);
}

typedef struct {
    slong exp;
    ulong prime;
    mp_limb_t gamma;
    mp_limb_t gammainv;
    mp_limb_t startingbeta;
    ulong co;
    ulong startinge;
    ulong idem;
    ulong cbound;
    ulong dbound;
    nmod_dlog_table_entry_struct * table; /* length cbound */
} nmod_dlog_entry_struct;

typedef struct {
    nmod_t mod;         /* p is mod.n */
    mp_limb_t alpha;    /* p.r. of p */
    mp_limb_t alphainv;
    slong num_factors;  /* factors of p - 1*/
    nmod_dlog_entry_struct * entries;
} nmod_dlog_env_struct;
typedef nmod_dlog_env_struct nmod_dlog_env_t[1];

/* assume that p is prime, don't check */
void nmod_dlog_env_init_prime(nmod_dlog_env_t L, mp_limb_t p)
{
    slong i;
    ulong c;
    nmod_dlog_entry_struct * Li;
    n_factor_t factors;

    n_factor_init(&factors);
    n_factor(&factors, p - 1, 1);

    nmod_init(&L->mod, p);
    L->entries = NULL;
    L->num_factors = factors.num;
    if (L->num_factors > 0)
    {
        L->entries = (nmod_dlog_entry_struct*) flint_malloc(L->num_factors*sizeof(nmod_dlog_entry_struct));
    }

    for (i = 0; i < L->num_factors; i++)
    {
        int success;
        fmpz_t pipow, pm1, temp, recp;

        Li = L->entries + i;

        Li->exp = factors.exp[i];
        Li->prime = factors.p[i];

        fmpz_init(recp);
        fmpz_init(temp);
        fmpz_init_set_ui(pipow, Li->prime);
        fmpz_pow_ui(pipow, pipow, Li->exp);
        fmpz_init_set_ui(pm1, p - 1);
        fmpz_divexact(recp, pm1, pipow);
        success = fmpz_invmod(temp, recp, pipow);
        FLINT_ASSERT(success);
        fmpz_mul(temp, temp, recp);

        Li->idem = fmpz_fdiv_ui(temp, p - 1);

        Li->co = fmpz_get_ui(recp);
        Li->startinge = fmpz_get_ui(pipow)/Li->prime;


        fmpz_clear(pipow);
        fmpz_clear(pm1);
        fmpz_clear(temp);
        fmpz_clear(recp);
    }

    L->alpha = 0;
try_alpha:
    L->alpha++;
    if (L->alpha >= p)
    {
        flint_printf("Exception in nmod_dlog_env_init: could not find primitive root\n");
        flint_abort();
    }
    for (i = 0; i < L->num_factors; i++)
    {
        Li = L->entries + i;
        Li->gamma = nmod_pow_ui(L->alpha, (p - 1) / Li->prime, L->mod);
        if (Li->gamma == 1)
        {
            goto try_alpha;
        }
    }

    L->alphainv = nmod_inv(L->alpha, L->mod);

    for (i = 0; i < L->num_factors; i++)
    {
        Li = L->entries + i;
        Li->gammainv = nmod_inv(Li->gamma, L->mod);

        Li->startingbeta = nmod_pow_ui(L->alphainv, Li->co, L->mod);

        Li->dbound = ceil(sqrt((double) Li->prime));
        Li->cbound = (Li->prime + Li->dbound - 1)/Li->dbound;
        while (Li->cbound > 100)
        {
            Li->dbound *= 2;
            Li->cbound = (Li->prime + Li->dbound - 1)/Li->dbound;
        }

        FLINT_ASSERT(Li->dbound > 0);
        FLINT_ASSERT(Li->cbound > 0);
        Li->table = (nmod_dlog_table_entry_struct *) flint_malloc(Li->cbound*sizeof(nmod_dlog_table_entry_struct));

        for (c = 0; c < Li->cbound; c++)
        {
            Li->table[c].cm = c*Li->dbound;
            Li->table[c].gammapow = nmod_pow_ui(Li->gamma, Li->table[c].cm, L->mod);
        }
        qsort(Li->table, Li->cbound, sizeof(nmod_dlog_table_entry_struct),
               (int(*)(const void*, const void*)) nmod_dlog_table_entry_struct_cmp);
        for (c = 1; c < Li->cbound; c++)
        {
            FLINT_ASSERT(Li->table[c - 1].gammapow < Li->table[c].gammapow);
            FLINT_ASSERT(Li->table[c].gammapow == nmod_pow_ui(Li->gamma, Li->table[c].cm, L->mod));
        }
    }    
}


void nmod_dlog_env_init(nmod_dlog_env_t L)
{
    L->num_factors = 0;
    L->entries = NULL;
    nmod_init(&L->mod, 2);
}

void nmod_dlog_env_clear(nmod_dlog_env_t L)
{
    slong i;
    nmod_dlog_entry_struct * Li;

    for (i = 0; i < L->num_factors; i++)
    {
        Li = L->entries + i;
        flint_free(Li->table);
    }

    if (L->entries)
    {
        flint_free(L->entries);
    }
}

/* return x such that x = alpha^y mod p, alpha is the p.r. L->alpha*/
ulong nmod_dlog_env_run(const nmod_dlog_env_t L, mp_limb_t y)
{
    slong i, j;
    ulong x, q, r, e, x0 = 0, x1 = 0, x2 = 0, pp0, pp1, acc, g, pipow;
    ulong lo, mid, hi, d;
    mp_limb_t beta, z, w;
    nmod_dlog_entry_struct * Li;

    FLINT_ASSERT(y != 0);
    FLINT_ASSERT(y < L->mod.n);

    i = 0;
    if (i < L->num_factors && L->entries[i].prime == 2)
    {
        Li = L->entries + i;
        FLINT_ASSERT(Li->prime == 2);

        z = nmod_pow_ui(y, Li->co, L->mod);
        beta = Li->startingbeta;
        e = Li->startinge;
        j = 0;
        pipow = 1; /* Li->prime^j */
        acc = 0;
        do {
            w = nmod_pow_ui(z, e, L->mod);
            /* solve Li->gamma ^ g == w mod p */

            if (w == 1)
            {
                g = 0;
            }
            else
            {
                FLINT_ASSERT(w == Li->gamma);
                g = 1;
                z = nmod_mul(z, beta, L->mod);
            }
            beta = nmod_mul(beta, beta, L->mod);
            acc += g*pipow;
            pipow = pipow*2;
            e = e / 2;
        } while (++j < Li->exp);

        umul_ppmm(pp1, pp0, acc, Li->idem);
        add_sssaaaaaa(x2, x1, x0, x2, x1, x0, WORD(0), pp1, pp0);
        i = 1;
    }

    for (; i < L->num_factors; i++)
    {
        Li = L->entries + i;

        z = nmod_pow_ui(y, Li->co, L->mod);
        beta = Li->startingbeta;
        e = Li->startinge;
        j = 0;
        pipow = 1; /* Li->prime^j */
        acc = 0;
        do {
            w = nmod_pow_ui(z, e, L->mod);
            /* solve Li->gamma ^ g == w mod p */
            d = 0;
            while (1)
            {
                lo = 0; hi = Li->cbound;
                while (hi - lo > 4)
                {
                    mid = lo + (hi - lo)/2;
                    if (Li->table[mid].gammapow == w)
                    {
                        g = Li->table[mid].cm + d;
                        goto found_g;
                    }
                    if (Li->table[mid].gammapow > w)
                        hi = mid;
                    else
                        lo = mid;
                }
                while (lo < hi)
                {
                    if (Li->table[lo].gammapow == w)
                    {
                        g = Li->table[lo].cm + d;
                        goto found_g;
                    }
                    lo++;
                }
                w = nmod_mul(w, Li->gammainv, L->mod);
                d++;
                /* should have found a solution if d is out of bounds */
                FLINT_ASSERT(d < Li->dbound);
            }
        found_g:
            FLINT_ASSERT(g < Li->prime);
            z = nmod_mul(z, nmod_pow_ui(beta, g, L->mod), L->mod);
            beta = nmod_pow_ui(beta, Li->prime, L->mod);
            acc += g*pipow;
            pipow = pipow*Li->prime;
            e = e / Li->prime;
        } while (++j < Li->exp);

        umul_ppmm(pp1, pp0, acc, Li->idem);
        add_sssaaaaaa(x2, x1, x0, x2, x1, x0, WORD(0), pp1, pp0);
    }

    udiv_qrnnd(q, r, x2, x1, L->mod.n - 1);
    udiv_qrnnd(q, x, r, x0, L->mod.n - 1);
    return x;
}




/*
    Assumption on fmpz_mod_dlog_env:
        p is prime.
        The prime factors of p - 1 all fit a ulong.
    The assumption p is prime can be removed, but then phi(p) needs to be calculated by someone somewhere.
    If the prime factors of p - 1 do not fit a ulong you do not want to calculate dlog mod p.

    so we have
    p - 1 = p1^e1 * ... * pn^en for ulong pi and ei
*/

typedef struct {
    fmpz_t gammapow;
    ulong cm;
} fmpz_mod_dlog_table_entry_struct;

static int fmpz_mod_dlog_table_entry_struct_cmp(const fmpz_mod_dlog_table_entry_struct * lhs,
                                       const fmpz_mod_dlog_table_entry_struct * rhs)
{
    return fmpz_cmp(lhs->gammapow, rhs->gammapow);
}


typedef struct {
    slong exp;
    ulong prime;
    fmpz_t gamma;
    fmpz_t gammainv;
    fmpz_t startingbeta;
    fmpz_t co;
    fmpz_t startinge;
    fmpz_t idem;
    ulong cbound;
    ulong dbound;
    fmpz_mod_dlog_table_entry_struct * table; /* length cbound */
} fmpz_mod_dlog_entry_struct;

typedef struct {
    fmpz_mod_ctx_t fpctx;
    fmpz_t pm1;      /* p - 1 */
    fmpz_t alpha;    /* p.r. of p */
    fmpz_t alphainv;
    slong num_factors;  /* factors of p - 1*/
    fmpz_mod_dlog_entry_struct * entries;
} fmpz_mod_dlog_env_struct;
typedef fmpz_mod_dlog_env_struct fmpz_mod_dlog_env_t[1];


/* assume that p is prime, don't check */
void fmpz_mod_dlog_env_init_prime(fmpz_mod_dlog_env_t L, const fmpz_t p)
{
    slong i;
    ulong c;
    fmpz_mod_dlog_entry_struct * Li;
    fmpz_factor_t factors;
    fmpz_t temp;

    fmpz_init(L->alpha);
    fmpz_init(L->alphainv);
    fmpz_init(L->pm1);
    fmpz_mod_ctx_init(L->fpctx, p);

    fmpz_init(temp);

    fmpz_factor_init(factors);
    fmpz_sub_ui(L->pm1, p, 1);
    fmpz_factor(factors, L->pm1);
    L->num_factors = factors->num;
    L->entries = NULL;
    if (L->num_factors > 0)
    {
        L->entries = (fmpz_mod_dlog_entry_struct*) flint_malloc(L->num_factors*sizeof(fmpz_mod_dlog_entry_struct));
    }
    for (i = 0; i < L->num_factors; i++)
    {
        int success;
        fmpz_t pipow, recp;

        Li = L->entries + i;

        fmpz_init(Li->idem);
        fmpz_init(Li->co);
        fmpz_init(Li->startinge);
        fmpz_init(Li->startingbeta);
        fmpz_init(Li->gamma);
        fmpz_init(Li->gammainv);

        FLINT_ASSERT(fmpz_abs_fits_ui(factors->p + i));
        Li->exp = factors->exp[i];
        Li->prime = fmpz_get_ui(factors->p + i);
/*
flint_printf("L[%wd].exp  : %wu\n", i, Li->exp);
flint_printf("L[%wd].prime: %wu\n", i, Li->prime);
*/
        fmpz_init(recp);
        fmpz_init_set_ui(pipow, Li->prime);
        fmpz_pow_ui(pipow, pipow, Li->exp);
        fmpz_divexact(recp, L->pm1, pipow);
        success = fmpz_invmod(temp, recp, pipow);
        FLINT_ASSERT(success);
        fmpz_mul(temp, temp, recp);

        fmpz_mod(Li->idem, temp, L->pm1);

        fmpz_set(Li->co, recp);
        fmpz_divexact_ui(Li->startinge, pipow, Li->prime);

        fmpz_clear(pipow);
        fmpz_clear(recp);
/*
flint_printf("L[%wd].idem: ", i); fmpz_print(Li->idem); printf("\n");
flint_printf("L[%wd].co  : ", i); fmpz_print(Li->co); printf("\n");
flint_printf("L[%wd].stae: ", i); fmpz_print(Li->startinge); printf("\n");
*/
    }
    fmpz_factor_clear(factors);

    fmpz_one(L->alpha);
try_alpha:
    fmpz_add_ui(L->alpha, L->alpha, 1);
/*
flint_printf("\nL.alpha   : ", i); fmpz_print(L->alpha); printf("\n");
*/

    if (fmpz_cmp(L->alpha, p) >= 0)
    {
        flint_printf("Exception in fmpz_mod_dlog_env_init: could not find primitive root\n");
        flint_abort();
    }
    for (i = 0; i < L->num_factors; i++)
    {
        Li = L->entries + i;
        fmpz_divexact_ui(temp, L->pm1, Li->prime);
/*
flint_printf("temp[%wd] : ", i); fmpz_print(temp); printf("\n");
*/

        fmpz_mod_pow_fmpz(Li->gamma, L->alpha, temp, L->fpctx);

/*
flint_printf("gamma[%wd]: ", i); fmpz_print(Li->gamma); printf("\n");
*/

        if (fmpz_is_one(Li->gamma))
        {
            goto try_alpha;
        }
    }

printf("prime "); fmpz_print(p); printf(" alpha "); fmpz_print(L->alpha); printf("\n");

    fmpz_mod_inv(L->alphainv, L->alpha, L->fpctx);
/*
flint_printf("L.alphainv: ", i); fmpz_print(L->alphainv); printf("\n");
*/
    for (i = 0; i < L->num_factors; i++)
    {
        Li = L->entries + i;
/*
flint_printf("L[%wd].p      : ", i); fmpz_print(L->fpctx->p); printf("\n");
flint_printf("L[%wd].gamma  : ", i); fmpz_print(Li->gamma); printf("\n");
*/
        fmpz_mod_inv(Li->gammainv, Li->gamma, L->fpctx);
/*
flint_printf("L[%wd].invgmma: ", i); fmpz_print(Li->gammainv); printf("\n");
*/
        fmpz_mod_pow_fmpz(Li->startingbeta, L->alphainv, Li->co, L->fpctx);
/*
flint_printf("L[%wd].stabeta: ", i); fmpz_print(Li->startinge); printf("\n");
*/

        Li->dbound = ceil(sqrt((double) Li->prime));
        Li->cbound = (Li->prime + Li->dbound - 1)/Li->dbound;
        while (Li->cbound > 100)
        {
            Li->dbound *= 2;
            Li->cbound = (Li->prime + Li->dbound - 1)/Li->dbound;
        }

        FLINT_ASSERT(Li->dbound > 0);
        FLINT_ASSERT(Li->cbound > 0);
        Li->table = (fmpz_mod_dlog_table_entry_struct *) flint_malloc(Li->cbound*sizeof(fmpz_mod_dlog_table_entry_struct));

        for (c = 0; c < Li->cbound; c++)
        {
            Li->table[c].cm = c*Li->dbound;
            fmpz_init(Li->table[c].gammapow);
            fmpz_mod_pow_ui(Li->table[c].gammapow, Li->gamma, Li->table[c].cm, L->fpctx);
        }
        qsort(Li->table, Li->cbound, sizeof(fmpz_mod_dlog_table_entry_struct),
               (int(*)(const void*, const void*)) fmpz_mod_dlog_table_entry_struct_cmp);
        for (c = 1; c < Li->cbound; c++)
        {
            FLINT_ASSERT(fmpz_cmp(Li->table[c - 1].gammapow, Li->table[c].gammapow) < 0);
        }
    }

    fmpz_clear(temp);
}

void fmpz_mod_dlog_env_init(fmpz_mod_dlog_env_t L)
{
    fmpz_t two;
    fmpz_init_set_ui(two, 2);

    L->num_factors = 0;
    L->entries = NULL;
    fmpz_init(L->alpha);
    fmpz_init(L->alphainv);
    fmpz_init(L->pm1);
    fmpz_mod_ctx_init(L->fpctx, two);

    fmpz_clear(two);
    return;
}

void fmpz_mod_dlog_env_clear(fmpz_mod_dlog_env_t L)
{
    slong i;
    ulong c;
    fmpz_mod_dlog_entry_struct * Li;

    for (i = 0; i < L->num_factors; i++)
    {
        Li = L->entries + i;
        fmpz_clear(Li->idem);
        fmpz_clear(Li->co);
        fmpz_clear(Li->startinge);
        fmpz_clear(Li->startingbeta);
        fmpz_clear(Li->gamma);
        fmpz_clear(Li->gammainv);
        for (c = 0; c < Li->cbound; c++)
        {
            fmpz_clear(Li->table[c].gammapow);
        }
        flint_free(Li->table);
    }

    if (L->entries)
    {
        flint_free(L->entries);
    }

    fmpz_clear(L->alpha);
    fmpz_clear(L->alphainv);
    fmpz_clear(L->pm1);

    fmpz_mod_ctx_clear(L->fpctx);

    return;
}

/* return x such that x = alpha^y mod p, alpha is the p.r. L->alpha*/
void fmpz_mod_dlog_env_run(const fmpz_mod_dlog_env_t L, fmpz_t xx, const fmpz_t y)
{
    slong i, j;
/*    ulong x, q, r, e, x0 = 0, x1 = 0, x2 = 0, pp0, pp1, acc, g, pipow;*/
    ulong g;
    fmpz_t x;
    fmpz_t pipow, e, acc;
    ulong lo, mid, hi, d;
    fmpz_t beta, z, w, temp;
    fmpz_mod_dlog_entry_struct * Li;

    fmpz_init(x);
    fmpz_init(acc);
    fmpz_init(pipow);
    fmpz_init(e);
    fmpz_init(beta);
    fmpz_init(z);
    fmpz_init(w);
    fmpz_init(temp);

    FLINT_ASSERT(!fmpz_is_zero(y));
    FLINT_ASSERT(fmpz_mod_is_canonical(y, L->fpctx));

    i = 0;
    if (i < L->num_factors && L->entries[i].prime == 2)
    {
        Li = L->entries + i;
        FLINT_ASSERT(Li->prime == 2);

        fmpz_mod_pow_fmpz(z, y, Li->co, L->fpctx);
        fmpz_set(beta, Li->startingbeta);
        fmpz_set(e, Li->startinge);
        j = 0;
        fmpz_one(pipow); /* Li->prime^j */
        fmpz_zero(acc);
        do {
            fmpz_mod_pow_fmpz(w, z, e, L->fpctx);
            /* solve Li->gamma ^ g == w mod p */
            if (fmpz_is_one(w))
            {
                g = 0;
            }
            else
            {
                if(!fmpz_equal(w, Li->gamma))
                {
                    flint_printf("Exception in fmpz_mod_dlog_env_run: could not find log\n");
                    flint_abort();
                }
                g = 1;
                fmpz_mod_mul(z, z, beta, L->fpctx);
            }
            fmpz_mod_mul(beta, beta, beta, L->fpctx);
            fmpz_addmul_ui(acc, pipow, g);
            fmpz_mul_2exp(pipow, pipow, 1);
            fmpz_tdiv_q_2exp(e, e, 1);
        } while (++j < Li->exp);

        fmpz_addmul(x, acc, Li->idem);
        i = 1;
    }

    for (; i < L->num_factors; i++)
    {
        Li = L->entries + i;

        fmpz_mod_pow_fmpz(z, y, Li->co, L->fpctx);
        fmpz_set(beta, Li->startingbeta);
        fmpz_set(e, Li->startinge);
        j = 0;
        fmpz_one(pipow); /* Li->prime^j */
        fmpz_zero(acc);
        do {
            fmpz_mod_pow_fmpz(w, z, e, L->fpctx);
            /* solve Li->gamma ^ g == w mod p */
            d = 0;
            while (1)
            {
                lo = 0; hi = Li->cbound;
                while (hi - lo > 4)
                {
                    int cmp;
                    mid = lo + (hi - lo)/2;
                    cmp = fmpz_cmp(Li->table[mid].gammapow, w);
                    if (cmp == 0)
                    {
                        g = Li->table[mid].cm + d;
                        goto found_g;
                    }
                    else if (cmp > 0)
                    {
                        hi = mid;
                    }
                    else
                    {
                        lo = mid;
                    }
                }
                while (lo < hi)
                {
                    if (fmpz_equal(Li->table[lo].gammapow, w))
                    {
                        g = Li->table[lo].cm + d;
                        goto found_g;
                    }
                    lo++;
                }
                fmpz_mod_mul(w, w, Li->gammainv, L->fpctx);
                d++;
                if (d >= Li->dbound)
                {
                    flint_printf("Exception in fmpz_mod_dlog_env_run: could not find log\n");
                    flint_abort();
                }
            }
        found_g:
            FLINT_ASSERT(g < Li->prime);
            fmpz_mod_pow_ui(temp, beta, g, L->fpctx);
            fmpz_mod_mul(z, z, temp, L->fpctx);
            fmpz_mod_pow_ui(beta, beta, Li->prime, L->fpctx);
            fmpz_addmul_ui(acc, pipow, g);
            fmpz_mul_ui(pipow, pipow, Li->prime);
            fmpz_divexact_ui(e, e, Li->prime);
        } while (++j < Li->exp);

        fmpz_addmul(x, acc, Li->idem);
    }

    fmpz_mod(xx, x, L->pm1);
    fmpz_clear(acc);
    fmpz_clear(pipow);
    fmpz_clear(e);
    fmpz_clear(beta);
    fmpz_clear(z);
    fmpz_clear(w);
    fmpz_clear(temp);
    fmpz_clear(x);
    return;
}




typedef struct
{
    slong * degbounds;
    ulong * subdegs;
    fmpz_mod_dlog_env_t dlogenv;
    nmod_dlog_env_t dlogenv_sp;
} mpoly_bma_interpolate_ctx_struct;
typedef mpoly_bma_interpolate_ctx_struct mpoly_bma_interpolate_ctx_t[1];

void mpoly_bma_interpolate_ctx_init(
    mpoly_bma_interpolate_ctx_t Ictx,
    slong nvars)
{
    Ictx->degbounds = (slong *) flint_malloc(nvars*sizeof(slong));
    Ictx->subdegs   = (ulong *) flint_malloc(nvars*sizeof(ulong));
    fmpz_mod_dlog_env_init(Ictx->dlogenv);
    nmod_dlog_env_init(Ictx->dlogenv_sp);
}

void mpoly_bma_interpolate_ctx_reset_prime(
    mpoly_bma_interpolate_ctx_t Ictx,
    const fmpz_t p)
{
    fmpz_mod_dlog_env_clear(Ictx->dlogenv);
    fmpz_mod_dlog_env_init_prime(Ictx->dlogenv, p);
}

void mpoly_bma_interpolate_ctx_reset_prime_sp(
    mpoly_bma_interpolate_ctx_t Ictx,
    mp_limb_t p)
{
    nmod_dlog_env_clear(Ictx->dlogenv_sp);
    nmod_dlog_env_init_prime(Ictx->dlogenv_sp, p);
}


void mpoly_bma_interpolate_ctx_clear(mpoly_bma_interpolate_ctx_t Ictx)
{
    flint_free(Ictx->degbounds);
    flint_free(Ictx->subdegs);
    fmpz_mod_dlog_env_clear(Ictx->dlogenv);
    nmod_dlog_env_clear(Ictx->dlogenv_sp);
}









typedef struct {
    slong half_point_count;
    nmod_poly_t R0, R1;
    nmod_poly_t V0, V1; /* V1 is our master polynomial */
    nmod_poly_t r; /* temporary */
    nmod_poly_t q; /* temporary also used as a queue of incoming points */
} nmod_bma_struct;
typedef nmod_bma_struct nmod_bma_t[1];

/*
    A = x^(2*n)       n = half_point_count
    deg(S) < 2*n
    U0*A + V0*S = R0   deg(R0) >= n
    U1*A + V1*S = R1   deg(R1) < n

    S is the reverse of the polynomial whose coefficients are the input points.
        S = a_0*x^(2n-1) + a_1*x^(2n-2) + ... + a_(2n-1)
    S can be updated with more points at any time via add_points.
    S can be reduced (compute V0, V1, R0, R1) at any time via reduce, which
        returns whether the reduction caused a change.

    The U0 and U1 are not stored.
*/
void nmod_bma_init(nmod_bma_t B, mp_limb_t p)
{
    B->half_point_count = 0;
    nmod_poly_init(B->V0, p);
    nmod_poly_init(B->R0, p);
    nmod_poly_one(B->R0);
    nmod_poly_init(B->V1, p);
    nmod_poly_one(B->V1);
    nmod_poly_init(B->R1, p);
    nmod_poly_init(B->r, p);
    nmod_poly_init(B->q, p);
    B->q->length = 0;
}


void nmod_bma_start_over(nmod_bma_t B)
{
    B->half_point_count = 0;
    nmod_poly_zero(B->V0);
    nmod_poly_one(B->R0);
    nmod_poly_one(B->V1);
    nmod_poly_zero(B->R1);
    B->q->length = 0;
}

void nmod_bma_clear(nmod_bma_t B)
{
    nmod_poly_clear(B->R0);
    nmod_poly_clear(B->R1);
    nmod_poly_clear(B->V0);
    nmod_poly_clear(B->V1);
    nmod_poly_clear(B->r);
    nmod_poly_clear(B->q);
}

void nmod_bma_reset_prime(nmod_bma_t B, mp_limb_t p)
{
    nmod_bma_clear(B);
    nmod_bma_init(B, p);
}

void nmod_bma_add_points(nmod_bma_t B, mp_limb_t * a, slong count)
{
    slong i;
    slong old_length = B->q->length;
    nmod_poly_fit_length(B->q, old_length + count);
    for (i = 0; i < count; i++)
    {
        B->q->coeffs[old_length + i] = a[i];
    }
    B->q->length = old_length + count;
}

void nmod_bma_add_point(nmod_bma_t B, mp_limb_t a)
{
    slong old_length = B->q->length;
    nmod_poly_fit_length(B->q, old_length + 1);
    B->q->coeffs[old_length] = a;
    B->q->length = old_length + 1;
}

int nmod_bma_reduce(nmod_bma_t B)
{
    int changed = 0;
    slong i, queue_length = B->q->length;

    if ((queue_length % 2) != 0)
    {
        flint_printf("Exception in nmod_bma_reduce: point count is not even\n");
        flint_abort();
    }

    /* reverse the queue into temp r */
    B->half_point_count += queue_length/2;
    nmod_poly_zero(B->r);
    for (i = 0; i < queue_length; i++)
    {
        nmod_poly_set_coeff_ui(B->r, queue_length - i - 1, B->q->coeffs[i]);
    }
    nmod_poly_mul(B->q, B->V0, B->r);
    nmod_poly_shift_left(B->R0, B->R0, queue_length);
    nmod_poly_add(B->R0, B->R0, B->q);
    nmod_poly_mul(B->q, B->V1, B->r);
    nmod_poly_shift_left(B->R1, B->R1, queue_length);
    nmod_poly_add(B->R1, B->R1, B->q);

    /* now reduce */
    while (B->half_point_count < nmod_poly_length(B->R1))
    {
        changed = 1;
        nmod_poly_divrem(B->q, B->r, B->R0, B->R1);
        nmod_poly_swap(B->R0, B->R1);
        nmod_poly_swap(B->R1, B->r);

        nmod_poly_mul(B->r, B->q, B->V1);
        nmod_poly_sub(B->q, B->V0, B->r);
        nmod_poly_swap(B->V0, B->V1);
        nmod_poly_swap(B->V1, B->q);
        FLINT_ASSERT(nmod_poly_degree(B->V1) > nmod_poly_degree(B->V0));
    }

    /* queue is empty now */
    B->q->length = 0;
    return changed;
}


/* split f assuming that f has degree(f) distinct nonzero roots in Fp */
static void _nmod_poly_rabinsplit(nmod_poly_t a, nmod_poly_t b, nmod_poly_t T, 
                                  const nmod_poly_t f, flint_rand_t randstate)
{
    mp_limb_t delta;

    FLINT_ASSERT(nmod_poly_degree(f) > 1);

try_again:

    delta = n_randint(randstate, f->mod.n);

    nmod_poly_zero(a);
    nmod_poly_set_coeff_ui(a, 1, 1);
    nmod_poly_set_coeff_ui(a, 0, delta);
    nmod_poly_powmod_ui_binexp(T, a, (f->mod.n - 1)/2, f);
    nmod_poly_zero(a);
    nmod_poly_set_coeff_ui(a, 0, 1);
    nmod_poly_sub(T, T, a);
    nmod_poly_gcd(a, T, f);
    FLINT_ASSERT(!nmod_poly_is_zero(a));
    if (0 >= nmod_poly_degree(a) || nmod_poly_degree(a) >= nmod_poly_degree(f))
    {
        goto try_again;
    }
    nmod_poly_div(b, f, a);
    /* deg a >= deg b */
    if (nmod_poly_degree(a) < nmod_poly_degree(b))
    {
        nmod_poly_swap(a, b);
    }
    return;
}

/* fill in roots with the t distinct nonzero roots of master, or fail */
static int _nmod_find_roots(mp_limb_t * roots, const nmod_poly_t master, slong t)
{
    mp_limb_t a0, a1;
    int success;
    slong i, roots_idx, sp;
    mp_limb_t delta;
    nmod_poly_struct * a , * b;
    nmod_poly_t f, T;
    nmod_poly_struct stack[FLINT_BITS + 1];
    flint_rand_t randstate;

    FLINT_ASSERT(t >= 0);
    FLINT_ASSERT(t == nmod_poly_degree(master));

    if (t == 0)
    {
        return 1;
    }
    else if (t == 1)
    {
        a0 = nmod_poly_get_coeff_ui(master, 0);
        a1 = nmod_poly_get_coeff_ui(master, 1);
        if (a0 == 0)
        {
            return 0;
        }
        roots[0] = nmod_mul(a0, nmod_inv(master->mod.n - a1, master->mod), master->mod);
        return 1;
    }

    flint_randinit(randstate);
    nmod_poly_init(T, master->mod.n);
    nmod_poly_init(f, master->mod.n);
    for (i = 0; i <= FLINT_BITS; i++)
    {
        nmod_poly_init(stack + i, master->mod.n);
    }

    roots_idx = 0;

    nmod_poly_make_monic(f, master);

    a = stack + 0;
    delta = n_randint(randstate, master->mod.n);
    nmod_poly_zero(a);
    nmod_poly_set_coeff_ui(a, 1, 1);
    nmod_poly_set_coeff_ui(a, 0, delta);
    nmod_poly_powmod_ui_binexp(T, a, (master->mod.n - 1)/2, f);
    nmod_poly_zero(a);
    nmod_poly_set_coeff_ui(a, 0, 1);
    nmod_poly_sub(T, T, a);
    nmod_poly_gcd(a, T, f);

    b = stack + 1;
    nmod_poly_zero(b);
    nmod_poly_set_coeff_ui(b, 0, 2);
    nmod_poly_add(T, T, b);
    nmod_poly_gcd(b, T, f);

    if (nmod_poly_degree(b) + nmod_poly_degree(a) != t)
    {
        success = 0;
        goto cleanup;
    }
    /* deg a >= deg b */
    if (nmod_poly_degree(a) < nmod_poly_degree(b))
    {
        nmod_poly_swap(a, b);
    }

    sp = nmod_poly_degree(b) > 0 ? 2 : 1;
    while (sp > 0)
    {
        FLINT_ASSERT(sp < FLINT_BITS);
        sp--;
        nmod_poly_swap(f, stack + sp);

        FLINT_ASSERT(nmod_poly_degree(f) > 0);
        if (nmod_poly_degree(f) == 1)
        {
            a0 = nmod_poly_get_coeff_ui(f, 0);
            a1 = nmod_poly_get_coeff_ui(f, 1);
            FLINT_ASSERT(a0 != 0);
            FLINT_ASSERT(a1 == 1);
            roots[roots_idx] = master->mod.n - a0;
            roots_idx++;
        }
        else
        {
            _nmod_poly_rabinsplit(stack + sp + 0, stack + sp + 1, T, f, randstate);
            FLINT_ASSERT(FLINT_BIT_COUNT(nmod_poly_degree(stack + sp + 1)) <= FLINT_BITS - sp - 1);
            sp += 2;
        }
    }

    success = 1;

cleanup:

    flint_randclear(randstate);
    nmod_poly_clear(T);
    nmod_poly_clear(f);
    for (i = 0; i <= FLINT_BITS; i++)
    {
        nmod_poly_clear(stack + i);
    }

    if (success)
    {
        FLINT_ASSERT(roots_idx == t);
    }

    return success;
}


typedef struct
{
    const nmod_mpoly_ctx_struct * polyctx;
    nmod_dlog_env_t dlogenv;
    slong * degbounds;
    ulong * subdegs;
    mp_bitcnt_t bits;
    ulong * inputexpmask;
} nmod_mpoly_bma_interpolate_ctx_struct;
typedef nmod_mpoly_bma_interpolate_ctx_struct nmod_mpoly_bma_interpolate_ctx_t[1];

typedef struct {
    nmod_poly_t roots;
    nmod_poly_t evals;
    nmod_bma_t bma;
} nmod_mpoly_bma_interpolate_struct;
typedef nmod_mpoly_bma_interpolate_struct nmod_mpoly_bma_interpolate_t[1];



void nmod_mpoly_bma_interpolate_ctx_init(nmod_mpoly_bma_interpolate_ctx_t Ictx, mp_bitcnt_t bits_,
                                                const nmod_mpoly_ctx_t pctx)
{
    slong N;

    Ictx->polyctx = pctx;
    Ictx->degbounds = (slong *) flint_malloc(pctx->minfo->nvars*sizeof(slong));
    Ictx->subdegs = (ulong *) flint_malloc(pctx->minfo->nvars*sizeof(ulong));
    nmod_dlog_env_init_prime(Ictx->dlogenv, pctx->ffinfo->mod.n);

    Ictx->bits = bits_;
    N = mpoly_words_per_exp_sp(bits_, pctx->minfo);
    Ictx->inputexpmask = (ulong *) flint_malloc(N*sizeof(ulong));
    mpoly_monomial_zero(Ictx->inputexpmask, N);
}

void nmod_mpoly_bma_interpolate_ctx_clear(nmod_mpoly_bma_interpolate_ctx_t Ictx)
{
    flint_free(Ictx->inputexpmask);
    nmod_dlog_env_clear(Ictx->dlogenv);
    flint_free(Ictx->degbounds);
    flint_free(Ictx->subdegs);
}


void nmod_mpoly_bma_interpolate_init(nmod_mpoly_bma_interpolate_t I, mp_limb_t p)
{
    nmod_poly_init(I->roots, p);
    nmod_poly_init(I->evals, p);
    nmod_bma_init(I->bma, p);
}

slong nmod_mpoly_bma_interpolate_pointcount(const nmod_mpoly_bma_interpolate_t I)
{
    return I->bma->q->length + 2*I->bma->half_point_count;
}

void nmod_mpoly_bma_interpolate_print(const nmod_mpoly_bma_interpolate_t I)
{
    flint_printf("(%wd) ", nmod_mpoly_bma_interpolate_pointcount(I));
    nmod_poly_print_pretty(I->bma->V1, "#");
}

void nmod_mpoly_bma_interpolate_reset_prime(
    nmod_mpoly_bma_interpolate_t I,
    mp_limb_t p)
{
    nmod_poly_clear(I->roots);
    nmod_poly_init(I->roots, p);
    nmod_poly_clear(I->evals);
    nmod_poly_init(I->evals, p);
    nmod_bma_reset_prime(I->bma, p);
}

void nmod_mpoly_bma_interpolate_clear(nmod_mpoly_bma_interpolate_t I)
{
    nmod_poly_clear(I->roots);
    nmod_poly_clear(I->evals);
    nmod_bma_clear(I->bma);
}

void nmod_mpoly_bma_interpolate_add_point(nmod_mpoly_bma_interpolate_t I, mp_limb_t a)
{
    nmod_poly_fit_length(I->evals, I->evals->length + 1);
    I->evals->coeffs[I->evals->length] = a;
    I->evals->length++;

    nmod_bma_add_point(I->bma, a);
}


void nmod_mpoly_bma_interpolate_zero(nmod_mpoly_bma_interpolate_t I, slong count)
{
    slong i;
    nmod_poly_fit_length(I->evals, count);
    nmod_bma_start_over(I->bma);
    for (i = 0; i < count; i++)
    {
        I->evals->coeffs[i] = 0;
        nmod_bma_add_point(I->bma, I->evals->coeffs[i]);
    }
    I->evals->length = count;
}

int nmod_mpoly_bma_interpolate_reduce(nmod_mpoly_bma_interpolate_t I)
{
    return nmod_bma_reduce(I->bma);
}

void nmod_mpoly_bma_interpolate_eval_init(nmod_mpoly_bma_interpolate_ctx_t Ictx, const nmod_mpoly_t A)
{    
    slong i, j, N = mpoly_words_per_exp_sp(Ictx->bits, Ictx->polyctx->minfo);
    FLINT_ASSERT(A->bits == Ictx->bits);

    for (i = 0; i < A->length; i++)
    {
        for (j = 0; j < N; j++)
        {
            Ictx->inputexpmask[j] |= (A->exps + N*i)[j];
        }
    }
}



/*
    put the evaluation of the monomials in A at alpha^w in the coeffs of E

    x_0     = alpha ^ (w * subdegs[n-1] * subdegs[n-2] * ... * * subdegs[1])
      ...
    x_(n-3) = alpha ^ (w * subdegs[n-1] * subdegs[n-2])
    x_(n-2) = alpha ^ (w * subdegs[n-1])
    x_(n-1) = alpha ^ (w)

    secret: subdegs[0] is not used
*/
void nmod_mpoly_bma_interpolate_eval_setskel(nmod_mpoly_t M, const nmod_mpoly_t A, ulong w, const nmod_mpoly_bma_interpolate_ctx_t Ictx)
{
    slong i, j;
    slong offset, shift;
    slong N = mpoly_words_per_exp_sp(A->bits, Ictx->polyctx->minfo);
    slong nvars = Ictx->polyctx->minfo->nvars;
    ulong * Aexp;
    mp_limb_t * Mcoeff;
    slong * LUToffset;
    ulong * LUTmask;
    mp_limb_t * LUTvalue;
    slong LUTlen;
    mp_limb_t xeval, xpoweval;
    TMP_INIT;

    TMP_START;

    LUToffset = (slong *) TMP_ALLOC(N*FLINT_BITS*sizeof(slong));
    LUTmask   = (ulong *) TMP_ALLOC(N*FLINT_BITS*sizeof(ulong));
    LUTvalue  = (mp_limb_t *) TMP_ALLOC(N*FLINT_BITS*sizeof(mp_limb_t));

    FLINT_ASSERT(M->bits == A->bits);

    nmod_mpoly_fit_length(M, A->length, Ictx->polyctx);
    M->length = A->length;

    Mcoeff = M->coeffs;
    Aexp = A->exps;

    LUTlen = 0;
    xeval = nmod_pow_ui(Ictx->dlogenv->alpha, w, Ictx->polyctx->ffinfo->mod);
    for (j = nvars - 1; j >= 0; j--)
    {
        mpoly_gen_offset_shift_sp(&offset, &shift, j, A->bits, Ictx->polyctx->minfo);

        xpoweval = xeval; /* xpoweval = xeval^(2^i) */
        for (i = 0; i < A->bits; i++)
        {
            LUToffset[LUTlen] = offset;
            LUTmask[LUTlen] = (UWORD(1) << (shift + i));
            LUTvalue[LUTlen] = xpoweval;
            if ((Ictx->inputexpmask[offset] & LUTmask[LUTlen]) != 0)
            {
                LUTlen++;
            }
            xpoweval = nmod_mul(xpoweval, xpoweval, Ictx->polyctx->ffinfo->mod);
        }
        xeval = nmod_pow_ui(xeval, Ictx->subdegs[j], Ictx->polyctx->ffinfo->mod);
    }
    FLINT_ASSERT(LUTlen < N*FLINT_BITS);

    for (i = 0; i < A->length; i++)
    {
        mp_limb_t t = 1;
        for (j = 0; j < LUTlen; j++)
        {
            if (((Aexp + N*i)[LUToffset[j]] & LUTmask[j]) != 0)
            {
                t = nmod_mul(t, LUTvalue[j], Ictx->polyctx->ffinfo->mod);
            }
        }
        Mcoeff[i] = t;
        mpoly_monomial_zero(M->exps + N*i, N);
    }

    TMP_END;
}

/* M = S */
void nmod_mpoly_bma_interpolate_eval_copyskel(nmod_mpoly_t M, const nmod_mpoly_t S, const nmod_mpoly_bma_interpolate_ctx_t Ictx)
{
    slong i, N;
    nmod_mpoly_fit_length(M, S->length, Ictx->polyctx);
    M->length = S->length;
    _nmod_vec_set(M->coeffs, S->coeffs, S->length);
    N = mpoly_words_per_exp(Ictx->bits, Ictx->polyctx->minfo);
    for (i = 0; i < M->length; i++)
    {
        mpoly_monomial_zero(M->exps + N*i, N);
    }
}

/* return A.M */
mp_limb_t nmod_mpoly_bma_interpolate_eval_useskel(const nmod_mpoly_t A, const nmod_mpoly_t M, const nmod_mpoly_bma_interpolate_ctx_t Ictx)
{
    slong i;
    mp_limb_t t = 0;

    FLINT_ASSERT(M->length == A->length);
    for (i = 0; i < A->length; i++)
    {
        t = nmod_add(t, nmod_mul(A->coeffs[i], M->coeffs[i], Ictx->polyctx->ffinfo->mod), Ictx->polyctx->ffinfo->mod);
    }
    return t;
}

/* return A.M and multiply M by S */
mp_limb_t nmod_mpoly_bma_interpolate_eval_useskelmul(const nmod_mpoly_t A, nmod_mpoly_t M, const nmod_mpoly_t S, const nmod_mpoly_bma_interpolate_ctx_t Ictx)
{
    slong i;
    mp_limb_t t = 0;

    for (i = 0; i < A->length; i++)
    {
        t = nmod_add(t, nmod_mul(A->coeffs[i], M->coeffs[i], Ictx->polyctx->ffinfo->mod), Ictx->polyctx->ffinfo->mod);
        M->coeffs[i] = nmod_mul(M->coeffs[i], S->coeffs[i], Ictx->polyctx->ffinfo->mod);
    }
    return t;
}



int nmod_mpoly_bma_interpolate_get_mpoly(
    fmpz_mpoly_t A,
    const fmpz_mpoly_ctx_t ctx,
    ulong alphashift,
    nmod_mpoly_bma_interpolate_t I,
    const mpoly_bma_interpolate_ctx_t Ictx,
    const nmodf_ctx_t fpctx)
{
    slong i, j, t, N;
    int success;
    ulong new_exp, this_exp;
    slong * shifts, * offsets;
    mp_limb_t * values, * roots;
    fmpz * Acoeff;
    ulong * Aexp;
    slong Alen;
    mp_limb_t T, S, V, V0, V1, V2, p0, p1, r;
    TMP_INIT;

    TMP_START;

    nmod_bma_reduce(I->bma);
    t = nmod_poly_degree(I->bma->V1);
    FLINT_ASSERT(t > 0);
    FLINT_ASSERT(I->evals->length >= t);

    nmod_poly_fit_length(I->roots, t);
    I->roots->length = t;
    success = _nmod_find_roots(I->roots->coeffs, I->bma->V1, t);
    if (!success)
    {
        goto cleanup;
    }

    roots = I->roots->coeffs;
    values = I->evals->coeffs;

    FLINT_ASSERT(ctx->minfo->ord == ORD_LEX);
    FLINT_ASSERT(A->bits <= FLINT_BITS);
    N = mpoly_words_per_exp_sp(A->bits, ctx->minfo);
    fmpz_mpoly_fit_length(A, t, ctx);
    Acoeff = A->coeffs;
    Aexp = A->exps;
    A->length = 0;

    shifts = (slong *) TMP_ALLOC(ctx->minfo->nvars);
    offsets = (slong *) TMP_ALLOC(ctx->minfo->nvars);
    for (j = 0; j < ctx->minfo->nvars; j++)
    {
        mpoly_gen_offset_shift_sp(offsets + j, shifts + j, j, A->bits, ctx->minfo);
    }

    Alen = 0;
    for (i = 0; i < t; i++)
    {
        /* coeffs[i] is (coeffs(P).values)/P(roots[i]) =: V/S
            where P(x) = master(x)/(x-roots[i])     */
        V0 = V1 = V2 = T = S = 0;
        r = roots[i];
        for (j = t; j > 0; j--)
        {
            T = nmod_add(nmod_mul(r, T, fpctx->mod), I->bma->V1->coeffs[j], fpctx->mod);
            S = nmod_add(nmod_mul(r, S, fpctx->mod), T, fpctx->mod);
            umul_ppmm(p1, p0, values[j - 1], T);
            add_sssaaaaaa(V2, V1, V0, V2, V1, V0, WORD(0), p1, p0);
        }
        /* roots[i] should be a root of master */
        FLINT_ASSERT(nmod_add(nmod_mul(r, T, fpctx->mod), I->bma->V1->coeffs[0], fpctx->mod) == 0);
        NMOD_RED3(V, V2, V1, V0, fpctx->mod);
        S = nmod_mul(S, nmod_pow_ui(r, alphashift, fpctx->mod), fpctx->mod);
        V0 = nmod_mul(V, nmod_inv(S, fpctx->mod), fpctx->mod);
        if (V0 == 0)
        {
            /* hmmm */
            continue;
        }

        fmpz_set_ui(Acoeff + Alen, V0);
        if (fpctx->mod.n - V0 < V0)
        {
            fmpz_sub_ui(Acoeff + Alen, Acoeff + Alen, fpctx->mod.n);
        }

        mpoly_monomial_zero(Aexp + N*Alen, N);
        new_exp = nmod_dlog_env_run(Ictx->dlogenv_sp, roots[i]);
        for (j = ctx->minfo->nvars - 1; j >= 0; j--)
        {
            this_exp = new_exp % Ictx->subdegs[j];
            new_exp = new_exp / Ictx->subdegs[j];
            if (this_exp >= Ictx->degbounds[j])
            {
                success = 0;
                goto cleanup;
            }
            (Aexp + N*Alen)[offsets[j]] |= this_exp << shifts[j];
        }
        if (new_exp != 0)
        {
            success = 0;
            goto cleanup;
        }
        Alen++;
    }
    A->length = Alen;

    fmpz_mpoly_sort_terms(A, ctx);

    success = 1;

cleanup:

    TMP_END;
    return success;
}











/***************************************************************

    start of fmpz !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

******************************************************************/





typedef struct
{
   mp_limb_t * coeffs;
   slong alloc;
   slong length;
} nmod_mpolyc_struct;

typedef nmod_mpolyc_struct nmod_mpolyc_t[1];

void nmod_mpolyc_init(nmod_mpolyc_t A)
{
    A->coeffs = NULL;
    A->alloc = 0;
    A->length = 0;
}

void nmod_mpolyc_clear(nmod_mpolyc_t A)
{
    if (A->coeffs)
        flint_free(A->coeffs);
}

void nmod_mpolyc_fit_length(nmod_mpolyc_t A, slong length)
{
    slong old_alloc = A->alloc;
    slong new_alloc = FLINT_MAX(length, A->alloc + A->alloc/2);

    if (length > old_alloc)
    {
        if (old_alloc == 0)
        {
            A->coeffs = (mp_limb_t *) flint_malloc(new_alloc*sizeof(mp_limb_t));
        }
        else
        {
            A->coeffs = (mp_limb_t *) flint_realloc(A->coeffs, new_alloc*sizeof(mp_limb_t));
        }

        A->alloc = new_alloc;
    }
}

typedef struct
{
   nmod_mpolyc_struct * coeffs;
   slong alloc;
   slong length;
} nmod_mpolycu_struct;

typedef nmod_mpolycu_struct nmod_mpolycu_t[1];

void nmod_mpolycu_init(nmod_mpolycu_t A)
{
    A->coeffs = NULL;
    A->alloc = 0;
    A->length = 0;
}

void nmod_mpolycu_clear(nmod_mpolycu_t A)
{
    slong i;
    for (i = 0; i < A->alloc; i++)
    {
        nmod_mpolyc_clear(A->coeffs + i);
    }
    if (A->coeffs)
        flint_free(A->coeffs);
}

void nmod_mpolycu_fit_length(nmod_mpolycu_t A, slong length)
{
    slong i;
    slong old_alloc = A->alloc;
    slong new_alloc = FLINT_MAX(length, A->alloc + A->alloc/2);

    if (length > old_alloc)
    {
        if (old_alloc == 0)
        {
            A->coeffs = (nmod_mpolyc_struct *) flint_malloc(new_alloc*sizeof(nmod_mpolyc_struct));
        }
        else
        {
            A->coeffs = (nmod_mpolyc_struct *) flint_realloc(A->coeffs, new_alloc*sizeof(nmod_mpolyc_struct));
        }

        for (i = old_alloc; i < new_alloc; i++)
        {
            nmod_mpolyc_init(A->coeffs + i);
        }
        A->alloc = new_alloc;
    }
}


typedef struct
{
   fmpz * coeffs;
   slong alloc;
   slong length;
} fmpz_mpolyc_struct;

typedef fmpz_mpolyc_struct fmpz_mpolyc_t[1];

void fmpz_mpolyc_init(fmpz_mpolyc_t A)
{
    A->coeffs = NULL;
    A->alloc = 0;
    A->length = 0;
}

void fmpz_mpolyc_clear(fmpz_mpolyc_t A)
{
    slong i;
    for (i = 0; i < A->alloc; i++)
    {
        fmpz_clear(A->coeffs + i);
    }
    if (A->coeffs)
        flint_free(A->coeffs);
}

void fmpz_mpolyc_fit_length(fmpz_mpolyc_t A, slong length)
{
    slong i;
    slong old_alloc = A->alloc;
    slong new_alloc = FLINT_MAX(length, A->alloc + A->alloc/2);

    if (length > old_alloc)
    {
        if (old_alloc == 0)
        {
            A->coeffs = (fmpz *) flint_malloc(new_alloc*sizeof(fmpz));
        }
        else
        {
            A->coeffs = (fmpz *) flint_realloc(A->coeffs, new_alloc*sizeof(fmpz));
        }

        for (i = old_alloc; i < new_alloc; i++)
        {
            fmpz_init(A->coeffs + i);
        }
        A->alloc = new_alloc;
    }
}

typedef struct
{
   fmpz_mpolyc_struct * coeffs;
   slong alloc;
   slong length;
} fmpz_mpolycu_struct;

typedef fmpz_mpolycu_struct fmpz_mpolycu_t[1];

void fmpz_mpolycu_init(fmpz_mpolycu_t A)
{
    A->coeffs = NULL;
    A->alloc = 0;
    A->length = 0;
}

void fmpz_mpolycu_clear(fmpz_mpolycu_t A)
{
    slong i;
    for (i = 0; i < A->alloc; i++)
    {
        fmpz_mpolyc_clear(A->coeffs + i);
    }
    if (A->coeffs)
        flint_free(A->coeffs);
}

void fmpz_mpolycu_fit_length(fmpz_mpolycu_t A, slong length)
{
    slong i;
    slong old_alloc = A->alloc;
    slong new_alloc = FLINT_MAX(length, A->alloc + A->alloc/2);

    if (length > old_alloc)
    {
        if (old_alloc == 0)
        {
            A->coeffs = (fmpz_mpolyc_struct *) flint_malloc(new_alloc*sizeof(fmpz_mpolyc_struct));
        }
        else
        {
            A->coeffs = (fmpz_mpolyc_struct *) flint_realloc(A->coeffs, new_alloc*sizeof(fmpz_mpolyc_struct));
        }

        for (i = old_alloc; i < new_alloc; i++)
        {
            fmpz_mpolyc_init(A->coeffs + i);
        }
        A->alloc = new_alloc;
    }
}













typedef struct
{
   fmpz_mod_poly_struct * coeffs;
   ulong * exps;
   slong alloc;
   slong length;
   slong bits;
} fmpz_mod_mpolyn_struct;
typedef fmpz_mod_mpolyn_struct fmpz_mod_mpolyn_t[1];

typedef struct
{
    fmpz_mod_mpolyn_struct * coeffs;
    ulong * exps;
    slong alloc;
    slong length;
    mp_bitcnt_t bits;   /* default bits to construct coeffs */
} fmpz_mod_mpolyun_struct;
typedef fmpz_mod_mpolyun_struct fmpz_mod_mpolyun_t[1];



void fmpz_mod_mpolyn_init(
    fmpz_mod_mpolyn_t A,
    mp_bitcnt_t bits,
    const fmpz_mpoly_ctx_t ctx,
    const fmpz_mod_ctx_t fpctx)
{
    A->coeffs = NULL;
    A->exps = NULL;
    A->alloc = 0;
    A->length = 0;
    A->bits = bits;
}

void fmpz_mod_mpolyn_clear(
    fmpz_mod_mpolyn_t A,
    const fmpz_mpoly_ctx_t ctx,
    const fmpz_mod_ctx_t fpctx)
{
    slong i;
    for (i = 0; i < A->alloc; i++)
        fmpz_mod_poly_clear(A->coeffs + i);
    if (A->coeffs)
        flint_free(A->coeffs);
    if (A->exps)
        flint_free(A->exps);
}

void fmpz_mod_mpolyn_swap(fmpz_mod_mpolyn_t A, fmpz_mod_mpolyn_t B)
{
   fmpz_mod_mpolyn_struct t = *A;
   *A = *B;
   *B = t;
}

void fmpz_mod_mpolyn_zero(
    fmpz_mod_mpolyn_t A,
    const fmpz_mpoly_ctx_t ctx,
    const fmpz_mod_ctx_t fpctx)
{
    A->length = 0;
}

int fmpz_mod_mpolyn_is_zero(
    fmpz_mod_mpolyn_t A,
    const fmpz_mpoly_ctx_t ctx,
    const fmpz_mod_ctx_t fpctx)
{
    return A->length == 0;
}

void fmpz_mod_mpolyn_print_pretty(
    const fmpz_mod_mpolyn_t A,
    const char ** x_in,
    const fmpz_mpoly_ctx_t ctx,
    const fmpz_mod_ctx_t fpctx)
{
    fmpz_mod_poly_struct * coeff = A->coeffs;
    slong len = A->length;
    ulong * exp = A->exps;
    slong bits = A->bits;
    slong i, j, N;
    fmpz * exponents;
    char ** x = (char **) x_in;
    TMP_INIT;

    if (len == 0)
    {
        flint_printf("0");
        return;
    }

    N = mpoly_words_per_exp(bits, ctx->minfo);

    TMP_START;

    if (x == NULL)
    {
        x = (char **) TMP_ALLOC(ctx->minfo->nvars*sizeof(char *));
        for (i = 0; i < ctx->minfo->nvars; i++)
        {
            x[i] = (char *) TMP_ALLOC(((FLINT_BITS+4)/3)*sizeof(char));
            flint_sprintf(x[i], "x%wd", i+1);
        }
    }

    exponents = (fmpz *) TMP_ALLOC(ctx->minfo->nvars*sizeof(ulong));
    for (i = 0; i < ctx->minfo->nvars; i++)
        fmpz_init(exponents + i);
   
    for (i = 0; i < len; i++)
    {
        if (i > 0)
        {
            printf(" + ");
        }

        printf("(");
        fmpz_mod_poly_print_pretty(coeff + i, "v");
        printf(")");

        mpoly_get_monomial_ffmpz(exponents, exp + N*i, bits, ctx->minfo);

        for (j = 0; j < ctx->minfo->nvars; j++)
        {
            int cmp = fmpz_cmp_ui(exponents + j, WORD(1));

            if (cmp > 0)
            {
                printf("*%s^", x[j]);
                fmpz_print(exponents + j);
            } else if (cmp == 0)
            {
                printf("*%s", x[j]);
            }
        }
    }

    for (i = 0; i < ctx->minfo->nvars; i++)
        fmpz_clear(exponents + i);

    TMP_END;
}

void fmpz_mod_mpolyn_fit_length(
    fmpz_mod_mpolyn_t A,
    slong length,
    const fmpz_mpoly_ctx_t ctx,
    const fmpz_mod_ctx_t fpctx)
{
    slong i;
    slong old_alloc = A->alloc;
    slong new_alloc = FLINT_MAX(length, 2*A->alloc);

    if (length > old_alloc)
    {
        slong N = mpoly_words_per_exp(A->bits, ctx->minfo);

        if (old_alloc == 0)
        {
            A->exps = (ulong *) flint_malloc(new_alloc*N*sizeof(ulong));
            A->coeffs = (fmpz_mod_poly_struct *) flint_malloc(new_alloc*sizeof(fmpz_mod_poly_struct));
        } else
        {
            A->exps = (ulong *) flint_realloc(A->exps, new_alloc*N*sizeof(ulong));
            A->coeffs = (fmpz_mod_poly_struct *) flint_realloc(A->coeffs, new_alloc*sizeof(fmpz_mod_poly_struct));
        }

        for (i = old_alloc; i < new_alloc; i++)
        {
            fmpz_mod_poly_init(A->coeffs + i, fpctx->p);
        }
        A->alloc = new_alloc;
    }
}




void fmpz_mod_mpolyun_init(
    fmpz_mod_mpolyun_t A,
    mp_bitcnt_t bits,
    const fmpz_mpoly_ctx_t ctx,
    const fmpz_mod_ctx_t fpctx)
{
    A->coeffs = NULL;
    A->exps = NULL;
    A->alloc = 0;
    A->length = 0;
    A->bits = bits;
}

void fmpz_mod_mpolyun_clear(
    fmpz_mod_mpolyun_t A,
    const fmpz_mpoly_ctx_t ctx,
    const fmpz_mod_ctx_t fpctx)
{
    slong i;
    for (i = 0; i < A->alloc; i++)
        fmpz_mod_mpolyn_clear(A->coeffs + i, ctx, fpctx);
    if (A->coeffs)
        flint_free(A->coeffs);
    if (A->exps)
        flint_free(A->exps);
}


void fmpz_mod_mpolyun_set_mod(fmpz_mod_mpolyun_t A, const fmpz_mod_ctx_t fpctx)
{
    slong i, j;

    for (i = 0; i < A->alloc; i++)
    {
        fmpz_mod_mpolyn_struct * Ac = A->coeffs + i;
        for (j = 0; j < Ac->alloc; j++)
        {
            fmpz_set(&(Ac->coeffs + j)->p, fpctx->p);
        }
    }
}

/*
    get the leading coeff in x_0,...,x_var
    A is in R[x_0, ... x_(var-1)][x_var]
*/
fmpz * fmpz_mod_mpolyn_leadcoeff_last_ref(
    const fmpz_mod_mpolyn_t A,
    const fmpz_mpoly_ctx_t ctx,
    const fmpz_mod_ctx_t fpctx)
{
    fmpz_mod_poly_struct * leadpoly;
    FLINT_ASSERT(A->length > 0);
    leadpoly = A->coeffs + 0;
    FLINT_ASSERT(leadpoly->length > 0);
    return leadpoly->coeffs + leadpoly->length - 1;
}

fmpz * fmpz_mod_mpolyun_leadcoeff_last_ref(
    const fmpz_mod_mpolyun_t A,
    const fmpz_mpoly_ctx_t ctx,
    const fmpz_mod_ctx_t fpctx)
{
    FLINT_ASSERT(A->length > 0);
    return fmpz_mod_mpolyn_leadcoeff_last_ref(A->coeffs + 0, ctx, fpctx);
}

fmpz_mod_poly_struct * fmpz_mod_mpolyn_leadcoeff_ref(
    fmpz_mod_mpolyn_t A,
    const fmpz_mpoly_ctx_t ctx,
    const fmpz_mod_ctx_t fpctx)
{
    FLINT_ASSERT(A->length > 0);
    return A->coeffs + 0;
}

fmpz_mod_poly_struct * fmpz_mod_mpolyun_leadcoeff_ref(
    fmpz_mod_mpolyun_t A,
    const fmpz_mpoly_ctx_t ctx,
    const fmpz_mod_ctx_t fpctx)
{
    FLINT_ASSERT(A->length > 0);
    return fmpz_mod_mpolyn_leadcoeff_ref(A->coeffs + 0, ctx, fpctx);
}

void fmpz_mod_mpolyun_swap(fmpz_mod_mpolyun_t A, fmpz_mod_mpolyun_t B)
{
   fmpz_mod_mpolyun_struct t = *A;
   *A = *B;
   *B = t;
}

void fmpz_mod_mpolyun_zero(
    fmpz_mod_mpolyun_t A,
    const fmpz_mpoly_ctx_t ctx,
    const fmpz_mod_ctx_t fpctx)
{
    slong i;
    for (i = 0; i < A->alloc; i++) {
        fmpz_mod_mpolyn_clear(A->coeffs + i, ctx, fpctx);
        fmpz_mod_mpolyn_init(A->coeffs + i, A->bits, ctx, fpctx);
    }
    A->length = 0;
}

void fmpz_mod_mpolyun_print_pretty(
    const fmpz_mod_mpolyun_t poly,
    const char ** x,
    const fmpz_mpoly_ctx_t ctx,
    const fmpz_mod_ctx_t fpctx)
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
        fmpz_mod_mpolyn_print_pretty(poly->coeffs + i,x,ctx, fpctx);
        flint_printf(")*X^%wd",poly->exps[i]);
    }
}

void fmpz_mod_mpolyun_fit_length(
    fmpz_mod_mpolyun_t A,
    slong length,
    const fmpz_mpoly_ctx_t ctx,
    const fmpz_mod_ctx_t fpctx)
{
    slong i;
    slong old_alloc = A->alloc;
    slong new_alloc = FLINT_MAX(length, 2*A->alloc);

    if (length > old_alloc)
    {
        if (old_alloc == 0)
        {
            A->exps = (ulong *) flint_malloc(new_alloc*sizeof(ulong));
            A->coeffs = (fmpz_mod_mpolyn_struct *) flint_malloc(
                                          new_alloc*sizeof(fmpz_mod_mpolyn_struct));
        } else
        {
            A->exps = (ulong *) flint_realloc(A->exps,
                                                      new_alloc*sizeof(ulong));
            A->coeffs = (fmpz_mod_mpolyn_struct *) flint_realloc(A->coeffs,
                                          new_alloc*sizeof(fmpz_mod_mpolyn_struct));
        }

        for (i = old_alloc; i < new_alloc; i++)
        {
            fmpz_mod_mpolyn_init(A->coeffs + i, A->bits, ctx, fpctx);
        }
        A->alloc = new_alloc;
    }
}



void fmpz_mod_mpolyun_content_last(
    fmpz_mod_poly_t a,
    const fmpz_mod_mpolyun_t B,
    const fmpz_mpoly_ctx_t ctx,
    const fmpz_mod_ctx_t fpctx)
{
    slong i, j;
    fmpz_mod_poly_t t;

    fmpz_mod_poly_init(t, fpctx->p);

    fmpz_mod_poly_zero(a);
    for (i = 0; i < B->length; i++)
    {
        for (j = 0; j < (B->coeffs + i)->length; j++)
        {
            fmpz_mod_poly_gcd(t, a, (B->coeffs + i)->coeffs + j);
            fmpz_mod_poly_swap(t, a);
            if (fmpz_mod_poly_degree(a) == 0)
                break;
        }
    }

    fmpz_mod_poly_clear(t);
}


void fmpz_mod_mpolyun_divexact_last(
    fmpz_mod_mpolyun_t A,
    const fmpz_mod_poly_t b,
    const fmpz_mpoly_ctx_t ctx,
    const fmpz_mod_ctx_t fpctx)
{
    slong i, j;
    fmpz_mod_poly_t r, t;

    fmpz_mod_poly_init(r, fpctx->p);
    fmpz_mod_poly_init(t, fpctx->p);

    for (i = 0; i < A->length; i++)
    {
        fmpz_mod_poly_struct * Ac = (A->coeffs + i)->coeffs;
        for (j = 0; j < (A->coeffs + i)->length; j++)
        {
            fmpz_mod_poly_divrem(t, r, Ac + j, b);
            FLINT_ASSERT(fmpz_mod_poly_is_zero(r));
            FLINT_ASSERT(!fmpz_mod_poly_is_zero(t));
            fmpz_mod_poly_swap(t, Ac + j);
        }
    }
    fmpz_mod_poly_clear(r);
    fmpz_mod_poly_clear(t);
}

void fmpz_mod_mpolyun_mul_last(
    fmpz_mod_mpolyun_t A,
    fmpz_mod_poly_t b,
    const fmpz_mpoly_ctx_t ctx,
    const fmpz_mod_ctx_t fpctx)
{
    slong i, j;
    fmpz_mod_poly_t t;

    fmpz_mod_poly_init(t, fpctx->p);

    for (i = 0; i < A->length; i++)
    {
        for (j = 0; j < (A->coeffs + i)->length; j++)
        {
            fmpz_mod_poly_mul(t, (A->coeffs + i)->coeffs + j, b);
            fmpz_mod_poly_swap(t, (A->coeffs + i)->coeffs + j);
        }
    }

    fmpz_mod_poly_clear(t);
}



slong fmpz_mod_mpolyun_lastdeg(
    const fmpz_mod_mpolyun_t A,
    const fmpz_mpoly_ctx_t ctx,
    const fmpz_mod_ctx_t fpctx)
{
    slong i, j;
    slong deg = -WORD(1);

    for (i = 0; i < A->length; i++)
    {
        for (j = 0; j < (A->coeffs + i)->length; j++)
        {
            slong newdeg = fmpz_mod_poly_degree((A->coeffs + i)->coeffs + j);
            deg = FLINT_MAX(deg, newdeg);
        }
    }

    return deg;
}

/*
    E = A(v = alpha)
    A is in R[X][v]
    E is in R[X]
*/
void fmpz_mod_mpolyun_eval_last_bivar(
    fmpz_mod_poly_t E,
    const fmpz_mod_mpolyun_t A,
    const fmpz_t alpha,
    const fmpz_mpoly_ctx_t ctx,
    const fmpz_mod_ctx_t fpctx)
{
    fmpz_t v;
    slong Ai, Alen;
    fmpz_mod_mpolyn_struct * Acoeff;
    ulong * Aexp;

    fmpz_init(v);

    Acoeff = A->coeffs;
    Aexp = A->exps;
    Alen = A->length;
    Ai = 0;
    fmpz_mod_poly_zero(E);
    for (Ai = 0; Ai < Alen; Ai++)
    {
        fmpz_mod_poly_evaluate_fmpz(v, (Acoeff + Ai)->coeffs + 0, alpha);
        fmpz_mod_poly_set_coeff_fmpz(E, Aexp[Ai], v);
    }

    fmpz_clear(v);
}



void fmpz_mod_mpolyn_one(
    fmpz_mod_mpolyn_t A,
    const fmpz_mpoly_ctx_t ctx,
    const fmpz_mod_ctx_t fpctx)
{
    fmpz_mod_poly_struct * Acoeff;
    ulong * Aexp;
    slong N;

    fmpz_mod_mpolyn_fit_length(A, 1, ctx, fpctx);
    Acoeff = A->coeffs;
    Aexp = A->exps;

    N = mpoly_words_per_exp(A->bits, ctx->minfo);

    fmpz_mod_poly_set_ui(Acoeff + 0, 1);
    mpoly_monomial_zero(Aexp + N*0, N);

    A->length = 1;
}

void fmpz_mod_mpolyun_one(
    fmpz_mod_mpolyun_t A,
    const fmpz_mpoly_ctx_t ctx,
    const fmpz_mod_ctx_t fpctx)
{
    fmpz_mod_mpolyun_fit_length(A, 1, ctx, fpctx);
    fmpz_mod_mpolyn_one(A->coeffs + 0, ctx, fpctx);
    A->exps[0] = 0;
    A->length = 1;
}





void fmpz_mod_mpolyn_scalar_mul_fmpz_mod(
    fmpz_mod_mpolyn_t A,
    const fmpz_t c,
    const fmpz_mpoly_ctx_t ctx,
    const fmpz_mod_ctx_t fpctx)
{
    slong i;
    for (i = 0; i < A->length; i++)
    {
        fmpz_mod_poly_scalar_mul_fmpz(A->coeffs + i, A->coeffs + i, c);
    }
}

void fmpz_mod_mpolyun_scalar_mul_fmpz_mod(
    fmpz_mod_mpolyun_t A,
    const fmpz_t c,
    const fmpz_mpoly_ctx_t ctx,
    const fmpz_mod_ctx_t fpctx)
{
    slong i;
    FLINT_ASSERT(!fmpz_is_zero(c));
    for (i = 0; i < A->length; i++)
    {
        fmpz_mod_mpolyn_scalar_mul_fmpz_mod(A->coeffs + i, c, ctx, fpctx);
    }
}


typedef struct {
    slong half_point_count;
    fmpz_mod_poly_t R0, R1;
    fmpz_mod_poly_t V0, V1; /* V1 is our master polynomial */
    fmpz_mod_poly_t r; /* temporary */
    fmpz_mod_poly_t q; /* temporary also used as a queue of incoming points */
} fmpz_mod_bma_struct;
typedef fmpz_mod_bma_struct fmpz_mod_bma_t[1];

/*
    A = x^(2*n)       n = half_point_count
    deg(S) < 2*n
    U0*A + V0*S = R0   deg(R0) >= n
    U1*A + V1*S = R1   deg(R1) < n

    S is the reverse of the polynomial whose coefficients are the input points.
        S = a_0*x^(2n-1) + a_1*x^(2n-2) + ... + a_(2n-1)
    S can be updated with more points at any time via add_points.
    S can be reduced (compute V0, V1, R0, R1) at any time via reduce, which
        returns whether the reduction caused a change.

    The U0 and U1 are not stored.
*/
void fmpz_mod_bma_init(fmpz_mod_bma_t B, const fmpz_t p)
{
    B->half_point_count = 0;
    fmpz_mod_poly_init(B->V0, p);
    fmpz_mod_poly_init(B->R0, p);
    fmpz_mod_poly_set_ui(B->R0, 1);
    fmpz_mod_poly_init(B->V1, p);
    fmpz_mod_poly_set_ui(B->V1, 1);
    fmpz_mod_poly_init(B->R1, p);
    fmpz_mod_poly_init(B->r, p);
    fmpz_mod_poly_init(B->q, p);
    B->q->length = 0;
}

void fmpz_mod_bma_start_over(fmpz_mod_bma_t B)
{
    B->half_point_count = 0;
    fmpz_mod_poly_zero(B->V0);
    fmpz_mod_poly_set_ui(B->R0, 1);
    fmpz_mod_poly_set_ui(B->V1, 1);
    fmpz_mod_poly_zero(B->R1);
    B->q->length = 0;
}

void fmpz_mod_bma_reset_prime(fmpz_mod_bma_t B, const fmpz_t p)
{
    fmpz_set(&B->R0->p, p);
    fmpz_set(&B->R1->p, p);
    fmpz_set(&B->V0->p, p);
    fmpz_set(&B->V1->p, p);
    fmpz_set(&B->r->p, p);
    fmpz_set(&B->q->p, p);
    fmpz_mod_bma_start_over(B);
}

void fmpz_mod_bma_clear(fmpz_mod_bma_t B)
{
    fmpz_mod_poly_clear(B->R0);
    fmpz_mod_poly_clear(B->R1);
    fmpz_mod_poly_clear(B->V0);
    fmpz_mod_poly_clear(B->V1);
    fmpz_mod_poly_clear(B->r);
    fmpz_mod_poly_clear(B->q);
}

void fmpz_mod_bma_add_points(fmpz_mod_bma_t B, const fmpz * a, slong count)
{
    slong i;
    slong old_length = B->q->length;
    fmpz_mod_poly_fit_length(B->q, old_length + count);
    for (i = 0; i < count; i++)
    {
        fmpz_set(B->q->coeffs + old_length + i, a + i);
    }
    B->q->length = old_length + count;
}

void fmpz_mod_bma_add_point(fmpz_mod_bma_t B, const fmpz_t a)
{
    slong old_length = B->q->length;
    fmpz_mod_poly_fit_length(B->q, old_length + 1);
    fmpz_set(B->q->coeffs + old_length, a);
    B->q->length = old_length + 1;
}

int fmpz_mod_bma_reduce(fmpz_mod_bma_t B)
{
    int changed = 0;
    slong i, queue_length = B->q->length;

    if ((queue_length % 2) != 0)
    {
        flint_printf("Exception in nmod_bma_reduce: point count is not even\n");
        flint_abort();
    }

    /* reverse the queue into temp r */
    B->half_point_count += queue_length/2;
    fmpz_mod_poly_zero(B->r);
    for (i = 0; i < queue_length; i++)
    {
        fmpz_mod_poly_set_coeff_fmpz(B->r, queue_length - i - 1, B->q->coeffs + i);
    }
    fmpz_mod_poly_mul(B->q, B->V0, B->r);
    fmpz_mod_poly_shift_left(B->R0, B->R0, queue_length);
    fmpz_mod_poly_add(B->R0, B->R0, B->q);
    fmpz_mod_poly_mul(B->q, B->V1, B->r);
    fmpz_mod_poly_shift_left(B->R1, B->R1, queue_length);
    fmpz_mod_poly_add(B->R1, B->R1, B->q);

    /* now reduce */
    while (B->half_point_count < fmpz_mod_poly_length(B->R1))
    {
        changed = 1;
        fmpz_mod_poly_divrem(B->q, B->r, B->R0, B->R1);
        fmpz_mod_poly_swap(B->R0, B->R1);
        fmpz_mod_poly_swap(B->R1, B->r);

        fmpz_mod_poly_mul(B->r, B->q, B->V1);
        fmpz_mod_poly_sub(B->q, B->V0, B->r);
        fmpz_mod_poly_swap(B->V0, B->V1);
        fmpz_mod_poly_swap(B->V1, B->q);
        FLINT_ASSERT(fmpz_mod_poly_degree(B->V1) > fmpz_mod_poly_degree(B->V0));
    }

    /* queue is empty now */
    B->q->length = 0;
    return changed;
}



/* split f assuming that f has degree(f) distinct nonzero roots in Fp */
static void _fmpz_mod_poly_rabinsplit(fmpz_mod_poly_t a, fmpz_mod_poly_t b, fmpz_mod_poly_t T,
                                  const fmpz_mod_poly_t f, flint_rand_t randstate)
{
    fmpz_t delta;

    fmpz_init(delta);

    FLINT_ASSERT(fmpz_mod_poly_degree(f) > 1);

try_again:

    fmpz_randm(delta, randstate, &f->p);

    fmpz_mod_poly_zero(a);
    fmpz_mod_poly_set_coeff_ui(a, 1, 1);
    fmpz_mod_poly_set_coeff_fmpz(a, 0, delta);
    fmpz_sub_ui(delta, &f->p, 1);
    fmpz_divexact_ui(delta, delta, 2);
    fmpz_mod_poly_powmod_fmpz_binexp(T, a, delta, f);
    fmpz_mod_poly_zero(a);
    fmpz_mod_poly_set_coeff_ui(a, 0, 1);
    fmpz_mod_poly_sub(T, T, a);
    fmpz_mod_poly_gcd(a, T, f);
    FLINT_ASSERT(!fmpz_mod_poly_is_zero(a));
    if (0 >= fmpz_mod_poly_degree(a) || fmpz_mod_poly_degree(a) >= fmpz_mod_poly_degree(f))
    {
        goto try_again;
    }
    fmpz_mod_poly_div_basecase(b, f, a);
    /* deg a >= deg b */
    if (fmpz_mod_poly_degree(a) < fmpz_mod_poly_degree(b))
    {
        fmpz_mod_poly_swap(a, b);
    }

    fmpz_clear(delta);
    return;
}

/* fill in roots with the t distinct nonzero roots of master, or fail */
static int _fmpz_mod_find_roots(fmpz * roots, const fmpz_mod_poly_t master, slong t, const fmpz_mod_ctx_t fpctx)
{
    fmpz_t a0, a1;
    int success;
    slong i, roots_idx, sp;
    fmpz_t delta;
    fmpz_mod_poly_struct * a , * b;
    fmpz_mod_poly_t f, T;
    fmpz_mod_poly_struct stack[FLINT_BITS + 1];
    flint_rand_t randstate;

    FLINT_ASSERT(t >= 0);
    FLINT_ASSERT(t == fmpz_mod_poly_degree(master));

    fmpz_init(a0);
    fmpz_init(a1);
    fmpz_init(delta);

    if (t == 0)
    {
        success = 1;
        goto cleanup1;
    }
    else if (t == 1)
    {
        fmpz_mod_poly_get_coeff_fmpz(a0, master, 0);
        fmpz_mod_poly_get_coeff_fmpz(a1, master, 1);
        if (fmpz_is_zero(a0))
        {
            success = 0;
            goto cleanup1;
        }
        fmpz_mod_inv(a1, a1, fpctx);
        fmpz_mod_neg(a1, a1, fpctx);
        fmpz_mod_mul(roots + 0, a0, a1, fpctx);
        success = 1;
        goto cleanup1;
    }

    flint_randinit(randstate);
    fmpz_mod_poly_init(T, &master->p);
    fmpz_mod_poly_init(f, &master->p);
    for (i = 0; i <= FLINT_BITS; i++)
    {
        fmpz_mod_poly_init(stack + i, &master->p);
    }

    roots_idx = 0;

    fmpz_mod_poly_make_monic(f, master);

    a = stack + 0;
    fmpz_randm(delta, randstate, &master->p);
    fmpz_mod_poly_zero(a);
    fmpz_mod_poly_set_coeff_ui(a, 1, 1);
    fmpz_mod_poly_set_coeff_fmpz(a, 0, delta);
    fmpz_sub_ui(delta, &master->p, 1);
    fmpz_divexact_ui(delta, delta, 2);
    fmpz_mod_poly_powmod_fmpz_binexp(T, a, delta, f);
    fmpz_mod_poly_zero(a);
    fmpz_mod_poly_set_coeff_ui(a, 0, 1);
    fmpz_mod_poly_sub(T, T, a);
    fmpz_mod_poly_gcd(a, T, f);

    b = stack + 1;
    fmpz_mod_poly_zero(b);
    fmpz_mod_poly_set_coeff_ui(b, 0, 2);
    fmpz_mod_poly_add(T, T, b);
    fmpz_mod_poly_gcd(b, T, f);

    if (fmpz_mod_poly_degree(b) + fmpz_mod_poly_degree(a) != t)
    {
        success = 0;
        goto cleanup;
    }
    /* deg a >= deg b */
    if (fmpz_mod_poly_degree(a) < fmpz_mod_poly_degree(b))
    {
        fmpz_mod_poly_swap(a, b);
    }

    sp = fmpz_mod_poly_degree(b) > 0 ? 2 : 1;
    while (sp > 0)
    {
        FLINT_ASSERT(sp < FLINT_BITS);
        sp--;
        fmpz_mod_poly_swap(f, stack + sp);

        FLINT_ASSERT(fmpz_mod_poly_degree(f) > 0);
        if (fmpz_mod_poly_degree(f) == 1)
        {
            fmpz_mod_poly_get_coeff_fmpz(a0, f, 0);
            fmpz_mod_poly_get_coeff_fmpz(a1, f, 1);
            FLINT_ASSERT(!fmpz_is_zero(a0));
            FLINT_ASSERT(fmpz_is_one(a1));
            fmpz_mod_neg(roots + roots_idx, a0, fpctx);
            roots_idx++;
        }
        else
        {
            _fmpz_mod_poly_rabinsplit(stack + sp + 0, stack + sp + 1, T, f, randstate);
            FLINT_ASSERT(FLINT_BIT_COUNT(fmpz_mod_poly_degree(stack + sp + 1)) <= FLINT_BITS - sp - 1);
            sp += 2;
        }
    }

    success = 1;

cleanup:

    flint_randclear(randstate);
    fmpz_mod_poly_clear(T);
    fmpz_mod_poly_clear(f);
    for (i = 0; i <= FLINT_BITS; i++)
    {
        fmpz_mod_poly_clear(stack + i);
    }

    if (success)
    {
        FLINT_ASSERT(roots_idx == t);
    }

cleanup1:

    fmpz_clear(a0);
    fmpz_clear(a1);
    fmpz_clear(delta);

    return success;
}


typedef struct {
    fmpz_mod_poly_t roots;
    fmpz_mod_poly_t evals;
    fmpz_mod_bma_t bma;
} fmpz_mod_mpoly_bma_interpolate_struct;
typedef fmpz_mod_mpoly_bma_interpolate_struct fmpz_mod_mpoly_bma_interpolate_t[1];


void fmpz_mod_mpoly_bma_interpolate_init(fmpz_mod_mpoly_bma_interpolate_t I, const fmpz_t p)
{
    fmpz_mod_poly_init(I->roots, p);
    fmpz_mod_poly_init(I->evals, p);
    fmpz_mod_bma_init(I->bma, p);
}

slong fmpz_mod_mpoly_bma_interpolate_pointcount(const fmpz_mod_mpoly_bma_interpolate_t I)
{
    return I->bma->q->length + 2*I->bma->half_point_count;
}

void fmpz_mod_mpoly_bma_interpolate_print(const fmpz_mod_mpoly_bma_interpolate_t I)
{
    flint_printf("(%wd) ", fmpz_mod_mpoly_bma_interpolate_pointcount(I));
    fmpz_mod_poly_print_pretty(I->bma->V1, "#");
}


void fmpz_mod_mpoly_bma_interpolate_reset_prime(
    fmpz_mod_mpoly_bma_interpolate_t I,
    const fmpz_t p)
{
    fmpz_set(&I->roots->p, p);
    fmpz_set(&I->evals->p, p);
    fmpz_mod_bma_reset_prime(I->bma, p);
}

void fmpz_mod_mpoly_bma_interpolate_clear(fmpz_mod_mpoly_bma_interpolate_t I)
{
    fmpz_mod_poly_clear(I->roots);
    fmpz_mod_poly_clear(I->evals);
    fmpz_mod_bma_clear(I->bma);
}

void fmpz_mod_mpoly_bma_interpolate_add_point(fmpz_mod_mpoly_bma_interpolate_t I, const fmpz_t a)
{
    fmpz_mod_poly_fit_length(I->evals, I->evals->length + 1);
    fmpz_set(I->evals->coeffs + I->evals->length, a);
    I->evals->length++;

    fmpz_mod_bma_add_point(I->bma, a);
}

void fmpz_mod_mpoly_bma_interpolate_zero(fmpz_mod_mpoly_bma_interpolate_t I, slong count)
{
    slong i;
    fmpz_mod_poly_fit_length(I->evals, count);
    fmpz_mod_bma_start_over(I->bma);
    for (i = 0; i < count; i++)
    {
        fmpz_zero(I->evals->coeffs + i);
        fmpz_mod_bma_add_point(I->bma, I->evals->coeffs + i);
    }
    I->evals->length = count;
}

int fmpz_mod_mpoly_bma_interpolate_reduce(fmpz_mod_mpoly_bma_interpolate_t I)
{
    return fmpz_mod_bma_reduce(I->bma);
}



void fmpz_mpolyuu_print_pretty(
    const fmpz_mpolyu_t poly,
    const char ** x,
    slong nmainvars,
    const fmpz_mpoly_ctx_t ctx)
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
        fmpz_mpoly_print_pretty(poly->coeffs + i, x, ctx);
        flint_printf(")");
        for (j = nmainvars - 1; j >= 0; j--)
        {
            flint_printf("*X%wd^%wd", nmainvars - 1 - j, 
                    mask & (poly->exps[i] >> (FLINT_BITS/nmainvars*j)));
        }
    }
}



/*
    set out to the evaluation of variables after ksub at alpha^w

    out[0]   = alpha ^ (w * subdegs[n-1] * subdegs[n-2] * ... * * subdegs[1])
      ...
    out[n-3] = alpha ^ (w * subdegs[n-1] * subdegs[n-2])
    out[n-2] = alpha ^ (w * subdegs[n-1])
    out[n-1] = alpha ^ (w)

    secret: subdegs[0] is not used
*/

void nmod_mpoly_bma_interpolate_alpha_powers(
    mp_limb_t * out,
    ulong w,
    const mpoly_bma_interpolate_ctx_t Ictx,
    const nmod_mpoly_ctx_t ctx)
{
    slong j = ctx->minfo->nvars - 1;
    out[j] = nmod_pow_ui(Ictx->dlogenv_sp->alpha, w, ctx->ffinfo->mod);
    for (; j > 0; j--)
    {
        out[j - 1] = nmod_pow_ui(out[j], Ictx->subdegs[j], ctx->ffinfo->mod);
    }
}

void fmpz_mod_mpoly_bma_interpolate_alpha_powers(
    fmpz * out,
    const fmpz_t w,
    const mpoly_bma_interpolate_ctx_t Ictx,
    const fmpz_mpoly_ctx_t ctx,
    const fmpz_mod_ctx_t fpctx)
{
    slong j = ctx->minfo->nvars - 1;
    fmpz_mod_pow_fmpz(out + j, Ictx->dlogenv->alpha, w, fpctx);
    for (; j > 0; j--)
    {
        fmpz_mod_pow_ui(out + j - 1, out + j, Ictx->subdegs[j], fpctx);
    }
}



/*
    Set the coefficients of E to the evaluation of coor monomials of A
    evaluation at x_i = alpha[i]
*/
void fmpz_mpoly_set_skel(
    fmpz_mpolyc_t M,
    const fmpz_mpoly_t A,
    const fmpz * alpha,
    const fmpz_mpoly_ctx_t ctx,
    const fmpz_mod_ctx_t fpctx)
{
    slong i, j;
    slong offset, shift;
    slong N = mpoly_words_per_exp_sp(A->bits, ctx->minfo);
    slong nvars = ctx->minfo->nvars;
    ulong * Aexp;
    fmpz * Mcoeff;
    slong * LUToffset;
    ulong * LUTmask;
    fmpz * LUTvalue;
    slong LUTlen;
    fmpz_t xpoweval;
    ulong * inputexpmask;
    TMP_INIT;

    TMP_START;

    fmpz_init(xpoweval);

    LUToffset = (slong *) TMP_ALLOC(N*FLINT_BITS*sizeof(slong));
    LUTmask   = (ulong *) TMP_ALLOC(N*FLINT_BITS*sizeof(ulong));
    LUTvalue  = (fmpz *) TMP_ALLOC(N*FLINT_BITS*sizeof(fmpz));
    for (i = 0; i < N*FLINT_BITS; i++)
    {
        fmpz_init(LUTvalue + i);
    }

    fmpz_mpolyc_fit_length(M, A->length);
    M->length = A->length;

    Mcoeff = M->coeffs;
    Aexp = A->exps;

    inputexpmask = (ulong *) TMP_ALLOC(N*sizeof(ulong));
    mpoly_monomial_zero(inputexpmask, N);
    for (i = 0; i < A->length; i++)
    {
        for (j = 0; j < N; j++)
        {
            inputexpmask[j] |= (Aexp + N*i)[j];
        }
    }

    LUTlen = 0;
    for (j = nvars - 1; j >= 0; j--)
    {
        mpoly_gen_offset_shift_sp(&offset, &shift, j, A->bits, ctx->minfo);

        fmpz_set(xpoweval, alpha + j); /* xpoweval = alpha[j]^(2^i) */
        for (i = 0; i < A->bits; i++)
        {
            LUToffset[LUTlen] = offset;
            LUTmask[LUTlen] = (UWORD(1) << (shift + i));
            fmpz_set(LUTvalue + LUTlen, xpoweval);
            if ((inputexpmask[offset] & LUTmask[LUTlen]) != 0)
            {
                LUTlen++;
            }
            fmpz_mod_mul(xpoweval, xpoweval, xpoweval, fpctx);
        }
    }
    FLINT_ASSERT(LUTlen < N*FLINT_BITS);

    for (i = 0; i < A->length; i++)
    {
        fmpz_one(xpoweval);
        for (j = 0; j < LUTlen; j++)
        {
            if (((Aexp + N*i)[LUToffset[j]] & LUTmask[j]) != 0)
            {
                fmpz_mod_mul(xpoweval, xpoweval, LUTvalue + j, fpctx);
            }
        }
        fmpz_set(Mcoeff + i, xpoweval);
    }

    fmpz_clear(xpoweval);
    for (i = 0; i < N*FLINT_BITS; i++)
    {
        fmpz_clear(LUTvalue + i);
    }
    TMP_END;
}

void fmpz_mpolyu_set_skel(
    fmpz_mpolycu_t M,
    const fmpz_mpolyu_t A,
    const fmpz * alpha,
    const fmpz_mpoly_ctx_t ctx,
    const fmpz_mod_ctx_t fpctx)
{
    slong i;
    fmpz_mpolycu_fit_length(M, A->length);
    M->length = A->length;
    for (i = 0; i < A->length; i++)
    {
        fmpz_mpoly_set_skel(M->coeffs + i, A->coeffs + i, alpha, ctx, fpctx);
    }
}



/* M = S */
void fmpz_mpoly_copy_skel(fmpz_mpolyc_t M, const fmpz_mpolyc_t S)
{
    fmpz_mpolyc_fit_length(M, S->length);
    M->length = S->length;
    _fmpz_vec_set(M->coeffs, S->coeffs, S->length);
}

void fmpz_mpolyu_copy_skel(fmpz_mpolycu_t M, const fmpz_mpolycu_t S)
{
    slong i;
    fmpz_mpolycu_fit_length(M, S->length);
    M->length = S->length;
    for (i = 0; i < S->length; i++)
    {
        fmpz_mpoly_copy_skel(M->coeffs + i, S->coeffs + i);
    }
}


/* M = S */
void fmpz_mpoly_red_skel(
    fmpz_mpolyc_t Ared,
    const fmpz_mpoly_t A,
    const fmpz_mod_ctx_t fpctx)
{
    fmpz_mpolyc_fit_length(Ared, A->length);
    Ared->length = A->length;
    _fmpz_vec_scalar_mod_fmpz(Ared->coeffs, A->coeffs, A->length, fpctx->p);
}

void fmpz_mpolyu_red_skel(
    fmpz_mpolycu_t Ared,
    const fmpz_mpolyu_t A,
    const fmpz_mod_ctx_t fpctx)
{
    slong i;
    fmpz_mpolycu_fit_length(Ared, A->length);
    Ared->length = A->length;
    for (i = 0; i < A->length; i++)
    {
        fmpz_mpoly_red_skel(Ared->coeffs + i, A->coeffs + i, fpctx);
    }
}


/*
    return A.M and multiply M by S
    the coefficients of A are not necesarily reduced mod fpctx->p
*/
void fmpz_mpoly_use_skel_mul(
    fmpz_t eval,
    fmpz_mpolyc_t Ared,
    fmpz_mpolyc_t M,
    const fmpz_mpolyc_t S,
    const fmpz_mod_ctx_t fpctx)
{
    slong i;
    fmpz_zero(eval);
    FLINT_ASSERT(M->length == Ared->length);
    for (i = 0; i < Ared->length; i++)
    {
        fmpz_addmul(eval, Ared->coeffs + i, M->coeffs + i);
        fmpz_mod_mul(M->coeffs + i, M->coeffs + i, S->coeffs + i, fpctx);
    }
    fmpz_mod(eval, eval, fpctx->p);
}


/*
    A is in ZZ[X,Y][x_0, ..., x_(n-1)]
    E is in Fp[X][Y]

    Set E to the evaluation of A. The evaluations of A's monomials are in M.
    M is then multiplied by S.
*/
void fmpz_mpolyuu_use_skel_mul(
    fmpz_mod_mpolyun_t E,
    const fmpz_mpolyu_t A,
    fmpz_mpolycu_t Ared,
    fmpz_mpolycu_t M,
    const fmpz_mpolycu_t S,
    const fmpz_mpoly_ctx_t ctx,
    const fmpz_mod_ctx_t fpctx)
{
    slong xexp, yexp;
    slong i;
    slong N = mpoly_words_per_exp_sp(E->bits, ctx->minfo);
    fmpz_t eval;

    fmpz_init(eval);

    FLINT_ASSERT(A->length == Ared->length);
    FLINT_ASSERT(A->length == M->length);
    FLINT_ASSERT(A->length == S->length);

    E->length = 0;
    for (i = 0; i < A->length; i++)
    {
        fmpz_mpoly_use_skel_mul(eval, Ared->coeffs + i, M->coeffs + i, S->coeffs + i, fpctx);
        if (fmpz_is_zero(eval))
        {
            continue;
        }

        xexp = A->exps[i] >> (FLINT_BITS/2);
        yexp = A->exps[i] & ((-UWORD(1)) >> (FLINT_BITS/2));

        if (E->length > 0 && E->exps[E->length - 1] == xexp)
        {
            fmpz_mod_poly_set_coeff_fmpz(E->coeffs[E->length-1].coeffs + 0, yexp, eval);
        }
        else
        {
            fmpz_mod_mpolyun_fit_length(E, E->length + 1, ctx, fpctx);
            fmpz_mod_mpolyn_fit_length(E->coeffs + E->length, 1, ctx, fpctx);
            (E->coeffs + E->length)->length = 1;
            mpoly_monomial_zero((E->coeffs + E->length)->exps, N);
            fmpz_mod_poly_zero((E->coeffs + E->length)->coeffs + 0);
            fmpz_mod_poly_set_coeff_fmpz((E->coeffs + E->length)->coeffs + 0, yexp, eval);
            E->exps[E->length] = xexp;
            E->length++;
        }
    }

    fmpz_clear(eval);
}



int fmpz_mod_mpoly_bma_interpolate_get_mpoly(
    fmpz_mpoly_t A,
    const fmpz_t alphashift,
    fmpz_mod_mpoly_bma_interpolate_t I,
    const mpoly_bma_interpolate_ctx_t Ictx,
    const fmpz_mpoly_ctx_t ctx,
    const fmpz_mod_ctx_t fpctx)
{
    slong i, j, t, N;
    int success;
    ulong this_exp;
    fmpz_t new_exp;
    slong * shifts, * offsets;
    fmpz * values, * roots;
    fmpz * Acoeff;
    ulong * Aexp;
    slong Alen;
    fmpz_t T, S, V, temp, halfp;
    TMP_INIT;

    TMP_START;

    fmpz_init(halfp);
    fmpz_init(T);
    fmpz_init(S);
    fmpz_init(V);
    fmpz_init(temp);
    fmpz_init(new_exp);

    fmpz_tdiv_q_2exp(halfp, fpctx->p, 1);

    fmpz_mod_bma_reduce(I->bma);
    t = fmpz_mod_poly_degree(I->bma->V1);
    FLINT_ASSERT(t > 0);
    FLINT_ASSERT(I->evals->length >= t);

    fmpz_mod_poly_fit_length(I->roots, t);
    I->roots->length = t;
    success = _fmpz_mod_find_roots(I->roots->coeffs, I->bma->V1, t, fpctx);
    if (!success)
    {
        goto cleanup;
    }

    roots = I->roots->coeffs;
    values = I->evals->coeffs;

    FLINT_ASSERT(ctx->minfo->ord == ORD_LEX);
    FLINT_ASSERT(A->bits <= FLINT_BITS);
    N = mpoly_words_per_exp_sp(A->bits, ctx->minfo);
    fmpz_mpoly_fit_length(A, t, ctx);
    Acoeff = A->coeffs;
    Aexp = A->exps;
    A->length = 0;

    shifts = (slong *) TMP_ALLOC(ctx->minfo->nvars);
    offsets = (slong *) TMP_ALLOC(ctx->minfo->nvars);
    for (j = 0; j < ctx->minfo->nvars; j++)
    {
        mpoly_gen_offset_shift_sp(offsets + j, shifts + j, j, A->bits, ctx->minfo);
    }

    Alen = 0;
    for (i = 0; i < t; i++)
    {
        /* coeffs[i] is (coeffs(P).values)/P(roots[i]) =: V/S
            where P(x) = master(x)/(x-roots[i])     */
        fmpz_zero(V);
        fmpz_zero(T);
        fmpz_zero(S);
        for (j = t; j > 0; j--)
        {
            fmpz_mod_mul(temp, roots + i, T, fpctx);
            fmpz_mod_add(T, temp, I->bma->V1->coeffs + j, fpctx);
            fmpz_mod_mul(temp, roots + i, S, fpctx);
            fmpz_mod_add(S, temp, T, fpctx);
            fmpz_mod_mul(temp, values + j - 1, T, fpctx);
            fmpz_mod_add(V, V, temp, fpctx);
        }
        /* roots[i] should be a root of master */
#if WANT_ASSERT
        fmpz_mod_mul(temp, roots + i, T, fpctx);
        fmpz_mod_add(temp, temp, I->bma->V1->coeffs + 0, fpctx);
        FLINT_ASSERT(fmpz_is_zero(temp));
#endif
        fmpz_mod_pow_fmpz(temp, roots + i, alphashift, fpctx);
        fmpz_mod_mul(S, S, temp, fpctx);
        fmpz_mod_inv(temp, S, fpctx);
        fmpz_mod_mul(Acoeff + Alen, V, temp, fpctx);
        if (fmpz_is_zero(Acoeff + Alen))
        {
            /* hmmm */
            continue;
        }

        if (fmpz_cmp(Acoeff + Alen, halfp) > 0)
        {
            fmpz_sub(Acoeff + Alen, Acoeff + Alen, fpctx->p);
        }


        mpoly_monomial_zero(Aexp + N*Alen, N);
        fmpz_mod_dlog_env_run(Ictx->dlogenv, new_exp, roots + i);
        for (j = ctx->minfo->nvars - 1; j >= 0; j--)
        {
            this_exp = fmpz_fdiv_ui(new_exp, Ictx->subdegs[j]);
            fmpz_fdiv_q_ui(new_exp, new_exp, Ictx->subdegs[j]);
            if (this_exp >= Ictx->degbounds[j])
            {
                success = 0;
                goto cleanup;
            }
            (Aexp + N*Alen)[offsets[j]] |= this_exp << shifts[j];
        }
        if (!fmpz_is_zero(new_exp))
        {
            success = 0;
            goto cleanup;
        }
        Alen++;
    }
    A->length = Alen;

    fmpz_mpoly_sort_terms(A, ctx);

    success = 1;

cleanup:

    fmpz_clear(T);
    fmpz_clear(S);
    fmpz_clear(V);
    fmpz_clear(temp);
    fmpz_clear(halfp);

    TMP_END;
    return success;
}





fmpz_mpoly_struct * _fmpz_mpolyu_get_coeff(fmpz_mpolyu_t A,
                             ulong pow, const fmpz_mpoly_ctx_t uctx);

/*
    Convert B to A using the variable permutation perm.
    The uctx should be the context of the coefficients of A.
    The ctx should be the context of B.

    operation on each term:

    for 0 <= k < m + 2
        l = perm[k]
        Aexp[k] = (Bexp[l] - shift[l])/stride[l]

    the most significant main variable uses Aexp[0]
    the least significant main variable uses Aexp[1]
    the coefficients of A use variables Aexp[2], ..., Aexp[m + 1]
*/
void fmpz_mpoly_to_mpolyuu_perm_deflate(
    fmpz_mpolyu_t A,
    const slong * perm,
    const ulong * shift,
    const ulong * stride,
    const fmpz_mpoly_ctx_t uctx,
    const fmpz_mpoly_t B,
    const fmpz_mpoly_ctx_t ctx)
{
    slong i, j, k, l;
    slong n = ctx->minfo->nvars;
    slong m = uctx->minfo->nvars;
    slong NA, NB;
    ulong * uexps;
    ulong * Bexps;
    fmpz_mpoly_struct * Ac;
    TMP_INIT;

    FLINT_ASSERT(A->bits <= FLINT_BITS);
    FLINT_ASSERT(B->bits <= FLINT_BITS);
    FLINT_ASSERT(m + 2 <= n);

    TMP_START;

    uexps = (ulong *) TMP_ALLOC((m + 2)*sizeof(ulong));
    Bexps = (ulong *) TMP_ALLOC(n*sizeof(ulong));

    fmpz_mpolyu_zero(A, uctx);

    NA = mpoly_words_per_exp(A->bits, uctx->minfo);
    NB = mpoly_words_per_exp(B->bits, ctx->minfo);

    for (j = 0; j < B->length; j++)
    {
        mpoly_get_monomial_ui(Bexps, B->exps + NB*j, B->bits, ctx->minfo);
        for (k = 0; k < m + 2; k++)
        {
            l = perm[k];
            FLINT_ASSERT(stride[l] != UWORD(0));
            FLINT_ASSERT(((Bexps[l] - shift[l]) % stride[l]) == UWORD(0));
            uexps[k] = (Bexps[l] - shift[l]) / stride[l];
        }
        FLINT_ASSERT(FLINT_BIT_COUNT(uexps[0]) < FLINT_BITS/2);
        FLINT_ASSERT(FLINT_BIT_COUNT(uexps[1]) < FLINT_BITS/2);
        Ac = _fmpz_mpolyu_get_coeff(A, (uexps[0] << (FLINT_BITS/2)) + uexps[1], uctx);
        FLINT_ASSERT(Ac->bits == A->bits);

        fmpz_mpoly_fit_length(Ac, Ac->length + 1, uctx);
        fmpz_set(Ac->coeffs + Ac->length, B->coeffs + j);
        mpoly_set_monomial_ui(Ac->exps + NA*Ac->length, uexps + 2, A->bits, uctx->minfo);
        Ac->length++;
    }

    for (i = 0; i < A->length; i++)
    {
        fmpz_mpoly_sort_terms(A->coeffs + i, uctx);
    }

    TMP_END;
}


/*
    Convert B to A using the variable permutation vector perm.
    A must be constructed with bits = Abits.

    operation on each term:

        for 0 <= l < n
            Aexp[l] = shift[l]

        for 0 <= k < m + 2
            l = perm[k]
            Aexp[l] += scale[l]*Bexp[k]
*/
void fmpz_mpoly_from_mpolyuu_perm_inflate( /* only for 2 main vars */
    fmpz_mpoly_t A,
    mp_bitcnt_t Abits,
    const fmpz_mpoly_ctx_t ctx,
    const fmpz_mpolyu_t B,
    const slong * perm,
    const ulong * shift,
    const ulong * stride,
    const fmpz_mpoly_ctx_t uctx)
{
    slong n = ctx->minfo->nvars;
    slong m = uctx->minfo->nvars;
    slong i, j, k, l;
    slong NA, NB;
    slong Alen;
    fmpz * Acoeff;
    ulong * Aexp;
    slong Aalloc;
    ulong * uexps;
    ulong * Aexps;
    TMP_INIT;

    FLINT_ASSERT(B->length > 0);
    FLINT_ASSERT(Abits <= FLINT_BITS);
    FLINT_ASSERT(B->bits <= FLINT_BITS);
    FLINT_ASSERT(m + 2 <= n);
    TMP_START;

    uexps = (ulong *) TMP_ALLOC((m + 2)*sizeof(ulong));
    Aexps = (ulong *) TMP_ALLOC(n*sizeof(ulong));

    NA = mpoly_words_per_exp(Abits, ctx->minfo);
    NB = mpoly_words_per_exp(B->bits, uctx->minfo);

    fmpz_mpoly_fit_bits(A, Abits, ctx);
    A->bits = Abits;

    Acoeff = A->coeffs;
    Aexp = A->exps;
    Aalloc = A->alloc;
    Alen = 0;
    for (i = 0; i < B->length; i++)
    {
        fmpz_mpoly_struct * Bc = B->coeffs + i;
        _fmpz_mpoly_fit_length(&Acoeff, &Aexp, &Aalloc, Alen + Bc->length, NA);
        FLINT_ASSERT(Bc->bits == B->bits);

        for (j = 0; j < Bc->length; j++)
        {
            fmpz_set(Acoeff + Alen + j, Bc->coeffs + j);
            mpoly_get_monomial_ui(uexps + 2, Bc->exps + NB*j, Bc->bits, uctx->minfo);
            uexps[0] = B->exps[i] >> (FLINT_BITS/2);
            uexps[1] = B->exps[i] & ((-UWORD(1)) >> (FLINT_BITS - FLINT_BITS/2));
            for (l = 0; l < n; l++)
            {
                Aexps[l] = shift[l];
            }
            for (k = 0; k < m + 2; k++)
            {
                l = perm[k];
                Aexps[l] += stride[l]*uexps[k];
            }
            mpoly_set_monomial_ui(Aexp + NA*(Alen + j), Aexps, Abits, ctx->minfo);
        }
        Alen += Bc->length;
    }
    A->coeffs = Acoeff;
    A->exps = Aexp;
    A->alloc = Aalloc;
    _fmpz_mpoly_set_length(A, Alen, ctx);

    fmpz_mpoly_sort_terms(A, ctx);
    TMP_END;
}



typedef struct {
    nmod_mpoly_bma_interpolate_struct * coeffs;
    ulong * exps;
    slong length;
    slong alloc;
    slong pointcount;
} nmod_bma_mpoly_struct;

typedef nmod_bma_mpoly_struct nmod_bma_mpoly_t[1];

void nmod_bma_mpoly_init(nmod_bma_mpoly_t A)
{
    A->length = 0;
    A->alloc = 0;
    A->exps = NULL;
    A->coeffs = NULL;
    A->pointcount = 0;
}

void nmod_bma_mpoly_reset_prime(
    nmod_bma_mpoly_t A,
    const nmodf_ctx_t fpctx)
{
    slong i;
    for (i = 0; i < A->alloc; i++)
    {
        nmod_mpoly_bma_interpolate_reset_prime(A->coeffs + i, fpctx->mod.n);
    }
}


void nmod_bma_mpoly_clear(nmod_bma_mpoly_t A)
{
    slong i;
    for (i = 0; i < A->alloc; i++)
    {
        nmod_mpoly_bma_interpolate_clear(A->coeffs + i);
    }

    if (A->exps)
        flint_free(A->exps);
    if (A->coeffs)
        flint_free(A->coeffs);
}

void nmod_bma_mpoly_print(nmod_bma_mpoly_t A)
{
    slong i;
    flint_printf("0");
    for (i = 0; i < A->length; i++)
    {
        flint_printf(" + [");
        nmod_mpoly_bma_interpolate_print(A->coeffs + i);
        flint_printf("]*X^%wd*Y^%wd", A->exps[i] >> (FLINT_BITS/2), A->exps[i] & ((-UWORD(1)) >> (FLINT_BITS/2)));
    }
}


void nmod_bma_mpoly_fit_length(
    nmod_bma_mpoly_t A,
    slong length,
    const nmodf_ctx_t fpctx)
{
    slong i;
    slong old_alloc = A->alloc;
    slong new_alloc = FLINT_MAX(length, A->alloc + A->alloc/2);

    if (length > old_alloc)
    {
        if (old_alloc == 0)
        {
            A->exps = (ulong *) flint_malloc(new_alloc*sizeof(ulong));
            A->coeffs = (nmod_mpoly_bma_interpolate_struct *) flint_malloc(new_alloc*sizeof(nmod_mpoly_bma_interpolate_struct));
        }
        else
        {
            A->exps = (ulong *) flint_realloc(A->exps, new_alloc*sizeof(ulong));
            A->coeffs = (nmod_mpoly_bma_interpolate_struct *) flint_realloc(A->coeffs, new_alloc*sizeof(nmod_mpoly_bma_interpolate_struct));
        }

        for (i = old_alloc; i < new_alloc; i++)
        {
            nmod_mpoly_bma_interpolate_init(A->coeffs + i, fpctx->mod.n);
        }
        A->alloc = new_alloc;
    }
}

void nmod_bma_mpoly_zero(nmod_bma_mpoly_t L)
{
    L->length = 0;
    L->pointcount = 0;
}





int nmod_bma_mpoly_reduce(nmod_bma_mpoly_t L)
{
    slong i;
    int changed;

    changed = 0;

    for (i = 0; i < L->length; i++)
    {
        FLINT_ASSERT(L->pointcount == nmod_mpoly_bma_interpolate_pointcount(L->coeffs + i));
        changed |= nmod_mpoly_bma_interpolate_reduce(L->coeffs + i);
    }

    return changed;
}

void nmod_bma_mpoly_add_point(
    nmod_bma_mpoly_t L,
    const nmod_mpolyun_t A,
    const nmodf_ctx_t fpctx)
{
    slong j;
    slong Alen = A->length;
    nmod_mpolyn_struct * Acoeff = A->coeffs;
    slong Li, Ai, ai;
    ulong Aexp;
    nmod_mpoly_bma_interpolate_struct * Lcoeff;
    slong Llen;
    ulong * Lexp;

    if (L->length == 0)
    {
        slong tot = 0;
        for (Ai = 0; Ai < Alen; Ai++)
        {
            for (ai = nmod_poly_degree((Acoeff + Ai)->coeffs + 0); ai >= 0; ai--)
            {
                tot += (0 != (Acoeff + Ai)->coeffs[0].coeffs[ai]);
            }
        }
        nmod_bma_mpoly_fit_length(L, tot, fpctx);
    }

    Lcoeff = L->coeffs;
    Llen = L->length;
    Lexp = L->exps;

    Li = 0;
    Ai = ai = 0;
    Aexp = 0;
    if (Ai < Alen)
    {
        ai = nmod_poly_degree((A->coeffs + Ai)->coeffs + 0);
        Aexp = (A->exps[Ai] << (FLINT_BITS/2)) + ai;
    }

    while (Li < Llen || Ai < Alen)
    {
        if (Li < Llen && Ai < Alen && Lexp[Li] == Aexp)
        {
add_same_exp:
            nmod_mpoly_bma_interpolate_add_point(Lcoeff + Li, (Acoeff + Ai)->coeffs[0].coeffs[ai]);
            Li++;
            do {
                ai--;
            } while (ai >= 0 && (0 == (Acoeff + Ai)->coeffs[0].coeffs[ai]));
            if (ai < 0)
            {
                Ai++;
                if (Ai < Alen)
                {
                    ai = nmod_poly_degree((A->coeffs + Ai)->coeffs + 0);        
                    Aexp = (A->exps[Ai] << (FLINT_BITS/2)) + ai;
                }
            }
            if (Ai < Alen)
            {
                Aexp = (A->exps[Ai] << (FLINT_BITS/2)) + ai;
            }
        }
        else if (Li < Llen && (Ai >= Alen || Lexp[Li] > Aexp))
        {
            nmod_mpoly_bma_interpolate_add_point(Lcoeff + Li, 0);            
        }
        else
        {
            FLINT_ASSERT(Li >= Llen || (Ai < Alen && Lexp[Li] < Aexp));
            {
                ulong texp;
                nmod_mpoly_bma_interpolate_struct tcoeff;

                nmod_bma_mpoly_fit_length(L, Llen + 1, fpctx);
                Lcoeff = L->coeffs;
                Lexp = L->exps;

                texp = Lexp[Llen];
                tcoeff = Lcoeff[Llen];
                for (j = Li; j < Llen; j++)
                {
                    Lexp[j + 1] = Lexp[j];
                    Lcoeff[j + 1] = Lcoeff[j];
                }
                Lexp[Li] = texp;
                Lcoeff[Li] = tcoeff;
            }

            nmod_mpoly_bma_interpolate_zero(Lcoeff + Li, L->pointcount);
            Lexp[Li] = Aexp;
            Llen++;
            L->length = Llen;
            goto add_same_exp;
        }
    }

    L->pointcount++;
}

int nmod_bma_mpoly_get_mpolyu(
    fmpz_mpolyu_t A,
    const fmpz_mpoly_ctx_t ctx,
    ulong alphashift,
    const nmod_bma_mpoly_t L,
    const mpoly_bma_interpolate_ctx_t Ictx,
    const nmodf_ctx_t fpctx)
{
    int success;
    slong i;

    fmpz_mpolyu_fit_length(A, L->length, ctx);
    A->length = 0;
    for (i = 0; i < L->length; i++)
    {
        A->exps[A->length] = L->exps[i];
        success = nmod_mpoly_bma_interpolate_get_mpoly(A->coeffs + A->length, ctx, alphashift, L->coeffs + i, Ictx, fpctx);
        if (!success)
        {
            return 0;
        }
        A->length += !fmpz_mpoly_is_zero(A->coeffs + A->length, ctx);
    }
    return 1;
}








typedef struct {
    fmpz_mod_mpoly_bma_interpolate_struct * coeffs;
    ulong * exps;
    slong length;
    slong alloc;
    slong pointcount;
} fmpz_mod_bma_mpoly_struct;

typedef fmpz_mod_bma_mpoly_struct fmpz_mod_bma_mpoly_t[1];

void fmpz_mod_bma_mpoly_init(fmpz_mod_bma_mpoly_t A)
{
    A->length = 0;
    A->alloc = 0;
    A->exps = NULL;
    A->coeffs = NULL;
    A->pointcount = 0;
}

void fmpz_mod_bma_mpoly_reset_prime(
    fmpz_mod_bma_mpoly_t A,
    const fmpz_mod_ctx_t fpctx)
{
    slong i;
    for (i = 0; i < A->alloc; i++)
    {
        fmpz_mod_mpoly_bma_interpolate_reset_prime(A->coeffs + i, fpctx->p);
    }
}


void fmpz_mod_bma_mpoly_clear(fmpz_mod_bma_mpoly_t A)
{
    slong i;
    for (i = 0; i < A->alloc; i++)
    {
        fmpz_mod_mpoly_bma_interpolate_clear(A->coeffs + i);
    }

    if (A->exps)
        flint_free(A->exps);
    if (A->coeffs)
        flint_free(A->coeffs);
}

void fmpz_mod_bma_mpoly_print(
    fmpz_mod_bma_mpoly_t A,
    const mpoly_bma_interpolate_ctx_t Ictx)
{
    slong i;
    flint_printf("0");
    for (i = 0; i < A->length; i++)
    {
        flint_printf(" + [");
        fmpz_mod_mpoly_bma_interpolate_print(A->coeffs + i);
        flint_printf("]*X^%wd*Y^%wd", A->exps[i] >> (FLINT_BITS/2), A->exps[i] & ((-UWORD(1)) >> (FLINT_BITS/2)));
    }
}


void fmpz_mod_bma_mpoly_fit_length(
    fmpz_mod_bma_mpoly_t A,
    slong length,
    const fmpz_mod_ctx_t fpctx)
{
    slong i;
    slong old_alloc = A->alloc;
    slong new_alloc = FLINT_MAX(length, A->alloc + A->alloc/2);

    if (length > old_alloc)
    {
        if (old_alloc == 0)
        {
            A->exps = (ulong *) flint_malloc(new_alloc*sizeof(ulong));
            A->coeffs = (fmpz_mod_mpoly_bma_interpolate_struct *) flint_malloc(new_alloc*sizeof(fmpz_mod_mpoly_bma_interpolate_struct));
        }
        else
        {
            A->exps = (ulong *) flint_realloc(A->exps, new_alloc*sizeof(ulong));
            A->coeffs = (fmpz_mod_mpoly_bma_interpolate_struct *) flint_realloc(A->coeffs, new_alloc*sizeof(fmpz_mod_mpoly_bma_interpolate_struct));
        }

        for (i = old_alloc; i < new_alloc; i++)
        {
            fmpz_mod_mpoly_bma_interpolate_init(A->coeffs + i, fpctx->p);
        }
        A->alloc = new_alloc;
    }
}

void fmpz_mod_bma_mpoly_zero(fmpz_mod_bma_mpoly_t L)
{
    L->length = 0;
    L->pointcount = 0;
}

int fmpz_mod_bma_mpoly_reduce(fmpz_mod_bma_mpoly_t L)
{
    slong i;
    int changed;

    changed = 0;

    for (i = 0; i < L->length; i++)
    {
        FLINT_ASSERT(L->pointcount == fmpz_mod_mpoly_bma_interpolate_pointcount(L->coeffs + i));
        changed |= fmpz_mod_mpoly_bma_interpolate_reduce(L->coeffs + i);
    }

    return changed;
}

void fmpz_mod_bma_mpoly_add_point(
    fmpz_mod_bma_mpoly_t L,
    const fmpz_mod_mpolyun_t A,
    const fmpz_mpoly_ctx_t ctx,
    const fmpz_mod_ctx_t fpctx)
{
    slong j;
    slong Alen = A->length;
    fmpz_mod_mpolyn_struct * Acoeff = A->coeffs;
    slong Li, Ai, ai;
    ulong Aexp;
    fmpz_mod_mpoly_bma_interpolate_struct * Lcoeff;
    slong Llen;
    ulong * Lexp;
    fmpz_t zero;

    fmpz_init(zero);

    if (L->length == 0)
    {
        slong tot = 0;
        for (Ai = 0; Ai < Alen; Ai++)
        {
            for (ai = fmpz_mod_poly_degree((Acoeff + Ai)->coeffs + 0); ai >= 0; ai--)
            {
                tot += !fmpz_is_zero((Acoeff + Ai)->coeffs[0].coeffs + ai);
            }
        }
        fmpz_mod_bma_mpoly_fit_length(L, tot, fpctx);
    }

    Lcoeff = L->coeffs;
    Llen = L->length;
    Lexp = L->exps;

    Li = 0;
    Ai = ai = 0;
    Aexp = 0;
    if (Ai < Alen)
    {
        ai = fmpz_mod_poly_degree((A->coeffs + Ai)->coeffs + 0);
        Aexp = (A->exps[Ai] << (FLINT_BITS/2)) + ai;
    }

    while (Li < Llen || Ai < Alen)
    {
/*
flint_printf("Li: %wd, Ai: %wd\n", Li, Ai);
usleep(1000000);
*/
        if (Li < Llen && Ai < Alen && Lexp[Li] == Aexp)
        {
add_same_exp:
            fmpz_mod_mpoly_bma_interpolate_add_point(Lcoeff + Li, (Acoeff + Ai)->coeffs[0].coeffs + ai);
            Li++;
            do {
                ai--;
            } while (ai >= 0 && fmpz_is_zero((Acoeff + Ai)->coeffs[0].coeffs + ai));
            if (ai < 0)
            {
                Ai++;
                if (Ai < Alen)
                {
                    ai = fmpz_mod_poly_degree((A->coeffs + Ai)->coeffs + 0);        
                    Aexp = (A->exps[Ai] << (FLINT_BITS/2)) + ai;
                }
            }
            if (Ai < Alen)
            {
                Aexp = (A->exps[Ai] << (FLINT_BITS/2)) + ai;
            }
        }
        else if (Li < Llen && (Ai >= Alen || Lexp[Li] > Aexp))
        {
            fmpz_mod_mpoly_bma_interpolate_add_point(Lcoeff + Li, zero);            
        }
        else
        {
            FLINT_ASSERT(Li >= Llen || (Ai < Alen && Lexp[Li] < Aexp));
            {
                ulong texp;
                fmpz_mod_mpoly_bma_interpolate_struct tcoeff;

                fmpz_mod_bma_mpoly_fit_length(L, Llen + 1, fpctx);
                Lcoeff = L->coeffs;
                Lexp = L->exps;

                texp = Lexp[Llen];
                tcoeff = Lcoeff[Llen];
                for (j = Li; j < Llen; j++)
                {
                    Lexp[j + 1] = Lexp[j];
                    Lcoeff[j + 1] = Lcoeff[j];
                }
                Lexp[Li] = texp;
                Lcoeff[Li] = tcoeff;
            }

            fmpz_mod_mpoly_bma_interpolate_zero(Lcoeff + Li, L->pointcount);
            Lexp[Li] = Aexp;
            Llen++;
            L->length = Llen;
            goto add_same_exp;
        }
    }

    fmpz_clear(zero);
    L->pointcount++;
}

int fmpz_mod_bma_mpoly_get_mpolyu(
    fmpz_mpolyu_t A,
    const fmpz_t alphashift,
    const fmpz_mod_bma_mpoly_t L,
    const mpoly_bma_interpolate_ctx_t Ictx,
    const fmpz_mpoly_ctx_t ctx,
    const fmpz_mod_ctx_t fpctx)
{
    int success;
    slong i;

    fmpz_mpolyu_fit_length(A, L->length, ctx);
    A->length = 0;
    for (i = 0; i < L->length; i++)
    {
        A->exps[A->length] = L->exps[i];
        success = fmpz_mod_mpoly_bma_interpolate_get_mpoly(A->coeffs + A->length, alphashift, L->coeffs + i, Ictx, ctx, fpctx);
        if (!success)
        {
            return 0;
        }
        A->length += !fmpz_mpoly_is_zero(A->coeffs + A->length, ctx);
    }
    return 1;
}




/*
    A = B
    A, B are in R[X]
*/
void fmpz_mod_mpolyun_set_poly(
    fmpz_mod_mpolyun_t A,
    const fmpz_mod_poly_t B,
    const fmpz_mpoly_ctx_t ctx,
    const fmpz_mod_ctx_t fpctx)
{
    slong N = mpoly_words_per_exp_sp(A->bits, ctx->minfo);
    slong Bexp;
    slong Blen = fmpz_mod_poly_length(B);
    fmpz * Bcoeff = B->coeffs;
    fmpz_mod_mpolyn_struct * Acoeff;
    ulong * Aexp;
    slong Ai;

    fmpz_mod_mpolyun_fit_length(A, Blen, ctx, fpctx);
    Acoeff = A->coeffs;
    Aexp = A->exps;

    Ai = 0;
    for (Bexp = Blen - 1; Bexp >= 0; Bexp--)
    {
        if (!fmpz_is_zero(Bcoeff + Bexp))
        {
            FLINT_ASSERT(Ai < A->alloc);

            fmpz_mod_mpolyn_fit_length(Acoeff + Ai, 1, ctx, fpctx);
            mpoly_monomial_zero((Acoeff + Ai)->exps + N*0, N);
            fmpz_mod_poly_zero((Acoeff + Ai)->coeffs + 0);
            fmpz_mod_poly_set_coeff_fmpz((Acoeff + Ai)->coeffs + 0, 0, Bcoeff + Bexp);
            Aexp[Ai] = Bexp;
            (Acoeff + Ai)->length = 1;
            Ai++;
        }
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
int fmpz_mod_mpolyun_addinterp_bivar(
    slong * lastdeg_,
    fmpz_mod_mpolyun_t F,
    fmpz_mod_mpolyun_t T,
    const fmpz_mod_poly_t A,
    const fmpz_mod_poly_t modulus,
    const fmpz_t alpha,
    const fmpz_mpoly_ctx_t ctx,
    const fmpz_mod_ctx_t fpctx)
{
    int changed = 0;
    slong lastdeg = -WORD(1);
    slong N = mpoly_words_per_exp_sp(F->bits, ctx->minfo);
    fmpz_t v;
    slong Fi, Toff, Aexp;
    fmpz * Acoeff = A->coeffs;
    slong Flen = F->length;
    fmpz_mod_mpolyn_struct * Fcoeff = F->coeffs;
    ulong * Fexp = F->exps;
    fmpz_mod_mpolyn_struct * Tcoeff;
    ulong * Texp;
    fmpz_mod_poly_t tp;
    
    Fi = 0;

    fmpz_init(v);

/*
flint_printf("\nnmod_mpolyun_addinterp_bivar: (alpha = %wu)\n", alpha);
flint_printf("A: "); nmod_poly_print_pretty(A, "X"); printf("\n");
flint_printf("F: "); nmod_mpolyun_print_pretty(F, NULL, ctx); printf("\n");
*/

    Aexp = fmpz_mod_poly_degree(A);

    fmpz_mod_poly_init(tp, fpctx->p);

    fmpz_mod_mpolyun_fit_length(T, Flen + Aexp + 1, ctx, fpctx);
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
            FLINT_ASSERT(!fmpz_mod_poly_is_zero((Fcoeff + Fi)->coeffs + 0));
            FLINT_ASSERT(fmpz_mod_poly_degree((Fcoeff + Fi)->coeffs + 0) < fmpz_mod_poly_degree(modulus));
        }

        if (Aexp >= 0)
        {
            FLINT_ASSERT(Acoeff[Aexp] != UWORD(0));
        }

        fmpz_mod_mpolyn_fit_length(Tcoeff + Toff, 1, ctx, fpctx);

        if (Fi < Flen && Aexp >= 0 && Fexp[Fi] == Aexp)
        {
            /* F term ok, A term ok */
            fmpz_mod_poly_evaluate_fmpz(v, (Fcoeff + Fi)->coeffs + 0, alpha);
            fmpz_mod_sub(v, Acoeff + Aexp, v, fpctx);
            if (!fmpz_is_zero(v))
            {
                changed = 1;
                fmpz_mod_poly_scalar_mul_fmpz(tp, modulus, v);
                fmpz_mod_poly_add((Tcoeff + Toff)->coeffs + 0, (Fcoeff + Fi)->coeffs + 0, tp);
            }
            else
            {
                fmpz_mod_poly_set((Tcoeff + Toff)->coeffs + 0, (Fcoeff + Fi)->coeffs + 0);
            }
            Texp[Toff] = Aexp;
            Fi++;
            do {
                Aexp--;
            } while (Aexp >= 0 && fmpz_is_zero(Acoeff + Aexp));
        }
        else if (Fi < Flen && (Aexp < 0 || Fexp[Fi] > Aexp))
        {
            /* F term ok, A term missing */
            fmpz_mod_poly_evaluate_fmpz(v, (Fcoeff + Fi)->coeffs + 0, alpha);
            if (!fmpz_is_zero(v))
            {
                changed = 1;
                fmpz_mod_poly_scalar_mul_fmpz(tp, modulus, v);
                fmpz_mod_poly_sub((Tcoeff + Toff)->coeffs + 0, (Fcoeff + Fi)->coeffs + 0, tp);
            }
            else
            {
                fmpz_mod_poly_set((Tcoeff + Toff)->coeffs + 0, (Fcoeff + Fi)->coeffs + 0);                
            }
            Texp[Toff] = Fexp[Fi];
            Fi++;
        }
        else if (Aexp >= 0 && (Fi >= Flen || Fexp[Fi] < Aexp))
        {
            /* F term missing, A term ok */
            changed = 1;
            fmpz_mod_poly_scalar_mul_fmpz((Tcoeff + Toff)->coeffs + 0, modulus, Acoeff + Aexp);
            Texp[Toff] = Aexp;
            do {
                Aexp--;
            } while (Aexp >= 0 && fmpz_is_zero(Acoeff + Aexp));
        }
        else
        {
            FLINT_ASSERT(0);
        }

        lastdeg = FLINT_MAX(lastdeg, fmpz_mod_poly_degree((Tcoeff + Toff)->coeffs + 0));
        mpoly_monomial_zero((Tcoeff + Toff)->exps + N*0, N);
        FLINT_ASSERT(!fmpz_mod_poly_is_zero((Tcoeff + Toff)->coeffs + 0));
        (Tcoeff + Toff)->length = 1;
        Toff++;
    }
    T->length = Toff;

    fmpz_mod_poly_clear(tp);

    if (changed)
    {
        fmpz_mod_mpolyun_swap(T, F);
    }

    fmpz_clear(v);

    *lastdeg_ = lastdeg;
    return changed;
}


int fmpz_mod_mpolyun_gcd_bivar(
    fmpz_mod_mpolyun_t G,
    fmpz_mod_mpolyun_t Abar,
    fmpz_mod_mpolyun_t Bbar,
    fmpz_mod_mpolyun_t A,
    fmpz_mod_mpolyun_t B,
    const fmpz_mpoly_ctx_t ctx,
    const fmpz_mod_ctx_t fpctx)
{
    int success;
    slong bound;
    fmpz_t alpha, temp, gammaeval;
    fmpz_mod_poly_t Aeval, Beval, Geval, Abareval, Bbareval;
    fmpz_mod_mpolyun_t T;
    slong deggamma, ldegG, ldegAbar, ldegBbar, ldegA, ldegB;
    fmpz_mod_poly_t cA, cB, cG, cAbar, cBbar, gamma, r;
    fmpz_mod_poly_t modulus, modulus2;
#if WANT_ASSERT
    fmpz_mod_poly_t leadA, leadB;
#endif

    fmpz_init(alpha);
    fmpz_init(temp);
    fmpz_init(gammaeval);
/*
printf("fmpz_mod_mpolyun_gcd_bivar called\n");
printf("A: "); fmpz_mod_mpolyun_print_pretty(A, NULL, ctx, fpctx); printf("\n");
printf("B: "); fmpz_mod_mpolyun_print_pretty(B, NULL, ctx, fpctx); printf("\n");
*/

#if WANT_ASSERT
    fmpz_mod_poly_init(leadA, fpctx->p);
    fmpz_mod_poly_init(leadB, fpctx->p);
    fmpz_mod_poly_set(leadA, fmpz_mod_mpolyun_leadcoeff_ref(A, ctx, fpctx));
    fmpz_mod_poly_set(leadB, fmpz_mod_mpolyun_leadcoeff_ref(B, ctx, fpctx));
#endif

    fmpz_mod_poly_init(r, fpctx->p);
    fmpz_mod_poly_init(cA, fpctx->p);
    fmpz_mod_poly_init(cB, fpctx->p);
    fmpz_mod_mpolyun_content_last(cA, A, ctx, fpctx);
    fmpz_mod_mpolyun_content_last(cB, B, ctx, fpctx);
    fmpz_mod_mpolyun_divexact_last(A, cA, ctx, fpctx);
    fmpz_mod_mpolyun_divexact_last(B, cB, ctx, fpctx);

    fmpz_mod_poly_init(cG, fpctx->p);
    fmpz_mod_poly_gcd_euclidean(cG, cA, cB);

    fmpz_mod_poly_init(cAbar, fpctx->p);
    fmpz_mod_poly_init(cBbar, fpctx->p);
    fmpz_mod_poly_divrem(cAbar, r, cA, cG);
    FLINT_ASSERT(fmpz_mod_poly_is_zero(r));
    fmpz_mod_poly_divrem(cBbar, r, cB, cG);
    FLINT_ASSERT(fmpz_mod_poly_is_zero(r));

    fmpz_mod_poly_init(gamma, fpctx->p);
    fmpz_mod_poly_gcd(gamma, fmpz_mod_mpolyun_leadcoeff_ref(A, ctx, fpctx),
                             fmpz_mod_mpolyun_leadcoeff_ref(B, ctx, fpctx));

    ldegA = fmpz_mod_mpolyun_lastdeg(A, ctx, fpctx);
    ldegB = fmpz_mod_mpolyun_lastdeg(B, ctx, fpctx);
    deggamma = fmpz_mod_poly_degree(gamma);
    bound = 1 + deggamma + FLINT_MAX(ldegA, ldegB);
/*
flint_printf("bound: %wd\n", bound);
*/
    fmpz_mod_poly_init(Aeval, fpctx->p);
    fmpz_mod_poly_init(Beval, fpctx->p);
    fmpz_mod_poly_init(Geval, fpctx->p);
    fmpz_mod_poly_init(Abareval, fpctx->p);
    fmpz_mod_poly_init(Bbareval, fpctx->p);

    fmpz_mod_mpolyun_init(T, A->bits, ctx, fpctx);

    fmpz_mod_poly_init(modulus, fpctx->p);
    fmpz_mod_poly_init(modulus2, fpctx->p);
    fmpz_mod_poly_set_ui(modulus, 1);

    fmpz_sub_ui(alpha, fpctx->p, 2);

choose_prime: /* prime is v - alpha */

    if (fmpz_cmp_ui(alpha, 2) < 1)
    {
        success = 0;
        goto cleanup;
    }

    fmpz_sub_ui(alpha, alpha, 1);

/*
printf("alpha: "); fmpz_print(alpha); printf("\n");
*/

    /* make sure evaluation point does not kill both lc(A) and lc(B) */
    fmpz_mod_poly_evaluate_fmpz(gammaeval, gamma, alpha);
    if (fmpz_is_zero(gammaeval))
    {
        goto choose_prime;
    }

    /* evaluation point should kill neither A nor B */
    fmpz_mod_mpolyun_eval_last_bivar(Aeval, A, alpha, ctx, fpctx);
    fmpz_mod_mpolyun_eval_last_bivar(Beval, B, alpha, ctx, fpctx);
    FLINT_ASSERT(Aeval->length > 0);
    FLINT_ASSERT(Beval->length > 0);
/*
printf("Aeval: "); fmpz_mod_poly_print_pretty(Aeval, "X"); printf("\n");
printf("Beval: "); fmpz_mod_poly_print_pretty(Beval, "X"); printf("\n");
*/
    fmpz_mod_poly_gcd(Geval, Aeval, Beval);
    fmpz_mod_poly_divrem(Abareval, r, Aeval, Geval);
    FLINT_ASSERT(fmpz_mod_poly_is_zero(r));
    fmpz_mod_poly_divrem(Bbareval, r, Beval, Geval);
    FLINT_ASSERT(fmpz_mod_poly_is_zero(r));

/*
printf("   Geval: "); fmpz_mod_poly_print_pretty(Geval, "X"); printf("\n");
printf("Abareval: "); fmpz_mod_poly_print_pretty(Abareval, "X"); printf("\n");
printf("Bbareval: "); fmpz_mod_poly_print_pretty(Bbareval, "X"); printf("\n");
*/

    FLINT_ASSERT(Geval->length > 0);
    FLINT_ASSERT(Abareval->length > 0);
    FLINT_ASSERT(Bbareval->length > 0);

    if (fmpz_mod_poly_degree(Geval) == 0)
    {
        fmpz_mod_mpolyun_one(G, ctx, fpctx);
        fmpz_mod_mpolyun_swap(Abar, A);
        fmpz_mod_mpolyun_swap(Bbar, B);
        goto successful_put_content;    
    }

    if (fmpz_mod_poly_degree(modulus) > 0)
    {
        FLINT_ASSERT(G->length > 0);
        if (fmpz_mod_poly_degree(Geval) > G->exps[0])
        {
            goto choose_prime;
        }
        else if (fmpz_mod_poly_degree(Geval) < G->exps[0])
        {
            fmpz_mod_poly_set_ui(modulus, 1);
        }
    }

    fmpz_mod_poly_scalar_mul_fmpz(Geval, Geval, gammaeval);

    if (fmpz_mod_poly_degree(modulus) > 0)
    {
        fmpz_mod_poly_evaluate_fmpz(temp, modulus, alpha);
        fmpz_mod_inv(temp, temp, fpctx);
        fmpz_mod_poly_scalar_mul_fmpz(modulus, modulus, temp);
        fmpz_mod_mpolyun_addinterp_bivar(&ldegG, G, T, Geval, modulus, alpha, ctx, fpctx);
        fmpz_mod_mpolyun_addinterp_bivar(&ldegAbar, Abar, T, Abareval, modulus, alpha, ctx, fpctx);
        fmpz_mod_mpolyun_addinterp_bivar(&ldegBbar, Bbar, T, Bbareval, modulus, alpha, ctx, fpctx);
    }
    else
    {
        fmpz_mod_mpolyun_set_poly(G, Geval, ctx, fpctx);
        fmpz_mod_mpolyun_set_poly(Abar, Abareval, ctx, fpctx);
        fmpz_mod_mpolyun_set_poly(Bbar, Bbareval, ctx, fpctx);
        ldegG = 0;
        ldegAbar = 0;
        ldegBbar = 0;
    }

    fmpz_mod_poly_scalar_mul_fmpz(modulus2, modulus, alpha);
    fmpz_mod_poly_shift_left(modulus, modulus, 1);
    fmpz_mod_poly_sub(modulus, modulus, modulus2);
/*
printf("modulus: "); fmpz_mod_poly_print_pretty(modulus, "v"); printf("\n");
printf("   G: "); fmpz_mod_mpolyun_print_pretty(G, NULL, ctx, fpctx); printf("\n");
printf("Abar: "); fmpz_mod_mpolyun_print_pretty(Abar, NULL, ctx, fpctx); printf("\n");
printf("Bbar: "); fmpz_mod_mpolyun_print_pretty(Bbar, NULL, ctx, fpctx); printf("\n");
*/
    if (fmpz_mod_poly_degree(modulus) < bound)
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

    fmpz_mod_poly_set_ui(modulus, 1);
    goto choose_prime;

successful:
/*
printf("successful\n");
printf("   G: "); fmpz_mod_mpolyun_print_pretty(G, NULL, ctx, fpctx); printf("\n");
printf("Abar: "); fmpz_mod_mpolyun_print_pretty(Abar, NULL, ctx, fpctx); printf("\n");
printf("Bbar: "); fmpz_mod_mpolyun_print_pretty(Bbar, NULL, ctx, fpctx); printf("\n");
*/

    fmpz_mod_mpolyun_content_last(modulus, G, ctx, fpctx);
    fmpz_mod_mpolyun_divexact_last(G, modulus, ctx, fpctx);
    fmpz_mod_mpolyun_divexact_last(Abar, fmpz_mod_mpolyun_leadcoeff_ref(G, ctx, fpctx), ctx, fpctx);
    fmpz_mod_mpolyun_divexact_last(Bbar, fmpz_mod_mpolyun_leadcoeff_ref(G, ctx, fpctx), ctx, fpctx);

successful_put_content:
/*
printf("put_content\n");
printf("   G: "); fmpz_mod_mpolyun_print_pretty(G, NULL, ctx, fpctx); printf("\n");
printf("Abar: "); fmpz_mod_mpolyun_print_pretty(Abar, NULL, ctx, fpctx); printf("\n");
printf("Bbar: "); fmpz_mod_mpolyun_print_pretty(Bbar, NULL, ctx, fpctx); printf("\n");
*/

    fmpz_mod_mpolyun_mul_last(G, cG, ctx, fpctx);
    fmpz_mod_mpolyun_mul_last(Abar, cAbar, ctx, fpctx);
    fmpz_mod_mpolyun_mul_last(Bbar, cBbar, ctx, fpctx);

    success = 1;

cleanup:

/*
printf("cleanup\n");
printf("   G: "); fmpz_mod_mpolyun_print_pretty(G, NULL, ctx, fpctx); printf("\n");
printf("Abar: "); fmpz_mod_mpolyun_print_pretty(Abar, NULL, ctx, fpctx); printf("\n");
printf("Bbar: "); fmpz_mod_mpolyun_print_pretty(Bbar, NULL, ctx, fpctx); printf("\n");
*/

#if WANT_ASSERT
    if (success)
    {
        FLINT_ASSERT(fmpz_is_one(fmpz_mod_mpolyun_leadcoeff_last_ref(G, ctx, fpctx)));
        fmpz_mod_poly_mul(modulus, fmpz_mod_mpolyun_leadcoeff_ref(G, ctx, fpctx),
                                   fmpz_mod_mpolyun_leadcoeff_ref(Abar, ctx, fpctx));
        FLINT_ASSERT(fmpz_mod_poly_equal(modulus, leadA));
        fmpz_mod_poly_mul(modulus, fmpz_mod_mpolyun_leadcoeff_ref(G, ctx, fpctx),
                                   fmpz_mod_mpolyun_leadcoeff_ref(Bbar, ctx, fpctx));
        FLINT_ASSERT(fmpz_mod_poly_equal(modulus, leadB));
    }
    fmpz_mod_poly_clear(leadA);
    fmpz_mod_poly_clear(leadB);
#endif

    fmpz_mod_poly_clear(r);
    fmpz_mod_poly_clear(cA);
    fmpz_mod_poly_clear(cB);
    fmpz_mod_poly_clear(cG);
    fmpz_mod_poly_clear(cAbar);
    fmpz_mod_poly_clear(cBbar);

    fmpz_mod_poly_clear(gamma);

    fmpz_mod_poly_clear(Aeval);
    fmpz_mod_poly_clear(Beval);
    fmpz_mod_poly_clear(Geval);
    fmpz_mod_poly_clear(Abareval);
    fmpz_mod_poly_clear(Bbareval);

    fmpz_mod_mpolyun_clear(T, ctx, fpctx);

    fmpz_mod_poly_clear(modulus);
    fmpz_mod_poly_clear(modulus2);

    return success;
}


ulong fmpz_mod_mpolyun_bidegree(
    const fmpz_mod_mpolyun_t A,
    const fmpz_mpoly_ctx_t ctx,
    const fmpz_mod_ctx_t fpctx)
{
    ulong degx, degy;

    FLINT_ASSERT(A->length > 0);

    degx = A->exps[0];
    degy = fmpz_mod_poly_degree((A->coeffs + 0)->coeffs + 0);

    FLINT_ASSERT(FLINT_BIT_COUNT(degx) < FLINT_BITS/2);
    FLINT_ASSERT(FLINT_BIT_COUNT(degy) < FLINT_BITS/2);

    return (degx << (FLINT_BITS/2)) + degy;
}

ulong nmod_mpolyun_bidegree(
    const nmod_mpolyun_t A,
    const nmod_mpoly_ctx_t ctx)
{
    ulong degx, degy;

    FLINT_ASSERT(A->length > 0);

    degx = A->exps[0];
    degy = nmod_poly_degree((A->coeffs + 0)->coeffs + 0);

    FLINT_ASSERT(FLINT_BIT_COUNT(degx) < FLINT_BITS/2);
    FLINT_ASSERT(FLINT_BIT_COUNT(degy) < FLINT_BITS/2);

    return (degx << (FLINT_BITS/2)) + degy;
}


void fmpz_mpoly_eval_fmpz_mod(
    fmpz_t eval,
    const fmpz_mpoly_t A,
    const fmpz * alpha,
    const fmpz_mpoly_ctx_t ctx,
    const fmpz_mod_ctx_t fpctx)
{
    slong i, j;
    slong offset, shift;
    slong N = mpoly_words_per_exp_sp(A->bits, ctx->minfo);
    slong nvars = ctx->minfo->nvars;
    ulong * Aexp;
    slong * LUToffset;
    ulong * LUTmask;
    fmpz * LUTvalue;
    slong LUTlen;
    fmpz_t xpoweval;
    ulong * inputexpmask;
    TMP_INIT;

    TMP_START;

    fmpz_init(xpoweval);

    LUToffset = (slong *) TMP_ALLOC(N*FLINT_BITS*sizeof(slong));
    LUTmask   = (ulong *) TMP_ALLOC(N*FLINT_BITS*sizeof(ulong));
    LUTvalue  = (fmpz *) TMP_ALLOC(N*FLINT_BITS*sizeof(fmpz));
    for (i = 0; i < N*FLINT_BITS; i++)
    {
        fmpz_init(LUTvalue + i);
    }

    Aexp = A->exps;

    inputexpmask = (ulong *) TMP_ALLOC(N*sizeof(ulong));
    mpoly_monomial_zero(inputexpmask, N);
    for (i = 0; i < A->length; i++)
    {
        for (j = 0; j < N; j++)
        {
            inputexpmask[j] |= (Aexp + N*i)[j];
        }
    }

    LUTlen = 0;
    for (j = nvars - 1; j >= 0; j--)
    {
        mpoly_gen_offset_shift_sp(&offset, &shift, j, A->bits, ctx->minfo);

        fmpz_set(xpoweval, alpha + j); /* xpoweval = alpha[i]^(2^i) */
        for (i = 0; i < A->bits; i++)
        {
            LUToffset[LUTlen] = offset;
            LUTmask[LUTlen] = (UWORD(1) << (shift + i));
            fmpz_set(LUTvalue + LUTlen, xpoweval);
            if ((inputexpmask[offset] & LUTmask[LUTlen]) != 0)
            {
                LUTlen++;
            }
            fmpz_mod_mul(xpoweval, xpoweval, xpoweval, fpctx);
        }
    }
    FLINT_ASSERT(LUTlen < N*FLINT_BITS);

    fmpz_zero(eval);
    for (i = 0; i < A->length; i++)
    {
        fmpz_mod(xpoweval, A->coeffs + i, fpctx->p);
        for (j = 0; j < LUTlen; j++)
        {
            if (((Aexp + N*i)[LUToffset[j]] & LUTmask[j]) != 0)
            {
                fmpz_mod_mul(xpoweval, xpoweval, LUTvalue + j, fpctx);
            }
        }
        fmpz_mod_add(eval, eval, xpoweval, fpctx);
    }

    fmpz_clear(xpoweval);
    for (i = 0; i < N*FLINT_BITS; i++)
    {
        fmpz_clear(LUTvalue + i);
    }
    TMP_END;

/*
    printf("fmpz_mpoly_eval_fmpz_mod returning\n");
*/
}

/* take coeffs from [0, p) to (-p/2, p/2]
*/
void fmpz_mpolyu_symmetrize_coeffs(
    fmpz_mpolyu_t A,
    const fmpz_mpoly_ctx_t ctx,
    const fmpz_mod_ctx_t fpctx)
{
    slong i, j;
    fmpz_mpoly_struct * Ac;

    for (i = 0; i < A->length; i++)
    {
        Ac = A->coeffs + i;
        for (j = 0; j < Ac->length; j++)
        {
            fmpz_mods(Ac->coeffs + j, Ac->coeffs + j, fpctx->p);
        }
    }
}


void fmpz_mpolyuu_eval_fmpz_mod(
    fmpz_mod_mpolyun_t E,
    const fmpz_mpolyu_t A,
    const fmpz * alpha,
    const fmpz_mpoly_ctx_t ctx,
    const fmpz_mod_ctx_t fpctx)
{
    slong N = mpoly_words_per_exp_sp(A->bits, ctx->minfo);
    slong i;
    slong xexp, yexp;
    fmpz_t eval;

    FLINT_ASSERT(E->bits == A->bits);

    fmpz_init(eval);

    E->length = 0;
    for (i = 0; i < A->length; i++)
    {
        fmpz_mpoly_eval_fmpz_mod(eval, A->coeffs + i, alpha, ctx, fpctx);
        if (fmpz_is_zero(eval))
        {
            continue;
        }

        xexp = A->exps[i] >> (FLINT_BITS/2);
        yexp = A->exps[i] & (-UWORD(1) >> (FLINT_BITS - FLINT_BITS/2));

/*
printf("eval: "); fmpz_print(eval); printf("\n");
flint_printf("exp: (%wd, %wd)\n", xexp, yexp);
*/

        if (E->length > 0 && E->exps[E->length - 1] == xexp)
        {
            fmpz_mod_poly_set_coeff_fmpz(E->coeffs[E->length-1].coeffs + 0, yexp, eval);
        }
        else
        {
            fmpz_mod_mpolyun_fit_length(E, E->length + 1, ctx, fpctx);
            fmpz_mod_mpolyn_fit_length(E->coeffs + E->length, 1, ctx, fpctx);
            (E->coeffs + E->length)->length = 1;
            mpoly_monomial_zero((E->coeffs + E->length)->exps, N);
            fmpz_mod_poly_zero((E->coeffs + E->length)->coeffs + 0);
            fmpz_mod_poly_set_coeff_fmpz((E->coeffs + E->length)->coeffs + 0, yexp, eval);
            E->exps[E->length] = xexp;
            E->length++;
        }
    }
/*
flint_printf("Elength: %wd\n", E->length);
*/
    fmpz_clear(eval);
}




int nmod_mpolyn_equal(
    const nmod_mpolyn_t A,
    const nmod_mpolyn_t B,
    const nmod_mpoly_ctx_t ctx)
{
    slong N = mpoly_words_per_exp(A->bits, ctx->minfo);
    slong i;

    FLINT_ASSERT(A->bits == B->bits);

    if (A->length != B->length)
    {
        return 0;
    }
    for (i = 0; i < A->length; i++)
    {
        if (!mpoly_monomial_equal(A->exps + N*i, B->exps + N*i, N))
        {
            return 0;
        }
        if (!nmod_poly_equal(A->coeffs + i, B->coeffs + i))
        {
            return 0;
        }
    }
    return 1;
}

int fmpz_mod_mpolyn_equal(
    const fmpz_mod_mpolyn_t A,
    const fmpz_mod_mpolyn_t B,
    const fmpz_mpoly_ctx_t ctx,
    const fmpz_mod_ctx_t fpctx)
{
    slong N = mpoly_words_per_exp(A->bits, ctx->minfo);
    slong i;

    FLINT_ASSERT(A->bits == B->bits);

    if (A->length != B->length)
    {
        return 0;
    }
    for (i = 0; i < A->length; i++)
    {
        if (!mpoly_monomial_equal(A->exps + N*i, B->exps + N*i, N))
        {
            return 0;
        }
        if (!fmpz_mod_poly_equal(A->coeffs + i, B->coeffs + i))
        {
            return 0;
        }
    }
    return 1;
}



int nmod_mpolyun_equal(
    const nmod_mpolyun_t A,
    const nmod_mpolyun_t B,
    const nmod_mpoly_ctx_t ctx)
{
    slong i;

    FLINT_ASSERT(A->bits == B->bits);

    if (A->length != B->length)
    {
        return 0;
    }
    for (i = 0; i < A->length; i++)
    {
        if (A->exps[i] != B->exps[i])
        {
            return 0;
        }
        if (!nmod_mpolyn_equal(A->coeffs + i, B->coeffs + i, ctx))
        {
            return 0;
        }
    }
    return 1;
}


int fmpz_mod_mpolyun_equal(
    const fmpz_mod_mpolyun_t A,
    const fmpz_mod_mpolyun_t B,
    const fmpz_mpoly_ctx_t ctx,
    const fmpz_mod_ctx_t fpctx)
{
    slong i;

    FLINT_ASSERT(A->bits == B->bits);

    if (A->length != B->length)
    {
        return 0;
    }
    for (i = 0; i < A->length; i++)
    {
        if (A->exps[i] != B->exps[i])
        {
            return 0;
        }
        if (!fmpz_mod_mpolyn_equal(A->coeffs + i, B->coeffs + i, ctx, fpctx))
        {
            return 0;
        }
    }
    return 1;
}







mp_limb_t fmpz_mpoly_eval_nmod(
    const nmodf_ctx_t fpctx,
    const fmpz_mpoly_t A,
    const mp_limb_t * alpha,
    const fmpz_mpoly_ctx_t ctx)
{
    mp_limb_t eval;
    slong i, j;
    slong offset, shift;
    slong N = mpoly_words_per_exp_sp(A->bits, ctx->minfo);
    slong nvars = ctx->minfo->nvars;
    ulong * Aexp;
    slong * LUToffset;
    ulong * LUTmask;
    mp_limb_t * LUTvalue;
    slong LUTlen;
    mp_limb_t xpoweval;
    ulong * inputexpmask;
    TMP_INIT;

    FLINT_ASSERT(A->bits <= FLINT_BITS);

    TMP_START;

    LUToffset = (slong *) TMP_ALLOC(N*FLINT_BITS*sizeof(slong));
    LUTmask   = (ulong *) TMP_ALLOC(N*FLINT_BITS*sizeof(ulong));
    LUTvalue  = (mp_limb_t *) TMP_ALLOC(N*FLINT_BITS*sizeof(mp_limb_t));

    Aexp = A->exps;

    inputexpmask = (ulong *) TMP_ALLOC(N*sizeof(ulong));
    mpoly_monomial_zero(inputexpmask, N);
    for (i = 0; i < A->length; i++)
    {
        for (j = 0; j < N; j++)
        {
            inputexpmask[j] |= (Aexp + N*i)[j];
        }
    }

    LUTlen = 0;
    for (j = nvars - 1; j >= 0; j--)
    {
        mpoly_gen_offset_shift_sp(&offset, &shift, j, A->bits, ctx->minfo);

        xpoweval = alpha[j]; /* xpoweval = alpha[i]^(2^i) */
        for (i = 0; i < A->bits; i++)
        {
            LUToffset[LUTlen] = offset;
            LUTmask[LUTlen] = (UWORD(1) << (shift + i));
            LUTvalue[LUTlen] = xpoweval;
            if ((inputexpmask[offset] & LUTmask[LUTlen]) != 0)
            {
                LUTlen++;
            }
            xpoweval = nmod_mul(xpoweval, xpoweval, fpctx->mod);
        }
    }
    FLINT_ASSERT(LUTlen < N*FLINT_BITS);

    eval = 0;
    for (i = 0; i < A->length; i++)
    {
        xpoweval = fmpz_fdiv_ui(A->coeffs + i, fpctx->mod.n);
        for (j = 0; j < LUTlen; j++)
        {
            if (((Aexp + N*i)[LUToffset[j]] & LUTmask[j]) != 0)
            {
                xpoweval = nmod_mul(xpoweval, LUTvalue[j], fpctx->mod);
            }
        }
        eval = nmod_add(eval, xpoweval, fpctx->mod);
    }

    TMP_END;
/*
    printf("fmpz_mpoly_eval_fmpz_mod returning\n");
*/
    return eval;
}




void fmpz_mpolyuu_eval_nmod(
    nmod_mpolyun_t E,
    const nmod_mpoly_ctx_t ctx_sp,
    const fmpz_mpolyu_t A,
    const mp_limb_t * alpha,
    const fmpz_mpoly_ctx_t ctx)
{
    slong N = mpoly_words_per_exp_sp(A->bits, ctx->minfo);
    slong i;
    slong xexp, yexp;
    mp_limb_t eval;

    FLINT_ASSERT(E->bits == A->bits);

    E->length = 0;
    for (i = 0; i < A->length; i++)
    {
        eval = fmpz_mpoly_eval_nmod(ctx_sp->ffinfo, A->coeffs + i, alpha, ctx);
        if (eval == 0)
        {
            continue;
        }

        xexp = A->exps[i] >> (FLINT_BITS/2);
        yexp = A->exps[i] & (-UWORD(1) >> (FLINT_BITS - FLINT_BITS/2));

        if (E->length > 0 && E->exps[E->length - 1] == xexp)
        {
            nmod_poly_set_coeff_ui(E->coeffs[E->length-1].coeffs + 0, yexp, eval);
        }
        else
        {
            nmod_mpolyun_fit_length(E, E->length + 1, ctx_sp);
            nmod_mpolyn_fit_length(E->coeffs + E->length, 1, ctx_sp);
            (E->coeffs + E->length)->length = 1;
            mpoly_monomial_zero((E->coeffs + E->length)->exps, N);
            nmod_poly_zero((E->coeffs + E->length)->coeffs + 0);
            nmod_poly_set_coeff_ui((E->coeffs + E->length)->coeffs + 0, yexp, eval);
            E->exps[E->length] = xexp;
            E->length++;
        }
    }
}




typedef struct {
    slong mlength;
    slong malloc;
    mp_limb_t * coeffs;
    mp_limb_t * monomials;
    slong ealloc;
    mp_limb_t * evals;
} nmod_zip_struct;
typedef nmod_zip_struct nmod_zip_t[1];


void nmod_zip_init(nmod_zip_t Z)
{
    Z->mlength = 0;
    Z->malloc = 0;
    Z->coeffs = NULL;
    Z->monomials = NULL;
    Z->ealloc = 0;
    Z->evals = NULL;
}

void nmod_zip_clear(nmod_zip_t Z)
{
    if (Z->coeffs)
        flint_free(Z->coeffs);
    if (Z->monomials)
        flint_free(Z->monomials);
    if (Z->evals)
        flint_free(Z->evals);
}

void nmod_zip_set_lengths(nmod_zip_t A, slong mlength, slong elength)
{
    slong old_malloc = A->malloc;
    slong new_malloc = FLINT_MAX(mlength, A->malloc + A->malloc/2);
    slong old_ealloc = A->ealloc;
    slong new_ealloc = FLINT_MAX(elength, A->ealloc + A->ealloc/2);

    FLINT_ASSERT(mlength > 0);
    FLINT_ASSERT(elength > 0);

    if (mlength > old_malloc)
    {
        if (old_malloc == 0)
        {
            A->coeffs    = (mp_limb_t *) flint_malloc(new_malloc*sizeof(mp_limb_t));
            A->monomials = (mp_limb_t *) flint_malloc(new_malloc*sizeof(mp_limb_t));
        }
        else
        {
            A->coeffs    = (mp_limb_t *) flint_realloc(A->coeffs,    new_malloc*sizeof(mp_limb_t));
            A->monomials = (mp_limb_t *) flint_realloc(A->monomials, new_malloc*sizeof(mp_limb_t));
        }
        A->malloc = new_malloc;
    }

    A->mlength = mlength;

    if (elength > old_ealloc)
    {
        if (old_ealloc == 0)
        {
            A->evals = (mp_limb_t *) flint_malloc(new_ealloc*sizeof(mp_limb_t));
        }
        else
        {
            A->evals = (mp_limb_t *) flint_realloc(A->evals, new_ealloc*sizeof(mp_limb_t));
        }
        A->ealloc = new_ealloc;
    }
}

void nmod_zip_print(const nmod_zip_t Z, slong elength)
{
    slong i;

    printf("m ");
    for (i = 0; i < Z->mlength; i++)
    {
        flint_printf("(%wu %wu) ", Z->coeffs[i], Z->monomials[i]);
    }
    printf("e ");
    for (i = 0; i < elength; i++)
    {
        flint_printf("%wu ", Z->evals[i]);
    }
}

typedef struct {
    nmod_zip_struct * coeffs;
    ulong * exps;
    slong length;
    slong alloc;
    slong pointcount;
} nmod_zip_mpolyu_struct;

typedef nmod_zip_mpolyu_struct nmod_zip_mpolyu_t[1];

void nmod_zip_mpolyu_init(nmod_zip_mpolyu_t Z)
{
    Z->alloc = 0;
    Z->length = 0;
    Z->exps = NULL;
    Z->coeffs = NULL;
    Z->pointcount = 0;
}

void nmod_zip_mpolyu_clear(nmod_zip_mpolyu_t Z)
{
    slong i;

    for (i = 0; i < Z->alloc; i++)
    {
        nmod_zip_clear(Z->coeffs + i);
    }

    if (Z->exps)
        flint_free(Z->exps);
    if (Z->coeffs)
        flint_free(Z->coeffs);
}

void nmod_zip_mpolyu_fit_length(
    nmod_zip_mpolyu_t A,
    slong length)
{
    slong i;
    slong old_alloc = A->alloc;
    slong new_alloc = FLINT_MAX(length, A->alloc + A->alloc/2);

    if (length > old_alloc)
    {
        if (old_alloc == 0)
        {
            A->exps = (ulong *) flint_malloc(new_alloc*sizeof(ulong));
            A->coeffs = (nmod_zip_struct *) flint_malloc(new_alloc*sizeof(nmod_zip_struct));
        }
        else
        {
            A->exps = (ulong *) flint_realloc(A->exps, new_alloc*sizeof(ulong));
            A->coeffs = (nmod_zip_struct *) flint_realloc(A->coeffs, new_alloc*sizeof(nmod_zip_struct));
        }

        for (i = old_alloc; i < new_alloc; i++)
        {
            nmod_zip_init(A->coeffs + i);
        }
        A->alloc = new_alloc;
    }
}

void nmod_zip_mpolyu_fit_poly(
    nmod_zip_mpolyu_t Z,
    fmpz_mpolyu_t H,
    slong eval_length)
{
    slong i;
    FLINT_ASSERT(H->length > 0);

    nmod_zip_mpolyu_fit_length(Z, H->length);

    for (i = 0; i < H->length; i++)
    {
        Z->exps[i] = H->exps[i];
        nmod_zip_set_lengths(Z->coeffs + i, (H->coeffs + i)->length, eval_length);
    }

    Z->length = H->length;
    Z->pointcount = 0;
}


void nmod_mpoly_set_skel(
    nmod_mpolyc_t S,
    const nmod_mpoly_ctx_t ctx_sp,
    const fmpz_mpoly_t A,
    const mp_limb_t * alpha,
    const fmpz_mpoly_ctx_t ctx);


void nmod_zip_mpolyu_set_skel(
    nmod_zip_mpolyu_t Z,
    const nmod_mpoly_ctx_t ctx_sp,
    const fmpz_mpolyu_t A,
    const mp_limb_t * alpha,
    const fmpz_mpoly_ctx_t ctx)
{
    slong i, j;

    nmod_mpolyc_t T;
    nmod_mpolyc_init(T);

    FLINT_ASSERT(Z->length == A->length);
    for (i = 0; i < Z->length; i++)
    {
        nmod_zip_struct * Zc = Z->coeffs + i;

        nmod_mpoly_set_skel(T, ctx_sp, A->coeffs + i, alpha, ctx);

        Z->exps[i] = A->exps[i];
        FLINT_ASSERT(Zc->mlength == T->length);
        for (j = 0; j < Zc->mlength; j++)
        {
            Zc->coeffs[j] = 0;
            Zc->monomials[j] = T->coeffs[j];
        }
    }
    Z->pointcount = 0;

    nmod_mpolyc_clear(T);
}

void nmod_zip_mpolyuu_print(const nmod_zip_mpolyu_t A)
{
    slong i;
    flint_printf("0");
    for (i = 0; i < A->length; i++)
    {
        flint_printf(" + [");
        nmod_zip_print(A->coeffs + i, A->pointcount);
        flint_printf("]*X^%wd*Y^%wd", A->exps[i] >> (FLINT_BITS/2), A->exps[i] & ((-UWORD(1)) >> (FLINT_BITS/2)));
    }
}

int nmod_zip_mpolyuu_add_point(
    nmod_zip_mpolyu_t L,
    const nmod_mpolyun_t A)
{
    slong Alen = A->length;
    nmod_mpolyn_struct * Acoeff = A->coeffs;
    slong Li, Ai, ai;
    ulong Aexp;
    nmod_zip_struct * Lcoeff;
    slong Llen;
    ulong * Lexp;
    slong pointcount = L->pointcount;

    Lcoeff = L->coeffs;
    Llen = L->length;
    Lexp = L->exps;

    Li = 0;
    Ai = ai = 0;
    Aexp = 0;
    if (Ai < Alen)
    {
        ai = nmod_poly_degree((A->coeffs + Ai)->coeffs + 0);
        Aexp = (A->exps[Ai] << (FLINT_BITS/2)) + ai;
    }

    for (Li = 0; Li < Llen; Li++)
    {
        nmod_zip_struct * Lc = Lcoeff + Li;

        if (Ai < Alen && Lexp[Li] == Aexp)
        {
            /* L present A present */
            Lc->evals[pointcount] = (Acoeff + Ai)->coeffs[0].coeffs[ai];
            do {
                ai--;
            } while (ai >= 0 && (0 == (Acoeff + Ai)->coeffs[0].coeffs[ai]));
            if (ai < 0)
            {
                Ai++;
                if (Ai < Alen)
                {
                    ai = nmod_poly_degree((A->coeffs + Ai)->coeffs + 0);        
                    Aexp = (A->exps[Ai] << (FLINT_BITS/2)) + ai;
                }
            }
            if (Ai < Alen)
            {
                Aexp = (A->exps[Ai] << (FLINT_BITS/2)) + ai;
            }
        }
        else if (Ai >= Alen || Lexp[Li] > Aexp)
        {
            /* L present A missing */
            Lc->evals[pointcount] = 0;
        }
        else
        {
            /* L missing A present */
            return 0;
        }
    }

    L->pointcount = pointcount + 1;
    return 1;
}




void nmod_mpoly_set_skel(
    nmod_mpolyc_t S,
    const nmod_mpoly_ctx_t ctx_sp,
    const fmpz_mpoly_t A,
    const mp_limb_t * alpha,
    const fmpz_mpoly_ctx_t ctx)
{
    slong i, j;
    slong offset, shift;
    slong N = mpoly_words_per_exp_sp(A->bits, ctx->minfo);
    slong nvars = ctx->minfo->nvars;
    ulong * Aexp;
    slong * LUToffset;
    ulong * LUTmask;
    mp_limb_t * LUTvalue;
    slong LUTlen;
    mp_limb_t xpoweval;
    ulong * inputexpmask;
    TMP_INIT;

    FLINT_ASSERT(A->bits <= FLINT_BITS);

    TMP_START;

    LUToffset = (slong *) TMP_ALLOC(N*FLINT_BITS*sizeof(slong));
    LUTmask   = (ulong *) TMP_ALLOC(N*FLINT_BITS*sizeof(ulong));
    LUTvalue  = (mp_limb_t *) TMP_ALLOC(N*FLINT_BITS*sizeof(mp_limb_t));

    Aexp = A->exps;

    inputexpmask = (ulong *) TMP_ALLOC(N*sizeof(ulong));
    mpoly_monomial_zero(inputexpmask, N);
    for (i = 0; i < A->length; i++)
    {
        for (j = 0; j < N; j++)
        {
            inputexpmask[j] |= (Aexp + N*i)[j];
        }
    }

    LUTlen = 0;
    for (j = nvars - 1; j >= 0; j--)
    {
        mpoly_gen_offset_shift_sp(&offset, &shift, j, A->bits, ctx->minfo);

        xpoweval = alpha[j]; /* xpoweval = alpha[i]^(2^i) */
        for (i = 0; i < A->bits; i++)
        {
            LUToffset[LUTlen] = offset;
            LUTmask[LUTlen] = (UWORD(1) << (shift + i));
            LUTvalue[LUTlen] = xpoweval;
            if ((inputexpmask[offset] & LUTmask[LUTlen]) != 0)
            {
                LUTlen++;
            }
            xpoweval = nmod_mul(xpoweval, xpoweval, ctx_sp->ffinfo->mod);
        }
    }
    FLINT_ASSERT(LUTlen < N*FLINT_BITS);

    nmod_mpolyc_fit_length(S, A->length);
    for (i = 0; i < A->length; i++)
    {
        xpoweval = 1;
        for (j = 0; j < LUTlen; j++)
        {
            if (((Aexp + N*i)[LUToffset[j]] & LUTmask[j]) != 0)
            {
                xpoweval = nmod_mul(xpoweval, LUTvalue[j], ctx_sp->ffinfo->mod);
            }
        }
        S->coeffs[i] = xpoweval;
    }
    S->length = A->length;

    TMP_END;
}

void nmod_mpolyu_set_skel(
    nmod_mpolycu_t S,
    const nmod_mpoly_ctx_t ctx_sp,
    const fmpz_mpolyu_t A,
    const mp_limb_t * alpha,
    const fmpz_mpoly_ctx_t ctx)
{
    slong i;
    nmod_mpolycu_fit_length(S, A->length);
    for (i = 0; i < A->length; i++)
    {
        nmod_mpoly_set_skel(S->coeffs + i, ctx_sp, A->coeffs + i, alpha, ctx);
    }
    S->length = A->length;
}



/* M = S */
void nmod_mpoly_red_skel(
    nmod_mpolyc_t Ared,
    const fmpz_mpoly_t A,
    const nmodf_ctx_t fpctx)
{
    nmod_mpolyc_fit_length(Ared, A->length);
    Ared->length = A->length;
    _fmpz_vec_get_nmod_vec(Ared->coeffs, A->coeffs, A->length, fpctx->mod);
}

void nmod_mpolyu_red_skel(
    nmod_mpolycu_t Ared,
    const fmpz_mpolyu_t A,
    const nmodf_ctx_t fpctx)
{
    slong i;
    nmod_mpolycu_fit_length(Ared, A->length);
    Ared->length = A->length;
    for (i = 0; i < A->length; i++)
    {
        nmod_mpoly_red_skel(Ared->coeffs + i, A->coeffs + i, fpctx);
    }
}



void nmod_mpoly_copy_skel(nmod_mpolyc_t M, const nmod_mpolyc_t S)
{
    slong i;
    nmod_mpolyc_fit_length(M, S->length);
    M->length = S->length;
    for (i = 0; i < S->length; i++)
    {
        M->coeffs[i] = S->coeffs[i];
    }
}

void nmod_mpolyu_copy_skel(nmod_mpolycu_t M, const nmod_mpolycu_t S)
{
    slong i;
    nmod_mpolycu_fit_length(M, S->length);
    M->length = S->length;
    for (i = 0; i < S->length; i++)
    {
        nmod_mpoly_copy_skel(M->coeffs + i, S->coeffs + i);
    }
}

mp_limb_t nmod_mpoly_use_skel_mul(
    const nmod_mpolyc_t Ared,
    nmod_mpolyc_t Acur,
    const nmod_mpolyc_t Ainc,
    const nmod_mpoly_ctx_t ctx_sp)
{
    slong i;
    mp_limb_t V, V0, V1, V2, p1, p0;
    FLINT_ASSERT(Ared->length == Acur->length);
    FLINT_ASSERT(Ared->length == Ainc->length);

    V0 = V1 = V2 = 0;
    for (i = 0; i < Ared->length; i++)
    {
        umul_ppmm(p1, p0, Ared->coeffs[i], Acur->coeffs[i]);
        add_sssaaaaaa(V2, V1, V0, V2, V1, V0, WORD(0), p1, p0);
        Acur->coeffs[i] = nmod_mul(Acur->coeffs[i], Ainc->coeffs[i], ctx_sp->ffinfo->mod);
    }
    
    NMOD_RED3(V, V2, V1, V0, ctx_sp->ffinfo->mod);
    return V;
}


void nmod_mpolyuu_use_skel_mul(
    nmod_mpolyun_t E,
    const fmpz_mpolyu_t A,
    const nmod_mpolycu_t Ared,
    nmod_mpolycu_t Acur,
    const nmod_mpolycu_t Ainc,
    const nmod_mpoly_ctx_t ctx_sp)
{
    slong xexp, yexp;
    slong i;
    slong N = mpoly_words_per_exp_sp(E->bits, ctx_sp->minfo);
    mp_limb_t eval;

    FLINT_ASSERT(A->length == Ared->length);
    FLINT_ASSERT(A->length == Acur->length);
    FLINT_ASSERT(A->length == Ainc->length);

    E->length = 0;
    for (i = 0; i < A->length; i++)
    {
        eval = nmod_mpoly_use_skel_mul(Ared->coeffs + i, Acur->coeffs + i, Ainc->coeffs + i, ctx_sp);
        if (eval == 0)
        {
            continue;
        }

        xexp = A->exps[i] >> (FLINT_BITS/2);
        yexp = A->exps[i] & ((-UWORD(1)) >> (FLINT_BITS - FLINT_BITS/2));
/*
flint_printf("xexp %wd  yexp %wd  eval "); fmpz_print(eval); printf("\n");
*/
        if (E->length > 0 && E->exps[E->length - 1] == xexp)
        {
            nmod_poly_set_coeff_ui(E->coeffs[E->length-1].coeffs + 0, yexp, eval);
        }
        else
        {
            nmod_mpolyun_fit_length(E, E->length + 1, ctx_sp);
            nmod_mpolyn_fit_length(E->coeffs + E->length, 1, ctx_sp);
            (E->coeffs + E->length)->length = 1;
            mpoly_monomial_zero((E->coeffs + E->length)->exps, N);
            nmod_poly_zero((E->coeffs + E->length)->coeffs + 0);
            nmod_poly_set_coeff_ui((E->coeffs + E->length)->coeffs + 0, yexp, eval);
            E->exps[E->length] = xexp;
            E->length++;
        }
    }
}


void nmod_mpolyun_scalar_mul_nmod(nmod_mpolyun_t A, mp_limb_t c, const nmod_mpoly_ctx_t ctx);


typedef enum {
    nmod_zip_find_coeffs_good,
    nmod_zip_find_coeffs_no_match,
    nmod_zip_find_coeffs_non_invertible
} nmod_zip_find_coeffs_ret_t;


nmod_zip_find_coeffs_ret_t nmod_zip_find_coeffs(
    nmod_zip_t Z,
    nmod_poly_t master,
    slong pointcount,
    const nmodf_ctx_t ffinfo)
{
    slong i, j;
    mp_limb_t V, V0, V1, V2, T, S, r, p0, p1;

    FLINT_ASSERT(pointcount >= Z->mlength);

    nmod_poly_product_roots_nmod_vec(master, Z->monomials, Z->mlength);

    for (i = 0; i < Z->mlength; i++)
    {
        /* coeffs[i] is (coeffs(P).values)/P(roots[i]) =: V/S
            where P(x) = master(x)/(x-roots[i])     */
        V0 = V1 = V2 = T = S = 0;
        r = Z->monomials[i];
        for (j = Z->mlength; j > 0; j--)
        {
            T = nmod_add(nmod_mul(r, T, ffinfo->mod), master->coeffs[j], ffinfo->mod);
            S = nmod_add(nmod_mul(r, S, ffinfo->mod), T, ffinfo->mod);
            umul_ppmm(p1, p0, Z->evals[j - 1], T);
            add_sssaaaaaa(V2, V1, V0, V2, V1, V0, WORD(0), p1, p0);
        }
        /* roots[i] should be a root of master */
        FLINT_ASSERT(nmod_add(nmod_mul(r, T, ffinfo->mod), master->coeffs[0], ffinfo->mod) == 0);
        NMOD_RED3(V, V2, V1, V0, ffinfo->mod);
        S = nmod_mul(S, r, ffinfo->mod); /* shift is one */
        if (S == 0)
        {
            return nmod_zip_find_coeffs_non_invertible;
        }
        Z->coeffs[i] = nmod_mul(V, nmod_inv(S, ffinfo->mod), ffinfo->mod);
    }

    /* use the coefficients of master as temp work space */
    for (i = 0; i < Z->mlength; i++)
    {
        master->coeffs[i] = nmod_pow_ui(Z->monomials[i], Z->mlength, ffinfo->mod);
    }

    /* check that the remaining points match */
    for (i = Z->mlength; i < pointcount; i++)
    {
        V0 = V1 = V2 = S = 0;
        for (j = 0; j < Z->mlength; j++)
        {
            master->coeffs[j] = nmod_mul(master->coeffs[j], Z->monomials[j], ffinfo->mod);
            umul_ppmm(p1, p0, Z->coeffs[j], master->coeffs[j]);
            add_sssaaaaaa(V2, V1, V0, V2, V1, V0, WORD(0), p1, p0);
        }
        NMOD_RED3(V, V2, V1, V0, ffinfo->mod);
        if (V != Z->evals[i])
        {
            return nmod_zip_find_coeffs_no_match;
        }
    }

    return nmod_zip_find_coeffs_good;
}

nmod_zip_find_coeffs_ret_t nmod_mpolyu_zip_find_coeffs(
    nmod_zip_mpolyu_t Z,
    const nmod_mpoly_ctx_t ctx_sp)
{
    slong i;
    nmod_zip_find_coeffs_ret_t ret;
    nmod_poly_t T;

    nmod_poly_init_mod(T, ctx_sp->ffinfo->mod);

    for (i = 0; i < Z->length; i++)
    {
        ret = nmod_zip_find_coeffs(Z->coeffs + i, T, Z->pointcount, ctx_sp->ffinfo);
        if (ret != nmod_zip_find_coeffs_good)
        {
            goto cleanup;
        }
    }

    ret = nmod_zip_find_coeffs_good;

cleanup:

    nmod_poly_clear(T);

    return ret;
}


int fmpz_mpolyu_addinterp_zip(
    fmpz_mpolyu_t H,
    const fmpz_t Hmodulus,
    const nmod_zip_mpolyu_t Z,
    const nmodf_ctx_t ffinfo)
{
    int changed = 0;
    slong i, j;
    fmpz_t t;

    fmpz_init(t);

    FLINT_ASSERT(H->length == Z->length);
    for (i = 0; i < H->length; i++)
    {
        fmpz_mpoly_struct * Hc = H->coeffs + i;
        nmod_zip_struct * Zc = Z->coeffs + i;

        FLINT_ASSERT(Hc->length == Zc->mlength);
        for (j = 0; j < Hc->length; j++)
        {
            fmpz_CRT_ui(t, Hc->coeffs + j, Hmodulus, Zc->coeffs[j], ffinfo->mod.n, 1);
            changed |= !fmpz_equal(t, Hc->coeffs + j);
            fmpz_swap(t, Hc->coeffs + j);
        }
    }

    fmpz_clear(t);
    return changed;
}



int fmpz_mpoly_repack_bits_inplace(
    fmpz_mpoly_t A,
    mp_bitcnt_t Abits,
    const fmpz_mpoly_ctx_t ctx)
{
    int success;
    ulong * texps;

    if (A->bits == Abits)
    {
        return 1;
    }

    texps = (ulong *) flint_malloc(A->alloc*mpoly_words_per_exp(Abits, ctx->minfo)*sizeof(ulong));
    success = mpoly_repack_monomials(texps, Abits, A->exps, A->bits, A->length, ctx->minfo);
    if (success)
    {
        ulong * t = A->exps;
        A->exps = texps;
        texps = t;
        A->bits = Abits;
    }
    flint_free(texps);
    return success;
}


int fmpz_mpolyu_repack_bits(
    fmpz_mpolyu_t A,
    mp_bitcnt_t Abits,
    const fmpz_mpoly_ctx_t ctx)
{
    mp_bitcnt_t org_bits = A->bits;
    int success;
    slong i, j;

    for (i = 0; i < A->length; i++)
    {
        FLINT_ASSERT((A->coeffs + i)->bits == A->bits);
        success = fmpz_mpoly_repack_bits_inplace(A->coeffs + i, Abits, ctx);
        if (!success)
        {
            /* repack changed coeffs */
            for (j = 0; j < i; j++)
            {
                success = fmpz_mpoly_repack_bits_inplace(A->coeffs + j, org_bits, ctx);
                FLINT_ASSERT(success);
            }
            return 0;
        }
    }
    return 1;
}


/* fit bits is not the best fxn to accomplish this */
void fmpz_mpoly_set_bits(
    fmpz_mpoly_t A,
    mp_bitcnt_t Abits,
    const fmpz_mpoly_ctx_t ctx)
{
    fmpz_mpoly_fit_bits(A, Abits, ctx);
    A->bits = Abits;
}

void fmpz_mpolyu_set_bits(
    fmpz_mpolyu_t A,
    mp_bitcnt_t Abits,
    const fmpz_mpoly_ctx_t ctx)
{
    slong i;
    for (i = 0; i < A->length; i++)
    {
        fmpz_mpoly_set_bits(A->coeffs + i, Abits, ctx);
    }
    A->bits = Abits;
}




/* A = D - B*C, D may be modified if saveD == 0 */
slong _fmpz_mpoly_mulsub(
                fmpz ** A_coeff, ulong ** A_exp, slong * A_alloc,
                 const fmpz * Dcoeff, const ulong * Dexp, slong Dlen, int saveD,
                 const fmpz * Bcoeff, const ulong * Bexp, slong Blen,
                 const fmpz * Ccoeff, const ulong * Cexp, slong Clen,
                              mp_bitcnt_t bits, slong N, const ulong * cmpmask)
{
    slong i, j;
    slong next_loc;
    slong heap_len = 2; /* heap zero index unused */
    mpoly_heap_s * heap;
    mpoly_heap_t * chain;
    slong * store, * store_base;
    mpoly_heap_t * x;
    slong Di;
    slong Alen;
    slong Aalloc = *A_alloc;
    fmpz * Acoeff = *A_coeff;
    ulong * Aexp = *A_exp;
    ulong * exp, * exps;
    ulong ** exp_list;
    slong exp_next;
    slong * hind;
    int small;
    ulong acc_sm[3], pp0, pp1;
    TMP_INIT;

   /* if exponent vectors fit in single word, call special version */
/*
   if (N == 1)
      return _fmpz_mpoly_mul_johnson1(poly1, exp1, alloc,
                             poly2, exp2, len2, poly3, exp3, len3, cmpmask[0]);
*/

    TMP_START;

    /* whether input coeffs are small, thus accumulation fit in three words */
    small =   _fmpz_mpoly_fits_small(Bcoeff, Blen)
           && _fmpz_mpoly_fits_small(Ccoeff, Clen);

    next_loc = Blen + 4; /* something bigger than heap can ever be */
    heap = (mpoly_heap_s *) TMP_ALLOC((Blen + 1)*sizeof(mpoly_heap_s));
    chain = (mpoly_heap_t *) TMP_ALLOC(Blen*sizeof(mpoly_heap_t));
    store = store_base = (slong *) TMP_ALLOC(2*Blen*sizeof(slong));
    exps = (ulong *) TMP_ALLOC(Blen*N*sizeof(ulong));
    exp_list = (ulong **) TMP_ALLOC(Blen*sizeof(ulong *));
    for (i = 0; i < Blen; i++)
        exp_list[i] = exps + i*N;

    /* space for heap indices */
    hind = (slong *) TMP_ALLOC(Blen*sizeof(slong));
    for (i = 0; i < Blen; i++)
        hind[i] = 1;

    /* start with no heap nodes and no exponent vectors in use */
    exp_next = 0;

    /* put (0, 0, Bexp[0] + Cexp[0]) on heap */
    x = chain + 0;
    x->i = 0;
    x->j = 0;
    x->next = NULL;

    heap[1].next = x;
    heap[1].exp = exp_list[exp_next++];

    if (bits <= FLINT_BITS)
        mpoly_monomial_add(heap[1].exp, Bexp + N*0, Cexp + N*0, N);
    else
        mpoly_monomial_add_mp(heap[1].exp, Bexp + N*0, Cexp + N*0, N);

    hind[0] = 2*1 + 0;

    Alen = 0;
    Di = 0;
    while (heap_len > 1)
    {
        exp = heap[1].exp;

        while (Di < Dlen && mpoly_monomial_gt(Dexp + N*Di, exp, N, cmpmask))
        {
            _fmpz_mpoly_fit_length(&Acoeff, &Aexp, &Aalloc, Alen + 1, N);
            mpoly_monomial_set(Aexp + N*Alen, Dexp + N*Di, N);
            if (saveD)
                fmpz_set(Acoeff + Alen, Dcoeff + Di);
            else
                fmpz_swap(Acoeff + Alen, (fmpz *)(Dcoeff + Di));
            Alen++;
            Di++;
        }

        _fmpz_mpoly_fit_length(&Acoeff, &Aexp, &Aalloc, Alen + 1, N);

        mpoly_monomial_set(Aexp + N*Alen, exp, N);

        if (small)
        {
            acc_sm[0] = acc_sm[1] = acc_sm[2] = 0;

            do
            {
                exp_list[--exp_next] = heap[1].exp;
                x = _mpoly_heap_pop(heap, &heap_len, N, cmpmask);

                hind[x->i] |= WORD(1);
                *store++ = x->i;
                *store++ = x->j;
                FLINT_ASSERT(!COEFF_IS_MPZ(Bcoeff[x->i]));
                FLINT_ASSERT(!COEFF_IS_MPZ(Ccoeff[x->j]));
                smul_ppmm(pp1, pp0, Bcoeff[x->i], Ccoeff[x->j]);
                sub_dddmmmsss(acc_sm[2], acc_sm[1], acc_sm[0],
                              acc_sm[2], acc_sm[1], acc_sm[0],
                              FLINT_SIGN_EXT(pp1), pp1, pp0);

                while ((x = x->next) != NULL)
                {
                    hind[x->i] |= WORD(1);
                    *store++ = x->i;
                    *store++ = x->j;
                    FLINT_ASSERT(!COEFF_IS_MPZ(Bcoeff[x->i]));
                    FLINT_ASSERT(!COEFF_IS_MPZ(Ccoeff[x->j]));
                    smul_ppmm(pp1, pp0, Bcoeff[x->i], Ccoeff[x->j]);
                    sub_dddmmmsss(acc_sm[2], acc_sm[1], acc_sm[0],
                                  acc_sm[2], acc_sm[1], acc_sm[0],
                                  FLINT_SIGN_EXT(pp1), pp1, pp0);
                }
            } while (heap_len > 1 && mpoly_monomial_equal(heap[1].exp, exp, N));

            fmpz_set_signed_uiuiui(Acoeff + Alen,
                                             acc_sm[2], acc_sm[1], acc_sm[0]);

            if (Di < Dlen && mpoly_monomial_equal(Dexp + N*Di, exp, N))
            {
                fmpz_add(Acoeff + Alen, Acoeff + Alen, Dcoeff + Di);
                Di++;
            }
        }
        else
        {
            if (Di < Dlen && mpoly_monomial_equal(Dexp + N*Di, exp, N))
            {
                fmpz_set(Acoeff + Alen, Dcoeff + Di);
                Di++;
            }
            else
            {
                fmpz_zero(Acoeff + Alen);
            }

            do
            {
                exp_list[--exp_next] = heap[1].exp;
                x = _mpoly_heap_pop(heap, &heap_len, N, cmpmask);

                hind[x->i] |= WORD(1);
                *store++ = x->i;
                *store++ = x->j;
                fmpz_submul(Acoeff + Alen, Bcoeff + x->i, Ccoeff + x->j);

                while ((x = x->next) != NULL)
                {
                    hind[x->i] |= WORD(1);
                    *store++ = x->i;
                    *store++ = x->j;
                    fmpz_submul(Acoeff + Alen, Bcoeff + x->i, Ccoeff + x->j);
                }
            } while (heap_len > 1 && mpoly_monomial_equal(heap[1].exp, exp, N));
        }

        Alen += !fmpz_is_zero(Acoeff + Alen);

        /* process nodes taken from the heap */
        while (store > store_base)
        {
            j = *--store;
            i = *--store;

            /* should we go right? */
            if (  (i + 1 < Blen)
                && (hind[i + 1] == 2*j + 1)
               )
            {
                x = chain + i + 1;
                x->i = i + 1;
                x->j = j;
                x->next = NULL;

                hind[x->i] = 2*(x->j+1) + 0;

                if (bits <= FLINT_BITS)
                    mpoly_monomial_add(exp_list[exp_next], Bexp + N*x->i,
                                                           Cexp + N*x->j, N);
                else
                    mpoly_monomial_add_mp(exp_list[exp_next], Bexp + N*x->i,
                                                              Cexp + N*x->j, N);

                exp_next += _mpoly_heap_insert(heap, exp_list[exp_next], x,
                                             &next_loc, &heap_len, N, cmpmask);
            }

            /* should we go up? */
            if (  (j + 1 < Clen)
               && ((hind[i] & 1) == 1)
               && ((i == 0) || (hind[i - 1] >= 2*(j + 2) + 1))
               )
            {
                x = chain + i;
                x->i = i;
                x->j = j + 1;
                x->next = NULL;

                hind[x->i] = 2*(x->j+1) + 0;

                if (bits <= FLINT_BITS)
                    mpoly_monomial_add(exp_list[exp_next], Bexp + N*x->i,
                                                           Cexp + N*x->j, N);
                else
                    mpoly_monomial_add_mp(exp_list[exp_next], Bexp + N*x->i,
                                                              Cexp + N*x->j, N);

                exp_next += _mpoly_heap_insert(heap, exp_list[exp_next], x,
                                             &next_loc, &heap_len, N, cmpmask);
            }
      }
   }

    FLINT_ASSERT(Di <= Dlen);
    _fmpz_mpoly_fit_length(&Acoeff, &Aexp, &Aalloc, Alen + Dlen - Di, N);
    if (saveD)
        _fmpz_vec_set(Acoeff + Alen, Dcoeff + Di, Dlen - Di);
    else
        _fmpz_vec_swap(Acoeff + Alen, (fmpz *)(Dcoeff + Di), Dlen - Di);
    mpoly_copy_monomials(Aexp + N*Alen, Dexp + N*Di, Dlen - Di, N);
    Alen += Dlen - Di;

    *A_coeff = Acoeff;
    *A_exp = Aexp;
    *A_alloc = Aalloc;

    TMP_END;

    return Alen;
}






int fmpz_mpolyuu_divides(
    fmpz_mpolyu_t Q,
    const fmpz_mpolyu_t A,
    const fmpz_mpolyu_t B,
    slong nmainvars,
    const fmpz_mpoly_ctx_t ctx)
{
    mp_bitcnt_t bits = A->bits;
    fmpz_mpoly_struct * Bcoeff = B->coeffs;
    slong Blen = B->length;
    ulong * Bexp = B->exps;
    fmpz_mpoly_struct * Acoeff = A->coeffs;
    slong Alen = A->length;
    ulong * Aexp = A->exps;
    fmpz_mpoly_struct * a, * b, * q;
    slong N;
    ulong * cmpmask;    /* cmp mask for lesser variables */
    fmpz_mpoly_t T, S;
    int success;
    ulong maskhi = 0;   /* main variables are in lex */
    int lt_divides;
    slong i, j, s;
    slong next_loc, heap_len;
    mpoly_heap1_s * heap;
    mpoly_heap_t * chain;
    slong * store, * store_base;
    mpoly_heap_t * x;
    slong * hind;
    ulong mask, exp, maxexp = Aexp[Alen - 1];
    TMP_INIT;

    FLINT_ASSERT(Alen > 0);
    FLINT_ASSERT(Blen > 0);
    FLINT_ASSERT(bits == A->bits);
    FLINT_ASSERT(bits == B->bits);
    FLINT_ASSERT(bits == Q->bits);

    TMP_START;

    N = mpoly_words_per_exp(bits, ctx->minfo);
    cmpmask = (ulong *) TMP_ALLOC(N*sizeof(ulong));
    mpoly_get_cmpmask(cmpmask, N, bits, ctx->minfo);

/*
printf("inside fmpz_mpolyuu_divides\n");
printf("A: "); fmpz_mpolyuu_print_pretty(A, NULL, nmainvars, ctx); printf("\n");
printf("B: "); fmpz_mpolyuu_print_pretty(B, NULL, nmainvars, ctx); printf("\n");
*/

    /* alloc array of heap nodes which can be chained together */
    next_loc = Blen + 4;   /* something bigger than heap can ever be */
    heap = (mpoly_heap1_s *) TMP_ALLOC((Blen + 1)*sizeof(mpoly_heap1_s));
    chain = (mpoly_heap_t *) TMP_ALLOC(Blen*sizeof(mpoly_heap_t));
    store = store_base = (slong *) TMP_ALLOC(2*Blen*sizeof(mpoly_heap_t *));

    /* space for flagged heap indicies */
    hind = (slong *) TMP_ALLOC(Blen*sizeof(slong));
    for (i = 0; i < B->length; i++)
        hind[i] = 1;

    /* mask with high bit set in each field of main exponent vector */
    mask = 0;
    for (i = 0; i < nmainvars; i++)
        mask = (mask << (FLINT_BITS/nmainvars)) + (UWORD(1) << (FLINT_BITS/nmainvars - 1));

    /* s is the number of terms * (latest quotient) we should put into heap */
    s = Blen;

    /* insert (-1, 0, exp2[0]) into heap */
    heap_len = 2;
    x = chain + 0;
    x->i = -WORD(1);
    x->j = 0;
    x->next = NULL;
    HEAP_ASSIGN(heap[1], Aexp[0], x);

    Q->length = 0;

    fmpz_mpoly_init3(T, 16, bits, ctx);
    fmpz_mpoly_init3(S, 16, bits, ctx);

    while (heap_len > 1)
    {
        exp = heap[1].exp;

        if (mpoly_monomial_overflows1(exp, mask))
            goto not_exact_division;

        fmpz_mpolyu_fit_length(Q, Q->length + 1, ctx);
        lt_divides = mpoly_monomial_divides1(Q->exps + Q->length, exp, Bexp[0], mask);

        T->length = 0;

        do
        {
            x = _mpoly_heap_pop1(heap, &heap_len, maskhi);
            do
            {
                *store++ = x->i;
                *store++ = x->j;
                if (x->i != -WORD(1))
                    hind[x->i] |= WORD(1);

                if (x->i == -WORD(1))
                {
                    a = Acoeff + x->j;
                    fmpz_mpoly_fit_length(S, T->length + a->length, ctx);
                    S->length = _fmpz_mpoly_add(
                                    S->coeffs, S->exps,
                                    T->coeffs, T->exps, T->length,
                                    a->coeffs, a->exps, a->length,
                                                              N, cmpmask);
                }
                else
                {
                    b = Bcoeff + x->i;
                    q = Q->coeffs + x->j;
                    S->length = _fmpz_mpoly_mulsub(
                                    &S->coeffs, &S->exps, &S->alloc,
                                    T->coeffs, T->exps, T->length, 0,
                                    b->coeffs, b->exps, b->length,
                                    q->coeffs, q->exps, q->length,
                                                         bits, N, cmpmask);
/*
        {
            fmpz_mpoly_ctx_t ctx;
            fmpz_mpoly_t SS;

            fmpz_mpoly_ctx_init(ctx, 2, ORD_LEX);


            SS->bits = bits;
            SS->coeffs = S->coeffs;
            SS->exps = S->exps;
            SS->length = S->length;
            SS->alloc = S->alloc;
            printf("S: "); fmpz_mpoly_print_pretty(SS, NULL, ctx); printf("\n");
            
            fmpz_mpoly_ctx_clear(ctx);
        }
*/
                }
                fmpz_mpoly_swap(S, T, ctx);

            } while ((x = x->next) != NULL);
        } while (heap_len > 1 && heap[1].exp == exp);

        /* process nodes taken from the heap */
        while (store > store_base)
        {
            j = *--store;
            i = *--store;

            if (i == -WORD(1))
            {
                /* take next dividend term */
                if (j + 1 < Alen)
                {
                    x = chain + 0;
                    x->i = i;
                    x->j = j + 1;
                    x->next = NULL;
                    _mpoly_heap_insert1(heap, Aexp[x->j], x,
                                                 &next_loc, &heap_len, maskhi);
                }
            }
            else
            {
                /* should we go right? */
                if (  (i + 1 < Blen)
                   && (hind[i + 1] == 2*j + 1)
                   )
                {
                    x = chain + i + 1;
                    x->i = i + 1;
                    x->j = j;
                    x->next = NULL;
                    hind[x->i] = 2*(x->j + 1) + 0;
                    _mpoly_heap_insert1(heap, Bexp[x->i] + Q->exps[x->j], x,
                                                 &next_loc, &heap_len, maskhi);
                }
                /* should we go up? */
                if (j + 1 == Q->length)
                {
                    s++;
                } else if (  ((hind[i] & 1) == 1)
                          && ((i == 1) || (hind[i - 1] >= 2*(j + 2) + 1))
                          )
                {
                    x = chain + i;
                    x->i = i;
                    x->j = j + 1;
                    x->next = NULL;
                    hind[x->i] = 2*(x->j + 1) + 0;
                    _mpoly_heap_insert1(heap, Bexp[x->i] + Q->exps[x->j], x,
                                                 &next_loc, &heap_len, maskhi);
                }
            }
        }

        if (T->length == 0)
        {
            continue;
        }

        if (mpoly_monomials_overflow_test(T->exps, T->length, bits, ctx->minfo))
        {
            goto not_exact_division;
        }

        q = Q->coeffs + Q->length;
        FLINT_ASSERT(q->bits == bits);
        b = Bcoeff + 0;

/*
        {
            fmpz_mpoly_ctx_t ctx;
            fmpz_mpoly_t TT, bb;

            fmpz_mpoly_ctx_init(ctx, 2, ORD_LEX);

            TT->bits = bits;
            TT->coeffs = T->coeffs;
            TT->exps = T->exps;
            TT->length = T->length;
            TT->alloc = T->alloc;
            printf("T: "); fmpz_mpoly_print_pretty(TT, NULL, ctx); printf("\n");

            bb->bits = bits;
            bb->coeffs = b->coeffs;
            bb->exps = b->exps;
            bb->length = b->length;
            bb->alloc = b->alloc;
            printf("b: "); fmpz_mpoly_print_pretty(bb, NULL, ctx); printf("\n");
            
            fmpz_mpoly_ctx_clear(ctx);
        }
*/

        q->length = _fmpz_mpoly_divides_monagan_pearce(
                            &q->coeffs, &q->exps, &q->alloc,
                            T->coeffs, T->exps, T->length,
                            b->coeffs, b->exps, b->length,
                                              bits, N, cmpmask);
        if (q->length == 0)
        {
            goto not_exact_division;
        }

        if (!lt_divides || (exp^maskhi) < (maxexp^maskhi))
        {
            goto not_exact_division;
        }
/*
        {
            fmpz_mpoly_ctx_t ctx;
            fmpz_mpoly_t qq;

            fmpz_mpoly_ctx_init(ctx, 2, ORD_LEX);


flint_printf("qexp: %016wx\n", Q->exps[Q->length]);
            qq->bits = bits;
            qq->coeffs = q->coeffs;
            qq->exps = q->exps;
            qq->length = q->length;
            qq->alloc = q->alloc;
            printf("q: "); fmpz_mpoly_print_pretty(qq, NULL, ctx); printf("\n");
            
            fmpz_mpoly_ctx_clear(ctx);
        }
*/
        /* put newly generated quotient term back into the heap if neccesary */
        if (s > 1)
        {
            i = 1;
            x = chain + i;
            x->i = i;
            x->j = Q->length;
            x->next = NULL;
            hind[x->i] = 2*(x->j + 1) + 0;
            _mpoly_heap_insert1(heap, Bexp[x->i] + Q->exps[x->j], x,
                                                 &next_loc, &heap_len, maskhi);
        }
        s = 1;
        Q->length++;
    }

    success = 1;

cleanup:

    fmpz_mpoly_clear(T, ctx);
    fmpz_mpoly_clear(S, ctx);

    TMP_END;

    return success;

not_exact_division:

    success = 0;
    Q->length = 0;
    goto cleanup;
}


/*
    A is in ZZ[x_0,..., x_(n-1)]
    out is in Fp[x_var]   (p is stored in out :( )

    x_0 = values[0]
    x_1 = values[1]
    ...
*/
slong fmpz_mpoly_eval_all_but_one_nmod(
    nmod_poly_t out,
    const fmpz_mpoly_t A,
    slong var,
    mp_limb_t * values,
    const fmpz_mpoly_ctx_t ctx)
{
    slong i, j;
    slong deg;
    ulong varexp, thisexp;
    mp_limb_t t, v;
    ulong mask;
    slong * offsets, * shifts;
    slong N = mpoly_words_per_exp_sp(A->bits, ctx->minfo);
    ulong * Aexp = A->exps;
    fmpz * Acoeff = A->coeffs;
    TMP_INIT;

    FLINT_ASSERT(A->bits <= FLINT_BITS);

    TMP_START;

    mask = (-UWORD(1)) >> (FLINT_BITS - A->bits);
    offsets = (slong *) TMP_ALLOC(ctx->minfo->nvars*sizeof(slong));
    shifts = (slong *) TMP_ALLOC(ctx->minfo->nvars*sizeof(slong));
    for (j = 0; j < ctx->minfo->nvars; j++)
    {
        mpoly_gen_offset_shift_sp(offsets + j, shifts + j, j, A->bits, ctx->minfo);
    }

    nmod_poly_zero(out);
    deg = -WORD(1);
    for (i = 0; i < A->length; i++)
    {
        varexp = ((Aexp + N*i)[offsets[var]]>>shifts[var])&mask;
        deg = FLINT_MAX(deg, (slong)(varexp));
        v = fmpz_fdiv_ui(Acoeff + i, out->mod.n);
        for (j = 0; j < ctx->minfo->nvars; j++)
        {
            thisexp = ((Aexp + N*i)[offsets[j]]>>shifts[j])&mask;
            if (j != var)
            {
                v = nmod_mul(v, nmod_pow_ui(values[j], thisexp, out->mod), out->mod);
            }
        }
        t = nmod_poly_get_coeff_ui(out, varexp);
        nmod_poly_set_coeff_ui(out, varexp, nmod_add(t, v, out->mod));
    }

    TMP_END;
    return deg;
}

/*
    A is in ZZ[X,Y][x_0,..., x_(n-1)]
    out is in Fp[x_var]   (p is stored in out :( )

    X = values[0]
    Y = values[1]
    x_0 = values[2]
    ...
*/
slong fmpz_mpolyuu_eval_all_but_one_nmod(
    nmod_poly_t out,
    const fmpz_mpolyu_t A,
    slong var,
    mp_limb_t * values,
    const fmpz_mpoly_ctx_t ctx)
{
    slong i;
    slong deg, thisdeg;
    ulong * Aexp = A->exps;
    const fmpz_mpoly_struct * Acoeff = A->coeffs;
    nmod_poly_t temp;
    mp_limb_t t, t1;

    FLINT_ASSERT(A->bits <= FLINT_BITS);

    nmod_poly_zero(out);
    nmod_poly_init(temp, out->mod.n);

    deg = -WORD(1);
    for (i = 0; i < A->length; i++)
    {
        t = nmod_pow_ui(values[0], Aexp[i] >> (FLINT_BITS/2), out->mod);
        t1 = nmod_pow_ui(values[1], Aexp[i] & ((-UWORD(1)) >> (FLINT_BITS/2)), out->mod);
        t = nmod_mul(t, t1, out->mod);
        thisdeg = fmpz_mpoly_eval_all_but_one_nmod(temp, Acoeff + i, var, values + 2, ctx);
        deg = FLINT_MAX(deg, thisdeg);
        nmod_poly_scalar_mul_nmod(temp, temp, t);
        nmod_poly_add(out, out, temp);
    }

    nmod_poly_clear(temp);   
    return deg; 
}

/*
    A, B are in ZZ[X,Y][x_0, ..., x_(n-1)]

    Set Adeg (resp Bdeg) to the degree of A (resp B) wrt x_var.
    Return an upper bound on the degree of gcd(A,B) wrt x_var.
*/
slong fmpz_mpolyuu_gcd_degree_bound_minor(
    slong * Adeg,
    slong * Bdeg,
    const fmpz_mpolyu_t A,
    const fmpz_mpolyu_t B,
    slong var,
    const fmpz_mpoly_ctx_t ctx,
    flint_rand_t state)
{
    slong i;
    int tries = 0;
    slong degA, degB, degRet;
    mp_limb_t p = UWORD(1) << (FLINT_BITS - 2);
    nmod_poly_t Geval, Aeval, Beval;
    mp_limb_t * values;
    TMP_INIT;

    TMP_START;

    values = (mp_limb_t *) TMP_ALLOC(ctx->minfo->nvars*sizeof(mp_limb_t));

    p = n_nextprime(p, 1);
    nmod_poly_init(Geval, p);
    nmod_poly_init(Aeval, p);
    nmod_poly_init(Beval, p);

try_again:

    for (i = 0; i < ctx->minfo->nvars + 2; i++)
    {
        values[i] = n_urandint(state, p);
    }

    degA = fmpz_mpolyuu_eval_all_but_one_nmod(Aeval, A, var, values, ctx);
    degB = fmpz_mpolyuu_eval_all_but_one_nmod(Beval, B, var, values, ctx);
    *Adeg = degA;
    *Bdeg = degB;

    if (degA != nmod_poly_degree(Aeval) || degB != nmod_poly_degree(Beval))
    {
        if (++tries > 100)
        {
            degRet = FLINT_MIN(degA, degB);
            goto cleanup;
        }
        p = n_nextprime(p, 1);
        nmod_poly_clear(Geval);
        nmod_poly_clear(Aeval);
        nmod_poly_clear(Beval);
        nmod_poly_init(Geval, p);
        nmod_poly_init(Aeval, p);
        nmod_poly_init(Beval, p);
        goto try_again;
    }

    nmod_poly_gcd(Geval, Aeval, Beval);
    degRet = nmod_poly_degree(Geval);

cleanup:

    nmod_poly_clear(Geval);
    nmod_poly_clear(Aeval);
    nmod_poly_clear(Beval);
    TMP_END;
    return degRet;
}

/*
    A is in ZZ[x_0, ..., x_(n-1)]

    After the substitutions
        x_0     = x ^ (sub[1] * sub[2] * ... * sub[n-1])

        x_(n-2) = x ^ (sub[n-1])
        x_(n-1) = x ^ (1)
    a univariate in ZZ[x] remainds. Return the content of this poly.
*/
void fmpz_mpoly_ksub_content(
    fmpz_t content,
    const fmpz_mpoly_t A,
    const ulong * subdegs,
    const fmpz_mpoly_ctx_t ctx)
{
    slong i, j;
    slong N = mpoly_words_per_exp_sp(A->bits, ctx->minfo);
    fmpz_mpoly_t T;
    fmpz_mpoly_ctx_t Tctx;
    fmpz_t e;
    fmpz * Acoeff = A->coeffs;
    ulong * Aexp = A->exps;
    ulong mask;
    slong * offsets, * shifts;
    TMP_INIT;

    TMP_START;
    fmpz_init(e);

    fmpz_mpoly_ctx_init(Tctx, 1, ORD_LEX);
    fmpz_mpoly_init(T, Tctx);

    mask = (-UWORD(1)) >> (FLINT_BITS - A->bits);
    offsets = (slong *) TMP_ALLOC(ctx->minfo->nvars*sizeof(slong));
    shifts = (slong *) TMP_ALLOC(ctx->minfo->nvars*sizeof(slong));
    for (j = 0; j < ctx->minfo->nvars; j++)
    {
        mpoly_gen_offset_shift_sp(offsets + j, shifts + j, j, A->bits, ctx->minfo);
    }

    for (i = 0; i < A->length; i++)
    {
        fmpz_zero(e);
        for (j = 0; j < ctx->minfo->nvars; j++)
        {
            fmpz_mul_ui(e, e, subdegs[j]);
            fmpz_add_ui(e, e, ((Aexp + N*i)[offsets[j]]>>shifts[j])&mask);
        }
        _fmpz_mpoly_push_exp_ffmpz(T, e, Tctx);
        fmpz_set(T->coeffs + T->length - 1, Acoeff + i);
    }
    fmpz_mpoly_sort_terms(T, Tctx);
    fmpz_mpoly_combine_like_terms(T, Tctx);
    _fmpz_vec_content(content, T->coeffs, T->length);
    fmpz_mpoly_clear(T, Tctx);
    fmpz_mpoly_ctx_clear(Tctx);

    fmpz_clear(e);
    TMP_END;
}

typedef enum {
    random_check_good,
    random_check_point_not_found,
    random_check_image_degree_high,
    random_check_image_degree_low,
    random_check_image_no_match    
} random_check_ret_t;



random_check_ret_t static _random_check_sp(
    ulong * GevaldegXY,
    ulong GdegboundXY,
    nmod_mpolyun_t Aeval_sp,
    nmod_mpolyun_t Beval_sp,
    nmod_mpolyun_t Geval_sp,
    nmod_mpolyun_t Abareval_sp,
    nmod_mpolyun_t Bbareval_sp,
    mp_limb_t * checkalpha_sp,
    const fmpz_mpolyu_t H,
    const fmpz_mpolyu_t A,
    const fmpz_mpolyu_t B,
    const fmpz_mpoly_t Gamma,
    const fmpz_mpoly_ctx_t ctx,
    const nmod_mpoly_ctx_t ctx_sp,
    flint_rand_t randstate)
{
    mp_limb_t Gammaeval_sp;
    int success;
    int point_try_count;
    slong i;

    /* try to test H at a random evaluation point */
    for (point_try_count = 0; point_try_count < 10; point_try_count++)
    {
        /* evaluate Gamma, A, and B at random point */
        for (i = 0; i < ctx->minfo->nvars; i++)
        {
            checkalpha_sp[i] = n_urandint(randstate, ctx_sp->ffinfo->mod.n);
/*
flint_printf("alpha_sp[%wu]: %wu\n", i, checkalpha_sp[i]);
*/
        }
        fmpz_mpolyuu_eval_nmod(Aeval_sp, ctx_sp, A, checkalpha_sp, ctx);
        fmpz_mpolyuu_eval_nmod(Beval_sp, ctx_sp, B, checkalpha_sp, ctx);

/*
printf("Aeval_sp: "); nmod_mpolyun_print_pretty(Aeval_sp, NULL, ctx_sp); printf("\n");
printf("Beval_sp: "); nmod_mpolyun_print_pretty(Beval_sp, NULL, ctx_sp); printf("\n");
*/
        /* make sure that evaluation did not kill either lc(A) or lc(B) */
        if ( Aeval_sp->length == 0 || Beval_sp->length == 0 
            || nmod_mpolyun_bidegree(Aeval_sp, ctx_sp) != A->exps[0]
            || nmod_mpolyun_bidegree(Beval_sp, ctx_sp) != B->exps[0])
        {
            continue;
        }

        /* Gamma is gcd(lc(A), lc(B)) so it evaluation should not be zero */
        Gammaeval_sp = fmpz_mpoly_eval_nmod(ctx_sp->ffinfo, Gamma, checkalpha_sp, ctx);
        FLINT_ASSERT(Gammaeval_sp != 0);
/*
flint_printf("Gammaeval_sp: %wu\n", Gammaeval_sp);
*/
        success = nmod_mpolyun_gcd_brown_smprime_bivar(Geval_sp, Abareval_sp, Bbareval_sp, Aeval_sp, Beval_sp, ctx_sp);
        if (!success)
        {
            continue;
        }
        nmod_mpolyun_scalar_mul_nmod(Geval_sp, Gammaeval_sp, ctx_sp);
/*
printf("Geval_sp: "); nmod_mpolyun_print_pretty(Geval_sp, NULL, ctx_sp); printf("\n");
*/
        FLINT_ASSERT(Geval_sp->length > 0);
        *GevaldegXY = nmod_mpolyun_bidegree(Geval_sp, ctx_sp);

        if (GdegboundXY < *GevaldegXY)
        {
            return random_check_image_degree_high;
        }
        else if (GdegboundXY > *GevaldegXY)
        {
            return random_check_image_degree_low;
        }

        /* reuse Bbareval for Heval */
        fmpz_mpolyuu_eval_nmod(Bbareval_sp, ctx_sp, H, checkalpha_sp, ctx);
/*
printf("Heval_sp: "); nmod_mpolyun_print_pretty(Bbareval_sp, NULL, ctx_sp); printf("\n");
*/

        if (!nmod_mpolyun_equal(Bbareval_sp, Geval_sp, ctx_sp))
        {
            return random_check_image_no_match;
        }

        return random_check_good;
    }

    return random_check_point_not_found;
}





random_check_ret_t static _random_check(
    ulong * GevaldegXY,
    ulong GdegboundXY,
    fmpz_mod_mpolyun_t Aeval,
    fmpz_mod_mpolyun_t Beval,
    fmpz_mod_mpolyun_t Geval,
    fmpz_mod_mpolyun_t Abareval,
    fmpz_mod_mpolyun_t Bbareval,
    fmpz_t Gammaeval,
    fmpz * checkalpha,
    const fmpz_mpolyu_t H,
    const fmpz_mpolyu_t A,
    const fmpz_mpolyu_t B,
    const fmpz_mpoly_t Gamma,
/*    const fmpz_mod_mpoly_bma_interpolate_ctx_t Ictx,*/
    const fmpz_mpoly_ctx_t ctx,
    const fmpz_mod_ctx_t fpctx,
    flint_rand_t randstate)
{
    int success;
    int point_try_count;
    slong i;

    /* try to test H at a random evaluation point */
    for (point_try_count = 0; point_try_count < 10; point_try_count++)
    {
        /* evaluate Gamma, A, and B at random point */
        for (i = 0; i < ctx->minfo->nvars; i++)
        {
            fmpz_randm(checkalpha + i, randstate, fpctx->p);
        }
        fmpz_mpolyuu_eval_fmpz_mod(Aeval, A, checkalpha, ctx, fpctx);
        fmpz_mpolyuu_eval_fmpz_mod(Beval, B, checkalpha, ctx, fpctx);

        /* make sure that evaluation did not kill either lc(A) or lc(B) */
        if ( Aeval->length == 0 || Beval->length == 0 
            || fmpz_mod_mpolyun_bidegree(Aeval, ctx, fpctx) != A->exps[0]
            || fmpz_mod_mpolyun_bidegree(Beval, ctx, fpctx) != B->exps[0])
        {
            continue;
        }

        /* Gamma is gcd(lc(A), lc(B)) so it evaluation should not be zero */
        fmpz_mpoly_eval_fmpz_mod(Gammaeval, Gamma, checkalpha, ctx, fpctx);
        FLINT_ASSERT(!fmpz_is_zero(Gammaeval));

        success = fmpz_mod_mpolyun_gcd_bivar(Geval, Abareval, Bbareval, Aeval, Beval, ctx, fpctx);
        if (!success)
        {
            continue;
        }
        fmpz_mod_mpolyun_scalar_mul_fmpz_mod(Geval, Gammaeval, ctx, fpctx);
/*    printf("Geval: "); fmpz_mod_mpolyun_print_pretty(Geval, NULL, Ictx->polyctx, Ictx->fpctx); printf("\n");*/

        FLINT_ASSERT(Geval->length > 0);
        *GevaldegXY = fmpz_mod_mpolyun_bidegree(Geval, ctx, fpctx);

        if (GdegboundXY < *GevaldegXY)
        {
            return random_check_image_degree_high;
        }
        else if (GdegboundXY > *GevaldegXY)
        {
            return random_check_image_degree_low;
        }

        /* reuse Bbareval for Heval */
        fmpz_mpolyuu_eval_fmpz_mod(Bbareval, H, checkalpha, ctx, fpctx);
        if (!fmpz_mod_mpolyun_equal(Bbareval, Geval, ctx, fpctx))
        {
            return random_check_image_no_match;
        }

        return random_check_good;
    }

    return random_check_point_not_found;
}






/*
Assumptions:

    A, B, G are in ZZ[X,Y][x_0, ..., x_(n-1)]    n = ctx->minfo->nvars
    Gamma is in ZZ[x_0, ..., x_(n-1)]
    The ordering of the variables is X > Y > x_0 > ... > x_(n-1).

    Since A and B are packed into fmpz_mpolyu_t, the main degrees auto satisfy
        deg_X(A) < 2^(FLINT_BITS/2 - 1)
        deg_Y(A) < 2^(FLINT_BITS/2 - 1)
    ditto for B.

    A, B are contentless wrt X, Y: the content in ZZ[x_0, ..., x_(n-1)] is 1.
    Gamma is the gcd of lc_XY(A) and lc_XY(B). Thus Gamma is a multiple of lc(G).

    The degree of things wrt the main pair of variables XY is always maintained
    as a bidegree, that is, as X^i*Y^j with comparisons using lex.

Procedure:

    We will construct H in ZZ[X,Y][x_0, ..., x_(n-1)] with lc_XY(H) = Gamma
    that is a ZZ[x_0, ..., x_(n-1)]-multiple of the true GCD. Then
    G is Hpp = H/content_XY(H) where content_XY(H) is in ZZ[x_0, ..., x_(n-1)]

    The kronecker substitution (ksub) in the lesser variables always takes the form
        x_0     = x ^ (subdegs[n-1] * subdegs[n-2] * ... * subdegs[1])
          ...
        x_(n-3) = x ^ (subdegs[n-1] * subdegs[n-2])
        x_(n-2) = x ^ (subdegs[n-1])
        x_(n-1) = x ^ (1)
    After such a ksub everything is in ZZ[X,Y][x].

    H is found from its image modulo a smooth prime p. This works by applying
    the ksub and interpolating coeffs in Fp[x] via bivariate images in Fp[X,Y]
    at values x = alpha^1, x = alpha^2, ... (where alpha is a gen of Fp*)
    then reversing the ksub to Fp[x_0, ..., x_(n-1)]. Then, we assume that
        (1) this first image of H is correct mod p, and
        (2) all terms present in H are present in this first image (zippel asmp)
    Zippel interpolation is used to find more images modulo machine primes p'.
    If (2) is false, i.e. p divides a coefficient of H, then Zippel interpolation
    will probably find a contradiction, in which case we need to start all over.

    H/content_XY(H) is the true gcd if it passes the divisibility test.
    This is tested after the candidate H stabilizes under chinese remaindering.
*/
int fmpz_mpolyuu_gcd_bma(
    fmpz_mpolyu_t G,
    const fmpz_mpolyu_t A,
    const fmpz_mpolyu_t B,
    const fmpz_mpoly_t Gamma,
    const fmpz_mpoly_ctx_t ctx)
{
    int changed, success, point_try_count;
    mp_bitcnt_t bits = A->bits, Hbits;
    mpoly_bma_interpolate_ctx_t Ictx;
    fmpz_mod_ctx_t fpctx;
    nmod_mpoly_ctx_t ctx_sp;
    fmpz_mpolyu_t H, Hpp, Abar, Bbar;
    fmpz_mpoly_t Hcontent;

    fmpz_mod_bma_mpoly_t Lambda;
    fmpz_mod_mpolyun_t Aeval, Beval, Geval, Abareval, Bbareval;
    fmpz_mpolycu_t Ainc, Acur, Binc, Bcur, Ared, Bred;
    fmpz_mpolyc_t Gammainc, Gammacur, Gammared;

    nmod_bma_mpoly_t Lambda_sp;
    nmod_mpolyun_t Aeval_sp, Beval_sp, Geval_sp, Abareval_sp, Bbareval_sp;
    nmod_mpolycu_t Ainc_sp, Acur_sp, Binc_sp, Bcur_sp, Ared_sp, Bred_sp;
    nmod_mpolyc_t Gammainc_sp, Gammacur_sp, Gammared_sp;

    mp_limb_t p_sp, sshift_sp, last_unlucky_sshift_plus_1_sp, image_count_sp;

    slong i, j;
    ulong GdegboundXY, GevaldegXY;
    slong * Gdegbounds, * Adegs, * Bdegs, * Gammadegs;
    flint_rand_t randstate;
    fmpz_t p, t, shift, subprod, cAksub, cBksub, sshift, last_unlucky_sshift_plus_1, image_count;
    fmpz_t Gammaeval;
    mp_limb_t Gammaeval_sp;
    fmpz * checkalpha;
    mp_limb_t * checkalpha_sp;
    int unlucky_count;
    fmpz_t Hmodulus;
    nmod_zip_mpolyu_t Z;
    slong zip_evals;

    FLINT_ASSERT(bits == A->bits);
    FLINT_ASSERT(bits == B->bits);
    FLINT_ASSERT(bits == G->bits);
    FLINT_ASSERT(bits == Gamma->bits);
/*
printf("\n****************\nfmpz_mpolyuu_gcd_bma called\n");
printf("Gamma: "); fmpz_mpoly_print_pretty(Gamma, NULL, ctx); printf("\n");
printf("    A: "); fmpz_mpolyuu_print_pretty(A, NULL, 2, ctx); printf("\n");
printf("    B: "); fmpz_mpolyuu_print_pretty(B, NULL, 2, ctx); printf("\n");
*/
    /* let's initialize everything at once to avoid complicated cleanup */

    flint_randinit(randstate);
    fmpz_init(p);
    fmpz_init(t);
    fmpz_init(shift);
    fmpz_init(image_count);
    fmpz_init(subprod);
    fmpz_init(cAksub);
    fmpz_init(cBksub);
    fmpz_init(sshift);
    fmpz_init(last_unlucky_sshift_plus_1);
    fmpz_init(Gammaeval);
    fmpz_init(Hmodulus);

    mpoly_bma_interpolate_ctx_init(Ictx, ctx->minfo->nvars);


    /* multiprecision workspace */
    fmpz_set_ui(p, 997);
    fmpz_mod_ctx_init(fpctx, p);
    fmpz_mod_bma_mpoly_init(Lambda);

    fmpz_mod_mpolyun_init(Aeval, bits, ctx, fpctx);
    fmpz_mod_mpolyun_init(Beval, bits, ctx, fpctx);
    fmpz_mod_mpolyun_init(Geval, bits, ctx, fpctx);
    fmpz_mod_mpolyun_init(Abareval, bits, ctx, fpctx);
    fmpz_mod_mpolyun_init(Bbareval, bits, ctx, fpctx);
    fmpz_mpolyc_init(Gammainc);
    fmpz_mpolyc_init(Gammacur);
    fmpz_mpolyc_init(Gammared);
    fmpz_mpolycu_init(Ainc);
    fmpz_mpolycu_init(Acur);
    fmpz_mpolycu_init(Ared);
    fmpz_mpolycu_init(Binc);
    fmpz_mpolycu_init(Bcur);
    fmpz_mpolycu_init(Bred);
    checkalpha = (fmpz *) flint_malloc(ctx->minfo->nvars*sizeof(fmpz));
    for (i = 0; i < ctx->minfo->nvars; i++)
    {
        fmpz_init(checkalpha + i);
    }

    /* machine precision workspace */
    nmod_mpoly_ctx_init_mpoly_ctx(ctx_sp, ctx->minfo, 2); /* mod no care */
    nmod_bma_mpoly_init(Lambda_sp);

    nmod_mpolyun_init(Aeval_sp, bits, ctx_sp);
    nmod_mpolyun_init(Beval_sp, bits, ctx_sp);
    nmod_mpolyun_init(Geval_sp, bits, ctx_sp);
    nmod_mpolyun_init(Abareval_sp, bits, ctx_sp);
    nmod_mpolyun_init(Bbareval_sp, bits, ctx_sp);
    nmod_mpolyc_init(Gammainc_sp);
    nmod_mpolyc_init(Gammacur_sp);
    nmod_mpolyc_init(Gammared_sp);
    nmod_mpolycu_init(Ainc_sp);
    nmod_mpolycu_init(Acur_sp);
    nmod_mpolycu_init(Ared_sp);
    nmod_mpolycu_init(Binc_sp);
    nmod_mpolycu_init(Bcur_sp);
    nmod_mpolycu_init(Bred_sp);
    checkalpha_sp = (mp_limb_t *) flint_malloc(ctx->minfo->nvars*sizeof(mp_limb_t));

    /* the zippler */
    nmod_zip_mpolyu_init(Z);

    Gdegbounds = (slong *) flint_malloc(ctx->minfo->nvars*sizeof(slong));
    Adegs      = (slong *) flint_malloc(ctx->minfo->nvars*sizeof(slong));
    Bdegs      = (slong *) flint_malloc(ctx->minfo->nvars*sizeof(slong));
    Gammadegs  = (slong *) flint_malloc(ctx->minfo->nvars*sizeof(slong));

    /* find a degree bound on G in the two main variables */
    GdegboundXY = FLINT_MIN(A->exps[0], B->exps[0]);
    p_sp = UWORD(1) << (FLINT_BITS - 2);
    for (point_try_count = 0; point_try_count < 10; point_try_count++)
    {
        p_sp = n_nextprime(p_sp, 1);
        nmod_mpoly_ctx_set_mod(ctx_sp, p_sp);
        /* unfortunate nmod_poly's need mod set */
        nmod_mpolyun_set_mod(Aeval_sp, ctx_sp->ffinfo->mod);
        nmod_mpolyun_set_mod(Beval_sp, ctx_sp->ffinfo->mod);
        nmod_mpolyun_set_mod(Geval_sp, ctx_sp->ffinfo->mod);
        nmod_mpolyun_set_mod(Abareval_sp, ctx_sp->ffinfo->mod);
        nmod_mpolyun_set_mod(Bbareval_sp, ctx_sp->ffinfo->mod);
        for (i = 0; i < ctx->minfo->nvars; i++)
        {
            checkalpha_sp[i] = n_urandint(randstate, p_sp);
        }
        fmpz_mpolyuu_eval_nmod(Aeval_sp, ctx_sp, A, checkalpha_sp, ctx);
        fmpz_mpolyuu_eval_nmod(Beval_sp, ctx_sp, B, checkalpha_sp, ctx);

        if (Aeval_sp->length == 0 || Beval_sp->length == 0
            || nmod_mpolyun_bidegree(Aeval_sp, ctx_sp) != A->exps[0]
            || nmod_mpolyun_bidegree(Beval_sp, ctx_sp) != B->exps[0])
        {
            /* evaluation killed at least one of lc(A) or lc(B) */
            continue;
        }
        success = nmod_mpolyun_gcd_brown_smprime_bivar(Geval_sp, Abareval_sp, Bbareval_sp, Aeval_sp, Beval_sp, ctx_sp);
        if (success)
        {
            FLINT_ASSERT(Geval_sp->length > 0);
            GdegboundXY = nmod_mpolyun_bidegree(Geval_sp, ctx_sp);
            break;
        }
    }
/*flint_printf("GdegboundXY: %016llx\n", GdegboundXY);*/

    /*
        Find degree bounds on G wrt lesser variables so that
            Gdegbounds[i] >= deg_(x_i)(G)
        Also fills in
            Adegs[i] = deg_(x_i)(A)
            Bdegs[i] = deg_(x_i)(B)
    */
    mpoly_degrees_si(Gammadegs, Gamma->exps, Gamma->length, bits, ctx->minfo);
    for (i = 0; i < ctx->minfo->nvars; i++)
    {
        Gdegbounds[i] = fmpz_mpolyuu_gcd_degree_bound_minor(Adegs + i, Bdegs + i, A, B, i, ctx, randstate);
/*flint_printf("Gdegbound[%wd]: %wd\n", i, Gdegbounds[i]);*/
    }

    /*
        Find bits into which H can be packed. The degrees satsify
            deg_(x_i)(H) <= deg_(x_i)(A)
            deg_(x_i)(H) <= deg_(x_i)(B)
            deg_(x_i)(H) <= deg_(x_i)(Gamma) + deg_(x_i)(G)
    */
    Hbits = bits;
    for (i = 0; i < ctx->minfo->nvars; i++)
    {
        mp_bitcnt_t Hibits;
        Ictx->degbounds[i] = FLINT_MIN(Adegs[i], Bdegs[i]);
        Ictx->degbounds[i] = FLINT_MIN(Ictx->degbounds[i], Gdegbounds[i] + Gammadegs[i]);
        Hibits = 1 + FLINT_BIT_COUNT(Ictx->degbounds[i]);
        Hbits = FLINT_MAX(Hbits, Hibits);

        /* degbounds[i] will be a strict degree bound on deg_(x_i)(H) */
        Ictx->degbounds[i]++;
    }

    fmpz_mpolyu_init(Abar, bits, ctx);
    fmpz_mpolyu_init(Bbar, bits, ctx);
    fmpz_mpolyu_init(H, Hbits, ctx);
    fmpz_mpolyu_init(Hpp, Hbits, ctx);
    fmpz_mpoly_init3(Hcontent, 0, Hbits, ctx);

    /* initialization done */

    if (GdegboundXY == 0)
    {
        fmpz_mpolyu_one(G, ctx);
        success = 1;
        goto cleanup;
    }

    if (Hbits > FLINT_BITS)
    {
        /* H cannot be guaranteed to be packed into FLINT_BITS - absolute falure */
        success = 0;
        goto cleanup;
    }

    /* initial choices for the ksub degrees are the strict degree bounds on H */
    for (i = 0; i < ctx->minfo->nvars; i++)
    {
        Ictx->subdegs[i] = Ictx->degbounds[i];
        if ((slong)(Ictx->subdegs[i]) <= 0)
        {
            /* ksub has overflown - absolute falure */
            success = 0;
            goto cleanup;
        }
    }
    goto got_ksub;

pick_ksub:

    if (ctx->minfo->nvars > 1)
    {
        /* just increment the smallest subdegs[j] */
        j = 1;
        for (i = 2; i < ctx->minfo->nvars; i++)
        {
            if (Ictx->subdegs[i] < Ictx->subdegs[j])
            {
                j = i;
            }
        }
        Ictx->subdegs[j]++;
        if ((slong)(Ictx->subdegs[j]) < 0)
        {
            /* ksub has overflown - absolute falure */
            success = 0;
            goto cleanup;
        }
    }

got_ksub:

    fmpz_one(subprod);
    for (i = 0; i < ctx->minfo->nvars; i++)
    {
        fmpz_mul_ui(subprod, subprod, Ictx->subdegs[i]);
/*
flint_printf("ksub[%wu]: %wu  ", i, Ictx->subdegs[i]);
*/
    }
/*
flint_printf("\n");
*/

    /* see if the ksub killed either lc(A) or lc(B) */
    fmpz_mpoly_ksub_content(cAksub, A->coeffs + 0, Ictx->subdegs, ctx);
    fmpz_mpoly_ksub_content(cBksub, B->coeffs + 0, Ictx->subdegs, ctx);
    if (fmpz_is_zero(cAksub) || fmpz_is_zero(cBksub))
    {
        /* try a new substitution if we killed either leading coefficient */
        goto pick_ksub;
    }

pick_bma_prime:
    /*
        Pick a prime p for first image. If p is large it should be smooth so
        that logs in Fp are possible. It should also be big enough so that the
        ksub is reversible.
    */
    if (fmpz_cmp(p, subprod) < 0)
        fmpz_set(p, subprod);
    fmpz_add_ui(p, p, 1);
    if (fmpz_is_probabprime(p) != 1)
        goto pick_bma_prime;    
    if (fmpz_is_prime(p) != 1)
        goto pick_bma_prime;

    /* make sure reduction does not kill either leading coeff after ksub */
    fmpz_gcd(t, p, cAksub);
    if (!fmpz_is_one(t))
        goto pick_bma_prime;
    fmpz_gcd(t, p, cBksub);
    if (!fmpz_is_one(t))
        goto pick_bma_prime;

    /* make sure p does not divide any coefficient of Gamma */
    for (i = 0; i < Gamma->length; i++)
    {
        if (fmpz_divisible(Gamma->coeffs + i, p))
            goto pick_bma_prime;        
    }

    if (fmpz_abs_fits_ui(p))
    {
        p_sp = fmpz_get_ui(p);
        sshift_sp = 1;

        unlucky_count = 0;
        last_unlucky_sshift_plus_1_sp = 0;
/*
    flint_printf("----- bma prime_sm: %wu\n", p_sp);
*/
        nmod_mpoly_ctx_set_mod(ctx_sp, p_sp);
        mpoly_bma_interpolate_ctx_reset_prime_sp(Ictx, p_sp);
        nmod_bma_mpoly_reset_prime(Lambda_sp, ctx_sp->ffinfo);
        nmod_bma_mpoly_zero(Lambda_sp);

        /* unfortunate nmod_poly's store their own ctx :( */
        nmod_mpolyun_set_mod(Aeval_sp, ctx_sp->ffinfo->mod);
        nmod_mpolyun_set_mod(Beval_sp, ctx_sp->ffinfo->mod);
        nmod_mpolyun_set_mod(Geval_sp, ctx_sp->ffinfo->mod);
        nmod_mpolyun_set_mod(Abareval_sp, ctx_sp->ffinfo->mod);
        nmod_mpolyun_set_mod(Bbareval_sp, ctx_sp->ffinfo->mod);
/*
    printf("Lambda: "); nmod_bma_mpoly_print(Lambda_sp); printf("\n");
*/
        FLINT_ASSERT(sshift_sp == 1);
        nmod_mpoly_bma_interpolate_alpha_powers(checkalpha_sp, sshift_sp, Ictx, ctx_sp);
/*
for (i = 0; i < ctx->minfo->nvars; i++)
{
    flint_printf("alpha_sp[%wd]: %wu\n",i,checkalpha_sp[i]);
}
*/

        /* set evaluation of monomials */
        nmod_mpoly_set_skel(Gammainc_sp, ctx_sp, Gamma, checkalpha_sp, ctx);
        nmod_mpolyu_set_skel(Ainc_sp, ctx_sp, A, checkalpha_sp, ctx);
        nmod_mpolyu_set_skel(Binc_sp, ctx_sp, B, checkalpha_sp, ctx);

        /* set reduction of coeffs */
        nmod_mpoly_red_skel(Gammared_sp, Gamma, ctx_sp->ffinfo);
        nmod_mpolyu_red_skel(Ared_sp, A, ctx_sp->ffinfo);
        nmod_mpolyu_red_skel(Bred_sp, B, ctx_sp->ffinfo);

        /* copy evaluation of monomials */
        nmod_mpoly_copy_skel(Gammacur_sp, Gammainc_sp);
        nmod_mpolyu_copy_skel(Acur_sp, Ainc_sp);
        nmod_mpolyu_copy_skel(Bcur_sp, Binc_sp);

        image_count_sp = 0;

    next_bma_image_sp:
/*
    flint_printf("next_bma_image_sp %wu sshift_sp %wu\n", image_count_sp, sshift_sp);
    printf("Lambda: "); nmod_bma_mpoly_print(Lambda_sp); printf("\n");
    usleep(100000);
*/
        /* image count is also the current power of alpha we are evaluating */
        image_count_sp++;
        FLINT_ASSERT(sshift_sp + Lambda_sp->pointcount == image_count_sp);

        if (image_count_sp >= p_sp - 1)
        {
            /* out of evaluation points alpha^image_count in Fp* */
            goto pick_bma_prime;
        }

        Gammaeval_sp = nmod_mpoly_use_skel_mul(Gammared_sp, Gammacur_sp, Gammainc_sp, ctx_sp);
        nmod_mpolyuu_use_skel_mul(Aeval_sp, A, Ared_sp, Acur_sp, Ainc_sp, ctx_sp);
        nmod_mpolyuu_use_skel_mul(Beval_sp, B, Bred_sp, Bcur_sp, Binc_sp, ctx_sp);
/*
flint_printf("Gammaeval_sp: %wu\n", Gammaeval_sp);
printf("Aeval_sp: "); nmod_mpolyun_print_pretty(Aeval_sp, NULL, ctx_sp); printf("\n");
printf("Beval_sp: "); nmod_mpolyun_print_pretty(Beval_sp, NULL, ctx_sp); printf("\n");
*/
        if (Aeval_sp->length == 0 || Beval_sp->length == 0
            || nmod_mpolyun_bidegree(Aeval_sp, ctx_sp) != A->exps[0]
            || nmod_mpolyun_bidegree(Beval_sp, ctx_sp) != B->exps[0])
        {
            /* evaluation killed either lc(A) or lc(B) */
            sshift_sp += Lambda_sp->pointcount + 1;
            nmod_bma_mpoly_zero(Lambda_sp);
            goto next_bma_image_sp;
        }

        /* the evaluation killed neither lc(A) nor lc(B) */
        FLINT_ASSERT(Gammaeval != 0);

        success = nmod_mpolyun_gcd_brown_smprime_bivar(Geval_sp, Abareval_sp, Bbareval_sp, Aeval_sp, Beval_sp, ctx_sp);
        if (!success)
        {
            sshift_sp += Lambda->pointcount + 1;
            nmod_bma_mpoly_zero(Lambda_sp);
            goto next_bma_image_sp;
        }

        FLINT_ASSERT(Geval_sp->length > 0);
        GevaldegXY = nmod_mpolyun_bidegree(Geval_sp, ctx_sp);
        nmod_mpolyun_scalar_mul_nmod(Geval_sp, Gammaeval_sp, ctx_sp);
/*
printf("Geval_sp: "); nmod_mpolyun_print_pretty(Geval_sp, NULL, ctx_sp); printf("\n");
*/
        FLINT_ASSERT(Gammaeval_sp == nmod_mpolyun_leadcoeff_last(Geval_sp, ctx_sp));

        if (GdegboundXY < GevaldegXY)
        {
    printf("bma image unlucky\n");

            /* this image in Fp[X,Y] was unlucky */
            if (sshift_sp == last_unlucky_sshift_plus_1_sp)
            {
                /* this ksub is probably unlucky */
                goto pick_ksub;
            }
            if (++unlucky_count > 2)
            {
                goto pick_bma_prime;
            }
            last_unlucky_sshift_plus_1_sp = sshift_sp + 1;
            sshift_sp += Lambda->pointcount + 1;
            nmod_bma_mpoly_zero(Lambda_sp);
            goto next_bma_image_sp;        
        }
        else if (GdegboundXY > GevaldegXY)
        {
    printf("bma image revealing\n");

            /* new bound on deg_XY(G) */
            sshift_sp += Lambda->pointcount;
            nmod_bma_mpoly_zero(Lambda_sp);
            nmod_bma_mpoly_add_point(Lambda_sp, Geval_sp, ctx_sp->ffinfo);
            GdegboundXY = GevaldegXY;
            if (GdegboundXY == 0)
            {
                fmpz_mpolyu_one(G, ctx);
                success = 1;
                goto cleanup;
            }
            goto next_bma_image_sp;
        }

        nmod_bma_mpoly_add_point(Lambda_sp, Geval_sp, ctx_sp->ffinfo);
/*
printf("Lambda_sp: "); nmod_bma_mpoly_print(Lambda_sp); printf("\n");
*/

        if ((Lambda_sp->pointcount & 1) != 0 || Gamma->length > Lambda_sp->pointcount/2)
        {
            goto next_bma_image_sp;
        }

        changed = nmod_bma_mpoly_reduce(Lambda_sp);
/*
printf("changed: %d\n", changed);
    printf("Lambda_sp: "); nmod_bma_mpoly_print(Lambda_sp); printf("\n");
*/
        if (changed)
        {
            goto next_bma_image_sp;
        }

        success = nmod_bma_mpoly_get_mpolyu(H, ctx, sshift_sp, Lambda_sp, Ictx, ctx_sp->ffinfo);
        if (!success)
        {
printf("could not get mpolyu!!!!!\n");
usleep(1000000);
            goto next_bma_image_sp;
        }

        if (H->length == 0 || (H->coeffs + 0)->length != Gamma->length)
        {
            goto next_bma_image_sp;
        }

        /* GdegboundXY should be the bidegree of H */
        FLINT_ASSERT(GdegboundXY == H->exps[0]);
/*
printf("H: "); fmpz_mpolyuu_print_pretty(H, NULL, 2, ctx); printf("\n");
*/
        switch (_random_check_sp(&GevaldegXY, GdegboundXY,
                    Aeval_sp, Beval_sp, Geval_sp, Abareval_sp, Bbareval_sp,
                            checkalpha_sp, H, A, B, Gamma, ctx, ctx_sp, randstate))
        {
            default:
                FLINT_ASSERT(0);
            case random_check_image_no_match:
            case random_check_image_degree_high:
                goto next_bma_image_sp;
            case random_check_image_degree_low:
                /* the random evaluation point gave us a better degree bound */
                sshift_sp += Lambda_sp->pointcount;
                nmod_bma_mpoly_zero(Lambda_sp);
                GdegboundXY = GevaldegXY;
                if (GdegboundXY == 0)
                {
                    fmpz_mpolyu_one(G, ctx);
                    success = 1;
                    goto cleanup;
                }
                goto next_bma_image_sp;
            case random_check_point_not_found:
                /* hmmm */
            case random_check_good:
                NULL;
        }
    }
    else
    {
        fmpz_one(sshift);

        unlucky_count = 0;
        fmpz_zero(last_unlucky_sshift_plus_1);

    printf("----- bma prime_lg: "); fmpz_print(p); printf("\n");

        fmpz_mod_ctx_set_mod(fpctx, p);
        mpoly_bma_interpolate_ctx_reset_prime(Ictx, p);
        fmpz_mod_bma_mpoly_reset_prime(Lambda, fpctx);
        fmpz_mod_bma_mpoly_zero(Lambda);

        /* unfortunate fmpz_mod_poly's store their own ctx :( */
        fmpz_mod_mpolyun_set_mod(Aeval, fpctx);
        fmpz_mod_mpolyun_set_mod(Beval, fpctx);
        fmpz_mod_mpolyun_set_mod(Geval, fpctx);
        fmpz_mod_mpolyun_set_mod(Abareval, fpctx);
        fmpz_mod_mpolyun_set_mod(Bbareval, fpctx);

    printf("Lambda: "); fmpz_mod_bma_mpoly_print(Lambda, Ictx); printf("\n");

        FLINT_ASSERT(fmpz_is_one(sshift));
        fmpz_mod_mpoly_bma_interpolate_alpha_powers(checkalpha, sshift, Ictx, ctx, fpctx);

        /* set evaluation of monomials */
        fmpz_mpoly_set_skel(Gammainc, Gamma, checkalpha, ctx, fpctx);
        fmpz_mpolyu_set_skel(Ainc, A, checkalpha, ctx, fpctx);
        fmpz_mpolyu_set_skel(Binc, B, checkalpha, ctx, fpctx);

        /* set reduction of coeffs */
        fmpz_mpoly_red_skel(Gammared, Gamma, fpctx);
        fmpz_mpolyu_red_skel(Ared, A, fpctx);
        fmpz_mpolyu_red_skel(Bred, B, fpctx);

        /* copy evaluation of monomials */
        fmpz_mpoly_copy_skel(Gammacur, Gammainc);
        fmpz_mpolyu_copy_skel(Acur, Ainc);
        fmpz_mpolyu_copy_skel(Bcur, Binc);

        fmpz_zero(image_count);

    next_bma_image:

    printf("next_bma_image "); fmpz_print(image_count); printf(" sshift: "); fmpz_print(sshift); printf("\n");
    printf("Lambda: "); fmpz_mod_bma_mpoly_print(Lambda, Ictx); printf("\n");
    usleep(1000000);

        /* image count is also the current power of alpha we are evaluating */
        fmpz_add_ui(image_count, image_count, 1);

    #if WANT_ASSERT
        /* image_count == sshift + Lambda->pointcount */
        fmpz_add_ui(t, sshift, Lambda->pointcount);
        FLINT_ASSERT(fmpz_equal(t, image_count));
    #endif

        fmpz_sub_ui(t, p, 1);
        if (fmpz_cmp(image_count, t) >= 0)
        {
            /* out of evaluation points alpha^image_count in Fp* */
            goto pick_bma_prime;
        }

        fmpz_mpoly_use_skel_mul(Gammaeval, Gammared, Gammacur, Gammainc, fpctx);
        fmpz_mpolyuu_use_skel_mul(Aeval, A, Ared, Acur, Ainc, ctx, fpctx);
        fmpz_mpolyuu_use_skel_mul(Beval, B, Bred, Bcur, Binc, ctx, fpctx);
    /*
    printf("Aeval: "); fmpz_mod_mpolyun_print_pretty(Aeval, NULL, Ictx->polyctx, Ictx->fpctx); printf("\n");
    printf("Beval: "); fmpz_mod_mpolyun_print_pretty(Beval, NULL, Ictx->polyctx, Ictx->fpctx); printf("\n");
    */
        if (Aeval->length == 0 || Beval->length == 0
            || fmpz_mod_mpolyun_bidegree(Aeval, ctx, fpctx) != A->exps[0]
            || fmpz_mod_mpolyun_bidegree(Beval, ctx, fpctx) != B->exps[0])
        {
            /* evaluation killed either lc(A) or lc(B) */
            fmpz_add_ui(sshift, sshift, Lambda->pointcount + 1);
            fmpz_mod_bma_mpoly_zero(Lambda);
            goto next_bma_image;
        }

        /* the evaluation killed neither lc(A) nor lc(B) */
        FLINT_ASSERT(!fmpz_is_zero(Gammaeval));

        success = fmpz_mod_mpolyun_gcd_bivar(Geval, Abareval, Bbareval, Aeval, Beval, ctx, fpctx);
        if (!success)
        {
            fmpz_add_ui(sshift, sshift, Lambda->pointcount + 1);
            fmpz_mod_bma_mpoly_zero(Lambda);
            goto next_bma_image;
        }

        FLINT_ASSERT(Geval->length > 0);
        GevaldegXY = fmpz_mod_mpolyun_bidegree(Geval, ctx, fpctx);
        fmpz_mod_mpolyun_scalar_mul_fmpz_mod(Geval, Gammaeval, ctx, fpctx);
    /*
    printf("Geval: "); fmpz_mod_mpolyun_print_pretty(Geval, NULL, Ictx->polyctx, Ictx->fpctx); printf("\n");
    */
        FLINT_ASSERT(fmpz_equal(Gammaeval, fmpz_mod_mpolyun_leadcoeff_last_ref(Geval, ctx, fpctx)));

        if (GdegboundXY < GevaldegXY)
        {
    printf("bma image unlucky\n");

            /* this image in Fp[X,Y] was unlucky */
            if (fmpz_equal(sshift, last_unlucky_sshift_plus_1))
            {
                /* this ksub is probably unlucky */
                goto pick_ksub;
            }
            if (++unlucky_count > 2)
            {
                goto pick_bma_prime;
            }
            fmpz_add_ui(last_unlucky_sshift_plus_1, sshift, 1);
            fmpz_add_ui(sshift, sshift, Lambda->pointcount + 1);
            fmpz_mod_bma_mpoly_zero(Lambda);
            goto next_bma_image;        
        }
        else if (GdegboundXY > GevaldegXY)
        {
    printf("bma image revealing\n");

            /* new bound on deg_XY(G) */
            fmpz_add_ui(sshift, sshift, Lambda->pointcount);
            fmpz_mod_bma_mpoly_zero(Lambda);
            fmpz_mod_bma_mpoly_add_point(Lambda, Geval, ctx, fpctx);
            GdegboundXY = GevaldegXY;
            if (GdegboundXY == 0)
            {
                fmpz_mpolyu_one(G, ctx);
                success = 1;
                goto cleanup;
            }
            goto next_bma_image;
        }

        fmpz_mod_bma_mpoly_add_point(Lambda, Geval, ctx, fpctx);
/*
    printf("Lambda: "); fmpz_mod_bma_mpoly_print(Lambda, Ictx); printf("\n");
*/
        if ((Lambda->pointcount & 1) != 0 || Gamma->length > Lambda->pointcount/2)
        {
            goto next_bma_image;
        }

        changed = fmpz_mod_bma_mpoly_reduce(Lambda);
/*
printf("changed: %d\n", changed);
    printf("Lambda: "); fmpz_mod_bma_mpoly_print(Lambda, Ictx); printf("\n");
*/
        if (changed)
        {
            goto next_bma_image;
        }

        success = fmpz_mod_bma_mpoly_get_mpolyu(H, sshift, Lambda, Ictx, ctx, fpctx);
        if (!success)
        {
printf("could not get mpolyu!!!!!\n");
usleep(1000000);
            goto next_bma_image;
        }

        if (H->length == 0 || (H->coeffs + 0)->length != Gamma->length)
        {
            goto next_bma_image;
        }

        /* GdegboundXY should be the bidegree of H */
        FLINT_ASSERT(GdegboundXY == H->exps[0]);

    printf("H: "); fmpz_mpolyuu_print_pretty(H, NULL, 2, ctx); printf("\n");

        switch (_random_check(&GevaldegXY, GdegboundXY,
                    Aeval, Beval, Geval, Abareval, Bbareval, Gammaeval,
                            checkalpha, H, A, B, Gamma, ctx, fpctx, randstate))
        {
            default:
                FLINT_ASSERT(0);
            case random_check_image_no_match:
            case random_check_image_degree_high:
                goto next_bma_image;
            case random_check_image_degree_low:
                /* the random evaluation point gave us a better degree bound */
                fmpz_add_ui(sshift, sshift, Lambda->pointcount);
                fmpz_mod_bma_mpoly_zero(Lambda);
                GdegboundXY = GevaldegXY;
                if (GdegboundXY == 0)
                {
                    fmpz_mpolyu_one(G, ctx);
                    success = 1;
                    goto cleanup;
                }
                goto next_bma_image;
            case random_check_point_not_found:
                /* hmmm */
            case random_check_good:
                NULL;
        }
    }

    /* assume that H is correct mod Hmodulus = p */
    fmpz_set(Hmodulus, p);
/*
printf("H("); fmpz_print(Hmodulus); printf("): "); fmpz_mpolyuu_print_pretty(H, NULL, 2, ctx); printf("\n");
*/
    /* find number of evals for zip interp */
    FLINT_ASSERT(H->length > 0);
    zip_evals = H->coeffs[0].length;
    for (i = 1; i < H->length; i++)
    {
        zip_evals = FLINT_MAX(zip_evals, H->coeffs[i].length);
    }
    zip_evals += 1; /* one extra check eval */
    nmod_zip_mpolyu_fit_poly(Z, H, zip_evals);

/*
flint_printf("zip_evals: %wd\n", zip_evals);
printf("Z: "); nmod_zip_mpolyuu_print(Z); printf("\n");
*/

    p_sp = UWORD(1) << (FLINT_BITS - 2);

pick_zip_prime:

    /*
        Get a new machine prime for zippel interpolation.
        H is currently interpolated modulo Hmodulus.
    */
    if (p_sp >= UWORD_MAX_PRIME)
    {
        /* ran out of machine primes - absolute failure */
        success = 0;
        goto cleanup;
    }
    p_sp = n_nextprime(p_sp, 1);

/*
flint_printf("zippel p: %wu\n", p_sp);
*/

    if (0 == fmpz_fdiv_ui(Hmodulus, p_sp))
    {
        goto pick_zip_prime;
    }

    nmod_mpoly_ctx_set_mod(ctx_sp, p_sp);
    /* unfortunate nmod_poly's need mod set */
    nmod_mpolyun_set_mod(Aeval_sp, ctx_sp->ffinfo->mod);
    nmod_mpolyun_set_mod(Beval_sp, ctx_sp->ffinfo->mod);
    nmod_mpolyun_set_mod(Geval_sp, ctx_sp->ffinfo->mod);
    nmod_mpolyun_set_mod(Abareval_sp, ctx_sp->ffinfo->mod);
    nmod_mpolyun_set_mod(Bbareval_sp, ctx_sp->ffinfo->mod);

    FLINT_ASSERT(p_sp > 3);
    for (i = 0; i < ctx->minfo->nvars; i++)
    {
        checkalpha_sp[i] = n_urandint(randstate, p_sp - 3) + 2;
/*flint_printf("alpha[%wd]: %wu\n", i, checkalpha_sp[i]);*/
    }

    nmod_zip_mpolyu_set_skel(Z, ctx_sp, H, checkalpha_sp, ctx);
/*
printf("Z: "); nmod_zip_mpolyuu_print(Z); printf("\n");
*/
    /* set evaluation of monomials */
    nmod_mpoly_set_skel(Gammainc_sp, ctx_sp, Gamma, checkalpha_sp, ctx);
    nmod_mpolyu_set_skel(Ainc_sp, ctx_sp, A, checkalpha_sp, ctx);
    nmod_mpolyu_set_skel(Binc_sp, ctx_sp, B, checkalpha_sp, ctx);

    /* set reduction of coeffs */
    nmod_mpoly_red_skel(Gammared_sp, Gamma, ctx_sp->ffinfo);
    nmod_mpolyu_red_skel(Ared_sp, A, ctx_sp->ffinfo);
    nmod_mpolyu_red_skel(Bred_sp, B, ctx_sp->ffinfo);

    /* copy evaluation of monomials */
    nmod_mpoly_copy_skel(Gammacur_sp, Gammainc_sp);
    nmod_mpolyu_copy_skel(Acur_sp, Ainc_sp);
    nmod_mpolyu_copy_skel(Bcur_sp, Binc_sp);

next_zip_image:

    Gammaeval_sp = nmod_mpoly_use_skel_mul(Gammared_sp, Gammacur_sp, Gammainc_sp, ctx_sp);
    nmod_mpolyuu_use_skel_mul(Aeval_sp, A, Ared_sp, Acur_sp, Ainc_sp, ctx_sp);
    nmod_mpolyuu_use_skel_mul(Beval_sp, B, Bred_sp, Bcur_sp, Binc_sp, ctx_sp);

    if (Aeval_sp->length == 0 || Beval_sp->length == 0
        || nmod_mpolyun_bidegree(Aeval_sp, ctx_sp) != A->exps[0]
        || nmod_mpolyun_bidegree(Beval_sp, ctx_sp) != B->exps[0])
    {
printf("zip image killed\n");
usleep(200000);
        /* evaluation point killed lc(A) or lc(B) */
        goto pick_zip_prime;
    }
/*
flint_printf("Gammaeval_sp: %wu\n", Gammaeval_sp);
printf("Aeval_sp: "); nmod_mpolyun_print_pretty(Aeval_sp, NULL, ctx_sp); printf("\n");
printf("Beval_sp: "); nmod_mpolyun_print_pretty(Beval_sp, NULL, ctx_sp); printf("\n");
*/
    /* the evaluation killed neither lc(A) nor lc(B) */
    FLINT_ASSERT(Gammaeval_sp != 0);

    success = nmod_mpolyun_gcd_brown_smprime_bivar(Geval_sp, Abareval_sp, Bbareval_sp, Aeval_sp, Beval_sp, ctx_sp);
    if (!success)
    {
printf("zip image failed\n");
usleep(200000);
        /* choose a bigger p */
        goto pick_zip_prime;
    }

    FLINT_ASSERT(Geval_sp->length > 0);
    GevaldegXY = nmod_mpolyun_bidegree(Geval_sp, ctx_sp);

    if (GevaldegXY > GdegboundXY)
    {
printf("zip image unlucky\n");
usleep(200000);
        /* this image in Fp'[X,Y] was unlucky */
        goto pick_zip_prime;        
    }
    else if (GevaldegXY < GdegboundXY)
    {
printf("zip image revealing\n");
usleep(200000);

        /* we have a new degree bound on deg_XY(G) */
        fmpz_add_ui(sshift, sshift, Lambda->pointcount);
        fmpz_mod_bma_mpoly_zero(Lambda);
        GdegboundXY = GevaldegXY;
        if (GdegboundXY == 0)
        {
            fmpz_mpolyu_one(G, ctx);
            success = 1;
            goto cleanup;
        }
        goto next_bma_image;
    }

    nmod_mpolyun_scalar_mul_nmod(Geval_sp, Gammaeval_sp, ctx_sp);
    FLINT_ASSERT(Gammaeval_sp == nmod_mpolyun_leadcoeff_last(Geval_sp, ctx_sp));
/*
printf("Geval_sp: "); nmod_mpolyun_print_pretty(Geval_sp, NULL, ctx_sp); printf("\n");
*/
    success = nmod_zip_mpolyuu_add_point(Z, Geval_sp);
    if (!success)
    {
printf("zip_image does not match\n");
usleep(200000);

        /*
            An image gcd in Fp'[X,Y] did not match the assumed formed in [X, Y].
            Start all over
        */
        goto pick_bma_prime;
    }
/*
printf("Z: "); nmod_zip_mpolyuu_print(Z); printf("\n");
*/
    if (Z->pointcount < zip_evals)
    {
        goto next_zip_image;
    }

    switch (nmod_mpolyu_zip_find_coeffs(Z, ctx_sp))
    {
        default:
            FLINT_ASSERT(0);
        case nmod_zip_find_coeffs_no_match:

printf("nmod_zip_find_coeffs_no_match\n");
usleep(200000);


            /*  The collection of image gcd's in Fp'[X,Y] could not be coerced
                into the assumed form in [X, Y][x_0, ..., x_(n-1)]. */
            goto pick_bma_prime;
        case nmod_zip_find_coeffs_non_invertible:

printf("nmod_zip_find_coeffs_non_invertible\n");
usleep(200000);


            /* The unlikely case where the evaluation points alpha produced
               a singular Vandermonde matrix. Assumed form is not nec wrong. */
            goto pick_zip_prime;
        case nmod_zip_find_coeffs_good:
            NULL;
    }

/*
printf("Z: "); nmod_zip_mpolyuu_print(Z); printf("\n");
*/
    FLINT_ASSERT(Hbits == H->bits);
    changed = fmpz_mpolyu_addinterp_zip(H, Hmodulus, Z, ctx_sp->ffinfo);
    fmpz_mul_ui(Hmodulus, Hmodulus, ctx_sp->ffinfo->mod.n);
/*
printf("H("); fmpz_print(Hmodulus); printf("): "); fmpz_mpolyuu_print_pretty(H, NULL, 2, ctx); printf("\n");
usleep(10000);
*/
    if (changed)
    {
        /* TODO if the coefficients of H are getting to large? */
        goto pick_zip_prime;
    }

/*
printf("zip unchanged\n");
usleep(100000);
*/

    /* compute content of H */
    fmpz_mpoly_set(Hcontent, H->coeffs + 0, ctx);
    FLINT_ASSERT(Hcontent->bits == Hbits);
    for (i = 1; i < H->length; i++)
    {
        success = _fmpz_mpoly_gcd(Hcontent, Hbits, Hcontent, H->coeffs + i, ctx);
        if (!success)
        {
            /* could not compute content - absolute failure */
            success = 0;
            goto cleanup;
        }
        FLINT_ASSERT(Hcontent->bits == Hbits);
    }
/*
printf("Hcontent: "); fmpz_mpoly_print_pretty(Hcontent, NULL, ctx); printf("\n");
*/
    /* upgrade Hpp to Hbits then try to pack down to bits */
    fmpz_mpolyu_set_bits(Hpp, Hbits, ctx);
    fmpz_mpolyu_divexact_mpoly(Hpp, H, Hcontent, ctx);
    success = fmpz_mpolyu_repack_bits(Hpp, bits, ctx);
    if (!success)
    {
        /* Hpp cannot be the GCD if it cannot be packed into bits */
        goto pick_zip_prime;
    }
    if (   !fmpz_mpolyuu_divides(Abar, A, Hpp, 2, ctx)
        || !fmpz_mpolyuu_divides(Bbar, B, Hpp, 2, ctx))
    {
        goto pick_zip_prime;
    }

    fmpz_mpolyu_swap(G, Hpp, ctx);
/*
printf("success!!!!!!!!!!!!!!!!!!!\n");
printf("G: "); fmpz_mpolyuu_print_pretty(G, NULL, 2, ctx); printf("\n");
*/
    success = 1;

cleanup:
/*
printf("cleanup\n");
*/
    fmpz_mpoly_clear(Hcontent, ctx);
    fmpz_mpolyu_clear(Hpp, ctx);
    fmpz_mpolyu_clear(H, ctx);
    fmpz_mpolyu_clear(Bbar, ctx);
    fmpz_mpolyu_clear(Abar, ctx);

    flint_free(Gdegbounds);
    flint_free(Adegs);
    flint_free(Bdegs);
    flint_free(Gammadegs);

    /* the zippler */
    nmod_zip_mpolyu_clear(Z);

    /* machine precision workspace */
    flint_free(checkalpha_sp);
    nmod_mpolyun_clear(Aeval_sp, ctx_sp);
    nmod_mpolyun_clear(Beval_sp, ctx_sp);
    nmod_mpolyun_clear(Geval_sp, ctx_sp);
    nmod_mpolyun_clear(Abareval_sp, ctx_sp);
    nmod_mpolyun_clear(Bbareval_sp, ctx_sp);
    nmod_mpolyc_clear(Gammainc_sp);
    nmod_mpolyc_clear(Gammacur_sp);
    nmod_mpolyc_clear(Gammared_sp);
    nmod_mpolycu_clear(Ainc_sp);
    nmod_mpolycu_clear(Acur_sp);
    nmod_mpolycu_clear(Ared_sp);
    nmod_mpolycu_clear(Binc_sp);
    nmod_mpolycu_clear(Bcur_sp);
    nmod_mpolycu_clear(Bred_sp);
    nmod_bma_mpoly_clear(Lambda_sp);
    nmod_mpoly_ctx_clear(ctx_sp);

    /* multiprecision workspace */
    for (i = 0; i < ctx->minfo->nvars; i++)
    {
        fmpz_clear(checkalpha + i);
    }
    flint_free(checkalpha);
    fmpz_mod_mpolyun_clear(Aeval, ctx, fpctx);
    fmpz_mod_mpolyun_clear(Beval, ctx, fpctx);
    fmpz_mod_mpolyun_clear(Geval, ctx, fpctx);
    fmpz_mod_mpolyun_clear(Abareval, ctx, fpctx);
    fmpz_mod_mpolyun_clear(Bbareval, ctx, fpctx);
    fmpz_mpolyc_clear(Gammainc);
    fmpz_mpolyc_clear(Gammacur);
    fmpz_mpolyc_clear(Gammared);
    fmpz_mpolycu_clear(Ainc);
    fmpz_mpolycu_clear(Acur);
    fmpz_mpolycu_clear(Ared);
    fmpz_mpolycu_clear(Binc);
    fmpz_mpolycu_clear(Bcur);
    fmpz_mpolycu_clear(Bred);
    fmpz_mod_bma_mpoly_clear(Lambda);
    fmpz_mod_ctx_clear(fpctx);

    mpoly_bma_interpolate_ctx_clear(Ictx);

    fmpz_clear(Hmodulus);
    fmpz_clear(Gammaeval);
    fmpz_clear(last_unlucky_sshift_plus_1);
    fmpz_clear(sshift);
    fmpz_clear(cBksub);
    fmpz_clear(cAksub);
    fmpz_clear(subprod);
    fmpz_clear(image_count);
    fmpz_clear(shift);
    fmpz_clear(t);
    fmpz_clear(p);
    flint_randclear(randstate);

    if (success)
    {
        FLINT_ASSERT(G->bits == bits);
    }

    return success;
}




int fmpz_mpoly_gcd_bma(
    fmpz_mpoly_t G,
    const fmpz_mpoly_t A,
    const fmpz_mpoly_t B,
    const fmpz_mpoly_ctx_t ctx)
{
    slong i, j;
    mp_bitcnt_t Gbits, ABbits;
    int success = 0;
    fmpz_mpoly_ctx_t uctx;
    fmpz_mpolyu_t Auu, Buu, Guu, Abar, Bbar, Gbar;
    fmpz_mpoly_t Acontent, Bcontent, Gamma;
    slong * Adegs, * Bdegs, * perm;
    ulong * shift, * stride;
    ulong max_main_degree, max_minor_degree;

timeit_t time;

    if (fmpz_mpoly_is_zero(A, ctx))
    {
        if (fmpz_mpoly_is_zero(B, ctx))
        {
            fmpz_mpoly_zero(G, ctx);
            return 1;
        }
        if (fmpz_sgn(B->coeffs + 0) < 0)
        {
            fmpz_mpoly_neg(G, B, ctx);
            return 1;
        } else
        {
            fmpz_mpoly_set(G, B, ctx);
            return 1;
        }
    }

    if (fmpz_mpoly_is_zero(B, ctx))
    {
        if (fmpz_sgn(A->coeffs + 0) < 0)
        {
            fmpz_mpoly_neg(G, A, ctx);
            return 1;
        } else
        {
            fmpz_mpoly_set(G, A, ctx);
            return 1;
        }
    }

    if (A->bits > FLINT_BITS || B->bits > FLINT_BITS)
    {
        return 0;
    }

    if (ctx->minfo->nvars < 3)
    {
        return fmpz_mpoly_gcd_zippel(G, A, B, ctx);
    }

    FLINT_ASSERT(A->bits <= FLINT_BITS);
    FLINT_ASSERT(B->bits <= FLINT_BITS);
    FLINT_ASSERT(ctx->minfo->nvars >= 3);
    FLINT_ASSERT(!fmpz_mpoly_is_zero(A, ctx));
    FLINT_ASSERT(!fmpz_mpoly_is_zero(B, ctx));


timeit_start(time);


    Adegs = (slong *) flint_malloc(ctx->minfo->nvars*sizeof(slong));
    Bdegs = (slong *) flint_malloc(ctx->minfo->nvars*sizeof(slong));
    perm = (slong *) flint_malloc((ctx->minfo->nvars)*sizeof(slong));
    shift = (ulong *) flint_malloc((ctx->minfo->nvars)*sizeof(ulong));
    stride = (ulong *) flint_malloc((ctx->minfo->nvars)*sizeof(ulong));

    mpoly_degrees_si(Adegs, A->exps, A->length, A->bits, ctx->minfo);
    mpoly_degrees_si(Bdegs, B->exps, B->length, B->bits, ctx->minfo);

    max_main_degree = 0;
    max_minor_degree = 0;
    for (i = 0; i < ctx->minfo->nvars; i++)
    {
        perm[i] = i;
        shift[i] = 0;
        stride[i] = 1;
        FLINT_ASSERT(Adegs[i] >= 0);
        FLINT_ASSERT(Bdegs[i] >= 0);
        if (i < 2)
        {
            max_main_degree = FLINT_MAX(max_main_degree, Adegs[i]);
            max_main_degree = FLINT_MAX(max_main_degree, Bdegs[i]);
        }
        else
        {
            max_minor_degree = FLINT_MAX(max_minor_degree, Adegs[i]);
            max_minor_degree = FLINT_MAX(max_minor_degree, Bdegs[i]);
        }
    }

    fmpz_mpoly_ctx_init(uctx, ctx->minfo->nvars - 2, ORD_LEX);

    /* ABbits is bits for intermediates in ZZ[x_0,x_1][x_2,...,x_(n-1)] */
    ABbits = 1 + FLINT_BIT_COUNT(max_minor_degree);
    ABbits = FLINT_MAX(MPOLY_MIN_BITS, ABbits);
    ABbits = mpoly_fix_bits(ABbits, uctx->minfo);
    FLINT_ASSERT(ABbits <= FLINT_BITS);

    fmpz_mpolyu_init(Auu, ABbits, uctx);
    fmpz_mpolyu_init(Buu, ABbits, uctx);
    fmpz_mpolyu_init(Guu, ABbits, uctx);
    fmpz_mpoly_init3(Acontent, 0, ABbits, uctx);
    fmpz_mpoly_init3(Bcontent, 0, ABbits, uctx);
    fmpz_mpoly_init3(Gamma, 0, ABbits, uctx);
    fmpz_mpolyu_init(Abar, ABbits, uctx);
    fmpz_mpolyu_init(Bbar, ABbits, uctx);
    fmpz_mpolyu_init(Gbar, ABbits, uctx);

    /* two main variables must be packed into bits = FLINT_BITS/2 */
    if (FLINT_BIT_COUNT(max_main_degree) >= FLINT_BITS/2)
    {
        success = 0;
        goto cleanup;
    }

    /* Gbits is bits for final answer in ZZ[x_0,...,x_(n-1)] */
    Gbits = FLINT_MIN(A->bits, B->bits);

    fmpz_mpoly_to_mpolyuu_perm_deflate(Auu, perm, shift, stride, uctx, A, ctx);
    fmpz_mpoly_to_mpolyuu_perm_deflate(Buu, perm, shift, stride, uctx, B, ctx);


timeit_stop(time);
flint_printf("bma setup time: %wu\n", time->cpu);

printf("Auu: "); fmpz_mpolyuu_print_pretty(Auu, NULL, 2, uctx); printf("\n");
printf("Buu: "); fmpz_mpolyuu_print_pretty(Buu, NULL, 2, uctx); printf("\n");


timeit_start(time);
    /* compute content of A */
    if (Auu->length < 1)
    {
        FLINT_ASSERT(Auu->length == 1);
        fmpz_mpoly_set(Acontent, Auu->coeffs + 0, uctx);
        
    }
    else
    {
        j = 0;
        for (i = 1; i < Auu->length; i++)
        {
            if ((Auu->coeffs + i)->length < (Auu->coeffs + j)->length)
            {
                j = i;
            }
        }

        if (j == 0)
        {
            success = _fmpz_mpoly_gcd(Acontent, ABbits, Auu->coeffs + 0, Auu->coeffs + 1, uctx);
            if (!success)
                goto cleanup;
            FLINT_ASSERT(Acontent->bits == ABbits);
            for (i = 2; i < Auu->length; i++)
            {
                success = _fmpz_mpoly_gcd(Acontent, ABbits, Acontent, Auu->coeffs + i, uctx);
                if (!success)
                    goto cleanup;
                FLINT_ASSERT(Acontent->bits == ABbits);
            }
        }
        else
        {
            success = _fmpz_mpoly_gcd(Acontent, ABbits, Auu->coeffs + 0, Auu->coeffs + j, uctx);
            if (!success)
                goto cleanup;
            FLINT_ASSERT(Acontent->bits == ABbits);
            for (i = 1; i < Auu->length; i++)
            {
                if (i == j)
                    continue;
                success = _fmpz_mpoly_gcd(Acontent, ABbits, Acontent, Auu->coeffs + i, uctx);
                if (!success)
                    goto cleanup;
                FLINT_ASSERT(Acontent->bits == ABbits);
            }
        }
    }

timeit_stop(time);
flint_printf("a content time: %wu\n", time->cpu);


timeit_start(time);
    /* compute content of B */
    fmpz_mpoly_set(Bcontent, Buu->coeffs + 0, uctx);
    FLINT_ASSERT(Bcontent->bits == ABbits);
    for (i = 1; i < Buu->length; i++)
    {
        success = _fmpz_mpoly_gcd(Bcontent, ABbits, Bcontent, Buu->coeffs + i, uctx);
        if (!success)
            goto cleanup;
        FLINT_ASSERT(Bcontent->bits == ABbits);
    }

timeit_stop(time);
flint_printf("b content time: %wu\n", time->cpu);


    /* remove content from A and B */
    fmpz_mpolyu_divexact_mpoly(Abar, Auu, Acontent, uctx);
    fmpz_mpolyu_divexact_mpoly(Bbar, Buu, Bcontent, uctx);

    /* compute GCD of leading coefficients */
    _fmpz_mpoly_gcd(Gamma, ABbits, Abar->coeffs + 0, Bbar->coeffs + 0, uctx);




timeit_start(time);
    success = fmpz_mpolyuu_gcd_bma(Gbar, Abar, Bbar, Gamma, uctx);
timeit_stop(time);
flint_printf("bma time: %wu\n", time->cpu);


    if (!success)
        goto cleanup;
    FLINT_ASSERT(Gamma->bits == ABbits);

/*
printf("Gbar: "); fmpz_mpolyuu_print_pretty(Gbar, NULL, 2, uctx); printf("\n");
printf("Acontent: "); fmpz_mpoly_print_pretty(Acontent, NULL, uctx); printf("\n");
printf("Bcontent: "); fmpz_mpoly_print_pretty(Bcontent, NULL, uctx); printf("\n");
*/

    /* put back content */
    success = _fmpz_mpoly_gcd(Acontent, ABbits, Acontent, Bcontent, uctx);
    if (!success)
        goto cleanup;
    FLINT_ASSERT(Acontent->bits == ABbits);
/*
printf("Gcontent: "); fmpz_mpoly_print_pretty(Acontent, NULL, uctx); printf("\n");
*/

    fmpz_mpolyu_mul_mpoly(Guu, Gbar, Acontent, uctx);

/*
printf("Guu: "); fmpz_mpolyuu_print_pretty(Guu, NULL, 2, uctx); printf("\n");
*/

    fmpz_mpoly_from_mpolyuu_perm_inflate(G, Gbits, ctx, Guu, perm, shift, stride, uctx);
/*
printf("G: "); fmpz_mpoly_print_pretty(G, NULL, ctx); printf("\n");
*/

    if (fmpz_sgn(G->coeffs + 0) < 0)
        fmpz_mpoly_neg(G, G, ctx);
    FLINT_ASSERT(G->bits == Gbits);

    success = 1;

cleanup:

    flint_free(Adegs);
    flint_free(Bdegs);
    flint_free(perm);
    flint_free(shift);
    flint_free(stride);

    fmpz_mpolyu_clear(Abar, uctx);
    fmpz_mpolyu_clear(Bbar, uctx);
    fmpz_mpolyu_clear(Gbar, uctx);
    fmpz_mpoly_clear(Acontent, uctx);
    fmpz_mpoly_clear(Bcontent, uctx);
    fmpz_mpoly_clear(Gamma, uctx);

    fmpz_mpolyu_clear(Auu, uctx);
    fmpz_mpolyu_clear(Buu, uctx);
    fmpz_mpolyu_clear(Guu, uctx);
    fmpz_mpoly_ctx_clear(uctx);

    return success;
}


int
main(void)
{
    slong i, j;
    FLINT_TEST_INIT(state);

    flint_printf("gcd_bma....\n");
    fflush(stdout);

    if (1)
    {
timeit_t time;
        fmpz_mpoly_ctx_t ctx;
        fmpz_mpoly_t G, A, B;
        const char* vars[] = {"y", "t", "x", "z"};

        fmpz_mpoly_ctx_init(ctx, 4, ORD_LEX);
        fmpz_mpoly_init(G, ctx);
        fmpz_mpoly_init(A, ctx);
        fmpz_mpoly_init(B, ctx);

        fmpz_mpoly_set_str_pretty(A, "39 - t*x - 7*x^2*y^3*z^11 + x^1000*y^3*z^20", vars, ctx);
        fmpz_mpoly_set_str_pretty(B, "1 + x^100 + x^3*y + 2*t^15*x^78*y^3*z^13", vars, ctx);
        fmpz_mpoly_set_str_pretty(G, "39 - t*x + 39*x^100 - t*x^101 + 39*x^3*y - t*x^4*y - 7*x^2*y^3*z^11"
      " - 7*x^102*y^3*z^11 - 7*x^5*y^4*z^11 + 78*t^15*x^78*y^3*z^13 - 2*t^16*x^79*y^3*z^13 + x^1000*y^3*z^20"
         " + x^1100*y^3*z^20 + x^1003*y^4*z^20 - 14*t^15*x^80*y^6*z^24 + 2*t^15*x^1078*y^6*z^33", vars, ctx);
printf("A: "); fmpz_mpoly_print_pretty(A, vars, ctx); printf("\n");
printf("B: "); fmpz_mpoly_print_pretty(B, vars, ctx); printf("\n");
printf("G: "); fmpz_mpoly_print_pretty(G, vars, ctx); printf("\n");
        fmpz_mpoly_mul(A, A, G, ctx);
        fmpz_mpoly_mul(B, B, G, ctx);

timeit_start(time);
        fmpz_mpoly_gcd_bma(G, A, B, ctx);
timeit_stop(time);
flint_printf("time: %wd\n", time->cpu);

timeit_start(time);
        fmpz_mpoly_gcd_bma(G, A, B, ctx);
timeit_stop(time);
flint_printf("time: %wd\n", time->cpu);


printf("A: "); fmpz_mpoly_print_pretty(A, vars, ctx); printf("\n");
printf("B: "); fmpz_mpoly_print_pretty(B, vars, ctx); printf("\n");
printf("G: "); fmpz_mpoly_print_pretty(G, vars, ctx); printf("\n");

        fmpz_mpoly_clear(A, ctx);
        fmpz_mpoly_clear(B, ctx);
        fmpz_mpoly_clear(G, ctx);
        fmpz_mpoly_ctx_clear(ctx);

    }

    for (i = 0; i < 0 * flint_test_multiplier(); i++)
    {
        fmpz_mpoly_ctx_t ctx;
        fmpz_mpoly_t a, b, g, ca, cb, cg, t;
        mp_bitcnt_t coeff_bits;
        slong len, len1, len2;
        slong degbound;
        int res;

        fmpz_mpoly_ctx_init_rand(ctx, state, 5);

        fmpz_mpoly_init(g, ctx);
        fmpz_mpoly_init(a, ctx);
        fmpz_mpoly_init(b, ctx);
        fmpz_mpoly_init(ca, ctx);
        fmpz_mpoly_init(cb, ctx);
        fmpz_mpoly_init(cg, ctx);
        fmpz_mpoly_init(t, ctx);

        len = n_randint(state, 10) + 1;
        len1 = n_randint(state, 20);
        len2 = n_randint(state, 20);

        degbound = 30/(2*ctx->minfo->nvars - 1);

        coeff_bits = n_randint(state, 100);

        for (j = 0; j < 4; j++)
        {
            do {
                fmpz_mpoly_randtest_bound(t, state, len, coeff_bits + 1, degbound, ctx);
            } while (t->length == 0);
            fmpz_mpoly_randtest_bound(a, state, len1, coeff_bits, degbound, ctx);
            fmpz_mpoly_randtest_bound(b, state, len2, coeff_bits, degbound, ctx);
            fmpz_mpoly_mul(a, a, t, ctx);
            fmpz_mpoly_mul(b, b, t, ctx);

            fmpz_mpoly_randtest_bits(g, state, len, coeff_bits, FLINT_BITS, ctx);

            res = fmpz_mpoly_gcd_bma(g, a, b, ctx);
            if (!res)
            {
                printf("FAIL\n");
                flint_printf("Check gcd could be computed\ni = %wd, j = %wd\n", i ,j);
                flint_abort();
            }

            fmpz_mpoly_assert_canonical(g, ctx);

            if (fmpz_mpoly_is_zero(g, ctx))
            {
                if (!fmpz_mpoly_is_zero(a, ctx) || !fmpz_mpoly_is_zero(b, ctx)) {
                    printf("FAIL\n");
                    flint_printf("Check zero gcd only results from zero inputs\ni = %wd, j = %wd\n", i ,j);
                    flint_abort();
                }
                continue;
            }

            if (fmpz_sgn(g->coeffs + 0) <= 0)
            {
                printf("FAIL\n");
                flint_printf("Check gcd has positive lc\ni = %wd, j = %wd\n", i ,j);
                flint_abort();
            }

            res = 1;
            res = res && fmpz_mpoly_divides(ca, a, g, ctx);
            res = res && fmpz_mpoly_divides(cb, b, g, ctx);
            if (!res)
            {
                printf("FAIL\n");
                flint_printf("Check divisibility\ni = %wd, j = %wd\n", i ,j);
                flint_abort();
            }

            res = fmpz_mpoly_gcd_bma(cg, ca, cb, ctx);

            if (!res)
            {
                printf("FAIL\n");
                flint_printf("Check cofactor gcd could be computed\ni = %wd, j = %wd\n", i ,j);
                flint_abort();
            }

            if (!fmpz_mpoly_equal_ui(cg, UWORD(1), ctx))
            {
                printf("FAIL\n");
                flint_printf("Check cofactors are relatively prime\ni = %wd, j = %wd\n", i ,j);                
                flint_abort();
            }
        }

        fmpz_mpoly_clear(g, ctx);
        fmpz_mpoly_clear(a, ctx);
        fmpz_mpoly_clear(b, ctx);
        fmpz_mpoly_clear(ca, ctx);
        fmpz_mpoly_clear(cb, ctx);
        fmpz_mpoly_clear(cg, ctx);
        fmpz_mpoly_clear(t, ctx);
        fmpz_mpoly_ctx_clear(ctx);
    }

    printf("PASS\n");
    FLINT_TEST_CLEANUP(state);

    return 0;
}
