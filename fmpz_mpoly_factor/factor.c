/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fmpz_mpoly_factor.h"
#include "fmpq_poly.h"
#include "fmpz_mod_mpoly.h"
#include "nmod_mpoly_factor.h"
#include "mpn_extras.h"
#include "nmod_vec.h"

#define LOW_THIRD_MASK ((-UWORD(1)) >> (FLINT_BITS - FLINT_BITS/3))

const char * ourvars[] = {"x", "y", "z", "w", "u", "v"};

void usleep(ulong);



/*************** n_polyun ***************************/

void n_polyun_sort_terms(n_polyun_t A)
{
    slong i, j;
    for (i = 1; i < A->length; i++)
    for (j = i; j > 0 && A->terms[j - 1].exp < A->terms[j].exp; j--)
        n_polyun_term_swap(A->terms + j - 1, A->terms + j);
    return;
}

void n_polyu_sort_terms(n_polyu_t A)
{
    slong i, j;
    for (i = 1; i < A->length; i++)
    for (j = i; j > 0 && A->terms[j - 1].exp < A->terms[j].exp; j--)
        n_polyu_term_swap(A->terms + j - 1, A->terms + j);
    return;
}

void n_polyun_mod_combine_like_terms(n_polyun_t A, nmod_t mod)
{
    slong in, out;
    n_polyun_term_struct * Aterms = A->terms;

    out = -1;

    for (in = 0; in < A->length; in++)
    {
        FLINT_ASSERT(in > out);

        if (out >= 0 && Aterms[out].exp == Aterms[in].exp)
        {
            n_poly_mod_add(Aterms[out].coeff, Aterms[out].coeff,
                                              Aterms[in].coeff, mod);
        }
        else
        {
            if (out < 0 || !n_poly_is_zero(Aterms[out].coeff))
                out++;

            if (out != in)
            {
                Aterms[out].exp = Aterms[in].exp;
                n_poly_swap(Aterms[out].coeff, Aterms[in].coeff);
            }
        }
    }

    if (out < 0 || !n_poly_is_zero(Aterms[out].coeff))
        out++;

    A->length = out;
}

void n_polyu_mod_combine_like_terms(n_polyu_t A, nmod_t mod)
{
    slong in, out;
    n_polyu_term_struct * Aterms = A->terms;

    out = -1;

    for (in = 0; in < A->length; in++)
    {
        FLINT_ASSERT(in > out);

        if (out >= 0 && Aterms[out].exp == Aterms[in].exp)
        {
            Aterms[out].coeff = nmod_add(Aterms[out].coeff,
                                         Aterms[in].coeff, mod);
        }
        else
        {
            if (out < 0 || Aterms[out].coeff != 0)
                out++;

            if (out != in)
            {
                Aterms[out].exp = Aterms[in].exp;
                Aterms[out].coeff = Aterms[in].coeff;
            }
        }
    }

    if (out < 0 || Aterms[out].coeff != 0)
        out++;

    A->length = out;
}

void n_polyu_get_n_polyun(n_polyu_t A, const n_polyun_t B)
{
    slong i, j;
    A->length = 0;
    for (i = 0; i < B->length; i++)
    {
        for (j = B->terms[i].coeff->length - 1; j >= 0; j--)
        {
            if (B->terms[i].coeff->coeffs[j] == 0)
                continue;
            n_polyu_fit_length(A, A->length + 1);
            A->terms[A->length].coeff = B->terms[i].coeff->coeffs[j];
            A->terms[A->length].exp = B->terms[i].exp + j;
            A->length++;
        }
    }
}

void n_polyun_mod_mul(
    n_polyun_t A,
    const n_polyun_t B,
    const n_polyun_t C,
    nmod_t mod)
{
    slong Ai, Bi, Ci;

    FLINT_ASSERT(A != B);
    FLINT_ASSERT(A != C);

    n_polyun_fit_length(A, B->length*C->length);
    Ai = 0;
    for (Bi = 0; Bi < B->length; Bi++)
    for (Ci = 0; Ci < C->length; Ci++)
    {
        A->terms[Ai].exp = B->terms[Bi].exp + C->terms[Ci].exp;
        n_poly_mod_mul(A->terms[Ai].coeff, B->terms[Bi].coeff,
                                           C->terms[Ci].coeff, mod);
        Ai++;
    }
    A->length = Ai;
    n_polyun_sort_terms(A);
    n_polyun_mod_combine_like_terms(A, mod);
}

void n_polyu_mod_mul(
    n_polyu_t A,
    const n_polyu_t B,
    const n_polyu_t C,
    nmod_t mod)
{
    slong Ai, Bi, Ci;

    FLINT_ASSERT(A != B);
    FLINT_ASSERT(A != C);

    n_polyu_fit_length(A, B->length*C->length);
    Ai = 0;
    for (Bi = 0; Bi < B->length; Bi++)
    for (Ci = 0; Ci < C->length; Ci++)
    {
        A->terms[Ai].exp = B->terms[Bi].exp + C->terms[Ci].exp;
        A->terms[Ai].coeff = nmod_mul(B->terms[Bi].coeff,
                                      C->terms[Ci].coeff, mod);
        Ai++;
    }
    A->length = Ai;
    n_polyu_sort_terms(A);
    n_polyu_mod_combine_like_terms(A, mod);
}

int n_polyu_equal(n_polyu_t A, n_polyu_t B)
{
    slong i;
    if (A->length != B->length)
        return 0;
    for (i = 0; i < A->length; i++)
    {
        if (A->terms[i].coeff != B->terms[i].coeff)
            return 0;
        if (A->terms[i].exp != B->terms[i].exp)
            return 0;
    }
    return 1;
}


int nmod_mpolyu_is_canonical(
    const nmod_mpolyu_t A,
    const nmod_mpoly_ctx_t ctx)
{
    slong i;
    for (i = 0; i < A->length; i++)
    {
        if (!nmod_mpoly_is_canonical(A->coeffs + i, ctx))
            return 0;
        if (nmod_mpoly_is_zero(A->coeffs + i, ctx))
            return 0;
        if (i > 0 && A->exps[i] >= A->exps[i - 1])
            return 0;
    }
    return 1;
}


void nmod_mpolyu2_print_pretty(
    nmod_mpolyu_t A,
    const char * var0,
    const char * var1,
    const char ** vars,
    const nmod_mpoly_ctx_t ctx)
{
    slong i;
    int first = 1;

    for (i = 0; i < A->length; i++)
    {
        if (!first)
            printf(" + ");
        first = 0;
        flint_printf("(");
        nmod_mpoly_print_pretty(A->coeffs + i, vars, ctx);
        flint_printf(")*%s^%wu*%s^%wu",
            var0, extract_exp(A->exps[i], 1, 2),
            var1, extract_exp(A->exps[i], 0, 2));
    }

    if (first)
        flint_printf("0");
}

void nmod_mpolyu3_print_pretty(
    const nmod_mpolyu_t A,
    const char * var0,
    const char * var1,
    const char * var2,
    const char ** vars,
    const nmod_mpoly_ctx_t ctx)
{
    slong i;
    int first = 1;

    for (i = 0; i < A->length; i++)
    {
        if (!first)
            printf(" + ");
        first = 0;
        flint_printf("(");
        nmod_mpoly_print_pretty(A->coeffs + i, vars, ctx);
        flint_printf(")*%s^%wu*%s^%wu*%s^%wu",
            var0, extract_exp(A->exps[i], 2, 3),
            var1, extract_exp(A->exps[i], 1, 3),
            var2, extract_exp(A->exps[i], 0, 3));
    }

    if (first)
        flint_printf("0");
}


static void nmod_mpoly_get_mpolyu2(
    nmod_mpolyu_t A,
    const nmod_mpoly_t B,
    slong var0,
    slong var1,
    const nmod_mpoly_ctx_t ctx)
{
    slong i;
    nmod_mpoly_struct * Ac;
    ulong * Bexps;
    slong NA, NB;

    FLINT_ASSERT(ctx->minfo->ord == ORD_LEX);

    FLINT_ASSERT(var0 < ctx->minfo->nvars);
    FLINT_ASSERT(var1 < ctx->minfo->nvars);

    Bexps = (ulong *) flint_malloc(ctx->minfo->nvars*sizeof(ulong));

    NA = mpoly_words_per_exp_sp(A->bits, ctx->minfo);
    NB = mpoly_words_per_exp(B->bits, ctx->minfo);

    A->length = 0;
    for (i = 0; i < B->length; i++)
    {
        mpoly_get_monomial_ui(Bexps, B->exps + NB*i, B->bits, ctx->minfo);
        FLINT_ASSERT(FLINT_BIT_COUNT(Bexps[var0]) < FLINT_BITS/2);
        FLINT_ASSERT(FLINT_BIT_COUNT(Bexps[var1]) < FLINT_BITS/2);
        Ac = _nmod_mpolyu_get_coeff(A, pack_exp2(Bexps[var0], Bexps[var1]), ctx);
        FLINT_ASSERT(Ac->bits == A->bits);
        nmod_mpoly_fit_length(Ac, Ac->length + 1, ctx);
        Ac->coeffs[Ac->length] = B->coeffs[i];
        Bexps[var0] = 0;
        Bexps[var1] = 0;
        mpoly_set_monomial_ui(Ac->exps + NA*Ac->length, Bexps, A->bits, ctx->minfo);
        Ac->length++;
    }

    flint_free(Bexps);
/*
printf("nmod_mpoly_get_mpolyu3 returning: "); nmod_mpolyu3_print_pretty(A, ourvars[var0], ourvars[var1], ourvars[var2], ourvars, ctx); printf("\n");
*/
    FLINT_ASSERT(nmod_mpolyu_is_canonical(A, ctx));
}


/*
    B vars: x0 x1 x2 x3 x4 xv           2 < v
    A vars: xv x0 x1 : 0 0 x2 x3 x4 0
*/
static void nmod_mpoly_get_mpolyu3(
    nmod_mpolyu_t A,
    const nmod_mpoly_t B,
    slong var0,
    slong var1,
    slong var2,
    const nmod_mpoly_ctx_t ctx)
{
    slong i;
    nmod_mpoly_struct * Ac;
    ulong * Bexps;
    slong NA, NB;

    FLINT_ASSERT(ctx->minfo->ord == ORD_LEX);

    FLINT_ASSERT(var0 < ctx->minfo->nvars);
    FLINT_ASSERT(var1 < ctx->minfo->nvars);
    FLINT_ASSERT(var2 < ctx->minfo->nvars);

    Bexps = (ulong *) flint_malloc(ctx->minfo->nvars*sizeof(ulong));

    NA = mpoly_words_per_exp_sp(A->bits, ctx->minfo);
    NB = mpoly_words_per_exp(B->bits, ctx->minfo);

    A->length = 0;
    for (i = 0; i < B->length; i++)
    {
        mpoly_get_monomial_ui(Bexps, B->exps + NB*i, B->bits, ctx->minfo);
        FLINT_ASSERT(FLINT_BIT_COUNT(Bexps[var0]) < FLINT_BITS/3);
        FLINT_ASSERT(FLINT_BIT_COUNT(Bexps[var1]) < FLINT_BITS/3);
        FLINT_ASSERT(FLINT_BIT_COUNT(Bexps[var2]) < FLINT_BITS/3);
        Ac = _nmod_mpolyu_get_coeff(A, pack_exp3(Bexps[var0], Bexps[var1], Bexps[var2]), ctx);
        FLINT_ASSERT(Ac->bits == A->bits);
        nmod_mpoly_fit_length(Ac, Ac->length + 1, ctx);
        Ac->coeffs[Ac->length] = B->coeffs[i];
        Bexps[var0] = 0;
        Bexps[var1] = 0;
        Bexps[var2] = 0;
        mpoly_set_monomial_ui(Ac->exps + NA*Ac->length, Bexps, A->bits, ctx->minfo);
        Ac->length++;
    }

    flint_free(Bexps);
/*
printf("nmod_mpoly_get_mpolyu3 returning: "); nmod_mpolyu3_print_pretty(A, ourvars[var0], ourvars[var1], ourvars[var2], ourvars, ctx); printf("\n");
*/
    FLINT_ASSERT(nmod_mpolyu_is_canonical(A, ctx));
}


/*
    fmpz_mod_mpoly_t
    sparse multivariates with fmpz_mod coeffs
*/
typedef struct
{
   fmpz * coeffs;
   ulong * exps;  
   slong alloc;
   slong length;
   flint_bitcnt_t bits;     /* number of bits per exponent */
} fmpz_mod_mpoly_struct;

typedef fmpz_mod_mpoly_struct fmpz_mod_mpoly_t[1];

void fmpz_mod_mpoly_clear(
    fmpz_mod_mpoly_t poly,
    const fmpz_mod_mpoly_ctx_t ctx)
{
   if (poly->coeffs != NULL)
   {
      slong i;

      for (i = 0; i < poly->alloc; i++)
         _fmpz_demote(poly->coeffs + i);

      flint_free(poly->coeffs);
      flint_free(poly->exps);
   }
}


void fmpz_mod_mpoly_init(
    fmpz_mod_mpoly_t A,
    const fmpz_mod_mpoly_ctx_t ctx)
{
    A->coeffs = NULL;
    A->exps = NULL;
    A->alloc = 0;
    A->length = 0;
    A->bits = MPOLY_MIN_BITS;
}

void fmpz_mod_mpoly_init3(
    fmpz_mod_mpoly_t A,
    slong alloc,
    flint_bitcnt_t bits,
    const fmpz_mod_mpoly_ctx_t ctx)
{
    if (alloc > 0)
    {
        slong N = mpoly_words_per_exp(bits, ctx->minfo);
        A->exps   = (ulong *) flint_malloc(alloc*N*sizeof(ulong));
        A->coeffs = (fmpz *) flint_calloc(alloc, sizeof(fmpz));
    }
    else
    {
        alloc = 0;
        A->coeffs = NULL;
        A->exps = NULL;
    }
    A->alloc = alloc;
    A->length = 0;
    A->bits = bits;
}

void fmpz_mod_mpoly_truncate(
    fmpz_mod_mpoly_t A,
    slong newlen,
    const fmpz_mod_mpoly_ctx_t ctx)
{
    if (A->length > newlen)
    {
        slong i;

        for (i = newlen; i < A->length; i++)
            _fmpz_demote(A->coeffs + i);

        A->length = newlen;
    }  
}

void fmpz_mod_mpoly_realloc(
    fmpz_mod_mpoly_t poly,
    slong alloc,
    const fmpz_mod_mpoly_ctx_t ctx)
{
    slong N;

    if (alloc <= 0)             /* Clear up, reinitialise */
    {
        fmpz_mod_mpoly_clear(poly, ctx);
        fmpz_mod_mpoly_init(poly, ctx);
        return;
    }

    N = mpoly_words_per_exp(poly->bits, ctx->minfo);

    if (poly->alloc != 0)            /* Realloc */
    {
        fmpz_mod_mpoly_truncate(poly, alloc, ctx);

        poly->coeffs = (fmpz *) flint_realloc(poly->coeffs, alloc*sizeof(fmpz));
        poly->exps = (ulong *) flint_realloc(poly->exps, alloc*N*sizeof(ulong));

        if (alloc > poly->alloc)
            memset(poly->coeffs + poly->alloc, 0,
                                           (alloc - poly->alloc)*sizeof(fmpz));
    }
    else                        /* Nothing allocated already so do it now */
    {
        poly->coeffs = (fmpz *) flint_calloc(alloc, sizeof(fmpz));
        poly->exps   = (ulong *) flint_malloc(alloc*N*sizeof(ulong));
    }

    poly->alloc = alloc;
}

void fmpz_mod_mpoly_fit_length(
    fmpz_mod_mpoly_t poly,
    slong len,
    const fmpz_mod_mpoly_ctx_t ctx)
{
    if (len > poly->alloc)
    {
        /* At least double number of allocated coeffs */
        if (len < 2 * poly->alloc)
            len = 2 * poly->alloc;
        fmpz_mod_mpoly_realloc(poly, len, ctx);
    }
}



    


/******* n_bpoly - (Z/nZ)[x,y] ****************************************/

void nmod_mpoly_to_bpoly(
    n_bpoly_t A,
    const nmod_mpoly_t B,
    slong varx,
    slong vary,
    const nmod_mpoly_ctx_t ctx);


/* M is a d x d array */
void computeM(
    mp_limb_t * M,
    slong d,
    const nmod_t ctx)
{
    slong i, j;
    mp_limb_t i2, i2j;

    for (i = 1; i <= d; i++)
    {
        i2 = nmod_mul(i, i, ctx);
        i2j = 1;
        for (j = 1; j <= d; j++)
        {
            i2j = nmod_mul(i2j, i2, ctx);
            M[j - 1 + d*(i - 1)] = i2j;
/*
flint_printf("M[%wd, %wd]: %wu\n", i - 1, j - 1, i2j);
*/
        }
    }
}

void coeffs2ptvals(
    mp_limb_t * Aptvals,
    slong d,
    const mp_limb_t * Acoeffs,
    slong Alength,
    const mp_limb_t * M,
    const nmod_t ctx)
{
    slong i, j;
    mp_limb_t vp0, vp1, vp2, vm0, vm1, vm2, pp1, pp0;

    if (Alength <= 1)
    {
        vp0 = Alength == 1 ? Acoeffs[0] : 0;
        for (i = 0; i <= 2*d; i++)
            Aptvals[i] = vp0;
        return;
    }

    Aptvals[0] = Acoeffs[0];

    for (i = 1; i <= d; i++)
    {
        vp2 = 0; vp1 = 0; vp0 = Acoeffs[0];
        vm2 = 0; vm1 = 0; vm0 = Acoeffs[1];
        for (j = 1; j < Alength/2; j++)
        {
            umul_ppmm(pp1, pp0, Acoeffs[2*j + 0], M[j - 1]);
            add_sssaaaaaa(vp2, vp1, vp0, vp2, vp1, vp0, 0, pp1, pp0);
            umul_ppmm(pp1, pp0, Acoeffs[2*j + 1], M[j - 1]);
            add_sssaaaaaa(vm2, vm1, vm0, vm2, vm1, vm0, 0, pp1, pp0);
        }
        NMOD_RED3(vm0, vm2, vm1, vm0, ctx);
        if ((Alength % 2) != 0)
        {
            umul_ppmm(pp1, pp0, Acoeffs[2*j + 0], M[j - 1]);
            add_sssaaaaaa(vp2, vp1, vp0, vp2, vp1, vp0, 0, pp1, pp0);
        }
        NMOD_RED3(vp0, vp2, vp1, vp0, ctx);
        vm0 = nmod_mul(i, vm0, ctx);
        Aptvals[2*i - 1] = nmod_add(vp0, vm0, ctx);
        Aptvals[2*i - 0] = nmod_sub(vp0, vm0, ctx);
        M += d;
    }
}

/* L is a (2*d + 1) x (d + 1) array */
void computeL(
    mp_limb_t * L,
    slong d,
    const nmod_t ctx)
{
    slong i, j;
    nmod_poly_t f, g;
    mp_limb_t c;

    nmod_poly_init_mod(f, ctx);
    nmod_poly_init_mod(g, ctx);

    nmod_poly_set_coeff_ui(f, 1, 1);
    nmod_poly_set_coeff_ui(g, 2, 1);
    for (i = 1; i <= d; i++)
    {
        nmod_poly_set_coeff_ui(g, 0, nmod_neg(nmod_mul(i, i, ctx), ctx));
        nmod_poly_mul(f, f, g);
    }
    for (i = 0; i <= d; i++)
    {
        nmod_poly_zero(g);
        nmod_poly_set_coeff_ui(g, 1, 1);
        nmod_poly_set_coeff_ui(g, 0, nmod_neg(i, ctx));
        nmod_poly_div(g, f, g);
        c = nmod_poly_evaluate_nmod(g, i);
        c = nmod_inv(c, ctx);
        nmod_poly_scalar_mul_nmod(g, g, c);
        for (j = 0; j <= 2*d; j++)
        {
            L[j*(d + 1) + i] = nmod_poly_get_coeff_ui(g, j);
/*
flint_printf("L[%wd, %wd]: %wu\n", j, i, nmod_poly_get_coeff_ui(g, j));
*/
        }
    }

    nmod_poly_clear(f);
    nmod_poly_clear(g);
}

void ptvals2coeffs(
    mp_limb_t * Acoeffs,
    slong Alength,
    const mp_limb_t * Aptvals,
    slong d,
    const mp_limb_t * L,
    mp_limb_t * t,
    const nmod_t ctx)
{
    slong i, j;
    mp_limb_t vp2, vp1, vp0, vm2, vm1, vm0, pp1, pp0;

    t[0] = Aptvals[0];
    for (i = 1; i <= d; i++)
    {
        vp0 = Aptvals[2*i - 1];
        vm0 = Aptvals[2*i - 0];
        t[2*i - 1] = nmod_add(vp0, vm0, ctx);
        t[2*i - 0] = nmod_sub(vp0, vm0, ctx);
    }

    for (i = 0; i + 2 <= Alength; i += 2)
    {
        vp2 = 0; umul_ppmm(vp1, vp0, L[0], t[0]);
        vm2 = 0; umul_ppmm(vm1, vm0, L[d + 1], t[0]);
        for (j = 1; j <= d; j++)
        {
            umul_ppmm(pp1, pp0, t[2*j - 1], L[j]);
            add_sssaaaaaa(vp2, vp1, vp0, vp2, vp1, vp0, 0, pp1, pp0);
            umul_ppmm(pp1, pp0, t[2*j - 0], L[j + d + 1]);
            add_sssaaaaaa(vm2, vm1, vm0, vm2, vm1, vm0, 0, pp1, pp0);
        }
        NMOD_RED3(Acoeffs[i + 0], vp2, vp1, vp0, ctx);
        NMOD_RED3(Acoeffs[i + 1], vm2, vm1, vm0, ctx);
        L += 2*(d + 1);
    }

    if (i < Alength)
    {
        vp2 = 0; umul_ppmm(vp1, vp0, L[0], t[0]);
        for (j = 1; j <= d; j++)
        {
            umul_ppmm(pp1, pp0, t[2*j - 1], L[j]);
            add_sssaaaaaa(vp2, vp1, vp0, vp2, vp1, vp0, 0, pp1, pp0);
        }
        NMOD_RED3(Acoeffs[i], vp2, vp1, vp0, ctx);
    }
}

void n_poly_mod_to_ptvals(
    n_poly_t v,
    slong d,
    const n_poly_t a,
    mp_limb_t * M,
    nmod_t mod)
{
    n_poly_fit_length(v, 2*d + 1);
    v->length = 2*d + 1;
    coeffs2ptvals(v->coeffs, d, a->coeffs, a->length, M, mod);
}

void nmod_ptvals_to_n_poly(
    n_poly_t a,
    const mp_limb_t * v,
    slong d,
    mp_limb_t * L,
    nmod_t mod)
{
    mp_limb_t * t = flint_malloc((2*d + 1)*sizeof(mp_limb_t));
    n_poly_fit_length(a, 2*d + 1);
    ptvals2coeffs(a->coeffs, 2*d + 1, v, d, L, t, mod);
    a->length = 2*d + 1;
    _n_poly_normalise(a);
    flint_free(t);
}


void _nmod_vec_addmul(
    mp_limb_t * a,
    const mp_limb_t * b,
    const mp_limb_t * c,
    slong length,
    nmod_t ctx)
{
    slong i;
    for (i = 0; i < length; i++)
        a[i] = nmod_add(a[i], nmod_mul(b[i], c[i], ctx), ctx);
}

void n_bpoly_mod_taylor_shift_outer(n_bpoly_t A, mp_limb_t c, nmod_t mod)
{
    slong n, i, j;
    n_poly_t t;

    FLINT_ASSERT(c < mod.n);

    if (c == 0)
        return;

    n_poly_init(t);
    n = A->length;

    for (i = n - 2; i >= 0; i--)
    {
        for (j = i; j < n - 1; j++)
        {
            n_poly_mod_scalar_mul_nmod(t, A->coeffs + j + 1, c, mod);
            n_poly_mod_add(A->coeffs + j, A->coeffs + j, t, mod);
        }
    }

    n_poly_clear(t);
}

slong n_bpoly_degree0(const n_bpoly_t A)
{
    return A->length - 1;
}

slong n_bpoly_degree1(const n_bpoly_t A)
{
    slong i, len = 0;
    for (i = 0; i < A->length; i++)
        len = FLINT_MAX(len, A->coeffs[i].length);
    return len - 1;    
}

slong n_bpoly_degree_inner(const n_bpoly_t A)
{
    slong i, len = 0;
    for (i = 0; i < A->length; i++)
        len = FLINT_MAX(len, A->coeffs[i].length);
    return len - 1;
}

void n_bpoly_mod_taylor_step(
    n_poly_t c,
    n_bpoly_t A,
    mp_limb_t alpha,
    nmod_t mod)
{
    slong i;

    FLINT_ASSERT(alpha == 0);

    if (A->length < 1)
    {
        n_poly_zero(c);
        return;
    }

    n_poly_swap(c, A->coeffs + 0);
    for (i = 0; i < A->length; i++)
        n_poly_swap(A->coeffs + i, A->coeffs + i + 1);

    A->length--;
}


void n_bpoly_mod_to_ptvals(
    n_bpoly_t A,
    slong d,
    const n_bpoly_t B,
    mp_limb_t * M,
    nmod_t mod)
{
    slong i;
    n_bpoly_fit_length(A, B->length);
    A->length = B->length;
    for (i = 0; i < B->length; i++)
        n_poly_mod_to_ptvals(A->coeffs + i, d, B->coeffs + i, M, mod);
}



int dense_bivar_disolve2(
    n_bpoly_t Delta1,
    n_bpoly_t Delta2,
    n_bpoly_t A,
    n_bpoly_t B1,
    n_bpoly_t B2,
    nmod_t mod,
    slong degy,
    slong Delta1degy,
    slong Delta2degy)
{
    slong d = 10, l, i, j;
    mp_limb_t alpha = 0;
    mp_limb_t * M, * L, * v;
    n_poly_t g, q, e, b2i, b1i, delta;
    n_bpoly_t beta1ptvals, beta2ptvals;
    n_bpoly_t delta1ptvals, delta2ptvals;

printf("dense_bivar_disolve2 called\n");

    printf("A: "); n_bpoly_print_pretty(A, "y", "x"); printf("\n");
    printf("B1: "); n_bpoly_print_pretty(B1, "y", "x"); printf("\n");
    printf("B2: "); n_bpoly_print_pretty(B2, "y", "x"); printf("\n");

    n_poly_init(delta);
    n_poly_init(g);
    n_poly_init(q);
    n_poly_init(e);
    n_poly_init(b1i);
    n_poly_init(b2i);
    n_bpoly_init(beta1ptvals);
    n_bpoly_init(beta2ptvals);
    n_bpoly_init(delta1ptvals);
    n_bpoly_init(delta2ptvals);
    v = (mp_limb_t *) flint_malloc((2*d + 1)*sizeof(mp_limb_t));
    M = (mp_limb_t *) flint_malloc(d*d*sizeof(mp_limb_t));
    L = (mp_limb_t *) flint_malloc((2*d + 1)*(d + 1)*sizeof(mp_limb_t));
    computeM(M, d, mod);
    computeL(L, d, mod);

try_alpha:

    n_bpoly_mod_taylor_shift_outer(B1, alpha, mod);
    n_bpoly_mod_taylor_shift_outer(B2, alpha, mod);

    FLINT_ASSERT(B1->length > 0);
    FLINT_ASSERT(B2->length > 0);

printf("B1[0]: "); n_poly_print_pretty(B1->coeffs + 0, "x"); printf("\n");
printf("B2[0]: "); n_poly_print_pretty(B2->coeffs + 0, "x"); printf("\n");

    n_poly_mod_xgcd(g, b2i, b1i, B1->coeffs + 0, B2->coeffs + 0, mod);
    if (!n_poly_is_one(g))
    {
        FLINT_ASSERT(0);
        alpha++;
        goto try_alpha;
    }

    n_bpoly_mod_to_ptvals(beta2ptvals, d, B1, M, mod);
    n_bpoly_mod_to_ptvals(beta1ptvals, d, B2, M, mod);

printf("beta1ptvals: "); n_bpoly_print_pretty(beta1ptvals, "y", "_"); printf("\n");
printf("beta2ptvals: "); n_bpoly_print_pretty(beta2ptvals, "y", "_"); printf("\n");

    n_bpoly_fit_length(delta1ptvals, degy + 1);
    n_bpoly_fit_length(delta2ptvals, degy + 1);
    delta1ptvals->length = 0;
    delta2ptvals->length = 0;

    n_bpoly_fit_length(Delta1, degy + 1);
    n_bpoly_fit_length(Delta2, degy + 1);

    for (l = 0; l <= degy; l++)
    {
flint_printf("------------ l: %wd ------------\n", l);
printf("A: "); n_bpoly_print_pretty(A, "y", "x"); printf("\n");

        _nmod_vec_zero(v, 2*d + 1);
        for (i = 0; i < l; i++)
        {
            j = l - i;
            if (i < delta1ptvals->length && j < beta1ptvals->length)
                _nmod_vec_addmul(v, delta1ptvals->coeffs[i].coeffs, beta1ptvals->coeffs[j].coeffs, 2*d + 1, mod);
            if (i < delta2ptvals->length && j < beta2ptvals->length)
                _nmod_vec_addmul(v, delta2ptvals->coeffs[i].coeffs, beta2ptvals->coeffs[j].coeffs, 2*d + 1, mod);
        }

        nmod_ptvals_to_n_poly(delta, v, d, L, mod);

printf("delta: "); n_poly_print_pretty(delta, "x"); printf("\n");

        n_bpoly_mod_taylor_step(e, A, alpha, mod);

flint_printf("A[%wd]: ", l); n_poly_print_pretty(e, "x"); printf("\n");


        n_poly_mod_sub(e, e, delta, mod);

printf("solving\n");
printf("e: "); n_poly_print_pretty(e, "x"); printf("\n");

printf("B1[0]: "); n_poly_print_pretty(B1->coeffs + 0, "x"); printf("\n");
printf("B2[0]: "); n_poly_print_pretty(B2->coeffs + 0, "x"); printf("\n");

        n_poly_mod_mul(g, e, b1i, mod);
        n_poly_mod_divrem(q, Delta1->coeffs + l, g, B1->coeffs + 0, mod);
        n_poly_mod_mul(g, e, b2i, mod);
        n_poly_mod_divrem(q, Delta2->coeffs + l, g, B2->coeffs + 0, mod);

flint_printf("Delta1[%wd]: ", l); n_poly_print_pretty(Delta1->coeffs + l, "x"); printf("\n");
flint_printf("Delta2[%wd]: ", l); n_poly_print_pretty(Delta2->coeffs + l, "x"); printf("\n");

        n_poly_mod_to_ptvals(delta1ptvals->coeffs + l, d, Delta1->coeffs + l, M, mod);
        n_poly_mod_to_ptvals(delta2ptvals->coeffs + l, d, Delta2->coeffs + l, M, mod);

        if (!n_poly_is_zero(Delta1->coeffs + l))
        {
            delta1ptvals->length = Delta1->length = l + 1;
            if (l > Delta1degy)
            {
                FLINT_ASSERT(0);
                return 0;
            }
        }

        if (!n_poly_is_zero(Delta2->coeffs + l))
        {
            delta2ptvals->length = Delta2->length = l + 1;
            if (l > Delta2degy)
            {
                FLINT_ASSERT(0);
                return 0;
            }
        }
    }

    n_bpoly_mod_taylor_shift_outer(Delta1, nmod_neg(alpha, mod), mod);
    n_bpoly_mod_taylor_shift_outer(Delta2, nmod_neg(alpha, mod), mod);

    printf("Delta1: "); n_bpoly_print_pretty(Delta1, "y", "x"); printf("\n");
    printf("Delta2: "); n_bpoly_print_pretty(Delta2, "y", "x"); printf("\n");

    n_poly_clear(delta);
    n_poly_clear(g);
    n_poly_clear(q);
    n_poly_clear(e);
    n_poly_clear(b1i);
    n_poly_clear(b2i);
    n_bpoly_clear(beta1ptvals);
    n_bpoly_clear(beta2ptvals);
    n_bpoly_clear(delta1ptvals);
    n_bpoly_clear(delta2ptvals);
    flint_free(v);
    flint_free(M);
    flint_free(L);

    return 1;
}




/*
    solve A/(B1*B2) = C1/B1 + C2/B2 in Fp[y][x]

    return:
        1: solution found with deg_y(Ci) <= ldegCibound
        0: no solution exists with deg_y(Ci) <= ldegCibound
       -1: could not find enough evaluation points where the Bi are pariwise prime
       -2: found no evaluation points where the Bi are pariwise prime
*/
int nmod_mpolyn_pfrac2_bivar(
    nmod_mpolyn_t C1, slong ldegC1bound,
    nmod_mpolyn_t C2, slong ldegC2bound,
    nmod_mpolyn_t A,
    nmod_mpolyn_t B1,
    nmod_mpolyn_t B2,
    const nmod_mpoly_ctx_t ctx,
    nmod_poly_stack_t Sp)
{
    int success;
    slong ldegA, ldegB1, ldegB2, ldegC1, ldegC2;
    slong bad_prime_count, bound;
    mp_limb_t alpha, temp;
    nmod_poly_struct * Aevalp, * B1evalp, * B2evalp, * C1evalp, * C2evalp;
    nmod_poly_struct * Aevalm, * B1evalm, * B2evalm, * C1evalm, * C2evalm;
    nmod_poly_struct * modulus, * modulus2, * alphapow, * t1, * t2;
    nmod_mpolyn_struct * T;
    slong N, off, shift;
    flint_bitcnt_t bits = A->bits;
    const char * vars[] = {"x", "y"};

    ldegA = nmod_mpolyn_lastdeg(A, ctx);
    ldegB1 = nmod_mpolyn_lastdeg(B1, ctx);
    ldegB2 = nmod_mpolyn_lastdeg(B2, ctx);

printf("nmod_mpolyn_pfrac2_bivar called\n");
flint_printf(" A(deg %wd): ", ldegA); nmod_mpolyn_print_pretty(A, vars, ctx); printf("\n");
flint_printf("B1(deg %wd): ", ldegB1); nmod_mpolyn_print_pretty(B1, vars, ctx); printf("\n");
flint_printf("B2(deg %wd): ", ldegB2); nmod_mpolyn_print_pretty(B2, vars, ctx); printf("\n");

    FLINT_ASSERT(Sp->ctx->ffinfo->mod.n == ctx->ffinfo->mod.n);
    FLINT_ASSERT(Sp->ctx->minfo->nvars == ctx->minfo->nvars);
    FLINT_ASSERT(Sp->bits == bits);
    FLINT_ASSERT(A->bits == bits);
    FLINT_ASSERT(B1->bits == bits);
    FLINT_ASSERT(B2->bits == bits);
    FLINT_ASSERT(C1->bits == bits);
    FLINT_ASSERT(C2->bits == bits);

    N = mpoly_words_per_exp_sp(bits, ctx->minfo);
    mpoly_gen_offset_shift_sp(&off, &shift, 0, bits, ctx->minfo);

    nmod_poly_stack_fit_request_poly(Sp, 15);
    Aevalp      = nmod_poly_stack_take_top_poly(Sp);
    B1evalp     = nmod_poly_stack_take_top_poly(Sp);
    B2evalp     = nmod_poly_stack_take_top_poly(Sp);
    C1evalp     = nmod_poly_stack_take_top_poly(Sp);
    C2evalp     = nmod_poly_stack_take_top_poly(Sp);
    Aevalm      = nmod_poly_stack_take_top_poly(Sp);
    B1evalm     = nmod_poly_stack_take_top_poly(Sp);
    B2evalm     = nmod_poly_stack_take_top_poly(Sp);
    C1evalm     = nmod_poly_stack_take_top_poly(Sp);
    C2evalm     = nmod_poly_stack_take_top_poly(Sp);
    modulus     = nmod_poly_stack_take_top_poly(Sp);
    modulus2    = nmod_poly_stack_take_top_poly(Sp);
    alphapow    = nmod_poly_stack_take_top_poly(Sp);
    t1          = nmod_poly_stack_take_top_poly(Sp);
    t2          = nmod_poly_stack_take_top_poly(Sp);

    nmod_poly_stack_fit_request_mpolyn(Sp, 1);
    T           = nmod_poly_stack_take_top_mpolyn(Sp);

    bound = ldegA;
    bound = FLINT_MAX(bound, ldegC1bound + ldegB2);
    bound = FLINT_MAX(bound, ldegB1 + ldegC2bound);
    bound += 1;

    bad_prime_count = 0;

    nmod_poly_fit_length(alphapow, FLINT_MAX(WORD(3), bound + 1));
    nmod_poly_one(modulus);

    if ((ctx->ffinfo->mod.n & UWORD(1)) == UWORD(0))
    {
        success = -1;
        goto cleanup;
    }

    alpha = (ctx->ffinfo->mod.n - UWORD(1))/UWORD(2);

    goto choose_prime;

bad_prime:

    if (bad_prime_count > bound)
    {
        success = nmod_poly_degree(modulus) > 0 ? -1 : -2;
        goto cleanup;
    }

    bad_prime_count++;

choose_prime: /* primes are v - alpha, v + alpha */

    if (alpha < 2)
    {
        success = -1;
        goto cleanup;
    }

    alpha--;

    FLINT_ASSERT(0 < alpha && alpha <= ctx->ffinfo->mod.n/2);
    FLINT_ASSERT(alphapow->alloc >= 2);
    alphapow->length = 2;
    alphapow->coeffs[0] = 1;
    alphapow->coeffs[1] = alpha;

    nmod_mpolyn_interp_reduce_2sm_poly(Aevalp, Aevalm, A, alphapow, ctx);
    nmod_mpolyn_interp_reduce_2sm_poly(B1evalp, B1evalm, B1, alphapow, ctx);
    nmod_mpolyn_interp_reduce_2sm_poly(B2evalp, B2evalm, B2, alphapow, ctx);

    /* make sure evaluation point did not drop the degree of a Bi */
    if (nmod_poly_degree(B1evalp) < ((B1->exps + N*0)[off]>>shift) ||
        nmod_poly_degree(B1evalm) < ((B1->exps + N*0)[off]>>shift) ||
        nmod_poly_degree(B2evalp) < ((B2->exps + N*0)[off]>>shift) ||
        nmod_poly_degree(B2evalm) < ((B2->exps + N*0)[off]>>shift))
    {
        goto choose_prime;
    }

    /* image pfrac's */
    if (!nmod_poly_invmod(t1, B2evalp, B1evalp))
        goto bad_prime;
    nmod_poly_mul(t2, Aevalp, t1);
    nmod_poly_rem(C1evalp, t2, B1evalp);
    nmod_poly_mul(t2, B2evalp, C1evalp);
    nmod_poly_sub(Aevalp, Aevalp, t2);
    nmod_poly_div(C2evalp, Aevalp, B1evalp);

    if (!nmod_poly_invmod(t1, B2evalm, B1evalm))
        goto bad_prime;
    nmod_poly_mul(t2, Aevalm, t1);
    nmod_poly_rem(C1evalm, t2, B1evalm);
    nmod_poly_mul(t2, B2evalm, C1evalm);
    nmod_poly_sub(Aevalm, Aevalm, t2);
    nmod_poly_div(C2evalm, Aevalm, B1evalm);

    /* update interpolants */
    if (nmod_poly_degree(modulus) > 0)
    {
        temp = nmod_poly_evaluate_nmod(modulus, alpha);
        FLINT_ASSERT(temp == nmod_poly_evaluate_nmod(modulus, ctx->ffinfo->mod.n - alpha));
        temp = nmod_mul(temp, alpha, ctx->ffinfo->mod);
        temp = nmod_add(temp, temp, ctx->ffinfo->mod);
        temp = n_invmod(temp, ctx->ffinfo->mod.n);
        nmod_poly_scalar_mul_nmod(modulus, modulus, temp);
        nmod_mpolyn_interp_crt_2sm_poly(&ldegC1, C1, T, C1evalp, C1evalm, modulus, alphapow, ctx);
        nmod_mpolyn_interp_crt_2sm_poly(&ldegC2, C2, T, C2evalp, C2evalm, modulus, alphapow, ctx);
    }
    else
    {
        nmod_mpolyn_interp_lift_2sm_poly(&ldegC1, C1, C1evalp, C1evalm, alpha, ctx);
        nmod_mpolyn_interp_lift_2sm_poly(&ldegC2, C2, C2evalp, C2evalm, alpha, ctx);
    }
    temp = nmod_mul(alpha, alpha, ctx->ffinfo->mod);
    nmod_poly_scalar_mul_nmod(modulus2, modulus, temp);
    nmod_poly_shift_left(modulus, modulus, 2);
    nmod_poly_sub(modulus, modulus, modulus2);

    if (ldegC1 > ldegC1bound || ldegC2 > ldegC2bound)
    {
        success = 0;
        goto cleanup;
    }

    if (nmod_poly_degree(modulus) < bound)
    {
        goto choose_prime;
    }

    success = 1;

cleanup:

    nmod_poly_stack_give_back_poly(Sp, 15);
    nmod_poly_stack_give_back_mpolyn(Sp, 1);

flint_printf("nmod_mpolyn_pfrac2_bivar returning success = %d\n", success);
printf("C1: "); nmod_mpolyn_print_pretty(C1, vars, ctx); printf("\n");
printf("C2: "); nmod_mpolyn_print_pretty(C2, vars, ctx); printf("\n");

    return success;
}



/*
    Try to solve A/(B[0]*...*B[r-1]) = C[0]/B[0] + ... + C[r-1]/B[r-1]
    for the C[i] in Fp[y][x].
    return:
        1: solution found with deg_y(Ci) <= ldegCibound
        0: no solution exists with deg_y(Ci) <= ldegCibound
       -1: could not find enough evaluation points where the Bi are pariwise prime
       -2: found no evaluation points where the Bi are pariwise prime
*/
int nmod_mpolyn_pfrac_bivar(
    slong r,
    nmod_mpolyn_struct * C,
    slong * ldegCbound,
    nmod_mpolyn_t A,
    nmod_mpolyn_struct * B,
    const nmod_mpoly_ctx_t ctx,
    nmod_poly_stack_t Sp)
{
    int success;
    slong ldegA, ldegBtotal, * ldegB, * ldegC;
    slong i, j, bad_prime_count, bound;
    mp_limb_t alpha, temp;
    nmod_poly_struct * Aevalp, ** Bevalp, ** Cevalp;
    nmod_poly_struct * Aevalm, ** Bevalm, ** Cevalm;
    nmod_poly_struct * modulus, * modulus2, * alphapow, * t1, * t2;
    nmod_mpolyn_struct * T;
    slong N, off, shift;
    flint_bitcnt_t bits = A->bits;
    const char * vars[] = {"x", "y"};
    TMP_INIT;

    TMP_START;

    ldegB = (slong *) TMP_ALLOC(2*r*sizeof(slong));
    ldegC = ldegB + r;

    Bevalp = (nmod_poly_struct **) TMP_ALLOC(4*r*sizeof(nmod_poly_struct *));
    Bevalm = Bevalp + r;
    Cevalp = Bevalm + r;
    Cevalm = Cevalp + r;

    ldegA = nmod_mpolyn_lastdeg(A, ctx);
    for (i = 0; i < r; i++)
        ldegB[i] = nmod_mpolyn_lastdeg(B + i, ctx);

printf("nmod_mpolyn_pfrac_bivar called\n");
flint_printf(" A(deg %wd): ", ldegA); nmod_mpolyn_print_pretty(A, vars, ctx); printf("\n");
for (i = 0; i < r; i++)
{
    flint_printf("B[%wd](deg %wd): ", i, ldegB[i]); nmod_mpolyn_print_pretty(B + i, vars, ctx); printf("\n");
}

    FLINT_ASSERT(r > 2);
    FLINT_ASSERT(Sp->ctx->ffinfo->mod.n == ctx->ffinfo->mod.n);
    FLINT_ASSERT(Sp->ctx->minfo->nvars == ctx->minfo->nvars);
    FLINT_ASSERT(Sp->bits == bits);
    FLINT_ASSERT(A->bits == bits);
    for (i = 0; i < r; i++)
    {
        FLINT_ASSERT(B[i].bits == bits);
        FLINT_ASSERT(C[i].bits == bits);
    }
    N = mpoly_words_per_exp_sp(bits, ctx->minfo);
    mpoly_gen_offset_shift_sp(&off, &shift, 0, bits, ctx->minfo);

    nmod_poly_stack_fit_request_poly(Sp, 7 + 4*r);
    Aevalp      = nmod_poly_stack_take_top_poly(Sp);
    Aevalm      = nmod_poly_stack_take_top_poly(Sp);
    modulus     = nmod_poly_stack_take_top_poly(Sp);
    modulus2    = nmod_poly_stack_take_top_poly(Sp);
    alphapow    = nmod_poly_stack_take_top_poly(Sp);
    t1          = nmod_poly_stack_take_top_poly(Sp);
    t2          = nmod_poly_stack_take_top_poly(Sp);
    for (i = 0; i < r; i++)
    {
        Bevalp[i] = nmod_poly_stack_take_top_poly(Sp);
        Bevalm[i] = nmod_poly_stack_take_top_poly(Sp);
        Cevalp[i] = nmod_poly_stack_take_top_poly(Sp);
        Cevalm[i] = nmod_poly_stack_take_top_poly(Sp);
    }

    nmod_poly_stack_fit_request_mpolyn(Sp, 1);
    T           = nmod_poly_stack_take_top_mpolyn(Sp);

    ldegBtotal = ldegB[0];
    for (i = 1; i < r; i++)
        ldegBtotal += ldegB[i];

    bound = ldegA;
    for (i = 0; i < r; i++)
        bound = FLINT_MAX(bound, ldegCbound[i] + ldegBtotal - ldegB[i]);
    bound += 1;

    bad_prime_count = 0;

    nmod_poly_fit_length(alphapow, FLINT_MAX(WORD(3), bound + 1));
    nmod_poly_one(modulus);

    if ((ctx->ffinfo->mod.n & UWORD(1)) == UWORD(0))
    {
        success = -1;
        goto cleanup;
    }

    alpha = (ctx->ffinfo->mod.n - UWORD(1))/UWORD(2);

    goto choose_prime;

bad_prime:

    if (bad_prime_count > 2*bound)
    {
        success = nmod_poly_degree(modulus) > 0 ? -1 : -2;
        goto cleanup;
    }

    bad_prime_count++;

choose_prime: /* primes are v - alpha, v + alpha */

    if (alpha < 2)
    {
        success = -1;
        goto cleanup;
    }

    alpha--;

    FLINT_ASSERT(0 < alpha && alpha <= ctx->ffinfo->mod.n/2);
    FLINT_ASSERT(alphapow->alloc >= 2);
    alphapow->length = 2;
    alphapow->coeffs[0] = 1;
    alphapow->coeffs[1] = alpha;

    /* generate images */
    nmod_mpolyn_interp_reduce_2sm_poly(Aevalp, Aevalm, A, alphapow, ctx);
    for (i = 0; i < r; i++)
    {
        nmod_mpolyn_interp_reduce_2sm_poly(Bevalp[i], Bevalm[i], B + i, alphapow, ctx);
        /* make sure evaluation point did not drop the degree of a Bi */
        if (nmod_poly_degree(Bevalp[i]) < ((B[i].exps + N*0)[off]>>shift) ||
            nmod_poly_degree(Bevalm[i]) < ((B[i].exps + N*0)[off]>>shift))
        {
            goto choose_prime;
        }
    }

    /* image pfrac's */
    for (i = 0; i < r; i++)
    {
        nmod_poly_one(t2);
        for (j = 0; j < r; j++)
        {
            if (j == i)
                continue;
            nmod_poly_mul(t1, t2, Bevalp[j]);
            nmod_poly_swap(t1, t2);
        }
        if (!nmod_poly_invmod(t1, t2, Bevalp[i]))
            goto bad_prime;
        nmod_poly_mul(t2, Aevalp, t1);
        nmod_poly_rem(Cevalp[i], t2, Bevalp[i]);

        nmod_poly_one(t2);
        for (j = 0; j < r; j++)
        {
            if (j == i)
                continue;
            nmod_poly_mul(t1, t2, Bevalm[j]);
            nmod_poly_swap(t1, t2);
        }
        if (!nmod_poly_invmod(t1, t2, Bevalm[i]))
            goto bad_prime;
        nmod_poly_mul(t2, Aevalm, t1);
        nmod_poly_rem(Cevalm[i], t2, Bevalm[i]);
    }

    /* update interpolants */
    if (nmod_poly_degree(modulus) > 0)
    {
        temp = nmod_poly_evaluate_nmod(modulus, alpha);
        FLINT_ASSERT(temp == nmod_poly_evaluate_nmod(modulus, ctx->ffinfo->mod.n - alpha));
        temp = nmod_mul(temp, alpha, ctx->ffinfo->mod);
        temp = nmod_add(temp, temp, ctx->ffinfo->mod);
        temp = n_invmod(temp, ctx->ffinfo->mod.n);
        nmod_poly_scalar_mul_nmod(modulus, modulus, temp);
        for (i = 0; i < r; i++)
        {
            nmod_mpolyn_interp_crt_2sm_poly(ldegC + i, C + i, T, Cevalp[i],
                                           Cevalm[i], modulus, alphapow, ctx);
        }
    }
    else
    {
        for (i = 0; i < r; i++)
        {
            nmod_mpolyn_interp_lift_2sm_poly(ldegC + i, C + i, Cevalp[i],
                                                       Cevalm[i], alpha, ctx);
        }
    }
    temp = nmod_mul(alpha, alpha, ctx->ffinfo->mod);
    nmod_poly_scalar_mul_nmod(modulus2, modulus, temp);
    nmod_poly_shift_left(modulus, modulus, 2);
    nmod_poly_sub(modulus, modulus, modulus2);

    for (i = 0; i < r; i++)
    {
        if (ldegC[i] > ldegCbound[i])
        {
            success = 0;
            goto cleanup;
        }
    }

    if (nmod_poly_degree(modulus) < bound)
    {
        goto choose_prime;
    }

    success = 1;

cleanup:

    nmod_poly_stack_give_back_poly(Sp, 7 + 4*r);
    nmod_poly_stack_give_back_mpolyn(Sp, 1);

    TMP_END;

flint_printf("nmod_mpolyn_pfrac2_bivar returning success = %d\n", success);
for (i = 0; i < r; i++)
{
flint_printf("C[%wd]: ", i); nmod_mpolyn_print_pretty(C + i, vars, ctx); printf("\n");
}

    return success;
}


/*
    Try to solve A/(B[0]*...*B[r-1]) = C[0]/B[0] + ... + C[r-1]/B[r-1]
    in Fp[x1,...,xn][X,Y]
*/
#if 0
int nmod_mpoly_pfrac(
    slong r,
    nmod_mpolyu_struct * C,
    nmod_zip_polyu_struct * Cassumed,
    const nmod_mpolyu_t B,
    const nmod_mpolyu_t A,
    const nmod_mpoly_t ctx)
{
    slong zip_evals;
    slong * Cdegbound;
    const nmod_mpoly_ctx_t ctx2;
    nmod_poly_stack_t Sp;

    nmod_mpoly_ctx_init(ctx2, 2, ORD_LEX, ctx->ffinfo->mod.n);

    Z->pointcount = 0;

    zip_evals = 0;
    for (i = 0; i < r; i++)
    {
        for (j = 0; j < Cassumed[i].length; j++)
            zip_evals = FLINT_MAX(zip_evals, Cassumed[i].coeffs[j].length);
    }
    zip_evals += 1;

    /* set evaluation of monomials */
    nmod_mpolyu_set_skel_nmod(Ainc, A, alpha, ctx);
    for (i = 0; i < r; i++)
        nmod_mpolyu_set_skel_nmod(Binc + i, B + i, alpha, ctx);

    /* set reduction of coeffs */
    nmod_mpolyu_red_skel_nmod(Ared, A, ctx_sp->ffinfo);
    for (i = 0; i < r; i++)
        nmod_mpolyu_red_skel_nmod(Bred + i, B + i, ctx_sp->ffinfo);

    /* copy evaluation of monomials */
    nmod_mpolyu_copy_skel(Acur, Ainc);
    for (i = 0; i < r; i++)
        nmod_mpolyu_copy_skel(Bcur + i, Binc + i);

next_zip_image:

    nmod_mpolyuu_use_skel_mul_nmod(Aeval, A, Ared, Acur, Ainc, ctx2);
    for (i = 0; i < r; i++)
        nmod_mpolyuu_use_skel_mul_nmod(Beval + i, B + i, Bred + i, Bcur + i, Binc + i, ctx2);

    success = nmod_mpolyn_pfrac_bivar(r, Ceval, Cdeg, Aeval, Beval, ctx, Sp);
    FLINT_ASSERT(success == 1);

    for (i = 0; i < r; i++)
    {
        success = nmod_zip_mpolyuu_add_point(Z + i, Ceval);
        FLINT_ASSERT(success == 1);
    }

    if (Z[0].pointcount < zip_evals)
        goto next_zip_image;

    

}
#endif

/*
int n_poly_mod_is_canonical(const n_poly_t A, nmod_t mod)
{
    slong i;

    if (A->length <= 0)
        return A->length == 0;

    for (i = 0; i < A->length; i++)
    {
        if (A->coeffs[i] >= mod.n)
            return 0;
    }

    return 0 != A->coeffs[A->length - 1];
}
*/

int n_bpoly_mod_is_canonical(const n_bpoly_t A, nmod_t mod)
{
    slong i;

    if (A->length <= 0)
        return A->length == 0;

    for (i = 0; i < A->length; i++)
    {
        if (!n_poly_mod_is_canonical(A->coeffs + i, mod))
        {
flint_printf("n_bpoly coeff i = %wd is bad\n", i);
            return 0;
        }
    }

    return !n_poly_is_zero(A->coeffs + A->length - 1);
}


void n_bpoly_mod_mul(n_bpoly_t A, const n_bpoly_t B, const n_bpoly_t C,
                                                                    nmod_t mod)
{
    slong i, j;
    n_poly_struct * t;

    FLINT_ASSERT(A != B);
    FLINT_ASSERT(A != C);

    if (B->length <= 0 || C->length <= 0)
    {
        A->length = 0;
        return;
    }

    n_bpoly_fit_length(A, B->length + C->length);
    for (i = 0; i < B->length + C->length - 1; i++)
        n_poly_zero(A->coeffs + i);

    t = A->coeffs + B->length + C->length - 1;

    for (i = 0; i < B->length; i++)
    {
        for (j = 0; j < C->length; j++)
        {
            _n_poly_mod_mul(t, B->coeffs + i, C->coeffs + j, mod);
            n_poly_mod_add(A->coeffs + i + j, A->coeffs + i + j, t, mod);
        }
    }

    A->length = B->length + C->length - 1;
}

int n_bpoly_equal(const n_bpoly_t A, const n_bpoly_t B)
{
    slong i;

    if (A->length != B->length)
        return 0;

    for (i = 0; i < A->length; i++)
    {
        if (!n_poly_equal(A->coeffs + i, B->coeffs + i))
            return 0;
    }

    return 1;
}

void n_bpoly_normalize(n_bpoly_t A)
{
    while (A->length > 0 && n_poly_is_zero(A->coeffs + A->length - 1))
        A->length--;
}

/*
    input A, B0, B1 with A(y,x) = B0(y,x) * B1(y,x) mod (y-alpha)
    return
       -1: B0(alpha,x) & B1(alpha,x) are not pairwise prime, or
           A(alpha,x) has wrong degree w.r.t x
        0: lift of B0 and B1 to true factors is impossible
        1: successfully lifted B0 and B1 to true factors without changing lc_x
*/
int n_bpoly_mod_hensel_lift2(
    n_bpoly_t A, /* clobbered (shifted by alpha) */
    n_bpoly_t B0,
    n_bpoly_t B1,
    mp_limb_t alpha,
    slong degree_inner, /* required degree in x */
    nmod_t mod)
{
    int success;
    slong i, j;
    n_poly_t c, s, t, u, v;
/*
flint_printf("n_bpoly_hensel_lift2 called: degree_inner = %wd, alpha = %wu\n", degree_inner, alpha);
printf(" A: "); n_bpoly_print_pretty(A, "Y", "X"); printf("\n"); fflush(stdout);
printf("B0: "); n_bpoly_print_pretty(B0, "Y", "X"); printf("\n"); fflush(stdout);
printf("B1: "); n_bpoly_print_pretty(B1, "Y", "X"); printf("\n"); fflush(stdout);
*/
    FLINT_ASSERT(n_bpoly_mod_is_canonical(A, mod));
    FLINT_ASSERT(n_bpoly_mod_is_canonical(B0, mod));
    FLINT_ASSERT(n_bpoly_mod_is_canonical(B1, mod));

    n_poly_init(c);
    n_poly_init(s);
    n_poly_init(t);
    n_poly_init(u);
    n_poly_init(v);

    n_bpoly_mod_taylor_shift_outer(A, alpha, mod);
    n_bpoly_mod_taylor_shift_outer(B0, alpha, mod);
    n_bpoly_mod_taylor_shift_outer(B1, alpha, mod);
/*
printf("shift A: "); n_bpoly_print_pretty(A, "Y", "X"); printf("\n"); fflush(stdout);
printf("shift B0: "); n_bpoly_print_pretty(B0, "Y", "X"); printf("\n"); fflush(stdout);
printf("shift B1: "); n_bpoly_print_pretty(B1, "Y", "X"); printf("\n"); fflush(stdout);
*/

    if (A->length <= 0 || B0->length <= 0 || B1->length <= 0)
    {
        success = -1;
        goto cleanup;
    }

    /* supposed to have A(alpha,x) = B0(alpha,x) * B1(alpha,x) */
    FLINT_ASSERT(n_poly_degree(A->coeffs + 0) == n_poly_degree(B0->coeffs + 0) +
                                                 n_poly_degree(B1->coeffs + 0));

    if (n_poly_degree(A->coeffs + 0) != degree_inner)
    {
        success = -1;
        goto cleanup;
    }

    /* the required degree in x is supposed to be deg_x(A) */
    FLINT_ASSERT(n_bpoly_degree_inner(A) == n_poly_degree(A->coeffs + 0));

    if (!n_poly_mod_invmod(s, B1->coeffs + 0, B0->coeffs + 0, mod))
    {
        success = -2;
        goto cleanup;
    }

    n_bpoly_fit_length(B0, A->length);
    n_bpoly_fit_length(B1, A->length);
    for (j = 1; j < A->length; j++)
    {
        n_poly_set(c, A->coeffs + j);
        for (i = 0; i <= j; i++)
        {
            if (i < B0->length && j - i < B1->length)
            {
                n_poly_mod_mul(t, B0->coeffs + i, B1->coeffs + j - i, mod);
                n_poly_mod_sub(c, c, t, mod);
            }
        }

        if (n_poly_is_zero(c))
            continue;

        n_poly_mod_mul(t, s, c, mod);
        n_poly_mod_rem(u, t, B0->coeffs + 0, mod);
        n_poly_mod_mul(t, u, B1->coeffs + 0, mod);
        n_poly_mod_sub(c, c, t, mod);
        n_poly_mod_div(v, c, B0->coeffs + 0, mod);

        if (j < B0->length)
            n_poly_mod_add(B0->coeffs + j, B0->coeffs + j, u, mod);
        else
            n_poly_set(B0->coeffs + j, u);

        if (j < B1->length)
            n_poly_mod_add(B1->coeffs + j, B1->coeffs + j, v, mod);
        else
            n_poly_set(B1->coeffs + j, v);

        if (!n_poly_is_zero(B0->coeffs + j))
            B0->length = FLINT_MAX(B0->length, j + 1);
        if (!n_poly_is_zero(B1->coeffs + j))
            B1->length = FLINT_MAX(B1->length, j + 1);

        if (B0->length - 1 + B1->length - 1 > A->length - 1)
        {
            success = 0;
            goto cleanup;
        }
    }

    n_bpoly_mod_taylor_shift_outer(B0, nmod_neg(alpha, mod), mod);
    n_bpoly_mod_taylor_shift_outer(B1, nmod_neg(alpha, mod), mod);

    success = 1;

cleanup:
/*
flint_printf("n_bpoly_hensel_lift2 returning %d\n", success);
flint_printf("B0: "); n_bpoly_print_pretty(B0, "Y", "X"); flint_printf("\n");
flint_printf("B1: "); n_bpoly_print_pretty(B1, "Y", "X"); flint_printf("\n");
*/

    if (success)
    {
        n_bpoly_t tp1, tp2;
        n_bpoly_init(tp1);
        n_bpoly_init(tp2);

        n_bpoly_mod_taylor_shift_outer(A, nmod_neg(alpha, mod), mod);
        n_bpoly_mod_mul(tp1, B0, B1, mod);
        FLINT_ASSERT(n_bpoly_equal(tp1, A));

        n_bpoly_clear(tp1);
        n_bpoly_clear(tp2);
    }

    n_poly_clear(c);
    n_poly_clear(s);
    n_poly_clear(t);
    n_poly_clear(u);
    n_poly_clear(v);

    return success;
}
/* r factor version */
int n_bpoly_mod_hensel_lift(
    slong r,
    n_bpoly_t A, /* clobbered (shifted by alpha) */
    n_bpoly_struct * B,
    mp_limb_t alpha,
    slong degree_inner, /* required degree in x */
    nmod_t mod)
{
    int success;
    slong i, j, k, tdeg;
    n_poly_struct * s, * v;
    n_poly_t c, t, u;
    n_bpoly_struct * U;
/*
flint_printf("n_bpoly_hensel_lift(%wd) called: degree_inner = %wd, alpha = %wu\n", r, degree_inner, alpha);
flint_printf(" A: ");
n_bpoly_print_pretty(A, "y", "x");
flint_printf("\n");
for (i = 0; i < r; i++)
{
flint_printf("B[%wd]: ", i);
n_bpoly_print_pretty(B + i, "y", "x");
flint_printf("\n");
}
*/
    FLINT_ASSERT(r > 2);
    FLINT_ASSERT(n_bpoly_mod_is_canonical(A, mod));

    for (i = 0; i < r; i++)
    {
        FLINT_ASSERT(B[i].length > 0);
        FLINT_ASSERT(n_bpoly_mod_is_canonical(B + i, mod));
    }

    U = (n_bpoly_struct *) flint_malloc(r*sizeof(n_bpoly_struct));
    for (i = 0; i < r; i++)
    {
        n_bpoly_init(U + i);
        n_bpoly_fit_length(U + i, A->length);
        for (j = 0; j < A->length; j++)
            n_poly_zero(U[i].coeffs + j);
        U[i].length = A->length;
        n_bpoly_fit_length(B + i, A->length);
    }

    s = (n_poly_struct *) flint_malloc(r*sizeof(n_poly_struct));
    v = (n_poly_struct *) flint_malloc(r*sizeof(n_poly_struct));
    for (i = 0; i < r; i++)
    {
        n_poly_init(s + i);
        n_poly_init(v + i);
    }

    n_poly_init(c);
    n_poly_init(t);
    n_poly_init(u);

    n_bpoly_mod_taylor_shift_outer(A, alpha, mod);
    for (i = 0; i < r; i++)
        n_bpoly_mod_taylor_shift_outer(B + i, alpha, mod);

    /* supposed to have A(alpha,x) = B0(alpha,x) * B1(alpha,x) * ... */
    j = 0;
    for (i = 0; i < r; i++)
        j += n_poly_degree(B[i].coeffs + 0);
    FLINT_ASSERT(j == n_poly_degree(A->coeffs + 0));

    if (n_poly_degree(A->coeffs + 0) != degree_inner)
    {
        success = -1;
        goto cleanup;
    }

    /* the required degree in x is supposed to be deg_x(A) */
    FLINT_ASSERT(n_bpoly_degree_inner(A) == n_poly_degree(A->coeffs + 0));

    for (i = 0; i < r; i++)
    {
        n_poly_one(t);
        for (j = 0; j < r; j++)
        {
            if (j != i)
                n_poly_mod_mul(t, t, B[j].coeffs + 0, mod);
        }
        if (!n_poly_mod_invmod(s + i, t, B[i].coeffs + 0, mod))
        {
            success = -1;
            goto cleanup;
        }
    }

    k = r - 2;
    n_poly_mod_mul(U[k].coeffs + 0, B[k].coeffs + 0, B[k + 1].coeffs + 0, mod);
    for (k--; k > 0; k--)
        n_poly_mod_mul(U[k].coeffs + 0, B[k].coeffs + 0, U[k + 1].coeffs + 0, mod);

    for (j = 1; j < A->length; j++)
    {
        for (k = 0; k < r; k++)
            n_poly_zero(U[k].coeffs + j);

        k = r - 2;
        n_poly_zero(U[k].coeffs + j);
        for (i = 0; i <= j; i++)
        {
            if (i < B[k].length && j - i < B[k + 1].length)
            {
                n_poly_mod_mul(t, B[k].coeffs + i, B[k + 1].coeffs + j - i, mod);
                n_poly_mod_add(U[k].coeffs + j, U[k].coeffs + j, t, mod);
            }
        }
        for (k--; k > 0; k--)
        {
            n_poly_zero(U[k].coeffs + j);
            for (i = 0; i <= j; i++)
            {
                if (i < B[k].length)
                {
                    n_poly_mod_mul(t, B[k].coeffs + i, U[k + 1].coeffs + j - i, mod);
                    n_poly_mod_add(U[k].coeffs + j, U[k].coeffs + j, t, mod);
                }
            }
        }

        n_poly_set(c, A->coeffs + j);

        for (i = 0; i <= j; i++)
        {
            if (i < B[0].length)
            {
                n_poly_mod_mul(t, B[0].coeffs + i, U[1].coeffs + j - i, mod);
                n_poly_mod_sub(c, c, t, mod);
            }
        }

        if (n_poly_is_zero(c))
            continue;

        tdeg = 0;
        for (i = 0; i < r; i++)
        {
            n_poly_mod_mul(t, s + i, c, mod);
            n_poly_mod_rem(v + i, t, B[i].coeffs + 0, mod);
            while (j >= B[i].length)
            {
                n_poly_zero(B[i].coeffs + B[i].length);
                B[i].length++;
            }
            n_poly_mod_add(B[i].coeffs + j, B[i].coeffs + j, v + i, mod);
            n_bpoly_normalize(B + i);
            tdeg += B[i].length - 1;
        }

        if (tdeg >= A->length)
        {
            success = 0;
            goto cleanup;
        }

        k = r - 2;
        n_poly_mod_mul(t, B[k].coeffs + 0, v + k + 1, mod);
        n_poly_mod_mul(u, B[k + 1].coeffs + 0, v + k, mod);
        n_poly_mod_add(t, t, u, mod);
        n_poly_mod_add(U[k].coeffs + j, U[k].coeffs + j, t, mod);
        for (k--; k > 0; k--)
        {
            n_poly_mod_mul(u, B[k].coeffs + 0, t, mod);
            n_poly_mod_mul(t, U[k + 1].coeffs + 0, v + k, mod);
            n_poly_mod_add(t, t, u, mod);
            n_poly_mod_add(U[k].coeffs + j, U[k].coeffs + j, t, mod);
        }
    }

    for (i = 0; i < r; i++)
        n_bpoly_mod_taylor_shift_outer(B + i, nmod_neg(alpha, mod), mod);

    success = 1;

cleanup:
/*
flint_printf("n_bpoly_hensel_lift(%wd) returning %d\n", r, success);
for (i = 0; i < r; i++)
{
flint_printf("B[%wd]: ", i);
n_bpoly_print_pretty(B + i, "y", "x");
flint_printf("\n");
}

FLINT_ASSERT(success);
*/
    if (success)
    {
        n_bpoly_t tp1, tp2;
        n_bpoly_init(tp1);
        n_bpoly_init(tp2);

        n_bpoly_mod_taylor_shift_outer(A, nmod_neg(alpha, mod), mod);
        n_bpoly_mod_mul(tp1, B + 0, B + 1, mod);
        for (i = 2; i < r; i++)
        {
            n_bpoly_mod_mul(tp2, tp1, B + i, mod);
            n_bpoly_swap(tp1, tp2);
        }
        FLINT_ASSERT(n_bpoly_equal(tp1, A));

        n_bpoly_clear(tp1);
        n_bpoly_clear(tp2);
    }

    for (i = 0; i < r; i++)
    {
        n_bpoly_clear(U + i);
        n_poly_clear(s + i);
        n_poly_clear(v + i);
    }
    flint_free(U);
    flint_free(s);
    flint_free(v);

    n_poly_clear(c);
    n_poly_clear(t);
    n_poly_clear(u);

    return success;
}




void nmod_mpolyn_cvtto_bpoly(
    n_bpoly_t A,
    nmod_mpolyn_t B,
    const nmod_mpoly_ctx_t ctx)
{
    slong i, x;

    FLINT_ASSERT(ctx->minfo->nvars == 3);
    FLINT_ASSERT(B->bits == FLINT_BITS/3);
    FLINT_ASSERT(B->length > 0);

    A->length = 0;
    for (i = 0; i < B->length; i++)
    {
        x = B->exps[i] >> (2*(FLINT_BITS/3));
        n_bpoly_fit_length(A, x + 1);
        while (A->length <= x)
        {
            n_poly_zero(A->coeffs + A->length);
            A->length++;
        }
        n_poly_set_nmod_poly(A->coeffs + x, B->coeffs + i);
    }
}

void n_bpoly_cvtto_mpolyn(
    nmod_mpolyn_t A,
    n_bpoly_t B,
    const nmod_mpoly_ctx_t ctx)
{
    slong i;

    FLINT_ASSERT(ctx->minfo->nvars == 3);
    FLINT_ASSERT(A->bits == FLINT_BITS/3);

    nmod_mpolyn_fit_length(A, B->length, ctx);
    A->length = 0;
    for (i = B->length - 1; i >= 0; i--)
    {
        if (!n_poly_is_zero(B->coeffs + i))
        {
            nmod_poly_set_n_poly(A->coeffs + A->length, B->coeffs + i);
            A->exps[A->length] = i << (2*(FLINT_BITS/3));
            A->length++;
        }
    }
}

void n_poly_fill_powers(
    n_poly_t alphapow,
    slong target,
    nmod_t mod)
{
    if (target + 1 > alphapow->length)
    {
        slong k;
        slong oldlength = alphapow->length;
        FLINT_ASSERT(2 <= oldlength);
        n_poly_fit_length(alphapow, target + 1);
        for (k = oldlength; k <= target; k++)
        {
            alphapow->coeffs[k] = nmod_mul(alphapow->coeffs[k - 1],
                                                alphapow->coeffs[1], mod);
        }
        alphapow->length = target + 1;
    }
}



/* multiply A by (x^k + c) */
void n_poly_mod_shift_left_scalar_addmul(n_poly_t A, slong k, mp_limb_t c,
                                                                    nmod_t mod)
{
    mp_limb_t * Acoeffs;
    slong i;
    slong Alen = A->length;

    n_poly_fit_length(A, Alen + k);

    Acoeffs = A->coeffs;

    flint_mpn_copyd(Acoeffs + k, Acoeffs, Alen);
    flint_mpn_zero(Acoeffs, k);

    for (i = 0; i < A->length; i++)
        Acoeffs[i] = nmod_addmul(Acoeffs[i], c, Acoeffs[i + k], mod);

    A->length = Alen + k;
}

/* A = B + C*(d1*x+d0) */
void n_poly_mod_addmul_linear(
    n_poly_t A,
    const n_poly_t B,
    const n_poly_t C,
    mp_limb_t d1, mp_limb_t d0,
    nmod_t mod)
{
    slong i;
    mp_limb_t * Acoeffs, * Bcoeffs, * Ccoeffs;
    slong Blen = B->length;
    slong Clen = C->length;
    slong Alen = FLINT_MAX(B->length, C->length + 1);
/*
flint_printf("n_poly_mod_addmul_linear called\n");
flint_printf("B: "); n_poly_print_pretty(B, "Z"); flint_printf("\n");
flint_printf("C: "); n_poly_print_pretty(C, "Z"); flint_printf("\n");
flint_printf("%wu*Z + %wu\n", d1, d0);
*/
    FLINT_ASSERT(A != B);
    FLINT_ASSERT(A != C);

    n_poly_fit_length(A, Alen);
    Acoeffs = A->coeffs;
    Bcoeffs = B->coeffs;
    Ccoeffs = C->coeffs;

    for (i = 0; i < Alen; i++)
    {
        ulong p1, p0, t0 = 0, t1 = 0, t2 = 0;

        if (i < Blen)
        {
            t0 = Bcoeffs[i];
        }
        if (i < Clen)
        {
            umul_ppmm(p1, p0, Ccoeffs[i], d0);
            add_ssaaaa(t1, t0, t1, t0, p1, p0);
        }
        if (0 < i && i - 1 < Clen)
        {
            umul_ppmm(p1, p0, Ccoeffs[i - 1], d1);
            add_sssaaaaaa(t2, t1, t0, t2, t1, t0, 0, p1, p0);
        }
        NMOD_RED3(Acoeffs[i], t2, t1, t0, mod);
    }

    A->length = Alen;
    _n_poly_normalise(A);
/*
flint_printf("n_poly_mod_addmul_linear returning\n");
flint_printf("A: "); n_poly_print_pretty(A, "Z"); flint_printf("\n");
*/
}

void _n_poly_mod_eval2_pow(
    mp_limb_t * vp,
    mp_limb_t * vm,
    const n_poly_t P, 
    n_poly_t alphapow,
    nmod_t mod)
{
    const mp_limb_t * Pcoeffs = P->coeffs;
    slong Plen = P->length;
    mp_limb_t * alpha_powers = alphapow->coeffs;
    mp_limb_t p1, p0, a0, a1, a2, q1, q0, b0, b1, b2;
    slong k;

    a0 = a1 = a2 = 0;
    b0 = b1 = b2 = 0;

    if (Plen > alphapow->length)
    {
        slong oldlength = alphapow->length;
        FLINT_ASSERT(2 <= oldlength);
        n_poly_fit_length(alphapow, Plen);
        for (k = oldlength; k < Plen; k++)
        {
            alphapow->coeffs[k] = nmod_mul(alphapow->coeffs[k - 1],
                                               alphapow->coeffs[1], mod);
        }
        alphapow->length = Plen;
        alpha_powers = alphapow->coeffs;
    }

    for (k = 0; k + 2 <= Plen; k += 2)
    {
        umul_ppmm(p1, p0, Pcoeffs[k + 0], alpha_powers[k + 0]);
        umul_ppmm(q1, q0, Pcoeffs[k + 1], alpha_powers[k + 1]);
        add_sssaaaaaa(a2, a1, a0, a2, a1, a0, 0, p1, p0);
        add_sssaaaaaa(b2, b1, b0, b2, b1, b0, 0, q1, q0);
    }

    if (k < Plen)
    {
        umul_ppmm(p1, p0, Pcoeffs[k + 0], alpha_powers[k + 0]);
        add_sssaaaaaa(a2, a1, a0, a2, a1, a0, 0, p1, p0);
        k++;
    }

    FLINT_ASSERT(k == Plen);

    NMOD_RED3(p0, a2, a1, a0, mod);
    NMOD_RED3(q0, b2, b1, b0, mod);

    *vp = nmod_add(p0, q0, mod);
    *vm = nmod_sub(p0, q0, mod);
}

void n_bpoly_mod_interp_reduce_2sm_poly(
    n_poly_t Ap,
    n_poly_t Am,
    const n_bpoly_t A,
    n_poly_t alphapow,
    nmod_t mod)
{
    slong i, Alen = A->length;
    const n_poly_struct * Ac = A->coeffs;
    mp_limb_t * Apc, * Amc;

    n_poly_fit_length(Ap, Alen);
    n_poly_fit_length(Am, Alen);

    Apc = Ap->coeffs;
    Amc = Am->coeffs;

    for (i = 0; i < Alen; i++)
        _n_poly_mod_eval2_pow(Apc + i, Amc + i, Ac + i, alphapow, mod);

    Ap->length = Alen;
    _n_poly_normalise(Ap);
    Am->length = Alen;
    _n_poly_normalise(Am);
}

void n_bpoly_mod_interp_lift_2sm_poly(
    slong * deg1,
    n_bpoly_t T,
    const n_poly_t A,
    const n_poly_t B,    
    mp_limb_t alpha,
    nmod_t mod)
{
    slong i;
    slong lastlength = 0;
    const mp_limb_t * Acoeffs = A->coeffs;
    const mp_limb_t * Bcoeffs = B->coeffs;
    n_poly_struct * Tcoeffs;
    slong Alen = A->length;
    slong Blen = B->length;
    slong Tlen = FLINT_MAX(Alen, Blen);
    mp_limb_t d0 = (1 + mod.n)/2;
    mp_limb_t d1 = nmod_inv(nmod_add(alpha, alpha, mod), mod);
    mp_limb_t Avalue, Bvalue, u, v;

    n_bpoly_fit_length(T, Tlen);

    Tcoeffs = T->coeffs;

    for (i = 0; i < Tlen; i++)
    {
        Avalue = (i < Alen) ? Acoeffs[i] : 0;
        Bvalue = (i < Blen) ? Bcoeffs[i] : 0;
        u = nmod_sub(Avalue, Bvalue, mod);
        v = nmod_add(Avalue, Bvalue, mod);
        u = nmod_mul(u, d1, mod);
        v = nmod_mul(v, d0, mod);
        if ((u | v) == 0)
        {
            n_poly_zero(Tcoeffs + i);
        }
        else
        {
            n_poly_fit_length(Tcoeffs + i, 2);
            Tcoeffs[i].coeffs[0] = v;
            Tcoeffs[i].coeffs[1] = u;
            Tcoeffs[i].length = 1 + (u != 0);
            lastlength = FLINT_MAX(lastlength, Tcoeffs[i].length);
        }
    }

    *deg1 = lastlength - 1;

    FLINT_ASSERT(Tlen <= 0 || !n_poly_is_zero(Tcoeffs + Tlen - 1));
    T->length = Tlen;
}

int n_bpoly_mod_interp_crt_2sm_poly(
    slong * deg1,
    n_bpoly_t F,
    n_bpoly_t T,
    n_poly_t A,
    n_poly_t B,
    const n_poly_t modulus,
    n_poly_t alphapow,
    nmod_t mod)
{
    int changed = 0;
    slong i, lastlength = 0;
    slong Alen = A->length;
    slong Blen = B->length;
    slong Flen = F->length;
    slong Tlen = FLINT_MAX(FLINT_MAX(Alen, Blen), Flen);
    n_poly_struct * Tcoeffs, * Fcoeffs;
    mp_limb_t * Acoeffs, * Bcoeffs;
    n_poly_t zero;
    mp_limb_t Avalue, Bvalue, FvalueA, FvalueB, u, v;
    n_poly_struct * Fvalue;
    mp_limb_t alpha = alphapow->coeffs[1];

    zero->alloc = 0;
    zero->length = 0;
    zero->coeffs = NULL;

    n_bpoly_fit_length(T, Tlen);
    Tcoeffs = T->coeffs;
    Acoeffs = A->coeffs;
    Bcoeffs = B->coeffs;
    Fcoeffs = F->coeffs;

    for (i = 0; i < Tlen; i++)
    {
        Fvalue = (i < Flen) ? Fcoeffs + i : zero;
        _n_poly_mod_eval2_pow(&FvalueA, &FvalueB, Fcoeffs + i, alphapow, mod);
        Avalue = (i < Alen) ? Acoeffs[i] : 0;
        Bvalue = (i < Blen) ? Bcoeffs[i] : 0;
        FvalueA = nmod_sub(FvalueA, Avalue, mod);
        FvalueB = nmod_sub(FvalueB, Bvalue, mod);
        u = nmod_sub(FvalueB, FvalueA, mod);
        v = nmod_mul(mod.n - alpha, nmod_add(FvalueB, FvalueA, mod), mod);
        if ((u | v) != 0)
        {
            changed = 1;
            n_poly_mod_addmul_linear(Tcoeffs + i, Fvalue, modulus, u, v, mod);
        }
        else
        {
            n_poly_set(Tcoeffs + i, Fvalue);
        }

        lastlength = FLINT_MAX(lastlength, Tcoeffs[i].length);
    }

    FLINT_ASSERT(Tlen <= 0 || !n_poly_is_zero(Tcoeffs + Tlen - 1));
    T->length = Tlen;

    if (changed)
        n_bpoly_swap(T, F);

    FLINT_ASSERT(n_bpoly_mod_is_canonical(F, mod));

    *deg1 = lastlength - 1;
    return changed;
}


void n_polyu3_mod_interp_reduce_2sm_bpoly(
    n_bpoly_t Ap,
    n_bpoly_t Am,
    const n_polyu_t A,
    n_poly_t alpha,
    nmod_t mod)
{
    slong i;
    slong cur0, cur1, e0, e1, e2;
    ulong tp0, tp1, tp2, tm0, tm1, tm2, p1, p0;
    n_polyu_term_struct * Aterms = A->terms;
/*
flint_printf("n_polyu3_mod_interp_reduce_2sm_bpoly called alpha = %wd\n", alpha->coeffs[1]);
flint_printf("A: "); n_polyu3_print_pretty(A, "Y", "X", "Z"); flint_printf("\n");
*/
    n_bpoly_zero(Ap);
    n_bpoly_zero(Am);

    FLINT_ASSERT(A->length > 0);

    i = 0;

    cur0 = extract_exp(Aterms[i].exp, 2, 3);
    cur1 = extract_exp(Aterms[i].exp, 1, 3);
    e2   = extract_exp(Aterms[i].exp, 0, 3);

    n_poly_fill_powers(alpha, e2, mod);

    tp2 = tp1 = tp0 = tm2 = tm1 = tm0 = 0;
    if (e2 & 1)
        umul_ppmm(tm1, tm0, alpha->coeffs[e2], Aterms[i].coeff);
    else
        umul_ppmm(tp1, tp0, alpha->coeffs[e2], Aterms[i].coeff);

    for (i = 1; i < A->length; i++)
    {
        e0 = extract_exp(Aterms[i].exp, 2, 3);
        e1 = extract_exp(Aterms[i].exp, 1, 3);
        e2 = extract_exp(Aterms[i].exp, 0, 3);

        FLINT_ASSERT(e0 <= cur0);
        if (e0 < cur0 || e1 < cur1)
        {
            NMOD_RED3(tp0, tp2, tp1, tp0, mod);
            NMOD_RED3(tm0, tm2, tm1, tm0, mod);

            n_bpoly_set_coeff(Ap, cur0, cur1, nmod_add(tp0, tm0, mod));
            n_bpoly_set_coeff(Am, cur0, cur1, nmod_sub(tp0, tm0, mod));

            tp2 = tp1 = tp0 = tm2 = tm1 = tm0 = 0;
        }
        else
        {
            FLINT_ASSERT(e0 == cur0);
            FLINT_ASSERT(e1 == cur1);
        }

        cur0 = e0;
        cur1 = e1;

        n_poly_fill_powers(alpha, e2, mod);

        FLINT_ASSERT(alpha->coeffs[e2] == nmod_pow_ui(alpha->coeffs[1], e2, mod));

        umul_ppmm(p1, p0, alpha->coeffs[e2], Aterms[i].coeff);
        if (e2 & 1)
            add_sssaaaaaa(tm2, tm1, tm0, tm2, tm1, tm0, 0, p1, p0);
        else
            add_sssaaaaaa(tp2, tp1, tp0, tp2, tp1, tp0, 0, p1, p0);
    }

    NMOD_RED3(tp0, tp2, tp1, tp0, mod);
    NMOD_RED3(tm0, tm2, tm1, tm0, mod);
    n_bpoly_set_coeff(Ap, cur0, cur1, nmod_add(tp0, tm0, mod));
    n_bpoly_set_coeff(Am, cur0, cur1, nmod_sub(tp0, tm0, mod));
/*
flint_printf("n_polyu3_mod_interp_reduce_2sm_bpoly returning\n");
flint_printf("Ap: "); n_bpoly_print_pretty(Ap, "Y", "X"); flint_printf("\n");
flint_printf("Am: "); n_bpoly_print_pretty(Am, "Y", "X"); flint_printf("\n");
*/
}


/*
    T(x0, x1, x2) is in F[x2][x0, x1]
    A(x0, x1) are B(x0, x1) are in F[x0, x1]
    set T so that
        T(x0, x1, x2) == A(x0, x1) mod (x2 - alpha)
        T(x0, x1, x2) == B(x0, x1) mod (x2 + alpha)
*/

void n_polyu3n_mod_interp_lift_2sm_bpoly(
    slong * lastdeg,
    n_polyun_t T,
    const n_bpoly_t A,
    const n_bpoly_t B,    
    mp_limb_t alpha,
    nmod_t mod)
{
    slong lastlength = 0;
    n_polyun_term_struct * Tterms;
    slong Ti;
    n_poly_struct * Acoeffs = A->coeffs;
    slong Ai, ai;
    n_poly_struct * Bcoeffs = B->coeffs;
    slong Bi, bi;
    mp_limb_t u, v, Avalue, Bvalue;
    mp_limb_t d0, d1;
/*
flint_printf("n_polyu3n_mod_interp_lift_2sm_bpoly called\n");
flint_printf("A: "); n_bpoly_print_pretty(A, "Y", "X"); flint_printf("\n");
flint_printf("B: "); n_bpoly_print_pretty(B, "Y", "X"); flint_printf("\n");
*/
    FLINT_ASSERT(2*alpha < mod.n);

    d0 = (1 + mod.n)/2;
    d1 = nmod_inv(nmod_add(alpha, alpha, mod), mod);

    n_polyun_fit_length(T, FLINT_MAX(A->length, B->length));
    Tterms = T->terms;

    Ti = 0;
    Ai = A->length - 1;
    Bi = B->length - 1;
    ai = (Ai < 0) ? 0 : n_poly_degree(A->coeffs + Ai);
    bi = (Bi < 0) ? 0 : n_poly_degree(B->coeffs + Bi);

    while (Ai >= 0 || Bi >= 0)
    {
        if (Ti >= T->alloc)
        {
            slong extra = FLINT_MAX(Ai, Bi);
            n_polyun_realloc(T, Ti + extra + 1);
            Tterms = T->terms;
        }

        FLINT_ASSERT(Ai < 0 || Acoeffs[Ai].coeffs[ai] != 0);
        FLINT_ASSERT(Bi < 0 || Bcoeffs[Bi].coeffs[bi] != 0);

        Avalue = 0;
        if (Ai >= 0)
        {
            Avalue = Acoeffs[Ai].coeffs[ai];
            Tterms[Ti].exp = pack_exp3(Ai, ai, 0);
        }

        Bvalue = 0;
        if (Bi >= 0)
        {
            ulong Bexp = pack_exp3(Bi, bi, 0);
            if (Avalue == 0)
            {
                Bvalue = Bcoeffs[Bi].coeffs[bi];
                Tterms[Ti].exp = Bexp;
            }
            else
            {
                if (Tterms[Ti].exp <= Bexp)
                {
                    Bvalue = Bcoeffs[Bi].coeffs[bi];
                }
                if (Tterms[Ti].exp < Bexp)
                {
                    Avalue = 0;
                    Tterms[Ti].exp = Bexp;
                }
            }
        }
/*
flint_printf("Texp: X^%wd*Y^%wd*?^%wd\n",
    extract_exp(Tterms[Ti].exp, 2, 3),
    extract_exp(Tterms[Ti].exp, 1, 3),
    extract_exp(Tterms[Ti].exp, 0, 3));
flint_printf("Avalue: %wu\n", Avalue);
flint_printf("Bvalue: %wu\n", Bvalue);
*/
        u = nmod_sub(Avalue, Bvalue, mod);
        v = nmod_add(Avalue, Bvalue, mod);
        u = nmod_mul(u, d1, mod);
        v = nmod_mul(v, d0, mod);
/*
flint_printf("u: %wu\n", u);
flint_printf("v: %wu\n", v);
*/
        FLINT_ASSERT(u != 0 || v != 0);

        n_poly_fit_length(Tterms[Ti].coeff, 2);
        Tterms[Ti].coeff->coeffs[0] = v;
        Tterms[Ti].coeff->coeffs[1] = u;
        Tterms[Ti].coeff->length = 1 + (u != 0);
        lastlength = FLINT_MAX(lastlength, Tterms[Ti].coeff->length);
        Ti++;

        if (Avalue != 0)
        {
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
        if (Bvalue != 0)
        {
            do {
                bi--;
            } while (bi >= 0 && Bcoeffs[Bi].coeffs[bi] == 0);
            if (bi < 0)
            {
                do {
                    Bi--;
                } while (Bi >= 0 && Bcoeffs[Bi].length == 0);
                if (Bi >= 0)
                    bi = n_poly_degree(Bcoeffs + Bi);
            }
        }
    }
    T->length = Ti;

    FLINT_ASSERT(n_polyun_mod_is_canonical(T, mod));

    *lastdeg = lastlength - 1;
/*
flint_printf("n_polyu3n_mod_interp_lift_2sm_bpoly returning lastdeg = %wd\n", lastlength - 1);
flint_printf("T: "); n_polyu3n_print_pretty(T, "Y", "X", "?", "Z"); flint_printf("\n");
*/

    return;
}


int n_polyu3n_mod_interp_crt_2sm_bpoly(
    slong * lastdeg,
    n_polyun_t F,
    n_polyun_t T,
    n_bpoly_t A,
    n_bpoly_t B,
    const n_poly_t modulus,
    n_poly_t alphapow,
    nmod_t mod)
{
    int changed = 0;
    slong lastlength = 0;
    n_poly_t zero;
    n_polyun_term_struct * Tterms;
    slong Ti;
    n_polyun_term_struct * Fterms = F->terms;
    slong Flen = F->length;
    slong Fi;
    n_poly_struct * Acoeffs = A->coeffs;
    slong Ai, ai;
    n_poly_struct * Bcoeffs = B->coeffs;
    slong Bi, bi;
    n_poly_struct * Fvalue;
    mp_limb_t u, v, Avalue, Bvalue, FvalueA, FvalueB;
    int texp_set, cmp;
    mp_limb_t alpha = alphapow->coeffs[1];
/*
flint_printf("n_polyu3n_mod_interp_crt_2sm_bpoly called alpha = %wu\n", alpha);
flint_printf("+: "); n_bpoly_print_pretty(A, "Y", "X"); flint_printf("\n");
flint_printf("-: "); n_bpoly_print_pretty(B, "Y", "X"); flint_printf("\n");
*/
#if WANT_ASSERT
    u = n_poly_mod_evaluate_nmod(modulus, alpha, mod);
    u = nmod_mul(u, alpha, mod);
    u = nmod_mul(u, 2, mod);
    FLINT_ASSERT(u == 1);
    u = n_poly_mod_evaluate_nmod(modulus, mod.n - alpha, mod);
    u = nmod_mul(u, alpha, mod);
    u = nmod_mul(u, 2, mod);
    FLINT_ASSERT(u == 1);
#endif

    FLINT_ASSERT(n_bpoly_mod_is_canonical(A, mod));
    FLINT_ASSERT(n_bpoly_mod_is_canonical(B, mod));
    FLINT_ASSERT(n_polyun_mod_is_canonical(F, mod));

    n_polyun_fit_length(T, FLINT_MAX(Flen, A->length));
    Tterms = T->terms;

    zero->alloc = 0;
    zero->length = 0;
    zero->coeffs = NULL;

    Ti = Fi = 0;
    Ai = A->length - 1;
    Bi = B->length - 1;
    ai = (Ai < 0) ? 0 : n_poly_degree(A->coeffs + Ai);
    bi = (Bi < 0) ? 0 : n_poly_degree(B->coeffs + Bi);

    while (Fi < Flen || Ai >= 0 || Bi >= 0)
    {
        if (Ti >= T->alloc)
        {
            slong extra = Flen - Fi;
            extra = FLINT_MAX(extra, Ai);
            extra = FLINT_MAX(extra, Bi);
            n_polyun_fit_length(T, Ti + extra + 1);
            Tterms = T->terms;
        }

        FLINT_ASSERT(Fi >= Flen || Fterms[Fi].coeff->length > 0);
        FLINT_ASSERT(Ai < 0 || Acoeffs[Ai].coeffs[ai] != 0);
        FLINT_ASSERT(Bi < 0 || Bcoeffs[Bi].coeffs[bi] != 0);

        Fvalue = zero;
        texp_set = 0;
        if (Fi < Flen)
        {
            Fvalue = Fterms[Fi].coeff;
            texp_set = 1;
            Tterms[Ti].exp = Fterms[Fi].exp;
        }

        Avalue = 0;
        if (Ai >= 0)
        {
            ulong Aexp = pack_exp3(Ai, ai, 0);
            cmp = (!texp_set) ? -1 :
                  Tterms[Ti].exp < Aexp ? -1 :
                  Tterms[Ti].exp > Aexp ? 1 : 0;
            if (cmp <= 0)
            {
                Avalue = Acoeffs[Ai].coeffs[ai];
            }
            if (cmp < 0)
            {
                Fvalue = zero;
                texp_set = 1;
                Tterms[Ti].exp = Aexp;
            }
        }

        Bvalue = 0;
        if (Bi >= 0)
        {
            ulong Bexp = pack_exp3(Bi, bi, 0);
            cmp = (!texp_set) ? -1 :
                  Tterms[Ti].exp < Bexp ? -1 :
                  Tterms[Ti].exp > Bexp ? 1 : 0;
            if (cmp <= 0)
            {
                Bvalue = Bcoeffs[Bi].coeffs[bi];
            }
            if (cmp < 0)
            {
                Fvalue = zero;
                Avalue = 0;
                texp_set = 1;
                Tterms[Ti].exp = Bexp;
            }
        }

        FLINT_ASSERT(texp_set);

/*
flint_printf("Texp: X^%wd*Y^%wd*?^%wd\n",
    extract_exp(Tterms[Ti].exp, 2, 3),
    extract_exp(Tterms[Ti].exp, 1, 3),
    extract_exp(Tterms[Ti].exp, 0, 3));
flint_printf("+value: %wu\n", Avalue);
flint_printf("-value: %wu\n", Bvalue);
flint_printf("Fvalue: "); n_poly_print_pretty(Fvalue, "Z"); flint_printf("\n");
*/

        _n_poly_mod_eval2_pow(&FvalueA, &FvalueB, Fvalue, alphapow, mod);
        FvalueA = nmod_sub(FvalueA, Avalue, mod);
        FvalueB = nmod_sub(FvalueB, Bvalue, mod);
        u = nmod_sub(FvalueB, FvalueA, mod);
        v = nmod_mul(mod.n - alpha, nmod_add(FvalueB, FvalueA, mod), mod);
        if (u != 0 || v != 0)
        {
            changed = 1;
            n_poly_mod_addmul_linear(Tterms[Ti].coeff, Fvalue, modulus,
                                                               u, v, mod);
        }
        else
        {
            FLINT_ASSERT(!n_poly_is_zero(Fvalue));
            n_poly_set(Tterms[Ti].coeff, Fvalue);
        }

        FLINT_ASSERT(Tterms[Ti].coeff->length >= Fvalue->length);
/*
flint_printf("Tvalue: "); n_poly_print_pretty(Tterms[Ti].coeff, "Z"); flint_printf("\n");
*/
        Fi += (Fvalue != zero);
        if (Avalue != 0)
        {
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
        if (Bvalue != 0)
        {
            do {
                bi--;
            } while (bi >= 0 && Bcoeffs[Bi].coeffs[bi] == 0);
            if (bi < 0)
            {
                do {
                    Bi--;
                } while (Bi >= 0 && Bcoeffs[Bi].length == 0);
                if (Bi >= 0)
                    bi = n_poly_degree(Bcoeffs + Bi);
            }
        }

        FLINT_ASSERT(!n_poly_is_zero(Tterms[Ti].coeff));
        lastlength = FLINT_MAX(lastlength, Tterms[Ti].coeff->length);
        Ti++;
    }
    T->length = Ti;

    if (changed)
        n_polyun_swap(T, F);

    FLINT_ASSERT(n_polyun_mod_is_canonical(F, mod));

    *lastdeg = lastlength - 1;
    return changed;
}



int n_bpoly_mod_disolve2(
    n_bpoly_t C1, n_bpoly_t C2,
    slong C1_deg1_bound, slong C2_deg1_bound,
    n_bpoly_t A,
    n_bpoly_t B1, n_bpoly_t B2,
    nmod_t mod)
{
    int success;
    slong A_deg1, B1_deg1, B2_deg1, C1_deg1, C2_deg1;
    slong bad_prime_count, bound;
    mp_limb_t alpha, c;
    n_poly_t Aevalp, B1evalp, B2evalp, C1evalp, C2evalp;
    n_poly_t Aevalm, B1evalm, B2evalm, C1evalm, C2evalm;
    n_poly_t modulus, alphapow, t1, t2;
    n_bpoly_t T;
/*
    const char * vars[] = {"x", "y"};
*/
    n_poly_init(Aevalp);
    n_poly_init(B1evalp);
    n_poly_init(B2evalp);
    n_poly_init(C1evalp);
    n_poly_init(C2evalp);
    n_poly_init(Aevalm);
    n_poly_init(B1evalm);
    n_poly_init(B2evalm);
    n_poly_init(C1evalm);
    n_poly_init(C2evalm);
    n_poly_init(modulus);
    n_poly_init(alphapow);
    n_poly_init(t1);
    n_poly_init(t2);

    n_bpoly_init(T);

    A_deg1 = n_bpoly_degree1(A);
    B1_deg1 = n_bpoly_degree1(B1);
    B2_deg1 = n_bpoly_degree1(B2);
/*
flint_printf("n_bpoly_mod_disolve2 called p = %wu\n", mod.n);
flint_printf(" A(deg1 %wd): ", A_deg1); n_bpoly_print_pretty(A, vars[0], vars[1]); flint_printf("\n");
flint_printf("B1(deg1 %wd): ", B1_deg1); n_bpoly_print_pretty(B1, vars[0], vars[1]); flint_printf("\n");
flint_printf("B2(deg1 %wd): ", B2_deg1); n_bpoly_print_pretty(B2, vars[0], vars[1]); flint_printf("\n");
*/
    if (B1_deg1 < 1 || B2_deg1 < 1)
    {
        success = -2;
        goto cleanup;
    }

    bound = A_deg1;
    bound = FLINT_MAX(bound, C1_deg1_bound + B2_deg1);
    bound = FLINT_MAX(bound, B1_deg1 + C2_deg1_bound);
    bound += 1;

    bad_prime_count = 0;

    n_poly_fit_length(alphapow, FLINT_MAX(WORD(3), bound + 1));
    n_poly_one(modulus);

    if ((mod.n & UWORD(1)) == UWORD(0))
    {
        success = -1;
        goto cleanup;
    }

    alpha = (mod.n - UWORD(1))/UWORD(2);

    goto choose_prime;

bad_prime:

    if (bad_prime_count > bound)
    {
        success = n_poly_degree(modulus) > 0 ? -1 : -2;
        goto cleanup;
    }

    bad_prime_count++;

choose_prime:

    if (alpha < 2)
    {
        success = -1;
        goto cleanup;
    }

    alpha--;

    FLINT_ASSERT(0 < alpha && alpha <= mod.n/2);
    FLINT_ASSERT(alphapow->alloc >= 2);
    alphapow->length = 2;
    alphapow->coeffs[0] = 1;
    alphapow->coeffs[1] = alpha;

    n_bpoly_mod_interp_reduce_2sm_poly(Aevalp, Aevalm, A, alphapow, mod);
    n_bpoly_mod_interp_reduce_2sm_poly(B1evalp, B1evalm, B1, alphapow, mod);
    n_bpoly_mod_interp_reduce_2sm_poly(B2evalp, B2evalm, B2, alphapow, mod);

    /* make sure evaluation point did not drop the degree of a Bi */
    if (n_poly_degree(B1evalp) < n_bpoly_degree0(B1) ||
        n_poly_degree(B1evalm) < n_bpoly_degree0(B1) ||
        n_poly_degree(B2evalp) < n_bpoly_degree0(B2) ||
        n_poly_degree(B2evalm) < n_bpoly_degree0(B2))
    {
        goto choose_prime;
    }

    /* image pfrac's */
    if (!n_poly_mod_invmod(t1, B2evalp, B1evalp, mod))
        goto bad_prime;
    _n_poly_mod_mul(t2, Aevalp, t1, mod);
    _n_poly_mod_rem(C1evalp, t2, B1evalp, mod);
    _n_poly_mod_mul(t2, B2evalp, C1evalp, mod);
    n_poly_mod_sub(Aevalp, Aevalp, t2, mod);
    _n_poly_mod_div(C2evalp, Aevalp, B1evalp, mod);

    if (!n_poly_mod_invmod(t1, B2evalm, B1evalm, mod))
        goto bad_prime;
    _n_poly_mod_mul(t2, Aevalm, t1, mod);
    _n_poly_mod_rem(C1evalm, t2, B1evalm, mod);
    _n_poly_mod_mul(t2, B2evalm, C1evalm, mod);
    n_poly_mod_sub(Aevalm, Aevalm, t2, mod);
    _n_poly_mod_div(C2evalm, Aevalm, B1evalm, mod);

    /* update interpolants */
    if (n_poly_degree(modulus) > 0)
    {
        c = n_poly_mod_evaluate_nmod(modulus, alpha, mod);
        FLINT_ASSERT(c == n_poly_mod_evaluate_nmod(modulus, mod.n - alpha, mod));
        c = nmod_mul(c, alpha, mod);
        c = nmod_add(c, c, mod);
        c = nmod_inv(c, mod);
        n_poly_mod_scalar_mul_nmod(modulus, modulus, c, mod);

        n_bpoly_mod_interp_crt_2sm_poly(&C1_deg1, C1, T, C1evalp, C1evalm, modulus, alphapow, mod);
        n_bpoly_mod_interp_crt_2sm_poly(&C2_deg1, C2, T, C2evalp, C2evalm, modulus, alphapow, mod);
    }
    else
    {
        n_bpoly_mod_interp_lift_2sm_poly(&C1_deg1, C1, C1evalp, C1evalm, alpha, mod);
        n_bpoly_mod_interp_lift_2sm_poly(&C2_deg1, C2, C2evalp, C2evalm, alpha, mod);
    }

    c = mod.n - nmod_mul(alpha, alpha, mod);
    n_poly_mod_shift_left_scalar_addmul(modulus, 2, c, mod);

    if (C1_deg1 > C1_deg1_bound || C2_deg1 > C2_deg1_bound)
    {
        success = 0;
        goto cleanup;
    }

    if (n_poly_degree(modulus) < bound)
    {
        goto choose_prime;
    }

    success = 1;

cleanup:

#if WANT_ASSERT
    if (success == 1)
    {
        n_bpoly_t T1, T2, T3;
        n_bpoly_init(T1);
        n_bpoly_init(T2);
        n_bpoly_init(T3);
        n_bpoly_set(T1, A);
        n_bpoly_mod_mul(T2, C1, B2, mod);
        n_bpoly_mod_sub(T3, T1, T2, mod);
        n_bpoly_swap(T1, T3);
        n_bpoly_mod_mul(T2, C2, B1, mod);
        n_bpoly_mod_sub(T3, T1, T2, mod);
        n_bpoly_swap(T1, T3);
        FLINT_ASSERT(n_bpoly_is_zero(T1));
        n_bpoly_clear(T1);
        n_bpoly_clear(T2);
        n_bpoly_clear(T3);
    }
#endif

    n_poly_clear(Aevalp);
    n_poly_clear(B1evalp);
    n_poly_clear(B2evalp);
    n_poly_clear(C1evalp);
    n_poly_clear(C2evalp);
    n_poly_clear(Aevalm);
    n_poly_clear(B1evalm);
    n_poly_clear(B2evalm);
    n_poly_clear(C1evalm);
    n_poly_clear(C2evalm);
    n_poly_clear(modulus);
    n_poly_clear(alphapow);
    n_poly_clear(t1);
    n_poly_clear(t2);

    n_bpoly_clear(T);
/*
flint_printf("n_bpoly_mod_disolve2 returning success = %d\n", success);
flint_printf("C1: "); n_bpoly_print_pretty(C1, vars[0], vars[1]); flint_printf("\n");
flint_printf("C2: "); n_bpoly_print_pretty(C2, vars[0], vars[1]); flint_printf("\n");
*/
    return success;
}


/*
    Try to solve A/(B[0]*...*B[r-1]) = C[0]/B[0] + ... + C[r-1]/B[r-1]
    for the C[i] in Fp[y][x].
    return:
        1: solution found with deg_y(Ci) <= ldegCibound
        0: no solution exists with deg_y(Ci) <= ldegCibound
       -1: could not find enough evaluation points where the Bi are pariwise prime
       -2: found no evaluation points where the Bi are pariwise prime
*/
int n_bpoly_mod_disolve(
    slong r,
    n_bpoly_struct * C,
    slong * C_deg1_bound,
    n_bpoly_t A,
    n_bpoly_struct * B,
    nmod_t mod)
{
    int success;
    slong i, j, bad_prime_count, bound;
    mp_limb_t alpha, c;
    n_poly_struct Aevalp[1], * Bevalp, * Cevalp;
    n_poly_struct Aevalm[1], * Bevalm, * Cevalm;
    n_poly_t modulus, alphapow, t1, t2;
    n_bpoly_t T;
    slong * B_deg1, * C_deg1, B_deg1_total, A_deg1;
/*    const char * vars[] = {"x", "y"};*/
    TMP_INIT;

    if (r < 3)
    {
        return n_bpoly_mod_disolve2(C + 0, C + 1, C_deg1_bound[0],
                                        C_deg1_bound[1], A, B + 0, B + 1, mod);
    }

    TMP_START;

    B_deg1 = (slong *) TMP_ALLOC(2*r*sizeof(slong));
    C_deg1 = B_deg1 + r;

    Bevalp = (n_poly_struct *) TMP_ALLOC(4*r*sizeof(n_poly_struct));
    Bevalm = Bevalp + r;
    Cevalp = Bevalm + r;
    Cevalm = Cevalp + r;

    n_poly_init(Aevalp);
    n_poly_init(Aevalm);
    n_poly_init(modulus);
    n_poly_init(alphapow);
    n_poly_init(t1);
    n_poly_init(t2);
    for (i = 0; i < r; i++)
    {
        n_poly_init(Bevalp + i);
        n_poly_init(Bevalm + i);
        n_poly_init(Cevalp + i);
        n_poly_init(Cevalm + i);
    }

    n_bpoly_init(T);

    A_deg1 = n_bpoly_degree1(A);
    for (i = 0; i < r; i++)
    {
        B_deg1[i] = n_bpoly_degree1(B + i);
        if (B_deg1[i] < 1)
        {
            success = -2;
            goto cleanup;
        }
    }

/*
printf("n_bpoly_mod_disolve called\n");
flint_printf(" A(deg %wd): ", A_deg1); n_bpoly_print_pretty(A, vars[0], vars[1]); flint_printf("\n");
for (i = 0; i < r; i++)
{
flint_printf("B[%wd](deg %wd): ", i, B_deg1[i]); n_bpoly_print_pretty(B + i, vars[0], vars[1]); flint_printf("\n");
}
*/

    B_deg1_total = B_deg1[0];
    for (i = 1; i < r; i++)
        B_deg1_total += B_deg1[i];

    bound = A_deg1;
    for (i = 0; i < r; i++)
        bound = FLINT_MAX(bound, C_deg1_bound[i] + B_deg1_total - B_deg1[i]);
    bound += 1;

    bad_prime_count = 0;

    n_poly_fit_length(alphapow, FLINT_MAX(WORD(3), bound + 1));
    n_poly_one(modulus);

    if ((mod.n & UWORD(1)) == UWORD(0))
    {
        success = -1;
        goto cleanup;
    }

    alpha = (mod.n - UWORD(1))/UWORD(2);

    goto choose_prime;

bad_prime:

    if (bad_prime_count > 2*bound)
    {
        success = n_poly_degree(modulus) > 0 ? -1 : -2;
        goto cleanup;
    }

    bad_prime_count++;

choose_prime:

    if (alpha < 2)
    {
        success = -1;
        goto cleanup;
    }

    alpha--;

    FLINT_ASSERT(0 < alpha && alpha <= mod.n/2);
    FLINT_ASSERT(alphapow->alloc >= 2);
    alphapow->length = 2;
    alphapow->coeffs[0] = 1;
    alphapow->coeffs[1] = alpha;

    n_bpoly_mod_interp_reduce_2sm_poly(Aevalp, Aevalm, A, alphapow, mod);
    for (i = 0; i < r; i++)
    {
        n_bpoly_mod_interp_reduce_2sm_poly(Bevalp + i, Bevalm + i, B + i, alphapow, mod);
        if (n_poly_degree(Bevalp + i) < n_bpoly_degree0(B + i) ||
            n_poly_degree(Bevalm + i) < n_bpoly_degree0(B + i))
        {
            goto choose_prime;
        }
    }

    for (i = 0; i < r; i++)
    {
        n_poly_one(t2);
        for (j = 0; j < r; j++)
        {
            if (j == i)
                continue;
            _n_poly_mod_mul(t1, t2, Bevalp + j, mod);
            n_poly_swap(t1, t2);
        }
        if (!n_poly_mod_invmod(t1, t2, Bevalp + i, mod))
            goto bad_prime;
        _n_poly_mod_mul(t2, Aevalp, t1, mod);
        _n_poly_mod_rem(Cevalp + i, t2, Bevalp + i, mod);

        n_poly_one(t2);
        for (j = 0; j < r; j++)
        {
            if (j == i)
                continue;
            _n_poly_mod_mul(t1, t2, Bevalm + j, mod);
            n_poly_swap(t1, t2);
        }
        if (!n_poly_mod_invmod(t1, t2, Bevalm + i, mod))
            goto bad_prime;
        _n_poly_mod_mul(t2, Aevalm, t1, mod);
        _n_poly_mod_rem(Cevalm + i, t2, Bevalm + i, mod);
    }

    if (n_poly_degree(modulus) > 0)
    {
        c = n_poly_mod_evaluate_nmod(modulus, alpha, mod);
        FLINT_ASSERT(c == n_poly_mod_evaluate_nmod(modulus, mod.n - alpha, mod));
        c = nmod_mul(c, alpha, mod);
        c = nmod_add(c, c, mod);
        c = nmod_inv(c, mod);
        n_poly_mod_scalar_mul_nmod(modulus, modulus, c, mod);

        for (i = 0; i < r; i++)
        {
            n_bpoly_mod_interp_crt_2sm_poly(C_deg1 + i, C + i, T, Cevalp + i,
                                           Cevalm + i, modulus, alphapow, mod);
        }
    }
    else
    {
        for (i = 0; i < r; i++)
        {
            n_bpoly_mod_interp_lift_2sm_poly(C_deg1 + i, C + i, Cevalp + i,
                                                       Cevalm + i, alpha, mod);
        }
    }

    c = mod.n - nmod_mul(alpha, alpha, mod);
    n_poly_mod_shift_left_scalar_addmul(modulus, 2, c, mod);

    for (i = 0; i < r; i++)
    {
        if (C_deg1[i] > C_deg1_bound[i])
        {
            success = 0;
            goto cleanup;
        }
    }

    if (n_poly_degree(modulus) < bound)
    {
        goto choose_prime;
    }

    success = 1;

cleanup:

#if WANT_ASSERT
    if (success == 1)
    {
        n_bpoly_t T1, T2, T3;
        n_bpoly_init(T1);
        n_bpoly_init(T2);
        n_bpoly_init(T3);
        n_bpoly_set(T1, A);
        for (i = 0; i < j; i++)
        {
            n_bpoly_one(T2);
            for (j = 0; j < r; j++)
            {
                n_bpoly_mod_mul(T3, T2, j == i ? C + j : B + j, mod);
                n_bpoly_swap(T2, T3);
            }
            n_bpoly_mod_sub(T3, T1, T2, mod);
            n_bpoly_swap(T1, T3);
        }
        FLINT_ASSERT(n_bpoly_is_zero(T1));
        n_bpoly_clear(T1);
        n_bpoly_clear(T2);
        n_bpoly_clear(T3);
    }
#endif

    n_poly_clear(Aevalp);
    n_poly_clear(Aevalm);
    n_poly_clear(modulus);
    n_poly_clear(alphapow);
    n_poly_clear(t1);
    n_poly_clear(t2);
    for (i = 0; i < r; i++)
    {
        n_poly_clear(Bevalp + i);
        n_poly_clear(Bevalm + i);
        n_poly_clear(Cevalp + i);
        n_poly_clear(Cevalm + i);
    }

    n_bpoly_clear(T);

    TMP_END;
/*
flint_printf("n_bpoly_mod_disolve returning success = %d\n", success);
for (i = 0; i < r; i++)
{
flint_printf("C[%wd]: ", i); n_bpoly_print_pretty(C + i, vars[0], vars[1]); flint_printf("\n");
}
*/
    return success;
}




/*
    Input A, B0, B1, with A(y,x,z) = B0(y,x,z)*B1(y,x,z) mod (y-beta)

    return
       -1: suspect that B0(beta,x,z) & B1(beta,x,z) are not pairwise prime, i.e. could not find an eval point for z making them pairwise prime, or
           A(y,x,z) has wrong degree in x, or
           could not find enough eval points for z with the required degree in x
        0: lift of B0*B1 to true factors is impossible
        1: successfully lifted B0*B1 to true factors BB0*BB1 without changing lc_x
*/
int n_polyu3_mod_factor_lift2(
    n_polyun_t BB0,
    n_polyun_t BB1,
    n_polyu_t A,
    n_polyu_t B0,
    n_polyu_t B1,
    mp_limb_t beta,
    slong degree_inner, /* required degree in x */
    const nmodf_ctx_t ctx)
{
    int success;
    n_polyun_t T;
    n_bpoly_t Ap, Am, B0p, B0m, B1p, B1m;
    n_poly_t modulus, alphapow, t1, t2;
    mp_limb_t alpha, c;
    slong ldegBB0, ldegBB1;
    slong Adegy, Adegz, Adegx;
    slong bad_primes_left;
/*
    const char * vars[] = {"Y", "X", "Z"};
flint_printf("+++++++++++++++++++++++\n");
flint_printf("n_polyu3_mod_factor_lift2 called: degree_inner = %wd, beta = %wd\n", degree_inner, beta);
flint_printf("A: "); n_polyu3_print_pretty(A, vars[0], vars[1], vars[2]); flint_printf("\n");
flint_printf("B0: "); n_polyu3_print_pretty(B0, vars[0], vars[1], vars[2]); flint_printf("\n");
flint_printf("B1: "); n_polyu3_print_pretty(B1, vars[0], vars[1], vars[2]); flint_printf("\n");
*/

    FLINT_ASSERT(n_polyu_mod_is_canonical(A, ctx->mod));
    FLINT_ASSERT(n_polyu_mod_is_canonical(B0, ctx->mod));
    FLINT_ASSERT(n_polyu_mod_is_canonical(B1, ctx->mod));

    n_polyun_init(T);
    n_bpoly_init(Ap);
    n_bpoly_init(Am);
    n_bpoly_init(B0p);
    n_bpoly_init(B0m);
    n_bpoly_init(B1p);
    n_bpoly_init(B1m);
    n_poly_init(modulus);
    n_poly_init(alphapow);
    n_poly_init(t1);
    n_poly_init(t2);

    n_polyu3_degrees(&Adegy, &Adegx, &Adegz, A);

    if (Adegx != degree_inner)
    {
        success = -1;
        goto cleanup;
    }

    n_poly_fit_length(alphapow, FLINT_MAX(WORD(3), Adegz + 2));
    n_poly_one(modulus);

    FLINT_ASSERT((ctx->mod.n & UWORD(1)) == UWORD(1));

    alpha = (ctx->mod.n - UWORD(1))/UWORD(2);

    bad_primes_left = FLINT_MAX(5, Adegz);

    goto choose_prime;

choose_prime:

    if (alpha < 2)
    {
        success = -1;
        goto cleanup;
    }

    alpha--;
/*
flint_printf("------ Z^2 - alpha^2: alpha = %wu ------\n", alpha);
usleep(1000000);
*/

    FLINT_ASSERT(0 < alpha && alpha <= ctx->mod.n/2);
    FLINT_ASSERT(alphapow->alloc >= 2);
    alphapow->length = 2;
    alphapow->coeffs[0] = 1;
    alphapow->coeffs[1] = alpha;

    n_polyu3_mod_interp_reduce_2sm_bpoly(Ap, Am, A, alphapow, ctx->mod);
    n_polyu3_mod_interp_reduce_2sm_bpoly(B0p, B0m, B0, alphapow, ctx->mod);
    n_polyu3_mod_interp_reduce_2sm_bpoly(B1p, B1m, B1, alphapow, ctx->mod);
/*
flint_printf(" Ap: "); n_bpoly_print_pretty(Ap, vars[0], vars[1]); flint_printf("\n");
flint_printf("B0p: "); n_bpoly_print_pretty(B0p, vars[0], vars[1]); flint_printf("\n");
flint_printf("B1p: "); n_bpoly_print_pretty(B1p, vars[0], vars[1]); flint_printf("\n");
*/
    success = n_bpoly_mod_hensel_lift2(Ap, B0p, B1p, beta, degree_inner, ctx->mod);
    if (success <= 0)
    {
        if (success == 0 || --bad_primes_left < 0)
            goto cleanup;
        goto choose_prime;
    }
/*
flint_printf("B0+: "); n_bpoly_print_pretty(B0p, vars[0], vars[1]); flint_printf("\n");
flint_printf("B1+: "); n_bpoly_print_pretty(B1p, vars[0], vars[1]); flint_printf("\n");
*/
/*
flint_printf(" Am: "); n_bpoly_print_pretty(Am, vars[0], vars[1]); flint_printf("\n");
flint_printf("B0m: "); n_bpoly_print_pretty(B0m, vars[0], vars[1]); flint_printf("\n");
flint_printf("B1m: "); n_bpoly_print_pretty(B1m, vars[0], vars[1]); flint_printf("\n");
*/
    success = n_bpoly_mod_hensel_lift2(Am, B0m, B1m, beta, degree_inner, ctx->mod);
    if (success <= 0)
    {
        if (success == 0 || --bad_primes_left < 0)
            goto cleanup;
        goto choose_prime;
    }
/*
flint_printf("B0-: "); n_bpoly_print_pretty(B0m, vars[0], vars[1]); flint_printf("\n");
flint_printf("B1-: "); n_bpoly_print_pretty(B1m, vars[0], vars[1]); flint_printf("\n");
*/
    if (n_poly_degree(modulus) > 0)
    {
        c = n_poly_mod_evaluate_nmod(modulus, alpha, ctx->mod);
        FLINT_ASSERT(c == n_poly_mod_evaluate_nmod(modulus,
                                                ctx->mod.n - alpha, ctx->mod));
        c = nmod_mul(c, alpha, ctx->mod);
        c = nmod_add(c, c, ctx->mod);
        c = n_invmod(c, ctx->mod.n);
        n_poly_mod_scalar_mul_nmod(modulus, modulus, c, ctx->mod);
        n_polyu3n_mod_interp_crt_2sm_bpoly(&ldegBB0, BB0, T, B0p, B0m,
                                                  modulus, alphapow, ctx->mod);
        n_polyu3n_mod_interp_crt_2sm_bpoly(&ldegBB1, BB1, T, B1p, B1m,
                                                  modulus, alphapow, ctx->mod);
    }
    else
    {
        n_polyu3n_mod_interp_lift_2sm_bpoly(&ldegBB0, BB0, B0p, B0m, alpha, ctx->mod);
        n_polyu3n_mod_interp_lift_2sm_bpoly(&ldegBB1, BB1, B1p, B1m, alpha, ctx->mod);
    }

    c = ctx->mod.n - nmod_mul(alpha, alpha, ctx->mod);
    n_poly_mod_shift_left_scalar_addmul(modulus, 2, c, ctx->mod);
/*
flint_printf("ldegBB0: %wd\n", ldegBB0);
flint_printf("ldegBB1: %wd\n", ldegBB1);
flint_printf("modulus: "); n_poly_print_pretty(modulus, "Z"); printf("\n");
flint_printf("BB0: "); n_polyu3n_print_pretty(BB0, vars[0], vars[1], "?", vars[2]); flint_printf("\n");
flint_printf("BB1: "); n_polyu3n_print_pretty(BB1, vars[0], vars[1], "?", vars[2]); flint_printf("\n");
*/
    if (ldegBB0 + ldegBB1 > Adegz)
    {
        success = 0;
        goto cleanup;
    }

    if (n_poly_degree(modulus) <= Adegz)
    {
        goto choose_prime;
    }

    success = 1;

cleanup:

/*
flint_printf("n_polyu3n_hensel_lift2 returning %d\n", success);
flint_printf("BB0: "); n_polyu3n_print_pretty(BB0, vars[0], vars[1], "?", vars[2]); flint_printf("\n");
flint_printf("BB1: "); n_polyu3n_print_pretty(BB1, vars[0], vars[1], "?", vars[2]); flint_printf("\n");
flint_printf("+++++++++++++++++++++++\n");
*/

#if WANT_ASSERT
    if (success == 1)
    {
        n_polyu_t T1, T2, T3;
        n_polyu_init(T1);
        n_polyu_init(T2);
        n_polyu_init(T3);
        n_polyu_get_n_polyun(T2, BB0);
        n_polyu_get_n_polyun(T3, BB1);
        n_polyu_mod_mul(T1, T2, T3, ctx->mod);
        FLINT_ASSERT(n_polyu_equal(A, T1));
        n_polyu_clear(T1);
        n_polyu_clear(T2);
        n_polyu_clear(T3);
    }
#endif

    n_polyun_clear(T);
    n_bpoly_clear(Ap);
    n_bpoly_clear(Am);
    n_bpoly_clear(B0p);
    n_bpoly_clear(B0m);
    n_bpoly_clear(B1p);
    n_bpoly_clear(B1m);
    n_poly_clear(modulus);
    n_poly_clear(alphapow);
    n_poly_clear(t1);
    n_poly_clear(t2);

    return success;
}
/* r factor version */
int n_polyu3_mod_factor_lift(
    slong r,
    n_polyun_struct * BB,
    n_polyu_t A,
    n_polyu_struct * B,
    mp_limb_t beta,
    slong degree_inner, /* required degree in x */
    const nmodf_ctx_t ctx)
{
    int success;
    slong i, j;
    n_polyun_t T;
    n_bpoly_struct * Bp, * Bm;
    n_bpoly_t Ap, Am;
    n_poly_t modulus, alphapow, t1, t2;
    mp_limb_t alpha, c;
    slong * BBdegZ;
    slong AdegY, AdegX, AdegZ;
    slong bad_primes_left;

/*
flint_printf("+++++++++++++++++++++++\n");
flint_printf("n_polyu3_mod_factor_lift called: degree_inner = %wd, beta = %wd\n", degree_inner, beta);
flint_printf("A: "); n_polyu3_print_pretty(A, "Y", "X", "Z"); printf("\n");
for (i = 0; i < r; i++)
{
flint_printf("B[%wd]: ", i); n_polyu3_print_pretty(B + i, "Y", "X", "Z"); flint_printf("\n");
}
*/

    if (r <= 2)
        return n_polyu3_mod_factor_lift2(BB + 0, BB + 1, A, B + 0, B + 1, beta, degree_inner, ctx);

    FLINT_ASSERT(n_polyu_mod_is_canonical(A, ctx->mod));
    for (i = 0; i < r; i++)
        FLINT_ASSERT(n_polyu_mod_is_canonical(B + i, ctx->mod));

    BBdegZ = (slong *) flint_malloc(r*sizeof(slong));
    Bp = (n_bpoly_struct *) flint_malloc(r*sizeof(n_bpoly_struct));
    Bm = (n_bpoly_struct *) flint_malloc(r*sizeof(n_bpoly_struct));
    for (i = 0; i < r; i++)
    {
        n_bpoly_init(Bp + i);
        n_bpoly_init(Bm + i);
    }
    n_polyun_init(T);
    n_bpoly_init(Ap);
    n_bpoly_init(Am);
    n_poly_init(modulus);
    n_poly_init(alphapow);
    n_poly_init(t1);
    n_poly_init(t2);

    n_polyu3_degrees(&AdegY, &AdegX, &AdegZ, A);
    if (AdegX != degree_inner)
    {
        success = -1;
        goto cleanup;
    }

    n_poly_fit_length(alphapow, FLINT_MAX(WORD(3), AdegZ + 2));
    n_poly_one(modulus);

    FLINT_ASSERT((ctx->mod.n & UWORD(1)) == UWORD(1));

    alpha = (ctx->mod.n - UWORD(1))/UWORD(2);

    bad_primes_left = FLINT_MAX(5, AdegZ);

    goto choose_prime;

choose_prime:

    if (alpha < 2)
    {
        success = -1;
        goto cleanup;
    }
    alpha--;
/*
flint_printf("------ alpha: %wu ------\n", alpha);
*/
    FLINT_ASSERT(0 < alpha && alpha <= ctx->mod.n/2);
    FLINT_ASSERT(alphapow->alloc >= 2);
    alphapow->length = 2;
    alphapow->coeffs[0] = 1;
    alphapow->coeffs[1] = alpha;

    n_polyu3_mod_interp_reduce_2sm_bpoly(Ap, Am, A, alphapow, ctx->mod);
/*flint_printf(" Ap: "); n_bpoly_print_pretty(Ap, "y", "x"); flint_printf("\n");*/

    for (i = 0; i < r; i++)
    {
        n_polyu3_mod_interp_reduce_2sm_bpoly(Bp + i, Bm + i,
                                                    B + i, alphapow, ctx->mod);
/*
flint_printf("Bp[%wd]: ", i);
n_bpoly_print_pretty(Bp + i, "Y", "X");
flint_printf("\n");
fflush(stdout);
*/
    }

    success = n_bpoly_mod_hensel_lift(r, Ap, Bp, beta, degree_inner, ctx->mod);
    if (success <= 0)
    {
        if (success == 0 || --bad_primes_left < 0)
            goto cleanup;
        goto choose_prime;
    }

    success = n_bpoly_mod_hensel_lift(r, Am, Bm, beta, degree_inner, ctx->mod);
    if (success <= 0)
    {
        if (success == 0 || --bad_primes_left < 0)
            goto cleanup;
        goto choose_prime;
    }
/*
for (i = 0; i < r; i++)
{
flint_printf("lifted Bp[%wd]: ", i);
n_bpoly_print_pretty(Bp + i, "y", "x");
flint_printf("\n");
}
*/
    if (n_poly_degree(modulus) > 0)
    {
        c = n_poly_mod_evaluate_nmod(modulus, alpha, ctx->mod);
        FLINT_ASSERT(c == n_poly_mod_evaluate_nmod(modulus,
                                                ctx->mod.n - alpha, ctx->mod));
        c = nmod_mul(c, alpha, ctx->mod);
        c = nmod_add(c, c, ctx->mod);
        c = nmod_inv(c, ctx->mod);
        n_poly_mod_scalar_mul_nmod(modulus, modulus, c, ctx->mod);
        for (i = 0; i < r; i++)
        {
            n_polyu3n_mod_interp_crt_2sm_bpoly(BBdegZ + i, BB + i, T,
                                  Bp + i, Bm + i, modulus, alphapow, ctx->mod);
        }
    }
    else
    {
        for (i = 0; i < r; i++)
        {
            n_polyu3n_mod_interp_lift_2sm_bpoly(BBdegZ + i, BB + i,
                                               Bp + i, Bm + i, alpha, ctx->mod);
        }
    }

    c = ctx->mod.n - nmod_mul(alpha, alpha, ctx->mod);
    n_poly_mod_shift_left_scalar_addmul(modulus, 2, c, ctx->mod);
/*
flint_printf("modulus: "); n_poly_print_pretty(modulus, "x"); flint_printf("\n");
for (i = 0; i < r; i++)
{
flint_printf("BB[%wd]: ", i);
n_polyu3n_print_pretty(BB + i, "Y", "X", "?", "Z");
flint_printf("\n");
}
*/
    j = BBdegZ[0];
    for (i = 1; i < r; i++)
        j += BBdegZ[i];

    if (j > AdegZ)
    {
        success = 0;
        goto cleanup;
    }

    if (n_poly_degree(modulus) <= AdegZ)
    {
        goto choose_prime;
    }

    success = 1;

cleanup:
/*
flint_printf("+++++++++++++++++++++++\n");
flint_printf("n_polyu3_mod_factor_lift returning %d\n", success);
for (i = 0; i < r; i++)
{
flint_printf("BB[%wd]: ", i);
n_polyu3n_print_pretty(BB + i, "Y", "X", "?", "Z");
printf("\n");
}
fflush(stdout);
*/
#if WANT_ASSERT
    if (success == 1)
    {
        n_polyu_t T1, T2, T3;
        n_polyu_init(T1);
        n_polyu_init(T2);
        n_polyu_init(T3);
        n_polyu_get_n_polyun(T2, BB + 0);
        n_polyu_get_n_polyun(T3, BB + 1);
        n_polyu_mod_mul(T1, T2, T3, ctx->mod);
        for (i = 2; i < r; i++)
        {
            n_polyu_get_n_polyun(T3, BB + i);
            n_polyu_mod_mul(T2, T1, T3, ctx->mod);
            n_polyu_swap(T2, T1);
        }
        FLINT_ASSERT(n_polyu_equal(A, T1));
        n_polyu_clear(T1);
        n_polyu_clear(T2);
        n_polyu_clear(T3);
    }
#endif

    n_polyun_clear(T);
    n_bpoly_clear(Ap);
    n_bpoly_clear(Am);

    for (i = 0; i < r; i++)
    {
        n_bpoly_clear(Bp + i);
        n_bpoly_clear(Bm + i);
    }

    flint_free(BBdegZ);
    flint_free(Bp);
    flint_free(Bm);

    n_poly_clear(modulus);
    n_poly_clear(alphapow);
    n_poly_clear(t1);
    n_poly_clear(t2);

    return success;
}






/* r factor version */
int fmpz_mod_tpoly_hensel_lift(
    slong r,
    fmpz_mod_mpolyn_struct * BB,
    fmpz_mod_mpoly_t A,
    fmpz_mod_mpoly_struct * B,
    fmpz_t beta,
    slong degree_inner, /* required degree in x */
    fmpz_mod_mpoly_ctx_t ctx)
{
    FLINT_ASSERT(0 && "fmpz_mod_tpoly_hensel_lift not implemented");
    return -1;

#if 0
    int success;
    slong i, j;
    nmod_mpolyn_t T;
    n_bpoly_struct * Bp, * Bm;
    n_bpoly_t Ap, Am;
    nmod_poly_t modulus, modulus2, alphapow, t1, t2;
    mp_limb_t alpha, temp;
    slong * BBdegz;
    slong degrees[3], Adegz, Adegx;
    slong bad_primes_left;
    const char * vars[] = {"y", "x", "z"};

    if (r <= 2)
        return nmod_tpoly_hensel_lift2(BB + 0, BB + 1, A, B + 0, B + 1, beta, degree_inner, ctx);

    FLINT_ASSERT(ctx->minfo->nvars == 3);
    FLINT_ASSERT(A->bits == FLINT_BITS/3);
    for (i = 0; i < r; i++)
    {
        FLINT_ASSERT(B[i].bits == FLINT_BITS/3);
        FLINT_ASSERT(BB[i].bits == FLINT_BITS/3);
        FLINT_ASSERT(nmod_mpoly_is_canonical(B + i, ctx));
    }
    FLINT_ASSERT(nmod_mpoly_is_canonical(A, ctx));

flint_printf("+++++++++++++++++++++++\n");
flint_printf("nmod_tpoly_hensel_lift called: degree_inner = %wd, beta = %wd\n", degree_inner, beta);
flint_printf("A: "); nmod_mpoly_print_pretty(A, vars, ctx); printf("\n");
for (i = 0; i < r; i++)
{
flint_printf("B[%wd]: ", i); nmod_mpoly_print_pretty(B + i, vars, ctx); printf("\n");
}

    BBdegz = (slong *) flint_malloc(r*sizeof(slong));
    Bp = (n_bpoly_struct *) flint_malloc(r*sizeof(n_bpoly_struct));
    Bm = (n_bpoly_struct *) flint_malloc(r*sizeof(n_bpoly_struct));
    for (i = 0; i < r; i++)
    {
        n_bpoly_init(Bp + i, ctx->ffinfo->mod);
        n_bpoly_init(Bm + i, ctx->ffinfo->mod);
    }
    nmod_mpolyn_init(T, FLINT_BITS/3, ctx);
    n_bpoly_init(Ap, ctx->ffinfo->mod);
    n_bpoly_init(Am, ctx->ffinfo->mod);
    nmod_poly_init_mod(modulus, ctx->ffinfo->mod);
    nmod_poly_init_mod(modulus2, ctx->ffinfo->mod);
    nmod_poly_init_mod(alphapow, ctx->ffinfo->mod);
    nmod_poly_init_mod(t1, ctx->ffinfo->mod);
    nmod_poly_init_mod(t2, ctx->ffinfo->mod);

    nmod_mpoly_degrees_si(degrees, A, ctx);
    Adegx = degrees[1];
    Adegz = degrees[2];

    if (Adegx != degree_inner)
    {
        success = -1;
        goto cleanup;
    }

    nmod_poly_fit_length(alphapow, FLINT_MAX(WORD(3), Adegz + 2));
    nmod_poly_one(modulus);

    FLINT_ASSERT((ctx->ffinfo->mod.n & UWORD(1)) == UWORD(1));

    alpha = (ctx->ffinfo->mod.n - UWORD(1))/UWORD(2);

    bad_primes_left = FLINT_MAX(5, Adegz);

    goto choose_prime;

choose_prime: /* primes are v - alpha, v + alpha */

    if (alpha < 2)
    {
        success = -1;
        goto cleanup;
    }
    alpha--;

flint_printf("------ alpha: %wu ------\n", alpha);

    FLINT_ASSERT(0 < alpha && alpha <= ctx->ffinfo->mod.n/2);
    FLINT_ASSERT(alphapow->alloc >= 2);
    alphapow->length = 2;
    alphapow->coeffs[0] = 1;
    alphapow->coeffs[1] = alpha;

    nmod_tpoly_interp_reduce_2sm_bpoly(Ap, Am, A, alphapow, ctx);
flint_printf(" Ap: "); n_bpoly_print(Ap, "y", "x"); printf("\n");

    for (i = 0; i < r; i++)
    {
        nmod_tpoly_interp_reduce_2sm_bpoly(Bp + i, Bm + i, B + i, alphapow, ctx);
flint_printf("Bp[%wd]: ", i); n_bpoly_print(Bp + i, "y", "x"); printf("\n");
    }

    success = n_bpoly_hensel_lift(r, Ap, Bp, beta, degree_inner, ctx->ffinfo->mod);
    if (success <= 0)
    {
        if (success == 0 || --bad_primes_left < 0)
            goto cleanup;
        goto choose_prime;
    }

    success = n_bpoly_hensel_lift(r, Am, Bm, beta, degree_inner, ctx->ffinfo->mod);
    if (success <= 0)
    {
        if (success == 0 || --bad_primes_left < 0)
            goto cleanup;
        goto choose_prime;
    }

for (i = 0; i < r; i++)
{
flint_printf("lifted Bp[%wd]: ", i); n_bpoly_print(Bp + i, "y", "x"); printf("\n");
}

    if (nmod_poly_degree(modulus) > 0)
    {
        temp = nmod_poly_evaluate_nmod(modulus, alpha);
        FLINT_ASSERT(temp == nmod_poly_evaluate_nmod(modulus, ctx->ffinfo->mod.n - alpha));
        temp = nmod_mul(temp, alpha, ctx->ffinfo->mod);
        temp = nmod_add(temp, temp, ctx->ffinfo->mod);
        temp = n_invmod(temp, ctx->ffinfo->mod.n);
        nmod_poly_scalar_mul_nmod(modulus, modulus, temp);
        for (i = 0; i < r; i++)
            nmod_tpolyn_interp_crt_2sm_bpoly(BBdegz + i, BB + i, T, Bp + i, Bm + i, modulus, alphapow, ctx);
    }
    else
    {
        for (i = 0; i < r; i++)
            nmod_tpolyn_interp_lift_2sm_bpoly(BBdegz + i, BB + i, Bp + i, Bm + i, alpha, ctx);
    }
    temp = nmod_mul(alpha, alpha, ctx->ffinfo->mod);
    nmod_poly_scalar_mul_nmod(modulus2, modulus, temp);
    nmod_poly_shift_left(modulus, modulus, 2);
    nmod_poly_sub(modulus, modulus, modulus2);

printf("modulus: "); nmod_poly_print_pretty(modulus, "x"); printf("\n");
for (i = 0; i < r; i++)
{
flint_printf("BB[%wd]: ", i); nmod_mpolyn_print_pretty(BB + i, vars, ctx); printf("\n");
}
fflush(stdout);

    j = BBdegz[0];
    for (i = 1; i < r; i++)
        j += BBdegz[i];

    if (j > Adegz)
    {
        success = 0;
        goto cleanup;
    }

    if (nmod_poly_degree(modulus) <= Adegz)
    {
        goto choose_prime;
    }

    success = 1;

cleanup:

flint_printf("+++++++++++++++++++++++\n");
flint_printf("nmod_tpoly_hensel_lift returning %d\n", success);
for (i = 0; i < r; i++)
{
flint_printf("BB[%wd]: ", i); nmod_mpolyn_print_pretty(BB + i, vars, ctx); printf("\n");
}
fflush(stdout);

    if (success == 1)
    {
        slong perm[3] = {0, 1, 2};
        ulong shift[3] = {0, 0, 0};
        ulong stride[3] = {1, 1, 1};
        nmod_mpoly_t tBB0, tBB1, tP;

        nmod_mpoly_init(tBB0, ctx);
        nmod_mpoly_init(tBB1, ctx);
        nmod_mpoly_init(tP, ctx);

        nmod_mpoly_from_mpolyn_perm_inflate(tBB0, FLINT_BITS/3, ctx, BB + 0, ctx, perm, shift, stride);
        nmod_mpoly_from_mpolyn_perm_inflate(tBB1, FLINT_BITS/3, ctx, BB + 1, ctx, perm, shift, stride);
        nmod_mpoly_mul(tP, tBB0, tBB1, ctx);
        for (i = 2; i < r; i++)
        {
            nmod_mpoly_from_mpolyn_perm_inflate(tBB0, FLINT_BITS/3, ctx, BB + i, ctx, perm, shift, stride);
            nmod_mpoly_mul(tBB1, tP, tBB0, ctx);
            nmod_mpoly_swap(tP, tBB1, ctx);
        }

        FLINT_ASSERT(nmod_mpoly_equal(tP, A, ctx));

        nmod_mpoly_clear(tBB0, ctx);
        nmod_mpoly_clear(tBB1, ctx);
        nmod_mpoly_clear(tP, ctx);
    }

    nmod_mpolyn_clear(T, ctx);
    n_bpoly_clear(Ap);
    n_bpoly_clear(Am);

    for (i = 0; i < r; i++)
    {
        n_bpoly_clear(Bp + i);
        n_bpoly_clear(Bm + i);
    }

    flint_free(BBdegz);
    flint_free(Bp);
    flint_free(Bm);

    nmod_poly_clear(modulus);
    nmod_poly_clear(modulus2);
    nmod_poly_clear(alphapow);
    nmod_poly_clear(t1);
    nmod_poly_clear(t2);

    return success;
#endif
}




fmpz_mpoly_struct * _fmpz_mpolyu_get_coeff(
    fmpz_mpolyu_t A,
    ulong pow,
    const fmpz_mpoly_ctx_t uctx);


/*
    B vars: x0 x1 x2 x3 x4 xv           2 < v
    A vars: xv x0 x1 : 0 0 x2 x3 x4 0
*/
static void fmpz_mpoly_to_fmpz_mpolyuuu(
    slong v,
    fmpz_mpolyu_t A,
    const fmpz_mpoly_ctx_t Actx,
    const fmpz_mpoly_t B,
    const fmpz_mpoly_ctx_t Bctx)
{
    slong i, j;
    fmpz_mpoly_struct * Ac;
    ulong * Bexps;
    slong NA, NB;


    FLINT_ASSERT(Bctx->minfo->nvars > v);
    FLINT_ASSERT(Actx->minfo->nvars == v - 2);

    Bexps = (ulong *) flint_malloc(Bctx->minfo->nvars*sizeof(ulong));

    NA = mpoly_words_per_exp_sp(A->bits, Actx->minfo);
    NB = mpoly_words_per_exp(B->bits, Bctx->minfo);

    A->length = 0;
    for (j = 0; j < B->length; j++)
    {
        mpoly_get_monomial_ui(Bexps, B->exps + NB*j, B->bits, Bctx->minfo);
        FLINT_ASSERT(FLINT_BIT_COUNT(Bexps[0]) < FLINT_BITS/3);
        FLINT_ASSERT(FLINT_BIT_COUNT(Bexps[1]) < FLINT_BITS/3);
        FLINT_ASSERT(FLINT_BIT_COUNT(Bexps[v]) < FLINT_BITS/3);
        Ac = _fmpz_mpolyu_get_coeff(A, (Bexps[v] << (2*(FLINT_BITS/3))) +
                                       (Bexps[0] << (1*(FLINT_BITS/3))) +
                                       (Bexps[1] << (0*(FLINT_BITS/3))), Actx);
        FLINT_ASSERT(Ac->bits == A->bits);
        fmpz_mpoly_fit_length(Ac, Ac->length + 1, Actx);
        fmpz_set(Ac->coeffs + Ac->length, B->coeffs + j);
        mpoly_set_monomial_ui(Ac->exps + NA*Ac->length, Bexps + 2, A->bits, Actx->minfo);
        Ac->length++;
    }

    for (i = 0; i < A->length; i++)
    {
        fmpz_mpoly_sort_terms(A->coeffs + i, Actx);
    }

    flint_free(Bexps);
}

void fmpz_mpoly_from_fmpz_mpolyuuu(
    slong v,
    fmpz_mpoly_t B,
    flint_bitcnt_t Bbits,
    const fmpz_mpoly_ctx_t Bctx,
    const fmpz_mpolyu_t A,
    const fmpz_mpoly_ctx_t Actx)
{
    slong i, j;
    ulong * exps;
    slong NA, NB;

    FLINT_ASSERT(Bctx->minfo->nvars > v);
    FLINT_ASSERT(Actx->minfo->nvars == v - 2);

    exps = (ulong *) flint_malloc(Bctx->minfo->nvars*sizeof(ulong));
    for (i = 0; i < Bctx->minfo->nvars; i++)
        exps[i] = 0;

    B->length = 0;
    fmpz_mpoly_fit_bits(B, Bbits, Bctx);
    B->bits = Bbits;

    NA = mpoly_words_per_exp_sp(A->bits, Actx->minfo);
    NB = mpoly_words_per_exp(B->bits, Bctx->minfo);

    for (i = 0; i < A->length; i++)
    {
        fmpz_mpoly_struct * Ac = A->coeffs + i;
        FLINT_ASSERT(Ac->bits == A->bits);
        fmpz_mpoly_fit_length(B, B->length + Ac->length, Bctx);
        for (j = 0; j < Ac->length; j++)
        {
            if (fmpz_is_zero(Ac->coeffs + j))
                continue;
            mpoly_get_monomial_ui(exps + 2, Ac->exps + NA*j, A->bits, Actx->minfo);
            exps[0] = (A->exps[i] >> (1*(FLINT_BITS/3))) & LOW_THIRD_MASK;
            exps[1] = (A->exps[i] >> (0*(FLINT_BITS/3))) & LOW_THIRD_MASK;
            exps[v] = (A->exps[i] >> (2*(FLINT_BITS/3))) & LOW_THIRD_MASK;
            fmpz_set(B->coeffs + B->length, Ac->coeffs + j);
            mpoly_set_monomial_ui(B->exps + NB*B->length, exps, B->bits, Bctx->minfo);
            B->length++;
        }
    }

    flint_free(exps);

    fmpz_mpoly_sort_terms(B, Bctx);
}


static void fmpz_mpolyuuu_print_pretty(
    const fmpz_mpolyu_t A,
    const char * var0,
    const char * var1,
    const char * var2,
    const char ** vars,
    const fmpz_mpoly_ctx_t ctx)
{
    slong i;
    int first = 1;

    for (i = 0; i < A->length; i++)
    {
        if (!first)
            flint_printf(" + ");
        first = 0;
        flint_printf("(");
        fmpz_mpoly_print_pretty(A->coeffs + i, vars, ctx);
        flint_printf(")*%s^%wu*%s^%wu*%s^%wu",
            var0, (A->exps[i] >> (2*(FLINT_BITS/3))) & LOW_THIRD_MASK,
            var1, (A->exps[i] >> (1*(FLINT_BITS/3))) & LOW_THIRD_MASK,
            var2, (A->exps[i] >> (0*(FLINT_BITS/3))) & LOW_THIRD_MASK);
    }

    if (first)
        flint_printf("0");
}

static void nmod_mpolyuuu_print_pretty(
    const nmod_mpolyu_t A,
    const char * var0,
    const char * var1,
    const char * var2,
    const char ** vars,
    const nmod_mpoly_ctx_t ctx)
{
    slong i;
    int first = 1;

    for (i = 0; i < A->length; i++)
    {
        if (!first)
            flint_printf(" + ");
        first = 0;
        flint_printf("(");
        nmod_mpoly_print_pretty(A->coeffs + i, vars, ctx);
        flint_printf(")*%s^%wu*%s^%wu*%s^%wu",
            var0, (A->exps[i] >> (2*(FLINT_BITS/3))) & LOW_THIRD_MASK,
            var1, (A->exps[i] >> (1*(FLINT_BITS/3))) & LOW_THIRD_MASK,
            var2, (A->exps[i] >> (0*(FLINT_BITS/3))) & LOW_THIRD_MASK);
    }

    if (first)
        flint_printf("0");
}



void nmod_tpoly_use_skel_mul(
    nmod_mpoly_t E,
    const fmpz_mpolyu_t A,
    const nmod_mpolycu_t Ared,
    nmod_mpolycu_t Acur,
    const nmod_mpolycu_t Ainc,
    const nmod_mpoly_ctx_t ctx_sp)
{
/*    slong xexp, yexp;*/
    slong i;
    mp_limb_t eval;
/*
flint_printf("nmod_tpoly_use_skel_mul called\n"); fflush(stdout);
*/

    FLINT_ASSERT(ctx_sp->minfo->nvars == 3);
    FLINT_ASSERT(E->bits == FLINT_BITS/3);
    FLINT_ASSERT(1 == mpoly_words_per_exp_sp(E->bits, ctx_sp->minfo));

    FLINT_ASSERT(A->length == Ared->length);
    FLINT_ASSERT(A->length == Acur->length);
    FLINT_ASSERT(A->length == Ainc->length);

    E->length = 0;
    for (i = 0; i < A->length; i++)
    {
/*
flint_printf("nmod_tpoly_use_skel_mul i = %wd\n", i); fflush(stdout);
*/
        eval = nmod_mpoly_use_skel_mul(Ared->coeffs + i, Acur->coeffs + i,
                                             Ainc->coeffs + i, ctx_sp->ffinfo);
        if (eval == 0)
        {
            continue;
        }

        nmod_mpoly_fit_length(E, E->length + 1, ctx_sp);
        E->exps[E->length] = A->exps[i];
        E->coeffs[E->length] = eval;
        E->length++;
    }
/*
flint_printf("nmod_tpoly_use_skel_mul returning\n"); fflush(stdout);
*/
}

void fmpz_mod_tpoly_use_skel_mul(
    fmpz_mod_mpoly_t E,
    const fmpz_mpolyu_t A,
    const fmpz_mpolycu_t Ared,
    fmpz_mpolycu_t Acur,
    const fmpz_mpolycu_t Ainc,
    const fmpz_mod_mpoly_ctx_t ctx_mp)
{
/*    slong xexp, yexp;*/
    slong i;
    fmpz_t eval;

    FLINT_ASSERT(ctx_mp->minfo->nvars == 3);
    FLINT_ASSERT(E->bits == FLINT_BITS/3);
    FLINT_ASSERT(1 == mpoly_words_per_exp_sp(E->bits, ctx_mp->minfo));

    FLINT_ASSERT(A->length == Ared->length);
    FLINT_ASSERT(A->length == Acur->length);
    FLINT_ASSERT(A->length == Ainc->length);

    fmpz_init(eval);

    E->length = 0;
    for (i = 0; i < A->length; i++)
    {
        fmpz_mod_mpoly_use_skel_mul(eval, Ared->coeffs + i, Acur->coeffs + i,
                                             Ainc->coeffs + i, ctx_mp->ffinfo);
        if (fmpz_is_zero(eval))
        {
            continue;
        }

        fmpz_mod_mpoly_fit_length(E, E->length + 1, ctx_mp);
        E->exps[E->length] = A->exps[i];
        fmpz_set(E->coeffs + E->length, eval);
        E->length++;
    }

    fmpz_clear(eval);
}


void nmod_bma_mpoly_add_point_tpoly(
    nmod_bma_mpoly_t L,
    const nmod_mpolyn_t A,
    const nmod_mpoly_ctx_t ctx_sp)
{
    slong j;
    slong Alen = A->length;
    nmod_poly_struct * Acoeff = A->coeffs;
    slong Li, Ai, ai;
    ulong Aexp;
    nmod_berlekamp_massey_struct * Lcoeff;
    slong Llen;
    ulong * Lexp;

    if (L->length == 0)
    {
        slong tot = 0;
        for (Ai = 0; Ai < Alen; Ai++)
        {
            for (ai = (Acoeff + Ai)->length - 1; ai >= 0; ai--)
            {
                tot += (0 != (Acoeff + Ai)->coeffs[ai]);
            }
        }
        nmod_bma_mpoly_fit_length(L, tot, ctx_sp->ffinfo);
    }

    Lcoeff = L->coeffs;
    Llen = L->length;
    Lexp = L->exps;

    Li = 0;
    Ai = ai = 0;
    Aexp = 0;
    if (Ai < Alen)
    {
        ai = nmod_poly_degree(A->coeffs + Ai);
        FLINT_ASSERT(FLINT_BIT_COUNT(ai) < FLINT_BITS/3);
        FLINT_ASSERT(0 == (A->exps[Ai] & LOW_THIRD_MASK));
        Aexp = A->exps[Ai] + ai;
    }

    while (Li < Llen || Ai < Alen)
    {
        if (Li < Llen && Ai < Alen && Lexp[Li] == Aexp)
        {
            /* L term present, A term present */
add_same_exp:
            nmod_berlekamp_massey_add_point(Lcoeff + Li, (Acoeff + Ai)->coeffs[ai]);
            Li++;
            do {
                ai--;
            } while (ai >= 0 && (0 == (Acoeff + Ai)->coeffs[ai]));
            if (ai < 0)
            {
                Ai++;
                if (Ai < Alen)
                {
                    ai = nmod_poly_degree(A->coeffs + Ai);
                    FLINT_ASSERT(FLINT_BIT_COUNT(ai) < FLINT_BITS/3);
                    FLINT_ASSERT(0 == (A->exps[Ai] & LOW_THIRD_MASK));
                    Aexp = A->exps[Ai] + ai;
                }
            }
            else
            {
                FLINT_ASSERT(Ai < Alen);
                FLINT_ASSERT(FLINT_BIT_COUNT(ai) < FLINT_BITS/3);
                FLINT_ASSERT(0 == (A->exps[Ai] & LOW_THIRD_MASK));
                Aexp = A->exps[Ai] + ai;
            }
        }
        else if (Li < Llen && (Ai >= Alen || Lexp[Li] > Aexp))
        {
            /* L term present, A term missing */
            nmod_berlekamp_massey_add_zeros(Lcoeff + Li, 1);
            Li++;
        }
        else
        {
            /* L term missing, A term present */
            FLINT_ASSERT(Ai < Alen && (Li >= Llen || Lexp[Li] < Aexp));
            {
                ulong texp;
                nmod_berlekamp_massey_struct tcoeff;

                nmod_bma_mpoly_fit_length(L, Llen + 1, ctx_sp->ffinfo);
                Lcoeff = L->coeffs;
                Lexp = L->exps;

                texp = Lexp[Llen];
                tcoeff = Lcoeff[Llen];
                for (j = Llen - 1; j >= Li; j--)
                {
                    Lexp[j + 1] = Lexp[j];
                    Lcoeff[j + 1] = Lcoeff[j];
                }
                Lexp[Li] = texp;
                Lcoeff[Li] = tcoeff;
            }

            nmod_berlekamp_massey_start_over(Lcoeff + Li);
            nmod_berlekamp_massey_add_zeros(Lcoeff + Li, L->pointcount);
            Lexp[Li] = Aexp;
            Llen++;
            L->length = Llen;

            goto add_same_exp;
        }
    }

    L->pointcount++;
}

void fmpz_mod_bma_mpoly_add_point_tpoly(
    fmpz_mod_bma_mpoly_t L,
    const fmpz_mod_mpolyn_t A,
    const fmpz_mod_mpoly_ctx_t ctx_mp)
{
    slong j;
    slong Alen = A->length;
    fmpz_mod_poly_struct * Acoeff = A->coeffs;
    slong Li, Ai, ai;
    ulong Aexp;
    fmpz_mod_berlekamp_massey_struct * Lcoeff;
    slong Llen;
    ulong * Lexp;

    if (L->length == 0)
    {
        slong tot = 0;
        for (Ai = 0; Ai < Alen; Ai++)
        {
            for (ai = (Acoeff + Ai)->length - 1; ai >= 0; ai--)
            {
                tot += !fmpz_is_zero((Acoeff + Ai)->coeffs + ai);
            }
        }
        fmpz_mod_bma_mpoly_fit_length(L, tot, ctx_mp->ffinfo);
    }

    Lcoeff = L->coeffs;
    Llen = L->length;
    Lexp = L->exps;

    Li = 0;
    Ai = ai = 0;
    Aexp = 0;
    if (Ai < Alen)
    {
        ai = fmpz_mod_poly_degree(A->coeffs + Ai);
        FLINT_ASSERT(FLINT_BIT_COUNT(ai) < FLINT_BITS/3);
        FLINT_ASSERT(0 == (A->exps[Ai] & LOW_THIRD_MASK));
        Aexp = A->exps[Ai] + ai;
    }

    while (Li < Llen || Ai < Alen)
    {
        if (Li < Llen && Ai < Alen && Lexp[Li] == Aexp)
        {
            /* L term present, A term present */
add_same_exp:
            fmpz_mod_berlekamp_massey_add_point(Lcoeff + Li,
                                                   (Acoeff + Ai)->coeffs + ai);
            Li++;
            do {
                ai--;
            } while (ai >= 0 && fmpz_is_zero((Acoeff + Ai)->coeffs + ai));
            if (ai < 0)
            {
                Ai++;
                if (Ai < Alen)
                {
                    ai = fmpz_mod_poly_degree(A->coeffs + Ai);
                    FLINT_ASSERT(FLINT_BIT_COUNT(ai) < FLINT_BITS/3);
                    FLINT_ASSERT(0 == (A->exps[Ai] & LOW_THIRD_MASK));
                    Aexp = A->exps[Ai] + ai;
                }
            }
            else
            {
                FLINT_ASSERT(Ai < Alen);
                FLINT_ASSERT(FLINT_BIT_COUNT(ai) < FLINT_BITS/3);
                FLINT_ASSERT(0 == (A->exps[Ai] & LOW_THIRD_MASK));
                Aexp = A->exps[Ai] + ai;
            }
        }
        else if (Li < Llen && (Ai >= Alen || Lexp[Li] > Aexp))
        {
            /* L term present, A term missing */
            fmpz_mod_berlekamp_massey_add_zeros(Lcoeff + Li, 1);
            Li++;
        }
        else
        {
            /* L term missing, A term present */
            FLINT_ASSERT(Ai < Alen && (Li >= Llen || Lexp[Li] < Aexp));
            {
                ulong texp;
                fmpz_mod_berlekamp_massey_struct tcoeff;

                fmpz_mod_bma_mpoly_fit_length(L, Llen + 1, ctx_mp->ffinfo);
                Lcoeff = L->coeffs;
                Lexp = L->exps;

                texp = Lexp[Llen];
                tcoeff = Lcoeff[Llen];
                for (j = Llen - 1; j >= Li; j--)
                {
                    Lexp[j + 1] = Lexp[j];
                    Lcoeff[j + 1] = Lcoeff[j];
                }
                Lexp[Li] = texp;
                Lcoeff[Li] = tcoeff;
            }

            fmpz_mod_berlekamp_massey_start_over(Lcoeff + Li);
            fmpz_mod_berlekamp_massey_add_zeros(Lcoeff + Li, L->pointcount);
            Lexp[Li] = Aexp;
            Llen++;
            L->length = Llen;
            goto add_same_exp;
        }
    }

    L->pointcount++;
}



int fmpz_mpolyu_equal(
    const fmpz_mpolyu_t A,
    const fmpz_mpolyu_t B,
    const fmpz_mpoly_ctx_t ctx)
{
    slong i;

    if (A->length != B->length)
        return 0;

    for (i = 0; i < A->length; i++)
    {
        if (A->exps[i] != B->exps[i])
            return 0;
        if (!fmpz_mpoly_equal(A->coeffs + i, B->coeffs + i, ctx))
            return 0;
    }

    return 1;
}


int nmod_zip_mpolyuuu_add_point(
    nmod_zip_mpolyu_t L,
    const nmod_mpolyn_t A)
{
    slong Alen = A->length;
    nmod_poly_struct * Acoeff = A->coeffs;
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
        ai = nmod_poly_degree(A->coeffs + Ai);
        FLINT_ASSERT(FLINT_BIT_COUNT(ai) < FLINT_BITS/3);
        FLINT_ASSERT(0 == (A->exps[Ai] & LOW_THIRD_MASK));
        Aexp = A->exps[Ai] + ai;
    }

    for (Li = 0; Li < Llen; Li++)
    {
        nmod_zip_struct * Lc = Lcoeff + Li;

        if (Ai < Alen && Lexp[Li] == Aexp)
        {
            /* L present A present */
            FLINT_ASSERT(pointcount <= Lc->ealloc);
            Lc->evals[pointcount] = (Acoeff + Ai)->coeffs[ai];
            do {
                ai--;
            } while (ai >= 0 && (0 == (Acoeff + Ai)->coeffs[ai]));
            if (ai < 0)
            {
                Ai++;
                if (Ai < Alen)
                {
                    ai = nmod_poly_degree(A->coeffs + Ai);        
                    FLINT_ASSERT(FLINT_BIT_COUNT(ai) < FLINT_BITS/3);
                    FLINT_ASSERT(0 == (A->exps[Ai] & LOW_THIRD_MASK));
                    Aexp = A->exps[Ai] + ai;
                }
            }
            else
            {
                FLINT_ASSERT(Ai < Alen);
                FLINT_ASSERT(FLINT_BIT_COUNT(ai) < FLINT_BITS/3);
                FLINT_ASSERT(0 == (A->exps[Ai] & LOW_THIRD_MASK));
                Aexp = A->exps[Ai] + ai;
            }
        }
        else if (Ai >= Alen || Lexp[Li] > Aexp)
        {
            /* L present A missing */
            FLINT_ASSERT(pointcount <= Lc->ealloc);
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



slong zip_form_new(
    ulong * Xdeg_,
    nmod_mpolyu_t M,
    const nmod_mpoly_ctx_t nctx,
    const fmpz_mpolyu_t B,
    const fmpz_mpoly_ctx_t ctx)
{
    ulong Xdeg;
    slong N = mpoly_words_per_exp(B->bits, ctx->minfo);
    slong i, j, max_length;
    const char * vars[] = {"x", "y", "z", "w", "?", "?", "?"};
    slong v = ctx->minfo->nvars + 2;

    FLINT_ASSERT(M->bits == B->bits);
    FLINT_ASSERT(N == mpoly_words_per_exp(M->bits, nctx->minfo));

    Xdeg = 0;
    for (i = 0; i < B->length; i++)
    {
        Xdeg = FLINT_MAX(Xdeg, (B->exps[i] >> (1*(FLINT_BITS/3))) & LOW_THIRD_MASK);
    }

    *Xdeg_ = Xdeg;

flint_printf("zip_form_new called Xdeg: %wu\n", Xdeg);
flint_printf("B: "); fmpz_mpolyuuu_print_pretty(B, vars[v], vars[0], vars[1], vars + 2, ctx); printf("\n");

    M->length = 0;
    max_length = 0;
    for (i = 0; i < B->length; i++)
    {
        fmpz_mpoly_struct * Bc = B->coeffs + i;
        nmod_mpoly_struct * Mc;
/*        ulong y = (B->exps[i] >> (2*(FLINT_BITS/3))) & LOW_THIRD_MASK;*/
        ulong x = (B->exps[i] >> (1*(FLINT_BITS/3))) & LOW_THIRD_MASK;
        ulong z = (B->exps[i] >> (0*(FLINT_BITS/3))) & LOW_THIRD_MASK;

        if (x >= Xdeg)
            continue;

        Mc = _nmod_mpolyu_get_coeff(M, (x << (1*(FLINT_BITS/3))) + (z << (0*(FLINT_BITS/3))), nctx);
        FLINT_ASSERT(Mc->length == 0);
        nmod_mpoly_fit_length(Mc, Bc->length, nctx);
        max_length = FLINT_MAX(max_length, Bc->length);
        Mc->length = Bc->length;
        for (j = 0; j < Bc->length; j++)
        {
            Mc->coeffs[j] = 1;
            mpoly_monomial_set(Mc->exps + N*j, Bc->exps + N*j, N);
        }
    }

flint_printf("M: "); nmod_mpolyuuu_print_pretty(M, vars[v], vars[0], vars[1], vars + 2, nctx); printf("\n");

flint_printf("returning max_length: %wd\n", max_length);

    return max_length;
}


/*
    return 
        -1: singular vandermonde matrix
        0:  inconsistent system
        1:  success
*/
int nmod_zip_find_coeffs_new(
    mp_limb_t * coeffs,             /* length mlength */
    const mp_limb_t * monomials,    /* length mlength */
    slong mlength,
    const mp_limb_t * evals,        /* length elength */
    slong elength,    
    nmod_poly_t master,
    const nmodf_ctx_t ffinfo)
{
    slong i, j;
    mp_limb_t V, V0, V1, V2, T, S, r, p0, p1;

    FLINT_ASSERT(elength >= mlength);

    nmod_poly_product_roots_nmod_vec(master, monomials, mlength);

    for (i = 0; i < mlength; i++)
    {
        /* coeffs[i] is (coeffs(P).values)/P(roots[i]) =: V/S
            where P(x) = master(x)/(x-roots[i])     */
        V0 = V1 = V2 = T = S = 0;
        r = monomials[i];
        for (j = mlength; j > 0; j--)
        {
            T = nmod_add(nmod_mul(r, T, ffinfo->mod), master->coeffs[j], ffinfo->mod);
            S = nmod_add(nmod_mul(r, S, ffinfo->mod), T, ffinfo->mod);
            umul_ppmm(p1, p0, evals[j - 1], T);
            add_sssaaaaaa(V2, V1, V0, V2, V1, V0, WORD(0), p1, p0);
        }
        /* roots[i] should be a root of master */
        FLINT_ASSERT(nmod_add(nmod_mul(r, T, ffinfo->mod), master->coeffs[0], ffinfo->mod) == 0);
        NMOD_RED3(V, V2, V1, V0, ffinfo->mod);
        S = nmod_mul(S, r, ffinfo->mod); /* shift is one */
        if (S == 0)
            return -1;
        coeffs[i] = nmod_mul(V, nmod_inv(S, ffinfo->mod), ffinfo->mod);
    }

    /* use the coefficients of master as temp work space */
    for (i = 0; i < mlength; i++)
    {
        master->coeffs[i] = nmod_pow_ui(monomials[i], mlength, ffinfo->mod);
    }

    /* check that the remaining points match */
    for (i = mlength; i < elength; i++)
    {
        V0 = V1 = V2 = S = 0;
        for (j = 0; j < mlength; j++)
        {
            master->coeffs[j] = nmod_mul(master->coeffs[j], monomials[j], ffinfo->mod);
            umul_ppmm(p1, p0, coeffs[j], master->coeffs[j]);
            add_sssaaaaaa(V2, V1, V0, V2, V1, V0, WORD(0), p1, p0);
        }
        NMOD_RED3(V, V2, V1, V0, ffinfo->mod);
        if (V != evals[i])
            return 0;
    }

    return 1;
}

int nmod_zip_find_coeffs_new2(
    mp_limb_t * coeffs,             /* length mlength */
    const mp_limb_t * monomials,    /* length mlength */
    slong mlength,
    const mp_limb_t * evals,        /* length elength */
    slong elength,
    const mp_limb_t * master,       /* length mlength + 1 */
    mp_limb_t * temp,               /* length mlength */
    nmod_t mod)
{
    slong i, j;
    mp_limb_t V, V0, V1, V2, T, S, r, p0, p1;

    FLINT_ASSERT(elength >= mlength);

    for (i = 0; i < mlength; i++)
    {
        /* coeffs[i] is (coeffs(P).values)/P(roots[i]) =: V/S
            where P(x) = master(x)/(x-roots[i])     */
        V0 = V1 = V2 = T = S = 0;
        r = monomials[i];
        for (j = mlength; j > 0; j--)
        {
            T = nmod_add(nmod_mul(r, T, mod), master[j], mod);
            S = nmod_add(nmod_mul(r, S, mod), T, mod);
            umul_ppmm(p1, p0, evals[j - 1], T);
            add_sssaaaaaa(V2, V1, V0, V2, V1, V0, 0, p1, p0);
        }
        /* roots[i] should be a root of master */
        FLINT_ASSERT(nmod_add(nmod_mul(r, T, mod), master[0], mod) == 0);
        NMOD_RED3(V, V2, V1, V0, mod);
        S = nmod_mul(S, r, mod); /* shift is one */
        if (S == 0)
            return -1;
        coeffs[i] = nmod_mul(V, nmod_inv(S, mod), mod);
    }

    /* check that the remaining points match */
    for (j = 0; j < mlength; j++)
        temp[j] = nmod_pow_ui(monomials[j], mlength, mod);

    for (i = mlength; i < elength; i++)
    {
        V0 = V1 = V2 = S = 0;
        for (j = 0; j < mlength; j++)
        {
            temp[j] = nmod_mul(temp[j], monomials[j], mod);
            umul_ppmm(p1, p0, coeffs[j], temp[j]);
            add_sssaaaaaa(V2, V1, V0, V2, V1, V0, 0, p1, p0);
        }
        NMOD_RED3(V, V2, V1, V0, mod);
        if (V != evals[i])
            return 0;
    }
    return 1;
}



void nmod_zip_fit_elength(nmod_zip_t A, slong elength)
{
    slong old_ealloc = A->ealloc;
    slong new_ealloc = FLINT_MAX(elength, A->ealloc + A->ealloc/2);

    FLINT_ASSERT(elength > 0);

    if (elength > old_ealloc)
    {
        if (old_ealloc == 0)
        {
            A->evals = (mp_limb_t *) flint_malloc(new_ealloc*sizeof(mp_limb_t));
        }
        else
        {
            A->evals = (mp_limb_t *) flint_realloc(A->evals,
                                                 new_ealloc*sizeof(mp_limb_t));
        }
        A->ealloc = new_ealloc;
    }
}

void nmod_zip_set_mlength(nmod_zip_t A, slong mlength)
{
    slong old_malloc = A->malloc;
    slong new_malloc = FLINT_MAX(mlength, A->malloc + A->malloc/2);

    FLINT_ASSERT(mlength > 0);

    if (mlength > old_malloc)
    {
        if (old_malloc == 0)
        {
            A->coeffs    = (mp_limb_t *) flint_malloc(
                                                 new_malloc*sizeof(mp_limb_t));
            A->monomials = (mp_limb_t *) flint_malloc(
                                                 new_malloc*sizeof(mp_limb_t));
        }
        else
        {
            A->coeffs    = (mp_limb_t *) flint_realloc(A->coeffs,
                                                 new_malloc*sizeof(mp_limb_t));
            A->monomials = (mp_limb_t *) flint_realloc(A->monomials,
                                                 new_malloc*sizeof(mp_limb_t));
        }
        A->malloc = new_malloc;
    }

    A->mlength = mlength;
}



void nmod_zip_mpolyu_fit_length(
    nmod_zip_mpolyu_t A,
    slong length);


void nmod_zip_mpolyuuu_add_point_allow_new(
    nmod_zip_mpolyu_t L,
    const nmod_mpolyn_t A)
{
    slong j;
    slong Alen = A->length;
    nmod_poly_struct * Acoeff = A->coeffs;
    slong Li, Ai, ai;
    ulong Aexp;
    nmod_zip_struct * Lcoeff;
    slong Llen;
    ulong * Lexp;
    slong pointcount = L->pointcount;

    if (L->length == 0)
    {
        slong tot = 0;
        for (Ai = 0; Ai < Alen; Ai++)
        {
            for (ai = (Acoeff + Ai)->length - 1; ai >= 0; ai--)
            {
                tot += (0 != (Acoeff + Ai)->coeffs[ai]);
            }
        }
        nmod_zip_mpolyu_fit_length(L, tot);
    }

    Lcoeff = L->coeffs;
    Llen = L->length;
    Lexp = L->exps;

    Li = 0;
    Ai = ai = 0;
    Aexp = 0;
    if (Ai < Alen)
    {
        ai = nmod_poly_degree(A->coeffs + Ai);
        FLINT_ASSERT(FLINT_BIT_COUNT(ai) < FLINT_BITS/3);
        FLINT_ASSERT(0 == (A->exps[Ai] & LOW_THIRD_MASK));
        Aexp = A->exps[Ai] + ai;
    }

    while (Li < Llen || Ai < Alen)
    {
        if (Li < Llen && Ai < Alen && Lexp[Li] == Aexp)
        {
            /* L term present, A term present */
add_same_exp:
            nmod_zip_fit_elength(Lcoeff + Li, pointcount + 1);
            (Lcoeff + Li)->evals[pointcount] = (Acoeff + Ai)->coeffs[ai];
            Li++;
            do {
                ai--;
            } while (ai >= 0 && (0 == (Acoeff + Ai)->coeffs[ai]));
            if (ai < 0)
            {
                Ai++;
                if (Ai < Alen)
                {
                    ai = nmod_poly_degree(A->coeffs + Ai);
                    FLINT_ASSERT(FLINT_BIT_COUNT(ai) < FLINT_BITS/3);
                    FLINT_ASSERT(0 == (A->exps[Ai] & LOW_THIRD_MASK));
                    Aexp = A->exps[Ai] + ai;
                }
            }
            else
            {
                FLINT_ASSERT(Ai < Alen);
                FLINT_ASSERT(FLINT_BIT_COUNT(ai) < FLINT_BITS/3);
                FLINT_ASSERT(0 == (A->exps[Ai] & LOW_THIRD_MASK));
                Aexp = A->exps[Ai] + ai;
            }
        }
        else if (Li < Llen && (Ai >= Alen || Lexp[Li] > Aexp))
        {
            /* L term present, A term missing */
            nmod_zip_fit_elength(Lcoeff + Li, pointcount + 1);
            (Lcoeff + Li)->evals[pointcount] = 0;
            Li++;
        }
        else
        {
            /* L term missing, A term present */
            FLINT_ASSERT(Ai < Alen && (Li >= Llen || Lexp[Li] < Aexp));
            {
                ulong texp;
                nmod_zip_struct tcoeff;

                nmod_zip_mpolyu_fit_length(L, Llen + 1);
                Lcoeff = L->coeffs;
                Lexp = L->exps;

                texp = Lexp[Llen];
                tcoeff = Lcoeff[Llen];
                for (j = Llen - 1; j >= Li; j--)
                {
                    Lexp[j + 1] = Lexp[j];
                    Lcoeff[j + 1] = Lcoeff[j];
                }
                Lexp[Li] = texp;
                Lcoeff[Li] = tcoeff;
            }

            nmod_zip_fit_elength(Lcoeff + Li, pointcount + 1);
            for (j = 0; j <= pointcount; j++)
                (Lcoeff + Li)->evals[j] = 0;

            Lexp[Li] = Aexp;
            Llen++;
            L->length = Llen;

            goto add_same_exp;
        }
    }

    L->pointcount = pointcount + 1;
}

void nmod_zip_mpolyuuu_print(const nmod_zip_mpolyu_t A)
{
    slong i;
    flint_printf("0");
    for (i = 0; i < A->length; i++)
    {
        flint_printf(" + [");
        nmod_zip_print(A->coeffs + i, A->pointcount);
        flint_printf("]*Y^%wd*X^%wd*Z^%wd",
                 (A->exps[i] >> (2*(FLINT_BITS/3))) & LOW_THIRD_MASK,
                 (A->exps[i] >> (1*(FLINT_BITS/3))) & LOW_THIRD_MASK,
                 (A->exps[i] >> (0*(FLINT_BITS/3))) & LOW_THIRD_MASK);
    }
}


void fmpz_mpolyu_sort_terms(
    fmpz_mpolyu_t A,
    const fmpz_mpoly_ctx_t ctx)
{
    slong i, j;
    for (i = 1; i < A->length; i++)
    {
        j = i;
        while (j > 0 && A->exps[j - 1] < A->exps[j])
        {
            ulong t = A->exps[j];
            A->exps[j] = A->exps[j - 1];
            A->exps[j - 1] = t;
            fmpz_mpoly_swap(A->coeffs + j - 1, A->coeffs + j, ctx);
            j--;
        }
    }
}


void nmod_mpoly_set_skel(
    nmod_mpolyc_t S,
    const nmod_mpoly_ctx_t ctx_sp,
    const fmpz_mpoly_t A,
    const mp_limb_t * alpha,
    const fmpz_mpoly_ctx_t ctx);


nmod_zip_find_coeffs_ret_t nmod_mpolyu_zip_find_coeffs_new(
    ulong Xdeg,
    fmpz_mpolyu_t H,
    const fmpz_mpolyu_t B,
    const fmpz_mpoly_ctx_t ctx,
    const mp_limb_t * alpha,
    nmod_zip_mpolyu_t Z,
    const nmod_mpoly_ctx_t nctx)
{
    slong i, j, k;
    nmod_poly_t master;
    nmod_mpolyc_t mon;
    slong N = mpoly_words_per_exp_sp(H->bits, ctx->minfo);
    slong lc_count;

printf("nmod_mpolyu_zip_find_coeffs_new called\n");

    FLINT_ASSERT(nctx->minfo->nvars == ctx->minfo->nvars);
    FLINT_ASSERT(H->bits == B->bits);

    nmod_mpolyc_init(mon);
    nmod_poly_init_mod(master, nctx->ffinfo->mod);

    fmpz_mpolyu_fit_length(H, Z->length, ctx);
    H->length = 0;

    lc_count = 0;

    for (i = 0; i < Z->length; i++)
    {
        fmpz_mpoly_struct * Hc = H->coeffs + H->length;
        nmod_zip_struct * Zc = Z->coeffs + i;
/*        ulong y = (Z->exps[i] >> (2*(FLINT_BITS/3))) & LOW_THIRD_MASK;*/
        ulong x = (Z->exps[i] >> (1*(FLINT_BITS/3))) & LOW_THIRD_MASK;
        ulong z = (Z->exps[i] >> (0*(FLINT_BITS/3))) & LOW_THIRD_MASK;

        /* deal with leading coeffs later */
        if (x >= Xdeg)
        {
            lc_count++;
            continue;
        }

        for (j = 0; j < B->length; j++)
        {
            fmpz_mpoly_struct * Bc = B->coeffs + j;
/*            ulong By = (B->exps[j] >> (2*(FLINT_BITS/3))) & LOW_THIRD_MASK;*/
            ulong Bx = (B->exps[j] >> (1*(FLINT_BITS/3))) & LOW_THIRD_MASK;
            ulong Bz = (B->exps[j] >> (0*(FLINT_BITS/3))) & LOW_THIRD_MASK;
            if (Bx == x && Bz == z)
            {
                H->exps[H->length] = Z->exps[i];

                nmod_mpoly_set_skel(mon, nctx, Bc, alpha, ctx);
                nmod_zip_set_mlength(Zc, Bc->length);
                for (k = 0; k < Bc->length; k++)
                    Zc->monomials[k] = mon->coeffs[k];

                switch (nmod_zip_find_coeffs_new(Zc->coeffs, Zc->monomials, Bc->length,
                                    Zc->evals, Z->pointcount, master, nctx->ffinfo))
                {
                    default:
                        FLINT_ASSERT(0);
                    case 0:
                        return 0;
                    case -1:
                        return -1;
                    case 1:
                        fmpz_mpoly_fit_length(Hc, Bc->length, ctx);
                        Hc->length = Bc->length;
                        for (k = 0; k < Bc->length; k++)
                        {
                            mp_limb_t V0 = Zc->coeffs[k];
                            fmpz_set_ui(Hc->coeffs + k, V0);
                            if (nctx->ffinfo->mod.n - V0 < V0)
                                fmpz_sub_ui(Hc->coeffs + k, Hc->coeffs + k, nctx->ffinfo->mod.n);
                            mpoly_monomial_set(Hc->exps + N*k, Bc->exps + N*k, N);
                        }
                        goto next;
                }
            }
        }
        return 0;
next:
        H->length++;
    }

    fmpz_mpolyu_fit_length(H, H->length + lc_count, ctx);

    for (j = 0; j < B->length; j++)
    {
        fmpz_mpoly_struct * Hc = H->coeffs + H->length;
        fmpz_mpoly_struct * Bc = B->coeffs + j;
/*        ulong By = (B->exps[j] >> (2*(FLINT_BITS/3))) & LOW_THIRD_MASK;*/
        ulong Bx = (B->exps[j] >> (1*(FLINT_BITS/3))) & LOW_THIRD_MASK;
/*        ulong Bz = (B->exps[j] >> (0*(FLINT_BITS/3))) & LOW_THIRD_MASK;*/

        if (Bx != Xdeg)
            continue;

        H->exps[H->length] = B->exps[j];

        fmpz_mpoly_fit_length(Hc, Bc->length, ctx);
        Hc->length = Bc->length;

        for (k = 0; k < Bc->length; k++)
        {
            mp_limb_t V0 = fmpz_fdiv_ui(Bc->coeffs + k, nctx->ffinfo->mod.n);
            fmpz_set_ui(Hc->coeffs + k, V0);
            if (nctx->ffinfo->mod.n - V0 < V0)
                fmpz_sub_ui(Hc->coeffs + k, Hc->coeffs + k, nctx->ffinfo->mod.n);
            mpoly_monomial_set(Hc->exps + N*k, Bc->exps + N*k, N);
        }

        H->length++;
    }

    fmpz_mpolyu_sort_terms(H, ctx);

    return nmod_zip_find_coeffs_good;
}



/*
    A and B are in ZZ[Y,X,Z](x1,...,xn) where n = ctx->minfo->nvars
    A = prod_{0<=i<r} B[i] mod (Y - beta)
    try to lift to a true factorization, return:
        -1: inconclusive, don't try again
        0: lift is impossible without changing the lc_X
        1: the B[i] were successfully lifted - yay!

    details:
        make assumed forms H[i] mod p0 for 0<=i<r
        Then use this assumed form to generate a bunch more lifts
        mod p1, p2, ... via zippel until p0*p1*... exceeds Mignotte.
*/
static int lift_uuu(
    slong r,                /* r = number of factors */
    fmpz_mpolyu_struct * B,
    const fmpz_mpolyu_t A,
    const fmpz_t beta,
    const slong * degs,     /* degs[i] = deg_xi(A) for 0 <= i < n */
    slong degree_inner,     /* degree_inner = deg_X(A) */
    const fmpz_mpoly_ctx_t ctx,
    fmpz_mpoly_factor_t fac_org,    /* in ZZ[X, Z, x1, ..., xn, Y] */
    const fmpz_mpoly_t A_org,       /* ditto */
    const fmpz_mpoly_ctx_t ctx_org) /* ctx->minfo->nvars = n + 3 */
{
    return -1;
#if 0
    int changed, success;
    flint_bitcnt_t bits = A->bits;
    mpoly_bma_interpolate_ctx_t Ictx;
    nmod_mpoly_ctx_t ctx_sp;
    fmpz_mod_mpoly_ctx_t ctx_mp;
    fmpz_mpolyu_struct * H;
    fmpz_mpoly_t T1, T2;
    /* multi precision workspace */
    fmpz_mod_bma_mpoly_struct * Lambda_mp;
    fmpz_mod_mpoly_t Aeval_mp;
    fmpz_mod_mpoly_struct * Beval_mp;
    fmpz_mod_mpolyn_struct * BBeval_mp;
    fmpz_mpolycu_t Ainc_mp, Acur_mp, Ared_mp;
    fmpz_mpolycu_struct * Binc_mp, * Bcur_mp, * Bred_mp;
    fmpz_t p, pm1, sshift_mp, image_count_mp, beta_mod_mp;
    fmpz * checkalpha_mp;
    /* single precision workspace */
    nmod_bma_mpoly_struct * Lambda_sp;
    nmod_mpoly_t Aeval_sp;
    nmod_mpoly_struct * Beval_sp;
    nmod_mpolyn_struct * BBeval_sp;
    nmod_mpolycu_t Ainc_sp, Acur_sp, Ared_sp;
    nmod_mpolycu_struct * Binc_sp, * Bcur_sp, * Bred_sp;
    mp_limb_t p_sp, sshift_sp, image_count_sp;
    mp_limb_t * checkalpha_sp;
    /* misc */
    slong i, j;
    flint_rand_t randstate;
    fmpz_t subprod;
    int unlucky_count;
    fmpz_t Hmodulus;
    nmod_zip_mpolyu_struct * Z;
    nmod_mpolyu_struct * M;
    slong zip_evals;
    ulong * Xdegs;
    nmod_mpoly_ctx_t nctx;
    slong v = ctx->minfo->nvars + 2;
    const char * vars[] = {"x", "y", "z", "w", "?", "?", "?"};
/*    const char * varsYXZ[] = {"Y", "X", "Z"};*/

    nmod_mpoly_ctx_init(nctx, ctx->minfo->nvars, ORD_LEX, 2);


flint_printf("lift_uuu called  bits: %wu\n", bits);
printf("A: "); fmpz_mpolyuuu_print_pretty(A, vars[v], vars[0], vars[1], vars + 2, ctx); printf("\n");
    for (i = 0; i < r; i++)
    {
        FLINT_ASSERT(B[i].bits == bits);
flint_printf("B[%wd]: ", i); fmpz_mpolyuuu_print_pretty(B + i, vars[v], vars[0], vars[1], vars + 2, ctx); printf("\n");
    }

    FLINT_ASSERT(bits == A->bits);
    for (i = 0; i < r; i++)
        FLINT_ASSERT(bits == B[i].bits);

    /* let's initialize everything at once to avoid complicated cleanup */

    flint_randinit(randstate);
    fmpz_init(p);
    fmpz_init(pm1); /* p - 1 */
    fmpz_init(image_count_mp);
    fmpz_init(subprod);
    fmpz_init(sshift_mp);
    fmpz_init(beta_mod_mp);
    fmpz_init(Hmodulus);

    mpoly_bma_interpolate_ctx_init(Ictx, ctx->minfo->nvars);

    /* multiprecision workspace "sp" */
    fmpz_set_ui(p, 2);    /* modulus no care */
    fmpz_mod_mpoly_ctx_init(ctx_mp, 3, ORD_LEX, p); /* modulus no care */
    Lambda_mp = (fmpz_mod_bma_mpoly_struct *) flint_malloc(r*sizeof(fmpz_mod_bma_mpoly_struct));
    for (i = 0; i < r; i++)
        fmpz_mod_bma_mpoly_init(Lambda_mp + i);

    fmpz_mod_mpoly_init3(Aeval_mp, 0, FLINT_BITS/3, ctx_mp);
    Beval_mp = (fmpz_mod_mpoly_struct *) flint_malloc(r*sizeof(fmpz_mod_mpoly_struct));
    BBeval_mp = (fmpz_mod_mpolyn_struct *) flint_malloc(r*sizeof(fmpz_mod_mpolyn_struct));
    for (i = 0; i < r; i++)
    {
        fmpz_mod_mpoly_init3(Beval_mp + i, 0, FLINT_BITS/3, ctx_mp);
        fmpz_mod_mpolyn_init(BBeval_mp + i, FLINT_BITS/3, ctx_mp);
    }

    fmpz_mpolycu_init(Ainc_mp);
    fmpz_mpolycu_init(Acur_mp);
    fmpz_mpolycu_init(Ared_mp);

    Binc_mp = (fmpz_mpolycu_struct *) flint_malloc(r*sizeof(fmpz_mpolycu_struct));
    Bcur_mp = (fmpz_mpolycu_struct *) flint_malloc(r*sizeof(fmpz_mpolycu_struct));
    Bred_mp = (fmpz_mpolycu_struct *) flint_malloc(r*sizeof(fmpz_mpolycu_struct));
    for (i = 0; i < r; i++)
    {
        fmpz_mpolycu_init(Binc_mp + i);
        fmpz_mpolycu_init(Bcur_mp + i);
        fmpz_mpolycu_init(Bred_mp + i);
    }

    checkalpha_mp = (fmpz *) flint_malloc(ctx->minfo->nvars*sizeof(fmpz));
    for (i = 0; i < ctx->minfo->nvars; i++)
        fmpz_init(checkalpha_mp + i);

    /* machine precision workspace "sp" */
    nmod_mpoly_ctx_init(ctx_sp, 3, ORD_LEX, 2); /* modulus no care */
    Lambda_sp = (nmod_bma_mpoly_struct *) flint_malloc(r*sizeof(nmod_bma_mpoly_struct));
    for (i = 0; i < r; i++)
        nmod_bma_mpoly_init(Lambda_sp + i);

    nmod_mpoly_init3(Aeval_sp, 0, FLINT_BITS/3, ctx_sp);
    Beval_sp = (nmod_mpoly_struct *) flint_malloc(r*sizeof(nmod_mpoly_struct));
    BBeval_sp = (nmod_mpolyn_struct *) flint_malloc(r*sizeof(nmod_mpolyn_struct));
    for (i = 0; i < r; i++)
    {
        nmod_mpoly_init3(Beval_sp + i, 0, FLINT_BITS/3, ctx_sp);
        nmod_mpolyn_init(BBeval_sp + i, FLINT_BITS/3, ctx_sp);
    }
    nmod_mpolycu_init(Ainc_sp);
    nmod_mpolycu_init(Acur_sp);
    nmod_mpolycu_init(Ared_sp);

    Binc_sp = (nmod_mpolycu_struct *) flint_malloc(r*sizeof(nmod_mpolycu_struct));
    Bcur_sp = (nmod_mpolycu_struct *) flint_malloc(r*sizeof(nmod_mpolycu_struct));
    Bred_sp = (nmod_mpolycu_struct *) flint_malloc(r*sizeof(nmod_mpolycu_struct));
    for (i = 0; i < r; i++)
    {
        nmod_mpolycu_init(Binc_sp + i);
        nmod_mpolycu_init(Bcur_sp + i);
        nmod_mpolycu_init(Bred_sp + i);
    }

    checkalpha_sp = (mp_limb_t *) flint_malloc(ctx->minfo->nvars*sizeof(mp_limb_t));

    /* the zippler */
    Z = (nmod_zip_mpolyu_struct *) flint_malloc(r*sizeof(nmod_zip_mpolyu_struct));
    for (i = 0; i < r; i++)
        nmod_zip_mpolyu_init(Z + i);

    H = (fmpz_mpolyu_struct *) flint_malloc(r*sizeof(fmpz_mpolyu_struct));
    for (i = 0; i < r; i++)
        fmpz_mpolyu_init(H + i, bits, ctx);

    M = (nmod_mpolyu_struct *) flint_malloc(r*sizeof(nmod_mpolyu_struct));
    for (i = 0; i < r; i++)
        nmod_mpolyu_init(M + i, bits, nctx);

    Xdegs = (ulong *) flint_malloc(r*sizeof(ulong));

    /*
        initial choices for the ksub degrees are the strict degree bounds on A
        (could be lower if we tried to get degree bounds on the factors ?)
    */
    for (i = 0; i < ctx->minfo->nvars; i++)
    {
        Ictx->degbounds[i] = degs[i] + 1;
        Ictx->subdegs[i] = Ictx->degbounds[i];
    }

    fmpz_mpoly_init3(T1, 0, A->bits, ctx_org);
    fmpz_mpoly_init3(T2, 0, A->bits, ctx_org);

    fmpz_one(subprod);
    for (i = 0; i < ctx->minfo->nvars; i++)
    {
        if ((slong)(Ictx->subdegs[i]) < 0)
        {
            /* ksub has overflown - absolute falure */
            success = 0;
            goto cleanup;
        }
        fmpz_mul_ui(subprod, subprod, Ictx->subdegs[i]);
    }
/*
    if (fmpz_is_zero(beta))
        goto pick_bma_prime;
*/
    /*
        Monagan's optional optimization:
        Use the exponents from the previous factorization instead of finding
        them from scratch. This is a bit tricky in the non-monic case where
        all lc_X have been predetermined. Recall we start with
            A = prod_i B[i] mod (Y - beta).
        where each B[i] has the form
            sum_j c_j(X,Z,x1,....,xn) Y^j
        We make the assumption that the monomials (sans Y) in the true factors
        are union_j monomials in the c_j.
        If this is not true, zip will catch this, and then we do indeed
        find them from scratch (or punt to another lifter)

        The following is done for each factor B.

        M holds terms of the form
            sum_{j,k} Y^0*X^j*Z^k*poly_{j,k}(x1,...,xn)
        for j < deg_X(B). Specifically, M holds
            1. the monomials of poly_{j,k}
            2. the evaluations of these monomials
            3. space for the poly with roots at these evaluations
        M is determined completely by B before any images are determined.

        We get a bunch of images of the form sum_{i,j,k} c_{i,j,k} Y^i*X^j*Z^k

        For c_{i,j,k} != 0 and j < deg_X(B):
            If X^j*Z^k does not appear in M, then punt to bma, else
            use the assumed form from M to make the corresponding entry in H

        For j = deg_X(B) and each i, k:
            copy the coefficient of Y^i*X^j*Z^k of B mods p into H
    */

    zip_evals = 0;
    for (i = 0; i < r; i++)
    {
        slong this_zip_evals;
        this_zip_evals = zip_form_new(Xdegs + i, M + i, nctx, B + i, ctx);
        zip_evals = FLINT_MAX(zip_evals, this_zip_evals);
    }
    zip_evals += 1;

    p_sp = UWORD(1) << (FLINT_BITS - 2);

pick_monagan_prime:

    if (p_sp >= UWORD_MAX_PRIME)
        goto pick_bma_prime;

    p_sp = n_nextprime(p_sp, 1);

flint_printf("monagan prime: %wu\n", p_sp); fflush(stdout);
usleep(1000000);

    nmod_mpoly_ctx_set_modulus(nctx, p_sp);
    nmod_mpoly_ctx_set_modulus(ctx_sp, p_sp);
    /* unfortunate nmod_poly's need mod set */
    for (i = 0; i < r; i++)
        nmod_mpolyn_set_mod(BBeval_sp + i, ctx_sp->ffinfo->mod);

    FLINT_ASSERT(p_sp > 3);
    for (i = 0; i < ctx->minfo->nvars; i++)
        checkalpha_sp[i] = n_urandint(randstate, p_sp - 3) + 2;

    /* set evaluation of monomials */
    nmod_mpolyu_set_skel(Ainc_sp, ctx_sp, A, checkalpha_sp, ctx);
    for (i = 0; i < r; i++)
        nmod_mpolyu_set_skel(Binc_sp + i, ctx_sp, B + i, checkalpha_sp, ctx);

    /* set reduction of coeffs */
    nmod_mpolyu_red_skel(Ared_sp, A, ctx_sp->ffinfo);
    for (i = 0; i < r; i++)
        nmod_mpolyu_red_skel(Bred_sp + i, B + i, ctx_sp->ffinfo);

    /* copy evaluation of monomials */
    nmod_mpolyu_copy_skel(Acur_sp, Ainc_sp);
    for (i = 0; i < r; i++)
        nmod_mpolyu_copy_skel(Bcur_sp + i, Binc_sp + i);

next_monagan_image:

printf("next_monagan_image\n"); fflush(stdout);

    nmod_tpoly_use_skel_mul(Aeval_sp, A, Ared_sp, Acur_sp, Ainc_sp, ctx_sp);
    for (i = 0; i < r; i++)
        nmod_tpoly_use_skel_mul(Beval_sp + i, B + i, Bred_sp + i, Bcur_sp + i, Binc_sp + i, ctx_sp);

    success = n_polyu3_hensel_lift(r, BBeval_sp, Aeval_sp, Beval_sp, fmpz_fdiv_ui(beta, p_sp), degree_inner, ctx_sp);
    if (success < 0)
    {
        goto pick_zip_prime;
    }
    else if (success == 0)
    {
        /* lifting is impossible - rare */
        goto cleanup;
    }

    /* update the zippler */
    for (i = 0; i < r; i++)
    {
/*
flint_printf("updating zippler[%wd] with: "); nmod_mpolyn_print_pretty(BBeval_sp, varsYXZ, ctx_sp); printf("\n");
*/
        nmod_zip_mpolyuuu_add_point_allow_new(Z + i, BBeval_sp + i);
    }
/*
flint_printf("Z[0].pointcount: %wd  ", Z[0].pointcount);
flint_printf("zip_evals: %wd\n", zip_evals);
*/
    if (Z[0].pointcount < zip_evals)
        goto next_monagan_image;
/*
printf("finding coeffs\n");
*/
    for (i = 0; i < r; i++)
    {
/*
flint_printf("Z[%wd]: ", i); nmod_zip_mpolyuuu_print(Z + i); printf("\n");
*/
        /* figure out the form of H + i from Z + i */
        switch (nmod_mpolyu_zip_find_coeffs_new(Xdegs[i], H + i, B + i, ctx,
                                                 checkalpha_sp, Z + i, nctx))
        {
            default:
                FLINT_ASSERT(0);
            case nmod_zip_find_coeffs_no_match:
                /*  The collection of images in Fp'[Y,X,Z] could not be coerced
                    into the assumed form in [Y,X,Z][x1, ..., xn]. */
                goto pick_bma_prime;
            case nmod_zip_find_coeffs_non_invertible:
                /* The unlikely case where the evaluation points alpha produced
                   a singular Vandermonde matrix. Assumed form is not nec wrong. */
                goto pick_monagan_prime;
            case nmod_zip_find_coeffs_good:
                (void)NULL;

flint_printf("H[%wd]: ", i);
fmpz_mpolyuuu_print_pretty(H + i, vars[v], vars[0], vars[1], vars + 2, ctx);
printf("\n");
fflush(stdout);

        }
    }

printf("monagan good\n");

    fmpz_set_ui(Hmodulus, p_sp);

    goto start_zip;

pick_bma_prime:

    if (fmpz_cmp(p, subprod) < 0)
        fmpz_set(p, subprod);

    success = fmpz_next_smooth_prime(p, p);
    fmpz_sub_ui(pm1, p, 1);
    if (!success)
    {
        /* ran out of smooth primes - absolute failure */
        success = 0;
        goto cleanup;
    }

printf("------- bma p: "); fmpz_print(p); printf(" --------\n"); fflush(stdout);
usleep(1000000);

    unlucky_count = 0;

    if (fmpz_abs_fits_ui(p))
    {
        p_sp = fmpz_get_ui(p);
        sshift_sp = 1;

        nmod_mpoly_ctx_set_modulus(ctx_sp, p_sp);
        nmod_discrete_log_pohlig_hellman_precompute_prime(Ictx->dlogenv_sp, p_sp);

        for (i = 0; i < r; i++)
        {
            nmod_bma_mpoly_reset_prime(Lambda_sp + i, ctx_sp->ffinfo);
            nmod_bma_mpoly_zero(Lambda_sp + i);
        }

        /* unfortunate nmod_poly's store their own ctx :( */
        for (i = 0; i < r; i++)
        {
            nmod_mpolyn_set_mod(BBeval_sp + i, ctx_sp->ffinfo->mod);
        }

        FLINT_ASSERT(sshift_sp == 1);
        nmod_mpoly_bma_interpolate_alpha_powers(checkalpha_sp, sshift_sp,
                                                    Ictx, ctx, ctx_sp->ffinfo);

        /* set evaluation of monomials */
        nmod_mpolyu_set_skel(Ainc_sp, ctx_sp, A, checkalpha_sp, ctx);
        for (i = 0; i < r; i++)
            nmod_mpolyu_set_skel(Binc_sp + i, ctx_sp, B + i, checkalpha_sp, ctx);

        /* set reduction of coeffs */
        nmod_mpolyu_red_skel(Ared_sp, A, ctx_sp->ffinfo);
        for (i = 0; i < r; i++)
            nmod_mpolyu_red_skel(Bred_sp + i, B + i, ctx_sp->ffinfo);

        /* copy evaluation of monomials */
        nmod_mpolyu_copy_skel(Acur_sp, Ainc_sp);
        for (i = 0; i < r; i++)
            nmod_mpolyu_copy_skel(Bcur_sp + i, Binc_sp + i);

        image_count_sp = 0;

    next_bma_image_sp:

        /* image count is also the current power of alpha we are evaluating */
        image_count_sp++;

        FLINT_ASSERT(sshift_sp + Lambda_sp->pointcount == image_count_sp);

        if (image_count_sp >= p_sp - 1)
        {
            /* out of evaluation points alpha^image_count in Fp* */
            goto pick_bma_prime;
        }

        nmod_tpoly_use_skel_mul(Aeval_sp, A, Ared_sp, Acur_sp, Ainc_sp, ctx_sp);
        for (i = 0; i < r; i++)
            nmod_tpoly_use_skel_mul(Beval_sp + i, B + i, Bred_sp + i, Bcur_sp + i, Binc_sp + i, ctx_sp);

        success = nmod_tpoly_hensel_lift(r, BBeval_sp, Aeval_sp, Beval_sp, fmpz_fdiv_ui(beta, p_sp), degree_inner, ctx_sp);
        if (success < 0)
        {
            if (++unlucky_count > 2)
                goto pick_bma_prime;
            sshift_sp += Lambda_sp[0].pointcount + 1;
            for (i = 0; i < r; i++)
                nmod_bma_mpoly_zero(Lambda_sp + i);
            goto next_bma_image_sp;
        }
        else if (success == 0)
        {
            /* lifting is impossible - rare */
            goto cleanup;
        }

        for (i = 0; i < r; i++)
            nmod_bma_mpoly_add_point_tpoly(Lambda_sp + i, BBeval_sp + i, ctx_sp);

        /* all of the lambda's should have the same point count */
        if ((Lambda_sp[0].pointcount & 1) != 0)
            goto next_bma_image_sp;

        for (i = 0; i < r; i++)
        {
            changed = nmod_bma_mpoly_reduce(Lambda_sp + i);
            if (changed)
                goto next_bma_image_sp;
        }

        for (i = 0; i < r; i++)
        {
            success = nmod_bma_mpoly_get_fmpz_mpolyu(H + i, ctx, sshift_sp,
                                              Lambda_sp + i, Ictx, ctx_sp->ffinfo);
        }
    }
    else
    {
        fmpz_one(sshift_mp);
        fmpz_mod_ctx_set_modulus(ctx_mp->ffinfo, p);
        fmpz_mod_discrete_log_pohlig_hellman_precompute_prime(Ictx->dlogenv, p);

        for (i = 0; i < r; i++)
        {
            fmpz_mod_bma_mpoly_reset_prime(Lambda_mp + i, ctx_mp->ffinfo);
            fmpz_mod_bma_mpoly_zero(Lambda_mp + i);
        }

        /* unfortunate fmpz_mod_poly's store their own ctx :( */
        for (i = 0; i < r; i++)
        {
            fmpz_mod_mpolyn_set_modulus(BBeval_mp + i, ctx_mp->ffinfo);
        }

        FLINT_ASSERT(fmpz_is_one(sshift_mp));
        fmpz_mod_mpoly_bma_interpolate_alpha_powers(checkalpha_mp, sshift_mp,
                                                    Ictx, ctx, ctx_mp->ffinfo);

        /* set evaluation of monomials */
        fmpz_mod_mpolyu_set_skel(Ainc_mp, ctx_mp, A, checkalpha_mp, ctx);
        for (i = 0; i < r; i++)
            fmpz_mod_mpolyu_set_skel(Binc_mp + i, ctx_mp, B + i, checkalpha_mp, ctx);

        /* set reduction of coeffs */
        fmpz_mod_mpolyu_red_skel(Ared_mp, A, ctx_mp->ffinfo);
        for (i = 0; i < r; i++)
            fmpz_mod_mpolyu_red_skel(Bred_mp + i, B + i, ctx_mp->ffinfo);

        /* copy evaluation of monomials */
        fmpz_mod_mpolyu_copy_skel(Acur_mp, Ainc_mp);
        for (i = 0; i < r; i++)
            fmpz_mod_mpolyu_copy_skel(Bcur_mp + i, Binc_mp + i);

        fmpz_zero(image_count_mp);

    next_bma_image_mp:

        /* image count is also the current power of alpha we are evaluating */
        fmpz_add_ui(image_count_mp, image_count_mp, 1);

    #if WANT_ASSERT
        /* image_count == sshift + Lambda->pointcount */
        {
            fmpz_t t;
            fmpz_init(t);
            fmpz_add_ui(t, sshift_mp, Lambda_mp[0].pointcount);
            FLINT_ASSERT(fmpz_equal(t, image_count_mp));
            fmpz_clear(t);
        }
    #endif

        if (fmpz_cmp(image_count_mp, pm1) >= 0)
        {
            /* out of evaluation points alpha^image_count in Fp* */
            goto pick_bma_prime;
        }

        fmpz_mod_tpoly_use_skel_mul(Aeval_mp, A, Ared_mp, Acur_mp, Ainc_mp, ctx_mp);
        for (i = 0; i < r; i++)
            fmpz_mod_tpoly_use_skel_mul(Beval_mp + i, B + i, Bred_mp + i, Bcur_mp + i, Binc_mp + i, ctx_mp);

        fmpz_mod(beta_mod_mp, beta, p);
        success = fmpz_mod_tpoly_hensel_lift(r, BBeval_mp, Aeval_mp, Beval_mp, beta_mod_mp, degree_inner, ctx_mp);
        if (success < 0)
        {
            if (++unlucky_count > 2)
                goto pick_bma_prime;
            fmpz_add_ui(sshift_mp, sshift_mp, Lambda_mp[0].pointcount + 1);
            for (i = 0; i < r; i++)
                fmpz_mod_bma_mpoly_zero(Lambda_mp + i);
            goto next_bma_image_mp;
        }
        else if (success == 0)
        {
            /* lifting is impossible - rare */
            goto cleanup;
        }

        for (i = 0; i < r; i++)
            fmpz_mod_bma_mpoly_add_point_tpoly(Lambda_mp + i, BBeval_mp + i, ctx_mp);

        if ((Lambda_mp[0].pointcount & 1) != 0)
            goto next_bma_image_mp;

        for (i = 0; i < r; i++)
        {
            changed = fmpz_mod_bma_mpoly_reduce(Lambda_mp + i);
            if (changed)
                goto next_bma_image_mp;
        }

        for (i = 0; i < r; i++)
        {
            success = fmpz_mod_bma_mpoly_get_fmpz_mpolyu(H + i, ctx, sshift_mp,
                                              Lambda_mp + i, Ictx, ctx_mp->ffinfo);
        }
    }

    /* assume that the H[i] are correct modulo Hmodulus = p */
    fmpz_set(Hmodulus, p);

start_zip:

    FLINT_ASSERT(fmpz_cmp_ui(Hmodulus, 3) > 0);

    /* find number of evals for zip interp */
    FLINT_ASSERT(H->length > 0);
    zip_evals = 2;
    for (i = 0; i < r; i++)
    {
        for (j = 0; j < H[i].length; j++)
            zip_evals = FLINT_MAX(zip_evals, H[i].coeffs[j].length);
    }

    for (i = 0; i < r; i++)
        nmod_zip_mpolyu_fit_poly(Z + i, H + i, zip_evals);

    p_sp = UWORD(1) << (FLINT_BITS - 2);

pick_zip_prime:
    /*
        Get a new machine prime for zippel interpolation.
        The H[i] are currently interpolated modulo Hmodulus.
    */
    if (p_sp >= UWORD_MAX_PRIME)
    {
        /* ran out of machine primes - absolute failure */
        success = -1;
        goto cleanup;
    }
    p_sp = n_nextprime(p_sp, 1);

    if (0 == fmpz_fdiv_ui(Hmodulus, p_sp))
    {
        goto pick_zip_prime;
    }

flint_printf("zip prime: %wu\n", p_sp); fflush(stdout);
usleep(1000000);

    nmod_mpoly_ctx_set_modulus(ctx_sp, p_sp);
    /* unfortunate nmod_poly's need mod set */
    for (i = 0; i < r; i++)
        nmod_mpolyn_set_mod(BBeval_sp + i, ctx_sp->ffinfo->mod);

    FLINT_ASSERT(p_sp > 3);
    for (i = 0; i < ctx->minfo->nvars; i++)
        checkalpha_sp[i] = n_urandint(randstate, p_sp - 3) + 2;

    /* set up the zippler */
    for (i = 0; i < r; i++)
        nmod_zip_mpolyu_set_skel(Z + i, ctx_sp, H + i, checkalpha_sp, ctx);

    /* set evaluation of monomials */
    nmod_mpolyu_set_skel(Ainc_sp, ctx_sp, A, checkalpha_sp, ctx);
    for (i = 0; i < r; i++)
        nmod_mpolyu_set_skel(Binc_sp + i, ctx_sp, B + i, checkalpha_sp, ctx);

    /* set reduction of coeffs */
    nmod_mpolyu_red_skel(Ared_sp, A, ctx_sp->ffinfo);
    for (i = 0; i < r; i++)
        nmod_mpolyu_red_skel(Bred_sp + i, B + i, ctx_sp->ffinfo);

    /* copy evaluation of monomials */
    nmod_mpolyu_copy_skel(Acur_sp, Ainc_sp);
    for (i = 0; i < r; i++)
        nmod_mpolyu_copy_skel(Bcur_sp + i, Binc_sp + i);

next_zip_image:

printf("next_zip_image\n"); fflush(stdout);

    nmod_tpoly_use_skel_mul(Aeval_sp, A, Ared_sp, Acur_sp, Ainc_sp, ctx_sp);
    for (i = 0; i < r; i++)
        nmod_tpoly_use_skel_mul(Beval_sp + i, B + i, Bred_sp + i, Bcur_sp + i, Binc_sp + i, ctx_sp);

    success = nmod_tpoly_hensel_lift(r, BBeval_sp, Aeval_sp, Beval_sp, fmpz_fdiv_ui(beta, p_sp), degree_inner, ctx_sp);
    if (success < 0)
    {
        goto pick_zip_prime;
    }
    else if (success == 0)
    {
        /* lifting is impossible - rare */
        goto cleanup;
    }

    /* update the zippler */
    for (i = 0; i < r; i++)
    {
/*
flint_printf("adding "); nmod_mpolyn_print_pretty(BBeval_sp + i, varsYXZ, ctx_sp); printf("\n");
flint_printf("to Z[%wd]: ", i); nmod_zip_mpolyuuu_print(Z + i); printf("\n");
*/
        success = nmod_zip_mpolyuuu_add_point(Z + i, BBeval_sp + i);
        if (!success)
        {
            /* An image in Fp'[Y,X,Z] did not match the assumed formed in [Y,X,Z] */
            goto pick_bma_prime;
        }
    }

#if WANT_ASSERT
    for (i = 1; i < r; i++)
        FLINT_ASSERT(Z[0].pointcount == Z[i].pointcount);
#endif

    /* all of the Z[i] have the same point count */
    if (Z[0].pointcount < zip_evals)
        goto next_zip_image;

    for (i = 0; i < r; i++)
    {
        switch (nmod_mpolyu_zip_find_coeffs(Z + i, ctx_sp))
        {
            default:
                FLINT_ASSERT(0);
            case nmod_zip_find_coeffs_no_match:
                /*  The collection of images in Fp'[Y,X,Z] could not be coerced
                    into the assumed form in [Y,X,Z][x1, ..., xn]. */
                goto pick_bma_prime;
            case nmod_zip_find_coeffs_non_invertible:
                /* The unlikely case where the evaluation points alpha produced
                   a singular Vandermonde matrix. Assumed form is not nec wrong. */
                goto pick_zip_prime;
            case nmod_zip_find_coeffs_good:
                NULL;
        }
    }

    changed = 0;
    for (i = 0; i < r; i++)
    {
        int this_changed;
        FLINT_ASSERT(bits == H[i].bits);
        this_changed = fmpz_mpolyu_addinterp_zip(H + i, Hmodulus, Z + i, ctx_sp->ffinfo);
        changed |= this_changed;

flint_printf("H[%wd](%d): ", i, success);
fmpz_mpolyuuu_print_pretty(H + i, vars[v], vars[0], vars[1], vars + 2, ctx);
printf("\n");
fflush(stdout);

    }

    fmpz_mul_ui(Hmodulus, Hmodulus, ctx_sp->ffinfo->mod.n);

    if (changed)
    {
        /* TODO: if the coefficients of the H[i] are getting too large? */
        goto pick_zip_prime;
    }

    /* write down the factorization */
    fmpz_mpoly_factor_fit_length(fac_org, r, ctx_org);
    fmpz_one(fac_org->constant);
    for (i = 0; i < r; i++)
    {
        fmpz_one(fac_org->exp + i);
        fmpz_mpoly_from_fmpz_mpolyuuu(ctx->minfo->nvars + 2, fac_org->poly + i, A->bits, ctx_org, H + i, ctx);
    }

    /* check the factorization */
    fmpz_mpoly_mul(T1, fac_org->poly + 0, fac_org->poly + 1, ctx_org);
    for (i = 2; i < r; i++)
    {
        fmpz_mpoly_mul(T2, T1, fac_org->poly + i, ctx_org);
        fmpz_mpoly_swap(T1, T2, ctx_org);
    }
    if (!fmpz_mpoly_equal(T1, A_org, ctx_org))
        goto pick_zip_prime;

    success = 1;

cleanup:

    return success;
#endif
}

static int play_lift(
    slong v,
    fmpz_mpoly_factor_t lfac,
    const fmpz * alpha,
    const fmpz_mpoly_t A,
    const slong * degs,
    const fmpz_mpoly_ctx_t ctx)
{
    slong i;
    slong r = lfac->num;
    const char * vars[] = {"x", "y", "z", "w", "?", "?", "?"};
    fmpz_mpolyu_t Au;
    fmpz_mpolyu_struct * Bu;
    flint_bitcnt_t bits;
    fmpz_mpoly_ctx_t ctx3, uctx;
    fmpz_mpoly_factor_t newlfac;

flint_printf(" play lift called v = %wd\n", v);
printf("lfac: "); fmpz_mpoly_factor_print_pretty(lfac, vars, ctx); printf("\n");
printf("A: "); fmpz_mpoly_print_pretty(A, vars, ctx); printf("\n");

    FLINT_ASSERT(2 < v);
    FLINT_ASSERT(v < ctx->minfo->nvars);

    fmpz_mpoly_ctx_init(ctx3, 3, ORD_LEX);
    fmpz_mpoly_ctx_init(uctx, v - 2, ORD_LEX);

    bits = 0;
    for (i = 0; i <= v; i++)
        bits = FLINT_MAX(bits, FLINT_BIT_COUNT(degs[i]));
    bits = mpoly_fix_bits(bits + 1, uctx->minfo);

    fmpz_mpolyu_init(Au, bits, uctx);
    fmpz_mpoly_to_fmpz_mpolyuuu(v, Au, uctx, A, ctx);
    Bu = (fmpz_mpolyu_struct *) flint_malloc(r*sizeof(fmpz_mpolyu_struct));
    for (i = 0; i < r; i++)
    {
        fmpz_mpolyu_init(Bu + i, bits, uctx);
        fmpz_mpoly_to_fmpz_mpolyuuu(v, Bu + i, uctx, lfac->poly + i, ctx);
    }

    fmpz_mpoly_factor_init(newlfac, ctx);
    lift_uuu(r, Bu, Au, alpha + r, degs + 2, degs[0], uctx,
             newlfac, A, ctx);

    return 0;
}


int fmpz_mpoly_factor_matches(
    const fmpz_mpoly_t A,
    const fmpz_mpoly_factor_t f,
    const fmpz_mpoly_ctx_t ctx)
{
    int matches;
    fmpz_mpoly_t T;
    fmpz_mpoly_init(T, ctx);
    fmpz_mpoly_factor_expand(T, f, ctx);
    matches = fmpz_mpoly_equal(T, A, ctx);
    fmpz_mpoly_clear(T, ctx);
    return matches;
}


/*********** mpoly_vec ************/

typedef struct
{
    fmpz_mpoly_struct * coeffs;
    slong alloc;
    slong length;
} fmpz_mpolyv_struct;

typedef fmpz_mpolyv_struct fmpz_mpolyv_t[1];

typedef struct
{
    nmod_mpoly_struct * coeffs;
    slong alloc;
    slong length;
} nmod_mpolyv_struct;

typedef nmod_mpolyv_struct nmod_mpolyv_t[1];


void fmpz_mpolyv_init(fmpz_mpolyv_t A, const fmpz_mpoly_ctx_t ctx)
{
    A->coeffs = NULL;
    A->alloc = 0;
    A->length = 0;
}

void nmod_mpolyv_init(nmod_mpolyv_t A, const nmod_mpoly_ctx_t ctx)
{
    A->coeffs = NULL;
    A->alloc = 0;
    A->length = 0;
}

void fmpz_mpolyv_clear(fmpz_mpolyv_t A, const fmpz_mpoly_ctx_t ctx)
{
    slong i;
    for (i = 0; i < A->alloc; i++)
        fmpz_mpoly_clear(A->coeffs + i, ctx);
    flint_free(A->coeffs);
}

void nmod_mpolyv_clear(nmod_mpolyv_t A, const nmod_mpoly_ctx_t ctx)
{
    slong i;
    for (i = 0; i < A->alloc; i++)
        nmod_mpoly_clear(A->coeffs + i, ctx);
    flint_free(A->coeffs);
}

void fmpz_mpolyv_swap(
    fmpz_mpolyv_t A,
    fmpz_mpolyv_t B,
    const fmpz_mpoly_ctx_t ctx)
{
   fmpz_mpolyv_struct t = *A;
   *A = *B;
   *B = t;
}

void nmod_mpolyv_swap(
    nmod_mpolyv_t A,
    nmod_mpolyv_t B,
    const nmod_mpoly_ctx_t ctx)
{
   nmod_mpolyv_struct t = *A;
   *A = *B;
   *B = t;
}

void fmpz_mpolyv_print_pretty(
    const fmpz_mpolyv_t poly,
    const char ** x,
    const fmpz_mpoly_ctx_t ctx)
{
    slong i;
    for (i = 0; i < poly->length; i++)
    {
        flint_printf("coeff[%wd]: ", i);
        fmpz_mpoly_print_pretty(poly->coeffs + i, x, ctx);
        flint_printf("\n");
    }
}

void nmod_mpolyv_print_pretty(
    const nmod_mpolyv_t poly,
    const char ** x,
    const nmod_mpoly_ctx_t ctx)
{
    slong i;
    for (i = 0; i < poly->length; i++)
    {
        flint_printf("coeff[%wd]: ", i);
        nmod_mpoly_print_pretty(poly->coeffs + i, x, ctx);
        flint_printf("\n");
    }
}

void fmpz_mpolyv_fit_length(
    fmpz_mpolyv_t A,
    slong length,
    const fmpz_mpoly_ctx_t ctx)
{
    slong i;
    slong old_alloc = A->alloc;
    slong new_alloc = FLINT_MAX(length, 2*A->alloc);

    if (length > old_alloc)
    {
        if (old_alloc == 0)
        {
            A->coeffs = (fmpz_mpoly_struct *) flint_malloc(
                                          new_alloc*sizeof(fmpz_mpoly_struct));
        }
        else
        {
            A->coeffs = (fmpz_mpoly_struct *) flint_realloc(A->coeffs,
                                          new_alloc*sizeof(fmpz_mpoly_struct));
        }

        for (i = old_alloc; i < new_alloc; i++)
        {
            fmpz_mpoly_init(A->coeffs + i, ctx);
        }
        A->alloc = new_alloc;
    }
}

void nmod_mpolyv_fit_length(
    nmod_mpolyv_t A,
    slong length,
    const nmod_mpoly_ctx_t ctx)
{
    slong i;
    slong old_alloc = A->alloc;
    slong new_alloc = FLINT_MAX(length, 2*A->alloc);

    if (length > old_alloc)
    {
        if (old_alloc == 0)
        {
            A->coeffs = (nmod_mpoly_struct *) flint_malloc(
                                          new_alloc*sizeof(nmod_mpoly_struct));
        }
        else
        {
            A->coeffs = (nmod_mpoly_struct *) flint_realloc(A->coeffs,
                                          new_alloc*sizeof(nmod_mpoly_struct));
        }

        for (i = old_alloc; i < new_alloc; i++)
        {
            nmod_mpoly_init(A->coeffs + i, ctx);
        }
        A->alloc = new_alloc;
    }
}

void fmpz_mpolyv_set_coeff(
    fmpz_mpolyv_t A,
    slong i,
    fmpz_mpoly_t c, /* clobbered */
    const fmpz_mpoly_ctx_t ctx)
{
    slong j;
    FLINT_ASSERT(!fmpz_mpoly_is_zero(c, ctx));
    fmpz_mpolyv_fit_length(A, i + 1, ctx);
    for (j = A->length; j < i; j++)
        fmpz_mpoly_zero(A->coeffs + j, ctx);
    fmpz_mpoly_swap(A->coeffs + i, c, ctx);
    A->length = FLINT_MAX(A->length, i + 1);
}


/************** subsets and tuples ****************************************/

void subset_first(fmpz_t a, slong n, slong r)
{
    FLINT_ASSERT(r >= 0);
    FLINT_ASSERT(r <= n);
    fmpz_one(a);
    fmpz_mul_2exp(a, a, r); 
    fmpz_sub_ui(a, a, 1);
}

int subset_next(fmpz_t a, const fmpz_t b, slong n)
{
    slong t1, t2, i;
    int r;
    if (a == b)
    {
        fmpz_t t;
        fmpz_init(t);
        r = subset_next(t, b, n);
        fmpz_swap(a, t);
        fmpz_clear(t);
        return r;
    }

    i = 0;
    while (i < n && fmpz_tstbit(b, i) == 0)
        i++;
    t1 = i;
/*flint_printf("t1: %wd\n", t1);*/
    while (i<n && fmpz_tstbit(b,i) == 1)
        i++;
    t2 = i;
/*flint_printf("t2: %wd\n", t2);*/
    if (t2 < n)
    {
        fmpz_t t;
        fmpz_init_set_ui(t, 1);
        fmpz_one(a);
        fmpz_mul_2exp(a, a, n - t2);
        fmpz_sub_ui(a, a, 1);
        fmpz_mul_2exp(a, a, t2);
        fmpz_and(a, b, a);
        fmpz_setbit(a, t2);
        if (t2 > t1)
            fmpz_mul_2exp(t, t, t2 - t1 - 1);
        fmpz_sub_ui(t, t, 1);
        fmpz_add(a, a, t);
        fmpz_clear(t);
        return 1;
    }
    else
    {
        return 0;
    }
}

void subset_print(const fmpz_t a, slong n)
{
    slong i;
    for (i = n - 1; i >= 0; i --)
    {
        flint_printf("%d",fmpz_tstbit(a, i));
    }
}

void subset_map_down(fmpz_t a, const fmpz_t b, const fmpz_t m)
{
    ulong i, j, bbits = fmpz_bits(b);

    FLINT_ASSERT(a != b);
    FLINT_ASSERT(a != m);

    j = 0;
    fmpz_zero(a);
    for (i = 0; i < bbits; i++)
    {
        if (fmpz_tstbit(b, i))
        {
            FLINT_ASSERT(!fmpz_tstbit(m, i));
            fmpz_setbit(a, j++);
        }
        else
        {
            j += !fmpz_tstbit(m, i);
        }
    }
}


void tuple_print(fmpz * alpha, slong n)
{
    slong j;
    for (j = 0; j < n; j++)
    {
        fmpz_print(alpha + j);
        flint_printf(j + 1 < n ? ", " : "\n");
    }
}


/* ensure that the first m values change upon the next call to tuple_next*/
void tuple_saturate(fmpz * alpha, slong n, slong m)
{
    slong i;

    for (i = m + 1; i < n; i++)
    {
        fmpz_add(alpha + m, alpha + m, alpha + i);
        fmpz_zero(alpha + i);
    }

    if (m < n && fmpz_is_zero(alpha + m))
    {
        for (i = 0; i < m; i++)
            if (!fmpz_is_zero(alpha + i))
                return;
        fmpz_one(alpha + m);
    }
}



void tuple_next(fmpz * alpha, slong n)
{
    slong i, t1, t2, t3;
    fmpz_t sum;

    fmpz_init(sum);
    for (i = 0; i < n; i++)
        fmpz_add(sum, sum, alpha + i);

    i = n - 1;
    while(i >= 0 && fmpz_is_zero(alpha + i))
        i--;
    t1 = i;
    while(i >= 0 && fmpz_cmp(alpha + i, sum) != 0)
        i--;
    t2 = i;
    while(i >= 0 && fmpz_cmp(alpha + i, sum) == 0)
        i--;
    t3 = i;

    if (t1 > 0 && t1 != t2)
    {
        fmpz_swap(alpha + t1, alpha + n - 1);
        fmpz_sub_ui(alpha + n - 1, alpha + n - 1, 1);
        fmpz_add_ui(alpha + t1 - 1, alpha + t1 - 1, 1);
    }
    else if (t1 > 0 && t1 == t2 && t3 >= 0)
    {
        fmpz_add_ui(alpha + t3, alpha + t3, 1);
        fmpz_zero(alpha + t3 + 1);
        fmpz_sub_ui(alpha + n - 1, sum, 1);
    }
    else
    {
        fmpz_add_ui(alpha + n - 1, alpha + 0, 1);
        if (n > 1)
            fmpz_zero(alpha + 0);
    }

    fmpz_clear(sum);
}


/********* univar ************************************************************/

int fmpz_mpoly_univar_content_mpoly(
    fmpz_mpoly_t g,
    const fmpz_mpoly_univar_t A,
    const fmpz_mpoly_ctx_t ctx)
{
    slong i;
    int success;

    fmpz_mpoly_zero(g, ctx);
    for (i = 0; i < A->length; i++)
    {
        success = fmpz_mpoly_gcd(g, g, A->coeffs + i, ctx);
        if (!success)
            return 0;
    }

    return 1;
}

void fmpz_mpoly_univar_divexact_mpoly(
    fmpz_mpoly_univar_t A,
    const fmpz_mpoly_t b,
    const fmpz_mpoly_ctx_t ctx)
{
    int success;
    slong i;

    for (i = 0; i < A->length; i++)
    {
        success = fmpz_mpoly_divides(A->coeffs + i, A->coeffs + i, b, ctx);
        FLINT_ASSERT(success);
    }
}

void _fmpz_mpoly_univar_shift_right(
    fmpz_mpoly_univar_t A,
    const fmpz_t shift,
    const fmpz_mpoly_ctx_t ctx)
{
    slong i;

    if (fmpz_is_zero(shift))
        return;

    for (i = 0; i < A->length; i++)
    {
        fmpz_sub(A->exps + i, A->exps + i, shift);
        FLINT_ASSERT(fmpz_sgn(A->exps + i) >= 0);
    }
}

/********* fmpz_mpoly_factor_t ***********************************************/

void fmpz_mpoly_factor_mul_mpoly_fmpz(
    fmpz_mpoly_factor_t f,
    const fmpz_mpoly_t A,
    const fmpz_t e,
    const fmpz_mpoly_ctx_t ctx)
{
    if (fmpz_mpoly_is_fmpz(A, ctx))
    {
        fmpz_t t;
        fmpz_init(t);
        fmpz_mpoly_get_fmpz(t, A, ctx);
        fmpz_pow_fmpz(t, t, e);
        fmpz_mul(f->constant, f->constant, t);
        fmpz_clear(t);
    }
    else
    {
        fmpz_mpoly_factor_append_fmpz(f, A, e, ctx);
    }
}

void fmpz_mpoly_factor_mul_mpoly_ui(
    fmpz_mpoly_factor_t f,
    const fmpz_mpoly_t A,
    ulong e,
    const fmpz_mpoly_ctx_t ctx)
{
    if (fmpz_mpoly_is_fmpz(A, ctx))
    {
        fmpz_t t;
        fmpz_init(t);
        fmpz_mpoly_get_fmpz(t, A, ctx);
        fmpz_pow_ui(t, t, e);
        fmpz_mul(f->constant, f->constant, t);
        fmpz_clear(t);
    }
    else
    {
        fmpz_mpoly_factor_append_ui(f, A, e, ctx);
    }
}

/* a *= b^e */
void fmpz_mpoly_factor_mul_factor_fmpz(
    fmpz_mpoly_factor_t a,
    const fmpz_mpoly_factor_t b,
    const fmpz_t e,
    const fmpz_mpoly_ctx_t ctx)
{
    slong i;
    fmpz_t t;

    fmpz_init(t);

    fmpz_pow_fmpz(t, b->constant, e);
    fmpz_mul(a->constant, a->constant, t);

    for (i = 0; i < b->num; i++)
    {
        fmpz_mul(t, b->exp + i, e);
        fmpz_mpoly_factor_append_fmpz(a, b->poly + i, t, ctx);
    }

    fmpz_clear(t);
}


/******* fmpz_mod_bpoly - (Z/nZ)[x,y] ****************************************/

typedef struct
{
    fmpz_mod_poly_struct * coeffs;
    slong alloc;
    slong length;
    fmpz_t modulus;
} fmpz_mod_bpoly_struct;

typedef fmpz_mod_bpoly_struct fmpz_mod_bpoly_t[1];


/******* fmpz_bpoly - Z[x,y] *************************************************/

typedef struct
{
    fmpz_poly_struct * coeffs;
    slong alloc;
    slong length;
} fmpz_bpoly_struct;

typedef fmpz_bpoly_struct fmpz_bpoly_t[1];


void fmpz_mod_bpoly_init(fmpz_mod_bpoly_t A, const fmpz_t modulus)
{
    A->coeffs = NULL;
    A->alloc = 0;
    A->length = 0;
    fmpz_init_set(A->modulus, modulus);
}

void fmpz_mod_bpoly_clear(fmpz_mod_bpoly_t A)
{
    fmpz_clear(A->modulus);
    if (A->coeffs)
    {
        slong i;
        for (i = 0; i < A->alloc; i++)
            fmpz_mod_poly_clear(A->coeffs + i);
        flint_free(A->coeffs);
    }
}

void fmpz_mod_bpoly_swap(fmpz_mod_bpoly_t A, fmpz_mod_bpoly_t B)
{
    fmpz_mod_bpoly_struct t = *A;
    *A = *B;
    *B = t;
}

void fmpz_bpoly_swap(fmpz_bpoly_t A, fmpz_bpoly_t B)
{
    fmpz_bpoly_struct t = *A;
    *A = *B;
    *B = t;
}

void fmpz_mod_bpoly_print(fmpz_mod_bpoly_t A, const char * xvar, const char * yvar)
{
    slong i;
    int first;

    first = 1;
    for (i = A->length - 1; i >= 0; i--)
    {
        if (fmpz_mod_poly_is_zero(A->coeffs + i))
            continue;

        if (!first)
            flint_printf(" + ");

        first = 0;

        flint_printf("(");
        fmpz_mod_poly_print_pretty(A->coeffs + i, yvar);
        flint_printf(")*%s^%wd", xvar, i);
    }

    if (first)
        flint_printf("0");
}

void fmpz_mod_bpoly_fit_length(fmpz_mod_bpoly_t A, slong len)
{
    slong i;

    if (len <= A->alloc)
        return;

    if (len < 2 * A->alloc)
        len = 2 * A->alloc;

    if (A->alloc == 0)
        A->coeffs = (fmpz_mod_poly_struct *) flint_malloc(len * sizeof(fmpz_mod_poly_struct));
    else
        A->coeffs = (fmpz_mod_poly_struct *) flint_realloc(A->coeffs, len * sizeof(fmpz_mod_poly_struct));

    for (i = A->alloc; i < len; i++)
        fmpz_mod_poly_init(A->coeffs + i, A->modulus);

    A->alloc = len;
}

void fmpz_mod_bpoly_set_fmpz_bpoly(fmpz_mod_bpoly_t A, const fmpz_bpoly_t B)
{
    slong i;
    fmpz_mod_bpoly_fit_length(A, B->length);
    for (i = 0; i < B->length; i++)
    {
        fmpz_mod_poly_set_fmpz_poly(A->coeffs + i, B->coeffs + i);
        if (!fmpz_mod_poly_is_zero(A->coeffs + i))
            A->length = i + 1;
    }
}

void fmpz_mod_bpoly_set_polyx(fmpz_mod_bpoly_t A, const fmpz_mod_poly_t B)
{
    slong i;
    fmpz_mod_bpoly_fit_length(A, B->length);
    A->length = 0;
    for (i = 0; i < B->length; i++)
    {
        fmpz_mod_poly_set_fmpz(A->coeffs + i, B->coeffs + i);
        if (!fmpz_mod_poly_is_zero(A->coeffs + i))
            A->length = i + 1;
    }    
}

void fmpz_mod_bpoly_set_polyy(fmpz_mod_bpoly_t A, const fmpz_mod_poly_t B)
{
    fmpz_mod_bpoly_fit_length(A, 1);
    fmpz_mod_poly_set(A->coeffs + 0, B);
    A->length = !fmpz_mod_poly_is_zero(A->coeffs + 0);
}


void fmpz_mod_bpoly_get_coeff(fmpz_t c, const fmpz_mod_bpoly_t A, slong xi, slong yi)
{
    if (xi >= A->length)
        fmpz_zero(c);

    fmpz_mod_poly_get_coeff_fmpz(c, A->coeffs + xi, yi);
}

void fmpz_mod_bpoly_make_monic(fmpz_mod_bpoly_t A, slong order)
{
    slong i;
    fmpz_mod_poly_t t, lcinv;

    FLINT_ASSERT(A->length > 0);

    fmpz_mod_poly_init(t, A->modulus);
    fmpz_mod_poly_init(lcinv, A->modulus);
    fmpz_mod_poly_inv_series(lcinv, A->coeffs + A->length - 1, order);

    for (i = 0; i < A->length; i++)
    {
        fmpz_mod_poly_mullow(t, A->coeffs + i, lcinv, order);
        fmpz_mod_poly_swap(A->coeffs + i, t);
    }

    fmpz_mod_poly_clear(t);
    fmpz_mod_poly_clear(lcinv);
}

void fmpz_mod_bpoly_mul(fmpz_mod_bpoly_t A, const fmpz_mod_bpoly_t B, const fmpz_mod_bpoly_t C, slong order)
{
    slong i, j;
    fmpz_mod_poly_t t;

    FLINT_ASSERT(A != B);
    FLINT_ASSERT(A != C);

    fmpz_mod_poly_init(t, B->modulus);

    fmpz_mod_bpoly_fit_length(A, B->length + C->length - 1);
    for (i = 0; i < B->length + C->length - 1; i++)
        fmpz_mod_poly_zero(A->coeffs + i);

    for (i = 0; i < B->length; i++)
    {
        for (j = 0; j < C->length; j++)
        {
            fmpz_mod_poly_mullow(t, B->coeffs + i, C->coeffs + j, order);
            fmpz_mod_poly_add(A->coeffs + i + j, A->coeffs + i + j, t);
        }
    }

    A->length = B->length + C->length - 1;

    fmpz_mod_poly_clear(t);
}

void fmpz_mod_bpoly_add_poly_shift(fmpz_mod_bpoly_t A, const fmpz_mod_poly_t B, slong yshift)
{
    slong i;
    fmpz_t c;

    FLINT_ASSERT(A->length > B->length);

    fmpz_init(c);

    for (i = 0; i < B->length; i++)
    {
        fmpz_mod_poly_get_coeff_fmpz(c, A->coeffs + i, yshift);
        FLINT_ASSERT(fmpz_is_zero(c));
        fmpz_mod_poly_set_coeff_fmpz(A->coeffs + i, yshift, B->coeffs + i);
    }

    fmpz_clear(c);
}

void fmpz_mod_bpoly_sub(fmpz_mod_bpoly_t A, const fmpz_mod_bpoly_t B, const fmpz_mod_bpoly_t C)
{
    slong i;

    FLINT_ASSERT(A != B);
    FLINT_ASSERT(A != C);

    fmpz_mod_bpoly_fit_length(A, FLINT_MAX(B->length, C->length));

    for (i = 0; i < FLINT_MAX(B->length, C->length) - 1; i++)
    {
        if (i < B->length)
        {
            if (i < C->length)
            {
                fmpz_mod_poly_sub(A->coeffs + i, B->coeffs + i, C->coeffs + i);
            }
            else
            {
                fmpz_mod_poly_set(A->coeffs + i, B->coeffs + i);
            }
        }
        else
        {
            FLINT_ASSERT(i < C->length);
            fmpz_mod_poly_neg(A->coeffs + i, C->coeffs + i);
        }

        if (!fmpz_mod_poly_is_zero(A->coeffs + i))
            A->length = i + 1;
    }
}

void fmpz_bpoly_init(fmpz_bpoly_t A)
{
    A->coeffs = NULL;
    A->alloc = 0;
    A->length = 0;
}

void fmpz_bpoly_clear(fmpz_bpoly_t A)
{
    if (A->coeffs)
    {
        slong i;
        for (i = 0; i < A->alloc; i++)
            fmpz_poly_clear(A->coeffs + i);
        flint_free(A->coeffs);
    }
}

void fmpz_bpoly_print(fmpz_bpoly_t A, const char * xvar, const char * yvar)
{
    slong i;
    int first = 1;

    for (i = A->length - 1; i >= 0; i--)
    {
        if (fmpz_poly_is_zero(A->coeffs + i))
            continue;

        if (first)
            flint_printf(" + ");
        first = 0;

        flint_printf("(");
        fmpz_poly_print_pretty(A->coeffs + i, yvar);
        flint_printf(")*%s^%wd", xvar, i);
    }

    if (first)
        flint_printf("0");
}

void fmpz_bpoly_fit_length(fmpz_bpoly_t A, slong len)
{
    slong i;

    if (len <= A->alloc)
        return;

    if (len < 2 * A->alloc)
        len = 2 * A->alloc;

    if (A->alloc == 0)
        A->coeffs = (fmpz_poly_struct *) flint_malloc(len * sizeof(fmpz_poly_struct));
    else
        A->coeffs = (fmpz_poly_struct *) flint_realloc(A->coeffs, len * sizeof(fmpz_poly_struct));

    for (i = A->alloc; i < len; i++)
        fmpz_poly_init(A->coeffs + i);

    A->alloc = len;
}

void fmpz_bpoly_set(fmpz_bpoly_t A, const fmpz_bpoly_t B)
{
    slong i;

    FLINT_ASSERT(A != B);

    fmpz_bpoly_fit_length(A, B->length);
    A->length = B->length;

    for (i = 0; i < B->length; i++)
        fmpz_poly_set(A->coeffs + i, B->coeffs + i);
}

void fmpz_bpoly_make_primitive(fmpz_bpoly_t A)
{
    slong i;
    fmpz_poly_t g, q;

    fmpz_poly_init(g);
    fmpz_poly_init(q);

    for (i = 0; i < A->length; i++)
    {
        fmpz_poly_gcd(q, g, A->coeffs + i);
        fmpz_poly_swap(g, q);
    }

    for (i = 0; i < A->length; i++)
    {
        fmpz_poly_div(q, A->coeffs + i, g);
        fmpz_poly_swap(A->coeffs + i, q);
    }

    fmpz_poly_clear(g);
    fmpz_poly_clear(q);
}

void fmpz_bpoly_set_coeff(fmpz_bpoly_t A, slong xi, slong yi, const fmpz_t c)
{
    slong i;

    FLINT_ASSERT(!fmpz_is_zero(c));

    if (xi >= A->length)
    {
        fmpz_bpoly_fit_length(A, xi + 1);
        for (i = A->length; i <= xi; i++)
            fmpz_poly_zero(A->coeffs + i);
        A->length = xi + 1;
    }

    fmpz_poly_set_coeff_fmpz(A->coeffs + xi, yi, c);
}

void fmpz_mod_bpoly_set_coeff(fmpz_mod_bpoly_t A, slong xi, slong yi, const fmpz_t c)
{
    slong i;

    FLINT_ASSERT(!fmpz_is_zero(c));

    if (xi >= A->length)
    {
        fmpz_mod_bpoly_fit_length(A, xi + 1);
        for (i = A->length; i <= xi; i++)
            fmpz_mod_poly_zero(A->coeffs + i);
        A->length = xi + 1;
    }

    fmpz_mod_poly_set_coeff_fmpz(A->coeffs + xi, yi, c);
}


void fmpz_bpoly_zero(fmpz_bpoly_t A)
{
    A->length = 0;
}

void fmpz_mod_bpoly_zero(fmpz_mod_bpoly_t A)
{
    A->length = 0;
}

int fmpz_bpoly_divides(fmpz_bpoly_t Q, fmpz_bpoly_t A, fmpz_bpoly_t B)
{
    slong i, qoff;
    int divides;
    fmpz_poly_t q, t;
    fmpz_bpoly_t R;

    FLINT_ASSERT(Q != A);
    FLINT_ASSERT(Q != B);
    FLINT_ASSERT(B->length > 0);

    fmpz_poly_init(q);
    fmpz_poly_init(t);
    fmpz_bpoly_init(R);
    fmpz_bpoly_set(R, A);

    Q->length = 0;

    while (R->length >= B->length)
    {
        divides = fmpz_poly_divides(q, R->coeffs + R->length - 1, B->coeffs + B->length - 1);
        if (!divides)
            goto cleanup;

        for (i = 0; i < B->length; i++)
        {
            fmpz_poly_mul(t, B->coeffs + i, q);
            fmpz_poly_sub(R->coeffs + i + R->length - B->length, R->coeffs + i + R->length - B->length, t);
        }

        qoff = R->length - B->length;

        FLINT_ASSERT(qoff >= 0);
        if (qoff >= Q->length)
        {
            fmpz_bpoly_fit_length(Q, qoff + 1);
            for (i = Q->length; i <= qoff; i++)
                fmpz_poly_zero(Q->coeffs + i);
            Q->length = qoff + 1;
        }

        fmpz_poly_set(Q->coeffs + qoff, q);

        while (R->length > 0 && fmpz_poly_is_zero(R->coeffs + R->length - 1))
            R->length--;
    }

    divides = (R->length == 0);

cleanup:

    fmpz_poly_clear(q);
    fmpz_poly_clear(t);
    fmpz_bpoly_clear(R);

    return divides;
}


void fmpz_mpoly_from_fmpz_bpoly(
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

void fmpz_bpoly_set_fmpz_mod_bpoly(fmpz_bpoly_t A, const fmpz_mod_bpoly_t B)
{
    slong i;

    fmpz_bpoly_fit_length(A, B->length);
    A->length = B->length;

    for (i = 0; i < B->length; i++)
    {
        fmpz_poly_fit_length(A->coeffs + i, (B->coeffs + i)->length);
        (A->coeffs + i)->length = (B->coeffs + i)->length;
        _fmpz_vec_scalar_smod_fmpz((A->coeffs + i)->coeffs, (B->coeffs + i)->coeffs, (B->coeffs + i)->length, B->modulus);
    }
}

void fmpz_mpoly_to_bpoly(fmpz_bpoly_t A, const fmpz_mpoly_t B, slong varx, slong vary, const fmpz_mpoly_ctx_t ctx)
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

void fmpz_bpoly_eval(fmpz_poly_t E, const fmpz_bpoly_t A, const fmpz_t alpha)
{
    slong i;
    fmpz_t t;

    fmpz_init(t);

    fmpz_poly_zero(E);
    for (i = A->length - 1; i >= 0; i--)
    {
        fmpz_poly_evaluate_fmpz(t, A->coeffs + i, alpha);
        fmpz_poly_set_coeff_fmpz(E, i, t);
    }

    fmpz_clear(t);
}

void fmpz_bpoly_taylor_shift(fmpz_bpoly_t A, const fmpz_t alpha)
{
    slong i;
    for (i = A->length - 1; i >= 0; i--)
        _fmpz_poly_taylor_shift((A->coeffs + i)->coeffs, alpha, (A->coeffs + i)->length);    
}

typedef struct {
    slong r; /* number of local factors */
    ulong k;
    slong lifting_prec;
    fmpz_t p;
    fmpz_t pk;
    fmpz_mod_bpoly_t Btilde;                /* mod p^k */
    fmpz_mod_bpoly_struct * newBitilde;     /* mod p^k */
    fmpz_mod_poly_struct * P;               /* mod p^k */
    fmpz_mod_poly_struct * d;               /* mod p^k */
    fmpz_mod_poly_struct * Bitilde;         /* mod p^k */
    fmpz_mod_poly_struct * d1;              /* mod p */
    fmpz_mod_poly_struct * Bitilde1;        /* mod p */
} bpoly_info_struct;

typedef bpoly_info_struct bpoly_info_t[1];

void bpoly_info_init(bpoly_info_t I, slong r, const fmpz_t p, ulong k)
{
    slong i;

    FLINT_ASSERT(r >= 2);

    I->r = r;

    I->lifting_prec = 0;

    I->k = k;
    fmpz_init_set(I->p, p);
    fmpz_init(I->pk);
    fmpz_pow_ui(I->pk, p, k);

    fmpz_mod_bpoly_init(I->Btilde, I->pk);

    I->newBitilde = (fmpz_mod_bpoly_struct *) flint_malloc(I->r*sizeof(fmpz_mod_bpoly_struct));
    I->P          = (fmpz_mod_poly_struct *) flint_malloc(I->r*sizeof(fmpz_mod_poly_struct));
    I->d          = (fmpz_mod_poly_struct *) flint_malloc(I->r*sizeof(fmpz_mod_poly_struct));
    I->Bitilde    = (fmpz_mod_poly_struct *) flint_malloc(I->r*sizeof(fmpz_mod_poly_struct));
    I->d1         = (fmpz_mod_poly_struct *) flint_malloc(I->r*sizeof(fmpz_mod_poly_struct));
    I->Bitilde1   = (fmpz_mod_poly_struct *) flint_malloc(I->r*sizeof(fmpz_mod_poly_struct));

    for (i = 0; i < I->r; i++)
    {
        fmpz_mod_bpoly_init(I->newBitilde + i, I->pk);
        fmpz_mod_poly_init(I->P + i, I->pk);
        fmpz_mod_poly_init(I->d + i, I->pk);
        fmpz_mod_poly_init(I->Bitilde + i, I->pk);
        fmpz_mod_poly_init(I->d1 + i, I->p);
        fmpz_mod_poly_init(I->Bitilde1 + i, I->p);
    }
}

void bpoly_info_clear(bpoly_info_t I)
{
    slong i;

    fmpz_clear(I->p);
    fmpz_clear(I->pk);

    fmpz_mod_bpoly_clear(I->Btilde);

    for (i = 0; i < I->r; i++)
    {
        fmpz_mod_bpoly_clear(I->newBitilde + i);
        fmpz_mod_poly_clear(I->P + i);
        fmpz_mod_poly_clear(I->d + i);
        fmpz_mod_poly_clear(I->Bitilde + i);
        fmpz_mod_poly_clear(I->d1 + i);
        fmpz_mod_poly_clear(I->Bitilde1 + i);
    }

    flint_free(I->newBitilde);
    flint_free(I->P);
    flint_free(I->d);
    flint_free(I->Bitilde);
    flint_free(I->d1);
    flint_free(I->Bitilde1);
}

/*
    set out[i] so that
    1/(f[0]*f[1]*...*f[n-1]) = out[0]/f[0] + ... + out[n-1]/f[n-1]
*/
int partial_fraction_coeffs(fmpz_mod_poly_struct * out, const fmpz_mod_poly_struct * f, const fmpz_t p, slong n)
{
    slong i;
    fmpz_mod_poly_t num, den, a, b, g, t;

    FLINT_ASSERT(n > 1);

    fmpz_mod_poly_init(num, p);
    fmpz_mod_poly_init(den, p);
    fmpz_mod_poly_init(a, p);
    fmpz_mod_poly_init(b, p);
    fmpz_mod_poly_init(g, p);
    fmpz_mod_poly_init(t, p);

    fmpz_mod_poly_set_ui(num, 1);
    fmpz_mod_poly_mul(den, f + 0, f + 1);
    for (i = 2; i < n; i++)
        fmpz_mod_poly_mul(den, den, f + i);

    for (i = 0; i < n; i++)
    {
        fmpz_mod_poly_divrem(den, t, den, f + i);
        FLINT_ASSERT(fmpz_mod_poly_is_zero(t));
        fmpz_mod_poly_xgcd(g, a, b, f + i, den);
        if (fmpz_mod_poly_degree(g) != 0)
            return 0;
        FLINT_ASSERT(fmpz_is_one(g->coeffs + 0));
        fmpz_mod_poly_mul(t, b, num);
        fmpz_mod_poly_rem(out + i, t, f + i);
        fmpz_mod_poly_mul(t, a, num);
        fmpz_mod_poly_rem(num, t, den);
    }

    fmpz_mod_poly_clear(num);
    fmpz_mod_poly_clear(den);
    fmpz_mod_poly_clear(a);
    fmpz_mod_poly_clear(b);
    fmpz_mod_poly_clear(g);
    fmpz_mod_poly_clear(t);
    return 1;
}


int bpoly_info_disolve(bpoly_info_t I)
{
    slong i, j;
    fmpz_t pj, t1;
    fmpz_mod_poly_t error, t, s, s1, s2;

    if (!partial_fraction_coeffs(I->d1, I->Bitilde1, I->p, I->r))
        return 0;

    fmpz_init(pj);
    fmpz_init(t1);
    fmpz_mod_poly_init(error, I->pk);
    fmpz_mod_poly_init(t, I->pk);
    fmpz_mod_poly_init(s, I->p);
    fmpz_mod_poly_init(s1, I->p);
    fmpz_mod_poly_init(s2, I->pk);

    for (i = 0; i < I->r; i++)
    {
        fmpz_mod_poly_set_ui(I->P + i, 1);
        for (j = 0; j < I->r; j++)
        {
            if (i == j)
                continue;
            fmpz_mod_poly_mul(I->P + i, I->P + i, I->Bitilde + j);
        }
    }

    fmpz_mod_poly_set_ui(error, 1);
    for (i = 0; i < I->r; i++)
    {
        fmpz_mod_poly_set(I->d + i, I->d1 + i); /* slight abuse because moduli are different */
        fmpz_mod_poly_mul(t, I->d + i, I->P + i);
        fmpz_mod_poly_sub(error, error, t);
    }

    fmpz_one(pj);
    for (j = 1; j < I->k; j++)
    {
        fmpz_mul(pj, pj, I->p);
        fmpz_mod_poly_zero(s);
        for (i = error->length - 1; i >= 0; i--)
        {
            FLINT_ASSERT(fmpz_divisible(error->coeffs + i, pj));
            fmpz_divexact(t1, error->coeffs + i, pj);
            fmpz_mod(t1, t1, I->p);
            fmpz_mod_poly_set_coeff_fmpz(s, i, t1);
        }

        for (i = 0; i < I->r; i++)
        {
            fmpz_mod_poly_mul(s1, s, I->d1 + i);
            fmpz_set(&s2->p, I->p);
            fmpz_mod_poly_rem(s2, s1, I->Bitilde1 + i);
            fmpz_set(&s2->p, I->pk);
            fmpz_mod_poly_scalar_mul_fmpz(s2, s2, pj);
            fmpz_mod_poly_add(I->d + i, I->d + i, s2);
        }

        fmpz_mod_poly_set_ui(error, 1);
        for (i = 0; i < I->r; i++)
        {
            fmpz_mod_poly_mul(t, I->d + i, I->P + i);
            fmpz_mod_poly_sub(error, error, t);
        }
    }

    FLINT_ASSERT(fmpz_mod_poly_is_zero(error));

    fmpz_clear(pj);
    fmpz_clear(t1);
    fmpz_mod_poly_clear(error);
    fmpz_mod_poly_clear(t);
    fmpz_mod_poly_clear(s);
    fmpz_mod_poly_clear(s1);
    fmpz_mod_poly_clear(s2);

    return 1;
}



static void _bivar_lift_quintic(bpoly_info_t I)
{
    slong i, j, k;
    fmpz_mod_bpoly_t tp, tp1, error;
    fmpz_mod_poly_t ss, tt;
/*
timeit_t timer;

timeit_start(timer);
*/
    fmpz_mod_poly_init(ss, I->pk);
    fmpz_mod_poly_init(tt, I->pk);
    fmpz_mod_bpoly_init(tp, I->pk);
    fmpz_mod_bpoly_init(tp1, I->pk);
    fmpz_mod_bpoly_init(error, I->pk);

    fmpz_mod_bpoly_mul(tp, I->newBitilde + 0, I->newBitilde + 1, I->lifting_prec);
    for (i = 2; i < I->r; i++)
    {
        fmpz_mod_bpoly_mul(tp1, tp, I->newBitilde + i, I->lifting_prec);
        fmpz_mod_bpoly_swap(tp1, tp);
    }
    fmpz_mod_bpoly_sub(error, I->Btilde, tp);

    for (j = 1; j < I->lifting_prec; j++)
    {
        fmpz_mod_poly_zero(ss);
        for (i = error->length - 1; i >= 0; i--)
        {
            fmpz_t ct;
            fmpz_init(ct);

            fmpz_mod_bpoly_get_coeff(ct, error, i, j);
            fmpz_mod_poly_set_coeff_fmpz(ss, i, ct);
            for (k = 0; k < j; k++)
            {
                fmpz_mod_bpoly_get_coeff(ct, error, i, k);
                FLINT_ASSERT(fmpz_is_zero(ct));
            }

            fmpz_clear(ct);
        }

        for (i = 0; i < I->r; i++)
        {
            fmpz_mod_poly_mul(tt, ss, I->d + i);
            fmpz_mod_poly_rem(tt, tt, I->Bitilde + i);
            fmpz_mod_bpoly_add_poly_shift(I->newBitilde + i, tt, j);
        }

        fmpz_mod_bpoly_mul(tp, I->newBitilde + 0, I->newBitilde + 1, I->lifting_prec);
        for (i = 2; i < I->r; i++)
        {
            fmpz_mod_bpoly_mul(tp1, tp, I->newBitilde + i, I->lifting_prec);
            fmpz_mod_bpoly_swap(tp1, tp);
        }
        fmpz_mod_bpoly_sub(error, I->Btilde, tp);
    }
/*
flint_printf("------------------\n");
for (k = 0; k < I->r; k++)
{
flint_printf("newBitilde[%wd]: ", k); fmpz_mod_bpoly_print(I->newBitilde + k, "x", "y"); printf("\n");
}

timeit_stop(timer);
flint_printf("_bivar_lift_quintic time: %wd\n", timer->wall);
*/
    fmpz_mod_poly_clear(ss);
    fmpz_mod_poly_clear(tt);
    fmpz_mod_bpoly_clear(tp);
    fmpz_mod_bpoly_clear(tp1);
    fmpz_mod_bpoly_clear(error);
}

void fmpz_mod_bpoly_reverse_vars(
    fmpz_mod_bpoly_t A,
    const fmpz_mod_bpoly_t B)
{
    slong i, j;
    fmpz_mod_bpoly_zero(A);
    for (i = 0; i < B->length; i++)
    {
        fmpz_mod_poly_struct * Bi = B->coeffs + i;
        for (j = 0; j < Bi->length; j++)
        {
            if (!fmpz_is_zero(Bi->coeffs + j))
            {
                fmpz_mod_bpoly_set_coeff(A, j, i, Bi->coeffs + j);
            }
        }
    }
}

static void _bivar_lift_quartic2(bpoly_info_t I)
{
    slong i, j, k;
    fmpz_mod_poly_t t, t1;
    fmpz_mod_bpoly_t btilde;
    fmpz_mod_bpoly_struct newbitilde[2];
/*
timeit_t timer;
timeit_start(timer);
*/
    FLINT_ASSERT(I->r == 2);

    fmpz_mod_poly_init(t, I->pk);
    fmpz_mod_poly_init(t1, I->pk);
    fmpz_mod_bpoly_init(btilde, I->pk);
    fmpz_mod_bpoly_reverse_vars(btilde, I->Btilde);
    for (k = 0; k < I->r; k++)
    {
        fmpz_mod_bpoly_init(newbitilde + k, I->pk);
        fmpz_mod_bpoly_reverse_vars(newbitilde + k, I->newBitilde + k);
        fmpz_mod_bpoly_fit_length(newbitilde + k, I->lifting_prec);
        FLINT_ASSERT((newbitilde + k)->length == 1);
    }

    for (j = 1; j < I->lifting_prec; j++)
    {
        if (j < btilde->length)
            fmpz_mod_poly_set(t, btilde->coeffs + j);
        else
            fmpz_mod_poly_zero(t);

        for (i = 1; i < j; i++)
        {
            fmpz_mod_poly_mul(t1, newbitilde[0].coeffs + i, newbitilde[1].coeffs + j - i);
            fmpz_mod_poly_sub(t, t, t1);
        }

        for (k = 0; k < I->r; k++)
        {
            fmpz_mod_poly_mul(t1, t, I->d + k);
            fmpz_mod_poly_rem(newbitilde[k].coeffs + j, t1, I->Bitilde + k);
            if (!fmpz_mod_poly_is_zero(newbitilde[k].coeffs + j))
                newbitilde[k].length = j + 1;
        }
    }

    for (k = 0; k < I->r; k++)
        fmpz_mod_bpoly_reverse_vars(I->newBitilde + k, newbitilde + k);
/*
flint_printf("------------------\n");
for (k = 0; k < I->r; k++)
{
flint_printf("newBitilde[%wd]: ", k); fmpz_mod_bpoly_print(I->newBitilde + k, "x", "y"); printf("\n");
}
timeit_stop(timer);
flint_printf("_bivar_lift_quartic2 time: %wd\n", timer->wall);
*/

    fmpz_mod_poly_clear(t);
    fmpz_mod_poly_clear(t1);
    fmpz_mod_bpoly_clear(btilde);
    for (k = 0; k < I->r; k++)
    {
        fmpz_mod_bpoly_clear(newbitilde + k);
    }
}

static void _bivar_lift_quartic(bpoly_info_t I)
{
    slong i, j, k;
    fmpz_mod_poly_t t, t1;
    fmpz_mod_bpoly_t btilde;
    fmpz_mod_bpoly_struct * newbitilde, * U;

    FLINT_ASSERT(I->r > 2);

    newbitilde = (fmpz_mod_bpoly_struct *) flint_malloc(I->r*sizeof(fmpz_mod_bpoly_struct));
    U = (fmpz_mod_bpoly_struct *) flint_malloc(I->r*sizeof(fmpz_mod_bpoly_struct));

    fmpz_mod_poly_init(t, I->pk);
    fmpz_mod_poly_init(t1, I->pk);
    fmpz_mod_bpoly_init(btilde, I->pk);
    fmpz_mod_bpoly_reverse_vars(btilde, I->Btilde);
    for (k = 0; k < I->r; k++)
    {
        fmpz_mod_bpoly_init(U + k, I->pk);
        fmpz_mod_bpoly_fit_length(U + k, I->lifting_prec);
        for (i = 0; i < I->lifting_prec; i++)
        {
            fmpz_mod_poly_zero(U[k].coeffs + i);
        }

        fmpz_mod_bpoly_init(newbitilde + k, I->pk);
        fmpz_mod_bpoly_reverse_vars(newbitilde + k, I->newBitilde + k);
        fmpz_mod_bpoly_fit_length(newbitilde + k, I->lifting_prec);
        FLINT_ASSERT(newbitilde[k].length == 1);
        for (i = 1; i < I->lifting_prec; i++)
        {
            fmpz_mod_poly_zero(newbitilde[k].coeffs + i);
        }
    }

    k = I->r - 2;
    fmpz_mod_poly_mul(U[k].coeffs + 0, newbitilde[k].coeffs + 0, newbitilde[k + 1].coeffs + 0);
    for (k--; k >= 1; k--)
        fmpz_mod_poly_mul(U[k].coeffs + 0, newbitilde[k].coeffs + 0, U[k + 1].coeffs + 0);

    for (j = 1; j < I->lifting_prec; j++)
    {
        k = I->r - 2;
        fmpz_mod_poly_zero(U[k].coeffs + j);
        for (i = 0; i <= j; i++)
        {
            fmpz_mod_poly_mul(t1, newbitilde[k].coeffs + i, newbitilde[k + 1].coeffs + j - i);
            fmpz_mod_poly_add(U[k].coeffs + j, U[k].coeffs + j, t1);
        }
        for (k--; k >= 1; k--)
        {
            fmpz_mod_poly_zero(U[k].coeffs + j);
            for (i = 0; i <= j; i++)
            {
                fmpz_mod_poly_mul(t1, newbitilde[k].coeffs + i, U[k + 1].coeffs + j - i);
                fmpz_mod_poly_add(U[k].coeffs + j, U[k].coeffs + j, t1);
            }
        }

        if (j < btilde->length)
            fmpz_mod_poly_set(t, btilde->coeffs + j);
        else
            fmpz_mod_poly_zero(t);

        for (i = 0; i <= j; i++)
        {
            fmpz_mod_poly_mul(t1, newbitilde[0].coeffs + i, U[1].coeffs + j - i);
            fmpz_mod_poly_sub(t, t, t1);
        }

        for (k = 0; k < I->r; k++)
        {
            fmpz_mod_poly_mul(t1, t, I->d + k);
            fmpz_mod_poly_rem(newbitilde[k].coeffs + j, t1, I->Bitilde + k);
            if (!fmpz_mod_poly_is_zero(newbitilde[k].coeffs + j))
                newbitilde[k].length = j + 1;
        }

        k = I->r - 2;
        fmpz_mod_poly_mul(t, newbitilde[k].coeffs + 0, newbitilde[k + 1].coeffs + j);
        fmpz_mod_poly_mul(t1, newbitilde[k].coeffs + j, newbitilde[k + 1].coeffs + 0);
        fmpz_mod_poly_add(t, t, t1);
        fmpz_mod_poly_add(U[k].coeffs + j, U[k].coeffs + j, t);
        for (k--; k >= 1; k--)
        {
            fmpz_mod_poly_mul(t1, newbitilde[k].coeffs + 0, t);
            fmpz_mod_poly_swap(t, t1);
            fmpz_mod_poly_mul(t1, newbitilde[k].coeffs + j, U[k + 1].coeffs + 0);
            fmpz_mod_poly_add(t, t, t1);
            fmpz_mod_poly_add(U[k].coeffs + j, U[k].coeffs + j, t);
        }
    }

    for (k = 0; k < I->r; k++)
        fmpz_mod_bpoly_reverse_vars(I->newBitilde + k, newbitilde + k);

    fmpz_mod_poly_clear(t);
    fmpz_mod_poly_clear(t1);
    fmpz_mod_bpoly_clear(btilde);
    for (k = 0; k < I->r; k++)
    {
        fmpz_mod_bpoly_clear(U + k);
        fmpz_mod_bpoly_clear(newbitilde + k);
    }

    flint_free(newbitilde);
    flint_free(U);
}


static void _recombine_naive(
    fmpz_mpoly_factor_t fac,
    flint_bitcnt_t bits,
    slong xvar,
    slong yvar,
    const fmpz_mpoly_ctx_t ctx,
    fmpz_bpoly_t B,
    fmpz_t alpha,
    bpoly_info_t I)
{
    fmpz_bpoly_t Q, R, trymez;
    fmpz_mod_bpoly_t tryme, trymet;
    fmpz_mod_poly_t leadB;
    slong i, j, r, len;
    fmpz_t subset, tsubset, test;
    fmpz_mpoly_t goodtry;
    slong * idx;

    fmpz_init(test);
    fmpz_init(tsubset);
    fmpz_init(subset);
    fmpz_mpoly_init(goodtry, ctx);

    fmpz_bpoly_init(Q);
    fmpz_bpoly_init(R);
    fmpz_bpoly_init(trymez);
    fmpz_mod_bpoly_init(tryme, I->pk);
    fmpz_mod_bpoly_init(trymet, I->pk);
    fmpz_mod_poly_init(leadB, I->pk);

    fmpz_mpoly_factor_one(fac, ctx);

    FLINT_ASSERT(B->length > 0);
    fmpz_mod_poly_set_fmpz_poly(leadB, B->coeffs + B->length - 1);

    len = I->r;
    idx = (slong *) flint_malloc(I->r * sizeof(slong));
    for (i = 0; i < len; i++)
        idx[i] = i;

    for (r = 1; r <= len/2; r++)
    {
        subset_first(subset, len, r);
        do {
try_subset:
            fmpz_mod_bpoly_set_polyy(tryme, leadB);

            for (i = 0; i < len; i++)
            {
                if (fmpz_tstbit(subset, i))
                {
                    fmpz_mod_bpoly_mul(trymet, tryme, I->newBitilde + idx[i], I->lifting_prec);
                    fmpz_mod_bpoly_swap(trymet, tryme);
                }
            }
            fmpz_bpoly_set_fmpz_mod_bpoly(trymez, tryme);
            fmpz_bpoly_make_primitive(trymez);

            if (fmpz_bpoly_divides(Q, B, trymez))
            {
                fmpz_neg(alpha, alpha);
                fmpz_bpoly_taylor_shift(trymez, alpha);
                fmpz_neg(alpha, alpha);
                fmpz_mpoly_from_fmpz_bpoly(goodtry, bits, trymez, xvar, yvar, ctx);
                fmpz_mpoly_factor_append_ui(fac, goodtry, 1, ctx);
                fmpz_bpoly_swap(B, Q);
                FLINT_ASSERT(B->length > 0);
                fmpz_mod_poly_set_fmpz_poly(leadB, B->coeffs + B->length - 1);

                /* fix indices */
                j = 0;
                for (i = 0; i < len; i++)
                {
                    if (!fmpz_tstbit(subset, i))
                        idx[j++] = idx[i];
                }
                len -= r;

                /* fix subsets */
                fmpz_set(tsubset, subset);
                do {
                    if (!subset_next(tsubset, tsubset, len + r))
                        goto rloop_continue;
                    fmpz_and(test, tsubset, subset);
                } while (!fmpz_is_zero(test));
                subset_map_down(test, tsubset, subset);
                fmpz_swap(test, subset);
                goto try_subset;
            }
        }
        while (subset_next(subset, subset, len));
rloop_continue:
        (void)(NULL);
    }

    fmpz_neg(alpha, alpha);
    fmpz_bpoly_taylor_shift(B, alpha);
    fmpz_neg(alpha, alpha);
    fmpz_mpoly_from_fmpz_bpoly(goodtry, bits, B, xvar, yvar, ctx);
    if (B->length > 1)
    {
        fmpz_mpoly_factor_append_ui(fac, goodtry, 1, ctx);
    }
    else
    {
        FLINT_ASSERT(fac->num > 0);
        fmpz_mpoly_mul(fac->poly + fac->num - 1,
                       fac->poly + fac->num - 1, goodtry, ctx);
    }

    fmpz_bpoly_clear(Q);
    fmpz_bpoly_clear(R);
    fmpz_bpoly_clear(trymez);
    fmpz_mod_bpoly_clear(tryme);
    fmpz_mod_bpoly_clear(trymet);
    fmpz_mod_poly_clear(leadB);

    fmpz_mpoly_clear(goodtry, ctx);
    fmpz_clear(subset);
    fmpz_clear(tsubset);
    fmpz_clear(test);
    flint_free(idx);
}


static int _irreducible_bivar_factors(
    fmpz_mpoly_factor_t fac,
    const fmpz_mpoly_t A,
    slong xvar,
    slong yvar,
    const fmpz_mpoly_ctx_t ctx)
{
    int success;
    slong i;
    fmpz_t alpha;
    fmpz_poly_t Beval;
    fmpz_bpoly_t B;
    fmpz_poly_factor_t Bevalfac;
    slong Blengthx, Blengthy;
    flint_bitcnt_t Bbits;
    ulong pkbits;
    ulong k;
    fmpz_t p;
    bpoly_info_t I;

    fmpz_mpoly_factor_one(fac, ctx);

    k = 1;
    fmpz_init_set_ui(p, UWORD(1) << (FLINT_BITS - 1));
    fmpz_init(alpha);
    fmpz_poly_init(Beval);
    fmpz_poly_factor_init(Bevalfac);
    fmpz_bpoly_init(B);
    bpoly_info_init(I, 2, p, k);

    fmpz_mpoly_to_bpoly(B, A, xvar, yvar, ctx);

    Blengthx = B->length;
    FLINT_ASSERT(Blengthx > 1);

    fmpz_zero(alpha);
    goto got_alpha;

next_alpha:

    fmpz_neg(alpha, alpha);
    fmpz_add_ui(alpha, alpha, fmpz_sgn(alpha) >= 0);

got_alpha:

    fmpz_bpoly_eval(Beval, B, alpha);

    /* if killed leading coeff, get new alpha */
    if (Beval->length != Blengthx)
        goto next_alpha;

    fmpz_one(&Bevalfac->c);
    Bevalfac->num = 0;
    fmpz_poly_factor(Bevalfac, Beval);

    /* if multiple factors, get new alpha */
    for (i = 0; i < Bevalfac->num; i++)
    {
        if (Bevalfac->exp[i] != 1)
            goto next_alpha;
    }

    /* if one factor, A is irreducible */
    if (Bevalfac->num == 1)
    {
        fmpz_mpoly_factor_append_ui(fac, A, 1, ctx);
        success = 1;
        goto cleanup;
    }

    fmpz_bpoly_taylor_shift(B, alpha);

    Blengthy = 0;
    Bbits = 0;
    for (i = 0; i < B->length; i++)
    {
        slong this_bits;
        Blengthy = FLINT_MAX(Blengthy, (B->coeffs + i)->length);
        this_bits = _fmpz_vec_max_bits((B->coeffs + i)->coeffs,
                                       (B->coeffs + i)->length);
        Bbits = FLINT_MAX(Bbits, FLINT_ABS(this_bits));
    }

    pkbits = (FLINT_BIT_COUNT(Blengthx*Blengthy) + 1)/2;
    pkbits += Blengthx + Blengthy + Bbits - 3;

next_prime:

    fmpz_nextprime(p, p, 1);

    FLINT_ASSERT(B->length > 0);
    FLINT_ASSERT((B->coeffs + B->length - 1)->length > 0);
    FLINT_ASSERT(!fmpz_is_zero((B->coeffs + B->length - 1)->coeffs + 0));

    if (fmpz_divisible((B->coeffs + B->length - 1)->coeffs + 0, p))
        goto next_prime;

    k = (pkbits + fmpz_bits(p))/fmpz_bits(p);

    bpoly_info_clear(I);
    bpoly_info_init(I, Bevalfac->num, p, k);
    I->lifting_prec = Blengthy + (B->coeffs + B->length - 1)->length;

    fmpz_mod_bpoly_set_fmpz_bpoly(I->Btilde, B);
    fmpz_mod_bpoly_make_monic(I->Btilde, I->lifting_prec);
    for (i = 0; i < I->r; i++)
    {
        fmpz_mod_poly_set_fmpz_poly(I->Bitilde1 + i, Bevalfac->p + i);
        fmpz_mod_poly_make_monic(I->Bitilde1 + i, I->Bitilde1 + i);
        fmpz_mod_poly_set_fmpz_poly(I->Bitilde + i, Bevalfac->p + i);
        fmpz_mod_poly_make_monic(I->Bitilde + i, I->Bitilde + i);
        fmpz_mod_bpoly_set_polyx(I->newBitilde + i, I->Bitilde + i);
    }

    FLINT_ASSERT(I->r > 1);

    if (!bpoly_info_disolve(I))
        goto next_prime;

/*
flint_printf("I->pk: "); fmpz_print(I->pk); printf("\n");
    for (k = 0; k < I->r; k++)
    {
flint_printf("newBitilde[%wd]: ", k); fmpz_mod_bpoly_print(I->newBitilde + k, "x", "y"); printf("\n");
    }
*/
    if (I->r == 2)
        _bivar_lift_quartic2(I);
    else if (I->r < 20)
        _bivar_lift_quartic(I);
    else
        _bivar_lift_quintic(I);
/*
printf("---------- lifted factors ----------\n");
    for (k = 0; k < I->r; k++)
    {
flint_printf("newBitilde[%wd]: ", k); fmpz_mod_bpoly_print(I->newBitilde + k, "x", "y"); printf("\n");
    }
*/

    _recombine_naive(fac, A->bits, xvar, yvar, ctx, B, alpha, I);
    success = 1;

cleanup:

    bpoly_info_clear(I);
    fmpz_bpoly_clear(B);
    fmpz_poly_factor_clear(Bevalfac);
    fmpz_poly_clear(Beval);
    fmpz_clear(alpha);
    fmpz_clear(p);

    FLINT_ASSERT(!success || fmpz_mpoly_factor_matches(A, fac, ctx));

    return success;
}





void fmpz_mpoly_convert_perm(
    fmpz_mpoly_t A,
    flint_bitcnt_t Abits,
    const fmpz_mpoly_ctx_t lctx,
    const fmpz_mpoly_t B,
    const fmpz_mpoly_ctx_t ctx,
    const slong * perm)
{
    slong n = ctx->minfo->nvars;
    slong m = lctx->minfo->nvars;
    slong i, k, l;
    slong NA, NB;
    ulong * Aexps;
    ulong * Bexps;
    TMP_INIT;

    FLINT_ASSERT(B->length > 0);
    FLINT_ASSERT(Abits <= FLINT_BITS);
    FLINT_ASSERT(B->bits <= FLINT_BITS);
    TMP_START;

    Aexps = (ulong *) TMP_ALLOC((m)*sizeof(ulong));
    Bexps = (ulong *) TMP_ALLOC(n*sizeof(ulong));

    NA = mpoly_words_per_exp(Abits, lctx->minfo);
    NB = mpoly_words_per_exp(B->bits, ctx->minfo);

    fmpz_mpoly_fit_bits(A, Abits, ctx);
    A->bits = Abits;
    fmpz_mpoly_fit_length(A, B->length, lctx);
    A->length = B->length;
    for (i = 0; i < B->length; i++)
    {        
        fmpz_set(A->coeffs + i, B->coeffs + i);
        mpoly_get_monomial_ui(Bexps, B->exps + NB*i, B->bits, ctx->minfo);
        for (k = 0; k < m; k++)
        {
            l = perm[k];
            Aexps[k] = l < 0 ? 0 : Bexps[l];
        }
        mpoly_set_monomial_ui(A->exps + NA*(i), Aexps, Abits, lctx->minfo);
     }  
    fmpz_mpoly_sort_terms(A, lctx);
    TMP_END;
}


static void _to_poly(fmpz_poly_t A, const fmpz_mpoly_t B, const fmpz_mpoly_ctx_t ctx)
{
    slong i;
    ulong * shift, * stride;

    if (B->length == 0)
    {
        fmpz_poly_zero(A);
        return;
    }

    shift = (ulong *) flint_malloc(ctx->minfo->nvars*sizeof(ulong));
    stride = (ulong *) flint_malloc(ctx->minfo->nvars*sizeof(ulong));
    for (i = 0; i < ctx->minfo->nvars; i++)
    {
        shift[i] = 0;
        stride[i] = 1;
    }

    _fmpz_mpoly_to_fmpz_poly_deflate(A, B, 0, shift, stride, ctx);

    flint_free(shift);
    flint_free(stride);
}

static void _to_polyq(
    fmpq_poly_t A,
    const fmpz_mpoly_t B,
    const fmpz_mpoly_ctx_t ctx)
{
    slong i;
    ulong mask;
    slong shift, off, N;
    slong Blen = B->length;
    fmpz * Bcoeff = B->coeffs;
    ulong * Bexp = B->exps;

    FLINT_ASSERT(B->bits <= FLINT_BITS);

    fmpq_poly_zero(A);

    N = mpoly_words_per_exp_sp(B->bits, ctx->minfo);
    mpoly_gen_offset_shift_sp(&off, &shift, 0, B->bits, ctx->minfo);

    mask = (-UWORD(1)) >> (FLINT_BITS - B->bits);
    for (i = 0; i < Blen; i++)
        fmpq_poly_set_coeff_fmpz(A, (Bexp[N*i + off] >> shift) & mask, Bcoeff + i);
}

static int _from_polyq(
    fmpz_mpoly_t A,
    flint_bitcnt_t Abits,
    const fmpq_poly_t B,
    const fmpz_mpoly_ctx_t ctx)
{
    slong var = 0;
    slong N;
    slong k;
    slong Alen;
    fmpz * Acoeff;
    ulong * Aexp;
    slong Aalloc;
    ulong * strideexp;
    TMP_INIT;

    if (B->length == 0)
    {
        fmpz_mpoly_zero(A, ctx);
        return 1;
    }

    if (!fmpz_is_one(fmpq_poly_denref(B)))
    {
        return 0;
    }

    TMP_START;

    FLINT_ASSERT(Abits <= FLINT_BITS);

    N = mpoly_words_per_exp(Abits, ctx->minfo);
    strideexp = (ulong*) TMP_ALLOC(N*sizeof(ulong));
    mpoly_gen_monomial_sp(strideexp, var, Abits, ctx->minfo);

    fmpz_mpoly_fit_bits(A, Abits, ctx);
    A->bits = Abits;

    Acoeff = A->coeffs;
    Aexp = A->exps;
    Aalloc = A->alloc;
    Alen = 0;
    for (k = B->length - 1; k >= 0; k--)
    {
        _fmpz_mpoly_fit_length(&Acoeff, &Aexp, &Aalloc, Alen + 1, N);
        if (!fmpz_is_zero(B->coeffs + k))
        {
            fmpz_swap(Acoeff + Alen, B->coeffs + k);
            mpoly_monomial_mul_ui(Aexp + N*Alen, strideexp, N, k);
            Alen++;
        }
    }
    A->coeffs = Acoeff;
    A->exps = Aexp;
    A->alloc = Aalloc;
    _fmpz_mpoly_set_length(A, Alen, ctx);

    TMP_END;

    return 1;
}


void fmpz_mpoly_to_mpolyv(
    fmpz_mpolyv_t A,
    const fmpz_mpoly_t B,
    const fmpz_mpoly_t xalpha,
    const fmpz_mpoly_ctx_t ctx)
{
    fmpz_mpoly_t Q, T;

    fmpz_mpoly_init(Q, ctx);
    fmpz_mpoly_init(T, ctx);

    fmpz_mpolyv_fit_length(A, 8, ctx);
    fmpz_mpoly_divrem(Q, A->coeffs + 0, B, xalpha, ctx);
    A->length = 1;

    while (!fmpz_mpoly_is_zero(Q, ctx))
    {
        fmpz_mpolyv_fit_length(A, A->length + 1, ctx);
        fmpz_mpoly_divrem(T, A->coeffs + A->length, Q, xalpha, ctx);
        fmpz_mpoly_swap(Q, T, ctx);
        A->length++;
    }

    while (A->length > 0 && fmpz_mpoly_is_zero(A->coeffs + A->length - 1, ctx))
        A->length--;

    fmpz_mpoly_clear(Q, ctx);
    fmpz_mpoly_clear(T, ctx);
}

void fmpz_mpoly_from_mpolyv(
    fmpz_mpoly_t A,
    const fmpz_mpolyv_t B,
    const fmpz_mpoly_t xalpha,
    const fmpz_mpoly_ctx_t ctx)
{
    slong i;
    fmpz_mpoly_t T;

    fmpz_mpoly_init(T, ctx);

    fmpz_mpoly_zero(A, ctx);
    for (i = B->length - 1; i >= 0; i--)
    {
        fmpz_mpoly_mul(T, A, xalpha, ctx);
        fmpz_mpoly_add(A, T, B->coeffs + i, ctx);
    }

    fmpz_mpoly_clear(T, ctx);
}


typedef struct {
    slong n;
    slong w;
    slong r;
    fmpq_poly_struct * inv_prod_dbetas;
    fmpq_poly_struct * dbetas;
    fmpz_mpoly_struct * prod_mbetas;
    fmpz_mpolyv_struct * prod_mbetas_coeffs;
    fmpz_mpoly_struct * mbetas;
    fmpz_mpoly_struct * deltas;

    fmpz_mpoly_struct * xalpha;
    fmpz_mpoly_struct * q;
    fmpz_mpoly_struct * qt;
    fmpz_mpoly_struct * newt;
    fmpz_mpolyv_struct * delta_coeffs;

    fmpq_poly_t dtq, S, R;

} fmpz_mpoly_disolve_struct;

typedef fmpz_mpoly_disolve_struct fmpz_mpoly_disolve_t[1];


int fmpz_mpoly_disolve_init(
    fmpz_mpoly_disolve_t I,
    slong r, slong w,
    const fmpz_mpoly_struct * betas,
    const fmpz * alpha,
    const fmpz_mpoly_ctx_t ctx)
{
    slong success = 1;
    slong i, j, k;
    fmpz_poly_t p;
    fmpq_poly_t G, S, pq;
/*
flint_printf("fmpz_mpoly_disolve_init called(l = %wd, w = %wd)\n",r,w);
*/
    I->r = r;
    I->w = w;

    I->dbetas = (fmpq_poly_struct *) flint_malloc(r*sizeof(fmpq_poly_struct));
    I->inv_prod_dbetas = (fmpq_poly_struct *) flint_malloc(
                                                   r*sizeof(fmpq_poly_struct));
    I->prod_mbetas = (fmpz_mpoly_struct *) flint_malloc(
                                          (w + 1)*r*sizeof(fmpz_mpoly_struct));
    I->prod_mbetas_coeffs = (fmpz_mpolyv_struct *) flint_malloc(
                                         (w + 1)*r*sizeof(fmpz_mpolyv_struct));
    I->mbetas = (fmpz_mpoly_struct *) flint_malloc(
                                          (w + 1)*r*sizeof(fmpz_mpoly_struct));
    I->deltas = (fmpz_mpoly_struct *) flint_malloc(
                                          (w + 1)*r*sizeof(fmpz_mpoly_struct));
    fmpq_poly_init(I->dtq);
    fmpq_poly_init(I->S);
    fmpq_poly_init(I->R);
    I->xalpha = (fmpz_mpoly_struct *) flint_malloc(
                                            (w + 1)*sizeof(fmpz_mpoly_struct));
    I->q = (fmpz_mpoly_struct *) flint_malloc(
                                            (w + 1)*sizeof(fmpz_mpoly_struct));
    I->qt = (fmpz_mpoly_struct *) flint_malloc(
                                            (w + 1)*sizeof(fmpz_mpoly_struct));
    I->newt = (fmpz_mpoly_struct *) flint_malloc(
                                            (w + 1)*sizeof(fmpz_mpoly_struct));
    I->delta_coeffs = (fmpz_mpolyv_struct *) flint_malloc(
                                         (w + 1)*r*sizeof(fmpz_mpolyv_struct));
    for (i = 0; i <= w; i++)
    {
        fmpz_mpoly_init(I->xalpha + i, ctx);
        fmpz_mpoly_init(I->q + i, ctx);
        fmpz_mpoly_init(I->qt + i, ctx);
        fmpz_mpoly_init(I->newt + i, ctx);
        for (j = 0; j < r; j++)
            fmpz_mpolyv_init(I->delta_coeffs + i*I->r + j, ctx);

        if (i < 1)
            continue;

        fmpz_mpoly_gen(I->xalpha + i, i, ctx);
        fmpz_mpoly_sub_fmpz(I->xalpha + i, I->xalpha + i, alpha + i - 1, ctx);
    }

    fmpz_poly_init(p);
    fmpq_poly_init(G);
    fmpq_poly_init(S);
    fmpq_poly_init(pq);

    /* initialize deltas */
    for (i = w; i >= 0; i--)
        for (j = 0; j < r; j++)
            fmpz_mpoly_init(I->deltas + i*r + j, ctx);

    /* set betas */
    i = w;
    for (j = 0; j < r; j++)
    {
        fmpz_mpoly_init(I->mbetas + i*r + j, ctx);
        fmpz_mpoly_set(I->mbetas + i*r + j, betas + j, ctx);
    }
    for (i--; i >= 0; i--)
    {
        for (j = 0; j < r; j++)
        {
            fmpz_mpoly_init(I->mbetas + i*r + j, ctx);
            fmpz_mpoly_evaluate_one_fmpz(I->mbetas + i*r + j,
                             I->mbetas + (i + 1)*r + j, i + 1, alpha + i, ctx);
        }
    }
    for (j = 0; j < r; j++)
    {
        _to_poly(p, I->mbetas + 0*r + j, ctx);
        fmpq_poly_init(I->dbetas + j);
        fmpq_poly_set_fmpz_poly(I->dbetas + j, p);
    }

    /* set product of betas */
    for (i = w; i >= 0; i--)
    {
        for (j = 0; j < r; j++)
        {
            fmpz_mpoly_init(I->prod_mbetas + i*r + j, ctx);
            fmpz_mpoly_one(I->prod_mbetas + i*r + j, ctx);
            for (k = 0; k < r; k++)
            {
                if (k == j)
                    continue;
                fmpz_mpoly_mul(I->prod_mbetas + i*r + j,
                           I->prod_mbetas + i*r + j, I->mbetas + i*r + k, ctx);
            }
            fmpz_mpolyv_init(I->prod_mbetas_coeffs + i*r + j, ctx);
            if (i > 0)
            {
                fmpz_mpoly_to_mpolyv(I->prod_mbetas_coeffs + i*r + j,
                                I->prod_mbetas + i*r + j, I->xalpha + i, ctx);
            }
        }        
    }

    for (j = 0; j < r; j++)
        fmpq_poly_init(I->inv_prod_dbetas + j);

    for (j = 0; success && j < r; j++)
    {
        if (fmpq_poly_degree(I->dbetas + j) !=
                 fmpz_mpoly_degree_si(betas + j, 0, ctx))
        {
            success = 0;
        }
    }

    for (j = 0; success && j < r; j++)
    {
        fmpq_poly_one(pq);
        for (k = 0; k < r; k++)
        {
            if (k == j)
                continue;
            fmpq_poly_mul(pq, pq, I->dbetas + k);
        }
        fmpq_poly_init(I->inv_prod_dbetas + j);
        fmpq_poly_xgcd(G, S, I->inv_prod_dbetas + j, I->dbetas + j, pq);
        if (!fmpq_poly_is_one(G))
            success = 0;
    }

    fmpz_poly_clear(p);
    fmpq_poly_clear(G);
    fmpq_poly_clear(S);
    fmpq_poly_clear(pq);

    FLINT_ASSERT(success == 1);

    return success;
}


static void fmpz_mpoly_disolve_clear(
    fmpz_mpoly_disolve_t I,
    const fmpz_mpoly_ctx_t ctx)
{
    slong i, j;

    fmpq_poly_clear(I->dtq);
    fmpq_poly_clear(I->S);
    fmpq_poly_clear(I->R);

    for (i = 0; i <= I->w; i++)
    {
        fmpz_mpoly_clear(I->xalpha + i, ctx);
        fmpz_mpoly_clear(I->q + i, ctx);
        fmpz_mpoly_clear(I->qt + i, ctx);
        fmpz_mpoly_clear(I->newt + i, ctx);
        for (j = 0; j < I->r; j++)
            fmpz_mpolyv_clear(I->delta_coeffs + i*I->r + j, ctx);
    }
    flint_free(I->xalpha);
    flint_free(I->q);
    flint_free(I->qt);
    flint_free(I->newt);
    flint_free(I->delta_coeffs);

    for (j = 0; j < I->r; j++)
    {
        fmpq_poly_clear(I->inv_prod_dbetas + j);
        fmpq_poly_clear(I->dbetas + j);
        for (i = 0; i <= I->w; i++)
        {
            fmpz_mpolyv_clear(I->prod_mbetas_coeffs + i*I->r + j, ctx);
            fmpz_mpoly_clear(I->prod_mbetas + i*I->r + j, ctx);
            fmpz_mpoly_clear(I->mbetas + i*I->r + j, ctx);
            fmpz_mpoly_clear(I->deltas + i*I->r + j, ctx);
        }
    }

    flint_free(I->inv_prod_dbetas);
    flint_free(I->dbetas);
    flint_free(I->prod_mbetas);
    flint_free(I->prod_mbetas_coeffs);
    flint_free(I->mbetas);
    flint_free(I->deltas);
}


static int fmpz_mpoly_disolve(
    flint_bitcnt_t bits,
    slong l,
    fmpz_mpoly_t t,
    const slong * degs,
    fmpz_mpoly_disolve_t I,
    const fmpz_mpoly_ctx_t ctx)
{
    slong i, j, k;
    int success;
    fmpz_mpoly_struct * deltas = I->deltas + l*I->r;
    fmpz_mpoly_struct * newdeltas = I->deltas + (l - 1)*I->r;
    fmpz_mpoly_struct * q = I->q + l;
    fmpz_mpoly_struct * qt = I->qt + l;
    fmpz_mpoly_struct * newt = I->newt + l;
    fmpz_mpolyv_struct * delta_coeffs = I->delta_coeffs + l*I->r;

    FLINT_ASSERT(bits <= FLINT_BITS);

    if (t->bits > FLINT_BITS)
        return -1;

/*
flint_printf("_mfactor_disolve(r = %wd, I->r = %wd) called:\n", r, I->r);
flint_printf("t: "); fmpz_mpoly_print_pretty(t, NULL, ctx); printf("\n");
for (i = 0; i < I->r; i++)
{
flint_printf("dbetas[%wd]: ", i);
fmpq_poly_print_pretty(I->dbetas + i, "x1");
flint_printf("\n");
}
*/

    if (l < 1)
    {
        _to_polyq(I->dtq, t, ctx);

        success = 1;
        for (i = 0; i < I->r; i++)
        {
            fmpq_poly_mul(I->S, I->dtq, I->inv_prod_dbetas + i);
            fmpq_poly_rem(I->R, I->S, I->dbetas + i);
            if (!_from_polyq(deltas + i, bits, I->R, ctx))
                return 0;
        }
        return 1;
    }

    for (i = 0; i < I->r; i++)
        delta_coeffs[i].length = 0;

    for (k = 0; k <= degs[l]; k++)
    {
        fmpz_mpoly_divrem(q, newt, t, I->xalpha + l, ctx);
        fmpz_mpoly_swap(t, q, ctx);
        for (j = 0; j < k; j++)
        for (i = 0; i < I->r; i++)
        {
            if (j >= delta_coeffs[i].length ||
                k - j >= I->prod_mbetas_coeffs[l*I->r + i].length)
            {
                continue;
            }

            fmpz_mpoly_mul(qt, delta_coeffs[i].coeffs + j,
                        I->prod_mbetas_coeffs[l*I->r + i].coeffs + k - j, ctx);
            fmpz_mpoly_sub(q, newt, qt, ctx);
            fmpz_mpoly_swap(newt, q, ctx);
        }

        success = fmpz_mpoly_disolve(bits, l - 1, newt, degs, I, ctx);
        if (success <= 0)
            return success;

        for (i = 0; i < I->r; i++)
        {
            if (fmpz_mpoly_is_zero(newdeltas + i, ctx))
                continue;

            if (k + I->prod_mbetas_coeffs[l*I->r + i].length - 1 > degs[l])
                return 0;

            fmpz_mpolyv_set_coeff(delta_coeffs + i, k, newdeltas + i, ctx);
        }
    }

    for (i = 0; i < I->r; i++)
        fmpz_mpoly_from_mpolyv(deltas + i, delta_coeffs + i, I->xalpha + l, ctx);

    return 1;
}



static int _mfactor_lift_quartic2(
    slong m,
    fmpz_mpoly_factor_t lfac,
    const fmpz * alpha,
    const fmpz_mpoly_t A,
    const slong * degs,
    const fmpz_mpoly_ctx_t ctx)
{
    int success;
    slong i, j;
    slong r = 2;
    fmpz_mpoly_t Aq, t, t2, t3, xalpha;
    fmpz_mpoly_struct * betas, * deltas;
    fmpz_mpoly_disolve_t I;
    fmpz_mpolyv_struct B[2];
    slong tdeg;
/*
flint_printf("_mfactor_lift_quartic2 called (degs[%wd] = %wd, alpha = %wd)\n", m, degs[m], *alpha);
flint_printf("lfac: "); fmpz_mpoly_factor_print_pretty(lfac, NULL, ctx); printf("\n");
flint_printf("A: "); fmpz_mpoly_print_pretty(A, NULL, ctx); printf("\n");
*/
    FLINT_ASSERT(r == lfac->num);

    fmpz_mpoly_init(t, ctx);
    fmpz_mpoly_init(t2, ctx);
    fmpz_mpoly_init(t3, ctx);
    fmpz_mpoly_init(xalpha, ctx);
    fmpz_mpoly_init(Aq, ctx);

    fmpz_mpoly_gen(xalpha, m, ctx);
    fmpz_mpoly_sub_fmpz(xalpha, xalpha, alpha + m - 1, ctx);

    betas  = (fmpz_mpoly_struct * ) flint_malloc(r*sizeof(fmpz_mpoly_struct));
    for (i = 0; i < r; i++)
    {
        fmpz_mpolyv_init(B + i, ctx);
        fmpz_mpoly_to_mpolyv(B + i, lfac->poly + i, xalpha, ctx);
        fmpz_mpolyv_fit_length(B + i, degs[m] + 1, ctx);
        for (j = B[i].length; j <= degs[m]; j++)
            fmpz_mpoly_zero(B[i].coeffs + j, ctx);
        betas[i] = B[i].coeffs[0];
    }

    success = fmpz_mpoly_disolve_init(I, lfac->num, m - 1, betas, alpha, ctx);
    FLINT_ASSERT(success == 1);

    deltas = I->deltas + (m - 1)*I->r;

    fmpz_mpoly_divrem(t2, t, A, xalpha, ctx);
    fmpz_mpoly_swap(Aq, t2, ctx);
#if WANT_ASSERT
    fmpz_mpoly_one(t2, ctx);
    for (i = 0; i < r; i++)
        fmpz_mpoly_mul(t2, t2, betas + i, ctx);
    FLINT_ASSERT(fmpz_mpoly_equal(t, t2, ctx));
#endif

    for (j = 1; j <= degs[m]; j++)
    {
        fmpz_mpoly_divrem(t2, t, Aq, xalpha, ctx);
        fmpz_mpoly_swap(Aq, t2, ctx);

        for (i = 0; i <= j; i++)
        {
            fmpz_mpoly_mul(t2, B[0].coeffs + i, B[1].coeffs + j - i, ctx);
            fmpz_mpoly_sub(t3, t, t2, ctx);
            fmpz_mpoly_swap(t, t3, ctx);
        }

        success = fmpz_mpoly_disolve(A->bits, m - 1, t, degs, I, ctx);
        if (success <= 0)
        {
            success = 0;
            goto cleanup;
        }

        tdeg = 0;
        for (i = 0; i < r; i++)
        {
            fmpz_mpoly_add(t3, B[i].coeffs + j, deltas + i, ctx);
            fmpz_mpoly_swap(B[i].coeffs + j, t3, ctx);
            if (!fmpz_mpoly_is_zero(B[i].coeffs + j, ctx))
                B[i].length = FLINT_MAX(B[i].length, j + 1);
            FLINT_ASSERT(B[i].length > 0);
            tdeg += B[i].length - 1;
        }

        if (tdeg > degs[m])
        {
            success = 0;
            goto cleanup;
        }
    }

    success = 1;

cleanup:

    fmpz_mpoly_disolve_clear(I, ctx);

    flint_free(betas);

    for (i = 0; i < r; i++)
    {
        if (success)
            fmpz_mpoly_from_mpolyv(lfac->poly + i, B + i, xalpha, ctx);
        fmpz_mpolyv_clear(B + i, ctx);
    }

    fmpz_mpoly_clear(t, ctx);
    fmpz_mpoly_clear(t2, ctx);
    fmpz_mpoly_clear(t3, ctx);
    fmpz_mpoly_clear(xalpha, ctx);
    fmpz_mpoly_clear(Aq, ctx);

    return success;
}

static int _mfactor_lift_quartic(
    slong m,
    fmpz_mpoly_factor_t lfac,
    const fmpz * alpha,
    const fmpz_mpoly_t A,
    const slong * degs,
    const fmpz_mpoly_ctx_t ctx)
{
    int success;
    slong i, j, k;
    slong r = lfac->num;
    fmpz_mpoly_t t, t1, t2, t3, xalpha;
    fmpz_mpoly_struct * betas, * deltas;
    fmpz_mpoly_disolve_t I;
    fmpz_mpolyv_t Av;
    fmpz_mpolyv_struct * B, * U;
    slong tdeg;
/*
flint_printf("_mfactor_lift_quartic called (degs[%wd] = %wd, alpha = %wd)\n", m, degs[m], *alpha);
flint_printf("lfac: "); fmpz_mpoly_factor_print_pretty(lfac, NULL, ctx); printf("\n");
flint_printf("A: "); fmpz_mpoly_print_pretty(A, NULL, ctx); printf("\n");
*/
    FLINT_ASSERT(r > 2);

    B = (fmpz_mpolyv_struct *) flint_malloc(r*sizeof(fmpz_mpolyv_struct));
    U = (fmpz_mpolyv_struct *) flint_malloc(r*sizeof(fmpz_mpolyv_struct));

    fmpz_mpoly_init(t, ctx);
    fmpz_mpoly_init(t1, ctx);
    fmpz_mpoly_init(t2, ctx);
    fmpz_mpoly_init(t3, ctx);
    fmpz_mpoly_init(xalpha, ctx);

    fmpz_mpoly_gen(xalpha, m, ctx);
    fmpz_mpoly_sub_fmpz(xalpha, xalpha, alpha + m - 1, ctx);

    fmpz_mpolyv_init(Av, ctx);
    fmpz_mpoly_to_mpolyv(Av, A, xalpha, ctx);
    fmpz_mpolyv_fit_length(Av, degs[m] + 1, ctx);
    for (j = Av->length; j <= degs[m]; j++)
        fmpz_mpoly_zero(Av->coeffs + j, ctx);

    for (k = 0; k < r; k++)
    {
        fmpz_mpolyv_init(U + k, ctx);
        fmpz_mpolyv_fit_length(U + k, degs[m] + 1, ctx);
        for (j = 0; j <= degs[m]; j++)
            fmpz_mpoly_zero(U[k].coeffs + j, ctx);

        fmpz_mpolyv_init(B + k, ctx);
        fmpz_mpoly_to_mpolyv(B + k, lfac->poly + k, xalpha, ctx);
        fmpz_mpolyv_fit_length(B + k, degs[m] + 1, ctx);
        for (j = Av->length; j <= degs[m]; j++)
            fmpz_mpoly_zero(B[k].coeffs + j, ctx);
    }

    betas  = (fmpz_mpoly_struct * ) flint_malloc(r*sizeof(fmpz_mpoly_struct));
    for (i = 0; i < r; i++)
        betas[i] = B[i].coeffs[0];

    fmpz_mpoly_disolve_init(I, lfac->num, m - 1, betas, alpha, ctx);
    deltas = I->deltas + (m - 1)*I->r;

    k = r - 2;
    fmpz_mpoly_mul(U[k].coeffs + 0, B[k].coeffs + 0, B[k + 1].coeffs + 0, ctx);
    for (k--; k >= 1; k--)
        fmpz_mpoly_mul(U[k].coeffs + 0, B[k].coeffs + 0, U[k + 1].coeffs + 0, ctx);

    for (j = 1; j <= degs[m]; j++)
    {
        k = r - 2;
        fmpz_mpoly_zero(U[k].coeffs + j, ctx);
        for (i = 0; i <= j; i++)
        {
            fmpz_mpoly_mul(t1, B[k].coeffs + i, B[k + 1].coeffs + j - i, ctx);
            fmpz_mpoly_add(U[k].coeffs + j, U[k].coeffs + j, t1, ctx);

        }
        for (k--; k >= 1; k--)
        {
            fmpz_mpoly_zero(U[k].coeffs + j, ctx);
            for (i = 0; i <= j; i++)
            {
                fmpz_mpoly_mul(t1, B[k].coeffs + i, U[k + 1].coeffs + j - i, ctx);
                fmpz_mpoly_add(U[k].coeffs + j, U[k].coeffs + j, t1, ctx);
            }
        }

        if (j < Av->length)
            fmpz_mpoly_set(t, Av->coeffs + j, ctx);
        else
            fmpz_mpoly_zero(t, ctx);

        for (i = 0; i <= j; i++)
        {
            fmpz_mpoly_mul(t2, B[0].coeffs + i, U[1].coeffs + j - i, ctx);
            fmpz_mpoly_sub(t3, t, t2, ctx);
            fmpz_mpoly_swap(t, t3, ctx);
        }

        if (fmpz_mpoly_is_zero(t, ctx))
            continue;

        success = fmpz_mpoly_disolve(A->bits, m - 1, t, degs, I, ctx);
        if (success <= 0)
        {
            success = 0;
            goto cleanup;
        }

        tdeg = 0;
        for (i = 0; i < r; i++)
        {
            fmpz_mpoly_add(t3, B[i].coeffs + j, deltas + i, ctx);
            fmpz_mpoly_swap(B[i].coeffs + j, t3, ctx);
            if (!fmpz_mpoly_is_zero(B[i].coeffs + j, ctx))
                B[i].length = FLINT_MAX(B[i].length, j + 1);
            FLINT_ASSERT(B[i].length > 0);
            tdeg += B[i].length - 1;
        }

        if (tdeg > degs[m])
        {
            success = 0;
            goto cleanup;
        }

        k = r - 2;
        fmpz_mpoly_mul(t, B[k].coeffs + 0, deltas + k + 1, ctx);
        fmpz_mpoly_mul(t1, deltas + k, B[k + 1].coeffs + 0, ctx);
        fmpz_mpoly_add(t, t, t1, ctx);
        fmpz_mpoly_add(U[k].coeffs + j, U[k].coeffs + j, t, ctx);
        for (k--; k >= 1; k--)
        {
            fmpz_mpoly_mul(t1, B[k].coeffs + 0, t, ctx);
            fmpz_mpoly_swap(t, t1, ctx);
            fmpz_mpoly_mul(t1, deltas + k, U[k + 1].coeffs + 0, ctx);
            fmpz_mpoly_add(t, t, t1, ctx);
            fmpz_mpoly_add(U[k].coeffs + j, U[k].coeffs + j, t, ctx);
        }

    }

    success = 1;

cleanup:

    fmpz_mpoly_disolve_clear(I, ctx);

    flint_free(betas);

    fmpz_mpolyv_clear(Av, ctx);
    for (i = 0; i < r; i++)
    {
        if (success)
            fmpz_mpoly_from_mpolyv(lfac->poly + i, B + i, xalpha, ctx);
        fmpz_mpolyv_clear(B + i, ctx);
        fmpz_mpolyv_clear(U + i, ctx);
    }

    flint_free(B);
    flint_free(U);
    fmpz_mpoly_clear(t, ctx);
    fmpz_mpoly_clear(t1, ctx);
    fmpz_mpoly_clear(t2, ctx);
    fmpz_mpoly_clear(t3, ctx);
    fmpz_mpoly_clear(xalpha, ctx);

    return success;
}


static int _mfactor_lift_quintic(
    slong m,
    fmpz_mpoly_factor_t lfac,
    const fmpz * alpha,
    const fmpz_mpoly_t A,
    const slong * degs,
    const fmpz_mpoly_ctx_t ctx)
{
    int success;
    slong i, j;
    slong r = lfac->num;
    fmpz_mpoly_t e, t, pow, xalpha, q;
    fmpz_mpoly_struct * betas, * deltas;
    fmpz_mpoly_disolve_t I;
/*
flint_printf("_mfactor_lift_quintic called (m = %wd, alpha = %wd)\n", m, *alpha);
flint_printf("lfac: "); fmpz_mpoly_factor_print_pretty(lfac, NULL, ctx); printf("\n");
flint_printf("A: "); fmpz_mpoly_print_pretty(A, NULL, ctx); printf("\n");
*/
    FLINT_ASSERT(r > 1);

    fmpz_mpoly_init(e, ctx);
    fmpz_mpoly_init(t, ctx);
    fmpz_mpoly_init(pow, ctx);
    fmpz_mpoly_init(xalpha, ctx);
    fmpz_mpoly_init(q, ctx);

    betas  = (fmpz_mpoly_struct * ) flint_malloc(r*sizeof(fmpz_mpoly_struct));
    for (i = 0; i < r; i++)
    {
        fmpz_mpoly_init(betas + i, ctx);
        fmpz_mpoly_evaluate_one_fmpz(betas + i, lfac->poly + i, m, alpha + m - 1, ctx);
    }

    fmpz_mpoly_mul(t, lfac->poly + 0, lfac->poly + 1, ctx);
    for (i = 2; i < r; i++)
        fmpz_mpoly_mul(t, t, lfac->poly + i, ctx);
    fmpz_mpoly_sub(e, A, t, ctx);

    fmpz_mpoly_one(pow, ctx);

    fmpz_mpoly_gen(xalpha, m, ctx);
    fmpz_mpoly_sub_fmpz(xalpha, xalpha, alpha + m - 1, ctx);

    fmpz_mpoly_disolve_init(I, lfac->num, m - 1, betas, alpha, ctx);
    deltas = I->deltas + (m - 1)*I->r;

    for (j = 1; j <= degs[m]; j++)
    {
        if (fmpz_mpoly_is_zero(e, ctx))
        {
            success = 1;
            goto cleanup;
        }

        fmpz_mpoly_mul(pow, pow, xalpha, ctx);
        success = fmpz_mpoly_divides(q, e, pow, ctx);
        FLINT_ASSERT(success);
        fmpz_mpoly_evaluate_one_fmpz(t, q, m, alpha + m - 1, ctx);

        success = fmpz_mpoly_disolve(A->bits, m - 1, t, degs, I, ctx);
        if (success <= 0)
        {
            success = 0;
            goto cleanup;
        }

        for (i = 0; i < r; i++)
        {
            fmpz_mpoly_mul(t, deltas + i, pow, ctx);
            fmpz_mpoly_add(lfac->poly + i, lfac->poly + i, t, ctx);
        }

        fmpz_mpoly_mul(t, lfac->poly + 0, lfac->poly + 1, ctx);
        for (i = 2; i < r; i++)
            fmpz_mpoly_mul(t, t, lfac->poly + i, ctx);
        fmpz_mpoly_sub(e, A, t, ctx);
    }

    success = fmpz_mpoly_is_zero(e, ctx);

cleanup:

    fmpz_mpoly_disolve_clear(I, ctx);

    fmpz_mpoly_clear(e, ctx);
    fmpz_mpoly_clear(t, ctx);
    fmpz_mpoly_clear(pow, ctx);
    fmpz_mpoly_clear(xalpha, ctx);
    fmpz_mpoly_clear(q, ctx);

    for (i = 0; i < r; i++)
        fmpz_mpoly_clear(betas + i, ctx);

    flint_free(betas);

    return success;
}



static void _fmpz_mpoly_get_lc(
    fmpz_mpoly_t c,
    const fmpz_mpoly_t A,
    const fmpz_mpoly_ctx_t ctx)
{
    slong dummyvars[] = {0};
    ulong dummydegs[] = {0};

    dummyvars[0] = 0;
    dummydegs[0] = fmpz_mpoly_degree_si(A, 0, ctx);
    fmpz_mpoly_get_coeff_vars_ui(c, A, dummyvars, dummydegs, 1, ctx);
}

static void _nmod_mpoly_get_lc(
    nmod_mpoly_t c,
    const nmod_mpoly_t A,
    const nmod_mpoly_ctx_t ctx)
{
    slong dummyvars[] = {0};
    ulong dummydegs[] = {0};

    dummyvars[0] = 0;
    dummydegs[0] = nmod_mpoly_degree_si(A, 0, ctx);
    nmod_mpoly_get_coeff_vars_ui(c, A, dummyvars, dummydegs, 1, ctx);
}

static void _fmpz_mpoly_set_lc(
    fmpz_mpoly_t A,
    const fmpz_mpoly_t B,
    const fmpz_mpoly_t c,
    const fmpz_mpoly_ctx_t ctx)
{
    slong deg;
    fmpz_mpoly_t t, g;

    fmpz_mpoly_init(t, ctx);
    fmpz_mpoly_init(g, ctx);

    deg = fmpz_mpoly_degree_si(B, 0, ctx);
    FLINT_ASSERT(deg >= 0);
    fmpz_mpoly_gen(g, 0, ctx);
    fmpz_mpoly_pow_ui(g, g, deg, ctx);
    _fmpz_mpoly_get_lc(t, B, ctx);
    fmpz_mpoly_sub(t, c, t, ctx);
    fmpz_mpoly_mul(t, t, g, ctx);
    fmpz_mpoly_add(A, B, t, ctx);

    fmpz_mpoly_clear(t, ctx);
    fmpz_mpoly_clear(g, ctx);
}

static void _nmod_mpoly_set_lc(
    nmod_mpoly_t A,
    const nmod_mpoly_t B,
    const nmod_mpoly_t c,
    const nmod_mpoly_ctx_t ctx)
{
    slong deg;
    nmod_mpoly_t t, g;

    nmod_mpoly_init(t, ctx);
    nmod_mpoly_init(g, ctx);

    deg = nmod_mpoly_degree_si(B, 0, ctx);
    FLINT_ASSERT(deg >= 0);
    nmod_mpoly_gen(g, 0, ctx);
    nmod_mpoly_pow_ui(g, g, deg, ctx);
    _nmod_mpoly_get_lc(t, B, ctx);
    nmod_mpoly_sub(t, c, t, ctx);
    nmod_mpoly_mul(t, t, g, ctx);
    nmod_mpoly_add(A, B, t, ctx);

    nmod_mpoly_clear(t, ctx);
    nmod_mpoly_clear(g, ctx);
}



/*
have
[pfac1]*[pfac2]*[pfac3] = p
is a factorization of primitivepart(q(x_m = alpha_m))

set
lcq = lc(q)
lcp = lc(q(x_m = alpha_m)), not nec the same as lc(p)

change to
[lcp/lc(pfac1)*pfac1]*[lcp/lc(pfac2)*pfac2]*[lcp/lc(pfac3)*pfac3] = lcp^2 p

now each lcp/lc(pfaci)*pfaci has leading coefficient lcp

replace the leading coefficient of each lcp/lc(pfaci)*pfaci with lcq

lift against lcq^2*q

remove content
*/
static int _try_lift(
    fmpz_mpoly_factor_t qfac,
    const fmpz_mpoly_t q,
    const fmpz_mpoly_factor_t pfac,
    const fmpz_mpoly_t p,
    slong m,
    fmpz * alpha,
    slong n,
    const fmpz_mpoly_ctx_t ctx)
{
    int success;
    slong i;
    slong * newdeg;
    fmpz_mpoly_t lcq, lcp, t, newq;
    fmpz_mpoly_univar_t u;

    FLINT_ASSERT(pfac->num > 1);

    newdeg = (slong *) flint_malloc((n + 1)*sizeof(slong));
    fmpz_mpoly_init(lcq, ctx);
    fmpz_mpoly_init(lcp, ctx);
    fmpz_mpoly_init(t, ctx);
    fmpz_mpoly_init(newq, ctx);
    fmpz_mpoly_univar_init(u, ctx);

    FLINT_ASSERT(fmpz_is_one(pfac->constant));
    FLINT_ASSERT(fmpz_mpoly_factor_matches(p, pfac, ctx));

    _fmpz_mpoly_get_lc(lcq, q, ctx);
    fmpz_mpoly_evaluate_one_fmpz(lcp, lcq, m, alpha + m - 1, ctx);

    FLINT_ASSERT(lcp->length > 0);

    fmpz_mpoly_pow_ui(t, lcq, pfac->num - 1, ctx);
    fmpz_mpoly_mul(newq, q, t, ctx);
    fmpz_mpoly_degrees_si(newdeg, newq, ctx);

    fmpz_set(qfac->constant, pfac->constant);
    fmpz_mpoly_factor_fit_length(qfac, pfac->num, ctx);
    qfac->num = pfac->num;
    for (i = 0; i < pfac->num; i++)
    {
        _fmpz_mpoly_get_lc(t, pfac->poly + i, ctx);
        success = fmpz_mpoly_divides(t, lcp, t, ctx);
        FLINT_ASSERT(success);
        fmpz_mpoly_mul(qfac->poly + i, pfac->poly + i, t, ctx);
        _fmpz_mpoly_set_lc(qfac->poly + i, qfac->poly + i, lcq, ctx);
        fmpz_one(qfac->exp + i);
    }

    if (qfac->num == 2)
        success = _mfactor_lift_quartic2(m, qfac, alpha, newq, newdeg, ctx);
    else if (qfac->num <  20)
        success = _mfactor_lift_quartic(m, qfac, alpha, newq, newdeg, ctx);
    else
        success = _mfactor_lift_quintic(m, qfac, alpha, newq, newdeg, ctx);

    if (!success)
        goto cleanup;

    for (i = 0; i < qfac->num; i++)
    {
        fmpz_mpoly_to_univar(u, qfac->poly + i, 0, ctx);
        success = fmpz_mpoly_univar_content_mpoly(t, u, ctx);
        if (!success)
            goto cleanup;
        success = fmpz_mpoly_divides(qfac->poly + i, qfac->poly + i, t, ctx);
        FLINT_ASSERT(success);
        FLINT_ASSERT(qfac->poly[i].length > 0);
        if (fmpz_sgn(qfac->poly[i].coeffs + 0) < 0)
            fmpz_mpoly_neg(qfac->poly + i, qfac->poly + i, ctx);
    }

    success = 1;

cleanup:

    flint_free(newdeg);
    fmpz_mpoly_clear(lcq, ctx);
    fmpz_mpoly_clear(lcp, ctx);
    fmpz_mpoly_clear(t, ctx);
    fmpz_mpoly_clear(newq, ctx);
    fmpz_mpoly_univar_clear(u, ctx);

    /* q and its factors are primitive with positive lc */
    FLINT_ASSERT(!success || (fmpz_is_one(qfac->constant) &&
                              fmpz_mpoly_factor_matches(q, qfac, ctx)));
    return success;
}


void fmpz_poly_factor_print_pretty(const fmpz_poly_factor_t f)
{
    slong i;
    fmpz_print(&f->c);
    for (i = 0; i < f->num; i++)
    {
        flint_printf("*(");
        fmpz_poly_print_pretty(f->p + i, "x");
        flint_printf(")^%wd", f->exp[i]);
    }
}



static int _wang_lcc(
    fmpz_mpolyv_t lc_divs,
    const fmpz_mpoly_factor_t lcAfac,
    const fmpz_poly_factor_t Aufac,
    const fmpz * alpha,
    const fmpz_mpoly_ctx_t ctx)
{
    int success;
    slong i, j, k;
    const slong n = ctx->minfo->nvars - 1;
    fmpz * lcAfaceval = _fmpz_vec_init(lcAfac->num);
    const fmpz ** salpha = (const fmpz **) flint_malloc((n + 1)*sizeof(fmpz *));
    fmpz * d = _fmpz_vec_init(1 + lcAfac->num);
    fmpz zero = 0;
    fmpz * dtilde = _fmpz_vec_init(Aufac->num);
    fmpz_t delta, q, r;
    fmpz_mpoly_t t;

    fmpz_init(delta);
    fmpz_init(q);
    fmpz_init(r);

    fmpz_mpoly_init(t, ctx);

    salpha[0] = &zero;
    for (i = 0; i < n; i++)
        salpha[i + 1] = alpha + i;

    for (j = 0; j < lcAfac->num; j++)
        fmpz_mpoly_evaluate_all_fmpz(lcAfaceval + j, lcAfac->poly + j,
                                                 (fmpz * const *) salpha, ctx);

    fmpz_mul(d + 0, &Aufac->c, lcAfac->constant);
    for (i = 0; i < lcAfac->num; i++)
    {
        fmpz_abs(q, lcAfaceval + i);
        if (fmpz_is_one(q))
        {
            success = 0;
            goto cleanup;
        }
        for (j = i; j >= 0; j--)
        {
            fmpz_set(r, d + j);
            while (!fmpz_is_one(r))
            {
                fmpz_gcd(r, r, q);
                fmpz_divexact(q, q, r);
                if (fmpz_is_one(q))
                {
                    success = 0;
                    goto cleanup;
                }
            }
        }
        fmpz_set(d + i + 1, q);
    }

    fmpz_mpolyv_fit_length(lc_divs, Aufac->num, ctx);
    lc_divs->length = Aufac->num;

    for (j = 0; j < Aufac->num; j++)
    {
        fmpz_mpoly_one(lc_divs->coeffs + j, ctx);
        fmpz_mul(r, Aufac->p[j].coeffs + Aufac->p[j].length - 1, &Aufac->c);
        fmpz_one(dtilde + j);
        for (i = lcAfac->num - 1; i >= 0; i--)
        {
            fmpz_abs(q, lcAfaceval + i);
            if (fmpz_cmp_ui(q, 2) < 0)
                continue;
            k = fmpz_remove(r, r, q);
            fmpz_mpoly_pow_ui(t, lcAfac->poly + i, k, ctx);
            fmpz_mpoly_mul(lc_divs->coeffs + j, lc_divs->coeffs + j, t, ctx);
            fmpz_pow_ui(q, lcAfaceval + i, k);
            fmpz_mul(dtilde + j, dtilde + j, q);
        }
    }

    fmpz_set(delta, &Aufac->c);
    for (j = 0; j < Aufac->num; j++)
    {
        FLINT_ASSERT(Aufac->p[j].length > 0);
        fmpz_gcd(r, Aufac->p[j].coeffs + Aufac->p[j].length - 1, dtilde + j);
        FLINT_ASSERT(fmpz_divisible(Aufac->p[j].coeffs + Aufac->p[j].length - 1, r));
        fmpz_divexact(q, Aufac->p[j].coeffs + Aufac->p[j].length - 1, r);
        fmpz_mpoly_scalar_mul_fmpz(lc_divs->coeffs + j, lc_divs->coeffs + j, q, ctx);
    }

    success = 1;

cleanup:

    fmpz_clear(delta);
    fmpz_clear(q);
    fmpz_clear(r);
    fmpz_mpoly_clear(t, ctx);
    _fmpz_vec_clear(lcAfaceval, lcAfac->num);
    _fmpz_vec_clear(d, 1 + lcAfac->num);
    _fmpz_vec_clear(dtilde, Aufac->num);
    flint_free(salpha);

    return success;
}


void nmod_poly_set_fmpz_poly(nmod_poly_t a, const fmpz_poly_t b)
{
    slong i;
    nmod_poly_fit_length(a, b->length);
    for (i = 0; i < b->length; i++)
        a->coeffs[i] = fmpz_fdiv_ui(b->coeffs + i, a->mod.n);
    a->length = b->length;
    _nmod_poly_normalise(a);
}

void _nmod_mpoly_set_fmpz_mpoly(
    nmod_mpoly_t Ap,
    const nmod_mpoly_ctx_t ctxp,
    const fmpz_mpoly_t A,
    const fmpz_mpoly_ctx_t ctx)
{
    slong N = mpoly_words_per_exp(A->bits, ctx->minfo);
    slong Ap_len, i;
    FLINT_ASSERT(ctx->minfo->nvars == ctx->minfo->nvars);
    FLINT_ASSERT(ctx->minfo->ord == ctx->minfo->ord);
    nmod_mpoly_fit_bits(Ap, A->bits, ctxp);
    Ap->bits = A->bits;
    nmod_mpoly_fit_length(Ap, A->length, ctxp);
    Ap_len = 0;
    for (i = 0; i < A->length; i++)
    {
        Ap->coeffs[Ap_len] = fmpz_fdiv_ui(A->coeffs + i, ctxp->ffinfo->mod.n);
        if (Ap->coeffs[Ap_len] == 0)
            continue;
        mpoly_monomial_set(Ap->exps + N*Ap_len, A->exps + N*i, N);
        Ap_len++;
    }
    Ap->length = Ap_len;
}

void fmpz_set_nmods(fmpz_t a, ulong n, ulong p)
{
    if (p - n < n)
        fmpz_neg_ui(a, p - n);
    else
        fmpz_set_ui(a, n);
}

/*
    Operation on each coeff:
        A = (A < 0 ? A + pk : A) + pk*mods(Ap, p)
*/
void _fmpz_mpoly_modpk_update(
    const fmpz_t pk,
    fmpz_mpoly_t A,
    const fmpz_mpoly_ctx_t ctx,
    const nmod_mpoly_t Ap,
    const nmod_mpoly_ctx_t ctxp)
{
    slong i;
    fmpz_mpoly_t T;
    slong N = mpoly_words_per_exp(Ap->bits, ctx->minfo);
/*
    for (i = 0; i < A->length; i++)
    {
        if (fmpz_sgn(A->coeffs + i) < 0)
            fmpz_add(A->coeffs + i, A->coeffs + i, pk);
    }

printf("positive A: "); fmpz_mpoly_print_pretty(A, NULL, ctx); printf("\n");
*/

    fmpz_mpoly_init3(T, Ap->length, Ap->bits, ctx);
    T->length = Ap->length;
    for (i = 0; i < Ap->length; i++)
    {
        fmpz_set_nmods(T->coeffs + i, Ap->coeffs[i], ctxp->ffinfo->mod.n);
        fmpz_mul(T->coeffs + i, T->coeffs + i, pk);
        mpoly_monomial_set(T->exps + N*i, Ap->exps + N*i, N);
    }

    fmpz_mpoly_add(A, A, T, ctx);
    fmpz_mpoly_clear(T, ctx);
}

void _fmpz_mpoly_modpk_taylor_coeff(
    const fmpz_t pk,
    nmod_mpoly_t T,
    const nmod_mpoly_ctx_t ctxp,
    const fmpz_mpoly_t E,
    const fmpz_mpoly_ctx_t ctx)
{
    slong i, Tlen;
    slong N = mpoly_words_per_exp(E->bits, ctx->minfo);
    fmpz_t t;

    fmpz_init(t);

    nmod_mpoly_fit_bits(T, E->bits, ctxp);
    T->bits = E->bits;
    nmod_mpoly_fit_length(T, E->length, ctxp);
    Tlen = 0;
    for (i = 0; i < E->length; i++)
    {
        FLINT_ASSERT(fmpz_divisible(E->coeffs + i, pk));
        fmpz_divexact(t, E->coeffs + i, pk);
        T->coeffs[Tlen] = fmpz_fdiv_ui(t, ctxp->ffinfo->mod.n);
        if (T->coeffs[Tlen] == 0)
            continue;
        mpoly_monomial_set(T->exps + N*Tlen, E->exps + N*i, N);
        Tlen++;
    }
    T->length = Tlen;

    fmpz_clear(t);
}

void _fmpz_mpoly_set_nmod_mpoly_smod(
    fmpz_mpoly_t A,
    const fmpz_mpoly_ctx_t ctx,
    const nmod_mpoly_t Ap,
    const nmod_mpoly_ctx_t ctxp)
{
    slong N = mpoly_words_per_exp(Ap->bits, ctxp->minfo);
/*
flint_printf("_fmpz_mpoly_set_nmod_mpoly_smod called bits = %wu\n", Ap->bits);
*/
    FLINT_ASSERT(ctx->minfo->ord == ctxp->minfo->ord);
    FLINT_ASSERT(ctx->minfo->nvars == ctxp->minfo->nvars);

    A->length = 0;
    fmpz_mpoly_fit_bits(A, Ap->bits, ctx);
    A->bits = Ap->bits;

    fmpz_mpoly_fit_length(A, Ap->length, ctx);
    A->length = Ap->length;

    mpoly_copy_monomials(A->exps, Ap->exps, Ap->length, N);
    _fmpz_vec_set_nmod_vec(A->coeffs, Ap->coeffs, Ap->length,
                                                      ctxp->ffinfo->mod);
/*
flint_printf("_fmpz_mpoly_set_nmod_mpoly_smod returning\n");
*/
}

/*
    The B[i] have the correct lc_x1 over ZZ.
    A = prod_i Bp[i] mod p
*/
int fmpz_mfactor_lift_prime_power(
    slong r,
    fmpz_mpoly_struct * B,
    const nmod_mpoly_struct * Bp,
    const fmpz_mpoly_t A,
    const slong * degs,
    const mp_limb_t * alphap,
    const fmpz_mpoly_ctx_t ctx,
    const nmod_mpoly_ctx_t ctxp,
    slong L)
{
    int success;
    slong i, k;
    nmod_disolve_t I;
    slong n = ctx->minfo->nvars - 1;
    fmpz_t pk;
    nmod_mpoly_struct * deltas;
    fmpz_mpoly_t e, t1, t2;
    nmod_mpoly_t tk;
    const char * vars [] = {"x", "y", "z", "w", "t", "u" ,"v"};

flint_printf("fmpz_mfactor_lift_prime_power\n");

    for (i = 0; i < r; i++)
    {
flint_printf("initial B[%wd]: ", i); fmpz_mpoly_print_pretty(B + i, vars, ctx); printf("\n");
    }


    FLINT_ASSERT(r > 1);

    nmod_disolve_init(I, r, n, Bp, alphap, ctxp);
    deltas = I->deltas + n*I->l;

    fmpz_init(pk);
    fmpz_one(pk);

    fmpz_mpoly_init(e, ctx);
    fmpz_mpoly_init(t1, ctx);
    fmpz_mpoly_init(t2, ctx);

    nmod_mpoly_init(tk, ctxp);

    k = 1;
    while (1)
    {
        fmpz_mul_ui(pk, pk, ctxp->ffinfo->mod.n);
        fmpz_mpoly_mul(t1, B + 0, B + 1, ctx);
        for (i = 2; i < r; i++)
        {
            fmpz_mpoly_mul(t2, t1, B + i, ctx);
            fmpz_mpoly_swap(t1, t2, ctx);
        }
        fmpz_mpoly_sub(e, A, t1, ctx);

flint_printf("e: "); fmpz_mpoly_print_pretty(e, vars, ctx); flint_printf("\n");

        if (fmpz_mpoly_is_zero(e, ctx))
        {
            success = 1;
            goto cleanup;
        }

        if (k > L)
        {
            success = 0;
            goto cleanup;
        }

        _fmpz_mpoly_modpk_taylor_coeff(pk, tk, ctxp, e, ctx);

flint_printf("calling nmod disolve\n");
        success = nmod_mfactor_disolve(A->bits, n, r, tk, alphap, degs, I, ctxp);
        FLINT_ASSERT(success);

        for (i = 0; i < r; i++)
        {
flint_printf("delta[%wd]: ", i); nmod_mpoly_print_pretty(deltas + i, vars, ctxp); flint_printf("\n");
            _fmpz_mpoly_modpk_update(pk, B + i, ctx, deltas + i, ctxp);
flint_printf("updated B[%wd]: ", i); fmpz_mpoly_print_pretty(B + i, vars, ctx); flint_printf("\n");

        }
    }

    FLINT_ASSERT(0 && "unreachable");

cleanup:

    nmod_disolve_clear(I, ctxp);

    return success;
}





void nmod_mpoly_monomial_evals(
    n_poly_t E,
    const nmod_mpoly_t A,
    const mp_limb_t * alpha,
    const nmod_mpoly_ctx_t ctx)
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
            xpoweval = nmod_mul(xpoweval, xpoweval, ctx->ffinfo->mod);
        }
    }
    FLINT_ASSERT(LUTlen < N*FLINT_BITS);

    n_poly_fit_length(E, A->length);
    E->length = A->length;
    for (i = 0; i < A->length; i++)
    {
        xpoweval = 1;
        for (j = 0; j < LUTlen; j++)
        {
            if (((Aexp + N*i)[LUToffset[j]] & LUTmask[j]) != 0)
            {
                xpoweval = nmod_mul(xpoweval, LUTvalue[j], ctx->ffinfo->mod);
            }
        }
        E->coeffs[i] = xpoweval;
    }

    TMP_END;
}

void _nmod_mpoly_monomial_evals(
    mp_limb_t * E,
    const ulong * Aexps,
    flint_bitcnt_t Abits,
    slong Alen,
    const mp_limb_t * alpha,
    slong vstart,
    const nmod_mpoly_ctx_t ctx)
{
    slong i, j;
    slong offset, shift;
    slong N = mpoly_words_per_exp_sp(Abits, ctx->minfo);
    slong nvars = ctx->minfo->nvars;
    slong * LUToffset;
    ulong * LUTmask;
    mp_limb_t * LUTvalue;
    slong LUTlen;
    mp_limb_t xpoweval;
    ulong * inputexpmask;
    TMP_INIT;

    TMP_START;
/*
flint_printf("_nmod_mpoly_monomial_evals called Alen: %wd\n", Alen);
*/
    LUToffset = (slong *) TMP_ALLOC(N*FLINT_BITS*sizeof(slong));
    LUTmask   = (ulong *) TMP_ALLOC(N*FLINT_BITS*sizeof(ulong));
    LUTvalue  = (mp_limb_t *) TMP_ALLOC(N*FLINT_BITS*sizeof(mp_limb_t));

    inputexpmask = (ulong *) TMP_ALLOC(N*sizeof(ulong));
    mpoly_monomial_zero(inputexpmask, N);
    for (i = 0; i < Alen; i++)
    {
        for (j = 0; j < N; j++)
        {
            inputexpmask[j] |= (Aexps + N*i)[j];
        }
    }

    LUTlen = 0;
    for (j = nvars - 1; j >= vstart; j--)
    {
        mpoly_gen_offset_shift_sp(&offset, &shift, j, Abits, ctx->minfo);

        xpoweval = alpha[j]; /* xpoweval = alpha[i]^(2^i) */
        for (i = 0; i < Abits; i++)
        {
            LUToffset[LUTlen] = offset;
            LUTmask[LUTlen] = (UWORD(1) << (shift + i));
            LUTvalue[LUTlen] = xpoweval;
            if ((inputexpmask[offset] & LUTmask[LUTlen]) != 0)
                LUTlen++;
            xpoweval = nmod_mul(xpoweval, xpoweval, ctx->ffinfo->mod);
        }
    }
    FLINT_ASSERT(LUTlen < N*FLINT_BITS);

    for (i = 0; i < Alen; i++)
    {
        xpoweval = 1;
        for (j = 0; j < LUTlen; j++)
        {
            if (((Aexps + N*i)[LUToffset[j]] & LUTmask[j]) != 0)
            {
                xpoweval = nmod_mul(xpoweval, LUTvalue[j], ctx->ffinfo->mod);
            }
        }
        E[i] = xpoweval;
    }

    TMP_END;
}

void _nmod_mpoly_monomial_evals_indirect(
    mp_limb_t * E,
    const ulong * Aexps,
    flint_bitcnt_t Abits,
    ulong * Aind,
    slong Alen,
    const mp_limb_t * alpha,
    slong vstart,
    slong vstop,
    const nmod_mpoly_ctx_t ctx)
{
    slong i, j;
    slong offset, shift;
    slong N = mpoly_words_per_exp_sp(Abits, ctx->minfo);
    slong * LUToffset;
    ulong * LUTmask;
    mp_limb_t * LUTvalue;
    slong LUTlen;
    mp_limb_t xpoweval;
    ulong * inputexpmask;
    const ulong * thisAexp;
    TMP_INIT;

    FLINT_ASSERT(0 <= vstart);
    FLINT_ASSERT(vstart < vstop);
    FLINT_ASSERT(vstop <= ctx->minfo->nvars);

    TMP_START;

    LUToffset = (slong *) TMP_ALLOC(N*FLINT_BITS*sizeof(slong));
    LUTmask   = (ulong *) TMP_ALLOC(N*FLINT_BITS*sizeof(ulong));
    LUTvalue  = (mp_limb_t *) TMP_ALLOC(N*FLINT_BITS*sizeof(mp_limb_t));

    inputexpmask = (ulong *) TMP_ALLOC(N*sizeof(ulong));
    mpoly_monomial_zero(inputexpmask, N);
    for (i = 0; i < Alen; i++)
    {
        thisAexp = Aexps + N*Aind[i];
        for (j = 0; j < N; j++)
            inputexpmask[j] |= thisAexp[j];
    }

    LUTlen = 0;
    for (j = vstop - 1; j >= vstart; j--)
    {
        mpoly_gen_offset_shift_sp(&offset, &shift, j, Abits, ctx->minfo);

        xpoweval = alpha[j]; /* xpoweval = alpha[i]^(2^i) */
        for (i = 0; i < Abits; i++)
        {
            LUToffset[LUTlen] = offset;
            LUTmask[LUTlen] = (UWORD(1) << (shift + i));
            LUTvalue[LUTlen] = xpoweval;
            if ((inputexpmask[offset] & LUTmask[LUTlen]) != 0)
                LUTlen++;
            xpoweval = nmod_mul(xpoweval, xpoweval, ctx->ffinfo->mod);
        }
    }
    FLINT_ASSERT(LUTlen < N*FLINT_BITS);

    for (i = 0; i < Alen; i++)
    {
        thisAexp = Aexps + N*Aind[i];
        xpoweval = 1;
        for (j = 0; j < LUTlen; j++)
        {
            if ((thisAexp[LUToffset[j]] & LUTmask[j]) != 0)
            {
                xpoweval = nmod_mul(xpoweval, LUTvalue[j], ctx->ffinfo->mod);
            }
        }
        E[i] = xpoweval;
    }

    TMP_END;
}




/*
    evaluation helper:
        3*i + 0 monomial evaluated at current power (changing)
        3*i + 1 coefficient (constant)
        3*i + 2 monomial evaluated at first power (constant)
*/

void nmod_mpolyu_set_eval_helper(
    n_polyun_t EH,
    const nmod_mpolyu_t A,
    const mp_limb_t * alpha,
    const nmod_mpoly_ctx_t ctx)
{
    slong i, j, n;
    n_polyun_term_struct * EHterms;
    mp_limb_t * p, * q;

    n_polyun_fit_length(EH, A->length);
    EH->length = A->length;
    EHterms = EH->terms;

    for (i = 0; i < A->length; i++)
    {
        EHterms[i].exp = A->exps[i];
        n = A->coeffs[i].length;
        n_poly_fit_length(EHterms[i].coeff, 3*n);
        nmod_mpoly_monomial_evals(EHterms[i].coeff, A->coeffs + i, alpha, ctx);
        FLINT_ASSERT(n == EHterms[i].coeff->length);
        p = EHterms[i].coeff->coeffs;
        q = A->coeffs[i].coeffs;
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


void nmod_mpoly_delete_duplicate_terms(
    nmod_mpoly_t A,
    const nmod_mpoly_ctx_t ctx)
{
    slong i, j;
    slong N = mpoly_words_per_exp(A->bits, ctx->minfo);
    j = -1;
    for (i = 0; i < A->length; i++)
    {
        if (j >= 0 && mpoly_monomial_equal(A->exps + N*j, A->exps + N*i, N))
        {
            FLINT_ASSERT(A->coeffs[j] == A->coeffs[i]);
            continue;
        }
        j++;
        A->coeffs[j] = A->coeffs[i];
        mpoly_monomial_set(A->exps + N*j, A->exps + N*i, N);
    }
    j++;
    A->length = j;
}

/*
    for each term Y^y*X^x*Z^z * pol(x1,...) in B with j < deg
    set Y^0*X^x*Z^z in H as the monomials with the monomial evals as coeffs
        merge monomial sets comming from different y's (shouldn't happen)
*/
slong nmod_mpolyu_set_eval_helper_and_zip_form(
    n_polyun_t EH,
    nmod_mpolyu_t H,
    ulong deg,
    const nmod_mpolyu_t B,
    const mp_limb_t * alpha,
    const nmod_mpoly_ctx_t ctx)
{
    slong i, j, n;
    ulong x, y, z;
    n_polyun_term_struct * EHterms;
    mp_limb_t * p, * q;
    nmod_mpoly_struct * Hc;
    slong old_len, zip_length = 0;
    slong N = mpoly_words_per_exp(B->bits, ctx->minfo);
/*
flint_printf("nmod_mpolyu_set_eval_helper_and_zip_form called\n");
flint_printf("deg = %wu\n", deg);
flint_printf("B: ");
nmod_mpolyu3_print_pretty(B, "Y", "X", "Z", NULL, ctx);
flint_printf("\n");
*/
    n_polyun_fit_length(EH, B->length);
    EH->length = B->length;
    EHterms = EH->terms;

    H->length = 0;

    for (i = 0; i < B->length; i++)
    {
        EHterms[i].exp = B->exps[i];
        y = extract_exp(EHterms[i].exp, 2, 3);
        x = extract_exp(EHterms[i].exp, 1, 3);
        z = extract_exp(EHterms[i].exp, 0, 3);
        n = B->coeffs[i].length;
        n_poly_fit_length(EHterms[i].coeff, 3*n);
        nmod_mpoly_monomial_evals(EHterms[i].coeff, B->coeffs + i, alpha, ctx);
        FLINT_ASSERT(n == EHterms[i].coeff->length);
        p = EHterms[i].coeff->coeffs;
        q = B->coeffs[i].coeffs;

        if (x < deg)
        {
            FLINT_ASSERT(y == 0 && "strange but ok");
            Hc = _nmod_mpolyu_get_coeff(H, pack_exp3(0, x, z), ctx);
            nmod_mpoly_fit_length(Hc, n, ctx);
            old_len = Hc->length;
            flint_mpn_copyi(Hc->coeffs + old_len, p, n);
            mpoly_copy_monomials(Hc->exps + N*old_len, B->coeffs[i].exps, n, N);
            Hc->length += n;
            zip_length = FLINT_MAX(zip_length, Hc->length);
            if (old_len > 0)
            {
                FLINT_ASSERT(0 && "strange but ok");
                nmod_mpoly_sort_terms(Hc, ctx);
                nmod_mpoly_delete_duplicate_terms(Hc, ctx);
            }
        }

        for (j = n - 1; j >= 0; j--)
        {
            mp_limb_t t1 = p[j];
            mp_limb_t t2 = q[j];
            p[3*j + 0] = t1;
            p[3*j + 1] = t2;
            p[3*j + 2] = t1;
        }
    }

    return zip_length;
}



/* if the coefficient doesn't exist, a new one is created (and set to zero) */
n_polyun_term_struct * n_polyun_get_term(n_polyun_t A, ulong k)
{
    slong i, j;
    n_polyun_term_struct * xk;

    for (i = 0; i < A->length && A->terms[i].exp >= k; i++)
    {
        if (A->terms[i].exp == k) 
        {
            return A->terms + i;
        }
    }

    n_polyun_fit_length(A, A->length + 1);

    for (j = A->length; j > i; j--)
        n_polyun_term_swap(A->terms + j, A->terms + j - 1);
    
    A->length++;

    xk = A->terms + i;
    xk->exp = k;
    xk->coeff->length = 0;
    return xk;
}

/*
    for each term Y^y*X^x*Z^z * pol(x1,...) in B with j < deg
    set Y^0*X^x*Z^z in H as the monomials with the monomial evals as coeffs
        merge monomial sets comming from different y's (shouldn't happen)
*/
static slong nmod_mpoly_set_eval_helper_and_zip_form(
    ulong * deg_,       /* deg_X(B), output */
    n_polyun_t EH,
    nmod_mpolyu_t H,
    const nmod_mpoly_t B,
    const mp_limb_t * alpha,
    slong yvar,         /* Y = gen(yvar) (X = gen(0), Z = gen(1))*/
    const nmod_mpoly_ctx_t ctx)
{
    slong xvar = 0;
    slong zvar = 1;
    slong i, j, n;
    ulong y, x, z;
    slong yoff, xoff, zoff;
    slong yshift, xshift, zshift;
    n_polyun_term_struct * EHterms;
    mp_limb_t * p;
    nmod_mpoly_struct * Hc;
    slong old_len, zip_length = 0;
    flint_bitcnt_t bits = B->bits;
    slong Blen = B->length;
    const ulong * Bexps = B->exps;
    const mp_limb_t * Bcoeffs = B->coeffs;
    slong N = mpoly_words_per_exp(bits, ctx->minfo);
    ulong mask = (-UWORD(1)) >> (FLINT_BITS - bits);
    ulong * ind;
    n_polyun_t T;
    n_polyun_term_struct * Tt;
    ulong deg;
/*
flint_printf("nmod_mpoly_set_eval_helper_and_zip_form called\n");
flint_printf("deg = %wu\n", deg);
flint_printf("B: ");
nmod_mpoly_print_pretty(B, NULL, ctx);
flint_printf("\n");
*/
    FLINT_ASSERT(bits <= FLINT_BITS);
    FLINT_ASSERT(bits == B->bits);
    FLINT_ASSERT(bits == H->bits);
    FLINT_ASSERT(Blen > 0);

    n_polyun_init(T);

    mpoly_gen_offset_shift_sp(&yoff, &yshift, yvar, bits, ctx->minfo);
    mpoly_gen_offset_shift_sp(&xoff, &xshift, xvar, bits, ctx->minfo);
    mpoly_gen_offset_shift_sp(&zoff, &zshift, zvar, bits, ctx->minfo);

    deg = (Bexps[N*0 + xoff] >> xshift) & mask;
    FLINT_ASSERT(deg == nmod_mpoly_degree_si(B, 0, ctx));

    /* TODO use a map here instead of this shit */
    for (i = 0; i < Blen; i++)
    {
        y = (Bexps[N*i + yoff] >> yshift) & mask;
        x = (Bexps[N*i + xoff] >> xshift) & mask;
        z = (Bexps[N*i + zoff] >> zshift) & mask;
        Tt = n_polyun_get_term(T, pack_exp3(y, x, z));
        FLINT_ASSERT(Tt->exp == pack_exp3(y, x, z));
        n_poly_fit_length(Tt->coeff, Tt->coeff->length + 1);
        Tt->coeff->coeffs[Tt->coeff->length] = i;
        Tt->coeff->length++;
    }

/*flint_printf("T:"); n_polyu3n_print_pretty(T, "Y", "X", "Z", "_"); flint_printf("\n");*/

    n_polyun_fit_length(EH, T->length);
    EH->length = T->length;
    EHterms = EH->terms;

    H->length = 0;

    for (i = 0; i < T->length; i++)
    {
        EHterms[i].exp = T->terms[i].exp;
        y = extract_exp(EHterms[i].exp, 2, 3);
        x = extract_exp(EHterms[i].exp, 1, 3);
        z = extract_exp(EHterms[i].exp, 0, 3);
        n = T->terms[i].coeff->length;
        n_poly_fit_length(EHterms[i].coeff, 3*n);
        EHterms[i].coeff->length = n;
        p = EHterms[i].coeff->coeffs;
        ind = T->terms[i].coeff->coeffs;
        _nmod_mpoly_monomial_evals_indirect(p, Bexps, bits, ind, n, alpha,
                                                                2, yvar, ctx);
        if (x < deg)
        {
            FLINT_ASSERT(y == 0 && "strange but ok");
            Hc = _nmod_mpolyu_get_coeff(H, pack_exp3(0, x, z), ctx);
            nmod_mpoly_fit_length(Hc, n, ctx);
            old_len = Hc->length;
            flint_mpn_copyi(Hc->coeffs + old_len, p, n);
            for (j = 0; j < n; j++)
            {
                mpoly_monomial_set(Hc->exps + N*(old_len + j),
                                   Bexps + N*ind[j], N);
            }
            Hc->length += n;
            zip_length = FLINT_MAX(zip_length, Hc->length);
            if (old_len > 0)
            {
                FLINT_ASSERT(0 && "strange but ok");
                nmod_mpoly_sort_terms(Hc, ctx);
                nmod_mpoly_delete_duplicate_terms(Hc, ctx);
            }
        }

        for (j = n - 1; j >= 0; j--)
        {
            mp_limb_t t1 = p[j];
            mp_limb_t t2 = Bcoeffs[ind[j]];
            p[3*j + 0] = t1;
            p[3*j + 1] = t2;
            p[3*j + 2] = t1;
        }
    }

    n_polyun_clear(T);

    *deg_ = deg;
/*
flint_printf("nmod_mpoly_set_eval_helper_and_zip_form returning\n");
*/
    return zip_length;
}



static mp_limb_t n_poly_mod_eval_step(n_poly_t A, nmod_t mod)
{
    slong i, Alen = A->length;
    mp_limb_t * Acoeffs = A->coeffs;
    ulong t0, t1, t2, p0, p1;

    FLINT_ASSERT(3*Alen <= A->alloc);

    t2 = t1 = t0 = 0;
    for (i = 0; i < Alen; i++)
    {
        umul_ppmm(p1, p0, Acoeffs[3*i + 0], Acoeffs[3*i + 1]);
        add_sssaaaaaa(t2, t1, t0, t2, t1, t0, 0, p1, p0);
        Acoeffs[3*i + 0] = nmod_mul(Acoeffs[3*i + 0], Acoeffs[3*i + 2], mod);
    }
    NMOD_RED3(t0, t2, t1, t0, mod);
    return t0;
}

void n_polyu_mod_eval_step(n_polyu_t E, n_polyun_t A, nmod_t mod)
{
    slong Ai, Ei;
    n_polyu_term_struct * Eterms;
    n_polyun_term_struct * Aterms;

    n_polyu_fit_length(E, A->length);

    Eterms = E->terms;
    Aterms = A->terms;
    Ei = 0;
    for (Ai = 0; Ai < A->length; Ai++)
    {
        FLINT_ASSERT(Ei < E->alloc);
        Eterms[Ei].exp = Aterms[Ai].exp;
        Eterms[Ei].coeff = n_poly_mod_eval_step(Aterms[Ai].coeff, mod);
        Ei += (Eterms[Ei].coeff != 0);
    }
    E->length = Ei;
}

void n_bpoly_mod_eval_step(n_bpoly_t E, n_polyun_t A, nmod_t mod)
{
    slong Ai;
    mp_limb_t c;
    ulong e0, e1;
    n_polyun_term_struct * Aterms = A->terms;

    n_bpoly_zero(E);
    for (Ai = 0; Ai < A->length; Ai++)
    {
        c = n_poly_mod_eval_step(Aterms[Ai].coeff, mod);
        e0 = extract_exp(Aterms[Ai].exp, 1, 2);
        e1 = extract_exp(Aterms[Ai].exp, 0, 2);
        if (c == 0)
            continue;
        n_bpoly_set_coeff_nonzero(E, e0, e1, c);
    }
}

void n_poly_eval_reset(n_poly_t A)
{
    slong i, Alen = A->length;
    mp_limb_t * Acoeffs = A->coeffs;

    FLINT_ASSERT(3*Alen <= A->alloc);

    for (i = 0; i < Alen; i++)
        Acoeffs[3*i + 0] = Acoeffs[3*i + 2];
}

void n_polyun_eval_reset(n_polyun_t A)
{
    slong Ai;
    for (Ai = 0; Ai < A->length; Ai++)
        n_poly_eval_reset(A->terms[Ai].coeff);
}

int n_polyu2_add_zip_must_match(
    n_polyun_t Z,
    const n_bpoly_t A,
    slong cur_length)
{
    slong i, Ai, ai;
    n_polyun_term_struct * Zt = Z->terms;
    const n_poly_struct * Acoeffs = A->coeffs;

    Ai = A->length - 1;
    ai = (Ai < 0) ? 0 : n_poly_degree(A->coeffs + Ai);

    for (i = 0; i < Z->length; i++)
    {
        if (Ai >= 0 && Zt[i].exp == pack_exp2(Ai, ai))
        {
            /* Z present, A present */
            Zt[i].coeff->coeffs[cur_length] = Acoeffs[Ai].coeffs[ai];
            Zt[i].coeff->length = cur_length + 1;
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
        else if (Ai < 0 || Zt[i].exp > pack_exp2(Ai, ai))
        {
            /* Z present, A missing */
            Zt[i].coeff->coeffs[cur_length] = 0;
            Zt[i].coeff->length = cur_length + 1;
        }
        else
        {
            /* Z missing, A present */
            return 0;
        }
    }

    return 1;
}

void n_polyu3_add_zip_limit1(
    n_polyun_t Z,
    const n_polyun_t A,
    const ulong deg1,
    slong cur_length,
    slong fit_length)
{
    const n_polyun_term_struct * At = A->terms;
    const n_polyun_term_struct * Ait;
    n_polyun_term_struct * Zit;
    slong Ai, ai, Zi, j;

    Ai = -1;
    ai = -1;
    do {
        Ai++;
        Ait = At + Ai;
    } while (Ai < A->length && extract_exp(Ait->exp, 1, 3) >= deg1);
    if (Ai < A->length)
        ai = n_poly_degree(Ait->coeff);

    Zi = 0;

    while (Ai < A->length && Zi < Z->length)
    {
        Zit = Z->terms + Zi;
        Ait = At + Ai;
        if (Ait->exp + ai > Zit->exp)
        {
            /* missing from Z */
            n_polyun_fit_length(Z, Z->length + 1);
            for (j = Z->length; j > Zi; j--)
                n_polyun_term_swap(Z->terms + j, Z->terms + j - 1);
            Z->length++;
            Zit = Z->terms + Zi;
            Zit->exp = Ait->exp + ai;
            n_poly_fit_length(Zit->coeff, fit_length);
            Zit->coeff->length = cur_length;
            mpn_zero(Zit->coeff->coeffs, cur_length);
            goto in_both;            
        }
        else if (Ait->exp + ai < Zit->exp)
        {
            /* missing from A */
            FLINT_ASSERT(cur_length + 1 <= Zit->coeff->alloc);
            Zit->coeff->coeffs[cur_length] = 0;
            Zit->coeff->length = cur_length + 1;
            Zi++;
        }
        else
        {
in_both:
            FLINT_ASSERT(cur_length == Zit->coeff->length);
            FLINT_ASSERT(cur_length + 1 <= Zit->coeff->alloc);
            Zit->coeff->coeffs[cur_length] = Ait->coeff->coeffs[ai];
            Zit->coeff->length = cur_length + 1;
            Zi++;
            do {
                ai--;
            } while (ai >= 0 && Ait->coeff->coeffs[ai] == 0);
            if (ai < 0)
            {
                do {
                    Ai++;
                    Ait = At + Ai;
                } while (Ai < A->length && extract_exp(Ait->exp, 1, 3) >= deg1);
                if (Ai < A->length)
                    ai = n_poly_degree(Ait->coeff);
            }
        }
    }

    /* everything in A must be put on the end of Z */
    while (Ai < A->length)
    {
        Zi = Z->length;
        n_polyun_fit_length(Z, Zi + A->length - Ai);
        Zit = Z->terms + Zi;
        Zit->exp = Ait->exp + ai;
        n_poly_fit_length(Zit->coeff, fit_length);
        Zit->coeff->length = cur_length;
        mpn_zero(Zit->coeff->coeffs, cur_length);
        Z->length = ++Zi;
        FLINT_ASSERT(cur_length == Zit->coeff->length);
        FLINT_ASSERT(cur_length + 1 <= Zit->coeff->alloc);
        Zit->coeff->coeffs[cur_length] = Ait->coeff->coeffs[ai];
        Zit->coeff->length = cur_length + 1;
        do {
            ai--;
        } while (ai >= 0 && Ait->coeff->coeffs[ai] == 0);
        if (ai < 0)
        {
            do {
                Ai++;
                Ait = At + Ai;
            } while (Ai < A->length && extract_exp(Ait->exp, 1, 3) >= deg1);
            if (Ai < A->length)
                ai = n_poly_degree(Ait->coeff);
        }
    }

    /* everything in Z must have a zero appended */
    while (Zi < Z->length)
    {
        Zit = Z->terms + Zi;
        FLINT_ASSERT(cur_length == Zit->coeff->length);
        FLINT_ASSERT(cur_length + 1 <= Zit->coeff->alloc);
        Zit->coeff->coeffs[cur_length] = 0;
        Zit->coeff->length = cur_length + 1;
        Zi++;
    }

    for (Zi = 0; Zi < Z->length; Zi++)
    {
        FLINT_ASSERT(Z->terms[Zi].coeff->length == cur_length + 1);
    }
}

slong nmod_mpolyu_find_term(const nmod_mpolyu_t A, ulong e)
{
    slong i;
    for (i = 0; i < A->length; i++)
        if (A->exps[i] == e)
            return i;
    return -1;
}

/*
    for each Y^y*X^x*Z^z in B with x = deg,
        keep the Y^y*X^x*Z^z*poly(x1,...) in B
    for each Y^y*X^x*Z^z in Z,
        assert that x < deg
        if there is no Y^0*X^x*Z^y in H, fail
        find coefficients of poly using this entry in H
        output Y^y*X^x*Z^z*poly(x1,...) to A
    sort A

    return
        -1: singular vandermonde matrix encountered
        0:  inconsistent system encountered
        1:  success
*/
static int nmod_mpoly_from_zip(
    nmod_mpoly_t B,
    const n_polyun_t Z,
    nmod_mpolyu_t H,
    ulong deg,
    slong yvar,     /* Y = gen(yvar) */
    const nmod_mpoly_ctx_t ctx)
{
    int success;
    slong Hi, Zi, Bi, i, j;
    slong xvar = 0;
    slong zvar = 1;
    ulong x, y, z;
    flint_bitcnt_t bits = B->bits;
    mp_limb_t * Bcoeffs;
    ulong * Bexps;
    slong N = mpoly_words_per_exp_sp(bits, ctx->minfo);
    ulong mask = (-UWORD(1)) >> (FLINT_BITS - bits);
    slong xoff, xshift, yoff, yshift, zoff, zshift;
    n_polyun_term_struct * Zt = Z->terms;
    nmod_mpoly_struct * Hc;
    slong Hlen = H->length;

    n_polyun_t M;  /* temp */
    n_polyun_init(M);
/*
flint_printf("-----------------");
flint_printf("nmod_mpoly_from_zip called vars %wd, %wd, %wd\n", yvar, xvar, zvar);
flint_printf("Z: "); n_polyu3n_print_pretty(Z, "Y", "X", "Z", "_"); printf("\n");
flint_printf("H: "); nmod_mpolyu3_print_pretty(H, "Y", "X", "Z", NULL, ctx); printf("\n");
flint_printf("deg: %wd\n", deg);
*/
    FLINT_ASSERT(bits == H->bits);

    n_polyun_fit_length(M, Hlen + 1);
    for (i = 0; i <= Hlen; i++)
        M->terms[i].coeff->length = 0;

    mpoly_gen_offset_shift_sp(&yoff, &yshift, yvar, bits, ctx->minfo);
    mpoly_gen_offset_shift_sp(&xoff, &xshift, xvar, bits, ctx->minfo);
    mpoly_gen_offset_shift_sp(&zoff, &zshift, zvar, bits, ctx->minfo);

    /* x is most significant in ctx, so keeping the lc_x in B is easy */
    FLINT_ASSERT(xvar == 0);
    
    for (Bi = 0; Bi < B->length; Bi++)
    {
        x = (((B->exps + N*Bi)[xoff] >> xshift) & mask);
        FLINT_ASSERT(x <= deg);
        if (x != deg)
            break;
    }

    for (Zi = 0; Zi < Z->length; Zi++)
    {
        y = extract_exp(Zt[Zi].exp, 2, 3);
        x = extract_exp(Zt[Zi].exp, 1, 3);
        z = extract_exp(Zt[Zi].exp, 0, 3);
        FLINT_ASSERT(x < deg);
        Hi = nmod_mpolyu_find_term(H, pack_exp3(0, x, z));
        if (Hi < 0)
            return 0;

        FLINT_ASSERT(Hi < Hlen);
        FLINT_ASSERT(H->exps[Hi] == pack_exp3(0, x, z));

        Hc = H->coeffs + Hi;
        FLINT_ASSERT(bits == Hc->bits);
        FLINT_ASSERT(Hc->length > 0);
        nmod_mpoly_fit_length(B, Bi + Hc->length, ctx);
        Bcoeffs = B->coeffs;

        if (M->terms[Hi].coeff->length < 1)
        {
            n_poly_mod_product_roots_nmod_vec(M->terms[Hi].coeff,
                                     Hc->coeffs, Hc->length, ctx->ffinfo->mod);
        }

        n_poly_fit_length(M->terms[Hlen].coeff, Hc->length);

        success = nmod_zip_find_coeffs_new2(Bcoeffs + Bi, Hc->coeffs,
                    Hc->length, Zt[Zi].coeff->coeffs, Zt[Zi].coeff->length,
                    M->terms[Hi].coeff->coeffs, M->terms[Hlen].coeff->coeffs,
                                                             ctx->ffinfo->mod);
        if (success <= 0)
            return success;

        Bexps = B->exps;
        for (j = Bi, i = 0; i < Hc->length; j++, i++)
        {
            if (Bcoeffs[j] == 0)
                continue;
            Bcoeffs[Bi] = Bcoeffs[j];
            FLINT_ASSERT(Bi < B->alloc);
            mpoly_monomial_set(Bexps + N*Bi, Hc->exps + N*i, N);
            (Bexps + N*Bi)[yoff] += y << yshift;
            Bi++;
        }
    }
    B->length = Bi;
    nmod_mpoly_sort_terms(B, ctx);
    FLINT_ASSERT(nmod_mpoly_is_canonical(B, ctx));
/*
flint_printf("nmod_mpoly_from_zip returning good\n");
flint_printf("B: "); nmod_mpoly_print_pretty(B, NULL, ctx); flint_printf("\n");
*/
    n_polyun_clear(M);

    return 1;
}


/*
    ctx is for Fp[x_0, ..., x_n] where n = ctx->minfo->nvars - 1

    We are lifting a factorization mod
        (x_m - alpha[m-1], ..., x_n - alpha[n-1])
    to a factorization mod
        (x_(m+1) - alpha[m], ..., x-n - alpha[n-1]),
    We are lifting variable x_m, and all polynomials live in Fp[x_0, ..., x_m].

    It is required that m >= 3. Rewrite polynomials as nmod_mpolyu's to
    live in Fp[x_2, ..., x_(m-1)][x_m, x_0, x_1] still using ctx for
    the Fp[x_2, ..., x_(m-1)] part. Evaluate inputs at
        (x_2, ..., x_(m-1)) = (beta_2^i, ..., beta_(m-1)^i) for i = 1, 2, ..., s
    The evaluations are trivariate n_polyu's in Fp[x_m, x_0, x_1], aka n_polyu3.
    Lift the evaluations in Fp[x_m, x_0, x_1] mod (x_m - alpha[m-1]) using
    n_polyu3_mod_factor_lift. The lifted factors come out as n_polyun's in
    Fp[x_1][x_m, x_0], i.e. x_1 is in dense storage. Interpolate the factors in
    Fp[x_2, ..., x_(m-1)][x_m, x_0, x_1] using Zippel's assumption.

    return:
       -1: function failed with no conclusion
        0: lifting is impossible without changing the lc_(x_0)
        1: B was successfully updated with the lifted factors
*/
int nmod_mpoly_factor_lift_tuncer(
    slong r,
    nmod_mpoly_struct * B,
    flint_rand_t state,
    const mp_limb_t * alpha,
    const nmod_mpoly_t A,
    const slong * degs,
    slong m,
    const nmod_mpoly_ctx_t ctx)
{
    int success;
    slong i;
    slong zip_fails_remaining;
    slong req_zip_images, cur_zip_image;
    nmod_mpolyu_struct Au[1], * H;
    n_polyun_struct Aeh[1], * Beh;
    n_polyu_struct Aeval[1], * Beval;
    n_polyun_struct * BBeval, * Z;
    mp_limb_t * beta;
    flint_bitcnt_t bits = A->bits;
    nmod_mpoly_t T1, T2;
    ulong * Bdegs;
    const slong degs0 = degs[0];
/*
    const char * vars[] = {"x", "y", "z", "w", "u", "v"};
flint_printf("nmod_mpoly_factor_lift_tuncer called (m = %wd, bits = %wu)\n", m, bits); fflush(stdout);
for (i = 0; i < n; i++)
{
flint_printf("alpha[%wd]: %wu\n", i, alpha[i]);
}
flint_printf("A: "); nmod_mpoly_print_pretty(A, vars, ctx); flint_printf("\n");
for (i = 0; i < r; i++)
{
flint_printf("B[%wd]: ", i); nmod_mpoly_print_pretty(B + i, vars, ctx); flint_printf("\n");
}
*/
    FLINT_ASSERT(m > 2);
    FLINT_ASSERT(r > 1);
    FLINT_ASSERT(bits <= FLINT_BITS);

    beta = (mp_limb_t *) flint_malloc(ctx->minfo->nvars*sizeof(mp_limb_t));
    Bdegs = (ulong *) flint_malloc(r*sizeof(ulong));
    H = (nmod_mpolyu_struct *) flint_malloc(r*sizeof(nmod_mpolyu_struct));
    Beh = (n_polyun_struct *) flint_malloc(r*sizeof(n_polyun_struct));
    Beval = (n_polyu_struct *) flint_malloc(r*sizeof(n_polyu_struct));
    BBeval = (n_polyun_struct *) flint_malloc(r*sizeof(n_polyun_struct));
    Z = (n_polyun_struct *) flint_malloc(r*sizeof(n_polyun_struct));

    nmod_mpolyu_init(Au, bits, ctx);
    n_polyun_init(Aeh);
    n_polyu_init(Aeval);
    for (i = 0; i < r; i++)
    {
        nmod_mpolyu_init(H + i, bits, ctx);
        n_polyun_init(Beh + i);
        n_polyu_init(Beval + i);
        n_polyun_init(BBeval + i);
        n_polyun_init(Z + i);
    }

    /* init done */

    for (i = 0; i < r; i++)
    {
        success = nmod_mpoly_repack_bits_inplace(B + i, bits, ctx);
        if (!success)
            goto cleanup;
    }

    zip_fails_remaining = 3;

choose_betas:

    /* only beta[2], beta[3], ..., beta[m - 1] will be used */
    FLINT_ASSERT(ctx->ffinfo->mod.n > 3);
    for (i = 0; i < ctx->minfo->nvars; i++)
        beta[i] = n_urandint(state, ctx->ffinfo->mod.n - 3) + 2;

    nmod_mpoly_get_mpolyu3(Au, A, m, 0, 1, ctx);
/*flint_printf("Au: ", i); nmod_mpolyu3_print_pretty(Au, vars[m], vars[0], vars[1], vars, ctx); flint_printf("\n");*/
    nmod_mpolyu_set_eval_helper(Aeh, Au, beta, ctx);
/*flint_printf("Aeh: "); n_polyu3n_print_pretty(Aeh, vars[m], vars[0], vars[1], "_"); flint_printf("\n");*/

    req_zip_images = 1;
    for (i = 0; i < r; i++)
    {
        slong this_zip_images;
        this_zip_images = nmod_mpoly_set_eval_helper_and_zip_form(Bdegs + i,
                                          Beh + i, H + i, B + i, beta, m, ctx);
        req_zip_images = FLINT_MAX(req_zip_images, this_zip_images);
        FLINT_ASSERT(Bdegs[i] > 0);

/*flint_printf("Beh[%wd]: ", i); n_polyu3n_print_pretty(Beh + i, vars[m], vars[0], vars[1], "_"); flint_printf("\n");
flint_printf("H[%wd]: ", i); nmod_mpolyu3_print_pretty(H + i, vars[m], vars[0], vars[1], vars, ctx); flint_printf("\n");
flint_printf("req_zip_images: %wd\n", req_zip_images);*/
    }

    cur_zip_image = 0;

next_zip_image:

    n_polyu_mod_eval_step(Aeval, Aeh, ctx->ffinfo->mod);
/*flint_printf("Aeval: "); n_polyu3_print_pretty(Aeval, vars[m], vars[0], vars[1]); flint_printf("\n");*/
    for (i = 0; i < r; i++)
    {
        n_polyu_mod_eval_step(Beval + i, Beh + i, ctx->ffinfo->mod);
/*flint_printf("Beval[%wd]: ", i); n_polyu3_print_pretty(Beval + i, vars[m], vars[0], vars[1]); flint_printf("\n");*/
    }

    success = n_polyu3_mod_factor_lift(r, BBeval, Aeval, Beval,
                                             alpha[m - 1], degs0, ctx->ffinfo);
    if (success < 1)
    {
        FLINT_ASSERT(0 && "spurious failure");
        if (--zip_fails_remaining >= 0)
            goto choose_betas;
        success = 0;
        goto cleanup;
    }

    for (i = 0; i < r; i++)
    {
/*flint_printf("BBeval[%wd]: ", i); n_polyu3n_print_pretty(BBeval + i, vars[m], vars[0], "?", vars[1]); flint_printf("\n");*/
        n_polyu3_add_zip_limit1(Z + i, BBeval + i, Bdegs[i],
                                                cur_zip_image, req_zip_images); 
/*flint_printf("Z[%wd]: ", i); n_polyu3n_print_pretty(Z + i, vars[m], vars[0], vars[1], "_"); flint_printf("\n");*/
    }

    cur_zip_image++;
    if (cur_zip_image < req_zip_images)
        goto next_zip_image;

    for (i = 0; i < r; i++)
    {
        success = nmod_mpoly_from_zip(B + i, Z + i, H + i, Bdegs[i], m, ctx);
        if (success < 1)
        {
            FLINT_ASSERT(0 && "spurious failure");
            success = 0;
            goto cleanup;
        }
    }

    nmod_mpoly_init3(T1, A->length, bits, ctx);
    nmod_mpoly_init3(T2, A->length, bits, ctx);
    nmod_mpoly_mul(T1, B + 0, B + 1, ctx);
    for (i = 2; i < r; i++)
    {
        nmod_mpoly_mul(T2, T1, B + i, ctx);
        nmod_mpoly_swap(T1, T2, ctx);
    }

    success = nmod_mpoly_equal(T1, A, ctx);
    nmod_mpoly_clear(T1, ctx);
    nmod_mpoly_clear(T2, ctx);

cleanup:

    nmod_mpolyu_clear(Au, ctx);
    n_polyun_clear(Aeh);
    n_polyu_clear(Aeval);
    for (i = 0; i < r; i++)
    {
        nmod_mpolyu_clear(H + i, ctx);
        n_polyun_clear(Beh + i);
        n_polyu_clear(Beval + i);
        n_polyun_clear(BBeval + i);
        n_polyun_clear(Z + i);
    }

    flint_free(beta);
    flint_free(Bdegs);
    flint_free(H);
    flint_free(Beh);
    flint_free(Beval);
    flint_free(BBeval);
    flint_free(Z);

    return success;
}


slong nmod_mpoly_set_eval_helper_and_zip_form2(
    slong * deg1_, /* degree of B wrt main var 1 */
    n_polyun_t EH,
    n_polyun_t H,
    n_polyun_t M,
    const nmod_mpoly_t B,
    const mp_limb_t * betas,
    const nmod_mpoly_ctx_t ctx)
{
    slong start, Bi, j, n;
    slong e0, e1, Hi, EHi;
    n_polyun_term_struct * EHterms, * Hterms, * Mterms;
    mp_limb_t * p;
    slong zip_length = 0;
    flint_bitcnt_t Bbits = B->bits;
    const mp_limb_t * Bcoeffs = B->coeffs;
    slong Blen = B->length;
    const ulong * Bexps = B->exps;
    ulong mask = (-UWORD(1)) >> (FLINT_BITS - Bbits);
    slong N = mpoly_words_per_exp_sp(Bbits, ctx->minfo);
    slong off0, off1, shift0, shift1;
    slong deg0, deg1 = -1;
/*
flint_printf("nmod_mpoly_set_eval_helper_and_zip_form2 called\n");
flint_printf("B: "); nmod_mpoly_print_pretty(B, NULL, ctx); flint_printf("\n");
*/
    FLINT_ASSERT(Blen > 0);

    mpoly_gen_offset_shift_sp(&off0, &shift0, 0, Bbits, ctx->minfo);
    mpoly_gen_offset_shift_sp(&off1, &shift1, 1, Bbits, ctx->minfo);

    Bi = 0;
    deg0 = (Bexps[N*Bi + off0] >> shift0) & mask;

    EHi = 0;
    Hi = 0;

    while (Bi < Blen)
    {
        start = Bi;
        e0 = (Bexps[N*Bi + off0] >> shift0) & mask;
        e1 = (Bexps[N*Bi + off1] >> shift1) & mask;
        deg1 = FLINT_MAX(deg1, e1);
        while (1)
        {
            Bi++;
            if (Bi >= Blen)
                break;
            if (((Bexps[N*Bi + off0] >> shift0) & mask) != e0)
                break;
            if (((Bexps[N*Bi + off1] >> shift1) & mask) != e1)
                break;
        }

        n = Bi - start;

        n_polyun_fit_length(EH, EHi + 1);
        EHterms = EH->terms;
        EHterms[EHi].exp = pack_exp2(e0, e1);
        n_poly_fit_length(EHterms[EHi].coeff, 3*n);
        EHterms[EHi].coeff->length = n;
        p = EHterms[EHi].coeff->coeffs;
        EHi++;

        _nmod_mpoly_monomial_evals(p, Bexps + N*start, Bbits, n, betas, 2, ctx);

        if (e0 < deg0)
        {
            n_polyun_fit_length(H, Hi + 1);
            n_polyun_fit_length(M, Hi + 1);
            Hterms = H->terms;
            Mterms = M->terms;
            Hterms[Hi].exp = pack_exp2(e0, e1);
            Mterms[Hi].exp = pack_exp2(e0, e1);
            n_poly_fit_length(Hterms[Hi].coeff, n);
            zip_length = FLINT_MAX(zip_length, n);
            Hterms[Hi].coeff->length = n;
            flint_mpn_copyi(Hterms[Hi].coeff->coeffs, p, n);
            n_poly_mod_product_roots_nmod_vec(Mterms[Hi].coeff, p, n, ctx->ffinfo->mod);
            Hi++;
        }

        for (j = n - 1; j >= 0; j--)
        {
            mp_limb_t t2 = Bcoeffs[start + j];
            mp_limb_t t1 = p[j];
            p[3*j + 0] = t1;
            p[3*j + 1] = t2;
            p[3*j + 2] = t1;
        }
    }

    EH->length = EHi;
    H->length = Hi;
    M->length = Hi;
/*
flint_printf("nmod_mpoly_set_eval_helper_and_zip_form2 returning deg1 = %wd\n", deg1);
*/
    *deg1_ = deg1;
    return zip_length;
}


int _fmpz_mpoly_modpk_update_zip(
    fmpz_t pk,
    fmpz_mpoly_t A,
    const fmpz_mpoly_ctx_t ctx,
    n_polyun_t Z,
    const n_polyun_t H,
    const n_polyun_t M,
    const nmod_mpoly_ctx_t ctxp)
{
    slong i, j, Ai, n;
    int success;
    slong off, shift;
    slong N = mpoly_words_per_exp_sp(A->bits, ctx->minfo);
    ulong start, mask = (-UWORD(1)) >> (FLINT_BITS - A->bits);
    n_poly_t c, t;
    mp_limb_t * ccoeffs;

    mpoly_gen_offset_shift_sp(&off, &shift, 0, A->bits, ctx->minfo);

    mask = mask << shift;

    Ai = 1;
    start = (A->exps + N*0)[off] & mask;

    while (Ai < A->length && start == ((A->exps + N*Ai)[off] & mask))
    {
        Ai++;
    }

    FLINT_ASSERT(Ai < A->length);

    if (Ai >= A->length)
        return 1;

    n_poly_init(c);
    n_poly_init(t);

    FLINT_ASSERT(Z->length == H->length);
    FLINT_ASSERT(Z->length == M->length);

    for (i = 0; i < Z->length; i++)
    {
        n = H->terms[i].coeff->length;
        FLINT_ASSERT(M->terms[i].coeff->length == n + 1);
        FLINT_ASSERT(Z->terms[i].coeff->length >= n);

        n_poly_fit_length(c, n);
        n_poly_fit_length(t, n);

        ccoeffs = c->coeffs;

        success = nmod_zip_find_coeffs_new2(c->coeffs,
                        H->terms[i].coeff->coeffs, n,
                        Z->terms[i].coeff->coeffs, Z->terms[i].coeff->length,
                        M->terms[i].coeff->coeffs, t->coeffs,
                                                            ctxp->ffinfo->mod);
        if (success <= 0)
        {
            n_poly_clear(t);
            n_poly_clear(c);
            return success;
        }

        FLINT_ASSERT(Ai + n <= A->length);

        for (j = 0; j < n; j++)
        {
            if (ctxp->ffinfo->mod.n - ccoeffs[j] < ccoeffs[j])
                fmpz_submul_ui(A->coeffs + Ai + j, pk, ctxp->ffinfo->mod.n - ccoeffs[j]);
            else
                fmpz_addmul_ui(A->coeffs + Ai + j, pk, ccoeffs[j]);
        }

        Ai += n;
    }

    FLINT_ASSERT(Ai == A->length);

    n_poly_clear(t);
    n_poly_clear(c);

    return 1;
}


/*
    n_polyun_t Beh has all x0^i*x1^j*poly(x2, ...) with
        coeffs triples suitable for sequential eval at x2,... = beta^i

    n_polyun_t H has x0^i*x1^j*poly(x2, ...) with i < deg_x0(B[i]) with
        the monomial evals at beta in the coeffs

    n_polyun_t M has x0^i*x1^j*poly(x2, ...) with i < deg_x0(B[i])
        the master poly of H[x0^i*x1^j] in the coeff

    n_polyun_t Z has x0^i*x1^j with i < deg_x0(B[i])
        ready to collect images
*/
int fmpz_mfactor_lift_prime_power_tuncer(
    slong r,
    fmpz_mpoly_struct * B,
    flint_rand_t state,
    const nmod_mpoly_struct * Bp,
    const fmpz_mpoly_t A,
    const mp_limb_t * alphap,
    const fmpz_mpoly_ctx_t ctx,
    const nmod_mpoly_ctx_t ctxp,
    slong L)
{
    slong req_zip_images, cur_zip_image;
    flint_bitcnt_t bits = A->bits;
    slong n = ctxp->minfo->nvars;
    int success;
    slong i, j, k;
    n_polyun_struct * H, * M, * Z, * Beh;
    n_bpoly_struct * Ceval, * Beval, Teval[1];
    slong * Cdegs1;
    nmod_mpoly_t Tp;
    nmod_mpolyu_t Tu;
    n_polyun_t Teh;
    fmpz_t pk;
    mp_limb_t * betas;
    fmpz_mpoly_t e, t1, t2;
/*
    const char * vars [] = {"x", "y", "z", "w", "t", "u" ,"v"};

flint_printf("fmpz_mfactor_lift_prime_power_tuncer called\n");
    for (i = 0; i < r; i++)
    {
flint_printf("initial B[%wd]: ", i); fmpz_mpoly_print_pretty(B + i, vars, ctx); flint_printf("\n");

    }
*/
    FLINT_ASSERT(r > 1);

    fmpz_init(pk);

    fmpz_mpoly_init(e, ctx);
    fmpz_mpoly_init(t1, ctx);
    fmpz_mpoly_init(t2, ctx);

    nmod_mpoly_init(Tp, ctxp);
    nmod_mpolyu_init(Tu, bits, ctxp);
    n_polyun_init(Teh);

    Cdegs1 = (slong *) flint_malloc(r*sizeof(slong));

    betas = (mp_limb_t *) flint_malloc(n*sizeof(mp_limb_t));

    FLINT_ASSERT(ctxp->ffinfo->mod.n > 3);
    for (i = 0; i < ctx->minfo->nvars; i++)
        betas[i] = n_urandint(state, ctxp->ffinfo->mod.n - 3) + 2;

    Beh = (n_polyun_struct *) flint_malloc(r*sizeof(n_polyun_struct));
    H = (n_polyun_struct *) flint_malloc(r*sizeof(n_polyun_struct));
    M = (n_polyun_struct *) flint_malloc(r*sizeof(n_polyun_struct));
    Z = (n_polyun_struct *) flint_malloc(r*sizeof(n_polyun_struct));
    Ceval = (n_bpoly_struct *) flint_malloc(r*sizeof(n_bpoly_struct));
    Beval = (n_bpoly_struct *) flint_malloc(r*sizeof(n_bpoly_struct));
    for (i = 0; i < r; i++)
    {
        n_polyun_init(Beh + i);
        n_polyun_init(H + i);
        n_polyun_init(M + i);
        n_polyun_init(Z + i);
        n_bpoly_init(Ceval + i);
        n_bpoly_init(Beval + i);
    }

    n_bpoly_init(Teval);

    /* init done */

    req_zip_images = 1;
    for (i = 0; i < r; i++)
    {
        slong this_images;
        this_images = nmod_mpoly_set_eval_helper_and_zip_form2(Cdegs1 + i,
                                   Beh + i, H + i, M + i, Bp + i, betas, ctxp);
        req_zip_images = FLINT_MAX(req_zip_images, this_images);
    }

    for (i = 0; i < r; i++)
    {
        n_polyun_fit_length(Z + i, H[i].length);
        Z[i].length = H[i].length;
        for (j = 0; j < H[i].length; j++)
        {
            Z[i].terms[j].exp = H[i].terms[j].exp;
            n_poly_fit_length(Z[i].terms[j].coeff, req_zip_images);
            Z[i].terms[j].coeff->length = 0;
        }
    }

    fmpz_one(pk);

    k = 1;

next_power:

    fmpz_mul_ui(pk, pk, ctxp->ffinfo->mod.n);

    fmpz_mpoly_mul(t1, B + 0, B + 1, ctx);
    for (i = 2; i < r; i++)
    {
        fmpz_mpoly_mul(t2, t1, B + i, ctx);
        fmpz_mpoly_swap(t1, t2, ctx);
    }
    fmpz_mpoly_sub(e, A, t1, ctx);

    if (fmpz_mpoly_is_zero(e, ctx))
    {
        success = 1;
        goto cleanup;
    }

    if (k > L)
    {
FLINT_ASSERT(0 && "spurious failure");
        success = 0;
        goto cleanup;
    }

    _fmpz_mpoly_modpk_taylor_coeff(pk, Tp, ctxp, e, ctx);

    nmod_mpoly_get_mpolyu2(Tu, Tp, 0, 1, ctxp);

    nmod_mpolyu_set_eval_helper(Teh, Tu, betas, ctxp);
    if (fmpz_cmp_ui(pk, ctxp->ffinfo->mod.n) > 0)
    {
        for (i = 0; i < r; i++)
            n_polyun_eval_reset(Beh + i);
    }

    cur_zip_image = 0;

next_zip_image:

    n_bpoly_mod_eval_step(Teval, Teh, ctxp->ffinfo->mod);
    for (i = 0; i < r; i++)
        n_bpoly_mod_eval_step(Beval + i, Beh + i, ctxp->ffinfo->mod);

    success = n_bpoly_mod_disolve(r, Ceval, Cdegs1, Teval, Beval,
                                                       ctxp->ffinfo->mod);
    if (success <= 0)
    {
        FLINT_ASSERT(0 && "spurious failure");
        success = 0;
        goto cleanup;
    }

    for (i = 0; i < r; i++)
    {
        success = n_polyu2_add_zip_must_match(Z + i, Ceval + i, cur_zip_image);
        if (!success)
        {
            FLINT_ASSERT(0 && "spurious failure");
            success = 0;
            goto cleanup;
        }
    }

    cur_zip_image++;
    if (cur_zip_image < req_zip_images)
        goto next_zip_image;

    for (i = 0; i < r; i++)
    {
        success = _fmpz_mpoly_modpk_update_zip(pk, B + i, ctx,
                                                    Z + i, H + i, M + i, ctxp);
        if (success <= 0)
        {
            FLINT_ASSERT(0 && "spurious failure");
            success = 0;
            goto cleanup;
        }
    }

    goto next_power;

cleanup:

    fmpz_clear(pk);

    fmpz_mpoly_clear(e, ctx);
    fmpz_mpoly_clear(t1, ctx);
    fmpz_mpoly_clear(t2, ctx);

    nmod_mpoly_clear(Tp, ctxp);
    nmod_mpolyu_clear(Tu, ctxp);
    n_polyun_clear(Teh);

    flint_free(Cdegs1);

    flint_free(betas);

    for (i = 0; i < r; i++)
    {
        n_polyun_clear(Beh + i);
        n_polyun_clear(H + i);
        n_polyun_clear(M + i);
        n_polyun_clear(Z + i);
        n_bpoly_clear(Ceval + i);
        n_bpoly_clear(Beval + i);
    }
    flint_free(Beh);
    flint_free(H);
    flint_free(M);
    flint_free(Z);
    flint_free(Ceval);
    flint_free(Beval);

    n_bpoly_clear(Teval);

    FLINT_ASSERT(success == 0 || success == 1);
    return success;
}


int fmpz_mpoly_is_fmpz_poly(
    const fmpz_mpoly_t A,
    slong var,
    const fmpz_mpoly_ctx_t ctx)
{
    int ret = 1;
    slong i, j;
    slong N = mpoly_words_per_exp(A->bits, ctx->minfo);
    fmpz * t;
    TMP_INIT;

    TMP_START;
    t = (fmpz *) TMP_ALLOC(ctx->minfo->nvars*sizeof(fmpz));
    for (i = 0; i < ctx->minfo->nvars; i++)
        fmpz_init(t + i);

    for (i = 0; i < A->length; i++)
    {
        mpoly_get_monomial_ffmpz(t, A->exps + N*i, A->bits, ctx->minfo);
        for (j = 0; j < ctx->minfo->nvars; j++)
        {
            if (j != var && !fmpz_is_zero(t + j))
            {
                ret = 0;
                goto cleanup;
            }
        }
    }

cleanup:

    for (i = 0; i < ctx->minfo->nvars; i++)
        fmpz_clear(t + i);

    TMP_END;

    return ret;
}


int fmpz_mpoly_get_fmpz_poly(
    fmpz_poly_t A,
    fmpz_mpoly_t B,
    slong var,
    const fmpz_mpoly_ctx_t ctx)
{
    slong Blen = B->length;
    const fmpz * Bcoeffs = B->coeffs;
    const ulong * Bexps = B->exps;
    flint_bitcnt_t Bbits = B->bits;
    slong i, N = mpoly_words_per_exp(Bbits, ctx->minfo);
    ulong k;

    if (B->length < 1)
    {
        fmpz_poly_zero(A);
        return 1;
    }

    if (Bbits <= FLINT_BITS)
    {
        ulong mask = (-UWORD(1)) >> (FLINT_BITS - Bbits);
        slong off, shift;

        mpoly_gen_offset_shift_sp(&off, &shift, var, Bbits, ctx->minfo);

        fmpz_poly_zero(A);
        for (i = 0; i < Blen; i++)
        {
            k = (Bexps[N*i + off] >> shift) & mask;
            fmpz_poly_set_coeff_fmpz(A, k, Bcoeffs + i);
        }
        return 1;
    }
    else
    {
        slong j, off;
        ulong check, wpf = Bbits/FLINT_BITS;

        off = mpoly_gen_offset_mp(var, Bbits, ctx->minfo);

        fmpz_poly_zero(A);
        for (i = 0; i < Blen; i++)
        {
            k = Bexps[N*i + off + 0];
            check = 0;
            for (j = 1; j < wpf; j++)
                check |= Bexps[N*i + off + j];

            if (check != 0 || (slong) k < 0)
                return 0;

            fmpz_poly_set_coeff_fmpz(A, k, Bcoeffs + i);
        }
        return 1;
    }
}

void fmpz_mpoly_fit_length_set_bits(
    fmpz_mpoly_t A,
    slong len,
    flint_bitcnt_t bits,
    const fmpz_mpoly_ctx_t ctx)
{
    slong N = mpoly_words_per_exp(bits, ctx->minfo);
    slong new_alloc;

    FLINT_ASSERT(len >= 0);

    if (A->alloc < len)
    {
        new_alloc = FLINT_MAX(len, 2*A->alloc);
        if (A->alloc > 0)
        {
            A->coeffs = (fmpz *) flint_realloc(A->coeffs, new_alloc*sizeof(fmpz));
            A->exps = (ulong *) flint_realloc(A->exps, new_alloc*N*sizeof(ulong));
            memset(A->coeffs + A->alloc, 0, (new_alloc - A->alloc)*sizeof(fmpz));
        }
        else
        {
            A->coeffs = (fmpz *) flint_calloc(new_alloc, sizeof(fmpz));
            A->exps   = (ulong *) flint_malloc(new_alloc*N*sizeof(ulong));
        }
        A->alloc = new_alloc;
    }
    else if (A->bits < bits)
    {
        if (A->alloc > 0)
            A->exps = (ulong *) flint_realloc(A->exps, A->alloc*N*sizeof(ulong));
    }

    A->bits = bits;
}

void _fmpz_mpoly_set_fmpz_poly(
    fmpz_mpoly_t A,
    flint_bitcnt_t Abits,
    fmpz * Bcoeffs,
    slong Blen,
    slong var,
    const fmpz_mpoly_ctx_t ctx)
{
    slong N = mpoly_words_per_exp(Abits, ctx->minfo);
    slong i, Alen;
    ulong * genexp;
    TMP_INIT;

    TMP_START;

    genexp = (ulong *) TMP_ALLOC(N*sizeof(ulong));
    if (Abits <= FLINT_BITS)
        mpoly_gen_monomial_sp(genexp, var, Abits, ctx->minfo);
    else
        mpoly_gen_monomial_offset_mp(genexp, var, Abits, ctx->minfo);

    fmpz_mpoly_fit_length_set_bits(A, 0, Abits, ctx);

    Alen = 0;
    for (i = Blen - 1; i >= 0; i--)
    {
        if (fmpz_is_zero(Bcoeffs + i))
            continue;
        fmpz_mpoly_fit_length(A, Alen + 1, ctx);
        fmpz_set(A->coeffs + Alen, Bcoeffs + i);
        if (Abits <= FLINT_BITS)
            mpoly_monomial_mul_ui(A->exps + N*Alen, genexp, N, i);
        else
            mpoly_monomial_mul_ui_mp(A->exps + N*Alen, genexp, N, i);
        Alen++;
    }
    A->length = Alen;

    TMP_END;
}

flint_bitcnt_t mpoly_gen_pow_exp_bits_required(
    slong v,
    ulong e,
    const mpoly_ctx_t mctx)
{
    return 1 + FLINT_BIT_COUNT(e); /* only lex and deg supported */
}


void fmpz_mpoly_set_fmpz_poly(
    fmpz_mpoly_t A,
    fmpz_poly_t B,
    slong v,
    const fmpz_mpoly_ctx_t ctx)
{
    flint_bitcnt_t bits;

    if (B->length < 1)
    {
        fmpz_mpoly_zero(A, ctx);
        return;
    }

    bits = mpoly_gen_pow_exp_bits_required(v, B->length - 1, ctx->minfo);
    bits = mpoly_fix_bits(bits, ctx->minfo);
    _fmpz_mpoly_set_fmpz_poly(A, bits, B->coeffs, B->length, v, ctx);
}


static int _try_tuncer(
    fmpz_mpoly_factor_t fac,
    const fmpz_mpoly_t A,
    const fmpz_mpoly_ctx_t ctx)
{
    int success;
    const slong n = ctx->minfo->nvars - 1;
    slong i, j, k;
    fmpz * alpha;
    fmpz_mpoly_struct * Aevals;
    slong * degs, * tdegs;
    fmpz_mpoly_factor_t tfac, lcAfac;
    fmpz_mpoly_t t, lcA, Acopy;
    fmpz_mpoly_struct * newA;
    fmpz_poly_t Au;
    fmpz_poly_factor_t Aufac;
    slong alpha_bits, alpha_count;
    flint_rand_t state;
    fmpz_mpoly_t m, mpow;
    fmpz_mpolyv_t Alc, lc_divs;
    fmpz_t q, facBound;
    mp_limb_t p;
    nmod_mpoly_ctx_t ctxp;
    nmod_mpoly_factor_t facp, tfacp;
    nmod_mpolyv_t Aevalp, Alcp;
    nmod_poly_t Aup;
    mp_limb_t * alphap;
    slong r, L;

    FLINT_ASSERT(A->bits <= FLINT_BITS);
    FLINT_ASSERT(n > 1);
/*
flint_printf("_try_tuncer called n = %wd\n", n);
flint_printf("A: "); fmpz_mpoly_print_pretty(A, NULL, ctx); printf("\n");
*/
    flint_randinit(state);

    fmpz_init(facBound);
    fmpz_init(q);

    fmpz_mpoly_init(Acopy, ctx);
    fmpz_mpoly_init(m, ctx);
    fmpz_mpoly_init(mpow, ctx);

    fmpz_mpolyv_init(Alc, ctx);
    fmpz_mpolyv_init(lc_divs, ctx);

    fmpz_poly_factor_init(Aufac);
    fmpz_poly_init(Au);

    fmpz_mpoly_init(lcA, ctx);
    fmpz_mpoly_factor_init(lcAfac, ctx);

    alpha = _fmpz_vec_init(n);
    alphap = (mp_limb_t *) flint_malloc(n*sizeof(mp_limb_t));

    degs  = (slong *) flint_malloc(2*(n + 1)*sizeof(slong));
    tdegs = degs + (n + 1);

    Aevals = (fmpz_mpoly_struct *) flint_malloc(n*sizeof(fmpz_mpoly_struct));
    for (i = 0; i < n; i++)
        fmpz_mpoly_init(Aevals + i, ctx);

    fmpz_mpoly_factor_init(tfac, ctx);
    fmpz_mpoly_init(t, ctx);

    nmod_mpoly_ctx_init(ctxp, n + 1, ORD_LEX, 2); /* modulus no care */
    nmod_mpoly_factor_init(facp, ctxp);
    nmod_mpoly_factor_init(tfacp, ctxp);
    nmod_mpolyv_init(Aevalp, ctxp);
    nmod_mpolyv_init(Alcp, ctxp);
    nmod_poly_init_mod(Aup, ctxp->ffinfo->mod);

    /* init done */

    _fmpz_mpoly_get_lc(lcA, A, ctx);
    success = fmpz_mpoly_factor(lcAfac, lcA, ctx);
    if (!success)
        goto cleanup;

    fmpz_mpoly_degrees_si(degs, A, ctx);

    alpha_count = 0;
    alpha_bits = 10;

next_alpha:

    /*
        choose the αi from +-[1, 2^alpha_bits].
        alpha_bits is incremented every so often until it gets too high.
    */

    alpha_count++;
    if (alpha_count >= alpha_bits)
    {
        alpha_count = 0;
        alpha_bits++;
        if (alpha_bits > FLINT_BITS/2)
        {
            success = 0;
            goto cleanup;
        }
    }

    for (i = 0; i < n; i++)
    {
        ulong l = n_randlimb(state);
        ulong mask = UWORD(1) << alpha_bits;
        if (l & mask)
            fmpz_neg_ui(alpha + i, 1 + (l & (mask - 1)));
        else
            fmpz_set_ui(alpha + i, 1 + (l & (mask - 1)));
    }
/*
usleep(1000);
printf("---------------------------\n");
printf("alpha = "); tuple_print(alpha, n);
*/
    /* ensure degrees do not drop under evaluation */
    for (i = n - 1; i >= 0; i--)
    {
        fmpz_mpoly_evaluate_one_fmpz(Aevals + i, i == n - 1 ? A : Aevals + i + 1,
                                                        i + 1, alpha + i, ctx);
        fmpz_mpoly_degrees_si(tdegs, Aevals + i, ctx);
        for (j = 0; j <= i; j++)
        {
            if (tdegs[j] != degs[j])
                goto next_alpha;
        }
    }

    /* make sure our univar is squarefree */
    success = fmpz_mpoly_get_fmpz_poly(Au, Aevals + 0, 0, ctx);
    FLINT_ASSERT(success);
    fmpz_poly_factor(Aufac, Au);
    r = Aufac->num;
    FLINT_ASSERT(r >= 1);
    for (j = 0; j < r; j++)
    {
        if (Aufac->exp[j] != 1)
            goto next_alpha;
    }

    /* if the univariate is irreducible, then A is irreducible */
    if (r < 2)
    {
        fmpz_mpoly_factor_one(fac, ctx);
        fmpz_mpoly_factor_append_ui(fac, A, 1, ctx);
        success = 1;
        goto cleanup;
    }

    if (lcAfac->num > 0)
    {
        success = _wang_lcc(lc_divs, lcAfac, Aufac, alpha, ctx);
        if (!success)
            goto next_alpha;
    }
    else
    {
        /* lcA is constant */
        fmpz_mpolyv_fit_length(lc_divs, r, ctx);
        lc_divs->length = r;
        for (i = 0; i < r; i++)
        {
            FLINT_ASSERT(Aufac->p[i].length > 0);
            fmpz_mpoly_set_fmpz(lc_divs->coeffs + i,
                          Aufac->p[i].coeffs + Aufac->p[i].length - 1, ctx);
        }
    }

    /*
        Assuming no extraneous factors, we have

            A(X, x1, ..., xn) = F1 * ... * Fr  where  r = Aufac->num

        and lead_divisor[i] is a divisor of lc_X(Fi). We also have the
        univariate factorization

            A(X, α1, ..., αn) = (c1 X^? + ... ) * ... * (cr X^? + ... )

        Set c(x1, ..., xn) = lc_X(A) and

            m(x1, ..., xn) = c/(prod_i lc_divs[i])   division is exact

        Lift the univariate factorization

            m(α)^(r-1) f(X, α) = (m(α)*lc_divs[0](α) X^? + ...) * ... *
                                 (m(α)*lc_divs[r-1](α) X^? + ...)

        against

            m(x1, ..., xn)^(r-1) f(X, x1, ..., xn)

        Note m(x1, ..., xn) is usually constant here, but it does not have to be.
    */

    FLINT_ASSERT(r > 1);
    success = fmpz_mpoly_divides(m, lcA, lc_divs->coeffs + 0, ctx);
    FLINT_ASSERT(success);
    for (i = 1; i < r; i++)
    {
        success = fmpz_mpoly_divides(m, m, lc_divs->coeffs + i, ctx);
        FLINT_ASSERT(success);
    }

    fmpz_mpoly_pow_ui(mpow, m, r - 1, ctx);
    if (fmpz_mpoly_is_one(mpow, ctx))
    {
        newA = (fmpz_mpoly_struct *) A;
    }
    else
    {
        newA = Acopy;
        fmpz_mpoly_mul(newA, A, mpow, ctx);
    }

    if (newA->bits > FLINT_BITS)
    {
        success = 0;
        goto cleanup;
    }

    fmpz_mpoly_degrees_si(degs, newA, ctx);
    _fmpz_vec_height(facBound, newA->coeffs, newA->length);
    if (!fmpz_mpoly_factor_bound_si(facBound, facBound, degs, n + 1))
    {
        success = 0;
        goto cleanup;
    }
    fmpz_mul_2exp(facBound, facBound, 1);
    fmpz_add_ui(facBound, facBound, 1);

/*
printf("mpow: "); fmpz_mpoly_print_pretty(mpow, NULL, ctx); printf("\n");

flint_printf("modified A: "); fmpz_mpoly_print_pretty(newA, NULL, ctx); printf("\n");
*/
    fmpz_mpoly_set(t, mpow, ctx);
    for (i = n - 1; i >= 0; i--)
    {
        fmpz_mpoly_evaluate_one_fmpz(t, mpow, i + 1, alpha + i, ctx);
        fmpz_mpoly_swap(t, mpow, ctx);
        fmpz_mpoly_mul(Aevals + i, Aevals + i, mpow, ctx);
    }

    fmpz_mpolyv_fit_length(Alc, (n + 1)*r, ctx);

    i = n;
    for (j = 0; j < r; j++)
    {
        fmpz_mpoly_mul(Alc->coeffs + i*r + j, lc_divs->coeffs + j, m, ctx);
    }
    for (i = n - 1; i >= 0; i--)
    {
        for (j = 0; j < r; j++)
        {
            fmpz_mpoly_evaluate_one_fmpz(Alc->coeffs + i*r + j,
                 Alc->coeffs + (i + 1)*r + j, i + 1, alpha + i, ctx);
        }
    }

    fmpz_mpoly_factor_fit_length(fac, r, ctx);
    fac->num = r;
    fmpz_one(fac->constant);
    for (i = 0; i < r; i++)
    {
        FLINT_ASSERT(fmpz_mpoly_is_fmpz(Alc->coeffs + 0*r + i, ctx));
        FLINT_ASSERT(fmpz_mpoly_length(Alc->coeffs + 0*r + i, ctx) == 1);
        FLINT_ASSERT(fmpz_divisible(Alc->coeffs[i].coeffs + 0,
                                 Aufac->p[i].coeffs + Aufac->p[i].length - 1));
        fmpz_divexact(q, Alc->coeffs[i].coeffs + 0,
                                  Aufac->p[i].coeffs + Aufac->p[i].length - 1);
        _fmpz_mpoly_set_fmpz_poly(fac->poly + i, newA->bits,
                               Aufac->p[i].coeffs, Aufac->p[i].length, 0, ctx);
        fmpz_mpoly_scalar_mul_fmpz(fac->poly + i, fac->poly + i, q, ctx);
        fmpz_one(fac->exp + i);
    }

    /*
        At this point, we have 

            Anew = m(x_1, ..., x_n)^(r-1)*A(X, x_1, ...., x_n)
            Aeval[n]   = Anew(X, x_1, ..., x_{n-1}, x_n)  (not stored)
            Aeval[n-1] = Anew(X, x_1, ..., x_{n-1}, α_n)
            ...
            Aeval[1]   = Anew(X, x_1, α_2, ..., α_n)
            Aeval[0]   = Anew(X, α_1, α_2, ..., α_n)

        and the leading coefficients for the factorization of
        Aeval[k] are in Alc[k][i] for 0 <= i < r, 0 <= k <= n

        Choose a prime p that
            1. Divides no coefficient of any Aeval[k], 0 <= k <= n
            2. Divides no coefficient of any Alc[k][i]
            3. Keeps Aeval[0] squarefree
            4. Divides no coefficient of any factor of Aeval[0]

        Try several p before going back an picking a new α.
    */

    p = UWORD(1) << (FLINT_BITS - 1);

next_tuncer_prime:

    if (p >= UWORD_MAX_PRIME)
    {
        success = 0;
        goto cleanup;
    }

    p = n_nextprime(p, 1);
/*
flint_printf("---------- next_tuncer_prime p = %wu -----------\n", p);
    usleep(1000);
*/
    nmod_mpoly_ctx_set_modulus(ctxp, p);
    nmod_mpoly_factor_fit_length(facp, r, ctxp);
    nmod_mpoly_factor_fit_length(tfacp, r, ctxp);
    facp->num = r;
    tfacp->num = r;

    nmod_mpolyv_fit_length(Aevalp, n + 1, ctxp);
    nmod_mpolyv_fit_length(Alcp, (n + 1)*r, ctxp);

    facp->constant = 1;
    for (i = 0; i < r; i++)
    {
        /* check condition 4 */
        _nmod_mpoly_set_fmpz_mpoly(facp->poly + i, ctxp, fac->poly + i, ctx);
        if (facp->poly[i].length != fac->poly[i].length)
            goto next_tuncer_prime;

        fmpz_one(facp->exp + i);
    }

    /* check condition 3 */
    Aup->mod = ctxp->ffinfo->mod;
    nmod_poly_set_fmpz_poly(Aup, Au);
    if (Aup->length != Au->length || !nmod_poly_is_squarefree(Aup))
        goto next_tuncer_prime;

    for (k = 0; k <= n; k++)
    {
        /* check condition 1 */
        _nmod_mpoly_set_fmpz_mpoly(Aevalp->coeffs + k, ctxp,
                                               k < n ? Aevals + k : newA, ctx);
        if (Aevalp->coeffs[k].length != (k < n ? Aevals + k : newA)->length)
            goto next_tuncer_prime;

        /* check condition 2 */
        for (i = 0; i < r; i++)
        {
            _nmod_mpoly_set_fmpz_mpoly(Alcp->coeffs + k*r + i, ctxp,
                                                   Alc->coeffs + k*r + i, ctx);
            if (Alcp->coeffs[k*r + i].length != Alc->coeffs[k*r + i].length)
                goto next_tuncer_prime;
        }
    }

    /* set alpha's mod p */
    for (i = 0; i < n; i++)
        alphap[i] = fmpz_fdiv_ui(alpha + i, ctxp->ffinfo->mod.n);

    for (k = 1; k <= n; k++)
    {
        /* LC fix */
        for (i = 0; i < r; i++)
        {
            fmpz_one(tfacp->exp + i);
            _nmod_mpoly_set_lc(tfacp->poly + i, facp->poly + i,
                                         Alcp->coeffs + k*r + i, ctxp);
        }
/*
flint_printf("tuncer step A[%wd]: ", k); nmod_mpoly_print_pretty(Aevalp + k, NULL, ctxp); printf("\n");
*/
        if (k >= 3)
        {
            success = nmod_mpoly_factor_lift_tuncer(r, tfacp->poly, state,
                                    alphap, Aevalp->coeffs + k, degs, k, ctxp);
            if (!success)
                goto next_tuncer_prime;
        }
        else
        {
            success = nmod_mfactor_lift(k, tfacp, alphap,
                                               Aevalp->coeffs + k, degs, ctxp);
            if (!success)
                goto next_tuncer_prime;
        }

        nmod_mpoly_factor_swap(tfacp, facp, ctxp);
    }

    /* LC fix */
    fmpz_mpoly_factor_fit_length(fac, r, ctx);
    fac->num = r;
    for (i = 0; i < r; i++)
    {
        fmpz_one(fac->exp + i);
        _fmpz_mpoly_set_nmod_mpoly_smod(fac->poly + i, ctx,
                                                         facp->poly + i, ctxp);
        _fmpz_mpoly_set_lc(fac->poly + i, fac->poly + i,
                                                   Alc->coeffs + n*r + i, ctx);
    }

    L = fmpz_clog_ui(facBound, ctxp->ffinfo->mod.n);

    if (0)
    {
        /* recursive approach */
        success = fmpz_mfactor_lift_prime_power(r, fac->poly,
                                 facp->poly, newA, degs, alphap, ctx, ctxp, L);
        FLINT_ASSERT(success == 1 && "don't know want to do");
    }
    else
    {
        /* more zippel :-( */
        success = fmpz_mfactor_lift_prime_power_tuncer(r, fac->poly, state,
                                       facp->poly, newA, alphap, ctx, ctxp, L);
        FLINT_ASSERT(success == 1 && "don't know want to do");
    }

    success = 1;

    if (fmpz_mpoly_is_fmpz(m, ctx))
    {
        for (i = 0; i < Aufac->num; i++)
        {
            _fmpz_vec_content(q, fac->poly[i].coeffs, fac->poly[i].length);
            if (fmpz_sgn(fac->poly[i].coeffs + 0) < 0)
                fmpz_neg(q, q);
            fmpz_mpoly_scalar_divexact_fmpz(fac->poly + i, fac->poly + i, q, ctx);
        }
    }
    else
    {
        fmpz_mpoly_univar_t u;
        fmpz_mpoly_univar_init(u, ctx);
        for (i = 0; i < Aufac->num; i++)
        {
            fmpz_mpoly_to_univar(u, fac->poly + i, 0, ctx);
            success = fmpz_mpoly_univar_content_mpoly(t, u, ctx);
            if (!success)
            {
                fmpz_mpoly_univar_clear(u, ctx);
                goto cleanup;
            }
            success = fmpz_mpoly_divides(fac->poly + i, fac->poly + i, t, ctx);
            FLINT_ASSERT(success);
            FLINT_ASSERT(fac->poly[i].length > 0);
            if (fmpz_sgn(fac->poly[i].coeffs + 0) < 0)
                fmpz_mpoly_neg(fac->poly + i, fac->poly + i, ctx);
        }
        fmpz_mpoly_univar_clear(u, ctx);
    }

cleanup:

    flint_randclear(state);

    fmpz_clear(facBound);
    fmpz_clear(q);

    fmpz_mpoly_clear(Acopy, ctx);
    fmpz_mpoly_clear(m, ctx);
    fmpz_mpoly_clear(mpow, ctx);

    fmpz_mpolyv_clear(Alc, ctx);
    fmpz_mpolyv_clear(lc_divs, ctx);

    fmpz_poly_factor_clear(Aufac);
    fmpz_poly_clear(Au);

    fmpz_mpoly_clear(lcA, ctx);
    fmpz_mpoly_factor_clear(lcAfac, ctx);

    _fmpz_vec_clear(alpha, n);
    flint_free(alphap);
    flint_free(degs);

    for (i = 0; i < n; i++)
        fmpz_mpoly_clear(Aevals + i, ctx);
    flint_free(Aevals);

    fmpz_mpoly_factor_clear(tfac, ctx);
    fmpz_mpoly_clear(t, ctx);

    nmod_mpoly_factor_clear(facp, ctxp);
    nmod_mpoly_factor_clear(tfacp, ctxp);
    nmod_mpolyv_clear(Aevalp, ctxp);
    nmod_mpolyv_clear(Alcp, ctxp);
    nmod_poly_clear(Aup);
    nmod_mpoly_ctx_clear(ctxp);

/*
flint_printf("_try_tuncer returning %d\n", success);
*/
    return success;
}




int _try_wang(
    fmpz_mpoly_factor_t fac,
    const fmpz_mpoly_t A,
    const fmpz_mpoly_ctx_t ctx)
{
    int success;
    const slong n = ctx->minfo->nvars - 1;
    slong i, j, k;
    fmpz * alpha;
    fmpz_mpoly_struct * Aevals;
    slong * degs, * degeval;
    fmpz_mpoly_factor_t tfac, lcAfac;
    fmpz_mpoly_t t, lcA, Acopy;
    fmpz_mpoly_struct * newA;
    fmpz_poly_t Au;
    fmpz_poly_factor_t Aufac;
    slong alpha_modulus, alpha_count;
    flint_rand_t state;
    fmpz_mpoly_t m, mpow;
    fmpz_mpolyv_t new_lcs, lc_divs;
    fmpz_t q;
    slong r;

    FLINT_ASSERT(A->bits <= FLINT_BITS);
    FLINT_ASSERT(n > 1);
/*
flint_printf("_try_wang called n = %wd\n", n);
flint_printf("A: "); fmpz_mpoly_print_pretty(A, NULL, ctx); flint_printf("\n");
*/
    fmpz_init(q);

    flint_randinit(state);

    fmpz_mpoly_init(Acopy, ctx);
    fmpz_mpoly_init(m, ctx);
    fmpz_mpoly_init(mpow, ctx);

    fmpz_mpolyv_init(new_lcs, ctx);
    fmpz_mpolyv_init(lc_divs, ctx);

    fmpz_poly_factor_init(Aufac);
    fmpz_poly_init(Au);

    fmpz_mpoly_init(lcA, ctx);
    fmpz_mpoly_factor_init(lcAfac, ctx);
    alpha = _fmpz_vec_init(n);

    degs    = (slong *) flint_malloc((n + 1)*sizeof(slong));
    degeval = (slong *) flint_malloc((n + 1)*sizeof(slong));
    Aevals    = (fmpz_mpoly_struct *) flint_malloc(n*sizeof(fmpz_mpoly_struct));
    for (i = 0; i < n; i++)
        fmpz_mpoly_init(Aevals + i, ctx);

    fmpz_mpoly_factor_init(tfac, ctx);
    fmpz_mpoly_init(t, ctx);

    _fmpz_mpoly_get_lc(lcA, A, ctx);
    success = fmpz_mpoly_factor(lcAfac, lcA, ctx);
    if (!success)
        goto cleanup;

    fmpz_mpoly_degrees_si(degs, A, ctx);

    alpha_count = 0;
    alpha_modulus = 1;
    goto got_alpha;

next_alpha:

    alpha_count++;
    if (alpha_count >= alpha_modulus)
    {
        alpha_count = 0;
        alpha_modulus++;
        if (alpha_modulus > 100)
        {
            success = 0;
            goto cleanup;
        }
    }

    for (i = 0; i < n; i++)
        fmpz_set_si(alpha + i, n_urandint(state, alpha_modulus) - alpha_modulus/2);

got_alpha:
/*
usleep(1000);
flint_printf("---------------------------\n");
flint_printf("(mod %wd) alpha = ", alpha_modulus); tuple_print(alpha, n);
*/
    /* ensure degrees do not drop under evaluation */
    for (i = n - 1; i >= 0; i--)
    {
        fmpz_mpoly_evaluate_one_fmpz(Aevals + i, i == n - 1 ? A : Aevals + i + 1, i + 1, alpha + i, ctx);
        fmpz_mpoly_degrees_si(degeval, Aevals + i, ctx);
        for (j = 0; j <= i; j++)
        {
            if (degeval[j] != degs[j])
                goto next_alpha;
        }
    }

    /* make sure our univar is squarefree */
    FLINT_ASSERT(fmpz_mpoly_is_fmpz_poly(Aevals + 0, 0, ctx));
    success = fmpz_mpoly_get_fmpz_poly(Au, Aevals + 0, 0, ctx);
    FLINT_ASSERT(success);
    fmpz_poly_factor(Aufac, Au);
    r = Aufac->num;
    for (j = 0; j < r; j++)
    {
        if (Aufac->exp[j] != 1)
            goto next_alpha;
    }

    /* if the univariate is irreducible, then A irreducible */
    if (r < 2)
    {
        fmpz_mpoly_factor_one(fac, ctx);
        fmpz_mpoly_factor_append_ui(fac, A, 1, ctx);
        success = 1;
        goto cleanup;
    }

    if (lcAfac->num > 0)
    {
        if (!_wang_lcc(lc_divs, lcAfac, Aufac, alpha, ctx))
            goto next_alpha;
    }
    else
    {
        /* lcA is constant */
        fmpz_mpolyv_fit_length(lc_divs, Aufac->num, ctx);
        lc_divs->length = Aufac->num;
        for (i = 0; i < Aufac->num; i++)
        {
            FLINT_ASSERT(Aufac->p[i].length > 0);
            fmpz_mpoly_set_fmpz(lc_divs->coeffs + i,
                             Aufac->p[i].coeffs + Aufac->p[i].length - 1, ctx);
        }
    }

    /*
        Assuming no extraneous divisors, we have

            A(X, x1, ..., xn) = F1 * ... * Fr  where  r = Aufac->num

        and lead_divisor[i] is a divisor of lc_X(Fi). We also have the
        univariate factorization

            A(X, α1, ..., αn) = (c1 X^? + ... ) * ... * (cr X^? + ... )

        Set c(x1, ..., xn) = lc_X(A) and

            m(x1, ..., xn) = c/(prod_i lc_divs[i])   division is exact

        Lift the univariate factorization

            m(α)^(r-1) f(X, α) = (m(α)*lc_divs[0](α) X^? + ...) * ... * (m(α)*lc_divs[r-1](α) X^? + ...)

        against

            m(x1, ..., xn)^(r-1) f(X, x1, ..., xn)

        Note m(x1, ..., xn) is usually constant here, but it certainly does not have to be.
    */

    FLINT_ASSERT(r > 1);
    success = fmpz_mpoly_divides(m, lcA, lc_divs->coeffs + 0, ctx);
    FLINT_ASSERT(success);
    for (i = 1; i < Aufac->num; i++)
    {
        success = fmpz_mpoly_divides(m, m, lc_divs->coeffs + i, ctx);
        FLINT_ASSERT(success);
    }
/*
printf("m: "); fmpz_mpoly_print_pretty(m, NULL, ctx); printf("\n");
*/
    fmpz_mpoly_pow_ui(mpow, m, r - 1, ctx);
    if (fmpz_mpoly_is_one(mpow, ctx))
    {
        newA = (fmpz_mpoly_struct *) A;
    }
    else
    {
        newA = Acopy;
        fmpz_mpoly_mul(newA, A, mpow, ctx);
    }

    fmpz_mpoly_degrees_si(degs, newA, ctx);

/*
printf("mpow: "); fmpz_mpoly_print_pretty(mpow, NULL, ctx); printf("\n");
flint_printf("modified A: "); fmpz_mpoly_print_pretty(newA, NULL, ctx); printf("\n");
*/
    fmpz_mpoly_set(t, mpow, ctx);
    for (i = n - 1; i >= 0; i--)
    {
        fmpz_mpoly_evaluate_one_fmpz(t, mpow, i + 1, alpha + i, ctx);
        fmpz_mpoly_swap(t, mpow, ctx);
        fmpz_mpoly_mul(Aevals + i, Aevals + i, mpow, ctx);
/*
flint_printf("Aeval[%wd]: ", i);  fmpz_mpoly_print_pretty(Aevals + i, NULL, ctx); printf("\n");
*/
    }

    fmpz_mpolyv_fit_length(new_lcs, (n + 1)*Aufac->num, ctx);

    i = n;
    for (j = 0; j < Aufac->num; j++)
    {
        fmpz_mpoly_mul(new_lcs->coeffs + i*Aufac->num + j,
                                            lc_divs->coeffs + j, m, ctx);
    }
    for (i = n - 1; i >= 0; i--)
    {
        for (j = 0; j < Aufac->num; j++)
        {
            fmpz_mpoly_evaluate_one_fmpz(new_lcs->coeffs + i*Aufac->num + j,
                 new_lcs->coeffs + (i + 1)*Aufac->num + j, i + 1, alpha + i, ctx);
        }
    }

    fmpz_mpoly_factor_fit_length(fac, r, ctx);
    fac->num = r;
    fmpz_one(fac->constant);
    for (i = 0; i < Aufac->num; i++)
    {
        FLINT_ASSERT(fmpz_mpoly_is_fmpz(new_lcs->coeffs + 0*r + i, ctx));
        FLINT_ASSERT(fmpz_mpoly_length(new_lcs->coeffs + 0*r + i, ctx) == 1);
        FLINT_ASSERT(fmpz_divisible(new_lcs->coeffs[i].coeffs + 0, Aufac->p[i].coeffs + Aufac->p[i].length - 1));
        fmpz_divexact(q, new_lcs->coeffs[i].coeffs + 0,
                                  Aufac->p[i].coeffs + Aufac->p[i].length - 1);
        _fmpz_mpoly_set_fmpz_poly(fac->poly + i, newA->bits,
                               Aufac->p[i].coeffs, Aufac->p[i].length, 0, ctx);
        fmpz_mpoly_scalar_mul_fmpz(fac->poly + i, fac->poly + i, q, ctx);
        fmpz_one(fac->exp + i);
    }

    fmpz_mpoly_factor_fit_length(tfac, Aufac->num, ctx);
    tfac->num = r;
    fmpz_one(tfac->constant);

    for (k = 1; k <= n; k++)
    {
        fmpz_mpoly_struct * target;
        for (i = 0; i < r; i++)
        {
            fmpz_one(tfac->exp + i);
            _fmpz_mpoly_set_lc(tfac->poly + i, fac->poly + i,
                                               new_lcs->coeffs + k*r + i, ctx);
        }

        target = k < n ? Aevals + k : newA;

/*
flint_printf("wang lifting variable %wd\n", k);
        if (k > 2)
            play_lift(k, tfac, alpha, target, degs, ctx);
*/

        if (tfac->num == 2)
            success = _mfactor_lift_quartic2(k, tfac, alpha, target, degs, ctx);
        else if (tfac->num <  20)
            success = _mfactor_lift_quartic(k, tfac, alpha, target, degs, ctx);
        else
            success = _mfactor_lift_quintic(k, tfac, alpha, target, degs, ctx);

        if (!success)
            goto next_alpha;

        fmpz_mpoly_factor_swap(tfac, fac, ctx);
    }

    success = 1;

    if (fmpz_mpoly_is_fmpz(m, ctx))
    {
        for (i = 0; i < r; i++)
        {
            _fmpz_vec_content(q, fac->poly[i].coeffs, fac->poly[i].length);
            if (fmpz_sgn(fac->poly[i].coeffs + 0) < 0)
                fmpz_neg(q, q);
            fmpz_mpoly_scalar_divexact_fmpz(fac->poly + i, fac->poly + i, q, ctx);
        }
    }
    else
    {
        fmpz_mpoly_univar_t u;
        fmpz_mpoly_univar_init(u, ctx);
        for (i = 0; i < r; i++)
        {
            fmpz_mpoly_to_univar(u, fac->poly + i, 0, ctx);
            success = fmpz_mpoly_univar_content_mpoly(t, u, ctx);
            if (!success)
            {
                fmpz_mpoly_univar_clear(u, ctx);
                goto cleanup;
            }
            success = fmpz_mpoly_divides(fac->poly + i, fac->poly + i, t, ctx);
            FLINT_ASSERT(success);
            FLINT_ASSERT(fac->poly[i].length > 0);
            if (fmpz_sgn(fac->poly[i].coeffs + 0) < 0)
                fmpz_mpoly_neg(fac->poly + i, fac->poly + i, ctx);
        }
        fmpz_mpoly_univar_clear(u, ctx);
    }

cleanup:

    fmpz_clear(q);

    fmpz_mpolyv_clear(new_lcs, ctx);
    fmpz_mpolyv_clear(lc_divs, ctx);

    fmpz_poly_factor_clear(Aufac);

    fmpz_mpoly_clear(lcA, ctx);
    fmpz_mpoly_factor_clear(lcAfac, ctx);
        
    _fmpz_vec_clear(alpha, n);

    for (i = 0; i < n; i++)
        fmpz_mpoly_clear(Aevals + i, ctx);
    flint_free(Aevals);

    flint_free(degs);
    flint_free(degeval);
    fmpz_mpoly_factor_clear(tfac, ctx);
    fmpz_mpoly_clear(t, ctx);

    fmpz_mpoly_clear(Acopy, ctx);
    fmpz_mpoly_clear(m, ctx);
    fmpz_mpoly_clear(mpow, ctx);

    fmpz_poly_clear(Au);

/*
flint_printf("wang returning %d\n", success);
*/
    return success;
}



/* A is square free and primitive w.r.t all variables */
static int _irreducible_mvar_factors(
    fmpz_mpoly_factor_t fac,
    fmpz_mpoly_t A,
    const fmpz_mpoly_ctx_t ctx)
{
    int success;
    const slong n = ctx->minfo->nvars - 1;
    slong i, j, m, r;
    fmpz_t subset;
    fmpz * alpha, * alphait;
    fmpz_mpoly_struct * Aevals;
    slong * deg, * degeval;
    fmpz_mpoly_factor_t qfac, pfac, tfac, dfac;
    fmpz_mpoly_t t, p, q;
    fmpz_mpoly_univar_t u;

    FLINT_ASSERT(n > 1);

    if (fmpz_sgn(A->coeffs + 0) < 0)
        fmpz_mpoly_neg(A, A, ctx);
/*
flint_printf("_irreducible_mvar_factors(n = %wd) called\nA: ", n); fmpz_mpoly_print_pretty(A, NULL, ctx); printf("\n");
*/

#if 0

/*flint_printf("trying wang\n");*/
    if (_try_wang(fac, A, ctx))
    {
/*flint_printf("wang success\n");*/
        return 1;
    }
flint_printf("wang failed\n");

#else

/*flint_printf("trying tuncer\n");*/
    if (_try_tuncer(fac, A, ctx))
    {
/*flint_printf("truncer success\n");*/
        return 1;
    }
flint_printf("truncer failed\n");

#endif

    fmpz_init(subset);
    alphait = _fmpz_vec_init(n);
    alpha   = _fmpz_vec_init(n);
    Aevals  = (fmpz_mpoly_struct *) flint_malloc(n*sizeof(fmpz_mpoly_struct));
    deg     = (slong *) flint_malloc((n + 1)*sizeof(slong));
    degeval = (slong *) flint_malloc((n + 1)*sizeof(slong));
    for (i = 0; i < n; i++)
        fmpz_mpoly_init(Aevals + i, ctx);
    fmpz_mpoly_factor_init(pfac, ctx);
    fmpz_mpoly_factor_init(qfac, ctx);
    fmpz_mpoly_factor_init(tfac, ctx);
    fmpz_mpoly_factor_init(dfac, ctx);
    fmpz_mpoly_init(t, ctx);
    fmpz_mpoly_init(p, ctx);
    fmpz_mpoly_init(q, ctx);
    fmpz_mpoly_univar_init(u, ctx);

    fmpz_mpoly_factor_one(fac, ctx);

    fmpz_mpoly_degrees_si(deg, A, ctx);
    goto got_alpha;

next_alpha:

    tuple_next(alphait, n);
    for (i = 0; i < n; i++)
    {
        j = n - 1 - i;
        fmpz_cdiv_q_2exp(alpha + j, alphait + i, 1);
        if (fmpz_is_even(alphait + i))
            fmpz_neg(alpha + j, alpha + j);
    }

got_alpha:

/*
usleep(1000000);
printf("alpha = "); tuple_print(alpha, n);
*/

    /* ensure degrees do not drop under evaluation */
    for (i = n - 1; i >= 0; i--)
    {
        fmpz_mpoly_evaluate_one_fmpz(Aevals + i, i == n - 1 ? A : Aevals + i + 1,
                                                        i + 1, alpha + i, ctx);
        fmpz_mpoly_degrees_si(degeval, Aevals + i, ctx);
        for (j = 0; j <= i; j++)
        {
            if (degeval[j] != deg[j])
            {
                tuple_saturate(alphait, n, n - i);
                goto next_alpha;
            }
        }
    }

    /* make sure our univar is squarefree */
    fmpz_mpoly_derivative(t, Aevals + 0, 0, ctx);
    fmpz_mpoly_gcd(t, t, Aevals + 0, ctx);
    if (!fmpz_mpoly_is_fmpz(t, ctx))
        goto next_alpha;

    /* make our evaluations primitive */
    for (i = n - 1; i > 0; i--)
    {
        fmpz_mpoly_to_univar(u, Aevals + i, 0, ctx);
        fmpz_mpoly_univar_content_mpoly(t, u, ctx);
        success = fmpz_mpoly_divides(Aevals + i, Aevals + i, t, ctx);
        FLINT_ASSERT(success);
        FLINT_ASSERT(Aevals[i].length > 0);
        if (fmpz_sgn(Aevals[i].coeffs + 0) < 0)
            fmpz_mpoly_neg(Aevals + i, Aevals + i, ctx);
    }

    _irreducible_bivar_factors(pfac, Aevals + 1, 0, 1, ctx);

    for (m = 2; m <= n; m++)
    {
        fmpz_mpoly_set(q, m < n ? Aevals + m : A, ctx);
        fmpz_mpoly_set(p, Aevals + m - 1, ctx);

        FLINT_ASSERT(fmpz_is_one(pfac->constant));
        FLINT_ASSERT(fmpz_mpoly_factor_matches(p, pfac, ctx));

        /* if pfac has only one factor, A must be irreducible */
        if (pfac->num < 2)
        {
            fmpz_mpoly_factor_append_ui(fac, A, 1, ctx);
            success = 1;
            goto cleanup;
        }

        success = _try_lift(qfac, q, pfac, p, m, alpha, n, ctx);
        if (success)
        {
            fmpz_mpoly_factor_swap(qfac, pfac, ctx);
            continue;
        }

        /* if we couldn't lift two local factors, A must be irreducible */
        if (pfac->num == 2)
        {
            fmpz_mpoly_factor_append_ui(fac, A, 1, ctx);
            success = 1;
            goto cleanup;
        }

        qfac->num = 0;

try_again:

        for (r = 1; r <= pfac->num/2; r++)
        {
            subset_first(subset, pfac->num, r);

            FLINT_ASSERT(fmpz_is_one(pfac->constant));
            FLINT_ASSERT(fmpz_mpoly_factor_matches(p, pfac, ctx));

            do {
                fmpz_mpoly_factor_fit_length(dfac, 2, ctx);
                dfac->num = 2;
                fmpz_one(dfac->constant);
                fmpz_mpoly_one(dfac->poly + 0, ctx);
                fmpz_mpoly_one(dfac->poly + 1, ctx);
                fmpz_one(dfac->exp + 0);
                fmpz_one(dfac->exp + 1);
                for (i = 0; i < pfac->num; i++)
                {
                    j = fmpz_tstbit(subset, i);
                    fmpz_mpoly_mul(dfac->poly + j, dfac->poly + j,
                                                          pfac->poly + i, ctx);
                }

                success = _try_lift(tfac, q, dfac, p, m, alpha, n, ctx);
                if (success)
                {
                    for (i = pfac->num - 1; i >= 0; i--)
                    {
                        if (fmpz_tstbit(subset, i))
                        {
                            fmpz_mpoly_swap(pfac->poly + i,
                                              pfac->poly + pfac->num - 1, ctx);
                            pfac->num--;
                        }
                    }
                    fmpz_mpoly_factor_append_ui(qfac, tfac->poly + 1, 1, ctx);
                    fmpz_mpoly_swap(q, tfac->poly + 0, ctx);
                    fmpz_mpoly_swap(p, dfac->poly + 0, ctx);
                    goto try_again;
                }
            }
            while (subset_next(subset, subset, pfac->num));
        }
        /* if pfac could not be combined, p must be irreducible */
        fmpz_mpoly_factor_append_ui(qfac, q, 1, ctx);
        fmpz_mpoly_factor_swap(qfac, pfac, ctx);
    }

    success = 1;

    for (i = 0; i < pfac->num; i++)
        fmpz_mpoly_factor_append_ui(fac, pfac->poly + i, 1, ctx);

cleanup:

    fmpz_clear(subset);
    _fmpz_vec_clear(alphait, n);
    _fmpz_vec_clear(alpha, n);
    for (i = 0; i < n; i++)
        fmpz_mpoly_clear(Aevals + i, ctx);
    flint_free(Aevals);
    flint_free(deg);
    flint_free(degeval);
    fmpz_mpoly_factor_clear(pfac, ctx);
    fmpz_mpoly_factor_clear(qfac, ctx);
    fmpz_mpoly_factor_clear(tfac, ctx);
    fmpz_mpoly_factor_clear(dfac, ctx);
    fmpz_mpoly_clear(t, ctx);
    fmpz_mpoly_clear(p, ctx);
    fmpz_mpoly_clear(q, ctx);
    fmpz_mpoly_univar_clear(u, ctx);

    FLINT_ASSERT(!success || (fmpz_is_one(fac->constant) &&
                              fmpz_mpoly_factor_matches(A, fac, ctx)));

    return success;
}


/*
    A is primitive w.r.t to any variable appearing in A.
    A is square free with positive lead coeff.

    return 1 for success, 0 for failure
*/
static int _irreducible_factors(
    fmpz_mpoly_factor_t f,
    const fmpz_mpoly_t A,
    const fmpz_mpoly_ctx_t ctx)
{
    int success;
    slong i, mvars;
    slong * Adegs, * perm, * iperm;
    ulong * shift, * stride;
    TMP_INIT;

    if (A->bits > FLINT_BITS)
        return 0;

    TMP_START;

    Adegs = (slong *) TMP_ALLOC(ctx->minfo->nvars*sizeof(slong));
    perm = (slong *) TMP_ALLOC(ctx->minfo->nvars*sizeof(slong));
    iperm = (slong *) TMP_ALLOC(ctx->minfo->nvars*sizeof(slong));
    shift = (ulong *) TMP_ALLOC(ctx->minfo->nvars*sizeof(slong));
    stride = (ulong *) TMP_ALLOC(ctx->minfo->nvars*sizeof(slong));
    fmpz_mpoly_degrees_si(Adegs, A, ctx);

    mvars = 0;

    for (i = 0; i < ctx->minfo->nvars; i++)
    {
        shift[i] = 0;
        stride[i] = 1;
        iperm[i] = -1;
        if (Adegs[i] > 0)
        {
            perm[mvars] = i;
            mvars++;
        }
    }
    /* TODO: figure out nice perm */


    /* invert perm */
    for (i = 0; i < mvars; i++)
        iperm[perm[i]] = i;

    if (mvars == 1)
    {
        fmpz_mpoly_t t;
        fmpz_poly_t Au;
        fmpz_poly_factor_t Aufac;

        fmpz_mpoly_init(t, ctx);
        fmpz_poly_init(Au);
        fmpz_poly_factor_init(Aufac);

        _fmpz_mpoly_to_fmpz_poly_deflate(Au, A, perm[0], shift, stride, ctx);
        fmpz_poly_factor(Aufac, Au);

        fmpz_mpoly_factor_one(f, ctx);
        fmpz_mul(f->constant, f->constant, &Aufac->c); /* Aufac->c should be 1 */
        for (i = 0; i < Aufac->num; i++)
        {
            _fmpz_mpoly_from_fmpz_poly_inflate(t, A->bits, Aufac->p + i, perm[0], shift, stride, ctx);
            fmpz_mpoly_factor_append_ui(f, t, Aufac->exp[i], ctx); /* Aufac->exp[i] should be 1 */
        }

        fmpz_mpoly_clear(t, ctx);
        fmpz_poly_clear(Au);
        fmpz_poly_factor_clear(Aufac);

        success = 1;
    }
    else if (mvars == 2)
    {
        success = _irreducible_bivar_factors(f, A, perm[0], perm[1], ctx);
        fmpz_mpoly_factor_fix_units(f, ctx);
    }
    else
    {
        fmpz_mpoly_ctx_t lctx;
        fmpz_mpoly_t Al, B;
        fmpz_mpoly_factor_t Amfactors;

        fmpz_mpoly_ctx_init(lctx, mvars, ORD_LEX);
        fmpz_mpoly_init(Al, lctx);
        fmpz_mpoly_init(B, ctx);
        fmpz_mpoly_factor_init(Amfactors, lctx);

        fmpz_mpoly_convert_perm(Al, A->bits, lctx, A, ctx, perm);

        success = _irreducible_mvar_factors(Amfactors, Al, lctx);
        FLINT_ASSERT(success);
        if (success)
        {
            fmpz_mpoly_factor_one(f, ctx);
            for (i = 0; i < Amfactors->num; i++)
            {
                fmpz_mpoly_convert_perm(B, A->bits, ctx, Amfactors->poly + i, lctx, iperm);
                FLINT_ASSERT(B->length > 0);
                if (fmpz_sgn(B->coeffs + 0) < 0)
                    fmpz_mpoly_neg(B, B, ctx);
                fmpz_mpoly_factor_append_fmpz(f, B, Amfactors->exp + i, ctx);
            }
        }

        fmpz_mpoly_factor_clear(Amfactors, lctx);
        fmpz_mpoly_clear(B, ctx);
        fmpz_mpoly_clear(Al, lctx);
        fmpz_mpoly_ctx_clear(lctx);
    }

    TMP_END;
/*
flint_printf("success: %d\n", success);
flint_printf("A: "); fmpz_mpoly_print_pretty(A, NULL, ctx); flint_printf("\n");
flint_printf("f: "); fmpz_mpoly_factor_print_pretty(f, NULL, ctx); flint_printf("\n");
*/
    FLINT_ASSERT(!success || fmpz_mpoly_factor_matches(A, f, ctx));

    return success;
}


/*
    A is primitive w.r.t each variable appearing in it
    return 1 for success, 0 for failure
*/
static int _squarefree_factors(
    fmpz_mpoly_factor_t f,
    const fmpz_mpoly_t A,
    const fmpz_mpoly_ctx_t ctx)
{
    int success;
    slong v;
    fmpz_t k;
    fmpz_mpoly_t S, Sp, Sm, Ss, Y, Z;

    fmpz_init(k);
    fmpz_mpoly_init(S, ctx);
    fmpz_mpoly_init(Sp, ctx);
    fmpz_mpoly_init(Sm, ctx);
    fmpz_mpoly_init(Ss, ctx);
    fmpz_mpoly_init(Y, ctx);
    fmpz_mpoly_init(Z, ctx);

    fmpz_mpoly_factor_one(f, ctx);

    for (v = 0; v < ctx->minfo->nvars; v++)
    {
        fmpz_mpoly_derivative(Sp, A, v, ctx);

        if (!fmpz_mpoly_is_zero(Sp, ctx))
        {
            success = fmpz_mpoly_gcd_cofactors(Sm, Ss, Y, A, Sp, ctx);
            if (!success)
                goto cleanup;

            for (fmpz_set_ui(k, 1); !(fmpz_mpoly_derivative(Sp, Ss, v, ctx),
                                      fmpz_mpoly_sub(Z, Y, Sp, ctx),
                                      fmpz_mpoly_is_zero(Z, ctx));
                                                         fmpz_add_ui(k, k, 1))
            {
                success = fmpz_mpoly_gcd_cofactors(S, Ss, Y, Ss, Z, ctx);
                if (!success)
                    goto cleanup;

                fmpz_mpoly_factor_mul_mpoly_fmpz(f, S, k, ctx);
            }

            fmpz_mpoly_factor_mul_mpoly_fmpz(f, Ss, k, ctx);

            success = 1;
            goto cleanup;
        }
    }

    FLINT_ASSERT(fmpz_mpoly_is_fmpz(A, ctx));

    fmpz_mpoly_factor_mul_mpoly_ui(f, A, 1, ctx);

    success = 1;

cleanup:

    fmpz_clear(k);
    fmpz_mpoly_clear(S, ctx);
    fmpz_mpoly_clear(Sp, ctx);
    fmpz_mpoly_clear(Sm, ctx);
    fmpz_mpoly_clear(Ss, ctx);
    fmpz_mpoly_clear(Y, ctx);
    fmpz_mpoly_clear(Z, ctx);

    return success;
}


int fmpz_mpoly_factor_squarefree(
    fmpz_mpoly_factor_t f,
    const fmpz_mpoly_t A,
    const fmpz_mpoly_ctx_t ctx)
{
    int success;
    slong j, v;
    fmpz_mpoly_t c;
    fmpz_mpoly_univar_t u;
    fmpz_mpoly_factor_t newf, tempf;
    fmpz * var_powers;

    /* 0. set trivial factorization */
    fmpz_mpoly_factor_one(f, ctx);
    if (fmpz_mpoly_is_fmpz(A, ctx))
    {
        fmpz_mpoly_get_fmpz(f->constant, A, ctx);
        return 1;
    }
    else
    {
        FLINT_ASSERT(A->length > 0);
        _fmpz_vec_content(f->constant, A->coeffs, A->length);
        if (fmpz_sgn(A->coeffs + 0) < 0)
            fmpz_neg(f->constant, f->constant);

        fmpz_mpoly_factor_fit_length(f, 1, ctx);
        fmpz_mpoly_scalar_divexact_fmpz(f->poly + 0, A, f->constant, ctx);
        fmpz_one(f->exp + 0);
        f->num = 1;
    }

    fmpz_mpoly_factor_init(newf, ctx);
    fmpz_mpoly_factor_init(tempf, ctx);
    fmpz_mpoly_univar_init(u, ctx);
    fmpz_mpoly_init(c, ctx);
    var_powers = _fmpz_vec_init(ctx->minfo->nvars);

    /* 1. ensure factors are primitive w.r.t any variable */
    FLINT_ASSERT(fmpz_mpoly_factor_matches(A, f, ctx));

    for (v = 0; v < ctx->minfo->nvars; v++)
    {
        fmpz_set(newf->constant, f->constant);
        newf->num = 0;
        for (j = 0; j < f->num; j++)
        {
            FLINT_ASSERT(fmpz_is_one(f->exp + j));

            fmpz_mpoly_to_univar(u, f->poly + j, v, ctx);
            FLINT_ASSERT(u->length > 0);

            success = fmpz_mpoly_univar_content_mpoly(c, u, ctx);
            if (!success)
                goto cleanup;

            fmpz_mpoly_univar_divexact_mpoly(u, c, ctx);

            fmpz_mpoly_factor_mul_mpoly_ui(newf, c, 1, ctx);

            fmpz_add(var_powers + v, var_powers + v, u->exps + u->length - 1);
            _fmpz_mpoly_univar_shift_right(u, u->exps + u->length - 1, ctx);

            if (u->length > 1)
            {
                fmpz_mpoly_from_univar_bits(c, A->bits, u, v, ctx);
                fmpz_mpoly_factor_append_ui(newf, c, 1, ctx);
            }
            else
            {
                FLINT_ASSERT(fmpz_mpoly_is_one(u->coeffs + 0, ctx));
            }
        }

        fmpz_mpoly_factor_swap(f, newf, ctx);
    }

    for (v = 0; v < ctx->minfo->nvars; v++)
    {
        if (fmpz_is_zero(var_powers + v))
            continue;
        fmpz_mpoly_gen(c, v, ctx);
        fmpz_mpoly_factor_append_fmpz(f, c, var_powers + v, ctx);
    }

    /* 2. ensure factors are squarefree */
    FLINT_ASSERT(fmpz_mpoly_factor_matches(A, f, ctx));

    fmpz_set(newf->constant, f->constant);
    newf->num = 0;
    for (j = 0; j < f->num; j++)
    {
        success = _squarefree_factors(tempf, f->poly + j, ctx);
        if (!success)
            goto cleanup;
        fmpz_mpoly_factor_mul_factor_fmpz(newf, tempf, f->exp + j, ctx);
    }
    fmpz_mpoly_factor_swap(f, newf, ctx);

    success = 1;

cleanup:

    fmpz_mpoly_factor_clear(newf, ctx);
    fmpz_mpoly_factor_clear(tempf, ctx);
    fmpz_mpoly_univar_clear(u, ctx);
    fmpz_mpoly_clear(c, ctx);
    _fmpz_vec_clear(var_powers, ctx->minfo->nvars);

    FLINT_ASSERT(fmpz_mpoly_factor_matches(A, f, ctx));

    return success;
}


int fmpz_mpoly_factor(
    fmpz_mpoly_factor_t f,
    const fmpz_mpoly_t A,
    const fmpz_mpoly_ctx_t ctx)
{
    int success;
    slong j, v;
    fmpz_mpoly_t c;
    fmpz_mpoly_univar_t u;
    fmpz_mpoly_factor_t newf, tempf;
    fmpz * var_powers;

    /* 0. set trivial factorization */
    fmpz_mpoly_factor_one(f, ctx);
    if (fmpz_mpoly_is_fmpz(A, ctx))
    {
        fmpz_mpoly_get_fmpz(f->constant, A, ctx);
        return 1;
    }
    else
    {
        FLINT_ASSERT(A->length > 0);
        _fmpz_vec_content(f->constant, A->coeffs, A->length);
        if (fmpz_sgn(A->coeffs + 0) < 0)
            fmpz_neg(f->constant, f->constant);

        fmpz_mpoly_factor_fit_length(f, 1, ctx);
        fmpz_mpoly_scalar_divexact_fmpz(f->poly + 0, A, f->constant, ctx);
        fmpz_one(f->exp + 0);
        f->num = 1;
    }

    if (A->bits > FLINT_BITS)
    {
        return 0;
    }

    fmpz_mpoly_factor_init(newf, ctx);
    fmpz_mpoly_factor_init(tempf, ctx);
    fmpz_mpoly_univar_init(u, ctx);
    fmpz_mpoly_init(c, ctx);
    var_powers = _fmpz_vec_init(ctx->minfo->nvars);

    /* 1. ensure factors are primitive w.r.t any variable */
    FLINT_ASSERT(fmpz_mpoly_factor_matches(A, f, ctx));

    for (v = 0; v < ctx->minfo->nvars; v++)
    {
        fmpz_set(newf->constant, f->constant);
        newf->num = 0;
        for (j = 0; j < f->num; j++)
        {
            FLINT_ASSERT(fmpz_is_one(f->exp + j));

            fmpz_mpoly_to_univar(u, f->poly + j, v, ctx);
            FLINT_ASSERT(u->length > 0);

            success = fmpz_mpoly_univar_content_mpoly(c, u, ctx);
            if (!success)
                goto cleanup;

            fmpz_mpoly_univar_divexact_mpoly(u, c, ctx);

            fmpz_mpoly_factor_mul_mpoly_ui(newf, c, 1, ctx);

            fmpz_add(var_powers + v, var_powers + v, u->exps + u->length - 1);
            _fmpz_mpoly_univar_shift_right(u, u->exps + u->length - 1, ctx);

            if (u->length > 1)
            {
                fmpz_mpoly_from_univar_bits(c, A->bits, u, v, ctx);
                fmpz_mpoly_factor_append_ui(newf, c, 1, ctx);
            }
            else
            {
                FLINT_ASSERT(fmpz_mpoly_is_one(u->coeffs + 0, ctx));
            }
        }

        fmpz_mpoly_factor_swap(f, newf, ctx);
    }

    for (v = 0; v < ctx->minfo->nvars; v++)
    {
        if (fmpz_is_zero(var_powers + v))
            continue;
        fmpz_mpoly_gen(c, v, ctx);
        fmpz_mpoly_factor_append_fmpz(f, c, var_powers + v, ctx);
    }

    /* 2. ensure factors are squarefree */
    FLINT_ASSERT(fmpz_mpoly_factor_matches(A, f, ctx));

    fmpz_set(newf->constant, f->constant);
    newf->num = 0;
    for (j = 0; j < f->num; j++)
    {
        success = _squarefree_factors(tempf, f->poly + j, ctx);
        if (!success)
            goto cleanup;
        fmpz_mpoly_factor_mul_factor_fmpz(newf, tempf, f->exp + j, ctx);
    }
    fmpz_mpoly_factor_swap(f, newf, ctx);

    /* 3. ensure factors are irreducible */
    FLINT_ASSERT(fmpz_mpoly_factor_matches(A, f, ctx));

    fmpz_set(newf->constant, f->constant);
    newf->num = 0;
    for (j = 0; j < f->num; j++)
    {
        success = _irreducible_factors(tempf, f->poly + j, ctx);
        if (!success)
            goto cleanup;
        fmpz_mpoly_factor_mul_factor_fmpz(newf, tempf, f->exp + j, ctx);
    }
    fmpz_mpoly_factor_swap(f, newf, ctx);

    success = 1;

cleanup:

    fmpz_mpoly_factor_clear(newf, ctx);
    fmpz_mpoly_factor_clear(tempf, ctx);
    fmpz_mpoly_univar_clear(u, ctx);
    fmpz_mpoly_clear(c, ctx);
    _fmpz_vec_clear(var_powers, ctx->minfo->nvars);

    FLINT_ASSERT(fmpz_mpoly_factor_matches(A, f, ctx));

    return success;
}

