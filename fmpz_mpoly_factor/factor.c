/*
    Copyright (C) 2019 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fmpz_mpoly_factor.h"
#include "fmpq_poly.h"
#include "fmpz_mod_mpoly.h"
#define LOW_THIRD_MASK ((-UWORD(1)) >> (FLINT_BITS - FLINT_BITS/3))

void usleep(ulong);

/******* nmod_bpoly - (Z/nZ)[x,y] ****************************************/

typedef struct
{
    nmod_poly_struct * coeffs;
    slong alloc;
    slong length;
    nmod_t modulus;
} nmod_bpoly_struct;

typedef nmod_bpoly_struct nmod_bpoly_t[1];

void nmod_bpoly_swap(nmod_bpoly_t A, nmod_bpoly_t B);

void nmod_bpoly_init(nmod_bpoly_t A, const nmod_t modulus);

void nmod_bpoly_clear(nmod_bpoly_t A);

void nmod_bpoly_fit_length(nmod_bpoly_t A, slong len);

void nmod_bpoly_set_coeff(nmod_bpoly_t A, slong xi, slong yi, mp_limb_t c);

void nmod_bpoly_print(nmod_bpoly_t A, const char * xvar, const char * yvar);

void nmod_mpoly_to_bpoly(
    nmod_bpoly_t A,
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

void nmod_poly_to_ptvals(
    nmod_poly_t v,
    slong d,
    const nmod_poly_t a,
    mp_limb_t * M,
    nmod_t ctx)
{
    nmod_poly_fit_length(v, 2*d + 1);
    v->length = 2*d + 1;
    coeffs2ptvals(v->coeffs, d, a->coeffs, a->length, M, ctx);
}

void nmod_ptvals_to_poly(
    nmod_poly_t a,
    const mp_limb_t * v,
    slong d,
    mp_limb_t * L,
    nmod_t ctx)
{
    mp_limb_t * t = flint_malloc((2*d + 1)*sizeof(mp_limb_t));
    nmod_poly_fit_length(a, 2*d + 1);
    ptvals2coeffs(a->coeffs, 2*d + 1, v, d, L, t, ctx);
    a->length = 2*d + 1;
    _nmod_poly_normalise(a);
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

void nmod_bpoly_taylor_shift_outer(
    nmod_bpoly_t A,
    mp_limb_t c)
{
    slong n, i, j;
    nmod_poly_t t;

    FLINT_ASSERT(c < A->modulus.n);

    if (c == 0)
        return;

    nmod_poly_init_mod(t, A->modulus);
    n = A->length;

    for (i = n - 2; i >= 0; i--)
    {
        for (j = i; j < n - 1; j++)
        {
            nmod_poly_scalar_mul_nmod(t, A->coeffs + j + 1, c);
            nmod_poly_add(A->coeffs + j, A->coeffs + j, t);
        }
    }

    nmod_poly_clear(t);
}

slong nmod_bpoly_degree_inner(const nmod_bpoly_t A)
{
    slong i, len = 0;
    for (i = 0; i < A->length; i++)
        len = FLINT_MAX(len, A->coeffs[i].length);
    return len - 1;
}

void nmod_bpoly_taylor_step(
    nmod_poly_t c,
    nmod_bpoly_t A,
    mp_limb_t alpha)
{
    slong i;

    FLINT_ASSERT(alpha == 0);

    if (A->length < 1)
    {
        nmod_poly_zero(c);
        return;
    }

    nmod_poly_swap(c, A->coeffs + 0);
    for (i = 0; i < A->length; i++)
        nmod_poly_swap(A->coeffs + i, A->coeffs + i + 1);

    A->length--;
}


void nmod_bpoly_to_ptvals(
    nmod_bpoly_t A,
    slong d,
    const nmod_bpoly_t B,
    mp_limb_t * M,
    nmod_t ctx)
{
    slong i;
    nmod_bpoly_fit_length(A, B->length);
    A->length = B->length;
    for (i = 0; i < B->length; i++)
        nmod_poly_to_ptvals(A->coeffs + i, d, B->coeffs + i, M, ctx);
}



int dense_bivar_disolve2(
    nmod_bpoly_t Delta1,
    nmod_bpoly_t Delta2,
    nmod_bpoly_t A,
    nmod_bpoly_t B1,
    nmod_bpoly_t B2,
    nmod_t ctx,
    slong degy,
    slong Delta1degy,
    slong Delta2degy)
{
    slong d = 10, l, i, j;
    mp_limb_t alpha = 0;
    mp_limb_t * M, * L, * v;
    nmod_poly_t g, q, e, b2i, b1i, delta;
    nmod_bpoly_t beta1ptvals, beta2ptvals;
    nmod_bpoly_t delta1ptvals, delta2ptvals;

printf("dense_bivar_disolve2 called\n");

    printf("A: "); nmod_bpoly_print(A, "y", "x"); printf("\n");
    printf("B1: "); nmod_bpoly_print(B1, "y", "x"); printf("\n");
    printf("B2: "); nmod_bpoly_print(B2, "y", "x"); printf("\n");

    nmod_poly_init_mod(delta, ctx);
    nmod_poly_init_mod(g, ctx);
    nmod_poly_init_mod(q, ctx);
    nmod_poly_init_mod(e, ctx);
    nmod_poly_init_mod(b1i, ctx);
    nmod_poly_init_mod(b2i, ctx);
    nmod_bpoly_init(beta1ptvals, ctx);
    nmod_bpoly_init(beta2ptvals, ctx);
    nmod_bpoly_init(delta1ptvals, ctx);
    nmod_bpoly_init(delta2ptvals, ctx);
    v = (mp_limb_t *) flint_malloc((2*d + 1)*sizeof(mp_limb_t));
    M = (mp_limb_t *) flint_malloc(d*d*sizeof(mp_limb_t));
    L = (mp_limb_t *) flint_malloc((2*d + 1)*(d + 1)*sizeof(mp_limb_t));
    computeM(M, d, ctx);
    computeL(L, d, ctx);

try_alpha:

    nmod_bpoly_taylor_shift_outer(B1, alpha);
    nmod_bpoly_taylor_shift_outer(B2, alpha);

    FLINT_ASSERT(B1->length > 0);
    FLINT_ASSERT(B2->length > 0);

printf("B1[0]: "); nmod_poly_print_pretty(B1->coeffs + 0, "x"); printf("\n");
printf("B2[0]: "); nmod_poly_print_pretty(B2->coeffs + 0, "x"); printf("\n");

    nmod_poly_xgcd(g, b2i, b1i, B1->coeffs + 0, B2->coeffs + 0);
    if (!nmod_poly_is_one(g))
    {
        FLINT_ASSERT(0);
        alpha++;
        goto try_alpha;
    }

    nmod_bpoly_to_ptvals(beta2ptvals, d, B1, M, ctx);
    nmod_bpoly_to_ptvals(beta1ptvals, d, B2, M, ctx);

printf("beta1ptvals: "); nmod_bpoly_print(beta1ptvals, "y", "_"); printf("\n");
printf("beta2ptvals: "); nmod_bpoly_print(beta2ptvals, "y", "_"); printf("\n");


    nmod_bpoly_fit_length(delta1ptvals, degy + 1);
    nmod_bpoly_fit_length(delta2ptvals, degy + 1);
    delta1ptvals->length = 0;
    delta2ptvals->length = 0;

    nmod_bpoly_fit_length(Delta1, degy + 1);
    nmod_bpoly_fit_length(Delta2, degy + 1);

    for (l = 0; l <= degy; l++)
    {
flint_printf("------------ l: %wd ------------\n", l);
printf("A: "); nmod_bpoly_print(A, "y", "x"); printf("\n");

        _nmod_vec_zero(v, 2*d + 1);
        for (i = 0; i < l; i++)
        {
            j = l - i;
            if (i < delta1ptvals->length && j < beta1ptvals->length)
                _nmod_vec_addmul(v, delta1ptvals->coeffs[i].coeffs, beta1ptvals->coeffs[j].coeffs, 2*d + 1, ctx);
            if (i < delta2ptvals->length && j < beta2ptvals->length)
                _nmod_vec_addmul(v, delta2ptvals->coeffs[i].coeffs, beta2ptvals->coeffs[j].coeffs, 2*d + 1, ctx);
        }

        nmod_ptvals_to_poly(delta, v, d, L, ctx);

printf("delta: "); nmod_poly_print_pretty(delta, "x"); printf("\n");

        nmod_bpoly_taylor_step(e, A, alpha);

flint_printf("A[%wd]: ", l); nmod_poly_print_pretty(e, "x"); printf("\n");


        nmod_poly_sub(e, e, delta);

printf("solving\n");
printf("e: "); nmod_poly_print_pretty(e, "x"); printf("\n");

printf("B1[0]: "); nmod_poly_print_pretty(B1->coeffs + 0, "x"); printf("\n");
printf("B2[0]: "); nmod_poly_print_pretty(B2->coeffs + 0, "x"); printf("\n");

        nmod_poly_mul(g, e, b1i);
        nmod_poly_divrem(q, Delta1->coeffs + l, g, B1->coeffs + 0);
        nmod_poly_mul(g, e, b2i);
        nmod_poly_divrem(q, Delta2->coeffs + l, g, B2->coeffs + 0);

flint_printf("Delta1[%wd]: ", l); nmod_poly_print_pretty(Delta1->coeffs + l, "x"); printf("\n");
flint_printf("Delta2[%wd]: ", l); nmod_poly_print_pretty(Delta2->coeffs + l, "x"); printf("\n");

        nmod_poly_to_ptvals(delta1ptvals->coeffs + l, d, Delta1->coeffs + l, M, ctx);
        nmod_poly_to_ptvals(delta2ptvals->coeffs + l, d, Delta2->coeffs + l, M, ctx);


        if (!nmod_poly_is_zero(Delta1->coeffs + l))
        {
            delta1ptvals->length = Delta1->length = l + 1;
            if (l > Delta1degy)
            {
                FLINT_ASSERT(0);
                return 0;
            }
        }

        if (!nmod_poly_is_zero(Delta2->coeffs + l))
        {
            delta2ptvals->length = Delta2->length = l + 1;
            if (l > Delta2degy)
            {
                FLINT_ASSERT(0);
                return 0;
            }
        }
    }

    nmod_bpoly_taylor_shift_outer(Delta1, -alpha);
    nmod_bpoly_taylor_shift_outer(Delta2, -alpha);

    printf("Delta1: "); nmod_bpoly_print(Delta1, "y", "x"); printf("\n");
    printf("Delta2: "); nmod_bpoly_print(Delta2, "y", "x"); printf("\n");

    nmod_poly_clear(delta);
    nmod_poly_clear(g);
    nmod_poly_clear(q);
    nmod_poly_clear(e);
    nmod_poly_clear(b1i);
    nmod_poly_clear(b2i);
    nmod_bpoly_clear(beta1ptvals);
    nmod_bpoly_clear(beta2ptvals);
    nmod_bpoly_clear(delta1ptvals);
    nmod_bpoly_clear(delta2ptvals);
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
    solve A/(B1*B2) = C1/B1 + C2/B2 in Fp[y][x]

    return:
        1: solution found with deg_y(Ci) <= ldegCibound
        0: no solution exists with deg_y(Ci) <= ldegCibound
       -1: could not find enough evaluation points where the Bi are pariwise prime
       -2: found no evaluation points where the Bi are pariwise prime
*/
int nmod_mpolyn_pfrac_bivar(
    nmod_mpolyn_struct * C,
    slong * ldegCbound,
    nmod_mpolyn_t A,
    nmod_mpolyn_struct * B,
    slong k,
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

    ldegB = (slong *) TMP_ALLOC(2*k*sizeof(slong));
    ldegC = ldegB + k;

    Bevalp = (nmod_poly_struct **) TMP_ALLOC(4*k*sizeof(nmod_poly_struct *));
    Bevalm = Bevalp + k;
    Cevalp = Bevalm + k;
    Cevalm = Cevalp + k;

    ldegA = nmod_mpolyn_lastdeg(A, ctx);
    for (i = 0; i < k; i++)
        ldegB[i] = nmod_mpolyn_lastdeg(B + i, ctx);

printf("nmod_mpolyn_pfrac_bivar called\n");
flint_printf(" A(deg %wd): ", ldegA); nmod_mpolyn_print_pretty(A, vars, ctx); printf("\n");
for (i = 0; i < k; i++)
{
    flint_printf("B[%wd](deg %wd): ", i, ldegB[i]); nmod_mpolyn_print_pretty(B + i, vars, ctx); printf("\n");
}

    FLINT_ASSERT(k > 2);
    FLINT_ASSERT(Sp->ctx->ffinfo->mod.n == ctx->ffinfo->mod.n);
    FLINT_ASSERT(Sp->ctx->minfo->nvars == ctx->minfo->nvars);
    FLINT_ASSERT(Sp->bits == bits);
    FLINT_ASSERT(A->bits == bits);
    for (i = 0; i < k; i++)
    {
        FLINT_ASSERT(B[i].bits == bits);
        FLINT_ASSERT(C[i].bits == bits);
    }
    N = mpoly_words_per_exp_sp(bits, ctx->minfo);
    mpoly_gen_offset_shift_sp(&off, &shift, 0, bits, ctx->minfo);

    nmod_poly_stack_fit_request_poly(Sp, 7 + 4*k);
    Aevalp      = nmod_poly_stack_take_top_poly(Sp);
    Aevalm      = nmod_poly_stack_take_top_poly(Sp);
    modulus     = nmod_poly_stack_take_top_poly(Sp);
    modulus2    = nmod_poly_stack_take_top_poly(Sp);
    alphapow    = nmod_poly_stack_take_top_poly(Sp);
    t1          = nmod_poly_stack_take_top_poly(Sp);
    t2          = nmod_poly_stack_take_top_poly(Sp);
    for (i = 0; i < k; i++)
    {
        Bevalp[i] = nmod_poly_stack_take_top_poly(Sp);
        Bevalm[i] = nmod_poly_stack_take_top_poly(Sp);
        Cevalp[i] = nmod_poly_stack_take_top_poly(Sp);
        Cevalm[i] = nmod_poly_stack_take_top_poly(Sp);
    }

    nmod_poly_stack_fit_request_mpolyn(Sp, 1);
    T           = nmod_poly_stack_take_top_mpolyn(Sp);

    ldegBtotal = ldegB[0];
    for (i = 1; i < k; i++)
        ldegBtotal += ldegB[i];

    bound = ldegA;
    for (i = 0; i < k; i++)
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
    for (i = 0; i < k; i++)
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
    for (i = 0; i < k; i++)
    {
        nmod_poly_one(t2);
        for (j = 0; j < k; j++)
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
        for (j = 0; j < k; j++)
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
        for (i = 0; i < k; i++)
        {
            nmod_mpolyn_interp_crt_2sm_poly(ldegC + i, C + i, T, Cevalp[i],
                                           Cevalm[i], modulus, alphapow, ctx);
        }
    }
    else
    {
        for (i = 0; i < k; i++)
        {
            nmod_mpolyn_interp_lift_2sm_poly(ldegC + i, C + i, Cevalp[i],
                                                       Cevalm[i], alpha, ctx);
        }
    }
    temp = nmod_mul(alpha, alpha, ctx->ffinfo->mod);
    nmod_poly_scalar_mul_nmod(modulus2, modulus, temp);
    nmod_poly_shift_left(modulus, modulus, 2);
    nmod_poly_sub(modulus, modulus, modulus2);

    for (i = 0; i < k; i++)
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

    nmod_poly_stack_give_back_poly(Sp, 7 + 4*k);
    nmod_poly_stack_give_back_mpolyn(Sp, 1);

    TMP_END;

flint_printf("nmod_mpolyn_pfrac2_bivar returning success = %d\n", success);
for (i = 0; i < k; i++)
{
flint_printf("C[%wd]: ", i); nmod_mpolyn_print_pretty(C + i, vars, ctx); printf("\n");
}

    return success;
}

int nmod_poly_is_canonical(const nmod_poly_t A)
{
    if (A->length <= 0)
        return A->length == 0;

    return 0 != A->coeffs[A->length - 1];
}

int nmod_bpoly_is_canonical(const nmod_bpoly_t A)
{
    slong i;

    if (A->length <= 0)
        return A->length == 0;

    for (i = 0; i < A->length; i++)
    {
        if (!nmod_poly_is_canonical(A->coeffs + i))
            return 0;
    }

    return !nmod_poly_is_zero(A->coeffs + A->length - 1);
}


void nmod_bpoly_mul(nmod_bpoly_t A, const nmod_bpoly_t B, const nmod_bpoly_t C)
{
    slong i, j;
    nmod_poly_t t;

    FLINT_ASSERT(A != B);
    FLINT_ASSERT(A != C);

    if (B->length <= 0 || C->length <= 0)
    {
        A->length = 0;
        return;
    }

    nmod_poly_init_mod(t, B->modulus);

    nmod_bpoly_fit_length(A, B->length + C->length - 1);
    for (i = 0; i < B->length + C->length - 1; i++)
        nmod_poly_zero(A->coeffs + i);

    for (i = 0; i < B->length; i++)
    {
        for (j = 0; j < C->length; j++)
        {
            nmod_poly_mul(t, B->coeffs + i, C->coeffs + j);
            nmod_poly_add(A->coeffs + i + j, A->coeffs + i + j, t);
        }
    }

    A->length = B->length + C->length - 1;

    nmod_poly_clear(t);
}

int nmod_bpoly_equal(const nmod_bpoly_t A, const nmod_bpoly_t B)
{
    slong i;

    if (A->length != B->length)
        return 0;

    for (i = 0; i < A->length; i++)
    {
        if (!nmod_poly_equal(A->coeffs + i, B->coeffs + i))
            return 0;
    }

    return 1;
}

/*
    input A, B0, B1 with A(y,x) = B0(y,x) * B1(y,x) mod (y-alpha)
    return
       -1: B0(alpha,x) & B1(alpha,x) are not pairwise prime, or
           A(alpha,x) has wrong degree w.r.t x
        0: lift of B0 and B1 to true factors is impossible
        1: successfully lifted B0 and B1 to true factors without changing lc_x
*/
int nmod_bpoly_hensel_lift2(
    nmod_bpoly_t A, /* clobbered (shifted by alpha) */
    nmod_bpoly_t B0,
    nmod_bpoly_t B1,
    mp_limb_t alpha,
    slong degree_inner, /* required degree in x */
    nmod_t ctx)
{
    int success;
    slong i, j;
    nmod_poly_t c, s, t, u, v;
/*
flint_printf("nmod_bpoly_hensel_lift2 called: degree_inner = %wd, alpha = %wu\n", degree_inner, alpha);
printf(" A: "); nmod_bpoly_print(A, "y", "x"); printf("\n");
printf("B0: "); nmod_bpoly_print(B0, "y", "x"); printf("\n");
printf("B1: "); nmod_bpoly_print(B1, "y", "x"); printf("\n");
*/
    FLINT_ASSERT(B0->length > 0);
    FLINT_ASSERT(B1->length > 0);
    FLINT_ASSERT(nmod_bpoly_is_canonical(A));
    FLINT_ASSERT(nmod_bpoly_is_canonical(B0));
    FLINT_ASSERT(nmod_bpoly_is_canonical(B1));

    nmod_poly_init_mod(c, ctx);
    nmod_poly_init_mod(s, ctx);
    nmod_poly_init_mod(t, ctx);
    nmod_poly_init_mod(u, ctx);
    nmod_poly_init_mod(v, ctx);

    nmod_bpoly_taylor_shift_outer(A, alpha);
    nmod_bpoly_taylor_shift_outer(B0, alpha);
    nmod_bpoly_taylor_shift_outer(B1, alpha);

    /* supposed to have A(alpha,x) = B0(alpha,x) * B1(alpha,x) */
    FLINT_ASSERT(nmod_poly_degree(A->coeffs + 0) == nmod_poly_degree(B0->coeffs + 0) +
                                                    nmod_poly_degree(B1->coeffs + 0));

    if (nmod_poly_degree(A->coeffs + 0) != degree_inner)
    {
        success = -1;
        goto cleanup;
    }

    /* the required degree in x is supposed to be deg_x(A) */
    FLINT_ASSERT(nmod_bpoly_degree_inner(A) == nmod_poly_degree(A->coeffs + 0));

    if (!nmod_poly_invmod(s, B1->coeffs + 0, B0->coeffs + 0))
    {
        success = -2;
        goto cleanup;
    }

    nmod_bpoly_fit_length(B0, A->length);
    nmod_bpoly_fit_length(B1, A->length);
    for (j = 1; j < A->length; j++)
    {
        nmod_poly_set(c, A->coeffs + j);
        for (i = 0; i <= j; i++)
        {
            if (i < B0->length && j - i < B1->length)
            {
                nmod_poly_mul(t, B0->coeffs + i, B1->coeffs + j - i);
                nmod_poly_sub(c, c, t);
            }
        }

        if (nmod_poly_is_zero(c))
            continue;

        nmod_poly_mul(t, s, c);
        nmod_poly_rem(u, t, B0->coeffs + 0);
        nmod_poly_mul(t, u, B1->coeffs + 0);
        nmod_poly_sub(c, c, t);
        nmod_poly_div(v, c, B0->coeffs + 0);
        nmod_poly_add(B0->coeffs + j, B0->coeffs + j, u);
        nmod_poly_add(B1->coeffs + j, B1->coeffs + j, v);

        if (!nmod_poly_is_zero(B0->coeffs + j))
            B0->length = FLINT_MAX(B0->length, j + 1);
        if (!nmod_poly_is_zero(B1->coeffs + j))
            B1->length = FLINT_MAX(B1->length, j + 1);

        if (B0->length - 1 + B1->length - 1 > A->length - 1)
        {
            success = 0;
            goto cleanup;
        }
    }

    nmod_bpoly_taylor_shift_outer(B0, nmod_neg(alpha, ctx));
    nmod_bpoly_taylor_shift_outer(B1, nmod_neg(alpha, ctx));

    success = 1;

cleanup:
/*
printf("nmod_bpoly_hensel_lift2 returning %d\n", success);
printf("B0: "); nmod_bpoly_print(B0, "y", "x"); printf("\n");
printf("B1: "); nmod_bpoly_print(B1, "y", "x"); printf("\n");
*/
    if (success)
    {
        nmod_bpoly_t tp1, tp2;
        nmod_bpoly_init(tp1, ctx);
        nmod_bpoly_init(tp2, ctx);

        nmod_bpoly_taylor_shift_outer(A, nmod_neg(alpha, ctx));
        nmod_bpoly_mul(tp1, B0, B1);
        FLINT_ASSERT(nmod_bpoly_equal(tp1, A));

        nmod_bpoly_clear(tp1);
        nmod_bpoly_clear(tp2);
    }

    nmod_poly_clear(c);
    nmod_poly_clear(s);
    nmod_poly_clear(t);
    nmod_poly_clear(u);
    nmod_poly_clear(v);

    return success;
}
/* r factor version */
int nmod_bpoly_hensel_lift(
    slong r,
    nmod_bpoly_t A, /* clobbered (shifted by alpha) */
    nmod_bpoly_struct * B,
    mp_limb_t alpha,
    slong degree_inner, /* required degree in x */
    nmod_t ctx)
{
    int success;
    slong i, j, k, tdeg;
    nmod_poly_struct * s, * v;
    nmod_poly_t c, t, u;
    nmod_bpoly_struct * U;
/*
flint_printf("nmod_bpoly_hensel_lift(%wd) called: degree_inner = %wd, alpha = %wu\n", r, degree_inner, alpha);
printf(" A: "); nmod_bpoly_print(A, "y", "x"); printf("\n");
for (i = 0; i < r; i++)
{
flint_printf("B[%wd]: ", i); nmod_bpoly_print(B + i, "y", "x"); printf("\n");
}
*/
    FLINT_ASSERT(r > 2);
    FLINT_ASSERT(nmod_bpoly_is_canonical(A));

    for (i = 0; i < r; i++)
    {
        FLINT_ASSERT(B[i].length > 0);
        FLINT_ASSERT(nmod_bpoly_is_canonical(B + i));
    }

    U = (nmod_bpoly_struct *) flint_malloc(r*sizeof(nmod_bpoly_struct));
    for (i = 0; i < r; i++)
    {
        nmod_bpoly_init(U + i, ctx);
        nmod_bpoly_fit_length(U + i, A->length + 1);
        U[i].length = A->length;
        nmod_bpoly_fit_length(B + i, A->length);
    }

    s = (nmod_poly_struct *) flint_malloc(r*sizeof(nmod_poly_struct));
    v = (nmod_poly_struct *) flint_malloc(r*sizeof(nmod_poly_struct));
    for (i = 0; i < r; i++)
    {
        nmod_poly_init_mod(s + i, ctx);
        nmod_poly_init_mod(v + i, ctx);
    }

    nmod_poly_init_mod(c, ctx);
    nmod_poly_init_mod(t, ctx);
    nmod_poly_init_mod(u, ctx);

    nmod_bpoly_taylor_shift_outer(A, alpha);
    for (i = 0; i < r; i++)
        nmod_bpoly_taylor_shift_outer(B + i, alpha);

    /* supposed to have A(alpha,x) = B0(alpha,x) * B1(alpha,x) * ... */
    j = 0;
    for (i = 0; i < r; i++)
        j += nmod_poly_degree(B[i].coeffs + 0);
    FLINT_ASSERT(j == nmod_poly_degree(A->coeffs + 0));

    if (nmod_poly_degree(A->coeffs + 0) != degree_inner)
    {
        success = -1;
        goto cleanup;
    }

    /* the required degree in x is supposed to be deg_x(A) */
    FLINT_ASSERT(nmod_bpoly_degree_inner(A) == nmod_poly_degree(A->coeffs + 0));

    for (i = 0; i < r; i++)
    {
        nmod_poly_one(t);
        for (j = 0; j < r; j++)
        {
            if (j != i)
                nmod_poly_mul(t, t, B[j].coeffs + 0);
        }
        if (!nmod_poly_invmod(s + i, t, B[i].coeffs + 0))
        {
            success = -1;
            goto cleanup;
        }
    }

    k = r - 2;
    nmod_poly_mul(U[k].coeffs + 0, B[k].coeffs + 0, B[k + 1].coeffs + 0);
    for (k--; k > 0; k--)
        nmod_poly_mul(U[k].coeffs + 0, B[k].coeffs + 0, U[k + 1].coeffs + 0);

    for (j = 1; j < A->length; j++)
    {
        for (k = 0; k < r; k++)
            nmod_poly_zero(U[k].coeffs + j);

        k = r - 2;
        nmod_poly_zero(U[k].coeffs + j);
        for (i = 0; i <= j; i++)
        {
            nmod_poly_mul(t, B[k].coeffs + i, B[k + 1].coeffs + j - i);
            nmod_poly_add(U[k].coeffs + j, U[k].coeffs + j, t);
        }
        for (k--; k > 0; k--)
        {
            nmod_poly_zero(U[k].coeffs + j);
            for (i = 0; i <= j; i++)
            {
                nmod_poly_mul(t, B[k].coeffs + i, U[k + 1].coeffs + j - i);
                nmod_poly_add(U[k].coeffs + j, U[k].coeffs + j, t);
            }
        }

        nmod_poly_set(c, A->coeffs + j);
        for (i = 0; i <= j; i++)
        {
            nmod_poly_mul(t, B[0].coeffs + i, U[1].coeffs + j - i);
            nmod_poly_sub(c, c, t);
        }

        if (nmod_poly_is_zero(c))
            continue;

        tdeg = 0;
        for (i = 0; i < r; i++)
        {
            nmod_poly_mul(t, s + i, c);
            nmod_poly_rem(v + i, t, B[i].coeffs + 0);
            nmod_poly_add(B[i].coeffs + j, B[i].coeffs + j, v + i);
            if (!nmod_poly_is_zero(B[i].coeffs + j))
                B[i].length = FLINT_MAX(B[i].length, j + 1);
            tdeg += B[i].length - 1;
        }

        if (tdeg >= A->length)
        {
            success = 0;
            goto cleanup;
        }

        k = r - 2;
        nmod_poly_mul(t, B[k].coeffs + 0, v + k + 1);
        nmod_poly_mul(u, B[k + 1].coeffs + 0, v + k);
        nmod_poly_add(t, t, u);
        nmod_poly_add(U[k].coeffs + j, U[k].coeffs + j, t);
        for (k--; k > 0; k--)
        {
            nmod_poly_mul(u, B[k].coeffs + 0, t);
            nmod_poly_mul(t, U[k + 1].coeffs + 0, v + k);
            nmod_poly_add(t, t, u);
            nmod_poly_add(U[k].coeffs + j, U[k].coeffs + j, t);
        }
    }

    for (i = 0; i < r; i++)
        nmod_bpoly_taylor_shift_outer(B + i, nmod_neg(alpha, ctx));

    success = 1;

cleanup:
/*
flint_printf("nmod_bpoly_hensel_lift(%wd) returning %d\n", r, success);
for (i = 0; i < r; i++)
{
flint_printf("B[%wd]: ", i); nmod_bpoly_print(B + i, "y", "x"); printf("\n");
}
*/
    if (success)
    {
        nmod_bpoly_t tp1, tp2;
        nmod_bpoly_init(tp1, ctx);
        nmod_bpoly_init(tp2, ctx);

        nmod_bpoly_taylor_shift_outer(A, nmod_neg(alpha, ctx));
        nmod_bpoly_mul(tp1, B + 0, B + 1);
        for (i = 2; i < r; i++)
        {
            nmod_bpoly_mul(tp2, tp1, B + i);
            nmod_bpoly_swap(tp1, tp2);
        }
        FLINT_ASSERT(nmod_bpoly_equal(tp1, A));

        nmod_bpoly_clear(tp1);
        nmod_bpoly_clear(tp2);
    }

    for (i = 0; i < r; i++)
    {
        nmod_bpoly_clear(U + i);
        nmod_poly_clear(s + i);
        nmod_poly_clear(v + i);
    }
    flint_free(U);
    flint_free(s);
    flint_free(v);

    nmod_poly_clear(c);
    nmod_poly_clear(t);
    nmod_poly_clear(u);

    return success;
}




void nmod_mpolyn_cvtto_bpoly(
    nmod_bpoly_t A,
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
        nmod_bpoly_fit_length(A, x + 1);
        while (A->length <= x)
        {
            nmod_poly_zero(A->coeffs + A->length);
            A->length++;
        }
        nmod_poly_swap(A->coeffs + x, B->coeffs + i);
    }
}

void nmod_bpoly_cvtto_mpolyn(
    nmod_mpolyn_t A,
    nmod_bpoly_t B,
    const nmod_mpoly_ctx_t ctx)
{
    slong i;

    FLINT_ASSERT(ctx->minfo->nvars == 3);
    FLINT_ASSERT(A->bits == FLINT_BITS/3);

    nmod_mpolyn_fit_length(A, B->length, ctx);
    A->length = 0;
    for (i = B->length - 1; i >= 0; i--)
    {
        if (!nmod_poly_is_zero(B->coeffs + i))
        {
            nmod_poly_swap(A->coeffs + A->length, B->coeffs + i);
            A->exps[A->length] = i << (2*(FLINT_BITS/3));
            A->length++;
        }
    }
}

void nmod_tpoly_interp_reduce_2sm_bpoly(
    nmod_bpoly_t Ap,
    nmod_bpoly_t Am,
    nmod_mpoly_t A,
    nmod_poly_t alphapow,
    const nmod_mpoly_ctx_t ctx)
{
    nmod_mpolyn_t Apn, Amn, An;

    nmod_mpolyn_init(An, FLINT_BITS/3, ctx);
    nmod_mpolyn_init(Apn, FLINT_BITS/3, ctx);
    nmod_mpolyn_init(Amn, FLINT_BITS/3, ctx);

    nmod_mpoly_cvtto_mpolyn(An, A, 2, ctx);
    nmod_mpolyn_interp_reduce_2sm_mpolyn(Apn, Amn, An, 2, alphapow, ctx);
    nmod_mpolyn_cvtto_bpoly(Ap, Apn, ctx);
    nmod_mpolyn_cvtto_bpoly(Am, Amn, ctx);

    nmod_mpolyn_clear(An, ctx);
    nmod_mpolyn_clear(Apn, ctx);
    nmod_mpolyn_clear(Amn, ctx);
}


void nmod_tpolyn_interp_lift_2sm_bpoly(
    slong * lastdeg,
    nmod_mpolyn_t A,
    nmod_bpoly_t Ap,
    nmod_bpoly_t Am,    
    mp_limb_t alpha,
    const nmod_mpoly_ctx_t ctx)
{
    nmod_mpolyn_t Apn, Amn;

    nmod_mpolyn_init(Apn, FLINT_BITS/3, ctx);
    nmod_mpolyn_init(Amn, FLINT_BITS/3, ctx);

    nmod_bpoly_cvtto_mpolyn(Apn, Ap, ctx);
    nmod_bpoly_cvtto_mpolyn(Amn, Am, ctx);

    nmod_mpolyn_interp_lift_2sm_mpolyn(lastdeg, A, Apn, Amn, 2, alpha, ctx);

    nmod_mpolyn_clear(Apn, ctx);
    nmod_mpolyn_clear(Amn, ctx);
}


int nmod_tpolyn_interp_crt_2sm_bpoly(
    slong * lastdeg,
    nmod_mpolyn_t F,
    nmod_mpolyn_t T,
    nmod_bpoly_t A,
    nmod_bpoly_t B,
    nmod_poly_t modulus,
    nmod_poly_t alphapow,
    const nmod_mpoly_ctx_t ctx)
{
    int changed;
    nmod_mpolyn_t An, Bn;

    nmod_mpolyn_init(An, FLINT_BITS/3, ctx);
    nmod_mpolyn_init(Bn, FLINT_BITS/3, ctx);

    nmod_bpoly_cvtto_mpolyn(An, A, ctx);
    nmod_bpoly_cvtto_mpolyn(Bn, B, ctx);

    changed = nmod_mpolyn_interp_crt_2sm_mpolyn(lastdeg, F, T, An, Bn, 2, modulus, alphapow, ctx);

    nmod_mpolyn_clear(An, ctx);
    nmod_mpolyn_clear(Bn, ctx);

    return changed;
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
int nmod_tpoly_hensel_lift2(
    nmod_mpolyn_t BB0,
    nmod_mpolyn_t BB1,
    nmod_mpoly_t A,
    nmod_mpoly_t B0,
    nmod_mpoly_t B1,
    mp_limb_t beta,
    slong degree_inner, /* required degree in x */
    nmod_mpoly_ctx_t ctx)
{
    int success;
    nmod_mpolyn_t T;
    nmod_bpoly_t Ap, Am, B0p, B0m, B1p, B1m;
    nmod_poly_t modulus, modulus2, alphapow, t1, t2;
    mp_limb_t alpha, temp;
    slong ldegBB0, ldegBB1;
    slong degrees[3], Adegz, Adegx;
    slong bad_primes_left;
    const char * vars[] = {"y", "x", "z"};

    FLINT_ASSERT(ctx->minfo->nvars == 3);
    FLINT_ASSERT(A->bits == FLINT_BITS/3);
    FLINT_ASSERT(B0->bits == FLINT_BITS/3);
    FLINT_ASSERT(B1->bits == FLINT_BITS/3);
    FLINT_ASSERT(BB0->bits == FLINT_BITS/3);
    FLINT_ASSERT(BB1->bits == FLINT_BITS/3);
    FLINT_ASSERT(nmod_mpoly_is_canonical(A, ctx));
    FLINT_ASSERT(nmod_mpoly_is_canonical(B0, ctx));
    FLINT_ASSERT(nmod_mpoly_is_canonical(B1, ctx));

flint_printf("+++++++++++++++++++++++\n");
flint_printf("nmod_tpoly_hensel_lift2 called: degree_inner = %wd, beta = %wd\n", degree_inner, beta);
flint_printf("A: "); nmod_mpoly_print_pretty(A, vars, ctx); printf("\n");
flint_printf("B0: "); nmod_mpoly_print_pretty(B0, vars, ctx); printf("\n");
flint_printf("B1: "); nmod_mpoly_print_pretty(B1, vars, ctx); printf("\n");

    nmod_mpolyn_init(T, FLINT_BITS/3, ctx);
    nmod_bpoly_init(Ap, ctx->ffinfo->mod);
    nmod_bpoly_init(Am, ctx->ffinfo->mod);
    nmod_bpoly_init(B0p, ctx->ffinfo->mod);
    nmod_bpoly_init(B0m, ctx->ffinfo->mod);
    nmod_bpoly_init(B1p, ctx->ffinfo->mod);
    nmod_bpoly_init(B1m, ctx->ffinfo->mod);
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
/*
flint_printf("------ alpha: %wu ------\n", alpha);
*/
    alpha--;

    FLINT_ASSERT(0 < alpha && alpha <= ctx->ffinfo->mod.n/2);
    FLINT_ASSERT(alphapow->alloc >= 2);
    alphapow->length = 2;
    alphapow->coeffs[0] = 1;
    alphapow->coeffs[1] = alpha;

    nmod_tpoly_interp_reduce_2sm_bpoly(Ap, Am, A, alphapow, ctx);
    nmod_tpoly_interp_reduce_2sm_bpoly(B0p, B0m, B0, alphapow, ctx);
    nmod_tpoly_interp_reduce_2sm_bpoly(B1p, B1m, B1, alphapow, ctx);
/*
printf(" Ap: "); nmod_bpoly_print(Ap, "y", "x"); printf("\n");
printf("B0p: "); nmod_bpoly_print(B0p, "y", "x"); printf("\n");
printf("B1p: "); nmod_bpoly_print(B1p, "y", "x"); printf("\n");
*/
    success = nmod_bpoly_hensel_lift2(Ap, B0p, B1p, beta, degree_inner, ctx->ffinfo->mod);
    if (success <= 0)
    {
        if (success == 0 || --bad_primes_left < 0)
            goto cleanup;
        goto choose_prime;
    }

    success = nmod_bpoly_hensel_lift2(Am, B0m, B1m, beta, degree_inner, ctx->ffinfo->mod);
    if (success <= 0)
    {
        if (success == 0 || --bad_primes_left < 0)
            goto cleanup;
        goto choose_prime;
    }
/*
printf("lifted B0p: "); nmod_bpoly_print(B0p, "y", "x"); printf("\n");
printf("       B1p: "); nmod_bpoly_print(B1p, "y", "x"); printf("\n");
*/
    if (nmod_poly_degree(modulus) > 0)
    {
        temp = nmod_poly_evaluate_nmod(modulus, alpha);
        FLINT_ASSERT(temp == nmod_poly_evaluate_nmod(modulus, ctx->ffinfo->mod.n - alpha));
        temp = nmod_mul(temp, alpha, ctx->ffinfo->mod);
        temp = nmod_add(temp, temp, ctx->ffinfo->mod);
        temp = n_invmod(temp, ctx->ffinfo->mod.n);
        nmod_poly_scalar_mul_nmod(modulus, modulus, temp);
        nmod_tpolyn_interp_crt_2sm_bpoly(&ldegBB0, BB0, T, B0p, B0m, modulus, alphapow, ctx);
        nmod_tpolyn_interp_crt_2sm_bpoly(&ldegBB1, BB1, T, B1p, B1m, modulus, alphapow, ctx);
    }
    else
    {
        nmod_tpolyn_interp_lift_2sm_bpoly(&ldegBB0, BB0, B0p, B0m, alpha, ctx);
        nmod_tpolyn_interp_lift_2sm_bpoly(&ldegBB1, BB1, B1p, B1m, alpha, ctx);
    }
    temp = nmod_mul(alpha, alpha, ctx->ffinfo->mod);
    nmod_poly_scalar_mul_nmod(modulus2, modulus, temp);
    nmod_poly_shift_left(modulus, modulus, 2);
    nmod_poly_sub(modulus, modulus, modulus2);
/*
printf("modulus: "); nmod_poly_print_pretty(modulus, "x"); printf("\n");
printf("BB0: "); nmod_mpolyn_print_pretty(BB0, vars, ctx); printf("\n");
printf("BB1: "); nmod_mpolyn_print_pretty(BB1, vars, ctx); printf("\n");
fflush(stdout);
*/
    if (ldegBB0 + ldegBB1 > Adegz)
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
flint_printf("nmod_tpoly_hensel_lift2 returning %d\n", success);
flint_printf("BB0: "); nmod_mpolyn_print_pretty(BB0, vars, ctx); printf("\n");
flint_printf("BB1: "); nmod_mpolyn_print_pretty(BB1, vars, ctx); printf("\n");
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

        nmod_mpoly_from_mpolyn_perm_inflate(tBB0, FLINT_BITS/3, ctx, BB0, ctx, perm, shift, stride);
        nmod_mpoly_from_mpolyn_perm_inflate(tBB1, FLINT_BITS/3, ctx, BB1, ctx, perm, shift, stride);
        nmod_mpoly_mul(tP, tBB0, tBB1, ctx);
        FLINT_ASSERT(nmod_mpoly_equal(tP, A, ctx));

        nmod_mpoly_clear(tBB0, ctx);
        nmod_mpoly_clear(tBB1, ctx);
        nmod_mpoly_clear(tP, ctx);
    }

    nmod_mpolyn_clear(T, ctx);
    nmod_bpoly_clear(Ap);
    nmod_bpoly_clear(Am);
    nmod_bpoly_clear(B0p);
    nmod_bpoly_clear(B0m);
    nmod_bpoly_clear(B1p);
    nmod_bpoly_clear(B1m);
    nmod_poly_clear(modulus);
    nmod_poly_clear(modulus2);
    nmod_poly_clear(alphapow);
    nmod_poly_clear(t1);
    nmod_poly_clear(t2);

    return success;
}
/* r factor version */
int nmod_tpoly_hensel_lift(
    slong r,
    nmod_mpolyn_struct * BB,
    nmod_mpoly_t A,
    nmod_mpoly_struct * B,
    mp_limb_t beta,
    slong degree_inner, /* required degree in x */
    nmod_mpoly_ctx_t ctx)
{
    int success;
    slong i, j;
    nmod_mpolyn_t T;
    nmod_bpoly_struct * Bp, * Bm;
    nmod_bpoly_t Ap, Am;
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
    Bp = (nmod_bpoly_struct *) flint_malloc(r*sizeof(nmod_bpoly_struct));
    Bm = (nmod_bpoly_struct *) flint_malloc(r*sizeof(nmod_bpoly_struct));
    for (i = 0; i < r; i++)
    {
        nmod_bpoly_init(Bp + i, ctx->ffinfo->mod);
        nmod_bpoly_init(Bm + i, ctx->ffinfo->mod);
    }
    nmod_mpolyn_init(T, FLINT_BITS/3, ctx);
    nmod_bpoly_init(Ap, ctx->ffinfo->mod);
    nmod_bpoly_init(Am, ctx->ffinfo->mod);
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
flint_printf(" Ap: "); nmod_bpoly_print(Ap, "y", "x"); printf("\n");

    for (i = 0; i < r; i++)
    {
        nmod_tpoly_interp_reduce_2sm_bpoly(Bp + i, Bm + i, B + i, alphapow, ctx);
flint_printf("Bp[%wd]: ", i); nmod_bpoly_print(Bp + i, "y", "x"); printf("\n");
    }

    success = nmod_bpoly_hensel_lift(r, Ap, Bp, beta, degree_inner, ctx->ffinfo->mod);
    if (success <= 0)
    {
        if (success == 0 || --bad_primes_left < 0)
            goto cleanup;
        goto choose_prime;
    }

    success = nmod_bpoly_hensel_lift(r, Am, Bm, beta, degree_inner, ctx->ffinfo->mod);
    if (success <= 0)
    {
        if (success == 0 || --bad_primes_left < 0)
            goto cleanup;
        goto choose_prime;
    }

for (i = 0; i < r; i++)
{
flint_printf("lifted Bp[%wd]: ", i); nmod_bpoly_print(Bp + i, "y", "x"); printf("\n");
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
    nmod_bpoly_clear(Ap);
    nmod_bpoly_clear(Am);

    for (i = 0; i < r; i++)
    {
        nmod_bpoly_clear(Bp + i);
        nmod_bpoly_clear(Bm + i);
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
}




void test_new_stuff()
{
    nmod_mpoly_ctx_t ctx;
    nmod_mpoly_t a;
    nmod_mpoly_struct b[3];
    nmod_mpolyn_struct bb[3];
    const char * vars[] = {"y", "x", "z"};

flint_printf("starting new stuff\n");

    nmod_mpoly_ctx_init(ctx, 3, ORD_LEX, 101);
    nmod_mpoly_init(a, ctx);
    nmod_mpoly_init(b + 0, ctx);
    nmod_mpoly_init(b + 1, ctx);
    nmod_mpoly_init(b + 2, ctx);
    nmod_mpolyn_init(bb + 0, FLINT_BITS/3, ctx);
    nmod_mpolyn_init(bb + 1, FLINT_BITS/3, ctx);
    nmod_mpolyn_init(bb + 2, FLINT_BITS/3, ctx);

    nmod_mpoly_set_str_pretty(a , "(y*x + y + z)"
                                  "*(z*x^2 + 1 + y + z)"
                                  "*(y^2*x + 2 + y + z)", vars, ctx);
    nmod_mpoly_set_str_pretty(b + 0, "(y*x + 2 + z)", vars, ctx);
    nmod_mpoly_set_str_pretty(b + 1, "(z*x^2 + 1 + 2 + z)", vars, ctx);
    nmod_mpoly_set_str_pretty(b + 2, "(y^2*x + 2 + 2 + z)", vars, ctx);

    nmod_tpoly_hensel_lift(3, bb, a, b, 2, 4, ctx);

    nmod_mpoly_clear(a, ctx);
    nmod_mpoly_clear(b + 0, ctx);
    nmod_mpoly_clear(b + 1, ctx);
    nmod_mpoly_clear(b + 2, ctx);
    nmod_mpolyn_clear(bb + 0, ctx);
    nmod_mpolyn_clear(bb + 1, ctx);
    nmod_mpolyn_clear(bb + 2, ctx);

    nmod_mpoly_ctx_clear(ctx);

#if 0
    nmod_mpoly_ctx_init(ctx, 2, ORD_LEX, 101);
    nmod_mpoly_init(a, ctx);
    nmod_mpoly_init(b0, ctx);
    nmod_mpoly_init(b1, ctx);
    nmod_mpoly_init(b2, ctx);
    nmod_bpoly_init(A, ctx->ffinfo->mod);
    nmod_bpoly_init(B + 0, ctx->ffinfo->mod);
    nmod_bpoly_init(B + 1, ctx->ffinfo->mod);
    nmod_bpoly_init(B + 2, ctx->ffinfo->mod);

    nmod_mpoly_set_str_pretty(a , "(x^5 + (y-2) + 1)*((y-2 + 1)*x^3 + (y-2)^2*x + (y-2))", vars, ctx);
    nmod_mpoly_set_str_pretty(b0, "(x^5 + 1)", vars, ctx);
    nmod_mpoly_set_str_pretty(b1, "((y-2 + 1)*x^3)", vars, ctx);

    nmod_mpoly_to_bpoly(A, a, 1, 0, ctx);
    nmod_mpoly_to_bpoly(B + 0, b0, 1, 0, ctx);
    nmod_mpoly_to_bpoly(B + 1, b1, 1, 0, ctx);

    nmod_bpoly_hensel_lift2(A, B + 0, B + 1, 2, ctx->ffinfo->mod);
#endif


#if 0
    slong i;
    slong d = 3;
    nmod_t ctx;
    mp_limb_t * M , * L, * ptvals, * coeffs, * temp;

    nmod_init(&ctx, 101);
    M = (mp_limb_t *) flint_malloc(d*d*sizeof(mp_limb_t));
    L = (mp_limb_t *) flint_malloc((2*d + 1)*(d + 1)*sizeof(mp_limb_t));
    ptvals = (mp_limb_t *) flint_malloc((2*d + 1)*sizeof(mp_limb_t));
    coeffs = (mp_limb_t *) flint_malloc((2*d + 1)*sizeof(mp_limb_t));
    temp = (mp_limb_t *) flint_malloc((2*d + 1)*sizeof(mp_limb_t));

    computeM(M, d, ctx);
    computeL(L, d, ctx);

    coeffs[0] = 1;
    coeffs[1] = 2;
    coeffs[2] = 3;
    coeffs[3] = 4;
    coeffs[4] = 5;
    coeffs[5] = 6;
    coeffs[6] = 0;

    coeffs2ptvals(ptvals, d, coeffs, 2*d+1, M, ctx);
for (i = 0; i <= 2*d; i++)
flint_printf("ptvals[%wd]: %wu\n", i, ptvals[i]);

    coeffs[6] = 7;

    ptvals2coeffs(coeffs, 2*d + 1, ptvals, d, L, temp, ctx);
for (i = 0; i <= 2*d; i++)
flint_printf("coeffs[%wd]: %wu\n", i, coeffs[i]);

    {
        nmod_bpoly_t A, B1, B2, Delta1, Delta2;
        nmod_bpoly_init(A, ctx);
        nmod_bpoly_init(B1, ctx);
        nmod_bpoly_init(B2, ctx);
        nmod_bpoly_init(Delta1, ctx);
        nmod_bpoly_init(Delta2, ctx);

        nmod_bpoly_set_coeff(A, 0, 1, 1);
        nmod_bpoly_set_coeff(A, 0, 4, 1);
        nmod_bpoly_set_coeff(A, 1, 0, 1);
        nmod_bpoly_set_coeff(A, 1, 3, 1);
        nmod_bpoly_set_coeff(A, 2, 1, 1);
        nmod_bpoly_set_coeff(A, 3, 0, 1);
        nmod_bpoly_set_coeff(A, 3, 3, 1);
        nmod_bpoly_set_coeff(A, 6, 0, 1);
        nmod_bpoly_set_coeff(B1, 0, 3, 1);
        nmod_bpoly_set_coeff(B1, 3, 0, 1);
        nmod_bpoly_set_coeff(B2, 0, 3, 1);
        nmod_bpoly_set_coeff(B2, 2, 0, 1);
        nmod_bpoly_set_coeff(B2, 0, 0, 1);
        dense_bivar_disolve2(Delta1, Delta2, A, B1, B2, ctx, 6, 2, 3);
    }
#endif

#if 0
    nmod_mpoly_ctx_t ctx;
    nmod_mpoly_t a, b0, b1, b2;
    nmod_mpolyn_t A;
    nmod_mpolyn_struct B[3], C[3];
    slong Cbounds[3] = {4, 6, 1};
    nmod_poly_stack_t Sp;
    const char * vars[] = {"x", "y"};
    slong perm[] = {0, 1};
    ulong shift[] = {0, 0};
    ulong stride[] = {1, 1};
    flint_bitcnt_t wbits;

    nmod_mpoly_ctx_init(ctx, 2, ORD_LEX, 101);
    nmod_mpoly_init(a, ctx);
    nmod_mpoly_init(b0, ctx);
    nmod_mpoly_init(b1, ctx);
    nmod_mpoly_init(b2, ctx);

    nmod_mpoly_set_str_pretty(a , "x^5 + x^8 + x*y + 3*x^4*y + 2*x^7*y + x*y^3 + x^2*y^3 + 2*x^4*y^3 + 2*x^5*y^3 + y^4 + 2*x*y^4 + 2*x^3*y^4 + 2*x^4*y^4 + x^6*y^4 + y^6 + 2*x*y^6 + x^2*y^6 + 2*x^3*y^6 + x^6*y^6 + y^7 + x^3*y^7 + 2*y^9 + 2*x^3*y^9 + y^12", vars, ctx);
    nmod_mpoly_set_str_pretty(b0, "(x^3 + y^3)", vars, ctx);
    nmod_mpoly_set_str_pretty(b1, "(x^3 + y^2 + 1)", vars, ctx);
    nmod_mpoly_set_str_pretty(b2, "(x^3 + y^3 + 1)", vars, ctx);

flint_printf("starting new stuff\n");

    wbits = a->bits;
    wbits = FLINT_MAX(wbits, b1->bits);
    wbits = FLINT_MAX(wbits, b2->bits);

    nmod_poly_stack_init(Sp, wbits, ctx);
    nmod_mpolyn_init(A, wbits, ctx);
    nmod_mpolyn_init(B + 0, wbits, ctx);
    nmod_mpolyn_init(B + 1, wbits, ctx);
    nmod_mpolyn_init(B + 2, wbits, ctx);
    nmod_mpolyn_init(C + 0, wbits, ctx);
    nmod_mpolyn_init(C + 1, wbits, ctx);
    nmod_mpolyn_init(C + 2, wbits, ctx);

    nmod_mpoly_to_mpolyn_perm_deflate_threaded_pool(A, ctx, a, ctx, perm, shift, stride, NULL, 0);
    nmod_mpoly_to_mpolyn_perm_deflate_threaded_pool(B + 0, ctx, b0, ctx, perm, shift, stride, NULL, 0);
    nmod_mpoly_to_mpolyn_perm_deflate_threaded_pool(B + 1, ctx, b1, ctx, perm, shift, stride, NULL, 0);
    nmod_mpoly_to_mpolyn_perm_deflate_threaded_pool(B + 2, ctx, b2, ctx, perm, shift, stride, NULL, 0);

    nmod_mpolyn_pfrac_bivar(C, Cbounds, A, B, 3, ctx, Sp);

#endif

flint_printf("finished new stuff\n");
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
    ulong mask = (UWORD(1) << (FLINT_BITS/3)) - 1;

    for (i = 0; i < A->length; i++)
    {
        if (!first)
            flint_printf(" + ");
        first = 0;
        flint_printf("(");
        fmpz_mpoly_print_pretty(A->coeffs + i, vars, ctx);
        flint_printf(")*%s^%wu*%s^%wu*%s^%wu",
            var0, (A->exps[i] >> (2*(FLINT_BITS/3))) & mask,
            var1, (A->exps[i] >> (1*(FLINT_BITS/3))) & mask,
            var2, (A->exps[i] >> (0*(FLINT_BITS/3))) & mask);
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
    slong xexp, yexp;
    slong i;
    mp_limb_t eval;

    FLINT_ASSERT(ctx_sp->minfo->nvars == 3);
    FLINT_ASSERT(E->bits == FLINT_BITS/3);
    FLINT_ASSERT(1 == mpoly_words_per_exp_sp(E->bits, ctx_sp->minfo));

    FLINT_ASSERT(A->length == Ared->length);
    FLINT_ASSERT(A->length == Acur->length);
    FLINT_ASSERT(A->length == Ainc->length);

    E->length = 0;
    for (i = 0; i < A->length; i++)
    {
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



/*
    everything is in ZZ[Y,X,Z](x1,...,xn) where n = ctx->minfo->nvars
    A = prod_{0<=i<r} B[i] mod (Y - beta)
    try to lift to a true factorization, return:
        -1: inconclusive, don't try again
        0: lift is impossible without changing the lc_X
        1: the B[i] were successfully lifted - yay!
*/
static int lift_uuu(
    fmpz_mpolyu_struct * B,
    slong r,                /* r = number of factors */
    const fmpz_mpolyu_t A,
    const fmpz_t beta,
    const slong * degs,     /* degs[i] = deg_xi(A) for 0 <= i < n */
    slong degree_inner,     /* degree_inner = deg_X(A) */
    const fmpz_mpoly_ctx_t ctx)
{
    int changed, success, point_try_count;
    flint_bitcnt_t bits = A->bits;
    mpoly_bma_interpolate_ctx_t Ictx;
    nmod_mpoly_ctx_t ctx_sp;
    fmpz_mod_mpoly_ctx_t ctx_mp;
    fmpz_mpolyu_struct * H;
    fmpz_mpolyu_t T1, T2;
    /* multi precision workspace */
/*
    fmpz_mod_bma_mpoly_t Lambda;
    fmpz_mod_mpoly_t Aeval;
    fmpz_mod_mpoly_struct * Beval;
    fmpz_mpolycu_t Ainc, Acur, Ared;
    fmpz_mpolycu_struct * Binc, * Bcur, * Ared;
*/
    fmpz_t p, pm1, sshift_mp, last_unlucky_sshift_plus_1_mp, image_count_mp;
    fmpz * checkalpha_mp;
    /* single precision workspace */
    nmod_bma_mpoly_struct * Lambda_sp;
    nmod_mpoly_t Aeval_sp;
    nmod_mpoly_struct * Beval_sp;
    nmod_mpolyn_struct * BBeval_sp;
    nmod_mpolycu_t Ainc_sp, Acur_sp, Ared_sp;
    nmod_mpolycu_struct * Binc_sp, * Bcur_sp, * Bred_sp;
    mp_limb_t p_sp, sshift_sp, last_unlucky_sshift_plus_1_sp, image_count_sp;
    mp_limb_t * checkalpha_sp;
    /* misc */
    slong i, j;
    flint_rand_t randstate;
    fmpz_t subprod;
    int unlucky_count;
    fmpz_t Hmodulus;
    nmod_zip_mpolyu_struct * Z;
    slong zip_evals;

    slong v = ctx->minfo->nvars + 2;
    const char * vars[] = {"x", "y", "z", "w", "?", "?", "?"};
    const char * varsYXZ[] = {"y", "x", "z"};

printf("lift_uuu called\n");
printf("A: "); fmpz_mpolyuuu_print_pretty(A, vars[v], vars[0], vars[1], vars + 2, ctx); printf("\n");
    for (i = 0; i < r; i++)
    {
printf("B[%wd]: ", i); fmpz_mpolyuuu_print_pretty(B + i, vars[v], vars[0], vars[1], vars + 2, ctx); printf("\n");
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
    fmpz_init(last_unlucky_sshift_plus_1_mp);
    fmpz_init(Hmodulus);

    mpoly_bma_interpolate_ctx_init(Ictx, ctx->minfo->nvars);

    /* multiprecision workspace */
    fmpz_set_ui(p, 2);    /* modulus no care */
    fmpz_mod_mpoly_ctx_init(ctx_mp, 3, ORD_LEX, p); /* modulus no care */

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

    fmpz_mpolyu_init(T1, bits, ctx);
    fmpz_mpolyu_init(T2, bits, ctx);
    H = (fmpz_mpolyu_struct *) flint_malloc(r*sizeof(fmpz_mpolyu_struct));
    for (i = 0; i < r; i++)
        fmpz_mpolyu_init(H + i, bits, ctx);

    /* initial choices for the ksub degrees are the strict degree bounds */
    for (i = 0; i < ctx->minfo->nvars; i++)
    {
        Ictx->degbounds[i] = degs[i] + 1;
        Ictx->subdegs[i] = Ictx->degbounds[i];
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
    }

got_ksub:

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

pick_bma_prime:

    if (fmpz_cmp(p, subprod) < 0)
        fmpz_set(p, subprod);

    success = fmpz_next_smooth_prime(p, p);
    fmpz_sub_ui(pm1, p, 1);
    if (!success)
    {
        /* ran out of smooth primes - absolute falure */
        success = 0;
        goto cleanup;
    }

printf("------- bma p: "); fmpz_print(p); printf(" --------\n"); fflush(stdout);
usleep(1000000);

    if (fmpz_abs_fits_ui(p))
    {
        p_sp = fmpz_get_ui(p);
        sshift_sp = 1;

        unlucky_count = 0;
        last_unlucky_sshift_plus_1_sp = 0;

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

        switch (nmod_tpoly_hensel_lift(r, BBeval_sp, Aeval_sp, Beval_sp, fmpz_fdiv_ui(beta, p_sp), degree_inner, ctx_sp))
        {
            case -1:
                flint_printf("lifting failed\n"); fflush(stdout);
                FLINT_ASSERT(0);
                goto pick_bma_prime; /* TODO: not correct - must check bad ksub */
            case 0:
                flint_printf("lifting is impossible\n"); fflush(stdout);
                FLINT_ASSERT(0);
                success = 0;
                goto cleanup;
            case 1:
                flint_printf("lifting succeeded\n"); fflush(stdout);
        }

        for (i = 0; i < r; i++)
            nmod_bma_mpoly_add_point_tpoly(Lambda_sp + i, BBeval_sp + i, ctx_sp);

        for (i = 0; i < r; i++)
        {
            if ((Lambda_sp->pointcount & 1) != 0)
                goto next_bma_image_sp;
        }

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
#if 0
        fmpz_one(sshift_sp);

        unlucky_count = 0;
        last_unlucky_sshift_plus_1_sp = 0;

        fmpz_mod_ctx_set_modulus(ctx_mp->ffinfo, p);
        fmpz_mod_discrete_log_pohlig_hellman_precompute_prime(Ictx->dlogenv_sp, p);

        for (i = 0; i < r; i++)
        {
            fmpz_mod_bma_mpoly_reset_prime(Lambda_sp + i, ctx_mp->ffinfo);
            fmpz_mod_bma_mpoly_zero(Lambda_sp + i);
        }

        /* unfortunate nmod_poly's store their own ctx :( */
        for (i = 0; i < r; i++)
        {
            fmpz_mod_mpolyn_set_modulus(BBeval_sp + i, ctx_sp->ffinfo);
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

    next_bma_image_mp:

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

        switch (nmod_tpoly_hensel_lift(r, BBeval_sp, Aeval_sp, Beval_sp, fmpz_fdiv_ui(beta, p_sp), degree_inner, ctx_sp))
        {
            case -1:
                flint_printf("lifting failed\n"); fflush(stdout);
                FLINT_ASSERT(0);
                goto pick_bma_prime; /* TODO: not correct - must check bad ksub */
            case 0:
                flint_printf("lifting is impossible\n"); fflush(stdout);
                FLINT_ASSERT(0);
                success = 0;
                goto cleanup;
            case 1:
                flint_printf("lifting succeeded\n"); fflush(stdout);
        }

        for (i = 0; i < r; i++)
            nmod_bma_mpoly_add_point_tpoly(Lambda_sp + i, BBeval_sp + i, ctx_sp);

        for (i = 0; i < r; i++)
        {
            if ((Lambda_sp->pointcount & 1) != 0)
                goto next_bma_image_sp;
        }

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
#endif
        FLINT_ASSERT(0);
    }

    /* assume that the H[i] are correct modulo Hmodulus = p */
    fmpz_set(Hmodulus, p);

    /* find number of evals for zip interp */
    FLINT_ASSERT(H->length > 0);
    for (i = 0; i < r; i++)
    {
        zip_evals = H[i].coeffs[0].length;
        for (j = 1; j < H[i].length; j++)
            zip_evals = FLINT_MAX(zip_evals, H[i].coeffs[j].length);
    }
    zip_evals += 1; /* one extra check eval */

    for (i = 0; i < r; i++)
        nmod_zip_mpolyu_fit_poly(Z + i, H + i, zip_evals);

    p_sp = UWORD(1) << (FLINT_BITS - 2);

pick_zip_prime:
    /*
        Get a new machine prime for zippel interpolation.
        H is currently interpolated modulo Hmodulus.
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

    switch (nmod_tpoly_hensel_lift(r, BBeval_sp, Aeval_sp, Beval_sp, fmpz_fdiv_ui(beta, p_sp), degree_inner, ctx_sp))
    {
        case -1:
            flint_printf("zip lifting failed\n"); fflush(stdout);
            FLINT_ASSERT(0);
            goto pick_zip_prime;
        case 0:
            flint_printf("zip lifting is impossible\n"); fflush(stdout);
            FLINT_ASSERT(0);
            success = 0;
            goto cleanup;
        case 1:
            flint_printf("zip lifting succeeded\n"); fflush(stdout);
    }

    /* update the zippler */
    for (i = 0; i < r; i++)
    {
        success = nmod_zip_mpolyuuu_add_point(Z + i, BBeval_sp + i);
        if (!success)
        {
            /*  An image in Fp'[Y,X,Z] did not match the assumed formed in [Y,X,Z].
                Start all over. */
            goto pick_bma_prime;
        }
    }

    for (i = 0; i < r; i++)
    {
        if (Z[i].pointcount < zip_evals)
            goto next_zip_image;
    }

printf("got enough zip evals\n");

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

printf("H[%wd](%d): ", i, success);
fmpz_mpolyuuu_print_pretty(H + i, vars[v], vars[0], vars[1], vars + 2, ctx);
printf("\n");
fflush(stdout);

    }

    fmpz_mul_ui(Hmodulus, Hmodulus, ctx_sp->ffinfo->mod.n);

    if (changed)
    {
        /* TODO if the coefficients of the H[i] are getting too large? */
        goto pick_zip_prime;
    }

    /*
        check A = prod_i H[i]. Since multiplication is not implemented,
        the check is implemented as A/H[1]/.../H[r-1] == H[0]
    */
    i = 1;
    if (!fmpz_mpolyuu_divides(T1, A, H + i, 3, ctx))
        goto pick_zip_prime;
    for (i++; i < r; i++)
    {
        if (!fmpz_mpolyuu_divides(T2, T1, H + i, 3, ctx))
            goto pick_zip_prime;
        fmpz_mpolyu_swap(T1, T2, ctx);
    }
    if (!fmpz_mpolyu_equal(T1, H + 0, ctx))
        goto pick_zip_prime;

printf("yay!!!!!!!!!!!!!!!!!!!\n"); fflush(stdout);

    /* good to go */
    success = 1;
    for (i = 0; i < r; i++)
        fmpz_mpolyu_swap(B + i, H + i, ctx);

cleanup:

    return success;
}

static int play_lift(
    slong v,
    fmpz_mpoly_factor_t lfac,
    const fmpz * alpha,
    const fmpz_mpoly_t A,
    const slong * degs,
    const fmpz_mpoly_ctx_t ctx)
{
    int success;
    slong i, j;
    slong r = lfac->length;
    fmpz_mpoly_t t, t2, t3;
    fmpz_mpoly_struct * betas, * deltas;
    slong tdeg;
    const char * vars[] = {"x", "y", "z", "w", "?", "?", "?"};
    fmpz_mpolyu_t Au;
    fmpz_mpolyu_struct * Bu;
    flint_bitcnt_t bits;
    fmpz_mpoly_ctx_t ctx3, uctx;

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

    lift_uuu(Bu, r, Au, alpha + r, degs + 2, degs[0], uctx);
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

void fmpz_mpolyv_init(fmpz_mpolyv_t A, const fmpz_mpoly_ctx_t ctx)
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


void fmpz_mpolyv_swap(
    fmpz_mpolyv_t A,
    fmpz_mpolyv_t B,
    const fmpz_mpoly_ctx_t ctx)
{
   fmpz_mpolyv_struct t = *A;
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
        printf(j + 1 < n ? ", " : "\n");
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

void fmpz_mpoly_univar_shift_right(
    fmpz_mpoly_univar_t A,
    ulong shift,
    const fmpz_mpoly_ctx_t ctx)
{
    slong i;
    for (i = 0; i < A->length; i++)
    {
        fmpz_sub_ui(A->exps + i, A->exps + i, shift);
        FLINT_ASSERT(fmpz_sgn(A->exps + i) >= 0);
    }
}

/********* fmpz_mpoly_factor_t ***********************************************/

void fmpz_mpoly_factor_one(
    fmpz_mpoly_factor_t f,
    const fmpz_mpoly_ctx_t ctx)
{
    fmpz_one(f->content);
    f->length = 0;
}

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
        fmpz_mul(f->content, f->content, t);
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
        fmpz_mul(f->content, f->content, t);
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

    fmpz_pow_fmpz(t, b->content, e);
    fmpz_mul(a->content, a->content, t);

    for (i = 0; i < b->length; i++)
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
    int first;

    first = 1;
    for (i = A->length - 1; i >= 0; i--)
    {
        if (fmpz_poly_is_zero(A->coeffs + i))
            continue;

        if (!first)
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
        FLINT_ASSERT(fac->length > 0);
        fmpz_mpoly_mul(fac->poly + fac->length - 1,
                       fac->poly + fac->length - 1, goodtry, ctx);
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

static void _to_polyq(fmpq_poly_t A, const fmpz_mpoly_t B, const fmpz_mpoly_ctx_t ctx)
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

/*
static void _from_poly(fmpz_mpoly_t A, flint_bitcnt_t Abits, const fmpz_poly_t B, const fmpz_mpoly_ctx_t ctx)
{
    slong i;
    ulong * shift, * stride;

    if (B->length == 0)
    {
        fmpz_mpoly_zero(A, ctx);
        return;
    }

    shift = (ulong *) flint_malloc(ctx->minfo->nvars*sizeof(ulong));
    stride = (ulong *) flint_malloc(ctx->minfo->nvars*sizeof(ulong));
    for (i = 0; i < ctx->minfo->nvars; i++)
    {
        shift[i] = 0;
        stride[i] = 1;
    }

    _fmpz_mpoly_from_fmpz_poly_inflate(A, Abits, B, 0, shift, stride, ctx);

    flint_free(shift);
    flint_free(stride);
}
*/

static int _from_polyq(fmpz_mpoly_t A, flint_bitcnt_t Abits, const fmpq_poly_t B, const fmpz_mpoly_ctx_t ctx)
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
    slong m,
    const fmpz_t alpha,
    const fmpz_mpoly_ctx_t ctx)
{
    fmpz_mpoly_t g, Q, T;

    fmpz_mpoly_init(g, ctx);
    fmpz_mpoly_init(Q, ctx);
    fmpz_mpoly_init(T, ctx);

    fmpz_mpoly_gen(g, m, ctx);
    fmpz_mpoly_sub_fmpz(g, g, alpha, ctx);

    fmpz_mpolyv_fit_length(A, 8, ctx);
    fmpz_mpoly_divrem(Q, A->coeffs + 0, B, g, ctx);
    A->length = 1;

    while (!fmpz_mpoly_is_zero(Q, ctx))
    {
        fmpz_mpolyv_fit_length(A, A->length + 1, ctx);
        fmpz_mpoly_divrem(T, A->coeffs + A->length, Q, g, ctx);
        fmpz_mpoly_swap(Q, T, ctx);
        A->length++;
    }

    while (A->length > 0 && fmpz_mpoly_is_zero(A->coeffs + A->length - 1, ctx))
        A->length--;

    fmpz_mpoly_clear(g, ctx);
    fmpz_mpoly_clear(Q, ctx);
    fmpz_mpoly_clear(T, ctx);
}

void fmpz_mpoly_from_mpolyv(
    fmpz_mpoly_t A,
    const fmpz_mpolyv_t B,
    slong m,
    const fmpz_t alpha,
    const fmpz_mpoly_ctx_t ctx)
{
    slong i;
    fmpz_mpoly_t g, T;

    fmpz_mpoly_init(g, ctx);
    fmpz_mpoly_init(T, ctx);

    fmpz_mpoly_gen(g, m, ctx);
    fmpz_mpoly_sub_fmpz(g, g, alpha, ctx);

    fmpz_mpoly_zero(A, ctx);
    for (i = B->length - 1; i >= 0; i--)
    {
        fmpz_mpoly_mul(T, A, g, ctx);
        fmpz_mpoly_add(A, T, B->coeffs + i, ctx);
    }

    fmpz_mpoly_clear(g, ctx);
    fmpz_mpoly_clear(T, ctx);
}


typedef struct {
    slong n;
    slong r;
    slong l;
    fmpq_poly_struct * inv_prod_dbetas;
    fmpq_poly_struct * dbetas;
    fmpz_mpoly_struct * prod_mbetas;
    fmpz_mpolyv_struct * prod_mbetas_coeffs;
    fmpz_mpoly_struct * mbetas;
    fmpz_mpoly_struct * deltas;

    fmpz_mpoly_struct * g;
    fmpz_mpoly_struct * q;
    fmpz_mpoly_struct * qt;
    fmpz_mpoly_struct * newt;
    fmpz_mpolyv_struct * delta_coeffs;

    fmpq_poly_t dtq, S, R;

} mfactor_disolve_struct;

typedef mfactor_disolve_struct mfactor_disolve_t[1];


static void _mfactor_disolve_init(
    mfactor_disolve_t I,
    slong l, slong r,
    const fmpz_mpoly_struct * betas,
    const fmpz * alpha,
    const fmpz_mpoly_ctx_t ctx)
{
    slong i, j, k;
    fmpz_poly_t p;
    fmpq_poly_t G, S, pq;
/*
flint_printf("_mfactor_disolve_init called(l = %wd, r = %wd)\n",l,r);
*/
    I->l = l;
    I->r = r;

    FLINT_ASSERT(l > 0);

    I->inv_prod_dbetas = (fmpq_poly_struct *) flint_malloc(l*sizeof(fmpq_poly_struct));
    I->dbetas = (fmpq_poly_struct *) flint_malloc(l*sizeof(fmpq_poly_struct));
    I->prod_mbetas = (fmpz_mpoly_struct *) flint_malloc((r + 1)*l*sizeof(fmpz_mpoly_struct));
    I->prod_mbetas_coeffs = (fmpz_mpolyv_struct *) flint_malloc((r + 1)*l*sizeof(fmpz_mpolyv_struct));
    I->mbetas = (fmpz_mpoly_struct *) flint_malloc((r + 1)*l*sizeof(fmpz_mpoly_struct));
    I->deltas = (fmpz_mpoly_struct *) flint_malloc((r + 1)*l*sizeof(fmpz_mpoly_struct));

    fmpq_poly_init(I->dtq);
    fmpq_poly_init(I->S);
    fmpq_poly_init(I->R);

    I->g = (fmpz_mpoly_struct *) flint_malloc((r + 1)*sizeof(fmpz_mpoly_struct));
    I->q = (fmpz_mpoly_struct *) flint_malloc((r + 1)*sizeof(fmpz_mpoly_struct));
    I->qt = (fmpz_mpoly_struct *) flint_malloc((r + 1)*sizeof(fmpz_mpoly_struct));
    I->newt = (fmpz_mpoly_struct *) flint_malloc((r + 1)*sizeof(fmpz_mpoly_struct));
    I->delta_coeffs = (fmpz_mpolyv_struct *) flint_malloc((r + 1)*l*sizeof(fmpz_mpolyv_struct));
    for (i = 0; i <= r; i++)
    {
        fmpz_mpoly_init(I->g + i, ctx);
        fmpz_mpoly_init(I->q + i, ctx);
        fmpz_mpoly_init(I->qt + i, ctx);
        fmpz_mpoly_init(I->newt + i, ctx);
        for (j = 0; j < l; j++)
            fmpz_mpolyv_init(I->delta_coeffs + i*I->l + j, ctx);
    }

    fmpz_poly_init(p);
    fmpq_poly_init(G);
    fmpq_poly_init(S);
    fmpq_poly_init(pq);

    /* initialize deltas */
    for (i = r; i >= 0; i--)
        for (j = 0; j < l; j++)
            fmpz_mpoly_init(I->deltas + i*l + j, ctx);

    /* set betas */
    i = r;
    for (j = 0; j < l; j++)
    {
        fmpz_mpoly_init(I->mbetas + i*l + j, ctx);
        fmpz_mpoly_set(I->mbetas + i*l + j, betas + j, ctx);
/*
flint_printf("mbetas[%wd][%wd]: ",i,j); fmpz_mpoly_print_pretty(I->mbetas + i*l + j, NULL, ctx); printf("\n");
*/
    }
    for (i--; i >= 0; i--)
    {
        for (j = 0; j < l; j++)
        {
            fmpz_mpoly_init(I->mbetas + i*l + j, ctx);
            fmpz_mpoly_evaluate_one_fmpz(I->mbetas + i*l + j, I->mbetas + (i + 1)*l + j, i + 1, alpha + i, ctx);
/*
flint_printf("mbetas[%wd][%wd]: ",i,j); fmpz_mpoly_print_pretty(I->mbetas + i*l + j, NULL, ctx); printf("\n");
*/
        }
    }
    for (j = 0; j < l; j++)
    {
        _to_poly(p, I->mbetas + 0*l + j, ctx);
        fmpq_poly_init(I->dbetas + j);
        fmpq_poly_set_fmpz_poly(I->dbetas + j, p);
/*
flint_printf("dbetas[%wd]: ",j); fmpq_poly_print_pretty(I->dbetas + j, "x1"); printf("\n");
*/
    }

    /* set product of betas */
    for (i = r; i >= 0; i--)
    {
        for (j = 0; j < l; j++)
        {
            fmpz_mpoly_init(I->prod_mbetas + i*l + j, ctx);
            fmpz_mpoly_one(I->prod_mbetas + i*l + j, ctx);
            for (k = 0; k < l; k++)
            {
                if (k == j)
                    continue;
                fmpz_mpoly_mul(I->prod_mbetas + i*l + j, I->prod_mbetas + i*l + j, I->mbetas + i*l + k, ctx);
            }
            fmpz_mpolyv_init(I->prod_mbetas_coeffs + i*l + j, ctx);
            if (i > 0)
            {
                fmpz_mpoly_to_mpolyv(I->prod_mbetas_coeffs + i*l + j,
                                     I->prod_mbetas + i*l + j, i, alpha + i - 1, ctx);
            }
/*
flint_printf("prod_mbetas[%wd][%wd]: ",i,j); fmpz_mpoly_print_pretty(I->prod_mbetas + i*l + j, NULL, ctx); printf("\n");
fmpz_mpolyv_print_pretty(I->prod_mbetas_coeffs + i*l + j, NULL, ctx);
*/
        }        
    }
    for (j = 0; j < l; j++)
    {
        fmpq_poly_one(pq);
        for (k = 0; k < l; k++)
        {
            if (k == j)
                continue;
            fmpq_poly_mul(pq, pq, I->dbetas + k);
        }
        fmpq_poly_init(I->inv_prod_dbetas + j);
        fmpq_poly_xgcd(G, S, I->inv_prod_dbetas + j, I->dbetas + j, pq);
        FLINT_ASSERT(fmpq_poly_is_one(G));
/*
flint_printf("inv_prod_dbetas[%wd]: ",j); fmpq_poly_print_pretty(I->inv_prod_dbetas + j, "x1"); printf("\n");
*/
    }

    fmpz_poly_clear(p);
    fmpq_poly_clear(G);
    fmpq_poly_clear(S);
    fmpq_poly_clear(pq);
}


static void _mfactor_disolve_clear(mfactor_disolve_t I, const fmpz_mpoly_ctx_t ctx)
{
    slong i, j;

    fmpq_poly_clear(I->dtq);
    fmpq_poly_clear(I->S);
    fmpq_poly_clear(I->R);

    for (i = 0; i <= I->r; i++)
    {
        fmpz_mpoly_clear(I->g + i, ctx);
        fmpz_mpoly_clear(I->q + i, ctx);
        fmpz_mpoly_clear(I->qt + i, ctx);
        fmpz_mpoly_clear(I->newt + i, ctx);
        for (j = 0; j < I->l; j++)
            fmpz_mpolyv_clear(I->delta_coeffs + i*I->l + j, ctx);
    }
    flint_free(I->g);
    flint_free(I->q);
    flint_free(I->qt);
    flint_free(I->newt);
    flint_free(I->delta_coeffs);

    for (j = 0; j < I->l; j++)
    {
        fmpq_poly_clear(I->inv_prod_dbetas + j);
        fmpq_poly_clear(I->dbetas + j);
        for (i = 0; i <= I->r; i++)
        {
            fmpz_mpolyv_clear(I->prod_mbetas_coeffs + i*I->l + j, ctx);
            fmpz_mpoly_clear(I->prod_mbetas + i*I->l + j, ctx);
            fmpz_mpoly_clear(I->mbetas + i*I->l + j, ctx);
            fmpz_mpoly_clear(I->deltas + i*I->l + j, ctx);
        }
    }

    flint_free(I->inv_prod_dbetas);
    flint_free(I->dbetas);
    flint_free(I->prod_mbetas);
    flint_free(I->prod_mbetas_coeffs);
    flint_free(I->mbetas);
    flint_free(I->deltas);
}


static int _mfactor_disolve(
    flint_bitcnt_t bits,
    slong r,
    slong num,
    fmpz_mpoly_t t,
    const fmpz * alpha,
    const slong * deg,
    mfactor_disolve_t I,
    const fmpz_mpoly_ctx_t ctx)
{
    slong i, j, s, l;
    int success;
    fmpz_mpoly_struct * deltas = I->deltas + r*I->l;
    fmpz_mpoly_struct * newdeltas = I->deltas + (r - 1)*I->l;

    FLINT_ASSERT(num == I->l);

/*
flint_printf("_mfactor_disolve(r = %wd, num = %wd) called:\n", r, num);
flint_printf("t: "); fmpz_mpoly_print_pretty(t, NULL, ctx); printf("\n");
for (i = 0; i < num; i++)
{
flint_printf("betas[%wd]: ", i); fmpz_mpoly_print_pretty(betas + i, NULL, ctx); printf("\n");
}
*/

    if (r == 0)
    {
        _to_polyq(I->dtq, t, ctx);

        success = 1;
        for (i = 0; i < num; i++)
        {
            fmpq_poly_mul(I->S, I->dtq, I->inv_prod_dbetas + i);
            fmpq_poly_rem(I->R, I->S, I->dbetas + i);
            if (!_from_polyq(deltas + i, bits, I->R, ctx))
                return 0;
        }
        return 1;
    }
    else if (1)
    {
        fmpz_mpoly_struct * g = I->g + r;
        fmpz_mpoly_struct * q = I->q + r;
        fmpz_mpoly_struct * qt = I->qt + r;
        fmpz_mpoly_struct * newt = I->newt + r;
        fmpz_mpolyv_struct * delta_coeffs = I->delta_coeffs + r*I->l;

        for (i = 0; i < I->l; i++)
            delta_coeffs[i].length = 0;

        fmpz_mpoly_gen(g, r, ctx);
        fmpz_mpoly_sub_fmpz(g, g, alpha + r - 1, ctx);

        for (l = 0; l <= deg[r]; l++)
        {
            fmpz_mpoly_divrem(q, newt, t, g, ctx);
            fmpz_mpoly_swap(t, q, ctx);
            for (s = 0; s < l; s++)
            {
                for (i = 0; i < I->l; i++)
                {
                    if ( s < delta_coeffs[i].length
                      && l - s < I->prod_mbetas_coeffs[r*I->l + i].length)
                    {
                        fmpz_mpoly_mul(qt, I->prod_mbetas_coeffs[r*I->l + i].coeffs + l - s,
                                           delta_coeffs[i].coeffs + s, ctx);
                        fmpz_mpoly_sub(q, newt, qt, ctx);
                        fmpz_mpoly_swap(newt, q, ctx);
                    }
                }
            }

            if (!_mfactor_disolve(bits, r - 1, num, newt, alpha, deg, I, ctx))
                return 0;

            for (i = 0; i < I->l; i++)
            {
                if (!fmpz_mpoly_is_zero(newdeltas + i, ctx))
                {
                    if (l + I->prod_mbetas_coeffs[r*I->l + i].length - 1 > deg[r])
                        return 0;
                    fmpz_mpolyv_set_coeff(delta_coeffs + i, l, newdeltas + i, ctx);
                }
            }
        }

        for (i = 0; i < I->l; i++)
            fmpz_mpoly_from_mpolyv(deltas + i, delta_coeffs + i, r, alpha + r - 1, ctx);

        return 1;
    }
    else
    {
        fmpz_mpoly_t e, tt, pow, g, q, newt;

FLINT_ASSERT(0);

        fmpz_mpoly_init(e, ctx);
        fmpz_mpoly_init(tt, ctx);
        fmpz_mpoly_init(pow, ctx);
        fmpz_mpoly_init(g, ctx);
        fmpz_mpoly_init(q, ctx);
        fmpz_mpoly_init(newt, ctx);
        for (i = 0; i < num; i++)
            fmpz_mpoly_zero(deltas + i, ctx);

        fmpz_mpoly_set(e, t, ctx);

        fmpz_mpoly_one(pow, ctx);
        fmpz_mpoly_gen(g, r, ctx);
        fmpz_mpoly_sub_fmpz(g, g, alpha + r - 1, ctx);
        for (j = 0; j <= deg[r]; j++)
        {
            success = fmpz_mpoly_divides(q, e, pow, ctx);
            FLINT_ASSERT(success);
            fmpz_mpoly_evaluate_one_fmpz(newt, q, r, alpha + r - 1, ctx);
            success = _mfactor_disolve(bits, r - 1, num, newt, alpha, deg, I, ctx);
            if (!success)
                goto cleanup;

            fmpz_mpoly_set(e, t, ctx);
            for (i = 0; i < num; i++)
            {
                fmpz_mpoly_mul(tt, newdeltas + i, pow, ctx);
                fmpz_mpoly_add(deltas + i, deltas + i, tt, ctx);
                fmpz_mpoly_mul(tt, deltas + i, I->prod_mbetas + r*I->l + i, ctx);
                fmpz_mpoly_sub(e, e, tt, ctx);
            }

            fmpz_mpoly_mul(pow, pow, g, ctx);
        }

        success = fmpz_mpoly_is_zero(e, ctx);

cleanup:

        fmpz_mpoly_clear(e, ctx);
        fmpz_mpoly_clear(tt, ctx);
        fmpz_mpoly_clear(pow, ctx);
        fmpz_mpoly_clear(g, ctx);
        fmpz_mpoly_clear(q, ctx);
        fmpz_mpoly_clear(newt, ctx);
    }

/*
flint_printf("_mfactor_disolve(r = %wd, num = %wd) returning %d:\n", r, num, success);
flint_printf("t: "); fmpz_mpoly_print_pretty(t, NULL, ctx); printf("\n");
for (i = 0; i < num; i++)
{
flint_printf("deltas[%wd]: ", i); fmpz_mpoly_print_pretty(deltas + i, NULL, ctx); printf("\n");
}
*/
    return success;
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
    slong r = lfac->length;
    fmpz_mpoly_t t, t2, t3;
    fmpz_mpoly_struct * betas, * deltas;
    mfactor_disolve_t I;
    fmpz_mpolyv_t Av;
    fmpz_mpolyv_struct B[2];
    slong tdeg;
/*
timeit_t timer;
slong timetot = 0;
*/
/*
flint_printf("_mfactor_lift_quartic2 called (degs[%wd] = %wd, alpha = %wd)\n", m, degs[m], *alpha);
flint_printf("lfac: "); fmpz_mpoly_factor_print_pretty(lfac, NULL, ctx); printf("\n");
flint_printf("A: "); fmpz_mpoly_print_pretty(A, NULL, ctx); printf("\n");
*/
    FLINT_ASSERT(r == 2);

    fmpz_mpoly_init(t, ctx);
    fmpz_mpoly_init(t2, ctx);
    fmpz_mpoly_init(t3, ctx);

    fmpz_mpolyv_init(Av, ctx);
    fmpz_mpoly_to_mpolyv(Av, A, m, alpha + m - 1, ctx);
    fmpz_mpolyv_fit_length(Av, degs[m] + 1, ctx);
    for (j = Av->length; j <= degs[m]; j++)
        fmpz_mpoly_zero(Av->coeffs + j, ctx);
/*
flint_printf("Av: "); fmpz_mpolyv_print_pretty(Av, NULL, ctx); printf("\n");
*/
    for (i = 0; i < r; i++)
    {
        fmpz_mpolyv_init(B + i, ctx);
        fmpz_mpoly_to_mpolyv(B + i, lfac->poly + i, m, alpha + m - 1, ctx);
        fmpz_mpolyv_fit_length(B + i, degs[m] + 1, ctx);
        for (j = Av->length; j <= degs[m]; j++)
            fmpz_mpoly_zero(B[i].coeffs + j, ctx);
/*
flint_printf("B + %wd: ", i); fmpz_mpolyv_print_pretty(B + i, NULL, ctx); printf("\n");
*/
    }
/*
timeit_start(timer);
*/
    betas  = (fmpz_mpoly_struct * ) flint_malloc(r*sizeof(fmpz_mpoly_struct));
    for (i = 0; i < r; i++)
    {
        betas[i] = B[i].coeffs[0];
    }
/*
timeit_start(timer);
*/
    _mfactor_disolve_init(I, lfac->length, m - 1, betas, alpha, ctx);
    deltas = I->deltas + (m - 1)*I->l;
/*
timeit_stop(timer);
timetot += timer->wall;
*/

    for (j = 1; j <= degs[m]; j++)
    {
        if (j < Av->length)
            fmpz_mpoly_set(t, Av->coeffs + j, ctx);
        else
            fmpz_mpoly_zero(t, ctx);

        for (i = 0; i <= j; i++)
        {
            fmpz_mpoly_mul(t2, B[0].coeffs + i, B[1].coeffs + j - i, ctx);
            fmpz_mpoly_sub(t3, t, t2, ctx);
            fmpz_mpoly_swap(t, t3, ctx);
        }
/*
flint_printf("disolve j = %wd, t = ", j); fmpz_mpoly_print_pretty(t, NULL, ctx); printf("\n");
*/
/*
timeit_start(timer);
*/
        success = _mfactor_disolve(A->bits, m - 1, r, t, alpha, degs, I, ctx);
/*
timeit_stop(timer);
timetot += timer->wall;
*/
        if (!success)
        {
            goto cleanup;
        }

        tdeg = -r;
        for (i = 0; i < r; i++)
        {
            fmpz_mpoly_add(t3, B[i].coeffs + j, deltas + i, ctx);
            fmpz_mpoly_swap(B[i].coeffs + j, t3, ctx);
            if (!fmpz_mpoly_is_zero(B[i].coeffs + j, ctx))
                B[i].length = FLINT_MAX(B[i].length, j + 1);
            FLINT_ASSERT(B[i].length > 0);
            tdeg += B[i].length;
        }

        if (tdeg > degs[m])
        {
            success = 0;
            goto cleanup;
        }
    }

    success = 1;

cleanup:

    _mfactor_disolve_clear(I, ctx);

    flint_free(betas);

    fmpz_mpolyv_clear(Av, ctx);
    for (i = 0; i < r; i++)
    {
        if (success)
            fmpz_mpoly_from_mpolyv(lfac->poly + i, B + i, m, alpha + m - 1, ctx);
        fmpz_mpolyv_clear(B + i, ctx);
    }

    fmpz_mpoly_clear(t, ctx);
    fmpz_mpoly_clear(t2, ctx);
    fmpz_mpoly_clear(t3, ctx);
/*
flint_printf("disolve total: %wd\n", timetot);
*/
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
    slong r = lfac->length;
    fmpz_mpoly_t t, t1, t2, t3;
    fmpz_mpoly_struct * betas, * deltas;
    mfactor_disolve_t I;
    fmpz_mpolyv_t Av;
    fmpz_mpolyv_struct * B, * U;
    slong tdeg;
/*
timeit_t timer;
slong timetot = 0;
*/
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

    fmpz_mpolyv_init(Av, ctx);
    fmpz_mpoly_to_mpolyv(Av, A, m, alpha + m - 1, ctx);
    fmpz_mpolyv_fit_length(Av, degs[m] + 1, ctx);
    for (j = Av->length; j <= degs[m]; j++)
        fmpz_mpoly_zero(Av->coeffs + j, ctx);
/*
flint_printf("Av: "); fmpz_mpolyv_print_pretty(Av, NULL, ctx); printf("\n");
flint_printf("r = %wd \n", r);
*/
    for (k = 0; k < r; k++)
    {
/*
flint_printf("init loop k = %wd\n", k);
*/
        fmpz_mpolyv_init(U + k, ctx);
        fmpz_mpolyv_fit_length(U + k, degs[m] + 1, ctx);
        for (j = 0; j <= degs[m]; j++)
            fmpz_mpoly_zero(U[k].coeffs + j, ctx);

        fmpz_mpolyv_init(B + k, ctx);
        fmpz_mpoly_to_mpolyv(B + k, lfac->poly + k, m, alpha + m - 1, ctx);
        fmpz_mpolyv_fit_length(B + k, degs[m] + 1, ctx);
        for (j = Av->length; j <= degs[m]; j++)
            fmpz_mpoly_zero(B[k].coeffs + j, ctx);
/*
flint_printf("B + %wd: ", k); fmpz_mpolyv_print_pretty(B + k, NULL, ctx); printf("\n");
*/
    }
/*
timeit_start(timer);
*/

    betas  = (fmpz_mpoly_struct * ) flint_malloc(r*sizeof(fmpz_mpoly_struct));
    for (i = 0; i < r; i++)
    {
        betas[i] = B[i].coeffs[0];
/*
flint_printf("betas[%wd]: ", i); fmpz_mpoly_print_pretty(betas + i, NULL, ctx); printf("\n");
*/
    }

/*
timeit_start(timer);
*/
    _mfactor_disolve_init(I, lfac->length, m - 1, betas, alpha, ctx);
    deltas = I->deltas + (m - 1)*I->l;
/*
timeit_stop(timer);
timetot += timer->wall;
*/
    k = r - 2;
    fmpz_mpoly_mul(U[k].coeffs + 0, B[k].coeffs + 0, B[k + 1].coeffs + 0, ctx);
    for (k--; k >= 1; k--)
        fmpz_mpoly_mul(U[k].coeffs + 0, B[k].coeffs + 0, U[k + 1].coeffs + 0, ctx);

    for (j = 1; j <= degs[m]; j++)
    {
/*
flint_printf("main loop j = %wd ", j);
*/
        k = r -2;
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
/*
flint_printf("disolve j = %wd, t = ", j); fmpz_mpoly_print_pretty(t, NULL, ctx); printf("\n");
*/
        if (fmpz_mpoly_is_zero(t, ctx))
        {
            continue;
        }
/*
timeit_start(timer);
*/
        success = _mfactor_disolve(A->bits, m - 1, r, t, alpha, degs, I, ctx);
/*
timeit_stop(timer);
timetot += timer->wall;
*/
        if (!success)
        {
/*
printf("dsolve failed\n");
*/
            goto cleanup;
        }

        tdeg = -r;
        for (i = 0; i < r; i++)
        {
            fmpz_mpoly_add(t3, B[i].coeffs + j, deltas + i, ctx);
            fmpz_mpoly_swap(B[i].coeffs + j, t3, ctx);
            if (!fmpz_mpoly_is_zero(B[i].coeffs + j, ctx))
                B[i].length = FLINT_MAX(B[i].length, j + 1);
            FLINT_ASSERT(B[i].length > 0);
            tdeg += B[i].length;
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
/*
printf("cleanup  success = %d\n", success);
*/
    _mfactor_disolve_clear(I, ctx);

    flint_free(betas);

    fmpz_mpolyv_clear(Av, ctx);
    for (i = 0; i < r; i++)
    {
        if (success)
            fmpz_mpoly_from_mpolyv(lfac->poly + i, B + i, m, alpha + m - 1, ctx);
        fmpz_mpolyv_clear(B + i, ctx);
        fmpz_mpolyv_clear(U + i, ctx);
    }

    flint_free(B);
    flint_free(U);
    fmpz_mpoly_clear(t, ctx);
    fmpz_mpoly_clear(t1, ctx);
    fmpz_mpoly_clear(t2, ctx);
    fmpz_mpoly_clear(t3, ctx);
/*
flint_printf("disolve total: %wd\n", timetot);
*/
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
    slong r = lfac->length;
    fmpz_mpoly_t e, t, pow, g, q;
    fmpz_mpoly_struct * betas, * deltas;
    mfactor_disolve_t I;
/*
timeit_t timer;
slong timetot = 0;
*/
/*
flint_printf("_mfactor_lift_quintic called (m = %wd, alpha = %wd)\n", m, *alpha);
flint_printf("lfac: "); fmpz_mpoly_factor_print_pretty(lfac, NULL, ctx); printf("\n");
flint_printf("A: "); fmpz_mpoly_print_pretty(A, NULL, ctx); printf("\n");
*/
    FLINT_ASSERT(r > 1);
/*
timeit_start(timer);
*/
    fmpz_mpoly_init(e, ctx);
    fmpz_mpoly_init(t, ctx);
    fmpz_mpoly_init(pow, ctx);
    fmpz_mpoly_init(g, ctx);
    fmpz_mpoly_init(q, ctx);

    betas  = (fmpz_mpoly_struct * ) flint_malloc(r*sizeof(fmpz_mpoly_struct));
    for (i = 0; i < r; i++)
    {
        fmpz_mpoly_init(betas + i, ctx);
        fmpz_mpoly_evaluate_one_fmpz(betas + i, lfac->poly + i, m, alpha + m - 1, ctx);
/*
flint_printf("betas[%wd]: ", i); fmpz_mpoly_print_pretty(betas + i, NULL, ctx); printf("\n");
*/
    }

    fmpz_mpoly_mul(t, lfac->poly + 0, lfac->poly + 1, ctx);
    for (i = 2; i < r; i++)
        fmpz_mpoly_mul(t, t, lfac->poly + i, ctx);
    fmpz_mpoly_sub(e, A, t, ctx);

    fmpz_mpoly_one(pow, ctx);
    fmpz_mpoly_gen(g, m, ctx);
    fmpz_mpoly_sub_fmpz(g, g, alpha + m - 1, ctx);
/*
timeit_start(timer);
*/
    _mfactor_disolve_init(I, lfac->length, m - 1, betas, alpha, ctx);
    deltas = I->deltas + (m - 1)*I->l;
/*
timeit_stop(timer);
timetot += timer->wall;
*/
    for (j = 1; j <= degs[m]; j++)
    {
/*
flint_printf("<_mfactor_lift> j = %wd, error (length: %wd)\n", j, e->length); fmpz_mpoly_print_pretty(e, NULL, ctx); printf("\n");
*/
        if (fmpz_mpoly_is_zero(e, ctx))
        {
            success = 1;
            goto cleanup;
        }

        fmpz_mpoly_mul(pow, pow, g, ctx);
        success = fmpz_mpoly_divides(q, e, pow, ctx);
        FLINT_ASSERT(success);
        fmpz_mpoly_evaluate_one_fmpz(t, q, m, alpha + m - 1, ctx);
/*
flint_printf("disolve j = %wd, t = ", j); fmpz_mpoly_print_pretty(t, NULL, ctx); printf("\n");
*/
/*
timeit_start(timer);
*/
        success = _mfactor_disolve(A->bits, m - 1, r, t, alpha, degs, I, ctx);
/*
timeit_stop(timer);
timetot += timer->wall;
*/
        if (!success)
            goto cleanup;

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

    _mfactor_disolve_clear(I, ctx);

    fmpz_mpoly_clear(e, ctx);
    fmpz_mpoly_clear(t, ctx);
    fmpz_mpoly_clear(pow, ctx);
    fmpz_mpoly_clear(g, ctx);
    fmpz_mpoly_clear(q, ctx);

    for (i = 0; i < r; i++)
    {
        fmpz_mpoly_clear(betas + i, ctx);
    }

    flint_free(betas);
/*
flint_printf("disolve total: %wd\n", timetot);
*/
/*
flint_printf("_mfactor_lift_quintic returning %d\n", success);
*/

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

    FLINT_ASSERT(pfac->length > 1);

    newdeg = (slong *) flint_malloc((n + 1)*sizeof(slong));
    fmpz_mpoly_init(lcq, ctx);
    fmpz_mpoly_init(lcp, ctx);
    fmpz_mpoly_init(t, ctx);
    fmpz_mpoly_init(newq, ctx);
    fmpz_mpoly_univar_init(u, ctx);

    FLINT_ASSERT(fmpz_is_one(pfac->content));
    FLINT_ASSERT(fmpz_mpoly_factor_matches(p, pfac, ctx));

    _fmpz_mpoly_get_lc(lcq, q, ctx);
    fmpz_mpoly_evaluate_one_fmpz(lcp, lcq, m, alpha + m - 1, ctx);

    FLINT_ASSERT(lcp->length > 0);

    fmpz_mpoly_pow_ui(t, lcq, pfac->length - 1, ctx);
    fmpz_mpoly_mul(newq, q, t, ctx);
    fmpz_mpoly_degrees_si(newdeg, newq, ctx);

    fmpz_set(qfac->content, pfac->content);
    fmpz_mpoly_factor_fit_length(qfac, pfac->length, ctx);
    qfac->length = pfac->length;
    for (i = 0; i < pfac->length; i++)
    {
        _fmpz_mpoly_get_lc(t, pfac->poly + i, ctx);
        success = fmpz_mpoly_divides(t, lcp, t, ctx);
        FLINT_ASSERT(success);
        fmpz_mpoly_mul(qfac->poly + i, pfac->poly + i, t, ctx);
        _fmpz_mpoly_set_lc(qfac->poly + i, qfac->poly + i, lcq, ctx);
        fmpz_one(qfac->exp + i);
    }

    if (qfac->length == 2)
        success = _mfactor_lift_quartic2(m, qfac, alpha, newq, newdeg, ctx);
    else if (qfac->length <  20)
        success = _mfactor_lift_quartic(m, qfac, alpha, newq, newdeg, ctx);
    else
        success = _mfactor_lift_quintic(m, qfac, alpha, newq, newdeg, ctx);

    if (!success)
        goto cleanup;

    for (i = 0; i < qfac->length; i++)
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
    FLINT_ASSERT(!success || (fmpz_is_one(qfac->content) &&
                              fmpz_mpoly_factor_matches(q, qfac, ctx)));
    return success;
}


static void fmpz_poly_factor_print_pretty(const fmpz_poly_factor_t f)
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
    fmpz * lcAfaceval = _fmpz_vec_init(lcAfac->length);
    const fmpz ** salpha = (const fmpz **) flint_malloc((n + 1)*sizeof(fmpz *));
    fmpz * d = _fmpz_vec_init(1 + lcAfac->length);
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

    for (j = 0; j < lcAfac->length; j++)
        fmpz_mpoly_evaluate_all_fmpz(lcAfaceval + j, lcAfac->poly + j,
                                                 (fmpz * const *) salpha, ctx);

    fmpz_mul(d + 0, &Aufac->c, lcAfac->content);
    for (i = 0; i < lcAfac->length; i++)
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
        for (i = lcAfac->length - 1; i >= 0; i--)
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
    _fmpz_vec_clear(lcAfaceval, lcAfac->length);
    _fmpz_vec_clear(d, 1 + lcAfac->length);
    _fmpz_vec_clear(dtilde, Aufac->num);
    flint_free(salpha);

    return success;
}


static int _try_wang(
    fmpz_mpoly_factor_t fac,
    const fmpz_mpoly_t A,
    const fmpz_mpoly_ctx_t ctx)
{
    int success;
    const slong n = ctx->minfo->nvars - 1;
    slong i, j, k;
    fmpz * alpha;
    fmpz_mpoly_struct * Aevals;
    slong * deg, * degeval;
    fmpz_mpoly_factor_t tfac, lcAfac;
    fmpz_mpoly_t t, lcA, Acopy;
    fmpz_mpoly_struct * newA;
    fmpz_poly_t Au;
    fmpz_poly_factor_t Aufac;
    ulong * shift, * stride;
    slong alpha_modulus, alpha_count;
    flint_rand_t state;
    fmpz_mpoly_t m, mpow;
    fmpz_mpolyv_t new_lcs, lc_divs;
    fmpz_t q;

    FLINT_ASSERT(A->bits <= FLINT_BITS);
    FLINT_ASSERT(n > 1);

flint_printf("_try_wang called n = %wd\n", n);
flint_printf("A: "); fmpz_mpoly_print_pretty(A, NULL, ctx); printf("\n");

    fmpz_init(q);

    flint_randinit(state);

    fmpz_mpoly_init(Acopy, ctx);
    fmpz_mpoly_init(m, ctx);
    fmpz_mpoly_init(mpow, ctx);

    fmpz_mpolyv_init(new_lcs, ctx);
    fmpz_mpolyv_init(lc_divs, ctx);

    fmpz_poly_factor_init(Aufac);
    fmpz_poly_init(Au);
    shift = (ulong *) flint_malloc((n + 1)*sizeof(ulong));
    stride = (ulong *) flint_malloc((n + 1)*sizeof(ulong));
    for (i = 0; i <= n; i++)
    {
        shift[i] = 0;
        stride[i] = 1;
    }

    fmpz_mpoly_init(lcA, ctx);
    fmpz_mpoly_factor_init(lcAfac, ctx);
    alpha = _fmpz_vec_init(n);

    deg     = (slong *) flint_malloc((n + 1)*sizeof(slong));
    degeval = (slong *) flint_malloc((n + 1)*sizeof(slong));
    Aevals    = (fmpz_mpoly_struct *) flint_malloc(n*sizeof(fmpz_mpoly_struct));
    for (i = 0; i < n; i++)
        fmpz_mpoly_init(Aevals + i, ctx);

    fmpz_mpoly_factor_init(tfac, ctx);
    fmpz_mpoly_init(t, ctx);

    _fmpz_mpoly_get_lc(lcA, A, ctx);
    success = fmpz_mpoly_factor(lcAfac, lcA, 1, ctx);
    if (!success)
        goto cleanup;

    fmpz_mpoly_degrees_si(deg, A, ctx);

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
printf("---------------------------\n");
printf("alpha = "); tuple_print(alpha, n);
*/
    /* ensure degrees do not drop under evaluation */
    for (i = n - 1; i >= 0; i--)
    {
        fmpz_mpoly_evaluate_one_fmpz(Aevals + i, i == n - 1 ? A : Aevals + i + 1, i + 1, alpha + i, ctx);
        fmpz_mpoly_degrees_si(degeval, Aevals + i, ctx);
        for (j = 0; j <= i; j++)
        {
            if (degeval[j] != deg[j])
                goto next_alpha;
        }
    }

    /* make sure our univar is squarefree */
    _fmpz_mpoly_to_fmpz_poly_deflate(Au, Aevals + 0, 0, shift, stride, ctx);
    fmpz_poly_factor(Aufac, Au);
    for (j = 0; j < Aufac->num; j++)
    {
        if (Aufac->exp[j] != 1)
            goto next_alpha;
    }

    /* if the univariate is irreducible, then A irreducible */
    if (Aufac->num == 1)
    {
        fmpz_mpoly_factor_append_ui(fac, A, 1, ctx);
        success = 1;
        goto cleanup;
    }

    if (lcAfac->length > 0)
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

    FLINT_ASSERT(Aufac->num > 1);
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
    fmpz_mpoly_pow_ui(mpow, m, Aufac->num - 1, ctx);
    if (fmpz_mpoly_is_one(mpow, ctx))
    {
        newA = (fmpz_mpoly_struct *) A;
    }
    else
    {
        newA = Acopy;
        fmpz_mpoly_mul(newA, A, mpow, ctx);
    }

    fmpz_mpoly_degrees_si(deg, newA, ctx);

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

    fmpz_mpoly_factor_fit_length(fac, Aufac->num, ctx);
    fac->length = Aufac->num;
    fmpz_one(fac->content);
    for (i = 0; i < Aufac->num; i++)
    {
        FLINT_ASSERT(fmpz_mpoly_is_fmpz(new_lcs->coeffs + 0*Aufac->num + i, ctx));
        FLINT_ASSERT(fmpz_mpoly_length(new_lcs->coeffs + 0*Aufac->num + i, ctx) == 1);
        FLINT_ASSERT(fmpz_divisible(new_lcs->coeffs[i].coeffs + 0, Aufac->p[i].coeffs + Aufac->p[i].length - 1));
        fmpz_divexact(q, new_lcs->coeffs[i].coeffs + 0,
                                  Aufac->p[i].coeffs + Aufac->p[i].length - 1);
        _fmpz_mpoly_from_fmpz_poly_inflate(fac->poly + i, newA->bits,
                                          Aufac->p + i, 0, shift, stride, ctx);
        fmpz_mpoly_scalar_mul_fmpz(fac->poly + i, fac->poly + i, q, ctx);
        fmpz_one(fac->exp + i);
    }

    fmpz_mpoly_factor_fit_length(tfac, Aufac->num, ctx);
    tfac->length = Aufac->num;
    fmpz_one(tfac->content);

    for (k = 1; k <= n; k++)
    {
        fmpz_mpoly_struct * target;
        for (i = 0; i < Aufac->num; i++)
        {
            fmpz_one(tfac->exp + i);
            _fmpz_mpoly_set_lc(tfac->poly + i, fac->poly + i,
                                         new_lcs->coeffs + Aufac->num*k + i, ctx);
        }

        target = k < n ? Aevals + k : newA;

flint_printf("wang lifting variable %wd\n", k);

        if (k > 2)
            play_lift(k, tfac, alpha, target, deg, ctx);

        if (tfac->length == 2)
            success = _mfactor_lift_quartic2(k, tfac, alpha, target, deg, ctx);
        else if (tfac->length <  20)
            success = _mfactor_lift_quartic(k, tfac, alpha, target, deg, ctx);
        else
            success = _mfactor_lift_quintic(k, tfac, alpha, target, deg, ctx);

        if (!success)
            goto next_alpha;

        fmpz_mpoly_factor_swap(tfac, fac, ctx);
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

    fmpz_clear(q);

    fmpz_mpolyv_clear(new_lcs, ctx);
    fmpz_mpolyv_clear(lc_divs, ctx);

    fmpz_poly_factor_clear(Aufac);
    flint_free(shift);
    flint_free(stride);

    fmpz_mpoly_clear(lcA, ctx);
    fmpz_mpoly_factor_clear(lcAfac, ctx);
        
    _fmpz_vec_clear(alpha, n);

    for (i = 0; i < n; i++)
        fmpz_mpoly_clear(Aevals + i, ctx);
    flint_free(Aevals);

    flint_free(deg);
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

    if (_try_wang(fac, A, ctx))
        return 1;

    fmpz_init(subset);
    alphait = _fmpz_vec_init(n);
    alpha = _fmpz_vec_init(n);
    Aevals    = (fmpz_mpoly_struct *) flint_malloc(n*sizeof(fmpz_mpoly_struct));
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
        fmpz_mpoly_evaluate_one_fmpz(Aevals + i, i == n - 1 ? A : Aevals + i + 1, i + 1, alpha + i, ctx);
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

        FLINT_ASSERT(fmpz_is_one(pfac->content));
        FLINT_ASSERT(fmpz_mpoly_factor_matches(p, pfac, ctx));

        /* if pfac has only one factor, A must be irreducible */
        if (pfac->length == 1)
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
        if (pfac->length == 2)
        {
            fmpz_mpoly_factor_append_ui(fac, A, 1, ctx);
            success = 1;
            goto cleanup;
        }

        qfac->length = 0;

try_again:

        for (r = 1; r <= pfac->length/2; r++)
        {
            subset_first(subset, pfac->length, r);

            FLINT_ASSERT(fmpz_is_one(pfac->content));
            FLINT_ASSERT(fmpz_mpoly_factor_matches(p, pfac, ctx));

            do {
                fmpz_mpoly_factor_fit_length(dfac, 2, ctx);
                dfac->length = 2;
                fmpz_one(dfac->content);
                fmpz_mpoly_one(dfac->poly + 0, ctx);
                fmpz_mpoly_one(dfac->poly + 1, ctx);
                fmpz_one(dfac->exp + 0);
                fmpz_one(dfac->exp + 1);
                for (i = 0; i < pfac->length; i++)
                {
                    j = fmpz_tstbit(subset, i);
                    fmpz_mpoly_mul(dfac->poly + j, dfac->poly + j, pfac->poly + i, ctx);
                }

                success = _try_lift(tfac, q, dfac, p, m, alpha, n, ctx);
                if (success)
                {
                    for (i = pfac->length - 1; i >= 0; i--)
                    {
                        if (fmpz_tstbit(subset, i))
                        {
                            fmpz_mpoly_swap(pfac->poly + i, pfac->poly + pfac->length - 1, ctx);
                            pfac->length--;
                        }
                    }
                    fmpz_mpoly_factor_append_ui(qfac, tfac->poly + 1, 1, ctx);
                    fmpz_mpoly_swap(q, tfac->poly + 0, ctx);
                    fmpz_mpoly_swap(p, dfac->poly + 0, ctx);
                    goto try_again;
                }
            }
            while (subset_next(subset, subset, pfac->length));
        }
        /* if pfac could not be combined, p must be irreducible */
        fmpz_mpoly_factor_append_ui(qfac, q, 1, ctx);
        fmpz_mpoly_factor_swap(qfac, pfac, ctx);
    }

    success = 1;

    for (i = 0; i < pfac->length; i++)
        fmpz_mpoly_factor_append_ui(fac, pfac->poly + i, 1, ctx);

cleanup:

    fmpz_clear(subset);
    _fmpz_vec_clear(alphait, n);
    _fmpz_vec_clear(alpha, n);
    for (i = 0; i < n; i++)
    {
        fmpz_mpoly_clear(Aevals + i, ctx);
    }
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

    FLINT_ASSERT(!success || (fmpz_is_one(fac->content) &&
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
        fmpz_mul(f->content, f->content, &Aufac->c); /* Aufac->c should be 1 */
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
        if (success)
        {
            fmpz_mpoly_factor_one(f, ctx);
            for (i = 0; i < Amfactors->length; i++)
            {
                fmpz_mpoly_convert_perm(B, A->bits, ctx, Amfactors->poly + i, lctx, iperm);
                FLINT_ASSERT(B->length > 0);
                if (fmpz_sgn(B->coeffs + 0) < 0)
                    fmpz_mpoly_neg(B, B, ctx);
                fmpz_mpoly_factor_append_fmpz(f, B, Amfactors->exp + i, ctx);
            }

            fmpz_mpoly_factor_clear(Amfactors, lctx);
            fmpz_mpoly_clear(B, ctx);
            fmpz_mpoly_clear(Al, lctx);
            fmpz_mpoly_ctx_clear(lctx);
        }
    }

    TMP_END;

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

int fmpz_mpoly_factor(
    fmpz_mpoly_factor_t f,
    const fmpz_mpoly_t A,
    int full,   /* 0 squarefree, 1 irreducible */
    const fmpz_mpoly_ctx_t ctx)
{
    int success;
    slong j, v;
    fmpz_mpoly_t c;
    fmpz_mpoly_univar_t u;
    fmpz_mpoly_factor_t newf, tempf;

    /* 0. set trivial factorization */
/*
{
timeit_t timer;
timeit_start(timer);
*/
    fmpz_mpoly_factor_one(f, ctx);
    if (fmpz_mpoly_is_fmpz(A, ctx))
    {
        fmpz_mpoly_get_fmpz(f->content, A, ctx);
        return 1;
    }
    else
    {
        _fmpz_vec_content(f->content, A->coeffs, A->length);
        if (fmpz_sgn(A->coeffs + 0) < 0)
            fmpz_neg(f->content, f->content);

        fmpz_mpoly_factor_fit_length(f, 1, ctx);
        fmpz_mpoly_scalar_divexact_fmpz(f->poly + 0, A, f->content, ctx);
        fmpz_one(f->exp + 0);
        f->length = 1;
    }

    if (A->bits > FLINT_BITS)
    {
        return 0;
    }

    fmpz_mpoly_factor_init(newf, ctx);
    fmpz_mpoly_factor_init(tempf, ctx);
    fmpz_mpoly_univar_init(u, ctx);
    fmpz_mpoly_init(c, ctx);
/*
timeit_stop(timer);
flint_printf("trivial time: %wd\n", timer->wall);
}
*/
    /* 1. ensure factors are primitive w.r.t any variable */
    FLINT_ASSERT(fmpz_mpoly_factor_matches(A, f, ctx));
/*
{
timeit_t timer;
timeit_start(timer);
*/
    for (v = 0; v < ctx->minfo->nvars; v++)
    {
        fmpz_set(newf->content, f->content);
        newf->length = 0;
        for (j = 0; j < f->length; j++)
        {
            if (fmpz_is_one(f->exp + j))
            {
                fmpz_mpoly_to_univar(u, f->poly + j, v, ctx);
                FLINT_ASSERT(u->length > 0);

                success = fmpz_mpoly_univar_content_mpoly(c, u, ctx);
                if (!success)
                    goto cleanup;

                fmpz_mpoly_univar_divexact_mpoly(u, c, ctx);

                fmpz_mpoly_factor_mul_mpoly_ui(newf, c, 1, ctx);

                if (u->exps[u->length - 1] != 0)
                {
                    fmpz_mpoly_gen(c, v, ctx);
                    fmpz_mpoly_factor_append_ui(newf, c, u->exps[u->length - 1], ctx);
                    fmpz_mpoly_univar_shift_right(u, u->exps[u->length - 1], ctx);
                }

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
            else
            {
                FLINT_ASSERT(fmpz_sgn(f->exp + j) > 0);
                fmpz_mpoly_factor_append_fmpz(newf, f->poly + j, f->exp + j, ctx);
            }
        }

        fmpz_mpoly_factor_swap(f, newf, ctx);
    }
/*
timeit_stop(timer);
flint_printf("primitive time: %wd\n", timer->wall);
}
*/
    /* 2. ensure factors are squarefree */
    FLINT_ASSERT(fmpz_mpoly_factor_matches(A, f, ctx));
/*
{
timeit_t timer;
timeit_start(timer);
*/
    fmpz_set(newf->content, f->content);
    newf->length = 0;
    for (j = 0; j < f->length; j++)
    {
        success = _squarefree_factors(tempf, f->poly + j, ctx);
        if (!success)
            goto cleanup;
        fmpz_mpoly_factor_mul_factor_fmpz(newf, tempf, f->exp + j, ctx);
    }
    fmpz_mpoly_factor_swap(f, newf, ctx);

    /* skip irreducible factorization if not wanted */
    if (!full)
    {
        success = 1;
        goto cleanup;
    }
/*
timeit_stop(timer);
flint_printf("squarefree time: %wd\n", timer->wall);
}
*/
    /* 3. ensure factors are irreducible */
    FLINT_ASSERT(fmpz_mpoly_factor_matches(A, f, ctx));
/*
{
timeit_t timer;
timeit_start(timer);
*/
    fmpz_set(newf->content, f->content);
    newf->length = 0;
    for (j = 0; j < f->length; j++)
    {
        success = _irreducible_factors(tempf, f->poly + j, ctx);
        if (!success)
            goto cleanup;
        fmpz_mpoly_factor_mul_factor_fmpz(newf, tempf, f->exp + j, ctx);
    }
    fmpz_mpoly_factor_swap(f, newf, ctx);
/*
timeit_stop(timer);
flint_printf("irreducible time: %wd\n", timer->wall);
}
*/
    success = 1;

cleanup:
/*
{
timeit_t timer;
timeit_start(timer);
*/
    FLINT_ASSERT(fmpz_mpoly_factor_matches(A, f, ctx));
    fmpz_mpoly_factor_clear(newf, ctx);
    fmpz_mpoly_factor_clear(tempf, ctx);
    fmpz_mpoly_univar_clear(u, ctx);
    fmpz_mpoly_clear(c, ctx);
/*
timeit_stop(timer);
flint_printf("check time: %wd\n", timer->wall);
}
*/

    return success;
}

