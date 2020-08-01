/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fq_nmod_mpoly_factor.h"


#if WANT_ASSERT
static void fq_nmod_polyu_get_fq_nmod_polyun(
    fq_nmod_polyu_t A,
    const fq_nmod_polyun_t B,
    const fq_nmod_ctx_t ctx)
{
    slong i, j;
    A->length = 0;
    for (i = 0; i < B->length; i++)
    {
        for (j = B->terms[i].coeff->length - 1; j >= 0; j--)
        {
            if (fq_nmod_is_zero(B->terms[i].coeff->coeffs + j, ctx))
                continue;
            fq_nmod_polyu_fit_length(A, A->length + 1, ctx);
            fq_nmod_set(A->coeffs + A->length, B->terms[i].coeff->coeffs + j, ctx);
            A->exps[A->length] = B->terms[i].exp + j;
            A->length++;
        }
    }
}

static void fq_nmod_polyu_sort_terms(fq_nmod_polyu_t A, const fq_nmod_ctx_t ctx)
{
    slong i, j;
    for (i = 1; i < A->length; i++)
    for (j = i; j > 0 && A->exps[j - 1] < A->exps[j]; j--)
    {
        ulong t = A->exps[j];
        A->exps[j] = A->exps[j - 1];
        A->exps[j - 1] = t;
        fq_nmod_swap(A->coeffs + j - 1, A->coeffs + j, ctx);
    }
    return;
}

static void fq_nmod_polyu_combine_like_terms(
    fq_nmod_polyu_t A,
    const fq_nmod_ctx_t ctx)
{
    slong in, out;

    out = -1;

    for (in = 0; in < A->length; in++)
    {
        FLINT_ASSERT(in > out);

        if (out >= 0 && A->exps[out] == A->exps[in])
        {
            fq_nmod_add(A->coeffs + out, A->coeffs + out, A->coeffs + in, ctx);
        }
        else
        {
            if (out < 0 || !fq_nmod_is_zero(A->coeffs + out, ctx))
                out++;

            if (out != in)
            {
                A->exps[out] = A->exps[in];
                fq_nmod_set(A->coeffs + out, A->coeffs + in, ctx);
            }
        }
    }

    if (out < 0 || !fq_nmod_is_zero(A->coeffs + out, ctx))
        out++;

    A->length = out;
}

static int fq_nmod_polyu_equal(
    fq_nmod_polyu_t A,
    fq_nmod_polyu_t B,
    const fq_nmod_ctx_t ctx)
{
    slong i;
    if (A->length != B->length)
        return 0;
    for (i = 0; i < A->length; i++)
    {
        if (!fq_nmod_equal(A->coeffs + i, B->coeffs + i, ctx))
            return 0;
        if (A->exps[i] != B->exps[i])
            return 0;
    }
    return 1;
}


static void fq_nmod_polyu_mul(
    fq_nmod_polyu_t A,
    const fq_nmod_polyu_t B,
    const fq_nmod_polyu_t C,
    const fq_nmod_ctx_t ctx)
{
    slong Ai, Bi, Ci;

    FLINT_ASSERT(A != B);
    FLINT_ASSERT(A != C);

    fq_nmod_polyu_fit_length(A, B->length*C->length, ctx);
    Ai = 0;
    for (Bi = 0; Bi < B->length; Bi++)
    for (Ci = 0; Ci < C->length; Ci++)
    {
        A->exps[Ai] = B->exps[Bi] + C->exps[Ci];
        fq_nmod_mul(A->coeffs + Ai, B->coeffs + Bi, C->coeffs + Ci, ctx);
        Ai++;
    }
    A->length = Ai;
    fq_nmod_polyu_sort_terms(A, ctx);
    fq_nmod_polyu_combine_like_terms(A, ctx);
}

#endif


void fq_nmod_polyu3_interp_reduce_bpoly(
    fq_nmod_bpoly_t Ap,
    const fq_nmod_polyu_t A,
    const fq_nmod_t alpha,
    const fq_nmod_ctx_t ctx)
{
    slong i;
    slong cur0, cur1, e0, e1, e2;
    fq_nmod_t p, t;
/*
    const fq_nmod_struct * Acoeffs = A->coeffs;
    const ulong * Aexps = A->exps;
flint_printf("n_polyu3_mod_interp_reduce_2sm_bpoly called alpha = %wd\n", alpha->coeffs[1]);
flint_printf("A: "); n_polyu3_print_pretty(A, "Y", "X", "Z"); flint_printf("\n");
*/
    fq_nmod_init(p, ctx);
    fq_nmod_init(t, ctx);

    fq_nmod_bpoly_zero(Ap, ctx);

    FLINT_ASSERT(A->length > 0);

    i = 0;

    cur0 = extract_exp(A->exps[i], 2, 3);
    cur1 = extract_exp(A->exps[i], 1, 3);
    e2   = extract_exp(A->exps[i], 0, 3);

    fq_nmod_pow_ui(t, alpha, e2, ctx);

    for (i = 1; i < A->length; i++)
    {
        e0 = extract_exp(A->exps[i], 2, 3);
        e1 = extract_exp(A->exps[i], 1, 3);
        e2 = extract_exp(A->exps[i], 0, 3);

        FLINT_ASSERT(e0 <= cur0);
        if (e0 < cur0 || e1 < cur1)
        {
            fq_nmod_bpoly_set_coeff(Ap, cur0, cur1, t, ctx);
            fq_nmod_zero(t, ctx);
        }
        else
        {
            FLINT_ASSERT(e0 == cur0);
            FLINT_ASSERT(e1 == cur1);
        }

        cur0 = e0;
        cur1 = e1;

        fq_nmod_pow_ui(p, alpha, e2, ctx);
        fq_nmod_add(t, t, p, ctx);
    }

    fq_nmod_bpoly_set_coeff(Ap, cur0, cur1, t, ctx);

    fq_nmod_clear(p, ctx);
    fq_nmod_clear(t, ctx);

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

void fq_nmod_polyu3n_interp_lift_sm_bpoly(
    slong * lastdeg,
    fq_nmod_polyun_t T,
    const fq_nmod_bpoly_t A,
    const fq_nmod_t alpha,
    const fq_nmod_ctx_t ctx)
{
    slong lastlength = 0;
    slong Ti;
    slong Ai, j;
/*
flint_printf("n_polyu3n_mod_interp_lift_2sm_bpoly called\n");
flint_printf("A: "); n_bpoly_print_pretty(A, "Y", "X"); flint_printf("\n");
flint_printf("B: "); n_bpoly_print_pretty(B, "Y", "X"); flint_printf("\n");
*/
    Ti = 0;

    for (Ai = A->length - 1; Ai >= 0; Ai--)
    {
        fq_nmod_poly_struct * Ac = A->coeffs + Ai;
        for (j = Ac->length - 1; j >= 0; j--)
        {
            if (fq_nmod_is_zero(Ac->coeffs + j, ctx))
                continue;
            fq_nmod_polyun_fit_length(T, Ti + 1, ctx);
            T->terms[Ti].exp = pack_exp3(Ai, j, 0);
            fq_nmod_poly_set_fq_nmod(T->terms[Ti].coeff, Ac->coeffs + j, ctx);
            lastlength = 1;
        }
    }

    T->length = Ti;

    *lastdeg = lastlength - 1;
/*
flint_printf("n_polyu3n_mod_interp_lift_2sm_bpoly returning lastdeg = %wd\n", lastlength - 1);
flint_printf("T: "); n_polyu3n_print_pretty(T, "Y", "X", "?", "Z"); flint_printf("\n");
*/

    return;
}

/*
    F is in Fq[x2][x0, x1]
    A is in Fq[x0, x1]
*/
int fq_nmod_polyu3n_interp_crt_sm_bpoly(
    slong * lastdeg,
    fq_nmod_polyun_t F,
    fq_nmod_polyun_t T,
    fq_nmod_bpoly_t A,
    const fq_nmod_poly_t modulus,
    const fq_nmod_t alpha,
    const fq_nmod_ctx_t ctx)
{
    int changed = 0;
    slong lastlength = 0;
    fq_nmod_polyun_term_struct * Tterms;
    slong Ti;
    fq_nmod_polyun_term_struct * Fterms = F->terms;
    slong Flen = F->length;
    slong Fi;
    fq_nmod_poly_struct * Acoeffs = A->coeffs;
    slong Ai, ai;
    fq_nmod_t v;
    fq_nmod_poly_t tp;

    fq_nmod_init(v, ctx);
    fq_nmod_poly_init(tp, ctx);
/*
flint_printf("n_polyu3n_mod_interp_crt_2sm_bpoly called alpha = %wu\n", alpha);
flint_printf("+: "); n_bpoly_print_pretty(A, "Y", "X"); flint_printf("\n");
flint_printf("-: "); n_bpoly_print_pretty(B, "Y", "X"); flint_printf("\n");
*/

    FLINT_ASSERT(fq_nmod_bpoly_is_canonical(A, ctx));
    FLINT_ASSERT(fq_nmod_polyun_is_canonical(F, ctx));

    fq_nmod_polyun_fit_length(T, FLINT_MAX(Flen, A->length), ctx);
    Tterms = T->terms;

    Ti = Fi = 0;
    Ai = A->length - 1;
    ai = (Ai < 0) ? 0 : fq_nmod_poly_degree(A->coeffs + Ai, ctx);

    while (Fi < Flen || Ai >= 0)
    {
        if (Ti >= T->alloc)
        {
            slong extra = Flen - Fi;
            extra = FLINT_MAX(extra, Ai);
            fq_nmod_polyun_fit_length(T, Ti + extra + 1, ctx);
            Tterms = T->terms;
        }

        FLINT_ASSERT(Fi >= Flen || Fterms[Fi].coeff->length > 0);
        FLINT_ASSERT(Ai < 0 || !fq_nmod_is_zero(Acoeffs[Ai].coeffs + ai, ctx));

        if (Fi < Flen && Ai >= 0 && Fterms[Fi].exp == pack_exp3(Ai, ai, 0))
        {
            /* F term ok, A term ok */
            fq_nmod_poly_evaluate_fq_nmod(v, Fterms[Fi].coeff, alpha, ctx);
            fq_nmod_sub(v, Acoeffs[Ai].coeffs + ai, v, ctx);
            if (!fq_nmod_is_zero(v, ctx))
            {
                changed = 1;
                fq_nmod_poly_scalar_mul_fq_nmod(tp, modulus, v, ctx);
                fq_nmod_poly_add(Tterms[Ti].coeff, Fterms[Fi].coeff, tp, ctx);
            }
            else
            {
                fq_nmod_poly_set(Tterms[Ti].coeff, Fterms[Fi].coeff, ctx);
            }
            Tterms[Ti].exp = Fterms[Fi].exp;

            Fi++;

            do {
                ai--;
            } while (ai >= 0 && fq_nmod_is_zero(Acoeffs[Ai].coeffs + ai, ctx));
            if (ai < 0)
            {
                do {
                    Ai--;
                } while (Ai >= 0 && Acoeffs[Ai].length == 0);
                if (Ai >= 0)
                    ai = fq_nmod_poly_degree(Acoeffs + Ai, ctx);
            }
        }
        else if (Fi < Flen && (Ai < 0 || Fterms[Fi].exp > pack_exp3(Ai, ai, 0)))
        {
            /* F term ok, A term missing */
            fq_nmod_poly_evaluate_fq_nmod(v, Fterms[Fi].coeff, alpha, ctx);
            if (!fq_nmod_is_zero(v, ctx))
            {
                changed = 1;
                fq_nmod_poly_scalar_mul_fq_nmod(tp, modulus, v, ctx);
                fq_nmod_poly_sub(Tterms[Ti].coeff, Fterms[Fi].coeff, tp, ctx);
            }
            else
            {
                fq_nmod_poly_set(Tterms[Ti].coeff, Fterms[Fi].coeff, ctx);
            }

            Tterms[Ti].exp = Fterms[Fi].exp;

            Fi++;
        }
        else
        {
            FLINT_ASSERT(Ai >= 0 && (Fi >= Flen || Fterms[Fi].exp < pack_exp3(Ai, ai, 0)));
            /* F term missing, Aterm ok */
            changed = 1;
            fq_nmod_poly_scalar_mul_fq_nmod(Tterms[Ti].coeff, modulus, Acoeffs[Ai].coeffs + ai, ctx);
            Tterms[Ti].exp = pack_exp3(Ai, ai, 0);

            do {
                ai--;
            } while (ai >= 0 && fq_nmod_is_zero(Acoeffs[Ai].coeffs + ai, ctx));
            if (ai < 0)
            {
                do {
                    Ai--;
                } while (Ai >= 0 && Acoeffs[Ai].length == 0);
                if (Ai >= 0)
                    ai = fq_nmod_poly_degree(Acoeffs + Ai, ctx);
            }
        }

        FLINT_ASSERT(!fq_nmod_poly_is_zero(Tterms[Ti].coeff, ctx));
        lastlength = FLINT_MAX(lastlength, Tterms[Ti].coeff->length);
        Ti++;
    }
    T->length = Ti;

    fq_nmod_clear(v, ctx);
    fq_nmod_poly_clear(tp, ctx);

    if (changed)
        fq_nmod_polyun_swap(T, F);

    FLINT_ASSERT(fq_nmod_polyun_is_canonical(F, ctx));

    *lastdeg = lastlength - 1;
    return changed;
}



/* r factor version */
int fq_nmod_polyu3_hlift(
    slong r,
    fq_nmod_polyun_struct * BB,
    fq_nmod_polyu_t A,
    fq_nmod_polyu_struct * B,
    const fq_nmod_t beta,
    slong degree_inner, /* required degree in x */
    const fq_nmod_ctx_t ctx)
{
    int success;
    slong i, j;
    fq_nmod_polyun_t T;
    fq_nmod_bpoly_struct * Bp;
    fq_nmod_bpoly_t Ap;
    fq_nmod_poly_t modulus, tempmod, t1, t2;
    fq_nmod_t alpha, c;
    slong * BBdegZ;
    slong AdegY, AdegX, AdegZ;
    slong bad_primes_left;


flint_printf("+++++++++++++++++++++++\n");
flint_printf("fq_nmod_polyu3_hlift called: degree_inner = %wd\n", degree_inner);
fq_nmod_ctx_print(ctx);

flint_printf("beta: "); fq_nmod_print_pretty(beta, ctx); flint_printf("\n");


flint_printf("A: "); fq_nmod_polyu3_print_pretty(A, "Y", "X", "Z", ctx); printf("\n");
for (i = 0; i < r; i++)
{
flint_printf("B[%wd]: ", i); fq_nmod_polyu3_print_pretty(B + i, "Y", "X", "Z", ctx); flint_printf("\n");
}


/*
    if (r < 3)
        return fq_nmod_polyu3_mod_hlift2(BB + 0, BB + 1, A, B + 0, B + 1,
                                                      beta, degree_inner, ctx);
*/
    fq_nmod_init(alpha, ctx);
    fq_nmod_init(c, ctx);

    FLINT_ASSERT(fq_nmod_polyu_is_canonical(A, ctx));
    for (i = 0; i < r; i++)
        FLINT_ASSERT(fq_nmod_polyu_is_canonical(B + i, ctx));

    BBdegZ = (slong *) flint_malloc(r*sizeof(slong));
    Bp = (fq_nmod_bpoly_struct *) flint_malloc(r*sizeof(fq_nmod_bpoly_struct));
    for (i = 0; i < r; i++)
        fq_nmod_bpoly_init(Bp + i, ctx);

    fq_nmod_polyun_init(T, ctx);
    fq_nmod_bpoly_init(Ap, ctx);
    fq_nmod_poly_init(modulus, ctx);
    fq_nmod_poly_init(tempmod, ctx);
    fq_nmod_poly_init(t1, ctx);
    fq_nmod_poly_init(t2, ctx);

    fq_nmod_polyu3_degrees(&AdegY, &AdegX, &AdegZ, A);
    if (AdegX != degree_inner)
    {
        success = -1;
        goto cleanup;
    }

    fq_nmod_poly_one(modulus, ctx);
    fq_nmod_poly_gen(tempmod, ctx);
    fq_nmod_poly_neg(tempmod, tempmod, ctx);

    fq_nmod_zero(alpha, ctx);

    bad_primes_left = FLINT_MAX(5, AdegZ);

choose_prime:

    if (fq_nmod_next(alpha, ctx) == 0)
    {
        success = -1;
        goto cleanup;
    }

flint_printf("------ alpha: "); fq_nmod_print_pretty(alpha, ctx); flint_printf(" -------\n");


    fq_nmod_polyu3_interp_reduce_bpoly(Ap, A, alpha, ctx);
/*flint_printf(" Ap: "); n_bpoly_print_pretty(Ap, "y", "x"); flint_printf("\n");*/

    for (i = 0; i < r; i++)
    {
        fq_nmod_polyu3_interp_reduce_bpoly(Bp + i, B + i, alpha, ctx);
/*
flint_printf("Bp[%wd]: ", i);
n_bpoly_print_pretty(Bp + i, "Y", "X");
flint_printf("\n");
fflush(stdout);
*/
    }

flint_printf("calling fq_nmod_bpoly_hlift r = %wd\n", r);
    if (r < 3)
        success = fq_nmod_bpoly_hlift2(Ap, Bp + 0, Bp + 1, beta, degree_inner, ctx);
    else
        success = fq_nmod_bpoly_hlift(r, Ap, Bp, beta, degree_inner, ctx);

flint_printf("returned from fq_nmod_bpoly_hlift success = %d\n", success);

    if (success < 1)
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
    if (fq_nmod_poly_degree(modulus, ctx) > 0)
    {
        fq_nmod_poly_evaluate_fq_nmod(c, modulus, alpha, ctx);
        fq_nmod_inv(c, c, ctx);
        fq_nmod_poly_scalar_mul_fq_nmod(modulus, modulus, c, ctx);
        for (i = 0; i < r; i++)
        {
            fq_nmod_polyu3n_interp_crt_sm_bpoly(BBdegZ + i, BB + i, T,
                                             Bp + i, modulus, alpha, ctx);
        }
    }
    else
    {
        for (i = 0; i < r; i++)
        {
            fq_nmod_polyu3n_interp_lift_sm_bpoly(BBdegZ + i, BB + i,
                                                       Bp + i, alpha, ctx);
        }
    }

    fq_nmod_poly_set_coeff(tempmod, 0, alpha, ctx);
    fq_nmod_poly_mul(modulus, modulus, tempmod, ctx);
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

    if (fq_nmod_poly_degree(modulus, ctx) <= AdegZ)
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
        fq_nmod_polyu_t T1, T2, T3;
        fq_nmod_polyu_init(T1, ctx);
        fq_nmod_polyu_init(T2, ctx);
        fq_nmod_polyu_init(T3, ctx);
        fq_nmod_polyu_get_fq_nmod_polyun(T2, BB + 0, ctx);
        fq_nmod_polyu_get_fq_nmod_polyun(T3, BB + 1, ctx);
        fq_nmod_polyu_mul(T1, T2, T3, ctx);
        for (i = 2; i < r; i++)
        {
            fq_nmod_polyu_get_fq_nmod_polyun(T3, BB + i, ctx);
            fq_nmod_polyu_mul(T2, T1, T3, ctx);
            fq_nmod_polyu_swap(T2, T1);
        }
        FLINT_ASSERT(fq_nmod_polyu_equal(A, T1, ctx));
        fq_nmod_polyu_clear(T1, ctx);
        fq_nmod_polyu_clear(T2, ctx);
        fq_nmod_polyu_clear(T3, ctx);
    }
#endif

    fq_nmod_polyun_clear(T, ctx);
    fq_nmod_bpoly_clear(Ap, ctx);

    for (i = 0; i < r; i++)
    {
        fq_nmod_bpoly_clear(Bp + i, ctx);
    }

    flint_free(BBdegZ);
    flint_free(Bp);

    fq_nmod_poly_clear(tempmod, ctx);
    fq_nmod_poly_clear(modulus, ctx);
    fq_nmod_poly_clear(t1, ctx);
    fq_nmod_poly_clear(t2, ctx);

    fq_nmod_clear(alpha, ctx);
    fq_nmod_clear(c, ctx);

FLINT_ASSERT(0);

    return success;
}
