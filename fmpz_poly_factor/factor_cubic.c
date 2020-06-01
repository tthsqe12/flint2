/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fmpz_poly_factor.h"

/* f(x) = x^3 - 3*a*x - b */
static void evalf(fmpz_t f, const fmpz_t x, const fmpz_t a, const fmpz_t b)
{
    fmpz_t t;
    fmpz_init(t);
    fmpz_mul(t, x, x);
    fmpz_submul_ui(t, a, 3);
    fmpz_mul(f, t, x);
    fmpz_sub(f, f, b);
    fmpz_clear(t);
}

/* f'(x) = 3*x^2 - 3*a */
static void evalfp(fmpz_t f, const fmpz_t x, const fmpz_t a, const fmpz_t b)
{
    fmpz_mul(f, x, x);
    fmpz_sub(f, f, a);
    fmpz_mul_ui(f, f, 3);
}

/*
    It is assumed that newton iteration (when performed over RR) never
    overshoots a root when starting from
            x = right (dir > 0) or
            x = left (dir < 0)
    and of course (left, right) contains a unique root.
*/
static int find_cubic_root(fmpz_t r,
    fmpz_t left, fmpz_t right,
    const fmpz_t a, const fmpz_t b,
    int dir)
{
    int ret;
    fmpz_t fleft, fright, fmid;
    fmpz_t fpleft, fpright, fpmid;
    fmpz_t t1, t2, d, mid;

    fmpz_init(t1);
    fmpz_init(t2);
    fmpz_init(d);
    fmpz_init(mid);
    fmpz_init(fleft);
    fmpz_init(fright);
    fmpz_init(fmid);
    fmpz_init(fpleft);
    fmpz_init(fpright);
    fmpz_init(fpmid);

    evalf(fleft, left, a, b);
    evalf(fright, right, a, b);

/*
printf("find_cubic_root\n");
printf("     a: "); fmpz_print(a); printf("\n");
printf("     b: "); fmpz_print(b); printf("\n");
printf("  left: "); fmpz_print(left); printf("\n");
printf(" right: "); fmpz_print(right); printf("\n");
printf(" fleft: "); fmpz_print(fleft); printf("\n");
printf("fright: "); fmpz_print(fright); printf("\n");
*/

    fmpz_mul(d, fright, fleft);
    FLINT_ASSERT(fmpz_sgn(d) < 0);

    fmpz_sub(d, right, left);
    FLINT_ASSERT(fmpz_sgn(d) > 0);

    if (fmpz_cmp_ui(d, 10) < 0)
    {
        while (fmpz_add_ui(left, left, 1), fmpz_cmp(left, right) < 0)
        {
            evalf(fleft, left, a, b);
            if (fmpz_is_zero(fleft))
            {
                fmpz_swap(r, left);
                ret = 1;
                goto cleanup;
            }
        }
        ret = 0;
        goto cleanup;
    }
    else if (dir > 0)
    {
        evalfp(fpright, right, a, b);
        fmpz_cdiv_q(mid, fright, fpright);
        FLINT_ASSERT(fmpz_sgn(mid) > 0);
        fmpz_sub(mid, right, mid);

        if (fmpz_cmp(mid, left) <= 0)
        {
            ret = 0;
            goto cleanup;
        }

        evalf(fmid, mid, a, b);

        if (fmpz_sgn(fmid) == 0)
        {
            fmpz_swap(r, mid);
            ret = 1;
            goto cleanup;
        }
        if (fmpz_sgn(fmid) == fmpz_sgn(fright))
        {
            ret = find_cubic_root(r, left, mid, a, b, dir);
            goto cleanup;
        }
        else
        {
            ret = 0;
            goto cleanup;
        }
    }
    else
    {
        evalfp(fpleft, left, a, b);
        fmpz_fdiv_q(mid, fleft, fpleft);
        FLINT_ASSERT(fmpz_sgn(mid) < 0);
        fmpz_sub(mid, left, mid);

        evalf(fmid, mid, a, b);

        if (fmpz_cmp(mid, right) >= 0)
        {
            ret = 0;
            goto cleanup;
        }

        if (fmpz_sgn(fmid) == 0)
        {
            fmpz_swap(r, mid);
            ret = 1;
            goto cleanup;
        }
        if (fmpz_sgn(fmid) == fmpz_sgn(fleft))
        {
            ret = find_cubic_root(r, mid, right, a, b, dir);
            goto cleanup;
        }
        else
        {
            ret = 0;
            goto cleanup;
        }
    }

cleanup:

    fmpz_clear(t1);
    fmpz_clear(t2);
    fmpz_clear(d);
    fmpz_clear(mid);
    fmpz_clear(fleft);
    fmpz_clear(fright);
    fmpz_clear(fmid);
    fmpz_clear(fpleft);
    fmpz_clear(fpright);
    fmpz_clear(fpmid);

    return ret;
}



/*
    write the integer roots of f(x) = x^3 - 3*a*x - b to r
    return 0: irreducible cubic
           1: (x - r[0])*(irreducible quadratic)
           2: (x - r[0])*(x - r[1])*(x - r[2])
           3: (x - r[0])*(x - r[1])^2
           4: (x - r[0])^3

    a and b are clobbered!
*/
int _fmpz_cubic_roots(fmpz * r, fmpz_t a, fmpz_t b)
{
    slong i;
    int ret, sign_a, sign_b, sign_d;
    fmpz_t d, t1, t2, t3, t4;
/*
printf("---------cubic roots called-----------\n");
printf("a: "); fmpz_print(a); printf("\n");
printf("b: "); fmpz_print(b); printf("\n");
*/
    FLINT_ASSERT(a != b);

    sign_a = fmpz_sgn(a);
    sign_b = fmpz_sgn(b);

    if (fmpz_is_pm1(b))
    {
        if (sign_a == 0)
        {
            fmpz_swap(r + 0, b);
            return 1;
        }
        else
        {
            return 0;
        }
    }

    if (fmpz_is_zero(b))
    {
        fmpz_zero(r + 0);

        if (sign_a <= 0)
            return sign_a == 0 ? 4 : 1;
        
        fmpz_mul_ui(a, a, 3);
        
        if (fmpz_is_square(a))
        {
            fmpz_sqrt(r + 1, a);
            fmpz_neg(r + 2, r + 1);
            return 2;
        }
        else
        {
            return 1;
        }
    }

    if (sign_b < 0)
        fmpz_neg(b, b);

    /* b >= 2 now, and sign_b has the sign of the original b */

    /* check irreducibility mod 2 */
    if (fmpz_is_odd(b) && fmpz_is_odd(a))
        return 0;

    /* check irreducibility mod some p */
    {
        nmod_t ctx;
        ulong A, B, F, G, H, P;
        nmod_init(&ctx, UWORD(5)*7*11*13*17*19*23*29);

        A = fmpz_fdiv_ui(a, ctx.n);
        B = fmpz_fdiv_ui(b, ctx.n);
        A = nmod_add(A, nmod_add(A, A, ctx), ctx);

        /*
            f(x) = B + A*x - x^3
            g(x) = f(x + 1) - f(x)
                 = A - 1 - 3*x - 3*x^2
            h(x) = g(x) - g(x + 1) 
                 = 6 + 6*x
        */

        F = B;
        G = nmod_sub(A, 1, ctx);
        H = 6;
        P = F;
        for (i = 1; i < 29; i++)
        {
            F = nmod_add(F, G, ctx);
            G = nmod_sub(G, H, ctx);
            H += 6;
            P = nmod_mul(P, F, ctx);
        }

        if (P != 0)
            return 0;
    }

    fmpz_init(d);
    fmpz_init(t1);
    fmpz_init(t2);
    fmpz_init(t3);
    fmpz_init(t4);

    /* d = b^2 - 4*a^3 = discriminant/-27 */
    fmpz_mul(d, b, b);
    fmpz_mul(t1, a, a);
    fmpz_mul(t2, t1, a);
    fmpz_submul_ui(d, t2, 4);

    sign_d = fmpz_sgn(d);

    if (sign_d == 0)
    {
        FLINT_ASSERT(fmpz_divisible(b, a));
        fmpz_divexact(r + 0, b, a);
        FLINT_ASSERT(fmpz_is_even(r + 0));
        fmpz_divexact_si(r + 1, r + 0, -2);
        ret = 3;
        goto cleanup;
    }

    if (sign_d > 0)
    {
        /*
            An interval containing the real root is (0, max(a,b)).
            Unfortunately newton converges too slowly starting from x = max(a,b).
            Luckily, we have a formula for the real root,

                cbrt(4b + sqrt(4*d)) + cbrt(4b - sqrt(4*d))
                -------------------------------------------,
                                     2

            which is calculated here by rounding intermediates to zero.
        */

        fmpz_sqrt(t3, d);
        fmpz_add(t4, b, t3);
        fmpz_mul_ui(t4, t4, 4);
        fmpz_root(t1, t4, 3);
        fmpz_sub(t4, b, t3);
        fmpz_mul_ui(t4, t4, 4);
        fmpz_root(t2, t4, 3);
        fmpz_add(t1, t1, t2);
        fmpz_fdiv_q_2exp(r + 0, t1, 1);

        /*
            If r is the actual real root in RR, then, experimentally,
            r[0] - 1 < r < r[0] + 2, so that r[0] - 1 cannot be a root. Test
            all three numbers r[0] - 1, r[0], r[0] + 1 anyways since these
            follow from easily-proven bounds.
        */

        fmpz_mul(t2, r + 0, r + 0);
        fmpz_sub(t3, t2, a);
        fmpz_submul_ui(t2, a, 3);
        fmpz_mul(t1, t2, r + 0);
        if (fmpz_equal(t1, b))
        {
            ret = 1;
            goto cleanup;
        }

        fmpz_submul_ui(b, r, 3);
        fmpz_mul_ui(t3, t3, 3);
        fmpz_add_ui(t3, t3, 1);
        fmpz_sub(t2, b, t3);
        if (fmpz_equal(t1, t2))
        {
            fmpz_add_ui(r + 0, r + 0, 1);
            ret = 1;
            goto cleanup;
        }

        fmpz_add(t2, b, t3);
        if (fmpz_equal(t1, t2))
        {
            fmpz_sub_ui(r + 0, r + 0, 1);
            ret = 1;
            goto cleanup;
        }

        ret = 0;
        goto cleanup;
    }

    /*
        Three real roots. The observations
                f(-2*sqrt(a)) < 0,   f(-sqrt(a)) > 0
                f(-b/(2*a)) > 0,     f(-b/(3*a)) < 0
                f(sqrt(a)) < 0,      f(2*sqrt(a)) > 0
        give intervals on which newton coverges quickly.
        However, we still do about log(a) full-precision computations, so it
        is asymptotically slower than fmpz_poly_hensel_lift_once
    */

    /* check irreducibility mod some more p */
    {
        nmod_t ctx;
        ulong A, B, F, G, H, P;
        nmod_init(&ctx, UWORD(31)*37*41*43*47);

        A = fmpz_fdiv_ui(a, ctx.n);
        B = fmpz_fdiv_ui(b, ctx.n);
        A = nmod_add(A, nmod_add(A, A, ctx), ctx);

        /*
            f(x) = B + A*x - x^3
            g(x) = f(x + 1) - f(x)
                 = A - 1 - 3*x - 3*x^2
            h(x) = g(x) - g(x + 1) 
                 = 6 + 6*x
        */

        F = B;
        G = nmod_sub(A, 1, ctx);
        H = 6;
        P = F;
        for (i = 1; i < 47; i++)
        {
            F = nmod_add(F, G, ctx);
            G = nmod_sub(G, H, ctx);
            H += 6;
            P = nmod_mul(P, F, ctx);
        }

        if (P != 0)
        {
            ret = 0;
            goto cleanup;
        }
    }

    FLINT_ASSERT(fmpz_cmp_ui(a, 2) > 0);
    fmpz_sqrt(t1, a);
    fmpz_add_ui(t2, t1, 1);
    fmpz_mul_ui(t2, t2, 2);

/*printf("finding positive real root\n");*/
    if (find_cubic_root(r + 0, t1, t2, a, b, +1))
    {
        fmpz_mul_ui(a, a, 4);
        fmpz_submul(a, r + 0, r + 0);
        fmpz_mul_ui(a, a, 3);
        if (fmpz_is_square(a))
        {
            fmpz_sqrt(t1, a);
            fmpz_sub(r + 1, t1, r + 0);
            fmpz_add(r + 2, t1, r + 0);
            FLINT_ASSERT(fmpz_is_even(r + 1));
            FLINT_ASSERT(fmpz_is_even(r + 2));
            fmpz_divexact_ui(r + 1, r + 1, 2);
            fmpz_divexact_si(r + 2, r + 2, -2);
            ret = 2;
            goto cleanup;
        }
        else
        {
            ret = 1;
            goto cleanup;
        }
    }

    /* the positive root is irrational, check the two negative ones */

    fmpz_sqrt(t1, a);
    fmpz_neg(t1, t1);

/*printf("t1: "); fmpz_print(t1); printf("\n");*/

    fmpz_sub_ui(t2, t1, 1);
    fmpz_mul_ui(t2, t2, 2);

    evalf(t3, t1, a, b);

/*printf("t3: "); fmpz_print(t3); printf("\n");*/

    if (fmpz_is_zero(t3))
    {
        fmpz_swap(r + 0, t1);
        ret = 1;
        goto cleanup;
    }
    else if (fmpz_sgn(t3) < 0)
    {
        fmpz_sub_ui(t1, t1, 1);
        evalf(t3, t1, a, b);
        if (fmpz_is_zero(t3))
        {
            fmpz_swap(r + 0, t1);
            ret = 1;
            goto cleanup;
        }
        else if (fmpz_sgn(t3) < 0)
        {
            ret = 0;
            goto cleanup;
        }
    }

    /* f(t1) > 0, f(t2) < 0 now */
/*printf("finding most negative real root\n");*/
    if (find_cubic_root(r + 0, t2, t1, a, b, -1))
    {
        ret = 1;
        goto cleanup;
    }

    fmpz_mul_ui(t3, a, 3);
    fmpz_fdiv_q(t4, b, t3);
    fmpz_neg(t4, t4);

    fmpz_mul_ui(t3, a, 2);
    fmpz_cdiv_q(t2, b, t3);
    fmpz_neg(t2, t2);

    if (fmpz_cmp(t2, t1) < 0)
        fmpz_swap(t2, t1);

    /* f(t2) > 0, f(t4) < 0 now */
/*printf("finding slightly negative real root\n");*/
    ret = find_cubic_root(r + 0, t2, t4, a, b, +1);
    goto cleanup;

cleanup:

    if (sign_b < 0)
    {
        fmpz_neg(r + 0, r + 0);
        fmpz_neg(r + 1, r + 1);
        fmpz_neg(r + 2, r + 2);
    }

/*
printf("returning %d\n", ret);
printf("r[0]: "); fmpz_print(r + 0); printf("\n");
printf("r[1]: "); fmpz_print(r + 1); printf("\n");
printf("r[2]: "); fmpz_print(r + 2); printf("\n");
*/

    fmpz_clear(d);
    fmpz_clear(t1);
    fmpz_clear(t2);
    fmpz_clear(t3);
    fmpz_clear(t4);

    return ret;
}

void _fmpz_poly_factor_cubic(fmpz_poly_factor_t fac, 
                                              const fmpz_poly_t f, slong exp)
{
    slong i;
    fmpz_t A,B,b2,T,ac;
    fmpz r[3];
    const fmpz *a, *b, *c, *d;
    fmpz_poly_t p;

    FLINT_ASSERT(f->length == 4);
    FLINT_ASSERT(fmpz_sgn(f->coeffs + 3) > 0);

    d = f->coeffs;
    c = f->coeffs + 1;
    b = f->coeffs + 2;
    a = f->coeffs + 3;

    fmpz_init(A);
    fmpz_init(B);
    fmpz_init(b2);
    fmpz_init(T);
    fmpz_init(ac);
    fmpz_init(r + 0);
    fmpz_init(r + 1);
    fmpz_init(r + 2);
    fmpz_poly_init2(p, 3);
    
    /* A = b^2 - 3*a*c */
    fmpz_mul(ac, a, c);
    fmpz_mul(b2, b, b);
    fmpz_set(A, b2);
    fmpz_submul_ui(A, ac, 3);

    /* B = (9*a*c-2*b^2)*b - 27*a^2*d */
    fmpz_mul_ui(b2, b2, 2);
    fmpz_mul_ui(ac, ac, 9);
    fmpz_sub(B, ac, b2);
    fmpz_mul(B, B, b);
    fmpz_mul(T, a, a);
    fmpz_mul(T, T, d);
    fmpz_submul_ui(B, T, 27);

    switch (_fmpz_cubic_roots(r, A, B))
    {
        case 4:
            fmpz_mul_ui(p->coeffs + 1, a, 3);
            fmpz_sub(p->coeffs + 0, b, r + 0);
            fmpz_gcd(T, p->coeffs + 0, p->coeffs + 1);
            fmpz_divexact(p->coeffs + 0, p->coeffs + 0, T);
            fmpz_divexact(p->coeffs + 1, p->coeffs + 1, T);
            _fmpz_poly_set_length(p, 2);
            fmpz_poly_factor_insert(fac, p, 3*exp);
            break;

        case 3:
            fmpz_mul_ui(p->coeffs + 1, a, 3);
            fmpz_sub(p->coeffs + 0, b, r + 0);
            fmpz_gcd(T, p->coeffs + 0, p->coeffs + 1);
            fmpz_divexact(p->coeffs + 0, p->coeffs + 0, T);
            fmpz_divexact(p->coeffs + 1, p->coeffs + 1, T);
            _fmpz_poly_set_length(p, 2);
            fmpz_poly_factor_insert(fac, p, 1*exp);

            fmpz_mul_ui(p->coeffs + 1, a, 3);
            fmpz_sub(p->coeffs + 0, b, r + 1);
            fmpz_gcd(T, p->coeffs + 0, p->coeffs + 1);
            fmpz_divexact(p->coeffs + 0, p->coeffs + 0, T);
            fmpz_divexact(p->coeffs + 1, p->coeffs + 1, T);
            _fmpz_poly_set_length(p, 2);
            fmpz_poly_factor_insert(fac, p, 2*exp);

            break;

        case 2:
            for (i = 0; i < 3; i++)
            {
                fmpz_mul_ui(p->coeffs + 1, a, 3);
                fmpz_sub(p->coeffs + 0, b, r + i);
                fmpz_gcd(T, p->coeffs + 0, p->coeffs + 1);
                fmpz_divexact(p->coeffs + 0, p->coeffs + 0, T);
                fmpz_divexact(p->coeffs + 1, p->coeffs + 1, T);
                _fmpz_poly_set_length(p, 2);
                fmpz_poly_factor_insert(fac, p, 1*exp);
            }

            break;

        case 1:
            fmpz_mul_ui(p->coeffs + 1, a, 3);
            fmpz_sub(p->coeffs + 0, b, r + 0);
            fmpz_gcd(T, p->coeffs + 0, p->coeffs + 1);
            fmpz_divexact(p->coeffs + 0, p->coeffs + 0, T);
            fmpz_divexact(p->coeffs + 1, p->coeffs + 1, T);
            _fmpz_poly_set_length(p, 2);
            fmpz_poly_factor_insert(fac, p, 1*exp);
            fmpz_poly_divides(p, f, p);
            fmpz_poly_factor_insert(fac, p, 1*exp);
            break;

        default:
            fmpz_poly_factor_insert(fac, f, exp);
    }
     
    fmpz_clear(A);
    fmpz_clear(B);
    fmpz_clear(b2);
    fmpz_clear(T);
    fmpz_clear(ac);
    fmpz_clear(r + 0);
    fmpz_clear(r + 1);
    fmpz_clear(r + 2);

    fmpz_poly_clear(p);
}
