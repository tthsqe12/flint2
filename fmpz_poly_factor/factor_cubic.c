/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fmpz_poly_factor.h"




static mp_limb_t fmpz_mod_64(const fmpz_t a)
{
    if (COEFF_IS_MPZ(*a))
    {
        const __mpz_struct * A = COEFF_TO_PTR(*a);
        return A->_mp_size > 0 ? A->_mp_d[0] : -A->_mp_d[0];
    }
    else
    {
        return *a;
    }
}


/*
    factor 2^e*y^3 + a*y + b = (y + r)*(2^e*y^2 - r*y + s) mod 2^p

    assuming:
        a and b are both odd and e > 0, or
        a is even and b is odd and e = 0.
*/

void binary_cubic_lift(
    fmpz_t r,
    fmpz_t s,
    const fmpz_t a,
    const fmpz_t b,
    slong e,
    slong p)
{
    slong n;
    fmpz_t r2, c, d, inv, t, t2;
    mp_limb_t A, B, C, D, INV, R, R2, S, E;
/*
printf("--- binary_cubic_root called ---\n");
flint_printf("      a: "); fmpz_print_binary(a); printf("\n");
flint_printf("      b: "); fmpz_print_binary(b); printf("\n");
flint_printf("      t: %wd\n", t);
flint_printf("   prec: %wd\n", p);
*/
    n = 1;

FLINT_ASSERT(p >= FLINT_BITS);

    A = fmpz_mod_64(a);
    B = fmpz_mod_64(b);
    R = 1;
    S = 1;
    INV = 1;
    R2 = R*R;
    E = (e < FLINT_BITS) ? (UWORD(1) << e) : 0;

    while (n <= FLINT_BITS/2)
    {
        mp_limb_t mask = (UWORD(1) << n);
        C = (A - (S - R2*E)) >> n;
        D = (B - (R*S)) >> n;
        R += (((D - C*R)*INV) % mask) << n;
        S += (((2*E*D*R + C*S)*INV) % mask) << n;
        n *= 2;
        R2 = R*R;
        INV = 2*INV - (INV*INV)*(2*R2*E + S);
    }

    fmpz_set_ui(r, R);
    fmpz_set_ui(s, S);

    if (n >= p)
        return;

    fmpz_init(t);
    fmpz_init(t2);
    fmpz_init(c);
    fmpz_init(d);
    fmpz_init_set_ui(inv, INV);
    fmpz_init_set_ui(r2, R);
    fmpz_mul_ui(r2, r2, R); 

    while (n < p)
    {
        fmpz_mul_2exp(c, r2, e);
        fmpz_add(c, c, a);
        fmpz_sub(c, c, s);
        fmpz_fdiv_q_2exp(c, c, n);

        fmpz_set(d, b);
        fmpz_submul(d, r, s);
        fmpz_fdiv_q_2exp(d, d, n);

        fmpz_mul(t, d, r);
        fmpz_mul_2exp(t, t, e + 1);
        fmpz_addmul(t, c, s);
        fmpz_fdiv_r_2exp(t, t, n);
        fmpz_mul(t, t, inv);
        fmpz_fdiv_r_2exp(t, t, n);
        fmpz_mul_2exp(t, t, n);
        fmpz_add(s, s, t);

        fmpz_submul(d, c, r);
        fmpz_fdiv_r_2exp(d, d, n);
        fmpz_mul(d, d, inv);
        fmpz_fdiv_r_2exp(d, d, n);
        fmpz_mul_2exp(d, d, n);
        fmpz_add(r, r, d);

        n *= 2;

        if (n < p)
        {
            fmpz_mul(r2, r, r);

            fmpz_mul(t, inv, inv);
            fmpz_mul_2exp(t2, r2, e + 1);
            fmpz_add(t2, t2, s);
            fmpz_mul_2exp(inv, inv, 1);
            fmpz_submul(inv, t, t2);
            fmpz_fdiv_r_2exp(inv, inv, n);
        }
    }
/*
    fmpz_mul(r2, r, r);
    fmpz_mul_2exp(c, r2, e);
    fmpz_add(c, c, a);
    fmpz_sub(c, c, s);
    fmpz_fdiv_r_2exp(t, c, p);
    FLINT_ASSERT(fmpz_is_zero(t));

    fmpz_set(d, b);
    fmpz_submul(d, r, s);
    fmpz_fdiv_r_2exp(t, d, p);
    FLINT_ASSERT(fmpz_is_zero(t));
*/
    fmpz_clear(t);
    fmpz_clear(t2);
    fmpz_clear(c);
    fmpz_clear(d);
    fmpz_clear(inv);
    fmpz_clear(r2);
/*
printf("--- binary_cubic_root returning ---\n");
flint_printf("      r: "); fmpz_print_binary(r); printf("\n");
flint_printf("      s: "); fmpz_print_binary(s); printf("\n");
*/
}

/* z = sqrt(x) clobber x */
slong binary_sqrt_lift(fmpz_t z, fmpz_t x, slong p)
{
    slong e, mp, n;
    fmpz_t t, t2, t3, t4, m;
    fmpz two = 2;

    FLINT_ASSERT(p > 0);
    fmpz_fdiv_r_2exp(x, x, p);

/*
printf("--- binary_sqrt called ---\n");
flint_printf("      x: "); fmpz_print_binary(x); printf("\n");
flint_printf("   prec: %wd\n", p);
*/

    if (fmpz_is_zero(x))
    {
        fmpz_zero(z);
        return p/2;
    }

    e = fmpz_remove(x, x, &two);

    if ((e % 2) != 0)
    {
        fmpz_zero(z);
/*
printf("--- binary_sqrt returning ---\n");
*/
        return -WORD(1);
    }

    mp = p - e - 3;
    if (mp < 1)
    {
        fmpz_one(z);
        fmpz_mul_2exp(z, z, e/2);
/*
printf("--- binary_sqrt returning ---\n");
*/
        return e/2 + 1;
    }

    if (fmpz_fdiv_ui(x, 8) != 1)
    {
        fmpz_zero(z);
/*
printf("--- binary_sqrt returning ---\n");
*/
        return -WORD(1);
    }

    fmpz_init(t);
    fmpz_init(t2);
    fmpz_init(t3);
    fmpz_init(t4);
    fmpz_init(m);

    fmpz_fdiv_q_ui(m, x, 8);
    fmpz_set(z, m);
/*
flint_printf("      z: "); fmpz_print_binary(z); printf("\n");
*/
    n = 4;
    while (n - 1 < mp + 2)
    {
        n = 2*n - 2;

        fmpz_mul(t2, z, z);

        fmpz_mul_2exp(t, z, 2);
        fmpz_sub_ui(t4, t, 3);
        fmpz_mul_si(t4, t4, -2);
        fmpz_mul(t4, t4, t2);

        fmpz_mul_2exp(t, z, 2);
        fmpz_sub_ui(t3, t, 1);
        fmpz_pow_ui(t3, t3, 3);

        fmpz_fdiv_r_2exp(t2, m, n);


        fmpz_mul(t3, t3, t2);

        fmpz_sub(z, t4, t3);
        fmpz_fdiv_r_2exp(z, z, n);
/*
flint_printf("      z: "); fmpz_print_binary(z); printf("\n");
*/
    }

    fmpz_mul_si(t, z, -4);
    fmpz_add_ui(t, t, 1);
    fmpz_mul(z, t, x);
    fmpz_fdiv_r_2exp(z, z, mp + 2);
    fmpz_mul_2exp(z, z, e/2);

    fmpz_clear(t);
    fmpz_clear(t2);
    fmpz_clear(t3);
    fmpz_clear(t4);
    fmpz_clear(m);
/*
printf("--- binary_sqrt returning ---\n");
flint_printf("sqrt(x): "); fmpz_print_binary(z); printf("\n");
flint_printf("   prec: %wd\n", mp + 2 + e/2);
*/
    return mp + 2 + e/2;
}

/* return f(0)*...*f(largest_prime - 1) mod prime_product */
static mp_limb_t eval_product_mod_n(
    const fmpz_t a,
    const fmpz_t b,
    mp_limb_t prime_product,
    mp_limb_t largest_prime)
{
    nmod_t ctx;
    ulong A, B, F, G, H, P, i;
    nmod_init(&ctx, prime_product);

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
    for (i = 1; i < largest_prime; i++)
    {
        F = nmod_add(F, G, ctx);
        G = nmod_sub(G, H, ctx);
        H += 6;
        P = nmod_mul(P, F, ctx);
    }

    return P;
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
int _fmpz_cubic_roots(fmpz * x, fmpz_t a, fmpz_t b)
{
    slong i;
    int ret, sign_a, sign_b, sign_d;
    fmpz_t d, t1, t2, t3, t4, ta, tb, r, s, z;
    slong prec;
    slong new_prec, sqrt_prec;
    ulong alpha, beta, alpha2, beta3;
    fmpz two = 2;



/*
printf("--------- cubic roots called -----------\n");
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
            fmpz_swap(x + 0, b);
            return 1;
        }
        else
        {
            return 0;
        }
    }

    if (fmpz_is_zero(b))
    {
        fmpz_zero(x + 0);

        if (sign_a <= 0)
            return sign_a == 0 ? 4 : 1;
        
        fmpz_mul_ui(a, a, 3);
        
        if (fmpz_is_square(a))
        {
            fmpz_sqrt(x + 1, a);
            fmpz_neg(x + 2, x + 1);
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
    if (0 != eval_product_mod_n(a, b, UWORD(5)*7*11*13*17*19*23*29, 29))
        return 0;

    fmpz_init(d);
    fmpz_init(t1);
    fmpz_init(t2);
    fmpz_init(t3);
    fmpz_init(t4);
    fmpz_init(ta);
    fmpz_init(tb);
    fmpz_init(r);
    fmpz_init(s);
    fmpz_init(z);

    /* d = b^2 - 4*a^3 = discriminant/-27 */
    fmpz_mul(d, b, b);
    fmpz_mul(t1, a, a);
    fmpz_mul(t2, t1, a);
    fmpz_submul_ui(d, t2, 4);

    sign_d = fmpz_sgn(d);

    if (sign_d == 0)
    {
        FLINT_ASSERT(fmpz_divisible(b, a));
        fmpz_divexact(x + 0, b, a);
        FLINT_ASSERT(fmpz_is_even(x + 0));
        fmpz_divexact_si(x + 1, x + 0, -2);
        ret = 3;
        goto cleanup;
    }

    if (sign_d > 0)
    {
        /*
            An interval containing the real root is (0, max(a,b)).
            Unfortunately newton converges too slowly starting from x = max(a,b).
            Luckily, we have a formula for the real root,

                cbrt(4*b + 4*sqrt(d)) + cbrt(4*b - 4*sqrt(d))
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
        fmpz_fdiv_q_2exp(x + 0, t1, 1);

        /*
            If r is the actual real root in RR, then, experimentally,
            r[0] - 1 < r < r[0] + 2, so that r[0] - 1 cannot be a root. Test
            all three numbers r[0] - 1, r[0], r[0] + 1 anyways since these
            follow from easily-proven bounds.
        */

        fmpz_mul(t2, x + 0, x + 0);
        fmpz_sub(t3, t2, a);
        fmpz_submul_ui(t2, a, 3);
        fmpz_mul(t1, t2, x + 0);
        if (fmpz_equal(t1, b))
        {
            ret = 1;
            goto cleanup;
        }

        fmpz_submul_ui(b, x + 0, 3);
        fmpz_mul_ui(t3, t3, 3);
        fmpz_add_ui(t3, t3, 1);
        fmpz_sub(t2, b, t3);
        if (fmpz_equal(t1, t2))
        {
            fmpz_add_ui(x + 0, x + 0, 1);
            ret = 1;
            goto cleanup;
        }

        fmpz_add(t2, b, t3);
        if (fmpz_equal(t1, t2))
        {
            fmpz_sub_ui(x + 0, x + 0, 1);
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
        give intervals on which newton coverges quickly over RR.
        However, we still do about log(a) full-precision computations, which
        kills the asymptotic performance. So, the lifting takes place in ZZ_2.
    */

    /* check irreducibility mod some more p */
    if (0 != eval_product_mod_n(a, b, UWORD(31)*37*41*43*47, 47))
    {
        ret = 0;
        goto cleanup;
    }

    /* 2*sqrt(a) is bound on absolute value of roots */
    prec = fmpz_bits(a)/2 + 3;
    prec = (prec + FLINT_BITS - 1)/FLINT_BITS * FLINT_BITS;

    fmpz_mul_si(ta, a, -3);
    fmpz_mul_si(tb, b, -1);

    FLINT_ASSERT(!fmpz_is_zero(ta));
    FLINT_ASSERT(!fmpz_is_zero(tb));

    alpha = fmpz_remove(ta, ta, &two);
    beta = fmpz_remove(tb, tb, &two);

    if (3*alpha == 2*beta)
    {
        ret = 0;
        goto cleanup;
    }
    else if (3*alpha > 2*beta)
    {
        /* only one root, which has valuation beta/3 */

        beta3 = beta/3;
        if ((beta % 3) != 0)
        {
            ret = 0;
            goto cleanup;
        }

        fmpz_mul_2exp(ta, ta, alpha - 2*beta3);
        binary_cubic_lift(r, s, ta, tb, 0, prec);
        fmpz_mul_2exp(x + 0, r, beta3);
        fmpz_neg(x + 0, x + 0);
        fmpz_fdiv_r_2exp(x + 0, x + 0, prec);

        if (fmpz_bits(x + 0) >= prec)
        {
            fmpz_neg(x + 0, x + 0);
            fmpz_fdiv_r_2exp(x + 0, x + 0, prec);
            fmpz_neg(x + 0, x + 0);
        }
        fmpz_mul(t1, x + 0, x + 0);
        fmpz_submul_ui(t1, a, 3);
        fmpz_mul(t2, t1, x + 0);
        ret = fmpz_equal(t2, b) ? 1 : 0;
        goto cleanup;
    }
    else
    {
        /* there is a root with valuation beta - alpha */

        alpha2 = alpha/2;
        if ((alpha % 2) != 0)
        {
            /* there are no other roots */

            binary_cubic_lift(r, s, ta, tb, 2*beta - 3*alpha, prec);
            fmpz_mul_2exp(x + 0, r, beta - alpha);
            fmpz_neg(x + 0, x + 0);
            fmpz_fdiv_r_2exp(x + 0, x + 0, prec);
            if (fmpz_bits(x + 0) >= prec)
            {
                fmpz_neg(x + 0, x + 0);
                fmpz_fdiv_r_2exp(x + 0, x + 0, prec);
                fmpz_neg(x + 0, x + 0);
            }

            fmpz_mul(t1, x + 0, x + 0);
            fmpz_submul_ui(t1, a, 3);
            fmpz_mul(t2, t1, x + 0);

            ret = fmpz_equal(t2, b) ? 1 : 0;
            goto cleanup;
        }

        new_prec = prec > FLINT_BITS ? prec + FLINT_BITS : prec;

        binary_cubic_lift(r, s, ta, tb, 2*beta - 3*alpha, new_prec);
        fmpz_mul_2exp(x + 0, r, beta - alpha);
        fmpz_neg(x + 0, x + 0);
        fmpz_fdiv_r_2exp(x + 0, x + 0, prec);
        if (fmpz_bits(x + 0) >= prec)
        {
            fmpz_neg(x + 0, x + 0);
            fmpz_fdiv_r_2exp(x + 0, x + 0, prec);
            fmpz_neg(x + 0, x + 0);
        }

        fmpz_mul(t1, x + 0, x + 0);
        fmpz_submul_ui(t1, a, 3);
        fmpz_mul(t2, t1, x + 0);

        /* there are two roots with valuation 2 */

        if (fmpz_equal(t2, b))
        {
            fmpz_mul_ui(a, a, 4);
            fmpz_submul(a, x + 0, x + 0);
            fmpz_mul_ui(a, a, 3);
            if (fmpz_is_square(a))
            {
                fmpz_sqrt(t1, a);
                fmpz_sub(x + 1, t1, x + 0);
                fmpz_add(x + 2, t1, x + 0);
                FLINT_ASSERT(fmpz_is_even(x + 1));
                FLINT_ASSERT(fmpz_is_even(x + 2));
                fmpz_divexact_ui(x + 1, x + 1, 2);
                fmpz_divexact_si(x + 2, x + 2, -2);
                ret = 2;
                goto cleanup;
            }
            else
            {
                ret = 1;
                goto cleanup;
            }
        }

        /* the root with valuation beta - alpha is irrational */

try_again:

        fmpz_mul(d, r, r);
        fmpz_mul_2exp(d, d, 2*beta - 3*alpha - 2);
        fmpz_sub(d, d, s);
        fmpz_fdiv_r_2exp(d, d, new_prec);
        sqrt_prec = binary_sqrt_lift(z, d, new_prec);

        if (sqrt_prec < 0)
        {
            /* roots with valuation alpha/2 are irrational*/
            ret = 0;
            goto cleanup;
        }

        if (sqrt_prec + alpha2 < prec)
        {
/*
flint_printf("%wd + %wd < %wd, beta: %wd, alpha: %wd\n", sqrt_prec, alpha2, prec, beta, alpha);
*/
            new_prec += FLINT_BITS + prec - (sqrt_prec + alpha2);
            binary_cubic_lift(r, s, ta, tb, 2*beta - 3*alpha, new_prec);
            goto try_again;
        }

        fmpz_mul_2exp(r, r, beta - alpha - 1);
        fmpz_mul_2exp(z, z, alpha2);
        fmpz_add(x + 1, r, z);
        fmpz_sub(x + 2, r, z);
        fmpz_fdiv_r_2exp(x + 1, x + 1, prec);
        fmpz_fdiv_r_2exp(x + 2, x + 2, prec);

        for (i = 1; i <= 2; i++)
        {
            if (fmpz_bits(x + i) >= prec)
            {
                fmpz_neg(x + i, x + i);
                fmpz_fdiv_r_2exp(x + i, x + i, prec);
                fmpz_neg(x + i, x + i);
            }

            fmpz_mul(t1, x + i, x + i);
            fmpz_submul_ui(t1, a, 3);
            fmpz_mul(t2, t1, x + i);
            if (fmpz_equal(t2, b))
            {
                fmpz_swap(x + i, x + 0);
                ret = 1;
                goto cleanup;
            }
        }

        ret = 0;
        goto cleanup;
    }

cleanup:

    if (sign_b < 0)
    {
        fmpz_neg(x + 0, x + 0);
        fmpz_neg(x + 1, x + 1);
        fmpz_neg(x + 2, x + 2);
    }
/*
printf("returning %d\n", ret);
printf("r[0]: "); fmpz_print(x + 0); printf("\n");
printf("r[1]: "); fmpz_print(x + 1); printf("\n");
printf("r[2]: "); fmpz_print(x + 2); printf("\n");
*/
    fmpz_clear(d);
    fmpz_clear(t1);
    fmpz_clear(t2);
    fmpz_clear(t3);
    fmpz_clear(t4);
    fmpz_clear(ta);
    fmpz_clear(tb);
    fmpz_clear(r);
    fmpz_clear(s);
    fmpz_clear(z);

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
