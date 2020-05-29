/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fmpz_poly_factor.h"

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
    int ret, sign_a, sign_b, sign_d;
    fmpz_t d, t1, t2, t3, t4;

printf("---------cubic roots called-----------\n");
printf("a: "); fmpz_print(a); printf("\n");
printf("b: "); fmpz_print(b); printf("\n");

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
        fmpz_divexact_si(r + 1, r + 0, -2);
        ret = 3;
        goto cleanup;
    }
    else if (sign_d > 0)
    {
        /* we have a formula for the real root */
        if (sign_a == 0)
        {
            fmpz_root(r + 0, b, 3);
        }
        else
        {
            fmpz_sqrt(t3, d);
            fmpz_add(t4, b, t3);
            fmpz_mul_ui(t4, t4, 4);
            fmpz_root(t1, t4, 3);
            if (sign_a > 0)
            {
                fmpz_sub(t4, b, t3);
                fmpz_mul_ui(t4, t4, 4);
                fmpz_root(t2, t4, 3);
                fmpz_add(t1, t1, t2);
            }
            else
            {
                fmpz_sub(t4, t3, b);
                fmpz_mul_ui(t4, t4, 4);
                fmpz_root(t2, t4, 3);                
                fmpz_sub(t1, t1, t2);
            }
            fmpz_fdiv_q_2exp(r + 0, t1, 1);
        }

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
    else
    {
        FLINT_ASSERT(0);
    }

cleanup:

    if (sign_b < 0)
    {
        fmpz_neg(r + 0, r + 0);
        fmpz_neg(r + 1, r + 1);
        fmpz_neg(r + 2, r + 2);
    }


printf("returning %d\n", ret);
printf("r[0]: "); fmpz_print(r + 0); printf("\n");
printf("r[1]: "); fmpz_print(r + 1); printf("\n");
printf("r[2]: "); fmpz_print(r + 2); printf("\n");


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
    
    fmpz_mul(ac, a, c);
    fmpz_mul(b2, b, b);
    fmpz_set(A, b2);
    fmpz_submul_ui(A, ac, 3);

    fmpz_mul_ui(b2, b2, 2);
    fmpz_mul_ui(ac, ac, 9);
    fmpz_sub(B, ac, b2);
    fmpz_mul(B, B, b);
    fmpz_mul(T, a, a);
    fmpz_mul(T, T, d);
    fmpz_mul_ui(T, T, 27);
    fmpz_sub(B, B, T);

    switch (_fmpz_cubic_roots(r, A, B))
    {
        case 4:
printf("case 4\n");
            fmpz_mul_ui(p->coeffs + 1, a, 3);
            fmpz_sub(p->coeffs + 0, b, r + 0);
            fmpz_gcd(T, p->coeffs + 0, p->coeffs + 1);
            fmpz_divexact(p->coeffs + 0, p->coeffs + 0, T);
            fmpz_divexact(p->coeffs + 1, p->coeffs + 1, T);
            _fmpz_poly_set_length(p, 2);
            fmpz_poly_factor_insert(fac, p, 3*exp);
            break;

        case 3:
printf("case 3\n");

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


        case 1:
printf("case 1\n");

            break;
        case 2: 
printf("case 2\n");

            break;


        default:
printf("default case\n");
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
}
