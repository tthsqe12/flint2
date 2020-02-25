/*
    Copyright (C) 2011 Andy Novocin
    Copyright (C) 2011 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include <stdlib.h>
#include "fmpz_poly.h"

void _fmpz_poly_factor_self_reciprocal(
    fmpz_poly_factor_t ffac,
    slong exp,
    const fmpz_poly_t f)
{
    slong i, j, k;
    slong half = (f->length - 1)/2;
    fmpz_poly_t g;
    fmpz_poly_factor_t gfac;
    fmpz * A, * B, * T;
    slong search_degrees[1];

    FLINT_ASSERT(half >= 2);
    FLINT_ASSERT(2*half == f->length - 1);

    A = _fmpz_vec_init(half + 2);
    B = _fmpz_vec_init(half + 2);

    fmpz_poly_init(g);
    fmpz_poly_factor_init(gfac);

    /* We need the solutions to a[0] = 2, a[1] = y, a[i] = y*a[i-1] - a[i-2].
       a[i] will be at A+1, a[i-1] will be at B+1
    */
    _fmpz_vec_zero(B, half + 2);
    _fmpz_vec_zero(A, half + 2);
    fmpz_set_si(B+1 + 0, 2);    /* a[0] = 2 */
    fmpz_set_si(A+1 + 1, 1);    /* a[1] = y */

    /* set g to the deflation of f */
    fmpz_poly_fit_length(g, f->length);
    g->length = half + 1;
    _fmpz_vec_zero(g->coeffs, g->length);
    fmpz_set(g->coeffs + 0, f->coeffs + half);
    fmpz_set(g->coeffs + 1, f->coeffs + half - 1);
    for (i = 2; i <= half; i++)
    {
        FLINT_ASSERT(fmpz_equal(f->coeffs + half - i, f->coeffs + half + i));
        _fmpz_vec_sub(B+1, A, B+1, i + 1); /* use A[0] = 0 */
        T = A; A = B; B = T;
        _fmpz_vec_scalar_addmul_fmpz(g->coeffs, A+1, i + 1, f->coeffs + half - i);
    }

    /* normalize g for the preconditions of _fmpz_poly_factor */
    if (fmpz_is_zero(g->coeffs + 0))
    {
        fmpz_poly_t t;
        fmpz_poly_shift_right(g, g, 1);
        fmpz_poly_init(t);
        fmpz_poly_set_coeff_si(t, 0, 1);
        fmpz_poly_set_coeff_si(t, 2, 1);
        fmpz_poly_factor_insert(ffac, t, exp);
        fmpz_poly_clear(t);
    }
    FLINT_ASSERT(!fmpz_is_zero(g->coeffs + 0));

    if (fmpz_sgn(g->coeffs + g->length - 1) < 0)
        fmpz_poly_neg(g, g);

    _fmpz_poly_factor_irreducible(gfac, 1, g,
                                FMPZ_POLY_FACTOR_USE_SELF_RECIPROCAL, NULL, 0);

    /* inflate each factor and factor again with restrictions */
    if (gfac->num > 1)
    {
        for (j = 0; j < gfac->num; j++)
        {
            fmpz_poly_struct * p = gfac->p + j;
            slong pdeg = p->length - 1;

            /* set g to the inflation of p */
            FLINT_ASSERT(pdeg >= 1);
            g->length = 2*pdeg + 1;
            _fmpz_vec_zero(g->coeffs, g->length);
            fmpz_one(A+1 + 0);
            fmpz_set(g->coeffs + pdeg, p->coeffs + 0);
            for (i = 1; i <= pdeg; i++)
            {
                fmpz_zero(A+1 + i);
                for (k = i; k >= 0; k--)
                {
                    fmpz_add(A+1 + k, A+1 + k, A+1 + k - 1); /* use A[0] = 0 */
                    fmpz_addmul(g->coeffs + pdeg - i + 2*k,
                                                 A+1 + k, p->coeffs + i);
                }
            }

            search_degrees[0] = pdeg;
            _fmpz_poly_factor_irreducible(ffac, exp, g, 0, search_degrees, 1);
        }
    }
    else
    {
        /* no need to reinflate if g was irreducible */
        search_degrees[0] = half;
        _fmpz_poly_factor_irreducible(ffac, exp, f, 0, search_degrees, 1);
    }

    _fmpz_vec_clear(B, half + 2);
    _fmpz_vec_clear(A, half + 2);

    fmpz_poly_clear(g);
    fmpz_poly_factor_clear(gfac);
}



static void add_possibilities(unsigned int * pos, slong len, int i,
                                                  const nmod_poly_factor_t fac)
{
    slong n, j, d;

    FLINT_ASSERT(i < 8*sizeof(int));

    pos[0] |= 1 << i;

    for (n = 0; n < fac->num; n++)
    {
        d = nmod_poly_degree(fac->p + n);
        for (j = len - d - 1; j >= 0; j--)
        {
            if (pos[j] & (1 << i))
                pos[j + d] |= (1 << i);
        }
    }
}


void _fmpz_poly_factor_irreducible(
    fmpz_poly_factor_t ffac,
    slong exp,
    const fmpz_poly_t f,
    unsigned int options,
    slong * search_degrees, slong num_search_degrees)
{
    slong cutoff = 10;
    slong num_possibilities;
    const slong lenF = f->length;
    slong i, idx;
    slong r = lenF;
    mp_limb_t p = 2;
    nmod_poly_t d, g, t;
    nmod_poly_factor_t fac;
    unsigned int * possible_degs;
    fmpz tmp[2];

    FLINT_ASSERT(lenF >= 2);
    FLINT_ASSERT(!fmpz_is_zero(f->coeffs + 0));
    FLINT_ASSERT(fmpz_sgn(f->coeffs + lenF - 1) > 0);

    if (lenF == 2)
    {
        fmpz_poly_factor_insert(ffac, f, exp);
        return;
    }

    /* check for a root at +-1 */
    fmpz_init(tmp + 0);
    fmpz_init(tmp + 1);
    for (i = 0; i < lenF; i++)
        fmpz_add(tmp + (i&1), tmp + (i&1), f->coeffs + i);

    if (fmpz_cmpabs(tmp + 0, tmp + 1) == 0)
    {
        int one_is_root = !fmpz_equal(tmp + 0, tmp + 1);
        fmpz_poly_t h;
        fmpz_poly_init(h);

        fmpz_poly_set_coeff_si(h, 1, 1);
        fmpz_poly_set_coeff_si(h, 0, one_is_root ? -1 : 1);
        fmpz_poly_factor_insert(ffac, h, exp);

        fmpz_poly_fit_length(h, lenF);
        h->length = lenF - 1;

        i = lenF - 1;
        fmpz_zero(h->coeffs + i);
        for ( ; i > 0; i--)
        {
            one_is_root ? fmpz_add(h->coeffs + i - 1, f->coeffs + i, h->coeffs + i)
                        : fmpz_sub(h->coeffs + i - 1, f->coeffs + i, h->coeffs + i);
        }

        _fmpz_poly_factor_irreducible(ffac, exp, h,
                                  options, search_degrees, num_search_degrees);

        fmpz_clear(tmp + 0);
        fmpz_clear(tmp + 1);
        fmpz_poly_clear(h);
        return;
    }
    fmpz_clear(tmp + 0);
    fmpz_clear(tmp + 1);

    /* check for self-reciprocal of even degree */
    if ((options & FMPZ_POLY_FACTOR_USE_SELF_RECIPROCAL) &&
        lenF > 4 && (lenF & 1) == 1)
    {
        slong half = (lenF - 1)/2;
        int ok = 1;
        for (i = 1; ok && i <= half; i++)
            ok = fmpz_equal(f->coeffs + half + i, f->coeffs + half - i);

        if (ok)
        {
            _fmpz_poly_factor_self_reciprocal(ffac, exp, f);
            return;
        }
    }

    possible_degs = (unsigned int *) flint_calloc(lenF + 1, sizeof(unsigned int));

    nmod_poly_factor_init(fac);
    nmod_poly_init_preinv(t, 1, 0);
    nmod_poly_init_preinv(d, 1, 0);
    nmod_poly_init_preinv(g, 1, 0);

    /* try three good primes */
    for (idx = 0; idx < 3; idx++)
    {
        for ( ; ; p = n_nextprime(p, 0))
        {
            nmod_init(&d->mod, p);
            g->mod = d->mod;
            t->mod = d->mod;

            fmpz_poly_get_nmod_poly(t, f);
            if (t->length == lenF && t->coeffs[0] != 0)
            {
                nmod_poly_derivative(d, t);
                nmod_poly_gcd(g, t, d);

                if (nmod_poly_is_one(g))
                {
                    nmod_poly_factor_t temp_fac;

                    nmod_poly_factor_init(temp_fac);
                    nmod_poly_factor(temp_fac, t);

                    add_possibilities(possible_degs, lenF, idx, temp_fac);

                    if (temp_fac->num <= r)
                    {
                        r = temp_fac->num;
                        nmod_poly_factor_set(fac, temp_fac);
                    }
                    nmod_poly_factor_clear(temp_fac);
                    break;
                }
            }
        }
        p = n_nextprime(p, 0);
    }
    nmod_poly_clear(d);
    nmod_poly_clear(g);
    nmod_poly_clear(t);

    for (i = 0; i < num_search_degrees; i++)
    {
        if (search_degrees[i] < lenF)
        {
            possible_degs[0] |= (1 << idx);
            possible_degs[search_degrees[i]] |= (1 << idx);
            idx++;
        }
    }

    num_possibilities = 0;
    for (i = 1; i < lenF - 1; i++)
    {
        num_possibilities += possible_degs[i] == possible_degs[0];
        #if TRACE_ZASSENHAUS == 1
        if (possible_degs[i] == possible_degs[0])
            flint_printf("degree %wd is possible\n", i);
        #endif
    }

    #if TRACE_ZASSENHAUS == 1
    flint_printf("num_possibilities: %wd\n", num_possibilities);
    #endif

    p = (fac->p + 0)->mod.n;
        
    if (r == 1 || num_possibilities == 0)
    {
        fmpz_poly_factor_insert(ffac, f, exp);
    }
    else if (r > cutoff)
    {
        fmpz_poly_factor_van_hoeij(ffac, fac, f, exp, p);
    }
    else
    {
        slong a;
        fmpz_t B;
        fmpz_poly_factor_t lifted_fac;

        fmpz_init(B);
        fmpz_poly_factor_init(lifted_fac);

        fmpz_poly_factor_mignotte(B, f);
        fmpz_mul_ui(B, B, 2);
        fmpz_add_ui(B, B, 1);
        a = fmpz_clog_ui(B, p);

        fmpz_poly_hensel_lift_once(lifted_fac, f, fac, a);

        fmpz_set_ui(B, p);
        fmpz_pow_ui(B, B, a);
        fmpz_poly_factor_zassenhaus_recombination(ffac,
                                         possible_degs, lifted_fac, f, B, exp);

        fmpz_clear(B);
        fmpz_poly_factor_clear(lifted_fac);
    }

    nmod_poly_factor_clear(fac);

    flint_free(possible_degs);
}



void fmpz_poly_factor(fmpz_poly_factor_t fac, const fmpz_poly_t G)
{
    const slong lenG = G->length;
    fmpz_poly_t g;

    fac->num = 0;

    if (lenG == 0)
    {
        fmpz_set_ui(&fac->c, 0);
        return;
    }
    if (lenG == 1)
    {
        fmpz_set(&fac->c, G->coeffs);
        return;
    }

    fmpz_poly_init(g);

    if (lenG == 2)
    {
        fmpz_poly_content(&fac->c, G);
        if (fmpz_sgn(fmpz_poly_lead(G)) < 0)
            fmpz_neg(&fac->c, &fac->c);
        fmpz_poly_scalar_divexact_fmpz(g, G, &fac->c);
        fmpz_poly_factor_insert(fac, g, 1);
    }
    else
    {
        slong j, k;
        fmpz_poly_factor_t sq_fr_fac;

        /* Does a presearch for a factor of form x^k */
        for (k = 0; fmpz_is_zero(G->coeffs + k); k++) ;

        if (k != 0)
        {
            fmpz_poly_t t;

            fmpz_poly_init(t);
            fmpz_poly_set_coeff_ui(t, 1, 1);
            fmpz_poly_factor_insert(fac, t, k);
            fmpz_poly_clear(t);
        }

        fmpz_poly_shift_right(g, G, k);

        fmpz_poly_factor_init(sq_fr_fac);
        fmpz_poly_factor_squarefree(sq_fr_fac, g);

        fmpz_set(&fac->c, &sq_fr_fac->c);

        /* Factor each square-free part */
        for (j = 0; j < sq_fr_fac->num; j++)
            _fmpz_poly_factor_irreducible(fac, sq_fr_fac->exp[j], sq_fr_fac->p + j,
                                FMPZ_POLY_FACTOR_USE_SELF_RECIPROCAL, NULL, 0);

        fmpz_poly_factor_clear(sq_fr_fac);
    }
    fmpz_poly_clear(g);
}

