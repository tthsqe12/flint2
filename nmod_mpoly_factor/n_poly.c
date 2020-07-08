/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "nmod_mpoly_factor.h"
#include "mpn_extras.h"
#include "nmod_vec.h"


void n_poly_realloc(n_poly_t A, slong len)
{
    slong old_alloc = A->alloc;
    slong new_alloc = FLINT_MAX(len, old_alloc + 1 + old_alloc/2);

    FLINT_ASSERT(A->alloc >= 0);
    if (len <= A->alloc)
        return;

    if (old_alloc > 0)
    {
        FLINT_ASSERT(A->coeffs != NULL);
        A->coeffs = (mp_limb_t *) flint_realloc(A->coeffs,
                                                  new_alloc*sizeof(mp_limb_t));
    }
    else
    {
        FLINT_ASSERT(A->coeffs == NULL);
        A->coeffs = (mp_limb_t *) flint_malloc(new_alloc*sizeof(mp_limb_t));
    }
    A->alloc = new_alloc;
}


void n_poly_print_pretty(const n_poly_t A, const char * x)
{
    slong i;
    int first = 1;

    for (i = A->length - 1; i >= 0; i--)
    {
        if (!first)
            flint_printf(" + ");
        first = 0;
        flint_printf("%wu*%s^%wd", A->coeffs[i], x, i);
    }

    if (first)
        flint_printf("0");
}

int n_poly_mod_is_canonical(const n_poly_t A, nmod_t mod)
{
    slong i;
    if (A->length < 0)
        return 0;
    for (i = 0; i < A->length; i++)
    {
        if (A->coeffs[i] >= mod.n)
            return 0;
        if (A->coeffs[i] == 0 && i + 1 == A->length)
            return 0;
    }
    return 1;
}

void n_poly_set_coeff(n_poly_t poly, slong j, ulong c)
{
    n_poly_fit_length(poly, j + 1);

    if (j + 1 < poly->length) /* interior */
    {
        poly->coeffs[j] = c;
    }
    else if (j + 1 == poly->length) /* leading coeff */
    {
        if (c != 0)
        {
            poly->coeffs[j] = c;
        }
        else
        {
            poly->length--;
            _n_poly_normalise(poly);
        }
    }
    else if (c != 0) /* extend polynomial */
    {
        flint_mpn_zero(poly->coeffs + poly->length, j - poly->length);
        poly->coeffs[j] = c;
        poly->length = j + 1;
    }
}


void n_poly_mod_set_coeff_ui(
    n_poly_t poly,
    slong j,
    ulong c,
    nmod_t mod)
{
    if (c >= mod.n)
        NMOD_RED(c, c, mod);

    n_poly_set_coeff(poly, j, c);
}



void n_poly_mod_mul(n_poly_t res, const n_poly_t poly1,
                 const n_poly_t poly2, nmod_t mod)
{
    slong len1, len2, len_out;
    
    len1 = poly1->length;
    len2 = poly2->length;

    if (len1 == 0 || len2 == 0)
    {
        n_poly_zero(res);

        return;
    }

    len_out = poly1->length + poly2->length - 1;

    if (res == poly1 || res == poly2)
    {
        n_poly_t temp;

        n_poly_init2(temp, len_out);

        if (len1 >= len2)
            _nmod_poly_mul(temp->coeffs, poly1->coeffs, len1,
                           poly2->coeffs, len2, mod);
        else
            _nmod_poly_mul(temp->coeffs, poly2->coeffs, len2,
                           poly1->coeffs, len1, mod);
        
        n_poly_swap(temp, res);
        n_poly_clear(temp);
    }
    else
    {
        n_poly_fit_length(res, len_out);
        
        if (len1 >= len2)
            _nmod_poly_mul(res->coeffs, poly1->coeffs, len1,
                           poly2->coeffs, len2, mod);
        else
            _nmod_poly_mul(res->coeffs, poly2->coeffs, len2,
                           poly1->coeffs, len1, mod);
    }

    res->length = len_out;
    _n_poly_normalise(res);
}

void n_poly_mod_mullow(
    n_poly_t res, 
    const n_poly_t poly1,
    const n_poly_t poly2,
    slong trunc,
    nmod_t mod)
{
    slong len1, len2, len_out;
    
    len1 = poly1->length;
    len2 = poly2->length;

    len_out = poly1->length + poly2->length - 1;
    if (trunc > len_out)
        trunc = len_out;
    
    if (len1 <= 0 || len2 <= 0 || trunc <= 0)
    {
        n_poly_zero(res);
        return;
    }

    if (res == poly1 || res == poly2)
    {
        n_poly_t temp;

        n_poly_init2(temp, trunc);

        if (len1 >= len2)
            _nmod_poly_mullow(temp->coeffs, poly1->coeffs, len1,
                           poly2->coeffs, len2, trunc, mod);
        else
            _nmod_poly_mullow(temp->coeffs, poly2->coeffs, len2,
                           poly1->coeffs, len1, trunc, mod);
        
        n_poly_swap(temp, res);
        n_poly_clear(temp);
    }
    else
    {
        n_poly_fit_length(res, trunc);
        
        if (len1 >= len2)
            _nmod_poly_mullow(res->coeffs, poly1->coeffs, len1,
                           poly2->coeffs, len2, trunc, mod);
        else
            _nmod_poly_mullow(res->coeffs, poly2->coeffs, len2,
                           poly1->coeffs, len1, trunc, mod);
    }

    res->length = trunc;
    _n_poly_normalise(res);
}



void n_poly_mod_div(n_poly_t Q, const n_poly_t A, const n_poly_t B, nmod_t mod)
{
    n_poly_t tQ;
    mp_ptr q;
    slong A_len, B_len;

    B_len = B->length;
    
    if (B_len == 0)
    {
        if (mod.n == 1)
        {
            n_poly_set(Q, A);
            return;
        }
        else
        {                                                                                
            flint_printf("Exception (n_poly_mod_div). Division by zero.\n");
            flint_abort();
        }
    }

    A_len = A->length;
    
    if (A_len < B_len)
    {
        n_poly_zero(Q);
        return;
    }

    if (Q == A || Q == B)
    {
        n_poly_init2(tQ, A_len - B_len + 1);
        q = tQ->coeffs;
    }
    else
    {
        n_poly_fit_length(Q, A_len - B_len + 1);
        q = Q->coeffs;
    }

    _nmod_poly_div(q, A->coeffs, A_len, B->coeffs, B_len, mod);

    if (Q == A || Q == B)
    {
        n_poly_swap(tQ, Q);
        n_poly_clear(tQ);
    }
    
    Q->length = A_len - B_len + 1;
}

void n_poly_mod_rem(n_poly_t R, const n_poly_t A, const n_poly_t B, nmod_t mod)
{
    const slong lenA = A->length, lenB = B->length;
    n_poly_t tR;
    mp_ptr r;

    if (lenB == 0)
    {
        flint_printf("Exception (nmod_poly_rem). Division by zero.\n");
        flint_abort();
    }
    if (lenA < lenB)
    {
        n_poly_set(R, A);
        return;
    }

    if (R == A || R == B)
    {
        n_poly_init2(tR, lenB - 1);
        r = tR->coeffs;
    }
    else
    {
        n_poly_fit_length(R, lenB - 1);
        r = R->coeffs;
    }

    _nmod_poly_rem(r, A->coeffs, lenA, B->coeffs, lenB, mod);

    if (R == A || R == B)
    {
        n_poly_swap(R, tR);
        n_poly_clear(tR);
    }
        
    R->length = lenB - 1;
    _n_poly_normalise(R);
}


void n_poly_mod_divrem(n_poly_t Q, n_poly_t R,
                                const n_poly_t A, const n_poly_t B, nmod_t mod)
{
    const slong lenA = A->length, lenB = B->length;
    n_poly_t tQ, tR;
    mp_ptr q, r;
    
    if (lenB == 0)
    {
        if (mod.n == 1)
        {
            n_poly_set(Q, A);
            n_poly_zero(R);
            return;
        }
        else
        {
            flint_printf("Exception (n_poly_mod_divrem). Division by zero.");
            flint_abort();
        }
    }

    if (lenA < lenB)
    {
        n_poly_set(R, A);
        n_poly_zero(Q);
        return;
    }

    if (Q == A || Q == B)
    {
        n_poly_init2(tQ, lenA - lenB + 1);
        q = tQ->coeffs;
    }
    else
    {
        n_poly_fit_length(Q, lenA - lenB + 1);
        q = Q->coeffs;
    }

    if (R == A || R == B)
    {
        n_poly_fit_length(tR, lenB - 1);
        r = tR->coeffs;
    }
    else
    {
        n_poly_fit_length(R, lenB - 1);
        r = R->coeffs;
    }

    _nmod_poly_divrem(q, r, A->coeffs, lenA, B->coeffs, lenB, mod);

    if (Q == A || Q == B)
    {
        n_poly_swap(Q, tQ);
        n_poly_clear(tQ);
    }
    if (R == A || R == B)
    {
        n_poly_swap(R, tR);
        n_poly_clear(tR);
    }
        
    Q->length = lenA - lenB + 1;
    R->length = lenB - 1;

    _n_poly_normalise(R);
}


int n_poly_mod_invmod(n_poly_t A, const n_poly_t B, const n_poly_t P, nmod_t mod)
{
    const slong lenB = B->length, lenP = P->length;
    mp_limb_t * a;
    n_poly_t tA;
    int ans;

    if (lenP < 2)
    {
        printf("Exception (nmod_poly_invmod). lenP < 2.\n");
        flint_abort();
    }
    if (lenB == 0)
    {
        n_poly_zero(A);
        return 0;
    }
    if (lenB >= lenP)
    {
        n_poly_t T;

        n_poly_init(T);
        n_poly_mod_rem(T, B, P, mod);
        ans = n_poly_mod_invmod(A, T, P, mod);
        n_poly_clear(T);
        return ans;
    }

    if (A != B && A != P)
    {
        n_poly_fit_length(A, lenP - 1);
        a = A->coeffs;
    }
    else
    {
        n_poly_init2(tA, lenP - 1);
        a = tA->coeffs;
    }

    ans = _nmod_poly_invmod(a, B->coeffs, lenB, P->coeffs, lenP, mod);

    if (A == B || A == P)
    {
        n_poly_swap(A, tA);
        n_poly_clear(tA);
    }

    A->length = lenP - 1;
    _n_poly_normalise(A);
    return ans;
}


void n_poly_mod_gcd(n_poly_t G, const n_poly_t A, const n_poly_t B, nmod_t mod)
{
    if (A->length < B->length)
    {
        n_poly_mod_gcd(G, B, A, mod);
    }
    else /* lenA >= lenB >= 0 */
    {
        slong lenA = A->length, lenB = B->length, lenG;
        n_poly_t tG;
        mp_ptr g;

        if (lenA == 0) /* lenA = lenB = 0 */
        {
            n_poly_zero(G);
        } 
        else if (lenB == 0) /* lenA > lenB = 0 */
        {
            n_poly_mod_make_monic(G, A, mod);
        }
        else /* lenA >= lenB >= 1 */
        {
            if (G == A || G == B)
            {
                n_poly_init2(tG, FLINT_MIN(lenA, lenB));
                g = tG->coeffs;
            }
            else
            {
                n_poly_fit_length(G, FLINT_MIN(lenA, lenB));
                g = G->coeffs;
            }

            lenG = _nmod_poly_gcd(g, A->coeffs, lenA, B->coeffs, lenB, mod);

            if (G == A || G == B)
            {
                n_poly_swap(tG, G);
                n_poly_clear(tG);
            }
            G->length = lenG;

            if (G->length == 1)
                G->coeffs[0] = 1;
            else
                n_poly_mod_make_monic(G, G, mod);
        }
    }
}

void n_poly_mod_xgcd(
    n_poly_t G,
    n_poly_t S,
    n_poly_t T,
    const n_poly_t A,
    const n_poly_t B,
    nmod_t mod)
{
    if (A->length < B->length)
    {
        n_poly_mod_xgcd(G, T, S, B, A, mod);
    }
    else  /* lenA >= lenB >= 0 */
    {
        const slong lenA = A->length, lenB = B->length;
        mp_limb_t inv;

        if (lenA == 0)  /* lenA = lenB = 0 */
        {
            n_poly_zero(G);
            n_poly_zero(S);
            n_poly_zero(T);
        }
        else if (lenB == 0)  /* lenA > lenB = 0 */
        {
            inv = n_invmod(A->coeffs[lenA - 1], mod.n);
            n_poly_mod_scalar_mul_nmod(G, A, inv, mod);
            n_poly_zero(T);
            n_poly_set_coeff(S, 0, inv);
            S->length = 1;
        }
        else if (lenB == 1)  /* lenA >= lenB = 1 */
        {
            n_poly_fit_length(T, 1);
            T->length = 1;
            T->coeffs[0] = n_invmod(B->coeffs[0], mod.n);
            n_poly_one(G);
            n_poly_zero(S);
        }
        else  /* lenA >= lenB >= 2 */
        {
            mp_ptr g, s, t;
            slong lenG;

            if (G == A || G == B)
            {
                g = _nmod_vec_init(FLINT_MIN(lenA, lenB));
            }
            else
            {
                n_poly_fit_length(G, FLINT_MIN(lenA, lenB));
                g = G->coeffs;
            }
            if (S == A || S == B)
            {
                s = _nmod_vec_init(lenB - 1);
            }
            else
            {
                n_poly_fit_length(S, lenB - 1);
                s = S->coeffs;
            }
            if (T == A || T == B)
            {
                t = _nmod_vec_init(lenA - 1);
            }
            else
            {
                n_poly_fit_length(T, lenA - 1);
                t = T->coeffs;
            }

            if (lenA >= lenB)
                lenG = _nmod_poly_xgcd(g, s, t, A->coeffs, lenA,
                                                         B->coeffs, lenB, mod);
            else
                lenG = _nmod_poly_xgcd(g, t, s, B->coeffs, lenB,
                                                         A->coeffs, lenA, mod);

            if (G == A || G == B)
            {
                flint_free(G->coeffs);
                G->coeffs = g;
                G->alloc  = FLINT_MIN(lenA, lenB);
            }
            if (S == A || S == B)
            {
                flint_free(S->coeffs);
                S->coeffs = s;
                S->alloc  = lenB - 1;
            }
            if (T == A || T == B)
            {
                flint_free(T->coeffs);
                T->coeffs = t;
                T->alloc  = lenA - 1;
            }

            G->length = lenG;
            S->length = FLINT_MAX(lenB - lenG, 1);
            T->length = FLINT_MAX(lenA - lenG, 1);
            MPN_NORM(S->coeffs, S->length);
            MPN_NORM(T->coeffs, T->length);

            if (G->coeffs[lenG - 1] != 1)
            {
                inv = nmod_inv(G->coeffs[lenG - 1], mod);
                n_poly_mod_scalar_mul_nmod(G, G, inv, mod);
                n_poly_mod_scalar_mul_nmod(S, S, inv, mod);
                n_poly_mod_scalar_mul_nmod(T, T, inv, mod);
            }
        }
    }
}

void n_poly_mod_inv_series(n_poly_t Qinv, const n_poly_t Q, slong n, nmod_t mod)
{
    slong Qlen = Q->length;

    Qlen = FLINT_MIN(Qlen, n);

    if (Qlen == 0)
    {
        flint_printf("Exception (nmod_poly_inv_series_newton). Division by zero.\n");
        flint_abort();
    }

    if (Qinv != Q)
    {
        n_poly_fit_length(Qinv, n);
        _nmod_poly_inv_series_newton(Qinv->coeffs, Q->coeffs, Qlen, n, mod);
    }
    else
    {
        n_poly_t t;
        n_poly_init(t);
        _nmod_poly_inv_series_newton(t->coeffs, Q->coeffs, Qlen, n, mod);
        n_poly_swap(Qinv, t);
        n_poly_clear(t);
    }

    Qinv->length = n;
    _n_poly_normalise(Qinv);
}