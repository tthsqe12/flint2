/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mpoly_factor.h"
#include "fmpq_poly.h"
#include "fmpz_mod_mpoly.h"
#include "nmod_mpoly_factor.h"

#define USE_ZAS 1
#define USE_WANG 2
#define USE_ZIP 4


void fmpz_mpoly_convert_perm(
    fmpz_mpoly_t A,
    flint_bitcnt_t Abits,
    const fmpz_mpoly_ctx_t Actx,
    const fmpz_mpoly_t B,
    const fmpz_mpoly_ctx_t Bctx,
    const slong * perm)
{
    slong n = Bctx->minfo->nvars;
    slong m = Actx->minfo->nvars;
    slong i, k, l;
    slong NA, NB;
    ulong * Aexps;
    ulong * Bexps;
    TMP_INIT;

    FLINT_ASSERT(A->bits <= FLINT_BITS);
    FLINT_ASSERT(B->bits <= FLINT_BITS);
    FLINT_ASSERT(B->length > 0);
    TMP_START;

    Aexps = (ulong *) TMP_ALLOC(m*sizeof(ulong));
    Bexps = (ulong *) TMP_ALLOC(n*sizeof(ulong));

    NA = mpoly_words_per_exp(Abits, Actx->minfo);
    NB = mpoly_words_per_exp(B->bits, Bctx->minfo);

    fmpz_mpoly_fit_length_set_bits(A, B->length, Abits, Actx);
    A->length = B->length;
    for (i = 0; i < B->length; i++)
    {        
        fmpz_set(A->coeffs + i, B->coeffs + i);
        mpoly_get_monomial_ui(Bexps, B->exps + NB*i, B->bits, Bctx->minfo);
        for (k = 0; k < m; k++)
        {
            l = perm[k];
            Aexps[k] = l < 0 ? 0 : Bexps[l];
        }
        mpoly_set_monomial_ui(A->exps + NA*i, Aexps, Abits, Actx->minfo);
     }  
    TMP_END;
    fmpz_mpoly_sort_terms(A, Actx);
}


slong _mpoly_compress(
    slong * V,      /* n*n matrix */
    slong * D,      /* n vector */
    slong * deg,    /* n vector */
    slong * S,      /* lxn matrix */
    slong n,
    slong l)
{
    slong mind, maxd;
    slong * minp = FLINT_ARRAY_ALLOC(n, slong);
    slong * maxp = FLINT_ARRAY_ALLOC(n, slong);
    slong * minn = FLINT_ARRAY_ALLOC(n, slong);
    slong * maxn = FLINT_ARRAY_ALLOC(n, slong);
    slong * perm = FLINT_ARRAY_ALLOC(n, slong);
    slong * tmp = maxn;
    slong i, j, k, m;
    slong best_prod, best_loc_i, best_loc_j, best_min, best_deg;

    FLINT_ASSERT(n > 1);
    FLINT_ASSERT(l > 1);

    for (i = 0; i < n; i++)
    for (j = 0; j < n; j++)
        V[i*n + j] = (i == j);

    for (j = 0; j < n; j++)
        minp[j] = maxp[j] = S[0*n + j];
    for (i = 1; i < l; i++)
    for (j = 0; j < n; j++)
    {
        minp[j] = FLINT_MIN(minp[j], S[i*n + j]);
        maxp[j] = FLINT_MAX(maxp[j], S[i*n + j]);
    }

    for (j = 0; j < n; j++)
    {
        D[j] = minp[j];
        deg[j] = 1 + maxp[j] - minp[j];
    }

    for (i = 0; i < l; i++)
    for (j = 0; j < n; j++)
        S[i*n + j] -= D[j];

again:
/*
    flint_printf("S:\n");
    for (k = 0; k < l; k++)
    {
        for (j = 0; j < n; j++)
            flint_printf("%05wd", S[k*n + j]);
        flint_printf("\n");
    }

    flint_printf("deg:\n");
        for (j = 0; j < n; j++)
            flint_printf("%05wd", deg[j] - 1);
        flint_printf("\n");
*/

    best_prod = 1;
    for (j = 0; j < n; j++)
        best_prod *= deg[j];

    best_loc_i = -1;
    best_loc_j = -1;

    for (i = 0; i < n; i++)
    {
        slong this_best_j, this_best_min = WORD_MAX, this_best_deg, this_prod;
/*
flint_printf("i = %wd\n", i);
*/
        mind = WORD_MAX;
        maxd = WORD_MIN;
        for (j = 0; j < n; j++)
        {
            minp[j] = minn[j] = WORD_MAX;
            maxp[j] = maxn[j] = WORD_MIN;
        }
        for (k = 0; k < l; k++)
        {
            slong * Sk = S + k*n;
            slong tot = 0;
            for (j = 0; j < n; j++)
            {
                tot += Sk[j];
                minp[j] = FLINT_MIN(minp[j], Sk[i] + Sk[j]);
                maxp[j] = FLINT_MAX(maxp[j], Sk[i] + Sk[j]);
                minn[j] = FLINT_MIN(minn[j], Sk[i] - Sk[j]);
                maxn[j] = FLINT_MAX(maxn[j], Sk[i] - Sk[j]);
            }
            mind = FLINT_MIN(mind, tot);
            maxd = FLINT_MAX(maxd, tot);
        }
/*
flint_printf("d: %wd - %wd\n", mind, maxd);
for (j = 0; j < n; j++)
{
flint_printf("p[%wd]: %wd - %wd\n", j, minp[j], maxp[j]);
flint_printf("n[%wd]: %wd - %wd\n", j, minn[j], maxn[j]);
}
*/
        this_best_deg = deg[i];
        this_best_j = n + 1; /* something > n */
        if (1 + maxd - mind < this_best_deg)
        {
            this_best_j = 0;
            this_best_min = mind;
            this_best_deg = 1 + maxd - mind;
        }
        for (j = 0; j < n; j++)
        {
            if (j == i)
                continue;
            if (1 + maxp[j] - minp[j] < this_best_deg)
            {
                this_best_j = 1 + j;
                this_best_min = minp[j];
                this_best_deg = 1 + maxp[j] - minp[j];
            }
            if (1 + maxn[j] - minn[j] < this_best_deg)
            {
                this_best_j = -1 - j;
                this_best_min = minn[j];
                this_best_deg = 1 + maxn[j] - minn[j];
            }
        }
/*
flint_printf("this_best_j: %wd\n", this_best_j);
*/
        if (this_best_j > n)
            continue;

        this_prod = this_best_deg;
        for (j = 0; j < n; j++)
            if (j != i)
                this_prod *= deg[j];
        if (this_prod < best_prod)
        {
            best_prod = this_prod;
            best_loc_i = i;
            best_loc_j = this_best_j;
            best_min = this_best_min;
            best_deg = this_best_deg;
        }
    }
/*
flint_printf("best_loc: %wd, %wd\n", best_loc_i, best_loc_j);
*/

    if (best_loc_i >= 0)
    {
        i = best_loc_i;
        j = best_loc_j;
        deg[i] = best_deg;
        if (j < 0)
        {
            j = -j - 1;
/*
flint_printf("x%wd -= x%wd\n", i+1, j+1);
*/
            for (k = 0; k < l; k++)
                S[k*n + i] += -S[k*n + j] - best_min;
            for (k = 0; k < n; k++)
            {
                D[k] += best_min*V[k*n + i];
                V[k*n + j] += V[k*n + i];
            }
        }
        else if (j > 0)
        {
            j = j - 1;
/*
flint_printf("x%wd += x%wd\n", i+1, j+1);
*/
            for (k = 0; k < l; k++)
                S[k*n + i] += S[k*n + j] - best_min;
            for (k = 0; k < n; k++)
            {
                D[k] += best_min*V[k*n + i];
                V[k*n + j] += -V[k*n + i];
            }         
        }
        else
        {
/*
flint_printf("x%wd = total\n", i+1);
*/
            for (k = 0; k < l; k++)
            {
                slong tot = 0;
                for (j = 0; j < n; j++)
                    tot += S[k*n + j];
                S[k*n + i] = tot - best_min;
            }
            for (k = 0; k < n; k++)
            {
                D[k] += best_min*V[k*n + i];
                for (j = 0; j < n; j++)
                    if (j != i)
                        V[k*n + j] += -V[k*n + i];
            }
        }

        goto again;
    }

    for (i = 0; i < n; i++)
        tmp[i] = i;
    for (i = 1; i < n; i++)
        for (j = i; j > 0 && deg[tmp[j]] < deg[tmp[j - 1]]; j--)
            SLONG_SWAP(tmp[j], tmp[j - 1]);
    m = 1;
    while (m < n && deg[tmp[n - (m + 1)]] > 1)
        m++;
    for (i = 0; i < n; i++)
        perm[i] = (i < m) ? tmp[n - m + i] : tmp[i - m];

    for (i = 0; i < n; i++)
        tmp[i] = deg[perm[i]] - 1;
    for (i = 0; i < n; i++)
        deg[i] = tmp[i];

    for (k = 0; k < l; k++)
    {
        slong * Sk = S + k*n;
        for (i = 0; i < n; i++)
            tmp[i] = Sk[perm[i]];
        for (i = 0; i < n; i++)
            Sk[i] = tmp[i];
    }

    for (k = 0; k < n; k++)
    {
        slong * Vk = V + k*n;
        for (i = 0; i < n; i++)
            tmp[i] = Vk[perm[i]];
        for (i = 0; i < n; i++)
            Vk[i] = tmp[i];
    }

    flint_free(minp);
    flint_free(maxp);
    flint_free(minn);
    flint_free(maxn);
    flint_free(perm);
/*
    flint_printf("S:\n");
    for (k = 0; k < l; k++)
    {
        for (j = 0; j < n; j++)
            flint_printf("%05wd", S[k*n + j]);
        flint_printf("\n");
    }

    flint_printf("V:\n");
    for (k = 0; k < n; k++)
    {
        for (j = 0; j < n; j++)
            flint_printf("%05wd", V[k*n + j]);
        flint_printf("\n");
    }

    flint_printf("delta:\n");
        for (j = 0; j < n; j++)
            flint_printf("%05wd", D[j]);
        flint_printf("\n");

    flint_printf("deg:\n");
        for (j = 0; j < n; j++)
            flint_printf("%05wd", deg[j]);
        flint_printf("\n");
*/
    return m;
}


typedef struct {
    fmpz_mpoly_ctx_t ctx;
    slong nvars;
    slong exps_alloc;
    slong * exps;
    slong * umat;
    slong * deltas;
    slong * degs;
    int is_id;
} fmpz_mpoly_perm_struct;

typedef fmpz_mpoly_perm_struct fmpz_mpoly_perm_t[1];



/*
    figure out an affine transformation
    if new nvars = old nvars and transformation is id
        swap(A, B)
        return
    end
    map 
    see if content needs to removed wrt new var0
*/
void fmpz_mpoly_perm_init(
    fmpz_mpoly_perm_t P,
    fmpz_mpoly_t L,
    fmpz_mpoly_t A,
    const fmpz_mpoly_ctx_t Actx)
{
    int same;
    slong i, j, max_deg;
    flint_bitcnt_t Lbits, Abits = A->bits;
    slong LN, AN = mpoly_words_per_exp_sp(Abits, Actx->minfo);
    slong Alen = A->length;
    slong mvars, nvars = Actx->minfo->nvars;

    P->nvars = nvars;
    P->umat = FLINT_ARRAY_ALLOC(nvars*nvars, slong);
    P->deltas = FLINT_ARRAY_ALLOC(nvars, slong);
    P->degs = FLINT_ARRAY_ALLOC(nvars, slong);
    P->exps = FLINT_ARRAY_ALLOC(Alen*nvars, slong);
    P->exps_alloc = Alen*nvars;
    for (i = 0; i < Alen; i++)
        mpoly_get_monomial_ui_sp((ulong *)P->exps + nvars*i,
                                           A->exps + AN*i, Abits, Actx->minfo);
    mvars = _mpoly_compress(P->umat, P->deltas, P->degs, P->exps, nvars, Alen);

    FLINT_ASSERT(mvars > 0);

    fmpz_mpoly_ctx_init(P->ctx, mvars, ORD_LEX);
    fmpz_mpoly_init(L, P->ctx);

    same = (mvars == nvars) && (Actx->minfo->ord == ORD_LEX);
    if (same)
    {
        for (i = 0; i < nvars; i++)
        {
            same = same && (P->deltas[i] == 0);
            for (j = 0; j < nvars; j++)
                same = same && (P->umat[i*nvars + j] == (i == j));
        }
    }

    P->is_id = same;
    if (same)
    {
        fmpz_mpoly_swap(L, A, Actx);
        return;
    }

    max_deg = P->degs[0];
    for (i = 1; i < mvars; i++)
        max_deg = FLINT_MAX(max_deg, P->degs[i]);
    Lbits = mpoly_fix_bits(1 + FLINT_BIT_COUNT(max_deg), P->ctx->minfo);

    fmpz_mpoly_fit_length_set_bits(L, Alen, Lbits, P->ctx);

    LN = mpoly_words_per_exp(Lbits, P->ctx->minfo);

    L->length = Alen;
    for (i = 0; i < Alen; i++)
    {
        fmpz_swap(L->coeffs + i, A->coeffs + i);
        mpoly_set_monomial_ui(L->exps + LN*i,
                             (ulong *)P->exps + nvars*i, Lbits, P->ctx->minfo);
    }

    fmpz_mpoly_sort_terms(L, P->ctx);
    fmpz_mpoly_unit_normalize(L, P->ctx);
}


void slong_array_fit_length(slong ** array, slong * alloc, slong len)
{
    if (len <= *alloc)
        return;
    len = FLINT_MAX(len, *alloc + *alloc/2 + 1);
    *array = flint_realloc(*array, len*sizeof(slong));
    *alloc = len;
}


void fmpz_mpoly_perm_expand_fmpz_poly(
    fmpz_mpoly_t A,
    flint_bitcnt_t Abits,
    const fmpz_mpoly_ctx_t Actx,
    fmpz_poly_t B,
    slong var,
    fmpz_mpoly_perm_t P)
{
    slong i, j, Alen;
    slong NA, nvars = Actx->minfo->nvars;
    slong * mins = FLINT_ARRAY_ALLOC(nvars, slong);
/*
flint_printf("fmpz_mpoly_perm_expand_fmpz_poly called\n");
flint_printf("B: ");
fmpz_poly_print_pretty(B, "X");
flint_printf("\n");
*/

    for (j = 0; j < nvars; j++)
        mins[j] = WORD_MAX;

    fmpz_mpoly_fit_length_set_bits(A, B->length, Abits, Actx);

    Alen = 0;
    for (i = B->length - 1; i >= 0; i--)
    {
        if (fmpz_is_zero(B->coeffs + i))
            continue;

        fmpz_swap(A->coeffs + Alen, B->coeffs + i);

        slong_array_fit_length(&P->exps, &P->exps_alloc, (Alen + 1)*nvars);

        for (j = 0; j < nvars; j++)
        {
            P->exps[Alen*nvars + j] = P->umat[j*nvars + var]*i +
                                      P->deltas[j];
            mins[j] = FLINT_MIN(mins[j], P->exps[Alen*nvars + j]);
        }
        Alen++;
    }

    A->length = Alen;
    NA = mpoly_words_per_exp(Abits, Actx->minfo);

    for (i = 0; i < Alen; i++)
    {
        for (j = 0; j < nvars; j++)
            P->exps[i*nvars + j] -= mins[j];

        mpoly_set_monomial_ui(A->exps + NA*i, (ulong *)P->exps + i*nvars, Abits, Actx->minfo);
    }

    fmpz_mpoly_sort_terms(A, Actx);
    fmpz_mpoly_unit_normalize(A, Actx);
/*
flint_printf("A: ");
fmpz_mpoly_print_pretty(A, NULL, Actx);
flint_printf("\n");
*/
}


void fmpz_mpoly_perm_expand_fmpz_bpoly(
    fmpz_mpoly_t A,
    flint_bitcnt_t Abits,
    const fmpz_mpoly_ctx_t Actx,
    fmpz_bpoly_t B,
    slong var0,
    slong var1,
    fmpz_mpoly_perm_t P)
{
    slong i, j, k, Ai, Alen;
    slong NA, nvars = Actx->minfo->nvars;
    slong * mins = FLINT_ARRAY_ALLOC(nvars, slong);
/*
flint_printf("fmpz_mpoly_perm_expand_fmpz_bpoly called\n");
flint_printf("B: ");
fmpz_bpoly_print_pretty(B, "X", "Y");
flint_printf("\n");
*/

    for (k = 0; k < nvars; k++)
        mins[k] = WORD_MAX;

    Alen = 0;
    for (i = B->length - 1; i >= 0; i--)
    {
        fmpz_poly_struct * Bi = B->coeffs + i;
        for (j = Bi->length - 1; j >= 0; j--)
            Alen += !fmpz_is_zero(Bi->coeffs + j);
    }

    fmpz_mpoly_fit_length_set_bits(A, Alen + 4, Abits, Actx);
    slong_array_fit_length(&P->exps, &P->exps_alloc, Alen*nvars);

    Ai = 0;
    for (i = B->length - 1; i >= 0; i--)
    {
        fmpz_poly_struct * Bi = B->coeffs + i;
        for (j = Bi->length - 1; j >= 0; j--)
        {
            if (fmpz_is_zero(Bi->coeffs + j))
                continue;
            fmpz_swap(A->coeffs + Ai, Bi->coeffs + j);
            for (k = 0; k < nvars; k++)
            {
                P->exps[Ai*nvars + k] = P->umat[k*nvars + var0]*i +
                                        P->umat[k*nvars + var1]*j +
                                        P->deltas[k];
                mins[k] = FLINT_MIN(mins[k], P->exps[Ai*nvars + k]);
            }
            Ai++;
        }
    }

    FLINT_ASSERT(Ai == Alen);

    A->length = Alen;
    NA = mpoly_words_per_exp(Abits, Actx->minfo);

    for (i = 0; i < Alen; i++)
    {
        for (k = 0; k < nvars; k++)
            P->exps[i*nvars + k] -= mins[k];

        mpoly_set_monomial_ui(A->exps + NA*i, (ulong *)P->exps + i*nvars, Abits, Actx->minfo);
    }

    fmpz_mpoly_sort_terms(A, Actx);
    fmpz_mpoly_unit_normalize(A, Actx);
/*
flint_printf("A: ");
fmpz_mpoly_print_pretty(A, NULL, Actx);
flint_printf("\n");
*/
}

void fmpz_mpoly_perm_expand_fmpz_mpoly(
    fmpz_mpoly_t A,
    flint_bitcnt_t Abits,
    const fmpz_mpoly_ctx_t Actx,
    fmpz_mpoly_t B,
    fmpz_mpoly_perm_t P)
{
    slong i, k, l;
    slong nvars = Actx->minfo->nvars;
    slong NA = mpoly_words_per_exp(Abits, Actx->minfo);
    slong mvars = P->ctx->minfo->nvars;
    flint_bitcnt_t Bbits = B->bits;
    slong NB = mpoly_words_per_exp(Bbits, P->ctx->minfo);
    slong * mins, * texps;
/*
flint_printf("fmpz_mpoly_perm_expand_fmpz_bpoly called\n");
flint_printf("B: ");
fmpz_bpoly_print_pretty(B, "X", "Y");
flint_printf("\n");
*/
    FLINT_ASSERT(fmpz_mpoly_degrees_fit_si(B, P->ctx));

    if (P->is_id)
    {
        fmpz_mpoly_swap(A, B, Actx);
        return;
    }

    texps = FLINT_ARRAY_ALLOC(nvars, slong);
    mins = FLINT_ARRAY_ALLOC(nvars, slong);
    for (k = 0; k < nvars; k++)
        mins[k] = WORD_MAX;

    slong_array_fit_length(&P->exps, &P->exps_alloc, B->length*nvars);
    fmpz_mpoly_fit_length_set_bits(A, B->length, Abits, Actx);
    _fmpz_mpoly_set_length(A, B->length, Actx);
    for (i = 0; i < B->length; i++)
    {
        fmpz_swap(A->coeffs + i, B->coeffs + i);
        mpoly_get_monomial_ui((ulong *)texps, B->exps + NB*i, Bbits, P->ctx->minfo);
        for (k = 0; k < nvars; k++)
        {
            slong tot = P->deltas[k];
            for (l = 0; l < mvars; l++)
                tot += P->umat[k*nvars + l]*texps[l];
            P->exps[i*nvars + k] = tot;
            mins[k] = FLINT_MIN(mins[k], tot);
        }
    }

    for (i = 0; i < B->length; i++)
    {
        for (k = 0; k < nvars; k++)
            P->exps[i*nvars + k] -= mins[k];
        mpoly_set_monomial_ui(A->exps + NA*i, (ulong *)P->exps + i*nvars, Abits, Actx->minfo);
    }

    fmpz_mpoly_sort_terms(A, Actx);
    fmpz_mpoly_unit_normalize(A, Actx);
/*
flint_printf("A: ");
fmpz_mpoly_print_pretty(A, NULL, Actx);
flint_printf("\n");
*/
}

/*
    A is primitive w.r.t to any variable appearing in A.
    A is square free with positive lead coeff.

    return 1 for success, 0 for failure
*/
static int _irreducible_factors(
    fmpz_mpolyv_t Af,
    fmpz_mpoly_t A,
    const fmpz_mpoly_ctx_t Actx,
    unsigned int algo)
{
    int success;
    slong i;
    flint_bitcnt_t Abits;
    double density;
    fmpz_mpoly_perm_t P;
    fmpz_mpoly_t L;
    fmpz_mpoly_ctx_struct *  Lctx;
    fmpz_poly_t u;
    fmpz_poly_factor_t uf;
    flint_rand_t state;
#if WANT_ASSERT
    fmpz_mpoly_t Aorg;

    fmpz_mpoly_init(Aorg, Actx);
    fmpz_mpoly_set(Aorg, A, Actx);
#endif

    FLINT_ASSERT(A->length == 0 || fmpz_sgn(A->coeffs + 0) > 0);

    if (A->length < 2)
    {
        FLINT_ASSERT(A->length == 1);
        FLINT_ASSERT(!fmpz_mpoly_is_fmpz(A, Actx));
        fmpz_mpolyv_fit_length(Af, 1, Actx);
        Af->length = 1;
        fmpz_mpoly_swap(Af->coeffs + 0, A, Actx);
        success = 1;
        goto cleanup_less;
    }

    if (!fmpz_mpoly_degrees_fit_si(A, Actx))
    {
        success = 0;
        goto cleanup_less;
    }

    if (A->bits > FLINT_BITS &&
        !fmpz_mpoly_repack_bits_inplace(A, FLINT_BITS, Actx))
    {
        success = 0;
        goto cleanup_less;
    }

    Abits = A->bits;

    flint_randinit(state);
    fmpz_poly_init(u);
    fmpz_poly_factor_init(uf);
    fmpz_mpoly_perm_init(P, L, A, Actx);
    Lctx = P->ctx;
/*
flint_printf("L: ");
fmpz_mpoly_print_pretty(L, NULL, Lctx);
flint_printf("\n");
*/
    density = L->length;
    for (i = 0; i < P->ctx->minfo->nvars; i++)
        density /= P->degs[i] + 1;

    if (!P->is_id)
    {
        fmpz_mpoly_univar_t U;
        fmpz_mpoly_t t;

flint_printf("it wasn't identity\n");

        fmpz_mpoly_init(t, Lctx);
        fmpz_mpoly_univar_init(U, Lctx);

        fmpz_mpoly_to_univar(U, L, 0, Lctx);
        success = _fmpz_mpoly_vec_content_mpoly(t, U->coeffs, U->length, Lctx);
        FLINT_ASSERT(success);
        FLINT_ASSERT(fmpz_mpoly_is_fmpz(t, Lctx));
        success = fmpz_mpoly_divides(L, L, t, Lctx);
        FLINT_ASSERT(success);
        fmpz_mpoly_unit_normalize(L, Lctx);

        fmpz_mpoly_clear(t, Lctx);
        fmpz_mpoly_univar_clear(U, Lctx);
    }

    if (P->degs[0] == 1)
    {
        fmpz_mpolyv_fit_length(Af, 1, Actx);
        Af->length = 1;
        fmpz_mpoly_perm_expand_fmpz_mpoly(Af->coeffs + 0, Abits, Actx, L, P);
        success = 1;
        goto cleanup_more;
    }
    else if (P->degs[0] == 2)
    {
        fmpz_mpoly_univar_t U;
        fmpz_mpoly_t aa, bb, cc, tt0, tt1, tt2;

flint_printf("hey quadratic\n");

        fmpz_mpoly_univar_init(U, Lctx);
        fmpz_mpoly_init(aa, Lctx);
        fmpz_mpoly_init(bb, Lctx);
        fmpz_mpoly_init(cc, Lctx);
        fmpz_mpoly_init(tt0, Lctx);
        fmpz_mpoly_init(tt1, Lctx);
        fmpz_mpoly_init(tt2, Lctx);

        fmpz_mpoly_to_univar(U, L, 0, Lctx);

flint_printf("U: ");
fmpz_mpoly_univar_print_pretty(U, NULL, Lctx);
flint_printf("\n");

        for (i = 0; i < U->length; i++)
        {
            if (fmpz_cmp_ui(U->exps + i, 2) == 0)
                fmpz_mpoly_swap(aa, U->coeffs + i, Lctx);
            else if (fmpz_cmp_ui(U->exps + i, 1) == 0)
                fmpz_mpoly_swap(bb, U->coeffs + i, Lctx);
            else if (fmpz_cmp_ui(U->exps + i, 0) == 0)
                fmpz_mpoly_swap(cc, U->coeffs + i, Lctx);
        }

flint_printf("aa: ");
fmpz_mpoly_print_pretty(aa, NULL, Lctx);
flint_printf("\n");

flint_printf("bb: ");
fmpz_mpoly_print_pretty(bb, NULL, Lctx);
flint_printf("\n");

flint_printf("cc: ");
fmpz_mpoly_print_pretty(cc, NULL, Lctx);
flint_printf("\n");

        fmpz_mpoly_mul(tt0, bb, bb, Lctx);
        fmpz_mpoly_mul(tt1, aa, cc, Lctx);
        fmpz_mpoly_scalar_mul_si(tt1, tt1, 4, Lctx);
        fmpz_mpoly_sub(tt2, tt0, tt1, Lctx);
        if (!fmpz_mpoly_sqrt(tt0, tt2, Lctx))
        {
            flint_printf("not a square");
            flint_printf("\n");
        }
        fmpz_mpoly_add(tt2, tt0, bb, Lctx);
        if (!fmpz_mpoly_scalar_divides_si(tt2, tt2, 2, Lctx))
        {
            flint_printf("should not happen\n");
            flint_abort();
        }

        fmpz_mpoly_gcd_cofactors(tt0, tt1, tt2, aa, tt2, Lctx);
        if (!fmpz_mpoly_divides(bb, cc, tt2, Lctx))
        {
            flint_printf("should not happen\n");
            flint_abort();
        }

        _mpoly_gen_shift_left(tt1->exps, tt1->bits, tt1->length, 0, 1, Lctx->minfo);
        fmpz_mpoly_add(tt1, tt1, tt2, Lctx);

        _mpoly_gen_shift_left(tt0->exps, tt0->bits, tt0->length, 0, 1, Lctx->minfo);
        fmpz_mpoly_add(tt0, tt0, bb, Lctx);

flint_printf("tt0: ");
fmpz_mpoly_print_pretty(tt0, NULL, Lctx);
flint_printf("\n");

flint_printf("tt1: ");
fmpz_mpoly_print_pretty(tt1, NULL, Lctx);
flint_printf("\n");

        fmpz_mpolyv_fit_length(Af, 2, Actx);
        Af->length = 2;
        fmpz_mpoly_perm_expand_fmpz_mpoly(Af->coeffs + 0, Abits, Actx, tt0, P);

flint_printf("Af[0]: ");
fmpz_mpoly_print_pretty(Af->coeffs + 0, NULL, Actx);
flint_printf("\n");

        fmpz_mpoly_perm_expand_fmpz_mpoly(Af->coeffs + 1, Abits, Actx, tt1, P);

flint_printf("Af[1]: ");
fmpz_mpoly_print_pretty(Af->coeffs + 1, NULL, Actx);
flint_printf("\n");


        success = 1;
        goto cleanup_less;
    }

    if (Lctx->minfo->nvars < 2)
    {
        FLINT_ASSERT(fmpz_mpoly_is_fmpz_poly(L, 0, Lctx));
        success = fmpz_mpoly_get_fmpz_poly(u, L, 0, Lctx);
        FLINT_ASSERT(success);
        fmpz_poly_factor(uf, u);

        FLINT_ASSERT(fmpz_is_pm1(&uf->c));
        for (i = 0; i < uf->num; i++)
            FLINT_ASSERT(uf->exp[i] == 1);

        fmpz_mpolyv_fit_length(Af, uf->num, Actx);
        Af->length = uf->num;
        for (i = 0; i < uf->num; i++)
        {
            fmpz_mpoly_perm_expand_fmpz_poly(Af->coeffs + i, Abits, Actx,
                                                              uf->p + i, 0, P);
        }

        success = 1;
    }
    else if (Lctx->minfo->nvars == 2)
    {
        fmpz_bpoly_t b;
        fmpz_tpoly_t bf;

        fmpz_bpoly_init(b);
        fmpz_tpoly_init(bf);

        fmpz_mpoly_get_bpoly(b, L, 0, 1, Lctx);
        fmpz_bpoly_factor(u, bf, b);
        fmpz_poly_factor(uf, u);

        fmpz_mpolyv_fit_length(Af, bf->length + uf->num, Actx);
        Af->length = 0;
        for (i = 0; i < bf->length; i++)
            fmpz_mpoly_perm_expand_fmpz_bpoly(Af->coeffs + (Af->length++),
                                         Abits, Actx, bf->coeffs + i, 0, 1, P);
        for (i = 0; i < uf->num; i++)
            fmpz_mpoly_perm_expand_fmpz_poly(Af->coeffs + (Af->length++),
                                                 Abits, Actx, uf->p + i, 1, P);
        fmpz_bpoly_clear(b);
        fmpz_tpoly_clear(bf);

        success = 1;
    }
    else
    {
        fmpz_mpoly_t lcL;
        fmpz_mpolyv_t Lf;
        fmpz_mpoly_factor_t lcLf;
        zassenhaus_prune_t Z;
        fmpz * alpha;
        int zero_ok, trying_zero, image_count, sqrfree;
        ulong alpha_modulus;
        slong main_degree = P->degs[0];

        fmpz_mpoly_init(lcL, Lctx);
        fmpz_mpolyv_init(Lf, Lctx);
        fmpz_mpoly_factor_init(lcLf, Lctx);
        zassenhaus_prune_init(Z);

        zassenhaus_prune_set_degree(Z, main_degree);

        /* some simple checks */

        alpha = _fmpz_vec_init(Lctx->minfo->nvars - 1);
        zero_ok = 0;
        trying_zero = 1;
        alpha_modulus = 5;
        image_count = 0;

        goto got_alpha;

next_alpha:

        trying_zero = 0;

        if (++alpha_modulus > 10)
            goto done_alpha;

        for (i = 0; i < Lctx->minfo->nvars - 1; i++)
        {
            slong a = n_urandint(state, alpha_modulus);
            a -= alpha_modulus/2;
            fmpz_set_si(alpha + i, a);
        }

got_alpha:

        _fmpz_mpoly_eval_rest_to_poly(u, L, alpha, Lctx);

        if (fmpz_poly_degree(u) != main_degree)
            goto next_alpha;

        fmpz_poly_factor(uf, u);

        zassenhaus_prune_start_add_factors(Z);                
        sqrfree = 1;
        for (i = 0; i < uf->num; i++)
        {
            if (uf->exp[i] != 1)
                sqrfree = 0;
            zassenhaus_prune_add_factor(Z, fmpz_poly_degree(uf->p + i), uf->exp[i]);
        }
        zassenhaus_prune_end_add_factors(Z);

        if (!sqrfree)
            goto next_alpha;

        zero_ok = zero_ok || trying_zero;
        if (++image_count < 3)
            goto next_alpha;

done_alpha:

        _fmpz_vec_clear(alpha, Lctx->minfo->nvars - 1);

        /* simple check done */

        if (zassenhaus_prune_must_be_irreducible(Z))
        {
            fmpz_mpolyv_fit_length(Af, 1, Actx);
            Af->length = 1;
            fmpz_mpoly_perm_expand_fmpz_mpoly(Af->coeffs + 0, Abits, Actx, L, P);
            success = 1;
            goto cleanup_more;
        }

        success = 0;

        if (algo & (USE_WANG | USE_ZIP))
        {
            _fmpz_mpoly_get_lead0(lcL, L, Lctx);

            if (fmpz_mpoly_factor_squarefree(lcLf, lcL, Lctx))
            {
                int irr_fac = 1;
                for (i = 0; i < lcLf->num; i++)
                    irr_fac = irr_fac && lcLf->poly[i].length < 4;

                irr_fac = irr_fac && fmpz_mpoly_factor_irred(lcLf, Lctx, algo);

                if (!(algo & USE_ZIP))
                {
                    success = fmpz_mpoly_factor_irred_wang(Lf, L,
                                        lcLf, irr_fac, lcL, Lctx, state, Z, 1);
                }
                else if (!(algo & USE_WANG))
                {
                    success = fmpz_mpoly_factor_irred_zippel(Lf, L,
                                           lcLf, irr_fac, lcL, Lctx, state, Z);
                }
                else
                {
                    if (density > 0.002 && zero_ok)
                    {
                        success = fmpz_mpoly_factor_irred_wang(Lf, L,
                                        lcLf, irr_fac, lcL, Lctx, state, Z, 0);
                    }

                    if (success == 0 && density > 0.04)
                    {
                        success = fmpz_mpoly_factor_irred_wang(Lf, L,
                                        lcLf, irr_fac, lcL, Lctx, state, Z, 1);
                    }

                    if (success == 0)
                    {
                        success = fmpz_mpoly_factor_irred_zippel(Lf, L,
                                           lcLf, irr_fac, lcL, Lctx, state, Z);
                    }

                    if (success == 0)
                    {
                        success = fmpz_mpoly_factor_irred_wang(Lf, L,
                                        lcLf, irr_fac, lcL, Lctx, state, Z, 1);
                    }
                }
            }
        }

        if (algo & USE_ZAS)
        {
            if (success == 0)
                success = fmpz_mpoly_factor_irred_zassenhaus(Lf, L, Lctx, Z);
        }

        success = (success > 0);
        if (success)
        {
            fmpz_mpolyv_fit_length(Af, Lf->length, Actx);
            Af->length = Lf->length;
            for (i = 0; i < Lf->length; i++)
                fmpz_mpoly_perm_expand_fmpz_mpoly(Af->coeffs + i, Abits, Actx,
                                                            Lf->coeffs + i, P);
        }

cleanup_more:

        fmpz_mpoly_clear(lcL, Lctx);
        fmpz_mpolyv_clear(Lf, Lctx);
        fmpz_mpoly_factor_clear(lcLf, Lctx);
        zassenhaus_prune_clear(Z);
    }

    fmpz_poly_clear(u);
    fmpz_poly_factor_clear(uf);
    flint_randclear(state);

cleanup_less:

#if WANT_ASSERT
    if (success)
    {
        fmpz_mpoly_t prod;
        fmpz_mpoly_init(prod, Actx);
        fmpz_mpoly_one(prod, Actx);
        for (i = 0; i < Af->length; i++)
            fmpz_mpoly_mul(prod, prod, Af->coeffs + i, Actx);
        FLINT_ASSERT(fmpz_mpoly_equal(prod, Aorg, Actx));
        fmpz_mpoly_clear(prod, Actx);
        fmpz_mpoly_clear(Aorg, Actx);
    }
#endif

    FLINT_ASSERT(success == 0 || success == 1);
    return success;
}


/* assume f is a square free factorization, replace it by an irreducible one */
int fmpz_mpoly_factor_irred(
    fmpz_mpoly_factor_t f,
    const fmpz_mpoly_ctx_t ctx,
    unsigned int algo)
{
    int success;
    slong i, j;
    fmpz_mpolyv_t t;
    fmpz_mpoly_factor_t g;

    fmpz_mpolyv_init(t, ctx);
    fmpz_mpoly_factor_init(g, ctx);

    fmpz_swap(g->constant, f->constant);
    g->num = 0;
    for (j = 0; j < f->num; j++)
    {
        success = _irreducible_factors(t, f->poly + j, ctx, algo);
        if (!success)
            goto cleanup;

        fmpz_mpoly_factor_fit_length(g, g->num + t->length, ctx);
        for (i = 0; i < t->length; i++)
        {
            fmpz_set(g->exp + g->num, f->exp + j);
            fmpz_mpoly_swap(g->poly + g->num, t->coeffs + i, ctx);
            g->num++;
        }
    }
    fmpz_mpoly_factor_swap(f, g, ctx);

    success = 1;

cleanup:

    fmpz_mpolyv_clear(t, ctx);
    fmpz_mpoly_factor_clear(g, ctx);

    return success;
}


int fmpz_mpoly_factor(
    fmpz_mpoly_factor_t f,
    const fmpz_mpoly_t A,
    const fmpz_mpoly_ctx_t ctx)
{
    int success;
    success = fmpz_mpoly_factor_squarefree(f, A, ctx) &&
              fmpz_mpoly_factor_irred(f, ctx, USE_ZAS | USE_WANG | USE_ZIP);
    FLINT_ASSERT(!success || fmpz_mpoly_factor_matches(A, f, ctx));
    return success;
}


int fmpz_mpoly_factor_zassenhaus(
    fmpz_mpoly_factor_t f,
    const fmpz_mpoly_t A,
    const fmpz_mpoly_ctx_t ctx)
{
    int success;
    success = fmpz_mpoly_factor_squarefree(f, A, ctx) &&
              fmpz_mpoly_factor_irred(f, ctx, USE_ZAS);
    FLINT_ASSERT(!success || fmpz_mpoly_factor_matches(A, f, ctx));
    return success;
}


int fmpz_mpoly_factor_wang(
    fmpz_mpoly_factor_t f,
    const fmpz_mpoly_t A,
    const fmpz_mpoly_ctx_t ctx)
{
    int success;
    success = fmpz_mpoly_factor_squarefree(f, A, ctx) &&
              fmpz_mpoly_factor_irred(f, ctx, USE_WANG);
    FLINT_ASSERT(!success || fmpz_mpoly_factor_matches(A, f, ctx));
    return success;
}


int fmpz_mpoly_factor_zippel(
    fmpz_mpoly_factor_t f,
    const fmpz_mpoly_t A,
    const fmpz_mpoly_ctx_t ctx)
{
    int success;
    success = fmpz_mpoly_factor_squarefree(f, A, ctx) &&
              fmpz_mpoly_factor_irred(f, ctx, USE_ZIP);
    FLINT_ASSERT(!success || fmpz_mpoly_factor_matches(A, f, ctx));
    return success;
}

