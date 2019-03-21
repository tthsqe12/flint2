/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fmpz_mpoly.h"
#include "nmod_mpoly.h"
#include "thread_pool.h"
#include "profiler.h"
int usleep(ulong usec);
void nmod_mpolyun_scalar_mul_nmod(nmod_mpolyun_t A, mp_limb_t c, const nmod_mpoly_ctx_t ctx);
mp_limb_t nmod_mpolyun_leadcoeff_last(nmod_mpolyun_t A, const nmod_mpoly_ctx_t ctx);
int nmod_mpolyun_is_canonical(const nmod_mpolyun_t A, const nmod_mpoly_ctx_t ctx);
int nmod_mpolyun_is_nonzero_nmod(const nmod_mpolyun_t A, const nmod_mpoly_ctx_t ctx);
int nmod_mpolyun_gcd_brown_smprime(nmod_mpolyun_t G,
                                nmod_mpolyun_t Abar, nmod_mpolyun_t Bbar,
                                    nmod_mpolyun_t A, nmod_mpolyun_t B,
                                         slong var, const nmod_mpoly_ctx_t ctx);


pthread_mutex_t iomutex;






/* instructions do A = B + I*(C - B) mod M */
typedef struct
{
    slong a_idx; /* index of A */
    slong b_idx; /* index of B */
    slong c_idx; /* index of C */
    fmpz_t idem;     /* I */
    fmpz_t modulus;  /* M */
} _fmpz_crt_prog_instr;

typedef struct
{
    _fmpz_crt_prog_instr * prog; /* straight line program */
    slong length; /* length of prog */
    slong alloc;  /* alloc of prog */
    slong localsize; /* length of outputs required in nmod_poly_crt_run */
    slong temp1loc; /* index of temporary used in run */
    slong temp2loc; /* index of another tempory used in run */
    int good;   /* the moduli are good for CRT, essentially relatively prime */
} fmpz_crt_struct;

typedef fmpz_crt_struct fmpz_crt_t[1];

/* general crt for nmod_poly_t - compile once, run many times ****************/

void fmpz_crt_init(fmpz_crt_t P);

int fmpz_crt_compile(fmpz_crt_t P, fmpz ** moduli, slong len);

void fmpz_crt_run(const fmpz_crt_t P, fmpz_t output, fmpz ** inputs);


void fmpz_crt_clear(fmpz_crt_t P);

void fmpz_crt_init(fmpz_crt_t P)
{
    P->prog = NULL;
    P->alloc = 0;
    P->length = 0;
    P->localsize = 1;
    P->temp1loc = 0;
    P->temp2loc = 0;
    P->good = 0;
}

static void fmpz_crt_fit_length(fmpz_crt_t P, slong k)
{
    k = FLINT_MAX(WORD(1), k);

    if (P->alloc == 0)
    {
        FLINT_ASSERT(P->prog == NULL);
        P->prog = (_fmpz_crt_prog_instr *) flint_malloc(k
                                           *sizeof(_fmpz_crt_prog_instr));
        P->alloc = k;
    }
    else if (k > P->alloc)
    {
        FLINT_ASSERT(P->prog != NULL);
        P->prog = (_fmpz_crt_prog_instr *) flint_realloc(P->prog, k
                                           *sizeof(_fmpz_crt_prog_instr));
        P->alloc = k;
    }
}

static void fmpz_crt_set_length(fmpz_crt_t P, slong k)
{
    slong i;

    FLINT_ASSERT(k <= P->length);

    for (i = k; i < P->length; i++)
    {
        fmpz_clear(P->prog[i].modulus);
        fmpz_clear(P->prog[i].idem);
    }
    P->length = k;
}

void fmpz_crt_clear(fmpz_crt_t P)
{
    fmpz_crt_set_length(P, 0);

    if (P->alloc > 0)
    {
        flint_free(P->prog);
    }
}

/*
    combine all moduli in [start, stop)
    return index of instruction that computes the result
*/
static slong _push_prog(fmpz_crt_t P,
                           fmpz ** moduli, slong * perm,
                                        slong ret_idx, slong start, slong stop)
{
    slong i, mid;
    slong b_idx, c_idx;
    mp_bitcnt_t lefttot, righttot;
    slong leftret, rightret;
    fmpz * leftmodulus, * rightmodulus;

    /* we should have at least 2 moduli */
    FLINT_ASSERT(start + 1 < stop);

    mid = start + (stop - start)/2;

    FLINT_ASSERT(start < mid);
    FLINT_ASSERT(mid < stop);

    lefttot = 0;
    for (i = start; i < mid; i++)
    {
        lefttot += fmpz_bits(moduli[perm[i]]);
    }

    righttot = 0;
    for (i = mid; i < stop; i++)
    {
        righttot += fmpz_bits(moduli[perm[i]]);
    }

    /* try to balance the total degree on left and right */
    while (lefttot < righttot
            && mid + 1 < stop
            && fmpz_bits(moduli[perm[mid]]) < righttot - lefttot)
    {
        lefttot += fmpz_bits(moduli[perm[mid]]);
        righttot -= fmpz_bits(moduli[perm[mid]]);
        mid++;
    }

    P->localsize = FLINT_MAX(P->localsize, 1 + ret_idx);

    /* compile left [start, mid) */
    if (start + 1 < mid)
    {
        b_idx = ret_idx + 1;
        leftret = _push_prog(P, moduli, perm, b_idx, start, mid);
        if (!P->good)
        {
            return -1;
        }
        leftmodulus = P->prog[leftret].modulus;
    }
    else
    {
        b_idx = -1 - perm[start];
        leftmodulus = moduli[perm[start]];
    }

    /* compile right [mid, end) */
    if (mid + 1 < stop)
    {
        c_idx = ret_idx + 2;
        rightret = _push_prog(P, moduli, perm, c_idx, mid, stop);
        if (!P->good)
        {
            return -1;
        }
        rightmodulus = P->prog[rightret].modulus;
    }
    else
    {
        c_idx = -1 - perm[mid];
        rightmodulus = moduli[perm[mid]];
    }

    /* check if fmpz_invmod is going to throw */
    if (fmpz_is_zero(leftmodulus) || fmpz_is_zero(rightmodulus))
    {
        P->good = 0;
        return -1;
    }

    /* compile [start, end) */
    i = P->length;
    fmpz_crt_fit_length(P, i + 1);
    fmpz_init(P->prog[i].modulus);
    fmpz_init(P->prog[i].idem);
    P->good = P->good && fmpz_invmod(P->prog[i].modulus, leftmodulus, rightmodulus);
    fmpz_mul(P->prog[i].idem, leftmodulus, P->prog[i].modulus);
    fmpz_mul(P->prog[i].modulus, leftmodulus, rightmodulus);
    P->prog[i].a_idx = ret_idx;
    P->prog[i].b_idx = b_idx;
    P->prog[i].c_idx = c_idx;
    P->length = i + 1;

    return i;
}

/*
    Return 1 if moduli can be CRT'ed, 0 otherwise.
    A return of 0 means that future calls to run will leave output undefined.
*/
int fmpz_crt_compile(fmpz_crt_t P, fmpz ** moduli, slong len)
{
    slong i, j;
    slong * perm;
    TMP_INIT;

    FLINT_ASSERT(len > 0);

    TMP_START;
    perm = (slong *) TMP_ALLOC(len * sizeof(slong));

    for (i = 0; i < len; i++)
    {
        perm[i] = i;
    }

    /* make perm sort the degs so that degs[perm[j-1]] <= degs[perm[j-0]] */
    for (i = 1; i < len; i++)
    {
        for (j = i; j > 0 && fmpz_bits(moduli[perm[j-1]])
                           > fmpz_bits(moduli[perm[j-0]]); j--)
        {
            slong temp = perm[j-1];
            perm[j-1] = perm[j-0];
            perm[j-0] = temp;
        }
    }

    fmpz_crt_fit_length(P, FLINT_MAX(WORD(1), len - 1));
    fmpz_crt_set_length(P, 0);
    P->localsize = 1;
    P->good = 1;

    if (1 < len)
    {
        _push_prog(P, moduli, perm, 0, 0, len);
    }
    else
    {
        /*
            There is only one modulus. Lets compute as
                output[0] = input[0] + 0*(input[0] - input[0]) mod moduli[0]
        */
        i = 0;
        fmpz_init(P->prog[i].modulus);
        fmpz_init(P->prog[i].idem);
        fmpz_set(P->prog[i].modulus, moduli[0]);
        P->prog[i].a_idx = 0;
        P->prog[i].b_idx = -WORD(1);
        P->prog[i].c_idx = -WORD(1);
        P->length = i + 1;

        P->good = !fmpz_is_zero(moduli[0]);
    }

    if (!P->good)
    {
        fmpz_crt_set_length(P, 0);
    }

    /* two more spots for temporaries */
    P->temp1loc = P->localsize++;
    P->temp2loc = P->localsize++;

    TMP_END;

    return P->good;
}

/*
    If P was set with a call to nmod_poly_crt_compile(P, m, len), return
    in outputs[0] polynomial r of smallest degree such that
        r = inputs[0] mod m[0]
        r = inputs[1] mod m[1]
            ...
        r = inputs[len-1] mod m[len-1]
    For thread safety "outputs" is expected to have enough space for all
    temporaries, thus should be at least as long as P->localsize.
*/

void fmpz_crt_run(const fmpz_crt_t P, fmpz_t output, fmpz ** inputs)
{
    slong i;
    slong a, b, c;
    fmpz * A, * B, * C, * t1, * t2, * outputs;
    TMP_INIT;

    TMP_START;
    outputs = (fmpz *) TMP_ALLOC(P->localsize*sizeof(fmpz));
    for (i = 0; i < P->localsize; i++)
    {
        fmpz_init(outputs + i);
    }

    fmpz_swap(output, outputs + 0);

    t1 = outputs + P->temp1loc;
    t2 = outputs + P->temp2loc;
    for (i = 0; i < P->length; i++)
    {
        a = P->prog[i].a_idx;
        b = P->prog[i].b_idx;
        c = P->prog[i].c_idx;
        A = outputs + a;
        B = b < 0 ? inputs[-b-1] : outputs + b;
        C = c < 0 ? inputs[-c-1] : outputs + c;

        /* A = B + I*(C - B) mod M */
        fmpz_sub(t1, B, C);
        fmpz_mul(t2, P->prog[i].idem, t1);
        fmpz_sub(t1, B, t2);
        fmpz_mods(A, t1, P->prog[i].modulus);
    }

    fmpz_swap(output, outputs + 0);

    for (i = 0; i < P->localsize; i++)
    {
        fmpz_clear(outputs + i);
    }

    TMP_END;
}















void fmpz_mpolyd_swap(fmpz_mpolyd_t A, fmpz_mpolyd_t B)
{
    fmpz_mpolyd_struct t = *A;
    *A = *B;
    *B = t;
}

void fmpz_mpolyd_ctx_init(fmpz_mpolyd_ctx_t dctx, slong nvars)
{
    slong i;

    dctx->nvars = nvars;
    dctx->perm = (slong *) flint_malloc(nvars*sizeof(slong));
    for (i = 0; i < nvars; i++)
    {
        dctx->perm[i] = i;
    }
}


int fmpz_mpolyd_ctx_init_version1(fmpz_mpolyd_ctx_t dctx,
                            const fmpz_mpoly_t A, const fmpz_mpoly_t B,
                                                    const fmpz_mpoly_ctx_t ctx)
{
    int success = 0;
    slong i, j, degb_prod;
    slong * Aexps, * Bexps, * deg_bounds;
    slong nvars = ctx->minfo->nvars;
    slong * perm = dctx->perm;
    TMP_INIT;

    TMP_START;
    Aexps = (slong *) TMP_ALLOC(nvars*sizeof(slong));
    Bexps = (slong *) TMP_ALLOC(nvars*sizeof(slong));
    if (A->bits > FLINT_BITS || B->bits > FLINT_BITS)
        goto cleanup;
    fmpz_mpoly_degrees_si(Aexps, A, ctx);
    fmpz_mpoly_degrees_si(Bexps, B, ctx);

    deg_bounds = (slong *) TMP_ALLOC(nvars*sizeof(slong));
    for (i = 0; i < nvars; i++)
    {
        dctx->perm[i] = i;
    }

    degb_prod = 1;
    for (i = 0; i < nvars; i++)
    {
        ulong hi;
        deg_bounds[i] = FLINT_MAX(Aexps[i] + 1, Bexps[i] + 1);
        umul_ppmm(hi, degb_prod, degb_prod, deg_bounds[i]);
        if (hi != WORD(0) || degb_prod < 0)
            goto cleanup;
    }

    success = 1;
    for (i = 1; i < nvars; i++)
    {
        for (j = i; (j > 0) && (deg_bounds[j-1] < deg_bounds[j-0]); j--)
        {
            slong t1, t2;
            t1 = deg_bounds[j-1];
            t2 = deg_bounds[j-0];
            deg_bounds[j-0] = t1;
            deg_bounds[j-1] = t2;
            t1 = perm[j-1];
            t2 = perm[j-0];
            perm[j-0] = t1;
            perm[j-1] = t2;
        }
    }

cleanup:
    TMP_END;
    return success;
}


void fmpz_mpolyd_ctx_clear(fmpz_mpolyd_ctx_t dctx)
{
    flint_free(dctx->perm);
}


void fmpz_mpolyd_init(fmpz_mpolyd_t poly, slong nvars)
{
    slong i;

    poly->nvars = nvars;
    poly->degb_alloc = nvars;
    poly->deg_bounds = (slong *) flint_malloc(poly->degb_alloc*sizeof(slong));
    for (i = 0; i < nvars; i++)
    {
        poly->deg_bounds[i] = WORD(1);
    }
    poly->coeff_alloc = WORD(16);
    poly->coeffs = (fmpz *) flint_malloc(poly->coeff_alloc*sizeof(fmpz));
    for (i = 0; i < poly->coeff_alloc; i++)
    {
        fmpz_init(poly->coeffs + i);
    }
}


void fmpz_mpolyd_fit_length(fmpz_mpolyd_t poly, slong len)
{
    if (poly->coeff_alloc < len) {
        slong i;
        poly->coeffs = (fmpz *) flint_realloc(poly->coeffs, len*sizeof(fmpz));
        for (i = poly->coeff_alloc; i < len; i++)
        {
            fmpz_init(poly->coeffs + i);
        }
        poly->coeff_alloc = len;
    }
}


void fmpz_mpolyd_set_nvars(fmpz_mpolyd_t poly, slong nvars)
{
    poly->nvars = nvars;
    if (poly->degb_alloc < nvars) {
        poly->deg_bounds = (slong *) flint_realloc(poly->deg_bounds, nvars*sizeof(slong));
        poly->degb_alloc = nvars;
    }
}


void fmpz_mpolyd_zero(fmpz_mpolyd_t poly)
{
    slong i;

    for (i = 0; i < poly->nvars; i++)
    {
        poly->deg_bounds[i] = WORD(1);
    }
    poly->coeffs[0] = UWORD(0);
}

void fmpz_mpolyd_set_fmpz(fmpz_mpolyd_t poly, fmpz_t num)
{
    slong i;

    for (i = 0; i < poly->nvars; i++)
    {
        poly->deg_bounds[i] = WORD(1);
    }
    fmpz_set(poly->coeffs + 0, num);
}


void fmpz_mpolyd_clear(fmpz_mpolyd_t poly)
{
    slong i;

    for (i = 0; i < poly->coeff_alloc; i++)
    {
        fmpz_clear(poly->coeffs + i);
    }

    flint_free(poly->deg_bounds);
    flint_free(poly->coeffs);
    poly->deg_bounds = NULL;
    poly->coeffs = NULL;
}


void fmpz_mpoly_convert_to_fmpz_mpolyd(
                            fmpz_mpolyd_t poly1, const fmpz_mpolyd_ctx_t dctx,
                          const fmpz_mpoly_t poly2, const fmpz_mpoly_ctx_t ctx)
{
    slong degb_prod;
    slong i, j, N;
    slong * exps;
    const slong * perm = dctx->perm;
    slong nvars = ctx->minfo->nvars;
    TMP_INIT;

    fmpz_mpolyd_set_nvars(poly1, ctx->minfo->nvars);

    FLINT_ASSERT(poly2->bits <= FLINT_BITS);

    if (poly2->length == 0)
    {
        fmpz_mpolyd_zero(poly1);
        return;
    }

    TMP_START;
    exps = (slong *) TMP_ALLOC(ctx->minfo->nvars*sizeof(slong));

    fmpz_mpoly_degrees_si(exps, poly2, ctx);
    degb_prod = WORD(1);
    for (i = 0; i < nvars; i++)
    {
        poly1->deg_bounds[i] = exps[perm[i]] + 1;
        degb_prod *= poly1->deg_bounds[i];
    }

    fmpz_mpolyd_fit_length(poly1, degb_prod);
    for (i = 0; i < degb_prod; i++)
    {
        fmpz_zero(poly1->coeffs + i);
    }

    N = mpoly_words_per_exp(poly2->bits, ctx->minfo);
    for (i = 0; i < poly2->length; i++)
    {
        slong off = 0;

        mpoly_get_monomial_ui((ulong *)exps, poly2->exps + N*i, poly2->bits, ctx->minfo);
        for (j = 0; j < nvars; j++)
        {
            off = exps[perm[j]] + poly1->deg_bounds[j]*off;
        }

        fmpz_set(poly1->coeffs + off, poly2->coeffs + i);
    }

    TMP_END;
}


/*
    m is the number of variables in A
*/
void fmpz_mpoly_to_fmpz_mpolyd_perm_deflate(fmpz_mpolyd_t A, slong m,
              const fmpz_mpoly_t B, const slong * perm, const ulong * shift,
        const ulong * stride, const ulong * degree, const fmpz_mpoly_ctx_t ctx)
{
    slong n = ctx->minfo->nvars;
    slong degb_prod;
    slong i, k, l, N;
    ulong * Bexp;
    TMP_INIT;

    FLINT_ASSERT(m <= n);
    FLINT_ASSERT(B->bits <= FLINT_BITS);
    FLINT_ASSERT(B->length > 0);

    fmpz_mpolyd_set_nvars(A, m);

    TMP_START;
    Bexp = (ulong *) TMP_ALLOC(n*sizeof(slong));

    degb_prod = WORD(1);
    for (k = 0; k < m; k++)
    {
        l = perm[k];
        FLINT_ASSERT(stride[l] != UWORD(0));
        FLINT_ASSERT((degree[l] - shift[l]) % stride[l] == UWORD(0));
        A->deg_bounds[k] = (degree[l] - shift[l])/stride[l] + 1;
        degb_prod *= A->deg_bounds[k];
        /* we should not be converting something whose dense size overflows */
        FLINT_ASSERT(degb_prod > 0);
        FLINT_ASSERT(degb_prod >= A->deg_bounds[k]);
    }

    fmpz_mpolyd_fit_length(A, degb_prod);
    for (i = 0; i < degb_prod; i++)
    {
        fmpz_zero(A->coeffs + i);
    }

    N = mpoly_words_per_exp(B->bits, ctx->minfo);
    for (i = 0; i < B->length; i++)
    {
        slong off = 0;
        mpoly_get_monomial_ui(Bexp, B->exps + N*i, B->bits, ctx->minfo);
        for (k = 0; k < m; k++)
        {
            l = perm[k];
            FLINT_ASSERT(stride[l] != UWORD(0));
            FLINT_ASSERT(((Bexp[l] - shift[l]) % stride[l]) == UWORD(0));
            FLINT_ASSERT((Bexp[l] - shift[l])/stride[l] < A->deg_bounds[k]);
            off = (Bexp[l] - shift[l])/stride[l] + A->deg_bounds[k]*off;
        }
        fmpz_set(A->coeffs + off, B->coeffs + i);
    }

    TMP_END;
}


void fmpz_mpoly_from_fmpz_mpolyd_perm_inflate(fmpz_mpoly_t A,
         mp_bitcnt_t Abits, const fmpz_mpoly_ctx_t ctx, const fmpz_mpolyd_t B,
                 const slong * perm, const ulong * shift, const ulong * stride)
{
    slong off;
    slong n = ctx->minfo->nvars;
    slong m = B->nvars;
    slong Alen;
    slong i, j, l, k, N;
    slong perm_nontrivial;
    ulong topmask;
    ulong * exps, * pcurexp, * pexps;
    TMP_INIT;

    FLINT_ASSERT(m <= n);
    FLINT_ASSERT(Abits <= FLINT_BITS);

    perm_nontrivial = n - m;

    /* we are going to push back terms manually */
    Alen = 0;
    fmpz_mpoly_zero(A, ctx);
    fmpz_mpoly_fit_bits(A, Abits, ctx);
    A->bits = Abits;

    N = mpoly_words_per_exp(Abits, ctx->minfo);

    TMP_START;

    /* find exponent vector for all variables in B */
    pexps = (ulong *) TMP_ALLOC(N*m*sizeof(ulong));
    for (k = 0; k < m; k++)
    {
        l = perm[k];
        perm_nontrivial |= l - k;
        mpoly_gen_monomial_sp(pexps + k*N, l, Abits, ctx->minfo);
        mpoly_monomial_mul_ui(pexps + k*N, pexps + k*N, N, stride[l]);
    }

    /* get most significant exponent in pcurexp and its vector in exps */
    pcurexp = (ulong *) TMP_ALLOC(N*sizeof(ulong));
    exps = (ulong *) TMP_ALLOC(m*sizeof(ulong));
    off = WORD(1);
    for (j = 0; j < m; j++)
    {
        off *= B->deg_bounds[j];
    }
    FLINT_ASSERT(off <= B->coeff_alloc);
    off--;
    mpoly_set_monomial_ui(pcurexp, shift, Abits, ctx->minfo);
    i = off;
    for (k = m - 1; k >= 0; k--) 
    {
        exps[k] = i % B->deg_bounds[k];
        i = i / B->deg_bounds[k];
        mpoly_monomial_madd(pcurexp, pcurexp, exps[k], pexps + N*k, N);
    }

    /* scan down through the exponents */
    topmask = 0;

    for (; off >= 0; off--)
    {
        if (!fmpz_is_zero(B->coeffs + off))
        {
            _fmpz_mpoly_fit_length(&A->coeffs, &A->exps, &A->alloc, Alen + 1, N);
            fmpz_set(A->coeffs + Alen, B->coeffs + off);
            mpoly_monomial_set(A->exps + N*Alen, pcurexp, N);
            topmask |= pcurexp[N - 1];
            Alen++;
        }

        k = m - 1;
        do {
            --exps[k];
            if ((slong)(exps[k]) < WORD(0))
            {
                FLINT_ASSERT(off == 0 || k > 0);
                FLINT_ASSERT(exps[k] == -UWORD(1));
                exps[k] = B->deg_bounds[k] - 1;
                mpoly_monomial_madd(pcurexp, pcurexp, exps[k], pexps + N*k, N);
            }
            else
            {
                mpoly_monomial_sub(pcurexp, pcurexp, pexps + N*k, N);
                break;
            }
        } while (--k >= 0);
    }
    _fmpz_mpoly_set_length(A, Alen, ctx);


    /* sort the terms if needed */
    if (ctx->minfo->ord != ORD_LEX || perm_nontrivial != WORD(0))
    {
        slong msb;
        mpoly_get_cmpmask(pcurexp, N, Abits, ctx->minfo);
        if (topmask != WORD(0))
        {
            count_leading_zeros(msb, topmask);
            msb = (FLINT_BITS - 1)^msb;
        }
        else
        {
            msb = -WORD(1);
        }
        if (N == 1)
        {
            if (msb >= WORD(0))
            {
                _fmpz_mpoly_radix_sort1(A, 0, A->length,
                                                   msb, pcurexp[0], topmask);
            }
        }
        else
        {
            _fmpz_mpoly_radix_sort(A, 0, A->length,
                                        (N - 1)*FLINT_BITS + msb, N, pcurexp);
        }
    }

    TMP_END;
}


void fmpz_mpolyd_print(fmpz_mpolyd_t poly, const char ** vars,
                                                  const fmpz_mpolyd_ctx_t dctx)
{
    int first = 0;
    slong i, j;
    slong degb_prod;

    degb_prod = WORD(1);
    for (j = 0; j < poly->nvars; j++) {
        degb_prod *= poly->deg_bounds[j];
    }

    first = 1;
    for (i = 0; i < degb_prod; i++) {
        ulong k = i;

        if (fmpz_is_zero(poly->coeffs + i))
            continue;

        if (!first)
            printf(" + ");

        fmpz_print(poly->coeffs + i);

        for (j = poly->nvars - 1; j >= 0; j--)
        {
            ulong m = poly->deg_bounds[j];
            ulong e = k % m;
            k = k / m;
            flint_printf("*%s^%wd", vars[dctx->perm[j]], e);
        }
        FLINT_ASSERT(k == 0);
        first = 0;
    }

    if (first)
        flint_printf("0");
}


void fmpz_mpolyd_print_simple(fmpz_mpolyd_t poly)
{
    int first = 0;
    slong i, j;
    slong degb_prod;

    degb_prod = WORD(1);
    for (j = 0; j < poly->nvars; j++) {
        degb_prod *= poly->deg_bounds[j];
    }

    first = 1;
    for (i = 0; i < degb_prod; i++) {
        ulong k = i;

        if (fmpz_is_zero(poly->coeffs + i))
            continue;

        if (!first)
            printf(" + ");

        fmpz_print(poly->coeffs + i);

        for (j = poly->nvars - 1; j >= 0; j--)
        {
            ulong m = poly->deg_bounds[j];
            ulong e = k % m;
            k = k / m;
            flint_printf("*x%d^%wd", j, e);
        }
        FLINT_ASSERT(k == 0);
        first = 0;
    }

    if (first)
        flint_printf("0");
}


void fmpz_mpoly_convert_from_fmpz_mpolyd(
                                  fmpz_mpoly_t A, const fmpz_mpoly_ctx_t ctx,
                           const fmpz_mpolyd_t B, const fmpz_mpolyd_ctx_t dctx)
{
    slong i, j;
    slong degb_prod;
    slong * perm = dctx->perm;
    ulong * exps;
    TMP_INIT;

    FLINT_ASSERT(ctx->minfo->nvars == B->nvars);

    degb_prod = WORD(1);
    for (j = 0; j < B->nvars; j++) {
        degb_prod *= B->deg_bounds[j];
    }

    TMP_START;
    exps = (ulong *) TMP_ALLOC(B->nvars*sizeof(ulong));

    fmpz_mpoly_zero(A, ctx);
    for (i = 0; i < degb_prod; i++) {
        ulong k = i;

        if (fmpz_is_zero(B->coeffs + i))
            continue;

        for (j = B->nvars - 1; j >= 0; j--) 
        {
            ulong m = B->deg_bounds[j];
            ulong e = k % m;
            k = k / m;
            exps[perm[j]] = e;
        }
        FLINT_ASSERT(k == 0);

        fmpz_mpoly_set_coeff_fmpz_ui(A, B->coeffs + i, exps, ctx);
    }

    TMP_END;
}

int fmpz_mpolyd_CRT_nmod(fmpz_mpolyd_t A,
                                 const fmpz_mpolyd_t B, const fmpz_t Bm,
                                 const nmod_mpolyd_t C, const nmodf_ctx_t fctx)
{
    int Bok, Cok;
    slong carry;
    slong Bind, Cind;
    slong i, j;
    slong * inds;
    slong nvars = B->nvars;
    slong degb_prod;
    ulong hi, c;
    slong diff;
    fmpz_t zero;
    fmpz_t Bmn;
    slong * temp_deg_bounds;
    TMP_INIT;

    FLINT_ASSERT(B->nvars == C->nvars);

    TMP_START;
    temp_deg_bounds = (slong *) TMP_ALLOC(nvars*sizeof(slong));

    degb_prod = 1;
    diff = 0;
    for (j = 0; j < nvars; j++)
    {
        diff |= B->deg_bounds[j] - C->deg_bounds[j];
        temp_deg_bounds[j] = FLINT_MAX(B->deg_bounds[j], C->deg_bounds[j]);
        umul_ppmm(hi, degb_prod, degb_prod, temp_deg_bounds[j]);
        if (hi != WORD(0) || degb_prod < 0)
            return 0;
    }

    fmpz_init_set_ui(zero, 0);
    fmpz_init(Bmn);

    fmpz_mul_ui(Bmn, Bm, fctx->mod.n);
    c = fmpz_fdiv_ui(Bm, fctx->mod.n);
    c = n_invmod(c, fctx->mod.n);

    if (diff == 0) {
        /* both polynomials are packed into the same bounds */

        fmpz_mpolyd_set_nvars(A, nvars);
        fmpz_mpolyd_fit_length(A, degb_prod);
        for (j = 0; j < nvars; j++)
            A->deg_bounds[j] = temp_deg_bounds[j];

        for (i = 0; i < degb_prod; i++)
        {
            _fmpz_CRT_ui_precomp(A->coeffs + i, B->coeffs + i, Bm, C->coeffs[i], fctx->mod.n, fctx->mod.ninv, Bmn, c, 1);
        }

    } else {
        /* different bounds for packing */

        fmpz_mpolyd_t temp;
        fmpz_mpolyd_struct * T;

        if (A == B)
        {
            T = temp;
            fmpz_mpolyd_init(T, nvars);
        } else {
            T = A;
        }

        fmpz_mpolyd_set_nvars(T, nvars);
        fmpz_mpolyd_fit_length(T, degb_prod);
        for (j = 0; j < nvars; j++)
            T->deg_bounds[j] = temp_deg_bounds[j];

        inds = (slong *) TMP_ALLOC(nvars*sizeof(slong));
        for (j = 0; j < nvars; j++)
            inds[j] = 0;
        Bok = 1;
        Cok = 1;
        Bind = 0;
        Cind = 0;
        for (i = 0; i < degb_prod; i++)
        {
                   if (Bok && Cok) {
                _fmpz_CRT_ui_precomp(T->coeffs + i, B->coeffs + Bind++, Bm, C->coeffs[Cind++], fctx->mod.n, fctx->mod.ninv, Bmn, c, 1);
            } else if (Bok && !Cok) {
                _fmpz_CRT_ui_precomp(T->coeffs + i, B->coeffs + Bind++, Bm, 0                , fctx->mod.n, fctx->mod.ninv, Bmn, c, 1);
            } else if (!Bok && Cok) {
                _fmpz_CRT_ui_precomp(T->coeffs + i, zero              , Bm, C->coeffs[Cind++], fctx->mod.n, fctx->mod.ninv, Bmn, c, 1);
            } else {
                fmpz_zero(T->coeffs + i);
            }

            Bok = 1;
            Cok = 1;
            carry = 1;
            for (j = nvars - 1; j >= 0; j--)
            {
                inds[j] += carry;
                if (inds[j] < T->deg_bounds[j])
                {
                    carry = 0;
                    Bok = Bok && (inds[j] < B->deg_bounds[j]);
                    Cok = Cok && (inds[j] < C->deg_bounds[j]);
                } else
                {
                    carry = 1;
                    inds[j] = 0;
                }
            }
        }

        if (A == B)
        {
            fmpz_mpolyd_swap(A, T);
            fmpz_mpolyd_clear(T);
        } else {

        }
    }

    fmpz_clear(zero);
    fmpz_clear(Bmn);

    TMP_END;
    return 1;
}


void fmpz_mpolyd_height(fmpz_t max, fmpz_mpolyd_t A)
{
    slong degb_prod, i, j;
    fmpz_t t;

    degb_prod = WORD(1);
    for (j = 0; j < A->nvars; j++)
    {
        degb_prod *= A->deg_bounds[j];
    }

    fmpz_init(t);
    fmpz_zero(max);
    for (i = 0; i < degb_prod; i++)
    {
        fmpz_abs(t, A->coeffs + i);
        if (fmpz_cmp(max, t) < 0)
            fmpz_set(max, t);
    }

    fmpz_clear(t);
}

void fmpz_mpolyd_heights(fmpz_t max, fmpz_t sum, fmpz_mpolyd_t A)
{
    slong degb_prod, i, j;
    fmpz_t t;

    degb_prod = WORD(1);
    for (j = 0; j < A->nvars; j++)
        degb_prod *= A->deg_bounds[j];

    fmpz_init(t);
    fmpz_zero(max);
    fmpz_zero(sum);
    for (i = 0; i < degb_prod; i++)
    {
        fmpz_abs(t, A->coeffs + i);
        fmpz_add(sum, sum, t);
        if (fmpz_cmp(max, t) < 0)
            fmpz_set(max, t);
    }

    fmpz_clear(t);
}

void fmpz_mpolyd_divexact_fmpz_inplace(fmpz_mpolyd_t A, fmpz_t c)
{
    slong degb_prod, i, j;

    degb_prod = WORD(1);
    for (j = 0; j < A->nvars; j++)
        degb_prod *= A->deg_bounds[j];

    for (i = 0; i < degb_prod; i++)
    {
        fmpz_divexact(A->coeffs + i, A->coeffs + i, c);
    }
}

void fmpz_mpolyd_mul_scalar_inplace(fmpz_mpolyd_t A, fmpz_t c)
{
    slong degb_prod, i, j;

    degb_prod = WORD(1);
    for (j = 0; j < A->nvars; j++)
        degb_prod *= A->deg_bounds[j];

    for (i = 0; i < degb_prod; i++)
    {
        fmpz_mul(A->coeffs + i, A->coeffs + i, c);
    }
}

void fmpz_mpolyd_content(fmpz_t c, const fmpz_mpolyd_t A)
{
    slong degb_prod, j;

    degb_prod = WORD(1);
    for (j = 0; j < A->nvars; j++)
        degb_prod *= A->deg_bounds[j];

    _fmpz_vec_content(c, A->coeffs, degb_prod);
}

void fmpz_mpolyd_to_nmod_mpolyd(nmod_mpolyd_t Ap, fmpz_mpolyd_t A, const nmodf_ctx_t fctx)
{
    slong j;
    slong degb_prod;

    nmod_mpolyd_set_nvars(Ap, A->nvars);

    degb_prod = WORD(1);
    for (j = 0; j < A->nvars; j++)
    {
        Ap->deg_bounds[j] = A->deg_bounds[j];
        degb_prod *= A->deg_bounds[j];
    }

    nmod_mpolyd_fit_length(Ap, degb_prod);
    _fmpz_vec_get_nmod_vec(Ap->coeffs, A->coeffs, degb_prod, fctx->mod);
}

void fmpz_mpolyd_set_nmod_mpolyd(fmpz_mpolyd_t A, nmod_mpolyd_t Ap, const nmodf_ctx_t fctx)
{
    slong j;
    slong degb_prod;

    fmpz_mpolyd_set_nvars(A, Ap->nvars);

    degb_prod = WORD(1);
    for (j = 0; j < Ap->nvars; j++)
    {
        A->deg_bounds[j] = Ap->deg_bounds[j];
        degb_prod *= Ap->deg_bounds[j];
    }

    fmpz_mpolyd_fit_length(A, degb_prod);
    _fmpz_vec_set_nmod_vec(A->coeffs, Ap->coeffs, degb_prod, fctx->mod);
}



void fmpz_mpolyd_set(fmpz_mpolyd_t A, const fmpz_mpolyd_t B)
{
    slong j;
    slong degb_prod;

    fmpz_mpolyd_set_nvars(A, B->nvars);

    degb_prod = WORD(1);
    for (j = 0; j < B->nvars; j++)
    {
        A->deg_bounds[j] = B->deg_bounds[j];
        degb_prod *= B->deg_bounds[j];
    }

    fmpz_mpolyd_fit_length(A, degb_prod);
    _fmpz_vec_set(A->coeffs, B->coeffs, degb_prod);
}


slong fmpz_mpolyd_leadmon(slong * exps, const fmpz_mpolyd_t A)
{

    slong i, j, k;
    slong degb_prod;

    degb_prod = WORD(1);
    for (j = 0; j < A->nvars; j++)
    {
        degb_prod *= A->deg_bounds[j];
    }

    for (i = degb_prod-1; i >= 0; i--)
    {
        if (!fmpz_is_zero(A->coeffs + i))
            break;
    }

    FLINT_ASSERT(i>=0);

    k = i;
    for (j = A->nvars - 1; j >= 0; j--) 
    {
        ulong m = A->deg_bounds[j];
        ulong e = k % m;
        k = k / m;
        exps[j] = e;
    }
    FLINT_ASSERT(k == 0);

    return i;
}

void fmpz_mpolyd_lc(fmpz_t a, const fmpz_mpolyd_t A)
{

    slong i, j;
    slong degb_prod;

    degb_prod = WORD(1);
    for (j = 0; j < A->nvars; j++)
    {
        degb_prod *= A->deg_bounds[j];
    }

    for (i = degb_prod-1; i >= 0; i--)
    {
        if (!fmpz_is_zero(A->coeffs + i))
        {
            fmpz_set(a, A->coeffs + i);
            return;
        }
    }
}



int fmpz_mpolyd_gcd_brown(fmpz_mpolyd_t G,
            fmpz_mpolyd_t Abar, fmpz_mpolyd_t Bbar,
                    fmpz_mpolyd_t A, fmpz_mpolyd_t B)
{
    int equal, success = 1;
    mp_limb_t p, old_p;
    slong j, nvars;
    slong lm_idx;
    slong * exp, * texp;
    fmpz_t gamma, m;
    fmpz_t gnm, gns, anm, ans, bnm, bns;
    fmpz_t lA, lB, cA, cB, cG, bound, temp, pp;
    nmod_mpolyd_t Gp, Apbar, Bpbar, Ap, Bp;
    nmodf_ctx_t fctx;
    TMP_INIT;

    TMP_START;

    nmodf_ctx_init(fctx, 2);
    nvars = A->nvars;

    nmod_mpolyd_init(Gp, nvars);
    nmod_mpolyd_init(Apbar, nvars);
    nmod_mpolyd_init(Bpbar, nvars);
    nmod_mpolyd_init(Ap, nvars);
    nmod_mpolyd_init(Bp, nvars);

    fmpz_init(cA);
    fmpz_init(cB);
    fmpz_init(cG);
    fmpz_init(lA);
    fmpz_init(lB);
    fmpz_init(gamma);
    fmpz_init(gnm);
    fmpz_init(gns);
    fmpz_init(anm);
    fmpz_init(ans);
    fmpz_init(bnm);
    fmpz_init(bns);
    fmpz_init(bound);
    fmpz_init(temp);
    fmpz_init_set_si(m, 1);
    fmpz_init(pp);

    fmpz_mpolyd_content(cA, A);
    fmpz_mpolyd_content(cB, B);
    fmpz_gcd(cG, cA, cB);
    fmpz_mpolyd_divexact_fmpz_inplace(A, cA);
    fmpz_mpolyd_divexact_fmpz_inplace(B, cB);
    fmpz_mpolyd_lc(lA, A);
    fmpz_mpolyd_lc(lB, B);
    fmpz_gcd(gamma, lA, lB);
    fmpz_mpolyd_height(bound, A);
    fmpz_mpolyd_height(temp, B);

    if (fmpz_cmp(bound, temp) < 0)
        fmpz_swap(bound, temp);
    fmpz_mul(bound, bound, gamma);
    fmpz_add(bound, bound, bound);

    exp = (slong *) TMP_ALLOC(nvars*sizeof(slong));
    texp = (slong *) TMP_ALLOC(nvars*sizeof(slong));

    fmpz_mpolyd_leadmon(exp, A);
    fmpz_mpolyd_leadmon(texp, B);
    for (j = 0; j < nvars; j++)
        exp[j] = FLINT_MIN(exp[j], texp[j]);

    p = UWORD(1) << (FLINT_BITS - 1);

choose_next_prime:

    old_p = p;
    p = n_nextprime(p, 1);
    if (p <= old_p) {
        /* ran out of primes */
        success = 0;
        goto done;
    }
    fmpz_set_ui(pp, p);
    if (fmpz_divisible(lA, pp) || fmpz_divisible(lB, pp))
        goto choose_next_prime;

    nmodf_ctx_reset(fctx, p);
    fmpz_mpolyd_to_nmod_mpolyd(Ap, A, fctx);
    fmpz_mpolyd_to_nmod_mpolyd(Bp, B, fctx);
    success = nmod_mpolyd_gcd_brown_smprime(Gp, Apbar, Bpbar, Ap, Bp, fctx);
    if (!success)
        goto choose_next_prime;

    lm_idx = nmod_mpolyd_leadmon(texp, Gp);
    if (lm_idx <= 0)
    {
        /* Gp is 1, which means A and B are r.p. */
        FLINT_ASSERT(lm_idx == 0);
        fmpz_mpolyd_set_fmpz(G, cG);
        fmpz_mpolyd_swap(Abar, A);
        fmpz_divexact(temp, cA, cG);
        fmpz_mpolyd_mul_scalar_inplace(Abar, temp);
        fmpz_mpolyd_swap(Bbar, B);
        fmpz_divexact(temp, cB, cG);
        fmpz_mpolyd_mul_scalar_inplace(Bbar, temp);
        goto done;
    }

    equal = 1;
    for (j = 0; j < nvars; j++)
    {
        if (texp[j] > exp[j])
        {
            goto choose_next_prime;
        } else if (texp[j] < exp[j])
        {
            equal = 0;
            break;
        }
    }

    nmod_mpolyd_mul_scalar(Gp, fmpz_fdiv_ui(gamma, p), fctx);

    if (fmpz_is_one(m) || !equal)
    {
        fmpz_mpolyd_set_nmod_mpolyd(G, Gp, fctx);
        fmpz_mpolyd_set_nmod_mpolyd(Abar, Apbar, fctx);
        fmpz_mpolyd_set_nmod_mpolyd(Bbar, Bpbar, fctx);
        fmpz_set_ui(m, p);
        for (j = 0; j < nvars; j++)
            exp[j] = texp[j];

        goto choose_next_prime;
    }

    success = 1;
    success = success && fmpz_mpolyd_CRT_nmod(G, G, m, Gp, fctx);
    success = success && fmpz_mpolyd_CRT_nmod(Abar, Abar, m, Apbar, fctx);
    success = success && fmpz_mpolyd_CRT_nmod(Bbar, Bbar, m, Bpbar, fctx);
    fmpz_mul(m, m, pp);
    if (!success)
        goto done;

    if (fmpz_cmp(m, bound) <= 0)
        goto choose_next_prime;

    fmpz_mpolyd_heights(gnm, gns, G);
    fmpz_mpolyd_heights(anm, ans, Abar);
    fmpz_mpolyd_heights(bnm, bns, Bbar);
    fmpz_mul(ans, ans, gnm);
    fmpz_mul(anm, anm, gns);
    fmpz_mul(bns, bns, gnm);
    fmpz_mul(bnm, bnm, gns);

    if (fmpz_cmp(ans, anm) > 0)
        fmpz_swap(ans, anm);
    if (fmpz_cmp(bns, bnm) > 0)
        fmpz_swap(bns, bnm);
    fmpz_add(ans, ans, ans);
    fmpz_add(bns, bns, bns);
    if (fmpz_cmp(ans, m) >= 0 || fmpz_cmp(bns, m) >= 0)
        goto choose_next_prime;

    fmpz_mpolyd_content(temp, G);
    fmpz_mpolyd_divexact_fmpz_inplace(G, temp);
    fmpz_mpolyd_lc(temp, G);
    fmpz_mpolyd_divexact_fmpz_inplace(Abar, temp);
    fmpz_mpolyd_divexact_fmpz_inplace(Bbar, temp);

    fmpz_mpolyd_mul_scalar_inplace(G, cG);
    fmpz_divexact(temp, cA, cG);
    fmpz_mpolyd_mul_scalar_inplace(Abar, temp);
    fmpz_divexact(temp, cB, cG);
    fmpz_mpolyd_mul_scalar_inplace(Bbar, temp);

done:

    fmpz_clear(cA);
    fmpz_clear(cB);
    fmpz_clear(cG);
    fmpz_clear(lA);
    fmpz_clear(lB);
    fmpz_clear(gamma);
    fmpz_clear(gnm);
    fmpz_clear(gns);
    fmpz_clear(anm);
    fmpz_clear(ans);
    fmpz_clear(bnm);
    fmpz_clear(bns);
    fmpz_clear(bound);
    fmpz_clear(temp);
    fmpz_clear(m);
    fmpz_clear(pp);

    nmod_mpolyd_clear(Gp);
    nmod_mpolyd_clear(Apbar);
    nmod_mpolyd_clear(Bpbar);
    nmod_mpolyd_clear(Ap);
    nmod_mpolyd_clear(Bp);

    nmodf_ctx_clear(fctx);

    TMP_END;
    return success;
}


int fmpz_mpoly_gcd_brown(fmpz_mpoly_t G,
                              const fmpz_mpoly_t A, const fmpz_mpoly_t B,
                                                    const fmpz_mpoly_ctx_t ctx)
{
    int success;
    fmpz_mpolyd_t Ad, Bd, Gd, Abar, Bbar;
    fmpz_mpolyd_ctx_t dctx;
    slong nvars = ctx->minfo->nvars;

    success = 1;

    if (fmpz_mpoly_is_zero(A, ctx)) {
        fmpz_mpoly_set(G, B, ctx);
        goto cleanup_stage0;
    }
    if (fmpz_mpoly_is_zero(B, ctx)) {
        fmpz_mpoly_set(G, A, ctx);
        goto cleanup_stage0;
    }

    fmpz_mpolyd_ctx_init(dctx, nvars);
    success = fmpz_mpolyd_ctx_init_version1(dctx, A, B, ctx);
    if (!success)
    {
        fmpz_mpoly_zero(G, ctx);
        goto cleanup_stage1;
    }

    fmpz_mpolyd_init(Ad, nvars);
    fmpz_mpolyd_init(Bd, nvars);
    fmpz_mpolyd_init(Gd, nvars);
    fmpz_mpolyd_init(Abar, nvars);
    fmpz_mpolyd_init(Bbar, nvars);

    fmpz_mpoly_convert_to_fmpz_mpolyd(Ad, dctx, A, ctx);
    fmpz_mpoly_convert_to_fmpz_mpolyd(Bd, dctx, B, ctx);

    success = fmpz_mpolyd_gcd_brown(Gd, Abar, Bbar, Ad, Bd);
    if (!success)
    {
        fmpz_mpoly_zero(G, ctx);
    } else
    {
        fmpz_mpoly_convert_from_fmpz_mpolyd(G, ctx, Gd, dctx);
    }

    fmpz_mpolyd_clear(Bbar);
    fmpz_mpolyd_clear(Abar);
    fmpz_mpolyd_clear(Gd);
    fmpz_mpolyd_clear(Bd);
    fmpz_mpolyd_clear(Ad);

cleanup_stage1:

    fmpz_mpolyd_ctx_clear(dctx);

cleanup_stage0:

    if (success && (G->length > 0) && (fmpz_sgn(G->coeffs + 0) < 0))
        fmpz_mpoly_neg(G, G, ctx);

    return success;
}





int fmpz_mpolyu_is_canonical(const fmpz_mpolyu_t A, const fmpz_mpoly_ctx_t ctx)
{
    slong i;

    if (A->length > A->alloc)
    {
        return 0;
    }

    for (i = 0; i < A->length; i++)
    {
        if (!fmpz_mpoly_is_canonical(A->coeffs + i, ctx))
        {
            return 0;
        }

        if (i > 0 && A->exps[i - 1] <= A->exps[i])
        {
            return 0;
        }
    }

    return 1;
}

void fmpz_mpolyu_content(fmpz_t g, const fmpz_mpolyu_t A, const fmpz_mpoly_ctx_t ctx)
{
    slong i, j;

    fmpz_zero(g);
    for (i = 0; i < A->length; i++)
    {
        fmpz_mpoly_struct * Ac = A->coeffs + i;
        for (j = 0; j < Ac->length; j++)
        {
            fmpz_gcd(g, g, Ac->coeffs + j);
            if (fmpz_is_one(g))
                return;
        }
    }
}

void fmpz_mpolyu_scalar_divexact(fmpz_mpolyu_t A, const fmpz_t c, const fmpz_mpoly_ctx_t ctx)
{
    slong i, j;

    if (fmpz_is_one(c))
        return;

    for (i = 0; i < A->length; i++)
    {
        fmpz_mpoly_struct * Ac = A->coeffs + i;
        for (j = 0; j < Ac->length; j++)
        {
            fmpz_divexact(Ac->coeffs + j, Ac->coeffs + j, c);
        }
    }
}


void fmpz_mpolyu_scalar_mul_fmpz(fmpz_mpolyu_t A, const fmpz_t c, const fmpz_mpoly_ctx_t ctx)
{
    slong i, j;

    if (fmpz_is_one(c))
        return;

    for (i = 0; i < A->length; i++)
    {
        fmpz_mpoly_struct * Ac = A->coeffs + i;
        for (j = 0; j < Ac->length; j++)
        {
            fmpz_mul(Ac->coeffs + j, Ac->coeffs + j, c);
        }
    }
}


void fmpz_mpolyu_height(fmpz_t max, const fmpz_mpolyu_t A, const fmpz_mpoly_ctx_t ctx)
{
    slong i, j;
    fmpz_t t;

    fmpz_init(t);
    fmpz_zero(max);

    for (i = 0; i < A->length; i++)
    {
        fmpz_mpoly_struct * Ac = A->coeffs + i;
        for (j = 0; j < Ac->length; j++)
        {
            fmpz_abs(t, Ac->coeffs + j);
            if (fmpz_cmp(max, t) < 0)
                fmpz_set(max, t);
        }
    }

    fmpz_clear(t);
}


void fmpz_mpolyu_heights(fmpz_t max, fmpz_t sum, const fmpz_mpolyu_t A, const fmpz_mpoly_ctx_t ctx)
{
    slong i, j;
    fmpz_t t;

    fmpz_init(t);
    fmpz_zero(max);
    fmpz_zero(sum);

    for (i = 0; i < A->length; i++)
    {
        fmpz_mpoly_struct * Ac = A->coeffs + i;
        for (j = 0; j < Ac->length; j++)
        {
            fmpz_abs(t, Ac->coeffs + j);
            fmpz_add(sum, sum, t);
            if (fmpz_cmp(max, t) < 0)
                fmpz_set(max, t);
        }
    }

    fmpz_clear(t);
}


/*
    E = A mod p
    E is in Fp[x_0,...,x_(m-2)][x_(m-1)]
    A is in ZZ[x_0,...,x_(m-2), x_(m-1)]
*/
void fmpz_mpoly_redto_nmod_mpolyn(nmod_mpolyn_t E, const nmod_mpoly_ctx_t pctx,
                              const fmpz_mpoly_t A, const fmpz_mpoly_ctx_t ctx)
{
    slong N = mpoly_words_per_exp_sp(A->bits, ctx->minfo);
    slong offset, shift, k;
    ulong mask;
    mp_limb_t v;
    fmpz * Acoeff = A->coeffs;
    ulong * Aexp = A->exps;
    slong Alen = A->length;
    slong Ai;
    nmod_poly_struct * Ecoeff;
    ulong * Eexp;
    slong Ei;
    slong m = ctx->minfo->nvars;

    mpoly_gen_offset_shift_sp(&offset, &shift, m - 1, A->bits, ctx->minfo);
    mask = (-UWORD(1)) >> (FLINT_BITS - A->bits);

    Ecoeff = E->coeffs;
    Eexp = E->exps;
    Ei = 0;
    for (Ai = 0; Ai < Alen; Ai++)
    {
        v = fmpz_fdiv_ui(Acoeff + Ai, pctx->ffinfo->mod.n);
        k = ((Aexp + N*Ai)[offset] >> shift) & mask;
        if (v == 0)
        {
            continue;
        }

        if (Ei > 0 && mpoly_monomial_equal_extra(Eexp + N*(Ei - 1), Aexp + N*Ai, N, offset, -(k << shift)))
        {
            /* append to previous */
FLINT_ASSERT((Ecoeff + Ei - 1)->mod.n == pctx->ffinfo->mod.n);
            nmod_poly_set_coeff_ui(Ecoeff + Ei - 1, k, v);
        }
        else
        {
            FLINT_ASSERT(Ei == 0 || mpoly_monomial_gt_nomask_extra(Eexp + N*(Ei - 1), Aexp + N*Ai, N, offset, -(k << shift)));

            /* create new */
            if (Ei >= E->alloc)
            {
                nmod_mpolyn_fit_length(E, Ei + 1, pctx);
                Ecoeff = E->coeffs;
                Eexp = E->exps;
            }

FLINT_ASSERT((Ecoeff + Ei)->mod.n == pctx->ffinfo->mod.n);

            mpoly_monomial_set_extra(Eexp + N*Ei, Aexp + N*Ai, N, offset, -(k << shift));
            nmod_poly_zero(Ecoeff + Ei);
            nmod_poly_set_coeff_ui(Ecoeff + Ei, k, v);
            Ei++;
        }
    }
    E->length = Ei;
}



/*
    E = A mod p
    E is in Fp[X][x_0,...,x_(m-2)][x_(m-1)]
    A is in ZZ[X][x_0,...,x_(m-2), x_(m-1)]
*/
void fmpz_mpolyu_redto_nmod_mpolyun(nmod_mpolyun_t E, const nmod_mpoly_ctx_t pctx,
                             const fmpz_mpolyu_t A, const fmpz_mpoly_ctx_t ctx)
{
    fmpz_mpoly_struct * Acoeff = A->coeffs;
    ulong * Aexp = A->exps;
    slong Alen = A->length;
    slong Ai;
    nmod_mpolyn_struct * Ecoeff;
    ulong * Eexp;
    slong Ei;
/*
printf("reduce starting\n");
printf("A: "); fmpz_mpolyu_print_pretty(A, NULL, ctx); printf("\n");
*/

    nmod_mpolyun_fit_length(E, Alen, pctx);
    Ecoeff = E->coeffs;
    Eexp = E->exps;

    Ei = 0;
    for (Ai = 0; Ai < Alen; Ai++)
    {
        fmpz_mpoly_redto_nmod_mpolyn(Ecoeff + Ei, pctx, Acoeff + Ai, ctx);
        Eexp[Ei] = Aexp[Ai];
        Ei += ((Ecoeff + Ei)->length != 0);
    }
    E->length = Ei;
/*
printf("reduce returnig\n");
printf("A: "); fmpz_mpolyu_print_pretty(A, NULL, ctx); printf("\n");
printf("E: "); nmod_mpolyun_print_pretty(E, NULL, pctx); printf("\n");
*/
    FLINT_ASSERT(nmod_mpolyun_is_canonical(E, pctx));
}






/*
    A = B
    A is in Fp[x_0, ..., x_(var-1), x_(var-1)][x_var]
    B is in Fp[x_0, ..., x_(var-2)][x_(var-1)]
*/
void fmpz_mpoly_startinterp_n(fmpz_mpoly_t A, const fmpz_mpoly_ctx_t ctx,
                            const nmod_mpolyn_t B, const nmod_mpoly_ctx_t pctx)
{
    slong N = mpoly_words_per_exp(B->bits, ctx->minfo);
    slong offset, shift;
    slong vi;

    nmod_poly_struct * Bcoeff = B->coeffs;
    ulong * Bexp = B->exps;
    slong Blen = B->length;
    slong Bi;

    fmpz * Acoeff;
    ulong * Aexp;
    slong Ai;
    slong var = ctx->minfo->nvars;

    FLINT_ASSERT(var = pctx->minfo->nvars);

    fmpz_mpoly_fit_length(A, Blen, ctx);
    Acoeff = A->coeffs;
    Aexp = A->exps;

    mpoly_gen_offset_shift_sp(&offset, &shift, var - 1, A->bits, ctx->minfo);

    Ai = 0;
    for (Bi = 0; Bi < Blen; Bi++)
    {
        if (Ai + (Bcoeff + Bi)->length >= A->alloc)
        {
            fmpz_mpoly_fit_length(A, Ai + (Bcoeff + Bi)->length, ctx);
            Acoeff = A->coeffs;
            Aexp = A->exps;
        }
        for (vi = (Bcoeff + Bi)->length - 1; vi >= 0; vi--)
        {
            if ((Bcoeff + Bi)->coeffs[vi] != 0)
            {
                mpoly_monomial_set_extra(Aexp + N*Ai, Bexp + N*Bi, N, offset, vi << shift);
                fmpz_set_ui_smod(Acoeff + Ai, (Bcoeff + Bi)->coeffs[vi], pctx->ffinfo->mod.n);
                Ai++;
            }
        }
    }
    A->length = Ai;
}


/*
    A = B using symmetric range
    A is in ZZ[X][x_0, ..., x_(var-2), x_(m-1)]
    B is in Fp[X][x_0, ..., x_(var-2)][x_(m-1)]
*/
void fmpz_mpolyu_startinterp_un(fmpz_mpolyu_t A, const fmpz_mpoly_ctx_t ctx,
                           const nmod_mpolyun_t B, const nmod_mpoly_ctx_t pctx)
{
    slong i;
    nmod_mpolyn_struct * Bcoeff = B->coeffs;
    ulong * Bexp = B->exps;
    slong Blen = B->length;

    fmpz_mpoly_struct * Acoeff;
    ulong * Aexp;

    fmpz_mpolyu_fit_length(A, Blen, ctx);
    Acoeff = A->coeffs;
    Aexp = A->exps;

    for (i = 0; i < Blen; i++)
    {
        Aexp[i] = Bexp[i];
        fmpz_mpoly_startinterp_n(Acoeff + i, ctx, Bcoeff + i, pctx);
    }
    A->length = Blen;

    FLINT_ASSERT(fmpz_mpolyu_is_canonical(A, ctx));
}



void nmod_mpolyun_change_nmod_poly_mod(nmod_mpolyun_t A, const nmod_mpoly_ctx_t ctx)
{
    slong i, j;

    for (i = 0; i < A->alloc; i++)
    {
        nmod_mpolyn_struct * Ac = A->coeffs + i;
        for (j = 0; j < Ac->alloc; j++)
        {
            (Ac->coeffs + j)->mod = ctx->ffinfo->mod;


            FLINT_ASSERT((Ac->coeffs + j)->mod.n == ctx->ffinfo->mod.n);

        }
    }
}




/*
    T = F + modulus*((A - F(mod p))/(modulus (mod p)))
    no assumptions about matching monomials
    F is in Fp[x_0, ..., x_(var-1), x_(var-1)][x_var]
    A is in Fp[x_0, ..., x_(var-2)][x_(var-1)]
*/
int fmpz_mpoly_addinterp_n(
            fmpz_mpoly_t T, const fmpz_mpoly_t F, const fmpz_mpoly_ctx_t ctx,
            fmpz_t modulus, const nmod_mpolyn_t A, const nmod_mpoly_ctx_t pctx)
{
    int changed = 0;
    slong N = mpoly_words_per_exp_sp(T->bits, ctx->minfo);
    slong offset, shift;
    slong vi;

    fmpz * Tcoeff;
    ulong * Texp;
    slong Ti;

    nmod_poly_struct * Acoeff = A->coeffs;
    slong Alen = A->length;
    ulong * Aexp = A->exps;
    slong Ai;

    fmpz * Fcoeff = F->coeffs;
    slong Flen = F->length;
    ulong * Fexp = F->exps;
    slong Fi;
    fmpz_t zero;
    slong var = ctx->minfo->nvars;

    FLINT_ASSERT(var == pctx->minfo->nvars);

    fmpz_init(zero);

    FLINT_ASSERT(T->bits == A->bits);
    FLINT_ASSERT(F->bits == A->bits);

    mpoly_gen_offset_shift_sp(&offset, &shift, var - 1, A->bits, ctx->minfo);

    Flen = F->length;

    fmpz_mpoly_fit_length(T, FLINT_MAX(Flen, Alen), ctx);
    Tcoeff = T->coeffs;
    Texp = T->exps;
    Ti = 0;

    Fi = Ai = vi = 0;
    if (Ai < Alen)
    {
        vi = nmod_poly_degree(A->coeffs + Ai);
    }
    while (Fi < Flen || Ai < Alen)
    {
        if (Ti >= T->alloc)
        {
            fmpz_mpoly_fit_length(T, Ti + FLINT_MAX(Flen - Fi, Alen - Ai), ctx);
            Tcoeff = T->coeffs;
            Texp = T->exps;
        }

        if (Fi < Flen && Ai < Alen && mpoly_monomial_equal_extra(Fexp + N*Fi, Aexp + N*Ai, N, offset, vi << shift))
        {
            /* F term ok, A term ok */
            fmpz_CRT_ui(Tcoeff + Ti, Fcoeff + Fi, modulus, (Acoeff + Ai)->coeffs[vi], pctx->ffinfo->mod.n, 1);
            changed |= !fmpz_equal(Tcoeff + Ti, Fcoeff + Fi);
            mpoly_monomial_set(Texp + N*Ti, Fexp + N*Fi, N);

            Fi++;
            do {
                vi--;
            } while (vi >= 0 && (Acoeff + Ai)->coeffs[vi] == 0);
            if (vi < 0)
            {
                Ai++;
                if (Ai < Alen)
                {
                    vi = nmod_poly_degree(A->coeffs + Ai);
                }
            }
        }
        else if (Fi < Flen && (Ai >= Alen || mpoly_monomial_gt_nomask_extra(Fexp + N*Fi, Aexp + N*Ai, N, offset, vi << shift)))
        {
            /* F term ok, A term missing */
            fmpz_CRT_ui(Tcoeff + Ti, Fcoeff + Fi, modulus, 0, pctx->ffinfo->mod.n, 1);
            changed |= !fmpz_equal(Tcoeff + Ti, Fcoeff + Fi);

            mpoly_monomial_set(Texp + N*Ti, Fexp + N*Fi, N);

            Fi++;
        }
        else
        {
            FLINT_ASSERT(Ai < Alen && (Fi >= Flen || mpoly_monomial_lt_nomask_extra(Fexp + N*Fi, Aexp + N*Ai, N, offset, vi << shift)));

            /* F term missing, A term ok */
            fmpz_CRT_ui(Tcoeff + Ti, zero, modulus, (Acoeff + Ai)->coeffs[vi], pctx->ffinfo->mod.n, 1);            
            FLINT_ASSERT(!fmpz_is_zero(Tcoeff + Ti));
            changed = 1;
            mpoly_monomial_set_extra(Texp + N*Ti, Aexp + N*Ai, N, offset, vi << shift);

            do {
                vi--;
            } while (vi >= 0 && (Acoeff + Ai)->coeffs[vi] == 0);
            if (vi < 0)
            {
                Ai++;
                if (Ai < Alen)
                {
                    vi = nmod_poly_degree(A->coeffs + Ai);
                }
            }
        }

        FLINT_ASSERT(!fmpz_is_zero(Tcoeff + Ti));
        Ti++;
    }
    T->length = Ti;

    fmpz_clear(zero);
    return changed;
}

/*
    F = F + modulus*((A - F(mod p))/(modulus (mod p)))
    no assumptions about matching monomials
    F is in ZZ[X][x_0, ..., x_(m-1), x_(m-1)]
    A is in Fp[X][x_0, ..., x_(m-2)][x_(m-1)]
*/
int fmpz_mpolyu_addinterp_un(
                fmpz_mpolyu_t F, fmpz_mpolyu_t T, const fmpz_mpoly_ctx_t ctx,
                 fmpz_t modulus, nmod_mpolyun_t A, const nmod_mpoly_ctx_t pctx)
{
    int changed = 0;
    slong i, j, k;
    ulong * Texp;
    ulong * Fexp;
    ulong * Aexp;
    slong Flen;
    slong Alen;
    fmpz_mpoly_struct * Tcoeff;
    fmpz_mpoly_struct * Fcoeff;
    nmod_mpolyn_struct  * Acoeff;
    fmpz_mpoly_t zero;
    nmod_mpolyn_t zerop;

    FLINT_ASSERT(F->bits == T->bits);
    FLINT_ASSERT(T->bits == A->bits);

    Flen = F->length;
    Alen = A->length;
    fmpz_mpolyu_fit_length(T, Flen + Alen, ctx);

    Tcoeff = T->coeffs;
    Fcoeff = F->coeffs;
    Acoeff = A->coeffs;
    Texp = T->exps;
    Fexp = F->exps;
    Aexp = A->exps;   

    fmpz_mpoly_init(zero, ctx);
    zero->bits = A->bits;
    zero->length = 0;

    nmod_mpolyn_init(zerop, A->bits, pctx);
    zero->length = 0;

    i = j = k = 0;
    while (i < Flen || j < Alen)
    {
        if (i < Flen && j < Alen && (Fexp[i] == Aexp[j]))
        {
            /* F term ok, A term ok */
            changed |= fmpz_mpoly_addinterp_n(Tcoeff + k, Fcoeff + i, ctx,
                                                    modulus, Acoeff + j, pctx);
            Texp[k] = Aexp[j];
            i++;
            j++;
        }
        else if (i < Flen && (j >= Alen || Fexp[i] > Aexp[j]))
        {
            /* F term ok, A term missing */
            changed |= fmpz_mpoly_addinterp_n(Tcoeff + k, Fcoeff + i, ctx,
                                                         modulus, zerop, pctx);
            Texp[k] = Fexp[i];
            i++;
        }
        else
        {
            FLINT_ASSERT(j < Alen && (i >= Flen || Aexp[j] > Fexp[i]));

            /* F term missing, A term ok */
            changed |= fmpz_mpoly_addinterp_n(Tcoeff + k, zero, ctx,
                                                    modulus, Acoeff + j, pctx);
            Texp[k] = Aexp[j];
            j++;
        }

        FLINT_ASSERT(!fmpz_mpoly_is_zero(Tcoeff + k, ctx));
        k++;
    }
    T->length = k;

    if (changed)
    {
        fmpz_mpolyu_swap(T, F, ctx);
    }

    fmpz_mpoly_clear(zero, ctx);
    nmod_mpolyn_clear(zerop, pctx);

    FLINT_ASSERT(fmpz_mpolyu_is_canonical(F, ctx));

    return changed;    
}




/*
    A and B are assumed to be primitive and therefore are marked const
*/
int fmpz_mpolyu_gcd_brown(fmpz_mpolyu_t G, fmpz_mpolyu_t Abar, fmpz_mpolyu_t Bbar,
      const fmpz_mpolyu_t A, const fmpz_mpolyu_t B, const fmpz_mpoly_ctx_t ctx)
{
    int success;
    fmpz_t bound;
    slong offset, shift;
    mp_limb_t p, t, gammared;
    fmpz_t gamma, modulus;
    fmpz_t gnm, gns, anm, ans, bnm, bns;
    fmpz_t temp;
    fmpz_mpolyu_t T;
    nmod_mpolyun_t Gp, Abarp, Bbarp, Ap, Bp;
    nmod_mpoly_ctx_t pctx;
    mp_bitcnt_t bits = G->bits;
    slong N = mpoly_words_per_exp_sp(bits, ctx->minfo);

    FLINT_ASSERT(bits <= FLINT_BITS);
    FLINT_ASSERT(bits == A->bits);
    FLINT_ASSERT(bits == B->bits);
    FLINT_ASSERT(bits == G->bits);
    FLINT_ASSERT(bits == Abar->bits);
    FLINT_ASSERT(bits == Bbar->bits);

/*
printf("fmpz_mpolyu_gcd_brown called\n");
printf("A: "); fmpz_mpolyu_print_pretty(A, NULL, ctx); printf("\n");
printf("B: "); fmpz_mpolyu_print_pretty(B, NULL, ctx); printf("\n");
*/

    mpoly_gen_offset_shift_sp(&offset, &shift, ctx->minfo->nvars - 1, G->bits, ctx->minfo);

    fmpz_init(gamma);
    fmpz_init(gnm);
    fmpz_init(gns);
    fmpz_init(anm);
    fmpz_init(ans);
    fmpz_init(bnm);
    fmpz_init(bns);
    fmpz_init(bound);
    fmpz_init(temp);
    fmpz_init_set_si(modulus, 1);

#if WANT_ASSERT
    fmpz_mpolyu_content(temp, A, ctx);
    FLINT_ASSERT(fmpz_is_one(temp));
    fmpz_mpolyu_content(temp, B, ctx);
    FLINT_ASSERT(fmpz_is_one(temp));
#endif

    fmpz_gcd(gamma, fmpz_mpolyu_leadcoeff_ref(A),
                    fmpz_mpolyu_leadcoeff_ref(B));

    fmpz_mpolyu_height(bound, A, ctx);
    fmpz_mpolyu_height(temp, B, ctx);
    if (fmpz_cmp(bound, temp) < 0)
        fmpz_swap(bound, temp);
    fmpz_mul(bound, bound, gamma);
    fmpz_add(bound, bound, bound);

    fmpz_mpolyu_init(T, bits, ctx);

    nmod_mpoly_ctx_init(pctx, ctx->minfo->nvars, ORD_LEX, 2);
    nmod_mpolyun_init(Ap, bits, pctx);
    nmod_mpolyun_init(Bp, bits, pctx);
    nmod_mpolyun_init(Gp, bits, pctx);
    nmod_mpolyun_init(Abarp, bits, pctx);
    nmod_mpolyun_init(Bbarp, bits, pctx);

    p = UWORD(1) << (FLINT_BITS - 1);

choose_prime:

    t = p;
    p = n_nextprime(p, 1);
    if (p <= t)
    {
        success = 0;
        goto cleanup;
    }

    /* make sure reduction does not kill both lc(A) and lc(B) */
    gammared = fmpz_fdiv_ui(gamma, p);
    if (gammared == 0)
    {
        goto choose_prime;
    }

    nmod_mpoly_ctx_clear(pctx);
    nmod_mpoly_ctx_init(pctx, ctx->minfo->nvars, ORD_LEX, p);
    /* the stupid nmod poly's store their own context :( */
    nmod_mpolyun_change_nmod_poly_mod(Ap, pctx);
    nmod_mpolyun_change_nmod_poly_mod(Bp, pctx);
    nmod_mpolyun_change_nmod_poly_mod(Gp, pctx);
    nmod_mpolyun_change_nmod_poly_mod(Abarp, pctx);
    nmod_mpolyun_change_nmod_poly_mod(Bbarp, pctx);

    /* reduction should kill neither A nor B */
    fmpz_mpolyu_redto_nmod_mpolyun(Ap, pctx, A, ctx);
    FLINT_ASSERT(Ap->length > 0);
    fmpz_mpolyu_redto_nmod_mpolyun(Bp, pctx, B, ctx);
    FLINT_ASSERT(Bp->length > 0);

    success = nmod_mpolyun_gcd_brown_smprime(Gp, Abarp, Bbarp,
                                          Ap, Bp, ctx->minfo->nvars - 1, pctx);
    if (!success)
    {
        goto choose_prime;
    }

    if (nmod_mpolyun_is_nonzero_nmod(Gp, pctx))
    {
        fmpz_mpolyu_one(G, ctx);
        fmpz_mpolyu_set(Abar, A, ctx);
        fmpz_mpolyu_set(Bbar, B, ctx);
        goto successful_put_content;
    }

    if (!fmpz_is_one(modulus))
    {
        int cmp = 0;
        FLINT_ASSERT(G->length > 0);
        if (G->exps[0] != Gp->exps[0])
        {
            cmp = G->exps[0] > Gp->exps[0] ? 1 : -1;
        }
        if (cmp == 0)
        {
            slong k = nmod_poly_degree((Gp->coeffs + 0)->coeffs + 0);
            cmp = mpoly_monomial_cmp_nomask_extra(
                        (G->coeffs + 0)->exps + N*0,
                       (Gp->coeffs + 0)->exps + N*0, N, offset, k << shift);
        }

        if (cmp < 0)
        {
            goto choose_prime;
        }
        else if (cmp > 0)
        {
            fmpz_one(modulus);
        }
    }

    FLINT_ASSERT(1 == nmod_mpolyun_leadcoeff_last(Gp, pctx));
    nmod_mpolyun_scalar_mul_nmod(Gp, gammared, pctx);

    if (!fmpz_is_one(modulus))
    {
        fmpz_mpolyu_addinterp_un(G, T, ctx, modulus, Gp, pctx);
        fmpz_mpolyu_addinterp_un(Abar, T, ctx, modulus, Abarp, pctx);
        fmpz_mpolyu_addinterp_un(Bbar, T, ctx, modulus, Bbarp, pctx);
    }
    else
    {
        fmpz_mpolyu_startinterp_un(G, ctx, Gp, pctx);
        fmpz_mpolyu_startinterp_un(Abar, ctx, Abarp, pctx);
        fmpz_mpolyu_startinterp_un(Bbar, ctx, Bbarp, pctx);
    }

    fmpz_mul_ui(modulus, modulus, p);

    if (fmpz_cmp(modulus, bound) <= 0)
    {
        goto choose_prime;
    }

    fmpz_mpolyu_heights(gnm, gns, G, ctx);
    fmpz_mpolyu_heights(anm, ans, Abar, ctx);
    fmpz_mpolyu_heights(bnm, bns, Bbar, ctx);
    fmpz_mul(ans, ans, gnm);
    fmpz_mul(anm, anm, gns);
    fmpz_mul(bns, bns, gnm);
    fmpz_mul(bnm, bnm, gns);

    if (fmpz_cmp(ans, anm) > 0)
        fmpz_swap(ans, anm);
    if (fmpz_cmp(bns, bnm) > 0)
        fmpz_swap(bns, bnm);
    fmpz_add(ans, ans, ans);
    fmpz_add(bns, bns, bns);
    if (fmpz_cmp(ans, modulus) < 0 && fmpz_cmp(bns, modulus) < 0)
    {
        goto successful;
    }

    /* do not reset modulus to 1 */
    goto choose_prime;

successful:

    fmpz_mpolyu_content(temp, G, ctx);
    fmpz_mpolyu_scalar_divexact(G, temp, ctx);
    fmpz_mpolyu_scalar_divexact(Abar, fmpz_mpolyu_leadcoeff_ref(G), ctx);
    fmpz_mpolyu_scalar_divexact(Bbar, fmpz_mpolyu_leadcoeff_ref(G), ctx);

successful_put_content:

    /* inputs were supposed to be primitive - nothing to do */
    success = 1;

cleanup:

#if WANT_ASSERT
    if (success)
    {
        fmpz_mul(temp, fmpz_mpolyu_leadcoeff_ref(G), fmpz_mpolyu_leadcoeff_ref(Abar));
        FLINT_ASSERT(fmpz_equal(temp, fmpz_mpolyu_leadcoeff_ref(A)));
        fmpz_mul(temp, fmpz_mpolyu_leadcoeff_ref(G), fmpz_mpolyu_leadcoeff_ref(Bbar));
        FLINT_ASSERT(fmpz_equal(temp, fmpz_mpolyu_leadcoeff_ref(B)));
    }
#endif

    fmpz_clear(gamma);
    fmpz_clear(gnm);
    fmpz_clear(gns);
    fmpz_clear(anm);
    fmpz_clear(ans);
    fmpz_clear(bnm);
    fmpz_clear(bns);
    fmpz_clear(bound);
    fmpz_clear(temp);
    fmpz_clear(modulus);

    nmod_mpolyun_clear(Gp, pctx);
    nmod_mpolyun_clear(Abarp, pctx);
    nmod_mpolyun_clear(Bbarp, pctx);
    nmod_mpolyun_clear(Ap, pctx);
    nmod_mpolyun_clear(Bp, pctx);

    nmod_mpoly_ctx_clear(pctx);

    fmpz_mpolyu_clear(T, ctx);

    return success;
}






/*
typedef struct
{
    volatile int idx;
    volatile slong G_exp, Abar_exp, Bbar_exp;
    pthread_mutex_t mutex;
    const nmod_mpoly_ctx_struct * ctx;
    nmod_poly_crt_struct * CRT;
    nmod_mpolyun_struct ** gptrs, ** abarptrs, ** bbarptrs;
    ulong numthreads;
}
_joinbase_struct;

typedef _joinbase_struct _joinbase_t[1];

typedef struct
{
    _joinbase_struct * base;
    nmod_mpolyun_t G, Abar, Bbar;
    slong G_lastdeg, Abar_lastdeg, Bbar_lastdeg;
}
_joinworker_arg_struct;
*/





/*
typedef struct
{
    volatile int gcd_is_one;
    pthread_mutex_t mutex;
    fmpz_t gamma;
    mp_limb_t p;
    const fmpz_mpoly_ctx_struct * ctx;
    fmpz_mpolyu_struct * A, * B;
}
_splitbase_struct;

typedef _splitbase_struct _splitbase_t[1];


typedef struct
{
    slong idx;
    _splitbase_struct * base;
    mp_limb_t p;
    nmod_mpolyu_t G;
    nmod_mpolyu_t Abar;
    nmod_mpolyu_t Bbar;
}
_splitworker_arg_struct;


typedef struct
{
    mp_limb_t p;
    nmod_mpolyu_t G;
    nmod_mpolyu_t Abar;
    nmod_mpolyu_t Bbar;
}
_parray_elem


typedef struct
{
    slong length;
    slong alloc;
    _parray_elem_struct * array;
}
_parray_struct

typedef _parray_struct _parray_t[1];

void _parray_init(_parray_t L)
{
    L->length = 0;
    L->alloc = 0;
    L->array = NULL;
}

void _parray_fit_length(_parray_t L, slong k)
{
    k = FLINT_MAX(WORD(1), k);

    if (L->alloc == 0)
    {
        FLINT_ASSERT(P->array == NULL);
        L->array = (_parray_elem *) flint_malloc(k*sizeof(_parray_elem));
        L->alloc = k;
    }
    else if (k > P->alloc)
    {
        FLINT_ASSERT(P->array != NULL);
        P->prog = (_parray_elem *) flint_realloc(P->prog, k*sizeof(_parray_elem));
        P->alloc = k;
    }
}

void _parray_push_back(_parray_t L, _splitworker_arg_struct * arg)
{
    _parray_elem_struct * top;
    _parray_fit_length(L, L->length + 1);

    top = L->array + L->length;
    top->p = arg->p;
    top->G = arg->G;
    top->Abar = arg->Abar;
    top->Bbar = arg->Bbar;
    L->length++;

    nmod_mpolyu_init(arg->G, NULL);
    nmod_mpolyu_init(arg->Abar, NULL);
    nmod_mpolyu_init(arg->Bbar, NULL);
}

void _parray_calculate_modulus(const _parray_t L, fmpz_t modulus)
{
    fmpz_one(modulus);
    for (i = 0; i < L->length; i++)
    {
        fmpz_mul_ui(modulus, modulus, L->array[i].p))
    }
}

void _parray_clear(_parray_t L)
{
    slong i;
    for (i = 0; i < L->length; i++)
    {
        nmod_mpolyu_clear(L->array[i].G, NULL);
        nmod_mpolyu_clear(L->array[i].Abar, NULL);
        nmod_mpolyu_clear(L->array[i].Bbar, NULL);
    }
    L->length = 0;
}
*/







typedef struct
{
    volatile int gcd_is_one;
    volatile mp_limb_t p;
    pthread_mutex_t mutex;
    fmpz_t gamma;
    const fmpz_mpoly_ctx_struct * ctx;
    fmpz_mpolyu_struct * A, * B;
    ulong numthreads;
    slong var;
}
_splitbase_struct;

typedef _splitbase_struct _splitbase_t[1];

typedef struct
{
    slong idx;
    _splitbase_struct * base;
    fmpz_mpolyu_t G, Abar, Bbar;
    fmpz_t modulus;
    mp_bitcnt_t required_bits;
}
_splitworker_arg_struct;


static void _splitworker(void * varg)
{
    _splitworker_arg_struct * arg = (_splitworker_arg_struct *) varg;
    _splitbase_struct * base = arg->base;
    const fmpz_mpoly_ctx_struct * ctx = base->ctx;
    mp_bitcnt_t bits = base->A->bits;
    slong N = mpoly_words_per_exp_sp(bits, ctx->minfo);
    slong offset, shift;
    nmod_mpolyun_t Ap, Bp, Gp, Abarp, Bbarp;
    fmpz_mpolyu_t T;
    int success;
    mp_limb_t p, gammared;
    nmod_mpoly_ctx_t pctx;

pthread_mutex_lock(&iomutex);
printf("splitworker started\n");
pthread_mutex_unlock(&iomutex);

    mpoly_gen_offset_shift_sp(&offset, &shift, ctx->minfo->nvars - 1, bits, ctx->minfo);

    fmpz_mpolyu_init(T, bits, ctx);

    nmod_mpoly_ctx_init(pctx, ctx->minfo->nvars, ORD_LEX, 2);
    nmod_mpolyun_init(Ap, bits, pctx);
    nmod_mpolyun_init(Bp, bits, pctx);
    nmod_mpolyun_init(Gp, bits, pctx);
    nmod_mpolyun_init(Abarp, bits, pctx);
    nmod_mpolyun_init(Bbarp, bits, pctx);

    while (fmpz_bits(arg->modulus) <= arg->required_bits)
    {
        /* get prime */
        pthread_mutex_lock(&base->mutex);
        p = n_nextprime(base->p, 1);
        base->p = p;        
        pthread_mutex_unlock(&base->mutex);


pthread_mutex_lock(&iomutex);
flint_printf("splitworker prime %wu, gamma: ", p); fmpz_print(base->gamma); printf("\n");
pthread_mutex_unlock(&iomutex);


        /* make sure reduction does not kill both lc(A) and lc(B) */
        gammared = fmpz_fdiv_ui(base->gamma, p);
        if (gammared == 0)
        {
            continue;
        }


pthread_mutex_lock(&iomutex);
flint_printf("splitworker gammared %wu\n", gammared);
pthread_mutex_unlock(&iomutex);


        nmod_mpoly_ctx_clear(pctx);
        nmod_mpoly_ctx_init(pctx, ctx->minfo->nvars, ORD_LEX, p);


pthread_mutex_lock(&iomutex);
flint_printf("splitworker changing mods\n");
pthread_mutex_unlock(&iomutex);


        /* the stupid nmod poly's store their own context :( */
        nmod_mpolyun_change_nmod_poly_mod(Ap, pctx);
        nmod_mpolyun_change_nmod_poly_mod(Bp, pctx);
        nmod_mpolyun_change_nmod_poly_mod(Gp, pctx);
        nmod_mpolyun_change_nmod_poly_mod(Abarp, pctx);
        nmod_mpolyun_change_nmod_poly_mod(Bbarp, pctx);

pthread_mutex_lock(&iomutex);
flint_printf("splitworker reducing now\n");
pthread_mutex_unlock(&iomutex);


        /* reduction should kill neither A nor B */
        fmpz_mpolyu_redto_nmod_mpolyun(Ap, pctx, base->A, ctx);
        FLINT_ASSERT(Ap->length > 0);
        fmpz_mpolyu_redto_nmod_mpolyun(Bp, pctx, base->B, ctx);
        FLINT_ASSERT(Bp->length > 0);


pthread_mutex_lock(&iomutex);
flint_printf("splitworker Ap:"); nmod_mpolyun_print_pretty(Ap, NULL, pctx); printf("\n");
flint_printf("splitworker Bp:"); nmod_mpolyun_print_pretty(Bp, NULL, pctx); printf("\n");
pthread_mutex_unlock(&iomutex);


        success = nmod_mpolyun_gcd_brown_smprime(Gp, Abarp, Bbarp,
                                              Ap, Bp, ctx->minfo->nvars - 1, pctx);
        if (!success)
        {
            continue;
        }

        FLINT_ASSERT(Gp->length > 0);
        FLINT_ASSERT(Abarp->length > 0);
        FLINT_ASSERT(Bbarp->length > 0);

        /* check up */
        if (base->gcd_is_one)
        {
FLINT_ASSERT(0);
            break;
        }

        if (nmod_mpolyun_is_nonzero_nmod(Gp, pctx))
        {
            base->gcd_is_one = 1;
FLINT_ASSERT(0);
            break;
        }

        if (!fmpz_is_one(arg->modulus))
        {
            int cmp = 0;
            FLINT_ASSERT(arg->G->length > 0);
            if (arg->G->exps[0] != Gp->exps[0])
            {
                cmp = arg->G->exps[0] > Gp->exps[0] ? 1 : -1;
            }
            if (cmp == 0)
            {
                slong k = nmod_poly_degree((Gp->coeffs + 0)->coeffs + 0);
                cmp = mpoly_monomial_cmp_nomask_extra(
                       (arg->G->coeffs + 0)->exps + N*0,
                           (Gp->coeffs + 0)->exps + N*0, N, offset, k << shift);
            }

            if (cmp < 0)
            {
                continue;
            }
            else if (cmp > 0)
            {
                fmpz_one(arg->modulus);
            }
        }

        FLINT_ASSERT(1 == nmod_mpolyun_leadcoeff_last(Gp, pctx));
        nmod_mpolyun_scalar_mul_nmod(Gp, gammared, pctx);

        if (!fmpz_is_one(arg->modulus))
        {
            fmpz_mpolyu_addinterp_un(arg->G, T, ctx, arg->modulus, Gp, pctx);
            fmpz_mpolyu_addinterp_un(arg->Abar, T, ctx, arg->modulus, Abarp, pctx);
            fmpz_mpolyu_addinterp_un(arg->Bbar, T, ctx, arg->modulus, Bbarp, pctx);
        }
        else
        {
            fmpz_mpolyu_startinterp_un(arg->G, ctx, Gp, pctx);
            fmpz_mpolyu_startinterp_un(arg->Abar, ctx, Abarp, pctx);
            fmpz_mpolyu_startinterp_un(arg->Bbar, ctx, Bbarp, pctx);
        }

        fmpz_mul_ui(arg->modulus, arg->modulus, p);
printf("arg->modulus: "); fmpz_print(arg->modulus); printf("\n");

    }


flint_printf("modulus bits %wd, required bits %wd\n", fmpz_bits(arg->modulus), arg->required_bits);


    fmpz_mpolyu_clear(T, ctx);

    nmod_mpolyun_clear(Ap, pctx);
    nmod_mpolyun_clear(Bp, pctx);
    nmod_mpolyun_clear(Gp, pctx);
    nmod_mpolyun_clear(Abarp, pctx);
    nmod_mpolyun_clear(Bbarp, pctx);

    nmod_mpoly_ctx_clear(pctx);
}










/* A = crt(B[0], ...., B[count-1]) wrt to P */
void fmpz_mpoly_crt(
    const fmpz_crt_t P,
    fmpz_t Amax, fmpz_t Asum,
    fmpz_mpoly_t A,
    fmpz_mpoly_struct * const * B,
    slong count,
    const fmpz_mpoly_ctx_t ctx)
{
    int cmp;
    slong N = mpoly_words_per_exp_sp(A->bits, ctx->minfo);
    fmpz ** input;
    slong * start;
    slong Ai;
    slong j, k;
    fmpz_t zero;
    TMP_INIT;

    TMP_START;

    fmpz_init(zero);

    input = (fmpz **) TMP_ALLOC(count * sizeof(fmpz *));
    start = (slong *) TMP_ALLOC(count * sizeof(slong));

    /* start[k] is the next available term in B[k] */
    for (k = 0; k < count; k++)
    {
        start[k] = 0;    
    }

    Ai = 0;
    while (1)
    {
        fmpz_mpoly_fit_length(A, Ai + 1, ctx);

        k = 0;
        do
        {
            input[k] = zero;
            if (start[k] < B[k]->length)
            {
                goto found_max;
            }
        } while (++k < count);

        break; /* all B[k] have been scanned completely */

    found_max:

        input[k] = B[k]->coeffs + start[k];
        mpoly_monomial_set(A->exps + N*Ai, B[k]->exps + N*start[k], N);
        start[k]++;

        for (k++; k < count; k++)
        {
            input[k] = zero;
            if (start[k] >= B[k]->length)
            {
                continue;
            }

            cmp = mpoly_monomial_cmp_nomask(B[k]->exps + N*start[k], A->exps + N*Ai, N);
            if (cmp == 0)
            {
                input[k] = B[k]->coeffs + start[k];
                start[k]++;
            }
            else if (cmp > 0)
            {
                /* undo previous max's */
                for (j = 0; j < k; j++)
                {
                    start[j] -= (input[j] != zero);
                    input[j] = zero;
                }
                goto found_max;
            }
        }

        fmpz_crt_run(P, A->coeffs + Ai, input);
pthread_mutex_lock(&iomutex);
printf("result: "); fmpz_print(A->coeffs + Ai); printf("\n");
pthread_mutex_unlock(&iomutex);

        if (fmpz_sgn(A->coeffs + Ai) > 0)
        {
            fmpz_add(Asum, Asum, A->coeffs + Ai);
        }
        else
        {
            fmpz_sub(Asum, Asum, A->coeffs + Ai);
        }


        if (fmpz_cmpabs(Amax, A->coeffs + Ai) < 0)
        {
            fmpz_set(Amax, A->coeffs + Ai);
            fmpz_abs(Amax, Amax);
        }

        Ai += !fmpz_is_zero(A->coeffs + Ai);

    }
    A->length = Ai;

    TMP_END;

    fmpz_clear(zero);
    return;
}

/*
    Append to A the result of crt'ing the coeff of X^exp
    Amax = max(Amax, abs(coeff0), abs(coeff1), ...)
    Asum = Asum + abs(coeffs0) + abs(coeffs1) + ...
*/
void fmpz_mpolyu_crt_exp(
    const fmpz_crt_t P,
    fmpz_t Amax, fmpz_t Asum,
    fmpz_mpolyu_t A,
    ulong exp,
    fmpz_mpolyu_struct * const * B,
    slong count,
    const fmpz_mpoly_ctx_t ctx)
{
    slong j, k;
    slong Ai;
    fmpz_mpoly_struct ** C;
    fmpz_mpoly_t zero;
    TMP_INIT;

    fmpz_mpoly_init(zero, ctx);
    zero->length = 0;

    TMP_START;
    C = (fmpz_mpoly_struct **) TMP_ALLOC(count * sizeof(fmpz_mpoly_struct *));
    for (k = 0; k < count; k++)
    {
        C[k] = zero;
        for (j = 0; j < B[k]->length; j++)
        {
            if (B[k]->exps[j] == exp)
            {
                C[k] = B[k]->coeffs + j;
                break;
            }
        }
    }

    Ai = A->length;
    fmpz_mpolyu_fit_length(A, Ai + 1, ctx);
    A->exps[Ai] = exp;
    fmpz_mpoly_crt(P, Amax, Asum, A->coeffs + Ai, C, count, ctx);
    A->length += (A->coeffs + Ai)->length != 0;

    TMP_END;
    fmpz_mpoly_clear(zero, ctx);
    return;
}



typedef struct
{
    volatile int idx;
    volatile slong G_exp, Abar_exp, Bbar_exp;
    pthread_mutex_t mutex;
    const fmpz_mpoly_ctx_struct * ctx;
    fmpz_crt_struct * CRT;
    fmpz_mpolyu_struct ** gptrs, ** abarptrs, ** bbarptrs;
    ulong numthreads;
}
_joinbase_struct;

typedef _joinbase_struct _joinbase_t[1];

typedef struct
{
    _joinbase_struct * base;
    fmpz_mpolyu_t G, Abar, Bbar;
    fmpz_t Gmax, Gsum, Abarmax, Abarsum, Bbarmax, Bbarsum;
}
_joinworker_arg_struct;

static void _joinworker(void * varg)
{
    _joinworker_arg_struct * arg = (_joinworker_arg_struct *) varg;
    _joinbase_struct * base = arg->base;
    slong our_G_exp, our_Abar_exp, our_Bbar_exp;

    while (1)
    {
        /* get exponent of either G, Abar, or Bbar to start working on */
        pthread_mutex_lock(&base->mutex);
        our_G_exp = base->G_exp;
        our_Abar_exp = base->Abar_exp;
        our_Bbar_exp = base->Bbar_exp;
        if (our_G_exp >= 0)
        {
            base->G_exp = our_G_exp - 1;
        }
        else if (our_Abar_exp >= 0)
        {
            base->Abar_exp = our_Abar_exp - 1;            
        }
        else if (our_Bbar_exp >= 0)
        {
            base->Bbar_exp = our_Bbar_exp - 1;            
        }
        pthread_mutex_unlock(&base->mutex);

        if (our_G_exp >= 0)
        {
            fmpz_mpolyu_crt_exp(base->CRT, arg->Gmax, arg->Gsum, arg->G, our_G_exp,
                                     base->gptrs, base->numthreads, base->ctx);
        }
        else if (our_Abar_exp >= 0)
        {
            fmpz_mpolyu_crt_exp(base->CRT, arg->Abarmax, arg->Abarsum, arg->Abar, our_Abar_exp,
                                  base->abarptrs, base->numthreads, base->ctx);
        }
        else if (our_Bbar_exp >= 0)
        {
            fmpz_mpolyu_crt_exp(base->CRT, arg->Bbarmax, arg->Bbarsum, arg->Bbar, our_Bbar_exp,
                                  base->bbarptrs, base->numthreads, base->ctx);
        }
        else
        {
            return;
        }
    }
}



int fmpz_mpolyu_gcd_brown_threaded(fmpz_mpolyu_t G, fmpz_mpolyu_t Abar,
                     fmpz_mpolyu_t Bbar, fmpz_mpolyu_t A, fmpz_mpolyu_t B,
                                const fmpz_mpoly_ctx_t ctx, slong thread_limit)
{
    slong i;
    mp_bitcnt_t bits = A->bits;
    slong N = mpoly_words_per_exp_sp(bits, ctx->minfo);
    ulong primes_needed;
    ulong numthreads, split_numthreads;
    int success;
    fmpz_t bound, modulus, temp;
    fmpz_t gnm, gns, anm, ans, bnm, bns;
    slong Gexp0, Abarexp0, Bbarexp0;
    fmpz ** mptrs;
    fmpz_mpolyu_struct ** gptrs, ** abarptrs, ** bbarptrs;
    fmpz_crt_t P;
    _splitworker_arg_struct * splitargs;
    _splitbase_t splitbase;
    _joinworker_arg_struct * joinargs;
    _joinbase_t joinbase;
    slong max_numworkers, numworkers;
    thread_pool_handle * handles;
    slong * starts;
    slong Gi;

    FLINT_ASSERT(global_thread_pool_initialized);
    max_numworkers = thread_pool_get_size(global_thread_pool);
    thread_limit = FLINT_MIN(thread_limit, max_numworkers + 1);
    if (thread_limit < 2)
    {
        return fmpz_mpolyu_gcd_brown(G, Abar, Bbar, A, B, ctx);
    }

    return 0;

    handles = (thread_pool_handle *) flint_malloc((thread_limit - 1)*sizeof(thread_pool_handle));

printf("gcd brown threaded called\n");
printf("A: "); fmpz_mpolyu_print_pretty(A, NULL, ctx); printf("\n");
printf("B: "); fmpz_mpolyu_print_pretty(B, NULL, ctx); printf("\n");


pthread_mutex_init(&iomutex, NULL);

    fmpz_init(gnm);
    fmpz_init(gns);
    fmpz_init(anm);
    fmpz_init(ans);
    fmpz_init(bnm);
    fmpz_init(bns);
    fmpz_init(bound);
    fmpz_init(temp);
    fmpz_init(modulus);

#if WANT_ASSERT
    fmpz_mpolyu_content(temp, A, ctx);
    FLINT_ASSERT(fmpz_is_one(temp));
    fmpz_mpolyu_content(temp, B, ctx);
    FLINT_ASSERT(fmpz_is_one(temp));
#endif

    fmpz_init(splitbase->gamma);
    fmpz_gcd(splitbase->gamma, fmpz_mpolyu_leadcoeff_ref(A),
                               fmpz_mpolyu_leadcoeff_ref(B));

    fmpz_mpolyu_height(bound, A, ctx);
    fmpz_mpolyu_height(temp, B, ctx);

    if (fmpz_cmp(bound, temp) < 0)
        fmpz_swap(bound, temp);
    fmpz_mul(bound, bound, splitbase->gamma);
    fmpz_add(bound, bound, bound);
printf("bound: "); fmpz_print(bound); printf("\n");


    splitbase->p = n_nextprime(UWORD(1) << (FLINT_BITS - 2), 1);


flint_printf("primes_needed: %wu\n", primes_needed);

    handles = (thread_pool_handle *) flint_malloc((thread_limit - 1)*sizeof(thread_pool_handle));
    numworkers = thread_pool_request(global_thread_pool, handles, thread_limit - 1);
    numthreads = numworkers + 1;

    primes_needed = (fmpz_bits(bound) + FLINT_BITS)/((ulong)(FLINT_BITS - 2));
    FLINT_ASSERT(primes_needed > 0);
    split_numthreads = FLINT_MIN(primes_needed, (ulong)(thread_limit));


    gptrs = (fmpz_mpolyu_struct **) flint_malloc(numthreads*sizeof(fmpz_mpolyu_struct *));
    abarptrs = (fmpz_mpolyu_struct **) flint_malloc(numthreads*sizeof(fmpz_mpolyu_struct *));
    bbarptrs = (fmpz_mpolyu_struct **) flint_malloc(numthreads*sizeof(fmpz_mpolyu_struct *));
    mptrs = (fmpz **) flint_malloc(numthreads*sizeof(fmpz *));
    splitargs = (_splitworker_arg_struct *) flint_malloc(numthreads*sizeof(_splitworker_arg_struct));
    for (i = 0; i < numthreads; i++)
    {
        fmpz_mpolyu_init(splitargs[i].G, bits, ctx);
        fmpz_mpolyu_init(splitargs[i].Abar, bits, ctx);
        fmpz_mpolyu_init(splitargs[i].Bbar, bits, ctx);
        fmpz_init_set_ui(splitargs[i].modulus, 1);
    }
    splitbase->numthreads = numthreads;
    splitbase->A = A;
    splitbase->B = B;
    splitbase->ctx = ctx;

    pthread_mutex_init(&splitbase->mutex, NULL);

compute_split:

    splitbase->gcd_is_one = 0;
    for (i = 0; i < numthreads; i++)
    {
        ulong ri = (fmpz_bits(bound) + numthreads)/numthreads;
        splitargs[i].idx = i;
        splitargs[i].base = splitbase;
        splitargs[i].required_bits = FLINT_MAX(ri, fmpz_bits(splitargs[i].modulus) + 1);

flint_printf("required bits[%wd]: %wd\n", i, splitargs[i].required_bits);

    }

    for (i = 0; i + 1 < numthreads; i++)
    {
        thread_pool_wake(global_thread_pool, handles[i], _splitworker, &splitargs[i]);
    }
    _splitworker(&splitargs[numthreads - 1]);
    for (i = 0; i + 1 < numthreads; i++)
    {
        thread_pool_wait(global_thread_pool, handles[i]);
    }

    if (splitbase->gcd_is_one)
    {
        FLINT_ASSERT(0);
    }

    fmpz_one(modulus);
    for (i = 0; i < numthreads; i++)
    {
        gptrs[i] = splitargs[i].G;
        abarptrs[i] = splitargs[i].Abar;
        bbarptrs[i] = splitargs[i].Bbar;
        mptrs[i] = splitargs[i].modulus;
        if (fmpz_bits(splitargs[i].modulus) <= splitargs[i].required_bits)
        {
            /* not enough primes - must fail (should not happen in practise) */
            FLINT_ASSERT(0);
        }
        FLINT_ASSERT(gptrs[i]->length > 0);
        FLINT_ASSERT(abarptrs[i]->length > 0);
        FLINT_ASSERT(bbarptrs[i]->length > 0);
        fmpz_mul(modulus, modulus, mptrs[i]);
    }

    /*
        Check for consistency in the leading monomial. All args have at least
        one image, so G, Abar, Bbar are defined and nonzero for each.
    */
    Gexp0 = gptrs[0]->exps[0];
    Abarexp0 = A->exps[0] - Gexp0;
    Bbarexp0 = B->exps[0] - Gexp0;
    for (i = 0; i < numthreads; i++)
    {
        if (gptrs[i]->exps[0] != Gexp0
                      || !mpoly_monomial_equal(
                                  (gptrs[i]->coeffs + 0)->exps + N*0,
                                  (gptrs[0]->coeffs + 0)->exps + N*0, N))
        {
            /* very unlucky - could try again or just fail */
            FLINT_ASSERT(0);
        }
        FLINT_ASSERT(splitargs[i].Abar->exps[0] == Abarexp0);
        FLINT_ASSERT(splitargs[i].Bbar->exps[0] == Bbarexp0);   


flint_printf("arg %wd:\n", i);
flint_printf("modulus: "); fmpz_print(mptrs[i]); printf("\n");
flint_printf("      G: "); fmpz_mpolyu_print_pretty(gptrs[i], NULL, ctx); printf("\n");
flint_printf("   Abar: "); fmpz_mpolyu_print_pretty(abarptrs[i], NULL, ctx); printf("\n");
flint_printf("   Bbar: "); fmpz_mpolyu_print_pretty(bbarptrs[i], NULL, ctx); printf("\n");

    }

    FLINT_ASSERT(fmpz_cmp(modulus, bound) > 0);

    fmpz_crt_init(P);
    success = fmpz_crt_compile(P, mptrs, numthreads);
    FLINT_ASSERT(success);

    joinbase->numthreads = numthreads;
    joinbase->gptrs = gptrs;
    joinbase->abarptrs = abarptrs;
    joinbase->bbarptrs = bbarptrs;
    joinbase->G_exp = Gexp0;
    joinbase->Abar_exp = Abarexp0;
    joinbase->Bbar_exp = Bbarexp0;
    joinbase->ctx = ctx;
    joinbase->CRT = P;
    pthread_mutex_init(&joinbase->mutex, NULL);

    joinargs = (_joinworker_arg_struct *) flint_malloc(numthreads*sizeof(_joinworker_arg_struct));

    for (i = 0; i < numthreads; i++)
    {
        joinargs[i].base = joinbase;
        fmpz_mpolyu_init(joinargs[i].G, bits, ctx);
        fmpz_mpolyu_init(joinargs[i].Abar, bits, ctx);
        fmpz_mpolyu_init(joinargs[i].Bbar, bits, ctx);
        fmpz_init(joinargs[i].Gmax);
        fmpz_init(joinargs[i].Gsum);
        fmpz_init(joinargs[i].Abarmax);
        fmpz_init(joinargs[i].Abarsum);
        fmpz_init(joinargs[i].Bbarmax);
        fmpz_init(joinargs[i].Bbarsum);
    }
    for (i = 0; i + 1 < numthreads; i++)
    {
        thread_pool_wake(global_thread_pool, handles[i], _joinworker, joinargs + i);
    }
    _joinworker(joinargs + numthreads - 1);
    for (i = 0; i + 1 < numthreads; i++)
    {
        thread_pool_wait(global_thread_pool, handles[i]);
    }
    pthread_mutex_destroy(&joinbase->mutex);

printf("joining G\n");
    for (i = 0; i < numthreads; i++)
    {
flint_printf("   arg[%wd].G:", i); fmpz_mpolyu_print_pretty(joinargs[i].G, NULL, ctx); printf("\n");
flint_printf("arg[%wd].Abar:", i); fmpz_mpolyu_print_pretty(joinargs[i].Abar, NULL, ctx); printf("\n");
flint_printf("arg[%wd].Bbar:", i); fmpz_mpolyu_print_pretty(joinargs[i].Bbar, NULL, ctx); printf("\n");
flint_printf("   G max[%wd]: ",i); fmpz_print(joinargs[i].Gmax); printf("\n");
flint_printf("   G sum[%wd]: ",i); fmpz_print(joinargs[i].Gsum); printf("\n");
flint_printf("Abar max[%wd]: ",i); fmpz_print(joinargs[i].Abarmax); printf("\n");
flint_printf("Abar sum[%wd]: ",i); fmpz_print(joinargs[i].Abarsum); printf("\n");
flint_printf("Bbar max[%wd]: ",i); fmpz_print(joinargs[i].Bbarmax); printf("\n");
flint_printf("Bbar sum[%wd]: ",i); fmpz_print(joinargs[i].Bbarsum); printf("\n");
    }


    fmpz_zero(gnm);
    fmpz_zero(gns);
    fmpz_zero(anm);
    fmpz_zero(ans);
    fmpz_zero(bnm);
    fmpz_zero(bns);
    for (i = 0; i < numthreads; i++)
    {
        FLINT_ASSERT(fmpz_sgn(joinargs[i].Gmax) >= 0);
        FLINT_ASSERT(fmpz_sgn(joinargs[i].Gsum) >= 0);
        if (fmpz_cmp(gnm, joinargs[i].Gmax) < 0)
            fmpz_set(gnm, joinargs[i].Gmax);
        fmpz_add(gns, gns, joinargs[i].Gsum);

        FLINT_ASSERT(fmpz_sgn(joinargs[i].Abarmax) >= 0);
        FLINT_ASSERT(fmpz_sgn(joinargs[i].Abarsum) >= 0);
        if (fmpz_cmp(anm, joinargs[i].Abarmax) < 0)
            fmpz_set(anm, joinargs[i].Abarmax);
        fmpz_add(ans, ans, joinargs[i].Abarsum);

        FLINT_ASSERT(fmpz_sgn(joinargs[i].Bbarmax) >= 0);
        FLINT_ASSERT(fmpz_sgn(joinargs[i].Bbarsum) >= 0);
        if (fmpz_cmp(bnm, joinargs[i].Bbarmax) < 0)
            fmpz_set(bnm, joinargs[i].Bbarmax);
        fmpz_add(bns, bns, joinargs[i].Bbarsum);
    }


printf("   G max: "); fmpz_print(gnm); printf("\n");
printf("   G sum: "); fmpz_print(gns); printf("\n");
printf("Abar max: "); fmpz_print(anm); printf("\n");
printf("Abar sum: "); fmpz_print(ans); printf("\n");
printf("Bbar max: "); fmpz_print(bnm); printf("\n");
printf("Bbar sum: "); fmpz_print(bns); printf("\n");

    /* quickndirty joiner of G */
    starts = (slong *) flint_malloc(numthreads*sizeof(slong));
    for (i = 0; i < numthreads; i++)
    {
        starts[i] = 0;
    }
    Gi = 0;
    while (1)
    {
        slong max_pos = -WORD(1);
        slong max_exp = -WORD(1);
        for (i = 0; i < numthreads; i++)
        {
            if (starts[i] < joinargs[i].G->length
                          && (slong)(joinargs[i].G->exps[starts[i]]) > max_exp)
            {
                max_pos = i;
                max_exp = joinargs[i].G->exps[starts[i]];
            }
        }
        if (max_pos < 0)
        {
            break;
        }
        fmpz_mpolyu_fit_length(G, Gi + 1, ctx);
        G->exps[Gi] = max_exp;
        fmpz_mpoly_swap(G->coeffs + Gi, joinargs[max_pos].G->coeffs + starts[max_pos], ctx);
        starts[max_pos]++;
        Gi++;
    }
    G->length = Gi;
    flint_free(starts);

    /* free join data */
    fmpz_crt_clear(P);
    for (i = 0; i < numthreads; i++)
    {
        fmpz_mpolyu_clear(joinargs[i].G, ctx);
        fmpz_mpolyu_clear(joinargs[i].Abar, ctx);
        fmpz_mpolyu_clear(joinargs[i].Bbar, ctx);
    }
    flint_free(joinargs);


    fmpz_mul(ans, ans, gnm);
    fmpz_mul(anm, anm, gns);
    fmpz_mul(bns, bns, gnm);
    fmpz_mul(bnm, bnm, gns);

    if (fmpz_cmp(ans, anm) > 0)
        fmpz_swap(ans, anm);
    if (fmpz_cmp(bns, bnm) > 0)
        fmpz_swap(bns, bnm);
    fmpz_add(ans, ans, ans);
    fmpz_add(bns, bns, bns);


printf("    ans: "); fmpz_print(ans); printf("\n");
printf("    bns: "); fmpz_print(bns); printf("\n");
printf("modulus: "); fmpz_print(modulus); printf("\n");

    if (fmpz_cmp(ans, modulus) >= 0 || fmpz_cmp(bns, modulus) >= 0)
    {
        FLINT_ASSERT(0);
    }

    success = 1;

cleanup_split:

    pthread_mutex_destroy(&splitbase->mutex);
    fmpz_clear(splitbase->gamma);

    for (i = 0; i < numthreads; i++)
    {
        fmpz_mpolyu_clear(splitargs[i].G, ctx);
        fmpz_mpolyu_clear(splitargs[i].Abar, ctx);
        fmpz_mpolyu_clear(splitargs[i].Bbar, ctx);
        fmpz_clear(splitargs[i].modulus);
    }
    for (i = 0; i + 1 < numthreads; i++)
    {
        thread_pool_give_back(global_thread_pool, handles[i]);
    }
    flint_free(handles);
    flint_free(gptrs);
    flint_free(abarptrs);
    flint_free(bbarptrs);
    flint_free(mptrs);
    flint_free(splitargs);


pthread_mutex_destroy(&iomutex);

    return 0;
}





void fmpz_mpoly_to_fmpz_poly_keepbits(fmpz_poly_t A, slong * Ashift,
               const fmpz_mpoly_t B, slong var, const fmpz_mpoly_ctx_t ctx);

void fmpz_mpoly_from_fmpz_poly_keepbits(fmpz_mpoly_t A, const fmpz_poly_t B,
                           slong Bshift, slong var, mp_bitcnt_t bits, const fmpz_mpoly_ctx_t ctx);

int fmpz_mpoly_gcd_brownnew(fmpz_mpoly_t G,
                              const fmpz_mpoly_t A, const fmpz_mpoly_t B,
                                                    const fmpz_mpoly_ctx_t ctx)
{
    int success;
    slong * perm;
    ulong * shift, * stride;
    slong i;
    mp_bitcnt_t new_bits;
    fmpz_mpoly_ctx_t uctx;
    fmpz_mpolyu_t Au, Bu, Gu, Abaru, Bbaru;
    timeit_t time;
    fmpz_t cA, cB, cG;

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

    if (ctx->minfo->nvars == 1)
    {
        slong shiftA, shiftB;
        fmpz_poly_t a, b, g;
        fmpz_poly_init(a);
        fmpz_poly_init(b);
        fmpz_poly_init(g);
        fmpz_mpoly_to_fmpz_poly_keepbits(a, &shiftA, A, 0, ctx);
        fmpz_mpoly_to_fmpz_poly_keepbits(b, &shiftB, B, 0, ctx);
        fmpz_poly_gcd(g, a, b);
        fmpz_mpoly_from_fmpz_poly_keepbits(G, g, FLINT_MIN(shiftA, shiftB), 0, A->bits, ctx);
        fmpz_poly_clear(a);
        fmpz_poly_clear(b);
        fmpz_poly_clear(g);
        return 1;
    }

    perm = (slong *) flint_malloc(ctx->minfo->nvars*sizeof(slong));
    shift = (ulong *) flint_malloc(ctx->minfo->nvars*sizeof(ulong));
    stride = (ulong *) flint_malloc(ctx->minfo->nvars*sizeof(ulong));
    for (i = 0; i < ctx->minfo->nvars; i++)
    {
        perm[i] = i + 1 < ctx->minfo->nvars ? i + 1 : 0;
        shift[i] = 0;
        stride[i] = 1;
    }

    new_bits = FLINT_MAX(A->bits, B->bits);

    fmpz_mpoly_ctx_init(uctx, ctx->minfo->nvars - 1, ORD_LEX);
    fmpz_mpolyu_init(Au, new_bits, uctx);
    fmpz_mpolyu_init(Bu, new_bits, uctx);
    fmpz_mpolyu_init(Gu, new_bits, uctx);
    fmpz_mpolyu_init(Abaru, new_bits, uctx);
    fmpz_mpolyu_init(Bbaru, new_bits, uctx);

    fmpz_init(cA);
    fmpz_init(cB);
    fmpz_init(cG);

/*
timeit_start(time);
*/
    fmpz_mpoly_to_mpolyu_perm_deflate(Au, A, perm, shift, stride, uctx, ctx);
    fmpz_mpoly_to_mpolyu_perm_deflate(Bu, B, perm, shift, stride, uctx, ctx);
/*
timeit_stop(time);
flint_printf(", (* convert in time *) %wd\n", time->wall);
*/

/*
timeit_start(time);
*/
    fmpz_mpolyu_content(cA, Au, uctx);
    fmpz_mpolyu_content(cB, Bu, uctx);
    fmpz_gcd(cG, cA, cB);
    fmpz_mpolyu_scalar_divexact(Au, cA, uctx);
    fmpz_mpolyu_scalar_divexact(Bu, cB, uctx);
    success = fmpz_mpolyu_gcd_brown_threaded(Gu, Abaru, Bbaru, Au, Bu, uctx, 10);
/*
timeit_stop(time);
flint_printf(", (*un gcd time *) %wd\n", time->wall);
*/
    if (success)
    {
        fmpz_mpoly_from_mpolyu_perm_inflate(G, new_bits, Gu, perm, shift, stride, uctx, ctx);
        if (fmpz_sgn(G->coeffs + 0) < 0)
            fmpz_neg(cG, cG);
        fmpz_mpoly_scalar_mul_fmpz(G, G, cG, ctx);
    }


    fmpz_clear(cA);
    fmpz_clear(cB);
    fmpz_clear(cG);

    fmpz_mpolyu_clear(Au, uctx);
    fmpz_mpolyu_clear(Bu, uctx);
    fmpz_mpolyu_clear(Gu, uctx);
    fmpz_mpolyu_clear(Abaru, uctx);
    fmpz_mpolyu_clear(Bbaru, uctx);
    fmpz_mpoly_ctx_clear(uctx);

    flint_free(perm);
    flint_free(shift);
    flint_free(stride);

    return success;
}
