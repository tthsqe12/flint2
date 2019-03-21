/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "nmod_mpoly.h"
#include "thread_pool.h"
#include "profiler.h"
int usleep(ulong usec);

pthread_mutex_t iomutex;


static __inline__
slong timeit_elapsed_wall(timeit_t t)
{
    struct timeval tv;
    gettimeofday(&tv, 0);
    return t->wall + tv.tv_sec * 1000 + tv.tv_usec / 1000;
}

/* instructions do A = B + I*(C - B) mod M */
typedef struct
{
    slong a_idx; /* index of A */
    slong b_idx; /* index of B */
    slong c_idx; /* index of C */
    nmod_poly_t idem;     /* I */
    nmod_poly_t modulus;  /* M */
} _nmod_poly_crt_prog_instr;

typedef struct
{
    _nmod_poly_crt_prog_instr * prog; /* straight line program */
    slong length; /* length of prog */
    slong alloc;  /* alloc of prog */
    slong localsize; /* length of outputs required in nmod_poly_crt_run */
    slong temp1loc; /* index of temporary used in run */
    slong temp2loc; /* index of another tempory used in run */
    int good;   /* the moduli are good for CRT, essentially relatively prime */
} nmod_poly_crt_struct;

typedef nmod_poly_crt_struct nmod_poly_crt_t[1];

/* general crt for nmod_poly_t - compile once, run many times ****************/

FLINT_DLL void nmod_poly_crt_init(nmod_poly_crt_t P);

FLINT_DLL int nmod_poly_crt_compile(nmod_poly_crt_t P,
                                 nmod_poly_struct * const * moduli, slong len);

NMOD_POLY_INLINE
slong nmod_poly_crt_local_size(nmod_poly_crt_t P)
{
   return P->localsize;
}

FLINT_DLL void nmod_poly_crt_run(const nmod_poly_crt_t P,
                        nmod_poly_t output, nmod_poly_struct * const * inputs);


FLINT_DLL void nmod_poly_crt_clear(nmod_poly_crt_t P);

FLINT_DLL void _nmod_poly_crt_run(const nmod_poly_crt_t P,
                                            nmod_poly_struct * const * outputs,
                                            nmod_poly_struct * const * inputs);

void nmod_poly_crt_init(nmod_poly_crt_t P)
{
    P->prog = NULL;
    P->alloc = 0;
    P->length = 0;
    P->localsize = 1;
    P->temp1loc = 0;
    P->temp2loc = 0;
    P->good = 0;
}

static void nmod_poly_crt_fit_length(nmod_poly_crt_t P, slong k)
{
    k = FLINT_MAX(WORD(1), k);

    if (P->alloc == 0)
    {
        FLINT_ASSERT(P->prog == NULL);
        P->prog = (_nmod_poly_crt_prog_instr *) flint_malloc(k
                                           *sizeof(_nmod_poly_crt_prog_instr));
        P->alloc = k;
    }
    else if (k > P->alloc)
    {
        FLINT_ASSERT(P->prog != NULL);
        P->prog = (_nmod_poly_crt_prog_instr *) flint_realloc(P->prog, k
                                           *sizeof(_nmod_poly_crt_prog_instr));
        P->alloc = k;
    }
}

static void nmod_poly_crt_set_length(nmod_poly_crt_t P, slong k)
{
    slong i;

    FLINT_ASSERT(k <= P->length);

    for (i = k; i < P->length; i++)
    {
        nmod_poly_clear(P->prog[i].modulus);
        nmod_poly_clear(P->prog[i].idem);
    }
    P->length = k;
}

void nmod_poly_crt_clear(nmod_poly_crt_t P)
{
    nmod_poly_crt_set_length(P, 0);

    if (P->alloc > 0)
    {
        flint_free(P->prog);
    }
}

/*
    combine all moduli in [start, stop)
    return index of instruction that computes the result
*/
static slong _push_prog(nmod_poly_crt_t P,
                           nmod_poly_struct * const * moduli, slong * perm,
                                        slong ret_idx, slong start, slong stop)
{
    slong i, mid;
    slong b_idx, c_idx;
    slong lefttot, righttot;
    slong leftret, rightret;
    nmod_poly_struct * leftmodulus, * rightmodulus;

    /* we should have at least 2 moduli */
    FLINT_ASSERT(start + 1 < stop);

    mid = start + (stop - start)/2;

    FLINT_ASSERT(start < mid);
    FLINT_ASSERT(mid < stop);

    lefttot = 0;
    for (i = start; i < mid; i++)
    {
        lefttot += nmod_poly_degree(moduli[perm[i]]);
    }

    righttot = 0;
    for (i = mid; i < stop; i++)
    {
        righttot += nmod_poly_degree(moduli[perm[i]]);
    }

    /* try to balance the total degree on left and right */
    while (lefttot < righttot
            && mid + 1 < stop
            && nmod_poly_degree(moduli[perm[mid]]) < righttot - lefttot)
    {
        lefttot += nmod_poly_degree(moduli[perm[mid]]);
        righttot -= nmod_poly_degree(moduli[perm[mid]]);
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
        leftmodulus = (nmod_poly_struct *) moduli[perm[start]];
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
        rightmodulus = (nmod_poly_struct *) moduli[perm[mid]];
    }

    /* check if nmod_poly_invmod is going to throw */
    if (nmod_poly_degree(leftmodulus) < 1 || nmod_poly_degree(rightmodulus) < 1)
    {
        P->good = 0;
        return -1;
    }

    /* compile [start, end) */
    i = P->length;
    nmod_poly_crt_fit_length(P, i + 1);
    nmod_poly_init(P->prog[i].modulus, rightmodulus->mod.n);
    nmod_poly_init(P->prog[i].idem, rightmodulus->mod.n);
    P->good = P->good && nmod_poly_invmod(P->prog[i].modulus, leftmodulus, rightmodulus);
    nmod_poly_mul(P->prog[i].idem, leftmodulus, P->prog[i].modulus);
    nmod_poly_mul(P->prog[i].modulus, leftmodulus, rightmodulus);
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
int nmod_poly_crt_compile(nmod_poly_crt_t P,
                                  nmod_poly_struct * const * moduli, slong len)
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
        for (j = i; j > 0 && nmod_poly_degree(moduli[perm[j-1]])
                           > nmod_poly_degree(moduli[perm[j-0]]); j--)
        {
            slong temp = perm[j-1];
            perm[j-1] = perm[j-0];
            perm[j-0] = temp;
        }
    }

    nmod_poly_crt_fit_length(P, FLINT_MAX(WORD(1), len - 1));
    nmod_poly_crt_set_length(P, 0);
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
        nmod_poly_init(P->prog[i].modulus, moduli[0]->mod.n);
        nmod_poly_init(P->prog[i].idem, moduli[0]->mod.n);
        nmod_poly_set(P->prog[i].modulus, moduli[0]);
        P->prog[i].a_idx = 0;
        P->prog[i].b_idx = -WORD(1);
        P->prog[i].c_idx = -WORD(1);
        P->length = i + 1;

        P->good = !nmod_poly_is_zero(moduli[0]);
    }

    if (!P->good)
    {
        nmod_poly_crt_set_length(P, 0);
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
void _nmod_poly_crt_run(const nmod_poly_crt_t P,
                                            nmod_poly_struct * const * outputs,
                                            nmod_poly_struct * const * inputs)
{
    slong i;
    slong a, b, c;
    nmod_poly_struct * A, * B, * C, * t1, * t2;

    t1 = outputs[P->temp1loc];
    t2 = outputs[P->temp2loc];

    for (i = 0; i < P->length; i++)
    {
        a = P->prog[i].a_idx;
        b = P->prog[i].b_idx;
        c = P->prog[i].c_idx;
        FLINT_ASSERT(a >= 0);
        A = outputs[a];
        B = b < 0 ? inputs[-b-1] : outputs[b];
        C = c < 0 ? inputs[-c-1] : outputs[c];

        /* A = B + I*(C - B) mod M */
        nmod_poly_sub(t1, B, C);
        nmod_poly_mul(t2, P->prog[i].idem, t1);
        nmod_poly_sub(t1, B, t2);

        if (nmod_poly_degree(t1) < nmod_poly_degree(P->prog[i].modulus))
        {
            nmod_poly_swap(A, t1);
        }
        else
        {
            nmod_poly_rem(A, t1, P->prog[i].modulus);
        }

        /* last calculation should write answer to outputs[0] */
        if (i + 1 >= P->length)
        {
            FLINT_ASSERT(A == outputs[0]);
        }
    }
}

void nmod_poly_crt_run(const nmod_poly_crt_t P, nmod_poly_t output,
                                             nmod_poly_struct * const * inputs)
{
    slong i;
    slong a, b, c;
    nmod_poly_struct * A, * B, * C, * t1, * t2, * outputs;
    TMP_INIT;

    TMP_START;
    outputs = (nmod_poly_struct *) TMP_ALLOC(P->localsize
                                                    *sizeof(nmod_poly_struct));
    for (i = 0; i < P->localsize; i++)
    {
        nmod_poly_init(outputs + i, inputs[0]->mod.n);
    }

    nmod_poly_swap(output, outputs + 0);

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
        nmod_poly_sub(t1, B, C);
        nmod_poly_mul(t2, P->prog[i].idem, t1);
        nmod_poly_sub(t1, B, t2);
        nmod_poly_rem(A, t1, P->prog[i].modulus);
    }

    nmod_poly_swap(output, outputs + 0);

    for (i = 0; i < P->localsize; i++)
    {
        nmod_poly_clear(outputs + i);
    }

    TMP_END;
}










/* return maxdegree crt(B[0], ...., B[l-1]) wrt to P */
slong nmod_mpolyn_crt1(
    const nmod_poly_crt_t P,
    nmod_mpolyn_struct * const * B,
    slong count,
    const nmod_mpoly_ctx_t ctx)
{
    int cmp;
    slong N = mpoly_words_per_exp_sp(B[0]->bits, ctx->minfo);
    slong lastdeg;
    nmod_poly_struct ** input;
    slong * start;
    slong j, k;
    nmod_poly_t zero, tcoeff;
    ulong * texp;
    TMP_INIT;

    TMP_START;

    nmod_poly_init(zero, ctx->ffinfo->mod.n);
    nmod_poly_init(tcoeff, ctx->ffinfo->mod.n);
    texp = (ulong *) TMP_ALLOC(N*sizeof(ulong));

    input = (nmod_poly_struct **) TMP_ALLOC(count * sizeof(nmod_poly_struct *));
    start = (slong *) TMP_ALLOC(count * sizeof(slong));

    /* start[k] is the next available term in B[k] */
    for (k = 0; k < count; k++)
    {
        start[k] = 0;    
    }

    lastdeg = -WORD(1);
    while (1)
    {
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
        mpoly_monomial_set(texp, B[k]->exps + N*start[k], N);
        start[k]++;

        for (k++; k < count; k++)
        {
            input[k] = zero;
            if (start[k] >= B[k]->length)
            {
                continue;
            }

            cmp = mpoly_monomial_cmp_nomask(B[k]->exps + N*start[k], texp, N);
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

        nmod_poly_crt_run(P, tcoeff, input);
        lastdeg = FLINT_MAX(lastdeg, nmod_poly_degree(tcoeff));
    }

    TMP_END;

    nmod_poly_clear(zero);
    nmod_poly_clear(tcoeff);
    return lastdeg;
}

/*
    return the last degree of the result of crt'ing the coeff of X^exp
*/
slong nmod_mpolyun_crt_exp1(
    const nmod_poly_crt_t P,
    ulong exp,
    nmod_mpolyun_struct * const * B,
    slong count,
    const nmod_mpoly_ctx_t ctx)
{
    slong j, k;
    slong lastdegree;
    nmod_mpolyn_struct ** C;
    nmod_mpolyn_t zero;
    TMP_INIT;

    nmod_mpolyn_init(zero, B[0]->bits, ctx);

    TMP_START;
    C = (nmod_mpolyn_struct **) TMP_ALLOC(count * sizeof(nmod_mpolyn_struct *));
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

    lastdegree = nmod_mpolyn_crt1(P, C, count, ctx);

    TMP_END;
    nmod_mpolyn_clear(zero, ctx);
    return lastdegree;
}







/* A = crt(B[0], ...., B[l-1]) wrt to P */
slong nmod_mpolyn_crt(
    const nmod_poly_crt_t P,
    nmod_mpolyn_t A,
    nmod_mpolyn_struct * const * B,
    slong count,
    const nmod_mpoly_ctx_t ctx)
{
    int cmp;
    slong N = mpoly_words_per_exp_sp(A->bits, ctx->minfo);
    slong lastdegree;
    nmod_poly_struct ** input;
    slong * start;
    slong Ai;
    slong j, k;
    nmod_poly_t zero;
    TMP_INIT;

    TMP_START;

    nmod_poly_init(zero, ctx->ffinfo->mod.n);

    input = (nmod_poly_struct **) TMP_ALLOC(count * sizeof(nmod_poly_struct *));
    start = (slong *) TMP_ALLOC(count * sizeof(slong));

    /* start[k] is the next available term in B[k] */
    for (k = 0; k < count; k++)
    {
        start[k] = 0;    
    }

    Ai = 0;
    lastdegree = -WORD(1);
    while (1)
    {
        nmod_mpolyn_fit_length(A, Ai + 1, ctx);

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

        nmod_poly_crt_run(P, A->coeffs + Ai, input);
        lastdegree = FLINT_MAX(lastdegree, nmod_poly_degree(A->coeffs + Ai));
        Ai += !nmod_poly_is_zero(A->coeffs + Ai);
    }
    A->length = Ai;

    TMP_END;

    nmod_poly_clear(zero);
    return lastdegree;
}

/*
    Append to A the result of crt'ing the coeff of X^exp
*/
slong nmod_mpolyun_crt_exp(
    const nmod_poly_crt_t P,
    nmod_mpolyun_t A,
    ulong exp,
    nmod_mpolyun_struct * const * B,
    slong count,
    const nmod_mpoly_ctx_t ctx)
{
    slong j, k;
    slong Ai;
    slong lastdegree;
    nmod_mpolyn_struct ** C;
    nmod_mpolyn_t zero;
    TMP_INIT;

    nmod_mpolyn_init(zero, A->bits, ctx);

    TMP_START;
    C = (nmod_mpolyn_struct **) TMP_ALLOC(count * sizeof(nmod_mpolyn_struct *));
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
    nmod_mpolyun_fit_length(A, Ai + 1, ctx);
    A->exps[Ai] = exp;
    lastdegree = nmod_mpolyn_crt(P, A->coeffs + Ai, C, count, ctx);
    A->length += (A->coeffs + Ai)->length != 0;

    TMP_END;
    nmod_mpolyn_clear(zero, ctx);
    return lastdegree;
}





int nmod_mpolyn_is_canonical(const nmod_mpolyn_t A, const nmod_mpoly_ctx_t ctx)
{
    slong i;

    if (!mpoly_monomials_valid_test(A->exps, A->length, A->bits, ctx->minfo))
    {
        return 0;
    }

    if (mpoly_monomials_overflow_test(A->exps, A->length, A->bits, ctx->minfo))
        return 0;

    if (!mpoly_monomials_inorder_test(A->exps, A->length, A->bits, ctx->minfo))
        return 0;

    for (i = 0; i < A->length; i++)
    {
        slong l = (A->coeffs + i)->length;

        if (l == 0)
        {
            printf("\n mpoln length zero!!!!!!!\n");
            return 0;
        }

        if ((A->coeffs + i)->coeffs[l - 1] == 0)
        {
            printf("\n mpolyn top coeff zero!!!!!!!\n");
            return 0;
        }
    }

    return 1;
}

int nmod_mpolyun_is_canonical(const nmod_mpolyun_t A, const nmod_mpoly_ctx_t ctx)
{
    slong i;

    if (A->length > A->alloc)
    {
        return 0;
    }

    for (i = 0; i < A->length; i++)
    {
        if (!nmod_mpolyn_is_canonical(A->coeffs + i, ctx))
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

/*
    F = F + modulus*(A - F(v = alpha))
    no assumptions about matching monomials
    F is in Fp[X][v]
    A is in Fp[X]
    it is expected that modulus(alpha) == 1
*/
int nmod_mpolyun_addinterp_bivar(slong * lastdeg_,
                   nmod_mpolyun_t F, nmod_mpolyun_t T, const nmod_poly_t A,
      const nmod_poly_t modulus,  mp_limb_t alpha,  const nmod_mpoly_ctx_t ctx)
{
    int changed = 0;
    slong lastdeg = -WORD(1);
    slong N = mpoly_words_per_exp_sp(F->bits, ctx->minfo);
    mp_limb_t v;
    slong Fi, Toff, Aexp;
    mp_limb_t * Acoeff = A->coeffs;
    slong Flen = F->length;
    nmod_mpolyn_struct * Fcoeff = F->coeffs;
    ulong * Fexp = F->exps;
    nmod_mpolyn_struct * Tcoeff;
    ulong * Texp;
    nmod_poly_t tp;
    
    Fi = 0;

/*
flint_printf("\nnmod_mpolyun_addinterp_bivar: (alpha = %wu)\n", alpha);
flint_printf("A: "); nmod_poly_print_pretty(A, "X"); printf("\n");
flint_printf("F: "); nmod_mpolyun_print_pretty(F, NULL, ctx); printf("\n");
*/

    Aexp = nmod_poly_degree(A);

    nmod_poly_init(tp, ctx->ffinfo->mod.n);

    nmod_mpolyun_fit_length(T, Flen + Aexp + 1, ctx);
    Tcoeff = T->coeffs;
    Texp = T->exps;
    Toff = 0;

    while (Fi < Flen || Aexp >= 0)
    {
        FLINT_ASSERT(Toff < T->alloc);
/*
flint_printf("Fi: %wd\n",Fi);
flint_printf("Aexp: %wd\n",Aexp);
*/
        if (Fi < Flen)
        {
            FLINT_ASSERT(!nmod_poly_is_zero((Fcoeff + Fi)->coeffs + 0));
            FLINT_ASSERT(nmod_poly_degree((Fcoeff + Fi)->coeffs + 0) < nmod_poly_degree(modulus));
        }

        if (Aexp >= 0)
        {
            FLINT_ASSERT(Acoeff[Aexp] != UWORD(0));
        }

        nmod_mpolyn_fit_length(Tcoeff + Toff, 1, ctx);

        if (Fi < Flen && Aexp >= 0 && Fexp[Fi] == Aexp)
        {
            /* F term ok, A term ok */
            v = nmod_poly_evaluate_nmod((Fcoeff + Fi)->coeffs + 0, alpha);
            v = nmod_sub(Acoeff[Aexp], v, ctx->ffinfo->mod);
            if (v != UWORD(0))
            {
                changed = 1;
                nmod_poly_scalar_mul_nmod(tp, modulus, v);
                nmod_poly_add((Tcoeff + Toff)->coeffs + 0, (Fcoeff + Fi)->coeffs + 0, tp);
            }
            else
            {
                nmod_poly_set((Tcoeff + Toff)->coeffs + 0, (Fcoeff + Fi)->coeffs + 0);                
            }
            Texp[Toff] = Aexp;
            Fi++;
            do {
                Aexp--;
            } while (Aexp >= 0 && Acoeff[Aexp] == UWORD(0));
        }
        else if (Fi < Flen && (Aexp < 0 || Fexp[Fi] > Aexp))
        {
            /* F term ok, A term missing */
            v = nmod_poly_evaluate_nmod((Fcoeff + Fi)->coeffs + 0, alpha);
            if (v != UWORD(0))
            {
                changed = 1;
                nmod_poly_scalar_mul_nmod(tp, modulus, v);
                nmod_poly_sub((Tcoeff + Toff)->coeffs + 0, (Fcoeff + Fi)->coeffs + 0, tp);
            }
            else
            {
                nmod_poly_set((Tcoeff + Toff)->coeffs + 0, (Fcoeff + Fi)->coeffs + 0);                
            }
            Texp[Toff] = Fexp[Fi];
            Fi++;
        }
        else if (Aexp >= 0 && (Fi >= Flen || Fexp[Fi] < Aexp))
        {
            /* F term missing, A term ok */
            changed = 1;
            nmod_poly_scalar_mul_nmod((Tcoeff + Toff)->coeffs + 0, modulus, Acoeff[Aexp]);
            Texp[Toff] = Aexp;
            do {
                Aexp--;
            } while (Aexp >= 0 && Acoeff[Aexp] == UWORD(0));
        }
        else
        {
            FLINT_ASSERT(0);
        }

        lastdeg = FLINT_MAX(lastdeg, nmod_poly_degree((Tcoeff + Toff)->coeffs + 0));
        mpoly_monomial_zero((Tcoeff + Toff)->exps + N*0, N);
        FLINT_ASSERT(!nmod_poly_is_zero((Tcoeff + Toff)->coeffs + 0));
        (Tcoeff + Toff)->length = 1;
        Toff++;
    }
    T->length = Toff;

    nmod_poly_clear(tp);

    if (changed)
    {
        nmod_mpolyun_swap(T, F);
    }

    *lastdeg_ = lastdeg;
    return changed;
}


/*
    set vp = P(alpha), vm = P(-alpha) given powers of alpha
*/
void _nmod_poly_eval2_pow(mp_limb_t * vp, mp_limb_t * vm, nmod_poly_t P, 
                                  nmod_poly_t alphapow, const nmodf_ctx_t fctx)
{
    mp_limb_t * Pcoeffs = P->coeffs;
    slong Plen = P->length;
    mp_limb_t * alpha_powers = alphapow->coeffs;
    mp_limb_t p1, p0, a0, a1, a2, q1, q0, b0, b1, b2;
    slong k;

    a0 = a1 = a2 = UWORD(0);
    b0 = b1 = b2 = UWORD(0);

    if (Plen > alphapow->length)
    {
        slong oldlength = alphapow->length;
        FLINT_ASSERT(2 <= oldlength);
        nmod_poly_fit_length(alphapow, Plen);
        for (k = oldlength; k < Plen; k++)
        {
            alphapow->coeffs[k] = nmod_mul(alphapow->coeffs[k - 1],
                                               alphapow->coeffs[1], fctx->mod);
        }
        alphapow->length = Plen;
    }

    for (k = 0; k + 2 <= Plen; k += 2)
    {
        umul_ppmm(p1, p0, Pcoeffs[k + 0], alpha_powers[k + 0]);
        umul_ppmm(q1, q0, Pcoeffs[k + 1], alpha_powers[k + 1]);
        add_sssaaaaaa(a2, a1, a0, a2, a1, a0, WORD(0), p1, p0);
        add_sssaaaaaa(b2, b1, b0, b2, b1, b0, WORD(0), q1, q0);
    }

    if (k < Plen)
    {
        umul_ppmm(p1, p0, Pcoeffs[k + 0], alpha_powers[k + 0]);
        add_sssaaaaaa(a2, a1, a0, a2, a1, a0, WORD(0), p1, p0);
        k++;
    }

    FLINT_ASSERT(k == Plen);

    NMOD_RED3(p0, a2, a1, a0, fctx->mod);
    NMOD_RED3(q0, b2, b1, b0, fctx->mod);

    vp[0] = nmod_add(p0, q0, fctx->mod);
    vm[0] = nmod_sub(p0, q0, fctx->mod);
}

/*
    update F from its value A at v = alpha and its value B at v = -alpha
    no assumptions about matching monomials
    F is in R[X][v]
    A is in R[X]
    B is in R[X]
    it is expected that modulus(alpha) == modulus(-alpha) == 1/(2*alpha)
*/
int nmod_mpolyun_addinterp2_bivar(slong * lastdeg_,
                   nmod_mpolyun_t F, nmod_mpolyun_t T, const nmod_poly_t A, const nmod_poly_t B,
      const nmod_poly_t modulus, nmod_poly_t alphapow, const nmod_mpoly_ctx_t ctx)
{
    int changed = 0, Finc;
    mp_limb_t alpha = nmod_poly_get_coeff_ui(alphapow, 1);
    slong lastdeg = -WORD(1);
    slong N = mpoly_words_per_exp_sp(F->bits, ctx->minfo);
    mp_limb_t u, v, FvalueA, FvalueB;
    slong Fi, Toff, Aexp, Bexp, e;
    mp_limb_t * Acoeff = A->coeffs;
    mp_limb_t * Bcoeff = B->coeffs;
    slong Flen = F->length;
    nmod_mpolyn_struct * Fcoeff = F->coeffs;
    ulong * Fexp = F->exps;
    nmod_mpolyn_struct * Tcoeff;
    ulong * Texp;
    nmod_poly_t tp;

    Fi = 0;
    Aexp = nmod_poly_degree(A);
    Bexp = nmod_poly_degree(B);

    nmod_poly_init(tp, ctx->ffinfo->mod.n);

    nmod_mpolyun_fit_length(T, Flen + FLINT_MAX(Aexp, Bexp) + 1, ctx);
    Tcoeff = T->coeffs;
    Texp = T->exps;
    Toff = 0;

#if WANT_ASSERT
    u = nmod_poly_evaluate_nmod(modulus, alpha);
    u = nmod_mul(u, alpha, ctx->ffinfo->mod);
    u = nmod_mul(u, 2, ctx->ffinfo->mod);
    FLINT_ASSERT(u == 1);
    u = nmod_poly_evaluate_nmod(modulus, ctx->ffinfo->mod.n - alpha);
    u = nmod_mul(u, alpha, ctx->ffinfo->mod);
    u = nmod_mul(u, 2, ctx->ffinfo->mod);
    FLINT_ASSERT(u == 1);
#endif
/*
pthread_mutex_lock(&iomutex);
flint_printf("addinterp2_bivar called  alpha = %wu\n", alpha);
flint_printf("A: "); nmod_poly_print_pretty(A, "X"); printf("\n");
flint_printf("B: "); nmod_poly_print_pretty(B, "X"); printf("\n");
flint_printf("F: "); nmod_mpolyun_print_pretty(F, NULL, ctx); printf("\n");
*/

    while (Fi < Flen || Aexp >= 0 || Bexp >= 0)
    {
        FLINT_ASSERT(Toff < T->alloc);

        e = -WORD(1);
        if (Fi < Flen)
        {
            e = Fexp[Fi];
            FLINT_ASSERT(!nmod_poly_is_zero((Fcoeff + Fi)->coeffs + 0));
            FLINT_ASSERT(nmod_poly_degree((Fcoeff + Fi)->coeffs + 0) < nmod_poly_degree(modulus));
        }
        if (Aexp >= 0)
        {
            e = FLINT_MAX(e, Aexp);
            FLINT_ASSERT(Acoeff[Aexp] != UWORD(0));
        }
        if (Bexp >= 0)
        {
            e = FLINT_MAX(e, Bexp);
            FLINT_ASSERT(Bcoeff[Bexp] != UWORD(0));
        }

        FLINT_ASSERT(e >= 0);
        nmod_mpolyn_fit_length(Tcoeff + Toff, 1, ctx);
        Texp[Toff] = e;

        FvalueA = FvalueB = 0;
        Finc = 0;
        if (Fi < Flen && e == Fexp[Fi])
        {
            Finc = 1;
            _nmod_poly_eval2_pow(&FvalueA, &FvalueB, (Fcoeff + Fi)->coeffs + 0, alphapow, ctx->ffinfo);
        }

        if (e == Aexp)
        {
            FvalueA = nmod_sub(FvalueA, Acoeff[Aexp], ctx->ffinfo->mod);
        }
        if (e == Bexp)
        {
            FvalueB = nmod_sub(FvalueB, Bcoeff[Bexp], ctx->ffinfo->mod);
        }

        u = nmod_sub(FvalueB, FvalueA, ctx->ffinfo->mod);
        v = nmod_mul(ctx->ffinfo->mod.n - alpha, nmod_add(FvalueB, FvalueA, ctx->ffinfo->mod), ctx->ffinfo->mod);

        if (u != 0 || v != 0)
        {
            changed = 1;
            nmod_poly_set_coeff_ui(tp, 0, v);
            nmod_poly_set_coeff_ui(tp, 1, u);
            nmod_poly_mul_classical((Tcoeff + Toff)->coeffs + 0, modulus, tp);
            if (Finc)
            {
                nmod_poly_add((Tcoeff + Toff)->coeffs + 0, (Tcoeff + Toff)->coeffs + 0, (Fcoeff + Fi)->coeffs + 0);
            }
        }
        else
        {
            FLINT_ASSERT(Finc == 1);
            nmod_poly_set((Tcoeff + Toff)->coeffs + 0, (Fcoeff + Fi)->coeffs + 0);
        }

        lastdeg = FLINT_MAX(lastdeg, nmod_poly_degree((Tcoeff + Toff)->coeffs + 0));
        mpoly_monomial_zero((Tcoeff + Toff)->exps + N*0, N);
        FLINT_ASSERT(!nmod_poly_is_zero((Tcoeff + Toff)->coeffs + 0));
        (Tcoeff + Toff)->length = 1;
        Toff++;

        Fi += Finc;
        if (e == Aexp)
        {
            do {
                Aexp--;
            } while (Aexp >= 0 && Acoeff[Aexp] == 0);
        }
        if (e == Bexp)
        {
            do {
                Bexp--;
            } while (Bexp >= 0 && Bcoeff[Bexp] == 0);
        }
    }
    T->length = Toff;

    nmod_poly_clear(tp);

    if (changed)
    {
        nmod_mpolyun_swap(T, F);
    }
/*
flint_printf("returning F: "); nmod_mpolyun_print_pretty(F, NULL, ctx); printf("\n");
pthread_mutex_unlock(&iomutex);
*/

    *lastdeg_ = lastdeg;
    return changed;
}




/*
    E = A(v = alpha), F = A(v = -alpha)
    A is in R[X][v]
    E is in R[X]
    F is in R[X]
*/
void nmod_mpolyun_eval2_last_bivar(nmod_poly_t E, nmod_poly_t F,
                      const nmod_mpolyun_t A, nmod_poly_t alphapow,
                                                    const nmod_mpoly_ctx_t ctx)
{
    mp_limb_t u, v;
    slong Ai, Alen;
    nmod_mpolyn_struct * Acoeff;
    ulong * Aexp;

    Acoeff = A->coeffs;
    Aexp = A->exps;
    Alen = A->length;
    Ai = 0;
    nmod_poly_zero(E);
    nmod_poly_zero(F);
    for (Ai = 0; Ai < Alen; Ai++)
    {
        _nmod_poly_eval2_pow(&u, &v, (Acoeff + Ai)->coeffs + 0, alphapow, ctx->ffinfo);
        nmod_poly_set_coeff_ui(E, Aexp[Ai], u);
        nmod_poly_set_coeff_ui(F, Aexp[Ai], v);
    }
}

/*
    set F from its value A at v = alpha and its value B at v = -alpha
    no assumptions about matching monomials
    F is in R[X][v]
    A is in R[X]
    B is in R[X]
*/
void nmod_mpolyun_startinterp2_bivar(slong * lastdeg_,
                   nmod_mpolyun_t F, const nmod_poly_t A, const nmod_poly_t B,
                                  mp_limb_t alpha,  const nmod_mpoly_ctx_t ctx)
{
    slong lastdeg = -WORD(1);
    slong N = mpoly_words_per_exp_sp(F->bits, ctx->minfo);
    mp_limb_t u, v, d0, d1, Avalue, Bvalue;
    slong Fi, Aexp, Bexp;
    mp_limb_t * Acoeff = A->coeffs;
    mp_limb_t * Bcoeff = B->coeffs;
    nmod_mpolyn_struct * Fcoeff;
    ulong * Fexp;
    slong e;

    Aexp = nmod_poly_degree(A);
    Bexp = nmod_poly_degree(B);

    nmod_mpolyun_fit_length(F, FLINT_MAX(Aexp, Bexp) + 1, ctx);
    Fcoeff = F->coeffs;
    Fexp = F->exps;

    d0 = n_invmod(UWORD(2), ctx->ffinfo->mod.n);
    d1 = n_invmod(nmod_add(alpha, alpha, ctx->ffinfo->mod), ctx->ffinfo->mod.n);

    Fi = 0;
    while (Aexp >= 0 || Bexp >= 0)
    {
        e = Aexp;
        Avalue = 0;
        Bvalue = 0;
        if (Aexp == Bexp)
        {
            Avalue = Acoeff[Aexp];
            Bvalue = Bcoeff[Bexp];
        }
        else if (Aexp > Bexp)
        {
            Avalue = Acoeff[Aexp];
        }
        else
        {
            FLINT_ASSERT(Bexp > Aexp);
            e = Bexp;
            Bvalue = Bcoeff[Bexp];
        }
        FLINT_ASSERT(Avalue != 0 || Bvalue != 0);
        u = nmod_add(Avalue, Bvalue, ctx->ffinfo->mod);
        v = nmod_sub(Avalue, Bvalue, ctx->ffinfo->mod);
        u = nmod_mul(u, d0, ctx->ffinfo->mod);
        v = nmod_mul(v, d1, ctx->ffinfo->mod);

        FLINT_ASSERT(Fi < F->alloc);
        nmod_mpolyn_fit_length(Fcoeff + Fi, 1, ctx);
        mpoly_monomial_zero((Fcoeff + Fi)->exps + N*0, N);
        nmod_poly_zero((Fcoeff + Fi)->coeffs + 0);
        nmod_poly_set_coeff_ui((Fcoeff + Fi)->coeffs + 0, 0, u);
        nmod_poly_set_coeff_ui((Fcoeff + Fi)->coeffs + 0, 1, v);
        lastdeg = FLINT_MAX(lastdeg, nmod_poly_degree((Fcoeff + Fi)->coeffs + 0));
        Fexp[Fi] = e;
        (Fcoeff + Fi)->length = 1;
        Fi++;

        if (e == Aexp)
        {
            do {
                Aexp--;
            } while (Aexp >= 0 && Acoeff[Aexp] == 0);
        }
        if (e == Bexp)
        {
            do {
                Bexp--;
            } while (Bexp >= 0 && Bcoeff[Bexp] == 0);
        }
    }
    F->length = Fi;
/*
pthread_mutex_lock(&iomutex);
flint_printf("startinterp2_bivar returning  alpha = %wu\n", alpha);
flint_printf("A: "); nmod_poly_print_pretty(A, "X"); printf("\n");
flint_printf("B: "); nmod_poly_print_pretty(B, "X"); printf("\n");
flint_printf("F: "); nmod_mpolyun_print_pretty(F, NULL, ctx); printf("\n");
pthread_mutex_unlock(&iomutex);
*/
    *lastdeg_ = lastdeg;
    return;
}




/*
    set T from
        value F modulo modulus
        value A at x_var = alpha
        value B at x_var = -alpha
    it is expected that modulus(alpha) == modulus(-alpha) == 1/(2*alpha)
    no assumptions about matching monomials
    F is in Fp[x_0, ..., x_(var-1), x_(var-1)][x_var]
    A is in Fp[x_0, ..., x_(var-2)][x_(var-1)]
    B is in Fp[x_0, ..., x_(var-2)][x_(var-1)]
*/
int nmod_mpolyn_addinterp2_n(slong * lastdeg_,
         nmod_mpolyn_t T, nmod_mpolyn_t F, nmod_mpolyn_t A, nmod_mpolyn_t B,
                   slong var, nmod_poly_t modulus, nmod_poly_t alphapow,
                                                    const nmod_mpoly_ctx_t ctx)
{
    int changed = 0;
    slong N = mpoly_words_per_exp_sp(A->bits, ctx->minfo);
    slong lastdeg = -WORD(1);
    slong offset, shift;
    nmod_poly_t tp;
    nmod_poly_t zero;

    nmod_poly_struct * Tcoeff;
    ulong * Texp;
    slong Ti;

    nmod_poly_struct * Fcoeff = F->coeffs;
    slong Flen = F->length;
    ulong * Fexp = F->exps;
    slong Fi;

    nmod_poly_struct * Acoeff = A->coeffs;
    slong Alen = A->length;
    ulong * Aexp = A->exps;
    slong Ai, ai;

    nmod_poly_struct * Bcoeff = B->coeffs;
    slong Blen = B->length;
    ulong * Bexp = B->exps;
    slong Bi, bi;

    nmod_poly_struct * Fvalue;
    mp_limb_t u, v, Avalue, Bvalue, FvalueA, FvalueB;
    int texp_set, cmp;
    mp_limb_t alpha = nmod_poly_get_coeff_ui(alphapow, 1);

#if WANT_ASSERT
    u = nmod_poly_evaluate_nmod(modulus, alpha);
    u = nmod_mul(u, alpha, ctx->ffinfo->mod);
    u = nmod_mul(u, 2, ctx->ffinfo->mod);
    FLINT_ASSERT(u == 1);
    u = nmod_poly_evaluate_nmod(modulus, ctx->ffinfo->mod.n - alpha);
    u = nmod_mul(u, alpha, ctx->ffinfo->mod);
    u = nmod_mul(u, 2, ctx->ffinfo->mod);
    FLINT_ASSERT(u == 1);
#endif

    FLINT_ASSERT(nmod_mpolyn_is_canonical(A, ctx));
    FLINT_ASSERT(nmod_mpolyn_is_canonical(B, ctx));
    FLINT_ASSERT(nmod_mpolyn_is_canonical(F, ctx));

    nmod_poly_init(tp, ctx->ffinfo->mod.n);
    nmod_poly_init(zero, ctx->ffinfo->mod.n);

    FLINT_ASSERT(var > 0);
    FLINT_ASSERT(T->bits == A->bits);
    FLINT_ASSERT(F->bits == A->bits);

    mpoly_gen_offset_shift_sp(&offset, &shift, var - 1, A->bits, ctx->minfo);

    Flen = F->length;

    nmod_mpolyn_fit_length(T, FLINT_MAX(Flen, Alen), ctx);
    Tcoeff = T->coeffs;
    Texp = T->exps;

    Ti = Fi = Ai = Bi = 0;
    ai = (Ai >= Alen) ? 0 : nmod_poly_degree(A->coeffs + Ai);
    bi = (Bi >= Blen) ? 0 : nmod_poly_degree(B->coeffs + Bi);

    while (Fi < Flen || Ai < Alen || Bi < Blen)
    {
        if (Ti >= T->alloc)
        {
            slong extra = Flen - Fi;
            extra = FLINT_MAX(extra, Alen - Ai);
            extra = FLINT_MAX(extra, Blen - Bi);
            nmod_mpolyn_fit_length(T, Ti + extra, ctx);
            Tcoeff = T->coeffs;
            Texp = T->exps;
        }

        FLINT_ASSERT(Fi >= Flen || (Fcoeff + Fi)->length != 0);
        FLINT_ASSERT(Ai >= Alen || (Acoeff + Ai)->coeffs[ai] != 0);
        FLINT_ASSERT(Bi >= Blen || (Bcoeff + Bi)->coeffs[bi] != 0);

        Fvalue = zero;
        texp_set = 0;
        if (Fi < Flen)
        {
            Fvalue = Fcoeff + N*Fi;
            texp_set = 1;
            mpoly_monomial_set(Texp + N*Ti, Fexp + N*Fi, N);
        }

        Avalue = 0;
        if (Ai < Alen)
        {
            cmp = (!texp_set) ? -1
                     : mpoly_monomial_cmp_nomask_extra(Texp + N*Ti,
                                      Aexp + N*Ai, N, offset, ai << shift);

            if (cmp <= 0)
            {
                Avalue = (Acoeff + Ai)->coeffs[ai];
            }
            if (cmp < 0)
            {
                Fvalue = zero;
                texp_set = 1;
                mpoly_monomial_set_extra(Texp + N*Ti,
                                         Aexp + N*Ai, N, offset, ai << shift);
            }
        }

        Bvalue = 0;
        if (Bi < Blen)
        {
            cmp = (!texp_set) ? -1 : mpoly_monomial_cmp_nomask_extra(
                             Texp + N*Ti, Bexp + N*Bi, N, offset, bi << shift);

            if (cmp <= 0)
            {
                Bvalue = (Bcoeff + Bi)->coeffs[bi];
            }
            if (cmp < 0)
            {
                Fvalue = zero;
                Avalue = 0;
                texp_set = 1;
                mpoly_monomial_set_extra(Texp + N*Ti,
                                         Bexp + N*Bi, N, offset, bi << shift);
            }
        }

        FLINT_ASSERT(texp_set);

        _nmod_poly_eval2_pow(&FvalueA, &FvalueB, Fvalue, alphapow, ctx->ffinfo);
        FvalueA = nmod_sub(FvalueA, Avalue, ctx->ffinfo->mod);
        FvalueB = nmod_sub(FvalueB, Bvalue, ctx->ffinfo->mod);
        u = nmod_sub(FvalueB, FvalueA, ctx->ffinfo->mod);
        v = nmod_mul(ctx->ffinfo->mod.n - alpha, nmod_add(FvalueB, FvalueA, ctx->ffinfo->mod), ctx->ffinfo->mod);
        if (u != 0 || v != 0)
        {
            changed = 1;
            nmod_poly_set_coeff_ui(tp, 0, v);
            nmod_poly_set_coeff_ui(tp, 1, u);
            nmod_poly_mul_classical(Tcoeff + Ti, modulus, tp);
            nmod_poly_add(Tcoeff + Ti, Tcoeff + Ti, Fvalue);

        }
        else
        {
            FLINT_ASSERT(!nmod_poly_is_zero(Fvalue));
            nmod_poly_set(Tcoeff + Ti, Fvalue);
        }

        Fi += (Fvalue != zero);
        if (Avalue != 0)
        {
            do {
                ai--;
            } while (ai >= 0 && (Acoeff + Ai)->coeffs[ai] == 0);
            if (ai < 0)
            {
                Ai++;
                if (Ai < Alen)
                {
                    ai = nmod_poly_degree(A->coeffs + Ai);
                }
            }
        }
        if (Bvalue != 0)
        {
            do {
                bi--;
            } while (bi >= 0 && (Bcoeff + Bi)->coeffs[bi] == 0);
            if (bi < 0)
            {
                Bi++;
                if (Bi < Blen)
                {
                    bi = nmod_poly_degree(B->coeffs + Bi);
                }
            }
        }

        FLINT_ASSERT(!nmod_poly_is_zero(Tcoeff + Ti));
        lastdeg = FLINT_MAX(lastdeg, nmod_poly_degree(Tcoeff + Ti));
        Ti++;
    }
    T->length = Ti;

    nmod_poly_clear(tp);
    nmod_poly_clear(zero);

    FLINT_ASSERT(nmod_mpolyn_is_canonical(T, ctx));


    *lastdeg_ = FLINT_MAX(*lastdeg_, lastdeg);
    return changed;
}

/*
    set F from
        value F modulo modulus
        value A at x_var = alpha
        value B at x_var = -alpha
    it is expected that modulus(alpha) == modulus(-alpha) == 1/(2*alpha)
    no assumptions about matching monomials
    F is in Fp[X][x_0, ..., x_(var-1), x_(var-1)][x_var]
    A is in Fp[X][x_0, ..., x_(var-2)][x_(var-1)]
    B is in Fp[X][x_0, ..., x_(var-2)][x_(var-1)]
*/
int nmod_mpolyun_addinterp2_un(slong * lastdeg,
         nmod_mpolyun_t F, nmod_mpolyun_t T, nmod_mpolyun_t A, nmod_mpolyun_t B,
                   slong var, nmod_poly_t modulus, nmod_poly_t alphapow,
                                                    const nmod_mpoly_ctx_t ctx)
{
    int changed = 0, Finc, Ainc, Binc;
    slong Ti, Fi, Ai, Bi;
    ulong * Texp;
    ulong * Fexp;
    ulong * Aexp;
    ulong * Bexp;
    slong Flen;
    slong Alen;
    slong Blen;
    slong e;
    nmod_mpolyn_struct * Tcoeff;
    nmod_mpolyn_struct * Fcoeff;
    nmod_mpolyn_struct  * Acoeff;
    nmod_mpolyn_struct  * Bcoeff;
    nmod_mpolyn_t zero;
    nmod_mpolyn_struct * Fvalue, * Avalue, * Bvalue;

    FLINT_ASSERT(var > 0);

    *lastdeg = -WORD(1);

    FLINT_ASSERT(F->bits == T->bits);
    FLINT_ASSERT(T->bits == A->bits);
    FLINT_ASSERT(T->bits == A->bits);

    Fcoeff = F->coeffs;
    Acoeff = A->coeffs;
    Bcoeff = B->coeffs;
    Fexp = F->exps;
    Aexp = A->exps;
    Bexp = B->exps;
    Flen = F->length;
    Alen = A->length;
    Blen = B->length;

    nmod_mpolyn_init(zero, A->bits, ctx);

    FLINT_ASSERT(nmod_mpolyun_is_canonical(A, ctx));
    FLINT_ASSERT(nmod_mpolyun_is_canonical(B, ctx));


    nmod_mpolyun_fit_length(T, FLINT_MAX(Alen, Blen), ctx);
    Tcoeff = T->coeffs;
    Texp = T->exps;

    Ti = Fi = Ai = Bi = 0;
    while (Fi < Flen || Ai < Alen || Bi < Blen)
    {
        if (Ti >= T->alloc)
        {
            slong extra = Flen - Fi;
            extra = FLINT_MAX(extra, Alen - Ai);
            extra = FLINT_MAX(extra, Blen - Blen);
            nmod_mpolyun_fit_length(T, Ti + extra, ctx);
            Tcoeff = T->coeffs;
            Texp = T->exps;
        }

        if (Fi < Flen && Ai < Alen && Bi < Blen &&
                                  Fexp[Fi] == Aexp[Ai] && Aexp[Ai] == Bexp[Bi])
        {
            /* F term ok, A term ok, B term ok */
            changed |= nmod_mpolyn_addinterp2_n(lastdeg, Tcoeff + Ti,
                                      Fcoeff + Fi, Acoeff + Ai, Bcoeff + Bi,
                                                  var, modulus, alphapow, ctx);
            Texp[Ti] = Fexp[Fi];
            Fi++;
            Ai++;
            Bi++;
        }
        else
        {
            /* at least one term is missing */
            e = -WORD(1);
            if (Fi < Flen)
            {
                e = Fexp[Fi];
            }
            if (Ai < Alen)
            {
                e = FLINT_MAX(e, Aexp[Ai]);
            }
            if (Bi < Blen)
            {
                e = FLINT_MAX(e, Bexp[Bi]);
            }

            Finc = Ainc = Binc = 0;
            Fvalue = Avalue = Bvalue = zero;
            if (Fi < Flen && e == Fexp[Fi])
            {
                Finc = 1;
                Fvalue = Fcoeff + Fi;
            }
            if (Ai < Alen && e == Aexp[Ai])
            {
                Ainc = 1;
                Avalue = Acoeff + Ai;
            }
            if (Bi < Blen && e == Bexp[Bi])
            {
                Binc = 1;
                Bvalue = Bcoeff + Bi;
            }
            FLINT_ASSERT(Finc || Ainc || Binc);

            Texp[Ti] = e;
            changed |= nmod_mpolyn_addinterp2_n(lastdeg, Tcoeff + Ti,
                                                  Fvalue, Avalue, Bvalue,
                                                  var, modulus, alphapow, ctx);
            Fi += Finc;
            Ai += Ainc;
            Bi += Binc;            
        }

        FLINT_ASSERT(!nmod_mpolyn_is_zero(Tcoeff + Ti, ctx));
        Ti++;
    }
    T->length = Ti;

    if (changed)
    {
        nmod_mpolyun_swap(T, F);
    }

    nmod_mpolyn_clear(zero, ctx);

    FLINT_ASSERT(nmod_mpolyun_is_canonical(F, ctx));

/*
pthread_mutex_lock(&iomutex);
flint_printf("nmod_mpolyun_addinterp2_un(%wd) alpha = %wu\n", var, nmod_poly_get_coeff_ui(alphapow,1));
printf("A: "); nmod_mpolyun_print_pretty(A, NULL, ctx); printf("\n");
printf("B: "); nmod_mpolyun_print_pretty(B, NULL, ctx); printf("\n");
printf("F: "); nmod_mpolyun_print_pretty(F, NULL, ctx); printf("\n");
pthread_mutex_unlock(&iomutex);
*/
    return changed;
}


/*
    E = A(x_var = alpha), F = A(x_var = -alpha)
    A is in [x_0, ..., x_(var-2), x_(var-1)][x_var]
    E is in [x_0, ..., x_(var-2)][x_(var-1)]
    F is in [x_0, ..., x_(var-2)][x_(var-1)]
*/
void nmod_mpolyn_eval2_last_n(nmod_mpolyn_t E, nmod_mpolyn_t F,
                    nmod_mpolyn_t A, slong var, nmod_poly_t alphapow,
                                                    const nmod_mpoly_ctx_t ctx)
{
    slong N = mpoly_words_per_exp(A->bits, ctx->minfo);
    slong offset, shift, k;
    ulong mask;
    mp_limb_t e, f;
    nmod_poly_struct * Acoeff = A->coeffs;
    ulong * Aexp = A->exps;
    slong Alen = A->length;
    slong Ai;
    nmod_poly_struct * Ecoeff;
    ulong * Eexp;
    slong Ei;
    nmod_poly_struct * Fcoeff;
    ulong * Fexp;
    slong Fi;

    FLINT_ASSERT(var > 0);
    mpoly_gen_offset_shift_sp(&offset, &shift, var - 1, A->bits, ctx->minfo);
    mask = (-UWORD(1)) >> (FLINT_BITS - A->bits);

    Ecoeff = E->coeffs;
    Eexp = E->exps;
    Ei = 0;
    Fcoeff = F->coeffs;
    Fexp = F->exps;
    Fi = 0;
    for (Ai = 0; Ai < Alen; Ai++)
    {
        _nmod_poly_eval2_pow(&e, &f, Acoeff + Ai, alphapow, ctx->ffinfo);
        k = ((Aexp + N*Ai)[offset] >> shift) & mask;

        if (e != 0)
        {
            if (Ei > 0 && mpoly_monomial_equal_extra(Eexp + N*(Ei - 1), Aexp + N*Ai, N, offset, -(k << shift)))
            {
                /* append to previous */
                nmod_poly_set_coeff_ui(Ecoeff + Ei - 1, k, e);
            }
            else
            {
                FLINT_ASSERT(Ei == 0 || mpoly_monomial_gt_nomask_extra(Eexp + N*(Ei - 1), Aexp + N*Ai, N, offset, -(k << shift)));

                /* create new */
                if (Ei >= E->alloc)
                {
                    nmod_mpolyn_fit_length(E, Ei + 1, ctx);
                    Ecoeff = E->coeffs;
                    Eexp = E->exps;
                }
                mpoly_monomial_set_extra(Eexp + N*Ei, Aexp + N*Ai, N, offset, -(k << shift));
                nmod_poly_zero(Ecoeff + Ei);
                nmod_poly_set_coeff_ui(Ecoeff + Ei, k, e);
                Ei++;
            }
        }

        if (f != 0)
        {
            if (Fi > 0 && mpoly_monomial_equal_extra(Fexp + N*(Fi - 1), Aexp + N*Ai, N, offset, -(k << shift)))
            {
                /* append to previous */
                nmod_poly_set_coeff_ui(Fcoeff + Fi - 1, k, f);
            }
            else
            {
                FLINT_ASSERT(Fi == 0 || mpoly_monomial_gt_nomask_extra(Fexp + N*(Fi - 1), Aexp + N*Ai, N, offset, -(k << shift)));

                /* create new */
                if (Fi >= F->alloc)
                {
                    nmod_mpolyn_fit_length(F, Fi + 1, ctx);
                    Fcoeff = F->coeffs;
                    Fexp = F->exps;
                }
                mpoly_monomial_set_extra(Fexp + N*Fi, Aexp + N*Ai, N, offset, -(k << shift));
                nmod_poly_zero(Fcoeff + Fi);
                nmod_poly_set_coeff_ui(Fcoeff + Fi, k, f);
                Fi++;
            }
        }
    }
    E->length = Ei;
    F->length = Fi;
}


/*
    E = A(x_var = alpha)
    A is in R[X][x_0, ..., x_(var-2), x_(var-1)][x_var]
    E is in R[X][x_0, ..., x_(var-2)][x_(var-1)]
*/
void nmod_mpolyun_eval2_last_un(nmod_mpolyun_t E, nmod_mpolyun_t F, 
                   nmod_mpolyun_t A, slong var, nmod_poly_t alphapow,
                                                    const nmod_mpoly_ctx_t ctx)
{
    nmod_mpolyn_struct * Acoeff = A->coeffs;
    ulong * Aexp = A->exps;
    slong Alen = A->length;
    slong Ai;
    nmod_mpolyn_struct * Ecoeff;
    ulong * Eexp;
    slong Ei;
    nmod_mpolyn_struct * Fcoeff;
    ulong * Fexp;
    slong Fi;

    nmod_mpolyun_fit_length(E, Alen, ctx);
    Ecoeff = E->coeffs;
    Eexp = E->exps;

    nmod_mpolyun_fit_length(F, Alen, ctx);
    Fcoeff = F->coeffs;
    Fexp = F->exps;

    Ei = Fi = 0;
    for (Ai = 0; Ai < Alen; Ai++)
    {
        nmod_mpolyn_eval2_last_n(Ecoeff + Ei, Fcoeff + Fi, Acoeff + Ai, var, alphapow, ctx);
        Eexp[Ei] = Aexp[Ai];
        Fexp[Fi] = Aexp[Ai];
        Ei += ((Ecoeff + Ei)->length != 0);
        Fi += ((Fcoeff + Fi)->length != 0);
    }
    E->length = Ei;
    F->length = Fi;

    FLINT_ASSERT(nmod_mpolyun_is_canonical(E, ctx));
    FLINT_ASSERT(nmod_mpolyun_is_canonical(F, ctx));

/*
pthread_mutex_lock(&iomutex);
flint_printf("nmod_mpolyun_eval2_last_un(%wd) alpha = %wu\n", var, nmod_poly_get_coeff_ui(alphapow,1));
printf("A: "); nmod_mpolyun_print_pretty(A, NULL, ctx); printf("\n");
printf("E: "); nmod_mpolyun_print_pretty(E, NULL, ctx); printf("\n");
printf("F: "); nmod_mpolyun_print_pretty(F, NULL, ctx); printf("\n");
pthread_mutex_unlock(&iomutex);
*/
}





/*
    T = A,B
    T is in [x_0, ..., x_(var-1), x_(var-1)][x_var]
    A is in [x_0, ..., x_(var-2)][x_(var-1)]
*/
void nmod_mpolyn_startinterp2_n(slong * lastdeg_,
              nmod_mpolyn_t T, nmod_mpolyn_t A, nmod_mpolyn_t B, slong var,
                                   mp_limb_t alpha, const nmod_mpoly_ctx_t ctx)
{
    slong N = mpoly_words_per_exp_sp(A->bits, ctx->minfo);
    slong lastdeg = -WORD(1);
    slong offset, shift;
    nmod_poly_t tp;
    nmod_poly_t zero;

    nmod_poly_struct * Tcoeff;
    ulong * Texp;
    slong Ti;

    nmod_poly_struct * Acoeff = A->coeffs;
    slong Alen = A->length;
    ulong * Aexp = A->exps;
    slong Ai, ai;

    nmod_poly_struct * Bcoeff = B->coeffs;
    slong Blen = B->length;
    ulong * Bexp = B->exps;
    slong Bi, bi;

    mp_limb_t u, v, Avalue, Bvalue, FvalueA, FvalueB;
    int cmp;
    mp_limb_t d0 = n_invmod(alpha + alpha, ctx->ffinfo->mod.n);
/*
flint_printf("nmod_mpolyn_startinterp2_n called (p = %wu)\n", ctx->ffinfo->mod.n);
printf("T: "); nmod_mpolyn_print_pretty(T, NULL,ctx); printf("\n");
printf("A: "); nmod_mpolyn_print_pretty(A, NULL,ctx); printf("\n");
printf("B: "); nmod_mpolyn_print_pretty(B, NULL,ctx); printf("\n");
*/

    nmod_poly_init(tp, ctx->ffinfo->mod.n);
    nmod_poly_init(zero, ctx->ffinfo->mod.n);

    FLINT_ASSERT(var > 0);
    FLINT_ASSERT(T->bits == A->bits);

    mpoly_gen_offset_shift_sp(&offset, &shift, var - 1, A->bits, ctx->minfo);

    nmod_mpolyn_fit_length(T, FLINT_MAX(Alen, Blen), ctx);
    Tcoeff = T->coeffs;
    Texp = T->exps;

    Ti = Ai = Bi = 0;
    ai = (Ai >= Alen) ? 0 : nmod_poly_degree(A->coeffs + Ai);
    bi = (Bi >= Blen) ? 0 : nmod_poly_degree(B->coeffs + Bi);

    while (Ai < Alen || Bi < Blen)
    {
/*
flint_printf("Ai: %wd ai: %wd, Bi: %wd bi: %wd\n", Ai, ai, Bi, bi);
*/
        if (Ti >= T->alloc)
        {
            slong extra = FLINT_MAX(Alen - Ai, Blen - Bi);
            nmod_mpolyn_fit_length(T, Ti + extra, ctx);
            Tcoeff = T->coeffs;
            Texp = T->exps;
        }

        FLINT_ASSERT(Ai >= Alen || (Acoeff + Ai)->coeffs[ai] != 0);
        FLINT_ASSERT(Bi >= Blen || (Bcoeff + Bi)->coeffs[bi] != 0);

        Avalue = 0;
        if (Ai < Alen)
        {
            Avalue = (Acoeff + Ai)->coeffs[ai];
            mpoly_monomial_set_extra(Texp + N*Ti,
                                     Aexp + N*Ai, N, offset, ai << shift);
        }

        Bvalue = 0;
        if (Bi < Blen)
        {
            cmp = (Avalue == 0) ? -1 : mpoly_monomial_cmp_nomask_extra(
                             Texp + N*Ti, Bexp + N*Bi, N, offset, bi << shift);
            if (cmp <= 0)
            {
                Bvalue = (Bcoeff + Bi)->coeffs[bi];
            }
            if (cmp < 0)
            {
                Avalue = 0;
                mpoly_monomial_set_extra(Texp + N*Ti,
                                         Bexp + N*Bi, N, offset, bi << shift);
            }
        }

        FvalueA = nmod_neg(Avalue, ctx->ffinfo->mod);
        FvalueB = nmod_neg(Bvalue, ctx->ffinfo->mod);
        u = nmod_sub(FvalueB, FvalueA, ctx->ffinfo->mod);
        v = nmod_mul(ctx->ffinfo->mod.n - alpha, nmod_add(FvalueB, FvalueA, ctx->ffinfo->mod), ctx->ffinfo->mod);
/*
flint_printf("u: %wu    v: %wu\n",u,v);
*/
        FLINT_ASSERT(u != 0 || v != 0);
        nmod_poly_zero(Tcoeff + Ti);
        u = nmod_mul(u, d0, ctx->ffinfo->mod);
        v = nmod_mul(v, d0, ctx->ffinfo->mod);
        nmod_poly_set_coeff_ui(Tcoeff + Ti, 0, v);
        nmod_poly_set_coeff_ui(Tcoeff + Ti, 1, u);
/*
flint_printf("u: %wu    v: %wu\n",u,v);
printf("Tcoeff+Ti: "); nmod_poly_print_pretty(Tcoeff + Ti, "v"); printf("\n");
*/
        if (Avalue != 0)
        {
            do {
                ai--;
            } while (ai >= 0 && (Acoeff + Ai)->coeffs[ai] == 0);
            if (ai < 0)
            {
                Ai++;
                if (Ai < Alen)
                {
                    ai = nmod_poly_degree(A->coeffs + Ai);
                }
            }
        }
        if (Bvalue != 0)
        {
            do {
                bi--;
            } while (bi >= 0 && (Bcoeff + Bi)->coeffs[bi] == 0);
            if (bi < 0)
            {
                Bi++;
                if (Bi < Blen)
                {
                    bi = nmod_poly_degree(B->coeffs + Bi);
                }
            }
        }

        FLINT_ASSERT(!nmod_poly_is_zero(Tcoeff + Ti));
        lastdeg = FLINT_MAX(lastdeg, nmod_poly_degree(Tcoeff + Ti));
        Ti++;
    }
    T->length = Ti;

    FLINT_ASSERT(nmod_mpolyn_is_canonical(T, ctx));

    *lastdeg_ = FLINT_MAX(*lastdeg_, lastdeg);
    return;
}


/*
    set F from
        A at x_var = alpha
        B at x_var = -alpha
    F is in R[X][x_0, ..., x_(var-2), x_(var-1)][x_var]
    A is in R[X][x_0, ..., x_(var-2)][x_(var-1)]
    B is in R[X][x_0, ..., x_(var-2)][x_(var-1)]
*/
void nmod_mpolyun_startinterp2_un(slong * lastdeg, nmod_mpolyun_t F,
                             nmod_mpolyun_t A, nmod_mpolyun_t B, slong var,
                                  mp_limb_t alpha, const nmod_mpoly_ctx_t ctx)
{
    nmod_mpolyn_struct * Bcoeff = B->coeffs;
    ulong * Bexp = B->exps;
    slong Blen = B->length;
    slong Bi;

    nmod_mpolyn_struct * Acoeff = A->coeffs;
    ulong * Aexp = A->exps;
    slong Alen = A->length;
    slong Ai;

    nmod_mpolyn_struct * Fcoeff;
    ulong * Fexp;
    slong Fi;

    nmod_mpolyn_t zero;

    nmod_mpolyn_init(zero, A->bits, ctx);

    FLINT_ASSERT(F->bits == A->bits);
    FLINT_ASSERT(F->bits == B->bits);

    nmod_mpolyun_fit_length(F, Alen + Blen, ctx);
    Fcoeff = F->coeffs;
    Fexp = F->exps;

    *lastdeg = -WORD(1);

    FLINT_ASSERT(nmod_mpolyun_is_canonical(A, ctx));
    FLINT_ASSERT(nmod_mpolyun_is_canonical(B, ctx));
/*
printf("nmod_mpolyun_startinterp2_un called\n");
printf("F: "); nmod_mpolyun_print_pretty(F, NULL,ctx); printf("\n");
printf("A: "); nmod_mpolyun_print_pretty(A, NULL,ctx); printf("\n");
printf("B: "); nmod_mpolyun_print_pretty(B, NULL,ctx); printf("\n");
*/

    Fi = Ai = Bi = 0;
    while (Ai < Alen || Bi < Blen)
    {
        FLINT_ASSERT(Fi < F->alloc);

        if (Ai < Alen && Bi < Blen && Aexp[Ai] == Bexp[Bi])
        {
            Fexp[Fi] = Aexp[Ai];
            nmod_mpolyn_startinterp2_n(lastdeg, Fcoeff + Fi, Acoeff + Ai, Bcoeff + Bi, var, alpha, ctx);
            Ai++;
            Bi++;
        }
        else if (Ai < Alen && (Bi >= Blen || Aexp[Ai] > Bexp[Bi]))
        {
            Fexp[Fi] = Aexp[Ai];
            nmod_mpolyn_startinterp2_n(lastdeg, Fcoeff + Fi, Acoeff + Ai, zero, var, alpha, ctx);
            Ai++;
        }
        else
        {
            FLINT_ASSERT(Bi < Blen && (Ai >= Alen || Bexp[Bi] > Aexp[Ai]));

            Fexp[Fi] = Bexp[Bi];
            nmod_mpolyn_startinterp2_n(lastdeg, Fcoeff + Fi, zero, Bcoeff + Bi, var, alpha, ctx);
            Bi++;
        }

        FLINT_ASSERT((Fcoeff + Fi)->length > 0);
        Fi++;
    }
    F->length = Fi;

    nmod_mpolyn_clear(zero, ctx);

    FLINT_ASSERT(nmod_mpolyun_is_canonical(F, ctx));

/*
pthread_mutex_lock(&iomutex);
flint_printf("nmod_mpolyun_startinterp2_un(%wd) alpha = %wu\n", var, alpha);
printf("A: "); nmod_mpolyun_print_pretty(A, NULL, ctx); printf("\n");
printf("B: "); nmod_mpolyun_print_pretty(B, NULL, ctx); printf("\n");
printf("F: "); nmod_mpolyun_print_pretty(F, NULL, ctx); printf("\n");
pthread_mutex_unlock(&iomutex);
*/

}














/*
    E = A(v = alpha)
    A is in R[X][v]
    E is in R[X]
*/
void nmod_mpolyun_eval_last_bivar(nmod_poly_t E, const nmod_mpolyun_t A,
                                  mp_limb_t alpha, const nmod_mpoly_ctx_t ctx)
{
    mp_limb_t v;
    slong Ai, Alen;
    nmod_mpolyn_struct * Acoeff;
    ulong * Aexp;

    Acoeff = A->coeffs;
    Aexp = A->exps;
    Alen = A->length;
    Ai = 0;
    nmod_poly_zero(E);
    for (Ai = 0; Ai < Alen; Ai++)
    {
        v = nmod_poly_evaluate_nmod((Acoeff + Ai)->coeffs + 0, alpha);
        nmod_poly_set_coeff_ui(E, Aexp[Ai], v);
    }
}

/*
    A = B
    A, B are in R[X]
*/
void nmod_mpolyun_set_poly(nmod_mpolyun_t A, const nmod_poly_t B,
                                                   const nmod_mpoly_ctx_t ctx)
{
    slong N = mpoly_words_per_exp_sp(A->bits, ctx->minfo);
    slong Bexp;
    slong Blen = nmod_poly_length(B);
    mp_limb_t * Bcoeff = B->coeffs;
    nmod_mpolyn_struct * Acoeff;
    ulong * Aexp;
    slong Ai;

    nmod_mpolyun_fit_length(A, Blen, ctx);
    Acoeff = A->coeffs;
    Aexp = A->exps;

    Ai = 0;
    for (Bexp = Blen - 1; Bexp >= 0; Bexp--)
    {
        if (Bcoeff[Bexp] != UWORD(0))
        {
            FLINT_ASSERT(Ai < A->alloc);

            nmod_mpolyn_fit_length(Acoeff + Ai, 1, ctx);
            mpoly_monomial_zero((Acoeff + Ai)->exps + N*0, N);
            nmod_poly_zero((Acoeff + Ai)->coeffs + 0);
            nmod_poly_set_coeff_ui((Acoeff + Ai)->coeffs + 0, 0, Bcoeff[Bexp]);
            Aexp[Ai] = Bexp;
            (Acoeff + Ai)->length = 1;
            Ai++;
        }
    }
    A->length = Ai;
}






/*
    T = F + modulus*(A - F(x_var = alpha))
    no assumptions about matching monomials
    F is in Fp[x_0, ..., x_(var-1), x_(var-1)][x_var]
    A is in Fp[x_0, ..., x_(var-2)][x_(var-1)]
    in order to fxn correctly, modulus(alpha) should be 1
*/
int nmod_mpolyn_addinterp_n(slong * lastdeg_,
             nmod_mpolyn_t T, nmod_mpolyn_t F, nmod_mpolyn_t A, slong var,
              nmod_poly_t modulus, mp_limb_t alpha, const nmod_mpoly_ctx_t ctx)
{
    int changed = 0;
    slong N = mpoly_words_per_exp_sp(A->bits, ctx->minfo);
    slong lastdeg = -WORD(1);
    slong offset, shift;
    slong vi;
    mp_limb_t v;
    nmod_poly_t tp;

    nmod_poly_struct * Tcoeff;
    ulong * Texp;
    slong Ti;

    nmod_poly_struct * Acoeff = A->coeffs;
    slong Alen = A->length;
    ulong * Aexp = A->exps;
    slong Ai;

    nmod_poly_struct * Fcoeff = F->coeffs;
    slong Flen = F->length;
    ulong * Fexp = F->exps;
    slong Fi;

    nmod_poly_init(tp, ctx->ffinfo->mod.n);

    FLINT_ASSERT(var > 0);
    FLINT_ASSERT(T->bits == A->bits);
    FLINT_ASSERT(F->bits == A->bits);

    mpoly_gen_offset_shift_sp(&offset, &shift, var - 1, A->bits, ctx->minfo);

    Flen = F->length;

    nmod_mpolyn_fit_length(T, FLINT_MAX(Flen, Alen), ctx);
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
            nmod_mpolyn_fit_length(T, Ti + FLINT_MAX(Flen - Fi, Alen - Ai), ctx);
            Tcoeff = T->coeffs;
            Texp = T->exps;
        }

        if (Fi < Flen && Ai < Alen && mpoly_monomial_equal_extra(Fexp + N*Fi, Aexp + N*Ai, N, offset, vi << shift))
        {
            /* F term ok, A term ok */
            v = nmod_poly_evaluate_nmod(Fcoeff + Fi, alpha);
            v = nmod_sub((Acoeff + Ai)->coeffs[vi], v, ctx->ffinfo->mod);
            if (v != UWORD(0))
            {
                changed = 1;
                nmod_poly_scalar_mul_nmod(tp, modulus, v);
                nmod_poly_add(Tcoeff + Ti, Fcoeff + Fi, tp);
            }
            else
            {
                nmod_poly_set(Tcoeff + Ti, Fcoeff + Fi);
            }
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
            v = nmod_poly_evaluate_nmod(Fcoeff + Fi, alpha);
            if (v != UWORD(0))
            {
                changed = 1;
                nmod_poly_scalar_mul_nmod(tp, modulus, v);
                nmod_poly_sub(Tcoeff + Ti, Fcoeff + Fi, tp);
            }
            else
            {
                nmod_poly_set(Tcoeff + Ti, Fcoeff + Fi);
            }
            mpoly_monomial_set(Texp + N*Ti, Fexp + N*Fi, N);

            Fi++;
        }
        else
        {
            FLINT_ASSERT(Ai < Alen && (Fi >= Flen || mpoly_monomial_lt_nomask_extra(Fexp + N*Fi, Aexp + N*Ai, N, offset, vi << shift)));

            /* F term missing, A term ok */
            changed = 1;
            nmod_poly_scalar_mul_nmod(Tcoeff + Ti, modulus, (Acoeff + Ai)->coeffs[vi]);
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

        FLINT_ASSERT(!nmod_poly_is_zero(Tcoeff + Ti));
        lastdeg = FLINT_MAX(lastdeg, nmod_poly_degree(Tcoeff + Ti));
        Ti++;
    }
    T->length = Ti;

    nmod_poly_clear(tp);

    *lastdeg_ = FLINT_MAX(*lastdeg_, lastdeg);
    return changed;
}

/*
    F = F + modulus*(A - F(x_var = alpha))
    no assumptions about matching monomials
    F is in Fp[X][x_0, ..., x_(var-1), x_(var-1)][x_var]
    A is in Fp[X][x_0, ..., x_(var-2)][x_(var-1)]
    in order to fxn correctly, modulus(alpha) should be 1
*/
int nmod_mpolyun_addinterp_un(slong * lastdeg,
             nmod_mpolyun_t F, nmod_mpolyun_t T, nmod_mpolyun_t A, slong var,
              nmod_poly_t modulus, mp_limb_t alpha, const nmod_mpoly_ctx_t ctx)
{
    int changed = 0;
    slong i, j, k;
    ulong * Texp;
    ulong * Fexp;
    ulong * Aexp;
    slong Flen;
    slong Alen;
    nmod_mpolyn_struct * Tcoeff;
    nmod_mpolyn_struct * Fcoeff;
    nmod_mpolyn_struct  * Acoeff;
    nmod_mpolyn_t zero;

    FLINT_ASSERT(var > 0);

    *lastdeg = -WORD(1);

    FLINT_ASSERT(F->bits == T->bits);
    FLINT_ASSERT(T->bits == A->bits);

    Flen = F->length;
    Alen = A->length;
    nmod_mpolyun_fit_length(T, Flen + Alen, ctx);

    Tcoeff = T->coeffs;
    Fcoeff = F->coeffs;
    Acoeff = A->coeffs;
    Texp = T->exps;
    Fexp = F->exps;
    Aexp = A->exps;   

    nmod_mpolyn_init(zero, A->bits, ctx);

    i = j = k = 0;
    while (i < Flen || j < Alen)
    {
        if (i < Flen && j < Alen && (Fexp[i] == Aexp[j]))
        {
            /* F term ok, A term ok */
            changed |= nmod_mpolyn_addinterp_n(lastdeg, Tcoeff + k, Fcoeff + i,
                                         Acoeff + j, var, modulus, alpha, ctx);
            Texp[k] = Aexp[j];
            i++;
            j++;
        }
        else if (i < Flen && (j >= Alen || Fexp[i] > Aexp[j]))
        {
            /* F term ok, A term missing */
            changed |= nmod_mpolyn_addinterp_n(lastdeg, Tcoeff + k, Fcoeff + i,
                                               zero, var, modulus, alpha, ctx);
            Texp[k] = Fexp[i];
            i++;
        }
        else
        {
            FLINT_ASSERT(j < Alen && (i >= Flen || Aexp[j] > Fexp[i]));

            /* F term missing, A term ok */
            changed |= nmod_mpolyn_addinterp_n(lastdeg, Tcoeff + k, zero,
                                         Acoeff + j, var, modulus, alpha, ctx);
            Texp[k] = Aexp[j];
            j++;
        }

        FLINT_ASSERT(!nmod_mpolyn_is_zero(Tcoeff + k, ctx));
        k++;
    }
    T->length = k;

    if (changed)
    {
        nmod_mpolyun_swap(T, F);
    }

    nmod_mpolyn_clear(zero, ctx);

    FLINT_ASSERT(nmod_mpolyun_is_canonical(F, ctx));


    return changed;    
}


/*
    A = B
    A is in Fp[x_0, ..., x_(var-1), x_(var-1)][x_var]
    B is in Fp[x_0, ..., x_(var-2)][x_(var-1)]
*/
void nmod_mpolyn_set_popup(nmod_mpolyn_t A, nmod_mpolyn_t B, slong var,
                                                    const nmod_mpoly_ctx_t ctx)
{
    slong N = mpoly_words_per_exp(B->bits, ctx->minfo);
    slong offset, shift;
    slong vi;

    nmod_poly_struct * Bcoeff = B->coeffs;
    ulong * Bexp = B->exps;
    slong Blen = B->length;
    slong Bi;

    nmod_poly_struct * Acoeff;
    ulong * Aexp;
    slong Ai;

    nmod_mpolyn_fit_length(A, Blen, ctx);
    Acoeff = A->coeffs;
    Aexp = A->exps;

    mpoly_gen_offset_shift_sp(&offset, &shift, var - 1, A->bits, ctx->minfo);

    Ai = 0;
    for (Bi = 0; Bi < Blen; Bi++)
    {
        if (Ai + (Bcoeff + Bi)->length >= A->alloc)
        {
            nmod_mpolyn_fit_length(A, Ai + (Bcoeff + Bi)->length, ctx);
            Acoeff = A->coeffs;
            Aexp = A->exps;
        }
        for (vi = (Bcoeff + Bi)->length - 1; vi >= 0; vi--)
        {
            if ((Bcoeff + Bi)->coeffs[vi] != 0)
            {
                mpoly_monomial_set_extra(Aexp + N*Ai, Bexp + N*Bi, N, offset, vi << shift);
                nmod_poly_zero(Acoeff + Ai);
                nmod_poly_set_coeff_ui(Acoeff + Ai, 0, (Bcoeff + Bi)->coeffs[vi]);
                Ai++;
            }
        }
    }
    A->length = Ai;
}


/*
    A = B
    A is in Fp[X][x_0, ..., x_(var-2), x_(var-1)][x_var]
    B is in Fp[X][x_0, ..., x_(var-2)][x_(var-1)]
*/
void nmod_mpolyun_set_popup(nmod_mpolyun_t A, nmod_mpolyun_t B, slong var,
                                                    const nmod_mpoly_ctx_t ctx)
{
    slong i;
    nmod_mpolyn_struct * Bcoeff = B->coeffs;
    ulong * Bexp = B->exps;
    slong Blen = B->length;

    nmod_mpolyn_struct * Acoeff;
    ulong * Aexp;

    nmod_mpolyun_fit_length(A, Blen, ctx);
    Acoeff = A->coeffs;
    Aexp = A->exps;

    for (i = 0; i < Blen; i++)
    {
        Aexp[i] = Bexp[i];
        nmod_mpolyn_set_popup(Acoeff + i, Bcoeff + i, var, ctx);
    }
    A->length = Blen;

    FLINT_ASSERT(nmod_mpolyun_is_canonical(A, ctx));

}

/*
    E = A(x_var = alpha)
    A is in Fp[x_0, ..., x_(var-2), x_(var-1)][x_var]
    E is in Fp[x_0, ..., x_(var-2)][x_(var-1)]
*/
void nmod_mpolyn_eval_last_n(nmod_mpolyn_t E, nmod_mpolyn_t A, slong var,
                                   mp_limb_t alpha, const nmod_mpoly_ctx_t ctx)
{
    slong N = mpoly_words_per_exp_sp(A->bits, ctx->minfo);
    slong offset, shift, k;
    ulong mask;
    mp_limb_t v;
    nmod_poly_struct * Acoeff = A->coeffs;
    ulong * Aexp = A->exps;
    slong Alen = A->length;
    slong Ai;
    nmod_poly_struct * Ecoeff;
    ulong * Eexp;
    slong Ei;

    FLINT_ASSERT(var > 0);
    mpoly_gen_offset_shift_sp(&offset, &shift, var - 1, A->bits, ctx->minfo);
    mask = (-UWORD(1)) >> (FLINT_BITS - A->bits);

    Ecoeff = E->coeffs;
    Eexp = E->exps;
    Ei = 0;
    for (Ai = 0; Ai < Alen; Ai++)
    {
        v = nmod_poly_evaluate_nmod(Acoeff + Ai, alpha);
        k = ((Aexp + N*Ai)[offset] >> shift) & mask;
        if (v == 0)
        {
            continue;
        }

        if (Ei > 0 && mpoly_monomial_equal_extra(Eexp + N*(Ei - 1), Aexp + N*Ai, N, offset, -(k << shift)))
        {
            /* append to previous */
            nmod_poly_set_coeff_ui(Ecoeff + Ei - 1, k, v);
        }
        else
        {
            FLINT_ASSERT(Ei == 0 || mpoly_monomial_gt_nomask_extra(Eexp + N*(Ei - 1), Aexp + N*Ai, N, offset, -(k << shift)));

            /* create new */
            if (Ei >= E->alloc)
            {
                nmod_mpolyn_fit_length(E, Ei + 1, ctx);
                Ecoeff = E->coeffs;
                Eexp = E->exps;
            }
            mpoly_monomial_set_extra(Eexp + N*Ei, Aexp + N*Ai, N, offset, -(k << shift));
            nmod_poly_zero(Ecoeff + Ei);
            nmod_poly_set_coeff_ui(Ecoeff + Ei, k, v);
            Ei++;
        }
    }
    E->length = Ei;
}


/*
    E = A(x_var = alpha)
    A is in Fp[X][x_0, ..., x_(var-2), x_(var-1)][x_var]
    E is in Fp[X][x_0, ..., x_(var-2)][x_(var-1)]
*/
void nmod_mpolyun_eval_last_un(nmod_mpolyun_t E, nmod_mpolyun_t A, slong var,
                                   mp_limb_t alpha, const nmod_mpoly_ctx_t ctx)
{
    nmod_mpolyn_struct * Acoeff = A->coeffs;
    ulong * Aexp = A->exps;
    slong Alen = A->length;
    slong Ai;
    nmod_mpolyn_struct * Ecoeff;
    ulong * Eexp;
    slong Ei;

    nmod_mpolyun_fit_length(E, Alen, ctx);
    Ecoeff = E->coeffs;
    Eexp = E->exps;

    Ei = 0;
    for (Ai = 0; Ai < Alen; Ai++)
    {
        nmod_mpolyn_eval_last_n(Ecoeff + Ei, Acoeff + Ai, var, alpha, ctx);
        Eexp[Ei] = Aexp[Ai];
        Ei += ((Ecoeff + Ei)->length != 0);
    }
    E->length = Ei;

    FLINT_ASSERT(nmod_mpolyun_is_canonical(E, ctx));
}



/*
    get the leading exponent in x_0,...,x_var
    A is in R[x_0, ... x_(var-1)][x_var]
*/
mp_limb_t nmod_mpolyn_leadcoeff_last(nmod_mpolyn_t A, const nmod_mpoly_ctx_t ctx)
{
    nmod_poly_struct * leadpoly;

    FLINT_ASSERT(A->length > 0);
    FLINT_ASSERT(nmod_poly_degree(A->coeffs + 0) >= 0);

    leadpoly = A->coeffs + 0;
    return leadpoly->coeffs[leadpoly->length - 1];
}

mp_limb_t nmod_mpolyun_leadcoeff_last(nmod_mpolyun_t A, const nmod_mpoly_ctx_t ctx)
{
    FLINT_ASSERT(A->length > 0);
    return nmod_mpolyn_leadcoeff_last(A->coeffs + 0, ctx);
}

/*
    get the leading exponent in x_0,...,x_var
    A is in R[x_0, ... x_(var-1)][x_var]
*/
/*
void nmod_mpolyn_leadmonomial(ulong * exp, nmod_mpolyn_t A, slong var, const nmod_mpoly_ctx_t ctx)
{
    slong N = mpoly_words_per_exp_sp(A->bits, ctx->minfo);
    slong offset, shift;

    FLINT_ASSERT(A->length > 0);
    FLINT_ASSERT(nmod_poly_degree(A->coeffs + 0) >= 0);

    mpoly_monomial_set(exp, A->exps + N*0, N);
    mpoly_gen_offset_shift_sp(&offset, &shift, var, A->bits, ctx->minfo);
    exp[offset] += nmod_poly_degree(A->coeffs + 0) << shift;
}

int nmod_mpolyn_leadmonomial_cmp(nmod_mpolyn_t A, slong var, ulong * rhs, const nmod_mpoly_ctx_t ctx)
{
    int r;
    slong N = mpoly_words_per_exp_sp(A->bits, ctx->minfo);
    ulong * lhs;
    TMP_INIT;
    TMP_START;
    lhs = TMP_ALLOC(N*sizeof(ulong));
    nmod_mpolyn_leadmonomial(lhs, A, var, ctx);

flint_printf("lhs: %wx  rhs: %wx\n", lhs[0], rhs[0]);

    r = mpoly_monomial_cmp_nomask(lhs, rhs, N);

flint_printf("leadmonomial_cmp returning %d\n", r);

    TMP_END;
    return r;
}
*/

int nmod_mpolyn_is_nonzero_nmod(const nmod_mpolyn_t A, const nmod_mpoly_ctx_t ctx)
{
    slong N;

    if (A->length != WORD(1))
    {
        return 0;
    }

    if (nmod_poly_degree(A->coeffs + 0) != 0)
    {
        return 0;
    }

    N = mpoly_words_per_exp(A->bits, ctx->minfo);
    return mpoly_monomial_is_zero(A->exps + N*0, N);
}

int nmod_mpolyun_is_nonzero_nmod(const nmod_mpolyun_t A, const nmod_mpoly_ctx_t ctx)
{
    if (A->length != 1 || A->exps[0] != UWORD(0))
    {
        return 0;
    }

    return nmod_mpolyn_is_nonzero_nmod(A->coeffs + 0, ctx);
}

void nmod_mpolyn_scalar_mul_nmod(nmod_mpolyn_t A, mp_limb_t c, const nmod_mpoly_ctx_t ctx)
{
    slong i;
    for (i = 0; i < A->length; i++)
    {
        nmod_poly_scalar_mul_nmod(A->coeffs + i, A->coeffs + i, c);
    }
}

void nmod_mpolyun_scalar_mul_nmod(nmod_mpolyun_t A, mp_limb_t c, const nmod_mpoly_ctx_t ctx)
{
    slong i;
    FLINT_ASSERT(c != 0);
    for (i = 0; i < A->length; i++)
    {
        nmod_mpolyn_scalar_mul_nmod(A->coeffs + i, c, ctx);
    }
}



/***********************************************
    gcd
**************************************************/



int nmod_mpolyun_gcd_brown_smprime_bivar_ref(nmod_mpolyun_t G,
  nmod_mpolyun_t Abar, nmod_mpolyun_t Bbar, nmod_mpolyun_t A, nmod_mpolyun_t B,
                       const nmod_mpoly_ctx_t ctx, nmod_poly_mpolyun_stack_t S)
{
    int success;
    slong bound;
    mp_limb_t alpha, temp, gammaeval;
    nmod_poly_t Aeval, Beval, Geval, Abareval, Bbareval;
    nmod_mpolyun_t T;
    slong deggamma, ldegG, ldegAbar, ldegBbar, ldegA, ldegB;
    nmod_poly_t cA, cB, cG, cAbar, cBbar, gamma;
    nmod_poly_t modulus, modulus2;
#if WANT_ASSERT
    nmod_poly_t leadA, leadB;
#endif

#if WANT_ASSERT
    nmod_poly_init(leadA, ctx->ffinfo->mod.n);
    nmod_poly_init(leadB, ctx->ffinfo->mod.n);
    nmod_poly_set(leadA, nmod_mpolyun_leadcoeff_ref(A, ctx));
    nmod_poly_set(leadB, nmod_mpolyun_leadcoeff_ref(B, ctx));
#endif

    nmod_poly_init(cA, ctx->ffinfo->mod.n);
    nmod_poly_init(cB, ctx->ffinfo->mod.n);
    nmod_mpolyun_content_last(cA, A, ctx);
    nmod_mpolyun_content_last(cB, B, ctx);
    nmod_mpolyun_divexact_last(A, cA, ctx);
    nmod_mpolyun_divexact_last(B, cB, ctx);

    nmod_poly_init(cG, ctx->ffinfo->mod.n);
    nmod_poly_gcd_euclidean(cG, cA, cB);

    nmod_poly_init(cAbar, ctx->ffinfo->mod.n);
    nmod_poly_init(cBbar, ctx->ffinfo->mod.n);
    nmod_poly_div(cAbar, cA, cG);
    nmod_poly_div(cBbar, cB, cG);

    nmod_poly_init(gamma, ctx->ffinfo->mod.n);
    nmod_poly_gcd(gamma, nmod_mpolyun_leadcoeff_ref(A, ctx),
                         nmod_mpolyun_leadcoeff_ref(B, ctx));

    ldegA = nmod_mpolyun_lastdeg(A, ctx);
    ldegB = nmod_mpolyun_lastdeg(B, ctx);
    deggamma = nmod_poly_degree(gamma);
    bound = 1 + deggamma + FLINT_MAX(ldegA, ldegB);

    nmod_poly_init(Aeval, ctx->ffinfo->mod.n);
    nmod_poly_init(Beval, ctx->ffinfo->mod.n);
    nmod_poly_init(Geval, ctx->ffinfo->mod.n);
    nmod_poly_init(Abareval, ctx->ffinfo->mod.n);
    nmod_poly_init(Bbareval, ctx->ffinfo->mod.n);

    nmod_mpolyun_init(T, A->bits, ctx);

    nmod_poly_init(modulus, ctx->ffinfo->mod.n);
    nmod_poly_init(modulus2, ctx->ffinfo->mod.n);
    nmod_poly_one(modulus);

    if ((ctx->ffinfo->mod.n & UWORD(1)) == UWORD(0))
    {
        success = 0;
        goto cleanup;
    }

    alpha = (ctx->ffinfo->mod.n - UWORD(1))/UWORD(2);

choose_prime: /* prime is v - alpha */

    if (alpha < 2)
    {
        success = 0;
        goto cleanup;
    }

    alpha--;

    /* make sure evaluation point does not kill both lc(A) and lc(B) */
    gammaeval = nmod_poly_evaluate_nmod(gamma, alpha);
    if (gammaeval == 0)
    {
        goto choose_prime;
    }

    /* evaluation point should kill neither A nor B */
    nmod_mpolyun_eval_last_bivar(Aeval, A, alpha, ctx);
    nmod_mpolyun_eval_last_bivar(Beval, B, alpha, ctx);
    FLINT_ASSERT(Aeval->length > 0);
    FLINT_ASSERT(Beval->length > 0);

    nmod_poly_gcd(Geval, Aeval, Beval);
    nmod_poly_div(Abareval, Aeval, Geval);
    nmod_poly_div(Bbareval, Beval, Geval);

    FLINT_ASSERT(Geval->length > 0);
    FLINT_ASSERT(Abareval->length > 0);
    FLINT_ASSERT(Bbareval->length > 0);

    if (nmod_poly_degree(Geval) == 0)
    {
        nmod_mpolyun_one(G, ctx);
        nmod_mpolyun_swap(Abar, A);
        nmod_mpolyun_swap(Bbar, B);
        goto successful_put_content;    
    }

    if (nmod_poly_degree(modulus) > 0)
    {
        FLINT_ASSERT(G->length > 0);
        if (nmod_poly_degree(Geval) > G->exps[0])
        {
            goto choose_prime;
        }
        else if (nmod_poly_degree(Geval) < G->exps[0])
        {
            nmod_poly_one(modulus);
        }
    }

    nmod_poly_scalar_mul_nmod(Geval, Geval, gammaeval);

    if (nmod_poly_degree(modulus) > 0)
    {
        temp = nmod_poly_evaluate_nmod(modulus, alpha);
        temp = n_invmod(temp, ctx->ffinfo->mod.n);
        nmod_poly_scalar_mul_nmod(modulus, modulus, temp);
        nmod_mpolyun_addinterp_bivar(&ldegG, G, T, Geval, modulus, alpha, ctx);
        nmod_mpolyun_addinterp_bivar(&ldegAbar, Abar, T, Abareval, modulus, alpha, ctx);
        nmod_mpolyun_addinterp_bivar(&ldegBbar, Bbar, T, Bbareval, modulus, alpha, ctx);
    }
    else
    {
        nmod_mpolyun_set_poly(G, Geval, ctx);
        nmod_mpolyun_set_poly(Abar, Abareval, ctx);
        nmod_mpolyun_set_poly(Bbar, Bbareval, ctx);
        ldegG = 0;
        ldegAbar = 0;
        ldegBbar = 0;
    }

    nmod_poly_scalar_mul_nmod(modulus2, modulus, alpha);
    nmod_poly_shift_left(modulus, modulus, 1);
    nmod_poly_sub(modulus, modulus, modulus2);

    if (nmod_poly_degree(modulus) < bound)
    {
        goto choose_prime;
    }

    FLINT_ASSERT(ldegG >= 0);
    FLINT_ASSERT(ldegAbar >= 0);
    FLINT_ASSERT(ldegBbar >= 0);

    if (   deggamma + ldegA == ldegG + ldegAbar
        && deggamma + ldegB == ldegG + ldegBbar )
    {
        goto successful;
    }

    nmod_poly_one(modulus);
    goto choose_prime;

successful:

    nmod_mpolyun_content_last(modulus, G, ctx);
    nmod_mpolyun_divexact_last(G, modulus, ctx);
    nmod_mpolyun_divexact_last(Abar, nmod_mpolyun_leadcoeff_ref(G, ctx), ctx);
    nmod_mpolyun_divexact_last(Bbar, nmod_mpolyun_leadcoeff_ref(G, ctx), ctx);

successful_put_content:

    nmod_mpolyun_mul_last(G, cG, ctx);
    nmod_mpolyun_mul_last(Abar, cAbar, ctx);
    nmod_mpolyun_mul_last(Bbar, cBbar, ctx);

    success = 1;

cleanup:

#if WANT_ASSERT
    if (success)
    {
        FLINT_ASSERT(1 == nmod_mpolyun_leadcoeff_last(G, ctx));
        nmod_poly_mul(modulus, nmod_mpolyun_leadcoeff_ref(G, ctx),
                               nmod_mpolyun_leadcoeff_ref(Abar, ctx));
        FLINT_ASSERT(nmod_poly_equal(modulus, leadA));
        nmod_poly_mul(modulus, nmod_mpolyun_leadcoeff_ref(G, ctx),
                               nmod_mpolyun_leadcoeff_ref(Bbar, ctx));
        FLINT_ASSERT(nmod_poly_equal(modulus, leadB));
    }
    nmod_poly_clear(leadA);
    nmod_poly_clear(leadB);
#endif

    nmod_poly_clear(cA);
    nmod_poly_clear(cB);
    nmod_poly_clear(cG);
    nmod_poly_clear(cAbar);
    nmod_poly_clear(cBbar);

    nmod_poly_clear(gamma);

    nmod_poly_clear(Aeval);
    nmod_poly_clear(Beval);
    nmod_poly_clear(Geval);
    nmod_poly_clear(Abareval);
    nmod_poly_clear(Bbareval);

    nmod_mpolyun_clear(T, ctx);

    nmod_poly_clear(modulus);
    nmod_poly_clear(modulus2);

    return success;
}



int nmod_mpolyun_gcd_brown_smprime_ref(nmod_mpolyun_t G,
                                nmod_mpolyun_t Abar, nmod_mpolyun_t Bbar,
                                    nmod_mpolyun_t A, nmod_mpolyun_t B,
                                         slong var, const nmod_mpoly_ctx_t ctx)
{
    int success;
    slong bound;
    slong offset, shift;
    mp_limb_t alpha, temp, gammaeval;
    nmod_mpolyun_t Aeval, Beval, Geval, Abareval, Bbareval;
    nmod_mpolyun_t T;
    slong deggamma, ldegG, ldegAbar, ldegBbar, ldegA, ldegB;
    nmod_poly_t cA, cB, cG, cAbar, cBbar, gamma;
    nmod_poly_t modulus, modulus2;
    mp_bitcnt_t bits = A->bits;
    slong N = mpoly_words_per_exp(bits, ctx->minfo);
#if WANT_ASSERT
    nmod_poly_t leadA, leadB;
#endif

    FLINT_ASSERT(var >= 0);
    if (var == 0)
    {
        /* bivariate is more comfortable separated */
        return nmod_mpolyun_gcd_brown_smprime_bivar_ref(G, Abar, Bbar, A, B, ctx);
    }

    mpoly_gen_offset_shift_sp(&offset, &shift, var - 1, G->bits, ctx->minfo);

#if WANT_ASSERT
    nmod_poly_init(leadA, ctx->ffinfo->mod.n);
    nmod_poly_init(leadB, ctx->ffinfo->mod.n);
    nmod_poly_set(leadA, nmod_mpolyun_leadcoeff_ref(A, ctx));
    nmod_poly_set(leadB, nmod_mpolyun_leadcoeff_ref(B, ctx));
#endif

    nmod_poly_init(cA, ctx->ffinfo->mod.n);
    nmod_poly_init(cB, ctx->ffinfo->mod.n);
    nmod_mpolyun_content_last(cA, A, ctx);
    nmod_mpolyun_content_last(cB, B, ctx);
    nmod_mpolyun_divexact_last(A, cA, ctx);
    nmod_mpolyun_divexact_last(B, cB, ctx);

    nmod_poly_init(cG, ctx->ffinfo->mod.n);
    nmod_poly_gcd_euclidean(cG, cA, cB);

    nmod_poly_init(cAbar, ctx->ffinfo->mod.n);
    nmod_poly_init(cBbar, ctx->ffinfo->mod.n);
    nmod_poly_div(cAbar, cA, cG);
    nmod_poly_div(cBbar, cB, cG);

    nmod_poly_init(gamma, ctx->ffinfo->mod.n);
    nmod_poly_gcd(gamma, nmod_mpolyun_leadcoeff_ref(A, ctx),
                         nmod_mpolyun_leadcoeff_ref(B, ctx));

    ldegA = nmod_mpolyun_lastdeg(A, ctx);
    ldegB = nmod_mpolyun_lastdeg(B, ctx);
    deggamma = nmod_poly_degree(gamma);
    bound = 1 + deggamma + FLINT_MAX(ldegA, ldegB);

    nmod_mpolyun_init(T, bits, ctx);
    nmod_poly_init(modulus, ctx->ffinfo->mod.n);
    nmod_poly_init(modulus2, ctx->ffinfo->mod.n);
    nmod_poly_one(modulus);

    nmod_mpolyun_init(Aeval, bits, ctx);
    nmod_mpolyun_init(Beval, bits, ctx);
    nmod_mpolyun_init(Geval, bits, ctx);
    nmod_mpolyun_init(Abareval, bits, ctx);
    nmod_mpolyun_init(Bbareval, bits, ctx);

    if ((ctx->ffinfo->mod.n & UWORD(1)) == UWORD(0))
    {
        success = 0;
        goto cleanup;
    }

    alpha = (ctx->ffinfo->mod.n - UWORD(1))/UWORD(2);

choose_prime: /* prime is v - alpha */

    if (alpha < 2)
    {
        success = 0;
        goto cleanup;
    }

    alpha--;

    /* make sure evaluation point does not kill both lc(A) and lc(B) */
    gammaeval = nmod_poly_evaluate_nmod(gamma, alpha);
    if (gammaeval == 0)
    {
        goto choose_prime;
    }

    /* evaluation point should kill neither A bor B */
    nmod_mpolyun_eval_last_un(Aeval, A, var, alpha, ctx);
    nmod_mpolyun_eval_last_un(Beval, B, var, alpha, ctx);
    FLINT_ASSERT(Aeval->length > 0);
    FLINT_ASSERT(Beval->length > 0);

    success = nmod_mpolyun_gcd_brown_smprime_ref(Geval, Abareval, Bbareval,
                                                   Aeval, Beval, var - 1, ctx);
    if (success == 0)
    {
        goto choose_prime;
    }

    FLINT_ASSERT(Geval->length > 0);
    FLINT_ASSERT(Abareval->length);
    FLINT_ASSERT(Bbareval->length > 0);

    if (nmod_mpolyun_is_nonzero_nmod(Geval, ctx))
    {
        nmod_mpolyun_one(G, ctx);
        nmod_mpolyun_swap(Abar, A);
        nmod_mpolyun_swap(Bbar, B);
        goto successful_put_content;    
    }

    if (nmod_poly_degree(modulus) > 0)
    {
        int cmp = 0;
        FLINT_ASSERT(G->length > 0);
        if (G->exps[0] != Geval->exps[0])
        {
            cmp = G->exps[0] > Geval->exps[0] ? 1 : -1;
        }
        if (cmp == 0)
        {
            slong k = nmod_poly_degree((Geval->coeffs + 0)->coeffs + 0);
            cmp = mpoly_monomial_cmp_nomask_extra(
                        (G->coeffs + 0)->exps + N*0,
                    (Geval->coeffs + 0)->exps + N*0, N, offset, k << shift);
        }

        if (cmp < 0)
        {
            goto choose_prime;
        }
        else if (cmp > 0)
        {
            nmod_poly_one(modulus);
        }
    }

    temp = nmod_mpolyn_leadcoeff_last(Geval->coeffs + 0, ctx);
    temp = n_invmod(temp, ctx->ffinfo->mod.n);
    temp = nmod_mul(gammaeval, temp, ctx->ffinfo->mod);
    nmod_mpolyun_scalar_mul_nmod(Geval, temp, ctx);

    if (nmod_poly_degree(modulus) > 0)
    {
        temp = nmod_poly_evaluate_nmod(modulus, alpha);
        temp = n_invmod(temp, ctx->ffinfo->mod.n);
        nmod_poly_scalar_mul_nmod(modulus, modulus, temp);
        nmod_mpolyun_addinterp_un(&ldegG, G, T, Geval, var, modulus, alpha, ctx);
        nmod_mpolyun_addinterp_un(&ldegAbar, Abar, T, Abareval, var, modulus, alpha, ctx);
        nmod_mpolyun_addinterp_un(&ldegBbar, Bbar, T, Bbareval, var, modulus, alpha, ctx);
    }
    else
    {
        nmod_mpolyun_set_popup(G, Geval, var, ctx);
        nmod_mpolyun_set_popup(Abar, Abareval, var, ctx);
        nmod_mpolyun_set_popup(Bbar, Bbareval, var, ctx);
        ldegG = 0;
        ldegAbar = 0;
        ldegBbar = 0;
    }

    nmod_poly_scalar_mul_nmod(modulus2, modulus, alpha);
    nmod_poly_shift_left(modulus, modulus, 1);
    nmod_poly_sub(modulus, modulus, modulus2);

    if (nmod_poly_degree(modulus) < bound)
    {
        goto choose_prime;
    }

    FLINT_ASSERT(ldegG >= 0);
    FLINT_ASSERT(ldegAbar >= 0);
    FLINT_ASSERT(ldegBbar >= 0);

    if (   deggamma + ldegA == ldegG + ldegAbar
        && deggamma + ldegB == ldegG + ldegBbar )
    {
        goto successful;
    }

    nmod_poly_one(modulus);
    goto choose_prime;

successful:

    nmod_mpolyun_content_last(modulus, G, ctx);
    nmod_mpolyun_divexact_last(G, modulus, ctx);
    nmod_mpolyun_divexact_last(Abar, nmod_mpolyun_leadcoeff_ref(G, ctx), ctx);
    nmod_mpolyun_divexact_last(Bbar, nmod_mpolyun_leadcoeff_ref(G, ctx), ctx);

successful_put_content:

    nmod_mpolyun_mul_last(G, cG, ctx);
    nmod_mpolyun_mul_last(Abar, cAbar, ctx);
    nmod_mpolyun_mul_last(Bbar, cBbar, ctx);

    success = 1;

cleanup:

#if WANT_ASSERT
    if (success)
    {
        FLINT_ASSERT(1 == nmod_mpolyun_leadcoeff_last(G, ctx));
        nmod_poly_mul(modulus, nmod_mpolyun_leadcoeff_ref(G, ctx),
                               nmod_mpolyun_leadcoeff_ref(Abar, ctx));
        FLINT_ASSERT(nmod_poly_equal(modulus, leadA));
        nmod_poly_mul(modulus, nmod_mpolyun_leadcoeff_ref(G, ctx),
                               nmod_mpolyun_leadcoeff_ref(Bbar, ctx));
        FLINT_ASSERT(nmod_poly_equal(modulus, leadB));
    }
    nmod_poly_clear(leadA);
    nmod_poly_clear(leadB);
#endif

    nmod_poly_clear(cA);
    nmod_poly_clear(cB);
    nmod_poly_clear(cG);
    nmod_poly_clear(cAbar);
    nmod_poly_clear(cBbar);

    nmod_poly_clear(gamma);

    nmod_mpolyun_clear(Aeval, ctx);
    nmod_mpolyun_clear(Beval, ctx);
    nmod_mpolyun_clear(Geval, ctx);
    nmod_mpolyun_clear(Abareval, ctx);
    nmod_mpolyun_clear(Bbareval, ctx);

    nmod_mpolyun_clear(T, ctx);

    nmod_poly_clear(modulus);
    nmod_poly_clear(modulus2);

    return success;
}




int nmod_mpolyun_gcd_brown_smprime_bivar(nmod_mpolyun_t G,
  nmod_mpolyun_t Abar, nmod_mpolyun_t Bbar, nmod_mpolyun_t A, nmod_mpolyun_t B,
                                                    const nmod_mpoly_ctx_t ctx)
{
    int success;
    slong bound;
    mp_limb_t alpha, temp, gammaevalp, gammaevalm;
    nmod_poly_t Aevalp, Bevalp, Gevalp, Abarevalp, Bbarevalp;
    nmod_poly_t Aevalm, Bevalm, Gevalm, Abarevalm, Bbarevalm;
    nmod_mpolyun_t T;
    slong deggamma, ldegG, ldegAbar, ldegBbar, ldegA, ldegB;
    nmod_poly_t cA, cB, cG, cAbar, cBbar, gamma;
    nmod_poly_t modulus, modulus2, alphapow, r;
    int gstab, astab, bstab, use_stab;
#if WANT_ASSERT
    nmod_poly_t leadA, leadB;
#endif

#if WANT_ASSERT
    nmod_poly_init(leadA, ctx->ffinfo->mod.n);
    nmod_poly_init(leadB, ctx->ffinfo->mod.n);
    nmod_poly_set(leadA, nmod_mpolyun_leadcoeff_ref(A, ctx));
    nmod_poly_set(leadB, nmod_mpolyun_leadcoeff_ref(B, ctx));
#endif

    nmod_poly_init(cA, ctx->ffinfo->mod.n);
    nmod_poly_init(cB, ctx->ffinfo->mod.n);
    nmod_mpolyun_content_last(cA, A, ctx);
    nmod_mpolyun_content_last(cB, B, ctx);
    nmod_mpolyun_divexact_last(A, cA, ctx);
    nmod_mpolyun_divexact_last(B, cB, ctx);

    nmod_poly_init(cG, ctx->ffinfo->mod.n);
    nmod_poly_gcd_euclidean(cG, cA, cB);

    nmod_poly_init(cAbar, ctx->ffinfo->mod.n);
    nmod_poly_init(cBbar, ctx->ffinfo->mod.n);
    nmod_poly_div(cAbar, cA, cG);
    nmod_poly_div(cBbar, cB, cG);

    nmod_poly_init(gamma, ctx->ffinfo->mod.n);
    nmod_poly_gcd(gamma, nmod_mpolyun_leadcoeff_ref(A, ctx),
                         nmod_mpolyun_leadcoeff_ref(B, ctx));

    ldegA = nmod_mpolyun_lastdeg(A, ctx);
    ldegB = nmod_mpolyun_lastdeg(B, ctx);
    deggamma = nmod_poly_degree(gamma);
    bound = 1 + deggamma + FLINT_MAX(ldegA, ldegB);

    nmod_poly_init(Aevalp, ctx->ffinfo->mod.n);
    nmod_poly_init(Bevalp, ctx->ffinfo->mod.n);
    nmod_poly_init(Gevalp, ctx->ffinfo->mod.n);
    nmod_poly_init(Abarevalp, ctx->ffinfo->mod.n);
    nmod_poly_init(Bbarevalp, ctx->ffinfo->mod.n);
    nmod_poly_init(Aevalm, ctx->ffinfo->mod.n);
    nmod_poly_init(Bevalm, ctx->ffinfo->mod.n);
    nmod_poly_init(Gevalm, ctx->ffinfo->mod.n);
    nmod_poly_init(Abarevalm, ctx->ffinfo->mod.n);
    nmod_poly_init(Bbarevalm, ctx->ffinfo->mod.n);

    nmod_mpolyun_init(T, A->bits, ctx);

    nmod_poly_init(r, ctx->ffinfo->mod.n);
    nmod_poly_init(alphapow, ctx->ffinfo->mod.n);
    nmod_poly_fit_length(alphapow, FLINT_MAX(WORD(3), bound + 1));
    nmod_poly_init(modulus, ctx->ffinfo->mod.n);
    nmod_poly_init(modulus2, ctx->ffinfo->mod.n);
    nmod_poly_one(modulus);

    if ((ctx->ffinfo->mod.n & UWORD(1)) == UWORD(0))
    {
        success = 0;
        goto cleanup;
    }

    use_stab = 1;
    gstab = bstab = astab = 0;

    alpha = (ctx->ffinfo->mod.n - UWORD(1))/UWORD(2);

choose_prime: /* primes are v - alpha, v + alpha */

    if (alpha < 2)
    {
        success = 0;
        goto cleanup;
    }

    alpha--;

    FLINT_ASSERT(0 < alpha && alpha <= ctx->ffinfo->mod.n/2);
    FLINT_ASSERT(alphapow->alloc >= 2);
    alphapow->length = 2;
    alphapow->coeffs[0] = 1;
    alphapow->coeffs[1] = alpha;

    /* make sure evaluation point does not kill both lc(A) and lc(B) */
    _nmod_poly_eval2_pow(&gammaevalp, &gammaevalm, gamma, alphapow, ctx->ffinfo);
    if (gammaevalp == 0 || gammaevalm == 0)
    {
        goto choose_prime;
    }

    /* evaluation point should kill neither A nor B */
    nmod_mpolyun_eval2_last_bivar(Aevalp, Aevalm, A, alphapow, ctx);
    nmod_mpolyun_eval2_last_bivar(Bevalp, Bevalm, B, alphapow, ctx);
    FLINT_ASSERT(Aevalp->length > 0);
    FLINT_ASSERT(Aevalm->length > 0);
    FLINT_ASSERT(Bevalp->length > 0);
    FLINT_ASSERT(Bevalm->length > 0);

    if (use_stab && gstab)
    {
        slong Gdeg;
        nmod_mpolyun_eval2_last_bivar(Gevalp, Gevalm, G, alphapow, ctx);
        Gdeg = G->exps[0];
        success = 1;
        success = success && nmod_poly_degree(Gevalp) == Gdeg;
        success = success && nmod_poly_degree(Gevalm) == Gdeg;
        success = success && Gevalp->coeffs[Gdeg] == gammaevalp;
        success = success && Gevalm->coeffs[Gdeg] == gammaevalm;
        nmod_poly_divrem_basecase(Abarevalp, r, Aevalp, Gevalp);
        success = success && (r->length == 0);
        nmod_poly_divrem_basecase(Abarevalm, r, Aevalm, Gevalm);
        success = success && (r->length == 0);
        nmod_poly_divrem_basecase(Bbarevalp, r, Bevalp, Gevalp);
        success = success && (r->length == 0);
        nmod_poly_divrem_basecase(Bbarevalm, r, Bevalm, Gevalm);
        success = success && (r->length == 0);

        if (!success)
        {
            use_stab = 0;
            nmod_poly_one(modulus);
            alpha = (ctx->ffinfo->mod.n - UWORD(1))/UWORD(2);
            goto choose_prime;
        }

        nmod_poly_scalar_mul_nmod(Abarevalp, Abarevalp, gammaevalp);
        nmod_poly_scalar_mul_nmod(Abarevalm, Abarevalm, gammaevalm);
        nmod_poly_scalar_mul_nmod(Bbarevalp, Bbarevalp, gammaevalp);
        nmod_poly_scalar_mul_nmod(Bbarevalm, Bbarevalm, gammaevalm);
    }
    else
    {
        nmod_poly_gcd(Gevalp, Aevalp, Bevalp);
        nmod_poly_div(Abarevalp, Aevalp, Gevalp);
        nmod_poly_div(Bbarevalp, Bevalp, Gevalp);
        nmod_poly_gcd(Gevalm, Aevalm, Bevalm);
        nmod_poly_div(Abarevalm, Aevalm, Gevalm);
        nmod_poly_div(Bbarevalm, Bevalm, Gevalm);
        gstab = astab = bstab = 0;
    }

    FLINT_ASSERT(Gevalp->length > 0);
    FLINT_ASSERT(Abarevalp->length > 0);
    FLINT_ASSERT(Bbarevalp->length > 0);
    FLINT_ASSERT(Gevalm->length > 0);
    FLINT_ASSERT(Abarevalm->length > 0);
    FLINT_ASSERT(Bbarevalm->length > 0);

    if (nmod_poly_degree(Gevalp) == 0 || nmod_poly_degree(Gevalm) == 0)
    {
        nmod_mpolyun_one(G, ctx);
        nmod_mpolyun_swap(Abar, A);
        nmod_mpolyun_swap(Bbar, B);
        goto successful_put_content;    
    }

    if (nmod_poly_degree(Gevalp) != nmod_poly_degree(Gevalm))
    {
        goto choose_prime;
    }

    /* the Geval have matching degrees */
    if (nmod_poly_degree(modulus) > 0)
    {
        FLINT_ASSERT(G->length > 0);
        if (nmod_poly_degree(Gevalp) > G->exps[0])
        {
            goto choose_prime;
        }
        else if (nmod_poly_degree(Gevalp) < G->exps[0])
        {
            nmod_poly_one(modulus);
        }
    }

    /* update interpolants */
    nmod_poly_scalar_mul_nmod(Gevalp, Gevalp, gammaevalp);
    nmod_poly_scalar_mul_nmod(Gevalm, Gevalm, gammaevalm);
    if (nmod_poly_degree(modulus) > 0)
    {
        temp = nmod_poly_evaluate_nmod(modulus, alpha);
        FLINT_ASSERT(temp == nmod_poly_evaluate_nmod(modulus, ctx->ffinfo->mod.n - alpha));
        temp = nmod_mul(temp, alpha, ctx->ffinfo->mod);
        temp = nmod_add(temp, temp, ctx->ffinfo->mod);
        nmod_poly_scalar_mul_nmod(modulus, modulus, n_invmod(temp, ctx->ffinfo->mod.n));
        if (!gstab)
        {
            gstab = !nmod_mpolyun_addinterp2_bivar(&ldegG, G, T, Gevalp, Gevalm, modulus, alphapow, ctx);
        }
        nmod_mpolyun_addinterp2_bivar(&ldegAbar, Abar, T, Abarevalp, Abarevalm, modulus, alphapow, ctx);
        nmod_mpolyun_addinterp2_bivar(&ldegBbar, Bbar, T, Bbarevalp, Bbarevalm, modulus, alphapow, ctx);
    }
    else
    {
        nmod_mpolyun_startinterp2_bivar(&ldegG, G, Gevalp, Gevalm, alpha, ctx);
        nmod_mpolyun_startinterp2_bivar(&ldegAbar, Abar, Abarevalp, Abarevalm, alpha, ctx);
        nmod_mpolyun_startinterp2_bivar(&ldegBbar, Bbar, Bbarevalp, Bbarevalm, alpha, ctx);
        gstab = astab = bstab = 0;
    }
    temp = nmod_mul(alpha, alpha, ctx->ffinfo->mod);
    nmod_poly_scalar_mul_nmod(modulus2, modulus, temp);
    nmod_poly_shift_left(modulus, modulus, 2);
    nmod_poly_sub(modulus, modulus, modulus2);

    if (nmod_poly_degree(modulus) < bound)
    {
        goto choose_prime;
    }

    FLINT_ASSERT(ldegG >= 0);
    FLINT_ASSERT(ldegAbar >= 0);
    FLINT_ASSERT(ldegBbar >= 0);

    if (   deggamma + ldegA == ldegG + ldegAbar
        && deggamma + ldegB == ldegG + ldegBbar )
    {
        goto successful;
    }

    nmod_poly_one(modulus);
    goto choose_prime;

successful:

    nmod_mpolyun_content_last(modulus, G, ctx);
    nmod_mpolyun_divexact_last(G, modulus, ctx);
    nmod_mpolyun_divexact_last(Abar, nmod_mpolyun_leadcoeff_ref(G, ctx), ctx);
    nmod_mpolyun_divexact_last(Bbar, nmod_mpolyun_leadcoeff_ref(G, ctx), ctx);

successful_put_content:

    nmod_mpolyun_mul_last(G, cG, ctx);
    nmod_mpolyun_mul_last(Abar, cAbar, ctx);
    nmod_mpolyun_mul_last(Bbar, cBbar, ctx);

    success = 1;

cleanup:

#if WANT_ASSERT
    if (success)
    {
        FLINT_ASSERT(1 == nmod_mpolyun_leadcoeff_last(G, ctx));
        nmod_poly_mul(modulus, nmod_mpolyun_leadcoeff_ref(G, ctx),
                               nmod_mpolyun_leadcoeff_ref(Abar, ctx));
        FLINT_ASSERT(nmod_poly_equal(modulus, leadA));
        nmod_poly_mul(modulus, nmod_mpolyun_leadcoeff_ref(G, ctx),
                               nmod_mpolyun_leadcoeff_ref(Bbar, ctx));
        FLINT_ASSERT(nmod_poly_equal(modulus, leadB));
    }
    nmod_poly_clear(leadA);
    nmod_poly_clear(leadB);
#endif

    nmod_poly_clear(cA);
    nmod_poly_clear(cB);
    nmod_poly_clear(cG);
    nmod_poly_clear(cAbar);
    nmod_poly_clear(cBbar);

    nmod_poly_clear(gamma);

    nmod_poly_clear(Aevalp);
    nmod_poly_clear(Bevalp);
    nmod_poly_clear(Gevalp);
    nmod_poly_clear(Abarevalp);
    nmod_poly_clear(Bbarevalp);
    nmod_poly_clear(Aevalm);
    nmod_poly_clear(Bevalm);
    nmod_poly_clear(Gevalm);
    nmod_poly_clear(Abarevalm);
    nmod_poly_clear(Bbarevalm);

    nmod_mpolyun_clear(T, ctx);

    nmod_poly_clear(r);
    nmod_poly_clear(alphapow);
    nmod_poly_clear(modulus);
    nmod_poly_clear(modulus2);

    return success;
}



int nmod_mpolyun_gcd_brown_smprime(nmod_mpolyun_t G,
                                nmod_mpolyun_t Abar, nmod_mpolyun_t Bbar,
                                    nmod_mpolyun_t A, nmod_mpolyun_t B,
                                         slong var, const nmod_mpoly_ctx_t ctx)
{
    int success;
    slong bound;
    slong offset, shift;
    mp_limb_t alpha, temp, gammaevalp, gammaevalm;
    nmod_mpolyun_t Aevalp, Bevalp, Gevalp, Abarevalp, Bbarevalp;
    nmod_mpolyun_t Aevalm, Bevalm, Gevalm, Abarevalm, Bbarevalm;
    nmod_mpolyun_t T;
    slong deggamma, ldegG, ldegAbar, ldegBbar, ldegA, ldegB;
    nmod_poly_t cA, cB, cG, cAbar, cBbar, gamma;
    nmod_poly_t modulus, modulus2, alphapow;
    mp_bitcnt_t bits = A->bits;
    slong N = mpoly_words_per_exp(bits, ctx->minfo);
#if WANT_ASSERT
    nmod_poly_t leadA, leadB;
#endif

    if (var == 0)
    {
        /* bivariate is more comfortable separated */
        return nmod_mpolyun_gcd_brown_smprime_bivar(G, Abar, Bbar, A, B, ctx);
    }

    mpoly_gen_offset_shift_sp(&offset, &shift, var - 1, G->bits, ctx->minfo);

#if WANT_ASSERT
    nmod_poly_init(leadA, ctx->ffinfo->mod.n);
    nmod_poly_init(leadB, ctx->ffinfo->mod.n);
    nmod_poly_set(leadA, nmod_mpolyun_leadcoeff_ref(A, ctx));
    nmod_poly_set(leadB, nmod_mpolyun_leadcoeff_ref(B, ctx));
#endif

    nmod_poly_init(cA, ctx->ffinfo->mod.n);
    nmod_poly_init(cB, ctx->ffinfo->mod.n);
    nmod_mpolyun_content_last(cA, A, ctx);
    nmod_mpolyun_content_last(cB, B, ctx);
    nmod_mpolyun_divexact_last(A, cA, ctx);
    nmod_mpolyun_divexact_last(B, cB, ctx);

    nmod_poly_init(cG, ctx->ffinfo->mod.n);
    nmod_poly_gcd_euclidean(cG, cA, cB);

    nmod_poly_init(cAbar, ctx->ffinfo->mod.n);
    nmod_poly_init(cBbar, ctx->ffinfo->mod.n);
    nmod_poly_div(cAbar, cA, cG);
    nmod_poly_div(cBbar, cB, cG);

    nmod_poly_init(gamma, ctx->ffinfo->mod.n);
    nmod_poly_gcd(gamma, nmod_mpolyun_leadcoeff_ref(A, ctx),
                         nmod_mpolyun_leadcoeff_ref(B, ctx));

    ldegA = nmod_mpolyun_lastdeg(A, ctx);
    ldegB = nmod_mpolyun_lastdeg(B, ctx);
    deggamma = nmod_poly_degree(gamma);
    bound = 1 + deggamma + FLINT_MAX(ldegA, ldegB);

    nmod_mpolyun_init(Aevalp, bits, ctx);
    nmod_mpolyun_init(Bevalp, bits, ctx);
    nmod_mpolyun_init(Gevalp, bits, ctx);
    nmod_mpolyun_init(Abarevalp, bits, ctx);
    nmod_mpolyun_init(Bbarevalp, bits, ctx);
    nmod_mpolyun_init(Aevalm, bits, ctx);
    nmod_mpolyun_init(Bevalm, bits, ctx);
    nmod_mpolyun_init(Gevalm, bits, ctx);
    nmod_mpolyun_init(Abarevalm, bits, ctx);
    nmod_mpolyun_init(Bbarevalm, bits, ctx);
    nmod_mpolyun_init(T, bits, ctx);

    nmod_poly_init(alphapow, ctx->ffinfo->mod.n);
    nmod_poly_fit_length(alphapow, FLINT_MAX(WORD(3), bound + 1));

    nmod_poly_init(modulus, ctx->ffinfo->mod.n);
    nmod_poly_init(modulus2, ctx->ffinfo->mod.n);
    nmod_poly_one(modulus);

    if ((ctx->ffinfo->mod.n & UWORD(1)) == UWORD(0))
    {
        success = 0;
        goto cleanup;
    }

    alpha = (ctx->ffinfo->mod.n - UWORD(1))/UWORD(2);

choose_prime:

    if (alpha < 2)
    {
        success = 0;
        goto cleanup;
    }

    alpha--;

    FLINT_ASSERT(0 < alpha && alpha <= ctx->ffinfo->mod.n/2);
    FLINT_ASSERT(alphapow->alloc >= 2);
    alphapow->length = 2;
    alphapow->coeffs[0] = 1;
    alphapow->coeffs[1] = alpha;

    /* make sure evaluation point does not kill both lc(A) and lc(B) */
    _nmod_poly_eval2_pow(&gammaevalp, &gammaevalm, gamma, alphapow, ctx->ffinfo);
    if (gammaevalp == 0 || gammaevalm == 0)
    {
        goto choose_prime;
    }

    /* evaluation point should kill neither A nor B */
    nmod_mpolyun_eval2_last_un(Aevalp, Aevalm, A, var, alphapow, ctx);
    nmod_mpolyun_eval2_last_un(Bevalp, Bevalm, B, var, alphapow, ctx);
    FLINT_ASSERT(Aevalp->length > 0);
    FLINT_ASSERT(Aevalm->length > 0);
    FLINT_ASSERT(Bevalp->length > 0);
    FLINT_ASSERT(Bevalm->length > 0);

    success = nmod_mpolyun_gcd_brown_smprime(Gevalp, Abarevalp, Bbarevalp,
                                               Aevalp, Bevalp, var - 1, ctx);
    success = success && nmod_mpolyun_gcd_brown_smprime(Gevalm, Abarevalm, Bbarevalm,
                                               Aevalm, Bevalm, var - 1, ctx);
    if (success == 0)
    {
        goto choose_prime;
    }

    FLINT_ASSERT(Gevalp->length > 0);
    FLINT_ASSERT(Abarevalp->length > 0);
    FLINT_ASSERT(Bbarevalp->length > 0);
    FLINT_ASSERT(Gevalm->length > 0);
    FLINT_ASSERT(Abarevalm->length > 0);
    FLINT_ASSERT(Bbarevalm->length > 0);

    if (   nmod_mpolyun_is_nonzero_nmod(Gevalp, ctx)
        || nmod_mpolyun_is_nonzero_nmod(Gevalm, ctx))
    {
        nmod_mpolyun_one(G, ctx);
        nmod_mpolyun_swap(Abar, A);
        nmod_mpolyun_swap(Bbar, B);
        goto successful_put_content;    
    }

    if (Gevalp->exps[0] != Gevalm->exps[0])
    {
        goto choose_prime;
    }
    if (   nmod_poly_degree((Gevalp->coeffs + 0)->coeffs + 0)
        != nmod_poly_degree((Gevalm->coeffs + 0)->coeffs + 0))
    {
        goto choose_prime;
    }
    if (!mpoly_monomial_equal((Gevalp->coeffs + 0)->exps + N*0,
                              (Gevalm->coeffs + 0)->exps + N*0, N))
    {
        goto choose_prime;
    }

    /* the Geval have matching degrees */
    if (nmod_poly_degree(modulus) > 0)
    {
        int cmp = 0;
        FLINT_ASSERT(G->length > 0);
        if (G->exps[0] != Gevalp->exps[0])
        {
            cmp = G->exps[0] > Gevalp->exps[0] ? 1 : -1;
        }
        if (cmp == 0)
        {
            slong k = nmod_poly_degree((Gevalp->coeffs + 0)->coeffs + 0);
            cmp = mpoly_monomial_cmp_nomask_extra(
                   (G->coeffs + 0)->exps + N*0,
                   (Gevalp->coeffs + 0)->exps + N*0, N, offset, k << shift);
        }

        if (cmp < 0)
        {
            goto choose_prime;
        }
        else if (cmp > 0)
        {
            nmod_poly_one(modulus);
        }
    }

    /* update interpolants */
    temp = nmod_mpolyn_leadcoeff_last(Gevalp->coeffs + 0, ctx);
    temp = n_invmod(temp, ctx->ffinfo->mod.n);
    temp = nmod_mul(gammaevalp, temp, ctx->ffinfo->mod);
    nmod_mpolyun_scalar_mul_nmod(Gevalp, temp, ctx);
    temp = nmod_mpolyn_leadcoeff_last(Gevalm->coeffs + 0, ctx);
    temp = n_invmod(temp, ctx->ffinfo->mod.n);
    temp = nmod_mul(gammaevalm, temp, ctx->ffinfo->mod);
    nmod_mpolyun_scalar_mul_nmod(Gevalm, temp, ctx);
    if (nmod_poly_degree(modulus) > 0)
    {
        temp = nmod_poly_evaluate_nmod(modulus, alpha);
        FLINT_ASSERT(temp == nmod_poly_evaluate_nmod(modulus, ctx->ffinfo->mod.n - alpha));
        temp = nmod_mul(temp, alpha, ctx->ffinfo->mod);
        temp = nmod_add(temp, temp, ctx->ffinfo->mod);
        nmod_poly_scalar_mul_nmod(modulus, modulus, n_invmod(temp, ctx->ffinfo->mod.n));
        nmod_mpolyun_addinterp2_un(&ldegG, G, T, Gevalp, Gevalm, var, modulus, alphapow, ctx);
        nmod_mpolyun_addinterp2_un(&ldegAbar, Abar, T, Abarevalp, Abarevalm, var, modulus, alphapow, ctx);
        nmod_mpolyun_addinterp2_un(&ldegBbar, Bbar, T, Bbarevalp, Bbarevalm, var, modulus, alphapow, ctx);
    }
    else
    {
        nmod_mpolyun_startinterp2_un(&ldegG, G, Gevalp, Gevalm, var, alpha, ctx);
        nmod_mpolyun_startinterp2_un(&ldegAbar, Abar, Abarevalp, Abarevalm, var, alpha, ctx);
        nmod_mpolyun_startinterp2_un(&ldegBbar, Bbar, Bbarevalp, Bbarevalm, var, alpha, ctx);
    }
    temp = nmod_mul(alpha, alpha, ctx->ffinfo->mod);
    nmod_poly_scalar_mul_nmod(modulus2, modulus, temp);
    nmod_poly_shift_left(modulus, modulus, 2);
    nmod_poly_sub(modulus, modulus, modulus2);

    if (nmod_poly_degree(modulus) < bound)
    {
        goto choose_prime;
    }

    FLINT_ASSERT(ldegG >= 0);
    FLINT_ASSERT(ldegAbar >= 0);
    FLINT_ASSERT(ldegBbar >= 0);

    if (   deggamma + ldegA == ldegG + ldegAbar
        && deggamma + ldegB == ldegG + ldegBbar )
    {
        goto successful;
    }

    nmod_poly_one(modulus);
    goto choose_prime;

successful:

    nmod_mpolyun_content_last(modulus, G, ctx);
    nmod_mpolyun_divexact_last(G, modulus, ctx);
    nmod_mpolyun_divexact_last(Abar, nmod_mpolyun_leadcoeff_ref(G, ctx), ctx);
    nmod_mpolyun_divexact_last(Bbar, nmod_mpolyun_leadcoeff_ref(G, ctx), ctx);

successful_put_content:

    nmod_mpolyun_mul_last(G, cG, ctx);
    nmod_mpolyun_mul_last(Abar, cAbar, ctx);
    nmod_mpolyun_mul_last(Bbar, cBbar, ctx);

    success = 1;

cleanup:

#if WANT_ASSERT
    if (success)
    {
        FLINT_ASSERT(1 == nmod_mpolyun_leadcoeff_last(G, ctx));
        nmod_poly_mul(modulus, nmod_mpolyun_leadcoeff_ref(G, ctx),
                               nmod_mpolyun_leadcoeff_ref(Abar, ctx));
        FLINT_ASSERT(nmod_poly_equal(modulus, leadA));
        nmod_poly_mul(modulus, nmod_mpolyun_leadcoeff_ref(G, ctx),
                               nmod_mpolyun_leadcoeff_ref(Bbar, ctx));
        FLINT_ASSERT(nmod_poly_equal(modulus, leadB));
    }
    nmod_poly_clear(leadA);
    nmod_poly_clear(leadB);
#endif

    nmod_poly_clear(cA);
    nmod_poly_clear(cB);
    nmod_poly_clear(cG);
    nmod_poly_clear(cAbar);
    nmod_poly_clear(cBbar);

    nmod_poly_clear(gamma);

    nmod_mpolyun_clear(Aevalp, ctx);
    nmod_mpolyun_clear(Bevalp, ctx);
    nmod_mpolyun_clear(Gevalp, ctx);
    nmod_mpolyun_clear(Abarevalp, ctx);
    nmod_mpolyun_clear(Bbarevalp, ctx);
    nmod_mpolyun_clear(Aevalm, ctx);
    nmod_mpolyun_clear(Bevalm, ctx);
    nmod_mpolyun_clear(Gevalm, ctx);
    nmod_mpolyun_clear(Abarevalm, ctx);
    nmod_mpolyun_clear(Bbarevalm, ctx);
    nmod_mpolyun_clear(T, ctx);

    nmod_poly_clear(alphapow);
    nmod_poly_clear(modulus);
    nmod_poly_clear(modulus2);

    return success;
}




/*****************************************************
    threaded
******************************************************/

typedef struct
{
    volatile int gcd_is_one;
    nmod_poly_struct * gamma;
    const nmod_mpoly_ctx_struct * ctx;
    nmod_mpolyun_struct * A, * B;
    ulong numthreads;
    slong var;
    slong bound;
}
_splitbase_struct;

typedef _splitbase_struct _splitbase_t[1];

typedef struct
{
    slong idx;
    _splitbase_struct * base;
    nmod_mpolyun_t G, Abar, Bbar;
    nmod_poly_t modulus;
    mp_limb_t alpha;
    slong required_images;
}
_splitworker_arg_struct;

static void _splitworker_bivar(void * varg)
{
    _splitworker_arg_struct * arg = (_splitworker_arg_struct *) varg;
    _splitbase_struct * base = arg->base;
    const nmod_mpoly_ctx_struct * ctx = base->ctx;
    nmod_poly_t Aeval, Beval, Geval, Abareval, Bbareval, modulus2;
    nmod_mpolyun_t T;
    mp_limb_t alpha, gammaeval, temp;
    ulong numthreads;
    slong ldeg;

    FLINT_ASSERT(base->var == 0);

    nmod_poly_init(Aeval, ctx->ffinfo->mod.n);
    nmod_poly_init(Beval, ctx->ffinfo->mod.n);
    nmod_poly_init(Geval, ctx->ffinfo->mod.n);
    nmod_poly_init(Abareval, ctx->ffinfo->mod.n);
    nmod_poly_init(Bbareval, ctx->ffinfo->mod.n);
    nmod_poly_init(modulus2, ctx->ffinfo->mod.n);
    nmod_mpolyun_init(T, base->A->bits, ctx);

    numthreads = base->numthreads;
    alpha = arg->alpha;

    nmod_poly_one(arg->modulus);
    while (nmod_poly_degree(arg->modulus) < arg->required_images)
    {
        /* get evaluation point */
        if (alpha <= numthreads)
        {
            break;
        }
        alpha -= numthreads;

        /* make sure evaluation point does not kill both lc(A) and lc(B) */
        gammaeval = nmod_poly_evaluate_nmod(base->gamma, alpha);
        if (gammaeval == 0)
        {
            continue;
        }

        /* make sure evaluation point does not kill either A or B */
        nmod_mpolyun_eval_last_bivar(Aeval, base->A, alpha, ctx);
        nmod_mpolyun_eval_last_bivar(Beval, base->B, alpha, ctx);
        if (Aeval->length == 0 || Beval->length == 0)
        {
            continue;
        }

        nmod_poly_gcd(Geval, Aeval, Beval);
        nmod_poly_div(Abareval, Aeval, Geval);
        nmod_poly_div(Bbareval, Beval, Geval);

        FLINT_ASSERT(Geval->length > 0);
        FLINT_ASSERT(Abareval->length > 0);
        FLINT_ASSERT(Bbareval->length > 0);

        /* check up */
        if (base->gcd_is_one)
        {
            break;
        }
        if (nmod_poly_degree(Geval) == 0)
        {
            base->gcd_is_one = 1;
            break;
        }
        if (nmod_poly_degree(arg->modulus) > 0)
        {
            FLINT_ASSERT(arg->G->length > 0);
            if (nmod_poly_degree(Geval) > arg->G->exps[0])
            {
                continue;
            }
            else if (nmod_poly_degree(Geval) < arg->G->exps[0])
            {
                nmod_poly_one(arg->modulus);
            }
        }

        /* update interpolants */
        nmod_poly_scalar_mul_nmod(Geval, Geval, gammaeval);
        if (nmod_poly_degree(arg->modulus) > 0)
        {
            temp = nmod_poly_evaluate_nmod(arg->modulus, alpha);
            temp = n_invmod(temp, ctx->ffinfo->mod.n);
            nmod_poly_scalar_mul_nmod(arg->modulus, arg->modulus, temp);
            nmod_mpolyun_addinterp_bivar(&ldeg, arg->G, T, Geval, arg->modulus, alpha, ctx);
            nmod_mpolyun_addinterp_bivar(&ldeg, arg->Abar, T, Abareval, arg->modulus, alpha, ctx);
            nmod_mpolyun_addinterp_bivar(&ldeg, arg->Bbar, T, Bbareval, arg->modulus, alpha, ctx);
        }
        else
        {
            nmod_mpolyun_set_poly(arg->G, Geval, ctx);
            nmod_mpolyun_set_poly(arg->Abar, Abareval, ctx);
            nmod_mpolyun_set_poly(arg->Bbar, Bbareval, ctx);
        }
        nmod_poly_scalar_mul_nmod(modulus2, arg->modulus, alpha);
        nmod_poly_shift_left(arg->modulus, arg->modulus, 1);
        nmod_poly_sub(arg->modulus, arg->modulus, modulus2);        
    }

    nmod_poly_clear(Aeval);
    nmod_poly_clear(Beval);
    nmod_poly_clear(Geval);
    nmod_poly_clear(Abareval);
    nmod_poly_clear(Bbareval);
    nmod_poly_clear(modulus2);
    nmod_mpolyun_clear(T, ctx);
}

static void _splitworker2_bivar(void * varg)
{
    _splitworker_arg_struct * arg = (_splitworker_arg_struct *) varg;
    _splitbase_struct * base = arg->base;
    const nmod_mpoly_ctx_struct * ctx = base->ctx;
    nmod_poly_t modulus2, alphapow, r;
    nmod_poly_t Aevalp, Bevalp, Gevalp, Abarevalp, Bbarevalp;
    nmod_poly_t Aevalm, Bevalm, Gevalm, Abarevalm, Bbarevalm;
    nmod_mpolyun_t T;
    mp_limb_t gammaevalp, alpha, temp;
    mp_limb_t gammaevalm;
    int gstab, astab, bstab, use_stab;
    ulong numthreads;
    slong ldeg;
    slong i;

    FLINT_ASSERT(base->var == 0);

    nmod_poly_init(r, ctx->ffinfo->mod.n);
    nmod_poly_init(modulus2, ctx->ffinfo->mod.n);
    nmod_poly_init(alphapow, ctx->ffinfo->mod.n);
    nmod_poly_fit_length(alphapow, base->bound + 1);

    nmod_poly_init(Aevalp, ctx->ffinfo->mod.n);
    nmod_poly_init(Bevalp, ctx->ffinfo->mod.n);
    nmod_poly_init(Gevalp, ctx->ffinfo->mod.n);
    nmod_poly_init(Abarevalp, ctx->ffinfo->mod.n);
    nmod_poly_init(Bbarevalp, ctx->ffinfo->mod.n);
    nmod_poly_init(Aevalm, ctx->ffinfo->mod.n);
    nmod_poly_init(Bevalm, ctx->ffinfo->mod.n);
    nmod_poly_init(Gevalm, ctx->ffinfo->mod.n);
    nmod_poly_init(Abarevalm, ctx->ffinfo->mod.n);
    nmod_poly_init(Bbarevalm, ctx->ffinfo->mod.n);
    nmod_mpolyun_init(T, base->A->bits, ctx);

    numthreads = base->numthreads;
    alpha = arg->alpha;

    use_stab = 1;
    gstab = bstab = astab = 0;

    nmod_poly_one(arg->modulus);
    while (nmod_poly_degree(arg->modulus) < arg->required_images)
    {
        /* get evaluation point */
        if (alpha <= numthreads)
        {
            break;
        }
        alpha -= numthreads;

        FLINT_ASSERT(0 < alpha && alpha <= ctx->ffinfo->mod.n/2);
        temp = 1;
        nmod_poly_set_coeff_ui(alphapow, 0, temp);
        for (i = 1; i <= base->bound; i++)
        {
            temp = nmod_mul(temp, alpha, ctx->ffinfo->mod);
            nmod_poly_set_coeff_ui(alphapow, i, temp);            
        }

        /* make sure evaluation point does not kill both lc(A) and lc(B) */
        _nmod_poly_eval2_pow(&gammaevalp, &gammaevalm, base->gamma, alphapow, ctx->ffinfo);
        if (gammaevalp == 0 || gammaevalm == 0)
        {
            continue;
        }

        /* make sure evaluation point does not kill either A or B */
        nmod_mpolyun_eval2_last_bivar(Aevalp, Aevalm, base->A, alphapow, ctx);
        nmod_mpolyun_eval2_last_bivar(Bevalp, Bevalm, base->B, alphapow, ctx);
        if (   Aevalp->length == 0 || Bevalp->length == 0
            || Aevalm->length == 0 || Bevalm->length == 0)
        {
            continue;
        }


        if (use_stab && gstab)
        {
            int success;
            slong Gdeg;

            nmod_mpolyun_eval2_last_bivar(Gevalp, Gevalm, arg->G, alphapow, ctx);
            Gdeg = arg->G->exps[0];
            success = 1;
            success = success && nmod_poly_degree(Gevalp) == Gdeg;
            success = success && nmod_poly_degree(Gevalm) == Gdeg;
            success = success && Gevalp->coeffs[Gdeg] == gammaevalp;
            success = success && Gevalm->coeffs[Gdeg] == gammaevalm;
            nmod_poly_divrem_basecase(Abarevalp, r, Aevalp, Gevalp);
            success = success && (r->length == 0);
            nmod_poly_divrem_basecase(Abarevalm, r, Aevalm, Gevalm);
            success = success && (r->length == 0);
            nmod_poly_divrem_basecase(Bbarevalp, r, Bevalp, Gevalp);
            success = success && (r->length == 0);
            nmod_poly_divrem_basecase(Bbarevalm, r, Bevalm, Gevalm);
            success = success && (r->length == 0);

            if (!success)
            {
                use_stab = 0;
                nmod_poly_one(arg->modulus);
                alpha = arg->alpha;
                continue;
            }

            nmod_poly_scalar_mul_nmod(Abarevalp, Abarevalp, gammaevalp);
            nmod_poly_scalar_mul_nmod(Abarevalm, Abarevalm, gammaevalm);
            nmod_poly_scalar_mul_nmod(Bbarevalp, Bbarevalp, gammaevalp);
            nmod_poly_scalar_mul_nmod(Bbarevalm, Bbarevalm, gammaevalm);
        }
        else
        {
            nmod_poly_gcd(Gevalp, Aevalp, Bevalp);
            nmod_poly_div(Abarevalp, Aevalp, Gevalp);
            nmod_poly_div(Bbarevalp, Bevalp, Gevalp);
            nmod_poly_gcd(Gevalm, Aevalm, Bevalm);
            nmod_poly_div(Abarevalm, Aevalm, Gevalm);
            nmod_poly_div(Bbarevalm, Bevalm, Gevalm);
        }

        FLINT_ASSERT(Gevalp->length > 0);
        FLINT_ASSERT(Abarevalp->length > 0);
        FLINT_ASSERT(Bbarevalp->length > 0);
        FLINT_ASSERT(Gevalm->length > 0);
        FLINT_ASSERT(Abarevalm->length > 0);
        FLINT_ASSERT(Bbarevalm->length > 0);

        /* check up */
        if (base->gcd_is_one)
        {
            break;
        }
        if (nmod_poly_degree(Gevalp) == 0 || nmod_poly_degree(Gevalm) == 0)
        {
            base->gcd_is_one = 1;
            break;
        }

        if (nmod_poly_degree(Gevalp) != nmod_poly_degree(Gevalm))
        {
            continue;
        }

        /* the Geval have matching degrees */
        if (nmod_poly_degree(arg->modulus) > 0)
        {
            FLINT_ASSERT(arg->G->length > 0);
            if (nmod_poly_degree(Gevalp) > arg->G->exps[0])
            {
                continue;
            }
            else if (nmod_poly_degree(Gevalp) < arg->G->exps[0])
            {
                nmod_poly_one(arg->modulus);
            }
        }
        /* update interpolants */
        nmod_poly_scalar_mul_nmod(Gevalp, Gevalp, gammaevalp);
        nmod_poly_scalar_mul_nmod(Gevalm, Gevalm, gammaevalm);
        if (nmod_poly_degree(arg->modulus) > 0)
        {
            temp = nmod_poly_evaluate_nmod(arg->modulus, alpha);
            FLINT_ASSERT(temp == nmod_poly_evaluate_nmod(arg->modulus, ctx->ffinfo->mod.n - alpha));
            temp = nmod_mul(temp, alpha, ctx->ffinfo->mod);
            temp = nmod_add(temp, temp, ctx->ffinfo->mod);
            nmod_poly_scalar_mul_nmod(arg->modulus, arg->modulus, n_invmod(temp, ctx->ffinfo->mod.n));
            if (!gstab)
            {
                gstab = !nmod_mpolyun_addinterp2_bivar(&ldeg, arg->G, T, Gevalp, Gevalm, arg->modulus, alphapow, ctx);
            }
            nmod_mpolyun_addinterp2_bivar(&ldeg, arg->Abar, T, Abarevalp, Abarevalm, arg->modulus, alphapow, ctx);
            nmod_mpolyun_addinterp2_bivar(&ldeg, arg->Bbar, T, Bbarevalp, Bbarevalm, arg->modulus, alphapow, ctx);
        }
        else
        {
            nmod_mpolyun_startinterp2_bivar(&ldeg, arg->G, Gevalp, Gevalm, alpha, ctx);
            nmod_mpolyun_startinterp2_bivar(&ldeg, arg->Abar, Abarevalp, Abarevalm, alpha, ctx);
            nmod_mpolyun_startinterp2_bivar(&ldeg, arg->Bbar, Bbarevalp, Bbarevalm, alpha, ctx);
            gstab = astab = bstab = 0;
        }
        temp = nmod_mul(alpha, alpha, ctx->ffinfo->mod);
        nmod_poly_scalar_mul_nmod(modulus2, arg->modulus, temp);
        nmod_poly_shift_left(arg->modulus, arg->modulus, 2);
        nmod_poly_sub(arg->modulus, arg->modulus, modulus2);
    }

    nmod_poly_clear(r);
    nmod_poly_clear(modulus2);
    nmod_poly_clear(alphapow);

    nmod_poly_clear(Aevalp);
    nmod_poly_clear(Bevalp);
    nmod_poly_clear(Gevalp);
    nmod_poly_clear(Abarevalp);
    nmod_poly_clear(Bbarevalp);
    nmod_poly_clear(Aevalm);
    nmod_poly_clear(Bevalm);
    nmod_poly_clear(Gevalm);
    nmod_poly_clear(Abarevalm);
    nmod_poly_clear(Bbarevalm);
    nmod_mpolyun_clear(T, ctx);
}


static void _splitworker(void * varg)
{
    _splitworker_arg_struct * arg = (_splitworker_arg_struct *) varg;
    _splitbase_struct * base = arg->base;
    const nmod_mpoly_ctx_struct * ctx = base->ctx;
    mp_bitcnt_t bits = base->A->bits;
    slong var = base->var;
    slong N = mpoly_words_per_exp(bits, ctx->minfo);
    slong offset, shift;
    nmod_poly_t modulus2;
    nmod_mpolyun_t Aeval, Beval, Geval, Abareval, Bbareval;
    nmod_mpolyun_t T;
    mp_limb_t alpha, gammaeval, temp;
    ulong numthreads;
    slong ldeg;
    int success;

    FLINT_ASSERT(var > 0);

    mpoly_gen_offset_shift_sp(&offset, &shift, var - 1, bits, ctx->minfo);

    nmod_mpolyun_init(Aeval, bits, ctx);
    nmod_mpolyun_init(Beval, bits, ctx);
    nmod_mpolyun_init(Geval, bits, ctx);
    nmod_mpolyun_init(Abareval, bits, ctx);
    nmod_mpolyun_init(Bbareval, bits, ctx);
    nmod_mpolyun_init(T, bits, ctx);
    nmod_poly_init(modulus2, ctx->ffinfo->mod.n);

    numthreads = base->numthreads;
    alpha = arg->alpha;

    nmod_poly_one(arg->modulus);
    while (nmod_poly_degree(arg->modulus) < arg->required_images)
    {
        /* get evaluation point */
        if (alpha <= numthreads)
        {
            break;
        }
        alpha -= numthreads;

        /* make sure evaluation point does not kill both lc(A) and lc(B) */
        gammaeval = nmod_poly_evaluate_nmod(base->gamma, alpha);
        if (gammaeval == 0)
        {
            continue;
        }

        /* make sure evaluation point does not kill either A or B */
        nmod_mpolyun_eval_last_un(Aeval, base->A, var, alpha, ctx);
        nmod_mpolyun_eval_last_un(Beval, base->B, var, alpha, ctx);
        if (Aeval->length == 0 || Beval->length == 0)
        {
            continue;
        }

        success = nmod_mpolyun_gcd_brown_smprime(Geval, Abareval, Bbareval,
                                                   Aeval, Beval, var - 1, ctx);
        if (success == 0)
        {
            continue;
        }

        FLINT_ASSERT(Geval->length > 0);
        FLINT_ASSERT(Abareval->length > 0);
        FLINT_ASSERT(Bbareval->length > 0);

        /* check up */
        if (base->gcd_is_one)
        {
            break;
        }

        if (nmod_mpolyun_is_nonzero_nmod(Geval, ctx))
        {
            base->gcd_is_one = 1;
            break;
        }

        if (nmod_poly_degree(arg->modulus) > 0)
        {
            int cmp = 0;
            FLINT_ASSERT(arg->G->length > 0);
            if (arg->G->exps[0] != Geval->exps[0])
            {
                cmp = arg->G->exps[0] > Geval->exps[0] ? 1 : -1;
            }
            if (cmp == 0)
            {
                slong k = nmod_poly_degree((Geval->coeffs + 0)->coeffs + 0);
                cmp = mpoly_monomial_cmp_nomask_extra(
                       (arg->G->coeffs + 0)->exps + N*0,
                        (Geval->coeffs + 0)->exps + N*0, N, offset, k << shift);
            }

            if (cmp < 0)
            {
                continue;
            }
            else if (cmp > 0)
            {
                nmod_poly_one(arg->modulus);
            }
        }

        /* update interpolants */
        temp = nmod_mpolyn_leadcoeff_last(Geval->coeffs + 0, ctx);
        temp = n_invmod(temp, ctx->ffinfo->mod.n);
        temp = nmod_mul(gammaeval, temp, ctx->ffinfo->mod);
        nmod_mpolyun_scalar_mul_nmod(Geval, temp, ctx);

        if (nmod_poly_degree(arg->modulus) > 0)
        {
            temp = nmod_poly_evaluate_nmod(arg->modulus, alpha);
            temp = n_invmod(temp, ctx->ffinfo->mod.n);
            nmod_poly_scalar_mul_nmod(arg->modulus, arg->modulus, temp);
            nmod_mpolyun_addinterp_un(&ldeg, arg->G, T, Geval, var, arg->modulus, alpha, ctx);
            nmod_mpolyun_addinterp_un(&ldeg, arg->Abar, T, Abareval, var, arg->modulus, alpha, ctx);
            nmod_mpolyun_addinterp_un(&ldeg, arg->Bbar, T, Bbareval, var, arg->modulus, alpha, ctx);
        }
        else
        {
            nmod_mpolyun_set_popup(arg->G, Geval, var, ctx);
            nmod_mpolyun_set_popup(arg->Abar, Abareval, var, ctx);
            nmod_mpolyun_set_popup(arg->Bbar, Bbareval, var, ctx);
        }

        nmod_poly_scalar_mul_nmod(modulus2, arg->modulus, alpha);
        nmod_poly_shift_left(arg->modulus, arg->modulus, 1);
        nmod_poly_sub(arg->modulus, arg->modulus, modulus2);
    }

    nmod_mpolyun_clear(Aeval, ctx);
    nmod_mpolyun_clear(Beval, ctx);
    nmod_mpolyun_clear(Geval, ctx);
    nmod_mpolyun_clear(Abareval, ctx);
    nmod_mpolyun_clear(Bbareval, ctx);
    nmod_mpolyun_clear(T, ctx);

    nmod_poly_clear(modulus2);
}

static void _splitworker2(void * varg)
{
    _splitworker_arg_struct * arg = (_splitworker_arg_struct *) varg;
    _splitbase_struct * base = arg->base;
    const nmod_mpoly_ctx_struct * ctx = base->ctx;
    mp_bitcnt_t bits = base->A->bits;
    slong var = base->var;
    slong N = mpoly_words_per_exp(bits, ctx->minfo);
    slong offset, shift;
    nmod_poly_t modulus2, alphapow;
    nmod_mpolyun_t Aevalp, Bevalp, Gevalp, Abarevalp, Bbarevalp;
    nmod_mpolyun_t Aevalm, Bevalm, Gevalm, Abarevalm, Bbarevalm;
    nmod_mpolyun_t T;
    mp_limb_t gammaevalp, alpha, temp;
    mp_limb_t gammaevalm;
    ulong numthreads;
    slong ldeg, i;
    int success;

    FLINT_ASSERT(var > 0);

    mpoly_gen_offset_shift_sp(&offset, &shift, var - 1, bits, ctx->minfo);

    nmod_poly_init(modulus2, ctx->ffinfo->mod.n);
    nmod_poly_init(alphapow, ctx->ffinfo->mod.n);
    nmod_poly_fit_length(alphapow, base->bound + 1);


    nmod_mpolyun_init(Aevalp, bits, ctx);
    nmod_mpolyun_init(Bevalp, bits, ctx);
    nmod_mpolyun_init(Gevalp, bits, ctx);
    nmod_mpolyun_init(Abarevalp, bits, ctx);
    nmod_mpolyun_init(Bbarevalp, bits, ctx);
    nmod_mpolyun_init(Aevalm, bits, ctx);
    nmod_mpolyun_init(Bevalm, bits, ctx);
    nmod_mpolyun_init(Gevalm, bits, ctx);
    nmod_mpolyun_init(Abarevalm, bits, ctx);
    nmod_mpolyun_init(Bbarevalm, bits, ctx);
    nmod_mpolyun_init(T, bits, ctx);

    numthreads = base->numthreads;
    alpha = arg->alpha;

    nmod_poly_one(arg->modulus);
    while (nmod_poly_degree(arg->modulus) < arg->required_images)
    {
        /* get evaluation point */
        if (alpha <= numthreads)
        {
            break;
        }
        alpha -= numthreads;

        FLINT_ASSERT(0 < alpha && alpha <= ctx->ffinfo->mod.n/2);
        temp = 1;
        nmod_poly_set_coeff_ui(alphapow, 0, temp);
        for (i = 1; i <= base->bound; i++)
        {
            temp = nmod_mul(temp, alpha, ctx->ffinfo->mod);
            nmod_poly_set_coeff_ui(alphapow, i, temp);            
        }

        /* make sure evaluation point does not kill both lc(A) and lc(B) */
        _nmod_poly_eval2_pow(&gammaevalp, &gammaevalm, base->gamma, alphapow, ctx->ffinfo);
        if (gammaevalp == 0 || gammaevalm == 0)
        {
            continue;
        }

        /* make sure evaluation point does not kill either A or B */
        nmod_mpolyun_eval2_last_un(Aevalp, Aevalm, base->A, var, alphapow, ctx);
        nmod_mpolyun_eval2_last_un(Bevalp, Bevalm, base->B, var, alphapow, ctx);
        if (   Aevalp->length == 0 || Bevalp->length == 0
            || Aevalm->length == 0 || Bevalm->length == 0)
        {
/*
printf("evaluation killed");
*/
            continue;
        }
/*
pthread_mutex_lock(&iomutex);
printf("Aevalp: "); nmod_mpolyun_print_pretty(Aevalp, NULL, ctx); printf("\n");
printf("Aevalm: "); nmod_mpolyun_print_pretty(Aevalm, NULL, ctx); printf("\n");
printf("Bevalp: "); nmod_mpolyun_print_pretty(Bevalp, NULL, ctx); printf("\n");
printf("Bevalm: "); nmod_mpolyun_print_pretty(Bevalm, NULL, ctx); printf("\n");
pthread_mutex_unlock(&iomutex);
*/


        success = nmod_mpolyun_gcd_brown_smprime(Gevalp, Abarevalp, Bbarevalp,
                                                   Aevalp, Bevalp, var - 1, ctx);
        success = success && nmod_mpolyun_gcd_brown_smprime(Gevalm, Abarevalm, Bbarevalm,
                                                   Aevalm, Bevalm, var - 1, ctx);
        if (success == 0)
        {
/*
printf("recursive gcd failed");
*/
            continue;
        }
/*
pthread_mutex_lock(&iomutex);
printf("Gevalp: "); nmod_mpolyun_print_pretty(Gevalp, NULL, ctx); printf("\n");
printf("Gevalm: "); nmod_mpolyun_print_pretty(Gevalm, NULL, ctx); printf("\n");
pthread_mutex_unlock(&iomutex);
*/
        FLINT_ASSERT(Gevalp->length > 0);
        FLINT_ASSERT(Abarevalp->length > 0);
        FLINT_ASSERT(Bbarevalp->length > 0);
        FLINT_ASSERT(Gevalm->length > 0);
        FLINT_ASSERT(Abarevalm->length > 0);
        FLINT_ASSERT(Bbarevalm->length > 0);

        /* check up */
        if (base->gcd_is_one)
        {
            break;
        }

        if (   nmod_mpolyun_is_nonzero_nmod(Gevalp, ctx)
            || nmod_mpolyun_is_nonzero_nmod(Gevalm, ctx))
        {
            base->gcd_is_one = 1;
            break;
        }

        if (Gevalp->exps[0] != Gevalm->exps[0])
        {
            continue;
        }
        if (   nmod_poly_degree((Gevalp->coeffs + 0)->coeffs + 0)
            != nmod_poly_degree((Gevalm->coeffs + 0)->coeffs + 0))
        {
            continue;
        }
        if (!mpoly_monomial_equal((Gevalp->coeffs + 0)->exps + N*0,
                                  (Gevalm->coeffs + 0)->exps + N*0, N))
        {
            continue;
        }

        /* the Geval have matching degrees */
        if (nmod_poly_degree(arg->modulus) > 0)
        {
            int cmp = 0;
            FLINT_ASSERT(arg->G->length > 0);
            if (arg->G->exps[0] != Gevalp->exps[0])
            {
                cmp = arg->G->exps[0] > Gevalp->exps[0] ? 1 : -1;
            }
            if (cmp == 0)
            {
                slong k = nmod_poly_degree((Gevalp->coeffs + 0)->coeffs + 0);
                cmp = mpoly_monomial_cmp_nomask_extra(
                       (arg->G->coeffs + 0)->exps + N*0,
                       (Gevalp->coeffs + 0)->exps + N*0, N, offset, k << shift);
            }

            if (cmp < 0)
            {
                continue;
            }
            else if (cmp > 0)
            {
                nmod_poly_one(arg->modulus);
            }
        }

        /* update interpolants */
        temp = nmod_mpolyn_leadcoeff_last(Gevalp->coeffs + 0, ctx);
        temp = n_invmod(temp, ctx->ffinfo->mod.n);
        temp = nmod_mul(gammaevalp, temp, ctx->ffinfo->mod);
        nmod_mpolyun_scalar_mul_nmod(Gevalp, temp, ctx);
        temp = nmod_mpolyn_leadcoeff_last(Gevalm->coeffs + 0, ctx);
        temp = n_invmod(temp, ctx->ffinfo->mod.n);
        temp = nmod_mul(gammaevalm, temp, ctx->ffinfo->mod);
        nmod_mpolyun_scalar_mul_nmod(Gevalm, temp, ctx);
        if (nmod_poly_degree(arg->modulus) > 0)
        {
            temp = nmod_poly_evaluate_nmod(arg->modulus, alpha);
            FLINT_ASSERT(temp == nmod_poly_evaluate_nmod(arg->modulus, ctx->ffinfo->mod.n - alpha));
            temp = nmod_mul(temp, alpha, ctx->ffinfo->mod);
            temp = nmod_add(temp, temp, ctx->ffinfo->mod);
            nmod_poly_scalar_mul_nmod(arg->modulus, arg->modulus, n_invmod(temp, ctx->ffinfo->mod.n));
            nmod_mpolyun_addinterp2_un(&ldeg, arg->G, T, Gevalp, Gevalm, var, arg->modulus, alphapow, ctx);
            nmod_mpolyun_addinterp2_un(&ldeg, arg->Abar, T, Abarevalp, Abarevalm, var, arg->modulus, alphapow, ctx);
            nmod_mpolyun_addinterp2_un(&ldeg, arg->Bbar, T, Bbarevalp, Bbarevalm, var, arg->modulus, alphapow, ctx);
        }
        else
        {
            nmod_mpolyun_startinterp2_un(&ldeg, arg->G, Gevalp, Gevalm, var, alpha, ctx);
            nmod_mpolyun_startinterp2_un(&ldeg, arg->Abar, Abarevalp, Abarevalm, var, alpha, ctx);
            nmod_mpolyun_startinterp2_un(&ldeg, arg->Bbar, Bbarevalp, Bbarevalm, var, alpha, ctx);
        }
        temp = nmod_mul(alpha, alpha, ctx->ffinfo->mod);
        nmod_poly_scalar_mul_nmod(modulus2, arg->modulus, temp);
        nmod_poly_shift_left(arg->modulus, arg->modulus, 2);
        nmod_poly_sub(arg->modulus, arg->modulus, modulus2);
    }

    nmod_poly_clear(modulus2);
    nmod_poly_clear(alphapow);

    nmod_mpolyun_clear(Aevalp, ctx);
    nmod_mpolyun_clear(Bevalp, ctx);
    nmod_mpolyun_clear(Gevalp, ctx);
    nmod_mpolyun_clear(Abarevalp, ctx);
    nmod_mpolyun_clear(Bbarevalp, ctx);
    nmod_mpolyun_clear(Aevalm, ctx);
    nmod_mpolyun_clear(Bevalm, ctx);
    nmod_mpolyun_clear(Gevalm, ctx);
    nmod_mpolyun_clear(Abarevalm, ctx);
    nmod_mpolyun_clear(Bbarevalm, ctx);
    nmod_mpolyun_clear(T, ctx);
}



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

static void _joinworker(void * varg)
{
    _joinworker_arg_struct * arg = (_joinworker_arg_struct *) varg;
    _joinbase_struct * base = arg->base;
    slong t, our_G_exp, our_Abar_exp, our_Bbar_exp;

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
            t = nmod_mpolyun_crt_exp(base->CRT, arg->G, our_G_exp,
                                 base->gptrs, base->numthreads, base->ctx);
            arg->G_lastdeg = FLINT_MAX(arg->G_lastdeg, t);
        }
        else if (our_Abar_exp >= 0)
        {
            t = nmod_mpolyun_crt_exp1(base->CRT, our_Abar_exp,
                              base->abarptrs, base->numthreads, base->ctx);
            arg->Abar_lastdeg = FLINT_MAX(arg->Abar_lastdeg, t);
        }
        else if (our_Bbar_exp >= 0)
        {
            t = nmod_mpolyun_crt_exp1(base->CRT, our_Bbar_exp,
                              base->bbarptrs, base->numthreads, base->ctx);
            arg->Bbar_lastdeg = FLINT_MAX(arg->Bbar_lastdeg, t);
        }
        else
        {
            return;
        }
    }
}



int nmod_mpolyun_gcd_brown_smprime_threaded(nmod_mpolyun_t G,
  nmod_mpolyun_t Abar, nmod_mpolyun_t Bbar, nmod_mpolyun_t A, nmod_mpolyun_t B,
                                         slong var, const nmod_mpoly_ctx_t ctx)
{
    slong i;
    mp_bitcnt_t bits = A->bits;
    slong N = mpoly_words_per_exp_sp(bits, ctx->minfo);
    ulong numthreads;
    int success;
    ulong bound;
    mp_limb_t alpha;
    slong deggamma, ldegGs, ldegAbars, ldegBbars, ldegA, ldegB;
    nmod_poly_t cA, cB, cG, cAbar, cBbar, gamma;
    nmod_poly_t cGs;
    slong Gexp0, Abarexp0, Bbarexp0;
    nmod_poly_struct ** mptrs;
    nmod_mpolyun_struct ** gptrs, ** abarptrs, ** bbarptrs;
    nmod_poly_crt_t P;
    _splitworker_arg_struct * splitargs;
    _splitbase_t splitbase;
    _joinworker_arg_struct * joinargs;
    _joinbase_t joinbase;
    slong max_numworkers, numworkers;
    thread_pool_handle * handles;
    slong * starts;
    slong Gi;
timeit_t time;

    FLINT_ASSERT(global_thread_pool_initialized);
    max_numworkers = thread_pool_get_size(global_thread_pool);
    if (max_numworkers == 0)
    {
        return nmod_mpolyun_gcd_brown_smprime(G, Abar, Bbar, A, B, var, ctx);
    }
/*
pthread_mutex_init(&iomutex, NULL);
*/


timeit_start(time);
    nmod_poly_init(cA, ctx->ffinfo->mod.n);
    nmod_poly_init(cB, ctx->ffinfo->mod.n);
    nmod_mpolyun_content_last(cA, A, ctx);
    nmod_mpolyun_content_last(cB, B, ctx);
    nmod_mpolyun_divexact_last(A, cA, ctx);
    nmod_mpolyun_divexact_last(B, cB, ctx);

    nmod_poly_init(cG, ctx->ffinfo->mod.n);
    nmod_poly_gcd_euclidean(cG, cA, cB);

    nmod_poly_init(cAbar, ctx->ffinfo->mod.n);
    nmod_poly_init(cBbar, ctx->ffinfo->mod.n);
    nmod_poly_div(cAbar, cA, cG);
    nmod_poly_div(cBbar, cB, cG);

    nmod_poly_init(gamma, ctx->ffinfo->mod.n);
    nmod_poly_gcd(gamma, nmod_mpolyun_leadcoeff_ref(A, ctx),
                         nmod_mpolyun_leadcoeff_ref(B, ctx));

    ldegA = nmod_mpolyun_lastdeg(A, ctx);
    ldegB = nmod_mpolyun_lastdeg(B, ctx);
    deggamma = nmod_poly_degree(gamma);
    bound = 1 + deggamma + FLINT_MAX(ldegA, ldegB);
/*
flint_printf("threaded var = %d\n", var);
flint_printf("A: "); nmod_mpolyun_print_pretty(A, NULL, ctx); printf("\n");
flint_printf("B: "); nmod_mpolyun_print_pretty(B, NULL, ctx); printf("\n");
flint_printf("bound: %wd\n", bound);
*/
/*
flint_printf("remove content: %wd\n", timeit_elapsed_wall(time));
*/
    if ((ctx->ffinfo->mod.n & UWORD(1)) == UWORD(0))
    {
        success = 0;
        goto cleanup;
    }

    alpha = (ctx->ffinfo->mod.n - UWORD(1))/UWORD(2);

    handles = (thread_pool_handle *) flint_malloc(max_numworkers *sizeof(thread_pool_handle));
    numworkers = thread_pool_request(global_thread_pool, handles, max_numworkers);
    numthreads = numworkers + 1;
    gptrs = (nmod_mpolyun_struct **) flint_malloc(numthreads*sizeof(nmod_mpolyun_struct *));
    abarptrs = (nmod_mpolyun_struct **) flint_malloc(numthreads*sizeof(nmod_mpolyun_struct *));
    bbarptrs = (nmod_mpolyun_struct **) flint_malloc(numthreads*sizeof(nmod_mpolyun_struct *));
    mptrs = (nmod_poly_struct **) flint_malloc(numthreads*sizeof(nmod_poly_struct *));
    splitargs = (_splitworker_arg_struct *) flint_malloc(numthreads*sizeof(_splitworker_arg_struct));
    for (i = 0; i < numthreads; i++)
    {
        nmod_mpolyun_init(splitargs[i].G, bits, ctx);
        nmod_mpolyun_init(splitargs[i].Abar, bits, ctx);
        nmod_mpolyun_init(splitargs[i].Bbar, bits, ctx);
        nmod_poly_init(splitargs[i].modulus, ctx->ffinfo->mod.n);
    }
    splitbase->numthreads = numthreads;
    splitbase->A = A;
    splitbase->B = B;
    splitbase->ctx = ctx;
    splitbase->gamma = gamma;
    splitbase->var = var;
    splitbase->bound = bound;

compute_split:

    if (alpha <= numthreads)
    {
        success = 0;
        goto cleanup_split;
    }
    alpha -= numthreads;

    splitbase->gcd_is_one = 0;
    for (i = 0; i < numthreads; i++)
    {
        slong ri = bound / numthreads + (i < (bound % numthreads));
        splitargs[i].idx = i;
        splitargs[i].base = splitbase;
        splitargs[i].alpha = alpha + i;
        splitargs[i].required_images = FLINT_MAX(WORD(1), ri);
    }

    for (i = 0; i + 1 < numthreads; i++)
    {
        thread_pool_wake(global_thread_pool, handles[i],
                  var == 0 ? _splitworker2_bivar : _splitworker2, &splitargs[i]);
    }
    (var == 0 ? _splitworker2_bivar : _splitworker2)(&splitargs[numthreads - 1]);
    for (i = 0; i + 1 < numthreads; i++)
    {
        thread_pool_wait(global_thread_pool, handles[i]);
    }

    if (splitbase->gcd_is_one)
    {
        success = 1;
        nmod_mpolyun_one(G, ctx);
        nmod_mpolyun_mul_last(G, cG, ctx);
        goto cleanup_split;
    }

    for (i = 0; i < numthreads; i++)
    {
        gptrs[i] = splitargs[i].G;
        abarptrs[i] = splitargs[i].Abar;
        bbarptrs[i] = splitargs[i].Bbar;
        mptrs[i] = splitargs[i].modulus;
        if (nmod_poly_degree(splitargs[i].modulus) < splitargs[i].required_images)
        {
            /* not enough evaluation points - must fail*/
            success = 0;
            goto cleanup_split;
        }
        FLINT_ASSERT(gptrs[i]->length > 0);
        FLINT_ASSERT(abarptrs[i]->length > 0);
        FLINT_ASSERT(bbarptrs[i]->length > 0);
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
            goto compute_split;
        }
        FLINT_ASSERT(splitargs[i].Abar->exps[0] == Abarexp0);
        FLINT_ASSERT(splitargs[i].Bbar->exps[0] == Bbarexp0);   
    }

/*
for (i = 0; i < numthreads; i++)
{
flint_printf("   m[%wd]: ", i); nmod_poly_print_pretty(mptrs[i], "v"); printf("\n");
flint_printf("   G[%wd]: ", i); nmod_mpolyun_print_pretty(gptrs[i], NULL, ctx); printf("\n");
flint_printf("Abar[%wd]: ", i); nmod_mpolyun_print_pretty(abarptrs[i], NULL, ctx); printf("\n");
flint_printf("Bbar[%wd]: ", i); nmod_mpolyun_print_pretty(bbarptrs[i], NULL, ctx); printf("\n");
}
*/

flint_printf("split: %wd\n", timeit_elapsed_wall(time));


    nmod_poly_crt_init(P);
    success = nmod_poly_crt_compile(P, mptrs, numthreads);
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
        joinargs[i].G_lastdeg = -WORD(1);
        joinargs[i].Abar_lastdeg = -WORD(1);
        joinargs[i].Bbar_lastdeg = -WORD(1);
        nmod_mpolyun_init(joinargs[i].G, bits, ctx);
        nmod_mpolyun_init(joinargs[i].Abar, bits, ctx);
        nmod_mpolyun_init(joinargs[i].Bbar, bits, ctx);
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

    ldegGs = ldegAbars = ldegBbars = -WORD(1);
    for (i = 0; i < numthreads; i++)
    {
        ldegGs = FLINT_MAX(ldegGs, joinargs[i].G_lastdeg);
        ldegAbars = FLINT_MAX(ldegAbars, joinargs[i].Abar_lastdeg);
        ldegBbars = FLINT_MAX(ldegBbars, joinargs[i].Bbar_lastdeg);
    }

/*
    flint_printf("ldegAbars: %wd", ldegAbars); printf("\n");
    flint_printf("ldegBbars: %wd", ldegBbars); printf("\n");
    flint_printf("ldegGs: %wd", ldegGs); printf("\n");

printf("joining G\n");
    for (i = 0; i < numthreads; i++)
    {
flint_printf("arg[%wd].G:", i); nmod_mpolyun_print_pretty(joinargs[i].G, NULL, ctx); printf("\n");
    }
*/

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
        nmod_mpolyun_fit_length(G, Gi + 1, ctx);
        G->exps[Gi] = max_exp;
        nmod_mpolyn_swap(G->coeffs + Gi, joinargs[max_pos].G->coeffs + starts[max_pos]);
        starts[max_pos]++;
        Gi++;
    }
    G->length = Gi;
    flint_free(starts);

    /* free join data */
    nmod_poly_crt_clear(P);
    for (i = 0; i < numthreads; i++)
    {
        nmod_mpolyun_clear(joinargs[i].G, ctx);
        nmod_mpolyun_clear(joinargs[i].Abar, ctx);
        nmod_mpolyun_clear(joinargs[i].Bbar, ctx);
    }
    flint_free(joinargs);

    if (   deggamma + ldegA != ldegGs + ldegAbars
        || deggamma + ldegB != ldegGs + ldegBbars
       )
    {
        /* divisibility test failed - could try again or just fail */
        goto compute_split;
    }

flint_printf("join: %wd\n", timeit_elapsed_wall(time));


    success = 1;

    nmod_poly_init(cGs, ctx->ffinfo->mod.n);
    nmod_mpolyun_content_last(cGs, G, ctx);
    nmod_mpolyun_divexact_last(G, cGs, ctx);
    nmod_mpolyun_mul_last(G, cG, ctx);
    nmod_poly_clear(cGs);

cleanup_split:

    for (i = 0; i < numthreads; i++)
    {
        nmod_mpolyun_clear(splitargs[i].G, ctx);
        nmod_mpolyun_clear(splitargs[i].Abar, ctx);
        nmod_mpolyun_clear(splitargs[i].Bbar, ctx);
        nmod_poly_clear(splitargs[i].modulus);
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

cleanup:

    nmod_poly_clear(cA);
    nmod_poly_clear(cB);
    nmod_poly_clear(cG);
    nmod_poly_clear(cAbar);
    nmod_poly_clear(cBbar);
    nmod_poly_clear(gamma);

flint_printf("finish: %wd\n", timeit_elapsed_wall(time));


/*
pthread_mutex_destroy(&iomutex);
*/
/*
if (!success)
{
printf("nmod_mpolyun_gcd_brown_smprime_threaded failed!!!\n");
}
*/
    return success;
}
