/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "threadpool.h"
#include "nmod_mpoly.h"
#include "fmpz_mpoly.h"
#include "profiler.h"

#define WORKER_COUNT 0

#define PROFILE_THIS TRUE

threadpool_t global_thread_pool;

const char * vars[] = {"x","y","z","t","u"};



#if defined(PROFILE_THIS)
typedef struct _vec_slong_struct
{
    slong * array;
    slong alloc;
    slong length;
} vec_slong_struct;

typedef vec_slong_struct vec_slong_t[1];

void vec_slong_init(vec_slong_t v)
{
    v->length = 0;
    v->alloc = 16;
    v->array = (slong *) flint_malloc(v->alloc*sizeof(slong));
}

void vec_slong_clear(vec_slong_t v)
{
    flint_free(v->array);
}

void vec_slong_push_back(vec_slong_t v, slong a)
{
    v->length++;
    if (v->length > v->alloc)
    {
        v->alloc = FLINT_MAX(v->length, 2*v->alloc);
        v->array = (slong *) flint_realloc(v->array, v->alloc*sizeof(slong));
    }
    v->array[v->length - 1] = a;
}

void vec_slong_print(const vec_slong_t v)
{
    slong i;
    flint_printf("[");
    for (i = 0; i < v->length; i++)
    {
        flint_printf("%wd",v->array[i]);
        if (i + 1 < v->length)
        {
            flint_printf(",",v->array[i]);
        }
    }
    flint_printf("]");
}
#endif


/*
    a thread safe mpoly supports three mutating operations
    - init from an array of terms
    - append an array of terms
    - clear out contents to a normal mpoly
*/
typedef struct _nmod_mpoly_ts_struct
{
    mp_limb_t * volatile coeffs; /* this is coeff_array[idx] */
    ulong * volatile exps;       /* this is exp_array[idx] */
    volatile slong length;
    slong alloc;
    mp_bitcnt_t bits;
    mp_bitcnt_t idx;    
    mp_limb_t * exp_array[FLINT_BITS];
    ulong * coeff_array[FLINT_BITS];
} nmod_mpoly_ts_struct;

typedef nmod_mpoly_ts_struct nmod_mpoly_ts_t[1];

void nmod_mpoly_ts_init(nmod_mpoly_ts_t A,
                              mp_limb_t * Bcoeff, ulong * Bexp, slong Blen,
                                                    mp_bitcnt_t bits, slong N)
{
    slong i;
    mp_bitcnt_t idx = FLINT_BIT_COUNT(Blen);
    idx = (idx <= 8) ? 0 : idx - 8;
    for (i = 0; i < FLINT_BITS; i++)
    {
        A->exp_array[i] = NULL;
        A->coeff_array[i] = NULL;
    }
    A->bits = bits;
    A->idx = idx;
    A->alloc = WORD(256) << idx;
    A->exps = A->exp_array[idx]
            = (ulong *) flint_malloc(N*A->alloc*sizeof(ulong));
    A->coeffs = A->coeff_array[idx]
              = (mp_limb_t *) flint_malloc(A->alloc*sizeof(mp_limb_t));
    A->length = Blen;
    for (i = 0; i < Blen; i++)
    {
        A->coeffs[i] = Bcoeff[i];
        mpoly_monomial_set(A->exps + N*i, Bexp + N*i, N);
    }
}

void nmod_mpoly_ts_print(const nmod_mpoly_ts_t B, const char ** x,
                                                    const nmod_mpoly_ctx_t ctx)
{
    nmod_mpoly_t A;
    A->length = B->length;
    A->alloc = B->alloc;
    A->coeffs = B->coeffs;
    A->exps = B->exps;
    A->bits = B->bits;
    nmod_mpoly_print_pretty(A, x, ctx);
}

void nmod_mpoly_ts_clear(nmod_mpoly_ts_t A)
{
    slong i;
    for (i = 0; i < FLINT_BITS; i++)
    {
        if (A->exp_array[i] != NULL)
        {
            FLINT_ASSERT(A->coeff_array[i] != NULL);
            flint_free(A->coeff_array[i]);    
            flint_free(A->exp_array[i]);
        }
    }
}

void nmod_mpoly_ts_clear_poly(nmod_mpoly_t Q, nmod_mpoly_ts_t A)
{
    if (Q->exps)
        flint_free(Q->exps);
    if (Q->coeffs)
        flint_free(Q->coeffs);

    Q->exps = A->exps;
    Q->coeffs = A->coeffs;
    Q->bits = A->bits;
    Q->alloc = A->alloc;
    Q->length = A->length;
    
    A->coeff_array[A->idx] = NULL;
    A->exp_array[A->idx] = NULL;
    nmod_mpoly_ts_clear(A);
}


/* put B on the end of A */
void nmod_mpoly_ts_append(nmod_mpoly_ts_t A,
                        mp_limb_t * Bcoeff, ulong * Bexps, slong Blen, slong N)
{
/* TODO: this needs barriers on non-x86 */

    slong i;
    ulong * oldexps = A->exps;
    ulong * oldcoeffs = A->coeffs;
    slong oldlength = A->length;
    slong newlength = A->length + Blen;

    /* we should be appending terms in order */
    FLINT_ASSERT(oldlength == 0 || mpoly_monomial_gt_nomask(
                                 oldexps + N*(oldlength - 1), Bexps + N*0, N));

    if (newlength <= A->alloc)
    {
        /* write new terms first */
        for (i = 0; i < Blen; i++)
        {
            oldcoeffs[oldlength + i] = Bcoeff[i];
            mpoly_monomial_set(oldexps + N*(oldlength + i), Bexps + N*i, N);
        }
    }
    else
    {
        slong newalloc;
        ulong * newexps;
        mp_limb_t * newcoeffs;
        mp_bitcnt_t newidx;
        newidx = FLINT_BIT_COUNT(newlength - 1);
        newidx = (newidx > 8) ? newidx - 8 : 0;
        FLINT_ASSERT(newidx > A->idx);

        newalloc = UWORD(256) << newidx;
        FLINT_ASSERT(newlength <= newalloc);
        newexps = A->exp_array[newidx]
                = (ulong *) flint_malloc(N*newalloc*sizeof(ulong));
        newcoeffs = A->coeff_array[newidx]
                  = (mp_limb_t *) flint_malloc(newalloc*sizeof(mp_limb_t));

        for (i = 0; i < oldlength; i++)
        {
            newcoeffs[i] = oldcoeffs[i];
            mpoly_monomial_set(newexps + N*i, oldexps + N*i, N);
        }
        for (i = 0; i < Blen; i++)
        {
            newcoeffs[oldlength + i] = Bcoeff[i];
            mpoly_monomial_set(newexps + N*(oldlength + i), Bexps + N*i, N);
        }

        A->alloc = newalloc;
        A->exps = newexps;
        A->coeffs = newcoeffs;
        A->idx = newidx;

        /* do not free oldcoeff/exps as other threads may be using them */
    }

    /* update length at the very end */
    A->length = newlength;
}

typedef struct _Lchunk_struct
{
    nmod_mpoly_t polyC;
    struct _Lchunk_struct * next;
    ulong * emin;
    ulong * emax;
    slong startidx;
    slong endidx;
    int upperclosed;
    volatile int lock;
    volatile int producer;
    volatile slong ma;
    volatile slong mq;
    int Cinited;
#if defined(PROFILE_THIS)
    slong idx;
#endif
} Lchunk_struct;

typedef Lchunk_struct Lchunk_t[1];



typedef struct
{
    pthread_mutex_t mutex;
    Lchunk_struct * head;
    Lchunk_struct * tail;
    Lchunk_struct * volatile cur;
    nmod_mpoly_t polyA;
    nmod_mpoly_t polyB;
    nmod_mpoly_ts_t polyQ;
    const nmod_mpoly_ctx_struct * ctx;
    slong length;
    slong N;
    mp_bitcnt_t bits;
    ulong * cmpmask;
    int failed;
#if defined(PROFILE_THIS)
    timeit_t timer;
#endif

} Lholder_struct;

typedef Lholder_struct Lholder_t[1];


static void Lholder_init(Lholder_t H)
{
    H->head = NULL;
    H->tail = NULL;
    H->cur = NULL;
    H->ctx = NULL;
    H->length = 0;
    H->N = 0;
    H->bits = 0;
    H->cmpmask = NULL;
}

static void Lchunk_clear(Lchunk_t L, Lholder_t H)
{
    if (L->Cinited)
    {
        nmod_mpoly_clear(L->polyC, H->ctx);
    }
}
static int Lholder_clear(nmod_mpoly_t Q, Lholder_t H)
{
    Lchunk_struct * L = H->head;
    while (L != NULL)
    {
        Lchunk_struct * nextL = L->next;
        Lchunk_clear(L, H);
        flint_free(L);
        L = nextL;
    }
    H->head = NULL;
    H->tail = NULL;
    H->cur = NULL;
    H->ctx = NULL;
    H->length = 0;
    H->N = 0;
    H->bits = 0;
    H->cmpmask = NULL;

    if (H->failed)
    {
        nmod_mpoly_zero(Q, H->ctx);
        nmod_mpoly_ts_clear(H->polyQ);
        return 0;
    }
    else
    {
        nmod_mpoly_ts_clear_poly(Q, H->polyQ);
        return 1;
    }
}

static void Lholder_add_chunk(Lholder_t H, Lchunk_t L)
{
    L->next = NULL;

    if (H->tail == NULL)
    {
        FLINT_ASSERT(H->head == NULL);
        H->tail = L;
        H->head = L;
    }
    else
    {
        Lchunk_struct * tail = H->tail;
        FLINT_ASSERT(tail->next == NULL);
        tail->next = L;
        H->tail = L;
    }
    H->length++;
}

static void Lchunk_print(const Lchunk_t L)
{
    flint_printf("********** L entry *********\n");
    flint_printf("emax: %016llx\n",L->emax[0]);
    flint_printf("emin: %016llx\n",L->emin[0]);
    flint_printf("ma: %wd\n", L->ma);
    flint_printf("mq: %wd\n", L->mq);
}

static void Lholder_print(const Lholder_t H)
{
    const Lchunk_struct * us = H->head;
    flint_printf("LHolder (%wd entries):\n", H->length);
    while (us != NULL)
    {
        Lchunk_print(us);
        us = us->next;
    }
}


int select_exps(fmpz_mpoly_t S, fmpz_mpoly_ctx_t zctx,
          ulong * Aexp, slong Alen, ulong * Bexp, slong Blen, mp_bitcnt_t bits)
{
    int failure;
    ulong mask;
    ulong * Sexp;
    slong Slen;
    fmpz * Scoeff;
    slong ns = 20;
    slong Astep = FLINT_MAX(Alen*2/ns, WORD(1));
    slong Bstep = FLINT_MAX(Blen*4/ns, WORD(1));
    slong tot;
    ulong * T0, * T1;
    slong i, N;
    TMP_INIT;

    TMP_START;

    N = mpoly_words_per_exp(bits, zctx->minfo);

    mask = 0; /* mask will be unused if bits > FLINT_BITS*/
    for (i = 0; i < FLINT_BITS/bits; i++)
        mask = (mask << bits) + (UWORD(1) << (bits - 1));

    FLINT_ASSERT(Alen > 0);
    FLINT_ASSERT(Blen > 0);

    tot = 10;
    tot += (Alen + Astep - 1)/Astep;
    tot += 2*((Blen + Bstep - 1)/Bstep);
    fmpz_mpoly_fit_bits(S, bits, zctx);
    S->bits = bits;
    fmpz_mpoly_fit_length(S, tot, zctx);
    Sexp = S->exps;
    Scoeff = S->coeffs;

    Slen = 0;
    for (i = 0; i < Alen; i += Astep)
    {
        mpoly_monomial_set(Sexp + N*Slen, Aexp + N*i, N);
        fmpz_one(Scoeff + Slen);
        Slen++;   
    }
    _fmpz_mpoly_set_length(S, Slen, zctx);

    T0 = (ulong *) TMP_ALLOC(N*sizeof(ulong));
    T1 = (ulong *) TMP_ALLOC(N*sizeof(ulong));
    mpoly_monomial_sub_mp(T0, Aexp + N*0, Bexp + N*0, N);
    mpoly_monomial_sub_mp(T1, Aexp + N*(Alen - 1), Bexp + N*(Blen - 1), N);
    if (bits <= FLINT_BITS)
    {
        if (   mpoly_monomial_overflows(T0, N, mask)
            || mpoly_monomial_overflows(T1, N, mask))
        {
            failure = 1;
            goto cleanup;
        }
    }
    else
    {
        if (   mpoly_monomial_overflows_mp(T0, N, bits)
            || mpoly_monomial_overflows_mp(T1, N, bits))
        {
            failure = 1;
            goto cleanup;
        }
    }

    for (i = 0; i < Blen; i += Bstep)
    {
        mpoly_monomial_sub_mp(Sexp + N*Slen, Aexp + N*0, Bexp + N*0, N);
        mpoly_monomial_add_mp(Sexp + N*Slen, Sexp + N*Slen, Bexp + N*i, N);
        fmpz_one(Scoeff + Slen);
        if (bits <= FLINT_BITS)
            Slen += !(mpoly_monomial_overflows(Sexp + N*Slen, N, mask));
        else
            Slen += !(mpoly_monomial_overflows_mp(Sexp + N*Slen, N, bits));

        mpoly_monomial_sub_mp(Sexp + N*Slen, Aexp + N*(Alen - 1), Bexp + N*(Blen - 1), N);
        mpoly_monomial_add_mp(Sexp + N*Slen, Sexp + N*Slen, Bexp + N*i, N);
        fmpz_one(Scoeff + Slen);
        if (bits <= FLINT_BITS)
            Slen += !(mpoly_monomial_overflows(Sexp + N*Slen, N, mask));
        else
            Slen += !(mpoly_monomial_overflows_mp(Sexp + N*Slen, N, bits));
    }

    mpoly_monomial_zero(Sexp + N*Slen, N);
    fmpz_one(Scoeff + Slen);
    Slen++;

    FLINT_ASSERT(Slen < tot);

    _fmpz_mpoly_set_length(S, Slen, zctx);

    fmpz_mpoly_sort_terms(S, zctx);
    fmpz_mpoly_combine_like_terms(S, zctx);

    failure = 0;

cleanup:
    TMP_END;
    return failure;
}


void nmod_mpoly_chop(nmod_mpoly_t C, const nmod_mpoly_t A,
                int upperclosed, ulong * emax, ulong * emin, const Lholder_t H)
{
    slong i, N;
    mp_limb_t * Ccoeff;
    ulong * Cexp;
    slong Calloc;
    slong Clen;

    FLINT_ASSERT(A->bits == H->bits);
    nmod_mpoly_fit_length(C, 16, H->ctx);
    nmod_mpoly_fit_bits(C, H->bits, H->ctx);
    C->bits = H->bits;
    Ccoeff = C->coeffs;
    Cexp = C->exps;
    Calloc =  C->alloc;
    Clen = 0;

    N = H->N;
    i = 0;
    while (i < A->length)
    {
        if (mpoly_monomial_gt(A->exps + N*i, emax, N, H->cmpmask))
        {
            break;
        }
        if (upperclosed && mpoly_monomial_equal(A->exps + N*i, emax, N))
        {
            break;
        }
        i++;
    }

    while (i < A->length)
    {
        if (mpoly_monomial_gt(A->exps + N*i, emin, N, H->cmpmask))
        {
            break;
        }
        _nmod_mpoly_fit_length(&Ccoeff, &Cexp, &Calloc, Clen + 1, N);
        Ccoeff[Clen] = A->coeffs[i];
        mpoly_monomial_set(Cexp + N*Clen, A->exps + N*i, N);
        Clen++;
        i++;
    }

    C->coeffs = Ccoeff;
    C->exps = Cexp;
    C->alloc = Calloc;
    C->length = Clen;
}

void nmod_mpoly_chopmulmonomial(nmod_mpoly_t T, const nmod_mpoly_t B,
           const nmod_mpoly_t Q, slong k, const Lchunk_t L, const Lholder_t H)
{
    slong i, N = H->N;
    mp_limb_t * Tcoeff = T->coeffs;
    ulong * Texp = T->exps;
    slong Talloc = T->alloc;
    slong Tlen;
    ulong * emin = L->emin;
    ulong * emax = L->emax;
    int upperclosed = L->upperclosed;

    FLINT_ASSERT(T->bits == H->bits);

    Tlen = 0;
    for (i = 0; i < B->length; i++)
    {
        _nmod_mpoly_fit_length(&Tcoeff, &Texp, &Talloc, Tlen + 1, N);
        mpoly_monomial_add_mp(Texp + N*Tlen, B->exps + N*i, Q->exps + N*k, N);
        if (mpoly_monomial_cmp(Texp + N*Tlen, emax, N, H->cmpmask) < upperclosed
            && mpoly_monomial_cmp(Texp + N*Tlen, emin, N, H->cmpmask) >= 0)
        {
            Tcoeff[Tlen] = nmod_mul(B->coeffs[i], Q->coeffs[k], H->ctx->ffinfo->mod);
            Tlen++;
        }
    }

    T->coeffs = Tcoeff;
    T->exps = Texp;
    T->alloc = Talloc;
    T->length = Tlen;
/*
flint_printf("B: "); nmod_mpoly_print_pretty(B, vars, H->ctx); flint_printf("\n");
flint_printf("Q: "); nmod_mpoly_print_pretty(Q, vars, H->ctx); flint_printf("\n");
flint_printf("Q index: %wd\n", k);
flint_printf("T: "); nmod_mpoly_print_pretty(T, vars, H->ctx); flint_printf("\n");
*/
}


int mpoly_monomial_cmp1(ulong a, ulong b, ulong cmpmask)
{
    if ((a^cmpmask) != (b^cmpmask))
    {
        if ((a^cmpmask) > (b^cmpmask))
            return 1;
        else
            return -1;
    }
    return 0;
}

int mpoly_monomial_gt1(ulong a, ulong b, ulong cmpmask)
{
    return (a^cmpmask) > (b^cmpmask);
}

int mpoly_monomial_ge1(ulong a, ulong b, ulong cmpmask)
{
    return (a^cmpmask) >= (b^cmpmask);
}


typedef struct _nmod_mpoly_stripe_struct
{
    char * big_mem;
    slong big_mem_alloc;
    const nmod_mpoly_ctx_struct * ctx;
    slong N;
    mp_bitcnt_t bits;
    ulong * cmpmask;
    slong * startidx;
    slong * endidx;
    ulong * emin;
    ulong * emax;
    int upperclosed;
} nmod_mpoly_stripe_struct;

typedef nmod_mpoly_stripe_struct nmod_mpoly_stripe_t[1];


typedef struct _worker_arg_struct
{
    Lholder_struct * H;
    nmod_mpoly_stripe_t S;
    nmod_mpoly_t polyT1;
    nmod_mpoly_t polyT2;
#if defined(PROFILE_THIS)
    vec_slong_t time_data;
#endif
} worker_arg_struct;

typedef worker_arg_struct worker_arg_t[1];


/*
A = B * C
*/
slong _nmod_mpoly_mul_stripe1(mp_limb_t ** Acoeff_, ulong ** Aexp_, slong * Aalloc_,
                 const mp_limb_t * Bcoeff, const ulong * Bexp, slong Blen,
                 const mp_limb_t * Ccoeff, const ulong * Cexp, slong Clen,
                                                   const nmod_mpoly_stripe_t S)
{
    int upperclosed;
    slong startidx, endidx;
    ulong prev_startidx;
    slong N = S->N;
    ulong maskhi = S->cmpmask[0];
    ulong emax = S->emax[0];
    ulong emin = S->emin[0];
    slong i, j;
    slong next_loc = Blen + 4;   /* something bigger than heap can ever be */
    slong Q_len = 0, heap_len = 1; /* heap zero index unused */
    mpoly_heap1_s * heap;
    mpoly_heap_t * chain;
    slong * Q;
    mpoly_heap_t * x;
    slong Alen;
    slong Aalloc = *Aalloc_;
    mp_limb_t * Acoeff = *Acoeff_;
    ulong * Aexp = *Aexp_;
    ulong acc0, acc1, acc2, pp0, pp1;
    ulong exp;
    slong * ends;
    ulong texp;
    slong * hind;
    slong toff;
    char * tbase = S->big_mem;

    FLINT_ASSERT(N == 1);

    toff = 0;
    hind = (slong *)(tbase + toff); toff += Blen*sizeof(slong);
    ends = (slong *)(tbase + toff); toff += Blen*sizeof(slong);
    Q = (slong *) (tbase + toff); toff += 2*Blen*sizeof(slong);
    heap = (mpoly_heap1_s *)(tbase + toff); toff += (Blen + 1)*sizeof(mpoly_heap1_s);
    chain = (mpoly_heap_t *)(tbase + toff); toff += Blen*sizeof(mpoly_heap_t);
    FLINT_ASSERT(toff <= S->big_mem_alloc);

/*
    hind = (slong *) TMP_ALLOC(Blen*sizeof(slong));
    ends = (slong *) TMP_ALLOC(Blen*sizeof(slong));
    heap = (mpoly_heap1_s *) TMP_ALLOC((Blen + 1)*sizeof(mpoly_heap1_s));
    chain = (mpoly_heap_t *) TMP_ALLOC(Blen*sizeof(mpoly_heap_t));
    Q = (slong *) TMP_ALLOC(2*Blen*sizeof(slong));
*/

    startidx = *S->startidx;
    endidx = *S->endidx;
    upperclosed = S->upperclosed;
    emax = S->emax[0];
    emin = S->emin[0];

    /* put all the starting nodes on the heap */
    prev_startidx = -UWORD(1);
    for (i = 0; i < Blen; i++)
    {
        if (startidx < Clen)
        {
            texp = Bexp[i] + Cexp[startidx];
            FLINT_ASSERT(mpoly_monomial_cmp1(emax, texp, maskhi) > -upperclosed);
        }
        while (startidx > 0)
        {
            texp = Bexp[i] + Cexp[startidx - 1];
            if (mpoly_monomial_cmp1(emax, texp, maskhi) <= -upperclosed)
            {
                break;
            }
            startidx--;
        }

        if (endidx < Clen)
        {
            texp = Bexp[i] + Cexp[endidx];
            FLINT_ASSERT(mpoly_monomial_cmp1(emin, texp, maskhi) > 0);
        }
        while (endidx > 0)
        {
            texp = Bexp[i] + Cexp[endidx - 1];
            if (mpoly_monomial_cmp1(emin, texp, maskhi) <= 0)
            {
                break;
            }
            endidx--;
        }

        ends[i] = endidx;

        hind[i] = 2*startidx + 1;

        if (  (startidx < endidx)
           && (((ulong)startidx) < prev_startidx) 
           )
        {
            x = chain + i;
            x->i = i;
            x->j = startidx;
            x->next = NULL;
            hind[x->i] = 2*(x->j + 1) + 0;
            _mpoly_heap_insert1(heap, Bexp[x->i] + Cexp[x->j], x,
                                      &next_loc, &heap_len, maskhi);
        }

        prev_startidx = startidx;
    }

    *S->startidx = startidx;
    *S->endidx = endidx;

    Alen = 0;
    while (heap_len > 1)
    {
        exp = heap[1].exp;

        _nmod_mpoly_fit_length(&Acoeff, &Aexp, &Aalloc, Alen + 1, N);

        Aexp[Alen] = exp;

        acc0 = acc1 = acc2 = 0;
        do
        {
            x = _mpoly_heap_pop1(heap, &heap_len, maskhi);

            hind[x->i] |= WORD(1);
            Q[Q_len++] = x->i;
            Q[Q_len++] = x->j;
            umul_ppmm(pp1, pp0, Bcoeff[x->i], Ccoeff[x->j]);
            add_sssaaaaaa(acc2, acc1, acc0, acc2, acc1, acc0, WORD(0), pp1, pp0);

            while ((x = x->next) != NULL)
            {
                hind[x->i] |= WORD(1);
                Q[Q_len++] = x->i;
                Q[Q_len++] = x->j;
                umul_ppmm(pp1, pp0, Bcoeff[x->i], Ccoeff[x->j]);
                add_sssaaaaaa(acc2, acc1, acc0, acc2, acc1, acc0, WORD(0), pp1, pp0);
            }
        } while (heap_len > 1 && heap[1].exp == exp);

        NMOD_RED3(Acoeff[Alen], acc2, acc1, acc0, S->ctx->ffinfo->mod);
        Alen += (Acoeff[Alen] != 0);

        /* for each node temporarily stored */
        while (Q_len > 0)
        {
            /* take node from store */
            j = Q[--Q_len];
            i = Q[--Q_len];

            /* should we go right? */
            if (  (i + 1 < Blen)
               && (j + 0 < ends[i + 1])
               && (hind[i + 1] == 2*j + 1)
               )
            {
                x = chain + i + 1;
                x->i = i + 1;
                x->j = j;
                x->next = NULL;

                hind[x->i] = 2*(x->j + 1) + 0;
                _mpoly_heap_insert1(heap, Bexp[x->i] + Cexp[x->j], x,
                                                 &next_loc, &heap_len, maskhi);
            }

            /* should we go up? */
            if (  (j + 1 < ends[i + 0])
               && ((hind[i] & 1) == 1)
               && (  (i == 0)
                  || (hind[i - 1] >= 2*(j + 2) + 1)
                  )
               )
            {
                x = chain + i;
                x->i = i;
                x->j = j + 1;
                x->next = NULL;

                hind[x->i] = 2*(x->j + 1) + 0;
                _mpoly_heap_insert1(heap, Bexp[x->i] + Cexp[x->j], x,
                                                 &next_loc, &heap_len, maskhi);
            }
        }
    }

    *Acoeff_ = Acoeff;
    *Aexp_ = Aexp;
    *Aalloc_ = Aalloc;

    return Alen;
}




slong _nmod_mpoly_mul_stripe(mp_limb_t ** Acoeff_, ulong ** Aexp_, slong * Aalloc_,
                 const mp_limb_t * Bcoeff, const ulong * Bexp, slong Blen,
                 const mp_limb_t * Ccoeff, const ulong * Cexp, slong Clen,
                                                   const nmod_mpoly_stripe_t S)
{
    int upperclosed;
    slong startidx, endidx;
    ulong prev_startidx;
    ulong * emax = S->emax;
    ulong * emin = S->emin;
    slong N = S->N;
    slong i, j;
    slong next_loc = Blen + 4;   /* something bigger than heap can ever be */
    slong Q_len = 0, heap_len = 1; /* heap zero index unused */
    mpoly_heap_s * heap;
    mpoly_heap_t * chain;
    slong * Q;
    mpoly_heap_t * x;
    slong Alen;
    slong Aalloc = *Aalloc_;
    mp_limb_t * Acoeff = *Acoeff_;
    ulong * Aexp = *Aexp_;
    ulong acc0, acc1, acc2, pp0, pp1;
    ulong * exp, * exps;
    ulong ** exp_list;
    slong exp_next;
    slong * ends;
    ulong * texp;
    slong * hind;
    slong toff;
    char * tbase = S->big_mem;

/*
    hind = (slong *) TMP_ALLOC(Blen*sizeof(slong));
    ends = (slong *) TMP_ALLOC(Blen*sizeof(slong));
    texp = (ulong *) TMP_ALLOC(N*sizeof(ulong));
    heap = (mpoly_heap_s *) TMP_ALLOC((Blen + 1)*sizeof(mpoly_heap_s));
    chain = (mpoly_heap_t *) TMP_ALLOC(Blen*sizeof(mpoly_heap_t));
    Q = (slong *) TMP_ALLOC(2*Blen*sizeof(slong));
    exps = (ulong *) TMP_ALLOC(Blen*N*sizeof(ulong));
    exp_list = (ulong **) TMP_ALLOC(Blen*sizeof(ulong *));
*/

    toff = 0;
    hind = (slong *)(tbase + toff);         toff += Blen*sizeof(slong);
    ends = (slong *)(tbase + toff);         toff += Blen*sizeof(slong);
    Q = (slong *) (tbase + toff);           toff += 2*Blen*sizeof(slong);
    heap = (mpoly_heap_s *)(tbase + toff);  toff += (Blen + 1)*sizeof(mpoly_heap_s);
    chain = (mpoly_heap_t *)(tbase + toff); toff += Blen*sizeof(mpoly_heap_t);
    exps = (ulong *)(tbase + toff);         toff +=  Blen*N*sizeof(ulong);
    exp_list = (ulong **)(tbase + toff);    toff +=  Blen*sizeof(ulong *);
    texp = (ulong *)(tbase + toff);         toff +=  N*sizeof(ulong);
    FLINT_ASSERT(toff <= S->big_mem_alloc);

    exp_next = 0;

    startidx = *S->startidx;
    endidx = *S->endidx;
    upperclosed = S->upperclosed;

    for (i = 0; i < Blen; i++)
        exp_list[i] = exps + i*N;

    /* put all the starting nodes on the heap */
    prev_startidx = -UWORD(1);
    for (i = 0; i < Blen; i++)
    {
        if (startidx < Clen)
        {
            mpoly_monomial_add_mp(texp, Bexp + N*i, Cexp + N*startidx, N);
            FLINT_ASSERT(mpoly_monomial_cmp(emax, texp, N, S->cmpmask) > -upperclosed);
        }
        while (startidx > 0)
        {
            mpoly_monomial_add_mp(texp, Bexp + N*i, Cexp + N*(startidx - 1), N);
            if (mpoly_monomial_cmp(emax, texp, N, S->cmpmask) <= -upperclosed)
            {
                break;
            }
            startidx--;
        }

        if (endidx < Clen)
        {
            mpoly_monomial_add_mp(texp, Bexp + N*i, Cexp + N*endidx, N);
            FLINT_ASSERT(mpoly_monomial_cmp(emin, texp, N, S->cmpmask) > 0);
        }
        while (endidx > 0)
        {
            mpoly_monomial_add_mp(texp, Bexp + N*i, Cexp + N*(endidx - 1), N);
            if (mpoly_monomial_cmp(emin, texp, N, S->cmpmask) <= 0)
            {
                break;
            }
            endidx--;
        }

        ends[i] = endidx;

        hind[i] = 2*startidx + 1;

        if (  (startidx < endidx)
           && (((ulong)startidx) < prev_startidx) 
           )
        {
            x = chain + i;
            x->i = i;
            x->j = startidx;
            x->next = NULL;
            hind[x->i] = 2*(x->j + 1) + 0;

            mpoly_monomial_add(exp_list[exp_next], Bexp + N*x->i, Cexp + N*x->j, N);

            if (!_mpoly_heap_insert(heap, exp_list[exp_next++], x,
                                      &next_loc, &heap_len, N, S->cmpmask))
               exp_next--;
        }

        prev_startidx = startidx;
    }

    *S->startidx = startidx;
    *S->endidx = endidx;

    Alen = 0;
    while (heap_len > 1)
    {
        exp = heap[1].exp;

        _nmod_mpoly_fit_length(&Acoeff, &Aexp, &Aalloc, Alen + 1, N);

        mpoly_monomial_set(Aexp + N*Alen, exp, N);

        acc0 = acc1 = acc2 = 0;
        do
        {
            exp_list[--exp_next] = heap[1].exp;

            x = _mpoly_heap_pop(heap, &heap_len, N, S->cmpmask);

            hind[x->i] |= WORD(1);
            Q[Q_len++] = x->i;
            Q[Q_len++] = x->j;
            umul_ppmm(pp1, pp0, Bcoeff[x->i], Ccoeff[x->j]);
            add_sssaaaaaa(acc2, acc1, acc0, acc2, acc1, acc0, WORD(0), pp1, pp0);

            while ((x = x->next) != NULL)
            {
                hind[x->i] |= WORD(1);
                Q[Q_len++] = x->i;
                Q[Q_len++] = x->j;
                umul_ppmm(pp1, pp0, Bcoeff[x->i], Ccoeff[x->j]);
                add_sssaaaaaa(acc2, acc1, acc0, acc2, acc1, acc0, WORD(0), pp1, pp0);
            }
        } while (heap_len > 1 && mpoly_monomial_equal(heap[1].exp, exp, N));

        NMOD_RED3(Acoeff[Alen], acc2, acc1, acc0, S->ctx->ffinfo->mod);
        Alen += (Acoeff[Alen] != 0);

        /* for each node temporarily stored */
        while (Q_len > 0)
        {
            /* take node from store */
            j = Q[--Q_len];
            i = Q[--Q_len];

            /* should we go right? */
            if (  (i + 1 < Blen)
               && (j + 0 < ends[i + 1])
               && (hind[i + 1] == 2*j + 1)
               )
            {
                x = chain + i + 1;
                x->i = i + 1;
                x->j = j;
                x->next = NULL;

                hind[x->i] = 2*(x->j + 1) + 0;

                mpoly_monomial_add_mp(exp_list[exp_next], Bexp + N*x->i, Cexp + N*x->j, N);

                if (!_mpoly_heap_insert(heap, exp_list[exp_next++], x,
                                      &next_loc, &heap_len, N, S->cmpmask))
                    exp_next--;
            }

            /* should we go up? */
            if (  (j + 1 < ends[i + 0])
               && ((hind[i] & 1) == 1)
               && (  (i == 0)
                  || (hind[i - 1] >= 2*(j + 2) + 1)
                  )
               )
            {
                x = chain + i;
                x->i = i;
                x->j = j + 1;
                x->next = NULL;

                hind[x->i] = 2*(x->j + 1) + 0;

                mpoly_monomial_add_mp(exp_list[exp_next], Bexp + N*x->i, Cexp + N*x->j, N);

                if (!_mpoly_heap_insert(heap, exp_list[exp_next++], x,
                                      &next_loc, &heap_len, N, S->cmpmask))
                    exp_next--;
            }
        }
    }

    *Acoeff_ = Acoeff;
    *Aexp_ = Aexp;
    *Aalloc_ = Aalloc;

    return Alen;
}


slong _nmod_mpoly_divides_stripe1(
                     mp_limb_t ** Qcoeff_,      ulong ** Qexp_, slong * Qalloc_,
                const mp_limb_t * Acoeff, const ulong * Aexp, slong Alen,
                const mp_limb_t * Bcoeff, const ulong * Bexp, slong Blen,
                                                   const nmod_mpoly_stripe_t S)
{
    mp_bitcnt_t bits = S->bits;
    ulong emin = S->emin[0];
    ulong cmpmask = S->cmpmask[0];
    ulong texp;
    int lt_divides;
    slong i, j, s;
    slong next_loc, heap_len;
    mpoly_heap1_s * heap;
    mpoly_heap_t * chain;
    slong * store, * store_base;
    mpoly_heap_t * x;

    slong Qlen;
    slong Qalloc = * Qalloc_;
    mp_limb_t * Qcoeff = * Qcoeff_;
    ulong * Qexp = * Qexp_;

    ulong exp;
    mp_limb_t lc_minus_inv, acc0, acc1, acc2, pp1, pp0;
    ulong mask;
    slong * hind;
    slong toff;
    char * tbase = S->big_mem;

    FLINT_ASSERT(Alen > 0);
    FLINT_ASSERT(Blen > 0);
    FLINT_ASSERT(S->N == 1);

    next_loc = Blen + 4;   /* something bigger than heap can ever be */
/*
    hind = (slong *) TMP_ALLOC(Blen*sizeof(slong));
    store = store_base = (slong *) TMP_ALLOC(2*Blen*sizeof(mpoly_heap_t *));
    heap = (mpoly_heap_s *) TMP_ALLOC((Blen + 1)*sizeof(mpoly_heap_s));
    chain = (mpoly_heap_t *) TMP_ALLOC(Blen*sizeof(mpoly_heap_t));
    exps = (ulong *) TMP_ALLOC(Blen*N*sizeof(ulong));
    exp_list = (ulong **) TMP_ALLOC(Blen*sizeof(ulong *));
*/

    toff = 0;
    hind = (slong *)(tbase + toff);         toff += Blen*sizeof(slong);
    store = store_base = (slong *) (tbase + toff);           toff += 2*Blen*sizeof(slong);
    heap = (mpoly_heap1_s *)(tbase + toff);  toff += (Blen + 1)*sizeof(mpoly_heap1_s);
    chain = (mpoly_heap_t *)(tbase + toff); toff += Blen*sizeof(mpoly_heap_t);
    FLINT_ASSERT(toff <= S->big_mem_alloc);

    for (i = 0; i < Blen; i++)
        hind[i] = 1;

    /* mask with high bit set in each word of each field of exponent vector */
    mask = 0;
    for (i = 0; i < FLINT_BITS/bits; i++)
        mask = (mask << bits) + (UWORD(1) << (bits - 1));

    Qlen = WORD(0);

    /* s is the number of terms * (latest quotient) we should put into heap */
    s = Blen;

    /* insert (-1, 0, exp2[0]) into heap */
    heap_len = 2;
    x = chain + 0;
    x->i = -WORD(1);
    x->j = 0;
    x->next = NULL;
    heap[1].next = x;
    HEAP_ASSIGN(heap[1], Aexp[0], x);

    FLINT_ASSERT(mpoly_monomial_cmp1(Aexp[0], emin, cmpmask) >= 0);

    /* precompute leading cofficient info */
    lc_minus_inv = S->ctx->ffinfo->mod.n - nmod_inv(Bcoeff[0], S->ctx->ffinfo->mod);

    while (heap_len > 1)
    {
        exp = heap[1].exp;

        if (mpoly_monomial_overflows1(exp, mask))
            goto not_exact_division;

        FLINT_ASSERT(mpoly_monomial_cmp1(exp, emin, cmpmask) >= 0);

        _nmod_mpoly_fit_length(&Qcoeff, &Qexp, &Qalloc, Qlen + 1, 1);

        lt_divides = mpoly_monomial_divides1(Qexp + Qlen, exp, Bexp[0], mask);

        acc0 = acc1 = acc2 = 0;
        do
        {
            x = _mpoly_heap_pop1(heap, &heap_len, cmpmask);
            do
            {
                *store++ = x->i;
                *store++ = x->j;
                if (x->i != -WORD(1))
                    hind[x->i] |= WORD(1);

                if (x->i == -WORD(1))
                {
                    add_sssaaaaaa(acc2, acc1, acc0, acc2, acc1, acc0,
                                 WORD(0), WORD(0), S->ctx->ffinfo->mod.n - Acoeff[x->j]);
                } else
                {
                    umul_ppmm(pp1, pp0, Bcoeff[x->i], Qcoeff[x->j]);
                    add_sssaaaaaa(acc2, acc1, acc0, acc2, acc1, acc0, WORD(0), pp1, pp0);                    
                }
            } while ((x = x->next) != NULL);
        } while (heap_len > 1 && heap[1].exp == exp);

        NMOD_RED3(Qcoeff[Qlen], acc2, acc1, acc0, S->ctx->ffinfo->mod);

        /* process nodes taken from the heap */
        while (store > store_base)
        {
            j = *--store;
            i = *--store;

            if (i == -WORD(1))
            {
                /* take next dividend term */
                if (j + 1 < Alen)
                {
                    x = chain + 0;
                    x->i = i;
                    x->j = j + 1;
                    x->next = NULL;

                    FLINT_ASSERT(mpoly_monomial_cmp1(Aexp[x->j], emin, cmpmask) >= 0);

                    _mpoly_heap_insert1(heap, Aexp[x->j], x,
                                                &next_loc, &heap_len, cmpmask);
                }
            } else
            {
                /* should we go up */
                if (  (i + 1 < Blen)
                   && (hind[i + 1] == 2*j + 1)
                   )
                {
                    x = chain + i + 1;
                    x->i = i + 1;
                    x->j = j;
                    x->next = NULL;
                    hind[x->i] = 2*(x->j + 1) + 0;

                    texp = Bexp[x->i] + Qexp[x->j];
                    if (mpoly_monomial_cmp1(texp, emin, cmpmask) >= 0)
                    {
                        _mpoly_heap_insert1(heap, texp, x,
                                                &next_loc, &heap_len, cmpmask);
                    }
                    else
                    {
                        hind[x->i] |= 1;
                    }
                }
                /* should we go up? */
                if (j + 1 == Qlen)
                {
                    s++;
                } else if (  ((hind[i] & 1) == 1)
                          && ((i == 1) || (hind[i - 1] >= 2*(j + 2) + 1))
                          )
                {
                    x = chain + i;
                    x->i = i;
                    x->j = j + 1;
                    x->next = NULL;
                    hind[x->i] = 2*(x->j + 1) + 0;

                    texp = Bexp[x->i] + Qexp[x->j];
                    if (mpoly_monomial_cmp1(texp, emin, cmpmask) >= 0)
                    {
                        _mpoly_heap_insert1(heap, texp, x,
                                                &next_loc, &heap_len, cmpmask);
                    }
                    else
                    {
                        hind[x->i] |= 1;
                    }
                }
            }
        }

        Qcoeff[Qlen] = nmod_mul(Qcoeff[Qlen], lc_minus_inv, S->ctx->ffinfo->mod);
        if (Qcoeff[Qlen] == 0)
        {
            continue;
        }

        if (!lt_divides)
        {
            goto not_exact_division;
        }

        if (s > 1)
        {
            i = 1;
            x = chain + i;
            x->i = i;
            x->j = Qlen;
            x->next = NULL;
            hind[x->i] = 2*(x->j + 1) + 0;
/*
            mpoly_monomial_add_mp(exp_list[exp_next], Bexp + N*x->i, Qexp + N*x->j, N);

            if (mpoly_monomial_cmp(exp_list[exp_next], S->emin, N, S->cmpmask) >= 0)
            {

                if (!_mpoly_heap_insert(heap, exp_list[exp_next++], x,
                                      &next_loc, &heap_len, N, S->cmpmask))
                    exp_next--;
            }
            else
            {
                hind[x->i] |= 1;                
            }
*/
            texp = Bexp[x->i] + Qexp[x->j];
            if (mpoly_monomial_cmp1(texp, emin, cmpmask) >= 0)
            {
                _mpoly_heap_insert1(heap, texp, x,
                                         &next_loc, &heap_len, cmpmask);
            }
            else
            {
                hind[x->i] |= 1;
            }
        }
        s = 1;      
        Qlen++;
    }


cleanup:

    *Qalloc_ = Qalloc;
    *Qcoeff_ = Qcoeff;
    *Qexp_ = Qexp;

    return Qlen;

not_exact_division:
    Qlen = 0;
    goto cleanup;
}

slong _nmod_mpoly_divides_stripe(
                     mp_limb_t ** Qcoeff_,      ulong ** Qexp_, slong * Qalloc_,
                const mp_limb_t * Acoeff, const ulong * Aexp, slong Alen,
                const mp_limb_t * Bcoeff, const ulong * Bexp, slong Blen,
                                                   const nmod_mpoly_stripe_t S)
{
    mp_bitcnt_t bits = S->bits;
    slong N = S->N;
    int lt_divides;
    slong i, j, s;
    slong next_loc, heap_len;
    mpoly_heap_s * heap;
    mpoly_heap_t * chain;
    slong * store, * store_base;
    mpoly_heap_t * x;

    slong Qlen;
    slong Qalloc = * Qalloc_;
    mp_limb_t * Qcoeff = * Qcoeff_;
    ulong * Qexp = * Qexp_;

    ulong * exp, * exps;
    ulong ** exp_list;
    slong exp_next;
    mp_limb_t lc_minus_inv, acc0, acc1, acc2, pp1, pp0;
    ulong mask;
    slong * hind;
    slong toff;
    char * tbase = S->big_mem;

    FLINT_ASSERT(Alen > 0);
    FLINT_ASSERT(Blen > 0);

    next_loc = Blen + 4;   /* something bigger than heap can ever be */
/*
    hind = (slong *) TMP_ALLOC(Blen*sizeof(slong));
    store = store_base = (slong *) TMP_ALLOC(2*Blen*sizeof(mpoly_heap_t *));
    heap = (mpoly_heap_s *) TMP_ALLOC((Blen + 1)*sizeof(mpoly_heap_s));
    chain = (mpoly_heap_t *) TMP_ALLOC(Blen*sizeof(mpoly_heap_t));
    exps = (ulong *) TMP_ALLOC(Blen*N*sizeof(ulong));
    exp_list = (ulong **) TMP_ALLOC(Blen*sizeof(ulong *));
*/

    toff = 0;
    hind = (slong *)(tbase + toff);         toff += Blen*sizeof(slong);
    store = store_base = (slong *) (tbase + toff);           toff += 2*Blen*sizeof(slong);
    heap = (mpoly_heap_s *)(tbase + toff);  toff += (Blen + 1)*sizeof(mpoly_heap_s);
    chain = (mpoly_heap_t *)(tbase + toff); toff += Blen*sizeof(mpoly_heap_t);
    exps = (ulong *)(tbase + toff);         toff +=  Blen*N*sizeof(ulong);
    exp_list = (ulong **)(tbase + toff);    toff +=  Blen*sizeof(ulong *);
    FLINT_ASSERT(toff <= S->big_mem_alloc);

    exp_next = 0;
    for (i = 0; i < Blen; i++)
        exp_list[i] = exps + i*N;

    for (i = 0; i < Blen; i++)
        hind[i] = 1;

    /* mask with high bit set in each word of each field of exponent vector */
    mask = 0;
    for (i = 0; i < FLINT_BITS/bits; i++)
        mask = (mask << bits) + (UWORD(1) << (bits - 1));

    Qlen = WORD(0);

    /* s is the number of terms * (latest quotient) we should put into heap */
    s = Blen;

    /* insert (-1, 0, exp2[0]) into heap */
    heap_len = 2;
    x = chain + 0;
    x->i = -WORD(1);
    x->j = 0;
    x->next = NULL;
    heap[1].next = x;
    heap[1].exp = exp_list[exp_next++];

    FLINT_ASSERT(mpoly_monomial_cmp(Aexp + N*0, S->emin, N, S->cmpmask) >= 0);

    mpoly_monomial_set(heap[1].exp, Aexp + N*0, N);

    /* precompute leading cofficient info */
    lc_minus_inv = S->ctx->ffinfo->mod.n - nmod_inv(Bcoeff[0], S->ctx->ffinfo->mod);

    while (heap_len > 1)
    {
        exp = heap[1].exp;

        if (bits <= FLINT_BITS)
        {
            if (mpoly_monomial_overflows(exp, N, mask))
            {
                goto not_exact_division;
            }
        } else
        {
            if (mpoly_monomial_overflows_mp(exp, N, bits))
            {
                goto not_exact_division;
            }
        }

        FLINT_ASSERT(mpoly_monomial_cmp(exp, S->emin, N, S->cmpmask) >= 0);

        _nmod_mpoly_fit_length(&Qcoeff, &Qexp, &Qalloc, Qlen + 1, N);

        if (bits <= FLINT_BITS)
            lt_divides = mpoly_monomial_divides(Qexp + N*Qlen, exp, Bexp + N*0, N, mask);
        else
            lt_divides = mpoly_monomial_divides_mp(Qexp + N*Qlen, exp, Bexp + N*0, N, bits);

        acc0 = acc1 = acc2 = 0;
        do
        {
            exp_list[--exp_next] = heap[1].exp;
            x = _mpoly_heap_pop(heap, &heap_len, N, S->cmpmask);
            do
            {
                *store++ = x->i;
                *store++ = x->j;
                if (x->i != -WORD(1))
                    hind[x->i] |= WORD(1);

                if (x->i == -WORD(1))
                {
                    add_sssaaaaaa(acc2, acc1, acc0, acc2, acc1, acc0,
                                 WORD(0), WORD(0), S->ctx->ffinfo->mod.n - Acoeff[x->j]);
                } else
                {
                    umul_ppmm(pp1, pp0, Bcoeff[x->i], Qcoeff[x->j]);
                    add_sssaaaaaa(acc2, acc1, acc0, acc2, acc1, acc0, WORD(0), pp1, pp0);                    
                }
            } while ((x = x->next) != NULL);
        } while (heap_len > 1 && mpoly_monomial_equal(heap[1].exp, exp, N));

        NMOD_RED3(Qcoeff[Qlen], acc2, acc1, acc0, S->ctx->ffinfo->mod);

        /* process nodes taken from the heap */
        while (store > store_base)
        {
            j = *--store;
            i = *--store;

            if (i == -WORD(1))
            {
                /* take next dividend term */
                if (j + 1 < Alen)
                {
                    x = chain + 0;
                    x->i = i;
                    x->j = j + 1;
                    x->next = NULL;
                    mpoly_monomial_set(exp_list[exp_next], Aexp + x->j*N, N);

                    FLINT_ASSERT(mpoly_monomial_cmp(exp_list[exp_next], S->emin, N, S->cmpmask) >= 0);

                    if (!_mpoly_heap_insert(heap, exp_list[exp_next++], x,
                                      &next_loc, &heap_len, N, S->cmpmask))
                        exp_next--;
                }
            } else
            {
                /* should we go up */
                if (  (i + 1 < Blen)
                   && (hind[i + 1] == 2*j + 1)
                   )
                {
                    x = chain + i + 1;
                    x->i = i + 1;
                    x->j = j;
                    x->next = NULL;
                    hind[x->i] = 2*(x->j + 1) + 0;

                    mpoly_monomial_add_mp(exp_list[exp_next], Bexp + N*x->i,
                                                              Qexp + N*x->j, N);

                    if (mpoly_monomial_cmp(exp_list[exp_next], S->emin, N, S->cmpmask) >= 0)
                    {
                        if (!_mpoly_heap_insert(heap, exp_list[exp_next++], x,
                                          &next_loc, &heap_len, N, S->cmpmask))
                            exp_next--;
                    }
                    else
                    {
                        hind[x->i] |= 1;
                    }
                }
                /* should we go up? */
                if (j + 1 == Qlen)
                {
                    s++;
                } else if (  ((hind[i] & 1) == 1)
                          && ((i == 1) || (hind[i - 1] >= 2*(j + 2) + 1))
                          )
                {
                    x = chain + i;
                    x->i = i;
                    x->j = j + 1;
                    x->next = NULL;
                    hind[x->i] = 2*(x->j + 1) + 0;

                    mpoly_monomial_add_mp(exp_list[exp_next], Bexp + N*x->i,
                                                              Qexp + N*x->j, N);

                    if (mpoly_monomial_cmp(exp_list[exp_next], S->emin, N, S->cmpmask) >= 0)
                    {
                        if (!_mpoly_heap_insert(heap, exp_list[exp_next++], x,
                                          &next_loc, &heap_len, N, S->cmpmask))
                            exp_next--;
                    }
                    else
                    {
                        hind[x->i] |= 1;
                    }
                }
            }
        }

        Qcoeff[Qlen] = nmod_mul(Qcoeff[Qlen], lc_minus_inv, S->ctx->ffinfo->mod);
        if (Qcoeff[Qlen] == 0)
        {
            continue;
        }

        if (!lt_divides)
        {
            goto not_exact_division;
        }

        if (s > 1)
        {
            i = 1;
            x = chain + i;
            x->i = i;
            x->j = Qlen;
            x->next = NULL;
            hind[x->i] = 2*(x->j + 1) + 0;

            mpoly_monomial_add_mp(exp_list[exp_next], Bexp + N*x->i, Qexp + N*x->j, N);

            if (mpoly_monomial_cmp(exp_list[exp_next], S->emin, N, S->cmpmask) >= 0)
            {

                if (!_mpoly_heap_insert(heap, exp_list[exp_next++], x,
                                      &next_loc, &heap_len, N, S->cmpmask))
                    exp_next--;
            }
            else
            {
                hind[x->i] |= 1;                
            }
        }
        s = 1;      
        Qlen++;
    }


cleanup:

    *Qalloc_ = Qalloc;
    *Qcoeff_ = Qcoeff;
    *Qexp_ = Qexp;

    return Qlen;

not_exact_division:
    Qlen = 0;
    goto cleanup;
}


slong Lfind_exp(ulong * exp, slong a, const Lholder_t H)
{
    slong N = H->N;
    slong b = H->polyA->length;
    const ulong * Aexp = H->polyA->exps;

try_again:
    FLINT_ASSERT(b >= a);

    FLINT_ASSERT(a > 0);
    FLINT_ASSERT(mpoly_monomial_cmp(Aexp + N*(a - 1), exp, N, H->cmpmask) >= 0);
    FLINT_ASSERT(b >= H->polyA->length
                  ||  mpoly_monomial_cmp(Aexp + N*b, exp, N, H->cmpmask) < 0);

    if (b - a < 5)
    {
        slong i = a;
        while (i < b
                && mpoly_monomial_cmp(Aexp + N*i, exp, N, H->cmpmask) >= 0)
        {
            i++;
        }
        return i;
    }
    else
    {
        slong c = a + (b - a)/2;
        if (mpoly_monomial_cmp(Aexp + N*c, exp, N, H->cmpmask) < 0)
        {
            b = c;
        }
        else
        {
            a = c;
        }
        goto try_again;
    }
}


slong _nmod_mpoly_sub(ulong * coeff1,       ulong * exp1,
                const ulong * coeff2, const ulong * exp2, slong len2,
                const ulong * coeff3, const ulong * exp3, slong len3,
                   slong N, const ulong * cmpmask, const nmodf_ctx_t fctx);

void mulsub(worker_arg_t W, Lchunk_t L, slong q_prev_length)
{
    Lholder_struct * H = W->H;
    slong N = H->N;
    nmod_mpoly_struct * C = L->polyC;
    const nmod_mpoly_struct * B = H->polyB;
    const nmod_mpoly_struct * A = H->polyA;
    nmod_mpoly_ts_struct * Q = H->polyQ;
    nmod_mpoly_struct * T1 = W->polyT1;
    nmod_mpoly_struct * T2 = W->polyT2;

    {
        nmod_mpoly_stripe_struct * S = W->S;
        S->startidx = &L->startidx;
        S->endidx = &L->endidx;
        S->emin = L->emin;
        S->emax = L->emax;
        S->upperclosed = L->upperclosed;
        if (N == 1)
        {
            T1->length = _nmod_mpoly_mul_stripe1(&T1->coeffs, &T1->exps, &T1->alloc,
                    Q->coeffs + L->mq, Q->exps + N*L->mq, q_prev_length - L->mq,
                    B->coeffs, B->exps, B->length,  S);
        }
        else
        {
            T1->length = _nmod_mpoly_mul_stripe(&T1->coeffs, &T1->exps, &T1->alloc,
                    Q->coeffs + L->mq, Q->exps + N*L->mq, q_prev_length - L->mq,
                    B->coeffs, B->exps, B->length,  S);
        }
    }

    if (T1->length > 0)
    {
        if (L->Cinited)
        {
            nmod_mpoly_fit_length(T2, C->length + T1->length, H->ctx);
            T2->length = _nmod_mpoly_sub(T2->coeffs, T2->exps, 
                            C->coeffs, C->exps, C->length,
                            T1->coeffs, T1->exps, T1->length,
                                        N, H->cmpmask, H->ctx->ffinfo);
            nmod_mpoly_swap(C, T2, H->ctx);
        }
        else
        {
            slong startidx, stopidx;
            if (L->upperclosed)
            {
                startidx = 0;
                stopidx = Lfind_exp(L->emin, 1, H);
            }
            else
            {
                startidx = Lfind_exp(L->emax, 1, H);
                stopidx = Lfind_exp(L->emin, startidx, H);
            }

            L->Cinited = 1;
            nmod_mpoly_init2(C, T1->length + stopidx - startidx, H->ctx);
            nmod_mpoly_fit_bits(C, H->bits, H->ctx);
            C->bits = H->bits;
            C->length = _nmod_mpoly_sub(C->coeffs, C->exps, 
                        A->coeffs + startidx, A->exps + N*startidx, stopidx - startidx,
                        T1->coeffs, T1->exps, T1->length,
                                        N, H->cmpmask, H->ctx->ffinfo);
        }
    }

    L->mq = q_prev_length;
}

void tryproc(worker_arg_t W, Lchunk_t L)
{
    Lholder_struct * H = W->H;
    slong i;
    slong N = H->N;
    nmod_mpoly_struct * C = L->polyC;
    slong q_prev_length;
    ulong mask;
    const nmod_mpoly_struct * B = H->polyB;
    const nmod_mpoly_struct * A = H->polyA;
    nmod_mpoly_ts_struct * Q = H->polyQ;
    nmod_mpoly_struct * T2 = W->polyT2;

    mask = 0;
    for (i = 0; i < FLINT_BITS/H->bits; i++)
        mask = (mask << H->bits) + (UWORD(1) << (H->bits - 1));

    /* return if this section has already finished processing */
    if (L->mq < 0)
    {
        return;
    }

    /* process more quotient terms if available */
    q_prev_length = Q->length;
    if (q_prev_length > L->mq)
    {
        if (L->producer == 0 && q_prev_length - L->mq < 20)
            return;

#if defined PROFILE_THIS
        vec_slong_push_back(W->time_data, 4*L->idx + 0);
        vec_slong_push_back(W->time_data, timeit_query_wall(H->timer));
#endif 
        mulsub(W, L, q_prev_length);
#if defined PROFILE_THIS
        vec_slong_push_back(W->time_data, 4*L->idx + 1);
        vec_slong_push_back(W->time_data, timeit_query_wall(H->timer));
#endif 

    }

    if (L->producer == 1)
    {
        Lchunk_struct * next;
        mp_limb_t * Rcoeff;
        ulong * Rexp;
        slong Rlen;

#if defined PROFILE_THIS
        vec_slong_push_back(W->time_data, 4*L->idx + 2);
        vec_slong_push_back(W->time_data, timeit_query_wall(H->timer));
#endif 

        /* process the remaining quotient terms */
        q_prev_length = Q->length;
        if (q_prev_length > L->mq)
        {
            mulsub(W, L, q_prev_length);
        }

        /* find location of remaining terms */
        if (L->Cinited)
        {
            Rlen = C->length;
            Rexp = C->exps;
            Rcoeff = C->coeffs;
        }
        else
        {
            slong startidx, stopidx;
            if (L->upperclosed)
            {
                startidx = 0;
                stopidx = Lfind_exp(L->emin, 1, H);
            }
            else
            {
                startidx = Lfind_exp(L->emax, 1, H);
                stopidx = Lfind_exp(L->emin, startidx, H);
            }
            Rlen = stopidx - startidx;
            Rcoeff = A->coeffs + startidx;
            Rexp = A->exps + N*startidx;
        }

        /* if we have remaining terms, add to quotient  */
        if (Rlen > 0)
        {
            nmod_mpoly_stripe_struct * S = W->S;
            S->startidx = &L->startidx;
            S->endidx = &L->endidx;
            S->emin = L->emin;
            S->emax = L->emax;
            S->upperclosed = L->upperclosed;
            if (N == 1)
            {
                T2->length = _nmod_mpoly_divides_stripe1(
                            &T2->coeffs, &T2->exps, &T2->alloc,
                               Rcoeff, Rexp, Rlen,
                               B->coeffs, B->exps, B->length,  S);
            }
            else
            {
                T2->length = _nmod_mpoly_divides_stripe(
                            &T2->coeffs, &T2->exps, &T2->alloc,
                               Rcoeff, Rexp, Rlen,
                               B->coeffs, B->exps, B->length,  S);
            }
            if (T2->length == 0)
            {
#if defined PROFILE_THIS
                vec_slong_push_back(W->time_data, 4*L->idx + 3);
                vec_slong_push_back(W->time_data, timeit_query_wall(H->timer));
#endif 
                H->failed = 1;
                return;
            }
            else
            {
                nmod_mpoly_ts_append(H->polyQ, T2->coeffs, T2->exps, T2->length, H->N);
            }
        }

#if defined PROFILE_THIS
        vec_slong_push_back(W->time_data, 4*L->idx + 3);
        vec_slong_push_back(W->time_data, timeit_query_wall(H->timer));
#endif 

        next = L->next;
        H->length--;
        H->cur = next;

        if (next != NULL)
        {
            next->producer = 1;
        }

        L->producer = 0;
        L->mq = -1;
    }

    return;
}


void worker_loop(void * varg)
{
    worker_arg_struct * W = (worker_arg_struct *) varg;
    Lholder_struct * H = W->H;
    nmod_mpoly_stripe_struct * S = W->S;
    const nmod_mpoly_struct * B = H->polyB;
    nmod_mpoly_struct * T1 = W->polyT1;
    nmod_mpoly_struct * T2 = W->polyT2;
    slong N = H->N;
    slong Blen = B->length;

    /* initialize stripe working memory */
    S->N = N;
    S->bits = H->bits;
    S->ctx = H->ctx;
    S->cmpmask = H->cmpmask;
    S->big_mem_alloc = 0;
    if (N == 1)
    {
        S->big_mem_alloc += Blen*sizeof(slong);
        S->big_mem_alloc += Blen*sizeof(slong);
        S->big_mem_alloc += 2*Blen*sizeof(slong);
        S->big_mem_alloc += (Blen + 1)*sizeof(mpoly_heap1_s);
        S->big_mem_alloc += Blen*sizeof(mpoly_heap_t);
    }
    else
    {
        S->big_mem_alloc = 0;
        S->big_mem_alloc += Blen*sizeof(slong);
        S->big_mem_alloc += Blen*sizeof(slong);
        S->big_mem_alloc += 2*Blen*sizeof(slong);
        S->big_mem_alloc += (Blen + 1)*sizeof(mpoly_heap_s);
        S->big_mem_alloc += Blen*sizeof(mpoly_heap_t);
        S->big_mem_alloc +=  Blen*N*sizeof(ulong);
        S->big_mem_alloc +=  Blen*sizeof(ulong *);
        S->big_mem_alloc +=  N*sizeof(ulong);
    }
    S->big_mem = (char *) flint_malloc(S->big_mem_alloc);
    nmod_mpoly_init2(T1, 16, H->ctx);
    nmod_mpoly_fit_bits(T1, H->bits, H->ctx);
    T1->bits = H->bits;
    nmod_mpoly_init2(T2, 16, H->ctx);
    nmod_mpoly_fit_bits(T2, H->bits, H->ctx);
    T2->bits = H->bits;

    while (!H->failed)
    {
        Lchunk_struct * L;
        L = H->cur;

        if (L == NULL)
        {
            break;
        }
        while (L != NULL)
        {
            pthread_mutex_lock(&H->mutex);
            if (L->lock != -1)
            {
                L->lock = -1;
                pthread_mutex_unlock(&H->mutex);
                tryproc(W, L);
                pthread_mutex_lock(&H->mutex);
                L->lock = 0;
                pthread_mutex_unlock(&H->mutex);
                break;
            }
            else
            {
                pthread_mutex_unlock(&H->mutex);
            }

            L = L->next;
        }
    }

    nmod_mpoly_clear(T1, H->ctx);
    nmod_mpoly_clear(T2, H->ctx);
    flint_free(S->big_mem);

    return;
}


/* return 1 if quotient is exact */
int nmod_mpoly_divides_heap_threaded(nmod_mpoly_t Q,
                  const nmod_mpoly_t A, const nmod_mpoly_t B,
                                                    const nmod_mpoly_ctx_t ctx)
{
    int ret = 1;
    fmpz_mpoly_ctx_t zctx;
    fmpz_mpoly_t S;
    slong i, k, N;
    mp_bitcnt_t exp_bits;
    ulong * cmpmask;
    ulong * Aexp, * Bexp;
    int freeAexp = 0, freeBexp = 0;
    TMP_INIT;
/*
printf("******************************\n");
printf("div A: "); nmod_mpoly_print_pretty(A, vars, ctx); printf("\n");
printf("div B: "); nmod_mpoly_print_pretty(B, vars, ctx); printf("\n");
*/
    FLINT_ASSERT(Q != A);
    FLINT_ASSERT(Q != B);

    if (B->length < 2 || A->length < 2)
    {
        if (B->length == 0)
            flint_throw(FLINT_DIVZERO, "Divide by zero in nmod_mpoly_divides_monagan_pearce");

        if (A->length == 0)
        {
            nmod_mpoly_zero(Q, ctx);
            return 1;
        }
        return nmod_mpoly_divides_monagan_pearce(Q, A, B, ctx);
    }
/*
time_prod = time_add = time_div = 0;
*/
    

    TMP_START;

    exp_bits = MPOLY_MIN_BITS;
    exp_bits = FLINT_MAX(exp_bits, A->bits);
    exp_bits = FLINT_MAX(exp_bits, B->bits);
    exp_bits = mpoly_fix_bits(exp_bits, ctx->minfo);

    N = mpoly_words_per_exp(exp_bits, ctx->minfo);
    cmpmask = (ulong*) TMP_ALLOC(N*sizeof(ulong));
    mpoly_get_cmpmask(cmpmask, N, exp_bits, ctx->minfo);

    /* ensure input exponents packed to same size as output exponents */
    Aexp = A->exps;
    if (exp_bits > A->bits)
    {
        freeAexp = 1;
        Aexp = (ulong *) flint_malloc(N*A->length*sizeof(ulong));
        mpoly_repack_monomials(Aexp, exp_bits, A->exps, A->bits,
                                                        A->length, ctx->minfo);
    }

    Bexp = B->exps;
    if (exp_bits > B->bits)
    {
        freeBexp = 1;
        Bexp = (ulong *) flint_malloc(N*B->length*sizeof(ulong));
        mpoly_repack_monomials(Bexp, exp_bits, B->exps, B->bits,
                                                    B->length, ctx->minfo);
    }

    fmpz_mpoly_ctx_init(zctx, ctx->minfo->nvars, ctx->minfo->ord);
    fmpz_mpoly_init(S, zctx);

    FLINT_ASSERT(exp_bits <= FLINT_BITS);

    if (select_exps(S, zctx, Aexp, A->length, Bexp, B->length, exp_bits))
    {
        ret = 0;
        nmod_mpoly_zero(Q, ctx);
        goto cleanup1;
    }

    {
#if defined PROFILE_THIS
        slong idx = 0;
#endif
        mp_limb_t qcoeff;
        ulong * texps, * qexps;
        mp_limb_t lc_inv;
        Lholder_t H;

        threadpool_threadhandle handles[WORKER_COUNT];
        slong num_workers;
        
        Lholder_init(H);

        H->polyA->coeffs = A->coeffs;
        H->polyA->exps = Aexp;
        H->polyA->bits = exp_bits;
        H->polyA->length = A->length;
        H->polyA->alloc = A->alloc;

        H->polyB->coeffs = B->coeffs;
        H->polyB->exps = Bexp;
        H->polyB->bits = exp_bits;
        H->polyB->length = B->length;
        H->polyB->alloc = B->alloc;

        H->ctx = ctx;
        H->bits = exp_bits;
        H->N = N;
        H->cmpmask = cmpmask;
        H->failed = 0;

        for (i = 0; i + 1 < S->length; i++)
        {
            Lchunk_struct * L;

            L = (Lchunk_struct *) malloc(sizeof(Lchunk_struct));
            L->ma = 0;
            L->mq = 0;
            L->emax = S->exps + N*i;
            L->emin = S->exps + N*(i + 1);
            L->upperclosed = 0;
            L->startidx = B->length;
            L->endidx = B->length;
            L->producer = 0;
            L->Cinited = 0;
            L->lock = -2;
#if defined PROFILE_THIS
            L->idx = idx++;
#endif
            Lholder_add_chunk(H, L);
        }

        H->head->upperclosed = 1;
        H->head->producer = 1;
        H->cur = H->head;

        texps = (ulong *) TMP_ALLOC(N*sizeof(ulong));
        qexps = (ulong *) TMP_ALLOC(N*sizeof(ulong));

        mpoly_monomial_sub_mp(qexps + N*0, Aexp + N*0, Bexp + N*0, N);
        lc_inv = nmod_inv(B->coeffs[0], ctx->ffinfo->mod);
        qcoeff = nmod_mul(lc_inv, A->coeffs[0], ctx->ffinfo->mod);

        nmod_mpoly_ts_init(H->polyQ, &qcoeff, qexps, 1, H->bits, H->N);

        mpoly_monomial_add_mp(texps, qexps + N*0, Bexp + N*1, N);

        k = 1;
        while (k < A->length && mpoly_monomial_gt(texps, Aexp + N*k, N, cmpmask))
        {
            mpoly_monomial_sub_mp(qexps, Aexp + N*k, Bexp + N*0, N);
            qcoeff = nmod_mul(lc_inv, A->coeffs[k], ctx->ffinfo->mod);
            nmod_mpoly_ts_append(H->polyQ, &qcoeff, qexps, 1, H->N);
            k++;
        }

        pthread_mutex_init(&H->mutex, NULL);

        threadpool_init(global_thread_pool, WORKER_COUNT);
        num_workers = threadpool_request(global_thread_pool, handles, WORKER_COUNT);

        {
            worker_arg_struct * worker_args;
            worker_args = (worker_arg_struct *) flint_malloc((num_workers + 1)*sizeof(worker_arg_t));

#if defined PROFILE_THIS
            for (i = 0; i < num_workers + 1; i++)
            {
                vec_slong_init((worker_args + i)->time_data);
            }
            timeit_start(H->timer);
#endif

            for (i = 0; i < num_workers; i++)
            {
                (worker_args + i)->H = H;
                threadpool_wake(global_thread_pool, handles[i], worker_loop, worker_args + i);
            }
            (worker_args + num_workers)->H = H;
            worker_loop(worker_args + num_workers);
            for (i = 0; i < num_workers; i++)
            {
                threadpool_wait(global_thread_pool, handles[i]);
                threadpool_giveback(global_thread_pool, handles[i]);
            }

#if defined PROFILE_THIS
            timeit_stop(H->timer);
            flint_printf("data = [");
            for (i = 0; i < num_workers + 1; i++)
            {
                flint_printf("[%wd,", i);
                vec_slong_print((worker_args + i)->time_data);
                flint_printf("],\n");
                vec_slong_clear((worker_args + i)->time_data);

            }
            flint_printf("%wd]\n", H->timer->wall);
#endif

            flint_free(worker_args);
        }


        threadpool_clear(global_thread_pool);

        pthread_mutex_destroy(&H->mutex);

        ret = Lholder_clear(Q, H);
/*
flint_printf("time_add: %wd\n", time_add);
flint_printf("time_prod: %wd\n", time_prod);
flint_printf("time_div: %wd\n", time_div);
*/
    }


cleanup1:
    fmpz_mpoly_clear(S, zctx);
    fmpz_mpoly_ctx_clear(zctx);

    if (freeAexp)
        flint_free(Aexp);

    if (freeBexp)
        flint_free(Bexp);

    TMP_END;

    return ret;
}
