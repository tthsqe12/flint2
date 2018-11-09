/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/


#include "nmod_mpoly.h"
#include "fmpz_mpoly.h"

const char * vars[] = {"x","y","z","t","u"};


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
} Lchunk_struct;

typedef Lchunk_struct Lchunk_t[1];



typedef struct
{
    Lchunk_struct * head;
    Lchunk_struct * tail;
    Lchunk_struct * cur;
    nmod_mpoly_t polyA;
    nmod_mpoly_t polyB;
    const nmod_mpoly_ctx_struct * ctx;
    slong length;
    slong N;
    mp_bitcnt_t bits;
    ulong * cmpmask;
    int failed;
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
static void Lholder_clear(Lholder_t H)
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
}

static void Lholder_add_chunk(Lholder_t H, Lchunk_t L)
{
/*
printf("add chunk called on %p\n",L);
*/
    L->next = NULL;

    if (H->tail == NULL)
    {
/*
printf("adding to empty\n");
*/
        FLINT_ASSERT(H->head == NULL);
        H->tail = L;
        H->head = L;
    }
    else
    {
        Lchunk_struct * tail = H->tail;
/*
printf("adding to end tail = %p\n",tail);
*/
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
    slong ns = 10;
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

/*
flint_printf("N: %wd   nvars: %wd\n", N, zctx->minfo->nvars);
flint_printf("S: "); fmpz_mpoly_print_pretty(S, vars, zctx); flint_printf("\n");
*/

/*
printf("here1\n");
*/
    Slen = 0;
    for (i = 0; i < Alen; i += Astep)
    {
        mpoly_monomial_set(Sexp + N*Slen, Aexp + N*i, N);
        fmpz_one(Scoeff + Slen);
        Slen++;   
    }
/*
flint_printf("Slen: %wd\n", Slen);
*/
    _fmpz_mpoly_set_length(S, Slen, zctx);
/*
flint_printf("S: "); fmpz_mpoly_print_pretty(S, vars, zctx); flint_printf("\n");
printf("here2\n");
*/
    T0 = (ulong *) TMP_ALLOC(N*sizeof(ulong));
    T1 = (ulong *) TMP_ALLOC(N*sizeof(ulong));
    mpoly_monomial_sub_mp(T0, Aexp + N*0, Bexp + N*0, N);
    mpoly_monomial_sub_mp(T1, Aexp + N*(Alen - 1), Bexp + N*(Blen - 1), N);
    if (bits <= FLINT_BITS)
    {
        if (   mpoly_monomial_overflows(T0, N, mask)
            || mpoly_monomial_overflows(T1, N, mask))
        {
/*
flint_printf(" first division failed sp\n");
*/
            failure = 1;
            goto cleanup;
        }
    }
    else
    {
        if (   mpoly_monomial_overflows_mp(T0, N, bits)
            || mpoly_monomial_overflows_mp(T1, N, bits))
        {
/*
flint_printf(" first division failed mp\n");
*/
            failure = 1;
            goto cleanup;
        }
    }
/*
printf("here3\n");
*/
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
/*
printf("here4\n");
*/
    mpoly_monomial_zero(Sexp + N*Slen, N);
    fmpz_one(Scoeff + Slen);
    Slen++;

    FLINT_ASSERT(Slen < tot);

    _fmpz_mpoly_set_length(S, Slen, zctx);
/*
flint_printf("before sort S: "); fmpz_mpoly_print_pretty(S, vars, zctx); flint_printf("\n");
*/
    fmpz_mpoly_sort_terms(S, zctx);
    fmpz_mpoly_combine_like_terms(S, zctx);
/*
flint_printf(" after sort S: "); fmpz_mpoly_print_pretty(S, vars, zctx); flint_printf("\n");
*/
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
/*
printf("A exp: %016llx   emax: %016llx\n", (A->exps + N*i)[0], L->emax[0]);
*/
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
/*
printf("A exp: %016llx   emin: %016llx\n", (A->exps + N*i)[0], L->emin[0]);
*/
        if (mpoly_monomial_gt(A->exps + N*i, emin, N, H->cmpmask))
        {
            break;
        }
/*
printf("pushing term\n");
*/
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
/*
flint_printf("C: "); nmod_mpoly_print_pretty(C, vars, H->ctx); flint_printf("\n");
*/
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

/*
A = B * C
*/
slong _nmod_mpoly_mul_stripe(mp_limb_t ** Acoeff_, ulong ** Aexp_, slong * Aalloc_,
                 const mp_limb_t * Bcoeff, const ulong * Bexp, slong Blen,
                 const mp_limb_t * Ccoeff, const ulong * Cexp, slong Clen,
                                                       Lchunk_t L, Lholder_t H)
{
    int upperclosed;
    slong startidx, endidx;
    slong N = H->N;
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
    slong * starts, * ends;
    ulong * texp, * emax, * emin;
    slong * hind;
    TMP_INIT;

    TMP_START;

    hind = (slong *) TMP_ALLOC(Blen*sizeof(slong));
    starts = (slong *) TMP_ALLOC(Blen*sizeof(slong));
    ends = (slong *) TMP_ALLOC(Blen*sizeof(slong));
    texp = (ulong *) TMP_ALLOC(N*sizeof(ulong));
    heap = (mpoly_heap_s *) TMP_ALLOC((Blen + 1)*sizeof(mpoly_heap_s));
    chain = (mpoly_heap_t *) TMP_ALLOC(Blen*sizeof(mpoly_heap_t));
    Q = (slong *) TMP_ALLOC(2*Blen*sizeof(slong));
    exps = (ulong *) TMP_ALLOC(Blen*N*sizeof(ulong));
    exp_list = (ulong **) TMP_ALLOC(Blen*sizeof(ulong *));
    exp_next = 0;

    startidx = L->startidx;
    endidx = L->endidx;
    upperclosed = L->upperclosed;
    emax = L->emax;
    emin = L->emin;

    for (i = 0; i < Blen; i++)
        exp_list[i] = exps + i*N;

    /* put all the starting nodes on the heap */
    for (i = 0; i < Blen; i++)
    {
        if (startidx < Clen)
        {
            mpoly_monomial_add_mp(texp, Bexp + N*i, Cexp + N*startidx, N);
            FLINT_ASSERT(mpoly_monomial_cmp(emax, texp, N, H->cmpmask) > -upperclosed);
        }
        while (startidx > 0)
        {
            mpoly_monomial_add_mp(texp, Bexp + N*i, Cexp + N*(startidx - 1), N);
            if (mpoly_monomial_cmp(emax, texp, N, H->cmpmask) <= -upperclosed)
            {
                break;
            }
            startidx--;
        }

        if (endidx < Clen)
        {
            mpoly_monomial_add_mp(texp, Bexp + N*i, Cexp + N*endidx, N);
            FLINT_ASSERT(mpoly_monomial_cmp(emin, texp, N, H->cmpmask) > 0);
        }
        while (endidx > 0)
        {
            mpoly_monomial_add_mp(texp, Bexp + N*i, Cexp + N*(endidx - 1), N);
            if (mpoly_monomial_cmp(emin, texp, N, H->cmpmask) <= 0)
            {
                break;
            }
            endidx--;
        }

        starts[i] = startidx;
        ends[i] = endidx;

        hind[i] = 2*startidx + 1;

        if (  (startidx < endidx)
           && (  (i == 0)
              || (startidx < starts[i - 1])
              )
           )
        {
            x = chain + i;
            x->i = i;
            x->j = startidx;
            x->next = NULL;
            hind[x->i] = 2*(x->j + 1) + 0;

            mpoly_monomial_add(exp_list[exp_next], Bexp + x->i*N, Cexp + x->j*N, N);

            if (!_mpoly_heap_insert(heap, exp_list[exp_next++], x,
                                      &next_loc, &heap_len, N, H->cmpmask))
               exp_next--;
        }
    }

    L->startidx = startidx;
    L->endidx = endidx;

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

            x = _mpoly_heap_pop(heap, &heap_len, N, H->cmpmask);

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

        NMOD_RED3(Acoeff[Alen], acc2, acc1, acc0, H->ctx->ffinfo->mod);
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
                                      &next_loc, &heap_len, N, H->cmpmask))
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
                                      &next_loc, &heap_len, N, H->cmpmask))
                    exp_next--;
            }
        }
    }

    *Acoeff_ = Acoeff;
    *Aexp_ = Aexp;
    *Aalloc_ = Aalloc;

    TMP_END;

    return Alen;
}


slong _nmod_mpoly_divides_stripe(
                     mp_limb_t ** Qcoeff_,      ulong ** Qexp_, slong * Qalloc_,
                const mp_limb_t * Acoeff, const ulong * Aexp, slong Alen,
                const mp_limb_t * Bcoeff, const ulong * Bexp, slong Blen,
                                                       Lchunk_t L, Lholder_t H)
{
    mp_bitcnt_t bits = H->bits;
    slong N = H->N;
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
    TMP_INIT;


    TMP_START;

    FLINT_ASSERT(Alen > 0);
    FLINT_ASSERT(Blen > 0);

    /* alloc array of heap nodes which can be chained together */
    next_loc = Blen + 4;   /* something bigger than heap can ever be */
    heap = (mpoly_heap_s *) TMP_ALLOC((Blen + 1)*sizeof(mpoly_heap_s));
    chain = (mpoly_heap_t *) TMP_ALLOC(Blen*sizeof(mpoly_heap_t));
    store = store_base = (slong *) TMP_ALLOC(2*Blen*sizeof(mpoly_heap_t *));

    /* array of exponent vectors, each of "N" words */
    exps = (ulong *) TMP_ALLOC(Blen*N*sizeof(ulong));
    /* list of pointers to available exponent vectors */
    exp_list = (ulong **) TMP_ALLOC(Blen*sizeof(ulong *));
    /* set up list of available exponent vectors */
    exp_next = 0;
    for (i = 0; i < Blen; i++)
        exp_list[i] = exps + i*N;

    /* space for flagged heap indicies */
    hind = (slong *) TMP_ALLOC(Blen*sizeof(slong));
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

    FLINT_ASSERT(mpoly_monomial_cmp(Aexp + N*0, L->emin, N, H->cmpmask) >= 0);

    mpoly_monomial_set(heap[1].exp, Aexp + N*0, N);

    /* precompute leading cofficient info */
    lc_minus_inv = H->ctx->ffinfo->mod.n - nmod_inv(Bcoeff[0], H->ctx->ffinfo->mod);

    while (heap_len > 1)
    {
        exp = heap[1].exp;

        if (bits <= FLINT_BITS)
        {
            if (mpoly_monomial_overflows(exp, N, mask))
                goto not_exact_division;
        } else
        {
            if (mpoly_monomial_overflows_mp(exp, N, bits))
                goto not_exact_division;
        }

        FLINT_ASSERT(mpoly_monomial_cmp(exp, L->emin, N, H->cmpmask) >= 0);

        _nmod_mpoly_fit_length(&Qcoeff, &Qexp, &Qalloc, Qlen + 1, N);

        if (bits <= FLINT_BITS)
            lt_divides = mpoly_monomial_divides(Qexp + N*Qlen, exp, Bexp + N*0, N, mask);
        else
            lt_divides = mpoly_monomial_divides_mp(Qexp + N*Qlen, exp, Bexp + N*0, N, bits);

        acc0 = acc1 = acc2 = 0;
        do
        {
            exp_list[--exp_next] = heap[1].exp;
            x = _mpoly_heap_pop(heap, &heap_len, N, H->cmpmask);
            do
            {
                *store++ = x->i;
                *store++ = x->j;
                if (x->i != -WORD(1))
                    hind[x->i] |= WORD(1);

                if (x->i == -WORD(1))
                {
                    add_sssaaaaaa(acc2, acc1, acc0, acc2, acc1, acc0,
                                 WORD(0), WORD(0), H->ctx->ffinfo->mod.n - Acoeff[x->j]);
                } else
                {
                    umul_ppmm(pp1, pp0, Bcoeff[x->i], Qcoeff[x->j]);
                    add_sssaaaaaa(acc2, acc1, acc0, acc2, acc1, acc0, WORD(0), pp1, pp0);                    
                }
            } while ((x = x->next) != NULL);
        } while (heap_len > 1 && mpoly_monomial_equal(heap[1].exp, exp, N));

        NMOD_RED3(Qcoeff[Qlen], acc2, acc1, acc0, H->ctx->ffinfo->mod);

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

                    FLINT_ASSERT(mpoly_monomial_cmp(exp_list[exp_next], L->emin, N, H->cmpmask) >= 0);

                    if (!_mpoly_heap_insert(heap, exp_list[exp_next++], x,
                                      &next_loc, &heap_len, N, H->cmpmask))
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

                    if (mpoly_monomial_cmp(exp_list[exp_next], L->emin, N, H->cmpmask) >= 0)
                    {
                        if (!_mpoly_heap_insert(heap, exp_list[exp_next++], x,
                                          &next_loc, &heap_len, N, H->cmpmask))
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

                    if (mpoly_monomial_cmp(exp_list[exp_next], L->emin, N, H->cmpmask) >= 0)
                    {
                        if (!_mpoly_heap_insert(heap, exp_list[exp_next++], x,
                                          &next_loc, &heap_len, N, H->cmpmask))
                            exp_next--;
                    }
                    else
                    {
                        hind[x->i] |= 1;
                    }
                }
            }
        }

        Qcoeff[Qlen] = nmod_mul(Qcoeff[Qlen], lc_minus_inv, H->ctx->ffinfo->mod);

        if (Qcoeff[Qlen] == 0)
        {
            continue;
        }

        if (!lt_divides ||
                mpoly_monomial_gt(exp, Aexp + N*(Alen - 1), N, H->cmpmask))
            goto not_exact_division;

        if (s > 1)
        {
            i = 1;
            x = chain + i;
            x->i = i;
            x->j = Qlen;
            x->next = NULL;
            hind[x->i] = 2*(x->j + 1) + 0;

            mpoly_monomial_add_mp(exp_list[exp_next], Bexp + N*x->i, Qexp + N*x->j, N);

            if (mpoly_monomial_cmp(exp_list[exp_next], L->emin, N, H->cmpmask) >= 0)
            {

                if (!_mpoly_heap_insert(heap, exp_list[exp_next++], x,
                                      &next_loc, &heap_len, N, H->cmpmask))
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

    TMP_END;

    return Qlen;

not_exact_division:
    Qlen = 0;
    goto cleanup;
}


slong Lfind_exp(ulong * exp, slong i, const Lholder_t H)
{
    slong N = H->N;
    slong Alen = H->polyA->length;
    const ulong * Aexp = H->polyA->exps;

    FLINT_ASSERT(i > 0);
    FLINT_ASSERT(mpoly_monomial_cmp(Aexp + N*(i - 1), exp, N, H->cmpmask) >= 0);

    while (i < Alen
            && mpoly_monomial_cmp(Aexp + N*i, exp, N, H->cmpmask) >= 0)
    {
        i++;
    }
    return i;
}

void tryproc(Lholder_t H, Lchunk_t L, nmod_mpoly_t Q)
{
    slong i;
    slong N = H->N;
    nmod_mpoly_struct * C = L->polyC;
    ulong mask;
    const nmod_mpoly_struct * B = H->polyB;
    const nmod_mpoly_struct * A = H->polyA;

    mask = 0;
    for (i = 0; i < FLINT_BITS/H->bits; i++)
        mask = (mask << H->bits) + (UWORD(1) << (H->bits - 1));


flint_printf("tryproducer called\n");


    FLINT_ASSERT(L->mq >= 0);


    if (Q->length > L->mq)
    {
        nmod_mpoly_t T;
        nmod_mpoly_init2(T, 16, H->ctx);
        T->length = _nmod_mpoly_mul_stripe(&T->coeffs, &T->exps, &T->alloc,
                Q->coeffs + L->mq, Q->exps + N*L->mq, Q->length - L->mq,
                B->coeffs, B->exps, B->length,
                        L, H);

        if (L->Cinited)
        {
            nmod_mpoly_sub(C, C, T, H->ctx);
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
            nmod_mpoly_init2(C, T->length + stopidx - startidx, H->ctx);
            nmod_mpoly_fit_bits(C, H->bits, H->ctx);
            C->bits = H->bits;

            C->length = _nmod_mpoly_sub(C->coeffs, C->exps, 
                        A->coeffs + startidx, A->exps + N*startidx, stopidx - startidx,
                        T->coeffs, T->exps, T->length,
                                        N, H->cmpmask, H->ctx->ffinfo);
        }
        nmod_mpoly_clear(T, H->ctx);
    }
/*
flint_printf("all C: "); nmod_mpoly_print_pretty(C, vars, H->ctx); flint_printf("\n");
flint_printf("all Q: "); nmod_mpoly_print_pretty(Q, vars, H->ctx); flint_printf("\n");
*/
    L->mq = Q->length;

    if (L->producer == 1)
    {
        FLINT_ASSERT(L->Cinited);

        if (C->length > 0)
        {
            nmod_mpoly_t T;
            nmod_mpoly_init2(T, 16, H->ctx);
            nmod_mpoly_fit_bits(T, H->bits, H->ctx);
            T->bits = H->bits;

            T->length = _nmod_mpoly_divides_stripe(
                        &T->coeffs, &T->exps, &T->alloc,
                           C->coeffs, C->exps, C->length,
                           B->coeffs, B->exps, B->length,
                                                       L, H);

            if (T->length == 0)
            {
                H->failed = 1;
                nmod_mpoly_clear(T, H->ctx);
                return;
            }
            else
            {
                mp_limb_t * Qcoeff = Q->coeffs;
                ulong * Qexp = Q->exps;
                slong Qalloc =  Q->alloc;
                slong Qlen = Q->length;

                _nmod_mpoly_fit_length(&Qcoeff, &Qexp, &Qalloc, Qlen + T->length, N);
                for (i = 0; i < T->length; i++)
                {
                    Qcoeff[Qlen] = T->coeffs[i];
                    mpoly_monomial_set(Qexp + N*Qlen, T->exps + N*i, N);
                    Qlen++;
                }
                Q->exps = Qexp;
                Q->coeffs = Qcoeff;
                Q->alloc = Qalloc;
                Q->length = Qlen;
                nmod_mpoly_clear(T, H->ctx);
            }
        }

        L->producer = 0;
        L->mq = -1;

        H->length--;
        H->cur = L->next;
        if (H->cur != NULL)
        {
            H->cur->producer = 1;
        }
    }

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

    if (B->length == 0)
        flint_throw(FLINT_DIVZERO, "Divide by zero in nmod_mpoly_divides_monagan_pearce");

    if (A->length == 0)
    {
        nmod_mpoly_zero(Q, ctx);
        return 1;
    }

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
/*
flint_printf("S: "); fmpz_mpoly_print_pretty(S, vars, zctx); flint_printf("\n");
*/
    {
        ulong * texps;
        mp_limb_t lc_inv;
        Lholder_t H;
        flint_rand_t randstate;
        
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
            Lholder_add_chunk(H, L);
        }

        H->head->upperclosed = 1;
        H->head->producer = 1;
        H->cur = H->head;
/*
Lholder_print(H);
*/
        texps = (ulong *) TMP_ALLOC(N*sizeof(ulong));

        nmod_mpoly_fit_length(Q, 1, ctx);
        nmod_mpoly_fit_bits(Q, exp_bits, ctx);
        Q->bits = exp_bits;

        mpoly_monomial_sub_mp(Q->exps + N*0, Aexp + N*0, Bexp + N*0, N);
        lc_inv = nmod_inv(B->coeffs[0], ctx->ffinfo->mod);
        Q->coeffs[0] = nmod_mul(lc_inv, A->coeffs[0], ctx->ffinfo->mod);
        Q->length = (Q->coeffs[0] != UWORD(0));
/*
printf("1 Q: "); nmod_mpoly_print_pretty(Q, vars, ctx); printf("\n");
*/
        mpoly_monomial_add_mp(texps, Q->exps + N*0, Bexp + N*1, N);

        k = 1;
        while (k < A->length && mpoly_monomial_gt(texps, Aexp + N*k, N, cmpmask))
        {
            nmod_mpoly_fit_length(Q, k + 1, ctx);
            mpoly_monomial_sub_mp(Q->exps + N*k, Aexp + N*k, Bexp + N*0, N);
            Q->coeffs[k] = nmod_mul(lc_inv, A->coeffs[k], ctx->ffinfo->mod);
            Q->length++;
            k++;
        }
/*
printf("2 Q: "); nmod_mpoly_print_pretty(Q, vars, ctx); printf("\n");
*/
        flint_randinit(randstate);

        while (!H->failed && H->cur != NULL)
        {
            slong i;
            Lchunk_struct * L = H->cur;

            /* pick random */
            i = n_randint(randstate, H->length);
/*
flint_printf("**** i: %wd *********************\n", i);
*/
            L = H->cur;
            while (i > 0)
            {
                L = L->next;
                i--;
            }
            tryproc(H, L, Q);
        }

        if (H->failed)
        {
            ret = 0;
            nmod_mpoly_zero(Q, ctx);
        }

        flint_randclear(randstate);

        Lholder_clear(H);
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
