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
    const nmod_mpoly_struct * polyA;
    const nmod_mpoly_ctx_struct * ctx;
    slong length;
    slong N;
    mp_bitcnt_t bits;
    ulong * cmpmask;
} Lholder_struct;

typedef Lholder_struct Lholder_t[1];


static void Lholder_init(Lholder_t H)
{
    H->head = NULL;
    H->tail = NULL;
    H->cur = NULL;
    H->polyA = NULL;
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
    H->polyA = NULL;
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
    slong ns = 8;
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


flint_printf("N: %wd   nvars: %wd\n", N, zctx->minfo->nvars);
flint_printf("S: "); fmpz_mpoly_print_pretty(S, vars, zctx); flint_printf("\n");



printf("here1\n");

    Slen = 0;
    for (i = 0; i < Alen; i += Astep)
    {
        mpoly_monomial_set(Sexp + N*Slen, Aexp + N*i, N);
        fmpz_one(Scoeff + Slen);
        Slen++;   
    }

flint_printf("Slen: %wd\n", Slen);

    _fmpz_mpoly_set_length(S, Slen, zctx);
flint_printf("S: "); fmpz_mpoly_print_pretty(S, vars, zctx); flint_printf("\n");


printf("here2\n");

    T0 = (ulong *) TMP_ALLOC(N*sizeof(ulong));
    T1 = (ulong *) TMP_ALLOC(N*sizeof(ulong));
    mpoly_monomial_sub_mp(T0, Aexp + N*0, Bexp + N*0, N);
    mpoly_monomial_sub_mp(T1, Aexp + N*(Alen - 1), Bexp + N*(Blen - 1), N);
    if (bits <= FLINT_BITS)
    {
        if (   mpoly_monomial_overflows(T0, N, mask)
            || mpoly_monomial_overflows(T1, N, mask))
        {
flint_printf(" first division failed sp\n");
            failure = 1;
            goto cleanup;
        }
    }
    else
    {
        if (   mpoly_monomial_overflows_mp(T0, N, bits)
            || mpoly_monomial_overflows_mp(T1, N, bits))
        {
flint_printf(" first division failed mp\n");
            failure = 1;
            goto cleanup;
        }
    }

printf("here3\n");

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

flint_printf("C: "); nmod_mpoly_print_pretty(C, vars, H->ctx); flint_printf("\n");

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

void tryproc(Lholder_t H, Lchunk_t L, nmod_mpoly_t Q, const nmod_mpoly_t B)
{
    slong i;
    slong N = H->N;
    nmod_mpoly_struct * C = L->polyC;
    ulong mask;

    mask = 0;
    for (i = 0; i < FLINT_BITS/H->bits; i++)
        mask = (mask << H->bits) + (UWORD(1) << (H->bits - 1));


/*
flint_printf("tryproducer called\n");
*/
    if (L->Cinited == 0)
    {
        nmod_mpoly_init(C, H->ctx);
        nmod_mpoly_chop(C, H->polyA, L->upperclosed, L->emax, L->emin, H);
        L->Cinited = 1;
    }

    FLINT_ASSERT(L->mq >= 0);

    for (i = L->mq; i < Q->length; i++)
    {
        nmod_mpoly_t T;
        nmod_mpoly_init2(T, 16, H->ctx);
        nmod_mpoly_fit_bits(T, H->bits, H->ctx);
        T->bits = H->bits;
        nmod_mpoly_chopmulmonomial(T, B, Q, i, L, H);
        nmod_mpoly_sub(C, C, T, H->ctx);
        nmod_mpoly_clear(T, H->ctx);
    }
flint_printf("all C: "); nmod_mpoly_print_pretty(C, vars, H->ctx); flint_printf("\n");
flint_printf("all Q: "); nmod_mpoly_print_pretty(Q, vars, H->ctx); flint_printf("\n");
    L->mq = Q->length;


    if (L->producer == 1)
    {
        mp_limb_t * Qcoeff = Q->coeffs;
        ulong * Qexp = Q->exps;
        slong Qalloc =  Q->alloc;
        slong Qlen = Q->length;
        mp_limb_t lc_inv = nmod_inv(B->coeffs[0], H->ctx->ffinfo->mod);
        int divides;
        while (C->length > 0)
        {
            _nmod_mpoly_fit_length(&Qcoeff, &Qexp, &Qalloc, Qlen + 1, N);
            Qcoeff[Qlen] = nmod_mul(lc_inv, C->coeffs[0], H->ctx->ffinfo->mod);
            divides = mpoly_monomial_divides(Qexp + N*Qlen, C->exps + N*0, B->exps + N*0, N, mask);
            FLINT_ASSERT(divides);
            Qlen++;
            Q->exps = Qexp;
            Q->coeffs = Qcoeff;
            Q->alloc = Qalloc;
            Q->length = Qlen;           

            {
                nmod_mpoly_t T;
                nmod_mpoly_init2(T, 16, H->ctx);
                nmod_mpoly_fit_bits(T, H->bits, H->ctx);
                T->bits = H->bits;
                nmod_mpoly_chopmulmonomial(T, B, Q, Qlen - 1, L, H);
                nmod_mpoly_sub(C, C, T, H->ctx);
                nmod_mpoly_clear(T, H->ctx);
            }
flint_printf("pro C: "); nmod_mpoly_print_pretty(L->polyC, vars, H->ctx); flint_printf("\n");
flint_printf("pro Q: "); nmod_mpoly_print_pretty(Q, vars, H->ctx); flint_printf("\n");
usleep(100000);

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
}

/* return 1 if quotient is exact */
int nmod_mpoly_divides_heap_threaded(nmod_mpoly_t Q,
                  const nmod_mpoly_t A, const nmod_mpoly_t B,
                                                    const nmod_mpoly_ctx_t ctx)
{
    int ret;
    fmpz_mpoly_ctx_t zctx;
    fmpz_mpoly_t S;
    slong i,j,k,N;
    mp_bitcnt_t exp_bits;
    ulong * cmpmask;
    ulong * Aexp, * Bexp;
    int freeAexp = 0, freeBexp = 0;
    TMP_INIT;

printf("******************************\n");
printf("div A: "); nmod_mpoly_print_pretty(A, vars, ctx); printf("\n");
printf("div B: "); nmod_mpoly_print_pretty(B, vars, ctx); printf("\n");

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

    if (select_exps(S, zctx, Aexp, A->length, Bexp, B->length, exp_bits))
    {
        ret = 0;
        nmod_mpoly_zero(Q, ctx);
        goto cleanup1;
    }

flint_printf("S: "); fmpz_mpoly_print_pretty(S, vars, zctx); flint_printf("\n");

    {
        ulong * texps;
        mp_limb_t lc_inv;
        Lholder_t H;
        flint_rand_t randstate;
        
        Lholder_init(H);
        H->polyA = A;
        H->ctx = ctx;
        H->bits = exp_bits;
        H->N = N;
        H->cmpmask = cmpmask;

        for (i = 0; i + 1 < S->length; i++)
        {
            Lchunk_struct * L;

            L = (Lchunk_struct *) malloc(sizeof(Lchunk_struct));
            L->ma = 0;
            L->mq = 0;
            L->emax = S->exps + N*i;
            L->emin = S->exps + N*(i + 1);
            L->upperclosed = 0;
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

        mpoly_monomial_sub_mp(Q->exps + N*0, A->exps + N*0, B->exps + N*0, N);
        lc_inv = nmod_inv(B->coeffs[0], ctx->ffinfo->mod);
        Q->coeffs[0] = nmod_mul(lc_inv, A->coeffs[0], ctx->ffinfo->mod);
        Q->length = (Q->coeffs[0] != UWORD(0));
printf("1 Q: "); nmod_mpoly_print_pretty(Q, vars, ctx); printf("\n");

        mpoly_monomial_add_mp(texps, Q->exps + N*0, B->exps + N*1, N);

        k = 1;
        while (k < A->length && mpoly_monomial_gt(texps, A->exps + N*k, N, cmpmask))
        {
            nmod_mpoly_fit_length(Q, k + 1, ctx);
            mpoly_monomial_sub_mp(Q->exps + N*k, A->exps + N*k, B->exps + N*0, N);
            Q->coeffs[k] = nmod_mul(lc_inv, A->coeffs[0], ctx->ffinfo->mod);
            Q->length += (Q->coeffs[k] != UWORD(0));
        }
printf("2 Q: "); nmod_mpoly_print_pretty(Q, vars, ctx); printf("\n");

        flint_randinit(randstate);

        j = 20;
        while (--j > 0 && H->cur != NULL)
        {
            slong i;
            Lchunk_struct * L = H->cur;

            /* pick random */
            i = n_randint(randstate, H->length);

flint_printf("**** i: %wd *********************\n", i);

            L = H->cur;
            while (i > 0)
            {
                L = L->next;
                i--;
            }
            tryproc(H, L, Q, B);
        }

        flint_randclear(randstate);

        Lholder_clear(H);

    }

    ret = 1;

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
