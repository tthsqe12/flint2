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
void usleep(ulong);


typedef struct {
    mp_limb_t gammapow;
    ulong cm;
} nmod_dlog_table_entry_struct;

static int nmod_dlog_table_entry_struct_cmp(const nmod_dlog_table_entry_struct * lhs,
                                       const nmod_dlog_table_entry_struct * rhs)
{
    return (lhs->gammapow < rhs->gammapow) ? -1 : (lhs->gammapow > rhs->gammapow);
}

typedef struct {
    slong exp;
    ulong prime;
    mp_limb_t gamma;
    mp_limb_t gammainv;
    mp_limb_t startingbeta;
    ulong co;
    ulong startinge;
    ulong idem;
    ulong cbound;
    ulong dbound;
    nmod_dlog_table_entry_struct * table; /* length cbound */
} nmod_dlog_entry_struct;

typedef struct {
    nmod_t mod;         /* p is mod.n */
    mp_limb_t alpha;    /* p.r. of p */
    mp_limb_t alphainv;
    slong num_factors;  /* factors of p - 1*/
    nmod_dlog_entry_struct * entries;
} nmod_dlog_env_struct;
typedef nmod_dlog_env_struct nmod_dlog_env_t[1];

/* assume that p is prime, don't check */
void nmod_dlog_env_init(nmod_dlog_env_t L, mp_limb_t p)
{
    slong i;
    ulong c;
    nmod_dlog_entry_struct * Li;
    n_factor_t factors;

    n_factor_init(&factors);
    n_factor(&factors, p - 1, 1);

    nmod_init(&L->mod, p);
    L->entries = NULL;
    L->num_factors = factors.num;
    if (L->num_factors > 0)
    {
        L->entries = (nmod_dlog_entry_struct*) flint_malloc(L->num_factors*sizeof(nmod_dlog_entry_struct));
    }

    for (i = 0; i < L->num_factors; i++)
    {
        int success;
        fmpz_t pipow, pm1, temp, recp;

        Li = L->entries + i;

        Li->exp = factors.exp[i];
        Li->prime = factors.p[i];

        fmpz_init(recp);
        fmpz_init(temp);
        fmpz_init_set_ui(pipow, Li->prime);
        fmpz_pow_ui(pipow, pipow, Li->exp);
        fmpz_init_set_ui(pm1, p - 1);
        fmpz_divexact(recp, pm1, pipow);
        success = fmpz_invmod(temp, recp, pipow);
        FLINT_ASSERT(success);
        fmpz_mul(temp, temp, recp);

        Li->idem = fmpz_fdiv_ui(temp, p - 1);

        Li->co = fmpz_get_ui(recp);
        Li->startinge = fmpz_get_ui(pipow)/Li->prime;


        fmpz_clear(pipow);
        fmpz_clear(pm1);
        fmpz_clear(temp);
        fmpz_clear(recp);
    }

    L->alpha = 0;
try_alpha:
    L->alpha++;
    if (L->alpha >= p)
    {
        flint_printf("Exception in nmod_dlog_env_init: could not find primitive root\n");
        flint_abort();
    }
    for (i = 0; i < L->num_factors; i++)
    {
        Li = L->entries + i;
        Li->gamma = nmod_pow_ui(L->alpha, (p - 1) / Li->prime, L->mod);
        if (Li->gamma == 1)
        {
            goto try_alpha;
        }
    }

    L->alphainv = nmod_inv(L->alpha, L->mod);

    for (i = 0; i < L->num_factors; i++)
    {
        Li = L->entries + i;
        Li->gammainv = nmod_inv(Li->gamma, L->mod);

        Li->startingbeta = nmod_pow_ui(L->alphainv, Li->co, L->mod);

        Li->dbound = ceil(sqrt((double) Li->prime));
        Li->cbound = (Li->prime + Li->dbound - 1)/Li->dbound;
        while (Li->cbound > 100)
        {
            Li->dbound *= 2;
            Li->cbound = (Li->prime + Li->dbound - 1)/Li->dbound;
        }

        FLINT_ASSERT(Li->dbound > 0);
        FLINT_ASSERT(Li->cbound > 0);
        Li->table = (nmod_dlog_table_entry_struct *) flint_malloc(Li->cbound*sizeof(nmod_dlog_table_entry_struct));

        for (c = 0; c < Li->cbound; c++)
        {
            Li->table[c].cm = c*Li->dbound;
            Li->table[c].gammapow = nmod_pow_ui(Li->gamma, Li->table[c].cm, L->mod);
        }
        qsort(Li->table, Li->cbound, sizeof(nmod_dlog_table_entry_struct),
               (int(*)(const void*, const void*)) nmod_dlog_table_entry_struct_cmp);
        for (c = 1; c < Li->cbound; c++)
        {
            FLINT_ASSERT(Li->table[c - 1].gammapow < Li->table[c].gammapow);
            FLINT_ASSERT(Li->table[c].gammapow == nmod_pow_ui(Li->gamma, Li->table[c].cm, L->mod));
        }
    }    
}

void nmod_dlog_env_clear(nmod_dlog_env_t L)
{
    slong i;
    nmod_dlog_entry_struct * Li;

    for (i = 0; i < L->num_factors; i++)
    {
        Li = L->entries + i;
        flint_free(Li->table);
    }

    if (L->entries)
    {
        flint_free(L->entries);
    }
}

/* return x such that x = alpha^y mod p, alpha is the p.r. L->alpha*/
ulong nmod_dlog_env_run(const nmod_dlog_env_t L, mp_limb_t y)
{
    slong i, j;
    ulong x, q, r, e, x0 = 0, x1 = 0, x2 = 0, pp0, pp1, acc, g, pipow;
    ulong lo, mid, hi, d;
    mp_limb_t beta, z, w;
    nmod_dlog_entry_struct * Li;

    FLINT_ASSERT(y != 0);
    FLINT_ASSERT(y < L->mod.n);

    i = 0;
    if (i < L->num_factors && L->entries[i].prime == 2)
    {
        Li = L->entries + i;
        FLINT_ASSERT(Li->prime == 2);

        z = nmod_pow_ui(y, Li->co, L->mod);
        beta = Li->startingbeta;
        e = Li->startinge;
        j = 0;
        pipow = 1; /* Li->prime^j */
        acc = 0;
        do {
            w = nmod_pow_ui(z, e, L->mod);
            /* solve Li->gamma ^ g == w mod p */

            if (w == 1)
            {
                g = 0;
            }
            else
            {
                FLINT_ASSERT(w == Li->gamma);
                g = 1;
                z = nmod_mul(z, beta, L->mod);
            }
            beta = nmod_mul(beta, beta, L->mod);
            acc += g*pipow;
            pipow = pipow*2;
            e = e / 2;
        } while (++j < Li->exp);

        umul_ppmm(pp1, pp0, acc, Li->idem);
        add_sssaaaaaa(x2, x1, x0, x2, x1, x0, WORD(0), pp1, pp0);
        i = 1;
    }

    for (; i < L->num_factors; i++)
    {
        Li = L->entries + i;

        z = nmod_pow_ui(y, Li->co, L->mod);
        beta = Li->startingbeta;
        e = Li->startinge;
        j = 0;
        pipow = 1; /* Li->prime^j */
        acc = 0;
        do {
            w = nmod_pow_ui(z, e, L->mod);
            /* solve Li->gamma ^ g == w mod p */
            d = 0;
            while (1)
            {
                lo = 0; hi = Li->cbound;
                while (hi - lo > 4)
                {
                    mid = lo + (hi - lo)/2;
                    if (Li->table[mid].gammapow == w)
                    {
                        g = Li->table[mid].cm + d;
                        goto found_g;
                    }
                    if (Li->table[mid].gammapow > w)
                        hi = mid;
                    else
                        lo = mid;
                }
                while (lo < hi)
                {
                    if (Li->table[lo].gammapow == w)
                    {
                        g = Li->table[lo].cm + d;
                        goto found_g;
                    }
                    lo++;
                }
                w = nmod_mul(w, Li->gammainv, L->mod);
                d++;
                /* should have found a solution if d is out of bounds */
                FLINT_ASSERT(d < Li->dbound);
            }
        found_g:
            FLINT_ASSERT(g < Li->prime);
            z = nmod_mul(z, nmod_pow_ui(beta, g, L->mod), L->mod);
            beta = nmod_pow_ui(beta, Li->prime, L->mod);
            acc += g*pipow;
            pipow = pipow*Li->prime;
            e = e / Li->prime;
        } while (++j < Li->exp);

        umul_ppmm(pp1, pp0, acc, Li->idem);
        add_sssaaaaaa(x2, x1, x0, x2, x1, x0, WORD(0), pp1, pp0);
    }

    udiv_qrnnd(q, r, x2, x1, L->mod.n - 1);
    udiv_qrnnd(q, x, r, x0, L->mod.n - 1);
    return x;
}



typedef struct {
    slong half_point_count;
    nmod_poly_t R0, R1;
    nmod_poly_t V0, V1; /* V1 is our master polynomial */
    nmod_poly_t r; /* temporary */
    nmod_poly_t q; /* temporary also used as a queue of incoming points */
} nmod_bma_struct;
typedef nmod_bma_struct nmod_bma_t[1];

/*
    A = x^(2*n)       n = half_point_count
    deg(S) < 2*n
    U0*A + V0*S = R0   deg(R0) >= n
    U1*A + V1*S = R1   deg(R1) < n

    S is the reverse of the polynomial whose coefficients are the input points.
        S = a_0*x^(2n-1) + a_1*x^(2n-2) + ... + a_(2n-1)
    S can be updated with more points at any time via add_points.
    S can be reduced (compute V0, V1, R0, R1) at any time via reduce, which
        returns whether the reduction caused a change.

    The U0 and U1 are not stored.
*/
void nmod_bma_init(nmod_bma_t B, mp_limb_t p)
{
    B->half_point_count = 0;
    nmod_poly_init(B->V0, p);
    nmod_poly_init(B->R0, p);
    nmod_poly_one(B->R0);
    nmod_poly_init(B->V1, p);
    nmod_poly_one(B->V1);
    nmod_poly_init(B->R1, p);
    nmod_poly_init(B->r, p);
    nmod_poly_init(B->q, p);
    B->q->length = 0;
}

void nmod_bma_start_over(nmod_bma_t B)
{
    B->half_point_count = 0;
    nmod_poly_zero(B->V0);
    nmod_poly_one(B->R0);
    nmod_poly_one(B->V1);
    nmod_poly_zero(B->R1);
    B->q->length = 0;
}

void nmod_bma_clear(nmod_bma_t B)
{
    nmod_poly_clear(B->R0);
    nmod_poly_clear(B->R1);
    nmod_poly_clear(B->V0);
    nmod_poly_clear(B->V1);
    nmod_poly_clear(B->r);
    nmod_poly_clear(B->q);
}

void nmod_bma_add_points(nmod_bma_t B, mp_limb_t * a, slong count)
{
    slong i;
    slong old_length = B->q->length;
    nmod_poly_fit_length(B->q, old_length + count);
    for (i = 0; i < count; i++)
    {
        B->q->coeffs[old_length + i] = a[i];
    }
    B->q->length = old_length + count;
}

void nmod_bma_add_point(nmod_bma_t B, mp_limb_t a)
{
    slong old_length = B->q->length;
    nmod_poly_fit_length(B->q, old_length + 1);
    B->q->coeffs[old_length] = a;
    B->q->length = old_length + 1;
}

int nmod_bma_reduce(nmod_bma_t B)
{
    int changed = 0;
    slong i, queue_length = B->q->length;

    if ((queue_length % 2) != 0)
    {
        flint_printf("Exception in nmod_bma_reduce: point count is not even\n");
        flint_abort();
    }

    /* reverse the queue into temp r */
    B->half_point_count += queue_length/2;
    nmod_poly_zero(B->r);
    for (i = 0; i < queue_length; i++)
    {
        nmod_poly_set_coeff_ui(B->r, queue_length - i - 1, B->q->coeffs[i]);
    }
    nmod_poly_mul(B->q, B->V0, B->r);
    nmod_poly_shift_left(B->R0, B->R0, queue_length);
    nmod_poly_add(B->R0, B->R0, B->q);
    nmod_poly_mul(B->q, B->V1, B->r);
    nmod_poly_shift_left(B->R1, B->R1, queue_length);
    nmod_poly_add(B->R1, B->R1, B->q);

    /* now reduce */
    while (B->half_point_count < nmod_poly_length(B->R1))
    {
        changed = 1;
        nmod_poly_divrem(B->q, B->r, B->R0, B->R1);
        nmod_poly_swap(B->R0, B->R1);
        nmod_poly_swap(B->R1, B->r);

        nmod_poly_mul(B->r, B->q, B->V1);
        nmod_poly_sub(B->q, B->V0, B->r);
        nmod_poly_swap(B->V0, B->V1);
        nmod_poly_swap(B->V1, B->q);
        FLINT_ASSERT(nmod_poly_degree(B->V1) > nmod_poly_degree(B->V0));
    }

    /* queue is empty now */
    B->q->length = 0;
    return changed;
}


/* split f assuming that f has degree(f) distinct nonzero roots in Fp */
static void _nmod_poly_rabinsplit(nmod_poly_t a, nmod_poly_t b, nmod_poly_t T, 
                                  const nmod_poly_t f, flint_rand_t randstate)
{
    mp_limb_t delta;

    FLINT_ASSERT(nmod_poly_degree(f) > 1);

try_again:

    delta = n_randint(randstate, f->mod.n);

    nmod_poly_zero(a);
    nmod_poly_set_coeff_ui(a, 1, 1);
    nmod_poly_set_coeff_ui(a, 0, delta);
    nmod_poly_powmod_ui_binexp(T, a, (f->mod.n - 1)/2, f);
    nmod_poly_zero(a);
    nmod_poly_set_coeff_ui(a, 0, 1);
    nmod_poly_sub(T, T, a);
    nmod_poly_gcd(a, T, f);
    FLINT_ASSERT(!nmod_poly_is_zero(a));
    if (0 >= nmod_poly_degree(a) || nmod_poly_degree(a) >= nmod_poly_degree(f))
    {
        goto try_again;
    }
    nmod_poly_div(b, f, a);
    /* deg a >= deg b */
    if (nmod_poly_degree(a) < nmod_poly_degree(b))
    {
        nmod_poly_swap(a, b);
    }
    return;
}

/* fill in roots with the t distinct nonzero roots of master, or fail */
static int _nmod_find_roots(mp_limb_t * roots, const nmod_poly_t master, slong t)
{
    mp_limb_t a0, a1;
    int success;
    slong i, roots_idx, sp;
    mp_limb_t delta;
    nmod_poly_struct * a , * b;
    nmod_poly_t f, T;
    nmod_poly_struct stack[FLINT_BITS + 1];
    flint_rand_t randstate;

    FLINT_ASSERT(t >= 0);
    FLINT_ASSERT(t == nmod_poly_degree(master));

    if (t == 0)
    {
        return 1;
    }
    else if (t == 1)
    {
        a0 = nmod_poly_get_coeff_ui(master, 0);
        a1 = nmod_poly_get_coeff_ui(master, 1);
        if (a0 == 0)
        {
            return 0;
        }
        roots[0] = nmod_mul(a0, nmod_inv(master->mod.n - a1, master->mod), master->mod);
        return 1;
    }

    flint_randinit(randstate);
    nmod_poly_init(T, master->mod.n);
    nmod_poly_init(f, master->mod.n);
    for (i = 0; i <= FLINT_BITS; i++)
    {
        nmod_poly_init(stack + i, master->mod.n);
    }

    roots_idx = 0;

    nmod_poly_make_monic(f, master);

    a = stack + 0;
    delta = n_randint(randstate, master->mod.n);
    nmod_poly_zero(a);
    nmod_poly_set_coeff_ui(a, 1, 1);
    nmod_poly_set_coeff_ui(a, 0, delta);
    nmod_poly_powmod_ui_binexp(T, a, (master->mod.n - 1)/2, f);
    nmod_poly_zero(a);
    nmod_poly_set_coeff_ui(a, 0, 1);
    nmod_poly_sub(T, T, a);
    nmod_poly_gcd(a, T, f);

    b = stack + 1;
    nmod_poly_zero(b);
    nmod_poly_set_coeff_ui(b, 0, 2);
    nmod_poly_add(T, T, b);
    nmod_poly_gcd(b, T, f);

    if (nmod_poly_degree(b) + nmod_poly_degree(a) != t)
    {
        success = 0;
        goto cleanup;
    }
    /* deg a >= deg b */
    if (nmod_poly_degree(a) < nmod_poly_degree(b))
    {
        nmod_poly_swap(a, b);
    }

    sp = nmod_poly_degree(b) > 0 ? 2 : 1;
    while (sp > 0)
    {
        FLINT_ASSERT(sp < FLINT_BITS);
        sp--;
        nmod_poly_swap(f, stack + sp);

        FLINT_ASSERT(nmod_poly_degree(f) > 0);
        if (nmod_poly_degree(f) == 1)
        {
            a0 = nmod_poly_get_coeff_ui(f, 0);
            a1 = nmod_poly_get_coeff_ui(f, 1);
            FLINT_ASSERT(a0 != 0);
            FLINT_ASSERT(a1 == 1);
            roots[roots_idx] = master->mod.n - a0;
            roots_idx++;
        }
        else
        {
            _nmod_poly_rabinsplit(stack + sp + 0, stack + sp + 1, T, f, randstate);
            FLINT_ASSERT(FLINT_BIT_COUNT(nmod_poly_degree(stack + sp + 1)) <= FLINT_BITS - sp - 1);
            sp += 2;
        }
    }

    success = 1;

cleanup:

    flint_randclear(randstate);
    nmod_poly_clear(T);
    nmod_poly_clear(f);
    for (i = 0; i <= FLINT_BITS; i++)
    {
        nmod_poly_clear(stack + i);
    }

    if (success)
    {
        FLINT_ASSERT(roots_idx == t);
    }

    return success;
}


typedef struct
{
    const nmod_mpoly_ctx_struct * polyctx;
    nmod_dlog_env_t dlogenv;
    slong * degbounds;
    ulong * subdegs;
    mp_bitcnt_t bits;
    ulong * inputexpmask;
} nmod_mpoly_bma_interpolate_ctx_struct;
typedef nmod_mpoly_bma_interpolate_ctx_struct nmod_mpoly_bma_interpolate_ctx_t[1];

typedef struct {
    nmod_poly_t roots;
    nmod_poly_t evals;
    nmod_bma_t bma;
} nmod_mpoly_bma_interpolate_struct;
typedef nmod_mpoly_bma_interpolate_struct nmod_mpoly_bma_interpolate_t[1];



void nmod_mpoly_bma_interpolate_ctx_init(nmod_mpoly_bma_interpolate_ctx_t Ictx, mp_bitcnt_t bits_,
                                                const nmod_mpoly_ctx_t pctx)
{
    slong N;

    Ictx->polyctx = pctx;
    Ictx->degbounds = (slong *) flint_malloc(pctx->minfo->nvars*sizeof(slong));
    Ictx->subdegs = (ulong *) flint_malloc(pctx->minfo->nvars*sizeof(ulong));
    nmod_dlog_env_init(Ictx->dlogenv, pctx->ffinfo->mod.n);

    Ictx->bits = bits_;
    N = mpoly_words_per_exp_sp(bits_, pctx->minfo);
    Ictx->inputexpmask = (ulong *) flint_malloc(N*sizeof(ulong));
    mpoly_monomial_zero(Ictx->inputexpmask, N);
}

void nmod_mpoly_bma_interpolate_ctx_clear(nmod_mpoly_bma_interpolate_ctx_t Ictx)
{
    flint_free(Ictx->inputexpmask);
    nmod_dlog_env_clear(Ictx->dlogenv);
    flint_free(Ictx->degbounds);
    flint_free(Ictx->subdegs);
}


void nmod_mpoly_bma_interpolate_init(nmod_mpoly_bma_interpolate_t I, mp_limb_t p)
{
    nmod_poly_init(I->roots, p);
    nmod_poly_init(I->evals, p);
    nmod_bma_init(I->bma, p);
}

void nmod_mpoly_bma_interpolate_clear(nmod_mpoly_bma_interpolate_t I)
{
    nmod_poly_clear(I->roots);
    nmod_poly_clear(I->evals);
    nmod_bma_clear(I->bma);
}

void nmod_mpoly_bma_interpolate_add_point(nmod_mpoly_bma_interpolate_t I, mp_limb_t a)
{
    nmod_poly_fit_length(I->evals, I->evals->length + 1);
    I->evals->coeffs[I->evals->length] = a;
    I->evals->length++;

    nmod_bma_add_point(I->bma, a);
}

int nmod_mpoly_bma_interpolate_reduce(nmod_mpoly_bma_interpolate_t I)
{
    return nmod_bma_reduce(I->bma);
}

void nmod_mpoly_bma_interpolate_eval_init(nmod_mpoly_bma_interpolate_ctx_t Ictx, const nmod_mpoly_t A)
{    
    slong i, j, N = mpoly_words_per_exp_sp(Ictx->bits, Ictx->polyctx->minfo);
    FLINT_ASSERT(A->bits == Ictx->bits);

    for (i = 0; i < A->length; i++)
    {
        for (j = 0; j < N; j++)
        {
            Ictx->inputexpmask[j] |= (A->exps + N*i)[j];
        }
    }
}



/*
    put the evaluation of the monomials in A at alpha^w in the coeffs of E

    x_0     = alpha ^ (w * subdegs[n-1] * subdegs[n-2] * ... * * subdegs[1])
      ...
    x_(n-3) = alpha ^ (w * subdegs[n-1] * subdegs[n-2])
    x_(n-2) = alpha ^ (w * subdegs[n-1])
    x_(n-1) = alpha ^ (w)

    secret: subdegs[0] is not used
*/
void nmod_mpoly_bma_interpolate_eval_setskel(nmod_mpoly_t M, const nmod_mpoly_t A, ulong w, const nmod_mpoly_bma_interpolate_ctx_t Ictx)
{
    slong i, j;
    slong offset, shift;
    slong N = mpoly_words_per_exp_sp(A->bits, Ictx->polyctx->minfo);
    slong nvars = Ictx->polyctx->minfo->nvars;
    ulong * Aexp;
    mp_limb_t * Mcoeff;
    slong * LUToffset;
    ulong * LUTmask;
    mp_limb_t * LUTvalue;
    slong LUTlen;
    mp_limb_t xeval, xpoweval;
    TMP_INIT;

    TMP_START;

    LUToffset = (slong *) TMP_ALLOC(N*FLINT_BITS*sizeof(slong));
    LUTmask   = (ulong *) TMP_ALLOC(N*FLINT_BITS*sizeof(ulong));
    LUTvalue  = (mp_limb_t *) TMP_ALLOC(N*FLINT_BITS*sizeof(mp_limb_t));

    FLINT_ASSERT(M->bits == A->bits);

    nmod_mpoly_fit_length(M, A->length, Ictx->polyctx);
    M->length = A->length;

    Mcoeff = M->coeffs;
    Aexp = A->exps;

    LUTlen = 0;
    xeval = nmod_pow_ui(Ictx->dlogenv->alpha, w, Ictx->polyctx->ffinfo->mod);
    for (j = nvars - 1; j >= 0; j--)
    {
        mpoly_gen_offset_shift_sp(&offset, &shift, j, A->bits, Ictx->polyctx->minfo);

        xpoweval = xeval; /* xpoweval = xeval^(2^i) */
        for (i = 0; i < A->bits; i++)
        {
            LUToffset[LUTlen] = offset;
            LUTmask[LUTlen] = (UWORD(1) << (shift + i));
            LUTvalue[LUTlen] = xpoweval;
            if ((Ictx->inputexpmask[offset] & LUTmask[LUTlen]) != 0)
            {
                LUTlen++;
            }
            xpoweval = nmod_mul(xpoweval, xpoweval, Ictx->polyctx->ffinfo->mod);
        }
        xeval = nmod_pow_ui(xeval, Ictx->subdegs[j], Ictx->polyctx->ffinfo->mod);
    }
    FLINT_ASSERT(LUTlen < N*FLINT_BITS);

    for (i = 0; i < A->length; i++)
    {
        mp_limb_t t = 1;
        for (j = 0; j < LUTlen; j++)
        {
            if (((Aexp + N*i)[LUToffset[j]] & LUTmask[j]) != 0)
            {
                t = nmod_mul(t, LUTvalue[j], Ictx->polyctx->ffinfo->mod);
            }
        }
        Mcoeff[i] = t;
        mpoly_monomial_zero(M->exps + N*i, N);
    }

    TMP_END;
}

/* M = S */
void nmod_mpoly_bma_interpolate_eval_copyskel(nmod_mpoly_t M, const nmod_mpoly_t S, const nmod_mpoly_bma_interpolate_ctx_t Ictx)
{
    slong i, N;
    nmod_mpoly_fit_length(M, S->length, Ictx->polyctx);
    M->length = S->length;
    _nmod_vec_set(M->coeffs, S->coeffs, S->length);
    N = mpoly_words_per_exp(Ictx->bits, Ictx->polyctx->minfo);
    for (i = 0; i < M->length; i++)
    {
        mpoly_monomial_zero(M->exps + N*i, N);
    }
}

/* return A.M */
mp_limb_t nmod_mpoly_bma_interpolate_eval_useskel(const nmod_mpoly_t A, const nmod_mpoly_t M, const nmod_mpoly_bma_interpolate_ctx_t Ictx)
{
    slong i;
    mp_limb_t t = 0;

    FLINT_ASSERT(M->length == A->length);
    for (i = 0; i < A->length; i++)
    {
        t = nmod_add(t, nmod_mul(A->coeffs[i], M->coeffs[i], Ictx->polyctx->ffinfo->mod), Ictx->polyctx->ffinfo->mod);
    }
    return t;
}

/* return A.M and multiply M by S */
mp_limb_t nmod_mpoly_bma_interpolate_eval_useskelmul(const nmod_mpoly_t A, nmod_mpoly_t M, const nmod_mpoly_t S, const nmod_mpoly_bma_interpolate_ctx_t Ictx)
{
    slong i;
    mp_limb_t t = 0;

    for (i = 0; i < A->length; i++)
    {
        t = nmod_add(t, nmod_mul(A->coeffs[i], M->coeffs[i], Ictx->polyctx->ffinfo->mod), Ictx->polyctx->ffinfo->mod);
        M->coeffs[i] = nmod_mul(M->coeffs[i], S->coeffs[i], Ictx->polyctx->ffinfo->mod);
    }
    return t;
}



int nmod_mpoly_bma_interpolate_get_mpoly(nmod_mpoly_t A, ulong alphashift, nmod_mpoly_bma_interpolate_t I,
                                        const nmod_mpoly_bma_interpolate_ctx_t Ictx)
{
    slong i, j, t, N;
    int success;
    ulong new_exp, this_exp;
    slong * shifts, * offsets;
    mp_limb_t * values, * roots;
    mp_limb_t * Acoeff;
    ulong * Aexp;
    slong Alen;
    mp_limb_t T, S, V, V0, V1, V2, p0, p1, r;
    TMP_INIT;

    TMP_START;

    nmod_bma_reduce(I->bma);
    t = nmod_poly_degree(I->bma->V1);
    FLINT_ASSERT(t > 0);
    FLINT_ASSERT(I->evals->length >= t);

    nmod_poly_fit_length(I->roots, t);
    I->roots->length = t;
    success = _nmod_find_roots(I->roots->coeffs, I->bma->V1, t);
    if (!success)
    {
        goto cleanup;
    }

    roots = I->roots->coeffs;
    values = I->evals->coeffs;

    FLINT_ASSERT(Ictx->polyctx->minfo->ord == ORD_LEX);
    FLINT_ASSERT(A->bits <= FLINT_BITS);
    N = mpoly_words_per_exp_sp(A->bits, Ictx->polyctx->minfo);
    nmod_mpoly_fit_length(A, t, Ictx->polyctx);
    Acoeff = A->coeffs;
    Aexp = A->exps;
    A->length = 0;

    shifts = (slong *) TMP_ALLOC(Ictx->polyctx->minfo->nvars);
    offsets = (slong *) TMP_ALLOC(Ictx->polyctx->minfo->nvars);
    for (j = 0; j < Ictx->polyctx->minfo->nvars; j++)
    {
        mpoly_gen_offset_shift_sp(offsets + j, shifts + j, j, A->bits, Ictx->polyctx->minfo);
    }

    Alen = 0;
    for (i = 0; i < t; i++)
    {
        /* coeffs[i] is (coeffs(P).values)/P(roots[i]) =: V/S
            where P(x) = master(x)/(x-roots[i])     */
        V0 = V1 = V2 = T = S = 0;
        r = roots[i];
        for (j = t; j > 0; j--)
        {
            T = nmod_add(nmod_mul(r, T, I->bma->V1->mod), I->bma->V1->coeffs[j], I->bma->V1->mod);
            S = nmod_add(nmod_mul(r, S, I->bma->V1->mod), T, I->bma->V1->mod);
            umul_ppmm(p1, p0, values[j - 1], T);
            add_sssaaaaaa(V2, V1, V0, V2, V1, V0, WORD(0), p1, p0);
        }
        /* roots[i] should be a root of master */
        FLINT_ASSERT(nmod_add(nmod_mul(r, T, I->bma->V1->mod), I->bma->V1->coeffs[0], I->bma->V1->mod) == 0);
        NMOD_RED3(V, V2, V1, V0, I->bma->V1->mod);
        S = nmod_mul(S, nmod_pow_ui(r, alphashift, I->bma->V1->mod), I->bma->V1->mod);
        Acoeff[Alen] = nmod_mul(V, nmod_inv(S, I->bma->V1->mod), I->bma->V1->mod);
        if (Acoeff[Alen] == 0)
        {
            /* hmmm */
            continue;
        }

        mpoly_monomial_zero(Aexp + N*Alen, N);
        new_exp = nmod_dlog_env_run(Ictx->dlogenv, roots[i]);
        for (j = Ictx->polyctx->minfo->nvars - 1; j >= 0; j--)
        {
            this_exp = new_exp % Ictx->subdegs[j];
            new_exp = new_exp / Ictx->subdegs[j];
            if (this_exp >= Ictx->degbounds[j])
            {
                success = 0;
                goto cleanup;
            }
            (Aexp + N*Alen)[offsets[j]] |= this_exp << shifts[j];
        }
        if (new_exp != 0)
        {
            success = 0;
            goto cleanup;
        }
        Alen++;
    }
    A->length = Alen;

    success = 1;

cleanup:

    TMP_END;
    return success;
}




/***************************************************************

    start of fmpz !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

******************************************************************/

typedef struct {
    fmpz_t p;
    fmpz_preinvn_t pinv;
} fmpz_mod_ctx_struct;
typedef fmpz_mod_ctx_struct fmpz_mod_ctx_t[1];

void fmpz_mod_ctx_init(fmpz_mod_ctx_t ctx, const fmpz_t p)
{
    fmpz_init_set(ctx->p, p);
    fmpz_preinvn_init(ctx->pinv, p);
}

void fmpz_mod_ctx_clear(fmpz_mod_ctx_t ctx)
{
    fmpz_preinvn_clear(ctx->pinv);
    fmpz_clear(ctx->p);
}

int fmpz_mod_is_canonical(const fmpz_t a, const fmpz_mod_ctx_t ctx)
{
    return fmpz_sgn(a) >= 0 && fmpz_cmp(a, ctx->p) < 0;
}

void fmpz_mod_add(fmpz_t a, const fmpz_t b, const fmpz_t c, const fmpz_mod_ctx_t ctx)
{
    FLINT_ASSERT(fmpz_mod_is_canonical(b, ctx));
    FLINT_ASSERT(fmpz_mod_is_canonical(c, ctx));
    fmpz_add(a, b, c);
    if (fmpz_cmpabs(a, ctx->p) >= 0)
    {
        fmpz_sub(a, a, ctx->p);
    }
    FLINT_ASSERT(fmpz_mod_is_canonical(a, ctx));
}

void fmpz_mod_sub(fmpz_t a, const fmpz_t b, const fmpz_t c, const fmpz_mod_ctx_t ctx)
{
    FLINT_ASSERT(fmpz_mod_is_canonical(b, ctx));
    FLINT_ASSERT(fmpz_mod_is_canonical(c, ctx));
    fmpz_sub(a, b, c);
    if (fmpz_sgn(a) < 0)
    {
        fmpz_add(a, a, ctx->p);
    }
    FLINT_ASSERT(fmpz_mod_is_canonical(a, ctx));
}

void fmpz_mod_neg(fmpz_t a, const fmpz_t b, const fmpz_mod_ctx_t ctx)
{
    FLINT_ASSERT(fmpz_mod_is_canonical(b, ctx));
    fmpz_neg(a, b);
    if (fmpz_sgn(a) < 0)
    {
        fmpz_add(a, a, ctx->p);
    }
    FLINT_ASSERT(fmpz_mod_is_canonical(a, ctx));
}

void fmpz_mod_mul(fmpz_t a, const fmpz_t b, const fmpz_t c, const fmpz_mod_ctx_t ctx)
{
/*
    fmpz_t q, t;
    FLINT_ASSERT(fmpz_mod_is_canonical(b, ctx));
    FLINT_ASSERT(fmpz_mod_is_canonical(c, ctx));
    fmpz_init(q);
    fmpz_init(t);
    fmpz_mul(t, b, c);
    fmpz_fdiv_qr_preinvn(q, a, t, ctx->p, ctx->pinv);
    fmpz_clear(q);
    fmpz_clear(t);
    FLINT_ASSERT(fmpz_mod_is_canonical(a, ctx));
*/

    fmpz_t t;
    FLINT_ASSERT(fmpz_mod_is_canonical(b, ctx));
    FLINT_ASSERT(fmpz_mod_is_canonical(c, ctx));
    fmpz_init(t);
    fmpz_mul(t, b, c);
    fmpz_mod(a, t, ctx->p);
    fmpz_clear(t);
    FLINT_ASSERT(fmpz_mod_is_canonical(a, ctx));

}

void fmpz_mod_inv(fmpz_t a, const fmpz_t b, const fmpz_mod_ctx_t ctx)
{
    fmpz_t d;
    FLINT_ASSERT(fmpz_mod_is_canonical(b, ctx));
    fmpz_init(d);
    fmpz_gcdinv(d, a, b, ctx->p);
/*
printf("inv p: "); fmpz_print(ctx->p); printf("\n");
printf("inv b: "); fmpz_print(b); printf("\n");
printf("inv a: "); fmpz_print(a); printf("\n");
printf("inv d: "); fmpz_print(d); printf("\n");
*/
    if (!fmpz_is_one(d))
    {
        flint_throw(FLINT_IMPINV, "Cannot invert in fmpz_mod_inv\n");        
    }
    fmpz_clear(d);

    fmpz_mod(a, a, ctx->p);

    FLINT_ASSERT(fmpz_mod_is_canonical(a, ctx));
}

void fmpz_mod_pow_ui(fmpz_t a, const fmpz_t b, ulong pow, const fmpz_mod_ctx_t ctx)
{
    fmpz_t r, x;

    fmpz_init_set_ui(r, 1);
    fmpz_init_set(x, b);

    while (pow != 0)
    {
        if ((pow & 1) != 0)
        {
            fmpz_mod_mul(r, r, x, ctx);
        }
        pow = pow >> 1;
        if (pow != 0)
        {
            fmpz_mod_mul(x, x, x, ctx);
        }
    }

    fmpz_swap(a, r);
    fmpz_clear(r);
    fmpz_clear(x);
    FLINT_ASSERT(fmpz_mod_is_canonical(a, ctx));
}

void fmpz_mod_pow_fmpz(fmpz_t a, const fmpz_t b, const fmpz_t pow, const fmpz_mod_ctx_t ctx)
{
    mp_bitcnt_t i, bits;
    fmpz_t r, x;

    fmpz_init_set_ui(r, 1);
    if (fmpz_sgn(pow) < 0)
    {
        fmpz_init(x);
        fmpz_mod_inv(x, b, ctx);
    }
    else
    {
        fmpz_init_set(x, b);
    }

    bits = fmpz_bits(pow);

    for (i = 0; i < bits; i++)
    {
        if (fmpz_tstbit(pow, i) != 0)
        {
            fmpz_mod_mul(r, r, x, ctx);
        }
        fmpz_mod_mul(x, x, x, ctx);
    }

    fmpz_swap(a, r);

    fmpz_clear(r);
    fmpz_clear(x);
    FLINT_ASSERT(fmpz_mod_is_canonical(a, ctx));
}




typedef struct
{
   fmpz_mod_poly_struct * coeffs;
   ulong * exps;
   slong alloc;
   slong length;
   slong bits;
} fmpz_mod_mpolyn_struct;
typedef fmpz_mod_mpolyn_struct fmpz_mod_mpolyn_t[1];

typedef struct
{
    fmpz_mod_mpolyn_struct * coeffs;
    ulong * exps;
    slong alloc;
    slong length;
    mp_bitcnt_t bits;   /* default bits to construct coeffs */
} fmpz_mod_mpolyun_struct;
typedef fmpz_mod_mpolyun_struct fmpz_mod_mpolyun_t[1];



void fmpz_mod_mpolyn_init(
    fmpz_mod_mpolyn_t A,
    mp_bitcnt_t bits,
    const fmpz_mpoly_ctx_t ctx,
    const fmpz_mod_ctx_t fpctx)
{
    A->coeffs = NULL;
    A->exps = NULL;
    A->alloc = 0;
    A->length = 0;
    A->bits = bits;
}

void fmpz_mod_mpolyn_clear(
    fmpz_mod_mpolyn_t A,
    const fmpz_mpoly_ctx_t ctx,
    const fmpz_mod_ctx_t fpctx)
{
    slong i;
    for (i = 0; i < A->alloc; i++)
        fmpz_mod_poly_clear(A->coeffs + i);
    if (A->coeffs)
        flint_free(A->coeffs);
    if (A->exps)
        flint_free(A->exps);
}

void fmpz_mod_mpolyn_swap(fmpz_mod_mpolyn_t A, fmpz_mod_mpolyn_t B)
{
   fmpz_mod_mpolyn_struct t = *A;
   *A = *B;
   *B = t;
}

void fmpz_mod_mpolyn_zero(
    fmpz_mod_mpolyn_t A,
    const fmpz_mpoly_ctx_t ctx,
    const fmpz_mod_ctx_t fpctx)
{
    A->length = 0;
}

int fmpz_mod_mpolyn_is_zero(
    fmpz_mod_mpolyn_t A,
    const fmpz_mpoly_ctx_t ctx,
    const fmpz_mod_ctx_t fpctx)
{
    return A->length == 0;
}

void fmpz_mod_mpolyn_print_pretty(
    const fmpz_mod_mpolyn_t A,
    const char ** x_in,
    const fmpz_mpoly_ctx_t ctx,
    const fmpz_mod_ctx_t fpctx)
{
    fmpz_mod_poly_struct * coeff = A->coeffs;
    slong len = A->length;
    ulong * exp = A->exps;
    slong bits = A->bits;
    slong i, j, N;
    fmpz * exponents;
    char ** x = (char **) x_in;
    TMP_INIT;

    if (len == 0)
    {
        flint_printf("0");
        return;
    }

    N = mpoly_words_per_exp(bits, ctx->minfo);

    TMP_START;

    if (x == NULL)
    {
        x = (char **) TMP_ALLOC(ctx->minfo->nvars*sizeof(char *));
        for (i = 0; i < ctx->minfo->nvars; i++)
        {
            x[i] = (char *) TMP_ALLOC(((FLINT_BITS+4)/3)*sizeof(char));
            flint_sprintf(x[i], "x%wd", i+1);
        }
    }

    exponents = (fmpz *) TMP_ALLOC(ctx->minfo->nvars*sizeof(ulong));
    for (i = 0; i < ctx->minfo->nvars; i++)
        fmpz_init(exponents + i);
   
    for (i = 0; i < len; i++)
    {
        if (i > 0)
        {
            printf(" + ");
        }

        printf("(");
        fmpz_mod_poly_print_pretty(coeff + i, "v");
        printf(")");

        mpoly_get_monomial_ffmpz(exponents, exp + N*i, bits, ctx->minfo);

        for (j = 0; j < ctx->minfo->nvars; j++)
        {
            int cmp = fmpz_cmp_ui(exponents + j, WORD(1));

            if (cmp > 0)
            {
                printf("*%s^", x[j]);
                fmpz_print(exponents + j);
            } else if (cmp == 0)
            {
                printf("*%s", x[j]);
            }
        }
    }

    for (i = 0; i < ctx->minfo->nvars; i++)
        fmpz_clear(exponents + i);

    TMP_END;
}

void fmpz_mod_mpolyn_fit_length(
    fmpz_mod_mpolyn_t A,
    slong length,
    const fmpz_mpoly_ctx_t ctx,
    const fmpz_mod_ctx_t fpctx)
{
    slong i;
    slong old_alloc = A->alloc;
    slong new_alloc = FLINT_MAX(length, 2*A->alloc);

    if (length > old_alloc)
    {
        slong N = mpoly_words_per_exp(A->bits, ctx->minfo);

        if (old_alloc == 0)
        {
            A->exps = (ulong *) flint_malloc(new_alloc*N*sizeof(ulong));
            A->coeffs = (fmpz_mod_poly_struct *) flint_malloc(new_alloc*sizeof(fmpz_mod_poly_struct));
        } else
        {
            A->exps = (ulong *) flint_realloc(A->exps, new_alloc*N*sizeof(ulong));
            A->coeffs = (fmpz_mod_poly_struct *) flint_realloc(A->coeffs, new_alloc*sizeof(fmpz_mod_poly_struct));
        }

        for (i = old_alloc; i < new_alloc; i++)
        {
            fmpz_mod_poly_init(A->coeffs + i, fpctx->p);
        }
        A->alloc = new_alloc;
    }
}




void fmpz_mod_mpolyun_init(
    fmpz_mod_mpolyun_t A,
    mp_bitcnt_t bits,
    const fmpz_mpoly_ctx_t ctx,
    const fmpz_mod_ctx_t fpctx)
{
    A->coeffs = NULL;
    A->exps = NULL;
    A->alloc = 0;
    A->length = 0;
    A->bits = bits;
}

void fmpz_mod_mpolyun_clear(
    fmpz_mod_mpolyun_t A,
    const fmpz_mpoly_ctx_t ctx,
    const fmpz_mod_ctx_t fpctx)
{
    slong i;
    for (i = 0; i < A->alloc; i++)
        fmpz_mod_mpolyn_clear(A->coeffs + i, ctx, fpctx);
    if (A->coeffs)
        flint_free(A->coeffs);
    if (A->exps)
        flint_free(A->exps);
}

/*
    get the leading coeff in x_0,...,x_var
    A is in R[x_0, ... x_(var-1)][x_var]
*/
fmpz * fmpz_mod_mpolyn_leadcoeff_last_ref(
    fmpz_mod_mpolyn_t A,
    const fmpz_mpoly_ctx_t ctx,
    const fmpz_mod_ctx_t fpctx)
{
    fmpz_mod_poly_struct * leadpoly;
    FLINT_ASSERT(A->length > 0);
    leadpoly = A->coeffs + 0;
    FLINT_ASSERT(leadpoly->length > 0);
    return leadpoly->coeffs + leadpoly->length - 1;
}

fmpz * fmpz_mod_mpolyun_leadcoeff_last_ref(
    fmpz_mod_mpolyun_t A,
    const fmpz_mpoly_ctx_t ctx,
    const fmpz_mod_ctx_t fpctx)
{
    FLINT_ASSERT(A->length > 0);
    return fmpz_mod_mpolyn_leadcoeff_last_ref(A->coeffs + 0, ctx, fpctx);
}

fmpz_mod_poly_struct * fmpz_mod_mpolyn_leadcoeff_ref(
    fmpz_mod_mpolyn_t A,
    const fmpz_mpoly_ctx_t ctx,
    const fmpz_mod_ctx_t fpctx)
{
    FLINT_ASSERT(A->length > 0);
    return A->coeffs + 0;
}

fmpz_mod_poly_struct * fmpz_mod_mpolyun_leadcoeff_ref(
    fmpz_mod_mpolyun_t A,
    const fmpz_mpoly_ctx_t ctx,
    const fmpz_mod_ctx_t fpctx)
{
    FLINT_ASSERT(A->length > 0);
    return fmpz_mod_mpolyn_leadcoeff_ref(A->coeffs + 0, ctx, fpctx);
}

void fmpz_mod_mpolyun_swap(fmpz_mod_mpolyun_t A, fmpz_mod_mpolyun_t B)
{
   fmpz_mod_mpolyun_struct t = *A;
   *A = *B;
   *B = t;
}

void fmpz_mod_mpolyun_zero(
    fmpz_mod_mpolyun_t A,
    const fmpz_mpoly_ctx_t ctx,
    const fmpz_mod_ctx_t fpctx)
{
    slong i;
    for (i = 0; i < A->alloc; i++) {
        fmpz_mod_mpolyn_clear(A->coeffs + i, ctx, fpctx);
        fmpz_mod_mpolyn_init(A->coeffs + i, A->bits, ctx, fpctx);
    }
    A->length = 0;
}

void fmpz_mod_mpolyun_print_pretty(
    const fmpz_mod_mpolyun_t poly,
    const char ** x,
    const fmpz_mpoly_ctx_t ctx,
    const fmpz_mod_ctx_t fpctx)
{
    slong i;
    if (poly->length == 0)
        flint_printf("0");
    for (i = 0; i < poly->length; i++)
    {
        if (i != 0)
            flint_printf(" + ");
        flint_printf("(");
        FLINT_ASSERT((poly->coeffs + i)->bits == poly->bits);
        fmpz_mod_mpolyn_print_pretty(poly->coeffs + i,x,ctx, fpctx);
        flint_printf(")*X^%wd",poly->exps[i]);
    }
}

void fmpz_mod_mpolyun_fit_length(
    fmpz_mod_mpolyun_t A,
    slong length,
    const fmpz_mpoly_ctx_t ctx,
    const fmpz_mod_ctx_t fpctx)
{
    slong i;
    slong old_alloc = A->alloc;
    slong new_alloc = FLINT_MAX(length, 2*A->alloc);

    if (length > old_alloc)
    {
        if (old_alloc == 0)
        {
            A->exps = (ulong *) flint_malloc(new_alloc*sizeof(ulong));
            A->coeffs = (fmpz_mod_mpolyn_struct *) flint_malloc(
                                          new_alloc*sizeof(fmpz_mod_mpolyn_struct));
        } else
        {
            A->exps = (ulong *) flint_realloc(A->exps,
                                                      new_alloc*sizeof(ulong));
            A->coeffs = (fmpz_mod_mpolyn_struct *) flint_realloc(A->coeffs,
                                          new_alloc*sizeof(fmpz_mod_mpolyn_struct));
        }

        for (i = old_alloc; i < new_alloc; i++)
        {
            fmpz_mod_mpolyn_init(A->coeffs + i, A->bits, ctx, fpctx);
        }
        A->alloc = new_alloc;
    }
}



void fmpz_mod_mpolyun_content_last(
    fmpz_mod_poly_t a,
    const fmpz_mod_mpolyun_t B,
    const fmpz_mpoly_ctx_t ctx,
    const fmpz_mod_ctx_t fpctx)
{
    slong i, j;
    fmpz_mod_poly_t t;

    fmpz_mod_poly_init(t, fpctx->p);

    fmpz_mod_poly_zero(a);
    for (i = 0; i < B->length; i++)
    {
        for (j = 0; j < (B->coeffs + i)->length; j++)
        {
            fmpz_mod_poly_gcd(t, a, (B->coeffs + i)->coeffs + j);
            fmpz_mod_poly_swap(t, a);
            if (fmpz_mod_poly_degree(a) == 0)
                break;
        }
    }

    fmpz_mod_poly_clear(t);
}


void fmpz_mod_mpolyun_divexact_last(
    fmpz_mod_mpolyun_t A,
    const fmpz_mod_poly_t b,
    const fmpz_mpoly_ctx_t ctx,
    const fmpz_mod_ctx_t fpctx)
{
    slong i, j;
    fmpz_mod_poly_t r, t;

    fmpz_mod_poly_init(r, fpctx->p);
    fmpz_mod_poly_init(t, fpctx->p);

    for (i = 0; i < A->length; i++)
    {
        fmpz_mod_poly_struct * Ac = (A->coeffs + i)->coeffs;
        for (j = 0; j < (A->coeffs + i)->length; j++)
        {
            fmpz_mod_poly_divrem(t, r, Ac + j, b);
            FLINT_ASSERT(fmpz_mod_poly_is_zero(r));
            FLINT_ASSERT(!fmpz_mod_poly_is_zero(t));
            fmpz_mod_poly_swap(t, Ac + j);
        }
    }
    fmpz_mod_poly_clear(r);
    fmpz_mod_poly_clear(t);
}

void fmpz_mod_mpolyun_mul_last(
    fmpz_mod_mpolyun_t A,
    fmpz_mod_poly_t b,
    const fmpz_mpoly_ctx_t ctx,
    const fmpz_mod_ctx_t fpctx)
{
    slong i, j;
    fmpz_mod_poly_t t;

    fmpz_mod_poly_init(t, fpctx->p);

    for (i = 0; i < A->length; i++)
    {
        for (j = 0; j < (A->coeffs + i)->length; j++)
        {
            fmpz_mod_poly_mul(t, (A->coeffs + i)->coeffs + j, b);
            fmpz_mod_poly_swap(t, (A->coeffs + i)->coeffs + j);
        }
    }

    fmpz_mod_poly_clear(t);
}



slong fmpz_mod_mpolyun_lastdeg(
    const fmpz_mod_mpolyun_t A,
    const fmpz_mpoly_ctx_t ctx,
    const fmpz_mod_ctx_t fpctx)
{
    slong i, j;
    slong deg = -WORD(1);

    for (i = 0; i < A->length; i++)
    {
        for (j = 0; j < (A->coeffs + i)->length; j++)
        {
            slong newdeg = fmpz_mod_poly_degree((A->coeffs + i)->coeffs + j);
            deg = FLINT_MAX(deg, newdeg);
        }
    }

    return deg;
}

/*
    E = A(v = alpha)
    A is in R[X][v]
    E is in R[X]
*/
void fmpz_mod_mpolyun_eval_last_bivar(
    fmpz_mod_poly_t E,
    const fmpz_mod_mpolyun_t A,
    const fmpz_t alpha,
    const fmpz_mpoly_ctx_t ctx,
    const fmpz_mod_ctx_t fpctx)
{
    fmpz_t v;
    slong Ai, Alen;
    fmpz_mod_mpolyn_struct * Acoeff;
    ulong * Aexp;

    fmpz_init(v);

    Acoeff = A->coeffs;
    Aexp = A->exps;
    Alen = A->length;
    Ai = 0;
    fmpz_mod_poly_zero(E);
    for (Ai = 0; Ai < Alen; Ai++)
    {
        fmpz_mod_poly_evaluate_fmpz(v, (Acoeff + Ai)->coeffs + 0, alpha);
        fmpz_mod_poly_set_coeff_fmpz(E, Aexp[Ai], v);
    }

    fmpz_clear(v);
}



void fmpz_mod_mpolyn_one(
    fmpz_mod_mpolyn_t A,
    const fmpz_mpoly_ctx_t ctx,
    const fmpz_mod_ctx_t fpctx)
{
    fmpz_mod_poly_struct * Acoeff;
    ulong * Aexp;
    slong N;

    fmpz_mod_mpolyn_fit_length(A, 1, ctx, fpctx);
    Acoeff = A->coeffs;
    Aexp = A->exps;

    N = mpoly_words_per_exp(A->bits, ctx->minfo);

    fmpz_mod_poly_set_ui(Acoeff + 0, 1);
    mpoly_monomial_zero(Aexp + N*0, N);

    A->length = 1;
}

void fmpz_mod_mpolyun_one(
    fmpz_mod_mpolyun_t A,
    const fmpz_mpoly_ctx_t ctx,
    const fmpz_mod_ctx_t fpctx)
{
    fmpz_mod_mpolyun_fit_length(A, 1, ctx, fpctx);
    fmpz_mod_mpolyn_one(A->coeffs + 0, ctx, fpctx);
    A->exps[0] = 0;
    A->length = 1;
}





void fmpz_mod_mpolyn_scalar_mul_fmpz_mod(
    fmpz_mod_mpolyn_t A,
    const fmpz_t c,
    const fmpz_mpoly_ctx_t ctx,
    const fmpz_mod_ctx_t fpctx)
{
    slong i;
    for (i = 0; i < A->length; i++)
    {
        fmpz_mod_poly_scalar_mul_fmpz(A->coeffs + i, A->coeffs + i, c);
    }
}

void fmpz_mod_mpolyun_scalar_mul_fmpz_mod(
    fmpz_mod_mpolyun_t A,
    const fmpz_t c,
    const fmpz_mpoly_ctx_t ctx,
    const fmpz_mod_ctx_t fpctx)
{
    slong i;
    FLINT_ASSERT(!fmpz_is_zero(c));
    for (i = 0; i < A->length; i++)
    {
        fmpz_mod_mpolyn_scalar_mul_fmpz_mod(A->coeffs + i, c, ctx, fpctx);
    }
}


/*
    Assumption on fmpz_mod_dlog_env:
        p is prime.
        The prime factors of p - 1 all fit a ulong.
    The assumption p is prime can be removed, but then phi(p) needs to be calculated by someone somewhere.
    If the prime factors of p - 1 do not fit a ulong you do not want to calculate dlog mod p.

    so we have
    p - 1 = p1^e1 * ... * pn^en for ulong pi and ei
*/

typedef struct {
    fmpz_t gammapow;
    ulong cm;
} fmpz_mod_dlog_table_entry_struct;

static int fmpz_mod_dlog_table_entry_struct_cmp(const fmpz_mod_dlog_table_entry_struct * lhs,
                                       const fmpz_mod_dlog_table_entry_struct * rhs)
{
    return fmpz_cmp(lhs->gammapow, rhs->gammapow);
}


typedef struct {
    slong exp;
    ulong prime;
    fmpz_t gamma;
    fmpz_t gammainv;
    fmpz_t startingbeta;
    fmpz_t co;
    fmpz_t startinge;
    fmpz_t idem;
    ulong cbound;
    ulong dbound;
    fmpz_mod_dlog_table_entry_struct * table; /* length cbound */
} fmpz_mod_dlog_entry_struct;

typedef struct {
    fmpz_mod_ctx_t fpctx;
    fmpz_t pm1;      /* p - 1 */
    fmpz_t alpha;    /* p.r. of p */
    fmpz_t alphainv;
    slong num_factors;  /* factors of p - 1*/
    fmpz_mod_dlog_entry_struct * entries;
} fmpz_mod_dlog_env_struct;
typedef fmpz_mod_dlog_env_struct fmpz_mod_dlog_env_t[1];

/* assume that p is prime, don't check */
void fmpz_mod_dlog_env_init(fmpz_mod_dlog_env_t L, const fmpz_t p)
{
    slong i;
    ulong c;
    fmpz_mod_dlog_entry_struct * Li;
    fmpz_factor_t factors;
    fmpz_t temp;

    fmpz_init(L->alpha);
    fmpz_init(L->alphainv);
    fmpz_init(L->pm1);
    fmpz_mod_ctx_init(L->fpctx, p);

    fmpz_init(temp);

    fmpz_factor_init(factors);
    fmpz_sub_ui(L->pm1, p, 1);
    fmpz_factor(factors, L->pm1);
    L->num_factors = factors->num;
    L->entries = NULL;
    if (L->num_factors > 0)
    {
        L->entries = (fmpz_mod_dlog_entry_struct*) flint_malloc(L->num_factors*sizeof(fmpz_mod_dlog_entry_struct));
    }
    for (i = 0; i < L->num_factors; i++)
    {
        int success;
        fmpz_t pipow, recp;

        Li = L->entries + i;

        fmpz_init(Li->idem);
        fmpz_init(Li->co);
        fmpz_init(Li->startinge);
        fmpz_init(Li->startingbeta);
        fmpz_init(Li->gamma);
        fmpz_init(Li->gammainv);

        FLINT_ASSERT(fmpz_abs_fits_ui(factors->p + i));
        Li->exp = factors->exp[i];
        Li->prime = fmpz_get_ui(factors->p + i);
/*
flint_printf("L[%wd].exp  : %wu\n", i, Li->exp);
flint_printf("L[%wd].prime: %wu\n", i, Li->prime);
*/
        fmpz_init(recp);
        fmpz_init_set_ui(pipow, Li->prime);
        fmpz_pow_ui(pipow, pipow, Li->exp);
        fmpz_divexact(recp, L->pm1, pipow);
        success = fmpz_invmod(temp, recp, pipow);
        FLINT_ASSERT(success);
        fmpz_mul(temp, temp, recp);

        fmpz_mod(Li->idem, temp, L->pm1);

        fmpz_set(Li->co, recp);
        fmpz_divexact_ui(Li->startinge, pipow, Li->prime);

        fmpz_clear(pipow);
        fmpz_clear(recp);
/*
flint_printf("L[%wd].idem: ", i); fmpz_print(Li->idem); printf("\n");
flint_printf("L[%wd].co  : ", i); fmpz_print(Li->co); printf("\n");
flint_printf("L[%wd].stae: ", i); fmpz_print(Li->startinge); printf("\n");
*/
    }
    fmpz_factor_clear(factors);

    fmpz_one(L->alpha);
try_alpha:
    fmpz_add_ui(L->alpha, L->alpha, 1);
/*
flint_printf("\nL.alpha   : ", i); fmpz_print(L->alpha); printf("\n");
*/

    if (fmpz_cmp(L->alpha, p) >= 0)
    {
        flint_printf("Exception in fmpz_mod_dlog_env_init: could not find primitive root\n");
        flint_abort();
    }
    for (i = 0; i < L->num_factors; i++)
    {
        Li = L->entries + i;
        fmpz_divexact_ui(temp, L->pm1, Li->prime);
/*
flint_printf("temp[%wd] : ", i); fmpz_print(temp); printf("\n");
*/

        fmpz_mod_pow_fmpz(Li->gamma, L->alpha, temp, L->fpctx);

/*
flint_printf("gamma[%wd]: ", i); fmpz_print(Li->gamma); printf("\n");
*/

        if (fmpz_is_one(Li->gamma))
        {
            goto try_alpha;
        }
    }

printf("prime "); fmpz_print(p); printf(" alpha "); fmpz_print(L->alpha); printf("\n");

    fmpz_mod_inv(L->alphainv, L->alpha, L->fpctx);
/*
flint_printf("L.alphainv: ", i); fmpz_print(L->alphainv); printf("\n");
*/
    for (i = 0; i < L->num_factors; i++)
    {
        Li = L->entries + i;
/*
flint_printf("L[%wd].p      : ", i); fmpz_print(L->fpctx->p); printf("\n");
flint_printf("L[%wd].gamma  : ", i); fmpz_print(Li->gamma); printf("\n");
*/
        fmpz_mod_inv(Li->gammainv, Li->gamma, L->fpctx);
/*
flint_printf("L[%wd].invgmma: ", i); fmpz_print(Li->gammainv); printf("\n");
*/
        fmpz_mod_pow_fmpz(Li->startingbeta, L->alphainv, Li->co, L->fpctx);
/*
flint_printf("L[%wd].stabeta: ", i); fmpz_print(Li->startinge); printf("\n");
*/

        Li->dbound = ceil(sqrt((double) Li->prime));
        Li->cbound = (Li->prime + Li->dbound - 1)/Li->dbound;
        while (Li->cbound > 100)
        {
            Li->dbound *= 2;
            Li->cbound = (Li->prime + Li->dbound - 1)/Li->dbound;
        }

        FLINT_ASSERT(Li->dbound > 0);
        FLINT_ASSERT(Li->cbound > 0);
        Li->table = (fmpz_mod_dlog_table_entry_struct *) flint_malloc(Li->cbound*sizeof(fmpz_mod_dlog_table_entry_struct));

        for (c = 0; c < Li->cbound; c++)
        {
            Li->table[c].cm = c*Li->dbound;
            fmpz_init(Li->table[c].gammapow);
            fmpz_mod_pow_ui(Li->table[c].gammapow, Li->gamma, Li->table[c].cm, L->fpctx);
        }
        qsort(Li->table, Li->cbound, sizeof(fmpz_mod_dlog_table_entry_struct),
               (int(*)(const void*, const void*)) fmpz_mod_dlog_table_entry_struct_cmp);
        for (c = 1; c < Li->cbound; c++)
        {
            FLINT_ASSERT(fmpz_cmp(Li->table[c - 1].gammapow, Li->table[c].gammapow) < 0);
        }
    }

    fmpz_clear(temp);
}

void fmpz_mod_dlog_env_clear(fmpz_mod_dlog_env_t L)
{
    slong i;
    ulong c;
    fmpz_mod_dlog_entry_struct * Li;

    for (i = 0; i < L->num_factors; i++)
    {
        Li = L->entries + i;
        fmpz_clear(Li->idem);
        fmpz_clear(Li->co);
        fmpz_clear(Li->startinge);
        fmpz_clear(Li->startingbeta);
        fmpz_clear(Li->gamma);
        fmpz_clear(Li->gammainv);
        for (c = 0; c < Li->cbound; c++)
        {
            fmpz_clear(Li->table[c].gammapow);
        }
        flint_free(Li->table);
    }

    if (L->entries)
    {
        flint_free(L->entries);
    }

    fmpz_clear(L->alpha);
    fmpz_clear(L->alphainv);
    fmpz_clear(L->pm1);

    fmpz_mod_ctx_clear(L->fpctx);

    return;
}

/* return x such that x = alpha^y mod p, alpha is the p.r. L->alpha*/
void fmpz_mod_dlog_env_run(const fmpz_mod_dlog_env_t L, fmpz_t xx, const fmpz_t y)
{
    slong i, j;
/*    ulong x, q, r, e, x0 = 0, x1 = 0, x2 = 0, pp0, pp1, acc, g, pipow;*/
    ulong g;
    fmpz_t x;
    fmpz_t pipow, e, acc;
    ulong lo, mid, hi, d;
    fmpz_t beta, z, w, temp;
    fmpz_mod_dlog_entry_struct * Li;

    fmpz_init(x);
    fmpz_init(acc);
    fmpz_init(pipow);
    fmpz_init(e);
    fmpz_init(beta);
    fmpz_init(z);
    fmpz_init(w);
    fmpz_init(temp);

    FLINT_ASSERT(!fmpz_is_zero(y));
    FLINT_ASSERT(fmpz_mod_is_canonical(y, L->fpctx));

    i = 0;
    if (i < L->num_factors && L->entries[i].prime == 2)
    {
        Li = L->entries + i;
        FLINT_ASSERT(Li->prime == 2);

        fmpz_mod_pow_fmpz(z, y, Li->co, L->fpctx);
        fmpz_set(beta, Li->startingbeta);
        fmpz_set(e, Li->startinge);
        j = 0;
        fmpz_one(pipow); /* Li->prime^j */
        fmpz_zero(acc);
        do {
            fmpz_mod_pow_fmpz(w, z, e, L->fpctx);
            /* solve Li->gamma ^ g == w mod p */
            if (fmpz_is_one(w))
            {
                g = 0;
            }
            else
            {
                if(!fmpz_equal(w, Li->gamma))
                {
                    flint_printf("Exception in fmpz_mod_dlog_env_run: could not find log\n");
                    flint_abort();
                }
                g = 1;
                fmpz_mod_mul(z, z, beta, L->fpctx);
            }
            fmpz_mod_mul(beta, beta, beta, L->fpctx);
            fmpz_addmul_ui(acc, pipow, g);
            fmpz_mul_2exp(pipow, pipow, 1);
            fmpz_tdiv_q_2exp(e, e, 1);
        } while (++j < Li->exp);

        fmpz_addmul(x, acc, Li->idem);
        i = 1;
    }

    for (; i < L->num_factors; i++)
    {
        Li = L->entries + i;

        fmpz_mod_pow_fmpz(z, y, Li->co, L->fpctx);
        fmpz_set(beta, Li->startingbeta);
        fmpz_set(e, Li->startinge);
        j = 0;
        fmpz_one(pipow); /* Li->prime^j */
        fmpz_zero(acc);
        do {
            fmpz_mod_pow_fmpz(w, z, e, L->fpctx);
            /* solve Li->gamma ^ g == w mod p */
            d = 0;
            while (1)
            {
                lo = 0; hi = Li->cbound;
                while (hi - lo > 4)
                {
                    int cmp;
                    mid = lo + (hi - lo)/2;
                    cmp = fmpz_cmp(Li->table[mid].gammapow, w);
                    if (cmp == 0)
                    {
                        g = Li->table[mid].cm + d;
                        goto found_g;
                    }
                    else if (cmp > 0)
                    {
                        hi = mid;
                    }
                    else
                    {
                        lo = mid;
                    }
                }
                while (lo < hi)
                {
                    if (fmpz_equal(Li->table[lo].gammapow, w))
                    {
                        g = Li->table[lo].cm + d;
                        goto found_g;
                    }
                    lo++;
                }
                fmpz_mod_mul(w, w, Li->gammainv, L->fpctx);
                d++;
                if (d >= Li->dbound)
                {
                    flint_printf("Exception in fmpz_mod_dlog_env_run: could not find log\n");
                    flint_abort();
                }
            }
        found_g:
            FLINT_ASSERT(g < Li->prime);
            fmpz_mod_pow_ui(temp, beta, g, L->fpctx);
            fmpz_mod_mul(z, z, temp, L->fpctx);
            fmpz_mod_pow_ui(beta, beta, Li->prime, L->fpctx);
            fmpz_addmul_ui(acc, pipow, g);
            fmpz_mul_ui(pipow, pipow, Li->prime);
            fmpz_divexact_ui(e, e, Li->prime);
        } while (++j < Li->exp);

        fmpz_addmul(x, acc, Li->idem);
    }

    fmpz_mod(xx, x, L->pm1);
    fmpz_clear(acc);
    fmpz_clear(pipow);
    fmpz_clear(e);
    fmpz_clear(beta);
    fmpz_clear(z);
    fmpz_clear(w);
    fmpz_clear(temp);
    fmpz_clear(x);
    return;
}

typedef struct {
    slong half_point_count;
    fmpz_mod_poly_t R0, R1;
    fmpz_mod_poly_t V0, V1; /* V1 is our master polynomial */
    fmpz_mod_poly_t r; /* temporary */
    fmpz_mod_poly_t q; /* temporary also used as a queue of incoming points */
} fmpz_mod_bma_struct;
typedef fmpz_mod_bma_struct fmpz_mod_bma_t[1];

/*
    A = x^(2*n)       n = half_point_count
    deg(S) < 2*n
    U0*A + V0*S = R0   deg(R0) >= n
    U1*A + V1*S = R1   deg(R1) < n

    S is the reverse of the polynomial whose coefficients are the input points.
        S = a_0*x^(2n-1) + a_1*x^(2n-2) + ... + a_(2n-1)
    S can be updated with more points at any time via add_points.
    S can be reduced (compute V0, V1, R0, R1) at any time via reduce, which
        returns whether the reduction caused a change.

    The U0 and U1 are not stored.
*/
void fmpz_mod_bma_init(fmpz_mod_bma_t B, const fmpz_t p)
{
    B->half_point_count = 0;
    fmpz_mod_poly_init(B->V0, p);
    fmpz_mod_poly_init(B->R0, p);
    fmpz_mod_poly_set_ui(B->R0, 1);
    fmpz_mod_poly_init(B->V1, p);
    fmpz_mod_poly_set_ui(B->V1, 1);
    fmpz_mod_poly_init(B->R1, p);
    fmpz_mod_poly_init(B->r, p);
    fmpz_mod_poly_init(B->q, p);
    B->q->length = 0;
}

void fmpz_mod_bma_start_over(fmpz_mod_bma_t B)
{
    B->half_point_count = 0;
    fmpz_mod_poly_zero(B->V0);
    fmpz_mod_poly_set_ui(B->R0, 1);
    fmpz_mod_poly_set_ui(B->V1, 1);
    fmpz_mod_poly_zero(B->R1);
    B->q->length = 0;
}

void fmpz_mod_bma_reset_prime(fmpz_mod_bma_t B, const fmpz_t p)
{
    fmpz_set(&B->R0->p, p);
    fmpz_set(&B->R1->p, p);
    fmpz_set(&B->V0->p, p);
    fmpz_set(&B->V1->p, p);
    fmpz_set(&B->r->p, p);
    fmpz_set(&B->q->p, p);
    fmpz_mod_bma_start_over(B);
}

void fmpz_mod_bma_clear(fmpz_mod_bma_t B)
{
    fmpz_mod_poly_clear(B->R0);
    fmpz_mod_poly_clear(B->R1);
    fmpz_mod_poly_clear(B->V0);
    fmpz_mod_poly_clear(B->V1);
    fmpz_mod_poly_clear(B->r);
    fmpz_mod_poly_clear(B->q);
}

void fmpz_mod_bma_add_points(fmpz_mod_bma_t B, const fmpz * a, slong count)
{
    slong i;
    slong old_length = B->q->length;
    fmpz_mod_poly_fit_length(B->q, old_length + count);
    for (i = 0; i < count; i++)
    {
        fmpz_set(B->q->coeffs + old_length + i, a + i);
    }
    B->q->length = old_length + count;
}

void fmpz_mod_bma_add_point(fmpz_mod_bma_t B, const fmpz_t a)
{
    slong old_length = B->q->length;
    fmpz_mod_poly_fit_length(B->q, old_length + 1);
    fmpz_set(B->q->coeffs + old_length, a);
    B->q->length = old_length + 1;
}

int fmpz_mod_bma_reduce(fmpz_mod_bma_t B)
{
    int changed = 0;
    slong i, queue_length = B->q->length;

    if ((queue_length % 2) != 0)
    {
        flint_printf("Exception in nmod_bma_reduce: point count is not even\n");
        flint_abort();
    }

    /* reverse the queue into temp r */
    B->half_point_count += queue_length/2;
    fmpz_mod_poly_zero(B->r);
    for (i = 0; i < queue_length; i++)
    {
        fmpz_mod_poly_set_coeff_fmpz(B->r, queue_length - i - 1, B->q->coeffs + i);
    }
    fmpz_mod_poly_mul(B->q, B->V0, B->r);
    fmpz_mod_poly_shift_left(B->R0, B->R0, queue_length);
    fmpz_mod_poly_add(B->R0, B->R0, B->q);
    fmpz_mod_poly_mul(B->q, B->V1, B->r);
    fmpz_mod_poly_shift_left(B->R1, B->R1, queue_length);
    fmpz_mod_poly_add(B->R1, B->R1, B->q);

    /* now reduce */
    while (B->half_point_count < fmpz_mod_poly_length(B->R1))
    {
        changed = 1;
        fmpz_mod_poly_divrem(B->q, B->r, B->R0, B->R1);
        fmpz_mod_poly_swap(B->R0, B->R1);
        fmpz_mod_poly_swap(B->R1, B->r);

        fmpz_mod_poly_mul(B->r, B->q, B->V1);
        fmpz_mod_poly_sub(B->q, B->V0, B->r);
        fmpz_mod_poly_swap(B->V0, B->V1);
        fmpz_mod_poly_swap(B->V1, B->q);
        FLINT_ASSERT(fmpz_mod_poly_degree(B->V1) > fmpz_mod_poly_degree(B->V0));
    }

    /* queue is empty now */
    B->q->length = 0;
    return changed;
}



/* split f assuming that f has degree(f) distinct nonzero roots in Fp */
static void _fmpz_mod_poly_rabinsplit(fmpz_mod_poly_t a, fmpz_mod_poly_t b, fmpz_mod_poly_t T,
                                  const fmpz_mod_poly_t f, flint_rand_t randstate)
{
    fmpz_t delta;

    fmpz_init(delta);

    FLINT_ASSERT(fmpz_mod_poly_degree(f) > 1);

try_again:

    fmpz_randm(delta, randstate, &f->p);

    fmpz_mod_poly_zero(a);
    fmpz_mod_poly_set_coeff_ui(a, 1, 1);
    fmpz_mod_poly_set_coeff_fmpz(a, 0, delta);
    fmpz_sub_ui(delta, &f->p, 1);
    fmpz_divexact_ui(delta, delta, 2);
    fmpz_mod_poly_powmod_fmpz_binexp(T, a, delta, f);
    fmpz_mod_poly_zero(a);
    fmpz_mod_poly_set_coeff_ui(a, 0, 1);
    fmpz_mod_poly_sub(T, T, a);
    fmpz_mod_poly_gcd(a, T, f);
    FLINT_ASSERT(!fmpz_mod_poly_is_zero(a));
    if (0 >= fmpz_mod_poly_degree(a) || fmpz_mod_poly_degree(a) >= fmpz_mod_poly_degree(f))
    {
        goto try_again;
    }
    fmpz_mod_poly_div_basecase(b, f, a);
    /* deg a >= deg b */
    if (fmpz_mod_poly_degree(a) < fmpz_mod_poly_degree(b))
    {
        fmpz_mod_poly_swap(a, b);
    }

    fmpz_clear(delta);
    return;
}

/* fill in roots with the t distinct nonzero roots of master, or fail */
static int _fmpz_mod_find_roots(fmpz * roots, const fmpz_mod_poly_t master, slong t, const fmpz_mod_ctx_t fpctx)
{
    fmpz_t a0, a1;
    int success;
    slong i, roots_idx, sp;
    fmpz_t delta;
    fmpz_mod_poly_struct * a , * b;
    fmpz_mod_poly_t f, T;
    fmpz_mod_poly_struct stack[FLINT_BITS + 1];
    flint_rand_t randstate;

    FLINT_ASSERT(t >= 0);
    FLINT_ASSERT(t == fmpz_mod_poly_degree(master));

    fmpz_init(a0);
    fmpz_init(a1);
    fmpz_init(delta);

    if (t == 0)
    {
        success = 1;
        goto cleanup1;
    }
    else if (t == 1)
    {
        fmpz_mod_poly_get_coeff_fmpz(a0, master, 0);
        fmpz_mod_poly_get_coeff_fmpz(a1, master, 1);
        if (fmpz_is_zero(a0))
        {
            success = 0;
            goto cleanup1;
        }
        fmpz_mod_inv(a1, a1, fpctx);
        fmpz_mod_neg(a1, a1, fpctx);
        fmpz_mod_mul(roots + 0, a0, a1, fpctx);
        success = 1;
        goto cleanup1;
    }

    flint_randinit(randstate);
    fmpz_mod_poly_init(T, &master->p);
    fmpz_mod_poly_init(f, &master->p);
    for (i = 0; i <= FLINT_BITS; i++)
    {
        fmpz_mod_poly_init(stack + i, &master->p);
    }

    roots_idx = 0;

    fmpz_mod_poly_make_monic(f, master);

    a = stack + 0;
    fmpz_randm(delta, randstate, &master->p);
    fmpz_mod_poly_zero(a);
    fmpz_mod_poly_set_coeff_ui(a, 1, 1);
    fmpz_mod_poly_set_coeff_fmpz(a, 0, delta);
    fmpz_sub_ui(delta, &master->p, 1);
    fmpz_divexact_ui(delta, delta, 2);
    fmpz_mod_poly_powmod_fmpz_binexp(T, a, delta, f);
    fmpz_mod_poly_zero(a);
    fmpz_mod_poly_set_coeff_ui(a, 0, 1);
    fmpz_mod_poly_sub(T, T, a);
    fmpz_mod_poly_gcd(a, T, f);

    b = stack + 1;
    fmpz_mod_poly_zero(b);
    fmpz_mod_poly_set_coeff_ui(b, 0, 2);
    fmpz_mod_poly_add(T, T, b);
    fmpz_mod_poly_gcd(b, T, f);

    if (fmpz_mod_poly_degree(b) + fmpz_mod_poly_degree(a) != t)
    {
        success = 0;
        goto cleanup;
    }
    /* deg a >= deg b */
    if (fmpz_mod_poly_degree(a) < fmpz_mod_poly_degree(b))
    {
        fmpz_mod_poly_swap(a, b);
    }

    sp = fmpz_mod_poly_degree(b) > 0 ? 2 : 1;
    while (sp > 0)
    {
        FLINT_ASSERT(sp < FLINT_BITS);
        sp--;
        fmpz_mod_poly_swap(f, stack + sp);

        FLINT_ASSERT(fmpz_mod_poly_degree(f) > 0);
        if (fmpz_mod_poly_degree(f) == 1)
        {
            fmpz_mod_poly_get_coeff_fmpz(a0, f, 0);
            fmpz_mod_poly_get_coeff_fmpz(a1, f, 1);
            FLINT_ASSERT(!fmpz_is_zero(a0));
            FLINT_ASSERT(fmpz_is_one(a1));
            fmpz_mod_neg(roots + roots_idx, a0, fpctx);
            roots_idx++;
        }
        else
        {
            _fmpz_mod_poly_rabinsplit(stack + sp + 0, stack + sp + 1, T, f, randstate);
            FLINT_ASSERT(FLINT_BIT_COUNT(fmpz_mod_poly_degree(stack + sp + 1)) <= FLINT_BITS - sp - 1);
            sp += 2;
        }
    }

    success = 1;

cleanup:

    flint_randclear(randstate);
    fmpz_mod_poly_clear(T);
    fmpz_mod_poly_clear(f);
    for (i = 0; i <= FLINT_BITS; i++)
    {
        fmpz_mod_poly_clear(stack + i);
    }

    if (success)
    {
        FLINT_ASSERT(roots_idx == t);
    }

cleanup1:

    fmpz_clear(a0);
    fmpz_clear(a1);
    fmpz_clear(delta);

    return success;
}


typedef struct
{
    const fmpz_mpoly_ctx_struct * polyctx;
    fmpz_mod_ctx_t fpctx;
    fmpz_mod_dlog_env_t dlogenv;
    slong * degbounds;
    ulong * subdegs;
    mp_bitcnt_t bits;
    ulong * inputexpmask;
} fmpz_mod_mpoly_bma_interpolate_ctx_struct;
typedef fmpz_mod_mpoly_bma_interpolate_ctx_struct fmpz_mod_mpoly_bma_interpolate_ctx_t[1];

typedef struct {
    fmpz_mod_poly_t roots;
    fmpz_mod_poly_t evals;
    fmpz_mod_bma_t bma;
} fmpz_mod_mpoly_bma_interpolate_struct;
typedef fmpz_mod_mpoly_bma_interpolate_struct fmpz_mod_mpoly_bma_interpolate_t[1];



void fmpz_mod_mpoly_bma_interpolate_ctx_init(
    fmpz_mod_mpoly_bma_interpolate_ctx_t Ictx,/* mp_bitcnt_t inputbits,*/
    const fmpz_mpoly_ctx_t pctx,
    const fmpz_t p)
{
/*
    slong N;
*/

    fmpz_mod_ctx_init(Ictx->fpctx, p);
    Ictx->polyctx = pctx;
    Ictx->degbounds = (slong *) flint_malloc(pctx->minfo->nvars*sizeof(slong));
    Ictx->subdegs = (ulong *) flint_malloc(pctx->minfo->nvars*sizeof(ulong));
    fmpz_mod_dlog_env_init(Ictx->dlogenv, p);

/*
    Ictx->inbits = inputbits;
    N = mpoly_words_per_exp_sp(inputbits, pctx->minfo);
    Ictx->inputexpmask = (ulong *) flint_malloc(N*sizeof(ulong));
    mpoly_monomial_zero(Ictx->inputexpmask, N);
*/
}

void fmpz_mod_mpoly_bma_interpolate_ctx_reset_prime(
    fmpz_mod_mpoly_bma_interpolate_ctx_t Ictx,
    const fmpz_t p)
{
    fmpz_mod_ctx_clear(Ictx->fpctx);
    fmpz_mod_dlog_env_clear(Ictx->dlogenv);

    fmpz_mod_ctx_init(Ictx->fpctx, p);
    fmpz_mod_dlog_env_init(Ictx->dlogenv, p);
}

void fmpz_mod_mpoly_bma_interpolate_ctx_clear(fmpz_mod_mpoly_bma_interpolate_ctx_t Ictx)
{
    flint_free(Ictx->inputexpmask);
    fmpz_mod_dlog_env_clear(Ictx->dlogenv);
    flint_free(Ictx->degbounds);
    flint_free(Ictx->subdegs);
    fmpz_mod_ctx_clear(Ictx->fpctx);
}


void fmpz_mod_mpoly_bma_interpolate_init(fmpz_mod_mpoly_bma_interpolate_t I, const fmpz_t p)
{
    fmpz_mod_poly_init(I->roots, p);
    fmpz_mod_poly_init(I->evals, p);
    fmpz_mod_bma_init(I->bma, p);
}

slong fmpz_mod_mpoly_bma_interpolate_pointcount(const fmpz_mod_mpoly_bma_interpolate_t I)
{
    return I->bma->q->length + 2*I->bma->half_point_count;
}

void fmpz_mod_mpoly_bma_interpolate_print(const fmpz_mod_mpoly_bma_interpolate_t I)
{
    flint_printf("(%wd) ", fmpz_mod_mpoly_bma_interpolate_pointcount(I));
    fmpz_mod_poly_print_pretty(I->bma->V1, "#");
}


void fmpz_mod_mpoly_bma_interpolate_reset_prime(
    fmpz_mod_mpoly_bma_interpolate_t I,
    const fmpz_t p)
{
    fmpz_set(&I->roots->p, p);
    fmpz_set(&I->evals->p, p);
    fmpz_mod_bma_reset_prime(I->bma, p);
}

void fmpz_mod_mpoly_bma_interpolate_clear(fmpz_mod_mpoly_bma_interpolate_t I)
{
    fmpz_mod_poly_clear(I->roots);
    fmpz_mod_poly_clear(I->evals);
    fmpz_mod_bma_clear(I->bma);
}

void fmpz_mod_mpoly_bma_interpolate_add_point(fmpz_mod_mpoly_bma_interpolate_t I, const fmpz_t a)
{
    fmpz_mod_poly_fit_length(I->evals, I->evals->length + 1);
    fmpz_set(I->evals->coeffs + I->evals->length, a);
    I->evals->length++;

    fmpz_mod_bma_add_point(I->bma, a);
}

void fmpz_mod_mpoly_bma_interpolate_zero(fmpz_mod_mpoly_bma_interpolate_t I, slong count)
{
    slong i;
    fmpz_mod_poly_fit_length(I->evals, count);
    fmpz_mod_bma_start_over(I->bma);
    for (i = 0; i < count; i++)
    {
        fmpz_zero(I->evals->coeffs + i);
        fmpz_mod_bma_add_point(I->bma, I->evals->coeffs + i);
    }
    I->evals->length = count;
}

int fmpz_mod_mpoly_bma_interpolate_reduce(fmpz_mod_mpoly_bma_interpolate_t I)
{
    return fmpz_mod_bma_reduce(I->bma);
}

void fmpz_mod_mpoly_bma_interpolate_eval_init(fmpz_mod_mpoly_bma_interpolate_ctx_t Ictx, const fmpz_mpoly_t A)
{    
    slong i, j, N = mpoly_words_per_exp_sp(Ictx->bits, Ictx->polyctx->minfo);
    FLINT_ASSERT(A->bits == Ictx->bits);

    for (i = 0; i < A->length; i++)
    {
        for (j = 0; j < N; j++)
        {
            Ictx->inputexpmask[j] |= (A->exps + N*i)[j];
        }
    }
}

/*
    put the evaluation of the monomials in A at alpha^w in the coeffs of E

    x_0     = alpha ^ (w * subdegs[n-1] * subdegs[n-2] * ... * * subdegs[1])
      ...
    x_(n-3) = alpha ^ (w * subdegs[n-1] * subdegs[n-2])
    x_(n-2) = alpha ^ (w * subdegs[n-1])
    x_(n-1) = alpha ^ (w)

    secret: subdegs[0] is not used
*/
void fmpz_mod_mpoly_bma_interpolate_eval_setskel(
    fmpz_mpoly_t M,
    const fmpz_mpoly_t A,
    const fmpz_t w,
    const fmpz_mod_mpoly_bma_interpolate_ctx_t Ictx)
{
    slong i, j;
    slong offset, shift;
    slong N = mpoly_words_per_exp_sp(A->bits, Ictx->polyctx->minfo);
    slong nvars = Ictx->polyctx->minfo->nvars;
    ulong * Aexp;
    fmpz * Mcoeff;
    slong * LUToffset;
    ulong * LUTmask;
    fmpz * LUTvalue;
    slong LUTlen;
    fmpz_t xeval, xpoweval;
    ulong * inputexpmask;
    TMP_INIT;

    TMP_START;

    fmpz_init(xeval);
    fmpz_init(xpoweval);

    LUToffset = (slong *) TMP_ALLOC(N*FLINT_BITS*sizeof(slong));
    LUTmask   = (ulong *) TMP_ALLOC(N*FLINT_BITS*sizeof(ulong));
    LUTvalue  = (fmpz *) TMP_ALLOC(N*FLINT_BITS*sizeof(fmpz));
    for (i = 0; i < N*FLINT_BITS; i++)
    {
        fmpz_init(LUTvalue + i);
    }

    FLINT_ASSERT(M->bits == A->bits);

    fmpz_mpoly_fit_length(M, A->length, Ictx->polyctx);
    M->length = A->length;

    Mcoeff = M->coeffs;
    Aexp = A->exps;

    inputexpmask = (ulong *) TMP_ALLOC(N*sizeof(ulong));
    mpoly_monomial_zero(inputexpmask, N);
    for (i = 0; i < A->length; i++)
    {
        for (j = 0; j < N; j++)
        {
            inputexpmask[j] |= (Aexp + N*i)[j];
        }
    }

    LUTlen = 0;
    fmpz_mod_pow_fmpz(xeval, Ictx->dlogenv->alpha, w, Ictx->fpctx);
    for (j = nvars - 1; j >= 0; j--)
    {
        mpoly_gen_offset_shift_sp(&offset, &shift, j, A->bits, Ictx->polyctx->minfo);

        fmpz_set(xpoweval, xeval); /* xpoweval = xeval^(2^i) */
        for (i = 0; i < A->bits; i++)
        {
            LUToffset[LUTlen] = offset;
            LUTmask[LUTlen] = (UWORD(1) << (shift + i));
            fmpz_set(LUTvalue + LUTlen, xpoweval);
            if ((inputexpmask[offset] & LUTmask[LUTlen]) != 0)
            {
                LUTlen++;
            }
            fmpz_mod_mul(xpoweval, xpoweval, xpoweval, Ictx->fpctx);
        }
        fmpz_mod_pow_ui(xeval, xeval, Ictx->subdegs[j], Ictx->fpctx);
    }
    FLINT_ASSERT(LUTlen < N*FLINT_BITS);

    for (i = 0; i < A->length; i++)
    {
        fmpz_one(xeval);
        for (j = 0; j < LUTlen; j++)
        {
            if (((Aexp + N*i)[LUToffset[j]] & LUTmask[j]) != 0)
            {
                fmpz_mod_mul(xeval, xeval, LUTvalue +j, Ictx->fpctx);
            }
        }
        fmpz_set(Mcoeff + i, xeval);
        mpoly_monomial_zero(M->exps + N*i, N);
    }

    fmpz_clear(xeval);
    fmpz_clear(xpoweval);
    for (i = 0; i < N*FLINT_BITS; i++)
    {
        fmpz_clear(LUTvalue + i);
    }
    TMP_END;
}

/* M = S */
void fmpz_mod_mpoly_bma_interpolate_eval_copyskel(
    fmpz_mpoly_t M,
    const fmpz_mpoly_t S,
    const fmpz_mod_mpoly_bma_interpolate_ctx_t Ictx)
{
    slong i, N;

    FLINT_ASSERT(M->bits == Ictx->bits);


    fmpz_mpoly_fit_length(M, S->length, Ictx->polyctx);
    M->length = S->length;
    _fmpz_vec_set(M->coeffs, S->coeffs, S->length);
    N = mpoly_words_per_exp(Ictx->bits, Ictx->polyctx->minfo);
    for (i = 0; i < M->length; i++)
    {
        mpoly_monomial_zero(M->exps + N*i, N);
    }
}

/* return A.M */
void fmpz_mod_mpoly_bma_interpolate_eval_useskel(fmpz_t eval, const fmpz_mpoly_t A, const fmpz_mpoly_t M, const fmpz_mod_mpoly_bma_interpolate_ctx_t Ictx)
{
    slong i;
    fmpz_t t;
    fmpz_init(t);
    fmpz_zero(eval);
    FLINT_ASSERT(M->length == A->length);
    for (i = 0; i < A->length; i++)
    {
        fmpz_mod_mul(t, A->coeffs + i, M->coeffs + i, Ictx->fpctx);
        fmpz_mod_add(eval, eval, t, Ictx->fpctx);
    }
    fmpz_clear(t);
}

/*
    return A.M and multiply M by S
    the coefficients of A are not necesarily reduced mod Ictx->fpctx->p
*/
void fmpz_mod_mpoly_bma_interpolate_eval_useskelmul(
    fmpz_t eval,
    const fmpz_mpoly_t A,
    fmpz_mpoly_t M,
    const fmpz_mpoly_t S,
    const fmpz_mod_mpoly_bma_interpolate_ctx_t Ictx)
{
    slong i;
    fmpz_zero(eval);
    FLINT_ASSERT(M->length == A->length);
    for (i = 0; i < A->length; i++)
    {
        fmpz_addmul(eval, A->coeffs + i, M->coeffs + i);
        fmpz_mod_mul(M->coeffs + i, M->coeffs + i, S->coeffs + i, Ictx->fpctx);
    }
    fmpz_mod(eval, eval, Ictx->fpctx->p);
}


void _fmpz_mpoly_set_skel(
    fmpz_mpoly_t M,
    const fmpz_mpoly_t A,
    const fmpz_mod_mpoly_bma_interpolate_ctx_t Ictx)
{
    fmpz_t one;
    fmpz_init_set_ui(one, 1);
    fmpz_mod_mpoly_bma_interpolate_eval_setskel(M, A, one, Ictx);
    fmpz_clear(one);
}

void _fmpz_mpoly_copy_skel(
    fmpz_mpoly_t M,
    const fmpz_mpoly_t S,
    const fmpz_mod_mpoly_bma_interpolate_ctx_t Ictx)
{
    fmpz_mod_mpoly_bma_interpolate_eval_copyskel(M, S, Ictx);
}

void _fmpz_mpolyuu_set_skel(
    fmpz_mpolyu_t M,
    const fmpz_mpolyu_t A,
    const fmpz_mod_mpoly_bma_interpolate_ctx_t Ictx)
{
    slong i;
    fmpz_t one;
    fmpz_init_set_ui(one, 1);
    fmpz_mpolyu_fit_length(M, A->length, Ictx->polyctx);
    for (i = 0; i < A->length; i++)
    {
        M->exps[i] = 0;
        fmpz_mod_mpoly_bma_interpolate_eval_setskel(M->coeffs + i, A->coeffs + i, one, Ictx);
    }
    M->length = A->length;
    fmpz_clear(one);
}

void _fmpz_mpolyuu_copy_skel(
    fmpz_mpolyu_t M,
    const fmpz_mpolyu_t S,
    const fmpz_mod_mpoly_bma_interpolate_ctx_t Ictx)
{
    slong i;
    fmpz_t one;
    fmpz_init_set_ui(one, 1);
    fmpz_mpolyu_fit_length(M, S->length, Ictx->polyctx);
    for (i = 0; i < S->length; i++)
    {
        M->exps[i] = 0;
        fmpz_mod_mpoly_bma_interpolate_eval_copyskel(M->coeffs + i, S->coeffs + i, Ictx);
    }
    M->length = S->length;
    fmpz_clear(one);
}


void _fmpz_mpoly_use_skel_mul(
    fmpz_t eval,
    const fmpz_mpoly_t A,
    fmpz_mpoly_t M,
    const fmpz_mpoly_t S,
    const fmpz_mod_mpoly_bma_interpolate_ctx_t Ictx)
{
    fmpz_mod_mpoly_bma_interpolate_eval_useskelmul(eval, A, M, S, Ictx);
}


void fmpz_mpolyuu_print_pretty(
    const fmpz_mpolyu_t poly,
    const char ** x,
    const fmpz_mpoly_ctx_t ctx);

void _fmpz_mpolyuu_use_skel_mul(
    fmpz_mod_mpolyun_t E,
    const fmpz_mpolyu_t A,
    fmpz_mpolyu_t M,
    const fmpz_mpolyu_t S,
    const fmpz_mod_mpoly_bma_interpolate_ctx_t Ictx)
{
    slong xexp, yexp;
    slong i;
    slong N = mpoly_words_per_exp_sp(E->bits, Ictx->polyctx->minfo);
    fmpz_t eval;

    fmpz_init(eval);
/*
printf("use called **********\n");
printf("A: "); fmpz_mpolyuu_print_pretty(A, NULL, Ictx->polyctx); printf("\n");
*/
    E->length = 0;
    for (i = 0; i < A->length; i++)
    {
/*
flint_printf("A[%wd]: ",i); fmpz_mpoly_print_pretty(A->coeffs + i, NULL, Ictx->polyctx); printf("\n");
flint_printf("M[%wd]: ",i); fmpz_mpoly_print_pretty(M->coeffs + i, NULL, Ictx->polyctx); printf("\n");
flint_printf("S[%wd]: ",i); fmpz_mpoly_print_pretty(S->coeffs + i, NULL, Ictx->polyctx); printf("\n");
*/

        fmpz_mod_mpoly_bma_interpolate_eval_useskelmul(eval, A->coeffs + i, M->coeffs + i, S->coeffs + i, Ictx);
        if (fmpz_is_zero(eval))
        {
            continue;
        }

        xexp = A->exps[i] >> (FLINT_BITS/2);
        yexp = A->exps[i] & ((-UWORD(1)) >> (FLINT_BITS/2));
/*
flint_printf("xexp %wd  yexp %wd  eval "); fmpz_print(eval); printf("\n");
*/
        if (E->length > 0 && E->exps[E->length - 1] == xexp)
        {
            fmpz_mod_poly_set_coeff_fmpz(E->coeffs[E->length-1].coeffs + 0, yexp, eval);
        }
        else
        {
            fmpz_mod_mpolyun_fit_length(E, E->length + 1, Ictx->polyctx, Ictx->fpctx);
            fmpz_mod_mpolyn_fit_length(E->coeffs + E->length, 1, Ictx->polyctx, Ictx->fpctx);
            (E->coeffs + E->length)->length = 1;
            mpoly_monomial_zero((E->coeffs + E->length)->exps, N);
            fmpz_mod_poly_zero((E->coeffs + E->length)->coeffs + 0);
            fmpz_mod_poly_set_coeff_fmpz((E->coeffs + E->length)->coeffs + 0, yexp, eval);
            E->exps[E->length] = xexp;
            E->length++;
        }
    }

    fmpz_clear(eval);
}



int fmpz_mod_mpoly_bma_interpolate_get_mpoly(fmpz_mpoly_t A, const fmpz_t alphashift, fmpz_mod_mpoly_bma_interpolate_t I,
                                        const fmpz_mod_mpoly_bma_interpolate_ctx_t Ictx)
{
    slong i, j, t, N;
    int success;
    ulong this_exp;
    fmpz_t new_exp;
    slong * shifts, * offsets;
    fmpz * values, * roots;
    fmpz * Acoeff;
    ulong * Aexp;
    slong Alen;
    fmpz_t T, S, V, temp;
    TMP_INIT;

    TMP_START;

    fmpz_init(T);
    fmpz_init(S);
    fmpz_init(V);
    fmpz_init(temp);
    fmpz_init(new_exp);

    fmpz_mod_bma_reduce(I->bma);
    t = fmpz_mod_poly_degree(I->bma->V1);
    FLINT_ASSERT(t > 0);
    FLINT_ASSERT(I->evals->length >= t);

    fmpz_mod_poly_fit_length(I->roots, t);
    I->roots->length = t;
    success = _fmpz_mod_find_roots(I->roots->coeffs, I->bma->V1, t, Ictx->fpctx);
    if (!success)
    {
        goto cleanup;
    }

    roots = I->roots->coeffs;
    values = I->evals->coeffs;

    FLINT_ASSERT(Ictx->polyctx->minfo->ord == ORD_LEX);
    FLINT_ASSERT(A->bits <= FLINT_BITS);
    N = mpoly_words_per_exp_sp(A->bits, Ictx->polyctx->minfo);
    fmpz_mpoly_fit_length(A, t, Ictx->polyctx);
    Acoeff = A->coeffs;
    Aexp = A->exps;
    A->length = 0;

    shifts = (slong *) TMP_ALLOC(Ictx->polyctx->minfo->nvars);
    offsets = (slong *) TMP_ALLOC(Ictx->polyctx->minfo->nvars);
    for (j = 0; j < Ictx->polyctx->minfo->nvars; j++)
    {
        mpoly_gen_offset_shift_sp(offsets + j, shifts + j, j, A->bits, Ictx->polyctx->minfo);
    }

    Alen = 0;
    for (i = 0; i < t; i++)
    {
        /* coeffs[i] is (coeffs(P).values)/P(roots[i]) =: V/S
            where P(x) = master(x)/(x-roots[i])     */
        fmpz_zero(V);
        fmpz_zero(T);
        fmpz_zero(S);
        for (j = t; j > 0; j--)
        {
            fmpz_mod_mul(temp, roots + i, T, Ictx->fpctx);
            fmpz_mod_add(T, temp, I->bma->V1->coeffs + j, Ictx->fpctx);
            fmpz_mod_mul(temp, roots + i, S, Ictx->fpctx);
            fmpz_mod_add(S, temp, T, Ictx->fpctx);
            fmpz_mod_mul(temp, values + j - 1, T, Ictx->fpctx);
            fmpz_mod_add(V, V, temp, Ictx->fpctx);
        }
        /* roots[i] should be a root of master */
#if WANT_ASSERT
        fmpz_mod_mul(temp, roots + i, T, Ictx->fpctx);
        fmpz_mod_add(temp, temp, I->bma->V1->coeffs + 0, Ictx->fpctx);
        FLINT_ASSERT(fmpz_is_zero(temp));
#endif
        fmpz_mod_pow_fmpz(temp, roots + i, alphashift, Ictx->fpctx);
        fmpz_mod_mul(S, S, temp, Ictx->fpctx);
        fmpz_mod_inv(temp, S, Ictx->fpctx);
        fmpz_mod_mul(Acoeff + Alen, V, temp, Ictx->fpctx);
        
        if (fmpz_is_zero(Acoeff + Alen))
        {
            /* hmmm */
            continue;
        }

        mpoly_monomial_zero(Aexp + N*Alen, N);
        fmpz_mod_dlog_env_run(Ictx->dlogenv, new_exp, roots + i);
        for (j = Ictx->polyctx->minfo->nvars - 1; j >= 0; j--)
        {
            this_exp = fmpz_fdiv_ui(new_exp, Ictx->subdegs[j]);
            fmpz_fdiv_q_ui(new_exp, new_exp, Ictx->subdegs[j]);
            if (this_exp >= Ictx->degbounds[j])
            {
                success = 0;
                goto cleanup;
            }
            (Aexp + N*Alen)[offsets[j]] |= this_exp << shifts[j];
        }
        if (!fmpz_is_zero(new_exp))
        {
            success = 0;
            goto cleanup;
        }
        Alen++;
    }
    A->length = Alen;

    success = 1;

cleanup:

    fmpz_clear(T);
    fmpz_clear(S);
    fmpz_clear(V);
    fmpz_clear(temp);

    TMP_END;
    return success;
}



void fmpz_mpolyuu_print_pretty(
    const fmpz_mpolyu_t poly,
    const char ** x,
    const fmpz_mpoly_ctx_t ctx)
{
    slong i;

    if (poly->length == 0)
        flint_printf("0");

    for (i = 0; i < poly->length; i++)
    {
        if (i != 0)
            flint_printf(" + ");
        flint_printf("(");
        fmpz_mpoly_print_pretty(poly->coeffs + i, x, ctx);
        flint_printf(")*X^%wd*Y^%wd", poly->exps[i] >> (FLINT_BITS/2), 
                                      poly->exps[i] & ((-UWORD(1)) >> (FLINT_BITS/2)));
    }
}

fmpz_mpoly_struct * _fmpz_mpolyu_get_coeff(fmpz_mpolyu_t A,
                             ulong pow, const fmpz_mpoly_ctx_t uctx);

/*
    Convert B to A using the variable permutation perm.
    The uctx should be the context of the coefficients of A.
    The ctx should be the context of B.

    operation on each term:

    for 0 <= k < m + 2
        l = perm[k]
        Aexp[k] = (Bexp[l] - shift[l])/stride[l]

    the most significant main variable uses Aexp[0]
    the least significant main variable uses Aexp[1]
    the coefficients of A use variables Aexp[2], ..., Aexp[m + 1]
*/
void fmpz_mpoly_to_mpolyuu_perm_deflate(
    fmpz_mpolyu_t A,
    const slong * perm,
    const ulong * shift,
    const ulong * stride,
    const fmpz_mpoly_ctx_t uctx,
    const fmpz_mpoly_t B,
    const fmpz_mpoly_ctx_t ctx)
{
    slong i, j, k, l;
    slong n = ctx->minfo->nvars;
    slong m = uctx->minfo->nvars;
    slong NA, NB;
    ulong * uexps;
    ulong * Bexps;
    fmpz_mpoly_struct * Ac;
    TMP_INIT;

    FLINT_ASSERT(A->bits <= FLINT_BITS);
    FLINT_ASSERT(B->bits <= FLINT_BITS);
    FLINT_ASSERT(m + 2 <= n);

    TMP_START;

    uexps = (ulong *) TMP_ALLOC((m + 2)*sizeof(ulong));
    Bexps = (ulong *) TMP_ALLOC(n*sizeof(ulong));

    fmpz_mpolyu_zero(A, uctx);

    NA = mpoly_words_per_exp(A->bits, uctx->minfo);
    NB = mpoly_words_per_exp(B->bits, ctx->minfo);

    for (j = 0; j < B->length; j++)
    {
        mpoly_get_monomial_ui(Bexps, B->exps + NB*j, B->bits, ctx->minfo);
        for (k = 0; k < m + 2; k++)
        {
            l = perm[k];
            FLINT_ASSERT(stride[l] != UWORD(0));
            FLINT_ASSERT(((Bexps[l] - shift[l]) % stride[l]) == UWORD(0));
            uexps[k] = (Bexps[l] - shift[l]) / stride[l];
        }
        FLINT_ASSERT(FLINT_BIT_COUNT(uexps[0]) < FLINT_BITS/2);
        FLINT_ASSERT(FLINT_BIT_COUNT(uexps[1]) < FLINT_BITS/2);
        Ac = _fmpz_mpolyu_get_coeff(A, (uexps[0] << (FLINT_BITS/2)) + uexps[1], uctx);
        FLINT_ASSERT(Ac->bits == A->bits);

        fmpz_mpoly_fit_length(Ac, Ac->length + 1, uctx);
        fmpz_set(Ac->coeffs + Ac->length, B->coeffs + j);
        mpoly_set_monomial_ui(Ac->exps + NA*Ac->length, uexps + 2, A->bits, uctx->minfo);
        Ac->length++;
    }

    for (i = 0; i < A->length; i++)
    {
        fmpz_mpoly_sort_terms(A->coeffs + i, uctx);
    }

    TMP_END;
}


/*
    Convert B to A using the variable permutation vector perm.
    A must be constructed with bits = Abits.

    operation on each term:

        for 0 <= l < n
            Aexp[l] = shift[l]

        for 0 <= k <= m
            l = perm[k]
            Aexp[l] += scale[l]*Bexp[k]
*/
void fmpz_mpoly_from_mpolyuu_perm_inflate(
    fmpz_mpoly_t A,
    mp_bitcnt_t Abits,
    const fmpz_mpoly_ctx_t ctx,
    const fmpz_mpolyu_t B,
    const slong * perm,
    const ulong * shift,
    const ulong * stride,
    const fmpz_mpoly_ctx_t uctx)
{
    slong n = ctx->minfo->nvars;
    slong m = uctx->minfo->nvars;
    slong i, j, k, l;
    slong NA, NB;
    slong Alen;
    fmpz * Acoeff;
    ulong * Aexp;
    slong Aalloc;
    ulong * uexps;
    ulong * Aexps;
    TMP_INIT;

    FLINT_ASSERT(B->length > 0);
    FLINT_ASSERT(Abits <= FLINT_BITS);
    FLINT_ASSERT(B->bits <= FLINT_BITS);
    FLINT_ASSERT(m + 2 <= n);
    TMP_START;

    uexps = (ulong *) TMP_ALLOC((m + 2)*sizeof(ulong));
    Aexps = (ulong *) TMP_ALLOC(n*sizeof(ulong));

    NA = mpoly_words_per_exp(Abits, ctx->minfo);
    NB = mpoly_words_per_exp(B->bits, uctx->minfo);

    fmpz_mpoly_fit_bits(A, Abits, ctx);
    A->bits = Abits;

    Acoeff = A->coeffs;
    Aexp = A->exps;
    Aalloc = A->alloc;
    Alen = 0;
    for (i = 0; i < B->length; i++)
    {
        fmpz_mpoly_struct * Bc = B->coeffs + i;
        _fmpz_mpoly_fit_length(&Acoeff, &Aexp, &Aalloc, Alen + Bc->length, NA);
        FLINT_ASSERT(Bc->bits == B->bits);

        for (j = 0; j < Bc->length; j++)
        {
            fmpz_set(Acoeff + Alen + j, Bc->coeffs + j);
            mpoly_get_monomial_ui(uexps + 2, Bc->exps + NB*j, Bc->bits, uctx->minfo);
            uexps[0] = B->exps[i] >> (FLINT_BITS/2);
            uexps[1] = B->exps[i] & ((-UWORD(1)) >> (FLINT_BITS/2));
            for (l = 0; l < n; l++)
            {
                Aexps[l] = shift[l];
            }
            for (k = 0; k <= m; k++)
            {
                l = perm[k];
                Aexps[l] += stride[l]*uexps[k];
            }
            mpoly_set_monomial_ui(Aexp + NA*(Alen + j), Aexps, Abits, ctx->minfo);
        }
        Alen += Bc->length;
    }
    A->coeffs = Acoeff;
    A->exps = Aexp;
    A->alloc = Aalloc;
    _fmpz_mpoly_set_length(A, Alen, ctx);

    fmpz_mpoly_sort_terms(A, ctx);
    TMP_END;
}

/* return the degree of A wrt variable var */
slong fmpz_mpoly_eval_all_but_one_nmod(
    nmod_poly_t out,
    const fmpz_mpoly_t A,
    slong var,
    mp_limb_t * values,
    const fmpz_mpoly_ctx_t ctx)
{
    slong i, j;
    slong deg;
    ulong varexp, thisexp;
    mp_limb_t t, v;
    ulong mask;
    slong * offsets, * shifts;
    slong N = mpoly_words_per_exp_sp(A->bits, ctx->minfo);
    ulong * Aexp = A->exps;
    fmpz * Acoeff = A->coeffs;
    TMP_INIT;

    FLINT_ASSERT(A->bits < FLINT_BITS);

    TMP_START;

    mask = (-UWORD(1)) >> (FLINT_BITS - A->bits);
    offsets = (slong *) TMP_ALLOC(ctx->minfo->nvars*sizeof(slong));
    shifts = (slong *) TMP_ALLOC(ctx->minfo->nvars*sizeof(slong));
    for (j = 0; j < ctx->minfo->nvars; j++)
    {
        mpoly_gen_offset_shift_sp(offsets + j, shifts + j, j, A->bits, ctx->minfo);
    }

    nmod_poly_zero(out);
    deg = -WORD(1);
    for (i = 0; i < A->length; i++)
    {
        varexp = ((Aexp + N*i)[offsets[var]]>>shifts[var])&mask;
        deg = FLINT_MAX(deg, (slong)(varexp));
        v = fmpz_fdiv_ui(Acoeff + i, out->mod.n);
        for (j = 0; j < ctx->minfo->nvars; j++)
        {
            thisexp = ((Aexp + N*i)[offsets[j]]>>shifts[j])&mask;
            if (j != var)
            {
                v = nmod_mul(v, nmod_pow_ui(values[j], thisexp, out->mod), out->mod);
            }
        }
        t = nmod_poly_get_coeff_ui(out, varexp);
        nmod_poly_set_coeff_ui(out, varexp, nmod_add(t, v, out->mod));
    }

    TMP_END;
    return deg;
}

slong fmpz_mpolyuu_eval_all_but_one_nmod(
    nmod_poly_t out,
    const fmpz_mpolyu_t A,
    slong var,
    mp_limb_t * values,
    const fmpz_mpoly_ctx_t ctx)
{
    slong i;
    slong deg, thisdeg;
    ulong * Aexp = A->exps;
    const fmpz_mpoly_struct * Acoeff = A->coeffs;
    nmod_poly_t temp;
    mp_limb_t t, t1;

    nmod_poly_zero(out);
    nmod_poly_init(temp, out->mod.n);

    deg = -WORD(1);
    for (i = 0; i < A->length; i++)
    {
        t = nmod_pow_ui(values[0], Aexp[i] >> (FLINT_BITS/2), out->mod);
        t1 = nmod_pow_ui(values[1], Aexp[i] & ((-UWORD(1)) >> (FLINT_BITS/2)), out->mod);
        t = nmod_mul(t, t1, out->mod);
        thisdeg = fmpz_mpoly_eval_all_but_one_nmod(temp, Acoeff + i, var, values + 2, ctx);
        deg = FLINT_MAX(deg, thisdeg);
        nmod_poly_scalar_mul_nmod(temp, temp, t);
        nmod_poly_add(out, out, temp);
    }

    nmod_poly_clear(temp);   
    return deg; 
}


slong _degree_bound(
    slong * Adeg,
    slong * Bdeg,
    const fmpz_mpolyu_t A,
    const fmpz_mpolyu_t B,
    slong var,
    const fmpz_mpoly_ctx_t ctx,
    flint_rand_t state)
{
    slong i;
    int tries = 0;
    slong degA, degB, degRet;
    mp_limb_t p = UWORD(1) << (FLINT_BITS - 2);
    nmod_poly_t Geval, Aeval, Beval;
    mp_limb_t * values;
    TMP_INIT;

    TMP_START;

    values = (mp_limb_t *) TMP_ALLOC(ctx->minfo->nvars*sizeof(mp_limb_t));

    p = n_nextprime(p, 1);
    nmod_poly_init(Geval, p);
    nmod_poly_init(Aeval, p);
    nmod_poly_init(Beval, p);

try_again:

usleep(100000);

    for (i = 0; i < ctx->minfo->nvars + 2; i++)
    {
        values[i] = n_urandint(state, p);
    }

    degA = fmpz_mpolyuu_eval_all_but_one_nmod(Aeval, A, var, values, ctx);
    degB = fmpz_mpolyuu_eval_all_but_one_nmod(Beval, B, var, values, ctx);
    *Adeg = degA;
    *Bdeg = degB;
/*
flint_printf("degA: %wd\n", degA);
flint_printf("degB: %wd\n", degB);
printf("Aeval: "); nmod_poly_print_pretty(Aeval, "v"); printf("\n");
printf("Beval: "); nmod_poly_print_pretty(Beval, "v"); printf("\n");
*/
    if (degA != nmod_poly_degree(Aeval) || degB != nmod_poly_degree(Beval))
    {
        if (++tries > 100)
        {
            degRet = FLINT_MIN(degA, degB);
            goto cleanup;
        }
        p = n_nextprime(p, 1);
        nmod_poly_clear(Geval);
        nmod_poly_clear(Aeval);
        nmod_poly_clear(Beval);
        nmod_poly_init(Geval, p);
        nmod_poly_init(Aeval, p);
        nmod_poly_init(Beval, p);
        goto try_again;
    }

    nmod_poly_gcd(Geval, Aeval, Beval);
    degRet = nmod_poly_degree(Geval);

cleanup:

    nmod_poly_clear(Geval);
    nmod_poly_clear(Aeval);
    nmod_poly_clear(Beval);
    TMP_END;
    return degRet;
}


void _ksub_content(
    fmpz_t content,
    const fmpz_mpoly_t A,
    const fmpz_mod_mpoly_bma_interpolate_ctx_t Ictx)
{
    slong i, j;
    slong N = mpoly_words_per_exp_sp(A->bits, Ictx->polyctx->minfo);
    fmpz_mpoly_t T;
    fmpz_mpoly_ctx_t Tctx;
    fmpz_t exp;
    fmpz * Acoeff = A->coeffs;
    ulong * Aexp = A->exps;
    ulong mask;
    slong * offsets, * shifts;
    TMP_INIT;

    TMP_START;

    fmpz_mpoly_ctx_init(Tctx, 1, ORD_LEX);
    fmpz_mpoly_init(T, Tctx);

    mask = (-UWORD(1)) >> (FLINT_BITS - A->bits);
    offsets = (slong *) TMP_ALLOC(Ictx->polyctx->minfo->nvars*sizeof(slong));
    shifts = (slong *) TMP_ALLOC(Ictx->polyctx->minfo->nvars*sizeof(slong));
    for (j = 0; j < Ictx->polyctx->minfo->nvars; j++)
    {
        mpoly_gen_offset_shift_sp(offsets + j, shifts + j, j, A->bits, Ictx->polyctx->minfo);
    }

    for (i = 0; i < A->length; i++)
    {
        fmpz_zero(exp);
        for (j = 0; j < Ictx->polyctx->minfo->nvars; j++)
        {
            fmpz_mul_ui(exp, exp, Ictx->subdegs[j]);
            fmpz_add_ui(exp, exp, ((Aexp + N*i)[offsets[j]]>>shifts[j])&mask);
        }
        _fmpz_mpoly_push_exp_ffmpz(T, exp, Tctx);
        fmpz_set(T->coeffs + T->length - 1, Acoeff + i);
    }
    fmpz_mpoly_sort_terms(T, Tctx);
    fmpz_mpoly_combine_like_terms(T, Tctx);

printf("after ksub: "); fmpz_mpoly_print_pretty(T, NULL, Tctx); printf("\n");

    _fmpz_vec_content(content, T->coeffs, T->length);
printf("content: "); fmpz_print(content); printf("\n");


    fmpz_mpoly_clear(T, Tctx);
    fmpz_mpoly_ctx_clear(Tctx);
    TMP_END;
}


typedef struct {
    fmpz_mod_mpoly_bma_interpolate_struct * coeffs;
    ulong * exps;
    slong length;
    slong alloc;
    slong pointcount;
} fmpz_mod_bma_mpoly_struct;

typedef fmpz_mod_bma_mpoly_struct fmpz_mod_bma_mpoly_t[1];

void fmpz_mod_bma_mpoly_init(
    fmpz_mod_bma_mpoly_t A,
    const fmpz_mod_mpoly_bma_interpolate_ctx_t Ictx)
{
    A->length = 0;
    A->alloc = 0;
    A->exps = NULL;
    A->coeffs = NULL;
    A->pointcount = 0;
}

void fmpz_mod_bma_mpoly_reset_prime(
    fmpz_mod_bma_mpoly_t A,
    const fmpz_mod_mpoly_bma_interpolate_ctx_t Ictx)
{
    slong i;
    for (i = 0; i < A->alloc; i++)
    {
        fmpz_mod_mpoly_bma_interpolate_reset_prime(A->coeffs + i, Ictx->fpctx->p);
    }
}


void fmpz_mod_bma_mpoly_clear(
    fmpz_mod_bma_mpoly_t A,
    const fmpz_mod_mpoly_bma_interpolate_ctx_t Ictx)
{
    slong i;
    for (i = 0; i < A->alloc; i++)
    {
        fmpz_mod_mpoly_bma_interpolate_clear(A->coeffs + i);
    }

    if (A->exps)
        flint_free(A->exps);
    if (A->coeffs)
        flint_free(A->coeffs);
}

void fmpz_mod_bma_mpoly_print(
    fmpz_mod_bma_mpoly_t A,
    const fmpz_mod_mpoly_bma_interpolate_ctx_t Ictx)
{
    slong i;
    flint_printf("0");
    for (i = 0; i < A->length; i++)
    {
        flint_printf(" + [");
        fmpz_mod_mpoly_bma_interpolate_print(A->coeffs + i);
        flint_printf("]*X^%wd*Y^%wd", A->exps[i] >> (FLINT_BITS/2), A->exps[i] & ((-UWORD(1)) >> (FLINT_BITS/2)));
    }
}


void fmpz_mod_bma_mpoly_fit_length(
    fmpz_mod_bma_mpoly_t A,
    slong length,
    const fmpz_mod_mpoly_bma_interpolate_ctx_t Ictx)
{
    slong i;
    slong old_alloc = A->alloc;
    slong new_alloc = FLINT_MAX(length, A->alloc + A->alloc/2);

    if (length > old_alloc)
    {
        if (old_alloc == 0)
        {
            A->exps = (ulong *) flint_malloc(new_alloc*sizeof(ulong));
            A->coeffs = (fmpz_mod_mpoly_bma_interpolate_struct *) flint_malloc(new_alloc*sizeof(fmpz_mod_mpoly_bma_interpolate_struct));
        }
        else
        {
            A->exps = (ulong *) flint_realloc(A->exps, new_alloc*sizeof(ulong));
            A->coeffs = (fmpz_mod_mpoly_bma_interpolate_struct *) flint_realloc(A->coeffs, new_alloc*sizeof(fmpz_mod_mpoly_bma_interpolate_struct));
        }

        for (i = old_alloc; i < new_alloc; i++)
        {
            fmpz_mod_mpoly_bma_interpolate_init(A->coeffs + i, Ictx->fpctx->p);
        }
        A->alloc = new_alloc;
    }
}

void fmpz_mod_bma_mpoly_zero(
    fmpz_mod_bma_mpoly_t L,
    const fmpz_mod_mpoly_bma_interpolate_ctx_t Ictx)
{
    L->length = 0;
    L->pointcount = 0;
}

int fmpz_mod_bma_mpoly_reduce(
    fmpz_mod_bma_mpoly_t L,
    const fmpz_mod_mpoly_bma_interpolate_ctx_t Ictx)
{
    slong i;
    int changed;

    changed = 0;

    for (i = 0; i < L->length; i++)
    {
        FLINT_ASSERT(L->pointcount == fmpz_mod_mpoly_bma_interpolate_pointcount(L->coeffs + i));
        changed |= fmpz_mod_mpoly_bma_interpolate_reduce(L->coeffs + i);
    }

    return changed;
}

void fmpz_mod_bma_mpoly_add_point(
    fmpz_mod_bma_mpoly_t L,
    const fmpz_mod_mpolyun_t A,
    const fmpz_mod_mpoly_bma_interpolate_ctx_t Ictx)
{
    slong j;
    slong Alen = A->length;
    fmpz_mod_mpolyn_struct * Acoeff = A->coeffs;
    slong Li, Ai, ai;
    ulong Aexp;
    fmpz_mod_mpoly_bma_interpolate_struct * Lcoeff;
    slong Llen;
    ulong * Lexp;
    fmpz_t zero;

    fmpz_init(zero);
/*
printf("inside add_point\n");
*/
    if (L->length == 0)
    {
        slong tot = 0;
        for (Ai = 0; Ai < Alen; Ai++)
        {
            for (ai = fmpz_mod_poly_degree((Acoeff + Ai)->coeffs + 0); ai >= 0; ai--)
            {
                tot += !fmpz_is_zero((Acoeff + Ai)->coeffs[0].coeffs + ai);
            }
        }
flint_printf("fitting to length %wd\n", tot);
        fmpz_mod_bma_mpoly_fit_length(L, tot, Ictx);
    }

    Lcoeff = L->coeffs;
    Llen = L->length;
    Lexp = L->exps;

    Li = 0;
    Ai = ai = 0;
    Aexp = 0;
    if (Ai < Alen)
    {
        ai = fmpz_mod_poly_degree((A->coeffs + Ai)->coeffs + 0);
        Aexp = (A->exps[Ai] << (FLINT_BITS/2)) + ai;
    }

    while (Li < Llen || Ai < Alen)
    {
/*
flint_printf("Li: %wd, Ai: %wd\n", Li, Ai);
usleep(1000000);
*/
        if (Li < Llen && Ai < Alen && Lexp[Li] == Aexp)
        {
add_same_exp:
            fmpz_mod_mpoly_bma_interpolate_add_point(Lcoeff + Li, (Acoeff + Ai)->coeffs[0].coeffs + ai);
            Li++;
            do {
                ai--;
            } while (ai >= 0 && fmpz_is_zero((Acoeff + Ai)->coeffs[0].coeffs + ai));
            if (ai < 0)
            {
                Ai++;
                if (Ai < Alen)
                {
                    ai = fmpz_mod_poly_degree((A->coeffs + Ai)->coeffs + 0);        
                    Aexp = (A->exps[Ai] << (FLINT_BITS/2)) + ai;
                }
            }
            if (Ai < Alen)
            {
                Aexp = (A->exps[Ai] << (FLINT_BITS/2)) + ai;
            }
        }
        else if (Li < Llen && (Ai >= Alen || Lexp[Li] > Aexp))
        {
            fmpz_mod_mpoly_bma_interpolate_add_point(Lcoeff + Li, zero);            
        }
        else
        {
            FLINT_ASSERT(Li >= Llen || (Ai < Alen && Lexp[Li] < Aexp));
            {
                ulong texp;
                fmpz_mod_mpoly_bma_interpolate_struct tcoeff;

                fmpz_mod_bma_mpoly_fit_length(L, Llen + 1, Ictx);
                Lcoeff = L->coeffs;
                Lexp = L->exps;

                texp = Lexp[Llen];
                tcoeff = Lcoeff[Llen];
                for (j = Li; j < Llen; j++)
                {
                    Lexp[j + 1] = Lexp[j];
                    Lcoeff[j + 1] = Lcoeff[j];
                }
                Lexp[Li] = texp;
                Lcoeff[Li] = tcoeff;
            }

            fmpz_mod_mpoly_bma_interpolate_zero(Lcoeff + Li, L->pointcount);
            Lexp[Li] = Aexp;
            Llen++;
            L->length = Llen;
            goto add_same_exp;
        }
    }

    fmpz_clear(zero);
    L->pointcount++;
}

int fmpz_mod_bma_mpoly_get_mpolyu(
    fmpz_mpolyu_t A,
    const fmpz_t alphashift,
    const fmpz_mod_bma_mpoly_t L,
    const fmpz_mod_mpoly_bma_interpolate_ctx_t Ictx)
{
    int success;
    slong i;

    fmpz_mpolyu_fit_length(A, L->length, Ictx->polyctx);
    A->length = 0;
    for (i = 0; i < L->length; i++)
    {
        A->exps[A->length] = L->exps[i];
        success = fmpz_mod_mpoly_bma_interpolate_get_mpoly(A->coeffs + A->length, alphashift, L->coeffs + i, Ictx);
        if (!success)
        {
            return 0;
        }
        A->length += !fmpz_mpoly_is_zero(A->coeffs + A->length, Ictx->polyctx);
    }
    return 1;
}




/*
    A = B
    A, B are in R[X]
*/
void fmpz_mod_mpolyun_set_poly(
    fmpz_mod_mpolyun_t A,
    const fmpz_mod_poly_t B,
    const fmpz_mpoly_ctx_t ctx,
    const fmpz_mod_ctx_t fpctx)
{
    slong N = mpoly_words_per_exp_sp(A->bits, ctx->minfo);
    slong Bexp;
    slong Blen = fmpz_mod_poly_length(B);
    fmpz * Bcoeff = B->coeffs;
    fmpz_mod_mpolyn_struct * Acoeff;
    ulong * Aexp;
    slong Ai;

    fmpz_mod_mpolyun_fit_length(A, Blen, ctx, fpctx);
    Acoeff = A->coeffs;
    Aexp = A->exps;

    Ai = 0;
    for (Bexp = Blen - 1; Bexp >= 0; Bexp--)
    {
        if (!fmpz_is_zero(Bcoeff + Bexp))
        {
            FLINT_ASSERT(Ai < A->alloc);

            fmpz_mod_mpolyn_fit_length(Acoeff + Ai, 1, ctx, fpctx);
            mpoly_monomial_zero((Acoeff + Ai)->exps + N*0, N);
            fmpz_mod_poly_zero((Acoeff + Ai)->coeffs + 0);
            fmpz_mod_poly_set_coeff_fmpz((Acoeff + Ai)->coeffs + 0, 0, Bcoeff + Bexp);
            Aexp[Ai] = Bexp;
            (Acoeff + Ai)->length = 1;
            Ai++;
        }
    }
    A->length = Ai;
}


/*
    F = F + modulus*(A - F(v = alpha))
    no assumptions about matching monomials
    F is in Fp[X][v]
    A is in Fp[X]
    it is expected that modulus(alpha) == 1
*/
int fmpz_mod_mpolyun_addinterp_bivar(
    slong * lastdeg_,
    fmpz_mod_mpolyun_t F,
    fmpz_mod_mpolyun_t T,
    const fmpz_mod_poly_t A,
    const fmpz_mod_poly_t modulus,
    const fmpz_t alpha,
    const fmpz_mpoly_ctx_t ctx,
    const fmpz_mod_ctx_t fpctx)
{
    int changed = 0;
    slong lastdeg = -WORD(1);
    slong N = mpoly_words_per_exp_sp(F->bits, ctx->minfo);
    fmpz_t v;
    slong Fi, Toff, Aexp;
    fmpz * Acoeff = A->coeffs;
    slong Flen = F->length;
    fmpz_mod_mpolyn_struct * Fcoeff = F->coeffs;
    ulong * Fexp = F->exps;
    fmpz_mod_mpolyn_struct * Tcoeff;
    ulong * Texp;
    fmpz_mod_poly_t tp;
    
    Fi = 0;

    fmpz_init(v);

/*
flint_printf("\nnmod_mpolyun_addinterp_bivar: (alpha = %wu)\n", alpha);
flint_printf("A: "); nmod_poly_print_pretty(A, "X"); printf("\n");
flint_printf("F: "); nmod_mpolyun_print_pretty(F, NULL, ctx); printf("\n");
*/

    Aexp = fmpz_mod_poly_degree(A);

    fmpz_mod_poly_init(tp, fpctx->p);

    fmpz_mod_mpolyun_fit_length(T, Flen + Aexp + 1, ctx, fpctx);
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
            FLINT_ASSERT(!fmpz_mod_poly_is_zero((Fcoeff + Fi)->coeffs + 0));
            FLINT_ASSERT(fmpz_mod_poly_degree((Fcoeff + Fi)->coeffs + 0) < fmpz_mod_poly_degree(modulus));
        }

        if (Aexp >= 0)
        {
            FLINT_ASSERT(Acoeff[Aexp] != UWORD(0));
        }

        fmpz_mod_mpolyn_fit_length(Tcoeff + Toff, 1, ctx, fpctx);

        if (Fi < Flen && Aexp >= 0 && Fexp[Fi] == Aexp)
        {
            /* F term ok, A term ok */
            fmpz_mod_poly_evaluate_fmpz(v, (Fcoeff + Fi)->coeffs + 0, alpha);
            fmpz_mod_sub(v, Acoeff + Aexp, v, fpctx);
            if (!fmpz_is_zero(v))
            {
                changed = 1;
                fmpz_mod_poly_scalar_mul_fmpz(tp, modulus, v);
                fmpz_mod_poly_add((Tcoeff + Toff)->coeffs + 0, (Fcoeff + Fi)->coeffs + 0, tp);
            }
            else
            {
                fmpz_mod_poly_set((Tcoeff + Toff)->coeffs + 0, (Fcoeff + Fi)->coeffs + 0);
            }
            Texp[Toff] = Aexp;
            Fi++;
            do {
                Aexp--;
            } while (Aexp >= 0 && fmpz_is_zero(Acoeff + Aexp));
        }
        else if (Fi < Flen && (Aexp < 0 || Fexp[Fi] > Aexp))
        {
            /* F term ok, A term missing */
            fmpz_mod_poly_evaluate_fmpz(v, (Fcoeff + Fi)->coeffs + 0, alpha);
            if (!fmpz_is_zero(v))
            {
                changed = 1;
                fmpz_mod_poly_scalar_mul_fmpz(tp, modulus, v);
                fmpz_mod_poly_sub((Tcoeff + Toff)->coeffs + 0, (Fcoeff + Fi)->coeffs + 0, tp);
            }
            else
            {
                fmpz_mod_poly_set((Tcoeff + Toff)->coeffs + 0, (Fcoeff + Fi)->coeffs + 0);                
            }
            Texp[Toff] = Fexp[Fi];
            Fi++;
        }
        else if (Aexp >= 0 && (Fi >= Flen || Fexp[Fi] < Aexp))
        {
            /* F term missing, A term ok */
            changed = 1;
            fmpz_mod_poly_scalar_mul_fmpz((Tcoeff + Toff)->coeffs + 0, modulus, Acoeff + Aexp);
            Texp[Toff] = Aexp;
            do {
                Aexp--;
            } while (Aexp >= 0 && fmpz_is_zero(Acoeff + Aexp));
        }
        else
        {
            FLINT_ASSERT(0);
        }

        lastdeg = FLINT_MAX(lastdeg, fmpz_mod_poly_degree((Tcoeff + Toff)->coeffs + 0));
        mpoly_monomial_zero((Tcoeff + Toff)->exps + N*0, N);
        FLINT_ASSERT(!fmpz_mod_poly_is_zero((Tcoeff + Toff)->coeffs + 0));
        (Tcoeff + Toff)->length = 1;
        Toff++;
    }
    T->length = Toff;

    fmpz_mod_poly_clear(tp);

    if (changed)
    {
        fmpz_mod_mpolyun_swap(T, F);
    }

    fmpz_clear(v);

    *lastdeg_ = lastdeg;
    return changed;
}


int fmpz_mod_mpolyun_gcd_bivar(
    fmpz_mod_mpolyun_t G,
    fmpz_mod_mpolyun_t Abar,
    fmpz_mod_mpolyun_t Bbar,
    fmpz_mod_mpolyun_t A,
    fmpz_mod_mpolyun_t B,
    const fmpz_mpoly_ctx_t ctx,
    const fmpz_mod_ctx_t fpctx)
{
    int success;
    slong bound;
    fmpz_t alpha, temp, gammaeval;
    fmpz_mod_poly_t Aeval, Beval, Geval, Abareval, Bbareval;
    fmpz_mod_mpolyun_t T;
    slong deggamma, ldegG, ldegAbar, ldegBbar, ldegA, ldegB;
    fmpz_mod_poly_t cA, cB, cG, cAbar, cBbar, gamma, r;
    fmpz_mod_poly_t modulus, modulus2;
#if WANT_ASSERT
    fmpz_mod_poly_t leadA, leadB;
#endif

    fmpz_init(alpha);
    fmpz_init(temp);
    fmpz_init(gammaeval);
/*
printf("fmpz_mod_mpolyun_gcd_bivar called\n");
printf("A: "); fmpz_mod_mpolyun_print_pretty(A, NULL, ctx, fpctx); printf("\n");
printf("B: "); fmpz_mod_mpolyun_print_pretty(B, NULL, ctx, fpctx); printf("\n");
*/

#if WANT_ASSERT
    fmpz_mod_poly_init(leadA, fpctx->p);
    fmpz_mod_poly_init(leadB, fpctx->p);
    fmpz_mod_poly_set(leadA, fmpz_mod_mpolyun_leadcoeff_ref(A, ctx, fpctx));
    fmpz_mod_poly_set(leadB, fmpz_mod_mpolyun_leadcoeff_ref(B, ctx, fpctx));
#endif

    fmpz_mod_poly_init(r, fpctx->p);
    fmpz_mod_poly_init(cA, fpctx->p);
    fmpz_mod_poly_init(cB, fpctx->p);
    fmpz_mod_mpolyun_content_last(cA, A, ctx, fpctx);
    fmpz_mod_mpolyun_content_last(cB, B, ctx, fpctx);
    fmpz_mod_mpolyun_divexact_last(A, cA, ctx, fpctx);
    fmpz_mod_mpolyun_divexact_last(B, cB, ctx, fpctx);

    fmpz_mod_poly_init(cG, fpctx->p);
    fmpz_mod_poly_gcd_euclidean(cG, cA, cB);

    fmpz_mod_poly_init(cAbar, fpctx->p);
    fmpz_mod_poly_init(cBbar, fpctx->p);
    fmpz_mod_poly_divrem(cAbar, r, cA, cG);
    FLINT_ASSERT(fmpz_mod_poly_is_zero(r));
    fmpz_mod_poly_divrem(cBbar, r, cB, cG);
    FLINT_ASSERT(fmpz_mod_poly_is_zero(r));

    fmpz_mod_poly_init(gamma, fpctx->p);
    fmpz_mod_poly_gcd(gamma, fmpz_mod_mpolyun_leadcoeff_ref(A, ctx, fpctx),
                             fmpz_mod_mpolyun_leadcoeff_ref(B, ctx, fpctx));

    ldegA = fmpz_mod_mpolyun_lastdeg(A, ctx, fpctx);
    ldegB = fmpz_mod_mpolyun_lastdeg(B, ctx, fpctx);
    deggamma = fmpz_mod_poly_degree(gamma);
    bound = 1 + deggamma + FLINT_MAX(ldegA, ldegB);
/*
flint_printf("bound: %wd\n", bound);
*/
    fmpz_mod_poly_init(Aeval, fpctx->p);
    fmpz_mod_poly_init(Beval, fpctx->p);
    fmpz_mod_poly_init(Geval, fpctx->p);
    fmpz_mod_poly_init(Abareval, fpctx->p);
    fmpz_mod_poly_init(Bbareval, fpctx->p);

    fmpz_mod_mpolyun_init(T, A->bits, ctx, fpctx);

    fmpz_mod_poly_init(modulus, fpctx->p);
    fmpz_mod_poly_init(modulus2, fpctx->p);
    fmpz_mod_poly_set_ui(modulus, 1);

    fmpz_sub_ui(alpha, fpctx->p, 2);

choose_prime: /* prime is v - alpha */

    if (fmpz_cmp_ui(alpha, 2) < 1)
    {
        success = 0;
        goto cleanup;
    }

    fmpz_sub_ui(alpha, alpha, 1);

/*
printf("alpha: "); fmpz_print(alpha); printf("\n");
*/

    /* make sure evaluation point does not kill both lc(A) and lc(B) */
    fmpz_mod_poly_evaluate_fmpz(gammaeval, gamma, alpha);
    if (fmpz_is_zero(gammaeval))
    {
        goto choose_prime;
    }

    /* evaluation point should kill neither A nor B */
    fmpz_mod_mpolyun_eval_last_bivar(Aeval, A, alpha, ctx, fpctx);
    fmpz_mod_mpolyun_eval_last_bivar(Beval, B, alpha, ctx, fpctx);
    FLINT_ASSERT(Aeval->length > 0);
    FLINT_ASSERT(Beval->length > 0);
/*
printf("Aeval: "); fmpz_mod_poly_print_pretty(Aeval, "X"); printf("\n");
printf("Beval: "); fmpz_mod_poly_print_pretty(Beval, "X"); printf("\n");
*/
    fmpz_mod_poly_gcd(Geval, Aeval, Beval);
    fmpz_mod_poly_divrem(Abareval, r, Aeval, Geval);
    FLINT_ASSERT(fmpz_mod_poly_is_zero(r));
    fmpz_mod_poly_divrem(Bbareval, r, Beval, Geval);
    FLINT_ASSERT(fmpz_mod_poly_is_zero(r));

/*
printf("   Geval: "); fmpz_mod_poly_print_pretty(Geval, "X"); printf("\n");
printf("Abareval: "); fmpz_mod_poly_print_pretty(Abareval, "X"); printf("\n");
printf("Bbareval: "); fmpz_mod_poly_print_pretty(Bbareval, "X"); printf("\n");
*/

    FLINT_ASSERT(Geval->length > 0);
    FLINT_ASSERT(Abareval->length > 0);
    FLINT_ASSERT(Bbareval->length > 0);

    if (fmpz_mod_poly_degree(Geval) == 0)
    {
        fmpz_mod_mpolyun_one(G, ctx, fpctx);
        fmpz_mod_mpolyun_swap(Abar, A);
        fmpz_mod_mpolyun_swap(Bbar, B);
        goto successful_put_content;    
    }

    if (fmpz_mod_poly_degree(modulus) > 0)
    {
        FLINT_ASSERT(G->length > 0);
        if (fmpz_mod_poly_degree(Geval) > G->exps[0])
        {
            goto choose_prime;
        }
        else if (fmpz_mod_poly_degree(Geval) < G->exps[0])
        {
            fmpz_mod_poly_set_ui(modulus, 1);
        }
    }

    fmpz_mod_poly_scalar_mul_fmpz(Geval, Geval, gammaeval);

    if (fmpz_mod_poly_degree(modulus) > 0)
    {
        fmpz_mod_poly_evaluate_fmpz(temp, modulus, alpha);
        fmpz_mod_inv(temp, temp, fpctx);
        fmpz_mod_poly_scalar_mul_fmpz(modulus, modulus, temp);
        fmpz_mod_mpolyun_addinterp_bivar(&ldegG, G, T, Geval, modulus, alpha, ctx, fpctx);
        fmpz_mod_mpolyun_addinterp_bivar(&ldegAbar, Abar, T, Abareval, modulus, alpha, ctx, fpctx);
        fmpz_mod_mpolyun_addinterp_bivar(&ldegBbar, Bbar, T, Bbareval, modulus, alpha, ctx, fpctx);
    }
    else
    {
        fmpz_mod_mpolyun_set_poly(G, Geval, ctx, fpctx);
        fmpz_mod_mpolyun_set_poly(Abar, Abareval, ctx, fpctx);
        fmpz_mod_mpolyun_set_poly(Bbar, Bbareval, ctx, fpctx);
        ldegG = 0;
        ldegAbar = 0;
        ldegBbar = 0;
    }

    fmpz_mod_poly_scalar_mul_fmpz(modulus2, modulus, alpha);
    fmpz_mod_poly_shift_left(modulus, modulus, 1);
    fmpz_mod_poly_sub(modulus, modulus, modulus2);
/*
printf("modulus: "); fmpz_mod_poly_print_pretty(modulus, "v"); printf("\n");
printf("   G: "); fmpz_mod_mpolyun_print_pretty(G, NULL, ctx, fpctx); printf("\n");
printf("Abar: "); fmpz_mod_mpolyun_print_pretty(Abar, NULL, ctx, fpctx); printf("\n");
printf("Bbar: "); fmpz_mod_mpolyun_print_pretty(Bbar, NULL, ctx, fpctx); printf("\n");
*/
    if (fmpz_mod_poly_degree(modulus) < bound)
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

    fmpz_mod_poly_set_ui(modulus, 1);
    goto choose_prime;

successful:
/*
printf("successful\n");
printf("   G: "); fmpz_mod_mpolyun_print_pretty(G, NULL, ctx, fpctx); printf("\n");
printf("Abar: "); fmpz_mod_mpolyun_print_pretty(Abar, NULL, ctx, fpctx); printf("\n");
printf("Bbar: "); fmpz_mod_mpolyun_print_pretty(Bbar, NULL, ctx, fpctx); printf("\n");
*/

    fmpz_mod_mpolyun_content_last(modulus, G, ctx, fpctx);
    fmpz_mod_mpolyun_divexact_last(G, modulus, ctx, fpctx);
    fmpz_mod_mpolyun_divexact_last(Abar, fmpz_mod_mpolyun_leadcoeff_ref(G, ctx, fpctx), ctx, fpctx);
    fmpz_mod_mpolyun_divexact_last(Bbar, fmpz_mod_mpolyun_leadcoeff_ref(G, ctx, fpctx), ctx, fpctx);

successful_put_content:
/*
printf("put_content\n");
printf("   G: "); fmpz_mod_mpolyun_print_pretty(G, NULL, ctx, fpctx); printf("\n");
printf("Abar: "); fmpz_mod_mpolyun_print_pretty(Abar, NULL, ctx, fpctx); printf("\n");
printf("Bbar: "); fmpz_mod_mpolyun_print_pretty(Bbar, NULL, ctx, fpctx); printf("\n");
*/

    fmpz_mod_mpolyun_mul_last(G, cG, ctx, fpctx);
    fmpz_mod_mpolyun_mul_last(Abar, cAbar, ctx, fpctx);
    fmpz_mod_mpolyun_mul_last(Bbar, cBbar, ctx, fpctx);

    success = 1;

cleanup:

/*
printf("cleanup\n");
printf("   G: "); fmpz_mod_mpolyun_print_pretty(G, NULL, ctx, fpctx); printf("\n");
printf("Abar: "); fmpz_mod_mpolyun_print_pretty(Abar, NULL, ctx, fpctx); printf("\n");
printf("Bbar: "); fmpz_mod_mpolyun_print_pretty(Bbar, NULL, ctx, fpctx); printf("\n");
*/

#if WANT_ASSERT
    if (success)
    {
        FLINT_ASSERT(fmpz_is_one(fmpz_mod_mpolyun_leadcoeff_last_ref(G, ctx, fpctx)));
        fmpz_mod_poly_mul(modulus, fmpz_mod_mpolyun_leadcoeff_ref(G, ctx, fpctx),
                                   fmpz_mod_mpolyun_leadcoeff_ref(Abar, ctx, fpctx));
        FLINT_ASSERT(fmpz_mod_poly_equal(modulus, leadA));
        fmpz_mod_poly_mul(modulus, fmpz_mod_mpolyun_leadcoeff_ref(G, ctx, fpctx),
                                   fmpz_mod_mpolyun_leadcoeff_ref(Bbar, ctx, fpctx));
        FLINT_ASSERT(fmpz_mod_poly_equal(modulus, leadB));
    }
    fmpz_mod_poly_clear(leadA);
    fmpz_mod_poly_clear(leadB);
#endif

    fmpz_mod_poly_clear(r);
    fmpz_mod_poly_clear(cA);
    fmpz_mod_poly_clear(cB);
    fmpz_mod_poly_clear(cG);
    fmpz_mod_poly_clear(cAbar);
    fmpz_mod_poly_clear(cBbar);

    fmpz_mod_poly_clear(gamma);

    fmpz_mod_poly_clear(Aeval);
    fmpz_mod_poly_clear(Beval);
    fmpz_mod_poly_clear(Geval);
    fmpz_mod_poly_clear(Abareval);
    fmpz_mod_poly_clear(Bbareval);

    fmpz_mod_mpolyun_clear(T, ctx, fpctx);

    fmpz_mod_poly_clear(modulus);
    fmpz_mod_poly_clear(modulus2);

    return success;
}


ulong fmpz_mod_mpolyun_bidegree(
    const fmpz_mod_mpolyun_t A,
    const fmpz_mpoly_ctx_t ctx,
    const fmpz_mod_ctx_t fpctx)
{
    ulong degx, degy;

    FLINT_ASSERT(A->length > 0);

    degx = A->exps[0];
    degy = fmpz_mod_poly_degree((A->coeffs + 0)->coeffs + 0);

    FLINT_ASSERT(FLINT_BIT_COUNT(degx) < FLINT_BITS/2);
    FLINT_ASSERT(FLINT_BIT_COUNT(degy) < FLINT_BITS/2);

    return (degx << (FLINT_BITS/2)) + degy;
}


void fmpz_mpoly_eval_fmpz_mod(
    fmpz_t eval,
    const fmpz_mpoly_t A,
    const fmpz * alpha,
    const fmpz_mpoly_ctx_t ctx,
    const fmpz_mod_ctx_t fpctx)
{
    slong i, j;
    slong offset, shift;
    slong N = mpoly_words_per_exp_sp(A->bits, ctx->minfo);
    slong nvars = ctx->minfo->nvars;
    ulong * Aexp;
    slong * LUToffset;
    ulong * LUTmask;
    fmpz * LUTvalue;
    slong LUTlen;
    fmpz_t xpoweval;
    ulong * inputexpmask;
    TMP_INIT;

    TMP_START;

    fmpz_init(xpoweval);

    LUToffset = (slong *) TMP_ALLOC(N*FLINT_BITS*sizeof(slong));
    LUTmask   = (ulong *) TMP_ALLOC(N*FLINT_BITS*sizeof(ulong));
    LUTvalue  = (fmpz *) TMP_ALLOC(N*FLINT_BITS*sizeof(fmpz));
    for (i = 0; i < N*FLINT_BITS; i++)
    {
        fmpz_init(LUTvalue + i);
    }

    Aexp = A->exps;

    inputexpmask = (ulong *) TMP_ALLOC(N*sizeof(ulong));
    mpoly_monomial_zero(inputexpmask, N);
    for (i = 0; i < A->length; i++)
    {
        for (j = 0; j < N; j++)
        {
            inputexpmask[j] |= (Aexp + N*i)[j];
        }
    }

    LUTlen = 0;
    for (j = nvars - 1; j >= 0; j--)
    {
        mpoly_gen_offset_shift_sp(&offset, &shift, j, A->bits, ctx->minfo);

        fmpz_set(xpoweval, alpha + j); /* xpoweval = alpha[i]^(2^i) */
        for (i = 0; i < A->bits; i++)
        {
            LUToffset[LUTlen] = offset;
            LUTmask[LUTlen] = (UWORD(1) << (shift + i));
            fmpz_set(LUTvalue + LUTlen, xpoweval);
            if ((inputexpmask[offset] & LUTmask[LUTlen]) != 0)
            {
                LUTlen++;
            }
            fmpz_mod_mul(xpoweval, xpoweval, xpoweval, fpctx);
        }
    }
    FLINT_ASSERT(LUTlen < N*FLINT_BITS);

    fmpz_zero(eval);
    for (i = 0; i < A->length; i++)
    {
        fmpz_mod(xpoweval, A->coeffs + i, fpctx->p);
        for (j = 0; j < LUTlen; j++)
        {
            if (((Aexp + N*i)[LUToffset[j]] & LUTmask[j]) != 0)
            {
                fmpz_mod_mul(xpoweval, xpoweval, LUTvalue + j, fpctx);
            }
        }
        fmpz_mod_add(eval, eval, xpoweval, fpctx);
    }

    fmpz_clear(xpoweval);
    for (i = 0; i < N*FLINT_BITS; i++)
    {
        fmpz_clear(LUTvalue + i);
    }
    TMP_END;

/*
    printf("fmpz_mpoly_eval_fmpz_mod returning\n");
*/
}

/* take coeffs from [0, p) to (-p/2, p/2]
*/
void fmpz_mpolyu_symmetrize_coeffs(
    fmpz_mpolyu_t A,
    const fmpz_mpoly_ctx_t ctx,
    const fmpz_mod_ctx_t fpctx)
{
    slong i, j;
    fmpz_mpoly_struct * Ac;

    for (i = 0; i < A->length; i++)
    {
        Ac = A->coeffs + i;
        for (j = 0; j < Ac->length; j++)
        {
            fmpz_mods(Ac->coeffs + j, Ac->coeffs + j, fpctx->p);
        }
    }
}


void fmpz_mpolyuu_eval_fmpz_mod(
    fmpz_mod_mpolyun_t E,
    const fmpz_mpolyu_t A,
    const fmpz * alpha,
    const fmpz_mpoly_ctx_t ctx,
    const fmpz_mod_ctx_t fpctx)
{
    slong N = mpoly_words_per_exp_sp(A->bits, ctx->minfo);
    slong i;
    slong xexp, yexp;
    fmpz_t eval;

    FLINT_ASSERT(E->bits == A->bits);

    fmpz_init(eval);

    E->length = 0;
    for (i = 0; i < A->length; i++)
    {
        fmpz_mpoly_eval_fmpz_mod(eval, A->coeffs + i, alpha, ctx, fpctx);
        if (fmpz_is_zero(eval))
        {
            continue;
        }

        xexp = A->exps[i] >> (FLINT_BITS/2);
        yexp = A->exps[i] & (-UWORD(1) >> (FLINT_BITS - FLINT_BITS/2));

/*
printf("eval: "); fmpz_print(eval); printf("\n");
flint_printf("exp: (%wd, %wd)\n", xexp, yexp);
*/

        if (E->length > 0 && E->exps[E->length - 1] == xexp)
        {
            fmpz_mod_poly_set_coeff_fmpz(E->coeffs[E->length-1].coeffs + 0, yexp, eval);
        }
        else
        {
            fmpz_mod_mpolyun_fit_length(E, E->length + 1, ctx, fpctx);
            fmpz_mod_mpolyn_fit_length(E->coeffs + E->length, 1, ctx, fpctx);
            (E->coeffs + E->length)->length = 1;
            mpoly_monomial_zero((E->coeffs + E->length)->exps, N);
            fmpz_mod_poly_zero((E->coeffs + E->length)->coeffs + 0);
            fmpz_mod_poly_set_coeff_fmpz((E->coeffs + E->length)->coeffs + 0, yexp, eval);
            E->exps[E->length] = xexp;
            E->length++;
        }
    }
/*
flint_printf("Elength: %wd\n", E->length);
*/
    fmpz_clear(eval);
}


int fmpz_mod_mpolyn_equal(
    const fmpz_mod_mpolyn_t A,
    const fmpz_mod_mpolyn_t B,
    const fmpz_mpoly_ctx_t ctx,
    const fmpz_mod_ctx_t fpctx)
{
    slong N = mpoly_words_per_exp(A->bits, ctx->minfo);
    slong i;

    FLINT_ASSERT(A->bits == B->bits);

    if (A->length != B->length)
    {
        return 0;
    }
    for (i = 0; i < A->length; i++)
    {
        if (!mpoly_monomial_equal(A->exps + N*i, B->exps + N*i, N))
        {
            return 0;
        }
        if (!fmpz_mod_poly_equal(A->coeffs + i, B->coeffs + i))
        {
            return 0;
        }
    }
    return 1;
}

int fmpz_mod_mpolyun_equal(
    const fmpz_mod_mpolyun_t A,
    const fmpz_mod_mpolyun_t B,
    const fmpz_mpoly_ctx_t ctx,
    const fmpz_mod_ctx_t fpctx)
{
    slong i;

    FLINT_ASSERT(A->bits == B->bits);

    if (A->length != B->length)
    {
        return 0;
    }
    for (i = 0; i < A->length; i++)
    {
        if (A->exps[i] != B->exps[i])
        {
flint_printf("un exponents differ at i = %wd\n", i);
            return 0;
        }
        if (!fmpz_mod_mpolyn_equal(A->coeffs + i, B->coeffs + i, ctx, fpctx))
        {
flint_printf("un coefficints differ at i = %wd\n", i);
            return 0;
        }
    }
    return 1;
}


int fmpz_mpolyuu_gcd_bma(fmpz_mpolyu_t G, const fmpz_mpolyu_t A, const fmpz_mpolyu_t B,
                                const fmpz_mpoly_t Gamma, const fmpz_mpoly_ctx_t ctx)
{
    int changed, success;
    mp_bitcnt_t bits = A->bits, Hbits;
    fmpz_mod_mpoly_bma_interpolate_ctx_t Ictx;
    fmpz_mod_bma_mpoly_t Lambda;
    fmpz_mpolyu_t H;
    fmpz_mpoly_t Gammainc, Gammacur;
    fmpz_mpolyu_t Ainc, Acur, Binc, Bcur;
    fmpz_mod_mpolyun_t Aeval, Beval, Geval, Abareval, Bbareval, Heval;
    slong i, j;
    ulong Gdegboundxy, Gevaldegxy;
    slong * Gdegbounds, * Adegs, * Bdegs, * Gammadegs;
    flint_rand_t randstate;
    fmpz_t p, t, shift, subprod, cAksub, cBksub, sshift, last_unlucky_sshift_plus_1, image_count;
    fmpz_t Gammaeval;
    fmpz * checkalpha;
    int unlucky_count;
    TMP_INIT;

    FLINT_ASSERT(bits == A->bits);
    FLINT_ASSERT(bits == B->bits);
    FLINT_ASSERT(bits == G->bits);
    FLINT_ASSERT(bits == Gamma->bits);

    TMP_START;
    flint_randinit(randstate);
    fmpz_init(p);
    fmpz_init(t);
    fmpz_init(shift);
    fmpz_init(image_count);
    fmpz_init(subprod);
    fmpz_init(cAksub);
    fmpz_init(cBksub);
    fmpz_init(sshift);
    fmpz_init(last_unlucky_sshift_plus_1);
    fmpz_init(Gammaeval);

    fmpz_mpoly_init3(Gammainc, Gamma->length, bits, ctx);
    fmpz_mpoly_init3(Gammacur, Gamma->length, bits, ctx);
    fmpz_mpolyu_init(Ainc, bits, ctx);
    fmpz_mpolyu_init(Acur, bits, ctx);
    fmpz_mpolyu_init(Binc, bits, ctx);
    fmpz_mpolyu_init(Bcur, bits, ctx);

printf("fmpz_mpolyuu_gcd_bma called\n");
printf("Gamma: "); fmpz_mpoly_print_pretty(Gamma, NULL, ctx); printf("\n");
printf("    A: "); fmpz_mpolyuu_print_pretty(A, NULL, ctx); printf("\n");
printf("    B: "); fmpz_mpolyuu_print_pretty(B, NULL, ctx); printf("\n");

    /* find a degree bound on the two main variables */
    Gdegboundxy = FLINT_MIN(A->exps[0], B->exps[0]);

flint_printf("Gdegboundxy: %016llx\n", Gdegboundxy);


    checkalpha = (fmpz *) TMP_ALLOC(ctx->minfo->nvars*sizeof(fmpz));
    for (i = 0; i < ctx->minfo->nvars; i++)
    {
        fmpz_init(checkalpha + i);
    }

    /* find degree bounds on lesser variables */
    Gdegbounds = (slong *) TMP_ALLOC(ctx->minfo->nvars*sizeof(slong));
    Adegs = (slong *) TMP_ALLOC(ctx->minfo->nvars*sizeof(slong));
    Bdegs = (slong *) TMP_ALLOC(ctx->minfo->nvars*sizeof(slong));
    Gammadegs = (slong *) TMP_ALLOC(ctx->minfo->nvars*sizeof(slong));
    mpoly_degrees_si(Gammadegs, Gamma->exps, Gamma->length, bits, ctx->minfo);
    for (i = 0; i < ctx->minfo->nvars; i++)
    {
        Gdegbounds[i] = _degree_bound(Adegs + i, Bdegs + i, A, B, i, ctx, randstate);
flint_printf("Gdegbound[%wd]: %wd\n", i, Gdegbounds[i]);
    }

    fmpz_set_ui(p, 997);
    fmpz_mod_mpoly_bma_interpolate_ctx_init(Ictx, ctx, p);
    Ictx->bits = bits;
    fmpz_mod_bma_mpoly_init(Lambda, Ictx);

    fmpz_mod_mpolyun_init(Aeval, bits, Ictx->polyctx, Ictx->fpctx);
    fmpz_mod_mpolyun_init(Beval, bits, Ictx->polyctx, Ictx->fpctx);
    fmpz_mod_mpolyun_init(Geval, bits, Ictx->polyctx, Ictx->fpctx);
    fmpz_mod_mpolyun_init(Abareval, bits, Ictx->polyctx, Ictx->fpctx);
    fmpz_mod_mpolyun_init(Bbareval, bits, Ictx->polyctx, Ictx->fpctx);
    fmpz_mod_mpolyun_init(Heval, bits, Ictx->polyctx, Ictx->fpctx);

    for (i = 0; i < ctx->minfo->nvars; i++)
    {        
        Ictx->degbounds[i] = FLINT_MIN(Adegs[i], Bdegs[i]);
        Ictx->degbounds[i] = FLINT_MIN(Ictx->degbounds[i], Gdegbounds[i] + Gammadegs[i]) + 1;
        Ictx->subdegs[i] = Ictx->degbounds[i];
        if ((slong)(Ictx->subdegs[i]) < 0)
        {
            /* absolute falure */
            FLINT_ASSERT(0);
        }
    }

    /* find bits into which H can be packed - degrees of H are < Ictx->degbounds */
    Hbits = bits;
    for (i = 0; i < ctx->minfo->nvars; i++)
    {
        mp_bitcnt_t asdfad = 1 + FLINT_BIT_COUNT(Ictx->degbounds[i] - 1);
        Hbits = FLINT_MAX(Hbits, asdfad);
    }

    if (Hbits > FLINT_BITS)
    {
        /* absolute falure */
        FLINT_ASSERT(0);
    }

    fmpz_mpolyu_init(H, Hbits, ctx);

    goto got_ksub;

pick_ksub:

    if (ctx->minfo->nvars > 2)
    {
        j = 1;
        for (i = 2; i < ctx->minfo->nvars; i++)
        {
            if (Ictx->subdegs[i] < Ictx->subdegs[j])
            {
                j = i;
            }
        }
        Ictx->subdegs[j]++;
        if ((slong)(Ictx->subdegs[j]) < 0)
        {
            /* absolute falure */
            FLINT_ASSERT(0);
        }
    }

got_ksub:

    fmpz_one(subprod);
    for (i = 0; i < ctx->minfo->nvars; i++)
    {
        fmpz_mul_ui(subprod, subprod, Ictx->subdegs[i]);
flint_printf("ksub[%wu]: %wu\n", i, Ictx->subdegs[i]);
    }

    _ksub_content(cAksub, A->coeffs + 0, Ictx);
    _ksub_content(cBksub, B->coeffs + 0, Ictx);

    if (fmpz_is_zero(cAksub) || fmpz_is_zero(cBksub))
    {
        /* try a new substitution if we killed either leading coefficient */
        goto pick_ksub;
    }

pick_prime:

    if (fmpz_cmp(p, subprod) < 0)
        fmpz_set(p, subprod);
    fmpz_add_ui(p, p, 1);
    if (fmpz_is_probabprime(p) != 1)
        goto pick_prime;    
    if (fmpz_is_prime(p) != 1)
        goto pick_prime;

    /* make sure reduction does not kill either leading coeff after ksub */
    fmpz_gcd(t, p, cAksub);
    if (!fmpz_is_one(t))
        goto pick_prime;
    fmpz_gcd(t, p, cBksub);
    if (!fmpz_is_one(t))
        goto pick_prime;

    fmpz_one(sshift);

    unlucky_count = 0;
    fmpz_zero(last_unlucky_sshift_plus_1);

    fmpz_mod_mpoly_bma_interpolate_ctx_reset_prime(Ictx, p);

    fmpz_mod_bma_mpoly_reset_prime(Lambda, Ictx);
    fmpz_mod_bma_mpoly_zero(Lambda, Ictx);

printf("Lambda: "); fmpz_mod_bma_mpoly_print(Lambda, Ictx); printf("\n");

    _fmpz_mpoly_set_skel(Gammainc, Gamma, Ictx);
/*printf("Gammainc: "); fmpz_mpoly_print_pretty(Gammainc, NULL, Ictx->polyctx); printf("\n");*/
    _fmpz_mpoly_copy_skel(Gammacur, Gammainc, Ictx);
/*printf("Gammacur: "); fmpz_mpoly_print_pretty(Gammacur, NULL, Ictx->polyctx); printf("\n");*/
    _fmpz_mpolyuu_set_skel(Ainc, A, Ictx);
/*printf("Ainc: "); fmpz_mpolyuu_print_pretty(Ainc, NULL, Ictx->polyctx); printf("\n");*/
    _fmpz_mpolyuu_copy_skel(Acur, Ainc, Ictx);
/*printf("Acur: "); fmpz_mpolyuu_print_pretty(Acur, NULL, Ictx->polyctx); printf("\n");*/
    _fmpz_mpolyuu_set_skel(Binc, B, Ictx);
/*printf("Binc: "); fmpz_mpolyuu_print_pretty(Binc, NULL, Ictx->polyctx); printf("\n");*/
    _fmpz_mpolyuu_copy_skel(Bcur, Binc, Ictx);
/*printf("Bcur: "); fmpz_mpolyuu_print_pretty(Bcur, NULL, Ictx->polyctx); printf("\n");*/

    fmpz_zero(image_count);

next_image:

    if (0)
    {
        /* out of evaluation points */
        goto pick_prime;
    }

printf("next_image  sshift: "); fmpz_print(sshift); printf("\n");
printf("Lambda: "); fmpz_mod_bma_mpoly_print(Lambda, Ictx); printf("\n");
usleep(1000000);

    /* image count is the current power of alpha we are evaluating */
    fmpz_add_ui(image_count, image_count, 1);
    /* image_count == sshift + Lambda->pointcount */
    fmpz_add_ui(t, sshift, Lambda->pointcount);
    FLINT_ASSERT(fmpz_equal(t, image_count));

    _fmpz_mpoly_use_skel_mul(Gammaeval, Gamma, Gammacur, Gammainc, Ictx);
/*printf("A2: "); fmpz_mpolyuu_print_pretty(A, NULL, Ictx->polyctx); printf("\n");*/

    _fmpz_mpolyuu_use_skel_mul(Aeval, A, Acur, Ainc, Ictx);

/*printf("A3: "); fmpz_mpolyuu_print_pretty(A, NULL, Ictx->polyctx); printf("\n");*/

    _fmpz_mpolyuu_use_skel_mul(Beval, B, Bcur, Binc, Ictx);
/*
printf("Gammaeval: "); fmpz_print(Gammaeval); printf("\n");
*/
printf("Aeval: "); fmpz_mod_mpolyun_print_pretty(Aeval, NULL, Ictx->polyctx, Ictx->fpctx); printf("\n");
printf("Beval: "); fmpz_mod_mpolyun_print_pretty(Beval, NULL, Ictx->polyctx, Ictx->fpctx); printf("\n");

/*
printf("Gammacur: "); fmpz_mpoly_print_pretty(Gammacur, NULL, Ictx->polyctx); printf("\n");
printf("Acur: "); fmpz_mpolyuu_print_pretty(Acur, NULL, Ictx->polyctx); printf("\n");
printf("Bcur: "); fmpz_mpolyuu_print_pretty(Bcur, NULL, Ictx->polyctx); printf("\n");
*/


    if (fmpz_is_zero(Gammaeval))
    {
        fmpz_add_ui(sshift, sshift, Lambda->pointcount + 1);
        fmpz_mod_bma_mpoly_zero(Lambda, Ictx);
        goto next_image;
    }

    FLINT_ASSERT(Aeval->length > 0);
    FLINT_ASSERT(Beval->length > 0);

    success = fmpz_mod_mpolyun_gcd_bivar(Geval, Abareval, Bbareval, Aeval, Beval, Ictx->polyctx, Ictx->fpctx);
    FLINT_ASSERT(success);

    fmpz_mod_mpolyun_scalar_mul_fmpz_mod(Geval, Gammaeval, Ictx->polyctx, Ictx->fpctx);
printf("Geval: "); fmpz_mod_mpolyun_print_pretty(Geval, NULL, Ictx->polyctx, Ictx->fpctx); printf("\n");

    FLINT_ASSERT(Geval->length > 0);
    Gevaldegxy = fmpz_mod_mpolyun_bidegree(Geval, Ictx->polyctx, Ictx->fpctx);

    FLINT_ASSERT(fmpz_equal(Gammaeval, (Geval->coeffs[0].coeffs + 0)->coeffs + fmpz_mod_poly_degree(Geval->coeffs[0].coeffs + 0)));

    if (Gdegboundxy < Gevaldegxy)
    {
printf("unlucky\n");
        if (fmpz_equal(sshift, last_unlucky_sshift_plus_1))
        {
            goto pick_ksub;
        }
        if (++unlucky_count > 2)
        {
            goto pick_prime;
        }
        fmpz_add_ui(last_unlucky_sshift_plus_1, sshift, 1);
        fmpz_add_ui(sshift, sshift, Lambda->pointcount + 1);
        fmpz_mod_bma_mpoly_zero(Lambda, Ictx);
        goto next_image;        
    }
    else if (Gdegboundxy > Gevaldegxy)
    {
printf("all previous were unlucky\n");

        Gdegboundxy = Gevaldegxy;
        fmpz_add_ui(sshift, sshift, Lambda->pointcount);
        fmpz_mod_bma_mpoly_zero(Lambda, Ictx);
        fmpz_mod_bma_mpoly_add_point(Lambda, Geval, Ictx);
        goto next_image;
    }

    fmpz_mod_bma_mpoly_add_point(Lambda, Geval, Ictx);

printf("Lambda: "); fmpz_mod_bma_mpoly_print(Lambda, Ictx); printf("\n");


    if ((Lambda->pointcount & 1) != 0)
    {
        goto next_image;
    }

    changed = fmpz_mod_bma_mpoly_reduce(Lambda, Ictx);

printf("changed: %d\n", changed);
printf("Lambda: "); fmpz_mod_bma_mpoly_print(Lambda, Ictx); printf("\n");

    if (changed)
    {
        goto next_image;
    }

    success = fmpz_mod_bma_mpoly_get_mpolyu(H, sshift, Lambda, Ictx);
    if (!success)
    {
        goto next_image;
    }

printf("H: "); fmpz_mpolyuu_print_pretty(H, NULL, ctx); printf("\n");

    /* evaluate Gamma, A, and B at random point  */

    for (i = 0; i < ctx->minfo->nvars; i++)
    {
        fmpz_randm(checkalpha + i, randstate, p);
flint_printf("checkalpha[%wd]: ", i); fmpz_print(checkalpha + i); printf("\n");
    }
    fmpz_mpolyuu_eval_fmpz_mod(Aeval, A, checkalpha, Ictx->polyctx, Ictx->fpctx);

printf("Aeval: "); fmpz_mod_mpolyun_print_pretty(Aeval, NULL, Ictx->polyctx, Ictx->fpctx); printf("\n");

    fmpz_mpolyuu_eval_fmpz_mod(Beval, B, checkalpha, Ictx->polyctx, Ictx->fpctx);


printf("Beval: "); fmpz_mod_mpolyun_print_pretty(Beval, NULL, Ictx->polyctx, Ictx->fpctx); printf("\n");


    FLINT_ASSERT(Aeval->length > 0 && fmpz_mod_mpolyun_bidegree(Aeval, Ictx->polyctx, Ictx->fpctx) == A->exps[0]);
    FLINT_ASSERT(Beval->length > 0 && fmpz_mod_mpolyun_bidegree(Beval, Ictx->polyctx, Ictx->fpctx) == B->exps[0]);
    fmpz_mpoly_eval_fmpz_mod(Gammaeval, Gamma, checkalpha, Ictx->polyctx, Ictx->fpctx);
    FLINT_ASSERT(!fmpz_is_zero(Gammaeval));

    success = fmpz_mod_mpolyun_gcd_bivar(Geval, Abareval, Bbareval, Aeval, Beval, Ictx->polyctx, Ictx->fpctx);
    FLINT_ASSERT(success);
    fmpz_mod_mpolyun_scalar_mul_fmpz_mod(Geval, Gammaeval, Ictx->polyctx, Ictx->fpctx);
printf("Geval: "); fmpz_mod_mpolyun_print_pretty(Geval, NULL, Ictx->polyctx, Ictx->fpctx); printf("\n");
    FLINT_ASSERT(Geval->length > 0 && fmpz_mod_mpolyun_bidegree(Geval, Ictx->polyctx, Ictx->fpctx) == H->exps[0]);

    fmpz_mpolyuu_eval_fmpz_mod(Heval, H, checkalpha, Ictx->polyctx, Ictx->fpctx);


printf("Heval: "); fmpz_mod_mpolyun_print_pretty(Heval, NULL, Ictx->polyctx, Ictx->fpctx); printf("\n");



    if (!fmpz_mod_mpolyun_equal(Heval, Geval, Ictx->polyctx, Ictx->fpctx))
    {
printf("not equal\n");
FLINT_ASSERT(0);
        goto next_image;
    }

printf("H: "); fmpz_mpolyuu_print_pretty(H, NULL, ctx); printf("\n");
    fmpz_mpolyu_symmetrize_coeffs(H, ctx, Ictx->fpctx);
printf("H: "); fmpz_mpolyuu_print_pretty(H, NULL, ctx); printf("\n");

    FLINT_ASSERT(0);

    fmpz_mod_mpolyun_clear(Heval, Ictx->polyctx, Ictx->fpctx);

    fmpz_mpolyu_clear(H, ctx);

    fmpz_mpoly_clear(Gammainc, ctx);
    fmpz_mpoly_clear(Gammacur, ctx);

    for (i = 0; i < ctx->minfo->nvars; i++)
    {
        fmpz_clear(checkalpha + i);
    }

    fmpz_clear(Gammaeval);
    fmpz_clear(p);
    fmpz_clear(image_count);
    fmpz_clear(t);
    fmpz_clear(shift);
    fmpz_clear(subprod);
    fmpz_clear(cAksub);
    fmpz_clear(cBksub);
    fmpz_clear(sshift);
    fmpz_clear(last_unlucky_sshift_plus_1);
    flint_randclear(randstate);

    TMP_END;
    return 0;
}




int fmpz_mpoly_gcd_bma(fmpz_mpoly_t G,
        const fmpz_mpoly_t A, const fmpz_mpoly_t B, const fmpz_mpoly_ctx_t ctx)
{
    slong i;
    mp_bitcnt_t Gbits, ABbits;
    int success = 0;
    fmpz_mpoly_ctx_t uctx;
    fmpz_mpolyu_t Auu, Buu, Guu, Abar, Bbar, Gbar;
    fmpz_mpoly_t Acontent, Bcontent, Gamma;
    slong * Adegs, * Bdegs, * perm;
    ulong * shift, * stride;
    TMP_INIT;

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

    if (ctx->minfo->nvars < 3)
    {
        return fmpz_mpoly_gcd_zippel(G, A, B, ctx);
    }

    FLINT_ASSERT(A->bits <= FLINT_BITS);
    FLINT_ASSERT(B->bits <= FLINT_BITS);

    FLINT_ASSERT(ctx->minfo->nvars >= WORD(3));
    FLINT_ASSERT(!fmpz_mpoly_is_zero(A, ctx));
    FLINT_ASSERT(!fmpz_mpoly_is_zero(B, ctx));

    TMP_START;

    Adegs = (slong *) TMP_ALLOC(ctx->minfo->nvars*sizeof(slong));
    Bdegs = (slong *) TMP_ALLOC(ctx->minfo->nvars*sizeof(slong));
    mpoly_degrees_si(Adegs, A->exps, A->length, A->bits, ctx->minfo);
    mpoly_degrees_si(Bdegs, B->exps, B->length, B->bits, ctx->minfo);

    for (i = 0; i < ctx->minfo->nvars; i++)
    {
        if (Adegs[i] >= 100 || Bdegs[i] >= 100)
        {
            success = 0;
            goto cleanup1;
        }
    }

    perm = (slong *) TMP_ALLOC((ctx->minfo->nvars)*sizeof(slong));
    shift = (ulong *) TMP_ALLOC((ctx->minfo->nvars)*sizeof(ulong));
    stride = (ulong *) TMP_ALLOC((ctx->minfo->nvars)*sizeof(ulong));
    for (i = 0; i < ctx->minfo->nvars; i++)
    {
        perm[i] = i;
        shift[i] = 0;
        stride[i] = 1;
    }
    /* ABbits is bits for itermediates in Z[X,Y][z_0,...,z_(n-3)] (in lex) */
    ABbits = 8;
    /* Gbits is bits for final answer in Z[x_0,...,x_(n-1)] */
    Gbits = FLINT_MIN(A->bits, B->bits);

    fmpz_mpoly_ctx_init(uctx, ctx->minfo->nvars - 2, ORD_LEX);
    fmpz_mpolyu_init(Auu, ABbits, uctx);
    fmpz_mpolyu_init(Buu, ABbits, uctx);
    fmpz_mpolyu_init(Guu, ABbits, uctx);

    fmpz_mpoly_to_mpolyuu_perm_deflate(Auu, perm, shift, stride, uctx, A, ctx);
    fmpz_mpoly_to_mpolyuu_perm_deflate(Buu, perm, shift, stride, uctx, B, ctx);

    fmpz_mpoly_init(Acontent, uctx);
    fmpz_mpoly_init(Bcontent, uctx);
    fmpz_mpoly_init(Gamma, uctx);
    fmpz_mpolyu_init(Abar, ABbits, uctx);
    fmpz_mpolyu_init(Bbar, ABbits, uctx);
    fmpz_mpolyu_init(Gbar, ABbits, uctx);

    /* compute content of A */
    fmpz_mpoly_set(Acontent, Auu->coeffs + 0, uctx);
    for (i = 1; i < Auu->length; i++)
    {
        success = _fmpz_mpoly_gcd(Acontent, ABbits, Acontent, Auu->coeffs + i, uctx);
        if (!success)
            goto cleanup;
        FLINT_ASSERT(Acontent->bits == ABbits);
    }

    /* compute content of B */
    fmpz_mpoly_set(Bcontent, Buu->coeffs + 0, uctx);
    for (i = 1; i < Buu->length; i++)
    {
        success = _fmpz_mpoly_gcd(Bcontent, ABbits, Bcontent, Buu->coeffs + i, uctx);
        if (!success)
            goto cleanup;
        FLINT_ASSERT(Bcontent->bits == ABbits);
    }

    /* remove content from A and B */
    fmpz_mpolyu_divexact_mpoly(Abar, Auu, Acontent, uctx);
    fmpz_mpolyu_divexact_mpoly(Bbar, Buu, Bcontent, uctx);

    /* compute GCD of leading coefficients */
    _fmpz_mpoly_gcd(Gamma, ABbits, Abar->coeffs + 0, Bbar->coeffs + 0, uctx);
    success = fmpz_mpolyuu_gcd_bma(Gbar, Abar, Bbar, Gamma, uctx);
    if (!success)
        goto cleanup;
    FLINT_ASSERT(Gamma->bits == ABbits);

    /* put back content */
    success = _fmpz_mpoly_gcd(Acontent, ABbits, Acontent, Bcontent, uctx);
    if (!success)
        goto cleanup;
    FLINT_ASSERT(Acontent->bits == ABbits);

    fmpz_mpolyu_mul_mpoly(Guu, Gbar, Acontent, uctx);
    fmpz_mpoly_from_mpolyuu_perm_inflate(G, Gbits, ctx, Guu, perm, shift, stride, uctx);
    if (fmpz_sgn(G->coeffs + 0) < 0)
        fmpz_mpoly_neg(G, G, ctx);
    FLINT_ASSERT(G->bits == Gbits);

    success = 1;

cleanup:

    fmpz_mpolyu_clear(Abar, uctx);
    fmpz_mpolyu_clear(Bbar, uctx);
    fmpz_mpolyu_clear(Gbar, uctx);
    fmpz_mpoly_clear(Acontent, uctx);
    fmpz_mpoly_clear(Bcontent, uctx);
    fmpz_mpoly_clear(Gamma, uctx);

    fmpz_mpolyu_clear(Auu, uctx);
    fmpz_mpolyu_clear(Buu, uctx);
    fmpz_mpolyu_clear(Guu, uctx);
    fmpz_mpoly_ctx_clear(uctx);


cleanup1:

    TMP_END;

    return success;
}


int
main(void)
{
    slong i;
    FLINT_TEST_INIT(state);

    flint_printf("gcd_brown....\n");
    fflush(stdout);

    {
        fmpz_mpoly_ctx_t ctx;
        fmpz_mpoly_t G, A, B;
        const char* vars[] = {"x", "y", "z", "t"};

        fmpz_mpoly_ctx_init(ctx, 4, ORD_LEX);
        fmpz_mpoly_init(G, ctx);
        fmpz_mpoly_init(A, ctx);
        fmpz_mpoly_init(B, ctx);
        fmpz_mpoly_set_str_pretty(A, "z*(z*x-2*t)", vars, ctx);
        fmpz_mpoly_set_str_pretty(B, "(t*y+z)", vars, ctx);
        fmpz_mpoly_set_str_pretty(G, "t*((z-3*t)*x+y+5*z-t)", vars, ctx);
printf("A: "); fmpz_mpoly_print_pretty(A, vars, ctx); printf("\n");
printf("B: "); fmpz_mpoly_print_pretty(B, vars, ctx); printf("\n");
printf("G: "); fmpz_mpoly_print_pretty(G, vars, ctx); printf("\n");
        fmpz_mpoly_mul(A, A, G, ctx);
        fmpz_mpoly_mul(B, B, G, ctx);
printf("A: "); fmpz_mpoly_print_pretty(A, vars, ctx); printf("\n");
printf("B: "); fmpz_mpoly_print_pretty(B, vars, ctx); printf("\n");
        fmpz_mpoly_gcd_bma(G, A, B, ctx);
        FLINT_ASSERT(0);
    }

    if (0)
    {
        fmpz_t pp, x, y;

        fmpz_init(pp);
        fmpz_init(x);
        fmpz_init(y);

        {
            ulong i;
            fmpz_mod_dlog_env_t L;

            fmpz_set_ui(pp, 4085);
            fmpz_mul_2exp(pp, pp, 117);
            fmpz_add_ui(pp, pp, 1);

            fmpz_mod_dlog_env_init(L, pp);
            for (i = 0; i < 1000; i++)
            {
                fmpz_mod_pow_ui(x, L->alpha, i, L->fpctx);
                fmpz_mod_dlog_env_run(L, y, x);
                FLINT_ASSERT(i == fmpz_get_ui(y));
            }
            fmpz_mod_dlog_env_clear(L);

        }

        fmpz_clear(pp);
        fmpz_clear(x);
        fmpz_clear(y);
    }

    if(0){
        int success;
        fmpz_t eval, p, shift;
        fmpz_mod_mpoly_bma_interpolate_ctx_t Ictx;
        fmpz_mod_mpoly_bma_interpolate_t I;
        fmpz_mpoly_ctx_t ctx;
        fmpz_mpoly_t A, S, M, R;
        const char * vars[] = {"x", "y", "z"};

        fmpz_init(eval);
        fmpz_init(p);
        fmpz_init_set_ui(shift, 1);

        fmpz_set_ui(p, 4085);
        fmpz_mul_2exp(p, p, 117);
        fmpz_add_ui(p, p, 1);

/*
        fmpz_set_ui(p, 1009);
*/
        fmpz_mpoly_ctx_init(ctx, 3, ORD_LEX);
        fmpz_mpoly_init(A, ctx);
        fmpz_mpoly_set_str_pretty(A, "10*z^3 + y + 3*x*y^4 + 2*x^2 + x^3*z^5", vars, ctx);

printf("A: "); fmpz_mpoly_print_pretty(A, vars, ctx); printf("\n");
FLINT_ASSERT(0);
/*
        fmpz_mod_mpoly_bma_interpolate_ctx_init(Ictx, A->bits, ctx, p);
*/
flint_printf("alpha: "); fmpz_print(Ictx->dlogenv->alpha); printf("\n");

        fmpz_mod_mpoly_bma_interpolate_init(I, p);

        fmpz_mod_mpoly_bma_interpolate_eval_init(Ictx, A);
        Ictx->degbounds[0] = 4;
        Ictx->subdegs[0] = 4;
        Ictx->degbounds[1] = 5;
        Ictx->subdegs[1] = 5;
        Ictx->degbounds[2] = 6;
        Ictx->subdegs[2] = 6;

        fmpz_mpoly_init3(S, A->length, A->bits, ctx);
        fmpz_mpoly_init3(M, A->length, A->bits, ctx);
        fmpz_mpoly_init3(R, A->length, A->bits, ctx);

        FLINT_ASSERT(0);
/*
        fmpz_mod_mpoly_bma_interpolate_eval_setskel(S, A, 1, Ictx);
*/
        fmpz_mod_mpoly_bma_interpolate_eval_copyskel(M, S, Ictx);

for (i = 0; i < 10; i++)
{
        fmpz_mod_mpoly_bma_interpolate_eval_useskelmul(eval, A, M, S, Ictx);
        fmpz_mod_mpoly_bma_interpolate_add_point(I, eval);
        fmpz_mod_mpoly_bma_interpolate_eval_useskelmul(eval, A, M, S, Ictx);
        fmpz_mod_mpoly_bma_interpolate_add_point(I, eval);
        fmpz_mod_mpoly_bma_interpolate_reduce(I);
flint_printf("master: "); fmpz_mod_poly_print_pretty(I->bma->V1, "#"); printf("\n");
}

        success = fmpz_mod_mpoly_bma_interpolate_get_mpoly(R, shift, I, Ictx);
        fmpz_mpoly_sort_terms(R, ctx);
printf("success: %d\n", success);
printf("R: "); fmpz_mpoly_print_pretty(R, vars, ctx); printf("\n");

        fmpz_mod_mpoly_bma_interpolate_clear(I);

        fmpz_mod_mpoly_bma_interpolate_ctx_clear(Ictx);

        fmpz_mpoly_clear(S, ctx);
        fmpz_mpoly_clear(M, ctx);
        fmpz_mpoly_clear(R, ctx);
        fmpz_mpoly_clear(A, ctx);
        fmpz_mpoly_ctx_clear(ctx);

        fmpz_clear(eval);
        fmpz_clear(p);
        fmpz_clear(shift);
    }




    if(0){
        int success;
        mp_limb_t eval;
        nmod_mpoly_bma_interpolate_ctx_t Ictx;
        nmod_mpoly_bma_interpolate_t I;
        nmod_mpoly_ctx_t ctx;
        nmod_mpoly_t A, S, M, R;
        const char * vars[] = {"x", "y", "z"};

        nmod_mpoly_ctx_init(ctx, 3, ORD_LEX, 1009);
        nmod_mpoly_init(A, ctx);
        nmod_mpoly_set_str_pretty(A, "10*z^3 + y + 3*x*y^4 + 2*x^2 + x^3*z^5", vars, ctx);

printf("A: "); nmod_mpoly_print_pretty(A, vars, ctx); printf("\n");

        nmod_mpoly_bma_interpolate_ctx_init(Ictx, A->bits, ctx);

flint_printf("alpha: %wu\n", Ictx->dlogenv->alpha);

        nmod_mpoly_bma_interpolate_init(I, ctx->ffinfo->mod.n);

        nmod_mpoly_bma_interpolate_eval_init(Ictx, A);
        Ictx->degbounds[0] = 4;
        Ictx->subdegs[0] = 4;
        Ictx->degbounds[1] = 5;
        Ictx->subdegs[1] = 5;
        Ictx->degbounds[2] = 6;
        Ictx->subdegs[2] = 6;

        nmod_mpoly_init3(S, A->length, A->bits, ctx);
        nmod_mpoly_init3(M, A->length, A->bits, ctx);
        nmod_mpoly_init3(R, A->length, A->bits, ctx);

        nmod_mpoly_bma_interpolate_eval_setskel(S, A, 1, Ictx);
        nmod_mpoly_bma_interpolate_eval_copyskel(M, S, Ictx);

for (i = 0; i < 10; i++)
{
        eval = nmod_mpoly_bma_interpolate_eval_useskelmul(A, M, S, Ictx);
        nmod_mpoly_bma_interpolate_add_point(I, eval);
        eval = nmod_mpoly_bma_interpolate_eval_useskelmul(A, M, S, Ictx);
        nmod_mpoly_bma_interpolate_add_point(I, eval);
        nmod_mpoly_bma_interpolate_reduce(I);
flint_printf("master: "); nmod_poly_print_pretty(I->bma->V1, "#"); printf("\n");
}

        success = nmod_mpoly_bma_interpolate_get_mpoly(R, 1, I, Ictx);
        nmod_mpoly_sort_terms(R, ctx);
printf("success: %d\n", success);
printf("R: "); nmod_mpoly_print_pretty(R, vars, ctx); printf("\n");

        nmod_mpoly_bma_interpolate_clear(I);

        nmod_mpoly_bma_interpolate_ctx_clear(Ictx);


        nmod_mpoly_clear(A, ctx);
        nmod_mpoly_ctx_clear(ctx);
    }


    if (0) {
        mp_limb_t p = 0;
        while (p<5)
        {
            ulong i;
            nmod_dlog_env_t L;
            p = n_nextprime(p, 1);
            nmod_dlog_env_init(L, p);
            for (i = 0; i < p - 1; i++)
            {
                FLINT_ASSERT(i == nmod_dlog_env_run(L, nmod_pow_ui(L->alpha, i, L->mod)));
            }
            nmod_dlog_env_clear(L);
        }

        p = (UWORD(29) << 57);

        {
            ulong i;
            nmod_dlog_env_t L;

flint_printf("     : %wu\n", p);
            p = n_nextprime(p, 1);
flint_printf("prime: %wu\n", p);

            nmod_dlog_env_init(L, p);
            for (i = 0; i < 10000; i++)
            {
                FLINT_ASSERT(i == nmod_dlog_env_run(L, nmod_pow_ui(L->alpha, i, L->mod)));
/*                dlog_env_run(L, i + 1);*/
            }
            nmod_dlog_env_clear(L);
        }   
    }

    printf("PASS\n");
    FLINT_TEST_CLEANUP(state);

    return 0;
}
