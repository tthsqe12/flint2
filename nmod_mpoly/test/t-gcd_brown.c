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

flint_printf("L[%wd].exp  : %wu\n", i, Li->exp);
flint_printf("L[%wd].prime: %wu\n", i, Li->prime);

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



void fmpz_mod_mpoly_bma_interpolate_ctx_init(fmpz_mod_mpoly_bma_interpolate_ctx_t Ictx, mp_bitcnt_t bits_,
                                                const fmpz_mpoly_ctx_t pctx, const fmpz_t p)
{
    slong N;

    fmpz_mod_ctx_init(Ictx->fpctx, p);
    Ictx->polyctx = pctx;
    Ictx->degbounds = (slong *) flint_malloc(pctx->minfo->nvars*sizeof(slong));
    Ictx->subdegs = (ulong *) flint_malloc(pctx->minfo->nvars*sizeof(ulong));
    fmpz_mod_dlog_env_init(Ictx->dlogenv, p);

    Ictx->bits = bits_;
    N = mpoly_words_per_exp_sp(bits_, pctx->minfo);
    Ictx->inputexpmask = (ulong *) flint_malloc(N*sizeof(ulong));
    mpoly_monomial_zero(Ictx->inputexpmask, N);
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
void fmpz_mod_mpoly_bma_interpolate_eval_setskel(fmpz_mpoly_t M, const fmpz_mpoly_t A, ulong w, const fmpz_mod_mpoly_bma_interpolate_ctx_t Ictx)
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

    LUTlen = 0;
    fmpz_mod_pow_ui(xeval, Ictx->dlogenv->alpha, w, Ictx->fpctx);
    for (j = nvars - 1; j >= 0; j--)
    {
        mpoly_gen_offset_shift_sp(&offset, &shift, j, A->bits, Ictx->polyctx->minfo);

        fmpz_set(xpoweval, xeval); /* xpoweval = xeval^(2^i) */
        for (i = 0; i < A->bits; i++)
        {
            LUToffset[LUTlen] = offset;
            LUTmask[LUTlen] = (UWORD(1) << (shift + i));
            fmpz_set(LUTvalue + LUTlen, xpoweval);
            if ((Ictx->inputexpmask[offset] & LUTmask[LUTlen]) != 0)
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
void fmpz_mod_mpoly_bma_interpolate_eval_copyskel(fmpz_mpoly_t M, const fmpz_mpoly_t S, const fmpz_mod_mpoly_bma_interpolate_ctx_t Ictx)
{
    slong i, N;
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

/* return A.M and multiply M by S */
void fmpz_mod_mpoly_bma_interpolate_eval_useskelmul(fmpz_t eval, const fmpz_mpoly_t A, fmpz_mpoly_t M, const fmpz_mpoly_t S, const fmpz_mod_mpoly_bma_interpolate_ctx_t Ictx)
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
        fmpz_mod_mul(M->coeffs + i, M->coeffs + i, S->coeffs + i, Ictx->fpctx);
    }
    fmpz_clear(t);
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







int
main(void)
{
    slong i, j;
    FLINT_TEST_INIT(state);

    flint_printf("gcd_brown....\n");
    fflush(stdout);

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

    if(1){
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

        fmpz_mod_mpoly_bma_interpolate_ctx_init(Ictx, A->bits, ctx, p);

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

        fmpz_mod_mpoly_bma_interpolate_eval_setskel(S, A, 1, Ictx);
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

    for (i = 0; i < 00 * flint_test_multiplier(); i++)
    {
        nmod_mpoly_ctx_t ctx;
        nmod_mpoly_t a, b, g, ca, cb, cg, t;
        slong len, len1, len2;
        slong degbound;
        mp_limb_t modulus;
        int res;

        modulus = n_randint(state, (i % 10 == 0) ? 4: FLINT_BITS - 1) + 1;
        modulus = n_randbits(state, modulus);
        modulus = n_nextprime(modulus, 1);

        nmod_mpoly_ctx_init_rand(ctx, state, modulus < 3000 ? 4 : 5, modulus);

        nmod_mpoly_init(g, ctx);
        nmod_mpoly_init(a, ctx);
        nmod_mpoly_init(b, ctx);
        nmod_mpoly_init(ca, ctx);
        nmod_mpoly_init(cb, ctx);
        nmod_mpoly_init(cg, ctx);
        nmod_mpoly_init(t, ctx);

        len = n_randint(state, 100) + 1;
        len1 = n_randint(state, 200);
        len2 = n_randint(state, 200);

        degbound = 100/ctx->minfo->nvars/ctx->minfo->nvars;

        for (j = 0; j < 4; j++)
        {
            do {
                nmod_mpoly_randtest_bound(t, state, len, degbound, ctx);
            } while (t->length == 0);
            nmod_mpoly_randtest_bound(a, state, len1, degbound, ctx);
            nmod_mpoly_randtest_bound(b, state, len2, degbound, ctx);
            nmod_mpoly_mul_johnson(a, a, t, ctx);
            nmod_mpoly_mul_johnson(b, b, t, ctx);

            nmod_mpoly_randtest_bits(g, state, len, FLINT_BITS, ctx);

            res = nmod_mpoly_gcd_brown(g, a, b, ctx);
            if (!res) {
                continue;
            }
            nmod_mpoly_assert_canonical(g, ctx);

            if (nmod_mpoly_is_zero(g, ctx))
            {
                if (!nmod_mpoly_is_zero(a, ctx) || !nmod_mpoly_is_zero(b, ctx))
                {
                    printf("FAIL\n");
                    flint_printf("Check zero gcd only results from zero inputs\ni = %wd, j = %wd\n", i ,j);
                    flint_abort();
                }
                continue;
            }

            if (g->coeffs[0] != UWORD(1))
            {
                printf("FAIL\n");
                flint_printf("Check gcd is monic\ni = %wd, j = %wd\n", i ,j);
                flint_abort();
            }

            res = 1;
            res = res && nmod_mpoly_divides_monagan_pearce(ca, a, g, ctx);
            res = res && nmod_mpoly_divides_monagan_pearce(cb, b, g, ctx);
            if (!res)
            {
                printf("FAIL\n");
                flint_printf("Check divisibility\ni = %wd, j = %wd\n", i ,j);
                flint_abort();
            }

            res = nmod_mpoly_gcd_brown(cg, ca, cb, ctx);

            if (!res)
                continue;

            if (!nmod_mpoly_equal_ui(cg, UWORD(1), ctx))
            {
                printf("FAIL\n");
                flint_printf("Check cofactors are relatively prime\ni = %wd, j = %wd\n", i ,j);                
                flint_abort();
            }
        }

        nmod_mpoly_clear(g, ctx);
        nmod_mpoly_clear(a, ctx);
        nmod_mpoly_clear(b, ctx);
        nmod_mpoly_clear(ca, ctx);
        nmod_mpoly_clear(cb, ctx);
        nmod_mpoly_clear(cg, ctx);
        nmod_mpoly_clear(t, ctx);
        nmod_mpoly_ctx_clear(ctx);
    }


    printf("PASS\n");
    FLINT_TEST_CLEANUP(state);

    return 0;
}
