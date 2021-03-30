/*
    Copyright (C) 2010, 2018 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mat.h"


static void _distribute_rows2(
    slong start, slong stop,
    slong * Astart, slong * Astop, slong Alen,
    slong * Bstart, slong * Bstop, slong Blen)
{
    FLINT_ASSERT(0 <= start);
    FLINT_ASSERT(start <= stop);
    FLINT_ASSERT(stop <= Alen + Blen);

    if (start >= Alen)
    {
        *Astart = 0;
        *Astop  = 0;
        *Bstart = start - Alen;
        *Bstop  = stop - Alen;
    }
    else if (stop <= Alen)
    {
        *Astart = start;
        *Astop  = stop;
        *Bstart = 0;
        *Bstop  = 0;
    }
    else
    {
        *Astart = start;
        *Astop  = Alen;
        *Bstart = 0;
        *Bstop  = stop - Alen;
    }
}

/*
    if f is continuous with
        f(0) = 0
        f'(x) = alpha for 0 < x < a
        f'(x) = beta for a < x < a + b
    return solution for x to f(x) = (a*alpha + b*beta)*yn/yd
    for 0 <= yn/yd <= 1
*/
static ulong _peicewise_linear_inverse2(
    ulong a, ulong alpha,
    ulong b, ulong beta,
    ulong yn, ulong yd)
{
    /* very low priority TODO: this can overflow only in very extreme cases */
    ulong y = yn*(a*alpha + b*beta)/yd;

    if (y <= a*alpha)
        return y/alpha;

    return a + (y - a*alpha)/beta;
}


mp_limb_t fmpz_get_nmod(
    const fmpz_t aa,
    nmod_t mod)
{
    fmpz A = *aa;
    mp_limb_t r, SA, UA;
    mpz_srcptr a;
    mp_limb_t * ad;
    slong an;

    if (!COEFF_IS_MPZ(A))
    {
        SA = FLINT_SIGN_EXT(A);
        UA = FLINT_ABS(A);
        NMOD_RED(r, UA, mod);
        return r + (SA&(mod.n - 2*r));
    }

    a = COEFF_TO_PTR(A);
    ad = a->_mp_d;
    an = a->_mp_size;

    if (an < 0)
    {
        SA = -UWORD(1);
        an = -an;
    }
    else
    {
        SA = 0;
    }

    if (an < 5)
    {
        r = 0;
        while (an > 0)
        {
            NMOD_RED2(r, r, ad[an - 1], mod);
            an--;
        }
    }
    else
    {
        r = mpn_mod_1(ad, an, mod.n);
    }

    return r + (SA&(mod.n - 2*r));
}


typedef struct {
    slong m;
    slong k;
    slong n;
    slong Astartrow;
    slong Astoprow;
    slong Bstartrow;
    slong Bstoprow;
    slong Cstartrow;
    slong Cstoprow;
    fmpz ** Arows;
    fmpz ** Brows;
    fmpz ** Crows;
    nmod_mat_t * mod_A;
    nmod_mat_t * mod_B;
    nmod_mat_t * mod_C;
    const fmpz_comb_struct * comb;
    slong num_primes;
    mp_ptr primes;
} _worker_arg;


static void _mod_worker(void * varg)
{
    _worker_arg * arg = (_worker_arg *) varg;
    slong i, j, l;
    slong k = arg->k;
    slong n = arg->n;
    slong Astartrow = arg->Astartrow;
    slong Astoprow = arg->Astoprow;
    slong Bstartrow = arg->Bstartrow;
    slong Bstoprow = arg->Bstoprow;
    fmpz ** Arows = arg->Arows;
    fmpz ** Brows = arg->Brows;
    nmod_mat_t * mod_A = arg->mod_A;
    nmod_mat_t * mod_B = arg->mod_B;
    slong num_primes = arg->num_primes;
    const fmpz_comb_struct * comb = arg->comb;

    if (comb != NULL)
    {
        mp_limb_t * residues;
        fmpz_comb_temp_t comb_temp;

        residues = FLINT_ARRAY_ALLOC(num_primes, mp_limb_t);
        fmpz_comb_temp_init(comb_temp, comb);

        for (i = Astartrow; i < Astoprow; i++)
        for (j = 0; j < k; j++)
        {
            fmpz_multi_mod_ui(residues, &Arows[i][j], comb, comb_temp);
            for (l = 0; l < num_primes; l++)
                mod_A[l]->rows[i][j] = residues[l];
        }

        for (i = Bstartrow; i < Bstoprow; i++)
        for (j = 0; j < n; j++)
        {
            fmpz_multi_mod_ui(residues, &Brows[i][j], comb, comb_temp);
            for (l = 0; l < num_primes; l++)
                mod_B[l]->rows[i][j] = residues[l];
        }

        flint_free(residues);
        fmpz_comb_temp_clear(comb_temp);
    }
    else
    {
        for (i = Astartrow; i < Astoprow; i++)
        for (j = 0; j < k; j++)
        {
            for (l = 0; l < num_primes; l++)
                nmod_mat_entry(mod_A[l], i, j) = fmpz_get_nmod(&Arows[i][j],
                                                                mod_A[l]->mod);
        }

        for (i = Bstartrow; i < Bstoprow; i++)
        for (j = 0; j < n; j++)
        {
            for (l = 0; l < num_primes; l++)
                nmod_mat_entry(mod_B[l], i, j) = fmpz_get_nmod(&Brows[i][j],
                                                                mod_A[l]->mod);
        }
    }
}

static void _crt_worker(void * varg)
{
    _worker_arg * arg = (_worker_arg *) varg;
    slong i, j, l;
    slong n = arg->n;
    slong Cstartrow = arg->Cstartrow;
    slong Cstoprow = arg->Cstoprow;
    fmpz ** Crows = arg->Crows;
    nmod_mat_t * mod_C = arg->mod_C;
    mp_limb_t * primes = arg->primes;
    slong num_primes = arg->num_primes;
    const fmpz_comb_struct * comb = arg->comb;

    if (comb != NULL)
    {
        mp_limb_t * residues;
        fmpz_comb_temp_t comb_temp;

        residues = FLINT_ARRAY_ALLOC(num_primes, mp_limb_t);
        fmpz_comb_temp_init(comb_temp, comb);

        for (i = Cstartrow; i < Cstoprow; i++)
        for (j = 0; j < n; j++)
        {
            for (l = 0; l < num_primes; l++)
                residues[l] = mod_C[l]->rows[i][j];

            fmpz_multi_CRT_ui(&Crows[i][j], residues, comb, comb_temp, 1);
        }

        flint_free(residues);
        fmpz_comb_temp_clear(comb_temp);
    }
    else if (num_primes == 1)
    {
        mp_limb_t r, t, p = primes[0];

        for (i = Cstartrow; i < Cstoprow; i++)
        for (j = 0; j < n; j++)
        {
            r = nmod_mat_entry(mod_C[0], i, j);
            t = p - r;
            if (t < r)
                fmpz_neg_ui(&Crows[i][j], t);
            else
                fmpz_set_ui(&Crows[i][j], r);
        }
    }
    else if (num_primes == 2)
    {
        mp_limb_t r0, r1, i0, i1, m0, m1, M[2], t[2], u[2];
        m0 = primes[0];
        m1 = primes[1];
        i0 = n_invmod(m1 % m0, m0);
        i1 = n_invmod(m0, m1);
        umul_ppmm(M[1], M[0], m0, m1);

        for (i = Cstartrow; i < Cstoprow; i++)
        for (j = 0; j < n; j++)
        {
            r0 = nmod_mul(i0, nmod_mat_entry(mod_C[0], i, j), mod_C[0]->mod);
            r1 = nmod_mul(i1, nmod_mat_entry(mod_C[1], i, j), mod_C[1]->mod);

            /*
                0 <= r0 <= m0 - 1
                0 <= r1 <= m1 - 1
                if m0*m1 < 2^(2*FLINT_BITS-1), then the sum
                    t = r0*m1 + r1*m0 <= 2*m0*m1 - m0 - m1
                does not overflow
            */
            FLINT_ASSERT(FLINT_BIT_COUNT(M[1]) < FLINT_BITS);
            umul_ppmm(t[1], t[0], r0, m1);
            umul_ppmm(u[1], u[0], r1, m0);
            add_ssaaaa(t[1], t[0], t[1], t[0], u[1], u[0]);

            /* t = t mod M */
            if (t[1] > M[1] || (t[1] == M[1] && t[0] > M[0]))
                sub_ddmmss(t[1], t[0], t[1], t[0], M[1], M[0]);

            FLINT_ASSERT(t[1] < M[1] || (t[1] == M[1] && t[0] < M[0]));

            sub_ddmmss(u[1], u[0], M[1], M[0], t[1], t[0]);
            if (u[1] < t[1] || (u[1] == t[1] && u[0] < t[0]))
                fmpz_neg_uiui(&Crows[i][j], u[1], u[0]);
            else
                fmpz_set_uiui(&Crows[i][j], t[1], t[0]);
        }
    }
    else
    {
        mp_ptr M, Ns, T, U;
        mp_size_t Msize, Nsize;
        mp_limb_t cy, ri;

        M = FLINT_ARRAY_ALLOC(num_primes + 1, mp_limb_t);

        M[0] = primes[0];
        Msize = 1;
        for (i = 1; i < num_primes; i++)
        {
            FLINT_ASSERT(Msize > 0);
            M[Msize] = cy = mpn_mul_1(M, M, Msize, primes[i]);
            Msize += (cy != 0);
        }

        /* We add terms with Msize + 1 limbs, with one extra limb for the
           carry accumulation. todo: reduce Nsize by 1 when the carries
           do not require an extra limb. */
        Nsize = Msize + 2;

        Ns = FLINT_ARRAY_ALLOC(Nsize*num_primes, mp_limb_t);
        T = FLINT_ARRAY_ALLOC(Nsize, mp_limb_t);
        U = FLINT_ARRAY_ALLOC(Nsize, mp_limb_t);

        for (i = 0; i < num_primes; i++)
        {
            Ns[i*Nsize + (Nsize - 1)] = 0;
            Ns[i*Nsize + (Nsize - 2)] = 0;
            mpn_divrem_1(Ns + i * Nsize, 0, M, Msize, primes[i]);
            ri = mpn_mod_1(Ns + i * Nsize, Msize, primes[i]);
            ri = n_invmod(ri, primes[i]);
            FLINT_ASSERT(Msize > 0);
            Ns[i*Nsize + Msize] = mpn_mul_1(Ns + i*Nsize, Ns + i*Nsize, Msize, ri);
        }

        for (i = Cstartrow; i < Cstoprow; i++)
        for (j = 0; j < n; j++)
        {
            ri = nmod_mat_entry(mod_C[0], i, j);
            FLINT_ASSERT(Nsize > 1);
            T[Nsize - 1] = mpn_mul_1(T, Ns, Nsize - 1, ri);

            for (l = 1; l < num_primes; l++)
            {
                ri = nmod_mat_entry(mod_C[l], i, j);
                T[Nsize - 1] += mpn_addmul_1(T, Ns + l*Nsize, Nsize - 1, ri);
            }

            mpn_tdiv_qr(U, T, 0, T, Nsize, M, Msize);
            mpn_sub_n(U, M, T, Msize);

            if (mpn_cmp(U, T, Msize) < 0)
            {
                fmpz_set_ui_array(&Crows[i][j], U, Msize);
                fmpz_neg(&Crows[i][j], &Crows[i][j]);
            }
            else
            {
                fmpz_set_ui_array(&Crows[i][j], T, Msize);
            }
        }

        flint_free(M);
        flint_free(Ns);
        flint_free(T);
        flint_free(U);
    }
}


void _fmpz_mat_mul_multi_mod(
    fmpz_mat_t C,
    const fmpz_mat_t A,
    const fmpz_mat_t B,
    flint_bitcnt_t bits)
{
    slong i, start, stop;
    slong m, k, n;
    flint_bitcnt_t primes_bits;
    _worker_arg mainarg;
    _worker_arg * args;
    fmpz_comb_t comb;
    slong limit, num_workers;
    thread_pool_handle * handles;

    mainarg.m = m = A->r;
    mainarg.k = k = A->c;
    mainarg.n = n = B->c;
    mainarg.Arows = A->rows;
    mainarg.Brows = B->rows;
    mainarg.Crows = C->rows;

    if (m < 1 || n < 1 || k < 1)
    {
        fmpz_mat_zero(C);
        return;
    }

    primes_bits = NMOD_MAT_OPTIMAL_MODULUS_BITS;

    if (bits < primes_bits)
    {
        primes_bits = bits;
        mainarg.num_primes = 1;
    }
    else
    {
        /* Round up in the division */
        mainarg.num_primes = (bits + primes_bits - 1) / primes_bits;
    }

    /* Initialize */
    mainarg.primes = FLINT_ARRAY_ALLOC(mainarg.num_primes, mp_limb_t);
    mainarg.primes[0] = n_nextprime(UWORD(1) << primes_bits, 0);
    for (i = 1; i < mainarg.num_primes; i++)
        mainarg.primes[i] = n_nextprime(mainarg.primes[i-1], 0);

    mainarg.mod_A = FLINT_ARRAY_ALLOC(mainarg.num_primes, nmod_mat_t);
    mainarg.mod_B = FLINT_ARRAY_ALLOC(mainarg.num_primes, nmod_mat_t);
    mainarg.mod_C = FLINT_ARRAY_ALLOC(mainarg.num_primes, nmod_mat_t);
    for (i = 0; i < mainarg.num_primes; i++)
    {
        nmod_mat_init(mainarg.mod_A[i], A->r, A->c, mainarg.primes[i]);
        nmod_mat_init(mainarg.mod_B[i], B->r, B->c, mainarg.primes[i]);
        nmod_mat_init(mainarg.mod_C[i], C->r, C->c, mainarg.primes[i]);
    }

    if (mainarg.num_primes > 200)
    {
        /* use comb */
        fmpz_comb_init(comb, mainarg.primes, mainarg.num_primes);
        mainarg.comb = comb;
    }
    else
    {
        /* don't use comb */
        mainarg.comb = NULL;
    }

    /* limit on the number of threads */
    limit = FLINT_MAX(k, n);
    limit = FLINT_MIN(limit, m);
    limit = (limit - 16 < 0) ? 0 : (limit - 16)/16;

    /* mod */
    if (limit < 2)
    {
mod_single:
        mainarg.Astartrow = 0;
        mainarg.Astoprow = m;
        mainarg.Bstartrow = 0;
        mainarg.Bstoprow = k;
        _mod_worker(&mainarg);
    }
    else
    {
        num_workers = flint_request_threads(&handles, limit);
        if (num_workers < 1)
        {
            flint_give_back_threads(handles, num_workers);
            goto mod_single;
        }

        args = FLINT_ARRAY_ALLOC(num_workers, _worker_arg);
        for (start = 0, i = 0; i < num_workers; start = stop, i++)
        {
            args[i] = mainarg;
            stop = _peicewise_linear_inverse2(m, k, k, n, i + 1, num_workers + 1);
            _distribute_rows2(start, stop,
                              &args[i].Astartrow, &args[i].Astoprow, m,
                              &args[i].Bstartrow, &args[i].Bstoprow, k);
        }

        _distribute_rows2(start, m + k,
                          &mainarg.Astartrow, &mainarg.Astoprow, m,
                          &mainarg.Bstartrow, &mainarg.Bstoprow, k);

        for (i = 0; i < num_workers; i++)
            thread_pool_wake(global_thread_pool, handles[i], 0, _mod_worker, &args[i]);
        _mod_worker(&mainarg);
        for (i = 0; i < num_workers; i++)
            thread_pool_wait(global_thread_pool, handles[i]);

        flint_give_back_threads(handles, num_workers);
        flint_free(args);
    }

    /* mul */
    for (i = 0; i < mainarg.num_primes; i++)
        nmod_mat_mul(mainarg.mod_C[i], mainarg.mod_A[i], mainarg.mod_B[i]);

    /* crt */
    if (limit < 2)
    {
crt_single:
        mainarg.Cstartrow = 0;
        mainarg.Cstoprow = m;
        _crt_worker(&mainarg);
    }
    else
    {
        num_workers = flint_request_threads(&handles, limit);
        if (num_workers < 1)
        {
            flint_give_back_threads(handles, num_workers);
            goto crt_single;
        }

        args = FLINT_ARRAY_ALLOC(num_workers, _worker_arg);
        for (start = 0, i = 0; i < num_workers; start = stop, i++)
        {
            args[i] = mainarg;
            stop = (i + 1)*m/(num_workers + 1);
            args[i].Cstartrow = start;
            args[i].Cstoprow = stop;
        }

        mainarg.Cstartrow = start;
        mainarg.Cstoprow = m;

        for (i = 0; i < num_workers; i++)
            thread_pool_wake(global_thread_pool, handles[i], 0, _crt_worker, &args[i]);
        _crt_worker(&mainarg);
        for (i = 0; i < num_workers; i++)
            thread_pool_wait(global_thread_pool, handles[i]);

        flint_give_back_threads(handles, num_workers);
        flint_free(args);
    }

    /* Cleanup */
    if (mainarg.comb != NULL)
        fmpz_comb_clear(comb);

    for (i = 0; i < mainarg.num_primes; i++)
    {
        nmod_mat_clear(mainarg.mod_A[i]);
        nmod_mat_clear(mainarg.mod_B[i]);
        nmod_mat_clear(mainarg.mod_C[i]);
    }

    flint_free(mainarg.mod_A);
    flint_free(mainarg.mod_B);
    flint_free(mainarg.mod_C);
    flint_free(mainarg.primes);
}

void
fmpz_mat_mul_multi_mod(fmpz_mat_t C, const fmpz_mat_t A, const fmpz_mat_t B)
{
    slong A_bits;
    slong B_bits;

    A_bits = fmpz_mat_max_bits(A);
    B_bits = fmpz_mat_max_bits(B);

    _fmpz_mat_mul_multi_mod(C, A, B, FLINT_ABS(A_bits) + FLINT_ABS(B_bits) \
        + FLINT_BIT_COUNT(A->c) + 1);
}
