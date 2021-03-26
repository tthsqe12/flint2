/*
    Copyright (C) 2010,2011,2018 Fredrik Johansson
    Copyright (C) 2021 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mat.h"

#if 0

static int
fmpz_get_sgnbit_mpn2(mp_ptr r, const fmpz_t x)
{
    if (!COEFF_IS_MPZ(*x))
    {
        slong v = *x;
        r[0] = FLINT_ABS(v);
        r[1] = 0;
        return v < 0;
    }
    else
    {
        __mpz_struct * p = COEFF_TO_PTR(*x);
        slong sz = p->_mp_size;
        r[0] = p->_mp_d[0];
        if (sz == 2 || sz == -2)
            r[1] = p->_mp_d[1];
        else
            r[1] = 0;
        return sz < 0;
    }
}

#define nn_mul_2x1(r2, r1, r0, a1, a0, b0)                  \
    do {                                                    \
        mp_limb_t t1;                                       \
        umul_ppmm(r1, r0, a0, b0);                          \
        umul_ppmm(r2, t1, a1, b0);                          \
        add_ssaaaa(r2, r1, r2, r1, 0, t1);                  \
    } while (0)

#define nn_mul_2x2(r3, r2, r1, r0, a1, a0, b1, b0)          \
    do {                                                    \
        mp_limb_t t1, t2, t3;                               \
        umul_ppmm(r1, r0, a0, b0);                          \
        umul_ppmm(r2, t1, a1, b0);                          \
        add_ssaaaa(r2, r1, r2, r1, 0, t1);                  \
        umul_ppmm(t1, t2, a0, b1);                          \
        umul_ppmm(r3, t3, a1, b1);                          \
        add_ssaaaa(r3, t1, r3, t1, 0, t3);                  \
        add_ssaaaa(r2, r1, r2, r1, t1, t2);                 \
        r3 += r2 < t1;                                      \
} while (0)

FLINT_DLL void
fmpz_mat_mul_4(fmpz_mat_t C, const fmpz_mat_t A, const fmpz_mat_t B)
{
    slong ar, ac, br, bc;
    slong i, j, k;
    mp_ptr AL, BL;
    char *AS, *BS;
    TMP_INIT;

    ar = fmpz_mat_nrows(A);
    ac = fmpz_mat_ncols(A);
    br = fmpz_mat_nrows(B);
    bc = fmpz_mat_ncols(B);

    TMP_START;

    AL = TMP_ALLOC(2 * sizeof(mp_limb_t) * ar * ac);
    BL = TMP_ALLOC(2 * sizeof(mp_limb_t) * br * bc);
    AS = TMP_ALLOC(sizeof(char) * ar * ac);
    BS = TMP_ALLOC(sizeof(char) * br * bc);

    for (i = 0; i < ar; i++)
        for (j = 0; j < ac; j++)
            AS[i * ac + j] = fmpz_get_sgnbit_mpn2(AL + 2 * i * ac + 2 * j,
                fmpz_mat_entry(A, i, j));

    for (i = 0; i < br; i++)
        for (j = 0; j < bc; j++)
            BS[j * br + i] = fmpz_get_sgnbit_mpn2(BL + 2 * j * br + 2 * i,
                fmpz_mat_entry(B, i, j));

    for (i = 0; i < ar; i++)
    {
        for (j = 0; j < bc; j++)
        {
            mp_limb_t s[4];
            mp_limb_t t[4];
            mp_limb_t u[4];

            flint_mpn_zero(s, 4);
            flint_mpn_zero(t, 4);
            flint_mpn_zero(u, 4);

            for (k = 0; k < br; k++)
            {
                mp_limb_t ah, al, bh, bl;
                mp_srcptr aptr, bptr;

                aptr = AL + 2 * i * ac + 2 * k;
                bptr = BL + 2 * j * br + 2 * k;

                al = aptr[0];
                ah = aptr[1];
                bl = bptr[0],
                bh = bptr[1];

                if (ah == 0 && bh == 0)
                {
                    if (al == 0 || bl == 0)
                        continue;

                    umul_ppmm(t[1], t[0], al, bl);

                    if (AS[i * ac + k] == BS[j * br + k])
                        add_sssaaaaaa(u[2], u[1], u[0], u[2], u[1], u[0], 0, t[1], t[0]);
                    else
                        sub_dddmmmsss(u[2], u[1], u[0], u[2], u[1], u[0], 0, t[1], t[0]);
                }
                else
                {
                    if (ah == 0)
                    {
                        nn_mul_2x1(t[2], t[1], t[0], bh, bl, al);
                        t[3] = 0;
                    }
                    else if (bh == 0)
                    {
                        nn_mul_2x1(t[2], t[1], t[0], ah, al, bl);
                        t[3] = 0;
                    }
                    else
                    {
                        nn_mul_2x2(t[3], t[2], t[1], t[0], ah, al, bh, bl);
                    }

                    if (AS[i * ac + k] == BS[j * br + k])
                        mpn_add_n(s, s, t, 4);
                    else
                        mpn_sub_n(s, s, t, 4);
                }
            }

            if (((mp_limb_signed_t) u[2]) >= 0)
            {
                s[3] += mpn_add_n(s, s, u, 3);
            }
            else
            {
                sub_dddmmmsss(u[2], u[1], u[0], 0, 0, 0, u[2], u[1], u[0]);
                s[3] -= mpn_sub_n(s, s, u, 3);
            }

            if (((mp_limb_signed_t) s[3]) >= 0)
            {
                fmpz_set_ui_array(fmpz_mat_entry(C, i, j), s, 4);
            }
            else
            {
                mpn_neg_n(s, s, 4);
                fmpz_set_ui_array(fmpz_mat_entry(C, i, j), s, 4);
                fmpz_neg(fmpz_mat_entry(C, i, j), fmpz_mat_entry(C, i, j));
            }
        }
    }

    TMP_END;
}


#else



/* loop "warmup" time */
#define FMPZ_MAT_MUL_4_BRANCHLESS_CUTOFF 16


void fmpz_get_signed_uiui(mp_limb_t * hi, mp_limb_t * lo, const fmpz_t x)
{
    ulong r0, r1, s;

    if (!COEFF_IS_MPZ(*x))
    {
        r0 = *x;
        r1 = FLINT_SIGN_EXT(r0);
    }
    else
    {
        __mpz_struct * p = COEFF_TO_PTR(*x);
        s = -(ulong)(p->_mp_size < 0);
        r0 = p->_mp_d[0];
        if (p->_mp_size == 2 || p->_mp_size == -2)
            r1 = p->_mp_d[1];
        else
            r1 = 0;

        sub_ddmmss(r1, r0, r1^s, r0^s, s, s);
    }

    *lo = r0;
    *hi = r1;
}


void _do_signed_row_branchfull(
    fmpz * CR,
    const mp_limb_t * AR,
    const mp_limb_t * B,
    slong br,
    slong bc)
{
    slong j, k, l;
    mp_limb_t s[4], t3, t2, t1, t0, w3, w2, w1, w0;
    mp_limb_t A0, A1, B0, B1;
    mp_limb_t u2, u1, u0;

    for (j = 0, l = 0; j < bc; j++)
    {
        t3 = t2 = t1 = t0 = 0;
        u2 = u1 = u0 = 0;

        for (k = 0; k < br; k++, l++)
        {
            A0 = AR[2*k + 0];
            A1 = AR[2*k + 1];
            B0 = B[2*l + 0];
            B1 = B[2*l + 1];

            if (FLINT_SIGN_EXT(A0) == A1 && FLINT_SIGN_EXT(B0) == B1)
            {
                smul_ppmm(w1, w0, B0, A0);
                add_sssaaaaaa(u2, u1, u0, u2, u1, u0,
                              FLINT_SIGN_EXT(w1), w1, w0);
            }
            else
            {
                sub_ddmmss(t3, t2, t3, t2, 0, FLINT_SIGN_EXT(A1)&B0);
                sub_ddmmss(t3, t2, t3, t2, 0, FLINT_SIGN_EXT(B1)&A0);

                smul_ppmm(w3, w2, B1, A1);
                add_ssaaaa(t3, t2, t3, t2, w3, w2);

                umul_ppmm(w1, w0, B0, A0);
                add_sssaaaaaa(u2, u1, u0, u2, u1, u0, UWORD(0), w1, w0);

                umul_ppmm(w2, w1, A1, B0);
                add_sssaaaaaa(t3, t2, t1, t3, t2, t1, UWORD(0), w2, w1);

                umul_ppmm(w2, w1, B1, A0);
                add_sssaaaaaa(t3, t2, t1, t3, t2, t1, UWORD(0), w2, w1);
            }
        }

        add_ssssaaaaaaaa(s[3], s[2], s[1], s[0], t3, t2, t1, t0,
                         FLINT_SIGN_EXT(u2), u2, u1, u0);

        fmpz_set_signed_ui_array(CR + j, s, 4);
    }
}

void _do_signed_row_branchless(
    fmpz * CR,
    const mp_limb_t * AR,
    const mp_limb_t * B,
    slong br,
    slong bc)
{
    slong j, k, l;
    mp_limb_t s[4], t3, t2, t1, t0, w3, w2, w1, w0;
    mp_limb_t A0, A1, B0, B1;
    mp_limb_t u3, u2, u1;

    for (j = 0, l = 0; j < bc; j++)
    {
        t3 = t2 = t1 = t0 = 0;
        u3 = u2 = u1 = 0;

        for (k = 0; k < br; k++, l++)
        {
            A0 = AR[2*k + 0];
            A1 = AR[2*k + 1];
            B0 = B[2*l + 0];
            B1 = B[2*l + 1];

            sub_ddmmss(u3, u2, u3, u2, 0, FLINT_SIGN_EXT(A1)&B0);
            sub_ddmmss(t3, t2, t3, t2, 0, FLINT_SIGN_EXT(B1)&A0);

            umul_ppmm(w1, w0, B0, A0);
            smul_ppmm(w3, w2, B1, A1);
            add_ssssaaaaaaaa(u3, u2, u1, t0, u3, u2, u1, t0, w3, w2, w1, w0);

            umul_ppmm(w2, w1, A1, B0);
            add_sssaaaaaa(t3, t2, t1, t3, t2, t1, 0, w2, w1);

            umul_ppmm(w2, w1, B1, A0);
            add_sssaaaaaa(u3, u2, u1, u3, u2, u1, 0, w2, w1);
        }

        s[0] = t0;
        add_sssaaaaaa(s[3], s[2], s[1], t3, t2, t1, u3, u2, u1);

        fmpz_set_signed_ui_array(CR + j, s, 4);
    }
}

static void _do_signed_row(
    fmpz * CR,
    const mp_limb_t * AR,
    const mp_limb_t * B,
    slong br,
    slong bc)
{
    if (br < FMPZ_MAT_MUL_4_BRANCHLESS_CUTOFF)
        _do_signed_row_branchfull(CR, AR, B, br, bc);
    else
        _do_signed_row_branchless(CR, AR, B, br, bc);
}


void _do_unsigned_row_branchfull(
    fmpz * CR,
    const mp_limb_t * AR,
    const mp_limb_t * B,
    slong br,
    slong bc)
{
    slong j, k, l;
    mp_limb_t s[4], t3, t2, t1, t0, w3, w2, w1, w0;
    mp_limb_t A0, A1, B0, B1;
    mp_limb_t u2, u1, u0;

    for (j = 0, l = 0; j < bc; j++)
    {
        t3 = t2 = t1 = t0 = 0;
        u2 = u1 = u0 = 0;

        for (k = 0; k < br; k++, l++)
        {
            A0 = AR[2*k + 0];
            A1 = AR[2*k + 1];
            B0 = B[2*l + 0];
            B1 = B[2*l + 1];

            umul_ppmm(w1, w0, B0, A0);
            add_sssaaaaaa(u2, u1, u0, u2, u1, u0, UWORD(0), w1, w0);

            if (0 != A1 || 0 != B1)
            {
                umul_ppmm(w3, w2, B1, A1);
                add_ssaaaa(t3, t2, t3, t2, w3, w2);

                add_sssaaaaaa(u2, u1, u0, u2, u1, u0, UWORD(0), w1, w0);

                umul_ppmm(w2, w1, A1, B0);
                add_sssaaaaaa(t3, t2, t1, t3, t2, t1, UWORD(0), w2, w1);

                umul_ppmm(w2, w1, B1, A0);
                add_sssaaaaaa(t3, t2, t1, t3, t2, t1, UWORD(0), w2, w1);
            }
        }

        add_ssssaaaaaaaa(s[3], s[2], s[1], s[0], t3, t2, t1, t0,
                         UWORD(0), u2, u1, u0);

        fmpz_set_ui_array(CR + j, s, 4);
    }
}

void _do_unsigned_row_branchless(
    fmpz * CR,
    const mp_limb_t * AR,
    const mp_limb_t * B,
    slong br,
    slong bc)
{
    slong j, k, l;
    mp_limb_t s[4], t3, t2, t1, t0, w3, w2, w1, w0;
    mp_limb_t A0, A1, B0, B1;
    mp_limb_t u3, u2, u1;

    for (j = 0, l = 0; j < bc; j++)
    {
        t3 = t2 = t1 = t0 = 0;
        u3 = u2 = u1 = 0;

        for (k = 0; k < br; k++, l++)
        {
            A0 = AR[2*k + 0];
            A1 = AR[2*k + 1];
            B0 = B[2*l + 0];
            B1 = B[2*l + 1];

            umul_ppmm(w1, w0, B0, A0);
            umul_ppmm(w3, w2, B1, A1);
            add_ssssaaaaaaaa(u3, u2, u1, t0, u3, u2, u1, t0, w3, w2, w1, w0);

            umul_ppmm(w2, w1, A1, B0);
            add_sssaaaaaa(t3, t2, t1, t3, t2, t1, UWORD(0), w2, w1);

            umul_ppmm(w2, w1, B1, A0);
            add_sssaaaaaa(u3, u2, u1, u3, u2, u1, UWORD(0), w2, w1);
        }

        s[0] = t0;
        add_sssaaaaaa(s[3], s[2], s[1], t3, t2, t1, u3, u2, u1);

        fmpz_set_ui_array(CR + j, s, 4);
    }
}

static void _do_unsigned_row(
    fmpz * CR,
    const mp_limb_t * AR,
    const mp_limb_t * B,
    slong br,
    slong bc)
{
    if (br < FMPZ_MAT_MUL_4_BRANCHLESS_CUTOFF)
        _do_unsigned_row_branchfull(CR, AR, B, br, bc);
    else
        _do_unsigned_row_branchless(CR, AR, B, br, bc);
}


typedef struct {
    slong Astartrow;
    slong Astoprow;
    slong Bstartcol;
    slong Bstopcol;
    slong br;
    slong bc;
    fmpz ** Crows;
    fmpz ** Arows;
    fmpz ** Brows;
    mp_limb_t * BL;
    int sign;
} _worker_arg;

static void _red_worker(void * varg)
{
    _worker_arg * arg = (_worker_arg *) varg;
    slong Bstartcol = arg->Bstartcol;
    slong Bstopcol = arg->Bstopcol;
    slong br = arg->br;
    fmpz ** Brows = arg->Brows;
    mp_limb_t * BL = arg->BL;
    int sign = arg->sign;
    slong i, j;

    if (sign)
    {
        for (j = Bstartcol; j < Bstopcol; j++)
            for (i = 0; i < br; i++)
                fmpz_get_signed_uiui(BL + 2*(j*br + i) + 1,
                                     BL + 2*(j*br + i) + 0, &Brows[i][j]);
    }
    else
    {
        for (j = Bstartcol; j < Bstopcol; j++)
            for (i = 0; i < br; i++)
                fmpz_get_uiui(BL + 2*(j*br + i) + 1,
                              BL + 2*(j*br + i) + 0, &Brows[i][j]);
    }
}

static void _mul_worker(void * varg)
{
    _worker_arg * arg = (_worker_arg *) varg;
    slong Astartrow = arg->Astartrow;
    slong Astoprow = arg->Astoprow;
    slong ac = arg->br;
    slong br = arg->br;
    slong bc = arg->bc;
    fmpz ** Crows = arg->Crows;
    fmpz ** Arows = arg->Arows;
    mp_limb_t * BL = arg->BL;
    int sign = arg->sign;
    mp_limb_t * AL;
    slong i, j;
    TMP_INIT;

    TMP_START;

    AL = TMP_ARRAY_ALLOC(2*ac, mp_limb_t);

    if (sign)
    {
        for (i = Astartrow; i < Astoprow; i++)
        {
            for (j = 0; j < ac; j++)
                fmpz_get_signed_uiui(AL + 2*j + 1, AL + 2*j, &Arows[i][j]);

            _do_signed_row(Crows[i], AL, BL, br, bc);
        }
    }
    else
    {
        for (i = Astartrow; i < Astoprow; i++)
        {
            for (j = 0; j < ac; j++)
                fmpz_get_uiui(AL + 2*j + 1, AL + 2*j, &Arows[i][j]);

            _do_unsigned_row(Crows[i], AL, BL, br, bc);
        }
    }

    TMP_END;
}


/*
    sign = 1:   max|A|, max|B| < 2^(2*FLINT_BITS - 1)
                max|C| < 2^(4*FLINT_BITS - 1)

    sign = 0:   all entries are >= 0 and
                max|A|, max|B| < 2^(2*FLINT_BITS)
                max|C| < 2^(4*FLINT_BITS)
*/
FLINT_DLL void fmpz_mat_mul_4(
    fmpz_mat_t C,
    const fmpz_mat_t A,
    const fmpz_mat_t B,
    int sign)
{
    slong i, ar, br, bc;
    _worker_arg mainarg;
    thread_pool_handle * handles;
    slong num_workers;
    _worker_arg * args;
    slong limit;
    TMP_INIT;

    TMP_START;

    br = fmpz_mat_nrows(B);
    bc = fmpz_mat_ncols(B);
    ar = fmpz_mat_nrows(A);

    /* limit on number of threads */
    limit = FLINT_MAX(bc, bc);
    limit = FLINT_MIN(limit, ar);
    limit = limit < 16 ? 0 : (ar - 16)/8;

    mainarg.Astartrow = 0;
    mainarg.Astoprow = ar;
    mainarg.Bstartcol = 0;
    mainarg.Bstopcol = bc;
    mainarg.br = br;
    mainarg.bc = bc;
    mainarg.Crows = C->rows;
    mainarg.Arows = A->rows;
    mainarg.Brows = B->rows;
    mainarg.BL = TMP_ARRAY_ALLOC(br*bc*2, mp_limb_t);
    mainarg.sign = sign;

    if (limit < 2)
    {
use_one_thread:
        _red_worker(&mainarg);
        _mul_worker(&mainarg);
        TMP_END;
        return;
    }

    num_workers = flint_request_threads(&handles, limit);
    if (num_workers < 1)
    {
        flint_give_back_threads(handles, num_workers);
        goto use_one_thread;
    }

    args = FLINT_ARRAY_ALLOC(num_workers, _worker_arg);

    for (i = 0; i < num_workers; i++)
    {
        args[i].Astartrow = (i + 0)*ar/(num_workers + 1);
        args[i].Astoprow = (i + 1)*ar/(num_workers + 1);
        args[i].Bstartcol = (i + 0)*bc/(num_workers + 1);
        args[i].Bstopcol = (i + 1)*bc/(num_workers + 1);
        args[i].br = mainarg.br;
        args[i].bc = mainarg.bc;
        args[i].Crows = mainarg.Crows;
        args[i].Arows = mainarg.Arows;
        args[i].Brows = mainarg.Brows;
        args[i].BL = mainarg.BL;
        args[i].sign = mainarg.sign;
    }

    i = num_workers;
    mainarg.Astartrow = (i + 0)*ar/(num_workers + 1);
    mainarg.Astoprow = (i + 1)*ar/(num_workers + 1);
    mainarg.Bstartcol = (i + 0)*bc/(num_workers + 1);
    mainarg.Bstopcol = (i + 1)*bc/(num_workers + 1);

    for (i = 0; i < num_workers; i++)
        thread_pool_wake(global_thread_pool, handles[i], 0, _red_worker, &args[i]);
    _red_worker(&mainarg);
    for (i = 0; i < num_workers; i++)
        thread_pool_wait(global_thread_pool, handles[i]);

    for (i = 0; i < num_workers; i++)
        thread_pool_wake(global_thread_pool, handles[i], 0, _mul_worker, &args[i]);
    _mul_worker(&mainarg);
    for (i = 0; i < num_workers; i++)
        thread_pool_wait(global_thread_pool, handles[i]);

    flint_give_back_threads(handles, num_workers);
    flint_free(args);

    TMP_END;
}

#endif
