/*
    Copyright (C) 2019 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mod_mpoly.h"
#include "profiler.h"

int
main(void)
{
    slong i, j, k;
    FLINT_TEST_INIT(state);

    flint_printf("crt....");
    fflush(stdout);

    /* test internal interface */
    {
        fmpz_mod_poly_multi_crt_t P;
        fmpz_mod_poly_struct * moduli, * inputs, * outputs;
        slong moduli_count = 3000;
        fmpz_mod_ctx_t ctx;
        fmpz_t modulus;
timeit_t timer;

        fmpz_init_set_ui(modulus, 3);
        fmpz_pow_ui(modulus, modulus, 8);
        fmpz_nextprime(modulus, modulus, 1);

        fmpz_mod_ctx_init(ctx, modulus);

        moduli = FLINT_ARRAY_ALLOC(moduli_count, fmpz_mod_poly_struct);
        inputs = FLINT_ARRAY_ALLOC(moduli_count, fmpz_mod_poly_struct);
        outputs = FLINT_ARRAY_ALLOC(moduli_count, fmpz_mod_poly_struct);

        for (k = 0; k < moduli_count; k++)
        {
            fmpz_mod_poly_init(moduli + k, ctx);
            fmpz_mod_poly_init(inputs + k, ctx);
            fmpz_mod_poly_init(outputs + k, ctx);

            fmpz_mod_poly_set_coeff_ui(moduli + k, 1, 1, ctx);
            fmpz_mod_poly_set_coeff_ui(moduli + k, 0, k, ctx);

            fmpz_mod_poly_set_coeff_ui(inputs + k, 0, (k*k)^k, ctx);
        }

        fmpz_mod_poly_multi_crt_init(P, ctx);
        if (!fmpz_mod_poly_multi_crt_precompute(P, moduli, moduli_count, ctx))
        {
            flint_printf("FAIL: Check simple example\n");
            flint_abort();            
        }

        FLINT_ASSERT(_fmpz_mod_poly_multi_crt_local_size(P) <= moduli_count);
timeit_start(timer);
        for (k = 0; k < 100; k++)
        {
            _fmpz_mod_poly_multi_crt_run(outputs, P, inputs, ctx);
        }
timeit_stop(timer);
flint_printf("time: %wd\n");

        for (k = 0; k < moduli_count; k++)
        {
            fmpz_mod_poly_clear(moduli + k, ctx);
            fmpz_mod_poly_clear(inputs + k, ctx);
            fmpz_mod_poly_clear(outputs + k, ctx);
        }

        flint_free(moduli);
        flint_free(inputs);
        flint_free(outputs);

        fmpz_mod_poly_multi_crt_clear(P, ctx);

        fmpz_mod_ctx_clear(ctx);
        fmpz_clear(modulus);
    }

    /* test flat interface */
    for (i = 0; i < 5 * flint_test_multiplier(); i++)
    {
        fmpz_mod_poly_multi_crt_t P;
        fmpz_mod_poly_t t, p;
        slong total_degree, moduli_length, moduli_count;
        fmpz_mod_poly_struct * moduli, * inputs, * rems;
        fmpz_mod_poly_t output;
        fmpz_mod_ctx_t ctx;

        fmpz_mod_ctx_init_rand_bits_prime(ctx, state, 100);

        fmpz_mod_poly_init(t, ctx);
        fmpz_mod_poly_init(p, ctx);
        fmpz_mod_poly_init(output, ctx);

        fmpz_mod_poly_multi_crt_init(P, ctx);

        for (j = 0; j < 4; j++)
        {
            moduli_length = n_randint(state, 20) + 1;
            moduli_count = n_randint(state, 50) + 1;

            moduli = FLINT_ARRAY_ALLOC(moduli_count, fmpz_mod_poly_struct);
            inputs = FLINT_ARRAY_ALLOC(moduli_count, fmpz_mod_poly_struct);
            rems = FLINT_ARRAY_ALLOC(moduli_count, fmpz_mod_poly_struct);
            for (k = 0; k < moduli_count; k++)
            {
                fmpz_mod_poly_init(moduli + k, ctx);
                fmpz_mod_poly_init(inputs + k, ctx);
                fmpz_mod_poly_init(rems + k, ctx);
                fmpz_mod_poly_randtest(moduli + k, state, moduli_length, ctx);
                fmpz_mod_poly_randtest(inputs + k, state, moduli_length, ctx);
            }

            if (fmpz_mod_poly_multi_crt_precompute(P, moduli, moduli_count, ctx))
            {
                fmpz_mod_poly_multi_crt_precomp(output, P, inputs, ctx);

                total_degree = 0;
                for (k = 0; k < moduli_count; k++)
                {
                    total_degree += fmpz_mod_poly_degree(moduli + k, ctx);

                    fmpz_mod_poly_sub(t, output, inputs + k, ctx);
                    fmpz_mod_poly_rem(t, t, moduli + k, ctx);
                    if (!fmpz_mod_poly_is_zero(t, ctx))
                    {
                        flint_printf("FAIL: Check remainder flat\n");
                        flint_printf("i = %wd, j = %wd, k = %wd\n", i, j, k);
                        flint_abort();
                    }
                }
                if (fmpz_mod_poly_degree(output, ctx) >= total_degree)
                {
                    flint_printf("FAIL: Check output degree flat\n");
                    flint_printf("i = %wd, j = %wd, k = %wd\n", i, j, k);
                    flint_abort();
                }

                fmpz_mod_poly_multi_mod(rems, moduli, output, moduli_count, ctx);
                for (k = 0; k < moduli_count; k++)
                {
                    fmpz_mod_poly_sub(t, rems + k, inputs + k, ctx);
                    fmpz_mod_poly_rem(t, t, moduli + k, ctx);
                    if (!fmpz_mod_poly_is_zero(t, ctx))
                    {
                        flint_printf("FAIL: Check multi remainder\n");
                        flint_printf("i = %wd, j = %wd, k = %wd\n", i, j, k);
                        flint_abort();
                    }
                }
            }
            else
            {
                /* check if it was ok to fail on these moduli */
                int ok = 0;

                fmpz_mod_poly_one(p, ctx);
                for (k = 0; !ok && k < moduli_count; k++)
                {
                    fmpz_mod_poly_mul(p, p, moduli + k, ctx);
                    if (fmpz_mod_poly_degree(moduli + k, ctx) < 1)
                        ok = 1;
                }

                for (k = 0; !ok && k < moduli_count; k++)
                {
                    fmpz_mod_poly_div(t, p, moduli + k, ctx);
                    fmpz_mod_poly_gcd(t, t, moduli + k, ctx);
                    if (fmpz_mod_poly_degree(t, ctx) > 0)
                        ok = 1;
                }

                if (!ok)
                {
                    flint_printf("FAIL: Check crt failure flat\n");
                    flint_printf("i = %wd, j = %wd\n", i, j);
                    flint_abort();                    
                }
            }

            for (k = 0; k < moduli_count; k++)
            {
                fmpz_mod_poly_clear(moduli + k, ctx);
                fmpz_mod_poly_clear(inputs + k, ctx);
                fmpz_mod_poly_clear(rems + k, ctx);
            }
            flint_free(moduli);
            flint_free(inputs);
            flint_free(rems);
        }

        fmpz_mod_poly_clear(t, ctx);
        fmpz_mod_poly_clear(p, ctx);
        fmpz_mod_poly_clear(output, ctx);

        fmpz_mod_poly_multi_crt_clear(P, ctx);

        fmpz_mod_ctx_clear(ctx);
    }

    /* test lazy interface */
    for (i = 0; i < 5 * flint_test_multiplier(); i++)
    {
        fmpz_mod_poly_t t, p;
        slong total_degree, moduli_length, moduli_count;
        fmpz_mod_poly_struct * moduli, * inputs;
        fmpz_mod_poly_t output;
        fmpz_mod_ctx_t ctx;

        fmpz_mod_ctx_init_rand_bits_prime(ctx, state, 100);

        fmpz_mod_poly_init(t, ctx);
        fmpz_mod_poly_init(p, ctx);
        fmpz_mod_poly_init(output, ctx);

        for (j = 0; j < 4; j++)
        {
            moduli_length = n_randint(state, 20) + 1;
            moduli_count = n_randint(state, 50) + 1;

            moduli = FLINT_ARRAY_ALLOC(moduli_count, fmpz_mod_poly_struct);
            inputs = FLINT_ARRAY_ALLOC(moduli_count, fmpz_mod_poly_struct);
            for (k = 0; k < moduli_count; k++)
            {
                fmpz_mod_poly_init(moduli + k, ctx);
                fmpz_mod_poly_init(inputs + k, ctx);
                fmpz_mod_poly_randtest(moduli + k, state, moduli_length, ctx);
                fmpz_mod_poly_randtest(inputs + k, state, moduli_length, ctx);
            }

            if (fmpz_mod_poly_multi_crt(output, moduli, inputs, moduli_count, ctx))
            {
                total_degree = 0;
                for (k = 0; k < moduli_count; k++)
                {
                    total_degree += fmpz_mod_poly_degree(moduli + k, ctx);

                    fmpz_mod_poly_sub(t, output, inputs + k, ctx);
                    fmpz_mod_poly_rem(t, t, moduli + k, ctx);
                    if (!fmpz_mod_poly_is_zero(t, ctx))
                    {
                        flint_printf("FAIL: Check remainder lazy\n");
                        flint_printf("i = %wd, j = %wd, k = %wd\n", i, j, k);
                        flint_abort();
                    }
                }
                if (fmpz_mod_poly_degree(output, ctx) >= total_degree)
                {
                    flint_printf("FAIL: Check output degree lazy\n");
                    flint_printf("i = %wd, j = %wd, k = %wd\n", i, j, k);
                    flint_abort();
                }
            }
            else
            {
                /* check if it was ok to fail on these moduli */
                int ok = 0;

                fmpz_mod_poly_one(p, ctx);
                for (k = 0; !ok && k < moduli_count; k++)
                {
                    fmpz_mod_poly_mul(p, p, moduli + k, ctx);
                    if (fmpz_mod_poly_degree(moduli + k, ctx) < 1)
                        ok = 1;
                }

                for (k = 0; !ok && k < moduli_count; k++)
                {
                    fmpz_mod_poly_div(t, p, moduli + k, ctx);
                    fmpz_mod_poly_gcd(t, t, moduli + k, ctx);
                    if (fmpz_mod_poly_degree(t, ctx) > 0)
                        ok = 1;
                }

                if (!ok)
                {
                    flint_printf("FAIL: Check crt failure lazy\n");
                    flint_printf(" i = %wd, j = %wd\n", i, j);
                    flint_abort();                    
                }
            }

            for (k = 0; k < moduli_count; k++)
            {
                fmpz_mod_poly_clear(moduli + k, ctx);
                fmpz_mod_poly_clear(inputs + k, ctx);
            }
            flint_free(moduli);
            flint_free(inputs);
        }

        fmpz_mod_poly_clear(t, ctx);
        fmpz_mod_poly_clear(p, ctx);
        fmpz_mod_poly_clear(output, ctx);

        fmpz_mod_ctx_clear(ctx);
    }

    flint_printf("PASS\n");
    FLINT_TEST_CLEANUP(state);

    return 0;
}
