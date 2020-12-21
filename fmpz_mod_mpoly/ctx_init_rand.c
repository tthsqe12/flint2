/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mod_mpoly.h"

void fmpz_mod_mpoly_ctx_init_rand(fmpz_mod_mpoly_ctx_t ctx,
                    flint_rand_t state, slong max_nvars, const fmpz_t modulus)
{
    mpoly_ctx_init_rand(ctx->minfo, state, max_nvars);
    fmpz_mod_ctx_init(ctx->ffinfo, modulus);
}

void fmpz_mod_mpoly_ctx_init_rand_bits_prime(fmpz_mod_mpoly_ctx_t ctx,
                  flint_rand_t state, slong max_nvars, flint_bitcnt_t max_bits)
{
    fmpz_t m;
    fmpz_init(m);
    fmpz_randtest_unsigned(m, state, max_bits);
    fmpz_nextprime(m, m, 0);
    fmpz_mod_mpoly_ctx_init_rand(ctx, state, max_nvars, m);
    fmpz_clear(m);
}

void fmpz_mod_mpoly_ctx_init_rand_bits(fmpz_mod_mpoly_ctx_t ctx,
                  flint_rand_t state, slong max_nvars, flint_bitcnt_t max_bits)
{
    fmpz_t m;
    fmpz_init(m);
    fmpz_randtest_unsigned(m, state, max_bits);
    fmpz_add_ui(m, m, 2);
    fmpz_mod_mpoly_ctx_init_rand(ctx, state, max_nvars, m);
    fmpz_clear(m);
}
