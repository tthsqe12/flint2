/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fq_nmod_mpoly_factor.h"

static void fq_nmod_mpoly_factor_mul_mpoly_fmpz(
	fq_nmod_mpoly_factor_t fac,
	const fq_nmod_mpoly_t a,
	const fmpz_t e,
	const fq_nmod_mpoly_ctx_t ctx)
{
	if (fq_nmod_mpoly_is_fq_nmod(a, ctx))
	{
        fq_nmod_t t;
        fq_nmod_init(t, ctx->fqctx);
		fq_nmod_mpoly_get_fq_nmod(t, a, ctx);
		fq_nmod_pow(t, t, e, ctx->fqctx);
		fq_nmod_mul(fac->constant, fac->constant, t, ctx->fqctx);
        fq_nmod_clear(t, ctx->fqctx);
		return;
	}
	else
	{
		fq_nmod_mpoly_factor_append_fmpz(fac, a, e, ctx);
	}
}

static void fq_nmod_mpoly_factor_mul_mpoly_ui(
	fq_nmod_mpoly_factor_t fac,
	const fq_nmod_mpoly_t a,
	ulong e,
	const fq_nmod_mpoly_ctx_t ctx)
{
	if (fq_nmod_mpoly_is_fq_nmod(a, ctx))
	{
        fq_nmod_t t;
        fq_nmod_init(t, ctx->fqctx);
		fq_nmod_mpoly_get_fq_nmod(t, a, ctx);
		fq_nmod_pow_ui(t, t, e, ctx->fqctx);
		fq_nmod_mul(fac->constant, fac->constant, t, ctx->fqctx);
        fq_nmod_clear(t, ctx->fqctx);
		return;
	}
	else
	{
		fq_nmod_mpoly_factor_append_ui(fac, a, e, ctx);
	}
}

static void _fq_nmod_mpoly_univar_shift_right(
    fq_nmod_mpoly_univar_t A,
    const fmpz_t shift,
    const fq_nmod_mpoly_ctx_t ctx)
{
    slong i;
    for (i = 0; i < A->length; i++)
    {
        fmpz_sub(A->exps + i, A->exps + i, shift);
        FLINT_ASSERT(fmpz_sgn(A->exps + i) >= 0);
    }
}

static int _squarefree_factors(
	fq_nmod_mpoly_factor_t f,
	const fq_nmod_mpoly_t A,
	const fq_nmod_mpoly_ctx_t ctx)
{
    return 0;
}


int fq_nmod_mpoly_factor_squarefree(
    fq_nmod_mpoly_factor_t f,
    const fq_nmod_mpoly_t A,
    const fq_nmod_mpoly_ctx_t ctx)
{
    int success;
    slong i, j, v;
    fq_nmod_mpoly_t c;
    fq_nmod_mpoly_univar_t u;
    fq_nmod_mpoly_factor_t g, t;
    fmpz * var_powers;

    f->num = 0;
    if (fq_nmod_mpoly_is_fq_nmod(A, ctx))
    {
		fq_nmod_mpoly_get_fq_nmod(f->constant, A, ctx);
        return 1;
    }

    FLINT_ASSERT(A->length > 0);
    fq_nmod_set(f->constant, A->coeffs + 0, ctx->fqctx);
	fq_nmod_mpoly_factor_fit_length(f, 1, ctx);
	fq_nmod_mpoly_make_monic(f->poly + 0, A, ctx);
	fmpz_one(f->exp + 0);
	f->num = 1;

    fq_nmod_mpoly_factor_init(g, ctx);
    fq_nmod_mpoly_factor_init(t, ctx);
    fq_nmod_mpoly_univar_init(u, ctx);
    fq_nmod_mpoly_init(c, ctx);
    var_powers = _fmpz_vec_init(ctx->minfo->nvars);

    FLINT_ASSERT(fq_nmod_mpoly_factor_matches(A, f, ctx));

    for (v = 0; v < ctx->minfo->nvars; v++)
    {
        fq_nmod_swap(g->constant, f->constant, ctx->fqctx);
        g->num = 0;
        for (j = 0; j < f->num; j++)
        {
            FLINT_ASSERT(fmpz_is_one(f->exp + j));

            fq_nmod_mpoly_to_univar(u, f->poly + j, v, ctx);
            FLINT_ASSERT(u->length > 0);

            success = fq_nmod_mpoly_univar_content_mpoly(c, u, ctx);
            if (!success)
                goto cleanup;

            fq_nmod_mpoly_univar_divexact_mpoly(u, c, ctx);

            fq_nmod_mpoly_factor_mul_mpoly_ui(g, c, 1, ctx);

            fmpz_add(var_powers + v, var_powers + v, u->exps + u->length - 1);
            _fq_nmod_mpoly_univar_shift_right(u, u->exps + u->length - 1, ctx);

            if (u->length > 1)
            {
                fq_nmod_mpoly_from_univar_bits(c, A->bits, u, v, ctx);
                fq_nmod_mpoly_factor_append_ui(g, c, 1, ctx);
            }
            else
            {
                FLINT_ASSERT(fq_nmod_mpoly_is_one(u->coeffs + 0, ctx));
            }
        }

        fq_nmod_mpoly_factor_swap(f, g, ctx);
    }
/*
	FLINT_ASSERT(fq_nmod_mpoly_factor_is_pairwise_prime(f, ctx) != 0);
*/
    fq_nmod_swap(g->constant, f->constant, ctx->fqctx);
    g->num = 0;
    for (j = 0; j < f->num; j++)
    {
		success = _squarefree_factors(t, f->poly + j, ctx);
		if (!success)
			goto cleanup;

        FLINT_ASSERT(fmpz_is_one(f->exp + j));
        FLINT_ASSERT(fq_nmod_is_one(t->constant, ctx->fqctx));

        fq_nmod_mpoly_factor_fit_length(g, g->num + t->num, ctx);
        for (i = 0; i < t->num; i++)
        {
            fmpz_swap(g->exp + g->num, t->exp + i);
            fq_nmod_mpoly_swap(g->poly + g->num, t->poly + i, ctx);
            g->num++;
        }
    }

    fq_nmod_mpoly_factor_swap(f, g, ctx);

    for (v = 0; v < ctx->minfo->nvars; v++)
    {
        if (fmpz_is_zero(var_powers + v))
            continue;

        fq_nmod_mpoly_factor_fit_length(f, f->num + 1, ctx);
        fq_nmod_mpoly_gen(f->poly + f->num, v, ctx);
        fmpz_swap(f->exp + f->num, var_powers + v);
        f->num++;
    }

    success = 1;

cleanup:

    fq_nmod_mpoly_factor_clear(t, ctx);
    fq_nmod_mpoly_factor_clear(g, ctx);
    fq_nmod_mpoly_univar_clear(u, ctx);
    fq_nmod_mpoly_clear(c, ctx);
    _fmpz_vec_clear(var_powers, ctx->minfo->nvars);

    FLINT_ASSERT(!success || fq_nmod_mpoly_factor_matches(A, f, ctx));

    return success;
}

