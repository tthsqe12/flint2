/*
    Copyright (C) 2019 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "nmod_mpoly_factor.h"

/*
    currently this function doesn't work if the coefficients depend on the main variable
    the assertion x->next == NULL would need to be removed and a loop put in place
    other asserts would need to be removed as well
*/
void nmod_mpoly_from_univar_bits(
	nmod_mpoly_t poly1,
	flint_bitcnt_t bits1,
	const nmod_mpoly_univar_t poly2,
	const nmod_mpoly_ctx_t ctx)
{
    slong i, bits, N;
    ulong k;
    slong next_loc, heap_len = 1;
    ulong * cmpmask;
    slong total_len, p_len;
    mp_limb_t * p_coeff;
    ulong * p_exp;
    slong p_alloc;
    slong var = poly2->var;
    mpoly_heap_s * heap;
    ulong ** poly2_exps;
    ulong * exp;
    ulong * one;
    mpoly_heap_t * chain, * x;
    TMP_INIT;

    bits = bits1;
    FLINT_ASSERT(bits <= FLINT_BITS);

    if (poly2->length == 0)
    {
        nmod_mpoly_fit_bits(poly1, bits, ctx);
        poly1->bits = bits;
		poly1->length = 0;
        return;
    }

    TMP_START;

    /* pack everything into bits */
    N = mpoly_words_per_exp_sp(bits, ctx->minfo);
    one = (ulong*) TMP_ALLOC(N*sizeof(ulong));
    cmpmask = (ulong*) TMP_ALLOC(N*sizeof(ulong));
    mpoly_gen_monomial_sp(one, var, bits, ctx->minfo);
    mpoly_get_cmpmask(cmpmask, N, bits, ctx->minfo);

    poly2_exps = (ulong **) TMP_ALLOC(poly2->length*sizeof(ulong*));
    total_len = 0;
    for (i = 0; i < poly2->length; i++)
    {
        total_len += (poly2->coeffs + i)->length;
        poly2_exps[i] = (poly2->coeffs + i)->exps;
        if (bits != (poly2->coeffs + i)->bits)
        {
            poly2_exps[i] = (ulong *) flint_malloc(
                                  N*(poly2->coeffs + i)->length*sizeof(ulong));
            if (!mpoly_repack_monomials(poly2_exps[i], bits,
                    (poly2->coeffs + i)->exps, (poly2->coeffs + i)->bits,
                                      (poly2->coeffs + i)->length, ctx->minfo))
            {
                FLINT_ASSERT(0 && "repack does not fit");
            }
        }
    }

    nmod_mpoly_fit_length(poly1, total_len, ctx);
    nmod_mpoly_fit_bits(poly1, bits, ctx);
    poly1->bits = bits;

    p_coeff = poly1->coeffs;
    p_exp = poly1->exps;
    p_alloc = poly1->alloc;

    next_loc = poly2->length + 2;
    heap = (mpoly_heap_s *) TMP_ALLOC((poly2->length + 1)*sizeof(mpoly_heap_s));
    exp = (ulong *) TMP_ALLOC(poly2->length*N*sizeof(ulong));
    chain = (mpoly_heap_t *) TMP_ALLOC(poly2->length*sizeof(mpoly_heap_t));

    for (i = 0; i < poly2->length; i++)
    {
        k = poly2->exps[i];
        x = chain + i;
        x->i = i;
        x->j = 0;
        x->next = NULL;
        mpoly_monomial_madd(exp + N*x->i, poly2_exps[x->i] + N*x->j, k, one, N);
        _mpoly_heap_insert(heap, exp + N*i, x, &next_loc, &heap_len, N,
                                                               cmpmask);
    }

    p_len = 0;
    while (heap_len > 1)
    {
        _nmod_mpoly_fit_length(&p_coeff, &p_exp, &p_alloc, p_len + 1, N);
        mpoly_monomial_set(p_exp + N*p_len, heap[1].exp, N);
        x = _mpoly_heap_pop(heap, &heap_len, N, cmpmask);
        p_coeff[p_len] = (poly2->coeffs + x->i)->coeffs[x->j];
        p_len++;

        FLINT_ASSERT(x->next == NULL);

        if (x->j + 1 < (poly2->coeffs + x->i)->length)
        {
            k = poly2->exps[x->i];
            x->j = x->j + 1;
            x->next = NULL;
            mpoly_monomial_madd(exp + N*x->i, poly2_exps[x->i] + N*x->j, k, one, N);
            _mpoly_heap_insert(heap, exp + N*x->i, x, &next_loc, &heap_len, N,
                                                               cmpmask);
        }
    }

    FLINT_ASSERT(total_len == p_len);
    poly1->coeffs = p_coeff;
    poly1->exps = p_exp;
    poly1->alloc = p_alloc;
    poly1->length = p_len;

    for (i = 0; i < poly2->length; i++)
    {
        if (poly2_exps[i] != (poly2->coeffs + i)->exps)
            flint_free(poly2_exps[i]);
    }

    TMP_END;
}




int nmod_mpoly_univar_content_mpoly(
    nmod_mpoly_t g,
    const nmod_mpoly_univar_t A,
    const nmod_mpoly_ctx_t ctx)
{
    slong i;

    nmod_mpoly_zero(g, ctx);
    for (i = 0; i < A->length; i++)
    {
		if (!nmod_mpoly_gcd(g, g, A->coeffs + i, ctx))
			return 0;
    }

    return 1;
}

void nmod_mpoly_univar_divexact_mpoly(
    nmod_mpoly_univar_t A,
    const nmod_mpoly_t b,
    const nmod_mpoly_ctx_t ctx)
{
    int success;
    slong i;

    for (i = 0; i < A->length; i++)
    {
        success = nmod_mpoly_divides(A->coeffs + i, A->coeffs + i, b, ctx);
        FLINT_ASSERT(success);
    }
}

void nmod_mpoly_univar_shift_right(
    nmod_mpoly_univar_t A,
    ulong shift,
    const nmod_mpoly_ctx_t ctx)
{
    slong i;
    for (i = 0; i < A->length; i++)
    {
        FLINT_ASSERT(A->exps[i] >= shift);
        A->exps[i] -= shift;
    }
}











void nmod_mpoly_factor_expand(nmod_mpoly_t a, const nmod_mpoly_factor_t f, const nmod_mpoly_ctx_t ctx)
{
    slong i;
    nmod_mpoly_t t1, t2;

    nmod_mpoly_init(t1, ctx);
    nmod_mpoly_init(t2, ctx);

    nmod_mpoly_set_ui(a, f->content, ctx);

    for (i = 0; i < f->length; i++)
    {
        nmod_mpoly_pow_fmpz(t1, f->poly + i, f->exp + i, ctx);
        nmod_mpoly_mul(t2, a, t1, ctx);
        nmod_mpoly_swap(a, t2, ctx);
    }

    nmod_mpoly_clear(t1, ctx);
    nmod_mpoly_clear(t2, ctx);
}

int nmod_mpoly_factor_matches(const nmod_mpoly_t a, const nmod_mpoly_factor_t f, const nmod_mpoly_ctx_t ctx)
{
    int matches;
    nmod_mpoly_t t;
    nmod_mpoly_init(t, ctx);
    nmod_mpoly_factor_expand(t, f, ctx);
    matches = nmod_mpoly_equal(t, a, ctx);
    nmod_mpoly_clear(t, ctx);
    return matches;
}



/* 	return 0 no, 1 yes, -1 don't know	*/
int nmod_mpoly_factor_is_pairwise_prime(
	const nmod_mpoly_factor_t f, 
	const nmod_mpoly_ctx_t ctx)
{
	int result;
	slong i, j;
	nmod_mpoly_t g;

	nmod_mpoly_init(g, ctx);

	for (i = 0; i + 1 < f->length; i++)
	for (j = i + 1; j < f->length; j++)
	{
		/* make sure factors are monic */
		if (   (f->poly + i)->length == 0
			|| (f->poly + j)->length == 0
			|| (f->poly + i)->coeffs[0] != 1
			|| (f->poly + j)->coeffs[0] != 1)
		{
			result = 0;
			goto cleanup;
		}

		if (nmod_mpoly_gcd(g, f->poly + i, f->poly + j, ctx))
		{
			if (!nmod_mpoly_is_one(g, ctx))
			{
				result = 0;
				goto cleanup;
			}
		}
		else
		{
			result = -1;
			goto cleanup;
		}
	}

	result = 1;

cleanup:

	nmod_mpoly_clear(g, ctx);

	return result;
}

void nmod_mpoly_factor_one(
	nmod_mpoly_factor_t a,
	const nmod_mpoly_ctx_t ctx)
{
	a->content = 1;
	a->length = 0;
}

void nmod_mpoly_factor_pow_fmpz(
	nmod_mpoly_factor_t a,
	const nmod_mpoly_factor_t b,
	const fmpz_t e,
	const nmod_mpoly_ctx_t ctx)
{
	if (a != b)
		nmod_mpoly_factor_set(a, b, ctx);
	a->content = nmod_pow_fmpz(a->content, e, ctx->ffinfo->mod);
	_fmpz_vec_scalar_mul_fmpz(a->exp, a->exp, a->length, e);
}

void nmod_mpoly_factor_push_fmpz(
	nmod_mpoly_factor_t fac,
	const nmod_mpoly_t a,
	const fmpz_t pow,
	const nmod_mpoly_ctx_t ctx)
{
	if (nmod_mpoly_is_ui(a, ctx))
	{
		ulong t = nmod_mpoly_get_ui(a, ctx);
		t = nmod_pow_fmpz(t, pow, ctx->ffinfo->mod);
		fac->content = nmod_mul(fac->content, t, ctx->ffinfo->mod);
		return;
	}
	else
	{
		nmod_mpoly_factor_append_fmpz(fac, a, pow, ctx);
	}
}

void nmod_mpoly_factor_push_ui(
	nmod_mpoly_factor_t fac,
	const nmod_mpoly_t a,
	ulong e,
	const nmod_mpoly_ctx_t ctx)
{
	if (nmod_mpoly_is_ui(a, ctx))
	{
		ulong t = nmod_mpoly_get_ui(a, ctx);
		t = nmod_pow_ui(t, e, ctx->ffinfo->mod);
		fac->content = nmod_mul(fac->content, t, ctx->ffinfo->mod);
		return;
	}
	else
	{
		nmod_mpoly_factor_append_ui(fac, a, e, ctx);
	}
}

/*
	b and c are pairwise prime
	produce a=b*c pairwise prime and return 1
		else return 0 with a undefined
*/
int nmod_mpoly_factor_mul_pairwise_prime(
	nmod_mpoly_factor_t a,
	const nmod_mpoly_factor_t b,
	const nmod_mpoly_factor_t c,
	const nmod_mpoly_ctx_t ctx)
{
	int success;
	slong i, j;
	nmod_mpoly_t T1, T2;
	nmod_mpoly_struct * g;
	fmpz_t t;

	if (a == b || a == c)
	{
		nmod_mpoly_factor_t ta;
		nmod_mpoly_factor_init(ta, ctx);
		success = nmod_mpoly_factor_mul_pairwise_prime(ta, b, c, ctx);
		nmod_mpoly_factor_swap(a, ta, ctx);
		nmod_mpoly_factor_clear(ta, ctx);
		return success;
	}

	fmpz_init(t);
	nmod_mpoly_init(T1, ctx);
	nmod_mpoly_init(T2, ctx);

	FLINT_ASSERT(nmod_mpoly_factor_is_pairwise_prime(b, ctx) == 1);
	FLINT_ASSERT(nmod_mpoly_factor_is_pairwise_prime(c, ctx) == 1);

	g = (nmod_mpoly_struct *) flint_malloc(b->length*c->length*sizeof(nmod_mpoly_struct));
	/* g[i,j] = gcd(b[i], c[j]) */
	for (i = 0; i < b->length; i++)
	for (j = 0; j < c->length; j++)
	{
		nmod_mpoly_init(g + i*c->length + j, ctx);
	}

	a->content = nmod_mul(b->content, c->content, ctx->ffinfo->mod);
	a->length = 0;

	for (i = 0; i < b->length; i++)
	for (j = 0; j < c->length; j++)
	{
		if (!nmod_mpoly_gcd(g + i*c->length + j, b->poly + i, c->poly + j, ctx))
		{
			success = 0;
			goto cleanup;
		}

		fmpz_add(t, b->exp + i, c->exp + j);
		nmod_mpoly_factor_push_fmpz(a, g + i*c->length + j, t, ctx);
	}

	for (i = 0; i < b->length; i++)
	{
		nmod_mpoly_set(T1, b->poly + i, ctx);
		for (j = 0; j < c->length; j++)
		{
			success = nmod_mpoly_divides(T1, T1, g + i*c->length + j, ctx); FLINT_ASSERT(success);
		}
		nmod_mpoly_factor_push_fmpz(a, T1, b->exp + i, ctx);
	}	

	for (j = 0; j < c->length; j++)
	{
		nmod_mpoly_set(T1, c->poly + j, ctx);
		for (i = 0; i < b->length; i++)
		{
			success = nmod_mpoly_divides(T1, T1, g + i*c->length + j, ctx); FLINT_ASSERT(success);
		}
		nmod_mpoly_factor_push_fmpz(a, T1, c->exp + j, ctx);
	}	

	success = 1;

cleanup:

	/* g[i,j] = gcd(b[i], c[j]) */
	for (i = 0; i < b->length; i++)
	for (j = 0; j < c->length; j++)
	{
		nmod_mpoly_clear(g + i*c->length + j, ctx);
	}
	flint_free(g);

	nmod_mpoly_clear(T1, ctx);
	nmod_mpoly_clear(T2, ctx);
	fmpz_clear(t);

	if (success)
	{
		nmod_mpoly_t ae, be, ce;
		nmod_mpoly_init(ae, ctx);
		nmod_mpoly_init(be, ctx);
		nmod_mpoly_init(ce, ctx);
		FLINT_ASSERT(nmod_mpoly_factor_is_pairwise_prime(a, ctx) == 1);
		nmod_mpoly_factor_expand(be, b, ctx);
		nmod_mpoly_factor_expand(ce, c, ctx);
		nmod_mpoly_mul(ae, be, ce, ctx);
		FLINT_ASSERT(nmod_mpoly_factor_matches(ae, a, ctx));
		nmod_mpoly_clear(ae, ctx);
		nmod_mpoly_clear(be, ctx);
		nmod_mpoly_clear(ce, ctx);
	}
	return success;
}


/* f *= g^e, no combination */
void nmod_mpoly_factor_append_factor_fmpz(
	nmod_mpoly_factor_t f,
	nmod_mpoly_factor_t g,
	fmpz_t e,
	const nmod_mpoly_ctx_t ctx)
{
	slong i;
	fmpz_t c;

	fmpz_init(c);

	f->content = nmod_mul(f->content, g->content, ctx->ffinfo->mod);
	for (i = 0; i < g->length; i++)
	{
		fmpz_mul(c, g->exp + i, e);
		nmod_mpoly_factor_append_fmpz(f, g->poly + i, c, ctx);
	}

	fmpz_clear(c);
}

void nmod_mpoly_factor_append_factor_ui(
	nmod_mpoly_factor_t f,
	nmod_mpoly_factor_t g,
	ulong e,
	const nmod_mpoly_ctx_t ctx)
{
	slong i;
	fmpz_t c;

	fmpz_init(c);

	f->content = nmod_mul(f->content, nmod_pow_ui(g->content, e, ctx->ffinfo->mod), ctx->ffinfo->mod);
	for (i = 0; i < g->length; i++)
	{
		fmpz_mul_si(c, g->exp + i, e);
		nmod_mpoly_factor_append_fmpz(f, g->poly + i, c, ctx);
	}

	fmpz_clear(c);
}


int _factor_sqrfd(
	nmod_mpoly_t c,
	nmod_mpoly_factor_t ans,
	slong var,
	const nmod_mpoly_ctx_t ctx)
{
	ulong i;
	int success;
	nmod_mpoly_t a, b, u, v, w, g;
/*
flint_printf("_factor_sqrfd(var = %wd) called c: ", var); nmod_mpoly_print_pretty(c, NULL, ctx); printf("\n");
*/
	nmod_mpoly_factor_one(ans, ctx);
	nmod_mpoly_init(a, ctx);
	nmod_mpoly_init(b, ctx);
	nmod_mpoly_init(u, ctx);
	nmod_mpoly_init(v, ctx);
	nmod_mpoly_init(w, ctx);
	nmod_mpoly_init(g, ctx);

	nmod_mpoly_set(a, c, ctx);
	nmod_mpoly_derivative(b, a, var, ctx);
	success = nmod_mpoly_gcd(c, a, b, ctx);
	if (!success)
		goto cleanup;

	success = nmod_mpoly_divides(w, a, c, ctx); FLINT_ASSERT(success);
	success = nmod_mpoly_divides(v, b, c, ctx); FLINT_ASSERT(success);
	nmod_mpoly_derivative(g, w, var, ctx);
	nmod_mpoly_sub(u, v, g, ctx);

	for (i = 1; i + 1 < ctx->ffinfo->mod.n && !nmod_mpoly_is_zero(u, ctx); i++)
	{
		success = nmod_mpoly_gcd(g, w, u, ctx);
		if (!success)
			goto cleanup;

		nmod_mpoly_factor_push_ui(ans, g, i, ctx);
		success = nmod_mpoly_divides(w, w, g, ctx); FLINT_ASSERT(success);
		success = nmod_mpoly_divides(c, c, w, ctx); FLINT_ASSERT(success);
		success = nmod_mpoly_divides(v, u, g, ctx); FLINT_ASSERT(success);
		nmod_mpoly_derivative(g, w, var, ctx);
		nmod_mpoly_sub(u, v, g, ctx);
	}

	nmod_mpoly_factor_push_ui(ans, w, i, ctx);

	success = 1;

cleanup:

	nmod_mpoly_clear(a, ctx);
	nmod_mpoly_clear(b, ctx);
	nmod_mpoly_clear(u, ctx);
	nmod_mpoly_clear(v, ctx);
	nmod_mpoly_clear(w, ctx);
	nmod_mpoly_clear(g, ctx);

	return success;
}

int _factor_sqrf(
	nmod_mpoly_factor_t r,
	const nmod_mpoly_t a,
	const nmod_mpoly_ctx_t ctx)
{
	int success;
	slong v;
	fmpz * shift, * stride;
	nmod_mpoly_t c;
	nmod_mpoly_factor_t ans, tfac;
	fmpz_t g, gr, p, pk;
/*
printf("_factor_sqrf called a: "); nmod_mpoly_print_pretty(a, NULL, ctx); printf("\n");
*/
	fmpz_init_set_ui(p, ctx->ffinfo->mod.n);
	fmpz_init(pk);
	fmpz_init(g);
	fmpz_init(gr);
	nmod_mpoly_init(c, ctx);
	nmod_mpoly_factor_init(ans, ctx);
	nmod_mpoly_factor_init(tfac, ctx);
	shift = _fmpz_vec_init(ctx->minfo->nvars);
	stride = _fmpz_vec_init(ctx->minfo->nvars);

	nmod_mpoly_factor_one(r, ctx);
	nmod_mpoly_set(c, a, ctx);

	for (v = 0; v < ctx->minfo->nvars; v++)
	{
		success = _factor_sqrfd(c, ans, v, ctx);
		if (!success)
			goto cleanup;

		nmod_mpoly_factor_append_factor_ui(r, ans, 1, ctx);
	}

	if (nmod_mpoly_is_ui(c, ctx))
	{
		r->content = nmod_mul(r->content, nmod_mpoly_get_ui(c, ctx), ctx->ffinfo->mod);
	}
	else
	{
		nmod_mpoly_deflation(shift, stride, c, ctx);
		fmpz_zero(g);
		for (v = 0; v < ctx->minfo->nvars; v++)
		{
			fmpz_gcd(g, g, stride + v);
			fmpz_gcd(g, g, shift + v);
		}

		FLINT_ASSERT(fmpz_remove(gr, g, p) > 0);
		fmpz_pow_ui(pk, p, fmpz_remove(gr, g, p));
		for (v = 0; v < ctx->minfo->nvars; v++)
		{
			fmpz_set(stride + v, pk);
			fmpz_zero(shift + v);
		}

		nmod_mpoly_deflate(c, c, shift, stride, ctx);

		success = _factor_sqrf(tfac, c, ctx);
		if (!success)
			goto cleanup;

		nmod_mpoly_factor_pow_fmpz(tfac, tfac, pk, ctx);
		nmod_mpoly_factor_mul_pairwise_prime(r, r, tfac, ctx);
	}

	success = 1;

cleanup:

	fmpz_clear(p);
	fmpz_clear(pk);
	fmpz_clear(g);
	fmpz_clear(gr);
	nmod_mpoly_clear(c, ctx);
	nmod_mpoly_factor_clear(ans, ctx);
	nmod_mpoly_factor_clear(tfac, ctx);
	_fmpz_vec_clear(shift, ctx->minfo->nvars);
	_fmpz_vec_clear(stride, ctx->minfo->nvars);

	return success;
}



int nmod_mpoly_factor(nmod_mpoly_factor_t fac, const nmod_mpoly_t A, int full, const nmod_mpoly_ctx_t ctx)
{
    int success;
    slong j, v;
    nmod_mpoly_factor_t newfac, sqrfac;
    fmpz_t cc;
    nmod_mpoly_t c;
    nmod_mpoly_univar_t u;
/*
flint_printf("nmod_mpoly_factor called\nA: "); nmod_mpoly_print_pretty(A, NULL, ctx); printf("\n");
*/
    /* 0. set trivial factorization */
    fac->length = 0;
    if (nmod_mpoly_is_ui(A, ctx))
    {
		fac->content = nmod_mpoly_get_ui(A, ctx);
        return 1;
    }
    else
    {
		fac->content = A->coeffs[0];
		nmod_mpoly_factor_fit_length(fac, 1, ctx);
		nmod_mpoly_make_monic(fac->poly + 0, A, ctx);
		fmpz_one(fac->exp + 0);
		fac->length = 1;
    }

    if (A->bits > FLINT_BITS)
    {
        return 0;
    }

    nmod_mpoly_factor_init(newfac, ctx);
    nmod_mpoly_factor_init(sqrfac, ctx);
    nmod_mpoly_univar_init(u, ctx);
    nmod_mpoly_init(c, ctx);
    fmpz_init(cc);

    /* 1. ensure factors are primitive w.r.t any variable */
/*
printf("******** step 1 fac: "); nmod_mpoly_factor_print_pretty(fac, NULL, ctx); printf("\n");
*/
    FLINT_ASSERT(nmod_mpoly_factor_matches(A, fac, ctx));

    for (v = 0; v < ctx->minfo->nvars; v++)
    {
        newfac->content = fac->content;
        newfac->length = 0;

        for (j = 0; j < fac->length; j++)
        {
            if (fmpz_is_one(fac->exp + j))
            {
                nmod_mpoly_to_univar(u, fac->poly + j, v, ctx);
                FLINT_ASSERT(u->length > 0);

                nmod_mpoly_univar_content_mpoly(c, u, ctx);
                nmod_mpoly_univar_divexact_mpoly(u, c, ctx);

                if (nmod_mpoly_is_ui(c, ctx))
                {
                    FLINT_ASSERT(nmod_mpoly_is_one(c, ctx));
                }
                else
                {
                    nmod_mpoly_factor_append_ui(newfac, c, 1, ctx);
                }

                FLINT_ASSERT(u->length > 0);

                if (u->exps[u->length - 1] != 0)
                {
                    nmod_mpoly_gen(c, v, ctx);
                    nmod_mpoly_factor_append_ui(newfac, c, u->exps[u->length - 1], ctx);
                    nmod_mpoly_univar_shift_right(u, u->exps[u->length - 1], ctx);
                }

                if (u->length > 1)
                {
                    nmod_mpoly_from_univar_bits(c, A->bits, u, ctx);
                    nmod_mpoly_factor_append_ui(newfac, c, 1, ctx);
                }
                else
                {
                    FLINT_ASSERT(nmod_mpoly_is_one(u->coeffs + 0, ctx));
                }
            }
            else
            {
                FLINT_ASSERT(fmpz_cmp_ui(fac->exp + j, 1) > 0);
                nmod_mpoly_factor_append_fmpz(newfac, fac->poly + j, fac->exp + j, ctx);
            }
        }

        nmod_mpoly_factor_swap(fac, newfac, ctx);
    }

    /* 2. ensure factors are squarefree */
/*
printf("******** step 2 fac: "); nmod_mpoly_factor_print_pretty(fac, NULL, ctx); printf("\n");
*/
    FLINT_ASSERT(nmod_mpoly_factor_matches(A, fac, ctx));
	FLINT_ASSERT(nmod_mpoly_factor_is_pairwise_prime(fac, ctx) != 0);

    newfac->content = fac->content;
    newfac->length = 0;
    for (j = 0; j < fac->length; j++)
    {
		success = _factor_sqrf(sqrfac, fac->poly + j, ctx);
/*
printf("_factor_sqrf returned %d: ", success); nmod_mpoly_factor_print_pretty(fac, NULL, ctx); printf("\n");
*/
		if (!success)
			goto cleanup;
		nmod_mpoly_factor_append_factor_fmpz(newfac, sqrfac, fac->exp + j, ctx);
    }
    nmod_mpoly_factor_swap(fac, newfac, ctx);

    success = 1;

cleanup:

    FLINT_ASSERT(nmod_mpoly_factor_matches(A, fac, ctx));

    nmod_mpoly_factor_clear(sqrfac, ctx);
    nmod_mpoly_factor_clear(newfac, ctx);
    nmod_mpoly_univar_clear(u, ctx);
    nmod_mpoly_clear(c, ctx);
    fmpz_clear(cc);

    return 1;
}
