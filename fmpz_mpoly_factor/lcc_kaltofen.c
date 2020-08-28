#include "fmpz_mpoly_factor.h"
#include "n_poly.h"


int fmpz_mpoly_evaluate_rest_except_one(
    fmpz_poly_t e,
    const fmpz_mpoly_t A,
    const fmpz * alphas,
    slong v,
    const fmpz_mpoly_ctx_t ctx)
{
    int success;
    slong i;
    fmpz_mpoly_t t;

    fmpz_mpoly_init(t, ctx);
    fmpz_mpoly_set(t, A, ctx);

    for (i = 1; i < ctx->minfo->nvars; i++)
    {
        if (i == v)
            continue;

        if (!fmpz_mpoly_evaluate_one_fmpz(t, t, i, alphas + i - 1, ctx))
        {
            success = 0;
            goto cleanup;
        }
    }

    success = fmpz_mpoly_is_fmpz_poly(t, v, ctx) &&
              fmpz_mpoly_get_fmpz_poly(e, t, v, ctx);

cleanup:

    fmpz_mpoly_clear(t, ctx);

    return success;
}



static void _make_bases_coprime(
    fmpz_poly_factor_t A,
    fmpz_poly_factor_t B)
{
    slong i, j;
    slong Alen = A->num;
    slong Blen = B->num;
    fmpz_poly_t g;

    fmpz_poly_init(g);

    for (i = 0; i < Alen; i++)
    for (j = 0; j < Blen; j++)
    {
        fmpz_poly_gcd(g, A->p + i, B->p + j);
        if (fmpz_poly_degree(g) > 0)
        {
            fmpz_poly_div(A->p + i, A->p + i, g);
            fmpz_poly_div(B->p + j, B->p + j, g);
            fmpz_poly_factor_fit_length(A, A->num + 1);
            fmpz_poly_set(A->p + A->num, g);
            A->exp[A->num] = A->exp[i];
            A->num++;
            fmpz_poly_factor_fit_length(B, B->num + 1);
            fmpz_poly_set(B->p + B->num, g);
            B->exp[B->num] = B->exp[j];
            B->num++;
        }
    }

    for (i = 0; i < A->num; i++)
    {
        if (fmpz_poly_degree(A->p + i) > 0)
            continue;
        A->num--;
        fmpz_poly_swap(A->p + i, A->p + A->num);
        SLONG_SWAP(A->exp[i], A->exp[A->num]);
        i--;
    }

    for (i = 0; i < B->num; i++)
    {
        if (fmpz_poly_degree(B->p + i) > 0)
            continue;
        B->num--;
        fmpz_poly_swap(B->p + i, B->p + B->num);
        SLONG_SWAP(B->exp[i], B->exp[B->num]);
        i--;
    }

    fmpz_poly_clear(g);
}

void fmpz_poly_vector_insert_poly(fmpz_bpoly_t v, const fmpz_poly_t a)
{
    slong i;

    for (i = 0; i < v->length; i++)
        if (fmpz_poly_equal(v->coeffs + i, a))
            return;

    fmpz_bpoly_fit_length(v, v->length + 1);
    fmpz_poly_set(v->coeffs + v->length, a);
    v->length++;
}


void fmpz_poly_factor_print_pretty(fmpz_poly_factor_t f, const char * x)
{
    slong i;
    fmpz_print(&f->c);
    for (i = 0; i < f->num; i++)
    {
        flint_printf("*(");
        fmpz_poly_print_pretty(f->p + i, x);
        flint_printf(")^%wd", f->exp[i]);
    }
}

static int _try_lift(
    fmpz_mpoly_struct * lifts,  /* length r */
    const fmpz_mpoly_t A,
    fmpz_poly_struct * Auf,     /* length r */
    slong r,
    const fmpz * alphas,
    slong v,
    const fmpz_mpoly_ctx_t Actx)
{
    int success;
    slong i, k;
    slong mvars, nvars = Actx->minfo->nvars;
    slong * Adegs, * Bdegs, * perm, * iperm;
    fmpz * Balphas;
    slong dummyvars[1];
    ulong dummydegs[1];
    fmpz_mpoly_t lcA, At, newA;
    fmpz_mpoly_ctx_t Bctx;
    flint_bitcnt_t Bbits;
    fmpz_mpoly_t Bt;
    fmpz_mpolyv_t fac, tfac;
    fmpz_mpoly_struct * Bevals, * Blcevals;
    fmpz_mpoly_univar_t u;
    fmpz_t q;
/*
flint_printf("try_lift called\n");
*/
    FLINT_ASSERT(0 < v && v < nvars);

    if (r < 2)
    {
        FLINT_ASSERT(r == 1);
        fmpz_mpoly_init(At, Actx);
        fmpz_mpoly_univar_init(u, Actx);
        fmpz_mpoly_to_univar(u, A, v, Actx);
        success = _fmpz_mpoly_vec_content_mpoly(At, u->coeffs, u->length, Actx);
        if (success)
        {
            success = fmpz_mpoly_divides(lifts + 0, A, At, Actx);
            FLINT_ASSERT(success);
            fmpz_mpoly_unit_normalize(lifts + 0, Actx);
        }
        fmpz_mpoly_clear(At, Actx);
        fmpz_mpoly_univar_clear(u, Actx);
        return success;
    }

    Adegs = FLINT_ARRAY_ALLOC(nvars, slong);
    Bdegs = FLINT_ARRAY_ALLOC(nvars, slong);
    perm = FLINT_ARRAY_ALLOC(nvars, slong);
    iperm = FLINT_ARRAY_ALLOC(nvars, slong);
    fmpz_init(q);
    fmpz_mpoly_init(lcA, Actx);
    fmpz_mpoly_init(At, Actx);
    fmpz_mpoly_init(newA, Actx);

    dummyvars[0] = v;
    dummydegs[0] = fmpz_mpoly_degree_si(A, v, Actx);
    fmpz_mpoly_get_coeff_vars_ui(lcA, A, dummyvars, dummydegs, 1, Actx);
    fmpz_mpoly_pow_ui(At, lcA, r - 1, Actx);
    fmpz_mpoly_mul(newA, A, At, Actx);

    if (newA->bits >= FLINT_BITS)
    {
        success = 0;
        goto cleanup_less;
    }

    fmpz_mpoly_degrees_si(Adegs, newA, Actx);

    perm[0] = v;
    mvars = 1;
    for (i = 0; i < nvars; i++)
    {
        if (i == v)
            continue;
        iperm[i] = -1;
        if (Adegs[i] > 0)
        {
            FLINT_ASSERT(i > 0);
            perm[mvars] = i;
            mvars++;
        }
    }

    /* TODO following code works if mvars = 1, but add special case? */

    fmpz_mpoly_ctx_init(Bctx, mvars, ORD_LEX);
    fmpz_mpoly_init(Bt, Bctx);
    fmpz_mpolyv_init(fac, Bctx);
    fmpz_mpolyv_init(tfac, Bctx);
    fmpz_mpoly_univar_init(u, Bctx);
    Balphas = _fmpz_vec_init(nvars);
    Bevals = FLINT_ARRAY_ALLOC(mvars, fmpz_mpoly_struct);
    Blcevals = FLINT_ARRAY_ALLOC(mvars, fmpz_mpoly_struct);
    for (i = 0; i < mvars; i++)
    {
        fmpz_mpoly_init(Bevals + i, Bctx);
        fmpz_mpoly_init(Blcevals + i, Bctx);
    }

    Bbits = mpoly_fix_bits(newA->bits, Bctx->minfo);

    /* invert perm */
    for (i = 0; i < mvars; i++)
    {
        iperm[perm[i]] = i;
        Bdegs[i] = Adegs[perm[i]];
        if (i > 0)
            fmpz_set(Balphas + i - 1, alphas + perm[i] - 1);
    }

    fmpz_mpoly_convert_perm(Bevals + mvars - 1, Bbits, Bctx, newA, Actx, perm);
    fmpz_mpoly_convert_perm(Blcevals + mvars - 1, Bbits, Bctx, lcA, Actx, perm);
/*
flint_printf("Bevals[%wd]: ", mvars - 1);
fmpz_mpoly_print_pretty(Bevals + mvars - 1, NULL, Bctx);
flint_printf("\n");

flint_printf("Blcevals[%wd]: ", mvars - 1);
fmpz_mpoly_print_pretty(Blcevals + mvars - 1, NULL, Bctx);
flint_printf("\n");
*/

    for (i = mvars - 2; i >= 0; i--)
    {
/*
flint_printf("Bevals[%d]: ", i); fmpz_print(Balphas + i); flint_printf("\n");
*/
        fmpz_mpoly_evaluate_one_fmpz(Bevals + i, Bevals + i + 1, i + 1, Balphas + i, Bctx);
        fmpz_mpoly_evaluate_one_fmpz(Blcevals + i, Blcevals + i + 1, i + 1, Balphas + i, Bctx);
/*
flint_printf("Bevals[%wd]: ", i);
fmpz_mpoly_print_pretty(Bevals + i, NULL, Bctx);
flint_printf("\n");

flint_printf("Blcevals[%wd]: ", i);
fmpz_mpoly_print_pretty(Blcevals + i, NULL, Bctx);
flint_printf("\n");
*/
    }

    fmpz_mpolyv_fit_length(fac, r, Bctx);
    fac->length = r;
    for (i = 0; i < r; i++)
    {
/*
flint_printf("Auf[%wd]: ", i);
fmpz_poly_print_pretty(Auf + i, "t");
flint_printf("\n");
*/
        FLINT_ASSERT(fmpz_mpoly_is_fmpz(Blcevals + 0, Bctx));
        FLINT_ASSERT(fmpz_mpoly_length(Blcevals + 0, Bctx) == 1);
        FLINT_ASSERT(fmpz_divisible(Blcevals[0].coeffs + 0, Auf[i].coeffs + Auf[i].length - 1));
        fmpz_divexact(q, Blcevals[0].coeffs + 0, Auf[i].coeffs + Auf[i].length - 1);
        _fmpz_mpoly_set_fmpz_poly(fac->coeffs + i, Bbits, Auf[i].coeffs, Auf[i].length, 0, Bctx);
        fmpz_mpoly_scalar_mul_fmpz(fac->coeffs + i, fac->coeffs + i, q, Bctx);
    }

    fmpz_mpolyv_fit_length(tfac, r, Bctx);
    tfac->length = r;
    for (k = 1; k <= mvars - 1; k++)
    {
        for (i = 0; i < r; i++)
            _fmpz_mpoly_set_lead0(tfac->coeffs + i, fac->coeffs + i, Blcevals + k, Bctx);

        success = fmpz_mpoly_hlift(k, tfac->coeffs, r, Balphas, Bevals + k, Bdegs, Bctx);
        if (!success)
            goto cleanup_more;

        fmpz_mpolyv_swap(tfac, fac, Bctx);
    }

    for (i = 0; i < r; i++)
    {
        fmpz_mpoly_to_univar(u, fac->coeffs + i, 0, Bctx);
        success = _fmpz_mpoly_vec_content_mpoly(Bt, u->coeffs, u->length, Bctx);
        if (!success)
            goto cleanup_more;
        success = fmpz_mpoly_divides(Bt, fac->coeffs + i, Bt, Bctx);
        FLINT_ASSERT(success);
        fmpz_mpoly_convert_perm(lifts + i, A->bits, Actx, Bt, Bctx, iperm);
        fmpz_mpoly_unit_normalize(lifts + i, Actx);
/*
flint_printf("lifts[%wd]: ", i);
fmpz_mpoly_print_pretty(lifts + i, NULL, Actx);
flint_printf("\n");
*/
    }

    success = 1;

cleanup_more:

    fmpz_mpoly_clear(Bt, Bctx);
    fmpz_mpolyv_clear(fac, Bctx);
    fmpz_mpolyv_clear(tfac, Bctx);
    fmpz_mpoly_univar_clear(u, Bctx);
    _fmpz_vec_clear(Balphas, nvars);
    for (i = 0; i < mvars; i++)
    {
        fmpz_mpoly_clear(Bevals + i, Bctx);
        fmpz_mpoly_clear(Blcevals + i, Bctx);
    }
    flint_free(Bevals);
    flint_free(Blcevals);

    fmpz_mpoly_ctx_clear(Bctx);

cleanup_less:

    flint_free(Adegs);
    flint_free(Bdegs);
    flint_free(perm);
    flint_free(iperm);

    fmpz_clear(q);
    fmpz_mpoly_clear(lcA, Actx);
    fmpz_mpoly_clear(At, Actx);
    fmpz_mpoly_clear(newA, Actx);

    return success;
}

/* assume content(b) = 1 for now */
void fmpz_mpoly_factor_divexact_mpoly_pow_ui(
    fmpz_mpoly_factor_t A,
    const fmpz_mpoly_t b_in,
    ulong e,
    const fmpz_mpoly_ctx_t ctx)
{
    int sgn;
    slong i;
    fmpz_mpoly_t b_copy;
    fmpz_mpoly_struct * b = (fmpz_mpoly_struct *) b_in;

    fmpz_mpoly_init(b_copy, ctx);

    i = 0; /* index strickly before which everything is coprime to b */

    while (i < A->num && !fmpz_mpoly_is_fmpz(b, ctx))
    {
        fmpz_mpoly_factor_fit_length(A, A->num + 1, ctx);
        fmpz_mpoly_gcd_cofactors(A->poly + A->num, A->poly + i, b_copy,
                                                          A->poly + i, b, ctx);
        b = b_copy;
        if (fmpz_mpoly_is_fmpz(A->poly + A->num, ctx))
        {
            i++;
            continue;
        }

        fmpz_sub_ui(A->exp + A->num, A->exp + i, e);
        sgn = fmpz_sgn(A->exp + A->num);
        if (sgn < 0)
        {
            flint_printf("non-exact division fmpz_mpoly_factor_divexact_mpoly_pow_ui");
            flint_abort();
        }
        else if (sgn > 0)
        {
            A->num++;
        }

        if (fmpz_mpoly_is_fmpz(A->poly + i, ctx))
        {
            A->num--;
            fmpz_mpoly_swap(A->poly + i, A->poly + A->num, ctx);
            fmpz_swap(A->exp + i, A->exp + A->num);
        }
        else
        {
            i++;
        }
    }

    if (!fmpz_mpoly_is_fmpz(b, ctx))
    {
        flint_printf("non-exact division fmpz_mpoly_factor_divexact_mpoly_pow_ui");
        flint_abort();
    }

    fmpz_mpoly_clear(b_copy, ctx);
}



int fmpz_mpoly_factor_lcc_kaltofen_step(
    fmpz_mpoly_struct * divs,   /* length r */
    slong r,
    fmpz_mpoly_factor_t Af, /* squarefree factorization of A */
    const fmpz_poly_struct * Au,
    slong v,                      /* minor bivar var*/
    const fmpz * alphas,
    const fmpz_mpoly_ctx_t ctx)
{
    int success;
    slong i, j, k, Af_deg_v;
    fmpz_poly_factor_struct * Auf;
    fmpz_mpoly_t Afp, tt;
    fmpz_mpolyv_t lfp;
    fmpz_bpoly_t f;
    fmpz_poly_t t, fp;

    fmpz_mpoly_init(Afp, ctx);
    fmpz_mpoly_init(tt, ctx);
    fmpz_bpoly_init(f);
    fmpz_poly_init(fp);
    fmpz_poly_init(t);
    fmpz_mpolyv_init(lfp, ctx);
    Auf = FLINT_ARRAY_ALLOC(r, fmpz_poly_factor_struct);
    for (i = 0; i < r; i++)
    {
        fmpz_poly_factor_init(Auf + i);
        fmpz_poly_factor_squarefree(Auf + i, Au + i);
    }

    Af_deg_v = 0;
    fmpz_mpoly_one(Afp, ctx);
    for (i = 0; i < Af->num; i++)
    {
        slong this_degree = fmpz_mpoly_degree_si(Af->poly + i, v, ctx);
        Af_deg_v += this_degree;
        /* we don't need to introduce needless content into Afp */
        if (this_degree != 0)
            fmpz_mpoly_mul(Afp, Afp, Af->poly + i, ctx);
    }

    if (!fmpz_mpoly_evaluate_rest_except_one(t, Afp, alphas, v, ctx))
    {
        success = 0;
        goto cleanup;
    }
    if (fmpz_poly_degree(t) < 1 || fmpz_poly_degree(t) != Af_deg_v)
    {
        success = 0;
        goto cleanup;
    }

    for (i = 0; i < r - 1; i++)
    for (j = i + 1; j < r; j++)
        _make_bases_coprime(Auf + i, Auf + j);

    f->length = 0;
    for (i = 0; i < r; i++)
    for (j = 0; j < Auf[i].num; j++)
        fmpz_poly_vector_insert_poly(f, Auf[i].p + j);

    fmpz_poly_primitive_part(t, t);

    fmpz_poly_one(fp);
    for (i = 0; i < f->length; i++)
        fmpz_poly_mul(fp, fp, f->coeffs + i);

    success = fmpz_poly_equal(fp, t);
    if (!success)
        goto cleanup;

    fmpz_mpolyv_fit_length(lfp, f->length, ctx);
    lfp->length = f->length;
    success = _try_lift(lfp->coeffs, Afp, f->coeffs, f->length, alphas, v, ctx);
    if (!success)
        goto cleanup;

    for (i = 0; i < r; i++)
    {
        fmpz_poly_factor_struct * auf = Auf + i;
        for (j = 0; j < auf->num; j++)
        {
            for (k = 0; k < f->length; k++)
            {
                if (fmpz_poly_equal(f->coeffs + k, auf->p + j))
                {
                    fmpz_mpoly_factor_divexact_mpoly_pow_ui(Af,
                                            lfp->coeffs + k, auf->exp[j], ctx);
                    fmpz_mpoly_pow_ui(tt, lfp->coeffs + k, auf->exp[j], ctx);
                    fmpz_mpoly_mul(divs + i, divs + i, tt, ctx);
                }
            }
        }
    }

    success = 1;

cleanup:

    fmpz_mpoly_clear(Afp, ctx);
    fmpz_mpoly_clear(tt, ctx);
    fmpz_bpoly_clear(f);
    fmpz_poly_clear(fp);
    fmpz_poly_clear(t);
    fmpz_mpolyv_clear(lfp, ctx);
    for (i = 0; i < r; i++)
        fmpz_poly_factor_clear(Auf + i);
    flint_free(Auf);

    return success;
}


