
#include "nmod_mpoly.h"


void
nmod_polydr_clear(nmod_polydr_t poly, const nmod_ctx_t ctx)
{
    if (poly->coeffs)
        flint_free(poly->coeffs);
}

void
nmod_polydr_realloc(nmod_polydr_t poly, slong alloc, const nmod_ctx_t ctx)
{
    if (alloc == 0)
    {
        nmod_polydr_clear(poly, ctx);
        poly->length = 0;
        poly->alloc = 0;
        poly->coeffs = NULL;

        return;
    }

    poly->coeffs = (mp_ptr) flint_realloc(poly->coeffs, alloc * sizeof(mp_limb_t));

    poly->alloc = alloc;

    /* truncate poly if necessary */
    if (poly->length > alloc)
    {
        poly->length = alloc;
        _nmod_polydr_normalise(poly);
    }
}

void
nmod_polydr_fit_length(nmod_polydr_t poly, slong alloc, const nmod_ctx_t ctx)
{
    if (alloc > poly->alloc)
    {
        if (alloc < 2 * poly->alloc)
            alloc = 2 * poly->alloc;

        nmod_polydr_realloc(poly, alloc, ctx);
    }
}


void nmod_polydr_set_coeff_ui(nmod_polydr_t poly, slong j, ulong c, const nmod_ctx_t ctx)
{
    if (c >= ctx->mod.n)
        NMOD_RED(c, c, ctx->mod);

    nmod_polydr_fit_length(poly, j + 1, ctx);

    if (j + 1 < poly->length) /* interior */
        poly->coeffs[j] = c;
    else if (j + 1 == poly->length) /* leading coeff */
    {
        if (c != 0)
            poly->coeffs[j] = c;
        else
        {
            poly->length--;
            _nmod_polydr_normalise(poly);
        }
    } else /* extend polynomial */
    {
        if (c == 0) return;
        else
        {
            flint_mpn_zero(poly->coeffs + poly->length, j - poly->length);

            poly->coeffs[j] = c;
            poly->length = j + 1;
        }
    }
}

void
nmod_polydr_add(nmod_polydr_t res, const nmod_polydr_t poly1,
                               const nmod_polydr_t poly2, const nmod_ctx_t ctx)
{
    slong max = FLINT_MAX(poly1->length, poly2->length);

    nmod_polydr_fit_length(res, max, ctx);

    _nmod_poly_add(res->coeffs, poly1->coeffs, poly1->length, poly2->coeffs,
                   poly2->length, ctx->mod);

    res->length = max;
    _nmod_polydr_normalise(res);  /* there may have been cancellation */
}


void
nmod_polydr_sub(nmod_polydr_t res, const nmod_polydr_t poly1,
                               const nmod_polydr_t poly2, const nmod_ctx_t ctx)
{
    slong max = FLINT_MAX(poly1->length, poly2->length);

    nmod_polydr_fit_length(res, max, ctx);

    _nmod_poly_sub(res->coeffs, poly1->coeffs, poly1->length, poly2->coeffs,
                   poly2->length, ctx->mod);

    res->length = max;
    _nmod_polydr_normalise(res);  /* there may have been cancellation */
}


void nmod_polydr_mul(nmod_polydr_t res, const nmod_polydr_t poly1,
                               const nmod_polydr_t poly2, const nmod_ctx_t ctx)
{
    slong len1, len2, len_out;
    
    len1 = poly1->length;
    len2 = poly2->length;

    if (len1 == 0 || len2 == 0)
    {
        nmod_polydr_zero(res, ctx);
        return;
    }

    len_out = poly1->length + poly2->length - 1;

    if (res == poly1 || res == poly2)
    {
        nmod_polydr_t temp;

        nmod_polydr_init2(temp, len_out, ctx);

        if (len1 >= len2)
            _nmod_poly_mul(temp->coeffs, poly1->coeffs, len1,
                           poly2->coeffs, len2, ctx->mod);
        else
            _nmod_poly_mul(temp->coeffs, poly2->coeffs, len2,
                           poly1->coeffs, len1, ctx->mod);
        
        nmod_polydr_swap(temp, res, ctx);
        nmod_polydr_clear(temp, ctx);
    }
    else
    {
        nmod_polydr_fit_length(res, len_out, ctx);
        
        if (len1 >= len2)
            _nmod_poly_mul(res->coeffs, poly1->coeffs, len1,
                           poly2->coeffs, len2, ctx->mod);
        else
            _nmod_poly_mul(res->coeffs, poly2->coeffs, len2,
                           poly1->coeffs, len1, ctx->mod);
    }

    res->length = len_out;
    _nmod_polydr_normalise(res);
}


void nmod_polydr_make_monic(nmod_polydr_t output, const nmod_polydr_t input,
                                                          const nmod_ctx_t ctx)
{
    if (input->length == 0)
    {
        flint_printf("Exception (nmod_poly_make_monic). Division by zero.\n");
        flint_abort();
    }

    nmod_polydr_fit_length(output, input->length, ctx);
    _nmod_poly_make_monic(output->coeffs, 
                            input->coeffs, input->length, ctx->mod);
    output->length = input->length;
}


void nmod_polydr_gcd(nmod_polydr_t G, 
                             const nmod_polydr_t A, const nmod_polydr_t B,
                                                          const nmod_ctx_t ctx)
{
    if (A->length < B->length)
    {
        nmod_polydr_gcd(G, B, A, ctx);
    }
    else /* lenA >= lenB >= 0 */
    {
        slong lenA = A->length, lenB = B->length, lenG;
        nmod_polydr_t tG;
        mp_ptr g;

        if (lenA == 0) /* lenA = lenB = 0 */
        {
            nmod_polydr_zero(G, ctx);
        } 
        else if (lenB == 0) /* lenA > lenB = 0 */
        {
            nmod_polydr_make_monic(G, A, ctx);
        }
        else /* lenA >= lenB >= 1 */
        {
            if (G == A || G == B)
            {
                nmod_polydr_init2(tG, FLINT_MIN(lenA, lenB), ctx);
                g = tG->coeffs;
            }
            else
            {
                nmod_polydr_fit_length(G, FLINT_MIN(lenA, lenB), ctx);
                g = G->coeffs;
            }

            lenG = _nmod_poly_gcd(g, A->coeffs, lenA,
                                               B->coeffs, lenB, ctx->mod);

            if (G == A || G == B)
            {
                nmod_polydr_swap(tG, G, ctx);
                nmod_polydr_clear(tG, ctx);
            }
            G->length = lenG;

            FLINT_ASSERT(ctx->mod.n > 1);

            if (G->length == 1)
                G->coeffs[0] = 1;
            else
                nmod_polydr_make_monic(G, G, ctx);
        }
    }
}

