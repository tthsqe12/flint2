/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "n_poly.h"
#include "ui_factor.h"


void n_poly_fq_add(
    n_poly_t A,
    const n_poly_t B,
    const n_poly_t C,
    const fq_nmod_ctx_t ctx)
{
    slong d = fq_nmod_ctx_degree(ctx);
    slong i;
    slong Blen = B->length;
    slong Clen = C->length;

    if (Blen > Clen)
    {
        n_poly_fit_length(A, d*Blen);
        _nmod_vec_add(A->coeffs, B->coeffs, C->coeffs, d*Clen, ctx->mod);
        if (A != B)
            for (i = d*Clen; i < d*Blen; i++)
                A->coeffs[i] = B->coeffs[i];
        A->length = Blen;
    }
    else if (Blen < Clen)
    {
        n_poly_fit_length(A, d*Clen);
        _nmod_vec_add(A->coeffs, B->coeffs, C->coeffs, d*Blen, ctx->mod);
        if (A != C)
            for (i = d*Blen; i < d*Clen; i++)
                A->coeffs[i] = C->coeffs[i];
        A->length = Clen;
    }
    else
    {
        n_poly_fit_length(A, d*Blen);
        _nmod_vec_add(A->coeffs, B->coeffs, C->coeffs, d*Clen, ctx->mod);
        A->length = Clen;
        _n_poly_fq_normalise(A, d);
    }
}
