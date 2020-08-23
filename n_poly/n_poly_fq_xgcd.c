/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "n_poly.h"
#undef WANT_ASSERT
#define WANT_ASSERT 0
#undef FLINT_ASSERT 
#define FLINT_ASSERT(param)
void n_poly_fq_xgcd(
    n_poly_t G,
    n_poly_t S,
    n_poly_t T,
    const n_poly_t A,
    const n_poly_t B,
    const fq_nmod_ctx_t ctx)
{
    fq_nmod_poly_t GG, SS, TT, AA, BB;
    fq_nmod_poly_init(GG, ctx);
    fq_nmod_poly_init(SS, ctx);
    fq_nmod_poly_init(TT, ctx);
    fq_nmod_poly_init(AA, ctx);
    fq_nmod_poly_init(BB, ctx);
    n_poly_fq_get_fq_nmod_poly(AA, A, ctx);
    n_poly_fq_get_fq_nmod_poly(BB, B, ctx);
    fq_nmod_poly_xgcd(GG, SS, TT, AA, BB, ctx);
    n_poly_fq_set_fq_nmod_poly(G, GG, ctx);
    n_poly_fq_set_fq_nmod_poly(S, SS, ctx);
    n_poly_fq_set_fq_nmod_poly(T, TT, ctx);
    fq_nmod_poly_clear(GG, ctx);
    fq_nmod_poly_clear(SS, ctx);
    fq_nmod_poly_clear(TT, ctx);
    fq_nmod_poly_clear(AA, ctx);
    fq_nmod_poly_clear(BB, ctx);
}
