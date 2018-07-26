/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "nmod_mpoly.h"



/*
    Try to set G to the gcd of A and B given the form f of G.
    return codes as enumerated in nmod_mpoly.h:

    nmod_mpoly_sgcd_success,
    nmod_mpoly_sgcd_form_wrong,
    nmod_mpoly_sgcd_no_solution,
    nmod_mpoly_sgcd_scales_not_found,
    nmod_mpoly_sgcd_eval_point_not_found,
    nmod_mpoly_sgcd_eval_gcd_deg_too_high
*/
nmod_sgcd_ret_t fq_nmod_mpolyu_sgcd_zippel(
    fq_nmod_mpolyu_t G,
    fq_nmod_mpolyu_t A,
    fq_nmod_mpolyu_t B,
    fq_nmod_mpolyu_t f,
    slong var,
    fq_nmod_mpoly_ctx_t ctx,
    flint_rand_t randstate,
    slong * degbound)
{
    return 0;
}



int fq_nmod_mpolyu_pgcd_zippel(
    fq_nmod_mpolyu_t G,
    fq_nmod_mpolyu_t A,
    fq_nmod_mpolyu_t B,
    slong var,
    fq_nmod_mpoly_ctx_t ctx,
    mpoly_zipinfo_t zinfo,
    flint_rand_t randstate)
{
    return 0;
}

