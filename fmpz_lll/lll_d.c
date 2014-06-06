/*=============================================================================

    This file is part of FLINT.

    FLINT is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    FLINT is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with FLINT; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301 USA

=============================================================================*/
/******************************************************************************

    Copyright (C) 2009, 2010 William Hart
    Copyright (C) 2009, 2010 Andy Novocin
    Copyright (C) 2014 Abhinav Baid

******************************************************************************/

#include "fmpz_lll.h"

#ifdef FUNC_NAME
#undef FUNC_NAME
#endif

#ifdef COMPUTE
#undef COMPUTE
#endif

#ifdef CALL_BABAI
#undef CALL_BABAI
#endif

#define FUNC_NAME fmpz_lll_d
#define COMPUTE(G, I, C)                                               \
do {                                                                   \
    d_mat_entry(G, I, I) = _d_vec_norm(appB->rows[I], C);              \
} while (0)
#define CALL_BABAI(NFF, BO, HF)                                        \
do {                                                                   \
    if (NFF < 50)                                                      \
    {                                                                  \
        BO =                                                           \
            fmpz_lll_check_babai(kappa, B, mu, r, s, appB, expo, A,    \
                                 alpha[kappa], zeros, kappamax,        \
                                 FLINT_MIN(kappamax + 1 + shift,       \
                                           n), fl);                    \
    }                                                                  \
    else                                                               \
    {                                                                  \
        BO = -1;                                                       \
    }                                                                  \
    if (BO == -1)                                                      \
    {                                                                  \
        NFF++;                                                         \
        HF =                                                           \
            fmpz_lll_check_babai_heuristic_d(kappa, B, mu, r, s,       \
                                             appB, expo, A,            \
                                             alpha[kappa], zeros,      \
                                             kappamax,                 \
                                             FLINT_MIN(kappamax +      \
                                                       1 + shift,      \
                                                       n), fl);        \
    }                                                                  \
} while (0)
#include "d_lll.c"
#undef FUNC_NAME
#undef COMPUTE
#undef CALL_BABAI
