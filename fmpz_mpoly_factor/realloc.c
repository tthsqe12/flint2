/*
    Copyright (C) 2019 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fmpz_mpoly_factor.h"

void fmpz_mpoly_factor_realloc(fmpz_mpoly_factor_t fac, slong alloc, const fmpz_mpoly_ctx_t ctx)
{
    if (alloc == 0)             /* Clear up, reinitialise */
    {
        fmpz_mpoly_factor_clear(fac, ctx);
        fmpz_mpoly_factor_init(fac, ctx);
    }
    else if (fac->alloc > 0)            /* Realloc */
    {
        if (fac->alloc > alloc)
        {
            slong i;

            for (i = alloc; i < fac->length; i++)
                fmpz_mpoly_clear(fac->poly + i, ctx);

            fac->poly = flint_realloc(fac->poly, alloc * sizeof(fmpz_mpoly_struct));
            fac->exp  = flint_realloc(fac->exp, alloc * sizeof(slong));
            fac->alloc     = alloc;
        }
        else if (fac->alloc < alloc)
        {
            slong i;

            fac->poly = flint_realloc(fac->poly, alloc * sizeof(fmpz_mpoly_struct));
            fac->exp  = flint_realloc(fac->exp, alloc * sizeof(slong));

            for (i = fac->alloc; i < alloc; i++)
            {
                fmpz_mpoly_init(fac->poly + i, ctx);
                fac->exp[i] = WORD(0);
            }
            fac->alloc = alloc;
        }
    }
    else                        /* Nothing allocated already so do it now */
    {
        slong i;

        fac->poly = flint_malloc(alloc * sizeof(fmpz_mpoly_struct));
        fac->exp  = flint_calloc(alloc, sizeof(slong));

        for (i = 0; i < alloc; i++)
        {
            fmpz_mpoly_init(fac->poly + i, ctx);
        }
        fac->length = 0;
        fac->alloc  = alloc;
    }
}

