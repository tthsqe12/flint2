/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "nmod_mpoly_factor.h"


void zassenhaus_prune_init(zassenhaus_prune_t zas)
{
    zas->deg = 0;
    zas->pos_degs = NULL;
    zas->new_length = 0;
    zas->new_total = 0;
    zas->new_degs = NULL;
    zas->alloc = 0;
}

void zassenhaus_prune_clear(zassenhaus_prune_t zas)
{
    if (zas->alloc > 0)
    {
        flint_free(zas->pos_degs);
        flint_free(zas->new_degs);
    }
}

void zassenhaus_prune_set_degree(zassenhaus_prune_t zas, slong d)
{
    slong i;

    if (d < 1)
    {
        flint_throw(FLINT_ERROR, "zassenhaus_prune_start");
        return;
    }

    if (zas->alloc > 0)
    {
        zas->pos_degs = (unsigned char *) flint_realloc(zas->pos_degs, (d + 1)*sizeof(unsigned char));
        zas->new_degs = (slong *) flint_realloc(zas->new_degs, (d + 1)*sizeof(slong));
    }
    else
    {
        zas->pos_degs = (unsigned char *) flint_malloc((d + 1)*sizeof(unsigned char));
        zas->new_degs = (slong *) flint_malloc((d + 1)*sizeof(slong));        
    }
    zas->alloc = d + 1;
    zas->deg = d;

    for (i = 0; i <= d; i++)
        zas->pos_degs[i] = 1;

    zas->new_length = 0;
    zas->new_total = 0;
}

void zassenhaus_prune_start_add_factors(zassenhaus_prune_t zas)
{
    zas->new_length = 0;
    zas->new_total = 0;
}


void zassenhaus_prune_add_factor(zassenhaus_prune_t zas, slong deg, slong exp)
{
    slong i;

    if (exp < 1 || deg < 1)
        return;

    for (i = 0; i < exp; i++)
    {
        if (zas->new_length >= zas->deg)
        {
            flint_throw(FLINT_ERROR, "zassenhaus_prune_add_factor");
            return;
        }
        zas->new_total += deg;
        zas->new_degs[zas->new_length] = deg;
        zas->new_length++;
    }
}

void zassenhaus_prune_finish_add_factors(zassenhaus_prune_t zas)
{
    slong i, j;
    unsigned char * a = zas->pos_degs;
    unsigned char pos_mask = 1;
    unsigned char new_mask = 2;

    if (zas->new_total != zas->deg)
    {
        flint_throw(FLINT_ERROR, "zassenhaus_prune_add_factor");
    }

    a[0] |= new_mask;
    for (j = 1; j <= zas->deg; j++)
        a[j] &= ~new_mask;

    for (i = 0; i < zas->new_length; i++)
    {
        slong d = zas->new_degs[i];

        for (j = zas->deg; j >= 0; j--)
        {
            if ((a[j] & new_mask) != 0)
            {
                if (j + d > zas->deg)
                {
                    flint_throw(FLINT_ERROR, "zassenhaus_prune_add_factor");
                }
                a[j + d] |= new_mask;
            }
        }
    }

    for (j = 0; j <= zas->deg; j++)
        a[j] &= a[j] >> 1;

    if (a[0] != pos_mask || a[zas->deg] != pos_mask)
    {
        flint_throw(FLINT_ERROR, "zassenhaus_prune_add_factor");
    }
}

int zassenhaus_prune_is_irreducible(const zassenhaus_prune_t zas)
{
    slong i;

    for (i = 1; i < zas->deg; i++)
    {
        if (zas->pos_degs[i] != 0)
            return 0;
    }

    return 1;
}

int zassenhaus_prune_degree_is_possible(const zassenhaus_prune_t zas, slong d)
{
    if (d <= 0)
        return d = 0;

    if (d >= zas->deg)
        return d == zas->deg;

    return zas->pos_degs[d] != 0;
}
