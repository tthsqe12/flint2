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
The conversion to Horner form can be stated as recursive. However, the call
stack has depth proportial to the length of the input polynomial in the worst
case. Therefore, we must convert it to an iterative algorithm.

The proceedure is

HornerForm(f):

    if f is simple to evaluate

        return eval(f)

    else
        choose a variable v and the smallest non zero exponent e appearing
            in the terms of f

        write f = q * v^e + r  where r is independent of the variable v

        return  HornerForm(q) * v^e + HornerForm(r)
*/

typedef struct
{
    slong f;
    slong r;
    slong v_var;
    ulong v_exp;
    int ret;
} stack_entry_struct;

typedef stack_entry_struct stack_entry_t[1];


/* A = A * X^pow */
void _nmod_mpoly_pmul(nmod_mpoly_t A, const nmod_mpoly_t X, slong pow,
                                          nmod_mpoly_t T, nmod_mpoly_ctx_t ctx)
{
    FLINT_ASSERT(pow >= 1);

    if (A->length <= X->length)
    {
        nmod_mpoly_pow_si(T, X, pow, ctx);
        nmod_mpoly_mul(A, A, T, ctx);
    } else
    {
        while (pow >= 1)
        {
            nmod_mpoly_mul(T, A, X, ctx);
            if (pow == 1)
            {
                nmod_mpoly_swap(A, T, ctx);
                break;
            }
            nmod_mpoly_mul(A, T, X, ctx);
            pow -= 2;
        }
    }
}

/*
    evaluate a f(xbar) at xbar = polys2,
*/
void _nmod_mpoly_compose_nmod_mpoly(nmod_mpoly_t res, nmod_mpoly_t poly1,
     nmod_mpoly_struct ** polys2, nmod_mpoly_ctx_t ctx1, nmod_mpoly_ctx_t ctx2)
{
    int ret;
    slong nvars = ctx1->minfo->nvars;
    slong i, j, k, N, cur, next, f, r, f_prev, r_prev, v, bits;
    slong sp, rp;
    stack_entry_struct * stack;
    nmod_mpoly_struct * regs;
    nmod_mpoly_t temp;
    slong * rtypes;
    ulong e, score;
    ulong totalcounts, maxcounts;
    ulong * counts, * mdegs;
    slong p_len;
    slong * p_list;
    mp_limb_t * p_coeff;
    ulong * p_exp, * p_uexp;

    TMP_INIT;

    p_len = poly1->length;
    p_coeff = poly1->coeffs;
    p_exp = poly1->exps;
    bits = poly1->bits;

    FLINT_ASSERT(res != poly1);
    FLINT_ASSERT(p_len > 0);

    TMP_START;

    /* unpack poly2 exponents */
    N = mpoly_words_per_exp(bits, ctx1->minfo);
    p_uexp = (ulong *) TMP_ALLOC(nvars*p_len*sizeof(ulong));
    for (i = 0; i < p_len; i++)
        mpoly_get_monomial_ui(p_uexp + nvars*i, p_exp + N*i, bits, ctx1->minfo);

    counts = (ulong *) TMP_ALLOC(nvars*sizeof(ulong));
    mdegs = (ulong * ) TMP_ALLOC(nvars*sizeof(ulong));

    /* stack */
    sp = -WORD(1); /* start with empty stack */
    stack = (stack_entry_struct *) TMP_ALLOC(nvars*(p_len + 1)*sizeof(stack_entry_struct));
    p_list = (slong *) TMP_ALLOC(p_len*sizeof(slong));

    /* registers of polynomials */
    rp = 0;
    rtypes = (slong *) TMP_ALLOC((nvars + 1)*sizeof(slong));
    regs   = (nmod_mpoly_struct *) TMP_ALLOC(nvars*sizeof(nmod_mpoly_struct));
    for (i = 0; i < nvars; i++)
        nmod_mpoly_init(regs + i, ctx2);
    nmod_mpoly_init(temp, ctx2);

    /* polynomials will be stored as link lists */
    for (i = 0; i + 1 < p_len; i++)
        p_list[i] = i + 1;
    p_list[i] = -WORD(1);

    sp++;
    (stack + sp)->ret = 0;
    (stack + sp)->f = 0;
HornerForm:

    f = (stack + sp)->f;

    FLINT_ASSERT(f != -WORD(1)); /* f is not supposed to be zero */

    /* obtain a count of the number of terms containing each variable */
    for (i = 0; i < nvars; i++)
    {
        counts[i] = 0;
        mdegs[i] = -WORD(1);
    }

    for (j = f; j != -WORD(1); j = p_list[j])
    {
        for (i = 0; i < nvars; i++)
        {
            e = p_uexp[nvars*j + i];
            if (e != 0)
            {
                counts[i]++;
                mdegs[i] = FLINT_MIN(mdegs[i], e);
            }
        }
    }

    totalcounts = 0;
    maxcounts = 0;
    v = -WORD(1);
    for (i = 0; i < nvars; i++)
    {
        maxcounts = FLINT_MAX(maxcounts, counts[i]);
        totalcounts += counts[i];
        if (counts[i] != 0)
            v = i;
    }


    /* handle simple cases */
    if (totalcounts == 0)
    {
        FLINT_ASSERT(p_list[f] == -WORD(1));    /* f should have had only one term */

        rtypes[rp] = f;

        goto HornerFormReturn;
        
    } else if (totalcounts == 1)
    {
        e = p_uexp[nvars*f + v];
        FLINT_ASSERT(e != 0);             /* this term should not be a scalar */
        nmod_mpoly_pow_si(regs + rp, polys2[v], e, ctx2);
        nmod_mpoly_scalar_mul_ui(regs + rp, regs + rp, p_coeff[f], ctx2);

        if (p_list[f] != -WORD(1)) /* if f has a second term */
        {
            FLINT_ASSERT(p_uexp[nvars*p_list[f] + v] == 0);   /* this term should be a scalar */
            nmod_mpoly_add_ui(regs + rp, regs + rp,  p_coeff[p_list[f]], ctx2);
        }   

        rtypes[rp] = -WORD(1);

        goto HornerFormReturn;
    }

    /* pick best power to pull out */
    k = 0;
    if (maxcounts == 1)
    {
        score = -WORD(1);
        for (i = 0; i < nvars; i++)
        {
            if (counts[i] == 1 && mdegs[i] < score)
            {
                score = mdegs[i];
                k = i;
            }
        }
    } else
    {
        score = 0;
        for (i = 0; i < nvars; i++)
        {
            if (counts[i] > 1 && mdegs[i]*(counts[i] - 1) > score)
            {
                score = mdegs[i]*(counts[i] - 1);
                k = i;
            }
        }
    }

    /* set variable power v */
    (stack + sp)->v_var = k;
    (stack + sp)->v_exp = mdegs[k];

    /* scan f and split into q and v with f = q*v + r then set f = q */
    r = -WORD(1);
    cur = f;
    f_prev = -WORD(1);
    r_prev = -WORD(1);
    while (cur != -WORD(1))
    {
        next = p_list[cur];
        e = p_uexp[nvars*cur + k];
        if (e == 0)
        {
            if (f_prev == -WORD(1))
                f = p_list[cur];
            else
                p_list[f_prev] = p_list[cur];

            if (r_prev == -WORD(1))
                r = cur;
            else
                p_list[r_prev] = cur;

            p_list[cur] = -WORD(1);
            r_prev = cur;
        } else
        {
            FLINT_ASSERT(e >= mdegs[k]);    /* mdegs[k] should be minimum non zero exponent */
            p_uexp[nvars*cur + k] = e - mdegs[k];
            f_prev = cur;
        }
        cur = next;
    }
    (stack + sp)->r = r;


    /* convert the quotient */
    sp++;
    (stack + sp)->ret = 1;
    (stack + sp)->f = f;
    goto HornerForm;
HornerForm1:

    /* convert the remainder */
    r = (stack + sp)->r;
    if (r != -WORD(1))
    {
        /* remainder is non zero */
        rp++;
        sp++;
        (stack + sp)->ret = 2;
        (stack + sp)->f = r;
        goto HornerForm;
HornerForm2:

        if (rtypes[rp - 1] == -WORD(1) && rtypes[rp] == -WORD(1))
        {
            /* both quotient and remainder are polynomials */
            _nmod_mpoly_pmul(regs + rp - 1, polys2[(stack + sp)->v_var], (stack + sp)->v_exp, temp, ctx2);
            nmod_mpoly_add(temp, regs + rp - 1, regs + rp, ctx2);
            nmod_mpoly_swap(temp, regs + rp - 1, ctx2);
        } else if (rtypes[rp - 1] == -WORD(1) && rtypes[rp] != -WORD(1))
        {
            /* quotient is a polynomial, remainder is a scalar */
            _nmod_mpoly_pmul(regs + rp - 1, polys2[(stack + sp)->v_var], (stack + sp)->v_exp, temp, ctx2);
            nmod_mpoly_add_ui(regs + rp - 1, regs + rp - 1, p_coeff[rtypes[rp]], ctx2);
        } else if (rtypes[rp - 1] != -WORD(1) && rtypes[rp] == -WORD(1))
        {
            /* quotient is a scalar, remainder is a polynomial */
            nmod_mpoly_pow_si(temp, polys2[(stack + sp)->v_var], (stack + sp)->v_exp, ctx2);
            nmod_mpoly_scalar_mul_ui(temp, temp, p_coeff[rtypes[rp - 1]], ctx2);
            nmod_mpoly_add(regs + rp - 1, temp, regs + rp, ctx2);
        } else
        {
            /* quotient is a scalar, remainder is a scalar */
            FLINT_ASSERT(0);    /* this should have been handled by simple case */
        }
        rp--;        

    } else
    {
        /* remainder is zero */
        if (rtypes[rp] == -WORD(1))
        {
            /* quotient is a polynomial */
            _nmod_mpoly_pmul(regs + rp, polys2[(stack + sp)->v_var], (stack + sp)->v_exp, temp, ctx2);

        } else
        {
            /* quotient is a scalar */
            FLINT_ASSERT(0);    /* this should have been handled by simple case */
        }
    }

    rtypes[rp] = -WORD(1);


HornerFormReturn:
    ret = (stack + sp)->ret;
    sp--;
    if (ret == 1) goto HornerForm1;
    if (ret == 2) goto HornerForm2;


    FLINT_ASSERT(rp == 0);
    FLINT_ASSERT(sp == -WORD(1));

    if (rtypes[rp] == -WORD(1))
    {
        nmod_mpoly_swap(res, regs + rp, ctx2);
    } else
    {
        nmod_mpoly_set_ui(res, p_coeff[rtypes[rp]], ctx2);
    }

    for (i = 0; i < nvars; i++)
        nmod_mpoly_clear(regs + i, ctx2);
    nmod_mpoly_clear(temp, ctx2);

    TMP_END;

    return;
}


void nmod_mpoly_compose_nmod_mpoly(nmod_mpoly_t A, nmod_mpoly_t B,
         nmod_mpoly_struct ** C, nmod_mpoly_ctx_t ctxB, nmod_mpoly_ctx_t ctxAC)
{
    FLINT_ASSERT(A != B);
    if (B->length == 0)
    {
        nmod_mpoly_zero(A, ctxAC);
        return;
    }
    else if (B->bits <= FLINT_BITS)
    {
        _nmod_mpoly_compose_nmod_mpoly(A, B, C, ctxB, ctxAC);
        return;
    }
    else
    {
        assert(0);
    }
}
