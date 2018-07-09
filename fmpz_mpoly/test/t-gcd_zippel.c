/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fmpz_mpoly.h"

int
main(void)
{
    slong i, j;
    FLINT_TEST_INIT(state);

    flint_printf("gcd_zippel....\n");
    fflush(stdout);

    if (0) {
        fmpz_mpoly_ctx_t ctx;
        fmpz_mpoly_t a, b, g;
        const char* vars[] = {"x0", "x1", "x2", "x3", "x4", "x5", "X", "x7", "x8", "x9", "x10"};

        fmpz_mpoly_ctx_init(ctx, 7, ORD_DEGREVLEX);
        fmpz_mpoly_init(a, ctx);
        fmpz_mpoly_init(b, ctx);
        fmpz_mpoly_init(g, ctx);
        fmpz_mpoly_set_str_pretty(a, "(83990204133285629858347513265339849451642301433702400*x1^2*x2^5*x3^3*x4^3*x5^2)*X^4 + (-176150975980870807675449573304781876534107448349067428495360*x0*x2^2*x3^2*x5^5)*X^3 + (163456000774756675154428230454675375458806918985*x0*x1^2*x2^6*x3^3*x4^6*x5^2)*X^2 + (-342813002581959632373050285628449851205102552695701504*x0^2*x2^3*x3^2*x4^3*x5^5)*X^1", vars, ctx);
        fmpz_mpoly_set_str_pretty(b, "(407430976282525688796438911588561033779200*x1*x2^2*x3^3*x4^3)*X^4 + (1083137323863694000171584191997897600521957187972407744100160815035617204961280*x0^2*x2^3*x4^2+792914348311528082358983902127160255*x0*x1*x2^3*x3^3*x4^6)*X^2 + (2107927907493537453702130932433991145186542215179572751824224810308429792*x0^3*x2^4*x4^5)*X^0", vars, ctx);

        fmpz_mpoly_gcd_zippel(g, a, b, ctx);
        fmpz_mpoly_assert_canonical(g, ctx);

printf("g: "); fmpz_mpoly_print_pretty(g, vars, ctx); printf("\n");


        fmpz_mpoly_clear(g, ctx);
        fmpz_mpoly_clear(a, ctx);
        fmpz_mpoly_clear(b, ctx);
        fmpz_mpoly_ctx_clear(ctx);
    }


    if (0) {
        fmpz_mpoly_ctx_t ctx;
        fmpz_mpoly_t a, b, g;
        const char* vars[] = {"x0", "x1", "x2", "x3", "x4", "x5", "x6", "x7", "x8", "x9", "x10"};

        fmpz_mpoly_ctx_init(ctx, 7, ORD_DEGREVLEX);
        fmpz_mpoly_init(a, ctx);
        fmpz_mpoly_init(b, ctx);
        fmpz_mpoly_init(g, ctx);
        fmpz_mpoly_set_str_pretty(a, "-342813002581959632373050285628449851205102552695701504*x0^2*x2^3*x3^2*x4^3*x5^5*x6+163456000774756675154428230454675375458806918985*x0*x1^2*x2^6*x3^3*x4^6*x5^2*x6^2-176150975980870807675449573304781876534107448349067428495360*x0*x2^2*x3^2*x5^5*x6^3+83990204133285629858347513265339849451642301433702400*x1^2*x2^5*x3^3*x4^3*x5^2*x6^4", vars, ctx);
        fmpz_mpoly_set_str_pretty(b, "6323783722480612361106392797301973435559626645538718255472674430925289376*x0^3*x2^6*x3*x4^5*x5^2+3249411971591082000514752575993692801565871563917223232300482445106851614883840*x0^2*x2^5*x3*x4^2*x5^2*x6^2+2378743044934584247076951706381480765*x0*x1*x2^5*x3^4*x4^6*x5^2*x6^2+1222292928847577066389316734765683101337600*x1*x2^4*x3^4*x4^3*x5^2*x6^4", vars, ctx);

        fmpz_mpoly_gcd_zippel(g, a, b, ctx);
        fmpz_mpoly_assert_canonical(g, ctx);

printf("g: "); fmpz_mpoly_print_pretty(g, vars, ctx); printf("\n");

/*
        fmpz_mpoly_gcd_brown(g, a, b, ctx);
        fmpz_mpoly_assert_canonical(g, ctx);
printf("g: "); fmpz_mpoly_print_pretty(g, vars, ctx); printf("\n");
*/

        fmpz_mpoly_clear(g, ctx);
        fmpz_mpoly_clear(a, ctx);
        fmpz_mpoly_clear(b, ctx);
        fmpz_mpoly_ctx_clear(ctx);
    }


    if (0) {
        fmpz_mpoly_ctx_t ctx;
        fmpz_mpoly_t a, b, g;
        const char* vars[] = {"x0", "x1", "x2", "x3", "x4", "x5", "x6", "x7", "x8", "x9", "x10"};

        fmpz_mpoly_ctx_init(ctx, 8, ORD_DEGREVLEX);
        fmpz_mpoly_init(a, ctx);
        fmpz_mpoly_init(b, ctx);
        fmpz_mpoly_init(g, ctx);
        fmpz_mpoly_set_str_pretty(a, "-7020808969457975896006939690168856095039316742934942973952*x0^4*x1^3*x2^6*x3^4*x4^4*x5^4*x6^2*x7^2-83531114730539976333193285975025691989061451433757*x0^3*x1^2*x2^6*x3^2*x4^5*x5^4*x6^3*x7^3+163456000774756675154428230454675375458806918985*x0*x1^2*x2^6*x3^3*x4^6*x5^2*x6^2*x7^5-3607571308062036012904959793294629717208968373244713112363335680*x0^3*x1^3*x2^5*x3^4*x4*x5^4*x6^4*x7^2-42921614039528361739448939843705473230946675109931386880*x0^2*x1^2*x2^5*x3^2*x4^2*x5^4*x6^5*x7^3+2378743044934584247076951706381480765*x0*x1*x2^5*x3^4*x4^6*x5^2*x6^2*x7^4+83990204133285629858347513265339849451642301433702400*x1^2*x2^5*x3^3*x4^3*x5^2*x6^4*x7^5+1222292928847577066389316734765683101337600*x1*x2^4*x3^4*x4^3*x5^2*x6^4*x7^4+6323783722480612361106392797301973435559626645538718255472674430925289376*x0^3*x2^6*x3*x4^5*x5^2*x7^4-342813002581959632373050285628449851205102552695701504*x0^2*x2^3*x3^2*x4^3*x5^5*x6*x7^5+3249411971591082000514752575993692801565871563917223232300482445106851614883840*x0^2*x2^5*x3*x4^2*x5^2*x6^2*x7^4-176150975980870807675449573304781876534107448349067428495360*x0*x2^2*x3^2*x5^5*x6^3*x7^5", vars, ctx);
        fmpz_mpoly_set_str_pretty(b, "290373907934093468189524695273411*x0^3*x1^2*x2^4*x3^4*x4^6*x5^2*x6^2*x7^4+149205680346798535031391328420602378240*x0^2*x1^2*x2^3*x3^4*x4^3*x5^2*x6^4*x7^4+706737904671899891757015097398939232284907987892413049546968399872*x0^2*x1^2*x2^3*x3^2*x4^6*x5^4*x7^3-423490795964240612469074109*x0^2*x1*x2^6*x3^3*x4^4*x5^3*x6*x7^2-38059888860737414552607525128334999552*x0*x2^5*x3^4*x4^3*x5^5*x7^2+363150086878245386240831811501920174903238965007028646131547769384468480*x0*x1^2*x2^2*x3^2*x4^3*x5^4*x6^2*x7^3-217606439855448202519761365698560*x0*x1*x2^5*x3^3*x4*x5^3*x6^3*x7^2-19556686934415534159066284250478551949639680*x2^4*x3^4*x5^5*x6^2*x7^2", vars, ctx);

        fmpz_mpoly_gcd_zippel(g, a, b, ctx);
        fmpz_mpoly_assert_canonical(g, ctx);

printf("g: "); fmpz_mpoly_print_pretty(g, vars, ctx); printf("\n");

/*
        fmpz_mpoly_gcd_brown(g, a, b, ctx);
        fmpz_mpoly_assert_canonical(g, ctx);
printf("g: "); fmpz_mpoly_print_pretty(g, vars, ctx); printf("\n");
*/

        fmpz_mpoly_clear(g, ctx);
        fmpz_mpoly_clear(a, ctx);
        fmpz_mpoly_clear(b, ctx);
        fmpz_mpoly_ctx_clear(ctx);
    }

    /* examples from Zippel's 1979 paper */
    if (0) {
        fmpz_mpoly_ctx_t ctx;
        fmpz_mpoly_t r, d, f, g;
        int success;
        const char* vars[] = {"x1", "x2", "x3", "x4", "x5", "x6", "x7", "x8", "x9", "x10"};

        const char * example[][3] =
        {{
            "x1^2 + x1 + 3",
            "2*x1^2 + 2*x1 + 1",
            "x1^2 + 2*x1 + 2"
        }, {
            "2*x1^2*x2^2 + x1*x2 + 2*x1",
            "x2^2 + 2*x1^2*x2 + x1^2 + 1",
            "x1^2*x2^2 + x1^2*x2 + x1*x2 + x1^2 + x1"
        }, {
            "x2^2*x3^2 + x2^2*x3 + 2*x1^2*x2*x3 + x1*x3",
            "x3^2 + x2^2*x3 + x1^2*x2*x3 + x1*x3 + x1^2*x2^2",
            "x2*x3 + 2*x1*x3 + x3 + x1"
        }, {
            "x1^2*x4^2 + x2^2*x3*x4 + x1^2*x2*x4 + x2*x4 + x1^2*x2*x3",
            "x1*x2*x3^2*x4^2 + x1*x3^2*x4^2 + x1*x4^2 + x4^2 + x1*x3*x4",
            "x1*x3^2*x4^2 + x3^2*x4^2 + x4^2 + x1*x2^2*x3*x4 + x1*x2^2"
        }, {
            "x1^3*x2^2*x3^2*x4*x5^2 + x1*x2^2*x5^2 + x1^3*x3*x4^2*x5"
                                    " + x1^3*x2*x3^2*x4*x5 + x1^2*x2*x3^2*x4^2"
            ,
            "x1*x2^2*x5^2 + x1*x2*x3^2*x4*x5 + x1*x2*x3^2*x4^2"
                                                          " + x1*x2^2*x4^2 + 1"
            ,
            "x1*x3^2*x4*x5^2 + x2*x5^2 + x1*x2*x4*x5 + x2*x5 + x1*x2*x3*x4^2"
        }, {
            "x1*x2*x4^2*x5^2*x6^2 + x1*x2^2*x3^2*x4*x5^2*x6^2  + x1^2*x3*x6^2"
                                   " + x1^2*x2*x3^2*x4*x5^2*x6 + x1^2*x3*x5*x6"
            ,
            "x1^2*x2*x4*x5^2*x6^2 + x1*x3*x5^2*x6^2 + x1*x2^2*x6^2"
                                      " + x1^2*x2^2*x3^2*x5*x6 + x1*x3^2*x4*x5"
            ,
            "x2^2*x3^2*x4*x5^2*x6 + x1*x4^2*x5*x6 + x2^2*x3^2*x4*x5*x6"
                                         " + x1*x2^2*x3*x4^2*x6 + x1^2*x3*x5^2"
        }, {
            "x1*x2^2*x4^2*x6^2*x7^2 + x1^2*x3*x4*x6^2*x7^2 + x3^2*x4^2*x7^2"
                                              " + x1^2*x2*x4^2*x6 + x3*x4*x5^2"
            ,
            "x1^2*x2*x4^2*x5*x6^2*x7^2 + x1*x2*x3*x6*x7 + x3*x4^2*x5^2*x7"
                                   " + x1*x4^2*x5^2*x7 + x1^2*x2*x3*x4^2+x5*x6"
            ,
            "x1*x3*x5*x6^2*x7^2 + x2^2*x3^2*x4^2*x5*x6*x7^2 + x4*x6*x7^2"
                                   " + x1^2*x2*x3*x5*x6*x7 + x1^2*x3^2*x4*x5^2"
        }, {
            "x2^2*x4*x5*x6*x7*x8^2 + x1^2*x2*x3^2*x4^2*x6^2*x7^2*x8"
                   " + x1^2*x3*x4^2*x6^2*x7^2 + x1^2*x2^2*x3^2*x4*x5^2*x6*x7^2"
                                                                " + x2^2*x4*x6"
            ,
            "x1^2*x2^2*x3*x4^2*x5*x6^2*x8^2 + x2*x5*x6^2*x8^2"
              " + x1^2*x2^2*x3^2*x4^2*x6^2*x7^2*x8 + x1^2*x3^2*x4*x5^2*x7^2*x8"
                                                      " + x1*x2^2*x3^2*x5^2*x7"
            ,
            "x1*x4^2*x5*x6*x7*x8^2 + x1*x2^2*x4^2*x5^2*x6^2*x8"
                       " + x1^2*x2*x3*x4^2*x6^2*x8 + x1^2*x2^2*x3^2*x4*x5^2*x8"
                                                           " + x1*x2*x4^2*x5^2"
        }, {
            "x1^2*x3^3*x4*x6*x8*x9^2 + x1*x2*x3*x4^2*x5^2*x8*x9"
                      " + x2*x3*x4*x5^2*x8*x9 + x1*x3^3*x4^2*x5^2*x6^2*x7*x8^2"
                                                  " + x2*x3*x4*x5^2*x6*x7*x8^2"
            ,
            "x1^2*x2^2*x3*x7^2*x8*x9 + x2^2*x9 + x1^2*x3*x4^2*x5^2*x6*x7^2"
                                            " + x4^2*x5^2*x7^2 + x3*x4^2*x6*x7"
            ,
            "x1^2*x2*x4*x5*x6*x7^2*x8^2*x9^2 + x1^2*x2*x3*x5*x6^2*x7^2*x8*x9^2"
                              " + x1^2*x3*x4*x6*x7^2*x8*x9 + x1^2*x2^2*x6*x8^2"
                                                        " + x2^2*x4*x5*x6^2*x7"
        }, {
            "x1*x2^2*x4^2*x8*x9^2*x10^2 + x2^2*x4*x5^2*x6*x7*x9*x10^2"
                        " + x1^2*x2*x3*x5^2*x7^2*x9^2 + x1*x3^2*x4^2*x7^2*x9^2"
                                                      " + x1^2*x3*x4*x7^2*x8^2"
            ,
            "x1*x2*x3^2*x4*x6*x7*x8*x9^2*x10^2 + x2^2*x3^2*x4^2*x6^2*x9*x10^2"
                                    " + x1*x2^2*x3^2*x4*x5*x6*x7*x8^2*x9^2*x10"
               " + x1^2*x2*x4^2*x5^2*x8^2*x9^2*x10 + x3*x4^2*x5*x6*x7^2*x9*x10"
            ,
            "x1*x2^2*x3^2*x5^2*x6^2*x7*x8*x9^2*x10^2 + x3*x8*x9^2*x10^2"
                  " + x1*x2^2*x3*x4*x5^2*x6^2*x8^2*x9*x10 + x1*x3*x6*x7*x8*x10"
                                                    " + x4^2*x5^2*x6^2*x7*x9^2"
        }};


        for (i = 1; i <= 10; i++)
        {
            fmpz_mpoly_ctx_init(ctx, i, ORD_DEGREVLEX);
            fmpz_mpoly_init(r, ctx);
            fmpz_mpoly_init(d, ctx);
            fmpz_mpoly_init(f, ctx);
            fmpz_mpoly_init(g, ctx);
            fmpz_mpoly_set_str_pretty(d, example[i - 1][0], vars, ctx);
            fmpz_mpoly_set_str_pretty(f, example[i - 1][1], vars, ctx);
            fmpz_mpoly_set_str_pretty(g, example[i - 1][2], vars, ctx);
            fmpz_mpoly_mul_johnson(f, f, d, ctx);
            fmpz_mpoly_mul_johnson(g, g, d, ctx);
            success = fmpz_mpoly_gcd_zippel(r, f, g, ctx);
            if (!success || !fmpz_mpoly_equal(r, d, ctx))
            {
                flint_printf("FAIL\ncheck example %wd\n",i);
                flint_abort();
            }
            fmpz_mpoly_clear(r, ctx);
            fmpz_mpoly_clear(d, ctx);
            fmpz_mpoly_clear(f, ctx);
            fmpz_mpoly_clear(g, ctx);
            fmpz_mpoly_ctx_clear(ctx);
        }

    }


    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {

flint_printf("TESTING i = %d\n", i);

        fmpz_mpoly_ctx_t ctx;
        fmpz_mpoly_t a, b, g, ca, cb, cg, t;
        mp_bitcnt_t coeff_bits;
        slong len, len1, len2;
        slong degbound;
        int res;

        fmpz_mpoly_ctx_init_rand(ctx, state, 10);

        fmpz_mpoly_init(g, ctx);
        fmpz_mpoly_init(a, ctx);
        fmpz_mpoly_init(b, ctx);
        fmpz_mpoly_init(ca, ctx);
        fmpz_mpoly_init(cb, ctx);
        fmpz_mpoly_init(cg, ctx);
        fmpz_mpoly_init(t, ctx);

        len = n_randint(state, 10) + 1;
        len1 = n_randint(state, 10);
        len2 = n_randint(state, 10);

        degbound = 70/(2*ctx->minfo->nvars - 1);

        coeff_bits = n_randint(state, 200);

        for (j = 0; j < 4; j++)
        {
            do {
                fmpz_mpoly_randtest_bound(t, state, len, coeff_bits + 1, degbound, ctx);
            } while (t->length == 0);
            fmpz_mpoly_randtest_bound(a, state, len1, coeff_bits, degbound, ctx);
            fmpz_mpoly_randtest_bound(b, state, len2, coeff_bits, degbound, ctx);
            fmpz_mpoly_mul_johnson(a, a, t, ctx);
            fmpz_mpoly_mul_johnson(b, b, t, ctx);

            fmpz_mpoly_randtest_bits(g, state, len, coeff_bits, FLINT_BITS, ctx);

            res = fmpz_mpoly_gcd_zippel(g, a, b, ctx);
            fmpz_mpoly_assert_canonical(g, ctx);

            if (!res) {
                continue;
            }

            if (fmpz_mpoly_is_zero(g, ctx))
            {
                if (!fmpz_mpoly_is_zero(a, ctx) || !fmpz_mpoly_is_zero(b, ctx)) {
                    printf("FAIL\n");
                    flint_printf("Check zero gcd only results from zero inputs\ni = %wd, j = %wd\n", i ,j);
                    flint_abort();
                }
                continue;
            }

            if (fmpz_sgn(g->coeffs + 0) <= 0)
            {
                printf("FAIL\n");
                flint_printf("Check gcd has positive lc\ni = %wd, j = %wd\n", i ,j);
                flint_abort();
            }

            res = 1;
            res = res && fmpz_mpoly_divides_monagan_pearce(ca, a, g, ctx);
            res = res && fmpz_mpoly_divides_monagan_pearce(cb, b, g, ctx);
            if (!res)
            {
                printf("FAIL\n");
                flint_printf("Check divisibility\ni = %wd, j = %wd\n", i ,j);
                flint_abort();
            }

            res = fmpz_mpoly_gcd_zippel(cg, ca, cb, ctx);

            if (!res)
                continue;

            if (!fmpz_mpoly_equal_ui(cg, UWORD(1), ctx))
            {
                printf("FAIL\n");
                flint_printf("Check cofactors are relatively prime\ni = %wd, j = %wd\n", i ,j);                
                flint_abort();
            }
        }

        fmpz_mpoly_clear(g, ctx);
        fmpz_mpoly_clear(a, ctx);
        fmpz_mpoly_clear(b, ctx);
        fmpz_mpoly_clear(ca, ctx);
        fmpz_mpoly_clear(cb, ctx);
        fmpz_mpoly_clear(cg, ctx);
        fmpz_mpoly_clear(t, ctx);
        fmpz_mpoly_ctx_clear(ctx);
    }


    printf("PASS\n");
    FLINT_TEST_CLEANUP(state);

    return 0;
}

