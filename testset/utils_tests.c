#include "test_set.h"
#include <stdint.h>

void cordic_test()
{
    int n = 16;
    int32_t t, s, c;
    printf("-=-=-=-=-=-= Cordic test =-=-=-=-=-=-\n");
    for (int i = -n; i <= n; i++)
    {
        t = i * (0x80000000l >> 4); // 1/16
        cordic_sin_cos(t, &s, &c);
        printf("sin(%8.5f) = %8.5f, cos(%8.5f) = %8.5f\n",
               (double)t / 0x80000000l * M_PI / 2,
               (double)s / 0x80000000l,
               (double)t / 0x80000000l * M_PI / 2,
               (double)c / 0x80000000l);
    }
}

void twiddle_test()
{
    int n = 64;
    float twiddle[n];
    printf("-=-=-=-=-=-= Twiddle test =-=-=-=-=-=-\n");
    precompute_twiddle_factors(twiddle, n);
    for (int i = 0; i < n / 2; i++)
    {
        printf("w[%4d] = %10.5f%+10.5fi\n", i, twiddle[2 * i], twiddle[2 * i + 1]);
    }
}

void fix_twiddle_test()
{
    int n = 128, l2n = 7;
    int32_t twiddle[n];
    float twiddle_f[n];
    printf("-=-=-=-=-=-= Fix twiddle test =-=-=-=-=-=-\n");
    precompute_twiddle_factors_fix(twiddle, l2n);
    for (int i = 0; i < n / 2; i++)
    {
        printf("w[%4d] = %10.5f%+10.5fi\n",
               i,
               (double)twiddle[2 * i] / 0x80000000l,
               (double)twiddle[2 * i + 1] / 0x80000000l);
    }
    printf("-=-=-=-=-=-= Fix twiddle error test =-=-=-=-=-=-\n");

    precompute_twiddle_factors(twiddle_f, n);
    for (int i = 0; i < n / 2; i++)
    {
        double absError = hypot((double)twiddle[2 * i] / 0x80000000l - twiddle_f[2 * i],
                                (double)twiddle[2 * i + 1] / 0x80000000l - twiddle_f[2 * i + 1]),
               relError = absError / hypot(twiddle_f[2 * i], twiddle_f[2 * i + 1]);

        printf("Abs Sq Error[%4d] = %11.5e, ", i, absError);
        printf("Rel Sq Error[%4d] = %10.5f%% \n", i, relError * 100);
    }
}
