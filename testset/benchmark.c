#include <fix_rfft.h>
#include <rfft.h>
#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <stdint.h>

void fixed_fft_benchmark()
{

    puts("-=-=-=-=-=-= Real FFT benchmark - fixed point =-=-=-=-=-=-");

    int l2n, n;
    int32_t twiddle[2048];
    int32_t data[2048];
    int bitrev[1024];

    puts("Generating sample data.");
    for (int i = 0; i < 2048; i++)
    {
        data[i] = random() << 1;
    }

    for (l2n = 3; l2n < 12; l2n++)
    {
        n = 1 << l2n;

        printf("%4d-points real-valued fixed-point FFT: ", n);
        precompute_bitrev_table_rfft_fix(bitrev, l2n);
        precompute_twiddle_factors_rfft_fix(twiddle, l2n);
        unsigned int start = clock();
        for (int i = 0; i < 5000; i++)
        {
            rfft_fix(data, l2n, twiddle, bitrev);
            irfft_fix(data, l2n, twiddle, bitrev);
        }
        unsigned int elapsed = clock() - start;
        printf("%12.5fus per transform\n", (float)elapsed / 10000);
    }
}

void fft_benchmark()
{

    puts("-=-=-=-=-=-= Real FFT benchmark - floating point =-=-=-=-=-=-");

    int l2n, n;
    float twiddle[2048];
    float data[2048];
    int bitrev[1024];

    puts("Generating sample data.");
    for (int i = 0; i < 2048; i++)
    {
        data[i] = (float)(random() - 0x40000000) / 0x80000000l;
    }

    for (l2n = 3; l2n < 12; l2n++)
    {
        n = 1 << l2n;

        printf("%4d-points real-valued fixed-point FFT: ", n);
        precompute_bitrev_table_rfft(bitrev, n);
        precompute_twiddle_factors_rfft(twiddle, n);
        unsigned int start = clock();
        for (int i = 0; i < 5000; i++)
        {
            rfft(data, n, twiddle, bitrev);
            irfft(data, n, twiddle, bitrev);
        }
        unsigned int elapsed = clock() - start;
        printf("%12.5fus per transform\n", (float)elapsed / 10000);
    }
}
