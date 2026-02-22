#include "test_set.h"
#include <stdint.h>

void print_rfft(float *data, int n)
{
    printf("%10.5f%+10.5fi\n", data[0], 0.0);
    for (int i = 1; i < n / 2; i++)
    {
        printf("%10.5f%+10.5fi\n", data[2 * i], data[2 * i + 1]);
    }
    printf("%10.5f%+10.5fi\n", data[1], 0.0);
    for (int i = n / 2 - 1; i > 0; i--)
    {
        printf("%10.5f%+10.5fi\n", data[2 * i], -data[2 * i + 1]);
    }
}

void real_fft_test()
{
    int n = 16;                                                                                              // Change this value for different sizes of the arrays (N)
    float data[] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, -8.0, -7.0, -6.0, -5.0, -4.0, -3.0, -2.0, -1.0}; // interleaved real and imaginary parts
    float twiddle[n];
    int bitrev[n / 2];

    printf("-=-=-=-=-=-= Real FFT test - floating point =-=-=-=-=-=-\n");
    n = 16;

    // Precompute twiddle factors
    precompute_twiddle_factors_rfft(twiddle, n);
    precompute_bitrev_table_rfft(bitrev, n);

    rfft(data, n, twiddle, bitrev);
    printf("FFT result:\n");
    print_rfft(data, n);

    irfft(data, n, twiddle, bitrev);
    printf("IFFT result:\n");
    for (int i = 0; i < n; i++)
    {
        printf("%10.5f\n", data[i]);
    }
}

void print_fix_rfft(int32_t *data, int n)
{
    printf("%10.5f%+10.5fi\n",
           (float)data[0] / 0x800000, 0.0);
    for (int i = 1; i < n / 2; i++)
    {
        printf("%10.5f%+10.5fi\n",
               (float)data[2 * i] / 0x800000,
               (float)data[2 * i + 1] / 0x800000);
    }
    printf("%10.5f%+10.5fi\n",
           (float)data[1] / 0x800000, 0.0);
    for (int i = n / 2 - 1; i > 0; i--)
    {
        printf("%10.5f%+10.5fi\n",
               (float)data[2 * i] / 0x800000,
               -(float)data[2 * i + 1] / 0x800000);
    }
}

void real_fft_test_fix()
{
    const int n = 16, l2n = 4;                                                   // Change this value for different sizes of the arrays (N)
    int32_t data[] = {1, 2, 3, 4, 5, 6, 7, 8, -8, -7, -6, -5, -4, -3, -2, -1}; // interleaved real and imaginary parts
    int32_t twiddle[n];
    int bitrev[n / 2];

    for (int i = 0; i < n; i++)
    {
        data[i] = data[i] << 23;
    }

    printf("-=-=-=-=-=-= Real FFT test - fixed point =-=-=-=-=-=-\n");

    // Precompute twiddle factors
    precompute_twiddle_factors_rfft_fix(twiddle, l2n);
    precompute_bitrev_table_rfft_fix(bitrev, l2n);

    rfft_fix(data, l2n, twiddle, bitrev);
    printf("FFT result:\n");

    print_fix_rfft(data, n);

    irfft_fix(data, l2n, twiddle, bitrev);
    printf("IFFT result:\n");
    for (int i = 0; i < n; i++)
    {
        printf("%10.5f\n", (float)data[i] / 0x800000);
    }
}
