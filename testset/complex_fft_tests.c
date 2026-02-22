#include "test_set.h"
#include <stdint.h>

void complex_fft_test()
{
    int n = 8;                                                                                               // Change this value for different sizes of the arrays (N)
    float data[] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, -8.0, -7.0, -6.0, -5.0, -4.0, -3.0, -2.0, -1.0}; // interleaved real and imaginary parts
    float twiddle[n];
    int bitrev[n];

    printf("-=-=-=-=-=-= Complex FFT test - floating point =-=-=-=-=-=-\n");
    // Precompute twiddle factors
    n = 8;
    precompute_twiddle_factors(twiddle, n);
    precompute_bitrev_table(bitrev, n);

    // Perform FFT
    fft(data, n, twiddle, bitrev);

    // Display the FFT results
    printf("FFT result:\n");
    for (int i = 0; i < n; i++)
    {
        printf("%10.5f%+10.5fi\n", data[2 * i], data[2 * i + 1]);
    }

    // Perform IFFT
    ifft(data, n, twiddle, bitrev);

    // Display the IFFT results
    printf("IFFT result:\n");
    for (int i = 0; i < n; i++)
    {
        printf("%10.5f%+10.5fi\n", data[2 * i], data[2 * i + 1]);
    }
}

void complex_fft_test_fix()
{
    const int n = 8, l2n = 3;                                                    // Change this value for different sizes of the arrays (N)
    int32_t data[] = {1, 2, 3, 4, 5, 6, 7, 8, -8, -7, -6, -5, -4, -3, -2, -1}; // interleaved real and imaginary parts
    int32_t twiddle[n];
    int bitrev[n];

    for (int i = 0; i < 2 * n; i++)
    {
        data[i] = data[i] << 23;
    }

    printf("-=-=-=-=-=-= Complex FFT test - fixed point =-=-=-=-=-=-\n");
    // Precompute twiddle factors
    precompute_twiddle_factors_fix(twiddle, l2n);
    precompute_bitrev_table_fix(bitrev, l2n);

    // Perform FFT
    fft_fix(data, l2n, twiddle, bitrev);

    // Display the FFT results
    printf("FFT result:\n");
    for (int i = 0; i < n; i++)
    {
        printf("%10.5f%+10.5fi\n",
               (float)data[2 * i] / 0x800000,
               (float)data[2 * i + 1] / 0x800000);
    }

    // Perform IFFT
    ifft_fix(data, l2n, twiddle, bitrev);

    // Display the IFFT results
    printf("IFFT result:\n");
    for (int i = 0; i < n; i++)
    {
        printf("%10.5f%+10.5fi\n",
               (float)data[2 * i] / 0x800000,
               (float)data[2 * i + 1] / 0x800000);
    }
}
