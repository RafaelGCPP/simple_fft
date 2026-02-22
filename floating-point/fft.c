#include <complex.h>
#include <math.h>
#include <stdio.h>
#include "fft.h"
#include "fft_core.h"


void fft(float* const data, int n, float* const twiddle, int* const bitrev)
{
    radix_2_fft(data, n, twiddle, bitrev, 1, 1);
}

void ifft(float* const data, int n, float* const twiddle, int* const bitrev)
{
    radix_2_fft(data, n, twiddle, bitrev, 1, -1);
}
