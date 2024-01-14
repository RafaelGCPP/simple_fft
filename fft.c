#include <complex.h>
#include <math.h>
#include <stdio.h>
#include "fft.h"
#include "fft_core.h"


void fft(float *data, int n, float *twiddle, int *bitrev)
{
    radix_2_dit_fft_float(data, n, twiddle, bitrev, 1, 1);
}

void ifft(float *data, int n, float *twiddle, int *bitrev)
{
    float complex *cdata = (float complex *)data;

    radix_2_dit_fft_float(data, n, twiddle, bitrev, 1, -1);

}

