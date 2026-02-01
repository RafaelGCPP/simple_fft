#include <complex.h>
#include <math.h>
#include <stdio.h>
#include "fix_fft.h"
#include "fix_fft_core.h"

void fft_fix(int* const data, int l2n, int* const twiddle, int* const bitrev)
{
    radix_2_dit_fft_fix(data, l2n, twiddle, bitrev, 1, 1);
}

void ifft_fix(int* const data, int l2n, int* const twiddle, int* const bitrev)
{
    radix_2_dit_fft_fix(data, l2n, twiddle, bitrev, 1, -1);
}
