#include <complex.h>
#include <math.h>
#include <stdio.h>
#include "fix_fft.h"
#include "fix_fft_core.h"
#include <stdint.h>

void fft_fix(int32_t* const data, int l2n, int32_t* const twiddle, int* const bitrev)
{
    radix_2_fft_fix(data, l2n, twiddle, bitrev, 1, 1);
}

void ifft_fix(int32_t* const data, int l2n, int32_t* const twiddle, int* const bitrev)
{
    radix_2_fft_fix(data, l2n, twiddle, bitrev, 1, -1);
}
