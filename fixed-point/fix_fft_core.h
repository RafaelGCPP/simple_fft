#pragma once
#include <stdint.h>

void radix_2_fft_fix(int32_t* const data, int l2n, int32_t* const twiddle, int* const bitrev, int ts, int direction);
