#pragma once
#include <stdint.h>

void fft_fix(int32_t* const data, int l2n, int32_t* const twiddle, int* const bitrev);
void ifft_fix(int32_t* const data, int l2n, int32_t* const twiddle, int* const bitrev);

void precompute_twiddle_factors_fix(int32_t* const twiddle, int n);
void precompute_bitrev_table_fix(int* const bitrev, int n);