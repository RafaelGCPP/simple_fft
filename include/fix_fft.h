#pragma once

void fft_fix(int* const data, int l2n, int* const twiddle, int* const bitrev);
void ifft_fix(int* const data, int l2n, int* const twiddle, int* const bitrev);

void precompute_twiddle_factors_fix(int* const twiddle, int n);
void precompute_bitrev_table_fix(int* const bitrev, int n);