#pragma once

void precompute_twiddle_factors(float* const twiddle, int n);
void precompute_twiddle_factors_sqrt(float* const twiddle, int n);
void precompute_bitrev_table(int* const bitrev, int n);

void fft(float* const data, int n, float* const twiddle, int* const bitrev);
void ifft(float* const data, int n, float* const twiddle, int* const bitrev);
