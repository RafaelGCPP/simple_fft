#pragma once

void rfft(float* const data, int n, float* const twiddle, int* const bitrev);
void irfft(float* const data, int n, float* const twiddle, int* const bitrev);

void precompute_twiddle_factors_rfft(float* const twiddle, int n);
void precompute_bitrev_table_rfft(int* const bitrev, int n);