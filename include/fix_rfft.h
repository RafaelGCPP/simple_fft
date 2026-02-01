#pragma once

void rfft_fix(int* const data, int l2n, int* const twiddle, int* const bitrev);
void irfft_fix(int* const data, int l2n, int* const twiddle, int* const bitrev);

void precompute_twiddle_factors_rfft_fix(int* const twiddle, int l2n);
void precompute_bitrev_table_rfft_fix(int* const bitrev, int l2n);