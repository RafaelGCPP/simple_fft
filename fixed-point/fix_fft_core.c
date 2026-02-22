#include "fix_fft_core.h"
#include <stdio.h>
#include <math.h>
#include "fix_cplx.h"
#include "fix_bitrev.h"

// This code assumes that n is a power of 2 and implements
// a radix-2 DIF FFT algorithm.
// data - input as real and imaginary interleaved.
// n - number of complex samples in the input data.
// twiddle - precomputed twiddle factors.
// ts - twiddle stride, used when computing the real FFT.
// This version uses no internal complex data types, and
// fixed point arithmetic. Numbers are represented as S8.23.
// The output is in bit reversed order, and the input is in normal order.

static void radix_2_dif_fft_fix_fwd(int *const data, int l2n, int *const twiddle, int ts)
{

    int n = 1 << l2n;

    fix_cplx *const cdata = (fix_cplx *)data;
    fix_cplx *const twd = (fix_cplx *)twiddle;

    // Perform the butterfly operations (DIF)
    for (int stride = n >> 1, tw_index = 1; stride > 0; stride >>= 1, tw_index <<= 1)
    {
        int jmax = n - stride;
        for (int j = 0; j < jmax; j += (stride << 1))
        {
            for (int i = 0; i < stride; i++)
            {
                fix_cplx a, b, w;
                int index = j + i;
                int ti = i * tw_index;

                w = twd[ti * ts];

                fix_cplx *const restrict x = &cdata[index];
                fix_cplx *const restrict y = &cdata[index + stride];

                a = *x;
                b = *y;

                // Butterfly for DIF:
                // x = a + b
                // y = (a - b) * W

                *x = fix_cplx_add(a, b);
                *y = fix_cplx_sub(a, b);

                if (ti != 0) // don't multiply by 1
                {
                    if (ti == (n >> 2))
                    {
                        // multiply by -j
                        fix_cplx_times_neg_j(y);
                    }
                    else
                    {
                        *y = fix_cplx_mul(*y, w, 31);
                    }
                }
            }
        }
    }
}

// This code assumes that n is a power of 2 and implements
// a radix-2 DIT IFFT algorithm.
// data - input as real and imaginary interleaved.
// n - number of complex samples in the input data.
// twiddle - precomputed twiddle factors.
// ts - twiddle stride, used when computing the real FFT.
// This version uses no internal complex data types, and
// fixed point arithmetic. Numbers are represented as S8.23.
// The input is in bit reversed order, and the output is in normal order.
static void radix_2_dit_fft_fix_rev(int *const data, int l2n, int *const twiddle, int ts)
{

    int n = 1 << l2n;

    fix_cplx *const cdata = (fix_cplx *)data;
    fix_cplx *const twd = (fix_cplx *)twiddle;

    // Perform the butterfly operations

    for (int k = 0, stride = 1, tw_index = n >> 1; stride < n; k++, stride <<= 1, tw_index >>= 1)
    {
        int jmax = n - stride;
        for (int j = 0; j < jmax; j += (stride << 1))
        {
            for (int i = 0; i < stride; i++)
            {
                fix_cplx a, b, w;
                int index = j + i;
                int ti = i * tw_index;

                w = twd[ti * ts];

                fix_cplx *const restrict x = &cdata[index];
                fix_cplx *const restrict y = &cdata[index + stride];

                a = *x;
                b = *y;

                if (stride == 1)
                {
                    a = fix_cplx_conj(a);
                    b = fix_cplx_conj(b);
                }
                if (ti != 0) // don't multiply by 1
                {
                    if (ti == (n >> 2))
                    {
                        // multiply by -j
                        fix_cplx_times_neg_j(&b);
                    }
                    else
                    {
                        b = fix_cplx_mul(b, w, 31);
                    }
                }

                *x = fix_cplx_add(a, b);
                *y = fix_cplx_sub(a, b);

                fix_cplx_half(x);
                fix_cplx_half(y);
                if (tw_index == 1)
                {
                    *x = fix_cplx_conj(*x);
                    *y = fix_cplx_conj(*y);
                }
            }
        }
    }
}

void radix_2_fft_fix(int *const data, int l2n, int *const twiddle, int *const bitrev, int ts, int direction)
{
    (void)bitrev;
    if (direction == 1)
    {
        radix_2_dif_fft_fix_fwd(data, l2n, twiddle, ts);
        fix_bitrev(data, l2n, bitrev);
    }
    else
    {
        fix_bitrev(data, l2n, bitrev);
        radix_2_dit_fft_fix_rev(data, l2n, twiddle, ts);
    }
}
