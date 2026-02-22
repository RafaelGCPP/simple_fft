#include "fft_core.h"
#include <stdio.h>
#include <math.h>
#include "cplx_float.h"
#include "bitrev.h"

// This code assumes that n is a power of 2 and implements
// a radix-2 DIF FFT algorithm.
// data - input as real and imaginary interleaved.
// n - number of complex samples in the input data.
// twiddle - precomputed twiddle factors.
// ts - twiddle stride, used when computing the real FFT.
// This version uses no internal complex data types.
// The output is in bit reversed order, and the input is in normal order.

static void radix_2_dif_fft_fwd(float *const data, int n, float *const twiddle, int ts)
{

    cplx *cdata = (cplx *)data;
    cplx *twd = (cplx *)twiddle;

    // Perform the butterfly operations

    for (int stride = n >> 1, tw_index = 1; stride > 0; stride >>= 1, tw_index <<= 1)
    {
        int jmax = n - stride;
        for (int j = 0; j < jmax; j += (stride << 1))
        {
            for (int i = 0; i < stride; i++)
            {
                cplx a, b, w;
                int index = j + i;

                w = twd[i * tw_index * ts];

                cplx *const restrict x = cdata + index;
                cplx *const restrict y = cdata + index + stride;

                a = *x;
                b = *y;

                *x = cplx_add(a, b);
                *y = cplx_sub(a, b);

                if (i * tw_index * ts != 0) // Optimization: skip mul by 1 (W^0)
                {
                    *y = cplx_mul(*y, w);
                }
            }
        }
    }
}

static void radix_2_dit_fft_rev(float *const data, int n, float *const twiddle, int ts)
{

    cplx *cdata = (cplx *)data;
    cplx *twd = (cplx *)twiddle;

    // Perform the butterfly operations

    for (int k = 0, stride = 1, tw_index = n >> 1; stride < n; k++, stride <<= 1, tw_index >>= 1)
    {
        int jmax = n - stride;
        for (int j = 0; j < jmax; j += (stride << 1))
        {
            for (int i = 0; i < stride; i++)
            {
                cplx a, b, w;
                int index = j + i;

                w = twd[i * tw_index * ts];

                cplx *const restrict x = cdata + index;
                cplx *const restrict y = cdata + index + stride;

                a = *x;
                b = *y;

                if (stride == 1)
                {
                    a = cplx_conj(a);
                    b = cplx_conj(b);
                }
                b = cplx_mul(b, w);

                *x = cplx_add(a, b);
                *y = cplx_sub(a, b);

                cplx_half(x);
                cplx_half(y);
                if (tw_index == 1)
                {
                    *x = cplx_conj(*x);
                    *y = cplx_conj(*y);
                }
            }
        }
    }
}

void radix_2_fft(float *const data, int n, float *const twiddle, int *const bitrev_table, int ts, int direction)
{
    if (direction == 1)
    {
        radix_2_dif_fft_fwd(data, n, twiddle, ts);
        bitrev(data, n, bitrev_table);
    }
    else
    {
        bitrev(data, n, bitrev_table);
        radix_2_dit_fft_rev(data, n, twiddle, ts);
    }
}