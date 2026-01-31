#include "fft_core.h"
#include <stdio.h>
#include <math.h>
#include "cplx_float.h"

// This code assumes that n is a power of 2 and implements
// a radix-2 DIT FFT algorithm.
// data - input as real and imaginary interleaved.
// n - number of complex samples in the input data.
// twiddle - precomputed twiddle factors.
// ts - twiddle stride, used when computing the real FFT.
// This version uses no internal complex data types.
void radix_2_dit_fft(float *data, int n, float *twiddle, int *bitrev, int ts, int direction)
{

    cplx *cdata = (cplx *)data;
    cplx *twd = (cplx *)twiddle;

    // Bit-reverse the input data

    for (int i = 1; i < n - 1; i++)
    {
        int j = bitrev[i];
        if (i < j)
        {
            cplx tmp = cdata[i];
            cdata[i] = cdata[j];
            cdata[j] = tmp;
        }
    }

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

                w=twd[i * tw_index * ts];

                a=cdata[index];
                b=cdata[index + stride];

                if (direction == -1 && stride == 1)
                {
                    cplx_conj(a, a);
                    cplx_conj(b, b);
                }
                cplx r;
                cplx_mul(r, b, w);
                b = r;

                cplx_add(cdata[index], a, b);
                cplx_sub(cdata[index + stride], a, b);

                if (direction == -1)
                {
                    cplx_half(cdata[index]);
                    cplx_half(cdata[index + stride]);
                    if (tw_index == 1)
                    {
                        cplx_conj(cdata[index], cdata[index]);
                        cplx_conj(cdata[index + stride], cdata[index + stride]);
                    }
                }
            }
        }
    }
}
