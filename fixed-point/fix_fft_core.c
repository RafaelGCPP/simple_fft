#include "fix_fft_core.h"
#include <stdio.h>
#include <math.h>
#include "fix_cplx.h"

// This code assumes that n is a power of 2 and implements
// a radix-2 DIT FFT algorithm.
// data - input as real and imaginary interleaved.
// n - number of complex samples in the input data.
// twiddle - precomputed twiddle factors.
// ts - twiddle stride, used when computing the real FFT.
// This version uses no internal complex data types, and
// fixed point arithmetic. Numbers are represented as S8.23.
void radix_2_dit_fft_fix(int *data, int l2n, int *twiddle, int *bitrev, int ts, int direction)
{

    int n = 1 << l2n;

    fix_cplx *cdata = (fix_cplx *)data;
    fix_cplx *twd = (fix_cplx *)twiddle;

    // Bit-reverse the input data

    for (int i = 1; i < n - 1; i++)
    {
        int j = bitrev[i];
        if (i < j)
        {
            fix_cplx tmp = cdata[i];
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
                fix_cplx a, b, w;
                int index = j + i;
                int ti = i * tw_index;

                w=twd[ti * ts];
                a=cdata[index];
                b=cdata[index + stride];

                if (direction == -1 && stride == 1)
                {
                    fix_cplx_conj(a, a);
                    fix_cplx_conj(b, b);
                }
                if (ti != 0) // don't multiply by 1
                {
                    if (ti == (n >> 2))
                    {
                        // multiply by -j
                        fix_cplx_times_neg_j(b);
                    }
                    else
                    {
                        fix_cplx r;
                        fix_cplx_mul(r, b, w, 31);
                        b = r;
                    }
                }

                fix_cplx_add(cdata[index], a, b);
                fix_cplx_sub(cdata[index + stride], a, b);

                if (direction == -1)
                {
                    fix_cplx_half(cdata[index]);
                    fix_cplx_half(cdata[index + stride]);
                    if (tw_index == 1)
                    {
                        fix_cplx_conj(cdata[index], cdata[index]);
                        fix_cplx_conj(cdata[index + stride], cdata[index + stride]);
                    }
                }
            }
        }
    }
}
