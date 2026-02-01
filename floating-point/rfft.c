#include <stdio.h>
#include <math.h>
#include "fft_core.h"
#include "cplx_float.h"

// This code assumes that n is a power of 2 and implements
// the real FFT algorithm. It is based on the fact that
// the FFT of a real signal is conjugate symmetric. Then you can
// consider the real input as interleaved complex data and perform
// a complex FFT. The even and odd parts are then computed.
// data - input as real data.
// n - number of real samples in the input data.
// twiddle - precomputed twiddle factors
// the result is stored in the input array, as half of the spectrum
// the N/2 sample is stored on the imaginary part of the first sample.
// Reference:
// Numerical Recipes in C, chapter 12.3
// this code uses no internal complex data types
static void rfft_core_forward(float *const data, int n, float *const twiddle, int *const bitrev)
{

    cplx *cdata = (cplx *)data;
    cplx *twd = (cplx *)twiddle;

    cplx even;
    cplx odd;

    radix_2_dit_fft(data, n / 2, twiddle, bitrev, 2, 1); // the twiddle stride is 2 so we can reuse them when computing the real FFT

    cplx tmp = cdata[0];
    cdata[0].real = tmp.real + tmp.imag;
    cdata[0].imag = tmp.real - tmp.imag;
    cplx_conj(cdata[n / 4], cdata[n / 4]);

    for (int i = 1; i < n / 4; i++)
    {
        // Symmetries!
        // even = (cdata[i] + conj(cdata[n / 2 - i])) / 2;
        // odd = (cdata[i] - conj(cdata[n / 2 - i])) / 2;
        // even1 = (cdata[n / 2 - i] + conj(cdata[i])) / 2; == conj(even1)
        // odd1 = (cdata[n / 2 - i] - conj(cdata[i])) / 2;  == -conj(odd1)

        cplx *const restrict x = &cdata[i];
        cplx *const restrict y = &cdata[n / 2 - i];

        cplx tmp;
        cplx_conj(tmp, *y); // tmp = conj(cdata[i])
        cplx_add(even, *x, tmp);
        cplx_half(even);
        cplx_sub(odd, *x, tmp);
        cplx_half(odd);

        // In forward transform
        // w  =  twd[i]
        // w1 =  twd[n/2 - i] == -conj(w)
        // cdata[i] = even - I * odd * w;
        // cdata[n / 2 - i] = even1 - I * odd1 * w1;
        //                  = conj(even) - i * conj (odd * w)
        // In inverse transform
        // w= conj(twd[i])
        // cdata[i] = even + I * odd * w;
        // cdata[n / 2 - i] = even1 + I * odd1 * w1;
        //                  = conj(even) + i * conj (odd * w)

        // we will build I * odd * w as tmp
        cplx_mul(tmp, odd, twd[i]);
        cplx_times_j(tmp); // tmp = I * odd * w

        // cdata[i] = even - I * odd * w;
        cplx_sub(*x, even, tmp);

        // tmp= I * odd1 * w1 = I * conj (tmp1)
        tmp.real = -tmp.real;

        // cdata[n / 2 - i] = even1 - I * odd1 * w1;
        // = conj(even) - i * conj (odd * w)
        cplx_conj(even, even);
        cplx_sub(*y, even, tmp);
    }
}

static void rfft_core_inverse(float *const data, int n, float *const twiddle, int *const bitrev)
{

    cplx *cdata = (cplx *)data;
    cplx *twd = (cplx *)twiddle;

    cplx even;
    cplx odd;

    cplx tmp;
    tmp.real = (cdata[0].real + cdata[0].imag) / 2;
    tmp.imag = (cdata[0].real - cdata[0].imag) / 2;
    cdata[0] = tmp;

    cplx_conj(cdata[n / 4], cdata[n / 4]);

    for (int i = 1; i < n / 4; i++)
    {
        // Symmetries!
        // even = (cdata[i] + conj(cdata[n / 2 - i])) / 2;
        // odd = (cdata[i] - conj(cdata[n / 2 - i])) / 2;
        // even1 = (cdata[n / 2 - i] + conj(cdata[i])) / 2; == conj(even1)
        // odd1 = (cdata[n / 2 - i] - conj(cdata[i])) / 2;  == -conj(odd1)

        cplx *const restrict x = &cdata[i];
        cplx *const restrict y = &cdata[n / 2 - i];

        cplx tmp;
        cplx_conj(tmp, *y); // tmp = conj(cdata[i])
        cplx_add(even, *x, tmp);
        cplx_half(even);
        cplx_sub(odd, *x, tmp);
        cplx_half(odd);

        // In forward transform
        // w  =  twd[i]
        // w1 =  twd[n/2 - i] == -conj(w)
        // cdata[i] = even - I * odd * w;
        // cdata[n / 2 - i] = even1 - I * odd1 * w1;
        //                  = conj(even) - i * conj (odd * w)
        // In inverse transform
        // w= conj(twd[i])
        // cdata[i] = even + I * odd * w;
        // cdata[n / 2 - i] = even1 + I * odd1 * w1;
        //                  = conj(even) + i * conj (odd * w)

        // we will build I * odd * w as tmp
        cplx w;
        cplx_conj(w, twd[i]);
        cplx_mul(tmp, odd, w);
        cplx_times_j(tmp); // tmp = I * odd * w
        // cdata[i] = even + I * odd * w;
        cplx_add(*x, even, tmp);

        // cdata[n / 2 - i] = even1 + I * odd1 * w1;
        // = conj(even) + i * conj (odd * w)
        cplx_conj(even, even);
        tmp.real = -tmp.real;
        cplx_add(*y, even, tmp);
    }
    radix_2_dit_fft(data, n / 2, twiddle, bitrev, 2, -1);
    
}

void rfft(float *const data, int n, float *const twiddle, int *const bitrev)
{
    rfft_core_forward(data, n, twiddle, bitrev);
}

void irfft(float *const data, int n, float *const twiddle, int *const bitrev)
{

    rfft_core_inverse(data, n, twiddle, bitrev);
}
