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
void rfft_core(float *data, int n, float *twiddle, int *bitrev, int direction)
{

    cplx *cdata = (cplx *)data;
    cplx *twd = (cplx *)twiddle;

    cplx even;
    cplx odd;

    if (direction == 1)
    {
        radix_2_dit_fft(data, n / 2, twiddle, bitrev, 2, 1); // the twiddle stride is 2 so we can reuse them when computing the real FFT

        cplx tmp=cdata[0];
        cdata[0].real = tmp.real + tmp.imag;
        cdata[0].imag = tmp.real - tmp.imag;    
    }
    else // (direction == -1)
    {
        // even = (creal(cdata[0]) + cimag(cdata[0])) / 2
        // odd = (creal(cdata[0]) - cimag(cdata[0])) / 2;
        // cdata[0] = even + I * odd; // since the input is real, the imaginary part of each is zero so

        cplx tmp;
        tmp.real = (cdata[0].real + cdata[0].imag) / 2;
        tmp.imag = (cdata[0].real - cdata[0].imag) / 2;        
        cdata[0]=tmp;
    }

    //if (direction == 1) {
        cplx tmp=cdata[n / 4];
        cplx_conj(cdata[n / 4], tmp);
    //} 

    for (int i = 1; i < n / 4; i++)
    {
        // Symmetries!
        // even = (cdata[i] + conj(cdata[n / 2 - i])) / 2;
        // odd = (cdata[i] - conj(cdata[n / 2 - i])) / 2;
        // even1 = (cdata[n / 2 - i] + conj(cdata[i])) / 2; == conj(even1)
        // odd1 = (cdata[n / 2 - i] - conj(cdata[i])) / 2;  == -conj(odd1)


        cplx tmp;
        cplx_conj(tmp, cdata[n/2 - i]); // tmp = conj(cdata[i])
        cplx_add(even, cdata[i], tmp);
        cplx_half(even);
        cplx_sub(odd, cdata[i], tmp);
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



        if (direction == 1)
        {
            // we will build I * odd * w as tmp
            cplx_mul(tmp, odd, twd[i]);
            cplx_times_j(tmp); // tmp = I * odd * w

            // cdata[i] = even - I * odd * w;
            cplx_sub(cdata[i], even, tmp);

            // tmp= I * odd1 * w1 = I * conj (tmp1)
            tmp.real = -tmp.real;
            
            // cdata[n / 2 - i] = even1 - I * odd1 * w1;
            // = conj(even) - i * conj (odd * w)
            cplx_conj(even, even);
            cplx_sub(cdata[n / 2 - i], even, tmp);
        }
        else
        {

            // we will build I * odd * w as tmp
            cplx w;
            cplx_conj(w, twd[i]);
            cplx_mul(tmp, odd, w);
            cplx_times_j(tmp); // tmp = I * odd * w
            // cdata[i] = even + I * odd * w;
            cplx_add(cdata[i], even, tmp);
           
            // cdata[n / 2 - i] = even1 + I * odd1 * w1;
            // = conj(even) + i * conj (odd * w)
            cplx_conj(even, even);
            tmp.real = -tmp.real;   
            cplx_add(cdata[n / 2 - i], even, tmp);

        }
    }
    if (direction == -1)
    {
        radix_2_dit_fft(data, n / 2, twiddle, bitrev, 2, -1);
    }
}

void rfft(float *data, int n, float *twiddle, int *bitrev)
{
    rfft_core(data, n, twiddle, bitrev, 1);
}

void irfft(float *data, int n, float *twiddle, int *bitrev)
{

    rfft_core(data, n, twiddle, bitrev, -1);
}
