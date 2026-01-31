#pragma once

// This is a fixed-point implementation of complex numbers.
// Fixed point numbers will be S8.23, i.e. 8 bits for the integer part and 23 bits for the fractional part.

typedef struct s_fix_complex
{
    int real;
    int imag;
} fix_cplx;

#define fix_mul(x, y, N) ((int)((long long)x) * ((long long)y) >> (N))

#define fix_cplx_half(x) (    \
    (x).real >>= 1, \
    (x).imag >>= 1)

#define fix_cplx_add(r, x, y) (     \
    (r).real = (x).real + (y).real, \
    (r).imag = (x).imag + (y).imag)

#define fix_cplx_sub(r, x, y) (     \
    (r).real = (x).real - (y).real, \
    (r).imag = (x).imag - (y).imag)

#define fix_cplx_mul(r, x, y, N) (                                                  \
    (r).real = fix_mul((x).real, (y).real, (N)) - fix_mul((x).imag, (y).imag, (N)), \
    (r).imag = fix_mul((x).real, (y).imag, (N)) + fix_mul((x).imag, (y).real, (N)))

#define fix_cplx_conj(r, x) ( \
    (r).real = (x).real,      \
    (r).imag = -(x).imag)

#define fix_cplx_times_j(x) {  int temp = (x).real;  (x).real = -(x).imag; (x).imag = temp;  }
#define fix_cplx_times_neg_j(x) {  int temp = (x).real;  (x).real = (x).imag; (x).imag = -temp;  }

