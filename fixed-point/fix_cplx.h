#pragma once
#if defined(__ARM_ACLE)
#include <arm_acle.h>
#endif

// This is a fixed-point implementation of complex numbers.
// Fixed point numbers will be S8.23, i.e. 8 bits for the integer part and 23 bits for the fractional part.

typedef struct s_fix_complex
{
    int real;
    int imag;
} fix_cplx;

// Fixed-point multiply: (x * y) >> N, using 64-bit intermediate to avoid overflow.
// The half-LSB bias (1 << (N-1)) is added before the shift to perform
// round-to-nearest instead of truncation, avoiding accumulated negative bias.
static inline int fix_mul(int x, int y, int N)
{
#if defined(__ARM_FEATURE_DSP)
    // __smull gives a true 64-bit result in a single instruction
    return (int)(__smull(x, y) >> N);
#else
    return (int)(((long long)x * (long long)y + (1LL << (N - 1))) >> N);
#endif
}

// Halve a complex number in place.
static inline void fix_cplx_half(fix_cplx *x)
{
    x->real >>= 1;
    x->imag >>= 1;
}

// Return x + y.
static inline fix_cplx fix_cplx_add(fix_cplx x, fix_cplx y)
{
#if defined(__ARM_FEATURE_DSP)
    return (fix_cplx){ __qadd(x.real, y.real), __qadd(x.imag, y.imag) };
#else
    return (fix_cplx){ x.real + y.real, x.imag + y.imag };
#endif
}

// Return x - y.
static inline fix_cplx fix_cplx_sub(fix_cplx x, fix_cplx y)
{
#if defined(__ARM_FEATURE_DSP)
    return (fix_cplx){ __qsub(x.real, y.real), __qsub(x.imag, y.imag) };
#else
    return (fix_cplx){x.real - y.real, x.imag - y.imag};
#endif
}

// Return x * y using fixed-point fractional bits N.
static inline fix_cplx fix_cplx_mul(fix_cplx x, fix_cplx y, int N)
{
    return (fix_cplx){
        fix_mul(x.real, y.real, N) - fix_mul(x.imag, y.imag, N),
        fix_mul(x.real, y.imag, N) + fix_mul(x.imag, y.real, N)};
}

// Return the complex conjugate of x.
static inline fix_cplx fix_cplx_conj(fix_cplx x)
{
    return (fix_cplx){x.real, -x.imag};
}

// Multiply x by j (i.e. rotate 90°) in place: x = j * x => (-imag, real).
static inline void fix_cplx_times_j(fix_cplx *x)
{
    int temp = x->real;
    x->real = -x->imag;
    x->imag = temp;
}

// Multiply x by -j (i.e. rotate -90°) in place: x = -j * x => (imag, -real).
static inline void fix_cplx_times_neg_j(fix_cplx *x)
{
    int temp = x->real;
    x->real = x->imag;
    x->imag = -temp;
}


#ifdef __ARM_FEATURE_DSP
static inline uint32_t rbit32(uint32_t x)
{
    uint32_t y;
    __asm__ volatile ("rbit %0, %1" : "=r"(y) : "r"(x));
    return y;
}
#endif