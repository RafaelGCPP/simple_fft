#pragma once

// Although most compilers implement complex numbers as a struct with two floats,
// we can't rely on this on embedded systems.  So we define our own complex type.

// Also, it will make it easier to port to fixed-point if we ever need to.

typedef struct s_complex_float
{
    float real;
    float imag;
} cplx;

// Scale x by scalar s in place.
static inline void cplx_scale(cplx *x, float s)
{
    x->real *= s;
    x->imag *= s;
}

// Halve x in place.
static inline void cplx_half(cplx *x)
{
    cplx_scale(x, 0.5f);
}

// Return x + y.
static inline cplx cplx_add(cplx x, cplx y)
{
    return (cplx){ x.real + y.real, x.imag + y.imag };
}

// Return x - y.
static inline cplx cplx_sub(cplx x, cplx y)
{
    return (cplx){ x.real - y.real, x.imag - y.imag };
}

// Return x * y.
static inline cplx cplx_mul(cplx x, cplx y)
{
    return (cplx){
        x.real * y.real - x.imag * y.imag,
        x.real * y.imag + x.imag * y.real
    };
}

// Return the complex conjugate of x.
static inline cplx cplx_conj(cplx x)
{
    return (cplx){ x.real, -x.imag };
}

// Multiply x by j (rotate 90°) in place: x = j * x => (-imag, real).
static inline void cplx_times_j(cplx *x)
{
    float temp = x->real;
    x->real    = -x->imag;
    x->imag    = temp;
}

// Multiply x by -j (rotate -90°) in place: x = -j * x => (imag, -real).
static inline void cplx_times_neg_j(cplx *x)
{
    float temp = x->real;
    x->real    = x->imag;
    x->imag    = -temp;
}
