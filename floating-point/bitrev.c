#include "bitrev.h"
#include "cplx_float.h"

#if defined(__ARM_ACLE)
#include <arm_acle.h>
#endif

// Perform bit-reversal permutation on n complex samples in data.
// l2n  - log2(n), used to shift the RBIT result on ARM.
// bitrev_table - precomputed bit-reversal table, or NULL if ARM RBIT
//                intrinsic is available (in which case no table is needed).
void bitrev(float *const data, int n, int *const bitrev_table)
{
    cplx *const cdata = (cplx *)data;

#if defined(__ARM_ACLE) && defined(__ARM_FEATURE_BIT_REVERSAL)
    // Single-cycle RBIT instruction available: table is not needed.
    (void)bitrev_table; // suppress unused-parameter warning
    for (int i = 1; i < n - 1; i++)
    {
        int j = (int)(__rbit((uint32_t)i) >> (32 - __builtin_ctz(n)));
        if (i < j)
        {
            cplx tmp = cdata[i];
            cdata[i] = cdata[j];
            cdata[j] = tmp;
        }
    }
#else
    // Fallback: precomputed table lookup. Table must not be NULL.
    for (int i = 1; i < n - 1; i++)
    {
        int j = bitrev_table[i];
        if (i < j)
        {
            cplx tmp = cdata[i];
            cdata[i] = cdata[j];
            cdata[j] = tmp;
        }
    }
#endif
}