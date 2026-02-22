#include "fix_bitrev.h"
#include "fix_cplx.h"

#if defined(__ARM_ACLE)
#include <arm_acle.h>
#endif

// Perform bit-reversal permutation on 2^l2n complex fixed-point samples.
// l2n  - log2(n).
// bitrev_table - precomputed bit-reversal table, or NULL if ARM RBIT
//                intrinsic is available (in which case no table is needed).
void fix_bitrev(int *const data, int l2n, int *const bitrev_table)
{
    int n = 1 << l2n;
    fix_cplx *const cdata = (fix_cplx *)data;

#if defined(__ARM_ACLE) && defined(__ARM_FEATURE_BIT_REVERSAL)
    // Single-cycle RBIT instruction available: table is not needed.
    (void)bitrev_table; // suppress unused-parameter warning
    for (int i = 1; i < n - 1; i++)
    {
        int j = (int)(__rbit((uint32_t)i) >> (32 - l2n));
        if (i < j)
        {
            fix_cplx tmp = cdata[i];
            cdata[i] = cdata[j];
            cdata[j] = tmp;
        }
    }
#else
    // Fallback: precomputed table lookup. Table must not be NULL.
    (void)l2n;
    for (int i = 1; i < n - 1; i++)
    {
        int j = bitrev_table[i];
        if (i < j)
        {
            fix_cplx tmp = cdata[i];
            cdata[i] = cdata[j];
            cdata[j] = tmp;
        }
    }
#endif
}