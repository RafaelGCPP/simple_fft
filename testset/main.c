
#include "test_set.h"
#include "benchmark.h"

int main()
{

    cordic_test();
    twiddle_test();
    fix_twiddle_test();
    complex_fft_test();
    real_fft_test();
    complex_fft_test_fix();
    real_fft_test_fix();
    fft_benchmark();
    fixed_fft_benchmark();

    return 0;
}
