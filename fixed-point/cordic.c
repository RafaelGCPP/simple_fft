#include "cordic.h"
#include <math.h>
#include <stdint.h>

#define ITERATIONS 31
#define K 0x4dba76d4 // 0.6072529350088812561694

// atan(2^-i) / (pi/2) * 2^31
const int32_t atan_table[] = {
    0x40000000, 0x25c80a3b, 0x13f670b6, 0x0a2223a8,
    0x05161a86, 0x028bafc2, 0x0145ec3c, 0x00a2f8aa,
    0x00517ca6, 0x0028be5d, 0x00145f30, 0x000a2f98,
    0x000517cc, 0x00028be6, 0x000145f3, 0x0000a2f9,
    0x0000517c, 0x000028be, 0x0000145f, 0x00000a2f,
    0x00000517, 0x0000028b, 0x00000145, 0x000000a2,
    0x00000051, 0x00000028, 0x00000014, 0x0000000a,
    0x00000005, 0x00000002, 0x00000001};

// theta is scaled so that [-1,1] maps to [-pi/2, pi/2]
// to use this function one must reduce the angle to the first or fourth quadrants
// and then scale it to the range [-1,1]
// 31 bits fractional representation means that 1 is represented as 0x7fffffff and -1 as 0x80000000
void cordic_sin_cos(int32_t theta, int32_t *sin, int32_t *cos)
{

    switch (theta)
    {
    case 0x00000000:
        *sin = 0;
        *cos = 0x7fffffff;
        return;
    case 0x7fffffff:
        *sin = 0x7fffffff;
        *cos = 0;
        return;
    case 0x80000000:
        *sin = 0x80000000;
        *cos = 0;
        return;
    }

    int32_t x = K;
    int32_t y = 0;
    int32_t z = theta;

    for (int i = 0; i < ITERATIONS; i++)
    {
        int32_t tx, ty;

        if (z >= 0)
        {
            tx = x - (y >> i);
            ty = y + (x >> i);
            z -= atan_table[i];
        }
        else
        {
            tx = x + (y >> i);
            ty = y - (x >> i);
            z += atan_table[i];
        }

        x = tx;
        y = ty;
    }

    *sin = y;
    *cos = x;
}
