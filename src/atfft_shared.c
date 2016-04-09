#include <math.h>
#include "../include/atfft.h"

int isEven (unsigned int x)
{
    return !(x % 2);
}

int isOdd (unsigned int x)
{
    return x % 2;
}

int isPowerOf2 (unsigned int x)
{
    return x && !(x & (x - 1));
}

void atfft_real (double *in, double *out, int size)
{
    int i = 0;

    for (i = 0; i < size; ++i)
    {
        out [i] = in [2 * i];
    }
}

void atfft_imag (double *in, double *out, int size)
{
    int i = 0;

    for (i = 0; i < size; ++i)
    {
        out [i] = in [2 * i + 1];
    }
}

void atfft_real_to_complex (double *in, double *out, int size)
{
    int i = 0;

    for (i = 0; i < size; ++i)
    {
        out [2 * i] = in [i];
        out [2 * i + 1] = 0;
    }
}

void atfft_halfcomplex_to_complex (double *in, double *out, int size)
{
    int i = 0;

    out [0] = in [0];
    out [1] = in [1];

    for (i = 1; i < size / 2; ++i)
    {
        out [2 * i] = in [2 * i];
        out [2 * i + 1] = in [2 * i + 1];

        out [2 * (size - i)] = in [2 * i];
        out [2 * (size - i) + 1] = -in [2 * i + 1];
    }

    out [2 * i] = in [2 * i];
    out [2 * i + 1] = in [2 * i + 1];

    if (isOdd (size))
    {
        out [2 * (size - i)] = in [2 * i];
        out [2 * (size - i) + 1] = -in [2 * i + 1];
    }
}
