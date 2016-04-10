/*
 * Copyright (C) 2016 Sean Enderby <sean.enderby@gmail.com>
 *
 * This work is free. You can redistribute it and/or modify it under the
 * terms of the Do What The Fuck You Want To Public License, Version 2,
 * as published by Sam Hocevar. See the COPYING file for more details.
 */

#include <math.h>
#include <string.h>
#include <atfft.h>

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

void atfft_real (struct atfft_complex *in, double *out, int size)
{
    int i = 0;

    for (i = 0; i < size; ++i)
    {
        out [i] = in [i].real;
    }
}

void atfft_imag (struct atfft_complex *in, double *out, int size)
{
    int i = 0;

    for (i = 0; i < size; ++i)
    {
        out [i] = in [i].imag;
    }
}

void atfft_real_to_complex (double *in, struct atfft_complex *out, int size)
{
    int i = 0;

    for (i = 0; i < size; ++i)
    {
        out [i].real = in [i];
        out [i].imag = 0;
    }
}

void atfft_halfcomplex_to_complex (struct atfft_complex *in, struct atfft_complex *out, int size)
{
    int i = 0;
    int lastBin = size / 2 + 1;

    memcpy (out, in, lastBin * sizeof (*out));

    if (isEven (size))
    {
        --lastBin;
    }

    for (i = 1; i < lastBin; ++i)
    {
        out [size - i].real = in [i].real;
        out [size - i].imag = -in [i].imag;
    }
}

void atfft_normalise_real (double *data, int size)
{
    int i = 0;

    for (i = 0; i < size; ++i)
    {
        data [i] /= size;
    }
}

void atfft_normalise_complex (struct atfft_complex *data, int size)
{
    int i = 0;

    for (i = 0; i < size; ++i)
    {
        data [i].real /= size;
        data [i].imag /= size;
    }
}
