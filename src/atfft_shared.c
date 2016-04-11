/*
 * Copyright (C) 2016 Sean Enderby <sean.enderby@gmail.com>
 *
 * This program is free software. It comes without any warranty, to
 * the extent permitted by applicable law. You can redistribute it 
 * and/or modify it under the terms of the Do What The Fuck You Want 
 * To Public License, Version 2, as published by Sam Hocevar. See 
 * the COPYING file for more details.
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

void atfft_real (atfft_complex *in, double *out, int size)
{
    int i = 0;

    for (i = 0; i < size; ++i)
    {
        out [i] = ATFFT_REAL (in [i]);
    }
}

void atfft_imag (atfft_complex *in, double *out, int size)
{
    int i = 0;

    for (i = 0; i < size; ++i)
    {
        out [i] = ATFFT_IMAG (in [i]);
    }
}

void atfft_real_to_complex (double *in, atfft_complex *out, int size)
{
    int i = 0;

    for (i = 0; i < size; ++i)
    {
        ATFFT_REAL (out [i]) = in [i];
        ATFFT_IMAG (out [i]) = 0;
    }
}

void atfft_halfcomplex_to_complex (atfft_complex *in, atfft_complex *out, int size)
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
        ATFFT_REAL (out [size - i]) = ATFFT_REAL (in [i]);
        ATFFT_IMAG (out [size - i]) = - ATFFT_IMAG (in [i]);
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

void atfft_normalise_complex (atfft_complex *data, int size)
{
    int i = 0;

    for (i = 0; i < size; ++i)
    {
        ATFFT_REAL (data [i]) /= size;
        ATFFT_IMAG (data [i]) /= size;
    }
}
