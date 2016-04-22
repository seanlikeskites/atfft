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

/* some handy integer tests */
int atfft_is_even (unsigned int x)
{
    return !(x % 2);
}

int atfft_is_odd (unsigned int x)
{
    return x % 2;
}

int atfft_is_power_of_2 (unsigned int x)
{
    return x && !(x & (x - 1));
}

/* normalisation */
void atfft_normalise_real (atfft_sample *data, int size)
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

atfft_sample atfft_abs (atfft_complex x)
{
    return sqrt (ATFFT_REAL (x) * ATFFT_REAL (x) + ATFFT_IMAG (x) * ATFFT_IMAG (x));
}

atfft_sample atfft_arg (atfft_complex x)
{
    return atan2 (ATFFT_IMAG (x), ATFFT_REAL (x));
}

/* conversion functions */
void atfft_real (atfft_complex *in, atfft_sample *out, int size)
{
    int i = 0;

    for (i = 0; i < size; ++i)
    {
        out [i] = ATFFT_REAL (in [i]);
    }
}

void atfft_imag (atfft_complex *in, atfft_sample *out, int size)
{
    int i = 0;

    for (i = 0; i < size; ++i)
    {
        out [i] = ATFFT_IMAG (in [i]);
    }
}

void atfft_real_to_complex (atfft_sample *in, atfft_complex *out, int size)
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

    if (atfft_is_even (size))
        --lastBin;

    for (i = 1; i < lastBin; ++i)
    {
        ATFFT_REAL (out [size - i]) = ATFFT_REAL (in [i]);
        ATFFT_IMAG (out [size - i]) = - ATFFT_IMAG (in [i]);
    }
}

#ifndef ATFFT_TYPE_FLOAT
void atfft_float_to_sample_real (float *in, atfft_sample *out, int size)
{
    int i = 0;

    for (i = 0; i < size; ++i)
    {
        out [i] = in [i];
    }
}

void atfft_sample_to_float_real (atfft_sample *in, float *out, int size)
{
    int i = 0;

    for (i = 0; i < size; ++i)
    {
        out [i] = in [i];
    }
}

void atfft_float_to_sample_complex (atfft_complex_f *in, atfft_complex *out, int size)
{
    int i = 0;

    for (i = 0; i < size; ++i)
    {
        ATFFT_REAL (out [i]) = ATFFT_REAL (in [i]);
        ATFFT_IMAG (out [i]) = ATFFT_IMAG (in [i]);
    }
}

void atfft_sample_to_float_complex (atfft_complex *in, atfft_complex_f *out, int size)
{
    int i = 0;

    for (i = 0; i < size; ++i)
    {
        ATFFT_REAL (out [i]) = ATFFT_REAL (in [i]);
        ATFFT_IMAG (out [i]) = ATFFT_IMAG (in [i]);
    }
}
#endif

#ifndef ATFFT_TYPE_DOUBLE
void atfft_double_to_sample_real (double *in, atfft_sample *out, int size)
{
    int i = 0;

    for (i = 0; i < size; ++i)
    {
        out [i] = in [i];
    }
}

void atfft_sample_to_double_real (atfft_sample *in, double *out, int size)
{
    int i = 0;

    for (i = 0; i < size; ++i)
    {
        out [i] = in [i];
    }
}

void atfft_double_to_sample_complex (atfft_complex_d *in, atfft_complex *out, int size)
{
    int i = 0;

    for (i = 0; i < size; ++i)
    {
        ATFFT_REAL (out [i]) = ATFFT_REAL (in [i]);
        ATFFT_IMAG (out [i]) = ATFFT_IMAG (in [i]);
    }
}

void atfft_sample_to_double_complex (atfft_complex *in, atfft_complex_d *out, int size)
{
    int i = 0;

    for (i = 0; i < size; ++i)
    {
        ATFFT_REAL (out [i]) = ATFFT_REAL (in [i]);
        ATFFT_IMAG (out [i]) = ATFFT_IMAG (in [i]);
    }
}
#endif

#ifndef ATFFT_TYPE_LONG_DOUBLE
void atfft_long_double_to_sample_real (long double *in, atfft_sample *out, int size)
{
    int i = 0;

    for (i = 0; i < size; ++i)
    {
        out [i] = in [i];
    }
}

void atfft_sample_to_long_double_real (atfft_sample *in, long double *out, int size)
{
    int i = 0;

    for (i = 0; i < size; ++i)
    {
        out [i] = in [i];
    }
}

void atfft_long_double_to_sample_complex (atfft_complex_l *in, atfft_complex *out, int size)
{
    int i = 0;

    for (i = 0; i < size; ++i)
    {
        ATFFT_REAL (out [i]) = ATFFT_REAL (in [i]);
        ATFFT_IMAG (out [i]) = ATFFT_IMAG (in [i]);
    }
}

void atfft_sample_to_long_double_complex (atfft_complex *in, atfft_complex_l *out, int size)
{
    int i = 0;

    for (i = 0; i < size; ++i)
    {
        ATFFT_REAL (out [i]) = ATFFT_REAL (in [i]);
        ATFFT_IMAG (out [i]) = ATFFT_IMAG (in [i]);
    }
}
#endif
