/*
 * Copyright (C) 2020 Sean Enderby <sean.enderby@gmail.com>
 *
 * This program is free software. It comes without any warranty, to
 * the extent permitted by applicable law. You can redistribute it
 * and/or modify it under the terms of the Do What The Fuck You Want
 * To Public License, Version 2, as published by Sam Hocevar. See
 * the COPYING file for more details.
 */

#include <math.h>
#include <string.h>
#include <atfft/atfft_shared.h>

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
void atfft_scale_real (atfft_sample *data, int size, atfft_sample scale_factor)
{
    for (int i = 0; i < size; ++i)
    {
        data [i] *= scale_factor;
    }
}

void atfft_normalise_real (atfft_sample *data, int size)
{
    atfft_scale_real (data, size, 1.0 / size);
}

void atfft_scale_complex (atfft_complex *data, int size, atfft_sample scale_factor)
{
    for (int i = 0; i < size; ++i)
    {
        ATFFT_RE (data [i]) *= scale_factor;
        ATFFT_IM (data [i]) *= scale_factor;
    }
}

void atfft_normalise_complex (atfft_complex *data, int size)
{
    atfft_scale_complex (data, size, 1.0 / size);
}

atfft_sample atfft_abs (const atfft_complex x)
{
    return sqrt (ATFFT_RE (x) * ATFFT_RE (x) + ATFFT_IM (x) * ATFFT_IM (x));
}

atfft_sample atfft_arg (const atfft_complex x)
{
    return atan2 (ATFFT_IM (x), ATFFT_RE (x));
}

/* conversion functions */
void atfft_real (atfft_complex *in, atfft_sample *out, int size)
{
    for (int i = 0; i < size; ++i)
    {
        out [i] = ATFFT_RE (in [i]);
    }
}

void atfft_real_stride (atfft_complex *in,
                        int in_stride,
                        atfft_sample *out,
                        int out_stride,
                        int size)
{
    for (int i = 0, o = 0; i < size * in_stride; i += in_stride, o += out_stride)
    {
        out [o] = ATFFT_RE (in [i]);
    }
}

void atfft_imag (atfft_complex *in, atfft_sample *out, int size)
{
    for (int i = 0; i < size; ++i)
    {
        out [i] = ATFFT_IM (in [i]);
    }
}

void atfft_imag_stride (atfft_complex *in,
                        int in_stride,
                        atfft_sample *out,
                        int out_stride,
                        int size)
{
    for (int i = 0, o = 0; i < size * in_stride; i += in_stride, o += out_stride)
    {
        out [o] = ATFFT_IM (in [i]);
    }
}

void atfft_real_to_complex (const atfft_sample *in, atfft_complex *out, int size)
{
    for (int i = 0; i < size; ++i)
    {
        ATFFT_RE (out [i]) = in [i];
        ATFFT_IM (out [i]) = 0;
    }
}

void atfft_real_to_complex_stride (const atfft_sample *in,
                                   int in_stride,
                                   atfft_complex *out,
                                   int out_stride,
                                   int size)
{
    for (int i = 0, o = 0; i < size * in_stride; i += in_stride, o += out_stride)
    {
        ATFFT_RE (out [o]) = in [i];
        ATFFT_IM (out [o]) = 0;
    }
}

void atfft_float_to_sample_real_stride (const float *in,
                                        int in_stride,
                                        atfft_sample *out,
                                        int out_stride,
                                        int size)
{
    for (int i = 0, o = 0; i < size * in_stride; i += in_stride, o += out_stride)
    {
        out [o] = in [i];
    }
}

void atfft_sample_to_float_real (const atfft_sample *in, float *out, int size)
{
    for (int i = 0; i < size; ++i)
    {
        out [i] = in [i];
    }
}

void atfft_sample_to_float_real_stride (const atfft_sample *in,
                                        int in_stride,
                                        float *out,
                                        int out_stride,
                                        int size)
{
    for (int i = 0, o = 0; i < size * in_stride; i += in_stride, o += out_stride)
    {
        out [o] = in [i];
    }
}

void atfft_float_to_sample_complex (atfft_complex_f *in, atfft_complex *out, int size)
{
    for (int i = 0; i < size; ++i)
    {
        ATFFT_RE (out [i]) = ATFFT_RE (in [i]);
        ATFFT_IM (out [i]) = ATFFT_IM (in [i]);
    }
}

void atfft_float_to_sample_complex_stride (atfft_complex_f *in,
                                           int in_stride,
                                           atfft_complex *out,
                                           int out_stride,
                                           int size)
{
    for (int i = 0, o = 0; i < size * in_stride; i += in_stride, o += out_stride)
    {
        ATFFT_RE (out [o]) = ATFFT_RE (in [i]);
        ATFFT_IM (out [o]) = ATFFT_IM (in [i]);
    }
}

void atfft_sample_to_float_complex (atfft_complex *in, atfft_complex_f *out, int size)
{
    for (int i = 0; i < size; ++i)
    {
        ATFFT_RE (out [i]) = ATFFT_RE (in [i]);
        ATFFT_IM (out [i]) = ATFFT_IM (in [i]);
    }
}

void atfft_sample_to_float_complex_stride (atfft_complex *in,
                                           int in_stride,
                                           atfft_complex_f *out,
                                           int out_stride,
                                           int size)
{
    for (int i = 0, o = 0; i < size * in_stride; i += in_stride, o += out_stride)
    {
        ATFFT_RE (out [o]) = ATFFT_RE (in [i]);
        ATFFT_IM (out [o]) = ATFFT_IM (in [i]);
    }
}

void atfft_double_to_sample_real (const double *in, atfft_sample *out, int size)
{
    for (int i = 0; i < size; ++i)
    {
        out [i] = in [i];
    }
}

void atfft_double_to_sample_real_stride (const double *in,
                                         int in_stride,
                                         atfft_sample *out,
                                         int out_stride,
                                         int size)
{
    for (int i = 0, o = 0; i < size * in_stride; i += in_stride, o += out_stride)
    {
        out [o] = in [i];
    }
}

void atfft_sample_to_double_real (const atfft_sample *in, double *out, int size)
{
    for (int i = 0; i < size; ++i)
    {
        out [i] = in [i];
    }
}

void atfft_sample_to_double_real_stride (const atfft_sample *in,
                                         int in_stride,
                                         double *out,
                                         int out_stride,
                                         int size)
{
    for (int i = 0, o = 0; i < size * in_stride; i += in_stride, o += out_stride)
    {
        out [o] = in [i];
    }
}

void atfft_double_to_sample_complex (atfft_complex_d *in, atfft_complex *out, int size)
{
    for (int i = 0; i < size; ++i)
    {
        ATFFT_RE (out [i]) = ATFFT_RE (in [i]);
        ATFFT_IM (out [i]) = ATFFT_IM (in [i]);
    }
}

void atfft_double_to_sample_complex_stride (atfft_complex_d *in,
                                            int in_stride,
                                            atfft_complex *out,
                                            int out_stride,
                                            int size)
{
    for (int i = 0, o = 0; i < size * in_stride; i += in_stride, o += out_stride)
    {
        ATFFT_RE (out [o]) = ATFFT_RE (in [i]);
        ATFFT_IM (out [o]) = ATFFT_IM (in [i]);
    }
}

void atfft_sample_to_double_complex (atfft_complex *in, atfft_complex_d *out, int size)
{
    for (int i = 0; i < size; ++i)
    {
        ATFFT_RE (out [i]) = ATFFT_RE (in [i]);
        ATFFT_IM (out [i]) = ATFFT_IM (in [i]);
    }
}

void atfft_sample_to_double_complex_stride (atfft_complex *in,
                                            int in_stride,
                                            atfft_complex_d *out,
                                            int out_stride,
                                            int size)
{
    for (int i = 0, o = 0; i < size * in_stride; i += in_stride, o += out_stride)
    {
        ATFFT_RE (out [o]) = ATFFT_RE (in [i]);
        ATFFT_IM (out [o]) = ATFFT_IM (in [i]);
    }
}

void atfft_long_double_to_sample_real (const long double *in, atfft_sample *out, int size)
{
    for (int i = 0; i < size; ++i)
    {
        out [i] = in [i];
    }
}

void atfft_long_double_to_sample_real_stride (const long double *in,
                                              int in_stride,
                                              atfft_sample *out,
                                              int out_stride,
                                              int size)
{
    for (int i = 0, o = 0; i < size * in_stride; i += in_stride, o += out_stride)
    {
        out [o] = in [i];
    }
}

void atfft_sample_to_long_double_real (const atfft_sample *in, long double *out, int size)
{
    for (int i = 0; i < size; ++i)
    {
        out [i] = in [i];
    }
}

void atfft_sample_to_long_double_real_stride (const atfft_sample *in,
                                              int in_stride,
                                              long double *out,
                                              int out_stride,
                                              int size)
{
    for (int i = 0, o = 0; i < size * in_stride; i += in_stride, o += out_stride)
    {
        out [o] = in [i];
    }
}

void atfft_long_double_to_sample_complex (atfft_complex_l *in, atfft_complex *out, int size)
{
    for (int i = 0; i < size; ++i)
    {
        ATFFT_RE (out [i]) = ATFFT_RE (in [i]);
        ATFFT_IM (out [i]) = ATFFT_IM (in [i]);
    }
}

void atfft_long_double_to_sample_complex_stride (atfft_complex_l *in,
                                                 int in_stride,
                                                 atfft_complex *out,
                                                 int out_stride,
                                                 int size)
{
    for (int i = 0, o = 0; i < size * in_stride; i += in_stride, o += out_stride)
    {
        ATFFT_RE (out [o]) = ATFFT_RE (in [i]);
        ATFFT_IM (out [o]) = ATFFT_IM (in [i]);
    }
}

void atfft_sample_to_long_double_complex (atfft_complex *in, atfft_complex_l *out, int size)
{
    for (int i = 0; i < size; ++i)
    {
        ATFFT_RE (out [i]) = ATFFT_RE (in [i]);
        ATFFT_IM (out [i]) = ATFFT_IM (in [i]);
    }
}

void atfft_sample_to_long_double_complex_stride (atfft_complex *in,
                                                 int in_stride,
                                                 atfft_complex_l *out,
                                                 int out_stride,
                                                 int size)
{
    for (int i = 0, o = 0; i < size * in_stride; i += in_stride, o += out_stride)
    {
        ATFFT_RE (out [o]) = ATFFT_RE (in [i]);
        ATFFT_IM (out [o]) = ATFFT_IM (in [i]);
    }
}
