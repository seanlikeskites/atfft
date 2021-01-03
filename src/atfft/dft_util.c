/*
 * Copyright (C) 2020 Sean Enderby <sean.enderby@gmail.com>
 *
 * This program is free software. It comes without any warranty, to
 * the extent permitted by applicable law. You can redistribute it
 * and/or modify it under the terms of the Do What The Fuck You Want
 * To Public License, Version 2, as published by Sam Hocevar. See
 * the COPYING file for more details.
 */

#include <string.h>
#include <atfft/dft_util.h>

int atfft_dft_halfcomplex_size (int size)
{
    return size / 2 + 1;
}

void atfft_halfcomplex_to_complex (atfft_complex *in, atfft_complex *out, int size)
{
    int last_bin = atfft_dft_halfcomplex_size (size);

    memcpy (out, in, last_bin * sizeof (*out));

    if (atfft_is_even (size))
        --last_bin;

    for (int i = 1; i < last_bin; ++i)
    {
        ATFFT_RE (out [size - i]) = ATFFT_RE (in [i]);
        ATFFT_IM (out [size - i]) = - ATFFT_IM (in [i]);
    }
}

void atfft_halfcomplex_to_complex_stride (atfft_complex *in,
                                          int in_stride,
                                          atfft_complex *out,
                                          int out_stride,
                                          int size)
{
    int last_bin = atfft_dft_halfcomplex_size (size);

    for (int i = 0, o = 0; i < last_bin * in_stride; i += in_stride, o += out_stride)
    {
        ATFFT_RE (out [o]) = ATFFT_RE (in [i]);
        ATFFT_IM (out [o]) = ATFFT_IM (in [i]);
    }

    if (atfft_is_even (size))
        --last_bin;

    for (int i = in_stride, o = (size - 1) * out_stride; i < last_bin * in_stride; i += in_stride, o -= out_stride)
    {
        ATFFT_RE (out [o]) = ATFFT_RE (in [i]);
        ATFFT_IM (out [o]) = - ATFFT_IM (in [i]);
    }
}


void atfft_complex_to_halfcomplex (atfft_complex *in, atfft_complex *out, int size)
{
    int last_bin = atfft_dft_halfcomplex_size (size);
    memcpy (out, in, last_bin * sizeof (*out));
}

void atfft_complex_to_halfcomplex_stride (atfft_complex *in,
                                          int in_stride,
                                          atfft_complex *out,
                                          int out_stride,
                                          int size)
{
    for (int i = 0, o = 0; i < atfft_dft_halfcomplex_size (size) * in_stride; i += in_stride, o += out_stride)
    {
        ATFFT_RE (out [o]) = ATFFT_RE (in [i]);
        ATFFT_IM (out [o]) = ATFFT_IM (in [i]);
    }
}

void atfft_float_to_sample_real (const float *in, atfft_sample *out, int size)
{
    for (int i = 0; i < size; ++i)
    {
        out [i] = in [i];
    }
}
