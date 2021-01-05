/*
 * Copyright (c) 2021 Sean Enderby <sean.enderby@gmail.com>
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

#include <string.h>
#include <atfft/dft_util.h>

int atfft_halfcomplex_size (int size)
{
    return size / 2 + 1;
}

void atfft_halfcomplex_to_complex (atfft_complex *in, atfft_complex *out, int size)
{
    int last_bin = atfft_halfcomplex_size (size);

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
    int last_bin = atfft_halfcomplex_size (size);

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
    int last_bin = atfft_halfcomplex_size (size);
    memcpy (out, in, last_bin * sizeof (*out));
}

void atfft_complex_to_halfcomplex_stride (atfft_complex *in,
                                          int in_stride,
                                          atfft_complex *out,
                                          int out_stride,
                                          int size)
{
    for (int i = 0, o = 0; i < atfft_halfcomplex_size (size) * in_stride; i += in_stride, o += out_stride)
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
