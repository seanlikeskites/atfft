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

#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <limits.h>
#include <atfft/dft.h>
#include <atfft/dft_util.h>
#include "atfft_internal.h"
#include "dft_cooley_tukey.h"
#include "dft_rader.h"
#include "dft_bluestein.h"
#include "dft_pfa.h"

typedef void (*complex_transform_function) (void*, atfft_complex*, int, atfft_complex*, int);
typedef void (*fft_destroy_function) (void*);

struct atfft_dft
{
    int size;
    enum atfft_direction direction;
    enum atfft_format format;

    void *fft;
    complex_transform_function complex_transform;
    fft_destroy_function fft_destroy;

    atfft_complex *real_in, *real_out;
};

int atfft_dft_is_supported_size (int size, enum atfft_format format)
{
    return size > 0;
}

struct atfft_dft* atfft_dft_create (int size, enum atfft_direction direction, enum atfft_format format)
{
    /* check size isn't negative */
    assert (atfft_dft_is_supported_size (size, format));

    struct atfft_dft *fft;

    if (!(fft = calloc (1, sizeof (*fft))))
        return NULL;

    fft->size = size;
    fft->direction = direction;
    fft->format = format;

    if (atfft_is_prime (size))
    {
        if (atfft_is_power_of_2 (size - 1))
        {
            /* Use Rader's algorithm */
            fft->fft = atfft_dft_rader_create (size, direction, ATFFT_COMPLEX);
            fft->complex_transform = atfft_dft_rader_complex_transform;
            fft->fft_destroy = atfft_dft_rader_destroy;
        }
        else
        {
            /* Use Bluestein's algorithm */
            fft->fft = atfft_dft_bluestein_create (size, direction, ATFFT_COMPLEX);
            fft->complex_transform = atfft_dft_bluestein_complex_transform;
            fft->fft_destroy = atfft_dft_bluestein_destroy;
        }
    }
    else
    {
        /* Use Cooley-Tukey */
        fft->fft = atfft_dft_cooley_tukey_create (size, direction, ATFFT_COMPLEX);
        fft->complex_transform = atfft_dft_cooley_tukey_complex_transform;
        fft->fft_destroy = atfft_dft_cooley_tukey_destroy;
    }

    if (!fft->fft)
        goto failed;

    if (format == ATFFT_REAL)
    {
        fft->real_in = malloc (size * sizeof (*(fft->real_in)));
        fft->real_out = malloc (size * sizeof (*(fft->real_out)));

        if (!(fft->real_in && fft->real_out))
            goto failed;
    }

    return fft;

failed:
    atfft_dft_destroy (fft);
    return NULL;
}

void atfft_dft_destroy (struct atfft_dft *fft)
{
    if (fft)
    {
        free (fft->real_out);
        free (fft->real_in);
        fft->fft_destroy (fft->fft);
        free (fft);
    }
}

void atfft_dft_complex_transform (struct atfft_dft *fft,
                                  atfft_complex *in,
                                  atfft_complex *out)
{
    atfft_dft_complex_transform_stride (fft, in, 1, out, 1);
}

void atfft_dft_complex_transform_stride (struct atfft_dft *fft,
                                         atfft_complex *in,
                                         int in_stride,
                                         atfft_complex *out,
                                         int out_stride)
{
    /* Only to be used with complex FFTs. */
    assert (fft->format == ATFFT_COMPLEX);

    fft->complex_transform (fft->fft, in, in_stride, out, out_stride);
}

void atfft_dft_real_forward_transform (struct atfft_dft *fft, const atfft_sample *in, atfft_complex *out)
{
    /* Only to be used for forward real FFTs. */
    assert ((fft->format == ATFFT_REAL) && (fft->direction == ATFFT_FORWARD));

    atfft_real_to_complex (in, fft->real_in, fft->size);
    fft->complex_transform (fft->fft, fft->real_in, 1, fft->real_out, 1);
    atfft_complex_to_halfcomplex (fft->real_out, out, fft->size);
}

void atfft_dft_real_forward_transform_stride (struct atfft_dft *fft,
                                              const atfft_sample *in,
                                              int in_stride,
                                              atfft_complex *out,
                                              int out_stride)
{
    /* Only to be used for forward real FFTs. */
    assert ((fft->format == ATFFT_REAL) && (fft->direction == ATFFT_FORWARD));

    atfft_real_to_complex_stride (in, in_stride, fft->real_in, 1, fft->size);
    fft->complex_transform (fft->fft, fft->real_in, 1, fft->real_out, 1);
    atfft_complex_to_halfcomplex_stride (fft->real_out, 1, out, out_stride, fft->size);
}

void atfft_dft_real_backward_transform (struct atfft_dft *fft, atfft_complex *in, atfft_sample *out)
{
    /* Only to be used for backward real FFTs. */
    assert ((fft->format == ATFFT_REAL) && (fft->direction == ATFFT_BACKWARD));

    atfft_halfcomplex_to_complex (in, fft->real_in, fft->size);
    fft->complex_transform (fft->fft, fft->real_in, 1, fft->real_out, 1);
    atfft_real (fft->real_out, out, fft->size);
}

void atfft_dft_real_backward_transform_stride (struct atfft_dft *fft,
                                               atfft_complex *in,
                                               int in_stride,
                                               atfft_sample *out,
                                               int out_stride)
{
    /* Only to be used for backward real FFTs. */
    assert ((fft->format == ATFFT_REAL) && (fft->direction == ATFFT_BACKWARD));

    atfft_halfcomplex_to_complex_stride (in, in_stride, fft->real_in, 1, fft->size);
    fft->complex_transform (fft->fft, fft->real_in, 1, fft->real_out, 1);
    atfft_real_stride (fft->real_out, 1, out, out_stride, fft->size);
}
