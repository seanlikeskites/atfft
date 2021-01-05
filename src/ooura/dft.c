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
#include <string.h>
#include <assert.h>
#include <math.h>
#include <atfft/dft.h>
#include "../ooura/ooura.h"

#ifdef ATFFT_TYPE_LONG_DOUBLE
#   ifdef _MSC_VER
#       pragma message(": warning: Ooura only supports double precision floating point, " \
                       "higher precision values will be demoted to double for FFT calculations.")
#   else
#       warning Ooura only supports double precision floating point, \
                higher precision values will be demoted to double for FFT calculations.
#   endif
#endif

struct atfft_dft
{
    int size;
    enum atfft_direction direction;
    int ooura_direction;
    enum atfft_format format;
    int data_size;
    double *data;
    int *work_area;
    double *tables;
};

struct atfft_dft* atfft_dft_create (int size, enum atfft_direction direction, enum atfft_format format)
{
    /* Ooura only supports sizes which are a power of 2. */
    assert (atfft_is_power_of_2 (size));

    struct atfft_dft *fft;

    if (!(fft = malloc (sizeof (*fft))))
        return NULL;

    fft->size = size;
    fft->direction = direction;
    fft->format = format;

    if (direction == ATFFT_FORWARD)
        fft->ooura_direction = -1;
    else
        fft->ooura_direction = 1;

    int work_size = 0;

    switch (format)
    {
        case ATFFT_COMPLEX:
            fft->data_size = 2 * size * sizeof (*(fft->data));
            work_size = (2 + (1 << (int) (log (size + 0.5) / log (2)) / 2)) * sizeof (*(fft->work_area));
            break;

        case ATFFT_REAL:
            fft->ooura_direction *= -1;
            fft->data_size = size * sizeof (*(fft->data));
            work_size = (2 + (1 << (int) (log (size / 2 + 0.5) / log (2)) / 2)) * sizeof (*(fft->work_area));
    }

    fft->data = malloc (fft->data_size);
    fft->work_area = malloc (work_size);
    fft->tables = malloc (size / 2 * sizeof (*(fft->tables)));

    /* clean up on failure */
    if (!(fft->data && fft->work_area && fft->tables))
    {
        atfft_dft_destroy (fft);
        fft = NULL;
    }
    else
    {
        fft->work_area [0] = 0;
    }

    return fft;
}

void atfft_dft_destroy (struct atfft_dft *fft)
{
    if (fft)
    {
        free (fft->tables);
        free (fft->work_area);
        free (fft->data);
        free (fft);
    }
}

void atfft_dft_complex_transform (struct atfft_dft *fft, atfft_complex *in, atfft_complex *out)
{
    /* Only to be used with complex FFTs. */
    assert (fft->format == ATFFT_COMPLEX);

#ifdef ATFFT_TYPE_DOUBLE
    memcpy (fft->data, in, fft->data_size);
#else
    atfft_sample_to_double_complex (in, (atfft_complex_d*) fft->data, fft->size);
#endif

    cdft (2 * fft->size, fft->ooura_direction, fft->data, fft->work_area, fft->tables);

#ifdef ATFFT_TYPE_DOUBLE
    memcpy (out, fft->data, fft->data_size);
#else
    atfft_double_to_sample_complex ((atfft_complex_d*) fft->data, out, fft->size);
#endif
}

void atfft_dft_complex_transform_stride (struct atfft_dft *fft,
                                         atfft_complex *in,
                                         int in_stride,
                                         atfft_complex *out,
                                         int out_stride)
{
    /* Only to be used with complex FFTs. */
    assert (fft->format == ATFFT_COMPLEX);

    atfft_sample_to_double_complex_stride (in,
                                           in_stride,
                                           (atfft_complex_d*) fft->data,
                                           1,
                                           fft->size);

    cdft (2 * fft->size, fft->ooura_direction, fft->data, fft->work_area, fft->tables);

    atfft_double_to_sample_complex_stride ((atfft_complex_d*) fft->data,
                                           1,
                                           out,
                                           out_stride,
                                           fft->size);
}

static void atfft_halfcomplex_ooura_to_fftw (const double *in,
                                             atfft_complex *out,
                                             int out_stride,
                                             int size)
{
    int half_size = size / 2;

    ATFFT_RE (out [0]) = in [0];
    ATFFT_IM (out [0]) = 0;

    int o = out_stride;

    for (int i = 1; i < half_size; ++i, o += out_stride)
    {
        ATFFT_RE (out [o]) = in [2 * i];
        ATFFT_IM (out [o]) = -in [2 * i + 1];
    }

    ATFFT_RE (out [o]) = in [1];
    ATFFT_IM (out [o]) = 0;
}

void atfft_dft_real_forward_transform (struct atfft_dft *fft, const atfft_sample *in, atfft_complex *out)
{
    /* Only to be used for forward real FFTs. */
    assert ((fft->format == ATFFT_REAL) && (fft->direction == ATFFT_FORWARD));

#ifdef ATFFT_TYPE_DOUBLE
    memcpy (fft->data, in, fft->data_size);
#else
    atfft_sample_to_double_real (in, fft->data, fft->size);
#endif

    rdft (fft->size, fft->ooura_direction, fft->data, fft->work_area, fft->tables);
    atfft_halfcomplex_ooura_to_fftw (fft->data, out, 1, fft->size);
}

void atfft_dft_real_forward_transform_stride (struct atfft_dft *fft,
                                              const atfft_sample *in,
                                              int in_stride,
                                              atfft_complex *out,
                                              int out_stride)
{
    /* Only to be used for forward real FFTs. */
    assert ((fft->format == ATFFT_REAL) && (fft->direction == ATFFT_FORWARD));

    atfft_sample_to_double_real_stride (in,
                                        in_stride,
                                        fft->data,
                                        1,
                                        fft->size);

    rdft (fft->size, fft->ooura_direction, fft->data, fft->work_area, fft->tables);
    atfft_halfcomplex_ooura_to_fftw (fft->data, out, out_stride, fft->size);
}

static void atfft_halfcomplex_fftw_to_ooura (atfft_complex *in,
                                             int in_stride,
                                             double *out,
                                             int size)
{
    int half_size = size / 2;

    out [0] = 2.0 * ATFFT_RE (in [0]);

    int i = in_stride;

    for (int o = 1; o < half_size; i += in_stride, ++o)
    {
        out [2 * o] = 2.0 * ATFFT_RE (in [i]);
        out [2 * o + 1] = 2.0 * - ATFFT_IM (in [i]);
    }

    out [1] = 2.0 * ATFFT_RE (in [i]);
}

void atfft_dft_real_backward_transform (struct atfft_dft *fft, atfft_complex *in, atfft_sample *out)
{
    /* Only to be used for backward real FFTs. */
    assert ((fft->format == ATFFT_REAL) && (fft->direction == ATFFT_BACKWARD));

    atfft_halfcomplex_fftw_to_ooura (in, 1, fft->data, fft->size);
    rdft (fft->size, fft->ooura_direction, fft->data, fft->work_area, fft->tables);

#ifdef ATFFT_TYPE_DOUBLE
    memcpy (out, fft->data, fft->data_size);
#else
    atfft_double_to_sample_real (fft->data, out, fft->size);
#endif
}

void atfft_dft_real_backward_transform_stride (struct atfft_dft *fft,
                                               atfft_complex *in,
                                               int in_stride,
                                               atfft_sample *out,
                                               int out_stride)
{
    /* Only to be used for backward real FFTs. */
    assert ((fft->format == ATFFT_REAL) && (fft->direction == ATFFT_BACKWARD));

    atfft_halfcomplex_fftw_to_ooura (in, in_stride, fft->data, fft->size);
    rdft (fft->size, fft->ooura_direction, fft->data, fft->work_area, fft->tables);
    atfft_double_to_sample_real_stride (fft->data,
                                        1,
                                        out,
                                        out_stride,
                                        fft->size);
}
