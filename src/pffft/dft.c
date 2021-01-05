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
#include <atfft/dft.h>
#include "../pffft/pffft.h"

#ifndef ATFFT_TYPE_FLOAT
#   ifdef _MSC_VER
#       pragma message(": warning: PFFFT only supports single precision floating point, " \
                       "higher precision values will be demoted to float for FFT calculations.")
#   else
#       warning PFFFT only supports single precision floating point, \
                higher precision values will be demoted to float for FFT calculations.
#   endif
#endif

struct atfft_dft
{
    int size;
    enum atfft_direction direction;
    pffft_direction_t pffft_direction;
    enum atfft_format format;
    int in_size, out_size;
    float *in, *out;
    float *work_area;
    PFFFT_Setup *plan;
};

static int atfft_is_supported_size_pffft (unsigned int size, enum atfft_format format)
{
    int min = 16;

    if (format == ATFFT_REAL)
        min = 32;

    if (!(size % min) && (size > 0))
        return 1;        
    else
        return 0;
}

struct atfft_dft* atfft_dft_create (int size, enum atfft_direction direction, enum atfft_format format)
{
    /* pffft only supports sizes which are a multiples of 32 (or 16 for complex transforms). */
    assert (atfft_is_supported_size_pffft (size, format));

    struct atfft_dft *fft;

    if (!(fft = malloc (sizeof (*fft))))
        return NULL;

    fft->size = size;
    fft->direction = direction;
    fft->format = format;

    if (direction == ATFFT_FORWARD)
        fft->pffft_direction = PFFFT_FORWARD;
    else
        fft->pffft_direction = PFFFT_BACKWARD;

    int work_size = 0;
    pffft_transform_t transform_type = 0;

    switch (format)
    {
        case ATFFT_COMPLEX:
            fft->in_size = 2 * size * sizeof (*(fft->in));
            fft->out_size = 2 * size * sizeof (*(fft->out));
            work_size = 2 * size * sizeof (*(fft->work_area));
            transform_type = PFFFT_COMPLEX;
            break;

        case ATFFT_REAL:
            fft->in_size = size * sizeof (*(fft->in));
            fft->out_size = size * sizeof (*(fft->out));
            work_size = size * sizeof (*(fft->work_area));
            transform_type = PFFFT_REAL;
    }

    fft->in = pffft_aligned_malloc (fft->in_size);
    fft->out = pffft_aligned_malloc (fft->out_size);
    fft->work_area = pffft_aligned_malloc (work_size);
    fft->plan = pffft_new_setup (size, transform_type);

    /* clean up on failure */
    if (!(fft->in && fft->out && fft->work_area && fft->plan))
    {
        atfft_dft_destroy (fft);
        fft = NULL;
    }

    return fft;
}

void atfft_dft_destroy (struct atfft_dft *fft)
{
    if (fft)
    {
        pffft_destroy_setup (fft->plan);
        pffft_aligned_free (fft->work_area);
        pffft_aligned_free (fft->out);
        pffft_aligned_free (fft->in);
        free (fft);
    }
}

void atfft_dft_complex_transform (struct atfft_dft *fft, atfft_complex *in, atfft_complex *out)
{
    /* Only to be used with complex FFTs. */
    assert (fft->format == ATFFT_COMPLEX);

#ifdef ATFFT_TYPE_FLOAT
    memcpy (fft->in, in, fft->in_size);
#else
    atfft_sample_to_float_complex (in, (atfft_complex_f*) fft->in, fft->size);
#endif

    pffft_transform_ordered (fft->plan, fft->in, fft->out, fft->work_area, fft->pffft_direction);

#ifdef ATFFT_TYPE_FLOAT
    memcpy (out, fft->out, fft->out_size);
#else
    atfft_float_to_sample_complex ((atfft_complex_f*) fft->out, out, fft->size);
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

    atfft_sample_to_float_complex_stride (in,
                                          in_stride,
                                          (atfft_complex_f*) fft->in,
                                          1,
                                          fft->size);

    pffft_transform_ordered (fft->plan, fft->in, fft->out, fft->work_area, fft->pffft_direction);

    atfft_float_to_sample_complex_stride ((atfft_complex_f*) fft->out,
                                          1,
                                          out,
                                          out_stride,
                                          fft->size);
}

static void atfft_halfcomplex_pffft_to_fftw (const float *in,
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
        ATFFT_IM (out [o]) = in [2 * i + 1];
    }

    ATFFT_RE (out [o]) = in [1];
    ATFFT_IM (out [o]) = 0;
}

void atfft_dft_real_forward_transform (struct atfft_dft *fft, const atfft_sample *in, atfft_complex *out)
{
    /* Only to be used for forward real FFTs. */
    assert ((fft->format == ATFFT_REAL) && (fft->direction == ATFFT_FORWARD));

#ifdef ATFFT_TYPE_FLOAT
    memcpy (fft->in, in, fft->in_size);
#else
    atfft_sample_to_float_real (in, fft->in, fft->size);
#endif

    pffft_transform_ordered (fft->plan, fft->in, fft->out, fft->work_area, fft->pffft_direction);
    atfft_halfcomplex_pffft_to_fftw (fft->out, out, 1, fft->size);
}

void atfft_dft_real_forward_transform_stride (struct atfft_dft *fft,
                                              const atfft_sample *in,
                                              int in_stride,
                                              atfft_complex *out,
                                              int out_stride)
{
    /* Only to be used for forward real FFTs. */
    assert ((fft->format == ATFFT_REAL) && (fft->direction == ATFFT_FORWARD));

    atfft_sample_to_float_real_stride (in,
                                       in_stride,
                                       fft->in,
                                       1,
                                       fft->size);

    pffft_transform_ordered (fft->plan, fft->in, fft->out, fft->work_area, fft->pffft_direction);
    atfft_halfcomplex_pffft_to_fftw (fft->out, out, out_stride, fft->size);
}

static void atfft_halfcomplex_fftw_to_pffft (atfft_complex *in,
                                             int in_stride,
                                             float *out,
                                             int size)
{
    int half_size = size / 2;

    out [0] = ATFFT_RE (in [0]);

    int i = in_stride;

    for (int o = 1; o < half_size; i += in_stride, ++o)
    {
        out [2 * o] = ATFFT_RE (in [i]);
        out [2 * o + 1] = ATFFT_IM (in [i]);
    }

    out [1] = ATFFT_RE (in [i]);
}

void atfft_dft_real_backward_transform (struct atfft_dft *fft, atfft_complex *in, atfft_sample *out)
{
    /* Only to be used for backward real FFTs. */
    assert ((fft->format == ATFFT_REAL) && (fft->direction == ATFFT_BACKWARD));

    atfft_halfcomplex_fftw_to_pffft (in, 1, fft->in, fft->size);
    pffft_transform_ordered (fft->plan, fft->in, fft->out, fft->work_area, fft->pffft_direction);

#ifdef ATFFT_TYPE_FLOAT
    memcpy (out, fft->out, fft->out_size);
#else
    atfft_float_to_sample_real (fft->out, out, fft->size);
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

    atfft_halfcomplex_fftw_to_pffft (in, in_stride, fft->in, fft->size);
    pffft_transform_ordered (fft->plan, fft->in, fft->out, fft->work_area, fft->pffft_direction);
    atfft_float_to_sample_real_stride (fft->out,
                                       1,
                                       out,
                                       out_stride,
                                       fft->size);
}
