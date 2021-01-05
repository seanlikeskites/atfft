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

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_fft_complex.h>
#include <gsl/gsl_fft_real.h>
#include <gsl/gsl_fft_halfcomplex.h>
#include <atfft/dft.h>

#ifdef ATFFT_TYPE_LONG_DOUBLE
#   ifdef _MSC_VER
#       pragma message(": warning: GSL only supports double precision floating point, " \
                       "higher precision values will be demoted to double for FFT calculations.")
#   else
#       warning GSL only supports double precision floating point, \
                higher precision values will be demoted to double for FFT calculations.
#   endif
#endif

struct atfft_dft
{
    int size;
    enum atfft_direction direction;
    gsl_fft_direction gsl_direction;
    enum atfft_format format;
    gsl_vector *data;
    void *work_area;
    void *tables;
};

struct atfft_dft* atfft_dft_create (int size, enum atfft_direction direction, enum atfft_format format)
{
    struct atfft_dft *fft;

    if (!(fft = malloc (sizeof (*fft))))
        return NULL;

    fft->size = size;
    fft->direction = direction;
    fft->format = format;

    switch (format)
    {
        case ATFFT_COMPLEX:
            fft->data = gsl_vector_alloc (2 * size);
            fft->work_area = gsl_fft_complex_workspace_alloc (size);
            fft->tables = gsl_fft_complex_wavetable_alloc (size);

            if (direction == ATFFT_FORWARD)
                fft->gsl_direction = -1;
            else
                fft->gsl_direction = 1;

            break;

        case ATFFT_REAL:
            fft->data = gsl_vector_alloc (size);
            fft->work_area = gsl_fft_real_workspace_alloc (size);
            fft->gsl_direction = 0; /* unused */

            if (direction == ATFFT_FORWARD)
                fft->tables = gsl_fft_real_wavetable_alloc (size);
            else
                fft->tables = gsl_fft_halfcomplex_wavetable_alloc (size);

            break;
    }

    /* clean up on failure (might be bad for old versions of GSL) */
    if (!(fft->data && fft->work_area && fft->tables))
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
        switch (fft->format)
        {
            case ATFFT_COMPLEX:
                gsl_fft_complex_wavetable_free (fft->tables);
                gsl_fft_complex_workspace_free (fft->work_area);
                break;

            case ATFFT_REAL:
                if (fft->direction == ATFFT_FORWARD)
                    gsl_fft_real_wavetable_free (fft->tables);
                else
                    gsl_fft_halfcomplex_wavetable_free (fft->tables);

                gsl_fft_real_workspace_free (fft->work_area);
                break;
        }

        gsl_vector_free (fft->data);

        free (fft);
    }
}

void atfft_dft_complex_transform (struct atfft_dft *fft, atfft_complex *in, atfft_complex *out)
{
    /* Only to be used with complex FFTs. */
    assert (fft->format == ATFFT_COMPLEX);

#ifdef ATFFT_TYPE_DOUBLE
    size_t n_bytes = fft->size * sizeof (*in);
#endif 

    gsl_complex_packed_array data = fft->data->data;

#ifdef ATFFT_TYPE_DOUBLE
    memcpy (data, in, n_bytes);
#else
    atfft_sample_to_double_complex (in, (atfft_complex_d*) data, fft->size);
#endif

    gsl_fft_complex_transform (data, fft->data->stride, fft->size, fft->tables, fft->work_area, fft->gsl_direction);

#ifdef ATFFT_TYPE_DOUBLE
    memcpy (out, data, n_bytes);
#else
    atfft_double_to_sample_complex ((atfft_complex_d*) data, out, fft->size);
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

    gsl_complex_packed_array data = fft->data->data;

    atfft_sample_to_double_complex_stride (in,
                                           in_stride,
                                           (atfft_complex_d*) data,
                                           1,
                                           fft->size);

    gsl_fft_complex_transform (data, fft->data->stride, fft->size, fft->tables, fft->work_area, fft->gsl_direction);

    atfft_double_to_sample_complex_stride ((atfft_complex_d*) data,
                                           1,
                                           out,
                                           out_stride,
                                           fft->size);
}

static void atfft_halfcomplex_gsl_to_fftw (const double *in,
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
        ATFFT_RE (out [o]) = in [2 * i - 1];
        ATFFT_IM (out [o]) = in [2 * i];
    }

    ATFFT_RE (out [o]) = in [size - 1];
    ATFFT_IM (out [o]) = 0;
}

void atfft_dft_real_forward_transform (struct atfft_dft *fft, const atfft_sample *in, atfft_complex *out)
{
    /* Only to be used for forward real FFTs. */
    assert ((fft->format == ATFFT_REAL) && (fft->direction == ATFFT_FORWARD));

    double *data = fft->data->data;

#ifdef ATFFT_TYPE_DOUBLE
    memcpy (data, in, fft->size * sizeof (*data));
#else
    atfft_sample_to_double_real (in, data, fft->size);
#endif

    gsl_fft_real_transform (data, fft->data->stride, fft->size, fft->tables, fft->work_area);
    atfft_halfcomplex_gsl_to_fftw (data, out, 1, fft->size);
}

void atfft_dft_real_forward_transform_stride (struct atfft_dft *fft,
                                              const atfft_sample *in,
                                              int in_stride,
                                              atfft_complex *out,
                                              int out_stride)
{
    /* Only to be used for forward real FFTs. */
    assert ((fft->format == ATFFT_REAL) && (fft->direction == ATFFT_FORWARD));

    double *data = fft->data->data;

    atfft_sample_to_double_real_stride (in,
                                        in_stride,
                                        data,
                                        1,
                                        fft->size);

    gsl_fft_real_transform (data, fft->data->stride, fft->size, fft->tables, fft->work_area);
    atfft_halfcomplex_gsl_to_fftw (data, out, out_stride, fft->size);
}

static void atfft_halfcomplex_fftw_to_gsl (atfft_complex *in,
                                           int in_stride,
                                           double *out,
                                           int size)
{
    int half_size = size / 2;

    out [0] = ATFFT_RE (in [0]);

    int i = in_stride;

    for (int o = 1; o < half_size; i += in_stride, ++o)
    {
        out [2 * o - 1] = ATFFT_RE (in [i]);
        out [2 * o] = ATFFT_IM (in [i]);
    }

    out [size - 1] = ATFFT_RE (in [i]);
}

void atfft_dft_real_backward_transform (struct atfft_dft *fft, atfft_complex *in, atfft_sample *out)
{
    /* Only to be used for backward real FFTs. */
    assert ((fft->format == ATFFT_REAL) && (fft->direction == ATFFT_BACKWARD));

    double *data = fft->data->data;

    atfft_halfcomplex_fftw_to_gsl (in, 1, data, fft->size);
    gsl_fft_halfcomplex_transform (data, fft->data->stride, fft->size, fft->tables, fft->work_area);   

#ifdef ATFFT_TYPE_DOUBLE
    memcpy (out, data, fft->size * sizeof (*data));
#else
    atfft_double_to_sample_real (data, out, fft->size);
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

    double *data = fft->data->data;

    atfft_halfcomplex_fftw_to_gsl (in, in_stride, data, fft->size);
    gsl_fft_halfcomplex_transform (data, fft->data->stride, fft->size, fft->tables, fft->work_area);   
    atfft_double_to_sample_real_stride (data,
                                        1,
                                        out,
                                        out_stride,
                                        fft->size);
}
