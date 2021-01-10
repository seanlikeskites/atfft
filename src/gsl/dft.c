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

    /* buffer for in place transform */
    size_t n_data_bytes;
    double *data;

    /* buffers for gsl internals */
    void *work_area;
    void *tables;
};

int atfft_dft_is_supported_size (int size, enum atfft_format format)
{
    return size > 0;
}

struct atfft_dft* atfft_dft_create (int size, enum atfft_direction direction, enum atfft_format format)
{
    /* gsl supports all sizes. */
    assert (atfft_dft_is_supported_size (size, format));

    struct atfft_dft *plan;

    if (!(plan = calloc (1, sizeof (*plan))))
        return NULL;

    plan->size = size;
    plan->direction = direction;
    plan->format = format;

    switch (format)
    {
        case ATFFT_COMPLEX:
            plan->n_data_bytes = 2 * size * sizeof (*(plan->data));
            plan->work_area = gsl_fft_complex_workspace_alloc (size);
            plan->tables = gsl_fft_complex_wavetable_alloc (size);

            if (direction == ATFFT_FORWARD)
                plan->gsl_direction = -1;
            else
                plan->gsl_direction = 1;

            break;

        case ATFFT_REAL:
            plan->n_data_bytes = size * sizeof (*(plan->data));
            plan->work_area = gsl_fft_real_workspace_alloc (size);

            if (direction == ATFFT_FORWARD)
                plan->tables = gsl_fft_real_wavetable_alloc (size);
            else
                plan->tables = gsl_fft_halfcomplex_wavetable_alloc (size);

            break;
    }

    plan->data = malloc (plan->n_data_bytes);

    /* clean up on failure (might be bad for old versions of GSL) */
    if (!(plan->data && plan->work_area && plan->tables))
    {
        atfft_dft_destroy (plan);
        plan = NULL;
    }

    return plan;
}

void atfft_dft_destroy (struct atfft_dft *plan)
{
    if (plan)
    {
        switch (plan->format)
        {
            case ATFFT_COMPLEX:
                gsl_fft_complex_wavetable_free (plan->tables);
                gsl_fft_complex_workspace_free (plan->work_area);
                break;

            case ATFFT_REAL:
                if (plan->direction == ATFFT_FORWARD)
                    gsl_fft_real_wavetable_free (plan->tables);
                else
                    gsl_fft_halfcomplex_wavetable_free (plan->tables);

                gsl_fft_real_workspace_free (plan->work_area);
                break;
        }

        free (plan->data);

        free (plan);
    }
}

void atfft_dft_complex_transform (struct atfft_dft *plan, atfft_complex *in, atfft_complex *out)
{
    /* Only to be used with complex FFTs. */
    assert (plan->format == ATFFT_COMPLEX);

    gsl_complex_packed_array data = plan->data;

#ifdef ATFFT_TYPE_DOUBLE
    memcpy (data, in, plan->n_data_bytes);
#else
    atfft_sample_to_double_complex (in, (atfft_complex_d*) data, plan->size);
#endif

    gsl_fft_complex_transform (data, 1, plan->size, plan->tables, plan->work_area, plan->gsl_direction);

#ifdef ATFFT_TYPE_DOUBLE
    memcpy (out, data, plan->n_data_bytes);
#else
    atfft_double_to_sample_complex ((atfft_complex_d*) data, out, plan->size);
#endif
}

void atfft_dft_complex_transform_stride (struct atfft_dft *plan,
                                         atfft_complex *in,
                                         int in_stride,
                                         atfft_complex *out,
                                         int out_stride)
{
    /* Only to be used with complex FFTs. */
    assert (plan->format == ATFFT_COMPLEX);

    gsl_complex_packed_array data = plan->data;

    atfft_sample_to_double_complex_stride (in,
                                           in_stride,
                                           (atfft_complex_d*) data,
                                           1,
                                           plan->size);

    gsl_fft_complex_transform (data, 1, plan->size, plan->tables, plan->work_area, plan->gsl_direction);

    atfft_double_to_sample_complex_stride ((atfft_complex_d*) data,
                                           1,
                                           out,
                                           out_stride,
                                           plan->size);
}

static void halfcomplex_gsl_to_atfft (const double *in,
                                      atfft_complex *out,
                                      int out_stride,
                                      int size)
{
    ATFFT_RE (out [0]) = in [0];
    ATFFT_IM (out [0]) = 0;

    int half_size = size / 2;
    int o = out_stride;
    int i = 1;

    for (; i < half_size; ++i, o += out_stride)
    {
        ATFFT_RE (out [o]) = in [2 * i - 1];
        ATFFT_IM (out [o]) = in [2 * i];
    }

    if (atfft_is_even (size))
    {
        ATFFT_RE (out [o]) = in [size - 1];
        ATFFT_IM (out [o]) = 0;
    }
    else if (size > 1)
    {
        ATFFT_RE (out [o]) = in [2 * i - 1];
        ATFFT_IM (out [o]) = in [2 * i];
    }
}

void atfft_dft_real_forward_transform (struct atfft_dft *plan, const atfft_sample *in, atfft_complex *out)
{
    /* Only to be used for forward real FFTs. */
    assert ((plan->format == ATFFT_REAL) && (plan->direction == ATFFT_FORWARD));

#ifdef ATFFT_TYPE_DOUBLE
    memcpy (plan->data, in, plan->n_data_bytes);
#else
    atfft_sample_to_double_real (in, plan->data, plan->size);
#endif

    gsl_fft_real_transform (plan->data, 1, plan->size, plan->tables, plan->work_area);
    halfcomplex_gsl_to_atfft (plan->data, out, 1, plan->size);
}

void atfft_dft_real_forward_transform_stride (struct atfft_dft *plan,
                                              const atfft_sample *in,
                                              int in_stride,
                                              atfft_complex *out,
                                              int out_stride)
{
    /* Only to be used for forward real FFTs. */
    assert ((plan->format == ATFFT_REAL) && (plan->direction == ATFFT_FORWARD));

    atfft_sample_to_double_real_stride (in,
                                        in_stride,
                                        plan->data,
                                        1,
                                        plan->size);

    gsl_fft_real_transform (plan->data, 1, plan->size, plan->tables, plan->work_area);
    halfcomplex_gsl_to_atfft (plan->data, out, out_stride, plan->size);
}

static void halfcomplex_atfft_to_gsl (atfft_complex *in,
                                      int in_stride,
                                      double *out,
                                      int size)
{
    out [0] = ATFFT_RE (in [0]);

    int half_size = size / 2;
    int i = in_stride;
    int o = 1;

    for (; o < half_size; i += in_stride, ++o)
    {
        out [2 * o - 1] = ATFFT_RE (in [i]);
        out [2 * o] = ATFFT_IM (in [i]);
    }

    if (atfft_is_even (size))
    {
        out [size - 1] = ATFFT_RE (in [i]);
    }
    else if (size > 1)
    {
        out [2 * o - 1] = ATFFT_RE (in [i]);
        out [2 * o] = ATFFT_IM (in [i]);
    }
}

void atfft_dft_real_backward_transform (struct atfft_dft *plan, atfft_complex *in, atfft_sample *out)
{
    /* Only to be used for backward real FFTs. */
    assert ((plan->format == ATFFT_REAL) && (plan->direction == ATFFT_BACKWARD));

    halfcomplex_atfft_to_gsl (in, 1, plan->data, plan->size);
    gsl_fft_halfcomplex_transform (plan->data, 1, plan->size, plan->tables, plan->work_area);   

#ifdef ATFFT_TYPE_DOUBLE
    memcpy (out, plan->data, plan->n_data_bytes);
#else
    atfft_double_to_sample_real (plan->data, out, plan->size);
#endif
}

void atfft_dft_real_backward_transform_stride (struct atfft_dft *plan,
                                               atfft_complex *in,
                                               int in_stride,
                                               atfft_sample *out,
                                               int out_stride)
{
    /* Only to be used for backward real FFTs. */
    assert ((plan->format == ATFFT_REAL) && (plan->direction == ATFFT_BACKWARD));

    halfcomplex_atfft_to_gsl (in, in_stride, plan->data, plan->size);
    gsl_fft_halfcomplex_transform (plan->data, 1, plan->size, plan->tables, plan->work_area);   
    atfft_double_to_sample_real_stride (plan->data,
                                        1,
                                        out,
                                        out_stride,
                                        plan->size);
}
