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
#include "ooura.h"

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

    /* buffer for in place transform */
    size_t n_data_bytes;
    double *data;

    /* buffers for ooura internals */
    int *work_area;
    double *tables;
};

int atfft_dft_is_supported_size (int size, enum atfft_format format)
{
    int min = 1;

    if (format == ATFFT_REAL)
        min = 2;

    return atfft_is_power_of_2 (size) && size >= min;
}

struct atfft_dft* atfft_dft_create (int size, enum atfft_direction direction, enum atfft_format format)
{
    /* ooura only supports sizes which are a power of 2. */
    assert (atfft_dft_is_supported_size (size, format));

    struct atfft_dft *plan;

    if (!(plan = calloc (1, sizeof (*plan))))
        return NULL;

    plan->size = size;
    plan->direction = direction;
    plan->format = format;

    if (direction == ATFFT_FORWARD)
        plan->ooura_direction = -1;
    else
        plan->ooura_direction = 1;

    int work_size = 0;

    if (format == ATFFT_COMPLEX)
    {
        plan->n_data_bytes = 2 * size * sizeof (*(plan->data));
        work_size = (2 + (1 << (int) (log (size + 0.5) / log (2)) / 2));
    }
    else
    {
        plan->ooura_direction *= -1;
        plan->n_data_bytes = size * sizeof (*(plan->data));
        work_size = (2 + (1 << (int) (log (size / 2 + 0.5) / log (2)) / 2));
    }

    plan->data = malloc (plan->n_data_bytes);
    plan->work_area = malloc (work_size * sizeof (*(plan->work_area)));
    plan->tables = malloc (size / 2 * sizeof (*(plan->tables)));

    if (!(plan->data && plan->work_area && plan->tables))
        goto failed;

    /* run a transform to initialise ooura state */
    plan->work_area [0] = 0;

    if (format == ATFFT_COMPLEX)
        cdft (2 * plan->size, plan->ooura_direction, plan->data, plan->work_area, plan->tables);
    else
        rdft (plan->size, plan->ooura_direction, plan->data, plan->work_area, plan->tables);

    return plan;

failed:
    atfft_dft_destroy (plan);
    return NULL;
}

void atfft_dft_destroy (struct atfft_dft *plan)
{
    if (plan)
    {
        free (plan->tables);
        free (plan->work_area);
        free (plan->data);
        free (plan);
    }
}

void atfft_dft_complex_transform (struct atfft_dft *plan, atfft_complex *in, atfft_complex *out)
{
    /* Only to be used with complex FFTs. */
    assert (plan->format == ATFFT_COMPLEX);

#ifdef ATFFT_TYPE_DOUBLE
    memcpy (plan->data, in, plan->n_data_bytes);
#else
    atfft_sample_to_double_complex (in, (atfft_complex_d*) plan->data, plan->size);
#endif

    cdft (2 * plan->size, plan->ooura_direction, plan->data, plan->work_area, plan->tables);

#ifdef ATFFT_TYPE_DOUBLE
    memcpy (out, plan->data, plan->n_data_bytes);
#else
    atfft_double_to_sample_complex ((atfft_complex_d*) plan->data, out, plan->size);
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

    atfft_sample_to_double_complex_stride (in,
                                           in_stride,
                                           (atfft_complex_d*) plan->data,
                                           1,
                                           plan->size);

    cdft (2 * plan->size, plan->ooura_direction, plan->data, plan->work_area, plan->tables);

    atfft_double_to_sample_complex_stride ((atfft_complex_d*) plan->data,
                                           1,
                                           out,
                                           out_stride,
                                           plan->size);
}

static void halfcomplex_ooura_to_atfft (const double *in,
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

void atfft_dft_real_forward_transform (struct atfft_dft *plan, const atfft_sample *in, atfft_complex *out)
{
    /* Only to be used for forward real FFTs. */
    assert ((plan->format == ATFFT_REAL) && (plan->direction == ATFFT_FORWARD));

#ifdef ATFFT_TYPE_DOUBLE
    memcpy (plan->data, in, plan->n_data_bytes);
#else
    atfft_sample_to_double_real (in, plan->data, plan->size);
#endif

    rdft (plan->size, plan->ooura_direction, plan->data, plan->work_area, plan->tables);
    halfcomplex_ooura_to_atfft (plan->data, out, 1, plan->size);
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

    rdft (plan->size, plan->ooura_direction, plan->data, plan->work_area, plan->tables);
    halfcomplex_ooura_to_atfft (plan->data, out, out_stride, plan->size);
}

static void halfcomplex_atfft_to_ooura (atfft_complex *in,
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

void atfft_dft_real_backward_transform (struct atfft_dft *plan, atfft_complex *in, atfft_sample *out)
{
    /* Only to be used for backward real FFTs. */
    assert ((plan->format == ATFFT_REAL) && (plan->direction == ATFFT_BACKWARD));

    halfcomplex_atfft_to_ooura (in, 1, plan->data, plan->size);
    rdft (plan->size, plan->ooura_direction, plan->data, plan->work_area, plan->tables);

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

    halfcomplex_atfft_to_ooura (in, in_stride, plan->data, plan->size);
    rdft (plan->size, plan->ooura_direction, plan->data, plan->work_area, plan->tables);
    atfft_double_to_sample_real_stride (plan->data,
                                        1,
                                        out,
                                        out_stride,
                                        plan->size);
}
