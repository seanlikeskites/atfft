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
#include "pffft.h"

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

    /* the pffft plan */
    PFFFT_Setup *plan;

    /* aligned input and output buffers for pffft transform */
    size_t n_in_out_bytes;
    float *in, *out;

    /* additional work buffer for pffft */
    float *work_area;
};

int atfft_dft_is_supported_size (int size, enum atfft_format format)
{
    if (size <= 0)
        return 0;

    /* remove any factors of 5 and 3 */
    int div = size;

    while (div % 5 == 0)
    {
        div /= 5;
    }

    while (div % 3 == 0)
    {
        div /= 3;
    }

    /* we should be left with a power of 2 */
    int min = 16;

    if (format == ATFFT_REAL)
        min = 32;

    return atfft_is_power_of_2 (div) && div >= min;
}

struct atfft_dft* atfft_dft_create (int size, enum atfft_direction direction, enum atfft_format format)
{
    /* pffft only supports sizes of the form:
     * 16 * 2^a * 3^b * 5^c for complex transforms
     * 32 * 2^a * 3^b * 5^c for real transforms. */
    assert (atfft_dft_is_supported_size (size, format));

    struct atfft_dft *plan;

    if (!(plan = calloc (1, sizeof (*plan))))
        return NULL;

    plan->size = size;
    plan->direction = direction;
    plan->format = format;

    if (direction == ATFFT_FORWARD)
        plan->pffft_direction = PFFFT_FORWARD;
    else
        plan->pffft_direction = PFFFT_BACKWARD;

    int work_size = 0;
    pffft_transform_t transform_type = 0;

    if (format == ATFFT_COMPLEX)
    {
        plan->n_in_out_bytes = 2 * size * sizeof (*(plan->in));
        work_size = 2 * size;
        transform_type = PFFFT_COMPLEX;
    }
    else
    {
        plan->n_in_out_bytes = size * sizeof (*(plan->in));
        work_size = size;
        transform_type = PFFFT_REAL;
    }

    plan->plan = pffft_new_setup (size, transform_type);
    plan->in = pffft_aligned_malloc (plan->n_in_out_bytes);
    plan->out = pffft_aligned_malloc (plan->n_in_out_bytes);
    plan->work_area = pffft_aligned_malloc (work_size * sizeof (*plan->work_area));

    /* clean up on failure */
    if (!(plan->plan && plan->in && plan->out && plan->work_area))
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
        pffft_destroy_setup (plan->plan);
        pffft_aligned_free (plan->work_area);
        pffft_aligned_free (plan->out);
        pffft_aligned_free (plan->in);
        free (plan);
    }
}

void atfft_dft_complex_transform (struct atfft_dft *plan, atfft_complex *in, atfft_complex *out)
{
    /* Only to be used with complex FFTs. */
    assert (plan->format == ATFFT_COMPLEX);

#ifdef ATFFT_TYPE_FLOAT
    memcpy (plan->in, in, plan->n_in_out_bytes);
#else
    atfft_sample_to_float_complex (in, (atfft_complex_f*) plan->in, plan->size);
#endif

    pffft_transform_ordered (plan->plan, plan->in, plan->out, plan->work_area, plan->pffft_direction);

#ifdef ATFFT_TYPE_FLOAT
    memcpy (out, plan->out, plan->n_in_out_bytes);
#else
    atfft_float_to_sample_complex ((atfft_complex_f*) plan->out, out, plan->size);
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

    atfft_sample_to_float_complex_stride (in,
                                          in_stride,
                                          (atfft_complex_f*) plan->in,
                                          1,
                                          plan->size);

    pffft_transform_ordered (plan->plan, plan->in, plan->out, plan->work_area, plan->pffft_direction);

    atfft_float_to_sample_complex_stride ((atfft_complex_f*) plan->out,
                                          1,
                                          out,
                                          out_stride,
                                          plan->size);
}

static void halfcomplex_pffft_to_atfft (const float *in,
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

void atfft_dft_real_forward_transform (struct atfft_dft *plan, const atfft_sample *in, atfft_complex *out)
{
    /* Only to be used for forward real FFTs. */
    assert ((plan->format == ATFFT_REAL) && (plan->direction == ATFFT_FORWARD));

#ifdef ATFFT_TYPE_FLOAT
    memcpy (plan->in, in, plan->n_in_out_bytes);
#else
    atfft_sample_to_float_real (in, plan->in, plan->size);
#endif

    pffft_transform_ordered (plan->plan, plan->in, plan->out, plan->work_area, plan->pffft_direction);
    halfcomplex_pffft_to_atfft (plan->out, out, 1, plan->size);
}

void atfft_dft_real_forward_transform_stride (struct atfft_dft *plan,
                                              const atfft_sample *in,
                                              int in_stride,
                                              atfft_complex *out,
                                              int out_stride)
{
    /* Only to be used for forward real FFTs. */
    assert ((plan->format == ATFFT_REAL) && (plan->direction == ATFFT_FORWARD));

    atfft_sample_to_float_real_stride (in,
                                       in_stride,
                                       plan->in,
                                       1,
                                       plan->size);

    pffft_transform_ordered (plan->plan, plan->in, plan->out, plan->work_area, plan->pffft_direction);
    halfcomplex_pffft_to_atfft (plan->out, out, out_stride, plan->size);
}

static void halfcomplex_atfft_to_pffft (atfft_complex *in,
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

void atfft_dft_real_backward_transform (struct atfft_dft *plan, atfft_complex *in, atfft_sample *out)
{
    /* Only to be used for backward real FFTs. */
    assert ((plan->format == ATFFT_REAL) && (plan->direction == ATFFT_BACKWARD));

    halfcomplex_atfft_to_pffft (in, 1, plan->in, plan->size);
    pffft_transform_ordered (plan->plan, plan->in, plan->out, plan->work_area, plan->pffft_direction);

#ifdef ATFFT_TYPE_FLOAT
    memcpy (out, plan->out, plan->n_in_out_bytes);
#else
    atfft_float_to_sample_real (plan->out, out, plan->size);
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

    halfcomplex_atfft_to_pffft (in, in_stride, plan->in, plan->size);
    pffft_transform_ordered (plan->plan, plan->in, plan->out, plan->work_area, plan->pffft_direction);
    atfft_float_to_sample_real_stride (plan->out,
                                       1,
                                       out,
                                       out_stride,
                                       plan->size);
}
