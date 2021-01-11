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

#define KFR_NO_C_COMPLEX_TYPES 1

#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <atfft/dft.h>
#include <atfft/dft_util.h>
#include <kfr/capi.h>
#include "kfr_definitions.h"

typedef void (*kfr_execute_function) (void*, atfft_kfr_sample*, const atfft_kfr_sample*, uint8_t*);

struct atfft_dft
{
    int size;
    enum atfft_direction direction;
    enum atfft_format format;

    /* the kfr plan */
    void *plan;
    kfr_execute_function transform_function;

    /* buffer for krf internals */
    uint8_t *work_area;

    /* input and output buffers for kfr transform */
    size_t in_size, out_size;
    atfft_kfr_sample *in, *out;
};

int atfft_dft_is_supported_size (int size, enum atfft_format format)
{
    int max = 16777216;

    if (format == ATFFT_COMPLEX)
        return size >= 2 && size <= max;
    else
        return size >= 4 && size <= max && atfft_is_even(size);
}

static void* init_kfr_plan (int size, enum atfft_format format)
{
    if (format == ATFFT_COMPLEX)
        return ATFFT_KFR_DFT_CREATE_PLAN (size);
    else
        return ATFFT_KFR_DFT_REAL_CREATE_PLAN (size, CCs);
}

static kfr_execute_function get_kfr_execute_function (enum atfft_direction direction,
                                                      enum atfft_format format)
{
    if (format == ATFFT_COMPLEX)
    {
        if (direction == ATFFT_FORWARD)
            return (kfr_execute_function) ATFFT_KFR_DFT_EXECUTE;
        else
            return (kfr_execute_function) ATFFT_KFR_DFT_EXECUTE_INVERSE;
    }
    else
    {
        if (direction == ATFFT_FORWARD)
            return (kfr_execute_function) ATFFT_KFR_DFT_REAL_EXECUTE;
        else
            return (kfr_execute_function) ATFFT_KFR_DFT_REAL_EXECUTE_INVERSE;
    }
}

struct atfft_dft* atfft_dft_create (int size, enum atfft_direction direction, enum atfft_format format)
{
    /* kfr most sizes for complex transforms, but real transforms must be even. */
    assert (atfft_dft_is_supported_size (size, format));

    struct atfft_dft *plan;

    if (!(plan = calloc (1, sizeof (*plan))))
        return NULL;

    plan->size = size;
    plan->direction = direction;
    plan->format = format;

    /* initialise the kfr plan */
    plan->plan = init_kfr_plan (size, format);

    if (!plan->plan)
        goto failed;

    plan->transform_function = get_kfr_execute_function (direction, format);

    /* once we have a plan we can find out how much work space it needs */
    size_t n_work_bytes = 0;

    if (format == ATFFT_COMPLEX)
        n_work_bytes = ATFFT_KFR_DFT_GET_TEMP_SIZE (plan->plan);
    else
        n_work_bytes = ATFFT_KFR_DFT_REAL_GET_TEMP_SIZE (plan->plan);

    plan->work_area = kfr_allocate (n_work_bytes);

    if (!plan->work_area)
        goto failed;

    /* allocate input and output buffers */
    if (format == ATFFT_COMPLEX)
    {
        plan->in_size = 2 * size;
        plan->out_size = plan->in_size;
    }
    else
    {
        if (direction == ATFFT_FORWARD)
        {
            plan->in_size = size;
            plan->out_size = 2 * atfft_halfcomplex_size (size);
        }
        else
        {
            plan->in_size = 2 * atfft_halfcomplex_size (size);
            plan->out_size = size;
        }
    }

    plan->in = kfr_allocate (plan->in_size * sizeof (*(plan->in)));
    plan->out = kfr_allocate (plan->out_size * sizeof (*(plan->out)));

    if (!(plan->in && plan->out))
        goto failed;


    return plan;

failed:
    atfft_dft_destroy (plan);
    return NULL;
}

void atfft_dft_destroy (struct atfft_dft *plan)
{
    if (plan)
    {
        kfr_deallocate (plan->out);
        kfr_deallocate (plan->in);
        kfr_deallocate (plan->work_area);

        if (plan->format == ATFFT_COMPLEX)
            ATFFT_KFR_DFT_DELETE_PLAN (plan->plan);
        else
            ATFFT_KFR_DFT_REAL_DELETE_PLAN (plan->plan);

        free (plan);
    }
}

static void apply_transform (struct atfft_dft *plan, const atfft_sample *in, atfft_sample *out)
{
#ifdef ATFFT_TYPE_LONG_DOUBLE
    atfft_sample_to_double_real (in, plan->in, plan->in_size);
    plan->transform_function (plan->plan, plan->out, plan->in, plan->work_area);
    atfft_double_to_sample_real (plan->out, out, plan->out_size);
#else
    plan->transform_function (plan->plan, out, in, plan->work_area);
#endif
}

void atfft_dft_complex_transform (struct atfft_dft *plan, atfft_complex *in, atfft_complex *out)
{
    /* Only to be used with complex FFTs. */
    assert (plan->format == ATFFT_COMPLEX);

    apply_transform (plan, (atfft_sample*) in, (atfft_sample*) out);
}

void atfft_dft_complex_transform_stride (struct atfft_dft *plan,
                                         atfft_complex *in,
                                         int in_stride,
                                         atfft_complex *out,
                                         int out_stride)
{
    /* Only to be used with complex FFTs. */
    assert (plan->format == ATFFT_COMPLEX);

#ifdef ATFFT_TYPE_FLOAT
    atfft_sample_to_float_complex_stride (in,
                                          in_stride,
                                          (atfft_complex*) plan->in,
                                          1,
                                          plan->in_size / 2);
#else
    atfft_sample_to_double_complex_stride (in,
                                           in_stride,
                                           (atfft_complex_d*) plan->in,
                                           1,
                                           plan->in_size / 2);
#endif

    plan->transform_function (plan->plan, plan->out, plan->in, plan->work_area);

#ifdef ATFFT_TYPE_FLOAT
    atfft_float_to_sample_complex_stride ((atfft_complex*) plan->out,
                                          1,
                                          out,
                                          out_stride,
                                          plan->out_size / 2);
#else
    atfft_double_to_sample_complex_stride ((atfft_complex_d*) plan->out,
                                           1,
                                           out,
                                           out_stride,
                                           plan->out_size / 2);
#endif
}

void atfft_dft_real_forward_transform (struct atfft_dft *plan, const atfft_sample *in, atfft_complex *out)
{
    /* Only to be used for forward real FFTs. */
    assert ((plan->format == ATFFT_REAL) && (plan->direction == ATFFT_FORWARD));

    apply_transform (plan, in, (atfft_sample*) out);
}

void atfft_dft_real_forward_transform_stride (struct atfft_dft *plan,
                                              const atfft_sample *in,
                                              int in_stride,
                                              atfft_complex *out,
                                              int out_stride)
{
    /* Only to be used for forward real FFTs. */
    assert ((plan->format == ATFFT_REAL) && (plan->direction == ATFFT_FORWARD));

#ifdef ATFFT_TYPE_FLOAT
    atfft_sample_to_float_real_stride (in,
                                       in_stride,
                                       plan->in,
                                       1,
                                       plan->in_size);
#else
    atfft_sample_to_double_real_stride (in,
                                        in_stride,
                                        plan->in,
                                        1,
                                        plan->in_size);
#endif

    plan->transform_function (plan->plan, plan->out, plan->in, plan->work_area);

#ifdef ATFFT_TYPE_FLOAT
    atfft_float_to_sample_complex_stride ((atfft_complex*) plan->out,
                                          1,
                                          out,
                                          out_stride,
                                          plan->out_size / 2);
#else
    atfft_double_to_sample_complex_stride ((atfft_complex_d*) plan->out,
                                           1,
                                           out,
                                           out_stride,
                                           plan->out_size / 2);
#endif
}

void atfft_dft_real_backward_transform (struct atfft_dft *plan, atfft_complex *in, atfft_sample *out)
{
    /* Only to be used for backward real FFTs. */
    assert ((plan->format == ATFFT_REAL) && (plan->direction == ATFFT_BACKWARD));

    apply_transform (plan, (atfft_sample*) in, out);
}

void atfft_dft_real_backward_transform_stride (struct atfft_dft *plan,
                                               atfft_complex *in,
                                               int in_stride,
                                               atfft_sample *out,
                                               int out_stride)
{
    /* Only to be used for backward real FFTs. */
    assert ((plan->format == ATFFT_REAL) && (plan->direction == ATFFT_BACKWARD));

#ifdef ATFFT_TYPE_FLOAT
    atfft_sample_to_float_complex_stride (in,
                                          in_stride,
                                          (atfft_complex*) plan->in,
                                          1,
                                          plan->in_size / 2);
#else
    atfft_sample_to_double_complex_stride (in,
                                           in_stride,
                                           (atfft_complex_d*) plan->in,
                                           1,
                                           plan->in_size / 2);
#endif

    plan->transform_function (plan->plan, plan->out, plan->in, plan->work_area);

#ifdef ATFFT_TYPE_FLOAT
    atfft_float_to_sample_real_stride (plan->out,
                                       1,
                                       out,
                                       out_stride,
                                       plan->out_size);
#else
    atfft_double_to_sample_real_stride (plan->out,
                                        1,
                                        out,
                                        out_stride,
                                        plan->out_size);
#endif
}
