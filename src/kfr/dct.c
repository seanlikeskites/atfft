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
#include <atfft/dct.h>
#include <kfr/capi.h>
#include "kfr_definitions.h"

typedef void (*kfr_execute_function) (void*, atfft_kfr_sample*, const atfft_kfr_sample*, uint8_t*);

struct atfft_dct
{
    int size;
    enum atfft_direction direction;

    /* the kfr plan */
    void *plan;
    kfr_execute_function transform_function;

    /* buffer for krf internals */
    uint8_t *work_area;

    /* input and output buffers for kfr transform */
    size_t in_size, out_size;
    atfft_kfr_sample *in, *out;
};

int atfft_dct_is_supported_size (int size)
{
    return size >= 4 && size <= 16777216;
}

struct atfft_dct* atfft_dct_create (int size, enum atfft_direction direction)
{
    /* kfr most sizes */
    assert (atfft_dct_is_supported_size (size));

    struct atfft_dct *plan;

    if (!(plan = calloc (1, sizeof (*plan))))
        return NULL;

    plan->size = size;
    plan->direction = direction;

    /* initialise the kfr plan */
    plan->plan = ATFFT_KFR_DCT_CREATE_PLAN (size);

    if (!plan->plan)
        goto failed;

    if (direction == ATFFT_FORWARD)
        plan->transform_function = (kfr_execute_function) ATFFT_KFR_DCT_EXECUTE;
    else
        plan->transform_function = (kfr_execute_function) ATFFT_KFR_DCT_EXECUTE_INVERSE;

    /* once we have a plan we can find out how much work space it needs */
    size_t n_work_bytes = 0;
    n_work_bytes = ATFFT_KFR_DCT_GET_TEMP_SIZE (plan->plan);
    plan->work_area = kfr_allocate (n_work_bytes);

    if (!plan->work_area)
        goto failed;

    /* allocate input and output buffers */
    plan->in_size = size;
    plan->out_size = size;
    plan->in = kfr_allocate (plan->in_size * sizeof (*(plan->in)));
    plan->out = kfr_allocate (plan->out_size * sizeof (*(plan->out)));

    if (!(plan->in && plan->out))
        goto failed;

    return plan;

failed:
    atfft_dct_destroy (plan);
    return NULL;
}

void atfft_dct_destroy (struct atfft_dct *plan)
{
    if (plan)
    {
        kfr_deallocate (plan->out);
        kfr_deallocate (plan->in);
        kfr_deallocate (plan->work_area);
        ATFFT_KFR_DCT_DELETE_PLAN (plan->plan);
        free (plan);
    }
}

void atfft_dct_transform (struct atfft_dct *plan, const atfft_sample *in, atfft_sample *out)
{
#ifdef ATFFT_TYPE_LONG_DOUBLE
    atfft_sample_to_double_real (in, plan->in, plan->in_size);
    plan->transform_function (plan->plan, plan->out, plan->in, plan->work_area);
    atfft_double_to_sample_real (plan->out, out, plan->out_size);
#else
    plan->transform_function (plan->plan, out, in, plan->work_area);
#endif

    if (plan->direction == ATFFT_BACKWARD)
        atfft_scale_real (out, plan->size, 2.0);
}
