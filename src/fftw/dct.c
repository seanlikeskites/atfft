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
#include <atfft/dct.h>
#include <fftw3.h>
#include "fftw_definitions.h"

struct atfft_dct
{
    int size;
    enum atfft_direction direction;

    /* the fftw plan */
    atfft_fftw_plan plan;

    /* aligned input and output buffers for fftw transform */
    size_t n_in_out_bytes;
    atfft_sample *in, *out;
};

int atfft_dct_is_supported_size (int size)
{
    return size > 0;
}

struct atfft_dct* atfft_dct_create (int size, enum atfft_direction direction)
{
    /* fftw supports all sizes. */
    assert (atfft_dct_is_supported_size (size));

    struct atfft_dct *plan;

    if (!(plan = calloc (1, sizeof (*plan))))
        return NULL;

    plan->size = size;
    plan->direction = direction;

    /* allocate input and output buffers */
    plan->n_in_out_bytes = size * sizeof (*(plan->in));
    plan->in = ATFFT_FFTW_MALLOC (plan->n_in_out_bytes);
    plan->out = ATFFT_FFTW_MALLOC (plan->n_in_out_bytes);

    if (!(plan->in && plan->out))
        goto failed;

    /* initialise the fftw plan */
    if (direction == ATFFT_FORWARD)
    {
        plan->plan = ATFFT_FFTW_PLAN_R2R_1D (size,
                                             plan->in,
                                             plan->out,
                                             FFTW_REDFT10,
                                             ATFFT_FFTW_PLANNING_METHOD);
    }
    else
    {
        plan->plan = ATFFT_FFTW_PLAN_R2R_1D (size,
                                             plan->in,
                                             plan->out,
                                             FFTW_REDFT01,
                                             ATFFT_FFTW_PLANNING_METHOD);
    }

    if (!plan->plan)
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
        ATFFT_FFTW_DESTROY_PLAN (plan->plan);
        ATFFT_FFTW_FREE (plan->out);
        ATFFT_FFTW_FREE (plan->in);
        free (plan);
    }
}

void atfft_dct_transform (struct atfft_dct *plan, const atfft_sample *in, atfft_sample *out)
{
    memcpy (plan->in, in, plan->n_in_out_bytes);
    ATFFT_FFTW_EXECUTE (plan->plan);
    memcpy (out, plan->out, plan->n_in_out_bytes);

    /* fftw multiplies DCT bins by 2 */
    atfft_scale_real (out, plan->size, 0.5);
}
