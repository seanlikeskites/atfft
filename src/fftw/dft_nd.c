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
#include <atfft/dft_nd.h>
#include <atfft/dft_nd_util.h>
#include <fftw3.h>
#include "fftw_definitions.h"

struct atfft_dft_nd
{
    enum atfft_direction direction;
    enum atfft_format format;

    /* the fftw plan */
    atfft_fftw_plan plan;

    /* aligned input and output buffers for fftw transform */
    size_t n_in_bytes, n_out_bytes;
    atfft_sample *in, *out;
};

static atfft_fftw_plan init_fftw_plan (const int *dims,
                                       int n_dims,
                                       enum atfft_direction direction,
                                       enum atfft_format format,
                                       atfft_sample *in,
                                       atfft_sample *out)
{
    if (format == ATFFT_COMPLEX)
    {
        if (direction == ATFFT_FORWARD)
            return ATFFT_FFTW_PLAN_DFT (n_dims,
                                        dims,
                                        (atfft_fftw_complex*) in,
                                        (atfft_fftw_complex*) out,
                                        FFTW_FORWARD,
                                        ATFFT_FFTW_PLANNING_METHOD);
        else
            return ATFFT_FFTW_PLAN_DFT (n_dims,
                                        dims,
                                        (atfft_fftw_complex*) in,
                                        (atfft_fftw_complex*) out,
                                        FFTW_BACKWARD,
                                        ATFFT_FFTW_PLANNING_METHOD);
    }
    else
    {
        if (direction == ATFFT_FORWARD)
            return ATFFT_FFTW_PLAN_DFT_R2C (n_dims,
                                            dims,
                                            in,
                                            (atfft_fftw_complex*) out,
                                            ATFFT_FFTW_PLANNING_METHOD);
        else
            return ATFFT_FFTW_PLAN_DFT_C2R (n_dims,
                                            dims,
                                            (atfft_fftw_complex*) in,
                                            out,
                                            ATFFT_FFTW_PLANNING_METHOD);
    }

    return NULL;
}

struct atfft_dft_nd* atfft_dft_nd_create (const int *dims,
                                          int n_dims,
                                          enum atfft_direction direction,
                                          enum atfft_format format)
{
    struct atfft_dft_nd *plan;

    if (!(plan = calloc (1, sizeof (*plan))))
        return NULL;

    plan->direction = direction;
    plan->format = format;

    int full_size = atfft_int_array_product (dims, n_dims);

    if (format == ATFFT_COMPLEX)
    {
        plan->n_in_bytes = 2 * full_size * sizeof (*(plan->in));
        plan->n_out_bytes = 2 * full_size * sizeof (*(plan->out));
    }
    else
    {
        int halfcomplex_size = atfft_nd_halfcomplex_size (dims, n_dims);

        if (direction == ATFFT_FORWARD)
        {
            plan->n_in_bytes = full_size * sizeof (*(plan->in));
            plan->n_out_bytes = 2 * halfcomplex_size * sizeof (*(plan->out));
        }
        else
        {
            plan->n_in_bytes = 2 * halfcomplex_size * sizeof (*(plan->in));
            plan->n_out_bytes = full_size * sizeof (*(plan->out));
        }
    }

    /* allocate input and output buffers */
    plan->in = ATFFT_FFTW_MALLOC (plan->n_in_bytes);
    plan->out = ATFFT_FFTW_MALLOC (plan->n_out_bytes);

    if (!(plan->in && plan->out))
        goto failed;

    /* initialise the fftw plan */
    plan->plan = init_fftw_plan (dims, n_dims, direction, format, plan->in, plan->out);

    if (!plan->plan)
        goto failed;

    return plan;

failed:
    atfft_dft_nd_destroy (plan);
    return NULL;
}

void atfft_dft_nd_destroy (struct atfft_dft_nd *plan)
{
    if (plan)
    {
        ATFFT_FFTW_DESTROY_PLAN (plan->plan);
        ATFFT_FFTW_FREE (plan->out);
        ATFFT_FFTW_FREE (plan->in);
        free (plan);
    }
}

static void apply_transform (struct atfft_dft_nd *plan, const atfft_sample *in, atfft_sample *out)
{
    memcpy (plan->in, in, plan->n_in_bytes);
    ATFFT_FFTW_EXECUTE (plan->plan);
    memcpy (out, plan->out, plan->n_out_bytes);
}

void atfft_dft_nd_complex_transform (struct atfft_dft_nd *plan, atfft_complex *in, atfft_complex *out)
{
    /* Only to be used with complex FFTs. */
    assert (plan->format == ATFFT_COMPLEX);

    apply_transform (plan, (atfft_sample*) in, (atfft_sample*) out);
}

void atfft_dft_nd_real_forward_transform (struct atfft_dft_nd *plan, const atfft_sample *in, atfft_complex *out)
{
    /* Only to be used for forward real FFTs. */
    assert ((plan->format == ATFFT_REAL) && (plan->direction == ATFFT_FORWARD));

    apply_transform (plan, in, (atfft_sample*) out);
}

void atfft_dft_nd_real_backward_transform (struct atfft_dft_nd *plan, atfft_complex *in, atfft_sample *out)
{
    /* Only to be used for backward real FFTs. */
    assert ((plan->format == ATFFT_REAL) && (plan->direction == ATFFT_BACKWARD));

    apply_transform (plan, (atfft_sample*) in, out);
}
