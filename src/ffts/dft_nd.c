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
#include <assert.h>
#include <atfft/dft_nd.h>
#include <atfft/dft_nd_util.h>
#include <ffts/ffts.h>

#ifndef ATFFT_TYPE_FLOAT
#   ifdef _MSC_VER
#       pragma message(": warning: FFTS only supports single precision floating point, " \
                       "higher precision values will be demoted to float for FFT calculations.")
#   else
#       warning FFTS only supports single precision floating point, \
                higher precision values will be demoted to float for FFT calculations.
#   endif
#endif

struct atfft_dft_nd
{
    size_t *dims;
    int n_dims;
    enum atfft_direction direction;
    enum atfft_format format;

    /* the ffts plan */
    ffts_plan_t *plan;

#ifndef ATFFT_TYPE_FLOAT
    /* input and output buffers for ffts transform */
    int in_size, out_size;
    float *in, *out;
#endif
};

static size_t* alloc_and_copy_dims_array (const int *dims,
                                          int n_dims)
{
    size_t *copy = malloc (n_dims * sizeof (*copy));

    if (!copy)
        return NULL;

    for (int i = 0; i < n_dims; ++i)
    {
        copy [i] = dims [i];
    }

    return copy;
}

static ffts_plan_t* init_ffts_plan (size_t *dims,
                                    int n_dims,
                                    enum atfft_direction direction,
                                    enum atfft_format format)
{
    int ffts_direction = 0;

    if (direction == ATFFT_FORWARD)
        ffts_direction = -1;
    else
        ffts_direction = 1;

    if (format == ATFFT_COMPLEX)
        return ffts_init_nd (n_dims, dims, ffts_direction);
    else
        return ffts_init_nd_real (n_dims, dims, ffts_direction);
}

struct atfft_dft_nd* atfft_dft_nd_create (const int *dims,
                                          int n_dims,
                                          enum atfft_direction direction,
                                          enum atfft_format format)
{
    /* FFTS only supports real transforms of even lengths */
    assert (format == ATFFT_COMPLEX || atfft_is_even (dims [n_dims - 1]));

    struct atfft_dft_nd *plan;

    if (!(plan = calloc (1, sizeof (*plan))))
        return NULL;

    plan->direction = direction;
    plan->format = format;
    plan->n_dims = n_dims;

    /* copy dims array */
    plan->dims = alloc_and_copy_dims_array (dims, n_dims);

    if (!plan->dims)
        goto failed;

#ifndef ATFFT_TYPE_FLOAT
    int full_size = atfft_int_array_product (dims, n_dims);

    if (format == ATFFT_COMPLEX)
    {
        plan->in_size = 2 * full_size;
        plan->out_size = 2 * full_size;
    }
    else
    {
        int halfcomplex_size = atfft_nd_halfcomplex_size (dims, n_dims);

        if (direction == ATFFT_FORWARD)
        {
            plan->in_size = full_size;
            plan->out_size = 2 * halfcomplex_size;
        }
        else
        {
            plan->in_size = 2 * halfcomplex_size;
            plan->out_size = full_size;
        }
    }

    plan->in = malloc (plan->in_size * sizeof (*(plan->in)));
    plan->out = malloc (plan->out_size * sizeof (*(plan->out)));

    if (!(plan->in && plan->out))
        goto failed;
#endif

    /* initialise the ffts plan */
    plan->plan = init_ffts_plan (plan->dims, plan->n_dims, direction, format);

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
        if (plan->plan) /* ffts can't hack freeing a null pointer */
            ffts_free (plan->plan);

#ifndef ATFFT_TYPE_FLOAT
        free (plan->out);
        free (plan->in);
#endif

        free (plan->dims);
        free (plan);
    }
}

static void apply_transform (struct atfft_dft_nd *plan, const atfft_sample *in, atfft_sample *out)
{
#ifdef ATFFT_TYPE_FLOAT
    ffts_execute (plan->plan, (const float*) in, (float*) out);
#else
    atfft_sample_to_float_real (in, plan->in, plan->in_size);
    ffts_execute (plan->plan, plan->in, plan->out);
    atfft_float_to_sample_real (plan->out, out, plan->out_size);
#endif
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
