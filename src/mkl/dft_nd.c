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
#include <atfft/dft_util.h>
#include "mkl_definitions.h"

struct atfft_dft_nd
{
    enum atfft_direction direction;
    enum atfft_format format;

    MKL_LONG *dims;
    MKL_LONG n_dims;

    /* the mkl plan */
    DFTI_DESCRIPTOR_HANDLE plan;

#ifdef ATFFT_TYPE_LONG_DOUBLE
    /* input and output buffers for mkl transform */
    int in_size, out_size;
    double *in, *out;
#endif
};

static MKL_LONG* alloc_and_copy_dims_array (const int *dims,
                                            int n_dims)
{
    MKL_LONG *copy = malloc (n_dims * sizeof (*copy));

    if (!copy)
        return NULL;

    for (int i = 0; i < n_dims; ++i)
    {
        copy [i] = dims [i];
    }

    return copy;
}

static MKL_LONG* init_halfcomplex_strides (const MKL_LONG *dims,
                                           MKL_LONG n_dims)
{
    MKL_LONG *strides = malloc ((n_dims + 1) * sizeof (*strides));

    if (!strides)
        return NULL;

    strides [0] = 0;
    strides [n_dims - 1] = atfft_halfcomplex_size (dims [n_dims - 1]);
    strides [n_dims] = 1;

    for (int i = n_dims - 2; i > 0; --i)
    {
        strides [i] = strides [i + 1] * dims [i];
    }

    return strides;
}

static int set_mkl_strides (DFTI_DESCRIPTOR_HANDLE plan,
                            const MKL_LONG *dims,
                            MKL_LONG n_dims,
                            enum atfft_direction direction)
{
    int ret = 0;
    MKL_LONG *halfcomplex_strides = init_halfcomplex_strides (dims, n_dims);

    if (!halfcomplex_strides)
        return -1;

    MKL_LONG status = DFTI_NO_ERROR;

    if (direction == ATFFT_FORWARD)
    {
        status = DftiSetValue(plan,
                              DFTI_OUTPUT_STRIDES,
                              halfcomplex_strides);

    }
    else
    {
        status = DftiSetValue(plan,
                              DFTI_INPUT_STRIDES,
                              halfcomplex_strides);

    }

    if (status != DFTI_NO_ERROR)
        ret = -1;

    free (halfcomplex_strides);
    return ret;
}

static DFTI_DESCRIPTOR_HANDLE init_mkl_plan (const MKL_LONG *dims,
                                             MKL_LONG n_dims,
                                             enum atfft_direction direction,
                                             enum atfft_format format)
{
    MKL_LONG status = DFTI_NO_ERROR;
    DFTI_DESCRIPTOR_HANDLE plan;

    /* create the plan */
    if (format == ATFFT_COMPLEX)
        status = DftiCreateDescriptor (&plan,
                                       ATFFT_MKL_PRECISION,
                                       DFTI_COMPLEX,
                                       n_dims,
                                       dims);
    else
        status = DftiCreateDescriptor (&plan,
                                       ATFFT_MKL_PRECISION,
                                       DFTI_REAL,
                                       n_dims,
                                       dims);

    if (status != DFTI_NO_ERROR)
        goto failed;

    /* set strides for real transforms */
    if (format == ATFFT_REAL)
    {
        if (set_mkl_strides (plan, dims, n_dims, direction) < 0)
            goto failed;
    }

    /* we aren't doing in place transforms */
    if (DftiSetValue(plan,
                     DFTI_PLACEMENT,
                     DFTI_NOT_INPLACE) != DFTI_NO_ERROR)
        goto failed;

    /* for real transforms we want an unpacked frequency domain */
    if (DftiSetValue(plan,
                     DFTI_CONJUGATE_EVEN_STORAGE,
                     DFTI_COMPLEX_COMPLEX) != DFTI_NO_ERROR)
        goto failed;
    
    /* commit the plan settings */
    if (DftiCommitDescriptor(plan) != DFTI_NO_ERROR)
        goto failed;

    return plan;

failed:
    DftiFreeDescriptor (&plan);
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
    plan->n_dims = n_dims;

    /* copy dims array */
    plan->dims = alloc_and_copy_dims_array (dims, n_dims);

    if (!plan->dims)
        goto failed;

    /* initialise the mkl plan */
    plan->plan = init_mkl_plan (plan->dims, plan->n_dims, direction, format);

    if (!plan->plan)
        goto failed;

#ifdef ATFFT_TYPE_LONG_DOUBLE
    /* allocate input and output buffers */
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

    return plan;

failed:
    atfft_dft_nd_destroy (plan);
    return NULL;
}

void atfft_dft_nd_destroy (struct atfft_dft_nd *plan)
{
    if (plan)
    {
#ifdef ATFFT_TYPE_LONG_DOUBLE
        free (plan->out);
        free (plan->in);
#endif

        DftiFreeDescriptor (&(plan->plan));
        free (plan->dims);
        free (plan);
    }
}

void atfft_dft_nd_complex_transform (struct atfft_dft_nd *plan, atfft_complex *in, atfft_complex *out)
{
    /* Only to be used with complex FFTs. */
    assert (plan->format == ATFFT_COMPLEX);

#ifdef ATFFT_TYPE_LONG_DOUBLE
    atfft_sample_to_double_real ((atfft_sample*) in, plan->in, plan->in_size);

    if (plan->direction == ATFFT_FORWARD)
        DftiComputeForward(plan->plan, plan->in, plan->out);
    else
        DftiComputeBackward(plan->plan, plan->in, plan->out);

    atfft_double_to_sample_real (plan->out, (atfft_sample*) out, plan->out_size);
#else
    if (plan->direction == ATFFT_FORWARD)
        DftiComputeForward(plan->plan, (atfft_sample*) in, (atfft_sample*) out);
    else
        DftiComputeBackward(plan->plan, (atfft_sample*) in, (atfft_sample*) out);
#endif
}

void atfft_dft_nd_real_forward_transform (struct atfft_dft_nd *plan, const atfft_sample *in, atfft_complex *out)
{
    /* Only to be used for forward real FFTs. */
    assert ((plan->format == ATFFT_REAL) && (plan->direction == ATFFT_FORWARD));

#ifdef ATFFT_TYPE_LONG_DOUBLE
    atfft_sample_to_double_real ((atfft_sample*) in, plan->in, plan->in_size);
    DftiComputeForward(plan->plan, plan->in, plan->out);
    atfft_double_to_sample_real (plan->out, (atfft_sample*) out, plan->out_size);
#else
    DftiComputeForward(plan->plan, (atfft_sample*) in, (atfft_sample*) out);
#endif
}

void atfft_dft_nd_real_backward_transform (struct atfft_dft_nd *plan, atfft_complex *in, atfft_sample *out)
{
    /* Only to be used for backward real FFTs. */
    assert ((plan->format == ATFFT_REAL) && (plan->direction == ATFFT_BACKWARD));

#ifdef ATFFT_TYPE_LONG_DOUBLE
    atfft_sample_to_double_real ((atfft_sample*) in, plan->in, plan->in_size);
    DftiComputeBackward(plan->plan, plan->in, plan->out);
    atfft_double_to_sample_real (plan->out, (atfft_sample*) out, plan->out_size);
#else
    DftiComputeBackward(plan->plan, (atfft_sample*) in, (atfft_sample*) out);
#endif
}
