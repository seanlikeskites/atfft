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
#include <atfft/dft.h>
#include <atfft/dft_util.h>
#include "mkl_definitions.h"

struct atfft_dft
{
    int size;
    enum atfft_direction direction;
    enum atfft_format format;

    /* the mkl plan */
    DFTI_DESCRIPTOR_HANDLE plan;

#ifdef ATFFT_TYPE_LONG_DOUBLE
    /* input and output buffers for mkl transform */
    int in_size, out_size;
    double *in, *out;
#endif
};

int atfft_dft_is_supported_size (int size, enum atfft_format format)
{
    return size > 0;
}

static DFTI_DESCRIPTOR_HANDLE init_mkl_plan (int size,
                                             enum atfft_format format)
{
    MKL_LONG status = DFTI_NO_ERROR;
    DFTI_DESCRIPTOR_HANDLE plan;

    /* create the plan */
    if (format == ATFFT_COMPLEX)
        status = DftiCreateDescriptor (&plan,
                                       ATFFT_MKL_PRECISION,
                                       DFTI_COMPLEX,
                                       1,
                                       size);
    else
        status = DftiCreateDescriptor (&plan,
                                       ATFFT_MKL_PRECISION,
                                       DFTI_REAL,
                                       1,
                                       size);

    if (status != DFTI_NO_ERROR)
        goto failed;

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

struct atfft_dft* atfft_dft_create (int size, enum atfft_direction direction, enum atfft_format format)
{
    /* mkl supports all sizes. */
    assert (atfft_dft_is_supported_size (size, format));

    struct atfft_dft *plan;

    if (!(plan = calloc (1, sizeof (*plan))))
        return NULL;

    plan->size = size;
    plan->direction = direction;
    plan->format = format;

    /* initialise the mkl plan */
    plan->plan = init_mkl_plan (size, format);

    if (!plan->plan)
        goto failed;

#ifdef ATFFT_TYPE_LONG_DOUBLE
    /* allocate input and output buffers */
    if (format == ATFFT_COMPLEX)
    {
        plan->in_size = 2 * size;
        plan->out_size = 2 * size;
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

    plan->in = malloc (plan->in_size * sizeof (*(plan->in)));
    plan->out = malloc (plan->out_size * sizeof (*(plan->out)));

    if (!(plan->in && plan->out))
        goto failed;
#endif

    return plan;

failed:
    atfft_dft_destroy (plan);
    return NULL;
}

void atfft_dft_destroy (struct atfft_dft *plan)
{
    if (plan)
    {
#ifdef ATFFT_TYPE_LONG_DOUBLE
        free (plan->out);
        free (plan->in);
#endif

        DftiFreeDescriptor (&(plan->plan));
        free (plan);
    }
}

void atfft_dft_complex_transform (struct atfft_dft *plan, atfft_complex *in, atfft_complex *out)
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

void atfft_dft_real_forward_transform (struct atfft_dft *plan, const atfft_sample *in, atfft_complex *out)
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

void atfft_dft_real_backward_transform (struct atfft_dft *plan, atfft_complex *in, atfft_sample *out)
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
