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
#include "mkl_definitions.h"
#include <atfft/dft_nd.h>
#include <atfft/dft_nd_util.h>
#include <atfft/dft_util.h>

struct atfft_dft_nd
{
    enum atfft_direction direction;
    enum atfft_format format;

    MKL_LONG *dims;
    MKL_LONG n_dims;
    MKL_LONG *halfcomplex_strides;

    DFTI_DESCRIPTOR_HANDLE plan;

    /* buffers to work in if working with long double */
#ifdef ATFFT_TYPE_LONG_DOUBLE
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

static MKL_LONG* init_halfcomplex_strides (const int *dims,
                                          int n_dims)
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

struct atfft_dft_nd* atfft_dft_nd_create (const int *dims,
                                          int n_dims,
                                          enum atfft_direction direction,
                                          enum atfft_format format)
{
    struct atfft_dft_nd *fft;

    if (!(fft = calloc (1, sizeof (*fft))))
        return NULL;

    fft->direction = direction;
    fft->format = format;
    fft->n_dims = n_dims;

    /* copy dims array */
    fft->dims = alloc_and_copy_dims_array (dims, n_dims);

    if (!fft->dims)
        goto failed;

    if (format == ATFFT_REAL)
    {
        /* init halfcomplex strides array */
        fft->halfcomplex_strides = init_halfcomplex_strides (dims, n_dims);

        if (!fft->halfcomplex_strides)
            goto failed;
    }

#ifdef ATFFT_TYPE_LONG_DOUBLE
    int full_size = atfft_int_array_product (dims, n_dims);
    int halfcomplex_size = atfft_nd_halfcomplex_size (dims, n_dims);
#endif

    /* create the plan */
    MKL_LONG status = DFTI_NO_ERROR;

    switch (format)
    {
        case ATFFT_COMPLEX:
            status = DftiCreateDescriptor (&(fft->plan),
                                           ATFFT_MKL_PRECISION,
                                           DFTI_COMPLEX,
                                           fft->n_dims,
                                           fft->dims);

#ifdef ATFFT_TYPE_LONG_DOUBLE
            fft->in_size = 2 * full_size;
            fft->out_size = 2 * full_size;
#endif
            break;

        case ATFFT_REAL:
            status = DftiCreateDescriptor (&(fft->plan),
                                           ATFFT_MKL_PRECISION,
                                           DFTI_REAL,
                                           fft->n_dims,
                                           fft->dims);

            if (direction == ATFFT_FORWARD)
            {
                status = DftiSetValue(fft->plan,
                                      DFTI_OUTPUT_STRIDES,
                                      fft->halfcomplex_strides);

#ifdef ATFFT_TYPE_LONG_DOUBLE
                fft->in_size = full_size;
                fft->out_size = 2 * halfcomplex_size;
#endif
            }
            else
            {
                status = DftiSetValue(fft->plan,
                                      DFTI_INPUT_STRIDES,
                                      fft->halfcomplex_strides);

#ifdef ATFFT_TYPE_LONG_DOUBLE
                fft->in_size = 2 * halfcomplex_size;
                fft->out_size = full_size;
#endif
            }

            break;
    }

#ifdef ATFFT_TYPE_LONG_DOUBLE
    fft->in = malloc (fft->in_size * sizeof (*(fft->in)));
    fft->out = malloc (fft->out_size * sizeof (*(fft->out)));
#endif

#ifdef ATFFT_TYPE_LONG_DOUBLE
    if (!(fft->in && fft->out) && status != DFTI_NO_ERROR)
#else
    if (status != DFTI_NO_ERROR)
#endif
        goto failed;

    /* set some other parameters on the plan */
    if (DftiSetValue(fft->plan,
                     DFTI_PLACEMENT,
                     DFTI_NOT_INPLACE) != DFTI_NO_ERROR)
        goto failed;

    if (DftiSetValue(fft->plan,
                     DFTI_CONJUGATE_EVEN_STORAGE,
                     DFTI_COMPLEX_COMPLEX) != DFTI_NO_ERROR)
        goto failed;

    /* commit the plan settings */
    status = DftiCommitDescriptor(fft->plan);

    if (status != DFTI_NO_ERROR)
        goto failed;

    return fft;

failed:
    atfft_dft_nd_destroy (fft);
    return NULL;
}

void atfft_dft_nd_destroy (struct atfft_dft_nd *fft)
{
    if (fft)
    {
#ifdef ATFFT_TYPE_LONG_DOUBLE
        free (fft->out);
        free (fft->in);
#endif

        DftiFreeDescriptor (&(fft->plan));
        free (fft->halfcomplex_strides);
        free (fft->dims);
        free (fft);
    }
}

void atfft_dft_nd_complex_transform (struct atfft_dft_nd *fft, atfft_complex *in, atfft_complex *out)
{
    /* Only to be used with complex FFTs. */
    assert (fft->format == ATFFT_COMPLEX);

#ifdef ATFFT_TYPE_LONG_DOUBLE
    atfft_sample_to_double_real ((atfft_sample*) in, fft->in, fft->in_size);

    if (fft->direction == ATFFT_FORWARD)
        DftiComputeForward(fft->plan, fft->in, fft->out);
    else
        DftiComputeBackward(fft->plan, fft->in, fft->out);

    atfft_double_to_sample_real (fft->out, (atfft_sample*) out, fft->out_size);
#else
    if (fft->direction == ATFFT_FORWARD)
        DftiComputeForward(fft->plan, (atfft_sample*) in, (atfft_sample*) out);
    else
        DftiComputeBackward(fft->plan, (atfft_sample*) in, (atfft_sample*) out);
#endif
}

void atfft_dft_nd_real_forward_transform (struct atfft_dft_nd *fft, const atfft_sample *in, atfft_complex *out)
{
    /* Only to be used for forward real FFTs. */
    assert ((fft->format == ATFFT_REAL) && (fft->direction == ATFFT_FORWARD));

#ifdef ATFFT_TYPE_LONG_DOUBLE
    atfft_sample_to_double_real ((atfft_sample*) in, fft->in, fft->in_size);
    DftiComputeForward(fft->plan, fft->in, fft->out);
    atfft_double_to_sample_real (fft->out, (atfft_sample*) out, fft->out_size);
#else
    DftiComputeForward(fft->plan, (atfft_sample*) in, (atfft_sample*) out);
#endif
}

void atfft_dft_nd_real_backward_transform (struct atfft_dft_nd *fft, atfft_complex *in, atfft_sample *out)
{
    /* Only to be used for backward real FFTs. */
    assert ((fft->format == ATFFT_REAL) && (fft->direction == ATFFT_BACKWARD));

#ifdef ATFFT_TYPE_LONG_DOUBLE
    atfft_sample_to_double_real ((atfft_sample*) in, fft->in, fft->in_size);
    DftiComputeBackward(fft->plan, fft->in, fft->out);
    atfft_double_to_sample_real (fft->out, (atfft_sample*) out, fft->out_size);
#else
    DftiComputeBackward(fft->plan, (atfft_sample*) in, (atfft_sample*) out);
#endif
}
