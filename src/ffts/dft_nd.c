/*
 * Copyright (C) 2020 Sean Enderby <sean.enderby@gmail.com>
 *
 * This program is free software. It comes without any warranty, to
 * the extent permitted by applicable law. You can redistribute it 
 * and/or modify it under the terms of the Do What The Fuck You Want 
 * To Public License, Version 2, as published by Sam Hocevar. See 
 * the COPYING file for more details.
 */

#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <ffts/ffts.h>
#include <atfft/atfft_dft_nd.h>

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
#ifndef ATFFT_TYPE_FLOAT
    int in_size, out_size;
    float *in, *out;
#endif
    ffts_plan_t *plan;
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

struct atfft_dft_nd* atfft_dft_nd_create (const int *dims,
                                          int n_dims,
                                          enum atfft_direction direction,
                                          enum atfft_format format)
{
    /* FFTS only supports real transforms of even lengths */
    //assert (format == ATFFT_COMPLEX || atfft_is_even (size));

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

    int ffts_direction = 0;

    if (direction == ATFFT_FORWARD)
        ffts_direction = -1;
    else
        ffts_direction = 1;

#ifndef ATFFT_TYPE_FLOAT
    int full_size = atfft_int_array_product (dims, n_dims);
    int halfcomplex_size = atfft_dft_nd_halfcomplex_size (dims, n_dims);
#endif

    switch (format)
    {
        case ATFFT_COMPLEX:
#ifndef ATFFT_TYPE_FLOAT
            fft->in_size = 2 * full_size;
            fft->out_size = 2 * full_size;
#endif
            fft->plan = ffts_init_nd (fft->n_dims, fft->dims, ffts_direction);
            break;

        case ATFFT_REAL:
#ifndef ATFFT_TYPE_FLOAT
            if (direction == ATFFT_FORWARD)
            {
                fft->in_size = full_size;
                fft->out_size = 2 * halfcomplex_size;
            }
            else
            {
                fft->in_size = 2 * halfcomplex_size;
                fft->out_size = full_size;
            }
#endif

            fft->plan = ffts_init_nd_real (fft->n_dims, fft->dims, ffts_direction);
            break;
    }

#ifndef ATFFT_TYPE_FLOAT
    fft->in = malloc (fft->in_size * sizeof (*(fft->in)));
    fft->out = malloc (fft->out_size * sizeof (*(fft->out)));
#endif

#ifndef ATFFT_TYPE_FLOAT
    if (!(fft->in && fft->out && fft->plan))
#else
    if (!fft->plan)
#endif
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
        if (fft->plan) /* ffts can't hack freeing a null pointer */
            ffts_free (fft->plan);

#ifndef ATFFT_TYPE_FLOAT
        free (fft->out);
        free (fft->in);
#endif

        free (fft->dims);
        free (fft);
    }
}

void atfft_dft_nd_complex_transform (struct atfft_dft_nd *fft, atfft_complex *in, atfft_complex *out)
{
    /* Only to be used with complex FFTs. */
    assert (fft->format == ATFFT_COMPLEX);

#ifdef ATFFT_TYPE_FLOAT
    ffts_execute (fft->plan, (const float*) in, (float*) out);
#else
    atfft_sample_to_float_real ((atfft_sample*) in, fft->in, fft->in_size);
    ffts_execute (fft->plan, fft->in, fft->out);
    atfft_float_to_sample_real (fft->out, (atfft_sample*) out, fft->out_size);
#endif
}

void atfft_dft_nd_real_forward_transform (struct atfft_dft_nd *fft, const atfft_sample *in, atfft_complex *out)
{
    /* Only to be used for forward real FFTs. */
    assert ((fft->format == ATFFT_REAL) && (fft->direction == ATFFT_FORWARD));

#ifdef ATFFT_TYPE_FLOAT
    ffts_execute (fft->plan, (const float*) in, (float*) out);
#else
    atfft_sample_to_float_real (in, fft->in, fft->in_size);
    ffts_execute (fft->plan, fft->in, fft->out);
    atfft_float_to_sample_real (fft->out, (atfft_sample*) out, fft->out_size);
#endif
}

void atfft_dft_nd_real_backward_transform (struct atfft_dft_nd *fft, atfft_complex *in, atfft_sample *out)
{
    /* Only to be used for backward real FFTs. */
    assert ((fft->format == ATFFT_REAL) && (fft->direction == ATFFT_BACKWARD));

#ifdef ATFFT_TYPE_FLOAT
    ffts_execute (fft->plan, (const float*) in, (float*) out);
#else
    atfft_sample_to_float_real ((atfft_sample*) in, fft->in, fft->in_size);
    ffts_execute (fft->plan, fft->in, fft->out);
    atfft_float_to_sample_real (fft->out, out, fft->out_size);
#endif
}
