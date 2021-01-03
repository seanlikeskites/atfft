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
#include <string.h>
#include <assert.h>
#include <atfft/dft_nd.h>
#include <atfft/dft_nd_util.h>
#include <atfft/dft.h>
#include <atfft/dft_util.h>
#include "atfft_internal.h"

struct atfft_dft_nd
{
    int *dims;
    int n_dims;
    enum atfft_direction direction;
    enum atfft_format format;

    /* plans for dimension sub-transforms */
    int n_sub_transforms;
    struct atfft_dft **sub_transforms;
    struct atfft_dft **dim_sub_transforms;

    /* additional plan for real transforms */
    struct atfft_dft *real_transform;

    atfft_complex *work_area, *real_backward_work_area;
    int *strides;
};

static void* alloc_and_copy_array (const void *array,
                                   size_t size)
{
    void *copy = malloc (size);

    if (!copy)
        return NULL;

    memcpy (copy, array, size);

    return copy;
}

static int* init_strides (const int *dims, int n_dims, int data_size, enum atfft_format format)
{
    int *strides = malloc (n_dims * sizeof (*strides));

    if (!strides)
        return NULL;

    int i = 0;

    for (; i < n_dims - 1; ++i)
    {
        strides [i] = data_size / dims [i];
    }

    if (format == ATFFT_REAL)
        strides [i] = data_size / atfft_dft_halfcomplex_size (dims [i]);
    else
        strides [i] = data_size / dims [i];

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
    fft->dims = alloc_and_copy_array (dims, n_dims * sizeof (*(fft->dims)));
    fft->dim_sub_transforms = calloc (n_dims, sizeof (*(fft->dim_sub_transforms)));

    if (!(fft->dims && fft->dim_sub_transforms))
        goto failed;

    /* allocate fft structs for each dimension */
    int n_complex_transforms = n_dims;
    int data_size = 0;

    if (format == ATFFT_REAL)
    {
        /* for real transforms we will use a 1d real transform
         * for the (n_dims - 1)th dimension */
        fft->real_transform = atfft_dft_create (dims [n_dims - 1],
                                                direction,
                                                ATFFT_REAL);

        n_complex_transforms = n_dims - 1;

        if (!fft->real_transform)
            goto failed;

        /* work space will be smaller for real transforms */
        data_size = atfft_dft_nd_halfcomplex_size (dims, n_dims);
    }
    else
    {
        data_size = atfft_int_array_product (dims, n_dims);
    }

    fft->sub_transforms = atfft_init_sub_transforms (dims,
                                                     n_complex_transforms,
                                                     &(fft->n_sub_transforms),
                                                     fft->dim_sub_transforms,
                                                     direction,
                                                     ATFFT_COMPLEX,
                                                     0);

    if (!fft->sub_transforms && n_complex_transforms > 0)
        goto failed;

    /* allocate work space */
    fft->work_area = malloc (data_size * sizeof (*(fft->work_area)));
    fft->strides = init_strides (dims, n_dims, data_size, format);

    if (!(fft->work_area && fft->strides))
        goto failed;

    /* for backward real transforms we need an additional working area,
     * to avoid overwriting the input signal */
    if (direction == ATFFT_BACKWARD && format == ATFFT_REAL)
    {
        fft->real_backward_work_area = malloc (data_size * sizeof (*(fft->work_area)));

        if (!fft->real_backward_work_area)
            goto failed;
    }

    return fft;

failed:
    atfft_dft_nd_destroy (fft);
    return NULL;
}

void atfft_dft_nd_destroy (struct atfft_dft_nd *fft)
{
    if (fft)
    {
        free (fft->strides);
        free (fft->real_backward_work_area);
        free (fft->work_area);
        atfft_dft_destroy (fft->real_transform);
        atfft_free_sub_transforms (fft->sub_transforms, fft->n_sub_transforms);
        free (fft->dim_sub_transforms);
        free (fft->dims);
        free (fft);
    }
}

static void complex_transform_and_transpose_left (struct atfft_dft *fft,
                                                  atfft_complex *in,
                                                  atfft_complex *out,
                                                  int size,
                                                  int stride)
{
    for (int i = 0; i < stride; ++i)
    {
        atfft_dft_complex_transform_stride (fft,
                                            in,
                                            1,
                                            out,
                                            stride);

        in += size;
        ++out;
    }
}

static void nd_complex_transform_left (int *dims,
                                       int *strides,
                                       struct atfft_dft **sub_transforms,
                                       int n_dims,
                                       atfft_complex *in,
                                       atfft_complex *work_area,
                                       atfft_complex *out)
{
    atfft_complex *work_areas[] = {work_area, out};
    int w = atfft_is_odd (n_dims);

    atfft_complex *current_in = in;
    atfft_complex *current_out = work_areas [w];

    for (int d = n_dims - 1; d >= 0; --d)
    {
        int size = dims [d];
        int stride = strides [d];
        struct atfft_dft *sub_transform = sub_transforms [d];

        complex_transform_and_transpose_left (sub_transform,
                                              current_in,
                                              current_out,
                                              size,
                                              stride);

        current_in = work_areas [w];
        w = 1 - w;
        current_out = work_areas [w];
    }
}

static void real_forward_transform_and_transpose_left (struct atfft_dft *fft,
                                                       const atfft_sample *in,
                                                       atfft_complex *out,
                                                       int size,
                                                       int stride)
{
    for (int i = 0; i < stride; ++i)
    {
        atfft_dft_real_forward_transform_stride (fft,
                                                 in,
                                                 1,
                                                 out,
                                                 stride);

        in += size;
        ++out;
    }
}

static void complex_transform_and_transpose_right (struct atfft_dft *fft,
                                                   atfft_complex *in,
                                                   atfft_complex *out,
                                                   int size,
                                                   int stride)
{
    for (int i = 0; i < stride; ++i)
    {
        atfft_dft_complex_transform_stride (fft,
                                            in,
                                            stride,
                                            out,
                                            1);

        ++in;
        out += size;
    }
}

static void nd_complex_transform_right (int *dims,
                                        int *strides,
                                        struct atfft_dft **sub_transforms,
                                        int n_dims,
                                        atfft_complex *in,
                                        atfft_complex *work_area,
                                        atfft_complex *out)
{
    atfft_complex *work_areas[] = {work_area, out};
    int w = atfft_is_odd (n_dims);

    atfft_complex *current_in = in;
    atfft_complex *current_out = work_areas [w];

    for (int d = 0; d < n_dims; ++d)
    {
        int size = dims [d];
        int stride = strides [d];
        struct atfft_dft *sub_transform = sub_transforms [d];

        complex_transform_and_transpose_right (sub_transform,
                                               current_in,
                                               current_out,
                                               size,
                                               stride);

        current_in = work_areas [w];
        w = 1 - w;
        current_out = work_areas [w];
    }
}

static void real_backward_transform_and_transpose_right (struct atfft_dft *fft,
                                                         atfft_complex *in,
                                                         atfft_sample *out,
                                                         int size,
                                                         int stride)
{
    for (int i = 0; i < stride; ++i)
    {
        atfft_dft_real_backward_transform_stride (fft,
                                                  in,
                                                  stride,
                                                  out,
                                                  1);

        ++in;
        out += size;
    }
}

void atfft_dft_nd_complex_transform (struct atfft_dft_nd *fft, atfft_complex *in, atfft_complex *out)
{
    /* Only to be used with complex FFTs. */
    assert (fft->format == ATFFT_COMPLEX);

    nd_complex_transform_right (fft->dims,
                                fft->strides,
                                fft->dim_sub_transforms,
                                fft->n_dims,
                                in,
                                fft->work_area,
                                out);
}

void atfft_dft_nd_real_forward_transform (struct atfft_dft_nd *fft, const atfft_sample *in, atfft_complex *out)
{
    /* Only to be used for forward real FFTs. */
    assert ((fft->format == ATFFT_REAL) && (fft->direction == ATFFT_FORWARD));

    /* perform a real transform on the last dimension */
    atfft_complex *real_transform_out = atfft_is_odd (fft->n_dims) ? out : fft->work_area;
    int last_dim = fft->n_dims - 1;

    real_forward_transform_and_transpose_left (fft->real_transform,
                                               in,
                                               real_transform_out,
                                               fft->dims [last_dim],
                                               fft->strides [last_dim]);

    /* do complex transforms for the remaining dimensions */
    nd_complex_transform_left (fft->dims,
                               fft->strides,
                               fft->dim_sub_transforms,
                               fft->n_dims - 1,
                               real_transform_out,
                               fft->work_area,
                               out);
}

void atfft_dft_nd_real_backward_transform (struct atfft_dft_nd *fft, atfft_complex *in, atfft_sample *out)
{
    /* Only to be used for backward real FFTs. */
    assert ((fft->format == ATFFT_REAL) && (fft->direction == ATFFT_BACKWARD));

    /* do complex transforms on the first n_dims - 1 dimensions */
    nd_complex_transform_right (fft->dims,
                                fft->strides,
                                fft->dim_sub_transforms,
                                fft->n_dims - 1,
                                in,
                                fft->work_area,
                                fft->real_backward_work_area);

    /* finally, perform the real transform on the last dimension */
    atfft_complex *real_transform_in = fft->real_backward_work_area;
    int last_dim = fft->n_dims - 1;

    real_backward_transform_and_transpose_right (fft->real_transform,
                                                 real_transform_in,
                                                 out,
                                                 fft->dims [last_dim],
                                                 fft->strides [last_dim]);
}
