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
#include <atfft/atfft_dft_nd.h>
#include <atfft/atfft_dft.h>
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

    atfft_complex *work_area;
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

static int int_array_prod (const int *array, int size)
{
    int prod = array [0];

    for (int i = 1; i < size; ++i)
    {
        prod *= array [i];
    }

    return prod;
}

static int* init_strides (const int *dims, int n_dims, int data_size)
{
    int *strides = malloc (n_dims * sizeof (*strides));

    if (!strides)
        return NULL;

    for (int i = 0; i < n_dims; ++i)
    {
        strides [i] = data_size / dims [i];
    }

    return strides;
}

struct atfft_dft_nd* atfft_dft_nd_create (const int *dims,
                                          int n_dims,
                                          enum atfft_direction direction,
                                          enum atfft_format format)
{
    /* only use this for data with 2 dimensions or more */
    assert (n_dims > 1);

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
    fft->sub_transforms = atfft_init_sub_transforms (dims,
                                                     n_dims,
                                                     &(fft->n_sub_transforms),
                                                     fft->dim_sub_transforms,
                                                     direction,
                                                     format,
                                                     0);

    if (!fft->sub_transforms)
        goto failed;

    /* allocate work space */
    int data_size = int_array_prod (dims, n_dims);
    fft->work_area = malloc (data_size * sizeof (*(fft->work_area)));
    fft->strides = init_strides (dims, n_dims, data_size);

    if (!(fft->work_area && fft->strides))
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
        free (fft->strides);
        free (fft->work_area);
        atfft_free_sub_transforms (fft->sub_transforms, fft->n_sub_transforms);
        free (fft->dim_sub_transforms);
        free (fft->dims);
        free (fft);
    }
}

void atfft_dft_nd_complex_transform (struct atfft_dft_nd *fft, atfft_complex *in, atfft_complex *out)
{
}
