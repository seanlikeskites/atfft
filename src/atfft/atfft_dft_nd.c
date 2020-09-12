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
#include <atfft/atfft_dft_nd.h>
#include <atfft/atfft_dft.h>

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

    int data_size;
    atfft_complex *work_area;
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

    /* copy dims array */
    fft->dims = alloc_and_copy_array (dims, n_dims * sizeof (*dims));

    if (!fft->dims)
        goto failed;

    fft->n_dims = n_dims;

    /* allocate fft structs for each dimension */

    /* allocate work space */
    

    return fft;

failed:
    atfft_dft_nd_destroy (fft);
    return NULL;
}

void atfft_dft_nd_destroy (struct atfft_dft_nd *fft)
{
    if (fft)
    {
        free (fft->dims);
        free (fft);
    }
}


void atfft_dft_nd_complex_transform (struct atfft_dft_nd *fft, atfft_complex *in, atfft_complex *out)
{
}
