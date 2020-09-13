/*
 * Copyright (C) 2020 Sean Enderby <sean.enderby@gmail.com>
 *
 * This program is free software. It comes without any warranty, to
 * the extent permitted by applicable law. You can redistribute it 
 * and/or modify it under the terms of the Do What The Fuck You Want 
 * To Public License, Version 2, as published by Sam Hocevar. See 
 * the COPYING file for more details.
 */

#include "atfft_internal.h"
#include <atfft/atfft_dft.h>
#include <math.h>
#include <stdlib.h>

/******************************************
 * External definitions of inline functions.
 ******************************************/
extern inline void atfft_copy_complex (const atfft_complex x, atfft_complex *y);
extern inline void atfft_swap_complex (const atfft_complex x, atfft_complex *y);
extern inline void atfft_sum_complex (const atfft_complex a,
                                      const atfft_complex b,
                                      atfft_complex *s);
extern inline void atfft_difference_complex (const atfft_complex a,
                                             const atfft_complex b,
                                             atfft_complex *d);
extern inline void atfft_product_complex (const atfft_complex a,
                                          const atfft_complex b,
                                          atfft_complex *p);
extern inline void atfft_multiply_by_complex (atfft_complex *a, const atfft_complex b);
extern inline void atfft_multiply_by_and_swap_complex (atfft_complex *a, const atfft_complex b);

/******************************************
 * Other functions.
 ******************************************/
int atfft_next_power_of_2 (int x)
{
    if (x <= 0)
        return 0;
    else
        return pow (2, (int) log2 (x) + 1);
}

void atfft_twiddle_factor (int n,
                           int N,
                           enum atfft_direction d,
                           atfft_complex *t)
{
    atfft_scaled_twiddle_factor (n, N, d, 1.0, t);
}

void atfft_scaled_twiddle_factor (int n,
                                  int N,
                                  enum atfft_direction d,
                                  atfft_sample s,
                                  atfft_complex *t)
{
    atfft_sample x = 2.0 * n * M_PI / N;
    ATFFT_RE (*t) = cos (x) / s;
    ATFFT_IM (*t) = sin (x) / s;

    if (d == ATFFT_FORWARD)
        ATFFT_IM (*t) *= -1.0;
}

int atfft_is_prime (int x)
{
    if (x <= 1 || (atfft_is_even (x) && x > 2))
        return 0;

    int sqrt_x = (int) sqrt ((double) x);

    for (int i = 2; i <= sqrt_x; ++i)
    {
        if ((x % i) == 0)
            return 0;
    }

    return 1;
}

/******************************************
 * Functions for allocating sub-transform
 * plans for use on transforms of large
 * prime sizes.
 ******************************************/
static int atfft_integer_is_in_array (const int *arr, int size, int member)
{
    for (int i = 0; i < size; ++i)
    {
        if (arr [i] == member)
            return 1;
    }

    return 0;
}

static int atfft_sub_transform_sizes (const int *sizes,
                                      int *sub_transform_sizes,
                                      int size,
                                      int threshold)
{
    int n_sub_transforms = 0;

    for (int i = 0; i < size; ++i)
    {
        if (sizes [i] > threshold && 
            !atfft_integer_is_in_array (sub_transform_sizes, n_sub_transforms, sizes [i]))
            sub_transform_sizes [n_sub_transforms++] = sizes [i];
    }

    return n_sub_transforms;
}

static struct atfft_dft* atfft_populate_sub_transform (int size,
                                                       const int *sizes,
                                                       int n_sizes,
                                                       struct atfft_dft **size_sub_transforms,
                                                       enum atfft_direction direction,
                                                       enum atfft_format format)
{
    struct atfft_dft *fft = atfft_dft_create (size, direction, format);

    if (!fft)
        return NULL;

    for (int i = 0; i < n_sizes; ++i)
    {
        if (sizes [i] == size)
            size_sub_transforms [i] = fft;
    }

    return fft;
}

struct atfft_dft** atfft_init_sub_transforms (const int *sizes,
                                              int n_sizes,
                                              int *n_sub_transforms,
                                              struct atfft_dft **size_sub_transforms,
                                              enum atfft_direction direction,
                                              enum atfft_format format,
                                              int threshold)
{
    struct atfft_dft **sub_transforms = NULL;
    int *sub_transform_sizes = malloc (n_sizes * sizeof (*sub_transform_sizes));

    if (!sub_transform_sizes)
        return NULL;

    /* get unique sizes for which sub-transforms are required */
    *n_sub_transforms = atfft_sub_transform_sizes (sizes, sub_transform_sizes, n_sizes, threshold);

    if (*n_sub_transforms == 0)
        goto succeeded; /* we don't need to allocate any sub-transforms */

    /* allocate some space for sub-transforms */
    sub_transforms = calloc (*n_sub_transforms, sizeof (*sub_transforms));

    if (!sub_transforms)
        goto failed;

    /* create the sub-transform dft structs */
    for (int i = 0; i < *n_sub_transforms; ++i)
    {
        struct atfft_dft *sub_transform = atfft_populate_sub_transform (sub_transform_sizes [i],
                                                                        sizes,
                                                                        n_sizes,
                                                                        size_sub_transforms,
                                                                        direction,
                                                                        format);

        if (!sub_transform)
            goto failed;

        sub_transforms [i] = sub_transform;
    }

    goto succeeded;

failed:
    atfft_free_sub_transforms (sub_transforms, *n_sub_transforms);
    sub_transforms = NULL;

succeeded:
    free (sub_transform_sizes);
    return sub_transforms;
}

void atfft_free_sub_transforms (struct atfft_dft **sub_transforms,
                                int size)
{
    if (sub_transforms)
    {
        for (int i = 0; i < size; ++i)
        {
            atfft_dft_destroy (sub_transforms [i]);
        }

        free (sub_transforms);
    }
}
