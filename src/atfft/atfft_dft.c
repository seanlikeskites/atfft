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
#include <limits.h>
#include <atfft/atfft_dft.h>
#include "atfft_internal.h"
#include "atfft_dft_cooley_tukey.h"
#include "atfft_dft_rader.h"
#include "atfft_dft_bluestein.h"

typedef void (*complex_transform_function) (void*, atfft_complex*, atfft_complex*, int);
typedef void (*fft_destroy_function) (void*);

struct atfft_dft
{
    int size;
    enum atfft_direction direction;
    enum atfft_format format;

    void *fft;
    complex_transform_function complex_transform;
    fft_destroy_function fft_destroy;
};

struct atfft_dft* atfft_dft_create (int size, enum atfft_direction direction, enum atfft_format format)
{
    struct atfft_dft *fft;

    if (!(fft = calloc (1, sizeof (*fft))))
        return NULL;

    fft->size = size;
    fft->direction = direction;
    fft->format = format;

    if (atfft_is_prime (size))
    {
        if (atfft_is_power_of_2 (size - 1))
        {
            /* Use Rader's algorithm */
            fft->fft = atfft_dft_rader_create (size, direction, format);
            fft->complex_transform = (void*) atfft_dft_rader_complex_transform;
            fft->fft_destroy = (void*) atfft_dft_rader_destroy;
        }
        else
        {
            /* Use Bluestein's algorithm */
            fft->fft = atfft_dft_bluestein_create (size, direction, format);
            fft->complex_transform = (void*) atfft_dft_bluestein_complex_transform;
            fft->fft_destroy = (void*) atfft_dft_bluestein_destroy;
        }
    }
    else
    {
        /* Use Cooley-Tukey */
        fft->fft = atfft_dft_cooley_tukey_create (size, direction, format);
        fft->complex_transform = (void*) atfft_dft_cooley_tukey_complex_transform;
        fft->fft_destroy = (void*) atfft_dft_cooley_tukey_destroy;
    }

    if (!fft->fft)
    {
        atfft_dft_destroy (fft);
        return NULL;
    }

    return fft;
}

void atfft_dft_destroy (struct atfft_dft *fft)
{
    if (fft)
    {
        fft->fft_destroy (fft->fft);
        free (fft);
    }
}

void atfft_dft_complex_transform (struct atfft_dft *fft, atfft_complex *in, atfft_complex *out, int stride)
{
    /* Only to be used with complex FFTs. */
    assert (fft->format == ATFFT_COMPLEX);

    fft->complex_transform (fft->fft, in, out, stride);
}

//void atfft_dft_calculate_bin_real (struct atfft_dft *fft,
//                                   const atfft_sample *in,
//                                   atfft_complex *out,
//                                   int k)
//{
//    ATFFT_RE (*out) = 0.0;
//    ATFFT_IM (*out) = 0.0;
//
//    for (int i = 0; i < fft->size; ++i)
//    {
//        atfft_complex *t_factor = &(fft->t_factors [(k * i) % fft->size]);
//        ATFFT_RE (*out) += in [i] * ATFFT_RE (*t_factor);
//        ATFFT_IM (*out) += in [i] * ATFFT_IM (*t_factor);
//    }
//}
//
//void atfft_dft_real_forward_transform (struct atfft_dft *fft, const atfft_sample *in, atfft_complex *out)
//{
//    /* Only to be used for forward real FFTs. */
//    assert ((fft->format == ATFFT_REAL) && (fft->direction == ATFFT_FORWARD));
//
//    for (int k = 0; k < (fft->size / 2) + 1; ++k)
//    {
//        atfft_dft_calculate_bin_real (fft, in, out + k, k);
//    }
//}
//
//void atfft_dft_calculate_sample_real (struct atfft_dft *fft,
//                                      atfft_complex *in,
//                                      atfft_sample *out,
//                                      int i)
//{
//    int k = 0;
//    atfft_complex *t_factor = NULL;
//
//    /* get DC component */
//    *out = ATFFT_RE (in [0]);
//
//    /* add other frequency components */
//    for (k = 1; k < (fft->size / 2) + 1; ++k)
//    {
//        t_factor = &(fft->t_factors [(k * i) % fft->size]);
//        *out += 2.0 * ATFFT_RE (in [k]) * ATFFT_RE (*t_factor);
//        *out -= 2.0 * ATFFT_IM (in [k]) * ATFFT_IM (*t_factor);
//    }
//
//    if (atfft_is_even (fft->size))
//        *out -= ATFFT_RE (in [k - 1]) * ATFFT_RE (*t_factor);
//}
//
//void atfft_dft_real_backward_transform (struct atfft_dft *fft, atfft_complex *in, atfft_sample *out)
//{
//    /* Only to be used for backward real FFTs. */
//    assert ((fft->format == ATFFT_REAL) && (fft->direction == ATFFT_BACKWARD));
//
//    for (int i = 0; i < fft->size; ++i)
//    {
//        atfft_dft_calculate_sample_real (fft, in, out + i, i);
//    }
//}
