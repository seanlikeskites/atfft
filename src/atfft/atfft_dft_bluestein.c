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
#include <atfft/atfft_dft.h>
#include "atfft_internal.h"
#include "atfft_dft_bluestein.h"

/**
 * Swap the real and imaginary parts of a and then multiply by b:
 *
 * p = j * conj(a) * b
 */
static void atfft_swap_and_product_complex (const atfft_complex a,
                                            const atfft_complex b,
                                            atfft_complex *p)
{
    ATFFT_RE (*p) = ATFFT_IM (a) * ATFFT_RE (b) -
                    ATFFT_RE (a) * ATFFT_IM (b);
    ATFFT_IM (*p) = ATFFT_IM (a) * ATFFT_IM (b) +
                    ATFFT_RE (a) * ATFFT_RE (b);
}


struct atfft_dft_bluestein
{
    int size;
    enum atfft_direction direction;
    enum atfft_format format;
    int conv_size;
    struct atfft_dft *fft;
    atfft_complex *sig, *sig_dft, *conv, *conv_dft, *factors;
};

static int atfft_bluestein_convolution_fft_size (int size)
{
    if (atfft_is_power_of_2 (size))
        return size;
    else
        return atfft_next_power_of_2 (2 * size - 1);
}

static int atfft_init_bluestein_convolution_dft (int size,
                                                 enum atfft_direction direction,
                                                 atfft_complex *conv_dft,
                                                 int conv_size,
                                                 atfft_complex *factors,
                                                 struct atfft_dft *fft)
{
    atfft_complex *sequence = calloc (conv_size, sizeof (*sequence));

    if (!sequence)
        return -1;

    /* produce convolution sequence */
    for (int i = 0; i < size; ++i)
    {
        atfft_twiddle_factor (-i * i, 2 * size, direction, sequence + i);
    }

    /* replicate samples for circular convolution */
    if (conv_size > size)
    {
        for (int i = 1; i < size; ++i)
        {
            atfft_copy_complex (sequence [i], &sequence [conv_size - i]);
        }
    }

    /* take DFT of sequence */
    atfft_dft_complex_transform (fft, sequence, conv_dft);
    atfft_normalise_complex (conv_dft, conv_size);

    /* take conjugate of the sequence for use later */
    for (int i = 0; i < size; ++i)
    {
        ATFFT_RE (factors [i]) = ATFFT_RE (sequence [i]);
        ATFFT_IM (factors [i]) = - ATFFT_IM (sequence [i]);
    }

    free (sequence);
    return 0;
}

struct atfft_dft_bluestein* atfft_dft_bluestein_create (int size,
                                                        enum atfft_direction direction,
                                                        enum atfft_format format)
{
    struct atfft_dft_bluestein *fft;

    if (!(fft = calloc (1, sizeof (*fft))))
        return NULL;

    fft->size = size;
    fft->direction = direction;
    fft->format = format;

    /* allocate some regular dft objects for performing the convolution */
    fft->conv_size = atfft_bluestein_convolution_fft_size (size);
    fft->fft = atfft_dft_create (fft->conv_size, ATFFT_FORWARD, ATFFT_COMPLEX);

    if (!fft->fft)
        goto failed;

    /* allocate work space for performing the dft */
    fft->sig = calloc (fft->conv_size, sizeof (*(fft->sig))); /* set up for zero padding */
    fft->sig_dft = malloc (fft->conv_size * sizeof (*(fft->sig_dft)));
    fft->conv = malloc (fft->conv_size * sizeof (*(fft->conv)));
    fft->conv_dft = malloc (fft->conv_size * sizeof (*(fft->conv_dft)));
    fft->factors = malloc (size * sizeof (*(fft->factors)));

    if (!(fft->sig && fft->sig_dft && fft->conv && fft->conv_dft && fft->factors))
        goto failed;

    /* calculate the convolution dft */
    if (atfft_init_bluestein_convolution_dft (size,
                                              direction,
                                              fft->conv_dft,
                                              fft->conv_size,
                                              fft->factors,
                                              fft->fft) < 0)
        goto failed;

    return fft;

failed:
    atfft_dft_bluestein_destroy (fft);
    return NULL;
}

void atfft_dft_bluestein_destroy (struct atfft_dft_bluestein *fft)
{
    if (fft)
    {
        free (fft->factors);
        free (fft->conv_dft);
        free (fft->conv);
        free (fft->sig_dft);
        free (fft->sig);
        atfft_dft_destroy (fft->fft);
        free (fft);
    }
}

void atfft_dft_bluestein_complex_transform (struct atfft_dft_bluestein *fft,
                                            atfft_complex *in,
                                            int in_stride,
                                            atfft_complex *out,
                                            int out_stride)
{
    /* multiply the input signal with the factors */
    for (int i = 0; i < fft->size; ++i)
    {
        atfft_product_complex (in [i * in_stride], fft->factors [i], fft->sig + i);
    }

    /* take DFT of the result */
    atfft_dft_complex_transform (fft->fft, fft->sig, fft->sig_dft);

    /* perform convolution in the frequency domain */
    for (int i = 0; i < fft->conv_size; ++i)
    {
        atfft_multiply_by_and_swap_complex (fft->sig_dft + i, fft->conv_dft [i]);
    }

    /* take the inverse DFT of the result */
    atfft_dft_complex_transform (fft->fft, fft->sig_dft, fft->conv);

    /* multiply the output transform with the factors */
    for (int i = 0; i < fft->size; ++i)
    {
        atfft_swap_and_product_complex (fft->conv [i], fft->factors [i], out + i * out_stride);
    }
}
