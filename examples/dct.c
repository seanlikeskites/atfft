/*
 * Copyright (C) 2016 Sean Enderby <sean.enderby@gmail.com>
 *
 * This program is free software. It comes without any warranty, to
 * the extent permitted by applicable law. You can redistribute it 
 * and/or modify it under the terms of the Do What The Fuck You Want 
 * To Public License, Version 2, as published by Sam Hocevar. See 
 * the COPYING file for more details.
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <atfft/atfft.h>

#ifndef M_PI
#   define M_PI 3.14159265358979323846
#endif

void printSampleArray (atfft_sample *data, int size)
{
    int i = 0;

#ifdef ATFFT_TYPE_LONG_DOUBLE
    printf ("%Lf", data [i]);
#else
    printf ("%f", data [i]);
#endif

    for (i = 1; i < size; ++i)
    {
#ifdef ATFFT_TYPE_LONG_DOUBLE
        printf (", %Lf", data [i]);
#else
        printf (", %f", data [i]);
#endif
    }

    printf ("\n");
}

int main()
{
    int nSamples = 32;
    atfft_sample *signal, *transform;
    struct atfft_dct *dctForward, *dctBackward;
    int i = 0;

    /* allocate some memory for the signals */
    signal = malloc (nSamples * sizeof (*signal));
    transform = malloc (nSamples * sizeof (*transform));

    /* construct some signal */
    for (i = 0; i < nSamples; ++i)
    {
        atfft_sample x = 2.0 * M_PI * i / nSamples;

        signal [i] = 0.3 + 0.6 * cos (2.0 * x - 0.3)
                         + 0.3 * cos (5.0 * x + 0.2)
                         + 0.1 * cos (8.0 * x - 0.8);
    }

    printf ("Original Signal:\n");
    printSampleArray (signal, nSamples);

    /* create some dcts */
    dctForward = atfft_dct_create (nSamples, ATFFT_FORWARD);
    dctBackward = atfft_dct_create (nSamples, ATFFT_BACKWARD);

    /* apply the forward transform */
    atfft_dct_transform (dctForward, signal, transform);
    printf ("\nFrequency Domain:\n");
    printSampleArray (transform, nSamples);

    /* apply the backward transform */
    atfft_dct_transform (dctBackward, transform, signal);
    printf ("\nReconstructed Signal:\n");
    printSampleArray (signal, nSamples);

    /* normalise the output */
    atfft_normalise_real (signal, nSamples);
    printf ("\nNormalised Signal:\n");
    printSampleArray (signal, nSamples);

    /* free everything */
    atfft_dct_destroy (dctBackward);
    atfft_dct_destroy (dctForward);
    free (transform);
    free (signal);

    return 0;
}

