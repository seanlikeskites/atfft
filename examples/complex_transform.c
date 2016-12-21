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

void printComplexArray (atfft_complex *data, int size)
{
    int i = 0;

#ifdef ATFFT_TYPE_LONG_DOUBLE
    printf ("(%Lf, %Lf)", ATFFT_REAL (data [i]), ATFFT_IMAG (data [i]));
#else
    printf ("(%f, %f)", ATFFT_REAL (data [i]), ATFFT_IMAG (data [i]));
#endif

    for (i = 1; i < size; ++i)
    {
#ifdef ATFFT_TYPE_LONG_DOUBLE
        printf (", (%Lf, %Lf)", ATFFT_REAL (data [i]), ATFFT_IMAG (data [i]));
#else
        printf (", (%f, %f)", ATFFT_REAL (data [i]), ATFFT_IMAG (data [i]));
#endif
    }

    printf ("\n");
}

int main()
{
    int nSamples = 32;
    atfft_sample *signal;
    atfft_complex *timeDomain, *freqDomain;
    struct atfft_dft *fftForward, *fftBackward;
    int i = 0;

    /* allocate some memory for the signals */
    signal = malloc (nSamples * sizeof (*signal));
    timeDomain = malloc (nSamples * sizeof (*timeDomain));
    freqDomain = malloc (nSamples * sizeof (*freqDomain));

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

    /* copy the signal into an array of complex numbers */
    atfft_real_to_complex (signal, timeDomain, nSamples);

    /* create some ffts */
    fftForward = atfft_dft_create (nSamples, ATFFT_FORWARD, ATFFT_COMPLEX);
    fftBackward = atfft_dft_create (nSamples, ATFFT_BACKWARD, ATFFT_COMPLEX);

    /* apply the forward transform */
    atfft_dft_complex_transform (fftForward, timeDomain, freqDomain);
    printf ("\nFrequency Domain:\n");
    printComplexArray (freqDomain, nSamples);

    /* apply the backward transform */
    atfft_dft_complex_transform (fftBackward, freqDomain, timeDomain);
    atfft_real (timeDomain, signal, nSamples);
    printf ("\nReconstructed Signal:\n");
    printSampleArray (signal, nSamples);

    /* normalise the output */
    atfft_normalise_real (signal, nSamples);
    printf ("\nNormalised Signal:\n");
    printSampleArray (signal, nSamples);

    /* free everything */
    atfft_dft_destroy (fftBackward);
    atfft_dft_destroy (fftForward);
    free (freqDomain);
    free (timeDomain);
    free (signal);

    return 0;
}
