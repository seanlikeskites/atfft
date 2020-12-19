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
    printf ("(%Lf, %Lf)", ATFFT_RE (data [i]), ATFFT_IM (data [i]));
#else
    printf ("(%f, %f)", ATFFT_RE (data [i]), ATFFT_IM (data [i]));
#endif

    for (i = 1; i < size; ++i)
    {
#ifdef ATFFT_TYPE_LONG_DOUBLE
        printf (", (%Lf, %Lf)", ATFFT_RE (data [i]), ATFFT_IM (data [i]));
#else
        printf (", (%f, %f)", ATFFT_RE (data [i]), ATFFT_IM (data [i]));
#endif
    }

    printf ("\n");
}

int main()
{
    int dims[] = {2, 2, 4};
    int nDims = 3;
    int nSamples = 16;

    atfft_complex *in = malloc (nSamples * sizeof (*in));
    atfft_complex *out = malloc (nSamples * sizeof (*in));

    for (int i = 0; i < nSamples; ++i)
    {
        ATFFT_RE (in [i]) = i;
        ATFFT_IM (in [i]) = nSamples - i;
    }

    printf ("Original Signal:\n");
    printComplexArray (in, nSamples);

    struct atfft_dft_nd *fft = atfft_dft_nd_create (dims, nDims, ATFFT_FORWARD, ATFFT_COMPLEX);
    struct atfft_dft_nd *ifft = atfft_dft_nd_create (dims, nDims, ATFFT_BACKWARD, ATFFT_COMPLEX);

    atfft_dft_nd_complex_transform (fft, in, out);
    printf ("\nFrequency Domain:\n");
    printComplexArray (out, nSamples);

    atfft_dft_nd_complex_transform (ifft, out, in);
    atfft_normalise_complex (in, nSamples);
    printf ("\nReconstructed Signal:\n");
    printComplexArray (in, nSamples);

    atfft_dft_nd_destroy (ifft);
    atfft_dft_nd_destroy (fft);
    free (out);
    free (in);

    return 0;
}

