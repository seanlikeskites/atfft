#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <atfft.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

void printDoubleArray (double *data, int size)
{
    int i = 0;

    printf ("%f", data [i]);

    for (i = 1; i < size; ++i)
    {
        printf (", %f", data [i]);
    }

    printf ("\n");
}

void printComplexArray (double *data, int size)
{
    int i = 0;

    printf ("(%f, %f)", data [i], data [i + 1]);

    for (i = 1; i < size; ++i)
    {
        printf (", (%f, %f)", data [2 * i], data [2 * i + 1]);
    }

    printf ("\n");
}

int main()
{
    int nSamples = 16;
    double *signal, *timeDomain, *freqDomain;
    struct atfft *fftForward, *fftBackward;
    int i = 0;

    /* allocate some memory for the signals */
    signal = malloc (nSamples * sizeof (*signal));
    timeDomain = malloc (2 * nSamples * sizeof (*timeDomain));
    freqDomain = malloc (2 * nSamples * sizeof (*freqDomain));

    /* construct some signal */
    for (i = 0; i < nSamples; ++i)
    {
        double x = 2.0 * M_PI * i / nSamples;

        signal [i] = 0.3 + 0.6 * cos (2.0 * x - 0.3)
                         + 0.3 * cos (5.0 * x + 0.2)
                         + 0.1 * cos (8.0 * x - 0.8);
    }

    printf ("Original Signal:\n");
    printDoubleArray (signal, nSamples);

    /* copy the signal into an array of complex numbers */
    atfft_real_to_complex (signal, timeDomain, nSamples);

    /* create some ffts */
    fftForward = atfft_create (nSamples, ATFFT_FORWARD, ATFFT_COMPLEX);
    fftBackward = atfft_create (nSamples, ATFFT_BACKWARD, ATFFT_COMPLEX);

    /* apply the forward transform */
    atfft_complex_transform (fftForward, timeDomain, freqDomain);
    printf ("\nFrequency Domain:\n");
    printComplexArray (freqDomain, nSamples);

    /* apply the backward transform */
    atfft_complex_transform (fftBackward, freqDomain, timeDomain);
    atfft_real (timeDomain, signal, nSamples);
    printf ("\nReconstructed Signal:\n");
    printDoubleArray (signal, nSamples);

    /* normalise the output */
    atfft_normalise_real (signal, nSamples);
    printf ("\nNormalised Signal:\n");
    printDoubleArray (signal, nSamples);

    /* free everything */
    atfft_free (fftBackward);
    atfft_free (fftForward);
    free (freqDomain);
    free (timeDomain);
    free (signal);

    return 0;
}
