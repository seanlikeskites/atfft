/*
 * Copyright (c) 2021 Sean Enderby <sean.enderby@gmail.com>
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
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
    int dims[] = {4, 4, 4};
    int nDims = 3;
    int nSamples = atfft_int_array_product (dims, nDims);

    atfft_complex *in = malloc (nSamples * sizeof (*in));
    atfft_complex *out = malloc (nSamples * sizeof (*out));

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
    atfft_normalise_dft_complex (in, nSamples);
    printf ("\nReconstructed Signal:\n");
    printComplexArray (in, nSamples);

    atfft_dft_nd_destroy (ifft);
    atfft_dft_nd_destroy (fft);
    free (out);
    free (in);

    return 0;
}
