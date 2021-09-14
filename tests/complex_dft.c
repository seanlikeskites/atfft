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

#include <atfft/atfft.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "test_functions.h"

int main()
{
    enum test_result ret = TEST_SUCCESS;
    enum test_result step = TEST_SUCCESS;
    int n_samples = 32;
    atfft_sample threshold = 10e-10;

    atfft_complex *input = malloc (n_samples * sizeof (*input));
    atfft_complex *expected = malloc (n_samples * sizeof (*expected));

    if (!(input && expected))
    {
        ret = TEST_FAILURE;
        goto cleanup;
    }

    /*************************************
     * Check the DFT of an impulse is DC
     *************************************/
    printf ("Checking DFT of impulse function:\n");

    atfft_complex complex_1 = {1.0, 0.0};
    generate_complex_impulse (input, n_samples, complex_1);
    generate_complex_dc (expected, n_samples, complex_1);

    step = test_complex_dft (input, expected, n_samples, ATFFT_FORWARD, threshold);

    if (step != TEST_SUCCESS)
        ret = TEST_FAILURE;

    /*************************************
     * Check the iDFT of DC is an impulse
     *************************************/
    printf ("\nChecking iDFT of DC:\n");

    atfft_complex complex_32 = {32.0, 0.0};
    generate_complex_dc (input, n_samples, complex_1);
    generate_complex_impulse (expected, n_samples, complex_32);

    step = test_complex_dft (input, expected, n_samples, ATFFT_BACKWARD, threshold);

    if (step != TEST_SUCCESS)
        ret = TEST_FAILURE;

    /*************************************
     * Check the DFT of a cosine
     *************************************/
    printf ("\nChecking iDFT of cosine wave:\n");

    atfft_complex complex_0 = {0.0, 0.0};
    generate_complex_cosine (input, n_samples, 5.0, 1.0, 0);
    generate_complex_dc (expected, n_samples, complex_0);
    ATFFT_RE (expected [5]) = 16.0;
    ATFFT_RE (expected [n_samples - 5]) = 16.0;

    step = test_complex_dft (input, expected, n_samples, ATFFT_FORWARD, threshold);

    if (step != TEST_SUCCESS)
        ret = TEST_FAILURE;

    /*************************************
     * Check the DFT of a sine
     *************************************/
    printf ("\nChecking iDFT of sine wave:\n");

    generate_complex_cosine (input, n_samples, 5.0, 1.0, - M_PI / 2.0);
    generate_complex_dc (expected, n_samples, complex_0);
    ATFFT_IM (expected [5]) = -16.0;
    ATFFT_IM (expected [n_samples - 5]) = 16.0;

    step = test_complex_dft (input, expected, n_samples, ATFFT_FORWARD, threshold);

    if (step != TEST_SUCCESS)
        ret = TEST_FAILURE;

cleanup:
    free (expected);
    free (input);

    return ret;
}
