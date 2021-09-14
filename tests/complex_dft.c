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
#include "test_functions.h"

int main()
{
    enum test_result ret = TEST_SUCCESS;
    int n_samples = 32;
    atfft_complex complex_one = {1.0, 0.0};
    atfft_sample threshold = 10e-10;

    /*************************************
     * Check the DFT of an impulse is DC
     *************************************/
    printf ("Checking DFT of impulse function:\n");
    atfft_complex *impulse = malloc (n_samples * sizeof (*impulse));
    atfft_complex *dc = malloc (n_samples * sizeof (*dc));

    if (!(impulse && dc))
    {
        ret = TEST_FAILURE;
        goto cleanup;
    }

    generate_complex_impulse (impulse, n_samples);
    generate_complex_dc (dc, n_samples, complex_one);

    ret = test_complex_dft (impulse, dc, n_samples, ATFFT_FORWARD, threshold);

    /*************************************
     * Check the iDFT of DC is an impulse
     *************************************/
    printf ("\nChecking iDFT of DC:\n");
    ATFFT_RE (impulse [0]) *= n_samples; /* scale impulse to correct magnitude */
    ret = test_complex_dft (dc, impulse, n_samples, ATFFT_BACKWARD, threshold);

cleanup:
    free (dc);
    free (impulse);

    return ret;
}
