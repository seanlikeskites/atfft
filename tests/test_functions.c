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

#include "test_functions.h"
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

void generate_real_dc (atfft_sample *sig, int size, atfft_sample amplitude)
{
    for (int i = 0; i < size; ++i)
    {
        sig [i] = amplitude;
    }
}

void generate_complex_dc (atfft_complex *sig, int size, atfft_complex amplitude)
{
    for (int i = 0; i < size; ++i)
    {
        ATFFT_RE (sig [i]) = ATFFT_RE (amplitude);
        ATFFT_IM (sig [i]) = ATFFT_IM (amplitude);
    }
}

void generate_real_impulse (atfft_sample *sig, int size, atfft_sample amplitude)
{
    sig [0] = amplitude;

    for (int i = 1; i < size; ++i)
    {
        sig [i] = 0.0;
    }
}

void generate_complex_impulse (atfft_complex *sig, int size, atfft_complex amplitude)
{
    ATFFT_RE (sig [0]) = ATFFT_RE (amplitude);
    ATFFT_IM (sig [0]) = ATFFT_IM (amplitude);

    for (int i = 1; i < size; ++i)
    {
        ATFFT_RE (sig [i]) = 0.0;
        ATFFT_IM (sig [i]) = 0.0;
    }
}

void generate_complex_cosine (atfft_complex *sig,
                              int size,
                              atfft_sample frequency,
                              atfft_sample amplitude,
                              atfft_sample phase)
{
    atfft_sample inc = 2 * M_PI * frequency / size;
    atfft_sample p = 0.0;

    for (int i = 0; i < size; ++i)
    {
        ATFFT_RE (sig [i]) = amplitude * cos (p + phase);
        ATFFT_IM (sig [i]) = 0.0;

        p = fmod (p + inc, 2 * M_PI);
    }
}

atfft_sample max_error_real (const atfft_sample *a, const atfft_sample *b, int size)
{
    atfft_sample error = 0.0;

    for (int i = 0; i < size; ++i)
    {
        atfft_sample e = fabs (a [i] - b [i]);

        if (e > error)
            error = e;
    }

    return error;
}

atfft_sample max_error_complex (atfft_complex *a, atfft_complex *b, int size)
{
    return max_error_real ((atfft_sample*) a, (atfft_sample*) b, size * 2);
}

enum test_result test_complex_dft (atfft_complex *in,
                                   atfft_complex *expected,
                                   int size,
                                   enum atfft_direction direction,
                                   atfft_sample threshold)
{
    enum test_result ret = TEST_SUCCESS;

    /* Allocate memory for output */
    atfft_complex *out = malloc (size * sizeof (*out));

    if (!out)
    {
        ret = TEST_FAILURE;
        goto cleanup;
    }

    /* make a DFT plan */
    struct atfft_dft *fft = atfft_dft_create (size, direction, ATFFT_COMPLEX);

    if (!fft)
    {
        ret = TEST_FAILURE;
        goto cleanup;
    }

    /* perform transform */
    atfft_dft_complex_transform (fft, in, out);

    /* verify output */
    atfft_sample error = max_error_complex (out, expected, size);

    printf ("Max Error: %e\n", error);

    if (error > threshold)
    {
        ret = TEST_FAILURE;
        goto cleanup;
    }

cleanup:
    atfft_dft_destroy (fft);
    free (out);

    return ret;
}
