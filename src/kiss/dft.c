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

#include <stdlib.h>
#include <assert.h>
#include <atfft/dft.h>

/* we need to make sure we are using kiss_fft with the correct type */
#define kiss_fft_scalar atfft_sample
#include <kiss_fft.h>
#include <kiss_fftr.h>

struct atfft_dft
{
    int size;
    enum atfft_direction direction;
    enum atfft_format format;

    /* the kissfft config */
    void *cfg;
};

int atfft_dft_is_supported_size (int size, enum atfft_format format)
{
    if (format == ATFFT_COMPLEX)
        return size > 0;
    else
        return size > 0 && atfft_is_even(size);
}

struct atfft_dft* atfft_dft_create (int size, enum atfft_direction direction, enum atfft_format format)
{
    /* kissfft supports all sizes for complex transforms, but real transforms must be even. */
    assert (atfft_dft_is_supported_size (size, format));

    struct atfft_dft *plan;

    if (!(plan = calloc (1, sizeof (*plan))))
        return NULL;

    plan->size = size;
    plan->direction = direction;
    plan->format = format;

    int inverse = direction == ATFFT_BACKWARD;

    if (format == ATFFT_COMPLEX)
        plan->cfg = kiss_fft_alloc (size, inverse, 0, 0);
    else
        plan->cfg = kiss_fftr_alloc (size, inverse, 0, 0);

    /* clean up on failure */
    if (!plan->cfg)
    {
        atfft_dft_destroy (plan);
        plan = NULL;
    }

    return plan;
}

void atfft_dft_destroy (struct atfft_dft *plan)
{
    if (plan)
    {
        free (plan->cfg);
        free (plan);
    }
}

void atfft_dft_complex_transform (struct atfft_dft *plan, atfft_complex *in, atfft_complex *out)
{
    /* Only to be used with complex FFTs. */
    assert (plan->format == ATFFT_COMPLEX);

    kiss_fft ((kiss_fft_cfg) plan->cfg, (kiss_fft_cpx*) in, (kiss_fft_cpx*) out);
}

void atfft_dft_real_forward_transform (struct atfft_dft *plan, const atfft_sample *in, atfft_complex *out)
{
    /* Only to be used for forward real FFTs. */
    assert ((plan->format == ATFFT_REAL) && (plan->direction == ATFFT_FORWARD));

    kiss_fftr (plan->cfg, in, (kiss_fft_cpx*) out);
}

void atfft_dft_real_backward_transform (struct atfft_dft *plan, atfft_complex *in, atfft_sample *out)
{
    /* Only to be used for backward real FFTs. */
    assert ((plan->format == ATFFT_REAL) && (plan->direction == ATFFT_BACKWARD));

    kiss_fftri (plan->cfg, (kiss_fft_cpx*) in, out);
}
