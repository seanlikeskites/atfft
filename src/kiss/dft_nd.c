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
#include <atfft/dft_nd.h>

/* we need to make sure we are using kiss_fft with the correct type */
#define kiss_fft_scalar atfft_sample
#include <kiss_fftnd.h>
#include <kiss_fftndr.h>

struct atfft_dft_nd
{
    enum atfft_direction direction;
    enum atfft_format format;
    void *cfg;
};

struct atfft_dft_nd* atfft_dft_nd_create (const int *dims,
                                          int n_dims,
                                          enum atfft_direction direction,
                                          enum atfft_format format)
{
    /* kiss_fft only does even length real transforms */
    assert ((format == ATFFT_COMPLEX) || atfft_is_even (dims [n_dims - 1]));

    struct atfft_dft_nd *fft;

    if (!(fft = malloc (sizeof (*fft))))
        return NULL;

    fft->direction = direction;
    fft->format = format;

    switch (format)
    {
        case ATFFT_COMPLEX:
            if (direction == ATFFT_FORWARD)
                fft->cfg = kiss_fftnd_alloc (dims, n_dims, 0, 0, 0);
            else
                fft->cfg = kiss_fftnd_alloc (dims, n_dims, 1, 0, 0);
            break;

        case ATFFT_REAL:
            if (direction == ATFFT_FORWARD)
                fft->cfg = kiss_fftndr_alloc (dims, n_dims, 0, 0, 0);
            else
                fft->cfg = kiss_fftndr_alloc (dims, n_dims, 1, 0, 0);
            break;
    };

    /* clean up on failure */
    if (!fft->cfg)
    {
        atfft_dft_nd_destroy (fft);
        fft = NULL;
    }

    return fft;
}

void atfft_dft_nd_destroy (struct atfft_dft_nd *fft)
{
    if (fft)
    {
        free (fft->cfg);
        free (fft);
    }
}

void atfft_dft_nd_complex_transform (struct atfft_dft_nd *fft, atfft_complex *in, atfft_complex *out)
{
    /* Only to be used with complex FFTs. */
    assert (fft->format == ATFFT_COMPLEX);

    kiss_fftnd ((kiss_fftnd_cfg) fft->cfg, (kiss_fft_cpx*) in, (kiss_fft_cpx*) out);
}

void atfft_dft_nd_real_forward_transform (struct atfft_dft_nd *fft, const atfft_sample *in, atfft_complex *out)
{
    /* Only to be used for forward real FFTs. */
    assert ((fft->format == ATFFT_REAL) && (fft->direction == ATFFT_FORWARD));

    kiss_fftndr (fft->cfg, in, (kiss_fft_cpx*) out);
}

void atfft_dft_nd_real_backward_transform (struct atfft_dft_nd *fft, atfft_complex *in, atfft_sample *out)
{
    /* Only to be used for backward real FFTs. */
    assert ((fft->format == ATFFT_REAL) && (fft->direction == ATFFT_BACKWARD));

    kiss_fftndri (fft->cfg, (kiss_fft_cpx*) in, out);
}
