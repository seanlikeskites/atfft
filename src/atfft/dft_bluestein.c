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
#include <atfft/dft.h>
#include <atfft/dft_util.h>
#include "atfft_internal.h"
#include "dft_bluestein.h"

/**
 * Swap the real and imaginary parts of a and then multiply by b:
 *
 * p = j * conj(a) * b
 */
static void atfft_swap_and_product_complex (const atfft_complex a,
                                            const atfft_complex b,
                                            atfft_complex *p)
{
    ATFFT_RE (*p) = ATFFT_IM (a) * ATFFT_RE (b) -
                    ATFFT_RE (a) * ATFFT_IM (b);
    ATFFT_IM (*p) = ATFFT_IM (a) * ATFFT_IM (b) +
                    ATFFT_RE (a) * ATFFT_RE (b);
}


struct atfft_dft_bluestein
{
    int size;
    enum atfft_direction direction;
    enum atfft_format format;
    int conv_size;
    struct atfft_dft *fft;
    atfft_complex *sig, *sig_dft, *conv, *conv_dft, *factors;
};

static int atfft_bluestein_convolution_fft_size (int size)
{
    if (atfft_is_power_of_2 (size))
        return size;
    else
        return atfft_next_power_of_2 (2 * size - 1);
}

static int atfft_init_bluestein_convolution_dft (int size,
                                                 enum atfft_direction direction,
                                                 atfft_complex *conv_dft,
                                                 int conv_size,
                                                 atfft_complex *factors,
                                                 struct atfft_dft *fft)
{
    atfft_complex *sin_table = malloc (2 * size * sizeof (*sin_table));
    atfft_complex *sequence = calloc (conv_size, sizeof (*sequence));
    int ret = 0;

    if (!(sequence && sin_table))
    {
        ret = -1;
        goto failed;
    }

    /* calculate sin table */
    for (int i = 0; i < 2 * size; ++i)
    {
        atfft_twiddle_factor (-i, 2 * size, direction, sin_table + i);
    }

    /* produce convolution sequence */
    for (int i = 0; i < size; ++i)
    {
        int index = (i * i) % (2 * size);
        atfft_copy_complex (sin_table [index], &sequence [i]);
    }

    /* replicate samples for circular convolution */
    if (conv_size > size)
    {
        for (int i = 1; i < size; ++i)
        {
            atfft_copy_complex (sequence [i], &sequence [conv_size - i]);
        }
    }

    /* take DFT of sequence */
    atfft_dft_complex_transform (fft, sequence, conv_dft);
    atfft_normalise_dft_complex (conv_dft, conv_size);

    /* take conjugate of the sequence for use later */
    for (int i = 0; i < size; ++i)
    {
        ATFFT_RE (factors [i]) = ATFFT_RE (sequence [i]);
        ATFFT_IM (factors [i]) = - ATFFT_IM (sequence [i]);
    }

failed:
    free (sequence);
    free (sin_table);
    return ret;
}

struct atfft_dft_bluestein* atfft_dft_bluestein_create (int size,
                                                        enum atfft_direction direction,
                                                        enum atfft_format format)
{
    struct atfft_dft_bluestein *fft;

    if (!(fft = calloc (1, sizeof (*fft))))
        return NULL;

    fft->size = size;
    fft->direction = direction;
    fft->format = format;

    /* allocate some regular dft objects for performing the convolution */
    fft->conv_size = atfft_bluestein_convolution_fft_size (size);
    fft->fft = atfft_dft_create (fft->conv_size, ATFFT_FORWARD, ATFFT_COMPLEX);

    if (!fft->fft)
        goto failed;

    /* allocate work space for performing the dft */
    fft->sig = calloc (fft->conv_size, sizeof (*(fft->sig))); /* set up for zero padding */
    fft->sig_dft = malloc (fft->conv_size * sizeof (*(fft->sig_dft)));
    fft->conv = malloc (fft->conv_size * sizeof (*(fft->conv)));
    fft->conv_dft = malloc (fft->conv_size * sizeof (*(fft->conv_dft)));
    fft->factors = malloc (size * sizeof (*(fft->factors)));

    if (!(fft->sig && fft->sig_dft && fft->conv && fft->conv_dft && fft->factors))
        goto failed;

    /* calculate the convolution dft */
    if (atfft_init_bluestein_convolution_dft (size,
                                              direction,
                                              fft->conv_dft,
                                              fft->conv_size,
                                              fft->factors,
                                              fft->fft) < 0)
        goto failed;

    return fft;

failed:
    atfft_dft_bluestein_destroy (fft);
    return NULL;
}

void atfft_dft_bluestein_destroy (void *fft)
{
    struct atfft_dft_bluestein *t = fft;

    if (t)
    {
        free (t->factors);
        free (t->conv_dft);
        free (t->conv);
        free (t->sig_dft);
        free (t->sig);
        atfft_dft_destroy (t->fft);
        free (t);
    }
}

void atfft_dft_bluestein_complex_transform (void *fft,
                                            atfft_complex *in,
                                            int in_stride,
                                            atfft_complex *out,
                                            int out_stride)
{
    struct atfft_dft_bluestein *t = fft;

    /* multiply the input signal with the factors */
    for (int i = 0; i < t->size; ++i)
    {
        atfft_product_complex (in [i * in_stride], t->factors [i], t->sig + i);
    }

    /* take DFT of the result */
    atfft_dft_complex_transform (t->fft, t->sig, t->sig_dft);

    /* perform convolution in the frequency domain */
    for (int i = 0; i < t->conv_size; ++i)
    {
        atfft_multiply_by_and_swap_complex (t->sig_dft + i, t->conv_dft [i]);
    }

    /* take the inverse DFT of the result */
    atfft_dft_complex_transform (t->fft, t->sig_dft, t->conv);

    /* multiply the output transform with the factors */
    for (int i = 0; i < t->size; ++i)
    {
        atfft_swap_and_product_complex (t->conv [i], t->factors [i], out + i * out_stride);
    }
}
