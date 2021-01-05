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
#include <math.h>
#include <libavutil/mem.h>
#include <libavcodec/avfft.h>
#include <atfft/dct.h>

#ifndef ATFFT_TYPE_FLOAT
#   ifdef _MSC_VER
#       pragma message(": warning: FFmpeg only supports single precision floating point, " \
                       "higher precision values will be demoted to float for FFT calculations.")
#   else
#       warning FFmpeg only supports single precision floating point, \
                higher precision values will be demoted to float for FFT calculations.
#   endif
#endif

struct atfft_dct
{
    int size;
    enum atfft_direction direction;
    FFTSample *data;
    DCTContext *context;
};

struct atfft_dct* atfft_dct_create (int size, enum atfft_direction direction)
{
    /* ffmpeg only supports sizes which are a power of 2. */
    assert (atfft_is_power_of_2 (size));

    struct atfft_dct *dct = NULL;

    if (!(dct = malloc (sizeof (*dct))))
        return NULL;

    dct->size = size;
    dct->direction = direction;

    size_t data_size = size * sizeof (*(dct->data));
    dct->data = av_malloc (data_size);

    if (direction == ATFFT_FORWARD)
        dct->context = av_dct_init (log2 (size), DCT_II);
    else
        dct->context = av_dct_init (log2 (size), DCT_III);

    /* clean up on failure */
    if (!(dct->data && dct->context))
    {
        atfft_dct_destroy (dct);
        dct = NULL;
    }

    return dct;
}

void atfft_dct_destroy (struct atfft_dct *dct)
{
    if (dct)
    {
        av_dct_end (dct->context);
        av_free (dct->data);
        free (dct);
    }
}

void atfft_dct_transform (struct atfft_dct *dct, const atfft_sample *in, atfft_sample *out)
{
#ifdef ATFFT_TYPE_FLOAT
    size_t n_bytes = dct->size * sizeof (*in);
#endif 

#ifdef ATFFT_TYPE_FLOAT
    memcpy (dct->data, in, n_bytes);
#else
    atfft_sample_to_float_real (in, dct->data, dct->size);
#endif

    av_dct_calc (dct->context, dct->data);

#ifdef ATFFT_TYPE_FLOAT
    memcpy (out, dct->data, n_bytes);
#else
    atfft_float_to_sample_real (dct->data, out, dct->size);
#endif

    if (dct->direction == ATFFT_BACKWARD)
        atfft_scale_real (out, dct->size, dct->size);
}
