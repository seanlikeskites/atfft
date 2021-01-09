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

    /* the libavcodec DCT context */
    DCTContext *context;

    /* buffer for libavcodec to operate on */
    size_t n_data_bytes;
    FFTSample *data;
};

int atfft_dct_is_supported_size (int size)
{
    int min = 1 << 4;
    int max = 1 << 16;

    return atfft_is_power_of_2 (size) && size >= min && size <= max;
}

struct atfft_dct* atfft_dct_create (int size, enum atfft_direction direction)
{
    /* ffmpeg only supports sizes which are a power of 2. */
    assert (atfft_dct_is_supported_size (size));

    struct atfft_dct *plan = NULL;

    if (!(plan = calloc (1, sizeof (*plan))))
        return NULL;

    plan->size = size;
    plan->direction = direction;

    if (direction == ATFFT_FORWARD)
        plan->context = av_dct_init (log2 (size), DCT_II);
    else
        plan->context = av_dct_init (log2 (size), DCT_III);

    plan->n_data_bytes = size * sizeof (*(plan->data));
    plan->data = av_malloc (plan->n_data_bytes);

    /* clean up on failure */
    if (!(plan->data && plan->context))
    {
        atfft_dct_destroy (plan);
        plan = NULL;
    }

    return plan;
}

void atfft_dct_destroy (struct atfft_dct *plan)
{
    if (plan)
    {
        av_free (plan->data);
        av_dct_end (plan->context);
        free (plan);
    }
}

void atfft_dct_transform (struct atfft_dct *plan, const atfft_sample *in, atfft_sample *out)
{
#ifdef ATFFT_TYPE_FLOAT
    memcpy (plan->data, in, plan->n_data_bytes);
#else
    atfft_sample_to_float_real (in, plan->data, plan->size);
#endif

    av_dct_calc (plan->context, plan->data);

#ifdef ATFFT_TYPE_FLOAT
    memcpy (out, plan->data, n_bytes);
#else
    atfft_float_to_sample_real (plan->data, out, plan->size);
#endif

    /* libavcodec normalises the backward transform */
    if (plan->direction == ATFFT_BACKWARD)
        atfft_scale_real (out, plan->size, plan->size);
}
