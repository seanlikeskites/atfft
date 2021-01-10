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
#include <string.h>
#include <assert.h>
#include <math.h>
#include <libavutil/mem.h>
#include <libavcodec/avfft.h>
#include <atfft/dft.h>

#ifndef ATFFT_TYPE_FLOAT
#   ifdef _MSC_VER
#       pragma message(": warning: FFmpeg only supports single precision floating point, " \
                       "higher precision values will be demoted to float for FFT calculations.")
#   else
#       warning FFmpeg only supports single precision floating point, \
                higher precision values will be demoted to float for FFT calculations.
#   endif
#endif

typedef void* (*init_context) (int, int);

struct atfft_dft
{
    int size;
    enum atfft_direction direction;
    enum atfft_format format;

    /* the libavcodec FFT context */
    void *context;

    /* buffer for in place transform */
    size_t n_data_bytes;
    FFTSample *data;
};

int atfft_dft_is_supported_size (int size, enum atfft_format format)
{
    int min = 1 << 2;
    int max = 1 << 16;

    if (format == ATFFT_REAL)
        min = 1 << 4;

    return atfft_is_power_of_2 (size) && size >= min && size <= max;
}

struct atfft_dft* atfft_dft_create (int size, enum atfft_direction direction, enum atfft_format format)
{
    /* ffmpeg only supports sizes which are a power of 2. */
    assert (atfft_dft_is_supported_size (size, format));

    struct atfft_dft *plan = NULL;

    if (!(plan = calloc (1, sizeof (*plan))))
        return NULL;

    plan->size = size;
    plan->direction = direction;
    plan->format = format;

    init_context init_ctx = NULL;
    int init_direction = 0;

    if (format == ATFFT_COMPLEX)
    {
        plan->n_data_bytes = 2 * size * sizeof (*(plan->data));
        init_ctx = (init_context) av_fft_init;

        if (direction == ATFFT_FORWARD)
            init_direction = 0;
        else
            init_direction = 1;
    }
    else
    {
        plan->n_data_bytes = size * sizeof (*(plan->data));
        init_ctx = (init_context) av_rdft_init;

        if (direction == ATFFT_FORWARD)
            init_direction = DFT_R2C;
        else
            init_direction = IDFT_C2R;
    }

    /* allocate libavcodec FFT context and data buffer */
    plan->context = init_ctx (log2 (size), init_direction);
    plan->data = av_malloc (plan->n_data_bytes);

    /* clean up on failure */
    if (!(plan->data && plan->context))
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
        av_free (plan->data);

        if (plan->format == ATFFT_COMPLEX)
            av_fft_end (plan->context);
        else
            av_rdft_end (plan->context);

        free (plan);
    }
}

void atfft_dft_complex_transform (struct atfft_dft *plan, atfft_complex *in, atfft_complex *out)
{
    /* Only to be used with complex FFTs. */
    assert (plan->format == ATFFT_COMPLEX);

#ifdef ATFFT_TYPE_FLOAT
    memcpy (plan->data, in, plan->n_data_bytes);
#else
    atfft_sample_to_float_complex (in, (atfft_complex_f*) plan->data, plan->size);
#endif

    av_fft_permute (plan->context, (FFTComplex*) plan->data);
    av_fft_calc (plan->context, (FFTComplex*) plan->data);

#ifdef ATFFT_TYPE_FLOAT
    memcpy (out, plan->data, plan->n_data_bytes);
#else
    atfft_float_to_sample_complex ((atfft_complex_f*) plan->data, out, plan->size);
#endif
}

void atfft_dft_complex_transform_stride (struct atfft_dft *plan,
                                         atfft_complex *in,
                                         int in_stride,
                                         atfft_complex *out,
                                         int out_stride)
{
    /* Only to be used with complex FFTs. */
    assert (plan->format == ATFFT_COMPLEX);

    atfft_sample_to_float_complex_stride (in,
                                          in_stride,
                                          (atfft_complex_f*) plan->data,
                                          1,
                                          plan->size);

    av_fft_permute (plan->context, (FFTComplex*) plan->data);
    av_fft_calc (plan->context, (FFTComplex*) plan->data);

    atfft_float_to_sample_complex_stride ((atfft_complex_f*) plan->data,
                                          1,
                                          out,
                                          out_stride,
                                          plan->size);
}

static void halfcomplex_ffmpeg_to_atfft (const FFTSample *in,
                                         atfft_complex *out,
                                         int out_stride,
                                         int size)
{
    int half_size = size / 2;

    ATFFT_RE (out [0]) = in [0];
    ATFFT_IM (out [0]) = 0;

    int o = out_stride;

    for (int i = 1; i < half_size; ++i, o += out_stride)
    {
        ATFFT_RE (out [o]) = in [2 * i];
        ATFFT_IM (out [o]) = in [2 * i + 1];
    }

    ATFFT_RE (out [o]) = in [1];
    ATFFT_IM (out [o]) = 0;
}

void atfft_dft_real_forward_transform (struct atfft_dft *plan, const atfft_sample *in, atfft_complex *out)
{
    /* Only to be used for forward real FFTs. */
    assert ((plan->format == ATFFT_REAL) && (plan->direction == ATFFT_FORWARD));

#ifdef ATFFT_TYPE_FLOAT
    memcpy (plan->data, in, plan->n_data_bytes);
#else
    atfft_sample_to_float_real (in, plan->data, plan->size);
#endif

    av_rdft_calc (plan->context, plan->data);
    halfcomplex_ffmpeg_to_atfft (plan->data, out, 1, plan->size);
}

void atfft_dft_real_forward_transform_stride (struct atfft_dft *plan,
                                              const atfft_sample *in,
                                              int in_stride,
                                              atfft_complex *out,
                                              int out_stride)
{
    /* Only to be used for forward real FFTs. */
    assert ((plan->format == ATFFT_REAL) && (plan->direction == ATFFT_FORWARD));

    atfft_sample_to_float_real_stride (in,
                                       in_stride,
                                       plan->data,
                                       1,
                                       plan->size);

    av_rdft_calc (plan->context, plan->data);
    halfcomplex_ffmpeg_to_atfft (plan->data, out, out_stride, plan->size);
}

static void halfcomplex_atfft_to_ffmpeg (atfft_complex *in,
                                         int in_stride,
                                         FFTSample *out,
                                         int size)
{
    int half_size = size / 2;

    out [0] = 2.0 * ATFFT_RE (in [0]);

    int i = in_stride;

    for (int o = 1; o < half_size; i += in_stride, ++o)
    {
        out [2 * o] = 2.0 * ATFFT_RE (in [i]);
        out [2 * o + 1] = 2.0 * ATFFT_IM (in [i]);
    }

    out [1] = 2.0 * ATFFT_RE (in [i]);
}

void atfft_dft_real_backward_transform (struct atfft_dft *plan, atfft_complex *in, atfft_sample *out)
{
    /* Only to be used for backward real FFTs. */
    assert ((plan->format == ATFFT_REAL) && (plan->direction == ATFFT_BACKWARD));

    halfcomplex_atfft_to_ffmpeg (in, 1, plan->data, plan->size);
    av_rdft_calc (plan->context, plan->data);

#ifdef ATFFT_TYPE_FLOAT
    memcpy (out, plan->data, plan->n_data_bytes);
#else
    atfft_float_to_sample_real (plan->data, out, plan->size);
#endif
}

void atfft_dft_real_backward_transform_stride (struct atfft_dft *plan,
                                               atfft_complex *in,
                                               int in_stride,
                                               atfft_sample *out,
                                               int out_stride)
{
    /* Only to be used for backward real FFTs. */
    assert ((plan->format == ATFFT_REAL) && (plan->direction == ATFFT_BACKWARD));

    halfcomplex_atfft_to_ffmpeg (in, in_stride, plan->data, plan->size);
    av_rdft_calc (plan->context, plan->data);
    atfft_float_to_sample_real_stride (plan->data,
                                       1,
                                       out,
                                       out_stride,
                                       plan->size);
}
