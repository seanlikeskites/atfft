/*
 * Copyright (C) 2016 Sean Enderby <sean.enderby@gmail.com>
 *
 * This program is free software. It comes without any warranty, to
 * the extent permitted by applicable law. You can redistribute it 
 * and/or modify it under the terms of the Do What The Fuck You Want 
 * To Public License, Version 2, as published by Sam Hocevar. See 
 * the COPYING file for more details.
 */

#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <libavutil/mem.h>
#include <libavcodec/avfft.h>
#include <atfft/atfft_dft.h>

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
    FFTSample *data;
    void *context;
};

struct atfft_dft* atfft_dft_create (int size, enum atfft_direction direction, enum atfft_format format)
{
    /* ffmpeg only supports sizes which are a power of 2. */
    assert (atfft_is_power_of_2 (size));

    struct atfft_dft *fft = NULL;

    if (!(fft = malloc (sizeof (*fft))))
        return NULL;

    fft->size = size;
    fft->direction = direction;
    fft->format = format;

    init_context init_ctx = NULL;
    size_t data_size = 0;
    int init_direction = 0;

    switch (format)
    {
        case ATFFT_COMPLEX:
            data_size = 2 * size * sizeof (*(fft->data));
            init_ctx = (init_context) av_fft_init;

            if (direction == ATFFT_FORWARD)
                init_direction = 0;
            else
                init_direction = 1;
            break;

        case ATFFT_REAL:
            data_size = size * sizeof (*(fft->data));
            init_ctx = (init_context) av_rdft_init;

            if (direction == ATFFT_FORWARD)
                init_direction = DFT_R2C;
            else
                init_direction = IDFT_C2R;
            break;
    }

    /* allocate some fft stuff */
    fft->data = av_malloc (data_size);
    fft->context = init_ctx (log2 (size), init_direction);

    /* clean up on failure */
    if (!(fft->data && fft->context))
    {
        atfft_dft_destroy (fft);
        fft = NULL;
    }

    return fft;
}

void atfft_dft_destroy (struct atfft_dft *fft)
{
    if (fft)
    {
        switch (fft->format)
        {
            case ATFFT_COMPLEX:
                av_fft_end (fft->context);
                break;

            case ATFFT_REAL:
                av_rdft_end (fft->context);
                break;
        }

        av_free (fft->data);
        free (fft);
    }
}

void atfft_dft_complex_transform (struct atfft_dft *fft, atfft_complex *in, atfft_complex *out)
{
#ifdef ATFFT_TYPE_FLOAT
    size_t n_bytes = fft->size * sizeof (*in);
#endif 

    /* Only to be used with complex FFTs. */
    assert (fft->format == ATFFT_COMPLEX);

#ifdef ATFFT_TYPE_FLOAT
    memcpy (fft->data, in, n_bytes);
#else
    atfft_sample_to_float_complex (in, (atfft_complex_f*) fft->data, fft->size);
#endif

    av_fft_permute (fft->context, (FFTComplex*) fft->data);
    av_fft_calc (fft->context, (FFTComplex*) fft->data);

#ifdef ATFFT_TYPE_FLOAT
    memcpy (out, fft->data, n_bytes);
#else
    atfft_float_to_sample_complex ((atfft_complex_f*) fft->data, out, fft->size);
#endif
}

static void atfft_halfcomplex_ffmpeg_to_fftw (const FFTSample *in, atfft_complex *out, int size)
{
    int half_size = size / 2;

    ATFFT_RE (out [0]) = in [0];
    ATFFT_IM (out [0]) = 0;

    for (int i = 1; i < half_size; ++i)
    {
        ATFFT_RE (out [i]) = in [2 * i];
        ATFFT_IM (out [i]) = in [2 * i + 1];
    }

    ATFFT_RE (out [half_size]) = in [1];
    ATFFT_IM (out [half_size]) = 0;
}

void atfft_dft_real_forward_transform (struct atfft_dft *fft, const atfft_sample *in, atfft_complex *out)
{
    /* Only to be used for forward real FFTs. */
    assert ((fft->format == ATFFT_REAL) && (fft->direction == ATFFT_FORWARD));

#ifdef ATFFT_TYPE_FLOAT
    memcpy (fft->data, in, fft->size * sizeof (*(fft->data)));
#else
    atfft_sample_to_float_real (in, fft->data, fft->size);
#endif

    av_rdft_calc (fft->context, fft->data);
    atfft_halfcomplex_ffmpeg_to_fftw (fft->data, out, fft->size);
}

static void atfft_halfcomplex_fftw_to_ffmpeg (atfft_complex *in, FFTSample *out, int size)
{
    int half_size = size / 2;

    out [0] = 2.0 * ATFFT_RE (in [0]);

    for (int i = 1; i < half_size; ++i)
    {
        out [2 * i] = 2.0 * ATFFT_RE (in [i]);
        out [2 * i + 1] = 2.0 * ATFFT_IM (in [i]);
    }

    out [1] = 2.0 * ATFFT_RE (in [half_size]);
}

void atfft_dft_real_backward_transform (struct atfft_dft *fft, atfft_complex *in, atfft_sample *out)
{
    /* Only to be used for backward real FFTs. */
    assert ((fft->format == ATFFT_REAL) && (fft->direction == ATFFT_BACKWARD));

    atfft_halfcomplex_fftw_to_ffmpeg (in, fft->data, fft->size);
    av_rdft_calc (fft->context, fft->data);

#ifdef ATFFT_TYPE_FLOAT
    memcpy (out, fft->data, fft->size * sizeof (*out));
#else
    atfft_float_to_sample_real (fft->data, out, fft->size);
#endif
}
