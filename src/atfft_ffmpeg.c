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
#include <string.h>
#include <libavcodec/avcodec.h>
#include <libavcodec/avfft.h>
#include <atfft.h>

struct atfft
{
    int size;
    int direction;
    enum atfft_format format;
    FFTSample *data;
    void *context;
};

struct atfft* atfft_create (int size, int direction, enum atfft_format format)
{
    struct atfft *fft;

    /* ffmpeg only supports sizes which are a power of 2. */
    assert (atfft_is_power_of_2 (size));

    fft = malloc (sizeof (*fft));
    fft->size = size;
    fft->direction = direction;
    fft->format = format;

    switch (format)
    {
        case ATFFT_COMPLEX:
            fft->data = av_malloc (2 * size * sizeof (*(fft->data)));

            if (direction == ATFFT_FORWARD)
            {
                fft->context = av_fft_init (log2 (size), 0);
            }
            else
            {
                fft->context = av_fft_init (log2 (size), 1);
            }

            break;

        case ATFFT_REAL:
            fft->data = av_malloc (size * sizeof (*(fft->data)));

            if (direction == ATFFT_FORWARD)
            {
                fft->context = av_rdft_init (log2 (size), DFT_R2C);
            }
            else
            {
                fft->context = av_rdft_init (log2 (size), IDFT_C2R);
            }

            break;
    }

    return fft;
}

void atfft_free (struct atfft *fft)
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

void atfft_complex_transform (struct atfft *fft, atfft_complex_double *in, atfft_complex_double *out)
{
    /* Only to be used with complex FFTs. */
    assert (fft->format == ATFFT_COMPLEX);

    atfft_double_to_float_complex (in, (atfft_complex_float*) fft->data, fft->size);
    av_fft_permute (fft->context, (FFTComplex*) fft->data);
    av_fft_calc (fft->context, (FFTComplex*) fft->data);
    atfft_float_to_double_complex ((atfft_complex_float*) fft->data, out, fft->size);
}

void atfft_halfcomplex_ffmpeg_to_fftw (FFTSample *in, atfft_complex_double *out, int size)
{
    int i = 0;
    int halfSize = size / 2;

    ATFFT_REAL (out [0]) = in [0];
    ATFFT_IMAG (out [0]) = 0;

    for (i = 1; i < halfSize; ++i)
    {
        ATFFT_REAL (out [i]) = in [2 * i];
        ATFFT_IMAG (out [i]) = in [2 * i + 1];
    }

    ATFFT_REAL (out [halfSize]) = in [1];
    ATFFT_IMAG (out [halfSize]) = 0;
}

void atfft_real_forward_transform (struct atfft *fft, double *in, atfft_complex_double *out)
{
    /* Only to be used for forward real FFTs. */
    assert ((fft->format == ATFFT_REAL) && (fft->direction == ATFFT_FORWARD));

    atfft_double_to_float_real (in, fft->data, fft->size);
    av_rdft_calc (fft->context, fft->data);
    atfft_halfcomplex_ffmpeg_to_fftw (fft->data, out, fft->size);
}

void atfft_halfcomplex_fftw_to_ffmpeg (atfft_complex_double *in, FFTSample *out, int size)
{
    int i = 0;
    int halfSize = size / 2;

    out [0] = 2.0 * ATFFT_REAL (in [0]);

    for (i = 1; i < halfSize; ++i)
    {
        out [2 * i] = 2.0 * ATFFT_REAL (in [i]);
        out [2 * i + 1] = 2.0 * ATFFT_IMAG (in [i]);
    }

    out [1] = 2.0 * ATFFT_REAL (in [halfSize]);
}

void atfft_real_backward_transform (struct atfft *fft, atfft_complex_double *in, double *out)
{
    /* Only to be used for backward real FFTs. */
    assert ((fft->format == ATFFT_REAL) && (fft->direction == ATFFT_BACKWARD));

    atfft_halfcomplex_fftw_to_ffmpeg (in, fft->data, fft->size);
    av_rdft_calc (fft->context, fft->data);
    atfft_float_to_double_real (fft->data, out, fft->size);
}
