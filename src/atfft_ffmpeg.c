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
#include <atfft.h>

typedef void* (*init_context) (int, int);

struct atfft
{
    int size;
    enum atfft_direction direction;
    enum atfft_format format;
    FFTSample *data;
    void *context;
};

struct atfft* atfft_create (int size, enum atfft_direction direction, enum atfft_format format)
{
    struct atfft *fft;

    init_context initContext;
    size_t dataSize;
    int initDirection;

    /* ffmpeg only supports sizes which are a power of 2. */
    assert (atfft_is_power_of_2 (size));

    if (!(fft = malloc (sizeof (*fft))))
        return NULL;

    fft->size = size;
    fft->direction = direction;
    fft->format = format;

    switch (format)
    {
        case ATFFT_COMPLEX:
            dataSize = 2 * size * sizeof (*(fft->data));
            initContext = (init_context) av_fft_init;

            if (direction == ATFFT_FORWARD)
                initDirection = 0;
            else
                initDirection = 1;
            break;

        case ATFFT_REAL:
            dataSize = size * sizeof (*(fft->data));
            initContext = (init_context) av_rdft_init;

            if (direction == ATFFT_FORWARD)
                initDirection = DFT_R2C;
            else
                initDirection = IDFT_C2R;
            break;
    }

    /* allocate some fft stuff */
    fft->data = av_malloc (dataSize);
    fft->context = initContext (log2 (size), initDirection);

    /* clean up on failure */
    if (!(fft->data && fft->context))
    {
        atfft_destroy (fft);
        fft = NULL;
    }

    return fft;
}

void atfft_destroy (struct atfft *fft)
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

void atfft_complex_transform (struct atfft *fft, atfft_complex *in, atfft_complex *out)
{
#ifdef ATFFT_TYPE_FLOAT
    size_t nBytes = fft->size * sizeof (*in);
#endif 

    /* Only to be used with complex FFTs. */
    assert (fft->format == ATFFT_COMPLEX);

#ifdef ATFFT_TYPE_FLOAT
    memcpy (fft->data, in, nBytes);
#else
    atfft_sample_to_float_complex (in, (atfft_complex_f*) fft->data, fft->size);
#endif

    av_fft_permute (fft->context, (FFTComplex*) fft->data);
    av_fft_calc (fft->context, (FFTComplex*) fft->data);

#ifdef ATFFT_TYPE_FLOAT
    memcpy (out, fft->data, nBytes);
#else
    atfft_float_to_sample_complex ((atfft_complex_f*) fft->data, out, fft->size);
#endif
}

void atfft_halfcomplex_ffmpeg_to_fftw (FFTSample *in, atfft_complex *out, int size)
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

void atfft_real_forward_transform (struct atfft *fft, atfft_sample *in, atfft_complex *out)
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

void atfft_halfcomplex_fftw_to_ffmpeg (atfft_complex *in, FFTSample *out, int size)
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

void atfft_real_backward_transform (struct atfft *fft, atfft_complex *in, atfft_sample *out)
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
