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
#include <string.h>
#include <assert.h>
#include <atfft.h>
#include "ooura/fft4g.c"

struct atfft
{
    int size;
    int direction;
    enum atfft_format format;
    int dataSize;
    double *data;
    int *workArea;
    double *tables;
};

struct atfft* atfft_create (int size, int direction, enum atfft_format format)
{
    struct atfft *fft;

    /* Ooura only supports sizes which are a power of 2. */
    assert (atfft_is_power_of_2 (size));

    fft = malloc (sizeof (*fft));
    fft->size = size;
    fft->direction = direction;
    fft->format = format;

    switch (format)
    {
        case ATFFT_COMPLEX:
            fft->dataSize = 2 * size * sizeof (*(fft->data));
            fft->data = malloc (fft->dataSize);
            fft->workArea = malloc ((2 + (1 << (int) (log (size + 0.5) / log (2)) / 2)) * sizeof (*(fft->workArea)));
            break;

        case ATFFT_REAL:
            fft->dataSize = size * sizeof (*(fft->data));
            fft->data = malloc (fft->dataSize);
            fft->workArea = malloc ((2 + (1 << (int) (log (size / 2 + 0.5) / log (2)) / 2)) * sizeof (*(fft->workArea)));
    }

    fft->tables = malloc (size / 2 * sizeof (*(fft->tables)));
    fft->workArea [0] = 0;

    return fft;
}

void atfft_free (struct atfft *fft)
{
    free (fft->tables);
    free (fft->workArea);
    free (fft->data);
    free (fft);
}

void atfft_complex_transform (struct atfft *fft, atfft_complex *in, atfft_complex *out)
{
    /* Only to be used with complex FFTs. */
    assert (fft->format == ATFFT_COMPLEX);

#ifdef ATFFT_TYPE_DOUBLE
    memcpy (fft->data, in, fft->dataSize);
#else
    atfft_sample_to_double_complex (in, (atfft_complex_d*) fft->data, fft->size);
#endif

    cdft (2 * fft->size, fft->direction, fft->data, fft->workArea, fft->tables);

#ifdef ATFFT_TYPE_DOUBLE
    memcpy (out, fft->data, fft->dataSize);
#else
    atfft_double_to_sample_complex ((atfft_complex_d*) fft->data, out, fft->size);
#endif
}

void atfft_halfcomplex_ooura_to_fftw (double *in, atfft_complex *out, int size)
{
    int i = 0;
    int halfSize = size / 2;

    ATFFT_REAL (out [0]) = in [0];
    ATFFT_IMAG (out [0]) = 0;

    for (i = 1; i < halfSize; ++i)
    {
        ATFFT_REAL (out [i]) = in [2 * i];
        ATFFT_IMAG (out [i]) = -in [2 * i + 1];
    }

    ATFFT_REAL (out [halfSize]) = in [1];
    ATFFT_IMAG (out [halfSize]) = 0;
}

void atfft_real_forward_transform (struct atfft *fft, atfft_sample *in, atfft_complex *out)
{
    /* Only to be used for forward real FFTs. */
    assert ((fft->format == ATFFT_REAL) && (fft->direction == ATFFT_FORWARD));

#ifdef ATFFT_TYPE_DOUBLE
    memcpy (fft->data, in, fft->dataSize);
#else
    atfft_sample_to_double_real (in, fft->data, fft->size);
#endif

    rdft (fft->size, 1, fft->data, fft->workArea, fft->tables);
    atfft_halfcomplex_ooura_to_fftw (fft->data, out, fft->size);
}

void atfft_halfcomplex_fftw_to_ooura (atfft_complex *in, double *out, int size)
{
    int i = 0;
    int halfSize = size / 2;

    out [0] = 2.0 * ATFFT_REAL (in [0]);

    for (i = 1; i < halfSize; ++i)
    {
        out [2 * i] = 2.0 * ATFFT_REAL (in [i]);
        out [2 * i + 1] = 2.0 * - ATFFT_IMAG (in [i]);
    }

    out [1] = 2.0 * ATFFT_REAL (in [halfSize]);
}

void atfft_real_backward_transform (struct atfft *fft, atfft_complex *in, atfft_sample *out)
{
    /* Only to be used for backward real FFTs. */
    assert ((fft->format == ATFFT_REAL) && (fft->direction == ATFFT_BACKWARD));

    atfft_halfcomplex_fftw_to_ooura (in, fft->data, fft->size);
    rdft (fft->size, -1, fft->data, fft->workArea, fft->tables);

#ifdef ATFFT_TYPE_DOUBLE
    memcpy (out, fft->data, fft->dataSize);
#else
    atfft_double_to_sample_real (fft->data, out, fft->size);
#endif
}