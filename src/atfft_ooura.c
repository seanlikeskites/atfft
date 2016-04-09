/*
 * Copyright (C) 2016 Sean Enderby <sean.enderby@gmail.com>
 *
 * This work is free. You can redistribute it and/or modify it under the
 * terms of the Do What The Fuck You Want To Public License, Version 2,
 * as published by Sam Hocevar. See the COPYING file for more details.
 */

#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "../include/atfft.h"
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
    assert (isPowerOf2 (size));

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

void atfft_complex_transform (struct atfft *fft, double *in, double *out)
{
    /* Only to be used with complex FFTs. */
    assert (fft->format == ATFFT_COMPLEX);

    memcpy (fft->data, in, fft->dataSize);
    cdft (2 * fft->size, fft->direction, fft->data, fft->workArea, fft->tables);
    memcpy (out, fft->data, fft->dataSize);
}

atfft_half_complex_ooura_to_fftw (double *in, double *out, int size)
{
    int i = 0;

    out [0] = in [0];
    out [1] = 0;

    for (i = 2; i < size; i += 2)
    {
        out [i] = in [i];
        out [i + 1] = -in [i + 1];
    }

    out [size] = in [1];
    out [size + 1] = 0;
}

void atfft_real_forward_transform (struct atfft *fft, double *in, double *out)
{
    /* Only to be used for forward real FFTs. */
    assert ((fft->format == ATFFT_REAL) && (fft->direction == ATFFT_FORWARD));

    memcpy (fft->data, in, fft->dataSize);
    rdft (fft->size, 1, fft->data, fft->workArea, fft->tables);
    atfft_half_complex_ooura_to_fftw (fft->data, out, fft->size);
}

atfft_half_complex_fftw_to_ooura (double *in, double *out, int size)
{
    int i = 0;

    out [0] = in [0];

    for (i = 2; i < size; i += 2)
    {
        out [i] = in [i];
        out [i + 1] = -in [i + 1] * 2;
    }

    out [1] = in [size];
}

void atfft_real_backward_transform (struct atfft *fft, double *in, double *out)
{
    /* Only to be used for backward real FFTs. */
    assert ((fft->format == ATFFT_REAL) && (fft->direction == ATFFT_BACKWARD));

    atfft_half_complex_fftw_to_ooura (in, fft->data, fft->size);
    rdft (fft->size, -1, fft->data, fft->workArea, fft->tables);
    memcpy (out, fft->data, fft->dataSize);
}
