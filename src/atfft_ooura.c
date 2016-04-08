/*
 * Copyright (C) 2016 Sean Enderby <sean.enderby@gmail.com>
 *
 * This work is free. You can redistribute it and/or modify it under the
 * terms of the Do What The Fuck You Want To Public License, Version 2,
 * as published by Sam Hocevar. See the COPYING file for more details.
 */

#include <stdlib.h>
#include <string.h>
#include "atfft.h"
#include "fft4g.c"

struct atfft
{
    int size;
    int direction;
    double *data;
    int *workArea;
    double *tables;
};

struct atfft* atfft_create (int size, int direction)
{
    struct atfft *fft;

    fft = malloc (sizeof (*fft));
    fft->size = size;
    fft->direction = direction;
    fft->data = malloc (2 * size * sizeof (*(fft->data)));
    fft->workArea = malloc ((2 + (1 << (int) (log (size + 0.5) / log (2)) / 2)) * sizeof (*(fft->workArea)));
    fft->tables = malloc (2 * size * sizeof (*(fft->tables)));

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
    size_t nBytes = 2 * fft->size * sizeof (*in);
    memcpy (fft->data, in, nBytes);
    cdft (2 * fft->size, fft->direction, fft->data, fft->workArea, fft->tables);
    memcpy (out, fft->data, nBytes);
}
