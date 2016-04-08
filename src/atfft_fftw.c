/*
 * Copyright (C) 2016 Sean Enderby <sean.enderby@gmail.com>
 *
 * This work is free. You can redistribute it and/or modify it under the
 * terms of the Do What The Fuck You Want To Public License, Version 2,
 * as published by Sam Hocevar. See the COPYING file for more details.
 */

#include <fftw3.h>
#include <stdlib.h>
#include <string.h>
#include "../include/atfft.h"

struct atfft
{
    int size;
    int direction;
    fftw_complex *in, *out;
    fftw_plan plan;
};

struct atfft* atfft_create (int size, int direction)
{
    struct atfft *fft;

    fft = malloc (sizeof (*fft));
    fft->size = size;
    fft->direction = direction;
    fft->in = fftw_malloc (size * sizeof (*(fft->in)));
    fft->out = fftw_malloc (size * sizeof (*(fft->out)));
    fft->plan = fftw_plan_dft_1d (size, fft->in, fft->out, direction, FFTW_ESTIMATE);

    return fft;
}

void atfft_free (struct atfft *fft)
{
    fftw_destroy_plan (fft->plan);
    fftw_free (fft->out);
    fftw_free (fft->in);
    free (fft);
}

void atfft_complex_transform (struct atfft *fft, double *in, double *out)
{
    size_t nBytes = fft->size * sizeof (*fft->in);
    memcpy (fft->in, in, nBytes);
    fftw_execute (fft->plan);
    memcpy (out, fft->out, nBytes);
}
