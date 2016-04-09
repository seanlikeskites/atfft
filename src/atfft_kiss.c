/*
 * Copyright (C) 2016 Sean Enderby <sean.enderby@gmail.com>
 *
 * This work is free. You can redistribute it and/or modify it under the
 * terms of the Do What The Fuck You Want To Public License, Version 2,
 * as published by Sam Hocevar. See the COPYING file for more details.
 */

/* we need to make sure we are using kiss_fft with doubles */
#define kiss_fft_scalar double

#include <stdlib.h>
#include <string.h>
#include "../include/atfft.h"
#include "kiss_fft.h"

struct atfft
{
    int size;
    int direction;
    enum atfft_format format;
    kiss_fft_cfg cfg;
    double *in, *out;
};

struct atfft* atfft_create (int size, int direction, enum atfft_format format)
{
    struct atfft *fft;

    fft = malloc (sizeof (*fft));
    fft->size = size;
    fft->direction = direction;
    fft->format = format;
    fft->cfg = kiss_fft_alloc (size, direction == 1, 0, 0);

    if (format == ATFFT_REAL)
    {
        fft->in = malloc (2 * size * sizeof (*(fft->in)));
        fft->out = malloc (2 * size * sizeof (*(fft->out)));
    }

    return fft;
}

void atfft_free (struct atfft *fft)
{
    if (fft->format == ATFFT_REAL)
    {
        free (fft->out);
        free (fft->in);
    }

    free (fft->cfg);
    free (fft);
}

void atfft_complex_transform (struct atfft *fft, double *in, double *out)
{
    kiss_fft (fft->cfg, (kiss_fft_cpx*) in, (kiss_fft_cpx*) out);
}

void atfft_complex_to_halfcomplex (double *in, double *out, int size)
{
    memcpy (out, in, 2 * (floor (size / 2) + 1) * sizeof (*out));
}

void atfft_real_forward_transform (struct atfft *fft, double *in, double *out)
{
    atfft_real_to_complex (in, fft->in, fft->size);
    atfft_complex_transform (fft, fft->in, fft->out);
    atfft_complex_to_halfcomplex (fft->out, out, fft->size);
}

void atfft_real_backward_transform (struct atfft *fft, double *in, double *out)
{
    atfft_halfcomplex_to_complex (in, fft->in, fft->size);
    atfft_complex_transform (fft, fft->in, fft->out);
    atfft_real (fft->out, out, fft->size);
}
