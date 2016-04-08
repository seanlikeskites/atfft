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
    kiss_fft_cfg cfg;
};

struct atfft* atfft_create (int size, int direction)
{
    struct atfft *fft;

    fft = malloc (sizeof (*fft));
    fft->size = size;
    fft->direction = direction;
    fft->cfg = kiss_fft_alloc (size, direction == 1, 0, 0);

    return fft;
}

void atfft_free (struct atfft *fft)
{
    free (fft->cfg);
    free (fft);
}

void atfft_complex_transform (struct atfft *fft, double *in, double *out)
{
    kiss_fft (fft->cfg, (kiss_fft_cpx*) in, (kiss_fft_cpx*) out);
}
