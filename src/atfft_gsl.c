/*
 * Copyright (C) 2016 Sean Enderby <sean.enderby@gmail.com>
 *
 * This work is free. You can redistribute it and/or modify it under the
 * terms of the Do What The Fuck You Want To Public License, Version 2,
 * as published by Sam Hocevar. See the COPYING file for more details.
 */

#include <math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_fft_complex.h>
#include <stdlib.h>
#include <string.h>
#include "../include/atfft.h"

struct atfft
{
    int size;
    int direction;
    gsl_vector_complex *data;
    gsl_fft_complex_workspace *workArea;
    gsl_fft_complex_wavetable *tables;
};

struct atfft* atfft_create (int size, int direction, enum atfft_format format)
{
    struct atfft *fft;

    fft = malloc (sizeof (*fft));
    fft->size = size;
    fft->direction = direction;
    fft->data = gsl_vector_complex_alloc (size);
    fft->workArea = gsl_fft_complex_workspace_alloc (size);
    fft->tables = gsl_fft_complex_wavetable_alloc (size);

    return fft;
}

void atfft_free (struct atfft *fft)
{
    gsl_fft_complex_wavetable_free (fft->tables);
    gsl_fft_complex_workspace_free (fft->workArea);
    gsl_vector_complex_free (fft->data);
    free (fft);
}

void atfft_complex_transform (struct atfft *fft, double *in, double *out)
{
    size_t nBytes = 2 * fft->size * sizeof (*in);
    memcpy (fft->data->data, in, nBytes);
    gsl_fft_complex_transform (fft->data->data, fft->data->stride, fft->data->size, fft->tables, fft->workArea, fft->direction);
    memcpy (out, fft->data->data, nBytes);
}
