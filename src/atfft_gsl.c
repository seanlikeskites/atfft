/*
 * Copyright (C) 2016 Sean Enderby <sean.enderby@gmail.com>
 *
 * This program is free software. It comes without any warranty, to
 * the extent permitted by applicable law. You can redistribute it 
 * and/or modify it under the terms of the Do What The Fuck You Want 
 * To Public License, Version 2, as published by Sam Hocevar. See 
 * the COPYING file for more details.
 */

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_fft_complex.h>
#include <gsl/gsl_fft_real.h>
#include <gsl/gsl_fft_halfcomplex.h>
#include <atfft.h>

struct atfft
{
    int size;
    int direction;
    enum atfft_format format;
    gsl_vector *data;
    void *workArea;
    void *tables;
};

struct atfft* atfft_create (int size, int direction, enum atfft_format format)
{
    struct atfft *fft;

    fft = malloc (sizeof (*fft));
    fft->size = size;
    fft->direction = direction;
    fft->format = format;

    switch (format)
    {
        case ATFFT_COMPLEX:
            fft->data = gsl_vector_alloc (2 * size);
            fft->workArea = gsl_fft_complex_workspace_alloc (size);
            fft->tables = gsl_fft_complex_wavetable_alloc (size);
            break;

        case ATFFT_REAL:
            fft->data = gsl_vector_alloc (size);
            fft->workArea = gsl_fft_real_workspace_alloc (size);

            if (direction == ATFFT_FORWARD)
            {
                fft->tables = gsl_fft_real_wavetable_alloc (size);
            }
            else
            {
                fft->tables = gsl_fft_halfcomplex_wavetable_alloc (size);
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
            gsl_fft_complex_wavetable_free (fft->tables);
            gsl_fft_complex_workspace_free (fft->workArea);
            break;

        case ATFFT_REAL:
            if (fft->direction == ATFFT_FORWARD)
            {
                gsl_fft_real_wavetable_free (fft->tables);
            }
            else
            {
                gsl_fft_halfcomplex_wavetable_free (fft->tables);
            }

            gsl_fft_real_workspace_free (fft->workArea);
            break;
    }

    gsl_vector_free (fft->data);

    free (fft);
}

void atfft_complex_transform (struct atfft *fft, struct atfft_complex *in, struct atfft_complex *out)
{
    size_t nBytes = fft->size * sizeof (*in);
    gsl_complex_packed_array data = fft->data->data;

    /* Only to be used with complex FFTs. */
    assert (fft->format == ATFFT_COMPLEX);

    memcpy (data, in, nBytes);
    gsl_fft_complex_transform (data, fft->data->stride, fft->size, fft->tables, fft->workArea, fft->direction);
    memcpy (out, data, nBytes);
}

void atfft_halfcomplex_gsl_to_fftw (double *in, struct atfft_complex *out, int size)
{
    out [0].real = in [0];
    out [0].imag = 0;

    memcpy (out + 1, in + 1, (size - 1) * sizeof (*in));

    if (isEven (size))
    {
        out [size / 2].imag = 0;
    }
}

void atfft_real_forward_transform (struct atfft *fft, double *in, struct atfft_complex *out)
{
    double *data = fft->data->data;

    /* Only to be used for forward real FFTs. */
    assert ((fft->format == ATFFT_REAL) && (fft->direction == ATFFT_FORWARD));

    memcpy (data, in, fft->size * sizeof (*data));
    gsl_fft_real_transform (data, fft->data->stride, fft->size, fft->tables, fft->workArea);
    atfft_halfcomplex_gsl_to_fftw (data, out, fft->size);
}

void atfft_halfcomplex_fftw_to_gsl (struct atfft_complex *in, double *out, int size)
{
    out [0] = in [0].real;
    memcpy (out + 1, in + 1, (size - 1) * sizeof (*out));
}

void atfft_real_backward_transform (struct atfft *fft, struct atfft_complex *in, double *out)
{
    double *data = fft->data->data;

    /* Only to be used for backward real FFTs. */
    assert ((fft->format == ATFFT_REAL) && (fft->direction == ATFFT_BACKWARD));

    atfft_halfcomplex_fftw_to_gsl (in, data, fft->size);
    gsl_fft_halfcomplex_transform (data, fft->data->stride, fft->size, fft->tables, fft->workArea);   
    memcpy (out, data, fft->size * sizeof (*data));
}
