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

void atfft_complex_transform (struct atfft *fft, atfft_complex *in, atfft_complex *out)
{
#ifdef ATFFT_TYPE_DOUBLE
    size_t nBytes = fft->size * sizeof (*in);
#endif 

    gsl_complex_packed_array data = fft->data->data;

    /* Only to be used with complex FFTs. */
    assert (fft->format == ATFFT_COMPLEX);

#ifdef ATFFT_TYPE_DOUBLE
    memcpy (data, in, nBytes);
#else
    atfft_sample_to_double_complex (in, (atfft_complex_d*) data, fft->size);
#endif

    gsl_fft_complex_transform (data, fft->data->stride, fft->size, fft->tables, fft->workArea, fft->direction);

#ifdef ATFFT_TYPE_DOUBLE
    memcpy (out, data, nBytes);
#else
    atfft_double_to_sample_complex ((atfft_complex_d*) data, out, fft->size);
#endif
}

void atfft_halfcomplex_gsl_to_fftw (double *in, atfft_complex *out, int size)
{
    int i = 0;
    int halfSize = size / 2;

    ATFFT_REAL (out [0]) = in [0];
    ATFFT_IMAG (out [0]) = 0;

    for (i = 1; i < halfSize; ++i)
    {
        ATFFT_REAL (out [i]) = in [2 * i - 1];
        ATFFT_IMAG (out [i]) = in [2 * i];
    }

    ATFFT_REAL (out [halfSize]) = in [size - 1];
    ATFFT_IMAG (out [halfSize]) = 0;
}

void atfft_real_forward_transform (struct atfft *fft, atfft_sample *in, atfft_complex *out)
{
    double *data = fft->data->data;

    /* Only to be used for forward real FFTs. */
    assert ((fft->format == ATFFT_REAL) && (fft->direction == ATFFT_FORWARD));

#ifdef ATFFT_TYPE_DOUBLE
    memcpy (data, in, fft->size * sizeof (*data));
#else
    atfft_sample_to_double_real (in, data, fft->size);
#endif

    gsl_fft_real_transform (data, fft->data->stride, fft->size, fft->tables, fft->workArea);
    atfft_halfcomplex_gsl_to_fftw (data, out, fft->size);
}

void atfft_halfcomplex_fftw_to_gsl (atfft_complex *in, double *out, int size)
{
    int i = 0;
    int halfSize = size / 2;

    out [0] = ATFFT_REAL (in [0]);

    for (i = 1; i < halfSize; ++i)
    {
        out [2 * i - 1] = ATFFT_REAL (in [i]);
        out [2 * i] = ATFFT_IMAG (in [i]);
    }

    out [size - 1] = ATFFT_REAL (in [halfSize]);
}

void atfft_real_backward_transform (struct atfft *fft, atfft_complex *in, atfft_sample *out)
{
    double *data = fft->data->data;

    /* Only to be used for backward real FFTs. */
    assert ((fft->format == ATFFT_REAL) && (fft->direction == ATFFT_BACKWARD));

    atfft_halfcomplex_fftw_to_gsl (in, data, fft->size);
    gsl_fft_halfcomplex_transform (data, fft->data->stride, fft->size, fft->tables, fft->workArea);   

#ifdef ATFFT_TYPE_DOUBLE
    memcpy (out, data, fft->size * sizeof (*data));
#else
    atfft_double_to_sample_real (data, out, fft->size);
#endif
}
