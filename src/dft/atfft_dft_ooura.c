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
#include <math.h>
#include <atfft/atfft_dft.h>
#include "../ooura/ooura.h"

#ifdef ATFFT_TYPE_LONG_DOUBLE
#   ifdef _MSC_VER
#       pragma message(": warning: Ooura only supports double precision floating point, " \
                       "higher precision values will be demoted to double for FFT calculations.")
#   else
#       warning Ooura only supports double precision floating point, \
                higher precision values will be demoted to double for FFT calculations.
#   endif
#endif

struct atfft_dft
{
    int size;
    enum atfft_direction direction;
    int oouraDirection;
    enum atfft_format format;
    int dataSize;
    double *data;
    int *workArea;
    double *tables;
};

struct atfft_dft* atfft_dft_create (int size, enum atfft_direction direction, enum atfft_format format)
{
    struct atfft_dft *fft;
    int workSize;

    /* Ooura only supports sizes which are a power of 2. */
    assert (atfft_is_power_of_2 (size));

    if (!(fft = malloc (sizeof (*fft))))
        return NULL;

    fft->size = size;
    fft->direction = direction;
    fft->format = format;

    if (direction == ATFFT_FORWARD)
        fft->oouraDirection = -1;
    else
        fft->oouraDirection = 1;

    switch (format)
    {
        case ATFFT_COMPLEX:
            fft->dataSize = 2 * size * sizeof (*(fft->data));
            workSize = (2 + (1 << (int) (log (size + 0.5) / log (2)) / 2)) * sizeof (*(fft->workArea));
            break;

        case ATFFT_REAL:
            fft->oouraDirection *= -1; /* dodgy inconsistent interface */
            fft->dataSize = size * sizeof (*(fft->data));
            workSize = (2 + (1 << (int) (log (size / 2 + 0.5) / log (2)) / 2)) * sizeof (*(fft->workArea));
    }

    fft->data = malloc (fft->dataSize);
    fft->workArea = malloc (workSize);
    fft->tables = malloc (size / 2 * sizeof (*(fft->tables)));

    /* clean up on failure */
    if (!(fft->data && fft->workArea && fft->tables))
    {
        atfft_dft_destroy (fft);
        fft = NULL;
    }
    else
    {
        fft->workArea [0] = 0;
    }

    return fft;
}

void atfft_dft_destroy (struct atfft_dft *fft)
{
    if (fft)
    {
        free (fft->tables);
        free (fft->workArea);
        free (fft->data);
        free (fft);
    }
}

void atfft_dft_complex_transform (struct atfft_dft *fft, atfft_complex *in, atfft_complex *out)
{
    /* Only to be used with complex FFTs. */
    assert (fft->format == ATFFT_COMPLEX);

#ifdef ATFFT_TYPE_DOUBLE
    memcpy (fft->data, in, fft->dataSize);
#else
    atfft_sample_to_double_complex (in, (atfft_complex_d*) fft->data, fft->size);
#endif

    cdft (2 * fft->size, fft->oouraDirection, fft->data, fft->workArea, fft->tables);

#ifdef ATFFT_TYPE_DOUBLE
    memcpy (out, fft->data, fft->dataSize);
#else
    atfft_double_to_sample_complex ((atfft_complex_d*) fft->data, out, fft->size);
#endif
}

void atfft_halfcomplex_ooura_to_fftw (const double *in, atfft_complex *out, int size)
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

void atfft_dft_real_forward_transform (struct atfft_dft *fft, const atfft_sample *in, atfft_complex *out)
{
    /* Only to be used for forward real FFTs. */
    assert ((fft->format == ATFFT_REAL) && (fft->direction == ATFFT_FORWARD));

#ifdef ATFFT_TYPE_DOUBLE
    memcpy (fft->data, in, fft->dataSize);
#else
    atfft_sample_to_double_real (in, fft->data, fft->size);
#endif

    rdft (fft->size, fft->oouraDirection, fft->data, fft->workArea, fft->tables);
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

void atfft_dft_real_backward_transform (struct atfft_dft *fft, atfft_complex *in, atfft_sample *out)
{
    /* Only to be used for backward real FFTs. */
    assert ((fft->format == ATFFT_REAL) && (fft->direction == ATFFT_BACKWARD));

    atfft_halfcomplex_fftw_to_ooura (in, fft->data, fft->size);
    rdft (fft->size, fft->oouraDirection, fft->data, fft->workArea, fft->tables);

#ifdef ATFFT_TYPE_DOUBLE
    memcpy (out, fft->data, fft->dataSize);
#else
    atfft_double_to_sample_real (fft->data, out, fft->size);
#endif
}
