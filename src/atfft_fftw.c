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
#include <fftw3.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <atfft.h>

struct atfft
{
    int size;
    int direction;
    enum atfft_format format;
    int inSize, outSize;
    double *in, *out;
    fftw_plan plan;
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
            fft->inSize = 2 * size * sizeof (*(fft->in));
            fft->in = fftw_malloc (fft->inSize);

            fft->outSize = 2 * size * sizeof (*(fft->out));
            fft->out = fftw_malloc (fft->outSize);

            fft->plan = fftw_plan_dft_1d (size, (fftw_complex*) fft->in, (fftw_complex*) fft->out, direction, FFTW_ESTIMATE);
            break;

        case ATFFT_REAL:
            if (direction == ATFFT_FORWARD)
            {
                fft->inSize = size * sizeof (*(fft->in));
                fft->in = fftw_malloc (fft->inSize);

                fft->outSize = 2 * (floor (size / 2) + 1) * sizeof (*(fft->out));
                fft->out = fftw_malloc (fft->outSize);

                fft->plan = fftw_plan_dft_r2c_1d (size, fft->in, (fftw_complex*) fft->out, FFTW_ESTIMATE);
            }
            else
            {
                fft->inSize = 2 * (floor (size / 2) + 1) * sizeof (*(fft->in));
                fft->in = fftw_malloc (fft->inSize);

                fft->outSize = size * sizeof (*(fft->out));
                fft->out = fftw_malloc (fft->outSize);

                fft->plan = fftw_plan_dft_c2r_1d (size, (fftw_complex*) fft->in, fft->out, FFTW_ESTIMATE);
            }
            break;
    }

    return fft;
}

void atfft_free (struct atfft *fft)
{
    fftw_destroy_plan (fft->plan);
    fftw_free (fft->out);
    fftw_free (fft->in);
    free (fft);
}

void atfft_fftw_apply_transform (struct atfft *fft, double *in, double *out)
{
    memcpy (fft->in, in, fft->inSize);
    fftw_execute (fft->plan);
    memcpy (out, fft->out, fft->outSize);
}

void atfft_complex_transform (struct atfft *fft, struct atfft_complex *in, struct atfft_complex *out)
{
    /* Only to be used with complex FFTs. */
    assert (fft->format == ATFFT_COMPLEX);

    atfft_fftw_apply_transform (fft, (double*) in, (double*) out);
}

void atfft_real_forward_transform (struct atfft *fft, double *in, struct atfft_complex *out)
{
    /* Only to be used for forward real FFTs. */
    assert ((fft->format == ATFFT_REAL) && (fft->direction == ATFFT_FORWARD));

    atfft_fftw_apply_transform (fft, in, (double*) out);
}

void atfft_real_backward_transform (struct atfft *fft, struct atfft_complex *in, double *out)
{
    /* Only to be used for backward real FFTs. */
    assert ((fft->format == ATFFT_REAL) && (fft->direction == ATFFT_BACKWARD));

    atfft_fftw_apply_transform (fft, (double*) in, out);
}
