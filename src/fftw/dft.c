/*
 * Copyright (C) 2020 Sean Enderby <sean.enderby@gmail.com>
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
#include <atfft/atfft_dft.h>
#include "fftw_definitions.h"

struct atfft_dft
{
    int size;
    enum atfft_direction direction;
    enum atfft_format format;
    int in_size, out_size;
    atfft_sample *in, *out;
    atfft_fftw_plan plan;
};

struct atfft_dft* atfft_dft_create (int size, enum atfft_direction direction, enum atfft_format format)
{
    struct atfft_dft *fft;

    if (!(fft = malloc (sizeof (*fft))))
        return NULL;

    fft->size = size;
    fft->direction = direction;
    fft->format = format;

    switch (format)
    {
        case ATFFT_COMPLEX:
            fft->in_size = 2 * size * sizeof (*(fft->in));
            fft->in = ATFFT_FFTW_MALLOC (fft->in_size);

            fft->out_size = 2 * size * sizeof (*(fft->out));
            fft->out = ATFFT_FFTW_MALLOC (fft->out_size);

            if (direction == ATFFT_FORWARD)
                fft->plan = ATFFT_FFTW_PLAN_DFT_1D (size,
                                                    (atfft_complex*) fft->in,
                                                    (atfft_complex*) fft->out, 
                                                    FFTW_FORWARD, 
                                                    FFTW_ESTIMATE);
            else
                fft->plan = ATFFT_FFTW_PLAN_DFT_1D (size,
                                                    (atfft_complex*) fft->in,
                                                    (atfft_complex*) fft->out, 
                                                    FFTW_BACKWARD, 
                                                    FFTW_ESTIMATE);
            break;

        case ATFFT_REAL:
            if (direction == ATFFT_FORWARD)
            {
                fft->in_size = size * sizeof (*(fft->in));
                fft->in = ATFFT_FFTW_MALLOC (fft->in_size);

                fft->out_size = 2 * (floor (size / 2) + 1) * sizeof (*(fft->out));
                fft->out = ATFFT_FFTW_MALLOC (fft->out_size);

                fft->plan = ATFFT_FFTW_PLAN_DFT_R2C_1D (size, 
                                                        fft->in,
                                                        (atfft_complex*) fft->out,
                                                        FFTW_ESTIMATE);
            }
            else
            {
                fft->in_size = 2 * (floor (size / 2) + 1) * sizeof (*(fft->in));
                fft->in = ATFFT_FFTW_MALLOC (fft->in_size);

                fft->out_size = size * sizeof (*(fft->out));
                fft->out = ATFFT_FFTW_MALLOC (fft->out_size);

                fft->plan = ATFFT_FFTW_PLAN_DFT_C2R_1D (size,
                                                        (atfft_complex*) fft->in,
                                                        fft->out,
                                                        FFTW_ESTIMATE);
            }
            break;
    }

    /* clean up on failure */
    if (!(fft->in && fft->out && fft->plan))
    {
        atfft_dft_destroy (fft);
        fft = NULL;
    }

    return fft;
}

void atfft_dft_destroy (struct atfft_dft *fft)
{
    if (fft)
    {
        ATFFT_FFTW_DESTROY_PLAN (fft->plan);
        ATFFT_FFTW_FREE (fft->out);
        ATFFT_FFTW_FREE (fft->in);
        free (fft);
    }
}

static void atfft_dft_fftw_apply_transform (struct atfft_dft *fft, const atfft_sample *in, atfft_sample *out)
{
    memcpy (fft->in, in, fft->in_size);
    ATFFT_FFTW_EXECUTE (fft->plan);
    memcpy (out, fft->out, fft->out_size);
}

void atfft_dft_complex_transform (struct atfft_dft *fft, atfft_complex *in, atfft_complex *out)
{
    /* Only to be used with complex FFTs. */
    assert (fft->format == ATFFT_COMPLEX);

    atfft_dft_fftw_apply_transform (fft, (atfft_sample*) in, (atfft_sample*) out);
}

void atfft_dft_real_forward_transform (struct atfft_dft *fft, const atfft_sample *in, atfft_complex *out)
{
    /* Only to be used for forward real FFTs. */
    assert ((fft->format == ATFFT_REAL) && (fft->direction == ATFFT_FORWARD));

    atfft_dft_fftw_apply_transform (fft, in, (atfft_sample*) out);
}

void atfft_dft_real_backward_transform (struct atfft_dft *fft, atfft_complex *in, atfft_sample *out)
{
    /* Only to be used for backward real FFTs. */
    assert ((fft->format == ATFFT_REAL) && (fft->direction == ATFFT_BACKWARD));

    atfft_dft_fftw_apply_transform (fft, (atfft_sample*) in, out);
}
