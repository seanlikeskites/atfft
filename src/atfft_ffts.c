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
#include <assert.h>
#include <ffts/ffts.h>
#include <atfft.h>

struct atfft
{
    int size;
    int direction;
    enum atfft_format format;
    float *in, *out;
    ffts_plan_t *plan;
};

struct atfft* atfft_create (int size, int direction, enum atfft_format format)
{
    struct atfft *fft;

    /* FFTS only supports sizes which are a power of 2. */
    assert (atfft_is_power_of_2 (size));

    fft = malloc (sizeof (*fft));
    fft->size = size;
    fft->direction = direction;
    fft->format = format;

    switch (format)
    {
        case ATFFT_COMPLEX:
            fft->in = malloc (2 * size * sizeof (*(fft->in)));
            fft->out = malloc (2 * size * sizeof (*(fft->out)));
            fft->plan = ffts_init_1d (size, direction);
            break;

        case ATFFT_REAL:
            if (direction == ATFFT_FORWARD)
            {
                fft->in = malloc (size * sizeof (*(fft->in)));
                fft->out = malloc (2 * (floor (size / 2) + 1) * sizeof (*(fft->out)));
            }
            else
            {
                fft->in = malloc (2 * (floor (size / 2) + 1) * sizeof (*(fft->in)));
                fft->out = malloc (size * sizeof (*(fft->out)));
            }

            fft->plan = ffts_init_1d_real (size, direction);
            break;
    }

    return fft;
}

void atfft_free (struct atfft *fft)
{
    ffts_free (fft->plan);
    free (fft->out);
    free (fft->in);
    free (fft);
}

void atfft_complex_transform (struct atfft *fft, atfft_complex_double *in, atfft_complex_double *out)
{
    /* Only to be used with complex FFTs. */
    assert (fft->format == ATFFT_COMPLEX);

    atfft_double_to_float_complex (in, (atfft_complex_float*) fft->in, fft->size);
    ffts_execute (fft->plan, fft->in, fft->out);
    atfft_float_to_double_complex ((atfft_complex_float*) fft->out, out, fft->size);
}

void atfft_real_forward_transform (struct atfft *fft, double *in, atfft_complex_double *out)
{
    /* Only to be used for forward real FFTs. */
    assert ((fft->format == ATFFT_REAL) && (fft->direction == ATFFT_FORWARD));

    atfft_double_to_float_real (in, fft->in, fft->size);
    ffts_execute (fft->plan, fft->in, fft->out);
    atfft_float_to_double_complex ((atfft_complex_float*) fft->out, out, fft->size / 2 + 1);
}

void atfft_real_backward_transform (struct atfft *fft, atfft_complex_double *in, double *out)
{
    /* Only to be used for backward real FFTs. */
    assert ((fft->format == ATFFT_REAL) && (fft->direction == ATFFT_BACKWARD));

    atfft_double_to_float_complex (in, (atfft_complex_float*) fft->in, fft->size / 2 + 1);
    ffts_execute (fft->plan, fft->in, fft->out);
    atfft_float_to_double_real (fft->out, out, fft->size);
}
