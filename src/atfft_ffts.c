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
#ifndef ATFFT_TYPE_FLOAT
    float *in, *out;
#endif
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
#ifndef ATFFT_TYPE_FLOAT
            fft->in = malloc (2 * size * sizeof (*(fft->in)));
            fft->out = malloc (2 * size * sizeof (*(fft->out)));
#endif
            fft->plan = ffts_init_1d (size, direction);
            break;

        case ATFFT_REAL:
#ifndef ATFFT_TYPE_FLOAT
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
#endif

            fft->plan = ffts_init_1d_real (size, direction);
            break;
    }

    return fft;
}

void atfft_destroy (struct atfft *fft)
{
    ffts_free (fft->plan);
#ifndef ATFFT_TYPE_FLOAT
    free (fft->out);
    free (fft->in);
#endif
    free (fft);
}

void atfft_complex_transform (struct atfft *fft, atfft_complex *in, atfft_complex *out)
{
    /* Only to be used with complex FFTs. */
    assert (fft->format == ATFFT_COMPLEX);

#ifdef ATFFT_TYPE_FLOAT
    ffts_execute (fft->plan, (float*) in, (float*) out);
#else
    atfft_sample_to_float_complex (in, (atfft_complex_f*) fft->in, fft->size);
    ffts_execute (fft->plan, fft->in, fft->out);
    atfft_float_to_sample_complex ((atfft_complex_f*) fft->out, out, fft->size);
#endif
}

void atfft_real_forward_transform (struct atfft *fft, atfft_sample *in, atfft_complex *out)
{
    /* Only to be used for forward real FFTs. */
    assert ((fft->format == ATFFT_REAL) && (fft->direction == ATFFT_FORWARD));

#ifdef ATFFT_TYPE_FLOAT
    ffts_execute (fft->plan, (float*) in, (float*) out);
#else
    atfft_sample_to_float_real (in, fft->in, fft->size);
    ffts_execute (fft->plan, fft->in, fft->out);
    atfft_float_to_sample_complex ((atfft_complex_f*) fft->out, out, fft->size / 2 + 1);
#endif
}

void atfft_real_backward_transform (struct atfft *fft, atfft_complex *in, atfft_sample *out)
{
    /* Only to be used for backward real FFTs. */
    assert ((fft->format == ATFFT_REAL) && (fft->direction == ATFFT_BACKWARD));

#ifdef ATFFT_TYPE_FLOAT
    ffts_execute (fft->plan, (float*) in, (float*) out);
#else
    atfft_sample_to_float_complex (in, (atfft_complex_f*) fft->in, fft->size / 2 + 1);
    ffts_execute (fft->plan, fft->in, fft->out);
    atfft_float_to_sample_real (fft->out, out, fft->size);
#endif
}
