/*
 * Copyright (C) 2020 Sean Enderby <sean.enderby@gmail.com>
 *
 * This program is free software. It comes without any warranty, to
 * the extent permitted by applicable law. You can redistribute it 
 * and/or modify it under the terms of the Do What The Fuck You Want 
 * To Public License, Version 2, as published by Sam Hocevar. See 
 * the COPYING file for more details.
 */

#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <ffts/ffts.h>
#include <atfft/atfft_dft.h>

#ifndef ATFFT_TYPE_FLOAT
#   ifdef _MSC_VER
#       pragma message(": warning: FFTS only supports single precision floating point, " \
                       "higher precision values will be demoted to float for FFT calculations.")
#   else
#       warning FFTS only supports single precision floating point, \
                higher precision values will be demoted to float for FFT calculations.
#   endif
#endif

struct atfft_dft
{
    int size;
    enum atfft_direction direction;
    enum atfft_format format;
#ifndef ATFFT_TYPE_FLOAT
    float *in, *out;
#endif
    ffts_plan_t *plan;
};

struct atfft_dft* atfft_dft_create (int size, enum atfft_direction direction, enum atfft_format format)
{
    struct atfft_dft *fft;

    if (!(fft = malloc (sizeof (*fft))))
        return NULL;

    fft->size = size;
    fft->direction = direction;
    fft->format = format;

    int ffts_direction = 0;

    if (direction == ATFFT_FORWARD)
        ffts_direction = -1;
    else
        ffts_direction = 1;

#ifndef ATFFT_TYPE_FLOAT
    int in_size = 0;
    int out_size = 0;
#endif

    switch (format)
    {
        case ATFFT_COMPLEX:
#ifndef ATFFT_TYPE_FLOAT
            in_size = 2 * size * sizeof (*(fft->in));
            out_size = 2 * size * sizeof (*(fft->out));
#endif
            fft->plan = ffts_init_1d (size, ffts_direction);
            break;

        case ATFFT_REAL:
#ifndef ATFFT_TYPE_FLOAT
            if (direction == ATFFT_FORWARD)
            {
                in_size = size * sizeof (*(fft->in));
                out_size = 2 * (floor (size / 2) + 1) * sizeof (*(fft->out));
            }
            else
            {
                in_size = 2 * (floor (size / 2) + 1) * sizeof (*(fft->in));
                out_size = size * sizeof (*(fft->out));
            }
#endif

            fft->plan = ffts_init_1d_real (size, ffts_direction);
            break;
    }

#ifndef ATFFT_TYPE_FLOAT
    fft->in = malloc (in_size);
    fft->out = malloc (out_size);
#endif

#ifndef ATFFT_TYPE_FLOAT
    if (!(fft->in && fft->out && fft->plan))
    {
        if (fft->plan) /* ffts can't hack freeing a null pointer */
            ffts_free (fft->plan);

        free (fft->out);
        free (fft->in);
        free (fft);
        fft = NULL;
    }
#else
    if (!fft->plan)
    {
        free (fft);
        fft = NULL;
    }
#endif

    return fft;
}

void atfft_dft_destroy (struct atfft_dft *fft)
{
    if (fft)
    {
        ffts_free (fft->plan);
#ifndef ATFFT_TYPE_FLOAT
        free (fft->out);
        free (fft->in);
#endif
        free (fft);
    }
}

void atfft_dft_complex_transform (struct atfft_dft *fft, atfft_complex *in, atfft_complex *out)
{
    /* Only to be used with complex FFTs. */
    assert (fft->format == ATFFT_COMPLEX);

#ifdef ATFFT_TYPE_FLOAT
    ffts_execute (fft->plan, (const float*) in, (float*) out);
#else
    atfft_sample_to_float_complex (in, (atfft_complex_f*) fft->in, fft->size);
    ffts_execute (fft->plan, fft->in, fft->out);
    atfft_float_to_sample_complex ((atfft_complex_f*) fft->out, out, fft->size);
#endif
}

void atfft_dft_real_forward_transform (struct atfft_dft *fft, const atfft_sample *in, atfft_complex *out)
{
    /* Only to be used for forward real FFTs. */
    assert ((fft->format == ATFFT_REAL) && (fft->direction == ATFFT_FORWARD));

#ifdef ATFFT_TYPE_FLOAT
    ffts_execute (fft->plan, (const float*) in, (float*) out);
#else
    atfft_sample_to_float_real (in, fft->in, fft->size);
    ffts_execute (fft->plan, fft->in, fft->out);
    atfft_float_to_sample_complex ((atfft_complex_f*) fft->out, out, fft->size / 2 + 1);
#endif
}

void atfft_dft_real_backward_transform (struct atfft_dft *fft, atfft_complex *in, atfft_sample *out)
{
    /* Only to be used for backward real FFTs. */
    assert ((fft->format == ATFFT_REAL) && (fft->direction == ATFFT_BACKWARD));

#ifdef ATFFT_TYPE_FLOAT
    ffts_execute (fft->plan, (const float*) in, (float*) out);
#else
    atfft_sample_to_float_complex (in, (atfft_complex_f*) fft->in, fft->size / 2 + 1);
    ffts_execute (fft->plan, fft->in, fft->out);
    atfft_float_to_sample_real (fft->out, out, fft->size);
#endif
}
