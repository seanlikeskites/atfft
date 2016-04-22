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
#include <atfft.h>

/* we need to make sure we are using kiss_fft with the correct type */
#define kiss_fft_scalar atfft_sample
#include <kiss_fft.h>

struct atfft
{
    int size;
    int direction;
    enum atfft_format format;
    kiss_fft_cfg cfg;
    atfft_complex *in, *out;
};

struct atfft* atfft_create (int size, int direction, enum atfft_format format)
{
    struct atfft *fft;
    int dataOk = 0;

    if (!(fft = malloc (sizeof (*fft))))
        return NULL;

    fft->size = size;
    fft->direction = direction;
    fft->format = format;
    fft->cfg = kiss_fft_alloc (size, direction == 1, 0, 0);

    if (format == ATFFT_REAL)
    {
        fft->in = malloc (size * sizeof (*(fft->in)));
        fft->out = malloc (size * sizeof (*(fft->out)));
        dataOk = fft->in && fft->out;
    }
    else
    {
        fft->in = NULL;
        fft->out = NULL;
        dataOk = 1;
    }

    /* clean up on failure */
    if (!(fft->cfg && dataOk))
    {
        free (fft->in);
        free (fft->cfg);
        free (fft);
        fft = NULL;
    }

    return fft;
}

void atfft_destroy (struct atfft *fft)
{
    if (fft)
    {
        free (fft->out);
        free (fft->in);
        free (fft->cfg);
        free (fft);
    }
}

void atfft_complex_transform (struct atfft *fft, atfft_complex *in, atfft_complex *out)
{
    /* Only to be used with complex FFTs. */
    assert (fft->format == ATFFT_COMPLEX);

    kiss_fft (fft->cfg, (kiss_fft_cpx*) in, (kiss_fft_cpx*) out);
}

void atfft_complex_to_halfcomplex (atfft_complex *in, atfft_complex *out, int size)
{
    memcpy (out, in, (floor (size / 2) + 1) * sizeof (*out));
}

void atfft_real_forward_transform (struct atfft *fft, atfft_sample *in, atfft_complex *out)
{
    /* Only to be used for forward real FFTs. */
    assert ((fft->format == ATFFT_REAL) && (fft->direction == ATFFT_FORWARD));

    atfft_real_to_complex (in, fft->in, fft->size);
    fft->format = ATFFT_COMPLEX;
    atfft_complex_transform (fft, fft->in, fft->out);
    fft->format = ATFFT_REAL;
    atfft_complex_to_halfcomplex (fft->out, out, fft->size);
}

void atfft_real_backward_transform (struct atfft *fft, atfft_complex *in, atfft_sample *out)
{
    /* Only to be used for backward real FFTs. */
    assert ((fft->format == ATFFT_REAL) && (fft->direction == ATFFT_BACKWARD));

    atfft_halfcomplex_to_complex (in, fft->in, fft->size);
    fft->format = ATFFT_COMPLEX;
    atfft_complex_transform (fft, fft->in, fft->out);
    fft->format = ATFFT_REAL;
    atfft_real (fft->out, out, fft->size);
}
