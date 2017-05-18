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
#include <math.h>
#include <atfft/atfft_dft.h>
#include <atfft/atfft_dct.h>

struct atfft_dct
{
    int size;
    enum atfft_direction direction;
    struct atfft_dft *dft;
    atfft_sample *cosins, *sins;
    atfft_complex *in, *out;
};

void atfft_dct_init_sins (atfft_sample *cosins, atfft_sample *sins, int size)
{
    int i = 0;

    for (i = 0; i < size; ++i)
    {
        atfft_sample x = i * M_PI / (2.0 * size);
        cosins [i] = cos (x);
        sins [i] = sin (x);
    }
}

struct atfft_dct* atfft_dct_create (int size, enum atfft_direction direction)
{
    struct atfft_dct *dct;

    if (!(dct = malloc (sizeof (*dct))))
        return NULL;

    dct->size = size;
    dct->direction = direction;
    dct->dft = atfft_dft_create (size, direction, ATFFT_COMPLEX);
    dct->cosins = malloc (size * sizeof (*(dct->cosins)));
    dct->sins = malloc (size * sizeof (*(dct->sins)));
    dct->in = malloc (size * sizeof (*(dct->in)));
    dct->out = malloc (size * sizeof (*(dct->out)));

    /* clean up on failure */
    if (!(dct->dft && dct->cosins && dct->sins && dct->in && dct->out))
    {
        atfft_dct_destroy (dct);
        dct = NULL;
    }
    else
    {
        atfft_dct_init_sins (dct->cosins, dct->sins, size);
    }

    return dct;
}

void atfft_dct_destroy (struct atfft_dct *dct)
{
    if (dct)
    {
        free (dct->out);
        free (dct->in);
        free (dct->sins);
        free (dct->cosins);
        atfft_dft_destroy (dct->dft);
        free (dct);
    }
}

void atfft_dct_rearrange_forward (const atfft_sample *in, atfft_complex *out, int size)
{
    int i = 0, j = 0;
    int start = 0;

    for (i = 0, j = 0; i < size; i+=2, ++j)
    {
        ATFFT_REAL (out [j]) = in [i];
        ATFFT_IMAG (out [j]) = 0.0;
    }

    if (atfft_is_even (size))
        start = size - 1;
    else
        start = size - 2;

    for (i = start; i > 0; i-=2, ++j)
    {
        ATFFT_REAL (out [j]) = in [i];
        ATFFT_IMAG (out [j]) = 0.0;
    }
}

void atfft_dct_scale_forward (struct atfft_dct *dct, atfft_sample *out)
{
    int i = 0;

    for (i = 0; i < dct->size; ++i)
    {
        atfft_sample cosComponent = ATFFT_REAL (dct->out [i]) * dct->cosins [i];
        atfft_sample sinComponent = ATFFT_IMAG (dct->out [i]) * dct->sins [i];
        out [i] = cosComponent + sinComponent;
    }
}

void atfft_dct_forward_transform (struct atfft_dct *dct, const atfft_sample *in, atfft_sample *out)
{
    atfft_dct_rearrange_forward (in, dct->in, dct->size);
    atfft_dft_complex_transform (dct->dft, dct->in, dct->out);
    atfft_dct_scale_forward (dct, out);
}

void atfft_dct_scale_backward (struct atfft_dct *dct, const atfft_sample *in)
{
    int i = 0;

    ATFFT_REAL (dct->in [0]) = in [0];
    ATFFT_IMAG (dct->in [0]) = 0.0;

    for (i = 1; i < dct->size; ++i)
    {
        atfft_sample realPart = in [i];
        atfft_sample imagPart = -in [dct->size - i];

        atfft_sample realCosComponent = dct->cosins [i] * realPart;
        atfft_sample realSinComponent = - dct->sins [i] * imagPart;

        atfft_sample imagCosComponent = dct->sins [i] * realPart;
        atfft_sample imagSinComponent = dct->cosins [i] * imagPart;

        ATFFT_REAL (dct->in [i]) = realCosComponent + realSinComponent;
        ATFFT_IMAG (dct->in [i]) = imagCosComponent + imagSinComponent;
    }
}

void atfft_dct_rearrange_backward (const atfft_complex *in, atfft_sample *out, int size)
{
    int i = 0, j = 0;
    int start = 0;

    for (i = 0, j = 0; i < size; i+=2, ++j)
    {
        out [i] = ATFFT_REAL (in [j]);
    }

    if (atfft_is_even (size))
        start = size - 1;
    else
        start = size - 2;

    for (i = start; i > 0; i-=2, ++j)
    {
        out [i] = ATFFT_REAL (in [j]);
    }
}

void atfft_dct_backward_transform (struct atfft_dct *dct, const atfft_sample *in, atfft_sample *out)
{
    atfft_dct_scale_backward (dct, in);
    atfft_dft_complex_transform (dct->dft, dct->in, dct->out);
    atfft_dct_rearrange_backward (dct->out, out, dct->size);
}

void atfft_dct_transform (struct atfft_dct *dct, const atfft_sample *in, atfft_sample *out)
{
    switch (dct->direction)
    {
        case ATFFT_FORWARD:
            atfft_dct_forward_transform (dct, in, out);
            break;

        case ATFFT_BACKWARD:
            atfft_dct_backward_transform (dct, in, out);
            break;
    };
}
