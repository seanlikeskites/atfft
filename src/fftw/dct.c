/*
 * Copyright (C) 2020 Sean Enderby <sean.enderby@gmail.com>
 *
 * This program is free software. It comes without any warranty, to
 * the extent permitted by applicable law. You can redistribute it 
 * and/or modify it under the terms of the Do What The Fuck You Want 
 * To Public License, Version 2, as published by Sam Hocevar. See 
 * the COPYING file for more details.
 */

#include <fftw3.h>
#include <stdlib.h>
#include <string.h>
#include <atfft/dct.h>
#include "fftw_definitions.h"

struct atfft_dct
{
    int size;
    enum atfft_direction direction;
    int data_size;
    atfft_sample *in, *out;
    atfft_fftw_plan plan;
};

struct atfft_dct* atfft_dct_create (int size, enum atfft_direction direction)
{
    struct atfft_dct *dct;

    if (!(dct = malloc (sizeof (*dct))))
        return NULL;

    dct->size = size;
    dct->direction = direction;
    dct->data_size = size * sizeof (*(dct->in));
    dct->in = ATFFT_FFTW_MALLOC (dct->data_size);
    dct->out = ATFFT_FFTW_MALLOC (dct->data_size);

    switch (direction)
    {
        case ATFFT_FORWARD:
            dct->plan = ATFFT_FFTW_PLAN_R2R_1D (size,
                                                dct->in,
                                                dct->out,
                                                FFTW_REDFT10,
                                                FFTW_ESTIMATE);
            break;

        case ATFFT_BACKWARD:
            dct->plan = ATFFT_FFTW_PLAN_R2R_1D (size,
                                                dct->in,
                                                dct->out,
                                                FFTW_REDFT01,
                                                FFTW_ESTIMATE);
            break;
    };

    /* clean up on failure */
    if (!(dct->in && dct->out && dct->plan))
    {
        atfft_dct_destroy (dct);
        dct = NULL;
    }

    return dct;
}

void atfft_dct_destroy (struct atfft_dct *dct)
{
    if (dct)
    {
        ATFFT_FFTW_DESTROY_PLAN (dct->plan);
        ATFFT_FFTW_FREE (dct->out);
        ATFFT_FFTW_FREE (dct->in);
        free (dct);
    }
}

void atfft_dct_transform (struct atfft_dct *dct, const atfft_sample *in, atfft_sample *out)
{
    memcpy (dct->in, in, dct->data_size);
    ATFFT_FFTW_EXECUTE (dct->plan);
    memcpy (out, dct->out, dct->data_size);

    if (dct->direction == ATFFT_FORWARD)
        atfft_scale_real (out, dct->size, 0.5);
}
