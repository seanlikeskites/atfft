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
#include <atfft/atfft_dft.h>
#include <atfft/atfft_dct.h>

struct atfft_dct
{
    int size;
    enum atfft_direction direction;
    struct atfft_dft *dft;
};

struct atfft_dct* atfft_dct_create (int size, enum atfft_direction direction)
{
    struct atfft_dct *dct;

    if (!(dct = malloc (sizeof (*dct))))
        return NULL;

    dct->size = size;
    dct->direction = direction;
    dct->dft = atfft_dft_create (size, direction, ATFFT_REAL);

    /* clean up on failure */
    if (!dct->dft)
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
        atfft_dft_destroy (dct->dft);
        free (dct);
    }
}

void atfft_dct_forward_transform (struct atfft_dct *dct, const atfft_sample *in, atfft_sample *out)
{
}

void atfft_dct_backward_transform (struct atfft_dct *dct, const atfft_sample *in, atfft_sample *out)
{
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
