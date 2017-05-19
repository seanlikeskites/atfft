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
#include <Accelerate/Accelerate.h>
#include <atfft/atfft_dct.h>

#ifndef ATFFT_TYPE_FLOAT
#   warning vDSP DCTs only support single precision floating point, \
            higher precision values will be demoted to float for DCT calculations.
#endif

struct atfft_dct
{
    int size;
    enum atfft_direction direction;
#ifndef ATFFT_TYPE_FLOAT
    float *in, *out;
#endif
    vDSP_DFT_Setup setup;
};

int atfft_is_supported_length_vdsp_dct (unsigned int length)
{
    int min = 16;

    if (atfft_is_power_of_2 (length))
        return 1;

    if ((!(length % 3) && atfft_is_power_of_2 (length / 3) && (length / 3 >= min)) ||
        (!(length % 5) && atfft_is_power_of_2 (length / 5) && (length / 5 >= min)))
        return 1;        
    else
        return 0;
}

struct atfft_dct* atfft_dct_create (int size, enum atfft_direction direction)
{
    struct atfft_dct *dct;

    /* vDSP only supports certain lengths */
    assert (atfft_is_supported_length_vdsp_dct (size));

    if (!(dct = malloc (sizeof (*dct))))
        return NULL;

    dct->size = size;
    dct->direction = direction;

#ifndef ATFFT_TYPE_FLOAT
    dct->in = malloc (size * sizeof (*(dct->in)));
    dct->out = malloc (size * sizeof (*(dct->out)));
#endif

    if (direction == ATFFT_FORWARD)
        dct->setup = vDSP_DCT_CreateSetup (NULL, size, vDSP_DCT_II);
    else
        dct->setup = vDSP_DCT_CreateSetup (NULL, size, vDSP_DCT_III);

#ifndef ATFFT_TYPE_FLOAT
    if (!(dct->in && dct->out && dct->setup))
#else
    if (!dct->setup)
#endif
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
        vDSP_DFT_DestroySetup (dct->setup);
#ifndef ATFFT_TYPE_FLOAT
        free (dct->out);
        free (dct->in);
#endif
        free (dct);
    }
}

void atfft_dct_transform (struct atfft_dct *dct, const atfft_sample *in, atfft_sample *out)
{
#ifdef ATFFT_TYPE_FLOAT
    vDSP_DCT_Execute (dct->setup, in, out);
#else
    atfft_sample_to_float_real (in, dct->in, dct->size);
    vDSP_DCT_Execute (dct->setup, dct->in, dct->out);
    atfft_float_to_sample_real (dct->out, out, dct->size);
#endif

    if (dct->direction == ATFFT_BACKWARD)
        atfft_scale_real (out, dct->size, 2.0);
}
