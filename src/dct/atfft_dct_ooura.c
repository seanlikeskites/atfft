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
#include <math.h>
#include <assert.h>
#include <atfft/atfft_dct.h>
#include "../ooura/ooura.h"

#ifdef ATFFT_TYPE_LONG_DOUBLE
#   ifdef _MSC_VER
#       pragma message(": warning: Ooura only supports double precision floating point, " \
                       "higher precision values will be demoted to double for FFT calculations.")
#   else
#       warning Ooura only supports double precision floating point, \
                higher precision values will be demoted to double for FFT calculations.
#   endif
#endif

struct atfft_dct
{
    int size;
    enum atfft_direction direction;
    int oouraDirection;
    int dataSize;
    double *data;
    int *workArea;
    double *tables;
};

struct atfft_dct* atfft_dct_create (int size, enum atfft_direction direction)
{
    struct atfft_dct *dct;
    int workSize;

    /* Ooura only supports sizes which are a power of 2. */
    assert (atfft_is_power_of_2 (size));

    if (!(dct = malloc (sizeof (*dct))))
        return NULL;

    dct->size = size;
    dct->direction = direction;

    if (direction == ATFFT_FORWARD)
        dct->oouraDirection = -1;
    else
        dct->oouraDirection = 1;

    dct->dataSize = size * sizeof (*(dct->data));
    dct->data = malloc (dct->dataSize);

    workSize = (2 + (1 << (int) (log (size / 2 + 0.5) / log (2)) / 2)) * sizeof (*(dct->workArea));
    dct->workArea = malloc (workSize);

    dct->tables = malloc ((size * 5 / 4) * sizeof (*(dct->tables)));

    /* clean up on failure */
    if (!(dct->data && dct->workArea && dct->tables))
    {
        atfft_dct_destroy (dct);
        dct = NULL;
    }
    else
    {
        dct->workArea [0] = 0;
    }

    return dct;
}

void atfft_dct_destroy (struct atfft_dct *dct)
{
    if (dct)
    {
        free (dct->tables);
        free (dct->workArea);
        free (dct->data);
        free (dct);
    }
}

void atfft_dct_transform (struct atfft_dct *dct, const atfft_sample *in, atfft_sample *out)
{
#ifdef ATFFT_TYPE_DOUBLE
    memcpy (dct->data, in, dct->dataSize);
#else
    atfft_sample_to_double_real (in, dct->data, dct->size);
#endif

    if (dct->direction == ATFFT_BACKWARD)
        dct->data [0] *= 0.5;

    ddct (dct->size, dct->oouraDirection, dct->data, dct->workArea, dct->tables);

#ifdef ATFFT_TYPE_DOUBLE
    memcpy (out, dct->data, dct->dataSize);
#else
    atfft_double_to_sample_real (dct->data, out, dct->size);
#endif

    if (dct->direction == ATFFT_BACKWARD)
        atfft_scale_real (out, dct->size, 2.0);
}
