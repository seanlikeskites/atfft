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
#include <atfft.h>

struct atfft
{
    int size;
    int direction;
    enum atfft_format format;
    double *inR, *inI, *outR, *outI;
    vDSP_DFT_SetupD setup;
};

int atfft_is_supported_length_vdsp (unsigned int length, enum atfft_format format)
{
    int min = 8;

    if (atfft_is_power_of_2 (length))
    {
        return 1;
    }

    if (format == ATFFT_REAL)
    {
        min = 16;
    }

    if ((!(length % 3) && atfft_is_power_of_2 (length / 3) && (length / 3 >= min)) ||
        (!(length % 5) && atfft_is_power_of_2 (length / 5) && (length / 5 >= min)))
    {
        return 1;        
    }
    else
    {
        return 0;
    }
}

struct atfft* atfft_create (int size, int direction, enum atfft_format format)
{
    struct atfft *fft;

    /* vDSP only supports certain lengths */
    assert (atfft_is_supported_length_vdsp (size, format));

    fft = malloc (sizeof (*fft));
    fft->size = size;
    fft->direction = direction;
    fft->format = format;

    fft->inR = malloc (size * sizeof (*(fft->inR)));
    fft->inI = malloc (size * sizeof (*(fft->inI)));
    fft->outR = malloc (size * sizeof (*(fft->outR)));
    fft->outI = malloc (size * sizeof (*(fft->outI)));

    switch (format)
    {
        case ATFFT_COMPLEX:
            fft->setup = vDSP_DFT_zop_CreateSetupD (NULL, size, -direction);
            break;

        case ATFFT_REAL:
            fft->setup = vDSP_DFT_zrop_CreateSetupD (NULL, size, -direction);
            break;
    }

    return fft;
}

void atfft_free (struct atfft *fft)
{
    vDSP_DFT_DestroySetupD (fft->setup);
    free (fft->outI);
    free (fft->outR);
    free (fft->inI);
    free (fft->inR);
    free (fft);
}

void atfft_complex_vdsp_to_fftw (double *real, double *imag, atfft_complex_double *out, int size)
{
    int i = 0;

    for (i = 0; i < size; ++i)
    {
        ATFFT_REAL (out [i]) = real [i];
        ATFFT_IMAG (out [i]) = imag [i];
    }
}

void atfft_complex_transform (struct atfft *fft, atfft_complex_double *in, atfft_complex_double *out)
{
    /* Only to be used with complex FFTs. */
    assert (fft->format == ATFFT_COMPLEX);

    atfft_real (in, fft->inR, fft->size);
    atfft_imag (in, fft->inI, fft->size);
    vDSP_DFT_ExecuteD (fft->setup, fft->inR, fft->inI, fft->outR, fft->outI);
    atfft_complex_vdsp_to_fftw (fft->outR, fft->outI, out, fft->size);
}

void atfft_real_fftw_to_vdsp (double *in, double *outE, double *outO, int size)
{
    int i = 0;

    for (i = 0; i < size / 2; ++i)
    {
        outE [i] = in [2 * i];
        outO [i] = in [2 * i + 1];
    }
}

void atfft_halfcomplex_vdsp_to_fftw (double *inE, double *inO, atfft_complex_double *out, int size)
{
    int i = 0;
    int halfSize = size / 2;

    ATFFT_REAL (out [0]) = inE [0] / 2.0;
    ATFFT_IMAG (out [0]) = 0;

    for (i = 1; i < halfSize; ++i)
    {
        ATFFT_REAL (out [i]) = inE [i] / 2.0;
        ATFFT_IMAG (out [i]) = inO [i] / 2.0;
    }

    ATFFT_REAL (out [halfSize]) = inO [0] / 2.0;
    ATFFT_IMAG (out [halfSize]) = 0;
}

void atfft_real_forward_transform (struct atfft *fft, double *in, atfft_complex_double *out)
{
    /* Only to be used for forward real FFTs. */
    assert ((fft->format == ATFFT_REAL) && (fft->direction == ATFFT_FORWARD));

    atfft_real_fftw_to_vdsp (in, fft->inR, fft->inI, fft->size);
    vDSP_DFT_ExecuteD (fft->setup, fft->inR, fft->inI, fft->outR, fft->outI);
    atfft_halfcomplex_vdsp_to_fftw (fft->outR, fft->outI, out, fft->size);
}

void atfft_halfcomplex_fftw_to_vdsp (atfft_complex_double *in, double *outE, double *outO, int size)
{
    int i = 0;
    int halfSize = size / 2;

    outE [0] = ATFFT_REAL (in [0]);

    for (i = 1; i < halfSize; ++i)
    {
        outE [i] = ATFFT_REAL (in [i]);
        outO [i] = ATFFT_IMAG (in [i]);
    }

    outO [0] = ATFFT_REAL (in [halfSize]);
}

void atfft_real_vdsp_to_fftw (double *inE, double *inO, double *out, int size)
{
    int i = 0;

    for (i = 0; i < size / 2; ++i)
    {
        out [2 * i] = inE [i];
        out [2 * i + 1] = inO [i];
    }
}

void atfft_real_backward_transform (struct atfft *fft, atfft_complex_double *in, double *out)
{
    /* Only to be used for backward real FFTs. */
    assert ((fft->format == ATFFT_REAL) && (fft->direction == ATFFT_BACKWARD));

    atfft_halfcomplex_fftw_to_vdsp (in, fft->inR, fft->inI, fft->size);
    vDSP_DFT_ExecuteD (fft->setup, fft->inR, fft->inI, fft->outR, fft->outI);
    atfft_real_vdsp_to_fftw (fft->outR, fft->outI, out, fft->size);
}

