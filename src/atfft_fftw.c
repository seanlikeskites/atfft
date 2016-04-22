/*
 * Copyright (C) 2016 Sean Enderby <sean.enderby@gmail.com>
 *
 * This program is free software. It comes without any warranty, to
 * the extent permitted by applicable law. You can redistribute it 
 * and/or modify it under the terms of the Do What The Fuck You Want 
 * To Public License, Version 2, as published by Sam Hocevar. See 
 * the COPYING file for more details.
 */

#include <math.h>
#include <fftw3.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <atfft.h>

#if defined(ATFFT_TYPE_FLOAT)
#   define ATFFT_FFTW_MALLOC fftwf_malloc
#   define ATFFT_FFTW_FREE fftwf_free
#   define ATFFT_FFTW_DESTROY_PLAN fftwf_destroy_plan
#   define ATFFT_FFTW_PLAN_DFT_1D fftwf_plan_dft_1d
#   define ATFFT_FFTW_PLAN_DFT_R2C_1D fftwf_plan_dft_r2c_1d
#   define ATFFT_FFTW_PLAN_DFT_C2R_1D fftwf_plan_dft_c2r_1d
#   define ATFFT_FFTW_EXECUTE fftwf_execute
    typedef fftwf_plan atfft_fftw_plan;
    
#elif defined (ATFFT_TYPE_DOUBLE)
#   define ATFFT_FFTW_MALLOC fftw_malloc
#   define ATFFT_FFTW_FREE fftw_free
#   define ATFFT_FFTW_DESTROY_PLAN fftw_destroy_plan
#   define ATFFT_FFTW_PLAN_DFT_1D fftw_plan_dft_1d
#   define ATFFT_FFTW_PLAN_DFT_R2C_1D fftw_plan_dft_r2c_1d
#   define ATFFT_FFTW_PLAN_DFT_C2R_1D fftw_plan_dft_c2r_1d
#   define ATFFT_FFTW_EXECUTE fftw_execute
    typedef fftw_plan atfft_fftw_plan;

#elif defined(ATFFT_TYPE_LONG_DOUBLE)
#   define ATFFT_FFTW_MALLOC fftwl_malloc
#   define ATFFT_FFTW_FREE fftwl_free
#   define ATFFT_FFTW_DESTROY_PLAN fftwl_destroy_plan
#   define ATFFT_FFTW_PLAN_DFT_1D fftwl_plan_dft_1d
#   define ATFFT_FFTW_PLAN_DFT_R2C_1D fftwl_plan_dft_r2c_1d
#   define ATFFT_FFTW_PLAN_DFT_C2R_1D fftwl_plan_dft_c2r_1d
#   define ATFFT_FFTW_EXECUTE fftwl_execute
    typedef fftwl_plan atfft_fftw_plan;
#endif

struct atfft
{
    int size;
    int direction;
    enum atfft_format format;
    int inSize, outSize;
    atfft_sample *in, *out;
    atfft_fftw_plan plan;
};

struct atfft* atfft_create (int size, int direction, enum atfft_format format)
{
    struct atfft *fft;

    if (!(fft = malloc (sizeof (*fft))))
        return NULL;

    fft->size = size;
    fft->direction = direction;
    fft->format = format;

    switch (format)
    {
        case ATFFT_COMPLEX:
            fft->inSize = 2 * size * sizeof (*(fft->in));
            fft->in = ATFFT_FFTW_MALLOC (fft->inSize);

            fft->outSize = 2 * size * sizeof (*(fft->out));
            fft->out = ATFFT_FFTW_MALLOC (fft->outSize);

            fft->plan = ATFFT_FFTW_PLAN_DFT_1D (size,
                                                (atfft_complex*) fft->in,
                                                (atfft_complex*) fft->out, 
                                                direction, 
                                                FFTW_ESTIMATE);
            break;

        case ATFFT_REAL:
            if (direction == ATFFT_FORWARD)
            {
                fft->inSize = size * sizeof (*(fft->in));
                fft->in = ATFFT_FFTW_MALLOC (fft->inSize);

                fft->outSize = 2 * (floor (size / 2) + 1) * sizeof (*(fft->out));
                fft->out = ATFFT_FFTW_MALLOC (fft->outSize);

                fft->plan = ATFFT_FFTW_PLAN_DFT_R2C_1D (size, 
                                                        fft->in,
                                                        (atfft_complex*) fft->out,
                                                        FFTW_ESTIMATE);
            }
            else
            {
                fft->inSize = 2 * (floor (size / 2) + 1) * sizeof (*(fft->in));
                fft->in = ATFFT_FFTW_MALLOC (fft->inSize);

                fft->outSize = size * sizeof (*(fft->out));
                fft->out = ATFFT_FFTW_MALLOC (fft->outSize);

                fft->plan = ATFFT_FFTW_PLAN_DFT_C2R_1D (size,
                                                        (atfft_complex*) fft->in,
                                                        fft->out,
                                                        FFTW_ESTIMATE);
            }
            break;
    }

    /* clean up on failure */
    if (!(fft->in && fft->out && fft->plan))
    {
        atfft_destroy (fft);
        fft = NULL;
    }

    return fft;
}

void atfft_destroy (struct atfft *fft)
{
    if (fft)
    {
        ATFFT_FFTW_DESTROY_PLAN (fft->plan);
        ATFFT_FFTW_FREE (fft->out);
        ATFFT_FFTW_FREE (fft->in);
        free (fft);
    }
}

void atfft_fftw_apply_transform (struct atfft *fft, atfft_sample *in, atfft_sample *out)
{
    memcpy (fft->in, in, fft->inSize);
    ATFFT_FFTW_EXECUTE (fft->plan);
    memcpy (out, fft->out, fft->outSize);
}

void atfft_complex_transform (struct atfft *fft, atfft_complex *in, atfft_complex *out)
{
    /* Only to be used with complex FFTs. */
    assert (fft->format == ATFFT_COMPLEX);

    atfft_fftw_apply_transform (fft, (atfft_sample*) in, (atfft_sample*) out);
}

void atfft_real_forward_transform (struct atfft *fft, atfft_sample *in, atfft_complex *out)
{
    /* Only to be used for forward real FFTs. */
    assert ((fft->format == ATFFT_REAL) && (fft->direction == ATFFT_FORWARD));

    atfft_fftw_apply_transform (fft, in, (atfft_sample*) out);
}

void atfft_real_backward_transform (struct atfft *fft, atfft_complex *in, atfft_sample *out)
{
    /* Only to be used for backward real FFTs. */
    assert ((fft->format == ATFFT_REAL) && (fft->direction == ATFFT_BACKWARD));

    atfft_fftw_apply_transform (fft, (atfft_sample*) in, out);
}
