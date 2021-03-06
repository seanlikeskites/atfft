/*
 * Copyright (c) 2021 Sean Enderby <sean.enderby@gmail.com>
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

#include <stdlib.h>
#include <assert.h>
#include <Accelerate/Accelerate.h>
#include <atfft/dft.h>

#ifdef ATFFT_TYPE_LONG_DOUBLE
#   warning vDSP only supports double precision floating point, \
            higher precision values will be demoted to double for FFT calculations.
#endif

#ifdef ATFFT_TYPE_FLOAT
#   define ATFFT_VDSP_DFT_SETUP vDSP_DFT_Setup
#   define ATFFT_VDSP_DFT_ZOP_CREATE_SETUP vDSP_DFT_zop_CreateSetup
#   define ATFFT_VDSP_DFT_ZROP_CREATE_SETUP vDSP_DFT_zrop_CreateSetup
#   define ATFFT_VDSP_DFT_DESTROY_SETUP vDSP_DFT_DestroySetup
#   define ATFFT_VDSP_DFT_EXECUTE vDSP_DFT_Execute
    typedef float atfft_vdsp_sample;

#else
#   define ATFFT_VDSP_DFT_SETUP vDSP_DFT_SetupD
#   define ATFFT_VDSP_DFT_ZOP_CREATE_SETUP vDSP_DFT_zop_CreateSetupD
#   define ATFFT_VDSP_DFT_ZROP_CREATE_SETUP vDSP_DFT_zrop_CreateSetupD
#   define ATFFT_VDSP_DFT_DESTROY_SETUP vDSP_DFT_DestroySetupD
#   define ATFFT_VDSP_DFT_EXECUTE vDSP_DFT_ExecuteD
    typedef double atfft_vdsp_sample;
#endif

struct atfft_dft
{
    int size;
    enum atfft_direction direction;
    enum atfft_format format;
    atfft_vdsp_sample *inR, *inI, *outR, *outI;
    ATFFT_VDSP_DFT_SETUP setup;
};

int atfft_is_supported_length_vdsp_dft (unsigned int length, enum atfft_format format)
{
    int min = 8;

    if (atfft_is_power_of_2 (length))
        return 1;

    if (format == ATFFT_REAL)
        min = 16;

    if ((!(length % 3) && atfft_is_power_of_2 (length / 3) && (length / 3 >= min)) ||
        (!(length % 5) && atfft_is_power_of_2 (length / 5) && (length / 5 >= min)))
        return 1;        
    else
        return 0;
}

struct atfft_dft* atfft_dft_create (int size, enum atfft_direction direction, enum atfft_format format)
{
    struct atfft_dft *fft;
    vDSP_DFT_Direction vdspDirection;

    /* vDSP only supports certain lengths */
    assert (atfft_is_supported_length_vdsp_dft (size, format));

    if (!(fft = malloc (sizeof (*fft))))
        return NULL;

    fft->size = size;
    fft->direction = direction;
    fft->format = format;

    fft->inR = malloc (size * sizeof (*(fft->inR)));
    fft->inI = malloc (size * sizeof (*(fft->inI)));
    fft->outR = malloc (size * sizeof (*(fft->outR)));
    fft->outI = malloc (size * sizeof (*(fft->outI)));

    if (direction == ATFFT_FORWARD)
        vdspDirection = vDSP_DFT_FORWARD;
    else
        vdspDirection = vDSP_DFT_INVERSE;

    switch (format)
    {
        case ATFFT_COMPLEX:
            fft->setup = ATFFT_VDSP_DFT_ZOP_CREATE_SETUP (NULL, size, vdspDirection);
            break;

        case ATFFT_REAL:
            fft->setup = ATFFT_VDSP_DFT_ZROP_CREATE_SETUP (NULL, size, vdspDirection);
            break;
    }

    if (!(fft->inR && fft->inI && fft->outR && fft->outI && fft->setup))
    {
        atfft_dft_destroy (fft);
        fft = NULL;
    }

    return fft;
}

void atfft_dft_destroy (struct atfft_dft *fft)
{
    if (fft)
    {
        ATFFT_VDSP_DFT_DESTROY_SETUP (fft->setup);
        free (fft->outI);
        free (fft->outR);
        free (fft->inI);
        free (fft->inR);
        free (fft);
    }
}

void atfft_complex_fftw_to_vdsp (atfft_complex *in,
                                 atfft_vdsp_sample *real, atfft_vdsp_sample *imag, int size)
{
    int i = 0;

    for (i = 0; i < size; ++i)
    {
        real [i] = ATFFT_RE (in [i]);
        imag [i] = ATFFT_IM (in [i]);
    }
}

void atfft_complex_vdsp_to_fftw (const atfft_vdsp_sample *real, const atfft_vdsp_sample *imag,
                                 atfft_complex *out, int size)
{
    int i = 0;

    for (i = 0; i < size; ++i)
    {
        ATFFT_RE (out [i]) = real [i];
        ATFFT_IM (out [i]) = imag [i];
    }
}

void atfft_dft_complex_transform (struct atfft_dft *fft, atfft_complex *in, atfft_complex *out)
{
    /* Only to be used with complex FFTs. */
    assert (fft->format == ATFFT_COMPLEX);

    atfft_complex_fftw_to_vdsp (in, fft->inR, fft->inI, fft->size);
    ATFFT_VDSP_DFT_EXECUTE (fft->setup, fft->inR, fft->inI, fft->outR, fft->outI);
    atfft_complex_vdsp_to_fftw (fft->outR, fft->outI, out, fft->size);
}

void atfft_real_fftw_to_vdsp (const atfft_sample *in, 
                              atfft_vdsp_sample *outE, atfft_vdsp_sample *outO, int size)
{
    int i = 0;

    for (i = 0; i < size / 2; ++i)
    {
        outE [i] = in [2 * i];
        outO [i] = in [2 * i + 1];
    }
}

void atfft_halfcomplex_vdsp_to_fftw (const atfft_vdsp_sample *inE, const atfft_vdsp_sample *inO,
                                     atfft_complex *out, int size)
{
    int i = 0;
    int halfSize = size / 2;

    ATFFT_RE (out [0]) = inE [0] / 2.0;
    ATFFT_IM (out [0]) = 0;

    for (i = 1; i < halfSize; ++i)
    {
        ATFFT_RE (out [i]) = inE [i] / 2.0;
        ATFFT_IM (out [i]) = inO [i] / 2.0;
    }

    ATFFT_RE (out [halfSize]) = inO [0] / 2.0;
    ATFFT_IM (out [halfSize]) = 0;
}

void atfft_dft_real_forward_transform (struct atfft_dft *fft, const atfft_sample *in, atfft_complex *out)
{
    /* Only to be used for forward real FFTs. */
    assert ((fft->format == ATFFT_REAL) && (fft->direction == ATFFT_FORWARD));

    atfft_real_fftw_to_vdsp (in, fft->inR, fft->inI, fft->size);
    ATFFT_VDSP_DFT_EXECUTE (fft->setup, fft->inR, fft->inI, fft->outR, fft->outI);
    atfft_halfcomplex_vdsp_to_fftw (fft->outR, fft->outI, out, fft->size);
}

void atfft_halfcomplex_fftw_to_vdsp (atfft_complex *in, 
                                     atfft_vdsp_sample *outE, atfft_vdsp_sample *outO, int size)
{
    int i = 0;
    int halfSize = size / 2;

    outE [0] = ATFFT_RE (in [0]);

    for (i = 1; i < halfSize; ++i)
    {
        outE [i] = ATFFT_RE (in [i]);
        outO [i] = ATFFT_IM (in [i]);
    }

    outO [0] = ATFFT_RE (in [halfSize]);
}

void atfft_real_vdsp_to_fftw (const atfft_vdsp_sample *inE, const atfft_vdsp_sample *inO,
                              atfft_sample *out, int size)
{
    int i = 0;

    for (i = 0; i < size / 2; ++i)
    {
        out [2 * i] = inE [i];
        out [2 * i + 1] = inO [i];
    }
}

void atfft_dft_real_backward_transform (struct atfft_dft *fft, atfft_complex *in, atfft_sample *out)
{
    /* Only to be used for backward real FFTs. */
    assert ((fft->format == ATFFT_REAL) && (fft->direction == ATFFT_BACKWARD));

    atfft_halfcomplex_fftw_to_vdsp (in, fft->inR, fft->inI, fft->size);
    ATFFT_VDSP_DFT_EXECUTE (fft->setup, fft->inR, fft->inI, fft->outR, fft->outI);
    atfft_real_vdsp_to_fftw (fft->outR, fft->outI, out, fft->size);
}

