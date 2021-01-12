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
#include <atfft/dct.h>

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
}
