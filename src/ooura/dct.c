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
#include <string.h>
#include <math.h>
#include <assert.h>
#include <atfft/dct.h>
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
    int ooura_direction;
    int data_size;
    double *data;
    int *work_area;
    double *tables;
};

struct atfft_dct* atfft_dct_create (int size, enum atfft_direction direction)
{
    /* Ooura only supports sizes which are a power of 2. */
    assert (atfft_is_power_of_2 (size));

    struct atfft_dct *dct;

    if (!(dct = malloc (sizeof (*dct))))
        return NULL;

    dct->size = size;
    dct->direction = direction;

    if (direction == ATFFT_FORWARD)
        dct->ooura_direction = -1;
    else
        dct->ooura_direction = 1;

    dct->data_size = size * sizeof (*(dct->data));
    dct->data = malloc (dct->data_size);

    int work_size = (2 + (1 << (int) (log (size / 2 + 0.5) / log (2)) / 2)) * sizeof (*(dct->work_area));
    dct->work_area = malloc (work_size);

    dct->tables = malloc ((size * 5 / 4) * sizeof (*(dct->tables)));

    /* clean up on failure */
    if (!(dct->data && dct->work_area && dct->tables))
    {
        atfft_dct_destroy (dct);
        dct = NULL;
    }
    else
    {
        dct->work_area [0] = 0;
    }

    return dct;
}

void atfft_dct_destroy (struct atfft_dct *dct)
{
    if (dct)
    {
        free (dct->tables);
        free (dct->work_area);
        free (dct->data);
        free (dct);
    }
}

void atfft_dct_transform (struct atfft_dct *dct, const atfft_sample *in, atfft_sample *out)
{
#ifdef ATFFT_TYPE_DOUBLE
    memcpy (dct->data, in, dct->data_size);
#else
    atfft_sample_to_double_real (in, dct->data, dct->size);
#endif

    if (dct->direction == ATFFT_BACKWARD)
        dct->data [0] *= 0.5;

    ddct (dct->size, dct->ooura_direction, dct->data, dct->work_area, dct->tables);

#ifdef ATFFT_TYPE_DOUBLE
    memcpy (out, dct->data, dct->data_size);
#else
    atfft_double_to_sample_real (dct->data, out, dct->size);
#endif

    if (dct->direction == ATFFT_BACKWARD)
        atfft_scale_real (out, dct->size, 2.0);
}
