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
#include <assert.h>
#include <math.h>
#include <atfft/dct.h>
#include "ooura.h"

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

    /* buffer for in place transform */
    size_t n_data_bytes;
    double *data;

    /* buffers for ooura state */
    int *work_area;
    double *tables;
};

int atfft_dct_is_supported_size (int size)
{
    int min = 2;
    return atfft_is_power_of_2 (size) && size >= min;
}

struct atfft_dct* atfft_dct_create (int size, enum atfft_direction direction)
{
    /* ooura only supports sizes which are a power of 2. */
    assert (atfft_dct_is_supported_size (size));

    struct atfft_dct *plan;

    if (!(plan = calloc (1, sizeof (*plan))))
        return NULL;

    plan->size = size;
    plan->direction = direction;

    if (direction == ATFFT_FORWARD)
        plan->ooura_direction = -1;
    else
        plan->ooura_direction = 1;

    plan->n_data_bytes = size * sizeof (*(plan->data));
    plan->data = malloc (plan->n_data_bytes);

    int work_size = (2 + (1 << (int) (log (size / 2 + 0.5) / log (2)) / 2));
    plan->work_area = malloc (work_size * sizeof (*plan->work_area));
    plan->tables = malloc ((size * 5 / 4) * sizeof (*(plan->tables)));

    /* clean up on failure */
    if (!(plan->data && plan->work_area && plan->tables))
        goto failed;

    /* run a transform to initialise ooura state */
    plan->work_area [0] = 0;
    ddct (plan->size, plan->ooura_direction, plan->data, plan->work_area, plan->tables);

    return plan;

failed:
    atfft_dct_destroy (plan);
    return NULL;
}

void atfft_dct_destroy (struct atfft_dct *plan)
{
    if (plan)
    {
        free (plan->tables);
        free (plan->work_area);
        free (plan->data);
        free (plan);
    }
}

void atfft_dct_transform (struct atfft_dct *plan, const atfft_sample *in, atfft_sample *out)
{
#ifdef ATFFT_TYPE_DOUBLE
    memcpy (plan->data, in, plan->n_data_bytes);
#else
    atfft_sample_to_double_real (in, plan->data, plan->size);
#endif

    if (plan->direction == ATFFT_BACKWARD)
        plan->data [0] *= 0.5;

    ddct (plan->size, plan->ooura_direction, plan->data, plan->work_area, plan->tables);

#ifdef ATFFT_TYPE_DOUBLE
    memcpy (out, plan->data, plan->n_data_bytes);
#else
    atfft_double_to_sample_real (plan->data, out, plan->size);
#endif
}
