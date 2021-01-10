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
#include <atfft/dft.h>
#include <atfft/dft_util.h>
#include <ffts/ffts.h>

#ifndef ATFFT_TYPE_FLOAT
#   ifdef _MSC_VER
#       pragma message(": warning: FFTS only supports single precision floating point, " \
                       "higher precision values will be demoted to float for FFT calculations.")
#   else
#       warning FFTS only supports single precision floating point, \
                higher precision values will be demoted to float for FFT calculations.
#   endif
#endif

struct atfft_dft
{
    int size;
    enum atfft_direction direction;
    enum atfft_format format;

    /* the ffts plan */
    ffts_plan_t *plan;

#ifndef ATFFT_TYPE_FLOAT
    /* input and output buffers for ffts transform */
    size_t in_size, out_size;
    float *in, *out;
#endif
};

int atfft_dft_is_supported_size (int size, enum atfft_format format)
{
    if (format == ATFFT_COMPLEX)
        return size > 0;
    else
        return size > 0 && atfft_is_even(size);
}

static ffts_plan_t* init_ffts_plan (int size,
                                    enum atfft_direction direction,
                                    enum atfft_format format)
{
    int ffts_direction = 0;

    if (direction == ATFFT_FORWARD)
        ffts_direction = -1;
    else
        ffts_direction = 1;

    if (format == ATFFT_COMPLEX)
        return ffts_init_1d (size, ffts_direction);
    else
        return ffts_init_1d_real (size, ffts_direction);
}

struct atfft_dft* atfft_dft_create (int size, enum atfft_direction direction, enum atfft_format format)
{
    /* fftw supports all sizes for complex transforms, but real transforms must be even. */
    assert (atfft_dft_is_supported_size (size, format));

    struct atfft_dft *plan;

    if (!(plan = calloc (1, sizeof (*plan))))
        return NULL;

    plan->size = size;
    plan->direction = direction;
    plan->format = format;

#ifndef ATFFT_TYPE_FLOAT
    /* allocate input and output buffers */
    if (format == ATFFT_COMPLEX)
    {
        plan->in_size = 2 * size;
        plan->out_size = 2 * size;
    }
    else
    {
        if (direction == ATFFT_FORWARD)
        {
            plan->in_size = size;
            plan->out_size = 2 * atfft_halfcomplex_size (size);
        }
        else
        {
            plan->in_size = 2 * atfft_halfcomplex_size (size);
            plan->out_size = size;
        }
    }

    plan->in = malloc (plan->in_size * sizeof (*(plan->in)));
    plan->out = malloc (plan->out_size * sizeof (*(plan->out)));

    if (!(plan->in && plan->out))
        goto failed;
#endif

    /* initialise the ffts plan */
    plan->plan = init_ffts_plan (size, direction, format);

    if (!plan->plan)
        goto failed;

    return plan;

failed:
    atfft_dft_destroy (plan);
    return NULL;
}

void atfft_dft_destroy (struct atfft_dft *plan)
{
    if (plan)
    {
        if (plan->plan) /* ffts can't hack freeing a null pointer */
            ffts_free (plan->plan);

#ifndef ATFFT_TYPE_FLOAT
        free (plan->out);
        free (plan->in);
#endif

        free (plan);
    }
}

static void apply_transform (struct atfft_dft *plan, const atfft_sample *in, atfft_sample *out)
{
#ifdef ATFFT_TYPE_FLOAT
    ffts_execute (plan->plan, (const float*) in, (float*) out);
#else
    atfft_sample_to_float_real (in, plan->in, plan->in_size);
    ffts_execute (plan->plan, plan->in, plan->out);
    atfft_float_to_sample_real (plan->out, out, plan->out_size);
#endif
}

void atfft_dft_complex_transform (struct atfft_dft *plan, atfft_complex *in, atfft_complex *out)
{
    /* Only to be used with complex FFTs. */
    assert (plan->format == ATFFT_COMPLEX);

    apply_transform (plan, (atfft_sample*) in, (atfft_sample*) out);
}

void atfft_dft_real_forward_transform (struct atfft_dft *plan, const atfft_sample *in, atfft_complex *out)
{
    /* Only to be used for forward real FFTs. */
    assert ((plan->format == ATFFT_REAL) && (plan->direction == ATFFT_FORWARD));

    apply_transform (plan, in, (atfft_sample*) out);
}

void atfft_dft_real_backward_transform (struct atfft_dft *plan, atfft_complex *in, atfft_sample *out)
{
    /* Only to be used for backward real FFTs. */
    assert ((plan->format == ATFFT_REAL) && (plan->direction == ATFFT_BACKWARD));

    apply_transform (plan, (atfft_sample*) in, out);
}
