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
#include <math.h>
#include <limits.h>
#include <atfft/dft.h>
#include "dft_cooley_tukey.h"
#include "atfft_internal.h"

/* Because an int is used to represent the size of the transform
 * valid sizes are anywhere between 0 and 2^(n - 1), where n is
 * the number of bits in an int. The maximum number of radices 
 * will therefore be (n - 1) as it will occur when the radices
 * are all 2s.
 */
#define MAX_RADICES (sizeof (int) * CHAR_BIT - 1)

#ifndef ATFFT_SUB_TRANSFORM_THRESHOLD
#define ATFFT_SUB_TRANSFORM_THRESHOLD 4
#endif /* ATFFT_SUB_TRANSFORM_THRESHOLD */

struct atfft_dft_cooley_tukey
{
    int size;
    enum atfft_direction direction;
    enum atfft_format format;

    /* radices and their associated sub transform sizes */
    int n_radices;
    int radices [MAX_RADICES];
    int sub_sizes [MAX_RADICES];

    /* complex sinusoids */
    atfft_complex *sinusoids;

    /* twiddle factors */
    atfft_complex **t_factors;

    /* constants */
    atfft_sample sin_2pi_on_3;

    /* working space for length-n butterflies */
    atfft_complex *work_space;

    /* plans for large prime factor sub-transforms */
    int n_sub_transforms;
    struct atfft_dft **sub_transforms;
    struct atfft_dft *radix_sub_transforms [MAX_RADICES];
};

/******************************************
 * Functions for decomposing transform
 * size into nice factors.
 ******************************************/
static int atfft_next_radix (int r)
{
    switch (r)
    {
        case 4:
            return 2;
        case 2:
            return 3;
        default:
            return r + 2;
    }
}

static int atfft_init_radices (int size, int *radices, int *sub_sizes, int *max_r)
{
    /* current radix */
    int r = 4;
    int n_radices = 0;
    int sqrt_size = (int) sqrt ((double) size);

    *max_r = 2;

    /* Factor out specific even radices first,
     * then any other prime factors.
     */
    do
    {
        while (size % r)
        {
            r = atfft_next_radix (r);

            /* a number will only have one prime factor greater than its square root */
            if (r > sqrt_size)
                r = size;
        }

        size /= r;

        radices [n_radices] = r;
        sub_sizes [n_radices] = size;

        if (r > *max_r)
            *max_r = r;

        ++n_radices;
    }
    while (size > 1);

    return n_radices;
}

/******************************************
 * Functions for pre-calculating complex
 * sinusoids and twiddle factors.
 ******************************************/
static void atfft_init_complex_sinusoids (atfft_complex *sinusoids,
                                          int size,
                                          enum atfft_direction direction)
{
    for (int i = 0; i < size; ++i)
    {
        atfft_twiddle_factor (i, size, direction, sinusoids + i);
    }
}

static void atfft_free_twiddle_factors (atfft_complex **factors,
                                        int n_radices)
{
    if (factors)
    {
        for (int i = 0; i < n_radices; ++i)
        {
            free (factors [i]);
        }

        free (factors);
    }
}

static atfft_complex* atfft_generate_twiddle_factors (int radix,
                                                      int sub_size,
                                                      enum atfft_direction direction)
{
    int size = radix * sub_size;
    int n_factors = size - sub_size;
    atfft_complex *f = malloc (n_factors * sizeof (*f));

    if (!f)
        return NULL;

    int n = 0;

    for (int k = 0; k < sub_size; ++k)
    {
        for (int r = 1; r < radix; ++r)
        {
            atfft_twiddle_factor (k * r, size, direction, f + n);
            ++n;
        }
    }

    return f;
}

static atfft_complex** atfft_init_twiddle_factors (int *radices,
                                                   int *sub_sizes,
                                                   int n_radices,
                                                   enum atfft_direction direction)
{
    atfft_complex **factors = calloc (n_radices, sizeof (*factors));

    if (!factors)
        return NULL;

    for (int i = 0; i < n_radices; ++i)
    {
        atfft_complex *f = atfft_generate_twiddle_factors (radices [i],
                                                           sub_sizes [i],
                                                           direction);

        if (!f)
            goto failed;

        factors [i] = f;
    }

    return factors;

failed:
    atfft_free_twiddle_factors (factors, n_radices);
    return NULL;
}

/******************************************
 * atfft_dft_cooley_tukey struct management
 ******************************************/
struct atfft_dft_cooley_tukey* atfft_dft_cooley_tukey_create (int size,
                                                              enum atfft_direction direction,
                                                              enum atfft_format format)
{
    struct atfft_dft_cooley_tukey *fft;

    if (!(fft = calloc (1, sizeof (*fft))))
        return NULL;

    fft->size = size;
    fft->direction = direction;
    fft->format = format;

    /* calculate radices */
    int max_r = 0;
    fft->n_radices = atfft_init_radices (size, fft->radices, fft->sub_sizes, &max_r);

    /* calculate twiddle factors */
    fft->sinusoids = malloc (size * sizeof (*(fft->sinusoids)));
    fft->t_factors = atfft_init_twiddle_factors (fft->radices,
                                                 fft->sub_sizes,
                                                 fft->n_radices,
                                                 direction);

    /* clean up on failure */
    if (!(fft->sinusoids && fft->t_factors))
        goto failed;
    else
        atfft_init_complex_sinusoids (fft->sinusoids, size, direction);

    /* initialise constants */
    if (direction == ATFFT_FORWARD)
        fft->sin_2pi_on_3 = - sin (2.0 * M_PI / 3.0);
    else
        fft->sin_2pi_on_3 = sin (2.0 * M_PI / 3.0);

    /* allocate some working space */
    fft->work_space = malloc (max_r * sizeof (*(fft->work_space)));

    if (!fft->work_space)
        goto failed;

    /* create any necessary sub-transform strucs */
    fft->sub_transforms = atfft_init_sub_transforms (fft->radices,
                                                     MAX_RADICES,
                                                     &(fft->n_sub_transforms),
                                                     fft->radix_sub_transforms,
                                                     direction,
                                                     format,
                                                     ATFFT_SUB_TRANSFORM_THRESHOLD);

    if (fft->n_sub_transforms > 0 && !fft->sub_transforms)
        goto failed;

    return fft;

failed:
    atfft_dft_cooley_tukey_destroy (fft);
    return NULL;
}

void atfft_dft_cooley_tukey_destroy (struct atfft_dft_cooley_tukey *fft)
{
    if (fft)
    {
        atfft_free_sub_transforms (fft->sub_transforms, fft->n_sub_transforms);
        free (fft->work_space);
        atfft_free_twiddle_factors (fft->t_factors, fft->n_radices);
        free (fft->sinusoids);
        free (fft);
    }
}

/******************************************
 * Optimised small size DFTs
 ******************************************/
static void atfft_dft_2 (atfft_complex *out,
                         int stride)
{
    /* The start of each block in the output DFT */
    atfft_complex *bin1 = out;
    atfft_complex *bin2 = bin1 + stride;
    atfft_complex t;

    atfft_copy_complex (*bin2, &t);
    atfft_difference_complex (*bin1, t, bin2);
    atfft_sum_complex (*bin1, t, bin1);
}

static void atfft_dft_3 (atfft_complex *out,
                         int stride,
                         atfft_sample sin_2pi_on_3)
{
    atfft_complex *bins [3];
    bins [0] = out;
    bins [1] = bins [0] + stride;
    bins [2] = bins [1] + stride;

    atfft_complex ts [3];

    atfft_sum_complex (*bins [1], *bins [2], &ts [0]);
    ATFFT_RE (ts [1]) = ATFFT_RE (*bins [0]) - ATFFT_RE (ts [0]) / 2.0;
    ATFFT_IM (ts [1]) = ATFFT_IM (*bins [0]) - ATFFT_IM (ts [0]) / 2.0;
    atfft_difference_complex (*bins [1], *bins [2], &ts [2]);
    ATFFT_RE (ts [2]) = sin_2pi_on_3 * ATFFT_RE (ts [2]);
    ATFFT_IM (ts [2]) = sin_2pi_on_3 * ATFFT_IM (ts [2]);

    atfft_sum_complex (*bins [0], ts [0], bins [0]);
    ATFFT_RE (*bins [1]) = ATFFT_RE (ts [1]) - ATFFT_IM (ts [2]);
    ATFFT_IM (*bins [1]) = ATFFT_IM (ts [1]) + ATFFT_RE (ts [2]);
    ATFFT_RE (*bins [2]) = ATFFT_RE (ts [1]) + ATFFT_IM (ts [2]);
    ATFFT_IM (*bins [2]) = ATFFT_IM (ts [1]) - ATFFT_RE (ts [2]);

}

static void atfft_dft_4 (atfft_complex **bins,
                         enum atfft_direction direction)
{
    atfft_complex ts [4];

    atfft_sum_complex (*bins [0], *bins [2], &ts [0]);
    atfft_sum_complex (*bins [1], *bins [3], &ts [1]);
    atfft_difference_complex (*bins [0], *bins [2], &ts [2]);

    if (direction == ATFFT_FORWARD)
        atfft_difference_complex (*bins [1], *bins [3], &ts [3]);
    else
        atfft_difference_complex (*bins [3], *bins [1], &ts [3]);

    atfft_sum_complex (ts [0], ts[1], bins [0]);
    ATFFT_RE (*bins [1]) = ATFFT_RE (ts [2]) + ATFFT_IM (ts [3]);
    ATFFT_IM (*bins [1]) = ATFFT_IM (ts [2]) - ATFFT_RE (ts [3]);
    atfft_difference_complex (ts [0], ts[1], bins [2]);
    ATFFT_RE (*bins [3]) = ATFFT_RE (ts [2]) - ATFFT_IM (ts [3]);
    ATFFT_IM (*bins [3]) = ATFFT_IM (ts [2]) + ATFFT_RE (ts [3]);
}

/******************************************
 * DFT Butterflies
 ******************************************/
static void atfft_butterfly_2 (atfft_complex *out,
                               int stride,
                               int sub_size,
                               atfft_complex *t_factors)
{
    int i = sub_size;
    int t = 0;
    int dft_stride = sub_size * stride;

    while (i--)
    {
        atfft_multiply_by_complex (out + dft_stride, t_factors [t]);
        atfft_dft_2 (out, dft_stride);

        out += stride;
        ++t;
    }
}

static void atfft_butterfly_3 (atfft_complex *out,
                               int stride,
                               int sub_size,
                               atfft_complex *t_factors,
                               atfft_sample sin_2pi_on_3)
{
    int i = sub_size;
    int radix = 3;
    int t = 0;
    int dft_stride = sub_size * stride;

    while (i--)
    {
        int m = dft_stride;

        for (int n = 1; n < radix; ++n)
        {
            atfft_multiply_by_complex (out + m, t_factors [t]);
            m += dft_stride;
            ++t;
        }

        atfft_dft_3 (out, dft_stride, sin_2pi_on_3);

        out += stride;
    }
}

static void atfft_butterfly_4 (atfft_complex *out,
                               int stride,
                               int sub_size,
                               enum atfft_direction direction,
                               atfft_complex *t_factors)
{
    int i = sub_size;
    int radix = 4;
    int t = 0;
    int dft_stride = sub_size * stride;

    atfft_complex *bins [radix];

    bins [0] = out;
    bins [1] = bins [0] + dft_stride;
    bins [2] = bins [1] + dft_stride;
    bins [3] = bins [2] + dft_stride;

    while (i--)
    {
        for (int n = 1; n < radix; ++n)
        {
            atfft_multiply_by_complex (bins [n], t_factors [t]);
            ++t;
        }

        atfft_dft_4 (bins, direction);

        for (int n = 0; n < radix; ++n)
        {
            bins [n] += stride;
        }
    }
}

static void atfft_butterfly_sub_transform (atfft_complex *out,
                                           int stride,
                                           int sub_size,
                                           int radix,
                                           atfft_complex *t_factors,
                                           struct atfft_dft *sub_transform)
{
    int i = sub_size;
    int t = 0;
    int dft_stride = sub_size * stride;

    while (i--)
    {
        int m = dft_stride;

        for (int n = 1; n < radix; ++n)
        {
            atfft_multiply_by_complex (out + m, t_factors [t]);
            m += dft_stride;
            ++t;
        }

        atfft_dft_complex_transform_stride (sub_transform, out, dft_stride, out, dft_stride);

        out += stride;
    }
}

static void atfft_butterfly_n (atfft_complex *out,
                               int stride,
                               int sub_size,
                               int radix,
                               atfft_complex *sinusoids,
                               int n_sinusoids,
                               int sin_stride,
                               atfft_complex *work_space)
{
    /* Combine radix DFTs of size sub_size,
     * into one DTF of size (radix * sub_size).
     */

    /* Loop over the bins in each of the sub-transforms. */
    for (int i = 0; i < sub_size; ++i)
    {
        /* Copy ith bin from each sub-transform into work_space. */
        for (int n = 0; n < radix; ++n)
        {
            int b = (n * sub_size + i) * stride;
            atfft_copy_complex (out [b], work_space + n);
        }

        /* Calculate the output bins. */
        for (int n = 0; n < radix; ++n)
        {
            /* In this iteration we are calculating the ith bin in the nth
             * block of the output DFT.
             * 
             *  ---------------------- Output DFT ----------------------
             *  ________________________________________________________
             * |                |                |     |                |
             * | 0, 1, ..., i-1 | 0, 1, ..., i-1 | ... | 0, 1, ..., i-1 |
             * |________________|________________|_____|________________|
             *  0th Block        1st Block              (radix - 1)th Block
             *
             *  k is the index of the current bin we are calculating in the
             *  output DFT.
             */
            int k = n * sub_size + i;
            int b = k * stride;

            /* copy the ith bin of the first sub-transform to the current
             * output bin.
             */
            atfft_copy_complex (work_space [0], out + b);

            /* Sum in the ith bins from the remaining sub-transforms,
             * multiplied by their respective twiddle factor.
             * out[b] += work_space [r] * sinusoids [(k * r * sin_stride) % n_sinusoids]
             */
            for (int r = 1; r < radix; ++r)
            {
                atfft_sample *bin = out [b];
                atfft_sample *in = work_space [r];
                atfft_sample *t = sinusoids [(k * r * sin_stride) % n_sinusoids];

                ATFFT_RE (bin) += ATFFT_RE (in) * ATFFT_RE (t) -
                                    ATFFT_IM (in) * ATFFT_IM (t);
                ATFFT_IM (bin) += ATFFT_RE (in) * ATFFT_IM (t) +
                                    ATFFT_IM (in) * ATFFT_RE (t);
            }
        }
    }
}

static void atfft_butterfly (const struct atfft_dft_cooley_tukey *fft,
                             atfft_complex *out,
                             int stride,
                             int sub_size,
                             int radix,
                             struct atfft_dft *sub_transform,
                             atfft_complex *t_factors,
                             int sin_stride)
{
    switch (radix)
    {
        case 2:
            atfft_butterfly_2 (out, stride, sub_size, t_factors);
            break;
        case 3:
            atfft_butterfly_3 (out, stride, sub_size, t_factors, fft->sin_2pi_on_3);
            break;
        case 4:
            atfft_butterfly_4 (out, stride, sub_size, fft->direction, t_factors);
            break;
        default:
            if (sub_transform)
                atfft_butterfly_sub_transform (out, stride, sub_size, radix, t_factors, sub_transform);
            else
                atfft_butterfly_n (out, stride, sub_size, radix, fft->sinusoids, fft->size, sin_stride, fft->work_space);
    }
}

/******************************************
 * DFT implementation
 ******************************************/
static void atfft_compute_dft_complex (const struct atfft_dft_cooley_tukey *fft,
                                       atfft_complex *in,
                                       atfft_complex *out,
                                       int stage,
                                       int in_stride,
                                       int out_stride,
                                       int sin_stride)
{
    /* Get the radix, R, for this stage of the transform.
     * We will split the transform into R sub-transforms
     * of size sub_size.
     */
    int R = fft->radices [stage];
    int sub_size = fft->sub_sizes [stage];

    if (stage < fft->n_radices - 1)
    {
        /* If there is another radix in the list
         * we recursively apply the transform, once
         * for each sub-transform in this stage. 
         */
        for (int r = 0; r < R; ++r)
        {
            atfft_compute_dft_complex (fft, 
                                       in + r * in_stride,
                                       out + r * sub_size * out_stride,
                                       stage + 1,
                                       in_stride * R,
                                       out_stride,
                                       sin_stride * R);
        }
    }
    else
    {
        /* If there are no more radices, we have
         * reached the first stage of our transform.
         * Here we apply the decimation in time.
         */
        for (int i = 0; i < sub_size * R; ++i)
        {
            atfft_copy_complex (in [i * in_stride], out + i * out_stride);
        }
    }

    /* Apply butterfly for this stage of the transform. */
    atfft_butterfly (fft,
                     out,
                     out_stride,
                     sub_size,
                     R,
                     fft->radix_sub_transforms [stage],
                     fft->t_factors [stage],
                     sin_stride);
}

void atfft_dft_cooley_tukey_complex_transform (struct atfft_dft_cooley_tukey *fft,
                                               atfft_complex *in,
                                               int in_stride,
                                               atfft_complex *out,
                                               int out_stride)
{
    /* Only to be used with complex FFTs. */
    assert (fft->format == ATFFT_COMPLEX);

    atfft_compute_dft_complex (fft,
                               in,
                               out,
                               0,
                               in_stride,
                               out_stride,
                               1);
}
