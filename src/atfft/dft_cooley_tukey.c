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
#include <atfft/dft.h>
#include "dft_cooley_tukey.h"
#include "atfft_internal.h"
#include "constants.h"
#include "dft_plan.h"

#ifndef ATFFT_SUB_TRANSFORM_THRESHOLD
#define ATFFT_SUB_TRANSFORM_THRESHOLD 8
#endif /* ATFFT_SUB_TRANSFORM_THRESHOLD */

static const atfft_sample one_on_root_two = 0.7071067811865475244008443621048;

struct atfft_dft_cooley_tukey
{
    enum atfft_dft_algorithm algorithm;

    int size;
    enum atfft_direction direction;
    enum atfft_format format;

    /* radices and their associated sub transform sizes */
    int n_radices;
    int radices [MAX_INT_FACTORS];
    int sub_sizes [MAX_INT_FACTORS];

    /* index permutation for decimation in time */
    int *permutation;

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
    struct atfft_dft *radix_sub_transforms [MAX_INT_FACTORS];
};

/******************************************
 * Functions for decomposing transform
 * size into nice factors.
 ******************************************/
static int atfft_next_radix (int r)
{
    switch (r)
    {
        case 8:
            return 4;
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
    int r = 8;
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
 * Functions for calculating decimation in
 * time permutation.
 ******************************************/
static void atfft_permute_index (int *permutation,
                                 int stage,
                                 int start,
                                 int step,
                                 const int *radices,
                                 const int *sub_sizes,
                                 int n_radices)
{
    /* Get the radix, R, for this stage of the transform.
     * We will split the transform into R sub-transforms
     * of size sub_size.
     */
    int R = radices [stage];
    int sub_size = sub_sizes [stage];

    if (stage < n_radices - 1)
    {
        /* If there is another radix in the list
         * we recursively apply the permutation, once
         * for each sub-transform in this stage. 
         */
        for (int r = 0; r < R; ++r)
        {
            atfft_permute_index (permutation + r * sub_size,
                                 stage + 1,
                                 start + r * step,
                                 step * R,
                                 radices,
                                 sub_sizes,
                                 n_radices);
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
            permutation [i] = start + i * step;
        }
    }
}

int* atfft_init_index_permutation (int size,
                                   const int *radices,
                                   const int *sub_sizes,
                                   int n_radices)
{
    int *permutation = malloc (size * sizeof (*permutation));
    int *test = malloc (size * sizeof (*test));

    if (!permutation)
        goto failed;

    atfft_permute_index (permutation, 0, 0, 1, radices, sub_sizes, n_radices);

failed:
    return permutation;
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

static atfft_complex** atfft_init_twiddle_factors (const int *radices,
                                                   const int *sub_sizes,
                                                   int n_radices,
                                                   enum atfft_direction direction)
{
    atfft_complex **factors = calloc (n_radices, sizeof (*factors));

    if (!factors)
        return NULL;

    /* The deepest stage has a sub-size of 1 so all twiddle factors
     * will be 1. A null pointer in the twiddle factors array
     * indicates no multiplication need be applied. */
    for (int i = 0; i < n_radices - 1; ++i)
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

    fft->algorithm = ATFFT_COOLEY_TUKEY;
    fft->size = size;
    fft->direction = direction;
    fft->format = format;

    /* calculate radices */
    int max_r = 0;
    fft->n_radices = atfft_init_radices (size, fft->radices, fft->sub_sizes, &max_r);

    /* calculate permutation */
    fft->permutation = atfft_init_index_permutation (size,
                                                     fft->radices,
                                                     fft->sub_sizes,
                                                     fft->n_radices);

    if (!fft->permutation)
        goto failed;

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
                                                     MAX_INT_FACTORS,
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

void atfft_dft_cooley_tukey_destroy (void *fft)
{
    struct atfft_dft_cooley_tukey *t = fft;

    if (t)
    {
        atfft_free_sub_transforms (t->sub_transforms, t->n_sub_transforms);
        free (t->work_space);
        atfft_free_twiddle_factors (t->t_factors, t->n_radices);
        free (t->sinusoids);
        free (t);
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

static void atfft_dft_8 (atfft_complex **bins,
                         enum atfft_direction direction)
{
    atfft_complex ts [8];
    atfft_complex qs [6];
    atfft_complex ss [4];

    for (int i = 0; i < 4; ++i)
    {
        atfft_sum_complex (*bins [i], *bins [i + 4], &ts [i]);
    }

    atfft_difference_complex (*bins [0], *bins [4], &ts [4]);
    atfft_difference_complex (*bins [1], *bins [5], &ts [5]);
    atfft_difference_complex (*bins [3], *bins [7], &ts [7]);

    for (int i = 0; i < 2; ++i)
    {
        atfft_sum_complex (ts [i], ts [i + 2], &qs [i]);
    }

    atfft_difference_complex (ts [0], ts [2], &qs [2]);
    ATFFT_RE (qs [5]) = one_on_root_two * (ATFFT_RE (ts [5]) - ATFFT_RE (ts [7]));
    ATFFT_IM (qs [5]) = one_on_root_two * (ATFFT_IM (ts [5]) - ATFFT_IM (ts [7]));

    if (direction == ATFFT_FORWARD)
    {
        atfft_difference_complex (*bins [2], *bins [6], &ts [6]);

        atfft_difference_complex (ts [1], ts [3], &qs [3]);

        ATFFT_RE (qs [4]) = one_on_root_two * (ATFFT_RE (ts [5]) + ATFFT_RE (ts [7]));
        ATFFT_IM (qs [4]) = one_on_root_two * (ATFFT_IM (ts [5]) + ATFFT_IM (ts [7]));
    }
    else
    {
        atfft_difference_complex (*bins [6], *bins [2], &ts [6]);

        atfft_difference_complex (ts [3], ts [1], &qs [3]);

        ATFFT_RE (qs [4]) = - one_on_root_two * (ATFFT_RE (ts [5]) + ATFFT_RE (ts [7]));
        ATFFT_IM (qs [4]) = - one_on_root_two * (ATFFT_IM (ts [5]) + ATFFT_IM (ts [7]));
    }

    atfft_difference_a_jb_complex (ts [4], qs [4], &ss [0]);
    atfft_sum_a_jb_complex (ts [4], qs [4], &ss [1]);
    atfft_difference_a_jb_complex (qs [5], ts [6], &ss [2]);
    atfft_sum_a_jb_complex (qs [5], ts [6], &ss [3]);

    atfft_sum_complex (qs [0], qs [1], bins [0]);
    atfft_sum_complex (ss [0], ss [2], bins [1]);
    atfft_difference_a_jb_complex (qs [2], qs [3], bins [2]);
    atfft_difference_complex (ss [0], ss [2], bins [3]);
    atfft_difference_complex (qs [0], qs [1], bins [4]);
    atfft_difference_complex (ss [1], ss [3], bins [5]);
    atfft_sum_a_jb_complex (qs [2], qs [3], bins [6]);
    atfft_sum_complex (ss [1], ss [3], bins [7]);
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
        if (t_factors)
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
        if (t_factors)
        {
            int m = dft_stride;

            for (int n = 1; n < radix; ++n)
            {
                atfft_multiply_by_complex (out + m, t_factors [t]);
                m += dft_stride;
                ++t;
            }
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

    atfft_complex *bins [4];

    bins [0] = out;
    bins [1] = bins [0] + dft_stride;
    bins [2] = bins [1] + dft_stride;
    bins [3] = bins [2] + dft_stride;

    while (i--)
    {
        if (t_factors)
        {
            for (int n = 1; n < radix; ++n)
            {
                atfft_multiply_by_complex (bins [n], t_factors [t]);
                ++t;
            }
        }

        atfft_dft_4 (bins, direction);

        for (int n = 0; n < radix; ++n)
        {
            bins [n] += stride;
        }
    }
}

static void atfft_butterfly_8 (atfft_complex *out,
                               int stride,
                               int sub_size,
                               enum atfft_direction direction,
                               atfft_complex *t_factors)
{
    int i = sub_size;
    const int radix = 8;
    int t = 0;
    int dft_stride = sub_size * stride;

    atfft_complex *bins [radix];

    bins [0] = out;

    for (int n = 1; n < radix; ++n)
    {
        bins [n] = bins [n - 1] + dft_stride;
    }

    while (i--)
    {
        if (t_factors)
        {
            for (int n = 1; n < radix; ++n)
            {
                atfft_multiply_by_complex (bins [n], t_factors [t]);
                ++t;
            }
        }

        atfft_dft_8 (bins, direction);

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
        if (t_factors)
        {
            int m = dft_stride;

            for (int n = 1; n < radix; ++n)
            {
                atfft_multiply_by_complex (out + m, t_factors [t]);
                m += dft_stride;
                ++t;
            }
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
        case 8:
            atfft_butterfly_8 (out, stride, sub_size, fft->direction, t_factors);
            break;
        default:
            if (sub_transform)
                atfft_butterfly_sub_transform (out, stride, sub_size, radix, t_factors, sub_transform);
            else
                atfft_butterfly_n (out, stride, sub_size, radix, fft->sinusoids, fft->size, sin_stride, fft->work_space);
    }
}

/******************************************
 * Recursive DFT implementation.
 ******************************************/
static void atfft_compute_dft_recursive (const struct atfft_dft_cooley_tukey *fft,
                                         atfft_complex *in,
                                         int in_stride,
                                         atfft_complex *out,
                                         int out_stride,
                                         int stage,
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
            atfft_compute_dft_recursive (fft, 
                                         in + r * in_stride,
                                         in_stride * R,
                                         out + r * sub_size * out_stride,
                                         out_stride,
                                         stage + 1,
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

/******************************************
 * Iterative DFT implementation.
 ******************************************/
static void atfft_decimate_in_time (const struct atfft_dft_cooley_tukey *fft,
                                    atfft_complex *in,
                                    int in_stride,
                                    atfft_complex *out,
                                    int out_stride)
{
    for (int i = 0; i < fft->size; ++i)
    {
        atfft_copy_complex (in [fft->permutation [i] * in_stride], out + i * out_stride);
    }
}

static void atfft_apply_dft_stage (const struct atfft_dft_cooley_tukey *fft,
                                   atfft_complex *out,
                                   int stride,
                                   int radix,
                                   int sub_size,
                                   atfft_complex *t_factors,
                                   struct atfft_dft *sub_transform,
                                   int sin_stride)
{
    for (int i = 0; i < fft->size; i += radix * sub_size)
    {
        /* Apply butterfly for this stage of the transform. */
        atfft_butterfly (fft,
                         out + i,
                         stride,
                         sub_size,
                         radix,
                         sub_transform,
                         t_factors,
                         sin_stride);
    }
}

static void atfft_compute_dft_iterative (const struct atfft_dft_cooley_tukey *fft,
                                         atfft_complex *in,
                                         int in_stride,
                                         atfft_complex *out,
                                         int out_stride)
{
    /* Apply decimation in time. */
    atfft_decimate_in_time (fft, in, in_stride, out, out_stride);

    int stage = fft->n_radices;
    int sine_stride = fft->size;

    while (stage--)
    {
        int R = fft->radices [stage];
        int sub_size = fft->sub_sizes [stage];
        atfft_complex *t_factors = fft->t_factors [stage];
        struct atfft_dft *sub_transform = fft->radix_sub_transforms [stage];

        sine_stride = sine_stride / R;

        atfft_apply_dft_stage (fft,
                               out,
                               out_stride,
                               R,
                               sub_size,
                               t_factors,
                               sub_transform,
                               sine_stride);
    }
}

/******************************************
 * Apply Transform
 ******************************************/
void atfft_dft_cooley_tukey_complex_transform (void *fft,
                                               atfft_complex *in,
                                               int in_stride,
                                               atfft_complex *out,
                                               int out_stride)
{
    struct atfft_dft_cooley_tukey *t = fft;

    /* Only to be used with complex FFTs. */
    assert (t->format == ATFFT_COMPLEX);

    if (0)
        atfft_compute_dft_recursive (t,
                                     in,
                                     in_stride,
                                     out,
                                     out_stride,
                                     0,
                                     1);
    else
        atfft_compute_dft_iterative (t,
                                     in,
                                     in_stride,
                                     out,
                                     out_stride);
}

/******************************************
 * Get plan info.
 ******************************************/
static cJSON* atfft_get_plan_stage (struct atfft_dft_cooley_tukey *fft, int stage_idx)
{
    cJSON *radix = NULL,
          *sub_size = NULL;

    cJSON *stage = cJSON_CreateObject();

    if (!stage)
        goto failed;

    radix = cJSON_AddNumberToObject (stage, "Radix", fft->radices [stage_idx]);
    sub_size = cJSON_AddNumberToObject (stage, "Sub-Size", fft->sub_sizes [stage_idx]);

    if (!(radix && sub_size))
        goto failed;

    if (fft->radix_sub_transforms [stage_idx])
    {
        cJSON *sub_plan = atfft_dft_get_plan (fft->radix_sub_transforms [stage_idx]);

        if (!sub_plan)
            goto failed;

        cJSON_AddItemToObject (stage, "Sub-Transform", sub_plan);
    }

    return stage;

failed:
    cJSON_Delete (stage);
    return NULL;
}

static cJSON* atfft_get_plan_stages (struct atfft_dft_cooley_tukey *fft)
{
    cJSON *stages = cJSON_CreateArray();

    if (!stages)
        goto failed;

    for (int i = 0; i < fft->n_radices; ++i)
    {
        cJSON* stage = atfft_get_plan_stage (fft, i);

        if (!stage)
            goto failed;

        cJSON_AddItemToArray (stages, stage);
    }

    return stages;

failed:
    cJSON_Delete (stages);
    return NULL;
}

cJSON* atfft_dft_cooley_tukey_get_plan (struct atfft_dft_cooley_tukey *fft)
{
    cJSON *alg = NULL,
          *size = NULL,
          *stages = NULL;

    cJSON *plan_structure = cJSON_CreateObject();

    if (!plan_structure)
        goto failed;

    alg = cJSON_AddStringToObject (plan_structure, "Algorithm", "Cooley-Tukey");
    size = cJSON_AddNumberToObject (plan_structure, "Size", fft->size);

    if (!(alg && size))
        goto failed;

    stages = atfft_get_plan_stages (fft);

    if (!stages)
        goto failed;

    cJSON_AddItemToObject (plan_structure, "Stages", stages);

    return plan_structure;

failed:
    cJSON_Delete (plan_structure);
    return NULL;
}
