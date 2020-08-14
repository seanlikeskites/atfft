/*
 * Copyright (C) 2020 Sean Enderby <sean.enderby@gmail.com>
 *
 * This program is free software. It comes without any warranty, to
 * the extent permitted by applicable law. You can redistribute it 
 * and/or modify it under the terms of the Do What The Fuck You Want 
 * To Public License, Version 2, as published by Sam Hocevar. See 
 * the COPYING file for more details.
 */

#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <limits.h>
#include <atfft/atfft_dft.h>
#include "../atfft_internal.h"
#include "atfft_dft_rader.h"

/* Because an int is used to represent the size of the transform
 * valid sizes are anywhere between 0 and 2^(n - 1), where n is
 * the number of bits in an int. The maximum number of radices 
 * will therefore be (n - 1) as it will occur when the radices
 * are all 2s.
 */
#define MAX_RADICES (sizeof (int) * CHAR_BIT - 1)

#ifndef ATFFT_RADER_THRESHOLD
#define ATFFT_RADER_THRESHOLD 4
#endif /* ATFFT_RADER_THRESHOLD */

struct atfft_dft
{
    int size;
    enum atfft_direction direction;
    enum atfft_format format;

    /* null terminated array of radices */
    int radices [MAX_RADICES + 1];

    /* twiddle factors */
    atfft_complex *t_factors;
    atfft_complex *stage_t_factors [MAX_RADICES];

    /* working space for butterflies */
    atfft_complex *work_space;

    /* sub fft objects for Rader's algorithm */
    int n_raders;
    struct atfft_dft_rader **raders;
    struct atfft_dft_rader *radix_raders [MAX_RADICES];
};

static void atfft_init_twiddle_factors (atfft_complex *factors,
                                        int size,
                                        enum atfft_direction direction)
{
    atfft_sample sin_factor = -1.0;

    if (direction == ATFFT_BACKWARD)
        sin_factor = 1.0;

    for (int i = 0; i < size; ++i)
    {
        atfft_sample x = 2.0 * i * M_PI / size;
        ATFFT_REAL (factors [i]) = cos (x);
        ATFFT_IMAG (factors [i]) = sin_factor * sin (x);
    }
}

static void atfft_free_twiddle_factors (atfft_complex **factors)
{
    if (factors)
    {
        for (int i = 0; i < MAX_RADICES; ++i)
        {
            free (factors [i]);
        }
    }
}

static void atfft_init_stage_twiddle_factors (atfft_complex **factors,
                                              int *radices,
                                              int *sub_sizes,
                                              int n_radices,
                                              enum atfft_direction direction)
{
    atfft_sample sin_factor = -1.0;

    if (direction == ATFFT_BACKWARD)
        sin_factor = 1.0;

    for (int i = 0; i < n_radices; ++i)
    {
        int radix = radices [i];
        int sub_size = sub_sizes [i];
        int size = radix * sub_size;
        int n_factors = size - sub_size;
        atfft_complex *f = malloc (n_factors * sizeof (*f));
        int n = 0;

        for (int k = 0; k < sub_size; ++k)
        {
            for (int r = 1; r < radix; ++r)
            {
                atfft_sample x = 2.0 * k * r * M_PI / size;
                ATFFT_REAL (f [n]) = cos (x);
                ATFFT_IMAG (f [n]) = sin_factor * sin (x);
                ++n;
            }
        }

        factors [i] = f;
    }
}

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

static int atfft_integer_is_in_array (const int *arr, int size, int member)
{
    for (int i = 0; i < size; ++i)
    {
        if (arr [i] == member)
            return 1;
    }

    return 0;
}

static int atfft_rader_radices (const int *radices, int *rader_radices, int size)
{
    int n_raders = 0;

    for (int i = 0; i < size; ++i)
    {
        if (radices [i] > ATFFT_RADER_THRESHOLD && 
            !atfft_integer_is_in_array (rader_radices, n_raders, radices [i]))
            rader_radices [n_raders++] = radices [i];
    }

    return n_raders;
}

static struct atfft_dft_rader* atfft_populate_rader (int radix,
                                                     const int *radices,
                                                     struct atfft_dft_rader **radix_raders,
                                                     enum atfft_direction direction,
                                                     enum atfft_format format)
{
    struct atfft_dft_rader *rader = atfft_dft_rader_create (radix, direction, format);

    if (!rader)
        return NULL;

    for (int i = 0; i < MAX_RADICES; ++i)
    {
        if (radices [i] == radix)
            radix_raders [i] = rader;
    }

    return rader;
}

static void atfft_free_raders (struct atfft_dft_rader **raders,
                               int size)
{
    if (raders)
    {
        for (int i = 0; i < size; ++i)
        {
            atfft_dft_rader_destroy (raders [i]);
        }

        free (raders);
    }
}

static struct atfft_dft_rader** atfft_init_raders (const int *radices,
                                                   int *n_raders,
                                                   struct atfft_dft_rader **radix_raders,
                                                   enum atfft_direction direction,
                                                   enum atfft_format format)
{
    int rader_radices [MAX_RADICES];

    /* get unique radices for which Rader's algorithm is required */
    *n_raders = atfft_rader_radices (radices, rader_radices, MAX_RADICES);

    if (*n_raders == 0)
        return NULL;

    /* allocate some space for them */
    struct atfft_dft_rader **raders = calloc (*n_raders, sizeof (*raders));

    if (!raders)
        return NULL;

    /* create the rader dft structs */
    for (int i = 0; i < *n_raders; ++i)
    {
        struct atfft_dft_rader *rader = atfft_populate_rader (rader_radices [i],
                                                              radices,
                                                              radix_raders,
                                                              direction,
                                                              format);

        if (!rader)
            goto failed;

        raders [i] = rader;
    }

    return raders;

failed:
    atfft_free_raders (raders, *n_raders);
    return NULL;
}

struct atfft_dft* atfft_dft_create (int size, enum atfft_direction direction, enum atfft_format format)
{
    struct atfft_dft *fft;

    if (!(fft = calloc (1, sizeof (*fft))))
        return NULL;

    fft->size = size;
    fft->direction = direction;
    fft->format = format;

    /* calculate radices */
    int max_r = 0;
    int sub_sizes [MAX_RADICES];
    int n_radices = atfft_init_radices (size, fft->radices, sub_sizes, &max_r);

    /* calculate twiddle factors */
    fft->t_factors = malloc (size * sizeof (*(fft->t_factors)));
    atfft_init_stage_twiddle_factors (fft->stage_t_factors,
                                      fft->radices,
                                      sub_sizes,
                                      n_radices,
                                      direction);

    /* clean up on failure */
    if (!fft->t_factors)
        goto failed;
    else
        atfft_init_twiddle_factors (fft->t_factors, size, direction);

    /* allocate some working space */
    fft->work_space = malloc (max_r * sizeof (*(fft->work_space)));

    if (!fft->work_space)
        goto failed;

    /* create any necessary Rader's algorithm strucs */
    fft->raders = atfft_init_raders (fft->radices,
                                     &(fft->n_raders),
                                     fft->radix_raders,
                                     direction,
                                     format);

    if (fft->n_raders > 0 && !fft->raders)
        goto failed;

    return fft;

failed:
    atfft_dft_destroy (fft);
    return NULL;
}

void atfft_dft_destroy (struct atfft_dft *fft)
{
    if (fft)
    {
        atfft_free_raders (fft->raders, fft->n_raders);
        free (fft->work_space);
        atfft_free_twiddle_factors (fft->stage_t_factors);
        free (fft->t_factors);
        free (fft);
    }
}

static void atfft_dft_2 (atfft_complex *out,
                         int sub_size)
{
    /* The start of each block in the output DFT */
    atfft_complex *bin1 = out;
    atfft_complex *bin2 = bin1 + sub_size;
    atfft_complex t;

    ATFFT_COPY_COMPLEX (*bin2, t);
    ATFFT_DIFFERENCE_COMPLEX (*bin1, t, *bin2);
    ATFFT_SUM_COMPLEX (*bin1, t, *bin1);
}

static void atfft_dft_3 (atfft_complex *out,
                         int sub_size,
                         atfft_sample sin_2pi_on_3)
{
    atfft_complex *bins [3];
    bins [0] = out;
    bins [1] = bins [0] + sub_size;
    bins [2] = bins [1] + sub_size;

    atfft_complex ts [3];

    ATFFT_SUM_COMPLEX (*bins [1], *bins [2], ts [0]);
    ATFFT_REAL (ts [1]) = ATFFT_REAL (*bins [0]) - ATFFT_REAL (ts [0]) / 2.0;
    ATFFT_IMAG (ts [1]) = ATFFT_IMAG (*bins [0]) - ATFFT_IMAG (ts [0]) / 2.0;
    ATFFT_DIFFERENCE_COMPLEX (*bins [1], *bins [2], ts [2]);
    ATFFT_REAL (ts [2]) = sin_2pi_on_3 * ATFFT_REAL (ts [2]);
    ATFFT_IMAG (ts [2]) = sin_2pi_on_3 * ATFFT_IMAG (ts [2]);

    ATFFT_SUM_COMPLEX (*bins [0], ts [0], *bins [0]);
    ATFFT_REAL (*bins [1]) = ATFFT_REAL (ts [1]) - ATFFT_IMAG (ts [2]);
    ATFFT_IMAG (*bins [1]) = ATFFT_IMAG (ts [1]) + ATFFT_REAL (ts [2]);
    ATFFT_REAL (*bins [2]) = ATFFT_REAL (ts [1]) + ATFFT_IMAG (ts [2]);
    ATFFT_IMAG (*bins [2]) = ATFFT_IMAG (ts [1]) - ATFFT_REAL (ts [2]);
}

static void atfft_dft_4 (atfft_complex **bins,
                         enum atfft_direction direction)
{
    atfft_complex ts [4];

    ATFFT_SUM_COMPLEX (*bins [0], *bins [2], ts [0]);
    ATFFT_SUM_COMPLEX (*bins [1], *bins [3], ts [1]);
    ATFFT_DIFFERENCE_COMPLEX (*bins [0], *bins [2], ts [2]);

    switch (direction)
    {
        case ATFFT_BACKWARD:
            ATFFT_DIFFERENCE_COMPLEX (*bins [3], *bins [1], ts [3]);
            break;
        default:
            ATFFT_DIFFERENCE_COMPLEX (*bins [1], *bins [3], ts [3]);
    }

    ATFFT_SUM_COMPLEX (ts [0], ts[1], *bins [0]);
    ATFFT_REAL (*bins [1]) = ATFFT_REAL (ts [2]) + ATFFT_IMAG (ts [3]);
    ATFFT_IMAG (*bins [1]) = ATFFT_IMAG (ts [2]) - ATFFT_REAL (ts [3]);
    ATFFT_DIFFERENCE_COMPLEX (ts [0], ts[1], *bins [2]);
    ATFFT_REAL (*bins [3]) = ATFFT_REAL (ts [2]) - ATFFT_IMAG (ts [3]);
    ATFFT_IMAG (*bins [3]) = ATFFT_IMAG (ts [2]) + ATFFT_REAL (ts [3]);
}

static void atfft_butterfly_n (atfft_complex *out,
                               int top_size,
                               int sub_size,
                               int stride,
                               int radix,
                               atfft_complex *t_factors,
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
            ATFFT_COPY_COMPLEX (out [n * sub_size + i], work_space [n]);
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

            /* copy the ith bin of the first sub-transform to the current
             * output bin.
             */
            ATFFT_COPY_COMPLEX (work_space [0], out [k]);

            /* Sum in the ith bins from the remaining sub-transforms,
             * multiplied by their respective twiddle factor.
             * out[k] += work_space [r] * t_factors [(k * r * stride) % top_size]
             */
            for (int r = 1; r < radix; ++r)
            {
                atfft_sample *bin = out [k];
                atfft_sample *in = work_space [r];
                atfft_sample *t = t_factors [(k * r * stride) % top_size];

                ATFFT_REAL (bin) += ATFFT_REAL (in) * ATFFT_REAL (t) -
                                    ATFFT_IMAG (in) * ATFFT_IMAG (t);
                ATFFT_IMAG (bin) += ATFFT_REAL (in) * ATFFT_IMAG (t) +
                                    ATFFT_IMAG (in) * ATFFT_REAL (t);
            }
        }
    }
}

static void atfft_butterfly_2 (atfft_complex *out,
                               int sub_size,
                               atfft_complex *t_factors)
{
    int i = sub_size;
    int t = 0;

    while (i--)
    {
        ATFFT_MULTIPLY_BY_COMPLEX (out [sub_size], t_factors [t]);
        atfft_dft_2 (out, sub_size);

        ++out;
        ++t;
    }
}

static void atfft_butterfly_3 (const struct atfft_dft *fft,
                               atfft_complex *out,
                               int sub_size,
                               int stride,
                               atfft_complex *t_factors)
{
    int i = sub_size;
    int radix = 3;
    int t = 0;

    atfft_sample sin_2pi_on_3 = ATFFT_IMAG (fft->t_factors [sub_size * stride]);

    while (i--)
    {
        int m = sub_size;

        for (int n = 1; n < radix; ++n)
        {
            ATFFT_MULTIPLY_BY_COMPLEX (out [m], t_factors [t]);
            m += sub_size;
            ++t;
        }

        atfft_dft_3 (out, sub_size, sin_2pi_on_3);

        ++out;
    }
}

static void atfft_butterfly_4 (atfft_complex *out,
                               int sub_size,
                               enum atfft_direction direction,
                               atfft_complex *t_factors)
{
    int i = sub_size;
    int radix = 4;
    int t = 0;

    atfft_complex *bins [radix];

    bins [0] = out;
    bins [1] = bins [0] + sub_size;
    bins [2] = bins [1] + sub_size;
    bins [3] = bins [2] + sub_size;

    while (i--)
    {
        for (int n = 1; n < radix; ++n)
        {
            ATFFT_MULTIPLY_BY_COMPLEX (*bins [n], t_factors [t]);
            ++t;
        }

        atfft_dft_4 (bins, direction);

        for (int n = 0; n < radix; ++n)
        {
            ++(bins [n]);
        }
    }
}

static void atfft_butterfly_rader (atfft_complex *out,
                                   int sub_size,
                                   int radix,
                                   atfft_complex *t_factors,
                                   struct atfft_dft_rader *rader)
{
    int i = sub_size;
    int t = 0;

    while (i--)
    {
        int m = sub_size;

        for (int n = 1; n < radix; ++n)
        {
            ATFFT_MULTIPLY_BY_COMPLEX (out [m], t_factors [t]);
            m += sub_size;
            ++t;
        }

        atfft_dft_rader_complex_transform (rader, out, out, sub_size);

        ++out;
    }
}

static void atfft_butterfly (const struct atfft_dft *fft,
                             atfft_complex *out,
                             int sub_size,
                             int stride,
                             int radix,
                             struct atfft_dft_rader *rader,
                             atfft_complex *t_factors)
{
    switch (radix)
    {
        case 2:
            atfft_butterfly_2 (out, sub_size, t_factors);
            break;
        case 3:
            atfft_butterfly_3 (fft, out, sub_size, stride, t_factors);
            break;
        case 4:
            atfft_butterfly_4 (out, sub_size, fft->direction, t_factors);
            break;
        default:
            if (rader)
                atfft_butterfly_rader (out, sub_size, radix, t_factors, rader);
            else
                atfft_butterfly_n (out, fft->size, sub_size, stride, radix, fft->t_factors, fft->work_space);
    }
}

static void atfft_compute_dft_complex (const struct atfft_dft *fft,
                                       atfft_complex *in,
                                       atfft_complex *out,
                                       int sub_size,
                                       int stride,
                                       const int *radices,
                                       struct atfft_dft_rader **raders,
                                       atfft_complex **t_factors)
{
    /* Get the radix, R, for this stage of the transform.
     * We will split the transform into R sub-transforms
     * of size next_size.
     */
    int R = *radices++;
    int next_size = sub_size / R;

    if (*radices)
    {
        /* If there is another radix in the list
         * we recursively apply the transform, once
         * for each sub-transform in this stage. 
         */
        for (int r = 0; r < R; ++r)
        {
            atfft_compute_dft_complex (fft, 
                                       in + r * stride,
                                       out + r * next_size,
                                       next_size,
                                       stride * R,
                                       radices,
                                       raders + 1,
                                       t_factors + 1);
        }
    }
    else
    {
        /* If there are no more radices, we have
         * reached the first stage of our transform.
         * Here we apply the decimation in time.
         */
        for (int i = 0; i < sub_size; ++i)
        {
            ATFFT_COPY_COMPLEX (in [i * stride], out [i]);
        }
    }

    /* Apply butterfly for this stage of the transform. */
    atfft_butterfly (fft, out, next_size, stride, R, *raders, *t_factors);
}

void atfft_dft_complex_transform (struct atfft_dft *fft, atfft_complex *in, atfft_complex *out)
{
    /* Only to be used with complex FFTs. */
    assert (fft->format == ATFFT_COMPLEX);

    atfft_compute_dft_complex (fft,
                               in,
                               out,
                               fft->size,
                               1,
                               fft->radices,
                               fft->radix_raders,
                               fft->stage_t_factors);
}

void atfft_dft_calculate_bin_real (struct atfft_dft *fft,
                                   const atfft_sample *in,
                                   atfft_complex *out,
                                   int k)
{
    ATFFT_REAL (*out) = 0.0;
    ATFFT_IMAG (*out) = 0.0;

    for (int i = 0; i < fft->size; ++i)
    {
        atfft_complex *t_factor = &(fft->t_factors [(k * i) % fft->size]);
        ATFFT_REAL (*out) += in [i] * ATFFT_REAL (*t_factor);
        ATFFT_IMAG (*out) += in [i] * ATFFT_IMAG (*t_factor);
    }
}

void atfft_dft_real_forward_transform (struct atfft_dft *fft, const atfft_sample *in, atfft_complex *out)
{
    /* Only to be used for forward real FFTs. */
    assert ((fft->format == ATFFT_REAL) && (fft->direction == ATFFT_FORWARD));

    for (int k = 0; k < (fft->size / 2) + 1; ++k)
    {
        atfft_dft_calculate_bin_real (fft, in, out + k, k);
    }
}

void atfft_dft_calculate_sample_real (struct atfft_dft *fft,
                                      atfft_complex *in,
                                      atfft_sample *out,
                                      int i)
{
    int k = 0;
    atfft_complex *t_factor = NULL;

    /* get DC component */
    *out = ATFFT_REAL (in [0]);

    /* add other frequency components */
    for (k = 1; k < (fft->size / 2) + 1; ++k)
    {
        t_factor = &(fft->t_factors [(k * i) % fft->size]);
        *out += 2.0 * ATFFT_REAL (in [k]) * ATFFT_REAL (*t_factor);
        *out -= 2.0 * ATFFT_IMAG (in [k]) * ATFFT_IMAG (*t_factor);
    }

    if (atfft_is_even (fft->size))
        *out -= ATFFT_REAL (in [k - 1]) * ATFFT_REAL (*t_factor);
}

void atfft_dft_real_backward_transform (struct atfft_dft *fft, atfft_complex *in, atfft_sample *out)
{
    /* Only to be used for backward real FFTs. */
    assert ((fft->format == ATFFT_REAL) && (fft->direction == ATFFT_BACKWARD));

    for (int i = 0; i < fft->size; ++i)
    {
        atfft_dft_calculate_sample_real (fft, in, out + i, i);
    }
}
