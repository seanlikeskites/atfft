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

extern void atfft_copy_complex (const atfft_complex x, atfft_complex y);
extern void atfft_sum_complex (const atfft_complex a,
                               const atfft_complex b,
                               atfft_complex s);
extern void atfft_difference_complex (const atfft_complex a,
                                      const atfft_complex b,
                                      atfft_complex d);
extern void atfft_product_complex (const atfft_complex a,
                                   const atfft_complex b,
                                   atfft_complex p);
extern void atfft_multiply_by_complex (atfft_complex a,
                                       const atfft_complex b);

/* Because an int is used to represent the size of the transform
 * valid sizes are anywhere between 0 and 2^(n - 1), where n is
 * the number of bits in an int. The maximum number of radices 
 * will therefore be (n - 1) as it will occur when the radices
 * are all 2s.
 */
#define MAX_RADICES (sizeof (int) * CHAR_BIT - 1)

struct atfft_dft
{
    int size;
    enum atfft_direction direction;
    enum atfft_format format;

    /* null terminated array of radices */
    int radices [MAX_RADICES + 1];

    /* twiddle factors */
    atfft_complex *tFactors;

    /* working space for length-n butterflies */
    atfft_complex *workSpace;
};

static void atfft_init_twiddle_factors (atfft_complex *factors,
                                        int size,
                                        enum atfft_direction direction)
{
    int i = 0;
    atfft_sample sinFactor = -1.0;

    if (direction == ATFFT_BACKWARD)
        sinFactor = 1.0;

    for (i = 0; i < size; ++i)
    {
        atfft_sample x = 2.0 * i * M_PI / size;
        ATFFT_REAL (factors [i]) = cos (x);
        ATFFT_IMAG (factors [i]) = sinFactor * sin (x);
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

static int atfft_init_radices (int size, int *radices)
{
    /* current radix */
    int r = 4;
    int maxR = 2;
    int sqrtSize = (int) sqrt ((double) size);

    /* factor out 2s, then other primes */
    do
    {
        while (size % r)
        {
            r = atfft_next_radix (r);

            /* a number will only have one prime factor greater than its square root */
            if (r > sqrtSize)
                r = size;
        }

        *radices++ = r;
        size /= r;

        if (r > maxR)
            maxR = r;
    }
    while (size > 1);

    return maxR;
}

struct atfft_dft* atfft_dft_create (int size, enum atfft_direction direction, enum atfft_format format)
{
    struct atfft_dft *fft;
    int maxR = 0;

    if (!(fft = calloc (1, sizeof (*fft))))
        return NULL;

    fft->size = size;
    fft->direction = direction;
    fft->format = format;

    /* calculate twiddle factors */
    fft->tFactors = malloc (size * sizeof (*(fft->tFactors)));

    /* clean up on failure */
    if (!fft->tFactors)
        goto failed;
    else
        atfft_init_twiddle_factors (fft->tFactors, size, direction);

    /* calculate radices */
    maxR = atfft_init_radices (size, fft->radices);

    /* allocate some working space */
    fft->workSpace = malloc (maxR * sizeof (*(fft->workSpace)));

    if (!fft->workSpace)
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
        free (fft->workSpace);
        free (fft->tFactors);
        free (fft);
    }
}

static void atfft_butterfly_2 (atfft_complex *out,
                               int subSize,
                               int stride,
                               atfft_complex *tFactors)
{
    /* The start of each block in the output DFT */
    atfft_complex *bin1 = out;
    atfft_complex *bin2 = bin1 + subSize;

    while (subSize--)
    {
        atfft_complex t;

        atfft_product_complex (*tFactors, *bin2, t);
        atfft_difference_complex (*bin1, t, *bin2);
        atfft_sum_complex (*bin1, t, *bin1);

        /* get the next twiddle factor and 
         * advance to the next output bins.
         */
        tFactors += stride;
        ++bin1;
        ++bin2;
    }
}

static void atfft_butterfly_3 (atfft_complex *out,
                               int subSize,
                               int stride,
                               atfft_complex *tFactors)
{
    int i = subSize;
    int n = 0;

    atfft_complex *bins [3];
    bins [0] = out;
    bins [1] = bins [0] + subSize;
    bins [2] = bins [1] + subSize;

    atfft_complex *tfs [2] = {tFactors, tFactors};
    atfft_complex zs [2];
    atfft_complex ts [3];

    atfft_sample sinPiOn3 = ATFFT_IMAG (tFactors [subSize * stride]);

    while (i--)
    {
        for (n = 0; n < 2; ++n)
        {
            atfft_product_complex (*tfs [n], *bins [n + 1], zs [n]);
            tfs [n] += (n + 1) * stride;
        }

        atfft_sum_complex (zs [0], zs [1], ts [0]);
        ATFFT_REAL (ts [1]) = ATFFT_REAL (*bins [0]) - ATFFT_REAL (ts [0]) / 2.0;
        ATFFT_IMAG (ts [1]) = ATFFT_IMAG (*bins [0]) - ATFFT_IMAG (ts [0]) / 2.0;
        atfft_difference_complex (zs [0], zs [1], ts [2]);
        ATFFT_REAL (ts [2]) = sinPiOn3 * ATFFT_REAL (ts [2]);
        ATFFT_IMAG (ts [2]) = sinPiOn3 * ATFFT_IMAG (ts [2]);

        atfft_sum_complex (*bins [0], ts [0], *bins [0]);
        ATFFT_REAL (*bins [1]) = ATFFT_REAL (ts [1]) - ATFFT_IMAG (ts [2]);
        ATFFT_IMAG (*bins [1]) = ATFFT_IMAG (ts [1]) + ATFFT_REAL (ts [2]);
        ATFFT_REAL (*bins [2]) = ATFFT_REAL (ts [1]) + ATFFT_IMAG (ts [2]);
        ATFFT_IMAG (*bins [2]) = ATFFT_IMAG (ts [1]) - ATFFT_REAL (ts [2]);

        for (n = 0; n < 3; ++n)
        {
            ++(bins [n]);
        }
    }
}

static void atfft_butterfly_4 (atfft_complex *out,
                               int subSize,
                               int stride,
                               enum atfft_direction direction,
                               atfft_complex *tFactors)
{
    int i = subSize;
    int n = 0;

    atfft_complex *bins [4];
    bins [0] = out;
    bins [1] = bins [0] + subSize;
    bins [2] = bins [1] + subSize;
    bins [3] = bins [2] + subSize;

    atfft_complex *tfs [3] = {tFactors, tFactors, tFactors};
    atfft_complex zs [3];
    atfft_complex ts [4];

    while (i--)
    {
        for (n = 0; n < 3; ++n)
        {
            atfft_product_complex (*tfs [n], *bins [n + 1], zs [n]);
            tfs [n] += (n + 1) * stride;
        }

        atfft_sum_complex (*bins [0], zs [1], ts [0]);
        atfft_sum_complex (zs [0], zs [2], ts [1]);
        atfft_difference_complex (*bins [0], zs [1], ts [2]);

        switch (direction)
        {
            case ATFFT_BACKWARD:
                atfft_difference_complex (zs [2], zs [0], ts [3]);
                break;
            default:
                atfft_difference_complex (zs [0], zs [2], ts [3]);
        }

        atfft_sum_complex (ts [0], ts[1], *bins [0]);
        ATFFT_REAL (*bins [1]) = ATFFT_REAL (ts [2]) + ATFFT_IMAG (ts [3]);
        ATFFT_IMAG (*bins [1]) = ATFFT_IMAG (ts [2]) - ATFFT_REAL (ts [3]);
        atfft_difference_complex (ts [0], ts[1], *bins [2]);
        ATFFT_REAL (*bins [3]) = ATFFT_REAL (ts [2]) - ATFFT_IMAG (ts [3]);
        ATFFT_IMAG (*bins [3]) = ATFFT_IMAG (ts [2]) + ATFFT_REAL (ts [3]);

        for (n = 0; n < 4; ++n)
        {
            ++(bins [n]);
        }
    }
}

static void atfft_butterfly_n (atfft_complex *out,
                               int topSize,
                               int subSize,
                               int stride,
                               int radix,
                               atfft_complex *tFactors,
                               atfft_complex *workSpace)
{
    /* Combine radix DFTs of size subSize,
     * into one DTF of size (radix * subSize).
     */
    int i = 0;
    int n = 0;

    /* Loop over the bins in each of the sub-transforms. */
    for (i = 0; i < subSize; ++i)
    {
        /* Copy ith bin from each sub-transform into workSpace. */
        for (n = 0; n < radix; ++n)
        {
            atfft_copy_complex (out [n * subSize + i], workSpace [n]);
        }

        /* Calculate the output bins. */
        for (n = 0; n < radix; ++n)
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
            int k = n * subSize + i;
            int r = 0;

            /* copy the ith bin of the first sub-transform to the current
             * output bin.
             */
            atfft_copy_complex (workSpace [0], out [k]);

            /* Sum in the ith bins from the remaining sub-transforms,
             * multiplied by their respective twiddle factor.
             * out[k] += workSpace [r] * tFactors [(k * r * stride) % topSize]
             */
            for (r = 1; r < radix; ++r)
            {
                atfft_sample *bin = out [k];
                atfft_sample *in = workSpace [r];
                atfft_sample *t = tFactors [(k * r * stride) % topSize];

                ATFFT_REAL (bin) += ATFFT_REAL (in) * ATFFT_REAL (t) -
                                    ATFFT_IMAG (in) * ATFFT_IMAG (t);
                ATFFT_IMAG (bin) += ATFFT_REAL (in) * ATFFT_IMAG (t) +
                                    ATFFT_IMAG (in) * ATFFT_REAL (t);
            }
        }
    }
}

void atfft_butterfly (const struct atfft_dft *fft,
                      atfft_complex *out,
                      int subSize,
                      int stride,
                      int radix)
{
    switch (radix)
    {
        case 2:
            atfft_butterfly_2 (out, subSize, stride, fft->tFactors);
            break;
        case 3:
            atfft_butterfly_3 (out, subSize, stride, fft->tFactors);
            break;
        case 4:
            atfft_butterfly_4 (out, subSize, stride, fft->direction, fft->tFactors);
            break;
        default:
            atfft_butterfly_n (out, fft->size, subSize, stride, radix, fft->tFactors, fft->workSpace);
    }
}

static void atfft_compute_dft_complex (const struct atfft_dft *fft,
                                       atfft_complex *in,
                                       atfft_complex *out,
                                       int subSize,
                                       int stride,
                                       const int *radices)
{
    /* Get the radix, R, for this stage of the transform.
     * We will split the transform into R sub-transforms
     * of size nextSize.
     */
    int R = *radices++;
    int nextSize = subSize / R;

    if (*radices)
    {
        /* If there is another radix in the list
         * we recursively apply the transform, once
         * for each sub-transform in this stage. 
         */
        int r = 0;

        for (r = 0; r < R; ++r)
        {
            atfft_compute_dft_complex (fft, 
                                       in + r * stride,
                                       out + r * nextSize,
                                       nextSize,
                                       stride * R,
                                       radices);
        }
    }
    else
    {
        /* If there are no more radices, we have
         * reached the first stage of our transform.
         * Here we apply the decimation in time.
         */
        int i = 0;

        for (i = 0; i < subSize; ++i)
        {
            atfft_copy_complex (in [i * stride], out [i]);
        }
    }

    /* Apply butterfly for this stage of the transform. */
    atfft_butterfly (fft, out, nextSize, stride, R);
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
                               fft->radices);
}

void atfft_dft_calculate_bin_real (struct atfft_dft *fft,
                                   const atfft_sample *in,
                                   atfft_complex *out,
                                   int k)
{
    int i = 0;
    ATFFT_REAL (*out) = 0.0;
    ATFFT_IMAG (*out) = 0.0;

    for (i = 0; i < fft->size; ++i)
    {
        atfft_complex *tFactor = &(fft->tFactors [(k * i) % fft->size]);
        ATFFT_REAL (*out) += in [i] * ATFFT_REAL (*tFactor);
        ATFFT_IMAG (*out) += in [i] * ATFFT_IMAG (*tFactor);
    }
}

void atfft_dft_real_forward_transform (struct atfft_dft *fft, const atfft_sample *in, atfft_complex *out)
{
    int k = 0;

    /* Only to be used for forward real FFTs. */
    assert ((fft->format == ATFFT_REAL) && (fft->direction == ATFFT_FORWARD));

    for (k = 0; k < (fft->size / 2) + 1; ++k)
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
    atfft_complex *tFactor = NULL;

    /* get DC component */
    *out = ATFFT_REAL (in [0]);

    /* add other frequency components */
    for (k = 1; k < (fft->size / 2) + 1; ++k)
    {
        tFactor = &(fft->tFactors [(k * i) % fft->size]);
        *out += 2.0 * ATFFT_REAL (in [k]) * ATFFT_REAL (*tFactor);
        *out -= 2.0 * ATFFT_IMAG (in [k]) * ATFFT_IMAG (*tFactor);
    }

    if (atfft_is_even (fft->size))
        *out -= ATFFT_REAL (in [k - 1]) * ATFFT_REAL (*tFactor);
}

void atfft_dft_real_backward_transform (struct atfft_dft *fft, atfft_complex *in, atfft_sample *out)
{
    int i = 0;

    /* Only to be used for backward real FFTs. */
    assert ((fft->format == ATFFT_REAL) && (fft->direction == ATFFT_BACKWARD));

    for (i = 0; i < fft->size; ++i)
    {
        atfft_dft_calculate_sample_real (fft, in, out + i, i);
    }
}
