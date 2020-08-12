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
    atfft_complex *tFactors;

    /* working space for butterflies */
    atfft_complex *workSpace;

    /* sub fft objects for Rader's algorithm */
    int nRaders;
    struct atfft_dft_rader **raders;
    struct atfft_dft_rader *radixRaders [MAX_RADICES];
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

static int atfft_integer_is_in_array (const int *arr, int size, int member)
{
    int i = 0;

    for (i = 0; i < size; ++i)
    {
        if (arr [i] == member)
            return 1;
    }

    return 0;
}

static int atfft_rader_radices (const int *radices, int *raderRadices, int size)
{
    int i = 0;
    int nRaders = 0;

    for (i = 0; i < size; ++i)
    {
        if (radices [i] > ATFFT_RADER_THRESHOLD && 
            !atfft_integer_is_in_array (raderRadices, nRaders, radices [i]))
            raderRadices [nRaders++] = radices [i];
    }

    return nRaders;
}

static struct atfft_dft_rader* atfft_populate_rader (int radix,
                                                     const int *radices,
                                                     struct atfft_dft_rader **radixRaders,
                                                     enum atfft_direction direction,
                                                     enum atfft_format format)
{
    struct atfft_dft_rader *rader = atfft_dft_rader_create (radix, direction, format);
    int i = 0;

    if (!rader)
        return NULL;

    for (i = 0; i < MAX_RADICES; ++i)
    {
        if (radices [i] == radix)
            radixRaders [i] = rader;
    }

    return rader;
}

static void atfft_free_raders (struct atfft_dft_rader **raders,
                               int size)
{
    if (raders)
    {
        int i = 0;

        for (i = 0; i < size; ++i)
        {
            atfft_dft_rader_destroy (raders [i]);
        }

        free (raders);
    }
}

static struct atfft_dft_rader** atfft_init_raders (const int *radices,
                                                   int *nRaders,
                                                   struct atfft_dft_rader **radixRaders,
                                                   enum atfft_direction direction,
                                                   enum atfft_format format)
{
    int raderRadices [MAX_RADICES];
    struct atfft_dft_rader **raders = NULL;
    int i = 0;

    /* get unique radices for which Rader's algorithm is required */
    *nRaders = atfft_rader_radices (radices, raderRadices, MAX_RADICES);

    if (*nRaders == 0)
        return NULL;

    /* allocate some space for them */
    raders = calloc (*nRaders, sizeof (*raders));

    if (!raders)
        return NULL;

    /* create the rader dft structs */
    for (i = 0; i < *nRaders; ++i)
    {
        struct atfft_dft_rader *rader = atfft_populate_rader (raderRadices [i],
                                                              radices,
                                                              radixRaders,
                                                              direction,
                                                              format);

        if (!rader)
            goto failed;

        raders [i] = rader;
    }

    return raders;

failed:
    atfft_free_raders (raders, *nRaders);
    return NULL;
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

    /* create any necessary Rader's algorithm strucs */
    fft->raders = atfft_init_raders (fft->radices,
                                     &(fft->nRaders),
                                     fft->radixRaders,
                                     direction,
                                     format);

    if (fft->nRaders > 0 && !fft->raders)
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
        atfft_free_raders (fft->raders, fft->nRaders);
        free (fft->workSpace);
        free (fft->tFactors);
        free (fft);
    }
}

static inline void atfft_dft_2 (atfft_complex *out,
                                int subSize)
{
    /* The start of each block in the output DFT */
    atfft_complex *bin1 = out;
    atfft_complex *bin2 = bin1 + subSize;
    atfft_complex t;

    ATFFT_COPY_COMPLEX (*bin2, t);
    ATFFT_DIFFERENCE_COMPLEX (*bin1, t, *bin2);
    ATFFT_SUM_COMPLEX (*bin1, t, *bin1);
}

static inline void atfft_dft_3 (atfft_complex *out,
                                int subSize,
                                int stride,
                                atfft_complex *tFactors)
{
    atfft_complex *bins [3];
    bins [0] = out;
    bins [1] = bins [0] + subSize;
    bins [2] = bins [1] + subSize;

    atfft_complex ts [3];

    atfft_sample sinPiOn3 = ATFFT_IMAG (tFactors [subSize * stride]);

    ATFFT_SUM_COMPLEX (*bins [1], *bins [2], ts [0]);
    ATFFT_REAL (ts [1]) = ATFFT_REAL (*bins [0]) - ATFFT_REAL (ts [0]) / 2.0;
    ATFFT_IMAG (ts [1]) = ATFFT_IMAG (*bins [0]) - ATFFT_IMAG (ts [0]) / 2.0;
    ATFFT_DIFFERENCE_COMPLEX (*bins [1], *bins [2], ts [2]);
    ATFFT_REAL (ts [2]) = sinPiOn3 * ATFFT_REAL (ts [2]);
    ATFFT_IMAG (ts [2]) = sinPiOn3 * ATFFT_IMAG (ts [2]);

    ATFFT_SUM_COMPLEX (*bins [0], ts [0], *bins [0]);
    ATFFT_REAL (*bins [1]) = ATFFT_REAL (ts [1]) - ATFFT_IMAG (ts [2]);
    ATFFT_IMAG (*bins [1]) = ATFFT_IMAG (ts [1]) + ATFFT_REAL (ts [2]);
    ATFFT_REAL (*bins [2]) = ATFFT_REAL (ts [1]) + ATFFT_IMAG (ts [2]);
    ATFFT_IMAG (*bins [2]) = ATFFT_IMAG (ts [1]) - ATFFT_REAL (ts [2]);
}

static inline void atfft_dft_4 (atfft_complex **bins,
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
            ATFFT_COPY_COMPLEX (out [n * subSize + i], workSpace [n]);
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
            ATFFT_COPY_COMPLEX (workSpace [0], out [k]);

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

static inline void atfft_butterfly_2 (const struct atfft_dft *fft,
                                      atfft_complex *out,
                                      int subSize,
                                      int stride)
{
    int i = subSize;
    int t = 0;

    atfft_dft_2 (out, subSize);

    while (--i)
    {
        ++out;
        t += stride;

        ATFFT_MULTIPLY_BY_COMPLEX (out [subSize], fft->tFactors [t]);
        atfft_dft_2 (out, subSize);
    }
}

static inline void atfft_butterfly_3 (const struct atfft_dft *fft,
                                      atfft_complex *out,
                                      int subSize,
                                      int stride)
{
    int radix = 3;
    int i = 0;
    int r = 0;

    atfft_dft_3 (out, subSize, stride, fft->tFactors);

    for (i = stride; i < subSize * stride; i += stride)
    {
        int m = subSize;
        int n = i;
        
        ++out;

        for (r = 1; r < radix; ++r)
        {
            ATFFT_MULTIPLY_BY_COMPLEX (out [m], fft->tFactors [n]);
            m += subSize;
            n += i;
        }

        atfft_dft_3 (out, subSize, stride, fft->tFactors);
    }
}

static inline void atfft_butterfly_4 (const struct atfft_dft *fft,
                                      atfft_complex *out,
                                      int subSize,
                                      int stride)
{
    int i = subSize;
    int radix = 4;

    atfft_complex *bins [4];
    atfft_complex *tfs [3] = {fft->tFactors, fft->tFactors, fft->tFactors};

    bins [0] = out;
    bins [1] = bins [0] + subSize;
    bins [2] = bins [1] + subSize;
    bins [3] = bins [2] + subSize;

    atfft_dft_4 (bins, fft->direction);

    while (--i)
    {
        int n = 0;

        ++(bins [0]);

        for (n = 1; n < radix; ++n)
        {
            ++(bins [n]);
            tfs [n - 1] += n * stride;
            ATFFT_MULTIPLY_BY_COMPLEX (*bins [n], *tfs [n - 1]);
        }

        atfft_dft_4 (bins, fft->direction);
    }
}

static inline void atfft_butterfly_rader (const struct atfft_dft *fft,
                                          atfft_complex *out,
                                          int subSize,
                                          int stride,
                                          int radix,
                                          struct atfft_dft_rader *rader)
{
    int i = 0;
    int r = 0;

    atfft_dft_rader_complex_transform (rader, out, out, subSize);

    for (i = stride; i < subSize * stride; i += stride)
    {
        int m = subSize;
        int n = i;
        
        ++out;

        for (r = 1; r < radix; ++r)
        {
            ATFFT_MULTIPLY_BY_COMPLEX (out [m], fft->tFactors [n]);
            m += subSize;
            n += i;
        }

        atfft_dft_rader_complex_transform (rader, out, out, subSize);
    }
}

void atfft_butterfly (const struct atfft_dft *fft,
                      atfft_complex *out,
                      int subSize,
                      int stride,
                      int radix,
                      struct atfft_dft_rader *rader)
{
    switch (radix)
    {
        case 2:
            atfft_butterfly_2 (fft, out, subSize, stride);
            break;
        case 3:
            atfft_butterfly_3 (fft, out, subSize, stride);
            break;
        case 4:
            atfft_butterfly_4 (fft, out, subSize, stride);
            break;
        default:
            if (rader)
                atfft_butterfly_rader (fft, out, subSize, stride, radix, rader);
            else
                atfft_butterfly_n (out, fft->size, subSize, stride, radix, fft->tFactors, fft->workSpace);
    }
}

static void atfft_compute_dft_complex (const struct atfft_dft *fft,
                                       atfft_complex *in,
                                       atfft_complex *out,
                                       int subSize,
                                       int stride,
                                       const int *radices,
                                       struct atfft_dft_rader **raders)
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
                                       radices,
                                       raders + 1);
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
            ATFFT_COPY_COMPLEX (in [i * stride], out [i]);
        }
    }

    /* Apply butterfly for this stage of the transform. */
    atfft_butterfly (fft, out, nextSize, stride, R, *raders);
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
                               fft->radixRaders);
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
