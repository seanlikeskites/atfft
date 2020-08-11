/*
 * Copyright (C) 2020 Sean Enderby <sean.enderby@gmail.com>
 *
 * This program is free software. It comes without any warranty, to
 * the extent permitted by applicable law. You can redistribute it 
 * and/or modify it under the terms of the Do What The Fuck You Want 
 * To Public License, Version 2, as published by Sam Hocevar. See 
 * the COPYING file for more details.
 */

#include <math.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <atfft/atfft_dft.h>
#include "atfft_dft_rader.h"
#include "../atfft_internal.h"

static int atfft_mod (int a, int n)
{
    int r = a % n;
    return r < 0 ? r + n : r;
}

static int atfft_is_prime (int x)
{
    int sqrtX = (int) sqrt ((double) x);
    int i = 2;

    if (x <= 1 || (atfft_is_even (x) && x > 2))
        return 0;

    for (i = 2; i <= sqrtX; ++i)
    {
        if ((x % i) == 0)
            return 0;
    }

    return 1;
}

static int atfft_next_power_of_2 (int x)
{
    if (x <= 0)
        return 0;
    else
        return pow (2, (int) log2 (x) + 1);
}

static int atfft_primitive_root_mod_n (int n)
{
    int g = 2;
    int isRoot = 1;

    /* this algorithm will only work for prime numbers */
    if (!atfft_is_prime (n))
        return -1;

    /* loop over potential primitive roots */
    for (g = 2; g < n; ++g)
    {
        int i = n - 2;
        int m = 1;
        isRoot = 1;

        while (i--)
        {
            m = (m * g) % n;

            if (m == 1)
            {
                isRoot = 0;
                break;
            }
        }

        if (isRoot)
            return g;
    }

    return -1;
}

static void atfft_gcd (int a, int b, int *gcd, int *x, int *y)
{
    int absA = abs (a);
    int absB = abs (b);
    int r0 = ATFFT_MAX (absA, absB);
    int r1 = ATFFT_MIN (absA, absB);
    int s0 = absA > absB ? 1 : 0;
    int s1 = 1 - s0;
    int t0 = s1;
    int t1 = s0;

    while (r1)
    {
        int q = r0 / r1;

        int r2 = r0 - q * r1;
        int s2 = s0 - q * s1;
        int t2 = t0 - q * t1;

        r0 = r1;
        r1 = r2;
        s0 = s1;
        s1 = s2;
        t0 = t1;
        t1 = t2;
    }

    *gcd = r0;
    *x = a < 0 ? -s0 : s0;
    *y = b < 0 ? -t0 : t0;
}

static int atfft_mult_inverse_mod_n (int a, int n)
{
    int gcd, x, y;

    atfft_gcd (a % n, n, &gcd, &x, &y);

    if (gcd == 1)
        return atfft_mod (x, n);
    else
        return -1;
}

static int atfft_rader_convolution_fft_size (int raderSize)
{
    if (atfft_is_power_of_2 (raderSize))
        return raderSize;
    else
        return atfft_next_power_of_2 (2 * raderSize - 1);
}

struct atfft_dft_rader
{
    int size;
    int raderSize;
    enum atfft_direction direction;
    enum atfft_format format;
    int pRoot1, pRoot2;
    int convSize;
    struct atfft_dft *convForward, *convBackward;
    int *perm1, *perm2;
    atfft_complex *sig, *sigDft, *conv, *convDft;
};

static void atfft_init_rader_permutations (int *perm, int size, int pRoot)
{
    int n = 0;
    int i = 1;

    for (n = 0; n < size - 1; ++n)
    {
        perm [n] = i;
        i = atfft_mod (i * pRoot, size);
    }
}

static int atfft_init_rader_convolution_dft (int size,
                                             enum atfft_direction direction,
                                             atfft_complex *convDft,
                                             int convSize,
                                             int *perm,
                                             int permSize,
                                             struct atfft_dft *fft)
{
    int i = 0;
    atfft_sample sinFactor = -1.0;
    atfft_complex *tFactors = calloc (convSize, sizeof (*tFactors));

    if (!tFactors)
        return -1;

    if (direction == ATFFT_BACKWARD)
        sinFactor = 1.0;

    /* produce rader twiddle factors */
    for (i = 0; i < permSize; ++i)
    {
        atfft_sample x = 2.0 * perm [i] * M_PI / size;
        ATFFT_REAL (tFactors [i]) = cos (x);
        ATFFT_IMAG (tFactors [i]) = sinFactor * sin (x);
    }

    /* replicate samples for circular convolution */
    if (convSize > permSize)
    {
        size_t nReplications = permSize - 1;
        atfft_complex *source = tFactors + 1;
        atfft_complex *dest = tFactors + convSize - nReplications;

        memcpy (dest, source, nReplications * sizeof (*dest));

    }

    /* take DFT of twiddle factors */
    atfft_dft_complex_transform (fft, tFactors, convDft);

    free (tFactors);
    return 0;
}

struct atfft_dft_rader* atfft_dft_rader_create (int size,
                                                enum atfft_direction direction,
                                                enum atfft_format format)
{
    struct atfft_dft_rader *fft;

    /* we can only find primitive roots for prime numbers */
    assert (atfft_is_prime (size));

    if (!(fft = calloc (1, sizeof (*fft))))
        return NULL;

    fft->size = size;
    fft->raderSize = size - 1;
    fft->direction = direction;
    fft->format = format;
    fft->pRoot1 = atfft_primitive_root_mod_n (size);
    fft->pRoot2 = atfft_mult_inverse_mod_n (fft->pRoot1, size);

    /* allocate some regular dft objects for performing the convolution */
    fft->convSize = atfft_rader_convolution_fft_size (fft->raderSize);
    fft->convForward = atfft_dft_create (fft->convSize, ATFFT_FORWARD, ATFFT_COMPLEX);
    fft->convBackward = atfft_dft_create (fft->convSize, ATFFT_BACKWARD, ATFFT_COMPLEX);

    if (!(fft->convForward && fft->convBackward))
        goto failed;

    /* generate permutation indices */
    fft->perm1 = malloc (fft->raderSize * sizeof (*(fft->perm1)));
    fft->perm2 = malloc (fft->raderSize * sizeof (*(fft->perm2)));

    if (!(fft->perm1 && fft->perm2))
        goto failed;

    atfft_init_rader_permutations (fft->perm1, size, fft->pRoot1);
    atfft_init_rader_permutations (fft->perm2, size, fft->pRoot2);

    /* allocate work space for performing the dft */
    fft->sig = calloc (fft->convSize, sizeof (*(fft->sig))); /* set up for zero padding */
    fft->sigDft = malloc (fft->convSize * sizeof (*(fft->sigDft)));
    fft->conv = malloc (fft->convSize * sizeof (*(fft->conv)));
    fft->convDft = malloc (fft->convSize * sizeof (*(fft->convDft)));

    if (!(fft->sig && fft->sigDft && fft->convDft))
        goto failed;

    /* calculate the convolution dft */
    if (atfft_init_rader_convolution_dft (size,
                                          direction,
                                          fft->convDft,
                                          fft->convSize,
                                          fft->perm2,
                                          fft->raderSize,
                                          fft->convForward) < 0)
        goto failed;

    return fft;

failed:
    atfft_dft_rader_destroy (fft);
    return NULL;
}

void atfft_dft_rader_destroy (struct atfft_dft_rader *fft)
{
    if (fft)
    {
        free (fft->convDft);
        free (fft->conv);
        free (fft->sigDft);
        free (fft->sig);
        free (fft->perm2);
        free (fft->perm1);
        atfft_dft_destroy (fft->convBackward);
        atfft_dft_destroy (fft->convForward);
        free (fft);
    }
}

static void atfft_rader_permute_input (int *perm,
                                       int size,
                                       atfft_complex *in, 
                                       int inStride,
                                       atfft_complex *out)
{
    int i = 0;

    for (i = 0; i < size; ++i)
    {
        atfft_copy_complex (in [inStride * perm [i]], out [i]);
    }
}

static void atfft_rader_permute_output (int *perm,
                                        int size,
                                        atfft_complex *in, 
                                        atfft_complex *out,
                                        int outStride)
{
    int i = 0;

    for (i = 0; i < size; ++i)
    {
        atfft_copy_complex (in [i], out [outStride * perm [i]]);
    }
}

void atfft_dft_rader_complex_transform (struct atfft_dft_rader *fft,
                                        atfft_complex *in,
                                        atfft_complex *out,
                                        int stride)
{
    int i = 0;
    atfft_complex in0, out0;

    atfft_copy_complex (in [0], in0);
    atfft_copy_complex (in0, out0);

    atfft_rader_permute_input (fft->perm1, 
                               fft->raderSize,
                               in,
                               stride,
                               fft->sig);

    atfft_dft_complex_transform (fft->convForward, fft->sig, fft->sigDft);

    for (i = 0; i < fft->convSize; ++i)
    {
        atfft_multiply_by_complex (fft->sigDft [i], fft->convDft [i]);
        atfft_sum_complex (out0, fft->sig [i], out0);
    }

    ATFFT_REAL (fft->sigDft [0]) += ATFFT_REAL (in0) * fft->convSize;
    ATFFT_IMAG (fft->sigDft [0]) += ATFFT_IMAG (in0) * fft->convSize;

    atfft_dft_complex_transform (fft->convBackward, fft->sigDft, fft->conv);
    atfft_normalise_complex (fft->conv, fft->convSize); /* can probably move this into the setup */

    atfft_rader_permute_output (fft->perm2, 
                                fft->raderSize,
                                fft->conv,
                                out,
                                stride);

    atfft_copy_complex (out0, out[0]);
}
