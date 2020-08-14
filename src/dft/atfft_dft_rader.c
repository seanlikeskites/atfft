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
    if (x <= 1 || (atfft_is_even (x) && x > 2))
        return 0;

    int sqrt_x = (int) sqrt ((double) x);

    for (int i = 2; i <= sqrt_x; ++i)
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
    /* this algorithm will only work for prime numbers */
    if (!atfft_is_prime (n))
        return -1;

    /* loop over potential primitive roots */
    for (int g = 2; g < n; ++g)
    {
        int i = n - 2;
        int m = 1;
        int is_root = 1;

        while (i--)
        {
            m = (m * g) % n;

            if (m == 1)
            {
                is_root = 0;
                break;
            }
        }

        if (is_root)
            return g;
    }

    return -1;
}

static void atfft_gcd (int a, int b, int *gcd, int *x, int *y)
{
    int abs_a = abs (a);
    int abs_b = abs (b);
    int r0 = ATFFT_MAX (abs_a, abs_b);
    int r1 = ATFFT_MIN (abs_a, abs_b);
    int s0 = abs_a > abs_b ? 1 : 0;
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

static int atfft_rader_convolution_fft_size (int rader_size)
{
    if (atfft_is_power_of_2 (rader_size))
        return rader_size;
    else
        return atfft_next_power_of_2 (2 * rader_size - 1);
}

struct atfft_dft_rader
{
    int size;
    int rader_size;
    enum atfft_direction direction;
    enum atfft_format format;
    int p_root1, p_root2;
    int conv_size;
    struct atfft_dft *conv_forward, *conv_backward;
    int *perm1, *perm2;
    atfft_complex *sig, *sig_dft, *conv, *conv_dft;
};

static void atfft_init_rader_permutations (int *perm, int size, int p_root)
{
    int i = 1;

    for (int n = 0; n < size - 1; ++n)
    {
        perm [n] = i;
        i = atfft_mod (i * p_root, size);
    }
}

static int atfft_init_rader_convolution_dft (int size,
                                             enum atfft_direction direction,
                                             atfft_complex *conv_dft,
                                             int conv_size,
                                             int *perm,
                                             int perm_size,
                                             struct atfft_dft *fft)
{
    atfft_complex *t_factors = calloc (conv_size, sizeof (*t_factors));

    if (!t_factors)
        return -1;

    atfft_sample sin_factor = -1.0;

    if (direction == ATFFT_BACKWARD)
        sin_factor = 1.0;

    /* produce rader twiddle factors */
    for (int i = 0; i < perm_size; ++i)
    {
        atfft_sample x = 2.0 * perm [i] * M_PI / size;
        ATFFT_REAL (t_factors [i]) = cos (x) / conv_size;
        ATFFT_IMAG (t_factors [i]) = sin_factor * sin (x) / conv_size;
    }

    /* replicate samples for circular convolution */
    if (conv_size > perm_size)
    {
        size_t n_replications = perm_size - 1;
        atfft_complex *source = t_factors + 1;
        atfft_complex *dest = t_factors + conv_size - n_replications;

        memcpy (dest, source, n_replications * sizeof (*dest));

    }

    /* take DFT of twiddle factors */
    atfft_dft_complex_transform (fft, t_factors, conv_dft);

    free (t_factors);
    return 0;
}

struct atfft_dft_rader* atfft_dft_rader_create (int size,
                                                enum atfft_direction direction,
                                                enum atfft_format format)
{
    /* we can only find primitive roots for prime numbers */
    assert (atfft_is_prime (size));

    struct atfft_dft_rader *fft;

    if (!(fft = calloc (1, sizeof (*fft))))
        return NULL;

    fft->size = size;
    fft->rader_size = size - 1;
    fft->direction = direction;
    fft->format = format;
    fft->p_root1 = atfft_primitive_root_mod_n (size);
    fft->p_root2 = atfft_mult_inverse_mod_n (fft->p_root1, size);

    /* allocate some regular dft objects for performing the convolution */
    fft->conv_size = atfft_rader_convolution_fft_size (fft->rader_size);
    fft->conv_forward = atfft_dft_create (fft->conv_size, ATFFT_FORWARD, ATFFT_COMPLEX);
    fft->conv_backward = atfft_dft_create (fft->conv_size, ATFFT_BACKWARD, ATFFT_COMPLEX);

    if (!(fft->conv_forward && fft->conv_backward))
        goto failed;

    /* generate permutation indices */
    fft->perm1 = malloc (fft->rader_size * sizeof (*(fft->perm1)));
    fft->perm2 = malloc (fft->rader_size * sizeof (*(fft->perm2)));

    if (!(fft->perm1 && fft->perm2))
        goto failed;

    atfft_init_rader_permutations (fft->perm1, size, fft->p_root1);
    atfft_init_rader_permutations (fft->perm2, size, fft->p_root2);

    /* allocate work space for performing the dft */
    fft->sig = calloc (fft->conv_size, sizeof (*(fft->sig))); /* set up for zero padding */
    fft->sig_dft = malloc (fft->conv_size * sizeof (*(fft->sig_dft)));
    fft->conv = malloc (fft->conv_size * sizeof (*(fft->conv)));
    fft->conv_dft = malloc (fft->conv_size * sizeof (*(fft->conv_dft)));

    if (!(fft->sig && fft->sig_dft && fft->conv_dft))
        goto failed;

    /* calculate the convolution dft */
    if (atfft_init_rader_convolution_dft (size,
                                          direction,
                                          fft->conv_dft,
                                          fft->conv_size,
                                          fft->perm2,
                                          fft->rader_size,
                                          fft->conv_forward) < 0)
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
        free (fft->conv_dft);
        free (fft->conv);
        free (fft->sig_dft);
        free (fft->sig);
        free (fft->perm2);
        free (fft->perm1);
        atfft_dft_destroy (fft->conv_backward);
        atfft_dft_destroy (fft->conv_forward);
        free (fft);
    }
}

static void atfft_rader_permute_input (int *perm,
                                       int size,
                                       atfft_complex *in, 
                                       int in_stride,
                                       atfft_complex *out)
{
    for (int i = 0; i < size; ++i)
    {
        ATFFT_COPY_COMPLEX (in [in_stride * perm [i]], out [i]);
    }
}

static void atfft_rader_permute_output (int *perm,
                                        int size,
                                        atfft_complex *in, 
                                        atfft_complex *out,
                                        int out_stride)
{
    for (int i = 0; i < size; ++i)
    {
        ATFFT_COPY_COMPLEX (in [i], out [out_stride * perm [i]]);
    }
}

void atfft_dft_rader_complex_transform (struct atfft_dft_rader *fft,
                                        atfft_complex *in,
                                        atfft_complex *out,
                                        int stride)
{
    atfft_complex in0, out0;

    ATFFT_COPY_COMPLEX (in [0], in0);
    ATFFT_COPY_COMPLEX (in0, out0);

    atfft_rader_permute_input (fft->perm1, 
                               fft->rader_size,
                               in,
                               stride,
                               fft->sig);

    atfft_dft_complex_transform (fft->conv_forward, fft->sig, fft->sig_dft);

    for (int i = 0; i < fft->conv_size; ++i)
    {
        ATFFT_MULTIPLY_BY_COMPLEX (fft->sig_dft [i], fft->conv_dft [i]);
        ATFFT_SUM_COMPLEX (out0, fft->sig [i], out0);
    }

    ATFFT_REAL (fft->sig_dft [0]) += ATFFT_REAL (in0);
    ATFFT_IMAG (fft->sig_dft [0]) += ATFFT_IMAG (in0);

    atfft_dft_complex_transform (fft->conv_backward, fft->sig_dft, fft->conv);

    atfft_rader_permute_output (fft->perm2, 
                                fft->rader_size,
                                fft->conv,
                                out,
                                stride);

    ATFFT_COPY_COMPLEX (out0, out[0]);
}
