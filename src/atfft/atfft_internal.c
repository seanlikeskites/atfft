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

#include "atfft_internal.h"
#include <atfft/dft.h>
#include <math.h>
#include <stdlib.h>

/******************************************
 * External definitions of inline functions.
 ******************************************/
extern inline void atfft_copy_complex (const atfft_complex x, atfft_complex *y);
extern inline void atfft_swap_complex (const atfft_complex x, atfft_complex *y);
extern inline void atfft_sum_complex (const atfft_complex a,
                                      const atfft_complex b,
                                      atfft_complex *s);
extern inline void atfft_difference_complex (const atfft_complex a,
                                             const atfft_complex b,
                                             atfft_complex *d);
extern inline void atfft_product_complex (const atfft_complex a,
                                          const atfft_complex b,
                                          atfft_complex *p);
extern inline void atfft_multiply_by_complex (atfft_complex *a, const atfft_complex b);
extern inline void atfft_multiply_by_and_swap_complex (atfft_complex *a, const atfft_complex b);

/******************************************
 * Other functions.
 ******************************************/
int atfft_next_power_of_2 (int x)
{
    if (x <= 0)
        return 0;
    else
        return pow (2, (int) log2 (x) + 1);
}

void atfft_twiddle_factor (int n,
                           int N,
                           enum atfft_direction d,
                           atfft_complex *t)
{
    atfft_scaled_twiddle_factor (n, N, d, 1.0, t);
}

void atfft_scaled_twiddle_factor (int n,
                                  int N,
                                  enum atfft_direction d,
                                  atfft_sample s,
                                  atfft_complex *t)
{
    atfft_sample x = 2.0 * n * M_PI / N;
    ATFFT_RE (*t) = cos (x) / s;
    ATFFT_IM (*t) = sin (x) / s;

    if (d == ATFFT_FORWARD)
        ATFFT_IM (*t) *= -1.0;
}

int atfft_is_prime (int x)
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

int atfft_mod (int a, int n)
{
    int r = a % n;
    return r < 0 ? r + n : r;
}

void atfft_gcd (int a, int b, int *gcd, int *x, int *y)
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

    if (x)
        *x = a < 0 ? -s0 : s0;

    if (y)
        *y = b < 0 ? -t0 : t0;
}

int atfft_mult_inverse_mod_n (int a, int n)
{
    int gcd, x, y;

    atfft_gcd (a % n, n, &gcd, &x, &y);

    if (gcd == 1)
        return atfft_mod (x, n);
    else
        return -1;
}

int atfft_prime_factors (int n, int *factors, int size)
{
    /* check array was provided */
    if (!factors)
        return -1;

    /* current factor */
    int f = 2;
    int n_factors = 0;
    int sqrt_size = (int) sqrt ((double) n);

    do
    {
        while (n % f)
        {
            f += f == 2 ? 1 : 2;

            /* a number will only have one prime factor greater than its square root */
            if (f > sqrt_size)
                f = n;
        }

        n /= f;
        factors [n_factors] = f;
        ++n_factors;
    }
    while (n > 1 && n_factors < size);

    return n_factors;
}

/******************************************
 * Functions for allocating sub-transform
 * plans for use on transforms of large
 * prime sizes.
 ******************************************/
static int atfft_integer_is_in_array (const int *arr, int size, int member)
{
    for (int i = 0; i < size; ++i)
    {
        if (arr [i] == member)
            return 1;
    }

    return 0;
}

static int atfft_sub_transform_sizes (const int *sizes,
                                      int *sub_transform_sizes,
                                      int size,
                                      int threshold)
{
    int n_sub_transforms = 0;

    for (int i = 0; i < size; ++i)
    {
        if (sizes [i] > threshold && 
            !atfft_integer_is_in_array (sub_transform_sizes, n_sub_transforms, sizes [i]))
            sub_transform_sizes [n_sub_transforms++] = sizes [i];
    }

    return n_sub_transforms;
}

static struct atfft_dft* atfft_populate_sub_transform (int size,
                                                       const int *sizes,
                                                       int n_sizes,
                                                       struct atfft_dft **size_sub_transforms,
                                                       enum atfft_direction direction,
                                                       enum atfft_format format)
{
    struct atfft_dft *fft = atfft_dft_create (size, direction, format);

    if (!fft)
        return NULL;

    for (int i = 0; i < n_sizes; ++i)
    {
        if (sizes [i] == size)
            size_sub_transforms [i] = fft;
    }

    return fft;
}

struct atfft_dft** atfft_init_sub_transforms (const int *sizes,
                                              int n_sizes,
                                              int *n_sub_transforms,
                                              struct atfft_dft **size_sub_transforms,
                                              enum atfft_direction direction,
                                              enum atfft_format format,
                                              int threshold)
{
    struct atfft_dft **sub_transforms = NULL;
    int *sub_transform_sizes = malloc (n_sizes * sizeof (*sub_transform_sizes));

    if (!sub_transform_sizes)
        return NULL;

    /* get unique sizes for which sub-transforms are required */
    *n_sub_transforms = atfft_sub_transform_sizes (sizes, sub_transform_sizes, n_sizes, threshold);

    if (*n_sub_transforms == 0)
        goto succeeded; /* we don't need to allocate any sub-transforms */

    /* allocate some space for sub-transforms */
    sub_transforms = calloc (*n_sub_transforms, sizeof (*sub_transforms));

    if (!sub_transforms)
        goto failed;

    /* create the sub-transform dft structs */
    for (int i = 0; i < *n_sub_transforms; ++i)
    {
        struct atfft_dft *sub_transform = atfft_populate_sub_transform (sub_transform_sizes [i],
                                                                        sizes,
                                                                        n_sizes,
                                                                        size_sub_transforms,
                                                                        direction,
                                                                        format);

        if (!sub_transform)
            goto failed;

        sub_transforms [i] = sub_transform;
    }

    goto succeeded;

failed:
    atfft_free_sub_transforms (sub_transforms, *n_sub_transforms);
    sub_transforms = NULL;

succeeded:
    free (sub_transform_sizes);
    return sub_transforms;
}

void atfft_free_sub_transforms (struct atfft_dft **sub_transforms,
                                int size)
{
    if (sub_transforms)
    {
        for (int i = 0; i < size; ++i)
        {
            atfft_dft_destroy (sub_transforms [i]);
        }

        free (sub_transforms);
    }
}
