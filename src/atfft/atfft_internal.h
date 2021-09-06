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

#ifndef ATFFT_INTERNAL_H_INCLUDED
#define ATFFT_INTERNAL_H_INCLUDED

#include <atfft/types.h>

/**
 * Copy the complex value x into the complex value y.
 */
inline void atfft_copy_complex (const atfft_complex x, atfft_complex *y)
{
    ATFFT_RE (*y) = ATFFT_RE (x);
    ATFFT_IM (*y) = ATFFT_IM (x);
}

/**
 * Swap the real and imaginary part of a complex number: y = j * conj (x)
 */
inline void atfft_swap_complex (const atfft_complex x, atfft_complex *y)
{
    ATFFT_RE (*y) = ATFFT_IM (x);
    ATFFT_IM (*y) = ATFFT_RE (x);
}

/**
 * Sum two complex numbers: s = a + b
 */
inline void atfft_sum_complex (const atfft_complex a,
                               const atfft_complex b,
                               atfft_complex *s)
{
    ATFFT_RE (*s) = ATFFT_RE (a) + ATFFT_RE (b);
    ATFFT_IM (*s) = ATFFT_IM (a) + ATFFT_IM (b);
}

/**
 * Take the difference of two complex numbers: d = a - b
 */
inline void atfft_difference_complex (const atfft_complex a,
                                      const atfft_complex b,
                                      atfft_complex *d)
{
    ATFFT_RE (*d) = ATFFT_RE (a) - ATFFT_RE (b);
    ATFFT_IM (*d) = ATFFT_IM (a) - ATFFT_IM (b);
}

/**
 * Find the product of two complex numbers: p = ab
 */
inline void atfft_product_complex (const atfft_complex a,
                                   const atfft_complex b,
                                   atfft_complex *p)
{
    ATFFT_RE (*p) = ATFFT_RE (a) * ATFFT_RE (b) -
                    ATFFT_IM (a) * ATFFT_IM (b);
    ATFFT_IM (*p) = ATFFT_RE (a) * ATFFT_IM (b) +
                    ATFFT_IM (a) * ATFFT_RE (b);
}

/**
 * Multiply a complex variable by another: a *= b
 */
inline void atfft_multiply_by_complex (atfft_complex *a, const atfft_complex b)
{
    atfft_sample temp;
    temp = ATFFT_RE (*a);

    ATFFT_RE (*a) = temp * ATFFT_RE (b) -
                    ATFFT_IM (*a) * ATFFT_IM (b);
    ATFFT_IM (*a) = temp * ATFFT_IM (b) +
                    ATFFT_IM (*a) * ATFFT_RE (b);
}

/**
 * Multiply a complex variable by another and swap real and imaginary parts:
 *
 * a = j * conj(ab)
 */
inline void atfft_multiply_by_and_swap_complex (atfft_complex *a, const atfft_complex b)
{
    atfft_sample temp;
    temp = ATFFT_RE (*a);

    ATFFT_RE (*a) = temp * ATFFT_IM (b) +
                    ATFFT_IM (*a) * ATFFT_RE (b);
    ATFFT_IM (*a) = temp * ATFFT_RE (b) -
                    ATFFT_IM (*a) * ATFFT_IM (b);
}

/**
 * Find the next power of 2 which is greater than x.
 */
int atfft_next_power_of_2 (int x);

/**
 * Compute a twiddle factor.
 *
 * When d == ATFFT_FORWARD:
 *     t = e^(-2*pi*j*n/N)
 *
 * When d == ATFFT_BACKWARD:
 *     t = e^(2*pi*j*n/N)
 */
void atfft_twiddle_factor (int n,
                           int N,
                           enum atfft_direction d,
                           atfft_complex *t);

/**
 * Compute a scaled twiddle factor.
 *
 * When d == ATFFT_FORWARD:
 *     t = e^(-2*pi*j*n/N) / s
 *
 * When d == ATFFT_BACKWARD:
 *     t = e^(2*pi*j*n/N) / s
 */
void atfft_scaled_twiddle_factor (int n,
                                  int N,
                                  enum atfft_direction d,
                                  atfft_sample s,
                                  atfft_complex *t);

/**
 * Return 1 if x is prime, 0 otherwise.
 */
int atfft_is_prime (int x);

/**
 * Return a mod n.
 */
int atfft_mod (int a, int n);

/**
 * Put the greatest common divisor of a and b in gcd.
 *
 * If x and y are not null they are set to values such that:
 *
 *   ax + by = gcd
 */
void atfft_gcd (int a, int b, int *gcd, int *x, int *y);

/**
 * Return the multiplicative inverse of a mod n.
 */
int atfft_mult_inverse_mod_n (int a, int n);

/******************************************
 * Functions for allocating sub-transform
 * plans for use on transforms of large
 * prime sizes.
 ******************************************/
struct atfft_dft** atfft_init_sub_transforms (const int *sizes,
                                              int n_sizes,
                                              int *n_sub_transforms,
                                              struct atfft_dft **size_sub_transforms,
                                              enum atfft_direction direction,
                                              enum atfft_format format,
                                              int threshold);

void atfft_free_sub_transforms (struct atfft_dft **sub_transforms,
                                int size);

#endif /* ATFFT_INTERNAL_H_INCLUDED */
