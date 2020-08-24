/*
 * Copyright (C) 2020 Sean Enderby <sean.enderby@gmail.com>
 *
 * This program is free software. It comes without any warranty, to
 * the extent permitted by applicable law. You can redistribute it 
 * and/or modify it under the terms of the Do What The Fuck You Want 
 * To Public License, Version 2, as published by Sam Hocevar. See 
 * the COPYING file for more details.
 */

#ifndef ATFFT_INTERNAL_H_INCLUDED
#define ATFFT_INTERNAL_H_INCLUDED

#include <atfft/atfft_shared.h>

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


#endif /* ATFFT_INTERNAL_H_INCLUDED */
