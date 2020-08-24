/*
 * Copyright (C) 2020 Sean Enderby <sean.enderby@gmail.com>
 *
 * This program is free software. It comes without any warranty, to
 * the extent permitted by applicable law. You can redistribute it 
 * and/or modify it under the terms of the Do What The Fuck You Want 
 * To Public License, Version 2, as published by Sam Hocevar. See 
 * the COPYING file for more details.
 */

#include "atfft_internal.h"
#include <math.h>

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
