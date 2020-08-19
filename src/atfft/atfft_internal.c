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
    ATFFT_REAL (*t) = cos (x) / s;
    ATFFT_IMAG (*t) = sin (x) / s;

    if (d == ATFFT_FORWARD)
        ATFFT_IMAG (*t) *= -1.0;
}
