/*
 * Copyright (C) 2016 Sean Enderby <sean.enderby@gmail.com>
 *
 * This program is free software. It comes without any warranty, to
 * the extent permitted by applicable law. You can redistribute it 
 * and/or modify it under the terms of the Do What The Fuck You Want 
 * To Public License, Version 2, as published by Sam Hocevar. See 
 * the COPYING file for more details.
 */

#include <math.h>
#include <atfft/atfft_windows.h>

void atfft_bartlett_window (atfft_sample *window, int size, enum atfft_window_symmetry symmetry)
{
    int i = 0;
    atfft_sample center;

    if (symmetry == ATFFT_SYMMETRIC)
        center = (size - 1.0) / 2.0;
    else
        center = size / 2.0;

    for (i = 0; i < size; ++i)
    {
        window [i] = 1.0 - fabs ((i - center) / center);
    }
}

void atfft_hann_window (atfft_sample *window, int size, enum atfft_window_symmetry symmetry)
{
    int i = 0;
    atfft_sample den;

    if (symmetry == ATFFT_SYMMETRIC)
        den = size - 1.0;
    else
        den = size;

    for (i = 0; i < size; ++i)
    {
        atfft_sample x = sin (M_PI * i / den);
        window [i] = x * x;
    }
}

void atfft_hamming_window (atfft_sample *window, int size, enum atfft_window_symmetry symmetry)
{
    int i = 0;
    atfft_sample den;

    if (symmetry == ATFFT_SYMMETRIC)
        den = size - 1.0;
    else
        den = size;

    for (i = 0; i < size; ++i)
    {
        atfft_sample x = cos (2 * M_PI * i / den);
        window [i] = 0.54 - 0.46 * x;
    }
}

void atfft_blackman_window (atfft_sample *window, int size, enum atfft_window_symmetry symmetry)
{
    int i = 0;
    atfft_sample den;

    if (symmetry == ATFFT_SYMMETRIC)
        den = size - 1.0;
    else
        den = size;

    for (i = 0; i < size; ++i)
    {
        atfft_sample x = 2 * M_PI * i / den;
        window [i] = 0.42 - 0.5 * cos (x) + 0.08 * cos (2.0 * x);
    }
}
