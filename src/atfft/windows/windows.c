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

#include <math.h>
#include <atfft/windows.h>
#include "constants.h"

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
        atfft_sample x = cos (2.0 * M_PI * i / den);
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
        atfft_sample x = 2.0 * M_PI * i / den;
        window [i] = 0.42 - 0.5 * cos (x) + 0.08 * cos (2.0 * x);
    }
}
