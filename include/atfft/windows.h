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

/** @file 
 * Functions for generating windows.
 */

#ifndef ATFFT_WINDOWS_H_INCLUDED
#define ATFFT_WINDOWS_H_INCLUDED

#include <atfft/types.h>

#ifdef __cplusplus
extern "C"
{
#endif

/** An enum to represent the symmetry of a window function. */
enum atfft_window_symmetry
{
    ATFFT_SYMMETRIC, /**< Create a symmetric window. */
    ATFFT_PERIODIC /**< Create a periodic window. */
};

/**
 * Generate a Bartlett window
 *
 * @param window an array to generate the window in
 * @param size the length of the window
 * @param size the symmetry of the window
 */
void atfft_bartlett_window (atfft_sample *window, int size, enum atfft_window_symmetry symmetry);

/**
 * Generate a Hann window
 *
 * @param window an array to generate the window in
 * @param size the length of the window
 * @param size the symmetry of the window
 */
void atfft_hann_window (atfft_sample *window, int size, enum atfft_window_symmetry symmetry);

/**
 * Generate a Hamming window
 *
 * @param window an array to generate the window in
 * @param size the length of the window
 * @param size the symmetry of the window
 */
void atfft_hamming_window (atfft_sample *window, int size, enum atfft_window_symmetry symmetry);

/**
 * Generate a Blackman window
 *
 * @param window an array to generate the window in
 * @param size the length of the window
 * @param size the symmetry of the window
 */
void atfft_blackman_window (atfft_sample *window, int size, enum atfft_window_symmetry symmetry);

#ifdef __cplusplus
}
#endif

#endif /* ATFFT_WINDOWS_H_INCLUDED */
