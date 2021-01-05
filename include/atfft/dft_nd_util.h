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
 * struct and functions for performing discrete fourier transforms.
 */

#ifndef ATFFT_DFT_ND_UTIL_H_INCLUDED
#define ATFFT_DFT_ND_UTIL_H_INCLUDED

#ifdef __cplusplus
extern "C"
{
#endif

/**
 * Take the product of an array of integers.
 */
int atfft_int_array_product (const int *array, int size);

/**
 * Return the size of the complex output when performing am
 * n-dimensional DFT on a real values signal.
 */
int atfft_nd_halfcomplex_size (const int *dims, int n_dims);

#ifdef __cplusplus
}
#endif

#endif /* ATFFT_DFT_ND_UTIL_H_INCLUDED */
