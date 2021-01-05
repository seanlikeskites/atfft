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
 * struct and functions for performing n-dimensional discrete fourier transforms.
 */

#ifndef ATFFT_DFT_ND_H_INCLUDED
#define ATFFT_DFT_ND_H_INCLUDED

#include <atfft/types.h>

#ifdef __cplusplus
extern "C"
{
#endif

/**
 * A Structure to hold internal n-dimensional FFT implementation.
 */
struct atfft_dft_nd;

/**
 * Create an n-dimensional fft structure.
 *
 * @param dims the length of each dimension
 * @param n_dims the number of dimensions
 * @param direction the direction of the transform
 * @param format the type of transform (real or complex)
 */
struct atfft_dft_nd* atfft_dft_nd_create (const int *dims,
                                          int n_dims,
                                          enum atfft_direction direction,
                                          enum atfft_format format);

/**
 * Free an n-dimensional fft structure.
 *
 * @param fft the structure to free
 */
void atfft_dft_nd_destroy (struct atfft_dft_nd *fft);

/**
 * Perform an n-dimensional complex DFT.
 *
 * Performs a forward or inverse transform depending on what the fft
 * structure passed was created for.
 *
 * @param fft a valid fft structure
 *            (should have been created with a format of ATFFT_COMPLEX)
 * @param in the input signal
 *           (should have the size and dimensionality the fft was created for)
 * @param out the output signal
 *           (should have the size and dimensionality the fft was created for)
 */
void atfft_dft_nd_complex_transform (struct atfft_dft_nd *fft, atfft_complex *in, atfft_complex *out);

/**
 * Perform an n-dimensional real forward DFT.
 *
 * @param fft a valid fft structure
 *            (should have been created with a format of ATFFT_COMPLEX)
 * @param in the input signal
 *           (should have the size and dimensionality the fft was created for)
 * @param out the output signal
 *           (should have the size and dimensionality the fft was created for)
 */
void atfft_dft_nd_real_forward_transform (struct atfft_dft_nd *fft, const atfft_sample *in, atfft_complex *out);

/**
 * Perform an n-dimensional real backward DFT.
 *
 * @param fft a valid fft structure
 *            (should have been created with a format of ATFFT_COMPLEX)
 * @param in the input signal
 *           (should have the size and dimensionality the fft was created for)
 * @param out the output signal
 *           (should have the size and dimensionality the fft was created for)
 */
void atfft_dft_nd_real_backward_transform (struct atfft_dft_nd *fft, atfft_complex *in, atfft_sample *out);

#ifdef __cplusplus
}
#endif

#endif /* ATFFT_DFT_ND_H_INCLUDED */
