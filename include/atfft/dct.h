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
 * struct and functions for performing discrete cosine transforms.
 */

#ifndef ATFFT_DCT_H_INCLUDED
#define ATFFT_DCT_H_INCLUDED

#include <atfft/types.h>

#ifdef __cplusplus
extern "C"
{
#endif

/** 
 * A Structure to hold internal FFT implementation.
 * 
 * When using atfft you will create one of these structures
 * using atfft_dct_create(), this structure is then passed 
 * to the calculation functions in order to compute DCTs.
 */
struct atfft_dct;

/**
 * Check whether a given signal length is supported by the current DCT implementation.
 *
 * @param size the signal length to check
 */
int atfft_dct_is_supported_size (int size);

/**
 * Create a dct structure.
 *
 * @param size the signal length the dct should operate on
 * @param direction the direction of the transform
 * @param format the type of transform (real or complex)
 */
struct atfft_dct* atfft_dct_create (int size, enum atfft_direction direction);

/**
 * Free a dct structure.
 *
 * @param dct the structure to free
 */
void atfft_dct_destroy (struct atfft_dct *dct);

/**
 * Perform a complex DCT.
 *
 * Performs a forward or inverse transform depending on what the dct
 * structure passed was created for.
 *
 * @param dct a valid dct structure 
 * @param in the input signal 
 *           (should have the number of samples the dct was created for)
 * @param out the output signal 
 *            (should have the number of samples the dct was created for)
 */
void atfft_dct_transform (struct atfft_dct *dct, const atfft_sample *in, atfft_sample *out);

#ifdef __cplusplus
}
#endif

#endif /* ATFFT_DCT_H_INCLUDED */
