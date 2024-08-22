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

#ifndef ATFFT_DFT_COOLEY_TUKEY_H_INCLUDED
#define ATFFT_DFT_COOLEY_TUKEY_H_INCLUDED

#ifdef __cplusplus
extern "C"
{
#endif

#include <atfft/types.h>
#include "../cJSON/cJSON.h"

enum atfft_dft_ct_method
{
    ATFFT_DFT_CT_RECURSIVE,
    ATFFT_DFT_CT_ITERATIVE
};

struct atfft_dft_ct;

struct atfft_dft_ct* atfft_dft_ct_create (int size,
                                          enum atfft_direction direction,
                                          enum atfft_format format,
                                          enum atfft_dft_ct_method);

void atfft_dft_ct_destroy (void *fft);

void atfft_dft_ct_complex_transform (void *fft,
                                     atfft_complex *in,
                                     int in_stride,
                                     atfft_complex *out,
                                     int out_stride);

int atfft_dft_ct_is_fast_size (int size);

cJSON* atfft_dft_ct_get_plan (struct atfft_dft_ct *fft);

#ifdef __cplusplus
}
#endif

#endif /* ATFFT_DFT_COOLEY_TUKEY_H_INCLUDED */
