/*
 * Copyright (C) 2020 Sean Enderby <sean.enderby@gmail.com>
 *
 * This program is free software. It comes without any warranty, to
 * the extent permitted by applicable law. You can redistribute it 
 * and/or modify it under the terms of the Do What The Fuck You Want 
 * To Public License, Version 2, as published by Sam Hocevar. See 
 * the COPYING file for more details.
 */

#ifndef ATFFT_DFT_BLUESTEIN_H_INCLUDED
#define ATFFT_DFT_BLUESTEIN_H_INCLUDED

#ifdef __cplusplus
extern "C"
{
#endif

#include <atfft/atfft_shared.h>

struct atfft_dft_bluestein;

struct atfft_dft_bluestein* atfft_dft_bluestein_create (int size,
                                                        enum atfft_direction direction,
                                                        enum atfft_format format);

void atfft_dft_bluestein_destroy (struct atfft_dft_bluestein *fft);

void atfft_dft_bluestein_complex_transform (struct atfft_dft_bluestein *fft,
                                            atfft_complex *in,
                                            int in_stride,
                                            atfft_complex *out,
                                            int out_stride);

#ifdef __cplusplus
}
#endif

#endif /* ATFFT_DFT_BLUESTEIN_H_INCLUDED */
