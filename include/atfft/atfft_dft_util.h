/*
 * Copyright (C) 2020 Sean Enderby <sean.enderby@gmail.com>
 *
 * This program is free software. It comes without any warranty, to
 * the extent permitted by applicable law. You can redistribute it
 * and/or modify it under the terms of the Do What The Fuck You Want
 * To Public License, Version 2, as published by Sam Hocevar. See
 * the COPYING file for more details.
 */

/** @file
 * struct and functions for performing discrete fourier transforms.
 */

#ifndef ATFFT_DFT_UTIL_H_INCLUDED
#define ATFFT_DFT_UTIL_H_INCLUDED

#include <atfft/atfft_shared.h>

#ifdef __cplusplus
extern "C"
{
#endif

/**
 * Return the size of the complex output when performing a DFT on
 * a real valued signal.
 */
int atfft_dft_halfcomplex_size (int size);

/**
 * Create a complex signal from a halfcomplex signal.
 *
 * @param in the input signal (should contain size / 2 + 1 elements)
 * @param out the output signal (should contain size elements)
 * @param size the length of the output signal
 */
void atfft_halfcomplex_to_complex (atfft_complex *in, atfft_complex *out, int size);
void atfft_halfcomplex_to_complex_stride (atfft_complex *in,
                                          int in_stride,
                                          atfft_complex *out,
                                          int out_stride,
                                          int size);

void atfft_complex_to_halfcomplex (atfft_complex *in, atfft_complex *out, int size);
void atfft_complex_to_halfcomplex_stride (atfft_complex *in,
                                          int in_stride,
                                          atfft_complex *out,
                                          int out_stride,
                                          int size);

#ifdef __cplusplus
}
#endif

#endif /* ATFFT_DFT_UTIL_H_INCLUDED */
