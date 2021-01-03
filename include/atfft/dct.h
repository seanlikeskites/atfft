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
