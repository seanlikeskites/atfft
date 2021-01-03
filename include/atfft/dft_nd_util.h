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
