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
