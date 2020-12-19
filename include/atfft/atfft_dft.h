/*
 * Copyright (C) 2016 Sean Enderby <sean.enderby@gmail.com>
 *
 * This program is free software. It comes without any warranty, to
 * the extent permitted by applicable law. You can redistribute it 
 * and/or modify it under the terms of the Do What The Fuck You Want 
 * To Public License, Version 2, as published by Sam Hocevar. See 
 * the COPYING file for more details.
 */

#ifndef ATFFT_DFT_H_INCLUDED
#define ATFFT_DFT_H_INCLUDED

#include <atfft/atfft_shared.h>

#ifdef __cplusplus
extern "C"
{
#endif

/** 
 * A Structure to hold internal FFT implementation.
 * 
 * When using atfft you will create one of these structures
 * using atfft_create(), this structure is then passed 
 * to the calculation functions in order to compute DFTs.
 */
struct atfft_dft;

/**
 * Create an fft structure.
 *
 * @param size the signal length the fft should operate on
 * @param direction the direction of the transform
 * @param format the type of transform (real or complex)
 */
struct atfft_dft* atfft_dft_create (int size, enum atfft_direction direction, enum atfft_format format);

/**
 * Free an fft structure.
 *
 * @param fft the structure to free
 */
void atfft_dft_destroy (struct atfft_dft *fft);

/**
 * Perform a complex DFT.
 *
 * Performs a forward or inverse transform depending on what the fft
 * structure passed was created for.
 *
 * @param fft a valid fft structure 
 *            (should have been created with a format of ATFFT_COMPLEX)
 * @param in the input signal 
 *           (should have the number of samples the fft was created for)
 * @param out the output signal 
 *            (should have the number of samples the fft was created for)
 */
void atfft_dft_complex_transform (struct atfft_dft *fft,
                                  atfft_complex *in,
                                  atfft_complex *out);

/**
 * Perform a complex DFT with a different stride for input and output.
 *
 * Performs a forward or inverse transform depending on what the fft
 * structure passed was created for.
 *
 * @param fft a valid fft structure 
 *            (should have been created with a format of ATFFT_COMPLEX)
 * @param in the input signal 
 *           (should have the number of samples the fft was created for)
 * @param in_stride the stride to take when reading the input signal
 * @param out the output signal 
 *            (should have the number of samples the fft was created for)
 * @param in_stride the stride to take when writing the output signal
 */
void atfft_dft_complex_transform_stride (struct atfft_dft *fft,
                                         atfft_complex *in,
                                         int in_stride,
                                         atfft_complex *out,
                                         int out_stride);

/**
 * Perform a real forward DFT.
 *
 * @param fft a valid fft structure 
 *            (should have been created with a direction of ATFFT_FORWARD
 *             and a format of ATFFT_REAL)
 * @param in the input signal 
 *           (should have the number of samples the fft was created for)
 * @param out the output signal 
 *            (should have size / 2 + 1 samples, where size if the
 *             size the fft was created for)
 */
void atfft_dft_real_forward_transform (struct atfft_dft *fft, const atfft_sample *in, atfft_complex *out);

/**
 * Perform a real inverse DFT.
 *
 * @param fft a valid fft structure 
 *            (should have been created with a direction of ATFFT_BACKWARD
 *             and a format of ATFFT_REAL)
 * @param in the input signal 
 *           (should have size / 2 + 1 samples, where size if the
 *            size the fft was created for)
 * @param out the output signal 
 *            (should have the number of samples the fft was created for)
 */
void atfft_dft_real_backward_transform (struct atfft_dft *fft, atfft_complex *in, atfft_sample *out);

#ifdef __cplusplus
}
#endif

#endif /* ATFFT_DFT_H_INCLUDED */
