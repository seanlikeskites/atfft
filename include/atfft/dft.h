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
 * This header provides definitions the necessary functions for creating, executing, and freeing Discrete Fourier Transform
 * (DFT) plans.
 *
 * The procedure for performing DFTs with ATFFT is as follows:
 *  - Create a plan for the type of transform you want to perform.
 *  - Use this plan with the relevant transform functions to perform the transform on a signal.
 *  - Free the plan once you are done with it.
 *
 * A DFT plan is created using the atfft_dft_create() function. When creating a plan you must specify the length of signal
 * it will operate on as well as the direction of the transform and the format of the time domain signal. Depending on the
 * direction and format specified, the plan can then be used with the relevant transform function. The following table
 * gives the correct transform function to be used for each combination of direction and format.
 *
 * <table>
 * <tr><td style="border:none"> <th>ATFFT_FORWARD                      <th>ATFFT_BACKWARD
 * <tr><th>ATFFT_COMPLEX        <td>atfft_dft_complex_transform()      <td>atfft_dft_complex_transform()
 * <tr><th>ATFFT_REAL           <td>atfft_dft_real_forward_transform() <td>atfft_dft_real_backward_transform()
 * </table>
 *
 * The following code provides a simple example of how to perform a DFT on complex valued data using ATFFT.
 * \include complex_dft.c
 *
 * The following code provides a simple example of how to perform a DFT on real valued data using ATFFT.
 * \include real_dft.c
 */

#ifndef ATFFT_DFT_H_INCLUDED
#define ATFFT_DFT_H_INCLUDED

#include <atfft/types.h>

#ifdef __cplusplus
extern "C"
{
#endif

/**
 * @struct atfft_dft atfft/dft.h
 *
 * An opaque struct holing the details of a DFT plan.
 *
 * To perform a DFT you will create one of these structures using atfft_dft_create(), this structure is then passed to the
 * calculation functions (atfft_dft_complex_transform(), atfft_dft_real_forward_transform(), etc.) in order to compute DFTs.
 *
 * The exact contents of this structure differs depending on which FFT implementation is being used.
 */
struct atfft_dft;

/**
 * Create a DFT plan.
 *
 * @param size the signal length the DFT should operate on
 * @param direction the direction of the transform
 * @param format the type of time domain signal the transform will apply to (real or complex)
 */
struct atfft_dft* atfft_dft_create (int size, enum atfft_direction direction, enum atfft_format format);

/**
 * Free a DFT plan.
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
 * @param out_stride the stride to take when writing the output signal
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
 *            (should have size / 2 + 1 samples, where size is the
 *             size the fft was created for)
 */
void atfft_dft_real_forward_transform (struct atfft_dft *fft, const atfft_sample *in, atfft_complex *out);

/**
 * Perform a real forward DFT.
 *
 * @param fft a valid fft structure
 *            (should have been created with a direction of ATFFT_FORWARD
 *             and a format of ATFFT_REAL)
 * @param in the input signal
 *           (should have (size - 1) * in_stride + 1 samples, where size is the
 *            size the fft was created for)
 * @param in_stride the stride to take when reading the input signal
 * @param out the output signal
 *           (should have out_stride * size / 2 + 1 samples, where size is the
 *            size the fft was created for)
 * @param out_stride the stride to take when writing the output signal
 */
void atfft_dft_real_forward_transform_stride (struct atfft_dft *fft,
                                              const atfft_sample *in,
                                              int in_stride,
                                              atfft_complex *out,
                                              int out_stride);

/**
 * Perform a real inverse DFT.
 *
 * @param fft a valid fft structure
 *            (should have been created with a direction of ATFFT_BACKWARD
 *             and a format of ATFFT_REAL)
 * @param in the input signal
 *           (should have size / 2 + 1 samples, where size is the
 *            size the fft was created for)
 * @param out the output signal
 *            (should have the number of samples the fft was created for)
 */
void atfft_dft_real_backward_transform (struct atfft_dft *fft, atfft_complex *in, atfft_sample *out);

/**
 * Perform a real inverse DFT.
 *
 * @param fft a valid fft structure
 *            (should have been created with a direction of ATFFT_BACKWARD
 *             and a format of ATFFT_REAL)
 * @param in the input signal
 *           (should have in_stride * size / 2 + 1 samples, where size is the
 *            size the fft was created for)
 * @param in_stride the stride to take when reading the input signal
 * @param out the output signal
 *           (should have (size - 1) * out_stride + 1 samples, where size is the
 *            size the fft was created for)
 * @param out_stride the stride to take when writing the output signal
 */
void atfft_dft_real_backward_transform_stride (struct atfft_dft *fft,
                                               atfft_complex *in,
                                               int in_stride,
                                               atfft_sample *out,
                                               int out_stride);

#ifdef __cplusplus
}
#endif

#endif /* ATFFT_DFT_H_INCLUDED */
