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
 * @param plan the plan to free
 */
void atfft_dft_destroy (struct atfft_dft *plan);

/**
 * Perform a DFT on complex data.
 *
 * Performs a forward or inverse transform depending on what the @p plan was created for.
 *
 * @param plan a valid DFT plan 
 *             (should have been created with a format of ATFFT_COMPLEX)
 * @param in the input signal
 *           (should contain at least as many elements as the signal size the @p plan was created for)
 * @param out the output signal
 *           (should contain at least as many elements as the signal size the @p plan was created for)
 */
void atfft_dft_complex_transform (struct atfft_dft *plan,
                                  atfft_complex *in,
                                  atfft_complex *out);

/**
 * Perform a DFT on complex data, taking strides of different length for input and output.
 *
 * Performs a forward or inverse transform depending on what the @p plan was created for.
 *
 * @param plan a valid DFT plan 
 *             (should have been created with a format of ATFFT_COMPLEX)
 * @param in the input signal
 *            (should contain at least <b>(@p in_stride * (size - 1) + 1)</b> elements,
 *             where size is the signal size the @p plan was created for)
 * @param in_stride the stride to take when reading the input
 * @param out the output signal
 *            (should contain at least <b>(@p out_stride * (size - 1) + 1)</b> elements,
 *             where size is the signal size the @p plan was created for)
 * @param out_stride the stride to take when writing the output
 */
void atfft_dft_complex_transform_stride (struct atfft_dft *plan,
                                         atfft_complex *in,
                                         int in_stride,
                                         atfft_complex *out,
                                         int out_stride);

/**
 * Perform a forward DFT on real data.
 *
 * @param plan a valid DFT plan 
 *            (should have been created with a direction of ATFFT_FORWARD and a format of ATFFT_REAL)
 * @param in the input signal
 *           (should contain at least as many elements as the signal size the @p plan was created for)
 * @param out the output signal
 *            (should contain at least <b>(\ref atfft_halfcomplex_size (size))</b> elements,
 *             where size is the signal size the @p plan was created for)
 */
void atfft_dft_real_forward_transform (struct atfft_dft *plan, const atfft_sample *in, atfft_complex *out);

/**
 * Perform a forward DFT on real data, taking strides of different length for input and output.
 *
 * @param plan a valid DFT plan 
 *            (should have been created with a direction of ATFFT_FORWARD and a format of ATFFT_REAL)
 * @param in the input signal
 *            (should contain at least <b>(@p in_stride * (size - 1) + 1)</b> elements,
 *             where size is the signal size the @p plan was created for)
 * @param in_stride the stride to take when reading the input
 * @param out the output signal
 *            (should contain at least <b>(@p out_stride * (\ref atfft_halfcomplex_size (size) - 1) + 1)</b> elements,
 *             where size is the signal size the @p plan was created for)
 * @param out_stride the stride to take when writing the output
 */
void atfft_dft_real_forward_transform_stride (struct atfft_dft *plan,
                                              const atfft_sample *in,
                                              int in_stride,
                                              atfft_complex *out,
                                              int out_stride);

/**
 * Perform a backward DFT on halfcomplex data.
 *
 * @param plan a valid DFT plan 
 *            (should have been created with a direction of ATFFT_BACKWARD and a format of ATFFT_REAL)
 * @param in the input signal
 *            (should contain at least <b>(\ref atfft_halfcomplex_size (size))</b> elements,
 *             where size is the signal size the @p plan was created for)
 * @param out the output signal
  *           (should contain at least as many elements as the signal size the @p plan was created for)
 */
void atfft_dft_real_backward_transform (struct atfft_dft *plan, atfft_complex *in, atfft_sample *out);

/**
 * Perform a backward DFT on halfcomplex data, taking strides of different length for input and output.
 *
 * @param plan a valid DFT plan 
 *            (should have been created with a direction of ATFFT_BACKWARD and a format of ATFFT_REAL)
 * @param in the input signal
 *           (should contain at least <b>(in_stride * (\ref atfft_halfcomplex_size (size) - 1) + 1)</b> elements,
 *            where size is the signal size the @p plan was created for)
 * @param in_stride the stride to take when reading the input
 * @param out the output signal
 *            (should contain at least <b>(@p out_stride * (size - 1) + 1)</b> elements,
 *             where size is the signal size the @p plan was created for)
 * @param out_stride the stride to take when writing the output
 */
void atfft_dft_real_backward_transform_stride (struct atfft_dft *plan,
                                               atfft_complex *in,
                                               int in_stride,
                                               atfft_sample *out,
                                               int out_stride);

#ifdef __cplusplus
}
#endif

#endif /* ATFFT_DFT_H_INCLUDED */
