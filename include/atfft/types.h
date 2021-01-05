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
 * Definitions of the types used by ATFFT and some
 * useful macros and functions for working with them.
 */

#ifndef ATFFT_SHARED_H_INCLUDED
#define ATFFT_SHARED_H_INCLUDED

#ifdef __cplusplus
extern "C"
{
#endif

/** An enum to represent the direction of a transform. */
enum atfft_direction
{
    ATFFT_FORWARD, /**< Create a plan for a forward transform (from time to frequency domain). */
    ATFFT_BACKWARD /**< Create a plan for a backward transform (from frequency to time domain). */
};

/** An enum to represent the type of values a transform will operate on. */
enum atfft_format
{
    ATFFT_COMPLEX, /**< Create a plan for operating on complex valued signals. */
    ATFFT_REAL /**< Create a plan for operating on real valued signals. */
};

/** A complex float type. */
typedef float atfft_complex_f [2];
/** A complex double type. */
typedef double atfft_complex_d [2];
/** A complex long double type. */
typedef long double atfft_complex_l [2];

/** A macro to access the real part of a complex value. */
#define ATFFT_RE(x) ((x) [0])
/** A macro to access the imaginary part of a complex value. */
#define ATFFT_IM(x) ((x) [1])

/** A macro to return the smallest of two values. */
#define ATFFT_MIN(x, y) ((x) < (y) ? (x) : (y))
/** A macro to return the largest of two values. */
#define ATFFT_MAX(x, y) ((x) > (y) ? (x) : (y))

/* typedefs for the types used in ATFFT */
#if defined(ATFFT_TYPE_FLOAT)
    typedef float atfft_sample;
    typedef atfft_complex_f atfft_complex;

#elif defined(ATFFT_TYPE_LONG_DOUBLE)
    typedef long double atfft_sample;
    typedef atfft_complex_l atfft_complex;

#else
    /** @cond */
#   if !defined(ATFFT_TYPE_DOUBLE)
#       define ATFFT_TYPE_DOUBLE
#   endif
    /** @endcond */
    /** A typedef defining the type used for real values in ATFFT.
     *
     * By default this will be double, but by providing additional definitions at compile time this can be changed.
     *
     * <table>
     * <tr><th>Definition                   <th>Sample Type
     * <tr><td>-D ATFFT_TYPE_FLOAT          <td>float
     * <tr><td>-D ATFFT_TYPE_LONG_DOUBLE    <td>long double
     * <tr><td>-D ATFFT_TYPE_DOUBLE         <td>double
     * </table>
     *
     * These definitions take precedence in the order given in the above table, i.e. if ATFFT_TYPE_FLOAT is defined, \ref
     * atfft_sample will always be typedefed to float, even if the other two definitions are present.
     *
     * In order for \ref atfft_sample to be double, either none of the above definitions must be present or only
     * ATFFT_TYPE_DOUBLE must be defined.
     */
    typedef double atfft_sample;
    /** A typedef defining the type used for complex values in ATFFT.
     *
     * Complex values are stored as an array of two values (a real and an imaginary part). Typedefs providing a complex type
     * for each precision of floating point value are described above (\ref atfft_complex_f, \ref atfft_complex_d, \ref
     * atfft_complex_l). \ref atfft_complex is typedefed to one of those three, such that its real and imaginary part have
     * the same type as \ref atfft_sample.
     *
     * By default atfft_complex will be an array of two doubles (\ref atfft_complex_d), but by providing additional
     * definitions at compile time this can be changed.
     *
     * <table>
     * <tr><th>Definition                   <th>Complex Type
     * <tr><td>-D ATFFT_TYPE_FLOAT          <td>\ref atfft_complex_f
     * <tr><td>-D ATFFT_TYPE_LONG_DOUBLE    <td>\ref atfft_complex_l
     * <tr><td>-D ATFFT_TYPE_DOUBLE         <td>\ref atfft_complex_d
     * </table>
     *
     * These definitions take precedence in the order given in the above table, i.e. if ATFFT_TYPE_FLOAT is defined,
     * \ref atfft_complex will always be typedefed to \ref atfft_complex_f, even if the other two definitions are present.
     *
     * In order for \ref atfft_complex to be \ref atfft_complex_d, either none of the above definitions must be present or
     * only ATFFT_TYPE_DOUBLE must be defined.
     */
    typedef atfft_complex_d atfft_complex;
#endif

/**
 * Check if an integer is even.
 *
 * @param x integer to check
 *
 * @return 1 if @p x is even, 0 otherwise
 */
int atfft_is_even (unsigned int x);

/**
 * Check if an integer is odd.
 *
 * @param x integer to check
 *
 * @return 1 if @p x is odd, 0 otherwise
 */
int atfft_is_odd (unsigned int x);

/**
 * Check if an integer is a power of 2
 *
 * @param x integer to check
 *
 * @return 1 if @p x is a power of 2, 0 otherwise
 */
int atfft_is_power_of_2 (unsigned int x);

/**
 * Scale a real valued signal.
 *
 * Applies a multiplicative scaling to a real valued signal.
 *
 * @param data the signal to scale (should contain at least @p size elements)
 * @param size the length of the signal
 * @param scale_factor the factor to scale by
 */
void atfft_scale_real (atfft_sample *data, int size, atfft_sample scale_factor);

/**
 * Normalise a real DFT output.
 *
 * Applies 1 / @p size scaling to a real valued signal.
 *
 * @param data the signal to normalise (should contain at least @p size elements)
 * @param size the length of the signal
 */
void atfft_normalise_real (atfft_sample *data, int size);

/**
 * Scale a complex valued signal.
 *
 * Applies a multiplicative scaling to a complex valued signal.
 *
 * @param data the signal to scale (should contain at least @p size elements)
 * @param size the length of the signal
 * @param scale_factor the factor to scale by
 */
void atfft_scale_complex (atfft_complex *data, int size, atfft_sample scale_factor);

/**
 * Normalise a complex DFT output.
 *
 * Applies 1 / @p size scaling to a complex valued signal.
 *
 * @param data the signal to normalise (should contain at least @p size elements)
 * @param size the length of the signal
 */
void atfft_normalise_complex (atfft_complex *data, int size);

/**
 * Find the magnitude of a complex value.
 *
 * @param x value to find the magnitude of
 *
 * @return The magnitude of @p x.
 *
 * @return Where @p x is some complex value \f$ a + bi \f$, atfft_abs() returns \f$ \sqrt{a^{2} + b^{2}} \f$.
 */
atfft_sample atfft_abs (const atfft_complex x);

/**
 * Find the argument (phase) of a complex value.
 *
 * @param x value to find the argument of
 *
 * @return The argument of @p x in the range \f$ -\pi \f$ to \f$ \pi \f$.
 *
 * @return Where @p x is some complex value \f$ a + bi \f$, atfft_arg() returns \f$ \mathop{\rm atan2}(b, a) \f$.
 */
atfft_sample atfft_arg (const atfft_complex x);

/**
 * Copy the real part of a complex valued signal into the output.
 *
 * @param in the signal to be copied from (should contain at least @p size elements)
 * @param out the signal to copy to (should contain at least @p size elements)
 * @param size the number of values to copy
 */
void atfft_real (atfft_complex *in, atfft_sample *out, int size);

/**
 * Copy the real part of a complex valued signal into the output, taking strides of different length for input and output.
 *
 * @param in the signal to be copied from
 *        (should contain at least (@p size - 1) * @p in_stride + 1 elements)
 * @param in_stride the stride to take when reading the input
 * @param out the signal to copy to
 *        (should contain at least (@p size - 1) * @p out_stride + 1 elements)
 * @param out_stride the stride to take when writing the output
 * @param size the number of values to copy
 */
void atfft_real_stride (atfft_complex *in,
                        int in_stride,
                        atfft_sample *out,
                        int out_stride,
                        int size);

/**
 * Copy the imaginary part of a complex valued signal into the output.
 *
 * @param in the signal to be copied from (should contain at least @p size elements)
 * @param out the signal to copy to (should contain at least @p size elements)
 * @param size the number of values to copy
 */
void atfft_imag (atfft_complex *in, atfft_sample *out, int size);

/**
 * Copy the imaginary part of a complex valued signal into the output, taking strides of different length for input and
 * output.
 *
 * @param in the signal to be copied from
 *        (should contain at least (@p size - 1) * @p in_stride + 1 elements)
 * @param in_stride the stride to take when reading the input
 * @param out the signal to copy to
 *        (should contain at least (@p size - 1) * @p out_stride + 1 elements)
 * @param out_stride the stride to take when writing the output
 * @param size the number of values to copy
 */
void atfft_imag_stride (atfft_complex *in,
                        int in_stride,
                        atfft_sample *out,
                        int out_stride,
                        int size);

/**
 * Create a complex signal from a real signal, setting all imaginary parts to 0.
 *
 * @param in the input signal (should contain at least @p size elements)
 * @param out the output signal (should contain at least @p size elements)
 * @param size the length of the signal
 */
void atfft_real_to_complex (const atfft_sample *in, atfft_complex *out, int size);

/**
 * Create a complex signal from a real signal, setting all imaginary parts to 0 and taking strides of different length for
 * input and output.
 *
 * @param in the input signal
 *        (should contain at least (@p size - 1) * @p in_stride + 1 elements)
 * @param in_stride the stride to take when reading the input
 * @param out the output signal
 *        (should contain at least (@p size - 1) * @p out_stride + 1 elements)
 * @param out_stride the stride to take when writing the output
 * @param size the length of the signal
 */
void atfft_real_to_complex_stride (const atfft_sample *in,
                                   int in_stride,
                                   atfft_complex *out,
                                   int out_stride,
                                   int size);

/**
 * Convert a real valued signal from single precision floats to \ref atfft_sample values.
 *
 * @param in the input signal (should contain at least @p size elements)
 * @param out the output signal (should contain at least @p size elements)
 * @param size the length of the signal
 */
void atfft_float_to_sample_real (const float *in, atfft_sample *out, int size);

/**
 * Convert a real valued signal from single precision floats to \ref atfft_sample values, taking strides of different length
 * for input and output.
 *
 * @param in the input signal
 *        (should contain at least (@p size - 1) * @p in_stride + 1 elements)
 * @param in_stride the stride to take when reading the input
 * @param out the output signal
 *        (should contain at least (@p size - 1) * @p out_stride + 1 elements)
 * @param out_stride the stride to take when writing the output
 * @param size the length of the signal
 */
void atfft_float_to_sample_real_stride (const float *in,
                                        int in_stride,
                                        atfft_sample *out,
                                        int out_stride,
                                        int size);

/**
 * Convert a real valued signal from \ref atfft_sample values to single precision floats.
 *
 * @param in the input signal (should contain at least @p size elements)
 * @param out the output signal (should contain at least @p size elements)
 * @param size the length of the signal
 */
void atfft_sample_to_float_real (const atfft_sample *in, float *out, int size);

/**
 * Convert a real valued signal from \ref atfft_sample values to single precision floats, taking strides of different length
 * for input and output.
 *
 * @param in the input signal
 *        (should contain at least (@p size - 1) * @p in_stride + 1 elements)
 * @param in_stride the stride to take when reading the input
 * @param out the output signal
 *        (should contain at least (@p size - 1) * @p out_stride + 1 elements)
 * @param out_stride the stride to take when writing the output
 * @param size the length of the signal
 */
void atfft_sample_to_float_real_stride (const atfft_sample *in,
                                        int in_stride,
                                        float *out,
                                        int out_stride,
                                        int size);

/**
 * Convert a complex valued signal from single precision floats to \ref atfft_sample values.
 *
 * @param in the input signal (should contain at least @p size elements)
 * @param out the output signal (should contain at least @p size elements)
 * @param size the length of the signal
 */
void atfft_float_to_sample_complex (atfft_complex_f *in, atfft_complex *out, int size);

/**
 * Convert a complex valued signal from single precision floats to \ref atfft_sample values, taking strides of different
 * length for input and output.
 *
 * @param in the input signal
 *        (should contain at least (@p size - 1) * @p in_stride + 1 elements)
 * @param in_stride the stride to take when reading the input
 * @param out the output signal
 *        (should contain at least (@p size - 1) * @p out_stride + 1 elements)
 * @param out_stride the stride to take when writing the output
 * @param size the length of the signal
 */
void atfft_float_to_sample_complex_stride (atfft_complex_f *in,
                                           int in_stride,
                                           atfft_complex *out,
                                           int out_stride,
                                           int size);

/**
 * Convert a complex valued signal from \ref atfft_sample values to single precision floats.
 *
 * @param in the input signal (should contain at least @p size elements)
 * @param out the output signal (should contain at least @p size elements)
 * @param size the length of the signal
 */
void atfft_sample_to_float_complex (atfft_complex *in, atfft_complex_f *out, int size);

/**
 * Convert a complex valued signal from \ref atfft_sample values to single precision floats, taking strides of different
 * length for input and output.
 *
 * @param in the input signal
 *        (should contain at least (@p size - 1) * @p in_stride + 1 elements)
 * @param in_stride the stride to take when reading the input
 * @param out the output signal
 *        (should contain at least (@p size - 1) * @p out_stride + 1 elements)
 * @param out_stride the stride to take when writing the output
 * @param size the length of the signal
 */
void atfft_sample_to_float_complex_stride (atfft_complex *in,
                                           int in_stride,
                                           atfft_complex_f *out,
                                           int out_stride,
                                           int size);

/**
 * Convert a real valued signal from double precision floats to \ref atfft_sample values.
 *
 * @param in the input signal (should contain at least @p size elements)
 * @param out the output signal (should contain at least @p size elements)
 * @param size the length of the signal
 */
void atfft_double_to_sample_real (const double *in, atfft_sample *out, int size);

/**
 * Convert a real valued signal from double precision floats to \ref atfft_sample values, taking strides of different length
 * for input and output.
 *
 * @param in the input signal
 *        (should contain at least (@p size - 1) * @p in_stride + 1 elements)
 * @param in_stride the stride to take when reading the input
 * @param out the output signal
 *        (should contain at least (@p size - 1) * @p out_stride + 1 elements)
 * @param out_stride the stride to take when writing the output
 * @param size the length of the signal
 */
void atfft_double_to_sample_real_stride (const double *in,
                                         int in_stride,
                                         atfft_sample *out,
                                         int out_stride,
                                         int size);

/**
 * Convert a real valued signal from \ref atfft_sample values to double precision floats.
 *
 * @param in the input signal (should contain at least @p size elements)
 * @param out the output signal (should contain at least @p size elements)
 * @param size the length of the signal
 */
void atfft_sample_to_double_real (const atfft_sample *in, double *out, int size);

/**
 * Convert a real valued signal from \ref atfft_sample values to double precision floats, taking strides of different length
 * for input and output.
 *
 * @param in the input signal
 *        (should contain at least (@p size - 1) * @p in_stride + 1 elements)
 * @param in_stride the stride to take when reading the input
 * @param out the output signal
 *        (should contain at least (@p size - 1) * @p out_stride + 1 elements)
 * @param out_stride the stride to take when writing the output
 * @param size the length of the signal
 */
void atfft_sample_to_double_real_stride (const atfft_sample *in,
                                         int in_stride,
                                         double *out,
                                         int out_stride,
                                         int size);

/**
 * Convert a complex valued signal from double precision floats to \ref atfft_sample values.
 *
 * @param in the input signal (should contain at least @p size elements)
 * @param out the output signal (should contain at least @p size elements)
 * @param size the length of the signal
 */
void atfft_double_to_sample_complex (atfft_complex_d *in, atfft_complex *out, int size);

/**
 * Convert a complex valued signal from double precision floats to \ref atfft_sample values, taking strides of different
 * length for input and output.
 *
 * @param in the input signal
 *        (should contain at least (@p size - 1) * @p in_stride + 1 elements)
 * @param in_stride the stride to take when reading the input
 * @param out the output signal
 *        (should contain at least (@p size - 1) * @p out_stride + 1 elements)
 * @param out_stride the stride to take when writing the output
 * @param size the length of the signal
 */
void atfft_double_to_sample_complex_stride (atfft_complex_d *in,
                                            int in_stride,
                                            atfft_complex *out,
                                            int out_stride,
                                            int size);

/**
 * Convert a complex valued signal from \ref atfft_sample values to double precision floats.
 *
 * @param in the input signal (should contain at least @p size elements)
 * @param out the output signal (should contain at least @p size elements)
 * @param size the length of the signal
 */
void atfft_sample_to_double_complex (atfft_complex *in, atfft_complex_d *out, int size);

/**
 * Convert a complex valued signal from \ref atfft_sample values to double precision floats, taking strides of different
 * length for input and output.
 *
 * @param in the input signal
 *        (should contain at least (@p size - 1) * @p in_stride + 1 elements)
 * @param in_stride the stride to take when reading the input
 * @param out the output signal
 *        (should contain at least (@p size - 1) * @p out_stride + 1 elements)
 * @param out_stride the stride to take when writing the output
 * @param size the length of the signal
 */
void atfft_sample_to_double_complex_stride (atfft_complex *in,
                                            int in_stride,
                                            atfft_complex_d *out,
                                            int out_stride,
                                            int size);

/**
 * Convert a real valued signal from long double precision floats to \ref atfft_sample values.
 *
 * @param in the input signal (should contain at least @p size elements)
 * @param out the output signal (should contain at least @p size elements)
 * @param size the length of the signal
 */
void atfft_long_double_to_sample_real (const long double *in, atfft_sample *out, int size);

/**
 * Convert a real valued signal from long double precision floats to \ref atfft_sample values, taking strides of different
 * length for input and output.
 *
 * @param in the input signal
 *        (should contain at least (@p size - 1) * @p in_stride + 1 elements)
 * @param in_stride the stride to take when reading the input
 * @param out the output signal
 *        (should contain at least (@p size - 1) * @p out_stride + 1 elements)
 * @param out_stride the stride to take when writing the output
 * @param size the length of the signal
 */
void atfft_long_double_to_sample_real_stride (const long double *in,
                                              int in_stride,
                                              atfft_sample *out,
                                              int out_stride,
                                              int size);

/**
 * Convert a real valued signal from \ref atfft_sample values to long double precision floats.
 *
 * @param in the input signal (should contain at least @p size elements)
 * @param out the output signal (should contain at least @p size elements)
 * @param size the length of the signal
 */
void atfft_sample_to_long_double_real (const atfft_sample *in, long double *out, int size);

/**
 * Convert a real valued signal from \ref atfft_sample values to long double precision floats, taking strides of different
 * length for input and output.
 *
 * @param in the input signal
 *        (should contain at least (@p size - 1) * @p in_stride + 1 elements)
 * @param in_stride the stride to take when reading the input
 * @param out the output signal
 *        (should contain at least (@p size - 1) * @p out_stride + 1 elements)
 * @param out_stride the stride to take when writing the output
 * @param size the length of the signal
 */
void atfft_sample_to_long_double_real_stride (const atfft_sample *in,
                                              int in_stride,
                                              long double *out,
                                              int out_stride,
                                              int size);

/**
 * Convert a complex valued signal from long double precision floats to \ref atfft_sample values.
 *
 * @param in the input signal (should contain at least @p size elements)
 * @param out the output signal (should contain at least @p size elements)
 * @param size the length of the signal
 */
void atfft_long_double_to_sample_complex (atfft_complex_l *in, atfft_complex *out, int size);

/**
 * Convert a complex valued signal from long double precision floats to \ref atfft_sample values, taking strides of
 * different length for input and output.
 *
 * @param in the input signal
 *        (should contain at least (@p size - 1) * @p in_stride + 1 elements)
 * @param in_stride the stride to take when reading the input
 * @param out the output signal
 *        (should contain at least (@p size - 1) * @p out_stride + 1 elements)
 * @param out_stride the stride to take when writing the output
 * @param size the length of the signal
 */
void atfft_long_double_to_sample_complex_stride (atfft_complex_l *in,
                                                int in_stride,
                                                atfft_complex *out,
                                                int out_stride,
                                                int size);

/**
 * Convert a complex valued signal from \ref atfft_sample values to long double precision floats.
 *
 * @param in the input signal (should contain at least @p size elements)
 * @param out the output signal (should contain at least @p size elements)
 * @param size the length of the signal
 */
void atfft_sample_to_long_double_complex (atfft_complex *in, atfft_complex_l *out, int size);

/**
 * Convert a complex valued signal from \ref atfft_sample values to long double precision floats, taking strides of
 * different length for input and output.
 *
 * @param in the input signal
 *        (should contain at least (@p size - 1) * @p in_stride + 1 elements)
 * @param in_stride the stride to take when reading the input
 * @param out the output signal
 *        (should contain at least (@p size - 1) * @p out_stride + 1 elements)
 * @param out_stride the stride to take when writing the output
 * @param size the length of the signal
 */
void atfft_sample_to_long_double_complex_stride (atfft_complex *in,
                                                int in_stride,
                                                atfft_complex_l *out,
                                                int out_stride,
                                                int size);

#ifdef __cplusplus
}
#endif

#endif /* ATFFT_SHARED_H_INCLUDED */
