/*
 * Copyright (C) 2016 Sean Enderby <sean.enderby@gmail.com>
 *
 * This program is free software. It comes without any warranty, to
 * the extent permitted by applicable law. You can redistribute it 
 * and/or modify it under the terms of the Do What The Fuck You Want 
 * To Public License, Version 2, as published by Sam Hocevar. See 
 * the COPYING file for more details.
 */

#ifndef ATFFT_H_INCLUDED
#define ATFFT_H_INCLUDED

#define ATFFT_FORWARD -1
#define ATFFT_BACKWARD 1

/**
 *
 *
 */
enum atfft_format
{
    ATFFT_COMPLEX,
    ATFFT_REAL
};

/** 
 * A Structure to hold internal FFT implementation.
 * 
 * When using atfft you will create one of these structures
 * using atfft_create(), this structure is then passed 
 * to the calculation functions in order to compute DFTs.
 */
struct atfft;

/** A complex float type. */
typedef float atfft_complex_f [2];
/** A complex double type. */
typedef double atfft_complex_d [2];
/** A complex long double type. */
typedef long double atfft_complex_l [2];

/** A macro to return the real part of the complex types. */
#define ATFFT_REAL(x) ((x) [0])
/** A macro to return the imaginary part of the complex types. */
#define ATFFT_IMAG(x) ((x) [1])

#if defined(ATFFT_TYPE_FLOAT)
    typedef float atfft_sample;
    typedef atfft_complex_f atfft_complex;

#elif defined(ATFFT_TYPE_LONG_DOUBLE)
    typedef long double atfft_sample;
    typedef atfft_complex_l atfft_complex;

#else
#   define ATFFT_TYPE_DOUBLE
    typedef double atfft_sample;
    typedef atfft_complex_d atfft_complex;
#endif

/**
 * Check if an integer is even.
 *
 * @param x integer to check
 */
int atfft_is_even (unsigned int x);

/**
 * Check if an integer is odd.
 *
 * @param x integer to check
 */
int atfft_is_odd (unsigned int x);

/**
 * Check if an integer is a power of 2
 *
 * @param x integer to check
 */
int atfft_is_power_of_2 (unsigned int x);

/**
 * Normalise a real DFT output.
 *
 * Applies 1 / size scaling to a real valued signal.
 *
 * @param data the signal to normalise (should contain size elements)
 * @param size the length of the signal
 */
void atfft_normalise_real (atfft_sample *data, int size);

/**
 * Normalise a complex DFT output.
 * 
 * Applies 1 / size scaling to a complex valued signal.
 *
 * @param data the signal to normalise (should contain size elements)
 * @param size the length of the signal
 */
void atfft_normalise_complex (atfft_complex *data, int size);

/**
 * Get the real part of a complex signal.
 *
 * @param in the input signal (should contain size elements)
 * @param out the output signal (should contain size elements)
 * @param size the length of the signals
 */
void atfft_real (atfft_complex *in, atfft_sample *out, int size);

/**
 * Get the imaginary part of a complex signal.
 *
 * @param in the input signal (should contain size elements)
 * @param out the output signal (should contain size elements)
 * @param size the length of the signals
 */
void atfft_imag (atfft_complex *in, atfft_sample *out, int size);

/**
 * Create a complex signal from a real signal.
 *
 * @param in the input signal (should contain size elements)
 * @param out the output signal (should contain size elements)
 * @param size the length of the signals
 */
void atfft_real_to_complex (atfft_sample *in, atfft_complex *out, int size);

/**
 * Create a complex signal from a halfcomplex signal.
 *
 * @param in the input signal (should contain size / 2 + 1 elements)
 * @param out the output signal (should contain size elements)
 * @param size the length of the output signal
 */
void atfft_halfcomplex_to_complex (atfft_complex *in, atfft_complex *out, int size);

#ifndef ATFFT_TYPE_FLOAT
void atfft_float_to_sample_real (float *in, atfft_sample *out, int size);
void atfft_sample_to_float_real (atfft_sample *in, float *out, int size);
void atfft_float_to_sample_complex (atfft_complex_f *in, atfft_complex *out, int size);
void atfft_sample_to_float_complex (atfft_complex *in, atfft_complex_f *out, int size);
#endif

#ifndef ATFFT_TYPE_DOUBLE
void atfft_double_to_sample_real (double *in, atfft_sample *out, int size);
void atfft_sample_to_double_real (atfft_sample *in, double *out, int size);
void atfft_double_to_sample_complex (atfft_complex_d *in, atfft_complex *out, int size);
void atfft_sample_to_double_complex (atfft_complex *in, atfft_complex_d *out, int size);
#endif

#ifndef ATFFT_TYPE_LONG_DOUBLE
void atfft_long_double_to_sample_real (long double *in, atfft_sample *out, int size);
void atfft_sample_to_long_double_real (atfft_sample *in, long double *out, int size);
void atfft_long_double_to_sample_complex (atfft_complex_l *in, atfft_complex *out, int size);
void atfft_sample_to_long_double_complex (atfft_complex *in, atfft_complex_l *out, int size);
#endif

/**
 * Create an fft structure.
 *
 * @param size the signal length the fft should operate on
 * @param direction the direction of the transform
 * @format the type of transform (real or complex)
 */
struct atfft* atfft_create (int size, int direction, enum atfft_format format);

/**
 * Free an fft structure.
 *
 * @param fft the structure to free
 */
void atfft_free (struct atfft *fft);

void atfft_complex_transform (struct atfft *fft, atfft_complex *in, atfft_complex *out);
void atfft_real_forward_transform (struct atfft *fft, atfft_sample *in, atfft_complex *out);
void atfft_real_backward_transform (struct atfft *fft, atfft_complex *in, atfft_sample *out);

#endif /* ATFFT_H_INCLUDED */
