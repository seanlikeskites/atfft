/*
 * Copyright (C) 2016 Sean Enderby <sean.enderby@gmail.com>
 *
 * This program is free software. It comes without any warranty, to
 * the extent permitted by applicable law. You can redistribute it 
 * and/or modify it under the terms of the Do What The Fuck You Want 
 * To Public License, Version 2, as published by Sam Hocevar. See 
 * the COPYING file for more details.
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
    ATFFT_FORWARD, /**< Perform a forward transform. */
    ATFFT_BACKWARD /**< Perform a backward transform. */
};

/** An enum to represent the type of values the transform will operate on. */
enum atfft_format
{
    ATFFT_COMPLEX, /**< Perform a transform on complex valued signals. */
    ATFFT_REAL /**< Perform a transform on real valued signals. */
};

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

/** Some typdefs for changing the type atfft works with. \cond */
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
/** \endcond */

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

/** Return the absolute value of a complex number. */
atfft_sample atfft_abs (const atfft_complex x);

/** Return the complex argument of a complex number. */
atfft_sample atfft_arg (const atfft_complex x);

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
void atfft_real_to_complex (const atfft_sample *in, atfft_complex *out, int size);

/**
 * Create a complex signal from a halfcomplex signal.
 *
 * @param in the input signal (should contain size / 2 + 1 elements)
 * @param out the output signal (should contain size elements)
 * @param size the length of the output signal
 */
void atfft_halfcomplex_to_complex (atfft_complex *in, atfft_complex *out, int size);

#ifndef ATFFT_TYPE_FLOAT
/**
 * Convert a real valued signal from single precision floats 
 * to the type atfft is using.
 *
 * @param in the input signal (should contain size elements)
 * @param out the output signal (should contain size elements)
 * @param size the length of the signals
 */
void atfft_float_to_sample_real (const float *in, atfft_sample *out, int size);

/**
 * Convert a real valued signal from the type atfft is using
 * to single precision floats
 *
 * @param in the input signal (should contain size elements)
 * @param out the output signal (should contain size elements)
 * @param size the length of the signals
 */
void atfft_sample_to_float_real (const atfft_sample *in, float *out, int size);

/**
 * Convert a complex valued signal from single precision floats 
 * to the type atfft is using.
 *
 * @param in the input signal (should contain size elements)
 * @param out the output signal (should contain size elements)
 * @param size the length of the signals
 */
void atfft_float_to_sample_complex (atfft_complex_f *in, atfft_complex *out, int size);

/**
 * Convert a complex valued signal from the type atfft is using
 * to single precision floats
 *
 * @param in the input signal (should contain size elements)
 * @param out the output signal (should contain size elements)
 * @param size the length of the signals
 */
void atfft_sample_to_float_complex (atfft_complex *in, atfft_complex_f *out, int size);
#endif

#ifndef ATFFT_TYPE_DOUBLE
/**
 * Convert a real valued signal from double precision floats 
 * to the type atfft is using.
 *
 * @param in the input signal (should contain size elements)
 * @param out the output signal (should contain size elements)
 * @param size the length of the signals
 */
void atfft_double_to_sample_real (const double *in, atfft_sample *out, int size);

/**
 * Convert a real valued signal from the type atfft is using
 * to double precision floats
 *
 * @param in the input signal (should contain size elements)
 * @param out the output signal (should contain size elements)
 * @param size the length of the signals
 */
void atfft_sample_to_double_real (const atfft_sample *in, double *out, int size);

/**
 * Convert a complex valued signal from double precision floats 
 * to the type atfft is using.
 *
 * @param in the input signal (should contain size elements)
 * @param out the output signal (should contain size elements)
 * @param size the length of the signals
 */
void atfft_double_to_sample_complex (atfft_complex_d *in, atfft_complex *out, int size);

/**
 * Convert a complex valued signal from the type atfft is using
 * to double precision floats
 *
 * @param in the input signal (should contain size elements)
 * @param out the output signal (should contain size elements)
 * @param size the length of the signals
 */
void atfft_sample_to_double_complex (atfft_complex *in, atfft_complex_d *out, int size);
#endif

#ifndef ATFFT_TYPE_LONG_DOUBLE
/**
 * Convert a real valued signal from long double precision floats 
 * to the type atfft is using.
 *
 * @param in the input signal (should contain size elements)
 * @param out the output signal (should contain size elements)
 * @param size the length of the signals
 */
void atfft_long_double_to_sample_real (const long double *in, atfft_sample *out, int size);

/**
 * Convert a real valued signal from the type atfft is using
 * to long double precision floats
 *
 * @param in the input signal (should contain size elements)
 * @param out the output signal (should contain size elements)
 * @param size the length of the signals
 */
void atfft_sample_to_long_double_real (const atfft_sample *in, long double *out, int size);

/**
 * Convert a complex valued signal from long double precision floats 
 * to the type atfft is using.
 *
 * @param in the input signal (should contain size elements)
 * @param out the output signal (should contain size elements)
 * @param size the length of the signals
 */
void atfft_long_double_to_sample_complex (atfft_complex_l *in, atfft_complex *out, int size);

/**
 * Convert a complex valued signal from the type atfft is using
 * to long double precision floats
 *
 * @param in the input signal (should contain size elements)
 * @param out the output signal (should contain size elements)
 * @param size the length of the signals
 */
void atfft_sample_to_long_double_complex (atfft_complex *in, atfft_complex_l *out, int size);
#endif

#ifdef __cplusplus
}
#endif

#endif /* ATFFT_SHARED_H_INCLUDED */
