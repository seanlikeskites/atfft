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

enum atfft_format
{
    ATFFT_COMPLEX,
    ATFFT_REAL
};

struct atfft;

typedef float atfft_complex_float [2];
typedef double atfft_complex_double [2];
typedef long double atfft_complex_long_double [2];

#define ATFFT_REAL(x) ((x) [0])
#define ATFFT_IMAG(x) ((x) [1])

#if defined(ATFFT_TYPE_FLOAT)
typedef float atfft_sample;
typedef atfft_complex_float atfft_complex;
#elif defined(ATFFT_TYPE_LONG_DOUBLE)
typedef long double atfft_sample;
typedef atfft_complex_long_double atfft_complex;
#else
#define ATFFT_TYPE_DOUBLE
typedef double atfft_sample;
typedef atfft_complex_double atfft_complex;
#endif

struct atfft* atfft_create (int size, int direction, enum atfft_format format);
void atfft_free (struct atfft *fft);

void atfft_complex_transform (struct atfft *fft, atfft_complex *in, atfft_complex *out);
void atfft_real_forward_transform (struct atfft *fft, double *in, atfft_complex *out);
void atfft_real_backward_transform (struct atfft *fft, atfft_complex *in, double *out);

int atfft_is_even (unsigned int x);
int atfft_is_odd (unsigned int x);
int atfft_is_power_of_2 (unsigned int x);

void atfft_real (atfft_complex *in, double *out, int size);
void atfft_imag (atfft_complex *in, double *out, int size);
void atfft_real_to_complex (double *in, atfft_complex *out, int size);
void atfft_halfcomplex_to_complex (atfft_complex *in, atfft_complex *out, int size);

#ifdef ATFFT_TYPE_FLOAT
void atfft_float_to_sample_real (float *in, atfft_sample *out, int size);
void atfft_sample_to_float_real (atfft_sample *in, float *out, int size);
void atfft_float_to_sample_complex (atfft_complex_float *in, atfft_complex *out, int size);
void atfft_sample_to_float_complex (atfft_complex *in, atfft_complex_float *out, int size);
#endif

#ifdef ATFFT_TYPE_DOUBLE
void atfft_double_to_sample_real (double *in, atfft_sample *out, int size);
void atfft_sample_to_double_real (atfft_sample *in, double *out, int size);
void atfft_double_to_sample_complex (atfft_complex_double *in, atfft_complex *out, int size);
void atfft_sample_to_double_complex (atfft_complex *in, atfft_complex_double *out, int size);
#endif

#ifdef ATFFT_TYPE_LONG_DOUBLE
void atfft_long_double_to_sample_real (long double *in, atfft_sample *out, int size);
void atfft_sample_to_long_double_real (atfft_sample *in, long double *out, int size);
void atfft_long_double_to_sample_complex (atfft_complex_long double *in, atfft_complex *out, int size);
void atfft_sample_to_long_double_complex (atfft_complex *in, atfft_complex_long double *out, int size);
#endif

void atfft_normalise_real (atfft_sample *data, int size);
void atfft_normalise_complex (atfft_complex *data, int size);

#endif /* ATFFT_H_INCLUDED */
