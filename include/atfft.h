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

typedef double atfft_complex_double [2];
typedef float atfft_complex_float [2];

#define ATFFT_REAL(x) ((x) [0])
#define ATFFT_IMAG(x) ((x) [1])

struct atfft* atfft_create (int size, int direction, enum atfft_format format);
void atfft_free (struct atfft *fft);

void atfft_complex_transform (struct atfft *fft, atfft_complex_double *in, atfft_complex_double *out);
void atfft_real_forward_transform (struct atfft *fft, double *in, atfft_complex_double *out);
void atfft_real_backward_transform (struct atfft *fft, atfft_complex_double *in, double *out);

int atfft_is_even (unsigned int x);
int atfft_is_odd (unsigned int x);
int atfft_is_power_of_2 (unsigned int x);

void atfft_real (atfft_complex_double *in, double *out, int size);
void atfft_imag (atfft_complex_double *in, double *out, int size);
void atfft_real_to_complex (double *in, atfft_complex_double *out, int size);
void atfft_halfcomplex_to_complex (atfft_complex_double *in, atfft_complex_double *out, int size);

void atfft_float_to_double_real (float *in, double *out, int size);
void atfft_double_to_float_real (double *in, float *out, int size);
void atfft_float_to_double_complex (atfft_complex_float *in, atfft_complex_double *out, int size);
void atfft_double_to_float_complex (atfft_complex_double *in, atfft_complex_float *out, int size);

void atfft_normalise_real (double *data, int size);
void atfft_normalise_complex (atfft_complex_double *data, int size);

#endif /* ATFFT_H_INCLUDED */
