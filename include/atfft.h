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

typedef double atfft_complex [2];

#define ATFFT_REAL(x) ((x) [0])
#define ATFFT_IMAG(x) ((x) [1])

struct atfft* atfft_create (int size, int direction, enum atfft_format format);
void atfft_free (struct atfft *fft);

void atfft_complex_transform (struct atfft *fft, atfft_complex *in, atfft_complex *out);
void atfft_real_forward_transform (struct atfft *fft, double *in, atfft_complex *out);
void atfft_real_backward_transform (struct atfft *fft, atfft_complex *in, double *out);

int isEven (unsigned int x);
int isOdd (unsigned int x);
int isPowerOf2 (unsigned int x);

void atfft_real (atfft_complex *in, double *out, int size);
void atfft_imag (atfft_complex *in, double *out, int size);
void atfft_real_to_complex (double *in, atfft_complex *out, int size);
void atfft_halfcomplex_to_complex (atfft_complex *in, atfft_complex *out, int size);

void atfft_normalise_real (double *data, int size);
void atfft_normalise_complex (atfft_complex *data, int size);

#endif /* ATFFT_H_INCLUDED */
