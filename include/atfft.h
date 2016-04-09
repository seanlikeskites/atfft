/*
 * Copyright (C) 2016 Sean Enderby <sean.enderby@gmail.com>
 *
 * This work is free. You can redistribute it and/or modify it under the
 * terms of the Do What The Fuck You Want To Public License, Version 2,
 * as published by Sam Hocevar. See the COPYING file for more details.
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

struct atfft* atfft_create (int size, int direction, enum atfft_format format);
void atfft_free (struct atfft *fft);

void atfft_complex_transform (struct atfft *fft, double *in, double *out);
void atfft_real_forward_transform (struct atfft *fft, double *in, double *out);
void atfft_real_backward_transform (struct atfft *fft, double *in, double *out);

int isPowerOf2 (unsigned int x);

#endif /* ATFFT_H_INCLUDED */
