/*
 * Copyright (C) 2020 Sean Enderby <sean.enderby@gmail.com>
 *
 * This program is free software. It comes without any warranty, to
 * the extent permitted by applicable law. You can redistribute it 
 * and/or modify it under the terms of the Do What The Fuck You Want 
 * To Public License, Version 2, as published by Sam Hocevar. See 
 * the COPYING file for more details.
 */

#include <math.h>
#include <stdlib.h>
#include <assert.h>
#include <atfft/atfft_dft.h>
#include "atfft_dft_bluestein.h"
#include "../atfft_internal.h"

struct atfft_dft_bluestein
{
    int size;
    enum atfft_direction direction;
    enum atfft_format format;
    int conv_size;
    struct atfft_dft *conv_forward, *conv_backward;
    atfft_complex *sig, *sig_dft, *conv, *conv_dft;
};

struct atfft_dft_bluestein* atfft_dft_bluestein_create (int size,
                                                        enum atfft_direction direction,
                                                        enum atfft_format format)
{
}

void atfft_dft_bluestein_destroy (struct atfft_dft_bluestein *fft)
{
}

void atfft_dft_bluestein_complex_transform (struct atfft_dft_bluestein *fft,
                                            atfft_complex *in,
                                            atfft_complex *out,
                                            int stride)
{
}
