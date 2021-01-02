/*
 * Copyright (C) 2020 Sean Enderby <sean.enderby@gmail.com>
 *
 * This program is free software. It comes without any warranty, to
 * the extent permitted by applicable law. You can redistribute it
 * and/or modify it under the terms of the Do What The Fuck You Want
 * To Public License, Version 2, as published by Sam Hocevar. See
 * the COPYING file for more details.
 */

#include <atfft/atfft_dft_util.h>

int atfft_int_array_product (const int *array, int size)
{
    int prod = array [0];

    for (int i = 1; i < size; ++i)
    {
        prod *= array [i];
    }

    return prod;
}

int atfft_dft_nd_halfcomplex_size (const int *dims, int n_dims)
{
    int prod = 1;

    for (int i = 0; i < n_dims - 1; ++i)
    {
        prod *= dims [i];
    }

    prod *= atfft_dft_halfcomplex_size (dims [n_dims - 1]);

    return prod;
}
