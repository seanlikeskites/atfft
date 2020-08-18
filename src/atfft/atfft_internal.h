/*
 * Copyright (C) 2020 Sean Enderby <sean.enderby@gmail.com>
 *
 * This program is free software. It comes without any warranty, to
 * the extent permitted by applicable law. You can redistribute it 
 * and/or modify it under the terms of the Do What The Fuck You Want 
 * To Public License, Version 2, as published by Sam Hocevar. See 
 * the COPYING file for more details.
 */

#ifndef ATFFT_INTERNAL_H_INCLUDED
#define ATFFT_INTERNAL_H_INCLUDED

#include <atfft/atfft_shared.h>

/**
 * Copy the complex value x into the complex value y.
 */
#define ATFFT_COPY_COMPLEX(x, y) \
    do \
    { \
        ATFFT_REAL (y) = ATFFT_REAL (x); \
        ATFFT_IMAG (y) = ATFFT_IMAG (x); \
    } \
    while (0)

/**
 * Sum two complex numbers: s = a + b
 */
#define ATFFT_SUM_COMPLEX(a, b, s) \
    do \
    { \
        ATFFT_REAL (s) = ATFFT_REAL (a) + ATFFT_REAL (b); \
        ATFFT_IMAG (s) = ATFFT_IMAG (a) + ATFFT_IMAG (b); \
    } \
    while (0)

/**
 * Take the difference of two complex numbers: d = a - b
 */
#define ATFFT_DIFFERENCE_COMPLEX(a, b, d) \
    do \
    { \
        ATFFT_REAL (d) = ATFFT_REAL (a) - ATFFT_REAL (b); \
        ATFFT_IMAG (d) = ATFFT_IMAG (a) - ATFFT_IMAG (b); \
    } \
    while (0)

/**
 * Find the product of two complex numbers: p = ab
 */
#define ATFFT_PRODUCT_COMPLEX(a, b, p) \
    do \
    { \
        ATFFT_REAL (p) = ATFFT_REAL (a) * ATFFT_REAL (b) - \
                         ATFFT_IMAG (a) * ATFFT_IMAG (b); \
        ATFFT_IMAG (p) = ATFFT_REAL (a) * ATFFT_IMAG (b) + \
                         ATFFT_IMAG (a) * ATFFT_REAL (b); \
    } \
    while (0)

/**
 * Multiply a complex variable by another: a *= b
 */
#define ATFFT_MULTIPLY_BY_COMPLEX(a, b) \
    do \
    { \
        atfft_sample temp; \
        temp = ATFFT_REAL (a); \
        \
        ATFFT_REAL (a) = temp * ATFFT_REAL (b) - \
                         ATFFT_IMAG (a) * ATFFT_IMAG (b); \
        ATFFT_IMAG (a) = temp * ATFFT_IMAG (b) + \
                         ATFFT_IMAG (a) * ATFFT_REAL (b); \
    } \
    while (0)

int atfft_next_power_of_2 (int x);

#endif /* ATFFT_INTERNAL_H_INCLUDED */
