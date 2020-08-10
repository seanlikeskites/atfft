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

/* Some inline functions for doing calculations with complex numbers */
inline void atfft_copy_complex (const atfft_complex x, atfft_complex y)
{
    ATFFT_REAL (y) = ATFFT_REAL (x);
    ATFFT_IMAG (y) = ATFFT_IMAG (x);
}

inline void atfft_sum_complex (const atfft_complex a,
                               const atfft_complex b,
                               atfft_complex s)
{
    ATFFT_REAL (s) = ATFFT_REAL (a) + ATFFT_REAL (b);
    ATFFT_IMAG (s) = ATFFT_IMAG (a) + ATFFT_IMAG (b);
}

inline void atfft_difference_complex (const atfft_complex a,
                                      const atfft_complex b,
                                      atfft_complex d)
{
    ATFFT_REAL (d) = ATFFT_REAL (a) - ATFFT_REAL (b);
    ATFFT_IMAG (d) = ATFFT_IMAG (a) - ATFFT_IMAG (b);
}

inline void atfft_product_complex (const atfft_complex a,
                                   const atfft_complex b,
                                   atfft_complex p)
{
    ATFFT_REAL (p) = ATFFT_REAL (a) * ATFFT_REAL (b) -
                     ATFFT_IMAG (a) * ATFFT_IMAG (b);
    ATFFT_IMAG (p) = ATFFT_REAL (a) * ATFFT_IMAG (b) +
                     ATFFT_IMAG (a) * ATFFT_REAL (b);
}

#endif /* ATFFT_INTERNAL_H_INCLUDED */
