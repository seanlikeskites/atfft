/*
 * Copyright (c) 2021 Sean Enderby <sean.enderby@gmail.com>
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

#ifndef ATFFT_MKL_DEFINITIONS_H_INCLUDED
#define ATFFT_MKL_DEFINITIONS_H_INCLUDED

#ifdef ATFFT_TYPE_LONG_DOUBLE
#   ifdef _MSC_VER
#       pragma message(": warning: PFFFT only supports single and double precision floating point, " \
                       "higher precision values will be demoted to double for FFT calculations.")
#   else
#       warning MKL only supports single and double precision floating point, \
                higher precision values will be demoted to double for FFT calculations.
#   endif
#endif

/* tell MKL about our complex types */
#include <atfft/types.h>
#define MKL_Complex8 atfft_complex_f
#define MKL_Complex16 atfft_complex_d
#include "mkl.h"

#if defined(ATFFT_TYPE_FLOAT)
#   define ATFFT_MKL_PRECISION DFTI_SINGLE
#else
#   define ATFFT_MKL_PRECISION DFTI_DOUBLE
#endif

#endif /* ATFFT_MKL_DEFINITIONS_H_INCLUDED */
