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

#ifndef TEST_FUNCTIONS_H_INCLUDED
#define TEST_FUNCTIONS_H_INCLUDED

#include <atfft/atfft.h>

enum test_result
{
    TEST_SUCCESS,
    TEST_FAILURE
};

void generate_real_dc (atfft_sample *sig, int size, atfft_sample offset);
void generate_complex_dc (atfft_complex *sig, int size, atfft_complex offset);

void generate_real_impulse (atfft_sample *sig, int size);
void generate_complex_impulse (atfft_complex *sig, int size);

atfft_sample max_error_real (const atfft_sample *a, const atfft_sample *b, int size);
atfft_sample max_error_complex (atfft_complex *a, atfft_complex *b, int size);

enum test_result test_complex_dft (atfft_complex *in,
                                   atfft_complex *expected,
                                   int size,
                                   enum atfft_direction,
                                   atfft_sample threshold);

#endif /* TEST_FUNCTIONS_H_INCLUDED */
