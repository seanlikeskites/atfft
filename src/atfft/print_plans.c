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

#include "atfft_internal.h"
#include "print_plans.h"
#include "dft_cooley_tukey.h"

void atfft_dft_print_plan_internal (void *fft, FILE *stream, int indent)
{
    enum atfft_dft_algorithm *algorithm = fft;

    switch (*algorithm) {
        case ATFFT_BASE:
            atfft_dft_base_print_plan (fft, stream, indent);
            break;

        case ATFFT_COOLEY_TUKEY:
            atfft_dft_cooley_tukey_print_plan (fft, stream, indent);
            break;

        case ATFFT_PFA:
            break;

        case ATFFT_RADER:
            break;

        case ATFFT_BLUESTEIN:
            break;
    }
}
