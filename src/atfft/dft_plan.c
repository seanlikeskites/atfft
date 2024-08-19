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

#include <stdlib.h>
#include "atfft_internal.h"
#include "dft_plan.h"
#include "dft_cooley_tukey.h"
#include "dft_pfa.h"
#include "dft_rader.h"
#include "dft_bluestein.h"

cJSON* atfft_dft_get_plan (void *fft)
{
    enum atfft_dft_algorithm *algorithm = fft;
    cJSON *plan_structure = NULL;

    switch (*algorithm) {
        case ATFFT_BASE:
            plan_structure = atfft_dft_base_get_plan (fft);
            break;

        case ATFFT_COOLEY_TUKEY:
            plan_structure = atfft_dft_cooley_tukey_get_plan (fft);
            break;

        case ATFFT_PFA:
            plan_structure = atfft_dft_pfa_get_plan (fft);
            break;

        case ATFFT_RADER:
            plan_structure = atfft_dft_rader_get_plan (fft);
            break;

        case ATFFT_BLUESTEIN:
            plan_structure = atfft_dft_bluestein_get_plan (fft);
            break;
    }

    return plan_structure;
}
