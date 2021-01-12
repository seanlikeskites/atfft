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

#ifndef ATFFT_KFR_DEFINITIONS_H_INCLUDED
#define ATFFT_KFR_DEFINITIONS_H_INCLUDED

#ifdef ATFFT_TYPE_LONG_DOUBLE
#   ifdef _MSC_VER
#       pragma message(": warning: KFR only supports single and double precision floating point, " \
                       "higher precision values will be demoted to double for FFT calculations.")
#   else
#       warning KFR only supports single and double precision floating point, \
                higher precision values will be demoted to double for FFT calculations.
#   endif
#endif

#if defined(ATFFT_TYPE_FLOAT)
#   define ATFFT_KFR_DFT_CREATE_PLAN kfr_dft_create_plan_f32
#   define ATFFT_KFR_DFT_GET_TEMP_SIZE kfr_dft_get_temp_size_f32
#   define ATFFT_KFR_DFT_DELETE_PLAN kfr_dft_delete_plan_f32
#   define ATFFT_KFR_DFT_EXECUTE kfr_dft_execute_f32
#   define ATFFT_KFR_DFT_EXECUTE_INVERSE kfr_dft_execute_inverse_f32
#   define ATFFT_KFR_DFT_REAL_CREATE_PLAN kfr_dft_real_create_plan_f32
#   define ATFFT_KFR_DFT_REAL_GET_TEMP_SIZE kfr_dft_real_get_temp_size_f32
#   define ATFFT_KFR_DFT_REAL_DELETE_PLAN kfr_dft_real_delete_plan_f32
#   define ATFFT_KFR_DFT_REAL_EXECUTE kfr_dft_real_execute_f32
#   define ATFFT_KFR_DFT_REAL_EXECUTE_INVERSE kfr_dft_real_execute_inverse_f32
#   define ATFFT_KFR_DCT_CREATE_PLAN kfr_dct_create_plan_f32
#   define ATFFT_KFR_DCT_GET_TEMP_SIZE kfr_dct_get_temp_size_f32
#   define ATFFT_KFR_DCT_DELETE_PLAN kfr_dct_delete_plan_f32
#   define ATFFT_KFR_DCT_EXECUTE kfr_dct_execute_f32
#   define ATFFT_KFR_DCT_EXECUTE_INVERSE kfr_dct_execute_inverse_f32
    typedef kfr_f32 atfft_kfr_sample;

#else
#   define ATFFT_KFR_DFT_CREATE_PLAN kfr_dft_create_plan_f64
#   define ATFFT_KFR_DFT_GET_TEMP_SIZE kfr_dft_get_temp_size_f64
#   define ATFFT_KFR_DFT_DELETE_PLAN kfr_dft_delete_plan_f64
#   define ATFFT_KFR_DFT_EXECUTE kfr_dft_execute_f64
#   define ATFFT_KFR_DFT_EXECUTE_INVERSE kfr_dft_execute_inverse_f64
#   define ATFFT_KFR_DFT_REAL_CREATE_PLAN kfr_dft_real_create_plan_f64
#   define ATFFT_KFR_DFT_REAL_GET_TEMP_SIZE kfr_dft_real_get_temp_size_f64
#   define ATFFT_KFR_DFT_REAL_DELETE_PLAN kfr_dft_real_delete_plan_f64
#   define ATFFT_KFR_DFT_REAL_EXECUTE kfr_dft_real_execute_f64
#   define ATFFT_KFR_DFT_REAL_EXECUTE_INVERSE kfr_dft_real_execute_inverse_f64
#   define ATFFT_KFR_DCT_CREATE_PLAN kfr_dct_create_plan_f64
#   define ATFFT_KFR_DCT_GET_TEMP_SIZE kfr_dct_get_temp_size_f64
#   define ATFFT_KFR_DCT_DELETE_PLAN kfr_dct_delete_plan_f64
#   define ATFFT_KFR_DCT_EXECUTE kfr_dct_execute_f64
#   define ATFFT_KFR_DCT_EXECUTE_INVERSE kfr_dct_execute_inverse_f64
    typedef kfr_f64 atfft_kfr_sample;
#endif

#endif /* ATFFT_KFR_DEFINITIONS_H_INCLUDED */
