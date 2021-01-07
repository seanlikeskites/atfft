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

#ifndef ATFFT_FFTW_DEFINITIONS_H_INCLUDED
#define ATFFT_FFTW_DEFINITIONS_H_INCLUDED

#if defined(ATFFT_TYPE_FLOAT)
#   define ATFFT_FFTW_MALLOC fftwf_malloc
#   define ATFFT_FFTW_FREE fftwf_free
#   define ATFFT_FFTW_DESTROY_PLAN fftwf_destroy_plan
#   define ATFFT_FFTW_PLAN_DFT_1D fftwf_plan_dft_1d
#   define ATFFT_FFTW_PLAN_DFT_R2C_1D fftwf_plan_dft_r2c_1d
#   define ATFFT_FFTW_PLAN_DFT_C2R_1D fftwf_plan_dft_c2r_1d
#   define ATFFT_FFTW_PLAN_DFT fftwf_plan_dft
#   define ATFFT_FFTW_PLAN_DFT_R2C fftwf_plan_dft_r2c
#   define ATFFT_FFTW_PLAN_DFT_C2R fftwf_plan_dft_c2r
#   define ATFFT_FFTW_PLAN_R2R_1D fftwf_plan_r2r_1d
#   define ATFFT_FFTW_EXECUTE fftwf_execute
    typedef fftwf_plan atfft_fftw_plan;
    
#elif defined (ATFFT_TYPE_DOUBLE)
#   define ATFFT_FFTW_MALLOC fftw_malloc
#   define ATFFT_FFTW_FREE fftw_free
#   define ATFFT_FFTW_DESTROY_PLAN fftw_destroy_plan
#   define ATFFT_FFTW_PLAN_DFT_1D fftw_plan_dft_1d
#   define ATFFT_FFTW_PLAN_DFT_R2C_1D fftw_plan_dft_r2c_1d
#   define ATFFT_FFTW_PLAN_DFT_C2R_1D fftw_plan_dft_c2r_1d
#   define ATFFT_FFTW_PLAN_DFT fftw_plan_dft
#   define ATFFT_FFTW_PLAN_DFT_R2C fftw_plan_dft_r2c
#   define ATFFT_FFTW_PLAN_DFT_C2R fftw_plan_dft_c2r
#   define ATFFT_FFTW_PLAN_R2R_1D fftw_plan_r2r_1d
#   define ATFFT_FFTW_EXECUTE fftw_execute
    typedef fftw_plan atfft_fftw_plan;

#elif defined(ATFFT_TYPE_LONG_DOUBLE)
#   define ATFFT_FFTW_MALLOC fftwl_malloc
#   define ATFFT_FFTW_FREE fftwl_free
#   define ATFFT_FFTW_DESTROY_PLAN fftwl_destroy_plan
#   define ATFFT_FFTW_PLAN_DFT_1D fftwl_plan_dft_1d
#   define ATFFT_FFTW_PLAN_DFT_R2C_1D fftwl_plan_dft_r2c_1d
#   define ATFFT_FFTW_PLAN_DFT_C2R_1D fftwl_plan_dft_c2r_1d
#   define ATFFT_FFTW_PLAN_DFT fftwl_plan_dft
#   define ATFFT_FFTW_PLAN_DFT_R2C fftwl_plan_dft_r2c
#   define ATFFT_FFTW_PLAN_DFT_C2R fftwl_plan_dft_c2r
#   define ATFFT_FFTW_PLAN_R2R_1D fftwl_plan_r2r_1d
#   define ATFFT_FFTW_EXECUTE fftwl_execute
    typedef fftwl_plan atfft_fftw_plan;
#endif

#endif /* ATFFT_FFTW_DEFINITIONS_H_INCLUDED */
