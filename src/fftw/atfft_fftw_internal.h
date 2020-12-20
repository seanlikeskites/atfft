/*
 * Copyright (C) 2020 Sean Enderby <sean.enderby@gmail.com>
 *
 * This program is free software. It comes without any warranty, to
 * the extent permitted by applicable law. You can redistribute it 
 * and/or modify it under the terms of the Do What The Fuck You Want 
 * To Public License, Version 2, as published by Sam Hocevar. See 
 * the COPYING file for more details.
 */

#ifndef ATFFT_FFTW_INTERNAL_H_INCLUDED
#define ATFFT_FFTW_INTERNAL_H_INCLUDED

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

#endif /* ATFFT_DFT_H_INCLUDED */
