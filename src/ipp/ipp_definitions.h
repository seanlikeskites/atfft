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

#ifndef ATFFT_IPP_DEFINITIONS_H_INCLUDED
#define ATFFT_IPP_DEFINITIONS_H_INCLUDED

#ifdef ATFFT_TYPE_LONG_DOUBLE
#   ifdef _MSC_VER
#       pragma message(": warning: IPP only supports single and double precision floating point, " \
                       "higher precision values will be demoted to double for FFT calculations.")
#   else
#       warning IPP only supports single and double precision floating point, \
                higher precision values will be demoted to double for FFT calculations.
#   endif
#endif

#if defined(ATFFT_TYPE_FLOAT)
#   define ATFFT_IPPS_DFT_SPEC_C IppsDFTSpec_C_32fc
#   define ATFFT_IPPS_DFT_GET_SIZE_C ippsDFTGetSize_C_32fc
#   define ATFFT_IPPS_DFT_INIT_C ippsDFTInit_C_32fc
#   define ATFFT_IPPS_DFT_FWD_CTOC ippsDFTFwd_CToC_32fc
#   define ATFFT_IPPS_DFT_INV_CTOC ippsDFTInv_CToC_32fc

#   define ATFFT_IPPS_DFT_SPEC_R IppsDFTSpec_R_32f
#   define ATFFT_IPPS_DFT_GET_SIZE_R ippsDFTGetSize_R_32f
#   define ATFFT_IPPS_DFT_INIT_R ippsDFTInit_R_32f
#   define ATFFT_IPPS_DFT_FWD_RTOCCS ippsDFTFwd_RToCCS_32f
#   define ATFFT_IPPS_DFT_INV_CCSTOR ippsDFTInv_CCSToR_32f
    typedef Ipp32f atfft_ipp_sample;
    typedef Ipp32fc atfft_ipp_complex;
    
#else
#   define ATFFT_IPPS_DFT_SPEC_C IppsDFTSpec_C_64fc
#   define ATFFT_IPPS_DFT_GET_SIZE_C ippsDFTGetSize_C_64fc
#   define ATFFT_IPPS_DFT_INIT_C ippsDFTInit_C_64fc
#   define ATFFT_IPPS_DFT_FWD_CTOC ippsDFTFwd_CToC_64fc
#   define ATFFT_IPPS_DFT_INV_CTOC ippsDFTInv_CToC_64fc

#   define ATFFT_IPPS_DFT_SPEC_R IppsDFTSpec_R_64f
#   define ATFFT_IPPS_DFT_GET_SIZE_R ippsDFTGetSize_R_64f
#   define ATFFT_IPPS_DFT_INIT_R ippsDFTInit_R_64f
#   define ATFFT_IPPS_DFT_FWD_RTOCCS ippsDFTFwd_RToCCS_64f
#   define ATFFT_IPPS_DFT_INV_CCSTOR ippsDFTInv_CCSToR_64f
    typedef Ipp64f atfft_ipp_sample;
    typedef Ipp64fc atfft_ipp_complex;

#endif

#endif /* ATFFT_IPP_DEFINITIONS_H_INCLUDED */
