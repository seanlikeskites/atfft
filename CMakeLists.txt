# the usual gumpf
cmake_minimum_required (VERSION 3.5)
project (atfft VERSION 1.0.0)

# configure sources shared between all implementations
add_library (atfft SHARED
	     src/atfft/dft_util.c
	     src/atfft/dft_nd_util.c
	     src/atfft/dct_util.c
	     src/atfft/types.c)

set_target_properties (atfft PROPERTIES VERSION ${PROJECT_VERSION}
	                                SOVERSION 1)

target_include_directories (atfft PUBLIC include)

# Work out which implementation we are building with
# and add necessary sources.
set (IMPL "atfft" CACHE STRING "The FFT implementation to use.")

if (IMPL STREQUAL "atfft")
	target_include_directories (atfft PRIVATE src/atfft)

	target_sources (atfft PRIVATE
		        src/atfft/atfft_internal.c
		        src/atfft/dct.c
		        src/atfft/dft.c
		        src/atfft/dft_bluestein.c
		        src/atfft/dft_cooley_tukey.c
		        src/atfft/dft_nd.c
		        src/atfft/dft_pfa.c
		        src/atfft/dft_rader.c)

elseif (IMPL STREQUAL "ffmpeg")
	find_library (AVUTIL avutil)
	find_library (AVCODEC avcodec)

	if (NOT AVUTIL)
		message (FATAL_ERROR "libavutil not found on system.")
	endif()

	if (NOT AVCODEC)
		message (FATAL_ERROR "libavcodec not found on system.")
	endif()

	target_sources (atfft PRIVATE
		        src/atfft/atfft_internal.c
		        src/atfft/dft_nd.c
		        src/ffmpeg/dct.c
		        src/ffmpeg/dft.c)

	target_link_libraries (atfft ${AVUTIL} ${AVCODEC})

elseif (IMPL STREQUAL "ffts")
	find_library (FFTS ffts)

	if (NOT FFTS)
		message (FATAL_ERROR "libffts not found on system.")
	endif()

	target_sources (atfft PRIVATE
		        src/atfft/dct.c
		        src/ffts/dft.c
		        src/ffts/dft_nd.c)

	target_link_libraries (atfft ${FFTS})

elseif (IMPL STREQUAL "fftw3")
	find_library (FFTW3 fftw3)

	if (NOT FFTW3)
		message (FATAL_ERROR "libfftw3 not found on system.")
	endif()

	target_sources (atfft PRIVATE
		        src/fftw/dct.c
		        src/fftw/dft.c
		        src/fftw/dft_nd.c)

	target_link_libraries (atfft ${FFTW3})

elseif (IMPL STREQUAL "gsl")
	find_library (GSL gsl)
	find_library (GSLCBLAS gslcblas)

	if (NOT GSL)
		message (FATAL_ERROR "libgsl not found on system.")
	endif()

	if (NOT GSLCBLAS)
		message (FATAL_ERROR "libgslcblas not found on system.")
	endif()

	target_sources (atfft PRIVATE
		        src/atfft/dct.c
		        src/atfft/atfft_internal.c
		        src/atfft/dft_nd.c
		        src/gsl/dft.c)

	target_link_libraries (atfft ${GSL} ${GSLCBLAS})

elseif (IMPL STREQUAL "ipp")
	string (REPLACE ":" ";" LIBRARY_DIRS $ENV{LD_LIBRARY_PATH})
	find_library (IPPS ipps PATHS ${LIBRARY_DIRS})
	find_library (IPPVM ippvm PATHS ${LIBRARY_DIRS})
	find_library (IPPCORE ippcore PATHS ${LIBRARY_DIRS})

	if (NOT IPPS)
		message (FATAL_ERROR "libipps not found on system.")
	endif()

	if (NOT IPPVM)
		message (FATAL_ERROR "libippvm not found on system.")
	endif()

	if (NOT IPPCORE)
		message (FATAL_ERROR "libippcore not found on system.")
	endif()

	target_include_directories (atfft PRIVATE src/ipp)

	target_sources (atfft PRIVATE
		        src/atfft/atfft_internal.c
		        src/atfft/dft_nd.c
		        src/ipp/dft.c
		        src/ipp/dct.c
		        src/ipp/ipp_definitions.c)

	target_link_libraries (atfft ${IPPS} ${IPPVM} ${IPPCORE})

elseif (IMPL STREQUAL "kfr")
	find_library (KFR kfr_capi)

	if (NOT KFR)
		message (FATAL_ERROR "libkfr_capi not found on system.")
	endif()

	target_include_directories (atfft PRIVATE src/kfr)

	target_sources (atfft PRIVATE
		        src/atfft/atfft_internal.c
		        src/atfft/dft_nd.c
		        src/kfr/dft.c
		        src/kfr/dct.c)

	target_link_libraries (atfft ${KFR})

elseif (IMPL STREQUAL "kiss")
	find_library (KISSFFT kissfft)

	if (NOT KISSFFT)
		message (FATAL_ERROR "libkissfft not found on system.")
	endif()

	target_sources (atfft PRIVATE
		        src/atfft/dct.c
		        src/kiss/dft.c
		        src/kiss/dft_nd.c)

		target_link_libraries (atfft ${KISSFFT})

elseif (IMPL STREQUAL "mkl")
	string (REPLACE ":" ";" LIBRARY_DIRS $ENV{LD_LIBRARY_PATH})
	find_library (MKL_RT mkl_rt PATHS ${LIBRARY_DIRS})

	if (NOT MKL_RT)
		message (FATAL_ERROR "libmkl_rt not found on system.")
	endif()

	target_include_directories (atfft PRIVATE src/mkl)

	target_sources (atfft PRIVATE
		        src/atfft/dct.c
		        src/mkl/dft.c
		        src/mkl/dft_nd.c)

	target_link_libraries (atfft mkl_rt)

elseif (IMPL STREQUAL "ooura")
	target_include_directories (atfft PRIVATE src/ooura)

	target_sources (atfft PRIVATE
		        src/atfft/atfft_internal.c
		        src/atfft/dft_nd.c
		        src/ooura/dft.c
		        src/ooura/dct.c
			src/ooura/fft4g.c)

elseif (IMPL STREQUAL "pffft")
	target_include_directories (atfft PRIVATE src/pffft)

	target_sources (atfft PRIVATE
		        src/atfft/atfft_internal.c
		        src/atfft/dft_nd.c
		        src/atfft/dct.c
		        src/pffft/dft.c
			src/pffft/pffft.c)

elseif (IMPL STREQUAL "vdsp")
	find_library (ACCELERATE Accelerate)

	if (NOT ACCELERATE)
		message (FATAL_ERROR "Accelerate framework not found on system.")
	endif()

	target_sources (atfft PRIVATE
		        src/atfft/atfft_internal.c
		        src/atfft/dft_nd.c
		        src/vdsp/dct.c
		        src/vdsp/dft.c)

	target_link_libraries (atfft ${ACCELERATE})

else()
	message (FATAL_ERROR "Unknown implementation: ${IMPL}")
endif()
