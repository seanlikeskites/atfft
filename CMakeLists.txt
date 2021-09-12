#############################################
# the usual gumpf
#############################################
cmake_minimum_required (VERSION 3.5)
project (atfft VERSION 1.0.0)

add_compile_options (-Wall -Wextra -pedantic)

#############################################
# Define Sources
#############################################
# public headers to be installed
set (HEADERS include/atfft/atfft.h
             include/atfft/dct.h
             include/atfft/dct_util.h
             include/atfft/dft.h
             include/atfft/dft_util.h
             include/atfft/dft_nd.h
             include/atfft/dft_nd_util.h
             include/atfft/types.h
             include/atfft/windows.h)

# sources shared between all implementations
set (COMMON_SRC src/atfft/dft_util.c
	        src/atfft/dft_nd_util.c
	        src/atfft/dct_util.c
	        src/atfft/types.c
	        src/atfft/windows/windows.c)

# Work out which implementation we are building with
# and add necessary sources.
set (IMPL "atfft" CACHE STRING "The FFT implementation to use.")

if (IMPL STREQUAL "atfft")
	set (IMPL_SRC src/atfft/atfft_internal.c
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

	set (IMPL_SRC src/atfft/atfft_internal.c
		      src/atfft/dft_nd.c
		      src/ffmpeg/dct.c
		      src/ffmpeg/dft.c)

	set (DOUBLE_LINK ${AVUTIL} ${AVCODEC})

elseif (IMPL STREQUAL "ffts")
	find_library (FFTS ffts)

	if (NOT FFTS)
		message (FATAL_ERROR "libffts not found on system.")
	endif()

	set (IMPL_SRC src/atfft/dct.c
		      src/ffts/dft.c
		      src/ffts/dft_nd.c)

	set (DOUBLE_LINK ${FFTS})

elseif (IMPL STREQUAL "fftw3")
	find_library (FFTW3 fftw3)

	if (NOT FFTW3)
		message (FATAL_ERROR "libfftw3 not found on system.")
	endif()

	set (IMPL_SRC src/fftw/dct.c
		      src/fftw/dft.c
		      src/fftw/dft_nd.c)

	set (DOUBLE_LINK ${FFTW3})

elseif (IMPL STREQUAL "gsl")
	find_library (GSL gsl)
	find_library (GSLCBLAS gslcblas)

	if (NOT GSL)
		message (FATAL_ERROR "libgsl not found on system.")
	endif()

	if (NOT GSLCBLAS)
		message (FATAL_ERROR "libgslcblas not found on system.")
	endif()

	set (IMPL_SRC src/atfft/dct.c
		      src/atfft/atfft_internal.c
		      src/atfft/dft_nd.c
		      src/gsl/dft.c)

	set (DOUBLE_LINK ${GSL} ${GSLCBLAS})

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

	set (IMPL_SRC src/atfft/atfft_internal.c
		      src/atfft/dft_nd.c
		      src/ipp/dft.c
		      src/ipp/dct.c
		      src/ipp/ipp_definitions.c)

	set (DOUBLE_LINK ${IPPS} ${IPPVM} ${IPPCORE})

elseif (IMPL STREQUAL "kfr")
	find_library (KFR kfr_capi)

	if (NOT KFR)
		message (FATAL_ERROR "libkfr_capi not found on system.")
	endif()

	set (IMPL_SRC src/atfft/atfft_internal.c
		      src/atfft/dft_nd.c
		      src/kfr/dft.c
		      src/kfr/dct.c)

	set (DOUBLE_LINK ${KFR})

elseif (IMPL STREQUAL "kiss")
	find_library (KISSFFT kissfft)

	if (NOT KISSFFT)
		message (FATAL_ERROR "libkissfft not found on system.")
	endif()

	set (IMPL_SRC src/atfft/dct.c
		      src/kiss/dft.c
		      src/kiss/dft_nd.c)

	set (DOUBLE_LINK ${KISSFFT})

elseif (IMPL STREQUAL "mkl")
	string (REPLACE ":" ";" LIBRARY_DIRS $ENV{LD_LIBRARY_PATH})
	find_library (MKL_RT mkl_rt PATHS ${LIBRARY_DIRS})

	if (NOT MKL_RT)
		message (FATAL_ERROR "libmkl_rt not found on system.")
	endif()

	set (IMPL_SRC src/atfft/dct.c
		      src/mkl/dft.c
		      src/mkl/dft_nd.c)

	set (DOUBLE_LINK mkl_rt)

elseif (IMPL STREQUAL "ooura")
	set (IMPL_SRC src/atfft/atfft_internal.c
		      src/atfft/dft_nd.c
		      src/ooura/dft.c
		      src/ooura/dct.c
		      src/ooura/fft4g.c)

elseif (IMPL STREQUAL "pffft")
	set (IMPL_SRC src/atfft/atfft_internal.c
		      src/atfft/dft_nd.c
		      src/atfft/dct.c
		      src/pffft/dft.c
		      src/pffft/pffft.c)

elseif (IMPL STREQUAL "vdsp")
	find_library (ACCELERATE Accelerate)

	if (NOT ACCELERATE)
		message (FATAL_ERROR "Accelerate framework not found on system.")
	endif()

	set (IMPL_SRC src/atfft/atfft_internal.c
		      src/atfft/dft_nd.c
		      src/vdsp/dct.c
		      src/vdsp/dft.c)

	set (DOUBLE_LINK ${ACCELERATE})

else()
	message (FATAL_ERROR "Unknown implementation: ${IMPL}")
endif()

#############################################
# Build for double type samples
#############################################
add_library (atfft SHARED "${COMMON_SRC}" "${IMPL_SRC}")
set_target_properties (atfft PROPERTIES VERSION ${PROJECT_VERSION}
	                                SOVERSION 1
					PUBLIC_HEADER "${HEADERS}")
target_include_directories (atfft PUBLIC include)
target_compile_definitions (atfft PUBLIC ATFFT_TYPE_DOUBLE)
target_link_libraries (atfft "${DOUBLE_LINK}" m)

##############################################
## Build for float type samples
##############################################
#add_library (atfftf SHARED "${COMMON_SRC}" "${IMPL_SRC}")
#set_target_properties (atfftf PROPERTIES VERSION ${PROJECT_VERSION}
#	                                 SOVERSION 1
#					 PUBLIC_HEADER "${HEADERS}")
#target_include_directories (atfftf PUBLIC include)
#target_compile_definitions (atfftf PUBLIC ATFFT_TYPE_FLOAT)
#target_link_libraries (atfftf "${FLOAT_LINK}")
#
##############################################
## Build for long double type samples
##############################################
#add_library (atfftl SHARED "${COMMON_SRC}" "${IMPL_SRC}")
#set_target_properties (atfftl PROPERTIES VERSION ${PROJECT_VERSION}
#	                                 SOVERSION 1
#					 PUBLIC_HEADER "${HEADERS}")
#target_include_directories (atfftl PUBLIC include)
#target_compile_definitions (atfftl PUBLIC ATFFT_TYPE_LONG_DOUBLE)
#target_link_libraries (atfftl "${LONG_DOUBLE_LINK}")

#############################################
# Installation Config
#############################################
include (GNUInstallDirs)
install (TARGETS atfft
	 LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
	 ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
	 PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/atfft)

#############################################
# Tests
#############################################
enable_testing()
#include (CTest)
add_subdirectory (tests tests EXCLUDE_FROM_ALL)
