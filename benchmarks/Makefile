CC = gcc
CPP = g++
CFLAGS = -Wall -O3
LIBS = -lm -lpthread -lbenchmark
INCLUDEDIRS = -I../include
FFT_SOURCES = ../src/atfft/atfft_shared.c

ifeq (FLOAT, $(TYPE))
	CFLAGS += -DATFFT_TYPE_FLOAT
else ifeq (LONG_DOUBLE, $(TYPE))
	CFLAGS += -DATFFT_TYPE_LONG_DOUBLE
endif

ifeq (FFTW3, $(IMPL))
	ifeq (FLOAT, $(TYPE))
		LIBS += -lfftw3f
	else ifeq (LONG_DOUBLE, $(TYPE))
		LIBS += -lfftw3l
	else
		LIBS += -lfftw3
	endif

	FFT_SOURCES += ../src/fftw/atfft_dft_fftw.c
	
else ifeq (FFTS, $(IMPL))
	LIBS += -lffts
	FFT_SOURCES += ../src/ffts/atfft_dft_ffts.c

else ifeq (GSL, $(IMPL))
	LIBS += -lgsl -lgslcblas
	FFT_SOURCES += ../src/gsl/atfft_dft_gsl.c

else ifeq (KISS, $(IMPL))
	ifeq (FLOAT, $(TYPE))
		LIBS += -lkissfftf
	else ifeq (LONG_DOUBLE, $(TYPE))
		LIBS += -lkissfftl
	else
		LIBS += -lkissfft
	endif

	FFT_SOURCES += ../src/kiss/atfft_dft_kiss.c
	
else ifeq (VDSP, $(IMPL))
	LIBS += -framework Accelerate
	FFT_SOURCES += ../src/vdsp/atfft_dft_vdsp.c

else ifeq (FFMPEG, $(IMPL))
	LIBS += -lavutil -lavcodec
	FFT_SOURCES += ../src/ffmpeg/atfft_dft_ffmpeg.c

else ifeq (PFFFT, $(IMPL))
	FFT_SOURCES += ../src/pffft/atfft_dft_pffft.c ../src/pffft/pffft.c
	CFLAGS += -msse

else ifeq (OOURA, $(IMPL))
	FFT_SOURCES += ../src/ooura/atfft_dft_ooura.c ../src/ooura/fft4g.c

else
	FFT_SOURCES += ../src/atfft/atfft_dft.c ../src/atfft/atfft_dft_rader.c ../src/atfft/atfft_internal.c
endif

FFT_TEMP1 = $(notdir $(FFT_SOURCES))
FFT_OBJECTS = $(subst .c,.o, $(FFT_TEMP1))

benchmark: benchmark.cpp $(FFT_SOURCES)
	$(CC) -c $(FFT_SOURCES) $(INCLUDEDIRS) $(CFLAGS)
	$(CPP) -o $@ $< $ $(FFT_OBJECTS) $(INCLUDEDIRS) $(LIBS) $(CFLAGS)

clean:
	rm -rf *.o benchmark
