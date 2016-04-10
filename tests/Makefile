CC = gcc
CFLAGS = -Wall
LIBS = -lm
INCLUDEDIRS = -I../include
FFT_SOURCES = ../src/atfft_shared.c

ifeq (FFTW3, $(IMPL))
	LIBS += -lfftw3
	FFT_SOURCES += ../src/atfft_fftw.c
else ifeq (GSL, $(IMPL))
	LIBS += -lgsl -lgslcblas
	FFT_SOURCES += ../src/atfft_gsl.c
else ifeq (KISS, $(IMPL))
	LIBS += -lkissfft
	FFT_SOURCES += ../src/atfft_kiss.c
else
	FFT_SOURCES += ../src/atfft_ooura.c
endif

all: complex_transform real_transform

complex_transform: complex_transform.c
	$(CC) -o $@ $< $(FFT_SOURCES) $(INCLUDEDIRS) $(LIBS) $(CFLAGS)

real_transform: real_transform.c
	$(CC) -o $@ $< $(FFT_SOURCES) $(INCLUDEDIRS) $(LIBS) $(CFLAGS)

clean:
	rm complex_transform real_transform
