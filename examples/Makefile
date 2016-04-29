CC = gcc -g
CFLAGS = -Wall
LIBS = -lm
INCLUDEDIRS = -I../include
FFT_SOURCES = ../src/atfft_shared.c

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

	FFT_SOURCES += ../src/atfft_dft_fftw.c ../src/atfft_dct_fftw.c
	
else ifeq (FFTS, $(IMPL))
	LIBS += -lffts
	FFT_SOURCES += ../src/atfft_dft_ffts.c

else ifeq (GSL, $(IMPL))
	LIBS += -lgsl -lgslcblas
	FFT_SOURCES += ../src/atfft_dft_gsl.c

else ifeq (KISS, $(IMPL))
	ifeq (FLOAT, $(TYPE))
		LIBS += -lkissfftf
	else ifeq (LONG_DOUBLE, $(TYPE))
		LIBS += -lkissfftl
	else
		LIBS += -lkissfft
	endif

	FFT_SOURCES += ../src/atfft_dft_kiss.c
	
else ifeq (VDSP, $(IMPL))
	LIBS += -framework Accelerate
	FFT_SOURCES += ../src/atfft_dft_vdsp.c ../src/atfft_dct_vdsp.c

else ifeq (FFMPEG, $(IMPL))
	LIBS += -lavutil -lavcodec
	FFT_SOURCES += ../src/atfft_dft_ffmpeg.c ../src/atfft_dct_ffmpeg.c

else
	FFT_SOURCES += ../src/atfft_dft_ooura.c ../src/atfft_dct_ooura.c ../src/ooura/fft4g.c
endif

all: complex_transform real_transform dct

complex_transform: complex_transform.c
	$(CC) -o $@ $< $(FFT_SOURCES) $(INCLUDEDIRS) $(LIBS) $(CFLAGS)

real_transform: real_transform.c
	$(CC) -o $@ $< $(FFT_SOURCES) $(INCLUDEDIRS) $(LIBS) $(CFLAGS)

dct: dct.c
	$(CC) -o $@ $< $(FFT_SOURCES) $(INCLUDEDIRS) $(LIBS) $(CFLAGS)

clean:
	rm -f complex_transform real_transform dct
