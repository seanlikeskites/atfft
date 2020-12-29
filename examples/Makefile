CC = gcc -g
CFLAGS = -Wall -O3
LIBS = -lm
INCLUDEDIRS = -I../include
DFT_SOURCES = ../src/atfft/atfft_shared.c
ND_DFT_SOURCES = ../src/atfft/atfft_shared.c
DCT_SOURCES = ../src/atfft/atfft_shared.c
PROGRAMS = complex_transform dct real_transform dct nd_complex nd_real

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

	DFT_SOURCES += ../src/fftw/dft.c
	ND_DFT_SOURCES += ../src/fftw/dft_nd.c
	DCT_SOURCES += ../src/fftw/dct.c
	
else ifeq (FFTS, $(IMPL))
	LIBS += -lffts
	DFT_SOURCES += ../src/ffts/dft.c
	DCT_SOURCES += ../src/atfft/dct.c ../src/ffts/dft.c

else ifeq (GSL, $(IMPL))
	LIBS += -lgsl -lgslcblas
	DFT_SOURCES += ../src/gsl/dft.c
	ND_DFT_SOURCES += ../src/gsl/dft.c ../src/atfft/dft_nd.c ../src/atfft/atfft_internal.c
	DCT_SOURCES += ../src/atfft/dct.c ../src/gsl/dft.c

else ifeq (KISS, $(IMPL))
	ifeq (FLOAT, $(TYPE))
		LIBS += -lkissfftf
	else ifeq (LONG_DOUBLE, $(TYPE))
		LIBS += -lkissfftl
	else
		LIBS += -lkissfft
	endif

	DFT_SOURCES += ../src/kiss/dft.c
	ND_DFT_SOURCES += ../src/kiss/dft_nd.c
	DCT_SOURCES += ../src/atfft/dct.c ../src/kiss/dft.c
	
else ifeq (VDSP, $(IMPL))
	LIBS += -framework Accelerate
	DFT_SOURCES += ../src/vdsp/dft.c
	DCT_SOURCES += ../src/vdsp/dct.c

else ifeq (FFMPEG, $(IMPL))
	LIBS += -lavutil -lavcodec
	DFT_SOURCES += ../src/ffmpeg/dft.c
	ND_DFT_SOURCES += ../src/ffmpeg/dft.c ../src/atfft/dft_nd.c ../src/atfft/atfft_internal.c
	DCT_SOURCES += ../src/ffmpeg/dct.c

else ifeq (PFFFT, $(IMPL))
	DFT_SOURCES += ../src/pffft/dft.c ../src/pffft/pffft.c
	ND_DFT_SOURCES += ../src/pffft/dft.c ../src/pffft/pffft.c ../src/atfft/dft_nd.c ../src/atfft/atfft_internal.c
	DCT_SOURCES += ../src/atfft/dct.c ../src/pffft/dft.c ../src/pffft/pffft.c
	CFLAGS += -msse

else ifeq (OOURA, $(IMPL))
	DFT_SOURCES += ../src/ooura/dft.c ../src/ooura/fft4g.c
	ND_DFT_SOURCES += ../src/ooura/dft.c ../src/ooura/fft4g.c ../src/atfft/dft_nd.c ../src/atfft/atfft_internal.c
	DCT_SOURCES += ../src/ooura/dct.c ../src/ooura/fft4g.c

else
	DFT_SOURCES += ../src/atfft/dft.c ../src/atfft/dft_rader.c ../src/atfft/atfft_internal.c ../src/atfft/dft_cooley_tukey.c ../src/atfft/dft_bluestein.c
	ND_DFT_SOURCES += ../src/atfft/dft.c ../src/atfft/dft_rader.c ../src/atfft/atfft_internal.c ../src/atfft/dft_cooley_tukey.c ../src/atfft/dft_bluestein.c ../src/atfft/dft_nd.c
	DCT_SOURCES += ../src/atfft/dct.c ../src/atfft/dft.c ../src/atfft/dft_rader.c ../src/atfft/atfft_internal.c ../src/atfft/dft_cooley_tukey.c ../src/atfft/dft_bluestein.c

endif

all: $(PROGRAMS)

complex_transform: complex_transform.c $(DFT_SOURCES)
	$(CC) -o $@ $< $(DFT_SOURCES) $(INCLUDEDIRS) $(LIBS) $(CFLAGS)

real_transform: real_transform.c $(DFT_SOURCES)
	$(CC) -o $@ $< $(DFT_SOURCES) $(INCLUDEDIRS) $(LIBS) $(CFLAGS)

nd_complex: nd_complex.c $(ND_DFT_SOURCES)
	$(CC) -o $@ $< $(ND_DFT_SOURCES) $(INCLUDEDIRS) $(LIBS) $(CFLAGS)

nd_real: nd_real.c $(ND_DFT_SOURCES)
	$(CC) -o $@ $< $(ND_DFT_SOURCES) $(INCLUDEDIRS) $(LIBS) $(CFLAGS)

dct: dct.c $(DCT_SOURCES)
	$(CC) -o $@ $< $(DCT_SOURCES) $(INCLUDEDIRS) $(LIBS) $(CFLAGS)

clean:
	rm -f $(PROGRAMS)
