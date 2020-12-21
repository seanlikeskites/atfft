CC = gcc -g
CFLAGS = -Wall -O3
LIBS = -lm
INCLUDEDIRS = -I../include
DFT_SOURCES = ../src/atfft/atfft_shared.c
ND_DFT_SOURCES = ../src/atfft/atfft_shared.c
DCT_SOURCES = ../src/atfft/atfft_shared.c
#PROGRAMS = complex_transform real_transform dct nd_complex nd_real
PROGRAMS = complex_transform dct nd_complex nd_real

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

	DFT_SOURCES += ../src/fftw/atfft_dft_fftw.c
	ND_DFT_SOURCES += ../src/fftw/atfft_dft_fftw.c ../src/fftw/atfft_dft_nd_fftw.c
	DCT_SOURCES += ../src/fftw/atfft_dct_fftw.c
	
else ifeq (FFTS, $(IMPL))
	LIBS += -lffts
	DFT_SOURCES += ../src/ffts/atfft_dft_ffts.c
	DCT_SOURCES += ../src/atfft/atfft_dct.c ../src/ffts/atfft_dft_ffts.c

else ifeq (GSL, $(IMPL))
	LIBS += -lgsl -lgslcblas
	DFT_SOURCES += ../src/gsl/atfft_dft_gsl.c
	ND_DFT_SOURCES += ../src/gsl/atfft_dft_gsl.c ../src/atfft/atfft_dft_nd.c ../src/atfft/atfft_internal.c
	DCT_SOURCES += ../src/atfft/atfft_dct.c ../src/gsl/atfft_dft_gsl.c

else ifeq (KISS, $(IMPL))
	ifeq (FLOAT, $(TYPE))
		LIBS += -lkissfftf
	else ifeq (LONG_DOUBLE, $(TYPE))
		LIBS += -lkissfftl
	else
		LIBS += -lkissfft
	endif

	DFT_SOURCES += ../src/kiss/atfft_dft_kiss.c
	DCT_SOURCES += ../src/atfft/atfft_dct.c ../src/kiss/atfft_dft_kiss.c
	
else ifeq (VDSP, $(IMPL))
	LIBS += -framework Accelerate
	DFT_SOURCES += ../src/vdsp/atfft_dft_vdsp.c
	DCT_SOURCES += ../src/vdsp/atfft_dct_vdsp.c

else ifeq (FFMPEG, $(IMPL))
	LIBS += -lavutil -lavcodec
	DFT_SOURCES += ../src/ffmpeg/atfft_dft_ffmpeg.c
	ND_DFT_SOURCES += ../src/ffmpeg/atfft_dft_ffmpeg.c ../src/atfft/atfft_dft_nd.c ../src/atfft/atfft_internal.c
	DCT_SOURCES += ../src/ffmpeg/atfft_dct_ffmpeg.c

else ifeq (PFFFT, $(IMPL))
	DFT_SOURCES += ../src/pffft/atfft_dft_pffft.c ../src/pffft/pffft.c
	ND_DFT_SOURCES += ../src/pffft/atfft_dft_pffft.c ../src/pffft/pffft.c ../src/atfft/atfft_dft_nd.c ../src/atfft/atfft_internal.c
	DCT_SOURCES += ../src/atfft/atfft_dct.c ../src/pffft/atfft_dft_pffft.c ../src/pffft/pffft.c
	CFLAGS += -msse

else ifeq (OOURA, $(IMPL))
	DFT_SOURCES += ../src/ooura/atfft_dft_ooura.c ../src/ooura/fft4g.c
	ND_DFT_SOURCES += ../src/ooura/atfft_dft_ooura.c ../src/ooura/fft4g.c ../src/atfft/atfft_dft_nd.c ../src/atfft/atfft_internal.c
	DCT_SOURCES += ../src/ooura/atfft_dct_ooura.c ../src/ooura/fft4g.c

else
	DFT_SOURCES += ../src/atfft/atfft_dft.c ../src/atfft/atfft_dft_rader.c ../src/atfft/atfft_internal.c ../src/atfft/atfft_dft_cooley_tukey.c ../src/atfft/atfft_dft_bluestein.c ../src/atfft/atfft_dft_nd.c
	ND_DFT_SOURCES += ../src/atfft/atfft_dft.c ../src/atfft/atfft_dft_rader.c ../src/atfft/atfft_internal.c ../src/atfft/atfft_dft_cooley_tukey.c ../src/atfft/atfft_dft_bluestein.c ../src/atfft/atfft_dft_nd.c
	DCT_SOURCES += ../src/atfft/atfft_dct.c ../src/atfft/atfft_dft.c ../src/atfft/atfft_dft_rader.c ../src/atfft/atfft_internal.c ../src/atfft/atfft_dft_cooley_tukey.c ../src/atfft/atfft_dft_bluestein.c

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
