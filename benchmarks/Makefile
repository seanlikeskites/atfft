CFLAGS := -Wall -O3
CXXFLAGS := -Wall -O3 -std=c++11 -DBENCHMARK_POWERS_OF_2
LIBS := -lm -lpthread -lbenchmark
INCLUDEDIRS := -I../include
FFT_SOURCES := ../src/atfft/atfft_shared.c

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

	FFT_SOURCES += ../src/fftw/dft.c
	CXXFLAGS += -DBENCHMARK_POWERS_OF_3 -DBENCHMARK_PRIMES
	
else ifeq (FFTS, $(IMPL))
	LIBS += -lffts
	FFT_SOURCES += ../src/ffts/dft.c
	CXXFLAGS += -DBENCHMARK_POWERS_OF_3 -DBENCHMARK_PRIMES

else ifeq (GSL, $(IMPL))
	LIBS += -lgsl -lgslcblas
	FFT_SOURCES += ../src/gsl/dft.c
	CXXFLAGS += -DBENCHMARK_POWERS_OF_3 -DBENCHMARK_PRIMES

else ifeq (KISS, $(IMPL))
	ifeq (FLOAT, $(TYPE))
		LIBS += -lkissfftf
	else ifeq (LONG_DOUBLE, $(TYPE))
		LIBS += -lkissfftl
	else
		LIBS += -lkissfft
	endif

	FFT_SOURCES += ../src/kiss/dft.c
	CXXFLAGS += -DBENCHMARK_POWERS_OF_3 -DBENCHMARK_PRIMES
	
else ifeq (VDSP, $(IMPL))
	LIBS += -framework Accelerate
	FFT_SOURCES += ../src/vdsp/dft.c

else ifeq (FFMPEG, $(IMPL))
	LIBS += -lavutil -lavcodec
	FFT_SOURCES += ../src/ffmpeg/dft.c

else ifeq (PFFFT, $(IMPL))
	FFT_SOURCES += ../src/pffft/dft.c ../src/pffft/pffft.c
	CFLAGS += -msse

else ifeq (OOURA, $(IMPL))
	FFT_SOURCES += ../src/ooura/dft.c ../src/ooura/fft4g.c

else
	FFT_SOURCES += ../src/atfft/dft.c ../src/atfft/dft_rader.c ../src/atfft/atfft_internal.c ../src/atfft/dft_cooley_tukey.c ../src/atfft/dft_bluestein.c
	CXXFLAGS += -DBENCHMARK_POWERS_OF_3 -DBENCHMARK_PRIMES
endif

FFT_TEMP1 = $(notdir $(FFT_SOURCES))
FFT_OBJECTS = $(subst .c,.o, $(FFT_TEMP1))

benchmark: benchmark.cpp $(FFT_SOURCES)
	$(CC) -c $(FFT_SOURCES) $(INCLUDEDIRS) $(CFLAGS)
	$(CXX) -o $@ $< $ $(FFT_OBJECTS) $(INCLUDEDIRS) $(LIBS) $(CXXFLAGS)

clean:
	rm -rf *.o benchmark
