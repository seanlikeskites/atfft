#include <cstdlib>
#include <ctime>
#include <benchmark/benchmark.h>
#include <atfft.h>

float randomFloat()
{
    return (float) rand() / (float) RAND_MAX;
}

static void complexTransform (benchmark::State &state)
{
    int size = state.range_x();
    atfft_complex *x = new atfft_complex [size];
    atfft_complex *y = new atfft_complex [size];
    atfft *fft = static_cast <atfft*> (atfft_create (size, ATFFT_FORWARD, ATFFT_COMPLEX));

    srand (time (NULL));

    for (int i = 0; i < size; ++i)
    {
        ATFFT_REAL (x [i]) = randomFloat();
        ATFFT_IMAG (x [i]) = randomFloat();
    }

    while (state.KeepRunning())
    {
        atfft_complex_transform (fft, x, y);
    }

    atfft_destroy (fft);
    delete[] y;
    delete[] x;
}

BENCHMARK (complexTransform)->Range (16, 1 << 16);

BENCHMARK_MAIN();
