#include <cstdlib>
#include <ctime>
#include <cmath>
#include <benchmark/benchmark.h>
#include <atfft/atfft.h>

float randomFloat()
{
    return (float) rand() / (float) RAND_MAX;
}

static void complexTransform (benchmark::State &state)
{
    int size = state.range (0);
    atfft_complex *x = new atfft_complex [size];
    atfft_complex *y = new atfft_complex [size];
    atfft_dft *fft = static_cast <atfft_dft*> (atfft_dft_create (size, ATFFT_FORWARD, ATFFT_COMPLEX));

    srand (time (NULL));

    for (int i = 0; i < size; ++i)
    {
        ATFFT_RE (x [i]) = randomFloat();
        ATFFT_IM (x [i]) = randomFloat();
    }

    for (auto _ : state)
    {
        atfft_dft_complex_transform (fft, x, y);
    }

    atfft_dft_destroy (fft);
    delete[] y;
    delete[] x;
}

static void realTransform (benchmark::State &state)
{
    int size = state.range (0);
    atfft_sample *x = new atfft_sample [size];
    atfft_complex *y = new atfft_complex [atfft_halfcomplex_size(size)];
    atfft_dft *fft = atfft_dft_create (size, ATFFT_FORWARD, ATFFT_REAL);

    srand (time (NULL));

    for (int i = 0; i < size; ++i)
    {
        x [i] = randomFloat();
    }

    for (auto _ : state)
    {
        atfft_dft_real_forward_transform (fft, x, y);
    }

    atfft_dft_destroy (fft);
    delete[] y;
    delete[] x;
}

#ifdef BENCHMARK_POWERS_OF_2
BENCHMARK (realTransform)->RangeMultiplier (2)->Range (1 << 5, 1 << 16);
BENCHMARK (complexTransform)->RangeMultiplier (2)->Range (1 << 5, 1 << 16);
#endif

#ifdef BENCHMARK_POWERS_OF_3
BENCHMARK (realTransform)->RangeMultiplier (3)->Range (27, std::pow (3, 10));
BENCHMARK (complexTransform)->RangeMultiplier (3)->Range (27, std::pow (3, 10));
#endif

#ifdef BENCHMARK_COMPOSITES
BENCHMARK (realTransform)->RangeMultiplier (30)->Range (30, std::pow (30, 4));
BENCHMARK (complexTransform)->RangeMultiplier (30)->Range (30, std::pow (30, 4));
#endif

#ifdef BENCHMARK_PRIMES
//BENCHMARK (realTransform)->Arg (37)
//                         ->Arg (67)
//                         ->Arg (131)
//                         ->Arg (257)
//                         ->Arg (521)
//                         ->Arg (1031)
//                         ->Arg (2053)
//                         ->Arg (4099)
//                         ->Arg (8191)
//                         ->Arg (16381)
//                         ->Arg (32771)
//                         ->Arg (65537);
BENCHMARK (complexTransform)->Arg (37)
                            ->Arg (67)
                            ->Arg (131)
                            ->Arg (257)
                            ->Arg (521)
                            ->Arg (1031)
                            ->Arg (2053)
                            ->Arg (4099)
                            ->Arg (8191)
                            ->Arg (16381)
                            ->Arg (32771)
                            ->Arg (65537);
#endif

BENCHMARK_MAIN();
