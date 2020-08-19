#include <cstdlib>
#include <ctime>
#include <benchmark/benchmark.h>
#include <atfft/atfft.h>

float randomFloat()
{
    return (float) rand() / (float) RAND_MAX;
}

static void complexTransform (benchmark::State &state)
{
    int size = state.range_x();
    atfft_complex *x = new atfft_complex [size];
    atfft_complex *y = new atfft_complex [size];
    atfft_dft *fft = static_cast <atfft_dft*> (atfft_dft_create (size, ATFFT_FORWARD, ATFFT_COMPLEX));

    srand (time (NULL));

    for (int i = 0; i < size; ++i)
    {
        ATFFT_RE (x [i]) = randomFloat();
        ATFFT_IM (x [i]) = randomFloat();
    }

    while (state.KeepRunning())
    {
        atfft_dft_complex_transform (fft, x, y);
    }

    atfft_dft_destroy (fft);
    delete[] y;
    delete[] x;
}

//BENCHMARK (complexTransform)->Arg(32)->Arg(64)->Arg(128)
//                            ->Arg(256)->Arg(512)->Arg(1024)->Arg(2048)->Arg(4096)->Arg(8192)
//                            ->Arg(16384)->Arg(32768)->Arg(65536);
BENCHMARK (complexTransform)->Arg(5)->Arg(7)->Arg(13)->Arg(31)->Arg(67)->Arg(127)
                            ->Arg(257)->Arg(509)->Arg(1021)->Arg(2053)->Arg(35987);

BENCHMARK_MAIN();
