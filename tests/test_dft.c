#include <atfft/atfft.h>
#include <stdlib.h>
#include <math.h>
#include <stdio.h>

void generate_real_dc (atfft_sample *sig, int size, atfft_sample offset)
{
    for (int i = 0; i < size; ++i)
    {
        sig [i] = offset;
    }
}

void generate_complex_dc (atfft_complex *sig, int size, atfft_complex offset)
{
    for (int i = 0; i < size; ++i)
    {
        ATFFT_RE (sig [i]) = ATFFT_RE (offset);
        ATFFT_IM (sig [i]) = ATFFT_IM (offset);
    }
}

void generate_real_impulse (atfft_sample *sig, int size)
{
    sig [0] = 1.0;

    for (int i = 1; i < size; ++i)
    {
        sig [i] = 0.0;
    }
}

void generate_complex_impulse (atfft_complex *sig, int size)
{
    generate_real_impulse ((atfft_sample*) sig, size * 2);
}

atfft_sample max_error_real (const atfft_sample *a, const atfft_sample *b, int size)
{
    atfft_sample error = 0.0;

    for (int i = 0; i < size; ++i)
    {
        atfft_sample e = fabs (a [i] - b [i]);

        if (e > error)
            error = e;
    }

    return error;
}

atfft_sample max_error_complex (atfft_complex *a, atfft_complex *b, int size)
{
    return max_error_real ((atfft_sample*) a, (atfft_sample*) b, size * 2);
}

int main()
{
    int ret = 0;
    int n_samples = 32;
    atfft_complex offset = {1.0, 0.0};

    atfft_complex *impulse = malloc (n_samples * sizeof (*impulse));
    atfft_complex *impulse_dft = malloc (n_samples * sizeof (*impulse_dft));
    atfft_complex *dc = malloc (n_samples * sizeof (*dc));

    if (!(impulse && impulse_dft && dc))
    {
        ret = 1;
        goto cleanup;
    }

    generate_complex_impulse (impulse, n_samples);
    generate_complex_dc (dc, n_samples, offset);

    struct atfft_dft *fft_forward = atfft_dft_create (n_samples, ATFFT_FORWARD, ATFFT_COMPLEX);

    if (!fft_forward)
    {
        ret = 1;
        goto cleanup;
    }

    atfft_dft_complex_transform (fft_forward, impulse, impulse_dft);

    atfft_sample error = max_error_complex (impulse_dft, dc, n_samples);

    printf ("Impulse error: %f\n", error);

    if (error > 10e-10)
    {
        ret = 1;
        goto cleanup;
    }

cleanup:
    atfft_dft_destroy (fft_forward);
    free (dc);
    free (impulse_dft);
    free (impulse);

    return ret;
}
