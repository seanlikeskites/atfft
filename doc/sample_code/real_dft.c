#include <atfft/atfft.h>
#define N_SAMPLES 32 /* Define the length of signals we are working with. */

int main()
{
    /* For a DFT on real valued data our time domain signal is real valued. */
    atfft_sample time_domain [N_SAMPLES];

    /* For a real valued time domain signal ATFFT will produce a halfcomplex
     * frequency domain signal. This is a complex valued signal with length
     * given by atfft_halfcomplex_size().
     */
    int freq_domain_length = atfft_halfcomplex_size (N_SAMPLES);
    atfft_complex *freq_domain = malloc (freq_domain_length * sizeof (*freq_domain));

    /* Here we would populate the timeDomain array with whatever signal we wish to
     * find the DFT of.
     */

    /* Create a DFT plan, specifying the length of the time domain signal, that we want
     * to do a forward transform (from time to frequency domain), and the fact that our
     * time domain signal is real valued.
     */
    struct atfft_dft *fft = atfft_dft_create (N_SAMPLES, ATFFT_FORWARD, ATFFT_REAL);

    /* Use the plan to perform a DFT on our signal. */
    atfft_dft_real_forward_transform (fft, time_domain, freq_domain);

    /* Now the freq_domain array contains the DFT of the signal we put in the time_domain
     * array.
     */

    /* Free the plan once we are done with it. */
    atfft_dft_destroy (fft);
    free (freq_domain);

    return 0;
}
