#include <atfft/atfft.h>
#define N_SAMPLES 32 /* Define the length of signals we are working with. */

int main()
{
    /* For a DFT on complex valued data, both the time and frequency domain signals
     * are complex valued and of the same length.
     */
    atfft_complex time_domain [N_SAMPLES];
    atfft_complex freq_domain [N_SAMPLES];

    /* Here we would populate the timeDomain array with whatever signal we wish to
     * find the DFT of.
     */

    /* Create a DFT plan, specifying the length of the time domain signal, that we want
     * to do a forward transform (from time to frequency domain), and the fact that our
     * time domain signal is complex valued.
     */
    struct atfft_dft *fft = atfft_dft_create (N_SAMPLES, ATFFT_FORWARD, ATFFT_COMPLEX);

    /* Use the plan to perform a DFT on our signal. */
    atfft_dft_complex_transform (fft, time_domain, freq_domain);

    /* Now the freq_domain array contains the DFT of the signal we put in the time_domain
     * array.
     */

    /* Free the plan once we are done with it. */
    atfft_dft_destroy (fft);

    return 0;
}
