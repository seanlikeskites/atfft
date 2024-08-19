/*
 * Copyright (c) 2021 Sean Enderby <sean.enderby@gmail.com>
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

#include <math.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <atfft/dft.h>
#include "atfft_internal.h"
#include "dft_rader.h"
#include "dft_plan.h"

static int atfft_primitive_root_mod_n (int n)
{
    /* this algorithm will only work for prime numbers */
    if (!atfft_is_prime (n))
        return -1;

    /* loop over potential primitive roots */
    for (int g = 2; g < n; ++g)
    {
        int i = n - 2;
        int m = 1;
        int is_root = 1;

        while (i--)
        {
            m = (m * g) % n;

            if (m == 1)
            {
                is_root = 0;
                break;
            }
        }

        if (is_root)
            return g;
    }

    return -1;
}

static int atfft_rader_convolution_fft_size (int rader_size)
{
    if (atfft_is_power_of_2 (rader_size))
        return rader_size;
    else
        return atfft_next_power_of_2 (2 * rader_size - 1);
}

struct atfft_dft_rader
{
    enum atfft_dft_algorithm algorithm;

    int size;
    int rader_size;
    enum atfft_direction direction;
    enum atfft_format format;
    int p_root1, p_root2;
    int conv_size;
    struct atfft_dft *fft;
    int *perm1, *perm2;
    atfft_complex *sig, *sig_dft, *conv, *conv_dft;
};

static void atfft_init_rader_permutations (int *perm, int size, int p_root)
{
    int i = 1;

    for (int n = 0; n < size - 1; ++n)
    {
        perm [n] = i;
        i = atfft_mod (i * p_root, size);
    }
}

static int atfft_init_rader_convolution_dft (int size,
                                             enum atfft_direction direction,
                                             atfft_complex *conv_dft,
                                             int conv_size,
                                             int *perm,
                                             int perm_size,
                                             struct atfft_dft *fft)
{
    atfft_complex *t_factors = calloc (conv_size, sizeof (*t_factors));

    if (!t_factors)
        return -1;

    /* produce rader twiddle factors */
    for (int i = 0; i < perm_size; ++i)
    {
        atfft_scaled_twiddle_factor (perm [i], size, direction, conv_size, t_factors + i);
    }

    /* replicate samples for circular convolution */
    if (conv_size > perm_size)
    {
        size_t n_replications = perm_size - 1;
        atfft_complex *source = t_factors + 1;
        atfft_complex *dest = t_factors + conv_size - n_replications;

        memcpy (dest, source, n_replications * sizeof (*dest));
    }

    /* take DFT of twiddle factors */
    atfft_dft_complex_transform (fft, t_factors, conv_dft);

    free (t_factors);
    return 0;
}

struct atfft_dft_rader* atfft_dft_rader_create (int size,
                                                enum atfft_direction direction,
                                                enum atfft_format format)
{
    /* we can only find primitive roots for prime numbers */
    assert (atfft_is_prime (size));

    struct atfft_dft_rader *fft;

    if (!(fft = calloc (1, sizeof (*fft))))
        return NULL;

    fft->algorithm = ATFFT_RADER;
    fft->size = size;
    fft->rader_size = size - 1;
    fft->direction = direction;
    fft->format = format;
    fft->p_root1 = atfft_primitive_root_mod_n (size);
    fft->p_root2 = atfft_mult_inverse_mod_n (fft->p_root1, size);

    /* allocate some regular dft objects for performing the convolution */
    fft->conv_size = atfft_rader_convolution_fft_size (fft->rader_size);
    fft->fft = atfft_dft_create (fft->conv_size, ATFFT_FORWARD, ATFFT_COMPLEX);

    if (!fft->fft)
        goto failed;

    /* generate permutation indices */
    fft->perm1 = malloc (fft->rader_size * sizeof (*(fft->perm1)));
    fft->perm2 = malloc (fft->rader_size * sizeof (*(fft->perm2)));

    if (!(fft->perm1 && fft->perm2))
        goto failed;

    atfft_init_rader_permutations (fft->perm1, size, fft->p_root1);
    atfft_init_rader_permutations (fft->perm2, size, fft->p_root2);

    /* allocate work space for performing the dft */
    fft->sig = calloc (fft->conv_size, sizeof (*(fft->sig))); /* set up for zero padding */
    fft->sig_dft = malloc (fft->conv_size * sizeof (*(fft->sig_dft)));
    fft->conv = malloc (fft->conv_size * sizeof (*(fft->conv)));
    fft->conv_dft = malloc (fft->conv_size * sizeof (*(fft->conv_dft)));

    if (!(fft->sig && fft->sig_dft && fft->conv_dft))
        goto failed;

    /* calculate the convolution dft */
    if (atfft_init_rader_convolution_dft (size,
                                          direction,
                                          fft->conv_dft,
                                          fft->conv_size,
                                          fft->perm2,
                                          fft->rader_size,
                                          fft->fft) < 0)
        goto failed;

    return fft;

failed:
    atfft_dft_rader_destroy (fft);
    return NULL;
}

void atfft_dft_rader_destroy (void *fft)
{
    struct atfft_dft_rader *t = fft;

    if (t)
    {
        free (t->conv_dft);
        free (t->conv);
        free (t->sig_dft);
        free (t->sig);
        free (t->perm2);
        free (t->perm1);
        atfft_dft_destroy (t->fft);
        free (t);
    }
}

static void atfft_rader_permute_input (int *perm,
                                       int size,
                                       atfft_complex *in, 
                                       int in_stride,
                                       atfft_complex *out)
{
    for (int i = 0; i < size; ++i)
    {
        atfft_copy_complex (in [in_stride * perm [i]], out + i);
    }
}

static void atfft_rader_permute_output (int *perm,
                                        int size,
                                        atfft_complex *in, 
                                        atfft_complex *out,
                                        int out_stride)
{
    for (int i = 0; i < size; ++i)
    {
        atfft_swap_complex (in [i], out + (out_stride * perm [i]));
    }
}

void atfft_dft_rader_complex_transform (void *fft,
                                        atfft_complex *in,
                                        int in_stride,
                                        atfft_complex *out,
                                        int out_stride)
{
    struct atfft_dft_rader *t = fft;

    atfft_complex in0, out0;

    atfft_copy_complex (in [0], &in0);
    atfft_copy_complex (in0, &out0);

    atfft_rader_permute_input (t->perm1, 
                               t->rader_size,
                               in,
                               in_stride,
                               t->sig);

    atfft_dft_complex_transform (t->fft, t->sig, t->sig_dft);

    for (int i = 0; i < t->conv_size; ++i)
    {
        atfft_multiply_by_and_swap_complex (t->sig_dft + i, t->conv_dft [i]);

        atfft_sum_complex (out0, t->sig [i], &out0);
    }

    /* add DC component to DC bin */
    ATFFT_RE (t->sig_dft [0]) += ATFFT_IM (in0);
    ATFFT_IM (t->sig_dft [0]) += ATFFT_RE (in0);

    atfft_dft_complex_transform (t->fft, t->sig_dft, t->conv);

    atfft_rader_permute_output (t->perm2, 
                                t->rader_size,
                                t->conv,
                                out,
                                out_stride);

    atfft_copy_complex (out0, &out[0]);
}

cJSON* atfft_dft_rader_get_plan (struct atfft_dft_rader *fft)
{
    cJSON *alg = NULL,
          *size = NULL,
          *conv_size = NULL,
          *conv_transform = NULL;

    cJSON *plan_structure = cJSON_CreateObject();

    if (!plan_structure)
        goto failed;

    alg = cJSON_AddStringToObject (plan_structure, "Algorithm", "Rader");
    size = cJSON_AddNumberToObject (plan_structure, "Size", fft->size);
    conv_size = cJSON_AddNumberToObject (plan_structure, "Convolution Transform Size", fft->conv_size);

    if (!(alg && size && conv_size))
        goto failed;

    conv_transform = atfft_dft_get_plan (fft->fft);

    if (!conv_transform)
        goto failed;

    cJSON_AddItemToObject (plan_structure, "Convolution Transform", conv_transform);

    return plan_structure;

failed:
    cJSON_Delete (plan_structure);
    return NULL;
}
