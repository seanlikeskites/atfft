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

#include <stdlib.h>
#include <assert.h>
#include <atfft/dft_nd.h>
#include "atfft_internal.h"
#include "dft_pfa.h"

struct atfft_dft_pfa
{
    enum atfft_dft_algorithm algorithm;

    int dims [2];
    int size;
    enum atfft_direction direction;
    enum atfft_format format;
    struct atfft_dft_nd *fft;
    atfft_complex *sig, *dft;
    int *i_perm, *o_perm;
};

static void atfft_init_pfa_permutations (int size1,
                                         int size2,
                                         int size,
                                         int *i_perm,
                                         int *o_perm)
{
    /* input permutation */
    for (int n1 = 0; n1 < size1; ++n1)
    {
        int row_start = n1 * size2;

        for (int n2 = 0; n2 < size2; ++n2)
        {
            int i = atfft_mod (row_start + n2 * size1, size);

            i_perm [row_start + n2] = i;
        }
    }

    /* output permutation */
    int minv_1_2 = atfft_mult_inverse_mod_n (size1, size2);
    int minv_2_1 = atfft_mult_inverse_mod_n (size2, size1);

    for (int k1 = 0; k1 < size1; ++k1)
    {
        int row_start = k1 * size2;

        for (int k2 = 0; k2 < size2; ++k2)
        {
            int o = atfft_mod (row_start * minv_2_1 + k2 * size1 * minv_1_2, size);

            o_perm [o] = row_start + k2;
        }
    }
}

struct atfft_dft_pfa* atfft_dft_pfa_create (int size1,
                                            int size2,
                                            enum atfft_direction direction,
                                            enum atfft_format format)
{
    /* size1 and size2 must be coprime */
    int gcd = 1;
    atfft_gcd (size1, size2, &gcd, NULL, NULL);

    assert (gcd == 1);

    struct atfft_dft_pfa *fft;

    if (!(fft = calloc (1, sizeof (*fft))))
        return NULL;

    fft->algorithm = ATFFT_PFA;
    fft->dims [0] = size1;
    fft->dims [1] = size2;
    fft->size = size1 * size2;
    fft->direction = direction;
    fft->format = format;

    /* allocate a 2D dft structure for performing dft */
    fft->fft = atfft_dft_nd_create (fft->dims, 2, direction, format);

    if (!fft->fft)
        goto failed;

    /* allocate working buffers */
    fft->sig = malloc (fft->size * sizeof (*(fft->sig)));
    fft->dft = malloc (fft->size * sizeof (*(fft->dft)));

    if (!(fft->sig && fft->dft))
        goto failed;

    /* allocate permutation maps */
    fft->i_perm = malloc (fft->size * sizeof (*(fft->i_perm)));
    fft->o_perm = malloc (fft->size * sizeof (*(fft->o_perm)));

    if (!(fft->i_perm && fft->o_perm))
        goto failed;

    atfft_init_pfa_permutations (size1, size2, fft->size, fft->i_perm, fft->o_perm);

    return fft;

failed:
    atfft_dft_pfa_destroy (fft);
    return NULL;
}

void atfft_dft_pfa_destroy (void *fft)
{
    struct atfft_dft_pfa *t = fft;

    if (t)
    {
        free (t->o_perm);
        free (t->i_perm);
        free (t->dft);
        free (t->sig);
        atfft_dft_nd_destroy (t->fft);
        free (t);
    }
}

static void atfft_pfa_permute (int *perm,
                               int size,
                               atfft_complex *in,
                               int in_stride,
                               atfft_complex *out,
                               int out_stride)
{
    for (int i = 0; i < size; ++i)
    {
        atfft_copy_complex (in [in_stride * perm [i]], out + out_stride * i);
    }
}

void atfft_dft_pfa_complex_transform (void *fft,
                                      atfft_complex *in,
                                      int in_stride,
                                      atfft_complex *out,
                                      int out_stride)
{
    struct atfft_dft_pfa *t = fft;

    /* permute intput */
    atfft_pfa_permute (t->i_perm, t->size, in, in_stride, t->sig, 1);

    /* peform 2d dft */
    atfft_dft_nd_complex_transform (t->fft, t->sig, t->dft);

    /* permute output */
    atfft_pfa_permute (t->o_perm, t->size, t->dft, 1, out, out_stride);
}

cJSON* atfft_dft_pfa_get_plan (struct atfft_dft_pfa *fft)
{
    cJSON *alg = NULL;

    cJSON *plan_structure = cJSON_CreateObject();

    if (!plan_structure)
        goto failed;

    alg = cJSON_AddStringToObject (plan_structure, "Algorithm", "Prime Factor Algorithm");

    if (!alg)
        goto failed;

    return plan_structure;

failed:
    cJSON_Delete (plan_structure);
    return NULL;
}
