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
#include <math.h>
#include <limits.h>
#include <atfft/dft.h>
#include <atfft/dft_util.h>
#include "atfft_internal.h"
#include "dft_cooley_tukey.h"
#include "dft_rader.h"
#include "dft_bluestein.h"
#include "dft_pfa.h"
#include "dft_plan.h"

static void atfft_init_even_real_sinusoids (atfft_complex *sinusoids,
                                            int sinusoids_size,
                                            int dft_size,
                                            enum atfft_direction direction)
{
    for (int i = 0; i < sinusoids_size; ++i)
    {
        atfft_twiddle_factor (i + 1, dft_size, direction, sinusoids + i);
    }
}

typedef void (*complex_transform_function) (void*, atfft_complex*, int, atfft_complex*, int);
typedef void (*fft_destroy_function) (void*);

struct atfft_dft
{
    enum atfft_dft_algorithm algorithm;

    int size;
    int internal_dft_size;
    enum atfft_direction direction;
    enum atfft_format format;

    void *fft;
    complex_transform_function complex_transform;
    fft_destroy_function fft_destroy;

    int is_even_real;
    atfft_complex *real_in, *real_out;
    atfft_complex *sinusoids;
};

struct atfft_dft* atfft_dft_create (int size, enum atfft_direction direction, enum atfft_format format)
{
    struct atfft_dft *fft;
    int sinusoids_size = 0;

    if (!(fft = calloc (1, sizeof (*fft))))
        return NULL;

    fft->algorithm = ATFFT_BASE;
    fft->size = size;
    fft->internal_dft_size = size;
    fft->direction = direction;
    fft->format = format;

    if (format == ATFFT_REAL)
    {
        /* Even length real transforms can be computed as a complex
         * transform of half the length. */
        if (atfft_is_even(size) && direction == ATFFT_FORWARD)
        {
            fft->is_even_real = 1;
            fft->internal_dft_size = size / 2;
            sinusoids_size = fft->internal_dft_size - 1;
        }
        else
        {
            fft->is_even_real = 0;
        }

        fft->real_in = malloc (fft->internal_dft_size * sizeof (*(fft->real_in)));
        fft->real_out = malloc (fft->internal_dft_size * sizeof (*(fft->real_out)));
        fft->sinusoids = malloc (sinusoids_size * sizeof (*(fft->sinusoids)));

        if (!(fft->real_in && fft->real_out && fft->sinusoids))
            goto failed;

        atfft_init_even_real_sinusoids (fft->sinusoids,
                                        sinusoids_size,
                                        size,
                                        direction);
    }

    if (atfft_is_prime (fft->internal_dft_size))
    {
        if (atfft_is_power_of_2 (fft->internal_dft_size - 1))
        {
            /* Use Rader's algorithm */
            fft->fft = atfft_dft_rader_create (fft->internal_dft_size, direction, ATFFT_COMPLEX);
            fft->complex_transform = atfft_dft_rader_complex_transform;
            fft->fft_destroy = atfft_dft_rader_destroy;
        }
        else
        {
            /* Use Bluestein's algorithm */
            fft->fft = atfft_dft_bluestein_create (fft->internal_dft_size, direction, ATFFT_COMPLEX);
            fft->complex_transform = atfft_dft_bluestein_complex_transform;
            fft->fft_destroy = atfft_dft_bluestein_destroy;
        }
    }
    else
    {
        /* Use Cooley-Tukey */
        fft->fft = atfft_dft_cooley_tukey_create (fft->internal_dft_size, direction, ATFFT_COMPLEX);
        fft->complex_transform = atfft_dft_cooley_tukey_complex_transform;
        fft->fft_destroy = atfft_dft_cooley_tukey_destroy;
    }

    if (!fft->fft)
        goto failed;

    return fft;

failed:
    atfft_dft_destroy (fft);
    return NULL;
}

void atfft_dft_destroy (struct atfft_dft *fft)
{
    if (fft)
    {
        free (fft->sinusoids);
        free (fft->real_out);
        free (fft->real_in);
        fft->fft_destroy (fft->fft);
        free (fft);
    }
}

void atfft_dft_complex_transform (struct atfft_dft *fft,
                                  atfft_complex *in,
                                  atfft_complex *out)
{
    atfft_dft_complex_transform_stride (fft, in, 1, out, 1);
}

void atfft_dft_complex_transform_stride (struct atfft_dft *fft,
                                         atfft_complex *in,
                                         int in_stride,
                                         atfft_complex *out,
                                         int out_stride)
{
    /* Only to be used with complex FFTs. */
    assert (fft->format == ATFFT_COMPLEX);

    fft->complex_transform (fft->fft, in, in_stride, out, out_stride);
}

static void atfft_dft_even_real_forward_transform (struct atfft_dft *fft,
                                                   const atfft_sample *in,
                                                   int in_stride,
                                                   atfft_complex *out,
                                                   int out_stride)
{
    for (int i = 0; i < fft->size / 2; ++i)
    {
        ATFFT_RE (fft->real_in [i]) = in [2 * i * in_stride];
        ATFFT_IM (fft->real_in [i]) = in [(2 * i + 1) * in_stride];
    }

    fft->complex_transform (fft->fft, fft->real_in, 1, fft->real_out, 1);

    ATFFT_RE (out [0]) = ATFFT_RE (fft->real_out [0]) + ATFFT_IM (fft->real_out [0]);

    for (int i = 1; i < fft->internal_dft_size; ++i)
    {
        atfft_complex E, O;

        ATFFT_RE (E) = (ATFFT_RE (fft->real_out [i]) +
                        ATFFT_RE (fft->real_out [fft->internal_dft_size - i])) / 2;
        ATFFT_IM (E) = (ATFFT_IM (fft->real_out [i]) -
                        ATFFT_IM (fft->real_out [fft->internal_dft_size - i])) / 2;

        ATFFT_RE (O) = (ATFFT_IM (fft->real_out [i]) +
                        ATFFT_IM (fft->real_out [fft->internal_dft_size - i])) / 2;
        ATFFT_IM (O) = (ATFFT_RE (fft->real_out [fft->internal_dft_size - i]) -
                        ATFFT_RE (fft->real_out [i])) / 2;

        atfft_multiply_by_complex (&O, fft->sinusoids [i - 1]);
        atfft_sum_complex (E, O, out + i);
    }

    ATFFT_RE (out [fft->size / 2]) = ATFFT_RE (fft->real_out [0]) - ATFFT_IM (fft->real_out [0]);
}

static void atfft_dft_trivial_real_forward_transform (struct atfft_dft *fft,
                                                      const atfft_sample *in,
                                                      int in_stride,
                                                      atfft_complex *out,
                                                      int out_stride)
{
    atfft_real_to_complex_stride (in, in_stride, fft->real_in, 1, fft->size);
    fft->complex_transform (fft->fft, fft->real_in, 1, fft->real_out, 1);
    atfft_complex_to_halfcomplex_stride (fft->real_out, 1, out, out_stride, fft->size);
}

static void atfft_perform_real_forward_dft (struct atfft_dft *fft,
                                            const atfft_sample *in,
                                            int in_stride,
                                            atfft_complex *out,
                                            int out_stride)
{
    if (fft->is_even_real)
        atfft_dft_even_real_forward_transform (fft, in, in_stride, out, out_stride);
    else
        atfft_dft_trivial_real_forward_transform (fft, in, in_stride, out, out_stride);

}

void atfft_dft_real_forward_transform (struct atfft_dft *fft, const atfft_sample *in, atfft_complex *out)
{
    /* Only to be used for forward real FFTs. */
    assert ((fft->format == ATFFT_REAL) && (fft->direction == ATFFT_FORWARD));

    atfft_perform_real_forward_dft (fft, in, 1, out, 1);
}

void atfft_dft_real_forward_transform_stride (struct atfft_dft *fft,
                                              const atfft_sample *in,
                                              int in_stride,
                                              atfft_complex *out,
                                              int out_stride)
{
    /* Only to be used for forward real FFTs. */
    assert ((fft->format == ATFFT_REAL) && (fft->direction == ATFFT_FORWARD));

    atfft_perform_real_forward_dft (fft, in, in_stride, out, out_stride);
}

void atfft_dft_real_backward_transform (struct atfft_dft *fft, atfft_complex *in, atfft_sample *out)
{
    /* Only to be used for backward real FFTs. */
    assert ((fft->format == ATFFT_REAL) && (fft->direction == ATFFT_BACKWARD));

    atfft_halfcomplex_to_complex (in, fft->real_in, fft->size);
    fft->complex_transform (fft->fft, fft->real_in, 1, fft->real_out, 1);
    atfft_real (fft->real_out, out, fft->size);
}

void atfft_dft_real_backward_transform_stride (struct atfft_dft *fft,
                                               atfft_complex *in,
                                               int in_stride,
                                               atfft_sample *out,
                                               int out_stride)
{
    /* Only to be used for backward real FFTs. */
    assert ((fft->format == ATFFT_REAL) && (fft->direction == ATFFT_BACKWARD));

    atfft_halfcomplex_to_complex_stride (in, in_stride, fft->real_in, 1, fft->size);
    fft->complex_transform (fft->fft, fft->real_in, 1, fft->real_out, 1);
    atfft_real_stride (fft->real_out, 1, out, out_stride, fft->size);
}

cJSON* atfft_dft_base_get_plan (struct atfft_dft *fft)
{
    cJSON *alg = NULL,
          *size = NULL,
          *direction = NULL,
          *format = NULL,
          *internal_transform = NULL;

    cJSON *plan_structure = cJSON_CreateObject();

    if (!plan_structure)
        goto failed;

    alg = cJSON_AddStringToObject (plan_structure, "Algorithm", "atfft Base Transform");
    size = cJSON_AddNumberToObject (plan_structure, "Size", fft->size);
    direction = cJSON_AddStringToObject (plan_structure, "Direction",
                                         fft->direction == ATFFT_FORWARD ? "forward" : "backward");
    format = cJSON_AddStringToObject (plan_structure, "Format",
                                      fft->format == ATFFT_COMPLEX ? "complex" : "real");

    if (!(alg && size && direction && format))
        goto failed;

    internal_transform = atfft_dft_get_plan (fft->fft);

    if (!internal_transform)
        goto failed;

    cJSON_AddItemToObject (plan_structure, "Internal Transform", internal_transform);
    
    return plan_structure;

failed:
    cJSON_Delete (plan_structure);
    return NULL;
}

void atfft_dft_print_plan (struct atfft_dft *fft, FILE *stream)
{
    char *plan_json = NULL;
    cJSON* plan_structure = atfft_dft_get_plan (fft);

    if (!plan_structure)
        goto failed;

    plan_json = cJSON_Print (plan_structure);

    if (!plan_json)
        goto failed;

    fprintf (stream, "%s\n", plan_json);

    goto succeeded;

failed:
    fprintf (stream, "{\"Error\": \"Failed to serialise plan.\"}\n");

succeeded:
    cJSON_Delete (plan_structure);
    free (plan_json);
}
