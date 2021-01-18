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
#include <atfft/dct.h>
#include "ipp_definitions.h"

struct atfft_dct
{
    int size;
    enum atfft_direction direction;

    /* the ipp plan and work area */
    void *plan;
    Ipp8u *plan_buffer;
    Ipp8u *work_area;

    /* scale factors and scaled buffer */
    atfft_sample factor0, factor1;
    atfft_ipp_sample *scaled_buffer;
};

int atfft_dct_is_supported_size (int size)
{
    return size > 0;
}

static void* init_ipp_plan (int size,
                            enum atfft_direction direction,
                            Ipp8u **plan_buffer,
                            Ipp8u **work_area)
{
    atfft_init_ipp();

    IppStatus status = ippStsNoErr;

    /* get the sizes of all the buffers we need */
    int plan_buffer_size = 0;
    int plan_init_size = 0;
    int work_area_size = 0;

    if (direction == ATFFT_FORWARD)
        status = ATFFT_IPPS_DCT_FWD_GET_SIZE (size,
                                              ippAlgHintNone,
                                              &plan_buffer_size,
                                              &plan_init_size,
                                              &work_area_size);
    else
        status = ATFFT_IPPS_DCT_INV_GET_SIZE (size,
                                              ippAlgHintNone,
                                              &plan_buffer_size,
                                              &plan_init_size,
                                              &work_area_size);

    if (status != ippStsNoErr)
        return NULL;

    /* allocate buffers */
    Ipp8u *plan_init_buffer = NULL;

    *plan_buffer = ippMalloc (plan_buffer_size);

    if (!(*plan_buffer))
        goto failed;

    if (plan_init_size > 0)
    {
        plan_init_buffer = ippMalloc (plan_init_size);

        if (!plan_init_buffer)
            goto failed;
    }

    if (work_area_size > 0)
    {
        *work_area = ippMalloc (work_area_size);

        if (!(*work_area))
            goto failed;
    }

    void *plan = NULL;

    /* initialise the ipp plan */
    if (direction == ATFFT_FORWARD)
        status = ATFFT_IPPS_DCT_FWD_INIT ((ATFFT_IPPS_DCT_FWD_SPEC**) &plan,
                                          size,
                                          ippAlgHintNone,
                                          *plan_buffer,
                                          plan_init_buffer);
    else
        status = ATFFT_IPPS_DCT_INV_INIT ((ATFFT_IPPS_DCT_INV_SPEC**) &plan,
                                          size,
                                          ippAlgHintNone,
                                          *plan_buffer,
                                          plan_init_buffer);

    if (status != ippStsNoErr)
        goto failed;

    ippFree (plan_init_buffer);

    return plan;

 failed:
    ippFree (*work_area);
    *work_area = NULL;

    ippFree (plan_init_buffer);

    ippFree (*plan_buffer);
    *plan_buffer = NULL;

    return NULL;
}

struct atfft_dct* atfft_dct_create (int size, enum atfft_direction direction)
{
    /* ipp supports all sizes. */
    assert (atfft_dct_is_supported_size (size));

    struct atfft_dct *plan;

    if (!(plan = calloc (1, sizeof (*plan))))
        return NULL;

    plan->size = size;
    plan->direction = direction;

    plan->plan = init_ipp_plan (size,
                                direction,
                                &(plan->plan_buffer),
                                &(plan->work_area));

    if (!plan->plan)
        goto failed;

    /* calculate scale factors and allocate the scaled buffer */
    if (direction == ATFFT_FORWARD)
        plan->factor0 = sqrt (size);
    else
        plan->factor0 = sqrt (size) / 2.0;

    plan->factor1 = sqrt(size / 2.0);
    plan->scaled_buffer = ippMalloc (size * sizeof (*(plan->scaled_buffer)));

    if (!plan->scaled_buffer)
        goto failed;

    return plan;

failed:
    atfft_dct_destroy (plan);
    return NULL;
}

void atfft_dct_destroy (struct atfft_dct *plan)
{
    if (plan)
    {
        ippFree (plan->scaled_buffer);
        ippFree (plan->work_area);
        ippFree (plan->plan);
        free (plan);
    }
}

static void scale_forward_transform (const atfft_ipp_sample* in,
                                     atfft_sample *out,
                                     int size,
                                     atfft_sample factor0,
                                     atfft_sample factor1)
{
    out [0] = factor0 * in [0];

    for (int i = 1; i < size; ++i)
    {
        out [i] = factor1 * in [i];
    }
}

static void scale_backward_transform (const atfft_sample* in,
                                      atfft_ipp_sample *out,
                                      int size,
                                      atfft_sample factor0,
                                      atfft_sample factor1)
{
    out [0] = factor0 * in [0];

    for (int i = 1; i < size; ++i)
    {
        out [i] = factor1 * in [i];
    }
}

void atfft_dct_transform (struct atfft_dct *plan, const atfft_sample *in, atfft_sample *out)
{
    if (plan->direction == ATFFT_FORWARD)
    {
        ATFFT_IPPS_DCT_FWD (in, plan->scaled_buffer, plan->plan, plan->work_area);
        scale_forward_transform (plan->scaled_buffer,
                                 out,
                                 plan->size,
                                 plan->factor0,
                                 plan->factor1);


    }
    else
    {
        scale_backward_transform (in,
                                  plan->scaled_buffer,
                                  plan->size,
                                  plan->factor0,
                                  plan->factor1);

        ATFFT_IPPS_DCT_INV (plan->scaled_buffer, out, plan->plan, plan->work_area);
    }

}
