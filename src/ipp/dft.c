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
#include <atfft/dft.h>
#include <atfft/dft_util.h>
#include <ipp.h>
#include "ipp_definitions.h"

static void atfft_init_ipp()
{
    static int ipp_initialised = 0;

    if (!ipp_initialised)
    {
        ippInit();
    }
}

struct atfft_dft
{
    int size;
    enum atfft_direction direction;
    enum atfft_format format;

    /* the ipp plan and work area */
    void *plan;
    Ipp8u *work_area;
};

int atfft_dft_is_supported_size (int size, enum atfft_format format)
{
#ifdef ATFFT_TYPE_FLOAT
    int max_complex_size = (1 << 27) - 1;
#else
    int max_complex_size = (1 << 26) - 1;
#endif

    if (format == ATFFT_COMPLEX)
        return size > 0 && size <= max_complex_size;
    else
        return size > 0;
}

static void* init_ipp_plan (int size,
                            enum atfft_format format,
                            Ipp8u **work_area)
{
    atfft_init_ipp();

    IppStatus status = ippStsNoErr;

    /* get the sizes of all the buffers we need */
    int plan_buffer_size = 0;
    int plan_init_size = 0;
    int work_area_size = 0;

    if (format == ATFFT_COMPLEX)
        status = ATFFT_IPPS_DFT_GET_SIZE_C (size,
                                            IPP_FFT_NODIV_BY_ANY,
                                            ippAlgHintNone,
                                            &plan_buffer_size,
                                            &plan_init_size,
                                            &work_area_size);
    else
        status = ATFFT_IPPS_DFT_GET_SIZE_R (size,
                                            IPP_FFT_NODIV_BY_ANY,
                                            ippAlgHintNone,
                                            &plan_buffer_size,
                                            &plan_init_size,
                                            &work_area_size);

    if (status != ippStsNoErr)
        return NULL;

    /* allocate buffers */
    void *plan = ippMalloc (plan_buffer_size);
    Ipp8u *plan_init_buffer = NULL;

    if (!plan)
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

    /* initialise the ipp plan */
    if (format == ATFFT_COMPLEX)
        status = ATFFT_IPPS_DFT_INIT_C (size,
                                        IPP_FFT_NODIV_BY_ANY,
                                        ippAlgHintNone,
                                        (ATFFT_IPPS_DFT_SPEC_C*) plan,
                                        plan_init_buffer);
    else
        status = ATFFT_IPPS_DFT_INIT_R (size,
                                        IPP_FFT_NODIV_BY_ANY,
                                        ippAlgHintNone,
                                        (ATFFT_IPPS_DFT_SPEC_R*) plan,
                                        plan_init_buffer);

    if (status != ippStsNoErr)
        goto failed;

    ippFree (plan_init_buffer);

    return plan;

 failed:
    ippFree (*work_area);
    *work_area = NULL;

    ippFree (plan_init_buffer);
    ippFree (plan);

    return NULL;
}

struct atfft_dft* atfft_dft_create (int size, enum atfft_direction direction, enum atfft_format format)
{
    /* ipp supports different lengths depending on a variety of parameters. */
    assert (atfft_dft_is_supported_size (size, format));

    struct atfft_dft *plan;

    if (!(plan = calloc (1, sizeof (*plan))))
        return NULL;

    plan->size = size;
    plan->direction = direction;
    plan->format = format;

    plan->plan = init_ipp_plan (size,
                                format,
                                &(plan->work_area));

    if (!plan->plan)
        goto failed;

    return plan;

failed:
    atfft_dft_destroy (plan);
    return NULL;
}

void atfft_dft_destroy (struct atfft_dft *plan)
{
    if (plan)
    {
        ippFree (plan->work_area);
        ippFree (plan->plan);
        free (plan);
    }
}

void atfft_dft_complex_transform (struct atfft_dft *plan, atfft_complex *in, atfft_complex *out)
{
    /* Only to be used with complex FFTs. */
    assert (plan->format == ATFFT_COMPLEX);

    if (plan->direction == ATFFT_FORWARD)
        ATFFT_IPPS_DFT_FWD_CTOC ((const atfft_ipp_complex*) in, (atfft_ipp_complex*) out, plan->plan, plan->work_area);
    else
        ATFFT_IPPS_DFT_INV_CTOC ((const atfft_ipp_complex*) in, (atfft_ipp_complex*) out, plan->plan, plan->work_area);
}

void atfft_dft_real_forward_transform (struct atfft_dft *plan, const atfft_sample *in, atfft_complex *out)
{
    /* Only to be used for forward real FFTs. */
    assert ((plan->format == ATFFT_REAL) && (plan->direction == ATFFT_FORWARD));

    ATFFT_IPPS_DFT_FWD_RTOCCS ((const atfft_ipp_sample*) in, (atfft_ipp_sample*) out, plan->plan, plan->work_area);
}

void atfft_dft_real_backward_transform (struct atfft_dft *plan, atfft_complex *in, atfft_sample *out)
{
    /* Only to be used for backward real FFTs. */
    assert ((plan->format == ATFFT_REAL) && (plan->direction == ATFFT_BACKWARD));

    ATFFT_IPPS_DFT_INV_CCSTOR ((const atfft_ipp_sample*) in, (atfft_ipp_sample*) out, plan->plan, plan->work_area);
}
