/*
  This file is part of the TGV package (https://github.com/chixindebaoyu/tgvnn).

  The MIT License (MIT)

  Copyright (c) Dong Wang

  Permission is hereby granted, free of charge, to any person obtaining a # copy
  of this software and associated documentation files (the "Software"), to # deal
  in the Software without restriction, including without limitation the # rights
  to use, copy, modify, merge, publish, distribute, sublicense, and/or # sell
  copies of the Software, and to permit persons to whom the Software is
  furnished to do so, subject to the following conditions:

  The above copyright notice and this permission notice shall be included # in all
  copies or substantial portions of the Software.

  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS # OR
  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL # THE
  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING # FROM,
  OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS # IN THE
  SOFTWARE.
*/

#ifndef _TGV_H
#define _TGV_H

#include <stdbool.h>
#include <stdlib.h>
#include <stdio.h>
#include <fcntl.h>
#include <unistd.h>
#include <err.h>
#include <errno.h>
#include <getopt.h>
#include <string.h>
#include <math.h>
#include <complex.h>
#include <time.h>
#include <stdint.h>

#include <cufft.h>
#include <cusolverDn.h>
#include <cuda_runtime.h>
#include <cublas_v2.h>

#include "float2math.h"
#include "ra.h"

#ifdef __cplusplus
#endif

void
gpuAssert(cudaError_t code, const char *file, int line);
void *
mallocAssert(const size_t n, const char *file, int line);
void
freeAssert(void *p, const char *file, int line);
static const char
*_cudaGetErrorEnum(cufftResult error);
void
print_usage();
void
save_rafile(float2 *h_out, const char *out_path,
  const size_t dim0, const size_t dim1, const size_t dim2, const size_t dim3);

__global__ void
scaledata(float2 *d_array, const size_t array_size, const float factor);
float
compute_maxmag(float2 *d_array, const size_t array_size);

void
fft2_init(const int rows, const int cols, const int ndyn);
void
forward(float2 *d_out, float2 *d_in, float2 *d_mask,
  const size_t N, const size_t rows, const size_t cols);
void
backward(float2 *d_out, float2 *d_in, float2 *d_mask,
  const size_t N, const size_t rows, const size_t cols);

__global__ void
arrayadd (float2 *array_c, float2 *array_a, float2 *array_b,
  const size_t array_size, float alpha, float beta);
__global__ void
arrayadd (float2 *array_c, float2 *array_a, float2 *array_b,
  const size_t array_size, float alpha, float beta);
__global__ void
arraydot (float2 *array_c, float2 *array_a, float2 *array_b,
  const size_t array_size, float alpha);
__global__ void
arrayreal (float2 *d_out, float2 *d_in, size_t array_size);

__global__ void
shrink (float *d_array, const float beta, const int array_size);
float
compute_alpha (float alpha0, float alpha1, int iter, int index);

__global__ void
grad_xx (float2 *d_out, float2 *d_in,
  const size_t N, const size_t rows, const size_t cols, const size_t ndyn);
__global__ void
grad_yy (float2 *d_out, float2 *d_in,
  const size_t N, const size_t rows, const size_t cols, const size_t ndyn);
__global__ void
grad_tt (float2 *d_out, float2 *d_in,
  const size_t N, const size_t rows, const size_t cols, const size_t ndyn,
  float mu);

__global__ void
grad_xx_bound(float2 *d_array,
  const size_t N, const size_t rows, const size_t cols, const size_t ndyn);
__global__ void
grad_yy_bound(float2 *d_array,
  const size_t N, const size_t rows, const size_t cols, const size_t ndyn);
__global__ void
grad_tt_bound(float2 *d_array,
  const size_t N, const size_t rows, const size_t cols, const size_t ndyn);

void
grad (float2 *d_out, float2 *d_in,
    const size_t N, const size_t rows, const size_t cols, const size_t ndyn,
    float mu, char mode);

__global__ void
grad_xx_adj (float2 *d_out, float2 *d_in,
  const size_t N, const size_t rows, const size_t cols, const size_t ndyn);
__global__ void
grad_yy_adj (float2 *d_out, float2 *d_in,
  const size_t N, const size_t rows, const size_t cols, const size_t ndyn);
__global__ void
grad_tt_adj (float2 *d_out, float2 *d_in,
  const size_t N, const size_t rows, const size_t cols, const size_t ndyn,
  float mu);

__global__ void
grad_xx_adj_zero (float2 *d_out, float2 *d_in,
  const size_t rows, const size_t cols, const size_t ndyn);
__global__ void
grad_yy_adj_zero (float2 *d_out, float2 *d_in,
  const size_t rows, const size_t cols, const size_t ndyn);
__global__ void
grad_tt_adj_zero (float2 *d_out, float2 *d_in,
  const size_t rows, const size_t cols, const size_t ndyn);

__global__ void
grad_xx_adj_bound (float2 *d_out, float2 *d_in,
  const size_t rows, const size_t cols, const size_t ndyn);
__global__ void
grad_yy_adj_bound (float2 *d_out, float2 *d_in,
  const size_t rows, const size_t cols, const size_t ndyn);
__global__ void
grad_tt_adj_bound (float2 *d_out, float2 *d_in,
  const size_t rows, const size_t cols, const size_t ndyn);

void
grad_adj (float2 *d_out, float2 *d_in,
  const size_t N, const size_t rows, const size_t cols,
  const size_t ndyn, float mu, char mode);

__global__ void
proj_p (float2 *d_p, float2 *d_tmp, float2 *d_wbar,
   const size_t N, float sigma, float alpha);
__global__ void
proj_q (float2 *d_q, float2 *d_tmp,
  const size_t N, float sigma, float alpha);

void
update_r (float2 *d_r, float2 *d_lbar, float2 *d_sbar,
  float2 *d_imgb, float2 *d_tmp, float2 *d_mask, float sigma,
  const size_t N, const size_t rows, const size_t cols);
__global__ void
update_s (float2 *d_imgs, float2 *d_tmp, float2 *d_imgz,
  const size_t N, float tau);
__global__ void
update_w (float2 *d_w, float2 *d_tmp, float2 *d_p, const size_t N, float tau);
__global__ void
update_r (float2 *d_r, float2 *d_tmp, float2 *d_imgb,
  size_t N, float sigma);

__global__ void
arrayabs(float *d_array, float2 *d_array2, const size_t array_size);
float
compute_ser(float2 *d_array, float2 *d_img, const size_t array_size);

void
tgv_cs(float2 *d_imgl, float2 *d_imgs, float2 *h_img, float2 *h_mask,
      const size_t N, const size_t rows, const size_t cols, const size_t ndyn,
      float alpha, float beta, float mu,
      float tau, float sigma, float reduction, int iter);

#ifdef __cplusplus
#endif

#ifndef MACROS_H_
#define MACROS_H_

#define cuTry(ans) gpuAssert((ans), __FILE__, __LINE__);
#define safe_malloc(n) mallocAssert(n, __FILE__, __LINE__);
#define safe_free(x) freeAssert(x, __FILE__, __LINE__);

#endif

#endif
