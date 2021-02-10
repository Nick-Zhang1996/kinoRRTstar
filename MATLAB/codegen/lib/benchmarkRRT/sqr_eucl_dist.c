/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * sqr_eucl_dist.c
 *
 * Code generation for function 'sqr_eucl_dist'
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "benchmarkRRT.h"
#include "sqr_eucl_dist.h"
#include "benchmarkRRT_emxutil.h"

/* Function Definitions */
void sqr_eucl_dist(const emxArray_real_T *array, emxArray_real_T *e_dist)
{
  emxArray_real_T *sqr_e_dist;
  int i1;
  int loop_ub;
  int vstride;
  emxInit_real_T(&sqr_e_dist, 2);
  i1 = sqr_e_dist->size[0] * sqr_e_dist->size[1];
  sqr_e_dist->size[0] = array->size[0];
  sqr_e_dist->size[1] = 2;
  emxEnsureCapacity_real_T(sqr_e_dist, i1);
  loop_ub = array->size[0] << 1;
  for (i1 = 0; i1 < loop_ub; i1++) {
    sqr_e_dist->data[i1] = 0.0;
  }

  loop_ub = array->size[0] - 1;
  for (i1 = 0; i1 <= loop_ub; i1++) {
    sqr_e_dist->data[i1] = array->data[i1] * array->data[i1];
  }

  loop_ub = array->size[0] - 1;
  for (i1 = 0; i1 <= loop_ub; i1++) {
    sqr_e_dist->data[i1 + sqr_e_dist->size[0]] = array->data[i1 + array->size[0]]
      * array->data[i1 + array->size[0]];
  }

  if (sqr_e_dist->size[0] == 0) {
    e_dist->size[0] = 0;
  } else {
    vstride = sqr_e_dist->size[0];
    i1 = e_dist->size[0];
    e_dist->size[0] = sqr_e_dist->size[0];
    emxEnsureCapacity_real_T(e_dist, i1);
    for (loop_ub = 0; loop_ub < vstride; loop_ub++) {
      e_dist->data[loop_ub] = sqr_e_dist->data[loop_ub];
    }

    for (loop_ub = 0; loop_ub < vstride; loop_ub++) {
      e_dist->data[loop_ub] += sqr_e_dist->data[vstride + loop_ub];
    }
  }

  emxFree_real_T(&sqr_e_dist);
}

/* End of code generation (sqr_eucl_dist.c) */
