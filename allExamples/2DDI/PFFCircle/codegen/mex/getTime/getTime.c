/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * getTime.c
 *
 * Code generation for function 'getTime'
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "getTime.h"
#include "getTime_data.h"
#include "matlab_c_timing_utils.h"

/* Function Definitions */
real_T getTime(const emlrtStack *sp)
{
  (void)sp;
  covrtLogFcn(&emlrtCoverageInstance, 0U, 0U);
  covrtLogBasicBlock(&emlrtCoverageInstance, 0U, 0U);
  return matlab_c_get_current_time_in_sec();
}

/* End of code generation (getTime.c) */
