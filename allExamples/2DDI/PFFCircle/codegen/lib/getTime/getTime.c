/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: getTime.c
 *
 * MATLAB Coder version            : 4.1
 * C/C++ source code generated on  : 22-Feb-2021 23:20:10
 */

/* Include Files */
#include "getTime.h"
#include "matlab_c_timing_utils.h"

/* Function Definitions */

/*
 * Arguments    : void
 * Return Type  : double
 */
double getTime(void)
{
  return matlab_c_get_current_time_in_sec();
}

/*
 * File trailer for getTime.c
 *
 * [EOF]
 */
