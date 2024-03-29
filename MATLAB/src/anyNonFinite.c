/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * anyNonFinite.c
 *
 * Code generation for function 'anyNonFinite'
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "benchmarkRRT.h"
#include "anyNonFinite.h"

/* Function Definitions */
boolean_T anyNonFinite(const creal_T x_data[], const int x_size[2])
{
  boolean_T p;
  int nx;
  int k;
  nx = x_size[0] * x_size[1];
  p = true;
  for (k = 0; k < nx; k++) {
    if (p && ((!rtIsInf(x_data[k].re)) && (!rtIsInf(x_data[k].im)) && ((!rtIsNaN
           (x_data[k].re)) && (!rtIsNaN(x_data[k].im))))) {
      p = true;
    } else {
      p = false;
    }
  }

  return !p;
}

/* End of code generation (anyNonFinite.c) */
