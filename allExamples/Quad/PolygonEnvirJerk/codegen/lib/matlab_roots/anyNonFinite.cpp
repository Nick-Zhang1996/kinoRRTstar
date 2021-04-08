/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * anyNonFinite.cpp
 *
 * Code generation for function 'anyNonFinite'
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "matlab_roots.h"
#include "anyNonFinite.h"

/* Function Definitions */
boolean_T anyNonFinite(const emxArray_creal_T *x)
{
  boolean_T p;
  int nx;
  int k;
  nx = x->size[0] * x->size[1];
  p = true;
  for (k = 0; k < nx; k++) {
    if (p && ((!rtIsInf(x->data[k].re)) && (!rtIsInf(x->data[k].im)) &&
              ((!rtIsNaN(x->data[k].re)) && (!rtIsNaN(x->data[k].im))))) {
      p = true;
    } else {
      p = false;
    }
  }

  return !p;
}

/* End of code generation (anyNonFinite.cpp) */
