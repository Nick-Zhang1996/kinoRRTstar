/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * linspace.c
 *
 * Code generation for function 'linspace'
 *
 */

/* Include files */
#include <math.h>
#include "rt_nonfinite.h"
#include "benchmarkRRT.h"
#include "linspace.h"

/* Function Definitions */
void linspace(const creal_T d2, creal_T y[12])
{
  double delta1_re;
  double delta1_im;
  int k;
  y[11] = d2;
  y[0].re = 0.0;
  y[0].im = 0.0;
  if (((d2.re < 0.0) && (fabs(d2.re) > 8.9884656743115785E+307)) || ((d2.im <
        0.0) && (fabs(d2.im) > 8.9884656743115785E+307))) {
    if (d2.im == 0.0) {
      delta1_re = d2.re / 11.0;
      delta1_im = 0.0;
    } else if (d2.re == 0.0) {
      delta1_re = 0.0;
      delta1_im = d2.im / 11.0;
    } else {
      delta1_re = d2.re / 11.0;
      delta1_im = d2.im / 11.0;
    }

    for (k = 0; k < 10; k++) {
      y[k + 1].re = delta1_re * (1.0 + (double)k);
      y[k + 1].im = delta1_im * (1.0 + (double)k);
    }
  } else {
    if (d2.im == 0.0) {
      delta1_re = d2.re / 11.0;
      delta1_im = 0.0;
    } else if (d2.re == 0.0) {
      delta1_re = 0.0;
      delta1_im = d2.im / 11.0;
    } else {
      delta1_re = d2.re / 11.0;
      delta1_im = d2.im / 11.0;
    }

    for (k = 0; k < 10; k++) {
      y[k + 1].re = (1.0 + (double)k) * delta1_re;
      y[k + 1].im = (1.0 + (double)k) * delta1_im;
    }
  }
}

/* End of code generation (linspace.c) */
