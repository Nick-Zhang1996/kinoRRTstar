/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * xnrm2.c
 *
 * Code generation for function 'xnrm2'
 *
 */

/* Include files */
#include <math.h>
#include "rt_nonfinite.h"
#include "benchmarkRRT.h"
#include "xnrm2.h"
#include "cost_eval.h"
#include "benchmarkRRT_rtwutil.h"

/* Function Definitions */
double xnrm2(int n, const creal_T x_data[], int ix0)
{
  double y;
  double scale;
  int kend;
  int k;
  double absxk;
  double t;
  y = 0.0;
  if (n >= 1) {
    if (n == 1) {
      y = rt_hypotd_snf(x_data[ix0 - 1].re, x_data[ix0 - 1].im);
    } else {
      scale = 3.3121686421112381E-170;
      kend = ix0 + 1;
      for (k = ix0; k <= kend; k++) {
        absxk = fabs(x_data[k - 1].re);
        if (absxk > scale) {
          t = scale / absxk;
          y = 1.0 + y * t * t;
          scale = absxk;
        } else {
          t = absxk / scale;
          y += t * t;
        }

        absxk = fabs(x_data[k - 1].im);
        if (absxk > scale) {
          t = scale / absxk;
          y = 1.0 + y * t * t;
          scale = absxk;
        } else {
          t = absxk / scale;
          y += t * t;
        }
      }

      y = scale * sqrt(y);
    }
  }

  return y;
}

/* End of code generation (xnrm2.c) */
