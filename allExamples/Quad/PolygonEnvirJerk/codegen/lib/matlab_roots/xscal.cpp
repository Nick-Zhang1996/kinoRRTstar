/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * xscal.cpp
 *
 * Code generation for function 'xscal'
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "matlab_roots.h"
#include "xscal.h"

/* Function Definitions */
void b_xscal(int n, const creal_T a, emxArray_creal_T *x, int ix0, int incx)
{
  int i4;
  int k;
  double x_re;
  double x_im;
  i4 = ix0 + incx * (n - 1);
  for (k = ix0; incx < 0 ? k >= i4 : k <= i4; k += incx) {
    x_re = x->data[k - 1].re;
    x_im = x->data[k - 1].im;
    x->data[k - 1].re = a.re * x_re - a.im * x_im;
    x->data[k - 1].im = a.re * x_im + a.im * x_re;
  }
}

void xscal(int n, const creal_T a, emxArray_creal_T *x, int ix0)
{
  int i3;
  int k;
  double x_re;
  double x_im;
  i3 = (ix0 + n) - 1;
  for (k = ix0; k <= i3; k++) {
    x_re = x->data[k - 1].re;
    x_im = x->data[k - 1].im;
    x->data[k - 1].re = a.re * x_re - a.im * x_im;
    x->data[k - 1].im = a.re * x_im + a.im * x_re;
  }
}

/* End of code generation (xscal.cpp) */
