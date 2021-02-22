/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * DI_costFreeVel.c
 *
 * Code generation for function 'DI_costFreeVel'
 *
 */

/* Include files */
#include <math.h>
#include "rt_nonfinite.h"
#include "benchmarkRRT.h"
#include "DI_costFreeVel.h"
#include "DI_cost.h"
#include "benchmarkRRT_rtwutil.h"

/* Function Definitions */
double DI_costFreeVel(const creal_T t_s, double x01, double x02, double x03,
                      double x04, double x11, double x12)
{
  double y_re;
  double y_im;
  double ar;
  double re;
  double r;
  double bim;
  double b_re;
  double b_t_s;
  double c_t_s;
  double b_y_re;
  double b_y_im;
  double c_y_re;
  double c_y_im;
  if ((t_s.im == 0.0) && (t_s.re >= 0.0)) {
    y_re = rt_powd_snf(t_s.re, 3.0);
    y_im = 0.0;
  } else if (t_s.re == 0.0) {
    y_re = 0.0;
    y_im = -rt_powd_snf(t_s.im, 3.0);
  } else {
    if (t_s.im == 0.0) {
      if (t_s.re < 0.0) {
        y_re = log(fabs(t_s.re));
        y_im = 3.1415926535897931;
      } else {
        y_re = log(t_s.re);
        y_im = 0.0;
      }
    } else if ((fabs(t_s.re) > 8.9884656743115785E+307) || (fabs(t_s.im) >
                8.9884656743115785E+307)) {
      y_re = log(rt_hypotd_snf(t_s.re / 2.0, t_s.im / 2.0)) +
        0.69314718055994529;
      y_im = rt_atan2d_snf(t_s.im, t_s.re);
    } else {
      y_re = log(rt_hypotd_snf(t_s.re, t_s.im));
      y_im = rt_atan2d_snf(t_s.im, t_s.re);
    }

    y_re *= 3.0;
    y_im *= 3.0;
    if (y_im == 0.0) {
      y_re = exp(y_re);
      y_im = 0.0;
    } else if (rtIsInf(y_im) && rtIsInf(y_re) && (y_re < 0.0)) {
      y_re = 0.0;
      y_im = 0.0;
    } else {
      r = exp(y_re / 2.0);
      y_re = r * (r * cos(y_im));
      y_im = r * (r * sin(y_im));
    }
  }

  ar = 3.0 * (x03 * x03) + 3.0 * (x04 * x04);
  if (t_s.im == 0.0) {
    re = ar / t_s.re;
  } else if (t_s.re == 0.0) {
    if (ar == 0.0) {
      re = 0.0 / t_s.im;
    } else {
      re = 0.0;
    }
  } else {
    r = fabs(t_s.re);
    bim = fabs(t_s.im);
    if (r > bim) {
      r = t_s.im / t_s.re;
      re = (ar + r * 0.0) / (t_s.re + r * t_s.im);
    } else if (bim == r) {
      if (t_s.re > 0.0) {
        b_t_s = 0.5;
      } else {
        b_t_s = -0.5;
      }

      if (t_s.im > 0.0) {
        c_t_s = 0.5;
      } else {
        c_t_s = -0.5;
      }

      re = (ar * b_t_s + 0.0 * c_t_s) / r;
    } else {
      r = t_s.re / t_s.im;
      re = r * ar / (t_s.im + r * t_s.re);
    }
  }

  ar = ((((3.0 * (x01 * x01) - 6.0 * x01 * x11) + 3.0 * (x02 * x02)) - 6.0 * x02
         * x12) + 3.0 * (x11 * x11)) + 3.0 * (x12 * x12);
  if (y_im == 0.0) {
    b_re = ar / y_re;
  } else if (y_re == 0.0) {
    if (ar == 0.0) {
      b_re = 0.0 / y_im;
    } else {
      b_re = 0.0;
    }
  } else {
    r = fabs(y_re);
    bim = fabs(y_im);
    if (r > bim) {
      r = y_im / y_re;
      b_re = (ar + r * 0.0) / (y_re + r * y_im);
    } else if (bim == r) {
      if (y_re > 0.0) {
        b_y_re = 0.5;
      } else {
        b_y_re = -0.5;
      }

      if (y_im > 0.0) {
        b_y_im = 0.5;
      } else {
        b_y_im = -0.5;
      }

      b_re = (ar * b_y_re + 0.0 * b_y_im) / r;
    } else {
      r = y_re / y_im;
      b_re = r * ar / (y_im + r * y_re);
    }
  }

  y_re = t_s.re * t_s.re - t_s.im * t_s.im;
  y_im = t_s.re * t_s.im + t_s.im * t_s.re;
  ar = ((6.0 * x01 * x03 + 6.0 * x02 * x04) - 6.0 * x03 * x11) - 6.0 * x04 * x12;
  if (y_im == 0.0) {
    r = ar / y_re;
  } else if (y_re == 0.0) {
    if (ar == 0.0) {
      r = 0.0 / y_im;
    } else {
      r = 0.0;
    }
  } else {
    r = fabs(y_re);
    bim = fabs(y_im);
    if (r > bim) {
      r = y_im / y_re;
      r = (ar + r * 0.0) / (y_re + r * y_im);
    } else if (bim == r) {
      if (y_re > 0.0) {
        c_y_re = 0.5;
      } else {
        c_y_re = -0.5;
      }

      if (y_im > 0.0) {
        c_y_im = 0.5;
      } else {
        c_y_im = -0.5;
      }

      r = (ar * c_y_re + 0.0 * c_y_im) / r;
    } else {
      r = y_re / y_im;
      r = r * ar / (y_im + r * y_re);
    }
  }

  return ((t_s.re + re) + b_re) + r;
}

/* End of code generation (DI_costFreeVel.c) */
