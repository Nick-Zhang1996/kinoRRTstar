/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * DI_state.c
 *
 * Code generation for function 'DI_state'
 *
 */

/* Include files */
#include <math.h>
#include "rt_nonfinite.h"
#include "benchmarkRRT.h"
#include "DI_state.h"
#include "DI_cost.h"
#include "benchmarkRRT_rtwutil.h"

/* Function Definitions */
void DI_state(const creal_T t, const creal_T t_s, double x01, double x02, double
              x03, double x04, double x11, double x12, double x13, double x14,
              double states[4])
{
  double a_re_tmp;
  double a_im_tmp;
  double y_re;
  double y_im;
  double a_re;
  double a_im;
  double b_y_re;
  double b_y_im;
  double x_re;
  double x_im;
  double c_y_re;
  double c_y_im;
  double r;
  double d_y_re;
  double d_y_im;
  double e_y_re;
  double e_y_im;
  double re_tmp;
  double b_re_tmp;
  double im_tmp;
  double b_im_tmp;
  double bim;
  double f_y_re;
  double c_re_tmp;
  double c_im_tmp;
  double f_y_im;
  double t_s_re_tmp;
  double t_s_im_tmp;
  double b_t_s_re_tmp;
  double d_re_tmp;
  double d_im_tmp;
  double b_t_s_im_tmp;
  double g_y_re;
  double g_y_im;
  double c_t_s_re_tmp;
  double c_t_s_im_tmp;
  double d_t_s_re_tmp;
  double d_t_s_im_tmp;
  double h_y_re;
  double h_y_im;
  double e_t_s_re_tmp;
  double e_t_s_im_tmp;
  double b_x_re;
  double b_x_im;
  a_re_tmp = t.re - t_s.re;
  a_im_tmp = t.im - t_s.im;
  if ((a_im_tmp == 0.0) && (a_re_tmp >= 0.0)) {
    y_re = rt_powd_snf(a_re_tmp, 3.0);
    y_im = 0.0;
  } else if (a_re_tmp == 0.0) {
    y_re = 0.0;
    y_im = -rt_powd_snf(a_im_tmp, 3.0);
  } else {
    if (a_im_tmp == 0.0) {
      if (a_re_tmp < 0.0) {
        a_re = log(fabs(a_re_tmp));
        a_im = 3.1415926535897931;
      } else {
        a_re = log(a_re_tmp);
        a_im = 0.0;
      }
    } else if ((fabs(a_re_tmp) > 8.9884656743115785E+307) || (fabs(a_im_tmp) >
                8.9884656743115785E+307)) {
      a_re = log(rt_hypotd_snf(a_re_tmp / 2.0, a_im_tmp / 2.0)) +
        0.69314718055994529;
      a_im = rt_atan2d_snf(a_im_tmp, a_re_tmp);
    } else {
      a_re = log(rt_hypotd_snf(a_re_tmp, a_im_tmp));
      a_im = rt_atan2d_snf(a_im_tmp, a_re_tmp);
    }

    y_re = 3.0 * a_re;
    y_im = 3.0 * a_im;
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

  if ((t_s.im == 0.0) && (t_s.re >= 0.0)) {
    b_y_re = rt_powd_snf(t_s.re, 3.0);
    b_y_im = 0.0;
  } else if (t_s.re == 0.0) {
    b_y_re = 0.0;
    b_y_im = -rt_powd_snf(t_s.im, 3.0);
  } else {
    if (t_s.im == 0.0) {
      if (t_s.re < 0.0) {
        x_re = log(fabs(t_s.re));
        x_im = 3.1415926535897931;
      } else {
        x_re = log(t_s.re);
        x_im = 0.0;
      }
    } else if ((fabs(t_s.re) > 8.9884656743115785E+307) || (fabs(t_s.im) >
                8.9884656743115785E+307)) {
      x_re = log(rt_hypotd_snf(t_s.re / 2.0, t_s.im / 2.0)) +
        0.69314718055994529;
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    } else {
      x_re = log(rt_hypotd_snf(t_s.re, t_s.im));
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    }

    b_y_re = 3.0 * x_re;
    b_y_im = 3.0 * x_im;
    if (b_y_im == 0.0) {
      b_y_re = exp(b_y_re);
      b_y_im = 0.0;
    } else if (rtIsInf(b_y_im) && rtIsInf(b_y_re) && (b_y_re < 0.0)) {
      b_y_re = 0.0;
      b_y_im = 0.0;
    } else {
      r = exp(b_y_re / 2.0);
      b_y_re = r * (r * cos(b_y_im));
      b_y_im = r * (r * sin(b_y_im));
    }
  }

  if ((a_im_tmp == 0.0) && (a_re_tmp >= 0.0)) {
    c_y_re = rt_powd_snf(a_re_tmp, 3.0);
    c_y_im = 0.0;
  } else if (a_re_tmp == 0.0) {
    c_y_re = 0.0;
    c_y_im = -rt_powd_snf(a_im_tmp, 3.0);
  } else {
    if (a_im_tmp == 0.0) {
      if (a_re_tmp < 0.0) {
        a_re = log(fabs(a_re_tmp));
        a_im = 3.1415926535897931;
      } else {
        a_re = log(a_re_tmp);
        a_im = 0.0;
      }
    } else if ((fabs(a_re_tmp) > 8.9884656743115785E+307) || (fabs(a_im_tmp) >
                8.9884656743115785E+307)) {
      a_re = log(rt_hypotd_snf(a_re_tmp / 2.0, a_im_tmp / 2.0)) +
        0.69314718055994529;
      a_im = rt_atan2d_snf(a_im_tmp, a_re_tmp);
    } else {
      a_re = log(rt_hypotd_snf(a_re_tmp, a_im_tmp));
      a_im = rt_atan2d_snf(a_im_tmp, a_re_tmp);
    }

    c_y_re = 3.0 * a_re;
    c_y_im = 3.0 * a_im;
    if (c_y_im == 0.0) {
      c_y_re = exp(c_y_re);
      c_y_im = 0.0;
    } else if (rtIsInf(c_y_im) && rtIsInf(c_y_re) && (c_y_re < 0.0)) {
      c_y_re = 0.0;
      c_y_im = 0.0;
    } else {
      r = exp(c_y_re / 2.0);
      c_y_re = r * (r * cos(c_y_im));
      c_y_im = r * (r * sin(c_y_im));
    }
  }

  if ((t_s.im == 0.0) && (t_s.re >= 0.0)) {
    d_y_re = rt_powd_snf(t_s.re, 3.0);
    d_y_im = 0.0;
  } else if (t_s.re == 0.0) {
    d_y_re = 0.0;
    d_y_im = -rt_powd_snf(t_s.im, 3.0);
  } else {
    if (t_s.im == 0.0) {
      if (t_s.re < 0.0) {
        x_re = log(fabs(t_s.re));
        x_im = 3.1415926535897931;
      } else {
        x_re = log(t_s.re);
        x_im = 0.0;
      }
    } else if ((fabs(t_s.re) > 8.9884656743115785E+307) || (fabs(t_s.im) >
                8.9884656743115785E+307)) {
      x_re = log(rt_hypotd_snf(t_s.re / 2.0, t_s.im / 2.0)) +
        0.69314718055994529;
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    } else {
      x_re = log(rt_hypotd_snf(t_s.re, t_s.im));
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    }

    d_y_re = 3.0 * x_re;
    d_y_im = 3.0 * x_im;
    if (d_y_im == 0.0) {
      d_y_re = exp(d_y_re);
      d_y_im = 0.0;
    } else if (rtIsInf(d_y_im) && rtIsInf(d_y_re) && (d_y_re < 0.0)) {
      d_y_re = 0.0;
      d_y_im = 0.0;
    } else {
      r = exp(d_y_re / 2.0);
      d_y_re = r * (r * cos(d_y_im));
      d_y_im = r * (r * sin(d_y_im));
    }
  }

  if ((t_s.im == 0.0) && (t_s.re >= 0.0)) {
    e_y_re = rt_powd_snf(t_s.re, 3.0);
    e_y_im = 0.0;
  } else if (t_s.re == 0.0) {
    e_y_re = 0.0;
    e_y_im = -rt_powd_snf(t_s.im, 3.0);
  } else {
    if (t_s.im == 0.0) {
      if (t_s.re < 0.0) {
        x_re = log(fabs(t_s.re));
        x_im = 3.1415926535897931;
      } else {
        x_re = log(t_s.re);
        x_im = 0.0;
      }
    } else if ((fabs(t_s.re) > 8.9884656743115785E+307) || (fabs(t_s.im) >
                8.9884656743115785E+307)) {
      x_re = log(rt_hypotd_snf(t_s.re / 2.0, t_s.im / 2.0)) +
        0.69314718055994529;
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    } else {
      x_re = log(rt_hypotd_snf(t_s.re, t_s.im));
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    }

    e_y_re = 3.0 * x_re;
    e_y_im = 3.0 * x_im;
    if (e_y_im == 0.0) {
      e_y_re = exp(e_y_re);
      e_y_im = 0.0;
    } else if (rtIsInf(e_y_im) && rtIsInf(e_y_re) && (e_y_re < 0.0)) {
      e_y_re = 0.0;
      e_y_im = 0.0;
    } else {
      r = exp(e_y_re / 2.0);
      e_y_re = r * (r * cos(e_y_im));
      e_y_im = r * (r * sin(e_y_im));
    }
  }

  if ((t_s.im == 0.0) && (t_s.re >= 0.0)) {
    x_re = rt_powd_snf(t_s.re, 3.0);
    x_im = 0.0;
  } else if (t_s.re == 0.0) {
    x_re = 0.0;
    x_im = -rt_powd_snf(t_s.im, 3.0);
  } else {
    if (t_s.im == 0.0) {
      if (t_s.re < 0.0) {
        x_re = log(fabs(t_s.re));
        x_im = 3.1415926535897931;
      } else {
        x_re = log(t_s.re);
        x_im = 0.0;
      }
    } else if ((fabs(t_s.re) > 8.9884656743115785E+307) || (fabs(t_s.im) >
                8.9884656743115785E+307)) {
      x_re = log(rt_hypotd_snf(t_s.re / 2.0, t_s.im / 2.0)) +
        0.69314718055994529;
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    } else {
      x_re = log(rt_hypotd_snf(t_s.re, t_s.im));
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    }

    x_re *= 3.0;
    x_im *= 3.0;
    if (x_im == 0.0) {
      x_re = exp(x_re);
      x_im = 0.0;
    } else if (rtIsInf(x_im) && rtIsInf(x_re) && (x_re < 0.0)) {
      x_re = 0.0;
      x_im = 0.0;
    } else {
      r = exp(x_re / 2.0);
      x_re = r * (r * cos(x_im));
      x_im = r * (r * sin(x_im));
    }
  }

  re_tmp = t_s.re * x03;
  b_re_tmp = ((2.0 * x01 - 2.0 * x11) + re_tmp) + t_s.re * x13;
  im_tmp = t_s.im * x03;
  b_im_tmp = im_tmp + t_s.im * x13;
  r = y_re * b_re_tmp - y_im * b_im_tmp;
  y_im = y_re * b_im_tmp + y_im * b_re_tmp;
  if (b_y_im == 0.0) {
    if (y_im == 0.0) {
      y_re = r / b_y_re;
    } else if (r == 0.0) {
      y_re = 0.0;
    } else {
      y_re = r / b_y_re;
    }
  } else if (b_y_re == 0.0) {
    if (r == 0.0) {
      y_re = y_im / b_y_im;
    } else if (y_im == 0.0) {
      y_re = 0.0;
    } else {
      y_re = y_im / b_y_im;
    }
  } else {
    a_re = fabs(b_y_re);
    bim = fabs(b_y_im);
    if (a_re > bim) {
      a_re = b_y_im / b_y_re;
      y_re = (r + a_re * y_im) / (b_y_re + a_re * b_y_im);
    } else if (bim == a_re) {
      if (b_y_re > 0.0) {
        f_y_re = 0.5;
      } else {
        f_y_re = -0.5;
      }

      if (b_y_im > 0.0) {
        f_y_im = 0.5;
      } else {
        f_y_im = -0.5;
      }

      y_re = (r * f_y_re + y_im * f_y_im) / a_re;
    } else {
      a_re = b_y_re / b_y_im;
      y_re = (a_re * r + y_im) / (b_y_im + a_re * b_y_re);
    }
  }

  a_re = a_re_tmp * a_re_tmp - a_im_tmp * a_im_tmp;
  a_im = a_re_tmp * a_im_tmp + a_im_tmp * a_re_tmp;
  c_re_tmp = 2.0 * t_s.re;
  c_im_tmp = 2.0 * t_s.im;
  y_im = ((3.0 * x01 - 3.0 * x11) + re_tmp) + c_re_tmp * x13;
  re_tmp = im_tmp + c_im_tmp * x13;
  r = a_re * y_im - a_im * re_tmp;
  a_im = a_re * re_tmp + a_im * y_im;
  t_s_re_tmp = t_s.re * t_s.re - t_s.im * t_s.im;
  t_s_im_tmp = t_s.re * t_s.im + t_s.im * t_s.re;
  if (t_s_im_tmp == 0.0) {
    if (a_im == 0.0) {
      a_re = r / t_s_re_tmp;
    } else if (r == 0.0) {
      a_re = 0.0;
    } else {
      a_re = r / t_s_re_tmp;
    }
  } else if (t_s_re_tmp == 0.0) {
    if (r == 0.0) {
      a_re = a_im / t_s_im_tmp;
    } else if (a_im == 0.0) {
      a_re = 0.0;
    } else {
      a_re = a_im / t_s_im_tmp;
    }
  } else {
    a_re = fabs(t_s_re_tmp);
    bim = fabs(t_s_im_tmp);
    if (a_re > bim) {
      a_re = t_s_im_tmp / t_s_re_tmp;
      a_re = (r + a_re * a_im) / (t_s_re_tmp + a_re * t_s_im_tmp);
    } else if (bim == a_re) {
      if (t_s_re_tmp > 0.0) {
        b_t_s_re_tmp = 0.5;
      } else {
        b_t_s_re_tmp = -0.5;
      }

      if (t_s_im_tmp > 0.0) {
        b_t_s_im_tmp = 0.5;
      } else {
        b_t_s_im_tmp = -0.5;
      }

      a_re = (r * b_t_s_re_tmp + a_im * b_t_s_im_tmp) / a_re;
    } else {
      a_re = t_s_re_tmp / t_s_im_tmp;
      a_re = (a_re * r + a_im) / (t_s_im_tmp + a_re * t_s_re_tmp);
    }
  }

  states[0] = ((x11 + x13 * a_re_tmp) + y_re) + a_re;
  re_tmp = t_s.re * x04;
  d_re_tmp = ((2.0 * x02 - 2.0 * x12) + re_tmp) + t_s.re * x14;
  im_tmp = t_s.im * x04;
  d_im_tmp = im_tmp + t_s.im * x14;
  y_re = c_y_re * d_re_tmp - c_y_im * d_im_tmp;
  y_im = c_y_re * d_im_tmp + c_y_im * d_re_tmp;
  if (d_y_im == 0.0) {
    if (y_im == 0.0) {
      y_re /= d_y_re;
    } else if (y_re == 0.0) {
      y_re = 0.0;
    } else {
      y_re /= d_y_re;
    }
  } else if (d_y_re == 0.0) {
    if (y_re == 0.0) {
      y_re = y_im / d_y_im;
    } else if (y_im == 0.0) {
      y_re = 0.0;
    } else {
      y_re = y_im / d_y_im;
    }
  } else {
    a_re = fabs(d_y_re);
    bim = fabs(d_y_im);
    if (a_re > bim) {
      a_re = d_y_im / d_y_re;
      y_re = (y_re + a_re * y_im) / (d_y_re + a_re * d_y_im);
    } else if (bim == a_re) {
      if (d_y_re > 0.0) {
        g_y_re = 0.5;
      } else {
        g_y_re = -0.5;
      }

      if (d_y_im > 0.0) {
        g_y_im = 0.5;
      } else {
        g_y_im = -0.5;
      }

      y_re = (y_re * g_y_re + y_im * g_y_im) / a_re;
    } else {
      a_re = d_y_re / d_y_im;
      y_re = (a_re * y_re + y_im) / (d_y_im + a_re * d_y_re);
    }
  }

  a_re = a_re_tmp * a_re_tmp - a_im_tmp * a_im_tmp;
  a_im = a_re_tmp * a_im_tmp + a_im_tmp * a_re_tmp;
  y_im = ((3.0 * x02 - 3.0 * x12) + re_tmp) + c_re_tmp * x14;
  re_tmp = im_tmp + c_im_tmp * x14;
  r = a_re * y_im - a_im * re_tmp;
  a_im = a_re * re_tmp + a_im * y_im;
  if (t_s_im_tmp == 0.0) {
    if (a_im == 0.0) {
      a_re = r / t_s_re_tmp;
    } else if (r == 0.0) {
      a_re = 0.0;
    } else {
      a_re = r / t_s_re_tmp;
    }
  } else if (t_s_re_tmp == 0.0) {
    if (r == 0.0) {
      a_re = a_im / t_s_im_tmp;
    } else if (a_im == 0.0) {
      a_re = 0.0;
    } else {
      a_re = a_im / t_s_im_tmp;
    }
  } else {
    a_re = fabs(t_s_re_tmp);
    bim = fabs(t_s_im_tmp);
    if (a_re > bim) {
      a_re = t_s_im_tmp / t_s_re_tmp;
      a_re = (r + a_re * a_im) / (t_s_re_tmp + a_re * t_s_im_tmp);
    } else if (bim == a_re) {
      if (t_s_re_tmp > 0.0) {
        c_t_s_re_tmp = 0.5;
      } else {
        c_t_s_re_tmp = -0.5;
      }

      if (t_s_im_tmp > 0.0) {
        c_t_s_im_tmp = 0.5;
      } else {
        c_t_s_im_tmp = -0.5;
      }

      a_re = (r * c_t_s_re_tmp + a_im * c_t_s_im_tmp) / a_re;
    } else {
      a_re = t_s_re_tmp / t_s_im_tmp;
      a_re = (a_re * r + a_im) / (t_s_im_tmp + a_re * t_s_re_tmp);
    }
  }

  states[1] = ((x12 + x14 * a_re_tmp) + y_re) + a_re;
  y_im = 2.0 * a_re_tmp;
  re_tmp = 2.0 * a_im_tmp;
  r = ((3.0 * x01 - 3.0 * x11) + t_s.re * x03) + c_re_tmp * x13;
  b_y_re = t_s.im * x03 + c_im_tmp * x13;
  b_y_im = y_im * r - re_tmp * b_y_re;
  re_tmp = y_im * b_y_re + re_tmp * r;
  if (t_s_im_tmp == 0.0) {
    if (re_tmp == 0.0) {
      y_im = b_y_im / t_s_re_tmp;
    } else if (b_y_im == 0.0) {
      y_im = 0.0;
    } else {
      y_im = b_y_im / t_s_re_tmp;
    }
  } else if (t_s_re_tmp == 0.0) {
    if (b_y_im == 0.0) {
      y_im = re_tmp / t_s_im_tmp;
    } else if (re_tmp == 0.0) {
      y_im = 0.0;
    } else {
      y_im = re_tmp / t_s_im_tmp;
    }
  } else {
    a_re = fabs(t_s_re_tmp);
    bim = fabs(t_s_im_tmp);
    if (a_re > bim) {
      a_re = t_s_im_tmp / t_s_re_tmp;
      y_im = (b_y_im + a_re * re_tmp) / (t_s_re_tmp + a_re * t_s_im_tmp);
    } else if (bim == a_re) {
      if (t_s_re_tmp > 0.0) {
        d_t_s_re_tmp = 0.5;
      } else {
        d_t_s_re_tmp = -0.5;
      }

      if (t_s_im_tmp > 0.0) {
        d_t_s_im_tmp = 0.5;
      } else {
        d_t_s_im_tmp = -0.5;
      }

      y_im = (b_y_im * d_t_s_re_tmp + re_tmp * d_t_s_im_tmp) / a_re;
    } else {
      a_re = t_s_re_tmp / t_s_im_tmp;
      y_im = (a_re * b_y_im + re_tmp) / (t_s_im_tmp + a_re * t_s_re_tmp);
    }
  }

  r = 3.0 * (a_re_tmp * a_re_tmp - a_im_tmp * a_im_tmp);
  re_tmp = 3.0 * (a_re_tmp * a_im_tmp + a_im_tmp * a_re_tmp);
  b_y_im = r * b_re_tmp - re_tmp * b_im_tmp;
  re_tmp = r * b_im_tmp + re_tmp * b_re_tmp;
  if (e_y_im == 0.0) {
    if (re_tmp == 0.0) {
      r = b_y_im / e_y_re;
    } else if (b_y_im == 0.0) {
      r = 0.0;
    } else {
      r = b_y_im / e_y_re;
    }
  } else if (e_y_re == 0.0) {
    if (b_y_im == 0.0) {
      r = re_tmp / e_y_im;
    } else if (re_tmp == 0.0) {
      r = 0.0;
    } else {
      r = re_tmp / e_y_im;
    }
  } else {
    a_re = fabs(e_y_re);
    bim = fabs(e_y_im);
    if (a_re > bim) {
      a_re = e_y_im / e_y_re;
      r = (b_y_im + a_re * re_tmp) / (e_y_re + a_re * e_y_im);
    } else if (bim == a_re) {
      if (e_y_re > 0.0) {
        h_y_re = 0.5;
      } else {
        h_y_re = -0.5;
      }

      if (e_y_im > 0.0) {
        h_y_im = 0.5;
      } else {
        h_y_im = -0.5;
      }

      r = (b_y_im * h_y_re + re_tmp * h_y_im) / a_re;
    } else {
      a_re = e_y_re / e_y_im;
      r = (a_re * b_y_im + re_tmp) / (e_y_im + a_re * e_y_re);
    }
  }

  states[2] = (x13 + y_im) + r;
  y_im = 2.0 * a_re_tmp;
  re_tmp = 2.0 * a_im_tmp;
  r = ((3.0 * x02 - 3.0 * x12) + t_s.re * x04) + c_re_tmp * x14;
  b_y_re = t_s.im * x04 + c_im_tmp * x14;
  b_y_im = y_im * r - re_tmp * b_y_re;
  re_tmp = y_im * b_y_re + re_tmp * r;
  if (t_s_im_tmp == 0.0) {
    if (re_tmp == 0.0) {
      y_im = b_y_im / t_s_re_tmp;
    } else if (b_y_im == 0.0) {
      y_im = 0.0;
    } else {
      y_im = b_y_im / t_s_re_tmp;
    }
  } else if (t_s_re_tmp == 0.0) {
    if (b_y_im == 0.0) {
      y_im = re_tmp / t_s_im_tmp;
    } else if (re_tmp == 0.0) {
      y_im = 0.0;
    } else {
      y_im = re_tmp / t_s_im_tmp;
    }
  } else {
    a_re = fabs(t_s_re_tmp);
    bim = fabs(t_s_im_tmp);
    if (a_re > bim) {
      a_re = t_s_im_tmp / t_s_re_tmp;
      y_im = (b_y_im + a_re * re_tmp) / (t_s_re_tmp + a_re * t_s_im_tmp);
    } else if (bim == a_re) {
      if (t_s_re_tmp > 0.0) {
        e_t_s_re_tmp = 0.5;
      } else {
        e_t_s_re_tmp = -0.5;
      }

      if (t_s_im_tmp > 0.0) {
        e_t_s_im_tmp = 0.5;
      } else {
        e_t_s_im_tmp = -0.5;
      }

      y_im = (b_y_im * e_t_s_re_tmp + re_tmp * e_t_s_im_tmp) / a_re;
    } else {
      a_re = t_s_re_tmp / t_s_im_tmp;
      y_im = (a_re * b_y_im + re_tmp) / (t_s_im_tmp + a_re * t_s_re_tmp);
    }
  }

  r = 3.0 * (a_re_tmp * a_re_tmp - a_im_tmp * a_im_tmp);
  re_tmp = 3.0 * (a_re_tmp * a_im_tmp + a_im_tmp * a_re_tmp);
  b_y_im = r * d_re_tmp - re_tmp * d_im_tmp;
  re_tmp = r * d_im_tmp + re_tmp * d_re_tmp;
  if (x_im == 0.0) {
    if (re_tmp == 0.0) {
      r = b_y_im / x_re;
    } else if (b_y_im == 0.0) {
      r = 0.0;
    } else {
      r = b_y_im / x_re;
    }
  } else if (x_re == 0.0) {
    if (b_y_im == 0.0) {
      r = re_tmp / x_im;
    } else if (re_tmp == 0.0) {
      r = 0.0;
    } else {
      r = re_tmp / x_im;
    }
  } else {
    a_re = fabs(x_re);
    bim = fabs(x_im);
    if (a_re > bim) {
      a_re = x_im / x_re;
      r = (b_y_im + a_re * re_tmp) / (x_re + a_re * x_im);
    } else if (bim == a_re) {
      if (x_re > 0.0) {
        b_x_re = 0.5;
      } else {
        b_x_re = -0.5;
      }

      if (x_im > 0.0) {
        b_x_im = 0.5;
      } else {
        b_x_im = -0.5;
      }

      r = (b_y_im * b_x_re + re_tmp * b_x_im) / a_re;
    } else {
      a_re = x_re / x_im;
      r = (a_re * b_y_im + re_tmp) / (x_im + a_re * x_re);
    }
  }

  states[3] = (x14 + y_im) + r;
}

/* End of code generation (DI_state.c) */
