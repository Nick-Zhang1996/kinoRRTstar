/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * state_traj.c
 *
 * Code generation for function 'state_traj'
 *
 */

/* Include files */
#include <math.h>
#include "rt_nonfinite.h"
#include "benchmarkRRT.h"
#include "state_traj.h"
#include "cost_eval.h"
#include "benchmarkRRT_rtwutil.h"

/* Function Definitions */
void state_traj(const creal_T t, const creal_T t_s, double x01, double x02,
                double x03, double x04, double x11, double x12, double x13,
                double x14, creal_T states[4])
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
  double brm;
  double s;
  double d;
  double c_re_tmp;
  double c_im_tmp;
  double re;
  double t_s_re_tmp;
  double t_s_im_tmp;
  double d_re_tmp;
  double d_im_tmp;
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
      y_im = 0.0;
    } else if (r == 0.0) {
      y_re = 0.0;
      y_im /= b_y_re;
    } else {
      y_re = r / b_y_re;
      y_im /= b_y_re;
    }
  } else if (b_y_re == 0.0) {
    if (r == 0.0) {
      y_re = y_im / b_y_im;
      y_im = 0.0;
    } else if (y_im == 0.0) {
      y_re = 0.0;
      y_im = -(r / b_y_im);
    } else {
      y_re = y_im / b_y_im;
      y_im = -(r / b_y_im);
    }
  } else {
    brm = fabs(b_y_re);
    a_re = fabs(b_y_im);
    if (brm > a_re) {
      s = b_y_im / b_y_re;
      d = b_y_re + s * b_y_im;
      y_re = (r + s * y_im) / d;
      y_im = (y_im - s * r) / d;
    } else if (a_re == brm) {
      if (b_y_re > 0.0) {
        s = 0.5;
      } else {
        s = -0.5;
      }

      if (b_y_im > 0.0) {
        d = 0.5;
      } else {
        d = -0.5;
      }

      y_re = (r * s + y_im * d) / brm;
      y_im = (y_im * s - r * d) / brm;
    } else {
      s = b_y_re / b_y_im;
      d = b_y_im + s * b_y_re;
      y_re = (s * r + y_im) / d;
      y_im = (s * y_im - r) / d;
    }
  }

  a_re = a_re_tmp * a_re_tmp - a_im_tmp * a_im_tmp;
  a_im = a_re_tmp * a_im_tmp + a_im_tmp * a_re_tmp;
  c_re_tmp = 2.0 * t_s.re;
  c_im_tmp = 2.0 * t_s.im;
  re = ((3.0 * x01 - 3.0 * x11) + re_tmp) + c_re_tmp * x13;
  im_tmp += c_im_tmp * x13;
  r = a_re * re - a_im * im_tmp;
  a_im = a_re * im_tmp + a_im * re;
  t_s_re_tmp = t_s.re * t_s.re - t_s.im * t_s.im;
  t_s_im_tmp = t_s.re * t_s.im + t_s.im * t_s.re;
  if (t_s_im_tmp == 0.0) {
    if (a_im == 0.0) {
      a_re = r / t_s_re_tmp;
      a_im = 0.0;
    } else if (r == 0.0) {
      a_re = 0.0;
      a_im /= t_s_re_tmp;
    } else {
      a_re = r / t_s_re_tmp;
      a_im /= t_s_re_tmp;
    }
  } else if (t_s_re_tmp == 0.0) {
    if (r == 0.0) {
      a_re = a_im / t_s_im_tmp;
      a_im = 0.0;
    } else if (a_im == 0.0) {
      a_re = 0.0;
      a_im = -(r / t_s_im_tmp);
    } else {
      a_re = a_im / t_s_im_tmp;
      a_im = -(r / t_s_im_tmp);
    }
  } else {
    brm = fabs(t_s_re_tmp);
    a_re = fabs(t_s_im_tmp);
    if (brm > a_re) {
      s = t_s_im_tmp / t_s_re_tmp;
      d = t_s_re_tmp + s * t_s_im_tmp;
      a_re = (r + s * a_im) / d;
      a_im = (a_im - s * r) / d;
    } else if (a_re == brm) {
      if (t_s_re_tmp > 0.0) {
        s = 0.5;
      } else {
        s = -0.5;
      }

      if (t_s_im_tmp > 0.0) {
        d = 0.5;
      } else {
        d = -0.5;
      }

      a_re = (r * s + a_im * d) / brm;
      a_im = (a_im * s - r * d) / brm;
    } else {
      s = t_s_re_tmp / t_s_im_tmp;
      d = t_s_im_tmp + s * t_s_re_tmp;
      a_re = (s * r + a_im) / d;
      a_im = (s * a_im - r) / d;
    }
  }

  states[0].re = ((x11 + x13 * a_re_tmp) + y_re) + a_re;
  states[0].im = (x13 * a_im_tmp + y_im) + a_im;
  re_tmp = t_s.re * x04;
  d_re_tmp = ((2.0 * x02 - 2.0 * x12) + re_tmp) + t_s.re * x14;
  im_tmp = t_s.im * x04;
  d_im_tmp = im_tmp + t_s.im * x14;
  y_re = c_y_re * d_re_tmp - c_y_im * d_im_tmp;
  y_im = c_y_re * d_im_tmp + c_y_im * d_re_tmp;
  if (d_y_im == 0.0) {
    if (y_im == 0.0) {
      b_y_re = y_re / d_y_re;
      y_im = 0.0;
    } else if (y_re == 0.0) {
      b_y_re = 0.0;
      y_im /= d_y_re;
    } else {
      b_y_re = y_re / d_y_re;
      y_im /= d_y_re;
    }
  } else if (d_y_re == 0.0) {
    if (y_re == 0.0) {
      b_y_re = y_im / d_y_im;
      y_im = 0.0;
    } else if (y_im == 0.0) {
      b_y_re = 0.0;
      y_im = -(y_re / d_y_im);
    } else {
      b_y_re = y_im / d_y_im;
      y_im = -(y_re / d_y_im);
    }
  } else {
    brm = fabs(d_y_re);
    a_re = fabs(d_y_im);
    if (brm > a_re) {
      s = d_y_im / d_y_re;
      d = d_y_re + s * d_y_im;
      b_y_re = (y_re + s * y_im) / d;
      y_im = (y_im - s * y_re) / d;
    } else if (a_re == brm) {
      if (d_y_re > 0.0) {
        s = 0.5;
      } else {
        s = -0.5;
      }

      if (d_y_im > 0.0) {
        d = 0.5;
      } else {
        d = -0.5;
      }

      b_y_re = (y_re * s + y_im * d) / brm;
      y_im = (y_im * s - y_re * d) / brm;
    } else {
      s = d_y_re / d_y_im;
      d = d_y_im + s * d_y_re;
      b_y_re = (s * y_re + y_im) / d;
      y_im = (s * y_im - y_re) / d;
    }
  }

  a_re = a_re_tmp * a_re_tmp - a_im_tmp * a_im_tmp;
  a_im = a_re_tmp * a_im_tmp + a_im_tmp * a_re_tmp;
  re = ((3.0 * x02 - 3.0 * x12) + re_tmp) + c_re_tmp * x14;
  im_tmp += c_im_tmp * x14;
  r = a_re * re - a_im * im_tmp;
  a_im = a_re * im_tmp + a_im * re;
  if (t_s_im_tmp == 0.0) {
    if (a_im == 0.0) {
      a_re = r / t_s_re_tmp;
      a_im = 0.0;
    } else if (r == 0.0) {
      a_re = 0.0;
      a_im /= t_s_re_tmp;
    } else {
      a_re = r / t_s_re_tmp;
      a_im /= t_s_re_tmp;
    }
  } else if (t_s_re_tmp == 0.0) {
    if (r == 0.0) {
      a_re = a_im / t_s_im_tmp;
      a_im = 0.0;
    } else if (a_im == 0.0) {
      a_re = 0.0;
      a_im = -(r / t_s_im_tmp);
    } else {
      a_re = a_im / t_s_im_tmp;
      a_im = -(r / t_s_im_tmp);
    }
  } else {
    brm = fabs(t_s_re_tmp);
    a_re = fabs(t_s_im_tmp);
    if (brm > a_re) {
      s = t_s_im_tmp / t_s_re_tmp;
      d = t_s_re_tmp + s * t_s_im_tmp;
      a_re = (r + s * a_im) / d;
      a_im = (a_im - s * r) / d;
    } else if (a_re == brm) {
      if (t_s_re_tmp > 0.0) {
        s = 0.5;
      } else {
        s = -0.5;
      }

      if (t_s_im_tmp > 0.0) {
        d = 0.5;
      } else {
        d = -0.5;
      }

      a_re = (r * s + a_im * d) / brm;
      a_im = (a_im * s - r * d) / brm;
    } else {
      s = t_s_re_tmp / t_s_im_tmp;
      d = t_s_im_tmp + s * t_s_re_tmp;
      a_re = (s * r + a_im) / d;
      a_im = (s * a_im - r) / d;
    }
  }

  states[1].re = ((x12 + x14 * a_re_tmp) + b_y_re) + a_re;
  states[1].im = (x14 * a_im_tmp + y_im) + a_im;
  re = 2.0 * a_re_tmp;
  im_tmp = 2.0 * a_im_tmp;
  r = ((3.0 * x01 - 3.0 * x11) + t_s.re * x03) + c_re_tmp * x13;
  b_y_im = t_s.im * x03 + c_im_tmp * x13;
  re_tmp = re * r - im_tmp * b_y_im;
  im_tmp = re * b_y_im + im_tmp * r;
  if (t_s_im_tmp == 0.0) {
    if (im_tmp == 0.0) {
      re = re_tmp / t_s_re_tmp;
      im_tmp = 0.0;
    } else if (re_tmp == 0.0) {
      re = 0.0;
      im_tmp /= t_s_re_tmp;
    } else {
      re = re_tmp / t_s_re_tmp;
      im_tmp /= t_s_re_tmp;
    }
  } else if (t_s_re_tmp == 0.0) {
    if (re_tmp == 0.0) {
      re = im_tmp / t_s_im_tmp;
      im_tmp = 0.0;
    } else if (im_tmp == 0.0) {
      re = 0.0;
      im_tmp = -(re_tmp / t_s_im_tmp);
    } else {
      re = im_tmp / t_s_im_tmp;
      im_tmp = -(re_tmp / t_s_im_tmp);
    }
  } else {
    brm = fabs(t_s_re_tmp);
    a_re = fabs(t_s_im_tmp);
    if (brm > a_re) {
      s = t_s_im_tmp / t_s_re_tmp;
      d = t_s_re_tmp + s * t_s_im_tmp;
      re = (re_tmp + s * im_tmp) / d;
      im_tmp = (im_tmp - s * re_tmp) / d;
    } else if (a_re == brm) {
      if (t_s_re_tmp > 0.0) {
        s = 0.5;
      } else {
        s = -0.5;
      }

      if (t_s_im_tmp > 0.0) {
        d = 0.5;
      } else {
        d = -0.5;
      }

      re = (re_tmp * s + im_tmp * d) / brm;
      im_tmp = (im_tmp * s - re_tmp * d) / brm;
    } else {
      s = t_s_re_tmp / t_s_im_tmp;
      d = t_s_im_tmp + s * t_s_re_tmp;
      re = (s * re_tmp + im_tmp) / d;
      im_tmp = (s * im_tmp - re_tmp) / d;
    }
  }

  r = 3.0 * (a_re_tmp * a_re_tmp - a_im_tmp * a_im_tmp);
  b_y_im = 3.0 * (a_re_tmp * a_im_tmp + a_im_tmp * a_re_tmp);
  re_tmp = r * b_re_tmp - b_y_im * b_im_tmp;
  b_y_im = r * b_im_tmp + b_y_im * b_re_tmp;
  if (e_y_im == 0.0) {
    if (b_y_im == 0.0) {
      r = re_tmp / e_y_re;
      b_y_im = 0.0;
    } else if (re_tmp == 0.0) {
      r = 0.0;
      b_y_im /= e_y_re;
    } else {
      r = re_tmp / e_y_re;
      b_y_im /= e_y_re;
    }
  } else if (e_y_re == 0.0) {
    if (re_tmp == 0.0) {
      r = b_y_im / e_y_im;
      b_y_im = 0.0;
    } else if (b_y_im == 0.0) {
      r = 0.0;
      b_y_im = -(re_tmp / e_y_im);
    } else {
      r = b_y_im / e_y_im;
      b_y_im = -(re_tmp / e_y_im);
    }
  } else {
    brm = fabs(e_y_re);
    a_re = fabs(e_y_im);
    if (brm > a_re) {
      s = e_y_im / e_y_re;
      d = e_y_re + s * e_y_im;
      r = (re_tmp + s * b_y_im) / d;
      b_y_im = (b_y_im - s * re_tmp) / d;
    } else if (a_re == brm) {
      if (e_y_re > 0.0) {
        s = 0.5;
      } else {
        s = -0.5;
      }

      if (e_y_im > 0.0) {
        d = 0.5;
      } else {
        d = -0.5;
      }

      r = (re_tmp * s + b_y_im * d) / brm;
      b_y_im = (b_y_im * s - re_tmp * d) / brm;
    } else {
      s = e_y_re / e_y_im;
      d = e_y_im + s * e_y_re;
      r = (s * re_tmp + b_y_im) / d;
      b_y_im = (s * b_y_im - re_tmp) / d;
    }
  }

  states[2].re = (x13 + re) + r;
  states[2].im = im_tmp + b_y_im;
  re = 2.0 * a_re_tmp;
  im_tmp = 2.0 * a_im_tmp;
  r = ((3.0 * x02 - 3.0 * x12) + t_s.re * x04) + c_re_tmp * x14;
  b_y_im = t_s.im * x04 + c_im_tmp * x14;
  re_tmp = re * r - im_tmp * b_y_im;
  im_tmp = re * b_y_im + im_tmp * r;
  if (t_s_im_tmp == 0.0) {
    if (im_tmp == 0.0) {
      re = re_tmp / t_s_re_tmp;
      im_tmp = 0.0;
    } else if (re_tmp == 0.0) {
      re = 0.0;
      im_tmp /= t_s_re_tmp;
    } else {
      re = re_tmp / t_s_re_tmp;
      im_tmp /= t_s_re_tmp;
    }
  } else if (t_s_re_tmp == 0.0) {
    if (re_tmp == 0.0) {
      re = im_tmp / t_s_im_tmp;
      im_tmp = 0.0;
    } else if (im_tmp == 0.0) {
      re = 0.0;
      im_tmp = -(re_tmp / t_s_im_tmp);
    } else {
      re = im_tmp / t_s_im_tmp;
      im_tmp = -(re_tmp / t_s_im_tmp);
    }
  } else {
    brm = fabs(t_s_re_tmp);
    a_re = fabs(t_s_im_tmp);
    if (brm > a_re) {
      s = t_s_im_tmp / t_s_re_tmp;
      d = t_s_re_tmp + s * t_s_im_tmp;
      re = (re_tmp + s * im_tmp) / d;
      im_tmp = (im_tmp - s * re_tmp) / d;
    } else if (a_re == brm) {
      if (t_s_re_tmp > 0.0) {
        s = 0.5;
      } else {
        s = -0.5;
      }

      if (t_s_im_tmp > 0.0) {
        d = 0.5;
      } else {
        d = -0.5;
      }

      re = (re_tmp * s + im_tmp * d) / brm;
      im_tmp = (im_tmp * s - re_tmp * d) / brm;
    } else {
      s = t_s_re_tmp / t_s_im_tmp;
      d = t_s_im_tmp + s * t_s_re_tmp;
      re = (s * re_tmp + im_tmp) / d;
      im_tmp = (s * im_tmp - re_tmp) / d;
    }
  }

  r = 3.0 * (a_re_tmp * a_re_tmp - a_im_tmp * a_im_tmp);
  b_y_im = 3.0 * (a_re_tmp * a_im_tmp + a_im_tmp * a_re_tmp);
  re_tmp = r * d_re_tmp - b_y_im * d_im_tmp;
  b_y_im = r * d_im_tmp + b_y_im * d_re_tmp;
  if (x_im == 0.0) {
    if (b_y_im == 0.0) {
      r = re_tmp / x_re;
      b_y_im = 0.0;
    } else if (re_tmp == 0.0) {
      r = 0.0;
      b_y_im /= x_re;
    } else {
      r = re_tmp / x_re;
      b_y_im /= x_re;
    }
  } else if (x_re == 0.0) {
    if (re_tmp == 0.0) {
      r = b_y_im / x_im;
      b_y_im = 0.0;
    } else if (b_y_im == 0.0) {
      r = 0.0;
      b_y_im = -(re_tmp / x_im);
    } else {
      r = b_y_im / x_im;
      b_y_im = -(re_tmp / x_im);
    }
  } else {
    brm = fabs(x_re);
    a_re = fabs(x_im);
    if (brm > a_re) {
      s = x_im / x_re;
      d = x_re + s * x_im;
      r = (re_tmp + s * b_y_im) / d;
      b_y_im = (b_y_im - s * re_tmp) / d;
    } else if (a_re == brm) {
      if (x_re > 0.0) {
        s = 0.5;
      } else {
        s = -0.5;
      }

      if (x_im > 0.0) {
        d = 0.5;
      } else {
        d = -0.5;
      }

      r = (re_tmp * s + b_y_im * d) / brm;
      b_y_im = (b_y_im * s - re_tmp * d) / brm;
    } else {
      s = x_re / x_im;
      d = x_im + s * x_re;
      r = (s * re_tmp + b_y_im) / d;
      b_y_im = (s * b_y_im - re_tmp) / d;
    }
  }

  states[3].re = (x14 + re) + r;
  states[3].im = im_tmp + b_y_im;
}

/* End of code generation (state_traj.c) */
