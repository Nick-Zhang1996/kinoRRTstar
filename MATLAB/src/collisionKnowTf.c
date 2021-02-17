/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * collisionKnowTf.c
 *
 * Code generation for function 'collisionKnowTf'
 *
 */

/* Include files */
#include <math.h>
#include "rt_nonfinite.h"
#include "benchmarkRRT.h"
#include "collisionKnowTf.h"
#include "state_traj.h"
#include "linspace.h"
#include "cost_eval.h"
#include "benchmarkRRT_rtwutil.h"
#include "benchmarkRRT_data.h"

/* Function Declarations */
static void __anon_fcn(const creal_T Tf, const double parent_data[], const
  creal_T node[8], const creal_T t, creal_T varargout_1[4]);

/* Function Definitions */
static void __anon_fcn(const creal_T Tf, const double parent_data[], const
  creal_T node[8], const creal_T t, creal_T varargout_1[4])
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
  double Tf_re_tmp;
  double Tf_im_tmp;
  double re_tmp;
  double re;
  double im_tmp;
  double im;
  double brm;
  double s;
  double d;
  double b_re_tmp;
  double b_im_tmp;
  double b_Tf_re_tmp;
  double b_Tf_im_tmp;
  double c_Tf_re_tmp;
  double c_Tf_im_tmp;
  a_re_tmp = t.re - Tf.re;
  a_im_tmp = t.im - Tf.im;
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

  if ((Tf.im == 0.0) && (Tf.re >= 0.0)) {
    b_y_re = rt_powd_snf(Tf.re, 3.0);
    b_y_im = 0.0;
  } else if (Tf.re == 0.0) {
    b_y_re = 0.0;
    b_y_im = -rt_powd_snf(Tf.im, 3.0);
  } else {
    if (Tf.im == 0.0) {
      if (Tf.re < 0.0) {
        x_re = log(fabs(Tf.re));
        x_im = 3.1415926535897931;
      } else {
        x_re = log(Tf.re);
        x_im = 0.0;
      }
    } else if ((fabs(Tf.re) > 8.9884656743115785E+307) || (fabs(Tf.im) >
                8.9884656743115785E+307)) {
      x_re = log(rt_hypotd_snf(Tf.re / 2.0, Tf.im / 2.0)) + 0.69314718055994529;
      x_im = rt_atan2d_snf(Tf.im, Tf.re);
    } else {
      x_re = log(rt_hypotd_snf(Tf.re, Tf.im));
      x_im = rt_atan2d_snf(Tf.im, Tf.re);
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

  if ((Tf.im == 0.0) && (Tf.re >= 0.0)) {
    d_y_re = rt_powd_snf(Tf.re, 3.0);
    d_y_im = 0.0;
  } else if (Tf.re == 0.0) {
    d_y_re = 0.0;
    d_y_im = -rt_powd_snf(Tf.im, 3.0);
  } else {
    if (Tf.im == 0.0) {
      if (Tf.re < 0.0) {
        x_re = log(fabs(Tf.re));
        x_im = 3.1415926535897931;
      } else {
        x_re = log(Tf.re);
        x_im = 0.0;
      }
    } else if ((fabs(Tf.re) > 8.9884656743115785E+307) || (fabs(Tf.im) >
                8.9884656743115785E+307)) {
      x_re = log(rt_hypotd_snf(Tf.re / 2.0, Tf.im / 2.0)) + 0.69314718055994529;
      x_im = rt_atan2d_snf(Tf.im, Tf.re);
    } else {
      x_re = log(rt_hypotd_snf(Tf.re, Tf.im));
      x_im = rt_atan2d_snf(Tf.im, Tf.re);
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

  if ((Tf.im == 0.0) && (Tf.re >= 0.0)) {
    e_y_re = rt_powd_snf(Tf.re, 3.0);
    e_y_im = 0.0;
  } else if (Tf.re == 0.0) {
    e_y_re = 0.0;
    e_y_im = -rt_powd_snf(Tf.im, 3.0);
  } else {
    if (Tf.im == 0.0) {
      if (Tf.re < 0.0) {
        x_re = log(fabs(Tf.re));
        x_im = 3.1415926535897931;
      } else {
        x_re = log(Tf.re);
        x_im = 0.0;
      }
    } else if ((fabs(Tf.re) > 8.9884656743115785E+307) || (fabs(Tf.im) >
                8.9884656743115785E+307)) {
      x_re = log(rt_hypotd_snf(Tf.re / 2.0, Tf.im / 2.0)) + 0.69314718055994529;
      x_im = rt_atan2d_snf(Tf.im, Tf.re);
    } else {
      x_re = log(rt_hypotd_snf(Tf.re, Tf.im));
      x_im = rt_atan2d_snf(Tf.im, Tf.re);
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

  if ((Tf.im == 0.0) && (Tf.re >= 0.0)) {
    x_re = rt_powd_snf(Tf.re, 3.0);
    x_im = 0.0;
  } else if (Tf.re == 0.0) {
    x_re = 0.0;
    x_im = -rt_powd_snf(Tf.im, 3.0);
  } else {
    if (Tf.im == 0.0) {
      if (Tf.re < 0.0) {
        x_re = log(fabs(Tf.re));
        x_im = 3.1415926535897931;
      } else {
        x_re = log(Tf.re);
        x_im = 0.0;
      }
    } else if ((fabs(Tf.re) > 8.9884656743115785E+307) || (fabs(Tf.im) >
                8.9884656743115785E+307)) {
      x_re = log(rt_hypotd_snf(Tf.re / 2.0, Tf.im / 2.0)) + 0.69314718055994529;
      x_im = rt_atan2d_snf(Tf.im, Tf.re);
    } else {
      x_re = log(rt_hypotd_snf(Tf.re, Tf.im));
      x_im = rt_atan2d_snf(Tf.im, Tf.re);
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

  Tf_re_tmp = Tf.re * node[2].re - Tf.im * node[2].im;
  Tf_im_tmp = Tf.re * node[2].im + Tf.im * node[2].re;
  re_tmp = Tf.re * parent_data[2];
  re = ((2.0 * parent_data[0] - 2.0 * node[0].re) + re_tmp) + Tf_re_tmp;
  im_tmp = Tf.im * parent_data[2];
  im = ((0.0 - 2.0 * node[0].im) + im_tmp) + Tf_im_tmp;
  r = y_re * re - y_im * im;
  y_im = y_re * im + y_im * re;
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
  b_re_tmp = 2.0 * Tf.re;
  b_im_tmp = 2.0 * Tf.im;
  re = ((3.0 * parent_data[0] - 3.0 * node[0].re) + re_tmp) + (b_re_tmp * node[2]
    .re - b_im_tmp * node[2].im);
  im = ((0.0 - 3.0 * node[0].im) + im_tmp) + (b_re_tmp * node[2].im + b_im_tmp *
    node[2].re);
  r = a_re * re - a_im * im;
  a_im = a_re * im + a_im * re;
  b_Tf_re_tmp = Tf.re * Tf.re - Tf.im * Tf.im;
  b_Tf_im_tmp = Tf.re * Tf.im + Tf.im * Tf.re;
  if (b_Tf_im_tmp == 0.0) {
    if (a_im == 0.0) {
      a_re = r / b_Tf_re_tmp;
      a_im = 0.0;
    } else if (r == 0.0) {
      a_re = 0.0;
      a_im /= b_Tf_re_tmp;
    } else {
      a_re = r / b_Tf_re_tmp;
      a_im /= b_Tf_re_tmp;
    }
  } else if (b_Tf_re_tmp == 0.0) {
    if (r == 0.0) {
      a_re = a_im / b_Tf_im_tmp;
      a_im = 0.0;
    } else if (a_im == 0.0) {
      a_re = 0.0;
      a_im = -(r / b_Tf_im_tmp);
    } else {
      a_re = a_im / b_Tf_im_tmp;
      a_im = -(r / b_Tf_im_tmp);
    }
  } else {
    brm = fabs(b_Tf_re_tmp);
    a_re = fabs(b_Tf_im_tmp);
    if (brm > a_re) {
      s = b_Tf_im_tmp / b_Tf_re_tmp;
      d = b_Tf_re_tmp + s * b_Tf_im_tmp;
      a_re = (r + s * a_im) / d;
      a_im = (a_im - s * r) / d;
    } else if (a_re == brm) {
      if (b_Tf_re_tmp > 0.0) {
        s = 0.5;
      } else {
        s = -0.5;
      }

      if (b_Tf_im_tmp > 0.0) {
        d = 0.5;
      } else {
        d = -0.5;
      }

      a_re = (r * s + a_im * d) / brm;
      a_im = (a_im * s - r * d) / brm;
    } else {
      s = b_Tf_re_tmp / b_Tf_im_tmp;
      d = b_Tf_im_tmp + s * b_Tf_re_tmp;
      a_re = (s * r + a_im) / d;
      a_im = (s * a_im - r) / d;
    }
  }

  varargout_1[0].re = ((node[0].re + (node[2].re * a_re_tmp - node[2].im *
    a_im_tmp)) + y_re) + a_re;
  varargout_1[0].im = ((node[0].im + (node[2].re * a_im_tmp + node[2].im *
    a_re_tmp)) + y_im) + a_im;
  c_Tf_re_tmp = Tf.re * node[3].re - Tf.im * node[3].im;
  c_Tf_im_tmp = Tf.re * node[3].im + Tf.im * node[3].re;
  re_tmp = Tf.re * parent_data[3];
  re = ((2.0 * parent_data[1] - 2.0 * node[1].re) + re_tmp) + c_Tf_re_tmp;
  im_tmp = Tf.im * parent_data[3];
  im = ((0.0 - 2.0 * node[1].im) + im_tmp) + c_Tf_im_tmp;
  y_re = c_y_re * re - c_y_im * im;
  y_im = c_y_re * im + c_y_im * re;
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
  re = ((3.0 * parent_data[1] - 3.0 * node[1].re) + re_tmp) + (b_re_tmp * node[3]
    .re - b_im_tmp * node[3].im);
  im = ((0.0 - 3.0 * node[1].im) + im_tmp) + (b_re_tmp * node[3].im + b_im_tmp *
    node[3].re);
  r = a_re * re - a_im * im;
  a_im = a_re * im + a_im * re;
  if (b_Tf_im_tmp == 0.0) {
    if (a_im == 0.0) {
      a_re = r / b_Tf_re_tmp;
      a_im = 0.0;
    } else if (r == 0.0) {
      a_re = 0.0;
      a_im /= b_Tf_re_tmp;
    } else {
      a_re = r / b_Tf_re_tmp;
      a_im /= b_Tf_re_tmp;
    }
  } else if (b_Tf_re_tmp == 0.0) {
    if (r == 0.0) {
      a_re = a_im / b_Tf_im_tmp;
      a_im = 0.0;
    } else if (a_im == 0.0) {
      a_re = 0.0;
      a_im = -(r / b_Tf_im_tmp);
    } else {
      a_re = a_im / b_Tf_im_tmp;
      a_im = -(r / b_Tf_im_tmp);
    }
  } else {
    brm = fabs(b_Tf_re_tmp);
    a_re = fabs(b_Tf_im_tmp);
    if (brm > a_re) {
      s = b_Tf_im_tmp / b_Tf_re_tmp;
      d = b_Tf_re_tmp + s * b_Tf_im_tmp;
      a_re = (r + s * a_im) / d;
      a_im = (a_im - s * r) / d;
    } else if (a_re == brm) {
      if (b_Tf_re_tmp > 0.0) {
        s = 0.5;
      } else {
        s = -0.5;
      }

      if (b_Tf_im_tmp > 0.0) {
        d = 0.5;
      } else {
        d = -0.5;
      }

      a_re = (r * s + a_im * d) / brm;
      a_im = (a_im * s - r * d) / brm;
    } else {
      s = b_Tf_re_tmp / b_Tf_im_tmp;
      d = b_Tf_im_tmp + s * b_Tf_re_tmp;
      a_re = (s * r + a_im) / d;
      a_im = (s * a_im - r) / d;
    }
  }

  varargout_1[1].re = ((node[1].re + (node[3].re * a_re_tmp - node[3].im *
    a_im_tmp)) + b_y_re) + a_re;
  varargout_1[1].im = ((node[1].im + (node[3].re * a_im_tmp + node[3].im *
    a_re_tmp)) + y_im) + a_im;
  re = 2.0 * a_re_tmp;
  im = 2.0 * a_im_tmp;
  a_re = ((3.0 * parent_data[0] - 3.0 * node[0].re) + Tf.re * parent_data[2]) +
    (b_re_tmp * node[2].re - b_im_tmp * node[2].im);
  a_im = ((0.0 - 3.0 * node[0].im) + Tf.im * parent_data[2]) + (b_re_tmp * node
    [2].im + b_im_tmp * node[2].re);
  re_tmp = re * a_re - im * a_im;
  im = re * a_im + im * a_re;
  if (b_Tf_im_tmp == 0.0) {
    if (im == 0.0) {
      re = re_tmp / b_Tf_re_tmp;
      im = 0.0;
    } else if (re_tmp == 0.0) {
      re = 0.0;
      im /= b_Tf_re_tmp;
    } else {
      re = re_tmp / b_Tf_re_tmp;
      im /= b_Tf_re_tmp;
    }
  } else if (b_Tf_re_tmp == 0.0) {
    if (re_tmp == 0.0) {
      re = im / b_Tf_im_tmp;
      im = 0.0;
    } else if (im == 0.0) {
      re = 0.0;
      im = -(re_tmp / b_Tf_im_tmp);
    } else {
      re = im / b_Tf_im_tmp;
      im = -(re_tmp / b_Tf_im_tmp);
    }
  } else {
    brm = fabs(b_Tf_re_tmp);
    a_re = fabs(b_Tf_im_tmp);
    if (brm > a_re) {
      s = b_Tf_im_tmp / b_Tf_re_tmp;
      d = b_Tf_re_tmp + s * b_Tf_im_tmp;
      re = (re_tmp + s * im) / d;
      im = (im - s * re_tmp) / d;
    } else if (a_re == brm) {
      if (b_Tf_re_tmp > 0.0) {
        s = 0.5;
      } else {
        s = -0.5;
      }

      if (b_Tf_im_tmp > 0.0) {
        d = 0.5;
      } else {
        d = -0.5;
      }

      re = (re_tmp * s + im * d) / brm;
      im = (im * s - re_tmp * d) / brm;
    } else {
      s = b_Tf_re_tmp / b_Tf_im_tmp;
      d = b_Tf_im_tmp + s * b_Tf_re_tmp;
      re = (s * re_tmp + im) / d;
      im = (s * im - re_tmp) / d;
    }
  }

  a_re = 3.0 * (a_re_tmp * a_re_tmp - a_im_tmp * a_im_tmp);
  a_im = 3.0 * (a_re_tmp * a_im_tmp + a_im_tmp * a_re_tmp);
  re_tmp = ((2.0 * parent_data[0] - 2.0 * node[0].re) + Tf.re * parent_data[2])
    + Tf_re_tmp;
  r = ((0.0 - 2.0 * node[0].im) + Tf.im * parent_data[2]) + Tf_im_tmp;
  b_y_im = a_re * re_tmp - a_im * r;
  a_im = a_re * r + a_im * re_tmp;
  if (e_y_im == 0.0) {
    if (a_im == 0.0) {
      a_re = b_y_im / e_y_re;
      a_im = 0.0;
    } else if (b_y_im == 0.0) {
      a_re = 0.0;
      a_im /= e_y_re;
    } else {
      a_re = b_y_im / e_y_re;
      a_im /= e_y_re;
    }
  } else if (e_y_re == 0.0) {
    if (b_y_im == 0.0) {
      a_re = a_im / e_y_im;
      a_im = 0.0;
    } else if (a_im == 0.0) {
      a_re = 0.0;
      a_im = -(b_y_im / e_y_im);
    } else {
      a_re = a_im / e_y_im;
      a_im = -(b_y_im / e_y_im);
    }
  } else {
    brm = fabs(e_y_re);
    a_re = fabs(e_y_im);
    if (brm > a_re) {
      s = e_y_im / e_y_re;
      d = e_y_re + s * e_y_im;
      a_re = (b_y_im + s * a_im) / d;
      a_im = (a_im - s * b_y_im) / d;
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

      a_re = (b_y_im * s + a_im * d) / brm;
      a_im = (a_im * s - b_y_im * d) / brm;
    } else {
      s = e_y_re / e_y_im;
      d = e_y_im + s * e_y_re;
      a_re = (s * b_y_im + a_im) / d;
      a_im = (s * a_im - b_y_im) / d;
    }
  }

  varargout_1[2].re = (node[2].re + re) + a_re;
  varargout_1[2].im = (node[2].im + im) + a_im;
  re = 2.0 * a_re_tmp;
  im = 2.0 * a_im_tmp;
  a_re = ((3.0 * parent_data[1] - 3.0 * node[1].re) + Tf.re * parent_data[3]) +
    (b_re_tmp * node[3].re - b_im_tmp * node[3].im);
  a_im = ((0.0 - 3.0 * node[1].im) + Tf.im * parent_data[3]) + (b_re_tmp * node
    [3].im + b_im_tmp * node[3].re);
  re_tmp = re * a_re - im * a_im;
  im = re * a_im + im * a_re;
  if (b_Tf_im_tmp == 0.0) {
    if (im == 0.0) {
      re = re_tmp / b_Tf_re_tmp;
      im = 0.0;
    } else if (re_tmp == 0.0) {
      re = 0.0;
      im /= b_Tf_re_tmp;
    } else {
      re = re_tmp / b_Tf_re_tmp;
      im /= b_Tf_re_tmp;
    }
  } else if (b_Tf_re_tmp == 0.0) {
    if (re_tmp == 0.0) {
      re = im / b_Tf_im_tmp;
      im = 0.0;
    } else if (im == 0.0) {
      re = 0.0;
      im = -(re_tmp / b_Tf_im_tmp);
    } else {
      re = im / b_Tf_im_tmp;
      im = -(re_tmp / b_Tf_im_tmp);
    }
  } else {
    brm = fabs(b_Tf_re_tmp);
    a_re = fabs(b_Tf_im_tmp);
    if (brm > a_re) {
      s = b_Tf_im_tmp / b_Tf_re_tmp;
      d = b_Tf_re_tmp + s * b_Tf_im_tmp;
      re = (re_tmp + s * im) / d;
      im = (im - s * re_tmp) / d;
    } else if (a_re == brm) {
      if (b_Tf_re_tmp > 0.0) {
        s = 0.5;
      } else {
        s = -0.5;
      }

      if (b_Tf_im_tmp > 0.0) {
        d = 0.5;
      } else {
        d = -0.5;
      }

      re = (re_tmp * s + im * d) / brm;
      im = (im * s - re_tmp * d) / brm;
    } else {
      s = b_Tf_re_tmp / b_Tf_im_tmp;
      d = b_Tf_im_tmp + s * b_Tf_re_tmp;
      re = (s * re_tmp + im) / d;
      im = (s * im - re_tmp) / d;
    }
  }

  a_re = 3.0 * (a_re_tmp * a_re_tmp - a_im_tmp * a_im_tmp);
  a_im = 3.0 * (a_re_tmp * a_im_tmp + a_im_tmp * a_re_tmp);
  re_tmp = ((2.0 * parent_data[1] - 2.0 * node[1].re) + Tf.re * parent_data[3])
    + c_Tf_re_tmp;
  r = ((0.0 - 2.0 * node[1].im) + Tf.im * parent_data[3]) + c_Tf_im_tmp;
  b_y_im = a_re * re_tmp - a_im * r;
  a_im = a_re * r + a_im * re_tmp;
  if (x_im == 0.0) {
    if (a_im == 0.0) {
      a_re = b_y_im / x_re;
      a_im = 0.0;
    } else if (b_y_im == 0.0) {
      a_re = 0.0;
      a_im /= x_re;
    } else {
      a_re = b_y_im / x_re;
      a_im /= x_re;
    }
  } else if (x_re == 0.0) {
    if (b_y_im == 0.0) {
      a_re = a_im / x_im;
      a_im = 0.0;
    } else if (a_im == 0.0) {
      a_re = 0.0;
      a_im = -(b_y_im / x_im);
    } else {
      a_re = a_im / x_im;
      a_im = -(b_y_im / x_im);
    }
  } else {
    brm = fabs(x_re);
    a_re = fabs(x_im);
    if (brm > a_re) {
      s = x_im / x_re;
      d = x_re + s * x_im;
      a_re = (b_y_im + s * a_im) / d;
      a_im = (a_im - s * b_y_im) / d;
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

      a_re = (b_y_im * s + a_im * d) / brm;
      a_im = (a_im * s - b_y_im * d) / brm;
    } else {
      s = x_re / x_im;
      d = x_im + s * x_re;
      a_re = (s * b_y_im + a_im) / d;
      a_im = (s * a_im - b_y_im) / d;
    }
  }

  varargout_1[3].re = (node[3].re + re) + a_re;
  varargout_1[3].im = (node[3].im + im) + a_im;
}

double b_collisionKnowTf(const double parent_data[], const double node_data[],
  const creal_T *Tf)
{
  double collision_flag;
  creal_T t[12];
  int i;
  int state_tmp;
  creal_T dcv1[4];
  int exitg3;
  double b_state[40];
  int exitg2;
  double c_state;
  int exitg1;
  int state_idx_0_tmp;
  static const signed char iv9[5] = { 6, 15, 10, 5, 15 };

  double state_idx_0;
  static const signed char iv10[5] = { 6, 14, 11, 15, 5 };

  collision_flag = 0.0;
  linspace(*Tf, t);
  for (i = 0; i < 10; i++) {
    state_traj(t[i + 1], *Tf, parent_data[0], parent_data[1], parent_data[2],
               parent_data[3], node_data[0], node_data[1], node_data[2],
               node_data[3], dcv1);
    state_tmp = i << 2;
    b_state[state_tmp] = dcv1[0].re;
    b_state[1 + state_tmp] = dcv1[1].re;
    b_state[2 + state_tmp] = dcv1[2].re;
    b_state[3 + state_tmp] = dcv1[3].re;
  }

  state_tmp = 0;
  do {
    exitg3 = 0;
    if (state_tmp < 10) {
      i = 0;
      do {
        exitg2 = 0;
        if (i < 2) {
          c_state = b_state[i + (state_tmp << 2)];
          if ((c_state > 20.0) || (c_state < 0.0)) {
            collision_flag = 1.0;
            exitg2 = 1;
          } else {
            i++;
          }
        } else {
          /*  check each obstacle */
          i = 0;
          exitg2 = 2;
        }
      } while (exitg2 == 0);

      if (exitg2 == 1) {
        exitg3 = 1;
      } else {
        do {
          exitg1 = 0;
          if (i < 5) {
            state_idx_0_tmp = state_tmp << 2;
            c_state = b_state[state_idx_0_tmp] - (double)iv9[i];
            c_state *= c_state;
            state_idx_0 = c_state;
            c_state = b_state[1 + state_idx_0_tmp] - (double)iv10[i];
            c_state *= c_state;
            if (state_idx_0 + c_state < (dv0[i] + 0.1) * (dv0[i] + 0.1)) {
              /*  (norm([p(1);p(2)]-[world.cx(i); world.cy(i)])<=1*world.radius(i)) */
              collision_flag = 1.0;
              exitg1 = 1;
            } else {
              i++;
            }
          } else {
            state_tmp++;
            exitg1 = 2;
          }
        } while (exitg1 == 0);

        if (exitg1 == 1) {
          exitg3 = 1;
        }
      }
    } else {
      /* %%%% dim=3 case has not been implemented yet %%%%% */
      exitg3 = 1;
    }
  } while (exitg3 == 0);

  return collision_flag;
}

double collisionKnowTf(const double parent_data[], const creal_T node[8], const
  creal_T *Tf)
{
  double collision_flag;
  creal_T t[12];
  int i;
  int state_tmp;
  creal_T dcv0[4];
  int exitg3;
  double b_state[40];
  int exitg2;
  double c_state;
  int exitg1;
  int state_idx_0_tmp;
  static const signed char iv7[5] = { 6, 15, 10, 5, 15 };

  double state_idx_0;
  static const signed char iv8[5] = { 6, 14, 11, 15, 5 };

  collision_flag = 0.0;
  linspace(*Tf, t);
  for (i = 0; i < 10; i++) {
    __anon_fcn(*Tf, parent_data, node, t[i + 1], dcv0);
    state_tmp = i << 2;
    b_state[state_tmp] = dcv0[0].re;
    b_state[1 + state_tmp] = dcv0[1].re;
    b_state[2 + state_tmp] = dcv0[2].re;
    b_state[3 + state_tmp] = dcv0[3].re;
  }

  state_tmp = 0;
  do {
    exitg3 = 0;
    if (state_tmp < 10) {
      i = 0;
      do {
        exitg2 = 0;
        if (i < 2) {
          c_state = b_state[i + (state_tmp << 2)];
          if ((c_state > 20.0) || (c_state < 0.0)) {
            collision_flag = 1.0;
            exitg2 = 1;
          } else {
            i++;
          }
        } else {
          /*  check each obstacle */
          i = 0;
          exitg2 = 2;
        }
      } while (exitg2 == 0);

      if (exitg2 == 1) {
        exitg3 = 1;
      } else {
        do {
          exitg1 = 0;
          if (i < 5) {
            state_idx_0_tmp = state_tmp << 2;
            c_state = b_state[state_idx_0_tmp] - (double)iv7[i];
            c_state *= c_state;
            state_idx_0 = c_state;
            c_state = b_state[1 + state_idx_0_tmp] - (double)iv8[i];
            c_state *= c_state;
            if (state_idx_0 + c_state < (dv0[i] + 0.1) * (dv0[i] + 0.1)) {
              /*  (norm([p(1);p(2)]-[world.cx(i); world.cy(i)])<=1*world.radius(i)) */
              collision_flag = 1.0;
              exitg1 = 1;
            } else {
              i++;
            }
          } else {
            state_tmp++;
            exitg1 = 2;
          }
        } while (exitg1 == 0);

        if (exitg1 == 1) {
          exitg3 = 1;
        }
      }
    } else {
      /* %%%% dim=3 case has not been implemented yet %%%%% */
      exitg3 = 1;
    }
  } while (exitg3 == 0);

  return collision_flag;
}

/* End of code generation (collisionKnowTf.c) */
