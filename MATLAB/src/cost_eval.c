/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * cost_eval.c
 *
 * Code generation for function 'cost_eval'
 *
 */

/* Include files */
#include <math.h>
#include "rt_nonfinite.h"
#include "benchmarkRRT.h"
#include "cost_eval.h"
#include "benchmarkRRT_rtwutil.h"

/* Function Definitions */
creal_T cost_eval(const creal_T t_s, double x01, double x02, double x03, double
                  x04, double x11, double x12, double x13, double x14)
{
  creal_T cost;
  double y_re;
  double y_im;
  double x_re;
  double x_im;
  double b_y_re;
  double b_y_im;
  double c_y_re;
  double c_y_im;
  double r;
  double d_y_re;
  double d_y_im;
  double e_y_re;
  double e_y_im;
  double ar;
  double re;
  double brm;
  double im;
  double bim;
  double d;
  double b_re;
  double b_im;
  double c_re;
  double c_im;
  double d_re;
  double d_im;
  double e_re;
  double e_im;
  double f_re;
  double f_im;
  double g_re;
  double g_im;
  double h_re;
  double h_im;
  double t_s_im_tmp;
  double i_re;
  double i_im;
  double j_re;
  double j_im;
  double k_re;
  double k_im;
  double l_re;
  double l_im;
  double m_re;
  double m_im;
  double n_re;
  double n_im;
  if ((t_s.im == 0.0) && (t_s.re >= 0.0)) {
    y_re = rt_powd_snf(t_s.re, 3.0);
    y_im = 0.0;
  } else if (t_s.re == 0.0) {
    y_re = 0.0;
    y_im = -rt_powd_snf(t_s.im, 3.0);
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

    y_re = 3.0 * x_re;
    y_im = 3.0 * x_im;
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

  if ((t_s.im == 0.0) && (t_s.re >= 0.0)) {
    c_y_re = rt_powd_snf(t_s.re, 3.0);
    c_y_im = 0.0;
  } else if (t_s.re == 0.0) {
    c_y_re = 0.0;
    c_y_im = -rt_powd_snf(t_s.im, 3.0);
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

    c_y_re = 3.0 * x_re;
    c_y_im = 3.0 * x_im;
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

  ar = 4.0 * (x03 * x03);
  if (t_s.im == 0.0) {
    re = ar / t_s.re;
    im = 0.0;
  } else if (t_s.re == 0.0) {
    if (ar == 0.0) {
      re = 0.0 / t_s.im;
      im = 0.0;
    } else {
      re = 0.0;
      im = -(ar / t_s.im);
    }
  } else {
    brm = fabs(t_s.re);
    bim = fabs(t_s.im);
    if (brm > bim) {
      bim = t_s.im / t_s.re;
      d = t_s.re + bim * t_s.im;
      re = (ar + bim * 0.0) / d;
      im = (0.0 - bim * ar) / d;
    } else if (bim == brm) {
      if (t_s.re > 0.0) {
        bim = 0.5;
      } else {
        bim = -0.5;
      }

      if (t_s.im > 0.0) {
        d = 0.5;
      } else {
        d = -0.5;
      }

      re = (ar * bim + 0.0 * d) / brm;
      im = (0.0 * bim - ar * d) / brm;
    } else {
      bim = t_s.re / t_s.im;
      d = t_s.im + bim * t_s.re;
      re = bim * ar / d;
      im = (bim * 0.0 - ar) / d;
    }
  }

  ar = 12.0 * (x01 * x01);
  if (y_im == 0.0) {
    b_re = ar / y_re;
    b_im = 0.0;
  } else if (y_re == 0.0) {
    if (ar == 0.0) {
      b_re = 0.0 / y_im;
      b_im = 0.0;
    } else {
      b_re = 0.0;
      b_im = -(ar / y_im);
    }
  } else {
    brm = fabs(y_re);
    bim = fabs(y_im);
    if (brm > bim) {
      bim = y_im / y_re;
      d = y_re + bim * y_im;
      b_re = (ar + bim * 0.0) / d;
      b_im = (0.0 - bim * ar) / d;
    } else if (bim == brm) {
      if (y_re > 0.0) {
        bim = 0.5;
      } else {
        bim = -0.5;
      }

      if (y_im > 0.0) {
        d = 0.5;
      } else {
        d = -0.5;
      }

      b_re = (ar * bim + 0.0 * d) / brm;
      b_im = (0.0 * bim - ar * d) / brm;
    } else {
      bim = y_re / y_im;
      d = y_im + bim * y_re;
      b_re = bim * ar / d;
      b_im = (bim * 0.0 - ar) / d;
    }
  }

  ar = 4.0 * (x04 * x04);
  if (t_s.im == 0.0) {
    c_re = ar / t_s.re;
    c_im = 0.0;
  } else if (t_s.re == 0.0) {
    if (ar == 0.0) {
      c_re = 0.0 / t_s.im;
      c_im = 0.0;
    } else {
      c_re = 0.0;
      c_im = -(ar / t_s.im);
    }
  } else {
    brm = fabs(t_s.re);
    bim = fabs(t_s.im);
    if (brm > bim) {
      bim = t_s.im / t_s.re;
      d = t_s.re + bim * t_s.im;
      c_re = (ar + bim * 0.0) / d;
      c_im = (0.0 - bim * ar) / d;
    } else if (bim == brm) {
      if (t_s.re > 0.0) {
        bim = 0.5;
      } else {
        bim = -0.5;
      }

      if (t_s.im > 0.0) {
        d = 0.5;
      } else {
        d = -0.5;
      }

      c_re = (ar * bim + 0.0 * d) / brm;
      c_im = (0.0 * bim - ar * d) / brm;
    } else {
      bim = t_s.re / t_s.im;
      d = t_s.im + bim * t_s.re;
      c_re = bim * ar / d;
      c_im = (bim * 0.0 - ar) / d;
    }
  }

  ar = 12.0 * (x02 * x02);
  if (b_y_im == 0.0) {
    d_re = ar / b_y_re;
    d_im = 0.0;
  } else if (b_y_re == 0.0) {
    if (ar == 0.0) {
      d_re = 0.0 / b_y_im;
      d_im = 0.0;
    } else {
      d_re = 0.0;
      d_im = -(ar / b_y_im);
    }
  } else {
    brm = fabs(b_y_re);
    bim = fabs(b_y_im);
    if (brm > bim) {
      bim = b_y_im / b_y_re;
      d = b_y_re + bim * b_y_im;
      d_re = (ar + bim * 0.0) / d;
      d_im = (0.0 - bim * ar) / d;
    } else if (bim == brm) {
      if (b_y_re > 0.0) {
        bim = 0.5;
      } else {
        bim = -0.5;
      }

      if (b_y_im > 0.0) {
        d = 0.5;
      } else {
        d = -0.5;
      }

      d_re = (ar * bim + 0.0 * d) / brm;
      d_im = (0.0 * bim - ar * d) / brm;
    } else {
      bim = b_y_re / b_y_im;
      d = b_y_im + bim * b_y_re;
      d_re = bim * ar / d;
      d_im = (bim * 0.0 - ar) / d;
    }
  }

  ar = 4.0 * (x13 * x13);
  if (t_s.im == 0.0) {
    e_re = ar / t_s.re;
    e_im = 0.0;
  } else if (t_s.re == 0.0) {
    if (ar == 0.0) {
      e_re = 0.0 / t_s.im;
      e_im = 0.0;
    } else {
      e_re = 0.0;
      e_im = -(ar / t_s.im);
    }
  } else {
    brm = fabs(t_s.re);
    bim = fabs(t_s.im);
    if (brm > bim) {
      bim = t_s.im / t_s.re;
      d = t_s.re + bim * t_s.im;
      e_re = (ar + bim * 0.0) / d;
      e_im = (0.0 - bim * ar) / d;
    } else if (bim == brm) {
      if (t_s.re > 0.0) {
        bim = 0.5;
      } else {
        bim = -0.5;
      }

      if (t_s.im > 0.0) {
        d = 0.5;
      } else {
        d = -0.5;
      }

      e_re = (ar * bim + 0.0 * d) / brm;
      e_im = (0.0 * bim - ar * d) / brm;
    } else {
      bim = t_s.re / t_s.im;
      d = t_s.im + bim * t_s.re;
      e_re = bim * ar / d;
      e_im = (bim * 0.0 - ar) / d;
    }
  }

  ar = 12.0 * (x11 * x11);
  if (c_y_im == 0.0) {
    f_re = ar / c_y_re;
    f_im = 0.0;
  } else if (c_y_re == 0.0) {
    if (ar == 0.0) {
      f_re = 0.0 / c_y_im;
      f_im = 0.0;
    } else {
      f_re = 0.0;
      f_im = -(ar / c_y_im);
    }
  } else {
    brm = fabs(c_y_re);
    bim = fabs(c_y_im);
    if (brm > bim) {
      bim = c_y_im / c_y_re;
      d = c_y_re + bim * c_y_im;
      f_re = (ar + bim * 0.0) / d;
      f_im = (0.0 - bim * ar) / d;
    } else if (bim == brm) {
      if (c_y_re > 0.0) {
        bim = 0.5;
      } else {
        bim = -0.5;
      }

      if (c_y_im > 0.0) {
        d = 0.5;
      } else {
        d = -0.5;
      }

      f_re = (ar * bim + 0.0 * d) / brm;
      f_im = (0.0 * bim - ar * d) / brm;
    } else {
      bim = c_y_re / c_y_im;
      d = c_y_im + bim * c_y_re;
      f_re = bim * ar / d;
      f_im = (bim * 0.0 - ar) / d;
    }
  }

  ar = 4.0 * (x14 * x14);
  if (t_s.im == 0.0) {
    g_re = ar / t_s.re;
    g_im = 0.0;
  } else if (t_s.re == 0.0) {
    if (ar == 0.0) {
      g_re = 0.0 / t_s.im;
      g_im = 0.0;
    } else {
      g_re = 0.0;
      g_im = -(ar / t_s.im);
    }
  } else {
    brm = fabs(t_s.re);
    bim = fabs(t_s.im);
    if (brm > bim) {
      bim = t_s.im / t_s.re;
      d = t_s.re + bim * t_s.im;
      g_re = (ar + bim * 0.0) / d;
      g_im = (0.0 - bim * ar) / d;
    } else if (bim == brm) {
      if (t_s.re > 0.0) {
        bim = 0.5;
      } else {
        bim = -0.5;
      }

      if (t_s.im > 0.0) {
        d = 0.5;
      } else {
        d = -0.5;
      }

      g_re = (ar * bim + 0.0 * d) / brm;
      g_im = (0.0 * bim - ar * d) / brm;
    } else {
      bim = t_s.re / t_s.im;
      d = t_s.im + bim * t_s.re;
      g_re = bim * ar / d;
      g_im = (bim * 0.0 - ar) / d;
    }
  }

  ar = 12.0 * (x12 * x12);
  if (d_y_im == 0.0) {
    h_re = ar / d_y_re;
    h_im = 0.0;
  } else if (d_y_re == 0.0) {
    if (ar == 0.0) {
      h_re = 0.0 / d_y_im;
      h_im = 0.0;
    } else {
      h_re = 0.0;
      h_im = -(ar / d_y_im);
    }
  } else {
    brm = fabs(d_y_re);
    bim = fabs(d_y_im);
    if (brm > bim) {
      bim = d_y_im / d_y_re;
      d = d_y_re + bim * d_y_im;
      h_re = (ar + bim * 0.0) / d;
      h_im = (0.0 - bim * ar) / d;
    } else if (bim == brm) {
      if (d_y_re > 0.0) {
        bim = 0.5;
      } else {
        bim = -0.5;
      }

      if (d_y_im > 0.0) {
        d = 0.5;
      } else {
        d = -0.5;
      }

      h_re = (ar * bim + 0.0 * d) / brm;
      h_im = (0.0 * bim - ar * d) / brm;
    } else {
      bim = d_y_re / d_y_im;
      d = d_y_im + bim * d_y_re;
      h_re = bim * ar / d;
      h_im = (bim * 0.0 - ar) / d;
    }
  }

  r = t_s.re * t_s.re - t_s.im * t_s.im;
  t_s_im_tmp = t_s.re * t_s.im + t_s.im * t_s.re;
  ar = 12.0 * x01 * x03;
  if (t_s_im_tmp == 0.0) {
    i_re = ar / r;
    i_im = 0.0;
  } else if (r == 0.0) {
    if (ar == 0.0) {
      i_re = 0.0 / t_s_im_tmp;
      i_im = 0.0;
    } else {
      i_re = 0.0;
      i_im = -(ar / t_s_im_tmp);
    }
  } else {
    brm = fabs(r);
    bim = fabs(t_s_im_tmp);
    if (brm > bim) {
      bim = t_s_im_tmp / r;
      d = r + bim * t_s_im_tmp;
      i_re = (ar + bim * 0.0) / d;
      i_im = (0.0 - bim * ar) / d;
    } else if (bim == brm) {
      if (r > 0.0) {
        bim = 0.5;
      } else {
        bim = -0.5;
      }

      if (t_s_im_tmp > 0.0) {
        d = 0.5;
      } else {
        d = -0.5;
      }

      i_re = (ar * bim + 0.0 * d) / brm;
      i_im = (0.0 * bim - ar * d) / brm;
    } else {
      bim = r / t_s_im_tmp;
      d = t_s_im_tmp + bim * r;
      i_re = bim * ar / d;
      i_im = (bim * 0.0 - ar) / d;
    }
  }

  ar = 12.0 * x02 * x04;
  if (t_s_im_tmp == 0.0) {
    j_re = ar / r;
    j_im = 0.0;
  } else if (r == 0.0) {
    if (ar == 0.0) {
      j_re = 0.0 / t_s_im_tmp;
      j_im = 0.0;
    } else {
      j_re = 0.0;
      j_im = -(ar / t_s_im_tmp);
    }
  } else {
    brm = fabs(r);
    bim = fabs(t_s_im_tmp);
    if (brm > bim) {
      bim = t_s_im_tmp / r;
      d = r + bim * t_s_im_tmp;
      j_re = (ar + bim * 0.0) / d;
      j_im = (0.0 - bim * ar) / d;
    } else if (bim == brm) {
      if (r > 0.0) {
        bim = 0.5;
      } else {
        bim = -0.5;
      }

      if (t_s_im_tmp > 0.0) {
        d = 0.5;
      } else {
        d = -0.5;
      }

      j_re = (ar * bim + 0.0 * d) / brm;
      j_im = (0.0 * bim - ar * d) / brm;
    } else {
      bim = r / t_s_im_tmp;
      d = t_s_im_tmp + bim * r;
      j_re = bim * ar / d;
      j_im = (bim * 0.0 - ar) / d;
    }
  }

  ar = 24.0 * x01 * x11;
  if (e_y_im == 0.0) {
    k_re = ar / e_y_re;
    k_im = 0.0;
  } else if (e_y_re == 0.0) {
    if (ar == 0.0) {
      k_re = 0.0 / e_y_im;
      k_im = 0.0;
    } else {
      k_re = 0.0;
      k_im = -(ar / e_y_im);
    }
  } else {
    brm = fabs(e_y_re);
    bim = fabs(e_y_im);
    if (brm > bim) {
      bim = e_y_im / e_y_re;
      d = e_y_re + bim * e_y_im;
      k_re = (ar + bim * 0.0) / d;
      k_im = (0.0 - bim * ar) / d;
    } else if (bim == brm) {
      if (e_y_re > 0.0) {
        bim = 0.5;
      } else {
        bim = -0.5;
      }

      if (e_y_im > 0.0) {
        d = 0.5;
      } else {
        d = -0.5;
      }

      k_re = (ar * bim + 0.0 * d) / brm;
      k_im = (0.0 * bim - ar * d) / brm;
    } else {
      bim = e_y_re / e_y_im;
      d = e_y_im + bim * e_y_re;
      k_re = bim * ar / d;
      k_im = (bim * 0.0 - ar) / d;
    }
  }

  ar = 12.0 * x01 * x13;
  if (t_s_im_tmp == 0.0) {
    l_re = ar / r;
    l_im = 0.0;
  } else if (r == 0.0) {
    if (ar == 0.0) {
      l_re = 0.0 / t_s_im_tmp;
      l_im = 0.0;
    } else {
      l_re = 0.0;
      l_im = -(ar / t_s_im_tmp);
    }
  } else {
    brm = fabs(r);
    bim = fabs(t_s_im_tmp);
    if (brm > bim) {
      bim = t_s_im_tmp / r;
      d = r + bim * t_s_im_tmp;
      l_re = (ar + bim * 0.0) / d;
      l_im = (0.0 - bim * ar) / d;
    } else if (bim == brm) {
      if (r > 0.0) {
        bim = 0.5;
      } else {
        bim = -0.5;
      }

      if (t_s_im_tmp > 0.0) {
        d = 0.5;
      } else {
        d = -0.5;
      }

      l_re = (ar * bim + 0.0 * d) / brm;
      l_im = (0.0 * bim - ar * d) / brm;
    } else {
      bim = r / t_s_im_tmp;
      d = t_s_im_tmp + bim * r;
      l_re = bim * ar / d;
      l_im = (bim * 0.0 - ar) / d;
    }
  }

  ar = 12.0 * x03 * x11;
  if (t_s_im_tmp == 0.0) {
    m_re = ar / r;
    m_im = 0.0;
  } else if (r == 0.0) {
    if (ar == 0.0) {
      m_re = 0.0 / t_s_im_tmp;
      m_im = 0.0;
    } else {
      m_re = 0.0;
      m_im = -(ar / t_s_im_tmp);
    }
  } else {
    brm = fabs(r);
    bim = fabs(t_s_im_tmp);
    if (brm > bim) {
      bim = t_s_im_tmp / r;
      d = r + bim * t_s_im_tmp;
      m_re = (ar + bim * 0.0) / d;
      m_im = (0.0 - bim * ar) / d;
    } else if (bim == brm) {
      if (r > 0.0) {
        bim = 0.5;
      } else {
        bim = -0.5;
      }

      if (t_s_im_tmp > 0.0) {
        d = 0.5;
      } else {
        d = -0.5;
      }

      m_re = (ar * bim + 0.0 * d) / brm;
      m_im = (0.0 * bim - ar * d) / brm;
    } else {
      bim = r / t_s_im_tmp;
      d = t_s_im_tmp + bim * r;
      m_re = bim * ar / d;
      m_im = (bim * 0.0 - ar) / d;
    }
  }

  ar = 4.0 * x03 * x13;
  if (t_s.im == 0.0) {
    n_re = ar / t_s.re;
    n_im = 0.0;
  } else if (t_s.re == 0.0) {
    if (ar == 0.0) {
      n_re = 0.0 / t_s.im;
      n_im = 0.0;
    } else {
      n_re = 0.0;
      n_im = -(ar / t_s.im);
    }
  } else {
    brm = fabs(t_s.re);
    bim = fabs(t_s.im);
    if (brm > bim) {
      bim = t_s.im / t_s.re;
      d = t_s.re + bim * t_s.im;
      n_re = (ar + bim * 0.0) / d;
      n_im = (0.0 - bim * ar) / d;
    } else if (bim == brm) {
      if (t_s.re > 0.0) {
        bim = 0.5;
      } else {
        bim = -0.5;
      }

      if (t_s.im > 0.0) {
        d = 0.5;
      } else {
        d = -0.5;
      }

      n_re = (ar * bim + 0.0 * d) / brm;
      n_im = (0.0 * bim - ar * d) / brm;
    } else {
      bim = t_s.re / t_s.im;
      d = t_s.im + bim * t_s.re;
      n_re = bim * ar / d;
      n_im = (bim * 0.0 - ar) / d;
    }
  }

  ar = 24.0 * x02 * x12;
  if (x_im == 0.0) {
    x_re = ar / x_re;
    y_im = 0.0;
  } else if (x_re == 0.0) {
    if (ar == 0.0) {
      x_re = 0.0 / x_im;
      y_im = 0.0;
    } else {
      x_re = 0.0;
      y_im = -(ar / x_im);
    }
  } else {
    brm = fabs(x_re);
    bim = fabs(x_im);
    if (brm > bim) {
      bim = x_im / x_re;
      d = x_re + bim * x_im;
      x_re = (ar + bim * 0.0) / d;
      y_im = (0.0 - bim * ar) / d;
    } else if (bim == brm) {
      if (x_re > 0.0) {
        bim = 0.5;
      } else {
        bim = -0.5;
      }

      if (x_im > 0.0) {
        d = 0.5;
      } else {
        d = -0.5;
      }

      x_re = (ar * bim + 0.0 * d) / brm;
      y_im = (0.0 * bim - ar * d) / brm;
    } else {
      bim = x_re / x_im;
      d = x_im + bim * x_re;
      x_re = bim * ar / d;
      y_im = (bim * 0.0 - ar) / d;
    }
  }

  ar = 12.0 * x02 * x14;
  if (t_s_im_tmp == 0.0) {
    b_y_re = ar / r;
    b_y_im = 0.0;
  } else if (r == 0.0) {
    if (ar == 0.0) {
      b_y_re = 0.0 / t_s_im_tmp;
      b_y_im = 0.0;
    } else {
      b_y_re = 0.0;
      b_y_im = -(ar / t_s_im_tmp);
    }
  } else {
    brm = fabs(r);
    bim = fabs(t_s_im_tmp);
    if (brm > bim) {
      bim = t_s_im_tmp / r;
      d = r + bim * t_s_im_tmp;
      b_y_re = (ar + bim * 0.0) / d;
      b_y_im = (0.0 - bim * ar) / d;
    } else if (bim == brm) {
      if (r > 0.0) {
        bim = 0.5;
      } else {
        bim = -0.5;
      }

      if (t_s_im_tmp > 0.0) {
        d = 0.5;
      } else {
        d = -0.5;
      }

      b_y_re = (ar * bim + 0.0 * d) / brm;
      b_y_im = (0.0 * bim - ar * d) / brm;
    } else {
      bim = r / t_s_im_tmp;
      d = t_s_im_tmp + bim * r;
      b_y_re = bim * ar / d;
      b_y_im = (bim * 0.0 - ar) / d;
    }
  }

  ar = 12.0 * x04 * x12;
  if (t_s_im_tmp == 0.0) {
    c_y_re = ar / r;
    c_y_im = 0.0;
  } else if (r == 0.0) {
    if (ar == 0.0) {
      c_y_re = 0.0 / t_s_im_tmp;
      c_y_im = 0.0;
    } else {
      c_y_re = 0.0;
      c_y_im = -(ar / t_s_im_tmp);
    }
  } else {
    brm = fabs(r);
    bim = fabs(t_s_im_tmp);
    if (brm > bim) {
      bim = t_s_im_tmp / r;
      d = r + bim * t_s_im_tmp;
      c_y_re = (ar + bim * 0.0) / d;
      c_y_im = (0.0 - bim * ar) / d;
    } else if (bim == brm) {
      if (r > 0.0) {
        bim = 0.5;
      } else {
        bim = -0.5;
      }

      if (t_s_im_tmp > 0.0) {
        d = 0.5;
      } else {
        d = -0.5;
      }

      c_y_re = (ar * bim + 0.0 * d) / brm;
      c_y_im = (0.0 * bim - ar * d) / brm;
    } else {
      bim = r / t_s_im_tmp;
      d = t_s_im_tmp + bim * r;
      c_y_re = bim * ar / d;
      c_y_im = (bim * 0.0 - ar) / d;
    }
  }

  ar = 4.0 * x04 * x14;
  if (t_s.im == 0.0) {
    d_y_re = ar / t_s.re;
    d_y_im = 0.0;
  } else if (t_s.re == 0.0) {
    if (ar == 0.0) {
      d_y_re = 0.0 / t_s.im;
      d_y_im = 0.0;
    } else {
      d_y_re = 0.0;
      d_y_im = -(ar / t_s.im);
    }
  } else {
    brm = fabs(t_s.re);
    bim = fabs(t_s.im);
    if (brm > bim) {
      bim = t_s.im / t_s.re;
      d = t_s.re + bim * t_s.im;
      d_y_re = (ar + bim * 0.0) / d;
      d_y_im = (0.0 - bim * ar) / d;
    } else if (bim == brm) {
      if (t_s.re > 0.0) {
        bim = 0.5;
      } else {
        bim = -0.5;
      }

      if (t_s.im > 0.0) {
        d = 0.5;
      } else {
        d = -0.5;
      }

      d_y_re = (ar * bim + 0.0 * d) / brm;
      d_y_im = (0.0 * bim - ar * d) / brm;
    } else {
      bim = t_s.re / t_s.im;
      d = t_s.im + bim * t_s.re;
      d_y_re = bim * ar / d;
      d_y_im = (bim * 0.0 - ar) / d;
    }
  }

  ar = 12.0 * x11 * x13;
  if (t_s_im_tmp == 0.0) {
    e_y_re = ar / r;
    e_y_im = 0.0;
  } else if (r == 0.0) {
    if (ar == 0.0) {
      e_y_re = 0.0 / t_s_im_tmp;
      e_y_im = 0.0;
    } else {
      e_y_re = 0.0;
      e_y_im = -(ar / t_s_im_tmp);
    }
  } else {
    brm = fabs(r);
    bim = fabs(t_s_im_tmp);
    if (brm > bim) {
      bim = t_s_im_tmp / r;
      d = r + bim * t_s_im_tmp;
      e_y_re = (ar + bim * 0.0) / d;
      e_y_im = (0.0 - bim * ar) / d;
    } else if (bim == brm) {
      if (r > 0.0) {
        bim = 0.5;
      } else {
        bim = -0.5;
      }

      if (t_s_im_tmp > 0.0) {
        d = 0.5;
      } else {
        d = -0.5;
      }

      e_y_re = (ar * bim + 0.0 * d) / brm;
      e_y_im = (0.0 * bim - ar * d) / brm;
    } else {
      bim = r / t_s_im_tmp;
      d = t_s_im_tmp + bim * r;
      e_y_re = bim * ar / d;
      e_y_im = (bim * 0.0 - ar) / d;
    }
  }

  ar = 12.0 * x12 * x14;
  if (t_s_im_tmp == 0.0) {
    y_re = ar / r;
    r = 0.0;
  } else if (r == 0.0) {
    if (ar == 0.0) {
      y_re = 0.0 / t_s_im_tmp;
      r = 0.0;
    } else {
      y_re = 0.0;
      r = -(ar / t_s_im_tmp);
    }
  } else {
    brm = fabs(r);
    bim = fabs(t_s_im_tmp);
    if (brm > bim) {
      bim = t_s_im_tmp / r;
      d = r + bim * t_s_im_tmp;
      y_re = (ar + bim * 0.0) / d;
      r = (0.0 - bim * ar) / d;
    } else if (bim == brm) {
      if (r > 0.0) {
        bim = 0.5;
      } else {
        bim = -0.5;
      }

      if (t_s_im_tmp > 0.0) {
        d = 0.5;
      } else {
        d = -0.5;
      }

      y_re = (ar * bim + 0.0 * d) / brm;
      r = (0.0 * bim - ar * d) / brm;
    } else {
      bim = r / t_s_im_tmp;
      d = t_s_im_tmp + bim * r;
      y_re = bim * ar / d;
      r = (bim * 0.0 - ar) / d;
    }
  }

  cost.re = (((((((((((((((((((t_s.re + re) + b_re) + c_re) + d_re) + e_re) +
    f_re) + g_re) + h_re) + i_re) + j_re) - k_re) + l_re) - m_re) + n_re) - x_re)
                + b_y_re) - c_y_re) + d_y_re) - e_y_re) - y_re;
  cost.im = (((((((((((((((((((t_s.im + im) + b_im) + c_im) + d_im) + e_im) +
    f_im) + g_im) + h_im) + i_im) + j_im) - k_im) + l_im) - m_im) + n_im) - y_im)
                + b_y_im) - c_y_im) + d_y_im) - e_y_im) - r;
  return cost;
}

/* End of code generation (cost_eval.c) */
