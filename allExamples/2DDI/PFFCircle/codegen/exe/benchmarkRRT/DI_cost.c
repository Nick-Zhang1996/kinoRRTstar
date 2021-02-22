/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * DI_cost.c
 *
 * Code generation for function 'DI_cost'
 *
 */

/* Include files */
#include <math.h>
#include "rt_nonfinite.h"
#include "benchmarkRRT.h"
#include "DI_cost.h"
#include "benchmarkRRT_rtwutil.h"

/* Function Definitions */
double DI_cost(const creal_T t_s, double x01, double x02, double x03, double x04,
               double x11, double x12, double x13, double x14)
{
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
  double bim;
  double b_re;
  double b_t_s;
  double c_re;
  double c_t_s;
  double f_y_re;
  double d_re;
  double f_y_im;
  double d_t_s;
  double e_re;
  double e_t_s;
  double g_y_re;
  double f_re;
  double g_y_im;
  double f_t_s;
  double g_re;
  double g_t_s;
  double h_y_re;
  double h_re;
  double h_y_im;
  double h_t_s;
  double t_s_im_tmp;
  double i_t_s;
  double i_y_re;
  double i_re;
  double i_y_im;
  double j_re;
  double b_r;
  double b_t_s_im_tmp;
  double c_r;
  double c_t_s_im_tmp;
  double j_y_re;
  double j_y_im;
  double d_r;
  double d_t_s_im_tmp;
  double e_r;
  double e_t_s_im_tmp;
  double j_t_s;
  double k_t_s;
  double b_x_re;
  double b_x_im;
  double f_r;
  double f_t_s_im_tmp;
  double g_r;
  double g_t_s_im_tmp;
  double l_t_s;
  double m_t_s;
  double h_r;
  double h_t_s_im_tmp;
  double i_r;
  double i_t_s_im_tmp;
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
  } else if (t_s.re == 0.0) {
    if (ar == 0.0) {
      re = 0.0 / t_s.im;
    } else {
      re = 0.0;
    }
  } else {
    brm = fabs(t_s.re);
    bim = fabs(t_s.im);
    if (brm > bim) {
      brm = t_s.im / t_s.re;
      re = (ar + brm * 0.0) / (t_s.re + brm * t_s.im);
    } else if (bim == brm) {
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

      re = (ar * b_t_s + 0.0 * c_t_s) / brm;
    } else {
      brm = t_s.re / t_s.im;
      re = brm * ar / (t_s.im + brm * t_s.re);
    }
  }

  ar = 12.0 * (x01 * x01);
  if (y_im == 0.0) {
    b_re = ar / y_re;
  } else if (y_re == 0.0) {
    if (ar == 0.0) {
      b_re = 0.0 / y_im;
    } else {
      b_re = 0.0;
    }
  } else {
    brm = fabs(y_re);
    bim = fabs(y_im);
    if (brm > bim) {
      brm = y_im / y_re;
      b_re = (ar + brm * 0.0) / (y_re + brm * y_im);
    } else if (bim == brm) {
      if (y_re > 0.0) {
        f_y_re = 0.5;
      } else {
        f_y_re = -0.5;
      }

      if (y_im > 0.0) {
        f_y_im = 0.5;
      } else {
        f_y_im = -0.5;
      }

      b_re = (ar * f_y_re + 0.0 * f_y_im) / brm;
    } else {
      brm = y_re / y_im;
      b_re = brm * ar / (y_im + brm * y_re);
    }
  }

  ar = 4.0 * (x04 * x04);
  if (t_s.im == 0.0) {
    c_re = ar / t_s.re;
  } else if (t_s.re == 0.0) {
    if (ar == 0.0) {
      c_re = 0.0 / t_s.im;
    } else {
      c_re = 0.0;
    }
  } else {
    brm = fabs(t_s.re);
    bim = fabs(t_s.im);
    if (brm > bim) {
      brm = t_s.im / t_s.re;
      c_re = (ar + brm * 0.0) / (t_s.re + brm * t_s.im);
    } else if (bim == brm) {
      if (t_s.re > 0.0) {
        d_t_s = 0.5;
      } else {
        d_t_s = -0.5;
      }

      if (t_s.im > 0.0) {
        e_t_s = 0.5;
      } else {
        e_t_s = -0.5;
      }

      c_re = (ar * d_t_s + 0.0 * e_t_s) / brm;
    } else {
      brm = t_s.re / t_s.im;
      c_re = brm * ar / (t_s.im + brm * t_s.re);
    }
  }

  ar = 12.0 * (x02 * x02);
  if (b_y_im == 0.0) {
    d_re = ar / b_y_re;
  } else if (b_y_re == 0.0) {
    if (ar == 0.0) {
      d_re = 0.0 / b_y_im;
    } else {
      d_re = 0.0;
    }
  } else {
    brm = fabs(b_y_re);
    bim = fabs(b_y_im);
    if (brm > bim) {
      brm = b_y_im / b_y_re;
      d_re = (ar + brm * 0.0) / (b_y_re + brm * b_y_im);
    } else if (bim == brm) {
      if (b_y_re > 0.0) {
        g_y_re = 0.5;
      } else {
        g_y_re = -0.5;
      }

      if (b_y_im > 0.0) {
        g_y_im = 0.5;
      } else {
        g_y_im = -0.5;
      }

      d_re = (ar * g_y_re + 0.0 * g_y_im) / brm;
    } else {
      brm = b_y_re / b_y_im;
      d_re = brm * ar / (b_y_im + brm * b_y_re);
    }
  }

  ar = 4.0 * (x13 * x13);
  if (t_s.im == 0.0) {
    e_re = ar / t_s.re;
  } else if (t_s.re == 0.0) {
    if (ar == 0.0) {
      e_re = 0.0 / t_s.im;
    } else {
      e_re = 0.0;
    }
  } else {
    brm = fabs(t_s.re);
    bim = fabs(t_s.im);
    if (brm > bim) {
      brm = t_s.im / t_s.re;
      e_re = (ar + brm * 0.0) / (t_s.re + brm * t_s.im);
    } else if (bim == brm) {
      if (t_s.re > 0.0) {
        f_t_s = 0.5;
      } else {
        f_t_s = -0.5;
      }

      if (t_s.im > 0.0) {
        g_t_s = 0.5;
      } else {
        g_t_s = -0.5;
      }

      e_re = (ar * f_t_s + 0.0 * g_t_s) / brm;
    } else {
      brm = t_s.re / t_s.im;
      e_re = brm * ar / (t_s.im + brm * t_s.re);
    }
  }

  ar = 12.0 * (x11 * x11);
  if (c_y_im == 0.0) {
    f_re = ar / c_y_re;
  } else if (c_y_re == 0.0) {
    if (ar == 0.0) {
      f_re = 0.0 / c_y_im;
    } else {
      f_re = 0.0;
    }
  } else {
    brm = fabs(c_y_re);
    bim = fabs(c_y_im);
    if (brm > bim) {
      brm = c_y_im / c_y_re;
      f_re = (ar + brm * 0.0) / (c_y_re + brm * c_y_im);
    } else if (bim == brm) {
      if (c_y_re > 0.0) {
        h_y_re = 0.5;
      } else {
        h_y_re = -0.5;
      }

      if (c_y_im > 0.0) {
        h_y_im = 0.5;
      } else {
        h_y_im = -0.5;
      }

      f_re = (ar * h_y_re + 0.0 * h_y_im) / brm;
    } else {
      brm = c_y_re / c_y_im;
      f_re = brm * ar / (c_y_im + brm * c_y_re);
    }
  }

  ar = 4.0 * (x14 * x14);
  if (t_s.im == 0.0) {
    g_re = ar / t_s.re;
  } else if (t_s.re == 0.0) {
    if (ar == 0.0) {
      g_re = 0.0 / t_s.im;
    } else {
      g_re = 0.0;
    }
  } else {
    brm = fabs(t_s.re);
    bim = fabs(t_s.im);
    if (brm > bim) {
      brm = t_s.im / t_s.re;
      g_re = (ar + brm * 0.0) / (t_s.re + brm * t_s.im);
    } else if (bim == brm) {
      if (t_s.re > 0.0) {
        h_t_s = 0.5;
      } else {
        h_t_s = -0.5;
      }

      if (t_s.im > 0.0) {
        i_t_s = 0.5;
      } else {
        i_t_s = -0.5;
      }

      g_re = (ar * h_t_s + 0.0 * i_t_s) / brm;
    } else {
      brm = t_s.re / t_s.im;
      g_re = brm * ar / (t_s.im + brm * t_s.re);
    }
  }

  ar = 12.0 * (x12 * x12);
  if (d_y_im == 0.0) {
    h_re = ar / d_y_re;
  } else if (d_y_re == 0.0) {
    if (ar == 0.0) {
      h_re = 0.0 / d_y_im;
    } else {
      h_re = 0.0;
    }
  } else {
    brm = fabs(d_y_re);
    bim = fabs(d_y_im);
    if (brm > bim) {
      brm = d_y_im / d_y_re;
      h_re = (ar + brm * 0.0) / (d_y_re + brm * d_y_im);
    } else if (bim == brm) {
      if (d_y_re > 0.0) {
        i_y_re = 0.5;
      } else {
        i_y_re = -0.5;
      }

      if (d_y_im > 0.0) {
        i_y_im = 0.5;
      } else {
        i_y_im = -0.5;
      }

      h_re = (ar * i_y_re + 0.0 * i_y_im) / brm;
    } else {
      brm = d_y_re / d_y_im;
      h_re = brm * ar / (d_y_im + brm * d_y_re);
    }
  }

  r = t_s.re * t_s.re - t_s.im * t_s.im;
  t_s_im_tmp = t_s.re * t_s.im + t_s.im * t_s.re;
  ar = 12.0 * x01 * x03;
  if (t_s_im_tmp == 0.0) {
    i_re = ar / r;
  } else if (r == 0.0) {
    if (ar == 0.0) {
      i_re = 0.0 / t_s_im_tmp;
    } else {
      i_re = 0.0;
    }
  } else {
    brm = fabs(r);
    bim = fabs(t_s_im_tmp);
    if (brm > bim) {
      brm = t_s_im_tmp / r;
      i_re = (ar + brm * 0.0) / (r + brm * t_s_im_tmp);
    } else if (bim == brm) {
      if (r > 0.0) {
        b_r = 0.5;
      } else {
        b_r = -0.5;
      }

      if (t_s_im_tmp > 0.0) {
        b_t_s_im_tmp = 0.5;
      } else {
        b_t_s_im_tmp = -0.5;
      }

      i_re = (ar * b_r + 0.0 * b_t_s_im_tmp) / brm;
    } else {
      brm = r / t_s_im_tmp;
      i_re = brm * ar / (t_s_im_tmp + brm * r);
    }
  }

  ar = 12.0 * x02 * x04;
  if (t_s_im_tmp == 0.0) {
    j_re = ar / r;
  } else if (r == 0.0) {
    if (ar == 0.0) {
      j_re = 0.0 / t_s_im_tmp;
    } else {
      j_re = 0.0;
    }
  } else {
    brm = fabs(r);
    bim = fabs(t_s_im_tmp);
    if (brm > bim) {
      brm = t_s_im_tmp / r;
      j_re = (ar + brm * 0.0) / (r + brm * t_s_im_tmp);
    } else if (bim == brm) {
      if (r > 0.0) {
        c_r = 0.5;
      } else {
        c_r = -0.5;
      }

      if (t_s_im_tmp > 0.0) {
        c_t_s_im_tmp = 0.5;
      } else {
        c_t_s_im_tmp = -0.5;
      }

      j_re = (ar * c_r + 0.0 * c_t_s_im_tmp) / brm;
    } else {
      brm = r / t_s_im_tmp;
      j_re = brm * ar / (t_s_im_tmp + brm * r);
    }
  }

  ar = 24.0 * x01 * x11;
  if (e_y_im == 0.0) {
    c_y_im = ar / e_y_re;
  } else if (e_y_re == 0.0) {
    if (ar == 0.0) {
      c_y_im = 0.0 / e_y_im;
    } else {
      c_y_im = 0.0;
    }
  } else {
    brm = fabs(e_y_re);
    bim = fabs(e_y_im);
    if (brm > bim) {
      brm = e_y_im / e_y_re;
      c_y_im = (ar + brm * 0.0) / (e_y_re + brm * e_y_im);
    } else if (bim == brm) {
      if (e_y_re > 0.0) {
        j_y_re = 0.5;
      } else {
        j_y_re = -0.5;
      }

      if (e_y_im > 0.0) {
        j_y_im = 0.5;
      } else {
        j_y_im = -0.5;
      }

      c_y_im = (ar * j_y_re + 0.0 * j_y_im) / brm;
    } else {
      brm = e_y_re / e_y_im;
      c_y_im = brm * ar / (e_y_im + brm * e_y_re);
    }
  }

  ar = 12.0 * x01 * x13;
  if (t_s_im_tmp == 0.0) {
    d_y_re = ar / r;
  } else if (r == 0.0) {
    if (ar == 0.0) {
      d_y_re = 0.0 / t_s_im_tmp;
    } else {
      d_y_re = 0.0;
    }
  } else {
    brm = fabs(r);
    bim = fabs(t_s_im_tmp);
    if (brm > bim) {
      brm = t_s_im_tmp / r;
      d_y_re = (ar + brm * 0.0) / (r + brm * t_s_im_tmp);
    } else if (bim == brm) {
      if (r > 0.0) {
        d_r = 0.5;
      } else {
        d_r = -0.5;
      }

      if (t_s_im_tmp > 0.0) {
        d_t_s_im_tmp = 0.5;
      } else {
        d_t_s_im_tmp = -0.5;
      }

      d_y_re = (ar * d_r + 0.0 * d_t_s_im_tmp) / brm;
    } else {
      brm = r / t_s_im_tmp;
      d_y_re = brm * ar / (t_s_im_tmp + brm * r);
    }
  }

  ar = 12.0 * x03 * x11;
  if (t_s_im_tmp == 0.0) {
    d_y_im = ar / r;
  } else if (r == 0.0) {
    if (ar == 0.0) {
      d_y_im = 0.0 / t_s_im_tmp;
    } else {
      d_y_im = 0.0;
    }
  } else {
    brm = fabs(r);
    bim = fabs(t_s_im_tmp);
    if (brm > bim) {
      brm = t_s_im_tmp / r;
      d_y_im = (ar + brm * 0.0) / (r + brm * t_s_im_tmp);
    } else if (bim == brm) {
      if (r > 0.0) {
        e_r = 0.5;
      } else {
        e_r = -0.5;
      }

      if (t_s_im_tmp > 0.0) {
        e_t_s_im_tmp = 0.5;
      } else {
        e_t_s_im_tmp = -0.5;
      }

      d_y_im = (ar * e_r + 0.0 * e_t_s_im_tmp) / brm;
    } else {
      brm = r / t_s_im_tmp;
      d_y_im = brm * ar / (t_s_im_tmp + brm * r);
    }
  }

  ar = 4.0 * x03 * x13;
  if (t_s.im == 0.0) {
    e_y_re = ar / t_s.re;
  } else if (t_s.re == 0.0) {
    if (ar == 0.0) {
      e_y_re = 0.0 / t_s.im;
    } else {
      e_y_re = 0.0;
    }
  } else {
    brm = fabs(t_s.re);
    bim = fabs(t_s.im);
    if (brm > bim) {
      brm = t_s.im / t_s.re;
      e_y_re = (ar + brm * 0.0) / (t_s.re + brm * t_s.im);
    } else if (bim == brm) {
      if (t_s.re > 0.0) {
        j_t_s = 0.5;
      } else {
        j_t_s = -0.5;
      }

      if (t_s.im > 0.0) {
        k_t_s = 0.5;
      } else {
        k_t_s = -0.5;
      }

      e_y_re = (ar * j_t_s + 0.0 * k_t_s) / brm;
    } else {
      brm = t_s.re / t_s.im;
      e_y_re = brm * ar / (t_s.im + brm * t_s.re);
    }
  }

  ar = 24.0 * x02 * x12;
  if (x_im == 0.0) {
    y_re = ar / x_re;
  } else if (x_re == 0.0) {
    if (ar == 0.0) {
      y_re = 0.0 / x_im;
    } else {
      y_re = 0.0;
    }
  } else {
    brm = fabs(x_re);
    bim = fabs(x_im);
    if (brm > bim) {
      brm = x_im / x_re;
      y_re = (ar + brm * 0.0) / (x_re + brm * x_im);
    } else if (bim == brm) {
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

      y_re = (ar * b_x_re + 0.0 * b_x_im) / brm;
    } else {
      brm = x_re / x_im;
      y_re = brm * ar / (x_im + brm * x_re);
    }
  }

  ar = 12.0 * x02 * x14;
  if (t_s_im_tmp == 0.0) {
    y_im = ar / r;
  } else if (r == 0.0) {
    if (ar == 0.0) {
      y_im = 0.0 / t_s_im_tmp;
    } else {
      y_im = 0.0;
    }
  } else {
    brm = fabs(r);
    bim = fabs(t_s_im_tmp);
    if (brm > bim) {
      brm = t_s_im_tmp / r;
      y_im = (ar + brm * 0.0) / (r + brm * t_s_im_tmp);
    } else if (bim == brm) {
      if (r > 0.0) {
        f_r = 0.5;
      } else {
        f_r = -0.5;
      }

      if (t_s_im_tmp > 0.0) {
        f_t_s_im_tmp = 0.5;
      } else {
        f_t_s_im_tmp = -0.5;
      }

      y_im = (ar * f_r + 0.0 * f_t_s_im_tmp) / brm;
    } else {
      brm = r / t_s_im_tmp;
      y_im = brm * ar / (t_s_im_tmp + brm * r);
    }
  }

  ar = 12.0 * x04 * x12;
  if (t_s_im_tmp == 0.0) {
    b_y_re = ar / r;
  } else if (r == 0.0) {
    if (ar == 0.0) {
      b_y_re = 0.0 / t_s_im_tmp;
    } else {
      b_y_re = 0.0;
    }
  } else {
    brm = fabs(r);
    bim = fabs(t_s_im_tmp);
    if (brm > bim) {
      brm = t_s_im_tmp / r;
      b_y_re = (ar + brm * 0.0) / (r + brm * t_s_im_tmp);
    } else if (bim == brm) {
      if (r > 0.0) {
        g_r = 0.5;
      } else {
        g_r = -0.5;
      }

      if (t_s_im_tmp > 0.0) {
        g_t_s_im_tmp = 0.5;
      } else {
        g_t_s_im_tmp = -0.5;
      }

      b_y_re = (ar * g_r + 0.0 * g_t_s_im_tmp) / brm;
    } else {
      brm = r / t_s_im_tmp;
      b_y_re = brm * ar / (t_s_im_tmp + brm * r);
    }
  }

  ar = 4.0 * x04 * x14;
  if (t_s.im == 0.0) {
    b_y_im = ar / t_s.re;
  } else if (t_s.re == 0.0) {
    if (ar == 0.0) {
      b_y_im = 0.0 / t_s.im;
    } else {
      b_y_im = 0.0;
    }
  } else {
    brm = fabs(t_s.re);
    bim = fabs(t_s.im);
    if (brm > bim) {
      brm = t_s.im / t_s.re;
      b_y_im = (ar + brm * 0.0) / (t_s.re + brm * t_s.im);
    } else if (bim == brm) {
      if (t_s.re > 0.0) {
        l_t_s = 0.5;
      } else {
        l_t_s = -0.5;
      }

      if (t_s.im > 0.0) {
        m_t_s = 0.5;
      } else {
        m_t_s = -0.5;
      }

      b_y_im = (ar * l_t_s + 0.0 * m_t_s) / brm;
    } else {
      brm = t_s.re / t_s.im;
      b_y_im = brm * ar / (t_s.im + brm * t_s.re);
    }
  }

  ar = 12.0 * x11 * x13;
  if (t_s_im_tmp == 0.0) {
    c_y_re = ar / r;
  } else if (r == 0.0) {
    if (ar == 0.0) {
      c_y_re = 0.0 / t_s_im_tmp;
    } else {
      c_y_re = 0.0;
    }
  } else {
    brm = fabs(r);
    bim = fabs(t_s_im_tmp);
    if (brm > bim) {
      brm = t_s_im_tmp / r;
      c_y_re = (ar + brm * 0.0) / (r + brm * t_s_im_tmp);
    } else if (bim == brm) {
      if (r > 0.0) {
        h_r = 0.5;
      } else {
        h_r = -0.5;
      }

      if (t_s_im_tmp > 0.0) {
        h_t_s_im_tmp = 0.5;
      } else {
        h_t_s_im_tmp = -0.5;
      }

      c_y_re = (ar * h_r + 0.0 * h_t_s_im_tmp) / brm;
    } else {
      brm = r / t_s_im_tmp;
      c_y_re = brm * ar / (t_s_im_tmp + brm * r);
    }
  }

  ar = 12.0 * x12 * x14;
  if (t_s_im_tmp == 0.0) {
    r = ar / r;
  } else if (r == 0.0) {
    if (ar == 0.0) {
      r = 0.0 / t_s_im_tmp;
    } else {
      r = 0.0;
    }
  } else {
    brm = fabs(r);
    bim = fabs(t_s_im_tmp);
    if (brm > bim) {
      brm = t_s_im_tmp / r;
      r = (ar + brm * 0.0) / (r + brm * t_s_im_tmp);
    } else if (bim == brm) {
      if (r > 0.0) {
        i_r = 0.5;
      } else {
        i_r = -0.5;
      }

      if (t_s_im_tmp > 0.0) {
        i_t_s_im_tmp = 0.5;
      } else {
        i_t_s_im_tmp = -0.5;
      }

      r = (ar * i_r + 0.0 * i_t_s_im_tmp) / brm;
    } else {
      brm = r / t_s_im_tmp;
      r = brm * ar / (t_s_im_tmp + brm * r);
    }
  }

  return (((((((((((((((((((t_s.re + re) + b_re) + c_re) + d_re) + e_re) + f_re)
                      + g_re) + h_re) + i_re) + j_re) - c_y_im) + d_y_re) -
                d_y_im) + e_y_re) - y_re) + y_im) - b_y_re) + b_y_im) - c_y_re)
    - r;
}

/* End of code generation (DI_cost.c) */
