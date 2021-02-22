/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * DI_stateFreeVel.c
 *
 * Code generation for function 'DI_stateFreeVel'
 *
 */

/* Include files */
#include <math.h>
#include "rt_nonfinite.h"
#include "benchmarkRRT.h"
#include "DI_stateFreeVel.h"
#include "DI_cost.h"
#include "benchmarkRRT_rtwutil.h"

/* Function Definitions */
void DI_stateFreeVel(const creal_T t, const creal_T t_s, double x01, double x02,
                     double x03, double x04, double x11, double x12, double
                     states[4])
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
  double t_re_tmp;
  double t_im_tmp;
  double b_t_re_tmp;
  double b_t_im_tmp;
  double t_im;
  double x01_re_tmp;
  double x01_im_tmp;
  double t_re;
  double br;
  double bi;
  double bim;
  double b_br;
  double x02_re_tmp;
  double x02_im_tmp;
  double b_bi;
  double c_br;
  double c_bi;
  double d_br;
  double d_bi;
  double e_br;
  double e_bi;
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

  t_re_tmp = t.re * t.re - t.im * t.im;
  t_im_tmp = t.re * t.im + t.im * t.re;
  b_t_re_tmp = t.re - 3.0 * t_s.re;
  b_t_im_tmp = t.im - 3.0 * t_s.im;
  r = t_re_tmp * b_t_re_tmp - t_im_tmp * b_t_im_tmp;
  t_im = t_re_tmp * b_t_im_tmp + t_im_tmp * b_t_re_tmp;
  x01_re_tmp = (x01 - x11) + t_s.re * x03;
  x01_im_tmp = t_s.im * x03;
  t_re = r * x01_re_tmp - t_im * x01_im_tmp;
  t_im = r * x01_im_tmp + t_im * x01_re_tmp;
  br = 2.0 * y_re;
  bi = 2.0 * y_im;
  if (bi == 0.0) {
    if (t_im == 0.0) {
      r = t_re / br;
    } else if (t_re == 0.0) {
      r = 0.0;
    } else {
      r = t_re / br;
    }
  } else if (br == 0.0) {
    if (t_re == 0.0) {
      r = t_im / bi;
    } else if (t_im == 0.0) {
      r = 0.0;
    } else {
      r = t_im / bi;
    }
  } else {
    r = fabs(br);
    bim = fabs(bi);
    if (r > bim) {
      r = bi / br;
      r = (t_re + r * t_im) / (br + r * bi);
    } else if (bim == r) {
      if (br > 0.0) {
        b_br = 0.5;
      } else {
        b_br = -0.5;
      }

      if (bi > 0.0) {
        b_bi = 0.5;
      } else {
        b_bi = -0.5;
      }

      r = (t_re * b_br + t_im * b_bi) / r;
    } else {
      r = br / bi;
      r = (r * t_re + t_im) / (bi + r * br);
    }
  }

  states[0] = (x01 + t.re * x03) + r;
  r = t_re_tmp * b_t_re_tmp - t_im_tmp * b_t_im_tmp;
  t_im = t_re_tmp * b_t_im_tmp + t_im_tmp * b_t_re_tmp;
  x02_re_tmp = (x02 - x12) + t_s.re * x04;
  x02_im_tmp = t_s.im * x04;
  t_re = r * x02_re_tmp - t_im * x02_im_tmp;
  t_im = r * x02_im_tmp + t_im * x02_re_tmp;
  br = 2.0 * b_y_re;
  bi = 2.0 * b_y_im;
  if (bi == 0.0) {
    if (t_im == 0.0) {
      r = t_re / br;
    } else if (t_re == 0.0) {
      r = 0.0;
    } else {
      r = t_re / br;
    }
  } else if (br == 0.0) {
    if (t_re == 0.0) {
      r = t_im / bi;
    } else if (t_im == 0.0) {
      r = 0.0;
    } else {
      r = t_im / bi;
    }
  } else {
    r = fabs(br);
    bim = fabs(bi);
    if (r > bim) {
      r = bi / br;
      r = (t_re + r * t_im) / (br + r * bi);
    } else if (bim == r) {
      if (br > 0.0) {
        c_br = 0.5;
      } else {
        c_br = -0.5;
      }

      if (bi > 0.0) {
        c_bi = 0.5;
      } else {
        c_bi = -0.5;
      }

      r = (t_re * c_br + t_im * c_bi) / r;
    } else {
      r = br / bi;
      r = (r * t_re + t_im) / (bi + r * br);
    }
  }

  states[1] = (x02 + t.re * x04) + r;
  y_re = 3.0 * t.re;
  y_im = 3.0 * t.im;
  t_re_tmp = t.re - 2.0 * t_s.re;
  t_im_tmp = t.im - 2.0 * t_s.im;
  r = y_re * t_re_tmp - y_im * t_im_tmp;
  b_t_re_tmp = y_re * t_im_tmp + y_im * t_re_tmp;
  b_t_im_tmp = r * x01_re_tmp - b_t_re_tmp * x01_im_tmp;
  b_t_re_tmp = r * x01_im_tmp + b_t_re_tmp * x01_re_tmp;
  br = 2.0 * c_y_re;
  bi = 2.0 * c_y_im;
  if (bi == 0.0) {
    if (b_t_re_tmp == 0.0) {
      r = b_t_im_tmp / br;
    } else if (b_t_im_tmp == 0.0) {
      r = 0.0;
    } else {
      r = b_t_im_tmp / br;
    }
  } else if (br == 0.0) {
    if (b_t_im_tmp == 0.0) {
      r = b_t_re_tmp / bi;
    } else if (b_t_re_tmp == 0.0) {
      r = 0.0;
    } else {
      r = b_t_re_tmp / bi;
    }
  } else {
    r = fabs(br);
    bim = fabs(bi);
    if (r > bim) {
      r = bi / br;
      r = (b_t_im_tmp + r * b_t_re_tmp) / (br + r * bi);
    } else if (bim == r) {
      if (br > 0.0) {
        d_br = 0.5;
      } else {
        d_br = -0.5;
      }

      if (bi > 0.0) {
        d_bi = 0.5;
      } else {
        d_bi = -0.5;
      }

      r = (b_t_im_tmp * d_br + b_t_re_tmp * d_bi) / r;
    } else {
      r = br / bi;
      r = (r * b_t_im_tmp + b_t_re_tmp) / (bi + r * br);
    }
  }

  states[2] = x03 + r;
  r = y_re * t_re_tmp - y_im * t_im_tmp;
  b_t_re_tmp = y_re * t_im_tmp + y_im * t_re_tmp;
  b_t_im_tmp = r * x02_re_tmp - b_t_re_tmp * x02_im_tmp;
  b_t_re_tmp = r * x02_im_tmp + b_t_re_tmp * x02_re_tmp;
  br = 2.0 * x_re;
  bi = 2.0 * x_im;
  if (bi == 0.0) {
    if (b_t_re_tmp == 0.0) {
      r = b_t_im_tmp / br;
    } else if (b_t_im_tmp == 0.0) {
      r = 0.0;
    } else {
      r = b_t_im_tmp / br;
    }
  } else if (br == 0.0) {
    if (b_t_im_tmp == 0.0) {
      r = b_t_re_tmp / bi;
    } else if (b_t_re_tmp == 0.0) {
      r = 0.0;
    } else {
      r = b_t_re_tmp / bi;
    }
  } else {
    r = fabs(br);
    bim = fabs(bi);
    if (r > bim) {
      r = bi / br;
      r = (b_t_im_tmp + r * b_t_re_tmp) / (br + r * bi);
    } else if (bim == r) {
      if (br > 0.0) {
        e_br = 0.5;
      } else {
        e_br = -0.5;
      }

      if (bi > 0.0) {
        e_bi = 0.5;
      } else {
        e_bi = -0.5;
      }

      r = (b_t_im_tmp * e_br + b_t_re_tmp * e_bi) / r;
    } else {
      r = br / bi;
      r = (r * b_t_im_tmp + b_t_re_tmp) / (bi + r * br);
    }
  }

  states[3] = x04 + r;
}

/* End of code generation (DI_stateFreeVel.c) */
