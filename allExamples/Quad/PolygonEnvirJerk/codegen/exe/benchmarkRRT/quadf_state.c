/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * quadf_state.c
 *
 * Code generation for function 'quadf_state'
 *
 */

/* Include files */
#include <math.h>
#include "rt_nonfinite.h"
#include "benchmarkRRT.h"
#include "quadf_state.h"
#include "quadf_cost.h"
#include "benchmarkRRT_rtwutil.h"

/* Function Definitions */
void quadf_state(const creal_T t, const creal_T t_s, double x01, double x02,
                 double x03, double x04, double x05, double x06, double x07,
                 double x08, double x09, double x11, double x12, double x13,
                 double x14, double x15, double x16, double x17, double x18,
                 double x19, double states[3])
{
  double a_re_tmp;
  double a_im_tmp;
  double y_re;
  double y_im;
  double r;
  double a_im;
  double b_y_re;
  double b_y_im;
  double x_re;
  double x_im;
  double c_y_re;
  double c_y_im;
  double d_y_re;
  double d_y_im;
  double e_y_re;
  double e_y_im;
  double f_y_re;
  double f_y_im;
  double g_y_re;
  double g_y_im;
  double h_y_re;
  double h_y_im;
  double i_y_re;
  double i_y_im;
  double j_y_re;
  double j_y_im;
  double k_y_re;
  double k_y_im;
  double l_y_re;
  double l_y_im;
  double m_y_re;
  double m_y_im;
  double n_y_re;
  double n_y_im;
  double o_y_re;
  double o_y_im;
  double p_y_re;
  double p_y_im;
  double q_y_re;
  double q_y_im;
  double re_tmp;
  double im_tmp;
  double t_s_re_tmp;
  double t_s_im_tmp;
  double re;
  double im;
  double br;
  double b_re_tmp;
  double b_im_tmp;
  double b_br;
  double c_re_tmp;
  double c_im_tmp;
  double b_im;
  double c_br;
  double d_im_tmp;
  double c_im;
  double d_br;
  double d_im;
  double e_br;
  double e_im;
  double f_br;
  double f_im;
  double g_br;
  double g_im;
  double h_br;
  double h_im;
  double i_br;
  double i_im;
  double j_br;
  double j_im;
  a_re_tmp = t.re - t_s.re;
  a_im_tmp = t.im - t_s.im;
  if ((a_im_tmp == 0.0) && (a_re_tmp >= 0.0)) {
    y_re = rt_powd_snf(a_re_tmp, 5.0);
    y_im = 0.0;
  } else if (a_re_tmp == 0.0) {
    y_re = 0.0;
    y_im = rt_powd_snf(a_im_tmp, 5.0);
  } else {
    if (a_im_tmp == 0.0) {
      if (a_re_tmp < 0.0) {
        r = log(fabs(a_re_tmp));
        a_im = 3.1415926535897931;
      } else {
        r = log(a_re_tmp);
        a_im = 0.0;
      }
    } else if ((fabs(a_re_tmp) > 8.9884656743115785E+307) || (fabs(a_im_tmp) >
                8.9884656743115785E+307)) {
      r = log(rt_hypotd_snf(a_re_tmp / 2.0, a_im_tmp / 2.0)) +
        0.69314718055994529;
      a_im = rt_atan2d_snf(a_im_tmp, a_re_tmp);
    } else {
      r = log(rt_hypotd_snf(a_re_tmp, a_im_tmp));
      a_im = rt_atan2d_snf(a_im_tmp, a_re_tmp);
    }

    y_re = 5.0 * r;
    y_im = 5.0 * a_im;
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
    b_y_re = rt_powd_snf(t_s.re, 5.0);
    b_y_im = 0.0;
  } else if (t_s.re == 0.0) {
    b_y_re = 0.0;
    b_y_im = rt_powd_snf(t_s.im, 5.0);
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

    b_y_re = 5.0 * x_re;
    b_y_im = 5.0 * x_im;
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
        r = log(fabs(a_re_tmp));
        a_im = 3.1415926535897931;
      } else {
        r = log(a_re_tmp);
        a_im = 0.0;
      }
    } else if ((fabs(a_re_tmp) > 8.9884656743115785E+307) || (fabs(a_im_tmp) >
                8.9884656743115785E+307)) {
      r = log(rt_hypotd_snf(a_re_tmp / 2.0, a_im_tmp / 2.0)) +
        0.69314718055994529;
      a_im = rt_atan2d_snf(a_im_tmp, a_re_tmp);
    } else {
      r = log(rt_hypotd_snf(a_re_tmp, a_im_tmp));
      a_im = rt_atan2d_snf(a_im_tmp, a_re_tmp);
    }

    c_y_re = 3.0 * r;
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

  if ((a_im_tmp == 0.0) && (a_re_tmp >= 0.0)) {
    e_y_re = rt_powd_snf(a_re_tmp, 4.0);
    e_y_im = 0.0;
  } else if (a_re_tmp == 0.0) {
    e_y_re = rt_powd_snf(a_im_tmp, 4.0);
    e_y_im = 0.0;
  } else {
    if (a_im_tmp == 0.0) {
      if (a_re_tmp < 0.0) {
        r = log(fabs(a_re_tmp));
        a_im = 3.1415926535897931;
      } else {
        r = log(a_re_tmp);
        a_im = 0.0;
      }
    } else if ((fabs(a_re_tmp) > 8.9884656743115785E+307) || (fabs(a_im_tmp) >
                8.9884656743115785E+307)) {
      r = log(rt_hypotd_snf(a_re_tmp / 2.0, a_im_tmp / 2.0)) +
        0.69314718055994529;
      a_im = rt_atan2d_snf(a_im_tmp, a_re_tmp);
    } else {
      r = log(rt_hypotd_snf(a_re_tmp, a_im_tmp));
      a_im = rt_atan2d_snf(a_im_tmp, a_re_tmp);
    }

    e_y_re = 4.0 * r;
    e_y_im = 4.0 * a_im;
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
    f_y_re = rt_powd_snf(t_s.re, 4.0);
    f_y_im = 0.0;
  } else if (t_s.re == 0.0) {
    f_y_re = rt_powd_snf(t_s.im, 4.0);
    f_y_im = 0.0;
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

    f_y_re = 4.0 * x_re;
    f_y_im = 4.0 * x_im;
    if (f_y_im == 0.0) {
      f_y_re = exp(f_y_re);
      f_y_im = 0.0;
    } else if (rtIsInf(f_y_im) && rtIsInf(f_y_re) && (f_y_re < 0.0)) {
      f_y_re = 0.0;
      f_y_im = 0.0;
    } else {
      r = exp(f_y_re / 2.0);
      f_y_re = r * (r * cos(f_y_im));
      f_y_im = r * (r * sin(f_y_im));
    }
  }

  if ((a_im_tmp == 0.0) && (a_re_tmp >= 0.0)) {
    g_y_re = rt_powd_snf(a_re_tmp, 5.0);
    g_y_im = 0.0;
  } else if (a_re_tmp == 0.0) {
    g_y_re = 0.0;
    g_y_im = rt_powd_snf(a_im_tmp, 5.0);
  } else {
    if (a_im_tmp == 0.0) {
      if (a_re_tmp < 0.0) {
        r = log(fabs(a_re_tmp));
        a_im = 3.1415926535897931;
      } else {
        r = log(a_re_tmp);
        a_im = 0.0;
      }
    } else if ((fabs(a_re_tmp) > 8.9884656743115785E+307) || (fabs(a_im_tmp) >
                8.9884656743115785E+307)) {
      r = log(rt_hypotd_snf(a_re_tmp / 2.0, a_im_tmp / 2.0)) +
        0.69314718055994529;
      a_im = rt_atan2d_snf(a_im_tmp, a_re_tmp);
    } else {
      r = log(rt_hypotd_snf(a_re_tmp, a_im_tmp));
      a_im = rt_atan2d_snf(a_im_tmp, a_re_tmp);
    }

    g_y_re = 5.0 * r;
    g_y_im = 5.0 * a_im;
    if (g_y_im == 0.0) {
      g_y_re = exp(g_y_re);
      g_y_im = 0.0;
    } else if (rtIsInf(g_y_im) && rtIsInf(g_y_re) && (g_y_re < 0.0)) {
      g_y_re = 0.0;
      g_y_im = 0.0;
    } else {
      r = exp(g_y_re / 2.0);
      g_y_re = r * (r * cos(g_y_im));
      g_y_im = r * (r * sin(g_y_im));
    }
  }

  if ((t_s.im == 0.0) && (t_s.re >= 0.0)) {
    h_y_re = rt_powd_snf(t_s.re, 5.0);
    h_y_im = 0.0;
  } else if (t_s.re == 0.0) {
    h_y_re = 0.0;
    h_y_im = rt_powd_snf(t_s.im, 5.0);
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

    h_y_re = 5.0 * x_re;
    h_y_im = 5.0 * x_im;
    if (h_y_im == 0.0) {
      h_y_re = exp(h_y_re);
      h_y_im = 0.0;
    } else if (rtIsInf(h_y_im) && rtIsInf(h_y_re) && (h_y_re < 0.0)) {
      h_y_re = 0.0;
      h_y_im = 0.0;
    } else {
      r = exp(h_y_re / 2.0);
      h_y_re = r * (r * cos(h_y_im));
      h_y_im = r * (r * sin(h_y_im));
    }
  }

  if ((a_im_tmp == 0.0) && (a_re_tmp >= 0.0)) {
    i_y_re = rt_powd_snf(a_re_tmp, 3.0);
    i_y_im = 0.0;
  } else if (a_re_tmp == 0.0) {
    i_y_re = 0.0;
    i_y_im = -rt_powd_snf(a_im_tmp, 3.0);
  } else {
    if (a_im_tmp == 0.0) {
      if (a_re_tmp < 0.0) {
        r = log(fabs(a_re_tmp));
        a_im = 3.1415926535897931;
      } else {
        r = log(a_re_tmp);
        a_im = 0.0;
      }
    } else if ((fabs(a_re_tmp) > 8.9884656743115785E+307) || (fabs(a_im_tmp) >
                8.9884656743115785E+307)) {
      r = log(rt_hypotd_snf(a_re_tmp / 2.0, a_im_tmp / 2.0)) +
        0.69314718055994529;
      a_im = rt_atan2d_snf(a_im_tmp, a_re_tmp);
    } else {
      r = log(rt_hypotd_snf(a_re_tmp, a_im_tmp));
      a_im = rt_atan2d_snf(a_im_tmp, a_re_tmp);
    }

    i_y_re = 3.0 * r;
    i_y_im = 3.0 * a_im;
    if (i_y_im == 0.0) {
      i_y_re = exp(i_y_re);
      i_y_im = 0.0;
    } else if (rtIsInf(i_y_im) && rtIsInf(i_y_re) && (i_y_re < 0.0)) {
      i_y_re = 0.0;
      i_y_im = 0.0;
    } else {
      r = exp(i_y_re / 2.0);
      i_y_re = r * (r * cos(i_y_im));
      i_y_im = r * (r * sin(i_y_im));
    }
  }

  if ((t_s.im == 0.0) && (t_s.re >= 0.0)) {
    j_y_re = rt_powd_snf(t_s.re, 3.0);
    j_y_im = 0.0;
  } else if (t_s.re == 0.0) {
    j_y_re = 0.0;
    j_y_im = -rt_powd_snf(t_s.im, 3.0);
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

    j_y_re = 3.0 * x_re;
    j_y_im = 3.0 * x_im;
    if (j_y_im == 0.0) {
      j_y_re = exp(j_y_re);
      j_y_im = 0.0;
    } else if (rtIsInf(j_y_im) && rtIsInf(j_y_re) && (j_y_re < 0.0)) {
      j_y_re = 0.0;
      j_y_im = 0.0;
    } else {
      r = exp(j_y_re / 2.0);
      j_y_re = r * (r * cos(j_y_im));
      j_y_im = r * (r * sin(j_y_im));
    }
  }

  if ((a_im_tmp == 0.0) && (a_re_tmp >= 0.0)) {
    k_y_re = rt_powd_snf(a_re_tmp, 4.0);
    k_y_im = 0.0;
  } else if (a_re_tmp == 0.0) {
    k_y_re = rt_powd_snf(a_im_tmp, 4.0);
    k_y_im = 0.0;
  } else {
    if (a_im_tmp == 0.0) {
      if (a_re_tmp < 0.0) {
        r = log(fabs(a_re_tmp));
        a_im = 3.1415926535897931;
      } else {
        r = log(a_re_tmp);
        a_im = 0.0;
      }
    } else if ((fabs(a_re_tmp) > 8.9884656743115785E+307) || (fabs(a_im_tmp) >
                8.9884656743115785E+307)) {
      r = log(rt_hypotd_snf(a_re_tmp / 2.0, a_im_tmp / 2.0)) +
        0.69314718055994529;
      a_im = rt_atan2d_snf(a_im_tmp, a_re_tmp);
    } else {
      r = log(rt_hypotd_snf(a_re_tmp, a_im_tmp));
      a_im = rt_atan2d_snf(a_im_tmp, a_re_tmp);
    }

    k_y_re = 4.0 * r;
    k_y_im = 4.0 * a_im;
    if (k_y_im == 0.0) {
      k_y_re = exp(k_y_re);
      k_y_im = 0.0;
    } else if (rtIsInf(k_y_im) && rtIsInf(k_y_re) && (k_y_re < 0.0)) {
      k_y_re = 0.0;
      k_y_im = 0.0;
    } else {
      r = exp(k_y_re / 2.0);
      k_y_re = r * (r * cos(k_y_im));
      k_y_im = r * (r * sin(k_y_im));
    }
  }

  if ((t_s.im == 0.0) && (t_s.re >= 0.0)) {
    l_y_re = rt_powd_snf(t_s.re, 4.0);
    l_y_im = 0.0;
  } else if (t_s.re == 0.0) {
    l_y_re = rt_powd_snf(t_s.im, 4.0);
    l_y_im = 0.0;
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

    l_y_re = 4.0 * x_re;
    l_y_im = 4.0 * x_im;
    if (l_y_im == 0.0) {
      l_y_re = exp(l_y_re);
      l_y_im = 0.0;
    } else if (rtIsInf(l_y_im) && rtIsInf(l_y_re) && (l_y_re < 0.0)) {
      l_y_re = 0.0;
      l_y_im = 0.0;
    } else {
      r = exp(l_y_re / 2.0);
      l_y_re = r * (r * cos(l_y_im));
      l_y_im = r * (r * sin(l_y_im));
    }
  }

  if ((a_im_tmp == 0.0) && (a_re_tmp >= 0.0)) {
    m_y_re = rt_powd_snf(a_re_tmp, 5.0);
    m_y_im = 0.0;
  } else if (a_re_tmp == 0.0) {
    m_y_re = 0.0;
    m_y_im = rt_powd_snf(a_im_tmp, 5.0);
  } else {
    if (a_im_tmp == 0.0) {
      if (a_re_tmp < 0.0) {
        x_re = log(fabs(a_re_tmp));
        x_im = 3.1415926535897931;
      } else {
        x_re = log(a_re_tmp);
        x_im = 0.0;
      }
    } else if ((fabs(a_re_tmp) > 8.9884656743115785E+307) || (fabs(a_im_tmp) >
                8.9884656743115785E+307)) {
      x_re = log(rt_hypotd_snf(a_re_tmp / 2.0, a_im_tmp / 2.0)) +
        0.69314718055994529;
      x_im = rt_atan2d_snf(a_im_tmp, a_re_tmp);
    } else {
      x_re = log(rt_hypotd_snf(a_re_tmp, a_im_tmp));
      x_im = rt_atan2d_snf(a_im_tmp, a_re_tmp);
    }

    m_y_re = 5.0 * x_re;
    m_y_im = 5.0 * x_im;
    if (m_y_im == 0.0) {
      m_y_re = exp(m_y_re);
      m_y_im = 0.0;
    } else if (rtIsInf(m_y_im) && rtIsInf(m_y_re) && (m_y_re < 0.0)) {
      m_y_re = 0.0;
      m_y_im = 0.0;
    } else {
      r = exp(m_y_re / 2.0);
      m_y_re = r * (r * cos(m_y_im));
      m_y_im = r * (r * sin(m_y_im));
    }
  }

  if ((t_s.im == 0.0) && (t_s.re >= 0.0)) {
    n_y_re = rt_powd_snf(t_s.re, 5.0);
    n_y_im = 0.0;
  } else if (t_s.re == 0.0) {
    n_y_re = 0.0;
    n_y_im = rt_powd_snf(t_s.im, 5.0);
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

    n_y_re = 5.0 * x_re;
    n_y_im = 5.0 * x_im;
    if (n_y_im == 0.0) {
      n_y_re = exp(n_y_re);
      n_y_im = 0.0;
    } else if (rtIsInf(n_y_im) && rtIsInf(n_y_re) && (n_y_re < 0.0)) {
      n_y_re = 0.0;
      n_y_im = 0.0;
    } else {
      r = exp(n_y_re / 2.0);
      n_y_re = r * (r * cos(n_y_im));
      n_y_im = r * (r * sin(n_y_im));
    }
  }

  if ((a_im_tmp == 0.0) && (a_re_tmp >= 0.0)) {
    o_y_re = rt_powd_snf(a_re_tmp, 3.0);
    o_y_im = 0.0;
  } else if (a_re_tmp == 0.0) {
    o_y_re = 0.0;
    o_y_im = -rt_powd_snf(a_im_tmp, 3.0);
  } else {
    if (a_im_tmp == 0.0) {
      if (a_re_tmp < 0.0) {
        x_re = log(fabs(a_re_tmp));
        x_im = 3.1415926535897931;
      } else {
        x_re = log(a_re_tmp);
        x_im = 0.0;
      }
    } else if ((fabs(a_re_tmp) > 8.9884656743115785E+307) || (fabs(a_im_tmp) >
                8.9884656743115785E+307)) {
      x_re = log(rt_hypotd_snf(a_re_tmp / 2.0, a_im_tmp / 2.0)) +
        0.69314718055994529;
      x_im = rt_atan2d_snf(a_im_tmp, a_re_tmp);
    } else {
      x_re = log(rt_hypotd_snf(a_re_tmp, a_im_tmp));
      x_im = rt_atan2d_snf(a_im_tmp, a_re_tmp);
    }

    o_y_re = 3.0 * x_re;
    o_y_im = 3.0 * x_im;
    if (o_y_im == 0.0) {
      o_y_re = exp(o_y_re);
      o_y_im = 0.0;
    } else if (rtIsInf(o_y_im) && rtIsInf(o_y_re) && (o_y_re < 0.0)) {
      o_y_re = 0.0;
      o_y_im = 0.0;
    } else {
      r = exp(o_y_re / 2.0);
      o_y_re = r * (r * cos(o_y_im));
      o_y_im = r * (r * sin(o_y_im));
    }
  }

  if ((t_s.im == 0.0) && (t_s.re >= 0.0)) {
    p_y_re = rt_powd_snf(t_s.re, 3.0);
    p_y_im = 0.0;
  } else if (t_s.re == 0.0) {
    p_y_re = 0.0;
    p_y_im = -rt_powd_snf(t_s.im, 3.0);
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

    p_y_re = 3.0 * x_re;
    p_y_im = 3.0 * x_im;
    if (p_y_im == 0.0) {
      p_y_re = exp(p_y_re);
      p_y_im = 0.0;
    } else if (rtIsInf(p_y_im) && rtIsInf(p_y_re) && (p_y_re < 0.0)) {
      p_y_re = 0.0;
      p_y_im = 0.0;
    } else {
      r = exp(p_y_re / 2.0);
      p_y_re = r * (r * cos(p_y_im));
      p_y_im = r * (r * sin(p_y_im));
    }
  }

  if ((a_im_tmp == 0.0) && (a_re_tmp >= 0.0)) {
    q_y_re = rt_powd_snf(a_re_tmp, 4.0);
    q_y_im = 0.0;
  } else if (a_re_tmp == 0.0) {
    q_y_re = rt_powd_snf(a_im_tmp, 4.0);
    q_y_im = 0.0;
  } else {
    if (a_im_tmp == 0.0) {
      if (a_re_tmp < 0.0) {
        x_re = log(fabs(a_re_tmp));
        x_im = 3.1415926535897931;
      } else {
        x_re = log(a_re_tmp);
        x_im = 0.0;
      }
    } else if ((fabs(a_re_tmp) > 8.9884656743115785E+307) || (fabs(a_im_tmp) >
                8.9884656743115785E+307)) {
      x_re = log(rt_hypotd_snf(a_re_tmp / 2.0, a_im_tmp / 2.0)) +
        0.69314718055994529;
      x_im = rt_atan2d_snf(a_im_tmp, a_re_tmp);
    } else {
      x_re = log(rt_hypotd_snf(a_re_tmp, a_im_tmp));
      x_im = rt_atan2d_snf(a_im_tmp, a_re_tmp);
    }

    q_y_re = 4.0 * x_re;
    q_y_im = 4.0 * x_im;
    if (q_y_im == 0.0) {
      q_y_re = exp(q_y_re);
      q_y_im = 0.0;
    } else if (rtIsInf(q_y_im) && rtIsInf(q_y_re) && (q_y_re < 0.0)) {
      q_y_re = 0.0;
      q_y_im = 0.0;
    } else {
      r = exp(q_y_re / 2.0);
      q_y_re = r * (r * cos(q_y_im));
      q_y_im = r * (r * sin(q_y_im));
    }
  }

  if ((t_s.im == 0.0) && (t_s.re >= 0.0)) {
    x_re = rt_powd_snf(t_s.re, 4.0);
    x_im = 0.0;
  } else if (t_s.re == 0.0) {
    x_re = rt_powd_snf(t_s.im, 4.0);
    x_im = 0.0;
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

    x_re *= 4.0;
    x_im *= 4.0;
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

  /*                       x14 + x17*(t - t_s) - (5*(t - t_s)^4*(12*x01 - 12*x11 + 6*t_s*x04 + 6*t_s*x14 + t_s^2*x07 - t_s^2*x17))/(2*t_s^5) - (3*(t - t_s)^2*(20*x01 - 20*x11 + 8*t_s*x04 + 12*t_s*x14 + t_s^2*x07 - 3*t_s^2*x17))/(2*t_s^3) - (2*(t - t_s)^3*(30*x01 - 30*x11 + 14*t_s*x04 + 16*t_s*x14 + 2*t_s^2*x07 - 3*t_s^2*x17))/t_s^4 */
  /*                       x15 + x18*(t - t_s) - (5*(t - t_s)^4*(12*x02 - 12*x12 + 6*t_s*x05 + 6*t_s*x15 + t_s^2*x08 - t_s^2*x18))/(2*t_s^5) - (3*(t - t_s)^2*(20*x02 - 20*x12 + 8*t_s*x05 + 12*t_s*x15 + t_s^2*x08 - 3*t_s^2*x18))/(2*t_s^3) - (2*(t - t_s)^3*(30*x02 - 30*x12 + 14*t_s*x05 + 16*t_s*x15 + 2*t_s^2*x08 - 3*t_s^2*x18))/t_s^4 */
  /*                       x16 + x19*(t - t_s) - (5*(t - t_s)^4*(12*x03 - 12*x13 + 6*t_s*x06 + 6*t_s*x16 + t_s^2*x09 - t_s^2*x19))/(2*t_s^5) - (3*(t - t_s)^2*(20*x03 - 20*x13 + 8*t_s*x06 + 12*t_s*x16 + t_s^2*x09 - 3*t_s^2*x19))/(2*t_s^3) - (2*(t - t_s)^3*(30*x03 - 30*x13 + 14*t_s*x06 + 16*t_s*x16 + 2*t_s^2*x09 - 3*t_s^2*x19))/t_s^4 */
  /*                               x17 - (3*(5000*t - 5000*t_s)*(20*x01 - 20*x11 + 8*t_s*x04 + 12*t_s*x14 + t_s^2*x07 - 3*t_s^2*x17))/(5000*t_s^3) - (10*(t - t_s)^3*(12*x01 - 12*x11 + 6*t_s*x04 + 6*t_s*x14 + t_s^2*x07 - t_s^2*x17))/t_s^5 - (6*(t - t_s)^2*(30*x01 - 30*x11 + 14*t_s*x04 + 16*t_s*x14 + 2*t_s^2*x07 - 3*t_s^2*x17))/t_s^4 */
  /*                               x18 - (3*(5000*t - 5000*t_s)*(20*x02 - 20*x12 + 8*t_s*x05 + 12*t_s*x15 + t_s^2*x08 - 3*t_s^2*x18))/(5000*t_s^3) - (10*(t - t_s)^3*(12*x02 - 12*x12 + 6*t_s*x05 + 6*t_s*x15 + t_s^2*x08 - t_s^2*x18))/t_s^5 - (6*(t - t_s)^2*(30*x02 - 30*x12 + 14*t_s*x05 + 16*t_s*x15 + 2*t_s^2*x08 - 3*t_s^2*x18))/t_s^4 */
  /*                               x19 - (3*(5000*t - 5000*t_s)*(20*x03 - 20*x13 + 8*t_s*x06 + 12*t_s*x16 + t_s^2*x09 - 3*t_s^2*x19))/(5000*t_s^3) - (10*(t - t_s)^3*(12*x03 - 12*x13 + 6*t_s*x06 + 6*t_s*x16 + t_s^2*x09 - t_s^2*x19))/t_s^5 - (6*(t - t_s)^2*(30*x03 - 30*x13 + 14*t_s*x06 + 16*t_s*x16 + 2*t_s^2*x09 - 3*t_s^2*x19))/t_s^4]; */
  a_im = x17 * (a_re_tmp * a_re_tmp - a_im_tmp * a_im_tmp);
  if (x17 * (a_re_tmp * a_im_tmp + a_im_tmp * a_re_tmp) == 0.0) {
    r = a_im / 2.0;
  } else if (a_im == 0.0) {
    r = 0.0;
  } else {
    r = a_im / 2.0;
  }

  re_tmp = 6.0 * t_s.re;
  im_tmp = 6.0 * t_s.im;
  t_s_re_tmp = t_s.re * t_s.re - t_s.im * t_s.im;
  t_s_im_tmp = t_s.re * t_s.im + t_s.im * t_s.re;
  re = ((((12.0 * x01 - 12.0 * x11) + re_tmp * x04) + re_tmp * x14) + t_s_re_tmp
        * x07) - t_s_re_tmp * x17;
  im = ((im_tmp * x04 + im_tmp * x14) + t_s_im_tmp * x07) - t_s_im_tmp * x17;
  a_im = y_re * re - y_im * im;
  y_im = y_re * im + y_im * re;
  br = 2.0 * b_y_re;
  im = 2.0 * b_y_im;
  if (im == 0.0) {
    if (y_im == 0.0) {
      y_re = a_im / br;
    } else if (a_im == 0.0) {
      y_re = 0.0;
    } else {
      y_re = a_im / br;
    }
  } else if (br == 0.0) {
    if (a_im == 0.0) {
      y_re = y_im / im;
    } else if (y_im == 0.0) {
      y_re = 0.0;
    } else {
      y_re = y_im / im;
    }
  } else {
    b_y_im = fabs(br);
    re = fabs(im);
    if (b_y_im > re) {
      b_y_im = im / br;
      y_re = (a_im + b_y_im * y_im) / (br + b_y_im * im);
    } else if (re == b_y_im) {
      if (br > 0.0) {
        b_br = 0.5;
      } else {
        b_br = -0.5;
      }

      if (im > 0.0) {
        b_im = 0.5;
      } else {
        b_im = -0.5;
      }

      y_re = (a_im * b_br + y_im * b_im) / b_y_im;
    } else {
      b_y_im = br / im;
      y_re = (b_y_im * a_im + y_im) / (im + b_y_im * br);
    }
  }

  b_re_tmp = 8.0 * t_s.re;
  b_im_tmp = 8.0 * t_s.im;
  c_re_tmp = 12.0 * t_s.re;
  c_im_tmp = 12.0 * t_s.im;
  re = ((((20.0 * x01 - 20.0 * x11) + b_re_tmp * x04) + c_re_tmp * x14) +
        t_s_re_tmp * x07) - 3.0 * t_s_re_tmp * x17;
  im = ((b_im_tmp * x04 + c_im_tmp * x14) + t_s_im_tmp * x07) - 3.0 * t_s_im_tmp
    * x17;
  b_y_re = c_y_re * re - c_y_im * im;
  y_im = c_y_re * im + c_y_im * re;
  br = 2.0 * d_y_re;
  im = 2.0 * d_y_im;
  if (im == 0.0) {
    if (y_im == 0.0) {
      b_y_re /= br;
    } else if (b_y_re == 0.0) {
      b_y_re = 0.0;
    } else {
      b_y_re /= br;
    }
  } else if (br == 0.0) {
    if (b_y_re == 0.0) {
      b_y_re = y_im / im;
    } else if (y_im == 0.0) {
      b_y_re = 0.0;
    } else {
      b_y_re = y_im / im;
    }
  } else {
    b_y_im = fabs(br);
    re = fabs(im);
    if (b_y_im > re) {
      b_y_im = im / br;
      b_y_re = (b_y_re + b_y_im * y_im) / (br + b_y_im * im);
    } else if (re == b_y_im) {
      if (br > 0.0) {
        c_br = 0.5;
      } else {
        c_br = -0.5;
      }

      if (im > 0.0) {
        c_im = 0.5;
      } else {
        c_im = -0.5;
      }

      b_y_re = (b_y_re * c_br + y_im * c_im) / b_y_im;
    } else {
      b_y_im = br / im;
      b_y_re = (b_y_im * b_y_re + y_im) / (im + b_y_im * br);
    }
  }

  c_y_im = 14.0 * t_s.re;
  d_y_re = 14.0 * t_s.im;
  d_y_im = 16.0 * t_s.re;
  d_im_tmp = 16.0 * t_s.im;
  re = ((((30.0 * x01 - 30.0 * x11) + c_y_im * x04) + d_y_im * x14) + 2.0 *
        t_s_re_tmp * x07) - 3.0 * t_s_re_tmp * x17;
  im = ((d_y_re * x04 + d_im_tmp * x14) + 2.0 * t_s_im_tmp * x07) - 3.0 *
    t_s_im_tmp * x17;
  c_y_re = e_y_re * re - e_y_im * im;
  y_im = e_y_re * im + e_y_im * re;
  br = 2.0 * f_y_re;
  im = 2.0 * f_y_im;
  if (im == 0.0) {
    if (y_im == 0.0) {
      c_y_re /= br;
    } else if (c_y_re == 0.0) {
      c_y_re = 0.0;
    } else {
      c_y_re /= br;
    }
  } else if (br == 0.0) {
    if (c_y_re == 0.0) {
      c_y_re = y_im / im;
    } else if (y_im == 0.0) {
      c_y_re = 0.0;
    } else {
      c_y_re = y_im / im;
    }
  } else {
    b_y_im = fabs(br);
    re = fabs(im);
    if (b_y_im > re) {
      b_y_im = im / br;
      c_y_re = (c_y_re + b_y_im * y_im) / (br + b_y_im * im);
    } else if (re == b_y_im) {
      if (br > 0.0) {
        d_br = 0.5;
      } else {
        d_br = -0.5;
      }

      if (im > 0.0) {
        d_im = 0.5;
      } else {
        d_im = -0.5;
      }

      c_y_re = (c_y_re * d_br + y_im * d_im) / b_y_im;
    } else {
      b_y_im = br / im;
      c_y_re = (b_y_im * c_y_re + y_im) / (im + b_y_im * br);
    }
  }

  states[0] = ((((x11 + x14 * a_re_tmp) + r) - y_re) - b_y_re) - c_y_re;
  a_im = x18 * (a_re_tmp * a_re_tmp - a_im_tmp * a_im_tmp);
  if (x18 * (a_re_tmp * a_im_tmp + a_im_tmp * a_re_tmp) == 0.0) {
    a_im /= 2.0;
  } else if (a_im == 0.0) {
    a_im = 0.0;
  } else {
    a_im /= 2.0;
  }

  re = ((((12.0 * x02 - 12.0 * x12) + re_tmp * x05) + re_tmp * x15) + t_s_re_tmp
        * x08) - t_s_re_tmp * x18;
  im = ((im_tmp * x05 + im_tmp * x15) + t_s_im_tmp * x08) - t_s_im_tmp * x18;
  y_re = g_y_re * re - g_y_im * im;
  y_im = g_y_re * im + g_y_im * re;
  br = 2.0 * h_y_re;
  im = 2.0 * h_y_im;
  if (im == 0.0) {
    if (y_im == 0.0) {
      y_re /= br;
    } else if (y_re == 0.0) {
      y_re = 0.0;
    } else {
      y_re /= br;
    }
  } else if (br == 0.0) {
    if (y_re == 0.0) {
      y_re = y_im / im;
    } else if (y_im == 0.0) {
      y_re = 0.0;
    } else {
      y_re = y_im / im;
    }
  } else {
    b_y_im = fabs(br);
    re = fabs(im);
    if (b_y_im > re) {
      b_y_im = im / br;
      y_re = (y_re + b_y_im * y_im) / (br + b_y_im * im);
    } else if (re == b_y_im) {
      if (br > 0.0) {
        e_br = 0.5;
      } else {
        e_br = -0.5;
      }

      if (im > 0.0) {
        e_im = 0.5;
      } else {
        e_im = -0.5;
      }

      y_re = (y_re * e_br + y_im * e_im) / b_y_im;
    } else {
      b_y_im = br / im;
      y_re = (b_y_im * y_re + y_im) / (im + b_y_im * br);
    }
  }

  re = ((((20.0 * x02 - 20.0 * x12) + b_re_tmp * x05) + c_re_tmp * x15) +
        t_s_re_tmp * x08) - 3.0 * t_s_re_tmp * x18;
  im = ((b_im_tmp * x05 + c_im_tmp * x15) + t_s_im_tmp * x08) - 3.0 * t_s_im_tmp
    * x18;
  b_y_re = i_y_re * re - i_y_im * im;
  y_im = i_y_re * im + i_y_im * re;
  br = 2.0 * j_y_re;
  im = 2.0 * j_y_im;
  if (im == 0.0) {
    if (y_im == 0.0) {
      b_y_re /= br;
    } else if (b_y_re == 0.0) {
      b_y_re = 0.0;
    } else {
      b_y_re /= br;
    }
  } else if (br == 0.0) {
    if (b_y_re == 0.0) {
      b_y_re = y_im / im;
    } else if (y_im == 0.0) {
      b_y_re = 0.0;
    } else {
      b_y_re = y_im / im;
    }
  } else {
    b_y_im = fabs(br);
    re = fabs(im);
    if (b_y_im > re) {
      b_y_im = im / br;
      b_y_re = (b_y_re + b_y_im * y_im) / (br + b_y_im * im);
    } else if (re == b_y_im) {
      if (br > 0.0) {
        f_br = 0.5;
      } else {
        f_br = -0.5;
      }

      if (im > 0.0) {
        f_im = 0.5;
      } else {
        f_im = -0.5;
      }

      b_y_re = (b_y_re * f_br + y_im * f_im) / b_y_im;
    } else {
      b_y_im = br / im;
      b_y_re = (b_y_im * b_y_re + y_im) / (im + b_y_im * br);
    }
  }

  re = ((((30.0 * x02 - 30.0 * x12) + c_y_im * x05) + d_y_im * x15) + 2.0 *
        t_s_re_tmp * x08) - 3.0 * t_s_re_tmp * x18;
  im = ((d_y_re * x05 + d_im_tmp * x15) + 2.0 * t_s_im_tmp * x08) - 3.0 *
    t_s_im_tmp * x18;
  c_y_re = k_y_re * re - k_y_im * im;
  y_im = k_y_re * im + k_y_im * re;
  br = 2.0 * l_y_re;
  im = 2.0 * l_y_im;
  if (im == 0.0) {
    if (y_im == 0.0) {
      c_y_re /= br;
    } else if (c_y_re == 0.0) {
      c_y_re = 0.0;
    } else {
      c_y_re /= br;
    }
  } else if (br == 0.0) {
    if (c_y_re == 0.0) {
      c_y_re = y_im / im;
    } else if (y_im == 0.0) {
      c_y_re = 0.0;
    } else {
      c_y_re = y_im / im;
    }
  } else {
    b_y_im = fabs(br);
    re = fabs(im);
    if (b_y_im > re) {
      b_y_im = im / br;
      c_y_re = (c_y_re + b_y_im * y_im) / (br + b_y_im * im);
    } else if (re == b_y_im) {
      if (br > 0.0) {
        g_br = 0.5;
      } else {
        g_br = -0.5;
      }

      if (im > 0.0) {
        g_im = 0.5;
      } else {
        g_im = -0.5;
      }

      c_y_re = (c_y_re * g_br + y_im * g_im) / b_y_im;
    } else {
      b_y_im = br / im;
      c_y_re = (b_y_im * c_y_re + y_im) / (im + b_y_im * br);
    }
  }

  states[1] = ((((x12 + x15 * a_re_tmp) + a_im) - y_re) - b_y_re) - c_y_re;
  a_im = x19 * (a_re_tmp * a_re_tmp - a_im_tmp * a_im_tmp);
  if (x19 * (a_re_tmp * a_im_tmp + a_im_tmp * a_re_tmp) == 0.0) {
    r = a_im / 2.0;
  } else if (a_im == 0.0) {
    r = 0.0;
  } else {
    r = a_im / 2.0;
  }

  re = ((((12.0 * x03 - 12.0 * x13) + re_tmp * x06) + re_tmp * x16) + t_s_re_tmp
        * x09) - t_s_re_tmp * x19;
  im = ((im_tmp * x06 + im_tmp * x16) + t_s_im_tmp * x09) - t_s_im_tmp * x19;
  y_re = m_y_re * re - m_y_im * im;
  y_im = m_y_re * im + m_y_im * re;
  br = 2.0 * n_y_re;
  im = 2.0 * n_y_im;
  if (im == 0.0) {
    if (y_im == 0.0) {
      y_re /= br;
    } else if (y_re == 0.0) {
      y_re = 0.0;
    } else {
      y_re /= br;
    }
  } else if (br == 0.0) {
    if (y_re == 0.0) {
      y_re = y_im / im;
    } else if (y_im == 0.0) {
      y_re = 0.0;
    } else {
      y_re = y_im / im;
    }
  } else {
    b_y_im = fabs(br);
    re = fabs(im);
    if (b_y_im > re) {
      b_y_im = im / br;
      y_re = (y_re + b_y_im * y_im) / (br + b_y_im * im);
    } else if (re == b_y_im) {
      if (br > 0.0) {
        h_br = 0.5;
      } else {
        h_br = -0.5;
      }

      if (im > 0.0) {
        h_im = 0.5;
      } else {
        h_im = -0.5;
      }

      y_re = (y_re * h_br + y_im * h_im) / b_y_im;
    } else {
      b_y_im = br / im;
      y_re = (b_y_im * y_re + y_im) / (im + b_y_im * br);
    }
  }

  re = ((((20.0 * x03 - 20.0 * x13) + b_re_tmp * x06) + c_re_tmp * x16) +
        t_s_re_tmp * x09) - 3.0 * t_s_re_tmp * x19;
  im = ((b_im_tmp * x06 + c_im_tmp * x16) + t_s_im_tmp * x09) - 3.0 * t_s_im_tmp
    * x19;
  b_y_re = o_y_re * re - o_y_im * im;
  y_im = o_y_re * im + o_y_im * re;
  br = 2.0 * p_y_re;
  im = 2.0 * p_y_im;
  if (im == 0.0) {
    if (y_im == 0.0) {
      b_y_re /= br;
    } else if (b_y_re == 0.0) {
      b_y_re = 0.0;
    } else {
      b_y_re /= br;
    }
  } else if (br == 0.0) {
    if (b_y_re == 0.0) {
      b_y_re = y_im / im;
    } else if (y_im == 0.0) {
      b_y_re = 0.0;
    } else {
      b_y_re = y_im / im;
    }
  } else {
    b_y_im = fabs(br);
    re = fabs(im);
    if (b_y_im > re) {
      b_y_im = im / br;
      b_y_re = (b_y_re + b_y_im * y_im) / (br + b_y_im * im);
    } else if (re == b_y_im) {
      if (br > 0.0) {
        i_br = 0.5;
      } else {
        i_br = -0.5;
      }

      if (im > 0.0) {
        i_im = 0.5;
      } else {
        i_im = -0.5;
      }

      b_y_re = (b_y_re * i_br + y_im * i_im) / b_y_im;
    } else {
      b_y_im = br / im;
      b_y_re = (b_y_im * b_y_re + y_im) / (im + b_y_im * br);
    }
  }

  re = ((((30.0 * x03 - 30.0 * x13) + c_y_im * x06) + d_y_im * x16) + 2.0 *
        t_s_re_tmp * x09) - 3.0 * t_s_re_tmp * x19;
  im = ((d_y_re * x06 + d_im_tmp * x16) + 2.0 * t_s_im_tmp * x09) - 3.0 *
    t_s_im_tmp * x19;
  c_y_re = q_y_re * re - q_y_im * im;
  y_im = q_y_re * im + q_y_im * re;
  br = 2.0 * x_re;
  im = 2.0 * x_im;
  if (im == 0.0) {
    if (y_im == 0.0) {
      c_y_re /= br;
    } else if (c_y_re == 0.0) {
      c_y_re = 0.0;
    } else {
      c_y_re /= br;
    }
  } else if (br == 0.0) {
    if (c_y_re == 0.0) {
      c_y_re = y_im / im;
    } else if (y_im == 0.0) {
      c_y_re = 0.0;
    } else {
      c_y_re = y_im / im;
    }
  } else {
    b_y_im = fabs(br);
    re = fabs(im);
    if (b_y_im > re) {
      b_y_im = im / br;
      c_y_re = (c_y_re + b_y_im * y_im) / (br + b_y_im * im);
    } else if (re == b_y_im) {
      if (br > 0.0) {
        j_br = 0.5;
      } else {
        j_br = -0.5;
      }

      if (im > 0.0) {
        j_im = 0.5;
      } else {
        j_im = -0.5;
      }

      c_y_re = (c_y_re * j_br + y_im * j_im) / b_y_im;
    } else {
      b_y_im = br / im;
      c_y_re = (b_y_im * c_y_re + y_im) / (im + b_y_im * br);
    }
  }

  states[2] = ((((x13 + x16 * a_re_tmp) + r) - y_re) - b_y_re) - c_y_re;
}

/* End of code generation (quadf_state.c) */
