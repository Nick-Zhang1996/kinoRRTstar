/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * quadf_costPFF.c
 *
 * Code generation for function 'quadf_costPFF'
 *
 */

/* Include files */
#include <math.h>
#include "rt_nonfinite.h"
#include "benchmarkRRT.h"
#include "quadf_costPFF.h"
#include "quadf_cost.h"
#include "benchmarkRRT_rtwutil.h"

/* Function Definitions */
double quadf_costPFF(const creal_T t_s, double x01, double x02, double x03,
                     double x04, double x05, double x06, double x07, double x08,
                     double x09, double x11, double x12, double x13)
{
  double cost;
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
  double r_y_re;
  double r_y_im;
  double s_y_re;
  double s_y_im;
  double t_y_re;
  double t_y_im;
  double u_y_re;
  double u_y_im;
  double v_y_re;
  double v_y_im;
  double w_y_re;
  double w_y_im;
  double x_y_re;
  double x_y_im;
  double y_y_re;
  double y_y_im;
  double ab_y_re;
  double ab_y_im;
  double bb_y_re;
  double bb_y_im;
  double cb_y_re;
  double cb_y_im;
  double db_y_re;
  double db_y_im;
  double eb_y_re;
  double eb_y_im;
  double fb_y_re;
  double fb_y_im;
  double gb_y_re;
  double gb_y_im;
  double hb_y_re;
  double hb_y_im;
  double ib_y_re;
  double ib_y_im;
  double jb_y_re;
  double jb_y_im;
  double kb_y_re;
  double kb_y_im;
  double lb_y_re;
  double lb_y_im;
  double mb_y_re;
  double mb_y_im;
  double nb_y_re;
  double nb_y_im;
  double ob_y_re;
  double ob_y_im;
  double pb_y_re;
  double pb_y_im;
  double qb_y_re;
  double qb_y_im;
  double rb_y_re;
  double rb_y_im;
  double sb_y_re;
  double sb_y_im;
  double tb_y_re;
  double tb_y_im;
  double a_re;
  double a_im;
  double ar_tmp;
  double ar;
  double re;
  double im;
  double b_ar_tmp;
  double b_re;
  double b_im;
  double c_ar_tmp;
  double c_re;
  double c_im;
  double d_re;
  double d_im;
  double e_re;
  double e_im;
  double f_im;
  double f_re;
  double g_im;
  double g_re;
  double h_im;
  double d_ar_tmp;
  double h_re;
  double i_im;
  double e_ar_tmp;
  double i_re;
  double j_im;
  double f_ar_tmp;
  double j_re;
  double k_im;
  double k_re;
  double l_im;
  double l_re;
  double m_im;
  double m_re;
  double n_im;
  double re_tmp;
  double im_tmp;
  double n_re;
  double o_im;
  double b_re_tmp;
  double b_im_tmp;
  double o_re;
  double p_im;
  double p_re;
  double q_im;
  double q_re;
  double r_im;
  double r_re;
  double s_im;
  double b_a_re;
  double b_x_im;
  if ((t_s.im == 0.0) && (t_s.re >= 0.0)) {
    y_re = rt_powd_snf(t_s.re, 8.0);
    y_im = 0.0;
  } else if (t_s.re == 0.0) {
    y_re = rt_powd_snf(t_s.im, 8.0);
    y_im = 0.0;
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

    y_re = 8.0 * x_re;
    y_im = 8.0 * x_im;
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
    b_y_re = rt_powd_snf(t_s.re, 7.0);
    b_y_im = 0.0;
  } else if (t_s.re == 0.0) {
    b_y_re = 0.0;
    b_y_im = -rt_powd_snf(t_s.im, 7.0);
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

    b_y_re = 7.0 * x_re;
    b_y_im = 7.0 * x_im;
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
    c_y_re = rt_powd_snf(t_s.re, 6.0);
    c_y_im = 0.0;
  } else if (t_s.re == 0.0) {
    c_y_re = -rt_powd_snf(t_s.im, 6.0);
    c_y_im = 0.0;
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

    c_y_re = 6.0 * x_re;
    c_y_im = 6.0 * x_im;
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
    d_y_re = rt_powd_snf(t_s.re, 6.0);
    d_y_im = 0.0;
  } else if (t_s.re == 0.0) {
    d_y_re = -rt_powd_snf(t_s.im, 6.0);
    d_y_im = 0.0;
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

    d_y_re = 6.0 * x_re;
    d_y_im = 6.0 * x_im;
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
    e_y_re = rt_powd_snf(t_s.re, 6.0);
    e_y_im = 0.0;
  } else if (t_s.re == 0.0) {
    e_y_re = -rt_powd_snf(t_s.im, 6.0);
    e_y_im = 0.0;
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

    e_y_re = 6.0 * x_re;
    e_y_im = 6.0 * x_im;
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
    f_y_re = rt_powd_snf(t_s.re, 6.0);
    f_y_im = 0.0;
  } else if (t_s.re == 0.0) {
    f_y_re = -rt_powd_snf(t_s.im, 6.0);
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

    f_y_re = 6.0 * x_re;
    f_y_im = 6.0 * x_im;
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

  if ((t_s.im == 0.0) && (t_s.re >= 0.0)) {
    g_y_re = rt_powd_snf(t_s.re, 5.0);
    g_y_im = 0.0;
  } else if (t_s.re == 0.0) {
    g_y_re = 0.0;
    g_y_im = rt_powd_snf(t_s.im, 5.0);
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

    g_y_re = 5.0 * x_re;
    g_y_im = 5.0 * x_im;
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

  if ((t_s.im == 0.0) && (t_s.re >= 0.0)) {
    i_y_re = rt_powd_snf(t_s.re, 5.0);
    i_y_im = 0.0;
  } else if (t_s.re == 0.0) {
    i_y_re = 0.0;
    i_y_im = rt_powd_snf(t_s.im, 5.0);
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

    i_y_re = 5.0 * x_re;
    i_y_im = 5.0 * x_im;
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
    j_y_re = rt_powd_snf(t_s.re, 5.0);
    j_y_im = 0.0;
  } else if (t_s.re == 0.0) {
    j_y_re = 0.0;
    j_y_im = rt_powd_snf(t_s.im, 5.0);
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

    j_y_re = 5.0 * x_re;
    j_y_im = 5.0 * x_im;
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

  if ((t_s.im == 0.0) && (t_s.re >= 0.0)) {
    k_y_re = rt_powd_snf(t_s.re, 5.0);
    k_y_im = 0.0;
  } else if (t_s.re == 0.0) {
    k_y_re = 0.0;
    k_y_im = rt_powd_snf(t_s.im, 5.0);
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

    k_y_re = 5.0 * x_re;
    k_y_im = 5.0 * x_im;
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
    l_y_re = rt_powd_snf(t_s.re, 5.0);
    l_y_im = 0.0;
  } else if (t_s.re == 0.0) {
    l_y_re = 0.0;
    l_y_im = rt_powd_snf(t_s.im, 5.0);
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

    l_y_re = 5.0 * x_re;
    l_y_im = 5.0 * x_im;
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

  if ((t_s.im == 0.0) && (t_s.re >= 0.0)) {
    m_y_re = rt_powd_snf(t_s.re, 4.0);
    m_y_im = 0.0;
  } else if (t_s.re == 0.0) {
    m_y_re = rt_powd_snf(t_s.im, 4.0);
    m_y_im = 0.0;
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

    m_y_re = 4.0 * x_re;
    m_y_im = 4.0 * x_im;
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
    n_y_re = rt_powd_snf(t_s.re, 4.0);
    n_y_im = 0.0;
  } else if (t_s.re == 0.0) {
    n_y_re = rt_powd_snf(t_s.im, 4.0);
    n_y_im = 0.0;
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

    n_y_re = 4.0 * x_re;
    n_y_im = 4.0 * x_im;
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

  if ((t_s.im == 0.0) && (t_s.re >= 0.0)) {
    o_y_re = rt_powd_snf(t_s.re, 4.0);
    o_y_im = 0.0;
  } else if (t_s.re == 0.0) {
    o_y_re = rt_powd_snf(t_s.im, 4.0);
    o_y_im = 0.0;
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

    o_y_re = 4.0 * x_re;
    o_y_im = 4.0 * x_im;
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
    p_y_re = rt_powd_snf(t_s.re, 4.0);
    p_y_im = 0.0;
  } else if (t_s.re == 0.0) {
    p_y_re = rt_powd_snf(t_s.im, 4.0);
    p_y_im = 0.0;
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

    p_y_re = 4.0 * x_re;
    p_y_im = 4.0 * x_im;
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

  if ((t_s.im == 0.0) && (t_s.re >= 0.0)) {
    q_y_re = rt_powd_snf(t_s.re, 4.0);
    q_y_im = 0.0;
  } else if (t_s.re == 0.0) {
    q_y_re = rt_powd_snf(t_s.im, 4.0);
    q_y_im = 0.0;
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
    r_y_re = rt_powd_snf(t_s.re, 4.0);
    r_y_im = 0.0;
  } else if (t_s.re == 0.0) {
    r_y_re = rt_powd_snf(t_s.im, 4.0);
    r_y_im = 0.0;
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

    r_y_re = 4.0 * x_re;
    r_y_im = 4.0 * x_im;
    if (r_y_im == 0.0) {
      r_y_re = exp(r_y_re);
      r_y_im = 0.0;
    } else if (rtIsInf(r_y_im) && rtIsInf(r_y_re) && (r_y_re < 0.0)) {
      r_y_re = 0.0;
      r_y_im = 0.0;
    } else {
      r = exp(r_y_re / 2.0);
      r_y_re = r * (r * cos(r_y_im));
      r_y_im = r * (r * sin(r_y_im));
    }
  }

  if ((t_s.im == 0.0) && (t_s.re >= 0.0)) {
    s_y_re = rt_powd_snf(t_s.re, 4.0);
    s_y_im = 0.0;
  } else if (t_s.re == 0.0) {
    s_y_re = rt_powd_snf(t_s.im, 4.0);
    s_y_im = 0.0;
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

    s_y_re = 4.0 * x_re;
    s_y_im = 4.0 * x_im;
    if (s_y_im == 0.0) {
      s_y_re = exp(s_y_re);
      s_y_im = 0.0;
    } else if (rtIsInf(s_y_im) && rtIsInf(s_y_re) && (s_y_re < 0.0)) {
      s_y_re = 0.0;
      s_y_im = 0.0;
    } else {
      r = exp(s_y_re / 2.0);
      s_y_re = r * (r * cos(s_y_im));
      s_y_im = r * (r * sin(s_y_im));
    }
  }

  if ((t_s.im == 0.0) && (t_s.re >= 0.0)) {
    t_y_re = rt_powd_snf(t_s.re, 4.0);
    t_y_im = 0.0;
  } else if (t_s.re == 0.0) {
    t_y_re = rt_powd_snf(t_s.im, 4.0);
    t_y_im = 0.0;
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

    t_y_re = 4.0 * x_re;
    t_y_im = 4.0 * x_im;
    if (t_y_im == 0.0) {
      t_y_re = exp(t_y_re);
      t_y_im = 0.0;
    } else if (rtIsInf(t_y_im) && rtIsInf(t_y_re) && (t_y_re < 0.0)) {
      t_y_re = 0.0;
      t_y_im = 0.0;
    } else {
      r = exp(t_y_re / 2.0);
      t_y_re = r * (r * cos(t_y_im));
      t_y_im = r * (r * sin(t_y_im));
    }
  }

  if ((t_s.im == 0.0) && (t_s.re >= 0.0)) {
    u_y_re = rt_powd_snf(t_s.re, 4.0);
    u_y_im = 0.0;
  } else if (t_s.re == 0.0) {
    u_y_re = rt_powd_snf(t_s.im, 4.0);
    u_y_im = 0.0;
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

    u_y_re = 4.0 * x_re;
    u_y_im = 4.0 * x_im;
    if (u_y_im == 0.0) {
      u_y_re = exp(u_y_re);
      u_y_im = 0.0;
    } else if (rtIsInf(u_y_im) && rtIsInf(u_y_re) && (u_y_re < 0.0)) {
      u_y_re = 0.0;
      u_y_im = 0.0;
    } else {
      r = exp(u_y_re / 2.0);
      u_y_re = r * (r * cos(u_y_im));
      u_y_im = r * (r * sin(u_y_im));
    }
  }

  if ((t_s.im == 0.0) && (t_s.re >= 0.0)) {
    v_y_re = rt_powd_snf(t_s.re, 4.0);
    v_y_im = 0.0;
  } else if (t_s.re == 0.0) {
    v_y_re = rt_powd_snf(t_s.im, 4.0);
    v_y_im = 0.0;
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

    v_y_re = 4.0 * x_re;
    v_y_im = 4.0 * x_im;
    if (v_y_im == 0.0) {
      v_y_re = exp(v_y_re);
      v_y_im = 0.0;
    } else if (rtIsInf(v_y_im) && rtIsInf(v_y_re) && (v_y_re < 0.0)) {
      v_y_re = 0.0;
      v_y_im = 0.0;
    } else {
      r = exp(v_y_re / 2.0);
      v_y_re = r * (r * cos(v_y_im));
      v_y_im = r * (r * sin(v_y_im));
    }
  }

  if ((t_s.im == 0.0) && (t_s.re >= 0.0)) {
    w_y_re = rt_powd_snf(t_s.re, 4.0);
    w_y_im = 0.0;
  } else if (t_s.re == 0.0) {
    w_y_re = rt_powd_snf(t_s.im, 4.0);
    w_y_im = 0.0;
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

    w_y_re = 4.0 * x_re;
    w_y_im = 4.0 * x_im;
    if (w_y_im == 0.0) {
      w_y_re = exp(w_y_re);
      w_y_im = 0.0;
    } else if (rtIsInf(w_y_im) && rtIsInf(w_y_re) && (w_y_re < 0.0)) {
      w_y_re = 0.0;
      w_y_im = 0.0;
    } else {
      r = exp(w_y_re / 2.0);
      w_y_re = r * (r * cos(w_y_im));
      w_y_im = r * (r * sin(w_y_im));
    }
  }

  if ((t_s.im == 0.0) && (t_s.re >= 0.0)) {
    x_y_re = rt_powd_snf(t_s.re, 4.0);
    x_y_im = 0.0;
  } else if (t_s.re == 0.0) {
    x_y_re = rt_powd_snf(t_s.im, 4.0);
    x_y_im = 0.0;
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

    x_y_re = 4.0 * x_re;
    x_y_im = 4.0 * x_im;
    if (x_y_im == 0.0) {
      x_y_re = exp(x_y_re);
      x_y_im = 0.0;
    } else if (rtIsInf(x_y_im) && rtIsInf(x_y_re) && (x_y_re < 0.0)) {
      x_y_re = 0.0;
      x_y_im = 0.0;
    } else {
      r = exp(x_y_re / 2.0);
      x_y_re = r * (r * cos(x_y_im));
      x_y_im = r * (r * sin(x_y_im));
    }
  }

  if ((t_s.im == 0.0) && (t_s.re >= 0.0)) {
    y_y_re = rt_powd_snf(t_s.re, 4.0);
    y_y_im = 0.0;
  } else if (t_s.re == 0.0) {
    y_y_re = rt_powd_snf(t_s.im, 4.0);
    y_y_im = 0.0;
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

    y_y_re = 4.0 * x_re;
    y_y_im = 4.0 * x_im;
    if (y_y_im == 0.0) {
      y_y_re = exp(y_y_re);
      y_y_im = 0.0;
    } else if (rtIsInf(y_y_im) && rtIsInf(y_y_re) && (y_y_re < 0.0)) {
      y_y_re = 0.0;
      y_y_im = 0.0;
    } else {
      r = exp(y_y_re / 2.0);
      y_y_re = r * (r * cos(y_y_im));
      y_y_im = r * (r * sin(y_y_im));
    }
  }

  if ((t_s.im == 0.0) && (t_s.re >= 0.0)) {
    ab_y_re = rt_powd_snf(t_s.re, 4.0);
    ab_y_im = 0.0;
  } else if (t_s.re == 0.0) {
    ab_y_re = rt_powd_snf(t_s.im, 4.0);
    ab_y_im = 0.0;
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

    ab_y_re = 4.0 * x_re;
    ab_y_im = 4.0 * x_im;
    if (ab_y_im == 0.0) {
      ab_y_re = exp(ab_y_re);
      ab_y_im = 0.0;
    } else if (rtIsInf(ab_y_im) && rtIsInf(ab_y_re) && (ab_y_re < 0.0)) {
      ab_y_re = 0.0;
      ab_y_im = 0.0;
    } else {
      r = exp(ab_y_re / 2.0);
      ab_y_re = r * (r * cos(ab_y_im));
      ab_y_im = r * (r * sin(ab_y_im));
    }
  }

  if ((t_s.im == 0.0) && (t_s.re >= 0.0)) {
    bb_y_re = rt_powd_snf(t_s.re, 4.0);
    bb_y_im = 0.0;
  } else if (t_s.re == 0.0) {
    bb_y_re = rt_powd_snf(t_s.im, 4.0);
    bb_y_im = 0.0;
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

    bb_y_re = 4.0 * x_re;
    bb_y_im = 4.0 * x_im;
    if (bb_y_im == 0.0) {
      bb_y_re = exp(bb_y_re);
      bb_y_im = 0.0;
    } else if (rtIsInf(bb_y_im) && rtIsInf(bb_y_re) && (bb_y_re < 0.0)) {
      bb_y_re = 0.0;
      bb_y_im = 0.0;
    } else {
      r = exp(bb_y_re / 2.0);
      bb_y_re = r * (r * cos(bb_y_im));
      bb_y_im = r * (r * sin(bb_y_im));
    }
  }

  if ((t_s.im == 0.0) && (t_s.re >= 0.0)) {
    cb_y_re = rt_powd_snf(t_s.re, 3.0);
    cb_y_im = 0.0;
  } else if (t_s.re == 0.0) {
    cb_y_re = 0.0;
    cb_y_im = -rt_powd_snf(t_s.im, 3.0);
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

    cb_y_re = 3.0 * x_re;
    cb_y_im = 3.0 * x_im;
    if (cb_y_im == 0.0) {
      cb_y_re = exp(cb_y_re);
      cb_y_im = 0.0;
    } else if (rtIsInf(cb_y_im) && rtIsInf(cb_y_re) && (cb_y_re < 0.0)) {
      cb_y_re = 0.0;
      cb_y_im = 0.0;
    } else {
      r = exp(cb_y_re / 2.0);
      cb_y_re = r * (r * cos(cb_y_im));
      cb_y_im = r * (r * sin(cb_y_im));
    }
  }

  if ((t_s.im == 0.0) && (t_s.re >= 0.0)) {
    db_y_re = rt_powd_snf(t_s.re, 3.0);
    db_y_im = 0.0;
  } else if (t_s.re == 0.0) {
    db_y_re = 0.0;
    db_y_im = -rt_powd_snf(t_s.im, 3.0);
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

    db_y_re = 3.0 * x_re;
    db_y_im = 3.0 * x_im;
    if (db_y_im == 0.0) {
      db_y_re = exp(db_y_re);
      db_y_im = 0.0;
    } else if (rtIsInf(db_y_im) && rtIsInf(db_y_re) && (db_y_re < 0.0)) {
      db_y_re = 0.0;
      db_y_im = 0.0;
    } else {
      r = exp(db_y_re / 2.0);
      db_y_re = r * (r * cos(db_y_im));
      db_y_im = r * (r * sin(db_y_im));
    }
  }

  if ((t_s.im == 0.0) && (t_s.re >= 0.0)) {
    eb_y_re = rt_powd_snf(t_s.re, 3.0);
    eb_y_im = 0.0;
  } else if (t_s.re == 0.0) {
    eb_y_re = 0.0;
    eb_y_im = -rt_powd_snf(t_s.im, 3.0);
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

    eb_y_re = 3.0 * x_re;
    eb_y_im = 3.0 * x_im;
    if (eb_y_im == 0.0) {
      eb_y_re = exp(eb_y_re);
      eb_y_im = 0.0;
    } else if (rtIsInf(eb_y_im) && rtIsInf(eb_y_re) && (eb_y_re < 0.0)) {
      eb_y_re = 0.0;
      eb_y_im = 0.0;
    } else {
      r = exp(eb_y_re / 2.0);
      eb_y_re = r * (r * cos(eb_y_im));
      eb_y_im = r * (r * sin(eb_y_im));
    }
  }

  if ((t_s.im == 0.0) && (t_s.re >= 0.0)) {
    fb_y_re = rt_powd_snf(t_s.re, 3.0);
    fb_y_im = 0.0;
  } else if (t_s.re == 0.0) {
    fb_y_re = 0.0;
    fb_y_im = -rt_powd_snf(t_s.im, 3.0);
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

    fb_y_re = 3.0 * x_re;
    fb_y_im = 3.0 * x_im;
    if (fb_y_im == 0.0) {
      fb_y_re = exp(fb_y_re);
      fb_y_im = 0.0;
    } else if (rtIsInf(fb_y_im) && rtIsInf(fb_y_re) && (fb_y_re < 0.0)) {
      fb_y_re = 0.0;
      fb_y_im = 0.0;
    } else {
      r = exp(fb_y_re / 2.0);
      fb_y_re = r * (r * cos(fb_y_im));
      fb_y_im = r * (r * sin(fb_y_im));
    }
  }

  if ((t_s.im == 0.0) && (t_s.re >= 0.0)) {
    gb_y_re = rt_powd_snf(t_s.re, 3.0);
    gb_y_im = 0.0;
  } else if (t_s.re == 0.0) {
    gb_y_re = 0.0;
    gb_y_im = -rt_powd_snf(t_s.im, 3.0);
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

    gb_y_re = 3.0 * x_re;
    gb_y_im = 3.0 * x_im;
    if (gb_y_im == 0.0) {
      gb_y_re = exp(gb_y_re);
      gb_y_im = 0.0;
    } else if (rtIsInf(gb_y_im) && rtIsInf(gb_y_re) && (gb_y_re < 0.0)) {
      gb_y_re = 0.0;
      gb_y_im = 0.0;
    } else {
      r = exp(gb_y_re / 2.0);
      gb_y_re = r * (r * cos(gb_y_im));
      gb_y_im = r * (r * sin(gb_y_im));
    }
  }

  if ((t_s.im == 0.0) && (t_s.re >= 0.0)) {
    hb_y_re = rt_powd_snf(t_s.re, 3.0);
    hb_y_im = 0.0;
  } else if (t_s.re == 0.0) {
    hb_y_re = 0.0;
    hb_y_im = -rt_powd_snf(t_s.im, 3.0);
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

    hb_y_re = 3.0 * x_re;
    hb_y_im = 3.0 * x_im;
    if (hb_y_im == 0.0) {
      hb_y_re = exp(hb_y_re);
      hb_y_im = 0.0;
    } else if (rtIsInf(hb_y_im) && rtIsInf(hb_y_re) && (hb_y_re < 0.0)) {
      hb_y_re = 0.0;
      hb_y_im = 0.0;
    } else {
      r = exp(hb_y_re / 2.0);
      hb_y_re = r * (r * cos(hb_y_im));
      hb_y_im = r * (r * sin(hb_y_im));
    }
  }

  if ((t_s.im == 0.0) && (t_s.re >= 0.0)) {
    ib_y_re = rt_powd_snf(t_s.re, 3.0);
    ib_y_im = 0.0;
  } else if (t_s.re == 0.0) {
    ib_y_re = 0.0;
    ib_y_im = -rt_powd_snf(t_s.im, 3.0);
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

    ib_y_re = 3.0 * x_re;
    ib_y_im = 3.0 * x_im;
    if (ib_y_im == 0.0) {
      ib_y_re = exp(ib_y_re);
      ib_y_im = 0.0;
    } else if (rtIsInf(ib_y_im) && rtIsInf(ib_y_re) && (ib_y_re < 0.0)) {
      ib_y_re = 0.0;
      ib_y_im = 0.0;
    } else {
      r = exp(ib_y_re / 2.0);
      ib_y_re = r * (r * cos(ib_y_im));
      ib_y_im = r * (r * sin(ib_y_im));
    }
  }

  if ((t_s.im == 0.0) && (t_s.re >= 0.0)) {
    jb_y_re = rt_powd_snf(t_s.re, 3.0);
    jb_y_im = 0.0;
  } else if (t_s.re == 0.0) {
    jb_y_re = 0.0;
    jb_y_im = -rt_powd_snf(t_s.im, 3.0);
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

    jb_y_re = 3.0 * x_re;
    jb_y_im = 3.0 * x_im;
    if (jb_y_im == 0.0) {
      jb_y_re = exp(jb_y_re);
      jb_y_im = 0.0;
    } else if (rtIsInf(jb_y_im) && rtIsInf(jb_y_re) && (jb_y_re < 0.0)) {
      jb_y_re = 0.0;
      jb_y_im = 0.0;
    } else {
      r = exp(jb_y_re / 2.0);
      jb_y_re = r * (r * cos(jb_y_im));
      jb_y_im = r * (r * sin(jb_y_im));
    }
  }

  if ((t_s.im == 0.0) && (t_s.re >= 0.0)) {
    kb_y_re = rt_powd_snf(t_s.re, 3.0);
    kb_y_im = 0.0;
  } else if (t_s.re == 0.0) {
    kb_y_re = 0.0;
    kb_y_im = -rt_powd_snf(t_s.im, 3.0);
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

    kb_y_re = 3.0 * x_re;
    kb_y_im = 3.0 * x_im;
    if (kb_y_im == 0.0) {
      kb_y_re = exp(kb_y_re);
      kb_y_im = 0.0;
    } else if (rtIsInf(kb_y_im) && rtIsInf(kb_y_re) && (kb_y_re < 0.0)) {
      kb_y_re = 0.0;
      kb_y_im = 0.0;
    } else {
      r = exp(kb_y_re / 2.0);
      kb_y_re = r * (r * cos(kb_y_im));
      kb_y_im = r * (r * sin(kb_y_im));
    }
  }

  if ((t_s.im == 0.0) && (t_s.re >= 0.0)) {
    lb_y_re = rt_powd_snf(t_s.re, 3.0);
    lb_y_im = 0.0;
  } else if (t_s.re == 0.0) {
    lb_y_re = 0.0;
    lb_y_im = -rt_powd_snf(t_s.im, 3.0);
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

    lb_y_re = 3.0 * x_re;
    lb_y_im = 3.0 * x_im;
    if (lb_y_im == 0.0) {
      lb_y_re = exp(lb_y_re);
      lb_y_im = 0.0;
    } else if (rtIsInf(lb_y_im) && rtIsInf(lb_y_re) && (lb_y_re < 0.0)) {
      lb_y_re = 0.0;
      lb_y_im = 0.0;
    } else {
      r = exp(lb_y_re / 2.0);
      lb_y_re = r * (r * cos(lb_y_im));
      lb_y_im = r * (r * sin(lb_y_im));
    }
  }

  if ((t_s.im == 0.0) && (t_s.re >= 0.0)) {
    mb_y_re = rt_powd_snf(t_s.re, 3.0);
    mb_y_im = 0.0;
  } else if (t_s.re == 0.0) {
    mb_y_re = 0.0;
    mb_y_im = -rt_powd_snf(t_s.im, 3.0);
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

    mb_y_re = 3.0 * x_re;
    mb_y_im = 3.0 * x_im;
    if (mb_y_im == 0.0) {
      mb_y_re = exp(mb_y_re);
      mb_y_im = 0.0;
    } else if (rtIsInf(mb_y_im) && rtIsInf(mb_y_re) && (mb_y_re < 0.0)) {
      mb_y_re = 0.0;
      mb_y_im = 0.0;
    } else {
      r = exp(mb_y_re / 2.0);
      mb_y_re = r * (r * cos(mb_y_im));
      mb_y_im = r * (r * sin(mb_y_im));
    }
  }

  if ((t_s.im == 0.0) && (t_s.re >= 0.0)) {
    nb_y_re = rt_powd_snf(t_s.re, 3.0);
    nb_y_im = 0.0;
  } else if (t_s.re == 0.0) {
    nb_y_re = 0.0;
    nb_y_im = -rt_powd_snf(t_s.im, 3.0);
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

    nb_y_re = 3.0 * x_re;
    nb_y_im = 3.0 * x_im;
    if (nb_y_im == 0.0) {
      nb_y_re = exp(nb_y_re);
      nb_y_im = 0.0;
    } else if (rtIsInf(nb_y_im) && rtIsInf(nb_y_re) && (nb_y_re < 0.0)) {
      nb_y_re = 0.0;
      nb_y_im = 0.0;
    } else {
      r = exp(nb_y_re / 2.0);
      nb_y_re = r * (r * cos(nb_y_im));
      nb_y_im = r * (r * sin(nb_y_im));
    }
  }

  if ((t_s.im == 0.0) && (t_s.re >= 0.0)) {
    ob_y_re = rt_powd_snf(t_s.re, 3.0);
    ob_y_im = 0.0;
  } else if (t_s.re == 0.0) {
    ob_y_re = 0.0;
    ob_y_im = -rt_powd_snf(t_s.im, 3.0);
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

    ob_y_re = 3.0 * x_re;
    ob_y_im = 3.0 * x_im;
    if (ob_y_im == 0.0) {
      ob_y_re = exp(ob_y_re);
      ob_y_im = 0.0;
    } else if (rtIsInf(ob_y_im) && rtIsInf(ob_y_re) && (ob_y_re < 0.0)) {
      ob_y_re = 0.0;
      ob_y_im = 0.0;
    } else {
      r = exp(ob_y_re / 2.0);
      ob_y_re = r * (r * cos(ob_y_im));
      ob_y_im = r * (r * sin(ob_y_im));
    }
  }

  if ((t_s.im == 0.0) && (t_s.re >= 0.0)) {
    pb_y_re = rt_powd_snf(t_s.re, 3.0);
    pb_y_im = 0.0;
  } else if (t_s.re == 0.0) {
    pb_y_re = 0.0;
    pb_y_im = -rt_powd_snf(t_s.im, 3.0);
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

    pb_y_re = 3.0 * x_re;
    pb_y_im = 3.0 * x_im;
    if (pb_y_im == 0.0) {
      pb_y_re = exp(pb_y_re);
      pb_y_im = 0.0;
    } else if (rtIsInf(pb_y_im) && rtIsInf(pb_y_re) && (pb_y_re < 0.0)) {
      pb_y_re = 0.0;
      pb_y_im = 0.0;
    } else {
      r = exp(pb_y_re / 2.0);
      pb_y_re = r * (r * cos(pb_y_im));
      pb_y_im = r * (r * sin(pb_y_im));
    }
  }

  if ((t_s.im == 0.0) && (t_s.re >= 0.0)) {
    qb_y_re = rt_powd_snf(t_s.re, 3.0);
    qb_y_im = 0.0;
  } else if (t_s.re == 0.0) {
    qb_y_re = 0.0;
    qb_y_im = -rt_powd_snf(t_s.im, 3.0);
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

    qb_y_re = 3.0 * x_re;
    qb_y_im = 3.0 * x_im;
    if (qb_y_im == 0.0) {
      qb_y_re = exp(qb_y_re);
      qb_y_im = 0.0;
    } else if (rtIsInf(qb_y_im) && rtIsInf(qb_y_re) && (qb_y_re < 0.0)) {
      qb_y_re = 0.0;
      qb_y_im = 0.0;
    } else {
      r = exp(qb_y_re / 2.0);
      qb_y_re = r * (r * cos(qb_y_im));
      qb_y_im = r * (r * sin(qb_y_im));
    }
  }

  if ((t_s.im == 0.0) && (t_s.re >= 0.0)) {
    rb_y_re = rt_powd_snf(t_s.re, 3.0);
    rb_y_im = 0.0;
  } else if (t_s.re == 0.0) {
    rb_y_re = 0.0;
    rb_y_im = -rt_powd_snf(t_s.im, 3.0);
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

    rb_y_re = 3.0 * x_re;
    rb_y_im = 3.0 * x_im;
    if (rb_y_im == 0.0) {
      rb_y_re = exp(rb_y_re);
      rb_y_im = 0.0;
    } else if (rtIsInf(rb_y_im) && rtIsInf(rb_y_re) && (rb_y_re < 0.0)) {
      rb_y_re = 0.0;
      rb_y_im = 0.0;
    } else {
      r = exp(rb_y_re / 2.0);
      rb_y_re = r * (r * cos(rb_y_im));
      rb_y_im = r * (r * sin(rb_y_im));
    }
  }

  if ((t_s.im == 0.0) && (t_s.re >= 0.0)) {
    sb_y_re = rt_powd_snf(t_s.re, 3.0);
    sb_y_im = 0.0;
  } else if (t_s.re == 0.0) {
    sb_y_re = 0.0;
    sb_y_im = -rt_powd_snf(t_s.im, 3.0);
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

    sb_y_re = 3.0 * x_re;
    sb_y_im = 3.0 * x_im;
    if (sb_y_im == 0.0) {
      sb_y_re = exp(sb_y_re);
      sb_y_im = 0.0;
    } else if (rtIsInf(sb_y_im) && rtIsInf(sb_y_re) && (sb_y_re < 0.0)) {
      sb_y_re = 0.0;
      sb_y_im = 0.0;
    } else {
      r = exp(sb_y_re / 2.0);
      sb_y_re = r * (r * cos(sb_y_im));
      sb_y_im = r * (r * sin(sb_y_im));
    }
  }

  if ((t_s.im == 0.0) && (t_s.re >= 0.0)) {
    tb_y_re = rt_powd_snf(t_s.re, 3.0);
    tb_y_im = 0.0;
  } else if (t_s.re == 0.0) {
    tb_y_re = 0.0;
    tb_y_im = -rt_powd_snf(t_s.im, 3.0);
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

    tb_y_re = 3.0 * x_re;
    tb_y_im = 3.0 * x_im;
    if (tb_y_im == 0.0) {
      tb_y_re = exp(tb_y_re);
      tb_y_im = 0.0;
    } else if (rtIsInf(tb_y_im) && rtIsInf(tb_y_re) && (tb_y_re < 0.0)) {
      tb_y_re = 0.0;
      tb_y_im = 0.0;
    } else {
      r = exp(tb_y_re / 2.0);
      tb_y_re = r * (r * cos(tb_y_im));
      tb_y_im = r * (r * sin(tb_y_im));
    }
  }

  if ((t_s.im == 0.0) && (t_s.re >= 0.0)) {
    x_re = rt_powd_snf(t_s.re, 5.0);
    x_im = 0.0;
  } else if (t_s.re == 0.0) {
    x_re = 0.0;
    x_im = rt_powd_snf(t_s.im, 5.0);
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

    x_re *= 5.0;
    x_im *= 5.0;
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

  a_re = 500.0 * t_s.re + 9.0;
  a_im = 500.0 * t_s.im;
  ar_tmp = x07 * x07;
  ar = 81.0 * v_y_re * ar_tmp;
  v_y_im = 81.0 * v_y_im * ar_tmp;
  if (v_y_im == 0.0) {
    re = ar / 1000.0;
    im = 0.0;
  } else if (ar == 0.0) {
    re = 0.0;
    im = v_y_im / 1000.0;
  } else {
    re = ar / 1000.0;
    im = v_y_im / 1000.0;
  }

  b_ar_tmp = x08 * x08;
  ar = 81.0 * x_y_re * b_ar_tmp;
  v_y_im = 81.0 * x_y_im * b_ar_tmp;
  if (v_y_im == 0.0) {
    b_re = ar / 1000.0;
    b_im = 0.0;
  } else if (ar == 0.0) {
    b_re = 0.0;
    b_im = v_y_im / 1000.0;
  } else {
    b_re = ar / 1000.0;
    b_im = v_y_im / 1000.0;
  }

  c_ar_tmp = x09 * x09;
  ar = 81.0 * ab_y_re * c_ar_tmp;
  v_y_im = 81.0 * ab_y_im * c_ar_tmp;
  if (v_y_im == 0.0) {
    c_re = ar / 1000.0;
    c_im = 0.0;
  } else if (ar == 0.0) {
    c_re = 0.0;
    c_im = v_y_im / 1000.0;
  } else {
    c_re = ar / 1000.0;
    c_im = v_y_im / 1000.0;
  }

  ar = x04 * (81.0 * jb_y_re) * x07;
  v_y_im = x04 * (81.0 * jb_y_im) * x07;
  if (v_y_im == 0.0) {
    d_re = ar / 250.0;
    d_im = 0.0;
  } else if (ar == 0.0) {
    d_re = 0.0;
    d_im = v_y_im / 250.0;
  } else {
    d_re = ar / 250.0;
    d_im = v_y_im / 250.0;
  }

  ar = x05 * (81.0 * mb_y_re) * x08;
  v_y_im = x05 * (81.0 * mb_y_im) * x08;
  if (v_y_im == 0.0) {
    e_re = ar / 250.0;
    e_im = 0.0;
  } else if (ar == 0.0) {
    e_re = 0.0;
    e_im = v_y_im / 250.0;
  } else {
    e_re = ar / 250.0;
    e_im = v_y_im / 250.0;
  }

  ar = x06 * (81.0 * pb_y_re) * x09;
  v_y_im = x06 * (81.0 * pb_y_im) * x09;
  if (v_y_im == 0.0) {
    jb_y_im = ar / 250.0;
    mb_y_re = 0.0;
  } else if (ar == 0.0) {
    jb_y_im = 0.0;
    mb_y_re = v_y_im / 250.0;
  } else {
    jb_y_im = ar / 250.0;
    mb_y_re = v_y_im / 250.0;
  }

  mb_y_im = t_s.re * t_s.re - t_s.im * t_s.im;
  pb_y_re = t_s.re * t_s.im + t_s.im * t_s.re;
  ar = x01 * (81.0 * mb_y_im) * x07;
  v_y_im = x01 * (81.0 * pb_y_re) * x07;
  if (v_y_im == 0.0) {
    pb_y_im = ar / 250.0;
    f_im = 0.0;
  } else if (ar == 0.0) {
    pb_y_im = 0.0;
    f_im = v_y_im / 250.0;
  } else {
    pb_y_im = ar / 250.0;
    f_im = v_y_im / 250.0;
  }

  ar = x02 * (81.0 * mb_y_im) * x08;
  v_y_im = x02 * (81.0 * pb_y_re) * x08;
  if (v_y_im == 0.0) {
    f_re = ar / 250.0;
    g_im = 0.0;
  } else if (ar == 0.0) {
    f_re = 0.0;
    g_im = v_y_im / 250.0;
  } else {
    f_re = ar / 250.0;
    g_im = v_y_im / 250.0;
  }

  ar = x03 * (81.0 * mb_y_im) * x09;
  v_y_im = x03 * (81.0 * pb_y_re) * x09;
  if (v_y_im == 0.0) {
    g_re = ar / 250.0;
    h_im = 0.0;
  } else if (ar == 0.0) {
    g_re = 0.0;
    h_im = v_y_im / 250.0;
  } else {
    g_re = ar / 250.0;
    h_im = v_y_im / 250.0;
  }

  d_ar_tmp = x04 * x04;
  ar = 81.0 * mb_y_im * d_ar_tmp;
  v_y_im = 81.0 * pb_y_re * d_ar_tmp;
  if (v_y_im == 0.0) {
    h_re = ar / 250.0;
    i_im = 0.0;
  } else if (ar == 0.0) {
    h_re = 0.0;
    i_im = v_y_im / 250.0;
  } else {
    h_re = ar / 250.0;
    i_im = v_y_im / 250.0;
  }

  e_ar_tmp = x05 * x05;
  ar = 81.0 * mb_y_im * e_ar_tmp;
  v_y_im = 81.0 * pb_y_re * e_ar_tmp;
  if (v_y_im == 0.0) {
    i_re = ar / 250.0;
    j_im = 0.0;
  } else if (ar == 0.0) {
    i_re = 0.0;
    j_im = v_y_im / 250.0;
  } else {
    i_re = ar / 250.0;
    j_im = v_y_im / 250.0;
  }

  f_ar_tmp = x06 * x06;
  ar = 81.0 * mb_y_im * f_ar_tmp;
  v_y_im = 81.0 * pb_y_re * f_ar_tmp;
  if (v_y_im == 0.0) {
    j_re = ar / 250.0;
    k_im = 0.0;
  } else if (ar == 0.0) {
    j_re = 0.0;
    k_im = v_y_im / 250.0;
  } else {
    j_re = ar / 250.0;
    k_im = v_y_im / 250.0;
  }

  ar = x07 * (81.0 * mb_y_im) * x11;
  v_y_im = x07 * (81.0 * pb_y_re) * x11;
  if (v_y_im == 0.0) {
    k_re = ar / 250.0;
    l_im = 0.0;
  } else if (ar == 0.0) {
    k_re = 0.0;
    l_im = v_y_im / 250.0;
  } else {
    k_re = ar / 250.0;
    l_im = v_y_im / 250.0;
  }

  ar = x08 * (81.0 * mb_y_im) * x12;
  v_y_im = x08 * (81.0 * pb_y_re) * x12;
  if (v_y_im == 0.0) {
    l_re = ar / 250.0;
    m_im = 0.0;
  } else if (ar == 0.0) {
    l_re = 0.0;
    m_im = v_y_im / 250.0;
  } else {
    l_re = ar / 250.0;
    m_im = v_y_im / 250.0;
  }

  ar = x09 * (81.0 * mb_y_im) * x13;
  v_y_im = x09 * (81.0 * pb_y_re) * x13;
  if (v_y_im == 0.0) {
    m_re = ar / 250.0;
    n_im = 0.0;
  } else if (ar == 0.0) {
    m_re = 0.0;
    n_im = v_y_im / 250.0;
  } else {
    m_re = ar / 250.0;
    n_im = v_y_im / 250.0;
  }

  re_tmp = 36.0 * t_s.re;
  im_tmp = 36.0 * t_s.im;
  r = 81.0 * t_s.re;
  v_y_re = 81.0 * t_s.im;
  ar = x01 * r * x04;
  v_y_im = x01 * v_y_re * x04;
  if (v_y_im == 0.0) {
    n_re = ar / 125.0;
    o_im = 0.0;
  } else if (ar == 0.0) {
    n_re = 0.0;
    o_im = v_y_im / 125.0;
  } else {
    n_re = ar / 125.0;
    o_im = v_y_im / 125.0;
  }

  b_re_tmp = 72.0 * t_s.re;
  b_im_tmp = 72.0 * t_s.im;
  ar = x02 * r * x05;
  v_y_im = x02 * v_y_re * x05;
  if (v_y_im == 0.0) {
    o_re = ar / 125.0;
    p_im = 0.0;
  } else if (ar == 0.0) {
    o_re = 0.0;
    p_im = v_y_im / 125.0;
  } else {
    o_re = ar / 125.0;
    p_im = v_y_im / 125.0;
  }

  ar = x03 * r * x06;
  v_y_im = x03 * v_y_re * x06;
  if (v_y_im == 0.0) {
    p_re = ar / 125.0;
    q_im = 0.0;
  } else if (ar == 0.0) {
    p_re = 0.0;
    q_im = v_y_im / 125.0;
  } else {
    p_re = ar / 125.0;
    q_im = v_y_im / 125.0;
  }

  ar = x04 * r * x11;
  v_y_im = x04 * v_y_re * x11;
  if (v_y_im == 0.0) {
    q_re = ar / 125.0;
    r_im = 0.0;
  } else if (ar == 0.0) {
    q_re = 0.0;
    r_im = v_y_im / 125.0;
  } else {
    q_re = ar / 125.0;
    r_im = v_y_im / 125.0;
  }

  ar = x05 * r * x12;
  v_y_im = x05 * v_y_re * x12;
  if (v_y_im == 0.0) {
    r_re = ar / 125.0;
    s_im = 0.0;
  } else if (ar == 0.0) {
    r_re = 0.0;
    s_im = v_y_im / 125.0;
  } else {
    r_re = ar / 125.0;
    s_im = v_y_im / 125.0;
  }

  ar = x06 * r * x13;
  v_y_im = x06 * v_y_re * x13;
  if (v_y_im == 0.0) {
    ab_y_im = ar / 125.0;
    jb_y_re = 0.0;
  } else if (ar == 0.0) {
    ab_y_im = 0.0;
    jb_y_re = v_y_im / 125.0;
  } else {
    ab_y_im = ar / 125.0;
    jb_y_re = v_y_im / 125.0;
  }

  r = a_re * a_re - a_im * a_im;
  a_im = a_re * a_im + a_im * a_re;
  a_re = x_re * r - x_im * a_im;
  x_im = x_re * a_im + x_im * r;
  r = x01 * x01;
  v_y_re = x02 * x02;
  v_y_im = x03 * x03;
  x_y_re = x11 * x11;
  x_y_im = x12 * x12;
  ab_y_re = x13 * x13;
  ar =
    (((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((
    ((((((((((((((250000.0 * y_re + 9000.0 * b_y_re) + 300.0 * c_y_re * ar_tmp)
    + 300.0 * d_y_re * b_ar_tmp) + 300.0 * e_y_re * c_ar_tmp) + 81.0 * f_y_re) +
    x04 * (1500.0 * g_y_re) * x07) + x05 * (1500.0 * h_y_re) * x08) + x06 *
    (1500.0 * i_y_re) * x09) + 9.0 * j_y_re * ar_tmp) + 9.0 * k_y_re * b_ar_tmp)
    + 9.0 * l_y_re * c_ar_tmp) + x01 * (1500.0 * m_y_re) * x07) + x02 * (1500.0 *
    n_y_re) * x08) + x03 * (1500.0 * o_y_re) * x09) + 2250.0 * p_y_re * d_ar_tmp)
    + x04 * (36.0 * q_y_re) * x07) + 2250.0 * r_y_re * e_ar_tmp) + x05 * (36.0 *
    s_y_re) * x08) + 2250.0 * t_y_re * f_ar_tmp) + x06 * (36.0 * u_y_re) * x09)
    + re) - x07 * (1500.0 * w_y_re) * x11) + b_re) - x08 * (1500.0 * y_y_re) *
    x12) + c_re) - x09 * (1500.0 * bb_y_re) * x13) + x01 * (4500.0 * cb_y_re) *
    x04) + x01 * (36.0 * db_y_re) * x07) + x02 * (4500.0 * eb_y_re) * x05) + x02
    * (36.0 * fb_y_re) * x08) + x03 * (4500.0 * gb_y_re) * x06) + x03 * (36.0 *
    hb_y_re) * x09) + 36.0 * ib_y_re * d_ar_tmp) + d_re) - x04 * (4500.0 *
    kb_y_re) * x11) + 36.0 * lb_y_re * e_ar_tmp) + e_re) - x05 * (4500.0 *
    nb_y_re) * x12) + 36.0 * ob_y_re * f_ar_tmp) + jb_y_im) - x06 * (4500.0 *
    qb_y_re) * x13) - x07 * (36.0 * rb_y_re) * x11) - x08 * (36.0 * sb_y_re) *
    x12) - x09 * (36.0 * tb_y_re) * x13) + 2250.0 * mb_y_im * r) + x01 * (72.0 *
    mb_y_im) * x04) + pb_y_im) - x01 * (4500.0 * mb_y_im) * x11) + 2250.0 *
    mb_y_im * v_y_re) + x02 * (72.0 * mb_y_im) * x05) + f_re) - x02 * (4500.0 *
    mb_y_im) * x12) + 2250.0 * mb_y_im * v_y_im) + x03 * (72.0 * mb_y_im) * x06)
    + g_re) - x03 * (4500.0 * mb_y_im) * x13) + h_re) - x04 * (72.0 * mb_y_im) *
    x11) + i_re) - x05 * (72.0 * mb_y_im) * x12) + j_re) - x06 * (72.0 * mb_y_im)
    * x13) - k_re) - l_re) - m_re) + 2250.0 * mb_y_im * x_y_re) + 2250.0 *
    mb_y_im * x_y_im) + 2250.0 * mb_y_im * ab_y_re) + re_tmp * r) + n_re) - x01 *
    b_re_tmp * x11) + re_tmp * v_y_re) + o_re) - x02 * b_re_tmp * x12) + re_tmp *
                     v_y_im) + p_re) - x03 * b_re_tmp * x13) - q_re) - r_re) -
                ab_y_im) + re_tmp * x_y_re) + re_tmp * x_y_im) + re_tmp *
             ab_y_re) + 81.0 * r / 250.0) - 81.0 * x01 * x11 / 125.0) + 81.0 *
          v_y_re / 250.0) - 81.0 * x02 * x12 / 125.0) + 81.0 * v_y_im / 250.0) -
       81.0 * x03 * x13 / 125.0) + 81.0 * x_y_re / 250.0) + 81.0 * x_y_im /
     250.0) + 81.0 * ab_y_re / 250.0;
  v_y_im =
    (((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((
    (((((250000.0 * y_im + 9000.0 * b_y_im) + 300.0 * c_y_im * ar_tmp) + 300.0 *
    d_y_im * b_ar_tmp) + 300.0 * e_y_im * c_ar_tmp) + 81.0 * f_y_im) + x04 *
    (1500.0 * g_y_im) * x07) + x05 * (1500.0 * h_y_im) * x08) + x06 * (1500.0 *
    i_y_im) * x09) + 9.0 * j_y_im * ar_tmp) + 9.0 * k_y_im * b_ar_tmp) + 9.0 *
    l_y_im * c_ar_tmp) + x01 * (1500.0 * m_y_im) * x07) + x02 * (1500.0 * n_y_im)
    * x08) + x03 * (1500.0 * o_y_im) * x09) + 2250.0 * p_y_im * d_ar_tmp) + x04 *
    (36.0 * q_y_im) * x07) + 2250.0 * r_y_im * e_ar_tmp) + x05 * (36.0 * s_y_im)
    * x08) + 2250.0 * t_y_im * f_ar_tmp) + x06 * (36.0 * u_y_im) * x09) + im) -
    x07 * (1500.0 * w_y_im) * x11) + b_im) - x08 * (1500.0 * y_y_im) * x12) +
    c_im) - x09 * (1500.0 * bb_y_im) * x13) + x01 * (4500.0 * cb_y_im) * x04) +
    x01 * (36.0 * db_y_im) * x07) + x02 * (4500.0 * eb_y_im) * x05) + x02 *
    (36.0 * fb_y_im) * x08) + x03 * (4500.0 * gb_y_im) * x06) + x03 * (36.0 *
    hb_y_im) * x09) + 36.0 * ib_y_im * d_ar_tmp) + d_im) - x04 * (4500.0 *
    kb_y_im) * x11) + 36.0 * lb_y_im * e_ar_tmp) + e_im) - x05 * (4500.0 *
    nb_y_im) * x12) + 36.0 * ob_y_im * f_ar_tmp) + mb_y_re) - x06 * (4500.0 *
    qb_y_im) * x13) - x07 * (36.0 * rb_y_im) * x11) - x08 * (36.0 * sb_y_im) *
    x12) - x09 * (36.0 * tb_y_im) * x13) + 2250.0 * pb_y_re * r) + x01 * (72.0 *
    pb_y_re) * x04) + f_im) - x01 * (4500.0 * pb_y_re) * x11) + 2250.0 * pb_y_re
    * v_y_re) + x02 * (72.0 * pb_y_re) * x05) + g_im) - x02 * (4500.0 * pb_y_re)
    * x12) + 2250.0 * pb_y_re * v_y_im) + x03 * (72.0 * pb_y_re) * x06) + h_im)
    - x03 * (4500.0 * pb_y_re) * x13) + i_im) - x04 * (72.0 * pb_y_re) * x11) +
    j_im) - x05 * (72.0 * pb_y_re) * x12) + k_im) - x06 * (72.0 * pb_y_re) * x13)
    - l_im) - m_im) - n_im) + 2250.0 * pb_y_re * x_y_re) + 2250.0 * pb_y_re *
                    x_y_im) + 2250.0 * pb_y_re * ab_y_re) + im_tmp * r) + o_im)
                - x01 * b_im_tmp * x11) + im_tmp * v_y_re) + p_im) - x02 *
             b_im_tmp * x12) + im_tmp * v_y_im) + q_im) - x03 * b_im_tmp * x13)
         - r_im) - s_im) - jb_y_re) + im_tmp * x_y_re) + im_tmp * x_y_im) +
    im_tmp * ab_y_re;
  if (x_im == 0.0) {
    if (v_y_im == 0.0) {
      cost = ar / a_re;
    } else if (ar == 0.0) {
      cost = 0.0;
    } else {
      cost = ar / a_re;
    }
  } else if (a_re == 0.0) {
    if (ar == 0.0) {
      cost = v_y_im / x_im;
    } else if (v_y_im == 0.0) {
      cost = 0.0;
    } else {
      cost = v_y_im / x_im;
    }
  } else {
    r = fabs(a_re);
    v_y_re = fabs(x_im);
    if (r > v_y_re) {
      r = x_im / a_re;
      cost = (ar + r * v_y_im) / (a_re + r * x_im);
    } else if (v_y_re == r) {
      if (a_re > 0.0) {
        b_a_re = 0.5;
      } else {
        b_a_re = -0.5;
      }

      if (x_im > 0.0) {
        b_x_im = 0.5;
      } else {
        b_x_im = -0.5;
      }

      cost = (ar * b_a_re + v_y_im * b_x_im) / r;
    } else {
      r = a_re / x_im;
      cost = (r * ar + v_y_im) / (x_im + r * a_re);
    }
  }

  return cost;
}

/* End of code generation (quadf_costPFF.c) */
