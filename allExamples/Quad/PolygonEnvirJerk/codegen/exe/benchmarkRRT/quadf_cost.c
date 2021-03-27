/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * quadf_cost.c
 *
 * Code generation for function 'quadf_cost'
 *
 */

/* Include files */
#include <math.h>
#include "rt_nonfinite.h"
#include "benchmarkRRT.h"
#include "quadf_cost.h"
#include "benchmarkRRT_rtwutil.h"

/* Function Definitions */
double quadf_cost(const creal_T t_s, double x01, double x02, double x03, double
                  x04, double x05, double x06, double x07, double x08, double
                  x09, double x11, double x12, double x13, double x14, double
                  x15, double x16, double x17, double x18, double x19)
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
  double ar;
  double br;
  double bi;
  double re;
  double brm;
  double bim;
  double b_br;
  double b_re;
  double b_bi;
  double c_br;
  double c_re;
  double c_bi;
  double d_br;
  double d_re;
  double d_bi;
  double qb_y_re;
  double e_re;
  double qb_y_im;
  double e_br;
  double f_re;
  double e_bi;
  double g_re;
  double f_br;
  double f_bi;
  double rb_y_re;
  double h_re;
  double rb_y_im;
  double i_re;
  double g_br;
  double g_bi;
  double sb_y_re;
  double j_re;
  double sb_y_im;
  double h_br;
  double k_re;
  double h_bi;
  double i_br;
  double l_re;
  double i_bi;
  double m_re;
  double j_br;
  double j_bi;
  double tb_y_re;
  double n_re;
  double tb_y_im;
  double k_br;
  double o_re;
  double k_bi;
  double p_re;
  double l_br;
  double l_bi;
  double ub_y_re;
  double q_re;
  double ub_y_im;
  double r_re;
  double m_br;
  double m_bi;
  double vb_y_re;
  double s_re;
  double vb_y_im;
  double n_br;
  double t_re;
  double n_bi;
  double o_br;
  double u_re;
  double o_bi;
  double t_s_im_tmp;
  double p_br;
  double v_re;
  double p_bi;
  double q_br;
  double w_re;
  double q_bi;
  double r_br;
  double ab_re;
  double r_bi;
  double s_br;
  double bb_re;
  double s_bi;
  double t_br;
  double t_bi;
  double u_br;
  double cb_re;
  double u_bi;
  double v_br;
  double v_bi;
  double w_br;
  double w_bi;
  double x_br;
  double x_bi;
  double y_br;
  double y_bi;
  double ab_br;
  double ab_bi;
  double bb_br;
  double bb_bi;
  double cb_br;
  double cb_bi;
  double db_br;
  double db_bi;
  double eb_br;
  double eb_bi;
  double fb_br;
  double fb_bi;
  double gb_br;
  double gb_bi;
  double hb_br;
  double hb_bi;
  double ib_br;
  double ib_bi;
  double jb_br;
  double jb_bi;
  double kb_br;
  double kb_bi;
  double lb_br;
  double lb_bi;
  double mb_br;
  double mb_bi;
  double nb_br;
  double nb_bi;
  double xb_y_re;
  double wb_y_im;
  double ob_br;
  double ob_bi;
  double pb_br;
  double pb_bi;
  double qb_br;
  double qb_bi;
  double rb_br;
  double rb_bi;
  double sb_br;
  double sb_bi;
  double yb_y_re;
  double xb_y_im;
  double tb_br;
  double tb_bi;
  double ub_br;
  double ub_bi;
  double ac_y_re;
  double yb_y_im;
  double vb_br;
  double vb_bi;
  double wb_br;
  double wb_bi;
  double xb_br;
  double xb_bi;
  double yb_br;
  double yb_bi;
  double ac_br;
  double ac_bi;
  double bc_br;
  double bc_bi;
  double cc_br;
  double cc_bi;
  double dc_br;
  double dc_bi;
  double ec_br;
  double ec_bi;
  if ((t_s.im == 0.0) && (t_s.re >= 0.0)) {
    y_re = rt_powd_snf(t_s.re, 5.0);
    y_im = 0.0;
  } else if (t_s.re == 0.0) {
    y_re = 0.0;
    y_im = rt_powd_snf(t_s.im, 5.0);
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

    y_re = 5.0 * x_re;
    y_im = 5.0 * x_im;
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
    c_y_re = rt_powd_snf(t_s.re, 5.0);
    c_y_im = 0.0;
  } else if (t_s.re == 0.0) {
    c_y_re = 0.0;
    c_y_im = rt_powd_snf(t_s.im, 5.0);
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

    c_y_re = 5.0 * x_re;
    c_y_im = 5.0 * x_im;
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
    e_y_re = rt_powd_snf(t_s.re, 5.0);
    e_y_im = 0.0;
  } else if (t_s.re == 0.0) {
    e_y_re = 0.0;
    e_y_im = rt_powd_snf(t_s.im, 5.0);
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

    e_y_re = 5.0 * x_re;
    e_y_im = 5.0 * x_im;
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
    f_y_re = rt_powd_snf(t_s.re, 3.0);
    f_y_im = 0.0;
  } else if (t_s.re == 0.0) {
    f_y_re = 0.0;
    f_y_im = -rt_powd_snf(t_s.im, 3.0);
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

    f_y_re = 3.0 * x_re;
    f_y_im = 3.0 * x_im;
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
    h_y_re = rt_powd_snf(t_s.re, 3.0);
    h_y_im = 0.0;
  } else if (t_s.re == 0.0) {
    h_y_re = 0.0;
    h_y_im = -rt_powd_snf(t_s.im, 3.0);
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

    h_y_re = 3.0 * x_re;
    h_y_im = 3.0 * x_im;
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
    l_y_re = rt_powd_snf(t_s.re, 3.0);
    l_y_im = 0.0;
  } else if (t_s.re == 0.0) {
    l_y_re = 0.0;
    l_y_im = -rt_powd_snf(t_s.im, 3.0);
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

    l_y_re = 3.0 * x_re;
    l_y_im = 3.0 * x_im;
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
    n_y_re = rt_powd_snf(t_s.re, 3.0);
    n_y_im = 0.0;
  } else if (t_s.re == 0.0) {
    n_y_re = 0.0;
    n_y_im = -rt_powd_snf(t_s.im, 3.0);
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

    n_y_re = 3.0 * x_re;
    n_y_im = 3.0 * x_im;
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
    r_y_re = rt_powd_snf(t_s.re, 3.0);
    r_y_im = 0.0;
  } else if (t_s.re == 0.0) {
    r_y_re = 0.0;
    r_y_im = -rt_powd_snf(t_s.im, 3.0);
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

    r_y_re = 3.0 * x_re;
    r_y_im = 3.0 * x_im;
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
    s_y_re = rt_powd_snf(t_s.re, 5.0);
    s_y_im = 0.0;
  } else if (t_s.re == 0.0) {
    s_y_re = 0.0;
    s_y_im = rt_powd_snf(t_s.im, 5.0);
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

    s_y_re = 5.0 * x_re;
    s_y_im = 5.0 * x_im;
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
    v_y_re = rt_powd_snf(t_s.re, 5.0);
    v_y_im = 0.0;
  } else if (t_s.re == 0.0) {
    v_y_re = 0.0;
    v_y_im = rt_powd_snf(t_s.im, 5.0);
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

    v_y_re = 5.0 * x_re;
    v_y_im = 5.0 * x_im;
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
    w_y_re = rt_powd_snf(t_s.re, 3.0);
    w_y_im = 0.0;
  } else if (t_s.re == 0.0) {
    w_y_re = 0.0;
    w_y_im = -rt_powd_snf(t_s.im, 3.0);
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

    w_y_re = 3.0 * x_re;
    w_y_im = 3.0 * x_im;
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
    x_y_re = rt_powd_snf(t_s.re, 3.0);
    x_y_im = 0.0;
  } else if (t_s.re == 0.0) {
    x_y_re = 0.0;
    x_y_im = -rt_powd_snf(t_s.im, 3.0);
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

    x_y_re = 3.0 * x_re;
    x_y_im = 3.0 * x_im;
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
    y_y_re = rt_powd_snf(t_s.re, 3.0);
    y_y_im = 0.0;
  } else if (t_s.re == 0.0) {
    y_y_re = 0.0;
    y_y_im = -rt_powd_snf(t_s.im, 3.0);
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

    y_y_re = 3.0 * x_re;
    y_y_im = 3.0 * x_im;
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
    cb_y_re = rt_powd_snf(t_s.re, 5.0);
    cb_y_im = 0.0;
  } else if (t_s.re == 0.0) {
    cb_y_re = 0.0;
    cb_y_im = rt_powd_snf(t_s.im, 5.0);
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

    cb_y_re = 5.0 * x_re;
    cb_y_im = 5.0 * x_im;
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
    gb_y_re = rt_powd_snf(t_s.re, 4.0);
    gb_y_im = 0.0;
  } else if (t_s.re == 0.0) {
    gb_y_re = rt_powd_snf(t_s.im, 4.0);
    gb_y_im = 0.0;
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

    gb_y_re = 4.0 * x_re;
    gb_y_im = 4.0 * x_im;
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
    hb_y_re = rt_powd_snf(t_s.re, 4.0);
    hb_y_im = 0.0;
  } else if (t_s.re == 0.0) {
    hb_y_re = rt_powd_snf(t_s.im, 4.0);
    hb_y_im = 0.0;
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

    hb_y_re = 4.0 * x_re;
    hb_y_im = 4.0 * x_im;
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
    lb_y_re = rt_powd_snf(t_s.re, 4.0);
    lb_y_im = 0.0;
  } else if (t_s.re == 0.0) {
    lb_y_re = rt_powd_snf(t_s.im, 4.0);
    lb_y_im = 0.0;
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

    lb_y_re = 4.0 * x_re;
    lb_y_im = 4.0 * x_im;
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
    nb_y_re = rt_powd_snf(t_s.re, 4.0);
    nb_y_im = 0.0;
  } else if (t_s.re == 0.0) {
    nb_y_re = rt_powd_snf(t_s.im, 4.0);
    nb_y_im = 0.0;
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

    nb_y_re = 4.0 * x_re;
    nb_y_im = 4.0 * x_im;
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
    pb_y_re = rt_powd_snf(t_s.re, 4.0);
    pb_y_im = 0.0;
  } else if (t_s.re == 0.0) {
    pb_y_re = rt_powd_snf(t_s.im, 4.0);
    pb_y_im = 0.0;
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

    pb_y_re = 4.0 * x_re;
    pb_y_im = 4.0 * x_im;
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

  ar = 18.0 * (x01 * x01);
  br = 125.0 * y_re;
  bi = 125.0 * y_im;
  if (bi == 0.0) {
    re = ar / br;
  } else if (br == 0.0) {
    if (ar == 0.0) {
      re = 0.0 / bi;
    } else {
      re = 0.0;
    }
  } else {
    brm = fabs(br);
    bim = fabs(bi);
    if (brm > bim) {
      brm = bi / br;
      re = (ar + brm * 0.0) / (br + brm * bi);
    } else if (bim == brm) {
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

      re = (ar * b_br + 0.0 * b_bi) / brm;
    } else {
      brm = br / bi;
      re = brm * ar / (bi + brm * br);
    }
  }

  ar = 24.0 * (x04 * x04);
  br = 625.0 * b_y_re;
  bi = 625.0 * b_y_im;
  if (bi == 0.0) {
    b_re = ar / br;
  } else if (br == 0.0) {
    if (ar == 0.0) {
      b_re = 0.0 / bi;
    } else {
      b_re = 0.0;
    }
  } else {
    brm = fabs(br);
    bim = fabs(bi);
    if (brm > bim) {
      brm = bi / br;
      b_re = (ar + brm * 0.0) / (br + brm * bi);
    } else if (bim == brm) {
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

      b_re = (ar * c_br + 0.0 * c_bi) / brm;
    } else {
      brm = br / bi;
      b_re = brm * ar / (bi + brm * br);
    }
  }

  ar = 18.0 * (x02 * x02);
  br = 125.0 * c_y_re;
  bi = 125.0 * c_y_im;
  if (bi == 0.0) {
    c_re = ar / br;
  } else if (br == 0.0) {
    if (ar == 0.0) {
      c_re = 0.0 / bi;
    } else {
      c_re = 0.0;
    }
  } else {
    brm = fabs(br);
    bim = fabs(bi);
    if (brm > bim) {
      brm = bi / br;
      c_re = (ar + brm * 0.0) / (br + brm * bi);
    } else if (bim == brm) {
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

      c_re = (ar * d_br + 0.0 * d_bi) / brm;
    } else {
      brm = br / bi;
      c_re = brm * ar / (bi + brm * br);
    }
  }

  ar = 9.0 * (x07 * x07);
  y_re = 5000.0 * t_s.re;
  y_im = 5000.0 * t_s.im;
  if (y_im == 0.0) {
    d_re = ar / y_re;
  } else if (y_re == 0.0) {
    if (ar == 0.0) {
      d_re = 0.0 / y_im;
    } else {
      d_re = 0.0;
    }
  } else {
    brm = fabs(y_re);
    bim = fabs(y_im);
    if (brm > bim) {
      brm = y_im / y_re;
      d_re = (ar + brm * 0.0) / (y_re + brm * y_im);
    } else if (bim == brm) {
      if (y_re > 0.0) {
        qb_y_re = 0.5;
      } else {
        qb_y_re = -0.5;
      }

      if (y_im > 0.0) {
        qb_y_im = 0.5;
      } else {
        qb_y_im = -0.5;
      }

      d_re = (ar * qb_y_re + 0.0 * qb_y_im) / brm;
    } else {
      brm = y_re / y_im;
      d_re = brm * ar / (y_im + brm * y_re);
    }
  }

  ar = 24.0 * (x05 * x05);
  br = 625.0 * d_y_re;
  bi = 625.0 * d_y_im;
  if (bi == 0.0) {
    e_re = ar / br;
  } else if (br == 0.0) {
    if (ar == 0.0) {
      e_re = 0.0 / bi;
    } else {
      e_re = 0.0;
    }
  } else {
    brm = fabs(br);
    bim = fabs(bi);
    if (brm > bim) {
      brm = bi / br;
      e_re = (ar + brm * 0.0) / (br + brm * bi);
    } else if (bim == brm) {
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

      e_re = (ar * e_br + 0.0 * e_bi) / brm;
    } else {
      brm = br / bi;
      e_re = brm * ar / (bi + brm * br);
    }
  }

  ar = 18.0 * (x03 * x03);
  br = 125.0 * e_y_re;
  bi = 125.0 * e_y_im;
  if (bi == 0.0) {
    f_re = ar / br;
  } else if (br == 0.0) {
    if (ar == 0.0) {
      f_re = 0.0 / bi;
    } else {
      f_re = 0.0;
    }
  } else {
    brm = fabs(br);
    bim = fabs(bi);
    if (brm > bim) {
      brm = bi / br;
      f_re = (ar + brm * 0.0) / (br + brm * bi);
    } else if (bim == brm) {
      if (br > 0.0) {
        f_br = 0.5;
      } else {
        f_br = -0.5;
      }

      if (bi > 0.0) {
        f_bi = 0.5;
      } else {
        f_bi = -0.5;
      }

      f_re = (ar * f_br + 0.0 * f_bi) / brm;
    } else {
      brm = br / bi;
      f_re = brm * ar / (bi + brm * br);
    }
  }

  ar = 9.0 * (x08 * x08);
  if (y_im == 0.0) {
    g_re = ar / y_re;
  } else if (y_re == 0.0) {
    if (ar == 0.0) {
      g_re = 0.0 / y_im;
    } else {
      g_re = 0.0;
    }
  } else {
    brm = fabs(y_re);
    bim = fabs(y_im);
    if (brm > bim) {
      brm = y_im / y_re;
      g_re = (ar + brm * 0.0) / (y_re + brm * y_im);
    } else if (bim == brm) {
      if (y_re > 0.0) {
        rb_y_re = 0.5;
      } else {
        rb_y_re = -0.5;
      }

      if (y_im > 0.0) {
        rb_y_im = 0.5;
      } else {
        rb_y_im = -0.5;
      }

      g_re = (ar * rb_y_re + 0.0 * rb_y_im) / brm;
    } else {
      brm = y_re / y_im;
      g_re = brm * ar / (y_im + brm * y_re);
    }
  }

  ar = 24.0 * (x06 * x06);
  br = 625.0 * f_y_re;
  bi = 625.0 * f_y_im;
  if (bi == 0.0) {
    h_re = ar / br;
  } else if (br == 0.0) {
    if (ar == 0.0) {
      h_re = 0.0 / bi;
    } else {
      h_re = 0.0;
    }
  } else {
    brm = fabs(br);
    bim = fabs(bi);
    if (brm > bim) {
      brm = bi / br;
      h_re = (ar + brm * 0.0) / (br + brm * bi);
    } else if (bim == brm) {
      if (br > 0.0) {
        g_br = 0.5;
      } else {
        g_br = -0.5;
      }

      if (bi > 0.0) {
        g_bi = 0.5;
      } else {
        g_bi = -0.5;
      }

      h_re = (ar * g_br + 0.0 * g_bi) / brm;
    } else {
      brm = br / bi;
      h_re = brm * ar / (bi + brm * br);
    }
  }

  ar = 9.0 * (x09 * x09);
  if (y_im == 0.0) {
    i_re = ar / y_re;
  } else if (y_re == 0.0) {
    if (ar == 0.0) {
      i_re = 0.0 / y_im;
    } else {
      i_re = 0.0;
    }
  } else {
    brm = fabs(y_re);
    bim = fabs(y_im);
    if (brm > bim) {
      brm = y_im / y_re;
      i_re = (ar + brm * 0.0) / (y_re + brm * y_im);
    } else if (bim == brm) {
      if (y_re > 0.0) {
        sb_y_re = 0.5;
      } else {
        sb_y_re = -0.5;
      }

      if (y_im > 0.0) {
        sb_y_im = 0.5;
      } else {
        sb_y_im = -0.5;
      }

      i_re = (ar * sb_y_re + 0.0 * sb_y_im) / brm;
    } else {
      brm = y_re / y_im;
      i_re = brm * ar / (y_im + brm * y_re);
    }
  }

  ar = 18.0 * (x11 * x11);
  br = 125.0 * g_y_re;
  bi = 125.0 * g_y_im;
  if (bi == 0.0) {
    j_re = ar / br;
  } else if (br == 0.0) {
    if (ar == 0.0) {
      j_re = 0.0 / bi;
    } else {
      j_re = 0.0;
    }
  } else {
    brm = fabs(br);
    bim = fabs(bi);
    if (brm > bim) {
      brm = bi / br;
      j_re = (ar + brm * 0.0) / (br + brm * bi);
    } else if (bim == brm) {
      if (br > 0.0) {
        h_br = 0.5;
      } else {
        h_br = -0.5;
      }

      if (bi > 0.0) {
        h_bi = 0.5;
      } else {
        h_bi = -0.5;
      }

      j_re = (ar * h_br + 0.0 * h_bi) / brm;
    } else {
      brm = br / bi;
      j_re = brm * ar / (bi + brm * br);
    }
  }

  ar = 24.0 * (x14 * x14);
  br = 625.0 * h_y_re;
  bi = 625.0 * h_y_im;
  if (bi == 0.0) {
    k_re = ar / br;
  } else if (br == 0.0) {
    if (ar == 0.0) {
      k_re = 0.0 / bi;
    } else {
      k_re = 0.0;
    }
  } else {
    brm = fabs(br);
    bim = fabs(bi);
    if (brm > bim) {
      brm = bi / br;
      k_re = (ar + brm * 0.0) / (br + brm * bi);
    } else if (bim == brm) {
      if (br > 0.0) {
        i_br = 0.5;
      } else {
        i_br = -0.5;
      }

      if (bi > 0.0) {
        i_bi = 0.5;
      } else {
        i_bi = -0.5;
      }

      k_re = (ar * i_br + 0.0 * i_bi) / brm;
    } else {
      brm = br / bi;
      k_re = brm * ar / (bi + brm * br);
    }
  }

  ar = 18.0 * (x12 * x12);
  br = 125.0 * i_y_re;
  bi = 125.0 * i_y_im;
  if (bi == 0.0) {
    l_re = ar / br;
  } else if (br == 0.0) {
    if (ar == 0.0) {
      l_re = 0.0 / bi;
    } else {
      l_re = 0.0;
    }
  } else {
    brm = fabs(br);
    bim = fabs(bi);
    if (brm > bim) {
      brm = bi / br;
      l_re = (ar + brm * 0.0) / (br + brm * bi);
    } else if (bim == brm) {
      if (br > 0.0) {
        j_br = 0.5;
      } else {
        j_br = -0.5;
      }

      if (bi > 0.0) {
        j_bi = 0.5;
      } else {
        j_bi = -0.5;
      }

      l_re = (ar * j_br + 0.0 * j_bi) / brm;
    } else {
      brm = br / bi;
      l_re = brm * ar / (bi + brm * br);
    }
  }

  ar = 9.0 * (x17 * x17);
  if (y_im == 0.0) {
    m_re = ar / y_re;
  } else if (y_re == 0.0) {
    if (ar == 0.0) {
      m_re = 0.0 / y_im;
    } else {
      m_re = 0.0;
    }
  } else {
    brm = fabs(y_re);
    bim = fabs(y_im);
    if (brm > bim) {
      brm = y_im / y_re;
      m_re = (ar + brm * 0.0) / (y_re + brm * y_im);
    } else if (bim == brm) {
      if (y_re > 0.0) {
        tb_y_re = 0.5;
      } else {
        tb_y_re = -0.5;
      }

      if (y_im > 0.0) {
        tb_y_im = 0.5;
      } else {
        tb_y_im = -0.5;
      }

      m_re = (ar * tb_y_re + 0.0 * tb_y_im) / brm;
    } else {
      brm = y_re / y_im;
      m_re = brm * ar / (y_im + brm * y_re);
    }
  }

  ar = 24.0 * (x15 * x15);
  br = 625.0 * j_y_re;
  bi = 625.0 * j_y_im;
  if (bi == 0.0) {
    n_re = ar / br;
  } else if (br == 0.0) {
    if (ar == 0.0) {
      n_re = 0.0 / bi;
    } else {
      n_re = 0.0;
    }
  } else {
    brm = fabs(br);
    bim = fabs(bi);
    if (brm > bim) {
      brm = bi / br;
      n_re = (ar + brm * 0.0) / (br + brm * bi);
    } else if (bim == brm) {
      if (br > 0.0) {
        k_br = 0.5;
      } else {
        k_br = -0.5;
      }

      if (bi > 0.0) {
        k_bi = 0.5;
      } else {
        k_bi = -0.5;
      }

      n_re = (ar * k_br + 0.0 * k_bi) / brm;
    } else {
      brm = br / bi;
      n_re = brm * ar / (bi + brm * br);
    }
  }

  ar = 18.0 * (x13 * x13);
  br = 125.0 * k_y_re;
  bi = 125.0 * k_y_im;
  if (bi == 0.0) {
    o_re = ar / br;
  } else if (br == 0.0) {
    if (ar == 0.0) {
      o_re = 0.0 / bi;
    } else {
      o_re = 0.0;
    }
  } else {
    brm = fabs(br);
    bim = fabs(bi);
    if (brm > bim) {
      brm = bi / br;
      o_re = (ar + brm * 0.0) / (br + brm * bi);
    } else if (bim == brm) {
      if (br > 0.0) {
        l_br = 0.5;
      } else {
        l_br = -0.5;
      }

      if (bi > 0.0) {
        l_bi = 0.5;
      } else {
        l_bi = -0.5;
      }

      o_re = (ar * l_br + 0.0 * l_bi) / brm;
    } else {
      brm = br / bi;
      o_re = brm * ar / (bi + brm * br);
    }
  }

  ar = 9.0 * (x18 * x18);
  if (y_im == 0.0) {
    p_re = ar / y_re;
  } else if (y_re == 0.0) {
    if (ar == 0.0) {
      p_re = 0.0 / y_im;
    } else {
      p_re = 0.0;
    }
  } else {
    brm = fabs(y_re);
    bim = fabs(y_im);
    if (brm > bim) {
      brm = y_im / y_re;
      p_re = (ar + brm * 0.0) / (y_re + brm * y_im);
    } else if (bim == brm) {
      if (y_re > 0.0) {
        ub_y_re = 0.5;
      } else {
        ub_y_re = -0.5;
      }

      if (y_im > 0.0) {
        ub_y_im = 0.5;
      } else {
        ub_y_im = -0.5;
      }

      p_re = (ar * ub_y_re + 0.0 * ub_y_im) / brm;
    } else {
      brm = y_re / y_im;
      p_re = brm * ar / (y_im + brm * y_re);
    }
  }

  ar = 24.0 * (x16 * x16);
  br = 625.0 * l_y_re;
  bi = 625.0 * l_y_im;
  if (bi == 0.0) {
    q_re = ar / br;
  } else if (br == 0.0) {
    if (ar == 0.0) {
      q_re = 0.0 / bi;
    } else {
      q_re = 0.0;
    }
  } else {
    brm = fabs(br);
    bim = fabs(bi);
    if (brm > bim) {
      brm = bi / br;
      q_re = (ar + brm * 0.0) / (br + brm * bi);
    } else if (bim == brm) {
      if (br > 0.0) {
        m_br = 0.5;
      } else {
        m_br = -0.5;
      }

      if (bi > 0.0) {
        m_bi = 0.5;
      } else {
        m_bi = -0.5;
      }

      q_re = (ar * m_br + 0.0 * m_bi) / brm;
    } else {
      brm = br / bi;
      q_re = brm * ar / (bi + brm * br);
    }
  }

  ar = 9.0 * (x19 * x19);
  if (y_im == 0.0) {
    r_re = ar / y_re;
  } else if (y_re == 0.0) {
    if (ar == 0.0) {
      r_re = 0.0 / y_im;
    } else {
      r_re = 0.0;
    }
  } else {
    brm = fabs(y_re);
    bim = fabs(y_im);
    if (brm > bim) {
      brm = y_im / y_re;
      r_re = (ar + brm * 0.0) / (y_re + brm * y_im);
    } else if (bim == brm) {
      if (y_re > 0.0) {
        vb_y_re = 0.5;
      } else {
        vb_y_re = -0.5;
      }

      if (y_im > 0.0) {
        vb_y_im = 0.5;
      } else {
        vb_y_im = -0.5;
      }

      r_re = (ar * vb_y_re + 0.0 * vb_y_im) / brm;
    } else {
      brm = y_re / y_im;
      r_re = brm * ar / (y_im + brm * y_re);
    }
  }

  ar = 18.0 * x01 * x04;
  br = 125.0 * m_y_re;
  bi = 125.0 * m_y_im;
  if (bi == 0.0) {
    s_re = ar / br;
  } else if (br == 0.0) {
    if (ar == 0.0) {
      s_re = 0.0 / bi;
    } else {
      s_re = 0.0;
    }
  } else {
    brm = fabs(br);
    bim = fabs(bi);
    if (brm > bim) {
      brm = bi / br;
      s_re = (ar + brm * 0.0) / (br + brm * bi);
    } else if (bim == brm) {
      if (br > 0.0) {
        n_br = 0.5;
      } else {
        n_br = -0.5;
      }

      if (bi > 0.0) {
        n_bi = 0.5;
      } else {
        n_bi = -0.5;
      }

      s_re = (ar * n_br + 0.0 * n_bi) / brm;
    } else {
      brm = br / bi;
      s_re = brm * ar / (bi + brm * br);
    }
  }

  ar = 3.0 * x01 * x07;
  br = 125.0 * n_y_re;
  bi = 125.0 * n_y_im;
  if (bi == 0.0) {
    t_re = ar / br;
  } else if (br == 0.0) {
    if (ar == 0.0) {
      t_re = 0.0 / bi;
    } else {
      t_re = 0.0;
    }
  } else {
    brm = fabs(br);
    bim = fabs(bi);
    if (brm > bim) {
      brm = bi / br;
      t_re = (ar + brm * 0.0) / (br + brm * bi);
    } else if (bim == brm) {
      if (br > 0.0) {
        o_br = 0.5;
      } else {
        o_br = -0.5;
      }

      if (bi > 0.0) {
        o_bi = 0.5;
      } else {
        o_bi = -0.5;
      }

      t_re = (ar * o_br + 0.0 * o_bi) / brm;
    } else {
      brm = br / bi;
      t_re = brm * ar / (bi + brm * br);
    }
  }

  ar = 18.0 * x02 * x05;
  br = 125.0 * o_y_re;
  bi = 125.0 * o_y_im;
  if (bi == 0.0) {
    u_re = ar / br;
  } else if (br == 0.0) {
    if (ar == 0.0) {
      u_re = 0.0 / bi;
    } else {
      u_re = 0.0;
    }
  } else {
    brm = fabs(br);
    bim = fabs(bi);
    if (brm > bim) {
      brm = bi / br;
      u_re = (ar + brm * 0.0) / (br + brm * bi);
    } else if (bim == brm) {
      if (br > 0.0) {
        p_br = 0.5;
      } else {
        p_br = -0.5;
      }

      if (bi > 0.0) {
        p_bi = 0.5;
      } else {
        p_bi = -0.5;
      }

      u_re = (ar * p_br + 0.0 * p_bi) / brm;
    } else {
      brm = br / bi;
      u_re = brm * ar / (bi + brm * br);
    }
  }

  r = t_s.re * t_s.re - t_s.im * t_s.im;
  t_s_im_tmp = t_s.re * t_s.im + t_s.im * t_s.re;
  ar = 9.0 * x04 * x07;
  br = 625.0 * r;
  bi = 625.0 * t_s_im_tmp;
  if (bi == 0.0) {
    v_re = ar / br;
  } else if (br == 0.0) {
    if (ar == 0.0) {
      v_re = 0.0 / bi;
    } else {
      v_re = 0.0;
    }
  } else {
    brm = fabs(br);
    bim = fabs(bi);
    if (brm > bim) {
      brm = bi / br;
      v_re = (ar + brm * 0.0) / (br + brm * bi);
    } else if (bim == brm) {
      if (br > 0.0) {
        q_br = 0.5;
      } else {
        q_br = -0.5;
      }

      if (bi > 0.0) {
        q_bi = 0.5;
      } else {
        q_bi = -0.5;
      }

      v_re = (ar * q_br + 0.0 * q_bi) / brm;
    } else {
      brm = br / bi;
      v_re = brm * ar / (bi + brm * br);
    }
  }

  ar = 3.0 * x02 * x08;
  br = 125.0 * p_y_re;
  bi = 125.0 * p_y_im;
  if (bi == 0.0) {
    w_re = ar / br;
  } else if (br == 0.0) {
    if (ar == 0.0) {
      w_re = 0.0 / bi;
    } else {
      w_re = 0.0;
    }
  } else {
    brm = fabs(br);
    bim = fabs(bi);
    if (brm > bim) {
      brm = bi / br;
      w_re = (ar + brm * 0.0) / (br + brm * bi);
    } else if (bim == brm) {
      if (br > 0.0) {
        r_br = 0.5;
      } else {
        r_br = -0.5;
      }

      if (bi > 0.0) {
        r_bi = 0.5;
      } else {
        r_bi = -0.5;
      }

      w_re = (ar * r_br + 0.0 * r_bi) / brm;
    } else {
      brm = br / bi;
      w_re = brm * ar / (bi + brm * br);
    }
  }

  ar = 18.0 * x03 * x06;
  br = 125.0 * q_y_re;
  bi = 125.0 * q_y_im;
  if (bi == 0.0) {
    ab_re = ar / br;
  } else if (br == 0.0) {
    if (ar == 0.0) {
      ab_re = 0.0 / bi;
    } else {
      ab_re = 0.0;
    }
  } else {
    brm = fabs(br);
    bim = fabs(bi);
    if (brm > bim) {
      brm = bi / br;
      ab_re = (ar + brm * 0.0) / (br + brm * bi);
    } else if (bim == brm) {
      if (br > 0.0) {
        s_br = 0.5;
      } else {
        s_br = -0.5;
      }

      if (bi > 0.0) {
        s_bi = 0.5;
      } else {
        s_bi = -0.5;
      }

      ab_re = (ar * s_br + 0.0 * s_bi) / brm;
    } else {
      brm = br / bi;
      ab_re = brm * ar / (bi + brm * br);
    }
  }

  ar = 9.0 * x05 * x08;
  br = 625.0 * r;
  bi = 625.0 * t_s_im_tmp;
  if (bi == 0.0) {
    bb_re = ar / br;
  } else if (br == 0.0) {
    if (ar == 0.0) {
      bb_re = 0.0 / bi;
    } else {
      bb_re = 0.0;
    }
  } else {
    brm = fabs(br);
    bim = fabs(bi);
    if (brm > bim) {
      brm = bi / br;
      bb_re = (ar + brm * 0.0) / (br + brm * bi);
    } else if (bim == brm) {
      if (br > 0.0) {
        t_br = 0.5;
      } else {
        t_br = -0.5;
      }

      if (bi > 0.0) {
        t_bi = 0.5;
      } else {
        t_bi = -0.5;
      }

      bb_re = (ar * t_br + 0.0 * t_bi) / brm;
    } else {
      brm = br / bi;
      bb_re = brm * ar / (bi + brm * br);
    }
  }

  ar = 3.0 * x03 * x09;
  br = 125.0 * r_y_re;
  bi = 125.0 * r_y_im;
  if (bi == 0.0) {
    r_y_im = ar / br;
  } else if (br == 0.0) {
    if (ar == 0.0) {
      r_y_im = 0.0 / bi;
    } else {
      r_y_im = 0.0;
    }
  } else {
    brm = fabs(br);
    bim = fabs(bi);
    if (brm > bim) {
      brm = bi / br;
      r_y_im = (ar + brm * 0.0) / (br + brm * bi);
    } else if (bim == brm) {
      if (br > 0.0) {
        u_br = 0.5;
      } else {
        u_br = -0.5;
      }

      if (bi > 0.0) {
        u_bi = 0.5;
      } else {
        u_bi = -0.5;
      }

      r_y_im = (ar * u_br + 0.0 * u_bi) / brm;
    } else {
      brm = br / bi;
      r_y_im = brm * ar / (bi + brm * br);
    }
  }

  ar = 9.0 * x06 * x09;
  br = 625.0 * r;
  bi = 625.0 * t_s_im_tmp;
  if (bi == 0.0) {
    cb_re = ar / br;
  } else if (br == 0.0) {
    if (ar == 0.0) {
      cb_re = 0.0 / bi;
    } else {
      cb_re = 0.0;
    }
  } else {
    brm = fabs(br);
    bim = fabs(bi);
    if (brm > bim) {
      brm = bi / br;
      cb_re = (ar + brm * 0.0) / (br + brm * bi);
    } else if (bim == brm) {
      if (br > 0.0) {
        v_br = 0.5;
      } else {
        v_br = -0.5;
      }

      if (bi > 0.0) {
        v_bi = 0.5;
      } else {
        v_bi = -0.5;
      }

      cb_re = (ar * v_br + 0.0 * v_bi) / brm;
    } else {
      brm = br / bi;
      cb_re = brm * ar / (bi + brm * br);
    }
  }

  ar = 36.0 * x01 * x11;
  br = 125.0 * s_y_re;
  bi = 125.0 * s_y_im;
  if (bi == 0.0) {
    r_y_re = ar / br;
  } else if (br == 0.0) {
    if (ar == 0.0) {
      r_y_re = 0.0 / bi;
    } else {
      r_y_re = 0.0;
    }
  } else {
    brm = fabs(br);
    bim = fabs(bi);
    if (brm > bim) {
      brm = bi / br;
      r_y_re = (ar + brm * 0.0) / (br + brm * bi);
    } else if (bim == brm) {
      if (br > 0.0) {
        w_br = 0.5;
      } else {
        w_br = -0.5;
      }

      if (bi > 0.0) {
        w_bi = 0.5;
      } else {
        w_bi = -0.5;
      }

      r_y_re = (ar * w_br + 0.0 * w_bi) / brm;
    } else {
      brm = br / bi;
      r_y_re = brm * ar / (bi + brm * br);
    }
  }

  ar = 18.0 * x01 * x14;
  br = 125.0 * t_y_re;
  bi = 125.0 * t_y_im;
  if (bi == 0.0) {
    q_y_im = ar / br;
  } else if (br == 0.0) {
    if (ar == 0.0) {
      q_y_im = 0.0 / bi;
    } else {
      q_y_im = 0.0;
    }
  } else {
    brm = fabs(br);
    bim = fabs(bi);
    if (brm > bim) {
      brm = bi / br;
      q_y_im = (ar + brm * 0.0) / (br + brm * bi);
    } else if (bim == brm) {
      if (br > 0.0) {
        x_br = 0.5;
      } else {
        x_br = -0.5;
      }

      if (bi > 0.0) {
        x_bi = 0.5;
      } else {
        x_bi = -0.5;
      }

      q_y_im = (ar * x_br + 0.0 * x_bi) / brm;
    } else {
      brm = br / bi;
      q_y_im = brm * ar / (bi + brm * br);
    }
  }

  ar = 18.0 * x04 * x11;
  br = 125.0 * u_y_re;
  bi = 125.0 * u_y_im;
  if (bi == 0.0) {
    q_y_re = ar / br;
  } else if (br == 0.0) {
    if (ar == 0.0) {
      q_y_re = 0.0 / bi;
    } else {
      q_y_re = 0.0;
    }
  } else {
    brm = fabs(br);
    bim = fabs(bi);
    if (brm > bim) {
      brm = bi / br;
      q_y_re = (ar + brm * 0.0) / (br + brm * bi);
    } else if (bim == brm) {
      if (br > 0.0) {
        y_br = 0.5;
      } else {
        y_br = -0.5;
      }

      if (bi > 0.0) {
        y_bi = 0.5;
      } else {
        y_bi = -0.5;
      }

      q_y_re = (ar * y_br + 0.0 * y_bi) / brm;
    } else {
      brm = br / bi;
      q_y_re = brm * ar / (bi + brm * br);
    }
  }

  ar = 36.0 * x02 * x12;
  br = 125.0 * v_y_re;
  bi = 125.0 * v_y_im;
  if (bi == 0.0) {
    p_y_im = ar / br;
  } else if (br == 0.0) {
    if (ar == 0.0) {
      p_y_im = 0.0 / bi;
    } else {
      p_y_im = 0.0;
    }
  } else {
    brm = fabs(br);
    bim = fabs(bi);
    if (brm > bim) {
      brm = bi / br;
      p_y_im = (ar + brm * 0.0) / (br + brm * bi);
    } else if (bim == brm) {
      if (br > 0.0) {
        ab_br = 0.5;
      } else {
        ab_br = -0.5;
      }

      if (bi > 0.0) {
        ab_bi = 0.5;
      } else {
        ab_bi = -0.5;
      }

      p_y_im = (ar * ab_br + 0.0 * ab_bi) / brm;
    } else {
      brm = br / bi;
      p_y_im = brm * ar / (bi + brm * br);
    }
  }

  ar = 3.0 * x01 * x17;
  br = 125.0 * w_y_re;
  bi = 125.0 * w_y_im;
  if (bi == 0.0) {
    p_y_re = ar / br;
  } else if (br == 0.0) {
    if (ar == 0.0) {
      p_y_re = 0.0 / bi;
    } else {
      p_y_re = 0.0;
    }
  } else {
    brm = fabs(br);
    bim = fabs(bi);
    if (brm > bim) {
      brm = bi / br;
      p_y_re = (ar + brm * 0.0) / (br + brm * bi);
    } else if (bim == brm) {
      if (br > 0.0) {
        bb_br = 0.5;
      } else {
        bb_br = -0.5;
      }

      if (bi > 0.0) {
        bb_bi = 0.5;
      } else {
        bb_bi = -0.5;
      }

      p_y_re = (ar * bb_br + 0.0 * bb_bi) / brm;
    } else {
      brm = br / bi;
      p_y_re = brm * ar / (bi + brm * br);
    }
  }

  ar = 42.0 * x04 * x14;
  br = 625.0 * x_y_re;
  bi = 625.0 * x_y_im;
  if (bi == 0.0) {
    o_y_im = ar / br;
  } else if (br == 0.0) {
    if (ar == 0.0) {
      o_y_im = 0.0 / bi;
    } else {
      o_y_im = 0.0;
    }
  } else {
    brm = fabs(br);
    bim = fabs(bi);
    if (brm > bim) {
      brm = bi / br;
      o_y_im = (ar + brm * 0.0) / (br + brm * bi);
    } else if (bim == brm) {
      if (br > 0.0) {
        cb_br = 0.5;
      } else {
        cb_br = -0.5;
      }

      if (bi > 0.0) {
        cb_bi = 0.5;
      } else {
        cb_bi = -0.5;
      }

      o_y_im = (ar * cb_br + 0.0 * cb_bi) / brm;
    } else {
      brm = br / bi;
      o_y_im = brm * ar / (bi + brm * br);
    }
  }

  ar = 3.0 * x07 * x11;
  br = 125.0 * y_y_re;
  bi = 125.0 * y_y_im;
  if (bi == 0.0) {
    o_y_re = ar / br;
  } else if (br == 0.0) {
    if (ar == 0.0) {
      o_y_re = 0.0 / bi;
    } else {
      o_y_re = 0.0;
    }
  } else {
    brm = fabs(br);
    bim = fabs(bi);
    if (brm > bim) {
      brm = bi / br;
      o_y_re = (ar + brm * 0.0) / (br + brm * bi);
    } else if (bim == brm) {
      if (br > 0.0) {
        db_br = 0.5;
      } else {
        db_br = -0.5;
      }

      if (bi > 0.0) {
        db_bi = 0.5;
      } else {
        db_bi = -0.5;
      }

      o_y_re = (ar * db_br + 0.0 * db_bi) / brm;
    } else {
      brm = br / bi;
      o_y_re = brm * ar / (bi + brm * br);
    }
  }

  ar = 18.0 * x02 * x15;
  br = 125.0 * ab_y_re;
  bi = 125.0 * ab_y_im;
  if (bi == 0.0) {
    n_y_im = ar / br;
  } else if (br == 0.0) {
    if (ar == 0.0) {
      n_y_im = 0.0 / bi;
    } else {
      n_y_im = 0.0;
    }
  } else {
    brm = fabs(br);
    bim = fabs(bi);
    if (brm > bim) {
      brm = bi / br;
      n_y_im = (ar + brm * 0.0) / (br + brm * bi);
    } else if (bim == brm) {
      if (br > 0.0) {
        eb_br = 0.5;
      } else {
        eb_br = -0.5;
      }

      if (bi > 0.0) {
        eb_bi = 0.5;
      } else {
        eb_bi = -0.5;
      }

      n_y_im = (ar * eb_br + 0.0 * eb_bi) / brm;
    } else {
      brm = br / bi;
      n_y_im = brm * ar / (bi + brm * br);
    }
  }

  ar = 18.0 * x05 * x12;
  br = 125.0 * bb_y_re;
  bi = 125.0 * bb_y_im;
  if (bi == 0.0) {
    n_y_re = ar / br;
  } else if (br == 0.0) {
    if (ar == 0.0) {
      n_y_re = 0.0 / bi;
    } else {
      n_y_re = 0.0;
    }
  } else {
    brm = fabs(br);
    bim = fabs(bi);
    if (brm > bim) {
      brm = bi / br;
      n_y_re = (ar + brm * 0.0) / (br + brm * bi);
    } else if (bim == brm) {
      if (br > 0.0) {
        fb_br = 0.5;
      } else {
        fb_br = -0.5;
      }

      if (bi > 0.0) {
        fb_bi = 0.5;
      } else {
        fb_bi = -0.5;
      }

      n_y_re = (ar * fb_br + 0.0 * fb_bi) / brm;
    } else {
      brm = br / bi;
      n_y_re = brm * ar / (bi + brm * br);
    }
  }

  ar = 36.0 * x03 * x13;
  br = 125.0 * cb_y_re;
  bi = 125.0 * cb_y_im;
  if (bi == 0.0) {
    l_y_im = ar / br;
  } else if (br == 0.0) {
    if (ar == 0.0) {
      l_y_im = 0.0 / bi;
    } else {
      l_y_im = 0.0;
    }
  } else {
    brm = fabs(br);
    bim = fabs(bi);
    if (brm > bim) {
      brm = bi / br;
      l_y_im = (ar + brm * 0.0) / (br + brm * bi);
    } else if (bim == brm) {
      if (br > 0.0) {
        gb_br = 0.5;
      } else {
        gb_br = -0.5;
      }

      if (bi > 0.0) {
        gb_bi = 0.5;
      } else {
        gb_bi = -0.5;
      }

      l_y_im = (ar * gb_br + 0.0 * gb_bi) / brm;
    } else {
      brm = br / bi;
      l_y_im = brm * ar / (bi + brm * br);
    }
  }

  ar = 6.0 * x04 * x17;
  br = 625.0 * r;
  bi = 625.0 * t_s_im_tmp;
  if (bi == 0.0) {
    m_y_re = ar / br;
  } else if (br == 0.0) {
    if (ar == 0.0) {
      m_y_re = 0.0 / bi;
    } else {
      m_y_re = 0.0;
    }
  } else {
    brm = fabs(br);
    bim = fabs(bi);
    if (brm > bim) {
      brm = bi / br;
      m_y_re = (ar + brm * 0.0) / (br + brm * bi);
    } else if (bim == brm) {
      if (br > 0.0) {
        hb_br = 0.5;
      } else {
        hb_br = -0.5;
      }

      if (bi > 0.0) {
        hb_bi = 0.5;
      } else {
        hb_bi = -0.5;
      }

      m_y_re = (ar * hb_br + 0.0 * hb_bi) / brm;
    } else {
      brm = br / bi;
      m_y_re = brm * ar / (bi + brm * br);
    }
  }

  ar = 6.0 * x07 * x14;
  br = 625.0 * r;
  bi = 625.0 * t_s_im_tmp;
  if (bi == 0.0) {
    m_y_im = ar / br;
  } else if (br == 0.0) {
    if (ar == 0.0) {
      m_y_im = 0.0 / bi;
    } else {
      m_y_im = 0.0;
    }
  } else {
    brm = fabs(br);
    bim = fabs(bi);
    if (brm > bim) {
      brm = bi / br;
      m_y_im = (ar + brm * 0.0) / (br + brm * bi);
    } else if (bim == brm) {
      if (br > 0.0) {
        ib_br = 0.5;
      } else {
        ib_br = -0.5;
      }

      if (bi > 0.0) {
        ib_bi = 0.5;
      } else {
        ib_bi = -0.5;
      }

      m_y_im = (ar * ib_br + 0.0 * ib_bi) / brm;
    } else {
      brm = br / bi;
      m_y_im = brm * ar / (bi + brm * br);
    }
  }

  ar = 3.0 * x02 * x18;
  br = 125.0 * db_y_re;
  bi = 125.0 * db_y_im;
  if (bi == 0.0) {
    l_y_re = ar / br;
  } else if (br == 0.0) {
    if (ar == 0.0) {
      l_y_re = 0.0 / bi;
    } else {
      l_y_re = 0.0;
    }
  } else {
    brm = fabs(br);
    bim = fabs(bi);
    if (brm > bim) {
      brm = bi / br;
      l_y_re = (ar + brm * 0.0) / (br + brm * bi);
    } else if (bim == brm) {
      if (br > 0.0) {
        jb_br = 0.5;
      } else {
        jb_br = -0.5;
      }

      if (bi > 0.0) {
        jb_bi = 0.5;
      } else {
        jb_bi = -0.5;
      }

      l_y_re = (ar * jb_br + 0.0 * jb_bi) / brm;
    } else {
      brm = br / bi;
      l_y_re = brm * ar / (bi + brm * br);
    }
  }

  ar = 42.0 * x05 * x15;
  br = 625.0 * eb_y_re;
  bi = 625.0 * eb_y_im;
  if (bi == 0.0) {
    k_y_im = ar / br;
  } else if (br == 0.0) {
    if (ar == 0.0) {
      k_y_im = 0.0 / bi;
    } else {
      k_y_im = 0.0;
    }
  } else {
    brm = fabs(br);
    bim = fabs(bi);
    if (brm > bim) {
      brm = bi / br;
      k_y_im = (ar + brm * 0.0) / (br + brm * bi);
    } else if (bim == brm) {
      if (br > 0.0) {
        kb_br = 0.5;
      } else {
        kb_br = -0.5;
      }

      if (bi > 0.0) {
        kb_bi = 0.5;
      } else {
        kb_bi = -0.5;
      }

      k_y_im = (ar * kb_br + 0.0 * kb_bi) / brm;
    } else {
      brm = br / bi;
      k_y_im = brm * ar / (bi + brm * br);
    }
  }

  ar = 3.0 * x08 * x12;
  br = 125.0 * fb_y_re;
  bi = 125.0 * fb_y_im;
  if (bi == 0.0) {
    k_y_re = ar / br;
  } else if (br == 0.0) {
    if (ar == 0.0) {
      k_y_re = 0.0 / bi;
    } else {
      k_y_re = 0.0;
    }
  } else {
    brm = fabs(br);
    bim = fabs(bi);
    if (brm > bim) {
      brm = bi / br;
      k_y_re = (ar + brm * 0.0) / (br + brm * bi);
    } else if (bim == brm) {
      if (br > 0.0) {
        lb_br = 0.5;
      } else {
        lb_br = -0.5;
      }

      if (bi > 0.0) {
        lb_bi = 0.5;
      } else {
        lb_bi = -0.5;
      }

      k_y_re = (ar * lb_br + 0.0 * lb_bi) / brm;
    } else {
      brm = br / bi;
      k_y_re = brm * ar / (bi + brm * br);
    }
  }

  ar = 18.0 * x03 * x16;
  br = 125.0 * gb_y_re;
  bi = 125.0 * gb_y_im;
  if (bi == 0.0) {
    j_y_im = ar / br;
  } else if (br == 0.0) {
    if (ar == 0.0) {
      j_y_im = 0.0 / bi;
    } else {
      j_y_im = 0.0;
    }
  } else {
    brm = fabs(br);
    bim = fabs(bi);
    if (brm > bim) {
      brm = bi / br;
      j_y_im = (ar + brm * 0.0) / (br + brm * bi);
    } else if (bim == brm) {
      if (br > 0.0) {
        mb_br = 0.5;
      } else {
        mb_br = -0.5;
      }

      if (bi > 0.0) {
        mb_bi = 0.5;
      } else {
        mb_bi = -0.5;
      }

      j_y_im = (ar * mb_br + 0.0 * mb_bi) / brm;
    } else {
      brm = br / bi;
      j_y_im = brm * ar / (bi + brm * br);
    }
  }

  ar = 18.0 * x06 * x13;
  br = 125.0 * hb_y_re;
  bi = 125.0 * hb_y_im;
  if (bi == 0.0) {
    h_y_im = ar / br;
  } else if (br == 0.0) {
    if (ar == 0.0) {
      h_y_im = 0.0 / bi;
    } else {
      h_y_im = 0.0;
    }
  } else {
    brm = fabs(br);
    bim = fabs(bi);
    if (brm > bim) {
      brm = bi / br;
      h_y_im = (ar + brm * 0.0) / (br + brm * bi);
    } else if (bim == brm) {
      if (br > 0.0) {
        nb_br = 0.5;
      } else {
        nb_br = -0.5;
      }

      if (bi > 0.0) {
        nb_bi = 0.5;
      } else {
        nb_bi = -0.5;
      }

      h_y_im = (ar * nb_br + 0.0 * nb_bi) / brm;
    } else {
      brm = br / bi;
      h_y_im = brm * ar / (bi + brm * br);
    }
  }

  ar = 3.0 * x07 * x17;
  y_re = 2500.0 * t_s.re;
  y_im = 2500.0 * t_s.im;
  if (y_im == 0.0) {
    i_y_re = ar / y_re;
  } else if (y_re == 0.0) {
    if (ar == 0.0) {
      i_y_re = 0.0 / y_im;
    } else {
      i_y_re = 0.0;
    }
  } else {
    brm = fabs(y_re);
    bim = fabs(y_im);
    if (brm > bim) {
      brm = y_im / y_re;
      i_y_re = (ar + brm * 0.0) / (y_re + brm * y_im);
    } else if (bim == brm) {
      if (y_re > 0.0) {
        xb_y_re = 0.5;
      } else {
        xb_y_re = -0.5;
      }

      if (y_im > 0.0) {
        wb_y_im = 0.5;
      } else {
        wb_y_im = -0.5;
      }

      i_y_re = (ar * xb_y_re + 0.0 * wb_y_im) / brm;
    } else {
      brm = y_re / y_im;
      i_y_re = brm * ar / (y_im + brm * y_re);
    }
  }

  ar = 6.0 * x05 * x18;
  br = 625.0 * r;
  bi = 625.0 * t_s_im_tmp;
  if (bi == 0.0) {
    i_y_im = ar / br;
  } else if (br == 0.0) {
    if (ar == 0.0) {
      i_y_im = 0.0 / bi;
    } else {
      i_y_im = 0.0;
    }
  } else {
    brm = fabs(br);
    bim = fabs(bi);
    if (brm > bim) {
      brm = bi / br;
      i_y_im = (ar + brm * 0.0) / (br + brm * bi);
    } else if (bim == brm) {
      if (br > 0.0) {
        ob_br = 0.5;
      } else {
        ob_br = -0.5;
      }

      if (bi > 0.0) {
        ob_bi = 0.5;
      } else {
        ob_bi = -0.5;
      }

      i_y_im = (ar * ob_br + 0.0 * ob_bi) / brm;
    } else {
      brm = br / bi;
      i_y_im = brm * ar / (bi + brm * br);
    }
  }

  ar = 6.0 * x08 * x15;
  br = 625.0 * r;
  bi = 625.0 * t_s_im_tmp;
  if (bi == 0.0) {
    j_y_re = ar / br;
  } else if (br == 0.0) {
    if (ar == 0.0) {
      j_y_re = 0.0 / bi;
    } else {
      j_y_re = 0.0;
    }
  } else {
    brm = fabs(br);
    bim = fabs(bi);
    if (brm > bim) {
      brm = bi / br;
      j_y_re = (ar + brm * 0.0) / (br + brm * bi);
    } else if (bim == brm) {
      if (br > 0.0) {
        pb_br = 0.5;
      } else {
        pb_br = -0.5;
      }

      if (bi > 0.0) {
        pb_bi = 0.5;
      } else {
        pb_bi = -0.5;
      }

      j_y_re = (ar * pb_br + 0.0 * pb_bi) / brm;
    } else {
      brm = br / bi;
      j_y_re = brm * ar / (bi + brm * br);
    }
  }

  ar = 3.0 * x03 * x19;
  br = 125.0 * ib_y_re;
  bi = 125.0 * ib_y_im;
  if (bi == 0.0) {
    h_y_re = ar / br;
  } else if (br == 0.0) {
    if (ar == 0.0) {
      h_y_re = 0.0 / bi;
    } else {
      h_y_re = 0.0;
    }
  } else {
    brm = fabs(br);
    bim = fabs(bi);
    if (brm > bim) {
      brm = bi / br;
      h_y_re = (ar + brm * 0.0) / (br + brm * bi);
    } else if (bim == brm) {
      if (br > 0.0) {
        qb_br = 0.5;
      } else {
        qb_br = -0.5;
      }

      if (bi > 0.0) {
        qb_bi = 0.5;
      } else {
        qb_bi = -0.5;
      }

      h_y_re = (ar * qb_br + 0.0 * qb_bi) / brm;
    } else {
      brm = br / bi;
      h_y_re = brm * ar / (bi + brm * br);
    }
  }

  ar = 42.0 * x06 * x16;
  br = 625.0 * jb_y_re;
  bi = 625.0 * jb_y_im;
  if (bi == 0.0) {
    g_y_im = ar / br;
  } else if (br == 0.0) {
    if (ar == 0.0) {
      g_y_im = 0.0 / bi;
    } else {
      g_y_im = 0.0;
    }
  } else {
    brm = fabs(br);
    bim = fabs(bi);
    if (brm > bim) {
      brm = bi / br;
      g_y_im = (ar + brm * 0.0) / (br + brm * bi);
    } else if (bim == brm) {
      if (br > 0.0) {
        rb_br = 0.5;
      } else {
        rb_br = -0.5;
      }

      if (bi > 0.0) {
        rb_bi = 0.5;
      } else {
        rb_bi = -0.5;
      }

      g_y_im = (ar * rb_br + 0.0 * rb_bi) / brm;
    } else {
      brm = br / bi;
      g_y_im = brm * ar / (bi + brm * br);
    }
  }

  ar = 3.0 * x09 * x13;
  br = 125.0 * kb_y_re;
  bi = 125.0 * kb_y_im;
  if (bi == 0.0) {
    e_y_im = ar / br;
  } else if (br == 0.0) {
    if (ar == 0.0) {
      e_y_im = 0.0 / bi;
    } else {
      e_y_im = 0.0;
    }
  } else {
    brm = fabs(br);
    bim = fabs(bi);
    if (brm > bim) {
      brm = bi / br;
      e_y_im = (ar + brm * 0.0) / (br + brm * bi);
    } else if (bim == brm) {
      if (br > 0.0) {
        sb_br = 0.5;
      } else {
        sb_br = -0.5;
      }

      if (bi > 0.0) {
        sb_bi = 0.5;
      } else {
        sb_bi = -0.5;
      }

      e_y_im = (ar * sb_br + 0.0 * sb_bi) / brm;
    } else {
      brm = br / bi;
      e_y_im = brm * ar / (bi + brm * br);
    }
  }

  ar = 3.0 * x08 * x18;
  if (y_im == 0.0) {
    f_y_re = ar / y_re;
  } else if (y_re == 0.0) {
    if (ar == 0.0) {
      f_y_re = 0.0 / y_im;
    } else {
      f_y_re = 0.0;
    }
  } else {
    brm = fabs(y_re);
    bim = fabs(y_im);
    if (brm > bim) {
      brm = y_im / y_re;
      f_y_re = (ar + brm * 0.0) / (y_re + brm * y_im);
    } else if (bim == brm) {
      if (y_re > 0.0) {
        yb_y_re = 0.5;
      } else {
        yb_y_re = -0.5;
      }

      if (y_im > 0.0) {
        xb_y_im = 0.5;
      } else {
        xb_y_im = -0.5;
      }

      f_y_re = (ar * yb_y_re + 0.0 * xb_y_im) / brm;
    } else {
      brm = y_re / y_im;
      f_y_re = brm * ar / (y_im + brm * y_re);
    }
  }

  ar = 6.0 * x06 * x19;
  br = 625.0 * r;
  bi = 625.0 * t_s_im_tmp;
  if (bi == 0.0) {
    f_y_im = ar / br;
  } else if (br == 0.0) {
    if (ar == 0.0) {
      f_y_im = 0.0 / bi;
    } else {
      f_y_im = 0.0;
    }
  } else {
    brm = fabs(br);
    bim = fabs(bi);
    if (brm > bim) {
      brm = bi / br;
      f_y_im = (ar + brm * 0.0) / (br + brm * bi);
    } else if (bim == brm) {
      if (br > 0.0) {
        tb_br = 0.5;
      } else {
        tb_br = -0.5;
      }

      if (bi > 0.0) {
        tb_bi = 0.5;
      } else {
        tb_bi = -0.5;
      }

      f_y_im = (ar * tb_br + 0.0 * tb_bi) / brm;
    } else {
      brm = br / bi;
      f_y_im = brm * ar / (bi + brm * br);
    }
  }

  ar = 6.0 * x09 * x16;
  br = 625.0 * r;
  bi = 625.0 * t_s_im_tmp;
  if (bi == 0.0) {
    g_y_re = ar / br;
  } else if (br == 0.0) {
    if (ar == 0.0) {
      g_y_re = 0.0 / bi;
    } else {
      g_y_re = 0.0;
    }
  } else {
    brm = fabs(br);
    bim = fabs(bi);
    if (brm > bim) {
      brm = bi / br;
      g_y_re = (ar + brm * 0.0) / (br + brm * bi);
    } else if (bim == brm) {
      if (br > 0.0) {
        ub_br = 0.5;
      } else {
        ub_br = -0.5;
      }

      if (bi > 0.0) {
        ub_bi = 0.5;
      } else {
        ub_bi = -0.5;
      }

      g_y_re = (ar * ub_br + 0.0 * ub_bi) / brm;
    } else {
      brm = br / bi;
      g_y_re = brm * ar / (bi + brm * br);
    }
  }

  ar = 3.0 * x09 * x19;
  if (y_im == 0.0) {
    e_y_re = ar / y_re;
  } else if (y_re == 0.0) {
    if (ar == 0.0) {
      e_y_re = 0.0 / y_im;
    } else {
      e_y_re = 0.0;
    }
  } else {
    brm = fabs(y_re);
    bim = fabs(y_im);
    if (brm > bim) {
      brm = y_im / y_re;
      e_y_re = (ar + brm * 0.0) / (y_re + brm * y_im);
    } else if (bim == brm) {
      if (y_re > 0.0) {
        ac_y_re = 0.5;
      } else {
        ac_y_re = -0.5;
      }

      if (y_im > 0.0) {
        yb_y_im = 0.5;
      } else {
        yb_y_im = -0.5;
      }

      e_y_re = (ar * ac_y_re + 0.0 * yb_y_im) / brm;
    } else {
      brm = y_re / y_im;
      e_y_re = brm * ar / (y_im + brm * y_re);
    }
  }

  ar = 18.0 * x11 * x14;
  br = 125.0 * lb_y_re;
  bi = 125.0 * lb_y_im;
  if (bi == 0.0) {
    d_y_im = ar / br;
  } else if (br == 0.0) {
    if (ar == 0.0) {
      d_y_im = 0.0 / bi;
    } else {
      d_y_im = 0.0;
    }
  } else {
    brm = fabs(br);
    bim = fabs(bi);
    if (brm > bim) {
      brm = bi / br;
      d_y_im = (ar + brm * 0.0) / (br + brm * bi);
    } else if (bim == brm) {
      if (br > 0.0) {
        vb_br = 0.5;
      } else {
        vb_br = -0.5;
      }

      if (bi > 0.0) {
        vb_bi = 0.5;
      } else {
        vb_bi = -0.5;
      }

      d_y_im = (ar * vb_br + 0.0 * vb_bi) / brm;
    } else {
      brm = br / bi;
      d_y_im = brm * ar / (bi + brm * br);
    }
  }

  ar = 3.0 * x11 * x17;
  br = 125.0 * mb_y_re;
  bi = 125.0 * mb_y_im;
  if (bi == 0.0) {
    d_y_re = ar / br;
  } else if (br == 0.0) {
    if (ar == 0.0) {
      d_y_re = 0.0 / bi;
    } else {
      d_y_re = 0.0;
    }
  } else {
    brm = fabs(br);
    bim = fabs(bi);
    if (brm > bim) {
      brm = bi / br;
      d_y_re = (ar + brm * 0.0) / (br + brm * bi);
    } else if (bim == brm) {
      if (br > 0.0) {
        wb_br = 0.5;
      } else {
        wb_br = -0.5;
      }

      if (bi > 0.0) {
        wb_bi = 0.5;
      } else {
        wb_bi = -0.5;
      }

      d_y_re = (ar * wb_br + 0.0 * wb_bi) / brm;
    } else {
      brm = br / bi;
      d_y_re = brm * ar / (bi + brm * br);
    }
  }

  ar = 18.0 * x12 * x15;
  br = 125.0 * nb_y_re;
  bi = 125.0 * nb_y_im;
  if (bi == 0.0) {
    c_y_re = ar / br;
  } else if (br == 0.0) {
    if (ar == 0.0) {
      c_y_re = 0.0 / bi;
    } else {
      c_y_re = 0.0;
    }
  } else {
    brm = fabs(br);
    bim = fabs(bi);
    if (brm > bim) {
      brm = bi / br;
      c_y_re = (ar + brm * 0.0) / (br + brm * bi);
    } else if (bim == brm) {
      if (br > 0.0) {
        xb_br = 0.5;
      } else {
        xb_br = -0.5;
      }

      if (bi > 0.0) {
        xb_bi = 0.5;
      } else {
        xb_bi = -0.5;
      }

      c_y_re = (ar * xb_br + 0.0 * xb_bi) / brm;
    } else {
      brm = br / bi;
      c_y_re = brm * ar / (bi + brm * br);
    }
  }

  ar = 9.0 * x14 * x17;
  br = 625.0 * r;
  bi = 625.0 * t_s_im_tmp;
  if (bi == 0.0) {
    c_y_im = ar / br;
  } else if (br == 0.0) {
    if (ar == 0.0) {
      c_y_im = 0.0 / bi;
    } else {
      c_y_im = 0.0;
    }
  } else {
    brm = fabs(br);
    bim = fabs(bi);
    if (brm > bim) {
      brm = bi / br;
      c_y_im = (ar + brm * 0.0) / (br + brm * bi);
    } else if (bim == brm) {
      if (br > 0.0) {
        yb_br = 0.5;
      } else {
        yb_br = -0.5;
      }

      if (bi > 0.0) {
        yb_bi = 0.5;
      } else {
        yb_bi = -0.5;
      }

      c_y_im = (ar * yb_br + 0.0 * yb_bi) / brm;
    } else {
      brm = br / bi;
      c_y_im = brm * ar / (bi + brm * br);
    }
  }

  ar = 3.0 * x12 * x18;
  br = 125.0 * ob_y_re;
  bi = 125.0 * ob_y_im;
  if (bi == 0.0) {
    b_y_im = ar / br;
  } else if (br == 0.0) {
    if (ar == 0.0) {
      b_y_im = 0.0 / bi;
    } else {
      b_y_im = 0.0;
    }
  } else {
    brm = fabs(br);
    bim = fabs(bi);
    if (brm > bim) {
      brm = bi / br;
      b_y_im = (ar + brm * 0.0) / (br + brm * bi);
    } else if (bim == brm) {
      if (br > 0.0) {
        ac_br = 0.5;
      } else {
        ac_br = -0.5;
      }

      if (bi > 0.0) {
        ac_bi = 0.5;
      } else {
        ac_bi = -0.5;
      }

      b_y_im = (ar * ac_br + 0.0 * ac_bi) / brm;
    } else {
      brm = br / bi;
      b_y_im = brm * ar / (bi + brm * br);
    }
  }

  ar = 18.0 * x13 * x16;
  br = 125.0 * pb_y_re;
  bi = 125.0 * pb_y_im;
  if (bi == 0.0) {
    y_im = ar / br;
  } else if (br == 0.0) {
    if (ar == 0.0) {
      y_im = 0.0 / bi;
    } else {
      y_im = 0.0;
    }
  } else {
    brm = fabs(br);
    bim = fabs(bi);
    if (brm > bim) {
      brm = bi / br;
      y_im = (ar + brm * 0.0) / (br + brm * bi);
    } else if (bim == brm) {
      if (br > 0.0) {
        bc_br = 0.5;
      } else {
        bc_br = -0.5;
      }

      if (bi > 0.0) {
        bc_bi = 0.5;
      } else {
        bc_bi = -0.5;
      }

      y_im = (ar * bc_br + 0.0 * bc_bi) / brm;
    } else {
      brm = br / bi;
      y_im = brm * ar / (bi + brm * br);
    }
  }

  ar = 9.0 * x15 * x18;
  br = 625.0 * r;
  bi = 625.0 * t_s_im_tmp;
  if (bi == 0.0) {
    b_y_re = ar / br;
  } else if (br == 0.0) {
    if (ar == 0.0) {
      b_y_re = 0.0 / bi;
    } else {
      b_y_re = 0.0;
    }
  } else {
    brm = fabs(br);
    bim = fabs(bi);
    if (brm > bim) {
      brm = bi / br;
      b_y_re = (ar + brm * 0.0) / (br + brm * bi);
    } else if (bim == brm) {
      if (br > 0.0) {
        cc_br = 0.5;
      } else {
        cc_br = -0.5;
      }

      if (bi > 0.0) {
        cc_bi = 0.5;
      } else {
        cc_bi = -0.5;
      }

      b_y_re = (ar * cc_br + 0.0 * cc_bi) / brm;
    } else {
      brm = br / bi;
      b_y_re = brm * ar / (bi + brm * br);
    }
  }

  ar = 3.0 * x13 * x19;
  br = 125.0 * x_re;
  bi = 125.0 * x_im;
  if (bi == 0.0) {
    y_re = ar / br;
  } else if (br == 0.0) {
    if (ar == 0.0) {
      y_re = 0.0 / bi;
    } else {
      y_re = 0.0;
    }
  } else {
    brm = fabs(br);
    bim = fabs(bi);
    if (brm > bim) {
      brm = bi / br;
      y_re = (ar + brm * 0.0) / (br + brm * bi);
    } else if (bim == brm) {
      if (br > 0.0) {
        dc_br = 0.5;
      } else {
        dc_br = -0.5;
      }

      if (bi > 0.0) {
        dc_bi = 0.5;
      } else {
        dc_bi = -0.5;
      }

      y_re = (ar * dc_br + 0.0 * dc_bi) / brm;
    } else {
      brm = br / bi;
      y_re = brm * ar / (bi + brm * br);
    }
  }

  ar = 9.0 * x16 * x19;
  br = 625.0 * r;
  bi = 625.0 * t_s_im_tmp;
  if (bi == 0.0) {
    r = ar / br;
  } else if (br == 0.0) {
    if (ar == 0.0) {
      r = 0.0 / bi;
    } else {
      r = 0.0;
    }
  } else {
    brm = fabs(br);
    bim = fabs(bi);
    if (brm > bim) {
      brm = bi / br;
      r = (ar + brm * 0.0) / (br + brm * bi);
    } else if (bim == brm) {
      if (br > 0.0) {
        ec_br = 0.5;
      } else {
        ec_br = -0.5;
      }

      if (bi > 0.0) {
        ec_bi = 0.5;
      } else {
        ec_bi = -0.5;
      }

      r = (ar * ec_br + 0.0 * ec_bi) / brm;
    } else {
      brm = br / bi;
      r = brm * ar / (bi + brm * br);
    }
  }

  return ((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((t_s.re +
    re) + b_re) + c_re) + d_re) + e_re) + f_re) + g_re) + h_re) + i_re) + j_re)
    + k_re) + l_re) + m_re) + n_re) + o_re) + p_re) + q_re) + r_re) + s_re) +
    t_re) + u_re) + v_re) + w_re) + ab_re) + bb_re) + r_y_im) + cb_re) - r_y_re)
    + q_y_im) - q_y_re) - p_y_im) - p_y_re) + o_y_im) - o_y_re) + n_y_im) -
    n_y_re) - l_y_im) - m_y_re) + m_y_im) - l_y_re) + k_y_im) - k_y_re) + j_y_im)
    - h_y_im) - i_y_re) - i_y_im) + j_y_re) - h_y_re) + g_y_im) - e_y_im) -
                     f_y_re) - f_y_im) + g_y_re) - e_y_re) - d_y_im) + d_y_re) -
               c_y_re) - c_y_im) + b_y_im) - y_im) - b_y_re) + y_re) - r;
}

/* End of code generation (quadf_cost.c) */
