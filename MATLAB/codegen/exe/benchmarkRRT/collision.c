/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * collision.c
 *
 * Code generation for function 'collision'
 *
 */

/* Include files */
#include <math.h>
#include "rt_nonfinite.h"
#include "benchmarkRRT.h"
#include "collision.h"
#include "state_traj.h"
#include "linspace.h"
#include "cost_eval.h"
#include "roots.h"
#include "opt_time.h"
#include "benchmarkRRT_rtwutil.h"
#include "benchmarkRRT_data.h"

/* Function Definitions */
void b_collision(const double parent_data[], double *collision_flag, creal_T *Tf)
{
  double dv4[5];
  creal_T tf_temp_data[4];
  int tf_temp_size[1];
  int nx;
  int k;
  double x_data[4];
  double y_data[4];
  boolean_T tmp_data[4];
  int trueCount;
  int i;
  signed char b_tmp_data[4];
  creal_T varargin_1_data[4];
  creal_T dc3;
  boolean_T SCALEA;
  double ma;
  creal_T t[12];
  double mb;
  boolean_T SCALEB;
  int exitg3;
  double b_state[40];
  double x;
  double br;
  int exitg2;
  double Mb;
  double bi;
  double Ma;
  int exitg1;
  static const signed char iv3[5] = { 6, 15, 10, 5, 15 };

  static const signed char iv4[5] = { 6, 14, 11, 15, 5 };

  *collision_flag = 0.0;
  opt_time(parent_data[0], parent_data[1], parent_data[2], parent_data[3], 18.0,
           18.0, 0.0, 0.0, dv4);
  roots(dv4, tf_temp_data, tf_temp_size);

  /*  DEBUG */
  /*          poly = opt_time(parent(1), parent(2), parent(3), parent(4), node(1), node(2), node(3), node(4)); */
  /*          fprintf("\n poly: "); */
  /*          for i = 1:size(poly,2) */
  /*              fprintf("%.5f, ", poly(1,i)); */
  /*          end */
  /*          fprintf("\n"); */
  nx = tf_temp_size[0];
  for (k = 0; k < nx; k++) {
    x_data[k] = tf_temp_data[k].im;
  }

  nx = tf_temp_size[0];
  for (k = 0; k < nx; k++) {
    y_data[k] = fabs(x_data[k]);
  }

  nx = (signed char)tf_temp_size[0];
  for (k = 0; k < nx; k++) {
    tmp_data[k] = (y_data[k] < 0.0001);
  }

  nx = (signed char)tf_temp_size[0] - 1;
  trueCount = 0;
  for (i = 0; i <= nx; i++) {
    if (tmp_data[i]) {
      trueCount++;
    }
  }

  k = 0;
  for (i = 0; i <= nx; i++) {
    if (tmp_data[i]) {
      tf_temp_data[k] = tf_temp_data[i];
      k++;
    }
  }

  /* NOTE */
  nx = trueCount - 1;
  trueCount = 0;
  for (i = 0; i <= nx; i++) {
    if (tf_temp_data[i].re >= 0.0) {
      trueCount++;
    }
  }

  k = 0;
  for (i = 0; i <= nx; i++) {
    if (tf_temp_data[i].re >= 0.0) {
      b_tmp_data[k] = (signed char)(i + 1);
      k++;
    }
  }

  for (k = 0; k < trueCount; k++) {
    varargin_1_data[k] = tf_temp_data[b_tmp_data[k] - 1];
  }

  *Tf = varargin_1_data[0];
  for (k = 2; k <= trueCount; k++) {
    dc3 = varargin_1_data[k - 1];
    if (rtIsNaN(dc3.re) || rtIsNaN(varargin_1_data[k - 1].im)) {
      SCALEA = false;
    } else if (rtIsNaN(Tf->re) || rtIsNaN(Tf->im)) {
      SCALEA = true;
    } else {
      ma = fabs(Tf->re);
      if ((ma > 8.9884656743115785E+307) || (fabs(Tf->im) >
           8.9884656743115785E+307)) {
        SCALEA = true;
      } else {
        SCALEA = false;
      }

      mb = fabs(varargin_1_data[k - 1].re);
      if ((mb > 8.9884656743115785E+307) || (fabs(varargin_1_data[k - 1].im) >
           8.9884656743115785E+307)) {
        SCALEB = true;
      } else {
        SCALEB = false;
      }

      if (SCALEA || SCALEB) {
        x = rt_hypotd_snf(Tf->re / 2.0, Tf->im / 2.0);
        br = rt_hypotd_snf(varargin_1_data[k - 1].re / 2.0, varargin_1_data[k -
                           1].im / 2.0);
      } else {
        x = rt_hypotd_snf(Tf->re, Tf->im);
        br = rt_hypotd_snf(varargin_1_data[k - 1].re, varargin_1_data[k - 1].im);
      }

      if (x == br) {
        Mb = fabs(Tf->im);
        bi = fabs(varargin_1_data[k - 1].im);
        if (ma > Mb) {
          Ma = ma;
          ma = Mb;
        } else {
          Ma = Mb;
        }

        if (mb > bi) {
          Mb = mb;
          mb = bi;
        } else {
          Mb = bi;
        }

        if (Ma > Mb) {
          if (ma < mb) {
            x = Ma - Mb;
            br = (ma / 2.0 + mb / 2.0) / (Ma / 2.0 + Mb / 2.0) * (mb - ma);
          } else {
            x = Ma;
            br = Mb;
          }
        } else if (Ma < Mb) {
          if (ma > mb) {
            br = Mb - Ma;
            x = (ma / 2.0 + mb / 2.0) / (Ma / 2.0 + Mb / 2.0) * (ma - mb);
          } else {
            x = Ma;
            br = Mb;
          }
        } else {
          x = ma;
          br = mb;
        }

        if (x == br) {
          x = rt_atan2d_snf(Tf->im, Tf->re);
          br = rt_atan2d_snf(varargin_1_data[k - 1].im, varargin_1_data[k - 1].
                             re);
          if (x == br) {
            br = varargin_1_data[k - 1].re;
            bi = varargin_1_data[k - 1].im;
            if (x > 0.78539816339744828) {
              if (x > 2.3561944901923448) {
                x = -Tf->im;
                br = -bi;
              } else {
                x = -Tf->re;
                br = -br;
              }
            } else if (x > -0.78539816339744828) {
              x = Tf->im;
              br = bi;
            } else if (x > -2.3561944901923448) {
              x = Tf->re;
            } else {
              x = -Tf->im;
              br = -bi;
            }

            if (x == br) {
              x = 0.0;
              br = 0.0;
            }
          }
        }
      }

      SCALEA = (x > br);
    }

    if (SCALEA) {
      *Tf = dc3;
    }
  }

  /* DEBUG */
  /*          for i = 1:size(tf_temp,1) */
  /*              fprintf("\n tf_temp %.0f, %.2f \n ",i,tf_temp(i,1)); */
  /*          end */
  /*          fprintf("collision Tf = %.2f \n",Tf); */
  /* Tf = sum(temp); */
  linspace(*Tf, t);
  for (i = 0; i < 10; i++) {
    /*  added real() */
    state_traj(t[i + 1], *Tf, parent_data[0], parent_data[1], parent_data[2],
               parent_data[3], 18.0, 18.0, 0.0, 0.0, tf_temp_data);
    nx = i << 2;
    b_state[nx] = tf_temp_data[0].re;
    b_state[1 + nx] = tf_temp_data[1].re;
    b_state[2 + nx] = tf_temp_data[2].re;
    b_state[3 + nx] = tf_temp_data[3].re;
  }

  nx = 0;
  do {
    exitg3 = 0;
    if (nx < 10) {
      i = 0;
      do {
        exitg2 = 0;
        if (i < 2) {
          Mb = b_state[i + (nx << 2)];
          if ((Mb > 20.0) || (Mb < 0.0)) {
            *collision_flag = 1.0;
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
            k = nx << 2;
            bi = b_state[k] - (double)iv3[i];
            bi *= bi;
            Mb = bi;
            bi = b_state[1 + k] - (double)iv4[i];
            bi *= bi;
            if (Mb + bi < (dv0[i] + 0.1) * (dv0[i] + 0.1)) {
              /*  (norm([p(1);p(2)]-[world.cx(i); world.cy(i)])<=1*world.radius(i)) */
              *collision_flag = 1.0;
              exitg1 = 1;
            } else {
              i++;
            }
          } else {
            nx++;
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
}

void collision(const double parent_data[], const double node[4], double
               *collision_flag, creal_T *Tf)
{
  double dv1[5];
  creal_T tf_temp_data[4];
  int tf_temp_size[1];
  int nx;
  int k;
  double x_data[4];
  double y_data[4];
  boolean_T tmp_data[4];
  int trueCount;
  int i;
  signed char b_tmp_data[4];
  creal_T varargin_1_data[4];
  creal_T dc0;
  boolean_T SCALEA;
  double ma;
  creal_T t[12];
  double mb;
  boolean_T SCALEB;
  int exitg3;
  double b_state[40];
  double x;
  double br;
  int exitg2;
  double Mb;
  double bi;
  double Ma;
  int exitg1;
  static const signed char iv1[5] = { 6, 15, 10, 5, 15 };

  static const signed char iv2[5] = { 6, 14, 11, 15, 5 };

  *collision_flag = 0.0;
  opt_time(parent_data[0], parent_data[1], parent_data[2], parent_data[3], node
           [0], node[1], node[2], node[3], dv1);
  roots(dv1, tf_temp_data, tf_temp_size);

  /*  DEBUG */
  /*          poly = opt_time(parent(1), parent(2), parent(3), parent(4), node(1), node(2), node(3), node(4)); */
  /*          fprintf("\n poly: "); */
  /*          for i = 1:size(poly,2) */
  /*              fprintf("%.5f, ", poly(1,i)); */
  /*          end */
  /*          fprintf("\n"); */
  nx = tf_temp_size[0];
  for (k = 0; k < nx; k++) {
    x_data[k] = tf_temp_data[k].im;
  }

  nx = tf_temp_size[0];
  for (k = 0; k < nx; k++) {
    y_data[k] = fabs(x_data[k]);
  }

  nx = (signed char)tf_temp_size[0];
  for (k = 0; k < nx; k++) {
    tmp_data[k] = (y_data[k] < 0.0001);
  }

  nx = (signed char)tf_temp_size[0] - 1;
  trueCount = 0;
  for (i = 0; i <= nx; i++) {
    if (tmp_data[i]) {
      trueCount++;
    }
  }

  k = 0;
  for (i = 0; i <= nx; i++) {
    if (tmp_data[i]) {
      tf_temp_data[k] = tf_temp_data[i];
      k++;
    }
  }

  /* NOTE */
  nx = trueCount - 1;
  trueCount = 0;
  for (i = 0; i <= nx; i++) {
    if (tf_temp_data[i].re >= 0.0) {
      trueCount++;
    }
  }

  k = 0;
  for (i = 0; i <= nx; i++) {
    if (tf_temp_data[i].re >= 0.0) {
      b_tmp_data[k] = (signed char)(i + 1);
      k++;
    }
  }

  for (k = 0; k < trueCount; k++) {
    varargin_1_data[k] = tf_temp_data[b_tmp_data[k] - 1];
  }

  *Tf = varargin_1_data[0];
  for (k = 2; k <= trueCount; k++) {
    dc0 = varargin_1_data[k - 1];
    if (rtIsNaN(dc0.re) || rtIsNaN(varargin_1_data[k - 1].im)) {
      SCALEA = false;
    } else if (rtIsNaN(Tf->re) || rtIsNaN(Tf->im)) {
      SCALEA = true;
    } else {
      ma = fabs(Tf->re);
      if ((ma > 8.9884656743115785E+307) || (fabs(Tf->im) >
           8.9884656743115785E+307)) {
        SCALEA = true;
      } else {
        SCALEA = false;
      }

      mb = fabs(varargin_1_data[k - 1].re);
      if ((mb > 8.9884656743115785E+307) || (fabs(varargin_1_data[k - 1].im) >
           8.9884656743115785E+307)) {
        SCALEB = true;
      } else {
        SCALEB = false;
      }

      if (SCALEA || SCALEB) {
        x = rt_hypotd_snf(Tf->re / 2.0, Tf->im / 2.0);
        br = rt_hypotd_snf(varargin_1_data[k - 1].re / 2.0, varargin_1_data[k -
                           1].im / 2.0);
      } else {
        x = rt_hypotd_snf(Tf->re, Tf->im);
        br = rt_hypotd_snf(varargin_1_data[k - 1].re, varargin_1_data[k - 1].im);
      }

      if (x == br) {
        Mb = fabs(Tf->im);
        bi = fabs(varargin_1_data[k - 1].im);
        if (ma > Mb) {
          Ma = ma;
          ma = Mb;
        } else {
          Ma = Mb;
        }

        if (mb > bi) {
          Mb = mb;
          mb = bi;
        } else {
          Mb = bi;
        }

        if (Ma > Mb) {
          if (ma < mb) {
            x = Ma - Mb;
            br = (ma / 2.0 + mb / 2.0) / (Ma / 2.0 + Mb / 2.0) * (mb - ma);
          } else {
            x = Ma;
            br = Mb;
          }
        } else if (Ma < Mb) {
          if (ma > mb) {
            br = Mb - Ma;
            x = (ma / 2.0 + mb / 2.0) / (Ma / 2.0 + Mb / 2.0) * (ma - mb);
          } else {
            x = Ma;
            br = Mb;
          }
        } else {
          x = ma;
          br = mb;
        }

        if (x == br) {
          x = rt_atan2d_snf(Tf->im, Tf->re);
          br = rt_atan2d_snf(varargin_1_data[k - 1].im, varargin_1_data[k - 1].
                             re);
          if (x == br) {
            br = varargin_1_data[k - 1].re;
            bi = varargin_1_data[k - 1].im;
            if (x > 0.78539816339744828) {
              if (x > 2.3561944901923448) {
                x = -Tf->im;
                br = -bi;
              } else {
                x = -Tf->re;
                br = -br;
              }
            } else if (x > -0.78539816339744828) {
              x = Tf->im;
              br = bi;
            } else if (x > -2.3561944901923448) {
              x = Tf->re;
            } else {
              x = -Tf->im;
              br = -bi;
            }

            if (x == br) {
              x = 0.0;
              br = 0.0;
            }
          }
        }
      }

      SCALEA = (x > br);
    }

    if (SCALEA) {
      *Tf = dc0;
    }
  }

  /* DEBUG */
  /*          for i = 1:size(tf_temp,1) */
  /*              fprintf("\n tf_temp %.0f, %.2f \n ",i,tf_temp(i,1)); */
  /*          end */
  /*          fprintf("collision Tf = %.2f \n",Tf); */
  /* Tf = sum(temp); */
  linspace(*Tf, t);
  for (i = 0; i < 10; i++) {
    /*  added real() */
    state_traj(t[i + 1], *Tf, parent_data[0], parent_data[1], parent_data[2],
               parent_data[3], node[0], node[1], node[2], node[3], tf_temp_data);
    nx = i << 2;
    b_state[nx] = tf_temp_data[0].re;
    b_state[1 + nx] = tf_temp_data[1].re;
    b_state[2 + nx] = tf_temp_data[2].re;
    b_state[3 + nx] = tf_temp_data[3].re;
  }

  nx = 0;
  do {
    exitg3 = 0;
    if (nx < 10) {
      i = 0;
      do {
        exitg2 = 0;
        if (i < 2) {
          Mb = b_state[i + (nx << 2)];
          if ((Mb > 20.0) || (Mb < 0.0)) {
            *collision_flag = 1.0;
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
            k = nx << 2;
            bi = b_state[k] - (double)iv1[i];
            bi *= bi;
            Mb = bi;
            bi = b_state[1 + k] - (double)iv2[i];
            bi *= bi;
            if (Mb + bi < (dv0[i] + 0.1) * (dv0[i] + 0.1)) {
              /*  (norm([p(1);p(2)]-[world.cx(i); world.cy(i)])<=1*world.radius(i)) */
              *collision_flag = 1.0;
              exitg1 = 1;
            } else {
              i++;
            }
          } else {
            nx++;
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
}

/* End of code generation (collision.c) */
