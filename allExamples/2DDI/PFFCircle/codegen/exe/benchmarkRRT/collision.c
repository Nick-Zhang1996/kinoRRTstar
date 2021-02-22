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
#include "DI_state.h"
#include "linspace.h"
#include "DI_cost.h"
#include "roots.h"
#include "DI_time.h"
#include "benchmarkRRT_rtwutil.h"
#include "benchmarkRRT_data.h"

/* Function Definitions */
void collision(const double parent_data[], double *collision_flag, creal_T *Tf)
{
  double dv2[5];
  creal_T Tf_data[4];
  int Tf_size[1];
  int nx;
  int k;
  double x_data[4];
  double y_data[4];
  boolean_T tmp_data[4];
  int trueCount;
  int i;
  signed char b_tmp_data[4];
  creal_T varargin_1_data[4];
  creal_T t[12];
  creal_T dc1;
  boolean_T SCALEA;
  double b_state[40];
  double ma;
  int exitg3;
  double mb;
  int exitg2;
  boolean_T SCALEB;
  double Mb;
  double x;
  double br;
  double bi;
  int exitg1;
  double Ma;
  static const signed char iv2[5] = { 6, 15, 10, 5, 15 };

  static const signed char iv3[5] = { 6, 14, 11, 15, 5 };

  *collision_flag = 0.0;

  /*      Tf = norm(parent(1 : 4) - node(1 : 4))/0.8; */
  DI_time(parent_data[0], parent_data[1], parent_data[2], parent_data[3], 18.0,
          18.0, 0.0, 0.0, dv2);
  roots(dv2, Tf_data, Tf_size);
  nx = Tf_size[0];
  for (k = 0; k < nx; k++) {
    x_data[k] = Tf_data[k].im;
  }

  nx = Tf_size[0];
  for (k = 0; k < nx; k++) {
    y_data[k] = fabs(x_data[k]);
  }

  nx = (signed char)Tf_size[0];
  for (k = 0; k < nx; k++) {
    tmp_data[k] = (y_data[k] < 0.0001);
  }

  nx = (signed char)Tf_size[0] - 1;
  trueCount = 0;
  for (i = 0; i <= nx; i++) {
    if (tmp_data[i]) {
      trueCount++;
    }
  }

  k = 0;
  for (i = 0; i <= nx; i++) {
    if (tmp_data[i]) {
      Tf_data[k] = Tf_data[i];
      k++;
    }
  }

  nx = trueCount - 1;
  trueCount = 0;
  for (i = 0; i <= nx; i++) {
    if (Tf_data[i].re >= 0.0) {
      trueCount++;
    }
  }

  k = 0;
  for (i = 0; i <= nx; i++) {
    if (Tf_data[i].re >= 0.0) {
      b_tmp_data[k] = (signed char)(i + 1);
      k++;
    }
  }

  for (k = 0; k < trueCount; k++) {
    varargin_1_data[k] = Tf_data[b_tmp_data[k] - 1];
  }

  *Tf = varargin_1_data[0];
  for (k = 2; k <= trueCount; k++) {
    dc1 = varargin_1_data[k - 1];
    if (rtIsNaN(dc1.re) || rtIsNaN(varargin_1_data[k - 1].im)) {
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
      *Tf = dc1;
    }
  }

  linspace(*Tf, t);
  for (i = 0; i < 10; i++) {
    DI_state(t[i + 1], *Tf, parent_data[0], parent_data[1], parent_data[2],
             parent_data[3], *(double (*)[4])&b_state[i << 2]);
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
            bi = b_state[k] - (double)iv2[i];
            bi *= bi;
            Mb = bi;
            bi = b_state[1 + k] - (double)iv3[i];
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
