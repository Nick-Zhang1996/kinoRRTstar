/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * segment_cost.c
 *
 * Code generation for function 'segment_cost'
 *
 */

/* Include files */
#include <math.h>
#include "rt_nonfinite.h"
#include "benchmarkRRT.h"
#include "segment_cost.h"
#include "cost_eval.h"
#include "roots.h"
#include "opt_time.h"
#include "benchmarkRRT_rtwutil.h"

/* Function Definitions */
void b_segment_cost(const double from_node[4], const double to_point_data[],
                    creal_T *cost, creal_T *Tf)
{
  double dv3[5];
  creal_T Tf_data[4];
  int Tf_size[1];
  int nx;
  int k;
  double x_data[4];
  double y_data[4];
  boolean_T tmp_data[4];
  int trueCount;
  int partialTrueCount;
  signed char b_tmp_data[4];
  creal_T varargin_1_data[4];
  creal_T dc2;
  boolean_T SCALEA;
  double ma;
  double mb;
  boolean_T SCALEB;
  double x;
  double br;
  double Mb;
  double Ma;
  opt_time(from_node[0], from_node[1], from_node[2], from_node[3],
           to_point_data[0], to_point_data[1], to_point_data[2], to_point_data[3],
           dv3);
  roots(dv3, Tf_data, Tf_size);
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
  for (k = 0; k <= nx; k++) {
    if (tmp_data[k]) {
      trueCount++;
    }
  }

  partialTrueCount = 0;
  for (k = 0; k <= nx; k++) {
    if (tmp_data[k]) {
      Tf_data[partialTrueCount] = Tf_data[k];
      partialTrueCount++;
    }
  }

  /* NOTE */
  nx = trueCount - 1;
  trueCount = 0;
  for (k = 0; k <= nx; k++) {
    if (Tf_data[k].re >= 0.0) {
      trueCount++;
    }
  }

  partialTrueCount = 0;
  for (k = 0; k <= nx; k++) {
    if (Tf_data[k].re >= 0.0) {
      b_tmp_data[partialTrueCount] = (signed char)(k + 1);
      partialTrueCount++;
    }
  }

  for (k = 0; k < trueCount; k++) {
    varargin_1_data[k] = Tf_data[b_tmp_data[k] - 1];
  }

  *Tf = varargin_1_data[0];
  for (k = 2; k <= trueCount; k++) {
    dc2 = varargin_1_data[k - 1];
    if (rtIsNaN(dc2.re) || rtIsNaN(varargin_1_data[k - 1].im)) {
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
        br = fabs(varargin_1_data[k - 1].im);
        if (ma > Mb) {
          Ma = ma;
          ma = Mb;
        } else {
          Ma = Mb;
        }

        if (mb > br) {
          Mb = mb;
          mb = br;
        } else {
          Mb = br;
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
            Mb = varargin_1_data[k - 1].im;
            if (x > 0.78539816339744828) {
              if (x > 2.3561944901923448) {
                x = -Tf->im;
                br = -Mb;
              } else {
                x = -Tf->re;
                br = -br;
              }
            } else if (x > -0.78539816339744828) {
              x = Tf->im;
              br = Mb;
            } else if (x > -2.3561944901923448) {
              x = Tf->re;
            } else {
              x = -Tf->im;
              br = -Mb;
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
      *Tf = dc2;
    }
  }

  *cost = cost_eval(*Tf, from_node[0], from_node[1], from_node[2], from_node[3],
                    to_point_data[0], to_point_data[1], to_point_data[2],
                    to_point_data[3]);
}

void c_segment_cost(const double from_node_data[], creal_T *cost, creal_T *Tf)
{
  double dv5[5];
  creal_T Tf_data[4];
  int Tf_size[1];
  int nx;
  int k;
  double x_data[4];
  double y_data[4];
  boolean_T tmp_data[4];
  int trueCount;
  int partialTrueCount;
  signed char b_tmp_data[4];
  creal_T varargin_1_data[4];
  creal_T dc4;
  boolean_T SCALEA;
  double ma;
  double mb;
  boolean_T SCALEB;
  double x;
  double br;
  double Mb;
  double Ma;
  opt_time(from_node_data[0], from_node_data[1], from_node_data[2],
           from_node_data[3], 18.0, 18.0, 0.0, 0.0, dv5);
  roots(dv5, Tf_data, Tf_size);
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
  for (k = 0; k <= nx; k++) {
    if (tmp_data[k]) {
      trueCount++;
    }
  }

  partialTrueCount = 0;
  for (k = 0; k <= nx; k++) {
    if (tmp_data[k]) {
      Tf_data[partialTrueCount] = Tf_data[k];
      partialTrueCount++;
    }
  }

  /* NOTE */
  nx = trueCount - 1;
  trueCount = 0;
  for (k = 0; k <= nx; k++) {
    if (Tf_data[k].re >= 0.0) {
      trueCount++;
    }
  }

  partialTrueCount = 0;
  for (k = 0; k <= nx; k++) {
    if (Tf_data[k].re >= 0.0) {
      b_tmp_data[partialTrueCount] = (signed char)(k + 1);
      partialTrueCount++;
    }
  }

  for (k = 0; k < trueCount; k++) {
    varargin_1_data[k] = Tf_data[b_tmp_data[k] - 1];
  }

  *Tf = varargin_1_data[0];
  for (k = 2; k <= trueCount; k++) {
    dc4 = varargin_1_data[k - 1];
    if (rtIsNaN(dc4.re) || rtIsNaN(varargin_1_data[k - 1].im)) {
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
        br = fabs(varargin_1_data[k - 1].im);
        if (ma > Mb) {
          Ma = ma;
          ma = Mb;
        } else {
          Ma = Mb;
        }

        if (mb > br) {
          Mb = mb;
          mb = br;
        } else {
          Mb = br;
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
            Mb = varargin_1_data[k - 1].im;
            if (x > 0.78539816339744828) {
              if (x > 2.3561944901923448) {
                x = -Tf->im;
                br = -Mb;
              } else {
                x = -Tf->re;
                br = -br;
              }
            } else if (x > -0.78539816339744828) {
              x = Tf->im;
              br = Mb;
            } else if (x > -2.3561944901923448) {
              x = Tf->re;
            } else {
              x = -Tf->im;
              br = -Mb;
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
      *Tf = dc4;
    }
  }

  *cost = cost_eval(*Tf, from_node_data[0], from_node_data[1], from_node_data[2],
                    from_node_data[3], 18.0, 18.0, 0.0, 0.0);
}

void segment_cost(const double from_node_data[], const double to_point[4],
                  creal_T *cost, creal_T *Tf)
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
  int partialTrueCount;
  signed char b_tmp_data[4];
  creal_T varargin_1_data[4];
  creal_T dc1;
  boolean_T SCALEA;
  double ma;
  double mb;
  boolean_T SCALEB;
  double x;
  double br;
  double Mb;
  double Ma;
  opt_time(from_node_data[0], from_node_data[1], from_node_data[2],
           from_node_data[3], to_point[0], to_point[1], to_point[2], to_point[3],
           dv2);
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
  for (k = 0; k <= nx; k++) {
    if (tmp_data[k]) {
      trueCount++;
    }
  }

  partialTrueCount = 0;
  for (k = 0; k <= nx; k++) {
    if (tmp_data[k]) {
      Tf_data[partialTrueCount] = Tf_data[k];
      partialTrueCount++;
    }
  }

  /* NOTE */
  nx = trueCount - 1;
  trueCount = 0;
  for (k = 0; k <= nx; k++) {
    if (Tf_data[k].re >= 0.0) {
      trueCount++;
    }
  }

  partialTrueCount = 0;
  for (k = 0; k <= nx; k++) {
    if (Tf_data[k].re >= 0.0) {
      b_tmp_data[partialTrueCount] = (signed char)(k + 1);
      partialTrueCount++;
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
        br = fabs(varargin_1_data[k - 1].im);
        if (ma > Mb) {
          Ma = ma;
          ma = Mb;
        } else {
          Ma = Mb;
        }

        if (mb > br) {
          Mb = mb;
          mb = br;
        } else {
          Mb = br;
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
            Mb = varargin_1_data[k - 1].im;
            if (x > 0.78539816339744828) {
              if (x > 2.3561944901923448) {
                x = -Tf->im;
                br = -Mb;
              } else {
                x = -Tf->re;
                br = -br;
              }
            } else if (x > -0.78539816339744828) {
              x = Tf->im;
              br = Mb;
            } else if (x > -2.3561944901923448) {
              x = Tf->re;
            } else {
              x = -Tf->im;
              br = -Mb;
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

  *cost = cost_eval(*Tf, from_node_data[0], from_node_data[1], from_node_data[2],
                    from_node_data[3], to_point[0], to_point[1], to_point[2],
                    to_point[3]);
}

/* End of code generation (segment_cost.c) */
