/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * collisionFreeVel.c
 *
 * Code generation for function 'collisionFreeVel'
 *
 */

/* Include files */
#include <math.h>
#include "rt_nonfinite.h"
#include "benchmarkRRT.h"
#include "collisionFreeVel.h"
#include "quadf_statePFFull.h"
#include "quadf_statePFF.h"
#include "linspace.h"
#include "quadf_costPFF.h"
#include "roots.h"
#include "quadf_timePFF.h"

/* Function Definitions */
void collisionFreeVel(const double parent_data[], const double node[9], double
                      *collision_flag, creal_T *Tf, double xf[9])
{
  double dv0[9];
  creal_T Tf_data[8];
  int Tf_size[1];
  int nx;
  int k;
  double x_data[8];
  double y_data[8];
  boolean_T tmp_data[8];
  int i;
  int partialTrueCount;
  double cost;
  creal_T Tf_data_tmp;
  creal_T t[12];
  double costnow;
  double b_state[30];
  boolean_T exitg3;
  int exitg2;
  static const signed char iv0[3] = { 20, 10, 10 };

  int exitg1;
  static const signed char iv1[4] = { 6, 6, 12, 12 };

  static const signed char iv2[4] = { 0, 5, 5, 0 };

  static const signed char iv3[4] = { 0, 0, 0, 5 };

  static const signed char iv4[4] = { 10, 5, 10, 5 };

  *collision_flag = 0.0;
  quadf_timePFF(parent_data[0], parent_data[1], parent_data[2], parent_data[3],
                parent_data[4], parent_data[5], parent_data[6], parent_data[7],
                parent_data[8], node[0], node[1], node[2], dv0);
  roots(dv0, Tf_data, Tf_size);
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
  k = 0;
  for (i = 0; i <= nx; i++) {
    if (tmp_data[i]) {
      k++;
    }
  }

  partialTrueCount = 0;
  for (i = 0; i <= nx; i++) {
    if (tmp_data[i]) {
      Tf_data[partialTrueCount] = Tf_data[i];
      partialTrueCount++;
    }
  }

  nx = k - 1;
  k = 0;
  for (i = 0; i <= nx; i++) {
    if (Tf_data[i].re >= 0.0) {
      k++;
    }
  }

  partialTrueCount = 0;
  for (i = 0; i <= nx; i++) {
    if (Tf_data[i].re >= 0.0) {
      Tf_data[partialTrueCount] = Tf_data[i];
      partialTrueCount++;
    }
  }

  *Tf = Tf_data[0];
  cost = quadf_costPFF(Tf_data[0], parent_data[0], parent_data[1], parent_data[2],
                       parent_data[3], parent_data[4], parent_data[5],
                       parent_data[6], parent_data[7], parent_data[8], node[0],
                       node[1], node[2]);
  for (i = 0; i <= k - 2; i++) {
    Tf_data_tmp = Tf_data[1 + i];
    costnow = quadf_costPFF(Tf_data_tmp, parent_data[0], parent_data[1],
      parent_data[2], parent_data[3], parent_data[4], parent_data[5],
      parent_data[6], parent_data[7], parent_data[8], node[0], node[1], node[2]);
    if (costnow < cost) {
      *Tf = Tf_data_tmp;
      cost = costnow;
    }
  }

  /*      Tf = norm(parent(1 : 3) - node(1 : 3))/4; */
  linspace(*Tf, t);
  for (i = 0; i < 10; i++) {
    quadf_statePFF(t[i + 1], *Tf, parent_data[0], parent_data[1], parent_data[2],
                   parent_data[3], parent_data[4], parent_data[5], parent_data[6],
                   parent_data[7], parent_data[8], node[0], node[1], node[2],
                   *(double (*)[3])&b_state[3 * i]);
  }

  quadf_statePFFull(*Tf, *Tf, parent_data[0], parent_data[1], parent_data[2],
                    parent_data[3], parent_data[4], parent_data[5], parent_data
                    [6], parent_data[7], parent_data[8], node[0], node[1], node
                    [2], xf);
  nx = 0;
  exitg3 = false;
  while ((!exitg3) && (nx < 10)) {
    i = 0;
    do {
      exitg2 = 0;
      if (i < 3) {
        cost = b_state[i + 3 * nx];
        if ((cost > iv0[i]) || (cost < 0.0)) {
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
      exitg3 = true;
    } else {
      do {
        exitg1 = 0;
        if (i < 4) {
          cost = b_state[3 * nx];
          if ((cost >= iv1[i]) && (cost <= (double)iv1[i] + 2.0)) {
            cost = b_state[1 + 3 * nx];
            if ((cost >= iv2[i]) && (cost <= (double)iv2[i] + 5.0)) {
              cost = b_state[2 + 3 * nx];
              if ((cost >= iv3[i]) && (cost <= iv3[i] + iv4[i])) {
                *collision_flag = 1.0;
                exitg1 = 2;
              } else {
                i++;
              }
            } else {
              i++;
            }
          } else {
            i++;
          }
        } else {
          nx++;
          exitg1 = 1;
        }
      } while (exitg1 == 0);

      if (exitg1 == 1) {
      } else {
        exitg3 = true;
      }
    }
  }
}

/* End of code generation (collisionFreeVel.c) */
