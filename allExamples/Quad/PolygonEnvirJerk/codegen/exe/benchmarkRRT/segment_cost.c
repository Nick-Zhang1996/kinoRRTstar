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
#include "quadf_cost.h"
#include "roots.h"
#include "quadf_time.h"

/* Function Definitions */
void segment_cost(const double from_node[9], const double to_point_data[],
                  double *cost, creal_T *Tf)
{
  double dv2[7];
  creal_T Tf_data[6];
  int Tf_size[1];
  int nx;
  int k;
  double x_data[8];
  double y_data[8];
  boolean_T tmp_data[6];
  int i;
  int partialTrueCount;
  creal_T Tf_data_tmp;
  double costnow;
  quadf_time(from_node[0], from_node[1], from_node[2], from_node[3], from_node[4],
             from_node[5], from_node[6], from_node[7], from_node[8],
             to_point_data[0], to_point_data[1], to_point_data[2],
             to_point_data[3], to_point_data[4], to_point_data[5],
             to_point_data[6], to_point_data[7], to_point_data[8], dv2);
  b_roots(dv2, Tf_data, Tf_size);
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
  *cost = quadf_cost(Tf_data[0], from_node[0], from_node[1], from_node[2],
                     from_node[3], from_node[4], from_node[5], from_node[6],
                     from_node[7], from_node[8], to_point_data[0],
                     to_point_data[1], to_point_data[2], to_point_data[3],
                     to_point_data[4], to_point_data[5], to_point_data[6],
                     to_point_data[7], to_point_data[8]);
  for (i = 0; i <= k - 2; i++) {
    Tf_data_tmp = Tf_data[1 + i];
    costnow = quadf_cost(Tf_data_tmp, from_node[0], from_node[1], from_node[2],
                         from_node[3], from_node[4], from_node[5], from_node[6],
                         from_node[7], from_node[8], to_point_data[0],
                         to_point_data[1], to_point_data[2], to_point_data[3],
                         to_point_data[4], to_point_data[5], to_point_data[6],
                         to_point_data[7], to_point_data[8]);
    if (costnow < *cost) {
      *Tf = Tf_data_tmp;
      *cost = costnow;
    }
  }

  /*      Tf = norm(from_node(1 : 3) - to_point(1 : 3))/4; */
  *cost = quadf_cost(*Tf, from_node[0], from_node[1], from_node[2], from_node[3],
                     from_node[4], from_node[5], from_node[6], from_node[7],
                     from_node[8], to_point_data[0], to_point_data[1],
                     to_point_data[2], to_point_data[3], to_point_data[4],
                     to_point_data[5], to_point_data[6], to_point_data[7],
                     to_point_data[8]);
}

/* End of code generation (segment_cost.c) */
