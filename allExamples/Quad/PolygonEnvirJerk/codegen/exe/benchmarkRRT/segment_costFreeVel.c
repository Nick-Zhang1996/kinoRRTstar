/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * segment_costFreeVel.c
 *
 * Code generation for function 'segment_costFreeVel'
 *
 */

/* Include files */
#include <math.h>
#include "rt_nonfinite.h"
#include "benchmarkRRT.h"
#include "segment_costFreeVel.h"
#include "quadf_costPFF.h"
#include "roots.h"
#include "quadf_timePFF.h"

/* Function Definitions */
void segment_costFreeVel(const double from_node_data[], const double to_point[9],
  double *cost, creal_T *Tf)
{
  double dv1[9];
  creal_T Tf_data[8];
  int Tf_size[1];
  int nx;
  int k;
  double x_data[8];
  double y_data[8];
  boolean_T tmp_data[8];
  int i;
  int partialTrueCount;
  creal_T Tf_data_tmp;
  double costnow;
  quadf_timePFF(from_node_data[0], from_node_data[1], from_node_data[2],
                from_node_data[3], from_node_data[4], from_node_data[5],
                from_node_data[6], from_node_data[7], from_node_data[8],
                to_point[0], to_point[1], to_point[2], dv1);
  roots(dv1, Tf_data, Tf_size);
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
  *cost = quadf_costPFF(Tf_data[0], from_node_data[0], from_node_data[1],
                        from_node_data[2], from_node_data[3], from_node_data[4],
                        from_node_data[5], from_node_data[6], from_node_data[7],
                        from_node_data[8], to_point[0], to_point[1], to_point[2]);
  for (i = 0; i <= k - 2; i++) {
    Tf_data_tmp = Tf_data[1 + i];
    costnow = quadf_costPFF(Tf_data_tmp, from_node_data[0], from_node_data[1],
      from_node_data[2], from_node_data[3], from_node_data[4], from_node_data[5],
      from_node_data[6], from_node_data[7], from_node_data[8], to_point[0],
      to_point[1], to_point[2]);
    if (costnow < *cost) {
      *Tf = Tf_data_tmp;
      *cost = costnow;
    }
  }

  /*      Tf = norm(from_node(1 : 3) - to_point(1 : 3))/4; */
  *cost = quadf_costPFF(*Tf, from_node_data[0], from_node_data[1],
                        from_node_data[2], from_node_data[3], from_node_data[4],
                        from_node_data[5], from_node_data[6], from_node_data[7],
                        from_node_data[8], to_point[0], to_point[1], to_point[2]);
}

/* End of code generation (segment_costFreeVel.c) */
