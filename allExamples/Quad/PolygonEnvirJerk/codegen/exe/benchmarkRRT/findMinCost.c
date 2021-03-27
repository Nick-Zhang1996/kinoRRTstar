/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * findMinCost.c
 *
 * Code generation for function 'findMinCost'
 *
 */

/* Include files */
#include <math.h>
#include "rt_nonfinite.h"
#include "benchmarkRRT.h"
#include "findMinCost.h"
#include "benchmarkRRT_emxutil.h"
#include "quadf_cost.h"
#include "roots.h"
#include "quadf_time.h"

/* Function Definitions */
double findMinCost(const emxArray_real_T *tree)
{
  double cost;
  emxArray_boolean_T *idx;
  int nx;
  int i2;
  int n;
  int i;
  emxArray_int32_T *r0;
  int partialTrueCount;
  emxArray_real_T *connectingNodes;
  emxArray_real_T *PathCost;
  int k;
  double dv4[7];
  double from_node_data[13];
  creal_T Tf_data[6];
  int Tf_size[1];
  double x_data[8];
  boolean_T exitg1;
  double y_data[8];
  boolean_T tmp_data[6];
  creal_T tf;
  double lcost;
  creal_T Tf_data_tmp;
  double costnow;
  emxInit_boolean_T(&idx, 1);
  nx = tree->size[0];
  i2 = idx->size[0];
  idx->size[0] = nx;
  emxEnsureCapacity_boolean_T(idx, i2);
  for (i2 = 0; i2 < nx; i2++) {
    idx->data[i2] = (tree->data[i2 + tree->size[0] * 9] == 1.0);
  }

  nx = idx->size[0] - 1;
  n = 0;
  for (i = 0; i <= nx; i++) {
    if (idx->data[i]) {
      n++;
    }
  }

  emxInit_int32_T(&r0, 1);
  i2 = r0->size[0];
  r0->size[0] = n;
  emxEnsureCapacity_int32_T(r0, i2);
  partialTrueCount = 0;
  for (i = 0; i <= nx; i++) {
    if (idx->data[i]) {
      r0->data[partialTrueCount] = i + 1;
      partialTrueCount++;
    }
  }

  emxFree_boolean_T(&idx);
  emxInit_real_T(&connectingNodes, 2);
  i2 = connectingNodes->size[0] * connectingNodes->size[1];
  connectingNodes->size[0] = r0->size[0];
  connectingNodes->size[1] = 13;
  emxEnsureCapacity_real_T(connectingNodes, i2);
  for (i2 = 0; i2 < 13; i2++) {
    nx = r0->size[0];
    for (n = 0; n < nx; n++) {
      connectingNodes->data[n + connectingNodes->size[0] * i2] = tree->data
        [(r0->data[n] + tree->size[0] * i2) - 1];
    }
  }

  emxInit_real_T(&PathCost, 1);
  i2 = r0->size[0];
  n = PathCost->size[0];
  PathCost->size[0] = r0->size[0];
  emxEnsureCapacity_real_T(PathCost, n);
  for (k = 0; k < i2; k++) {
    for (n = 0; n < 13; n++) {
      from_node_data[n] = connectingNodes->data[k + connectingNodes->size[0] * n];
    }

    quadf_time(connectingNodes->data[k], connectingNodes->data[k +
               connectingNodes->size[0]], connectingNodes->data[k +
               (connectingNodes->size[0] << 1)], connectingNodes->data[k +
               connectingNodes->size[0] * 3], connectingNodes->data[k +
               (connectingNodes->size[0] << 2)], connectingNodes->data[k +
               connectingNodes->size[0] * 5], connectingNodes->data[k +
               connectingNodes->size[0] * 6], connectingNodes->data[k +
               connectingNodes->size[0] * 7], connectingNodes->data[k +
               (connectingNodes->size[0] << 3)], 18.0, 8.0, 8.0, 0.0, 0.0, 0.0,
               0.0, 0.0, 0.0, dv4);
    b_roots(dv4, Tf_data, Tf_size);
    nx = Tf_size[0];
    for (n = 0; n < nx; n++) {
      x_data[n] = Tf_data[n].im;
    }

    nx = Tf_size[0];
    for (n = 0; n < nx; n++) {
      y_data[n] = fabs(x_data[n]);
    }

    nx = (signed char)Tf_size[0];
    for (n = 0; n < nx; n++) {
      tmp_data[n] = (y_data[n] < 0.0001);
    }

    nx = (signed char)Tf_size[0] - 1;
    n = 0;
    for (i = 0; i <= nx; i++) {
      if (tmp_data[i]) {
        n++;
      }
    }

    partialTrueCount = 0;
    for (i = 0; i <= nx; i++) {
      if (tmp_data[i]) {
        Tf_data[partialTrueCount] = Tf_data[i];
        partialTrueCount++;
      }
    }

    nx = n - 1;
    n = 0;
    for (i = 0; i <= nx; i++) {
      if (Tf_data[i].re >= 0.0) {
        n++;
      }
    }

    partialTrueCount = 0;
    for (i = 0; i <= nx; i++) {
      if (Tf_data[i].re >= 0.0) {
        Tf_data[partialTrueCount] = Tf_data[i];
        partialTrueCount++;
      }
    }

    tf = Tf_data[0];
    lcost = quadf_cost(Tf_data[0], connectingNodes->data[k],
                       connectingNodes->data[k + connectingNodes->size[0]],
                       connectingNodes->data[k + (connectingNodes->size[0] << 1)],
                       connectingNodes->data[k + connectingNodes->size[0] * 3],
                       connectingNodes->data[k + (connectingNodes->size[0] << 2)],
                       connectingNodes->data[k + connectingNodes->size[0] * 5],
                       connectingNodes->data[k + connectingNodes->size[0] * 6],
                       connectingNodes->data[k + connectingNodes->size[0] * 7],
                       connectingNodes->data[k + (connectingNodes->size[0] << 3)],
                       18.0, 8.0, 8.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
    for (i = 0; i <= n - 2; i++) {
      Tf_data_tmp = Tf_data[1 + i];
      costnow = quadf_cost(Tf_data_tmp, from_node_data[0], from_node_data[1],
                           from_node_data[2], from_node_data[3], from_node_data
                           [4], from_node_data[5], from_node_data[6],
                           from_node_data[7], from_node_data[8], 18.0, 8.0, 8.0,
                           0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
      if (costnow < lcost) {
        tf = Tf_data_tmp;
        lcost = costnow;
      }
    }

    /*      Tf = norm(from_node(1 : 3) - to_point(1 : 3))/4; */
    lcost = quadf_cost(tf, connectingNodes->data[k], connectingNodes->data[k +
                       connectingNodes->size[0]], connectingNodes->data[k +
                       (connectingNodes->size[0] << 1)], connectingNodes->data[k
                       + connectingNodes->size[0] * 3], connectingNodes->data[k
                       + (connectingNodes->size[0] << 2)], connectingNodes->
                       data[k + connectingNodes->size[0] * 5],
                       connectingNodes->data[k + connectingNodes->size[0] * 6],
                       connectingNodes->data[k + connectingNodes->size[0] * 7],
                       connectingNodes->data[k + (connectingNodes->size[0] << 3)],
                       18.0, 8.0, 8.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
    PathCost->data[k] = tree->data[(r0->data[k] + tree->size[0] * 10) - 1] +
      lcost;
  }

  emxFree_int32_T(&r0);
  emxFree_real_T(&connectingNodes);

  /*  find minimum cost last node */
  n = PathCost->size[0];
  if (PathCost->size[0] <= 2) {
    if (PathCost->size[0] == 1) {
      cost = PathCost->data[0];
    } else if ((PathCost->data[0] > PathCost->data[1]) || (rtIsNaN
                (PathCost->data[0]) && (!rtIsNaN(PathCost->data[1])))) {
      cost = PathCost->data[1];
    } else {
      cost = PathCost->data[0];
    }
  } else {
    if (!rtIsNaN(PathCost->data[0])) {
      nx = 1;
    } else {
      nx = 0;
      k = 2;
      exitg1 = false;
      while ((!exitg1) && (k <= PathCost->size[0])) {
        if (!rtIsNaN(PathCost->data[k - 1])) {
          nx = k;
          exitg1 = true;
        } else {
          k++;
        }
      }
    }

    if (nx == 0) {
      cost = PathCost->data[0];
    } else {
      cost = PathCost->data[nx - 1];
      i2 = nx + 1;
      for (k = i2; k <= n; k++) {
        if (cost > PathCost->data[k - 1]) {
          cost = PathCost->data[k - 1];
        }
      }
    }
  }

  emxFree_real_T(&PathCost);
  return cost;
}

/* End of code generation (findMinCost.c) */
