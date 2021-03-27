/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * RRTstar3D.c
 *
 * Code generation for function 'RRTstar3D'
 *
 */

/* Include files */
#include <math.h>
#include <string.h>
#include "rt_nonfinite.h"
#include "benchmarkRRT.h"
#include "RRTstar3D.h"
#include "collisionKnowTf.h"
#include "segment_cost.h"
#include "benchmarkRRT_emxutil.h"
#include "collisionFreeVelKnowTf.h"
#include "segment_costFreeVel.h"
#include "sqr_eucl_dist.h"
#include "quadf_cost.h"
#include "quadf_costPFF.h"
#include "collisionFreeVel.h"
#include "norm.h"
#include "rand.h"
#include "collision.h"
#include "benchmarkRRT_rtwutil.h"

/* Function Definitions */
void extendTree(const emxArray_real_T *tree, double GChild[360000],
                emxArray_real_T *new_tree, double *flag)
{
  double new_node_data[14];
  int i3;
  int flag1;
  emxArray_real_T *tmp_dist;
  emxArray_real_T *sqr_dist;
  emxArray_boolean_T *x;
  emxArray_int32_T *ii;
  cell_wrap_2 reshapes[2];
  double randomPoint[9];
  int n;
  double Vect[3];
  double near_cost;
  creal_T unusedU5;
  int nx;
  int i;
  int new_node_idx;
  double cost_near;
  int iidx;
  boolean_T exitg1;
  int min_parent_idx;
  double d2;
  double new_point[9];
  double tree_data[13];
  static const signed char iv10[4] = { 6, 6, 12, 12 };

  static const signed char iv11[4] = { 0, 5, 5, 0 };

  static const signed char iv12[4] = { 0, 0, 0, 5 };

  static const signed char iv13[4] = { 10, 5, 10, 5 };

  double min_cost;
  int i4;
  double rnewcost;
  boolean_T b0;
  boolean_T bv0[600];
  short tmp_data[600];

  /*  segmentLength: maximum stepsize, r: neighbor radius */
  memset(&new_node_data[0], 0, 14U * sizeof(double));
  i3 = new_tree->size[0] * new_tree->size[1];
  new_tree->size[0] = 1;
  new_tree->size[1] = 13;
  emxEnsureCapacity_real_T(new_tree, i3);
  for (i3 = 0; i3 < 13; i3++) {
    new_tree->data[i3] = 0.0;
  }

  flag1 = 0;
  emxInit_real_T(&tmp_dist, 2);
  emxInit_real_T(&sqr_dist, 1);
  emxInit_boolean_T(&x, 1);
  emxInit_int32_T(&ii, 1);
  emxInitMatrix_cell_wrap_2(reshapes);
  while (flag1 == 0) {
    /*  select a random point */
    b_rand();
    memset(&randomPoint[0], 0, 9U * sizeof(double));
    randomPoint[0] = 20.0 * b_rand();
    randomPoint[1] = 10.0 * b_rand();
    randomPoint[2] = 10.0 * b_rand();

    /*  find node that is closest to randomPoint (Eucl. dist. between positions).  */
    /* fprintf("find node that is closest to randomPoint \n"); */
    /*  NOTE */
    /*      tmp = tree(:, 1 : dim) - randomPoint(1 : dim); */
    i3 = tmp_dist->size[0] * tmp_dist->size[1];
    tmp_dist->size[0] = tree->size[0];
    tmp_dist->size[1] = 3;
    emxEnsureCapacity_real_T(tmp_dist, i3);
    nx = tree->size[0] * 3;
    for (i3 = 0; i3 < nx; i3++) {
      tmp_dist->data[i3] = 0.0;
    }

    i3 = tree->size[0];
    for (i = 0; i < i3; i++) {
      tmp_dist->data[i] = tree->data[i] - randomPoint[0];
      tmp_dist->data[i + tmp_dist->size[0]] = tree->data[i + tree->size[0]] -
        randomPoint[1];
      tmp_dist->data[i + (tmp_dist->size[0] << 1)] = tree->data[i + (tree->size
        [0] << 1)] - randomPoint[2];
    }

    sqr_eucl_dist(tmp_dist, sqr_dist);
    n = sqr_dist->size[0];
    if (sqr_dist->size[0] <= 2) {
      if (sqr_dist->size[0] == 1) {
        cost_near = sqr_dist->data[0];
        iidx = 0;
      } else if ((sqr_dist->data[0] > sqr_dist->data[1]) || (rtIsNaN
                  (sqr_dist->data[0]) && (!rtIsNaN(sqr_dist->data[1])))) {
        cost_near = sqr_dist->data[1];
        iidx = 1;
      } else {
        cost_near = sqr_dist->data[0];
        iidx = 0;
      }
    } else {
      if (!rtIsNaN(sqr_dist->data[0])) {
        new_node_idx = 1;
      } else {
        new_node_idx = 0;
        nx = 2;
        exitg1 = false;
        while ((!exitg1) && (nx <= sqr_dist->size[0])) {
          if (!rtIsNaN(sqr_dist->data[nx - 1])) {
            new_node_idx = nx;
            exitg1 = true;
          } else {
            nx++;
          }
        }
      }

      if (new_node_idx == 0) {
        cost_near = sqr_dist->data[0];
        iidx = 0;
      } else {
        cost_near = sqr_dist->data[new_node_idx - 1];
        iidx = new_node_idx - 1;
        i3 = new_node_idx + 1;
        for (nx = i3; nx <= n; nx++) {
          if (cost_near > sqr_dist->data[nx - 1]) {
            cost_near = sqr_dist->data[nx - 1];
            iidx = nx - 1;
          }
        }
      }
    }

    min_parent_idx = iidx;
    Vect[0] = randomPoint[0] - tree->data[iidx];
    Vect[1] = randomPoint[1] - tree->data[iidx + tree->size[0]];
    Vect[2] = randomPoint[2] - tree->data[iidx + (tree->size[0] << 1)];
    d2 = b_norm(Vect);
    Vect[0] /= d2;
    Vect[1] /= d2;
    Vect[2] /= d2;

    /*  NOTE */
    memset(&new_point[0], 0, 9U * sizeof(double));

    /*  find new_point that is within the range of min_parent_idx in terms of segmentLength (Eucl. dist. between positions). */
    /* fprintf("find new_point that is within the range of min_parent_idx in terms of segmentLength (Eucl. dist. between positions)\n"); */
    if (cost_near > 16.0) {
      /*  generate a new point that is closest to randomPoint, segmentLength away from tree(idx,1:dim) */
      new_point[0] = tree->data[iidx] + Vect[0] * 4.0;
      new_point[1] = tree->data[iidx + tree->size[0]] + Vect[1] * 4.0;
      new_point[2] = tree->data[iidx + (tree->size[0] << 1)] + Vect[2] * 4.0;
    } else {
      new_point[0] = randomPoint[0];
      new_point[1] = randomPoint[1];
      new_point[2] = randomPoint[2];
    }

    /*  check if the new_point is in collision */
    /* fprintf("check if the new_point is in collision\n"); */
    /*  check if a point is in collision */
    n = 0;
    if ((new_point[0] > 20.0) || (new_point[0] < 0.0)) {
      n = 1;
    }

    if ((new_point[1] > 10.0) || (new_point[1] < 0.0)) {
      n = 1;
    }

    if ((new_point[2] > 10.0) || (new_point[2] < 0.0)) {
      n = 1;
    }

    if (n == 0) {
      /*  check each obstacle */
      i = 0;
      exitg1 = false;
      while ((!exitg1) && (i < 4)) {
        if ((new_point[0] >= iv10[i]) && (new_point[0] <= (double)iv10[i] + 2.0)
            && (new_point[1] >= iv11[i]) && (new_point[1] <= (double)iv11[i] +
             5.0) && (new_point[2] >= iv12[i]) && (new_point[2] <= iv12[i] +
             iv13[i])) {
          n = 1;
          exitg1 = true;
        } else {
          i++;
        }
      }
    }

    if (n == 0) {
      for (i3 = 0; i3 < 13; i3++) {
        tree_data[i3] = tree->data[iidx + tree->size[0] * i3];
      }

      collisionFreeVel(tree_data, new_point, &near_cost, &unusedU5, randomPoint);

      /*  this collision checking includes a steering function and forward simulation */
      if (near_cost == 0.0) {
        /* fprintf("no collision \n"); */
        for (i3 = 0; i3 < 6; i3++) {
          new_point[3 + i3] = randomPoint[3 + i3];
        }

        /* calculate the cost from root to to_point through parent from_node */
        min_cost = quadf_costPFF(unusedU5, tree->data[iidx], tree->data[iidx +
          tree->size[0]], tree->data[iidx + (tree->size[0] << 1)], tree->
          data[iidx + tree->size[0] * 3], tree->data[iidx + (tree->size[0] << 2)],
          tree->data[iidx + tree->size[0] * 5], tree->data[iidx + tree->size[0] *
          6], tree->data[iidx + tree->size[0] * 7], tree->data[iidx +
          (tree->size[0] << 3)], new_point[0], new_point[1], new_point[2]);
        min_cost += tree->data[iidx + tree->size[0] * 10];

        /*  total cost from root to new_point through its parent tree(idx,:) */
        /*        tmp_dist = tree(:, 1 : dim) - new_point(1 : dim); */
        i3 = tmp_dist->size[0] * tmp_dist->size[1];
        tmp_dist->size[0] = tree->size[0];
        tmp_dist->size[1] = 3;
        emxEnsureCapacity_real_T(tmp_dist, i3);
        nx = tree->size[0] * 3;
        for (i3 = 0; i3 < nx; i3++) {
          tmp_dist->data[i3] = 0.0;
        }

        i3 = tree->size[0];
        for (i = 0; i < i3; i++) {
          tmp_dist->data[i] = tree->data[i] - new_point[0];
          tmp_dist->data[i + tmp_dist->size[0]] = tree->data[i + tree->size[0]]
            - new_point[1];
          tmp_dist->data[i + (tmp_dist->size[0] << 1)] = tree->data[i +
            (tree->size[0] << 1)] - new_point[2];
        }

        /*  find near neighbors    */
        cost_near = fmin(40.0 * rt_powd_snf(log((double)tree->size[0] + 1.0) /
          (double)tree->size[0], 0.33333333333333331), 4.0);
        cost_near *= cost_near;
        sqr_eucl_dist(tmp_dist, sqr_dist);
        i3 = x->size[0];
        x->size[0] = sqr_dist->size[0];
        emxEnsureCapacity_boolean_T(x, i3);
        nx = sqr_dist->size[0];
        for (i3 = 0; i3 < nx; i3++) {
          x->data[i3] = (sqr_dist->data[i3] <= cost_near);
        }

        nx = x->size[0];
        new_node_idx = 0;
        i3 = ii->size[0];
        ii->size[0] = x->size[0];
        emxEnsureCapacity_int32_T(ii, i3);
        n = 0;
        exitg1 = false;
        while ((!exitg1) && (n <= nx - 1)) {
          if (x->data[n]) {
            new_node_idx++;
            ii->data[new_node_idx - 1] = n + 1;
            if (new_node_idx >= nx) {
              exitg1 = true;
            } else {
              n++;
            }
          } else {
            n++;
          }
        }

        if (x->size[0] == 1) {
          if (new_node_idx == 0) {
            ii->size[0] = 0;
          }
        } else if (1 > new_node_idx) {
          ii->size[0] = 0;
        } else {
          i3 = ii->size[0];
          ii->size[0] = new_node_idx;
          emxEnsureCapacity_int32_T(ii, i3);
        }

        i3 = sqr_dist->size[0];
        sqr_dist->size[0] = ii->size[0];
        emxEnsureCapacity_real_T(sqr_dist, i3);
        nx = ii->size[0];
        for (i3 = 0; i3 < nx; i3++) {
          sqr_dist->data[i3] = ii->data[i3];
        }

        /* fprintf("choosing new parent\n") */
        if (sqr_dist->size[0] >= 1) {
          i3 = sqr_dist->size[0];
          for (i = 0; i < i3; i++) {
            /*  choose parent node */
            if ((int)sqr_dist->data[i] != iidx + 1) {
              n = (int)sqr_dist->data[i];
              for (i4 = 0; i4 < 13; i4++) {
                tree_data[i4] = tree->data[(n + tree->size[0] * i4) - 1];
              }

              segment_costFreeVel(tree_data, new_point, &cost_near, &unusedU5);
              cost_near += tree->data[((int)sqr_dist->data[i] + tree->size[0] *
                10) - 1];
              if (cost_near + 0.001 < min_cost) {
                n = (int)sqr_dist->data[i];
                for (i4 = 0; i4 < 13; i4++) {
                  tree_data[i4] = tree->data[(n + tree->size[0] * i4) - 1];
                }

                collisionFreeVelKnowTf(tree_data, new_point, &unusedU5,
                  &near_cost, randomPoint);
                if (near_cost == 0.0) {
                  min_cost = cost_near;
                  min_parent_idx = (int)sqr_dist->data[i] - 1;
                  for (i4 = 0; i4 < 6; i4++) {
                    new_point[3 + i4] = randomPoint[3 + i4];
                  }
                }
              }
            }
          }
        }

        /* fprintf("assembling new node \n"); */
        memcpy(&tree_data[0], &new_point[0], 9U * sizeof(double));
        tree_data[9] = 0.0;
        tree_data[10] = min_cost;
        tree_data[11] = min_parent_idx + 1;
        tree_data[12] = 0.0;
        memcpy(&new_node_data[0], &tree_data[0], 13U * sizeof(double));
        if (tree->size[0] != 0) {
          n = tree->size[0];
        } else {
          n = 0;
        }

        i3 = reshapes[1].f1->size[0] * reshapes[1].f1->size[1];
        reshapes[1].f1->size[0] = 1;
        reshapes[1].f1->size[1] = 13;
        emxEnsureCapacity_real_T(reshapes[1].f1, i3);
        for (i3 = 0; i3 < 13; i3++) {
          reshapes[1].f1->data[i3] = new_node_data[i3];
        }

        i3 = new_tree->size[0] * new_tree->size[1];
        new_tree->size[0] = n + reshapes[1].f1->size[0];
        new_tree->size[1] = 13;
        emxEnsureCapacity_real_T(new_tree, i3);
        for (i3 = 0; i3 < 13; i3++) {
          for (i4 = 0; i4 < n; i4++) {
            new_tree->data[i4 + new_tree->size[0] * i3] = tree->data[i4 +
              tree->size[0] * i3];
          }
        }

        for (i3 = 0; i3 < 13; i3++) {
          nx = reshapes[1].f1->size[0];
          for (i4 = 0; i4 < nx; i4++) {
            new_tree->data[(i4 + n) + new_tree->size[0] * i3] = reshapes[1]
              .f1->data[i4 + reshapes[1].f1->size[0] * i3];
          }
        }

        new_node_idx = new_tree->size[0] - 1;
        new_tree->data[min_parent_idx + new_tree->size[0] * 12]++;

        /*  ChildNum + 1 */
        /* fprintf("updating GChild \n"); */
        d2 = new_node_idx + 1;
        GChild[min_parent_idx + 600 * ((int)new_tree->data[min_parent_idx +
          new_tree->size[0] * 12] - 1)] = d2;

        /*  update GChild matrix */
        /* fprintf("rewiring \n"); */
        if (sqr_dist->size[0] >= 1) {
          /*  rewire */
          i3 = sqr_dist->size[0];
          for (flag1 = 0; flag1 < i3; flag1++) {
            if ((int)sqr_dist->data[flag1] != min_parent_idx + 1) {
              near_cost = new_tree->data[((int)sqr_dist->data[flag1] +
                new_tree->size[0] * 10) - 1];
              n = (int)sqr_dist->data[flag1];
              for (i4 = 0; i4 < 13; i4++) {
                tree_data[i4] = new_tree->data[(n + new_tree->size[0] * i4) - 1];
              }

              segment_cost(new_point, tree_data, &cost_near, &unusedU5);
              rnewcost = min_cost + cost_near;
              if (near_cost > rnewcost + 0.001) {
                n = (int)sqr_dist->data[flag1];
                for (i4 = 0; i4 < 13; i4++) {
                  tree_data[i4] = new_tree->data[(n + new_tree->size[0] * i4) -
                    1];
                }

                near_cost = collisionKnowTf(new_node_data, tree_data, &unusedU5);
                if (near_cost == 0.0) {
                  cost_near = sqr_dist->data[flag1];
                  n = (int)new_tree->data[((int)sqr_dist->data[flag1] +
                    new_tree->size[0] * 11) - 1];
                  nx = 0;
                  for (i = 0; i < 600; i++) {
                    b0 = (GChild[(n + 600 * i) - 1] == cost_near);
                    bv0[i] = b0;
                    if (b0) {
                      nx++;
                    }
                  }

                  n = 0;
                  for (i = 0; i < 600; i++) {
                    if (bv0[i]) {
                      tmp_data[n] = (short)(i + 1);
                      n++;
                    }
                  }

                  n = (int)new_tree->data[((int)sqr_dist->data[flag1] +
                    new_tree->size[0] * 11) - 1];
                  for (i4 = 0; i4 < nx; i4++) {
                    GChild[(n + 600 * (tmp_data[i4] - 1)) - 1] = -1.0;
                  }

                  /*  parent of reduced_idx(j) before rewire, change its child list. */
                  new_tree->data[((int)sqr_dist->data[flag1] + new_tree->size[0]
                                  * 10) - 1] = rnewcost;

                  /*  update the cost and parent information of the node being rewired, reduced_idx(j) */
                  new_tree->data[((int)sqr_dist->data[flag1] + new_tree->size[0]
                                  * 11) - 1] = d2;
                  new_tree->data[new_node_idx + new_tree->size[0] * 12]++;

                  /*  add the node being rewired to the child list of the new added node, new_node_idx. */
                  GChild[new_node_idx + 600 * ((int)new_tree->data[new_node_idx
                    + new_tree->size[0] * 12] - 1)] = sqr_dist->data[flag1];

                  /*  update all cost of the descendant of the node being rewired */
                }
              }
            }
          }
        }

        /* fprintf("rewiring done \n"); */
        flag1 = 1;
      }
    }
  }

  emxFreeMatrix_cell_wrap_2(reshapes);
  emxFree_int32_T(&ii);
  emxFree_boolean_T(&x);
  emxFree_real_T(&sqr_dist);
  emxFree_real_T(&tmp_dist);
  n = 0;

  /*  check to see if new node connects directly to end_node */
  Vect[0] = new_node_data[0] - 18.0;
  Vect[1] = new_node_data[1] - 8.0;
  Vect[2] = new_node_data[2] - 8.0;
  if (b_norm(Vect) < 4.0) {
    collision(new_node_data, &near_cost, &unusedU5);
    if (near_cost == 0.0) {
      n = 1;
      new_tree->data[(new_tree->size[0] + new_tree->size[0] * 9) - 1] = 1.0;

      /*  mark node as connecting to end. */
    }
  }

  *flag = n;
}

/* End of code generation (RRTstar3D.c) */
