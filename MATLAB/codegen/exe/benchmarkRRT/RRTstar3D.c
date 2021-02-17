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
#include "sqr_eucl_dist.h"
#include "cost_eval.h"
#include "collision.h"
#include "rand.h"
#include "norm.h"
#include "benchmarkRRT_data.h"

/* Function Definitions */
void extendTree(const emxArray_real_T *tree, double GChild[1200000],
                emxArray_real_T *new_tree, double *flag)
{
  int i3;
  double new_node_data[8];
  int flag1;
  emxArray_real_T *tmp;
  emxArray_real_T *sqr_dist;
  emxArray_boolean_T *x;
  emxArray_int32_T *ii;
  cell_wrap_3 reshapes[2];
  double randomPoint_idx_0;
  double randomPoint_idx_1;
  int n;
  int nx;
  double Vect[2];
  creal_T cost_near;
  int new_node_idx;
  double ex;
  int iidx;
  boolean_T exitg1;
  int min_parent_idx;
  double d2;
  double d3;
  double new_point[4];
  int i;
  double tree_data[8];
  creal_T Tf;
  static const signed char iv5[5] = { 6, 15, 10, 5, 15 };

  static const signed char iv6[5] = { 6, 14, 11, 15, 5 };

  creal_T min_cost;
  creal_T new_node[8];
  int i4;
  boolean_T b0;
  boolean_T bv0[300];
  short tmp_data[300];

  /*  segmentLength: maximum stepsize, r: neighbor radius */
  /*  NOTE */
  for (i3 = 0; i3 < 5; i3++) {
    new_node_data[i3] = 0.0;
  }

  i3 = new_tree->size[0] * new_tree->size[1];
  new_tree->size[0] = 1;
  new_tree->size[1] = 8;
  emxEnsureCapacity_real_T(new_tree, i3);
  for (i3 = 0; i3 < 8; i3++) {
    new_tree->data[i3] = 0.0;
  }

  flag1 = 0;
  emxInit_real_T(&tmp, 2);
  emxInit_real_T(&sqr_dist, 1);
  emxInit_boolean_T(&x, 1);
  emxInit_int32_T(&ii, 1);
  emxInitMatrix_cell_wrap_3(reshapes);
  while (flag1 == 0) {
    /*  select a random point */
    b_rand();
    randomPoint_idx_0 = 20.0 * b_rand();
    randomPoint_idx_1 = 20.0 * b_rand();

    /*  find node that is closest to randomPoint (Eucl. dist. between positions).  */
    /* tmp = tree(:, 1 : dim) - randomPoint(1 : dim); */
    nx = tree->size[0];
    i3 = tmp->size[0] * tmp->size[1];
    tmp->size[0] = nx;
    tmp->size[1] = 2;
    emxEnsureCapacity_real_T(tmp, i3);
    for (i3 = 0; i3 < nx; i3++) {
      tmp->data[i3] = tree->data[i3];
    }

    for (i3 = 0; i3 < nx; i3++) {
      tmp->data[i3 + tmp->size[0]] = tree->data[i3 + tree->size[0]];
    }

    i3 = tree->size[0];
    for (n = 0; n < i3; n++) {
      tmp->data[n] -= randomPoint_idx_0;
      tmp->data[n + tmp->size[0]] -= randomPoint_idx_1;
    }

    sqr_eucl_dist(tmp, sqr_dist);
    n = sqr_dist->size[0];
    if (sqr_dist->size[0] <= 2) {
      if (sqr_dist->size[0] == 1) {
        ex = sqr_dist->data[0];
        iidx = 0;
      } else if ((sqr_dist->data[0] > sqr_dist->data[1]) || (rtIsNaN
                  (sqr_dist->data[0]) && (!rtIsNaN(sqr_dist->data[1])))) {
        ex = sqr_dist->data[1];
        iidx = 1;
      } else {
        ex = sqr_dist->data[0];
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
        ex = sqr_dist->data[0];
        iidx = 0;
      } else {
        ex = sqr_dist->data[new_node_idx - 1];
        iidx = new_node_idx - 1;
        i3 = new_node_idx + 1;
        for (nx = i3; nx <= n; nx++) {
          if (ex > sqr_dist->data[nx - 1]) {
            ex = sqr_dist->data[nx - 1];
            iidx = nx - 1;
          }
        }
      }
    }

    min_parent_idx = iidx;
    Vect[0] = randomPoint_idx_0 - tree->data[iidx];
    Vect[1] = randomPoint_idx_1 - tree->data[iidx + tree->size[0]];
    d2 = b_norm(Vect);
    Vect[0] /= d2;
    Vect[1] /= d2;
    d2 = b_rand();
    d3 = b_rand();
    new_point[2] = -1.0 + 2.0 * d2;
    new_point[3] = -1.0 + 2.0 * d3;

    /*  vmin = -1£¬ vmax = 1 */
    /*          if (size(new_point(dim + 1 : 2 * dim),1) ~= 1) */
    /*              display("xx") */
    /*          end */
    /*  find new_point that is within the range of min_parent_idx in terms of segmentLength (Eucl. dist. between positions). */
    if (ex > 16.0) {
      /*  generate a new point that is closest to randomPoint, segmentLength away from tree(idx,1:dim) */
      new_point[0] = tree->data[iidx] + Vect[0] * 4.0;
      new_point[1] = tree->data[iidx + tree->size[0]] + Vect[1] * 4.0;
    } else {
      new_point[0] = randomPoint_idx_0;
      new_point[1] = randomPoint_idx_1;
    }

    /*  check if the new_point is in collision */
    /*  check if a point is in collision */
    n = 0;
    if ((new_point[0] > 20.0) || (new_point[0] < 0.0)) {
      n = 1;
    }

    if ((new_point[1] > 20.0) || (new_point[1] < 0.0)) {
      n = 1;
    }

    if (n == 0) {
      /*  check each obstacle */
      i = 0;
      exitg1 = false;
      while ((!exitg1) && (i < 5)) {
        d2 = new_point[0] - (double)iv5[i];
        d2 *= d2;
        Vect[0] = d2;
        d2 = new_point[1] - (double)iv6[i];
        d2 *= d2;
        if (sqrt(Vect[0] + d2) < dv0[i] + 0.1) {
          /*  (norm([p(1);p(2)]-[world.cx(i); world.cy(i)])<world.radius(i)+0.1) */
          n = 1;
          exitg1 = true;
        } else {
          i++;
        }
      }
    }

    if (n == 0) {
      for (i3 = 0; i3 < 8; i3++) {
        tree_data[i3] = tree->data[iidx + tree->size[0] * i3];
      }

      collision(tree_data, new_point, &randomPoint_idx_0, &Tf);

      /*  this collision checking includes a steering function */
      if (randomPoint_idx_0 == 0.0) {
        /* calculate the cost from root to to_point through parent from_node */
        /*      if ~exist('Tf','var') || isempty(Tf)                 */
        /*          Tf = roots(obj.eval_arrival_internal(from_node(1), from_node(2), from_node(3), from_node(4), to_point(1), to_point(2), to_point(3), to_point(4))); */
        /*          Tf = Tf(imag(Tf) == 0); */
        /*          Tf = min(Tf(Tf >= 0)); */
        /*      end */
        min_cost = cost_eval(Tf, tree->data[iidx], tree->data[iidx + tree->size
                             [0]], tree->data[iidx + (tree->size[0] << 1)],
                             tree->data[iidx + tree->size[0] * 3], new_point[0],
                             new_point[1], new_point[2], new_point[3]);
        min_cost.re += tree->data[iidx + tree->size[0] * 5];

        /*  total cost from root to new_point through its parent tree(idx,:) */
        new_node[0].re = new_point[0];
        new_node[0].im = 0.0;
        new_node[1].re = new_point[1];
        new_node[1].im = 0.0;
        new_node[2].re = new_point[2];
        new_node[2].im = 0.0;
        new_node[3].re = new_point[3];
        new_node[3].im = 0.0;
        new_node[4].re = 0.0;
        new_node[4].im = 0.0;
        new_node[5] = min_cost;
        new_node[6].re = iidx + 1;
        new_node[6].im = 0.0;
        new_node[7].re = 0.0;
        new_node[7].im = 0.0;

        /*  new node candidate */
        /* tmp_dist = tree(:, 1 : dim) - new_point(1 : dim); */
        nx = tree->size[0];
        i3 = tmp->size[0] * tmp->size[1];
        tmp->size[0] = nx;
        tmp->size[1] = 2;
        emxEnsureCapacity_real_T(tmp, i3);
        for (i3 = 0; i3 < nx; i3++) {
          tmp->data[i3] = tree->data[i3];
        }

        for (i3 = 0; i3 < nx; i3++) {
          tmp->data[i3 + tmp->size[0]] = tree->data[i3 + tree->size[0]];
        }

        i3 = tree->size[0];
        for (n = 0; n < i3; n++) {
          tmp->data[n] -= new_point[0];
          tmp->data[n + tmp->size[0]] -= new_point[1];
        }

        /*  find near neighbors    */
        randomPoint_idx_0 = fmin(40.0 * sqrt(log((double)tree->size[0] + 1.0) /
          (double)tree->size[0]), 3.0);
        randomPoint_idx_0 *= randomPoint_idx_0;
        sqr_eucl_dist(tmp, sqr_dist);
        i3 = x->size[0];
        x->size[0] = sqr_dist->size[0];
        emxEnsureCapacity_boolean_T(x, i3);
        nx = sqr_dist->size[0];
        for (i3 = 0; i3 < nx; i3++) {
          x->data[i3] = (sqr_dist->data[i3] <= randomPoint_idx_0);
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

        if (sqr_dist->size[0] >= 1) {
          i3 = sqr_dist->size[0];
          for (i = 0; i < i3; i++) {
            /*  choose parent node */
            if ((int)sqr_dist->data[i] != iidx + 1) {
              n = (int)sqr_dist->data[i];
              for (i4 = 0; i4 < 8; i4++) {
                tree_data[i4] = tree->data[(n + tree->size[0] * i4) - 1];
              }

              segment_cost(tree_data, new_point, &cost_near, &Tf);
              cost_near.re += tree->data[((int)sqr_dist->data[i] + tree->size[0]
                * 5) - 1];
              if (cost_near.re + 0.01 < min_cost.re) {
                n = (int)sqr_dist->data[i];
                for (i4 = 0; i4 < 8; i4++) {
                  tree_data[i4] = tree->data[(n + tree->size[0] * i4) - 1];
                }

                randomPoint_idx_0 = collisionKnowTf(tree_data, new_node, &Tf);
                if (randomPoint_idx_0 == 0.0) {
                  min_cost = cost_near;
                  min_parent_idx = (int)sqr_dist->data[i] - 1;
                }
              }
            }
          }
        }

        tree_data[0] = new_point[0];
        tree_data[1] = new_point[1];
        tree_data[2] = new_point[2];
        tree_data[3] = new_point[3];
        tree_data[4] = 0.0;
        tree_data[5] = min_cost.re;
        tree_data[6] = min_parent_idx + 1;
        tree_data[7] = 0.0;
        memcpy(&new_node_data[0], &tree_data[0], sizeof(double) << 3);
        if (tree->size[0] != 0) {
          n = tree->size[0];
        } else {
          n = 0;
        }

        i3 = reshapes[1].f1->size[0] * reshapes[1].f1->size[1];
        reshapes[1].f1->size[0] = 1;
        reshapes[1].f1->size[1] = 8;
        emxEnsureCapacity_real_T(reshapes[1].f1, i3);
        for (i3 = 0; i3 < 8; i3++) {
          reshapes[1].f1->data[i3] = new_node_data[i3];
        }

        i3 = new_tree->size[0] * new_tree->size[1];
        new_tree->size[0] = n + reshapes[1].f1->size[0];
        new_tree->size[1] = 8;
        emxEnsureCapacity_real_T(new_tree, i3);
        for (i3 = 0; i3 < 8; i3++) {
          for (i4 = 0; i4 < n; i4++) {
            new_tree->data[i4 + new_tree->size[0] * i3] = tree->data[i4 +
              tree->size[0] * i3];
          }
        }

        for (i3 = 0; i3 < 8; i3++) {
          nx = reshapes[1].f1->size[0];
          for (i4 = 0; i4 < nx; i4++) {
            new_tree->data[(i4 + n) + new_tree->size[0] * i3] = reshapes[1]
              .f1->data[i4 + reshapes[1].f1->size[0] * i3];
          }
        }

        /*  DEBUG check new node cost */
        /* fprintf("new node cost %.2f",real(min_cost)); */
        new_node_idx = new_tree->size[0] - 1;
        new_tree->data[min_parent_idx + new_tree->size[0] * 7]++;

        /*  ChildNum + 1 */
        d2 = new_node_idx + 1;
        GChild[min_parent_idx + 4000 * ((int)new_tree->data[min_parent_idx +
          new_tree->size[0] * 7] - 1)] = d2;

        /*  update GChild matrix */
        if (sqr_dist->size[0] >= 1) {
          /*  rewire */
          i3 = sqr_dist->size[0];
          for (flag1 = 0; flag1 < i3; flag1++) {
            if ((int)sqr_dist->data[flag1] != min_parent_idx + 1) {
              randomPoint_idx_0 = new_tree->data[((int)sqr_dist->data[flag1] +
                new_tree->size[0] * 5) - 1];
              n = (int)sqr_dist->data[flag1];
              for (i4 = 0; i4 < 8; i4++) {
                tree_data[i4] = new_tree->data[(n + new_tree->size[0] * i4) - 1];
              }

              b_segment_cost(new_point, tree_data, &cost_near, &Tf);
              cost_near.re += min_cost.re;
              if (randomPoint_idx_0 > cost_near.re + 0.01) {
                n = (int)sqr_dist->data[flag1];
                for (i4 = 0; i4 < 8; i4++) {
                  tree_data[i4] = new_tree->data[(n + new_tree->size[0] * i4) -
                    1];
                }

                randomPoint_idx_0 = b_collisionKnowTf(new_node_data, tree_data,
                  &Tf);
                if (randomPoint_idx_0 == 0.0) {
                  randomPoint_idx_0 = sqr_dist->data[flag1];
                  n = (int)new_tree->data[((int)sqr_dist->data[flag1] +
                    new_tree->size[0] * 6) - 1];
                  nx = 0;
                  for (i = 0; i < 300; i++) {
                    b0 = (GChild[(n + 4000 * i) - 1] == randomPoint_idx_0);
                    bv0[i] = b0;
                    if (b0) {
                      nx++;
                    }
                  }

                  n = 0;
                  for (i = 0; i < 300; i++) {
                    if (bv0[i]) {
                      tmp_data[n] = (short)(i + 1);
                      n++;
                    }
                  }

                  n = (int)new_tree->data[((int)sqr_dist->data[flag1] +
                    new_tree->size[0] * 6) - 1];
                  for (i4 = 0; i4 < nx; i4++) {
                    GChild[(n + 4000 * (tmp_data[i4] - 1)) - 1] = -1.0;
                  }

                  /*  parent of reduced_idx(j) before rewire, change its child list. */
                  new_tree->data[((int)sqr_dist->data[flag1] + new_tree->size[0]
                                  * 5) - 1] = cost_near.re;

                  /*  update the cost and parent information of the node being rewired, reduced_idx(j) */
                  new_tree->data[((int)sqr_dist->data[flag1] + new_tree->size[0]
                                  * 6) - 1] = d2;
                  new_tree->data[new_node_idx + new_tree->size[0] * 7]++;

                  /*  add the node being rewired to the child list of the new added node, new_node_idx. */
                  GChild[new_node_idx + 4000 * ((int)new_tree->data[new_node_idx
                    + new_tree->size[0] * 7] - 1)] = sqr_dist->data[flag1];

                  /*  update all cost of the descendant of the node being rewired */
                }
              }
            }
          }
        }

        flag1 = 1;
      }
    }
  }

  emxFreeMatrix_cell_wrap_3(reshapes);
  emxFree_int32_T(&ii);
  emxFree_boolean_T(&x);
  emxFree_real_T(&sqr_dist);
  emxFree_real_T(&tmp);
  n = 0;

  /*  check to see if new node connects directly to end_node */
  Vect[0] = new_node_data[0] - 18.0;
  Vect[1] = new_node_data[1] - 18.0;
  if (b_norm(Vect) < 0.2) {
    n = 1;
    new_tree->data[(new_tree->size[0] + (new_tree->size[0] << 2)) - 1] = 1.0;

    /*  mark node as connecting to end. */
  } else {
    Vect[0] = new_node_data[0] - 18.0;
    Vect[1] = new_node_data[1] - 18.0;
    if (b_norm(Vect) < 4.0) {
      b_collision(new_node_data, &randomPoint_idx_0, &cost_near);
      if (randomPoint_idx_0 == 0.0) {
        n = 1;
        new_tree->data[(new_tree->size[0] + (new_tree->size[0] << 2)) - 1] = 1.0;

        /*  mark node as connecting to end. */
      }
    }
  }

  *flag = n;
}

/* End of code generation (RRTstar3D.c) */
