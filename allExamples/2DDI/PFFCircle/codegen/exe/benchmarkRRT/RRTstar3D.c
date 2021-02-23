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
#include "DI_costFreeVel.h"
#include "collisionFreeVel.h"
#include "norm.h"
#include "rand.h"
#include "collision.h"
#include "benchmarkRRT_data.h"

/* Function Definitions */
void extendTree(const emxArray_real_T *tree, double GChild[16000000],
                emxArray_real_T *new_tree, double *flag)
{
  double new_node[8];
  int i3;
  int flag1;
  emxArray_real_T *tmp_dist;
  emxArray_real_T *sqr_dist;
  emxArray_boolean_T *x;
  emxArray_int32_T *ii;
  double randomPoint[4];
  int n;
  double Vect[2];
  double near_cost;
  creal_T unusedU4;
  int i;
  int new_node_idx;
  double cost_near;
  int nx;
  int iidx;
  boolean_T exitg1;
  int min_parent_idx;
  double d2;
  double new_point[4];
  double tree_data[8];
  static const signed char iv5[5] = { 6, 15, 10, 5, 15 };

  static const signed char iv6[5] = { 6, 14, 11, 15, 5 };

  double min_cost;
  int i4;
  double rnewcost;
  boolean_T b0;
  boolean_T bv0[4000];
  short tmp_data[4000];

  /*  segmentLength: maximum stepsize, r: neighbor radius */
  /*  NOTE */
  /* fprintf("sampling node %.0f \n",size(tree,1)+1); */
  memset(&new_node[0], 0, sizeof(double) << 3);
  i3 = new_tree->size[0] * new_tree->size[1];
  new_tree->size[0] = 1;
  new_tree->size[1] = 8;
  emxEnsureCapacity_real_T(new_tree, i3);
  for (i3 = 0; i3 < 8; i3++) {
    new_tree->data[i3] = 0.0;
  }

  flag1 = 0;
  emxInit_real_T(&tmp_dist, 2);
  emxInit_real_T(&sqr_dist, 1);
  emxInit_boolean_T(&x, 1);
  emxInit_int32_T(&ii, 1);
  while (flag1 == 0) {
    /*  select a random point */
    b_rand();
    randomPoint[0] = 20.0 * b_rand();
    randomPoint[1] = 20.0 * b_rand();

    /*  find node that is closest to randomPoint (Eucl. dist. between positions).  */
    /* tmp = tree(:, 1 : dim) - randomPoint(1 : dim); */
    i3 = tmp_dist->size[0] * tmp_dist->size[1];
    tmp_dist->size[0] = tree->size[0];
    tmp_dist->size[1] = 2;
    emxEnsureCapacity_real_T(tmp_dist, i3);
    n = tree->size[0] << 1;
    for (i3 = 0; i3 < n; i3++) {
      tmp_dist->data[i3] = 0.0;
    }

    i3 = tree->size[0];
    for (i = 0; i < i3; i++) {
      tmp_dist->data[i] = tree->data[i] - randomPoint[0];
      tmp_dist->data[i + tmp_dist->size[0]] = tree->data[i + tree->size[0]] -
        randomPoint[1];
    }

    /* fprintf("  sampled pos: %.1f, %.1f\n",randomPoint(1), randomPoint(2)); */
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
    d2 = b_norm(Vect);
    Vect[0] /= d2;
    Vect[1] /= d2;
    new_point[2] = 0.0;
    new_point[3] = 0.0;

    /*  find new_point that is within the range of min_parent_idx in terms of segmentLength (Eucl. dist. between positions). */
    if (cost_near > 16.0) {
      /*  generate a new point that is closest to randomPoint, segmentLength away from tree(idx,1:dim) */
      new_point[0] = tree->data[iidx] + Vect[0] * 4.0;
      new_point[1] = tree->data[iidx + tree->size[0]] + Vect[1] * 4.0;
    } else {
      new_point[0] = randomPoint[0];
      new_point[1] = randomPoint[1];
    }

    /* fprintf("adjusted pos: %.1f, %.1f\n",new_point(1), new_point(2)); */
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
      /* fprintf("no collision\n"); */
      for (i3 = 0; i3 < 8; i3++) {
        tree_data[i3] = tree->data[iidx + tree->size[0] * i3];
      }

      collisionFreeVel(tree_data, new_point, &near_cost, &unusedU4, randomPoint);

      /*  this collision checking includes a steering function and forward simulation */
      if (near_cost == 0.0) {
        new_point[2] = randomPoint[2];
        new_point[3] = randomPoint[3];

        /* fprintf("    adjusted pos: %.1f, %.1f\n",new_point(1), new_point(2)); */
        /* print tree */
        /*        fprintf("    tree =  "); */
        /*        for k = 1:size(tree,1) */
        /*          for j = 1:size(tree,2) */
        /*              fprintf(" %.1f, ",tree(k,j)) */
        /*          end */
        /*          fprintf("\n"); */
        /*        end */
        /* calculate the cost from root to to_point through parent from_node */
        min_cost = DI_costFreeVel(unusedU4, tree->data[iidx], tree->data[iidx +
          tree->size[0]], tree->data[iidx + (tree->size[0] << 1)], tree->
          data[iidx + tree->size[0] * 3], new_point[0], new_point[1]);

        /* fprintf("DI_costFreeVel %.1f\n", cost); */
        /*      fprintf("    from node "); */
        /*      for i = 1:6 */
        /*          fprintf(" %.2f, ",from_node(i)); */
        /*      end */
        /*      fprintf("\n"); */
        /* fprintf("    parent cost %.1f\n", from_node( 2 * dim + 2)); */
        /*  NOTE do not use broadcasting */
        min_cost += tree->data[iidx + tree->size[0] * 5];

        /*  total cost from root to new_point through its parent tree(idx,:) */
        /* fprintf("    base cost from %.0f = %.2f\n",idx, min_cost); */
        /* fprintf("tree size = %.0f\n",size(tree,1)); */
        /*  NOTE, broadcast does not translate to C */
        /* tmp_dist = tree(:, 1 : dim) - new_point(1 : dim); */
        i3 = tmp_dist->size[0] * tmp_dist->size[1];
        tmp_dist->size[0] = tree->size[0];
        tmp_dist->size[1] = 2;
        emxEnsureCapacity_real_T(tmp_dist, i3);
        n = tree->size[0] << 1;
        for (i3 = 0; i3 < n; i3++) {
          tmp_dist->data[i3] = 0.0;
        }

        i3 = tree->size[0];
        for (i = 0; i < i3; i++) {
          tmp_dist->data[i] = tree->data[i] - new_point[0];
          tmp_dist->data[i + tmp_dist->size[0]] = tree->data[i + tree->size[0]]
            - new_point[1];
        }

        /* fprintf("tmp_dist size = %.0f\n",size(tmp_dist,1)); */
        /*  find near neighbors    */
        cost_near = fmin(60.0 * sqrt(log((double)tree->size[0] + 1.0) / (double)
          tree->size[0]), 3.0);

        /*        fprintf("    dist to node: "); */
        /*        for index = 1:size(dist_sqr,1) */
        /*            fprintf(" %.1f, ",dist_sqr(index).^0.5); */
        /*        end */
        /*        fprintf("\n"); */
        cost_near *= cost_near;
        sqr_eucl_dist(tmp_dist, sqr_dist);
        i3 = x->size[0];
        x->size[0] = sqr_dist->size[0];
        emxEnsureCapacity_boolean_T(x, i3);
        n = sqr_dist->size[0];
        for (i3 = 0; i3 < n; i3++) {
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
        n = ii->size[0];
        for (i3 = 0; i3 < n; i3++) {
          sqr_dist->data[i3] = ii->data[i3];
        }

        /* fprintf("    near idx count = %.0f\n",size(near_idx,1)); */
        if (sqr_dist->size[0] >= 1) {
          i3 = sqr_dist->size[0];
          for (i = 0; i < i3; i++) {
            /*  choose parent node */
            /* fprintf("    checking index %.0f\n",near_idx(i)); */
            if ((int)sqr_dist->data[i] != iidx + 1) {
              n = (int)sqr_dist->data[i];
              for (i4 = 0; i4 < 8; i4++) {
                tree_data[i4] = tree->data[(n + tree->size[0] * i4) - 1];
              }

              segment_costFreeVel(tree_data, new_point, &cost_near, &unusedU4);
              cost_near += tree->data[((int)sqr_dist->data[i] + tree->size[0] *
                5) - 1];

              /* fprintf("    checking %.0f, cost = %.2f + %.2f(seg) = %.2f \n",near_idx(i), tree(near_idx(i), 2 * dim + 2), segcost,cost_near); */
              if (cost_near + 0.001 < min_cost) {
                n = (int)sqr_dist->data[i];
                for (i4 = 0; i4 < 8; i4++) {
                  tree_data[i4] = tree->data[(n + tree->size[0] * i4) - 1];
                }

                collisionFreeVelKnowTf(tree_data, new_point, &unusedU4,
                  &near_cost, randomPoint);
                if (near_cost == 0.0) {
                  min_cost = cost_near;
                  min_parent_idx = (int)sqr_dist->data[i] - 1;
                  new_point[2] = randomPoint[2];
                  new_point[3] = randomPoint[3];
                }
              }
            }
          }
        }

        new_node[0] = new_point[0];
        new_node[1] = new_point[1];
        new_node[2] = new_point[2];
        new_node[3] = new_point[3];
        new_node[4] = 0.0;
        new_node[5] = min_cost;
        new_node[6] = min_parent_idx + 1;
        new_node[7] = 0.0;

        /* fprintf("    adding new code %.0f, parent = %.0f, cost = %.2f\n",size(tree,1)+1,min_parent_idx,min_cost); */
        if (tree->size[0] != 0) {
          n = tree->size[0];
        } else {
          n = 0;
        }

        i3 = new_tree->size[0] * new_tree->size[1];
        new_tree->size[0] = n + 1;
        new_tree->size[1] = 8;
        emxEnsureCapacity_real_T(new_tree, i3);
        for (i3 = 0; i3 < 8; i3++) {
          for (i4 = 0; i4 < n; i4++) {
            new_tree->data[i4 + new_tree->size[0] * i3] = tree->data[i4 +
              tree->size[0] * i3];
          }
        }

        for (i3 = 0; i3 < 8; i3++) {
          new_tree->data[n + new_tree->size[0] * i3] = new_node[i3];
        }

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
            /* fprintf("    checking index %.0f\n",reduced_idx(j)); */
            if ((int)sqr_dist->data[flag1] != min_parent_idx + 1) {
              near_cost = new_tree->data[((int)sqr_dist->data[flag1] +
                new_tree->size[0] * 5) - 1];
              n = (int)sqr_dist->data[flag1];
              for (i4 = 0; i4 < 8; i4++) {
                tree_data[i4] = new_tree->data[(n + new_tree->size[0] * i4) - 1];
              }

              segment_cost(new_point, tree_data, &cost_near, &unusedU4);
              rnewcost = min_cost + cost_near;

              /* fprintf("checking %.0f, cost = %.2f + %.2f(seg) = %.2f \n",reduced_idx(j), near_cost, rnewcost); */
              if (near_cost > rnewcost + 0.001) {
                n = (int)sqr_dist->data[flag1];
                for (i4 = 0; i4 < 8; i4++) {
                  tree_data[i4] = new_tree->data[(n + new_tree->size[0] * i4) -
                    1];
                }

                near_cost = collisionKnowTf(new_node, tree_data, &unusedU4);
                if (near_cost == 0.0) {
                  cost_near = sqr_dist->data[flag1];
                  n = (int)new_tree->data[((int)sqr_dist->data[flag1] +
                    new_tree->size[0] * 6) - 1];
                  nx = 0;
                  for (i = 0; i < 4000; i++) {
                    b0 = (GChild[(n + 4000 * i) - 1] == cost_near);
                    bv0[i] = b0;
                    if (b0) {
                      nx++;
                    }
                  }

                  n = 0;
                  for (i = 0; i < 4000; i++) {
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
                                  * 5) - 1] = rnewcost;

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

  emxFree_int32_T(&ii);
  emxFree_boolean_T(&x);
  emxFree_real_T(&sqr_dist);
  emxFree_real_T(&tmp_dist);
  n = 0;

  /*  check to see if new node connects directly to end_node */
  Vect[0] = new_node[0] - 18.0;
  Vect[1] = new_node[1] - 18.0;
  if (b_norm(Vect) < 0.2) {
    n = 1;
    new_tree->data[(new_tree->size[0] + (new_tree->size[0] << 2)) - 1] = 1.0;

    /*  mark node as connecting to end. */
  } else {
    Vect[0] = new_node[0] - 18.0;
    Vect[1] = new_node[1] - 18.0;
    if (b_norm(Vect) < 4.0) {
      collision(new_node, &near_cost, &unusedU4);
      if (near_cost == 0.0) {
        n = 1;
        new_tree->data[(new_tree->size[0] + (new_tree->size[0] << 2)) - 1] = 1.0;

        /*  mark node as connecting to end. */
      }
    }
  }

  *flag = n;
}

/* End of code generation (RRTstar3D.c) */
