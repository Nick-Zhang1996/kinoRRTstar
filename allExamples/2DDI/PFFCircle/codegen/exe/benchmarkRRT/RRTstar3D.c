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
#include "rt_nonfinite.h"
#include "benchmarkRRT.h"
#include "RRTstar3D.h"
#include "benchmarkRRT_emxutil.h"
#include "segment_cost.h"
#include "norm.h"
#include "DI_costFreeVel.h"
#include "collisionFreeVel.h"
#include "rand.h"
#include "collision.h"
#include "benchmarkRRT_data.h"
#include <stdio.h>

/* Function Declarations */
static void extendTree(const emxArray_real_T *tree, emxArray_real_T *GChild,
  emxArray_real_T *new_tree, double *flag);

/* Function Definitions */
static void extendTree(const emxArray_real_T *tree, emxArray_real_T *GChild,
  emxArray_real_T *new_tree, double *flag)
{
  int i3;
  double new_node_data[8];
  int flag1;
  emxArray_real_T *result;
  cell_wrap_2 reshapes[2];
  emxArray_real_T *b_tree;
  double randomPoint[4];
  int input_sizes_idx_0;
  double tmp[2];
  int loop_ub;
  double collision_flag;
  creal_T unusedU4;
  double sqr_e_dist_idx_0;
  double d2;
  double new_point[4];
  int i;
  boolean_T exitg1;
  double tree_data[8];
  static const signed char iv4[5] = { 6, 15, 10, 5, 15 };

  static const signed char iv5[5] = { 6, 14, 11, 15, 5 };

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
  emxInit_real_T(&result, 2);
  emxInitMatrix_cell_wrap_2(reshapes);
  emxInit_real_T(&b_tree, 2);
  while (flag1 == 0) {
    /*  select a random point */
    b_rand();
    randomPoint[0] = 20.0 * b_rand();
    randomPoint[1] = 20.0 * b_rand();

    /*  find node that is closest to randomPoint (Eucl. dist. between positions).  */
    loop_ub = tree->size[0];
    i3 = b_tree->size[0] * b_tree->size[1];
    b_tree->size[0] = loop_ub;
    b_tree->size[1] = 2;
    emxEnsureCapacity_real_T(b_tree, i3);
    for (i3 = 0; i3 < loop_ub; i3++) {
      b_tree->data[i3] = tree->data[i3];
    }

    for (i3 = 0; i3 < loop_ub; i3++) {
      b_tree->data[i3 + b_tree->size[0]] = tree->data[i3 + tree->size[0]];
    }

    collision_flag = b_tree->data[0] - randomPoint[0];
    sqr_e_dist_idx_0 = collision_flag * collision_flag;
    collision_flag = b_tree->data[1] - randomPoint[1];
    tmp[0] = randomPoint[0] - tree->data[0];
    tmp[1] = randomPoint[1] - tree->data[tree->size[0]];
    d2 = b_norm(tmp);
    tmp[0] /= d2;
    tmp[1] /= d2;
    new_point[2] = 0.0;
    new_point[3] = 0.0;

    /*  find new_point that is within the range of min_parent_idx in terms of segmentLength (Eucl. dist. between positions). */
    if (sqr_e_dist_idx_0 + collision_flag * collision_flag > 16.0) {
      /*  generate a new point that is closest to randomPoint, segmentLength away from tree(idx,1:dim) */
      new_point[0] = tree->data[0] + tmp[0] * 4.0;
      new_point[1] = tree->data[tree->size[0]] + tmp[1] * 4.0;
    } else {
      new_point[0] = randomPoint[0];
      new_point[1] = randomPoint[1];
    }

    /*  check if the new_point is in collision */
    /*  check if a point is in collision */
    input_sizes_idx_0 = 0;
    if ((new_point[0] > 20.0) || (new_point[0] < 0.0)) {
      input_sizes_idx_0 = 1;
    }

    if ((new_point[1] > 20.0) || (new_point[1] < 0.0)) {
      input_sizes_idx_0 = 1;
    }

    if (input_sizes_idx_0 == 0) {
      /*  check each obstacle */
      i = 0;
      exitg1 = false;
      while ((!exitg1) && (i < 5)) {
        collision_flag = new_point[0] - (double)iv4[i];
        collision_flag *= collision_flag;
        tmp[0] = collision_flag;
        collision_flag = new_point[1] - (double)iv5[i];
        collision_flag *= collision_flag;
        if (sqrt(tmp[0] + collision_flag) < dv0[i] + 0.1) {
          /*  (norm([p(1);p(2)]-[world.cx(i); world.cy(i)])<world.radius(i)+0.1) */
          input_sizes_idx_0 = 1;
          exitg1 = true;
        } else {
          i++;
        }
      }
    }

    if (input_sizes_idx_0 == 0) {
      for (i3 = 0; i3 < 8; i3++) {
        tree_data[i3] = tree->data[tree->size[0] * i3];
      }

      collisionFreeVel(tree_data, new_point, &collision_flag, &unusedU4,
                       randomPoint);

      /*  this collision checking includes a steering function and forward simulation */
      if (collision_flag == 0.0) {
        /* calculate the cost from root to to_point through parent from_node */
        collision_flag = DI_costFreeVel(unusedU4, tree->data[0], tree->data
          [tree->size[0]], tree->data[tree->size[0] << 1], tree->data[tree->
          size[0] * 3], new_point[0], new_point[1]);
        collision_flag += tree->data[tree->size[0] * 5];

        /*  total cost from root to new_point through its parent tree(idx,:) */
        /*  find near neighbors    */
        new_node_data[0] = new_point[0];
        new_node_data[1] = new_point[1];
        new_node_data[2] = randomPoint[2];
        new_node_data[3] = randomPoint[3];
        new_node_data[4] = 0.0;
        new_node_data[5] = collision_flag;
        new_node_data[6] = 1.0;
        new_node_data[7] = 0.0;
        if (tree->size[0] != 0) {
          input_sizes_idx_0 = tree->size[0];
        } else {
          input_sizes_idx_0 = 0;
        }

        i3 = reshapes[0].f1->size[0] * reshapes[0].f1->size[1];
        reshapes[0].f1->size[0] = input_sizes_idx_0;
        reshapes[0].f1->size[1] = 8;
        emxEnsureCapacity_real_T(reshapes[0].f1, i3);
        loop_ub = input_sizes_idx_0 << 3;
        for (i3 = 0; i3 < loop_ub; i3++) {
          reshapes[0].f1->data[i3] = tree->data[i3];
        }

        i3 = reshapes[1].f1->size[0] * reshapes[1].f1->size[1];
        reshapes[1].f1->size[0] = 1;
        reshapes[1].f1->size[1] = 8;
        emxEnsureCapacity_real_T(reshapes[1].f1, i3);
        for (i3 = 0; i3 < 8; i3++) {
          reshapes[1].f1->data[i3] = new_node_data[i3];
        }

        i3 = result->size[0] * result->size[1];
        result->size[0] = reshapes[0].f1->size[0] + reshapes[1].f1->size[0];
        result->size[1] = reshapes[0].f1->size[1];
        emxEnsureCapacity_real_T(result, i3);
        loop_ub = reshapes[0].f1->size[1];
        for (i3 = 0; i3 < loop_ub; i3++) {
          i = reshapes[0].f1->size[0];
          for (flag1 = 0; flag1 < i; flag1++) {
            result->data[flag1 + result->size[0] * i3] = reshapes[0].f1->
              data[flag1 + reshapes[0].f1->size[0] * i3];
          }
        }

        loop_ub = reshapes[1].f1->size[1];
        for (i3 = 0; i3 < loop_ub; i3++) {
          i = reshapes[1].f1->size[0];
          for (flag1 = 0; flag1 < i; flag1++) {
            result->data[(flag1 + reshapes[0].f1->size[0]) + result->size[0] *
              i3] = reshapes[1].f1->data[flag1 + reshapes[1].f1->size[0] * i3];
          }
        }

        i3 = new_tree->size[0] * new_tree->size[1];
        new_tree->size[0] = input_sizes_idx_0 + 1;
        new_tree->size[1] = 8;
        emxEnsureCapacity_real_T(new_tree, i3);
        for (i3 = 0; i3 < 8; i3++) {
          for (flag1 = 0; flag1 < input_sizes_idx_0; flag1++) {
            new_tree->data[flag1 + new_tree->size[0] * i3] = tree->data[flag1 +
              input_sizes_idx_0 * i3];
          }
        }

        for (i3 = 0; i3 < 8; i3++) {
          for (flag1 = 0; flag1 < 1; flag1++) {
            new_tree->data[input_sizes_idx_0 + new_tree->size[0] * i3] =
              new_node_data[i3];
          }
        }

        new_tree->data[new_tree->size[0] * 7] = result->data[result->size[0] * 7]
          + 1.0;

        /*  ChildNum + 1 */
        GChild->data[4000 * ((int)new_tree->data[new_tree->size[0] * 7] - 1)] =
          result->size[0];

        /*  update GChild matrix */
        flag1 = 1;
      }
    }
  }

  emxFree_real_T(&b_tree);
  emxFreeMatrix_cell_wrap_2(reshapes);
  emxFree_real_T(&result);
  input_sizes_idx_0 = 0;

  /*  check to see if new node connects directly to end_node */
  tmp[0] = new_node_data[0] - 18.0;
  tmp[1] = new_node_data[1] - 18.0;
  if (b_norm(tmp) < 0.2) {
    input_sizes_idx_0 = 1;
    new_tree->data[(new_tree->size[0] + (new_tree->size[0] << 2)) - 1] = 1.0;

    /*  mark node as connecting to end. */
  } else {
    tmp[0] = new_node_data[0] - 18.0;
    tmp[1] = new_node_data[1] - 18.0;
    if (b_norm(tmp) < 4.0) {
      collision(new_node_data, &collision_flag, &unusedU4);
      if (collision_flag == 0.0) {
        input_sizes_idx_0 = 1;
        new_tree->data[(new_tree->size[0] + (new_tree->size[0] << 2)) - 1] = 1.0;

        /*  mark node as connecting to end. */
      }
    }
  }

  *flag = input_sizes_idx_0;
}

void RRTstar3D(emxArray_real_T *Its, emxArray_real_T *time, double Cost_data[],
               int Cost_size[2])
{
  emxArray_real_T *tree;
  int firstSol;
  int i0;
  emxArray_real_T *GChild;
  static const signed char start_node[8] = { 2, 2, 0, 0, 0, 0, 0, 0 };

  double numPaths;
  int i;
  emxArray_real_T *b;
  emxArray_real_T *connectingNodes;
  int b_i;
  emxArray_boolean_T *idx;
  int i1;
  int loop_ub;
  double c;
  emxArray_int32_T *r0;
  double b_connectingNodes[2];
  double connectingNodes_data[8];
  creal_T unusedU0;
  int exitg1;
  emxInit_real_T(&tree, 2);
  time->size[0] = 1;
  time->size[1] = 0;
  Its->size[0] = 1;
  Its->size[1] = 0;
  firstSol = 0;

  /*  planning in state space */
  /*  create random world */
  /*  Each Row Contains States, ConnectToEnd flag, Costz, ParentNodeIdx, and ChildNum */
  /*  establish tree starting with the start node */
  i0 = tree->size[0] * tree->size[1];
  tree->size[0] = 1;
  tree->size[1] = 8;
  emxEnsureCapacity_real_T(tree, i0);
  for (i0 = 0; i0 < 8; i0++) {
    tree->data[i0] = start_node[i0];
  }

  emxInit_real_T(&GChild, 2);
  i0 = GChild->size[0] * GChild->size[1];
  GChild->size[0] = 4000;
  GChild->size[1] = 4000;
  emxEnsureCapacity_real_T(GChild, i0);
  for (i0 = 0; i0 < 16000000; i0++) {
    GChild->data[i0] = 0.0;
  }

  /* tic */
  /*  check to see if start_node connects directly to end_node */
  numPaths = 0.0;
  i = 1;
  emxInit_real_T(&b, 2);
  emxInit_real_T(&connectingNodes, 2);
  for (b_i = 0; b_i < 4000; b_i++) {
    i = b_i + 1;
    i0 = b->size[0] * b->size[1];
    b->size[0] = GChild->size[0];
    b->size[1] = GChild->size[1];
    emxEnsureCapacity_real_T(b, i0);
    loop_ub = GChild->size[0] * GChild->size[1];
    for (i0 = 0; i0 < loop_ub; i0++) {
      b->data[i0] = GChild->data[i0];
    }

    extendTree(tree, b, connectingNodes, &c);

    /* [tree,GChild,flag] = [a,b,c]; */
    i0 = tree->size[0] * tree->size[1];
    tree->size[0] = connectingNodes->size[0];
    tree->size[1] = 8;
    emxEnsureCapacity_real_T(tree, i0);
    loop_ub = connectingNodes->size[0] * connectingNodes->size[1];
    for (i0 = 0; i0 < loop_ub; i0++) {
      tree->data[i0] = connectingNodes->data[i0];
    }

    i0 = GChild->size[0] * GChild->size[1];
    GChild->size[0] = 4000;
    GChild->size[1] = 4000;
    emxEnsureCapacity_real_T(GChild, i0);
    for (i0 = 0; i0 < 16000000; i0++) {
      GChild->data[i0] = b->data[i0];
    }

    if (fmod(1.0 + (double)b_i, 100.0) == 0.0) {
      printf("nodes = %.0f\n", 1.0 + (double)b_i);
      fflush(stdout);
    }

    numPaths += c;
    if ((numPaths == 1.0) && (firstSol == 0)) {
      firstSol = 1;
      i0 = Its->size[1];
      i1 = Its->size[0] * Its->size[1];
      Its->size[1] = i0 + 1;
      emxEnsureCapacity_real_T(Its, i1);
      Its->data[i0] = 1.0 + (double)b_i;
      i0 = time->size[1];
      i1 = time->size[0] * time->size[1];
      time->size[1] = i0 + 1;
      emxEnsureCapacity_real_T(time, i1);
      time->data[i0] = 0.0;
    }

    if (1 + b_i == 200) {
      i0 = Its->size[1];
      i1 = Its->size[0] * Its->size[1];
      Its->size[1] = i0 + 1;
      emxEnsureCapacity_real_T(Its, i1);
      Its->data[i0] = 200.0;
      i0 = time->size[1];
      i1 = time->size[0] * time->size[1];
      time->size[1] = i0 + 1;
      emxEnsureCapacity_real_T(time, i1);
      time->data[i0] = 0.0;
    }

    if (1 + b_i == 400) {
      i0 = Its->size[1];
      i1 = Its->size[0] * Its->size[1];
      Its->size[1] = i0 + 1;
      emxEnsureCapacity_real_T(Its, i1);
      Its->data[i0] = 400.0;
      i0 = time->size[1];
      i1 = time->size[0] * time->size[1];
      time->size[1] = i0 + 1;
      emxEnsureCapacity_real_T(time, i1);
      time->data[i0] = 0.0;
    }

    if (1 + b_i == 600) {
      i0 = Its->size[1];
      i1 = Its->size[0] * Its->size[1];
      Its->size[1] = i0 + 1;
      emxEnsureCapacity_real_T(Its, i1);
      Its->data[i0] = 600.0;
      i0 = time->size[1];
      i1 = time->size[0] * time->size[1];
      time->size[1] = i0 + 1;
      emxEnsureCapacity_real_T(time, i1);
      time->data[i0] = 0.0;
    }

    if (1 + b_i == 800) {
      i0 = Its->size[1];
      i1 = Its->size[0] * Its->size[1];
      Its->size[1] = i0 + 1;
      emxEnsureCapacity_real_T(Its, i1);
      Its->data[i0] = 800.0;
      i0 = time->size[1];
      i1 = time->size[0] * time->size[1];
      time->size[1] = i0 + 1;
      emxEnsureCapacity_real_T(time, i1);
      time->data[i0] = 0.0;
    }

    if (1 + b_i == 1000) {
      i0 = Its->size[1];
      i1 = Its->size[0] * Its->size[1];
      Its->size[1] = i0 + 1;
      emxEnsureCapacity_real_T(Its, i1);
      Its->data[i0] = 1000.0;
      i0 = time->size[1];
      i1 = time->size[0] * time->size[1];
      time->size[1] = i0 + 1;
      emxEnsureCapacity_real_T(time, i1);
      time->data[i0] = 0.0;
    }

    if (1 + b_i == 1200) {
      i0 = Its->size[1];
      i1 = Its->size[0] * Its->size[1];
      Its->size[1] = i0 + 1;
      emxEnsureCapacity_real_T(Its, i1);
      Its->data[i0] = 1200.0;
      i0 = time->size[1];
      i1 = time->size[0] * time->size[1];
      time->size[1] = i0 + 1;
      emxEnsureCapacity_real_T(time, i1);
      time->data[i0] = 0.0;
    }

    if (1 + b_i == 1400) {
      i0 = Its->size[1];
      i1 = Its->size[0] * Its->size[1];
      Its->size[1] = i0 + 1;
      emxEnsureCapacity_real_T(Its, i1);
      Its->data[i0] = 1400.0;
      i0 = time->size[1];
      i1 = time->size[0] * time->size[1];
      time->size[1] = i0 + 1;
      emxEnsureCapacity_real_T(time, i1);
      time->data[i0] = 0.0;
    }

    if (1 + b_i == 1600) {
      i0 = Its->size[1];
      i1 = Its->size[0] * Its->size[1];
      Its->size[1] = i0 + 1;
      emxEnsureCapacity_real_T(Its, i1);
      Its->data[i0] = 1600.0;
      i0 = time->size[1];
      i1 = time->size[0] * time->size[1];
      time->size[1] = i0 + 1;
      emxEnsureCapacity_real_T(time, i1);
      time->data[i0] = 0.0;
    }

    if (1 + b_i == 1800) {
      i0 = Its->size[1];
      i1 = Its->size[0] * Its->size[1];
      Its->size[1] = i0 + 1;
      emxEnsureCapacity_real_T(Its, i1);
      Its->data[i0] = 1800.0;
      i0 = time->size[1];
      i1 = time->size[0] * time->size[1];
      time->size[1] = i0 + 1;
      emxEnsureCapacity_real_T(time, i1);
      time->data[i0] = 0.0;
    }

    if (1 + b_i == 2000) {
      i0 = Its->size[1];
      i1 = Its->size[0] * Its->size[1];
      Its->size[1] = i0 + 1;
      emxEnsureCapacity_real_T(Its, i1);
      Its->data[i0] = 2000.0;
      i0 = time->size[1];
      i1 = time->size[0] * time->size[1];
      time->size[1] = i0 + 1;
      emxEnsureCapacity_real_T(time, i1);
      time->data[i0] = 0.0;
    }

    if (1 + b_i == 2500) {
      i0 = Its->size[1];
      i1 = Its->size[0] * Its->size[1];
      Its->size[1] = i0 + 1;
      emxEnsureCapacity_real_T(Its, i1);
      Its->data[i0] = 2500.0;
      i0 = time->size[1];
      i1 = time->size[0] * time->size[1];
      time->size[1] = i0 + 1;
      emxEnsureCapacity_real_T(time, i1);
      time->data[i0] = 0.0;
    }

    if (1 + b_i == 3000) {
      i0 = Its->size[1];
      i1 = Its->size[0] * Its->size[1];
      Its->size[1] = i0 + 1;
      emxEnsureCapacity_real_T(Its, i1);
      Its->data[i0] = 3000.0;
      i0 = time->size[1];
      i1 = time->size[0] * time->size[1];
      time->size[1] = i0 + 1;
      emxEnsureCapacity_real_T(time, i1);
      time->data[i0] = 0.0;
    }

    if (1 + b_i == 3500) {
      i0 = Its->size[1];
      i1 = Its->size[0] * Its->size[1];
      Its->size[1] = i0 + 1;
      emxEnsureCapacity_real_T(Its, i1);
      Its->data[i0] = 3500.0;
      i0 = time->size[1];
      i1 = time->size[0] * time->size[1];
      time->size[1] = i0 + 1;
      emxEnsureCapacity_real_T(time, i1);
      time->data[i0] = 0.0;
    }
  }

  emxFree_real_T(&b);
  emxFree_real_T(&GChild);
  emxInit_boolean_T(&idx, 1);
  i0 = time->size[1];
  i1 = time->size[0] * time->size[1];
  time->size[1] = i0 + 1;
  emxEnsureCapacity_real_T(time, i1);
  time->data[i0] = 0.0;
  i0 = Its->size[1];
  i1 = Its->size[0] * Its->size[1];
  Its->size[1] = i0 + 1;
  emxEnsureCapacity_real_T(Its, i1);
  Its->data[i0] = i;

  /*  figure; */
  /*  [path_first,cost] = findMinimumPath(tree_small,end_node,dim); */
  /*  Cost=[Cost, cost]; */
  /*  plotExpandedTree(tree_small,dim); */
  /*  plotWorld(world,path_first,dim); */
  /*  plotTraj(path_first,dim); */
  /*  % figure; */
  /*  [path_200,cost] = findMinimumPath(tree_200,end_node,dim); */
  /*  Cost=[Cost, cost]; */
  /*  % plotExpandedTree(tree_200,dim); */
  /*  % plotWorld(world,path_200,dim); */
  /*  % plotTraj(path_200,dim); */
  /*  figure; */
  /*  [path_400,cost] = findMinimumPath(tree_400,end_node,dim); */
  /*  Cost=[Cost, cost]; */
  /*  plotExpandedTree(tree_400,dim); */
  /*  plotWorld(world,path_400,dim); */
  /*  plotTraj(path_400,dim); */
  /*  % figure; */
  /*  [path_600,cost] = findMinimumPath(tree_600,end_node,dim); */
  /*  Cost=[Cost, cost]; */
  /*  % plotExpandedTree(tree_600,dim); */
  /*  % plotWorld(world,path_600,dim); */
  /*  % plotTraj(path_600,dim); */
  /*  % figure; */
  /*  [path_800,cost] = findMinimumPath(tree_800,end_node,dim); */
  /*  Cost=[Cost, cost]; */
  /*  % plotExpandedTree(tree_800,dim); */
  /*  % plotWorld(world,path_800,dim); */
  /*  % plotTraj(path_800,dim); */
  /*  % figure; */
  /*  [path_1000,cost] = findMinimumPath(tree_1000,end_node,dim); */
  /*  Cost=[Cost, cost]; */
  /*  % plotExpandedTree(tree_1000,dim); */
  /*  % plotWorld(world,path_1000,dim); */
  /*  % plotTraj(path_1000,dim); */
  /*  % figure; */
  /*  [path_1200,cost] = findMinimumPath(tree_1200,end_node,dim); */
  /*  Cost=[Cost, cost]; */
  /*  % plotExpandedTree(tree_1200,dim); */
  /*  % plotWorld(world,path_1200,dim); */
  /*  % plotTraj(path_1200,dim); */
  /*  % figure; */
  /*  [path_1400,cost] = findMinimumPath(tree_1400,end_node,dim); */
  /*  Cost=[Cost, cost]; */
  /*  % plotExpandedTree(tree_1400,dim); */
  /*  % plotWorld(world,path_1400,dim); */
  /*  % plotTraj(path_1400,dim); */
  /*  % figure; */
  /*  [path_1600,cost] = findMinimumPath(tree_1600,end_node,dim); */
  /*  Cost=[Cost, cost]; */
  /*  % plotExpandedTree(tree_1600,dim); */
  /*  % plotWorld(world,path_1600,dim); */
  /*  % plotTraj(path_1600,dim); */
  /*  % figure; */
  /*  [path_1800,cost] = findMinimumPath(tree_1800,end_node,dim); */
  /*  Cost=[Cost, cost]; */
  /*  % plotExpandedTree(tree_1800,dim); */
  /*  % plotWorld(world,path_1800,dim); */
  /*  % plotTraj(path_1800,dim); */
  /*  figure; */
  /*  [path_2000,cost] = findMinimumPath(tree_2000,end_node,dim); */
  /*  Cost=[Cost, cost]; */
  /*  plotExpandedTree(tree_2000,dim); */
  /*  plotWorld(world,path_2000,dim); */
  /*  plotTraj(path_2000,dim); */
  /*  % % figure; */
  /*  [path_2500,cost] = findMinimumPath(tree_2500,end_node,dim); */
  /*  Cost=[Cost, cost]; */
  /*  % % plotExpandedTree(tree_2500,dim); */
  /*  % % plotWorld(world,path_2500,dim); */
  /*  % % plotTraj(path_2500,dim); */
  /*  % % figure; */
  /*  [path_3000,cost] = findMinimumPath(tree_3000,end_node,dim); */
  /*  Cost=[Cost, cost]; */
  /*  % % plotExpandedTree(tree_3000,dim); */
  /*  % % plotWorld(world,path_3000,dim); */
  /*  % % plotTraj(path_3000,dim); */
  /*  % % figure; */
  /*  [path_3500,cost] = findMinimumPath(tree_3500,end_node,dim); */
  /*  Cost=[Cost, cost]; */
  /*  % % plotExpandedTree(tree_3500,dim); */
  /*  % % plotWorld(world,path_3500,dim); */
  /*  % % plotTraj(path_3500,dim); */
  /*  figure; */
  /* NOTE */
  /*  find nodes that connect to end_node */
  /*      connectingNodes = []; */
  /*      for i=1:size(tree,1) */
  /*          if tree(i,2*dim+1)==1 */
  /*              tree(i,2*dim+2)=tree(i,2*dim+2)+segment_cost(tree(i,:),end_node(1:2*dim),dim); */
  /*              connectingNodes = [connectingNodes ; tree(i,:)]; */
  /*          end */
  /*      end */
  loop_ub = tree->size[0];
  i0 = idx->size[0];
  idx->size[0] = loop_ub;
  emxEnsureCapacity_boolean_T(idx, i0);
  for (i0 = 0; i0 < loop_ub; i0++) {
    idx->data[i0] = (tree->data[i0 + (tree->size[0] << 2)] == 1.0);
  }

  b_i = idx->size[0] - 1;
  firstSol = 0;
  for (i = 0; i <= b_i; i++) {
    if (idx->data[i]) {
      firstSol++;
    }
  }

  emxInit_int32_T(&r0, 1);
  i0 = r0->size[0];
  r0->size[0] = firstSol;
  emxEnsureCapacity_int32_T(r0, i0);
  firstSol = 0;
  for (i = 0; i <= b_i; i++) {
    if (idx->data[i]) {
      r0->data[firstSol] = i + 1;
      firstSol++;
    }
  }

  emxFree_boolean_T(&idx);
  i0 = connectingNodes->size[0] * connectingNodes->size[1];
  connectingNodes->size[0] = r0->size[0];
  connectingNodes->size[1] = 8;
  emxEnsureCapacity_real_T(connectingNodes, i0);
  for (i0 = 0; i0 < 8; i0++) {
    loop_ub = r0->size[0];
    for (i1 = 0; i1 < loop_ub; i1++) {
      connectingNodes->data[i1 + connectingNodes->size[0] * i0] = tree->data
        [(r0->data[i1] + tree->size[0] * i0) - 1];
    }
  }

  emxFree_real_T(&tree);

  /* NOTE */
  i0 = r0->size[0];
  emxFree_int32_T(&r0);
  for (firstSol = 0; firstSol < i0; firstSol++) {
    b_connectingNodes[0] = connectingNodes->data[firstSol] - 18.0;
    b_connectingNodes[1] = connectingNodes->data[firstSol +
      connectingNodes->size[0]] - 18.0;
    if (b_norm(b_connectingNodes) >= 0.2) {
      for (i1 = 0; i1 < 8; i1++) {
        connectingNodes_data[i1] = connectingNodes->data[firstSol +
          connectingNodes->size[0] * i1];
      }

      segment_cost(connectingNodes_data, &numPaths, &unusedU0);
      connectingNodes->data[firstSol + connectingNodes->size[0] * 5] += numPaths;
    }
  }

  /*  find minimum cost last node */
  i0 = connectingNodes->size[0];
  i1 = connectingNodes->size[0];
  if (i1 <= 2) {
    i0 = connectingNodes->size[0];
    if (i0 == 1) {
      c = connectingNodes->data[connectingNodes->size[0] * 5];
      b_i = 1;
    } else if ((connectingNodes->data[connectingNodes->size[0] * 5] >
                connectingNodes->data[1 + connectingNodes->size[0] * 5]) ||
               (rtIsNaN(connectingNodes->data[connectingNodes->size[0] * 5]) &&
                (!rtIsNaN(connectingNodes->data[1 + connectingNodes->size[0] * 5]))))
    {
      c = connectingNodes->data[1 + connectingNodes->size[0] * 5];
      b_i = 2;
    } else {
      c = connectingNodes->data[connectingNodes->size[0] * 5];
      b_i = 1;
    }
  } else {
    if (!rtIsNaN(connectingNodes->data[connectingNodes->size[0] * 5])) {
      b_i = 1;
    } else {
      b_i = 0;
      firstSol = 2;
      do {
        exitg1 = 0;
        i1 = connectingNodes->size[0];
        if (firstSol <= i1) {
          if (!rtIsNaN(connectingNodes->data[(firstSol + connectingNodes->size[0]
                * 5) - 1])) {
            b_i = firstSol;
            exitg1 = 1;
          } else {
            firstSol++;
          }
        } else {
          exitg1 = 1;
        }
      } while (exitg1 == 0);
    }

    if (b_i == 0) {
      c = connectingNodes->data[connectingNodes->size[0] * 5];
      b_i = 1;
    } else {
      c = connectingNodes->data[(b_i + connectingNodes->size[0] * 5) - 1];
      i1 = b_i + 1;
      for (firstSol = i1; firstSol <= i0; firstSol++) {
        if (c > connectingNodes->data[(firstSol + connectingNodes->size[0] * 5)
            - 1]) {
          c = connectingNodes->data[(firstSol + connectingNodes->size[0] * 5) -
            1];
          b_i = firstSol;
        }
      }
    }
  }

  /*  construct lowest cost path */
  b_connectingNodes[0] = connectingNodes->data[b_i - 1] - 18.0;
  b_connectingNodes[1] = connectingNodes->data[(b_i + connectingNodes->size[0])
    - 1] - 18.0;
  if (b_norm(b_connectingNodes) >= 0.2) {
    for (i0 = 0; i0 < 8; i0++) {
      connectingNodes_data[i0] = connectingNodes->data[(b_i +
        connectingNodes->size[0] * i0) - 1];
    }

    segment_cost(connectingNodes_data, &numPaths, &unusedU0);
  }

  emxFree_real_T(&connectingNodes);
  Cost_size[0] = 1;
  Cost_size[1] = 1;
  Cost_data[0] = c;

  /*  plotExpandedTree(tree,dim); */
  /*  plotWorld(world,path,dim); hold on */
  /*  plotTraj(path,dim); */
}

/* End of code generation (RRTstar3D.c) */
