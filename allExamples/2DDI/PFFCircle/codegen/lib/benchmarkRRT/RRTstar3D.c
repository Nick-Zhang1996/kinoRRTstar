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
#include "collisionFreeVel.h"
#include "rand.h"
#include "benchmarkRRT_data.h"

/* Function Declarations */
static void extendTree(const emxArray_real_T *tree, emxArray_real_T *new_tree);

/* Function Definitions */
static void extendTree(const emxArray_real_T *tree, emxArray_real_T *new_tree)
{
  int collision_flag;
  emxArray_real_T *b_tree;
  double randomPoint[4];
  int loop_ub;
  double b_collision_flag;
  double sqr_e_dist_idx_0;
  double tmp[2];
  double d2;
  double new_point[4];
  boolean_T exitg1;
  double tree_data[8];
  creal_T unusedU4;
  static const signed char iv2[5] = { 6, 15, 10, 5, 15 };

  static const signed char iv3[5] = { 6, 14, 11, 15, 5 };

  /*  segmentLength: maximum stepsize, r: neighbor radius */
  /*  NOTE */
  collision_flag = new_tree->size[0] * new_tree->size[1];
  new_tree->size[0] = 1;
  new_tree->size[1] = 8;
  emxEnsureCapacity_real_T(new_tree, collision_flag);
  for (collision_flag = 0; collision_flag < 8; collision_flag++) {
    new_tree->data[collision_flag] = 0.0;
  }

  emxInit_real_T(&b_tree, 2);
  while (1) {
    /*  select a random point */
    b_rand();
    randomPoint[0] = 20.0 * b_rand();
    randomPoint[1] = 20.0 * b_rand();

    /*  find node that is closest to randomPoint (Eucl. dist. between positions).  */
    loop_ub = tree->size[0];
    collision_flag = b_tree->size[0] * b_tree->size[1];
    b_tree->size[0] = loop_ub;
    b_tree->size[1] = 2;
    emxEnsureCapacity_real_T(b_tree, collision_flag);
    for (collision_flag = 0; collision_flag < loop_ub; collision_flag++) {
      b_tree->data[collision_flag] = tree->data[collision_flag];
    }

    for (collision_flag = 0; collision_flag < loop_ub; collision_flag++) {
      b_tree->data[collision_flag + b_tree->size[0]] = tree->data[collision_flag
        + tree->size[0]];
    }

    b_collision_flag = b_tree->data[0] - randomPoint[0];
    sqr_e_dist_idx_0 = b_collision_flag * b_collision_flag;
    b_collision_flag = b_tree->data[1] - randomPoint[1];
    tmp[0] = randomPoint[0] - tree->data[0];
    tmp[1] = randomPoint[1] - tree->data[tree->size[0]];
    d2 = b_norm(tmp);
    tmp[0] /= d2;
    tmp[1] /= d2;
    new_point[2] = 0.0;
    new_point[3] = 0.0;

    /*  find new_point that is within the range of min_parent_idx in terms of segmentLength (Eucl. dist. between positions). */
    if (sqr_e_dist_idx_0 + b_collision_flag * b_collision_flag > 16.0) {
      /*  generate a new point that is closest to randomPoint, segmentLength away from tree(idx,1:dim) */
      new_point[0] = tree->data[0] + tmp[0] * 4.0;
      new_point[1] = tree->data[tree->size[0]] + tmp[1] * 4.0;
    } else {
      new_point[0] = randomPoint[0];
      new_point[1] = randomPoint[1];
    }

    /*  check if the new_point is in collision */
    /*  check if a point is in collision */
    collision_flag = 0;
    if ((new_point[0] > 20.0) || (new_point[0] < 0.0)) {
      collision_flag = 1;
    }

    if ((new_point[1] > 20.0) || (new_point[1] < 0.0)) {
      collision_flag = 1;
    }

    if (collision_flag == 0) {
      /*  check each obstacle */
      loop_ub = 0;
      exitg1 = false;
      while ((!exitg1) && (loop_ub < 5)) {
        b_collision_flag = new_point[0] - (double)iv2[loop_ub];
        b_collision_flag *= b_collision_flag;
        tmp[0] = b_collision_flag;
        b_collision_flag = new_point[1] - (double)iv3[loop_ub];
        b_collision_flag *= b_collision_flag;
        if (sqrt(tmp[0] + b_collision_flag) < dv0[loop_ub] + 0.1) {
          /*  (norm([p(1);p(2)]-[world.cx(i); world.cy(i)])<world.radius(i)+0.1) */
          collision_flag = 1;
          exitg1 = true;
        } else {
          loop_ub++;
        }
      }
    }

    if (collision_flag == 0) {
      for (collision_flag = 0; collision_flag < 8; collision_flag++) {
        tree_data[collision_flag] = tree->data[tree->size[0] * collision_flag];
      }

      collisionFreeVel(tree_data, new_point, &b_collision_flag, &unusedU4,
                       randomPoint);

      /*  this collision checking includes a steering function and forward simulation */
      if (b_collision_flag == 0.0) {
        /* calculate the cost from root to to_point through parent from_node */
        /*  total cost from root to new_point through its parent tree(idx,:) */
        /*  find near neighbors    */
        /*  ChildNum + 1 */
        abort();
      }
    }
  }
}

void RRTstar3D(emxArray_real_T *Its, emxArray_real_T *time, double Cost_data[],
               int Cost_size[2])
{
  emxArray_real_T *tree;
  int firstSol;
  int i0;
  static const signed char start_node[8] = { 2, 2, 0, 0, 0, 0, 0, 0 };

  int i;
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
  double ex;
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

  /*  GChild = zeros(400, 500); */
  /* tic */
  /*  check to see if start_node connects directly to end_node */
  i = 1;
  emxInit_real_T(&connectingNodes, 2);
  for (b_i = 0; b_i < 4000; b_i++) {
    i = b_i + 1;
    extendTree(tree, connectingNodes);

    /* [tree,GChild,flag] = [a,b,c]; */
    i0 = tree->size[0] * tree->size[1];
    tree->size[0] = connectingNodes->size[0];
    tree->size[1] = 8;
    emxEnsureCapacity_real_T(tree, i0);
    loop_ub = connectingNodes->size[0] * connectingNodes->size[1];
    for (i0 = 0; i0 < loop_ub; i0++) {
      tree->data[i0] = connectingNodes->data[i0];
    }

    if (((int)c == 1.0) && (firstSol == 0)) {
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

      segment_cost(connectingNodes_data, &c, &unusedU0);
      connectingNodes->data[firstSol + connectingNodes->size[0] * 5] += c;
    }
  }

  /*  find minimum cost last node */
  i0 = connectingNodes->size[0];
  i1 = connectingNodes->size[0];
  if (i1 <= 2) {
    i0 = connectingNodes->size[0];
    if (i0 == 1) {
      ex = connectingNodes->data[connectingNodes->size[0] * 5];
      b_i = 1;
    } else if ((connectingNodes->data[connectingNodes->size[0] * 5] >
                connectingNodes->data[1 + connectingNodes->size[0] * 5]) ||
               (rtIsNaN(connectingNodes->data[connectingNodes->size[0] * 5]) &&
                (!rtIsNaN(connectingNodes->data[1 + connectingNodes->size[0] * 5]))))
    {
      ex = connectingNodes->data[1 + connectingNodes->size[0] * 5];
      b_i = 2;
    } else {
      ex = connectingNodes->data[connectingNodes->size[0] * 5];
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
      ex = connectingNodes->data[connectingNodes->size[0] * 5];
      b_i = 1;
    } else {
      ex = connectingNodes->data[(b_i + connectingNodes->size[0] * 5) - 1];
      i1 = b_i + 1;
      for (firstSol = i1; firstSol <= i0; firstSol++) {
        if (ex > connectingNodes->data[(firstSol + connectingNodes->size[0] * 5)
            - 1]) {
          ex = connectingNodes->data[(firstSol + connectingNodes->size[0] * 5) -
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

    segment_cost(connectingNodes_data, &c, &unusedU0);
  }

  emxFree_real_T(&connectingNodes);
  Cost_size[0] = 1;
  Cost_size[1] = 1;
  Cost_data[0] = ex;

  /*  plotExpandedTree(tree,dim); */
  /*  plotWorld(world,path,dim); hold on */
  /*  plotTraj(path,dim); */
}

/* End of code generation (RRTstar3D.c) */
