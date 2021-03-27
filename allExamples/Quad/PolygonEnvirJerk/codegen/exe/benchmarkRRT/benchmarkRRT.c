/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * benchmarkRRT.c
 *
 * Code generation for function 'benchmarkRRT'
 *
 */

/* Include files */
#include <math.h>
#include <string.h>
#include "rt_nonfinite.h"
#include "benchmarkRRT.h"
#include "benchmarkRRT_emxutil.h"
#include "findMinCost.h"
#include "RRTstar3D.h"

/* Custom Source Code */
#include "main.h"
#include <stdio.h>

/* Function Definitions */
double benchmarkRRT(void)
{
  emxArray_real_T *tree;
  int k;
  static const signed char start_node[13] = { 2, 2, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0 };

  static double GChild[360000];
  double numPaths;
  emxArray_real_T *b_tree;
  int i;
  double flag;
  int loop_ub;
  boolean_T tf;
  boolean_T exitg1;
  int exponent;
  emxInit_real_T(&tree, 2);

  /*  clc; */
  /*  close all; */
  /*  clear all; */
  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
  printf("starting run\n");
  fflush(stdout);

  /*  planning in state space */
  /*  create random world */
  printf("world created... \n");
  fflush(stdout);

  /*  Each Row Contains States, ConnectToEnd flag, Costz, ParentNodeIdx, and ChildNum */
  /*  establish tree starting with the start node */
  k = tree->size[0] * tree->size[1];
  tree->size[0] = 1;
  tree->size[1] = 13;
  emxEnsureCapacity_real_T(tree, k);
  for (k = 0; k < 13; k++) {
    tree->data[k] = start_node[k];
  }

  /* GChild  = []; */
  /* coder.varsize('GChild') */
  memset(&GChild[0], 0, 360000U * sizeof(double));

  /*  tic */
  /*  check to see if start_node connects directly to end_node */
  numPaths = 0.0;
  emxInit_real_T(&b_tree, 2);
  for (i = 0; i < 600; i++) {
    /* fprintf("trying node... \n"); */
    extendTree(tree, GChild, b_tree, &flag);
    k = tree->size[0] * tree->size[1];
    tree->size[0] = b_tree->size[0];
    tree->size[1] = 13;
    emxEnsureCapacity_real_T(tree, k);
    loop_ub = b_tree->size[0] * b_tree->size[1];
    for (k = 0; k < loop_ub; k++) {
      tree->data[k] = b_tree->data[k];
    }

    /* fprintf("nodes added"); */
    numPaths += flag;
    if ((flag == 1.0) && (numPaths == 1.0)) {
      flag = findMinCost(b_tree);
      printf("nodes: %.0f, min cost %.3f \n", 1.0 + (double)i, flag);
      fflush(stdout);
    }

    tf = false;
    k = 0;
    exitg1 = false;
    while ((!exitg1) && (k < 4)) {
      frexp((200.0 * (double)k + 400.0) / 2.0, &exponent);
      if (fabs((double)(200 * k - i) + 399.0) < ldexp(1.0, exponent - 53)) {
        tf = true;
        exitg1 = true;
      } else {
        k++;
      }
    }

    if (tf) {
      flag = findMinCost(b_tree);
      printf("nodes: %.0f, min cost %.3f \n", 1.0 + (double)i, flag);
      fflush(stdout);
    }
  }

  emxFree_real_T(&b_tree);
  emxFree_real_T(&tree);

  /*  if show_output == 1 */
  /*  figure; */
  /*  [path_first,cost] = findMinimumPath(tree_small,end_node,dim,state_dim); */
  /*  Cost=[Cost, cost]; */
  /*  plotExpandedTree(tree_small,dim,state_dim); */
  /*  plotWorld(world,dim); */
  /*  plotTraj(path_first,dim); */
  /*  figure; */
  /*  [path_400,cost] = findMinimumPath(tree_400,end_node,dim,state_dim); */
  /*  Cost=[Cost, cost]; */
  /*  plotExpandedTree(tree_400,dim,state_dim); */
  /*  plotWorld(world,dim); */
  /*  plotTraj(path_400,dim); */
  /*  % figure; */
  /*  [path_600,cost] = findMinimumPath(tree_600,end_node,dim,state_dim); */
  /*  Cost=[Cost, cost]; */
  /*  % plotExpandedTree(tree_600,dim,state_dim); */
  /*  % plotWorld(world,dim); */
  /*  % plotTraj(path_600,dim); */
  /*  % figure; */
  /*  [path_800,cost] = findMinimumPath(tree_800,end_node,dim,state_dim); */
  /*  Cost=[Cost, cost]; */
  /*  % plotExpandedTree(tree_800,dim,state_dim); */
  /*  % plotWorld(world,dim); */
  /*  % plotTraj(path_800,dim); */
  /*  % figure; */
  /*  [path_1000,cost] = findMinimumPath(tree_1000,end_node,dim,state_dim); */
  /*  Cost=[Cost, cost]; */
  /*  % % plotExpandedTree(tree_1000,dim,state_dim); */
  /*  % plotWorld(world,dim); */
  /*  % plotTraj(path_1000,dim); */
  /*  % figure; */
  /*  [path_1200,cost] = findMinimumPath(tree_1200,end_node,dim,state_dim); */
  /*  Cost=[Cost, cost]; */
  /*  % plotExpandedTree(tree_1200,dim,state_dim); */
  /*  % plotWorld(world,dim); */
  /*  % plotTraj(path_1200,dim); */
  /*  % figure; */
  /*  [path_1400,cost] = findMinimumPath(tree_1400,end_node,dim,state_dim); */
  /*  Cost=[Cost, cost]; */
  /*  % plotExpandedTree(tree_1400,dim,state_dim); */
  /*  % plotWorld(world,dim); */
  /*  % plotTraj(path_1400,dim); */
  /*  % figure; */
  /*  [path_1600,cost] = findMinimumPath(tree_1600,end_node,dim,state_dim); */
  /*  Cost=[Cost, cost]; */
  /*  % plotExpandedTree(tree_1600,dim,state_dim); */
  /*  % plotWorld(world,dim); */
  /*  % plotTraj(path_1600,dim); */
  /*  % figure; */
  /*  [path_1800,cost] = findMinimumPath(tree_1800,end_node,dim,state_dim); */
  /*  Cost=[Cost, cost]; */
  /*  % plotExpandedTree(tree_1800,dim,state_dim); */
  /*  % plotWorld(world,dim); */
  /*  % plotTraj(path_1800,dim); */
  /*  figure; */
  /*  [path,~] = findMinimumPath(tree,end_node,dim,state_dim); */
  /*  % Cost=[Cost, cost]; */
  /*  % plotExpandedTree(tree,dim,state_dim); */
  /*  plotWorld(world,dim); hold on */
  /*  plotTraj(path,dim); */
  /*  end */
  return 0.0;
}

/* End of code generation (benchmarkRRT.c) */
