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
#include <string.h>
#include "rt_nonfinite.h"
#include "benchmarkRRT.h"
#include "benchmarkRRT_emxutil.h"
#include "segment_cost.h"
#include "norm.h"
#include "RRTstar3D.h"

/* Function Definitions */
double benchmarkRRT(void)
{
  double retval;
  emxArray_real_T *tree;
  int i0;
  static const signed char start_node[8] = { 2, 2, 0, 0, 0, 0, 0, 0 };

  static double GChild[400000];
  emxArray_real_T *connectingNodes;
  int i;
  emxArray_boolean_T *idx;
  double flag;
  int loop_ub;
  int end;
  emxArray_int32_T *r0;
  double b_connectingNodes[2];
  double connectingNodes_data[8];
  creal_T lcost;
  creal_T unusedU0;
  int exitg1;
  emxInit_real_T(&tree, 2);

  /*  clc; */
  /*  close all; */
  /*  clear all; */
  /* tic */
  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
  /*  firstSol=0; */
  /*  Planning in state space */
  /*  Create random world */
  /*  Each Row Contains States, ConnectToEnd flag, Costz, ParentNodeIdx, and ChildNum */
  /*  Establish tree starting with the start node */
  i0 = tree->size[0] * tree->size[1];
  tree->size[0] = 1;
  tree->size[1] = 8;
  emxEnsureCapacity_real_T(tree, i0);
  for (i0 = 0; i0 < 8; i0++) {
    tree->data[i0] = start_node[i0];
  }

  /* coder.varsize('GChild') */
  memset(&GChild[0], 0, 400000U * sizeof(double));

  /* GChild = [0]; */
  /*  tic */
  /*  check to see if start_node connects directly to end_node */
  emxInit_real_T(&connectingNodes, 2);
  
  // NOTE timing
  struct timeval stop, start;
  gettimeofday(&start, NULL);
  for (i = 0; i < 2000; i++) {
    extendTree(tree, GChild, connectingNodes, &flag);
    i0 = tree->size[0] * tree->size[1];
    tree->size[0] = connectingNodes->size[0];
    tree->size[1] = 8;
    emxEnsureCapacity_real_T(tree, i0);
    loop_ub = connectingNodes->size[0] * connectingNodes->size[1];
    for (i0 = 0; i0 < loop_ub; i0++) {
      tree->data[i0] = connectingNodes->data[i0];
    }

    /*        if numPaths==1 && firstSol==0 */
    /*          toc   */
    /*          firstSol=1; */
    /*          tree_small = tree; */
    /*          its */
    /*        end */
    /*        if its==500 */
    /*          toc */
    /*          tree_500 = tree; */
    /*        end */
    /*        if its==1000 */
    /*          toc */
    /*          tree_1000 = tree; */
    /*        end */
    /*        if its==1500 */
    /*          toc */
    /*          tree_1500 = tree; */
    /*        end */
  }
  gettimeofday(&stop, NULL);
  printf("sampling 2000 nodes took %lu us\n", (stop.tv_sec - start.tv_sec) * 1000000 + stop.tv_usec - start.tv_usec); 
  emxInit_boolean_T(&idx, 1);

  /*  toc */
  /*  run_time = toc */
  /*  figure; */
  /*  path_first = findMinimumPath(tree_small,end_node,dim); */
  /*  plotExpandedTree(world,tree_small,dim); */
  /*  plotWorld(world,path_first,dim); */
  /*  plotTraj(path_first,dim); */
  /*  figure; */
  /*  path_500 = findMinimumPath(tree_500,end_node,dim); */
  /*  plotExpandedTree(world,tree_500,dim); */
  /*  plotWorld(world,path_500,dim); */
  /*  plotTraj(path_500,dim); */
  /*  figure; */
  /*  path_1000 = findMinimumPath(tree_1000,end_node,dim); */
  /*  plotExpandedTree(world,tree_1000,dim); */
  /*  plotWorld(world,path_1000,dim); */
  /*  plotTraj(path_1000,dim); */
  /*  figure; */
  /*  path_1500 = findMinimumPath(tree_1500,end_node,dim); */
  /*  plotExpandedTree(world,tree_1500,dim); */
  /*  plotWorld(world,path_1500,dim); */
  /*  plotTraj(path_1500,dim); */
  /*  figure; */
  /*  figure; */
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

  end = idx->size[0] - 1;
  loop_ub = 0;
  for (i = 0; i <= end; i++) {
    if (idx->data[i]) {
      loop_ub++;
    }
  }

  emxInit_int32_T(&r0, 1);
  i0 = r0->size[0];
  r0->size[0] = loop_ub;
  emxEnsureCapacity_int32_T(r0, i0);
  loop_ub = 0;
  for (i = 0; i <= end; i++) {
    if (idx->data[i]) {
      r0->data[loop_ub] = i + 1;
      loop_ub++;
    }
  }

  emxFree_boolean_T(&idx);
  i0 = connectingNodes->size[0] * connectingNodes->size[1];
  connectingNodes->size[0] = r0->size[0];
  connectingNodes->size[1] = 8;
  emxEnsureCapacity_real_T(connectingNodes, i0);
  for (i0 = 0; i0 < 8; i0++) {
    loop_ub = r0->size[0];
    for (i = 0; i < loop_ub; i++) {
      connectingNodes->data[i + connectingNodes->size[0] * i0] = tree->data
        [(r0->data[i] + tree->size[0] * i0) - 1];
    }
  }

  emxFree_real_T(&tree);

  /* NOTE */
  i0 = r0->size[0];
  emxFree_int32_T(&r0);
  for (loop_ub = 0; loop_ub < i0; loop_ub++) {
    b_connectingNodes[0] = connectingNodes->data[loop_ub] - 18.0;
    b_connectingNodes[1] = connectingNodes->data[loop_ub + connectingNodes->
      size[0]] - 18.0;
    if (b_norm(b_connectingNodes) >= 0.2) {
      for (i = 0; i < 8; i++) {
        connectingNodes_data[i] = connectingNodes->data[loop_ub +
          connectingNodes->size[0] * i];
      }

      c_segment_cost(connectingNodes_data, &lcost, &unusedU0);
      connectingNodes->data[loop_ub + connectingNodes->size[0] * 5] += lcost.re;
    }
  }

  /*  find minimum cost last node */
  i0 = connectingNodes->size[0];
  i = connectingNodes->size[0];
  if (i <= 2) {
    i0 = connectingNodes->size[0];
    if (i0 == 1) {
      end = 1;
    } else if ((connectingNodes->data[connectingNodes->size[0] * 5] >
                connectingNodes->data[1 + connectingNodes->size[0] * 5]) ||
               (rtIsNaN(connectingNodes->data[connectingNodes->size[0] * 5]) &&
                (!rtIsNaN(connectingNodes->data[1 + connectingNodes->size[0] * 5]))))
    {
      end = 2;
    } else {
      end = 1;
    }
  } else {
    if (!rtIsNaN(connectingNodes->data[connectingNodes->size[0] * 5])) {
      end = 1;
    } else {
      end = 0;
      loop_ub = 2;
      do {
        exitg1 = 0;
        i = connectingNodes->size[0];
        if (loop_ub <= i) {
          if (!rtIsNaN(connectingNodes->data[(loop_ub + connectingNodes->size[0]
                * 5) - 1])) {
            end = loop_ub;
            exitg1 = 1;
          } else {
            loop_ub++;
          }
        } else {
          exitg1 = 1;
        }
      } while (exitg1 == 0);
    }

    if (end == 0) {
      end = 1;
    } else {
      flag = connectingNodes->data[(end + connectingNodes->size[0] * 5) - 1];
      i = end + 1;
      for (loop_ub = i; loop_ub <= i0; loop_ub++) {
        if (flag > connectingNodes->data[(loop_ub + connectingNodes->size[0] * 5)
            - 1]) {
          flag = connectingNodes->data[(loop_ub + connectingNodes->size[0] * 5)
            - 1];
          end = loop_ub;
        }
      }
    }
  }

  /*      cost */
  /*  construct lowest cost path */
  b_connectingNodes[0] = connectingNodes->data[end - 1] - 18.0;
  b_connectingNodes[1] = connectingNodes->data[(end + connectingNodes->size[0])
    - 1] - 18.0;
  if (b_norm(b_connectingNodes) >= 0.2) {
    for (i0 = 0; i0 < 8; i0++) {
      connectingNodes_data[i0] = connectingNodes->data[(end +
        connectingNodes->size[0] * i0) - 1];
    }

    c_segment_cost(connectingNodes_data, &lcost, &unusedU0);
  }

  emxFree_real_T(&connectingNodes);

  /*  plotExpandedTree(world,tree,dim); */
  /*  plotWorld(world,path,dim); hold on */
  /*  plotTraj(path,dim); */
  retval = 0.0;

  /* toc */
  return retval;
}

/* End of code generation (benchmarkRRT.c) */
