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
#include "rt_nonfinite.h"
#include "benchmarkRRT.h"
#include "findMinCost.h"
#include "benchmarkRRT_emxutil.h"
#include "segment_cost.h"
#include "norm.h"

/* Function Definitions */
double findMinCost(const emxArray_real_T *tree)
{
  double min_cost;
  emxArray_boolean_T *close_to_goal_idx;
  int loop_ub;
  int i2;
  int end;
  int i;
  emxArray_int32_T *r1;
  emxArray_real_T *connectingNodes;
  double b_connectingNodes[2];
  double connectingNodes_data[8];
  creal_T lcost;
  creal_T unusedU0;
  double total_cost;
  emxInit_boolean_T(&close_to_goal_idx, 1);

  /*  find nodes that connect to end_node */
  /*      connectingNodes = []; */
  /*      for i=1:size(tree,1) */
  /*          if tree(i,2*dim+1)==1 */
  /*              tree(i,2*dim+2)=tree(i,2*dim+2)+segment_cost(tree(i,:),end_node(1:2*dim),dim); */
  /*              connectingNodes = [connectingNodes ; tree(i,:)]; */
  /*          end */
  /*      end */
  loop_ub = tree->size[0];
  i2 = close_to_goal_idx->size[0];
  close_to_goal_idx->size[0] = loop_ub;
  emxEnsureCapacity_boolean_T(close_to_goal_idx, i2);
  for (i2 = 0; i2 < loop_ub; i2++) {
    close_to_goal_idx->data[i2] = (tree->data[i2 + (tree->size[0] << 2)] == 1.0);
  }

  /* fprintf("close to goal index %.1f \n",sum(nonzeros(close_to_goal_idx))); */
  end = close_to_goal_idx->size[0] - 1;
  loop_ub = 0;
  for (i = 0; i <= end; i++) {
    if (close_to_goal_idx->data[i]) {
      loop_ub++;
    }
  }

  emxInit_int32_T(&r1, 1);
  i2 = r1->size[0];
  r1->size[0] = loop_ub;
  emxEnsureCapacity_int32_T(r1, i2);
  loop_ub = 0;
  for (i = 0; i <= end; i++) {
    if (close_to_goal_idx->data[i]) {
      r1->data[loop_ub] = i + 1;
      loop_ub++;
    }
  }

  emxFree_boolean_T(&close_to_goal_idx);
  emxInit_real_T(&connectingNodes, 2);
  i2 = connectingNodes->size[0] * connectingNodes->size[1];
  connectingNodes->size[0] = r1->size[0];
  connectingNodes->size[1] = 8;
  emxEnsureCapacity_real_T(connectingNodes, i2);
  for (i2 = 0; i2 < 8; i2++) {
    loop_ub = r1->size[0];
    for (i = 0; i < loop_ub; i++) {
      connectingNodes->data[i + connectingNodes->size[0] * i2] = tree->data
        [(r1->data[i] + tree->size[0] * i2) - 1];
    }
  }

  /* DEBUG */
  /* fprintf("connecting Nodes  : %.1f, %.1f, %.1f, %.1f, %.1f, %.1f, %.1f, %.1f, \n",connectingNodes(i,1),connectingNodes(i,2),connectingNodes(i,3),connectingNodes(i,4),connectingNodes(i,5),connectingNodes(i,6),connectingNodes(i,7),connectingNodes(i,8)); */
  min_cost = 10000.0;

  /* fprintf("connecting Nodes size %.0f",size(connectingNodes,1)); */
  i2 = r1->size[0];
  for (loop_ub = 0; loop_ub < i2; loop_ub++) {
    b_connectingNodes[0] = connectingNodes->data[loop_ub] - 18.0;
    b_connectingNodes[1] = connectingNodes->data[loop_ub + connectingNodes->
      size[0]] - 18.0;
    if (b_norm(b_connectingNodes) >= 0.2) {
      for (i = 0; i < 8; i++) {
        connectingNodes_data[i] = connectingNodes->data[loop_ub +
          connectingNodes->size[0] * i];
      }

      c_segment_cost(connectingNodes_data, &lcost, &unusedU0);

      /* fprintf("evaluate node cost, local to end = %.2f \n",lcost); */
      total_cost = tree->data[(r1->data[loop_ub] + tree->size[0] * 5) - 1] +
        lcost.re;

      /* fprintf("evaluate node cost, total = %.2f \n",total_cost); */
      if (total_cost < min_cost) {
        min_cost = total_cost;
      }
    }
  }

  emxFree_int32_T(&r1);
  emxFree_real_T(&connectingNodes);
  return min_cost;
}

/* End of code generation (findMinCost.c) */
