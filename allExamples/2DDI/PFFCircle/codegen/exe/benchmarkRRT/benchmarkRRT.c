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
#include "findMinCost.h"
#include "RRTstar3D.h"

/* Custom Source Code */
#include "main.h"
#include <stdio.h>

/* Function Definitions */
double benchmarkRRT(void)
{
  emxArray_real_T *tree;
  int imid;
  static const signed char start_node[8] = { 2, 2, 0, 0, 0, 0, 0, 0 };

  static double GChild[16000000];
  double numPaths;
  emxArray_real_T *connectingNodes;
  int i;
  emxArray_boolean_T *idx;
  double flag;
  int n;
  boolean_T tf;
  int ilo;
  int ihi;
  boolean_T exitg1;
  static const short iv0[5] = { 400, 1000, 2000, 3000, 4000 };

  emxArray_int32_T *r0;
  double b_connectingNodes[2];
  double connectingNodes_data[8];
  creal_T unusedU0;
  int exitg2;
  emxInit_real_T(&tree, 2);

  /*  clc; */
  /*  close all; */
  /*  clear all; */
  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
  /*  firstSol=0; */
  /*  Planning in state space */
  /*  Create random world */
  /*  Each Row Contains States, ConnectToEnd flag, Costz, ParentNodeIdx, and ChildNum */
  /*  Establish tree starting with the start node */
  imid = tree->size[0] * tree->size[1];
  tree->size[0] = 1;
  tree->size[1] = 8;
  emxEnsureCapacity_real_T(tree, imid);
  for (imid = 0; imid < 8; imid++) {
    tree->data[imid] = start_node[imid];
  }

  /* coder.varsize('GChild') */
  memset(&GChild[0], 0, 16000000U * sizeof(double));

  /* GChild = [0]; */
  /* tic */
  /*  check to see if start_node connects directly to end_node */
  numPaths = 0.0;
  emxInit_real_T(&connectingNodes, 2);
  for (i = 0; i < 4000; i++) {
    extendTree(tree, GChild, connectingNodes, &flag);
    imid = tree->size[0] * tree->size[1];
    tree->size[0] = connectingNodes->size[0];
    tree->size[1] = 8;
    emxEnsureCapacity_real_T(tree, imid);
    n = connectingNodes->size[0] * connectingNodes->size[1];
    for (imid = 0; imid < n; imid++) {
      tree->data[imid] = connectingNodes->data[imid];
    }

    numPaths += flag;

    /*        if mod(i,100) == 0 */
    /*            fprintf("%.0f nodes \n",i); */
    /*        end */
    /*  report first solution */
    /*  its, time, cost */
    if ((flag == 1.0) && (numPaths == 0.0)) {
      flag = findMinCost(connectingNodes);
      printf("nodes: %.0f, min cost %.3f \n", 1.0 + (double)i, flag);
      fflush(stdout);
    }

    tf = false;
    n = -1;
    ilo = 1;
    ihi = 5;
    exitg1 = false;
    while ((!exitg1) && (ihi >= ilo)) {
      imid = ((ilo >> 1) + (ihi >> 1)) - 1;
      if (((ilo & 1) == 1) && ((ihi & 1) == 1)) {
        imid++;
      }

      if (1 + i == iv0[imid]) {
        n = imid;
        exitg1 = true;
      } else if (1 + i < iv0[imid]) {
        ihi = imid;
      } else {
        ilo = imid + 2;
      }
    }

    if (n + 1 > 0) {
      while ((n > 0) && (1 + i == iv0[n - 1])) {
        n--;
      }
    }

    if (n + 1 > 0) {
      tf = true;
    }

    if (tf) {
      flag = findMinCost(connectingNodes);
      printf("nodes: %.0f, min cost %.3f \n", 1.0 + (double)i, flag);
      fflush(stdout);
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
  /* NOTE */
  /*  find nodes that connect to end_node */
  /*      connectingNodes = []; */
  /*      for i=1:size(tree,1) */
  /*          if tree(i,2*dim+1)==1 */
  /*              tree(i,2*dim+2)=tree(i,2*dim+2)+segment_cost(tree(i,:),end_node(1:2*dim),dim); */
  /*              connectingNodes = [connectingNodes ; tree(i,:)]; */
  /*          end */
  /*      end */
  n = tree->size[0];
  imid = idx->size[0];
  idx->size[0] = n;
  emxEnsureCapacity_boolean_T(idx, imid);
  for (imid = 0; imid < n; imid++) {
    idx->data[imid] = (tree->data[imid + (tree->size[0] << 2)] == 1.0);
  }

  ilo = idx->size[0] - 1;
  n = 0;
  for (i = 0; i <= ilo; i++) {
    if (idx->data[i]) {
      n++;
    }
  }

  emxInit_int32_T(&r0, 1);
  imid = r0->size[0];
  r0->size[0] = n;
  emxEnsureCapacity_int32_T(r0, imid);
  n = 0;
  for (i = 0; i <= ilo; i++) {
    if (idx->data[i]) {
      r0->data[n] = i + 1;
      n++;
    }
  }

  emxFree_boolean_T(&idx);
  imid = connectingNodes->size[0] * connectingNodes->size[1];
  connectingNodes->size[0] = r0->size[0];
  connectingNodes->size[1] = 8;
  emxEnsureCapacity_real_T(connectingNodes, imid);
  for (imid = 0; imid < 8; imid++) {
    n = r0->size[0];
    for (ihi = 0; ihi < n; ihi++) {
      connectingNodes->data[ihi + connectingNodes->size[0] * imid] = tree->data
        [(r0->data[ihi] + tree->size[0] * imid) - 1];
    }
  }

  emxFree_real_T(&tree);

  /* NOTE */
  imid = r0->size[0];
  emxFree_int32_T(&r0);
  for (n = 0; n < imid; n++) {
    b_connectingNodes[0] = connectingNodes->data[n] - 18.0;
    b_connectingNodes[1] = connectingNodes->data[n + connectingNodes->size[0]] -
      18.0;
    if (b_norm(b_connectingNodes) >= 0.2) {
      for (ihi = 0; ihi < 8; ihi++) {
        connectingNodes_data[ihi] = connectingNodes->data[n +
          connectingNodes->size[0] * ihi];
      }

      b_segment_cost(connectingNodes_data, &flag, &unusedU0);
      connectingNodes->data[n + connectingNodes->size[0] * 5] += flag;
    }
  }

  /*  find minimum cost last node */
  imid = connectingNodes->size[0];
  ihi = connectingNodes->size[0];
  if (ihi <= 2) {
    imid = connectingNodes->size[0];
    if (imid == 1) {
      ilo = 1;
    } else if ((connectingNodes->data[connectingNodes->size[0] * 5] >
                connectingNodes->data[1 + connectingNodes->size[0] * 5]) ||
               (rtIsNaN(connectingNodes->data[connectingNodes->size[0] * 5]) &&
                (!rtIsNaN(connectingNodes->data[1 + connectingNodes->size[0] * 5]))))
    {
      ilo = 2;
    } else {
      ilo = 1;
    }
  } else {
    if (!rtIsNaN(connectingNodes->data[connectingNodes->size[0] * 5])) {
      ilo = 1;
    } else {
      ilo = 0;
      n = 2;
      do {
        exitg2 = 0;
        ihi = connectingNodes->size[0];
        if (n <= ihi) {
          if (!rtIsNaN(connectingNodes->data[(n + connectingNodes->size[0] * 5)
                       - 1])) {
            ilo = n;
            exitg2 = 1;
          } else {
            n++;
          }
        } else {
          exitg2 = 1;
        }
      } while (exitg2 == 0);
    }

    if (ilo == 0) {
      ilo = 1;
    } else {
      flag = connectingNodes->data[(ilo + connectingNodes->size[0] * 5) - 1];
      ihi = ilo + 1;
      for (n = ihi; n <= imid; n++) {
        if (flag > connectingNodes->data[(n + connectingNodes->size[0] * 5) - 1])
        {
          flag = connectingNodes->data[(n + connectingNodes->size[0] * 5) - 1];
          ilo = n;
        }
      }
    }
  }

  /*  construct lowest cost path */
  b_connectingNodes[0] = connectingNodes->data[ilo - 1] - 18.0;
  b_connectingNodes[1] = connectingNodes->data[(ilo + connectingNodes->size[0])
    - 1] - 18.0;
  if (b_norm(b_connectingNodes) >= 0.2) {
    for (imid = 0; imid < 8; imid++) {
      connectingNodes_data[imid] = connectingNodes->data[(ilo +
        connectingNodes->size[0] * imid) - 1];
    }

    b_segment_cost(connectingNodes_data, &flag, &unusedU0);
  }

  emxFree_real_T(&connectingNodes);

  /*  plotExpandedTree(world,tree,dim); */
  /*  plotWorld(world,path,dim); hold on */
  /*  plotTraj(path,dim); */
  return 0.0;
}

/* End of code generation (benchmarkRRT.c) */
