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
#include "rt_nonfinite.h"
#include "benchmarkRRT.h"
#include "benchmarkRRT_emxutil.h"
#include "RRTstar3D.h"

/* Custom Source Code */
#include "main.h"

/* Function Definitions */
double benchmarkRRT(void)
{
  double retval;
  emxArray_real_T *Its;
  emxArray_real_T *time;
  double cost_data[1];
  int cost_size[2];
  emxInit_real_T(&Its, 2);
  emxInit_real_T(&time, 2);

  /*  clc; */
  /*  close all; */
  /*  clear all; */
  /* tic */
  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
  RRTstar3D(Its, time, cost_data, cost_size);
  retval = 0.0;

  /* toc */
  emxFree_real_T(&time);
  emxFree_real_T(&Its);
  return retval;
}

/* End of code generation (benchmarkRRT.c) */
