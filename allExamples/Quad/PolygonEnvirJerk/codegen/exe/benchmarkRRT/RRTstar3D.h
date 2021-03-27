/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * RRTstar3D.h
 *
 * Code generation for function 'RRTstar3D'
 *
 */

#ifndef RRTSTAR3D_H
#define RRTSTAR3D_H

/* Include files */
#include <stddef.h>
#include <stdlib.h>
#include "rtwtypes.h"
#include "benchmarkRRT_types.h"

/* Function Declarations */
extern void extendTree(const emxArray_real_T *tree, double GChild[360000],
  emxArray_real_T *new_tree, double *flag);

#endif

/* End of code generation (RRTstar3D.h) */
