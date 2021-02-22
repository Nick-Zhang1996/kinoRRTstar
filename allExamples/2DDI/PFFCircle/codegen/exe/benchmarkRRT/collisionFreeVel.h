/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * collisionFreeVel.h
 *
 * Code generation for function 'collisionFreeVel'
 *
 */

#ifndef COLLISIONFREEVEL_H
#define COLLISIONFREEVEL_H

/* Include files */
#include <stddef.h>
#include <stdlib.h>
#include "rtwtypes.h"
#include "benchmarkRRT_types.h"

/* Function Declarations */
extern void collisionFreeVel(const double parent_data[], const double node[4],
  double *collision_flag, creal_T *Tf, double xf[4]);

#endif

/* End of code generation (collisionFreeVel.h) */
