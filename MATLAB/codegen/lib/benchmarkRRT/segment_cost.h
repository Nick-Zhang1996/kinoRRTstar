/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * segment_cost.h
 *
 * Code generation for function 'segment_cost'
 *
 */

#ifndef SEGMENT_COST_H
#define SEGMENT_COST_H

/* Include files */
#include <stddef.h>
#include <stdlib.h>
#include "rtwtypes.h"
#include "benchmarkRRT_types.h"

/* Function Declarations */
extern void b_segment_cost(const double from_node[4], const double
  to_point_data[], creal_T *cost, creal_T *Tf);
extern void c_segment_cost(const double from_node_data[], creal_T *cost, creal_T
  *Tf);
extern void segment_cost(const double from_node_data[], const double to_point[4],
  creal_T *cost, creal_T *Tf);

#endif

/* End of code generation (segment_cost.h) */
