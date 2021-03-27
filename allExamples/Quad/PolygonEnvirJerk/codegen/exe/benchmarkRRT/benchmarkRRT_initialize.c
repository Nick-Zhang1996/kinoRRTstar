/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * benchmarkRRT_initialize.c
 *
 * Code generation for function 'benchmarkRRT_initialize'
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "benchmarkRRT.h"
#include "benchmarkRRT_initialize.h"
#include "eml_rand_mt19937ar_stateful.h"

/* Function Definitions */
void benchmarkRRT_initialize(void)
{
  rt_InitInfAndNaN(8U);
  c_eml_rand_mt19937ar_stateful_i();
}

/* End of code generation (benchmarkRRT_initialize.c) */
