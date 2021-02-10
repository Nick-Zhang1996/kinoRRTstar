/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * _coder_benchmarkRRT_api.h
 *
 * Code generation for function '_coder_benchmarkRRT_api'
 *
 */

#ifndef _CODER_BENCHMARKRRT_API_H
#define _CODER_BENCHMARKRRT_API_H

/* Include files */
#include "tmwtypes.h"
#include "mex.h"
#include "emlrt.h"
#include <stddef.h>
#include <stdlib.h>
#include "_coder_benchmarkRRT_api.h"

/* Variable Declarations */
extern emlrtCTX emlrtRootTLSGlobal;
extern emlrtContext emlrtContextGlobal;

/* Function Declarations */
extern real_T benchmarkRRT(void);
extern void benchmarkRRT_api(int32_T nlhs, const mxArray *plhs[1]);
extern void benchmarkRRT_atexit(void);
extern void benchmarkRRT_initialize(void);
extern void benchmarkRRT_terminate(void);
extern void benchmarkRRT_xil_terminate(void);

#endif

/* End of code generation (_coder_benchmarkRRT_api.h) */
