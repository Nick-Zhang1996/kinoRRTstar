/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: _coder_getTime_api.h
 *
 * MATLAB Coder version            : 4.1
 * C/C++ source code generated on  : 22-Feb-2021 23:20:10
 */

#ifndef _CODER_GETTIME_API_H
#define _CODER_GETTIME_API_H

/* Include Files */
#include "tmwtypes.h"
#include "mex.h"
#include "emlrt.h"
#include <stddef.h>
#include <stdlib.h>
#include "_coder_getTime_api.h"

/* Variable Declarations */
extern emlrtCTX emlrtRootTLSGlobal;
extern emlrtContext emlrtContextGlobal;

/* Function Declarations */
extern real_T getTime(void);
extern void getTime_api(int32_T nlhs, const mxArray *plhs[1]);
extern void getTime_atexit(void);
extern void getTime_initialize(void);
extern void getTime_terminate(void);
extern void getTime_xil_terminate(void);

#endif

/*
 * File trailer for _coder_getTime_api.h
 *
 * [EOF]
 */
