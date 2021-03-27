/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * _coder_getTime_api.c
 *
 * Code generation for function '_coder_getTime_api'
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "getTime.h"
#include "_coder_getTime_api.h"
#include "getTime_data.h"

/* Function Declarations */
static const mxArray *emlrt_marshallOut(const real_T u);

/* Function Definitions */
static const mxArray *emlrt_marshallOut(const real_T u)
{
  const mxArray *y;
  const mxArray *m0;
  y = NULL;
  m0 = emlrtCreateDoubleScalar(u);
  emlrtAssign(&y, m0);
  return y;
}

void getTime_api(int32_T nlhs, const mxArray *plhs[1])
{
  real_T current_time_in_sec;
  emlrtStack st = { NULL,              /* site */
    NULL,                              /* tls */
    NULL                               /* prev */
  };

  (void)nlhs;
  st.tls = emlrtRootTLSGlobal;

  /* Invoke the target function */
  current_time_in_sec = getTime(&st);

  /* Marshall function outputs */
  plhs[0] = emlrt_marshallOut(current_time_in_sec);
}

/* End of code generation (_coder_getTime_api.c) */
