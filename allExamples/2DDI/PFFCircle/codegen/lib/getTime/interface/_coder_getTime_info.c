/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: _coder_getTime_info.c
 *
 * MATLAB Coder version            : 4.1
 * C/C++ source code generated on  : 22-Feb-2021 23:20:10
 */

/* Include Files */
#include "_coder_getTime_info.h"

/* Function Definitions */

/*
 * Arguments    : void
 * Return Type  : const mxArray *
 */
const mxArray *emlrtMexFcnResolvedFunctionsInfo(void)
{
  const mxArray *nameCaptureInfo;
  const char * data[8] = {
    "789ced584d4f1341189e9a8278c0f4a2674dbc2849075a01c168fab9a5dad602052a04e37699d2a1b3bb65774b0bf160c21ff00778e8d1a3ffc078f30fe8c11f"
    "61bc7b71b7db699795814acb12b63b49337dfb74e779e799779ebc59e04b677d0080dbc01cf71e9af364270e74e61be0e4b0e3bece7cd316d33106fceddf266d",
    "f8c7ce2cc892869a9a19102ca15c5d2c21450f245e44dd657664114bbca4150e6b08284895c901da6923654c50018b28235b8225ac07226781ba810119dfe315"
    "245457eb22502a6a2f5d620d80459f4f8cfdfbfbd427c5d0c7f89f062a5d7c2bb99d58840959a88b48d25498c2da52bd04ab589257560aaac62b902724d9e4c5",
    "1a412a0c25126998e7b838560482e02ed28c3d06c56ede6f2d79f92c798d9f92b715f731663a6e81094bf4fb0be56b32d63b4da7d3f8ee30f8a84e1417e41da4"
    "04b15e3c8ac493a0aa2958da05e0dc7df79b07ab9e033a62e45197aa92dc90ba7cad01f9224c3e73df14df4a678af1ede422cc460b99680cae84a6679e94a026",
    "cba42437211249fb33d5d6074e5181a029505074ae2e228585773fa2ced6e1f9e7ef16be41efd95d061fad378a87e74af865484c37a6abe1cd6c741e1f9667d5"
    "782f8ffc393ce7e50118b153eb5b7dfd223a8e8aafb3f8fad569dc16d311d01143276ae0c3e2f3dbe21e9fbfcd275478a3c9e8f2d506e4b3f7373d3eb30e28be",
    "558cc7b6a91983d1f12bb7f339e5c76b9bb9f5b9d46a622fb7b4bfb6b77b98e7c37921e51e3f1e8d7bf838e2967ed93e5879d0e1f5cb2787d72fbbcb9fddde2f"
    "b718cf8feebd35488e23c3ea5b276c311d011d31f4c1aa5ae315155d559fdc1a902fc2e4fbf7fc93679f3fc12528f21ae14bd0d0a45c972055c751df9e78f6c6",
    "f3edebeedb8d8db4ca91ea3a172a2babe460efd50be968d9f36d77fbf6fb63af0fbf209f7bebc2ebc32f93cfebc387b3fe67c6f3fdea9863ac4f75a4781f7dd8"
    "99f7f6be242b224ff0118ab73b4947ebeedb87a1f5e5630cbe808e987db92237c0d5bdbf6e0dc8f794c967d603c5ffaf2f4744ff02dbdab4dfb23979fe913fbf",
    "7e7ef77cfc92f89cf2f1f9cc7c28ba50e4361b5a6e3abebfaf352ba16517bdef6e319e1fdd7b6b8c070ebe4fd177aae0e675f5ede74c3ef3fc297eb1f337b531"
    "4ac041dffefa28e5f5dfd7ddb773427523ad704299d38a6161667936f63a1b73816fff053a9820b0",
    "" };

  nameCaptureInfo = NULL;
  emlrtNameCaptureMxArrayR2016a(data, 10552U, &nameCaptureInfo);
  return nameCaptureInfo;
}

/*
 * Arguments    : void
 * Return Type  : mxArray *
 */
mxArray *emlrtMexFcnProperties(void)
{
  mxArray *xResult;
  mxArray *xEntryPoints;
  const char * fldNames[6] = { "Name", "NumberOfInputs", "NumberOfOutputs",
    "ConstantInputs", "FullPath", "TimeStamp" };

  mxArray *xInputs;
  const char * b_fldNames[4] = { "Version", "ResolvedFunctions", "EntryPoints",
    "CoverageInfo" };

  xEntryPoints = emlrtCreateStructMatrix(1, 1, 6, fldNames);
  xInputs = emlrtCreateLogicalMatrix(1, 0);
  emlrtSetField(xEntryPoints, 0, "Name", emlrtMxCreateString("getTime"));
  emlrtSetField(xEntryPoints, 0, "NumberOfInputs", emlrtMxCreateDoubleScalar(0.0));
  emlrtSetField(xEntryPoints, 0, "NumberOfOutputs", emlrtMxCreateDoubleScalar
                (1.0));
  emlrtSetField(xEntryPoints, 0, "ConstantInputs", xInputs);
  emlrtSetField(xEntryPoints, 0, "FullPath", emlrtMxCreateString(
    "D:\\Documents\\GitHub\\kinoRRTstar\\allExamples\\2DDI\\PFFCircle\\getTime.m"));
  emlrtSetField(xEntryPoints, 0, "TimeStamp", emlrtMxCreateDoubleScalar
                (738209.97052083339));
  xResult = emlrtCreateStructMatrix(1, 1, 4, b_fldNames);
  emlrtSetField(xResult, 0, "Version", emlrtMxCreateString(
    "9.5.0.1033004 (R2018b) Update 2"));
  emlrtSetField(xResult, 0, "ResolvedFunctions", (mxArray *)
                emlrtMexFcnResolvedFunctionsInfo());
  emlrtSetField(xResult, 0, "EntryPoints", xEntryPoints);
  return xResult;
}

/*
 * File trailer for _coder_getTime_info.c
 *
 * [EOF]
 */
