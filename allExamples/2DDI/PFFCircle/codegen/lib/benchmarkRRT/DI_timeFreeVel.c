/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * DI_timeFreeVel.c
 *
 * Code generation for function 'DI_timeFreeVel'
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "benchmarkRRT.h"
#include "DI_timeFreeVel.h"

/* Function Definitions */
void DI_timeFreeVel(double x01, double x02, double x03, double x04, double x11,
                    double x12, double time_equ[5])
{
  time_equ[0] = 1.0;
  time_equ[1] = 0.0;
  time_equ[2] = -3.0 * (x03 * x03) - 3.0 * (x04 * x04);
  time_equ[3] = ((12.0 * x03 * x11 - 12.0 * x02 * x04) - 12.0 * x01 * x03) +
    12.0 * x04 * x12;
  time_equ[4] = ((((-9.0 * (x01 * x01) + 18.0 * x01 * x11) - 9.0 * (x02 * x02))
                  + 18.0 * x02 * x12) - 9.0 * (x11 * x11)) - 9.0 * (x12 * x12);
}

/* End of code generation (DI_timeFreeVel.c) */
