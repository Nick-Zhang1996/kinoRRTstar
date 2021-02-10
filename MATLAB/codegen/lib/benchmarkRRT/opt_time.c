/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * opt_time.c
 *
 * Code generation for function 'opt_time'
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "benchmarkRRT.h"
#include "opt_time.h"

/* Function Definitions */
void opt_time(double x01, double x02, double x03, double x04, double x11, double
              x12, double x13, double x14, double time_equ[5])
{
  time_equ[0] = 1.0;
  time_equ[1] = 0.0;
  time_equ[2] = ((((-4.0 * (x03 * x03) - 4.0 * x03 * x13) - 4.0 * (x04 * x04)) -
                  4.0 * x04 * x14) - 4.0 * (x13 * x13)) - 4.0 * (x14 * x14);
  time_equ[3] = ((((((24.0 * x03 * x11 - 24.0 * x02 * x04) - 24.0 * x01 * x13) -
                    24.0 * x01 * x03) - 24.0 * x02 * x14) + 24.0 * x04 * x12) +
                 24.0 * x11 * x13) + 24.0 * x12 * x14;
  time_equ[4] = ((((-36.0 * (x01 * x01) + 72.0 * x01 * x11) - 36.0 * (x02 * x02))
                  + 72.0 * x02 * x12) - 36.0 * (x11 * x11)) - 36.0 * (x12 * x12);
}

/* End of code generation (opt_time.c) */
