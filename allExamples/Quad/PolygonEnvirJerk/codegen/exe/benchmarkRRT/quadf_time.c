/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * quadf_time.c
 *
 * Code generation for function 'quadf_time'
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "benchmarkRRT.h"
#include "quadf_time.h"

/* Function Definitions */
void quadf_time(double x01, double x02, double x03, double x04, double x05,
                double x06, double x07, double x08, double x09, double x11,
                double x12, double x13, double x14, double x15, double x16,
                double x17, double x18, double x19, double time_equ[7])
{
  time_equ[0] = 1.0;
  time_equ[1] = 0.0;
  time_equ[2] = (((((((3.0 * x07 * x17 / 2500.0 + 3.0 * x08 * x18 / 2500.0) +
                      3.0 * x09 * x19 / 2500.0) - 9.0 * (x07 * x07) / 5000.0) -
                    9.0 * (x08 * x08) / 5000.0) - 9.0 * (x09 * x09) / 5000.0) -
                  9.0 * (x17 * x17) / 5000.0) - 9.0 * (x18 * x18) / 5000.0) -
    9.0 * (x19 * x19) / 5000.0;
  time_equ[3] = ((((((((((12.0 * x04 * x17 / 625.0 - 18.0 * x05 * x08 / 625.0) -
    18.0 * x06 * x09 / 625.0) - 18.0 * x04 * x07 / 625.0) - 12.0 * x07 * x14 /
                       625.0) + 12.0 * x05 * x18 / 625.0) - 12.0 * x08 * x15 /
                     625.0) + 12.0 * x06 * x19 / 625.0) - 12.0 * x09 * x16 /
                   625.0) + 18.0 * x14 * x17 / 625.0) + 18.0 * x15 * x18 / 625.0)
    + 18.0 * x16 * x19 / 625.0;
  time_equ[4] = (((((((((((((((((((9.0 * x01 * x17 / 125.0 - 9.0 * x02 * x08 /
    125.0) - 9.0 * x03 * x09 / 125.0) - 9.0 * x01 * x07 / 125.0) - 126.0 * x04 *
    x14 / 625.0) + 9.0 * x07 * x11 / 125.0) + 9.0 * x02 * x18 / 125.0) - 126.0 *
    x05 * x15 / 625.0) + 9.0 * x08 * x12 / 125.0) + 9.0 * x03 * x19 / 125.0) -
    126.0 * x06 * x16 / 625.0) + 9.0 * x09 * x13 / 125.0) - 9.0 * x11 * x17 /
                        125.0) - 9.0 * x12 * x18 / 125.0) - 9.0 * x13 * x19 /
                      125.0) - 72.0 * (x04 * x04) / 625.0) - 72.0 * (x05 * x05) /
                    625.0) - 72.0 * (x06 * x06) / 625.0) - 72.0 * (x14 * x14) /
                  625.0) - 72.0 * (x15 * x15) / 625.0) - 72.0 * (x16 * x16) /
    625.0;
  time_equ[5] = ((((((((((72.0 * x04 * x11 / 125.0 - 72.0 * x02 * x05 / 125.0) -
    72.0 * x03 * x06 / 125.0) - 72.0 * x01 * x14 / 125.0) - 72.0 * x01 * x04 /
                       125.0) - 72.0 * x02 * x15 / 125.0) + 72.0 * x05 * x12 /
                     125.0) - 72.0 * x03 * x16 / 125.0) + 72.0 * x06 * x13 /
                   125.0) + 72.0 * x11 * x14 / 125.0) + 72.0 * x12 * x15 / 125.0)
    + 72.0 * x13 * x16 / 125.0;
  time_equ[6] = (((((((36.0 * x01 * x11 / 25.0 + 36.0 * x02 * x12 / 25.0) + 36.0
                      * x03 * x13 / 25.0) - 18.0 * (x01 * x01) / 25.0) - 18.0 *
                    (x02 * x02) / 25.0) - 18.0 * (x03 * x03) / 25.0) - 18.0 *
                  (x11 * x11) / 25.0) - 18.0 * (x12 * x12) / 25.0) - 18.0 * (x13
    * x13) / 25.0;
}

/* End of code generation (quadf_time.c) */
