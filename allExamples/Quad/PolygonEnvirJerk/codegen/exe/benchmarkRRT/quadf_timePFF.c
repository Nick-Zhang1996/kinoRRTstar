/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * quadf_timePFF.c
 *
 * Code generation for function 'quadf_timePFF'
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "benchmarkRRT.h"
#include "quadf_timePFF.h"

/* Function Definitions */
void quadf_timePFF(double x01, double x02, double x03, double x04, double x05,
                   double x06, double x07, double x08, double x09, double x11,
                   double x12, double x13, double time_equ[9])
{
  double time_equ_tmp;
  double b_time_equ_tmp;
  double c_time_equ_tmp;
  double d_time_equ_tmp;
  double e_time_equ_tmp;
  double f_time_equ_tmp;
  double g_time_equ_tmp;
  double h_time_equ_tmp;
  double i_time_equ_tmp;
  time_equ[0] = 2.5E+8;
  time_equ[1] = 9.0E+6;
  time_equ_tmp = x07 * x07;
  b_time_equ_tmp = x08 * x08;
  c_time_equ_tmp = x09 * x09;
  time_equ[2] = ((-300000.0 * time_equ_tmp - 300000.0 * b_time_equ_tmp) -
                 300000.0 * c_time_equ_tmp) + 81000.0;
  time_equ[3] = ((((-9000.0 * time_equ_tmp - 3.0E+6 * x04 * x07) - 9000.0 *
                   b_time_equ_tmp) - 3.0E+6 * x05 * x08) - 9000.0 *
                 c_time_equ_tmp) - 3.0E+6 * x06 * x09;
  d_time_equ_tmp = x04 * x04;
  e_time_equ_tmp = x05 * x05;
  f_time_equ_tmp = x06 * x06;
  time_equ[4] = (((((((((((((4.5E+6 * x07 * x11 - 4.5E+6 * x02 * x08) - 81000.0 *
    x04 * x07) - 4.5E+6 * x03 * x09) - 81000.0 * x05 * x08) - 81000.0 * x06 *
    x09) - 4.5E+6 * x01 * x07) + 4.5E+6 * x08 * x12) + 4.5E+6 * x09 * x13) -
                     6.75E+6 * d_time_equ_tmp) - 6.75E+6 * e_time_equ_tmp) -
                   6.75E+6 * f_time_equ_tmp) - 81.0 * time_equ_tmp) - 81.0 *
                 b_time_equ_tmp) - 81.0 * c_time_equ_tmp;
  time_equ[5] = ((((((((((((((((1.8E+7 * x04 * x11 - 1.8E+7 * x02 * x05) -
    126000.0 * x01 * x07) - 1.8E+7 * x03 * x06) - 126000.0 * x02 * x08) - 648.0 *
    x04 * x07) - 126000.0 * x03 * x09) - 648.0 * x05 * x08) - 1.8E+7 * x01 * x04)
                        - 648.0 * x06 * x09) + 1.8E+7 * x05 * x12) + 126000.0 *
                      x07 * x11) + 1.8E+7 * x06 * x13) + 126000.0 * x08 * x12) +
                   126000.0 * x09 * x13) - 153000.0 * d_time_equ_tmp) - 153000.0
                 * e_time_equ_tmp) - 153000.0 * f_time_equ_tmp;
  time_equ_tmp = x01 * x01;
  b_time_equ_tmp = x02 * x02;
  c_time_equ_tmp = x03 * x03;
  g_time_equ_tmp = x11 * x11;
  h_time_equ_tmp = x12 * x12;
  i_time_equ_tmp = x13 * x13;
  time_equ[6] = ((((((((((((((((((((((-1.125E+7 * time_equ_tmp - 423000.0 * x01 *
    x04) + 2.25E+7 * x01 * x11) - 972.0 * x07 * x01) - 1.125E+7 * b_time_equ_tmp)
    - 423000.0 * x02 * x05) + 2.25E+7 * x02 * x12) - 972.0 * x08 * x02) -
    1.125E+7 * c_time_equ_tmp) - 423000.0 * x03 * x06) + 2.25E+7 * x03 * x13) -
    972.0 * x09 * x03) - 972.0 * d_time_equ_tmp) + 423000.0 * x04 * x11) - 972.0
    * e_time_equ_tmp) + 423000.0 * x05 * x12) - 972.0 * f_time_equ_tmp) +
                      423000.0 * x06 * x13) - 1.125E+7 * g_time_equ_tmp) + 972.0
                    * x07 * x11) - 1.125E+7 * h_time_equ_tmp) + 972.0 * x08 *
                  x12) - 1.125E+7 * i_time_equ_tmp) + 972.0 * x09 * x13;
  time_equ[7] = (((((((((((((-270000.0 * time_equ_tmp + 540000.0 * x01 * x11) -
    2592.0 * x04 * x01) - 270000.0 * b_time_equ_tmp) + 540000.0 * x02 * x12) -
    2592.0 * x05 * x02) - 270000.0 * c_time_equ_tmp) + 540000.0 * x03 * x13) -
                      2592.0 * x06 * x03) - 270000.0 * g_time_equ_tmp) + 2592.0 *
                    x04 * x11) - 270000.0 * h_time_equ_tmp) + 2592.0 * x05 * x12)
                 - 270000.0 * i_time_equ_tmp) + 2592.0 * x06 * x13;
  time_equ[8] = (((((((-1620.0 * time_equ_tmp + 3240.0 * x01 * x11) - 1620.0 *
                      b_time_equ_tmp) + 3240.0 * x02 * x12) - 1620.0 *
                    c_time_equ_tmp) + 3240.0 * x03 * x13) - 1620.0 *
                  g_time_equ_tmp) - 1620.0 * h_time_equ_tmp) - 1620.0 *
    i_time_equ_tmp;
}

/* End of code generation (quadf_timePFF.c) */
