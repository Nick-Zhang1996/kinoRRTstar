/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * eig.c
 *
 * Code generation for function 'eig'
 *
 */

/* Include files */
#include <math.h>
#include "rt_nonfinite.h"
#include "benchmarkRRT.h"
#include "eig.h"
#include "schur.h"
#include "xzgeev.h"
#include "anyNonFinite.h"

/* Function Definitions */
void eig(const creal_T A_data[], const int A_size[2], creal_T V_data[], int
         V_size[1])
{
  boolean_T p;
  int info;
  int i;
  boolean_T exitg2;
  creal_T beta1_data[8];
  int beta1_size[1];
  creal_T T_data[64];
  int T_size[2];
  double V_data_re;
  int exitg1;
  double brm;
  double bim;
  double d;
  if (anyNonFinite(A_data, A_size)) {
    if ((A_size[0] == 1) && (A_size[1] == 1)) {
      V_size[0] = 1;
      V_data[0].re = rtNaN;
      V_data[0].im = 0.0;
    } else {
      V_size[0] = A_size[0];
      info = A_size[0];
      for (i = 0; i < info; i++) {
        V_data[i].re = rtNaN;
        V_data[i].im = 0.0;
      }
    }
  } else if ((A_size[0] == 1) && (A_size[1] == 1)) {
    V_size[0] = 1;
    V_data[0] = A_data[0];
  } else {
    p = (A_size[0] == A_size[1]);
    if (p) {
      info = 0;
      exitg2 = false;
      while ((!exitg2) && (info <= A_size[1] - 1)) {
        i = 0;
        do {
          exitg1 = 0;
          if (i <= info) {
            if ((!(A_data[i + A_size[0] * info].re == A_data[info + A_size[0] *
                   i].re)) || (!(A_data[i + A_size[0] * info].im == -A_data[info
                                 + A_size[0] * i].im))) {
              p = false;
              exitg1 = 1;
            } else {
              i++;
            }
          } else {
            info++;
            exitg1 = 2;
          }
        } while (exitg1 == 0);

        if (exitg1 == 1) {
          exitg2 = true;
        }
      }
    }

    if (p) {
      schur(A_data, A_size, T_data, T_size);
      info = T_size[0];
      V_size[0] = T_size[0];
      for (i = 0; i < info; i++) {
        V_data[i] = T_data[i + T_size[0] * i];
      }
    } else {
      xzgeev(A_data, A_size, &info, V_data, V_size, beta1_data, beta1_size);
      info = V_size[0];
      for (i = 0; i < info; i++) {
        V_data_re = V_data[i].re;
        if (beta1_data[i].im == 0.0) {
          if (V_data[i].im == 0.0) {
            V_data[i].re /= beta1_data[i].re;
            V_data[i].im = 0.0;
          } else if (V_data[i].re == 0.0) {
            V_data[i].re = 0.0;
            V_data[i].im /= beta1_data[i].re;
          } else {
            V_data[i].re /= beta1_data[i].re;
            V_data[i].im /= beta1_data[i].re;
          }
        } else if (beta1_data[i].re == 0.0) {
          if (V_data[i].re == 0.0) {
            V_data[i].re = V_data[i].im / beta1_data[i].im;
            V_data[i].im = 0.0;
          } else if (V_data[i].im == 0.0) {
            V_data[i].re = 0.0;
            V_data[i].im = -(V_data_re / beta1_data[i].im);
          } else {
            V_data[i].re = V_data[i].im / beta1_data[i].im;
            V_data[i].im = -(V_data_re / beta1_data[i].im);
          }
        } else {
          brm = fabs(beta1_data[i].re);
          bim = fabs(beta1_data[i].im);
          if (brm > bim) {
            bim = beta1_data[i].im / beta1_data[i].re;
            d = beta1_data[i].re + bim * beta1_data[i].im;
            V_data[i].re = (V_data[i].re + bim * V_data[i].im) / d;
            V_data[i].im = (V_data[i].im - bim * V_data_re) / d;
          } else if (bim == brm) {
            if (beta1_data[i].re > 0.0) {
              bim = 0.5;
            } else {
              bim = -0.5;
            }

            if (beta1_data[i].im > 0.0) {
              d = 0.5;
            } else {
              d = -0.5;
            }

            V_data[i].re = (V_data[i].re * bim + V_data[i].im * d) / brm;
            V_data[i].im = (V_data[i].im * bim - V_data_re * d) / brm;
          } else {
            bim = beta1_data[i].re / beta1_data[i].im;
            d = beta1_data[i].im + bim * beta1_data[i].re;
            V_data[i].re = (bim * V_data[i].re + V_data[i].im) / d;
            V_data[i].im = (bim * V_data[i].im - V_data_re) / d;
          }
        }
      }
    }
  }
}

/* End of code generation (eig.c) */
