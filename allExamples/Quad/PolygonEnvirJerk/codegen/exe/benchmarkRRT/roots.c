/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * roots.c
 *
 * Code generation for function 'roots'
 *
 */

/* Include files */
#include <math.h>
#include <string.h>
#include "rt_nonfinite.h"
#include "benchmarkRRT.h"
#include "roots.h"
#include "eig.h"

/* Function Definitions */
void b_roots(const double c[7], creal_T r_data[], int r_size[1])
{
  int k1;
  int k2;
  int nTrailingZeros;
  int companDim;
  boolean_T exitg1;
  int j;
  int a_size[2];
  boolean_T exitg2;
  creal_T a_data[36];
  double ctmp[7];
  creal_T eiga_data[8];
  int eiga_size[1];
  memset(&r_data[0], 0, 6U * sizeof(creal_T));
  k1 = 1;
  while ((k1 <= 7) && (!(c[k1 - 1] != 0.0))) {
    k1++;
  }

  k2 = 7;
  while ((k2 >= k1) && (!(c[k2 - 1] != 0.0))) {
    k2--;
  }

  nTrailingZeros = 6 - k2;
  if (k1 < k2) {
    companDim = k2 - k1;
    exitg1 = false;
    while ((!exitg1) && (companDim > 0)) {
      j = 0;
      exitg2 = false;
      while ((!exitg2) && (j + 1 <= companDim)) {
        ctmp[j] = c[k1 + j] / c[k1 - 1];
        if (rtIsInf(fabs(ctmp[j]))) {
          exitg2 = true;
        } else {
          j++;
        }
      }

      if (j + 1 > companDim) {
        exitg1 = true;
      } else {
        k1++;
        companDim--;
      }
    }

    if (companDim < 1) {
      if (1 > 7 - k2) {
        r_size[0] = 0;
      } else {
        r_size[0] = 7 - k2;
      }
    } else {
      a_size[0] = companDim;
      a_size[1] = companDim;
      memset(&a_data[0], 0, (unsigned int)(companDim * companDim * (int)sizeof
              (creal_T)));
      for (k1 = 0; k1 <= companDim - 2; k1++) {
        j = companDim * k1;
        a_data[j].re = -ctmp[k1];
        a_data[j].im = 0.0;
        j = (k1 + j) + 1;
        a_data[j].re = 1.0;
        a_data[j].im = 0.0;
      }

      j = companDim * (companDim - 1);
      a_data[j].re = -ctmp[companDim - 1];
      a_data[j].im = 0.0;
      for (k1 = 0; k1 <= nTrailingZeros; k1++) {
        r_data[k1].re = 0.0;
        r_data[k1].im = 0.0;
      }

      eig(a_data, a_size, eiga_data, eiga_size);
      for (k1 = 0; k1 < companDim; k1++) {
        r_data[(k1 - k2) + 7] = eiga_data[k1];
      }

      r_size[0] = (companDim - k2) + 7;
    }
  } else if (1 > 7 - k2) {
    r_size[0] = 0;
  } else {
    r_size[0] = 7 - k2;
  }
}

void roots(const double c[9], creal_T r_data[], int r_size[1])
{
  int k1;
  int k2;
  int nTrailingZeros;
  int companDim;
  boolean_T exitg1;
  int j;
  int a_size[2];
  boolean_T exitg2;
  creal_T a_data[64];
  double ctmp[9];
  creal_T eiga_data[8];
  int eiga_size[1];
  memset(&r_data[0], 0, sizeof(creal_T) << 3);
  k1 = 1;
  while ((k1 <= 9) && (!(c[k1 - 1] != 0.0))) {
    k1++;
  }

  k2 = 9;
  while ((k2 >= k1) && (!(c[k2 - 1] != 0.0))) {
    k2--;
  }

  nTrailingZeros = 8 - k2;
  if (k1 < k2) {
    companDim = k2 - k1;
    exitg1 = false;
    while ((!exitg1) && (companDim > 0)) {
      j = 0;
      exitg2 = false;
      while ((!exitg2) && (j + 1 <= companDim)) {
        ctmp[j] = c[k1 + j] / c[k1 - 1];
        if (rtIsInf(fabs(ctmp[j]))) {
          exitg2 = true;
        } else {
          j++;
        }
      }

      if (j + 1 > companDim) {
        exitg1 = true;
      } else {
        k1++;
        companDim--;
      }
    }

    if (companDim < 1) {
      if (1 > 9 - k2) {
        r_size[0] = 0;
      } else {
        r_size[0] = 9 - k2;
      }
    } else {
      a_size[0] = companDim;
      a_size[1] = companDim;
      memset(&a_data[0], 0, (unsigned int)(companDim * companDim * (int)sizeof
              (creal_T)));
      for (k1 = 0; k1 <= companDim - 2; k1++) {
        j = companDim * k1;
        a_data[j].re = -ctmp[k1];
        a_data[j].im = 0.0;
        j = (k1 + j) + 1;
        a_data[j].re = 1.0;
        a_data[j].im = 0.0;
      }

      j = companDim * (companDim - 1);
      a_data[j].re = -ctmp[companDim - 1];
      a_data[j].im = 0.0;
      for (k1 = 0; k1 <= nTrailingZeros; k1++) {
        r_data[k1].re = 0.0;
        r_data[k1].im = 0.0;
      }

      eig(a_data, a_size, eiga_data, eiga_size);
      for (k1 = 0; k1 < companDim; k1++) {
        r_data[(k1 - k2) + 9] = eiga_data[k1];
      }

      r_size[0] = (companDim - k2) + 9;
    }
  } else if (1 > 9 - k2) {
    r_size[0] = 0;
  } else {
    r_size[0] = 9 - k2;
  }
}

/* End of code generation (roots.c) */
