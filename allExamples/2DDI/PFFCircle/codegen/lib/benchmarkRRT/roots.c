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
#include "xzhseqr.h"
#include "xgehrd.h"
#include "anyNonFinite.h"
#include "xzgeev.h"

/* Function Definitions */
void roots(const double c[5], creal_T r_data[], int r_size[1])
{
  int k1;
  int k2;
  int nTrailingZeros;
  int companDim;
  boolean_T exitg1;
  int j;
  int a_size[2];
  boolean_T exitg2;
  creal_T a_data[16];
  int m;
  double ctmp[5];
  int i2;
  boolean_T p;
  creal_T eiga_data[4];
  int eiga_size[1];
  creal_T beta1_data[4];
  int beta1_size[1];
  int exitg3;
  double eiga_data_re;
  double brm;
  double bim;
  int jend;
  double d;
  memset(&r_data[0], 0, sizeof(creal_T) << 2);
  k1 = 1;
  while ((k1 <= 5) && (!(c[k1 - 1] != 0.0))) {
    k1++;
  }

  k2 = 5;
  while ((k2 >= k1) && (!(c[k2 - 1] != 0.0))) {
    k2--;
  }

  nTrailingZeros = 4 - k2;
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
      if (1 > 5 - k2) {
        r_size[0] = 0;
      } else {
        r_size[0] = 5 - k2;
      }
    } else {
      a_size[0] = companDim;
      a_size[1] = companDim;
      memset(&a_data[0], 0, (unsigned int)(companDim * companDim * (int)sizeof
              (creal_T)));
      for (m = 0; m <= companDim - 2; m++) {
        i2 = companDim * m;
        a_data[i2].re = -ctmp[m];
        a_data[i2].im = 0.0;
        i2 = (m + i2) + 1;
        a_data[i2].re = 1.0;
        a_data[i2].im = 0.0;
      }

      i2 = companDim * (companDim - 1);
      a_data[i2].re = -ctmp[companDim - 1];
      a_data[i2].im = 0.0;
      for (m = 0; m <= nTrailingZeros; m++) {
        r_data[m].re = 0.0;
        r_data[m].im = 0.0;
      }

      if (anyNonFinite(a_data, a_size)) {
        if (companDim == 1) {
          eiga_data[0].re = rtNaN;
          eiga_data[0].im = 0.0;
        } else {
          for (i2 = 0; i2 < companDim; i2++) {
            eiga_data[i2].re = rtNaN;
            eiga_data[i2].im = 0.0;
          }
        }
      } else if (companDim == 1) {
        eiga_data[0] = a_data[0];
      } else {
        p = true;
        j = 0;
        exitg1 = false;
        while ((!exitg1) && (j <= companDim - 1)) {
          k1 = 0;
          do {
            exitg3 = 0;
            if (k1 <= j) {
              if ((!(a_data[k1 + companDim * j].re == a_data[j + companDim * k1]
                     .re)) || (!(a_data[k1 + companDim * j].im == -a_data[j +
                                 companDim * k1].im))) {
                p = false;
                exitg3 = 1;
              } else {
                k1++;
              }
            } else {
              j++;
              exitg3 = 2;
            }
          } while (exitg3 == 0);

          if (exitg3 == 1) {
            exitg1 = true;
          }
        }

        if (p) {
          if (anyNonFinite(a_data, a_size)) {
            a_size[0] = (signed char)companDim;
            k1 = (signed char)companDim * (signed char)companDim;
            for (i2 = 0; i2 < k1; i2++) {
              a_data[i2].re = rtNaN;
              a_data[i2].im = 0.0;
            }

            m = (signed char)companDim;
            if (1 < (signed char)companDim) {
              nTrailingZeros = 2;
              if ((signed char)companDim - 2 < (signed char)companDim - 1) {
                jend = (signed char)companDim - 1;
              } else {
                jend = (signed char)companDim;
              }

              for (j = 0; j < jend; j++) {
                for (k1 = nTrailingZeros; k1 <= m; k1++) {
                  i2 = (k1 + (signed char)companDim * j) - 1;
                  a_data[i2].re = 0.0;
                  a_data[i2].im = 0.0;
                }

                nTrailingZeros++;
              }
            }
          } else {
            xgehrd(a_data, a_size);
            eml_zlahqr(a_data, a_size);
            if (3 < a_size[0]) {
              for (k1 = 4; k1 < 5; k1++) {
                a_data[3].re = 0.0;
                a_data[3].im = 0.0;
              }
            }
          }

          k1 = a_size[0];
          for (m = 0; m < k1; m++) {
            eiga_data[m] = a_data[m + a_size[0] * m];
          }
        } else {
          xzgeev(a_data, a_size, &k1, eiga_data, eiga_size, beta1_data,
                 beta1_size);
          k1 = eiga_size[0];
          for (i2 = 0; i2 < k1; i2++) {
            eiga_data_re = eiga_data[i2].re;
            if (beta1_data[i2].im == 0.0) {
              if (eiga_data[i2].im == 0.0) {
                eiga_data[i2].re /= beta1_data[i2].re;
                eiga_data[i2].im = 0.0;
              } else if (eiga_data[i2].re == 0.0) {
                eiga_data[i2].re = 0.0;
                eiga_data[i2].im /= beta1_data[i2].re;
              } else {
                eiga_data[i2].re /= beta1_data[i2].re;
                eiga_data[i2].im /= beta1_data[i2].re;
              }
            } else if (beta1_data[i2].re == 0.0) {
              if (eiga_data[i2].re == 0.0) {
                eiga_data[i2].re = eiga_data[i2].im / beta1_data[i2].im;
                eiga_data[i2].im = 0.0;
              } else if (eiga_data[i2].im == 0.0) {
                eiga_data[i2].re = 0.0;
                eiga_data[i2].im = -(eiga_data_re / beta1_data[i2].im);
              } else {
                eiga_data[i2].re = eiga_data[i2].im / beta1_data[i2].im;
                eiga_data[i2].im = -(eiga_data_re / beta1_data[i2].im);
              }
            } else {
              brm = fabs(beta1_data[i2].re);
              bim = fabs(beta1_data[i2].im);
              if (brm > bim) {
                bim = beta1_data[i2].im / beta1_data[i2].re;
                d = beta1_data[i2].re + bim * beta1_data[i2].im;
                eiga_data[i2].re = (eiga_data[i2].re + bim * eiga_data[i2].im) /
                  d;
                eiga_data[i2].im = (eiga_data[i2].im - bim * eiga_data_re) / d;
              } else if (bim == brm) {
                if (beta1_data[i2].re > 0.0) {
                  bim = 0.5;
                } else {
                  bim = -0.5;
                }

                if (beta1_data[i2].im > 0.0) {
                  d = 0.5;
                } else {
                  d = -0.5;
                }

                eiga_data[i2].re = (eiga_data[i2].re * bim + eiga_data[i2].im *
                                    d) / brm;
                eiga_data[i2].im = (eiga_data[i2].im * bim - eiga_data_re * d) /
                  brm;
              } else {
                bim = beta1_data[i2].re / beta1_data[i2].im;
                d = beta1_data[i2].im + bim * beta1_data[i2].re;
                eiga_data[i2].re = (bim * eiga_data[i2].re + eiga_data[i2].im) /
                  d;
                eiga_data[i2].im = (bim * eiga_data[i2].im - eiga_data_re) / d;
              }
            }
          }
        }
      }

      for (m = 0; m < companDim; m++) {
        r_data[(m - k2) + 5] = eiga_data[m];
      }

      r_size[0] = (companDim - k2) + 5;
    }
  } else if (1 > 5 - k2) {
    r_size[0] = 0;
  } else {
    r_size[0] = 5 - k2;
  }
}

/* End of code generation (roots.c) */
