/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * matlab_roots.cpp
 *
 * Code generation for function 'matlab_roots'
 *
 */

/* Include files */
#include <cmath>
#include "rt_nonfinite.h"
#include "matlab_roots.h"
#include "matlab_roots_emxutil.h"
#include "xzhseqr.h"
#include "xgehrd.h"
#include "anyNonFinite.h"
#include "xzgeev.h"

/* Function Definitions */
void matlab_roots(const emxArray_real_T *p, emxArray_creal_T *res)
{
  int istart;
  int k2;
  int k1;
  int nTrailingZeros;
  emxArray_real_T *ctmp;
  int companDim;
  boolean_T exitg1;
  int j;
  emxArray_creal_T *a;
  boolean_T exitg2;
  emxArray_creal_T *eiga;
  boolean_T b_p;
  emxArray_creal_T *beta1;
  int exitg3;
  double a_re;
  double eiga_re;
  int jend;
  double a_im;
  double eiga_im;
  boolean_T b_a;
  double beta1_re;
  double beta1_im;
  double brm;
  istart = res->size[0];
  res->size[0] = p->size[1] - 1;
  emxEnsureCapacity_creal_T(res, istart);
  k2 = p->size[1];
  for (istart = 0; istart <= k2 - 2; istart++) {
    res->data[istart].re = 0.0;
    res->data[istart].im = 0.0;
  }

  k1 = 1;
  while ((k1 <= p->size[1]) && (!(p->data[k1 - 1] != 0.0))) {
    k1++;
  }

  k2 = p->size[1];
  while ((k2 >= k1) && (!(p->data[k2 - 1] != 0.0))) {
    k2--;
  }

  nTrailingZeros = p->size[1] - k2;
  if (k1 < k2) {
    emxInit_real_T(&ctmp, 2);
    companDim = k2 - k1;
    istart = ctmp->size[0] * ctmp->size[1];
    ctmp->size[0] = p->size[0];
    ctmp->size[1] = p->size[1];
    emxEnsureCapacity_real_T(ctmp, istart);
    exitg1 = false;
    while ((!exitg1) && (companDim > 0)) {
      j = 0;
      exitg2 = false;
      while ((!exitg2) && (j + 1 <= companDim)) {
        ctmp->data[j] = p->data[k1 + j] / p->data[k1 - 1];
        if (rtIsInf(std::abs(ctmp->data[j]))) {
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
      if (1 > nTrailingZeros) {
        res->size[0] = 0;
      } else {
        istart = res->size[0];
        res->size[0] = nTrailingZeros;
        emxEnsureCapacity_creal_T(res, istart);
      }
    } else {
      emxInit_creal_T(&a, 2);
      istart = a->size[0] * a->size[1];
      a->size[0] = companDim;
      a->size[1] = companDim;
      emxEnsureCapacity_creal_T(a, istart);
      k2 = companDim * companDim;
      for (istart = 0; istart < k2; istart++) {
        a->data[istart].re = 0.0;
        a->data[istart].im = 0.0;
      }

      for (k2 = 0; k2 <= companDim - 2; k2++) {
        a->data[a->size[0] * k2].re = -ctmp->data[k2];
        a->data[a->size[0] * k2].im = 0.0;
        a->data[(k2 + a->size[0] * k2) + 1].re = 1.0;
        a->data[(k2 + a->size[0] * k2) + 1].im = 0.0;
      }

      a->data[a->size[0] * (companDim - 1)].re = -ctmp->data[companDim - 1];
      a->data[a->size[0] * (companDim - 1)].im = 0.0;
      for (k2 = 0; k2 < nTrailingZeros; k2++) {
        res->data[k2].re = 0.0;
        res->data[k2].im = 0.0;
      }

      emxInit_creal_T(&eiga, 1);
      if (anyNonFinite(a)) {
        if ((a->size[0] == 1) && (a->size[1] == 1)) {
          istart = eiga->size[0];
          eiga->size[0] = 1;
          emxEnsureCapacity_creal_T(eiga, istart);
          eiga->data[0].re = rtNaN;
          eiga->data[0].im = 0.0;
        } else {
          istart = eiga->size[0];
          eiga->size[0] = a->size[0];
          emxEnsureCapacity_creal_T(eiga, istart);
          k2 = a->size[0];
          for (istart = 0; istart < k2; istart++) {
            eiga->data[istart].re = rtNaN;
            eiga->data[istart].im = 0.0;
          }
        }
      } else if ((a->size[0] == 1) && (a->size[1] == 1)) {
        istart = eiga->size[0];
        eiga->size[0] = 1;
        emxEnsureCapacity_creal_T(eiga, istart);
        eiga->data[0] = a->data[0];
      } else {
        b_p = (a->size[0] == a->size[1]);
        if (b_p) {
          j = 0;
          exitg1 = false;
          while ((!exitg1) && (j <= a->size[1] - 1)) {
            k2 = 0;
            do {
              exitg3 = 0;
              if (k2 <= j) {
                a_re = a->data[j + a->size[0] * k2].re;
                a_im = -a->data[j + a->size[0] * k2].im;
                b_a = ((a->data[k2 + a->size[0] * j].re == a_re) && (a->data[k2
                        + a->size[0] * j].im == a_im));
                if (!b_a) {
                  b_p = false;
                  exitg3 = 1;
                } else {
                  k2++;
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
        }

        if (b_p) {
          if (anyNonFinite(a)) {
            k2 = a->size[0];
            k1 = a->size[1];
            istart = a->size[0] * a->size[1];
            a->size[0] = k2;
            a->size[1] = k1;
            emxEnsureCapacity_creal_T(a, istart);
            k2 *= k1;
            for (istart = 0; istart < k2; istart++) {
              a->data[istart].re = rtNaN;
              a->data[istart].im = 0.0;
            }

            k1 = a->size[0];
            if (1 < a->size[0]) {
              istart = 2;
              if (a->size[0] - 2 < a->size[1] - 1) {
                jend = a->size[0] - 1;
              } else {
                jend = a->size[1];
              }

              for (j = 0; j < jend; j++) {
                for (k2 = istart; k2 <= k1; k2++) {
                  a->data[(k2 + a->size[0] * j) - 1].re = 0.0;
                  a->data[(k2 + a->size[0] * j) - 1].im = 0.0;
                }

                istart++;
              }
            }
          } else {
            xgehrd(a);
            eml_zlahqr(a);
            k1 = a->size[0];
            if (3 < a->size[0]) {
              istart = 4;
              if (a->size[0] - 4 < a->size[1] - 1) {
                jend = a->size[0] - 3;
              } else {
                jend = a->size[1];
              }

              for (j = 0; j < jend; j++) {
                for (k2 = istart; k2 <= k1; k2++) {
                  a->data[(k2 + a->size[0] * j) - 1].re = 0.0;
                  a->data[(k2 + a->size[0] * j) - 1].im = 0.0;
                }

                istart++;
              }
            }
          }

          k1 = a->size[0];
          istart = eiga->size[0];
          eiga->size[0] = a->size[0];
          emxEnsureCapacity_creal_T(eiga, istart);
          for (k2 = 0; k2 < k1; k2++) {
            eiga->data[k2] = a->data[k2 + a->size[0] * k2];
          }
        } else {
          emxInit_creal_T(&beta1, 1);
          xzgeev(a, &k2, eiga, beta1);
          istart = eiga->size[0];
          emxEnsureCapacity_creal_T(eiga, istart);
          k2 = eiga->size[0];
          for (istart = 0; istart < k2; istart++) {
            eiga_re = eiga->data[istart].re;
            eiga_im = eiga->data[istart].im;
            beta1_re = beta1->data[istart].re;
            beta1_im = beta1->data[istart].im;
            if (beta1_im == 0.0) {
              if (eiga_im == 0.0) {
                eiga->data[istart].re = eiga_re / beta1_re;
                eiga->data[istart].im = 0.0;
              } else if (eiga_re == 0.0) {
                eiga->data[istart].re = 0.0;
                eiga->data[istart].im = eiga_im / beta1_re;
              } else {
                eiga->data[istart].re = eiga_re / beta1_re;
                eiga->data[istart].im = eiga_im / beta1_re;
              }
            } else if (beta1_re == 0.0) {
              if (eiga_re == 0.0) {
                eiga->data[istart].re = eiga_im / beta1_im;
                eiga->data[istart].im = 0.0;
              } else if (eiga_im == 0.0) {
                eiga->data[istart].re = 0.0;
                eiga->data[istart].im = -(eiga_re / beta1_im);
              } else {
                eiga->data[istart].re = eiga_im / beta1_im;
                eiga->data[istart].im = -(eiga_re / beta1_im);
              }
            } else {
              brm = std::abs(beta1_re);
              a_re = std::abs(beta1_im);
              if (brm > a_re) {
                a_im = beta1_im / beta1_re;
                a_re = beta1_re + a_im * beta1_im;
                eiga->data[istart].re = (eiga_re + a_im * eiga_im) / a_re;
                eiga->data[istart].im = (eiga_im - a_im * eiga_re) / a_re;
              } else if (a_re == brm) {
                if (beta1_re > 0.0) {
                  a_im = 0.5;
                } else {
                  a_im = -0.5;
                }

                if (beta1_im > 0.0) {
                  a_re = 0.5;
                } else {
                  a_re = -0.5;
                }

                eiga->data[istart].re = (eiga_re * a_im + eiga_im * a_re) / brm;
                eiga->data[istart].im = (eiga_im * a_im - eiga_re * a_re) / brm;
              } else {
                a_im = beta1_re / beta1_im;
                a_re = beta1_im + a_im * beta1_re;
                eiga->data[istart].re = (a_im * eiga_re + eiga_im) / a_re;
                eiga->data[istart].im = (a_im * eiga_im - eiga_re) / a_re;
              }
            }
          }

          emxFree_creal_T(&beta1);
        }
      }

      emxFree_creal_T(&a);
      for (k2 = 0; k2 < companDim; k2++) {
        res->data[k2 + nTrailingZeros] = eiga->data[k2];
      }

      emxFree_creal_T(&eiga);
      k2 = nTrailingZeros + companDim;
      if (1 > k2) {
        res->size[0] = 0;
      } else {
        istart = res->size[0];
        res->size[0] = k2;
        emxEnsureCapacity_creal_T(res, istart);
      }
    }

    emxFree_real_T(&ctmp);
  } else if (1 > nTrailingZeros) {
    res->size[0] = 0;
  } else {
    istart = res->size[0];
    res->size[0] = nTrailingZeros;
    emxEnsureCapacity_creal_T(res, istart);
  }
}

/* End of code generation (matlab_roots.cpp) */
