/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * xgehrd.cpp
 *
 * Code generation for function 'xgehrd'
 *
 */

/* Include files */
#include <cmath>
#include "rt_nonfinite.h"
#include "matlab_roots.h"
#include "xgehrd.h"
#include "matlab_roots_emxutil.h"
#include "xscal.h"
#include "recip.h"
#include "xdlapy3.h"
#include "xnrm2.h"

/* Function Definitions */
void xgehrd(emxArray_creal_T *a)
{
  emxArray_creal_T *tau;
  emxArray_creal_T *work;
  int n;
  int i;
  int i0;
  int b_i;
  int im1n_tmp;
  int in;
  creal_T alpha1;
  int n_tmp;
  int lastc;
  double xnorm;
  double beta1;
  int iv0_tmp;
  boolean_T b_tau;
  int knt;
  double beta1_im;
  int i1;
  int lastv;
  int b_lastc;
  boolean_T exitg1;
  int k;
  creal_T c;
  double temp_im;
  int ix;
  int exitg5;
  int exitg2;
  int i2;
  int exitg4;
  emxInit_creal_T(&tau, 1);
  emxInit_creal_T(&work, 1);
  n = a->size[0];
  i = a->size[0] - 1;
  i0 = tau->size[0];
  tau->size[0] = i;
  emxEnsureCapacity_creal_T(tau, i0);
  i = a->size[0];
  i0 = work->size[0];
  work->size[0] = i;
  emxEnsureCapacity_creal_T(work, i0);
  for (i0 = 0; i0 < i; i0++) {
    work->data[i0].re = 0.0;
    work->data[i0].im = 0.0;
  }

  i0 = a->size[0];
  for (b_i = 0; b_i <= i0 - 2; b_i++) {
    im1n_tmp = b_i * n;
    in = (b_i + 1) * n;
    alpha1 = a->data[(b_i + a->size[0] * b_i) + 1];
    i = b_i + 3;
    if (i >= n) {
      i = n;
    }

    i += im1n_tmp;
    n_tmp = n - b_i;
    lastc = n_tmp - 2;
    tau->data[b_i].re = 0.0;
    tau->data[b_i].im = 0.0;
    if (lastc + 1 > 0) {
      xnorm = xnrm2(lastc, a, i);
      if ((xnorm != 0.0) || (a->data[(b_i + a->size[0] * b_i) + 1].im != 0.0)) {
        beta1 = xdlapy3(a->data[(b_i + a->size[0] * b_i) + 1].re, a->data[(b_i +
          a->size[0] * b_i) + 1].im, xnorm);
        if (a->data[(b_i + a->size[0] * b_i) + 1].re >= 0.0) {
          beta1 = -beta1;
        }

        if (std::abs(beta1) < 1.0020841800044864E-292) {
          knt = -1;
          i1 = (i + lastc) - 1;
          do {
            knt++;
            for (k = i; k <= i1; k++) {
              xnorm = a->data[k - 1].re;
              beta1_im = a->data[k - 1].im;
              a->data[k - 1].re = 9.9792015476736E+291 * xnorm - 0.0 * beta1_im;
              a->data[k - 1].im = 9.9792015476736E+291 * beta1_im + 0.0 * xnorm;
            }

            beta1 *= 9.9792015476736E+291;
            alpha1.re *= 9.9792015476736E+291;
            alpha1.im *= 9.9792015476736E+291;
          } while (!(std::abs(beta1) >= 1.0020841800044864E-292));

          beta1 = xdlapy3(alpha1.re, alpha1.im, xnrm2(lastc, a, i));
          if (alpha1.re >= 0.0) {
            beta1 = -beta1;
          }

          xnorm = beta1 - alpha1.re;
          if (0.0 - alpha1.im == 0.0) {
            tau->data[b_i].re = xnorm / beta1;
            tau->data[b_i].im = 0.0;
          } else if (xnorm == 0.0) {
            tau->data[b_i].re = 0.0;
            tau->data[b_i].im = (0.0 - alpha1.im) / beta1;
          } else {
            tau->data[b_i].re = xnorm / beta1;
            tau->data[b_i].im = (0.0 - alpha1.im) / beta1;
          }

          c.re = alpha1.re - beta1;
          c.im = alpha1.im;
          xscal(lastc, recip(c), a, i);
          for (k = 0; k <= knt; k++) {
            beta1 *= 1.0020841800044864E-292;
          }

          alpha1.re = beta1;
          alpha1.im = 0.0;
        } else {
          xnorm = beta1 - a->data[(b_i + a->size[0] * b_i) + 1].re;
          beta1_im = 0.0 - a->data[(b_i + a->size[0] * b_i) + 1].im;
          if (beta1_im == 0.0) {
            tau->data[b_i].re = xnorm / beta1;
            tau->data[b_i].im = 0.0;
          } else if (xnorm == 0.0) {
            tau->data[b_i].re = 0.0;
            tau->data[b_i].im = beta1_im / beta1;
          } else {
            tau->data[b_i].re = xnorm / beta1;
            tau->data[b_i].im = beta1_im / beta1;
          }

          c.re = a->data[(b_i + a->size[0] * b_i) + 1].re - beta1;
          c.im = a->data[(b_i + a->size[0] * b_i) + 1].im;
          xscal(lastc, recip(c), a, i);
          alpha1.re = beta1;
          alpha1.im = 0.0;
        }
      }
    }

    a->data[(b_i + a->size[0] * b_i) + 1].re = 1.0;
    a->data[(b_i + a->size[0] * b_i) + 1].im = 0.0;
    n_tmp -= 3;
    iv0_tmp = (b_i + im1n_tmp) + 2;
    im1n_tmp = in + 1;
    b_tau = ((tau->data[b_i].re != 0.0) || (tau->data[b_i].im != 0.0));
    if (b_tau) {
      lastv = n_tmp + 2;
      i = iv0_tmp + n_tmp;
      exitg1 = false;
      while ((!exitg1) && (lastv > 0)) {
        b_tau = ((a->data[i].re == 0.0) && (a->data[i].im == 0.0));
        if (b_tau) {
          lastv--;
          i--;
        } else {
          exitg1 = true;
        }
      }

      b_lastc = n;
      exitg1 = false;
      while ((!exitg1) && (b_lastc > 0)) {
        i = in + b_lastc;
        knt = i + (lastv - 1) * n;
        do {
          exitg2 = 0;
          if (((n > 0) && (i <= knt)) || ((n < 0) && (i >= knt))) {
            b_tau = ((a->data[i - 1].re != 0.0) || (a->data[i - 1].im != 0.0));
            if (b_tau) {
              exitg2 = 1;
            } else {
              i += n;
            }
          } else {
            b_lastc--;
            exitg2 = 2;
          }
        } while (exitg2 == 0);

        if (exitg2 == 1) {
          exitg1 = true;
        }
      }
    } else {
      lastv = 0;
      b_lastc = 0;
    }

    if (lastv > 0) {
      if (b_lastc != 0) {
        for (i = 0; i < b_lastc; i++) {
          work->data[i].re = 0.0;
          work->data[i].im = 0.0;
        }

        ix = iv0_tmp;
        i1 = (in + n * (lastv - 1)) + 1;
        for (knt = im1n_tmp; n < 0 ? knt >= i1 : knt <= i1; knt += n) {
          c.re = a->data[ix - 1].re - 0.0 * a->data[ix - 1].im;
          c.im = a->data[ix - 1].im + 0.0 * a->data[ix - 1].re;
          i = 0;
          i2 = (knt + b_lastc) - 1;
          for (k = knt; k <= i2; k++) {
            xnorm = a->data[k - 1].re * c.re - a->data[k - 1].im * c.im;
            beta1_im = a->data[k - 1].re * c.im + a->data[k - 1].im * c.re;
            work->data[i].re += xnorm;
            work->data[i].im += beta1_im;
            i++;
          }

          ix++;
        }
      }

      c.re = -tau->data[b_i].re;
      c.im = -tau->data[b_i].im;
      if ((!(c.re == 0.0)) || (!(c.im == 0.0))) {
        i = in;
        knt = iv0_tmp - 1;
        for (k = 0; k < lastv; k++) {
          b_tau = ((a->data[knt].re != 0.0) || (a->data[knt].im != 0.0));
          if (b_tau) {
            beta1 = a->data[knt].re * c.re + a->data[knt].im * c.im;
            temp_im = a->data[knt].re * c.im - a->data[knt].im * c.re;
            ix = 0;
            i1 = i + 1;
            i2 = b_lastc + i;
            for (im1n_tmp = i1; im1n_tmp <= i2; im1n_tmp++) {
              xnorm = work->data[ix].re * beta1 - work->data[ix].im * temp_im;
              beta1_im = work->data[ix].re * temp_im + work->data[ix].im * beta1;
              a->data[im1n_tmp - 1].re += xnorm;
              a->data[im1n_tmp - 1].im += beta1_im;
              ix++;
            }
          }

          knt++;
          i += n;
        }
      }
    }

    im1n_tmp = (b_i + in) + 2;
    beta1 = tau->data[b_i].re;
    temp_im = -tau->data[b_i].im;
    if ((beta1 != 0.0) || (temp_im != 0.0)) {
      lastv = n_tmp + 2;
      i = iv0_tmp + n_tmp;
      do {
        exitg5 = 0;
        if (lastv > 0) {
          b_tau = ((a->data[i].re == 0.0) && (a->data[i].im == 0.0));
          if (b_tau) {
            lastv--;
            i--;
          } else {
            exitg5 = 1;
          }
        } else {
          exitg5 = 2;
        }
      } while (exitg5 == 0);

      do {
        exitg4 = 0;
        if (lastc + 1 > 0) {
          i = im1n_tmp + lastc * n;
          k = i;
          do {
            exitg2 = 0;
            if (k <= (i + lastv) - 1) {
              b_tau = ((a->data[k - 1].re != 0.0) || (a->data[k - 1].im != 0.0));
              if (b_tau) {
                exitg2 = 1;
              } else {
                k++;
              }
            } else {
              lastc--;
              exitg2 = 2;
            }
          } while (exitg2 == 0);

          if (exitg2 == 1) {
            exitg4 = 1;
          }
        } else {
          exitg4 = 1;
        }
      } while (exitg4 == 0);
    } else {
      lastv = 0;
      lastc = -1;
    }

    if (lastv > 0) {
      if (lastc + 1 != 0) {
        for (i = 0; i <= lastc; i++) {
          work->data[i].re = 0.0;
          work->data[i].im = 0.0;
        }

        i = 0;
        i1 = im1n_tmp + n * lastc;
        for (knt = im1n_tmp; n < 0 ? knt >= i1 : knt <= i1; knt += n) {
          ix = iv0_tmp - 1;
          c.re = 0.0;
          c.im = 0.0;
          i2 = (knt + lastv) - 1;
          for (k = knt; k <= i2; k++) {
            c.re += a->data[k - 1].re * a->data[ix].re + a->data[k - 1].im *
              a->data[ix].im;
            c.im += a->data[k - 1].re * a->data[ix].im - a->data[k - 1].im *
              a->data[ix].re;
            ix++;
          }

          work->data[i].re += c.re - 0.0 * c.im;
          work->data[i].im += c.im + 0.0 * c.re;
          i++;
        }
      }

      c.re = -beta1;
      c.im = -temp_im;
      if ((!(-beta1 == 0.0)) || (!(-temp_im == 0.0))) {
        i = im1n_tmp - 1;
        knt = 0;
        for (k = 0; k <= lastc; k++) {
          b_tau = ((work->data[knt].re != 0.0) || (work->data[knt].im != 0.0));
          if (b_tau) {
            beta1 = work->data[knt].re * c.re + work->data[knt].im * c.im;
            temp_im = work->data[knt].re * c.im - work->data[knt].im * c.re;
            ix = iv0_tmp;
            i1 = i + 1;
            i2 = lastv + i;
            for (im1n_tmp = i1; im1n_tmp <= i2; im1n_tmp++) {
              xnorm = a->data[ix - 1].re * beta1 - a->data[ix - 1].im * temp_im;
              beta1_im = a->data[ix - 1].re * temp_im + a->data[ix - 1].im *
                beta1;
              a->data[im1n_tmp - 1].re += xnorm;
              a->data[im1n_tmp - 1].im += beta1_im;
              ix++;
            }
          }

          knt++;
          i += n;
        }
      }
    }

    a->data[(b_i + a->size[0] * b_i) + 1] = alpha1;
  }

  emxFree_creal_T(&work);
  emxFree_creal_T(&tau);
}

/* End of code generation (xgehrd.cpp) */
