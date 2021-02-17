/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * xgehrd.c
 *
 * Code generation for function 'xgehrd'
 *
 */

/* Include files */
#include <math.h>
#include "rt_nonfinite.h"
#include "benchmarkRRT.h"
#include "xgehrd.h"
#include "recip.h"
#include "xdlapy3.h"
#include "xnrm2.h"

/* Function Definitions */
void xgehrd(creal_T a_data[], int a_size[2])
{
  int n;
  int i;
  int i5;
  creal_T work_data[4];
  int b_i;
  int im1n_tmp;
  int in;
  int alpha1_tmp;
  creal_T alpha1;
  int n_tmp_tmp;
  int n_tmp;
  creal_T tau_data[3];
  double c_re;
  double beta1;
  int iv0_tmp;
  int lastv;
  int knt;
  int lastc;
  double c_im;
  int i6;
  int jy;
  boolean_T exitg2;
  int ix;
  double a_data_im;
  int exitg1;
  int i7;
  double temp_im;
  creal_T b_alpha1;
  n = a_size[0];
  i = (signed char)a_size[0];
  for (i5 = 0; i5 < i; i5++) {
    work_data[i5].re = 0.0;
    work_data[i5].im = 0.0;
  }

  i5 = a_size[0];
  for (b_i = 0; b_i <= i5 - 2; b_i++) {
    im1n_tmp = b_i * n;
    in = (b_i + 1) * n;
    alpha1_tmp = (b_i + a_size[0] * b_i) + 1;
    alpha1 = a_data[alpha1_tmp];
    i = b_i + 3;
    if (i >= n) {
      i = n;
    }

    i += im1n_tmp;
    n_tmp_tmp = n - b_i;
    n_tmp = n_tmp_tmp - 3;
    tau_data[b_i].re = 0.0;
    tau_data[b_i].im = 0.0;
    if (n_tmp + 2 > 0) {
      c_re = xnrm2(n_tmp + 1, a_data, i);
      if ((c_re != 0.0) || (a_data[(b_i + a_size[0] * b_i) + 1].im != 0.0)) {
        beta1 = xdlapy3(a_data[(b_i + a_size[0] * b_i) + 1].re, a_data[(b_i +
          a_size[0] * b_i) + 1].im, c_re);
        if (a_data[(b_i + a_size[0] * b_i) + 1].re >= 0.0) {
          beta1 = -beta1;
        }

        if (fabs(beta1) < 1.0020841800044864E-292) {
          knt = -1;
          i6 = i + n_tmp;
          do {
            knt++;
            for (jy = i; jy <= i6; jy++) {
              c_re = a_data[jy - 1].re;
              a_data_im = a_data[jy - 1].im;
              a_data[jy - 1].re = 9.9792015476736E+291 * c_re - 0.0 * a_data_im;
              a_data[jy - 1].im = 9.9792015476736E+291 * a_data_im + 0.0 * c_re;
            }

            beta1 *= 9.9792015476736E+291;
            alpha1.re *= 9.9792015476736E+291;
            alpha1.im *= 9.9792015476736E+291;
          } while (!(fabs(beta1) >= 1.0020841800044864E-292));

          beta1 = xdlapy3(alpha1.re, alpha1.im, xnrm2(n_tmp + 1, a_data, i));
          if (alpha1.re >= 0.0) {
            beta1 = -beta1;
          }

          c_re = beta1 - alpha1.re;
          if (0.0 - alpha1.im == 0.0) {
            tau_data[b_i].re = c_re / beta1;
            tau_data[b_i].im = 0.0;
          } else if (c_re == 0.0) {
            tau_data[b_i].re = 0.0;
            tau_data[b_i].im = (0.0 - alpha1.im) / beta1;
          } else {
            tau_data[b_i].re = c_re / beta1;
            tau_data[b_i].im = (0.0 - alpha1.im) / beta1;
          }

          b_alpha1.re = alpha1.re - beta1;
          b_alpha1.im = alpha1.im;
          alpha1 = recip(b_alpha1);
          i6 = i + n_tmp;
          for (jy = i; jy <= i6; jy++) {
            c_re = a_data[jy - 1].re;
            a_data_im = a_data[jy - 1].im;
            a_data[jy - 1].re = alpha1.re * c_re - alpha1.im * a_data_im;
            a_data[jy - 1].im = alpha1.re * a_data_im + alpha1.im * c_re;
          }

          for (jy = 0; jy <= knt; jy++) {
            beta1 *= 1.0020841800044864E-292;
          }

          alpha1.re = beta1;
          alpha1.im = 0.0;
        } else {
          c_re = beta1 - a_data[(b_i + a_size[0] * b_i) + 1].re;
          c_im = 0.0 - a_data[(b_i + a_size[0] * b_i) + 1].im;
          if (c_im == 0.0) {
            tau_data[b_i].re = c_re / beta1;
            tau_data[b_i].im = 0.0;
          } else if (c_re == 0.0) {
            tau_data[b_i].re = 0.0;
            tau_data[b_i].im = c_im / beta1;
          } else {
            tau_data[b_i].re = c_re / beta1;
            tau_data[b_i].im = c_im / beta1;
          }

          alpha1.re = a_data[(b_i + a_size[0] * b_i) + 1].re - beta1;
          alpha1.im = a_data[(b_i + a_size[0] * b_i) + 1].im;
          alpha1 = recip(alpha1);
          i6 = i + n_tmp;
          for (jy = i; jy <= i6; jy++) {
            c_re = a_data[jy - 1].re;
            a_data_im = a_data[jy - 1].im;
            a_data[jy - 1].re = alpha1.re * c_re - alpha1.im * a_data_im;
            a_data[jy - 1].im = alpha1.re * a_data_im + alpha1.im * c_re;
          }

          alpha1.re = beta1;
          alpha1.im = 0.0;
        }
      }
    }

    a_data[alpha1_tmp].re = 1.0;
    a_data[alpha1_tmp].im = 0.0;
    iv0_tmp = (b_i + im1n_tmp) + 2;
    im1n_tmp = in + 1;
    if ((tau_data[b_i].re != 0.0) || (tau_data[b_i].im != 0.0)) {
      lastv = n_tmp + 2;
      i = iv0_tmp + n_tmp;
      while ((lastv > 0) && ((a_data[i].re == 0.0) && (a_data[i].im == 0.0))) {
        lastv--;
        i--;
      }

      lastc = n;
      exitg2 = false;
      while ((!exitg2) && (lastc > 0)) {
        i = in + lastc;
        knt = i;
        do {
          exitg1 = 0;
          if ((n > 0) && (knt <= i + (lastv - 1) * n)) {
            if ((a_data[knt - 1].re != 0.0) || (a_data[knt - 1].im != 0.0)) {
              exitg1 = 1;
            } else {
              knt += n;
            }
          } else {
            lastc--;
            exitg1 = 2;
          }
        } while (exitg1 == 0);

        if (exitg1 == 1) {
          exitg2 = true;
        }
      }
    } else {
      lastv = 0;
      lastc = 0;
    }

    if (lastv > 0) {
      if (lastc != 0) {
        for (i = 0; i < lastc; i++) {
          work_data[i].re = 0.0;
          work_data[i].im = 0.0;
        }

        ix = iv0_tmp;
        i6 = (in + n * (lastv - 1)) + 1;
        for (jy = im1n_tmp; n < 0 ? jy >= i6 : jy <= i6; jy += n) {
          c_re = a_data[ix - 1].re - 0.0 * a_data[ix - 1].im;
          c_im = a_data[ix - 1].im + 0.0 * a_data[ix - 1].re;
          i = 0;
          i7 = (jy + lastc) - 1;
          for (knt = jy; knt <= i7; knt++) {
            work_data[i].re += a_data[knt - 1].re * c_re - a_data[knt - 1].im *
              c_im;
            work_data[i].im += a_data[knt - 1].re * c_im + a_data[knt - 1].im *
              c_re;
            i++;
          }

          ix++;
        }
      }

      c_re = -tau_data[b_i].re;
      c_im = -tau_data[b_i].im;
      if ((!(-tau_data[b_i].re == 0.0)) || (!(-tau_data[b_i].im == 0.0))) {
        i = in;
        jy = iv0_tmp - 1;
        for (knt = 0; knt < lastv; knt++) {
          if ((a_data[jy].re != 0.0) || (a_data[jy].im != 0.0)) {
            beta1 = a_data[jy].re * c_re + a_data[jy].im * c_im;
            temp_im = a_data[jy].re * c_im - a_data[jy].im * c_re;
            ix = 0;
            i6 = i + 1;
            i7 = lastc + i;
            for (im1n_tmp = i6; im1n_tmp <= i7; im1n_tmp++) {
              a_data[im1n_tmp - 1].re += work_data[ix].re * beta1 - work_data[ix]
                .im * temp_im;
              a_data[im1n_tmp - 1].im += work_data[ix].re * temp_im +
                work_data[ix].im * beta1;
              ix++;
            }
          }

          jy++;
          i += n;
        }
      }
    }

    im1n_tmp = (b_i + in) + 2;
    if ((tau_data[b_i].re != 0.0) || (-tau_data[b_i].im != 0.0)) {
      lastv = n_tmp + 2;
      i = iv0_tmp + n_tmp;
      while ((lastv > 0) && ((a_data[i].re == 0.0) && (a_data[i].im == 0.0))) {
        lastv--;
        i--;
      }

      lastc = n_tmp_tmp - 2;
      exitg2 = false;
      while ((!exitg2) && (lastc + 1 > 0)) {
        i = im1n_tmp + lastc * n;
        knt = i;
        do {
          exitg1 = 0;
          if (knt <= (i + lastv) - 1) {
            if ((a_data[knt - 1].re != 0.0) || (a_data[knt - 1].im != 0.0)) {
              exitg1 = 1;
            } else {
              knt++;
            }
          } else {
            lastc--;
            exitg1 = 2;
          }
        } while (exitg1 == 0);

        if (exitg1 == 1) {
          exitg2 = true;
        }
      }
    } else {
      lastv = 0;
      lastc = -1;
    }

    if (lastv > 0) {
      if (lastc + 1 != 0) {
        for (i = 0; i <= lastc; i++) {
          work_data[i].re = 0.0;
          work_data[i].im = 0.0;
        }

        i = 0;
        i6 = im1n_tmp + n * lastc;
        for (jy = im1n_tmp; n < 0 ? jy >= i6 : jy <= i6; jy += n) {
          ix = iv0_tmp - 1;
          c_re = 0.0;
          c_im = 0.0;
          i7 = (jy + lastv) - 1;
          for (knt = jy; knt <= i7; knt++) {
            c_re += a_data[knt - 1].re * a_data[ix].re + a_data[knt - 1].im *
              a_data[ix].im;
            c_im += a_data[knt - 1].re * a_data[ix].im - a_data[knt - 1].im *
              a_data[ix].re;
            ix++;
          }

          work_data[i].re += c_re - 0.0 * c_im;
          work_data[i].im += c_im + 0.0 * c_re;
          i++;
        }
      }

      c_re = -tau_data[b_i].re;
      c_im = tau_data[b_i].im;
      if ((!(-tau_data[b_i].re == 0.0)) || (!(tau_data[b_i].im == 0.0))) {
        i = im1n_tmp - 1;
        jy = 0;
        for (knt = 0; knt <= lastc; knt++) {
          if ((work_data[jy].re != 0.0) || (work_data[jy].im != 0.0)) {
            beta1 = work_data[jy].re * c_re + work_data[jy].im * c_im;
            temp_im = work_data[jy].re * c_im - work_data[jy].im * c_re;
            ix = iv0_tmp;
            i6 = i + 1;
            i7 = lastv + i;
            for (im1n_tmp = i6; im1n_tmp <= i7; im1n_tmp++) {
              a_data_im = a_data[ix - 1].re * temp_im + a_data[ix - 1].im *
                beta1;
              a_data[im1n_tmp - 1].re += a_data[ix - 1].re * beta1 - a_data[ix -
                1].im * temp_im;
              a_data[im1n_tmp - 1].im += a_data_im;
              ix++;
            }
          }

          jy++;
          i += n;
        }
      }
    }

    a_data[alpha1_tmp] = alpha1;
  }
}

/* End of code generation (xgehrd.c) */
