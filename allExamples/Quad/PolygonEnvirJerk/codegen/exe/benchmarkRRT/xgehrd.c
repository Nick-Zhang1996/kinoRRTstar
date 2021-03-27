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
#include "xzlarf.h"
#include "recip.h"
#include "xdlapy3.h"
#include "xnrm2.h"

/* Function Definitions */
void xgehrd(creal_T a_data[], int a_size[2])
{
  int n;
  int loop_ub;
  int i5;
  creal_T work_data[8];
  int i;
  int im1n;
  int in;
  int alpha1_tmp;
  creal_T alpha1;
  int u0;
  creal_T tau_data[7];
  double xnorm_tmp;
  double beta1;
  creal_T b_a_data;
  int knt;
  double ai;
  int i6;
  int k;
  n = a_size[0];
  loop_ub = (signed char)a_size[0];
  for (i5 = 0; i5 < loop_ub; i5++) {
    work_data[i5].re = 0.0;
    work_data[i5].im = 0.0;
  }

  i5 = a_size[0];
  for (i = 0; i <= i5 - 2; i++) {
    loop_ub = i * n;
    im1n = loop_ub + 2;
    in = (i + 1) * n;
    alpha1_tmp = (i + a_size[0] * i) + 1;
    alpha1 = a_data[alpha1_tmp];
    u0 = i + 3;
    if (u0 >= n) {
      u0 = n;
    }

    loop_ub += u0;
    u0 = (n - i) - 3;
    tau_data[i].re = 0.0;
    tau_data[i].im = 0.0;
    if (u0 + 2 > 0) {
      xnorm_tmp = xnrm2(u0 + 1, a_data, loop_ub);
      if ((xnorm_tmp != 0.0) || (a_data[(i + a_size[0] * i) + 1].im != 0.0)) {
        beta1 = xdlapy3(a_data[(i + a_size[0] * i) + 1].re, a_data[(i + a_size[0]
          * i) + 1].im, xnorm_tmp);
        if (a_data[(i + a_size[0] * i) + 1].re >= 0.0) {
          beta1 = -beta1;
        }

        if (fabs(beta1) < 1.0020841800044864E-292) {
          knt = -1;
          i6 = loop_ub + u0;
          do {
            knt++;
            for (k = loop_ub; k <= i6; k++) {
              xnorm_tmp = a_data[k - 1].re;
              ai = a_data[k - 1].im;
              a_data[k - 1].re = 9.9792015476736E+291 * xnorm_tmp - 0.0 * ai;
              a_data[k - 1].im = 9.9792015476736E+291 * ai + 0.0 * xnorm_tmp;
            }

            beta1 *= 9.9792015476736E+291;
            alpha1.re *= 9.9792015476736E+291;
            alpha1.im *= 9.9792015476736E+291;
          } while (!(fabs(beta1) >= 1.0020841800044864E-292));

          beta1 = xdlapy3(alpha1.re, alpha1.im, xnrm2(u0 + 1, a_data, loop_ub));
          if (alpha1.re >= 0.0) {
            beta1 = -beta1;
          }

          xnorm_tmp = beta1 - alpha1.re;
          if (0.0 - alpha1.im == 0.0) {
            tau_data[i].re = xnorm_tmp / beta1;
            tau_data[i].im = 0.0;
          } else if (xnorm_tmp == 0.0) {
            tau_data[i].re = 0.0;
            tau_data[i].im = (0.0 - alpha1.im) / beta1;
          } else {
            tau_data[i].re = xnorm_tmp / beta1;
            tau_data[i].im = (0.0 - alpha1.im) / beta1;
          }

          b_a_data.re = alpha1.re - beta1;
          b_a_data.im = alpha1.im;
          alpha1 = recip(b_a_data);
          i6 = loop_ub + u0;
          for (k = loop_ub; k <= i6; k++) {
            xnorm_tmp = a_data[k - 1].re;
            ai = a_data[k - 1].im;
            a_data[k - 1].re = alpha1.re * xnorm_tmp - alpha1.im * ai;
            a_data[k - 1].im = alpha1.re * ai + alpha1.im * xnorm_tmp;
          }

          for (k = 0; k <= knt; k++) {
            beta1 *= 1.0020841800044864E-292;
          }

          alpha1.re = beta1;
          alpha1.im = 0.0;
        } else {
          xnorm_tmp = beta1 - a_data[(i + a_size[0] * i) + 1].re;
          ai = 0.0 - a_data[(i + a_size[0] * i) + 1].im;
          if (ai == 0.0) {
            tau_data[i].re = xnorm_tmp / beta1;
            tau_data[i].im = 0.0;
          } else if (xnorm_tmp == 0.0) {
            tau_data[i].re = 0.0;
            tau_data[i].im = ai / beta1;
          } else {
            tau_data[i].re = xnorm_tmp / beta1;
            tau_data[i].im = ai / beta1;
          }

          b_a_data.re = a_data[(i + a_size[0] * i) + 1].re - beta1;
          b_a_data.im = a_data[(i + a_size[0] * i) + 1].im;
          alpha1 = recip(b_a_data);
          i6 = loop_ub + u0;
          for (k = loop_ub; k <= i6; k++) {
            xnorm_tmp = a_data[k - 1].re;
            ai = a_data[k - 1].im;
            a_data[k - 1].re = alpha1.re * xnorm_tmp - alpha1.im * ai;
            a_data[k - 1].im = alpha1.re * ai + alpha1.im * xnorm_tmp;
          }

          alpha1.re = beta1;
          alpha1.im = 0.0;
        }
      }
    }

    a_data[alpha1_tmp].re = 1.0;
    a_data[alpha1_tmp].im = 0.0;
    xzlarf(n, (n - i) - 1, i + im1n, tau_data[i], a_data, in + 1, n, work_data);
    b_a_data.re = tau_data[i].re;
    b_a_data.im = -tau_data[i].im;
    b_xzlarf((n - i) - 1, (n - i) - 1, i + im1n, b_a_data, a_data, (i + in) + 2,
             n, work_data);
    a_data[alpha1_tmp] = alpha1;
  }
}

/* End of code generation (xgehrd.c) */
