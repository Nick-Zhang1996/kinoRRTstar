/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * xzgeev.c
 *
 * Code generation for function 'xzgeev'
 *
 */

/* Include files */
#include <string.h>
#include "rt_nonfinite.h"
#include "benchmarkRRT.h"
#include "xzgeev.h"
#include "xzlartg.h"
#include "xzhgeqz.h"
#include "xzggbal.h"
#include "quadf_cost.h"
#include "benchmarkRRT_rtwutil.h"

/* Function Definitions */
void xzgeev(const creal_T A_data[], const int A_size[2], int *info, creal_T
            alpha1_data[], int alpha1_size[1], creal_T beta1_data[], int
            beta1_size[1])
{
  int At_size[2];
  int loop_ub_tmp;
  creal_T At_data[64];
  double anrm;
  int k;
  boolean_T exitg1;
  double absxk;
  boolean_T ilascl;
  double anrmto;
  int ilo;
  int ihi;
  int rscale_data[8];
  int rscale_size[1];
  double ctoc;
  int n;
  boolean_T notdone;
  int jcol;
  double stemp_im;
  double cto1;
  int jcolp1;
  int jrow;
  double a;
  creal_T s;
  int stemp_re_tmp;
  At_size[0] = A_size[0];
  At_size[1] = A_size[1];
  loop_ub_tmp = A_size[0] * A_size[1];
  if (0 <= loop_ub_tmp - 1) {
    memcpy(&At_data[0], &A_data[0], (unsigned int)(loop_ub_tmp * (int)sizeof
            (creal_T)));
  }

  *info = 0;
  anrm = 0.0;
  k = 0;
  exitg1 = false;
  while ((!exitg1) && (k <= loop_ub_tmp - 1)) {
    absxk = rt_hypotd_snf(A_data[k].re, A_data[k].im);
    if (rtIsNaN(absxk)) {
      anrm = rtNaN;
      exitg1 = true;
    } else {
      if (absxk > anrm) {
        anrm = absxk;
      }

      k++;
    }
  }

  if (rtIsInf(anrm) || rtIsNaN(anrm)) {
    alpha1_size[0] = A_size[0];
    loop_ub_tmp = A_size[0];
    for (k = 0; k < loop_ub_tmp; k++) {
      alpha1_data[k].re = rtNaN;
      alpha1_data[k].im = 0.0;
    }

    beta1_size[0] = A_size[0];
    loop_ub_tmp = A_size[0];
    for (k = 0; k < loop_ub_tmp; k++) {
      beta1_data[k].re = rtNaN;
      beta1_data[k].im = 0.0;
    }
  } else {
    ilascl = false;
    anrmto = anrm;
    if ((anrm > 0.0) && (anrm < 6.7178761075670888E-139)) {
      anrmto = 6.7178761075670888E-139;
      ilascl = true;
    } else {
      if (anrm > 1.4885657073574029E+138) {
        anrmto = 1.4885657073574029E+138;
        ilascl = true;
      }
    }

    if (ilascl) {
      absxk = anrm;
      ctoc = anrmto;
      notdone = true;
      while (notdone) {
        stemp_im = absxk * 2.0041683600089728E-292;
        cto1 = ctoc / 4.9896007738368E+291;
        if ((stemp_im > ctoc) && (ctoc != 0.0)) {
          a = 2.0041683600089728E-292;
          absxk = stemp_im;
        } else if (cto1 > absxk) {
          a = 4.9896007738368E+291;
          ctoc = cto1;
        } else {
          a = ctoc / absxk;
          notdone = false;
        }

        loop_ub_tmp = At_size[0] * At_size[1] - 1;
        for (k = 0; k <= loop_ub_tmp; k++) {
          At_data[k].re *= a;
          At_data[k].im *= a;
        }
      }
    }

    xzggbal(At_data, At_size, &ilo, &ihi, rscale_data, rscale_size);
    n = At_size[0];
    if ((At_size[0] > 1) && (ihi >= ilo + 2)) {
      for (jcol = ilo - 1; jcol + 1 < ihi - 1; jcol++) {
        jcolp1 = jcol + 2;
        for (jrow = ihi - 1; jrow + 1 > jcol + 2; jrow--) {
          loop_ub_tmp = jrow + At_size[0] * jcol;
          xzlartg(At_data[loop_ub_tmp - 1], At_data[loop_ub_tmp], &absxk, &s,
                  &At_data[(jrow + At_size[0] * jcol) - 1]);
          At_data[loop_ub_tmp].re = 0.0;
          At_data[loop_ub_tmp].im = 0.0;
          for (loop_ub_tmp = jcolp1; loop_ub_tmp <= n; loop_ub_tmp++) {
            k = jrow + At_size[0] * (loop_ub_tmp - 1);
            stemp_re_tmp = k - 1;
            ctoc = absxk * At_data[stemp_re_tmp].re + (s.re * At_data[k].re -
              s.im * At_data[jrow + At_size[0] * (loop_ub_tmp - 1)].im);
            stemp_im = absxk * At_data[(jrow + At_size[0] * (loop_ub_tmp - 1)) -
              1].im + (s.re * At_data[jrow + At_size[0] * (loop_ub_tmp - 1)].im
                       + s.im * At_data[jrow + At_size[0] * (loop_ub_tmp - 1)].
                       re);
            cto1 = At_data[(jrow + At_size[0] * (loop_ub_tmp - 1)) - 1].re;
            At_data[k].re = absxk * At_data[jrow + At_size[0] * (loop_ub_tmp - 1)]
              .re - (s.re * At_data[(jrow + At_size[0] * (loop_ub_tmp - 1)) - 1]
                     .re + s.im * At_data[(jrow + At_size[0] * (loop_ub_tmp - 1))
                     - 1].im);
            At_data[k].im = absxk * At_data[k].im - (s.re * At_data[(jrow +
              At_size[0] * (loop_ub_tmp - 1)) - 1].im - s.im * cto1);
            At_data[stemp_re_tmp].re = ctoc;
            At_data[stemp_re_tmp].im = stemp_im;
          }

          s.re = -s.re;
          s.im = -s.im;
          for (loop_ub_tmp = 1; loop_ub_tmp <= ihi; loop_ub_tmp++) {
            k = (loop_ub_tmp + At_size[0] * (jrow - 1)) - 1;
            stemp_re_tmp = (loop_ub_tmp + At_size[0] * jrow) - 1;
            ctoc = absxk * At_data[stemp_re_tmp].re + (s.re * At_data[k].re -
              s.im * At_data[(loop_ub_tmp + At_size[0] * (jrow - 1)) - 1].im);
            stemp_im = absxk * At_data[(loop_ub_tmp + At_size[0] * jrow) - 1].im
              + (s.re * At_data[(loop_ub_tmp + At_size[0] * (jrow - 1)) - 1].im
                 + s.im * At_data[(loop_ub_tmp + At_size[0] * (jrow - 1)) - 1].
                 re);
            cto1 = At_data[(loop_ub_tmp + At_size[0] * jrow) - 1].re;
            At_data[k].re = absxk * At_data[(loop_ub_tmp + At_size[0] * (jrow -
              1)) - 1].re - (s.re * At_data[(loop_ub_tmp + At_size[0] * jrow) -
                             1].re + s.im * At_data[(loop_ub_tmp + At_size[0] *
              jrow) - 1].im);
            At_data[k].im = absxk * At_data[k].im - (s.re * At_data[(loop_ub_tmp
              + At_size[0] * jrow) - 1].im - s.im * cto1);
            At_data[stemp_re_tmp].re = ctoc;
            At_data[stemp_re_tmp].im = stemp_im;
          }
        }
      }
    }

    xzhgeqz(At_data, At_size, ilo, ihi, info, alpha1_data, alpha1_size,
            beta1_data, beta1_size);
    if ((*info == 0) && ilascl) {
      notdone = true;
      while (notdone) {
        stemp_im = anrmto * 2.0041683600089728E-292;
        cto1 = anrm / 4.9896007738368E+291;
        if ((stemp_im > anrm) && (anrm != 0.0)) {
          a = 2.0041683600089728E-292;
          anrmto = stemp_im;
        } else if (cto1 > anrmto) {
          a = 4.9896007738368E+291;
          anrm = cto1;
        } else {
          a = anrm / anrmto;
          notdone = false;
        }

        loop_ub_tmp = alpha1_size[0];
        for (k = 0; k < loop_ub_tmp; k++) {
          alpha1_data[k].re *= a;
          alpha1_data[k].im *= a;
        }
      }
    }
  }
}

/* End of code generation (xzgeev.c) */
