/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * xzgeev.cpp
 *
 * Code generation for function 'xzgeev'
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "matlab_roots.h"
#include "xzgeev.h"
#include "matlab_roots_emxutil.h"
#include "xzlartg.h"
#include "xzhgeqz.h"
#include "xzhseqr.h"
#include "matlab_roots_rtwutil.h"

/* Function Definitions */
void xzgeev(const emxArray_creal_T *A, int *info, emxArray_creal_T *alpha1,
            emxArray_creal_T *beta1)
{
  emxArray_creal_T *At;
  int nzcount;
  int ii;
  double anrm;
  boolean_T exitg1;
  double absxk;
  boolean_T ilascl;
  double anrmto;
  int ilo;
  double ctoc;
  int ihi;
  boolean_T notdone;
  int exitg3;
  double cfrom1;
  int i;
  int n;
  double cto1;
  int j;
  double At_re;
  int jrow;
  creal_T b_At;
  creal_T c_At;
  boolean_T exitg4;
  double c;
  creal_T atmp;
  int exitg2;
  boolean_T d_At;
  double stemp_re;
  emxInit_creal_T(&At, 2);
  nzcount = At->size[0] * At->size[1];
  At->size[0] = A->size[0];
  At->size[1] = A->size[1];
  emxEnsureCapacity_creal_T(At, nzcount);
  ii = A->size[0] * A->size[1];
  for (nzcount = 0; nzcount < ii; nzcount++) {
    At->data[nzcount] = A->data[nzcount];
  }

  *info = 0;
  anrm = 0.0;
  nzcount = At->size[0] * At->size[1];
  ii = 0;
  exitg1 = false;
  while ((!exitg1) && (ii <= nzcount - 1)) {
    absxk = rt_hypotd_snf(At->data[ii].re, At->data[ii].im);
    if (rtIsNaN(absxk)) {
      anrm = rtNaN;
      exitg1 = true;
    } else {
      if (absxk > anrm) {
        anrm = absxk;
      }

      ii++;
    }
  }

  if (rtIsInf(anrm) || rtIsNaN(anrm)) {
    nzcount = alpha1->size[0];
    alpha1->size[0] = At->size[0];
    emxEnsureCapacity_creal_T(alpha1, nzcount);
    ii = At->size[0];
    for (nzcount = 0; nzcount < ii; nzcount++) {
      alpha1->data[nzcount].re = rtNaN;
      alpha1->data[nzcount].im = 0.0;
    }

    nzcount = beta1->size[0];
    beta1->size[0] = At->size[0];
    emxEnsureCapacity_creal_T(beta1, nzcount);
    ii = At->size[0];
    for (nzcount = 0; nzcount < ii; nzcount++) {
      beta1->data[nzcount].re = rtNaN;
      beta1->data[nzcount].im = 0.0;
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
        cfrom1 = absxk * 2.0041683600089728E-292;
        cto1 = ctoc / 4.9896007738368E+291;
        if ((cfrom1 > ctoc) && (ctoc != 0.0)) {
          At_re = 2.0041683600089728E-292;
          absxk = cfrom1;
        } else if (cto1 > absxk) {
          At_re = 4.9896007738368E+291;
          ctoc = cto1;
        } else {
          At_re = ctoc / absxk;
          notdone = false;
        }

        nzcount = At->size[0] * At->size[1];
        ii = At->size[0] * At->size[1];
        emxEnsureCapacity_creal_T(At, ii);
        ii = nzcount - 1;
        for (nzcount = 0; nzcount <= ii; nzcount++) {
          At->data[nzcount].re *= At_re;
          At->data[nzcount].im *= At_re;
        }
      }
    }

    ilo = 1;
    ihi = At->size[0];
    if (At->size[0] <= 1) {
      ihi = 1;
    } else {
      do {
        exitg3 = 0;
        i = 0;
        j = -1;
        notdone = false;
        ii = ihi;
        exitg1 = false;
        while ((!exitg1) && (ii > 0)) {
          nzcount = 0;
          i = ii;
          j = ihi - 1;
          jrow = 0;
          exitg4 = false;
          while ((!exitg4) && (jrow <= ihi - 1)) {
            d_At = ((At->data[(ii + At->size[0] * jrow) - 1].re != 0.0) ||
                    (At->data[(ii + At->size[0] * jrow) - 1].im != 0.0));
            if (d_At || (ii == jrow + 1)) {
              if (nzcount == 0) {
                j = jrow;
                nzcount = 1;
                jrow++;
              } else {
                nzcount = 2;
                exitg4 = true;
              }
            } else {
              jrow++;
            }
          }

          if (nzcount < 2) {
            notdone = true;
            exitg1 = true;
          } else {
            ii--;
          }
        }

        if (!notdone) {
          exitg3 = 2;
        } else {
          n = At->size[0];
          if (i != ihi) {
            for (ii = 1; ii <= n; ii++) {
              atmp = At->data[(i + At->size[0] * (ii - 1)) - 1];
              At->data[(i + At->size[0] * (ii - 1)) - 1] = At->data[(ihi +
                At->size[0] * (ii - 1)) - 1];
              At->data[(ihi + At->size[0] * (ii - 1)) - 1] = atmp;
            }
          }

          if (j + 1 != ihi) {
            for (ii = 0; ii < ihi; ii++) {
              atmp = At->data[ii + At->size[0] * j];
              At->data[ii + At->size[0] * j] = At->data[ii + At->size[0] * (ihi
                - 1)];
              At->data[ii + At->size[0] * (ihi - 1)] = atmp;
            }
          }

          ihi--;
          if (ihi == 1) {
            exitg3 = 1;
          }
        }
      } while (exitg3 == 0);

      if (exitg3 == 1) {
      } else {
        do {
          exitg2 = 0;
          i = 0;
          j = 0;
          notdone = false;
          jrow = ilo;
          exitg1 = false;
          while ((!exitg1) && (jrow <= ihi)) {
            nzcount = 0;
            i = ihi;
            j = jrow;
            ii = ilo;
            exitg4 = false;
            while ((!exitg4) && (ii <= ihi)) {
              d_At = ((At->data[(ii + At->size[0] * (jrow - 1)) - 1].re != 0.0) ||
                      (At->data[(ii + At->size[0] * (jrow - 1)) - 1].im != 0.0));
              if (d_At || (ii == jrow)) {
                if (nzcount == 0) {
                  i = ii;
                  nzcount = 1;
                  ii++;
                } else {
                  nzcount = 2;
                  exitg4 = true;
                }
              } else {
                ii++;
              }
            }

            if (nzcount < 2) {
              notdone = true;
              exitg1 = true;
            } else {
              jrow++;
            }
          }

          if (!notdone) {
            exitg2 = 1;
          } else {
            n = At->size[0];
            if (i != ilo) {
              for (ii = ilo; ii <= n; ii++) {
                atmp = At->data[(i + At->size[0] * (ii - 1)) - 1];
                At->data[(i + At->size[0] * (ii - 1)) - 1] = At->data[(ilo +
                  At->size[0] * (ii - 1)) - 1];
                At->data[(ilo + At->size[0] * (ii - 1)) - 1] = atmp;
              }
            }

            if (j != ilo) {
              for (ii = 0; ii < ihi; ii++) {
                atmp = At->data[ii + At->size[0] * (j - 1)];
                At->data[ii + At->size[0] * (j - 1)] = At->data[ii + At->size[0]
                  * (ilo - 1)];
                At->data[ii + At->size[0] * (ilo - 1)] = atmp;
              }
            }

            ilo++;
            if (ilo == ihi) {
              exitg2 = 1;
            }
          }
        } while (exitg2 == 0);
      }
    }

    n = At->size[0];
    if ((At->size[0] > 1) && (ihi >= ilo + 2)) {
      for (ii = ilo - 1; ii + 1 < ihi - 1; ii++) {
        nzcount = ii + 2;
        for (jrow = ihi - 1; jrow + 1 > ii + 2; jrow--) {
          b_At = At->data[(jrow + At->size[0] * ii) - 1];
          c_At = At->data[jrow + At->size[0] * ii];
          xzlartg(b_At, c_At, &c, &atmp, &At->data[(jrow + At->size[0] * ii) - 1]);
          At->data[jrow + At->size[0] * ii].re = 0.0;
          At->data[jrow + At->size[0] * ii].im = 0.0;
          for (j = nzcount; j <= n; j++) {
            absxk = atmp.re * At->data[jrow + At->size[0] * (j - 1)].re -
              atmp.im * At->data[jrow + At->size[0] * (j - 1)].im;
            ctoc = atmp.re * At->data[jrow + At->size[0] * (j - 1)].im + atmp.im
              * At->data[jrow + At->size[0] * (j - 1)].re;
            stemp_re = c * At->data[(jrow + At->size[0] * (j - 1)) - 1].re +
              absxk;
            absxk = c * At->data[(jrow + At->size[0] * (j - 1)) - 1].im + ctoc;
            ctoc = At->data[(jrow + At->size[0] * (j - 1)) - 1].re;
            cfrom1 = At->data[(jrow + At->size[0] * (j - 1)) - 1].im;
            cto1 = At->data[(jrow + At->size[0] * (j - 1)) - 1].im;
            At_re = At->data[(jrow + At->size[0] * (j - 1)) - 1].re;
            At->data[jrow + At->size[0] * (j - 1)].re = c * At->data[jrow +
              At->size[0] * (j - 1)].re - (atmp.re * ctoc + atmp.im * cfrom1);
            At->data[jrow + At->size[0] * (j - 1)].im = c * At->data[jrow +
              At->size[0] * (j - 1)].im - (atmp.re * cto1 - atmp.im * At_re);
            At->data[(jrow + At->size[0] * (j - 1)) - 1].re = stemp_re;
            At->data[(jrow + At->size[0] * (j - 1)) - 1].im = absxk;
          }

          atmp.re = -atmp.re;
          atmp.im = -atmp.im;
          for (i = 1; i <= ihi; i++) {
            absxk = atmp.re * At->data[(i + At->size[0] * (jrow - 1)) - 1].re -
              atmp.im * At->data[(i + At->size[0] * (jrow - 1)) - 1].im;
            ctoc = atmp.re * At->data[(i + At->size[0] * (jrow - 1)) - 1].im +
              atmp.im * At->data[(i + At->size[0] * (jrow - 1)) - 1].re;
            stemp_re = c * At->data[(i + At->size[0] * jrow) - 1].re + absxk;
            absxk = c * At->data[(i + At->size[0] * jrow) - 1].im + ctoc;
            ctoc = At->data[(i + At->size[0] * jrow) - 1].re;
            cfrom1 = At->data[(i + At->size[0] * jrow) - 1].im;
            cto1 = At->data[(i + At->size[0] * jrow) - 1].im;
            At_re = At->data[(i + At->size[0] * jrow) - 1].re;
            At->data[(i + At->size[0] * (jrow - 1)) - 1].re = c * At->data[(i +
              At->size[0] * (jrow - 1)) - 1].re - (atmp.re * ctoc + atmp.im *
              cfrom1);
            At->data[(i + At->size[0] * (jrow - 1)) - 1].im = c * At->data[(i +
              At->size[0] * (jrow - 1)) - 1].im - (atmp.re * cto1 - atmp.im *
              At_re);
            At->data[(i + At->size[0] * jrow) - 1].re = stemp_re;
            At->data[(i + At->size[0] * jrow) - 1].im = absxk;
          }
        }
      }
    }

    xzhgeqz(At, ilo, ihi, info, alpha1, beta1);
    if ((*info == 0) && ilascl) {
      notdone = true;
      while (notdone) {
        cfrom1 = anrmto * 2.0041683600089728E-292;
        cto1 = anrm / 4.9896007738368E+291;
        if ((cfrom1 > anrm) && (anrm != 0.0)) {
          At_re = 2.0041683600089728E-292;
          anrmto = cfrom1;
        } else if (cto1 > anrmto) {
          At_re = 4.9896007738368E+291;
          anrm = cto1;
        } else {
          At_re = anrm / anrmto;
          notdone = false;
        }

        nzcount = alpha1->size[0];
        emxEnsureCapacity_creal_T(alpha1, nzcount);
        ii = alpha1->size[0];
        for (nzcount = 0; nzcount < ii; nzcount++) {
          alpha1->data[nzcount].re *= At_re;
          alpha1->data[nzcount].im *= At_re;
        }
      }
    }
  }

  emxFree_creal_T(&At);
}

/* End of code generation (xzgeev.cpp) */
