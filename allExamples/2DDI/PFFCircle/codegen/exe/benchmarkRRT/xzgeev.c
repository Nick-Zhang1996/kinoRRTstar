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
#include "DI_cost.h"
#include "benchmarkRRT_rtwutil.h"

/* Function Definitions */
void xzgeev(const creal_T A_data[], const int A_size[2], int *info, creal_T
            alpha1_data[], int alpha1_size[1], creal_T beta1_data[], int
            beta1_size[1])
{
  int At_size[2];
  int ii;
  creal_T At_data[16];
  double anrm;
  int jcol;
  boolean_T exitg1;
  double absxk;
  boolean_T ilascl;
  double anrmto;
  int nzcount;
  int ilo;
  double ctoc;
  int ihi;
  boolean_T notdone;
  int exitg3;
  double stemp_im;
  int i;
  int n;
  double cto1;
  int j;
  double a;
  int jcolp1;
  int jrow;
  int At_data_tmp;
  creal_T atmp;
  boolean_T exitg4;
  int exitg2;
  At_size[0] = A_size[0];
  At_size[1] = A_size[1];
  ii = A_size[0] * A_size[1];
  if (0 <= ii - 1) {
    memcpy(&At_data[0], &A_data[0], (unsigned int)(ii * (int)sizeof(creal_T)));
  }

  *info = 0;
  anrm = 0.0;
  jcol = 0;
  exitg1 = false;
  while ((!exitg1) && (jcol <= ii - 1)) {
    absxk = rt_hypotd_snf(At_data[jcol].re, At_data[jcol].im);
    if (rtIsNaN(absxk)) {
      anrm = rtNaN;
      exitg1 = true;
    } else {
      if (absxk > anrm) {
        anrm = absxk;
      }

      jcol++;
    }
  }

  if (rtIsInf(anrm) || rtIsNaN(anrm)) {
    alpha1_size[0] = A_size[0];
    ii = A_size[0];
    for (nzcount = 0; nzcount < ii; nzcount++) {
      alpha1_data[nzcount].re = rtNaN;
      alpha1_data[nzcount].im = 0.0;
    }

    beta1_size[0] = A_size[0];
    ii = A_size[0];
    for (nzcount = 0; nzcount < ii; nzcount++) {
      beta1_data[nzcount].re = rtNaN;
      beta1_data[nzcount].im = 0.0;
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

        ii = At_size[0] * At_size[1] - 1;
        for (nzcount = 0; nzcount <= ii; nzcount++) {
          At_data[nzcount].re *= a;
          At_data[nzcount].im *= a;
        }
      }
    }

    ilo = 1;
    ihi = A_size[0];
    if (A_size[0] <= 1) {
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
          jcol = 0;
          exitg4 = false;
          while ((!exitg4) && (jcol <= ihi - 1)) {
            if ((At_data[(ii + At_size[0] * jcol) - 1].re != 0.0) || (At_data
                 [(ii + At_size[0] * jcol) - 1].im != 0.0) || (ii == jcol + 1))
            {
              if (nzcount == 0) {
                j = jcol;
                nzcount = 1;
                jcol++;
              } else {
                nzcount = 2;
                exitg4 = true;
              }
            } else {
              jcol++;
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
          n = At_size[0];
          if (i != ihi) {
            for (jcol = 1; jcol <= n; jcol++) {
              ii = At_size[0] * (jcol - 1);
              nzcount = (i + ii) - 1;
              atmp = At_data[nzcount];
              At_data_tmp = (ihi + ii) - 1;
              At_data[nzcount] = At_data[At_data_tmp];
              At_data[At_data_tmp] = atmp;
            }
          }

          if (j + 1 != ihi) {
            for (jcol = 0; jcol < ihi; jcol++) {
              ii = jcol + At_size[0] * j;
              atmp = At_data[ii];
              At_data_tmp = jcol + At_size[0] * (ihi - 1);
              At_data[ii] = At_data[At_data_tmp];
              At_data[At_data_tmp] = atmp;
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
          jcol = ilo;
          exitg1 = false;
          while ((!exitg1) && (jcol <= ihi)) {
            nzcount = 0;
            i = ihi;
            j = jcol;
            ii = ilo;
            exitg4 = false;
            while ((!exitg4) && (ii <= ihi)) {
              if ((At_data[(ii + At_size[0] * (jcol - 1)) - 1].re != 0.0) ||
                  (At_data[(ii + At_size[0] * (jcol - 1)) - 1].im != 0.0) || (ii
                   == jcol)) {
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
              jcol++;
            }
          }

          if (!notdone) {
            exitg2 = 1;
          } else {
            n = At_size[0];
            if (i != ilo) {
              for (jcol = ilo; jcol <= n; jcol++) {
                ii = At_size[0] * (jcol - 1);
                nzcount = (i + ii) - 1;
                atmp = At_data[nzcount];
                At_data_tmp = (ilo + ii) - 1;
                At_data[nzcount] = At_data[At_data_tmp];
                At_data[At_data_tmp] = atmp;
              }
            }

            if (j != ilo) {
              for (jcol = 0; jcol < ihi; jcol++) {
                ii = jcol + At_size[0] * (j - 1);
                atmp = At_data[ii];
                At_data_tmp = jcol + At_size[0] * (ilo - 1);
                At_data[ii] = At_data[At_data_tmp];
                At_data[At_data_tmp] = atmp;
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

    n = A_size[0];
    if ((A_size[0] > 1) && (ihi >= ilo + 2)) {
      for (jcol = ilo - 1; jcol + 1 < ihi - 1; jcol++) {
        jcolp1 = jcol + 2;
        for (jrow = ihi - 1; jrow + 1 > jcol + 2; jrow--) {
          At_data_tmp = jrow + At_size[0] * jcol;
          xzlartg(At_data[At_data_tmp - 1], At_data[At_data_tmp], &absxk, &atmp,
                  &At_data[(jrow + At_size[0] * jcol) - 1]);
          At_data[At_data_tmp].re = 0.0;
          At_data[At_data_tmp].im = 0.0;
          for (j = jcolp1; j <= n; j++) {
            ii = jrow + At_size[0] * (j - 1);
            nzcount = ii - 1;
            ctoc = absxk * At_data[nzcount].re + (atmp.re * At_data[ii].re -
              atmp.im * At_data[jrow + At_size[0] * (j - 1)].im);
            stemp_im = absxk * At_data[(jrow + At_size[0] * (j - 1)) - 1].im +
              (atmp.re * At_data[jrow + At_size[0] * (j - 1)].im + atmp.im *
               At_data[jrow + At_size[0] * (j - 1)].re);
            cto1 = At_data[(jrow + At_size[0] * (j - 1)) - 1].re;
            At_data[ii].re = absxk * At_data[jrow + At_size[0] * (j - 1)].re -
              (atmp.re * At_data[(jrow + At_size[0] * (j - 1)) - 1].re + atmp.im
               * At_data[(jrow + At_size[0] * (j - 1)) - 1].im);
            At_data[ii].im = absxk * At_data[ii].im - (atmp.re * At_data[(jrow +
              At_size[0] * (j - 1)) - 1].im - atmp.im * cto1);
            At_data[nzcount].re = ctoc;
            At_data[nzcount].im = stemp_im;
          }

          atmp.re = -atmp.re;
          atmp.im = -atmp.im;
          for (i = 1; i <= ihi; i++) {
            ii = (i + At_size[0] * (jrow - 1)) - 1;
            nzcount = (i + At_size[0] * jrow) - 1;
            ctoc = absxk * At_data[nzcount].re + (atmp.re * At_data[ii].re -
              atmp.im * At_data[(i + At_size[0] * (jrow - 1)) - 1].im);
            stemp_im = absxk * At_data[(i + At_size[0] * jrow) - 1].im +
              (atmp.re * At_data[(i + At_size[0] * (jrow - 1)) - 1].im + atmp.im
               * At_data[(i + At_size[0] * (jrow - 1)) - 1].re);
            cto1 = At_data[(i + At_size[0] * jrow) - 1].re;
            At_data[ii].re = absxk * At_data[(i + At_size[0] * (jrow - 1)) - 1].
              re - (atmp.re * At_data[(i + At_size[0] * jrow) - 1].re + atmp.im *
                    At_data[(i + At_size[0] * jrow) - 1].im);
            At_data[ii].im = absxk * At_data[ii].im - (atmp.re * At_data[(i +
              At_size[0] * jrow) - 1].im - atmp.im * cto1);
            At_data[nzcount].re = ctoc;
            At_data[nzcount].im = stemp_im;
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

        ii = alpha1_size[0];
        for (nzcount = 0; nzcount < ii; nzcount++) {
          alpha1_data[nzcount].re *= a;
          alpha1_data[nzcount].im *= a;
        }
      }
    }
  }
}

/* End of code generation (xzgeev.c) */
