/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * xzhgeqz.c
 *
 * Code generation for function 'xzhgeqz'
 *
 */

/* Include files */
#include <math.h>
#include <string.h>
#include "rt_nonfinite.h"
#include "benchmarkRRT.h"
#include "xzhgeqz.h"
#include "xzlartg.h"
#include "sqrt.h"

/* Function Definitions */
void xzhgeqz(const creal_T A_data[], const int A_size[2], int ilo, int ihi, int *
             info, creal_T alpha1_data[], int alpha1_size[1], creal_T
             beta1_data[], int beta1_size[1])
{
  int A_size_idx_0;
  int jp1;
  creal_T b_A_data[64];
  int n;
  int ad22_re_tmp;
  double eshift_re;
  double eshift_im;
  creal_T ctemp;
  double anorm;
  double scale;
  double reAij;
  double sumsq;
  double b_atol;
  boolean_T firstNonZero;
  int j;
  double ascale;
  int i;
  double bscale;
  double imAij;
  boolean_T guard1 = false;
  boolean_T guard2 = false;
  double temp2;
  int ifirst;
  int istart;
  int ilast;
  int ilastm1;
  int ifrstm;
  int ilastm;
  int iiter;
  boolean_T goto60;
  boolean_T goto70;
  boolean_T goto90;
  int jiter;
  int exitg1;
  boolean_T b_guard1 = false;
  boolean_T guard3 = false;
  boolean_T exitg2;
  creal_T b_ascale;
  creal_T shift;
  double ascale_re;
  double ad22_re;
  double ad22_im;
  A_size_idx_0 = A_size[0];
  jp1 = A_size[0] * A_size[1];
  if (0 <= jp1 - 1) {
    memcpy(&b_A_data[0], &A_data[0], (unsigned int)(jp1 * (int)sizeof(creal_T)));
  }

  *info = 0;
  if ((A_size[0] == 1) && (A_size[1] == 1)) {
    ihi = 1;
  }

  n = A_size[0];
  alpha1_size[0] = A_size[0];
  if (0 <= A_size[0] - 1) {
    memset(&alpha1_data[0], 0, (unsigned int)(A_size[0] * (int)sizeof(creal_T)));
  }

  beta1_size[0] = A_size[0];
  jp1 = A_size[0];
  for (ad22_re_tmp = 0; ad22_re_tmp < jp1; ad22_re_tmp++) {
    beta1_data[ad22_re_tmp].re = 1.0;
    beta1_data[ad22_re_tmp].im = 0.0;
  }

  eshift_re = 0.0;
  eshift_im = 0.0;
  ctemp.re = 0.0;
  ctemp.im = 0.0;
  anorm = 0.0;
  if (ilo <= ihi) {
    scale = 0.0;
    sumsq = 0.0;
    firstNonZero = true;
    for (j = ilo; j <= ihi; j++) {
      ad22_re_tmp = j + 1;
      if (ihi < j + 1) {
        ad22_re_tmp = ihi;
      }

      for (i = ilo; i <= ad22_re_tmp; i++) {
        reAij = A_data[(i + A_size[0] * (j - 1)) - 1].re;
        imAij = A_data[(i + A_size[0] * (j - 1)) - 1].im;
        if (reAij != 0.0) {
          anorm = fabs(reAij);
          if (firstNonZero) {
            sumsq = 1.0;
            scale = anorm;
            firstNonZero = false;
          } else if (scale < anorm) {
            temp2 = scale / anorm;
            sumsq = 1.0 + sumsq * temp2 * temp2;
            scale = anorm;
          } else {
            temp2 = anorm / scale;
            sumsq += temp2 * temp2;
          }
        }

        if (imAij != 0.0) {
          anorm = fabs(imAij);
          if (firstNonZero) {
            sumsq = 1.0;
            scale = anorm;
            firstNonZero = false;
          } else if (scale < anorm) {
            temp2 = scale / anorm;
            sumsq = 1.0 + sumsq * temp2 * temp2;
            scale = anorm;
          } else {
            temp2 = anorm / scale;
            sumsq += temp2 * temp2;
          }
        }
      }
    }

    anorm = scale * sqrt(sumsq);
  }

  reAij = 2.2204460492503131E-16 * anorm;
  b_atol = 2.2250738585072014E-308;
  if (reAij > 2.2250738585072014E-308) {
    b_atol = reAij;
  }

  reAij = 2.2250738585072014E-308;
  if (anorm > 2.2250738585072014E-308) {
    reAij = anorm;
  }

  ascale = 1.0 / reAij;
  bscale = 1.0 / sqrt(A_size[0]);
  firstNonZero = true;
  ad22_re_tmp = ihi + 1;
  for (j = ad22_re_tmp; j <= n; j++) {
    alpha1_data[j - 1] = A_data[(j + A_size[0] * (j - 1)) - 1];
  }

  guard1 = false;
  guard2 = false;
  if (ihi >= ilo) {
    ifirst = ilo;
    istart = ilo;
    ilast = ihi - 1;
    ilastm1 = ihi - 2;
    ifrstm = ilo;
    ilastm = ihi;
    iiter = 0;
    goto60 = false;
    goto70 = false;
    goto90 = false;
    jiter = 0;
    do {
      exitg1 = 0;
      if (jiter <= 30 * ((ihi - ilo) + 1) - 1) {
        b_guard1 = false;
        if (ilast + 1 == ilo) {
          goto60 = true;
          b_guard1 = true;
        } else {
          ad22_re_tmp = ilast + A_size_idx_0 * ilastm1;
          if (fabs(b_A_data[ad22_re_tmp].re) + fabs(b_A_data[ilast +
               A_size_idx_0 * ilastm1].im) <= b_atol) {
            b_A_data[ad22_re_tmp].re = 0.0;
            b_A_data[ad22_re_tmp].im = 0.0;
            goto60 = true;
            b_guard1 = true;
          } else {
            j = ilastm1;
            guard3 = false;
            exitg2 = false;
            while ((!exitg2) && (j + 1 >= ilo)) {
              if (j + 1 == ilo) {
                guard3 = true;
                exitg2 = true;
              } else {
                ad22_re_tmp = j + A_size_idx_0 * (j - 1);
                if (fabs(b_A_data[ad22_re_tmp].re) + fabs(b_A_data[j +
                     A_size_idx_0 * (j - 1)].im) <= b_atol) {
                  b_A_data[ad22_re_tmp].re = 0.0;
                  b_A_data[ad22_re_tmp].im = 0.0;
                  guard3 = true;
                  exitg2 = true;
                } else {
                  j--;
                  guard3 = false;
                }
              }
            }

            if (guard3) {
              ifirst = j + 1;
              goto70 = true;
            }

            if (goto70) {
              b_guard1 = true;
            } else {
              jp1 = alpha1_size[0];
              for (ad22_re_tmp = 0; ad22_re_tmp < jp1; ad22_re_tmp++) {
                alpha1_data[ad22_re_tmp].re = rtNaN;
                alpha1_data[ad22_re_tmp].im = 0.0;
              }

              jp1 = beta1_size[0];
              for (ad22_re_tmp = 0; ad22_re_tmp < jp1; ad22_re_tmp++) {
                beta1_data[ad22_re_tmp].re = rtNaN;
                beta1_data[ad22_re_tmp].im = 0.0;
              }

              *info = 1;
              exitg1 = 1;
            }
          }
        }

        if (b_guard1) {
          if (goto60) {
            goto60 = false;
            alpha1_data[ilast] = b_A_data[ilast + A_size_idx_0 * ilast];
            ilast = ilastm1;
            ilastm1--;
            if (ilast + 1 < ilo) {
              firstNonZero = false;
              guard2 = true;
              exitg1 = 1;
            } else {
              iiter = 0;
              eshift_re = 0.0;
              eshift_im = 0.0;
              ilastm = ilast + 1;
              if (ifrstm > ilast + 1) {
                ifrstm = ilo;
              }

              jiter++;
            }
          } else {
            if (goto70) {
              goto70 = false;
              iiter++;
              ifrstm = ifirst;
              if (iiter - iiter / 10 * 10 != 0) {
                anorm = ascale * b_A_data[ilastm1 + A_size_idx_0 * ilastm1].re;
                reAij = ascale * b_A_data[ilastm1 + A_size_idx_0 * ilastm1].im;
                if (reAij == 0.0) {
                  shift.re = anorm / bscale;
                  shift.im = 0.0;
                } else if (anorm == 0.0) {
                  shift.re = 0.0;
                  shift.im = reAij / bscale;
                } else {
                  shift.re = anorm / bscale;
                  shift.im = reAij / bscale;
                }

                jp1 = A_size_idx_0 * ilast;
                anorm = ascale * b_A_data[ilast + jp1].re;
                reAij = ascale * b_A_data[ilast + A_size_idx_0 * ilast].im;
                if (reAij == 0.0) {
                  ad22_re = anorm / bscale;
                  ad22_im = 0.0;
                } else if (anorm == 0.0) {
                  ad22_re = 0.0;
                  ad22_im = reAij / bscale;
                } else {
                  ad22_re = anorm / bscale;
                  ad22_im = reAij / bscale;
                }

                imAij = 0.5 * (shift.re + ad22_re);
                temp2 = 0.5 * (shift.im + ad22_im);
                anorm = ascale * b_A_data[ilastm1 + jp1].re;
                reAij = ascale * b_A_data[ilastm1 + A_size_idx_0 * ilast].im;
                if (reAij == 0.0) {
                  ascale_re = anorm / bscale;
                  sumsq = 0.0;
                } else if (anorm == 0.0) {
                  ascale_re = 0.0;
                  sumsq = reAij / bscale;
                } else {
                  ascale_re = anorm / bscale;
                  sumsq = reAij / bscale;
                }

                anorm = ascale * b_A_data[ilast + A_size_idx_0 * ilastm1].re;
                reAij = ascale * b_A_data[ilast + A_size_idx_0 * ilastm1].im;
                if (reAij == 0.0) {
                  scale = anorm / bscale;
                  anorm = 0.0;
                } else if (anorm == 0.0) {
                  scale = 0.0;
                  anorm = reAij / bscale;
                } else {
                  scale = anorm / bscale;
                  anorm = reAij / bscale;
                }

                reAij = shift.re * ad22_im + shift.im * ad22_re;
                shift.re = ((imAij * imAij - temp2 * temp2) + (ascale_re * scale
                  - sumsq * anorm)) - (shift.re * ad22_re - shift.im * ad22_im);
                shift.im = ((imAij * temp2 + temp2 * imAij) + (ascale_re * anorm
                  + sumsq * scale)) - reAij;
                b_sqrt(&shift);
                if ((imAij - ad22_re) * shift.re + (temp2 - ad22_im) * shift.im <=
                    0.0) {
                  shift.re += imAij;
                  shift.im += temp2;
                } else {
                  shift.re = imAij - shift.re;
                  shift.im = temp2 - shift.im;
                }
              } else {
                anorm = ascale * b_A_data[ilast + A_size_idx_0 * ilastm1].re;
                reAij = ascale * b_A_data[ilast + A_size_idx_0 * ilastm1].im;
                if (reAij == 0.0) {
                  ascale_re = anorm / bscale;
                  sumsq = 0.0;
                } else if (anorm == 0.0) {
                  ascale_re = 0.0;
                  sumsq = reAij / bscale;
                } else {
                  ascale_re = anorm / bscale;
                  sumsq = reAij / bscale;
                }

                eshift_re += ascale_re;
                eshift_im += sumsq;
                shift.re = eshift_re;
                shift.im = eshift_im;
              }

              j = ilastm1;
              jp1 = ilastm1 + 1;
              exitg2 = false;
              while ((!exitg2) && (j + 1 > ifirst)) {
                istart = j + 1;
                n = A_size_idx_0 * j;
                ctemp.re = ascale * b_A_data[j + n].re - shift.re * bscale;
                ctemp.im = ascale * b_A_data[j + A_size_idx_0 * j].im - shift.im
                  * bscale;
                anorm = fabs(ctemp.re) + fabs(ctemp.im);
                temp2 = ascale * (fabs(b_A_data[jp1 + n].re) + fabs(b_A_data[jp1
                  + A_size_idx_0 * j].im));
                reAij = anorm;
                if (temp2 > anorm) {
                  reAij = temp2;
                }

                if ((reAij < 1.0) && (reAij != 0.0)) {
                  anorm /= reAij;
                  temp2 /= reAij;
                }

                if ((fabs(b_A_data[j + A_size_idx_0 * (j - 1)].re) + fabs
                     (b_A_data[j + A_size_idx_0 * (j - 1)].im)) * temp2 <= anorm
                    * b_atol) {
                  goto90 = true;
                  exitg2 = true;
                } else {
                  jp1 = j;
                  j--;
                }
              }

              if (!goto90) {
                istart = ifirst;
                ctemp.re = ascale * b_A_data[(ifirst + A_size_idx_0 * (ifirst -
                  1)) - 1].re - shift.re * bscale;
                ctemp.im = ascale * b_A_data[(ifirst + A_size_idx_0 * (ifirst -
                  1)) - 1].im - shift.im * bscale;
                goto90 = true;
              }
            }

            if (goto90) {
              goto90 = false;
              b_ascale.re = ascale * b_A_data[istart + A_size_idx_0 * (istart -
                1)].re;
              b_ascale.im = ascale * b_A_data[istart + A_size_idx_0 * (istart -
                1)].im;
              b_xzlartg(ctemp, b_ascale, &anorm, &shift);
              j = istart;
              n = istart - 2;
              while (j < ilast + 1) {
                if (j > istart) {
                  jp1 = j + A_size_idx_0 * n;
                  xzlartg(b_A_data[jp1 - 1], b_A_data[jp1], &anorm, &shift,
                          &b_A_data[(j + A_size_idx_0 * n) - 1]);
                  b_A_data[jp1].re = 0.0;
                  b_A_data[jp1].im = 0.0;
                }

                for (n = j; n <= ilastm; n++) {
                  jp1 = j + A_size_idx_0 * (n - 1);
                  ad22_re_tmp = jp1 - 1;
                  ad22_re = anorm * b_A_data[ad22_re_tmp].re + (shift.re *
                    b_A_data[jp1].re - shift.im * b_A_data[j + A_size_idx_0 * (n
                    - 1)].im);
                  ad22_im = anorm * b_A_data[(j + A_size_idx_0 * (n - 1)) - 1].
                    im + (shift.re * b_A_data[j + A_size_idx_0 * (n - 1)].im +
                          shift.im * b_A_data[j + A_size_idx_0 * (n - 1)].re);
                  reAij = b_A_data[(j + A_size_idx_0 * (n - 1)) - 1].re;
                  b_A_data[jp1].re = anorm * b_A_data[j + A_size_idx_0 * (n - 1)]
                    .re - (shift.re * b_A_data[(j + A_size_idx_0 * (n - 1)) - 1]
                           .re + shift.im * b_A_data[(j + A_size_idx_0 * (n - 1))
                           - 1].im);
                  b_A_data[jp1].im = anorm * b_A_data[jp1].im - (shift.re *
                    b_A_data[(j + A_size_idx_0 * (n - 1)) - 1].im - shift.im *
                    reAij);
                  b_A_data[ad22_re_tmp].re = ad22_re;
                  b_A_data[ad22_re_tmp].im = ad22_im;
                }

                shift.re = -shift.re;
                shift.im = -shift.im;
                n = j;
                if (ilast + 1 < j + 2) {
                  n = ilast - 1;
                }

                for (i = ifrstm; i <= n + 2; i++) {
                  jp1 = (i + A_size_idx_0 * (j - 1)) - 1;
                  ad22_re_tmp = (i + A_size_idx_0 * j) - 1;
                  ad22_re = anorm * b_A_data[ad22_re_tmp].re + (shift.re *
                    b_A_data[jp1].re - shift.im * b_A_data[(i + A_size_idx_0 *
                    (j - 1)) - 1].im);
                  ad22_im = anorm * b_A_data[(i + A_size_idx_0 * j) - 1].im +
                    (shift.re * b_A_data[(i + A_size_idx_0 * (j - 1)) - 1].im +
                     shift.im * b_A_data[(i + A_size_idx_0 * (j - 1)) - 1].re);
                  reAij = b_A_data[(i + A_size_idx_0 * j) - 1].re;
                  b_A_data[jp1].re = anorm * b_A_data[(i + A_size_idx_0 * (j - 1))
                    - 1].re - (shift.re * b_A_data[(i + A_size_idx_0 * j) - 1].
                               re + shift.im * b_A_data[(i + A_size_idx_0 * j) -
                               1].im);
                  b_A_data[jp1].im = anorm * b_A_data[jp1].im - (shift.re *
                    b_A_data[(i + A_size_idx_0 * j) - 1].im - shift.im * reAij);
                  b_A_data[ad22_re_tmp].re = ad22_re;
                  b_A_data[ad22_re_tmp].im = ad22_im;
                }

                n = j - 1;
                j++;
              }
            }

            jiter++;
          }
        }
      } else {
        guard2 = true;
        exitg1 = 1;
      }
    } while (exitg1 == 0);
  } else {
    guard1 = true;
  }

  if (guard2) {
    if (firstNonZero) {
      *info = ilast + 1;
      for (jp1 = 0; jp1 <= ilast; jp1++) {
        alpha1_data[jp1].re = rtNaN;
        alpha1_data[jp1].im = 0.0;
        beta1_data[jp1].re = rtNaN;
        beta1_data[jp1].im = 0.0;
      }
    } else {
      guard1 = true;
    }
  }

  if (guard1) {
    for (j = 0; j <= ilo - 2; j++) {
      alpha1_data[j] = b_A_data[j + A_size_idx_0 * j];
    }
  }
}

/* End of code generation (xzhgeqz.c) */
