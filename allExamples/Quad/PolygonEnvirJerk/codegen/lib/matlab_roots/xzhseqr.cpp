/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * xzhseqr.cpp
 *
 * Code generation for function 'xzhseqr'
 *
 */

/* Include files */
#include <cmath>
#include "rt_nonfinite.h"
#include "matlab_roots.h"
#include "xzhseqr.h"
#include "xscal.h"
#include "xzlarfg.h"
#include "sqrt.h"
#include "matlab_roots_rtwutil.h"

/* Function Definitions */
int eml_zlahqr(emxArray_creal_T *h)
{
  int info;
  int n;
  int u1;
  double itmax;
  int ldh;
  int j;
  int i;
  double SMLNUM;
  double tst;
  boolean_T exitg1;
  double aa;
  double htmp1;
  int L;
  creal_T sc;
  boolean_T goto140;
  int its;
  boolean_T exitg2;
  int k;
  boolean_T exitg3;
  creal_T x2;
  double t_re;
  double t_im;
  double ab;
  boolean_T goto70;
  double ba;
  int m;
  double u_re;
  double u_im;
  double s;
  int b_k;
  creal_T v[2];
  double b_SMLNUM;
  n = h->size[0];
  u1 = h->size[0];
  if (10 > u1) {
    u1 = 10;
  }

  itmax = 30.0 * (double)u1;
  ldh = h->size[0];
  info = 0;
  if (1 != h->size[0]) {
    u1 = h->size[0];
    for (j = 0; j <= u1 - 4; j++) {
      h->data[(j + h->size[0] * j) + 2].re = 0.0;
      h->data[(j + h->size[0] * j) + 2].im = 0.0;
      h->data[(j + h->size[0] * j) + 3].re = 0.0;
      h->data[(j + h->size[0] * j) + 3].im = 0.0;
    }

    if (1 <= n - 2) {
      h->data[(n + h->size[0] * (n - 3)) - 1].re = 0.0;
      h->data[(n + h->size[0] * (n - 3)) - 1].im = 0.0;
    }

    for (i = 2; i <= n; i++) {
      if (h->data[(i + h->size[0] * (i - 2)) - 1].im != 0.0) {
        tst = h->data[(i + h->size[0] * (i - 2)) - 1].re;
        aa = h->data[(i + h->size[0] * (i - 2)) - 1].im;
        htmp1 = std::abs(h->data[(i + h->size[0] * (i - 2)) - 1].re) + std::abs
          (h->data[(i + h->size[0] * (i - 2)) - 1].im);
        if (aa == 0.0) {
          sc.re = tst / htmp1;
          sc.im = 0.0;
        } else if (tst == 0.0) {
          sc.re = 0.0;
          sc.im = aa / htmp1;
        } else {
          sc.re = tst / htmp1;
          sc.im = aa / htmp1;
        }

        htmp1 = rt_hypotd_snf(sc.re, sc.im);
        if (-sc.im == 0.0) {
          sc.re /= htmp1;
          sc.im = 0.0;
        } else if (sc.re == 0.0) {
          sc.re = 0.0;
          sc.im = -sc.im / htmp1;
        } else {
          sc.re /= htmp1;
          sc.im = -sc.im / htmp1;
        }

        tst = h->data[(i + h->size[0] * (i - 2)) - 1].re;
        aa = h->data[(i + h->size[0] * (i - 2)) - 1].im;
        h->data[(i + h->size[0] * (i - 2)) - 1].re = rt_hypotd_snf(tst, aa);
        h->data[(i + h->size[0] * (i - 2)) - 1].im = 0.0;
        b_xscal((n - i) + 1, sc, h, i + (i - 1) * ldh, ldh);
        x2.re = sc.re;
        x2.im = -sc.im;
        u1 = i + 1;
        if (n < u1) {
          u1 = n;
        }

        xscal(u1, x2, h, 1 + (i - 1) * ldh);
      }
    }

    SMLNUM = 2.2250738585072014E-308 * ((double)n / 2.2204460492503131E-16);
    i = n - 1;
    exitg1 = false;
    while ((!exitg1) && (i + 1 >= 1)) {
      L = -1;
      goto140 = false;
      its = 0;
      exitg2 = false;
      while ((!exitg2) && (its <= (int)itmax)) {
        k = i;
        exitg3 = false;
        while ((!exitg3) && ((k + 1 > L + 2) && (!(std::abs(h->data[k + h->size
                   [0] * (k - 1)].re) + std::abs(h->data[k + h->size[0] * (k - 1)]
                   .im) <= SMLNUM)))) {
          tst = (std::abs(h->data[(k + h->size[0] * (k - 1)) - 1].re) + std::abs
                 (h->data[(k + h->size[0] * (k - 1)) - 1].im)) + (std::abs
            (h->data[k + h->size[0] * k].re) + std::abs(h->data[k + h->size[0] *
            k].im));
          if (tst == 0.0) {
            if (k - 1 >= 1) {
              tst = std::abs(h->data[(k + h->size[0] * (k - 2)) - 1].re);
            }

            if (k + 2 <= n) {
              tst += std::abs(h->data[(k + h->size[0] * k) + 1].re);
            }
          }

          if (std::abs(h->data[k + h->size[0] * (k - 1)].re) <=
              2.2204460492503131E-16 * tst) {
            htmp1 = std::abs(h->data[k + h->size[0] * (k - 1)].re) + std::abs
              (h->data[k + h->size[0] * (k - 1)].im);
            tst = std::abs(h->data[(k + h->size[0] * k) - 1].re) + std::abs
              (h->data[(k + h->size[0] * k) - 1].im);
            if (htmp1 > tst) {
              ab = htmp1;
              ba = tst;
            } else {
              ab = tst;
              ba = htmp1;
            }

            htmp1 = std::abs(h->data[k + h->size[0] * k].re) + std::abs(h->
              data[k + h->size[0] * k].im);
            t_re = h->data[(k + h->size[0] * (k - 1)) - 1].re - h->data[k +
              h->size[0] * k].re;
            t_im = h->data[(k + h->size[0] * (k - 1)) - 1].im - h->data[k +
              h->size[0] * k].im;
            tst = std::abs(t_re) + std::abs(t_im);
            if (htmp1 > tst) {
              aa = htmp1;
              htmp1 = tst;
            } else {
              aa = tst;
            }

            s = aa + ab;
            tst = 2.2204460492503131E-16 * (htmp1 * (aa / s));
            if ((SMLNUM > tst) || rtIsNaN(tst)) {
              b_SMLNUM = SMLNUM;
            } else {
              b_SMLNUM = tst;
            }

            if (ba * (ab / s) <= b_SMLNUM) {
              exitg3 = true;
            } else {
              k--;
            }
          } else {
            k--;
          }
        }

        L = k - 1;
        if (k + 1 > 1) {
          h->data[k + h->size[0] * (k - 1)].re = 0.0;
          h->data[k + h->size[0] * (k - 1)].im = 0.0;
        }

        if (k + 1 >= i + 1) {
          goto140 = true;
          exitg2 = true;
        } else {
          if (its == 10) {
            t_re = 0.75 * std::abs(h->data[(k + h->size[0] * k) + 1].re) +
              h->data[k + h->size[0] * k].re;
            t_im = h->data[k + h->size[0] * k].im;
          } else if (its == 20) {
            t_re = 0.75 * std::abs(h->data[i + h->size[0] * (i - 1)].re) +
              h->data[i + h->size[0] * i].re;
            t_im = h->data[i + h->size[0] * i].im;
          } else {
            t_re = h->data[i + h->size[0] * i].re;
            t_im = h->data[i + h->size[0] * i].im;
            x2 = h->data[(i + h->size[0] * i) - 1];
            b_sqrt(&x2);
            sc = h->data[i + h->size[0] * (i - 1)];
            b_sqrt(&sc);
            u_re = x2.re * sc.re - x2.im * sc.im;
            u_im = x2.re * sc.im + x2.im * sc.re;
            s = std::abs(u_re) + std::abs(u_im);
            if (s != 0.0) {
              tst = h->data[(i + h->size[0] * (i - 1)) - 1].re - h->data[i +
                h->size[0] * i].re;
              aa = h->data[(i + h->size[0] * (i - 1)) - 1].im - h->data[i +
                h->size[0] * i].im;
              t_re = 0.5 * tst;
              t_im = 0.5 * aa;
              tst = std::abs(t_re) + std::abs(t_im);
              if ((!(s > tst)) && (!rtIsNaN(tst))) {
                s = tst;
              }

              if (t_im == 0.0) {
                x2.re = t_re / s;
                x2.im = 0.0;
              } else if (t_re == 0.0) {
                x2.re = 0.0;
                x2.im = t_im / s;
              } else {
                x2.re = t_re / s;
                x2.im = t_im / s;
              }

              aa = x2.re;
              htmp1 = x2.re;
              x2.re = x2.re * x2.re - x2.im * x2.im;
              x2.im = aa * x2.im + x2.im * htmp1;
              if (u_im == 0.0) {
                sc.re = u_re / s;
                sc.im = 0.0;
              } else if (u_re == 0.0) {
                sc.re = 0.0;
                sc.im = u_im / s;
              } else {
                sc.re = u_re / s;
                sc.im = u_im / s;
              }

              x2.re += sc.re * sc.re - sc.im * sc.im;
              x2.im += sc.re * sc.im + sc.im * sc.re;
              b_sqrt(&x2);
              sc.re = s * x2.re;
              sc.im = s * x2.im;
              if (tst > 0.0) {
                if (t_im == 0.0) {
                  x2.re = t_re / tst;
                  x2.im = 0.0;
                } else if (t_re == 0.0) {
                  x2.re = 0.0;
                  x2.im = t_im / tst;
                } else {
                  x2.re = t_re / tst;
                  x2.im = t_im / tst;
                }

                if (x2.re * sc.re + x2.im * sc.im < 0.0) {
                  sc.re = -sc.re;
                  sc.im = -sc.im;
                }
              }

              htmp1 = t_re + sc.re;
              aa = t_im + sc.im;
              if (aa == 0.0) {
                if (u_im == 0.0) {
                  ba = u_re / htmp1;
                  tst = 0.0;
                } else if (u_re == 0.0) {
                  ba = 0.0;
                  tst = u_im / htmp1;
                } else {
                  ba = u_re / htmp1;
                  tst = u_im / htmp1;
                }
              } else if (htmp1 == 0.0) {
                if (u_re == 0.0) {
                  ba = u_im / aa;
                  tst = 0.0;
                } else if (u_im == 0.0) {
                  ba = 0.0;
                  tst = -(u_re / aa);
                } else {
                  ba = u_im / aa;
                  tst = -(u_re / aa);
                }
              } else {
                ab = std::abs(htmp1);
                tst = std::abs(aa);
                if (ab > tst) {
                  s = aa / htmp1;
                  tst = htmp1 + s * aa;
                  ba = (u_re + s * u_im) / tst;
                  tst = (u_im - s * u_re) / tst;
                } else if (tst == ab) {
                  if (htmp1 > 0.0) {
                    htmp1 = 0.5;
                  } else {
                    htmp1 = -0.5;
                  }

                  if (aa > 0.0) {
                    tst = 0.5;
                  } else {
                    tst = -0.5;
                  }

                  ba = (u_re * htmp1 + u_im * tst) / ab;
                  tst = (u_im * htmp1 - u_re * tst) / ab;
                } else {
                  s = htmp1 / aa;
                  tst = aa + s * htmp1;
                  ba = (s * u_re + u_im) / tst;
                  tst = (s * u_im - u_re) / tst;
                }
              }

              t_re = h->data[i + h->size[0] * i].re - (u_re * ba - u_im * tst);
              t_im = h->data[i + h->size[0] * i].im - (u_re * tst + u_im * ba);
            }
          }

          goto70 = false;
          m = i;
          exitg3 = false;
          while ((!exitg3) && (m > k + 1)) {
            sc.re = h->data[(m + h->size[0] * (m - 1)) - 1].re - t_re;
            sc.im = h->data[(m + h->size[0] * (m - 1)) - 1].im - t_im;
            tst = h->data[m + h->size[0] * (m - 1)].re;
            s = (std::abs(sc.re) + std::abs(sc.im)) + std::abs(tst);
            if (sc.im == 0.0) {
              sc.re /= s;
              sc.im = 0.0;
            } else if (sc.re == 0.0) {
              sc.re = 0.0;
              sc.im /= s;
            } else {
              sc.re /= s;
              sc.im /= s;
            }

            tst /= s;
            v[0] = sc;
            v[1].re = tst;
            v[1].im = 0.0;
            if (std::abs(h->data[(m + h->size[0] * (m - 2)) - 1].re) * std::abs
                (tst) <= 2.2204460492503131E-16 * ((std::abs(sc.re) + std::abs
                  (sc.im)) * ((std::abs(h->data[(m + h->size[0] * (m - 1)) - 1].
                    re) + std::abs(h->data[(m + h->size[0] * (m - 1)) - 1].im))
                              + (std::abs(h->data[m + h->size[0] * m].re) + std::
                                 abs(h->data[m + h->size[0] * m].im))))) {
              goto70 = true;
              exitg3 = true;
            } else {
              m--;
            }
          }

          if (!goto70) {
            sc.re = h->data[k + h->size[0] * k].re - t_re;
            sc.im = h->data[k + h->size[0] * k].im - t_im;
            tst = h->data[(k + h->size[0] * k) + 1].re;
            s = (std::abs(sc.re) + std::abs(sc.im)) + std::abs(tst);
            if (sc.im == 0.0) {
              v[0].re = sc.re / s;
              v[0].im = 0.0;
            } else if (sc.re == 0.0) {
              v[0].re = 0.0;
              v[0].im = sc.im / s;
            } else {
              v[0].re = sc.re / s;
              v[0].im = sc.im / s;
            }

            tst /= s;
            v[1].re = tst;
            v[1].im = 0.0;
          }

          for (b_k = m; b_k <= i; b_k++) {
            if (b_k > m) {
              v[0] = h->data[(b_k + h->size[0] * (b_k - 2)) - 1];
              v[1] = h->data[b_k + h->size[0] * (b_k - 2)];
            }

            sc = xzlarfg(&v[0], &v[1]);
            if (b_k > m) {
              h->data[(b_k + h->size[0] * (b_k - 2)) - 1] = v[0];
              h->data[b_k + h->size[0] * (b_k - 2)].re = 0.0;
              h->data[b_k + h->size[0] * (b_k - 2)].im = 0.0;
            }

            t_re = v[1].re;
            t_im = v[1].im;
            tst = sc.re * v[1].re - sc.im * v[1].im;
            for (j = b_k; j <= n; j++) {
              aa = sc.re * h->data[(b_k + h->size[0] * (j - 1)) - 1].re - -sc.im
                * h->data[(b_k + h->size[0] * (j - 1)) - 1].im;
              htmp1 = sc.re * h->data[(b_k + h->size[0] * (j - 1)) - 1].im +
                -sc.im * h->data[(b_k + h->size[0] * (j - 1)) - 1].re;
              x2.re = aa + tst * h->data[b_k + h->size[0] * (j - 1)].re;
              x2.im = htmp1 + tst * h->data[b_k + h->size[0] * (j - 1)].im;
              h->data[(b_k + h->size[0] * (j - 1)) - 1].re -= x2.re;
              h->data[(b_k + h->size[0] * (j - 1)) - 1].im -= x2.im;
              h->data[b_k + h->size[0] * (j - 1)].re -= x2.re * t_re - x2.im *
                t_im;
              h->data[b_k + h->size[0] * (j - 1)].im -= x2.re * t_im + x2.im *
                t_re;
            }

            if (b_k + 2 < i + 1) {
              u1 = b_k + 1;
            } else {
              u1 = i;
            }

            for (j = 0; j <= u1; j++) {
              aa = sc.re * h->data[j + h->size[0] * (b_k - 1)].re - sc.im *
                h->data[j + h->size[0] * (b_k - 1)].im;
              htmp1 = sc.re * h->data[j + h->size[0] * (b_k - 1)].im + sc.im *
                h->data[j + h->size[0] * (b_k - 1)].re;
              x2.re = aa + tst * h->data[j + h->size[0] * b_k].re;
              x2.im = htmp1 + tst * h->data[j + h->size[0] * b_k].im;
              h->data[j + h->size[0] * (b_k - 1)].re -= x2.re;
              h->data[j + h->size[0] * (b_k - 1)].im -= x2.im;
              h->data[j + h->size[0] * b_k].re -= x2.re * t_re - x2.im * -t_im;
              h->data[j + h->size[0] * b_k].im -= x2.re * -t_im + x2.im * t_re;
            }

            if ((b_k == m) && (m > k + 1)) {
              sc.re = 1.0 - sc.re;
              sc.im = 0.0 - sc.im;
              htmp1 = rt_hypotd_snf(sc.re, sc.im);
              if (sc.im == 0.0) {
                sc.re /= htmp1;
                sc.im = 0.0;
              } else if (sc.re == 0.0) {
                sc.re = 0.0;
                sc.im /= htmp1;
              } else {
                sc.re /= htmp1;
                sc.im /= htmp1;
              }

              tst = h->data[m + h->size[0] * (m - 1)].re;
              aa = h->data[m + h->size[0] * (m - 1)].im;
              h->data[m + h->size[0] * (m - 1)].re = tst * sc.re - aa * -sc.im;
              h->data[m + h->size[0] * (m - 1)].im = tst * -sc.im + aa * sc.re;
              if (m + 2 <= i + 1) {
                tst = h->data[(m + h->size[0] * m) + 1].re;
                aa = h->data[(m + h->size[0] * m) + 1].im;
                h->data[(m + h->size[0] * m) + 1].re = tst * sc.re - aa * sc.im;
                h->data[(m + h->size[0] * m) + 1].im = tst * sc.im + aa * sc.re;
              }

              for (j = m; j <= i + 1; j++) {
                if (j != m + 1) {
                  if (n > j) {
                    b_xscal(n - j, sc, h, j + j * ldh, ldh);
                  }

                  x2.re = sc.re;
                  x2.im = -sc.im;
                  xscal(j - 1, x2, h, 1 + (j - 1) * ldh);
                }
              }
            }
          }

          sc = h->data[i + h->size[0] * (i - 1)];
          if (h->data[i + h->size[0] * (i - 1)].im != 0.0) {
            tst = rt_hypotd_snf(h->data[i + h->size[0] * (i - 1)].re, h->data[i
                                + h->size[0] * (i - 1)].im);
            h->data[i + h->size[0] * (i - 1)].re = tst;
            h->data[i + h->size[0] * (i - 1)].im = 0.0;
            if (sc.im == 0.0) {
              sc.re /= tst;
              sc.im = 0.0;
            } else if (sc.re == 0.0) {
              sc.re = 0.0;
              sc.im /= tst;
            } else {
              sc.re /= tst;
              sc.im /= tst;
            }

            if (n > i + 1) {
              x2.re = sc.re;
              x2.im = -sc.im;
              b_xscal((n - i) - 1, x2, h, (i + (i + 1) * ldh) + 1, ldh);
            }

            xscal(i, sc, h, 1 + i * ldh);
          }

          its++;
        }
      }

      if (!goto140) {
        info = i + 1;
        exitg1 = true;
      } else {
        i = L;
      }
    }
  }

  return info;
}

/* End of code generation (xzhseqr.cpp) */
