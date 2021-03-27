/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * xzhseqr.c
 *
 * Code generation for function 'xzhseqr'
 *
 */

/* Include files */
#include <math.h>
#include "rt_nonfinite.h"
#include "benchmarkRRT.h"
#include "xzhseqr.h"
#include "quadf_cost.h"
#include "xzlarfg.h"
#include "sqrt.h"
#include "benchmarkRRT_rtwutil.h"

/* Function Definitions */
int eml_zlahqr(creal_T h_data[], int h_size[2])
{
  int info;
  int n;
  int ldh;
  int i11;
  int j;
  int u1;
  int ix0_tmp;
  int i;
  double SMLNUM;
  boolean_T exitg1;
  double tst;
  double ab;
  double bb;
  int L;
  boolean_T goto140;
  creal_T sc;
  int its;
  boolean_T exitg2;
  int k;
  boolean_T exitg3;
  double ba;
  int ix0;
  double t_re;
  boolean_T goto70;
  creal_T x2;
  int m;
  double u_re;
  double u_im;
  int b_k;
  double s;
  double aa;
  creal_T v[2];
  double b_u_re;
  n = h_size[0];
  ldh = h_size[0];
  info = 0;
  if (1 != h_size[0]) {
    i11 = h_size[0];
    for (j = 0; j <= i11 - 4; j++) {
      u1 = j + h_size[0] * j;
      ix0_tmp = u1 + 2;
      h_data[ix0_tmp].re = 0.0;
      h_data[ix0_tmp].im = 0.0;
      u1 += 3;
      h_data[u1].re = 0.0;
      h_data[u1].im = 0.0;
    }

    if (1 <= n - 2) {
      i11 = (n + h_size[0] * (n - 3)) - 1;
      h_data[i11].re = 0.0;
      h_data[i11].im = 0.0;
    }

    for (i = 2; i <= n; i++) {
      i11 = (i + h_size[0] * (i - 2)) - 1;
      if (h_data[i11].im != 0.0) {
        tst = h_data[(i + h_size[0] * (i - 2)) - 1].re;
        ab = h_data[(i + h_size[0] * (i - 2)) - 1].im;
        bb = fabs(h_data[(i + h_size[0] * (i - 2)) - 1].re) + fabs(h_data[(i +
          h_size[0] * (i - 2)) - 1].im);
        if (ab == 0.0) {
          sc.re = tst / bb;
          sc.im = 0.0;
        } else if (tst == 0.0) {
          sc.re = 0.0;
          sc.im = ab / bb;
        } else {
          sc.re = tst / bb;
          sc.im = ab / bb;
        }

        bb = rt_hypotd_snf(sc.re, sc.im);
        if (-sc.im == 0.0) {
          sc.re /= bb;
          sc.im = 0.0;
        } else if (sc.re == 0.0) {
          sc.re = 0.0;
          sc.im = -sc.im / bb;
        } else {
          sc.re /= bb;
          sc.im = -sc.im / bb;
        }

        h_data[i11].re = rt_hypotd_snf(h_data[(i + h_size[0] * (i - 2)) - 1].re,
          h_data[(i + h_size[0] * (i - 2)) - 1].im);
        h_data[i11].im = 0.0;
        ix0_tmp = (i - 1) * ldh;
        ix0 = i + ix0_tmp;
        i11 = ix0 + ldh * (n - i);
        for (k = ix0; ldh < 0 ? k >= i11 : k <= i11; k += ldh) {
          ab = h_data[k - 1].re;
          tst = h_data[k - 1].im;
          h_data[k - 1].re = sc.re * ab - sc.im * tst;
          h_data[k - 1].im = sc.re * tst + sc.im * ab;
        }

        ix0 = ix0_tmp + 1;
        sc.im = -sc.im;
        u1 = i + 1;
        if (n < u1) {
          u1 = n;
        }

        i11 = ix0_tmp + u1;
        for (k = ix0; k <= i11; k++) {
          ab = h_data[k - 1].re;
          tst = h_data[k - 1].im;
          h_data[k - 1].re = sc.re * ab - sc.im * tst;
          h_data[k - 1].im = sc.re * tst + sc.im * ab;
        }
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
      while ((!exitg2) && (its < 301)) {
        k = i;
        exitg3 = false;
        while ((!exitg3) && (k + 1 > L + 2)) {
          i11 = k + h_size[0] * (k - 1);
          ab = fabs(h_data[i11].re);
          ba = ab + fabs(h_data[k + h_size[0] * (k - 1)].im);
          if (ba <= SMLNUM) {
            exitg3 = true;
          } else {
            u1 = k + h_size[0] * k;
            bb = fabs(h_data[u1].re) + fabs(h_data[k + h_size[0] * k].im);
            tst = (fabs(h_data[i11 - 1].re) + fabs(h_data[(k + h_size[0] * (k -
                      1)) - 1].im)) + bb;
            if (tst == 0.0) {
              if (k - 1 >= 1) {
                tst = fabs(h_data[(k + h_size[0] * (k - 2)) - 1].re);
              }

              if (k + 2 <= n) {
                tst += fabs(h_data[u1 + 1].re);
              }
            }

            if (ab <= 2.2204460492503131E-16 * tst) {
              tst = fabs(h_data[u1 - 1].re) + fabs(h_data[(k + h_size[0] * k) -
                1].im);
              if (ba > tst) {
                ab = ba;
                ba = tst;
              } else {
                ab = tst;
              }

              tst = fabs(h_data[(k + h_size[0] * (k - 1)) - 1].re - h_data[k +
                         h_size[0] * k].re) + fabs(h_data[(k + h_size[0] * (k -
                1)) - 1].im - h_data[k + h_size[0] * k].im);
              if (bb > tst) {
                aa = bb;
                bb = tst;
              } else {
                aa = tst;
              }

              s = aa + ab;
              if (ba * (ab / s) <= fmax(SMLNUM, 2.2204460492503131E-16 * (bb *
                    (aa / s)))) {
                exitg3 = true;
              } else {
                k--;
              }
            } else {
              k--;
            }
          }
        }

        L = k - 1;
        if (k + 1 > 1) {
          h_data[k + h_size[0] * (k - 1)].re = 0.0;
          h_data[k + h_size[0] * (k - 1)].im = 0.0;
        }

        if (k + 1 >= i + 1) {
          goto140 = true;
          exitg2 = true;
        } else {
          if (its == 10) {
            t_re = 0.75 * fabs(h_data[(k + h_size[0] * k) + 1].re) + h_data[k +
              h_size[0] * k].re;
            ba = h_data[k + h_size[0] * k].im;
          } else if (its == 20) {
            t_re = 0.75 * fabs(h_data[i + h_size[0] * (i - 1)].re) + h_data[i +
              h_size[0] * i].re;
            ba = h_data[i + h_size[0] * i].im;
          } else {
            ix0_tmp = i + h_size[0] * i;
            t_re = h_data[ix0_tmp].re;
            ba = h_data[i + h_size[0] * i].im;
            x2 = h_data[ix0_tmp - 1];
            b_sqrt(&x2);
            u1 = i + h_size[0] * (i - 1);
            sc = h_data[u1];
            b_sqrt(&sc);
            u_re = x2.re * sc.re - x2.im * sc.im;
            u_im = x2.re * sc.im + x2.im * sc.re;
            s = fabs(u_re) + fabs(u_im);
            if (s != 0.0) {
              t_re = 0.5 * (h_data[u1 - 1].re - h_data[i + h_size[0] * i].re);
              ba = 0.5 * (h_data[(i + h_size[0] * (i - 1)) - 1].im - h_data[i +
                          h_size[0] * i].im);
              tst = fabs(t_re) + fabs(ba);
              s = fmax(s, tst);
              if (ba == 0.0) {
                x2.re = t_re / s;
                x2.im = 0.0;
              } else if (t_re == 0.0) {
                x2.re = 0.0;
                x2.im = ba / s;
              } else {
                x2.re = t_re / s;
                x2.im = ba / s;
              }

              ab = x2.re;
              aa = x2.re;
              x2.re = x2.re * x2.re - x2.im * x2.im;
              x2.im = ab * x2.im + x2.im * aa;
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
                if (ba == 0.0) {
                  x2.re = t_re / tst;
                  x2.im = 0.0;
                } else if (t_re == 0.0) {
                  x2.re = 0.0;
                  x2.im = ba / tst;
                } else {
                  x2.re = t_re / tst;
                  x2.im = ba / tst;
                }

                if (x2.re * sc.re + x2.im * sc.im < 0.0) {
                  sc.re = -sc.re;
                  sc.im = -sc.im;
                }
              }

              bb = t_re + sc.re;
              aa = ba + sc.im;
              if (aa == 0.0) {
                if (u_im == 0.0) {
                  b_u_re = u_re / bb;
                  tst = 0.0;
                } else if (u_re == 0.0) {
                  b_u_re = 0.0;
                  tst = u_im / bb;
                } else {
                  b_u_re = u_re / bb;
                  tst = u_im / bb;
                }
              } else if (bb == 0.0) {
                if (u_re == 0.0) {
                  b_u_re = u_im / aa;
                  tst = 0.0;
                } else if (u_im == 0.0) {
                  b_u_re = 0.0;
                  tst = -(u_re / aa);
                } else {
                  b_u_re = u_im / aa;
                  tst = -(u_re / aa);
                }
              } else {
                ba = fabs(bb);
                tst = fabs(aa);
                if (ba > tst) {
                  s = aa / bb;
                  tst = bb + s * aa;
                  b_u_re = (u_re + s * u_im) / tst;
                  tst = (u_im - s * u_re) / tst;
                } else if (tst == ba) {
                  if (bb > 0.0) {
                    ab = 0.5;
                  } else {
                    ab = -0.5;
                  }

                  if (aa > 0.0) {
                    tst = 0.5;
                  } else {
                    tst = -0.5;
                  }

                  b_u_re = (u_re * ab + u_im * tst) / ba;
                  tst = (u_im * ab - u_re * tst) / ba;
                } else {
                  s = bb / aa;
                  tst = aa + s * bb;
                  b_u_re = (s * u_re + u_im) / tst;
                  tst = (s * u_im - u_re) / tst;
                }
              }

              t_re = h_data[i + h_size[0] * i].re - (u_re * b_u_re - u_im * tst);
              ba = h_data[i + h_size[0] * i].im - (u_re * tst + u_im * b_u_re);
            }
          }

          goto70 = false;
          m = i;
          exitg3 = false;
          while ((!exitg3) && (m > k + 1)) {
            u1 = m + h_size[0] * (m - 1);
            sc.re = h_data[u1 - 1].re - t_re;
            sc.im = h_data[(m + h_size[0] * (m - 1)) - 1].im - ba;
            ab = h_data[u1].re;
            s = (fabs(sc.re) + fabs(sc.im)) + fabs(ab);
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

            ab /= s;
            v[0] = sc;
            v[1].re = ab;
            v[1].im = 0.0;
            if (fabs(h_data[(m + h_size[0] * (m - 2)) - 1].re) * fabs(ab) <=
                2.2204460492503131E-16 * ((fabs(sc.re) + fabs(sc.im)) * ((fabs
                   (h_data[(m + h_size[0] * (m - 1)) - 1].re) + fabs(h_data[(m +
                     h_size[0] * (m - 1)) - 1].im)) + (fabs(h_data[m + h_size[0]
                    * m].re) + fabs(h_data[m + h_size[0] * m].im))))) {
              goto70 = true;
              exitg3 = true;
            } else {
              m--;
            }
          }

          if (!goto70) {
            sc.re = h_data[k + h_size[0] * k].re - t_re;
            sc.im = h_data[k + h_size[0] * k].im - ba;
            ab = h_data[(k + h_size[0] * k) + 1].re;
            s = (fabs(sc.re) + fabs(sc.im)) + fabs(ab);
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

            ab /= s;
            v[1].re = ab;
            v[1].im = 0.0;
          }

          for (b_k = m; b_k <= i; b_k++) {
            if (b_k > m) {
              u1 = b_k + h_size[0] * (b_k - 2);
              v[0] = h_data[u1 - 1];
              v[1] = h_data[u1];
            }

            sc = xzlarfg(&v[0], &v[1]);
            if (b_k > m) {
              h_data[(b_k + h_size[0] * (b_k - 2)) - 1] = v[0];
              h_data[b_k + h_size[0] * (b_k - 2)].re = 0.0;
              h_data[b_k + h_size[0] * (b_k - 2)].im = 0.0;
            }

            t_re = v[1].re;
            ba = v[1].im;
            tst = sc.re * v[1].re - sc.im * v[1].im;
            for (j = b_k; j <= n; j++) {
              ix0_tmp = b_k + h_size[0] * (j - 1);
              u1 = ix0_tmp - 1;
              x2.re = (sc.re * h_data[u1].re - -sc.im * h_data[(b_k + h_size[0] *
                        (j - 1)) - 1].im) + tst * h_data[ix0_tmp].re;
              x2.im = (sc.re * h_data[(b_k + h_size[0] * (j - 1)) - 1].im +
                       -sc.im * h_data[(b_k + h_size[0] * (j - 1)) - 1].re) +
                tst * h_data[b_k + h_size[0] * (j - 1)].im;
              h_data[u1].re = h_data[(b_k + h_size[0] * (j - 1)) - 1].re - x2.re;
              h_data[u1].im = h_data[(b_k + h_size[0] * (j - 1)) - 1].im - x2.im;
              h_data[ix0_tmp].re = h_data[b_k + h_size[0] * (j - 1)].re - (x2.re
                * t_re - x2.im * ba);
              h_data[ix0_tmp].im = h_data[b_k + h_size[0] * (j - 1)].im - (x2.re
                * ba + x2.im * t_re);
            }

            if (b_k + 2 < i + 1) {
              i11 = b_k + 1;
            } else {
              i11 = i;
            }

            for (j = 0; j <= i11; j++) {
              ix0_tmp = j + h_size[0] * (b_k - 1);
              u1 = j + h_size[0] * b_k;
              x2.re = (sc.re * h_data[ix0_tmp].re - sc.im * h_data[j + h_size[0]
                       * (b_k - 1)].im) + tst * h_data[u1].re;
              x2.im = (sc.re * h_data[j + h_size[0] * (b_k - 1)].im + sc.im *
                       h_data[j + h_size[0] * (b_k - 1)].re) + tst * h_data[j +
                h_size[0] * b_k].im;
              h_data[ix0_tmp].re = h_data[j + h_size[0] * (b_k - 1)].re - x2.re;
              h_data[ix0_tmp].im = h_data[j + h_size[0] * (b_k - 1)].im - x2.im;
              h_data[u1].re = h_data[j + h_size[0] * b_k].re - (x2.re * t_re -
                x2.im * -ba);
              h_data[u1].im = h_data[j + h_size[0] * b_k].im - (x2.re * -ba +
                x2.im * t_re);
            }

            if ((b_k == m) && (m > k + 1)) {
              bb = rt_hypotd_snf(1.0 - sc.re, 0.0 - sc.im);
              if (0.0 - sc.im == 0.0) {
                t_re = (1.0 - sc.re) / bb;
                ba = 0.0;
              } else if (1.0 - sc.re == 0.0) {
                t_re = 0.0;
                ba = (0.0 - sc.im) / bb;
              } else {
                t_re = (1.0 - sc.re) / bb;
                ba = (0.0 - sc.im) / bb;
              }

              ab = h_data[m + h_size[0] * (m - 1)].re;
              tst = h_data[m + h_size[0] * (m - 1)].im;
              h_data[m + h_size[0] * (m - 1)].re = ab * t_re - tst * -ba;
              h_data[m + h_size[0] * (m - 1)].im = ab * -ba + tst * t_re;
              if (m + 2 <= i + 1) {
                u1 = (m + h_size[0] * m) + 1;
                ab = h_data[u1].re;
                tst = h_data[(m + h_size[0] * m) + 1].im;
                h_data[u1].re = ab * t_re - tst * ba;
                h_data[u1].im = ab * ba + tst * t_re;
              }

              for (j = m; j <= i + 1; j++) {
                if (j != m + 1) {
                  if (n > j) {
                    ix0 = j + j * ldh;
                    i11 = ix0 + ldh * ((n - j) - 1);
                    for (u1 = ix0; ldh < 0 ? u1 >= i11 : u1 <= i11; u1 += ldh) {
                      ab = h_data[u1 - 1].re;
                      tst = h_data[u1 - 1].im;
                      h_data[u1 - 1].re = t_re * ab - ba * tst;
                      h_data[u1 - 1].im = t_re * tst + ba * ab;
                    }
                  }

                  ix0_tmp = (j - 1) * ldh;
                  ix0 = ix0_tmp + 1;
                  i11 = (ix0_tmp + j) - 1;
                  for (u1 = ix0; u1 <= i11; u1++) {
                    ab = h_data[u1 - 1].re;
                    tst = h_data[u1 - 1].im;
                    h_data[u1 - 1].re = t_re * ab - -ba * tst;
                    h_data[u1 - 1].im = t_re * tst + -ba * ab;
                  }
                }
              }
            }
          }

          t_re = h_data[i + h_size[0] * (i - 1)].re;
          ba = h_data[i + h_size[0] * (i - 1)].im;
          if (h_data[i + h_size[0] * (i - 1)].im != 0.0) {
            tst = rt_hypotd_snf(h_data[i + h_size[0] * (i - 1)].re, h_data[i +
                                h_size[0] * (i - 1)].im);
            h_data[i + h_size[0] * (i - 1)].re = tst;
            h_data[i + h_size[0] * (i - 1)].im = 0.0;
            if (ba == 0.0) {
              t_re /= tst;
              ba = 0.0;
            } else if (t_re == 0.0) {
              t_re = 0.0;
              ba /= tst;
            } else {
              t_re /= tst;
              ba /= tst;
            }

            if (n > i + 1) {
              ix0 = (i + (i + 1) * ldh) + 1;
              i11 = ix0 + ldh * ((n - i) - 2);
              for (k = ix0; ldh < 0 ? k >= i11 : k <= i11; k += ldh) {
                ab = h_data[k - 1].re;
                tst = h_data[k - 1].im;
                h_data[k - 1].re = t_re * ab - -ba * tst;
                h_data[k - 1].im = t_re * tst + -ba * ab;
              }
            }

            ix0_tmp = i * ldh;
            ix0 = ix0_tmp + 1;
            i11 = ix0_tmp + i;
            for (k = ix0; k <= i11; k++) {
              ab = h_data[k - 1].re;
              tst = h_data[k - 1].im;
              h_data[k - 1].re = t_re * ab - ba * tst;
              h_data[k - 1].im = t_re * tst + ba * ab;
            }
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

/* End of code generation (xzhseqr.c) */
