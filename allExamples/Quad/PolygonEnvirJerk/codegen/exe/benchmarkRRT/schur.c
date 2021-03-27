/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * schur.c
 *
 * Code generation for function 'schur'
 *
 */

/* Include files */
#include <string.h>
#include "rt_nonfinite.h"
#include "benchmarkRRT.h"
#include "schur.h"
#include "xzhseqr.h"
#include "xgehrd.h"
#include "anyNonFinite.h"

/* Function Definitions */
void schur(const creal_T A_data[], const int A_size[2], creal_T V_data[], int
           V_size[2])
{
  signed char unnamed_idx_0;
  signed char unnamed_idx_1;
  int m;
  int i1;
  int istart;
  int jend;
  int j;
  int i;
  if (anyNonFinite(A_data, A_size)) {
    unnamed_idx_0 = (signed char)A_size[0];
    unnamed_idx_1 = (signed char)A_size[1];
    V_size[0] = unnamed_idx_0;
    V_size[1] = unnamed_idx_1;
    m = unnamed_idx_0 * unnamed_idx_1;
    for (i1 = 0; i1 < m; i1++) {
      V_data[i1].re = rtNaN;
      V_data[i1].im = 0.0;
    }

    m = unnamed_idx_0;
    if (1 < unnamed_idx_0) {
      istart = 2;
      if (unnamed_idx_0 - 2 < unnamed_idx_1 - 1) {
        jend = unnamed_idx_0 - 1;
      } else {
        jend = unnamed_idx_1;
      }

      for (j = 0; j < jend; j++) {
        for (i = istart; i <= m; i++) {
          i1 = (i + unnamed_idx_0 * j) - 1;
          V_data[i1].re = 0.0;
          V_data[i1].im = 0.0;
        }

        istart++;
      }
    }
  } else {
    V_size[0] = A_size[0];
    V_size[1] = A_size[1];
    m = A_size[0] * A_size[1];
    if (0 <= m - 1) {
      memcpy(&V_data[0], &A_data[0], (unsigned int)(m * (int)sizeof(creal_T)));
    }

    xgehrd(V_data, V_size);
    eml_zlahqr(V_data, V_size);
    m = V_size[0];
    if (3 < V_size[0]) {
      istart = 4;
      if (V_size[0] - 4 < V_size[1] - 1) {
        jend = V_size[0] - 3;
      } else {
        jend = V_size[1];
      }

      for (j = 0; j < jend; j++) {
        for (i = istart; i <= m; i++) {
          i1 = (i + V_size[0] * j) - 1;
          V_data[i1].re = 0.0;
          V_data[i1].im = 0.0;
        }

        istart++;
      }
    }
  }
}

/* End of code generation (schur.c) */
