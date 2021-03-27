/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * collisionFreeVelKnowTf.c
 *
 * Code generation for function 'collisionFreeVelKnowTf'
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "benchmarkRRT.h"
#include "collisionFreeVelKnowTf.h"
#include "quadf_statePFFull.h"
#include "quadf_statePFF.h"
#include "linspace.h"

/* Function Definitions */
void collisionFreeVelKnowTf(const double parent_data[], const double node[9],
  const creal_T *Tf, double *collision_flag, double xf[9])
{
  creal_T t[12];
  int i;
  double b_state[30];
  int j;
  boolean_T exitg3;
  int exitg2;
  double d3;
  static const signed char iv14[3] = { 20, 10, 10 };

  int exitg1;
  static const signed char iv15[4] = { 6, 6, 12, 12 };

  static const signed char iv16[4] = { 0, 5, 5, 0 };

  static const signed char iv17[4] = { 0, 0, 0, 5 };

  static const signed char iv18[4] = { 10, 5, 10, 5 };

  *collision_flag = 0.0;
  linspace(*Tf, t);
  for (i = 0; i < 10; i++) {
    quadf_statePFF(t[i + 1], *Tf, parent_data[0], parent_data[1], parent_data[2],
                   parent_data[3], parent_data[4], parent_data[5], parent_data[6],
                   parent_data[7], parent_data[8], node[0], node[1], node[2],
                   *(double (*)[3])&b_state[3 * i]);
  }

  quadf_statePFFull(*Tf, *Tf, parent_data[0], parent_data[1], parent_data[2],
                    parent_data[3], parent_data[4], parent_data[5], parent_data
                    [6], parent_data[7], parent_data[8], node[0], node[1], node
                    [2], xf);
  j = 0;
  exitg3 = false;
  while ((!exitg3) && (j < 10)) {
    i = 0;
    do {
      exitg2 = 0;
      if (i < 3) {
        d3 = b_state[i + 3 * j];
        if ((d3 > iv14[i]) || (d3 < 0.0)) {
          *collision_flag = 1.0;
          exitg2 = 1;
        } else {
          i++;
        }
      } else {
        /*  check each obstacle */
        i = 0;
        exitg2 = 2;
      }
    } while (exitg2 == 0);

    if (exitg2 == 1) {
      exitg3 = true;
    } else {
      do {
        exitg1 = 0;
        if (i < 4) {
          d3 = b_state[3 * j];
          if ((d3 >= iv15[i]) && (d3 <= (double)iv15[i] + 2.0)) {
            d3 = b_state[1 + 3 * j];
            if ((d3 >= iv16[i]) && (d3 <= (double)iv16[i] + 5.0)) {
              d3 = b_state[2 + 3 * j];
              if ((d3 >= iv17[i]) && (d3 <= iv17[i] + iv18[i])) {
                *collision_flag = 1.0;
                exitg1 = 2;
              } else {
                i++;
              }
            } else {
              i++;
            }
          } else {
            i++;
          }
        } else {
          j++;
          exitg1 = 1;
        }
      } while (exitg1 == 0);

      if (exitg1 == 1) {
      } else {
        exitg3 = true;
      }
    }
  }
}

/* End of code generation (collisionFreeVelKnowTf.c) */
