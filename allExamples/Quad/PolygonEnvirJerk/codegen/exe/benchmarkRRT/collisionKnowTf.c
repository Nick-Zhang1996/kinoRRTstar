/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * collisionKnowTf.c
 *
 * Code generation for function 'collisionKnowTf'
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "benchmarkRRT.h"
#include "collisionKnowTf.h"
#include "quadf_state.h"
#include "linspace.h"

/* Function Definitions */
double collisionKnowTf(const double parent_data[], const double node_data[],
  const creal_T *Tf)
{
  double collision_flag;
  creal_T t[12];
  int i;
  int j;
  double b_state[30];
  boolean_T exitg3;
  int exitg2;
  double d4;
  static const signed char iv19[3] = { 20, 10, 10 };

  int exitg1;
  static const signed char iv20[4] = { 6, 6, 12, 12 };

  static const signed char iv21[4] = { 0, 5, 5, 0 };

  static const signed char iv22[4] = { 0, 0, 0, 5 };

  static const signed char iv23[4] = { 10, 5, 10, 5 };

  collision_flag = 0.0;
  linspace(*Tf, t);
  for (i = 0; i < 10; i++) {
    quadf_state(t[i + 1], *Tf, parent_data[0], parent_data[1], parent_data[2],
                parent_data[3], parent_data[4], parent_data[5], parent_data[6],
                parent_data[7], parent_data[8], node_data[0], node_data[1],
                node_data[2], node_data[3], node_data[4], node_data[5],
                node_data[6], node_data[7], node_data[8], *(double (*)[3])&
                b_state[3 * i]);
  }

  j = 0;
  exitg3 = false;
  while ((!exitg3) && (j < 10)) {
    i = 0;
    do {
      exitg2 = 0;
      if (i < 3) {
        d4 = b_state[i + 3 * j];
        if ((d4 > iv19[i]) || (d4 < 0.0)) {
          collision_flag = 1.0;
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
          d4 = b_state[3 * j];
          if ((d4 >= iv20[i]) && (d4 <= (double)iv20[i] + 2.0)) {
            d4 = b_state[1 + 3 * j];
            if ((d4 >= iv21[i]) && (d4 <= (double)iv21[i] + 5.0)) {
              d4 = b_state[2 + 3 * j];
              if ((d4 >= iv22[i]) && (d4 <= iv22[i] + iv23[i])) {
                collision_flag = 1.0;
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

  return collision_flag;
}

/* End of code generation (collisionKnowTf.c) */
