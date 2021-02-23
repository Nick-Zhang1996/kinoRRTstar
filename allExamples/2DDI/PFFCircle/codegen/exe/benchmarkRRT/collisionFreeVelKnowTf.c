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
#include "DI_stateFreeVel.h"
#include "linspace.h"
#include "benchmarkRRT_data.h"

/* Function Definitions */
void collisionFreeVelKnowTf(const double parent_data[], const double node[4],
  const creal_T *Tf, double *collision_flag, double xf[4])
{
  creal_T t[12];
  int i;
  double b_state[40];
  int j;
  int exitg3;
  int exitg2;
  double c_state;
  int exitg1;
  int state_idx_0_tmp;
  static const signed char iv7[5] = { 6, 15, 10, 5, 15 };

  double state_idx_0;
  static const signed char iv8[5] = { 6, 14, 11, 15, 5 };

  *collision_flag = 0.0;
  linspace(*Tf, t);
  for (i = 0; i < 10; i++) {
    DI_stateFreeVel(t[i + 1], *Tf, parent_data[0], parent_data[1], parent_data[2],
                    parent_data[3], node[0], node[1], *(double (*)[4])&b_state[i
                    << 2]);
  }

  DI_stateFreeVel(*Tf, *Tf, parent_data[0], parent_data[1], parent_data[2],
                  parent_data[3], node[0], node[1], xf);
  j = 0;
  do {
    exitg3 = 0;
    if (j < 10) {
      i = 0;
      do {
        exitg2 = 0;
        if (i < 2) {
          c_state = b_state[i + (j << 2)];
          if ((c_state > 20.0) || (c_state < 0.0)) {
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
        exitg3 = 1;
      } else {
        do {
          exitg1 = 0;
          if (i < 5) {
            state_idx_0_tmp = j << 2;
            c_state = b_state[state_idx_0_tmp] - (double)iv7[i];
            c_state *= c_state;
            state_idx_0 = c_state;
            c_state = b_state[1 + state_idx_0_tmp] - (double)iv8[i];
            c_state *= c_state;
            if (state_idx_0 + c_state < (dv0[i] + 0.1) * (dv0[i] + 0.1)) {
              /*  (norm([p(1);p(2)]-[world.cx(i); world.cy(i)])<=1*world.radius(i)) */
              *collision_flag = 1.0;
              exitg1 = 1;
            } else {
              i++;
            }
          } else {
            j++;
            exitg1 = 2;
          }
        } while (exitg1 == 0);

        if (exitg1 == 1) {
          exitg3 = 1;
        }
      }
    } else {
      /* %%%% dim=3 case has not been implemented yet %%%%% */
      exitg3 = 1;
    }
  } while (exitg3 == 0);
}

/* End of code generation (collisionFreeVelKnowTf.c) */
