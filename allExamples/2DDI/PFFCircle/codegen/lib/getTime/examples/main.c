/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: main.c
 *
 * MATLAB Coder version            : 4.1
 * C/C++ source code generated on  : 22-Feb-2021 23:20:10
 */

/*************************************************************************/
/* This automatically generated example C main file shows how to call    */
/* entry-point functions that MATLAB Coder generated. You must customize */
/* this file for your application. Do not modify this file directly.     */
/* Instead, make a copy of this file, modify it, and integrate it into   */
/* your development environment.                                         */
/*                                                                       */
/* This file initializes entry-point function arguments to a default     */
/* size and value before calling the entry-point functions. It does      */
/* not store or use any values returned from the entry-point functions.  */
/* If necessary, it does pre-allocate memory for returned values.        */
/* You can use this file as a starting point for a main function that    */
/* you can deploy in your application.                                   */
/*                                                                       */
/* After you copy the file, and before you deploy it, you must make the  */
/* following changes:                                                    */
/* * For variable-size function arguments, change the example sizes to   */
/* the sizes that your application requires.                             */
/* * Change the example values of function arguments to the values that  */
/* your application requires.                                            */
/* * If the entry-point functions return values, store these values or   */
/* otherwise use them as required by your application.                   */
/*                                                                       */
/*************************************************************************/
/* Include Files */
#include "getTime.h"
#include "main.h"
#include "getTime_terminate.h"
#include "getTime_initialize.h"

/* Function Declarations */
static void main_getTime(void);

/* Function Definitions */

/*
 * Arguments    : void
 * Return Type  : void
 */
static void main_getTime(void)
{
  double current_time_in_sec;

  /* Call the entry-point 'getTime'. */
  current_time_in_sec = getTime();
}

/*
 * Arguments    : int argc
 *                const char * const argv[]
 * Return Type  : int
 */
int main(int argc, const char * const argv[])
{
  (void)argc;
  (void)argv;

  /* Initialize the application.
     You do not need to do this more than one time. */
  getTime_initialize();

  /* Invoke the entry-point functions.
     You can call entry-point functions multiple times. */
  main_getTime();

  /* Terminate the application.
     You do not need to do this more than one time. */
  getTime_terminate();
  return 0;
}

/*
 * File trailer for main.c
 *
 * [EOF]
 */
