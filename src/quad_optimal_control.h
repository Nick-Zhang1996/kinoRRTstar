// optimal control solver for quadcopter case
#ifndef QUAD_OPTIMAL_CONTROL_H_
#define QUAD_OPTIMAL_CONTROL_H_

#include <assert.h>
#include <stdlib.h>
#include <algorithm>
#include <exception>
#include <vector>
#include <math.h>
#include <iostream>
#include <limits>
#include <unsupported/Eigen/Polynomials>


using std::vector;
using std::cout;

class err_NoValidSolution : public std::exception{};
class err_Unexcepted : public std::exception{};

class QuadOptimalControl{
  public:
  const int interior_point_count;
  double *buffer_interior_states, *buffer_interior_pos;
  // N: no of interior points when evaluating states between initial and final for collision checking
  QuadOptimalControl(const int N) : interior_point_count(N) {
    buffer_interior_states= new double[9*interior_point_count];
    buffer_interior_pos= new double[3*interior_point_count];
  }

  ~QuadOptimalControl(){
    delete[] buffer_interior_states;
    delete[] buffer_interior_pos;
  }

  // following group of functions deals with real root finding for a polynomial
  // nth degree polynomial is represented as a std::vector<double> with coeffs as elements
  // <an, an-1, an-2, ... a0> (size:n+1)
  // an*x^n + an-1 * x*n-1 + ... + a0
  //
  // find the smallest positive real root of a polynomial
  double minPositiveRoot(vector<double> coeffs);

  // find real roots of a <deg> degree polynomial 
  // in ascending order
  // OBSOLETE
  //vector<double> roots(vector<double> coeffs);

  // find derivative of a polynomial
  //vector<double> derivative(vector<double> coeffs);

  // find root of polynomial in interval
  // must contain one and only one root in interval
  //double falsePosition(vector<double> coeffs, double low_bound, double high_bound, double err=1e-4);

  // evaluate polynomial 
  /*
  inline double peval( vector<double> coeffs, double x){
    double val = 0.0;
    int deg = coeffs.size()-1;
    for (int i=0; i<(deg+1); i++){
      val += coeffs[i] * pow(x,deg-i);
    }
    return val;
  }
  */

    
  // optimal control related functions for quadcopter
  // x0i : initial state
  // x1i : final state
  // xi1-9: x,y,z,vx,vy,vz,ax,xy,az
  //
  // given initial and final full state, give cost
  double cost(double t_s, double x01, double  x02, double  x03, double  x04, double  x05, double  x06, double  x07, double  x08, double  x09, double  x11, double  x12, double  x13, double  x14, double  x15, double  x16, double  x17, double  x18, double  x19);

  //  find cost given partial state (xyz) only
  double costPartialFreeFinalState(double t_s, double x01, double x02, double x03, double x04, double x05, double x06, double x07, double x08, double x09, double x11, double x12, double x13);

  // find interior position vector, double  including endpoints
  // t_s : final time / segment time. initial state @ t=0, double  final state at t = t_s
  // return: dim (N, double n), double  N being number of interior points
  // Old name
  //double state(double t_s, double x01, double  x02, double  x03, double  x04, double  x05, double  x06, double  x07, double  x08, double  x09, double  x11, double  x12, double  x13, double  x14, double  x15, double  x16, double  x17, double  x18, double  x19);
  double* interiorPosition(double t_s, double x01, double  x02, double  x03, double  x04, double  x05, double  x06, double  x07, double  x08, double  x09, double  x11, double  x12, double  x13, double  x14, double  x15, double  x16, double  x17, double  x18, double  x19);

  // calculate only xyz
  //double statePartialFinalState(double t_s, double x01, double x02, double x03, double x04, double x05, double x06, double x07, double x08, double x09, double x11, double x12, double x13);
  double* interiorPositionPartialFinalState(double t_s, double x01, double x02, double x03, double x04, double x05, double x06, double x07, double x08, double x09, double x11, double x12, double x13);
  // calculate full state vector
  //double statePartialFinalStateFull(double t_s, double x01, double x02, double x03, double x04, double x05, double x06, double x07, double x08, double x09, double x11, double x12, double x13);
  double* interiorStatePartialFinalState(double t_s, double x01, double x02, double x03, double x04, double x05, double x06, double x07, double x08, double x09, double x11, double x12, double x13);

  // find time to go from initial to final state
  // this requires solving a polynomial
  double time(double x01, double  x02, double  x03, double  x04, double  x05, double  x06, double  x07, double  x08, double  x09, double  x11, double  x12, double  x13, double  x14, double  x15, double  x16, double  x17, double  x18, double  x19);
  double timePartialFinalState(double x01, double  x02, double  x03, double  x04, double  x05, double  x06, double  x07, double  x08, double  x09, double  x11, double  x12, double  x13);


  // test scripts
  void test_cost();
  void test_costPartialFreeFinalState();
  void test_interiorPosition();
  void test_interiorPositionPartialFinalState();
  void test_interiorStatePartialFinalState();
  void test_time();
  void test_timePartialFinalState();

};



#endif
