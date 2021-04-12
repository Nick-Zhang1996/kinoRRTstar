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


using std::vector;
using std::cout;

class err_NoValidSolution : public std::exception{};
class err_Unexcepted : public std::exception{};

class QuadOptimalControl{
  public:
  const int interior_point_count;
  float *buffer_interior_states, *buffer_interior_pos;
  // N: no of interior points when evaluating states between initial and final for collision checking
  QuadOptimalControl(const int N) : interior_point_count(N) {
    buffer_interior_states= new float[9*interior_point_count];
    buffer_interior_pos= new float[3*interior_point_count];
  }

  ~QuadOptimalControl(){
    delete[] buffer_interior_states;
    delete[] buffer_interior_pos;
  }

  // following group of functions deals with real root finding for a polynomial
  // nth degree polynomial is represented as a std::vector<float> with coeffs as elements
  // <an, an-1, an-2, ... a0> (size:n+1)
  // an*x^n + an-1 * x*n-1 + ... + a0
  //
  // find the smallest positive real root of a polynomial
  float minPositiveRoot(vector<float> coeffs);

  // find real roots of a <deg> degree polynomial 
  // in ascending order
  vector<float> roots(vector<float> coeffs);

  // find derivative of a polynomial
  vector<float> derivative(vector<float> coeffs);

  // find root of polynomial in interval
  // must contain one and only one root in interval
  float falsePosition(vector<float> coeffs, float low_bound, float high_bound, float err=1e-4);

  // evaluate polynomial 
  inline float peval( vector<float> coeffs, float x){
    float val = 0.0;
    int deg = coeffs.size()-1;
    for (int i=0; i<(deg+1); i++){
      val += coeffs.at(i) * pow(x,deg-i);
    }
    return val;
  }

    
  /*
  // optimal control related functions for quadcopter
  // x0i : initial state
  // x1i : final state
  // xi1-9: x,y,z,vx,vy,vz,ax,xy,az
  //
  // given initial and final full state, give cost
  float cost(float t_s, float  float x01, float  float  x02, float  float  x03, float  float  x04, float  float  x05, float  float  x06, float  float  x07, float  float  x08, float  float  x09, float  float  x11, float  float  x12, float  float  x13, float  float  x14, float  float  x15, float  float  x16, float  float  x17, float  float  x18, float  float  x19);

  //  find cost given partial state (xyz) only
  float costPartialFreeFinalState(t_s, float x01, float x02, float x03, float x04, float x05, float x06, float x07, float x08, float x09, float x11, float x12, float x13);

  // find interior position vector, float  including endpoints
  // t_s : final time / segment time. initial state @ t=0, float  final state at t = t_s
  // return: dim (N, float n), float  N being number of interior points
  // Old name
  //float state(float t_s, float x01, float  x02, float  x03, float  x04, float  x05, float  x06, float  x07, float  x08, float  x09, float  x11, float  x12, float  x13, float  x14, float  x15, float  x16, float  x17, float  x18, float  x19);
  float* interiorPosition(float t_s, float x01, float  x02, float  x03, float  x04, float  x05, float  x06, float  x07, float  x08, float  x09, float  x11, float  x12, float  x13, float  x14, float  x15, float  x16, float  x17, float  x18, float  x19);

  // calculate only xyz
  //float statePartialFinalState(float t_s, float x01, float x02, float x03, float x04, float x05, float x06, float x07, float x08, float x09, float x11, float x12, float x13);
  float* interiorPositionPartialFinalState(float t_s, float x01, float x02, float x03, float x04, float x05, float x06, float x07, float x08, float x09, float x11, float x12, float x13);
  // calculate full state vector
  //float statePartialFinalStateFull(float t_s, float x01, float x02, float x03, float x04, float x05, float x06, float x07, float x08, float x09, float x11, float x12, float x13);
  float* interiorStatePartialFinalState(float t_s, float x01, float x02, float x03, float x04, float x05, float x06, float x07, float x08, float x09, float x11, float x12, float x13);

  // find time to go from initial to final state
  float time(float x01, float  x02, float  x03, float  x04, float  x05, float  x06, float  x07, float  x08, float  x09, float  x11, float  x12, float  x13, float  x14, float  x15, float  x16, float  x17, float  x18, float  x19);
  float timePartialFinalState(float float x01, float  x02, float  x03, float  x04, float  x05, float  x06, float  x07, float  x08, float  x09, float  x11, float  x12, float  x13);


  */
};



#endif
