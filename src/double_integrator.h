// optimal control solver for doble integrator
#ifndef DOUBLE_INTEGRATOR_H
#define DOUBLE_INTEGRATOR_H

#include <assert.h>
#include <stdlib.h>
#include <algorithm>
#include <exception>
#include <vector>
#include <math.h>
#include <iostream>
#include <limits>
//#include <unsupported/Eigen/Polynomials>
#include <chrono>
#include <ctime>    
#include "common.h"
#include "tree.h"


using std::vector;
using std::cout;

class err_NoSolution: public std::exception{};
class err_Unexcepted : public std::exception{};
class err_TooManyIteration : public std::exception{};

class DoubleIntegrator{
  private:
    double total_time;
    int total_count;
  public:
  const int interior_point_count;
  double *buffer_interior_states, *buffer_interior_pos;
  double *buffer_interior_state;
  // N: no of interior points when evaluating states between initial and final for collision checking
  DoubleIntegrator(const int N) : 
    total_time(0.0),
    total_count(0),
    interior_point_count(N){
    buffer_interior_states= new double[6*interior_point_count];
    buffer_interior_pos= new double[3*interior_point_count];
    buffer_interior_state= new double[6];
  }

  ~DoubleIntegrator(){
    delete[] buffer_interior_state;
    delete[] buffer_interior_states;
    delete[] buffer_interior_pos;
  }

  void printTotalTime(){
    cout << "total time in roots " << total_time << "sec \n";
    cout << "total count " << total_count << "\n";
    cout << "avg time " << total_time/total_count << "\n";
  }

  // following group of functions deals with real root finding for a polynomial
  // nth degree polynomial is represented as a std::vector<double> with coeffs as elements
  // <an, an-1, an-2, ... a0> (size:n+1)
  // an*x^n + an-1 * x*n-1 + ... + a0
  //
  // find the smallest positive real root of a polynomial
  double minPositiveRoot(vector<double> coeffs);
  //double minPositiveRootEigen(vector<double> coeffs);
  // find root of polynomial coeff in low, high
  double bisect( vector<double> coeffs, double low, double high);

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
  inline double peval( vector<double> coeffs, double x){
    double val = 0.0;
    int deg = coeffs.size()-1;
    for (int i=0; i<(deg+1); i++){
      val += coeffs[i] * pow(x,deg-i);
    }
    return val;
  }

    
  // optimal control related functions for double integrator
  // x0i : initial state
  // x1i : final state
  // xi1-6: x,y,z,vx,vy,vz
  //
  // given initial and final full state, give cost
  double cost(double t_s, double x01, double  x02, double  x03, double  x04, double  x05, double  x06,  double  x11, double  x12, double  x13, double  x14, double  x15, double  x16);
  double cost(double t_s, Node& node_i, Node& node_f){ return cost(t_s, node_i.x, node_i.y, node_i.z, node_i.vx, node_i.vy, node_i.vz, node_f.x, node_f.y, node_f.z, node_f.vx, node_f.vy, node_f.vz);}

  //  find cost given partial state (xyz) only
  double costPartialFreeFinalState(double t_s, double x01, double x02, double x03, double x04, double x05, double x06,  double x11, double x12, double x13);
  double costPartialFreeFinalState(double t_s, Node& node_i, Node& node_f){ return costPartialFreeFinalState(t_s, node_i.x, node_i.y, node_i.z, node_i.vx, node_i.vy, node_i.vz,  node_f.x, node_f.y, node_f.z);}

  // find interior position vector, double  including endpoints
  // t_s : final time / segment time. initial state @ t=0, double  final state at t = t_s
  // return: dim (N, double n), double  N being number of interior points
  double* interiorPosition(double t_s, double x01, double  x02, double  x03, double  x04, double  x05, double  x06, double  x11, double  x12, double  x13, double  x14, double  x15, double  x16);
  double* interiorPosition(double t_s, Node& node_i, Node& node_f){ return interiorPosition(t_s, node_i.x, node_i.y, node_i.z, node_i.vx, node_i.vy, node_i.vz, node_f.x, node_f.y, node_f.z, node_f.vx, node_f.vy, node_f.vz);}

  // calculate only xyz
  //double statePartialFinalState(double t_s, double x01, double x02, double x03, double x04, double x05, double x06, double x11, double x12, double x13);
  double* interiorPositionPartialFinalState(double t_s, double x01, double x02, double x03, double x04, double x05, double x06, double x11, double x12, double x13);
  double* interiorPositionPartialFinalState(double t_s, Node& node_i, Node& node_f){ return interiorPositionPartialFinalState(t_s, node_i.x, node_i.y, node_i.z, node_i.vx, node_i.vy, node_i.vz, node_f.x, node_f.y, node_f.z);}
  // calculate full state vector
  //double statePartialFinalStateFull(double t_s, double x01, double x02, double x03, double x04, double x05, double x06, double x11, double x12, double x13);
  double* interiorStatePartialFinalState(double t_s, double x01, double x02, double x03, double x04, double x05, double x06, double x11, double x12, double x13);
  double* interiorStatePartialFinalState(double t_s, Node& node_i, Node& node_f){ return interiorStatePartialFinalState(t_s, node_i.x, node_i.y, node_i.z, node_i.vx, node_i.vy, node_i.vz, node_f.x, node_f.y, node_f.z);}

  // leave v in node_f, solve optimal control problem and fill in 
  // v for node_f
  void setFullStatePartialFinalState(double t_s, Node& node_i, Node& node_f);

  // given section time between two end states calculate full state at t and store in buffer_interior_state
  void state(double t,double t_s, double x01, double  x02, double  x03, double  x04, double  x05, double  x06, double  x11, double  x12, double  x13, double  x14, double  x15, double  x16);
  void state(double t,double t_s, Node& node_i, Node& node_f){ return state(t, t_s,  node_i.x, node_i.y, node_i.z, node_i.vx, node_i.vy, node_i.vz,  node_f.x, node_f.y, node_f.z,node_f.vx, node_f.vy, node_f.vz);}
  void state(double t,double t_s, Waypoint& node_i, Waypoint& node_f){ return state(t, t_s,  node_i.x, node_i.y, node_i.z, node_i.vx, node_i.vy, node_i.vz, node_f.x, node_f.y, node_f.z,node_f.vx, node_f.vy, node_f.vz);}


  // find time to go from initial to final state
  // this requires solving a polynomial
  double time(double x01, double  x02, double  x03, double  x04, double  x05, double  x06,  double  x11, double  x12, double  x13, double  x14, double  x15, double  x16);
  double time( Node& node_i, Node& node_f){ return time( node_i.x, node_i.y, node_i.z, node_i.vx, node_i.vy, node_i.vz, node_f.x, node_f.y, node_f.z, node_f.vx, node_f.vy, node_f.vz);}
  double time( Waypoint& node_i, Waypoint& node_f){ return time( node_i.x, node_i.y, node_i.z, node_i.vx, node_i.vy, node_i.vz,  node_f.x, node_f.y, node_f.z, node_f.vx, node_f.vy, node_f.vz);}

  double timePartialFinalState(double x01, double  x02, double  x03, double  x04, double  x05, double  x06, double  x11, double  x12, double  x13);
  double timePartialFinalState( Node& node_i, Node& node_f){ return timePartialFinalState( node_i.x, node_i.y, node_i.z, node_i.vx, node_i.vy, node_i.vz, node_f.x, node_f.y, node_f.z);}
  double timePartialFinalState( Waypoint& node_i, Waypoint& node_f){ return timePartialFinalState( node_i.x, node_i.y, node_i.z, node_i.vx, node_i.vy, node_i.vz, node_f.x, node_f.y, node_f.z);}


};



#endif
