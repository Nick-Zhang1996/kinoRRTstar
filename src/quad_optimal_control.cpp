#include "quad_optimal_control.h"

// find the smallest positive real root of a polynomial
/*
double QuadOptimalControl::minPositiveRootEigen(vector<double> coeffs){


  // input ordering is highest order coeff first
  // converting it to lowest order coeff first for Eigen
  std::reverse(coeffs.begin(), coeffs.end());
  int deg = coeffs.size()-1;
  Eigen::VectorXd poly(deg+1);

  for (int i=0; i<=deg; i++){
    poly(i) = coeffs[i];
  }
  Eigen::PolynomialSolver<double, Eigen::Dynamic> solver;

  auto start = std::chrono::system_clock::now();
  solver.compute(poly);
  auto end = std::chrono::system_clock::now();
  const Eigen::PolynomialSolver<double, Eigen::Dynamic>::RootsType &r = solver.roots();
  double res = std::numeric_limits<double>::max();

  for (int i=0; i<r.rows(); i++){
    if ( abs(r(i,0).imag()) < 0.0001 && r(i,0).real() > 0.0 && r(i,0).real()<res ){
      res = r(i,0).real();
    }
  }


  std::chrono::duration<double> elapsed_seconds = end-start;
  total_time += elapsed_seconds.count();
  total_count++;

  if (res > 0.0){
    return res;
  }
  throw err_NoSolution();

}
*/

// highest order coeff first
// very very sketch, don't ask don't tell shhhh.....
double QuadOptimalControl::minPositiveRoot(vector<double> coeffs){
  auto start = std::chrono::system_clock::now();

  bool sign_left = peval(coeffs, 0.0) > 0;
  double step = 0.1;
  double max_range = 5.0;
  double lower_bound = 0.0;
  double upper_bound = 0.0;
  double retval;

  bool found_solution_interval = false;
  while ( upper_bound < max_range ){
    if (peval(coeffs,upper_bound) > 0 != sign_left){
      found_solution_interval = true;
      break;
    } else { 
      lower_bound = upper_bound;
    }
    upper_bound += step;
  }

  if (!found_solution_interval){ throw new err_NoSolution; }
  
  try{
    retval = bisect(coeffs, lower_bound, upper_bound);
  } catch (err_TooManyIteration) {
    throw new err_NoSolution;
  }

  auto end = std::chrono::system_clock::now();
  std::chrono::duration<double> elapsed_seconds = end-start;
  total_time += elapsed_seconds.count();
  total_count++;

  return retval;

}

// low = left = a
// high = right = b
double QuadOptimalControl::bisect( vector<double> coeffs, double low, double high){
  assert (high> low);
  bool sign_a = peval( coeffs, low) > 0;
  bool sign_b = peval( coeffs, high) > 0;

  assert (sign_a != sign_b);

  int loop_count = 0;
  while (high-low > 1e-3){
    double mid = (low+high)/2.0;
    if (peval( coeffs, mid) > 0 == sign_a) { low = mid; }
    else { high = mid; }
    loop_count++;
    if (loop_count > 100){ 
      cout << "bisect: max iter exceeded \n";
      throw new err_TooManyIteration; 
    }
  }
  return (high+low)/2.0;

}
/*
// find real roots of a polynomial
vector<double> QuadOptimalControl::roots(vector<double> coeffs){
  int deg = coeffs.size()-1;
  //cout << "deg = " << deg;
  vector<double> sol;
  if (deg == 1){
    if ( abs(coeffs.front()) < 0.00001 ){
      return sol;
    } else {
      // FIXME what if front() = 0
      sol.push_back(-coeffs[1]/coeffs.front());
      return sol;
    }
  }

  // find critical points
  auto der = derivative(coeffs);
  auto critical_points = roots(der);
  // evaluate polynomial at each critical point
  // find interval of sign change
  auto p_point = critical_points.begin();
  // true means positive
  bool last_sign = peval(coeffs, *p_point) > 0;
  double last_point = *p_point;
  p_point++;
  
  for (p_point; p_point!=critical_points.end(); p_point++){
    // find sign change
    if (peval(coeffs, *p_point) > 0 != last_sign){
      sol.push_back(falsePosition(coeffs, last_point, *p_point));
      last_sign = !last_sign;
      last_point = *p_point;
    }
  }

  // handle -inf~1st critical point 
  // true means positive
  bool sign_neg_inf = (deg % 2 == 0)?(coeffs.front() > 0):(coeffs.front() < 0);
  bool sign_pos_inf = (deg % 2 == 0)?(coeffs.front() > 0):(coeffs.front() > 0);
  // sign change in -inf~1st critical point 
  bool sol_at_left  = sign_neg_inf != peval(coeffs, critical_points.front()) > 0;
  bool sol_at_right = sign_pos_inf != peval(coeffs, critical_points.front()) > 0;

  if (sol_at_left || sol_at_right){
    auto der2 = derivative(der);
    // TODO not sure about this, finding a sufficiently large value for bound on root finding
    double a_max = *std::max_element(coeffs.begin(), coeffs.end());
    double a_min = *std::min_element(coeffs.begin(), coeffs.end());
    double a_max_abs = std::max(abs(a_max), abs(a_min));
    if (sol_at_left){
      double low_bound = - (a_max_abs / abs(coeffs.front()))*10;
      auto temp = falsePosition(coeffs, low_bound, critical_points.front());
      sol.insert(sol.begin(),temp);
    }
    if (sol_at_right){
      double high_bound = (a_max_abs / abs(coeffs.front()))*10;
      sol.push_back(falsePosition(coeffs, critical_points.back(), high_bound));
    }
  }

  return sol;

}

// find derivative of a polynomial
vector<double> QuadOptimalControl::derivative(vector<double> coeffs){
  int deg = coeffs.size()-1;
  assert (deg > 0);

  vector<double> derivative_polynomial;
  for (int i=0; i<deg; i++){
    derivative_polynomial.push_back( (deg-i)*coeffs[i] );
  }

  return derivative_polynomial;

}

double QuadOptimalControl::falsePosition(vector<double> coeffs, double low_bound, double high_bound, double err){
  // false position method:
  // draw a line between f(a) and f(b)
  // find intercept with x
  // use x as new bound
  assert (high_bound > low_bound);
  //terminal condition 1
  if (high_bound - low_bound < err){ return 0.5*(high_bound + low_bound); }

  double fa = peval( coeffs, low_bound);
  double fb = peval( coeffs, high_bound);
  assert (fa>0 != fb>0);
  double intercept = (-fa) * (high_bound - low_bound) / (fb - fa) + low_bound;
  double f_intercept = peval( coeffs, intercept);
  // terminal condition 2
  if ( abs(f_intercept) < err ) { return intercept; }

  if (f_intercept > 0 == fa > 0){
    return falsePosition(coeffs, intercept, high_bound);
  } else {
    return falsePosition(coeffs,  low_bound, intercept);
  }
  throw err_Unexcepted();
}

*/


double QuadOptimalControl::cost(double t_s, double x01, double  x02, double  x03, double  x04, double  x05, double  x06, double  x07, double  x08, double  x09, double  x11, double  x12, double  x13, double  x14, double  x15, double  x16, double  x17, double  x18, double  x19){
  return t_s + (18*pow(x01,2))/(125*pow(t_s,5)) + (24*pow(x04,2))/(625*pow(t_s,3)) + (18*pow(x02,2))/(125*pow(t_s,5)) + (9*pow(x07,2))/(5000*t_s) + (24*pow(x05,2))/(625*pow(t_s,3)) + (18*pow(x03,2))/(125*pow(t_s,5)) + (9*pow(x08,2))/(5000*t_s) + (24*pow(x06,2))/(625*pow(t_s,3)) + (9*pow(x09,2))/(5000*t_s) + (18*pow(x11,2))/(125*pow(t_s,5)) + (24*pow(x14,2))/(625*pow(t_s,3)) + (18*pow(x12,2))/(125*pow(t_s,5)) + (9*pow(x17,2))/(5000*t_s) + (24*pow(x15,2))/(625*pow(t_s,3)) + (18*pow(x13,2))/(125*pow(t_s,5)) + (9*pow(x18,2))/(5000*t_s) + (24*pow(x16,2))/(625*pow(t_s,3)) + (9*pow(x19,2))/(5000*t_s) + (18*x01*x04)/(125*pow(t_s,4)) + (3*x01*x07)/(125*pow(t_s,3)) + (18*x02*x05)/(125*pow(t_s,4)) + (9*x04*x07)/(625*pow(t_s,2)) + (3*x02*x08)/(125*pow(t_s,3)) + (18*x03*x06)/(125*pow(t_s,4)) + (9*x05*x08)/(625*pow(t_s,2)) + (3*x03*x09)/(125*pow(t_s,3)) + (9*x06*x09)/(625*pow(t_s,2)) - (36*x01*x11)/(125*pow(t_s,5)) + (18*x01*x14)/(125*pow(t_s,4)) - (18*x04*x11)/(125*pow(t_s,4)) - (36*x02*x12)/(125*pow(t_s,5)) - (3*x01*x17)/(125*pow(t_s,3)) + (42*x04*x14)/(625*pow(t_s,3)) - (3*x07*x11)/(125*pow(t_s,3)) + (18*x02*x15)/(125*pow(t_s,4)) - (18*x05*x12)/(125*pow(t_s,4)) - (36*x03*x13)/(125*pow(t_s,5)) - (6*x04*x17)/(625*pow(t_s,2)) + (6*x07*x14)/(625*pow(t_s,2)) - (3*x02*x18)/(125*pow(t_s,3)) + (42*x05*x15)/(625*pow(t_s,3)) - (3*x08*x12)/(125*pow(t_s,3)) + (18*x03*x16)/(125*pow(t_s,4)) - (18*x06*x13)/(125*pow(t_s,4)) - (3*x07*x17)/(2500*t_s) - (6*x05*x18)/(625*pow(t_s,2)) + (6*x08*x15)/(625*pow(t_s,2)) - (3*x03*x19)/(125*pow(t_s,3)) + (42*x06*x16)/(625*pow(t_s,3)) - (3*x09*x13)/(125*pow(t_s,3)) - (3*x08*x18)/(2500*t_s) - (6*x06*x19)/(625*pow(t_s,2)) + (6*x09*x16)/(625*pow(t_s,2)) - (3*x09*x19)/(2500*t_s) - (18*x11*x14)/(125*pow(t_s,4)) + (3*x11*x17)/(125*pow(t_s,3)) - (18*x12*x15)/(125*pow(t_s,4)) - (9*x14*x17)/(625*pow(t_s,2)) + (3*x12*x18)/(125*pow(t_s,3)) - (18*x13*x16)/(125*pow(t_s,4)) - (9*x15*x18)/(625*pow(t_s,2)) + (3*x13*x19)/(125*pow(t_s,3)) - (9*x16*x19)/(625*pow(t_s,2));

}

double QuadOptimalControl::costPartialFreeFinalState(double t_s, double x01, double x02, double x03, double x04, double x05, double x06, double x07, double x08, double x09, double x11, double x12, double x13){
    return (250000*pow(t_s,8) + 9000*pow(t_s,7) + 300*pow(t_s,6)*pow(x07,2) + 300*pow(t_s,6)*pow(x08,2) + 300*pow(t_s,6)*pow(x09,2) + 81*pow(t_s,6) + 1500*pow(t_s,5)*x04*x07 + 1500*pow(t_s,5)*x05*x08 + 1500*pow(t_s,5)*x06*x09 + 9*pow(t_s,5)*pow(x07,2) + 9*pow(t_s,5)*pow(x08,2) + 9*pow(t_s,5)*pow(x09,2) + 1500*pow(t_s,4)*x01*x07 + 1500*pow(t_s,4)*x02*x08 + 1500*pow(t_s,4)*x03*x09 + 2250*pow(t_s,4)*pow(x04,2) + 36*pow(t_s,4)*x04*x07 + 2250*pow(t_s,4)*pow(x05,2) + 36*pow(t_s,4)*x05*x08 + 2250*pow(t_s,4)*pow(x06,2) + 36*pow(t_s,4)*x06*x09 + (81*pow(t_s,4)*pow(x07,2))/1000 - 1500*pow(t_s,4)*x07*x11 + (81*pow(t_s,4)*pow(x08,2))/1000 - 1500*pow(t_s,4)*x08*x12 + (81*pow(t_s,4)*pow(x09,2))/1000 - 1500*pow(t_s,4)*x09*x13 + 4500*pow(t_s,3)*x01*x04 + 36*pow(t_s,3)*x01*x07 + 4500*pow(t_s,3)*x02*x05 + 36*pow(t_s,3)*x02*x08 + 4500*pow(t_s,3)*x03*x06 + 36*pow(t_s,3)*x03*x09 + 36*pow(t_s,3)*pow(x04,2) + (81*pow(t_s,3)*x04*x07)/250 - 4500*pow(t_s,3)*x04*x11 + 36*pow(t_s,3)*pow(x05,2) + (81*pow(t_s,3)*x05*x08)/250 - 4500*pow(t_s,3)*x05*x12 + 36*pow(t_s,3)*pow(x06,2) + (81*pow(t_s,3)*x06*x09)/250 - 4500*pow(t_s,3)*x06*x13 - 36*pow(t_s,3)*x07*x11 - 36*pow(t_s,3)*x08*x12 - 36*pow(t_s,3)*x09*x13 + 2250*pow(t_s,2)*pow(x01,2) + 72*pow(t_s,2)*x01*x04 + (81*pow(t_s,2)*x01*x07)/250 - 4500*pow(t_s,2)*x01*x11 + 2250*pow(t_s,2)*pow(x02,2) + 72*pow(t_s,2)*x02*x05 + (81*pow(t_s,2)*x02*x08)/250 - 4500*pow(t_s,2)*x02*x12 + 2250*pow(t_s,2)*pow(x03,2) + 72*pow(t_s,2)*x03*x06 + (81*pow(t_s,2)*x03*x09)/250 - 4500*pow(t_s,2)*x03*x13 + (81*pow(t_s,2)*pow(x04,2))/250 - 72*pow(t_s,2)*x04*x11 + (81*pow(t_s,2)*pow(x05,2))/250 - 72*pow(t_s,2)*x05*x12 + (81*pow(t_s,2)*pow(x06,2))/250 - 72*pow(t_s,2)*x06*x13 - (81*pow(t_s,2)*x07*x11)/250 - (81*pow(t_s,2)*x08*x12)/250 - (81*pow(t_s,2)*x09*x13)/250 + 2250*pow(t_s,2)*pow(x11,2) + 2250*pow(t_s,2)*pow(x12,2) + 2250*pow(t_s,2)*pow(x13,2) + 36*t_s*pow(x01,2) + (81*t_s*x01*x04)/125 - 72*t_s*x01*x11 + 36*t_s*pow(x02,2) + (81*t_s*x02*x05)/125 - 72*t_s*x02*x12 + 36*t_s*pow(x03,2) + (81*t_s*x03*x06)/125 - 72*t_s*x03*x13 - (81*t_s*x04*x11)/125 - (81*t_s*x05*x12)/125 - (81*t_s*x06*x13)/125 + 36*t_s*pow(x11,2) + 36*t_s*pow(x12,2) + 36*t_s*pow(x13,2) + (81*pow(x01,2))/250 - (81*x01*x11)/125 + (81*pow(x02,2))/250 - (81*x02*x12)/125 + (81*pow(x03,2))/250 - (81*x03*x13)/125 + (81*pow(x11,2))/250 + (81*pow(x12,2))/250 + (81*pow(x13,2))/250)/(pow(t_s,5)*pow((500*t_s + 9),2));
}


double* QuadOptimalControl::interiorPosition(double t_s, double x01, double  x02, double  x03, double  x04, double  x05, double  x06, double  x07, double  x08, double  x09, double  x11, double  x12, double  x13, double  x14, double  x15, double  x16, double  x17, double  x18, double  x19){
  // prepare time vector
  double t_step = t_s/(interior_point_count+1);
  double t = t_step;
  for (int i=0;i<interior_point_count;i++){
    *(buffer_interior_pos + i * 3 + 0) = x11 + x14*(t - t_s) + (x17*pow( (t - t_s),2))/2 - (pow( (t - t_s),5)*(12*x01 - 12*x11 + 6*t_s*x04 + 6*t_s*x14 + pow(t_s,2)*x07 - pow(t_s,2)*x17))/(2*pow(t_s,5)) - (pow( (t - t_s),3)*(20*x01 - 20*x11 + 8*t_s*x04 + 12*t_s*x14 + pow(t_s,2)*x07 - 3*pow(t_s,2)*x17))/(2*pow(t_s,3)) - (pow( (t - t_s),4)*(30*x01 - 30*x11 + 14*t_s*x04 + 16*t_s*x14 + 2*pow(t_s,2)*x07 - 3*pow(t_s,2)*x17))/(2*pow(t_s,4));
    *(buffer_interior_pos + i * 3 + 1) = x12 + x15*(t - t_s) + (x18*pow( (t - t_s),2))/2 - (pow( (t - t_s),5)*(12*x02 - 12*x12 + 6*t_s*x05 + 6*t_s*x15 + pow(t_s,2)*x08 - pow(t_s,2)*x18))/(2*pow(t_s,5)) - (pow( (t - t_s),3)*(20*x02 - 20*x12 + 8*t_s*x05 + 12*t_s*x15 + pow(t_s,2)*x08 - 3*pow(t_s,2)*x18))/(2*pow(t_s,3)) - (pow( (t - t_s),4)*(30*x02 - 30*x12 + 14*t_s*x05 + 16*t_s*x15 + 2*pow(t_s,2)*x08 - 3*pow(t_s,2)*x18))/(2*pow(t_s,4));
    *(buffer_interior_pos + i * 3 + 2) = x13 + x16*(t - t_s) + (x19*pow( (t - t_s),2))/2 - (pow( (t - t_s),5)*(12*x03 - 12*x13 + 6*t_s*x06 + 6*t_s*x16 + pow(t_s,2)*x09 - pow(t_s,2)*x19))/(2*pow(t_s,5)) - (pow( (t - t_s),3)*(20*x03 - 20*x13 + 8*t_s*x06 + 12*t_s*x16 + pow(t_s,2)*x09 - 3*pow(t_s,2)*x19))/(2*pow(t_s,3)) - (pow( (t - t_s),4)*(30*x03 - 30*x13 + 14*t_s*x06 + 16*t_s*x16 + 2*pow(t_s,2)*x09 - 3*pow(t_s,2)*x19))/(2*pow(t_s,4));

    t += t_step;
  }
  return buffer_interior_pos;
}

double* QuadOptimalControl::interiorPositionPartialFinalState(double t_s, double x01, double x02, double x03, double x04, double x05, double x06, double x07, double x08, double x09, double x11, double x12, double x13){
  // prepare time vector
  double t_step = t_s/(interior_point_count+1);
  double t = t_step;
  for (int i=0;i<interior_point_count;i++){
    *(buffer_interior_pos + i * 3 + 0) = x01 + t*x04 + (pow(t,2)*x07)/2 + (625*pow(t,3)*x01)/(500*pow(t_s,3) + 9*pow(t_s,2)) - (15*pow(t,3)*x01)/(500*pow(t_s,4) + 9*pow(t_s,3)) + (15*pow(t,4)*x01)/(2*(500*pow(t_s,5) + 9*pow(t_s,4))) - (3*pow(t,5)*x01)/(2*(500*pow(t_s,6) + 9*pow(t_s,5))) - (625*pow(t,3)*x11)/(500*pow(t_s,3) + 9*pow(t_s,2)) + (15*pow(t,3)*x11)/(500*pow(t_s,4) + 9*pow(t_s,3)) - (15*pow(t,4)*x11)/(2*(500*pow(t_s,5) + 9*pow(t_s,4))) + (3*pow(t,5)*x11)/(2*(500*pow(t_s,6) + 9*pow(t_s,5))) - (1875*pow(t,3)*t_s*x01)/(500*pow(t_s,4) + 9*pow(t_s,3)) + (625*pow(t,3)*t_s*x04)/(500*pow(t_s,3) + 9*pow(t_s,2)) - (15*pow(t,3)*t_s*x04)/(500*pow(t_s,4) + 9*pow(t_s,3)) + (1875*pow(t,4)*t_s*x01)/(2*(500*pow(t_s,5) + 9*pow(t_s,4))) + (15*pow(t,4)*t_s*x04)/(2*(500*pow(t_s,5) + 9*pow(t_s,4))) - (375*pow(t,5)*t_s*x01)/(2*(500*pow(t_s,6) + 9*pow(t_s,5))) - (3*pow(t,5)*t_s*x04)/(2*(500*pow(t_s,6) + 9*pow(t_s,5))) + (1875*pow(t,3)*t_s*x11)/(500*pow(t_s,4) + 9*pow(t_s,3)) - (1875*pow(t,4)*t_s*x11)/(2*(500*pow(t_s,5) + 9*pow(t_s,4))) + (375*pow(t,5)*t_s*x11)/(2*(500*pow(t_s,6) + 9*pow(t_s,5))) - (1875*pow(t,3)*pow(t_s,2)*x04)/(500*pow(t_s,4) + 9*pow(t_s,3)) + (125*pow(t,3)*pow(t_s,2)*x07)/(500*pow(t_s,3) + 9*pow(t_s,2)) - (15*pow(t,3)*pow(t_s,2)*x07)/(2*(500*pow(t_s,4) + 9*pow(t_s,3))) + (1875*pow(t,4)*pow(t_s,2)*x04)/(2*(500*pow(t_s,5) + 9*pow(t_s,4))) - (625*pow(t,3)*pow(t_s,3)*x07)/(500*pow(t_s,4) + 9*pow(t_s,3)) + (15*pow(t,4)*pow(t_s,2)*x07)/(4*(500*pow(t_s,5) + 9*pow(t_s,4))) - (375*pow(t,5)*pow(t_s,2)*x04)/(2*(500*pow(t_s,6) + 9*pow(t_s,5))) + (625*pow(t,4)*pow(t_s,3)*x07)/(2*(500*pow(t_s,5) + 9*pow(t_s,4))) - (3*pow(t,5)*pow(t_s,2)*x07)/(4*(500*pow(t_s,6) + 9*pow(t_s,5))) - (125*pow(t,5)*pow(t_s,3)*x07)/(2*(500*pow(t_s,6) + 9*pow(t_s,5)));
    *(buffer_interior_pos + i * 3 + 1) =  x02 + t*x05 + (pow(t,2)*x08)/2 + (625*pow(t,3)*x02)/(500*pow(t_s,3) + 9*pow(t_s,2)) - (15*pow(t,3)*x02)/(500*pow(t_s,4) + 9*pow(t_s,3)) + (15*pow(t,4)*x02)/(2*(500*pow(t_s,5) + 9*pow(t_s,4))) - (3*pow(t,5)*x02)/(2*(500*pow(t_s,6) + 9*pow(t_s,5))) - (625*pow(t,3)*x12)/(500*pow(t_s,3) + 9*pow(t_s,2)) + (15*pow(t,3)*x12)/(500*pow(t_s,4) + 9*pow(t_s,3)) - (15*pow(t,4)*x12)/(2*(500*pow(t_s,5) + 9*pow(t_s,4))) + (3*pow(t,5)*x12)/(2*(500*pow(t_s,6) + 9*pow(t_s,5))) - (1875*pow(t,3)*t_s*x02)/(500*pow(t_s,4) + 9*pow(t_s,3)) + (625*pow(t,3)*t_s*x05)/(500*pow(t_s,3) + 9*pow(t_s,2)) - (15*pow(t,3)*t_s*x05)/(500*pow(t_s,4) + 9*pow(t_s,3)) + (1875*pow(t,4)*t_s*x02)/(2*(500*pow(t_s,5) + 9*pow(t_s,4))) + (15*pow(t,4)*t_s*x05)/(2*(500*pow(t_s,5) + 9*pow(t_s,4))) - (375*pow(t,5)*t_s*x02)/(2*(500*pow(t_s,6) + 9*pow(t_s,5))) - (3*pow(t,5)*t_s*x05)/(2*(500*pow(t_s,6) + 9*pow(t_s,5))) + (1875*pow(t,3)*t_s*x12)/(500*pow(t_s,4) + 9*pow(t_s,3)) - (1875*pow(t,4)*t_s*x12)/(2*(500*pow(t_s,5) + 9*pow(t_s,4))) + (375*pow(t,5)*t_s*x12)/(2*(500*pow(t_s,6) + 9*pow(t_s,5))) - (1875*pow(t,3)*pow(t_s,2)*x05)/(500*pow(t_s,4) + 9*pow(t_s,3)) + (125*pow(t,3)*pow(t_s,2)*x08)/(500*pow(t_s,3) + 9*pow(t_s,2)) - (15*pow(t,3)*pow(t_s,2)*x08)/(2*(500*pow(t_s,4) + 9*pow(t_s,3))) + (1875*pow(t,4)*pow(t_s,2)*x05)/(2*(500*pow(t_s,5) + 9*pow(t_s,4))) - (625*pow(t,3)*pow(t_s,3)*x08)/(500*pow(t_s,4) + 9*pow(t_s,3)) + (15*pow(t,4)*pow(t_s,2)*x08)/(4*(500*pow(t_s,5) + 9*pow(t_s,4))) - (375*pow(t,5)*pow(t_s,2)*x05)/(2*(500*pow(t_s,6) + 9*pow(t_s,5))) + (625*pow(t,4)*pow(t_s,3)*x08)/(2*(500*pow(t_s,5) + 9*pow(t_s,4))) - (3*pow(t,5)*pow(t_s,2)*x08)/(4*(500*pow(t_s,6) + 9*pow(t_s,5))) - (125*pow(t,5)*pow(t_s,3)*x08)/(2*(500*pow(t_s,6) + 9*pow(t_s,5)));
    *(buffer_interior_pos + i * 3 + 2) =  x03 + t*x06 + (pow(t,2)*x09)/2 + (625*pow(t,3)*x03)/(500*pow(t_s,3) + 9*pow(t_s,2)) - (15*pow(t,3)*x03)/(500*pow(t_s,4) + 9*pow(t_s,3)) + (15*pow(t,4)*x03)/(2*(500*pow(t_s,5) + 9*pow(t_s,4))) - (3*pow(t,5)*x03)/(2*(500*pow(t_s,6) + 9*pow(t_s,5))) - (625*pow(t,3)*x13)/(500*pow(t_s,3) + 9*pow(t_s,2)) + (15*pow(t,3)*x13)/(500*pow(t_s,4) + 9*pow(t_s,3)) - (15*pow(t,4)*x13)/(2*(500*pow(t_s,5) + 9*pow(t_s,4))) + (3*pow(t,5)*x13)/(2*(500*pow(t_s,6) + 9*pow(t_s,5))) - (1875*pow(t,3)*t_s*x03)/(500*pow(t_s,4) + 9*pow(t_s,3)) + (625*pow(t,3)*t_s*x06)/(500*pow(t_s,3) + 9*pow(t_s,2)) - (15*pow(t,3)*t_s*x06)/(500*pow(t_s,4) + 9*pow(t_s,3)) + (1875*pow(t,4)*t_s*x03)/(2*(500*pow(t_s,5) + 9*pow(t_s,4))) + (15*pow(t,4)*t_s*x06)/(2*(500*pow(t_s,5) + 9*pow(t_s,4))) - (375*pow(t,5)*t_s*x03)/(2*(500*pow(t_s,6) + 9*pow(t_s,5))) - (3*pow(t,5)*t_s*x06)/(2*(500*pow(t_s,6) + 9*pow(t_s,5))) + (1875*pow(t,3)*t_s*x13)/(500*pow(t_s,4) + 9*pow(t_s,3)) - (1875*pow(t,4)*t_s*x13)/(2*(500*pow(t_s,5) + 9*pow(t_s,4))) + (375*pow(t,5)*t_s*x13)/(2*(500*pow(t_s,6) + 9*pow(t_s,5))) - (1875*pow(t,3)*pow(t_s,2)*x06)/(500*pow(t_s,4) + 9*pow(t_s,3)) + (125*pow(t,3)*pow(t_s,2)*x09)/(500*pow(t_s,3) + 9*pow(t_s,2)) - (15*pow(t,3)*pow(t_s,2)*x09)/(2*(500*pow(t_s,4) + 9*pow(t_s,3))) + (1875*pow(t,4)*pow(t_s,2)*x06)/(2*(500*pow(t_s,5) + 9*pow(t_s,4))) - (625*pow(t,3)*pow(t_s,3)*x09)/(500*pow(t_s,4) + 9*pow(t_s,3)) + (15*pow(t,4)*pow(t_s,2)*x09)/(4*(500*pow(t_s,5) + 9*pow(t_s,4))) - (375*pow(t,5)*pow(t_s,2)*x06)/(2*(500*pow(t_s,6) + 9*pow(t_s,5))) + (625*pow(t,4)*pow(t_s,3)*x09)/(2*(500*pow(t_s,5) + 9*pow(t_s,4))) - (3*pow(t,5)*pow(t_s,2)*x09)/(4*(500*pow(t_s,6) + 9*pow(t_s,5))) - (125*pow(t,5)*pow(t_s,3)*x09)/(2*(500*pow(t_s,6) + 9*pow(t_s,5)));

    t += t_step;
  }
  return buffer_interior_pos;
}


double* QuadOptimalControl::interiorStatePartialFinalState(double t_s, double x01, double x02, double x03, double x04, double x05, double x06, double x07, double x08, double x09, double x11, double x12, double x13){
  // prepare time vector
  double t_step = t_s/(interior_point_count+1);
  double t = t_step;
  for (int i=0;i<interior_point_count;i++){
    *(buffer_interior_states + i * 9 + 0) = x01 + t*x04 + (pow(t,2)*x07)/2 + (625*pow(t,3)*x01)/(500*pow(t_s,3) + 9*pow(t_s,2)) - (15*pow(t,3)*x01)/(500*pow(t_s,4) + 9*pow(t_s,3)) + (15*pow(t,4)*x01)/(2*(500*pow(t_s,5) + 9*pow(t_s,4))) - (3*pow(t,5)*x01)/(2*(500*pow(t_s,6) + 9*pow(t_s,5))) - (625*pow(t,3)*x11)/(500*pow(t_s,3) + 9*pow(t_s,2)) + (15*pow(t,3)*x11)/(500*pow(t_s,4) + 9*pow(t_s,3)) - (15*pow(t,4)*x11)/(2*(500*pow(t_s,5) + 9*pow(t_s,4))) + (3*pow(t,5)*x11)/(2*(500*pow(t_s,6) + 9*pow(t_s,5))) - (1875*pow(t,3)*t_s*x01)/(500*pow(t_s,4) + 9*pow(t_s,3)) + (625*pow(t,3)*t_s*x04)/(500*pow(t_s,3) + 9*pow(t_s,2)) - (15*pow(t,3)*t_s*x04)/(500*pow(t_s,4) + 9*pow(t_s,3)) + (1875*pow(t,4)*t_s*x01)/(2*(500*pow(t_s,5) + 9*pow(t_s,4))) + (15*pow(t,4)*t_s*x04)/(2*(500*pow(t_s,5) + 9*pow(t_s,4))) - (375*pow(t,5)*t_s*x01)/(2*(500*pow(t_s,6) + 9*pow(t_s,5))) - (3*pow(t,5)*t_s*x04)/(2*(500*pow(t_s,6) + 9*pow(t_s,5))) + (1875*pow(t,3)*t_s*x11)/(500*pow(t_s,4) + 9*pow(t_s,3)) - (1875*pow(t,4)*t_s*x11)/(2*(500*pow(t_s,5) + 9*pow(t_s,4))) + (375*pow(t,5)*t_s*x11)/(2*(500*pow(t_s,6) + 9*pow(t_s,5))) - (1875*pow(t,3)*pow(t_s,2)*x04)/(500*pow(t_s,4) + 9*pow(t_s,3)) + (125*pow(t,3)*pow(t_s,2)*x07)/(500*pow(t_s,3) + 9*pow(t_s,2)) - (15*pow(t,3)*pow(t_s,2)*x07)/(2*(500*pow(t_s,4) + 9*pow(t_s,3))) + (1875*pow(t,4)*pow(t_s,2)*x04)/(2*(500*pow(t_s,5) + 9*pow(t_s,4))) - (625*pow(t,3)*pow(t_s,3)*x07)/(500*pow(t_s,4) + 9*pow(t_s,3)) + (15*pow(t,4)*pow(t_s,2)*x07)/(4*(500*pow(t_s,5) + 9*pow(t_s,4))) - (375*pow(t,5)*pow(t_s,2)*x04)/(2*(500*pow(t_s,6) + 9*pow(t_s,5))) + (625*pow(t,4)*pow(t_s,3)*x07)/(2*(500*pow(t_s,5) + 9*pow(t_s,4))) - (3*pow(t,5)*pow(t_s,2)*x07)/(4*(500*pow(t_s,6) + 9*pow(t_s,5))) - (125*pow(t,5)*pow(t_s,3)*x07)/(2*(500*pow(t_s,6) + 9*pow(t_s,5)));
    *(buffer_interior_states + i * 9 + 1) =  x02 + t*x05 + (pow(t,2)*x08)/2 + (625*pow(t,3)*x02)/(500*pow(t_s,3) + 9*pow(t_s,2)) - (15*pow(t,3)*x02)/(500*pow(t_s,4) + 9*pow(t_s,3)) + (15*pow(t,4)*x02)/(2*(500*pow(t_s,5) + 9*pow(t_s,4))) - (3*pow(t,5)*x02)/(2*(500*pow(t_s,6) + 9*pow(t_s,5))) - (625*pow(t,3)*x12)/(500*pow(t_s,3) + 9*pow(t_s,2)) + (15*pow(t,3)*x12)/(500*pow(t_s,4) + 9*pow(t_s,3)) - (15*pow(t,4)*x12)/(2*(500*pow(t_s,5) + 9*pow(t_s,4))) + (3*pow(t,5)*x12)/(2*(500*pow(t_s,6) + 9*pow(t_s,5))) - (1875*pow(t,3)*t_s*x02)/(500*pow(t_s,4) + 9*pow(t_s,3)) + (625*pow(t,3)*t_s*x05)/(500*pow(t_s,3) + 9*pow(t_s,2)) - (15*pow(t,3)*t_s*x05)/(500*pow(t_s,4) + 9*pow(t_s,3)) + (1875*pow(t,4)*t_s*x02)/(2*(500*pow(t_s,5) + 9*pow(t_s,4))) + (15*pow(t,4)*t_s*x05)/(2*(500*pow(t_s,5) + 9*pow(t_s,4))) - (375*pow(t,5)*t_s*x02)/(2*(500*pow(t_s,6) + 9*pow(t_s,5))) - (3*pow(t,5)*t_s*x05)/(2*(500*pow(t_s,6) + 9*pow(t_s,5))) + (1875*pow(t,3)*t_s*x12)/(500*pow(t_s,4) + 9*pow(t_s,3)) - (1875*pow(t,4)*t_s*x12)/(2*(500*pow(t_s,5) + 9*pow(t_s,4))) + (375*pow(t,5)*t_s*x12)/(2*(500*pow(t_s,6) + 9*pow(t_s,5))) - (1875*pow(t,3)*pow(t_s,2)*x05)/(500*pow(t_s,4) + 9*pow(t_s,3)) + (125*pow(t,3)*pow(t_s,2)*x08)/(500*pow(t_s,3) + 9*pow(t_s,2)) - (15*pow(t,3)*pow(t_s,2)*x08)/(2*(500*pow(t_s,4) + 9*pow(t_s,3))) + (1875*pow(t,4)*pow(t_s,2)*x05)/(2*(500*pow(t_s,5) + 9*pow(t_s,4))) - (625*pow(t,3)*pow(t_s,3)*x08)/(500*pow(t_s,4) + 9*pow(t_s,3)) + (15*pow(t,4)*pow(t_s,2)*x08)/(4*(500*pow(t_s,5) + 9*pow(t_s,4))) - (375*pow(t,5)*pow(t_s,2)*x05)/(2*(500*pow(t_s,6) + 9*pow(t_s,5))) + (625*pow(t,4)*pow(t_s,3)*x08)/(2*(500*pow(t_s,5) + 9*pow(t_s,4))) - (3*pow(t,5)*pow(t_s,2)*x08)/(4*(500*pow(t_s,6) + 9*pow(t_s,5))) - (125*pow(t,5)*pow(t_s,3)*x08)/(2*(500*pow(t_s,6) + 9*pow(t_s,5)));
    *(buffer_interior_states + i * 9 + 2) =   x03 + t*x06 + (pow(t,2)*x09)/2 + (625*pow(t,3)*x03)/(500*pow(t_s,3) + 9*pow(t_s,2)) - (15*pow(t,3)*x03)/(500*pow(t_s,4) + 9*pow(t_s,3)) + (15*pow(t,4)*x03)/(2*(500*pow(t_s,5) + 9*pow(t_s,4))) - (3*pow(t,5)*x03)/(2*(500*pow(t_s,6) + 9*pow(t_s,5))) - (625*pow(t,3)*x13)/(500*pow(t_s,3) + 9*pow(t_s,2)) + (15*pow(t,3)*x13)/(500*pow(t_s,4) + 9*pow(t_s,3)) - (15*pow(t,4)*x13)/(2*(500*pow(t_s,5) + 9*pow(t_s,4))) + (3*pow(t,5)*x13)/(2*(500*pow(t_s,6) + 9*pow(t_s,5))) - (1875*pow(t,3)*t_s*x03)/(500*pow(t_s,4) + 9*pow(t_s,3)) + (625*pow(t,3)*t_s*x06)/(500*pow(t_s,3) + 9*pow(t_s,2)) - (15*pow(t,3)*t_s*x06)/(500*pow(t_s,4) + 9*pow(t_s,3)) + (1875*pow(t,4)*t_s*x03)/(2*(500*pow(t_s,5) + 9*pow(t_s,4))) + (15*pow(t,4)*t_s*x06)/(2*(500*pow(t_s,5) + 9*pow(t_s,4))) - (375*pow(t,5)*t_s*x03)/(2*(500*pow(t_s,6) + 9*pow(t_s,5))) - (3*pow(t,5)*t_s*x06)/(2*(500*pow(t_s,6) + 9*pow(t_s,5))) + (1875*pow(t,3)*t_s*x13)/(500*pow(t_s,4) + 9*pow(t_s,3)) - (1875*pow(t,4)*t_s*x13)/(2*(500*pow(t_s,5) + 9*pow(t_s,4))) + (375*pow(t,5)*t_s*x13)/(2*(500*pow(t_s,6) + 9*pow(t_s,5))) - (1875*pow(t,3)*pow(t_s,2)*x06)/(500*pow(t_s,4) + 9*pow(t_s,3)) + (125*pow(t,3)*pow(t_s,2)*x09)/(500*pow(t_s,3) + 9*pow(t_s,2)) - (15*pow(t,3)*pow(t_s,2)*x09)/(2*(500*pow(t_s,4) + 9*pow(t_s,3))) + (1875*pow(t,4)*pow(t_s,2)*x06)/(2*(500*pow(t_s,5) + 9*pow(t_s,4))) - (625*pow(t,3)*pow(t_s,3)*x09)/(500*pow(t_s,4) + 9*pow(t_s,3)) + (15*pow(t,4)*pow(t_s,2)*x09)/(4*(500*pow(t_s,5) + 9*pow(t_s,4))) - (375*pow(t,5)*pow(t_s,2)*x06)/(2*(500*pow(t_s,6) + 9*pow(t_s,5))) + (625*pow(t,4)*pow(t_s,3)*x09)/(2*(500*pow(t_s,5) + 9*pow(t_s,4))) - (3*pow(t,5)*pow(t_s,2)*x09)/(4*(500*pow(t_s,6) + 9*pow(t_s,5))) - (125*pow(t,5)*pow(t_s,3)*x09)/(2*(500*pow(t_s,6) + 9*pow(t_s,5)));
    *(buffer_interior_states + i * 9 + 3) = x04 + t*x07 - ((5*pow(t,2)*(6*pow(t,2)*x01 - 6*pow(t,2)*x11))/4 + (5*pow(t,2)*pow(t_s,3)*(3000*x01 + 36*x04 - 3000*x11 - 3000*t*x04 - 12*t*x07 + 250*pow(t,2)*x07))/4 + (5*pow(t,2)*pow(t_s,2)*(36*x01 - 36*x11 - 3000*t*x01 - 24*t*x04 + 3000*t*x11 + 750*pow(t,2)*x04 + 3*pow(t,2)*x07))/4 + (5*pow(t,2)*pow(t_s,4)*(3000*x04 + 18*x07 - 1000*t*x07))/4 + 1500*pow(t,2)*pow(t_s,5)*x07 + (5*pow(t,2)*t_s*(24*t*x11 - 24*t*x01 + 750*pow(t,2)*x01 + 6*pow(t,2)*x04 - 750*pow(t,2)*x11))/4)/(pow(t_s,5)*(500*t_s + 9));
    *(buffer_interior_states + i * 9 + 4) = x05 + t*x08 - ((5*pow(t,2)*(6*pow(t,2)*x02 - 6*pow(t,2)*x12))/4 + (5*pow(t,2)*pow(t_s,3)*(3000*x02 + 36*x05 - 3000*x12 - 3000*t*x05 - 12*t*x08 + 250*pow(t,2)*x08))/4 + (5*pow(t,2)*pow(t_s,2)*(36*x02 - 36*x12 - 3000*t*x02 - 24*t*x05 + 3000*t*x12 + 750*pow(t,2)*x05 + 3*pow(t,2)*x08))/4 + (5*pow(t,2)*pow(t_s,4)*(3000*x05 + 18*x08 - 1000*t*x08))/4 + 1500*pow(t,2)*pow(t_s,5)*x08 + (5*pow(t,2)*t_s*(24*t*x12 - 24*t*x02 + 750*pow(t,2)*x02 + 6*pow(t,2)*x05 - 750*pow(t,2)*x12))/4)/(pow(t_s,5)*(500*t_s + 9));
    *(buffer_interior_states + i * 9 + 5) = x06 + t*x09 - ((5*pow(t,2)*(6*pow(t,2)*x03 - 6*pow(t,2)*x13))/4 + (5*pow(t,2)*pow(t_s,3)*(3000*x03 + 36*x06 - 3000*x13 - 3000*t*x06 - 12*t*x09 + 250*pow(t,2)*x09))/4 + (5*pow(t,2)*pow(t_s,2)*(36*x03 - 36*x13 - 3000*t*x03 - 24*t*x06 + 3000*t*x13 + 750*pow(t,2)*x06 + 3*pow(t,2)*x09))/4 + (5*pow(t,2)*pow(t_s,4)*(3000*x06 + 18*x09 - 1000*t*x09))/4 + 1500*pow(t,2)*pow(t_s,5)*x09 + (5*pow(t,2)*t_s*(24*t*x13 - 24*t*x03 + 750*pow(t,2)*x03 + 6*pow(t,2)*x06 - 750*pow(t,2)*x13))/4)/(pow(t_s,5)*(500*t_s + 9));
    *(buffer_interior_states + i * 9 + 6) = x07 - (3000*t*x07*pow(t_s,5) + 5*t*(1500*x04 + 9*x07 - 750*t*x07)*pow(t_s,4) + 5*t*(1500*x01 + 18*x04 - 1500*x11 - 2250*t*x04 - 9*t*x07 + 250*pow(t,2)*x07)*pow(t_s,3) + 5*t*(18*x01 - 18*x11 - 2250*t*x01 - 18*t*x04 + 2250*t*x11 + 750*pow(t,2)*x04 + 3*pow(t,2)*x07)*pow(t_s,2) + 5*t*(18*t*x11 - 18*t*x01 + 750*pow(t,2)*x01 + 6*pow(t,2)*x04 - 750*pow(t,2)*x11)*t_s + 5*t*(6*pow(t,2)*x01 - 6*pow(t,2)*x11))/(pow(t_s,5)*(500*t_s + 9));
    *(buffer_interior_states + i * 9 + 7) = x08 - (3000*t*x08*pow(t_s,5) + 5*t*(1500*x05 + 9*x08 - 750*t*x08)*pow(t_s,4) + 5*t*(1500*x02 + 18*x05 - 1500*x12 - 2250*t*x05 - 9*t*x08 + 250*pow(t,2)*x08)*pow(t_s,3) + 5*t*(18*x02 - 18*x12 - 2250*t*x02 - 18*t*x05 + 2250*t*x12 + 750*pow(t,2)*x05 + 3*pow(t,2)*x08)*pow(t_s,2) + 5*t*(18*t*x12 - 18*t*x02 + 750*pow(t,2)*x02 + 6*pow(t,2)*x05 - 750*pow(t,2)*x12)*t_s + 5*t*(6*pow(t,2)*x02 - 6*pow(t,2)*x12))/(pow(t_s,5)*(500*t_s + 9));
    *(buffer_interior_states + i * 9 + 8) = x09 - (3000*t*x09*pow(t_s,5) + 5*t*(1500*x06 + 9*x09 - 750*t*x09)*pow(t_s,4) + 5*t*(1500*x03 + 18*x06 - 1500*x13 - 2250*t*x06 - 9*t*x09 + 250*pow(t,2)*x09)*pow(t_s,3) + 5*t*(18*x03 - 18*x13 - 2250*t*x03 - 18*t*x06 + 2250*t*x13 + 750*pow(t,2)*x06 + 3*pow(t,2)*x09)*pow(t_s,2) + 5*t*(18*t*x13 - 18*t*x03 + 750*pow(t,2)*x03 + 6*pow(t,2)*x06 - 750*pow(t,2)*x13)*t_s + 5*t*(6*pow(t,2)*x03 - 6*pow(t,2)*x13))/(pow(t_s,5)*(500*t_s + 9));

    t += t_step;
  }
  return buffer_interior_pos;

}


void QuadOptimalControl::setFullStatePartialFinalState(double t_s, Node& node_i, Node& node_f){
  double x01 = node_i.x;
  double x02 = node_i.y;
  double x03 = node_i.z;
  double x04 = node_i.vx;
  double x05 = node_i.vy;
  double x06 = node_i.vz;
  double x07 = node_i.ax;
  double x08 = node_i.ay;
  double x09 = node_i.az;
  double x11 = node_f.x;
  double x12 = node_f.y;
  double x13 = node_f.z;

  double t = t_s;
  node_f.vx = x04 + t*x07 - ((5*pow(t,2)*(6*pow(t,2)*x01 - 6*pow(t,2)*x11))/4 + (5*pow(t,2)*pow(t_s,3)*(3000*x01 + 36*x04 - 3000*x11 - 3000*t*x04 - 12*t*x07 + 250*pow(t,2)*x07))/4 + (5*pow(t,2)*pow(t_s,2)*(36*x01 - 36*x11 - 3000*t*x01 - 24*t*x04 + 3000*t*x11 + 750*pow(t,2)*x04 + 3*pow(t,2)*x07))/4 + (5*pow(t,2)*pow(t_s,4)*(3000*x04 + 18*x07 - 1000*t*x07))/4 + 1500*pow(t,2)*pow(t_s,5)*x07 + (5*pow(t,2)*t_s*(24*t*x11 - 24*t*x01 + 750*pow(t,2)*x01 + 6*pow(t,2)*x04 - 750*pow(t,2)*x11))/4)/(pow(t_s,5)*(500*t_s + 9));
  node_f.vy = x05 + t*x08 - ((5*pow(t,2)*(6*pow(t,2)*x02 - 6*pow(t,2)*x12))/4 + (5*pow(t,2)*pow(t_s,3)*(3000*x02 + 36*x05 - 3000*x12 - 3000*t*x05 - 12*t*x08 + 250*pow(t,2)*x08))/4 + (5*pow(t,2)*pow(t_s,2)*(36*x02 - 36*x12 - 3000*t*x02 - 24*t*x05 + 3000*t*x12 + 750*pow(t,2)*x05 + 3*pow(t,2)*x08))/4 + (5*pow(t,2)*pow(t_s,4)*(3000*x05 + 18*x08 - 1000*t*x08))/4 + 1500*pow(t,2)*pow(t_s,5)*x08 + (5*pow(t,2)*t_s*(24*t*x12 - 24*t*x02 + 750*pow(t,2)*x02 + 6*pow(t,2)*x05 - 750*pow(t,2)*x12))/4)/(pow(t_s,5)*(500*t_s + 9));
  node_f.vz = x06 + t*x09 - ((5*pow(t,2)*(6*pow(t,2)*x03 - 6*pow(t,2)*x13))/4 + (5*pow(t,2)*pow(t_s,3)*(3000*x03 + 36*x06 - 3000*x13 - 3000*t*x06 - 12*t*x09 + 250*pow(t,2)*x09))/4 + (5*pow(t,2)*pow(t_s,2)*(36*x03 - 36*x13 - 3000*t*x03 - 24*t*x06 + 3000*t*x13 + 750*pow(t,2)*x06 + 3*pow(t,2)*x09))/4 + (5*pow(t,2)*pow(t_s,4)*(3000*x06 + 18*x09 - 1000*t*x09))/4 + 1500*pow(t,2)*pow(t_s,5)*x09 + (5*pow(t,2)*t_s*(24*t*x13 - 24*t*x03 + 750*pow(t,2)*x03 + 6*pow(t,2)*x06 - 750*pow(t,2)*x13))/4)/(pow(t_s,5)*(500*t_s + 9));

  node_f.ax = x07 - (3000*t*x07*pow(t_s,5) + 5*t*(1500*x04 + 9*x07 - 750*t*x07)*pow(t_s,4) + 5*t*(1500*x01 + 18*x04 - 1500*x11 - 2250*t*x04 - 9*t*x07 + 250*pow(t,2)*x07)*pow(t_s,3) + 5*t*(18*x01 - 18*x11 - 2250*t*x01 - 18*t*x04 + 2250*t*x11 + 750*pow(t,2)*x04 + 3*pow(t,2)*x07)*pow(t_s,2) + 5*t*(18*t*x11 - 18*t*x01 + 750*pow(t,2)*x01 + 6*pow(t,2)*x04 - 750*pow(t,2)*x11)*t_s + 5*t*(6*pow(t,2)*x01 - 6*pow(t,2)*x11))/(pow(t_s,5)*(500*t_s + 9));
  node_f.ay = x08 - (3000*t*x08*pow(t_s,5) + 5*t*(1500*x05 + 9*x08 - 750*t*x08)*pow(t_s,4) + 5*t*(1500*x02 + 18*x05 - 1500*x12 - 2250*t*x05 - 9*t*x08 + 250*pow(t,2)*x08)*pow(t_s,3) + 5*t*(18*x02 - 18*x12 - 2250*t*x02 - 18*t*x05 + 2250*t*x12 + 750*pow(t,2)*x05 + 3*pow(t,2)*x08)*pow(t_s,2) + 5*t*(18*t*x12 - 18*t*x02 + 750*pow(t,2)*x02 + 6*pow(t,2)*x05 - 750*pow(t,2)*x12)*t_s + 5*t*(6*pow(t,2)*x02 - 6*pow(t,2)*x12))/(pow(t_s,5)*(500*t_s + 9));
  node_f.az= x09 - (3000*t*x09*pow(t_s,5) + 5*t*(1500*x06 + 9*x09 - 750*t*x09)*pow(t_s,4) + 5*t*(1500*x03 + 18*x06 - 1500*x13 - 2250*t*x06 - 9*t*x09 + 250*pow(t,2)*x09)*pow(t_s,3) + 5*t*(18*x03 - 18*x13 - 2250*t*x03 - 18*t*x06 + 2250*t*x13 + 750*pow(t,2)*x06 + 3*pow(t,2)*x09)*pow(t_s,2) + 5*t*(18*t*x13 - 18*t*x03 + 750*pow(t,2)*x03 + 6*pow(t,2)*x06 - 750*pow(t,2)*x13)*t_s + 5*t*(6*pow(t,2)*x03 - 6*pow(t,2)*x13))/(pow(t_s,5)*(500*t_s + 9));

}

double QuadOptimalControl::time(double x01, double  x02, double  x03, double  x04, double  x05, double  x06, double  x07, double  x08, double  x09, double  x11, double  x12, double  x13, double  x14, double  x15, double  x16, double  x17, double  x18, double  x19){
  vector<double> coeffs{1, 0, (3*x07*x17)/2500 + (3*x08*x18)/2500 + (3*x09*x19)/2500 - (9*pow(x07,2))/5000 - (9*pow(x08,2))/5000 - (9*pow(x09,2))/5000 - (9*pow(x17,2))/5000 - (9*pow(x18,2))/5000 - (9*pow(x19,2))/5000, (12*x04*x17)/625 - (18*x05*x08)/625 - (18*x06*x09)/625 - (18*x04*x07)/625 - (12*x07*x14)/625 + (12*x05*x18)/625 - (12*x08*x15)/625 + (12*x06*x19)/625 - (12*x09*x16)/625 + (18*x14*x17)/625 + (18*x15*x18)/625 + (18*x16*x19)/625, (9*x01*x17)/125 - (9*x02*x08)/125 - (9*x03*x09)/125 - (9*x01*x07)/125 - (126*x04*x14)/625 + (9*x07*x11)/125 + (9*x02*x18)/125 - (126*x05*x15)/625 + (9*x08*x12)/125 + (9*x03*x19)/125 - (126*x06*x16)/625 + (9*x09*x13)/125 - (9*x11*x17)/125 - (9*x12*x18)/125 - (9*x13*x19)/125 - (72*pow(x04,2))/625 - (72*pow(x05,2))/625 - (72*pow(x06,2))/625 - (72*pow(x14,2))/625 - (72*pow(x15,2))/625 - (72*pow(x16,2))/625, (72*x04*x11)/125 - (72*x02*x05)/125 - (72*x03*x06)/125 - (72*x01*x14)/125 - (72*x01*x04)/125 - (72*x02*x15)/125 + (72*x05*x12)/125 - (72*x03*x16)/125 + (72*x06*x13)/125 + (72*x11*x14)/125 + (72*x12*x15)/125 + (72*x13*x16)/125, (36*x01*x11)/25 + (36*x02*x12)/25 + (36*x03*x13)/25 - (18*pow(x01,2))/25 - (18*pow(x02,2))/25 - (18*pow(x03,2))/25 - (18*pow(x11,2))/25 - (18*pow(x12,2))/25 - (18*pow(x13,2))/25};
  return minPositiveRoot(coeffs);
 }

double QuadOptimalControl::timePartialFinalState(double x01, double  x02, double  x03, double  x04, double  x05, double  x06, double  x07, double  x08, double  x09, double  x11, double  x12, double  x13){
    vector<double> coeffs{250000000, 9000000, - 300000*pow(x07,2) - 300000*pow(x08,2) - 300000*pow(x09,2) + 81000, - 9000*pow(x07,2) - 3000000*x04*x07 - 9000*pow(x08,2) - 3000000*x05*x08 - 9000*pow(x09,2) - 3000000*x06*x09, 4500000*x07*x11 - 4500000*x02*x08 - 81000*x04*x07 - 4500000*x03*x09 - 81000*x05*x08 - 81000*x06*x09 - 4500000*x01*x07 + 4500000*x08*x12 + 4500000*x09*x13 - 6750000*pow(x04,2) - 6750000*pow(x05,2) - 6750000*pow(x06,2) - 81*pow(x07,2) - 81*pow(x08,2) - 81*pow(x09,2), 18000000*x04*x11 - 18000000*x02*x05 - 126000*x01*x07 - 18000000*x03*x06 - 126000*x02*x08 - 648*x04*x07 - 126000*x03*x09 - 648*x05*x08 - 18000000*x01*x04 - 648*x06*x09 + 18000000*x05*x12 + 126000*x07*x11 + 18000000*x06*x13 + 126000*x08*x12 + 126000*x09*x13 - 153000*pow(x04,2) - 153000*pow(x05,2) - 153000*pow(x06,2), - 11250000*pow(x01,2) - 423000*x01*x04 + 22500000*x01*x11 - 972*x07*x01 - 11250000*pow(x02,2) - 423000*x02*x05 + 22500000*x02*x12 - 972*x08*x02 - 11250000*pow(x03,2) - 423000*x03*x06 + 22500000*x03*x13 - 972*x09*x03 - 972*pow(x04,2) + 423000*x04*x11 - 972*pow(x05,2) + 423000*x05*x12 - 972*pow(x06,2) + 423000*x06*x13 - 11250000*pow(x11,2) + 972*x07*x11 - 11250000*pow(x12,2) + 972*x08*x12 - 11250000*pow(x13,2) + 972*x09*x13, - 270000*pow(x01,2) + 540000*x01*x11 - 2592*x04*x01 - 270000*pow(x02,2) + 540000*x02*x12 - 2592*x05*x02 - 270000*pow(x03,2) + 540000*x03*x13 - 2592*x06*x03 - 270000*pow(x11,2) + 2592*x04*x11 - 270000*pow(x12,2) + 2592*x05*x12 - 270000*pow(x13,2) + 2592*x06*x13, - 1620*pow(x01,2) + 3240*x01*x11 - 1620*pow(x02,2) + 3240*x02*x12 - 1620*pow(x03,2) + 3240*x03*x13 - 1620*pow(x11,2) - 1620*pow(x12,2) - 1620*pow(x13,2)};
  return minPositiveRoot(coeffs);
}

