# include "double_integrator.h"
// highest order coeff first
// very very sketch, don't ask don't tell shhhh.....
double DoubleIntegrator::minPositiveRoot(vector<double> coeffs){
  auto start = std::chrono::system_clock::now();

  bool sign_left = peval(coeffs, 0.0) > 0;
  double step = 0.1;
  double max_range = 15.0;
  double lower_bound = 0.0;
  double upper_bound = 0.0;
  double retval;

  bool found_solution_interval = false;
  while ( upper_bound < max_range ){
    if ( (peval(coeffs,upper_bound) > 0) != sign_left){
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
  } catch (err_TooManyIteration&) {
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
double DoubleIntegrator::bisect( vector<double> coeffs, double low, double high){
  assert (high> low);
  bool sign_a = peval( coeffs, low) > 0;
  bool sign_b = peval( coeffs, high) > 0;

  assert (sign_a != sign_b);

  int loop_count = 0;
  while (high-low > 1e-3){
    double mid = (low+high)/2.0;
    if ( (peval( coeffs, mid) > 0) == sign_a) { low = mid; }
    else { high = mid; }
    loop_count++;
    if (loop_count > 100){ 
      cout << "bisect: max iter exceeded \n";
      throw new err_TooManyIteration; 
    }
  }
  return (high+low)/2.0;

}
double DoubleIntegrator::cost(double t_s, double x01, double x02, double x03, double x04, double x05, double x06, double x11, double x12, double x13, double x14, double x15, double x16)
{
  return t_s + (12*pow(x01,2))/pow(t_s,3) + (4*pow(x04,2))/t_s + (12*pow(x02,2))/pow(t_s,3) + (4*pow(x05,2))/t_s + (12*pow(x03,2))/pow(t_s,3) + (4*pow(x06,2))/t_s + (12*pow(x11,2))/pow(t_s,3) + (4*pow(x14,2))/t_s + (12*pow(x12,2))/pow(t_s,3) + (4*pow(x15,2))/t_s + (12*pow(x13,2))/pow(t_s,3) + (4*pow(x16,2))/t_s      + (12*x01*x04)/pow(t_s,2) + (12*x02*x05)/pow(t_s,2) + (12*x03*x06)/pow(t_s,2) - (24*x01*x11)/pow(t_s,3) + (12*x01*x14)/pow(t_s,2)          - (12*x04*x11)/pow(t_s,2) - (24*x02*x12)/pow(t_s,3) + (4*x04*x14)/t_s + (12*x02*x15)/pow(t_s,2) - (12*x05*x12)/pow(t_s,2)          - (24*x03*x13)/pow(t_s,3) + (4*x05*x15)/t_s + (12*x03*x16)/pow(t_s,2) - (12*x06*x13)/pow(t_s,2) + (4*x06*x16)/t_s - (12*x11*x14)/pow(t_s,2) - (12*x12*x15)/pow(t_s,2) - (12*x13*x16)/pow(t_s,2);
}


double DoubleIntegrator::costPartialFreeFinalState(double t_s, double x01, double x02, double x03, double x04, double x05, double x06, double x11, double x12, double x13)
{
  return t_s + (3*pow(x04,2) + 3*pow(x05,2) + 3*pow(x06,2))/t_s + (6*x01*x04 + 6*x02*x05 + 6*x03*x06 - 6*x04*x11 - 6*x05*x12 - 6*x06*x13)/pow(t_s,2) + (3*pow(x01,2) - 6*x01*x11 + 3*pow(x02,2) - 6*x02*x12 + 3*pow(x03,2) - 6*x03*x13 + 3*pow(x11,2) + 3*pow(x12,2) + 3*pow(x13,2))/pow(t_s,3);
}


double* DoubleIntegrator::interiorPositionPartialFinalState(double t_s, double x01, double x02, double x03, double x04, double x05, double x06, double x11, double x12, double x13)
{
  // prepare time vector
  double t_step = t_s/(interior_point_count+1);
  double t = t_step;
  for (int i=0;i<interior_point_count;i++){
    *(buffer_interior_pos + i*3 + 0) = x01 + t*x04 + (pow(t,2)*(t - 3*t_s)*(x01 - x11 + t_s*x04))/(2*pow(t_s,3));
    *(buffer_interior_pos + i*3 + 1) = x02 + t*x05 + (pow(t,2)*(t - 3*t_s)*(x02 - x12 + t_s*x05))/(2*pow(t_s,3));
    *(buffer_interior_pos + i*3 + 2) = x03 + t*x06 + (pow(t,2)*(t - 3*t_s)*(x03 - x13 + t_s*x06))/(2*pow(t_s,3));
    //*(buffer_interior_pos + i*r + 3) = x04 + (3*t*(t - 2*t_s)*(x01 - x11 + t_s*x04))/(2*pow(t_s,3));
    //*(buffer_interior_pos + i*r + 4) = x05 + (3*t*(t - 2*t_s)*(x02 - x12 + t_s*x05))/(2*pow(t_s,3));
    //*(buffer_interior_pos + i*r + 5) = x06 + (3*t*(t - 2*t_s)*(x03 - x13 + t_s*x06))/(2*pow(t_s,3));
    t += t_step;

  }
  return buffer_interior_pos;
}



double* DoubleIntegrator::interiorPosition( double t_s, double x01, double x02, double x03, double x04, double x05, double x06, double x11, double x12, double x13, double x14, double x15, double x16)
{
  // prepare time vector
  double t_step = t_s/(interior_point_count+1);
  double t = t_step;
  for (int i=0;i<interior_point_count;i++){
    *(buffer_interior_pos + i*3 + 0) = x11 + x14*(t - t_s) + (pow(t-t_s,3)*(2*x01 - 2*x11 + t_s*x04 + t_s*x14))/pow(t_s,3) + (pow(t-t_s,2)*(3*x01 - 3*x11 + t_s*x04 + 2*t_s*x14))/pow(t_s,2);
    *(buffer_interior_pos + i*3 + 1) = x12 + x15*(t - t_s) + (pow(t-t_s,3)*(2*x02 - 2*x12 + t_s*x05 + t_s*x15))/pow(t_s,3) + (pow(t-t_s,2)*(3*x02 - 3*x12 + t_s*x05 + 2*t_s*x15))/pow(t_s,2);
    *(buffer_interior_pos + i*3 + 2) = x13 + x16*(t - t_s) + (pow(t-t_s,3)*(2*x03 - 2*x13 + t_s*x06 + t_s*x16))/pow(t_s,3) + (pow(t-t_s,2)*(3*x03 - 3*x13 + t_s*x06 + 2*t_s*x16))/pow(t_s,2);
    //*(buffer_interior_pos + i*3 + 3) = x14 + (2*(t - t_s)*(3*x01 - 3*x11 + t_s*x04 + 2*t_s*x14))/pow(t_s,2) + (3*pow(t-t_s,2)*(2*x01 - 2*x11 + t_s*x04 + t_s*x14))/pow(t_s,3);
    //*(buffer_interior_pos + i*3 + 4) = x15 + (2*(t - t_s)*(3*x02 - 3*x12 + t_s*x05 + 2*t_s*x15))/pow(t_s,2) + (3*pow(t-t_s,2)*(2*x02 - 2*x12 + t_s*x05 + t_s*x15))/pow(t_s,3);
    //*(buffer_interior_pos + i*3 + 5) = x16 + (2*(t - t_s)*(3*x03 - 3*x13 + t_s*x06 + 2*t_s*x16))/pow(t_s,2) + (3*pow(t-t_s,2)*(2*x03 - 2*x13 + t_s*x06 + t_s*x16))/pow(t_s,3);
    t += t_step;
  }
  return buffer_interior_pos;
}
void DoubleIntegrator::state(double t, double t_s, double x01, double x02, double x03, double x04, double x05, double x06, double x11, double x12, double x13, double x14, double x15, double x16)
{
  *(buffer_interior_state  + 0) = x11 + x14*(t - t_s) + (pow(t-t_s,3)*(2*x01 - 2*x11 + t_s*x04 + t_s*x14))/pow(t_s,3) + (pow(t-t_s,2)*(3*x01 - 3*x11 + t_s*x04 + 2*t_s*x14))/pow(t_s,2);
  *(buffer_interior_state  + 1) = x12 + x15*(t - t_s) + (pow(t-t_s,3)*(2*x02 - 2*x12 + t_s*x05 + t_s*x15))/pow(t_s,3) + (pow(t-t_s,2)*(3*x02 - 3*x12 + t_s*x05 + 2*t_s*x15))/pow(t_s,2);
  *(buffer_interior_state  + 2) = x13 + x16*(t - t_s) + (pow(t-t_s,3)*(2*x03 - 2*x13 + t_s*x06 + t_s*x16))/pow(t_s,3) + (pow(t-t_s,2)*(3*x03 - 3*x13 + t_s*x06 + 2*t_s*x16))/pow(t_s,2);
  *(buffer_interior_state  + 3) = x14 + (2*(t - t_s)*(3*x01 - 3*x11 + t_s*x04 + 2*t_s*x14))/pow(t_s,2) + (3*pow(t-t_s,2)*(2*x01 - 2*x11 + t_s*x04 + t_s*x14))/pow(t_s,3);
  *(buffer_interior_state  + 4) = x15 + (2*(t - t_s)*(3*x02 - 3*x12 + t_s*x05 + 2*t_s*x15))/pow(t_s,2) + (3*pow(t-t_s,2)*(2*x02 - 2*x12 + t_s*x05 + t_s*x15))/pow(t_s,3);
  *(buffer_interior_state  + 5) = x16 + (2*(t - t_s)*(3*x03 - 3*x13 + t_s*x06 + 2*t_s*x16))/pow(t_s,2) + (3*pow(t-t_s,2)*(2*x03 - 2*x13 + t_s*x06 + t_s*x16))/pow(t_s,3);
}

void DoubleIntegrator::setFullStatePartialFinalState(double t_s, Node& node_i, Node& node_f)
{
  double x01 = node_i.x;
  double x02 = node_i.y;
  double x03 = node_i.z;
  double x04 = node_i.vx;
  double x05 = node_i.vy;
  double x06 = node_i.vz;
  double x11 = node_f.x;
  double x12 = node_f.y;
  double x13 = node_f.z;
  double t = t_s;
  node_f.vx = x04 + (3*t*(t - 2*t_s)*(x01 - x11 + t_s*x04))/(2*pow(t_s,3));
  node_f.vy = x05 + (3*t*(t - 2*t_s)*(x02 - x12 + t_s*x05))/(2*pow(t_s,3));
  node_f.vz = x06 + (3*t*(t - 2*t_s)*(x03 - x13 + t_s*x06))/(2*pow(t_s,3));

}

double DoubleIntegrator::timePartialFinalState(double x01, double  x02, double  x03, double  x04, double  x05, double  x06, double  x11, double  x12, double  x13)
{
    vector<double> coeffs{1, 0,
        - 3*pow(x04,2) - 3*pow(x05,2) - 3*pow(x06,2), 12*x04*x11 - 12*x02*x05 - 12*x03*x06 - 12*x01*x04 + 12*x05*x12 + 12*x06*x13,
        - 9*pow(x01,2) + 18*x01*x11 - 9*pow(x02,2) + 18*x02*x12 - 9*pow(x03,2) + 18*x03*x13 - 9*pow(x11,2) - 9*pow(x12,2) - 9*pow(x13,2)};
    return minPositiveRoot(coeffs);
}

double DoubleIntegrator::time(double x01, double  x02, double  x03, double  x04, double  x05, double  x06, double  x11, double  x12, double  x13, double  x14, double  x15, double  x16)
{
  vector<double> coeffs{1, 0,
      - 4*pow(x04,2) - 4*x04*x14 - 4*pow(x05,2) - 4*x05*x15 - 4*pow(x06,2) - 4*x06*x16 - 4*pow(x14,2) - 4*pow(x15,2) - 4*pow(x16,2), 
      24*x04*x11 - 24*x02*x05 - 24*x03*x06 - 24*x01*x14 - 24*x01*x04 - 24*x02*x15 + 24*x05*x12 - 24*x03*x16 + 24*x06*x13 + 24*x11*x14 + 24*x12*x15 + 24*x13*x16,
      - 36*pow(x01,2) + 72*x01*x11 - 36*pow(x02,2) + 72*x02*x12 - 36*pow(x03,2) + 72*x03*x13 - 36*pow(x11,2) - 36*pow(x12,2) - 36*pow(x13,2)};
    return minPositiveRoot(coeffs);
}
