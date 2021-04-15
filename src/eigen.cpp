#include <unsupported/Eigen/Polynomials>
#include <iostream>

int main(){
  // lowest order coeff first
  Eigen::Vector3d coeff(-0.2151830138973625, -0.3111717537041549, 0.708563939215852);
  Eigen::PolynomialSolver<double, Eigen::Dynamic> solver;

  solver.compute(coeff);
  const Eigen::PolynomialSolver<double, Eigen::Dynamic>::RootsType &r = solver.roots();
  std::cout << "rows=" << r.rows() << " , cols="  << r.cols();
  std::cout << r(0,0).real();
  std::cout << r(1,0).imag();

}
