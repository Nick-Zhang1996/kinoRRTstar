#include <unsupported/Eigen/Polynomials>
#include <iostream>
#include <vector>

int main(){
  // high order coeffs first
  std::vector<double> coeffs{1.0000, 0, -0.0002, -0.0670, -6.9903, 1.8660, -8.2350};
  std::reverse(coeffs.begin(), coeffs.end());

  Eigen::VectorXd poly(7);

  for (int i=0; i<coeffs.size(); i++){
    poly(i) = coeffs[i];
  }
  // lowest order coeff first
  Eigen::PolynomialSolver<double, Eigen::Dynamic> solver;

  solver.compute(poly);
  const Eigen::PolynomialSolver<double, Eigen::Dynamic>::RootsType &r = solver.roots();
  std::cout << "rows=" << r.rows() << " , cols="  << r.cols() << "\n";
  for (int i=0; i<r.rows(); i++){
    std::cout << r(i,0).real() << " +" << r(i,0).imag() << "i\n";
  }

}
