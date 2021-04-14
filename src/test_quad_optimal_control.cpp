#include <iostream>
#include "quad_optimal_control.h"

int main(){
  QuadOptimalControl oc(10);
  // peval
  auto val = oc.peval(vector<float>{1, -2, 1}, 4);
  std::cout << val << "\n"; // expect 9

  val = oc.peval(vector<float>{1, 1, 0, 0}, 2);
  std::cout << val << "\n"; // expect 12

  // derivative
  auto sol = oc.derivative(vector<float>{1, -2, 1});
  for (auto i=sol.begin(); i<sol.end(); i++){
    std::cout << *i << ",";
  }
  std::cout << "\n"; // expect 2, -2

  sol = oc.derivative(vector<float>{2, 1, 0, 0});
  for (auto i=sol.begin(); i<sol.end(); i++){
    std::cout << *i << ",";
  }
  std::cout << "\n"; // expect 6 2 0
  // falsePosition
  vector<float> coeffs{1,-2,0};
  auto temp = oc.falsePosition(coeffs, -10, 1);

  // roots
  vector<float> vec{1, -2, 0};
  auto root = oc.roots(vec);
  for (auto i=root.begin(); i<root.end(); i++){
    std::cout << *i << ",";
  }
  std::cout << "\n"; // expect 0, 2

  val = oc.minPositiveRoot(vector<float>{1, -2, 0});
  std::cout << val << "\n"; // expect 2
}