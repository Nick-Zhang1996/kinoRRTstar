#include <iostream>
#include "quad_optimal_control.h"

int main(){
  oc = QuadOptimalControl(10);
  auto val = oc.peval(4,vector<float>{1, -2, 1});
  std::cout << val << "\n"; // expect 9

  val = oc.peval(4,vector<float>{1, 1, 0, 0});
  std::cout << val << "\n"; // expect 12

  auto sol = oc.peval(4,vector<float>{1, -2, 1});
  for (auto i=sol.begin(); i<sol.end(); i++){
    std::cout << i << ",";
  }
  std::cout << "\n"; // expect 2, -2, 1

  sol = oc.peval(4,vector<float>{2, 1, 0, 0});
  for (auto i=sol.begin(); i<sol.end(); i++){
    std::cout << i << ",";
  }
  std::cout << "\n"; // expect 6 2 0

  sol = oc.roots(vector<float>{1, -2, 0});
  for (auto i=sol.begin(); i<sol.end(); i++){
    std::cout << i << ",";
  }
  std::cout << "\n"; // expect 0, 2

  sol = oc.minPositiveRoot(vector<float>{1, -2, 0});
  for (auto i=sol.begin(); i<sol.end(); i++){
    std::cout << i << ",";
  }
  std::cout << "\n"; // expect 2
}
