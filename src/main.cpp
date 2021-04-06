// Main for kino RRT
//

#include <iostream>
#include "world.h"

int main(){
  // defind world
  World world(10,10,10);
  Point p1{1,1,1};
  Point p2{5,5,5};

  Box obs1(3,6,3,6,3,6);
  world.addObstacle(&obs1);
  cout << world.checkNoCollision(p1) << "\n" ;
  cout << world.checkNoCollision(p2) << "\n";
  
  return 0;
}
