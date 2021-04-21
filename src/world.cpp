#include "world.h"

bool World::checkNoPathCollision(double * pos_buffer){
  for (int i=0; i<interior_point_n; i++){
    Point p;
    // x
    p.x = *(pos_buffer + i*3 + 0);
    p.y = *(pos_buffer + i*3 + 1);
    p.z = *(pos_buffer + i*3 + 2);
    if (!checkNoCollision( p )){ return false; }
  }
  return true;
}
