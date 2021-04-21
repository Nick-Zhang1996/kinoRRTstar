#include "main.h"
#include "common.h"
#include "tree.h"
#include "world.h"
#include "kino_rrt_star.h"
#include "quad_optimal_control.h"

int main(){
  World world(10,10,10);
  Box obs1(3,6,0,4,0,10);
  Box obs2(3,6,6,10,0,10);
  Node start_node(0,0,0);
  Node goal_node( 8, 5, 5 );

  world.addObstacle(obs1);
  world.addObstacle(obs2);

  KinoRrtStar rrt(world, start_node, goal_node, 6000, 10);
  rrt.run();


}
