#include "main.h"
#include "common.h"
#include "tree.h"
#include "world.h"
#include "kino_rrt_star.h"
#include "quad_optimal_control.h"

int main(){
  World world(20,10,10);
  Box obs1(6,0,0,6+2,5,10);
  Box obs2(6,5,0,6+2,5+5,5);
  Box obs3(12,5,0,12+2,5+5,10);
  Box obs4(12,0,5,12+2,5,5+5);


  Node start_node(2,2,2);
  Node goal_node( 18, 8, 8 );

  // TODO check obstacle size
  world.addObstacle(obs1);
  world.addObstacle(obs2);
  world.addObstacle(obs3);
  world.addObstacle(obs4);

  KinoRrtStar rrt(world, start_node, goal_node, 600, 10);
  rrt.run();
  cout << "program complete\n";


}
