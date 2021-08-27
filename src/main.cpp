#include "main.h"
#include "common.h"
#include "tree.h"
#include "world.h"
#include "kino_rrt_star.h"

int main(){
  /*
  World world(0,20,0,10,0,10);
  // x_l, x_h, y_l,y_h, z_l, z_h
  Box obs1(6,8, 0,5, 0,10);
  Box obs2(6,8, 5,10, 0,5);
  Box obs3(12,14, 5,10, 0,10);
  Box obs4(12,14, 0,5, 5,10);
  Node start_node(2,2,2);
  Node goal_node( 18, 8, 8 );
  */

  World world(-10,10,-10,10,-10,10);

  Box obs1(-4,-2, -10,0, -10,10);
  Box obs2(-4,-2, 0,10, -10,0);
  Box obs3(2,4, 0,10, -10,10);
  Box obs4(2,4, -10,0, 0,10);

  Node start_node(-8, 2, -3);
  Node goal_node(8, -2, 3);

  // TODO check obstacle size
  world.addObstacle(obs1);
  world.addObstacle(obs2);
  world.addObstacle(obs3);
  world.addObstacle(obs4);


  KinoRrtStar rrt(world, start_node, goal_node, 10);

  auto start = std::chrono::system_clock::now();

  rrt.run(2000);
  int retval = rrt.prepareSolution();
  if (!retval){
    cout << "RRT cannot find solution" << endl;
    return 1;
  }
  while (true){
    auto p = rrt.getNextWaypoint();
    if (!p.valid){break;}
  }

      
  auto end = std::chrono::system_clock::now();
  std::chrono::duration<double> elapsed_seconds = end-start;

  cout << "rrt complete in " << elapsed_seconds.count() << "s" << endl;

  // find trajectory to ddebug getTrajectory;
  double traj_t = rrt.getTrajectoryTime();
  for (double t=0.0; t<traj_t; t+=0.01){
    rrt.getTrajectory(t);
  }


}

