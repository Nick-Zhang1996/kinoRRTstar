#include "sst.h"
#include "world.h"
SST::SST():
  duration(10),
  waypoint(),
  world(0,0,0,0,0,0),
  start_node(),
  goal_node()
{

}

SST::SST(World& in_world, Node& in_start_node, Node& in_end_node, double in_duration):
  duration(in_duration),
  waypoint(),
  world(in_world),
  start_node(in_start_node),
  goal_node(in_end_node)
{

  // state space
  // x,y,z, vx,vy,vz
  int state_dim = 6;
  int control_dim = 3;
  auto space(std::make_shared<ob::RealVectorStateSpace>(state_dim));
  ob::RealVectorBounds bounds(state_dim);
  // NOTE assuming world to be a cube centered at origin
  bounds.setLow(world.x_l);
  bounds.setHigh(world.x_h);
  space->setBounds(bounds);

  // control space
  auto cspace(std::make_shared<oc::RealVectorControlSpace>(space,control_dim));
  ob::RealVectorBounds cbounds(control_dim);
  cbounds.setLow(-10);
  cbounds.setHigh(10);
  cspace->setBounds(cbounds);

  ss = new oc::SimpleSetup(cspace);
  ss->setStatePropagator(propagate);
  oc::SimpleSetup *local_ss = ss;

  ss->setStateValidityChecker(
      [&in_world, &local_ss](const ob::State *state) { return SST::isStateValid(in_world, local_ss->getSpaceInformation().get(), state); });

  // start state
  ob::ScopedState<ob::RealVectorStateSpace> start(space);
  for (int i=0; i<state_dim; i++)
  {
    start[i] = 0.0;
  }
  start[0] = start_node.x;
  start[1] = start_node.y;
  start[2] = start_node.z;


  // goal state
  ob::ScopedState<ob::RealVectorStateSpace> goal(space);
  for (int i=0; i<state_dim; i++)
  {
    goal[i] = 0.0;
  }
  goal[0] = goal_node.x;
  goal[1] = goal_node.y;
  goal[2] = goal_node.z;

  ss->setStartAndGoalStates(start, goal);
  // set planner
  ss->setPlanner(std::make_shared<oc::SST>(ss->getSpaceInformation()));

}

Waypoint SST::getWaypoint(int index)
{
  ob::State* _state = ss->getSolutionPath().getState(index);
  double* states = _state->as<ob::RealVectorStateSpace::StateType>()->values;

  waypoint.t = 0.0;
  waypoint.valid = true;
  waypoint.x = states[0];
  waypoint.y = states[1];
  waypoint.z = states[2];
  waypoint.vx = states[3];
  waypoint.vy = states[4];
  waypoint.vz = states[5];
  waypoint.ax = states[6];
  waypoint.ay = states[7];
  waypoint.az = states[8];
  return waypoint;
}

int SST::getWaypointCount()
{
  return (int)ss->getSolutionPath().getStateCount();
}

// passing a class method to ompl is tricky(idk how to)
// so instead leave it out of class
bool SST::isStateValid(World &world, const oc::SpaceInformation *si, const ob::State *state)
{
  const auto* local_state = state->as<ob::RealVectorStateSpace::StateType>()->values;
  return world.checkNoCollision(local_state[0], local_state[1], local_state[2]);
}

// system dynamics
void SST::propagate(const ob::State *start, const oc::Control *control, const double duration, ob::State *result)
{
  const auto* state = start->as<ob::RealVectorStateSpace::StateType>()->values;
  const auto* ctrl = control->as<oc::RealVectorControlSpace::ControlType>()->values;
  auto* out = result->as<ob::RealVectorStateSpace::StateType>()->values; // x+ = x + vx * dt + 1/2 * ax * dt*dt

  out[0] = state[0] + state[3] * duration + ctrl[0] * 0.5 * duration * duration;
  out[1] = state[1] + state[4] * duration + ctrl[1] * 0.5 * duration * duration;
  out[2] = state[2] + state[5] * duration + ctrl[2] * 0.5 * duration * duration;
  // vx+ = vx + ax * dt
  out[3] = state[3] + ctrl[0] * duration;
  out[4] = state[4] + ctrl[1] * duration;
  out[5] = state[5] + ctrl[2] * duration;

}

void SST::planWithSimpleSetup()
{
  // create world 
  World world(-10,10,-10,10,-10,10);
  // x_l, x_h, y_l,y_h, z_l, z_h
  Box obs1(-4,-2, -10,0, -10,10);
  Box obs2(-4,-2, 0,10, -10,0);
  Box obs3(2,4, 0,10, -10,10);
  Box obs4(2,4, -10,0, 0,10);

  world.addObstacle(obs1);
  world.addObstacle(obs2);
  world.addObstacle(obs3);
  world.addObstacle(obs4);
  
  // state space
  // x,y,z, vx,vy,vz
  int state_dim = 6;
  int control_dim = 3;
  auto space(std::make_shared<ob::RealVectorStateSpace>(state_dim));
  ob::RealVectorBounds bounds(state_dim);
  bounds.setLow(-15);
  bounds.setHigh(15);
  space->setBounds(bounds);

  // control space
  auto cspace(std::make_shared<oc::RealVectorControlSpace>(space,control_dim));
  ob::RealVectorBounds cbounds(control_dim);
  cbounds.setLow(-10);
  cbounds.setHigh(10);
  cspace->setBounds(cbounds);

  ss = new oc::SimpleSetup(cspace);
  ss->setStatePropagator(propagate);
  oc::SimpleSetup *local_ss = ss;

  ss->setStateValidityChecker(
      [&world, &local_ss](const ob::State *state) { return isStateValid(world, local_ss->getSpaceInformation().get(), state); });

  // start state
  ob::ScopedState<ob::RealVectorStateSpace> start(space);
  for (int i=0; i<state_dim; i++)
  {
    start[i] = 0.0;
  }
  start[0] = -8;
  start[1] = 2;
  start[2] = -3;


  // goal state
  ob::ScopedState<ob::RealVectorStateSpace> goal(space);
  for (int i=0; i<state_dim; i++)
  {
    goal[i] = 0.0;
  }
  goal[0] = 8;
  goal[1] = -2;
  goal[2] = 3;

  ss->setStartAndGoalStates(start, goal);
  // set planner
  ss->setPlanner(std::make_shared<oc::SST>(ss->getSpaceInformation()));
  ob::PlannerStatus solved = ss->solve(30.0);
  if(solved)
  {
      std::cout << "Found solution:" << std::endl;
      // print the path to screen

      ss->getSolutionPath().printAsMatrix(std::cout);
      ss->getSolutionPath().print(std::cout);
  }
  else
      std::cout << "No solution found" << std::endl;


}

bool SST::solve()
{
  std::cout << "sst solving with time limit: " << duration << std::endl;
  ob::PlannerStatus solved = ss->solve(duration);
  if(solved)
  {
    std::cout << "Found solution:" << std::endl;
    // print the path to screen

    ss->getSolutionPath().printAsMatrix(std::cout);
    ss->getSolutionPath().print(std::cout);
    return true;
  }
  else
  {
    std::cout << "No solution found" << std::endl;
    return false;
  }

}

// --------- main --------------

int main(int, char**)
{
  SST sst;
  sst.planWithSimpleSetup();
  return 0;

}
