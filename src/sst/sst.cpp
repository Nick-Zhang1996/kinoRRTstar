#include "sst.h"
#include "world.h"
  
namespace ob = ompl::base;
namespace oc = ompl::control;

bool isStateValid(World &world, const oc::SpaceInformation *si, const ob::State *state)
{
  const auto* local_state = state->as<ob::RealVectorStateSpace::StateType>()->values;
  return world.checkNoCollision(local_state[0], local_state[1], local_state[2]);
}

// system dynamics
void propagate(const ob::State *start, const oc::Control *control, const double duration, ob::State *result)
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

void planWithSimpleSetup()
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

  oc::SimpleSetup ss(cspace);
  ss.setStatePropagator(propagate);

  ss.setStateValidityChecker(
      [&world, &ss](const ob::State *state) { return isStateValid(world, ss.getSpaceInformation().get(), state); });

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

  ss.setStartAndGoalStates(start, goal);
  // set planner
  ss.setPlanner(std::make_shared<oc::SST>(ss.getSpaceInformation()));
  ob::PlannerStatus solved = ss.solve(30.0);
  if(solved)
  {
      std::cout << "Found solution:" << std::endl;
      // print the path to screen

      ss.getSolutionPath().printAsMatrix(std::cout);
      ss.getSolutionPath().print(std::cout);
  }
  else
      std::cout << "No solution found" << std::endl;


}

int main(int, char**)
{
   planWithSimpleSetup();
   return 0;

}
