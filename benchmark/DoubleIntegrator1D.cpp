#include <ompl/control/SpaceInformation.h>
 #include <ompl/base/goals/GoalState.h>
 #include <ompl/base/spaces/SE2StateSpace.h>
 #include <ompl/control/spaces/RealVectorControlSpace.h>
 #include <ompl/control/planners/kpiece/KPIECE1.h>
 #include <ompl/control/planners/rrt/RRT.h>
 #include <ompl/control/planners/est/EST.h>
 #include <ompl/control/planners/sst/SST.h>
 #include <ompl/control/planners/syclop/SyclopRRT.h>
 #include <ompl/control/planners/syclop/SyclopEST.h>
 #include <ompl/control/planners/pdst/PDST.h>
 #include <ompl/control/planners/syclop/GridDecomposition.h>
 #include <ompl/control/SimpleSetup.h>
 #include <ompl/config.h>
 #include <iostream>
  
namespace ob = ompl::base;
namespace oc = ompl::control;

//TODO
bool isStateValid(const oc::SpaceInformation *si, const ob::State *state)
{
  return true;
}

void propagate(const ob::State *start, const oc::Control *control, const double duration, ob::State *result)
{
  const auto* state = start->as<ob::RealVectorStateSpace::StateType>()->values;
  const auto* ctrl = control->as<oc::RealVectorControlSpace::ControlType>()->values;
  auto* out = result->as<ob::RealVectorStateSpace::StateType>()->values; // x+ = x + vx * dt + 1/2 * ax * dt*dt

  out[0] = state[0] + state[1] * duration + ctrl[0] * 0.5 * duration * duration;
  // vx+ = vx + ax * dt
  out[1] = state[1] + ctrl[0] * duration;

}

void planWithSimpleSetup()
{
  // state space
  // x,y,z, vx,vy,vz
  int state_dim = 2;
  int control_dim = 1;
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
      [&ss](const ob::State *state) { return isStateValid(ss.getSpaceInformation().get(), state); });

  // start state
  ob::ScopedState<ob::RealVectorStateSpace> start(space);
  for (int i=0; i<state_dim; i++)
  {
    start[i] = 0.0;
  }

  ob::ScopedState<ob::RealVectorStateSpace> goal(space);
  goal[0] = 10.0;
  goal[1] = 0.0;

  ss.setStartAndGoalStates(start, goal, 0.01);
  // set planner
  ss.setPlanner(std::make_shared<oc::SST>(ss.getSpaceInformation()));
  ob::PlannerStatus solved = ss.solve(10.0);
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
