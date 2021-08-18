#include "sst.h"
OmplBenchmark::OmplBenchmark():
  duration(10),
  waypoint(),
  world(0,0,0,0,0,0),
  start_node(),
  goal_node()
{

}

OmplBenchmark::OmplBenchmark(World& in_world, Node& in_start_node, Node& in_end_node):
  duration(0),
  waypoint(),
  world(in_world),
  start_node(in_start_node),
  goal_node(in_end_node)
{

  // state space
  // x,y,z, vx,vy,vz
  #ifdef SINGLE_INTEGRATOR_DYNAMICS
  int state_dim = 3;
  #endif
  #ifdef DOUBLE_INTEGRATOR_DYNAMICS
  int state_dim = 6;
  #endif
  #ifdef TRIPLE_INTEGRATOR_DYNAMICS
  int state_dim = 9;
  #endif

  int control_dim = 3;
  space = std::make_shared<ob::RealVectorStateSpace>(state_dim);
  bounds = new ob::RealVectorBounds(state_dim);
  // NOTE assuming world to be a cube centered at origin
  bounds->setLow(world.x_l);
  bounds->setHigh(world.x_h);
  space->setBounds(*bounds);
  space->setLongestValidSegmentFraction(0.02/space->getMaximumExtent());

  // control space
  auto cspace(std::make_shared<oc::RealVectorControlSpace>(space,control_dim));
  cbounds = new ob::RealVectorBounds(control_dim);
  cbounds->setLow(-10);
  cbounds->setHigh(10);
  cspace->setBounds(*cbounds);

  ss = new oc::SimpleSetup(cspace);
  ss->setStatePropagator(propagate);
  //ss->setOptimizationObjective(getMyPathLengthObjective(ss->getSpaceInformation()));
  optimization_objective = getMixedObjective(ss->getSpaceInformation());
  ss->setOptimizationObjective(optimization_objective);


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
  /*
  ob::ScopedState<ob::RealVectorStateSpace> goal(space);
  for (int i=0; i<state_dim; i++)
  {
    goal[i] = 0.0;
  }
  goal[0] = goal_node.x;
  goal[1] = goal_node.y;
  goal[2] = goal_node.z;
  */

  //ss->setStartAndGoalStates(start, goal);
  //auto goal = std::make_shared<myGoal>(ss->getSpaceInformation(), 8, -2, 3);
  auto goal = std::make_shared<myGoal>(ss->getSpaceInformation());
  goal->setPos(8, -2, 3);
  ss->addStartState(start);
  ss->setGoal(goal);
  // set planner
  ss->setPlanner(std::make_shared<oc::mySST>(ss->getSpaceInformation()));
  //ss->setPlanner(std::make_shared<oc::RRT>(ss->getSpaceInformation()));

  //oc::SimpleSetup *local_ss = ss;
  //World local_world = world;

  ss->setStateValidityChecker(
      [&local_world = world, &local_ss = ss](const ob::State *state) { return OmplBenchmark::isStateValid(local_world, local_ss->getSpaceInformation().get(), state); });
}

bool OmplBenchmark::solve(double in_duration)
{

  duration = in_duration;
  std::cout << "sst solving with time limit: " << duration << std::endl;
  ob::PlannerStatus solved = ss->solve(duration);
  if(solved)
  {
    std::cout << "Found solution:" << std::endl;
    // print the path to screen

    ss->getSolutionPath().printAsMatrix(std::cout);
    //ss->getSolutionPath().print(std::cout);
    double cost = ss->getSolutionPath().asGeometric().cost(optimization_objective).value();
    std::cout << "cost = " << cost << std::endl;
    return true;
  }
  else
  {
    std::cout << "No solution found" << std::endl;
    return false;
  }
}

bool OmplBenchmark::solveIncrementally(double in_duration, double step)
{

  duration = in_duration;
  std::cout << "sst solving with time limit: " << duration << std::endl;

  double total_time = 0;
  auto start = std::chrono::system_clock::now();
  //std::chrono::duration<double> elapsed_seconds = start-start;
  ob::PlannerStatus solved;
  int iter = (int) (in_duration/step);
  for (int i=0; i<iter; i++)
  {
    solved = ss->solve(step);

    // cost
    double min_cost_path = ss->getSolutionPath().asGeometric().cost(optimization_objective).value();
    //double min_cost_planner = ss->getPlanner()->as<oc::mySST>()->getBestCost().value();
    //double min_cost_planner = ss->getPlanner()->as<oc::RRTstar>()->bestCost().value();
    min_cost_hist.append(min_cost_path);

    // node count ?
    ob::PlannerData data(ss->getSpaceInformation());
    ss->getPlanner()->getPlannerData(data);
    int node_count = data.numVertices();
    node_count_hist.append(node_count);

    // solution count ?
    int solution_count = data.numGoalVertices();
    solution_count_hist.append(solution_count);

    //elapsed_seconds = std::chrono::system_clock::now()-start;
    //cout << "runtime: " << i*step << " nodes: " << node_count << " cost: " << min_cost << " solutions: " << solution_count << endl;
    cout << "runtime: " << i*step << " nodes: " << node_count << " cost(path): " << min_cost_path <<" cost(planenr): " << 0  << " solutions: " << solution_count << endl;

  }



  //ob::PlannerStatus solved = ss->solve(duration);
  if(solved)
  {
    std::cout << solved.asString() << std::endl;
    std::cout << "Found solution:" << std::endl;
    // print the path to screen

    ss->getSolutionPath().printAsMatrix(std::cout);
    //ss->getSolutionPath().print(std::cout);
    double cost = ss->getSolutionPath().asGeometric().cost(optimization_objective).value();
    std::cout << "cost = " << cost << std::endl;
    return true;
  }
  else
  {
    std::cout << "No solution found" << std::endl;
    return false;
  }
}

Waypoint OmplBenchmark::getWaypoint(int index)
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

int OmplBenchmark::getWaypointCount()
{
  return (int)ss->getSolutionPath().getStateCount();
}

// passing a class method to ompl is tricky(idk how to)
// so instead leave it out of class
bool OmplBenchmark::isStateValid(World &world, const oc::SpaceInformation *si, const ob::State *state)
{
  const auto* local_state = state->as<ob::RealVectorStateSpace::StateType>()->values;
  return world.checkNoCollision(local_state[0], local_state[1], local_state[2]);
}

// system dynamics
// single integrator
#ifdef SINGLE_INTEGRATOR_DYNAMICS
void OmplBenchmark::propagate(const ob::State *start, const oc::Control *control, const double duration, ob::State *result)
{
  const auto* state = start->as<ob::RealVectorStateSpace::StateType>()->values;
  const auto* ctrl = control->as<oc::RealVectorControlSpace::ControlType>()->values;
  auto* out = result->as<ob::RealVectorStateSpace::StateType>()->values; 
  // x+ = x + vx * dt
  out[0] = state[0] + ctrl[0] * duration;
  out[1] = state[1] + ctrl[1] * duration;
  out[2] = state[2] + ctrl[2] * duration;

}
#endif
// double integrator
#ifdef DOUBLE_INTEGRATOR_DYNAMICS
void OmplBenchmark::propagate(const ob::State *start, const oc::Control *control, const double duration, ob::State *result)
{
  const auto* state = start->as<ob::RealVectorStateSpace::StateType>()->values;
  const auto* ctrl = control->as<oc::RealVectorControlSpace::ControlType>()->values;
  auto* out = result->as<ob::RealVectorStateSpace::StateType>()->values; 
  // x+ = x + vx * dt + 1/2 * ax * dt*dt
  out[0] = state[0] + state[3] * duration + ctrl[0] * 0.5 * duration * duration;
  out[1] = state[1] + state[4] * duration + ctrl[1] * 0.5 * duration * duration;
  out[2] = state[2] + state[5] * duration + ctrl[2] * 0.5 * duration * duration;
  // vx+ = vx + ax * dt
  out[3] = state[3] + ctrl[0] * duration;
  out[4] = state[4] + ctrl[1] * duration;
  out[5] = state[5] + ctrl[2] * duration;

}
#endif
// triple integrator
#ifdef TRIPLE_INTEGRATOR_DYNAMICS
void OmplBenchmark::propagate(const ob::State *start, const oc::Control *control, const double duration, ob::State *result)
{
  const auto* state = start->as<ob::RealVectorStateSpace::StateType>()->values;
  const auto* ctrl = control->as<oc::RealVectorControlSpace::ControlType>()->values;
  auto* out = result->as<ob::RealVectorStateSpace::StateType>()->values; 
  // r+ = r + v * dt + 1/2 * a * dt2 + 1/6 * j * dt3
  out[0] = state[0] + state[3] * duration + 0.5 * state[6] * duration * duration + 1.0/6.0 * ctrl[0] * duration * duration * duration;
  out[1] = state[1] + state[3] * duration + 0.5 * state[7] * duration * duration + 1.0/6.0 * ctrl[1] * duration * duration * duration;
  out[2] = state[2] + state[3] * duration + 0.5 * state[8] * duration * duration + 1.0/6.0 * ctrl[2] * duration * duration * duration;
  // v+ = v + a*dt + 1/2*j*dt2 
  out[3] = state[3] + state[6] * duration + 0.5* ctrl[0] * duration * duration;
  out[4] = state[4] + state[7] * duration + 0.5* ctrl[1] * duration * duration;
  out[5] = state[5] + state[8] * duration + 0.5* ctrl[2] * duration * duration;
  // a+ = a + j * dt
  out[6] = state[6] + ctrl[0] * duration;
  out[7] = state[7] + ctrl[1] * duration;
  out[8] = state[8] + ctrl[2] * duration;

}
#endif

boost::python::list OmplBenchmark::getNodeCountHistPy()
{
  return node_count_hist;
}

boost::python::list OmplBenchmark::getMinCostHistPy()
{
  return min_cost_hist;
}

boost::python::list OmplBenchmark::getSolutionCountHistPy()
{
  return solution_count_hist;
}

void OmplBenchmark::planWithSimpleSetup()
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


ob::OptimizationObjectivePtr OmplBenchmark::getMyPathLengthObjective(const ob::SpaceInformationPtr& si)
{
    return ob::OptimizationObjectivePtr(new myPathLengthOptimizationObjective(si));
}

ob::OptimizationObjectivePtr OmplBenchmark::getMixedObjective(const ob::SpaceInformationPtr& si)
{
    return ob::OptimizationObjectivePtr(new mixedOptimizationObjective(si));
}
// -------------  myPathLengthOptimizationObjective --------

myPathLengthOptimizationObjective::myPathLengthOptimizationObjective(const ob::SpaceInformationPtr &si)
  : ob::StateCostIntegralObjective(si,true)
{
    description_ = "My Path Length";
}

ob::Cost myPathLengthOptimizationObjective::stateCost(const ob::State *) const
{
    return identityCost();
}

ob::Cost myPathLengthOptimizationObjective::motionCost(const ob::State *s1, const ob::State *s2) const
{
 // TODO
  double* _s1 = s1->as<ob::RealVectorStateSpace::StateType>()->values;
  double* _s2 = s2->as<ob::RealVectorStateSpace::StateType>()->values;
  auto sqr = [](double a){return a*a;};
  double distance = sqr(_s1[0]-_s2[0]) + sqr(_s1[1]-_s2[1]) + sqr(_s1[2]-_s2[2]);
  distance = sqrt(distance);
  return ob::Cost(distance);
}

ob::Cost myPathLengthOptimizationObjective::motionCostHeuristic(const ob::State *s1,
                                                                                  const ob::State *s2) const
{
    return motionCost(s1, s2);
}

// -------------  mixedOptimizationObjective --------
//  integrate: (1 + u'Ru) dt, R = lambda * I (uniform cost)

mixedOptimizationObjective::mixedOptimizationObjective(const ob::SpaceInformationPtr &si)
  : ob::StateCostIntegralObjective(si,true)
{
    description_ = "Mixed Time + Control Effort";
}

// is this really used?
ob::Cost mixedOptimizationObjective::stateCost(const ob::State *) const
{
  std::cout << "stateCost called" << std::endl;
  return identityCost();
}
// 3 dof state
#ifdef SINGLE_INTEGRATOR_DYNAMICS
ob::Cost mixedOptimizationObjective::motionCost(const ob::State *s1, const ob::State *s2) const
{
  double* _s1 = s1->as<ob::RealVectorStateSpace::StateType>()->values;
  double* _s2 = s2->as<ob::RealVectorStateSpace::StateType>()->values;
  auto sqr = [](double a){return a*a;};
  double distance = sqr(_s1[0]-_s2[0]) + sqr(_s1[1]-_s2[1]) + sqr(_s1[2]-_s2[2]);
  distance = sqrt(distance);
  return ob::Cost(distance);
}
#endif

// 6 dof state
#ifdef DOUBLE_INTEGRATOR_DYNAMICS
ob::Cost mixedOptimizationObjective::motionCost(const ob::State *s1, const ob::State *s2) const
{
  double* _s1 = s1->as<ob::RealVectorStateSpace::StateType>()->values;
  double* _s2 = s2->as<ob::RealVectorStateSpace::StateType>()->values;
  // calculate time interval
  double dt = (_s2[0] - _s1[0]) / (_s2[3] + _s1[3]) * 2; 
  double dt_alt = (_s2[1] - _s1[1]) / (_s2[4] + _s1[4]) * 2; 
  if (dt/dt_alt > 1.01 || dt/dt_alt < 0.99 || dt < 0)
  {
    std::cout << " --- bu hao le ---" << std::endl;
  }
  // calculate  control effort
  double ax = (_s2[3] - _s2[3]) / dt;
  double ay = (_s2[4] - _s2[4]) / dt;
  double az = (_s2[5] - _s2[5]) / dt;
  double cost = dt * (ax*ax + ay*ay + az*az) + dt;

  return ob::Cost(cost);
}
#endif

// 9 dof state
#ifdef TRIPLE_INTEGRATOR_DYNAMICS
ob::Cost mixedOptimizationObjective::motionCost(const ob::State *s1, const ob::State *s2) const
{
  double* _s1 = s1->as<ob::RealVectorStateSpace::StateType>()->values;
  double* _s2 = s2->as<ob::RealVectorStateSpace::StateType>()->values;
  // calculate time interval
  double dt = (_s2[3] - _s1[3]) / (_s2[6] + _s1[6]) * 2; 
  double dt_alt = (_s2[4] - _s1[4]) / (_s2[7] + _s1[7]) * 2; 
  if (dt/dt_alt > 1.01 || dt/dt_alt < 0.99 || dt < 0)
  {
    std::cout << " --- bu hao le ---" << std::endl;
  }
  // calculate  control effort
  double ax = (_s2[3] - _s2[3]) / dt;
  double ay = (_s2[4] - _s2[4]) / dt;
  double az = (_s2[5] - _s2[5]) / dt;
  double cost = dt * (ax*ax + ay*ay + az*az)*0.0002 + dt;

  return ob::Cost(cost);
}
#endif

ob::Cost mixedOptimizationObjective::motionCostHeuristic(const ob::State *s1,
                                                                                  const ob::State *s2) const
{
    return motionCost(s1, s2);
}

bool myGoal::isSatisfied(const ob::State* st) const
{
  double* state = st->as<ob::RealVectorStateSpace::StateType>()->values;
  if (std::abs(_x-state[0]) < 1.5f && std::abs(_y-state[1]) < 1.5f && std::abs(_z-state[2]) < 1.5f)
  {
    return true;
  } else 
  {
    return false;
  }
}


// --------- main --------------

int main(int, char**)
{
  OmplBenchmark sst;
  sst.planWithSimpleSetup();
  return 0;

}
