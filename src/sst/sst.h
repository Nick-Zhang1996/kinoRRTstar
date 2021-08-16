#ifndef SST_H
#define SST_H
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
#include <ompl/base/objectives/StateCostIntegralObjective.h>

#include <iostream>
#include <cmath>
#include <chrono>
#include <boost/python.hpp>
#include "common.h"
#include "world.h"
#define SINGLE_INTEGRATOR_DYNAMICS


namespace ob = ompl::base;
namespace oc = ompl::control;

class myPathLengthOptimizationObjective : public ob::StateCostIntegralObjective
{
  public:
    myPathLengthOptimizationObjective(const ob::SpaceInformationPtr &si);
    ob::Cost stateCost(const ob::State *s) const override;
    ob::Cost motionCost(const ob::State *s1, const ob::State *s2) const override;
    ob::Cost motionCostHeuristic(const ob::State *s1, const ob::State *s2) const override;
    //ob::InformedSamplerPtr allocInformedStateSampler(const ob::ProblemDefinitionPtr &probDefn,
                                                  //unsigned int maxNumberCalls) const override;
};

class mixedOptimizationObjective : public ob::StateCostIntegralObjective
{
  public:
    mixedOptimizationObjective(const ob::SpaceInformationPtr &si);
    ob::Cost stateCost(const ob::State *s) const override;
    ob::Cost motionCost(const ob::State *s1, const ob::State *s2) const override;
    ob::Cost motionCostHeuristic(const ob::State *s1, const ob::State *s2) const override;
    //ob::InformedSamplerPtr allocInformedStateSampler(const ob::ProblemDefinitionPtr &probDefn,
                                                  //unsigned int maxNumberCalls) const override;
};

namespace ompl
{
  namespace control
  { 
    class mySST : public SST
    {
      using SST::SST;
      public:
        base::Cost getBestCost(){ return prevSolutionCost_;}
    };
  }
}

class myGoal : public ompl::base::Goal
{
  private:
    double _x,_y,_z;
  public:
    myGoal(const ob::SpaceInformationPtr &si, double x, double y, double z) : 
      ob::Goal(si),
      _x(x),
      _y(y),
      _z(z)
    {
    }

    virtual bool isSatisfied(const ob::State* st) const;
};



class OmplBenchmark
{
  private:
    // solution time
    double duration;
    Waypoint waypoint;
    World world;
    Node start_node,goal_node;
    std::shared_ptr<ob::RealVectorStateSpace> space;
    ob::RealVectorBounds* bounds;
    ob::RealVectorBounds* cbounds;
    oc::SimpleSetup *ss;
    ob::OptimizationObjectivePtr optimization_objective;
    boost::python::list node_count_hist;
    boost::python::list min_cost_hist;
    boost::python::list solution_count_hist;

    

    //static callback functions
    static bool isStateValid(World &world, const oc::SpaceInformation *si, const ob::State *state);
    static void propagate(const ob::State *start, const oc::Control *control, const double duration, ob::State *result);
    static ob::OptimizationObjectivePtr getMyPathLengthObjective(const ob::SpaceInformationPtr& si);
    static ob::OptimizationObjectivePtr getMixedObjective(const ob::SpaceInformationPtr& si);
    
  public:
    OmplBenchmark(World& in_world, Node& in_start_node, Node& in_end_node);
    OmplBenchmark();
    bool solve(double duration);
    bool solveIncrementally(double in_duration, double step);
    Waypoint getWaypoint(int index);
    int getWaypointCount();
    void planWithSimpleSetup();
    boost::python::list getNodeCountHistPy();
    boost::python::list getMinCostHistPy();
    boost::python::list getSolutionCountHistPy();
};

#endif
