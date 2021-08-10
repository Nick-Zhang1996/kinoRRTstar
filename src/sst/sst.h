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
#include <iostream>

#include "common.h"
#include "world.h"


namespace ob = ompl::base;
namespace oc = ompl::control;

class mySST{
  private:
    // solution time
    double duration;
    Waypoint waypoint;
    World world;
    Node start_node,goal_node;
    oc::SimpleSetup *ss;

    //static callback functions
    static bool isStateValid(World &world, const oc::SpaceInformation *si, const ob::State *state);
    static void propagate(const ob::State *start, const oc::Control *control, const double duration, ob::State *result);
    
  public:
    mySST(World& in_world, Node& in_start_node, Node& in_end_node, double in_duration);
    mySST();
    bool solve();
    Waypoint getWaypoint(int index);
    int getWaypointCount();
    void planWithSimpleSetup();

};

#endif
