// kinodynamic RRT* implementation
#ifndef KINO_RRT_STAR_H
#define KINO_RRT_STAR_H

#include "common.h"
#include "tree.h"
#include <iostream>
#include "world.h"
#include <cstdlib>
#include <ctime>
#ifdef TRIPLE_INTEGRATOR_DYNAMICS
#include "quad_optimal_control.h"
#elif defined DOUBLE_INTEGRATOR_DYNAMICS
#include "double_integrator.h"
#endif

#include <limits>
#include <math.h>
#include <algorithm>
#include <chrono>
#include <boost/python.hpp>
using std::cout;
using std::rand;
using std::min;

class KinoRrtStar{
  private:
    // tree
    Tree tree;
    Node start_node, end_node;
    int target_node_count;
    World world;
    #ifdef TRIPLE_INTEGRATOR_DYNAMICS
    QuadOptimalControl oc;
    #elif defined DOUBLE_INTEGRATOR_DYNAMICS
    DoubleIntegrator oc;
    #endif
    double overall_lowest_cost;
    int overall_lowest_cost_id;
    int interior_point_count;
    int rewire_count;
    int neighbour_total_count;
    int neighbour_count;
    list<Waypoint>::iterator waypoints_iter;
    list<Waypoint> waypoints;

    // progress statistics 
    // log various progress metrics by time interval log_interval_s
    // 0: 0s, 1:log_interval_s, 2: 2*log_interval_s
    double log_interval_s;
    vector<int> node_count_hist;
    vector<double> min_cost_hist;
    vector<int> solution_count_hist;


  public:
    KinoRrtStar(World& in_world, Node& in_start_node, Node& in_end_node, int in_interior_point_count );
    // n_nodes: number of nodes to add to tree, if 0 then stop after first solution
    void buildTree(int n_nodes);
    bool sampleSpace();
    void run(int in_target_node_count);
    void runWithTimeLimit(double duration);
    void buildTreeTillNodeCount();
    void buildTreeTillFirstSolution();
    void sampleNode();
    void showResult();

    // prepare waypoints for optimal solution, including interior points
    // return number of waypoints
    int prepareSolutionWithInteriorPoints();
    int prepareSolution();
    double getTrajectoryTime();
    Waypoint getTrajectory(double t);
    Waypoint waypoint;
    // after prepareSolution() is called
    // repeatedly call this function to get all waypoints
    // this function should be called <number of waypoints> times
    // if called afterwards, Waypoint.valid will be false;
    Waypoint getNextWaypoint();

    // generate a random point in world, only x,y,z are concerned
    // TODO test
    Node getRandomNode();
    // move node to target if its further than radius
    // TODO test
    void moveToVicinity(Node& node, Node& target, double radius);
    bool connectToGoal(Node& node);
    double getNeighnourRadius();
    double getSampleRadius();

    double sqr(double a) { return a*a; }
    double dist(Node& a, Node& b) { return sqrt(sqr(a.x-b.x) + sqr(a.y-b.y) + sqr(a.z-b.z)); }

    vector<int> getNodeCountHist();
    vector<double> getMinCostHist();
    vector<int> getSolutionCountHist();

    boost::python::list getNodeCountHistPy();
    boost::python::list getMinCostHistPy();
    boost::python::list getSolutionCountHistPy();

};

#endif
