// kinodynamic RRT* implementation
#ifndef KINO_RRT_STAR_H
#define KINO_RRT_STAR_H

#include "common.h"
#include "tree.h"
#include <iostream>
#include "world.h"
#include <cstdlib>
#include <ctime>
#include "quad_optimal_control.h"
using std::cout;
using std::rand;

class KinoRrtStar{
  private:
    // tree
    Tree tree;
    Node start_node, end_node;
    int target_node_count;
    World world;
    QuadOptimalControl oc;
  public:
    KinoRrtStar(World& in_world, Node& in_start_node, Node& in_end_node, int in_target_node_count, int interior_point_count );
    // n_nodes: number of nodes to add to tree, if 0 then stop after first solution
    void buildTree(int n_nodes);
    bool sampleSpace();
    void run();
    void buildTreeTillNodeCount();
    void buildTreeTillFirstSolution();
    void sampleNode();

    // generate a random point in world, only x,y,z are concerned
    // TODO test
    Node getRandomNode();
    // move node to target if its further than radius
    // TODO test
    void moveToVicinity(Node& node, Node& target, double radius);
    bool connectToGoal(Node& node);
    double getNeighnourRadius(){ return 0.5; };

    double sqr(double a) { return a*a; }
    double dist(Node& a, Node& b) { return sqrt(sqr(a.x-b.x) + sqr(a.y-b.y) + sqr(a.z-b.z)); }

};

#endif
