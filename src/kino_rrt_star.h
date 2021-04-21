// kinodynamic RRT* implementation
#ifndef KINO_RRT_STAR_H
#define KINO_RRT_STAR_H

#include "tree.h"
#include <iostream>
#include "world.h"
using std::cout;

class KinoRrtStar{
  private:
    // tree
    Tree tree;
    Node start_node, end_node;
    int target_node_count;
    World world;
  public:
    // n_nodes: number of nodes to add to tree, if 0 then stop after first solution
    void buildTree(int n_nodes);
    bool sampleSpace();
    bool rewire();



};

#endif
