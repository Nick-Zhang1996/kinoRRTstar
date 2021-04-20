// tree structure
//
#ifndef TREE_H
#define TREE_H
#include <iostream>
#include <list>
#include <vector>
#include <limits>
#include <assert.h>
#include <math.h>
using std::cout;
using std::list;
using std::vector;
using std::sqrt;

// represent one node
struct Node {
    double x,y,z,vx,vy,vz,ax,ay,az;
    list<int> id_children;
    int id_parent;
    // cost from start
    double cost;
    // connected to end-goal
    bool is_end;
    // may not be necessary
    int id;
};
class err_cant_find_child : public std::exception{};

class Tree{
  private:
    // node count
    int n_nodes = 0;
    vector<Node> tree;
    inline double sqr(double x){return x*x;};
    // euclidean distance
    inline double dist(Node& a, Node& b){return sqrt(sqr(a.x-b.x)+sqr(a.y-b.y)+sqr(a.z-b.z));}
    // nodes that connect to goal
    list<int> id_end_nodes;
    int n_solutions = 0;

  public:
    Tree(Node& root);
    int getSolutionCount(){ return n_solutions; }
    int getNodeCount(){ return n_nodes; }

    void addNode(Node& new_node, Node& parent);
    void addNode(Node& new_node, int id_parent);

    list<int>& getChildrenId(Node& node);
    list<int>& getChildrenId(int id_node);

    void transferChild(int id_child, int id_new_parent);

    bool isEnd(Node& node);
    bool isEnd(int id_node);

    // get neighbour list within radius
    list<int> getNeighbourId(Node& node, double radius);
    list<int> getNeighbourId(int id_node, double radius);

    // get closest node
    int getClosestId(int id_node);
    Node& getClosest(int id_node);

    void updateCost(Node& node, double cost_delta);
    void updateCost(int id_node, double cost_delta);

    Node& node(int id){ return tree[id]; }

    void printInfo();
};


#endif
