// tree structure
//
#ifndef TREE_H
#define TREE_H
#include <list>
#include <vector>
#include <limits>
#include <assert.h>
#include <math.h>
using std::list;
using std::vector;
using std::sqrt;

// represent one node
struct Node {
    double x,y,z,vx,vy,vz,ax,ay,az;
    list<Node*> p_children;
    Node* p_parent;
    // cost from start
    double cost;
    // connected to end-goal
    bool is_end;
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
    vector<Node*> end_nodes;
    int n_solutions = 0;

  public:
    Tree(Node& root);
    int getSolutionCount(){ return n_solutions; }
    int getNodeCount(){ return n_nodes; }

    void addNode(Node& new_node, Node& parent);
    list<Node*>& getChildrenP(Node& node);
    void transferChild(Node& child, Node& new_parent);
    bool isEnd(Node& node);
    // get neighbour list within radius
    list<Node*> getNeighbour(Node& node, double radius);
    Node& getParent(Node& node){ return *node.p_parent; }

    // get closest node
    Node& getClosest(Node& node);
    void updateCost(Node& node, double cost_delta);
};


#endif
