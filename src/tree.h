// tree structure
//
#ifndef TREE_H
#define TREE_H
#include <list>
#include <vector>
using std::list;
using std::vector;

// represent one node
struct Node {
    double x,y,z,vx,vy,vz,ax,ay,az;
    list<Node*> children;
    Node* parent;
    // cost from start
    double cost;
    // connected to end-goal
    bool is_end;
};

class Tree{
  private:
    // node count
    int n_nodes = 0;
    vector<Node> tree;
    // euclidean distance
    double dist(Node& a, Node& b);

  public:
    Tree(Node& root);
    int getSolutionCount();
    int getNodeCount();

    void addNode(Node& new_node, Node& parent);
    list<Node*>& getChildren(Node& node);
    void transferChild(Node& child, Node& new_parent);
    bool isEnd(Node& node);

    // get neighbour list within radius
    list<Node*> getNeighbour(Node& node, double radius);
    // get closest node
    Node& getClosest(Node& node);
    void updateCost(Node& node, double cost_delta);


    

};


#endif
