#include "tree.h"
#include <iostream>
using std::cout;

int main(){
  Node root_node;
  root_node.x = 0;
  root_node.cost = 0.0;
  root_node.is_end = false;

  Tree tree(root_node);
  tree.printInfo();

  Node node1;
  node1.x = 1;
  node1.cost = 1.0;
  node1.is_end = false;
  tree.addNode(node1, root_node);
  tree.printInfo();

  Node node2;
  node2.x = 2;
  node2.cost = 2.0;
  node2.is_end = false;
  tree.addNode(node2, node1);
  tree.printInfo();

  Node node3;
  node3.x = 3;
  node3.cost = 3.0;
  node3.is_end = true;
  tree.addNode(tree.getParent(node2), node1);
  tree.printInfo();


  // 0->1 , 1->2,3
  auto children = tree.getChildrenP(node1);
  for (auto i=children.begin(); i!=children.end(); i++){
    cout << (*i)->x << ", ";
  }
  cout << "\n";
  // expect 2,3

  auto neighbour = tree.getNeighbour(node2,1.4);
  for (auto i=neighbour.begin(); i!=neighbour.end(); i++){
    cout << (*i)->x << ", ";
  }
  cout << "\n";
  // expect 1,3

  auto closest = tree.getClosest(node3);
  cout << closest.x << "\n";
  // expect 2
  

  tree.updateCost(node1,1.0);
  cout << node3.cost << "\n";
  // expect 4.0
  cout << node2.cost << "\n";
  // expect 3.0
  cout << root_node.cost << "\n";
  // expect 0.0
  //
  tree.transferChild(node3, root_node);
  children = tree.getChildrenP(root_node);
  for (auto i=children.begin(); i!=children.end(); i++){
    cout << (*i)->x << ", ";
  }
  cout << "\n";
  // expect 1,3
  cout << node3.p_parent->x << "\n";
  // expect 0
  children = tree.getChildrenP(node1);
  for (auto i=children.begin(); i!=children.end(); i++){
    cout << (*i)->x << ", ";
  }
  cout << "\n";
  // expect 2
    




}
