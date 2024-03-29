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
  tree.addNode(node1, 0);
  tree.printInfo();

  Node node2;
  node2.x = 2;
  node2.cost = 2.0;
  node2.is_end = false;
  tree.addNode(node2, 1);
  tree.printInfo();

  Node node3;
  node3.x = 3;
  node3.cost = 3.0;
  node3.is_end = true;
  tree.addNode(node3, tree.node(2).id_parent);
  tree.printInfo();


  // 0->1 , 1->2,3
  auto children = tree.getChildrenId(1);
  for (auto i=children.begin(); i!=children.end(); i++){
    cout << *i << ", ";
  }
  cout << "\n";
  // expect 2,3

  auto neighbour = tree.getNeighbourId(2,1.4);
  for (auto i=neighbour.begin(); i!=neighbour.end(); i++){
    cout << *i << ", ";
  }
  cout << "\n";
  // expect 1,3

  auto closest = tree.getClosest(3);
  cout << closest.x << "\n";
  // expect 2
  

  tree.updateCost(node1,1.0);
  cout << tree.node(3).cost << "\n";
  // expect 4.0
  cout << tree.node(2).cost << "\n";
  // expect 3.0
  cout << tree.node(0).cost << "\n";
  // expect 0.0
  //
  tree.transferChild(3, 0);
  // 0->1,3  1->2
  auto children2 = tree.getChildrenId(root_node);
  for (auto i=children2.begin(); i!=children2.end(); i++){
    cout << tree.node(*i).id << ", ";
  }
  cout << "\n";
  // expect 1,3
  //
  cout << tree.node(3).id_parent << "\n";
  // expect 0
  auto children3 = tree.getChildrenId(1);
  for (auto i=children3.begin(); i!=children3.end(); i++){
    cout << (*i) << ", ";
  }
  cout << "\n";
  // expect 2

}
