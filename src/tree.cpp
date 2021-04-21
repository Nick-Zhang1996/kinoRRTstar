#include "tree.h"

Tree::Tree(Node& root){
  // is this valid?
  tree.push_back(root);
  n_nodes++;
}

void Tree::addNode(Node& new_node, Node& parent){
  addNode(new_node, parent.id);
}
void Tree::addNode(Node& new_node, int id_parent){
  assert (id_parent < n_nodes);
  new_node.id = n_nodes;
  new_node.id_parent = id_parent;
  tree.push_back(new_node);
  tree[id_parent].id_children.push_back(new_node.id);
  n_nodes++;
}

list<int>& Tree::getChildrenId(int id_node){
  assert (id_node < n_nodes);
  return tree[id_node].id_children;
}
list<int>& Tree::getChildrenId(Node& node){
  return getChildrenId(node.id);
}

void Tree::transferChild(int id_child, int id_new_parent){
  assert (tree[id_child].id != 0);
  bool flag_found_old_children = false;

  // remove child from old parent child list
  Node& old_parent = tree[tree[id_child].id_parent];
  for (auto p_id_old_children = old_parent.id_children.begin(); p_id_old_children != old_parent.id_children.end(); p_id_old_children++){
    if (*p_id_old_children == id_child){ 
      old_parent.id_children.erase(p_id_old_children); 
      flag_found_old_children = true;
      break;
    }
  }
  if (!flag_found_old_children){ throw err_cant_find_child(); }
  // set child parent
  tree[id_child].id_parent = id_new_parent;
  tree[id_new_parent].id_children.push_back(id_child);
}

bool Tree::isEnd(int id_node){
  return tree[id_node].is_end;
}
bool Tree::isEnd(Node& node){
  return isEnd(node.id);
}

// TODO quadtree
list<int> Tree::getNeighbourId(Node& node, double radius){
  getNeighbourId(node.id, radius);
}

list<int> Tree::getNeighbourId(int id_node, double radius){
  Node& node = tree[id_node];
  list<int> res;
  for (auto i=tree.begin(); i!=tree.end(); i++){
    if (i->id == id_node){ continue; }
    if (dist(node,*i) < radius){
      res.push_back(i->id);
    }
  }
  return res;
}

Node& Tree::getClosest(int id_node){
  return tree[getClosestId(id_node)];

}

int Tree::getClosestId(int id_node){
  double min_dist = std::numeric_limits<double>::max();
  int min_node_id= -1;
  Node& node = tree[id_node];
  for (auto i=tree.begin(); i!=tree.end(); i++){
    if (i->id == id_node){ continue; }
    if (dist(node,*i) < min_dist){
      min_dist = dist(node,*i);
      min_node_id = i->id;
    }
  }
  assert (min_node_id != -1);
  return min_node_id;
}

void Tree::updateCost(Node& node, double cost_delta){
  updateCost(node.id, cost_delta);
}

void Tree::updateCost(int id_node, double cost_delta){
  tree[id_node].cost += cost_delta;
  for (auto p_child_id=tree[id_node].id_children.begin(); p_child_id!=tree[id_node].id_children.end(); p_child_id++){
    updateCost(*p_child_id,cost_delta);
  }
}

void Tree::printInfo(){
  int node_id = 0;
  for (auto node=tree.begin(); node!=tree.end(); node++){
    cout << "node_id : " << node->id << "\n";
    if (node->id == 0){
      cout << "root node \n";
    } else {
      cout << "parent_id : " << tree[node->id_parent].id << "\n";
    }
    cout << "children: ";
    for (auto child=node->id_children.begin(); child!=node->id_children.end(); child++){
      cout << (*child) << ", ";
    }
    cout << "\n \n";
    node_id++;
  }
  cout << " ---------------- \n";


  

}
