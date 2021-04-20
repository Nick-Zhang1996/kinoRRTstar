#include "tree.h"

Tree::Tree(Node& root){
  // is this valid?
  tree.push_back(root);
  n_nodes++;
}

void Tree::addNode(Node& new_node, Node& parent){
  new_node.id = n_nodes;
  new_node.p_parent = &parent;
  tree.push_back(new_node);
  // new_node is copied into tree, need pointer to the new copied one in tree
  parent.p_children.push_back(&tree.back());
  n_nodes++;
}

list<Node*>& Tree::getChildrenP(Node& node){
  return node.p_children;

}

void Tree::transferChild(Node& child, Node& new_parent){
  assert (child.id != 0);
  bool flag_found_old_children = false;
  for (auto pp_old_children = child.p_parent->p_children.begin();pp_old_children!=child.p_parent->p_children.end();pp_old_children++){
    if ((*pp_old_children) == &child){
      child.p_parent->p_children.erase(pp_old_children);
      flag_found_old_children = true;
      break;
    }
  }
  if (!flag_found_old_children){ throw err_cant_find_child(); }
  child.p_parent = &new_parent;
  new_parent.p_children.push_back(&child);
}

bool Tree::isEnd(Node& node){
  return node.is_end;
}

// TODO quadtree
list<Node*> Tree::getNeighbour(Node& node, double radius){
  list<Node*> res;
  for (auto i=tree.begin(); i!=tree.end(); i++){
    if (dist(node,*i) < radius){
      res.push_back(&(*i));
    }
  }
  return res;
}

Node& Tree::getClosest(Node& node){
  double min_dist = std::numeric_limits<double>::max();
  Node* min_node=NULL;
  for (auto i=tree.begin(); i!=tree.end(); i++){
    if (dist(node,*i) < min_dist){
      min_dist = dist(node,*i);
      min_node = &(*i);
    }
  }
  assert (min_node != NULL);
  return *min_node;
}

void Tree::updateCost(Node& node, double cost_delta){
  node.cost += cost_delta;
  for (auto child=node.p_children.begin(); child!=node.p_children.end(); child++){
    updateCost(**child,cost_delta);
  }
}

void Tree::printInfo(){
  int node_id = 0;
  for (auto node=tree.begin(); node!=tree.end(); node++){
    cout << "node_id : " << node_id << "\n";
    if (node->id == 0){
      cout << "root node \n";
    } else {
      cout << "parent_id : " << node->p_parent->id << "\n";
    }
    cout << "children: ";
    for (auto child=node->p_children.begin(); child!=node->p_children.end(); child++){
      cout << (*child)->id << ", ";
    }
    cout << "\n \n";
    node_id++;
  }
  cout << " ---------------- \n";


  

}
