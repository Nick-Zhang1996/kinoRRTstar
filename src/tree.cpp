#include "tree.h"

Tree::Tree(Node& root){
  // is this valid?
  tree.push_back(root);
  n_nodes++;
}

void Tree::addNode(Node& new_node, Node& parent){
  tree.push_back(new_node);
  new_node.p_parent = &parent;
  // new_node is copied into tree, need pointer to the one in tree
  parent.p_children.push_back(&tree.back());
}

list<Node*>& Tree::getChildrenP(Node& node){
  return node.p_children;

}

void Tree::transferChild(Node& child, Node& new_parent){
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
