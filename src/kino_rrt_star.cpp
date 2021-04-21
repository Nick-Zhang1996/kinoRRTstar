#include "kino_rrt_star.h"

KinoRrtStar::KinoRrtStar(World& in_world, Node& in_start_node, Node& in_end_node, int in_target_node_count, int interior_point_count ) : 
    tree(in_start_node),
    start_node(in_start_node),
    end_node(in_end_node),
    target_node_count(in_target_node_count),
    oc(interior_point_count),
    world(in_world) {
    std::srand(std::time(nullptr)); // use current time as seed for random generator
    world.setInteriorPointCount(interior_point_count);
    cout << "Using only position coordinate for start/end \n";
}

void KinoRrtStar::run(){
  // build tree
  if (target_node_count <= 0){
    buildTreeTillFirstSolution();
  } else {
    buildTreeTillNodeCount();
  }

  if (tree.getSolutionCount() > 0){
    cout << "found " << tree.getSolutionCount() << " solutions \n";
  }

}

void KinoRrtStar::buildTreeTillFirstSolution(){
  while (tree.getSolutionCount() == 0 && tree.getNodeCount() < 10000) {
    sampleNode();
  }
}

void KinoRrtStar::buildTreeTillNodeCount(){
  // keep sampling node until terminal condition is met
  //
  while (tree.getNodeCount() < target_node_count ){
    sampleNode();
  }

}

void KinoRrtStar::sampleNode(){
  // sample node
  Node new_node = getRandomNode();

  Node& closest_node = tree.getClosest(new_node);
  // move to closest existing node
  // NOTE in place operation
  // TODO calculate radius properly
  double radius = getNeighnourRadius();
  moveToVicinity(new_node, closest_node, radius);

  //   check node collision
  if (!world.checkNoCollision(new_node)){ return; }

  // find traj
  double t = oc.timePartialFinalState(closest_node, new_node);
  double* interior_pos = oc.interiorPositionPartialFinalState(t, closest_node, new_node);
  //   check traj collision
  if (!world.checkNoPathCollision(interior_pos)){ return; }
  // search for better parent (lower segment cost)
  list<int> id_parent_candidates = tree.getNeighbourId(new_node, radius);
  int id_best_parent = closest_node.id;
  double lowest_cost = closest_node.cost + oc.costPartialFreeFinalState(t, closest_node, new_node);
  
  for (auto i=id_parent_candidates.begin(); i!=id_parent_candidates.end(); i++){
    Node& candidate_node = tree.node(*i);
    double this_t = oc.timePartialFinalState(candidate_node, new_node);
    double cost = candidate_node.cost + oc.costPartialFreeFinalState(this_t, candidate_node, new_node);
    if (cost < lowest_cost){
      lowest_cost = cost;
      id_best_parent = *i;
      t = this_t;
    }
  }

  // TODO set full state in new_node
  oc.setFullStatePartialFinalState(t, tree.node(id_best_parent), new_node);
  
  // check if new node is close to finish
  if (connectToGoal(new_node)){
    new_node.is_end = true;
    // solution count will be updated when node is added in tree
  }
  // add to tree
  tree.addNode(new_node,id_best_parent);

  //
  // rewire
  // can new node be a parent to neighbour nodes?
  
  // If yes then update cost of neighbour's children

}

Node KinoRrtStar::getRandomNode(){
  Node node;
  node.x = (double)rand()/(RAND_MAX + 1u) * world.getXSize();
  node.y = (double)rand()/(RAND_MAX + 1u) * world.getYSize();
  node.z = (double)rand()/(RAND_MAX + 1u) * world.getZSize();
  return node;
}

void KinoRrtStar::moveToVicinity(Node& node, Node& target, double radius){
  assert (radius > 0);
  double ori_dist = dist(node, target);
  if ( ori_dist < radius ) { return; }
  double scale = radius / ori_dist;
  node.x += (target.x - node.x) * scale;
  node.y += (target.y - node.y) * scale;
  node.z += (target.z - node.z) * scale;
  assert ( dist(node,target) < radius + 0.0001 );
  return;
}

bool KinoRrtStar::connectToGoal(Node& node){
  // find traj
  double t = oc.timePartialFinalState(node, end_node);
  double* interior_pos = oc.interiorPositionPartialFinalState(t, node, end_node);
  //   check traj collision
  if (world.checkNoPathCollision(interior_pos)){ 
    return true; 
  } else {
    return false;
  }
}

