#include "kino_rrt_star.h"

KinoRrtStar::KinoRrtStar(World& in_world, Node& in_start_node, Node& in_end_node, int in_target_node_count, int interior_point_count ) : 
    tree(in_start_node),
    start_node(in_start_node),
    end_node(in_end_node),
    target_node_count(in_target_node_count),
    oc(interior_point_count),
    world(in_world) {
    overall_lowest_cost = std::numeric_limits<double>::max();
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
    showResult();
  } else {
    cout << "no results find with " << target_node_count << "nodes \n";
  }


}

void KinoRrtStar::showResult(){
  int this_node_id = overall_lowest_cost_id;
  assert (tree.node(this_node_id).is_end);

  while (this_node_id != 0){
    Node& node = tree.node(this_node_id);
    cout << "id: " << node.id;
    cout << " x: " << node.x;
    cout << " y: " << node.y;
    cout << " z: " << node.z;
    cout << " cost: " << node.cost;
    cout << "\n";
    this_node_id = node.id_parent;
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
    if (tree.getNodeCount() % 1000 == 0){
      cout << tree.getNodeCount() << "\n";
    }
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
  // sets x,y,z
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
  // sets vxyz and axyz
  oc.setFullStatePartialFinalState(t, tree.node(id_best_parent), new_node);
  new_node.cost = lowest_cost;
  
  // check if new node is close to finish
  if (connectToGoal(new_node)){
    new_node.is_end = true;
    // solution count will be updated when node is added in tree
    if (new_node.cost < overall_lowest_cost){
      overall_lowest_cost = new_node.cost;
      overall_lowest_cost_id = tree.getNodeCount();
      cout << "cost: " << overall_lowest_cost << "\n";
    }
  }
  // add to tree
  // sets id_parent, id
  tree.addNode(new_node,id_best_parent);

  // rewire
  list<int> id_neighnour = tree.getNeighbourId(new_node, radius);
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
  node.x = target.x + (node.x - target.x) * scale;
  node.y = target.y + (node.y - target.y) * scale;
  node.z = target.z + (node.z - target.z) * scale;
  assert ( dist(node,target) < radius + 0.0001 );
  return;
}

bool KinoRrtStar::connectToGoal(Node& node){
  //TODO handle this properly
  if (dist(node, end_node) > 1.0){ return false; }
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

