#include "kino_rrt_star.h"

KinoRrtStar::KinoRrtStar(World& in_world, Node& in_start_node, Node& in_end_node, int in_target_node_count, int in_interior_point_count ) : 
    tree(in_start_node),
    start_node(in_start_node),
    end_node(in_end_node),
    target_node_count(in_target_node_count),
    world(in_world),
    oc(in_interior_point_count),
    interior_point_count(in_interior_point_count),
    rewire_count(0)
{
    overall_lowest_cost = std::numeric_limits<double>::max();
    std::srand(std::time(nullptr)); // use current time as seed for random generator
    world.setInteriorPointCount(in_interior_point_count);
    cout << "Using only position coordinate for start/end \n";
}

void KinoRrtStar::run(){
  cout << "KinoRrtStar.run()" << endl;
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

  oc.printTotalTime();

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
  cout << "rewire: " << rewire_count << "\n";

}

void KinoRrtStar::buildTreeTillFirstSolution(){
  while (tree.getSolutionCount() == 0 && tree.getNodeCount() < 10000) {
    sampleNode();
    if (tree.getNodeCount() % 100 == 0){
      cout << tree.getNodeCount() << "\n";
    }
  }
}

void KinoRrtStar::buildTreeTillNodeCount(){
  // keep sampling node until terminal condition is met
  //
  while (tree.getNodeCount() < target_node_count ){
    sampleNode();
    if (tree.getNodeCount() % 100 == 0){
      cout << tree.getNodeCount() << "\n";
    }
  }

}

void KinoRrtStar::sampleNode(){
  // sample node
  Node new_node = getRandomNode();

  Node& closest_node = tree.getClosest(new_node);
  double radius = getNeighnourRadius();
  // move to closest existing node
  // NOTE in place operation
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
  
  neighbour_total_count += id_parent_candidates.size();
  neighbour_count += 1;
  for (auto i=id_parent_candidates.begin(); i!=id_parent_candidates.end(); i++){
    Node& candidate_node = tree.node(*i);
    double this_t = oc.timePartialFinalState(candidate_node, new_node);
    double cost = candidate_node.cost + oc.costPartialFreeFinalState(this_t, candidate_node, new_node);
    if (cost < lowest_cost){
      // check for path collision
      double* interior_pos = oc.interiorPositionPartialFinalState(this_t, candidate_node, new_node);
      if (!world.checkNoPathCollision(interior_pos)){ continue; }
      // no collision
      lowest_cost = cost;
      id_best_parent = *i;
      t = this_t;
    }
  }

  // up to now only x,y,z in new_node is used
  // now we know the final parent
  // set full state in new_node
  // this sets vxyz and axyz
  oc.setFullStatePartialFinalState(t, tree.node(id_best_parent), new_node);
  new_node.cost = lowest_cost;
  
  // check if new node is close to finish
  if (connectToGoal(new_node)){
    new_node.is_end = true;
    // solution count will be updated when node is added in tree
    if (new_node.cost < overall_lowest_cost){
      overall_lowest_cost = new_node.cost;
      overall_lowest_cost_id = tree.getNodeCount();
      //cout << "cost: " << overall_lowest_cost << "\n";
      //cout << "avg neighbor count: " << (double) neighbour_total_count / neighbour_count << "\n";
      neighbour_total_count = 0;
      neighbour_count = 0;
    }
  }
  // add to tree
  // sets id_parent, id
  int id_new_node = tree.addNode(new_node,id_best_parent);
  // rewire
  list<int> id_neighbour = tree.getNeighbourId(id_new_node, radius);
  // can new node be a parent to neighbour nodes?
  for (auto i=id_neighbour.begin(); i!=id_neighbour.end(); i++){
    try {
      double new_t = oc.time(new_node, tree.node(*i));
      double new_cost = new_node.cost + oc.cost(new_t, new_node, tree.node(*i));
      // If yes then let new_node be new parent
      // and update cost of neighbour and its descendents
      if (new_cost > tree.node(*i).cost){ continue; }
      double* interior_pos = oc.interiorPosition(new_t, new_node, tree.node(*i));
      //   check traj collision
      if (!world.checkNoPathCollision(interior_pos)){ continue; }
      tree.transferChild(*i, id_new_node);
      tree.updateCost(*i, new_cost - tree.node(*i).cost);
      //cout << "rewiring... \n";
      rewire_count++;
    } catch (err_NoSolution&){ continue; }
  }
  

}

Node KinoRrtStar::getRandomNode(){
  Node node;
  node.x = (double)rand()/(RAND_MAX + 1u) * (world.x_h - world.x_l) + world.x_l;
  node.y = (double)rand()/(RAND_MAX + 1u) * (world.y_h - world.y_l) + world.y_l;
  node.z = (double)rand()/(RAND_MAX + 1u) * (world.z_h - world.z_l) + world.z_l;
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
  if (dist(node, end_node) > getNeighnourRadius()){ return false; }
  // find traj
  double t = oc.time(node, end_node);
  double* interior_pos = oc.interiorPosition(t, node, end_node);
  //   check traj collision
  if (world.checkNoPathCollision(interior_pos)){ 
    return true; 
  } else {
    return false;
  }
}

// per original implementation from MATLAB
double KinoRrtStar::getNeighnourRadius(){
  double nun = (double) tree.getNodeCount();
  double ner = 40.0 * pow( ( log(nun + 1) / nun ), 1.0/3.0 );
  return min(ner,4.0);

}

int KinoRrtStar::prepareSolutionWithInteriorPoints(){
  // first retrieve the optimal path, this is done in reverse order from goal to start
  // key_waypoints_reversed includes index of lowest cost node connected to end node ... root node(start) 
  list<int> key_waypoints_reversed;
  int this_node_id = overall_lowest_cost_id;
  assert (tree.node(this_node_id).is_end);

  key_waypoints_reversed.push_back(this_node_id);
  while (this_node_id != 0){
    Node& node = tree.node(this_node_id);
    this_node_id = node.id_parent;
    key_waypoints_reversed.push_back(this_node_id);

  }


  /*
  // checking items in key_waypoints_reversed
  cout << " showing items in key_waypoints_reversed " << endl;
  for (auto i=key_waypoints_reversed.begin(); i!=key_waypoints_reversed.end(); i++){
    Node& node = tree.node(*i);
    cout << "id: " << *i << " " << node.x << ","  << node.y << "," << node.z << endl;
  }
  cout << "done" << endl;
  */

  // now reverse iterate on the key waypoints, adding in interior points
  double waypoint_t = 0.0;
  waypoints.clear();
  // start from second
  list<int>::reverse_iterator i=key_waypoints_reversed.rbegin();
  ++i;

  for ( ; i != key_waypoints_reversed.rend(); ++i){

    Node& late_node = tree.node(*i);
    Node& early_node = tree.node(late_node.id_parent);

    // add early_node
    Waypoint p;
    p.valid = true;
    p.t = waypoint_t;
    p.x = early_node.x;
    p.y = early_node.y;
    p.z = early_node.z;
    waypoints.push_back(p);
    // DEBUG
    //cout << "examining path from node x:" << p.x << ", y:" << p.y << ", z:" << p.z << endl;
    //cout << "                 to node x:" << late_node.x << ", y:" << late_node.y << ", z:" << late_node.z << endl;

    double section_time = oc.timePartialFinalState(early_node, late_node);
    double* pos_buffer = oc.interiorPositionPartialFinalState(section_time, early_node, late_node);
    double step_time = section_time / (interior_point_count+1);

    waypoint_t += step_time;
    for (int i=0; i<interior_point_count; i++){
      Waypoint p;
      p.valid = true;
      p.t = waypoint_t;
      p.x = *(pos_buffer + i*3 + 0);
      p.y = *(pos_buffer + i*3 + 1);
      p.z = *(pos_buffer + i*3 + 2);
      waypoints.push_back(p);
      assert (world.checkNoCollision(p.x,p.y,p.z));
      waypoint_t += step_time;
    }

    /*
    if (tree.isRoot(late_node)){
      // now add the very last node
      Waypoint p;
      p.valid = true;
      p.t = waypoint_t;
      p.x = late_node.x;
      p.y = late_node.y;
      p.z = late_node.z;
      waypoints.push_back(p);
    }
    */

  }

  // now add in last waypoint to goal_node
  Node& late_node = end_node;
  Node& early_node = tree.node(overall_lowest_cost_id);

  Waypoint p;
  p.valid = true;
  p.t = waypoint_t;
  p.x = early_node.x;
  p.y = early_node.y;
  p.z = early_node.z;
  waypoints.push_back(p);

  double section_time = oc.timePartialFinalState(early_node, late_node);
  double* pos_buffer = oc.interiorPositionPartialFinalState(section_time, early_node, late_node);
  double step_time = section_time / (interior_point_count+1);

  waypoint_t += step_time;
  for (int i=0; i<interior_point_count; i++){
    Waypoint p;
    p.valid = true;
    p.t = waypoint_t;
    p.x = *(pos_buffer + i*3 + 0);
    p.y = *(pos_buffer + i*3 + 1);
    p.z = *(pos_buffer + i*3 + 2);
    waypoints.push_back(p);
    assert (world.checkNoCollision(p.x,p.y,p.z));
    waypoint_t += step_time;
  }

  // add in end node (goal)
  p.t = waypoint_t;
  p.x = end_node.x;
  p.y = end_node.y;
  p.z = end_node.z;
  waypoints.push_back(p);

  waypoints_iter = waypoints.begin();
  return waypoints.size();

}


int KinoRrtStar::prepareSolution(){
  // first retrieve the optimal path, this is done in reverse order from goal to start
  // key_waypoints_reversed includes index of lowest cost node connected to end node ... root node(start) 
  list<int> key_waypoints_reversed;
  int this_node_id = overall_lowest_cost_id;
  if (!tree.node(this_node_id).is_end){
    return 0;
  }

  key_waypoints_reversed.push_back(this_node_id);
  while (this_node_id != 0){
    Node& node = tree.node(this_node_id);
    this_node_id = node.id_parent;
    key_waypoints_reversed.push_back(this_node_id);
  }

  // now reverse iterate on the key waypoints construct list and add timestamp
  double waypoint_t = 0.0;
  waypoints.clear();
  list<int>::reverse_iterator i=key_waypoints_reversed.rbegin();
  ++i;
  for ( ; i != key_waypoints_reversed.rend(); ++i){
    Node& late_node = tree.node(*i);
    Node& early_node = tree.node(late_node.id_parent);
    double section_time = oc.timePartialFinalState(early_node, late_node);
    // add early_node
    Waypoint p;
    p.t = waypoint_t;
    p = early_node;

    waypoints.push_back(p);
    waypoint_t += section_time;
    // DEBUG verify no collision
    if (!world.checkNoCollision(p.x,p.y,p.z)){
      cout << "[prepareSolution()] found invalid waypoint()" << endl;
    }
  }

  // now add in last waypoint to goal_node
  Node& late_node = end_node;
  Node& early_node = tree.node(overall_lowest_cost_id);
  double section_time = oc.timePartialFinalState(early_node, late_node);

  Waypoint p;
  p.t = waypoint_t;
  p = early_node;
  waypoints.push_back(p);
  waypoint_t += section_time;

  // add in end node (goal)
  p.t = waypoint_t;
  p = end_node;
  waypoints.push_back(p);

  waypoints_iter = waypoints.begin();
  return waypoints.size();
}

Waypoint KinoRrtStar::getNextWaypoint(){
  if (waypoints_iter == waypoints.end()){
    waypoint.valid = false;
  } else {
    // copy to persistent variable so we don't return a local var
    auto retval = *waypoints_iter;
    waypoint.valid = true;
    waypoint = retval;
    waypoints_iter++;
  }
  return waypoint;
}

// this requires prepareSolution to be called
// get trajectory time
double KinoRrtStar::getTrajectoryTime(){
  return waypoints.back().t;

}

// this invalidates Waypoint related functions
Waypoint KinoRrtStar::getTrajectory(double t){
  try {
    Waypoint early_waypoint,late_waypoint;
    for (auto iter = waypoints.begin(); iter != waypoints.end(); iter++){
      if ( (*iter).t <= t) {
        early_waypoint = *iter;
        continue;
      }
      late_waypoint = *(iter);
      if (iter == waypoints.begin()){
        // if t < (*iter).t, then late_waypoint = first waypoint and early_waypoint is never set
        early_waypoint =  *(iter);
        late_waypoint = *(++iter);
      }
      break;
    }



    //cout << "assigning late waypoint" << endl;
    //cout << "done assigning late waypoint" << endl;

    //double section_time = oc.time(early_waypoint, late_waypoint);
    double section_time = late_waypoint.t - early_waypoint.t;
    double dt = t - early_waypoint.t;
    //cout << "dt = " << dt << endl;
    //cout << "section_t = " << section_time << endl;
    //cout << " late_waypoint.t = " << late_waypoint.t << endl;
    //cout << " late-early = " << late_waypoint.t-early_waypoint.t << endl;

    //assert (dt > -0.0001);
    //assert (dt < section_time-0.0001);
    //
    //cout << "early: " << early_waypoint.t << "," << early_waypoint.x << "," << early_waypoint.y << "," << early_waypoint.z << endl;
    //cout << "late: " << late_waypoint.t << "," << late_waypoint.x << "," << late_waypoint.y << "," << late_waypoint.z << endl;
    //cout << "state(" << dt << "," << section_time << "," << early_waypoint.x << "," << early_waypoint.y << "," << early_waypoint.z << "," << early_waypoint.vx << "," << early_waypoint.vy << "," << early_waypoint.vz << "," << early_waypoint.ax << "," << early_waypoint.ay << "," << early_waypoint.az << "," << late_waypoint.x << "," << late_waypoint.y << "," << late_waypoint.z << endl;<< "," << late_waypoint.vx << "," << late_waypoint.vy << "," << late_waypoint.vz << endl;<< "," << late_waypoint.ax << "," << late_waypoint.ay << "," << late_waypoint.az << ")" << endl;


    oc.state(dt, section_time, early_waypoint, late_waypoint);
    waypoint.t = t;
    waypoint.valid = true;
    waypoint.x = *(oc.buffer_interior_state + 0);
    waypoint.y = *(oc.buffer_interior_state + 1);
    waypoint.z = *(oc.buffer_interior_state + 2);
    waypoint.vx = *(oc.buffer_interior_state + 3);
    waypoint.vy = *(oc.buffer_interior_state + 4);
    waypoint.vz = *(oc.buffer_interior_state + 5);
    waypoint.ax = *(oc.buffer_interior_state + 6);
    waypoint.ay = *(oc.buffer_interior_state + 7);
    waypoint.az = *(oc.buffer_interior_state + 8);
    return waypoint;
  } catch ( const std::exception &e) {
    cout << "error: ";
    cout << e.what() << endl;
  } 

  return waypoint;

}
