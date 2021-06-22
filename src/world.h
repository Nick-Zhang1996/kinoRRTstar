// describe world obstacles
#ifndef WORLD_H
#define WORLD_H

#include <iostream>
#include <list>
#include <iterator>
#include "common.h"

class World {
  private:
    int interior_point_n;
    list <Box> obstacles;

  public:

    //x,y,z, low, high
    double x_l,x_h,y_l,y_h,z_l,z_h;


    // constructor with x,y,z dimension 
    World(double dim_x_l, double dim_x_h, double dim_y_l, double dim_y_h, double dim_z_l, double dim_z_h) : x_l(dim_x_l), x_h(dim_x_h), y_l(dim_y_l), y_h(dim_y_h), z_l(dim_z_l), z_h(dim_z_h) {}
    // add obstacle to world
    // [obstacle] can be  Box
    int addObstacle( Box& obstacle ) {
      obstacles.push_back(obstacle);
      cout << "obstacle " << obstacles.size() << endl;
      cout << " x: " << obstacle.x_l << " - " << obstacle.x_h << endl;
      cout << " y: " << obstacle.y_l << " - " << obstacle.y_h << endl;
      cout << " z: " << obstacle.z_l << " - " << obstacle.z_h << endl;
      
      return 0;
    }
    // check collision
    // check if point is in collision with registered obstacle
    // also check if point is outside of world
    // return: Trus if no collision
    bool checkNoCollision( Point& point){
      if (point.x < x_l || point.x > x_h || point.y < y_l || point.y > y_h ||point.z < z_l || point.z > z_h) { return false; }
      for(list <Box>::iterator it = obstacles.begin(); it != obstacles.end(); it++){
        if ( it->isInside(point) ){ return false; }
      }
      return true;
    }

    bool checkNoCollision( Node& node){
      Point point = {node.x, node.y, node.z};
      return checkNoCollision( point );
    }
    
    bool checkNoCollision( double x, double y, double z){
      Point point = {x, y, z};
      return checkNoCollision( point );
    }

    void setInteriorPointCount(int val){ interior_point_n = val; }
    // TODO test
    bool checkNoPathCollision(double* pos_buffer);

};



#endif
