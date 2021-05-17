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

    double x,y,z;


    // constructor with x,y,z dimension 
    World( double dim_x, double dim_y, double dim_z ) : x(dim_x), y(dim_y), z(dim_z) {}
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
      if (point.x < 0.0 || point.x > x || point.y < 0.0 || point.y > y ||point.z < 0.0 || point.z > z) { return false; }
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
    double getXSize() { return x; }
    double getYSize() { return y; }
    double getZSize() { return z; }

   

};



#endif
