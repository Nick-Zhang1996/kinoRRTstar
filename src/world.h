// describe world obstacles
#ifndef WORLD_H
#define WORLD_H

#include <iostream>
#include <list>
#include <iterator>
using namespace std;

struct Point {
  float x,y,z;
};

class Body {
  public:
    virtual bool isInside( Point point) const = 0;
};

class Box : public Body {
  private:
    float x_l, x_h, y_l, y_h, z_l, z_h;
  public:
    Box(float x_l, float x_h, float y_l, float y_h, float z_l, float z_h) :
      x_l(x_l), x_h(x_h), y_l(y_l), y_h(y_h), z_l(z_l), z_h(z_h) {}
    bool isInside( Point point) const {
      if (x_l < point.x and point.x < x_h
          and y_l < point.y and point.y < y_h
          and z_l < point.z and point.z < z_h){
        return true;
      } else { return false; }

    }

};

class World {
  private:
    float x,y,z;
    list <Body*> obstacles;

  public:
    // constructor with x,y,z dimension 
    World( float dim_x, float dim_y, float dim_z ) : x(dim_x), y(dim_y), z(dim_z) {}
    // add obstacle to world
    // [obstacle] can be Sphere, Box, or other subclass of Body
    int addObstacle( Body* p_obstacle ) {
      obstacles.push_back(p_obstacle);
      return 0;
    }
    // check collision
    // check if point is in collision with registered obstacle
    // return: Trus if no collision
    bool checkNoCollision( Point point){
      for(list <Body*>::iterator it = obstacles.begin(); it != obstacles.end(); it++){
        if ( (*it)->isInside(point) ){ return false; }
      }
      return true;
    }

};



#endif
