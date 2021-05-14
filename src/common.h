#ifndef COMMON_H_
#define COMMON_H_

#include <iostream>
#include <list>
#include <vector>
using std::list;
using std::vector;
using std::cout;
using std::endl;

struct Waypoint {
    double x,y,z,vx,vy,vz,ax,ay,az,dummy;
    bool valid;
};

// represent one node
struct Node {
    double x,y,z,vx,vy,vz,ax,ay,az;
    list<int> id_children;
    int id_parent;
    // cost from start
    double cost;
    // connected to end-goal
    bool is_end;
    // may not be necessary
    int id;
    Node():is_end(false){};
    Node(double xx, double yy, double zz):is_end(false),x(xx),y(yy),z(zz){};
};

struct Point {
  double x,y,z;
};

class Body {
  public:
    virtual bool isInside( Point point) const = 0;
};

class Box : public Body {
  private:
    double x_l, x_h, y_l, y_h, z_l, z_h;
  public:
    Box(double x_l, double x_h, double y_l, double y_h, double z_l, double z_h) :
      x_l(x_l), x_h(x_h), y_l(y_l), y_h(y_h), z_l(z_l), z_h(z_h) {}
    bool isInside( Point point) const {
      if (x_l < point.x and point.x < x_h
          and y_l < point.y and point.y < y_h
          and z_l < point.z and point.z < z_h){
        return true;
      } else { return false; }

    }

};
#endif
