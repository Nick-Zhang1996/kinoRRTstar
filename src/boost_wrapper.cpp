#include <boost/python.hpp>
#include "common.h"
#include "world.h"
#include "kino_rrt_star.h"

using namespace boost::python;
BOOST_PYTHON_MODULE(kinoRRT){

  class_<Box>("Box", init<double, double, double, double, double, double>())
    .def_readonly("x_l", &Box::x_l)
    .def_readonly("x_h", &Box::x_h)
    .def_readonly("y_l", &Box::y_l)
    .def_readonly("y_h", &Box::y_h)
    .def_readonly("z_l", &Box::z_l)
    .def_readonly("z_h", &Box::z_h)
    ;

  class_<Node>("Node",init<double, double, double>());

  class_<Waypoint>("Waypoint")
    .def_readonly("t", &Waypoint::t)
    .def_readonly("x", &Waypoint::x)
    .def_readonly("y", &Waypoint::y)
    .def_readonly("z", &Waypoint::z)
    .def_readonly("vx", &Waypoint::vx)
    .def_readonly("vy", &Waypoint::vy)
    .def_readonly("vz", &Waypoint::vz)
    .def_readonly("ax", &Waypoint::ax)
    .def_readonly("ay", &Waypoint::ay)
    .def_readonly("az", &Waypoint::az)
    .def_readonly("valid", &Waypoint::valid)
    ;

  class_<World>("World", init<double, double, double, double, double, double>())
    .def("addObstacle", &World::addObstacle)
    .def("setInteriorPointCount", &World::setInteriorPointCount)
    .def_readonly("x_l", &World::x_l)
    .def_readonly("x_h", &World::x_h)
    .def_readonly("y_l", &World::y_l)
    .def_readonly("y_h", &World::y_h)
    .def_readonly("z_l", &World::z_l)
    .def_readonly("z_h", &World::z_h)
    ;

  class_<KinoRrtStar>("KinoRrtStar",init<World&, Node&, Node&, int, int>())
    .def("run", &KinoRrtStar::run)
    .def("prepareSolution", &KinoRrtStar::prepareSolution)
    .def("prepareSolutionWithInteriorPoints", &KinoRrtStar::prepareSolutionWithInteriorPoints)
    .def("getNextWaypoint", &KinoRrtStar::getNextWaypoint)
    .def("getTrajectoryTime", &KinoRrtStar::getTrajectoryTime)
    .def("getTrajectory", &KinoRrtStar::getTrajectory)
    .def_readonly("waypoint", &KinoRrtStar::waypoint)
    ;
}

