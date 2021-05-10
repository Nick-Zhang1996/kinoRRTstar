#include <boost/python.hpp>
#include "common.h"
#include "world.h"
#include "kino_rrt_star.h"

using namespace boost::python;
BOOST_PYTHON_MODULE(kinoRRT){

  class_<Box>("Box", init<double, double, double, double, double, double>());

  class_<Node>("Node",init<double, double, double>());

  class_<Waypoint>("Waypoint")
    .def_readonly("x", &Waypoint::x)
    .def_readonly("y", &Waypoint::y)
    .def_readonly("z", &Waypoint::z)
    .def_readonly("valid", &Waypoint::valid)
    ;

  class_<World>("World", init<double, double, double>())
    .def("addObstacle", &World::addObstacle)
    .def("setInteriorPointCount", &World::setInteriorPointCount)
    ;

  class_<KinoRrtStar>("KinoRrtStar",init<World&, Node&, Node&, int, int>())
    .def("run", &KinoRrtStar::run)
    .def("prepareSolution", &KinoRrtStar::prepareSolution)
    .def("getNextWaypoint", &KinoRrtStar::getNextWaypoint)
    ;
}

