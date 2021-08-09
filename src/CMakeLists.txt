cmake_minimum_required(VERSION 3.10)
project(kinoRrt VERSION 1.0)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)
#include_directories(SYSTEM ${Boost_INCLULDE_DIR})

find_package(PythonLibs 3.6 EXACT REQUIRED)
find_package(Boost 1.58 QUIET REQUIRED COMPONENTS serialization filesystem system program_options python)

include_directories(${PROJECT_SOURCE_DIR})

add_executable(main main.cpp kino_rrt_star.cpp quad_optimal_control.cpp world.cpp tree.cpp)

add_library(kinoRrt MODULE boost_wrapper world kino_rrt_star quad_optimal_control tree) 
target_include_directories(kinoRrt PRIVATE
  /usr/include/python3.9)

target_link_libraries(kinoRrt 
  ${PYTHON_LIBRARIES}
  ${Boost_LIBRARIES})


#  ${Boost_FILESYSTEM_LIBRARY}
#  ${Boost_SYSTEM_LIBRARY}
#  ${Boost_PROGRAM_OPTIONS_LIBRARY})

add_subdirectory(tests)
