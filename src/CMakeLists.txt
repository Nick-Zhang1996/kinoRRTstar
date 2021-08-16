cmake_minimum_required(VERSION 3.10)
project(kinoRRT VERSION 1.0)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-undef")


find_package(PythonLibs 3.6 REQUIRED)
#find_package(Boost 1.58 QUIET REQUIRED COMPONENTS serialization filesystem system program_options python3)
find_package(Boost 1.58 QUIET REQUIRED COMPONENTS serialization filesystem system program_options python)
find_package(ompl)

include_directories(SYSTEM ${Boost_INCLULDE_DIR})
include_directories(${PROJECT_SOURCE_DIR} ${PROJECT_SOURCE_DIR}/sst)
link_directories($ENV{HOME}/.local/lib)

# main
add_executable(main main.cpp kino_rrt_star.cpp quad_optimal_control.cpp world.cpp tree.cpp)
target_include_directories(main PRIVATE
  ${PYTHON_INCLUDE_DIRS})
target_link_libraries(main
  ${PYTHON_LIBRARIES}
  ${Boost_LIBRARIES})
message(STATUS "[main]Boost_LIBRARIES: ${Boost_LIBRARIES}")

# boost library kinoRRT
add_library(kinoRRT SHARED boost_wrapper.cpp world.cpp kino_rrt_star.cpp quad_optimal_control.cpp tree.cpp sst/sst.cpp) 
set_target_properties(kinoRRT PROPERTIES PREFIX "")
target_include_directories(kinoRRT PRIVATE
  ${PYTHON_INCLUDE_DIRS}
  ${OMPL_INCLUDE_DIRS})

target_link_libraries(kinoRRT 
  ompl
  ${PYTHON_LIBRARIES}
  ${Boost_LIBRARIES})

add_subdirectory(tests)
add_subdirectory(sst)
