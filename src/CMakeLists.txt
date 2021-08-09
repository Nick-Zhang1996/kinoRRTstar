cmake_minimum_required(VERSION 3.10)
project(kinoRRT VERSION 1.0)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED True)
include_directories(SYSTEM ${Boost_INCLULDE_DIR})

find_package(PythonLibs 3.6 EXACT REQUIRED)
find_package(Boost 1.58 QUIET REQUIRED COMPONENTS serialization filesystem system program_options python3)

include_directories(${PROJECT_SOURCE_DIR})

add_executable(main main.cpp kino_rrt_star.cpp quad_optimal_control.cpp world.cpp tree.cpp)

add_library(kinoRRT SHARED boost_wrapper world kino_rrt_star quad_optimal_control tree) 
set_target_properties(kinoRRT PROPERTIES PREFIX "")
target_include_directories(kinoRRT PRIVATE
  ${PYTHON_INCLUDE_DIRS})

target_link_libraries(kinoRRT 
  ${PYTHON_LIBRARIES}
  ${Boost_LIBRARIES})

message(STATUS "PYTHON_LIBRARIES: ${PYTHON_LIBRARIES}")
message(STATUS "Boost_LIBRARIES: ${Boost_LIBRARIES}")


#  ${Boost_FILESYSTEM_LIBRARY}
#  ${Boost_SYSTEM_LIBRARY}
#  ${Boost_PROGRAM_OPTIONS_LIBRARY})

add_subdirectory(tests)
add_subdirectory(sst)
