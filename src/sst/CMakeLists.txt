include_directories(SYSTEM ${Boost_INCLULDE_DIR})

find_package(Boost 1.58 QUIET REQUIRED COMPONENTS serialization filesystem system program_options)
find_package(ompl)
find_package(Eigen3 3.3 REQUIRED)

add_executable(RigidBodyPlanning RigidBodyPlanning.cpp)
target_include_directories(RigidBodyPlanning PRIVATE ${OMPL_INCLUDE_DIRS})

target_link_libraries(RigidBodyPlanning
  ompl
  Eigen3::Eigen
  ${Boost_FILESYSTEM_LIBRARY}
  ${Boost_SYSTEM_LIBRARY}
  ${Boost_PROGRAM_OPTIONS_LIBRARY})
