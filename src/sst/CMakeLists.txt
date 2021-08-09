add_executable(RigidBodyPlanning RigidBodyPlanning.cpp)

target_link_libraries(RigidBodyPlanning
  ompl
  Eigen3::Eigen
  ${Boost_FILESYSTEM_LIBRARY}
  ${Boost_SYSTEM_LIBRARY}
  ${Boost_PROGRAM_OPTIONS_LIBRARY})
