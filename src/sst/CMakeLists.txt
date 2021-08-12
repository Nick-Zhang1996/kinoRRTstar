include_directories(SYSTEM ${Boost_INCLULDE_DIR})

add_executable(sst sst.cpp)
target_include_directories(sst PRIVATE 
  ${PYTHON_INCLUDE_DIRS}
  ${OMPL_INCLUDE_DIRS})
message(STATUS "ompl_include_dir: ${OMPL_INCLUDE_DIRS}")
message(STATUS "[sst]Boost_LIBRARIES: ${Boost_LIBRARIES}")

target_link_libraries(sst
  ompl
#  Eigen3::Eigen
  ${PYTHON_LIBRARIES}
  ${Boost_LIBRARIES})
#  ${Boost_FILESYSTEM_LIBRARY}
#  ${Boost_SYSTEM_LIBRARY}
#  ${Boost_PROGRAM_OPTIONS_LIBRARY})
