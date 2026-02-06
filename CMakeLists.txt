# Check cmake version meets our requirements
cmake_minimum_required(VERSION 3.3)

# Declare project, which will configure compiler for us
project(FalaiseCimrmanPlugin)

# #Finds ROOT and Falaise packages and defines cmake variables
find_package(Falaise REQUIRED)
find_package(ROOT    REQUIRED)

include_directories(include ${Bayeux_INCLUDE_DIRS})

include_directories(${CMAKE_CURRENT_BINARY_DIR}/source)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/source)

list(APPEND Cimrman_HEADERS
  source/module/Cimrman.h
  source/reconstruction/Algos.h
  source/reconstruction/KinkFinder.h
  source/reconstruction/Sinogram.h
  source/reconstruction/TrajectoryBuilder.h
  source/geometry/Geometry.h
  source/visualization/Visu.h
  source/datamodel/Association.h
  source/datamodel/Cluster.h
  source/datamodel/DelayedCluster.h
  source/datamodel/Event.h
  source/datamodel/Likelihood.h
  source/datamodel/LinearFit.h
  source/datamodel/OMHit.h
  source/datamodel/Point.h
  source/datamodel/PreclusterSolution.h
  source/datamodel/Precluster.h
  source/datamodel/Solution.h
  source/datamodel/Track.h
  source/datamodel/TrackerHit.h
  source/datamodel/Trajectory.h
  )

list(APPEND Cimrman_SOURCES
  source/module/Cimrman.cpp
  source/reconstruction/Algos.cpp
  source/reconstruction/KinkFinder.cpp
  source/reconstruction/Sinogram.cpp
  source/reconstruction/TrajectoryBuilder.cpp
  source/visualization/Visu.cpp
  source/datamodel/Association.cpp
  source/datamodel/Cluster.cpp
  source/datamodel/DelayedCluster.cpp
  source/datamodel/Event.cpp
  source/datamodel/Likelihood.cpp
  source/datamodel/LinearFit.cpp
  source/datamodel/OMHit.cpp
  source/datamodel/Point.cpp
  source/datamodel/PreclusterSolution.cpp
  source/datamodel/Precluster.cpp
  source/datamodel/Solution.cpp
  source/datamodel/Track.cpp
  source/datamodel/TrackerHit.cpp
  source/datamodel/Trajectory.cpp
  )

add_library(Falaise_Cimrman SHARED
  ${Cimrman_SOURCES}
  ${Cimrman_HEADERS} )


target_link_libraries(Falaise_Cimrman Falaise::Falaise ${ROOT_LIBRARIES} ROOT::Geom)
# PUBLIC Falaise::Falaise FalaiseModule ${ROOT_LIBRARIES} ROOT::Geom)  

# Apple linker requires dynamic lookup of symbols, so we
# add link flags on this platform
if(APPLE)
  set_target_properties(Falaise_Cimrman
    PROPERTIES
      LINK_FLAGS "-undefined dynamic_lookup"
      INSTALL_RPATH "@loader_path"
    )
endif()

# Install Targets
install(TARGETS Falaise_Cimrman DESTINATION ${CMAKE_INSTALL_PLUGINDIR})


# Test support:
option(FalaiseCimrmanPlugin_ENABLE_TESTING "Build unit testing system for FalaiseCimrmanPlugin" OFF)
if(FalaiseCimrmanPlugin_ENABLE_TESTING)
 enable_testing()
 add_subdirectory(testing)
endif()

# configure_file("CimrmanPipeline.conf.in" "CimrmanPipeline.conf" @ONLY)
