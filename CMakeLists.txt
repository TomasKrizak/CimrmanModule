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
  source/tkrec/Visu.h
  source/tkrec/Algos.h
  source/tkrec/Geometry.h
  source/tkrec/Cluster.h
  source/tkrec/Event.h
  source/tkrec/OMHit.h
  source/tkrec/Point.h
  source/tkrec/Cimrman.h
  source/tkrec/Track.h
  source/tkrec/Trajectory.h
  source/tkrec/TrackerHit.h
  source/tkrec/Solution.h
  source/tkrec/PreclusterSolution.h
  source/tkrec/Precluster.h
  source/tkrec/LinearFit.h
  source/tkrec/Likelihood.h
  )

list(APPEND Cimrman_SOURCES
  source/tkrec/Cimrman.cpp
  source/tkrec/Visu.cpp
  source/tkrec/Algos.cpp
  source/tkrec/Event.cpp
  source/tkrec/Trajectory.cpp
  source/tkrec/Track.cpp
  source/tkrec/Cluster.cpp
  source/tkrec/OMHit.cpp
  source/tkrec/TrackerHit.cpp
  source/tkrec/Point.cpp
  source/tkrec/Solution.cpp
  source/tkrec/PreclusterSolution.cpp
  source/tkrec/Precluster.cpp
  source/tkrec/LinearFit.cpp
  source/tkrec/Likelihood.cpp
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
