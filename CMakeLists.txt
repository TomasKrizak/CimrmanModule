# Check cmake version meets our requirements
cmake_minimum_required(VERSION 3.3)

# Declare project, which will configure compiler for us
project(FalaiseTKReconstructPlugin)

# #Finds ROOT and Falaise packages and defines cmake variables
find_package(Falaise REQUIRED)
find_package(ROOT    REQUIRED)

include_directories(include ${Bayeux_INCLUDE_DIRS})

include_directories(${CMAKE_CURRENT_BINARY_DIR}/source)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/source)

list(APPEND TKReconstruct_HEADERS
  source/tkrec/Visu.h
  source/tkrec/Algos.h
  source/tkrec/Geometry.h
  source/tkrec/Cluster.h
  source/tkrec/Event.h
  source/tkrec/OMHit.h
  source/tkrec/Point.h
  source/tkrec/TKReconstruct.h
  source/tkrec/Track.h
  source/tkrec/Trajectory.h
  source/tkrec/TrackerHit.h
  source/tkrec/Solution.h
  source/tkrec/PreclusterSolution.h
  source/tkrec/Precluster.h
  source/tkrec/LinearFit.h
  source/tkrec/Likelihood.h
  )

list(APPEND TKReconstruct_SOURCES
  source/tkrec/TKReconstruct.cpp
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

add_library(Falaise_TKReconstruct SHARED
  ${TKReconstruct_SOURCES}
  ${TKReconstruct_HEADERS} )


target_link_libraries(Falaise_TKReconstruct Falaise::Falaise ${ROOT_LIBRARIES} ROOT::Geom)
# PUBLIC Falaise::Falaise FalaiseModule ${ROOT_LIBRARIES} ROOT::Geom)  

# Apple linker requires dynamic lookup of symbols, so we
# add link flags on this platform
if(APPLE)
  set_target_properties(Falaise_TKReconstruct
    PROPERTIES
      LINK_FLAGS "-undefined dynamic_lookup"
      INSTALL_RPATH "@loader_path"
    )
endif()

# Install Targets
install(TARGETS Falaise_TKReconstruct DESTINATION ${CMAKE_INSTALL_PLUGINDIR})


# Test support:
option(FalaiseTKReconstructPlugin_ENABLE_TESTING "Build unit testing system for FalaiseTKReconstructPlugin" OFF)
if(FalaiseTKReconstructPlugin_ENABLE_TESTING)
 enable_testing()
 add_subdirectory(testing)
endif()

# configure_file("TKReconstructPipeline.conf.in" "TKReconstructPipeline.conf" @ONLY)
