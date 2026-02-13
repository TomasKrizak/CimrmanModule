#ifndef FALAISE_CIMRMAN_TRAJECTORYBUILDER_H
#define FALAISE_CIMRMAN_TRAJECTORYBUILDER_H

// Standard headers
#include <vector>
#include <memory>
#include <tuple>

// Cimrman headers
#include "reconstruction/KinkFinder.h"
#include "datamodel/Trajectory.h"

namespace cimrman {
  
  // Forward declarations
  namespace datamodel {
    class PreclusterSolution;
    using PreclusterSolutionHdl = std::shared_ptr<PreclusterSolution>;
  
    class Trajectory;
    using TrajectoryHdl = std::shared_ptr<Trajectory>;
    
    class Point;
    using PointHdl = std::shared_ptr<Point>;
  }
  
  struct Geometry;
  struct PolylinesConfig;
  
  
  class TrajectoryBuilder
  {
  private:
  
    datamodel::PreclusterSolutionHdl precluster_solution = nullptr;
    std::vector<KinkFinder::ConnectionStrategy> connection_strategies;
    const PolylinesConfig & config;
    const Geometry & geom;
                                                
    static std::pair<datamodel::TrackHdl, datamodel::PointHdl> get_trajectory_end(
                                                          datamodel::TrajectoryHdl trajectory, 
                                                          const datamodel::Trajectory::EndPoint endpoint);
                                                          
    KinkFinder::Status try_connecting(datamodel::TrajectoryHdl base_trajectory,
                                      datamodel::TrajectoryHdl other_trajectory);                         
  
    // merges "other_trajectory" into "base_trajectory" based on the chosen ends
    // and provided kink point, leaving other_trajectory empty! 
    static void merge_trajectories_into_base(datamodel::TrajectoryHdl base_trajectory,
                                             datamodel::Trajectory::EndPoint base_endpoint, 
                                             datamodel::TrajectoryHdl other_trajectory, 
                                             datamodel::Trajectory::EndPoint other_endpoint,
                                             datamodel::PointHdl kink_point);

  public:
  
    TrajectoryBuilder() = default;
    
    TrajectoryBuilder(datamodel::PreclusterSolutionHdl _precluster_solution,
                      const PolylinesConfig & _config,
                      const Geometry & _geom);
  
    void process();
    
    void add_connection_strategy(KinkFinder::ConnectionStrategy _connection_strategy);
  
  };

} //  end of namespace cimrman

#endif // FALAISE_CIMRMAN_TRAJECTORYBUILDER_H
