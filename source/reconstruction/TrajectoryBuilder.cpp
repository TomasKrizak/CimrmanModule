// Cimrman headers
#include "reconstruction/TrajectoryBuilder.h"
#include "datamodel/PreclusterSolution.h"
#include "datamodel/Point.h"

// Standard headers
//#include <iostream>

using namespace cimrman::datamodel;

namespace cimrman {

  TrajectoryBuilder::TrajectoryBuilder(datamodel::PreclusterSolutionHdl _precluster_solution, 
                                      const PolylinesConfig & _config,
                                      const Geometry & _geom)
    : precluster_solution(_precluster_solution),
      config(_config),
      geom(_geom)
  {
    return;
  }

  void TrajectoryBuilder::add_connection_strategy(KinkFinder::ConnectionStrategy _connection_strategy)
  {
    connection_strategies.push_back(_connection_strategy);
  }
  
  // merges "other_trajectory" into "base_trajectory" based on the chosen ends
  // and provided kink point, leaving other_trajectory empty! 
  void TrajectoryBuilder::merge_trajectories_into_base(TrajectoryHdl base_trajectory, Trajectory::EndPoint base_endpoint, 
                                                       TrajectoryHdl other_trajectory, Trajectory::EndPoint other_endpoint,
                                                       PointHdl kink_point)
  {
    // Trajectory contains vector of points of the polyline and the vector of segments 
    auto & base_trajectory_points  =  base_trajectory->get_trajectory_points();
    auto & other_trajectory_points = other_trajectory->get_trajectory_points();

    auto & base_segments  =  base_trajectory->get_segments();
    auto & other_segments = other_trajectory->get_segments();

    DT_THROW_IF(base_segments.empty(), 
      std::logic_error, "Base trajectory is empty");

    DT_THROW_IF(other_segments.empty(),
      std::logic_error, "Other trajectory is empty");
    
    // connection is always done by inserting other_XXX vector to the end of base_XXX vector
    // ( that means BACK-FRONT connection ), so the trajectories is first reoriented accordingly
    
    // Orienting base_trajectory
    if( base_endpoint == Trajectory::EndPoint::FRONT )
    {
      // reverses the sotrage order of trajectory_points and segments
      base_trajectory->reverse_representation();
    }
    
    // Orienting other_trajectory
    if( other_endpoint == Trajectory::EndPoint::BACK )
    {
      // reverses the sotrage order of trajectory_points and segments
      other_trajectory->reverse_representation();
    }
    
    // last point of base_trajectory is replaced with the kink point
    base_trajectory_points.back() = kink_point;
    
    // connecting trajectory points (first point of other_trajectory is omitted)
    base_trajectory_points.insert(base_trajectory_points.end(),
                                  other_trajectory_points.begin() + 1, 
                                  other_trajectory_points.end());      
    
    // connecting segments
    base_segments.insert(base_segments.end(), 
                         other_segments.begin(), 
                         other_segments.end());
    
    // in case the base_trajectory was not a polyline before
    base_trajectory->mark_as_kinked();
    
    other_trajectory_points.clear();
    other_segments.clear();
  }
  
  
  std::pair<TrackHdl, PointHdl> TrajectoryBuilder::get_trajectory_end(TrajectoryHdl trajectory, 
                                                                      const Trajectory::EndPoint endpoint)
  {
    auto & trajectory_points = trajectory->get_trajectory_points();
    auto & segments = trajectory->get_segments();
    if(endpoint == Trajectory::EndPoint::FRONT)
    {
      return { segments.front() , trajectory_points.front()};
    }
    else
    {
      return { segments.back() , trajectory_points.back()};
    }
  }
  
  struct EndPointPair {
    datamodel::Trajectory::EndPoint endpoint1;
    datamodel::Trajectory::EndPoint endpoint2;
  };
  
  static constexpr std::array<EndPointPair ,4> endpoint_order = {{
    {datamodel::Trajectory::EndPoint::BACK,  datamodel::Trajectory::EndPoint::FRONT},  
    {datamodel::Trajectory::EndPoint::FRONT, datamodel::Trajectory::EndPoint::BACK},
    {datamodel::Trajectory::EndPoint::BACK,  datamodel::Trajectory::EndPoint::BACK},
    {datamodel::Trajectory::EndPoint::FRONT, datamodel::Trajectory::EndPoint::FRONT}
  }};
  
  KinkFinder::Status TrajectoryBuilder::try_connecting(TrajectoryHdl base_trajectory, TrajectoryHdl other_trajectory)
  {
    // trying all combinations of the trajectories endpoints
    for(auto [base_endpoint, other_endpoint] : endpoint_order)
    {
      auto [segment1, point1] = get_trajectory_end(base_trajectory, base_endpoint);
      auto [segment2, point2] = get_trajectory_end(other_trajectory, other_endpoint);
            
      KinkFinder finder(segment1, point1,
                        segment2, point2,
                        config, geom);
                  
      for(const auto & connection_strategy : connection_strategies)
      {
        finder.process(connection_strategy);
        
        if(finder.get_status() == KinkFinder::Status::KINK_FOUND)
        {
          Point kink_point = finder.get_kink_point();
          PointHdl kink_point_hdl = std::make_shared<Point>(kink_point);
          merge_trajectories_into_base(base_trajectory, base_endpoint, 
                                       other_trajectory, other_endpoint,
                                       kink_point_hdl);
          
          return KinkFinder::Status::KINK_FOUND;
        }
      }
    }
    return KinkFinder::Status::KINK_NOT_FOUND;
  }
      

  void TrajectoryBuilder::process()
  {
    // algorithm going through pairs of trajectories, trying to connecting them
    // greedy algorithm - connecting as as possible to one trajectory, then moving to another one
    auto i = 0u;
    auto & trajectories = precluster_solution->get_trajectories(); 
    while(i < trajectories.size())
    {       
      bool found_connection = false;
      auto & base_trajectory = trajectories[i];
      
      for(auto j = i + 1; j < trajectories.size(); ++j)
      {
        auto & other_trajectory = trajectories[j];
        
        // try_connecting investigates all possible combinations of ends and all allowed connection strategies
        KinkFinder::Status status = try_connecting( base_trajectory, other_trajectory);
        if( status == KinkFinder::Status::KINK_FOUND )
        {
          trajectories.erase(trajectories.begin() + j);
          found_connection = true;
          break;
        }
      }
      if(not found_connection)
      {
        i++;
      }
    }
  }      
        
      

} //  end of namespace cimrman
