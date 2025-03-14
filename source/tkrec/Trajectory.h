#ifndef FALAISE_TKRECONSTRUCT_TRAJECTORY_H
#define FALAISE_TKRECONSTRUCT_TRAJECTORY_H

// Standard headers
#include <iostream>
#include <vector>

#include "tkrec/TrackerHit.h"
#include "tkrec/Track.h"
#include "tkrec/Point.h"

namespace tkrec {

  class Trajectory
  {
  private:
		
    std::vector<PointHdl> trajectory_points; // startpoint, kink points and endpoint of trajectory
    std::vector<TrackHdl> segments; // linear parts of the polyline trajectory
    bool kinked_trajectory = false;
    
		
  public:
		
    Trajectory() = default;
    Trajectory(TrackHdl & segment);
    Trajectory(std::vector<TrackHdl> & _segments);
    virtual ~Trajectory() = default;
		
    std::vector<PointHdl> & get_trajectory_points();
    std::vector<ConstPointHdl> get_trajectory_points() const;
    
    std::vector<TrackHdl> & get_segments();
    std::vector<ConstTrackHdl> get_segments() const;
    
    const std::vector<Association> get_associations() const;

    void mark_as_kinked();
    void mark_as_straight();
    
    bool has_kink() const;

    void update_segments();

    void print(std::ostream & out_ = std::clog) const;
		    
  };
  
  typedef std::shared_ptr<Trajectory> TrajectoryHdl;
  typedef std::shared_ptr<const Trajectory> ConstTrajectoryHdl;

} //  end of namespace tkrec

#endif // FALAISE_TKRECONSTRUCT_TRAJECTORY_H
