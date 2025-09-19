#ifndef FALAISE_CIMRMAN_PRECLUSTERSOLUTION_H
#define FALAISE_CIMRMAN_PRECLUSTERSOLUTION_H

// Standard headers
#include <iostream>
#include <memory>

#include "tkrec/Trajectory.h"
#include "tkrec/Track.h"
#include "tkrec/TrackerHit.h"

#include <datatools/logger.h>

namespace tkrec {

  class PreclusterSolution
  {
  private:
		
    std::vector<TrackHdl> unprocessed_tracks;
    std::vector<TrajectoryHdl> trajectories;
    std::vector<ConstTrackerHitHdl> unclustered_tracker_hits;
    
  public:
    
    datatools::logger::priority verbosity = datatools::logger::PRIO_FATAL;
    
  public:
    
    PreclusterSolution() = default;		
    virtual ~PreclusterSolution() = default;    

    std::vector<TrajectoryHdl> & get_trajectories();
    std::vector<ConstTrajectoryHdl> get_trajectories() const;
    
    std::vector<TrackHdl> & get_unprocessed_tracks();
    std::vector<ConstTrackHdl> get_unprocessed_tracks() const;
    
    std::vector<ConstTrackerHitHdl> & get_unclustered_tracker_hits();
    const std::vector<ConstTrackerHitHdl> & get_unclustered_tracker_hits() const;
    
    void add_track(TrackHdl track);
   
    void print(std::ostream & out_ = std::cout) const;
  };

  typedef std::shared_ptr<PreclusterSolution> PreclusterSolutionHdl;
  typedef std::shared_ptr<const PreclusterSolution> ConstPreclusterSolutionHdl;

} //  end of namespace tkrec

#endif // FALAISE_CIMRMAN_PRECLUSTERSOLUTION_H
