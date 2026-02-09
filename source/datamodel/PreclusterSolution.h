#ifndef FALAISE_CIMRMAN_PRECLUSTERSOLUTION_H
#define FALAISE_CIMRMAN_PRECLUSTERSOLUTION_H

// Standard headers
#include <iostream>
#include <memory>
#include <vector>

namespace cimrman::datamodel {

  // Forward declarations
  class TrackerHit;
  using TrackerHitHdl = std::shared_ptr<TrackerHit>;
  using ConstTrackerHitHdl = std::shared_ptr<const TrackerHit>;
  
  class Track;
  using TrackHdl = std::shared_ptr<Track>;
  using ConstTrackHdl = std::shared_ptr<const Track>;
  
  class Trajectory;
  using TrajectoryHdl = std::shared_ptr<Trajectory>;
  using ConstTrajectoryHdl = std::shared_ptr<const Trajectory>;

 
  class PreclusterSolution
  {
  private:
		
    std::vector<TrackHdl> unprocessed_tracks;
    std::vector<TrajectoryHdl> trajectories;
    std::vector<ConstTrackerHitHdl> unclustered_tracker_hits;
    
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

} //  end of namespace cimrman::datamodel

#endif // FALAISE_CIMRMAN_PRECLUSTERSOLUTION_H
