// Cimrman headers
#include "datamodel/PreclusterSolution.h"
#include "datamodel/Trajectory.h"
#include "datamodel/Track.h"
#include "datamodel/TrackerHit.h"

// Bayeux:
#include "datatools/exception.h"

namespace cimrman::datamodel {  
  
  void PreclusterSolution::add_track(TrackHdl track)
  {
    unprocessed_tracks.push_back(track);	
    return;
  }
  
  std::vector<TrajectoryHdl> & PreclusterSolution::get_trajectories()
  {
    return trajectories;
  }  
  
  std::vector<ConstTrajectoryHdl> PreclusterSolution::get_trajectories() const
  {
    std::vector<ConstTrajectoryHdl> traj;
    for(const auto & tr : trajectories)
    {
      traj.push_back(tr);
    }
      return traj;
  }
  
  std::vector<TrackHdl> & PreclusterSolution::get_unprocessed_tracks()
  {
    return unprocessed_tracks;
  }  
  
  std::vector<ConstTrackHdl> PreclusterSolution::get_unprocessed_tracks() const
  {
    std::vector<ConstTrackHdl> tracks;
    for(const auto & track : unprocessed_tracks)
    {
      tracks.push_back(track);
    }
    return tracks;
  }

  std::vector<ConstTrackerHitHdl> & PreclusterSolution::get_unclustered_tracker_hits()
  {
    return unclustered_tracker_hits;
  }

  const std::vector<ConstTrackerHitHdl> & PreclusterSolution::get_unclustered_tracker_hits() const
  {
    return unclustered_tracker_hits;
  }

  void PreclusterSolution::print(std::ostream & out_) const
  {
    out_ <<"Precluster solution: " << std::endl;
    out_ << "	" << trajectories.size() << " precluster trajectories: " << std::endl;
    for(const auto & traj : trajectories)
    {
      out_ << "	";
      traj->print();
    }
    return;
  }

} //  end of namespace cimrman::datamodel
