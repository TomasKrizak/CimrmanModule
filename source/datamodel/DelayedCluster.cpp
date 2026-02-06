// Cimrman headers
#include "datamodel/DelayedCluster.h"
#include "datamodel/TrackerHit.h"

// Standard library:
#include <limits>  // infinity
#include <cmath>   // isnan
#include <numeric> // accumulate

// Bayeux headers
#include "datatools/exception.h"

namespace cimrman::datamodel {

  DelayedCluster::DelayedCluster(const std::vector<TrackerHitHdl> & _tracker_hits,
		       double _phi_min,
		       double _phi_max)
    : Cluster()
  {
    DT_THROW_IF(_tracker_hits.size() < 1, std::logic_error, "No tracker hits");
    
    tracker_hits = _tracker_hits;
    phi_min = _phi_min;
    phi_max = _phi_max;
    return;
  }

  void DelayedCluster::set_phi_min(double _phi_min)
  {
    phi_min = _phi_min;
  }

  void DelayedCluster::set_phi_max(double _phi_max)
  {
    phi_max = _phi_max;
  }

  double DelayedCluster::get_phi_min() const
  {
    return phi_min;
  }

  double DelayedCluster::get_phi_max() const
  {
    return phi_max;
  }
  
  void DelayedCluster::set_center_x(double _center_x)
  {
    center_x = _center_x;
  }

  void DelayedCluster::set_center_y(double _center_y)
  {
    center_y = _center_y;
  }

  double DelayedCluster::get_center_x() const
  {
    return center_x;
  }

  double DelayedCluster::get_center_y() const
  {
    return center_y;
  }
  
  void DelayedCluster::set_reference_time_min(double _reference_time_min)
  {
    reference_time_min = _reference_time_min;
  }
  
  double DelayedCluster::get_reference_time_min() const
  {
    return reference_time_min;
  }

  void DelayedCluster::set_reference_time_max(double _reference_time_max)
  {
    reference_time_max = _reference_time_max;
  }
  
  double DelayedCluster::get_reference_time_max() const
  {
    return reference_time_max;
  }

  void DelayedCluster::set_reference_time(double _reference_time)
  {
    reference_time = _reference_time;
  }
  
  double DelayedCluster::get_reference_time() const
  {
    return reference_time;
  }
  
  // calculates the average x,y as the center of the cluster
  void DelayedCluster::find_center()
  {
    if(tracker_hits.empty()) return;
    
    double sum_x = std::accumulate(tracker_hits.begin(), tracker_hits.end(), 0.0,
        [](double acc, const auto& hit){ return acc + hit->get_x();});
        
    double sum_y = std::accumulate(tracker_hits.begin(), tracker_hits.end(), 0.0,
        [](double acc, const auto& hit){ return acc + hit->get_y();});

    double no_hits = static_cast<double>(tracker_hits.size());     
    center_x = sum_x / no_hits;
    center_y = sum_y / no_hits;
  }
  
    
  void DelayedCluster::update_drift_radii()
  { 
    DT_THROW_IF(reference_time == datatools::invalid_real(), std::logic_error, "no reference time to recalculate drift radii!");
    for(auto & hit : tracker_hits)
    {
      //TrackerHit& hit = const_cast<TrackerHit&>(*hit);
      double drift_time = hit->get_delayed_time() - reference_time;
      hit->update_drift_radius( drift_time );
    }
  }

  void DelayedCluster::print(std::ostream & out_) const
  {
    out_ << "Delayed Cluster |"
	 << " size "           << tracker_hits.size()
	 << ", phi min: " << phi_min 
	 << ", phi max: "   << phi_max
	 << ", center x: " << center_x 
	 << ", center y: "   << center_y << std::endl;
  }

} //  end of namespace cimrman::datamodel
