#ifndef FALAISE_CIMRMAN_DELAYEDCLUSTER_H
#define FALAISE_CIMRMAN_DELAYEDCLUSTER_H

// Standard library:
#include <cmath> // M_PI

// Cimrman headers
#include "tkrec/Cluster.h"

namespace tkrec {

  class DelayedCluster : public Cluster
  {
  private:
		
    double phi_min = datatools::invalid_real();
    double phi_max = datatools::invalid_real();
    double center_x = datatools::invalid_real();
    double center_y = datatools::invalid_real();

    double reference_time_min = datatools::invalid_real();
    double reference_time_max = datatools::invalid_real();
    double reference_time = datatools::invalid_real();

  public:
    
    datatools::logger::priority verbosity = datatools::logger::PRIO_FATAL;
    
  public:
    
    DelayedCluster() = default;		
    DelayedCluster(const std::vector<TrackerHitHdl> & _tracker_hits,
		        double _phi_min,
		        double _phi_max,
		        double _r_min,
		        double _r_max);
		        
    DelayedCluster(const std::vector<TrackerHitHdl> & _tracker_hits,
		        double _phi_min = 0.0,
		        double _phi_max = M_PI);

    // calculates the average x,y as the center of the cluster
    void find_center();
    
    // recalculates drift radii of cluster hits based on set reference time
    void update_drift_radii();

    void set_phi_min(double _phi_min);
    double get_phi_min() const;
		
    void set_phi_max(double _phi_max);
    double get_phi_max() const;
    
    void set_center_x(double _center_x);
    double get_center_x() const;
    
    void set_center_y(double _center_y);
    double get_center_y() const;
    
    void set_reference_time_min(double _reference_time_min);
    double get_reference_time_min() const;
    
    void set_reference_time_max(double _reference_time_max);
    double get_reference_time_max() const;
    
    void set_reference_time(double _reference_time);
    double get_reference_time() const;

    void print(std::ostream & out_ = std::cout) const;
  };

  typedef std::shared_ptr<DelayedCluster> DelayedClusterHdl;
  typedef std::shared_ptr<const DelayedCluster> ConstDelayedClusterHdl;

} //  end of namespace tkrec

#endif // FALAISE_CIMRMAN_DELAYEDCLUSTER_H
