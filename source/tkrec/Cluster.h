#ifndef FALAISE_CIMRMAN_CLUSTER_H
#define FALAISE_CIMRMAN_CLUSTER_H

// Standard headers
#include <iostream>
#include <memory>

// Cimrman headers
#include "tkrec/TrackerHit.h"
#include "tkrec/LinearFit.h"

#include <datatools/logger.h>

namespace tkrec {

  class Cluster
  {
  protected:
		
    // 0 == no ambiguity
    // 1 == mirror image along line x = x0 
    // 2 == mirror image along line y = y0 
    // 3 == mirror image along line y = x + (y0-x0) 
    // 4 == mirror image along line y = -x + (y0-x0) 
    int ambiguity_type = 0;
    
    std::vector<TrackerHitHdl> tracker_hits;
    std::vector<LinearFitHdl> linear_fits;
    
    double phi_estimate = datatools::invalid_real();
    double r_estimate = datatools::invalid_real();

  public:
    
    datatools::logger::priority verbosity = datatools::logger::PRIO_FATAL;
    
  public:
    
    Cluster() = default;		
    Cluster(const std::vector<TrackerHitHdl> & _tracker_hits,
		        double _phi_estimate,
		        double _r_estimate);
    virtual ~Cluster() = default;
		
    std::vector<TrackerHitHdl> & get_tracker_hits();
    const std::vector<ConstTrackerHitHdl> get_const_tracker_hits() const;
    
    std::vector<LinearFitHdl> & get_linear_fits();
    std::vector<ConstLinearFitHdl> get_linear_fits() const;
     
    void set_phi_estimate(double _phi_estimate);
    double get_phi_estimate() const;
		
    void set_r_estimate(double _r_estimate);
    double get_r_estimate() const;

    void set_ambiguity_type(int _ambiguity_type);
    int get_ambiguity_type() const;

    void print(std::ostream & out_ = std::cout) const;
		
  };

  typedef std::shared_ptr<Cluster> ClusterHdl;
  typedef std::shared_ptr<const Cluster> ConstClusterHdl;

} //  end of namespace tkrec

#endif // FALAISE_CIMRMAN_CLUSTER_H
