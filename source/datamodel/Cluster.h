#ifndef FALAISE_CIMRMAN_CLUSTER_H
#define FALAISE_CIMRMAN_CLUSTER_H

// Standard headers
#include <iostream>
#include <memory>
#include <vector>

// Bayeux headers
#include <datatools/utils.h>

namespace cimrman::datamodel {

  // forward declarations
  class TrackerHit;
  using TrackerHitHdl = std::shared_ptr<TrackerHit>;
  using ConstTrackerHitHdl = std::shared_ptr<const TrackerHit>;

  class LinearFit;
  using LinearFitHdl = std::shared_ptr<LinearFit>;
  using ConstLinearFitHdl = std::shared_ptr<const LinearFit>;
  


  // 0 == no ambiguity
  // 1 == mirror image along line x = x0 
  // 2 == mirror image along line y = y0 
  // 3 == mirror image along line y = x + (y0-x0) 
  // 4 == mirror image along line y = -x + (y0-x0) 
  enum class Ambiguity
  {
    None = 0,
    Vertical = 1,
    Horizontal = 2,
    AntiDiagonal = 3,
    MainDiagonal = 4
  };


  class Cluster
  {
  protected:
		
    Ambiguity ambiguity_type = Ambiguity::None;
    
    std::vector<TrackerHitHdl> tracker_hits;
    std::vector<LinearFitHdl> linear_fits;
    
    double phi_estimate = datatools::invalid_real();
    double r_estimate = datatools::invalid_real();
    
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

    void set_ambiguity_type(Ambiguity _ambiguity_type);
    Ambiguity get_ambiguity_type() const;

    void print(std::ostream & out_ = std::cout) const;
		
  };

  typedef std::shared_ptr<Cluster> ClusterHdl;
  typedef std::shared_ptr<const Cluster> ConstClusterHdl;

} //  end of namespace cimrman::datamodel

#endif // FALAISE_CIMRMAN_CLUSTER_H
