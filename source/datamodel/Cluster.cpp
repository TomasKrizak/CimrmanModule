// Cimrman headers
#include "datamodel/Cluster.h"
#include "datamodel/LinearFit.h"
#include "datamodel/TrackerHit.h"

namespace cimrman::datamodel {

  Cluster::Cluster(const std::vector<TrackerHitHdl> & _tracker_hits,
		       double _phi_estimate,
		       double _r_estimate)
    : Cluster()
    
  {
    DT_THROW_IF(_tracker_hits.size() < 1, std::logic_error, "No tracker hits");
    tracker_hits = _tracker_hits;
    phi_estimate = _phi_estimate; 
    r_estimate = _r_estimate;
    return;
  }

  std::vector<TrackerHitHdl> & Cluster::get_tracker_hits()
  {
    return tracker_hits;
  }
  
  const std::vector<ConstTrackerHitHdl> Cluster::get_const_tracker_hits() const 
  {
    std::vector<ConstTrackerHitHdl> hits;
    for (const auto & hit : tracker_hits)
    {
      hits.push_back(hit);
    }
    return hits;
  }

  std::vector<LinearFitHdl> & Cluster::get_linear_fits()
  {
    return linear_fits;
  }

  std::vector<ConstLinearFitHdl> Cluster::get_linear_fits() const
  {
    std::vector<ConstLinearFitHdl> fits;
    for (const auto & fit : linear_fits)
    {
      fits.push_back(fit);
    }
    return fits;
  }

  void Cluster::set_phi_estimate(double _phi_estimate)
  {
    phi_estimate = _phi_estimate;
  }

  void Cluster::set_r_estimate(double _r_estimate)
  {
    r_estimate = _r_estimate;
  }

  void Cluster::set_ambiguity_type(Ambiguity _ambiguity_type)
  {
    ambiguity_type = _ambiguity_type;
  }

  Ambiguity Cluster::get_ambiguity_type() const
  {
    return ambiguity_type;
  }

  double Cluster::get_phi_estimate() const
  {
    return phi_estimate;
  }

  double Cluster::get_r_estimate() const
  {
    return r_estimate;
  }

  void Cluster::print(std::ostream & out_) const
  {
    out_ << "Cluster |"
	 << " size: "           << tracker_hits.size()
	 << ", phi estimate: " << phi_estimate 
	 << ", r estimate: "   << r_estimate 
	 << ", ambiguity type: " << int(ambiguity_type) << std::endl;
  }

} //  end of namespace cimrman::datamodel
