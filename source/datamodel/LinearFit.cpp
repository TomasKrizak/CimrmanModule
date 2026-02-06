// Cimrman headers
#include "datamodel/LinearFit.h"
#include "datamodel/Cluster.h"

// Standard headers
#include <cmath>

namespace cimrman::datamodel {
  
  LinearFit::LinearFit(double _phi, double _r, double _theta, double _h)
  {
    phi   = _phi;
    r     = _r;
    theta = _theta;
    h     = _h;
    
    a = std::tan(phi);
    b = -r / std::cos(phi);
    c = std::tan(theta) / std::cos(phi);
    d = h - r * std::tan(phi) * std::tan(theta);
    return;
  }

  LinearFit::LinearFit(ConstClusterHdl cluster)
  {
    phi   = cluster->get_phi_estimate();
    r     = cluster->get_r_estimate();
    theta = 0.0;
    h     = 0.0; 
    
    a = std::tan(phi);
    b = -r / std::cos(phi);
    c = 0.0;
    d = 0.0;
    
    likelihood = Likelihood(cluster);
    
    origin_cluster = cluster;
    return;
  }

  LinearFit::LinearFit(ConstLinearFitHdl & fit)
	: LinearFit(*fit)
  {
    return;
  }

  void LinearFit::print(std::ostream & out_) const
  {
    out_ << "Track: "
	 << ", a = " << a 
	 << ", b = " << b 
	 << ", c = " << c 
	 << ", d = " << d << std::endl; 
    out_ << "		 phi = " << phi 
	 << ", r = " << r 
	 << ", theta = " << theta 
	 << ", h = " << h << std::endl; 
    out_ << "	chi squared: " << chi_squared << std::endl;
  }

} //  end of namespace cimrman::datamodel
