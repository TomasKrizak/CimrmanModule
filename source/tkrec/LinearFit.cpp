// TK headers
#include "tkrec/LinearFit.h"
#include "tkrec/Cluster.h" 

// ClassImp(tkrec::LinearFit);

namespace tkrec {
  
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
    return;
  }

  LinearFit::LinearFit(ConstLinearFitHdl & fit)
	: LinearFit(*fit)
  {
    return;
  }

/*
  void TKtrack::set_chi_squared(double _chi_squared)
  {
    chi_squared = _chi_squared;
    return;
  }

  void TKtrack::set_chi_squared_R(double _chi_squared_R)
  {
    chi_squared_R = _chi_squared_R;
    return;
  }

  void TKtrack::set_chi_squared_Z(double _chi_squared_Z)
  {
    chi_squared_Z = _chi_squared_Z;
    return;
  }

  void TKtrack::set_quality(double _quality)
  {
    quality = _quality;
    return;
  }

  void TKtrack::set_quality_R(double _quality_R)
  {
    quality_R = _quality_R;
    return;
  }

  void TKtrack::set_quality_Z(double _quality_Z)
  {
    quality_Z = _quality_Z;
    return;
  }

  void TKtrack::set_likelihood(double _likelihood)
  {
    likelihood = _likelihood;
    return;
  }

  void TKtrack::set_likelihood_R(double _likelihood_R)
  {
    likelihood_R = _likelihood_R;
    return;
  }

  void TKtrack::set_likelihood_Z(double _likelihood_Z)
  {
    likelihood_Z = _likelihood_Z;
    return;
  }
*/

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
    out_ << "	chi squared R: " << chi_squared_R << std::endl;
    out_ << "	chi squared Z: " << chi_squared_Z << std::endl;
  }

} //  end of namespace tkrec

