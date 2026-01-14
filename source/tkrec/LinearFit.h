#ifndef FALAISE_CIMRMAN_LINEARFIT_H
#define FALAISE_CIMRMAN_LINEARFIT_H

// Standard headers
#include <iostream>
#include <memory>

// Bayeux headers
#include "datatools/utils.h"

// Cimrman headers
#include "Likelihood.h"

namespace tkrec {

  // Forward declarations
  class Cluster;
  using ConstClusterHdl = std::shared_ptr<const Cluster>;

  // line in the form:
  //	y = ax + b
  //	z = cx + d

  // parametrically:
  //	x(s) = s
  //	y(s) = a*s + b
  //	z(s) = c*s + d

  // or in the form:
  // 	x(t) = cos(phi)*t + r*sin(phi)
  // 	y(t) = sin(phi)*t - r*cos(phi)
  // 	z(t) = tan(theta)*t + h

  // set of transformations:
  //	a = tan(phi)
  //	b = -r/cos(phi)
  //	c = tan(theta)/cos(phi)
  //	d = h - r*tan(phi)*tan(theta)

  //	phi = atan(a) 
  //	r = -b/sqrt(a*a+1.0)
  //	theta = atan(c/sqrt(a*a+1.0))
  //	h = d - a*b*c/(a*a+1.0)
  
  struct LinearFit
  {	
    double phi = datatools::invalid_real();
    double r = datatools::invalid_real();
    double theta = datatools::invalid_real();
    double h = datatools::invalid_real();

    double a = datatools::invalid_real();
    double b = datatools::invalid_real();
    double c = datatools::invalid_real();
    double d = datatools::invalid_real();	

    double chi_squared = datatools::invalid_real();
    
    std::weak_ptr<const Cluster> origin_cluster;

    Likelihood likelihood;
        
    LinearFit(ConstClusterHdl cluster);
    LinearFit(double _phi, double _r, double _theta, double _h);    

    // copy constructor needed for branching into different solutions
    LinearFit(const LinearFit & fit) = default;
    LinearFit(std::shared_ptr<const LinearFit> & fit);

    void print(std::ostream & out_ = std::cout) const;
  };

  typedef std::shared_ptr<LinearFit> LinearFitHdl;
  typedef std::shared_ptr<const LinearFit> ConstLinearFitHdl;

} //  end of namespace tkrec

#endif // FALAISE_CIMRMAN_LINEARFIT_H
