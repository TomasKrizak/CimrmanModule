#ifndef FALAISE_CIMRMAN_LIKELIHOOD_H
#define FALAISE_CIMRMAN_LIKELIHOOD_H

// Standard headers
#include <iostream>
#include <memory>
#include <vector>

// Bayeux headers
#include <datatools/utils.h>

namespace cimrman::datamodel {

  // forward definitions
  class Cluster;
  using ConstClusterHdl = std::shared_ptr<const Cluster>;
  
  class TrackerHit;
  using TrackerHitHdl = std::shared_ptr<TrackerHit>;


  struct Likelihood
  {	
    //averages weighted by sigma_R
    double R = 0.0;
    double Rr = 0.0;
    double Rx = 0.0;
    double Ry = 0.0;

    double Rrr = 0.0;
    double Rxx = 0.0;
    double Ryy = 0.0;
    double Rrx = 0.0;
    double Rry = 0.0;
    double Rxy = 0.0;

    double Cov_Rrr = 0.0;
    double Cov_Rxx = 0.0;
    double Cov_Ryy = 0.0;
    double Cov_Rrx = 0.0;
    double Cov_Rry = 0.0;
    double Cov_Rxy = 0.0;

    //averages weighted by sigma_Z
    double Z = 0.0;
    double Zz = 0.0;
    double Zx = 0.0;
    double Zy = 0.0;

    double Zzz = 0.0;
    double Zxx = 0.0;
    double Zyy = 0.0;
    double Zzx = 0.0;
    double Zzy = 0.0;
    double Zxy = 0.0;

    double Cov_Zzz = 0.0;
    double Cov_Zxx = 0.0;
    double Cov_Zyy = 0.0;
    double Cov_Zzx = 0.0;
    double Cov_Zzy = 0.0;
    double Cov_Zxy = 0.0;

    int no_R = 0;
    int no_Z = 0;

    Likelihood() = default;
    Likelihood(ConstClusterHdl & cluster);
    Likelihood(const std::vector<TrackerHitHdl> & hits, std::vector<bool> & signs);
    
    double log_likelihood_value(double phi) const;
    double log_likelihood_derivative(double phi) const;
  };

  typedef std::shared_ptr<Likelihood> LikelihoodHdl;
  typedef std::shared_ptr<const Likelihood> ConstLikelihoodHdl;

} //  end of namespace cimrman::datamodel

#endif // FALAISE_CIMRMAN_LIKELIHOOD_H
