// TK headers
#include "tkrec/Likelihood.h"
#include <tkrec/Cluster.h>

// ClassImp(tkrec::Likelihood);

namespace tkrec 
{

  Likelihood::Likelihood(ConstClusterHdl & cluster)
  {
  
    double x, y, z, r;
    double const_R, const_Z;
    no_R = cluster->get_tracker_hits().size();
    for(const auto & hit : cluster->get_tracker_hits())
    {
      x = hit->get_x();
      y = hit->get_y();
      r = hit->get_R();
      z = hit->get_Z();
      
      // sign of drift radius differs in calculations based on which side of the anode wire the track is
      if(cluster->get_r_estimate() - x * std::sin(cluster->get_phi_estimate()) + y * std::cos(cluster->get_phi_estimate()) < 0.0) 
      {
	      r *= -1.0;
      }

      const_R = std::pow(hit->get_sigma_R(), -2.0);

      R  += const_R;
      Rr += r * const_R;
      Rx += x * const_R;
      Ry += y * const_R;

      Rrr += r * r * const_R;
      Rrx += r * x * const_R;
      Rry += r * y * const_R;
      Rxx += x * x * const_R;
      Ryy += y * y * const_R;
      Rxy += x * y * const_R;

      if(hit->has_valid_Z())
      {
	    no_Z++;
	    const_Z = std::pow(hit->get_sigma_Z(), -2.0);

	    Z  += const_Z;
	    Zz += z * const_Z;
	    Zx += x * const_Z;
	    Zy += y * const_Z;
	    
	    Zzz += z * z * const_Z;
	    Zzx += z * x * const_Z;
	    Zzy += z * y * const_Z;
	    Zxx += x * x * const_Z;
	    Zyy += y * y * const_Z;
	    Zxy += x * y * const_Z;
      }
    }

    Cov_Rrr = Rrr / R - std::pow(Rr / R, 2.0);
    Cov_Rrx = Rrx / R - (Rr / R) * (Rx / R);
    Cov_Rry = Rry / R - (Rr / R) * (Ry / R);
    Cov_Rxx = Rxx / R - std::pow(Rx / R, 2.0);
    Cov_Ryy = Ryy / R - std::pow(Ry / R, 2.0);
    Cov_Rxy = Rxy / R - (Rx / R) * (Ry / R);

    if(no_Z > 0)
    {
      Cov_Zzz = Zzz / Z - std::pow(Zz / Z, 2.0);
      Cov_Zzx = Zzx / Z - (Zz / Z) * (Zx / Z);
      Cov_Zzy = Zzy / Z - (Zz / Z) * (Zy / Z);
      Cov_Zxx = Zxx / Z - std::pow(Zx / Z, 2.0);
      Cov_Zyy = Zyy / Z - std::pow(Zy / Z, 2.0);
      Cov_Zxy = Zxy / Z - (Zx / Z) * (Zy / Z);
    }	  
  	return;
  }

  double Likelihood::log_likelihood_value(double phi) const
  {
    double sin_phi = std::sin(phi);
    double cos_phi = std::cos(phi);

    double lik_R = Cov_Rrr
    + Cov_Rxx * std::pow(sin_phi, 2.0) 
    + Cov_Ryy * std::pow(cos_phi, 2.0) 
    - 2.0 * Cov_Rry * cos_phi 
    + 2.0 * Cov_Rrx * sin_phi
    - 2.0 * Cov_Rxy * sin_phi * cos_phi;

    double lik_Z = 0;
    if(no_Z > 0)
    {		
      lik_Z = std::pow(Cov_Zzy * sin_phi + Cov_Zzx * cos_phi, 2.0);
      lik_Z /= Cov_Zxx * std::pow(cos_phi, 2.0) 
      + Cov_Zyy * std::pow(sin_phi, 2.0) 
      + 2.0 * sin_phi * cos_phi * Cov_Zxy;
      lik_Z = Cov_Zzz - lik_Z;
    }

    return -0.5 * (R * lik_R + Z * lik_Z);
  }

  double Likelihood::log_likelihood_derivative(double phi) const
  {
    double sin_phi = std::sin(phi);
    double cos_phi = std::cos(phi);

    double der_R = (Cov_Rxx - Cov_Ryy) * sin_phi * cos_phi
    - Cov_Rxy * (2.0 * cos_phi * cos_phi - 1.0)
    + Cov_Rry * sin_phi
    + Cov_Rrx * cos_phi;

    double der_Z = 0;
    if(no_Z > 0)
    {
      der_Z = -( Cov_Zzy * sin_phi + Cov_Zzx * cos_phi );
      der_Z *= sin_phi * (Cov_Zzy * Cov_Zxy - Cov_Zzx * Cov_Zyy)
      + cos_phi * (Cov_Zzy * Cov_Zxx - Cov_Zzx * Cov_Zxy);
      double temp = Cov_Zyy * sin_phi * sin_phi 
      + 2.0 * Cov_Zxy * sin_phi * cos_phi
      + Cov_Zxx * cos_phi * cos_phi;
      der_Z /= std::pow( temp, 2.0 );
    }

    return -(R * der_R + Z * der_Z);
  }

} //  end of namespace tkrec
