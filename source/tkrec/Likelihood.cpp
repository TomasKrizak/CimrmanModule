// Cimrman headers
#include "tkrec/Likelihood.h"
#include <tkrec/Cluster.h>

// ClassImp(tkrec::Likelihood);

namespace tkrec 
{

  Likelihood::Likelihood(ConstClusterHdl & cluster)
  {
    double x, y, z, r;
    double const_R, const_Z;
    const auto hits = cluster->get_const_tracker_hits();
    no_R = hits.size();
    for(const auto & hit : hits)
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
      const_R = 1.0 / (hit->get_sigma_R() * hit->get_sigma_R());

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
	      const_Z = 1.0 / (hit->get_sigma_Z() * hit->get_sigma_Z());

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

    Cov_Rrr = (Rrr - Rr * Rr / R) / R;
    Cov_Rrx = (Rrx - Rr * Rx / R) / R;
    Cov_Rry = (Rry - Rr * Ry / R) / R;
    Cov_Rxx = (Rxx - Rx * Rx / R) / R;
    Cov_Ryy = (Ryy - Ry * Ry / R) / R;
    Cov_Rxy = (Rxy - Rx * Ry / R) / R;

    if(no_Z > 0)
    {
      Cov_Zzz = (Zzz - Zz * Zz / Z) / Z;
      Cov_Zzx = (Zzx - Zz * Zx / Z) / Z;
      Cov_Zzy = (Zzy - Zz * Zy / Z) / Z;
      Cov_Zxx = (Zxx - Zx * Zx / Z) / Z;
      Cov_Zyy = (Zyy - Zy * Zy / Z) / Z;
      Cov_Zxy = (Zxy - Zx * Zy / Z) / Z;
    }	  
    
  	return;
  }

  double Likelihood::log_likelihood_value(double phi) const
  {
    double sin_phi = std::sin(phi);
    double cos_phi = std::cos(phi);

    double lik_R = Cov_Rrr
    + Cov_Rxx * sin_phi * sin_phi 
    + Cov_Ryy * cos_phi * cos_phi 
    - 2.0 * Cov_Rry * cos_phi 
    + 2.0 * Cov_Rrx * sin_phi
    - 2.0 * Cov_Rxy * sin_phi * cos_phi;

    double lik_Z = 0;
    if(no_Z > 0)
    {		
      double temp = Cov_Zzy * sin_phi + Cov_Zzx * cos_phi;
      lik_Z = temp * temp;
      lik_Z /= Cov_Zxx * cos_phi * cos_phi 
            + Cov_Zyy * sin_phi * sin_phi 
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
      der_Z /= (temp * temp);
    }

    return -(R * der_R + Z * der_Z);
  }

} //  end of namespace tkrec
