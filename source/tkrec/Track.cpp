// TK headers
#include "tkrec/Track.h"

// ClassImp(tkrec::Track);

namespace tkrec {
  
  using namespace std;

  Track::Track(const ConstLinearFitHdl & linear_fit, const std::vector<ConstTrackerHitHdl> & tracker_hits)
    : fit(std::make_shared<LinearFit>(*linear_fit))
  {
    double sin_phi = std::sin(fit->phi);
    double cos_phi = std::cos(fit->phi);
    double tan_theta = std::tan(fit->theta);
    for(auto & hit : tracker_hits)
    {
      double t0 = hit->get_x() * cos_phi + hit->get_y() * sin_phi;

      double x0 = t0 * cos_phi + fit->r * sin_phi;
      double y0 = t0 * sin_phi - fit->r * cos_phi;
      double z0 = t0 * tan_theta + fit->h;

       // creating Point obejcts to represent reconstructed avalanche origin points
      Association association(hit);
      association.parameter_t = t0;
      association.point = std::make_shared<Point>(x0, y0, z0);
      associations.push_back( association );
    }
    return;
  }
  
  void Track::add_associated_tracker_hit(ConstTrackerHitHdl & tracker_hit)
  {
    Association association(tracker_hit);

    double t0 = tracker_hit->get_x() * std::cos(fit->phi) + tracker_hit->get_y() * std::sin(fit->phi);
    association.parameter_t = t0;
    
    double x0 = t0 * std::cos(fit->phi) + fit->r * std::sin(fit->phi);
    double y0 = t0 * std::sin(fit->phi) - fit->r * std::cos(fit->phi);
    double z0 = t0 * std::tan(fit->theta) + fit->h;

     // creating Point obejcts to represent reconstructed avalanche origin points
    association.point = std::make_shared<Point>(x0, y0, z0);

    associations.push_back( association );
    return;
  }
  
  Association::Association(const ConstTrackerHitHdl & hit)
    : tracker_hit(hit)
  {
    return;
  }
  
  // sorting the track associations based on position along the track
  void Track::sort_associations()
  {
     std::sort(associations.begin(), associations.end(), 
              [&](Association assoc1, Association assoc2)
              {
                return assoc1.parameter_t < assoc2.parameter_t;
              });
     return;
  }
  
  // TODO can I add static storage for sin(phi), cos(phi)?
  double Track::calculate_t(double x, double y) const      
  {
    return x * std::cos(fit->phi) + y * std::sin(fit->phi); 
  }
  
  double Track::calculate_t(ConstPointHdl & point) const
  {
    return point->x * std::cos(fit->phi) + point->y * std::sin(fit->phi);
  }
  
  double Track::calculate_t(const Point & point) const
  {
    return point.x * std::cos(fit->phi) + point.y * std::sin(fit->phi);
  }


  double Track::horizontal_distance_to_line(const Point & point) const
  {
    return std::abs(fit->r - point.x * std::sin(fit->phi) + point.y * std::cos(fit->phi)); 
  }

  void Track::update_associations()
  {
    double sin_phi = std::sin(fit->phi);
    double cos_phi = std::cos(fit->phi);
    double tan_theta = std::tan(fit->theta);
    for(auto & association : associations)
    {
      ConstTrackerHitHdl & hit = association.tracker_hit;
      double t0 = hit->get_x() * cos_phi + hit->get_y() * sin_phi;
      
      association.parameter_t = t0;
      association.point->x = t0 * cos_phi + fit->r * sin_phi;
      association.point->y = t0 * sin_phi - fit->r * cos_phi;
      association.point->z = t0 * tan_theta + fit->h;
    }
  }


  void Track::set_a(double _a)
  {
    fit->a = _a;
  }

  void Track::set_b(double _b)
  {
    fit->b = _b;
  }

  void Track::set_c(double _c)
  {
    fit->c = _c;
  }

  void Track::set_d(double _d)
  {
    fit->d = _d;
  }

  void Track::set_phi(double _phi)
  {
    fit->phi = _phi;
  }

  void Track::set_r(double _r)
  {
    fit->r = _r;
  }

  void Track::set_theta(double _theta)
  {
    fit->theta = _theta;
  }

  void Track::set_h(double _h)
  {
    fit->h = _h;
  }

  double Track::get_a() const
  {
    return fit->a;
  } 

  double Track::get_b() const
  {
    return fit->b;
  } 

  double Track::get_c() const
  {
    return fit->c;
  } 

  double Track::get_d() const
  {
    return fit->d;
  } 

  double Track::get_phi() const
  {
    return fit->phi;
  } 

  double Track::get_r() const
  {
    return fit->r;
  } 

  double Track::get_theta() const
  {
    return fit->theta;
  } 

  double Track::get_h() const
  {
    return fit->h;
  } 
  
  std::vector<Association> & Track::get_associations()
  {
    return associations;
  }

  const std::vector<Association> & Track::get_associations() const
  {
    return associations;
  }

  int Track::hit_split_counter(double t_reference) const
  {
    int count = 0;
    for(const auto & association : associations)
    {
      if(association.parameter_t > t_reference)
      {
        count++;
      }
      else
      {
        count--;
      }
    }
    return count;
  }
  
  Point get_intersection(const TrackHdl & track1, const TrackHdl & track2)
  {
      double a1 = track1->get_a();
      double b1 = track1->get_b();

      double a2 = track2->get_a();
      double b2 = track2->get_b();
        
      if(a1 == a2) return Point();
        
      double x = (b2 - b1) / (a1 - a2);
      double y = a1 * x + b1;

      double z1 = track1->get_c() * x + track1->get_d();
      double z2 = track2->get_c() * x + track2->get_d();

      return Point(x, y, (z1 + z2) / 2.0);	
  }
  
  
  

  double Track::get_chi_squared() const
  {
    return fit->chi_squared;
  }

  double Track::get_chi_squared_R() const
  {
    return fit->chi_squared_R;
  }

  double Track::get_chi_squared_Z() const
  {
    return fit->chi_squared_Z;
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



  double TKtrack::get_quality() const
  {
    return quality;
  } 

  double TKtrack::get_quality_R() const
  {
    return quality_R;
  } 

  double TKtrack::get_quality_Z() const
  {
    return quality_Z;
  } 

  double TKtrack::get_likelihood() const
  {
    return likelihood;
  } 

  double TKtrack::get_likelihood_R() const
  {
    return likelihood_R;
  } 

  double TKtrack::get_likelihood_Z() const
  {
    return likelihood_Z;
  } 
 
*/
  void Track::print(std::ostream & out_) const
  {
    out_ << "Track: "
	 << "a = " << fit->a 
	 << ", b = " << fit->b 
	 << ", c = " << fit->c 
	 << ", d = " << fit->d << std::endl; 
    out_ << "		 phi = " << fit->phi 
	 << ", r = " << fit->r 
	 << ", theta = " << fit->theta 
	 << ", h = " << fit->h << std::endl; 
   /* out_ << "	chi squared: " << fit->chi_squared << std::endl;
    out_ << "	chi squared R: " << fit->chi_squared_R << std::endl;
    out_ << "	chi squared Z: " << fit->chi_squared_Z << std::endl;*/
    out_ << "	number of associated tracker hits: " << associations.size() << std::endl;
  }

} //  end of namespace tkrec

