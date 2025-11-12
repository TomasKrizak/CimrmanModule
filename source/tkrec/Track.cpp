// Cimrman headers
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
    double sin_phi = std::sin(fit->phi);
    double cos_phi = std::cos(fit->phi);

    double t0 = tracker_hit->get_x() * cos_phi + tracker_hit->get_y() * sin_phi;
    association.parameter_t = t0;
    
    double x0 = t0 * cos_phi + fit->r * sin_phi;
    double y0 = t0 * sin_phi - fit->r * cos_phi;
    double z0 = t0 * std::tan(fit->theta) + fit->h;

     // creating Point obejcts to represent reconstructed avalanche origin points
    association.point = std::make_shared<Point>(x0, y0, z0);

    associations.push_back( association );
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
  
  const LinearFitHdl & Track::get_fit() const
  {
    return fit;
  }
  
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

  void Track::evaluate()
  {
    square_error = 0.0;
    square_error_R = 0.0;
    square_error_Z = 0.0;
    
    chi_squared = 0.0;
    chi_squared_R = 0.0;
    chi_squared_Z = 0.0;
    
    int no_Z = 0;
    int no_R = associations.size();
    
    double sin_phi = std::sin(fit->phi);
    double cos_phi = std::cos(fit->phi);
    for(const auto & association : associations)
    {
      const ConstTrackerHitHdl & hit = association.tracker_hit;
      
      double distance_R2 =  std::abs(fit->r - hit->get_x() * sin_phi + hit->get_y() * cos_phi); 
      distance_R2 = std::pow(distance_R2 - hit->get_R(), 2.0);
      
      square_error_R += distance_R2;
      chi_squared_R += (distance_R2 / std::pow(hit->get_sigma_R(), 2.0));

      if(hit->has_valid_Z()) 
      {
        const ConstPointHdl point = association.point;
        no_Z++;
        double distance_Z2 = std::pow( (hit->get_Z() - point->z) , 2.0);
        square_error_Z += distance_Z2;
        chi_squared_Z += (distance_Z2 / std::pow(hit->get_sigma_Z(), 2.0));  
      }
    }
    square_error = square_error_R + square_error_Z;
    chi_squared = chi_squared_R + chi_squared_Z;     
    
    if(no_Z == 0)
    {
      square_error_Z = datatools::invalid_real();
    }
    return;
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
    return chi_squared;
  }
  
  double Track::get_chi_squared_R() const
  {
    return chi_squared_R;
  }
  
  double Track::get_chi_squared_Z() const
  {
    return chi_squared_Z;
  }
    
  void Track::set_chi_squared(double _chi_squared)
  {
    chi_squared = _chi_squared;
  }
  
  void Track::set_chi_squared_R(double _chi_squared_R)
  {
    chi_squared_R = _chi_squared_R;
  }
  
  void Track::set_chi_squared_Z(double _chi_squared_Z)
  {
    chi_squared_Z = _chi_squared_Z;
  }
  
  double Track::get_square_error() const
  {
    return square_error;
  }
  
  double Track::get_square_error_R() const
  {
    return square_error_R;
  }
  
  double Track::get_square_error_Z() const
  {
    return square_error_Z;
  }
  
  void Track::set_square_error(double _square_error)
  {
    square_error = _square_error;
  }
  
  void Track::set_square_error_R(double _square_error_R)
  {
    square_error_R = _square_error_R;
  }
  
  void Track::set_square_error_Z(double _square_error_Z)
  {
    square_error_Z = _square_error_Z;
  }
  
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
    out_ << "	chi squared: " << chi_squared << std::endl;
    out_ << "	chi squared R: " << chi_squared_R << std::endl;
    out_ << "	chi squared Z: " << chi_squared_Z << std::endl;
    out_ << "	number of associated tracker hits: " << associations.size() << std::endl;
  }

} //  end of namespace tkrec

