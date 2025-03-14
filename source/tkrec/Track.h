#ifndef FALAISE_TKRECONSTRUCT_TRACK_H
#define FALAISE_TKRECONSTRUCT_TRACK_H

// Standard headers
#include <iostream>
#include <vector>
#include <cmath>
#include <memory>
#include <algorithm>

#include "tkrec/TrackerHit.h"
#include "tkrec/LinearFit.h"
#include "tkrec/Point.h"

#include <datatools/utils.h>

namespace tkrec {

  // Track parametrised by LinearFit in two alternative ways:
  // line in the form: (does not describe lines parallel to source foil)
  //	y = ax + b
  //	z = cx + d

  // or in the form:
  // 	x(t) = cos(phi)*cos(theta)*t + r*sin(phi)
  // 	y(t) = sin(phi)*cos(theta)*t - r*cos(phi)
  // 	z(t) = sin(theta)*t + h

  // set of transformations:
  //	a = tan(phi)
  //	b = -r/cos(phi)
  //	c = tan(theta)/cos(phi)
  //	d = h - r*tan(phi)*tan(theta)

  //	phi = atan(a) 
  //	r = -b/sqrt(a*a+1.0)
  //	theta = atan(c/sqrt(a*a+1.0))
  //	h = d - a*b*c/(a*a+1.0)
  
  struct Association
  {
    ConstTrackerHitHdl tracker_hit;
    PointHdl point;
    double parameter_t;
    
    Association(const ConstTrackerHitHdl & hit);
  };
  
  // Association stores relevant information about Track - TrackerHit connection
    
  // association point: 
  // closest point on the track to a tracekr hit
  // (tracker hit - anode wire (xi,yi), vertical position zi)
  // 	x(t) = (xi*cos(phi)*+yi*sin(phi))*cos(phi) + r*sin(phi)
  // 	y(t) = (xi*cos(phi)*+yi*sin(phi))*sin(phi) - r*cos(phi)
  // 	z(t) = (xi*cos(phi)*+yi*sin(phi))*tan(theta) + h - zi

  // association parameter t:
  // t describes the position of the association along the track
  // it is calculated as a projection of any point (x,y) onto the track
  // t = x * cos(phi) + y * sin(phi)
  
  class Track
  {
  private:
  
    LinearFitHdl fit;
    std::vector<Association> associations;
	
  public:
	
    Track() = default;
    Track(const ConstLinearFitHdl & linear_fit, const std::vector<ConstTrackerHitHdl> & tracker_hits);
    virtual ~Track() = default;
		
    void add_associated_tracker_hit(ConstTrackerHitHdl & tracker_hit);
				    
    std::vector<Association> & get_associations(); 
    const std::vector<Association> & get_associations() const; 
    
    void sort_associations();
    double calculate_t(double x, double y) const;   
    double calculate_t(const Point & point) const;   
    double calculate_t(ConstPointHdl & point) const;
    int hit_split_counter(double t_reference) const;
    double horizontal_distance_to_line(const Point & point) const;
  
    void update_associations();
  
    void set_a(double _a);
    void set_b(double _b);
    void set_c(double _c);
    void set_d(double _d);
		
    void set_phi(double _phi);
    void set_r(double _r);
    void set_theta(double _theta);
    void set_h(double _h);

    double get_a() const;
    double get_b() const;
    double get_c() const;
    double get_d() const;
		
    double get_phi() const;
    double get_r() const;
    double get_theta() const;
    double get_h() const;

    double get_chi_squared() const;
    double get_chi_squared_R() const;
    double get_chi_squared_Z() const;
		
	/*
    void calculate_tr_hit_points();
			
    void set_chi_squared(double _chi_squared);
    void set_chi_squared_R(double _chi_squared_R);
    void set_chi_squared_Z(double _chi_squared_Z);
		
    void set_quality(double _quality);
    void set_quality_R(double _quality_R);
    void set_quality_Z(double _quality_Z);
		
    void set_likelihood(double _likelihood);
    void set_likelihood_R(double _likelihood_R);
    void set_likelihood_Z(double _likelihood_Z);
		
    double get_quality() const;
    double get_quality_R() const;
    double get_quality_Z() const;
		
    double get_likelihood() const;
    double get_likelihood_R() const;
    double get_likelihood_Z() const;
*/
    void print(std::ostream & out_ = std::cout) const;
    
  };
  
  typedef std::shared_ptr<Track> TrackHdl;
  typedef std::shared_ptr<const Track> ConstTrackHdl;
  
  Point get_intersection(const TrackHdl & track1, const TrackHdl & track2); 

} //  end of namespace tkrec

#endif // FALAISE_TKRECONSTRUCT_TRACK_H
