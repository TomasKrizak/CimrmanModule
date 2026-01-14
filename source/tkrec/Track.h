#ifndef FALAISE_CIMRMAN_TRACK_H
#define FALAISE_CIMRMAN_TRACK_H

// Standard headers
#include <iostream>
#include <vector>
#include <memory>

// Bayeux headers
#include <datatools/utils.h>

// Cimrman headers
#include "tkrec/Association.h"

namespace tkrec {

  // Forward declaration 
  class TrackerHit;
  using TrackerHitHdl = std::shared_ptr<TrackerHit>;
  using ConstTrackerHitHdl = std::shared_ptr<const TrackerHit>;

  class Point;
  using PointHdl = std::shared_ptr<Point>;
  using ConstPointHdl = std::shared_ptr<const Point>;
  
  class LinearFit;
  using LinearFitHdl = std::shared_ptr<LinearFit>;
  using ConstLinearFitHdl = std::shared_ptr<const LinearFit>;


  class Track
  {
  private:
  
    // Track geometrically described in LinearFit class
    LinearFitHdl fit;
    std::vector<Association> associations;
    
    double chi_squared = datatools::invalid_real();
    double chi_squared_R = datatools::invalid_real();
    double chi_squared_Z = datatools::invalid_real();
    
    double square_error = datatools::invalid_real();
    double square_error_R = datatools::invalid_real();
    double square_error_Z = datatools::invalid_real();

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
    void evaluate();
    
    const LinearFitHdl & get_fit() const;
    
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
    
    void set_chi_squared(double _chi_squared);
    void set_chi_squared_R(double _chi_squared_R);
    void set_chi_squared_Z(double _chi_squared_Z);
    
    double get_square_error() const;
    double get_square_error_R() const;
    double get_square_error_Z() const;
    
    void set_square_error(double _square_error);
    void set_square_error_R(double _square_error_R);
    void set_square_error_Z(double _square_error_Z);
    
    void print(std::ostream & out_ = std::cout) const;
    
  };
  
  typedef std::shared_ptr<Track> TrackHdl;
  typedef std::shared_ptr<const Track> ConstTrackHdl;
  
  Point get_intersection(const TrackHdl & track1, const TrackHdl & track2); 

} //  end of namespace tkrec

#endif // FALAISE_CIMRMAN_TRACK_H
