#ifndef FALAISE_CIMRMAN_TRAJECTORY_H
#define FALAISE_CIMRMAN_TRAJECTORY_H

// Standard headers
#include <iostream>
#include <vector>
#include <memory>

// Bayeux headers
#include <datatools/utils.h>

// Cimrman headers
#include "tkrec/Association.h"

namespace tkrec {

  // forward definitions
  class Track;
  using TrackHdl = std::shared_ptr<Track>;
  using ConstTrackHdl = std::shared_ptr<const Track>;
  
  class Point;
  using PointHdl = std::shared_ptr<Point>;
  using ConstPointHdl = std::shared_ptr<const Point>;
  
  
  class Trajectory
  {
  public:
  
    enum EndPoint{ FRONT, BACK };
    enum Type{ ELECTRON, ALPHA, ANY };

  private:
		
    std::vector<PointHdl> trajectory_points; // startpoint, kink points and endpoint of trajectory
    std::vector<TrackHdl> segments; // linear parts of the polyline trajectory
    bool kinked_trajectory = false;
    
    // chi_squared over NDF (number of measurements - number of parameters of the polyline)
    // number of measurements: number of associated tracker hits with available R
    //                       + number of associated tracker hits with available Z
    // parameters of trajectory: 4 for a straight trajectory
    //                           4 + 3k for kinked trajectory where k is number of kinks
    double chi_squared = datatools::invalid_real();
    double chi_squared_R = datatools::invalid_real();
    double chi_squared_Z = datatools::invalid_real();
  
    // mean square error  
    double MSE = datatools::invalid_real();
    double MSE_R = datatools::invalid_real();
    double MSE_Z = datatools::invalid_real();
    
  public:
		
    Trajectory() = default;
    Trajectory(TrackHdl & segment);
    Trajectory(std::vector<TrackHdl> & _segments);
    virtual ~Trajectory() = default;
		
    std::vector<PointHdl> & get_trajectory_points();
    std::vector<ConstPointHdl> get_trajectory_points() const;
    
    std::vector<TrackHdl> & get_segments();
    std::vector<ConstTrackHdl> get_segments() const;
    
    const std::vector<Association> get_associations() const;

    void mark_as_kinked();
    void mark_as_straight();
    bool has_kink() const;

    void update_segments();

    double get_chi_squared() const;
    double get_chi_squared_R() const;
    double get_chi_squared_Z() const;
    
    void set_chi_squared(double _chi_squared);
    void set_chi_squared_R(double _chi_squared_R);
    void set_chi_squared_Z(double _chi_squared_Z);
    
    double get_MSE() const;
    double get_MSE_R() const;
    double get_MSE_Z() const;
    
    void set_MSE(double _MSE);
    void set_MSE_R(double _MSE_R);
    void set_MSE_Z(double _MSE_Z);

    void print(std::ostream & out_ = std::clog) const;
		    
  };
  
  typedef std::shared_ptr<Trajectory> TrajectoryHdl;
  typedef std::shared_ptr<const Trajectory> ConstTrajectoryHdl;

} //  end of namespace tkrec

#endif // FALAISE_CIMRMAN_TRAJECTORY_H
