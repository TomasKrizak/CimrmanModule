// Cimrman headers
#include "tkrec/Trajectory.h"

// Bayeux:
#include <datatools/exception.h>

// ClassImp(tkrec::Trajectory);

namespace tkrec {


  Trajectory::Trajectory(TrackHdl & segment)
    : Trajectory()
  {
    segments.push_back(segment);
    return;
  }


  Trajectory::Trajectory(std::vector<TrackHdl> & _segments)
    : Trajectory()
  {
    segments = _segments;
    if( segments.size() > 1)
    {
      kinked_trajectory = true;
    }
    return;
  }

  std::vector<TrackHdl> & Trajectory::get_segments()
  {
    return segments;
  }

  std::vector<ConstTrackHdl> Trajectory::get_segments() const
  {
    std::vector<ConstTrackHdl> segs;
    for(const auto & track : segments) {
      segs.push_back(track);
    }
    return segs;
  }

  std::vector<PointHdl> & Trajectory::get_trajectory_points()
  {
     return trajectory_points;
  }

  std::vector<ConstPointHdl> Trajectory::get_trajectory_points() const
  {
    std::vector<ConstPointHdl> traj_points;
    for(const auto & point : trajectory_points)
    {
      traj_points.push_back(point);
    }
    return traj_points;
  }

  const std::vector<Association> Trajectory::get_associations() const
  {
    std::vector<Association> all_associations;
    for(const auto & track : segments)
    {
      const std::vector<Association> associations = track->get_associations();
      all_associations.insert(all_associations.end(), associations.begin(), associations.end());
    }
    return all_associations;
  }

  void Trajectory::mark_as_kinked()
  {
    kinked_trajectory = true;
  }
  
  void Trajectory::mark_as_straight()
  {
    kinked_trajectory = false; 
  }

  bool Trajectory::has_kink() const
  {
    return kinked_trajectory;
  }
  
  void Trajectory::update_segments()
  {
    for(auto i = 0u; i < segments.size(); i++)
    {
      ConstPointHdl point1 = trajectory_points[i];
      ConstPointHdl point2 = trajectory_points[i+1];
      
      TrackHdl & segment = segments[i];
      double a = (point2->y - point1->y) / (point2->x - point1->x);
      double b = point1->y - a * point1->x;
      double c = (point2->z - point1->z) / (point2->x - point1->x);
      double d = point1->z - c * point1->x;
      
      segment->set_a(a);
      segment->set_b(b);
      segment->set_c(c);
      segment->set_d(d);
      
      segment->set_phi( std::atan(a) );
      segment->set_r( -b / std::sqrt(a * a + 1.0) );
      segment->set_theta( std::atan(c / std::sqrt(a * a + 1.0)) );
      segment->set_h( d - (a * b * c) / (a * a + 1.0) );
      
      segment->update_associations();
    }
    
    return;
  }
  
  double Trajectory::get_chi_squared() const
  {
    return chi_squared;
  }
  
  double Trajectory::get_chi_squared_R() const
  {
    return chi_squared_R;
  }
  
  double Trajectory::get_chi_squared_Z() const
  {
    return chi_squared_Z;
  }
    
  void Trajectory::set_chi_squared(double _chi_squared)
  {
    chi_squared = _chi_squared;
  }
  
  void Trajectory::set_chi_squared_R(double _chi_squared_R)
  {
    chi_squared_R = _chi_squared_R;
  }
  
  void Trajectory::set_chi_squared_Z(double _chi_squared_Z)
  {
    chi_squared_Z = _chi_squared_Z;
  }
  
  double Trajectory::get_MSE() const
  {
    return MSE;
  }
  
  double Trajectory::get_MSE_R() const
  {
    return MSE_R;
  }
  
  double Trajectory::get_MSE_Z() const
  {
    return MSE_Z;
  }
  
  void Trajectory::set_MSE(double _MSE)
  {
    MSE = _MSE;
  }
  
  void Trajectory::set_MSE_R(double _MSE_R)
  {
    MSE_R = _MSE_R;
  }
  
  void Trajectory::set_MSE_Z(double _MSE_Z)
  {
    MSE_Z = _MSE_Z;
  }
  
  void Trajectory::print(std::ostream & out_) const
  {
    out_ <<"Trajectory: " << std::endl;
    out_ << "	" << segments.size() << " segments: " << std::endl;
    for(auto i = 0u; i < segments.size(); i++)
    {
      out_ << "	";
      segments[i]->print();
    }
    out_ << "	" << trajectory_points.size() << " trajectory points: " << std::endl;
    for(auto i = 0u; i < trajectory_points.size(); i++)
    {
      out_ << "	" << i+1 << ". ";
      trajectory_points[i]->print(); 
    }
    return;
  }

} //  end of namespace tkrec
