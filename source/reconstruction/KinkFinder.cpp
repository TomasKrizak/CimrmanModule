// Cimrman headers
#include "reconstruction/KinkFinder.h"
#include "reconstruction/Algos.h"
#include "geometry/Geometry.h"
#include "datamodel/TrackerHit.h"
#include "datamodel/Track.h"

// Standard headers
#include <cmath>
#include <tuple>
//#include <iostream>

using namespace cimrman::datamodel;

namespace cimrman{
    
  KinkFinder::KinkFinder( TrackHdl _track1, PointHdl _endpoint1,
                          TrackHdl _track2, PointHdl _endpoint2,
                          const PolylinesConfig & _config,
                          const Geometry & _geom)
    : track1( _track1),
      endpoint1(_endpoint1),
      track2( _track2),
      endpoint2(_endpoint2),
      config(_config),
      geom(_geom) 
  {
    return;
  }
      
  // passes if the 2D angular deviation is in bounds
  bool KinkFinder::check_angle_2D_after( const double min_angle, const double max_angle ) const
  {
    if( status == Status::UNPROCESSED ) return false;
       
    // we need to calculate angle between the kink point and the opposite ends of the tracks
    Point opposite_end1;
    const auto & associations1 = track1->get_associations();
    
    // XXX I have no idea what is a good tolerance for points mismatch - but 1mm works fine
    if( Point::distance_3D( associations1.front().point, endpoint1 ) < 1.0 )
    {
      opposite_end1 = *(associations1.back().point);
    }
    else
    {
      opposite_end1 = *(associations1.front().point);
    }
    
    Point opposite_end2;
    const auto & associations2 = track2->get_associations();
    if( Point::distance_3D( associations2.front().point, endpoint2 ) < 1.0 )
    {
      opposite_end2 = *(associations2.back().point);
    }
    else
    {
      opposite_end2 = *(associations2.front().point);
    }
       
    const double angle = Point::calculate_angle_2D(opposite_end1, kink_point, opposite_end2);
    const double angular_deviation = M_PI - angle;
    
    return (min_angle <= angular_deviation && angular_deviation <= max_angle);
  }
  
  // passes if the 3D angular deviation is in bounds
  bool KinkFinder::check_angle_3D_after( const double min_angle, const double max_angle ) const
  {
    if( status == Status::UNPROCESSED ) return false;
    
    // we need to calculate angle between the kink point and the opposite ends of the tracks
    Point opposite_end1;
    const auto & associations1 = track1->get_associations();
    if( Point::distance_3D( associations1.front().point, endpoint1 ) < 1.0 )
    {
      opposite_end1 = *(associations1.back().point);
    }
    else
    {
      opposite_end1 = *(associations1.front().point);
    }
    
    Point opposite_end2;
    const auto & associations2 = track2->get_associations();
    if( Point::distance_3D( associations2.front().point, endpoint2 ) < 1.0 )
    {
      opposite_end2 = *(associations2.back().point);
    }
    else
    {
      opposite_end2 = *(associations2.front().point);
    }

    const double angle = Point::calculate_angle_3D(opposite_end1, kink_point, opposite_end2);
    const double angular_deviation = M_PI - angle;
    
    return (min_angle <= angular_deviation && angular_deviation <= max_angle);
  } 

  // passes if both tracks have sufficiently close hits to the kink point 
  bool KinkFinder::check_close_hits_2D( const double max_distance ) const
  {
    if( status == Status::UNPROCESSED ) return false;
    
    const auto & associations1 = track1->get_associations();
    const auto & associations2 = track2->get_associations();
    
    const auto & kink_point = this->kink_point;
    auto is_close = [&kink_point, max_distance](const auto & association)
    {
      const double distance = Point::distance_2D( kink_point, *(association.point) ); 
      return distance < max_distance;
    };
    
    const auto close_hit1 = (std::find_if(associations1.begin(), associations1.end(), is_close));
    bool has_close_hit1 = (close_hit1 != associations1.end()); 
    
    const auto close_hit2 = std::find_if(associations2.begin(), associations2.end(), is_close);
    bool has_close_hit2 = (close_hit2 != associations2.end()); 
   
    return (has_close_hit1 && has_close_hit2);
  }
  
  // passes if the 3D angular deviation is in bounds
  bool KinkFinder::check_angle_3D_before( const double min_angle, const double max_angle ) const
  { 
    // we need to calculate angle between the kink point and the opposite ends of the tracks
    Point opposite_end1;
    const auto & associations1 = track1->get_associations();
    if( Point::distance_3D( associations1.front().point, endpoint1 ) < 1.0 )
    {
      opposite_end1 = *(associations1.back().point);
    }
    else
    {
      opposite_end1 = *(associations1.front().point);
    }
    
    Point opposite_end2;
    const auto & associations2 = track2->get_associations();
    if( Point::distance_3D( associations2.front().point, endpoint2 ) < 1.0 )
    {
      opposite_end2 = *(associations2.back().point);
    }
    else
    {
      opposite_end2 = *(associations2.front().point);
    }
    
    Point direction1 = {endpoint1->x - opposite_end1.x,
                        endpoint1->y - opposite_end1.y,
                        endpoint1->z - opposite_end1.z};
                        
    Point direction2 = {endpoint2->x - opposite_end2.x,
                        endpoint2->y - opposite_end2.y,
                        endpoint2->z - opposite_end2.z};
    
    const double angle = Point::calculate_angle_3D(direction1, {0.0, 0.0, 0.0}, direction2);
    const double angular_deviation = M_PI - angle;
    
    return (min_angle <= angular_deviation && angular_deviation <= max_angle);
  } 
  
  
  // passes if the 3D angular deviation is in bounds
  bool KinkFinder::check_angle_2D_before( const double min_angle, const double max_angle ) const
  { 
    // we need to calculate angle between the kink point and the opposite ends of the tracks
    Point opposite_end1;
    const auto & associations1 = track1->get_associations();
    if( Point::distance_3D( associations1.front().point, endpoint1 ) < 1.0 )
    {
      opposite_end1 = *(associations1.back().point);
    }
    else
    {
      opposite_end1 = *(associations1.front().point);
    }
    
    Point opposite_end2;
    const auto & associations2 = track2->get_associations();
    if( Point::distance_3D( associations2.front().point, endpoint2 ) < 1.0 )
    {
      opposite_end2 = *(associations2.back().point);
    }
    else
    {
      opposite_end2 = *(associations2.front().point);
    }
    
    Point direction1 = {endpoint1->x - opposite_end1.x,
                        endpoint1->y - opposite_end1.y,
                        endpoint1->z - opposite_end1.z};
                        
    Point direction2 = {endpoint2->x - opposite_end2.x,
                        endpoint2->y - opposite_end2.y,
                        endpoint2->z - opposite_end2.z};
    
    const double angle = Point::calculate_angle_2D(direction1, {0.0, 0.0, 0.0}, direction2);
    const double angular_deviation = M_PI - angle;
    
    return (min_angle <= angular_deviation && angular_deviation <= max_angle);
  } 
  
  
  // passes if the kink point is inside the tracker volume
  bool KinkFinder::check_is_inside_tracker() const
  {
    if( status == Status::UNPROCESSED ) return false;
    
    return geom.is_inside_tracker(kink_point);
  }
  
  // passes if kink point is more than min_distance from mainwall
  bool KinkFinder::check_distance_to_MW( const double min_distance ) const
  { 
    if( status == Status::UNPROCESSED ) return false;
    
    const double distance = geom.distance_to_MW(kink_point);
    return distance > min_distance;
  }
  
  // passes if kink point is more than min_distance from Xwall
  bool KinkFinder::check_distance_to_XW( const double min_distance ) const
  {
    if( status == Status::UNPROCESSED ) return false;
    
    double distance = geom.distance_to_XW(kink_point);
    return distance > min_distance;
  }
  
  // passes if kink point is more than min_distance from source foil
  bool KinkFinder::check_distance_to_SF( const double min_distance ) const
  {
    if( status == Status::UNPROCESSED ) return false;
    
    const double distance = geom.distance_to_SF(kink_point);
    return distance > min_distance;
  }

  // passes if the ends of the tracks are close enough in 2D
  bool KinkFinder::check_are_ends_close_2D( const double max_distance ) const
  {
    if( status == Status::UNPROCESSED ) return false;
    
    const double distance_2D = Point::distance_2D( endpoint1, endpoint2 ); 
    return distance_2D < max_distance;
  }
  
  // passes if the ends of the tracks are close enough in 3D
  bool KinkFinder::check_are_ends_close_3D( const double max_distance ) const
  {
    if( status == Status::UNPROCESSED ) return false;
    
    const double distance_3D = Point::distance_3D( endpoint1, endpoint2 ); 
    return distance_3D < max_distance;
  }
  
  // designed for large kink procedure
  // passes if the vertical mismatch of the track in the kink point is small enough 
  bool KinkFinder::check_vertical_distance( const double max_distance) const
  {
    if( status == Status::UNPROCESSED ) return false;
    
    const double z1 = track1->get_c() * kink_point.x + track1->get_d();
    const double z2 = track2->get_c() * kink_point.x + track2->get_d();
    
    return std::abs(z2 - z1) < max_distance;		
  }
  
  // designed for small kink procedure
  // passes if the kink point is close enough to both lineear fits (not segments)
  bool KinkFinder::check_middle_point_distance( const double max_distance) const
  {
    if( status == Status::UNPROCESSED ) return false;
       
    double distance1 = track1->horizontal_distance_to_line(kink_point); 
    double distance2 = track2->horizontal_distance_to_line(kink_point);  
               
    return (distance1 < max_distance) && (distance2 < max_distance); 
  }
   
   
  void KinkFinder::process(ConnectionStrategy strategy)
  {
    if( strategy == ConnectionStrategy::VERTICAL_ALIGNMENT ) 
    {
      vertical_alignment_procedure();
    }
    else if( strategy == ConnectionStrategy::ENDPOINTS_MIDDLE ) 
    {
      endpoints_middle_procedure();
    }
  }
  
  bool KinkFinder::check_criteria(std::initializer_list<Criterion> criteria) const
  {
    for(const auto & criterium : criteria)
    {
      bool passed = criterium();
      if( not passed )
      {
        return false;
      }
    }
    return true;
  }
  
  void KinkFinder::vertical_alignment_procedure()
  {
    DT_THROW_IF(endpoint1 == nullptr || endpoint2 == nullptr, std::logic_error, 
      "Missing endpoints for kink investigation!");
    
    // connection strategy
    kink_point = Track::get_intersection_2D(track1, track2);    
    status = Status::CANDIDATE_LOCATED;
    
    // constructing a list of criteria
    std::initializer_list<Criterion> criteria = {
      [&]{ return check_are_ends_close_2D( config.max_trajectory_endpoints_distance ); },
      [&]{ return check_close_hits_2D( config.max_tracker_hits_distance ); },
      [&]{ return check_angle_3D_after( config.min_kink_angle, config.max_kink_angle ); },
      [&]{ return check_angle_3D_before( config.min_kink_angle, config.max_kink_angle ); },
      [&]{ return check_vertical_distance( config.max_vertical_distance ); },
      [&]{ return check_is_inside_tracker(); },
      [&]{ return check_distance_to_SF( config.min_distance_from_foil ); },
      [&]{ return check_distance_to_MW( config.min_distance_from_main_walls ); },
      [&]{ return check_distance_to_XW( config.min_distance_from_X_walls ); }
    };
    
    // applying criteria
    bool passed = check_criteria( criteria );
    if( passed ) 
    {
      status = Status::KINK_FOUND;
    }
    else 
    {
      status = Status::KINK_NOT_FOUND;
    }
  }
  
  void KinkFinder::endpoints_middle_procedure()
  {
    DT_THROW_IF(endpoint1 == nullptr || endpoint2 == nullptr, std::logic_error,
      "Missing endpoints for kink investigation!");

    // connection strategy
    kink_point = Point::get_middle_point(endpoint1, endpoint2);
    status = Status::CANDIDATE_LOCATED;

    // constructing a list of criteria
    std::initializer_list<Criterion> criteria = {
      [&]{ return check_are_ends_close_3D( config.max_trajectory_endpoints_distance ); },
      [&]{ return check_close_hits_2D( config.max_tracker_hits_distance ); },
      [&]{ return check_middle_point_distance( config.max_trajectories_middlepoint_distance ); },
      [&]{ return check_angle_3D_after( config.min_trajectory_connection_angle, config.max_trajectory_connection_angle ); },
    };
    
    // applying criteria
    bool passed = check_criteria( criteria );
    if( passed )
    {
      status = Status::KINK_FOUND;
    }
    else
    {
      status = Status::KINK_NOT_FOUND;
    }
  }
  
  KinkFinder::Status KinkFinder::get_status() const
  {
    return status;
  }
  
  Point KinkFinder::get_kink_point() const
  {
    DT_THROW_IF(status == Status::UNPROCESSED, std::logic_error, 
      "Kink point not constructed!");
    
    return kink_point;
  }

} //  end of namespace cimrman
