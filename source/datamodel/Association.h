#ifndef FALAISE_CIMRMAN_ASSOCIATION_H
#define FALAISE_CIMRMAN_ASSOCIATION_H

// Standard headers
#include<memory>

namespace cimrman::datamodel {
  
  // Forward declarations
  class TrackerHit;
  using ConstTrackerHitHdl = std::shared_ptr<const TrackerHit>;

  class Point;
  using PointHdl = std::shared_ptr<Point>;


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
  
  struct Association
  {
    ConstTrackerHitHdl tracker_hit;
    PointHdl point;
    double parameter_t;
    
    Association(const ConstTrackerHitHdl & hit);
  };

} //  end of namespace cimrman::datamodel

#endif // FALAISE_CIMRMAN_ASSOCIATION_H
