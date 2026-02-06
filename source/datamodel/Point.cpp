// Cimrman headers
#include "datamodel/Point.h"

// Standard headers
#include <cmath>
#include <tuple>

namespace cimrman::datamodel {

  Point::Point(double _x, double _y, double _z)
    : x(_x), y(_y), z(_z)
  {
    return;
  }
  
  Point::Point(const std::shared_ptr<const Point>& point)
     : x(point->x), y(point->y), z(point->z) 
  {
    return;
  }

  void Point::print(std::ostream & out_) const
  {
    out_ << "point (" << x << ", " << y << ", " << z << ")" << std::endl;
  }

  double Point::distance_2D(const PointHdl point1, const PointHdl point2)
  {
    return std::hypot(point1->x - point2->x, point1->y - point2->y);
  }

  double Point::distance_3D(const PointHdl point1, const PointHdl point2)
  {
    return std::hypot(point1->x - point2->x, point1->y - point2->y, point1->z - point2->z);
  }

  double Point::distance_2D(const Point & point1, const Point & point2)
  {
    return std::hypot(point1.x - point2.x, point1.y - point2.y);
  }

  double Point::distance_3D(const Point & point1, const Point & point2)
  {
    return std::hypot(point1.x - point2.x, point1.y - point2.y, point1.z - point2.z);
  }
 
  double Point::calculate_angle_3D(const Point & point1, const Point & point2, const Point & point3)
  {
    // calculating 3D angle between 2 vectors
    const Point vec1 = {point1.x - point2.x,
                        point1.y - point2.y,
                        point1.z - point2.z };
                  
    const Point vec2 = {point3.x - point2.x,
                        point3.y - point2.y,
                        point3.z - point2.z };

    const double norm1 = std::hypot(vec1.x, vec1.y, vec1.z);
    const double norm2 = std::hypot(vec2.x, vec2.y, vec2.z);
    const double scalar_prod = (vec1.x * vec2.x) + (vec1.y * vec2.y) + (vec1.z * vec2.z);
    double angle = scalar_prod / (norm1 * norm2);
    angle = std::max(-1.0, std::min(1.0, angle));
    return std::acos(angle);
  } 
  
  double Point::calculate_angle_3D(const PointHdl point1, const PointHdl point2, const PointHdl point3)
  {
    return calculate_angle_3D( *point1, *point2, *point3 );
  }
  
  double Point::calculate_angle_2D(const Point & point1, const Point & point2, const Point & point3)
  {
     // calculating 2D angle between 2 vectors (in the horizontal plane)
    const std::pair<double, double> vec1 = {point1.x - point2.x,
                                            point1.y - point2.y};
                  
    const std::pair<double, double> vec2 = {point3.x - point2.x,
                                            point3.y - point2.y};
                  
    const double norm1 = std::hypot(vec1.first, vec1.second);
    const double norm2 = std::hypot(vec2.first, vec2.second);
    
    const double scalar_prod = (vec1.first * vec2.second) + (vec1.first * vec2.second);
    double angle = scalar_prod / (norm1 * norm2);
    angle = std::max(-1.0, std::min(1.0, angle));
    return angle = std::acos(angle);
  }
  
  double Point::calculate_angle_2D(const PointHdl point1, const PointHdl point2, const PointHdl point3)
  {
    return calculate_angle_2D( *point1, *point2, *point3 );
  }
  
  Point Point::get_middle_point(const Point & point1, const Point & point2)
  {
    Point middle_point;
    middle_point.x = (point1.x + point2.x) / 2.0;
    middle_point.y = (point1.y + point2.y) / 2.0;
    middle_point.z = (point1.z + point2.z) / 2.0;
    return middle_point;
  }   
  
  Point Point::get_middle_point(const PointHdl point1, const PointHdl point2)
  {
    return get_middle_point(*point1, *point2);
  }   

} //  end of namespace cimrman::datamodel
