#ifndef FALAISE_CIMRMAN_POINT_H
#define FALAISE_CIMRMAN_POINT_H

// Standard headers
#include <iostream>
#include <memory>

// Bayeux
#include <datatools/utils.h>

namespace cimrman::datamodel {

  struct Point;
  typedef std::shared_ptr<Point> PointHdl;
  typedef std::shared_ptr<const Point> ConstPointHdl;

  struct Point
  {		
    double x = datatools::invalid_real();
    double y = datatools::invalid_real();
    double z = datatools::invalid_real();
    
    Point() = default;
    Point(double _x, double _y, double _z);
    Point(const std::shared_ptr<const Point>& point);

    void print(std::ostream & out_ = std::clog) const;
  
    // calculates the distance between two points in the horizontal plane
    static double distance_2D(const Point & point1, const Point & point2);
    static double distance_2D(const PointHdl point1, const PointHdl point2);

    // calculates the distance between two points in the 3D space
    static double distance_3D(const Point & point1, const Point & point2);
    static double distance_3D(const PointHdl point1, const PointHdl point2);
    
    // calculates angle (in radians) created by three points (in the horizontal plane)
    static double calculate_angle_2D(const Point & point1, const Point & point2, const Point & point3);
    static double calculate_angle_2D(const PointHdl point1, const PointHdl point2, const PointHdl point3);
    
    // calculates angle (in radians) created by three points (in full 3D space)
    static double calculate_angle_3D(const Point & point1, const Point & point2, const Point & point3);
    static double calculate_angle_3D(const PointHdl point1, const PointHdl point2, const PointHdl point3);
    
    static Point get_middle_point(const Point & point1, const Point & point2);
    static Point get_middle_point(const PointHdl point1, const PointHdl point2);
  };

} //  end of namespace cimrman::datamodel

#endif // FALAISE_CIMRMAN_POINT_H
