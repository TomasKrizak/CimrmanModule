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
    
    // calculates the point between the two points
    static Point get_middle_point(const Point & point1, const Point & point2);
    static Point get_middle_point(const PointHdl point1, const PointHdl point2);
    
    // vector helpers:
    
    // treats the 3D points as vectors and calculates their cross product
    static Point cross_product(const Point & point1, const Point & point2);
    static Point cross_product(const PointHdl point1, const PointHdl point2);
    
    // treats the 3D points as vectors and calculates their dot product
    static double dot_product(const Point & point1, const Point & point2);
    static double dot_product(const PointHdl point1, const PointHdl point2);
    
    Point operator+(const Point & point1) const;
    Point operator-(const Point & point1) const;
    Point operator*(const double factor) const;
    Point operator/(const double factor) const;
    
  };

} //  end of namespace cimrman::datamodel

#endif // FALAISE_CIMRMAN_POINT_H
