#ifndef FALAISE_CIMRMAN_POINT_H
#define FALAISE_CIMRMAN_POINT_H

// Standard headers
#include <iostream>
#include <cmath>
#include <memory>

#include <datatools/utils.h>

namespace tkrec {

  struct Point
  {		
    double x = datatools::invalid_real();
    double y = datatools::invalid_real();
    double z = datatools::invalid_real();
    
    Point() = default;
    Point(double _x, double _y, double _z);
    Point(const std::shared_ptr<const Point>& point);

    void print(std::ostream & out_ = std::clog) const;

  };

  typedef std::shared_ptr<Point> PointHdl;
  typedef std::shared_ptr<const Point> ConstPointHdl;

  double distance_2D(const PointHdl & point1, const PointHdl & point2);
  double distance_3D(const PointHdl & point1, const PointHdl & point2);
  
  double distance_2D(const Point & point1, const Point & point2);
  double distance_3D(const Point & point1, const Point & point2);

} //  end of namespace tkrec

#endif // FALAISE_CIMRMAN_POINT_H
