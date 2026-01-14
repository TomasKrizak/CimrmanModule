// Cimrman headers
#include "tkrec/Point.h"

// Standard headers
#include <cmath>

// ClassImp(tkrec::Point);

namespace tkrec {

  using namespace std;

  Point::Point(double _x, double _y, double _z)
    : x(_x), y(_y), z(_z)
  {
  }
  
  Point::Point(const std::shared_ptr<const Point>& point)
     : x(point->x), y(point->y), z(point->z) 
  {
  }

  void Point::print(std::ostream & out_) const
  {
    out_ << "point (" << x << ", " << y << ", " << z << ")" << std::endl;
  }

  double distance_2D(const PointHdl & point1, const PointHdl & point2)
  {
    return std::hypot(point1->x - point2->x, point1->y - point2->y);
  }

  double distance_3D(const PointHdl & point1, const PointHdl & point2)
  {
    return std::hypot(point1->x - point2->x, point1->y - point2->y, point1->z - point2->z);
  }

  double distance_2D(const Point & point1, const Point & point2)
  {
    return std::hypot(point1.x - point2.x, point1.y - point2.y);
  }

  double distance_3D(const Point & point1, const Point & point2)
  {
    return std::hypot(point1.x - point2.x, point1.y - point2.y, point1.z - point2.z);
  }


} //  end of namespace tkrec
