// TK headers
#include "tkrec/Point.h"

// ClassImp(tkrec::Point);

namespace tkrec {

  using namespace std;

  Point::Point(double _x, double _y, double _z)
    : x(_x)
    , y(_y)
    , z(_z)
  {
  }
  
  Point::Point(const std::shared_ptr<const Point>& point)
     : x(point->x)
     , y(point->y)
     , z(point->z) 
  {
  }

  void Point::print(std::ostream & out_) const
  {
    out_ << "point (" << x << ", " << y << ", " << z << ")" << std::endl;
  }

  double distance_2D(const PointHdl & point1, const PointHdl & point2)
  {
    double temp = pow(point1->x - point2->x, 2) + pow(point1->y - point2->y, 2);
    return sqrt(temp);
  }

  double distance_3D(const PointHdl & point1, const PointHdl & point2)
  {
    double temp = pow(point1->x - point2->x, 2) + pow(point1->y - point2->y, 2) + pow(point1->z - point2->z, 2);
    return sqrt(temp);
  }

  double distance_2D(const Point & point1, const Point & point2)
  {
    double temp = pow(point1.x - point2.x, 2) + pow(point1.y - point2.y, 2);
    return sqrt(temp);
  }

  double distance_3D(const Point & point1, const Point & point2)
  {
    double temp = pow(point1.x - point2.x, 2) + pow(point1.y - point2.y, 2) + pow(point1.z - point2.z, 2);
    return sqrt(temp);
  }


} //  end of namespace tkrec
