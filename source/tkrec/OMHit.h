#ifndef FALAISE_CIMRMAN_OMHIT_H
#define FALAISE_CIMRMAN_OMHIT_H

// Standard headers
#include <iostream>
#include <memory>

#include <datatools/utils.h>

namespace tkrec {

  class OMHit
  {
  
    int OM_num = -1;
    int SWCR[4]; // 0 = Side, 1 = Wall, 2 = Column, 3 = Row
    double x = datatools::invalid_real();
    double y = datatools::invalid_real();
    double z = datatools::invalid_real();

  public:
  		
    OMHit(int _OM_num);
    OMHit(int _SWCR[4]);
    virtual ~OMHit() = default;
		
    void set_OM_num(int _OM_num);	
    int get_OM_num() const;
    
    int* get_SWCR();
    const int* get_SWCR() const;
    
    double get_x() const;
    double get_y() const;
    double get_z() const;
    
    void set_x(double _x);
    void set_y(double _y);
    void set_z(double _z);		

    void print(std::ostream & out_ = std::cout) const;
  };

  typedef std::shared_ptr<OMHit> OMHitHdl;
  typedef std::shared_ptr<const OMHit> ConstOMHitHdl;
  
} //  end of namespace tkrec

#endif // FALAISE_CIMRMAN_OMHIT_H
