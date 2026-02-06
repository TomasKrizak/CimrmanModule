// Cimrman headers
#include "datamodel/OMHit.h"

namespace cimrman::datamodel {

  OMHit::OMHit(int _OM_num)
  { 
    OM_num = _OM_num;
    //mainwall IT
    if(OM_num < 260) 
    {
      SWCR[0] = 0;
      SWCR[1] = -1;
      SWCR[2] = OM_num / 13;
      SWCR[3] = OM_num % 13;
    }
    //mainwall FR
    else if(OM_num < 520)
    {
      SWCR[0] = 1;
      SWCR[1] = -1;
      SWCR[2] = (OM_num - 260) / 13;
      SWCR[3] = (OM_num - 260) % 13;
    }
    //Xcalo IT	
    else if(OM_num < 584)
    {
      SWCR[0] = 0;
      SWCR[1] = (OM_num - 520) / 32;
      SWCR[2] = ((OM_num - 520) / 16) % 2;
      SWCR[3] = (OM_num -520) % 16;
    }
    //Xcalo FR
    else if(OM_num < 648)
    {
      SWCR[0] = 1;
      SWCR[1] = (OM_num - 520 - 64) / 32;
      SWCR[2] = ((OM_num - 520 - 64) / 16) % 2;
      SWCR[3] = (OM_num -520 - 64) % 16;
    }
    //GVeto IT
    else if(OM_num < 680)
    {
      SWCR[0] = 0;
      SWCR[1] = (OM_num - 520 - 128) / 16;
      SWCR[2] = (OM_num - 520 - 128) % 16;
      SWCR[3] = -1;
    }
    //GVeto FR
    else if(OM_num < 712)
    {
      SWCR[0] = 1;
      SWCR[1] = (OM_num - 520 - 128 - 32) / 16;
      SWCR[2] = (OM_num - 520 - 128 - 32) % 16;
      SWCR[3] = -1;
    }
    return;
  }
  

  OMHit::OMHit(int _SWCR[4])
  {
    for(int i = 0; i < 4; i++)
    {
      SWCR[i] = _SWCR[i];    
    }
    // auto detect MW
    if (SWCR[0] != -1 && 
	      SWCR[1] == -1 && 
	      SWCR[2] != -1 && 
	      SWCR[3] != -1 )
    {
      OM_num = 260*SWCR[0] + 13*SWCR[2] + SWCR[3];
    }
    // auto detect XW
    else if (SWCR[0] != -1 && 
	           SWCR[1] != -1 && 
	           SWCR[2] != -1 && 
	           SWCR[3] != -1 )
    {
      OM_num = 520 + 64*SWCR[0] + 32*SWCR[1] + 16*SWCR[2] + SWCR[3];
    }
    // auto detect GV
    else if (SWCR[0] != -1 && 
	           SWCR[1] != -1 && 
	           SWCR[2] != -1 && 
	           SWCR[3] == -1 )
    {
      OM_num = 520 + 128 + 32*SWCR[0] + 16*SWCR[1] + SWCR[2];
    }
    return;
  }

  void OMHit::set_OM_num(int _OM_num)
  {
    OM_num = _OM_num;
    return;
  }

  int OMHit::get_OM_num() const
  {
    return OM_num;
  }

  int* OMHit::get_SWCR()
  {	
    return SWCR;
  }

  const int* OMHit::get_SWCR() const
  {	
    return SWCR;
  }

  double OMHit::get_x() const
  {	
    return x;
  }
  
  double OMHit::get_y() const
  {	
    return y;
  }
  
  double OMHit::get_z() const
  {	
    return z;
  }
  
  
  void OMHit::set_x(double _x)
  {
    x = _x;
  }
      
  void OMHit::set_y(double _y)
  {
    y = _y;
  }
    
  void OMHit::set_z(double _z)
  {
    z = _z;
  }

  void OMHit::print(std::ostream & out_) const
  {
    out_ << "	OM " << SWCR[0] << "." << SWCR[1] << "." << SWCR[2] << "." << SWCR[3] << "." << std::endl;
    return;
  }

} //  end of namespace cimrman::datamodel
