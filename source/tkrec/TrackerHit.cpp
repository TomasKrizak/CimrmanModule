// TK headers
#include "tkrec/TrackerHit.h"

#include <datatools/exception.h>

//ClassImp(tkrec::TrackerHit);

namespace tkrec {


  // dimensions in mm
  // origin in the center of detector
  const double TrackerHit::default_sigma_R = 2.0; // in mm //_config_.default_sigma_r; TODO??
  const double TrackerHit::default_sigma_Z = 17.0; // in mm
  const double TrackerHit::default_R = 22.0; // in mm
  const double TrackerHit::default_Z = 0.0; // in mm


  TrackerHit::TrackerHit()
  {
    SRL[0] = -1;
    SRL[1] = -1;
    SRL[2] = -1;
  }
  
  TrackerHit::TrackerHit(int _SRL[3])
  {
    cell_num = 113 * 9 * _SRL[0] + 9 * _SRL[1] + _SRL[2];
    DT_THROW_IF(cell_num > 2033 || cell_num < 0, std::logic_error,
	"Not valid tracker cell SRL combination (" << _SRL[0] << ", " << _SRL[1] << ", " << _SRL[2] << ")"); 
    for(int i = 0; i < 3; i++)
    {
      SRL[i] = _SRL[i];    
    }
  }
  
  TrackerHit::TrackerHit(int _cell_num)
  {
    DT_THROW_IF(_cell_num > 2033 || _cell_num < 0, std::logic_error,
		"Not valid tracker cell number " << _cell_num); 
    cell_num = _cell_num;
    
    SRL[0] =  cell_num / 1017; 	// compute side
    SRL[1] = (cell_num % 1017)/ 9; // compute row
    SRL[2] =  cell_num % 9; 	// compute layer
  }
 
  void TrackerHit::set_CDbank_tr_hit(const datatools::handle<snemo::datamodel::calibrated_tracker_hit> & _CDbank_tr_hit)
  {
    CDbank_tr_hit = _CDbank_tr_hit;
  }
  
  const datatools::handle<snemo::datamodel::calibrated_tracker_hit> & TrackerHit::get_CDbank_tr_hit() const
  {
    return CDbank_tr_hit;
  }

  double TrackerHit::get_x() const
  {
    return x;
  }

  double TrackerHit::get_y() const
  {
    return y;
  }
  
  double TrackerHit::get_R() const
  {
    return R;
  }
  
  double TrackerHit::get_Z() const
  {
    return Z;
  }
  
  double TrackerHit::get_sigma_R() const
  {
    return sigma_R;
  }
  
  double TrackerHit::get_sigma_Z() const
  {
    return sigma_Z;
  }
  
  bool TrackerHit::has_valid_Z() const
  {
    return valid_Z;
  }

  bool TrackerHit::has_valid_R() const
  {
    return valid_R;
  }

  double TrackerHit::get_delayed_time() const
  {
    return delayed_time;
  }
  
  void TrackerHit::set_delayed_time(double _delayed_time)
  {
    delayed_time = _delayed_time;
  }

  int* TrackerHit::get_SRL()
  {
    return SRL;
  }
  
  const int* TrackerHit::get_SRL() const
  {
    return SRL;
  }
  
  int TrackerHit::get_cell_num() const
  {
    return cell_num;
  }

  void TrackerHit::set_x(double _x)
  {
    x = _x;
  }
  
  void TrackerHit::set_y(double _y)
  {
    y = _y;
  }
  void TrackerHit::set_R(double _R)
  {
    R = _R;
  }
  
  void TrackerHit::set_Z(double _Z)
  {
    Z = _Z;
  }
  
  void TrackerHit::set_sigma_R(double _sigma_R)
  {
    sigma_R = _sigma_R;
  }

  void TrackerHit::set_sigma_Z(double _sigma_Z)
  {
    sigma_Z = _sigma_Z;
  }

  void TrackerHit::set_valid_R()
  {
    valid_R = true;
  }

  void TrackerHit::set_valid_Z()
  {
    valid_Z = true;
  }
  
  void TrackerHit::set_invalid_R()
  {
    valid_R = false;
  }

  void TrackerHit::set_invalid_Z()
  {
    valid_Z = false;
  }
  
  bool TrackerHit::is_prompt() const
  {
  	return prompt_hit;
  }
  
  void TrackerHit::set_as_prompt()
  {
  	prompt_hit = true;
  }

  void TrackerHit::update_drift_radius(double drift_time)
  {
    const double time_usec = drift_time / 1000.0;
    const double _tracker_drift_model_manu_params_[5] = {0.263223, -0.030965, -0.571594, 6.01392e-02, 1.13142e+03};
    double radius;
    if (time_usec < 0)
    {
      radius = 0;
      valid_R = false;
    }
    else
    {
      //if (time_usec > 30.0)
        //DT_LOG_WARNING(get_logging_priority(), " tracker hit with too large anode drift time = " << time_usec << " us");
      const double r2 = _tracker_drift_model_manu_params_[3] * std::log(1 + time_usec * _tracker_drift_model_manu_params_[4]);
      if (time_usec > 10.0)
        radius = r2;
      else
      {
        const double r1 = _tracker_drift_model_manu_params_[0] * std::exp(_tracker_drift_model_manu_params_[1] * time_usec) / std::pow(time_usec, _tracker_drift_model_manu_params_[2]);
        radius = std::min(r1, r2);
      }

      valid_R = true;
      R = radius * 44.0;
      sigma_R = 1.1;
    }
  }

  void TrackerHit::print(std::ostream & out_) const
  {
    out_ << "	SRL: " << SRL[0] << "." << SRL[1] << "." << SRL[2] << std::endl
	 << "	x = " << x << " mm"
	 << " ,y = " << y << " mm"
	 << " ,R = " << R << " mm"
	 << " ,Z = " << Z << " mm" << std::endl;
  }

} //  end of namespace tkrec
