#ifndef FALAISE_CIMRMAN_TRACKERHIT_H
#define FALAISE_CIMRMAN_TRACKERHIT_H

// Standard headers
#include <iostream>
#include <limits>
#include <memory>

#include "falaise/snemo/datamodels/calibrated_tracker_hit.h"
#include <bayeux/datatools/handle.h>

#include <datatools/utils.h>

namespace tkrec {

  class TrackerHit
  {    
  private:

    int SRL[3]; // side, column, layer
    int cell_num = -1;
  
    static const double default_R; // in mm
    static const double default_Z; // in mm
    static const double default_sigma_R; // in mm
    static const double default_sigma_Z; // in mm
    		
    double x = datatools::invalid_real();
    double y = datatools::invalid_real(); 
    double R = default_R;
    double Z = default_Z;
    double sigma_R = default_sigma_R;
    double sigma_Z = default_sigma_Z;
    bool valid_R = false;
    bool valid_Z = false;
    
    bool prompt_hit = false; 
    double delayed_time = datatools::invalid_real(); 
    
    datatools::handle<snemo::datamodel::calibrated_tracker_hit> CDbank_tr_hit;
    
  public:
    
    TrackerHit();
    TrackerHit(int _SRL[3]);
    TrackerHit(int _cell_num);
    
    void set_CDbank_tr_hit(const datatools::handle<snemo::datamodel::calibrated_tracker_hit> & _CDbank_tr_hit);
    const datatools::handle<snemo::datamodel::calibrated_tracker_hit> & get_CDbank_tr_hit() const;
    
    double get_x() const;
    double get_y() const;
    double get_R() const;
    double get_Z() const;
    double get_sigma_R() const;
    double get_sigma_Z() const;
    bool has_valid_R() const;
    bool has_valid_Z() const;
    
    void set_x(double _x);
    void set_y(double _y);
    void set_R(double _R);
    void set_Z(double _Z);
    void set_sigma_R(double _sigma_R = default_sigma_R);
    void set_sigma_Z(double _sigma_Z = default_sigma_Z);
    void set_valid_R();
    void set_valid_Z();
    void set_invalid_R();
    void set_invalid_Z();
    
    const int* get_SRL() const;
    int* get_SRL();
    int get_cell_num() const;
    
    double get_delayed_time() const;
    void set_delayed_time(double _delayed_time);
    
    bool is_prompt() const;
    void set_as_prompt();
    
    void update_drift_radius(double drift_time);

    void print(std::ostream & out_ = std::cout) const;

  };
  
  typedef std::shared_ptr<TrackerHit> TrackerHitHdl;
  typedef std::shared_ptr<const TrackerHit> ConstTrackerHitHdl;

} //  end of namespace tkrec

#endif // FALAISE_CIMRMAN_TRACKERHIT_H
